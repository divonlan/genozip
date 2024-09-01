// ------------------------------------------------------------------
//   gencomp_piz.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include "zfile.h"
#include "buffer.h"
#include "dispatcher.h"
#include "compressor.h"
#include "file.h"
#include "strings.h"
#include "writer_private.h"
#include "piz.h"
#include "threads.h"

// PIZ main thread, genocat --show-plan (note: genounzip --show-plan uses a different code path)
void noreturn gencomp_piz_genocat_show_plan (void)
{
    for (int64_t i=0; i < z_file->recon_plan.len; i++) {
        ReconPlanItemP p = B(ReconPlanItem, z_file->recon_plan, i);

        if (p->flavor == PLAN_VB_PLAN) {
            evb->vblock_i = p->vb_i;                
            piz_read_vb_header (evb);
            sam_piz_after_vb_header (evb);    // adjust the vb->vb_plan Buffer
            gencomp_piz_vb_to_plan (evb, &i, true); // also displays vb_plan for this MAIN VB
            
            buf_free (evb->z_data);
            buf_free (evb->vb_plan);
        }
    }

    exit_ok;
}

// PIZ main thread: For each MAIN VB, set v->first_gc_vb/v->first_gc_vb_line to starting point of the 
// PRIM/DEPN data to be integrated back into this MAIN VB. The detailed vb integration plan (vb_plan)
// is in the payload of each VB_HEADER section
void gencomp_piz_initialize_vb_info (void) 
{
    START_TIMER;

    // get SEC_GENCOMP data
    ASSERTNOTINUSE (evb->scratch);
    zfile_get_global_section (SectionHeader, sections_first_sec (SEC_GENCOMP, HARD_FAIL), &evb->scratch, "scratch");
    
    evb->scratch.len /= sizeof (uint32_t);
    BGEN_u32_buf (&evb->scratch, NULL);

    evb->scratch.len /= (sizeof (GencompSecItem) / sizeof (uint32_t)); // to units of GencompSecItem;

    // Section of PRIM and DEPN VbHeader sections, ordered by vb_i
    for (CompIType comp_i=SAM_COMP_PRIM; comp_i <= SAM_COMP_DEPN; comp_i++) {

        Section gc_vb_header_sec=NULL;
        bool first_sec = true;
        uint32_t remaining_gc_lines_in_gc_vb=0, consume=0;
        
        // iterate on MAIN VBs in order of their absorption in ZIP (evb->scratch (=SEC_GENCOMP) is ordered by order of absorption)
        for_buf (GencompSecItem, main_vb_ent, evb->scratch) {
            VbInfo *v = VBINFO(main_vb_ent->vb_i);                    
            v->num_gc_lines[comp_i-1] = main_vb_ent->num_gc_lines[comp_i-1];

            uint32_t remaining_gc_lines_in_main = main_vb_ent->num_gc_lines[comp_i-1];

            if (remaining_gc_lines_in_gc_vb) // continuing a gc VB that start in an earlier MAIN VB
                v->num_gc_vbs++;

            // progress to next MAIN VB by "consuming" all gc data for this MAIN VB
            bool v_data_set = false;
            while (remaining_gc_lines_in_main) { // iterate to "consume" PRIM/DEPN lines
                if (!remaining_gc_lines_in_gc_vb) {
                    if (first_sec || gc_vb_header_sec) // make sure we don't circle back to the first section after already exhausting all sections
                        sections_get_next_vb_header_sec (comp_i, &gc_vb_header_sec); 
                    
                    ASSERT (gc_vb_header_sec, "We still need to consume %u %s lines in MAIN/%u, but %s data has exhausted prematurely", 
                            remaining_gc_lines_in_main, comp_name (comp_i), main_vb_ent->vb_i, comp_name (comp_i));

                    v->num_gc_vbs++;
                    remaining_gc_lines_in_gc_vb = gc_vb_header_sec->num_lines;
                    first_sec = false;
                }

                // set data at first gc line of this main VBs 
                if (!v_data_set) {
                    v->first_gc_vb       [comp_i-1] = gc_vb_header_sec->vblock_i;
                    v->first_gc_vb_line_i[comp_i-1] = gc_vb_header_sec->num_lines - remaining_gc_lines_in_gc_vb;

                    v_data_set = true;
                }

                consume = MIN_(remaining_gc_lines_in_main, remaining_gc_lines_in_gc_vb);
                remaining_gc_lines_in_gc_vb -= consume;
                remaining_gc_lines_in_main  -= consume;

                MAXIMIZE (VBINFO(gc_vb_header_sec->vblock_i)->ending_vb, v->vblock_i);
            }
        }

        // we finished all MAIN VBs, expecting gc data to be exhausted too:
        ASSERT (!remaining_gc_lines_in_gc_vb && 
                (  (consume && !sections_get_next_vb_header_sec (comp_i, &gc_vb_header_sec)) // either : gc lines were consumed in the last iteration, thus exhausting current gc VB and there are no more gc VBs
                || (!consume && !gc_vb_header_sec)), // or: final gc line consumed in a previous iteration and we already know that there are no more gc VBs
                "Finished assigning %s lines to all MAIN VBs, expecting %s lines to be all used up, but they are not", comp_name (comp_i), comp_name (comp_i));
    }

    if (flag.show_sec_gencomp) {
        iprint0 ("SEC_GENCOMP:\n");
        for_buf (GencompSecItem, gci, evb->scratch) {
            VbInfo *v = VBINFO(gci->vb_i);                    
            iprintf ("MAIN vb_i=%-4u num_prim_lines=%-4u num_depn_lines=%-4u Calculated start: PRIM/%u/%u DEPN/%u/%u\n", 
                     gci->vb_i, gci->num_gc_lines[0], gci->num_gc_lines[1], v->first_gc_vb[0], v->first_gc_vb_line_i[0], v->first_gc_vb[1], v->first_gc_vb_line_i[1]);
        }

        if (is_genocat) exit_ok;
    }

    buf_free (evb->scratch);

    // case: we're asked to output a single PRIM or DEPN VB 
    if (flag.one_vb && flag.one_vb <= z_file->num_vbs && IN_RANGX (VBINFO(flag.one_vb)->comp_i, SAM_COMP_PRIM, SAM_COMP_DEPN)) 
        writer_add_single_vb_plan (flag.one_vb);
    
    // case: plan of all MAIN VBs with embedded PRIM / DEPN lines as presribed
    else 
        writer_add_trivial_plan (SAM_COMP_MAIN, PLAN_VB_PLAN);

    if (is_genocat && flag.show_recon_plan && !flag_has_head && !flag_has_tail)
        gencomp_piz_genocat_show_plan (); // doesn't return

    COPY_TIMER_EVB (gencomp_piz_initialize_vb_info);
}

// PIZ main thread:used for expanding PLAN_VB_PLAN for calculation of --head, --tail or --lines
void gencomp_piz_expand_PLAN_VB_PLAN (ReconPlanItemP *p, int64_t *i, UpdatedIType updated_i_type)
{
    VBINFO((*p)->vb_i)->vb_plan_expanded = true;

    evb->vblock_i = (*p)->vb_i;                
    piz_read_vb_header (evb);
    sam_piz_after_vb_header (evb); // adjust the vb->vb_plan Buffer
    
    int64_t plan_len_before = z_file->recon_plan.len;

    CLEAR_FLAG (show_recon_plan);
    gencomp_piz_vb_to_plan (evb, i, false); // if --show-plan, also displays vb_plan for this MAIN VB
    RESTORE_FLAG (show_recon_plan);

    // *i is now at one before beginning. change to one after end if needed
    if (updated_i_type == AFTER_LAST_EXPANDED_ITEM)
        *i += z_file->recon_plan.len - plan_len_before + 2;

    // update p as buffer was realloced
    *p = B(ReconPlanItem, z_file->recon_plan, *i);

    buf_free (evb->z_data);
}

// PIZ writer thread: revise z_file->recon_plan to include plan sent from zip through vb->vb_plan.
// also remove already-executed steps, to keep memory footprint small.
void gencomp_piz_vb_to_plan (VBlockP vb, int64_t *i, bool remove_previous)
{
    START_TIMER;

    // if vb->vb_plan is empty, it means entire VB (if this MAIN, it also means this VB does not include PRIM/DEPN lines)
    // possibly not a MAIN VB, if requested a PRIM/DEPN VB with genocat --one-vb 
    if (!vb->vb_plan.len) {
        B(ReconPlanItem, z_file->recon_plan, *i)->flavor = PLAN_FULL_VB; // just replace PLAN_VB_PLAN
        
        if (flag.show_recon_plan)
            recon_plan_show (vb->vblock_i, *i, 1);

        (*i)--;
        goto done;
    }

    VbInfo *v = VBINFO(vb->vblock_i);                    

    // move the remaining part of recon_plan aside, and discard the already-executed part
    VBlockP scratch_vb = threads_am_i_main_thread() ? evb : wvb;
    ASSERTNOTINUSE (scratch_vb->scratch);
    buf_copy (scratch_vb, &scratch_vb->scratch, &z_file->recon_plan, ReconPlanItem, *i+1, 0, "scratch");

    // note: we calculated num_gc_vbs in gencomp_piz_initialize_vb_info because it is not possible to calculate the exact
    // num_ending_gc_vbs as gencomp_piz_initialize_vb_info iterates in absoption order, while PLAN_END_OF_VB appears in reconstruction order 
    buf_alloc (NULL, &z_file->recon_plan, 0, vb->vb_plan.len + (1 + v->num_gc_vbs)/*upper bound on number of PLAN_END_OF_VBs*/ + (remove_previous ? scratch_vb->scratch.len : (z_file->recon_plan.len - 1)), ReconPlanItem, 0, NULL);
    z_file->recon_plan.len = remove_previous ? 0 : *i;
    
    // iterators
    uint32_t main_next_line = 0;
    Section gc_vb_sec[2]    = {}; // VB_HEADER section of PRIM/DEPN vb from which we are currently consuming
    uint32_t gc_line_i[2]   = {}; // next line within PRIM/DEPN current vb

    CompIType prev_comp_i = SAM_COMP_PRIM; // initialize (same as in sam_zip_init_vb). note: we don't use vb->vb_plan.prev_comp_i bc we process the entire VB here, and this function might be called multiple times

    for_buf (VbPlanItem, vpi, vb->vb_plan) {
        CompIType comp_i = VB_PLAN_COMP_PIZ (vpi->comp, prev_comp_i);
        prev_comp_i = comp_i;

        if (!vpi->n_lines) // an entry only for transitioning prev_comp_i - no plan entry
            continue;

        if (comp_i == SAM_COMP_MAIN) {
            ASSERT (vpi->n_lines + main_next_line <= vb->lines.len32,
                    "%s: expecting vpi->n_lines=%u + main_next_line=%u <= vb->lines.len32=%u", VB_NAME, vpi->n_lines, main_next_line, vb->lines.len32);
        
            BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
                .vb_i       = vb->vblock_i,
                .start_line = main_next_line, 
                .num_lines  = vpi->n_lines,
                .flavor     = PLAN_RANGE 
            };

            main_next_line += vpi->n_lines;
        }

        else {
            int c = comp_i - 1; // 0=PRIM 1=DEPN

            // consume lines from possibly multiple VBs of comp_i
            uint32_t vpi_n_lines = vpi->n_lines; // a copy so to not destroy the original 
            while (vpi_n_lines) {
                // case: current MAIN vb has not yet consumed any lines from this component 
                if (!gc_vb_sec[c]) { 
                    gc_vb_sec[c] = sections_vb_header (v->first_gc_vb[c]);
                    gc_line_i[c] = v->first_gc_vb_line_i[c];
                }

                // case: we have consumed all lines from current gc vb - move to next gc vb
                else if (gc_line_i[c] == gc_vb_sec[c]->num_lines) {
                    // case: this MAIN vb is the last VB (by vb_i, not order of absorption!) that consumes the gc VB, so add a PLAN_END_OF_VB
                    if (VBINFO(gc_vb_sec[c]->vblock_i)->ending_vb == vb->vblock_i) 
                        BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
                            .vb_i   = gc_vb_sec[c]->vblock_i,
                            .flavor = PLAN_END_OF_VB 
                        };

                    sections_get_next_vb_header_sec (comp_i, &gc_vb_sec[c]);
                    gc_line_i[c] = 0;
                }
                
                uint32_t n_consumed = MIN_(gc_vb_sec[c]->num_lines - gc_line_i[c], vpi_n_lines);
                        
                BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
                    .vb_i       = gc_vb_sec[c]->vblock_i,
                    .start_line = gc_line_i[c],
                    .num_lines  = n_consumed,
                    .flavor     = PLAN_RANGE 
                };
            
                vpi_n_lines  -= n_consumed;
                gc_line_i[c] += n_consumed;                  
            }
        }
    }

    for (int c=SAM_COMP_PRIM-1; c <= SAM_COMP_DEPN-1; c++)
        // case: this MAIN vb is last VB (by vb_i, not order of absorption!) that consumes gc VB is last in its plan (gc VBs that were switched out were considered earlier)
        if (gc_vb_sec[c] && VBINFO(gc_vb_sec[c]->vblock_i)->ending_vb == vb->vblock_i) 
            BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
                .vb_i   = gc_vb_sec[c]->vblock_i,
                .flavor = PLAN_END_OF_VB 
            };

    // our MAIN VB also ends here
    BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
        .vb_i   = vb->vblock_i,
        .flavor = PLAN_END_OF_VB 
    };

    if (flag.show_recon_plan)
        recon_plan_show (vb->vblock_i, 0, z_file->recon_plan.len); // recon plan for this MAIN VB only

    // return subsequent MAIN VBs to plan
    buf_append_buf (NULL, &z_file->recon_plan, &scratch_vb->scratch, ReconPlanItem, NULL);
    buf_free (scratch_vb->scratch);

    ASSERT (main_next_line == vb->lines.len32, 
            "Only %u of %u MAIN lines were consumed when creating plan for %s", main_next_line, vb->lines.len32, VB_NAME);

    *i = remove_previous ? -1 : (*i - 1); // main loop to continue to beginning of updated recon_plan 

done:
    COPY_TIMER (gencomp_piz_vb_to_plan);
}

static void gencomp_piz_add_one_gc_vb_to_rl (VBIType gc_vb_i, bool *initialized)
{
    VbInfo *v = VBINFO(gc_vb_i); // vb_info of gc vb         

    if (!v->encountered) { // first encounter with this VB - add it to the reading list
        if (! *initialized) { // defer initialization until we know that we actually need to change the reading list
            // save yet-read suffix of readling list - we will append them back when we're done expanding the placeholder
            buf_copy (evb, &evb->scratch, &z_file->piz_reading_list, SectionEnt, z_file->piz_reading_list.next, 0, "scratch");
            z_file->piz_reading_list.len = 0; // discard all sections that were already read
            *initialized = true;
        }

        sections_reading_list_add_vb_header (gc_vb_i);
        v->encountered = true;
    }
}

// PIZ main thread: called to expand the reading list to include all of a MAIN VB's dependent 
// PRIM and DEPN VBs that have not been read yet
void gencomp_piz_update_reading_list (VBlockP vb)
{
    START_TIMER;

    VbInfo *v = VBINFO(vb->vblock_i); // vb_info of MAIN vb               

    // VB was already expanded and PRIM/DEPN VBs were already added to the reading list (or maybe filtered out) 
    // (this happens when PLAN_VB_PLAN was expanded during z_file->recon_plan creation, in order to filter out lines with --head/tail/lines)
    if (v->vb_plan_expanded) return;

    // iterators
    Section gc_vb_sec[2]  = {}; // VB_HEADER section of PRIM/DEPN vb from which we are currently consuming
    uint32_t gc_line_i[2] = {}; // next line within PRIM/DEPN current vb
    bool initialized      = false;

    CompIType prev_comp_i = SAM_COMP_PRIM; // initialize (same as in sam_zip_init_vb). note: we don't use vb->vb_plan.prev_comp_i bc we process the entire VB here

    for_buf (VbPlanItem, vpi, vb->vb_plan) {
        CompIType comp_i = VB_PLAN_COMP_PIZ (vpi->comp, prev_comp_i);
        prev_comp_i = comp_i;

        if (comp_i != SAM_COMP_MAIN) {
            int c = comp_i - 1;

            // consume lines from possibly multiple VBs of comp_i
            uint32_t vpi_n_lines = vpi->n_lines; // a copy so to not destroy the original 
            while (vpi_n_lines) {
                // case: current MAIN vb has not yet consumed any lines from this component 
                if (!gc_vb_sec[c]) { 
                    ASSERT (v->first_gc_vb[c], "v->first_gc_vb[%s]=0", comp_name(c+1));

                    gc_vb_sec[c] = sections_vb_header (v->first_gc_vb[c]);
                    gc_line_i[c] = v->first_gc_vb_line_i[c];
                    gencomp_piz_add_one_gc_vb_to_rl (gc_vb_sec[c]->vblock_i, &initialized);
                }

                // case: we have consumed all lines from current gc vb - move to next gc vb
                else if (gc_line_i[c] == gc_vb_sec[c]->num_lines) {
                    sections_get_next_vb_header_sec (comp_i, &gc_vb_sec[c]);
                    gc_line_i[c] = 0;
                    gencomp_piz_add_one_gc_vb_to_rl (gc_vb_sec[c]->vblock_i, &initialized);
                }
                
                uint32_t n_consumed = MIN_(gc_vb_sec[c]->num_lines - gc_line_i[c], vpi_n_lines);
                                    
                vpi_n_lines  -= n_consumed;
                gc_line_i[c] += n_consumed;
            }
        }
    }

    // if we made any changes: append the yet-unread sections back to the reading list
    if (initialized) {
        if (flag.show_reading_list) {
            iprintf ("Inserting to reading_list after reading SEC_VB_HEADER of vblock_i=%u :\n", vb->vblock_i);
            sections_show_section_list (z_file->data_type, &z_file->piz_reading_list);
        }

        buf_append_buf (evb, &z_file->piz_reading_list, &evb->scratch, SectionEnt, NULL);
        buf_free (evb->scratch);

        z_file->piz_reading_list.next = 0; // piz_one_txt_file to begin the updated header list at its beginning
    }

    COPY_TIMER (gencomp_piz_update_reading_list);
}
