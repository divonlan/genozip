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

// For each MAIN VB, set v->first_gc_vb/v->first_gc_vb_line to starting point of the 
// PRIM/DEPN data to be integrated back into this MAIN VB. The detailed vb integration plan (vb_plan)
// is in the payload of each VB_HEADER section
void gencomp_piz_initialize_vb_info (void) 
{
    // get SEC_GENCOMP data
    ASSERTNOTINUSE (evb->scratch);
    zfile_get_global_section (SectionHeader, sections_get_gencomp_sec(), &evb->scratch, "scratch");
    
    evb->scratch.len /= sizeof (uint32_t);
    BGEN_u32_buf (&evb->scratch, NULL);

    evb->scratch.len /= (sizeof (GencompSecItem) / sizeof (uint32_t)); // to units of GencompSecItem;

    // Section of PRIM and DEPN VbHeader sections, ordered by vb_i
    for (CompIType comp_i=SAM_COMP_PRIM; comp_i <= SAM_COMP_DEPN; comp_i++) {

        Section gc_vb_header_sec=NULL;
        bool first_sec = true;
        uint32_t remaining_gc_lines_in_gc_vb=0, consume=0;
        
        // iterate on MAIN VBs in order of their absorption in ZIP (evb->scratch (=SEC_GENCOMP) is ordered by order of absorption)
        for_buf (GencompSecItem, main_vb_info, evb->scratch) {
            VbInfo *v = VBINFO(main_vb_info->vb_i);                    

            uint32_t remaining_gc_lines_in_main = main_vb_info->num_gc_lines[comp_i-1];

            // progress to next MAIN VB by "consuming" all gc data for this MAIN VB
            bool v_data_set = false;
            while (remaining_gc_lines_in_main) { // iterate to "consume" PRIM/DEPN lines
                if (!remaining_gc_lines_in_gc_vb) {
                    if (first_sec || gc_vb_header_sec) // make sure we don't circle back to the first section after already exhausting all sections
                        sections_get_next_vb_header_sec (comp_i, &gc_vb_header_sec); 
                    
                    ASSERT (gc_vb_header_sec, "We still need to consume %u %s lines in MAIN/%u, but %s data has exhausted prematurely", 
                            remaining_gc_lines_in_main, comp_name (comp_i), main_vb_info->vb_i, comp_name (comp_i));

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

    writer_add_trivial_plan (SAM_COMP_MAIN, PLAN_VB_PLAN);
}

// revise z_file->recon_plan to include plan sent from zip through vb->vb_plan.
// also remove already-executed steps, to keep memory footprint small.
void gencomp_piz_vb_to_plan (VBlockP vb, int64_t *i)
{
    ASSERT0 (Z_DT(SAM) && vb->comp_i == SAM_COMP_MAIN, "Invalid dt or comp_i");

    // if vb->vb_plan is empty, it means entire VB is MAIN
    if (!vb->vb_plan.len) {
        B(ReconPlanItem, z_file->recon_plan, *i)->flavor = PLAN_FULL_VB; // just replace PLAN_VB_PLAN
        (*i)--;
        return;
    }

    VbInfo *v = VBINFO(vb->vblock_i);                    

    // move the remaining part of recon_plan aside, and discard the already-executed part
    ASSERTNOTINUSE (wvb->scratch);
    buf_copy (wvb, &wvb->scratch, &z_file->recon_plan, ReconPlanItem, *i+1, 0, "scratch");

    buf_alloc (NULL, &z_file->recon_plan, 0, vb->vb_plan.len + v->num_gc_vbs/*potential PLAN_END_OF_VBs*/ + wvb->scratch.len, ReconPlanItem, 0, NULL);
    z_file->recon_plan.len = 0;
    
    // iterators
    uint32_t main_next_line = 0;
    Section gc_vb_sec[2]    = {}; // VB_HEADER section of PRIM/DEPN vb from which we are currently consuming
    uint32_t gc_line_i[2]   = {}; // next line within PRIM/DEPN current vb

    for_buf (VbPlanItem, vpi, vb->vb_plan) {
        if (vpi->comp_i == SAM_COMP_MAIN) {
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
            int c = vpi->comp_i - 1;

            // consume lines from possibly multiple VBs of vpi->comp_i
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

                    sections_get_next_vb_header_sec (vpi->comp_i, &gc_vb_sec[c]);
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
        // case: this MAIN vb is last VB (by vb_i, not order of absorption!) that consumes gc VB in its is last in its plan (gc VBs that were switched out were considered earlier)
        if (VBINFO(gc_vb_sec[c]->vblock_i)->ending_vb == vb->vblock_i) 
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
        recon_plan_show (vb->vblock_i); // recon plan for this MAIN VB only

    // return subsequent MAIN VBs to plan
    buf_append_buf (NULL, &z_file->recon_plan, &wvb->scratch, ReconPlanItem, NULL);
    buf_free (wvb->scratch);

    ASSERT (main_next_line == vb->lines.len32, 
            "Only %u of %u MAIN lines were consumed when creating plan for %s", main_next_line, vb->lines.len32, VB_NAME);

    *i = -1; // main loop to continue to beginning of updated recon_plan 
//xxx verify that if no subsetting, all prim/depnn was consumed
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

        sections_txt_list_add_vb_header (gc_vb_i);
        v->encountered = true;
    }
}

// called by main thread to expand a MAIN Vb to include all its dependent PRIM and DEPN 
// VBs that have not been read yet
void gencomp_piz_update_piz_reading_list (VBlockP vb)
{
    VbInfo *v = VBINFO(vb->vblock_i); // vb_info of MAIN vb               

    // iterators
    Section gc_vb_sec[2]  = {}; // VB_HEADER section of PRIM/DEPN vb from which we are currently consuming
    uint32_t gc_line_i[2] = {}; // next line within PRIM/DEPN current vb
    bool initialized      = false;

    for_buf (VbPlanItem, vpi, vb->vb_plan) {
        if (vpi->comp_i != SAM_COMP_MAIN) {
            int c = vpi->comp_i - 1;

            // consume lines from possibly multiple VBs of vpi->comp_i
            uint32_t vpi_n_lines = vpi->n_lines; // a copy so to not destroy the original 
            while (vpi_n_lines) {
                // case: current MAIN vb has not yet consumed any lines from this component 
                if (!gc_vb_sec[c]) { 
                    gc_vb_sec[c] = sections_vb_header (v->first_gc_vb[c]);
                    gc_line_i[c] = v->first_gc_vb_line_i[c];
                    gencomp_piz_add_one_gc_vb_to_rl (gc_vb_sec[c]->vblock_i, &initialized);
                }

                // case: we have consumed all lines from current gc vb - move to next gc vb
                else if (gc_line_i[c] == gc_vb_sec[c]->num_lines) {
                    sections_get_next_vb_header_sec (vpi->comp_i, &gc_vb_sec[c]);
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
        buf_append_buf (evb, &z_file->piz_reading_list, &evb->scratch, SectionEnt, NULL);
        buf_free (evb->scratch);

        z_file->piz_reading_list.next = 0; // piz_one_txt_file to begin the updated header list at its beginning
    }
}
