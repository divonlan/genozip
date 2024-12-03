// ------------------------------------------------------------------
//   writer.c
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "context.h"
#include "profiler.h"
#include "strings.h"
#include "threads.h"
#include "sections.h"
#include "random_access.h"
#include "writer.h"
#include "writer_private.h"
#include "zfile.h"
#include "mgzip.h"
#include "mutex.h"
#include "codec.h" 
#include "endianness.h"
#include "bits.h"
#include "version.h"
#include "gencomp.h"
#include "piz.h"
#include "dispatcher.h"
#include "progress.h"
#include "buf_list.h"

// ---------------
// Data structures
// ---------------

rom recon_plan_flavors[8] = PLAN_FLAVOR_NAMES;

#define IS_DROPPED_BUF_NAME "is_dropped"

#define TXTINFO(comp_i) ({ ASSERTINRANGE ((comp_i), 0, txt_header_info.len); B(VbInfo, txt_header_info, (comp_i)); })

// #define vb_info         z_file->vb_info[0] // we use only [0] on the PIZ side
#define txt_header_info z_file->txt_header_info

VBlockP wvb = NULL;

static ThreadId writer_thread = THREAD_ID_NONE;

static void writer_release_vb (VbInfo *v)
{
    if (v->vb) vb_release_vb (&v->vb, PIZ_TASK_NAME);
    buf_destroy (v->is_dropped_buf); // destroy and not free as this memory is not going to be recycled - its a buffer inside a buffer
}

// --------------------------------
// VBlock and Txt Header properties
// --------------------------------

// Note: txtheader lines are used only for error reporting, using writer_get_txt_line_i. They are not
// included in recon_plan, as they are not counted for --head, --downsample etc
void writer_set_num_txtheader_lines (CompIType comp_i, uint32_t num_txtheader_lines)
{
    TXTINFO(comp_i)->num_lines = num_txtheader_lines;
}

// calculate 1-based textual line in txt_file (including txt_header), based on vb_i, vb->line_i and recon_plan
// or record in binary file (eg alignment in BAM). In case of a text file, it also includes the header lines.
// For now: returns 0 if reconstruction modifies the file, and therefore recon_plan doesn't reconstruct the original file
// To do: return line according to recon plan even if modified. challenge: lines dropped by reconstructor and not known yet to writer
int64_t writer_get_txt_line_i (VBlockP vb, LineIType line_in_vb/*0-based (if MAIN VB, excludes gencomp lines) */)
{
    if (flag.maybe_lines_dropped_by_reconstructor) return -1;
    if (!vb->lines.len32) return 0;

    // since data is unmodified, we have only FULL_VB, RANGE, END_OF_VB and TXTHEADER plan items 
    int64_t txt_num_lines  = 0;
    LineIType vb_num_lines = 0;
    int64_t factor = VB_DT(FASTQ) ? 4 : 1; // in FASTQ, each plan line is 4 txt lines

    // if vb is PRIM/DEPN - find the MAIN VB embedding the desired line and expand it
    if (z_has_gencomp && vb->comp_i != COMP_MAIN) {
        // TODO (bug 1137): traverse MAIN VBs in VbInfo in order of absorption(!), looking for the VB embedding the requested prim/depn line:
        // 1. set vb_num_lines to the sum of num_gc_lines[vb->comp_i] of all MAIN VBs prior (in absorption order!) 
        //    to the VB embedding the the requested PRIM/DEPN line  
        // 2. find the VB in which the requested prim/depn line: call a version of gencomp_piz_expand_PLAN_VB_PLAN 
        //    that can run from an arbitrary thread to expand the PLAN_VB_PLAN in recon_plan
        // -. We probably will need update gencomp_piz_initialize_vb_info to add more fields to VbInfo to do this
        return -999; // not implemented
    }
    
    for (int64_t i=0; i < z_file->recon_plan.len; i++) {
        ReconPlanItemP p = B(ReconPlanItem, z_file->recon_plan, i);
    
        // case: interleaved plan  
        if (p->flavor == PLAN_INTERLEAVE) {
            if (p->vb_i == vb->vblock_i && p->num_lines/2 > line_in_vb) // R1
                return txt_num_lines + 2 * line_in_vb * factor + 1; // +1 because 1-based 

            else if (p->vb2_i == vb->vblock_i && p->num_lines/2 > line_in_vb) // R2
                return txt_num_lines + (2 * line_in_vb + 1) * factor + 1; // +1 because 1-based
        }

        // case: MAIN VB containing the requested line in a gencomp file - expand PLAN_VB_PLAN first
        else if (p->flavor == PLAN_VB_PLAN && p->vb_i == vb->vblock_i) {
            return -999; // not implemented 

            // TODO (bug 1137): we need a version of gencomp_piz_expand_PLAN_VB_PLAN that can run from an arbitrary 
            // thread (not just main or writer) - this will require using a copy of z_file->recon_plan, 
            // because we can't change z_file->recon_plan itself under the feet of the writer thread
            
            // gencomp_piz_expand_PLAN_VB_PLAN (&p, &i, BEFORE_FIRST_EXPANDED_ITEM);
            // continue;
        } 

        else if (p->vb_i == vb->vblock_i) {
            // case: non-interleaved plan item containing this line
            if (vb_num_lines + p->num_lines > line_in_vb) 
                return txt_num_lines + (line_in_vb - vb_num_lines) * factor + 1; // +1 because 1-based 

            // case: plan item is from current VB, but we have not yet reached the current line
            else
                vb_num_lines += p->num_lines;
        }

        // note: txtheader lines are not included in the plan item, bc they are not counted for --header, --downsample etc
        bool no_write = (p->vb_i && !writer_does_vb_need_write (p->vb_i)) ||
                        (p->flavor == PLAN_TXTHEADER && !writer_does_txtheader_need_write (p->comp_i));

        if (!no_write)
            txt_num_lines += (p->flavor == PLAN_TXTHEADER) ? TXTINFO(p->comp_i)->num_lines
                                                           : (factor * p->num_lines); 
    }

    WARN_ONCE ("While Unexpectedly, unable to find current line %s/%u in recon_plan", VB_NAME, line_in_vb);
    return 0;
}

// returns my own number in the pair (1 or 2) and pair_vb_i. return 0 if file is not paired.
bool writer_am_i_pair_2 (VBIType vb_i, uint32_t *pair_1_vb_i)
{
    VbInfo *v = VBINFO(vb_i);
    VbInfo *t = TXTINFO(v->comp_i);

    if (!t->pair) return 0; // not paired

    if (v->pair_vb_i < vb_i) {
        *pair_1_vb_i = v->pair_vb_i;
        return true;
    }
    else    
        return false;
}

bool writer_does_txtheader_need_write (CompIType comp_i)
{
    bool needs_write = !flag.no_writer_thread && TXTINFO(comp_i)->needs_write;
    return needs_write;                             
}

bool writer_does_txtheader_need_recon (CompIType comp_i)
{
    bool needs_recon = TXTINFO(comp_i)->needs_recon;
    return needs_recon;                             
}

bool writer_get_fasta_contig_grepped_out (VBIType vb_i)
{
    return VBINFO(vb_i)->fasta_contig_grepped_out;
}

void writer_set_fasta_contig_grepped_out (VBIType vb_i)
{
    VBINFO(vb_i)->fasta_contig_grepped_out = true;
}

bool writer_does_vb_need_recon (VBIType vb_i)
{
    bool needs_recon = !z_file->vb_info[0].len || VBINFO(vb_i)->needs_recon;
    return needs_recon;
}

bool writer_does_vb_need_write (VBIType vb_i)
{
    bool needs_write = !z_file->vb_info[0].len || VBINFO(vb_i)->needs_write;
    return needs_write;
}

// called by: 1. PIZ main thread, when writer creates plan and filters by lines/head/tail
// 2. PIZ compute thread
Bits *writer_get_is_dropped (VBIType vb_i)
{
    VbInfo *v = VBINFO(vb_i);

    // allocate if needed. buffer was put in wvb buffer_list by writer_z_initialize
    if (!v->is_dropped && v->needs_write) {
        buf_alloc_bits_exact (NULL, &v->is_dropped_buf, v->num_lines, CLEAR, 0, IS_DROPPED_BUF_NAME);
        v->is_dropped = (BitsP)&v->is_dropped_buf;
    }

    return v->is_dropped;
}   

// -----------------------------
// PIZ reconstruction plan stuff
// -----------------------------

// PIZ main thread
static void writer_init_txt_header_info (void)
{
    buf_alloc_exact_zero (evb, txt_header_info, z_file->num_components, VbInfo, "txt_header_info"); // info on txt_headers

    Section sec = NULL;

    for (CompIType comp_i=0; comp_i < txt_header_info.len; comp_i++) {

        VbInfo *comp = TXTINFO(comp_i);
        comp->comp_i = comp_i;

        sec = sections_get_comp_txt_header_sec (comp_i);
        if (!sec) continue; // no SEC_TXT_HEADER section for this component (eg PRIM/DEPN components of SAM/BAM)

        // conditions entire txt header should be read 
        comp->needs_recon = 
            (  !Z_DT(SAM) // This clause only limits SAM/BAM files  
            || comp_i == SAM_COMP_MAIN // SAM file is reconstructed
            || (comp_i >= SAM_COMP_FQ00 && flag.deep)); // FQ files are reconstructed only if flag.deep set in flags_update_piz_one_z_file

        // conditions we write the txt header (doesn't affect the VBs of this component)
        comp->needs_write =
           !flag_loading_auxiliary
        && comp->needs_recon
        && !flag.no_header
        && DTPZ(txt_header_required) // not HDR_NONE
        && (!Z_DT(SAM) || comp_i == SAM_COMP_MAIN); // For SAM/BAM, show the txt_header only of MAIN, not PRIM/DEPN and (deep) FASTQ components           

        if (comp->needs_write && !flag.no_writer_thread) {
            // mutex: locked:    here (at initialization)
            //        waited on: writer thread, wanting the data
            //        unlocked:  by main thread after txtheader data is handed over.
            mutex_initialize (comp->wait_for_data);
            mutex_lock (comp->wait_for_data);
        } 

        if (VER(15) && (Z_DT(SAM) || Z_DT(FASTQ)))
            comp->pair = sec->flags.txt_header.pair;

        else if (!VER(15) && flag.pair) 
            comp->pair = comp_i ? PAIR_R2 : PAIR_R1; // only happens for FASTQ
    }
}

// PIZ main thread - initialize vb, component, txt_file info. this is run once per z_file
void writer_z_initialize (void)
{
    if (txt_header_info.len) return; // already initialized

    wvb = vb_initialize_nonpool_vb (VB_ID_WRITER, DT_NONE, WRITER_TASK_NAME);

    writer_init_txt_header_info(); 

    buf_alloc_exact_zero (evb, z_file->vb_info[0], z_file->num_vbs+1, VbInfo, "z_file->vb_info"); // +1 as first vb_i=1 (entry 0 will be unused) so we have num_vbs+1 entries

    for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {

        Section sec = sections_vb_header (vb_i);
            
        VbInfo *v = VBINFO(vb_i); 
        VbInfo *t = TXTINFO(sec->comp_i);

        v->vblock_i = vb_i;
        v->comp_i   = sec->comp_i;
        
        // since v14, num_lines is carried by Section. 
        if (VER(14))
            v->num_lines = sec->num_lines;

        // up to v13, num_lines was in SectionHeaderVbHeader
        else {  
            SectionHeaderVbHeader header = zfile_read_section_header (evb, sec, SEC_VB_HEADER).vb_header;
            v->num_lines = BGEN32 (header.v13_top_level_repeats);
        }

        ASSERT (!t->pair || flag.pair == PAIRED, "t->pair=%s inconsistent with flag.pair=%s", pair_type_name (t->pair), pair_type_name (flag.pair));

        // set pairs (used even if not interleaving)
        if (t->pair == PAIR_R1) 
            v->pair_vb_i = vb_i + sections_get_num_vbs (v->comp_i);
        
        else if (t->pair == PAIR_R2) {
            v->pair_vb_i = vb_i - sections_get_num_vbs (v->comp_i);
            ASSERTINRANGX (v->pair_vb_i, 1, z_file->num_vbs);

            // verify that this VB and its pair have the same number of lines (test when initiazing the second one)
            uint32_t R1_num_lines = VBINFO(v->pair_vb_i)->num_lines;
            ASSERT (v->num_lines == R1_num_lines, "vb=%u has %u lines, but vb=%u, its pair, has %u lines", 
                    vb_i, v->num_lines, v->pair_vb_i, R1_num_lines);
        }

        // conditions this VB should not be read or reconstructed 
        v->needs_recon = true; // default
        #define DROP v->needs_recon = false

        // --one-vb: user only wants to see a single VB, and this is not it
        if (flag.one_vb && flag.one_vb != vb_i && 
            !(flag.deep_fq_only && v->comp_i <= SAM_COMP_PRIM)) // but: if deep and requested VB is a FASTQ VB, don't drop MAIN/PRIM SAM VBs required for creating deep_ents 
            DROP; 

        // --header-only: drop all VBs (except VCF - handled separately below, 
        // and except FASTQ and FASTA, for which --header-only sets flag.header_only_fast, not flag.header_only)
        else if (flag.header_only && z_file->data_type != DT_VCF) DROP; 

        else switch (z_file->data_type) {

            case DT_FASTQ:
                // --R1 or --R2: we don't show the other component
                if (flag.one_component && flag.one_component-1 != v->comp_i) DROP;

                break;

            case DT_VCF:
                // --header-only VCF - drop VBs
                if (flag.header_only) DROP;

                break;

            case DT_SAM:
                // Deep file, with --sam or --bam - we don't need the FQ data
                if (flag.deep_sam_only && (v->comp_i >= SAM_COMP_FQ00)) DROP; 

                // Deep, reconstructing only FASTQ data: We don't need DEPN as it doesn't particpate in Deep data
                else if (flag.deep_fq_only && v->comp_i == SAM_COMP_DEPN) DROP;

                // Deep: --R1 --R2: we don't need the non-requested FASTQ components
                else if (flag.deep && flag.one_component-1 >= SAM_COMP_FQ00 &&
                    v->comp_i >= SAM_COMP_FQ00 && v->comp_i != flag.one_component-1) DROP;

                break;

            case DT_FASTA:
                if (!fasta_piz_is_vb_needed (vb_i)) DROP;
                break;
            
            default: {} // no special handling of other data types
        }

        #undef DROP

        // we add is_dropped_buf to the wvb buffer list. allocation will occur in main thread when writer 
        // create the plan, for VBs on the boundary of head/tail/lines, otherwise in reconstruction compute thread.
        if (v->needs_recon && (flag.maybe_lines_dropped_by_reconstructor || flag.maybe_lines_dropped_by_writer))
            buf_set_promiscuous_do (wvb, &v->is_dropped_buf, "is_dropped_buf", __FUNCLINE);
        
        // conditions in which VB should be written. if false, but needs_recon, VB is still read, but not written (eg reading aux files)
        v->needs_write = 
            v->needs_recon
        &&  !flag_loading_auxiliary // we're ingesting, but not reconstructing, an auxiliary file
        &&  !(flag.deep && flag.deep_fq_only && v->comp_i <= SAM_COMP_DEPN) // in Deep: --R1/--R2/--interleaved/--fq: reconstruct MAIN/PRIM/DEPN components but don't write them
        &&  !(flag.deep && flag.no_writer && !flag.show_recon_plan && v->comp_i >= SAM_COMP_FQ00); // in Deep with no_writer (eg --test): no need to write FASTQ components as unlike gencomp SAM components, these are digested in the compute thread, not in writer

        if (v->needs_write && !flag.no_writer_thread) {
            // mutex: locked:    here (at initialization)
            //        waited on: writer thread, wanting the data
            //        unlocked:  by main thread after vb data is handed over.
            mutex_initialize (v->wait_for_data);
            mutex_lock (v->wait_for_data);
        }

        else if (!v->needs_write)
            writer_release_vb (v);            
    }
}

static int64_t writer_get_plan_num_lines (void)
{
    int64_t num_lines = 0;
    for_buf (ReconPlanItem, p, z_file->recon_plan) 
        // note: in PLAN_TXT_HEADER num_lines is always 0, see writer_add_txtheader_plan
        if (!p->vb_i || VBINFO(p->vb_i)->needs_write) // exclude VBs that are not written (e.g. Deep SAM VB if writing a FASTQ component)
            num_lines += p->num_lines; // note: for VB_PLAN_VB, num_lines equals the sum of lines when expanded

    return num_lines;
}

static inline void writer_drop_lines (VbInfo *v, uint32_t start_line, uint32_t num_lines)
{
    if (!v->is_dropped)
        v->is_dropped = writer_get_is_dropped (v->vblock_i);

    #define DROPL(l) bits_set (v->is_dropped, start_line + (l))
    
    // shortcuts for the most common cases
    switch (num_lines) {
        case 5  : DROPL(4); // fall through
        case 4  : DROPL(3); 
        case 3  : DROPL(2); 
        case 2  : DROPL(1); 
        case 1  : DROPL(0); break;
        default : bits_set_region (v->is_dropped, start_line, num_lines);
    }

    #undef DROPL
}

// PIZ main thread
static void writer_drop_plan_item (ReconPlanItemP p, uint32_t drop_flavor)
{
    VbInfo *v = VBINFO(p->vb_i);

    switch (p->flavor) {
        case PLAN_RANGE:
            if (!v->is_dropped)
                v->is_dropped = writer_get_is_dropped (p->vb_i);

            writer_drop_lines (v, p->start_line, p->num_lines);
            p->flavor = drop_flavor;
            
            break;

        case PLAN_FULL_VB:
        case PLAN_VB_PLAN:
            v->needs_recon = false;
            p->flavor = drop_flavor;
            break;
        
        case PLAN_INTERLEAVE:
            v->needs_recon = VBINFO(p->vb2_i)->needs_recon = false;
            p->flavor = drop_flavor;
            break;
    
        default: {}  // TXTHEADER is never removed by head/tail/lines, 
                     // END_OF_VB will be removed by writer_cleanup_recon_plan_after_genocat_filtering if needed
    }
}

static void writer_drop_item_lines_from_start (ReconPlanItemP  p, uint32_t drop_lines)
{
    VbInfo *v = VBINFO(p->vb_i);

    switch (p->flavor) {
        case PLAN_RANGE:
            writer_drop_lines (v, p->start_line, drop_lines);
            break;

        case PLAN_FULL_VB:
            writer_drop_lines (v, 0, drop_lines);
            break;
        
        case PLAN_INTERLEAVE:
            writer_drop_lines (v, 0, (drop_lines+1) / 2);            // round up in case of odd number
            writer_drop_lines (VBINFO(p->vb2_i), 0, drop_lines / 2); // round down in case of odd number

            // case: pair-1 is removed in its entirety and pair-2 only the last line surviving: update plan item to single line
            if (drop_lines + 1 == p->num_lines) {             
                    uint32_t plan_i = BNUM (z_file->recon_plan, p);
                    ReconPlanItem old_p = *p;

                    // We still need to keep pair-1 in the section list, because we need it for reconstructing pair-2. 
                    // So we insert a zero-line PLAN_RANGE for pair-1
                    buf_insert (NULL, z_file->recon_plan, ReconPlanItem, plan_i, 
                                ((ReconPlanItem[]){ { .flavor = PLAN_RANGE,     .vb_i = old_p.vb_i,  .num_lines = 0, .start_line = 0                       },  
                                                    { .flavor = PLAN_END_OF_VB, .vb_i = old_p.vb_i                                                         },
                                                    { .flavor = PLAN_RANGE,     .vb_i = old_p.vb2_i, .num_lines = 1, .start_line = old_p.num_lines / 2 - 1 } } ), // each of the two VBs has (old_p.num_lines / 2) lines
                                3, NULL);
                                
                    *B(ReconPlanItem, z_file->recon_plan, plan_i + 3) = 
                                 (ReconPlanItem)    { .flavor = PLAN_END_OF_VB, .vb_i = old_p.vb2_i                                                        };
            }

            break;

        default: ABORT0 ("invalid flavor");
    }
}

static void writer_drop_item_lines_from_end (ReconPlanItemP p, uint32_t drop_lines)
{
    VbInfo *v = VBINFO(p->vb_i);

    switch (p->flavor) {
        case PLAN_RANGE:
            writer_drop_lines (v, p->start_line + p->num_lines - drop_lines, drop_lines);
            break;

        case PLAN_FULL_VB:
            writer_drop_lines (v, v->num_lines - drop_lines, drop_lines);
            break;
        
        case PLAN_INTERLEAVE: {
            VbInfo *v2 = VBINFO(p->vb2_i);
            writer_drop_lines (v, v->num_lines - drop_lines/2, drop_lines/2);           // round down in case of odd number
            writer_drop_lines (v2, v2->num_lines - (drop_lines+1)/2, (drop_lines+1)/2); // round up in case of odd number

            // case: pair-2 is removed in its entirety and pair-1 has only the first line surviving: update plan item to single line
            if (drop_lines + 1 == p->num_lines) {       
                VBINFO (p->vb2_i)->needs_recon = VBINFO (p->vb2_i)->needs_write = false;
                
                *p = (ReconPlanItem){ .flavor = PLAN_RANGE, .vb_i = p->vb_i, .start_line = 0, .num_lines = 1 };

                uint32_t plan_i = BNUM (z_file->recon_plan, p);

                buf_insert (NULL, z_file->recon_plan, ReconPlanItem, plan_i + 1, 
                            (&(ReconPlanItem){ .flavor = PLAN_END_OF_VB, .vb_i = p->vb_i }), 
                            1, NULL);
            }

            break;
        }
        
        default: ABORT ("invalid flavor: %s", display_plan_item (p).s);
    }
}

static void writer_do_tail (int64_t lines_to_trim)
{
    // remove plan items before first_line
    int64_t lines_so_far = 0;
    for (int64_t i=0; i < z_file->recon_plan.len/*modified in PLAN_VB_PLAN*/ && lines_so_far < lines_to_trim; i++) {
        ReconPlanItemP p = B(ReconPlanItem, z_file->recon_plan, i); // note: recon_plan might be realloaced if PLAN_VB_PLAN
        
        if (!p->num_lines/*always for PLAN_TXT_HEADER*/ || !VBINFO(p->vb_i)->needs_write)
            continue; 

        lines_so_far += p->num_lines;
        
        // case: p is fully included in lines to be removed
        if (lines_so_far <= lines_to_trim)  
            writer_drop_plan_item (p, PLAN_REMOVE_ME);
    
        else if (p->flavor == PLAN_VB_PLAN) {
            lines_so_far -= p->num_lines; // undo
            gencomp_piz_expand_PLAN_VB_PLAN (&p, &i, BEFORE_FIRST_EXPANDED_ITEM);
            continue;
        }
        
        // case: some but not all of the lines of p need to be removed:
        // keep the plan item, and use v->is_dropped bitmap to mark some of the lines for dropping
        else 
            writer_drop_item_lines_from_start (p, p->num_lines - (lines_so_far - lines_to_trim));
    }
}

static void writer_do_head (int64_t num_plan_lines, int64_t lines_to_trim)
{
    // remove plan items after last_line
    int64_t i; for (i=z_file->recon_plan.len-1; i >= 0 && lines_to_trim > 0; i--) {
        ReconPlanItemP p = B(ReconPlanItem, z_file->recon_plan, i); // note: recon_plan might be realloaced if PLAN_VB_PLAN
        if (!p->num_lines) continue; 
        
        // case: p is fully included in lines to be removed
        else if (p->num_lines <= lines_to_trim) {
            writer_drop_plan_item (p, PLAN_REMOVE_ME);
            lines_to_trim -= p->num_lines;
        }

        else if (p->flavor == PLAN_VB_PLAN) {
            gencomp_piz_expand_PLAN_VB_PLAN (&p, &i, AFTER_LAST_EXPANDED_ITEM); // p and i now correpond item after last expanded item
            continue;
        }

        // case: p is partially included in lines to be removed
        else {
            writer_drop_item_lines_from_end (p, lines_to_trim);
            lines_to_trim = 0;
        }
    }
}

static void writer_downsample_plan (void)
{
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);

    uint64_t lines_so_far = 0;

    for (uint64_t i=0; i < plan_len; i++) {
        ReconPlanItemP  p = &plan[i];

        uint64_t num_lines_to_next = (flag.downsample - (lines_so_far % flag.downsample) + flag.shard) % flag.downsample;

        // case: we don't need any line from this plan item
        if (p->num_lines && num_lines_to_next >= p->num_lines)
            writer_drop_plan_item (p, PLAN_DOWNSAMPLE); // we count these lines for the purpose of downsampling, but don't write them

        lines_so_far += p->num_lines;
    }
}

// PIZ main thread: after all ordinal filtering (where all lines where available for the ordinal filters), 
// we now remove VBs that are entirely filtered out by --regions, and the reconstructor can further filter
// out lines from the surviving VBs
static void writer_filter_regions (void)
{
    for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) 
        VBINFO(vb_i)->needs_recon &= random_access_is_vb_included (vb_i); // --regions: this VB is excluded
}

// PIZ main thread
static void writer_cleanup_recon_plan_after_genocat_filtering (void)
{
    // VBs that have all their lines dropped - make sure they're marked with !need_recon
    for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {
        VbInfo *v = VBINFO(vb_i);
        if (v->needs_recon && v->is_dropped && bits_is_fully_set (v->is_dropped) &&
            (v->pair_vb_i < vb_i/*not paired, or vb_i is pair-2*/ || bits_is_fully_set (VBINFO(v->pair_vb_i)->is_dropped))) { // if paired, drop a pair_1 only if its pair_2 is dropped as well (for example, in --tail=1 we need the last read of pair_2, so we can't drop pair_1 either) 

            v->needs_recon = false;
            writer_release_vb (v);
        }
    }

    // remove PLAN_REMOVE_ME items or !need_recon from the recon_plan
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);
    uint64_t new_i = 0;


    for (uint64_t old_i=0; old_i < plan_len; old_i++) {
        ReconPlanItemP old_p = &plan[old_i];

        if (!old_p->vb_i ||
            old_p->flavor == PLAN_DOWNSAMPLE || // keep it in plan - just for downsample counting purposes
            (VBINFO(old_p->vb_i)->needs_recon && old_p->flavor != PLAN_REMOVE_ME)) {
            
            if (old_i != new_i) // no need to copy if already in place
                plan[new_i] = *old_p;

            // possibly merge (compact) RANGE with previous line
            if (new_i && 
                plan[new_i].vb_i == plan[new_i-1].vb_i &&
                plan[new_i].flavor == PLAN_RANGE && plan[new_i-1].flavor == PLAN_RANGE &&
                plan[new_i].start_line == plan[new_i-1].start_line + plan[new_i-1].num_lines)

                plan[new_i-1].num_lines += plan[new_i].num_lines; // just add the lines to the previous plan item

            else // no merge - keep added item
                new_i++;
        }
    }

    z_file->recon_plan.len = new_i;
}

// generates z_file->piz_reading_list consisting of the TXT_HEADER and VB_HEADER sections needed for
// reconstruction, in the order of reconstruction
static void writer_create_piz_reading_list (CompIType comp_i)
{
    // alloate in case we need to replace the section_list
    buf_alloc (evb, &z_file->piz_reading_list, z_file->num_vbs + sections_txt_header_get_num_fragments(), 0, SectionEnt, 1, "z_file->piz_reading_list"); 
    z_file->piz_reading_list.next = 0;
    
    // add all TXT_HEADERs and VB_HEADERs according to the order of their first appearance in the recon_plan
    for_buf (ReconPlanItem, pi, z_file->recon_plan)
        if (pi->flavor == PLAN_TXTHEADER)
            sections_reading_list_add_txt_header (pi->comp_i); // add all fragments

        else {
            VbInfo *v = VBINFO(pi->vb_i); 
            if (!v->encountered) { // first encounter with this VB
                sections_reading_list_add_vb_header (v->vblock_i);
                v->encountered = true;
            }

            if (pi->flavor == PLAN_INTERLEAVE) {
                VbInfo *v2 = VBINFO(pi->vb2_i); 
                if (!v2->encountered) { // first encounter with this VB
                    sections_reading_list_add_vb_header (v2->vblock_i);
                    v2->encountered = true;
                }
            }
        }
}

// PIZ main thread: add txtheader entry for the component
static void writer_add_txtheader_plan (CompIType comp_i)
{
    if (!TXTINFO(comp_i)->needs_recon) return; 

    buf_alloc (wvb, &z_file->recon_plan, 1, 1000, ReconPlanItem, 1.5, "recon_plan");

    BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
        .flavor = PLAN_TXTHEADER,
        .comp_i = comp_i
        // num_lines remains 0, and txt_header data is not counted for --head, --tail, --downsample etc
    };

    if (flag.show_recon_plan && z_has_gencomp && comp_i == SAM_COMP_MAIN)
        recon_plan_show (0, z_file->recon_plan.len-1, 1); // recon plan for this TXT_HEADER only
}

// plan for --one-vb of a PRIM or DEPN VB
void writer_add_one_vb_plan_prim_or_depn (VBIType vblock_i)
{
    BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
        .flavor    = PLAN_VB_PLAN, // all lines
        .vb_i      = vblock_i,
        .num_lines = VBINFO(vblock_i)->num_lines
    };
}

// PIZ main thread: add "full vb" entry for each VB of the component
void writer_add_trivial_plan (CompIType comp_i, PlanFlavor flavor)
{
    VBIType num_vbs = sections_get_num_vbs (comp_i);
    if (!num_vbs) return;

    buf_alloc (wvb, &z_file->recon_plan, num_vbs + 2, 1000, ReconPlanItem, 1, "recon_plan"); // +2 for potentially 2 TXT_HEADERs in case paired FASTQ

    Section vb_header_sec = NULL;
    while (sections_get_next_vb_header_sec (comp_i, &vb_header_sec)) { // in order of vb_i of this component
        VbInfo *v = VBINFO(vb_header_sec->vblock_i);

        if (v->needs_recon) 
            BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
                .flavor    = flavor, // all lines
                .vb_i      = vb_header_sec->vblock_i,
                .num_lines = v->num_lines + (flavor == PLAN_VB_PLAN ? (v->num_gc_lines[SAM_COMP_PRIM-1] + v->num_gc_lines[SAM_COMP_DEPN-1]) : 0)
            };
    }
}

// PIZ main thread: add interleave entry for each VB of the component
// also interleaves VBs of the two components and eliminates the second TXT_HEADER entry
static void writer_add_interleaved_plan (CompIType comp_1, CompIType comp_2)
{
    // expecting both componenets to have the same number of VBs
    VBIType R1_num_vbs = sections_get_num_vbs (comp_1);
    VBIType R2_num_vbs = sections_get_num_vbs (comp_2);
    
    ASSERT (R1_num_vbs == R2_num_vbs, "R1 has %u VBs, but R2 has %u VBs - expecting them to be equal when --interleaved",
             R1_num_vbs, R2_num_vbs);

    ASSERT0 (R1_num_vbs, "Component has no VBs");

    writer_add_txtheader_plan (comp_1); // enough to add one txt header - just to open the txt file

    buf_alloc (wvb, &z_file->recon_plan, R1_num_vbs, 0, ReconPlanItem, 0, "recon_plan");

    for (VBIType i=0; i < R1_num_vbs; i++) {

        VBIType vb1_i = i + sections_get_first_vb_i (comp_1);
        VBIType vb2_2 = i + sections_get_first_vb_i (comp_2);
        
        VbInfo *v1 = VBINFO(vb1_i);
        VbInfo *v2 = VBINFO(vb2_2);

        // if one of them is filtered out, both are, and if no reconstructing then obviously no writing either            
        if (!v1->needs_recon || !v2->needs_recon) 
            v1->needs_recon = v2->needs_recon = v1->needs_write = v2->needs_write = false;

        if (v1->needs_write && v2->needs_write)           
            BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
                .flavor    = PLAN_INTERLEAVE, 
                .vb_i      = vb1_i,
                .vb2_i     = vb2_2,
                .num_lines = v1->num_lines * 2
            };
        
        else // if one of them is filtered out, both are
            v1->needs_write = v2->needs_write = false;
    }     
}

static void writer_add_SAM_plan (void)
{
    writer_add_txtheader_plan (SAM_COMP_MAIN);

    if (z_has_gencomp && VER2(15,65))                 
        gencomp_piz_initialize_vb_info();
    
    else if (z_has_gencomp && !VER2(15,65)) 
        recon_plan_add_prescribed_by_recon_plan_section();
    
    else
        writer_add_trivial_plan (SAM_COMP_MAIN, PLAN_FULL_VB);
}

static void noreturn writer_report_count (void)
{
    uint64_t count=0;
    for_buf (ReconPlanItem, p, z_file->recon_plan) 
        if (p->flavor==PLAN_RANGE || p->flavor==PLAN_FULL_VB || p->flavor==PLAN_INTERLEAVE || p->flavor==PLAN_DOWNSAMPLE || p->flavor==PLAN_VB_PLAN)
            if (VBINFO(p->vb_i)->needs_write)
                count += p->num_lines;

    if (flag.downsample) {
        bool rounddown = flag.shard >= (count % flag.downsample);
        count = (count + (flag.downsample-1)) / flag.downsample - rounddown;
    }

    iprintf ("%"PRIu64"\n", count);

    exit_ok;
}

static void writer_apply_genocat_flags_to_recon_plan (void)
{
    bool has_regions_filter = random_access_has_filter();
    
    // filtering
    if (flag.maybe_lines_dropped_by_writer || has_regions_filter) {

        int64_t num_lines = (flag.lines_first != NO_LINE || flag.tail) ? writer_get_plan_num_lines() : 0; // calculate if needed
        int64_t lines_to_trim;

        // first - filters based on line ordinal position (--head, --lines, --tail)
        if ((flag.lines_first != NO_LINE || flag.tail) 
            && (lines_to_trim = flag.tail ? MIN_(num_lines - flag.tail, num_lines) : flag.lines_first)) 
            writer_do_tail (lines_to_trim);

        if (flag.lines_last != NO_LINE 
            && (lines_to_trim = num_lines - MIN_(flag.lines_last+1, num_lines)))
            writer_do_head (num_lines, lines_to_trim);

        // second - filters based on contents of line - remove VBs by random_access here, more in reconstructor
        if (has_regions_filter) 
            writer_filter_regions (); // note: only filters out whole VBs, does not mark lines for dropping

        // if --downsample is the only line dropping filter, we can remove some plan items now 
        // note: this is an optimization, lines which are not removed now, will be discarded by the writer thread,
        // and lines that are moved, will still be counted for downsampling purposes by the writer thread.
        if (flag.downsample && !flag.tail && flag.lines_first == NO_LINE && !flag.maybe_lines_dropped_by_reconstructor)
            writer_downsample_plan(); 
    }

    // mark VBs that are fully dropped as !needs_recon + remove REMOVE_ME items from recon_plan + compact plan
    if (flag.maybe_lines_dropped_by_writer || has_regions_filter || flag.one_vb) 
        writer_cleanup_recon_plan_after_genocat_filtering();

    // if user wants just genocat --count and reconstructor doesn't drop lines - we know the answer from the plan
    if (flags_writer_counts()) 
        writer_report_count(); // doesn't return
}

// PIZ main thread: 
// (1) create reconstruction plan for one txt_file
// (2) create z_file->piz_reading_list consisting of TXT_HEADER, VB_HEADER sections in the order needed for recontructing this txt file
void writer_create_plan (CompIType comp_i)
{  
    START_TIMER;

    // initialize as might be set by previous components. Note: we use z_file and not txt_file, bc txt_file is not allocated yet.    
    z_file->recon_plan.len = z_file->piz_reading_list.len = 0; 

    ASSINP (!flag.interleaved || txt_header_info.len == 2 || flag.deep_fq_only, // in case of Deep, we verify the conditions for interleaved in flags_piz_set_flags_for_deep
            "--interleave cannot be used because %s was not compressed with --pair", z_name);

    // case: Single FASTQ component: --R1, --R2 or --one-vb of a FASTQ VB
    if (flag.one_component && (Z_DT(FASTQ) || (Z_DT(SAM) && flag.one_component-1 >= SAM_COMP_FQ00))) {
        if (Z_DT(SAM)) writer_add_SAM_plan();  // one Deep FASTQ component - need to get SAM first
            
        writer_add_txtheader_plan (flag.one_component-1);
        writer_add_trivial_plan (flag.one_component-1, PLAN_FULL_VB); // one_component is copm_i+1
    }

    // case: interleaved - Pair
    else if ((Z_DT(FASTQ) && flag.interleaved))  
        writer_add_interleaved_plan (FQ_COMP_R1, FQ_COMP_R2); // single txt_file, recon_plan with PLAN_INTERLEAVE to interleave pairs of VBs

    // case: interleaved - Deep
    else if (Z_DT(SAM) && flag.deep && flag.interleaved) { 
        writer_add_SAM_plan();  
        writer_add_interleaved_plan (SAM_COMP_FQ00, SAM_COMP_FQ01); // single txt_file, recon_plan with PLAN_INTERLEAVE to interleave pairs of VBs
    }

    // case: SAM with PRIM/DEPN (i.e. gencomp) - either non-deep or deep and --sam or --bam specified
    else if (Z_DT(SAM) && comp_i == SAM_COMP_MAIN) 
        writer_add_SAM_plan();

    // case: paried FASTQ with --one-vb : include both components - some are dropped in VBINFO
    else if (Z_DT(FASTQ) && flag.pair && flag.one_vb) {
        writer_add_txtheader_plan (FQ_COMP_R1);
        writer_add_trivial_plan (FQ_COMP_R1, PLAN_FULL_VB);
        writer_add_txtheader_plan (FQ_COMP_R2);
        writer_add_trivial_plan (FQ_COMP_R2, PLAN_FULL_VB);
    }

    // normal file, not one of the above
    else {
        writer_add_txtheader_plan (comp_i);
        writer_add_trivial_plan (comp_i, PLAN_FULL_VB);
    }

    // apply genocat flags
    if (is_genocat) 
        writer_apply_genocat_flags_to_recon_plan();

    // create z_file->piz_reading_list recon_plan (only TXT_HEADER and VB_HEADER sections) 
    // note: if we have PRIM/DEPN data in 15.0.64 or later, this will be just the MAIN VBs for now and expanded later in gencomp_piz_update_reading_list
    writer_create_piz_reading_list (comp_i); 

    if (flag.show_reading_list) 
        sections_show_section_list (z_file->data_type, &z_file->piz_reading_list, SEC_NONE);

    // actual number of buffers - the maximum of reported by ZIP, but not less than 3 or more than num_vbs
    z_file->max_conc_writing_vbs = MIN_(z_file->num_vbs, MAX_(3, z_file->max_conc_writing_vbs));

    // show recon plan, except if genocat with gencomp and no head/tail/lines (in which case we've already shown it in )
    // to do: in case of SAM gencomp with --head/--tail/--lines: we show the plan here - with the one VB from which we cut lines expanded, and others still in VB_PLAN_VB format. Ideally we should expand them too.
    if (flag.show_recon_plan) {
        recon_plan_show_all(); // recon_plan for entire file    
        if (is_genocat) exit_ok;
    }
    
#if defined _WIN32 || defined APPLE
    ASSERTW ((uint64_t)z_file->max_conc_writing_vbs * segconf.vb_size < MEMORY_WARNING_THREASHOLD,
             "\nWARNING: This file's output will be re-ordered in-memory, which will consume %"PRIu64" MB.\n"
             "Alternatively, use the --unsorted option to avoid in-memory sorting", (segconf.vb_size * z_file->max_conc_writing_vbs) >> 20);
#endif
    COPY_TIMER_EVB (writer_create_plan);
}

// -------------------
// Writer thread stuff
// -------------------

#define BGZF_FLUSH_THRESHOLD  (32 MB) 
#define PLAIN_FLUSH_THRESHOLD (4 MB) 
#define FLUSH_THRESHOLD (TXT_IS_BGZF ? BGZF_FLUSH_THRESHOLD : PLAIN_FLUSH_THRESHOLD)

static void writer_write (BufferP buf, uint64_t txt_data_len)
{
    START_TIMER;

    if (!buf->len) return;
    
    txtfile_fwrite (STRb(*buf));
    
    txt_file->disk_so_far += buf->len;

    buf_free (*buf);

    COPY_TIMER_FULL (wvb, write, false);
}

uint32_t writer_get_max_bgzf_threads (void)
{
    return MIN_(global_max_threads, 10); // note: we don't need more threads than the ratio between BGZF compression speed and disk write speed
}

// if bgzf pool is full, wait for one thread to complete
static bool writer_output_one_processed_bgzf (Dispatcher dispatcher, bool blocking)
{
    VBlockP bgzf_vb = dispatcher_get_processed_vb (dispatcher, NULL, blocking);

    if (bgzf_vb) {
        writer_write (&bgzf_vb->comp_txt_data, bgzf_vb->txt_data.len);
        dispatcher_recycle_vbs (dispatcher, true);  // also release VB
    }

    return !!bgzf_vb;
}

static void writer_flush_vb (Dispatcher dispatcher, VBlockP vb, bool is_txt_header, bool is_last)
{
    ASSERTNOTNULL(vb);

    if (!Ltxt && !is_last) return; // no data to flush

    txt_file->txt_data_so_far_single += vb->txt_data.len;

    if (!flag.no_writer) { // note: we might have a writer thread despite no writing - eg for calculated the digest if we have SAM gencomp
        // case: BGZF compression - offload compression to compute threads (a max of BGZF_FLUSH_THRESHOLD per thread)
        if (dispatcher) {
            // in case of wvb, write the entire buffer, if its "native" vb, chop it up
            uint32_t chunk_size = (vb==wvb) ? Ltxt : BGZF_FLUSH_THRESHOLD;

            for (uint32_t i=0; i < Ltxt || is_last; i += chunk_size) { // if txt_data is empty, we still do one iteration in case of is_last

                if (vb_pool_is_full (POOL_BGZF))
                    writer_output_one_processed_bgzf (dispatcher, true); 

                bgzf_dispatch_compress (dispatcher, Btxt (i), MIN_(chunk_size, Ltxt - i), vb->comp_i,
                                        is_last && (i+chunk_size >= Ltxt));
            
                is_last = false; // if is is_last, we enter the loop once if txt_data.len=0
            }

            buf_free (vb->txt_data);
        }

        // case: no BGZF, so we can write to disk right away
        else if (Ltxt) 
            writer_write (&vb->txt_data, Ltxt);
    }

    // case this is wvb: free, so we can use it for more lines (for non-wvb, we will release the VB after verifying digest)
    if (vb == wvb)
        buf_free (vb->txt_data);
}

// final filter of output lines, based on downsample and shard. returns true if line is to be output
static inline bool writer_line_survived_downsampling (VbInfo *v)
{
    if (!flag.downsample) return true;

    uint64_t line_i = txt_file->lines_written_so_far / ((flag.interleaved || flag.sequential) ? 2ULL : 1ULL); // 2 lines at a time if interleave (FASTQ) or sequential (FASTA)

    // show (or not) the line based on our downsampling rate and shard value
    return (line_i % (uint64_t)flag.downsample) == flag.shard;
}

// write lines one at a time, watching for dropped lines
static void writer_write_line_range (VbInfo *v, uint32_t start_line, uint32_t num_lines)
{
    ASSERTNOTNULL (v);
    ASSERTNOTNULL (v->vb);
    
    ASSERT (!v->vb->scratch.len, "expecting vb=%s/%u data to be BGZF-compressed by writer at flush, but it is already compressed by reconstructor: txt_file->effective_codec=%s compressed.len=%"PRIu64" txt_data.len=%"PRIu64,
            comp_name (v->vb->comp_i), v->vb->vblock_i, codec_name (txt_file->effective_codec), v->vb->scratch.len, v->vb->txt_data.len);

    if (!v->vb->txt_data.len) return; // no data in the VB 

    ARRAY (uint32_t, lines, v->vb->lines); 

    ASSERT (start_line + num_lines <= lines_len, "vb_i=%u invalid range: least start_line=%u + num_lines=%u == %u but lines_len=%"PRIu64, 
            VBINFO_NUM(v), start_line, num_lines, start_line + num_lines, lines_len);

    for (uint32_t line_i=start_line; line_i < start_line + num_lines; line_i++) {
        
        rom start = Bc (v->vb->txt_data, lines[line_i]);
        rom after = Bc (v->vb->txt_data, lines[line_i+1]);  // note: lines has one extra entry so this is always correct
        uint32_t line_len = (uint32_t)(after - start);
        
        bool is_dropped = v->is_dropped && bits_get (v->is_dropped, line_i);

        ASSERT (after >= start, "vb_i=%u Writing line %i (start_line=%u num_lines=%u): expecting start=%p <= after=%p", 
                VBINFO_NUM(v), line_i, start_line, num_lines, start, after);
 
        if (!is_dropped) { // don't output lines dropped in container_  reconstruct_do due to vb->drop_curr_line
            
            if (writer_line_survived_downsampling(v))
                buf_add_more (wvb, &wvb->txt_data, start, line_len, "txt_data");
            
            txt_file->lines_written_so_far++; // increment even if downsampled-out, but not if filtered out during reconstruction (for downsampling accounting)
            v->vb->num_nondrop_lines--;       // update lines remaining to be written (initialized during reconstruction in container_reconstruct)
        }
    }
}

static void writer_write_lines_add_pair (VbInfo *v, rom start, unsigned len, unsigned pair /* 1 or 2 */)
{
    static rom suffixes[3] = { "", "/1", "/2" }; // suffixes for pair 1 and pair 2 reads

    SAFE_NUL (&start[len]);
    unsigned qname_len = strcspn (start, " \n\t");
    SAFE_RESTORE;

    rom sep = &start[qname_len];
    ASSERT (qname_len < len, "Expected to find a newline in the 4-line read:\n%.*s\n", len, start);

    // write up to the separator
    buf_add_more (wvb, &wvb->txt_data, start, qname_len, "txt_data");

    // write suffix if requested, and suffix is not already present
    if (qname_len < 3 || sep[-2] != '/' || sep[-1] != '0' + pair)
        buf_add_more (wvb, &wvb->txt_data, suffixes[pair], 2, "txt_data");

    buf_add_more (wvb, &wvb->txt_data, sep, len - qname_len, "txt_data");
}

// write one fastq "line" (actually 4 textual lines) at time, interleaved from v1 and v2
static void writer_write_lines_interleaves (Dispatcher dispatcher, VbInfo *v1, VbInfo *v2)
{
    if (!v1->vb->txt_data.len || !v2->vb->txt_data.len) return; // no data in the either of the VBs 

    ASSERT (v1->num_lines == v2->num_lines, "when interleaving, expecting number of lines of vb_i=%u (%u) to be the same as vb_i=%u (%u)",
            v1->vb->vblock_i, v1->num_lines, v2->vb->vblock_i, v2->num_lines);
    
    ARRAY (uint32_t, lines1, v1->vb->lines);
    ARRAY (uint32_t, lines2, v2->vb->lines);
    
    for (uint32_t line_i=0; line_i < v1->num_lines; line_i++) {

        rom start1 = Bc (v1->vb->txt_data, lines1[line_i]);
        rom start2 = Bc (v2->vb->txt_data, lines2[line_i]);
        unsigned len1 = lines1[line_i+1] - lines1[line_i];
        unsigned len2 = lines2[line_i+1] - lines2[line_i];
        
        bool is_dropped1 = v1->is_dropped && bits_get (v1->is_dropped, line_i);
        bool is_dropped2 = v2->is_dropped && bits_get (v2->is_dropped, line_i);

        // skip lines dropped in container_reconstruct due to vb->drop_curr_line for either leaf
        if ((flag.interleaved == INTERLEAVE_BOTH && !is_dropped1 && !is_dropped2) || 
            (flag.interleaved == INTERLEAVE_EITHER && (!is_dropped1 || !is_dropped2))) {
            // TODO - INTERLEAVE_EITHER doesn't work yet, bug 396
            if (writer_line_survived_downsampling(v1)) { 
                if (len1) writer_write_lines_add_pair (v1, start1, len1, 1); // also adds /1 and /2 if paired
                if (len2) writer_write_lines_add_pair (v2, start2, len2, 2);
            }

            txt_file->lines_written_so_far += 2; // increment even if downsampled-out, but not if filtered out during reconstruction (for downsampling accounting)
        }

        if (wvb->txt_data.len > FLUSH_THRESHOLD) 
            writer_flush_vb (dispatcher, wvb, false, false);
    }
}

// writer thread: waiting for data from a VB and loading it
static void writer_load_vb (Dispatcher dispatcher, VbInfo *v, bool is_txtheader)
{
    if (flag_show_threads && !is_txtheader) 
        iprintf ("writer: vb=%s/%u WAITING FOR VB\n", comp_name(v->comp_i), VBINFO_NUM(v));

    if (flag_is_show_vblocks (PIZ_TASK_NAME) && !is_txtheader) 
        iprintf ("WRITER_WAITING_FOR vb=%s/%u\n", comp_name(v->comp_i), VBINFO_NUM(v));

    else if (flag_show_threads && is_txtheader) 
        iprintf ("writer: comp=%s WAITING FOR TXT_HEADER\n", comp_name(v->comp_i));

    bool got_data;
    do {
        bool computing_bgzf = dispatcher && !vb_pool_is_empty (POOL_BGZF);

        got_data = mutex_wait (v->wait_for_data, !computing_bgzf); // blocking if there is no BGZF currently processing 

        // if data is ready yet, use this waiting time to write processed bgzf to disk
        if (!got_data) {
            // output to txt_file all BGZF vbs that completed compression
            bool any_data_flushed = false;
            while (writer_output_one_processed_bgzf (dispatcher, false)) // non-blocking
                any_data_flushed = true;

            if (!any_data_flushed) 
                usleep (15000); // neither the VB, or processed data is ready - take a nap
        } 

    } while (!got_data);

    threads_log_by_vb (v->vb, "writer", v->vb->vblock_i ? "WRITER LOADED VB" : "WRITER LOADED TXT_HEADER", 0);

    if (flag_is_show_vblocks (PIZ_TASK_NAME) && !is_txtheader) 
        iprintf ("WRITER_LOADED vb=%s/%u\n", comp_name(v->comp_i), VBINFO_NUM(v));

    v->is_loaded = true;

    ASSERT (v->vblock_i == v->vb->vblock_i, "Wrong VB: v->vblock_i=%u but v->vb->vblock_i=%u", v->vblock_i, v->vb->vblock_i);
}

// Thread entry point for writer thread - at this point, the reconstruction plan is ready and unmutable
static void writer_main_loop (VBlockP wvb) // same as wvb global variable
{
    START_TIMER;

    ASSERTNOTNULL (wvb);

    threads_set_writer_thread();

    // if we need to BGZF-compress, we will dispatch the compression workload to compute threads
    Dispatcher dispatcher = (!flag.no_writer && TXT_IS_BGZF && txt_file->mgzip_flags.library != BGZF_EXTERNAL_LIB) ? 
        dispatcher_init ("bgzf", NULL, POOL_BGZF, writer_get_max_bgzf_threads(), 0, false, false, NULL, 0, NULL) : NULL;

    // normally, we digest in the compute thread but in case gencomp lines can be inserted into the vb we digest here.
    bool do_digest = piz_need_digest && z_has_gencomp;

    // execute reconstruction plan
    for (int64_t i=0; i < z_file->recon_plan.len; i++) { // note: recon_plan.len maybe 0 if everything is filtered out
        ReconPlanItemP p = B(ReconPlanItem, z_file->recon_plan, i);

        VbInfo *v  = p->vb_i ? VBINFO(p->vb_i) 
                   : p->flavor == PLAN_TXTHEADER ? TXTINFO(p->comp_i)
                   : NULL;

        VbInfo *v2 = (p->vb_i && p->flavor == PLAN_INTERLEAVE) ? VBINFO(p->vb2_i) : NULL;

        if ((!v || !v->needs_write) && (!v2 || !v2->needs_write)) continue;

        // if data for this VB is not ready yet, wait for it
        if (v && v->needs_write && !v->is_loaded && p->flavor != PLAN_DOWNSAMPLE)       
            writer_load_vb (dispatcher, v, p->flavor == PLAN_TXTHEADER);

        if (v2 && v2->needs_write && !v2->is_loaded) 
            writer_load_vb (dispatcher, v2, false);

        // Note: A Deep FASTQ, even in a gencomp SAM file, is digested in the compute thread.
        bool do_digest_v = do_digest && v->vb->data_type != DT_FASTQ;

        switch (p->flavor) {

            case PLAN_TXTHEADER:   
                writer_flush_vb (dispatcher, v->vb, true, false); // write the txt header in its entirety
                writer_release_vb (v);
                break;

            case PLAN_VB_PLAN: // gencomp VB starting 15.0.64
                gencomp_piz_vb_to_plan (v->vb, &i, true);
                continue; // items already executed are removed and VB is inserted at beginging of modified plan ; i reset to 0
            
            case PLAN_FULL_VB:   
                // note: normally digest calculation is done in the compute thread in piz_reconstruct_one_vb, but in
                // case of an unmodified VB that inserts lines from gencomp VBs, we do it here, as we re-assemble the original VB
                if (do_digest_v) digest_one_vb (v->vb, false, NULL); 

                if (!flag.downsample) {

                    if (flag_is_show_vblocks (PIZ_TASK_NAME)) // only displayed for entire VBs, not line ranges etc 
                        iprintf ("VB_FLUSH_FULL_VB(id=%d) vb=%s/%d txt_data.len=%u\n", v->vb->id, comp_name (v->vb->comp_i), v->vb->vblock_i, v->vb->txt_data.len32);

                    writer_flush_vb (dispatcher, wvb, false, false); // flush any remaining unflushed wvb lines from previous VBs
                    writer_flush_vb (dispatcher, v->vb, false, false); // write entire VB
                    txt_file->lines_written_so_far += v->vb->num_nondrop_lines;
                }
                else 
                    writer_write_line_range (v, 0, v->num_lines); // this skips downsampled-out lines

                dispatcher_increment_progress ("txt_write", 1); // done writing VB
                writer_release_vb (v);
                break;

            case PLAN_DOWNSAMPLE:
                txt_file->lines_written_so_far += p->num_lines;
                break;

            case PLAN_END_OF_VB: { // done with VB - free the memory (happens after a series of "default" line range entries)
                v->vb->lines.len32 = wvb->lines.len32; // for correct digest error printing
                bool needs_flush = do_digest_v ? digest_one_vb (v->vb, false, &wvb->txt_data) : true;

                // note: if we're digesting gencomp VBs, we dont flush wvb when finishing a non-digestable VB (eg PRIM/DEPN in SAM/BAM)
                if (needs_flush || i == z_file->recon_plan.len-1/*final VB*/)  
                    writer_flush_vb (dispatcher, wvb, false, false); // flush remaining unflushed lines of this VB

                writer_release_vb (v);
                wvb->lines.len32 = 0;
                
                v->range_vb_in_use = false; 

                dispatcher_increment_progress ("txt_write", 1); // done writing VB
                break;
            }
            case PLAN_INTERLEAVE:
                writer_write_lines_interleaves (dispatcher, v, v2);
                writer_release_vb (v);
                writer_release_vb (v2);
                dispatcher_increment_progress ("txt_writeX2", 2); // done writing both VBs
                break;

            default:   // PLAN_RANGE
                wvb->lines.len32 += p->num_lines;
                writer_write_line_range (v, p->start_line, p->num_lines);

                v->range_vb_in_use = true; 
                break;
        }

        if (!do_digest_v && wvb->txt_data.len > FLUSH_THRESHOLD)
            writer_flush_vb (dispatcher, wvb, false, false);
    }

    ASSERT (!do_digest || !wvb->txt_data.len, "expecting wvb->txt_data.len=%u to be 0. Perhaps a dropped END_OF_VB?", wvb->txt_data.len32); // sanity

    // end a VB if PLAN_END_OF_VB went missing. This happens when a PLAN_VB_PLAN containing the PLAN_END_OF_VB is eliminated in its entirety (e.g. due to --head),
    // but an earlier PLAN_VB_PLAN contains PLAN_RANGE of that PRIM/DEPN VB.
    for_buf (VbInfo, v, z_file->vb_info[0])
        if (v->range_vb_in_use) {
            ASSERT (is_genocat && z_has_gencomp && (v->comp_i == SAM_COMP_DEPN || v->comp_i == SAM_COMP_PRIM),
                    "PLAN_END_OF_VB is missing for %s/%u, but this is only possible for SAM DEPN or PRIM vbs.", comp_name (v->comp_i), v->vblock_i);

            writer_release_vb (v);
            dispatcher_increment_progress ("txt_write", 1); // done writing VB
        }
            
    // this might have data (eg with flag.downsample or flag.interleave) or not
    writer_flush_vb (dispatcher, wvb, false, true);

    // output for all remain BGZF (if any)
    if (dispatcher) {
        while (!vb_pool_is_empty (POOL_BGZF)) 
            writer_output_one_processed_bgzf (dispatcher, true);

        bgzf_write_finalize(); // write final data to wvb->comp_txt_data
        writer_write (&wvb->comp_txt_data, 0);

        dispatcher_finish (&dispatcher, NULL, false, false);
    }
  
    // The recon plan is supposed to end with all VBs flushed due to END_OF_VB, FULL_VB etc
    ASSERT (!wvb->txt_data.len, "Expected wvb to be flushed, but wvb->txt_data.len=%"PRIu64, wvb->txt_data.len);

    threads_unset_writer_thread();

    COPY_TIMER_FULL (wvb, writer_main_loop, false);
}

//---------------------------------------------------
// PIZ interaction with writer: start, stop, handover
//---------------------------------------------------

// PIZ main thread: launch writer thread to write to a single txt_file 
static void writer_start_writing (CompIType unbind_comp_i)
{
    ASSERTNOTNULL (txt_file);

    if (writer_thread != THREAD_ID_NONE || flag.no_writer_thread) return;

    writer_thread = threads_create (writer_main_loop, wvb);
}

// PIZ main thread: wait for writer thread to finish writing a single txt_file
void writer_finish_writing (bool is_last_txt_file)
{
    if (writer_thread != THREAD_ID_NONE) 
        // wait for thread to complete (possibly it completed already)
        threads_join (&writer_thread, WRITER_TASK_NAME); // also sets writer_thread=THREAD_ID_NONE
    
    // all mutexes destroyed by main thread, that created them (because we will proceed to freeing z_file in which they live)    
    if (is_last_txt_file) {
        for_buf (VbInfo, comp, txt_header_info)
            mutex_destroy (comp->wait_for_data);

        for_buf (VbInfo, v, z_file->vb_info[0])
            mutex_destroy (v->wait_for_data);

        // destroy all is_dropped bufs allocated by writer in evb (more efficient multiple destroys)
        buflist_destroy_bufs_by_name (IS_DROPPED_BUF_NAME, true); 
    }
}

// PIZ main thread. returns true if indeed handed over
static bool writer_handover (VbInfo *v, VBlockP vb)
{
    if (flag.no_writer_thread || !v->needs_write) return false;

    writer_start_writing (flag.unbind ? v->comp_i : COMP_NONE); // start writer thread if not already started
    
    // writer thread now owns this VB, and is the only thread that will modify it, until finally destroying it
    v->vb = vb;

    threads_log_by_vb (vb, vb->compute_task, vb->vblock_i ? "HANDING OVER VB" : "HANDING OVER TXT_HEADER", 0);

    mutex_unlock (v->wait_for_data); 

    if (flag_is_show_vblocks (PIZ_TASK_NAME)) 
        iprintf ("HANDED_OVER(task=%s id=%d) vb_i=%s/%u txt_data.len=%u\n", PIZ_TASK_NAME, vb->id, comp_name(v->comp_i), vb->vblock_i, Ltxt);

    return true;
}

// PIZ main thread: hand over data from a txtheader whose reconstruction main thread has completed, to the writer thread
bool writer_handover_txtheader (VBlockP *txt_header_vb_p)
{
    bool is_handed_over = writer_does_txtheader_need_write ((*txt_header_vb_p)->comp_i)
                       && writer_handover (TXTINFO((*txt_header_vb_p)->comp_i), *txt_header_vb_p);
    
    if (is_handed_over) *txt_header_vb_p = NULL;

    return is_handed_over;
}

// PIZ main thread: hand over data from a VB whose reconstruction compute thread has completed, to the writer thread
bool writer_handover_data (VBlockP *vb_p)
{
    bool is_handed_over = writer_handover (VBINFO((*vb_p)->vblock_i), *vb_p);
    
    if (is_handed_over) *vb_p = NULL;

    return is_handed_over;
}

