// ------------------------------------------------------------------
//   writer.c
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "context.h"
#include "profiler.h"
#include "strings.h"
#include "file.h"
#include "threads.h"
#include "sections.h"
#include "random_access.h"
#include "writer.h"
#include "zfile.h"
#include "bgzf.h"
#include "mutex.h"
#include "codec.h" 
#include "endianness.h"
#include "bits.h"
#include "version.h"
#include "gencomp.h"
#include "piz.h"
#include "recon_plan_io.h"
#include "dispatcher.h"
#include "progress.h"

// ---------------
// Data structures
// ---------------

#define WRITER_TASK_NAME "writer"

typedef struct {
    Mutex wait_for_data;   // initialized locked, unlocked when data is moved to txt_data and is_loaded is set
    VBlockP vb;            // data handed over main -> compute -> main -> writer
    VBIType vblock_i;      // this is the same as the index of this entry, and also of vb->vblock_i (after vb is loaded). for convenience.
    VBIType pair_vb_i;     // vb_i of pair vb, or 0 if there isn't any. if pair_vb_i > vblock_i then we are pair_1
    uint32_t num_lines;    // number of lines in this VB - readed from VB header only if needed
    CompIType comp_i;      // comp_i of this VB (0-3)
    bool is_loaded;        // has data moving to txt_data completed
    bool needs_recon;      // VB/TXT_HEADER should be reconstructed, but writing it depends on needs_write
    bool needs_write;      // VB/TXT_HEADER is written to output file, and hence is in recon_plan
    bool no_ordinal_filter;// Never filter out lines from this VB due to head/tail/lines/downsample
    bool full_vb;
    bool encountered;      // used by writer_update_section_list
    bool fasta_contig_grepped_out; // FASTA: a VB that starts with SEQ data, inherits this from the last contig of the previous VB 
    Buffer is_dropped_buf; // a bitarray with a bit set is the line is marked for dropping
    BitsP is_dropped;  // pointer into is_dropped_buf
} VbInfo; // used both for VBs and txt headers

#define VBINFO(vb_i) B(VbInfo, vb_info, (vb_i))

#define vb_info         z_file->vb_info[0] // we use only [0] on the PIZ side
#define txt_header_info z_file->txt_header_info

VBlockP wvb = NULL;

static ThreadId writer_thread = THREAD_ID_NONE;

// --------------------------------
// VBlock and Txt Header properties
// --------------------------------

// Note: txtheader lines are used only for error reporting, using writer_get_txt_line_i. They are not
// included in recon_plan, as they are not counted for --head, --downsample etc
void writer_set_num_txtheader_lines (CompIType comp_i, uint32_t num_txtheader_lines)
{
    B(VbInfo, txt_header_info, comp_i)->num_lines = num_txtheader_lines;
}

// calculate 1-based textual line in txt_file (including txt_header), based on vb_i, vb->line_i and recon_plan
// or record in binary file (eg alignment in BAM). In case of a text file, it also includes the header lines.
// For now: returns 0 if reconstruction modifies the file, and therefore recon_plan doesn't reconstruct the original file
// To do: return line according to recon plan even if modified. challenge: lines dropped by reconstructor and not known yet to writer
uint64_t writer_get_txt_line_i (VBlockP vb)
{
    if (flag.data_modified) return 0;

    // since data is unmodified, we have only FULL_VB, RANGE, END_OF_VB and TXTHEADER plan items 
    int64_t txt_num_lines  = 0;
    LineIType vb_num_lines = 0;
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);
    
    for (uint64_t i=0; i < plan_len; i++) {
        if (plan[i].vb_i == vb->vblock_i) {
            // case: plan item containing this line
            if (vb_num_lines + plan[i].num_lines > vb->line_i) 
                return txt_num_lines + (vb->line_i - vb_num_lines) + 1; // +1 because 1-based 

            // case: plan item is from current VB, but we have not yet reached the current line
            else
                vb_num_lines += plan[i].num_lines;
        }

        // note: txtheader lines are not included in the plan item, bc they are not counted for --header, --downsample etc
        txt_num_lines += (plan[i].flavor == PLAN_TXTHEADER) ? B(VbInfo, txt_header_info, plan[i].comp_i)->num_lines
                                                            : VB_DT(DT_FASTQ) ? (4 * plan[i].num_lines) : plan[i].num_lines; // in FASTQ, each plan line is 4 txt lines
        //xxx plan_item_lines = VB_DT(DT_FASTQ) ? (4 * plan[i].num_lines) : plan[i].num_lines; // in FASTQ, each plan line is 4 txt lines
        
        // if (plan[i].vb_i == vb->vblock_i) {
        //     // case: plan item containing this line
        //     if (vb_num_lines + plan_item_lines > vb->line_i) 
        //         return txt_num_lines + (vb->line_i - vb_num_lines) + 1; // +1 because 1-based 

        //     // case: plan item is from current VB, but we have not yet reached the current line
        //     else
        //         vb_num_lines += plan_item_lines;
        // }

        // // note: txtheader lines are not included in the plan item, bc they are not counted for --header, --downsample etc
        // txt_num_lines += (plan[i].flavor == PLAN_TXTHEADER) ? B(VbInfo, txt_header_info, plan[i].comp_i)->num_lines
        //                                                     : plan_item_lines;
    }

    WARN_ONCE0 ("Unexpectedly, unable to find current vb/line_i in recon_plan");
    return 0;
}

// returns my own number in the pair (1 or 2) and pair_vb_i. return 0 if file is not paired.
unsigned writer_get_pair (VBIType vb_i, uint32_t *pair_vb_i)
{
    if (!z_file->z_flags.dts_paired) return 0; // not paired

    VbInfo *v = VBINFO(vb_i);

    *pair_vb_i = v->pair_vb_i;
    return *pair_vb_i > vb_i ? 1 : 2;
}

bool writer_does_txtheader_need_write (Section sec)
{
    ASSERT (sec->st == SEC_TXT_HEADER, "sec->st=%s is not SEC_TXT_HEADER", st_name (sec->st));
    ASSERT (sec->comp_i < txt_header_info.len, "sec->comp_i=%u out of range [0,%d]", sec->comp_i, (int)txt_header_info.len-1);
    
    bool needs_write = B(VbInfo, txt_header_info, sec->comp_i)->needs_write;
    return needs_write;                             
}

bool writer_does_txtheader_need_recon (Section sec)
{
    ASSERT (sec->st == SEC_TXT_HEADER, "sec->st=%s is not SEC_TXT_HEADER", st_name (sec->st));
    ASSERT (sec->comp_i < txt_header_info.len, "comp_i=%u out of range [0,%d]", sec->comp_i, (int)txt_header_info.len-1);

    bool needs_recon = B(VbInfo, txt_header_info, sec->comp_i)->needs_recon;
    return needs_recon;                             
}

bool writer_get_fasta_contig_grepped_out (VBIType vb_i)
{
    ASSERT (vb_i >= 1 && vb_i < vb_info.len, "vb_i=%u out of range [1,%d]", vb_i, (int)vb_info.len-1);

    return VBINFO(vb_i)->fasta_contig_grepped_out;
}

void writer_set_fasta_contig_grepped_out (VBIType vb_i)
{
    ASSERT (vb_i >= 1 && vb_i < vb_info.len, "vb_i=%u out of range [1,%d]", vb_i, (int)vb_info.len-1);

    VBINFO(vb_i)->fasta_contig_grepped_out = true;
}

bool writer_does_vb_need_recon (VBIType vb_i)
{
    ASSERT (vb_i >= 1 && vb_i < vb_info.len, "vb_i=%u out of range [1,%d]", vb_i, (int)vb_info.len-1);

    bool needs_recon = !vb_info.len || VBINFO(vb_i)->needs_recon;
    return needs_recon;
}

// called by: 1. PIZ main thread, when writer creates plan and filters by lines/head/tail
// 2. PIZ compute thread
Bits *writer_get_is_dropped (VBIType vb_i)
{
    ASSERT (vb_i >= 1 && vb_i < vb_info.len, "vb_i=%u out of range [1,%d]", vb_i, (int)vb_info.len-1);
    VbInfo *v = VBINFO(vb_i);

    // allocate if needed. buffer was put in evb buffer_list by writer_init_vb_info
    if (!v->is_dropped) {
        buf_alloc_bits (evb, &v->is_dropped_buf, v->num_lines, "is_dropped");
        buf_zero (&v->is_dropped_buf);
        v->is_dropped = (BitsP)&v->is_dropped_buf;
    }

    return v->is_dropped;
}

// -----------------------------
// PIZ reconstruction plan stuff
// -----------------------------

// PIZ main thread
static VBIType writer_init_txt_header_info (void)
{
    buf_alloc_exact_zero (evb, txt_header_info, sections_get_num_comps(), VbInfo, "txt_header_info"); // info on txt_headers

    VBIType num_vbs = 0;
    Section sec = NULL;

    for (CompIType comp_i=0; comp_i < txt_header_info.len; comp_i++) {

        VbInfo *comp = B(VbInfo, txt_header_info, comp_i);
        comp->comp_i = comp_i;

        num_vbs += sections_get_num_vbs (comp_i);

        sec = sections_get_comp_txt_header_sec (comp_i);
        if (!sec) continue; // no SEC_TXT_HEADER section for this component (eg PRIM/DEPN components of SAM/BAM)

        // conditions entire txt header should be read 
        comp->needs_recon = 
            (  !Z_DT(DT_VCF) // This clause only limits VCF files  
            || ( flag.luft && (comp_i == VCF_COMP_MAIN || (comp_i == VCF_COMP_PRIM_ONLY && !flag.header_one))) 
            || (!flag.luft && (comp_i == VCF_COMP_MAIN || (comp_i == VCF_COMP_LUFT_ONLY && !flag.header_one)))); // --header-one - we don't need the ##primary_only / ##luft_only lines 
            
        // conditions we write the txt header (doesn't affect the VBs of this component)
        comp->needs_write =
           !flag_loading_auxiliary
        &&
           comp->needs_recon
        &&
           !flag.no_header
        &&
           (  !Z_DT(DT_VCF) // This clause only limits VCF files 
           || ( flag.luft && comp_i == VCF_COMP_PRIM_ONLY)  // if luft rendition, show ##primary_only rejects components (appears in the vcf header)
           || (!flag.luft && comp_i == VCF_COMP_LUFT_ONLY)  // if primary rendtion, show ##luft_only rejects components (appears in the vcf header)
           || comp_i == VCF_COMP_MAIN)
        &&
           (  !Z_DT(DT_SAM) // This clause only limits SAM/BAM
           || comp_i == SAM_COMP_MAIN);               // Show only the txt header of the MAIN component

        if (comp->needs_write) {
            // mutex: locked:    here (at initialization)
            //        waited on: writer thread, wanting the data
            //        unlocked:  by main thread after txtheader data is handed over.
            mutex_initialize (comp->wait_for_data);
            mutex_lock (comp->wait_for_data);
        } 
    }
    
    return num_vbs;
}

// PIZ main thread - initialize vb, component, txt_file info. this is run once per z_file
static void writer_init_vb_info (void)
{
    if (txt_header_info.len) return; // already initialized

    vb_info.len = writer_init_txt_header_info() + 1; // +1 as first vb_i=1 (entry 0 will be unused) so we have num_vb+1 entries

    buf_alloc_zero (evb, &vb_info, 0, vb_info.len, VbInfo, 1, "z_file->vb_info");

    uint32_t num_vbs_R=0;
    if (Z_DT(DT_FASTQ) && z_file->z_flags.dts_paired) {
        num_vbs_R = sections_get_num_vbs (FQ_COMP_R1);
        uint32_t num_vbs_R2 = sections_get_num_vbs (FQ_COMP_R2);
        ASSERT (num_vbs_R == num_vbs_R2 || (!VER(10) && !num_vbs_R2), // dts_paired introduced 9.0.13
                "R1 and R2 have a different number of VBs (%u vs %u respecitvely)", num_vbs_R, num_vbs_R2);
    
        // v8-9: case not paired after all
        if (!num_vbs_R2) {
            num_vbs_R=0;  
            z_file->z_flags.dts_paired = false;
        }
    } 

    for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {

        Section sec = sections_vb_header (vb_i, false);
            
        VbInfo *v = VBINFO(vb_i); 
        v->vblock_i  = vb_i;
        v->comp_i    = sec->comp_i;

        // since v14, num_lines is carried by Section. 
        // Note: in DVCF, this is different than the number of lines in the default reconstruction which drops luft_only lines.
        if (VER(14))
            v->num_lines = sec->num_lines;

        // up until v13, num_lines was in SectionHeaderVbHeader
        else {  
            SectionHeaderVbHeader header = zfile_read_section_header (evb, sec->offset, vb_i, sec->st).vb_header;
            v->num_lines = BGEN32 (header.v13_top_level_repeats);
        }

        // we add is_dropped_buf to the evb buffer list. allocation will occur in main thread whenw writer 
        // create the plan, for VBs on the boundary of head/tail/lines, otherwise in reconstruction compute thread.
        buf_add_to_buffer_list (evb, &v->is_dropped_buf);

        // head/tail/lines/downsample don't apply to DVCF reject components which are part of the VCF header
        if (z_is_dvcf && v->comp_i != VCF_COMP_MAIN)
            v->no_ordinal_filter = true;

        // set pairs (used even if not interleaving)
        if (num_vbs_R) {
            v->pair_vb_i = vb_i + num_vbs_R * (v->comp_i % 2 ? -1 : 1);
        
            // verify that this VB and its pair have the same number of lines (test when initiazing the second one)
            uint32_t pair_num_lines = VBINFO(v->pair_vb_i)->num_lines;
            ASSERT (v->pair_vb_i > vb_i || v->num_lines == pair_num_lines, "vb=%u has %u lines, but vb=%u, its pair, has %u lines", 
                    vb_i, v->num_lines, v->pair_vb_i, pair_num_lines);
        }

        // conditions this VB should not be read or reconstructed 
        v->needs_recon = true; // default
        #define DROP v->needs_recon = false

        if (!v->num_lines && // can happen if all lines where sent to gencomp
            !(piz_need_digest && z_has_gencomp)) DROP; // keep anyway if we need to digest

        // --one-vb: user only wants to see a single VB, and this is not it
        else if (flag.one_vb && flag.one_vb != vb_i) DROP; 

        // --header-only: drop all VBs (except VCF - handled separately below)
        else if (flag.header_only && z_file->data_type != DT_VCF) DROP; 

        else switch (z_file->data_type) {

            case DT_FASTQ:
                // --R1 or --R2: we don't show the other component
                if (flag.one_component && flag.one_component-1 != v->comp_i) DROP;

                break;

            case DT_VCF:
                // --header-one or --no-header: we don't show ##primary_only/##luft_only components as they are part of theader            
                if ((flag.header_one || flag.no_header) && v->comp_i != VCF_COMP_MAIN) DROP;

                // in DVCF primary rendition, we don't need ##primary_only lines
                else if (!flag.luft && v->comp_i == VCF_COMP_PRIM_ONLY) DROP;

                // in DVCF luft rendition, we don't need ##luft_only lines
                else if (flag.luft && v->comp_i == VCF_COMP_LUFT_ONLY) DROP;

                // --single-coord:  we don't need the ##primary_only / ##luft_only lines (we do need the TXT_HEADER tough as it contains the ##fileformat line)
                else if (flag.single_coord && v->comp_i != VCF_COMP_MAIN) DROP;

                // --header-only VCF - drop MAIN comp VBs (we DO still show reject VBs as they are part of the VCF header)
                else if (flag.header_only && v->comp_i == VCF_COMP_MAIN) DROP;

                break;

            case DT_FASTA:
                if (!fasta_piz_is_vb_needed (vb_i)) DROP;
                break;
            
            default: {} // no special handling of other data types
        }

        #undef DROP
        
        // conditions in which VB should be written. if false, but needs_recon, VB is still read, but not written (eg reading aux files)
        v->needs_write = 
            v->needs_recon
        &&  !flag_loading_auxiliary; // we're ingesting, but not reconstructing, an auxiliary file

        if (v->needs_write) {
            // mutex: locked:    here (at initialization)
            //        waited on: writer thread, wanting the data
            //        unlocked:  by main thread after txtheader data is handed over.
            mutex_initialize (v->wait_for_data);
            mutex_lock (v->wait_for_data);
        }
    }
}

static int64_t writer_get_plan_num_lines (void)
{
    int64_t num_lines = 0;
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);
    
    for (uint64_t i=0; i < plan_len; i++) 
        num_lines += plan[i].num_lines;

    return num_lines;
}

static inline void writer_drop_lines (VbInfo *v, uint32_t start_line, uint32_t num_lines)
{
    if (!v->is_dropped)
        v->is_dropped = writer_get_is_dropped (v->vblock_i);

    #define DROPL(l) bits_set (v->is_dropped, start_line + (l))
    
    // shortcuts for the most common cases
    if      (num_lines == 1) { DROPL(0); }
    else if (num_lines == 2) { DROPL(0); DROPL(1); } 
    else if (num_lines == 3) { DROPL(0); DROPL(1); DROPL(2); } 
    else if (num_lines == 4) { DROPL(0); DROPL(1); DROPL(2); DROPL(3); } 
    else if (num_lines == 5) { DROPL(0); DROPL(1); DROPL(2); DROPL(3); DROPL(4); } 
    else bits_set_region (v->is_dropped, start_line, num_lines);
    
    #undef DROPL
}

// PIZ main thread
static void writer_drop_plan_item (ReconPlanItem *p, uint32_t drop_flavor)
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
            v->needs_recon = false;
            p->flavor = drop_flavor;
            break;
        
        case PLAN_INTERLEAVE:
            v->needs_recon = VBINFO(p->vb2_i)->needs_recon = false;
            p->flavor = drop_flavor;
            break;
    
        default: {}  // TXTHEADER is never removed by head/tail/lines, 
                     // END_OF_VB will be removed by writer_cleanup_recon_plan_after_filtering if needed
    }
}

static void writer_drop_item_lines_from_start (ReconPlanItem *p, uint32_t drop_lines)
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
            writer_drop_lines (v, 0, (drop_lines+1) / 2);           // round up in case of odd number
            writer_drop_lines (VBINFO(p->vb2_i), 0, drop_lines / 2); // round down in case of odd number
            break;

        default: ABORT0 ("invalid flavor");
    }
}

static void writer_drop_item_lines_from_end (ReconPlanItem *p, uint32_t drop_lines)
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
            writer_drop_lines (v2, v2->num_lines + (drop_lines+1)/2, (drop_lines+1)/2); // round up in case of odd number
            break;
        }
        
        default: ABORT ("invalid flavor: %s", recon_plan_flavors[p->flavor]);
    }
}

static void writer_trim_lines_from_plan_start (int64_t lines_to_trim)
{
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);

    // remove plan items before first_line
    int64_t lines_so_far = 0;
    for (uint64_t i=0; i < plan_len && lines_so_far < lines_to_trim; i++) {
        ReconPlanItem *p = &plan[i];
        if (!p->num_lines) continue; 

        lines_so_far += p->num_lines;
        
        // case: plan[i] is fully included in lines to be removed
        if (lines_so_far <= lines_to_trim)  
            writer_drop_plan_item (p, PLAN_REMOVE_ME);
    
        // case: plan[i] is partially included in lines to be removed
        else 
            writer_drop_item_lines_from_start (p, p->num_lines - (lines_so_far - lines_to_trim));
    }
}

static void writer_trim_lines_from_plan_end (int64_t num_plan_lines, int64_t lines_to_trim)
{
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);

    // remove plan items after last_line
    int64_t lines_so_far = 0;
    for (int64_t i=plan_len-1; i >= 0 && lines_so_far <= lines_to_trim; i--) {
        ReconPlanItem *p = &plan[i];
        if (!p->num_lines) continue; 

        lines_so_far += p->num_lines;
        
        // case: plan[i] is fully included in lines to be removed
        if (lines_so_far <= lines_to_trim)  
            writer_drop_plan_item (p, PLAN_REMOVE_ME);
    
        // case: plan[i] is partially included in lines to be removed
        else 
            writer_drop_item_lines_from_end (p, p->num_lines - (lines_so_far - lines_to_trim));
    }
}

static void writer_downsample_plan (void)
{
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);

    uint64_t lines_so_far = 0;

    for (uint64_t i=0; i < plan_len; i++) {
        ReconPlanItem *p = &plan[i];

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

static void writer_cleanup_recon_plan_after_filtering (void)
{
    // VBs that have all their lines dropped - make sure they're marked with !need_recon
    for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {
        VbInfo *v = VBINFO(vb_i);
        if (v->needs_recon && v->is_dropped && bits_is_fully_set (v->is_dropped))
            v->needs_recon = false;
    }

    // remove PLAN_REMOVE_ME items or !need_recon from the recon_plan
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);
    uint64_t new_i = 0;

    for (uint64_t old_i=0; old_i < plan_len; old_i++) {
        ReconPlanItem *old_p = &plan[old_i];
        if (!old_p->vb_i ||
            old_p->flavor == PLAN_DOWNSAMPLE || // keep it in plan - just for downsample counting purposes
            (VBINFO(old_p->vb_i)->needs_recon && old_p->flavor != PLAN_REMOVE_ME)) {
                
            if (old_i != new_i) // no need to copy if already in place
                plan[new_i] = *old_p;
            new_i++;
        }
    }

    z_file->recon_plan.len = new_i;
}

static void writer_update_section_list (void)
{
    // alloate in case we need to replace the section_list
    #define new_list (&evb->z_section_headers)
    ASSERTNOTINUSE (*new_list);
    buf_alloc (evb, new_list, z_file->section_list_buf.len, 0, SectionEnt, 1, "z_section_headers"); // note: the data will be moved to z_file->section_list_buf upon commit
    
    // add all TXT_HEADERs and VBs according to the order of their first appearance in the recon_plan
    ARRAY (ReconPlanItem, plan, z_file->recon_plan);
    for (uint64_t i=0; i < plan_len; i++) 
        if (plan[i].flavor == PLAN_TXTHEADER)
            sections_new_list_add_txt_header (new_list, plan[i].comp_i);

        else {
            VbInfo *v = VBINFO(plan[i].vb_i); 
            if (!v->encountered) { // first encounter with this VB
                sections_new_list_add_vb (new_list, v->vblock_i);
                v->encountered = true;
            }

            if (plan[i].flavor == PLAN_INTERLEAVE) {
                VbInfo *v2 = VBINFO(plan[i].vb2_i); 
                if (!v2->encountered) { // first encounter with this VB
                    sections_new_list_add_vb (new_list, v2->vblock_i);
                    v2->encountered = true;
                }
            }
        }

    // case: SAM - add all PRIM VBs that need to be loaded
    // case: FASTQ with --R2 - add R1 VBs needed for pair lookup 
    if ((Z_DT (DT_SAM) && sections_get_num_vbs (SAM_COMP_PRIM)) ||
        (Z_DT(DT_FASTQ) && flag.one_component == 2)) {
        Section vb_header = NULL;
        while (sections_get_next_vb_of_comp_sec ((Z_DT (DT_SAM) ? SAM_COMP_PRIM : FQ_COMP_R1), &vb_header)) 
            if (!VBINFO(vb_header->vblock_i)->encountered)  // VB not already add bc it is in recon_plan
                sections_new_list_add_vb (new_list, vb_header->vblock_i);
    }

    // add SEC_BGZF of all components that have it
    sections_new_list_add_bgzf (new_list);

    // copy all global sections from current section list to new one
    sections_new_list_add_global_sections (new_list);

    // replace section list with new list
    sections_commit_new_list (new_list); 
    #undef new_list
}

// PIZ main thread: add txtheader entry for the component
static void writer_add_txtheader_plan (CompIType comp_i)
{
    if (!B(VbInfo, txt_header_info, comp_i)->needs_write) return; 

    buf_alloc (evb, &z_file->recon_plan, 1, 1000, ReconPlanItem, 1.5, "recon_plan");

    BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
        .flavor    = PLAN_TXTHEADER,
        .comp_i = comp_i
        // num_lines remains 0, and txt_header data is not counted for --head, --tail, --downsample etc
    };
}

// PIZ main thread: add "full vb" entry for each VB of the component
static void writer_add_trivial_plan (CompIType comp_i)
{
    VBIType num_vbs = sections_get_num_vbs (comp_i);
    if (!num_vbs) return;

    buf_alloc (evb, &z_file->recon_plan, num_vbs, 1000, ReconPlanItem, 1.5, "recon_plan");
    buf_alloc_exact_zero (evb, z_file->recon_plan_index, z_file->num_vbs+1, BufWord, "z_file->recon_plan_index");

    if (z_file->recon_plan.len) // we have a txt header item
        *B1ST(BufWord, z_file->recon_plan_index) = (BufWord){ .len = 1 };

    Section vb_header = NULL;
    while (sections_get_next_vb_of_comp_sec (comp_i, &vb_header)) 
        if (VBINFO(vb_header->vblock_i)->needs_recon) {
            BNXT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
                .flavor    = PLAN_FULL_VB, // all lines
                .vb_i      = vb_header->vblock_i,
                .num_lines = B(VbInfo, vb_info, vb_header->vblock_i)->num_lines
            };

            VBINFO(vb_header->vblock_i)->full_vb = true;

            *B(BufWord, z_file->recon_plan_index, vb_header->vblock_i) = 
                (BufWord){ .index = z_file->recon_plan.len - 1, .len = 1};
        }

    // sort by VB as VBs in z_file might be out of order (in DVCF we don't allow out-of-order so no need to sort)
    if (!z_is_dvcf)
        recon_plan_sort_by_vb (z_file);
}

// PIZ main thread: add interleave entry for each VB of the component
// also interleaves VBs of the two components and eliminates the second TXT_HEADER entry
static void writer_add_interleaved_plan (void)
{
    // expecting both componenets to have the same number of VBs
    VBIType R1_num_vbs = sections_get_num_vbs (FQ_COMP_R1);
    VBIType R2_num_vbs = sections_get_num_vbs (FQ_COMP_R2);
    
    ASSERT (R1_num_vbs == R2_num_vbs, "R1 has %u VBs, but R2 has %u VBs - expecting them to be equal when --interleaved",
             R1_num_vbs, R2_num_vbs);

    ASSERT0 (R1_num_vbs, "Component has no VBs");

    buf_alloc (evb, &z_file->recon_plan, R1_num_vbs, 0, ReconPlanItem, 0, "recon_plan");

    for (VBIType i=0; i < R1_num_vbs; i++) {

        VBIType vb1_i = i + sections_get_first_vb_i (FQ_COMP_R1);
        VBIType vb2_2 = i + sections_get_first_vb_i (FQ_COMP_R2);
        
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

// PIZ main thread: for each VB in the component pointed by sec, sort VBs according to first appearance in recon plan
static void writer_add_plan_from_recon_section (CompIType comp_i, bool is_luft,
                                                VBIType *conc_writing_vbs, uint32_t *vblock_mb) // out
{
    Section recon_plan_sec = sections_get_comp_recon_plan_sec (comp_i, is_luft);
    if (!recon_plan_sec) {
        writer_add_trivial_plan (comp_i);
        return;
    }

    recon_plan_uncompress (recon_plan_sec, conc_writing_vbs, vblock_mb); // into evb->scratch 

    if (evb->scratch.len) { // this component has a non-empty RECON_PLAN 
        ARRAY (ReconPlanItem, plan, evb->scratch);

        buf_alloc (evb, &z_file->recon_plan, plan_len, 0, ReconPlanItem, 0, "recon_plan");      

        // convert v12/13
        if (!VER(14)) {
            for (uint64_t i=0; i < plan_len; i++) 
                if (plan[i].flavor == PLAN_END_OF_VB)
                    plan[i].num_lines = 0; // in v12/13, this was all 1s 
        }

        for (uint64_t i=0; i < plan_len; i++) {
            VbInfo *v = VBINFO(plan[i].vb_i);
            if (v->needs_write) {
                if (plan[i].flavor == PLAN_FULL_VB) {
                    plan[i].num_lines = v->num_lines; // note: in the file format num_lines=0 - improves compression of RECON_PLAN section
                    v->full_vb = true;
                }

                BNXT (ReconPlanItem, z_file->recon_plan) = plan[i];
            }
        }
    }
    
    buf_free (evb->scratch);
}

// PIZ main thread: 
// (1) create reconstruction plan for this z_file, taking account RECON_PLAN sections and interleaving.
// (2) re-writes z_file->section_list to include only TXT_HEADERs/VBs which need reading, and order them in the order
// they will be needed by writer.
void writer_create_plan (void)
{  
    writer_init_vb_info();

    uint32_t vblock_mb=0;
    z_file->max_conc_writing_vbs = 0;

    ASSINP (!flag.interleaved || txt_header_info.len == 2, "--interleave cannot be used because %s was not compressed with --pair", z_name);

    // case: SAM with PRIM/DEPN
    if (z_sam_gencomp) {        
        writer_add_txtheader_plan (SAM_COMP_MAIN);

        writer_add_plan_from_recon_section (SAM_COMP_MAIN, false, &z_file->max_conc_writing_vbs, &vblock_mb);

        z_file->num_components  = 1; // one txt_file
    }

    // case: DVCF
    else if (z_is_dvcf) {
        // rejected variants that are reconstructed to be the first part of the VCF header
        CompIType reject_comp_i = flag.luft ? VCF_COMP_PRIM_ONLY : VCF_COMP_LUFT_ONLY;
        writer_add_txtheader_plan (reject_comp_i);
        writer_add_trivial_plan (reject_comp_i);

        // main TXT_HEADER is reconstructed after the rejects
        writer_add_txtheader_plan (VCF_COMP_MAIN);

        // main component (lines with rejected variants are dropped by reconstructor)
        writer_add_plan_from_recon_section (VCF_COMP_MAIN, flag.luft, &z_file->max_conc_writing_vbs, &vblock_mb);
    }

    // case: VCF with --sort, other than DVCF
    else if (Z_DT(DT_VCF) && flag.sort) 
        writer_add_plan_from_recon_section (VCF_COMP_MAIN, false, &z_file->max_conc_writing_vbs, &vblock_mb);

    // case: paired FASTQ to be written interleaved 
    else if (Z_DT(DT_FASTQ) && flag.interleaved) 
        writer_add_interleaved_plan(); // single txt_file, recon_plan with PLAN_INTERLEAVE to interleave pairs of VBs

    // case: paired FASTQ with --R1 or --R2
    else if (Z_DT(DT_FASTQ) && flag.one_component) 
        writer_add_trivial_plan (flag.one_component-1); // one_component is 1 or 2

    // cases: paired FASTQ with --unbind (genounzip only - no filtering/modifications)
    else if (Z_DT(DT_FASTQ) && flag.unbind) {
        writer_add_txtheader_plan (FQ_COMP_R1); // reading the txt_header is required when unbinding (it contains the filename, digest...)
        writer_add_trivial_plan (FQ_COMP_R1);
        z_file->recon_plan.count = z_file->recon_plan.len; // count - length of plan of R1

        writer_add_txtheader_plan (FQ_COMP_R2);
        writer_add_trivial_plan (FQ_COMP_R2);
    }

    // normal file, not one of the above
    else {
        writer_add_txtheader_plan (COMP_MAIN);
        writer_add_trivial_plan (COMP_MAIN);
    }

    bool has_regions_filter = random_access_has_filter();

    // filtering (these are genocat flags, not available in genounzip): 
    if (flag.maybe_lines_dropped_by_writer || has_regions_filter) {

        int64_t num_lines = (flag.lines_first != NO_LINE || flag.tail) ? writer_get_plan_num_lines() : 0; // calculate if needed
        int64_t lines_to_trim;

        // first - filters based on line ordinal position (--head, --lines, --tail)
        if ((flag.lines_first != NO_LINE || flag.tail) 
            && (lines_to_trim = flag.tail ? MIN_(num_lines - flag.tail, num_lines) : flag.lines_first))
            writer_trim_lines_from_plan_start (lines_to_trim);

        if (flag.lines_last != NO_LINE 
            && (lines_to_trim = num_lines - MIN_(flag.lines_last+1, num_lines)))
            writer_trim_lines_from_plan_end (num_lines, lines_to_trim);

        // second - filters based on contents of line - remove VBs by random_access here, more in reconstructor
        if (has_regions_filter) 
            writer_filter_regions (); // note: only filters out whole VBs, does not mark lines for dropping

        // if --downsample is the only line dropping filter, we can remove some plan items now 
        // note: this is an optimization, lines which are not removed now, will be discarded by the writer thread,
        // and lines that are moved, will still be counted for downsampling purposes by the writer thread.
        if (flag.downsample && !flag.tail && flag.lines_first == NO_LINE && !flag.maybe_lines_dropped_by_reconstructor)
            writer_downsample_plan(); 

        // mark VBs that are fully dropped as !needs_recon + remove REMOVE_ME items from recon_plan
        writer_cleanup_recon_plan_after_filtering();
    }

    // build and commit new section list from recon_plan 
    writer_update_section_list(); 

    if (flag.show_gheader == 2) { // --show-gheader=2 : show modified section list 
        sections_show_section_list (z_file->data_type);
        if (is_genocat) exit_ok();
    }

    // actual number of buffers - the maximum of any recon_plan in this z_file, but not less than 3 or more than num_vbs
    z_file->max_conc_writing_vbs = MIN_(vb_info.len, MAX_(3, z_file->max_conc_writing_vbs));

    if (flag.show_recon_plan) {
        recon_plan_show (z_file, flag.luft, z_file->max_conc_writing_vbs, vblock_mb);    
        if (is_genocat) exit_ok();
    }

#if defined _WIN32 || defined APPLE
    ASSERTW ((uint64_t)z_file->max_conc_writing_vbs * (((uint64_t)vblock_mb) << 20) < MEMORY_WARNING_THREASHOLD,
             "\nWARNING: This file's output will be re-ordered in-memory, which will consume %u MB.\n"
             "Alternatively, use the --unsorted option to avoid in-memory sorting", vblock_mb * z_file->max_conc_writing_vbs);
#endif
}

// -------------------
// Writer thread stuff
// -------------------

#define BGZF_FLUSH_THRESHOLD  (32*1024*1024) 
#define PLAIN_FLUSH_THRESHOLD (4 *1024*1024) 
#define FLUSH_THRESHOLD ((txt_file->codec == CODEC_BGZF) ? BGZF_FLUSH_THRESHOLD : PLAIN_FLUSH_THRESHOLD)

static void writer_write (BufferP buf, uint64_t txt_data_len)
{
    START_TIMER;

    if (!buf->len) return;
    
    file_write (txt_file, STRb(*buf));
    
    txt_file->txt_data_so_far_single += txt_data_len;
    txt_file->disk_so_far            += buf->len;

    buf_free (*buf);

    COPY_TIMER_FULL (wvb, write);
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
        writer_write (&bgzf_vb->z_data, bgzf_vb->txt_data.len);
        dispatcher_recycle_vbs (dispatcher, true);  // also release VB
    }

    return !!bgzf_vb;
}

static void writer_flush_vb (Dispatcher dispatcher, VBlockP vb, bool is_txt_header, bool is_last)
{
    if (!vb->txt_data.len32 && !is_last) return; // no data to flush

    if (!flag.no_writer) { // note: we might have a writer thread despite no writing - eg for calculated the digest if we have SAM gencomp
        // case: BGZF compression - offload compression to compute threads (a max of BGZF_FLUSH_THRESHOLD per thread)
        if (dispatcher) {
            // in case of wvb, write the entire buffer, if its "native" vb, chop it up
            uint32_t chunk_size = (vb==wvb) ? vb->txt_data.len32 : BGZF_FLUSH_THRESHOLD;

            for (uint32_t i=0; i < vb->txt_data.len32 || is_last; i += chunk_size) { // if txt_data is empty, we still do one iteration in case of is_last

                if (vb_pool_is_full (POOL_BGZF))
                    writer_output_one_processed_bgzf (dispatcher, true); 

                bgzf_dispatch_compress (dispatcher, Btxt(i), MIN_(chunk_size, vb->txt_data.len32 - i), 
                                        is_last && (i+chunk_size >= vb->txt_data.len32));
            
                is_last = false; // if is is_last, we enter the loop once if txt_data.len=0
            }

            buf_free (vb->txt_data);
        }

        // case: no BGZF, so we can write to disk right away
        else if (vb->txt_data.len) 
            writer_write (&vb->txt_data, vb->txt_data.len);
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
    ASSERT (!v->vb->scratch.len, "expecting vb=%s/%u data to be BGZF-compressed by writer at flush, but it is already compressed by reconstructor: txt_file->codec=%s compressed.len=%"PRIu64" txt_data.len=%"PRIu64,
            comp_name (v->vb->comp_i), v->vb->vblock_i, codec_name (txt_file->codec), v->vb->scratch.len, v->vb->txt_data.len);

    if (!v->vb->txt_data.len) return; // no data in the VB 

    ARRAY (uint32_t, lines, v->vb->lines); 

    ASSERT (start_line + num_lines <= lines_len, "vb_i=%u invalid range: least start_line=%u + num_lines=%u == %u but lines_len=%"PRIu64, 
            BNUM (vb_info, v), start_line, num_lines, start_line + num_lines, lines_len);

    for (uint32_t line_i=start_line; line_i < start_line + num_lines; line_i++) {
        
        rom start = Bc (v->vb->txt_data, lines[line_i]);
        rom after = Bc (v->vb->txt_data, lines[line_i+1]);  // note: lines has one extra entry so this is always correct
        uint32_t line_len = (uint32_t)(after - start);
        
        bool is_dropped = bits_get (v->is_dropped, line_i);

        ASSERT (after >= start, "vb_i=%u Writing line %i (start_line=%u num_lines=%u): expecting start=%p <= after=%p", 
                BNUM (vb_info, v), line_i, start_line, num_lines, start, after);
 
        if (!is_dropped) { // don't output lines dropped in container_  reconstruct_do due to vb->drop_curr_line
            
            if (writer_line_survived_downsampling(v))
                buf_add_more (wvb, &wvb->txt_data, start, line_len, "txt_data");
            
            txt_file->lines_written_so_far++; // increment even if downsampled-out, but not if filtered out during reconstruction (for downsampling accounting)
            v->vb->num_nondrop_lines--;       // update lines remaining to be written (initialized during reconstruction in container_reconstruct_do)
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
        
        bool is_dropped1 = bits_get (v1->is_dropped, line_i);
        bool is_dropped2 = bits_get (v2->is_dropped, line_i);

        // skip lines dropped in container_reconstruct_do due to vb->drop_curr_line for either leaf
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
static void writer_load_vb (Dispatcher dispatcher, VbInfo *v)
{
    bool is_comp = (v < B1ST (VbInfo, vb_info) || v > BLST (VbInfo, vb_info));

    if (flag.show_threads && !is_comp) 
        iprintf ("writer: vb=%s/%u WAITING FOR VB\n", comp_name(v->comp_i), BNUM (vb_info, v));

    if (flag.show_vblocks && !is_comp) 
        iprintf ("WRITER_WAITING_FOR vb=%s/%u\n", comp_name(v->comp_i), BNUM (vb_info, v));

    else if (flag.show_threads && is_comp) 
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

    if (flag.show_vblocks && !is_comp) 
        iprintf ("WRITER_LOADED vb=%s/%u\n", comp_name(v->comp_i), BNUM (vb_info, v));

    v->is_loaded = true;

    ASSERT (v->vblock_i == v->vb->vblock_i, "Wrong VB: v->vblock_i=%u but v->vb->vblock_i=%u", v->vblock_i, v->vb->vblock_i);
}

// Thread entry point for writer thread - at this point, the reconstruction plan is ready and unmutable
static void writer_main_loop (VBlockP wvb) // same as wvb global variable
{
    ASSERTNOTNULL (wvb);

    threads_set_writer_thread();

    // if we need to BGZF-compress, we will dispatch the compression workload to compute threads
    Dispatcher dispatcher = (!flag.no_writer && txt_file->codec == CODEC_BGZF) ? 
        dispatcher_init ("bgzf", POOL_BGZF, writer_get_max_bgzf_threads(), 0, false, false, NULL, 0, NULL) : NULL;
        
    bool do_digest = piz_need_digest && z_has_gencomp;

    // execute reconstruction plan
    for (uint64_t i=0; i < txt_file->recon_plan.len; i++) { // note: recon_plan.len maybe 0 if everything is filtered out

        ReconPlanItem *p = B(ReconPlanItem, txt_file->recon_plan, i);

        ASSERT (p->vb_i >= 0 && p->vb_i <= vb_info.len, 
                "plan[%"PRIu64"].vb_i=%u expected to be in range [1,%u] ", i, p->vb_i, vb_info.len32);

        VbInfo *v  = p->vb_i ? VBINFO(p->vb_i) 
                   : p->flavor == PLAN_TXTHEADER ? B(VbInfo, txt_header_info, p->comp_i)
                   : NULL;

        VbInfo *v2 = (p->vb_i && p->flavor == PLAN_INTERLEAVE) ? VBINFO(p->vb2_i) : NULL;

        // if data for this VB is not ready yet, wait for it
        if (v && !v->is_loaded && p->flavor != PLAN_DOWNSAMPLE)
            writer_load_vb (dispatcher, v);

        if (v2 && !v2->is_loaded) 
            writer_load_vb (dispatcher, v2);

        switch (p->flavor) {

            case PLAN_TXTHEADER:   
                writer_flush_vb (dispatcher, v->vb, true, false); // write the txt header in its entirety
                vb_release_vb (&v->vb, PIZ_TASK_NAME);
                break;

            case PLAN_FULL_VB:   
                // note: normally digest calculation is done in the compute thread in piz_reconstruct_one_vb, but in
                // case of an unmodified VB that inserts lines from gencomp VBs, we do it here, as we re-assemble the original VB
                if (do_digest) digest_one_vb (v->vb, false, NULL); 

                if (!flag.downsample) {

                    if (flag.show_vblocks) // only displayed for entire VBs, not line ranges etc 
                        iprintf ("VB_FLUSH_FULL_VB(id=%d) vb=%s/%d txt_data.len=%u\n", v->vb->id, comp_name (v->vb->comp_i), v->vb->vblock_i, v->vb->txt_data.len32);

                    writer_flush_vb (dispatcher, wvb, false, false); // flush any remaining unflushed wvb lines from previous VBs
                    writer_flush_vb (dispatcher, v->vb, false, false); // write entire VB
                    txt_file->lines_written_so_far += v->vb->num_nondrop_lines;
                }
                else 
                    writer_write_line_range (v, 0, v->num_lines); // this skips downsampled-out lines

                dispatcher_increment_progress ("txt_write", 1); // done writing VB
                vb_release_vb (&v->vb, PIZ_TASK_NAME);
                break;

            case PLAN_DOWNSAMPLE:
                txt_file->lines_written_so_far += p->num_lines;
                break;

            case PLAN_END_OF_VB: { // done with VB - free the memory (happens after a series of "default" line range entries)
                bool needs_flush = do_digest ? digest_one_vb (v->vb, false, &wvb->txt_data) : true;

                // note: if we're digesting gencomp VBs, we dont flush wvb when finishing a non-digestable VB (eg PRIM/DEPN in SAM/BAM)
                if (needs_flush)
                    writer_flush_vb (dispatcher, wvb, false, false); // flush remaining unflushed lines of this VB

                vb_release_vb (&v->vb, PIZ_TASK_NAME);
                dispatcher_increment_progress ("txt_write", 1); // done writing VB
                break;
            }
            case PLAN_INTERLEAVE:
                writer_write_lines_interleaves (dispatcher, v, v2);
                vb_release_vb (&v->vb, PIZ_TASK_NAME);
                vb_release_vb (&v2->vb, PIZ_TASK_NAME);
                dispatcher_increment_progress ("txt_writeX2", 2); // done writing both VBs
                break;

            default:   // PLAN_RANGE
                writer_write_line_range (v, p->start_line, p->num_lines);
                break;
        }

        if (!do_digest && wvb->txt_data.len > FLUSH_THRESHOLD)
            writer_flush_vb (dispatcher, wvb, false, false);
    }

    ASSERT (!do_digest || !wvb->txt_data.len, "expecting wvb->txt_data.len=%u to be 0. Perhaps a dropped END_OF_VB?", wvb->txt_data.len32); // sanity

    // this might have data (eg with flag.downsample or flag.interleave) or not
    writer_flush_vb (dispatcher, wvb, false, true);

    // output for all remain BGZF (if any)
    if (dispatcher) {
        while (!vb_pool_is_empty (POOL_BGZF)) 
            writer_output_one_processed_bgzf (dispatcher, true);

        bgzf_write_finalize(); // write final data to wvb->z_data
        writer_write (&wvb->z_data, 0);

        dispatcher_finish (&dispatcher, NULL, false, false);
    }
  
    // The recon plan is supposed to end with all VBs flushed due to END_OF_VB, FULL_VB etc
    ASSERT (!wvb->txt_data.len, "Expected wvb to be flushed, but wvb->txt_data.len=%"PRIu64, wvb->txt_data.len);

    vb_release_vb (&wvb, WRITER_TASK_NAME); 

    threads_unset_writer_thread();
}

//---------------------------------------------------
// PIZ interaction with writer: start, stop, handover
//---------------------------------------------------

// PIZ main thread: launch writer thread to write to a single txt_file 
static void writer_start_writing (CompIType unbind_comp_i)
{
    ASSERTNOTNULL (txt_file);

    if (writer_thread != THREAD_ID_NONE || flag.no_writer_thread) return;

    // copy the portion of the reconstruction plan (note: recon_plan may exists only in cases of a single txt file)
    if (!flag.unbind)
        buf_copy (evb, &txt_file->recon_plan, &z_file->recon_plan, ReconPlanItem, 0, 0, 0);
    else if (unbind_comp_i == FQ_COMP_R1) 
        buf_copy (evb, &txt_file->recon_plan, &z_file->recon_plan, ReconPlanItem, 0, z_file->recon_plan.count, 0); // param = length of R1
    else // FQ_COMP_R2
        buf_copy (evb, &txt_file->recon_plan, &z_file->recon_plan, ReconPlanItem, z_file->recon_plan.count, 0, 0);

    wvb = vb_get_vb (POOL_MAIN, WRITER_TASK_NAME, 0, COMP_NONE);
    writer_thread = threads_create (writer_main_loop, wvb);
}

// PIZ main thread: wait for writer thread to finish writing a single txt_file
void writer_finish_writing (bool is_last_txt_file)
{
    if (flag.no_writer_thread || !txt_file || writer_thread == THREAD_ID_NONE) return;

    // wait for thread to complete (possibly it completed already)
    threads_join (&writer_thread); // also sets writer_thread=THREAD_ID_NONE
    
    // all mutexes destroyed by main thread, that created them (because we will proceed to freeing z_file in which they live)    
    if (is_last_txt_file) {
        for (CompIType comp_i=0; comp_i < txt_header_info.len; comp_i++) {
            VbInfo *comp = B(VbInfo, txt_header_info, comp_i);
            mutex_destroy (comp->wait_for_data);
        }
        
        for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {
            VbInfo *v = VBINFO(vb_i);
            mutex_destroy (v->wait_for_data);
            buf_destroy (v->is_dropped_buf);
        }
    }                
}

// PIZ main thread
static void writer_handover (VbInfo *v, VBlockP vb)
{
    if (flag.no_writer_thread) return;

    if (!v->needs_write) return; // we don't need this data as we are not going to write any of it

    writer_start_writing (flag.unbind ? v->comp_i : COMP_NONE); // start writer thread if not already started
    
    // writer thread now owns this VB, and is the only thread that will modify it, until finally destroying it
    v->vb = vb;

    threads_log_by_vb (vb, vb->compute_task, vb->vblock_i ? "HANDING OVER VB" : "HANDING OVER TXT_HEADER", 0);

    mutex_unlock (v->wait_for_data); 

    if (flag.show_vblocks) 
        iprintf ("HANDED_OVER(task=%s id=%u) vb_i=%s/%u txt_data.len=%u\n", PIZ_TASK_NAME, vb->id, comp_name(v->comp_i), vb->vblock_i, vb->txt_data.len32);
}

// PIZ main thread: hand over data from a txtheader whose reconstruction main thread has completed, to the writer thread
void writer_handover_txtheader (VBlockP *txt_header_vb_p)
{
    writer_handover (B(VbInfo, txt_header_info, (*txt_header_vb_p)->comp_i), *txt_header_vb_p);
    *txt_header_vb_p = NULL;
}

// PIZ main thread: hand over data from a VB whose reconstruction compute thread has completed, to the writer thread
void writer_handover_data (VBlockP *vb_p)
{
    writer_handover (VBINFO((*vb_p)->vblock_i), *vb_p);
    *vb_p = NULL;
}

