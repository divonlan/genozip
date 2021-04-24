// ------------------------------------------------------------------
//   writer.c
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "context.h"
#include "profiler.h"
#include "strings.h"
#include "file.h"
#include "threads.h"
#include "linesorter.h"
#include "sections.h"
#include "random_access.h"
#include "writer.h"
#include "zfile.h"
#include "bgzf.h"
#include "mutex.h"
#include "codec.h" 

// ---------------
// Data structures
// ---------------

typedef struct {
    uint32_t comp_i;              // component_i of this VB
    bool is_loaded;               // has data moving to txt_data completed
    Mutex wait_for_data;          // initialized locked, unlocked when data is moved to txt_data and is_loaded is set
    uint32_t pair_vb_i;           // vb_i of pair vb, or 0 if there isn't any. if pair_vb_i > vblock_i then we are pair_1
    bool in_plan;                 // this VB is used in the reconstruction plan
    bool no_read;                 // entire VB filtered out and should not be read or reconstructed (despite maybe being in the plan)
                                  // if no_read=false, VB should be reconstructed, but writing it depends on in_plan
    VBlock *vb;                   // data handed over main -> compute -> main -> writer
} VbInfo;

typedef struct {
    VbInfo info;                  // data in here (in_plan, no_read etc) refers to entire component, vb is txt_header
    const SecLiEnt *txt_header_sl;
    uint32_t first_vb_i;
    uint32_t num_vbs;
    bool liftover_rejects;        // this component is the "liftover rejects" component
} CompInfo;

typedef struct {
    uint32_t txt_file_i;
    uint32_t first_comp_i, num_comps;
    uint32_t first_plan_i, plan_len;
} TxtFileInfo;

#define vb_info       z_file->vb_info[0] // we use only [0] on the PIZ side
#define comp_info     z_file->comp_info
#define txt_file_info z_file->txt_file_info

static ThreadId writer_thread = THREAD_ID_NONE;

// --------------------------------
// VBlock and Txt Header properties
// --------------------------------

// returns my own number in the pair (1 or 2) and pair_vb_i. return 0 if file is not paired.
unsigned writer_get_pair (uint32_t vb_i, uint32_t *pair_vb_i)
{
    if (!z_file->z_flags.dts_paired) return 0; // not paired

    VbInfo *v = ENT (VbInfo, vb_info, vb_i);

    *pair_vb_i = v->pair_vb_i;
    return *pair_vb_i > vb_i ? 1 : 2;
}

bool writer_is_liftover_rejects (uint32_t component_i)
{
    bool is_liftover_rejects = comp_info.len && // will fail if loading an auxiliary file
                               ENT (CompInfo, comp_info, component_i)->liftover_rejects;
    return is_liftover_rejects;
}

bool writer_is_txtheader_in_plan (uint32_t component_i)
{
    bool is_txtheader_in_plan = comp_info.len && // will fail if loading an auxiliary file
                                ENT (CompInfo, comp_info, component_i)->info.in_plan;
    return is_txtheader_in_plan;                             
}

bool writer_is_component_no_read (uint32_t component_i)
{
    bool is_component_skipped = comp_info.len && // will be false if loading an auxiliary file
                                ENT (CompInfo, comp_info, component_i)->info.no_read;
    return is_component_skipped;                             
}

bool writer_is_vb_no_read (uint32_t vb_i)
{
    bool no_read = !vb_info.len || ENT (VbInfo, vb_info, vb_i)->no_read;
    return no_read;
}

void writer_get_txt_file_info (uint32_t *first_comp_i, uint32_t *num_comps, ConstSecLiEntP *start_sl) // out
{
    const TxtFileInfo *tf = ENT (TxtFileInfo, txt_file_info, z_file->num_txt_components_so_far);
    const CompInfo *first_comp = ENT (CompInfo, comp_info, tf->first_comp_i);

    *first_comp_i = tf->first_comp_i;
    *num_comps    = tf->num_comps;
    *start_sl     = first_comp->txt_header_sl->offset ? first_comp->txt_header_sl - 1
                                                      : NULL; /* before start of array */
}

// -----------------------------
// PIZ reconstruction plan stuff
// -----------------------------

// PIZ main thread
static void sort_piz_init_txt_file_info (void)
{
    buf_alloc (evb, &txt_file_info, 0, flag.unbind || flag.one_component ? comp_info.len : 1, TxtFileInfo, 1, "txt_file_info"); // maximal size

    if (flag.unbind || flag.one_component) {
        uint32_t txt_file_i = 0;
        for (uint32_t comp_i=0; comp_i < comp_info.len; comp_i++) {
            
            bool two_components = flag.interleave || ENT (CompInfo, comp_info, comp_i)->liftover_rejects;

            NEXTENT (TxtFileInfo, txt_file_info) = (TxtFileInfo){
                .txt_file_i   = txt_file_i++,
                .first_comp_i = comp_i,
                .num_comps    = 1 + two_components
            };

            if (two_components) comp_i++;
        }
    }
    else // concatenating all components to a single txt file
        NEXTENT (TxtFileInfo, txt_file_info) = (TxtFileInfo){ .num_comps = comp_info.len };
}

// PIZ main thread
static uint32_t writer_init_comp_info (void)
{
    buf_alloc_zero (evb, &comp_info, 0, z_file->num_components, CompInfo, 0, "comp_info");
    comp_info.len = z_file->num_components;

    uint32_t num_vbs = 0;
    const SecLiEnt *sl = NULL;

    for (unsigned comp_i=0; comp_i < z_file->num_components; comp_i++) {

        CompInfo *comp = ENT (CompInfo, comp_info, comp_i);

        ASSERT (sections_next_sec (&sl, SEC_TXT_HEADER, 0, 0),
                "Expecting %s to have %u components, but found only %u", z_name, z_file->num_components, comp_i);

        *comp = (CompInfo){ 
            .info.comp_i      = comp_i, 
            .txt_header_sl    = sl, 
            .first_vb_i       = 0xffffffff,
            .liftover_rejects = sl->flags.txt_header.liftover_rejects,
        };

        // conditions entire component (header and all VBs) should be skipped (i.e. not even read and inspected)
        comp->info.no_read =
           (!flag.luft && comp->liftover_rejects && !flag.one_component)
        ||
           (flag.one_component && flag.one_component-1 != comp_i); // --component specifies a single component, and this is not it 

        // conditions we write the txt header
        comp->info.in_plan =
           !flag_loading_auxiliary
        &&
           !comp->info.no_read
        &&
           !flag.no_header
        &&
           (  !flag.interleave                      // when interleaving, show every header of the first component of every two
           || comp_i % 2 == 0) 
        &&                                          // whether to show this section considering concatenating or unbinding
           (  flag.unbind                           // unbinding: we show all components
           || (flag.luft && comp->liftover_rejects) // concatenating: if --luft, show all rejects components
           || flag.one_component-1 == comp_i        // single component: this is it (note: if dual coord, this will pointing a rejects components since we switched their order)
           || (!flag.one_component && comp_i==0)    // concatenating: show first header (0 if no rejects)
           || (!flag.one_component && (comp-1)->liftover_rejects)); // concatenating: show first header (first after rejects)

        if (comp->info.in_plan) {
            // mutex: locked:    here (at initialization)
            //        waited on: writer thread, wanting the data
            //        unlocked:  by main thread after txtheader data is handed over.
            mutex_initialize (comp->info.wait_for_data);
            mutex_lock (comp->info.wait_for_data);
        } 
        else
            z_file->disk_size_minus_skips -= sections_get_section_size (sl); // remove txt header here, VBs will be removed in writer_init_vb_info

        sections_count_component_vbs (sl, &comp->first_vb_i, &comp->num_vbs);
        num_vbs += comp->num_vbs;
    }

    return num_vbs;
}

// PIZ main thread - called from writer_init_vb_info
static void write_set_first_vb_for_tail (void)
{
    const SecLiEnt *sl = NULL;
    int64_t count_lines = flag.tail;

    while (sections_prev_sec (&sl, SEC_VB_HEADER)) {
        zfile_read_section_header (evb, sl->offset, sl->vblock_i, SEC_VB_HEADER); // reads into evb->compressed

        SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)evb->compressed.data;
        uint32_t vb_num_lines = BGEN32 (header->num_lines); 

        // case: this is the first VB (from the file end) we need
        if (vb_num_lines > count_lines) { 
            txt_file->tail_1st_vb = sl->vblock_i;
            txt_file->tail_1st_line_1st_vb = vb_num_lines - (uint32_t)count_lines;
            if (flag.out_dt != DT_BAM) flag.no_header = true;
            return;
        }

        count_lines -= vb_num_lines;
    }

    // all lines are included - we cancel --tail for this txt_file
    flag.tail = 0;
}

// PIZ main thread - initialize vb, component, txt_file info. this is run once per z_file
static void writer_init_vb_info (void)
{
    if (comp_info.len) return; // already initialized

    // if we have --tail - find the first VB (near the end) that is needed
    if (flag.tail)
        write_set_first_vb_for_tail(); // sets txt_file->tail_1st_vb,tail_1st_line_1st_vb which are hereinafter immutable

    vb_info.len = writer_init_comp_info() + 1; // +1 as first vb_i=1 (entry 0 will be unused) so we have num_vb+1 entries

    sort_piz_init_txt_file_info();

    buf_alloc_zero (evb, &vb_info, 0, vb_info.len, VbInfo, 1, "z_file->vb_info");

    const SecLiEnt *sl = NULL;

    for (unsigned comp_i=0; comp_i < comp_info.len; comp_i++) {

        CompInfo *comp = ENT (CompInfo, comp_info, comp_i);

        for (uint32_t v_comp_i=0; v_comp_i < comp->num_vbs; v_comp_i++) {

            ASSERT0 (sections_next_sec (&sl, SEC_VB_HEADER, 0, 0), "Unexpected end of section list");
            
            VbInfo *v = ENT (VbInfo, vb_info, sl->vblock_i); // note: VBs are out of order in --luft bc writer_move_liftover_rejects_to_front()
            v->comp_i = comp_i;

            // set pairs (used even if not interleaving)
            if (z_file->z_flags.dts_paired) 
                v->pair_vb_i = sl->vblock_i + comp->num_vbs * (comp_i % 2 ? -1 : 1);  

            // conditions this VB should not be read or reconstructed 
            v->no_read = 
                comp->info.no_read                            // entire component is skipped
            ||  (flag.tail && sl->vblock_i < txt_file->tail_1st_vb) // --tail: this VB is too early, not needed
            ||  (flag.one_vb && flag.one_vb != sl->vblock_i)  // --one-vb: user only wants to see a single VB, and this is not it
            ||  (flag.no_header && comp->liftover_rejects)    // --no-header: this a rejects VB which is displayed as a header
            ||  (flag.header_only && !comp->liftover_rejects) // --header-only (except rejects VB)
            ||  !random_access_is_vb_included (sl->vblock_i); // --regions: this VB is excluded

            // conditions in which VB should be written 
            v->in_plan = 
                !v->no_read
            &&  !flag_loading_auxiliary;                      // we're ingesting, but not reconstructing, an auxiliary file

            if (v->in_plan) {
                // mutex: locked:    here (at initialization)
                //        waited on: writer thread, wanting the data
                //        unlocked:  by main thread after txtheader data is handed over.
                mutex_initialize (v->wait_for_data);
                mutex_lock (v->wait_for_data);
            }
            
            if (v->no_read)
                z_file->disk_size_minus_skips -= sections_get_vb_size (sl);
            else 
                z_file->disk_size_minus_skips -= sections_get_vb_skipped_sections_size (sl);
        }
    }
}

// PIZ with --luft: 
// concatenating - all rejects components move to the front
// unbinding - each reject component switches place with its primary component
// note: currently we support only a single primary+rejects dual component. bug 333.
void writer_move_liftover_rejects_to_front (void)
{
    ASSERT0 (flag.luft, "Expecting flag.luft=true");

    const SecLiEnt *primary = (SecLiEnt *)sections_get_first_section_of_type (SEC_TXT_HEADER, false);
    const SecLiEnt *rejects = sections_component_last (primary) + 1;

    sections_pull_component_up (NULL, rejects);
}

// PIZ main thread: add txtheader entry for the component
static void writer_add_txtheader_plan (CompInfo *comp)
{
    if (!comp->info.in_plan) return; 

    buf_alloc (evb, &z_file->recon_plan, 1, 1000, ReconPlanItem, 1.5, "recon_plan");

    NEXTENT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
        .vb_i       = 0, // component, not VB
        .rp_comp_i  = comp->info.comp_i,
        .plan_type  = PLAN_TXTHEADER
    };
}

// PIZ main thread: add "all vb" entry for each VB of the component
static void writer_add_trival_plan (const CompInfo *comp)
{
    buf_alloc (evb, &z_file->recon_plan, comp->num_vbs, 1000, ReconPlanItem, 1.5, "recon_plan");

    for (uint32_t i=0; i < comp->num_vbs; i++) {
        uint32_t vb_i = comp->first_vb_i + i;
        VbInfo *v = ENT (VbInfo, vb_info, vb_i);

        if (!v->in_plan) continue;

        NEXTENT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
            .vb_i       = vb_i,
            .plan_type  = PLAN_FULL_VB // all lines
        };
    }
}

// PIZ main thread: add interleave entry for each VB of the component
// also interleaves VBs of the two components and eliminates the second TXT_HEADER entry
static void writer_add_interleave_plan (const CompInfo *comp)
{
    const SecLiEnt *sl1 = comp->txt_header_sl;
    const SecLiEnt *sl2 = (comp+1)->txt_header_sl;

    unsigned num_vbs_per_component=0;

    while (sections_next_sec2 (&sl1, SEC_VB_HEADER, SEC_TXT_HEADER, false, false) && 
           sl1->st == SEC_VB_HEADER) {

        ASSERT0 (sections_next_sec2 (&sl2, SEC_VB_HEADER, SEC_TXT_HEADER, false, false) &&
                  sl2->st == SEC_VB_HEADER, "Failed to find matching VB in second component, when --interleave");

        VbInfo *v1 = ENT (VbInfo, vb_info, sl1->vblock_i);
        VbInfo *v2 = ENT (VbInfo, vb_info, sl2->vblock_i);

        if (v1->in_plan && v2->in_plan) {            
            buf_alloc (evb, &z_file->recon_plan, 2, 1000, ReconPlanItem, 1.5, "recon_plan");
            NEXTENT (ReconPlanItem, z_file->recon_plan) = (ReconPlanItem){ 
                .vb_i      = sl1->vblock_i,
                .vb2_i     = sl2->vblock_i,
                .plan_type = PLAN_INTERLEAVE // all lines
            };
        }
        else // if one of them is filtered out, both are
            v1->in_plan = v2->in_plan = false;

        // get last section index, before pull up
        const SecLiEnt *last_1 = sections_vb_last (sl1);
        const SecLiEnt *last_2 = sections_vb_last (sl2);

        sl1 = sections_pull_vb_up (sl2->vblock_i, last_1); // move the 2nd component VB to be after the corresponding 1st component VB
        sl2 = last_2; 
        num_vbs_per_component++;
    }

    ASSERT0 (num_vbs_per_component, "Component has no VBs");
    
    // we've iterated to the end of component1 VBs, make sure component2 doesn't have any additional VB
    ASSERT (!sections_next_sec2 (&sl2, SEC_VB_HEADER, SEC_TXT_HEADER, false, false) ||
            sl2->st == SEC_TXT_HEADER, // either no more VB/TXT headers, or the next header is TXT
            "First component has %u num_vbs_per_component VBs, but second component has more, when --interleave", num_vbs_per_component);
}

// PIZ main thread: for each VB in the component pointed by sl, sort VBs according to first appearance in recon plan
static void writer_add_plan_from_recon_section (const CompInfo *comp, const SecLiEnt *recon_plan_sl,
                                                    uint32_t *conc_writing_vbs, uint32_t *vblock_mb) // out
{
    zfile_get_global_section (SectionHeaderReconPlan, SEC_RECON_PLAN, recon_plan_sl, 
                              &evb->compressed, "compressed");

    // assign outs
    *conc_writing_vbs = MAX (*conc_writing_vbs, BGEN32 (header.conc_writing_vbs));
    *vblock_mb = BGEN32 (header.vblock_mb);

    evb->compressed.len /= sizeof (uint32_t); // len to units of uint32_t
    BGEN_u32_buf (&evb->compressed, 0);
    evb->compressed.len /= 3; // len to units of ReconPlanItem

    // concatenate reconstruction plan
    buf_add_buf (evb, &z_file->recon_plan, &evb->compressed, ReconPlanItem, "recon_plan");
            
    // track if VB has been pulled up already
    uint8_t *vb_is_pulled_up = CALLOC (comp->num_vbs); // dynamic allocation as number of VBs in a component is unbound

    ARRAY (const ReconPlanItem, plan, evb->compressed);
    
    const SecLiEnt *sl = comp->txt_header_sl;
    for (uint64_t i=0; i < plan_len; i++)
        if (!vb_is_pulled_up[plan[i].vb_i - comp->first_vb_i]) { // first encounter with this VB in the plan
            // move all sections of vb_i to be immediately after sl ; returns last section of vb_i after move
            sl = sections_pull_vb_up (plan[i].vb_i, sl); 
            vb_is_pulled_up[plan[i].vb_i - comp->first_vb_i] = true;
        }

    FREE (vb_is_pulled_up);
    buf_free (&evb->compressed);
}

// returns SL if this component has a SEC_RECON_PLAN section which matchs flag.luft
static const SecLiEnt *writer_get_recon_plan_sl (const CompInfo *comp)
{
    const SecLiEnt *sl = comp->txt_header_sl;

    if (sections_next_sec2 (&sl, SEC_RECON_PLAN, SEC_TXT_HEADER, 0, 0) && sl->st == SEC_RECON_PLAN) {

        if (sl->flags.recon_plan.luft == flag.luft) 
            return sl;  // first recon plan matches flag.luft
                
        if ((sl+1)->st == SEC_RECON_PLAN && (sl+1)->flags.recon_plan.luft == flag.luft) 
            return sl+1;  // 2nd recon plan matches flag.luft
    }
    
    return NULL; // no RECON_PLAN matching flag.luft was found in this component
}

// PIZ main thread: create reconstruction plan for this z_file, taking account RECON_PLAN sections and interleaving.
// re-sort all VBs within each component in sections according to VB appearance order in reconstruction plan
void writer_create_plan (void)
{  
    // when showing in liftover coordinates, bring the rejects component(s) forward before the primary component(s)
    if (z_file->z_flags.dual_coords && flag.luft) 
        writer_move_liftover_rejects_to_front(); // must be before writer_init_vb_info as components order changes

    writer_init_vb_info();

    if (flag.no_writer && !flag.show_recon_plan) return; // we prepare the *_info_* data even if not writing

    uint32_t conc_writing_vbs=0, vblock_mb=0;

    ASSERT (!flag.interleave || !(z_file->num_components % 2), "%s has %u components, but --interleave expects an even number", 
            z_name, z_file->num_components);

    TxtFileInfo *tf_ent = FIRSTENT (TxtFileInfo, txt_file_info);

    for (unsigned comp_i = 0; comp_i < z_file->num_components; comp_i += 1 + flag.interleave) {

        CompInfo *comp = ENT (CompInfo, comp_info, comp_i);
        const SecLiEnt *recon_plan_sl;

        // start plan for this txt_file (unless already started in previous comp)
        if (!tf_ent->plan_len) 
            tf_ent->first_plan_i = (uint32_t)z_file->recon_plan.len;

        // add this txt_header to the plan if needed
        if (comp->info.in_plan)
            writer_add_txtheader_plan (comp);
                    
        // case: interleave the VBs of two components - we ignore RECON_PLAN sections 
        if (flag.interleave) 
            writer_add_interleave_plan (comp);

        // case: we have a SEC_RECON section (occurs in dual-coordinates files)
        else if (flag.sort && (recon_plan_sl = writer_get_recon_plan_sl (comp)))
            writer_add_plan_from_recon_section (comp, recon_plan_sl, &conc_writing_vbs, &vblock_mb);

        // case: just a plain-vanilla VB
        else 
            writer_add_trival_plan (comp);

        tf_ent->plan_len = (uint32_t)z_file->recon_plan.len - tf_ent->first_plan_i;
        
        if (flag.unbind && !comp->liftover_rejects) tf_ent++; // note: if this is a rejects component, we stay on the same txt file one more component 
    }

    // actual number of buffers - compute threads: conc_writing_vbs ; writer thread: conc_writing_vbs (at least 3) ; but no more than num_vbs
    conc_writing_vbs = MIN (vb_info.len, MAX (3, conc_writing_vbs) + global_max_threads);

    if (flag.show_recon_plan) {
        linesorter_show_recon_plan (z_file, flag.luft, conc_writing_vbs, vblock_mb);    
        if (exe_type == EXE_GENOCAT) exit_ok;
    }

#if defined _WIN32 || defined APPLE
    ASSERTW ((uint64_t)conc_writing_vbs * (((uint64_t)vblock_mb) << 20) < MEMORY_WARNING_THREASHOLD,
             "\nWARNING: This file will be output sorted. Sorting is done in-memory and will consume %u MB.\n"
             "Alternatively, use the --unsorted option to avoid in-memory sorting", vblock_mb * conc_writing_vbs);
#endif
}

// -------------------
// Writer thread stuff
// -------------------

static void writer_flush_vb (VBlockP vb)
{
    START_TIMER;

    if (txt_file->codec == CODEC_BGZF) {
        // compress now if not compressed yet (--interleave, --downsample, sorting)
        if (flag.maybe_vb_modified_by_writer) 
            bgzf_compress_vb (vb); // compress data into vb->compressed (using BGZF blocks from source file or new ones)

        bgzf_write_to_disk (vb); 
    }

    else if (vb->txt_data.len) {
        file_write (txt_file, vb->txt_data.data, vb->txt_data.len);

        txt_file->txt_data_so_far_single += vb->txt_data.len;
        txt_file->disk_so_far            += vb->txt_data.len;
    }

    // free, so in case this is wvb we can use it for more lines
    buf_free (&vb->txt_data);
    buf_free (&vb->compressed);
    buf_free (&vb->bgzf_blocks);

    COPY_TIMER (write);
}

// final filter of output lines, based on downsample and shard. returns true if line is to be output
static inline bool writer_line_survived_downsampling (VbInfo *v)
{
    if (!flag.downsample) return true;

    uint64_t line_i = txt_file->lines_written_so_far / (flag.interleave ? 2ULL : 1ULL); // 2 lines at a time if interleave

    // show (or not) the line based on our downsampling rate and shard value
    return (line_i % (uint64_t)flag.downsample) == flag.shard;
}

// write lines one at a time, watching for dropped lines
static void writer_write_line_range (VBlock *wvb, VbInfo *v, uint32_t start_line, uint32_t num_lines)
{
    ASSERT (!v->vb->compressed.len, "expecting vb_i=%u data to be BGZF-compressed by writer at flush, but it is already compressed by reconstructor: txt_file->codec=%s compressed.len=%"PRIu64" txt_data.len=%"PRIu64,
            v->vb->vblock_i, codec_name (txt_file->codec), v->vb->compressed.len, v->vb->txt_data.len);

    if (!v->vb->txt_data.len) return; // no data in the VB 

    ARRAY (const char *, lines, v->vb->lines); 

    ASSERT (start_line + num_lines <= lines_len, "vb_i=%u Expecting lines_len=%u to be at least %u", 
            ENTNUM (vb_info, v), (unsigned)lines_len, start_line + num_lines);

    for (uint32_t line_i=start_line; line_i < start_line + num_lines; line_i++) {
        
        const char *start = lines[line_i];
        const char *after = lines[line_i+1];  // note: lines has one extra entry so this is always correct
        uint32_t line_len = (uint32_t)(after - start);

        ASSERT (after >= start, "vb_i=%u Writing line %i (start_line=%u num_lines=%u): expecting start=%p <= after=%p", 
                ENTNUM (vb_info, v), line_i, start_line, num_lines, start, after);

        if (line_len) { // don't output lines dropped in container_reconstruct_do due to vb->drop_curr_line
            
            if (writer_line_survived_downsampling(v))
                buf_add_more (wvb, &wvb->txt_data, start, line_len, "txt_data");
            
            txt_file->lines_written_so_far++; // increment even if downsampled-out, but not if filtered out during reconstruction (for downsampling accounting)
            v->vb->num_nondrop_lines--;       // update lines remaining to be written (initialized during reconstruction in container_reconstruct_do)
        }
    }
}

static void writer_write_lines_add_pair (VBlock *wvb, VbInfo *v, const char *start, unsigned len, unsigned pair /* 1 or 2 */)
{
    static const char *suffixes[3] = { "", "/1", "/2" }; // suffixes for pair 1 and pair 2 reads

    SAFE_NUL (&start[len]);
    unsigned qname_len = strcspn (start, " \n\t");
    SAFE_RESTORE;

    const char *sep = &start[qname_len];
    ASSERT (qname_len < len, "Expected to find a newline in the 4-line read:\n%.*s\n", len, start);

    // write up to the separator
    buf_add_more (wvb, &wvb->txt_data, start, qname_len, "txt_data");

    // write suffix if requested, and suffix is not already present
    if (qname_len < 3 || sep[-2] != '/' || sep[-1] != '0' + pair)
        buf_add_more (wvb, &wvb->txt_data, suffixes[pair], 2, "txt_data");

    buf_add_more (wvb, &wvb->txt_data, sep, len - qname_len, "txt_data");
}

// write one fastq "line" (actually 4 textual lines) at time, interleaved from v1 and v2
static void writer_write_lines_interleaves (VBlock *wvb, VbInfo *v1, VbInfo *v2)
{
    if (!v1->vb->txt_data.len || !v2->vb->txt_data.len) return; // no data in the either of the VBs 

    ASSERT (v1->vb->lines.len == v2->vb->lines.len, "when interleaving, expecting number of lines of vb_i=%u (%"PRIu64") to be the same as vb_i=%u (%"PRIu64")",
            v1->vb->vblock_i, v1->vb->lines.len, v2->vb->vblock_i, v2->vb->lines.len);

    for (uint32_t line_i=0; line_i < v1->vb->lines.len; line_i++) {

        const char **start1 = ENT (const char *, v1->vb->lines, line_i);
        const char **start2 = ENT (const char *, v2->vb->lines, line_i);
        unsigned len1 = (unsigned)(*(start1+1) - *start1);
        unsigned len2 = (unsigned)(*(start2+1) - *start2);

        // skip lines dropped in container_reconstruct_do due to vb->drop_curr_line for either leaf
        if (len1 && len2) {

            if (writer_line_survived_downsampling(v1)) { 
                writer_write_lines_add_pair (wvb, v1, *start1, len1, 1); // also adds /1 and /2 if paired
                writer_write_lines_add_pair (wvb, v2, *start2, len2, 2);
            }

            txt_file->lines_written_so_far += 2; // increment even if downsampled-out, but not if filtered out during reconstruction (for downsampling accounting)
        }
    }
}

// writer thread: waiting for data from a VB and loading it
static void writer_load_vb (VbInfo *v)
{
    bool is_comp = (v < FIRSTENT (VbInfo, vb_info) || v > LASTENT (VbInfo, vb_info));

    if (flag.show_threads && !is_comp) 
        iprintf ("writer: vb_i=%u WAITING FOR VB\n", ENTNUM (vb_info, v));
    
    else if (flag.show_threads && is_comp) 
        iprintf ("writer: component_i=%u WAITING FOR TXT_HEADER\n", v->comp_i);

    mutex_wait (v->wait_for_data);

    threads_log_by_vb (v->vb, v->vb->compute_task, v->vb->vblock_i ? "WRITER LOADED VB" : "WRITER LOADED TXT_HEADER", 0);

    v->is_loaded = true;
}

// Thread entry point for writer thread - at this point, the reconstruction plan is ready and unmutable
static void writer_main_loop (VBlockP wvb)
{
    ASSERTNOTNULL (wvb);

    // execute reconstruction plan
    for (uint64_t i=0; i < txt_file->recon_plan.len; i++) { // note: recon_plan.len maybe 0 if everything is filtered out

        ReconPlanItem *p = ENT (ReconPlanItem, txt_file->recon_plan, i);

        ASSERT (p->vb_i >= 0 && p->vb_i <= vb_info.len, 
                "plan[%u].vb_i=%u expected to be in range [1,%u] ", (unsigned)i, p->vb_i, (unsigned)vb_info.len);

        VbInfo *v  = p->vb_i                        ? ENT (VbInfo, vb_info, p->vb_i) 
                   : p->plan_type == PLAN_TXTHEADER ? &ENT (CompInfo, comp_info, p->rp_comp_i)->info
                   :                                  NULL;

        VbInfo *v2 = p->vb_i && flag.interleave ? ENT (VbInfo, vb_info, p->vb2_i) : NULL;

        // if data for this VB is not ready yet, wait for it
        if (v && !v->is_loaded) 
            writer_load_vb (v);

        if (v2 && !v2->is_loaded) 
            writer_load_vb (v2);

        // case: free data if this is the last line of the VB
        switch (p->num_lines) {

            case PLAN_TXTHEADER:   
                writer_flush_vb (v->vb); // write the txt header in its entirety
                vb_release_vb (&v->vb);
                break;

            case PLAN_FULL_VB:   
                if (!flag.downsample) {
                    writer_flush_vb (v->vb); // write entire VB
                    txt_file->lines_written_so_far += v->vb->num_nondrop_lines;
                }
                else 
                    writer_write_line_range (wvb, v, 0, v->vb->lines.len);

                vb_release_vb (&v->vb);
                break;

            case PLAN_END_OF_VB: // done with VB - free the memory (happens after a series of "default" line range entries)
                vb_release_vb (&v->vb);
                break;

            case PLAN_INTERLEAVE:
                writer_write_lines_interleaves (wvb, v, v2);
                vb_release_vb (&v->vb);
                vb_release_vb (&v2->vb);
                break;

            default: 
                writer_write_line_range (wvb, v, p->start_line, p->num_lines);
                break;
        }

        if (wvb->txt_data.len > 4*1024*1024) // flush VB every 4MB
            writer_flush_vb (wvb);
    }

    writer_flush_vb (wvb);
    vb_release_vb (&wvb); 
}

// PIZ main thread: launch writer thread to write to a single txt_file (one or more components)
static void writer_start_writing (uint32_t txt_file_i)
{
    ASSERTNOTNULL (txt_file);

    if (writer_thread != THREAD_ID_NONE || flag.no_writer) return;

    ASSERT (txt_file_i < txt_file_info.len, "txt_file_i=%u is out of range txt_file_info.len=%u", 
            txt_file_i, (unsigned)txt_file_info.len);

    TxtFileInfo *tf = ENT (TxtFileInfo, txt_file_info, txt_file_i);

    // copy the portion of the reconstruction plan related to this txt file 
    buf_copy (evb, &txt_file->recon_plan, &z_file->recon_plan, ReconPlanItem, tf->first_plan_i, tf->plan_len, 0);

    VBlockP wvb = vb_get_vb ("writer", 0);
    writer_thread = threads_create (writer_main_loop, wvb);
}

// PIZ main thread: wait for writer thread to finish writing a single txt_file (one or more components)
void writer_finish_writing (bool is_last_txt_file)
{
    if (flag.no_writer || !txt_file || writer_thread == THREAD_ID_NONE) return;

    // wait for thread to complete (possibly it completed already)
    threads_join (&writer_thread, true); // also sets writer_thread=THREAD_ID_NONE
    
    // all mutexes destroyed by main thread, that created them (not sure its important)    
    if (is_last_txt_file) {
        for (unsigned comp_i=0; comp_i < comp_info.len; comp_i++) {
            CompInfo *comp = ENT (CompInfo, comp_info, comp_i);
            mutex_destroy (comp->info.wait_for_data);
        }
        
        for (uint32_t vb_i=0; vb_i < vb_info.len; vb_i++) {
            VbInfo *v = ENT (VbInfo, vb_info, vb_i);
            mutex_destroy (v->wait_for_data);
        }
    }                
}

// PIZ main thread
static void writer_handover (VbInfo *v, VBlock *vb)
{
    if (flag.no_writer) return;

    if (!v->in_plan) return; // we don't need this data as we are not going to write any of it

    writer_start_writing (flag.unbind ? v->comp_i : 0); // start writer thread if not already started

    v->vb = vb; // writer thread now owns this VB, and is the only thread that will modify it, until finally destroying it

    threads_log_by_vb (vb, vb->compute_task, vb->vblock_i ? "HANDING OVER VB" : "HANDING OVER TXT_HEADER", 0);

    mutex_unlock (v->wait_for_data); 
}

// PIZ main thread: hand over data from a txtheader whose reconstruction main thread has completed, to the writer thread
void writer_handover_txtheader (VBlockP *comp_vb_p, uint32_t component_i)
{
    writer_handover (&ENT (CompInfo, comp_info, component_i)->info, *comp_vb_p);
    *comp_vb_p = NULL;
}

// PIZ main thread: hand over data from a VB whose reconstruction compute thread has completed, to the writer thread
void writer_handover_data (VBlockP *vb_p)
{
    writer_handover (ENT (VbInfo, vb_info, (*vb_p)->vblock_i), *vb_p);
    *vb_p = NULL;
}

