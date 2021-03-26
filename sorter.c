// ------------------------------------------------------------------
//   sorter.c
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "vblock.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "context.h"
#include "codec.h"
#include "zfile.h"
#include "compressor.h"
#include "profiler.h"
#include "strings.h"
#include "file.h"
#include "threads.h"
#include "sorter.h"
#include "sections.h"

// an entry per VB
typedef struct {
    PosType last_line_pos;        // of last line in VB
    WordIndex last_line_chrom_wi; // of last line in VB
    uint32_t start_i, len;        // index and length of this VB's data in txt_file->line_info
    uint32_t num_lines;           // number of lines in this VB
    uint32_t final_index_i;       // final index entry in which this vb appears
    bool is_unsorted;             // copy of vb->is_unsorted
    bool accessed;                // used for accounting for threads
} ZipVbInfo;

typedef struct {
    // sychronization between the main thread and the Writer thread
    uint32_t vblock_i;
    bool is_loaded;               // has data moving to txt_data completed
    Mutex wait_for_data;          // initialized locked, unlocked when data is moved to txt_data and is_loaded is set
    uint32_t pair_vb_i;           // vb_i of pair vb, or 0 if there isn't any. if pair_vb_i > vblock_i then we are pair_1

    // data handed over from the VB
    Buffer txt_data;              // data to be written to disk
    Buffer bgzf_blocks;           // details of BGZF blocks to be written to disk
    Buffer compressed;            // data already BGZF-compressed in compute thread, except perhaps the flanking regions to be compressed in bgzf_write_to_disk
    Buffer line_start;            // array of num_lines x (char *) - pointer to with txt_data - start of each line

    uint32_t num_lines;           // number of textual lines (not data-type lines) in txt_data. 1 if data type is not textual.
    uint32_t vb_data_size;        // data size of VB in default reconstruction
    uint32_t first_line;          // 1-based line number of first line of this VB in the resulting txt file
    bool is_header;               // should this data be dropped if flag.no_header
} PizVbInfo;

typedef struct {
    unsigned comp_i;
    const SectionListEntry *txt_header_sl;
    uint32_t first_vb_i;
    uint32_t num_vbs;
} PizCompInfo;

// an entry per line of Primary and per line of Luft
typedef struct {
    uint32_t vblock_i;
    uint32_t start_line, num_lines;
    WordIndex chrom_wi;
    PosType start_pos, end_pos;
} LineInfo;

static int writer_thread=0;

static void sorter_show_recon_plan (bool is_luft, uint32_t num_txt_data_bufs, uint32_t vblock_mb)
{
    ARRAY (ReconPlanItem, plan, txt_file->recon_plan);

    iprintf ("\nReconstruction plan%s: entries=%u num_txt_data_bufs=%u x %u MB\n", 
             !z_file->z_flags.dual_coords ? "" : is_luft ? " of LUFT view" : " of PRIMARY view", 
             (unsigned)plan_len, num_txt_data_bufs, vblock_mb);

    for (uint32_t i=0; i < plan_len; i++)
        iprintf ("vb=%u\t%s=%u\t%s=%s\n", plan[i].vb_i, 
                 plan[i].plan_type == PLAN_INTERLEAVE ? "vb2" : "start_line", 
                 plan[i].start_line, 
                 plan[i].num_lines < 0xfffffff0       ? "num_lines"  : "plan_type",
                 plan[i].plan_type == PLAN_END_OF_VB  ? "END_OF_VB"  : 
                 plan[i].plan_type == PLAN_FULL_VB    ? "FULL_VB"    :
                 plan[i].plan_type == PLAN_INTERLEAVE ? "INTERLEAVE" :
                                                        str_int_s (plan[i].num_lines).s); 
}

// --------
// ZIP side
// --------

static void sorter_zip_merge_vb_do (VBlock *vb, DidIType chrom_did_i)
{
    ASSERT_DT_FUNC (vb, sizeof_zip_dataline);
    const unsigned dl_size = DT_FUNC (vb, sizeof_zip_dataline)();
    bool is_luft = !!chrom_did_i; // 0=PRIMARY 1=LUFT

    mutex_lock (txt_file->recon_plan_mutex[is_luft]);

    buf_alloc (evb, &txt_file->line_info[is_luft], vb->lines.len, 0, LineInfo, CTX_GROWTH, 0); // added to evb buf_list in file_initialize_txt_file_data

    WordIndex last_chrom_word_index = WORD_INDEX_NONE;
    PosType last_pos = -1;
    uint32_t start_i = txt_file->line_info[is_luft].len;

    for (uint32_t line_i=0 ; line_i < vb->lines.len ; line_i++) {
        ZipDataLine *dl = (ZipDataLine *)(&vb->lines.data[dl_size * line_i]);
        PosType pos = dl->pos[is_luft];
        WordIndex chrom_word_index = node_word_index (vb, chrom_did_i, dl->chrom_index[is_luft]);

        // special case (eg GVCF) - this line is exactly one position after previous line - just add to previous record
        if (is_luft && chrom_word_index == WORD_INDEX_NONE) {
            // In Luft, exclude rejected lines (we exclude here if sorted, and in vcf_piz_TOPLEVEL_cb_drop_line_if_bad_oSTATUS_or_no_header if not sorted)
        }
        else if (last_pos >= 1 && chrom_word_index == last_chrom_word_index && pos == last_pos+1) {
            LineInfo *last = LASTENT (LineInfo, txt_file->line_info[is_luft]);
            last->end_pos  = pos;
            last->num_lines++; 
        }
        else
            NEXTENT (LineInfo, txt_file->line_info[is_luft]) = (LineInfo){
                .vblock_i   = vb->vblock_i,
                .start_line = line_i,
                .num_lines  = 1,
                .chrom_wi   = chrom_word_index,
                .start_pos  = pos,
                .end_pos    = pos
            };
        
        last_pos = pos;
        last_chrom_word_index = chrom_word_index;
    }

    // store information about this VB, that will help us figure out if the entire file is sorted or not
    buf_alloc_zero (evb, &txt_file->vb_info[is_luft], 0, MIN (1000, vb->vblock_i), ZipVbInfo, CTX_GROWTH, 0); // added to evb buf_list in file_initialize_txt_file_data
    
    *ENT (ZipVbInfo, txt_file->vb_info[is_luft], vb->vblock_i-1) = (ZipVbInfo){ // -1 as vblock_i is 1-based
        .is_unsorted        = vb->is_unsorted[is_luft],
        .last_line_chrom_wi = last_chrom_word_index,
        .last_line_pos      = last_pos,
        .start_i            = start_i,
        .len                = txt_file->line_info[is_luft].len - start_i,
        .num_lines          = (uint32_t)vb->lines.len
    };

    txt_file->vb_info[is_luft].len = MAX (txt_file->vb_info[is_luft].len, vb->vblock_i); // note: vblock_i is 1-based

    mutex_unlock (txt_file->recon_plan_mutex[is_luft]);
}

// ZIP compute thread: merge data from one VB (possibly VBs are out of order) into txt_file->vb_info/line_info
void sorter_zip_merge_vb (VBlock *vb)
{
    sorter_zip_merge_vb_do (vb, CHROM);
    if (z_file->z_flags.dual_coords) // Luft coordinates exist
        sorter_zip_merge_vb_do (vb, DTF(ochrom));
}

static bool sorter_zip_is_file_sorted (bool is_luft)
{
    ZipVbInfo *last_v = NULL;

    for (uint32_t vb_i=0; vb_i < txt_file->num_vbs; vb_i++) {
        ZipVbInfo *v = ENT (ZipVbInfo, txt_file->vb_info[is_luft], vb_i);
        
        // cases file line order is NOT sorted 
        if (v->is_unsorted ||  // this vb is not sorted internally
            (vb_i && (v->last_line_chrom_wi < last_v->last_line_chrom_wi || // this vb's first chrom is smaller than last vb's last chrom
                     (v->last_line_chrom_wi == last_v->last_line_chrom_wi && v->last_line_pos < last_v->last_line_pos)))) // same chrom but decreasing pos 
            return false;

        last_v = v;
    }

    return true; // no evidence of unsortedness found, meaning it is sorted
}

static void sorter_zip_create_index (Buffer *index_buf, bool is_luft)
{
    buf_alloc (evb, index_buf, 0, txt_file->line_info[is_luft].len, uint32_t, 1, "compressed");
    
    ZipVbInfo *v = FIRSTENT (ZipVbInfo, txt_file->vb_info[is_luft]);

    for (uint32_t vb_i=0; vb_i < txt_file->num_vbs; vb_i++, v++) 
        for (uint32_t line_i=v->start_i; line_i < v->start_i + v->len; line_i++)
            NEXTENT (uint32_t, *index_buf) = line_i;
}

static bool is_luft_sorter;
static int sorter_line_cmp (const void *a_, const void *b_)
{
    LineInfo *a = ENT (LineInfo, txt_file->line_info[is_luft_sorter], *(uint32_t *)a_);
    LineInfo *b = ENT (LineInfo, txt_file->line_info[is_luft_sorter], *(uint32_t *)b_);

    // case: different chrome: order by chrom
    if (a->chrom_wi != b->chrom_wi) return a->chrom_wi - b->chrom_wi;

    // case: same chrom, different pos - order by pos (note: at this point pos ranges are gapless)
    if (a->end_pos < b->start_pos) return -1; // a is smaller (don't use substraction as POS is 64b and return is 32b)
    
    if (b->end_pos < a->start_pos) return +1; // b is smaller

    // same chrome and pos ranges overlap (not expected in normal VCFs) - sort by vb_i
    return a->vblock_i - b->vblock_i;
}

static void sorter_zip_get_final_index_i (Buffer *line_info, Buffer *vb_info, Buffer *index_buf)
{
    ARRAY (LineInfo, lorder, *line_info);
    ARRAY (ZipVbInfo, v, *vb_info);
    ARRAY (uint32_t, index, *index_buf);

    uint32_t count_vbs=0;
    uint32_t prev_vb_i=(uint32_t)-1;

    // initialize
    for (unsigned vb_i=0; vb_i < txt_file->num_vbs; vb_i++)
        v[vb_i].final_index_i = 0xffffffff; 

    // scan backwards
    for (int32_t index_i=index_len-1; index_i >= 0 ; index_i--) {
        
        uint32_t i = index[index_i];
        uint32_t vb_i = lorder[i].vblock_i - 1; // -1 bc vblock_i is 1-based

        // case: new VB that was not encountered before (scanning backwards) - this is its final occurance
        if (vb_i != prev_vb_i) {

            if (v[vb_i].final_index_i == 0xffffffff) {
                v[vb_i].final_index_i = index_i;
         
                count_vbs++;
                if (count_vbs == txt_file->num_vbs) break; // we're done
            }
            prev_vb_i = vb_i;
        }
    }
}

// reconstruction plan format: txt_file->recon_plan is an array of ReconPlanItem
// returns - max number of concurrent txt_data buffers in memory needed for execution of the plan
static uint32_t sorter_zip_plan_reconstruction (const Buffer *line_info, const Buffer *vb_info, const Buffer *index_buf)
{
    ARRAY (uint32_t, index, *index_buf);

    buf_alloc (evb, &txt_file->recon_plan, 0, 100000, ReconPlanItem, 1, "txt_file->recon_plan");

    LineInfo *prev_li = NULL;
    uint32_t num_txt_data_bufs=0, max_num_txt_data_bufs=0;
    
    for (uint32_t index_i=0; index_i < index_len; index_i++) {
        
        LineInfo *li = ENT (LineInfo, *line_info, index[index_i]);

        // verify sort
        #ifdef DEBUG
        ASSERTE0 (!prev_li || 
                  li->chrom_wi > prev_li->chrom_wi || 
                  (li->chrom_wi == prev_li->chrom_wi && li->start_pos >= prev_li->start_pos),
                  "line_info not sorted correctly");
        #endif
    
        // same vb and consecutive lines - collapse to a single line
        if (prev_li && 
            li->vblock_i == prev_li->vblock_i &&  
            prev_li->start_line + prev_li->num_lines == li->start_line)
            LASTENT (ReconPlanItem, txt_file->recon_plan)->num_lines += li->num_lines;

        // different vb or non-consecutive lines - start new entry
        else {
            buf_alloc (evb, &txt_file->recon_plan, 1, 0, ReconPlanItem, 2, 0);
            NEXTENT (ReconPlanItem, txt_file->recon_plan) = (ReconPlanItem){
                .vb_i       = li->vblock_i,
                .start_line = li->start_line,
                .num_lines  = li->num_lines
            };
        }
        
        // account for number of concurrent threads - a thread lives between the first and last access to its VB
        ZipVbInfo *v = ENT (ZipVbInfo, *vb_info, li->vblock_i - 1); // -1 as vblock_i is 1-based 
        
        if (!v->accessed) { // case: first access to a VB (in PIZ - we would wait on VB data here)
            v->accessed = true;
            num_txt_data_bufs++;

            if (num_txt_data_bufs > max_num_txt_data_bufs) 
                max_num_txt_data_bufs = num_txt_data_bufs;
        }

        if (v->final_index_i == index_i) { // case: last access to a VB (in PIZ - we would free VB data here)
            num_txt_data_bufs--;

            buf_alloc (evb, &txt_file->recon_plan, 1, 0, ReconPlanItem, 2, 0);
            NEXTENT (ReconPlanItem, txt_file->recon_plan) = (ReconPlanItem){ 
                .vb_i      = li->vblock_i,
                .num_lines = PLAN_END_OF_VB,
            }; 
        }

        prev_li = li;
    }

    return max_num_txt_data_bufs;
}

static void sorter_compress_recon_plan_do (bool is_luft)
{
    #define index_buf evb->compressed // an array of uint32_t - indices into txt_file->line_info

    // case: no need for a thread plan, as file is sorted
    if (sorter_zip_is_file_sorted (is_luft)) return; 

    // build index sorted by VB - indices into txt_file->line_info 
    sorter_zip_create_index (&index_buf, is_luft);

    // sort lines
    is_luft_sorter = is_luft; // ugly

    START_TIMER
    qsort (index_buf.data, index_buf.len, sizeof (uint32_t), sorter_line_cmp);
    COPY_TIMER_VB (evb, sorter_compress_qsort);

    // find final entry in index_buf of each VB - needed to calculate how many concurrent threads will be needed for piz
    sorter_zip_get_final_index_i (&txt_file->line_info[is_luft], &txt_file->vb_info[is_luft], &index_buf);

    // create txt_file->reconstruction_plan
    uint32_t num_txt_data_bufs = sorter_zip_plan_reconstruction (&txt_file->line_info[is_luft], &txt_file->vb_info[is_luft], &index_buf);
    buf_free (&index_buf)

    // get best codec for the reconstruction plan data
    Codec codec = codec_assign_best_codec (evb, 0, &txt_file->recon_plan, SEC_RECON_PLAN);
    if (codec == CODEC_UNKNOWN) codec = CODEC_NONE; // too small for compression

    if (flag.show_recon_plan)
        sorter_show_recon_plan (is_luft, num_txt_data_bufs, (uint32_t)(flag.vblock_memory >> 20));

    // prepare section header and compress
    SectionHeaderReconPlan header = (SectionHeaderReconPlan){
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_RECON_PLAN,
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderReconPlan)),
        .h.data_uncompressed_len = BGEN32 (txt_file->recon_plan.len * sizeof (ReconPlanItem)),
        .h.codec                 = codec,
        .h.flags.recon_plan.luft = is_luft,
        .num_txt_data_bufs       = BGEN32 (num_txt_data_bufs),
        .vblock_mb               = BGEN32 ((uint32_t)(flag.vblock_memory >> 20))
    };

    txt_file->recon_plan.len *= 3; // each ReconPlanItem is 3xuint32_t
    BGEN_u32_buf (&txt_file->recon_plan, 0);

    comp_compress (evb, &evb->z_data, (SectionHeader*)&header, txt_file->recon_plan.data, NULL);
    buf_free (&txt_file->recon_plan);
    buf_free (&txt_file->vb_info[is_luft]);
    buf_free (&txt_file->line_info[is_luft]);
}

void sorter_compress_recon_plan (void)
{
    START_TIMER;

    // reconstruction plans for primary and luft (might not exist)
    if (txt_file->vb_info[0].len) sorter_compress_recon_plan_do (false);
    if (txt_file->vb_info[1].len) sorter_compress_recon_plan_do (true );

    COPY_TIMER_VB (evb, sorter_compress_recon_plan);
}

// -----------------------------
// PIZ reconstruction plan stuff
// -----------------------------

// PIZ with --luft: switch the order of the primary and reject txt files in the section list, to convince
// piz_one_file to piz them in reverse order
void sorter_move_liftover_rejects_to_front (void)
{
    ASSERTE0 (flag.luft, "Expecting flag.luft=true");

    SectionListEntry *primary = (SectionListEntry *)sections_get_first_section_of_type (SEC_TXT_HEADER, false);

    SectionListEntry *rejects = primary;
    sections_get_next_section_of_type ((const SectionListEntry **)&rejects, SEC_TXT_HEADER, false, true);

    SectionListEntry *after = rejects;
    sections_get_next_section_of_type ((const SectionListEntry **)&after, SEC_DICT, false, true); // first non-VB section type

    unsigned num_primary_secs = rejects - primary;
    unsigned num_rejects_sec  = after - rejects;

    // now swap them...
    SectionListEntry *temp = MALLOC (num_rejects_sec * sizeof (SectionListEntry));
    memcpy (temp, rejects, num_rejects_sec * sizeof (SectionListEntry));
    memmove (primary + num_rejects_sec, primary, num_primary_secs * sizeof (SectionListEntry));
    memcpy (primary, temp, num_rejects_sec * sizeof (SectionListEntry));
    FREE (temp);
}

// PIZ without --luft: remove rejects
void sorter_remove_liftover_rejects (void)
{
    ASSERTE0 (!flag.luft, "Expecting flag.luft=false");

    SectionListEntry *primary = (SectionListEntry *)sections_get_first_section_of_type (SEC_TXT_HEADER, false);

    SectionListEntry *rejects = primary;
    sections_get_next_section_of_type ((const SectionListEntry **)&rejects, SEC_TXT_HEADER, false, true);

    SectionListEntry *after = rejects;
    sections_get_next_section_of_type ((const SectionListEntry **)&after, SEC_DICT, false, true); // first non-VB section type

    unsigned num_rejects_sec  = after - rejects;
    unsigned num_after_sec    = AFTERENT (SectionListEntry, z_file->section_list_buf) - after;

    // now remove rejects
    memmove (rejects, after, num_after_sec * sizeof (SectionListEntry));
    z_file->section_list_buf.len -= num_rejects_sec;
    z_file->num_components--;
}


// returns my own number in the pair (1 or 2) and pair_vb_i. return 0 if file is not paired.
unsigned sorter_piz_get_pair (uint32_t vb_i, uint32_t *pair_vb_i)
{
    if (!z_file->z_flags.dts_paired) return 0; // not paired

    PizVbInfo *v = ENT (PizVbInfo, z_file->vb_info[0], vb_i - 1);

    *pair_vb_i = v->pair_vb_i;
    return *pair_vb_i > vb_i ? 1 : 2;
}

static uint32_t sorter_piz_initialize_components (void)
{
    buf_alloc (evb, &z_file->comp_info, 0, z_file->num_components, PizCompInfo, 0, "z_file->comp_info");
    z_file->comp_info.len = z_file->num_components;

    uint32_t num_vbs = 0;
    const SectionListEntry *sl = NULL;

    for (unsigned i=0; i < z_file->num_components; i++) {

        PizCompInfo *comp = ENT (PizCompInfo, z_file->comp_info, i);

        ASSERTE (sections_get_next_section_of_type (&sl, SEC_TXT_HEADER, 0, 0),
                 "Expecting %s to have %u components, but found only %u", z_name, z_file->num_components, i);

        *comp = (PizCompInfo){ .comp_i = i, .txt_header_sl = sl, .first_vb_i = 0xffffffff };

        sections_count_component_vbs (sl, &comp->first_vb_i, &comp->num_vbs);
        num_vbs += comp->num_vbs;
    }

    return num_vbs;
}

// initialize component info. this is run once per z_file
static void sorter_piz_init_vb_info (void)
{
    if (z_file->comp_info.len) return; // already initialized

    z_file->vb_info[0].len = sorter_piz_initialize_components();
    buf_alloc_zero (evb, &z_file->vb_info[0], 0, z_file->vb_info[0].len, PizVbInfo, 1, "z_file->vb_info");

    const SectionListEntry *sl = NULL;
    uint32_t v_i=0;
    for (unsigned comp_i=0; comp_i < z_file->comp_info.len; comp_i++) {

        uint32_t comp_num_vbs = ENT (PizCompInfo, z_file->comp_info, comp_i)->num_vbs;
        for (uint32_t v_comp_i=0; v_comp_i < comp_num_vbs; v_comp_i++, v_i++) {

            ASSERTE0 (sections_get_next_section_of_type (&sl, SEC_VB_HEADER, 0, 0), "Unexpected end of section list");
            
            PizVbInfo *v = ENT (PizVbInfo, z_file->vb_info[0], v_i);
            v->vblock_i = sl->vblock_i;

            // set pairs (used even if not interleaving)
            if (z_file->z_flags.dts_paired) 
                v->pair_vb_i = v->vblock_i + comp_num_vbs * (comp_i % 2 ? -1 : 1);  

            // we will unlock this when VB is reconstructed and data handed over for writer thread to consume.
            // writer thread will wait on this mutex
            mutex_initialize (v->wait_for_data);
            mutex_lock (v->wait_for_data); // to be unlocked when VB compute thread is joined
        }
    }
}

// PIZ main thread: add "all vb" entry for each VB of the component
static void sorter_piz_add_trival_plan (const PizCompInfo *comp)
{
    buf_alloc (evb, &txt_file->recon_plan, comp->num_vbs, 1000, ReconPlanItem, 1.5, "recon_plan");

    for (uint32_t i=0; i < comp->num_vbs; i++) {
        NEXTENT (ReconPlanItem, txt_file->recon_plan) = (ReconPlanItem){ 
            .vb_i       = comp->first_vb_i + i,
            .plan_type  = PLAN_FULL_VB // all lines
        };
    }
}

// PIZ main thread: add interleave entry for each VB of the component
// also interleaves VBs of the two components and eliminates the second TXT_HEADER entry
static void sorter_piz_add_interleave_plan (const PizCompInfo *comp)
{
    const SectionListEntry *sl1 = comp->txt_header_sl;
    const SectionListEntry *sl2 = (comp+1)->txt_header_sl;

    ((SectionListEntry *)sl2)->section_type = SEC_NONE; // effectively merge the components by neutralizing the 2nd TXT_HEADER

    unsigned num_vbs_per_component=0;

    while (sections_get_next_section_of_type2 (&sl1, SEC_VB_HEADER, SEC_NONE, false, false) && // 2nd component TXT_HEADER is not SEC_NONE
           sl1->section_type == SEC_VB_HEADER) {

        ASSERTE0 (sections_get_next_section_of_type2 (&sl2, SEC_VB_HEADER, SEC_TXT_HEADER, false, false) &&
                  sl2->section_type == SEC_VB_HEADER, "Failed to find matching VB in second component, when --interleave");

        const SectionListEntry *last_1 = sections_get_last_section_of_type2 (sl1, SEC_B250, SEC_LOCAL);
        const SectionListEntry *last_2 = sections_get_last_section_of_type2 (sl2, SEC_B250, SEC_LOCAL);
        
        buf_alloc (evb, &txt_file->recon_plan, 2, 1000, ReconPlanItem, 1.5, "recon_plan");
        NEXTENT (ReconPlanItem, txt_file->recon_plan) = (ReconPlanItem){ 
            .vb_i      = sl1->vblock_i,
            .vb2_i     = sl2->vblock_i,
            .plan_type = PLAN_INTERLEAVE // all lines
        };

        sl1 = sections_pull_vb_up (sl2->vblock_i, last_1); // move the 2nd component VB to be after the corresponding 1st component VB
        sl2 = last_2; 
        num_vbs_per_component++;
    }

    ASSERTE0 (num_vbs_per_component, "Component has no VBs");
    
    // we've iterated to the end of component1 VBs, make sure component2 doesn't have any additional VB
    ASSERTE (!sections_get_next_section_of_type2 (&sl2, SEC_VB_HEADER, SEC_TXT_HEADER, false, false) ||
             sl2->section_type == SEC_TXT_HEADER, // either no more VB/TXT headers, or the next header is TXT
             "First component has %u num_vbs_per_component VBs, but second component has more, when --interleave", num_vbs_per_component);
}

// returns true if this component has a SEC_RECON_PLAN section
static bool sorter_piz_is_comp_have_recon_plan (const PizCompInfo *comp)
{
    const SectionListEntry *sl = comp->txt_header_sl;

    return sections_get_next_section_of_type2 (&sl, SEC_RECON_PLAN, SEC_TXT_HEADER, 0, 0) &&
           sl->section_type == SEC_RECON_PLAN;
}

// for each VB in the component pointed by sl, sort VBs according to first appearance in recon plan
static void sorter_piz_add_plan_from_recon_section (const PizCompInfo *comp, 
                                                    uint32_t *num_txt_data_bufs, uint32_t *vblock_mb) // out
{
    zfile_get_global_section (SectionHeaderReconPlan, SEC_RECON_PLAN, comp->txt_header_sl, 
                              &evb->compressed, "compressed");

    // assign outs
    *num_txt_data_bufs = MAX (*num_txt_data_bufs, BGEN32 (header.num_txt_data_bufs));
    *vblock_mb = BGEN32 (header.vblock_mb);

    evb->compressed.len /= sizeof (uint32_t); // len to units of uint32_t
    BGEN_u32_buf (&evb->compressed, 0);
    evb->compressed.len /= 3; // len to units of ReconPlanItem

    // concatenate reconstruction plan
    buf_add_buf (evb, &txt_file->recon_plan, &evb->compressed, ReconPlanItem, "recon_plan");
            
    // track if VB has been pulled up already
    uint8_t *vb_is_pulled_up = CALLOC (comp->num_vbs); // dynamic allocation as number of VBs in a component is unbound

    ARRAY (const ReconPlanItem, plan, evb->compressed);
    
    const SectionListEntry *sl = comp->txt_header_sl;
    for (uint64_t i=0; i < plan_len; i++)
        if (!vb_is_pulled_up[plan[i].vb_i - comp->first_vb_i]) {
            // move all sections of vb_i to be immediately after sl ; returns last section of vb_i after move
            sl = sections_pull_vb_up (plan[i].vb_i, sl); 
            vb_is_pulled_up[plan[i].vb_i - comp->first_vb_i] = true;
        }

    FREE (vb_is_pulled_up);
    buf_free (&evb->compressed);
}

// PIZ main thread: if --sort: read one (if unbind) or all (if concat) reconstruction plans (for primary OR luft) from z_file
// if a component has no plan (because not --sort or not in z_file) - create a trival reconstruction plan
// re-sort all VBs within each component in sections according to VB appearance order in reconstruction plan
static void sorter_piz_create_plan (uint32_t component_i /* 0 if not unbinding */)
{
    if (flag.genocat_no_reconstruct_output && !flag.show_recon_plan) return;

    sorter_piz_init_vb_info();

    unsigned first_comp_i = flag.unbind ? component_i : 0;
    unsigned last_comp_i  = flag.unbind ? component_i : z_file->num_components-1;
    unsigned num_comps = last_comp_i - first_comp_i + 1;
    uint32_t num_txt_data_bufs=0, vblock_mb=0;

    ASSERTE (!flag.interleave || !(num_comps % 2), "%s has %u components, but --interleave expects an even number", 
             z_name, num_comps);

    for (unsigned comp_i = first_comp_i; comp_i <= last_comp_i; comp_i += 1 + flag.interleave) {

        PizCompInfo *comp = ENT (PizCompInfo, z_file->comp_info, comp_i);

        // case: interleave the VBs of two components - we ignore RECON_PLAN sections 
        if (flag.interleave) 
            sorter_piz_add_interleave_plan (comp);

        else if (flag.sort && sorter_piz_is_comp_have_recon_plan (comp))
            sorter_piz_add_plan_from_recon_section (comp, &num_txt_data_bufs, &vblock_mb);

        else 
            sorter_piz_add_trival_plan (comp);
    }

    // actual number of buffers - compute threads: num_txt_data_bufs ; writer thread: num_txt_data_bufs (at least 3) ; but no more than num_vbs
    num_txt_data_bufs = MIN (z_file->vb_info[0].len, MAX (3, num_txt_data_bufs) + global_max_threads);

    if (flag.show_recon_plan) {
        sorter_show_recon_plan (flag.luft, num_txt_data_bufs, vblock_mb);    
        if (exe_type == EXE_GENOCAT) exit_ok;
    }

#if defined _WIN32 || defined APPLE
    ASSERTW ((uint64_t)num_txt_data_bufs * (((uint64_t)vblock_mb) << 20) < MEMORY_WARNING_THREASHOLD,
                "\nWARNING: This file will be output sorted. Sorting is done in-memory and will consume %u MB.\n"
                "Alternatively, use the --unsorted option to avoid in-memory sorting", vblock_mb * num_txt_data_bufs);
#endif
}

// ----------------
// PIZ writer stuff
// ----------------

static void destroy_one_vb_info (PizVbInfo *v)
{
    buf_destroy (&v->txt_data);
    buf_destroy (&v->bgzf_blocks);
    buf_destroy (&v->compressed);
    buf_destroy (&v->line_start);
}

// final filter of output lines, based on downsample and shard. returns true if line is to be output
static inline bool sorter_piz_line_survived_downsampling (void)
{
    if (!flag.downsample) return true;

    uint32_t line_i = txt_file->lines_so_far / (DTPT (line_height) << !!flag.interleave); // double if interleave

    // show (or not) the line based on our downsampling rate and shard value
    return (line_i % flag.downsample) == flag.shard;
}

// write lines one at a time, watching for dropped lines
static void sorter_piz_write_line_range (PizVbInfo *v, uint32_t start_line, uint32_t num_lines)
{
    ARRAY (const char *, lines, v->line_start); 

    for (uint32_t line_i=start_line; line_i < start_line + num_lines; line_i++) {
        const char *start = lines[line_i];
        const char *after = lines[line_i+1];  // note: lines has one extra entry so this is always correct
        uint32_t line_len = after - start;
        if (!(line_len==2 && start[0]=='\b')) { // don't output lines dropped in container_reconstruct_do due to vb->dont_show_curr_line
            
            if (sorter_piz_line_survived_downsampling())
                txtfile_write_to_disk (0, start, line_len);
            
            txt_file->lines_so_far++; // increment even if downsampled-out, but not if filtered out during reconstruction
        }
    }
    txtfile_flush_if_stdout();
}

// write 4 lines at time, interleaved from v and v2
static void sorter_piz_write_lines_interleaves (PizVbInfo *v1, PizVbInfo *v2)
{
    ASSERTE (v1->num_lines % 4 == 0, "vb_i=%u has %u lines, which is not a multiple of 4", v1->vblock_i, v1->num_lines);
    ASSERTE (v1->num_lines == v2->num_lines, "when interleaving, expecting number of lines of vb_i=%u (%u) to be the same as vb_i=%u (%u)",
             v1->vblock_i, v1->num_lines, v2->vblock_i, v2->num_lines);

    for (uint32_t line_i=0; line_i < v1->num_lines; line_i += 4) {

        const char **start1 = ENT (const char *, v1->line_start, line_i);
        const char **start2 = ENT (const char *, v2->line_start, line_i);
        
        // skip lines dropped in container_reconstruct_do due to vb->dont_show_curr_line for either leaf
        if (start1[0][0] != '\b' && start2[0][0] != '\b') { 

            if (sorter_piz_line_survived_downsampling()) { 
                txtfile_write_4_lines (&v1->txt_data, start1, 1);
                txtfile_write_4_lines (&v2->txt_data, start2, 2);
            }

            txt_file->lines_so_far += 8; // increment even if downsampled-out, but not if filtered out during reconstruction
        }
    }

    txtfile_flush_if_stdout();
}

// Thread entry point for writer thread - at this point, the reconstruction plan is ready and unmutable
static void *sorter_piz_writer (void *unused)
{
    ARRAY (ReconPlanItem, plan, txt_file->recon_plan);
    ASSERTE0 (plan_len, "txt_file has no reconstruction plan");

    for (uint64_t i=0; i < plan_len; i++) {

        ASSERTE (plan[i].vb_i >= 1 && plan[i].vb_i <= z_file->vb_info[0].len, 
                 "plan[%u].vb_i=%u expected to be in range [1,%u] ", (unsigned)i, plan[i].vb_i, (unsigned)z_file->vb_info[0].len);

        PizVbInfo *v  =                   ENT (PizVbInfo, z_file->vb_info[0], plan[i].vb_i  - 1),
                  *v2 = flag.interleave ? ENT (PizVbInfo, z_file->vb_info[0], plan[i].vb2_i - 1) : NULL;

        // if data for this VB is not ready yet, wait for it
        if (!v->is_loaded) {
            if (flag.show_threads) iprintf ("Writer thread: waiting for data from vb_i=%u\n", plan[i].vb_i);
            mutex_wait (v->wait_for_data);
        }

        if (v2 && !v2->is_loaded) {
            if (flag.show_threads) iprintf ("Writer thread: waiting for data from vb_i=%u\n", plan[i].vb2_i);
            mutex_wait (v2->wait_for_data);
        }

        // case: free data if this is the last line of the VB
        switch (plan[i].num_lines) {

            case PLAN_FULL_VB:   // write entire VB
                if (!(v->is_header && flag.no_header)) {

                    if (!flag.may_drop_lines) {
                        txtfile_write_one_vblock (&v->txt_data, &v->bgzf_blocks, &v->compressed, 
                                                  plan[i].vb_i, v->vb_data_size, v->num_lines, v->first_line);
                        
                        txt_file->lines_so_far += v->num_lines; // textual, not data-type, line
                    }

                    else  // expecting some lines to be dropped - we write lines one at a time, watching for dropped lines
                        sorter_piz_write_line_range (v, 0, v->num_lines);

                    txtfile_flush_if_stdout();
                }
                destroy_one_vb_info (v);
                break;

            case PLAN_END_OF_VB: // done with VB - free the memory (happens after a series of "default" line range entries)
                destroy_one_vb_info (v);
                break;

            case PLAN_INTERLEAVE:
                sorter_piz_write_lines_interleaves (v, v2);
                destroy_one_vb_info (v);
                destroy_one_vb_info (v2);
                break;

            default: 
                sorter_piz_write_line_range (v, plan->start_line, plan->num_lines);
                break;
        }
    }
    return NULL;
}

// PIZ main thread: hand over data from a VB whose reconstruction compute thread has completed, to the writer thread
void sorter_piz_handover_data (VBlock *vb)
{
    if (flag.genocat_no_reconstruct_output) return;

    PizVbInfo *v = ENT (PizVbInfo, z_file->vb_info[0], vb->vblock_i - 1);

    // move data from the vb to vb_info
    buf_grab (&v->txt_data,    "vb_info->txt_data",    &vb->txt_data);
    buf_grab (&v->bgzf_blocks, "vb_info->bgzf_blocks", &vb->bgzf_blocks);
    buf_grab (&v->compressed,  "vb_info->compressed",  &vb->compressed);
    
    v->num_lines    = DTPT (line_height) ? vb->lines.len * DTPT (line_height) : 1; // textual lines, not Seg lines
    v->vb_data_size = vb->vb_data_size;
    v->first_line   = vb->first_line;
    v->is_header    = vb->is_rejects_vb;
    
    buf_alloc (evb, &v->line_start, 0, v->num_lines + 1, char *, 0, "vb_info->line_start");
    v->line_start.len = v->num_lines + 1;
    *LASTENT (const char *, v->line_start) = AFTERENT (char, v->txt_data); // one extra "start" for after the txt data - makes len calculations easier

    if (DTPT (line_height)) 
        str_split (v->txt_data.data, v->txt_data.len, v->num_lines, '\n', FIRSTENT (const char *, v->line_start), NULL);
    else // binary data types - the whole txt_data is a single "line"
        *FIRSTENT (const char *, v->line_start) = v->txt_data.data;

    v->is_loaded = true;
    mutex_unlock (v->wait_for_data);
}

// PIZ main thread: create reconstruction plan for the requested component(s) and launch writer thread
void sorter_piz_start_writing (uint32_t component_i /* 0 if not unbinding */)
{
    sorter_piz_create_plan (component_i);

    writer_thread = threads_create (sorter_piz_writer, NULL, "writer_thread", THREADS_NO_VB);
}

void sorter_piz_finish_writing (void)
{
    if (flag.genocat_no_reconstruct_output || !txt_file) return;

    // wait for thread to complete (possibly it completed already)
    threads_join (writer_thread);

    // sanity check
    for (uint32_t i=0; i < z_file->vb_info[0].len; i++) {
        PizVbInfo *v = ENT (PizVbInfo, z_file->vb_info[0], i);
        ASSERTNOTINUSE (v->txt_data);
        ASSERTNOTINUSE (v->bgzf_blocks);
        ASSERTNOTINUSE (v->compressed);
        ASSERTNOTINUSE (v->line_start);
    }
}