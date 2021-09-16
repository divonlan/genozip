// ------------------------------------------------------------------
//   linesorter.c
//   Copyright (C) 2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "context.h"
#include "codec.h"
#include "zfile.h"
#include "compressor.h"
#include "strings.h"
#include "file.h"
#include "linesorter.h"
#include "sections.h"
#include "profiler.h"
#include "endianness.h"

// an entry per VB
typedef struct {
    PosType last_line_pos;        // of last line in VB
    WordIndex last_line_chrom_wi; // of last line in VB
    uint32_t start_i, len;        // index and length of this VB's data in txt_file->line_info
    uint32_t num_lines;           // number of lines in this VB
    uint32_t final_index_i;       // final index entry in which this vb appears
    bool is_unsorted;             // copy of vb->is_unsorted
    bool accessed;                // used for accounting for threads
} LsVbInfo;

// a txt_file->line_info entry - one per a group of consecutive lines of Primary or Luft
typedef struct {
    uint32_t vblock_i;
    uint32_t start_line, num_lines;
    WordIndex chrom_wi, prim_chrom_wi;
    PosType start_pos, prim_start_pos, end_pos;
    uint32_t tie_breaker;
} LineInfo;

void linesorter_show_recon_plan (File *file, bool is_luft, uint32_t conc_writing_vbs, uint32_t vblock_mb)
{
    ARRAY (ReconPlanItem, plan, file->recon_plan);

    iprintf ("\nReconstruction plan%s: entries=%u conc_writing_vbs=%u x %u MB\n", 
             !z_dual_coords ? "" : is_luft ? " of LUFT view" : " of PRIMARY view", 
             (unsigned)plan_len, conc_writing_vbs, vblock_mb);

    for (uint32_t i=0; i < plan_len; i++)
        iprintf ("vb=%u\t%s=%u\t%s=%s\n", plan[i].vb_i, 
                 plan[i].plan_type == PLAN_INTERLEAVE ? "vb2"        
               : plan[i].plan_type == PLAN_TXTHEADER  ? "component"
               :                                        "start_line",
                 plan[i].start_line, 
                 plan[i].num_lines < 0xfffffff0        ? "num_lines"   
               :                                         "plan_type",
                 plan[i].plan_type == PLAN_END_OF_VB   ? "END_OF_VB"    
               : plan[i].plan_type == PLAN_FULL_VB     ? "FULL_VB"     
               : plan[i].plan_type == PLAN_INTERLEAVE  ? "INTERLEAVE"  
               : plan[i].plan_type == PLAN_TXTHEADER   ? "TXTHEADER"   
               :                                         str_int_s (plan[i].num_lines).s); 
}

static void linesorter_merge_vb_do (VBlock *vb, bool is_luft)
{
    ASSERT_DT_FUNC (vb, sizeof_zip_dataline);
    const unsigned dl_size = DT_FUNC (vb, sizeof_zip_dataline)();
    
    mutex_lock (txt_file->recon_plan_mutex[is_luft]);

    buf_alloc (evb, &txt_file->line_info[is_luft], vb->lines.len, 0, LineInfo, CTX_GROWTH, 0); // added to evb buf_list in file_initialize_txt_file_data

    WordIndex last_chrom_word_index = WORD_INDEX_NONE;
    PosType last_pos = -1;
    uint32_t start_i = txt_file->line_info[is_luft].len;

    for (uint32_t line_i=0 ; line_i < vb->lines.len ; line_i++) {
        ZipDataLine *dl = (ZipDataLine *)(&vb->lines.data[dl_size * line_i]);
        PosType pos = dl->pos[is_luft];
        WordIndex chrom_word_index = node_word_index (vb, (is_luft ? DTF(ochrom) : CHROM), dl->chrom[is_luft]);

        // exclude rejected lines (we exclude here if sorted, and in vcf_lo_piz_TOPLEVEL_cb_filter_line if not sorted)
        if (chrom_word_index != WORD_INDEX_NONE) 
            NEXTENT (LineInfo, txt_file->line_info[is_luft]) = (LineInfo){
                .vblock_i        = vb->vblock_i,
                .start_line      = line_i,
                .num_lines       = 1,
                .chrom_wi        = chrom_word_index,
                .start_pos       = pos,
                .end_pos         = pos + MAX_(0, dl->end_delta),
                .tie_breaker     = dl->tie_breaker,

                // Primary coord data (for duplicate detection in Luft)
                .prim_chrom_wi  = is_luft ? node_word_index (vb, CHROM, dl->chrom[0]) : WORD_INDEX_NONE,
                .prim_start_pos = is_luft ? dl->pos[0] : 0,
            };
        
        last_pos = pos;
        last_chrom_word_index = chrom_word_index;
    }

    // store information about this VB, that will help us figure out if the entire file is sorted or not
    buf_alloc_zero (evb, &txt_file->vb_info[is_luft], 0, MIN_(1000, vb->vblock_i), LsVbInfo, CTX_GROWTH, 0); // added to evb buf_list in file_initialize_txt_file_data
    
    *ENT (LsVbInfo, txt_file->vb_info[is_luft], vb->vblock_i-1) = (LsVbInfo){ // -1 as vblock_i is 1-based
        .is_unsorted        = vb->is_unsorted[is_luft],
        .last_line_chrom_wi = last_chrom_word_index,
        .last_line_pos      = last_pos,
        .start_i            = start_i,
        .len                = txt_file->line_info[is_luft].len - start_i,
        .num_lines          = (uint32_t)vb->lines.len
    };

    txt_file->vb_info[is_luft].len = MAX_(txt_file->vb_info[is_luft].len, vb->vblock_i); // note: vblock_i is 1-based

    mutex_unlock (txt_file->recon_plan_mutex[is_luft]);
}

// ZIP compute thread: merge data from one VB (possibly VBs are out of order) into txt_file->vb_info/line_info
void linesorter_merge_vb (VBlock *vb)
{
    linesorter_merge_vb_do (vb, false);
    if (z_dual_coords) // Luft coordinates exist
        linesorter_merge_vb_do (vb, true);
}

static bool linesorter_is_file_sorted (bool is_luft)
{
    LsVbInfo *last_v = NULL;

    for (uint32_t vb_i=0; vb_i < txt_file->num_vbs; vb_i++) {
        LsVbInfo *v = ENT (LsVbInfo, txt_file->vb_info[is_luft], vb_i);
        
        // cases file line order is NOT sorted 
        if (v->is_unsorted ||  // this vb is not sorted internally
            (vb_i && (v->last_line_chrom_wi < last_v->last_line_chrom_wi || // this vb's first chrom is smaller than last vb's last chrom
                     (v->last_line_chrom_wi == last_v->last_line_chrom_wi && v->last_line_pos < last_v->last_line_pos)))) // same chrom but decreasing pos 
            return false;

        last_v = v;
    }

    return true; // no evidence of unsortedness found, meaning it is sorted
}

static void linesorter_create_index (Buffer *index_buf, bool is_luft)
{
    buf_alloc (evb, index_buf, 0, txt_file->line_info[is_luft].len, uint32_t, 1, "compressed");
    
    LsVbInfo *v = FIRSTENT (LsVbInfo, txt_file->vb_info[is_luft]);

    for (uint32_t vb_i=0; vb_i < txt_file->num_vbs; vb_i++, v++) 
        for (uint32_t line_i=v->start_i; line_i < v->start_i + v->len; line_i++)
            NEXTENT (uint32_t, *index_buf) = line_i;
}

// ZIP: detect duplicate variants in Luft, that are *not* duplicate in Primary (i.e. the duplication is due to many-to-one mapping in the chain file)
// and output to the a dups file, which can be used with genocat --regions-file or bcftools --regions-file
static void line_sorter_detect_duplicates (const Buffer *index_buf)
{
    static Buffer dups = { .name = "dups" };

    ARRAY (uint32_t, index, *index_buf);
    ARRAY (const LineInfo, li, txt_file->line_info[1]);

    int64_t last_dup=-1;
    for (int64_t index_i=1; index_i < index_len; index_i++) {

        const LineInfo *this_li = &li[index[index_i]];
        const LineInfo *prev_li = &li[index[index_i-1]];

        bool is_dup_luft = (this_li->chrom_wi      == prev_li->chrom_wi)      && (this_li->start_pos      == prev_li->start_pos);
        bool is_dup_prim = (this_li->prim_chrom_wi == prev_li->prim_chrom_wi) && (this_li->prim_start_pos == prev_li->prim_start_pos);
        if (is_dup_luft && !is_dup_prim) {
            if (last_dup != index_i-1) // not duplicate duplicate (i.e. 3 or more variants are the same)
                bufprintf (evb, &dups, "%s\t%"PRIu64"\t%"PRIu64"\n", ctx_get_zf_nodes_snip (ZCTX(DTFZ(ochrom)), this_li->chrom_wi), this_li->start_pos, this_li->start_pos);

            last_dup = index_i;
        }
    }

    if (dups.len) { 
        char fn[strlen(z_name) + 30];
        sprintf (fn, "%s.overlaps", z_name);
        buf_dump_to_file (fn, &dups, 1, false, true, false);
        
        buf_free (&dups);

        WARN ("FYI: Genozip detected cases of two or more variants with different Primary coordinates, mapped to the same Luft coordinates.\n"
              "These duplicates (in Luft coordinates) were output to %s, and may be filtered out with:\ngenocat %s --luft --regions-file ^%s\n",
              fn, z_name, fn);
    }
}

static bool is_luft_sorter;
static int linesorter_line_cmp (const void *a_, const void *b_)
{
    LineInfo *a = ENT (LineInfo, txt_file->line_info[is_luft_sorter], *(uint32_t *)a_);
    LineInfo *b = ENT (LineInfo, txt_file->line_info[is_luft_sorter], *(uint32_t *)b_);

    // case: different chrome: order by chrom
    if (a->chrom_wi != b->chrom_wi) return a->chrom_wi - b->chrom_wi;

    // case: same chrom, different pos - order by start_pos
    if (a->start_pos < b->start_pos) return -1; // a is smaller (don't use substraction as POS is 64b and return is 32b)
    if (b->start_pos < a->start_pos) return +1; // b is smaller

    // case: same chrom, same start_pos - order by end pos
    if (a->end_pos < b->end_pos) return -1; // a is smaller (don't use substraction as POS is 64b and return is 32b)
    if (b->end_pos < a->end_pos) return +1; // b is smaller

    // same chrom and pos ranges - sort by tie_breaker
    if (a->tie_breaker < b->tie_breaker) return -1; // a is smaller (don't use substraction as its unsigned 32b)
    if (b->tie_breaker < a->tie_breaker) return +1; // b is smaller
    
    return 0; // same chrom, pos, tie-breaker -> undeterministic sorting, hopefully this doesn't happen too often
}

static void linesorter_get_final_index_i (Buffer *line_info, Buffer *vb_info, Buffer *index_buf)
{
    ARRAY (LineInfo, lorder, *line_info);
    ARRAY (LsVbInfo, v, *vb_info);
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

// ZIP: reconstruction plan format: txt_file->recon_plan is an array of ReconPlanItem
// returns - max number of concurrent txt_data buffers in memory needed for execution of the plan
static uint32_t linesorter_plan_reconstruction (const Buffer *line_info, const Buffer *vb_info, const Buffer *index_buf)
{
    ARRAY (uint32_t, index, *index_buf);

    buf_alloc (evb, &txt_file->recon_plan, 0, 100000, ReconPlanItem, 1, "txt_file->recon_plan");

    LineInfo *prev_li = NULL;
    uint32_t conc_writing_vbs=0, max_num_txt_data_bufs=0;
    
    for (uint64_t index_i=0; index_i < index_len; index_i++) {
        
        LineInfo *li = ENT (LineInfo, *line_info, index[index_i]);

        // verify sort
        #ifdef DEBUG
        ASSERT (!prev_li || 
                li->chrom_wi > prev_li->chrom_wi || 
                (li->chrom_wi == prev_li->chrom_wi && li->start_pos >= prev_li->start_pos),
                "line_info sorting: [%"PRIu64"]={chrom=%d,pos=%"PRId64"-%"PRId64",nlines=%u} [%"PRIu64"]={chrom=%d,pos=%"PRId64"-%"PRId64",nlines=%u}",
                index_i-1, prev_li->chrom_wi, prev_li->start_pos, prev_li->end_pos, prev_li->num_lines, 
                index_i, li->chrom_wi, li->start_pos, li->end_pos, li->num_lines);
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
        LsVbInfo *v = ENT (LsVbInfo, *vb_info, li->vblock_i - 1); // -1 as vblock_i is 1-based 
        
        if (!v->accessed) { // case: first access to a VB (in PIZ - we would wait on VB data here)
            v->accessed = true;
            conc_writing_vbs++;

            if (conc_writing_vbs > max_num_txt_data_bufs) 
                max_num_txt_data_bufs = conc_writing_vbs;
        }

        if (v->final_index_i == index_i) { // case: last access to a VB (in PIZ - we would free VB data here)
            conc_writing_vbs--;

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

#define index_buf evb->compressed // an array of uint32_t - indices into txt_file->line_info

// ZIP
static void linesorter_compress_recon_plan_do (bool is_luft)
{
    // case: no need for a thread plan, as file is sorted
    if (linesorter_is_file_sorted (is_luft)) return; 

    // build index sorted by VB - indices into txt_file->line_info 
    linesorter_create_index (&index_buf, is_luft);

    // sort lines
    is_luft_sorter = is_luft; // ugly

    START_TIMER
    qsort (index_buf.data, index_buf.len, sizeof (uint32_t), linesorter_line_cmp);
    COPY_TIMER_VB (evb, linesorter_compress_qsort);

    // detect cases where distinct variants in Primary got mapped to the same coordinates in Luft
    if (is_luft) line_sorter_detect_duplicates(&index_buf);

    // find final entry in index_buf of each VB - needed to calculate how many concurrent threads will be needed for piz
    linesorter_get_final_index_i (&txt_file->line_info[is_luft], &txt_file->vb_info[is_luft], &index_buf);

    // create txt_file->recon_plan
    uint32_t conc_writing_vbs = linesorter_plan_reconstruction (&txt_file->line_info[is_luft], &txt_file->vb_info[is_luft], &index_buf);
    buf_free (&index_buf)

    // get best codec for the reconstruction plan data
    Codec codec = codec_assign_best_codec (evb, 0, &txt_file->recon_plan, SEC_RECON_PLAN);
    if (codec == CODEC_UNKNOWN) codec = CODEC_NONE; // too small for compression

    if (flag.show_recon_plan)
        linesorter_show_recon_plan (txt_file, is_luft, conc_writing_vbs, (uint32_t)(flag.vblock_memory >> 20));

    // prepare section header and compress
    SectionHeaderReconPlan header = (SectionHeaderReconPlan){
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_RECON_PLAN,
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderReconPlan)),
        .h.data_uncompressed_len = BGEN32 (txt_file->recon_plan.len * sizeof (ReconPlanItem)),
        .h.codec                 = codec,
        .h.flags.recon_plan.luft = is_luft,
        .conc_writing_vbs       = BGEN32 (conc_writing_vbs),
        .vblock_mb               = BGEN32 ((uint32_t)(flag.vblock_memory >> 20))
    };

    txt_file->recon_plan.len *= 3; // each ReconPlanItem is 3xuint32_t
    BGEN_u32_buf (&txt_file->recon_plan, 0);

    comp_compress (evb, &evb->z_data, (SectionHeader*)&header, txt_file->recon_plan.data, NULL);
    buf_free (&txt_file->recon_plan);
    buf_free (&txt_file->vb_info[is_luft]);
    buf_free (&txt_file->line_info[is_luft]);
}

// ZIP: main thread, called by zip_one_file()
void linesorter_compress_recon_plan (void)
{
    START_TIMER;

    // reconstruction plans for primary and luft (might not exist)
    if (txt_file->vb_info[0].len) linesorter_compress_recon_plan_do (false);
    if (txt_file->vb_info[1].len) linesorter_compress_recon_plan_do (true);

    COPY_TIMER_VB (evb, linesorter_compress_recon_plan);
}

// PIZ: get the maximum value of conc_writing_vbs across all SEC_RECON_PLAN sections in the z_file
uint32_t linesorter_get_max_conc_writing_vbs (void)
{
    Section sl = NULL;
    uint32_t max_conc_writing_vbs=0;

    while (sections_next_sec (&sl, SEC_RECON_PLAN)) {
        uint32_t conc_writing_vbs = zfile_read_section_header (evb, sl->offset, 0, SEC_RECON_PLAN).recon_plan.conc_writing_vbs;
        max_conc_writing_vbs = MAX_(max_conc_writing_vbs, BGEN32 (conc_writing_vbs));
    }

    return max_conc_writing_vbs;
}
