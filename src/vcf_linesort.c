// ------------------------------------------------------------------
//   linesorter.c
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

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
#include "sections.h"
#include "profiler.h"
#include "endianness.h"
#include "segconf.h"
#include "version.h"
#include "recon_plan_io.h"
#include "vcf_private.h"

// an entry per VB
typedef struct {
    LineCmpInfo last_line;        // data from the last line of this VB
    uint32_t start_i, len;        // index and length of this VB's data in txt_file->line_info
    uint32_t num_lines;           // number of lines in this VB
    uint32_t final_index_i;       // final index entry in which this vb appears
    bool is_unsorted;             // copy of vb->is_unsorted
    bool accessed;                // used for accounting for threads
} LsVbInfo;

typedef struct {
    VBIType vblock_i;
    uint32_t start_line;
    LineCmpInfo info;
    WordIndex prim_chrom_wi;
    PosType32 prim_start_pos;
} LineInfo;

bool vcf_is_sorting (CompIType comp_i)
{
    if (chain_is_loaded && (comp_i==VCF_COMP_PRIM_ONLY || comp_i==VCF_COMP_LUFT_ONLY)) // rejects components during DVCF creation
        return false;
    else if (flag.unsorted)
        return false;
    else if (txt_file->coords || z_has_gencomp) // dvcf
        return true;
    else if (flag.sort)
        return true;
    else
        return false;
}

static void vcf_linesort_merge_vb_do (VBlockP vb, bool is_luft)
{
    mutex_lock (txt_file->recon_plan_mutex[is_luft]);

    buf_alloc (evb, &txt_file->line_info[is_luft], vb->lines.len, 0, LineInfo, CTX_GROWTH, 0); // added to evb buf_list in file_initialize_txt_file_data

    uint32_t start_i = txt_file->line_info[is_luft].len32;

    for (uint32_t line_i=0 ; line_i < vb->lines.len32 ; line_i++) {
        ZipDataLineVCF *dl = DATA_LINE(line_i);
        PosType64 pos = dl->pos[is_luft];
        WordIndex chrom_word_index = node_word_index (vb, (is_luft ? DTF(luft_chrom) : DTF(prim_chrom)), dl->chrom[is_luft]);

        // exclude rejected lines (we exclude here if sorted, and in vcf_lo_piz_TOPLEVEL_cb_filter_line if not sorted)
        if (chrom_word_index != WORD_INDEX_NONE) 
            BNXT (LineInfo, txt_file->line_info[is_luft]) = (LineInfo){
                .vblock_i         = vb->vblock_i,
                .start_line       = line_i,
                .info.chrom_wi    = chrom_word_index,
                .info.start_pos   = pos,
                .info.end_pos     = pos + MAX_(0, dl->end_delta),
                .info.tie_breaker = dl->tie_breaker,

                // Primary coord data (for duplicate detection in Luft)
                .prim_chrom_wi    = is_luft ? node_word_index (vb, CHROM, dl->chrom[0]) : WORD_INDEX_NONE,
                .prim_start_pos   = is_luft ? dl->pos[0] : 0,
            };
    }

    // store information about this VB, that will help us figure out if the entire file is sorted or not
    buf_alloc_zero (evb, &txt_file->vb_info[is_luft], 0, MAX_(1000, vb->vblock_i), LsVbInfo, CTX_GROWTH, 0); // added to evb buf_list in file_initialize_txt_file_data

    LineInfo *last_line = BLST (LineInfo, txt_file->line_info[is_luft]);

    *B(LsVbInfo, txt_file->vb_info[is_luft], vb->vblock_i-1) = (LsVbInfo){ // -1 as vblock_i is 1-based
        .is_unsorted = VB_VCF->is_unsorted[is_luft],
        .last_line   = last_line->info,
        .start_i     = start_i,
        .len         = txt_file->line_info[is_luft].len - start_i,
        .num_lines   = vb->lines.len32
    };

    txt_file->vb_info[is_luft].len = MAX_(txt_file->vb_info[is_luft].len, vb->vblock_i); // note: vblock_i is 1-based

    mutex_unlock (txt_file->recon_plan_mutex[is_luft]);
}

// ZIP compute thread (DVCF or sorting): merge data from one VB (possibly VBs are out of order) into txt_file->vb_info/line_info
// (currently used only for VCF - SAM merge happens in main thread in sam_zip_gc_after_compute)
void vcf_linesort_merge_vb (VBlockP vb)
{
    START_TIMER;

    vcf_linesort_merge_vb_do (vb, false);
    if (z_is_dvcf) // Luft coordinates exist
        vcf_linesort_merge_vb_do (vb, true);

    COPY_TIMER(vcf_linesort_merge_vb);
}

static bool vcf_linesort_is_file_sorted (bool is_luft)
{
    LsVbInfo *last_v = NULL;

    for (VBIType vb_i=0; vb_i < txt_file->num_vbs; vb_i++) {
        LsVbInfo *v = B(LsVbInfo, txt_file->vb_info[is_luft], vb_i);

        if (v->is_unsorted ||  // this vb is not sorted internally
            (vb_i && vcf_linesort_cmp (last_v->last_line, v->last_line) > 0)) // the previous and this VBs are unsorted between them
            return false; 

        last_v = v;
    }

    return true; // no evidence of unsortedness found, meaning it is sorted
}

static void vcf_linesort_create_index (BufferP index_buf, bool is_luft)
{
    buf_alloc (evb, index_buf, 0, txt_file->line_info[is_luft].len, uint32_t, 1, "scratch");
    
    LsVbInfo *v = B1ST (LsVbInfo, txt_file->vb_info[is_luft]);

    for (VBIType vb_i=0; vb_i < txt_file->num_vbs; vb_i++, v++) 
        for (uint32_t line_i=v->start_i; line_i < v->start_i + v->len; line_i++)
            BNXT32 (*index_buf) = line_i;
}

// ZIP: detect duplicate variants in Luft, that are *not* duplicate in Primary (i.e. the duplication is due to many-to-one mapping in the chain file)
// and output to the a dups file, which can be used with genocat --regions-file or bcftools --regions-file
static void line_sorter_detect_duplicates (ConstBufferP index_buf)
{
    static Buffer dups = { .name = "dups" };

    ARRAY (uint32_t, index, *index_buf);
    ARRAY (const LineInfo, li, txt_file->line_info[1]);

    char overlaps_fn[strlen(z_name) + 30];
    sprintf (overlaps_fn, "%s.overlaps", z_name);

    bufprintf (evb, &dups, "##fileformat=GENOZIP-REGIONS-FILE\n"
                           "##contents=Luft coordinates of variants that have different coordinates in the Primary reference and the same coordinates in the Luft reference\n"
                           "##reference=%s\n"
                           "##usage_see_variants=genocat %s --luft --no-header --regions-file %s\n"
                           "##usage_filter_out_variants=genocat %s --luft --regions-file ^%s\n"
                           "##documentation=" WEBSITE_GENOCAT " and " WEBSITE_DVCF "\n"
                           "#CHROM\tSTART\tEND\n", ref_get_filename(gref), z_name, overlaps_fn, z_name, overlaps_fn);

    uint64_t empty_len = dups.len;

    int64_t last_dup=-1;
    for (int64_t index_i=1; index_i < index_len; index_i++) {

        const LineInfo *this_li = &li[index[index_i]];
        const LineInfo *prev_li = &li[index[index_i-1]];

        bool is_dup_luft = (this_li->info.chrom_wi == prev_li->info.chrom_wi) && (this_li->info.start_pos == prev_li->info.start_pos);
        bool is_dup_prim = (this_li->prim_chrom_wi == prev_li->prim_chrom_wi) && (this_li->prim_start_pos == prev_li->prim_start_pos);
        if (is_dup_luft && !is_dup_prim) {
            if (last_dup != index_i-1) // not duplicate duplicate (i.e. 3 or more variants are the same)
                bufprintf (evb, &dups, "%s\t%d\t%d\n", 
                           ctx_snip_from_zf_nodes (ZCTX(DTFZ(luft_chrom)), this_li->info.chrom_wi, 0, 0), 
                           this_li->info.start_pos, this_li->info.start_pos);

            last_dup = index_i;
        }
    }

    if (dups.len > empty_len) { 
        buf_dump_to_file (overlaps_fn, &dups, 1, false, false, false, false);
        
        WARN ("FYI: Genozip detected cases of two or more variants with different Primary coordinates, mapped to the same Luft coordinates.\n"
              "These duplicates (in Luft coordinates) were output to %s, and may be filtered out with:\ngenocat %s --luft --regions-file ^%s\n",
              overlaps_fn, z_name, overlaps_fn);
    }

    buf_free (dups);
}

#ifdef DEBUG
typedef struct { char s[200]; } LcText;
static LcText vcf_linesort_dis_info (LineCmpInfo a)
{
    LcText s;
    sprintf (s.s, "wi=%d start_pos=%d end_pos=%d tie_breaker=%u", a.chrom_wi, a.start_pos, a.end_pos, a.tie_breaker);
    return s;
}
#endif

int vcf_linesort_cmp (LineCmpInfo a, LineCmpInfo b)
{
    if (a.chrom_wi != b.chrom_wi) return a.chrom_wi - b.chrom_wi;

    // case: same chrom, different pos - order by start_pos
    if (a.start_pos < b.start_pos) return -1; // a is smaller (don't use substraction as POS is 64b and return is 32b)
    if (b.start_pos < a.start_pos) return +1; // b is smaller

    // case: same chrom, same start_pos - order by end pos
    if (a.end_pos < b.end_pos) return -1; // a is smaller (don't use substraction as POS is 64b and return is 32b)
    if (b.end_pos < a.end_pos) return +1; // b is smaller

    // same chrom and pos ranges - sort by tie_breaker
    if (a.tie_breaker < b.tie_breaker) return -1; // a is smaller (don't use substraction as its unsigned 32b)
    if (b.tie_breaker < a.tie_breaker) return +1; // b is smaller
    
    return 0; // same chrom, pos, tie-breaker -> undeterministic sorting, hopefully this doesn't happen too often
}

static bool is_luft_sorter;
static int vcf_linesort_line_cmp (const void *a_, const void *b_)
{
    LineInfo *a = B(LineInfo, txt_file->line_info[is_luft_sorter], *(uint32_t *)a_);
    LineInfo *b = B(LineInfo, txt_file->line_info[is_luft_sorter], *(uint32_t *)b_);

    return vcf_linesort_cmp (a->info, b->info);
}

static void vcf_linesort_get_final_index_i (ConstBufferP line_info, BufferP vb_info, ConstBufferP index_buf)
{
    ARRAY (const LineInfo, lorder, *line_info);
    ARRAY (LsVbInfo, v, *vb_info);
    ARRAY (const uint32_t, index, *index_buf);

    uint32_t count_vbs=0;
    uint32_t prev_vb_i=(uint32_t)-1;

    // initialize
    for (VBIType vb_i=0; vb_i < txt_file->num_vbs; vb_i++)
        v[vb_i].final_index_i = 0xffffffff; 

    // scan backwards
    for (int32_t index_i=index_len-1; index_i >= 0 ; index_i--) {
        
        uint32_t i = index[index_i];
        VBIType vb_i = lorder[i].vblock_i - 1; // -1 bc vblock_i is 1-based

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
static uint32_t vcf_linesort_plan_reconstruction (bool is_luft)
{
    ConstBufferP line_info = &txt_file->line_info[is_luft];
    BufferP vb_info        = &txt_file->vb_info[is_luft];
    BufferP index_buf      = &evb->scratch; // an array of uint32_t - indices into line_info
    BufferP recon_plan     = &txt_file->recon_plan;

    // case: no need for a reconstruction plan, as file is sorted
    if (vcf_linesort_is_file_sorted (is_luft)) return 0; 

    // build index sorted by VB - indices into line_info 
    vcf_linesort_create_index (index_buf, is_luft);

    // sort lines
    is_luft_sorter = is_luft; // ugly

    START_TIMER
    qsort (STRb(*index_buf), sizeof (uint32_t), vcf_linesort_line_cmp);
    ARRAY (uint32_t, index, *index_buf);
    COPY_TIMER_EVB (vcf_linesort_compress_qsort);

    // detect cases where distinct variants in Primary got mapped to the same coordinates in Luft
    if (is_luft) line_sorter_detect_duplicates(index_buf);

    // find final entry in index_buf of each VB - needed to calculate how many concurrent threads will be needed for piz
    vcf_linesort_get_final_index_i (line_info, vb_info, index_buf);

    buf_alloc (evb, recon_plan, 0, 100000, ReconPlanItem, 1, "txt_file->recon_plan");

    LineInfo *prev_li = NULL;
    VBIType conc_writing_vbs=0;
    uint32_t max_num_txt_data_bufs=0;
    
    for (uint64_t index_i=0; index_i < index_len; index_i++) {
        
        LineInfo *li = B(LineInfo, *line_info, index[index_i]);

        // verify sort
        #ifdef DEBUG
        ASSERT (!prev_li || vcf_linesort_cmp (prev_li->info, li->info) <= 0,
                "line_info sorting: [%"PRIu64"]={%s} [%"PRIu64"]={%s}", 
                index_i-1, vcf_linesort_dis_info (prev_li->info).s, index_i, vcf_linesort_dis_info (li->info).s);
        #endif
    
        // same vb and consecutive lines - collapse to a single entry
        if (prev_li && 
            li->vblock_i == prev_li->vblock_i &&  
            prev_li->start_line + 1 == li->start_line)
            
            BLST (ReconPlanItem, *recon_plan)->num_lines++;

        // different vb or non-consecutive lines - start new entry
        else {
            buf_alloc (evb, recon_plan, 1, 0, ReconPlanItem, 2, 0);
            BNXT (ReconPlanItem, *recon_plan) = (ReconPlanItem){
                .vb_i       = li->vblock_i,
                .start_line = li->start_line,
                .num_lines  = 1
            };
        }
        
        // account for number of concurrent threads - a thread lives between the first and last access to its VB
        LsVbInfo *v = B(LsVbInfo, *vb_info, li->vblock_i - 1); // -1 as vblock_i is 1-based 
        
        if (!v->accessed) { // case: first access to a VB (in PIZ - we would wait on VB data here)
            v->accessed = true;
            conc_writing_vbs++;

            if (conc_writing_vbs > max_num_txt_data_bufs) 
                max_num_txt_data_bufs = conc_writing_vbs;
        }

        if (v->final_index_i == index_i) { // case: last access to a VB (in PIZ - we would free VB data here)
            conc_writing_vbs--;

            buf_alloc (evb, recon_plan, 1, 0, ReconPlanItem, 2, 0);
            BNXT (ReconPlanItem, txt_file->recon_plan) = (ReconPlanItem){ 
                .vb_i   = li->vblock_i,
                .flavor = PLAN_END_OF_VB,
            }; 
        }

        prev_li = li;
    }

    buf_free (*index_buf);
    
    return max_num_txt_data_bufs;
}

// ZIP main thread
static void vcf_zip_generate_recon_plan_do (bool is_luft)
{
    // create txt_file->recon_plan
    uint32_t conc_writing_vbs = vcf_linesort_plan_reconstruction (is_luft);
    if (!conc_writing_vbs) return; // no need for a recon plan

    recon_plan_compress (conc_writing_vbs, is_luft);
    
    buf_free (txt_file->vb_info[is_luft]);
    buf_free (txt_file->line_info[is_luft]);
}

// ZIP: main thread, called by zip_one_file()
void vcf_zip_generate_recon_plan (void)
{
    START_TIMER;

    // we only need a linesorter recon plan if we're sorting (note: we always sort for DVCF)
    if (!flag.sort && !z_is_dvcf) return;

    // reconstruction plans for primary and luft (might not exist)
    if (txt_file->vb_info[0].len) vcf_zip_generate_recon_plan_do (false);
    if (txt_file->vb_info[1].len) vcf_zip_generate_recon_plan_do (true);

    COPY_TIMER_EVB (generate_recon_plan);
}

