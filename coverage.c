// ------------------------------------------------------------------
//   coverage.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "coverage.h"
#include "vblock.h"
#include "file.h"
#include "strings.h"
#include "txtfile.h"
#include "contigs.h"
#include "sam_private.h"

void coverage_initialize (VBlock *vb)
{
    vb->coverage.len = CTX(CHROM)->word_list.len + NUM_COVER_TYPES;
    buf_alloc (vb, &vb->coverage, 0, vb->coverage.len, uint64_t, 1, "coverage");
    buf_zero (&vb->coverage); // zero every VB even if already allocated from prev VB

    vb->read_count.len = CTX(CHROM)->word_list.len + NUM_COVER_TYPES;
    buf_alloc (vb, &vb->read_count, 0, vb->read_count.len, uint64_t, 1, "read_count");
    buf_zero (&vb->read_count); // zero every VB even if already allocated from prev VB
    
    vb->unmapped_read_count.len = CTX(CHROM)->word_list.len;
    buf_alloc (vb, &vb->unmapped_read_count, 0, vb->unmapped_read_count.len, uint64_t, 1, "unmapped_read_count");
    buf_zero (&vb->unmapped_read_count); // zero every VB even if already allocated from prev VB
}

void coverage_add_one_vb (VBlockP vb)
{
    buf_add_vector (evb, uint64_t, txt_file->coverage, vb->coverage, "txt_file->coverage");
    buf_add_vector (evb, uint64_t, txt_file->read_count, vb->read_count, "txt_file->read_count");
    buf_add_vector (evb, uint64_t, txt_file->unmapped_read_count, vb->unmapped_read_count, "txt_file->unmapped_read_count");
}

static WordIndex coverage_get_chrom_index (char chrom_char)
{
    const char chrom_names[3][5] = { { chrom_char, 0 },
                                   { 'c', 'h', 'r', chrom_char, 0 },
                                   { 'C', 'h', 'r', chrom_char, 0 } };

    for (unsigned i=0; i < sizeof (chrom_names) / sizeof (chrom_names[0]); i++) {
        WordIndex chrom_word_index = ctx_get_word_index_by_snip (ZCTX(CHROM), chrom_names[i]);
        if (chrom_word_index != WORD_INDEX_NONE) return chrom_word_index;
    }

    return WORD_INDEX_NONE;
}

static ConstContigPkgP coverage_get_contigs (const char *option_name)
{
    ConstContigPkgP contigs = ((Z_DT(DT_SAM) || Z_DT(DT_BAM)) && sam_hdr_contigs) ? sam_hdr_contigs : ref_get_ctgs (gref);

    ASSINP (!Z_DT(DT_FASTQ) || contigs->contigs.len, "%s: %s for FASTQ only works on files compressed with a reference", z_name, option_name);

    return contigs;
}

// get coverage and length of primary autosome contigs
static float coverage_get_autosome_depth (ConstContigPkgP contigs, WordIndex index_chrX, WordIndex index_chrY)
{
    uint64_t coverage_AS=0, len_AS=0;
    
    for (uint64_t i=0; i < txt_file->coverage.len - NUM_COVER_TYPES; i++) {

        if (i == index_chrX || i == index_chrY) continue; // not autosome

        uint64_t coverage = *ENT (uint64_t, txt_file->coverage, i);
        if (!coverage) continue;

        unsigned cn_len;
        ctx_get_snip_by_word_index (ZCTX(CHROM), i, 0, &cn_len);
        
        PosType len = contigs_get_LN (contigs, i);

        if (cn_len > 5 || len <= (1>>20)) continue; // not primary autosome (note: in SAM/BAM compression with REF_INTERNAL the minimum length of a contig 1MB)

        coverage_AS += coverage;
        len_AS      += len;
    }

    ASSERT0 (len_AS, "Cannot calculate autosome data, because autosome contigs are not loaded");

    return (float)coverage_AS / (float)len_AS;
}

// output of genocat --sex, called from piz_one_txt_file
void coverage_sex_classifier (bool is_first_z_file)
{    
    ConstContigPkgP contigs = coverage_get_contigs ("--sex");

    bool is_sam   = (Z_DT(DT_SAM));
    bool is_fastq = (Z_DT(DT_FASTQ));

    WordIndex index_chr1 = coverage_get_chrom_index ('1');
    WordIndex index_chrX = coverage_get_chrom_index ('X');
    WordIndex index_chrY = coverage_get_chrom_index ('Y');

    float len_chr1 = index_chr1 == WORD_INDEX_NONE ? 1 : contigs_get_LN (contigs, index_chr1);
    float len_chrX = index_chrX == WORD_INDEX_NONE ? 1 : contigs_get_LN (contigs, index_chrX);
    float len_chrY = index_chrY == WORD_INDEX_NONE ? 1 : contigs_get_LN (contigs, index_chrY);

    ARRAY (uint64_t, coverage, txt_file->coverage);

    if (!coverage_len || !coverage[index_chr1] || !coverage[index_chrX]) {
        WARN ("%s: %s: --sex doesn't work for this file, because it is missing some of the autosomal or sex chromosomes (chr1_bases=%"PRIu64" chrX_bases=%"PRIu64")", 
               global_cmd, z_name, coverage_len ? coverage[index_chr1] : 0, coverage_len ? coverage[index_chrX] : 0);
        return;
    }

    //iprintf ("chrY: index=%u cov=%"PRIu64" len=%f\n", index_chrY, *ENT (uint64_t, txt_file->coverage, index_chrY), len_chrY);
    
    // in SAM, it is sufficient to look at chr1 as it is already mapped. this allows as to run with --regions, 10 times faster
    // in FASTQ, we rely on our own approximate aligner, so we compare against all autosomes
    float depth_AS   = is_sam ? (len_chr1 ? (float)coverage[index_chr1] / len_chr1 : 0)
                               : coverage_get_autosome_depth (contigs, index_chrX, index_chrY);

    float depth_chrX = len_chrX ? (float)coverage[index_chrX] / len_chrX : 0;
    float depth_chrY = len_chrY ? (float)coverage[index_chrY] / len_chrY : 0;

    // correct for Genozip Aligner's bias in favour X and Y (vs autosomes) in humans (in FASTQ)
    float correction = is_fastq ? 1.333 : 1;

    float ratio_AS_X = !depth_AS   ? 0 
                      : !depth_chrX ? 1000 
                      :               correction * depth_AS / depth_chrX;
    
    float ratio_X_Y  = !depth_chrX ? 0 
                      : !depth_chrY ? 1000 
                      :               depth_chrX / depth_chrY;
    
    typedef enum { XTEST_X, XTEST_XX_MINUS, XTEST_XX, XTEST_NA } XTest;

    XTest by_as_x = !depth_AS || !depth_chrX ? XTEST_NA
                  : ratio_AS_X > 1.75        ? XTEST_X
                  : ratio_AS_X < 1.1         ? XTEST_XX
                  : ratio_AS_X < 1.3         ? XTEST_XX_MINUS // not enough X for a female or XXY - possibly XY/XXY mosaic 
                  :                            XTEST_NA;
                 
    typedef enum { XYTEST_HAS_Y, XYTEST_HAS_Y_MINUS, XYTEST_NO_Y, XYTEST_NA } XYTest;

    XYTest by_x_y = !depth_chrX              ? XYTEST_NA
                  : ratio_X_Y < 1.8          ? XYTEST_HAS_Y
                  : ratio_X_Y < 5            ? XYTEST_HAS_Y_MINUS // has Y, but too little vs X for a male - possible XXY
                  : ratio_X_Y > 9            ? XYTEST_NO_Y
                  :                            XYTEST_NA;
                 
    //  by X/Y ratio →                XYTEST_HAS_Y         XYTEST_HAS_Y_MINUS           XYTEST_NO_Y   XYTEST_NA          ↓ by AS/X ratio ↓
    static char *sam_call[4][4] = { { "Male",              "Male",                      "Unassigned", "Male",      }, // XTEST_X
                                    { "Unassigned",        "Male-XXY or XY/XXY mosaic", "Female",     "Unassigned" }, // XTEST_XX_MINUS
                                    { "Male-XXY or XXYY",  "Male-XXY",                  "Female",     "Female"     }, // XTEST_XX
                                    { "Unassigned",        "Unassigned",                "Unassigned", "Unassigned" }  // XTEST_NA
                                  };

    //  by X/Y ratio →                XYTEST_HAS_Y         XYTEST_HAS_Y_MINUS           XYTEST_NO_Y   XYTEST_NA          ↓ by AS/X ratio ↓
    static char *fq_call[4][4] =  { { "Male",              "Male",                      "Unassigned", "Male",      }, // XTEST_X
                                    { "Unassigned",        "Unassigned",                "Female",     "Female"     }, // XTEST_XX_MINUS
                                    { "Male-XXY or XXYY",  "Unassigned",                "Female",     "Female"     }, // XTEST_XX
                                    { "Unassigned",        "Unassigned",                "Unassigned", "Unassigned" }  // XTEST_NA
                                  };

    if (!flag.multiple_files) 
        iprintf (is_info_stream_terminal ? "\n--sex for: %s\n%-10s  %-*s  %-6s  %-6s  %-6s  %-4s  %-4s\n" 
                                          : "\n--sex for: %s\n%s\t%*s\t%s\t%s\t%s\t%s\t%s\n",
                z_name, "Sex", flag.longest_filename, "File", is_sam ? "DP_1" : "DP_AS", "DP_X", "DP_Y", 
                is_sam ? "1/X" : "AS/X", "X/Y");
                
    iprintf (is_info_stream_terminal ? "%-10s  %-*s  %-6.3f  %-6.3f  %-6.3f  %-4.1f  %-4.1f\n" : "%s\t%*s\t%f\t%f\t%f\t%f\t%f\n",  
            is_sam ? sam_call[by_as_x][by_x_y] : fq_call[by_as_x][by_x_y], 
            flag.longest_filename, z_name, 
            depth_AS, depth_chrX, depth_chrY, ratio_AS_X, ratio_X_Y);

    fflush (info_stream); // in case output is redirected
}

// output of genocat --coverage, called from piz_one_txt_file
void coverage_show_coverage (void)
{
    ASSERTNOTEMPTY (txt_file->coverage);

    ConstContigPkgP contigs = coverage_get_contigs ("--coverage");

    unsigned chr_width = (flag.show_coverage == COV_ALL ? 26 : 13);

    if (flag.show_coverage != COV_ONE)
        iprintf (is_info_stream_terminal ? "\n--coverage for: %s\n%-*s  %-8s  %-11s  %-10s  %-5s   %s\n" 
                                         : "\n--coverage for: %s\n%*s\t%s\t%s\t%s\t%s\t%s\n", 
                 z_name, chr_width, "Contig", "LN", "Reads", "Bases", "Bases", "Depth");

    txt_file->coverage.len -= NUM_COVER_TYPES; // real contigs only
    ARRAY (uint64_t, coverage, txt_file->coverage);
    ARRAY (uint64_t, read_count, txt_file->read_count);

    uint64_t *coverage_special   = &coverage[coverage_len];
    uint64_t *read_count_special = &read_count[coverage_len];

    // calculate CVR_TOTAL - all bases in the file
    for (uint64_t i=0; i < coverage_len + NUM_COVER_TYPES - 1; i++) {
        coverage_special  [CVR_TOTAL] += coverage[i];
        read_count_special[CVR_TOTAL] += read_count[i];
    }

    // calculate CVR_PRIMARY - all bases in the contig (i.e. excluding unmapped, secondary, supplementary, duplicate, soft-clipped)
    for (uint64_t i=0; i < coverage_len; i++) {
        coverage_special  [CVR_ALL_CONTIGS] += coverage[i];
        read_count_special[CVR_ALL_CONTIGS] += read_count[i];
    }

    PosType genome_nbases = 0;

    for (uint64_t i=0; i < coverage_len; i++) {
        if (!coverage[i]) continue;

        unsigned cn_len;
        const char *chrom_name = ctx_get_snip_by_word_index (ZCTX(CHROM), i, 0, &cn_len);

        PosType len = (cn_len==1 && chrom_name[0]=='*') ? 0 : contigs_get_LN (contigs, i);

        if (flag.show_coverage == COV_ALL || (flag.show_coverage == COV_CHROM && cn_len <= 5))
            iprintf (is_info_stream_terminal ? "%-*s  %-8s  %-11s  %-10s  %-4.1f%%  %6.2f\n" : "%*s\t%s\t%s\t%s\t%4.1f\t%6.2f\n", 
                     chr_width, chrom_name, str_bases(len).s, str_uint_commas (read_count[i]).s, str_bases(coverage[i]).s, 
                     100.0 * (float)coverage[i] / (float)coverage_special[CVR_TOTAL], 
                     len ? (float)coverage[i] / (float)len : 0);

        else {
            coverage_special[CVR_OTHER_CONTIGS]   += coverage[i]; // other non-chromosome contigs
            read_count_special[CVR_OTHER_CONTIGS] += read_count[i];
        }

        if (cn_len <= 5) genome_nbases += len;
    }

    char all_coverage[7] = "0";
    if (genome_nbases) // avoid division by 0
        sprintf (all_coverage, "%*.2f", (int)sizeof(all_coverage)-1, (float)coverage_special[CVR_ALL_CONTIGS] / (float)genome_nbases);

    if (flag.show_coverage == COV_ONE) 
        iprintf ("%s:\t%s\n", z_name, all_coverage);

    else {
        char *cvr_names[NUM_COVER_TYPES] = { "Other contigs", "All contigs", "Soft clip", "Unmapped", "Secondary", "Supplementary", "Failed filters", "Duplicate", "TOTAL"};

        for (uint64_t i=0; i < NUM_COVER_TYPES; i++) {

            if (coverage_special[i])
                iprintf (is_info_stream_terminal ? "%-*s  %-8s  %-11s  %-10s  %-4.1f%%  %s\n" : "%*s\t%s\t%s\t%s\t%4.1f%%\t%s\n", 
                        chr_width, cvr_names[i], "",
                        (i == CVR_SOFT_CLIP ? "" : str_uint_commas (read_count_special[i]).s),
                        str_bases(coverage_special[i]).s,
                        100.0 * (float)coverage_special[i] / (float)coverage_special[CVR_TOTAL],
                        i == CVR_ALL_CONTIGS ? all_coverage : "");
            
            if (i == CVR_OTHER_CONTIGS && is_info_stream_terminal)
                iprint0 ("-----\n");
        }
    }
}

// output of genocat --idxstats - designed to identical to samtools idxstats
void coverage_show_idxstats (void)
{
    ASSERTNOTEMPTY (txt_file->coverage);

    ConstContigPkgP contigs = coverage_get_contigs ("--idxstats");

    txt_file->read_count.len -= NUM_COVER_TYPES; // real contigs only
    ARRAY (uint64_t, read_count, txt_file->read_count);
    ARRAY (uint64_t, unmapped_read_count, txt_file->unmapped_read_count);

    for (uint64_t i=0; i < read_count_len; i++) {
        unsigned cn_len;
        const char *chrom_name = ctx_get_snip_by_word_index (ZCTX(CHROM), i, 0, &cn_len);

        PosType len = (cn_len==1 && chrom_name[0]=='*') ? 0 : contigs_get_LN (contigs, i);

        iprintf ("%.*s\t%"PRIu64"\t%"PRId64"\t%"PRIu64"\n", cn_len, chrom_name, len, read_count[i], unmapped_read_count[i]);
    }

    // FASTQ (but not SAM) unmapped reads
    uint64_t unmapped_unknown = AFTERENT(uint64_t, txt_file->read_count)[CVR_UNMAPPED]; 
    if (unmapped_unknown)
        iprintf ("*\t0\t0\t%"PRIu64"\n", unmapped_unknown);

    fflush (info_stream); // in case output is redirected
}
