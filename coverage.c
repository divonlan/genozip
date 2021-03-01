// ------------------------------------------------------------------
//   coverage.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "coverage.h"
#include "vblock.h"
#include "file.h"

void coverage_initialize (VBlock *vb)
{
    vb->coverage.len = vb->contexts[CHROM].word_list.len;
    buf_alloc_more (vb, &vb->coverage, 0, vb->coverage.len, uint64_t, 1, "coverage");
    buf_zero (&vb->coverage); // zero every VB even if already allocated from prev VB
}

void coverage_add_one_vb (VBlockP vb)
{
    if (!vb->coverage.len) return;

    if (!txt_file->coverage.len)
        buf_copy (evb, &txt_file->coverage, &vb->coverage, sizeof (uint64_t), 0, 0, "txt_file->coverage");

    else
        for (uint64_t i=0; i < vb->coverage.len; i++)
            *ENT (uint64_t, txt_file->coverage, i) += *ENT (uint64_t, vb->coverage, i);
}

static WordIndex coverage_get_chrom_index (char chrom_char)
{
    const char chrom_names[3][5] = { { chrom_char, 0 },
                                   { 'c', 'h', 'r', chrom_char, 0 },
                                   { 'C', 'h', 'r', chrom_char, 0 } };

    for (unsigned i=0; i < sizeof (chrom_names) / sizeof (chrom_names[0]); i++) {
        WordIndex chrom_word_index = ctx_get_word_index_by_snip (&z_file->contexts[CHROM], chrom_names[i]);
        if (chrom_word_index != WORD_INDEX_NONE) return chrom_word_index;
    }

    return WORD_INDEX_NONE;
}

// get coverage and length of primary autosome contigs
static double coverage_get_autosome_depth (WordIndex index_chrX, WordIndex index_chrY)
{
    uint64_t coverage_AS=0, len_AS=0;

    for (uint64_t i=0; i < txt_file->coverage.len; i++) {

        if (i == index_chrX || i == index_chrY) continue; // not autosome

        uint64_t coverage = *ENT (uint64_t, txt_file->coverage, i);
        if (!coverage) continue;

        unsigned cn_len;
        ctx_get_snip_by_word_index (&z_file->contexts[CHROM], i, 0, &cn_len);
        
        PosType len = has_header_contigs ? ENT (RefContig, header_contigs, i)->max_pos
                    :                      ENT (RefContig, loaded_contigs, i)->max_pos;

        if (cn_len > 5 || len <= (1>>20)) continue; // not primary autosome (note: in SAM/BAM compression with REF_INTERNAL the minimum length of a contig 1MB)

        coverage_AS += coverage;
        len_AS      += len;
    }

    ASSERTE0 (len_AS, "Cannot calculate autosome data, because autosome contigs are not loaded");

    return (double)coverage_AS / (double)len_AS;
}

// output of genocat --show-sex, called from piz_one_file
void coverage_sex_classifier (bool is_first_z_file)
{    
    bool is_sam   = (z_file->data_type == DT_SAM);
    bool is_fastq = (z_file->data_type == DT_FASTQ);

    WordIndex index_chr1 = coverage_get_chrom_index ('1');
    WordIndex index_chrX = coverage_get_chrom_index ('X');
    WordIndex index_chrY = coverage_get_chrom_index ('Y');

    double len_chr1 = index_chr1 == WORD_INDEX_NONE ? 1
                    : has_header_contigs            ? ENT (RefContig, header_contigs, index_chr1)->max_pos
                    :                                 ENT (RefContig, loaded_contigs, index_chr1)->max_pos;

    double len_chrX = index_chrX == WORD_INDEX_NONE ? 1
                    : has_header_contigs            ? ENT (RefContig, header_contigs, index_chrX)->max_pos
                    :                                 ENT (RefContig, loaded_contigs, index_chrX)->max_pos;

    double len_chrY = index_chrY == WORD_INDEX_NONE ? 1
                    : has_header_contigs            ? ENT (RefContig, header_contigs, index_chrY)->max_pos
                    :                                 ENT (RefContig, loaded_contigs, index_chrY)->max_pos;


    //printf ("chrY: index=%u cov=%"PRIu64" len=%f\n", index_chrY, *ENT (uint64_t, txt_file->coverage, index_chrY), len_chrY);
    
    // in SAM, it is sufficient to look at chr1 as it is already mapped. this allows as to run with --regions, 10 times faster
    // in FASTQ, we rely on our own approximate aligner, so we compare against all autosomes
    double depth_AS   = is_sam ? (len_chr1 ? (double)*ENT (uint64_t, txt_file->coverage, index_chr1) / len_chr1 : 0)
                               : coverage_get_autosome_depth (index_chrX, index_chrY);

    double depth_chrX = len_chrX ? (double)*ENT (uint64_t, txt_file->coverage, index_chrX) / len_chrX : 0;
    double depth_chrY = len_chrY ? (double)*ENT (uint64_t, txt_file->coverage, index_chrY) / len_chrY : 0;

    // correct for Genozip Aligner's bias in favour X and Y (vs autosomes) in humans (in FASTQ)
    double correction = is_fastq ? 1.333 : 1;

    double ratio_AS_X = !depth_AS   ? 0 
                      : !depth_chrX ? 1000 
                      :               correction * depth_AS / depth_chrX;
    
    double ratio_X_Y  = !depth_chrX ? 0 
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

    if (is_first_z_file) 
        printf ("File\t%s_Depth\tX_Depth\tY_Depth\t%s\tX_Depth/Y_Depth\tSex\n", 
                is_sam ? "1" : "AS", 
                is_sam ? "1_Depth/X_Depth" : "AS_Depth/X_Depth (corrected)"); 
        
    printf ("%s\t%8.5f\t%8.5f\t%8.5f\t%4.1f\t%4.1f\t%s\n", 
            z_name, depth_AS, depth_chrX, depth_chrY, ratio_AS_X, ratio_X_Y, 
            is_sam ? sam_call[by_as_x][by_x_y] : fq_call[by_as_x][by_x_y]);

    fflush (stdout); // in case output is redirected
}

// output of genocat --show-coverage, called from piz_one_file
void coverage_show_coverage (void)
{
    printf ("Contig\tLN\tCoverage\tDepth\n");
        
    for (uint64_t i=0; i < txt_file->coverage.len; i++) {
        uint64_t coverage = *ENT (uint64_t, txt_file->coverage, i);
        if (!coverage) continue;

        PosType len = has_header_contigs ? ENT (RefContig, header_contigs, i)->max_pos
                    :                      ENT (RefContig, loaded_contigs, i)->max_pos;

        unsigned cn_len;
        const char *chrom_name = ctx_get_snip_by_word_index (&z_file->contexts[CHROM], i, 0, &cn_len);

        if (flag.show_coverage==1 || cn_len <= 5) 
            printf ("%s\t%"PRIu64"\t%"PRIu64"\t%6.2f\n", chrom_name, len, coverage, len ? (double)coverage / (double)len : 0);
    }

    fflush (stdout); // in case output is redirected
}
