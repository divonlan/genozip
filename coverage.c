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

// output of genocat --show-sex, called from piz_one_file
void coverage_sex_classifier (void)
{    
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

    //printf ("chr1: index=%u cov=%"PRIu64" len=%f\n", index_chr1, *ENT (uint64_t, txt_file->coverage, index_chr1), len_chr1);
    //printf ("chrX: index=%u cov=%"PRIu64" len=%f\n", index_chrX, *ENT (uint64_t, txt_file->coverage, index_chrX), len_chrX);
    //printf ("chrY: index=%u cov=%"PRIu64" len=%f\n", index_chrY, *ENT (uint64_t, txt_file->coverage, index_chrY), len_chrY);
    
    double depth_chr1 = (double)*ENT (uint64_t, txt_file->coverage, index_chr1) / len_chr1;
    double depth_chrX = (double)*ENT (uint64_t, txt_file->coverage, index_chrX) / len_chrX;
    double depth_chrY = (double)*ENT (uint64_t, txt_file->coverage, index_chrY) / len_chrY;

    double ratio_X_Y  =   !depth_chrX ? 0 
                        : !depth_chrY ? 1000 
                        :               depth_chrX / depth_chrY;
    
    double ratio_1_X  =   !depth_chr1 ? 0 
                        : !depth_chrX ? 1000 
                        :               depth_chr1 / depth_chrX;
    
    typedef enum { MALE, BORDERLINE, FEMALE, UNASSIGNED } Sexes;

    Sexes by_1_x = !depth_chr1 || !depth_chrX ? UNASSIGNED
                 : ratio_1_X > 1.75           ? MALE
                 : ratio_1_X < 1.1            ? FEMALE
                 : ratio_1_X < 1.3            ? BORDERLINE // not enough X for a female or XXY - possibly XY/XXY mosaic 
                 :                              UNASSIGNED;
                 
    Sexes by_x_y = !depth_chrX                ? UNASSIGNED
                 : ratio_X_Y < 1.8            ? MALE
                 : ratio_X_Y < 5              ? BORDERLINE // too much X for a male - possible XXY
                 : ratio_X_Y > 9              ? FEMALE
                 :                              UNASSIGNED;
                 
    //  by X/Y ratio →                Male                 Borderline                   Female        Unassigned      ↓ by 1/X ratio ↓
    static char *decision[4][4] = { { "Male",              "Male",                      "Unassigned", "Male",      }, // Male
                                    { "Unassigned",        "Male-XXY or XY/XXY mosaic", "Female",     "Unassigned" }, // Borderline
                                    { "Male-XXY or XXYY",  "Male-XXY",                  "Female",     "Female"     }, // Female
                                    { "Unassigned",        "Unassigned",                "Unassigned", "Unassigned" }  // Unassigned
                                  };

    printf ("%s\t%8.5f\t%8.5f\t%8.5f\t%4.1f\t%4.1f\t%s\n", 
            z_name, depth_chr1, depth_chrX, depth_chrY, ratio_1_X, ratio_X_Y, decision[by_1_x][by_x_y]);
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

        const char *chrom_name = ctx_get_snip_by_word_index (&z_file->contexts[CHROM], i, 0, 0);

        if (flag.show_coverage==1 || strlen (chrom_name) <= 5)
            printf ("%s\t%"PRIu64"\t%"PRIu64"\t%6.2f\n", chrom_name, len, coverage, len ? (double)coverage / (double)len : 0);
    }
}
