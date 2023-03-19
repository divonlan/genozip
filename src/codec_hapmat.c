// ------------------------------------------------------------------
//   codec_hapmat.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// *** RETIRED CODEC. USED JUST FOR DECOMPRESSING OLD VCF FILES (retired officially in v13, but hasn't been in use since some earlier version) ***

#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"
#include "strings.h"
#include "file.h"
#include "piz.h"
#include "endianness.h"
#include "context.h"
#include "dict_id.h"
#include "reconstruct.h"
#include "vcf_private.h"
#include "data_types.h"

//--------------
// PIZ side
//--------------

// PIZ: for each haplotype column, retrieve its it address in the haplotype sections. Note that since the haplotype sections are
// transposed, each column will be a row, or a contiguous array, in the section data. This function returns an array
// of pointers, each pointer being a beginning of column data within the section array
void codec_hapmat_piz_calculate_columns (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    vb->ht_matrix_ctx = ECTX (_PBWT_HT_MATRIX);

    uint32_t ht_per_line = vb->ht_per_line = CTX(FORMAT_GT_HT)->local.len32 / vb->lines.len32;

    vb->hapmat_one_array.len32 = ht_per_line + 7; // +7 because depermuting_loop works on a word (32/64 bit) boundary
    buf_alloc (vb, &vb->hapmat_one_array, 0, vb->hapmat_one_array.len, char, 1, "hapmat_one_array");
    buf_alloc (vb, &vb->hapmat_columns_data, 0, vb->hapmat_one_array.len, char*, 1, "hapmat_columns_data"); // realloc for exact size (+15 is padding for 64b operations)

    // each entry is a pointer to the beginning of haplotype column located in vb->haplotype_sections_data
    // note: haplotype columns are permuted only within their own sample block
    ARRAY (rom , hapmat_columns_data, vb->hapmat_columns_data); 

    ContextP perm_ctx = ECTX (_PBWT_GT_HT_INDEX);
    ARRAY (const unsigned, permutatation_index, perm_ctx->local);
    
    ASSERT (permutatation_index, "%s.local is empty", perm_ctx->tag_name);

    // provide 7 extra zero-columns for the convenience of the permuting loop (supporting 64bit assignments)
    // note: txt_file->max_lines_per_vb will be zero if genozip file was created by redirecting output
    buf_alloc (vb, &vb->hapmat_column_of_zeros, 0, MAX_(txt_file->max_lines_per_vb, vb->lines.len), char, 1, "hapmat_column_of_zeros");
    buf_zero (&vb->hapmat_column_of_zeros);

    for (uint32_t ht_i = 0; ht_i < ht_per_line; ht_i++) 
        hapmat_columns_data[ht_i] = Bc (vb->ht_matrix_ctx->local, permutatation_index[ht_i] * vb->lines.len);

    for (unsigned ht_i=ht_per_line; ht_i < ht_per_line + 7; ht_i++)
        hapmat_columns_data[ht_i] = vb->hapmat_column_of_zeros.data;
}

// PIZ: build haplotype for a line - reversing the permutation and the transposal.
static inline void codec_hapmat_piz_get_one_line (VBlockVCFP vb)
{
    START_TIMER;

    if (flag.samples) buf_zero (&vb->hapmat_one_array); // if we're not filling in all samples, initialize to 0;

    ARRAY (rom , hapmat_columns_data, vb->hapmat_columns_data); 
    uint32_t ht_i_after = vb->ht_per_line; // automatic variable - faster
    uint64_t *next = B1ST64 (vb->hapmat_one_array);
    // LineIType vb_line_i = vb->line_i - vb->first_line;

    // this loop can consume up to 25-50% of the entire decompress compute time (tested with 1KGP data)
    // note: we do memory assignment 64 bit at time (its about 10% faster than byte-by-byte)
    
    for (uint32_t ht_i=0; ht_i < ht_i_after; ht_i += 8) {

#ifdef __LITTLE_ENDIAN__
        *(next++) = ((uint64_t)(uint8_t)hapmat_columns_data[ht_i    ][vb->line_i]      ) |  // this is LITTLE ENDIAN order
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 1][vb->line_i] << 8 ) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 2][vb->line_i] << 16) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 3][vb->line_i] << 24) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 4][vb->line_i] << 32) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 5][vb->line_i] << 40) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 6][vb->line_i] << 48) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 7][vb->line_i] << 56) ;  // no worries if ht_per_line is not a multiple of 4 - we have extra columns of zero
#else
        *(next++) = ((uint64_t)(uint8_t)hapmat_columns_data[ht_i    ][vb->line_i] << 56) |  // this is BIG ENDIAN order
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 1][vb->line_i] << 48) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 2][vb->line_i] << 40) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 3][vb->line_i] << 32) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 4][vb->line_i] << 24) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 5][vb->line_i] << 16) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 6][vb->line_i] << 8 ) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 7][vb->line_i]      ) ;  
#endif
    }

    if (flag.show_alleles)
        iprintf ("Line %-2d : %.*s\n", vb->line_i, vb->hapmat_one_array.len32, vb->hapmat_one_array.data);

    COPY_TIMER (codec_hapmat_piz_get_one_line);
}

// Explanation of the reconstruction process of VCF haplotype matrix data compressed with the HT codec:
// 1) The GT_HT and GT_HT_INDEX sections as decompressed normally using the sub-codecs
// 2) When reconstructing, this function codec_hapmat_reconstruct is called for every haplotype value needed (eg the "1" in "1/0"),
//    by PIZ as the LT_CODEC reconstructor for CODEC_HAPM. It takes data from GT_HT, consulting GT_HT_INDEX
//    to find the correct column in the matrix, and skips '*'s (missing haplotypes due to mixed ploidy, missing samples,
//    or lines missing GT in FORMAT) until it finds a valid haplotype value.
CODEC_RECONSTRUCT (codec_hapmat_reconstruct)
{
    // get one row of the haplotype matrix for this line into vb->hapmat_one_array if we don't have it already
    if (VB_VCF->hapmat_one_array.param != vb->line_i + 1) { // we store the line_i in param (+1 to distguish from "not set" value of 0)
        codec_hapmat_piz_get_one_line (VB_VCF);
        VB_VCF->hapmat_one_array.len = 0; // length of data consumed
        VB_VCF->hapmat_one_array.param = vb->line_i + 1;
    }

    // find next allele - skipping unused spots ('*')
    uint8_t ht = '*';
    while (ht == '*' && VB_VCF->hapmat_one_array.len32 < vb->ht_per_line)
        ht = *B8(VB_VCF->hapmat_one_array, VB_VCF->hapmat_one_array.len++);

    if (vb->drop_curr_line) return;

    switch (ht) {
        case '.': case '0' ... '9':
            RECONSTRUCT1 (ht); break;

        case 58 ... 147:  // allele 10 to 99 (ascii 58 to 147)
            RECONSTRUCT_INT (ht - '0'); break;
        
        case '-': // ploidy padding (starting 10.0.2) - appears as the 2nd+ HT - unlike '*', these are counted in GT.repeats
            Ltxt--; // remove previous phase character;
            break;

        // "%|%" is used instead of "./." in case the previous sample had a | phase. This is meant to reduce entroy of GT.b250
        // in the case we have phased data, but missing samples appear as "./.".
        case '%': 
            if (*BLSTtxt == '|' || *BLSTtxt == '/') { // second % in "%|%" and act on the 2nd
                Ltxt -= 1;  // remove | or /
                RECONSTRUCT ("/.", 2);
            }
            else
                RECONSTRUCT1 ('.'); // first % in "%|%" 
            break;

        default: 
            ASSPIZ (false, "Invalid character found in decompressed HT array: '%c' (ASCII %u)",  ht, ht);
    }
}

