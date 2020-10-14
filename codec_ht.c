// ------------------------------------------------------------------
//   codec_ht.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"
#include "strings.h"
#include "file.h"
#include "endianness.h"
#include "context.h"
#include "dict_id.h"
#include "piz.h"

//--------------
// ZIP side
//--------------

typedef struct {
    int num_alt_alleles;
    uint32_t index_in_original_line;
    uint32_t index_in_sorted_line;
} HaploTypeSortHelperIndex;

void codec_ht_comp_init (VBlock *vb)
{
    vb->ht_ctx         = mtf_get_ctx (vb, dict_id_FORMAT_GT_HT);
    vb->ht_ctx->ltype  = LT_HT;
    vb->ht_ctx->lcodec = CODEC_HT;

    vb->ht_index_ctx         = mtf_get_ctx (vb, dict_id_FORMAT_GT_HT_INDEX);
    vb->ht_index_ctx->ltype  = LT_UINT32;
    vb->ht_index_ctx->lcodec = codec_args[CODEC_HT].sub_codec2; // BSC is 4-5% better than LZMA, only slightly slower
}

static int sort_by_alt_allele_comparator(const void *p, const void *q)  
{ 
    int l = ((HaploTypeSortHelperIndex *)p)->num_alt_alleles; 
    int r = ((HaploTypeSortHelperIndex *)q)->num_alt_alleles;  
    return (l - r); 
}

static int sort_by_original_index_comparator(const void *p, const void *q)  
{ 
    int l = ((HaploTypeSortHelperIndex *)p)->index_in_original_line; 
    int r = ((HaploTypeSortHelperIndex *)q)->index_in_original_line;  
    return (l - r); 
}

static HaploTypeSortHelperIndex *codec_ht_count_alt_alleles (VBlock *vb)
{
    START_TIMER; 

    buf_alloc (vb, &vb->helper_index_buf, vb->num_haplotypes_per_line * sizeof(HaploTypeSortHelperIndex), 0,
               "helper_index_buf", vb->vblock_i);
    buf_zero (&vb->helper_index_buf);
    ARRAY (HaploTypeSortHelperIndex, helper_index, vb->helper_index_buf);

    // build index array 
    for (uint32_t ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i++) 
        helper_index[ht_i].index_in_original_line = ht_i;

    for (uint32_t line_i=0; line_i < (uint32_t)vb->lines.len; line_i++) {
        for (uint32_t ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i++) {

            // we count as alt alleles : 1 - 99 (ascii 49 to 147)
            //             ref alleles : 0 . (unknown) * (missing)
            uint8_t one_ht = *ENT (uint8_t, vb->ht_ctx->local, line_i * vb->num_haplotypes_per_line + ht_i);
            if (one_ht >= '1')
                helper_index[ht_i].num_alt_alleles++;
        }
    }
    COPY_TIMER (codec_ht_count_alt_alleles);

    return helper_index;
}

static void codec_ht_compress_one_array (VBlockP vb, uint32_t ht_i, 
                                         char **line_data_1, uint32_t *line_data_len_1,
                                         char **line_data_2, uint32_t *line_data_len_2)
{
    uint32_t orig_col_i = ENT (HaploTypeSortHelperIndex, vb->helper_index_buf, ht_i)->index_in_original_line; 
    const uint32_t increment = vb->num_haplotypes_per_line; // automatic var that can be optimized to a register

    uint8_t *ht_data_ptr = ENT (uint8_t, vb->ht_ctx->local, orig_col_i);
    ARRAY (uint8_t, column, vb->ht_one_array);

    for (uint32_t line_i=0; line_i < vb->ht_one_array.len; line_i++, ht_data_ptr += increment) 
        column[line_i] = *ht_data_ptr;

    *line_data_1 = (char *)column;
    *line_data_len_1 = vb->ht_one_array.len;

    *line_data_2 = NULL;
    *line_data_len_2 = 0;

    if (flag_show_alleles)
        printf ("Col %-2u : %.*s\n", orig_col_i, (int)vb->ht_one_array.len, column);
}

// sort haplogroups by alt allele count within the variant group, create an index for it, and split
// it to sample groups. for each sample a haplotype is just a string of 1 and 0 etc (could be other alleles too)
// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
bool codec_ht_compress (VBlock *vb, 
                       Codec *codec,             // out
                       const char *uncompressed, // option 1 - compress contiguous data
                       uint32_t *uncompressed_len, 
                       LocalGetLineCB callback,  // option 2 - not supported
                       char *compressed, uint32_t *compressed_len /* in/out */, 
                       bool soft_fail)
{
    START_TIMER;
    
    ASSERT0 (uncompressed && !callback, "Error in codec_ht_compress: callback option not supported");

    HaploTypeSortHelperIndex *helper_index = codec_ht_count_alt_alleles (vb);
    
    // sort the helper index by number of alt alleles
    qsort (helper_index, vb->num_haplotypes_per_line, sizeof (HaploTypeSortHelperIndex), sort_by_alt_allele_comparator);

    vb->ht_one_array.len = vb->lines.len;
    buf_alloc (vb, &vb->ht_one_array, vb->ht_one_array.len, 1, "ht_one_array", 0);
    
    if (flag_show_alleles) 
        printf ("\nAfter transpose and sorting:\n");

    // compress the matrix one column at a time, by the order of helper index
    uint64_t save_lines_len = vb->lines.len;
    vb->lines.len = vb->num_haplotypes_per_line; // temporarily set vb->lines.len to number of columns, as this is the number of time the callback will be called
    
    // compress as a normal sub-codec section. all piz-side logic will happen during reconstruction.
    *codec = codec_args[CODEC_HT].sub_codec1;
    Compressor compress = codec_args[*codec].compress;

    PAUSE_TIMER; //  don't include sub-codec compressor - it accounts for itself

    bool success = compress (vb, codec, 0, uncompressed_len, codec_ht_compress_one_array, compressed, compressed_len, soft_fail);
    
    RESUME_TIMER (compressor_ht);

    vb->lines.len = save_lines_len;

    if (!success) return false; // soft fail (if it was hard fail, compress() already failed)

    // final step - build the reverse index that will allow access by the original index to the sorted array
    // this will be included in the genozip file
    for (uint32_t ht_i=0; ht_i < vb->num_haplotypes_per_line ; ht_i++)
        helper_index[ht_i].index_in_sorted_line = ht_i;

    // sort array back to its original order
    qsort (helper_index, vb->num_haplotypes_per_line, sizeof (HaploTypeSortHelperIndex), sort_by_original_index_comparator);

    // create a permutation index for the vblock
    // we populate the ht_index_ctx local, and it will be written after us, as the context is create after the ht_ctx 
    buf_alloc (vb, &vb->ht_index_ctx->local, vb->num_haplotypes_per_line * sizeof(uint32_t), 
               0, "context->local", vb->ht_index_ctx->did_i);
    vb->ht_index_ctx->local.len = vb->num_haplotypes_per_line;
    
    ARRAY (uint32_t, hp_index, vb->ht_index_ctx->local);
    for (uint32_t ht_i=0; ht_i < vb->num_haplotypes_per_line ; ht_i++)
        hp_index[ht_i] = BGEN32 (helper_index[ht_i].index_in_sorted_line);

    COPY_TIMER (compressor_ht);

    return true;
}

//--------------
// PIZ side
//--------------

// PIZ: for each haplotype column, retrieve its it address in the haplotype sections. Note that since the haplotype sections are
// transposed, each column will be a row, or a contiguous array, in the section data. This function returns an array
// of pointers, each pointer being a beginning of column data within the section array
void codec_ht_piz_calculate_columns (VBlock *vb)
{
    uint32_t num_haplotypes_per_line = vb->num_haplotypes_per_line = (uint32_t)(vb->ht_ctx->local.len / vb->lines.len);

    vb->ht_one_array.len = num_haplotypes_per_line + 7; // +7 because depermuting_loop works on a word (32/64 bit) boundary
    buf_alloc (vb, &vb->ht_one_array, vb->ht_one_array.len, 1, "ht_one_array", vb->vblock_i);
    buf_alloc (vb, &vb->ht_columns_data, sizeof (char *) * vb->ht_one_array.len, 1, "ht_columns_data", 0); // realloc for exact size (+15 is padding for 64b operations)

    // each entry is a pointer to the beginning of haplotype column located in vb->haplotype_sections_data
    // note: haplotype columns are permuted only within their own sample block
    ARRAY (const char *, ht_columns_data, vb->ht_columns_data); 
    ARRAY (const unsigned, permutatation_index, vb->ht_index_ctx->local);
    
    ASSERT0 (permutatation_index, "Error in codec_ht_piz_calculate_columns: ht_index_ctx.local is empty");

    // provide 7 extra zero-columns for the convenience of the permuting loop (supporting 64bit assignments)
    // note: txt_file->max_lines_per_vb will be zero if genozip file was created by redirecting output
    buf_alloc (vb, &vb->column_of_zeros, MAX (txt_file->max_lines_per_vb, vb->lines.len), 1, "column_of_zeros", 0);
    buf_zero (&vb->column_of_zeros);

    for (uint32_t ht_i = 0; ht_i < num_haplotypes_per_line; ht_i++) 
        ht_columns_data[ht_i] = ENT (char, vb->ht_ctx->local, permutatation_index[ht_i] * vb->lines.len);

    for (unsigned ht_i=num_haplotypes_per_line; ht_i < num_haplotypes_per_line + 7; ht_i++)
        ht_columns_data[ht_i] = vb->column_of_zeros.data;
}

// PIZ: build haplotype for a line - reversing the permutation and the transposal.
static inline void codec_ht_piz_get_one_line (VBlock *vb)
{
    START_TIMER;

    if (flag_samples) buf_zero (&vb->ht_one_array); // if we're not filling in all samples, initialize to 0;

    ARRAY (const char *, ht_columns_data, vb->ht_columns_data); 
    uint32_t ht_i_after = vb->num_haplotypes_per_line; // automatic variable - faster
    uint64_t *next = FIRSTENT (uint64_t, vb->ht_one_array);
    uint32_t vb_line_i = vb->line_i - vb->first_line;

    // this loop can consume up to 25-50% of the entire decompress compute time (tested with 1KGP data)
    // note: we do memory assignment 64 bit at time (its about 10% faster than byte-by-byte)
    
    for (uint32_t ht_i=0; ht_i < ht_i_after; ht_i += 8) {

#ifdef __LITTLE_ENDIAN__
        *(next++) = ((uint64_t)(uint8_t)ht_columns_data[ht_i    ][vb_line_i]      ) |  // this is LITTLE ENDIAN order
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 1][vb_line_i] << 8 ) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 2][vb_line_i] << 16) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 3][vb_line_i] << 24) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 4][vb_line_i] << 32) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 5][vb_line_i] << 40) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 6][vb_line_i] << 48) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 7][vb_line_i] << 56) ;  // no worries if num_haplotypes_per_line is not a multiple of 4 - we have extra columns of zero
#else
        *(next++) = ((uint64_t)(uint8_t)ht_columns_data[ht_i    ][vb_line_i] << 56) |  // this is BIG ENDIAN order
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 1][vb_line_i] << 48) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 2][vb_line_i] << 40) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 3][vb_line_i] << 32) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 4][vb_line_i] << 24) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 5][vb_line_i] << 16) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 6][vb_line_i] << 8 ) |
                    ((uint64_t)(uint8_t)ht_columns_data[ht_i + 7][vb_line_i]      ) ;  
#endif
    }

    if (flag_show_alleles)
        printf ("Line %-2u : %.*s\n", vb_line_i, (int)vb->ht_one_array.len, vb->ht_one_array.data);

    COPY_TIMER (codec_ht_piz_get_one_line);
}

// Explanation of the reconstruction process of VCF haplotype matrix data compressed with the HT codec:
// 1) The GT_HT and GT_HT_INDEX sections as decompressed normally using the sub-codecs
// 2) When reconstructing, this function codec_ht_reconstruct is called for every haplotype value needed (eg the "1" in "1/0"),
//    by PIZ as the LT_HT reconstructor. It takes data from GT_HT, consulting GT_HT_INDEX
//    to find the correct column in the matrix, and skips '*'s (missing haplotypes due to mixed ploidy, missing samples,
//    or lines missing GT in FORMAT) until it finds a valid haplotype value.
void codec_ht_reconstruct (VBlock *vb, Context *ctx)
{
    if (vb->dont_show_curr_line) return;

    // get one row of the haplotype matrix for this line into vb->ht_one_array if we don't have it already
    if (vb->ht_one_array_line_i != vb->line_i) {
        codec_ht_piz_get_one_line (vb);
        vb->ht_one_array.len = 0; // length of data consumed
        vb->ht_one_array_line_i = vb->line_i;
    }

    // find next allele - skipping unused spots ('*')
    uint8_t ht = '*';
    while (ht == '*' && vb->ht_one_array.len < vb->num_haplotypes_per_line)
        ht = *ENT(uint8_t, vb->ht_one_array, vb->ht_one_array.len++);

    if (ht == '.' || IS_DIGIT(ht)) 
        RECONSTRUCT1 (ht);
    
    else if (ht == '*') 
        ABORT ("Error in codec_ht_reconstruct: reconstructing txt_line=%u vb_i=%u: unexpected end of ctx->local data in %s (len=%u)", 
               vb->line_i, vb->vblock_i, ctx->name, (uint32_t)ctx->local.len)
    
    else { // allele 10 to 99 (ascii 58 to 147)
        RECONSTRUCT_INT (ht - '0');
    }
}
