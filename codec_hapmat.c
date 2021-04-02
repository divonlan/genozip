// ------------------------------------------------------------------
//   codec_hapmat.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
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
#include "reconstruct.h"
#include "vcf_private.h"

//--------------
// ZIP side
//--------------

typedef struct {
    int num_alt_alleles;
    uint32_t index_in_original_line;
    uint32_t index_in_sorted_line;
} HaploTypeSortHelperIndex;

void codec_hapmat_comp_init (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    vb->ht_matrix_ctx->ltype  = LT_CODEC;
    vb->ht_matrix_ctx->lcodec = CODEC_HAPM;

    vb->hapmat_index_ctx         = ctx_get_ctx (vb, dict_id_FORMAT_GT_HT_INDEX);
    vb->hapmat_index_ctx->ltype  = LT_UINT32;

    // in --stats, consolidate stats into GT
    vb->ht_matrix_ctx->st_did_i = vb->hapmat_index_ctx->st_did_i = ctx_get_ctx (vb, dict_id_FORMAT_GT)->did_i;
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

static HaploTypeSortHelperIndex *codec_hapmat_count_alt_alleles (VBlockVCF *vb)
{
    START_TIMER; 

    buf_alloc_zero (vb, &vb->hapmat_helper_index_buf, 0, vb->ht_per_line, HaploTypeSortHelperIndex, 0, "hapmat_helper_index_buf");
    ARRAY (HaploTypeSortHelperIndex, helper_index, vb->hapmat_helper_index_buf);

    // build index array 
    for (uint32_t ht_i=0; ht_i < vb->ht_per_line; ht_i++) 
        helper_index[ht_i].index_in_original_line = ht_i;

    for (uint32_t line_i=0; line_i < (uint32_t)vb->lines.len; line_i++) {
        for (uint32_t ht_i=0; ht_i < vb->ht_per_line; ht_i++) {

            // we count as alt alleles : 1 - 99 (ascii 49 to 147)
            //             ref alleles : 0 . (unknown) * (missing)
            uint8_t one_ht = *ENT (uint8_t, vb->ht_matrix_ctx->local, line_i * vb->ht_per_line + ht_i);
            if (one_ht >= '1')
                helper_index[ht_i].num_alt_alleles++;
        }
    }
    COPY_TIMER (codec_hapmat_count_alt_alleles);

    return helper_index;
}

static void codec_hapmat_compress_one_array (VBlockP vb_, uint32_t ht_i, 
                                             char **line_data_1, uint32_t *line_data_len_1,
                                             uint32_t unused_maximum_len)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    uint32_t orig_col_i = ENT (HaploTypeSortHelperIndex, vb->hapmat_helper_index_buf, ht_i)->index_in_original_line; 
    const uint32_t increment = vb->ht_per_line; // automatic var that can be optimized to a register

    uint8_t *ht_data_ptr = ENT (uint8_t, vb->ht_matrix_ctx->local, orig_col_i);
    ARRAY (uint8_t, column, vb->hapmat_one_array);

    for (uint32_t line_i=0; line_i < vb->hapmat_one_array.len; line_i++, ht_data_ptr += increment) 
        column[line_i] = *ht_data_ptr;

    *line_data_1 = (char *)column;
    *line_data_len_1 = vb->hapmat_one_array.len;

    if (flag.show_alleles)
        iprintf ("Col %-2u : %.*s\n", orig_col_i, (int)vb->hapmat_one_array.len, column);
}

// build the reverse index that will allow access by the original index to the sorted array, to be included in the genozip file
static void codec_hapmap_compress_build_index (VBlockVCF *vb, HaploTypeSortHelperIndex *helper_index)
{
    for (uint32_t ht_i=0; ht_i < vb->ht_per_line ; ht_i++)
        helper_index[ht_i].index_in_sorted_line = ht_i;

    // sort array back to its original order
    qsort (helper_index, vb->ht_per_line, sizeof (HaploTypeSortHelperIndex), sort_by_original_index_comparator);

    // create a permutation index for the vblock
    // we populate the hapmat_hapmat_index_ctx local, and it will be written after us, as the context is create after the hapmat_ctx 
    buf_alloc (vb, &vb->hapmat_index_ctx->local, 0, vb->ht_per_line, uint32_t, 0, "contexts->local");
    vb->hapmat_index_ctx->local.len = vb->ht_per_line;
    
    ARRAY (uint32_t, hp_index, vb->hapmat_index_ctx->local);
    for (uint32_t ht_i=0; ht_i < vb->ht_per_line ; ht_i++)
        hp_index[ht_i] = BGEN32 (helper_index[ht_i].index_in_sorted_line);    
}

// sort haplogroups by alt allele count within the variant group, create an index for it, and split
// it to sample groups. for each sample a haplotype is just a string of 1 and 0 etc (could be other alleles too)
// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
bool codec_hapmat_compress (VBlock *vb_, 
                            SectionHeader *header,    
                            const char *uncompressed, // option 1 - compress contiguous data
                            uint32_t *uncompressed_len, 
                            LocalGetLineCB callback,  // option 2 - not supported
                            char *compressed, uint32_t *compressed_len /* in/out */, 
                            bool soft_fail)
{
    START_TIMER;
    
    ASSERT0 (uncompressed && !callback, "callback option not supported");

    VBlockVCF *vb = (VBlockVCF *)vb_;
    HaploTypeSortHelperIndex *helper_index = codec_hapmat_count_alt_alleles (vb);
    
    // sort the helper index by number of alt alleles
    qsort (helper_index, vb->ht_per_line, sizeof (HaploTypeSortHelperIndex), sort_by_alt_allele_comparator);

    vb->hapmat_one_array.len = vb->lines.len;
    buf_alloc_old (vb, &vb->hapmat_one_array, vb->hapmat_one_array.len, 1, "hapmat_one_array");
    
    if (flag.show_alleles) 
        iprint0 ("\nAfter transpose and sorting:\n");

    // compress the matrix one column at a time, by the order of helper index
    uint64_t save_lines_len = vb->lines.len;
    vb->lines.len = vb->ht_per_line; // temporarily set vb->lines.len to number of columns, as this is the number of time the callback will be called
    
    Codec sub_codec = codec_args[CODEC_HAPM].sub_codec;
    CodecCompress *compress = codec_args[sub_codec].compress;

    PAUSE_TIMER; //  don't include sub-codec compressor - it accounts for itself
    bool success = compress ((VBlockP)vb, header, 0, uncompressed_len, codec_hapmat_compress_one_array, compressed, compressed_len, soft_fail);    
    RESUME_TIMER (compressor_hapmat);

    // build data in hapmat_index_ctx 
    codec_hapmap_compress_build_index (vb, helper_index);

    vb->lines.len = save_lines_len;

    if (!success) return false; // soft fail (if it was hard fail, compress() already failed)

    // since codecs were already assigned to contexts before compression of all contexts begun, but
    // we just created this context now, we assign a codec manually
    // BUG: I can't this to work... reconstruction failed. no idea why. bug 215.
    //codec_assign_best_codec (vb, vb->hapmat_index_ctx, NULL, SEC_LOCAL, vb->ht_per_line * sizeof(uint32_t));
    vb->hapmat_index_ctx->lcodec = CODEC_BSC; // in the mean time, until bug 215 is fixed

    COPY_TIMER (compressor_hapmat);

    return true;
}

//--------------
// PIZ side
//--------------

// PIZ: for each haplotype column, retrieve its it address in the haplotype sections. Note that since the haplotype sections are
// transposed, each column will be a row, or a contiguous array, in the section data. This function returns an array
// of pointers, each pointer being a beginning of column data within the section array
void codec_hapmat_piz_calculate_columns (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    uint32_t ht_per_line = vb->ht_per_line = (uint32_t)(vb->ht_matrix_ctx->local.len / vb->lines.len);

    vb->hapmat_one_array.len = ht_per_line + 7; // +7 because depermuting_loop works on a word (32/64 bit) boundary
    buf_alloc_old (vb, &vb->hapmat_one_array, vb->hapmat_one_array.len, 1, "hapmat_one_array");
    buf_alloc_old (vb, &vb->hapmat_columns_data, sizeof (char *) * vb->hapmat_one_array.len, 1, "hapmat_columns_data"); // realloc for exact size (+15 is padding for 64b operations)

    // each entry is a pointer to the beginning of haplotype column located in vb->haplotype_sections_data
    // note: haplotype columns are permuted only within their own sample block
    ARRAY (const char *, hapmat_columns_data, vb->hapmat_columns_data); 
    ARRAY (const unsigned, permutatation_index, vb->hapmat_index_ctx->local);
    
    ASSERT0 (permutatation_index, "hapmat_hapmat_index_ctx.local is empty");

    // provide 7 extra zero-columns for the convenience of the permuting loop (supporting 64bit assignments)
    // note: txt_file->max_lines_per_vb will be zero if genozip file was created by redirecting output
    buf_alloc_old (vb, &vb->hapmat_column_of_zeros, MAX (txt_file->max_lines_per_vb, vb->lines.len), 1, "hapmat_column_of_zeros");
    buf_zero (&vb->hapmat_column_of_zeros);

    for (uint32_t ht_i = 0; ht_i < ht_per_line; ht_i++) 
        hapmat_columns_data[ht_i] = ENT (char, vb->ht_matrix_ctx->local, permutatation_index[ht_i] * vb->lines.len);

    for (unsigned ht_i=ht_per_line; ht_i < ht_per_line + 7; ht_i++)
        hapmat_columns_data[ht_i] = vb->hapmat_column_of_zeros.data;
}

// PIZ: build haplotype for a line - reversing the permutation and the transposal.
static inline void codec_hapmat_piz_get_one_line (VBlockVCF *vb)
{
    START_TIMER;

    if (flag.samples) buf_zero (&vb->hapmat_one_array); // if we're not filling in all samples, initialize to 0;

    ARRAY (const char *, hapmat_columns_data, vb->hapmat_columns_data); 
    uint32_t ht_i_after = vb->ht_per_line; // automatic variable - faster
    uint64_t *next = FIRSTENT (uint64_t, vb->hapmat_one_array);
    uint32_t vb_line_i = vb->line_i - vb->first_line;

    // this loop can consume up to 25-50% of the entire decompress compute time (tested with 1KGP data)
    // note: we do memory assignment 64 bit at time (its about 10% faster than byte-by-byte)
    
    for (uint32_t ht_i=0; ht_i < ht_i_after; ht_i += 8) {

#ifdef __LITTLE_ENDIAN__
        *(next++) = ((uint64_t)(uint8_t)hapmat_columns_data[ht_i    ][vb_line_i]      ) |  // this is LITTLE ENDIAN order
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 1][vb_line_i] << 8 ) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 2][vb_line_i] << 16) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 3][vb_line_i] << 24) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 4][vb_line_i] << 32) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 5][vb_line_i] << 40) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 6][vb_line_i] << 48) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 7][vb_line_i] << 56) ;  // no worries if ht_per_line is not a multiple of 4 - we have extra columns of zero
#else
        *(next++) = ((uint64_t)(uint8_t)hapmat_columns_data[ht_i    ][vb_line_i] << 56) |  // this is BIG ENDIAN order
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 1][vb_line_i] << 48) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 2][vb_line_i] << 40) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 3][vb_line_i] << 32) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 4][vb_line_i] << 24) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 5][vb_line_i] << 16) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 6][vb_line_i] << 8 ) |
                    ((uint64_t)(uint8_t)hapmat_columns_data[ht_i + 7][vb_line_i]      ) ;  
#endif
    }

    if (flag.show_alleles)
        iprintf ("Line %-2u : %.*s\n", vb_line_i, (int)vb->hapmat_one_array.len, vb->hapmat_one_array.data);

    COPY_TIMER (codec_hapmat_piz_get_one_line);
}

// Explanation of the reconstruction process of VCF haplotype matrix data compressed with the HT codec:
// 1) The GT_HT and GT_HT_INDEX sections as decompressed normally using the sub-codecs
// 2) When reconstructing, this function codec_hapmat_reconstruct is called for every haplotype value needed (eg the "1" in "1/0"),
//    by PIZ as the LT_CODEC reconstructor for CODEC_HAPM. It takes data from GT_HT, consulting GT_HT_INDEX
//    to find the correct column in the matrix, and skips '*'s (missing haplotypes due to mixed ploidy, missing samples,
//    or lines missing GT in FORMAT) until it finds a valid haplotype value.
void codec_hapmat_reconstruct (VBlock *vb_, Codec codec, Context *ctx)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    // get one row of the haplotype matrix for this line into vb->hapmat_one_array if we don't have it already
    if (vb->hapmat_one_array.param != vb->line_i) { // we store the line_i in param
        codec_hapmat_piz_get_one_line (vb);
        vb->hapmat_one_array.len = 0; // length of data consumed
        vb->hapmat_one_array.param = vb->line_i;
    }

    // find next allele - skipping unused spots ('*')
    uint8_t ht = '*';
    while (ht == '*' && vb->hapmat_one_array.len < vb->ht_per_line)
        ht = *ENT(uint8_t, vb->hapmat_one_array, vb->hapmat_one_array.len++);

    if (vb->dont_show_curr_line) return;

    switch (ht) {
        case '.': case '0' ... '9':
            RECONSTRUCT1 (ht); break;

        case 58 ... 147: { // allele 10 to 99 (ascii 58 to 147)
            RECONSTRUCT_INT (ht - '0'); break;
        }
        case '-': // ploidy padding (starting 10.0.2) - appears as the 2nd+ HT - unlike '*', these are counted in GT.repeats
            vb->txt_data.len--; // remove previous phase character;
            break;

        // "%|%" is used instead of "./." in case the previous sample had a | phase. This is meant to reduce entroy of GT.b250
        // in the case we have phased data, but missing samples appear as "./.".
        case '%': 
            if (*LASTENT (char, vb->txt_data) == '|' || *LASTENT (char, vb->txt_data) == '/') { // second % in "%|%" and act on the 2nd
                vb->txt_data.len -= 1;  // remove | or /
                RECONSTRUCT ("/.", 2);
            }
            else
                RECONSTRUCT1 ('.'); // first % in "%|%" 
            break;

        default: 
            ABORT ("Error in codec_hapmat_reconstruct: reconstructing txt_line=%u vb_i=%u: Invalid character found in decompressed HT array: '%c' (ASCII %u)", 
                   vb->line_i, vb->vblock_i, ht, ht);
    }
}

