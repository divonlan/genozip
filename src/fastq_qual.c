// ------------------------------------------------------------------
//   fastq_qual.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"
#include "deep.h"
#include "codec.h"

//---------------
// ZIP
//---------------

static void fastq_seg_QUAL_segconf (VBlockFASTQP vb, STRp(qual))
{
    for (uint32_t i=0; i < qual_len; i++)
        if (IS_QUAL_SCORE(qual[i]))
            segconf.qual_histo[QHT_QUAL][qual[i]-33].count++;
}

void fastq_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qual))
{
    START_TIMER;

    decl_ctx (FASTQ_QUAL);

    if (dl->deep_qual)
        fastq_deep_seg_QUAL (vb, dl, ctx, qual_len);

    else if (str_is_monochar (STRa(qual))) {
        seg_special1 (VB, FASTQ_SPECIAL_monochar_QUAL, qual[0], ctx, qual_len);
        dl->dont_compress_QUAL = true;
    }

    else {
        seg_simple_lookup (VB, ctx, qual_len);
        ctx->local.len32 += qual_len;

        if (segconf_running) segconf.nontrivial_qual = true; // not monochar
    }

    if (segconf_running)
        fastq_seg_QUAL_segconf (vb, STRa(qual));

    COPY_TIMER (fastq_seg_QUAL);
}

// callback function for compress to get data of one line
COMPRESSOR_CALLBACK (fastq_zip_qual) 
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb_line_i);
    char *qual = Btxt (dl->qual.index);

    // Deep: modify QUAL to contain only trim data (i.e. without the middle part copied from SAM)
    // careful to do this only once, as this function can be called multiple times (e.g. in CODEC_HOMP)
    if (dl->deep_qual && !dl->qual_is_trims_only && !dl->dont_compress_QUAL) {
        uint32_t trim_1_len = dl->sam_seq_offset;
        uint32_t trim_2_len = dl->qual.len - dl->sam_seq_offset - dl->sam_seq_len;
        memmove (qual + trim_1_len, qual + dl->qual.len - trim_2_len, trim_2_len);

        dl->qual.len = trim_1_len + trim_2_len;

        dl->qual_is_trims_only = true;
    }

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = dl->dont_compress_QUAL ? 0 : MIN_(maximum_size, dl->qual.len);
    
    if (__builtin_expect (!line_data, false)) return; // only lengths were requested

    *line_data = qual;
    
    if (is_rev) *is_rev = 0;
}

void fastq_update_qual_len (VBlockP vb, uint32_t line_i, uint32_t new_len) 
{ 
    DATA_LINE (line_i)->qual.len = new_len;
}

//---------------
// PIZ
//---------------

SPECIAL_RECONSTRUCTOR (fastq_special_monochar_QUAL)
{
    START_TIMER;

    char monochar = snip[0];
    char *c = BAFTtxt;
    char *aft = c + vb->seq_len;

    while (c < aft) *c++ = monochar;
    Ltxt += vb->seq_len;

    COPY_TIMER (fastq_special_monochar_QUAL);

    return NO_NEW_VALUE;
}
