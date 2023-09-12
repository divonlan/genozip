// ------------------------------------------------------------------
//   fastq_qual.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"
#include "seg.h"
#include "piz.h"
#include "deep.h"
#include "codec.h"
#include "segconf.h"

//---------------
// ZIP
//---------------

void fastq_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qual), bool deep)
{
    START_TIMER;

    decl_ctx (FASTQ_QUAL);

    if (deep)
        fastq_deep_seg_QUAL (vb, dl, ctx, qual_len);

    else if (str_is_monochar (STRa(qual))) {
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_monochar_QUAL, qual[0] }, 3, ctx, qual_len);
        dl->dont_compress_QUAL = true;
    }

    else {
        seg_simple_lookup (VB, ctx, qual_len);
        ctx->local.len32 += qual_len;

        if (segconf.running) segconf.nontrivial_qual = true; // not monochar
    }

    COPY_TIMER (fastq_seg_QUAL);
}

// callback function for compress to get data of one line
COMPRESSOR_CALLBACK (fastq_zip_qual) 
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    // note: in sam_seq_len (i.e. Deep and we copied the prefix of QUAL from SAM) - compress only the uncopied suffix
    *line_data_len  = dl->dont_compress_QUAL ? 0 : MIN_(maximum_size, dl->qual.len - dl->sam_seq_len);
    
    if (!line_data) return; // only lengths were requested

    *line_data = Btxt (dl->qual.index) + dl->sam_seq_len;
    
    if (is_rev) *is_rev = 0;
}

void fastq_update_qual_len (VBlockP vb, uint32_t line_i, uint32_t new_len) 
{ 
    ZipDataLineFASTQ *dl = DATA_LINE (line_i);

    // note: add back sam_seq_len bc it will be substracted again when sub-codec calls fastq_zip_qual
    dl->qual.len = new_len + dl->sam_seq_len;
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
