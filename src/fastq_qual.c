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

void fastq_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qual))
{
    START_TIMER;

    decl_ctx (FASTQ_QUAL);

    if (dl->deep_qual)
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

    // note: if we're copy QUAL from Deep (fully or partially), compress only the uncopied suffix starting at sam_seq_len
    uint32_t offset = dl->deep_qual ? dl->sam_seq_len : 0;

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = dl->dont_compress_QUAL ? 0 : MIN_(maximum_size, dl->qual.len - offset);
    
    if (!line_data) return; // only lengths were requested

    *line_data = Btxt (dl->qual.index) + offset;
    
    if (is_rev) *is_rev = 0;
}

void fastq_update_qual_len (VBlockP vb, uint32_t line_i, uint32_t new_len) 
{ 
    ZipDataLineFASTQ *dl = DATA_LINE (line_i);

    uint32_t offset = dl->deep_qual ? dl->sam_seq_len : 0; // as in fastq_zip_qual

    // note: new_len is based line_data_len returned from fastq_zip_qual and hence does not include offset - we need to add it back.
    dl->qual.len = new_len + offset;
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
