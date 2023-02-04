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

#define QUAL_RECON_FROM_SAM 'S'
#define QUAL_RECON_FROM_CTX 'X'

void fastq_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, uint32_t qual_len, bool deep)
{
    if (!flag.deep) {
        CTX(FASTQ_QUAL)->local.len32 += qual_len;
        CTX(FASTQ_QUAL)->txt_len     += qual_len;
    }

    // note: in case of --deep, we use a special reconstructor. otherwise, we just reconstruct directly from local
    else if (flag.deep && deep) {
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_QUAL, QUAL_RECON_FROM_SAM }, 3, FASTQ_QUAL, qual_len);
        dl->dont_compress_QUAL = true;
    }

    else if (flag.deep && !deep) {
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_QUAL, QUAL_RECON_FROM_CTX }, 3, FASTQ_QUAL, qual_len);
        CTX(FASTQ_QUAL)->local.len32 += qual_len;
    }
}

// called to reconstruct QUAL in --deep 
SPECIAL_RECONSTRUCTOR (fastq_special_QUAL)
{
    if (snip_len == 1 && *snip == QUAL_RECON_FROM_SAM) {
        // xxx reconstruct from SAM
    }

    else if (ctx->ltype == LT_CODEC)
        codec_args[ctx->lcodec].reconstruct (VB, ctx->lcodec, ctx, NULL, 0); 

    else if (ctx->ltype == LT_SEQUENCE) 
        reconstruct_from_local_sequence (VB, ctx, NULL, 0, reconstruct); 

    else 
        ASSPIZ (false, "Invalid ltype=%s for %s", lt_name (ctx->ltype), ctx->tag_name);

    return NO_NEW_VALUE;
}

// callback function for compress to get data of one line (called by codec_bz2_compress)
COMPRESSOR_CALLBACK (fastq_zip_qual) 
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = dl->dont_compress_QUAL ? 0 : MIN_(dl->seq.len, maximum_size);
    
    if (!line_data) return; // only lengths were requested

    *line_data = Btxt (dl->qual_index);

    if (is_rev) *is_rev = 0;
}
