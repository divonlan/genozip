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

void fastq_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qual), bool deep)
{
    START_TIMER;

    decl_ctx (FASTQ_QUAL);

    if (segconf.running)
        segconf_update_qual (STRa(qual)); // get stats on qual scores

    if (!flag.deep) {
        ctx->local.len32 += qual_len;
        ctx->txt_len     += qual_len;
    }

    // note: in case of --deep, we use a special reconstructor. otherwise, we just reconstruct directly from local
    else if (flag.deep && deep) {
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_deep_copy_QUAL }, 2, ctx, qual_len);
        dl->dont_compress_QUAL = true;
    }

    else if (flag.deep && !deep) {
        seg_simple_lookup (VB, ctx, qual_len);
        CTX(FASTQ_QUAL)->local.len32 += qual_len;
    }

    COPY_TIMER (fastq_seg_QUAL);
}

// callback function for compress to get data of one line
COMPRESSOR_CALLBACK (fastq_zip_qual) 
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = dl->dont_compress_QUAL ? 0 : MIN_(dl->seq.len, maximum_size);
    
    if (!line_data) return; // only lengths were requested

    *line_data = Btxt (dl->qual_index);

    if (is_rev) *is_rev = 0;
}
