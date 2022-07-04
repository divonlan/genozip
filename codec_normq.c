// ------------------------------------------------------------------
//   codec_normq.c
//   Copyright (C) 2022-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

// fallback compression algorithm for SAM QUAL values LONGR and DOMQ are not applicable

#include "vblock.h"
#include "data_types.h"
#include "profiler.h"
#include "codec.h"
#include "sam_private.h" // this is a SAM-only codec

//--------------
// ZIP side
//--------------

COMPRESS(codec_normq_compress)
{
    START_TIMER;
    ASSERT0 (!uncompressed && get_line_cb, "only callback option is supported");

    ContextP qual_ctx = ECTX (((SectionHeaderCtx *)header)->dict_id); // may be SAM_QUAL or OPTION_U2_Z  
    BufferP qual_buf = &qual_ctx->local;

    // case: this is our second entry, after soft-failing. Just continue from where we stopped
    if (!soft_fail) goto do_compress;

    buf_alloc_exact (vb, *qual_buf, qual_buf->len,  char, "contexts->local"); 
    
    uint32_t next = 0;
    for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) {
        STRw(qual);
        bool is_rev;
        get_line_cb (vb, line_i, pSTRa (qual), CALLBACK_NO_SIZE_LIMIT, &is_rev);

        if (qual_len) { 
            if (is_rev) str_reverse (Bc(*qual_buf, next), qual, qual_len);
            else        memcpy      (Bc(*qual_buf, next), qual, qual_len);

            next += qual_len;
        }
    }

    qual_ctx->lcodec = CODEC_UNKNOWN;
    header->sub_codec = codec_assign_best_codec (vb, qual_ctx, &qual_ctx->local, SEC_LOCAL); // provide BufferP to override callback
    if (header->sub_codec == CODEC_UNKNOWN) header->sub_codec = CODEC_NONE; // really small

do_compress: ({});
    CodecCompress *compress = codec_args[header->sub_codec].compress;
    *uncompressed_len = qual_buf->len32;

    // make sure we have enough memory
    uint32_t min_required_compressed_len = codec_args[header->sub_codec].est_size (header->sub_codec, qual_buf->len);
    if (*compressed_len < min_required_compressed_len) {
        if (soft_fail) return false; // call me again with more memory
        ABORT ("%s: Compressing %s with %s need %u bytes, but allocated only %u", VB_NAME, qual_ctx->tag_name, codec_name(header->sub_codec), min_required_compressed_len, *compressed_len);
    }

    COPY_TIMER_COMPRESS (compressor_normq); // don't account for sub-codec compressor, it accounts for itself

    return compress (vb, header, B1STc(*qual_buf), uncompressed_len, NULL, STRa(compressed), false, name);
}

//--------------
// PIZ side
//--------------

CODEC_RECONSTRUCT (codec_normq_reconstruct)
{   
    if (!ctx->is_loaded) return;

    bool reconstruct = true;

    rom next_qual = Bc(ctx->local, ctx->next_local);

    if (*next_qual==' ') { // this is QUAL="*"
        sam_reconstruct_missing_quality (vb, reconstruct);
        ctx->next_local++;
    }

    else {
        if (reconstruct) {
            if (last_flags.rev_comp) str_reverse (BAFTtxt, next_qual, vb->seq_len);
            else                     memcpy      (BAFTtxt, next_qual, vb->seq_len);
        }

        ctx->next_local += vb->seq_len;
        vb->txt_data.len     += vb->seq_len;
    }
}
