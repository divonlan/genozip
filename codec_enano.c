// ------------------------------------------------------------------
//   codec_enano.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

// This codec is for quality scores of Nanopore data, based on https://pubmed.ncbi.nlm.nih.gov/32470109/. This file contains Genozip code, and the
// algorithm itself, derived from ENANO source code, is located in codec_enano_alg.c

#include "genozip.h"
#include "codec.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "dict_id.h"
#include "reconstruct.h"
#include "strings.h"
#include "compressor.h"
#include "profiler.h"
#include "context.h"
#include "endianness.h"
#include "piz.h"
#include "strings.h"
#include "stats.h"

#include "codec_enano_alg.c"

#define enano_state con_cache

// must be called before merge, so that new contexts are added to z_file->contexts
void codec_enano_seg_init (VBlock *vb, ContextP qual_ctx)
{
    qual_ctx->lcodec = CODEC_ENANO;

    buf_alloc (vb, &qual_ctx->enano_state, 1, 0, EnanoState, 0, "enano_state");
    EnanoState *state = FIRSTENT (EnanoState, qual_ctx->enano_state);

    memset (state, 0, sizeof(EnanoState));

    for (int i=0; i < ENANO_NUM_CTXS; i++) {
        // dict_id carefully chosen to be unique in respect to ALT_KEY
        DictId dict_id = dict_id_make ((char[]){ 'a' + (i >> 10), (i >> 5) & 31, i & 31, 'E','N','A','N','O'}, 8, DTYPE_FIELD);
        state->ctxs[i] = ctx_get_ctx (vb, dict_id);
        state->ctxs[i]->local_dep = DEP_L1;    // compress the sub contexts in the second round
        state->ctxs[i]->lcodec = flag.fast ? CODEC_RANS8 : CODEC_ARITH8; // ARITH slightly better and slower than RANS, and no other alg comes close, so we can skip testing 8K contexts...
        stats_set_consolidation (vb, qual_ctx->did_i, 1, state->ctxs[i]->did_i);
    }
} 

bool codec_enano_compress (VBlock *vb, 
                           SectionHeader *header,    
                           const char *uncompressed,     // option 1 - not supported
                           uint32_t *uncompressed_len, 
                           LocalGetLineCB qual_callback, // option 2 - get one line
                           char *compressed, uint32_t *compressed_len, // in/out 
                           bool soft_fail)               // soft fail not supported
{
    //xxx 1. doesn't work - compare context selection 1 line of FQ with enano (print context_i in both)
    //xxx 2. need to support BAM (binary SEQ?)
    //xxx 3. decompressor
    START_TIMER;

    ContextP qual_ctx = CTX(VB_DT(DT_FASTQ) ? FASTQ_QUAL : SAM_QUAL);
    LocalGetLineCB *seq_callback = (VB_DT(DT_FASTQ) ? fastq_zip_seq : sam_zip_seq);
    
    EnanoState *state = FIRSTENT (EnanoState, qual_ctx->enano_state);

    for (uint64_t line_i=0; line_i < vb->lines.len; line_i++) {

        char *seq; uint32_t seq_len;
        bool is_rev;
        seq_callback (vb, line_i, pSTRa(seq), 0, &is_rev);

        char *qual; uint32_t qual_len; 
        qual_callback (vb, line_i, pSTRa(qual), vb->txt_data.len, NULL);
        
        ASSERT (seq_len == qual_len, "In line=%"PRIu64", expected seq_len=%u == qual_len=%u", line_i, seq_len, qual_len);
        ASSERT (seq_len >= ENANO_MIN_READ_LEN, "enano codec requires that all reads have at least %u bases, but line_i=%"PRIu64" read has only %u", ENANO_MIN_READ_LEN, line_i+1, seq_len);

        if (is_rev)
            codec_enano_encode_qual_rev (vb, state, STRa(seq), qual);
        else
            codec_enano_encode_qual (vb, state, STRa(seq), qual);
    }

    *compressed_len = 0;

    COPY_TIMER (compressor_enano);

    return true;
}

//--------------
// PIZ side
//--------------

// two options: 1. the length maybe given (textually) in snip/snip_len. in that case, it is used and vb->seq_len is updated.
// if snip_len==0, then the length is taken from seq_len.
void codec_enano_uncompress (VBlock *vb, Codec codec, uint8_t param,
                             STRp(compressed), Buffer *uncompressed_buf, uint64_t uncompressed_len,
                             Codec sub_codec)
{

}