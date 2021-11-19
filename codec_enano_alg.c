// ------------------------------------------------------------------
//   codec_enano_alg.c
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt
//   Copyright claimed on additions and modifications vs EnanoFASTQ.
//
// This is code is derived from: https://github.com/guilledufort/EnanoFASTQ
// which has been modified extensively. The unmodified terms if the license EnanoFASTQ are as follows:
// 
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
// MIT License
//
// Copyright (c) 2020 Guillermo Dufort y Álvarez
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#include "genozip.h"
#include "buffer.h"
#include "codec.h"

// algorithm parameters
#define ENANO_PARAM_L  6
#define ENANO_NUM_CTXS 8192 //7168
#define B_HEAD_OF_Q (1 + ENANO_PARAM_L/2) // 4

#if B_HEAD_OF_Q != ENANO_MIN_READ_LEN
#error B_HEAD_OF_Q != ENANO_MIN_READ_LEN
#endif

#define ROUND_UP(x)    (x == 0? 0 : 1 << ((x) - 1))
/* Quantize quality scores q' = ((q + ROUND_UP(QUANT_Q)) >> QUANT_Q) */

#define DIV_ROUND(q, SHIFT) ((q + ROUND_UP(SHIFT)) >> SHIFT)

/* Keep as a power of 2 */
#define QMAX 128
#define A_LOG 2

#define QUANT_Q 3

#define LOG_AVGS 4
#define LOG_AVGS_ERRS 2

#define DIF_CTX_LOG 3
#define DIF_CANT (1 << DIF_CTX_LOG)

#define QUANT_D_LOG 5
#define QUANT_D_CANT (1 << QUANT_D_LOG)
#define QUANT_D_MASK (QUANT_D_CANT - 1)
#define QUANT_D_MAX QUANT_D_MASK

#define Q_LOG (7 - QUANT_Q)
#define TOTAL_Q_LOG (DIF_CTX_LOG + Q_LOG)

#define Q_CTX (1 << TOTAL_Q_LOG)

#define BC_CTX (1 << (LOG_AVGS + LOG_AVGS_ERRS))

#define CTX_CNT (BC_CTX * Q_CTX)

#define Q_LOG_CANT (1 << Q_LOG)

#define AVG_SHIFT 4
#define TOTAL_ERR_SHIFT 4

#define ROUND_UP(x)    (x == 0? 0 : 1 << ((x) - 1))
/* Quantize quality scores q' = ((q + ROUND_UP(QUANT_Q)) >> QUANT_Q) */

#define DIV_ROUND(q, SHIFT) ((q + ROUND_UP(SHIFT)) >> SHIFT)

#define AVG_SHIFT 4

#define B_CTX_LEN ENANO_PARAM_L
#define B_CTX     (1 << (B_CTX_LEN * A_LOG))
#define AVG_CANT  (B_CTX * Q_CTX)
#define B_MASK    (B_CTX - 1)

static const unsigned char QDif[8][14] = {
    {0, 1, 2, 3, 4, 5, 5, 6, 6, 7, 7},
    {1, 0, 2, 3, 4, 5, 5, 6, 6, 6, 7},
    {2, 1, 0, 3, 4, 5, 5, 6, 6, 6, 7},
    {3, 2, 1, 0, 3, 4, 5, 5, 6, 6, 7},
    {3, 3, 2, 1, 0, 4, 5, 5, 6, 6, 7},
    {3, 3, 3, 2, 1, 0, 4, 5, 6, 6, 7},
    {3, 3, 3, 2, 2, 1, 0, 4, 5, 6, 7},
    {0, 1, 2, 3, 4, 3, 4, 3, 4, 5, 6}
};

static const unsigned char QBin[128] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14,
    14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
};

typedef struct {
    uint16_t ctx_avgs_sums[AVG_CANT];
    uint16_t ctx_avgs_err_sums[AVG_CANT];
    uint32_t ctx_err_avgs_total[Q_CTX];
    uint32_t ctx_i, B_prev_ctx_i, Q_prev_ctx_i;
    uint8_t prev_q;
    ContextP ctxs[ENANO_NUM_CTXS];
} EnanoState;

static inline uint32_t codec_enano_get_ctx_i (EnanoState *state, uint8_t b /*0 to 3*/, uint8_t q1, uint8_t q2) 
{
    // Update averages
    uint32_t nctx = (state->B_prev_ctx_i << TOTAL_Q_LOG) + state->Q_prev_ctx_i;
    int32_t err = q1 - (DIV_ROUND(state->ctx_avgs_sums[nctx], AVG_SHIFT));
    state->ctx_avgs_sums[nctx] += err;
    uint32_t abs_err = ABS(err);
    state->ctx_avgs_err_sums[nctx] += abs_err - (DIV_ROUND(state->ctx_avgs_err_sums[nctx], AVG_SHIFT));
    state->ctx_err_avgs_total[state->Q_prev_ctx_i] += abs_err - (DIV_ROUND(state->ctx_err_avgs_total[state->Q_prev_ctx_i], TOTAL_ERR_SHIFT));

    state->B_prev_ctx_i = ((state->B_prev_ctx_i << A_LOG) + b) & B_MASK;

    uint32_t dif_ctx = 0;
    if (q1 < 7) 
        dif_ctx = QDif[q1][MIN_(q2, 10)];
    
    else {
        int dif = q2 - q1;
        if      (dif == 0) dif_ctx = 0;
        else if (dif < 0)  dif_ctx = QDif[7][MIN_(-2 * dif - 1, 9)];
        else               dif_ctx = QDif[7][MIN_(2 * dif, 10)];
    }

    state->Q_prev_ctx_i = ((dif_ctx << Q_LOG) + (uint32_t)QBin[q1]);

    nctx = (state->B_prev_ctx_i << TOTAL_Q_LOG) + state->Q_prev_ctx_i;

    uint32_t avg = QBin[DIV_ROUND(state->ctx_avgs_sums[nctx], AVG_SHIFT)];
    uint32_t total_err_avg = (state->ctx_err_avgs_total[state->Q_prev_ctx_i] >> (TOTAL_ERR_SHIFT - AVG_SHIFT));
    uint32_t avg_err = state->ctx_avgs_err_sums[nctx];
    uint32_t err_c = 0;

    if      (avg_err < (total_err_avg >> 1)) err_c = 0;
    else if (avg_err < total_err_avg)        err_c = 1;
    else if (avg_err < (total_err_avg << 1)) err_c = 2;
    else                                     err_c = 3;

    return (err_c << (TOTAL_Q_LOG + LOG_AVGS)) + (avg << (TOTAL_Q_LOG)) + state->Q_prev_ctx_i;
}

static void codec_enano_alg_init (EnanoState *state, STRp(seq)) // seq must be at least 4 bases
{
    for (uint32_t s_ctx=0; s_ctx < B_CTX; s_ctx++) 
        for (uint32_t dif=0; dif < DIF_CANT; dif++) 
            for (uint32_t q_quant = 0; q_quant < Q_LOG_CANT; q_quant++) {
                uint32_t avg_ctx = (s_ctx << TOTAL_Q_LOG) + (dif << Q_LOG) + q_quant;
                state->ctx_avgs_sums[avg_ctx] = q_quant << AVG_SHIFT;
            }

    memset (state, 0, (char*)state->ctxs - (char*)state); // clear all except for contexts
    memset (state->ctx_err_avgs_total, 1 << TOTAL_ERR_SHIFT, Q_CTX * sizeof(uint16_t));

    // Get first context
    for (int32_t i=0; i < B_HEAD_OF_Q; i++) 
        state->ctx_i = codec_enano_get_ctx_i (state, acgt_encode[i < seq_len ? seq[i] : 'A'], 0, 0);
}


static void codec_enano_encode_qual (VBlockP vb, EnanoState *state, STRp(seq), const char *qual) 
{
    codec_enano_alg_init (state, seq, seq_len);

    for (uint32_t i=0; i < seq_len; i++) {
        uint8_t q = (uint8_t)qual[i] - '!'; 

        buf_alloc (vb, &state->ctxs[state->ctx_i]->local, 1, 4096, uint8_t, CTX_GROWTH, "contexts->local");
        NEXTENT(uint8_t, state->ctxs[state->ctx_i]->local) = q;   

        // next b - from the position 4 ahead of q (so bases for the last for quality scores, bases are taken from the next line)   
        char next_b = (i + B_HEAD_OF_Q < seq_len) ? seq[i + B_HEAD_OF_Q] : 'A';

        state->ctx_i  = codec_enano_get_ctx_i (state, acgt_encode[(uint8_t)next_b], q, state->prev_q);
        state->prev_q = q;
    }
}

// this function is called if SAM/BAM read has been reverse complemented by the aligner, we analyze it in its original state - 
// revcomp of SEQ and reverse of QUAL.
static void codec_enano_encode_qual_rev (VBlockP vb, EnanoState *state, STRp(seq), const char *qual) 
{
    codec_enano_alg_init (state, seq, seq_len);

    for (int32_t i=seq_len-1; i >= 0; i--) {
        uint8_t q = (uint8_t)qual[i] - '!'; 

        buf_alloc (vb, &state->ctxs[state->ctx_i]->local, 1, 4096, uint8_t, CTX_GROWTH, "contexts->local");
        NEXTENT(uint8_t, state->ctxs[state->ctx_i]->local) = q;   

        // next b - from the position 4 ahead of q (so bases for the last for quality scores, bases are taken from the next line)   
        char next_b = (i - B_HEAD_OF_Q >= 0) ? seq[i - B_HEAD_OF_Q] : 'A';

        state->ctx_i  = codec_enano_get_ctx_i (state, 3 - (acgt_encode[(uint8_t)next_b]/*revcomp*/), q, state->prev_q);
        state->prev_q = q;
    }
}

void codec_enano_decode_qual (EnanoState *state, STRp(seq), char *qual) 
{
    codec_enano_alg_init (state, seq, seq_len);

    for (uint32_t i=0; i < seq_len; i++) {
        
        uint8_t q = NEXTLOCAL (uint8_t, state->ctxs[state->ctx_i]);

        qual[i] = q + '!';

        // update state     
        state->ctx_i  = codec_enano_get_ctx_i (state, ((i+B_HEAD_OF_Q < seq_len) ? seq[i+B_HEAD_OF_Q] : 'A'), q, state->prev_q);
        state->prev_q = q;
    }
}
