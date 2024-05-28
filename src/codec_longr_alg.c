// ------------------------------------------------------------------
//   codec_longr_alg.c
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.
//   Copyright claimed on additions and modifications vs EnanoFASTQ.
//
// The code in this file is partially derived from: https://github.com/guilledufort/EnanoFASTQ/tree/d6119dbccd19485c8222ec6d72516114401d1a08
// and has been modified extensively. To the extent there are still parts of this file that are
// subject to the copyright asserted in the EnanoFASTQ source code, those parts only are licensed according to
// the EnanoFASTQ license ("MIT License") below. All other parts are subject to the Genozip license in LICENSE.txt.
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

#include "compressor.h"

// algorithm parameters
#define N_BASES      6
#define B_AHEAD_OF_Q 3 

#define ROUND_UP(x) ((x)==0 ? 0 : (1 << ((x) - 1)))
#define DIV_ROUND(q, SHIFT) (((q) + ROUND_UP(SHIFT)) >> SHIFT) // divides by 2^SHIFT, rounding to the nearest

#define AVG_SHIFT   4  
#define TOTAL_ERR_SHIFT 4 

#define B_BITS      (2*N_BASES)
#define DIFQ_BITS   4 
#define QUAL_BITS   5 
#define AVG_BITS    QUAL_BITS
#define ERR_BITS    2
#define UNUSED_BITS (32 - B_BITS - DIFQ_BITS - QUAL_BITS - AVG_BITS - ERR_BITS)
#define CHAN_BITS   (DIFQ_BITS + QUAL_BITS + AVG_BITS + ERR_BITS)
#define Q_BITS      (QUAL_BITS + DIFQ_BITS)
#define NCTX_BITS   (B_BITS + DIFQ_BITS + QUAL_BITS)

#if (CHAN_BITS > 16)
#error need to increase width of LongrState.base_chan
#endif

#define LONGR_NUM_CHANNELS (1 << CHAN_BITS)

typedef union { // 32 bit
    struct {
        uint32_t B      : B_BITS;      // 6 neighboring bases - ACGT encoding
        uint32_t difq   : DIFQ_BITS;   // quantized delta vs previous qual
        uint32_t qbin   : QUAL_BITS;   // quantized current qual
        uint32_t avg    : AVG_BITS;
        uint32_t err_c  : ERR_BITS;
        uint32_t unused : UNUSED_BITS;
    } bits;

    struct {
        uint32_t na     : B_BITS;      // not included in the channel
        uint32_t n      : CHAN_BITS;   // qbin, difq, avg, err_c
        uint32_t unused : UNUSED_BITS;
    } channel;

    struct {
        uint32_t na1    : B_BITS;      // not included in the Q
        uint32_t n      : Q_BITS;      // qbin and difq
        uint32_t na2    : AVG_BITS + ERR_BITS; 
        uint32_t unused : UNUSED_BITS;
    } Q;

    struct {
        uint32_t n      : NCTX_BITS; // B, qbin, difq
        uint32_t na     : AVG_BITS + ERR_BITS;  // not included in nctx
        uint32_t unused : UNUSED_BITS;
    } nctx;

    uint32_t value;
} LongrChannel;

typedef struct {
    uint8_t *value_to_bin;
    uint16_t chan_avgs_sums[1 << NCTX_BITS];
    uint16_t chan_avgs_err_sums[1 << NCTX_BITS];
    uint32_t chan_err_avgs_total[1 << Q_BITS];
    LongrChannel chan;
    uint32_t chan_num_bases[1 << CHAN_BITS];
    uint16_t *base_chan; // channel of each base quality
    uint32_t next_base;  // index into base_chan
} LongrState;

static inline void codec_longr_update_state (LongrState *s, uint8_t b, int32_t q1, int32_t q2) 
{
    LongrChannel chan = s->chan; // automatic var for effeciency
    
    // update averages
    int32_t err = q1 - (DIV_ROUND(s->chan_avgs_sums[chan.nctx.n], AVG_SHIFT)); 
    s->chan_avgs_sums[chan.nctx.n] += err;
    int32_t abs_err = ABS(err);
    s->chan_avgs_err_sums[chan.nctx.n] += abs_err - (DIV_ROUND(s->chan_avgs_err_sums[chan.nctx.n], AVG_SHIFT));
    s->chan_err_avgs_total[chan.Q.n]   += abs_err - (DIV_ROUND(s->chan_err_avgs_total[chan.Q.n], TOTAL_ERR_SHIFT));

    chan.bits.B = (chan.bits.B << 2) | b; // shift and add new base 

    chan.bits.difq = MIN_((1L<<DIFQ_BITS)-1, INTERLACE(int32_t, q1-q2)); // difference might be [-93,93], interlaced to [0,186] and capped 

    chan.bits.qbin = s->value_to_bin[q1];
    
    chan.bits.avg = s->value_to_bin[DIV_ROUND(s->chan_avgs_sums[chan.nctx.n], AVG_SHIFT)];

    uint32_t total_err_avg = (s->chan_err_avgs_total[chan.Q.n] >> (TOTAL_ERR_SHIFT - AVG_SHIFT));
    uint32_t avg_err = s->chan_avgs_err_sums[chan.nctx.n];    
    
    chan.bits.err_c = (avg_err < (total_err_avg >> 1)) ? 0
                    : (avg_err < total_err_avg)        ? 1
                    : (avg_err < (total_err_avg << 1)) ? 2
                    :                                    3;
    
    s->chan = chan; // copy back to state
}

static void codec_longr_alg_init (LongrState *state) 
{
    ASSERT (sizeof (LongrChannel) == 4, "Expecting sizeof (LongrChannel)=%d == 4", (int)sizeof (LongrChannel));

    memset (state->chan_err_avgs_total, 1 << TOTAL_ERR_SHIFT, (1 << Q_BITS) * sizeof(uint32_t));

    for (uint32_t n=0 ; n < (1 << NCTX_BITS); n++) 
        state->chan_avgs_sums[n] = ((LongrChannel){ .nctx.n = n }).bits.qbin << AVG_SHIFT;
}

static void codec_longr_alg_init_read (LongrState *state, STRp(seq), bool is_rev)
{
    ASSERTNOTNULL (seq);
    state->chan = (LongrChannel){};

    // get first channel
    for (int32_t i=0; i < B_AHEAD_OF_Q; i++) 
        codec_longr_update_state (state, 
                                  is_rev ? acgt_encode_comp[(seq_len-1-i) >= 0 ? seq[seq_len-1-i] : 'T']
                                         : acgt_encode[i < seq_len ? seq[i] : 'A'],
                                  0, 0);
}

static inline char codec_longr_next_base (STRp(seq), uint32_t base_i)
{
    return (base_i + B_AHEAD_OF_Q < seq_len) ? seq[base_i + B_AHEAD_OF_Q] : 'A'; 
}

static inline char codec_longr_next_base_rev (STRp(seq), uint32_t base_i)
{
    return (base_i >= B_AHEAD_OF_Q) ? seq[base_i - B_AHEAD_OF_Q] : 'T';
}
