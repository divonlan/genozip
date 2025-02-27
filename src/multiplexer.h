// ------------------------------------------------------------------
//   multiplexer.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// note: separated from seg.h, to avoid including seg.h everywhere

#pragma once

#include "genozip.h"

#define BASE64_DICT_ID_LEN 14
#define MULTIPLEXER(n_channels)                         \
struct {                                                \
    /* all fields ara aligned - padding only at end  */ \
    ContextP ctx;                                       \
    bool no_stons;                                      \
    uint8_t special_code;                               \
    uint16_t num_channels;                              \
    uint32_t snip_len;                                  \
    DictId dict_ids[n_channels];                        \
    ContextP channel_ctx[n_channels];                   \
    char snip[BASE64_DICT_ID_LEN * (n_channels) + 6];   \
}                               
#define MUX ((MultiplexerP)mux)
#define MUX_CAPACITY(mux)    (sizeof((mux).dict_ids)/sizeof(DictId)) // max number of channels this mux can contain
#define MUX_CHANNEL_CTX(mux) ((ContextP *)((mux)->dict_ids + (mux)->num_channels))
#define MUX_SNIP(mux)        ((char *)(MUX_CHANNEL_CTX(mux) + (mux)->num_channels))
#define MUX_SIZEOF_SNIP(mux) (BASE64_DICT_ID_LEN * ((mux)->num_channels) + 6)

typedef MULTIPLEXER(1000) *MultiplexerP;
typedef const MULTIPLEXER(1000) *ConstMultiplexerP;

typedef MULTIPLEXER(2)  Multiplexer2,  *Multiplexer2P;
typedef MULTIPLEXER(3)  Multiplexer3,  *Multiplexer3P;
typedef MULTIPLEXER(4)  Multiplexer4,  *Multiplexer4P;
typedef MULTIPLEXER(5)  Multiplexer5,  *Multiplexer5P;
typedef MULTIPLEXER(6)  Multiplexer6,  *Multiplexer6P;
typedef MULTIPLEXER(7)  Multiplexer7,  *Multiplexer7P;
typedef MULTIPLEXER(10) Multiplexer10, *Multiplexer10P;
