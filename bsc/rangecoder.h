/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Range coder                                               */
/*-----------------------------------------------------------*/

/*--

This file is a part of bsc and/or libbsc, a program and a library for
lossless, block-sorting data compression.

   Copyright (c) 2009-2012 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Please see the file LICENSE for full copyright information and file AUTHORS
for full list of contributors.

See also the bsc and libbsc web site:
  http://libbsc.com/ for more information.

--*/


// ------------------------------------------------------------------
//   This file was extensively modified to adapt it to genozip. All modifications:
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef _LIBBSC_CODER_RANGECODER_H
#define _LIBBSC_CODER_RANGECODER_H

#include "platform.h"

typedef struct 
{
    union ari
    {
        struct u
        {
            unsigned int low32;
            unsigned int carry;
        } u;
        unsigned long long low;
    } ari;

    unsigned int ari_code;
    unsigned int ari_ffnum;
    unsigned int ari_cache;
    unsigned int ari_range;

    const unsigned short * restrict ari_input;
          unsigned short * restrict ari_output;
          unsigned short * restrict ari_outputEOB;
          unsigned short * restrict ari_outputStart;
} RangeCoder;

// private functions

static inline void OutputShort (RangeCoder *coder, unsigned short s)
{
    *coder->ari_output++ = s;
};

static inline unsigned short InputShort (RangeCoder *coder)
{
    return *coder->ari_input++;
};

static inline void ShiftLow (RangeCoder *coder)
{
    if (coder->ari.u.low32 < 0xffff0000U || coder->ari.u.carry)
    {
        OutputShort(coder, coder->ari_cache + coder->ari.u.carry);
        if (coder->ari_ffnum)
        {
            unsigned short s = coder->ari.u.carry - 1;
            do { OutputShort(coder, s); } while (--coder->ari_ffnum);
        }
        coder->ari_cache = coder->ari.u.low32 >> 16; coder->ari.u.carry = 0;
    } else coder->ari_ffnum++;
    coder->ari.u.low32 <<= 16;
}

// public functions

static inline bool CheckEOB (RangeCoder *coder)
{
    return coder->ari_output >= coder->ari_outputEOB;
}

static inline void InitEncoder (RangeCoder *coder, unsigned char * output, int outputSize)
{
    coder->ari_outputStart = (unsigned short *)output;
    coder->ari_output      = (unsigned short *)output;
    coder->ari_outputEOB   = (unsigned short *)(output + outputSize - 16);
    coder->ari.low         = 0;
    coder->ari_ffnum       = 0;
    coder->ari_cache       = 0;
    coder->ari_range       = 0xffffffff;
};

static inline int FinishEncoder (RangeCoder *coder)
{
    ShiftLow(coder); ShiftLow(coder); ShiftLow(coder);
    return (int)(coder->ari_output - coder->ari_outputStart) * sizeof(coder->ari_output[0]);
}

static inline void EncodeBit0 (RangeCoder *coder, int probability)
{
    coder->ari_range = (coder->ari_range >> 12) * probability;
    if (coder->ari_range < 0x10000)
    {
        coder->ari_range <<= 16; ShiftLow(coder);
    }
}

static inline void EncodeBit1 (RangeCoder *coder, int probability)
{
    unsigned int range = (coder->ari_range >> 12) * probability;
    coder->ari.low += range; coder->ari_range -= range;
    if (coder->ari_range < 0x10000)
    {
        coder->ari_range <<= 16; ShiftLow(coder);
    }
}

static inline void EncodeBit (RangeCoder *coder, unsigned int bit)
{
    if (bit) EncodeBit1(coder, 2048); else EncodeBit0(coder, 2048);
};

static inline void EncodeByte (RangeCoder *coder, unsigned int byte)
{
    for (int bit = 7; bit >= 0; --bit)
    {
        EncodeBit(coder, byte & (1 << bit));
    }
};

static inline void EncodeWord (RangeCoder *coder, unsigned int word)
{
    for (int bit = 31; bit >= 0; --bit)
    {
        EncodeBit(coder, word & (1 << bit));
    }
};

static inline void InitDecoder (RangeCoder *coder, const unsigned char * input)
{
    coder->ari_input = (unsigned short *)input;
    coder->ari_code  = 0;
    coder->ari_range = 0xffffffff;
    coder->ari_code  = (coder->ari_code << 16) | InputShort(coder);
    coder->ari_code  = (coder->ari_code << 16) | InputShort(coder);
    coder->ari_code  = (coder->ari_code << 16) | InputShort(coder);
};

static inline int DecodeBitProb (RangeCoder *coder, int probability)
{
    unsigned int range = (coder->ari_range >> 12) * probability;
    if (coder->ari_code >= range)
    {
        coder->ari_code -= range; coder->ari_range -= range;
        if (coder->ari_range < 0x10000)
        {
            coder->ari_range <<= 16; coder->ari_code = (coder->ari_code << 16) | InputShort(coder);
        }
        return 1;
    }
    coder->ari_range = range;
    if (coder->ari_range < 0x10000)
    {
        coder->ari_range <<= 16; coder->ari_code = (coder->ari_code << 16) | InputShort(coder);
    }
    return 0;
}

static inline unsigned int DecodeBit (RangeCoder *coder)
{
    return DecodeBitProb(coder, 2048);
}

static inline unsigned int DecodeByte (RangeCoder *coder)
{
    unsigned int byte = 0;
    for (int bit = 7; bit >= 0; --bit)
    {
        byte += byte + DecodeBit(coder);
    }
    return byte;
}

static inline unsigned int DecodeWord (RangeCoder *coder)
{
    unsigned int word = 0;
    for (int bit = 31; bit >= 0; --bit)
    {
        word += word + DecodeBit (coder);
    }
    return word;
}

#endif
