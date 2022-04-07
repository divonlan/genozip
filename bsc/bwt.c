/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Burrows Wheeler Transform                                 */
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
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <stdlib.h>
#include <memory.h>

#include "bwt.h"
#include "platform.h"
#include "libbsc.h"
#include "divsufsort.h"

int bsc_bwt_encode(void *vb, unsigned char * T, int n, unsigned char * num_indexes, int * indexes, int features)
{
    int index = divbwt (vb, T, T, NULL, n, num_indexes, indexes, 0);
    switch (index)
    {
        case -1 : return LIBBSC_BAD_PARAMETER;
        case -2 : return LIBBSC_NOT_ENOUGH_MEMORY;
    }
    return index;
}

static int bsc_unbwt_mergedTL_serial(void *vb, unsigned char * restrict T, unsigned int * restrict P, int n, int index)
{
    unsigned int bucket[ALPHABET_SIZE] = { 0 };

    for (int i = 0; i < index; ++i)
    {
        unsigned char c = T[i];
        P[i] = ((bucket[c]++) << 8) | c;
    }
    for (int i = index; i < n; ++i)
    {
        unsigned char c = T[i];
        P[i + 1] = ((bucket[c]++) << 8) | c;
    }

    for (int sum = 1, i = 0; i < ALPHABET_SIZE; ++i)
    {
        int tmp = sum; sum += bucket[i]; bucket[i] = tmp;
    }

    for (int p = 0, i = n - 1; i >= 0; --i)
    {
        unsigned int  u = P[p];
        unsigned char c = u & 0xff;
        T[i] = c; p = (u >> 8) + bucket[c];
    }

    return LIBBSC_NO_ERROR;
}

#define BWT_NUM_FASTBITS (17)

static int bsc_unbwt_biPSI_serial(void *vb, unsigned char * restrict T, unsigned int * restrict P, int n, int index)
{
    int *bucket;
    if (!(bucket = (int *)bsc_zero_malloc(vb, ALPHABET_SIZE * ALPHABET_SIZE * sizeof(int))))
        return LIBBSC_NOT_ENOUGH_MEMORY;

    unsigned short *fastbits;
    if (!(fastbits = (unsigned short *)bsc_malloc(vb,(1 + (1 << BWT_NUM_FASTBITS)) * sizeof(unsigned short)))) {
        return LIBBSC_NOT_ENOUGH_MEMORY;
        bsc_free(vb, bucket);
    }

    int count[ALPHABET_SIZE] = { 0 };

    int shift = 0; while ((n >> shift) > (1 << BWT_NUM_FASTBITS)) shift++;

    for (int i = 0; i < n; ++i) count[T[i]]++;

    for (int sum = 1, c = 0; c < ALPHABET_SIZE; ++c)
    {
        int tmp = sum; sum += count[c]; count[c] = tmp;
        if (count[c] != sum)
        {
            int * restrict bucket_p = &bucket[c << 8];

            int hi = index; if (sum < hi) hi = sum;
            for (int i = count[c]; i < hi; ++i) bucket_p[T[i]]++;

            int lo = index + 1; if (count[c] > lo) lo = count[c];
            for (int i = lo; i < sum; ++i) bucket_p[T[i - 1]]++;
        }
    }

    int lastc = T[0];
    for (int v = 0, sum = 1, c = 0; c < ALPHABET_SIZE; ++c)
    {
        if (c == lastc) sum++;

        int * restrict bucket_p = &bucket[c];
        for (int d = 0; d < ALPHABET_SIZE; ++d)
        {
            int tmp = sum; sum += bucket_p[d << 8]; bucket_p[d << 8] = tmp;
            if (bucket_p[d << 8] != sum)
            {
                for (; v <= ((sum - 1) >> shift); ++v) fastbits[v] = (c << 8) | d;
            }
        }
    }

    for (int i = 0; i < index; ++i)
    {
        unsigned char c = T[i];
        int           p = count[c]++;

        if (p < index) P[bucket[(c << 8) | T[p    ]]++] = i; else
        if (p > index) P[bucket[(c << 8) | T[p - 1]]++] = i;
    }

    for (int i = index; i < n; ++i)
    {
        unsigned char c = T[i];
        int           p = count[c]++;

        if (p < index) P[bucket[(c << 8) | T[p    ]]++] = i + 1; else
        if (p > index) P[bucket[(c << 8) | T[p - 1]]++] = i + 1;
    }

    for (int c = 0; c < ALPHABET_SIZE; ++c)
    {
        for (int d = 0; d < c; ++d)
        {
            int tmp = bucket[(d << 8) | c]; bucket[(d << 8) | c] = bucket[(c << 8) | d]; bucket[(c << 8) | d] = tmp;
        }
    }

    for (int p = index, i = 1; i < n; i += 2)
    {
        int c = fastbits[p >> shift]; while (bucket[c] <= p) c++;
        T[i - 1] = (unsigned char)(c >> 8); T[i] = (unsigned char)(c & 0xff);
        p = P[p];
    }

    T[n - 1] = (unsigned char)lastc;

    bsc_free(vb, fastbits); bsc_free(vb, bucket);
    return LIBBSC_NO_ERROR;
}

static int bsc_unbwt_reconstruct_serial(void *vb, unsigned char * T, unsigned int * P, int n, int index)
{
    if (n < 3 * 1024 * 1024) return bsc_unbwt_mergedTL_serial(vb, T, P, n, index);
    return bsc_unbwt_biPSI_serial(vb, T, P, n, index);
}

int bsc_bwt_decode(void *vb, unsigned char * T, int n, int index, unsigned char num_indexes, int * indexes, int features)
{
    ASSERT0 (T && n>=0 && index>0 && index <= n, "bad parameter");

    if (n <= 1)
    {
        return LIBBSC_NO_ERROR;
    }
    unsigned int * P = (unsigned int *)bsc_malloc (vb, (n + 1) * sizeof(unsigned int));
    if (P)
    {
        int result = bsc_unbwt_reconstruct_serial(vb, T, P, n, index);
        bsc_free (vb, P);
        return result;
    };
    return LIBBSC_NOT_ENOUGH_MEMORY;
}
