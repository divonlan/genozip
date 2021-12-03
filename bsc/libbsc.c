/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Compression/decompression functions                       */
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

#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "../libdeflate/libdeflate.h"

#include "platform.h"
#include "libbsc.h"
#include "bwt.h"
#include "lzp.h"
#include "coder.h"

void* (* bsc_malloc_do)(void *vb, size_t size, const char *,uint32_t) = NULL;
void  (* bsc_free_do)(void *vb, void* address, const char *,uint32_t) = NULL;

int bsc_init_full(int features, 
                  void* (* malloc)(void *vb, size_t size, const char *,uint32_t), 
                  void (* free)(void *vb, void* address, const char *,uint32_t))
{
    bsc_malloc_do = malloc;
    bsc_free_do   = free;
    return bsc_coder_init(features);
}

void *bsc_zero_malloc (void *vb, size_t size)
{
    void *address = bsc_malloc (vb, size);
    if (address) memset (address, 0, size);
    return address;
}

int bsc_store(void *vb, const unsigned char * input, unsigned char * output, int n, int features)
{
    unsigned int adler32_data = libdeflate_adler32 (1, input, n);

    memmove(output + LIBBSC_HEADER_SIZE, input, n);
    *(int *)(output +  0) = n + LIBBSC_HEADER_SIZE;
    *(int *)(output +  4) = n;
    *(int *)(output +  8) = 0;
    *(int *)(output + 12) = 0;
    *(int *)(output + 16) = adler32_data;
    *(int *)(output + 20) = adler32_data;
    *(int *)(output + 24) = libdeflate_adler32(1, output, 24);
    return n + LIBBSC_HEADER_SIZE;
}

static int bsc_compress_inplace (void *vb, unsigned char * data, int n, int lzpHashSize, int lzpMinLen, int blockSorter, int coder, int features)
{
    int             indexes[256];
    unsigned char   num_indexes;

    int mode = 0;

    switch (blockSorter)
    {
        case LIBBSC_BLOCKSORTER_BWT : mode = LIBBSC_BLOCKSORTER_BWT; break;

        default : return LIBBSC_BAD_PARAMETER;
    }

    switch (coder)
    {
        case LIBBSC_CODER_QLFC_STATIC   : mode += (LIBBSC_CODER_QLFC_STATIC   << 5); break;
        case LIBBSC_CODER_QLFC_ADAPTIVE : mode += (LIBBSC_CODER_QLFC_ADAPTIVE << 5); break;

        default : return LIBBSC_BAD_PARAMETER;
    }

    if (lzpMinLen != 0 || lzpHashSize != 0)
    {
        if (lzpMinLen < 4 || lzpMinLen > 255) return LIBBSC_BAD_PARAMETER;
        if (lzpHashSize < 10 || lzpHashSize > 28) return LIBBSC_BAD_PARAMETER;
        mode += (lzpMinLen << 8);
        mode += (lzpHashSize << 16);
    }
    if (n < 0 || n > 1073741824) return LIBBSC_BAD_PARAMETER;
    if (n <= LIBBSC_HEADER_SIZE)
    {
        return bsc_store(vb, data, data, n, features);
    }

    unsigned int adler32_data = libdeflate_adler32 (1, data, n);

    int lzSize = n;
    if (mode != (mode & 0xff))
    {
        unsigned char * buffer = (unsigned char *)bsc_malloc (vb, n);
        if (buffer == NULL) return LIBBSC_NOT_ENOUGH_MEMORY;

        lzSize = bsc_lzp_compress(vb, data, buffer, n, lzpHashSize, lzpMinLen, features);
        if (lzSize < LIBBSC_NO_ERROR)
        {
            lzSize = n; mode &= 0xff;
        }
        else
        {
            memcpy(data, buffer, lzSize);
        }

        bsc_free (vb, buffer);
    }

    if (lzSize <= LIBBSC_HEADER_SIZE)
    {
        blockSorter = LIBBSC_BLOCKSORTER_BWT;
        mode = (mode & 0xffffffe0) | LIBBSC_BLOCKSORTER_BWT;
    }

    int index = LIBBSC_BAD_PARAMETER; num_indexes = 0;
    switch (blockSorter)
    {
        case LIBBSC_BLOCKSORTER_BWT : index = bsc_bwt_encode(vb, data, lzSize, &num_indexes, indexes, features); break;
        default : return LIBBSC_BAD_PARAMETER;
    }

    if (n < 64 * 1024) num_indexes = 0;

    if (index < LIBBSC_NO_ERROR)
    {
        return index;
    }

    unsigned char * buffer = (unsigned char *)bsc_malloc (vb, lzSize + 4096);
    if (buffer)
    {
        int result = bsc_coder_compress(vb, data, buffer, lzSize, coder, features);
        if (result >= LIBBSC_NO_ERROR) memcpy(data + LIBBSC_HEADER_SIZE, buffer, result);
        bsc_free (vb, buffer);
        if ((result < LIBBSC_NO_ERROR) || (result + 1 + 4 * num_indexes >= n))
        {
            return LIBBSC_NOT_COMPRESSIBLE;
        }
        {
            if (num_indexes > 0)
            {
                memcpy(data + LIBBSC_HEADER_SIZE + result, indexes, 4 * num_indexes);
            }
            data[LIBBSC_HEADER_SIZE + result + 4 * num_indexes] = num_indexes;
            result += 1 + 4 * num_indexes;
        }
        *(int *)(data +  0) = result + LIBBSC_HEADER_SIZE;
        *(int *)(data +  4) = n;
        *(int *)(data +  8) = mode;
        *(int *)(data + 12) = index;
        *(int *)(data + 16) = adler32_data;
        *(int *)(data + 20) = libdeflate_adler32 (1, data + LIBBSC_HEADER_SIZE, result);
        *(int *)(data + 24) = libdeflate_adler32 (1, data, 24);
        return result + LIBBSC_HEADER_SIZE;
    }

    return LIBBSC_NOT_ENOUGH_MEMORY;
}

int bsc_compress (void *vb, const unsigned char * input, unsigned char * output, int n, int lzpHashSize, int lzpMinLen, int blockSorter, int coder, int features)
{
    if (input == output)
        return bsc_compress_inplace (vb, output, n, lzpHashSize, lzpMinLen, blockSorter, coder, features);

    int             indexes[256];
    unsigned char   num_indexes;

    int mode = 0;

    switch (blockSorter)
    {
        case LIBBSC_BLOCKSORTER_BWT : mode = LIBBSC_BLOCKSORTER_BWT; break;
        default : return LIBBSC_BAD_PARAMETER;
    }

    switch (coder)
    {
        case LIBBSC_CODER_QLFC_STATIC   : mode += (LIBBSC_CODER_QLFC_STATIC   << 5); break;
        case LIBBSC_CODER_QLFC_ADAPTIVE : mode += (LIBBSC_CODER_QLFC_ADAPTIVE << 5); break;

        default : return LIBBSC_BAD_PARAMETER;
    }

    if (lzpMinLen != 0 || lzpHashSize != 0)
    {
        if (lzpMinLen < 4 || lzpMinLen > 255) return LIBBSC_BAD_PARAMETER;
        if (lzpHashSize < 10 || lzpHashSize > 28) return LIBBSC_BAD_PARAMETER;
        mode += (lzpMinLen << 8);
        mode += (lzpHashSize << 16);
    }
    if (n < 0 || n > 1073741824) return LIBBSC_BAD_PARAMETER;
    if (n <= LIBBSC_HEADER_SIZE)
    {
        return bsc_store(vb, input, output, n, features);
    }
    int lzSize = 0;
    if (mode != (mode & 0xff))
    {
        lzSize = bsc_lzp_compress(vb, input, output, n, lzpHashSize, lzpMinLen, features);
        if (lzSize < LIBBSC_NO_ERROR)
        {
            mode &= 0xff;
        }
    }
    if (mode == (mode & 0xff))
    {
        lzSize = n; memcpy(output, input, n);
    }

    if (lzSize <= LIBBSC_HEADER_SIZE)
    {
        blockSorter = LIBBSC_BLOCKSORTER_BWT;
        mode = (mode & 0xffffffe0) | LIBBSC_BLOCKSORTER_BWT;
    }

    int index = LIBBSC_BAD_PARAMETER; num_indexes = 0;
    switch (blockSorter)
    {
        case LIBBSC_BLOCKSORTER_BWT : index = bsc_bwt_encode(vb, output, lzSize, &num_indexes, indexes, features); break;
        default : return LIBBSC_BAD_PARAMETER;
    }

    if (n < 64 * 1024) num_indexes = 0;

    if (index < LIBBSC_NO_ERROR)
    {
        return index;
    }

    unsigned char * buffer = (unsigned char *)bsc_malloc (vb, lzSize + 4096);
    if (buffer)
    {
        int result = bsc_coder_compress(vb, output, buffer, lzSize, coder, features);
        if (result >= LIBBSC_NO_ERROR) memcpy(output + LIBBSC_HEADER_SIZE, buffer, result);
        bsc_free (vb, buffer);

        if ((result < LIBBSC_NO_ERROR) || (result + 1 + 4 * num_indexes >= n))
            return bsc_store(vb, input, output, n, features);

        if (num_indexes > 0)
            memcpy(output + LIBBSC_HEADER_SIZE + result, indexes, 4 * num_indexes);

        output[LIBBSC_HEADER_SIZE + result + 4 * num_indexes] = num_indexes;
        result += 1 + 4 * num_indexes;

        *(int *)(output +  0) = result + LIBBSC_HEADER_SIZE;
        *(int *)(output +  4) = n;
        *(int *)(output +  8) = mode;
        *(int *)(output + 12) = index;
        *(int *)(output + 16) = libdeflate_adler32 (1, input, n);
        *(int *)(output + 20) = libdeflate_adler32 (1, output + LIBBSC_HEADER_SIZE, result);
        *(int *)(output + 24) = libdeflate_adler32 (1, output, 24);
        return result + LIBBSC_HEADER_SIZE;
    }

    return LIBBSC_NOT_ENOUGH_MEMORY;
}

int bsc_block_info (void *vb, const unsigned char * blockHeader, int headerSize, int * pBlockSize, int * pDataSize, int features)
{
    if (headerSize < LIBBSC_HEADER_SIZE)
    {
        return LIBBSC_UNEXPECTED_EOB;
    }

    if (*(unsigned int *)(blockHeader + 24) != libdeflate_adler32 (1, blockHeader, 24))
    {
        return LIBBSC_DATA_CORRUPT;
    }

    int blockSize    = *(int *)(blockHeader +  0);
    int dataSize     = *(int *)(blockHeader +  4);
    int mode         = *(int *)(blockHeader +  8);
    int index        = *(int *)(blockHeader + 12);

    int lzpHashSize  = (mode >> 16) & 0xff;
    int lzpMinLen    = (mode >>  8) & 0xff;
    int coder        = (mode >>  5) & 0x7;
    int blockSorter  = (mode >>  0) & 0x1f;

    int test_mode = 0;

    switch (blockSorter)
    {
        case LIBBSC_BLOCKSORTER_BWT : test_mode = LIBBSC_BLOCKSORTER_BWT; break;
        default : if (blockSorter > 0) return LIBBSC_DATA_CORRUPT;
    }

    switch (coder)
    {
        case LIBBSC_CODER_QLFC_STATIC   : test_mode += (LIBBSC_CODER_QLFC_STATIC   << 5); break;
        case LIBBSC_CODER_QLFC_ADAPTIVE : test_mode += (LIBBSC_CODER_QLFC_ADAPTIVE << 5); break;

        default : if (coder > 0) return LIBBSC_DATA_CORRUPT;
    }

    if (lzpMinLen != 0 || lzpHashSize != 0)
    {
        if (lzpMinLen < 4 || lzpMinLen > 255) return LIBBSC_DATA_CORRUPT;
        if (lzpHashSize < 10 || lzpHashSize > 28) return LIBBSC_DATA_CORRUPT;
        test_mode += (lzpMinLen << 8);
        test_mode += (lzpHashSize << 16);
    }

    if (test_mode != mode)
    {
        return LIBBSC_DATA_CORRUPT;
    }

    if (blockSize < LIBBSC_HEADER_SIZE || blockSize > LIBBSC_HEADER_SIZE + dataSize)
    {
        return LIBBSC_DATA_CORRUPT;
    }

    if (index < 0 || index > dataSize)
    {
        return LIBBSC_DATA_CORRUPT;
    }

    if (pBlockSize != NULL) *pBlockSize = blockSize;
    if (pDataSize != NULL) *pDataSize = dataSize;

    return LIBBSC_NO_ERROR;
}

static int bsc_decompress_inplace (void *vb, unsigned char * data, int inputSize, int outputSize, int features)
{
    int             indexes[256];
    unsigned char   num_indexes;

    int blockSize = 0, dataSize = 0;

    int info = bsc_block_info (vb, data, inputSize, &blockSize, &dataSize, features);
    if (info != LIBBSC_NO_ERROR)
    {
        return info;
    }

    if (inputSize < blockSize || outputSize < dataSize)
    {
        return LIBBSC_UNEXPECTED_EOB;
    }

    if (*(unsigned int *)(data + 20) != libdeflate_adler32 (1, data + LIBBSC_HEADER_SIZE, blockSize - LIBBSC_HEADER_SIZE))
    {
        return LIBBSC_DATA_CORRUPT;
    }

    int mode = *(int *)(data + 8);
    if (mode == 0)
    {
        memmove(data, data + LIBBSC_HEADER_SIZE, dataSize);
        return LIBBSC_NO_ERROR;
    }

    int             index           = *(int *)(data + 12);
    unsigned int    adler32_data    = *(int *)(data + 16);

    num_indexes = data[blockSize - 1];
    if (num_indexes > 0)
    {
        memcpy(indexes, data + blockSize - 1 - 4 * num_indexes, 4 * num_indexes);
    }

    int lzpHashSize  = (mode >> 16) & 0xff;
    int lzpMinLen    = (mode >>  8) & 0xff;
    int coder        = (mode >>  5) & 0x7;
    int blockSorter  = (mode >>  0) & 0x1f;

    int lzSize = LIBBSC_NO_ERROR;
    {
        unsigned char * buffer = (unsigned char *)bsc_malloc (vb, blockSize);
        if (buffer == NULL) return LIBBSC_NOT_ENOUGH_MEMORY;

        memcpy(buffer, data, blockSize);

        lzSize = bsc_coder_decompress(vb, buffer + LIBBSC_HEADER_SIZE, data, coder, features);

        bsc_free (vb, buffer);
    }
    if (lzSize < LIBBSC_NO_ERROR)
    {
        return lzSize;
    }

    int result;
    switch (blockSorter)
    {
        case LIBBSC_BLOCKSORTER_BWT : result = bsc_bwt_decode(vb, data, lzSize, index, num_indexes, indexes, features); break;
        default : return LIBBSC_DATA_CORRUPT;
    }
    if (result < LIBBSC_NO_ERROR)
    {
        return result;
    }

    if (mode != (mode & 0xff))
    {
        unsigned char * buffer = (unsigned char *)bsc_malloc (vb, lzSize);
        if (buffer)
        {
            memcpy(buffer, data, lzSize);
            result = bsc_lzp_decompress(vb, buffer, data, lzSize, lzpHashSize, lzpMinLen, features);
            bsc_free (vb, buffer);
            if (result < LIBBSC_NO_ERROR)
            {
                return result;
            }
            return result == dataSize ? (adler32_data == libdeflate_adler32 (1, data, dataSize) ? LIBBSC_NO_ERROR : LIBBSC_DATA_CORRUPT) : LIBBSC_DATA_CORRUPT;
        }
        return LIBBSC_NOT_ENOUGH_MEMORY;
    }

    return lzSize == dataSize ? (adler32_data == libdeflate_adler32 (1, data, dataSize) ? LIBBSC_NO_ERROR : LIBBSC_DATA_CORRUPT) : LIBBSC_DATA_CORRUPT;
}

int bsc_decompress (void *vb, const unsigned char * input, int inputSize, unsigned char * output, int outputSize, int features)
{
    int             indexes[256];
    unsigned char   num_indexes;

    if (input == output)
    {
        return bsc_decompress_inplace(vb, output, inputSize, outputSize, features);
    }

    int blockSize = 0, dataSize = 0;

    int info = bsc_block_info (vb, input, inputSize, &blockSize, &dataSize, features);
    if (info != LIBBSC_NO_ERROR)
    {
        return info;
    }

    if (inputSize < blockSize || outputSize < dataSize)
    {
        return LIBBSC_UNEXPECTED_EOB;
    }

    if (*(unsigned int *)(input + 20) != libdeflate_adler32 (1, input + LIBBSC_HEADER_SIZE, blockSize - LIBBSC_HEADER_SIZE))
    {
        return LIBBSC_DATA_CORRUPT;
    }

    int mode = *(int *)(input + 8);
    if (mode == 0)
    {
        memcpy(output, input + LIBBSC_HEADER_SIZE, dataSize);
        return LIBBSC_NO_ERROR;
    }

    int             index           = *(int *)(input + 12);
    unsigned int    adler32_data    = *(int *)(input + 16);

    num_indexes = input[blockSize - 1];
    if (num_indexes > 0)
    {
        memcpy(indexes, input + blockSize - 1 - 4 * num_indexes, 4 * num_indexes);
    }

    int lzpHashSize  = (mode >> 16) & 0xff;
    int lzpMinLen    = (mode >>  8) & 0xff;
    int coder        = (mode >>  5) & 0x7;
    int blockSorter  = (mode >>  0) & 0x1f;

    int lzSize = bsc_coder_decompress(vb, input + LIBBSC_HEADER_SIZE, output, coder, features);
    if (lzSize < LIBBSC_NO_ERROR)
    {
        return lzSize;
    }

    int result;
    switch (blockSorter)
    {
        case LIBBSC_BLOCKSORTER_BWT : result = bsc_bwt_decode(vb, output, lzSize, index, num_indexes, indexes, features); break;

        default : return LIBBSC_DATA_CORRUPT;
    }
    if (result < LIBBSC_NO_ERROR)
    {
        return result;
    }

    if (mode != (mode & 0xff))
    {
        unsigned char * buffer = (unsigned char *)bsc_malloc (vb, lzSize);
        if (buffer)
        {
            memcpy(buffer, output, lzSize);
            result = bsc_lzp_decompress (vb, buffer, output, lzSize, lzpHashSize, lzpMinLen, features);
            bsc_free (vb, buffer); 
            if (result < LIBBSC_NO_ERROR)
            {
                return result;
            }
            return result == dataSize ? (adler32_data == libdeflate_adler32 (1, output, dataSize) ? LIBBSC_NO_ERROR : LIBBSC_DATA_CORRUPT) : LIBBSC_DATA_CORRUPT;
        }
        return LIBBSC_NOT_ENOUGH_MEMORY;
    }

    return lzSize == dataSize ? (adler32_data == libdeflate_adler32 (1, output, dataSize) ? LIBBSC_NO_ERROR : LIBBSC_DATA_CORRUPT) : LIBBSC_DATA_CORRUPT;
}

/*-------------------------------------------------*/
/* End                                  libbsc.cpp */
/*-------------------------------------------------*/
