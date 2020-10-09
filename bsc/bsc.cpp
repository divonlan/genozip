/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Block Sorting Compressor                                  */
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

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <memory.h>

#if defined(_WIN32)
  #include <windows.h>
  SIZE_T g_LargePageSize = 0;
#endif

#include "libbsc.h"
#include "filters.h"
#include "platform.h"

#pragma pack(push, 1)

#define LIBBSC_CONTEXTS_AUTODETECT   3

unsigned char bscFileSign[4] = {'b', 's', 'c', 0x31};

typedef struct BSC_BLOCK_HEADER
{
    long long       blockOffset;
    signed char     recordSize;
    signed char     sortingContexts;
} BSC_BLOCK_HEADER;

#pragma pack(pop)

int paramBlockSize                = 25 * 1024 * 1024;
int paramBlockSorter              = LIBBSC_BLOCKSORTER_BWT;
int paramCoder                    = LIBBSC_CODER_QLFC_STATIC;
int paramSortingContexts          = LIBBSC_CONTEXTS_FOLLOWING;

int paramEnableParallelProcessing = 1;
int paramEnableMultiThreading     = 1;
int paramEnableFastMode           = 0;
int paramEnableLargePages         = 0;
int paramEnableCUDA               = 0;
int paramEnableSegmentation       = 0;
int paramEnableReordering         = 0;
int paramEnableLZP                = 1;
int paramLZPHashSize              = 16;
int paramLZPMinLen                = 128;

int paramFeatures()
{
    int features =
        (paramEnableFastMode       ? LIBBSC_FEATURE_FASTMODE       : LIBBSC_FEATURE_NONE) |
        (paramEnableLargePages     ? LIBBSC_FEATURE_LARGEPAGES     : LIBBSC_FEATURE_NONE)
    ;

    return features;
}

#if defined(__GNUC__) && (defined(_GLIBCXX_USE_LFS) || defined(__MINGW32__))
    #define BSC_FSEEK fseeko64
    #define BSC_FTELL ftello64
    #define BSC_FILEOFFSET off64_t
#elif defined(_MSC_VER) && _MSC_VER >= 1400
    #define BSC_FSEEK _fseeki64
    #define BSC_FTELL _ftelli64
    #define BSC_FILEOFFSET __int64
#else
    #define BSC_FSEEK fseek
    #define BSC_FTELL ftell
    #define BSC_FILEOFFSET long
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__CYGWIN__) || defined(__MINGW32__) || defined(__BORLANDC__) || defined(_MSC_VER)
  #include <windows.h>
  double BSC_CLOCK() { return 0.001 * GetTickCount(); }
#elif defined (__unix) || defined (__linux__) || defined (__QNX__) || defined (_AIX)  || defined (__NetBSD__) || defined(macintosh) || defined (_MAC)
  #include <sys/time.h>
  double BSC_CLOCK() { timeval tv; gettimeofday(&tv, 0); return tv.tv_sec + tv.tv_usec * 0.000001; }
#else
  double BSC_CLOCK() { return (double)clock() / CLOCKS_PER_SEC; }
#endif

int segmentedBlock[256];

static void * bsc_default_malloc(void *vb, size_t size)
{
#if defined(_WIN32)
    if ((g_LargePageSize != 0) && (size >= 256 * 1024))
    {
        void * address = VirtualAlloc(0, (size + g_LargePageSize - 1) & (~(g_LargePageSize - 1)), MEM_COMMIT | MEM_LARGE_PAGES, PAGE_READWRITE);
        if (address != NULL) return address;
    }
    return VirtualAlloc(0, size, MEM_COMMIT, PAGE_READWRITE);
#else
    return malloc(size);
#endif
}

static void * bsc_default_zero_malloc(void *vb, size_t size)
{
#if defined(_WIN32)
    if ((g_LargePageSize != 0) && (size >= 256 * 1024))
    {
        void * address = VirtualAlloc(0, (size + g_LargePageSize - 1) & (~(g_LargePageSize - 1)), MEM_COMMIT | MEM_LARGE_PAGES, PAGE_READWRITE);
        if (address != NULL) return address;
    }
    return VirtualAlloc(0, size, MEM_COMMIT, PAGE_READWRITE);
#else
    return calloc(1, size);
#endif
}

static void bsc_default_free(void *vb, void * address)
{
#if defined(_WIN32)
    VirtualFree(address, 0, MEM_RELEASE);
#else
    free(address);
#endif
}

void Compression(char * argv[])
{
    if (!paramEnableLZP)
    {
        paramLZPHashSize = 0;
        paramLZPMinLen = 0;
    }

    FILE * fInput = fopen(argv[2], "rb");
    if (fInput == NULL)
    {
        fprintf(stderr, "Can't open input file: %s!\n", argv[2]);
        exit(1);
    }

    FILE * fOutput = fopen(argv[3], "wb");
    if (fOutput == NULL)
    {
        fprintf(stderr, "Can't create output file: %s!\n", argv[3]);
        exit(1);
    }

    if (BSC_FSEEK(fInput, 0, SEEK_END))
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
        exit(1);
    }

    BSC_FILEOFFSET fileSize = BSC_FTELL(fInput);
    if (fileSize < 0)
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
        exit(1);
    }

    if (BSC_FSEEK(fInput, 0, SEEK_SET))
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
        exit(1);
    }

    if (paramBlockSize > fileSize)
    {
        paramBlockSize = (int)fileSize;
    }

    if (fwrite(bscFileSign, sizeof(bscFileSign), 1, fOutput) != 1)
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[3]);
        exit(1);
    }

    int nBlocks = paramBlockSize > 0 ? (int)((fileSize + paramBlockSize - 1) / paramBlockSize) : 0;
    if (fwrite(&nBlocks, sizeof(nBlocks), 1, fOutput) != 1)
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[3]);
        exit(1);
    }

    double startTime = BSC_CLOCK();

    int segmentationStart = 0, segmentationEnd = 0;

    {
        unsigned char * buffer = (unsigned char *)bsc_malloc(NULL, paramBlockSize + LIBBSC_HEADER_SIZE);
        if (buffer == NULL)
        {
            {
                fprintf(stderr, "Not enough memory! Please check README file for more information.\n");
                exit(2);
            }
        }

        while (true)
        {
            BSC_FILEOFFSET  blockOffset     = 0;
            int             dataSize        = 0;
            {
                if ((feof(fInput) == 0) && (BSC_FTELL(fInput) != fileSize))
                {
                    {
                        double progress = (100.0 * (double)BSC_FTELL(fInput)) / fileSize;
                        fprintf(stdout, "\rCompressing %.55s(%02d%%)", argv[2], (int)progress);
                        fflush(stdout);
                    }

                    blockOffset = BSC_FTELL(fInput);

                    int currentBlockSize = paramBlockSize;
                    if (paramEnableSegmentation)
                    {
                        if (segmentationEnd - segmentationStart > 1) currentBlockSize = segmentedBlock[segmentationStart];
                    }

                    dataSize = (int)fread(buffer, 1, currentBlockSize, fInput);
                    if (dataSize <= 0)
                    {
                        fprintf(stderr, "\nIO error on file: %s!\n", argv[2]);
                        exit(1);
                    }

                    if (paramEnableSegmentation)
                    {
                        bool bSegmentation = false;

                        if (segmentationStart == segmentationEnd) bSegmentation = true;
                        if ((segmentationEnd - segmentationStart == 1) && (dataSize != segmentedBlock[segmentationStart])) bSegmentation = true;

                        if (bSegmentation)
                        {
                            segmentationStart = 0; segmentationEnd = bsc_detect_segments(NULL, buffer, dataSize, segmentedBlock, 256, paramFeatures());
                            if (segmentationEnd <= LIBBSC_NO_ERROR)
                            {
                                switch (segmentationEnd)
                                {
                                    case LIBBSC_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                                    default                         : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                                }
                                exit(2);
                            }
                        }

                        int newDataSize = segmentedBlock[segmentationStart++];
                        if (dataSize != newDataSize)
                        {
                            BSC_FILEOFFSET pos = BSC_FTELL(fInput) - dataSize + newDataSize;
                            BSC_FSEEK(fInput, pos, SEEK_SET);
                            dataSize = newDataSize;
                        }
                    }
                }
            }

            if (dataSize == 0) break;

            signed char recordSize = 1;
            if (paramEnableReordering)
            {
                recordSize = bsc_detect_recordsize(NULL, buffer, dataSize, paramFeatures());
                if (recordSize < LIBBSC_NO_ERROR)
                {
                    {
                        switch (recordSize)
                        {
                            case LIBBSC_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                            default                         : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                        }
                        exit(2);
                    }
                }
                if (recordSize > 1)
                {
                    int result = bsc_reorder_forward(NULL, buffer, dataSize, recordSize, paramFeatures());
                    if (result != LIBBSC_NO_ERROR)
                    {
                        {
                            switch (result)
                            {
                                case LIBBSC_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                                default                         : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                            }
                            exit(2);
                        }
                    }
                }
            }

            signed char sortingContexts = paramSortingContexts;
            if (paramSortingContexts == LIBBSC_CONTEXTS_AUTODETECT)
            {
                sortingContexts = bsc_detect_contextsorder(NULL, buffer, dataSize, paramFeatures());
                if (sortingContexts < LIBBSC_NO_ERROR)
                {
                    {
                        switch (sortingContexts)
                        {
                            case LIBBSC_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough memory!\n"); break;
                            default                         : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                        }
                        exit(2);
                    }
                }
            }
            if (sortingContexts == LIBBSC_CONTEXTS_PRECEDING)
            {
                int result = bsc_reverse_block(NULL, buffer, dataSize, paramFeatures());
                if (result != LIBBSC_NO_ERROR)
                {
                    {
                        fprintf(stderr, "\nInternal program error, please contact the author!\n");
                        exit(2);
                    }
                }
            }

            int blockSize = bsc_compress(NULL, buffer, buffer, dataSize, paramLZPHashSize, paramLZPMinLen, paramBlockSorter, paramCoder, paramFeatures());
            if (blockSize == LIBBSC_NOT_COMPRESSIBLE)
            {
                {
                    sortingContexts = LIBBSC_CONTEXTS_FOLLOWING; recordSize = 1;

                    BSC_FILEOFFSET pos = BSC_FTELL(fInput);
                    {
                        BSC_FSEEK(fInput, blockOffset, SEEK_SET);
                        if (dataSize != (int)fread(buffer, 1, dataSize, fInput))
                        {
                            fprintf(stderr, "\nInternal program error, please contact the author!\n");
                            exit(2);
                        }
                    }
                    BSC_FSEEK(fInput, pos, SEEK_SET);
                }

                blockSize = bsc_store(NULL, buffer, buffer, dataSize, paramFeatures());
            }
            if (blockSize < LIBBSC_NO_ERROR)
            {
                {
                    switch (blockSize)
                    {
                        case LIBBSC_NOT_ENOUGH_MEMORY       : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                        case LIBBSC_NOT_SUPPORTED           : fprintf(stderr, "\nSpecified compression method is not supported on this platform!\n"); break;
                        case LIBBSC_GPU_ERROR               : fprintf(stderr, "\nGeneral GPU failure! Please check README file for more information.\n"); break;
                        case LIBBSC_GPU_NOT_SUPPORTED       : fprintf(stderr, "\nYour GPU is not supported! Please check README file for more information.\n"); break;
                        case LIBBSC_GPU_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough GPU memory! Please check README file for more information.\n"); break;

                        default                             : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                    }
                    exit(2);
                }
            }
            {
                BSC_BLOCK_HEADER header = {blockOffset, recordSize, sortingContexts};

                if (fwrite(&header, sizeof(BSC_BLOCK_HEADER), 1, fOutput) != 1)
                {
                    fprintf(stderr, "\nIO error on file: %s!\n", argv[3]);
                    exit(1);
                }

                if ((int)fwrite(buffer, 1, blockSize, fOutput) != blockSize)
                {
                    fprintf(stderr, "\nIO error on file: %s!\n", argv[3]);
                    exit(1);
                }
            }

        }

        bsc_free (vb, NULL, buffer);
    }

    fprintf(stdout, "\r%.55s compressed %.0f into %.0f in %.3f seconds.\n", argv[2], (double)fileSize, (double)BSC_FTELL(fOutput), BSC_CLOCK() - startTime);

    fclose(fInput); fclose(fOutput);
}

void Decompression(char * argv[])
{
    FILE * fInput = fopen(argv[2], "rb");
    if (fInput == NULL)
    {
        fprintf(stderr, "Can't open input file: %s!\n", argv[2]);
        exit(1);
    }

    FILE * fOutput = fopen(argv[3], "wb");
    if (fOutput == NULL)
    {
        fprintf(stderr, "Can't create output file: %s!\n", argv[3]);
        exit(1);
    }

    if (BSC_FSEEK(fInput, 0, SEEK_END))
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
        exit(1);
    }

    BSC_FILEOFFSET fileSize = BSC_FTELL(fInput);
    if (fileSize < 0)
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
        exit(1);
    }

    if (BSC_FSEEK(fInput, 0, SEEK_SET))
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[2]);
        exit(1);
    }

    unsigned char inputFileSign[sizeof(bscFileSign)];

    if (fread(inputFileSign, sizeof(bscFileSign), 1, fInput) != 1)
    {
        fprintf(stderr, "This is not bsc archive!\n");
        exit(1);
    }

    if (memcmp(inputFileSign, bscFileSign, sizeof(bscFileSign)) != 0)
    {
        fprintf(stderr, "This is not bsc archive or invalid compression method!\n");
        exit(2);
    }

    int nBlocks = 0;
    if (fread(&nBlocks, sizeof(nBlocks), 1, fInput) != 1)
    {
        fprintf(stderr, "This is not bsc archive!\n");
        exit(1);
    }

    double startTime = BSC_CLOCK();

    {
        int bufferSize = -1; unsigned char * buffer = NULL;

        while (true)
        {
            BSC_FILEOFFSET  blockOffset     = 0;

            signed char     sortingContexts = 0;
            signed char     recordSize      = 0;
            int             blockSize       = 0;
            int             dataSize        = 0;

            {
                if ((feof(fInput) == 0) && (BSC_FTELL(fInput) != fileSize))
                {
                    {
                        double progress = (100.0 * (double)BSC_FTELL(fInput)) / fileSize;
                        fprintf(stdout, "\rDecompressing %.55s(%02d%%)", argv[2], (int)progress);
                        fflush(stdout);
                    }

                    BSC_BLOCK_HEADER header = {0, 0, 0};
                    if (fread(&header, sizeof(BSC_BLOCK_HEADER), 1, fInput) != 1)
                    {
                        fprintf(stderr, "\nUnexpected end of file: %s!\n", argv[2]);
                        exit(1);
                    }

                    recordSize = header.recordSize;
                    if (recordSize < 1)
                    {
                        fprintf(stderr, "\nThis is not bsc archive or invalid compression method!\n");
                        exit(2);
                    }

                    sortingContexts = header.sortingContexts;
                    if ((sortingContexts != LIBBSC_CONTEXTS_FOLLOWING) && (sortingContexts != LIBBSC_CONTEXTS_PRECEDING))
                    {
                        fprintf(stderr, "\nThis is not bsc archive or invalid compression method!\n");
                        exit(2);
                    }

                    blockOffset = (BSC_FILEOFFSET)header.blockOffset;

                    unsigned char bscBlockHeader[LIBBSC_HEADER_SIZE];

                    if (fread(bscBlockHeader, LIBBSC_HEADER_SIZE, 1, fInput) != 1)
                    {
                        fprintf(stderr, "\nUnexpected end of file: %s!\n", argv[2]);
                        exit(1);
                    }

                    if (bsc_block_info(NULL, bscBlockHeader, LIBBSC_HEADER_SIZE, &blockSize, &dataSize, paramFeatures()) != LIBBSC_NO_ERROR)
                    {
                        fprintf(stderr, "\nThis is not bsc archive or invalid compression method!\n");
                        exit(2);
                    }

                    if ((blockSize > bufferSize) || (dataSize > bufferSize))
                    {
                        if (blockSize > bufferSize) bufferSize = blockSize;
                        if (dataSize  > bufferSize) bufferSize = dataSize;

                        if (buffer != NULL) bsc_free (vb, NULL, buffer); 
                        buffer = (unsigned char *)bsc_malloc(NULL, bufferSize);
                    }

                    if (buffer == NULL)
                    {
                        fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n");
                        exit(2);
                    }

                    memcpy(buffer, bscBlockHeader, LIBBSC_HEADER_SIZE);

                    if (fread(buffer + LIBBSC_HEADER_SIZE, blockSize - LIBBSC_HEADER_SIZE, 1, fInput) != 1)
                    {
                        fprintf(stderr, "\nUnexpected end of file: %s!\n", argv[2]);
                        exit(1);
                    }
                }
            }

            if (dataSize == 0) break;

            int result = bsc_decompress(NULL, buffer, blockSize, buffer, dataSize, paramFeatures());
            if (result < LIBBSC_NO_ERROR)
            {
                {
                    switch (result)
                    {
                        case LIBBSC_DATA_CORRUPT            : fprintf(stderr, "\nThe compressed data is corrupted!\n"); break;
                        case LIBBSC_NOT_ENOUGH_MEMORY       : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                        case LIBBSC_GPU_ERROR               : fprintf(stderr, "\nGeneral GPU failure! Please check README file for more information.\n"); break;
                        case LIBBSC_GPU_NOT_SUPPORTED       : fprintf(stderr, "\nYour GPU is not supported! Please check README file for more information.\n"); break;
                        case LIBBSC_GPU_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough GPU memory! Please check README file for more information.\n"); break;

                        default                             : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                    }
                    exit(2);
                }
            }

            if (sortingContexts == LIBBSC_CONTEXTS_PRECEDING)
            {
                result = bsc_reverse_block(NULL, buffer, dataSize, paramFeatures());
                if (result != LIBBSC_NO_ERROR)
                {
                    {
                        fprintf(stderr, "\nInternal program error, please contact the author!\n");
                        exit(2);
                    }
                }
            }

            if (recordSize > 1)
            {
                result = bsc_reorder_reverse(NULL, buffer, dataSize, recordSize, paramFeatures());
                if (result != LIBBSC_NO_ERROR)
                {
                    {
                        switch (result)
                        {
                            case LIBBSC_NOT_ENOUGH_MEMORY   : fprintf(stderr, "\nNot enough memory! Please check README file for more information.\n"); break;
                            default                         : fprintf(stderr, "\nInternal program error, please contact the author!\n");
                        }
                        exit(2);
                    }
                }
            }
            {
                if (BSC_FSEEK(fOutput, blockOffset, SEEK_SET))
                {
                    fprintf(stderr, "\nIO error on file: %s!\n", argv[3]);
                    exit(1);
                }

                if ((int)fwrite(buffer, 1, dataSize, fOutput) != dataSize)
                {
                    fprintf(stderr, "\nIO error on file: %s!\n", argv[3]);
                    exit(1);
                }
            }
        }

        if (buffer != NULL) bsc_free (vb, NULL, buffer);
    }

    if (BSC_FSEEK(fOutput, 0, SEEK_END))
    {
        fprintf(stderr, "IO error on file: %s!\n", argv[3]);
        exit(1);
    }

    fprintf(stdout, "\r%.55s decompressed %.0f into %.0f in %.3f seconds.\n", argv[2], (double)fileSize, (double)BSC_FTELL(fOutput), BSC_CLOCK() - startTime);

    fclose(fInput); fclose(fOutput);
}

void ShowUsage(void)
{
    fprintf(stdout, "Usage: bsc <e|d> inputfile outputfile <options>\n\n");

    fprintf(stdout, "Block sorting options:\n");
    fprintf(stdout, "  -b<size> Block size in megabytes, default: -b25\n");
    fprintf(stdout, "             minimum: -b1, maximum: -b1024\n");
    fprintf(stdout, "  -m<algo> Block sorting algorithm, default: -m0\n");
    fprintf(stdout, "             -m0 Burrows Wheeler Transform (default)\n");
    fprintf(stdout, "  -c<ctx>  Contexts for sorting, default: -cf\n");
    fprintf(stdout, "             -cf Following contexts (default)\n");
    fprintf(stdout, "             -cp Preceding contexts\n");
    fprintf(stdout, "             -ca Autodetect (experimental)\n");
    fprintf(stdout, "  -e<algo> Entropy encoding algorithm, default: -e1\n");
    fprintf(stdout, "             -e1 Static Quantized Local Frequency Coding (default)\n");
    fprintf(stdout, "             -e2 Adaptive Quantized Local Frequency Coding (best compression)\n");
   
    fprintf(stdout, "\nPreprocessing options:\n");
    fprintf(stdout, "  -p       Disable all preprocessing techniques\n");
    fprintf(stdout, "  -s       Enable segmentation (adaptive block size), default: disable\n");
    fprintf(stdout, "  -r       Enable structured data reordering, default: disable\n");
    fprintf(stdout, "  -l       Enable Lempel-Ziv preprocessing, default: enable\n");
    fprintf(stdout, "  -H<size> LZP dictionary size in bits, default: -H16\n");
    fprintf(stdout, "             minimum: -H10, maximum: -H28\n");
    fprintf(stdout, "  -M<size> LZP minimum match length, default: -M128\n");
    fprintf(stdout, "             minimum: -M4, maximum: -M255\n");

    fprintf(stdout, "\nPlatform specific options:\n");

#ifdef _WIN32
    fprintf(stdout, "  -P       Enable large 2MB RAM pages, default: disable\n");
#endif

    fprintf(stdout,"\nOptions may be combined into one, like -b128p -m5e1\n");
    exit(0);
}

void ProcessSwitch(char * s)
{
    if (*s == 0)
    {
        ShowUsage();
    }

    for (; *s != 0; )
    {
        switch (*s++)
        {
            case 'b':
            {
                char * strNum = s; while ((*s >= '0') && (*s <= '9')) s++;
                paramBlockSize = atoi(strNum) * 1024 * 1024;
                if ((paramBlockSize < 1024 * 1024) || (paramBlockSize > 1024 * 1024 * 1024)) ShowUsage();
                break;
            }

            case 'm':
            {
                char * strNum = s; while ((*s >= '0') && (*s <= '9')) s++;
                switch (atoi(strNum))
                {
                    case 0   : paramBlockSorter = LIBBSC_BLOCKSORTER_BWT; break;

                    default  : ShowUsage();
                }
                break;
            }

            case 'c':
            {
                switch (*s++)
                {
                    case 'f' : paramSortingContexts = LIBBSC_CONTEXTS_FOLLOWING;  break;
                    case 'p' : paramSortingContexts = LIBBSC_CONTEXTS_PRECEDING;  break;
                    case 'a' : paramSortingContexts = LIBBSC_CONTEXTS_AUTODETECT; break;
                    default  : ShowUsage();
                }
                break;
            }

            case 'e':
            {
                switch (*s++)
                {
                    case '1' : paramCoder = LIBBSC_CODER_QLFC_STATIC;   break;
                    case '2' : paramCoder = LIBBSC_CODER_QLFC_ADAPTIVE; break;
                    default  : ShowUsage();
                }
                break;
            }

            case 'H':
            {
                char * strNum = s; while ((*s >= '0') && (*s <= '9')) s++;
                paramLZPHashSize = atoi(strNum);
                if ((paramLZPHashSize < 10) || (paramLZPHashSize > 28)) ShowUsage();
                break;
            }

            case 'M':
            {
                char * strNum = s; while ((*s >= '0') && (*s <= '9')) s++;
                paramLZPMinLen = atoi(strNum);
                if ((paramLZPMinLen < 4) || (paramLZPMinLen > 255)) ShowUsage();
                break;
            }

            case 'l': paramEnableLZP            = 1; break;
            case 's': paramEnableSegmentation   = 1; break;
            case 'r': paramEnableReordering     = 1; break;

            case 'p': paramEnableLZP = paramEnableSegmentation = paramEnableReordering = 0; break;

#ifdef _WIN32
            case 'P': paramEnableLargePages     = 1; break;
#endif

            default : ShowUsage();
        }
    }
}

void ProcessCommandline(int argc, char * argv[])
{
    if (argc < 4 || strlen(argv[1]) != 1)
    {
        ShowUsage();
    }

    for (int i = 4; i < argc; ++i)
    {
        if (argv[i][0] == '-')
        {
            ProcessSwitch(&argv[i][1]);
        }
        else
        {
            ShowUsage();
        }
    }
}

int main(int argc, char * argv[])
{
    fprintf(stdout, "This is bsc, Block Sorting Compressor. Version 3.1.0. 8 July 2012.\n");
    fprintf(stdout, "Copyright (c) 2009-2012 Ilya Grebnov <Ilya.Grebnov@gmail.com>.\n\n");

    ProcessCommandline(argc, argv);

    if (bsc_init_full (paramFeatures(), bsc_default_malloc, bsc_default_zero_malloc, bsc_default_free) != LIBBSC_NO_ERROR)
    {
        fprintf(stderr, "\nInternal program error, please contact the author!\n");
        exit(2);
    }

    switch (*argv[1])
    {
        case 'e' : case 'E' : Compression(argv); break;
        case 'd' : case 'D' : Decompression(argv); break;
        default  : ShowUsage();
    }

    return 0;
}

/*-----------------------------------------------------------*/
/* End                                               bsc.cpp */
/*-----------------------------------------------------------*/
