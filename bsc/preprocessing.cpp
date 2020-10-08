/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Data preprocessing functions                              */
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

#include <stdlib.h>
#include <memory.h>

#include "filters.h"
#include "platform.h"
#include "libbsc.h"

int bsc_reverse_block(unsigned char * T, int n, int features)
{
    for (int i = 0, j = n - 1; i < j; ++i, --j)
    {
        unsigned char tmp = T[i]; T[i] = T[j]; T[j] = tmp;
    }
 
    return LIBBSC_NO_ERROR;
}

int bsc_reorder_forward (void *vb, unsigned char * T, int n, int recordSize, int features)
{
    if (recordSize <= 0) return LIBBSC_BAD_PARAMETER;
    if (recordSize == 1) return LIBBSC_NO_ERROR;

    if (unsigned char * buffer = (unsigned char *)bsc_malloc (vb, n))
    {
        memcpy(buffer, T, n);

        unsigned char * RESTRICT S = buffer;
        unsigned char * RESTRICT D = T;

        int chunk = (n / recordSize);

        switch (recordSize)
        {
            case 2: for (int i = 0; i < chunk; ++i) { D[0] = S[0]; D[chunk] = S[1]; D++; S += 2; } break;
            case 3: for (int i = 0; i < chunk; ++i) { D[0] = S[0]; D[chunk] = S[1]; D[chunk * 2] = S[2]; D++; S += 3; } break;
            case 4: for (int i = 0; i < chunk; ++i) { D[0] = S[0]; D[chunk] = S[1]; D[chunk * 2] = S[2]; D[chunk * 3] = S[3]; D++; S += 4; } break;
            default:
                for (int i = 0; i < chunk; ++i) { for (int j = 0; j < recordSize; ++j) D[j * chunk] = S[j]; D++; S += recordSize; }
        }

        bsc_free (vb, buffer); return LIBBSC_NO_ERROR;
    }

    return LIBBSC_NOT_ENOUGH_MEMORY;
}

int bsc_reorder_reverse (void *vb, unsigned char * T, int n, int recordSize, int features)
{
    if (recordSize <= 0) return LIBBSC_BAD_PARAMETER;
    if (recordSize == 1) return LIBBSC_NO_ERROR;

    if (unsigned char * buffer = (unsigned char *)bsc_malloc (vb, n))
    {
        memcpy(buffer, T, n);

        unsigned char * RESTRICT S = buffer;
        unsigned char * RESTRICT D = T;

        int chunk = (n / recordSize);

        switch (recordSize)
        {
            case 2: for (int i = 0; i < chunk; ++i) { D[0] = S[0]; D[1] = S[chunk]; D += 2; S++; } break;
            case 3: for (int i = 0; i < chunk; ++i) { D[0] = S[0]; D[1] = S[chunk]; D[2] = S[chunk * 2]; D += 3; S++; } break;
            case 4: for (int i = 0; i < chunk; ++i) { D[0] = S[0]; D[1] = S[chunk]; D[2] = S[chunk * 2]; D[3] = S[chunk * 3]; D += 4; S++; } break;
            default:
                for (int i = 0; i < chunk; ++i) { for (int j = 0; j < recordSize; ++j) D[j] = S[j * chunk]; D += recordSize; S++; }
        }

        bsc_free (vb, buffer); return LIBBSC_NO_ERROR;
    }

    return LIBBSC_NOT_ENOUGH_MEMORY;
}

/*-------------------------------------------------*/
/* End                           preprocessing.cpp */
/*-------------------------------------------------*/
