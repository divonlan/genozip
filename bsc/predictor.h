/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Probability counter and logistic mixer                    */
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
//   All modifications:
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef _LIBBSC_CODER_PREDICTOR_H
#define _LIBBSC_CODER_PREDICTOR_H

#include "platform.h"
#include "tables.h"

static inline void ProbCounter_UpdateBit0(short *probability, const int threshold, const int adaptationRate)
{
    *probability += (((4096 - threshold - *probability) * adaptationRate) >> 12);
};

static inline void ProbCounter_UpdateBit1(short *probability, const int threshold, const int adaptationRate)
{
    *probability -= (((*probability - threshold) * adaptationRate) >> 12);
};

typedef struct {
    short stretchedProbability0;
    short stretchedProbability1;
    short stretchedProbability2;
    int   mixedProbability;
    int   index;

    short probabilityMap[17];

    int weight0;
    int weight1;
    int weight2;
} ProbabilityMixer;

static inline void MixerInit(ProbabilityMixer *mixer)
{
    mixer->weight0 = mixer->weight1 = 2048 << 5; mixer->weight2 = 0;
    for (int p = 0; p < 17; ++p)
    {
        mixer->probabilityMap[p] = bsc_squash((p - 8) * 256);
    }
}

static inline int Mixup(ProbabilityMixer *mixer, const int probability0, const int probability1, const int probability2)
{
    mixer->stretchedProbability0 = bsc_stretch(probability0);
    mixer->stretchedProbability1 = bsc_stretch(probability1);
    mixer->stretchedProbability2 = bsc_stretch(probability2);

    short stretchedProbability = (mixer->stretchedProbability0 * mixer->weight0 + 
                                  mixer->stretchedProbability1 * mixer->weight1 + 
                                  mixer->stretchedProbability2 * mixer->weight2) >> 17;

    if (stretchedProbability < -2047) stretchedProbability = -2047;
    if (stretchedProbability >  2047) stretchedProbability =  2047;

    mixer->index                = (stretchedProbability + 2048) >> 8;
    const int weight            = stretchedProbability & 255;
    const int probability       = bsc_squash(stretchedProbability);
    const int mappedProbability = mixer->probabilityMap[mixer->index] + (((mixer->probabilityMap[mixer->index + 1] - mixer->probabilityMap[mixer->index]) * weight) >> 8);

    return mixer->mixedProbability = (3 * probability + mappedProbability) >> 2;
};

static inline int MixupAndUpdateBit0(ProbabilityMixer * restrict mixer, 
                                     const int probability0,  const int probability1,  const int probability2,
                                     const int learningRate0, const int learningRate1, const int learningRate2,
                                     const int threshold,     const int adaptationRate
)
{
    const short stretchedProbability0 = bsc_stretch(probability0);
    const short stretchedProbability1 = bsc_stretch(probability1);
    const short stretchedProbability2 = bsc_stretch(probability2);

    short stretchedProbability = (stretchedProbability0 * mixer->weight0 + stretchedProbability1 * mixer->weight1 + stretchedProbability2 * mixer->weight2) >> 17;

    if (stretchedProbability < -2047) stretchedProbability = -2047;
    if (stretchedProbability >  2047) stretchedProbability =  2047;

    const int weight            = stretchedProbability & 255;
    const int index             = (stretchedProbability + 2048) >> 8;
    const int probability       = bsc_squash(stretchedProbability);
    const int mappedProbability = mixer->probabilityMap[index] + (((mixer->probabilityMap[index + 1] - mixer->probabilityMap[index]) * weight) >> 8);
    const int mixedProbability  = (3 * probability + mappedProbability) >> 2;

    ProbCounter_UpdateBit0(&mixer->probabilityMap[index], threshold, adaptationRate);
    ProbCounter_UpdateBit0(&mixer->probabilityMap[index + 1], threshold, adaptationRate);

    const int eps = mixedProbability - 4095;

    mixer->weight0 -= (learningRate0 * eps * stretchedProbability0) >> 16;
    mixer->weight1 -= (learningRate1 * eps * stretchedProbability1) >> 16;
    mixer->weight2 -= (learningRate2 * eps * stretchedProbability2) >> 16;

    return mixedProbability;
};

static inline int MixupAndUpdateBit1(ProbabilityMixer *mixer, 
                                     const int probability0,  const int probability1,  const int probability2,
                                     const int learningRate0, const int learningRate1, const int learningRate2,
                                     const int threshold,     const int adaptationRate
)
{
    const short stretchedProbability0 = bsc_stretch(probability0);
    const short stretchedProbability1 = bsc_stretch(probability1);
    const short stretchedProbability2 = bsc_stretch(probability2);

    short stretchedProbability = (stretchedProbability0 * mixer->weight0 + stretchedProbability1 * mixer->weight1 + stretchedProbability2 * mixer->weight2) >> 17;

    if (stretchedProbability < -2047) stretchedProbability = -2047;
    if (stretchedProbability >  2047) stretchedProbability =  2047;

    const int weight            = stretchedProbability & 255;
    const int index             = (stretchedProbability + 2048) >> 8;
    const int probability       = bsc_squash(stretchedProbability);
    const int mappedProbability = mixer->probabilityMap[index] + (((mixer->probabilityMap[index + 1] - mixer->probabilityMap[index]) * weight) >> 8);
    const int mixedProbability  = (3 * probability + mappedProbability) >> 2;

    ProbCounter_UpdateBit1(&mixer->probabilityMap[index], threshold, adaptationRate);
    ProbCounter_UpdateBit1(&mixer->probabilityMap[index + 1], threshold, adaptationRate);

    const int eps = mixedProbability - 1;

    mixer->weight0 -= (learningRate0 * eps * stretchedProbability0) >> 16;
    mixer->weight1 -= (learningRate1 * eps * stretchedProbability1) >> 16;
    mixer->weight2 -= (learningRate2 * eps * stretchedProbability2) >> 16;

    return mixedProbability;
};

static inline void UpdateBit0(ProbabilityMixer *mixer, 
                              const int learningRate0, const int learningRate1, const int learningRate2,
                              const int threshold,     const int adaptationRate
)
{
    ProbCounter_UpdateBit0(&mixer->probabilityMap[mixer->index], threshold, adaptationRate);
    ProbCounter_UpdateBit0(&mixer->probabilityMap[mixer->index + 1], threshold, adaptationRate);

    const int eps = mixer->mixedProbability - 4095;

    mixer->weight0 -= (learningRate0 * eps * mixer->stretchedProbability0) >> 16;
    mixer->weight1 -= (learningRate1 * eps * mixer->stretchedProbability1) >> 16;
    mixer->weight2 -= (learningRate2 * eps * mixer->stretchedProbability2) >> 16;
};

static inline void UpdateBit1(ProbabilityMixer *mixer, 
                              const int learningRate0, const int learningRate1, const int learningRate2,
                              const int threshold,     const int adaptationRate
)
{
    ProbCounter_UpdateBit1(&mixer->probabilityMap[mixer->index], threshold, adaptationRate);
    ProbCounter_UpdateBit1(&mixer->probabilityMap[mixer->index + 1], threshold, adaptationRate);

    const int eps = mixer->mixedProbability - 1;

    mixer->weight0 -= (learningRate0 * eps * mixer->stretchedProbability0) >> 16;
    mixer->weight1 -= (learningRate1 * eps * mixer->stretchedProbability1) >> 16;
    mixer->weight2 -= (learningRate2 * eps * mixer->stretchedProbability2) >> 16;
};

#endif
