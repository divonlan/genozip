/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Quantized Local Frequency Coding functions                */
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
//   This file was extensively modified to adapt it to genozip. 

#include <stdlib.h>
#include <memory.h>

#include "qlfc.h"
#include "libbsc.h"
#include "platform.h"
#include "rangecoder.h"
#include "tables.h"
#include "predictor.h"
#include "qlfc_model.h"

int bsc_qlfc_init(int features)
{
    return bsc_qlfc_init_static_model();
}

static unsigned char * bsc_qlfc_transform(void *vb, const unsigned char * restrict input, unsigned char * restrict buffer, int n, unsigned char * restrict MTFTable)
{
    unsigned char Flag[ALPHABET_SIZE];

    for (int i = 0; i < ALPHABET_SIZE; ++i) Flag[i] = 0;
    for (int i = 0; i < ALPHABET_SIZE; ++i) MTFTable[i] = i;

    if (input[n - 1] == 0)
    {
        MTFTable[0] = 1; MTFTable[1] = 0;
    }

    int index = n, nSymbols = 0;
    for (int i = n - 1; i >= 0;)
    {
        unsigned char currentChar = input[i--];
        for (; (i >= 0) && (input[i] == currentChar); --i) ;

        unsigned char previousChar = MTFTable[0], rank = 1; MTFTable[0] = currentChar;
        while (true)
        {
            unsigned char temporaryChar0 = MTFTable[rank + 0]; MTFTable[rank + 0] = previousChar;
            if (temporaryChar0 == currentChar) { rank += 0; break; }

            unsigned char temporaryChar1 = MTFTable[rank + 1]; MTFTable[rank + 1] = temporaryChar0;
            if (temporaryChar1 == currentChar) { rank += 1; break; }

            unsigned char temporaryChar2 = MTFTable[rank + 2]; MTFTable[rank + 2] = temporaryChar1;
            if (temporaryChar2 == currentChar) { rank += 2; break; }

            unsigned char temporaryChar3 = MTFTable[rank + 3]; MTFTable[rank + 3] = temporaryChar2;
            if (temporaryChar3 == currentChar) { rank += 3; break; }

            rank += 4; previousChar = temporaryChar3;
        }

        if (Flag[currentChar] == 0)
        {
            Flag[currentChar] = 1;
            rank = nSymbols++;
        }

        buffer[--index] = rank;
    }

    buffer[n - 1] = 1;

    for (int rank = 1; rank < ALPHABET_SIZE; ++rank)
    {
        if (Flag[MTFTable[rank]] == 0)
        {
            MTFTable[rank] = MTFTable[rank - 1];
            break;
        }
    }

    return buffer + index;
}

static int bsc_qlfc_adaptive_encode(void *vb, const unsigned char * input, unsigned char * output, unsigned char * buffer, int inputSize, int outputSize, QlfcStatisticalModel * model)
{
    unsigned char MTFTable[ALPHABET_SIZE];

    bsc_qlfc_init_model(model);

    int contextRank0 = 0;
    int contextRank4 = 0;
    int contextRun   = 0;
    int maxRank      = 7;
    int avgRank      = 0;

    unsigned char rankHistory[ALPHABET_SIZE], runHistory[ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; ++i)
    {
        rankHistory[i] = runHistory[i] = 0;
    }

    unsigned char * rankArray = bsc_qlfc_transform(vb, input, buffer, inputSize, MTFTable);

    RangeCoder coder;

    InitEncoder(&coder, output, outputSize);
    EncodeWord(&coder, (unsigned int)inputSize);

    unsigned char usedChar[ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; ++i) usedChar[i] = 0;

    int prevChar = -1;
    for (int rank = 0; rank < ALPHABET_SIZE; ++rank)
    {
        int currentChar = MTFTable[rank];

        for (int bit = 7; bit >= 0; --bit)
        {
            bool bit0 = false, bit1 = false;

            for (int c = 0; c < ALPHABET_SIZE; ++c)
            {
                if (c == prevChar || usedChar[c] == 0)
                {
                    if ((currentChar >> (bit + 1)) == (c >> (bit + 1)))
                    {
                        if (c & (1 << bit)) bit1 = true; else bit0 = true;
                        if (bit0 && bit1) break;
                    }
                }
            }

            if (bit0 && bit1)
            {
                EncodeBit(&coder, currentChar & (1 << bit));
            }
        }

        if (currentChar == prevChar)
        {
            maxRank = bsc_log2_256(rank - 1);
            break;
        }

        prevChar = currentChar; usedChar[currentChar] = 1;
    }

    for (const unsigned char * inputEnd = input + inputSize; input < inputEnd;)
    {
        if (CheckEOB(&coder))
        {
            return LIBBSC_NOT_COMPRESSIBLE;
        }

        int currentChar = *input, runSize;
        {
            const unsigned char * inputStart = input++;
            while (true)
            {
                if (input <= inputEnd - 4)
                {
                    if (input[0] != currentChar) { input += 0; break; }
                    if (input[1] != currentChar) { input += 1; break; }
                    if (input[2] != currentChar) { input += 2; break; }
                    if (input[3] != currentChar) { input += 3; break; }

                    input += 4;
                }
                else
                {
                    while ((input < inputEnd) && (*input == currentChar)) ++input;
                    break;
                }
            }

            runSize = (int)(input - inputStart);
        }

        int                 rank            =   *rankArray++;
        int                 history         =   rankHistory[currentChar];
        int                 state           =   model_rank_state(contextRank4, contextRun, history);

        short *            restrict statePredictor  = & model->Rank.StateModel[state];
        short *            restrict charPredictor   = & model->Rank.CharModel[currentChar];
        short *            restrict staticPredictor = & model->Rank.StaticModel;
        ProbabilityMixer * restrict mixer           = & model->mixerOfRank[currentChar];

        if (avgRank < 32)
        {
            if (rank == 1)
            {
                rankHistory[currentChar] = 0;

                int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                ProbCounter_UpdateBit0(statePredictor,  M_RANK_TS_TH0, M_RANK_TS_AR0);
                ProbCounter_UpdateBit0(charPredictor,   M_RANK_TC_TH0, M_RANK_TC_AR0);
                ProbCounter_UpdateBit0(staticPredictor, M_RANK_TP_TH0, M_RANK_TP_AR0);

                EncodeBit0(&coder, MixupAndUpdateBit0 (mixer, probability0, probability1, probability2, M_RANK_TM_LR0, M_RANK_TM_LR1, M_RANK_TM_LR2, M_RANK_TM_TH0, M_RANK_TM_AR0));
            }
            else
            {
                {
                    int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                    ProbCounter_UpdateBit1(statePredictor,  M_RANK_TS_TH1, M_RANK_TS_AR1);
                    ProbCounter_UpdateBit1(charPredictor,   M_RANK_TC_TH1, M_RANK_TC_AR1);
                    ProbCounter_UpdateBit1(staticPredictor, M_RANK_TP_TH1, M_RANK_TP_AR1);

                    EncodeBit1(&coder, MixupAndUpdateBit1 (mixer, probability0, probability1, probability2, M_RANK_TM_LR0, M_RANK_TM_LR1, M_RANK_TM_LR2, M_RANK_TM_TH1, M_RANK_TM_AR1));
                }

                int bitRankSize = bsc_log2_256(rank); rankHistory[currentChar] = bitRankSize;

                statePredictor  = & model->Rank.Exponent.StateModel[state][0];
                charPredictor   = & model->Rank.Exponent.CharModel[currentChar][0];
                staticPredictor = & model->Rank.Exponent.StaticModel[0];
                mixer           = & model->mixerOfRankExponent[history < 1 ? 1 : history][1];

                for (int bit = 1; bit < bitRankSize; ++bit, ++statePredictor, ++charPredictor, ++staticPredictor)
                {
                    int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                    ProbCounter_UpdateBit1(statePredictor,  M_RANK_ES_TH1, M_RANK_ES_AR1);
                    ProbCounter_UpdateBit1(charPredictor,   M_RANK_EC_TH1, M_RANK_EC_AR1);
                    ProbCounter_UpdateBit1(staticPredictor, M_RANK_EP_TH1, M_RANK_EP_AR1);

                    EncodeBit1(&coder, MixupAndUpdateBit1 (mixer, probability0, probability1, probability2, M_RANK_EM_LR0, M_RANK_EM_LR1, M_RANK_EM_LR2, M_RANK_EM_TH1, M_RANK_EM_AR1));

                    mixer = & model->mixerOfRankExponent[history <= bit ? bit + 1 : history][bit + 1];
                }
                if (bitRankSize < maxRank)
                {
                    int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                    ProbCounter_UpdateBit0(statePredictor,  M_RANK_ES_TH0, M_RANK_ES_AR0);
                    ProbCounter_UpdateBit0(charPredictor,   M_RANK_EC_TH0, M_RANK_EC_AR0);
                    ProbCounter_UpdateBit0(staticPredictor, M_RANK_EP_TH0, M_RANK_EP_AR0);

                    EncodeBit0(&coder, MixupAndUpdateBit0 (mixer, probability0, probability1, probability2, M_RANK_EM_LR0, M_RANK_EM_LR1, M_RANK_EM_LR2, M_RANK_EM_TH0, M_RANK_EM_AR0));
                }

                statePredictor  = & model->Rank.Mantissa[bitRankSize].StateModel[state][0];
                charPredictor   = & model->Rank.Mantissa[bitRankSize].CharModel[currentChar][0];
                staticPredictor = & model->Rank.Mantissa[bitRankSize].StaticModel[0];
                mixer           = & model->mixerOfRankMantissa[bitRankSize];

                for (int context = 1, bit = bitRankSize - 1; bit >= 0; --bit)
                {
                    if (rank & (1 << bit))
                    {
                        int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                        ProbCounter_UpdateBit1(&statePredictor[context],  M_RANK_MS_TH1, M_RANK_MS_AR1);
                        ProbCounter_UpdateBit1(&charPredictor[context],   M_RANK_MC_TH1, M_RANK_MC_AR1);
                        ProbCounter_UpdateBit1(&staticPredictor[context], M_RANK_MP_TH1, M_RANK_MP_AR1);

                        EncodeBit1(&coder, MixupAndUpdateBit1 (mixer, probability0, probability1, probability2, M_RANK_MM_LR0, M_RANK_MM_LR1, M_RANK_MM_LR2, M_RANK_MM_TH1, M_RANK_MM_AR1));

                        context += context + 1;
                    }
                    else
                    {
                        int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                        ProbCounter_UpdateBit0(&statePredictor[context],  M_RANK_MS_TH0, M_RANK_MS_AR0);
                        ProbCounter_UpdateBit0(&charPredictor[context],   M_RANK_MC_TH0, M_RANK_MC_AR0);
                        ProbCounter_UpdateBit0(&staticPredictor[context], M_RANK_MP_TH0, M_RANK_MP_AR0);

                        EncodeBit0(&coder, MixupAndUpdateBit0 (mixer, probability0, probability1, probability2, M_RANK_MM_LR0, M_RANK_MM_LR1, M_RANK_MM_LR2, M_RANK_MM_TH0, M_RANK_MM_AR0));

                        context += context;
                    }
                }
            }
        }
        else
        {
            rankHistory[currentChar] = bsc_log2_256(rank);

            statePredictor  = & model->Rank.Escape.StateModel[state][0];
            charPredictor   = & model->Rank.Escape.CharModel[currentChar][0];
            staticPredictor = & model->Rank.Escape.StaticModel[0];

            for (int context = 1, bit = maxRank; bit >= 0; --bit)
            {
                mixer = & model->mixerOfRankEscape[context];

                if (rank & (1 << bit))
                {
                    int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                    ProbCounter_UpdateBit1(&statePredictor[context],  M_RANK_PS_TH1, M_RANK_PS_AR1);
                    ProbCounter_UpdateBit1(&charPredictor[context],   M_RANK_PC_TH1, M_RANK_PC_AR1);
                    ProbCounter_UpdateBit1(&staticPredictor[context], M_RANK_PP_TH1, M_RANK_PP_AR1);

                    EncodeBit1(&coder, MixupAndUpdateBit1 (mixer, probability0, probability1, probability2, M_RANK_PM_LR0, M_RANK_PM_LR1, M_RANK_PM_LR2, M_RANK_PM_TH1, M_RANK_PM_AR1));

                    context += context + 1;
                }
                else
                {
                    int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                    ProbCounter_UpdateBit0(&statePredictor[context],  M_RANK_PS_TH0, M_RANK_PS_AR0);
                    ProbCounter_UpdateBit0(&charPredictor[context],   M_RANK_PC_TH0, M_RANK_PC_AR0);
                    ProbCounter_UpdateBit0(&staticPredictor[context], M_RANK_PP_TH0, M_RANK_PP_AR0);

                    EncodeBit0(&coder, MixupAndUpdateBit0 (mixer, probability0, probability1, probability2, M_RANK_PM_LR0, M_RANK_PM_LR1, M_RANK_PM_LR2, M_RANK_PM_TH0, M_RANK_PM_AR0));

                    context += context;
                }
            }
        }

        avgRank         =   (avgRank * 124 + rank * 4) >> 7;
        rank            =   rank - 1;
        history         =   runHistory[currentChar];
        state           =   model_run_state(contextRank0, contextRun, rank, history);
        statePredictor  = & model->Run.StateModel[state];
        charPredictor   = & model->Run.CharModel[currentChar];
        staticPredictor = & model->Run.StaticModel;
        mixer           = & model->mixerOfRun[currentChar];

        if (runSize == 1)
        {
            runHistory[currentChar] = (runHistory[currentChar] + 2) >> 2;

            int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

            ProbCounter_UpdateBit0(statePredictor,  M_RUN_TS_TH0, M_RUN_TS_AR0);
            ProbCounter_UpdateBit0(charPredictor,   M_RUN_TC_TH0, M_RUN_TC_AR0);
            ProbCounter_UpdateBit0(staticPredictor, M_RUN_TP_TH0, M_RUN_TP_AR0);

            EncodeBit0(&coder, MixupAndUpdateBit0 (mixer, probability0, probability1, probability2, M_RUN_TM_LR0, M_RUN_TM_LR1, M_RUN_TM_LR2, M_RUN_TM_TH0, M_RUN_TM_AR0));
        }
        else
        {
            {
                int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                ProbCounter_UpdateBit1(statePredictor,  M_RUN_TS_TH1, M_RUN_TS_AR1);
                ProbCounter_UpdateBit1(charPredictor,   M_RUN_TC_TH1, M_RUN_TC_AR1);
                ProbCounter_UpdateBit1(staticPredictor, M_RUN_TP_TH1, M_RUN_TP_AR1);

                EncodeBit1(&coder, MixupAndUpdateBit1 (mixer, probability0, probability1, probability2, M_RUN_TM_LR0, M_RUN_TM_LR1, M_RUN_TM_LR2, M_RUN_TM_TH1, M_RUN_TM_AR1));
            }

            int bitRunSize = bsc_log2(runSize); runHistory[currentChar] = (runHistory[currentChar] + 3 * bitRunSize + 3) >> 2;

            statePredictor  = & model->Run.Exponent.StateModel[state][0];
            charPredictor   = & model->Run.Exponent.CharModel[currentChar][0];
            staticPredictor = & model->Run.Exponent.StaticModel[0];
            mixer           = & model->mixerOfRunExponent[history < 1 ? 1 : history][1];

            for (int bit = 1; bit < bitRunSize; ++bit, ++statePredictor, ++charPredictor, ++staticPredictor)
            {
                int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                ProbCounter_UpdateBit1(statePredictor,  M_RUN_ES_TH1, M_RUN_ES_AR1);
                ProbCounter_UpdateBit1(charPredictor,   M_RUN_EC_TH1, M_RUN_EC_AR1);
                ProbCounter_UpdateBit1(staticPredictor, M_RUN_EP_TH1, M_RUN_EP_AR1);

                EncodeBit1(&coder, MixupAndUpdateBit1 (mixer, probability0, probability1, probability2, M_RUN_EM_LR0, M_RUN_EM_LR1, M_RUN_EM_LR2, M_RUN_EM_TH1, M_RUN_EM_AR1));

                mixer = & model->mixerOfRunExponent[history <= bit ? bit + 1 : history][bit + 1];
            }
            {
                int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                ProbCounter_UpdateBit0(statePredictor,  M_RUN_ES_TH0, M_RUN_ES_AR0);
                ProbCounter_UpdateBit0(charPredictor,   M_RUN_EC_TH0, M_RUN_EC_AR0);
                ProbCounter_UpdateBit0(staticPredictor, M_RUN_EP_TH0, M_RUN_EP_AR0);

                EncodeBit0(&coder, MixupAndUpdateBit0 (mixer, probability0, probability1, probability2, M_RUN_EM_LR0, M_RUN_EM_LR1, M_RUN_EM_LR2, M_RUN_EM_TH0, M_RUN_EM_AR0));
            }

            statePredictor  = & model->Run.Mantissa[bitRunSize].StateModel[state][0];
            charPredictor   = & model->Run.Mantissa[bitRunSize].CharModel[currentChar][0];
            staticPredictor = & model->Run.Mantissa[bitRunSize].StaticModel[0];
            mixer           = & model->mixerOfRunMantissa[bitRunSize];

            for (int context = 1, bit = bitRunSize - 1; bit >= 0; --bit)
            {
                if (runSize & (1 << bit))
                {
                    int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                    ProbCounter_UpdateBit1(&statePredictor[context],  M_RUN_MS_TH1, M_RUN_MS_AR1);
                    ProbCounter_UpdateBit1(&charPredictor[context],   M_RUN_MC_TH1, M_RUN_MC_AR1);
                    ProbCounter_UpdateBit1(&staticPredictor[context], M_RUN_MP_TH1, M_RUN_MP_AR1);

                    EncodeBit1(&coder, MixupAndUpdateBit1 (mixer, probability0, probability1, probability2, M_RUN_MM_LR0, M_RUN_MM_LR1, M_RUN_MM_LR2, M_RUN_MM_TH1, M_RUN_MM_AR1));

                    if (bitRunSize <= 5) context += context + 1; else context++;
                }
                else
                {
                    int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                    ProbCounter_UpdateBit0(&statePredictor[context],  M_RUN_MS_TH0, M_RUN_MS_AR0);
                    ProbCounter_UpdateBit0(&charPredictor[context],   M_RUN_MC_TH0, M_RUN_MC_AR0);
                    ProbCounter_UpdateBit0(&staticPredictor[context], M_RUN_MP_TH0, M_RUN_MP_AR0);

                    EncodeBit0(&coder, MixupAndUpdateBit0 (mixer, probability0, probability1, probability2, M_RUN_MM_LR0, M_RUN_MM_LR1, M_RUN_MM_LR2, M_RUN_MM_TH0, M_RUN_MM_AR0));

                    if (bitRunSize <= 5) context += context + 0; else context++;
                }
            }
        }

        contextRank0 = ((contextRank0 << 1) | (rank == 0   ? 1    : 0)) & 0x7;
        contextRank4 = ((contextRank4 << 2) | (rank < 3    ? rank : 3)) & 0xff;
        contextRun   = ((contextRun   << 1) | (runSize < 3 ? 1    : 0)) & 0xf;
    }

    return FinishEncoder(&coder);
}

static int bsc_qlfc_static_encode(void *vb, const unsigned char * input, unsigned char * output, unsigned char * buffer, int inputSize, int outputSize, QlfcStatisticalModel * model)
{
    unsigned char MTFTable[ALPHABET_SIZE];

    bsc_qlfc_init_model(model);

    int contextRank0 = 0;
    int contextRank4 = 0;
    int contextRun   = 0;
    int maxRank      = 7;
    int avgRank      = 0;

    unsigned char rankHistory[ALPHABET_SIZE], runHistory[ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; ++i)
    {
        rankHistory[i] = runHistory[i] = 0;
    }

    unsigned char * rankArray = bsc_qlfc_transform (vb, input, buffer, inputSize, MTFTable);

    RangeCoder coder;

    InitEncoder(&coder, output, outputSize);
    EncodeWord(&coder, (unsigned int)inputSize);

    unsigned char usedChar[ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; ++i) usedChar[i] = 0;

    int prevChar = -1;
    for (int rank = 0; rank < ALPHABET_SIZE; ++rank)
    {
        int currentChar = MTFTable[rank];

        for (int bit = 7; bit >= 0; --bit)
        {
            bool bit0 = false, bit1 = false;

            for (int c = 0; c < ALPHABET_SIZE; ++c)
            {
                if (c == prevChar || usedChar[c] == 0)
                {
                    if ((currentChar >> (bit + 1)) == (c >> (bit + 1)))
                    {
                        if (c & (1 << bit)) bit1 = true; else bit0 = true;
                        if (bit0 && bit1) break;
                    }
                }
            }

            if (bit0 && bit1)
            {
                EncodeBit(&coder,currentChar & (1 << bit));
            }
        }

        if (currentChar == prevChar)
        {
            maxRank = bsc_log2_256(rank - 1);
            break;
        }

        prevChar = currentChar; usedChar[currentChar] = 1;
    }

    for (const unsigned char * inputEnd = input + inputSize; input < inputEnd;)
    {
        if (CheckEOB(&coder))
        {
            return LIBBSC_NOT_COMPRESSIBLE;
        }

        int currentChar = *input, runSize;
        {
            const unsigned char * inputStart = input++;
            while (true)
            {
                if (input <= inputEnd - 4)
                {
                    if (input[0] != currentChar) { input += 0; break; }
                    if (input[1] != currentChar) { input += 1; break; }
                    if (input[2] != currentChar) { input += 2; break; }
                    if (input[3] != currentChar) { input += 3; break; }

                    input += 4;
                }
                else
                {
                    while ((input < inputEnd) && (*input == currentChar)) ++input;
                    break;
                }
            }

            runSize = (int)(input - inputStart);
        }

        int                 rank            =   *rankArray++;
        int                 history         =   rankHistory[currentChar];
        int                 state           =   model_rank_state(contextRank4, contextRun, history);

        short * restrict    statePredictor  = & model->Rank.StateModel[state];
        short * restrict    charPredictor   = & model->Rank.CharModel[currentChar];
        short * restrict    staticPredictor = & model->Rank.StaticModel;

        if (avgRank < 32)
        {
            if (rank == 1)
            {
                rankHistory[currentChar] = 0;

                int probability = ((*charPredictor) * F_RANK_TM_LR0 + (*statePredictor) * F_RANK_TM_LR1 + (*staticPredictor) * F_RANK_TM_LR2) >> 5;

                ProbCounter_UpdateBit0(statePredictor,  F_RANK_TS_TH0, F_RANK_TS_AR0);
                ProbCounter_UpdateBit0(charPredictor,   F_RANK_TC_TH0, F_RANK_TC_AR0);
                ProbCounter_UpdateBit0(staticPredictor, F_RANK_TP_TH0, F_RANK_TP_AR0);

                EncodeBit0(&coder,probability);
            }
            else
            {
                {
                    int probability = ((*charPredictor) * F_RANK_TM_LR0 + (*statePredictor) * F_RANK_TM_LR1 + (*staticPredictor) * F_RANK_TM_LR2) >> 5;

                    ProbCounter_UpdateBit1(statePredictor,  F_RANK_TS_TH1, F_RANK_TS_AR1);
                    ProbCounter_UpdateBit1(charPredictor,   F_RANK_TC_TH1, F_RANK_TC_AR1);
                    ProbCounter_UpdateBit1(staticPredictor, F_RANK_TP_TH1, F_RANK_TP_AR1);

                    EncodeBit1(&coder,probability);
                }

                int bitRankSize = bsc_log2_256(rank); rankHistory[currentChar] = bitRankSize;

                statePredictor  = & model->Rank.Exponent.StateModel[state][0];
                charPredictor   = & model->Rank.Exponent.CharModel[currentChar][0];
                staticPredictor = & model->Rank.Exponent.StaticModel[0];

                for (int bit = 1; bit < bitRankSize; ++bit, ++statePredictor, ++charPredictor, ++staticPredictor)
                {
                    int probability = ((*charPredictor) * F_RANK_EM_LR0 + (*statePredictor) * F_RANK_EM_LR1 + (*staticPredictor) * F_RANK_EM_LR2) >> 5;

                    ProbCounter_UpdateBit1(statePredictor,  F_RANK_ES_TH1, F_RANK_ES_AR1);
                    ProbCounter_UpdateBit1(charPredictor,   F_RANK_EC_TH1, F_RANK_EC_AR1);
                    ProbCounter_UpdateBit1(staticPredictor, F_RANK_EP_TH1, F_RANK_EP_AR1);

                    EncodeBit1(&coder,probability);
                }
                if (bitRankSize < maxRank)
                {
                    int probability = ((*charPredictor) * F_RANK_EM_LR0 + (*statePredictor) * F_RANK_EM_LR1 + (*staticPredictor) * F_RANK_EM_LR2) >> 5;

                    ProbCounter_UpdateBit0(statePredictor,  F_RANK_ES_TH0, F_RANK_ES_AR0);
                    ProbCounter_UpdateBit0(charPredictor,   F_RANK_EC_TH0, F_RANK_EC_AR0);
                    ProbCounter_UpdateBit0(staticPredictor, F_RANK_EP_TH0, F_RANK_EP_AR0);

                    EncodeBit0(&coder,probability);
                }

                statePredictor  = & model->Rank.Mantissa[bitRankSize].StateModel[state][0];
                charPredictor   = & model->Rank.Mantissa[bitRankSize].CharModel[currentChar][0];
                staticPredictor = & model->Rank.Mantissa[bitRankSize].StaticModel[0];

                for (int context = 1, bit = bitRankSize - 1; bit >= 0; --bit)
                {
                    int probability = (charPredictor[context] * F_RANK_MM_LR0 + statePredictor[context] * F_RANK_MM_LR1 + staticPredictor[context] * F_RANK_MM_LR2) >> 5;

                    if (rank & (1 << bit))
                    {
                        ProbCounter_UpdateBit1(&statePredictor[context],  F_RANK_MS_TH1, F_RANK_MS_AR1);
                        ProbCounter_UpdateBit1(&charPredictor[context],   F_RANK_MC_TH1, F_RANK_MC_AR1);
                        ProbCounter_UpdateBit1(&staticPredictor[context], F_RANK_MP_TH1, F_RANK_MP_AR1);

                        EncodeBit1(&coder,probability); context += context + 1;
                    }
                    else
                    {
                        ProbCounter_UpdateBit0(&statePredictor[context],  F_RANK_MS_TH0, F_RANK_MS_AR0);
                        ProbCounter_UpdateBit0(&charPredictor[context],   F_RANK_MC_TH0, F_RANK_MC_AR0);
                        ProbCounter_UpdateBit0(&staticPredictor[context], F_RANK_MP_TH0, F_RANK_MP_AR0);

                        EncodeBit0(&coder,probability); context += context;
                    }
                }
            }
        }
        else
        {
            rankHistory[currentChar] = bsc_log2_256(rank);

            statePredictor  = & model->Rank.Escape.StateModel[state][0];
            charPredictor   = & model->Rank.Escape.CharModel[currentChar][0];
            staticPredictor = & model->Rank.Escape.StaticModel[0];

            for (int context = 1, bit = maxRank; bit >= 0; --bit)
            {
                int probability = (charPredictor[context] * F_RANK_PM_LR0 + statePredictor[context] * F_RANK_PM_LR1 + staticPredictor[context] * F_RANK_PM_LR2) >> 5;

                if (rank & (1 << bit))
                {
                    ProbCounter_UpdateBit1(&statePredictor[context],  F_RANK_PS_TH1, F_RANK_PS_AR1);
                    ProbCounter_UpdateBit1(&charPredictor[context],   F_RANK_PC_TH1, F_RANK_PC_AR1);
                    ProbCounter_UpdateBit1(&staticPredictor[context], F_RANK_PP_TH1, F_RANK_PP_AR1);

                    EncodeBit1(&coder,probability); context += context + 1;
                }
                else
                {
                    ProbCounter_UpdateBit0(&statePredictor[context],  F_RANK_PS_TH0, F_RANK_PS_AR0);
                    ProbCounter_UpdateBit0(&charPredictor[context],   F_RANK_PC_TH0, F_RANK_PC_AR0);
                    ProbCounter_UpdateBit0(&staticPredictor[context], F_RANK_PP_TH0, F_RANK_PP_AR0);

                    EncodeBit0(&coder,probability); context += context;
                }
            }
        }

        avgRank         =   (avgRank * 124 + rank * 4) >> 7;
        rank            =   rank - 1;
        history         =   runHistory[currentChar];
        state           =   model_run_state(contextRank0, contextRun, rank, history);
        statePredictor  = & model->Run.StateModel[state];
        charPredictor   = & model->Run.CharModel[currentChar];
        staticPredictor = & model->Run.StaticModel;

        if (runSize == 1)
        {
            runHistory[currentChar] = (runHistory[currentChar] + 2) >> 2;

            int probability = ((*charPredictor) * F_RUN_TM_LR0 + (*statePredictor) * F_RUN_TM_LR1 + (*staticPredictor) * F_RUN_TM_LR2) >> 5;

            ProbCounter_UpdateBit0(statePredictor,  F_RUN_TS_TH0, F_RUN_TS_AR0);
            ProbCounter_UpdateBit0(charPredictor,   F_RUN_TC_TH0, F_RUN_TC_AR0);
            ProbCounter_UpdateBit0(staticPredictor, F_RUN_TP_TH0, F_RUN_TP_AR0);

            EncodeBit0(&coder,probability);
        }
        else
        {
            {
                int probability = ((*charPredictor) * F_RUN_TM_LR0 + (*statePredictor) * F_RUN_TM_LR1 + (*staticPredictor) * F_RUN_TM_LR2) >> 5;

                ProbCounter_UpdateBit1(statePredictor,  F_RUN_TS_TH1, F_RUN_TS_AR1);
                ProbCounter_UpdateBit1(charPredictor,   F_RUN_TC_TH1, F_RUN_TC_AR1);
                ProbCounter_UpdateBit1(staticPredictor, F_RUN_TP_TH1, F_RUN_TP_AR1);

                EncodeBit1(&coder,probability);
            }

            int bitRunSize = bsc_log2(runSize); runHistory[currentChar] = (runHistory[currentChar] + 3 * bitRunSize + 3) >> 2;

            statePredictor  = & model->Run.Exponent.StateModel[state][0];
            charPredictor   = & model->Run.Exponent.CharModel[currentChar][0];
            staticPredictor = & model->Run.Exponent.StaticModel[0];

            for (int bit = 1; bit < bitRunSize; ++bit, ++statePredictor, ++charPredictor, ++staticPredictor)
            {
                int probability = ((*charPredictor) * F_RUN_EM_LR0 + (*statePredictor) * F_RUN_EM_LR1 + (*staticPredictor) * F_RUN_EM_LR2) >> 5;

                ProbCounter_UpdateBit1(statePredictor,  F_RUN_ES_TH1, F_RUN_ES_AR1);
                ProbCounter_UpdateBit1(charPredictor,   F_RUN_EC_TH1, F_RUN_EC_AR1);
                ProbCounter_UpdateBit1(staticPredictor, F_RUN_EP_TH1, F_RUN_EP_AR1);

                EncodeBit1(&coder,probability);
            }
            {
                int probability = ((*charPredictor) * F_RUN_EM_LR0 + (*statePredictor) * F_RUN_EM_LR1 + (*staticPredictor) * F_RUN_EM_LR2) >> 5;

                ProbCounter_UpdateBit0(statePredictor,  F_RUN_ES_TH0, F_RUN_ES_AR0);
                ProbCounter_UpdateBit0(charPredictor,   F_RUN_EC_TH0, F_RUN_EC_AR0);
                ProbCounter_UpdateBit0(staticPredictor, F_RUN_EP_TH0, F_RUN_EP_AR0);

                EncodeBit0(&coder,probability);
            }

            statePredictor  = & model->Run.Mantissa[bitRunSize].StateModel[state][0];
            charPredictor   = & model->Run.Mantissa[bitRunSize].CharModel[currentChar][0];
            staticPredictor = & model->Run.Mantissa[bitRunSize].StaticModel[0];

            for (int context = 1, bit = bitRunSize - 1; bit >= 0; --bit)
            {
                int probability = (charPredictor[context] * F_RUN_MM_LR0 + statePredictor[context] * F_RUN_MM_LR1 + staticPredictor[context] * F_RUN_MM_LR2) >> 5;
                if (runSize & (1 << bit))
                {
                    ProbCounter_UpdateBit1(&statePredictor[context],  F_RUN_MS_TH1, F_RUN_MS_AR1);
                    ProbCounter_UpdateBit1(&charPredictor[context],   F_RUN_MC_TH1, F_RUN_MC_AR1);
                    ProbCounter_UpdateBit1(&staticPredictor[context], F_RUN_MP_TH1, F_RUN_MP_AR1);

                    EncodeBit1(&coder,probability); if (bitRunSize <= 5) context += context + 1; else context++;
                }
                else
                {
                    ProbCounter_UpdateBit0(&statePredictor[context],  F_RUN_MS_TH0, F_RUN_MS_AR0);
                    ProbCounter_UpdateBit0(&charPredictor[context],   F_RUN_MC_TH0, F_RUN_MC_AR0);
                    ProbCounter_UpdateBit0(&staticPredictor[context], F_RUN_MP_TH0, F_RUN_MP_AR0);

                    EncodeBit0(&coder,probability); if (bitRunSize <= 5) context += context + 0; else context++;
                }
            }
        }

        contextRank0 = ((contextRank0 << 1) | (rank == 0   ? 1    : 0)) & 0x7;
        contextRank4 = ((contextRank4 << 2) | (rank < 3    ? rank : 3)) & 0xff;
        contextRun   = ((contextRun   << 1) | (runSize < 3 ? 1    : 0)) & 0xf;
    }

    return FinishEncoder(&coder);
}

int bsc_qlfc_adaptive_decode(void *vb, const unsigned char * input, unsigned char * output, QlfcStatisticalModel * model)
{
    RangeCoder coder;

    unsigned char MTFTable[ALPHABET_SIZE];

    bsc_qlfc_init_model(model);

    int contextRank0 = 0;
    int contextRank4 = 0;
    int contextRun   = 0;
    int maxRank      = 7;
    int avgRank      = 0;

    unsigned char rankHistory[ALPHABET_SIZE], runHistory[ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; ++i)
    {
        rankHistory[i] = runHistory[i] = 0;
    }

    InitDecoder(&coder, input);
    int n = (int)DecodeWord(&coder);

    unsigned char usedChar[ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; ++i) usedChar[i] = 0;

    int prevChar = -1;
    for (int rank = 0; rank < ALPHABET_SIZE; ++rank)
    {
        int currentChar = 0;

        for (int bit = 7; bit >= 0; --bit)
        {
            bool bit0 = false, bit1 = false;

            for (int c = 0; c < ALPHABET_SIZE; ++c)
            {
                if (c == prevChar || usedChar[c] == 0)
                {
                    if (currentChar == (c >> (bit + 1)))
                    {
                        if (c & (1 << bit)) bit1 = true; else bit0 = true;
                        if (bit0 && bit1) break;
                    }
                }
            }

            if (bit0 && bit1)
            {
                currentChar += currentChar + DecodeBit (&coder);
            }
            else
            {
                if (bit0) currentChar += currentChar + 0;
                if (bit1) currentChar += currentChar + 1;
            }
        }

        MTFTable[rank] =  currentChar;

        if (currentChar == prevChar)
        {
            maxRank = bsc_log2_256(rank - 1);
            break;
        }

        prevChar = currentChar; usedChar[currentChar] = 1;
    }

    for (int i = 0; i < n;)
    {
        int                 currentChar     =   MTFTable[0];
        int                 history         =   rankHistory[currentChar];
        int                 state           =   model_rank_state(contextRank4, contextRun, history);

        short *            restrict statePredictor  = & model->Rank.StateModel[state];
        short *            restrict charPredictor   = & model->Rank.CharModel[currentChar];
        short *            restrict staticPredictor = & model->Rank.StaticModel;
        ProbabilityMixer * restrict mixer           = & model->mixerOfRank[currentChar];

        int rank = 1;
        if (avgRank < 32)
        {
            if (DecodeBitProb (&coder, Mixup (mixer, *charPredictor, *statePredictor, *staticPredictor)))
            {
                ProbCounter_UpdateBit1(statePredictor,  M_RANK_TS_TH1, M_RANK_TS_AR1);
                ProbCounter_UpdateBit1(charPredictor,   M_RANK_TC_TH1, M_RANK_TC_AR1);
                ProbCounter_UpdateBit1(staticPredictor, M_RANK_TP_TH1, M_RANK_TP_AR1);
                UpdateBit1 (mixer, M_RANK_TM_LR0, M_RANK_TM_LR1, M_RANK_TM_LR2, M_RANK_TM_TH1, M_RANK_TM_AR1);

                statePredictor  = & model->Rank.Exponent.StateModel[state][0];
                charPredictor   = & model->Rank.Exponent.CharModel[currentChar][0];
                staticPredictor = & model->Rank.Exponent.StaticModel[0];
                mixer           = & model->mixerOfRankExponent[history < 1 ? 1 : history][1];

                int bitRankSize = 1;
                while (true)
                {
                    if (bitRankSize == maxRank) break;
                    if (DecodeBitProb (&coder, Mixup (mixer, *charPredictor, *statePredictor, *staticPredictor)))
                    {
                        ProbCounter_UpdateBit1(statePredictor,  M_RANK_ES_TH1, M_RANK_ES_AR1); statePredictor++;
                        ProbCounter_UpdateBit1(charPredictor,   M_RANK_EC_TH1, M_RANK_EC_AR1); charPredictor++;
                        ProbCounter_UpdateBit1(staticPredictor, M_RANK_EP_TH1, M_RANK_EP_AR1); staticPredictor++;
                        UpdateBit1 (mixer, M_RANK_EM_LR0, M_RANK_EM_LR1, M_RANK_EM_LR2, M_RANK_EM_TH1, M_RANK_EM_AR1);
                        bitRankSize++;
                        mixer = & model->mixerOfRankExponent[history < bitRankSize ? bitRankSize : history][bitRankSize];
                    }
                    else
                    {
                        ProbCounter_UpdateBit0(statePredictor,  M_RANK_ES_TH0, M_RANK_ES_AR0);
                        ProbCounter_UpdateBit0(charPredictor,   M_RANK_EC_TH0, M_RANK_EC_AR0);
                        ProbCounter_UpdateBit0(staticPredictor, M_RANK_EP_TH0, M_RANK_EP_AR0);
                        UpdateBit0 (mixer, M_RANK_EM_LR0, M_RANK_EM_LR1, M_RANK_EM_LR2, M_RANK_EM_TH0, M_RANK_EM_AR0);
                        break;
                    }
                }

                rankHistory[currentChar] = bitRankSize;

                statePredictor  = & model->Rank.Mantissa[bitRankSize].StateModel[state][0];
                charPredictor   = & model->Rank.Mantissa[bitRankSize].CharModel[currentChar][0];
                staticPredictor = & model->Rank.Mantissa[bitRankSize].StaticModel[0];
                mixer           = & model->mixerOfRankMantissa[bitRankSize];

                for (int bit = bitRankSize - 1; bit >= 0; --bit)
                {
                    if (DecodeBitProb (&coder, Mixup (mixer, charPredictor[rank], statePredictor[rank], staticPredictor[rank])))
                    {
                        ProbCounter_UpdateBit1(&statePredictor[rank],  M_RANK_MS_TH1, M_RANK_MS_AR1);
                        ProbCounter_UpdateBit1(&charPredictor[rank],   M_RANK_MC_TH1, M_RANK_MC_AR1);
                        ProbCounter_UpdateBit1(&staticPredictor[rank], M_RANK_MP_TH1, M_RANK_MP_AR1);
                        UpdateBit1 (mixer, M_RANK_MM_LR0, M_RANK_MM_LR1, M_RANK_MM_LR2, M_RANK_MM_TH1, M_RANK_MM_AR1);
                        rank += rank + 1;
                    }
                    else
                    {
                        ProbCounter_UpdateBit0(&statePredictor[rank],  M_RANK_MS_TH0, M_RANK_MS_AR0);
                        ProbCounter_UpdateBit0(&charPredictor[rank],   M_RANK_MC_TH0, M_RANK_MC_AR0);
                        ProbCounter_UpdateBit0(&staticPredictor[rank], M_RANK_MP_TH0, M_RANK_MP_AR0);
                        UpdateBit0 (mixer, M_RANK_MM_LR0, M_RANK_MM_LR1, M_RANK_MM_LR2, M_RANK_MM_TH0, M_RANK_MM_AR0);
                        rank += rank;
                    }
                }
            }
            else
            {
                rankHistory[currentChar] = 0;
                ProbCounter_UpdateBit0(statePredictor, M_RANK_TS_TH0,  M_RANK_TS_AR0);
                ProbCounter_UpdateBit0(charPredictor, M_RANK_TC_TH0,   M_RANK_TC_AR0);
                ProbCounter_UpdateBit0(staticPredictor, M_RANK_TP_TH0, M_RANK_TP_AR0);
                UpdateBit0 (mixer, M_RANK_TM_LR0, M_RANK_TM_LR1, M_RANK_TM_LR2, M_RANK_TM_TH0, M_RANK_TM_AR0);
            }
        }
        else
        {
            statePredictor  = & model->Rank.Escape.StateModel[state][0];
            charPredictor   = & model->Rank.Escape.CharModel[currentChar][0];
            staticPredictor = & model->Rank.Escape.StaticModel[0];

            rank = 0;
            for (int context = 1, bit = maxRank; bit >= 0; --bit)
            {
                mixer = & model->mixerOfRankEscape[context];

                if (DecodeBitProb (&coder, Mixup (mixer, charPredictor[context], statePredictor[context], staticPredictor[context])))
                {
                    ProbCounter_UpdateBit1(&statePredictor[context],  M_RANK_PS_TH1, M_RANK_PS_AR1);
                    ProbCounter_UpdateBit1(&charPredictor[context],   M_RANK_PC_TH1, M_RANK_PC_AR1);
                    ProbCounter_UpdateBit1(&staticPredictor[context], M_RANK_PP_TH1, M_RANK_PP_AR1);
                    UpdateBit1 (mixer, M_RANK_PM_LR0, M_RANK_PM_LR1, M_RANK_PM_LR2, M_RANK_PM_TH1, M_RANK_PM_AR1);
                    context += context + 1; rank += rank + 1;
                }
                else
                {
                    ProbCounter_UpdateBit0(&statePredictor[context],  M_RANK_PS_TH0, M_RANK_PS_AR0);
                    ProbCounter_UpdateBit0(&charPredictor[context],   M_RANK_PC_TH0, M_RANK_PC_AR0);
                    ProbCounter_UpdateBit0(&staticPredictor[context], M_RANK_PP_TH0, M_RANK_PP_AR0);
                    UpdateBit0 (mixer, M_RANK_PM_LR0, M_RANK_PM_LR1, M_RANK_PM_LR2, M_RANK_PM_TH0, M_RANK_PM_AR0);
                    context += context; rank += rank;
                }
            }

            rankHistory[currentChar] = bsc_log2_256(rank);
        }

        {
            for (int r = 0; r < rank; ++r)
            {
                MTFTable[r] = MTFTable[r + 1];
            }
            MTFTable[rank] = currentChar;
        }

        avgRank         =   (avgRank * 124 + rank * 4) >> 7;
        rank            =   rank - 1;
        history         =   runHistory[currentChar];
        state           =   model_run_state(contextRank0, contextRun, rank, history);
        statePredictor  = & model->Run.StateModel[state];
        charPredictor   = & model->Run.CharModel[currentChar];
        staticPredictor = & model->Run.StaticModel;
        mixer           = & model->mixerOfRun[currentChar];

        int runSize = 1;
        if (DecodeBitProb (&coder, Mixup (mixer, *charPredictor, *statePredictor, *staticPredictor)))
        {
            ProbCounter_UpdateBit1(statePredictor,  M_RUN_TS_TH1, M_RUN_TS_AR1);
            ProbCounter_UpdateBit1(charPredictor,   M_RUN_TC_TH1, M_RUN_TC_AR1);
            ProbCounter_UpdateBit1(staticPredictor, M_RUN_TP_TH1, M_RUN_TP_AR1);
            UpdateBit1 (mixer, M_RUN_TM_LR0, M_RUN_TM_LR1, M_RUN_TM_LR2, M_RUN_TM_TH1, M_RUN_TM_AR1);

            statePredictor  = & model->Run.Exponent.StateModel[state][0];
            charPredictor   = & model->Run.Exponent.CharModel[currentChar][0];
            staticPredictor = & model->Run.Exponent.StaticModel[0];
            mixer           = & model->mixerOfRunExponent[history < 1 ? 1 : history][1];

            int bitRunSize = 1;
            while (true)
            {
                if (DecodeBitProb (&coder, Mixup (mixer, *charPredictor, *statePredictor, *staticPredictor)))
                {
                    ProbCounter_UpdateBit1(statePredictor,  M_RUN_ES_TH1, M_RUN_ES_AR1); statePredictor++;
                    ProbCounter_UpdateBit1(charPredictor,   M_RUN_EC_TH1, M_RUN_EC_AR1); charPredictor++;
                    ProbCounter_UpdateBit1(staticPredictor, M_RUN_EP_TH1, M_RUN_EP_AR1); staticPredictor++;
                    UpdateBit1 (mixer, M_RUN_EM_LR0, M_RUN_EM_LR1, M_RUN_EM_LR2, M_RUN_EM_TH1, M_RUN_EM_AR1);
                    bitRunSize++; mixer = & model->mixerOfRunExponent[history < bitRunSize ? bitRunSize : history][bitRunSize];
                }
                else
                {
                    ProbCounter_UpdateBit0(statePredictor,  M_RUN_ES_TH0, M_RUN_ES_AR0);
                    ProbCounter_UpdateBit0(charPredictor,   M_RUN_EC_TH0, M_RUN_EC_AR0);
                    ProbCounter_UpdateBit0(staticPredictor, M_RUN_EP_TH0, M_RUN_EP_AR0);
                    UpdateBit0 (mixer, M_RUN_EM_LR0, M_RUN_EM_LR1, M_RUN_EM_LR2, M_RUN_EM_TH0, M_RUN_EM_AR0);
                    break;
                }
            }

            runHistory[currentChar] = (runHistory[currentChar] + 3 * bitRunSize + 3) >> 2;

            statePredictor  = & model->Run.Mantissa[bitRunSize].StateModel[state][0];
            charPredictor   = & model->Run.Mantissa[bitRunSize].CharModel[currentChar][0];
            staticPredictor = & model->Run.Mantissa[bitRunSize].StaticModel[0];
            mixer           = & model->mixerOfRunMantissa[bitRunSize];

            for (int context = 1, bit = bitRunSize - 1; bit >= 0; --bit)
            {
                if (DecodeBitProb (&coder, Mixup (mixer, charPredictor[context], statePredictor[context], staticPredictor[context])))
                {
                    ProbCounter_UpdateBit1(&statePredictor[context],  M_RUN_MS_TH1, M_RUN_MS_AR1);
                    ProbCounter_UpdateBit1(&charPredictor[context],   M_RUN_MC_TH1, M_RUN_MC_AR1);
                    ProbCounter_UpdateBit1(&staticPredictor[context], M_RUN_MP_TH1, M_RUN_MP_AR1);
                    UpdateBit1 (mixer, M_RUN_MM_LR0, M_RUN_MM_LR1, M_RUN_MM_LR2, M_RUN_MM_TH1, M_RUN_MM_AR1);
                    runSize += runSize + 1; if (bitRunSize <= 5) context += context + 1; else context++;
                }
                else
                {
                    ProbCounter_UpdateBit0(&statePredictor[context],  M_RUN_MS_TH0, M_RUN_MS_AR0);
                    ProbCounter_UpdateBit0(&charPredictor[context],   M_RUN_MC_TH0, M_RUN_MC_AR0);
                    ProbCounter_UpdateBit0(&staticPredictor[context], M_RUN_MP_TH0, M_RUN_MP_AR0);
                    UpdateBit0 (mixer, M_RUN_MM_LR0, M_RUN_MM_LR1, M_RUN_MM_LR2, M_RUN_MM_TH0, M_RUN_MM_AR0);
                    runSize += runSize; if (bitRunSize <= 5) context += context; else context++;
                }
            }

        }
        else
        {
            runHistory[currentChar] = (runHistory[currentChar] + 2) >> 2;
            ProbCounter_UpdateBit0(statePredictor,  M_RUN_TS_TH0, M_RUN_TS_AR0);
            ProbCounter_UpdateBit0(charPredictor,   M_RUN_TC_TH0, M_RUN_TC_AR0);
            ProbCounter_UpdateBit0(staticPredictor, M_RUN_TP_TH0, M_RUN_TP_AR0);
            UpdateBit0 (mixer, M_RUN_TM_LR0, M_RUN_TM_LR1, M_RUN_TM_LR2, M_RUN_TM_TH0, M_RUN_TM_AR0);
        }

        contextRank0 = ((contextRank0 << 1) | (rank == 0   ? 1    : 0)) & 0x7;
        contextRank4 = ((contextRank4 << 2) | (rank < 3    ? rank : 3)) & 0xff;
        contextRun   = ((contextRun   << 1) | (runSize < 3 ? 1    : 0)) & 0xf;

        for (; runSize > 0; --runSize) output[i++] = currentChar;
    }

    return n;
}

int bsc_qlfc_static_decode(void *vb, const unsigned char * input, unsigned char * output, QlfcStatisticalModel * model)
{
    RangeCoder coder;

    unsigned char MTFTable[ALPHABET_SIZE];

    bsc_qlfc_init_model(model);

    int contextRank0 = 0;
    int contextRank4 = 0;
    int contextRun   = 0;
    int maxRank      = 7;
    int avgRank      = 0;

    unsigned char rankHistory[ALPHABET_SIZE], runHistory[ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; ++i)
    {
        rankHistory[i] = runHistory[i] = 0;
    }

    InitDecoder(&coder, input);
    int n = (int)DecodeWord(&coder);

    unsigned char usedChar[ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; ++i) usedChar[i] = 0;

    int prevChar = -1;
    for (int rank = 0; rank < ALPHABET_SIZE; ++rank)
    {
        int currentChar = 0;

        for (int bit = 7; bit >= 0; --bit)
        {
            bool bit0 = false, bit1 = false;

            for (int c = 0; c < ALPHABET_SIZE; ++c)
            {
                if (c == prevChar || usedChar[c] == 0)
                {
                    if (currentChar == (c >> (bit + 1)))
                    {
                        if (c & (1 << bit)) bit1 = true; else bit0 = true;
                        if (bit0 && bit1) break;
                    }
                }
            }

            if (bit0 && bit1)
            {
                currentChar += currentChar + DecodeBit(&coder);
            }
            else
            {
                if (bit0) currentChar += currentChar + 0;
                if (bit1) currentChar += currentChar + 1;
            }
        }

        MTFTable[rank] =  currentChar;

        if (currentChar == prevChar)
        {
            maxRank = bsc_log2_256(rank - 1);
            break;
        }

        prevChar = currentChar; usedChar[currentChar] = 1;
    }

    for (int i = 0; i < n;)
    {
        int                 currentChar     =   MTFTable[0];
        int                 history         =   rankHistory[currentChar];
        int                 state           =   model_rank_state(contextRank4, contextRun, history);

        short * restrict    statePredictor  = & model->Rank.StateModel[state];
        short * restrict    charPredictor   = & model->Rank.CharModel[currentChar];
        short * restrict    staticPredictor = & model->Rank.StaticModel;

        int rank = 1;
        if (avgRank < 32)
        {
            if (DecodeBitProb(&coder, (*charPredictor * F_RANK_TM_LR0 + *statePredictor * F_RANK_TM_LR1 + *staticPredictor * F_RANK_TM_LR2) >> 5))
            {
                ProbCounter_UpdateBit1(statePredictor,  F_RANK_TS_TH1, F_RANK_TS_AR1);
                ProbCounter_UpdateBit1(charPredictor,   F_RANK_TC_TH1, F_RANK_TC_AR1);
                ProbCounter_UpdateBit1(staticPredictor, F_RANK_TP_TH1, F_RANK_TP_AR1);

                statePredictor  = & model->Rank.Exponent.StateModel[state][0];
                charPredictor   = & model->Rank.Exponent.CharModel[currentChar][0];
                staticPredictor = & model->Rank.Exponent.StaticModel[0];

                int bitRankSize = 1;
                while (true)
                {
                    if (bitRankSize == maxRank) break;
                    if (DecodeBitProb (&coder, (*charPredictor * F_RANK_EM_LR0 + *statePredictor * F_RANK_EM_LR1 + *staticPredictor * F_RANK_EM_LR2) >> 5))
                    {
                        ProbCounter_UpdateBit1(statePredictor,  F_RANK_ES_TH1, F_RANK_ES_AR1); statePredictor++;
                        ProbCounter_UpdateBit1(charPredictor,   F_RANK_EC_TH1, F_RANK_EC_AR1); charPredictor++;
                        ProbCounter_UpdateBit1(staticPredictor, F_RANK_EP_TH1, F_RANK_EP_AR1); staticPredictor++;
                        bitRankSize++;
                    }
                    else
                    {
                        ProbCounter_UpdateBit0(statePredictor,  F_RANK_ES_TH0, F_RANK_ES_AR0);
                        ProbCounter_UpdateBit0(charPredictor,   F_RANK_EC_TH0, F_RANK_EC_AR0);
                        ProbCounter_UpdateBit0(staticPredictor, F_RANK_EP_TH0, F_RANK_EP_AR0);
                        break;
                    }
                }

                rankHistory[currentChar] = bitRankSize;

                statePredictor  = & model->Rank.Mantissa[bitRankSize].StateModel[state][0];
                charPredictor   = & model->Rank.Mantissa[bitRankSize].CharModel[currentChar][0];
                staticPredictor = & model->Rank.Mantissa[bitRankSize].StaticModel[0];

                for (int bit = bitRankSize - 1; bit >= 0; --bit)
                {
                    if (DecodeBitProb (&coder, (charPredictor[rank] * F_RANK_MM_LR0 + statePredictor[rank] * F_RANK_MM_LR1 + staticPredictor[rank] * F_RANK_MM_LR2) >> 5))
                    {
                        ProbCounter_UpdateBit1(&statePredictor[rank],  F_RANK_MS_TH1, F_RANK_MS_AR1);
                        ProbCounter_UpdateBit1(&charPredictor[rank],   F_RANK_MC_TH1, F_RANK_MC_AR1);
                        ProbCounter_UpdateBit1(&staticPredictor[rank], F_RANK_MP_TH1, F_RANK_MP_AR1);
                        rank += rank + 1;
                    }
                    else
                    {
                        ProbCounter_UpdateBit0(&statePredictor[rank],  F_RANK_MS_TH0, F_RANK_MS_AR0);
                        ProbCounter_UpdateBit0(&charPredictor[rank],   F_RANK_MC_TH0, F_RANK_MC_AR0);
                        ProbCounter_UpdateBit0(&staticPredictor[rank], F_RANK_MP_TH0, F_RANK_MP_AR0);
                        rank += rank;
                    }
                }
            }
            else
            {
                rankHistory[currentChar] = 0;
                ProbCounter_UpdateBit0(statePredictor,  F_RANK_TS_TH0, F_RANK_TS_AR0);
                ProbCounter_UpdateBit0(charPredictor,   F_RANK_TC_TH0, F_RANK_TC_AR0);
                ProbCounter_UpdateBit0(staticPredictor, F_RANK_TP_TH0, F_RANK_TP_AR0);
            }
        }
        else
        {
            statePredictor  = & model->Rank.Escape.StateModel[state][0];
            charPredictor   = & model->Rank.Escape.CharModel[currentChar][0];
            staticPredictor = & model->Rank.Escape.StaticModel[0];

            rank = 0;
            for (int context = 1, bit = maxRank; bit >= 0; --bit)
            {
                if (DecodeBitProb (&coder, (charPredictor[context] * F_RANK_PM_LR0 + statePredictor[context] * F_RANK_PM_LR1 + staticPredictor[context] * F_RANK_PM_LR2) >> 5))
                {
                    ProbCounter_UpdateBit1(&statePredictor[context],  F_RANK_PS_TH1, F_RANK_PS_AR1);
                    ProbCounter_UpdateBit1(&charPredictor[context],   F_RANK_PC_TH1, F_RANK_PC_AR1);
                    ProbCounter_UpdateBit1(&staticPredictor[context], F_RANK_PP_TH1, F_RANK_PP_AR1);
                    context += context + 1; rank += rank + 1;
                }
                else
                {
                    ProbCounter_UpdateBit0(&statePredictor[context],  F_RANK_PS_TH0, F_RANK_PS_AR0);
                    ProbCounter_UpdateBit0(&charPredictor[context],   F_RANK_PC_TH0, F_RANK_PC_AR0);
                    ProbCounter_UpdateBit0(&staticPredictor[context], F_RANK_PP_TH0, F_RANK_PP_AR0);
                    context += context; rank += rank;
                }
            }

            rankHistory[currentChar] = bsc_log2_256(rank);
        }

        {
            for (int r = 0; r < rank; ++r)
            {
                MTFTable[r] = MTFTable[r + 1];
            }
            MTFTable[rank] = currentChar;
        }

        avgRank         =   (avgRank * 124 + rank * 4) >> 7;
        rank            =   rank - 1;
        history         =   runHistory[currentChar];
        state           =   model_run_state(contextRank0, contextRun, rank, history);
        statePredictor  = & model->Run.StateModel[state];
        charPredictor   = & model->Run.CharModel[currentChar];
        staticPredictor = & model->Run.StaticModel;

        int runSize = 1;
        if (DecodeBitProb (&coder, (*charPredictor * F_RUN_TM_LR0 + *statePredictor * F_RUN_TM_LR1 + *staticPredictor * F_RUN_TM_LR2) >> 5))
        {
            ProbCounter_UpdateBit1(statePredictor,  F_RUN_TS_TH1, F_RUN_TS_AR1);
            ProbCounter_UpdateBit1(charPredictor,   F_RUN_TC_TH1, F_RUN_TC_AR1);
            ProbCounter_UpdateBit1(staticPredictor, F_RUN_TP_TH1, F_RUN_TP_AR1);

            statePredictor  = & model->Run.Exponent.StateModel[state][0];
            charPredictor   = & model->Run.Exponent.CharModel[currentChar][0];
            staticPredictor = & model->Run.Exponent.StaticModel[0];

            int bitRunSize = 1;
            while (true)
            {
                if (DecodeBitProb (&coder, (*charPredictor * F_RUN_EM_LR0 + *statePredictor * F_RUN_EM_LR1 + *staticPredictor * F_RUN_EM_LR2) >> 5))
                {
                    ProbCounter_UpdateBit1(statePredictor,  F_RUN_ES_TH1, F_RUN_ES_AR1); statePredictor++;
                    ProbCounter_UpdateBit1(charPredictor,   F_RUN_EC_TH1, F_RUN_EC_AR1); charPredictor++;
                    ProbCounter_UpdateBit1(staticPredictor, F_RUN_EP_TH1, F_RUN_EP_AR1); staticPredictor++;
                    bitRunSize++;
                }
                else
                {
                    ProbCounter_UpdateBit0(statePredictor,  F_RUN_ES_TH0, F_RUN_ES_AR0);
                    ProbCounter_UpdateBit0(charPredictor,   F_RUN_EC_TH0, F_RUN_EC_AR0);
                    ProbCounter_UpdateBit0(staticPredictor, F_RUN_EP_TH0, F_RUN_EP_AR0);
                    break;
                }
            }

            runHistory[currentChar] = (runHistory[currentChar] + 3 * bitRunSize + 3) >> 2;

            statePredictor  = & model->Run.Mantissa[bitRunSize].StateModel[state][0];
            charPredictor   = & model->Run.Mantissa[bitRunSize].CharModel[currentChar][0];
            staticPredictor = & model->Run.Mantissa[bitRunSize].StaticModel[0];

            for (int context = 1, bit = bitRunSize - 1; bit >= 0; --bit)
            {
                if (DecodeBitProb (&coder, (charPredictor[context] * F_RUN_MM_LR0 + statePredictor[context] * F_RUN_MM_LR1 + staticPredictor[context] * F_RUN_MM_LR2) >> 5))
                {
                    ProbCounter_UpdateBit1(&statePredictor[context],  F_RUN_MS_TH1, F_RUN_MS_AR1);
                    ProbCounter_UpdateBit1(&charPredictor[context],   F_RUN_MC_TH1, F_RUN_MC_AR1);
                    ProbCounter_UpdateBit1(&staticPredictor[context], F_RUN_MP_TH1, F_RUN_MP_AR1);
                    runSize += runSize + 1; if (bitRunSize <= 5) context += context + 1; else context++;
                }
                else
                {
                    ProbCounter_UpdateBit0(&statePredictor[context],  F_RUN_MS_TH0, F_RUN_MS_AR0);
                    ProbCounter_UpdateBit0(&charPredictor[context],   F_RUN_MC_TH0, F_RUN_MC_AR0);
                    ProbCounter_UpdateBit0(&staticPredictor[context], F_RUN_MP_TH0, F_RUN_MP_AR0);
                    runSize += runSize; if (bitRunSize <= 5) context += context; else context++;
                }
            }

        }
        else
        {
            runHistory[currentChar] = (runHistory[currentChar] + 2) >> 2;
            ProbCounter_UpdateBit0(statePredictor,  F_RUN_TS_TH0, F_RUN_TS_AR0);
            ProbCounter_UpdateBit0(charPredictor,   F_RUN_TC_TH0, F_RUN_TC_AR0);
            ProbCounter_UpdateBit0(staticPredictor, F_RUN_TP_TH0, F_RUN_TP_AR0);
        }

        contextRank0 = ((contextRank0 << 1) | (rank == 0   ? 1    : 0)) & 0x7;
        contextRank4 = ((contextRank4 << 2) | (rank < 3    ? rank : 3)) & 0xff;
        contextRun   = ((contextRun   << 1) | (runSize < 3 ? 1    : 0)) & 0xf;

        for (; runSize > 0; --runSize) output[i++] = currentChar;
    }

    return n;
}

int bsc_qlfc_static_encode_block(void *vb, const unsigned char * input, unsigned char * output, int inputSize, int outputSize)
{
    QlfcStatisticalModel * model = (QlfcStatisticalModel *)bsc_malloc (vb, sizeof(QlfcStatisticalModel));
    if (model)
    {
        unsigned char * buffer = (unsigned char *)bsc_malloc (vb, inputSize * sizeof(unsigned char));
        if (buffer)
        {
            int result = bsc_qlfc_static_encode (vb, input, output, buffer, inputSize, outputSize, model);

            bsc_free (vb, buffer); bsc_free (vb, model);

            return result;
        };
        bsc_free (vb, model);
    };
    return LIBBSC_NOT_ENOUGH_MEMORY;
}

int bsc_qlfc_adaptive_encode_block (void *vb, const unsigned char * input, unsigned char * output, int inputSize, int outputSize)
{
    QlfcStatisticalModel * model = (QlfcStatisticalModel *)bsc_malloc (vb, sizeof(QlfcStatisticalModel));
    if (model)
    {
        unsigned char * buffer = (unsigned char *)bsc_malloc (vb, inputSize * sizeof(unsigned char));
        if (buffer)
        {
            int result = bsc_qlfc_adaptive_encode (vb, input, output, buffer, inputSize, outputSize, model);

            bsc_free (vb, buffer); bsc_free (vb, model);

            return result;
        };
        bsc_free (vb, model);
    };
    return LIBBSC_NOT_ENOUGH_MEMORY;
}

int bsc_qlfc_static_decode_block(void *vb, const unsigned char * input, unsigned char * output)
{
    QlfcStatisticalModel * model = (QlfcStatisticalModel *)bsc_malloc (vb, sizeof(QlfcStatisticalModel));
    if (model)
    {
        int result = bsc_qlfc_static_decode(vb, input, output, model);

        bsc_free (vb, model);

        return result;
    };
    return LIBBSC_NOT_ENOUGH_MEMORY;
}

int bsc_qlfc_adaptive_decode_block(void *vb, const unsigned char * input, unsigned char * output)
{
    QlfcStatisticalModel * model = (QlfcStatisticalModel *)bsc_malloc (vb, sizeof(QlfcStatisticalModel));
    if (model)
    {
        int result = bsc_qlfc_adaptive_decode(vb, input, output, model);

        bsc_free (vb, model);

        return result;
    };
    return LIBBSC_NOT_ENOUGH_MEMORY;
}
