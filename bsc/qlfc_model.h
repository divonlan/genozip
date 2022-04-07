/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Statistical data compression model for QLFC               */
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

#pragma once

#include "predictor.h"

static const int M_RANK_TS_TH0 =    1; static const int M_RANK_TS_AR0 =   57;
static const int M_RANK_TS_TH1 = -111; static const int M_RANK_TS_AR1 =   31;
static const int M_RANK_TC_TH0 =  291; static const int M_RANK_TC_AR0 =  250;
static const int M_RANK_TC_TH1 =  154; static const int M_RANK_TC_AR1 =  528;
static const int M_RANK_TP_TH0 =  375; static const int M_RANK_TP_AR0 =  163;
static const int M_RANK_TP_TH1 =  313; static const int M_RANK_TP_AR1 =  639;
static const int M_RANK_TM_TH0 =  -41; static const int M_RANK_TM_AR0 =   96;
static const int M_RANK_TM_TH1 =   53; static const int M_RANK_TM_AR1 =   49;
static const int M_RANK_TM_LR0 =   20; static const int M_RANK_TM_LR1 =   47;
static const int M_RANK_TM_LR2 =   27;

static const int M_RANK_ES_TH0 = -137; static const int M_RANK_ES_AR0 =   17;
static const int M_RANK_ES_TH1 =  482; static const int M_RANK_ES_AR1 =   40;
static const int M_RANK_EC_TH0 =   61; static const int M_RANK_EC_AR0 =  192;
static const int M_RANK_EC_TH1 =  200; static const int M_RANK_EC_AR1 =  133;
static const int M_RANK_EP_TH0 =   54; static const int M_RANK_EP_AR0 = 1342;
static const int M_RANK_EP_TH1 =  578; static const int M_RANK_EP_AR1 = 1067;
static const int M_RANK_EM_TH0 =  -11; static const int M_RANK_EM_AR0 =  318;
static const int M_RANK_EM_TH1 =  144; static const int M_RANK_EM_AR1 =  848;
static const int M_RANK_EM_LR0 =   49; static const int M_RANK_EM_LR1 =   41;
static const int M_RANK_EM_LR2 =   40;

static const int M_RANK_MS_TH0 = -145; static const int M_RANK_MS_AR0 =   18;
static const int M_RANK_MS_TH1 =  114; static const int M_RANK_MS_AR1 =   24;
static const int M_RANK_MC_TH0 =  -43; static const int M_RANK_MC_AR0 =   69;
static const int M_RANK_MC_TH1 =  -36; static const int M_RANK_MC_AR1 =   78;
static const int M_RANK_MP_TH0 =   -2; static const int M_RANK_MP_AR0 = 1119;
static const int M_RANK_MP_TH1 =   11; static const int M_RANK_MP_AR1 = 1181;
static const int M_RANK_MM_TH0 = -203; static const int M_RANK_MM_AR0 =   20;
static const int M_RANK_MM_TH1 = -271; static const int M_RANK_MM_AR1 =   15;
static const int M_RANK_MM_LR0 =  263; static const int M_RANK_MM_LR1 =  175;
static const int M_RANK_MM_LR2 =   17;

static const int M_RANK_PS_TH0 =  -99; static const int M_RANK_PS_AR0 =   32;
static const int M_RANK_PS_TH1 =  318; static const int M_RANK_PS_AR1 =   42;
static const int M_RANK_PC_TH0 =   17; static const int M_RANK_PC_AR0 =  101;
static const int M_RANK_PC_TH1 = 1116; static const int M_RANK_PC_AR1 =  246;
static const int M_RANK_PP_TH0 =   22; static const int M_RANK_PP_AR0 =  964;
static const int M_RANK_PP_TH1 =   -2; static const int M_RANK_PP_AR1 = 1110;
static const int M_RANK_PM_TH0 = -194; static const int M_RANK_PM_AR0 =   21;
static const int M_RANK_PM_TH1 = -129; static const int M_RANK_PM_AR1 =   20;
static const int M_RANK_PM_LR0 =  480; static const int M_RANK_PM_LR1 =  202;
static const int M_RANK_PM_LR2 =   17;

static const int M_RUN_TS_TH0 =  -93; static const int M_RUN_TS_AR0 =   34;
static const int M_RUN_TS_TH1 =   -4; static const int M_RUN_TS_AR1 =   51;
static const int M_RUN_TC_TH0 =  139; static const int M_RUN_TC_AR0 =  423;
static const int M_RUN_TC_TH1 =  244; static const int M_RUN_TC_AR1 =  162;
static const int M_RUN_TP_TH0 =  275; static const int M_RUN_TP_AR0 =  450;
static const int M_RUN_TP_TH1 =   -6; static const int M_RUN_TP_AR1 =  579;
static const int M_RUN_TM_TH0 =  -68; static const int M_RUN_TM_AR0 =   25;
static const int M_RUN_TM_TH1 =    1; static const int M_RUN_TM_AR1 =   64;
static const int M_RUN_TM_LR0 =   15; static const int M_RUN_TM_LR1 =   50;
static const int M_RUN_TM_LR2 =   78;

static const int M_RUN_ES_TH0 = -116; static const int M_RUN_ES_AR0 =   31;
static const int M_RUN_ES_TH1 =   43; static const int M_RUN_ES_AR1 =   45;
static const int M_RUN_EC_TH0 =  165; static const int M_RUN_EC_AR0 =  222;
static const int M_RUN_EC_TH1 =   30; static const int M_RUN_EC_AR1 =  324;
static const int M_RUN_EP_TH0 =  315; static const int M_RUN_EP_AR0 =  857;
static const int M_RUN_EP_TH1 =  109; static const int M_RUN_EP_AR1 =  867;
static const int M_RUN_EM_TH0 =  -14; static const int M_RUN_EM_AR0 =  215;
static const int M_RUN_EM_TH1 =   61; static const int M_RUN_EM_AR1 =   73;
static const int M_RUN_EM_LR0 =   35; static const int M_RUN_EM_LR1 =   37;
static const int M_RUN_EM_LR2 =   42;

static const int M_RUN_MS_TH0 = -176; static const int M_RUN_MS_AR0 =   14;
static const int M_RUN_MS_TH1 = -141; static const int M_RUN_MS_AR1 =   21;
static const int M_RUN_MC_TH0 =   84; static const int M_RUN_MC_AR0 =  172;
static const int M_RUN_MC_TH1 =   37; static const int M_RUN_MC_AR1 =  263;
static const int M_RUN_MP_TH0 =    2; static const int M_RUN_MP_AR0 =   15;
static const int M_RUN_MP_TH1 = -197; static const int M_RUN_MP_AR1 =   20;
static const int M_RUN_MM_TH0 =  -27; static const int M_RUN_MM_AR0 =  142;
static const int M_RUN_MM_TH1 = -146; static const int M_RUN_MM_AR1 =   27;
static const int M_RUN_MM_LR0 =   51; static const int M_RUN_MM_LR1 =   44;
static const int M_RUN_MM_LR2 =   80;

static const int F_RANK_TS_TH0 = -116; static const int F_RANK_TS_AR0 =   33;
static const int F_RANK_TS_TH1 =  -78; static const int F_RANK_TS_AR1 =   34;
static const int F_RANK_TC_TH0 =   -2; static const int F_RANK_TC_AR0 =  282;
static const int F_RANK_TC_TH1 =   12; static const int F_RANK_TC_AR1 =  274;
static const int F_RANK_TP_TH0 =    4; static const int F_RANK_TP_AR0 =  697;
static const int F_RANK_TP_TH1 =   55; static const int F_RANK_TP_AR1 = 1185;
static const int F_RANK_TM_LR0 =   17; static const int F_RANK_TM_LR1 =   14;
static const int F_RANK_TM_LR2 =    1;

static const int F_RANK_ES_TH0 = -177; static const int F_RANK_ES_AR0 =   23;
static const int F_RANK_ES_TH1 = -370; static const int F_RANK_ES_AR1 =   11;
static const int F_RANK_EC_TH0 =  -14; static const int F_RANK_EC_AR0 =  271;
static const int F_RANK_EC_TH1 =    3; static const int F_RANK_EC_AR1 =  308;
static const int F_RANK_EP_TH0 =   -3; static const int F_RANK_EP_AR0 =  788;
static const int F_RANK_EP_TH1 =  135; static const int F_RANK_EP_AR1 = 1364;
static const int F_RANK_EM_LR0 =   22; static const int F_RANK_EM_LR1 =    6;
static const int F_RANK_EM_LR2 =    4;

static const int F_RANK_MS_TH0 = -254; static const int F_RANK_MS_AR0 =   16;
static const int F_RANK_MS_TH1 = -177; static const int F_RANK_MS_AR1 =   20;
static const int F_RANK_MC_TH0 =  -55; static const int F_RANK_MC_AR0 =   73;
static const int F_RANK_MC_TH1 =  -54; static const int F_RANK_MC_AR1 =   74;
static const int F_RANK_MP_TH0 =   -6; static const int F_RANK_MP_AR0 =  575;
static const int F_RANK_MP_TH1 = 1670; static const int F_RANK_MP_AR1 = 1173;
static const int F_RANK_MM_LR0 =   15; static const int F_RANK_MM_LR1 =   10;
static const int F_RANK_MM_LR2 =    7;

static const int F_RANK_PS_TH0 = -126; static const int F_RANK_PS_AR0 =   32;
static const int F_RANK_PS_TH1 = -126; static const int F_RANK_PS_AR1 =   32;
static const int F_RANK_PC_TH0 =  -33; static const int F_RANK_PC_AR0 =  120;
static const int F_RANK_PC_TH1 =  -25; static const int F_RANK_PC_AR1 =  157;
static const int F_RANK_PP_TH0 =   -6; static const int F_RANK_PP_AR0 =  585;
static const int F_RANK_PP_TH1 =  150; static const int F_RANK_PP_AR1 =  275;
static const int F_RANK_PM_LR0 =   16; static const int F_RANK_PM_LR1 =   11;
static const int F_RANK_PM_LR2 =    5;

static const int F_RUN_TS_TH0 =  -68; static const int F_RUN_TS_AR0 =   38;
static const int F_RUN_TS_TH1 = -112; static const int F_RUN_TS_AR1 =   36;
static const int F_RUN_TC_TH0 =   -4; static const int F_RUN_TC_AR0 =  221;
static const int F_RUN_TC_TH1 =  -13; static const int F_RUN_TC_AR1 =  231;
static const int F_RUN_TP_TH0 =    0; static const int F_RUN_TP_AR0 =    0;
static const int F_RUN_TP_TH1 =    0; static const int F_RUN_TP_AR1 =    0;
static const int F_RUN_TM_LR0 =   14; static const int F_RUN_TM_LR1 =   18;
static const int F_RUN_TM_LR2 =    0;

static const int F_RUN_ES_TH0 =  -90; static const int F_RUN_ES_AR0 =   45;
static const int F_RUN_ES_TH1 =  -92; static const int F_RUN_ES_AR1 =   44;
static const int F_RUN_EC_TH0 =   -3; static const int F_RUN_EC_AR0 =  325;
static const int F_RUN_EC_TH1 =  -11; static const int F_RUN_EC_AR1 =  341;
static const int F_RUN_EP_TH0 =   24; static const int F_RUN_EP_AR0 =  887;
static const int F_RUN_EP_TH1 =   -4; static const int F_RUN_EP_AR1 =  765;
static const int F_RUN_EM_LR0 =   14; static const int F_RUN_EM_LR1 =   15;
static const int F_RUN_EM_LR2 =    3;

static const int F_RUN_MS_TH0 = -275; static const int F_RUN_MS_AR0 =   14;
static const int F_RUN_MS_TH1 = -185; static const int F_RUN_MS_AR1 =   22;
static const int F_RUN_MC_TH0 =  -18; static const int F_RUN_MC_AR0 =  191;
static const int F_RUN_MC_TH1 =  -15; static const int F_RUN_MC_AR1 =  241;
static const int F_RUN_MP_TH0 =  -73; static const int F_RUN_MP_AR0 =   54;
static const int F_RUN_MP_TH1 = -214; static const int F_RUN_MP_AR1 =   19;
static const int F_RUN_MM_LR0 =    7; static const int F_RUN_MM_LR1 =   15;
static const int F_RUN_MM_LR2 =   10;

typedef struct {
    ProbabilityMixer mixerOfRank[ALPHABET_SIZE];
    ProbabilityMixer mixerOfRankExponent[8][8];
    ProbabilityMixer mixerOfRankMantissa[8];
    ProbabilityMixer mixerOfRankEscape[ALPHABET_SIZE];
    ProbabilityMixer mixerOfRun[ALPHABET_SIZE];
    ProbabilityMixer mixerOfRunExponent[32][32];
    ProbabilityMixer mixerOfRunMantissa[32];

    struct 
    {
        short StaticModel;
        short StateModel[ALPHABET_SIZE];
        short CharModel[ALPHABET_SIZE];

        struct
        {
            short StaticModel[8];
            short StateModel[ALPHABET_SIZE][8];
            short CharModel[ALPHABET_SIZE][8];
        } Exponent;

        struct
        {
            short StaticModel[ALPHABET_SIZE];
            short StateModel[ALPHABET_SIZE][ALPHABET_SIZE];
            short CharModel[ALPHABET_SIZE][ALPHABET_SIZE];
        } Mantissa[8];

        struct
        {
            short StaticModel[ALPHABET_SIZE];
            short StateModel[ALPHABET_SIZE][ALPHABET_SIZE];
            short CharModel[ALPHABET_SIZE][ALPHABET_SIZE];
        } Escape;

    } Rank;

    struct 
    {
        short StaticModel;
        short StateModel[ALPHABET_SIZE];
        short CharModel[ALPHABET_SIZE];

        struct
        {
            short StaticModel[32];
            short StateModel[ALPHABET_SIZE][32];
            short CharModel[ALPHABET_SIZE][32];
        } Exponent;

        struct
        {
            short StaticModel[32];
            short StateModel[ALPHABET_SIZE][32];
            short CharModel[ALPHABET_SIZE][32];
        } Mantissa[32];

    } Run;
} QlfcStatisticalModel;

int  bsc_qlfc_init_static_model();
void bsc_qlfc_init_model(QlfcStatisticalModel * model);


/*-----------------------------------------------------------*/
/* End                                          qlfc_model.h */
/*-----------------------------------------------------------*/
