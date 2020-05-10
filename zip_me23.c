// ------------------------------------------------------------------
//   zip_me23.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "vblock.h"
#include "buffer.h"
#include "zfile.h"
#include "zip.h"

void zip_me23_compress_one_vb (VBlockP vb_)
{ 
    // compress ID and Genotype data
    VBlockME23 *vb = (VBlockME23 *)vb_;
    COMPRESS_DATA_SECTION (SEC_NUMERIC_ID_DATA, id_numeric_data, char, COMP_LZMA, false);
    COMPRESS_DATA_SECTION (SEC_HT_DATA, genotype_data, char, COMP_BZ2, false);    
    COMPRESS_DATA_SECTION (SEC_RANDOM_POS_DATA, random_pos_data, uint32_t, COMP_LZMA, true);    
  }
