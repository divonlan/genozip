// ------------------------------------------------------------------
//   vcf_vblock.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// vb stands for VBlock - it started its life as VBlockVCF when genozip could only compress VCFs, but now
// it means a block of lines from the text file. 

#include "vcf_private.h"

unsigned vcf_vb_size (DataType dt) { return sizeof (VBlockVCF); }
unsigned vcf_vb_zip_dl_size (void) { return sizeof (ZipDataLineVCF); }

// ZIP/PIZ: called before segging / reconstructing each line
void vcf_reset_line (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    vb->sample_i = 0;
    vb->deferred_q_len = 0;
    vb->mate_line_i = NO_LINE;
    ii_buf.len32 = 0; 

    N_ALTS = 0; // = ALT not parsed yet
    CTX(FORMAT_GT_HT)->use_HT_matrix = false; 
    CTX(FORMAT_RGQ)->line_has_RGQ = unknown;

    Did reset_ctx_specific[] = { 
        FORMAT_SB, INFO_AN, INFO_DP, INFO_QD, INFO_RU, VCF_QUAL,
    };
    for (int i=0; i < ARRAY_LEN(reset_ctx_specific); i++)
        CTX(reset_ctx_specific[i])->ctx_specific = 0;
    
    if (IS_ZIP) {
        memset (B1ST16(CTX(VCF_FORMAT)->sf_i), 0xff, MAX_DICTS * sizeof (uint16_t));  // set all sf_i to NO_SF_I
        memset (&vb->first_idx, 0xff, (char*)&vb->after_idx - (char*)&vb->first_idx); // set all idx's to -1
    }

    else {
        memset (vb->is_deferred_, 0, sizeof (vb->is_deferred_));
    }
}
