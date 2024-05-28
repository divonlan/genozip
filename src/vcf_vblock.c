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
    vb->n_alts = 0; // = ALT not parsed yet
    vb->deferred_q_len = 0;
    vb->mate_line_i = NO_LINE;

    CTX(FORMAT_GT_HT)->use_HT_matrix = false; 
    CTX(FORMAT_RGQ)->line_has_RGQ = unknown;
    CTX(FORMAT_SB)->ctx_specific = 0;
    CTX(INFO_BaseCounts)->ctx_specific = 0;
    CTX(INFO_AN)->ctx_specific = 0;
    CTX(INFO_DP)->ctx_specific = 0;
    CTX(INFO_QD)->ctx_specific = 0;
    CTX(INFO_RU)->ctx_specific = 0;
    CTX(VCF_QUAL)->ctx_specific = 0;
    
    if (IS_ZIP) {
        for (Did did_i=0; did_i < NUM_VCF_FIELDS; did_i++)
            CTX(did_i)->sf_i = -1; // initialize

        memset (&vb->first_idx, 0xff, (char*)&vb->after_idx - (char*)&vb->first_idx); // set all idx's to -1

        CTX(VCF_INFO)->info_items.len32 = 0; // info_items 
    }
}
