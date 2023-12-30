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
#include "strings.h"

unsigned vcf_vb_size (DataType dt) { return sizeof (VBlockVCF); }
unsigned vcf_vb_zip_dl_size (void) { return sizeof (ZipDataLineVCF); }
bool vcf_vb_has_haplotype_data (VBlockP vb) { return !!VB_VCF->ht_matrix_ctx; }

// this is also used to determine whether we should reconstruct this VB in LUFT coordinates - the is_translation callback defined in TRANSLATIONS
bool vcf_vb_is_luft    (VBlockP vb) { return vb && VB_VCF->vb_coords == DC_LUFT; }
bool vcf_vb_is_primary (VBlockP vb) { return vb && VB_VCF->vb_coords == DC_PRIMARY; }

int32_t vcf_vb_get_reject_bytes (VBlockP vb) { return VB_VCF->reject_bytes; }

rom vcf_coords_name (int coord)
{
    static rom coords_names[4] = { "NONE", "PRIM", "LUFT", "BOTH" };
    
    return (coord < 0 || coord >= NUM_COORDS) ? "(invalid coord)" : coords_names[coord];
}

// ZIP/PIZ: called before segging / reconstructing each line
void vcf_reset_line (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    vb->sample_i = 0;
    vb->n_alts = 0; // = ALT not parsed yet
    
    CTX(FORMAT_RGQ)->line_has_RGQ = unknown;
    CTX(INFO_DP)->ctx_specific = 0;
    CTX(INFO_QD)->ctx_specific = 0;

    if (IS_ZIP) {
        for (Did did_i=0; did_i < NUM_VCF_FIELDS; did_i++)
            CTX(did_i)->sf_i = -1; // initialize
    }
}
