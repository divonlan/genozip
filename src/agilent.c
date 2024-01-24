// ------------------------------------------------------------------
//   agilent.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "reconstruct.h"
#include "seg.h"
#include "sam_private.h"

void agilent_seg_initialize (VBlockP vb)
{
    ctx_set_ltype (VB, LT_BLOB, OPTION_QX_Z, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_QX_Z, OPTION_QX_DOMQRUNS, OPTION_QX_QUALMPLX, OPTION_QX_DIVRQUAL, DID_EOL);

    CTX(OPTION_QX_Z)->no_callback = true; // QX is normally compressed with a callback, but not with AGeNT Trimmer
}

//---------------------------------------------------------------------------------------------------
// RX:Z (FASTQ and SAM/BAM): first half + second of molecular barcode, format: "TTA-TTC" (3-hyphen-3)
//---------------------------------------------------------------------------------------------------

void agilent_seg_RX (VBlockP vb_, ContextP ctx, STRp(rx), unsigned add_bytes)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    bool use_special = false;

    if (rx_len != 7 || rx[3] != '-') goto do_seg;

    if (VB_DT(FASTQ)) 
        use_special = ctx_encountered_in_line (VB, OPTION_ZA_Z)   && ctx_encountered_in_line (VB, OPTION_ZB_Z)   &&
                      vb->last_txt_len (OPTION_ZA_Z) >= 3         && vb->last_txt_len (OPTION_ZB_Z) >= 3         &&
                      !memcmp (last_txt (VB, OPTION_ZA_Z), rx, 3) && !memcmp (last_txt (VB, OPTION_ZB_Z), rx+4, 3);
    
    else { // SAM/BAM
        STR(za); STR(zb);
        bool is_bam = IS_BAM_ZIP;
        sam_seg_get_aux_Z (vb, VB_SAM->idx_ZA_Z, pSTRa(za), is_bam);
        sam_seg_get_aux_Z (vb, VB_SAM->idx_ZB_Z, pSTRa(zb), is_bam);

        use_special = has(ZA_Z) && has(ZB_Z) && za_len >= 3 && zb_len >= 3 && !memcmp (rx, za, 3) && !memcmp (rx+4, zb, 3);
    }

do_seg:
    if (use_special)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VB_DT(FASTQ) ? FASTQ_SPECIAL_AGENT_RX : SAM_SPECIAL_AGENT_RX }, 2, ctx, add_bytes);
    else
        seg_by_ctx (VB, STRa(rx), ctx, add_bytes);
}

// e.g. ZA="CCGT" + ZB="GTAGT" --> RX="CCG-GTA" (used for SAM and FASTQ)
SPECIAL_RECONSTRUCTOR (agilent_special_AGENT_RX)
{    
    if (reconstruct) {
        STR(snip); 
        reconstruct_peek (vb, CTX(OPTION_ZA_Z), pSTRa(snip));
        RECONSTRUCT_SEP (snip, 3, '-');

        reconstruct_peek (vb, CTX(OPTION_ZB_Z), pSTRa(snip));
        RECONSTRUCT (snip, 3);
    }

    return NO_NEW_VALUE;
}

//---------------------------------------------------------------------------------------------------
// QX:Z (FASTQ and SAM/BAM): base qualities of RX, format:"DDD DDA"
//---------------------------------------------------------------------------------------------------

void agilent_seg_QX (VBlockP vb, ContextP ctx, STRp(qx), unsigned add_bytes)
{
    if (qx_len == 7 && qx[3] == ' ') {
        seg_add_to_local_fixed (VB, ctx, &qx[0], 3, LOOKUP_NONE, 0);
        seg_add_to_local_fixed (VB, ctx, &qx[4], 3, LOOKUP_NONE, 0);
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VB_DT(FASTQ) ? FASTQ_SPECIAL_AGENT_QX : SAM_SPECIAL_AGENT_QX }, 2, ctx, add_bytes);
    }

    else 
        seg_by_ctx (VB, STRa(qx), ctx, add_bytes);
}

// e.g. "FFF FDF" (used for SAM and FASTQ)
SPECIAL_RECONSTRUCTOR (agilent_special_AGENT_QX)
{
    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, (char[]){ SNIP_LOOKUP, '3'}, 2, reconstruct, __FUNCLINE);
    if (reconstruct) RECONSTRUCT1 (' ');
    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, (char[]){ SNIP_LOOKUP, '3'}, 2, reconstruct, __FUNCLINE);
    
    return NO_NEW_VALUE;
}
