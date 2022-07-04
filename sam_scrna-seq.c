// ------------------------------------------------------------------
//   sam_star.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "chrom.h"
#include "codec.h"
#include "profiler.h"
#include "container.h"

// fields used by STARsolo and 10xGenomics cellranger. Some are SAM-standard and some not.

//------------------------------------------------------------------------------------------------------------------------
// CR:Z - "Cellular barcode. The uncorrected sequence bases of the cellular barcode as reported by the sequencing machine"
//------------------------------------------------------------------------------------------------------------------------

void sam_seg_CR_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(cr), unsigned add_bytes)
{
    if (segconf.running && !segconf.CR_CB_seperator) {
        if      (memchr (cr, '_', cr_len)) segconf.CR_CB_seperator = '_'; // as observed in STARsolo-generated files
        else if (memchr (cr, '-', cr_len)) segconf.CR_CB_seperator = '-'; // as recommended by the SAM spec
    
        segconf.n_CR_CB_seps = str_count_char (STRa(cr), segconf.CR_CB_seperator);
    }

    dl->CR = TXTWORD(cr);
    seg_set_last_txt (VB, CTX(OPTION_CR_Z), STRa(cr));

    if (sam_is_depn_vb && vb->sag && IS_SAG_SOLO && str_issame_(STRa(cr), STRBw(z_file->solo_data, vb->solo_aln->CR)))
        sam_seg_against_sa_group (vb, CTX(OPTION_CR_Z), add_bytes);

    else
        sam_seg_buddied_Z_fields (vb, dl, BUDDIED_CR, OPTION_CR_CB, STRa(cr), segconf.CR_CB_seperator, _OPTION_CR_CB, add_bytes);

}

//------------------------------------------------------------------------------------------------------------------
// CB:Z - "Cell identifier, consisting of the optionally-corrected cellular barcode sequence and an optional suffix"
//------------------------------------------------------------------------------------------------------------------

void sam_seg_CB_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(cb), unsigned add_bytes)
{
    dl->CB = TXTWORD(cb);

    if (sam_is_depn_vb && vb->sag && IS_SAG_SOLO && str_issame_(STRa(cb), STRBw(z_file->solo_data, vb->solo_aln->CB)))
        sam_seg_against_sa_group (vb, CTX(OPTION_CB_Z), add_bytes);

    else if (ctx_encountered_in_line (VB, OPTION_CR_Z) &&
        str_issame_(STRa(cb), last_txt(VB, OPTION_CR_Z), vb->last_txt_len(OPTION_CR_Z)))
    
        seg_by_did (VB, STRa(copy_CR_snip), OPTION_CB_Z, add_bytes); // copy from CR

    else
        sam_seg_buddied_Z_fields (vb, dl, BUDDIED_CB, OPTION_CR_CB, STRa(cb), segconf.CR_CB_seperator, _OPTION_CR_CB, add_bytes);
}

//-----------------------------------------------------------------------------------------------------
// CY:Z - "Phred quality of the cellular barcode sequence in the CR tag"
//-----------------------------------------------------------------------------------------------------

void sam_seg_CY_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(cy), unsigned add_bytes)
{
    dl->CY = TXTWORD(cy);

    if (segconf.running && !segconf.CY_con_snip_len) {
        char sep = 0;
        if      (str_count_char (STRa(cy), ' ') == segconf.n_CR_CB_seps) sep = ' '; // as recommended by the SAM spec (not valid QUAL phred score)
        // note: '_' might be a valid phred score - hence the additional check of verifying n_CR_CB_seps
        else if (str_count_char (STRa(cy), '_') == segconf.n_CR_CB_seps) sep = '_'; // as observed in STARsolo-generated files
    
        if (sep) {
            str_split (cy, cy_len, 0, sep, item, false);
            
            segconf.CY_con = (MiniContainer){
                .nitems_lo         = 1,
                .repeats           = n_items,
                .repsep            = { sep },
                .drop_final_repsep = true,
                .items[0]          = { .dict_id = (DictId)_OPTION_CY_ARR }
            };

            segconf.CY_con_snip_len = sizeof (segconf.CY_con_snip);
            container_prepare_snip ((ContainerP)&segconf.CY_con, 0, 0, segconf.CY_con_snip, &segconf.CY_con_snip_len);
        }
    }

    if (sam_is_depn_vb && vb->sag && IS_SAG_SOLO && str_issame_(STRa(cy), STRBw(z_file->solo_data, vb->solo_aln->CY)))
        sam_seg_against_sa_group (vb, CTX(OPTION_CY_Z), add_bytes);

    else if (segconf.CY_con_snip_len) {

        str_split (cy, cy_len, segconf.CY_con.repeats, segconf.CY_con.repsep[0], item, true);
        if (!n_items) goto fallback;

        for (int i=0; i < n_items; i++) {
            seg_add_to_local_fixed (VB, CTX(OPTION_CY_ARR), STRi(item, i));
            
            SNIPi1 (SNIP_LOOKUP, item_lens[i]);
            seg_by_did (VB, STRa(snip), OPTION_CY_ARR, item_lens[i]);
        }

        seg_by_did (VB, STRa(segconf.CY_con_snip), OPTION_CY_Z, (add_bytes - cy_len) + (n_items-1)); // account for separators
    }

    else fallback: 
        seg_by_did (VB, STRa(cy), OPTION_CY_Z, add_bytes);
}

//------------------------------------------------------------------------------------------------------------------------
// UR:Z - "Chromium molecular barcode sequence as reported by the sequencer"
// UB:Z - "Chromium molecular barcode sequence that is error-corrected among other molecular barcodes with the same cellular barcode and gene alignment"
//------------------------------------------------------------------------------------------------------------------------

void sam_seg_UR_UB_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(value), Did source_did, unsigned add_bytes)
{
    segconf_set_has (source_did);
    if (source_did == OPTION_UR_Z) dl->UR = TXTWORD(value);
    else                           dl->UB = TXTWORD(value);

    // note: UR and UB are segged together - UB is an alias of UR

    if (sam_is_depn_vb && vb->sag && IS_SAG_SOLO && 
          ((source_did == OPTION_UR_Z && str_issame_(STRa(value), STRBw(z_file->solo_data, vb->solo_aln->UR))) ||
           (source_did == OPTION_UB_Z && str_issame_(STRa(value), STRBw(z_file->solo_data, vb->solo_aln->UB)))))
        
        sam_seg_against_sa_group (vb, CTX(OPTION_UR_Z), add_bytes);

    else
        sam_seg_buddied_Z_fields (vb, dl, BUDDIED_UR, OPTION_UR_Z, STRa(value), 0, DICT_ID_NONE, add_bytes);

    seg_set_last_txt (VB, CTX(source_did), STRa(value));
}

//-----------------------------------------------------------------------------------------------------
// UY:Z - "Chromium molecular barcode read quality. Phred scores as reported by sequencer"
//-----------------------------------------------------------------------------------------------------

void sam_seg_UY_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(uy), unsigned add_bytes)
{
    segconf_set_has (OPTION_UY_Z);
    dl->UY = TXTWORD(uy);

    if (sam_is_depn_vb && vb->sag && IS_SAG_SOLO && str_issame_(STRa(uy), STRBw(z_file->solo_data, vb->solo_aln->UY))) {
        sam_seg_against_sa_group (vb, CTX(OPTION_UY_Z), add_bytes);
        dl->dont_compress_UY = true;
    }
    
    else {
        CTX(OPTION_UY_Z)->local.len32 += uy_len;

        SNIPi1 (SNIP_LOOKUP, uy_len);
        seg_by_did (VB, STRa(snip), OPTION_UY_Z, add_bytes);
    }
}
