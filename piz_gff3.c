// ------------------------------------------------------------------
//   piz_gff3.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// this is the compute thread entry point. It receives all data of a variant block and processes it
// in memory to the uncompressed format. This thread then terminates the I/O thread writes the output.

#include "vblock.h"
#include "piz.h"
#include "seg.h"
#include "endianness.h"
#include "buffer.h"
#include "move_to_front.h"
#include "regions.h"
#include "file.h"

static void piz_gff3_reconstruct_array_of_struct (VBlockGFF3 *vb, uint8_t did_i, DictIdType dict_id, 
                                                  unsigned num_items_in_struct)
{
    DECLARE_SNIP;
    snip = AFTERENT (const char, vb->txt_data);
    snip_len = piz_reconstruct_from_ctx ((VBlockP)vb, did_i, "", 0);
    
    if (snip_len <= 1 || snip[0] != SNIP_STRUCTURED) return;  // not AoS - just copy
        
    vb->txt_data.len -= snip_len; // roll back

    unsigned num_entries = atoi (snip+1);

    MtfContext *ctxs[MAX_AoS_ITEMS]; // an array of length num_items_in_struct (pointer to start of sub-array in vb->mtf_ctx)
    MtfContext *enst_ctx;
    seg_gff3_array_of_struct_ctxs (vb, dict_id, num_items_in_struct, ctxs, &enst_ctx);

    for (unsigned entry_i=0; entry_i < num_entries; entry_i++) {

        for (unsigned item_i=0; item_i < num_items_in_struct; item_i++) 
            piz_reconstruct_from_ctx ((VBlockP)vb, ctxs[item_i]->did_i, " ", 1);

        piz_reconstruct_from_ctx ((VBlockP)vb, enst_ctx->did_i, "", 0);

        if (entry_i < num_entries-1) RECONSTRUCT1 (',');
    }
}

static bool piz_gff3_reconstruct_special_info_subfields (VBlock *vb_, uint8_t did_i, DictIdType dict_id)
{
    VBlockGFF3 *vb = (VBlockGFF3 *)vb_;

    if (dict_id.num == dict_id_ATTR_ID) {
        DECLARE_SNIP;
        snip = AFTERENT (const char, vb->txt_data);
        snip_len = piz_reconstruct_from_ctx ((VBlockP)vb, did_i, "", 0);
        RECONSTRUCT1 (SNIP_SEP);
        vb->txt_data.len -= snip_len+1; // roll back

        RECONSTRUCT_FROM_DICT_POS (DID_I_NONE, vb->last_id, true, NULL, false);
        return false;
    }

    if (dict_id.num == dict_id_ATTR_Dbxref) {
        piz_reconstruct_from_ctx ((VBlockP)vb, did_i, "", 0);
        return false;
    }

    if (dict_id.num == dict_id_ATTR_Variant_seq   ||
        dict_id.num == dict_id_ATTR_Reference_seq ||
        dict_id.num == dict_id_ATTR_ancestral_allele) {
        
        // note: all three are stored together in dict_id_ATTR_Variant_seq as they are correlated
        MtfContext *ctx = mtf_get_ctx (vb, (DictIdType)dict_id_ATTR_Variant_seq); 
        piz_reconstruct_from_local_text ((VBlockP)vb, ctx);
        return false;
    }

    if (dict_id.num == dict_id_ATTR_Variant_effect      ||
        dict_id.num == dict_id_ATTR_sift_prediction     ||
        dict_id.num == dict_id_ATTR_polyphen_prediction ||
        dict_id.num == dict_id_ATTR_variant_peptide) {
        unsigned num_item_in_struct = (dict_id.num == dict_id_ATTR_variant_peptide ? 2 : 3); // num items excluding ENST Id
        piz_gff3_reconstruct_array_of_struct (vb, did_i, dict_id, num_item_in_struct);
        return false;
    }

    return true;// proceed with normal reconstruction
}

void piz_gff3_reconstruct_vb (VBlockGFF3 *vb)
{
    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_data_start = vb->txt_data.len;
        vb->line_i = vb->first_line + vb_line_i;

        DECLARE_SNIP;
        uint32_t seqid_word_index = RECONSTRUCT_FROM_DICT (GFF3_SEQID, true);
        piz_reconstruct_from_ctx ((VBlockP)vb, GFF3_SOURCE, "\t", 1);
        piz_reconstruct_from_ctx ((VBlockP)vb, GFF3_TYPE,   "\t", 1);
        RECONSTRUCT_FROM_DICT_POS (GFF3_START, vb->last_pos, true,  NULL, true); // delta vs. previous line START
        RECONSTRUCT_FROM_DICT_POS (GFF3_END,   vb->last_pos, false, NULL, true); // delta vs. START
        piz_reconstruct_from_ctx ((VBlockP)vb, GFF3_SCORE,  "\t", 1);
        piz_reconstruct_from_ctx ((VBlockP)vb, GFF3_STRAND, "\t", 1);
        piz_reconstruct_from_ctx ((VBlockP)vb, GFF3_PHASE,  "\t", 1);

        uint32_t iname_word_index = LOAD_SNIP (GFF3_ATTRS);        
        bool has_13;
        piz_reconstruct_info ((VBlockP)vb, iname_word_index, snip, snip_len, piz_gff3_reconstruct_special_info_subfields, 
                              &has_13);
        // add the end-of-line
        RECONSTRUCT (has_13 ? "\r\n" : "\n" , 1 + has_13);

        // after consuming sections' data, if this line is not to be outputed - shorten txt_data back to start of line
        if (flag_regions && !regions_is_site_included (seqid_word_index, vb->last_pos))
            vb->txt_data.len = txt_data_start; // remove excluded line
    }
}

bool piz_gff3_read_one_vb (VBlock *vb, SectionListEntry *sl)
{ 
    if (vb->vblock_i == 1) piz_map_iname_subfields(vb);
    return true;
}
