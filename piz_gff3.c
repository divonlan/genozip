// ------------------------------------------------------------------
//   piz_gff3.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// this is the compute thread entry point. It receives all data of a variant block and processes it
// in memory to the uncompressed format. This thread then terminates the I/O thread writes the output.

#include "vblock.h"
#include "piz.h"
#include "seg.h"
#include "profiler.h"
#include "endianness.h"
#include "zfile.h"
#include "buffer.h"
#include "header.h"
#include "move_to_front.h"
#include "regions.h"
#include "file.h"

static void piz_gff3_reconstruct_array_of_struct (VBlockGFF3 *vb, uint8_t did_i, DictIdType dict_id, 
                                                  unsigned num_items_in_struct, uint32_t txt_line_i)
{
    DECLARE_SNIP;
    LOAD_SNIP (did_i);

    if (snip_len <= 1 || snip[0] != AOS_NUM_ENTRIES) { // not AoS - just copy
        RECONSTRUCT (snip, snip_len);
        return;
    }

    unsigned num_entries = atoi (snip+1);

    MtfContext *ctxs; // an array of length num_items_in_struct (pointer to start of sub-array in vb->mtf_ctx)
    MtfContext *enst_ctx;
    seg_gff3_array_of_struct_ctxs (vb, dict_id, num_items_in_struct, &ctxs, &enst_ctx);

    for (unsigned entry_i=0; entry_i < num_entries; entry_i++) {

        for (unsigned item_i=0; item_i < num_items_in_struct; item_i++) {
            RECONSTRUCT_FROM_DICT (ctxs[item_i].did_i, false);
            RECONSTRUCT1 (" ");
        }

        RECONSTRUCT_ID (enst_ctx->did_i, &vb->enst_data, &vb->next_enst, NULL, false);

        if (entry_i < num_entries-1) RECONSTRUCT1 (",");
    }
}

static bool piz_gff3_reconstruct_special_info_subfields (VBlock *vb_, uint8_t did_i, DictIdType dict_id, uint32_t txt_line_i)
{
    VBlockGFF3 *vb = (VBlockGFF3 *)vb_;

    if (dict_id.num == dict_id_ATTR_ID) {
        DECLARE_SNIP;
        RECONSTRUCT_FROM_DICT_POS (did_i, vb->last_id, true, NULL, false);
        return false;
    }

    if (dict_id.num == dict_id_ATTR_Dbxref) {
        RECONSTRUCT_ID (did_i, &vb->dbxref_numeric_data, &vb->next_dbxref_numeric_data, NULL, false)
        return false;
    }

    if (dict_id.num == dict_id_ATTR_Variant_seq   ||
        dict_id.num == dict_id_ATTR_Reference_seq ||
        dict_id.num == dict_id_ATTR_ancestral_allele) {
        RECONSTRUCT_FROM_BUF (vb->seq_data, vb->next_seq, "SEQ", "", 0);
        return false;
    }

    if (dict_id.num == dict_id_ATTR_Variant_effect      ||
        dict_id.num == dict_id_ATTR_sift_prediction     ||
        dict_id.num == dict_id_ATTR_polyphen_prediction ||
        dict_id.num == dict_id_ATTR_variant_peptide) {
        unsigned num_item_in_struct = (dict_id.num == dict_id_ATTR_variant_peptide ? 2 : 3); // num items excluding ENST Id
        piz_gff3_reconstruct_array_of_struct (vb, did_i, dict_id, num_item_in_struct, txt_line_i);
        return false;
    }

    return true;// proceed with normal reconstruction
}

static void piz_gff3_reconstruct_vb (VBlockGFF3 *vb)
{
    START_TIMER;

    buf_alloc (vb, &vb->txt_data, vb->vb_data_size, 1.1, "txt_data", vb->vblock_i);
    
    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_data_start = vb->txt_data.len;
        uint32_t txt_line_i = vb->first_line + vb_line_i;

        DECLARE_SNIP;
        uint32_t seqid_word_index = RECONSTRUCT_FROM_DICT (GFF3_SEQID, true);
        RECONSTRUCT_FROM_DICT (GFF3_SOURCE, true);
        RECONSTRUCT_FROM_DICT (GFF3_TYPE, true);
        RECONSTRUCT_FROM_DICT_POS (GFF3_START, vb->last_pos, true, NULL, true); // delta vs. previous line START
        RECONSTRUCT_FROM_DICT_POS (GFF3_END, vb->last_pos, false, NULL, true);  // delta vs. START
        RECONSTRUCT_FROM_DICT (GFF3_SCORE, true);
        RECONSTRUCT_FROM_DICT (GFF3_STRAND, true);
        RECONSTRUCT_FROM_DICT (GFF3_PHASE, true);

        uint32_t iname_word_index = LOAD_SNIP (GFF3_ATTRS);        
        bool has_13;
        piz_reconstruct_info ((VBlockP)vb, iname_word_index, snip, snip_len, piz_gff3_reconstruct_special_info_subfields, 
                              txt_line_i, &has_13);
        // add the end-of-line
        RECONSTRUCT (has_13 ? "\r\n" : "\n" , 1 + has_13);

        // after consuming sections' data, if this line is not to be outputed - shorten txt_data back to start of line
        if (flag_regions && !regions_is_site_included (seqid_word_index, vb->last_pos))
            vb->txt_data.len = txt_data_start; // remove excluded line
    }

    COPY_TIMER(vb->profile.piz_reconstruct_vb);
}

void piz_gff3_uncompress_one_vb (VBlock *vb_)
{
    UNCOMPRESS_HEADER_AND_FIELDS (VBlockGFF3, true);
    //UNCOMPRESS_SUBFIELDS (vb->num_info_subfields, SEC_GFF3_ATTRS_SF_B250);
    UNCOMPRESS_DATA_SECTION (SEC_RANDOM_POS_DATA, random_pos_data, uint32_t, true);    
    UNCOMPRESS_DATA_SECTION (SEC_SEQ_DATA, seq_data, char, true);    
    UNCOMPRESS_DATA_SECTION (SEC_NUMERIC_ID_DATA, dbxref_numeric_data, char, true);    
    UNCOMPRESS_DATA_SECTION (SEC_ENST_DATA, enst_data, char, true);    

    piz_gff3_reconstruct_vb (vb);
    UNCOMPRESS_DONE;
}

bool piz_gff3_read_one_vb (VBlock *vb, SectionListEntry *sl)
{ 
    if (vb->vblock_i == 1) piz_map_iname_subfields();
    
    READ_DATA_SECTION (SEC_RANDOM_POS_DATA, true); // POS data that failed delta
    READ_DATA_SECTION (SEC_SEQ_DATA, true); // Data of Variant_seq, Reference_seq and ancestral_allele
    READ_DATA_SECTION (SEC_NUMERIC_ID_DATA, true); // Data of Dbxref
    READ_DATA_SECTION (SEC_ENST_DATA, true); // Data of ENST IDs coming from Variant_effect, sift_prediction, polyphen_prediction, variant_peptide

    return true;
}
