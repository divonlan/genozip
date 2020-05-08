// ------------------------------------------------------------------
//   piz_gff3.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// this is the compute thread entry point. It receives all data of a variant block and processes it
// in memory to the uncompressed format. This thread then terminates the I/O thread writes the output.

#include "vblock.h"
#include "piz.h"
#include "profiler.h"
#include "endianness.h"
#include "zfile.h"
#include "buffer.h"
#include "header.h"
#include "move_to_front.h"
#include "regions.h"
#include "file.h"

static bool piz_gff3_reconstruct_special_info_subfields (VBlock *vb_, uint8_t did_i, DictIdType dict_id, uint32_t txt_line_i)
{
    VBlockGFF3 *vb = (VBlockGFF3 *)vb_;

    if (dict_id.num == dict_id_ATTR_ID) {
        DECLARE_SNIP;
        RECONSTRUCT_FROM_DICT_POS (did_i, vb->last_id, true, NULL, false);
        return false;
    }

    if (dict_id.num == dict_id_ATTR_Dbxref) {
        if (!flag_strip) 
            RECONSTRUCT_ID (did_i, &vb->dbxref_numeric_data, &vb->next_dbxref_numeric_data, NULL, false)
        else 
            RECONSTRUCT1 (".");

        return false;
    }

    if (dict_id.num == dict_id_ATTR_Variant_seq   ||
        dict_id.num == dict_id_ATTR_Reference_seq ||
        dict_id.num == dict_id_ATTR_ancestral_allele) {
        RECONSTRUCT_FROM_BUF (vb->seq_data, vb->next_seq, "SEQ", '\t', "", 0);
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

static void piz_gff3_uncompress_all_sections (VBlockGFF3 *vb)
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
    vb->first_line       = BGEN32 (header->first_line);
    vb->lines.len        = BGEN32 (header->num_lines);
    vb->vb_data_size     = BGEN32 (header->vb_data_size);
    vb->longest_line_len = BGEN32 (header->longest_line_len);

    // in case of --split, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every sam component
    if (flag_split) vb->vblock_i = BGEN32 (header->h.vblock_i);
    
    unsigned section_i=1;

    // uncompress the fields     
    piz_uncompress_fields ((VBlockP)vb, section_index, &section_i);

    // uncompress ATTRIBUTE subfields
    for (uint8_t sf_i=0; sf_i < vb->num_info_subfields ; sf_i++) {
        
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);

        MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, &vb->num_info_subfields, 
                                                  header->dict_id, SEC_GFF3_ATTRS_SF_DICT);

        zfile_uncompress_section ((VBlockP)vb, header, &ctx->b250, "mtf_ctx.b250", SEC_GFF3_ATTRS_SF_B250);    
    }

    // uncompress the seq data
    SectionHeader *seq_header  = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, seq_header, &vb->seq_data, "seq_data", SEC_SEQ_DATA);    

    // uncompress the Dbxref data
    SectionHeader *dbxref_header  = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, dbxref_header, &vb->dbxref_numeric_data, "dbxref_numeric_data", SEC_NUMERIC_ID_DATA);    
}

void piz_gff3_uncompress_one_vb (VBlock *vb_)
{
    START_TIMER;

    VBlockGFF3 *vb = (VBlockGFF3 *)vb_;

    piz_gff3_uncompress_all_sections (vb);

    piz_gff3_reconstruct_vb (vb);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway
    
    COPY_TIMER (vb->profile.compute);
}
