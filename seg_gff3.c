// ------------------------------------------------------------------
//   seg_gff3.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "move_to_front.h"
#include "header.h"
#include "random_access.h"
#include "file.h"
#include "strings.h"
#include "dict_id.h"

#define DATA_LINE(i) ENT (ZipDataLineGFF3, vb->lines, i)

// called from seg_all_data_lines
void seg_gff3_initialize (VBlock *vb_)
{
    VBlockGFF3 *vb = (VBlockGFF3 *)vb_;

    buf_alloc (vb, &vb->dbxref_numeric_data, sizeof(uint32_t) * vb->lines.len, 1, "dbxref_numeric_data", vb->vblock_i);    
    buf_alloc (vb, &vb->seq_data, 20 * vb->lines.len, 1, "seq_data", vb->vblock_i); // should normally be more than enough, but if not, seg_add_to_data_buf will realloc
}

static bool seg_gff3_special_info_subfields(VBlockP vb_, MtfContextP ctx, const char **this_value, unsigned *this_value_len, char *optimized_snip)
{
    VBlockGFF3 *vb = (VBlockGFF3 *)vb_;

    // ID - this is a sequential number (at least in GRCh37/38)
    if (ctx->dict_id.num == dict_id_ATTR_ID) {
        vb->last_id = seg_pos_field ((VBlockP)vb, vb->last_id, NULL, true, ctx->did_i, SEC_GFF3_ATTRS_SF_B250, 
                                     *this_value, *this_value_len, "ID");
        vb->txt_section_bytes[SEC_GFF3_ATTRS_SF_B250]--; // exclude the separator included by default by seg_pos_field

        return false; // do not add to dictionary/b250 - we already did it
    }

    // Dbxref (example: "dbSNP_151:rs1307114892") - we divide to the non-numeric part which we store
    // in a dictionary and the numeric part which store in a NUMERICAL_ID_DATA section
    if (ctx->dict_id.num == dict_id_ATTR_Dbxref) {
        seg_id_field (vb_, &vb->dbxref_numeric_data, ctx->dict_id, SEC_GFF3_ATTRS_SF_B250, 
                      *this_value, *this_value_len, false, false); // discard the const as seg_id_field modifies

        return false; // do not add to dictionary/b250 - we already did it
    }

    // we store these 3 in one dictionary, as they are correlated and will compress better together
    if (ctx->dict_id.num == dict_id_ATTR_Variant_seq   ||
        ctx->dict_id.num == dict_id_ATTR_Reference_seq ||
        ctx->dict_id.num == dict_id_ATTR_ancestral_allele) {

        seg_add_to_data_buf (vb_, &vb->seq_data, SEC_SEQ_DATA, *this_value, *this_value_len, '\t', *this_value_len);
        return false; // do not add to dictionary/b250 - we already did it
    }

    return true; // all other cases -  procedue with adding to dictionary/b250
}

const char *seg_gff3_data_line (VBlock *vb_,   
                                const char *field_start_line)     // index in vb->txt_data where this line starts
{
    VBlockGFF3 *vb = (VBlockGFF3 *)vb_;
    ZipDataLineGFF3 *dl = DATA_LINE (vb->line_i);

    const char *next_field, *field_start;
    unsigned field_len=0;
    char separator;
    bool has_13 = false; // does this line end in Windows-style \r\n rather than Unix-style \n

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    // SEQID
    field_start = field_start_line;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, &has_13, "SEQID");
    seg_chrom_field (vb_, field_start, field_len);

    // SOURCE
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, &has_13, "SOURCE");
    seg_one_field (vb, field_start, field_len, GFF3_SOURCE);

    // TYPE
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, &has_13, "TYPE");
    seg_one_field (vb, field_start, field_len, GFF3_TYPE);

    // START - delta vs previous line 
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, &has_13, "START");
    vb->last_pos = seg_pos_field (vb_, vb->last_pos, NULL, false, GFF3_START, SEC_GFF3_START_B250, field_start, field_len, "START");
    random_access_update_pos (vb_, vb->last_pos);

    // END - delta vs START
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, &has_13, "END");
    seg_pos_field (vb_, vb->last_pos, NULL, false, GFF3_END, SEC_GFF3_END_B250, field_start, field_len, "END");

    // SCORE
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, &has_13, "SCORE");
    seg_one_field (vb, field_start, field_len, GFF3_SCORE);

    // STRAND
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, &has_13, "STRAND");
    seg_one_field (vb, field_start, field_len, GFF3_STRAND);

    // PHASE
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, &has_13, "PHASE");
    seg_one_field (vb, field_start, field_len, GFF3_PHASE);

    // ATTRIBUTES
    field_start = next_field; 
    next_field = seg_get_next_item (vb, field_start, &len, true, false, false, &field_len, &separator, &has_13, 
                                    field_names[DT_GFF3][GFF3_ATTRS] /* pointer to string to allow pointer comparison */); 

    seg_info_field (vb_, &dl->attrs_mtf_i, &vb->iname_mapper_buf, &vb->num_info_subfields, seg_gff3_special_info_subfields,
                    (char*)field_start, field_len, has_13); // we break the const bc seg_vcf_info_field might add a :#

    vb->txt_section_bytes[SEC_STATS_HT_SEPERATOR] -= has_13; // the \r in case of Windows \r\n line ending (WHY IS THIS?)

    return next_field;
}
