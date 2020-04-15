// ------------------------------------------------------------------
//   piz_sam.c
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
#include "seg.h"
#include "header.h"

static void piz_sam_map_optional_subfields (VBlockSAM *vb)
{
    // terminology: we call a list of INFO subfield names, an "oname". An oname looks something like
    // this: "MX:Z:ab:i:". Each oname consists of OPTIONAL subfields. These subfields are not unique to this
    // oname and can appear in other onames. Each optional subfield is made of a segment of the template eg "MX:Z:"
    // and a value eg "abcded" which is stored in the b250 of that optional subfield.

    // note: an oname might end with a "#" indicating Windows-style newline. We ignore it for now.

    const MtfContext *optional_ctx = &vb->mtf_ctx[SAM_OPTIONAL];
    vb->optional_mapper_buf.len = optional_ctx->word_list.len;
    buf_alloc (vb, &vb->optional_mapper_buf, sizeof (SubfieldMapper) * vb->optional_mapper_buf.len,
               1, "optional_mapper_buf", 0);
    buf_zero (&vb->optional_mapper_buf);

    ARRAY (SubfieldMapper, all_optional_mappers, vb->optional_mapper_buf);

    ARRAY (const MtfWord, all_optionals, optional_ctx->word_list);

    for (unsigned optional_i=0; optional_i < vb->optional_mapper_buf.len; optional_i++) {

        char *optional = ENT (char, optional_ctx->dict, all_optionals[optional_i].char_index); // e.g. "MX:Z:ab:i:" - pointer into the OPTIONAL dictionary
        SubfieldMapper *optional_mapper = &all_optional_mappers[optional_i]; // optional_mapper of this specific set of names "I1=I2=I3="
        optional_mapper->num_subfields = all_optionals[optional_i].snip_len / 5; 

        // traverse the subfields of one optional. E.g. if the optional is "MX:Z:ab:i:" then we traverse MX, ab
        // note - optional 
        for (unsigned i=0; i < optional_mapper->num_subfields; i++) {
                        
            DictIdType dict_id = dict_id_sam_optnl_sf (dict_id_make (&optional[i*5], 2)); // get dict_id from first 2 characters eg "MX"
            optional_mapper->did_i[i] = mtf_get_existing_did_i_by_dict_id ((VBlockP)vb, dict_id); 
        }
    }
}

static void piz_sam_reconstruct_qname (VBlockSAM *vb, const char *template, unsigned template_len, uint32_t sam_line_i)
{
    for (unsigned i=0; i <= template_len; i++) {
        
        const char *snip;
        unsigned snip_len;
        mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[vb->qname_mapper.did_i[i]], NULL, &snip, &snip_len, sam_line_i);

        buf_add (&vb->txt_data, snip, snip_len);

        if (i < template_len) buf_add (&vb->txt_data, &template[i], 1); // add separator
    }
}

static void piz_sam_reconstruct_optional_fields (VBlockSAM *vb, const char *oname, unsigned oname_len, uint32_t opt_word_index, 
                                                 uint32_t sam_line_i)
{
    SubfieldMapper *opt_map = ENT (SubfieldMapper, vb->optional_mapper_buf, opt_word_index);

    ASSERT (opt_map->num_subfields == oname_len / 5, "Error: opt_map->num_subfields=%u but oname=%*.s indicates %u optional fields. sam_line=%u", 
            opt_map->num_subfields, oname_len, oname, oname_len/5, sam_line_i);

    for (unsigned i=0; i < opt_map->num_subfields; i++) {

        const char *snip;
        unsigned snip_len;
        mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[opt_map->did_i[i]], NULL, &snip, &snip_len, sam_line_i);

        buf_add (&vb->txt_data, &oname[i*5], 5)
        buf_add (&vb->txt_data, snip, snip_len);

        if (i != opt_map->num_subfields-1)
            buf_add (&vb->txt_data, "\t", 1);
    }
    
    if (oname_len % 5 == 1 && oname[oname_len-1] == '#') // Windows style end-of-line \r\n
        buf_add (&vb->txt_data, "\r", 1);
}

static void piz_sam_reconstruct_seq_qual (VBlockSAM *vb, uint32_t cigar_seq_len, 
                                          const Buffer *data, uint32_t *next,
                                          const char *field_name, uint32_t sam_line_i)
{
    // seq and qual are expected to be either of length cigar_seq_len, or "*" 
    uint32_t len = (*next >= data->len || data->data[*next] == '*') ? 1 : cigar_seq_len;
    
    ASSERT (*next + len <= data->len, "Error reading sam_line=%u: unexpected end of %s data", sam_line_i, field_name);

    buf_add (&vb->txt_data, &data->data[*next], len);
    buf_add (&vb->txt_data, "\t", 1);
    
    *next += len;
}

static void piz_sam_reconstruct_vb (VBlockSAM *vb)
{
    START_TIMER;

    buf_alloc (vb, &vb->txt_data, vb->vb_data_size, 1.1, "piz_sam_reconstruct_vb", vb->vblock_i);

    uint32_t cigar_seq_len=0, snip_len, rname_snip_len;
    const char *snip, *rname_snip;
    char pos_str[50], pnext_str[50];

    for (uint32_t vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        for (SamFields f=SAM_QNAME; f <= SAM_OPTIONAL; f++) {
                
            uint32_t word_index = mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[f], NULL, &snip, &snip_len, vb->first_line + vb_line_i);

            switch (f) {

            // reconstruct QNAME from its subfield components
            case SAM_QNAME:
                piz_sam_reconstruct_qname (vb, snip, snip_len, vb->first_line + vb_line_i);
                break;    

            // reconstruct pos from delta vs. previous line
            case SAM_POS:
                vb->last_pos = piz_decode_pos (vb->last_pos, snip, snip_len, pos_str, &snip_len); 
                snip = pos_str;
                break;

            // rnext & pnext
            case SAM_PNEXT:

                // rnext - get from RNAME dictionary
                mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[SAM_RNAME], NULL, &rname_snip, &rname_snip_len, vb->first_line + vb_line_i);                
                buf_add (&vb->txt_data, rname_snip, rname_snip_len);
                buf_add (&vb->txt_data, "\t", 1);

                // pnext - reconstruct from delta vs. pos, unless its "*" which means the original was "0"
                if (snip[0] == '*' && snip_len == 1)
                    snip = "0";
                else {
                    piz_decode_pos (vb->last_pos, snip, snip_len, pnext_str, &snip_len); 
                    snip = pnext_str;
                }
                break;

            // cigar - get length of SEQ and QUAL, and if original CIGAR was "*" - recover it
            case SAM_CIGAR:
                cigar_seq_len = seg_sam_seq_len_from_cigar (snip, snip_len);
                
                if (snip[snip_len-1] == '*') {
                    snip = "*";
                    snip_len = 1;
                }
                break;

            case SAM_OPTIONAL:

                // SEQ & QUAL data
                piz_sam_reconstruct_seq_qual (vb, cigar_seq_len, &vb->seq_data,  &vb->next_seq,  "SEQ",  vb->first_line + vb_line_i);
                piz_sam_reconstruct_seq_qual (vb, cigar_seq_len, &vb->qual_data, &vb->next_qual, "QUAL", vb->first_line + vb_line_i);

                // add OPTIONAL fields, and Windows-style \r if needed
                piz_sam_reconstruct_optional_fields (vb, snip, snip_len, word_index, vb->first_line + vb_line_i);
                break;

            default: break; // FLAG, RNAME, MAPQ, TLEN
            }
            
            if (f != SAM_QNAME && f != SAM_OPTIONAL)
                buf_add (&vb->txt_data, snip, snip_len);

            buf_add (&vb->txt_data, f==SAM_OPTIONAL ? "\n" : "\t", 1);
        }
    }

    COPY_TIMER(vb->profile.piz_reconstruct_vb);
}

static void piz_sam_uncompress_all_sections (VBlockSAM *vb)
{
    // The VB is read from disk in zfile_sam_read_one_vb(), in the I/O thread, and is decompressed here in the 
    // Compute thread, with the exception of dictionaries that are processed by the I/O thread
    // Order of sections in a V2 VB:
    // 1. SEC_SAM_VB_HEADER - no data, header only
    // 2. (the dictionaries were here in the file, but they are omitted from vb->z_data)
    // 3. b250 of SAM_QNAME to SAM_OPTIONAL
    // 4. b250 of all QNAME subfields
    // 5. b250 of all OPTIONAL subfields
    // 6. SEQ data
    // 7. QUAL data

    unsigned *section_index = (unsigned *)vb->z_section_headers.data;

    SectionHeaderVbHeaderSAM *header = (SectionHeaderVbHeaderSAM *)(vb->z_data.data + section_index[0]);
    vb->first_line              = BGEN32 (header->first_line);
    vb->num_lines               = BGEN32 (header->num_lines);
    vb->vb_data_size            = BGEN32 (header->vb_data_size);

    // in case of --split, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every sam component
    if (flag_split) vb->vblock_i = BGEN32 (header->h.vblock_i);
    
    unsigned section_i=1;

    // uncompress the fields     
    for (SamFields f=SAM_QNAME; f <= SAM_OPTIONAL; f++) {

        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);

        zfile_uncompress_section ((VBlockP)vb, header, &vb->mtf_ctx[f].b250, "mtf_ctx.b250", SEC_SAM_QNAME_B250 + f*2);
    }

    // QNAME subfields
    for (uint8_t sf_i=0; sf_i < vb->qname_mapper.num_subfields ; sf_i++) {
        
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);

        MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, NULL, header->dict_id, SEC_SAM_QNAME_SF_DICT);

        vb->qname_mapper.did_i[sf_i] = ctx->did_i;

        zfile_uncompress_section ((VBlockP)vb, header, &ctx->b250, "mtf_ctx.b250", SEC_SAM_QNAME_SF_B250);    
    }

    // OPTIONAL subfields
    for (uint8_t sf_i=0; sf_i < vb->num_optional_subfield_b250s ; sf_i++) {
        
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);

        MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, NULL, header->dict_id, SEC_SAM_OPTNL_SF_DICT);

        zfile_uncompress_section ((VBlockP)vb, header, &ctx->b250, "mtf_ctx.b250", SEC_SAM_OPTNL_SF_B250);    
    }

    piz_sam_map_optional_subfields (vb);

    // SEQ and QUAL
    SectionHeader *seq_header  = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, seq_header, &vb->seq_data, "seq_data", SEC_SAM_SEQ_DATA);    
    
    SectionHeader *qual_header = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, qual_header, &vb->qual_data, "qual_data", SEC_SAM_QUAL_DATA);    
}

void piz_sam_uncompress_one_vb (VBlock *vb_)
{
    START_TIMER;

    VBlockSAM *vb = (VBlockSAM *)vb_;

    piz_sam_uncompress_all_sections (vb);

    piz_sam_reconstruct_vb (vb);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway
    
    COPY_TIMER (vb->profile.compute);
}
