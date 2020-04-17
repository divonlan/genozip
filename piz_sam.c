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
#include "move_to_front.h"
#include "regions.h"
#include "zfile.h"

#   define LOAD_SNIP(did_i) mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[(did_i)], NULL, &snip, &snip_len, sam_line_i); 
#   define RECONSTRUCT_FROM_DICT(f)  \
        LOAD_SNIP(f) \
        buf_add (&vb->reconstructed_line, snip, snip_len); \
        buf_add (&vb->reconstructed_line, "\t", 1); 

static void piz_sam_reconstruct_pos (VBlockSAM *vb, uint32_t sam_line_i, char separator, bool skip)
{
    unsigned i=vb->next_random_pos; for (; i < vb->random_pos_data.len && vb->random_pos_data.data[i] != '\t'; i++);
    
    ASSERT (i < vb->random_pos_data.len, "Error reading sam_line=%u: unexpected end of RANDOM_POS data", sam_line_i);

    if (!skip) {
        buf_add (&vb->reconstructed_line, &vb->random_pos_data.data[vb->next_random_pos], i - vb->next_random_pos + (separator=='\t'));
        if (separator != '\t') buf_add (&vb->reconstructed_line, &separator, 1);
    }

    vb->next_random_pos = i+1; // skip the trailing \t too
}

static void piz_sam_reconstruct_seq_qual (VBlockSAM *vb, uint32_t cigar_seq_len, 
                                          const Buffer *data, uint32_t *next,
                                          const char *field_name, uint32_t sam_line_i, bool skip)
{
    // seq and qual are expected to be either of length cigar_seq_len, or "*" 
    uint32_t len = (*next >= data->len || data->data[*next] == '*') ? 1 : cigar_seq_len;
    ASSERT (*next + len <= data->len, "Error reading sam_line=%u: unexpected end of %s data", sam_line_i, field_name);

    if (!skip) buf_add (&vb->reconstructed_line, &data->data[*next], len);
    
    *next += len;
}

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

    // store the did_i of the NM field, if we have it
    vb->nm_did_i     = mtf_get_existing_did_i_by_dict_id ((VBlockP)vb, (DictIdType)dict_id_OPTION_NM); 
    vb->strand_did_i = mtf_get_existing_did_i_by_dict_id ((VBlockP)vb, (DictIdType)dict_id_OPTION_STRAND); 
}

static void piz_sam_reconstruct_qname (VBlockSAM *vb, const char *template, unsigned template_len, uint32_t sam_line_i)
{
    for (unsigned i=0; i <= template_len; i++) {
        
        const char *snip;
        unsigned snip_len;
        mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[vb->qname_mapper.did_i[i]], NULL, &snip, &snip_len, sam_line_i);

        buf_add (&vb->reconstructed_line, snip, snip_len);

        buf_add (&vb->reconstructed_line, (i < template_len) ? &template[i] : "\t", 1); // add separator
    }
}

static void piz_sam_reconstruct_sa_oa_xa (VBlockSAM *vb, bool is_xa, uint32_t sam_line_i)
{
    // XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
    // OA and SA format is: (rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ . in OA - NM is optional (but its , is not)

#   define ADD_SNIP(start,sep) { if (!flag_strip) { \
                                    buf_add (&vb->reconstructed_line, &snip[start], snip_len-start); \
                                    buf_add (&vb->reconstructed_line, sep, 1); \
                                 }\
                               }

#   define GET_ADD_SNIP(did_i, sep) { mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[did_i], NULL, &snip, &snip_len, sam_line_i);\
                                      ADD_SNIP(0, sep); }

    const char *snip; unsigned snip_len;

    GET_ADD_SNIP (SAM_RNAME, ","); // rname - consume but don't output - shares the same dictionary as main RNAME

    // strand
    if (!flag_strip) {
        mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[vb->strand_did_i], NULL, &snip, &snip_len, sam_line_i);\

        if (is_xa) buf_add (&vb->reconstructed_line, snip, 1); // XA: add strand concatenated with pos
    }

    // pos
    piz_sam_reconstruct_pos (vb, sam_line_i, ',', flag_strip); // if flag_strip - consume but don't output
    
    if (!is_xa) ADD_SNIP(0, ","); // SA and OA: add strand now

    GET_ADD_SNIP (SAM_CIGAR, ","); // cigar (read even if flag_snip because consumes same CIGAR as main)

    if (!flag_strip) {
        if (!is_xa) GET_ADD_SNIP (SAM_MAPQ, ","); // SA and OA only: mapq

        // nm
        GET_ADD_SNIP (vb->nm_did_i, ";");
    }
}

static void piz_sam_reconstruct_optional_fields (VBlockSAM *vb, uint32_t cigar_seq_len, 
                                                 const char *oname, unsigned oname_len, uint32_t opt_word_index, 
                                                 uint32_t sam_line_i)
{

    SubfieldMapper *opt_map = ENT (SubfieldMapper, vb->optional_mapper_buf, opt_word_index);

    ASSERT (opt_map->num_subfields == oname_len / 5, "Error: opt_map->num_subfields=%u but oname=%.*s indicates %u optional fields. sam_line=%u", 
            opt_map->num_subfields, oname_len, oname, oname_len/5, sam_line_i);

    const char *snip; unsigned snip_len;
    for (unsigned sf_i=0; sf_i < opt_map->num_subfields; sf_i++) {

        if (!flag_strip) buf_add (&vb->reconstructed_line, &oname[sf_i*5], 5)

        DictIdType dict_id = dict_id_sam_optnl_sf (dict_id_make (&oname[sf_i*5], 2));
        MtfContext *ctx = &vb->mtf_ctx[opt_map->did_i[sf_i]];

        // E2 doesn't have a dictionary - its data is stored in SEQ
        if (dict_id.num == dict_id_OPTION_E2)
            piz_sam_reconstruct_seq_qual (vb, cigar_seq_len, &vb->seq_data,  &vb->next_seq, "E2", sam_line_i, flag_strip);

        // U2 doesn't have a dictionary - its data is stored in QUAL
        else if (dict_id.num == dict_id_OPTION_U2) {
            if (!flag_strip) piz_sam_reconstruct_seq_qual (vb, cigar_seq_len, &vb->qual_data, &vb->next_qual, "U2", sam_line_i, false);
        }

        // MC and OC are stored in the CIGAR dictionary
        else if (dict_id.num == dict_id_OPTION_MC || dict_id.num == dict_id_OPTION_OC) {
            if (!flag_strip) {
                mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[SAM_CIGAR], NULL, &snip, &snip_len, sam_line_i);
                buf_add (&vb->reconstructed_line, snip, snip_len);
            }
        }

        else if (dict_id.num == dict_id_OPTION_mc && oname[sf_i*5 + 3] == 'i') {
            if (!flag_strip) {
                LOAD_SNIP (ctx->did_i);
                char mc_str[30];
                piz_decode_pos (vb->last_pos, snip, snip_len, mc_str, &snip_len); 
                buf_add (&vb->reconstructed_line, mc_str, snip_len);
            }
        }
        
        // SA, XA and OA have subsubfields IF snip starts with ascii 255
        else if (dict_id.num == dict_id_OPTION_SA || dict_id.num == dict_id_OPTION_OA ||
                 (dict_id.num == dict_id_OPTION_XA && oname[sf_i*5 + 3] == 'Z')) {

            mtf_get_next_snip ((VBlockP)vb, ctx, NULL, &snip, &snip_len, sam_line_i);

            if (snip[0] == -1) {
                unsigned repeats = atoi (&snip[1]);
                for (unsigned r=0; r < repeats; r++)             
                    piz_sam_reconstruct_sa_oa_xa (vb, (dict_id.num == dict_id_OPTION_XA), sam_line_i);
            }
            else if (!flag_strip)
                buf_add (&vb->reconstructed_line, snip, snip_len);
        }

        // other optional fields - get each for its own dictionary
        else {        
            if (!flag_strip) {
                mtf_get_next_snip ((VBlockP)vb, ctx, NULL, &snip, &snip_len, sam_line_i);
                buf_add (&vb->reconstructed_line, snip, snip_len);
            }
        }

        if (!flag_strip && sf_i != opt_map->num_subfields-1)
            buf_add (&vb->reconstructed_line, "\t", 1);
    }
    
    if (oname_len % 5 == 1 && oname[oname_len-1] == '#') // Windows style end-of-line \r\n
        buf_add (&vb->reconstructed_line, "\r", 1);
}

static void piz_sam_reconstruct_vb (VBlockSAM *vb)
{
#   define IFNOTSTRIP(def,len)  if (flag_strip) { \
                                buf_add (&vb->reconstructed_line, def, len)  \
                                buf_add (&vb->reconstructed_line, "\t", 1) \
                            } else 
    START_TIMER;

    buf_alloc (vb, &vb->txt_data, vb->vb_data_size, 1.1, "piz_sam_reconstruct_vb", vb->vblock_i);
    buf_alloc (vb, &vb->reconstructed_line, vb->longest_line_len, 1, "reconstructed_line", vb->vblock_i);

    uint32_t cigar_seq_len=0, snip_len;
    const char *snip;
    uint32_t last_rname_word_index = (uint32_t)-1;
    char pos_str[50], pnext_str[50];
    
    for (uint32_t vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        vb->reconstructed_line.len = 0; // initialize for a new line

        uint32_t sam_line_i = vb->first_line + vb_line_i;

        // QNAME - reconstruct from its subfield components
        IFNOTSTRIP("*",1) {
            LOAD_SNIP (SAM_QNAME);
            piz_sam_reconstruct_qname (vb, snip, snip_len, vb->first_line + vb_line_i);
        }

        // FLAG, RNAME - from their dictionaries
        IFNOTSTRIP("0",1) { RECONSTRUCT_FROM_DICT (SAM_FLAG);}

        uint32_t this_rname_word_index = RECONSTRUCT_FROM_DICT (SAM_RNAME);

        // POS - reconstruct from pos_data
        if (this_rname_word_index == last_rname_word_index) {
            LOAD_SNIP (SAM_POS);
            vb->last_pos = piz_decode_pos (vb->last_pos, snip, snip_len, pos_str, &snip_len); 
            buf_add (&vb->reconstructed_line, pos_str, snip_len);
            buf_add (&vb->reconstructed_line, "\t", 1);
        }
        else { // different rname - get from random_pos
            vb->last_pos = seg_pos_snip_to_int (&vb->random_pos_data.data[vb->next_random_pos], vb_line_i, "POS");
            piz_sam_reconstruct_pos (vb, vb->first_line + vb_line_i, '\t', false);
        }
        last_rname_word_index = this_rname_word_index;

        // MAPQ - from its dictionary
        IFNOTSTRIP("255",3) { RECONSTRUCT_FROM_DICT (SAM_MAPQ); }

        // CIGAR - get length of SEQ and QUAL, and if original CIGAR was "*" - recover it
        LOAD_SNIP (SAM_CIGAR);
        cigar_seq_len = seg_sam_seq_len_from_cigar (snip, snip_len);

        IFNOTSTRIP("*",1) {         
            if (snip[snip_len-1] == '*') {  // if its something like "151*" make it "*"
                buf_add (&vb->reconstructed_line, "*\t", 2);
            }
            else {
                buf_add (&vb->reconstructed_line, snip, snip_len); 
                buf_add (&vb->reconstructed_line, "\t", 1);
            }
        }
        
        // RNEXT - from RNAME dictionary
        uint32_t this_rnext_word_index = LOAD_SNIP (SAM_RNAME) // always read, because it consumes the same b250 as RNAME
        IFNOTSTRIP("*",1) {         
            buf_add (&vb->reconstructed_line, snip, snip_len); 
            buf_add (&vb->reconstructed_line, "\t", 1); 
        }
        
        // PNEXT - 2.5 options:
        // 1. if RNEXT is the same as RNAME or unknown (i.e. either it is "=" or "*" or the word index is the same) - 
        //    we decode the delta from the RNEXT b250
        //    1A. If 1, but PNEXT is "*" we recover the "0"
        // 2. If RNEXT and RNAME are different - get the pos data from from vb->random_pos_data 
        bool get_from_dictionary = (this_rnext_word_index == this_rname_word_index) ||
                                   (snip_len == 1 && snip[0] == '=') || // this refers to the RNEXT snip
                                   (snip_len == 1 && snip[0] == '*');

        if (get_from_dictionary) {
            IFNOTSTRIP("0",1) {         
                LOAD_SNIP (SAM_PNEXT);

                if (snip_len == 1 && snip[0] == '*') {  // this refers to the PNEXT snip
                    buf_add (&vb->reconstructed_line, "0\t", 2);
                }
                else { 
                    piz_decode_pos (vb->last_pos, snip, snip_len, pnext_str, &snip_len); 
                    buf_add (&vb->reconstructed_line, pnext_str, snip_len);
                    buf_add (&vb->reconstructed_line, "\t", 1);
                }
            }
        }
        else {
            // note: we always consume (but possible skip outputting), even if flag_strip because it consumes the same random_pos as POS
            piz_sam_reconstruct_pos (vb, vb->first_line + vb_line_i, '\t', flag_strip);
            IFNOTSTRIP("0",1) {}; // add the default if we're stripping
        }
        // TLEN - from its dictionary
        IFNOTSTRIP("0",1) { RECONSTRUCT_FROM_DICT (SAM_TLEN); }
        
        // SEQ & QUAL data
        piz_sam_reconstruct_seq_qual (vb, cigar_seq_len, &vb->seq_data,  &vb->next_seq,  "SEQ",  vb->first_line + vb_line_i, false);
        buf_add (&vb->reconstructed_line, "\t", 1);

        IFNOTSTRIP("*",1) {
            piz_sam_reconstruct_seq_qual (vb, cigar_seq_len, &vb->qual_data, &vb->next_qual, "QUAL", vb->first_line + vb_line_i, false);
        }

        // OPTIONAL fields, and Windows-style \r if needed
        uint32_t word_index = LOAD_SNIP (SAM_OPTIONAL);
        if (snip_len != 1 || (snip[0] != '#' && snip[0] != '*')) buf_add (&vb->reconstructed_line, "\t", 1);
        
        piz_sam_reconstruct_optional_fields (vb, cigar_seq_len, snip, snip_len, word_index, vb->first_line + vb_line_i);

        buf_add (&vb->reconstructed_line, "\n", 1);

        // output the reconstructed line - unless it needs to be skipped according to --regions
        if (!flag_regions || regions_is_site_included (last_rname_word_index, vb->last_pos)) 
            buf_add (&vb->txt_data, vb->reconstructed_line.data, vb->reconstructed_line.len);
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
    // 6. RANDOM POS data
    // 7. SEQ data
    // 8. QUAL data

    unsigned *section_index = (unsigned *)vb->z_section_headers.data;

    SectionHeaderVbHeaderSAM *header = (SectionHeaderVbHeaderSAM *)(vb->z_data.data + section_index[0]);
    vb->first_line       = BGEN32 (header->first_line);
    vb->num_lines        = BGEN32 (header->num_lines);
    vb->vb_data_size     = BGEN32 (header->vb_data_size);
    vb->longest_line_len = BGEN32 (header->longest_line_len);

    // in case of --split, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every sam component
    if (flag_split) vb->vblock_i = BGEN32 (header->h.vblock_i);
    
    unsigned section_i=1;

    // uncompress the fields     
    for (SamFields f=SAM_QNAME; f <= SAM_OPTIONAL; f++) {

        SectionType b250_sec = SEC_SAM_QNAME_B250 + f*2;

        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);
        if (zfile_skip_section_by_flags (b250_sec, DICT_ID_NONE)) continue;

        zfile_uncompress_section ((VBlockP)vb, header, &vb->mtf_ctx[f].b250, "mtf_ctx.b250", b250_sec);
    }

    // QNAME subfields
    for (uint8_t sf_i=0; sf_i < vb->qname_mapper.num_subfields ; sf_i++) {
        
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);
        if (flag_strip) continue;

        MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, NULL, header->dict_id, SEC_SAM_QNAME_SF_DICT);

        vb->qname_mapper.did_i[sf_i] = ctx->did_i;

        zfile_uncompress_section ((VBlockP)vb, header, &ctx->b250, "mtf_ctx.b250", SEC_SAM_QNAME_SF_B250);    
    }

    // OPTIONAL subfields
    for (uint8_t sf_i=0; sf_i < vb->num_optional_subfield_b250s ; sf_i++) {
        
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);
        if (zfile_skip_section_by_flags (SEC_SAM_OPTNL_SF_B250, header->dict_id)) continue;

        MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, NULL, header->dict_id, SEC_SAM_OPTNL_SF_DICT);

        zfile_uncompress_section ((VBlockP)vb, header, &ctx->b250, "mtf_ctx.b250", SEC_SAM_OPTNL_SF_B250);    
    }

    piz_sam_map_optional_subfields (vb);

    // POS, SEQ and QUAL (data also contains E2 and U2 respectively, if they exist)
    SectionHeader *pos_header  = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, pos_header, &vb->random_pos_data, "random_pos_data", SEC_SAM_RAND_POS_DATA);    
    
    SectionHeader *seq_header  = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, seq_header, &vb->seq_data, "seq_data", SEC_SAM_SEQ_DATA);    

    SectionHeader *qual_header = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    if (!flag_strip) zfile_uncompress_section ((VBlockP)vb, qual_header, &vb->qual_data, "qual_data", SEC_SAM_QUAL_DATA);    
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
