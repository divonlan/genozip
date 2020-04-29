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
#include "file.h"
#include "strings.h"

// reconstructs POS from random_pos_data, and returns a pointer to the reconstructed pos string (terminated by separator) 
static uint32_t piz_sam_reconstruct_random_pos (VBlockSAM *vb, uint32_t txt_line_i, char separator, bool skip)
{
    ASSERT (vb->next_random_pos < vb->random_pos_data.len, "Error reading sam_line=%u: unexpected end of RANDOM_POS data", txt_line_i);
    
    uint32_t pos = BGEN32 (*ENT (uint32_t, vb->random_pos_data, vb->next_random_pos++));
    
    if (!skip) {
        char pos_str[20];
        unsigned pos_str_len;
        str_uint (pos, pos_str, &pos_str_len);

        buf_add (&vb->txt_data, pos_str, pos_str_len);
        buf_add (&vb->txt_data, &separator, 1);
    }

    return pos;
}

static void piz_sam_reconstruct_tlen (VBlockSAM *vb, const char *tlen, unsigned tlen_len)
{
    ASSERT0 (tlen, "Error in piz_sam_reconstruct_tlen: tlen=NULL");
    ASSERT0 (tlen_len, "Error in piz_sam_reconstruct_tlen: tlen_len=0");

    // case - tlen is the negative of the previous line
    if ((tlen_len == 1) && (tlen[0] == '*')) {
        ASSERT0 (vb->last_tlen_abs && vb->last_tlen_abs_len, "Error in piz_sam_reconstruct_tlen: tlen=\"*\" but this is the first line in the vb");

        // if previous line was a positive, this line is negative
        if (vb->last_tlen_is_positive) 
            buf_add (&vb->txt_data, "-", 1); 

        buf_add (&vb->txt_data, vb->last_tlen_abs, vb->last_tlen_abs_len); 

        // flip the sign of last_tlen
        vb->last_tlen_is_positive = !vb->last_tlen_is_positive;
    }

    // case - tlen is not a negative of previous line - output as-is
    else {
        buf_add (&vb->txt_data, tlen, tlen_len); 
     
        // update last_tlen_*
        vb->last_tlen_is_positive = tlen[0] != '-';
        vb->last_tlen_abs         = &tlen[!vb->last_tlen_is_positive];
        vb->last_tlen_abs_len     = tlen_len - !vb->last_tlen_is_positive;
    }

    buf_add (&vb->txt_data, "\t", 1); 
}

static inline void piz_sam_reconstruct_AS (VBlockSAM *vb, const char *snip, unsigned snip_len,
                                           uint32_t txt_line_i, uint32_t cigar_seq_len)
{
    if (snip[0] == '*') { // delta vs seq_len
        unsigned value = cigar_seq_len - atoi (&snip[1]);
        char value_str[20];
        unsigned value_str_len;
        str_uint (value, value_str, &value_str_len);
        buf_add (&vb->txt_data, value_str, value_str_len);
    }
    else // not delta encoded
        buf_add (&vb->txt_data, snip, snip_len);
}

static inline void piz_sam_reconstruct_MD (VBlockSAM *vb, uint32_t txt_line_i, uint32_t cigar_seq_len)
{
    const char *snip; 
    unsigned snip_len;
    LOAD_SNIP_FROM_BUF (vb->md_data, vb->next_md, "MD", '\t');

    char reconstruced_md_str[MAX_SAM_MD_LEN]; 
    unsigned reconstruced_md_str_len;

    // case: MD is an empty string - reconstruct the original MD that is the sequence length
    if (!snip_len) {
        str_uint (cigar_seq_len, reconstruced_md_str, &reconstruced_md_str_len);
        buf_add (&vb->txt_data, reconstruced_md_str, reconstruced_md_str_len);
    }
    
    // case: MD ends with a * eg "119C*" - we reconstruct the original, eg "119C31" - using the sequence length 
    else if (snip[snip_len-1] == '*') {
        unsigned partial_seq_len_by_md_field = seg_sam_get_seq_len_by_MD_field (snip, snip_len-1, NULL);

        memcpy (reconstruced_md_str, snip, snip_len-1);
        str_uint (cigar_seq_len - partial_seq_len_by_md_field, 
                                    &reconstruced_md_str[snip_len-1], &reconstruced_md_str_len);
        buf_add (&vb->txt_data, reconstruced_md_str, snip_len-1 + reconstruced_md_str_len);
    }
    
    // case: MD is stored as-is - just copy it
    else 
        buf_add (&vb->txt_data, snip, snip_len);
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
                        
            DictIdType dict_id = dict_id_sam_optnl_sf (dict_id_make (&optional[i*5], 4)); // get dict_id from first 4 characters eg "MX:i"
            optional_mapper->did_i[i] = mtf_get_existing_did_i_by_dict_id (dict_id); 
        }
    }

    // store the did_i of the NM field, if we have it
    vb->nm_did_i     = mtf_get_existing_did_i_by_dict_id ((DictIdType)dict_id_OPTION_NM); 
    vb->strand_did_i = mtf_get_existing_did_i_by_dict_id ((DictIdType)dict_id_OPTION_STRAND); 
}

static void piz_sam_reconstruct_SA_OA_XA (VBlockSAM *vb, bool is_xa, uint32_t txt_line_i)
{
    // XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
    // OA and SA format is: (rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ . in OA - NM is optional (but its , is not)

#   define ADD_SNIP(start,sep) { if (!flag_strip) { \
                                    buf_add (&vb->txt_data, &snip[start], snip_len-start); \
                                    buf_add (&vb->txt_data, sep, 1); \
                                 }\
                               }

#   define GET_ADD_SNIP(did_i, sep) { mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[did_i], NULL, &snip, &snip_len, txt_line_i);\
                                      ADD_SNIP(0, sep); }

    const char *snip; unsigned snip_len;

    GET_ADD_SNIP (SAM_RNAME, ","); // rname - consume but don't output - shares the same dictionary as main RNAME

    // strand
    if (!flag_strip) {
        mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[vb->strand_did_i], NULL, &snip, &snip_len, txt_line_i);\

        if (is_xa) buf_add (&vb->txt_data, snip, 1); // XA: add strand concatenated with pos
    }

    // pos
    piz_sam_reconstruct_random_pos (vb, txt_line_i, ',', flag_strip); // if flag_strip - consume but don't output
    
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
                                                 uint32_t txt_line_i)
{

    SubfieldMapper *opt_map = ENT (SubfieldMapper, vb->optional_mapper_buf, opt_word_index);

    ASSERT (opt_map->num_subfields == oname_len / 5, "Error: opt_map->num_subfields=%u but oname=%.*s indicates %u optional fields. sam_line=%u", 
            opt_map->num_subfields, oname_len, oname, oname_len/5, txt_line_i);

    const char *snip; unsigned snip_len;
    for (unsigned sf_i=0; sf_i < opt_map->num_subfields; sf_i++) {

        if (!flag_strip) buf_add (&vb->txt_data, &oname[sf_i*5], 5)

        DictIdType dict_id = dict_id_sam_optnl_sf (dict_id_make (&oname[sf_i*5], 4));
        uint8_t did_i = opt_map->did_i[sf_i];
        MtfContext *ctx = did_i != DID_I_NONE ? &vb->mtf_ctx[did_i] : NULL; // NULL if this subfield has no dictionary (eg MD, E2, U2)

        // E2 doesn't have a dictionary - its data is stored in SEQ
        if (dict_id.num == dict_id_OPTION_E2)
            piz_reconstruct_seq_qual ((VBlockP)vb, cigar_seq_len, &vb->seq_data, &vb->next_seq, SEC_SEQ_DATA, txt_line_i, false);

        // U2 doesn't have a dictionary - its data is stored in QUAL
        else if (dict_id.num == dict_id_OPTION_U2) {
            if (!flag_strip) piz_reconstruct_seq_qual ((VBlockP)vb, cigar_seq_len, &vb->qual_data, &vb->next_qual, SEC_QUAL_DATA, txt_line_i, false);
        }

        else if (dict_id.num == dict_id_OPTION_MD) {
            if (!flag_strip) 
                piz_sam_reconstruct_MD (vb, txt_line_i, cigar_seq_len);
        }

        // AS is normally stored as a delta vs seq_len
        else if (dict_id.num == dict_id_OPTION_AS) {
            if (!flag_strip) {
                LOAD_SNIP (did_i);
                piz_sam_reconstruct_AS (vb, snip, snip_len, txt_line_i, cigar_seq_len);
            }
        }

        // MC and OC are stored in the CIGAR dictionary
        else if (dict_id.num == dict_id_OPTION_MC || dict_id.num == dict_id_OPTION_OC) {
            if (!flag_strip) {
                mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[SAM_CIGAR], NULL, &snip, &snip_len, txt_line_i);
                buf_add (&vb->txt_data, snip, snip_len);
            }
        }

        else if (dict_id.num == dict_id_OPTION_mc) {
            if (!flag_strip) 
                RECONSTRUCT_FROM_DICT_POS (did_i, false, NULL, false);
        }
        
        // SA, XA and OA have subsubfields IF snip starts with ascii 255
        else if (dict_id.num == dict_id_OPTION_SA || dict_id.num == dict_id_OPTION_OA || dict_id.num == dict_id_OPTION_XA) {

            mtf_get_next_snip ((VBlockP)vb, ctx, NULL, &snip, &snip_len, txt_line_i);

            if (snip[0] == -1) {
                unsigned repeats = atoi (&snip[1]);
                for (unsigned r=0; r < repeats; r++)             
                    piz_sam_reconstruct_SA_OA_XA (vb, (dict_id.num == dict_id_OPTION_XA), txt_line_i);
            }
            else if (!flag_strip)
                buf_add (&vb->txt_data, snip, snip_len);
        }

        // other optional fields - get each for its own dictionary
        else {        
            if (!flag_strip) {
                mtf_get_next_snip ((VBlockP)vb, ctx, NULL, &snip, &snip_len, txt_line_i);
                buf_add (&vb->txt_data, snip, snip_len);
            }
        }

        if (!flag_strip && sf_i != opt_map->num_subfields-1)
            buf_add (&vb->txt_data, "\t", 1);
    }
    
    if (oname_len % 5 == 1 && oname[oname_len-1] == '#') // Windows style end-of-line \r\n
        buf_add (&vb->txt_data, "\r", 1);
}

static void piz_sam_reconstruct_vb (VBlockSAM *vb)
{
    START_TIMER;

    buf_alloc (vb, &vb->txt_data, vb->vb_data_size, 1.1, "txt_data", vb->vblock_i);

    uint32_t cigar_seq_len=0, snip_len;
    const char *snip;
    
    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_data_start = vb->txt_data.len;
        uint32_t txt_line_i     = vb->first_line + vb_line_i;

        // QNAME - reconstruct from its subfield components
        IFNOTSTRIP("*",1) {
            LOAD_SNIP (SAM_QNAME);
            piz_reconstruct_compound_field ((VBlockP)vb, &vb->qname_mapper, "\t", 1, snip, snip_len, txt_line_i);
        }

        // FLAG, RNAME - from their dictionaries
        IFNOTSTRIP("0",1) { RECONSTRUCT_FROM_DICT (SAM_FLAG);}

        uint32_t rname_word_index = RECONSTRUCT_FROM_DICT (SAM_RNAME);

        // POS - reconstruct from pos_data
        // chr5 5000000   rname_minus_2=chr5          rname_minus_3=chr5
        // chr1 1000001   rname_minus_1=chr1          rname_minus_2=chr1
        // chr1 1000002   rname_minus_1=chr1          rname=chr1 <--- DO delta
        // chr1 2000001   rname=chr1 <-- DONT delta
        // chr1 2000002
        if (rname_word_index == vb->rname_index_minus_1 && 
            (rname_word_index != vb->rname_index_minus_2 || rname_word_index == vb->rname_index_minus_3)) {
            RECONSTRUCT_FROM_DICT_POS (SAM_POS, true, NULL, true); // same rname - reconstruct from delta 
        }
        else  // different rname - get from random_pos
            vb->last_pos = piz_sam_reconstruct_random_pos (vb, txt_line_i, '\t', false);

        vb->rname_index_minus_3 = vb->rname_index_minus_2;
        vb->rname_index_minus_2 = vb->rname_index_minus_1;
        vb->rname_index_minus_1 = rname_word_index;

        // MAPQ - from its dictionary
        IFNOTSTRIP("255",3) { RECONSTRUCT_FROM_DICT (SAM_MAPQ); }

        // CIGAR - get length of SEQ and QUAL, and if original CIGAR was "*" - recover it
        LOAD_SNIP (SAM_CIGAR);
        cigar_seq_len = seg_sam_seq_len_from_cigar (snip, snip_len);

        IFNOTSTRIP("*",1) {         
            if (snip[snip_len-1] == '*') {  // if its something like "151*" make it "*"
                buf_add (&vb->txt_data, "*\t", 2);
            }
            else {
                buf_add (&vb->txt_data, snip, snip_len); 
                buf_add (&vb->txt_data, "\t", 1);
            }
        }
        
        // RNEXT - from RNAME dictionary
        uint32_t this_rnext_word_index = LOAD_SNIP (SAM_RNAME) // always read, because it consumes the same b250 as RNAME
        IFNOTSTRIP("*",1) {         
            buf_add (&vb->txt_data, snip, snip_len); 
            buf_add (&vb->txt_data, "\t", 1); 
        }
        
        // PNEXT - 2.5 options:
        // 1. if RNEXT is the same as RNAME or unknown (i.e. either it is "=" or "*" or the word index is the same) - 
        //    we decode the delta from the RNEXT b250
        //    1A. If 1, but PNEXT is "*" we recover the "0"
        // 2. If RNEXT and RNAME are different - get the pos data from from vb->random_pos_data 
        bool get_from_dictionary = (this_rnext_word_index == rname_word_index) ||
                                   (snip_len == 1 && snip[0] == '=') || // this refers to the RNEXT snip
                                   (snip_len == 1 && snip[0] == '*');

        if (get_from_dictionary) {
            IFNOTSTRIP("0",1) {         
                LOAD_SNIP (SAM_PNEXT);

                if (snip_len == 1 && snip[0] == '*') {  // this refers to the PNEXT snip
                    buf_add (&vb->txt_data, "0\t", 2);
                }
                else { 
                    RECONSTRUCT_FROM_DICT_POS (DID_I_NONE, false, &vb->last_pnext_delta, true);
                }
            }
        }
        else {
            // note: we always consume (but possible skip outputting), even if flag_strip because it consumes the same random_pos as POS
            piz_sam_reconstruct_random_pos (vb, txt_line_i, '\t', flag_strip);
            IFNOTSTRIP("0",1) {}; // add the default if we're stripping
        }
        // TLEN - from its dictionary
        IFNOTSTRIP("0",1) {         
            LOAD_SNIP (SAM_TLEN);
            piz_sam_reconstruct_tlen (vb, snip, snip_len);
        }
        
        // SEQ & QUAL data
        piz_reconstruct_seq_qual ((VBlockP)vb, cigar_seq_len, &vb->seq_data, &vb->next_seq, SEC_SEQ_DATA, txt_line_i, false);
        buf_add (&vb->txt_data, "\t", 1);

        IFNOTSTRIP("*",1) {
            piz_reconstruct_seq_qual ((VBlockP)vb, cigar_seq_len, &vb->qual_data, &vb->next_qual, SEC_QUAL_DATA, txt_line_i, false);
        }

        // OPTIONAL fields, and Windows-style \r if needed
        uint32_t word_index = LOAD_SNIP (SAM_OPTIONAL);
        if (snip_len != 1 || (snip[0] != '#' && snip[0] != '*')) buf_add (&vb->txt_data, "\t", 1);
        
        piz_sam_reconstruct_optional_fields (vb, cigar_seq_len, snip, snip_len, word_index, txt_line_i);

        buf_add (&vb->txt_data, "\n", 1);

        // after consuming sections' data, if this line is not to be outputed - shorten txt_data back to start of line
        if (flag_regions && !regions_is_site_included (vb->rname_index_minus_1, vb->last_pos))
            vb->txt_data.len = txt_data_start; // remove excluded line
    }

    COPY_TIMER(vb->profile.piz_reconstruct_vb);
}

static void piz_sam_uncompress_all_sections (VBlockSAM *vb)
{
    // The VB is read from disk in zfile_sam_read_one_vb(), in the I/O thread, and is decompressed here in the 
    // Compute thread, with the exception of dictionaries that are processed by the I/O thread
    // Order of sections in a V2 VB:
    // 1. SEC_VB_HEADER - no data, header only
    // 2. (the dictionaries were here in the file, but they are omitted from vb->z_data)
    // 3. b250 of SAM_QNAME to SAM_OPTIONAL
    // 4. b250 of all QNAME subfields
    // 5. b250 of all OPTIONAL subfields
    // 6. RANDOM POS data
    // 7. SEQ data
    // 8. QUAL data

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
    
    // QNAME subfields
    piz_uncompress_compound_field ((VBlockP)vb, SEC_SAM_QNAME_B250, SEC_SAM_QNAME_SF_B250, &vb->qname_mapper, &section_i);

    // OPTIONAL subfields
    for (uint8_t sf_i=0; sf_i < vb->num_optional_subfield_b250s ; sf_i++) {
        
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]);
        if (zfile_is_skip_section (vb, SEC_SAM_OPTNL_SF_B250, header->dict_id)) continue;

        MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, NULL, header->dict_id, SEC_SAM_OPTNL_SF_DICT);

        zfile_uncompress_section ((VBlockP)vb, header, &ctx->b250, "mtf_ctx.b250", SEC_SAM_OPTNL_SF_B250);    
    }

    piz_sam_map_optional_subfields (vb);

    // RAND_POS, MD, SEQ and QUAL (data also contains E2 and U2 respectively, if they exist)
    SectionHeader *pos_header  = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, pos_header, &vb->random_pos_data, "random_pos_data", SEC_SAM_RAND_POS_DATA);    
    vb->random_pos_data.len /= 4; // its an array of uint32_t

    SectionHeader *md_header  = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, md_header, &vb->md_data, "md_data", SEC_SAM_MD_DATA);    
    
    SectionHeader *seq_header  = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, seq_header, &vb->seq_data, "seq_data", SEC_SEQ_DATA);    

    SectionHeader *qual_header = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    if (!flag_strip) zfile_uncompress_section ((VBlockP)vb, qual_header, &vb->qual_data, "qual_data", SEC_QUAL_DATA);    
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
