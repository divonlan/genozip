// ------------------------------------------------------------------
//   piz_sam.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// this is the compute thread entry point. It receives all data of a variant block and processes it
// in memory to the uncompressed format. This thread then terminates the I/O thread writes the output.

#include "vblock.h"
#include "piz.h"
#include "endianness.h"
#include "zfile.h"
#include "buffer.h"
#include "seg.h"
#include "move_to_front.h"
#include "regions.h"
#include "file.h"
#include "strings.h"

// reconstructs POS from random_pos_data, and returns a pointer to the reconstructed pos string (terminated by separator) 
static uint32_t piz_sam_reconstruct_random_pos (VBlockSAM *vb, char separator, bool skip)
{
    MtfContext *pos_ctx = &vb->mtf_ctx[DTF(pos)];
    ASSERT (pos_ctx->next_local < pos_ctx->local.len, "Error in piz_sam_reconstruct_random_pos: reading vb->line_i=%u: unexpected end of %s local data", vb->line_i, pos_ctx->name);
    uint32_t pos = BGEN32 (*ENT (uint32_t, pos_ctx->local, pos_ctx->next_local++));
    
    if (!skip) {
        char pos_str[20];
        unsigned pos_str_len;
        str_int (pos, pos_str, &pos_str_len);

        RECONSTRUCT (pos_str, pos_str_len);
        RECONSTRUCT1 (separator);
    }

    return pos;
}

static void piz_sam_reconstruct_tlen (VBlockSAM *vb, const char *tlen, unsigned tlen_len, int32_t pnext_pos_delta, int32_t cigar_seq_len)
{
    ASSERT0 (tlen, "Error in piz_sam_reconstruct_tlen: tlen=NULL");
    ASSERT0 (tlen_len, "Error in piz_sam_reconstruct_tlen: tlen_len=0");

    // case - tlen is the negative of the previous line
    if ((tlen_len == 1) && (tlen[0] == '*')) {
        ASSERT0 (vb->last_tlen_abs && vb->last_tlen_abs_len, "Error in piz_sam_reconstruct_tlen: tlen=\"*\" but this is the first line in the vb");

        // if previous line was a positive, this line is negative
        if (vb->last_tlen_is_positive) 
            RECONSTRUCT1 ('-'); 

        RECONSTRUCT (vb->last_tlen_abs, vb->last_tlen_abs_len); 

        // flip the sign of last_tlen
        vb->last_tlen_is_positive = !vb->last_tlen_is_positive;
    }

    else if ((tlen_len > 1) && (tlen[0] == '.')) {
        int32_t tlen_by_calc = atoi (&tlen[1]);
        int32_t tlen_val = tlen_by_calc + pnext_pos_delta + cigar_seq_len;

        // update last_tlen_*
        vb->last_tlen_is_positive = tlen_val >= 0;
        vb->last_tlen_abs = AFTERENT (char, vb->txt_data) + !vb->last_tlen_is_positive;

        unsigned reconstructed_tlen_len;
        str_int (tlen_val, AFTERENT (char, vb->txt_data), &reconstructed_tlen_len);
        vb->txt_data.len += reconstructed_tlen_len;
        
        vb->last_tlen_abs_len = reconstructed_tlen_len - !vb->last_tlen_is_positive;
    }

    // case - tlen is not a negative of previous line - output as-is
    else {
        RECONSTRUCT (tlen, tlen_len); 
     
        // update last_tlen_*
        vb->last_tlen_is_positive = tlen[0] != '-';
        vb->last_tlen_abs         = &tlen[!vb->last_tlen_is_positive];
        vb->last_tlen_abs_len     = tlen_len - !vb->last_tlen_is_positive;
    }

    RECONSTRUCT1 ('\t'); 
}

static inline void piz_sam_reconstruct_AS (VBlockSAM *vb, uint8_t did_i, uint32_t cigar_seq_len)
{
    const char *snip = AFTERENT (char, vb->txt_data);
    unsigned snip_len = piz_reconstruct_from_ctx ((VBlockP)vb, did_i, "", 0);

    if (snip[0] == '*') { // delta vs seq_len
        vb->txt_data.len -= snip_len; // roll back
        unsigned value = cigar_seq_len - atoi (&snip[1]);
        char value_str[20];
        unsigned value_str_len;
        str_int (value, value_str, &value_str_len);
        RECONSTRUCT (value_str, value_str_len);
    }
}

static inline void piz_sam_reconstruct_MD (VBlockSAM *vb, uint32_t cigar_seq_len)
{
    MtfContext *ctx = mtf_get_ctx (vb, (DictIdType)dict_id_OPTION_MD);

    DECLARE_SNIP;
    LOAD_SNIP_FROM_BUF (ctx->local, ctx->next_local, ctx->name);

    char reconstruced_md_str[MAX_SAM_MD_LEN]; 
    unsigned reconstruced_md_str_len;

    // case: MD is an empty string - reconstruct the original MD that is the sequence length
    if (!snip_len) {
        str_int (cigar_seq_len, reconstruced_md_str, &reconstruced_md_str_len);
        RECONSTRUCT (reconstruced_md_str, reconstruced_md_str_len);
    }
    
    // case: MD ends with a * eg "119C*" - we reconstruct the original, eg "119C31" - using the sequence length 
    else if (snip[snip_len-1] == '*') {
        unsigned partial_seq_len_by_md_field = seg_sam_get_seq_len_by_MD_field (snip, snip_len-1, NULL);

        memcpy (reconstruced_md_str, snip, snip_len-1);
        str_int (cigar_seq_len - partial_seq_len_by_md_field, 
                                    &reconstruced_md_str[snip_len-1], &reconstruced_md_str_len);
        RECONSTRUCT (reconstruced_md_str, snip_len-1 + reconstruced_md_str_len);
    }
    
    // case: MD is stored as-is - just copy it
    else 
        RECONSTRUCT (snip, snip_len);
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
//            optional_mapper->did_i[i] = mtf_get_existing_did_i_by_dict_id (dict_id); 
            optional_mapper->did_i[i] = mtf_get_ctx (vb, dict_id)->did_i; 
        }
    }
}

static void piz_sam_reconstruct_SA_OA_XA (VBlockSAM *vb, bool is_xa)
{
    // XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
    // OA and SA format is: (rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ . in OA - NM is optional (but its , is not)

    DECLARE_SNIP;
    piz_reconstruct_from_ctx ((VBlockP)vb, SAM_RNAME, ",", 1); // rname - shares the same dictionary as main RNAME

    // strand
    MtfContext *strand_ctx = mtf_get_ctx (vb, (DictIdType)dict_id_OPTION_STRAND);
    mtf_get_next_snip ((VBlockP)vb, strand_ctx, NULL, &snip, &snip_len);\

    if (is_xa) RECONSTRUCT (snip, 1); // XA: add strand concatenated with pos

    // pos
    piz_sam_reconstruct_random_pos (vb, ',', false); 
    
    if (!is_xa) RECONSTRUCT_SEP (snip, snip_len, ",");

    piz_reconstruct_from_ctx ((VBlockP)vb, SAM_CIGAR, ",", 1); // cigar (read even if flag_snip because consumes same CIGAR as main)

    if (!is_xa) piz_reconstruct_from_ctx ((VBlockP)vb, SAM_MAPQ, ",", 1); // SA and OA only: mapq

    // nm
    MtfContext *nm_ctx = mtf_get_ctx (vb, (DictIdType)dict_id_OPTION_NM);
    piz_reconstruct_from_ctx ((VBlockP)vb, nm_ctx->did_i, ";", 1);
}

static void piz_sam_reconstruct_optional_fields (VBlockSAM *vb, uint32_t cigar_seq_len, 
                                                 const char *oname, unsigned oname_len, uint32_t opt_word_index)
{

    SubfieldMapper *opt_map = ENT (SubfieldMapper, vb->optional_mapper_buf, opt_word_index);

    ASSERT (opt_map->num_subfields == oname_len / 5, "Error: opt_map->num_subfields=%u but oname=%.*s indicates %u optional fields. sam_line=%u", 
            opt_map->num_subfields, oname_len, oname, oname_len/5, vb->line_i);

    DECLARE_SNIP;
    char *this_bi=NULL, *this_bd=NULL;
    for (unsigned sf_i=0; sf_i < opt_map->num_subfields; sf_i++) {

        RECONSTRUCT (&oname[sf_i*5], 5)

        DictIdType dict_id = dict_id_sam_optnl_sf (dict_id_make (&oname[sf_i*5], 4));
        uint8_t did_i = opt_map->did_i[sf_i];
        MtfContext *ctx = did_i != DID_I_NONE ? &vb->mtf_ctx[did_i] : NULL; // NULL if this subfield has no dictionary (eg MD, E2, U2)

        // E2 doesn't have a dictionary - its data is stored in SEQ
        if (dict_id.num == dict_id_OPTION_E2)
            piz_reconstruct_seq_qual ((VBlockP)vb, &vb->mtf_ctx[SAM_SEQ], cigar_seq_len, false);

        // U2 doesn't have a dictionary - its data is stored in QUAL
        else if (dict_id.num == dict_id_OPTION_U2) 
            piz_reconstruct_seq_qual ((VBlockP)vb, &vb->mtf_ctx[SAM_QUAL], cigar_seq_len, false);

        else if (dict_id.num == dict_id_OPTION_MD) 
            piz_sam_reconstruct_MD (vb, cigar_seq_len);

        else if (dict_id.num == dict_id_OPTION_BD) {
            this_bd = AFTERENT (char, vb->txt_data);
            piz_reconstruct_seq_qual ((VBlockP)vb, ctx, cigar_seq_len, false);
        }

        else if (dict_id.num == dict_id_OPTION_BI) {
            this_bi = AFTERENT (char, vb->txt_data);
            piz_reconstruct_seq_qual ((VBlockP)vb, ctx, cigar_seq_len, false);
        }

        // AS is normally stored as a delta vs seq_len
        else if (dict_id.num == dict_id_OPTION_AS) 
            piz_sam_reconstruct_AS (vb, did_i, cigar_seq_len);

        // MC and OC are stored in the CIGAR dictionary
        else if (dict_id.num == dict_id_OPTION_MC || dict_id.num == dict_id_OPTION_OC) {
            mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[SAM_CIGAR], NULL, &snip, &snip_len);
            RECONSTRUCT (snip, snip_len);
        }

        else if (dict_id.num == dict_id_OPTION_mc)
            RECONSTRUCT_FROM_DICT_POS (did_i, vb->last_pos, false, NULL, false)
        
        // SA, XA and OA have subsubfields IF snip starts with ascii 255
        else if (dict_id.num == dict_id_OPTION_SA || dict_id.num == dict_id_OPTION_OA || dict_id.num == dict_id_OPTION_XA) {

            // temporary ugliness until we get Structed into piz_reconstruct_from_ctx

            const char *snip  = AFTERENT (char, vb->txt_data);
            uint32_t snip_len = piz_reconstruct_from_ctx ((VBlockP)vb, did_i, "", 0);

            if (snip[0] == SNIP_STRUCTURED) {
                vb->txt_data.len -= snip_len; // roll back

                unsigned repeats = atoi (&snip[1]);
                for (unsigned r=0; r < repeats; r++)             
                    piz_sam_reconstruct_SA_OA_XA (vb, (dict_id.num == dict_id_OPTION_XA));
            }
        }

        // other optional fields - get each for its own dictionary
        else piz_reconstruct_from_ctx ((VBlockP)vb, did_i, "", 0);        

        if (sf_i != opt_map->num_subfields-1)
            RECONSTRUCT1 ('\t');
    }
    
    if (oname_len % 5 == 1 && oname[oname_len-1] == '#') // Windows style end-of-line \r\n
        RECONSTRUCT1 ('\r');

    // if this line has both BD and BI data, then the BI was calculated as a delta vs. the BD. We apply the delta now
    if (this_bd && this_bi)
        for (uint32_t i=0; i < cigar_seq_len; i++) 
            *(this_bi++) += *(this_bd++);
}

void piz_sam_reconstruct_vb (VBlockSAM *vb)
{
    piz_map_compound_field ((VBlockP)vb, dict_id_is_sam_qname_sf, &vb->qname_mapper);
    piz_sam_map_optional_subfields (vb);

    uint32_t cigar_seq_len=0;
    DECLARE_SNIP;
    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_data_start = vb->txt_data.len;
        vb->line_i = vb->first_line + vb_line_i;

        // QNAME - reconstruct from its subfield components
        LOAD_SNIP (SAM_QNAME);
        piz_reconstruct_compound_field ((VBlockP)vb, &vb->qname_mapper, "\t", 1, snip, snip_len);

        // FLAG, RNAME - from their dictionaries
        piz_reconstruct_from_ctx ((VBlockP)vb, SAM_FLAG, "\t", 1);

        uint32_t rname_word_index = RECONSTRUCT_FROM_DICT (SAM_RNAME, true);

        // POS - reconstruct from pos_data
        // chr5 5000000   rname_minus_2=chr5          rname_minus_3=chr5
        // chr1 1000001   rname_minus_1=chr1          rname_minus_2=chr1
        // chr1 1000002   rname_minus_1=chr1          rname=chr1 <--- DO delta
        // chr1 2000001   rname=chr1 <-- DONT delta
        // chr1 2000002
        if (rname_word_index == vb->rname_index_minus_1 && 
            (rname_word_index != vb->rname_index_minus_2 || rname_word_index == vb->rname_index_minus_3)) {
            RECONSTRUCT_FROM_DICT_POS (SAM_POS, vb->last_pos, true, NULL, true); // same rname - reconstruct from delta 
        }
        else  // different rname - get from random_pos
            vb->last_pos = piz_sam_reconstruct_random_pos (vb, '\t', false);

        vb->rname_index_minus_3 = vb->rname_index_minus_2;
        vb->rname_index_minus_2 = vb->rname_index_minus_1;
        vb->rname_index_minus_1 = rname_word_index;

        // MAPQ - from its dictionary
        piz_reconstruct_from_ctx ((VBlockP)vb, SAM_MAPQ, "\t", 1);

        // CIGAR - get length of SEQ and QUAL, and if original CIGAR was "*" - recover it
        LOAD_SNIP (SAM_CIGAR);
        cigar_seq_len = seg_sam_seq_len_from_cigar (snip, snip_len);

        if (snip[snip_len-1] == '*')   // if its something like "151*" make it "*"
            RECONSTRUCT ("*\t", 2)
        else 
            RECONSTRUCT_TABBED (snip, snip_len); 
        
        
        // RNEXT - from RNAME dictionary
        uint32_t this_rnext_word_index = LOAD_SNIP (SAM_RNAME) // always read, because it consumes the same b250 as RNAME
        RECONSTRUCT_TABBED (snip, snip_len); 
        
        // PNEXT - 2.5 options:
        // 1. if RNEXT is the same as RNAME or unknown (i.e. either it is "=" or "*" or the word index is the same) - 
        //    we decode the delta from the RNEXT b250
        //    1A. If 1, but PNEXT is "*" we recover the "0"
        // 2. If RNEXT and RNAME are different - get the pos data from from vb->random_pos_data 
        bool get_from_dictionary = (this_rnext_word_index == rname_word_index) ||
                                   (snip_len == 1 && snip[0] == '=') || // this refers to the RNEXT snip
                                   (snip_len == 1 && snip[0] == '*');

        if (get_from_dictionary) {
            LOAD_SNIP (SAM_PNEXT);

            if (snip_len == 1 && snip[0] == '*')   // this refers to the PNEXT snip
                RECONSTRUCT ("0\t", 2)
            else 
                RECONSTRUCT_FROM_DICT_POS (DID_I_NONE, vb->last_pos, false, &vb->last_pnext_delta, true);
        }
        else {
            piz_sam_reconstruct_random_pos (vb, '\t', false);
        }
        // TLEN - from its dictionary
        LOAD_SNIP (SAM_TLEN);
        piz_sam_reconstruct_tlen (vb, snip, snip_len, vb->last_pnext_delta, cigar_seq_len);
        
        // SEQ & QUAL data
        piz_reconstruct_seq_qual ((VBlockP)vb, &vb->mtf_ctx[SAM_SEQ], cigar_seq_len, false);
        RECONSTRUCT1 ('\t');

        piz_reconstruct_seq_qual ((VBlockP)vb, &vb->mtf_ctx[SAM_QUAL], cigar_seq_len, false);

        // OPTIONAL fields, and Windows-style \r if needed
        uint32_t word_index = LOAD_SNIP (SAM_OPTIONAL);
        if (snip_len != 1 || (snip[0] != '#' && snip[0] != '*')) RECONSTRUCT1 ('\t');
        
        piz_sam_reconstruct_optional_fields (vb, cigar_seq_len, snip, snip_len, word_index);

        RECONSTRUCT1 ('\n');

        // after consuming sections' data, if this line is not to be outputed - shorten txt_data back to start of line
        if (flag_regions && !regions_is_site_included (vb->rname_index_minus_1, vb->last_pos))
            vb->txt_data.len = txt_data_start; // remove excluded line
    }
}
