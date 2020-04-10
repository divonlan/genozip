// ------------------------------------------------------------------
//   v2v3.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// in genozip v2 and v3 we had a bug where when there was an INFO subfield with a value following one or more
// subfields without a value, then the dictionary name of included the names of the valueless subfields
//
// Consider the file:
// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
// 13      2072    rs2     T       T       100     PASS    I2=a;N1;N2;I3=x
//
// --show-sections shows the dictionaries as:
// I2       FIELD             1            1          509       2 B
// N1;N2;I3 INFO              1            1          509       2 B
//
// We fixed this bug in 4.0.0 incrementing the major version to reflect the resulting change
// in the file format, and included this code below for decompressing files that were compressed with the
// buggy v2/v3 format
//-----------------------------------------------------------------------------------------------------------

// for each unique type of INFO fields (each one containing multiple names), create a unique mapping
// info field node index (i.e. in b250) -> list of names, lengths and the context of the subfields
static void v2v3_piz_vcf_map_iname_subfields (VBlockVCF *vb)
{
    // terminology: we call a list of INFO subfield names, an "iname". An iname looks something like
    // this: "I1=I2=I3=". Each iname consists of info subfields. These fields are not unique to this
    // iname and can appear in other inames. The INFO field contains the iname, and values of the subfields.
    // iname _mapper maps these subfields. This function creates an iname_mapper for every unique iname.

    const MtfContext *info_ctx = &vb->mtf_ctx[VCF_INFO];
    vb->iname_mapper_buf.len = info_ctx->word_list.len;
    buf_alloc (vb, &vb->iname_mapper_buf, sizeof (SubfieldMapper) * vb->iname_mapper_buf.len,
               1, "iname_mapper_buf", 0);
    buf_zero (&vb->iname_mapper_buf);

    SubfieldMapper *all_iname_mappers = (SubfieldMapper*)vb->iname_mapper_buf.data;

    const MtfWord *all_inames = (const MtfWord *)info_ctx->word_list.data;

    for (unsigned iname_i=0; iname_i < vb->iname_mapper_buf.len; iname_i++) {

        const char *iname = (const char *)&info_ctx->dict.data[all_inames[iname_i].char_index]; // e.g. "I1=I2=I3=" - pointer into the INFO dictionary
        unsigned iname_len = all_inames[iname_i].snip_len; 
        SubfieldMapper *iname_mapper = &all_iname_mappers[iname_i]; // iname_mapper of this specific set of names "I1=I2=I3="

        // get INFO subfield snips - which are the values of the INFO subfield, where the names are
        // in the INFO snip in the format "info1=info2=info3="). 
        DictIdType dict_id;
        iname_mapper->num_subfields = 0;

        // traverse the subfields of one iname. E.g. if the iname is "I1=I2=I3=" then we traverse I1, I2, I3
        for (unsigned i=0; i < iname_len; i++) {
            
            // traverse the iname, and get the dict_id for each subfield name (using only the first 8 characers)
            dict_id.num = 0;
            unsigned j=0; 
            while (iname[i] != '=' && iname[i] != '\t') { // value-less INFO names can be terminated by the end-of-word \t in the dictionary
                if (j < DICT_ID_LEN) 
                    dict_id.id[j] = iname[i]; // scan the whole name, but copy only the first 8 bytes to dict_id
                i++, j++;
            }
            dict_id = dict_id_vcf_info_sf (dict_id);

            iname_mapper->did_i[iname_mapper->num_subfields] = mtf_get_existing_did_i_by_dict_id ((VBlockP)vb, dict_id); // it will be NIL if this is an INFO name without values            
            iname_mapper->num_subfields++;
        }
    }
}

static inline void v2v3_piz_vcf_reconstruct_info (VBlockVCF *vb, const SubfieldMapper *iname_mapper, const char *snip,
                                              const char **info_sf_value_snip, uint32_t *info_sf_value_snip_len, bool dummy)
{
    for (unsigned sf_i=0; sf_i < iname_mapper->num_subfields ; sf_i++) {

        // get the name eg "AF="
        const char *start = snip;
        for (; *snip != '=' && *snip != '\t'; snip++);
        if (*snip == '=') snip++; // move past the '='

        buf_add (&vb->line_variant_data, start, (unsigned)(snip-start)); // name inc. '=' e.g. "Info1="

        if (iname_mapper->did_i[sf_i] != DID_I_NONE)  // some info names can be without values, in which case there will be no ctx
            buf_add (&vb->line_variant_data, info_sf_value_snip[sf_i], info_sf_value_snip_len[sf_i]); // value e.g "value1"

        if (sf_i != iname_mapper->num_subfields-1)
            buf_add (&vb->line_variant_data, ";", 1); // seperator between each two name=value pairs e.g "name1=value;name2=value2"
    }
}