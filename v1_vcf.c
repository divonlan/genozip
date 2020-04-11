// ------------------------------------------------------------------
//   v1.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// read section header - called from the I/O thread, but for a specific VB
// returns offset of header within data, EOF if end of file

// NOTE: this file contains all the functions that have a different version for v1 vs v2 genozip.
// it is #included at the end of the affected source files

#ifdef V1_ZFILE

// note: the first section (i.e. the first VCF header) is always read by v1_zfile_read_section() (the current version)
// if header.genozip_version is 1, then subsequent sections, including subsequent VCF headers, will be read by this function
int v1_zfile_read_section (VBlock *vb,
                           Buffer *data, const char *buf_name, /* buffer to append */
                           unsigned header_size, SectionType expected_sec_type,
                           bool allow_eof)
{
    bool is_encrypted = crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    unsigned header_offset = expected_sec_type == SEC_TXT_HEADER ? 0 : data->len;
    
    unsigned requested_bytes = header_size;
    if (expected_sec_type == SEC_TXT_HEADER) requested_bytes -= data->len;

    buf_alloc (vb, data, header_offset + header_size, 2, buf_name, 1);

    SectionHeader *header = zfile_read_from_disk (vb, data, requested_bytes, false); // note: header in file can be shorter than header_size if its an earlier version
    ASSERT (header || allow_eof, "Error: genozip v1 file - Failed to read header, section_type=%s", st_name(expected_sec_type));
    
    // if this is the first VCF header, part of it was already read by zfile_read_genozip_header() and placed in data (which is z_file->v1_next_vcf_header)
    if (expected_sec_type == SEC_TXT_HEADER) header = (SectionHeader *)data->data;

    if (!header) return EOF; 

    // decrypt header
    if (is_encrypted) {
        ASSERT (BGEN32 (header->magic) != GENOZIP_MAGIC, 
                "Error: genozip v1 file - password provided, but file %s is not encrypted", z_name);

        v1_crypt_do (vb, (uint8_t*)header, header_size, vb->vblock_i, --vb->z_next_header_i); // negative section_i for a header
    }
    bool is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;
    if (!is_magical && !is_encrypted && expected_sec_type == SEC_TXT_HEADER) {

        // file appears to be encrypted but user hasn't provided a password    
        crypt_prompt_for_password();

        unsigned padding;
        is_encrypted = crypt_get_encrypted_len (&header_size, &padding); // update header size if encrypted
            
        if (padding) {
            char *header_extra_bytes = zfile_read_from_disk (vb, data, padding, false);
            ASSERT0 (header_extra_bytes, "Error: genozip v1 file - Failed to read header padding");
        }

        v1_crypt_do (vb, (uint8_t*)header, header_size, vb->vblock_i, --vb->z_next_header_i);
        is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;
    }

    if (!is_magical && is_encrypted && expected_sec_type == SEC_TXT_HEADER) {
        ABORT ("Error: genozip v1 file - password is wrong for file %s", z_name); // mostly likely its because of a wrong password
    }

    // case: encryption failed because this is actually a SEC_TXT_HEADER of a new vcf component of a concatenated file
    bool new_vcf_component = false;
    unsigned new_header_size=0;

    if (!is_magical && is_encrypted && expected_sec_type == SEC_VCF_VB_HEADER) {
    
        // reverse failed decryption
        v1_crypt_do (vb, (uint8_t*)header, header_size, vb->vblock_i, vb->z_next_header_i);

        new_header_size = sizeof (SectionHeaderTxtHeader);
        unsigned padding;
        crypt_get_encrypted_len (&new_header_size, &padding); // adjust header size with encryption block size

        // read additional bytes - if we need to
        int additional_bytes = new_header_size - header_size;

        bool success;
        if (additional_bytes > 0) {
            buf_alloc (vb, data, data->len + additional_bytes, 1, 0, 0);
            header = (SectionHeader *)&data->data[header_offset]; // update after realloc
            success = (bool)zfile_read_from_disk (vb, data, new_header_size - header_size, true);
        }
        else 
            success = true; // we don't need any extra bytes

        if (success) { // success
            // attempt to re-decrypt, with a key for the vcf header
            v1_crypt_do (vb, (uint8_t*)header, new_header_size, 0, -1); // vb_i=0 and sec_i=-1 always, for all VCFHeader section headers
            is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;

            vb->z_next_header_i++; // roll back
            new_vcf_component = true; 
        }
    } 

    ASSERT (is_magical, "Error: genozip v1 file - corrupt data (magic is wrong) when reading file %s", z_name);

    unsigned compressed_offset   = BGEN32 (header->compressed_offset);
    ASSERT (compressed_offset, "Error: genozip v1 file - header.compressed_offset is 0 when reading section_type=%s", st_name(expected_sec_type));

    unsigned data_compressed_len = BGEN32 (header->data_compressed_len);
    ASSERT (data_compressed_len, "Error: genozip v1 file - header.data_compressed_len is 0 when reading section_type=%s", st_name(expected_sec_type));

    unsigned data_encrypted_len  = BGEN32 (header->data_encrypted_len);

    unsigned data_len = MAX (data_compressed_len, data_encrypted_len);

    // We found a VCF header - possibly a result of concatenating files:
    // if regular mode (no split) - we will just read the data from disk to skip this header (unless its expected)
    // if in split mode - we will end processing this output file here and store the header for the next output vcf file
    bool found_a_vcf_header = header->section_type == SEC_TXT_HEADER;
    
    // check that we received the section type we expect, 
    ASSERT (header->section_type == expected_sec_type || (found_a_vcf_header && expected_sec_type == SEC_VCF_VB_HEADER),
            "Error: genozip v1 file - unexpected section type: expecting %s, found %s", st_name(expected_sec_type), st_name(header->section_type));

    unsigned expected_header_size = new_vcf_component ? new_header_size : header_size;

    // case: this is actually not a v1 file, but it got truncated so the v2 reader didn't find the genozip header section
    // and referred it to v1
    ASSERT (!(header->section_type == SEC_TXT_HEADER && compressed_offset == sizeof(SectionHeaderTxtHeader)), 
            "Error: failed to read file %s - it appears to be truncated or corrupted", z_name);

    ASSERT (compressed_offset == crypt_padded_len (expected_header_size) || expected_sec_type == SEC_VCF_VB_HEADER, // for variant data, we also have the permutation index
            "Error: genozip v1 file - invalid header - expecting compressed_offset to be %u but found %u", expected_header_size, compressed_offset);

    // allocate more memory for the rest of the header + data (note: after this realloc, header pointer is no longer valid)
    buf_alloc (vb, data, header_offset + compressed_offset + data_len, 2, "v1_zfile_read_section", 2);
    header = (SectionHeader *)&data->data[header_offset]; // update after realloc

    // in case we're expecting SEC_VCF_VB_HEADER - read the rest of the header: 
    // if we found SEC_VCF_VB_HEADER, we need to read haplotype index, and if we found a SEC_TXT_HEADER of a concatenated
    // file componented - we need to read the read of the header as it is bigger than Variant Data header that was read
    if (expected_sec_type == SEC_VCF_VB_HEADER && !new_vcf_component) {

        int bytes_left_over = compressed_offset - header_size;
        ASSERT (bytes_left_over >= 0, "Error: genozip v1 file - expected bytes_left_over=%d to be >=0", bytes_left_over)

        if (bytes_left_over) { // there will be an Index only if this VCF has samples
            
            uint8_t *left_over_data = zfile_read_from_disk (vb, data, bytes_left_over, false);
            ASSERT (left_over_data, "Error: genozip v1 file - Failed to read left over bytes. bytes_left_over=%u", bytes_left_over);

            if (is_encrypted) { // this is just for the ht index - we've already handle encrypted vcf header and set new_vcf_component=true 
                ASSERT (bytes_left_over == crypt_padded_len (bytes_left_over), "Error: genozip v1 file - bad length of bytes_left_over=%u", bytes_left_over); // we expected it to be aligned, bc the total encryption block is aligned and header_size is aligned 
                
                // for the haplotype index - it is part of the header so we just continue the encryption stream
                crypt_continue (vb, left_over_data, bytes_left_over);
            }
            if (header->section_type == SEC_TXT_HEADER) 
                new_vcf_component = true;
        }
    }

    // read section data
    if (data_len) {
        ASSERT (zfile_read_from_disk (vb, data, data_len, false), 
                "Error: genozip v1 file - failed to read section data, section_type=%s", st_name(header->section_type));
    }

    // deal with a VCF header that was encountered while expecting a SEC_VCF_VB_HEADER (i.e. 2nd+ component of a concatenated file)
    if (new_vcf_component) {

        if (flag_split) {
            // move the header from vb->z_data to z_file->v1_next_vcf_header
            buf_copy (vb, &z_file->v1_next_vcf_header, data, 1, header_offset, data->len - header_offset, "z_file->v1_next_vcf_header", 0);
            data->len = header_offset; // shrink - remove the vcf header that doesn't belong to this component

            return EOF; // end of VCF component in flag_split mode
        }
        else {
            // since we're not in split mode, so we can skip this vcf header section - just return the next section instead
            return v1_zfile_read_section (vb, data, buf_name, header_size, expected_sec_type, allow_eof);
        }
    }

    return header_offset;
}


bool v1_zfile_vcf_read_one_vb (VBlockVCF *vb)
{ 
    START_TIMER;

    int vardata_header_offset = v1_zfile_read_section ((VBlockP)vb, &vb->z_data, "z_data",
                                                       sizeof(v1_SectionHeaderVariantData), SEC_VCF_VB_HEADER, true);
    if (vardata_header_offset == EOF) {

        // update size - in case they were not known (pipe, gzip etc)
        if (!z_file->disk_size) 
            z_file->disk_size = z_file->disk_so_far;
            
        COPY_TIMER (vb->profile.zfile_read_one_vb);
        return false; // end of file
    }

    // note - copying values here z_data.data can get reallocated each call to v1_zfile_read_section
    v1_SectionHeaderVariantData *vardata_header = (v1_SectionHeaderVariantData *)&vb->z_data.data[vardata_header_offset];
    unsigned num_sample_blocks       = BGEN32 (vardata_header->num_sample_blocks);
    bool has_genotype_data           = vardata_header->has_genotype_data;
    PhaseType phase_type             = (PhaseType)vardata_header->phase_type;
    unsigned num_haplotypes_per_line = BGEN32 (vardata_header->num_haplotypes_per_line);
    unsigned num_dictionary_sections = BGEN16 (vardata_header->num_dictionary_sections);

    // dictionaries are processed right here by the dispatcher thread - the compute
    // thread only access the dictionaries on the z_file->mtf_ctx
    if (num_dictionary_sections) {

        unsigned start_dictionary_sections = vb->z_data.len;

        // read all sections into memory
        for (unsigned did_i=0; did_i < num_dictionary_sections; did_i++) {

            unsigned start_i = vb->z_data.len; // vb->z_data.len is updated next, by v1_zfile_read_section()
            v1_zfile_read_section ((VBlockP)vb, &vb->z_data, "z_data", sizeof(SectionHeaderDictionary), SEC_VCF_FRMT_SF_DICT, false);    

            // update dictionaries in z_file->mtf_ctx with dictionary data from this VB
            mtf_integrate_dictionary_fragment ((VBlockP)vb, &vb->z_data.data[start_i]);
        }

        vb->z_data.len = start_dictionary_sections; // shrink z_data back
    }

    // overlay all available dictionaries (not just those that have fragments in this variant block) to the vb
    mtf_overlay_dictionaries_to_vb ((VBlockP)vb);

    // read the other sections

    buf_alloc (vb, &vb->z_section_headers, 1000 /* arbitrary initial value */ * sizeof(char*), 0, "z_section_headers", 1);
    
    ((unsigned *)vb->z_section_headers.data)[0] = vardata_header_offset; // variant data header is at index 0

    // is_sb_included was introduced in version 4.0.x, we just set it "all included" for v1
    buf_alloc (vb, &vb->is_sb_included, num_sample_blocks * sizeof(bool), 1, "is_sb_included", vb->vblock_i);

    unsigned section_i=1;

    for (unsigned sb_i=0; sb_i < num_sample_blocks; sb_i++) {

        // make sure we have enough space for the section pointers
        buf_alloc (vb, &vb->z_section_headers, sizeof (unsigned) * (section_i + 3), 2, "z_section_headers", 2);

        if (has_genotype_data) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            v1_zfile_read_section ((VBlockP)vb, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_VCF_GT_DATA, false);
        }

        if (phase_type == PHASE_MIXED_PHASED) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            v1_zfile_read_section ((VBlockP)vb, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_VCF_PHASE_DATA, false);
        }

        if (num_haplotypes_per_line) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            v1_zfile_read_section ((VBlockP)vb, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_VCF_HT_DATA , false);    
        }

        *ENT (bool, &vb->is_sb_included, sb_i) = true;
    }
    
    COPY_TIMER (vb->profile.zfile_read_one_vb);

    return true; 
}

#endif // V1_ZFILE

#ifdef V1_PIZ

#define V1_SAMPLES_PER_BLOCK 4096

// decode the delta-encoded value of the POS field
static inline void v1_piz_vcf_decode_pos (VBlockVCF *vb, const char *str,
                                      char *pos_str, const char **delta_pos_start, unsigned *delta_pos_len, int *add_len /* out */)
{
    START_TIMER;

    *delta_pos_start = str;

    int32_t delta=0;

    // we parse the string ourselves - this is hugely faster than sscanf.
    unsigned negative = *str == '-'; // 1 or 0

    const char *s; for (s=(str + negative); *s != '\t' ; s++)
        delta = delta*10 + (*s - '0');

    if (negative) delta = -delta;

    vb->last_pos += delta;
    
    ASSERT (vb->last_pos >= 0, "Error: vb->last_pos=%d is negative", vb->last_pos);

    // create number string without calling slow sprintf

    // create reverse string
    char reverse_pos_str[50];
    uint32_t n = vb->last_pos;
    unsigned len=0; 
    while (n) {
        reverse_pos_str[len++] = '0' + (n % 10);
        n /= 10;
    }

    // reverse it
    for (unsigned i=0; i < len; i++) pos_str[i] = reverse_pos_str[len-i-1];
    pos_str[len] = '\0';

    *delta_pos_len = s - str;
    *add_len = len - *delta_pos_len;

    COPY_TIMER(vb->profile.piz_vcf_decode_pos);
}

static void v1_piz_get_line_get_num_subfields (VBlockVCF *vb, unsigned line_i, // line in vcf file
                                               const char **line, unsigned *remaining_len,
                                               const char **subfields_start, unsigned *subfields_len, int *num_subfields)
{    
    START_TIMER;

    const char *after = *line + *remaining_len;

    unsigned column=1, i=0; for (; i < *remaining_len && column < 9; i++)
        if      ((*line)[i] == '\t') column++;
        else if ((*line)[i] == '\n') break;

    ASSERT (column==9, "Error: line %u is missing a FORMAT field", line_i); 

    *line += i; // past the tab
    
    if ((*line)[0] == '\n') { // this row has no genotype data, but maybe elsewhere in the VB there is
        *num_subfields = (unsigned)vb->has_genotype_data; // if we have genotype, then zip inserted BASE250_MISSING_SF 
        (*line)++; // past the newline
        goto cleanup;
    }
    
    // Check for GT field in FORMAT columns - must always appear first per VCF spec (if at appears)
    if ((*line)[0] == 'G' && (*line)[1] == 'T' && ((*line)[2] == ':' || (*line)[2] == '\n')) { 
        *line += 3; // past the GT and : or \n
        if ((*line)[-1] == '\n') { // no subfields in this line
            *num_subfields = (unsigned)vb->has_genotype_data; // if we have genotype, then zip inserted BASE250_MISSING_SF 
            goto cleanup;
        }
    }

    *subfields_start = *line;

    // count the number of subfields in the line
    for (*num_subfields = 1; *num_subfields <= MAX_SUBFIELDS; (*num_subfields)++) {
        uint32_t len = after - *line;
        seg_get_format_subfield (line, &len, line_i);
        if ((*line)[-1] == '\n') break;
    } 

    *subfields_len = (unsigned)(*line - *subfields_start); // length including separator

    ASSERT ((*line)[-1] == '\n', "Error: number of subfields declared in FORMAT exceeds maximum allowed of %u (excluding GT), line=%u", MAX_SUBFIELDS, vb->first_line+line_i);

cleanup:
    *remaining_len = after - *line;

    COPY_TIMER (vb->profile.piz_vcf_get_format_info)
}

static void v1_piz_vcf_get_variant_data_line (VBlockVCF *vb, 
                                          unsigned line_i, // line in vcf file
                                          unsigned *length_remaining, // for safety
                                          const char **line_start) // out
{
    START_TIMER;

    // decoding the delta-encouded value of the POS field
    const char *delta_pos_start=0; // start of delta-POS field
    unsigned delta_pos_len=0; // length of encoded delta-pos
    int add_len=0; // how much more (or less) characters do we need in the string after decoding POS?
    char pos_str[22]; // decoded POS value

    const char *next = *line_start;
    const char *after = *line_start + *length_remaining;
    unsigned column = 1;

    while (next < after) {

        // starting a new field
        if (*next == '\t') {
            column++;

            // if we're at the POS field, decode the delta encoding
            if (column == 2)  
                v1_piz_vcf_decode_pos (vb, next+1, pos_str, &delta_pos_start, &delta_pos_len, &add_len);
        }

        else if (*next == '\n') {
            next++; // past the newline
            unsigned line_strlen = next - *line_start; // inc. the newline, not including \0

            buf_alloc (vb, &vb->line_variant_data, line_strlen + add_len + 1 /* +1 sprintf adds a \0 */, 
                       1.2, "line_variant_data", 0);

            const char *after_delta = delta_pos_start + delta_pos_len;

            sprintf (vb->line_variant_data.data, "%.*s%s%.*s",
                     (int)(delta_pos_start - *line_start), *line_start,    // substring until \t preceding delta
                     pos_str,                           // decoded pos string
                     (int)(next - after_delta), after_delta);  // substring starting \t after delta
            
            vb->line_variant_data.len = line_strlen + add_len;

            *line_start = next;
            *length_remaining = after - next;

            goto cleanup;
        }

        next++;
    }
    
    ABORT0 ("Error: corrupt genozip file - at end of variant_data buffer, and no newline was found");

cleanup:
    COPY_TIMER(vb->profile.piz_vcf_get_variant_data_line);
    return;
}

static void v1_piz_get_line_subfields (VBlockVCF *vb, unsigned line_i, // line in vcf file
                                       const char *subfields_start, unsigned subfields_len,
                                       int *line_subfields) // out
{
    START_TIMER;

    for (unsigned i=0; i < MAX_SUBFIELDS; i++) line_subfields[i] = NIL;

    // case: this line has no subfields, despite other lines in the VB having
    if (!subfields_len) {
        line_subfields[0] = -2;
        return;
    }

    const char *after = subfields_start + subfields_len;
    for (unsigned i=0; i < MAX_SUBFIELDS; i++) {
        uint32_t len = after - subfields_start;
        DictIdType subfield = seg_get_format_subfield (&subfields_start, &len, line_i);

        // the dictionaries were already read, so all subfields are expected to have a ctx
        unsigned did_i=0 ; for (; did_i < z_file->num_dict_ids; did_i++) 
            if (subfield.num == z_file->mtf_ctx[did_i].dict_id.num) {
                // entry i corresponds to subfield i in FORMAT (excluding GT), and contains the index in mtf_ctx of this subfield
                line_subfields[i] = did_i;
                break;
            }
#ifdef DEBUG
        ASSERTW (did_i < z_file->num_dict_ids, 
                 "Warning: subfield %.*s not found in dictionaries, line=%u. This can happen legitimately if the subfield is declared in FORMAT, but a value is never provided in any sample", DICT_ID_LEN, subfield.id, line_i);
#endif
        if (subfields_start[-1] == '\t' || subfields_start[-1] == '\n') break;
    } 

    COPY_TIMER (vb->profile.piz_get_line_subfields)
}

// convert genotype data from sample block format of indices in base-250 to line format
// of tab-separated genotypes
static void v1_piz_vcf_get_genotype_data_line (VBlockVCF *vb, unsigned line_i, int *line_subfields)
{
    START_TIMER;

    SnipIterator *sample_iterator = (SnipIterator *)vb->sample_iterator.data; // for convenience

    char *next = vb->line_gt_data.data;
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned first_sample = sb_i*V1_SAMPLES_PER_BLOCK;
        unsigned num_samples_in_sb = vb_vcf_num_samples_in_sb (vb, sb_i);

        for (unsigned sample_i=first_sample; 
             sample_i < first_sample + num_samples_in_sb; 
             sample_i++) {

            const char *snip = NULL; // will be set to a pointer into a dictionary

            // case: this line has no non-GT subfields defined in format, despite other lines in the VB having
            // we have theremore, only one "fake" subfield in this line, which we shall skip
            if (line_subfields[0] == -2) {
                ASSERT (*sample_iterator[sample_i].next_b250 == BASE250_MISSING_SF, 
                        "Error in vcf line %u - line has no subfields, expecting BASE250_MISSING_SF but not seeting it", vb->first_line + line_i);
                
                sample_iterator[sample_i].next_b250++; // skip
            }
            else {
                for (uint8_t sf_i=0; sf_i < vb->num_format_subfields; sf_i++) {

                    if (line_subfields[sf_i] != NIL) {  // this line has this subfield (according to its FORMAT field)

                        // add a colon before, if needed
                        if (snip) *(next++) = ':'; // this works for empty "" snip too

                        unsigned snip_len;
                        mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[line_subfields[sf_i]], // note: line_subfields[sf_i] maybe -2 (set in piz_get_line_subfields()), and this is an invalid value. this is ok, bc in this case sample_iterator[sample_i] will be a control character
                                        &sample_iterator[sample_i], &snip, &snip_len, vb->first_line + line_i);

                        if (snip && snip_len) { // it can be a valid empty subfield if snip="" and snip_len=0
                            memcpy (next, snip, snip_len);
                            next += snip_len;
                        }
                    }
                }
            }

            // if we ended with a : - remove it
            next -= (next[-1] == ':');

            // add sample terminator - \t
            *(next++) = '\t';

            // safety
            ASSERT (next <= vb->line_gt_data.data + vb->line_gt_data.size, 
                    "Error: line_gt_data buffer overflow. vblock_i=%u vcf_line=%u sb_i=%u sample_i=%u",
                    vb->vblock_i, line_i + vb->first_line, sb_i, sample_i);
        }
    }

    // change last terminator to a \n
    next[-1] = '\n';

    vb->line_gt_data.len = next - vb->line_gt_data.data;

    DATA_LINE(vb, line_i)->has_genotype_data = vb->line_gt_data.len > global_vcf_num_samples; // not all just \t

    COPY_TIMER(vb->profile.piz_vcf_get_genotype_data_line);
}

static void v1_piz_initialize_next_gt_in_sample (VBlockVCF *vb, int *num_subfields)
{
    START_TIMER;
    
    buf_alloc (vb, &vb->sample_iterator, sizeof(SnipIterator) * global_vcf_num_samples, 1, "sample_iterator", 0);
    SnipIterator *sample_iterator = (SnipIterator *)vb->sample_iterator.data; 
    
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sb = vb_vcf_num_samples_in_sb (vb, sb_i);
        
        const uint8_t *next = (const uint8_t *)vb->genotype_sections_data[sb_i].data;
        const uint8_t *after = next + vb->genotype_sections_data[sb_i].len;

        unsigned sample_after = sb_i * V1_SAMPLES_PER_BLOCK + num_samples_in_sb;
        
        unsigned sample_i = sb_i * V1_SAMPLES_PER_BLOCK; 
        for (;sample_i < sample_after && next < after; sample_i++) {
            
            sample_iterator[sample_i].next_b250 = next; // line=0 of each sample_i (column)
            sample_iterator[sample_i].prev_word_index = -1;

            // now skip all remaining genotypes in this column, arriving at the beginning of the next column
            // (gt data is stored transposed - i.e. column by column)
            for (unsigned line_i=0; line_i < vb->num_lines; line_i++)
                for (unsigned sf=0; sf < num_subfields[line_i]; sf++) 
                    next += v1_base250_len (next); 
        }

        // sanity checks to see we read the correct amount of genotypes
        ASSERT (sample_i == sample_after, "Error: expected to find %u genotypes in sb_i=%u of vblock_i=%u, but found only %u",
                vb->num_lines * num_samples_in_sb, sb_i, vb->vblock_i, vb->num_lines * (sample_i - sb_i * V1_SAMPLES_PER_BLOCK));

        ASSERT (next == after, "Error: expected to find %u genotypes in sb_i=%u of vblock_i=%u, but found more. ",
                vb->num_lines * num_samples_in_sb, sb_i, vb->vblock_i);
    }

    COPY_TIMER (vb->profile.piz_vcf_initialize_sample_iterators)
}

// combine all the sections of a variant block to regenerate the variant_data, haplotype_data,
// genotype_data and phase_data for each row of the variant block
void v1_piz_vcf_reconstruct_line_components (VBlockVCF *vb)
{
    START_TIMER;

    if (!vb->data_lines) 
        vb->data_lines = calloc (txt_file->max_lines_per_vb, sizeof (PizDataLineVCF));

    // initialize phase data if needed
    if (vb->phase_type == PHASE_MIXED_PHASED) 
        buf_alloc (vb, &vb->line_phase_data, global_vcf_num_samples, 1, "line_phase_data", 0);

    // initialize haplotype stuff
    const char **ht_columns_data=NULL;
    if (vb->has_haplotype_data) {

        //  memory - realloc for exact size, add 7 because depermuting_loop works on a word (32/64 bit) boundary
        buf_alloc (vb, &vb->line_ht_data, vb->num_haplotypes_per_line + 7, 1, "line_ht_data", 0);

        ht_columns_data = piz_vcf_get_ht_columns_data (vb);
    }

    // traverse the variant data first, only processing the FORMAT field - populate
    // the subfield data needed by v1_piz_initialize_next_gt_in_sample

    const char *variant_data_next_line = vb->v1_variant_data_section_data.data;
    unsigned variant_data_length_remaining = vb->v1_variant_data_section_data.len;
    
    buf_alloc (vb, &vb->v1_subfields_start_buf, vb->num_lines * sizeof (char *),   1, "v1_subfields_start_buf", 0);
    buf_zero (&vb->v1_subfields_start_buf);

    buf_alloc (vb, &vb->v1_subfields_len_buf,   vb->num_lines * sizeof (unsigned), 1, "v1_subfields_len_buf", 0);
    buf_zero (&vb->v1_subfields_len_buf);

    buf_alloc (vb, &vb->v1_num_subfields_buf,   vb->num_lines * sizeof (int),      1, "v1_num_subfields_buf", 0);
    buf_zero (&vb->v1_num_subfields_buf);
    
    const char **subfields_start = (const char **) vb->v1_subfields_start_buf.data; // pointer within the FORMAT field
    unsigned *subfields_len      = (unsigned *)    vb->v1_subfields_len_buf.data;   // length of the FORMAT field, excluding GT, including the separator
    int *num_subfields           = (int *)         vb->v1_num_subfields_buf.data;   // number of subfields excluding GT
        
    // initialize genotype stuff
    if (vb->has_genotype_data) {
        for (unsigned line_i=0; line_i < vb->num_lines; line_i++) 
            // get subfields info from the FORMAT field
            v1_piz_get_line_get_num_subfields (vb, vb->first_line + line_i, 
                                            &variant_data_next_line, &variant_data_length_remaining,
                                            &subfields_start[line_i], &subfields_len[line_i], &num_subfields[line_i]);

        v1_piz_initialize_next_gt_in_sample(vb, num_subfields);

        buf_alloc (vb, &vb->line_gt_data, vb->max_gt_line_len, 1, "line_gt_data", 0);
    }

    // initialize again - for piz_vcf_get_variant_data_line
    variant_data_next_line = vb->v1_variant_data_section_data.data;
    variant_data_length_remaining = vb->v1_variant_data_section_data.len;

    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {

        // de-permute variant data into vb->line_variant_data
        v1_piz_vcf_get_variant_data_line (vb, vb->first_line + line_i, &variant_data_length_remaining, &variant_data_next_line);

        // reset len for next line - no need to realloc as we have realloced what is needed already
        vb->line_ht_data.len = vb->line_gt_data.len = vb->line_phase_data.len = 0;

        // transform sample blocks (each block: n_lines x s_samples) into line components (each line: 1 line x ALL_samples)
        if (vb->has_genotype_data)  {
            int line_subfields[MAX_SUBFIELDS]; // entry i corresponds to subfield i in FORMAT (excluding GT), and contains the index in mtf_ctx of this subfield
            v1_piz_get_line_subfields (vb, vb->first_line + line_i,
                                       subfields_start[line_i], subfields_len[line_i], line_subfields);

            v1_piz_vcf_get_genotype_data_line (vb, line_i, line_subfields);
        }
        if (vb->phase_type == PHASE_MIXED_PHASED) 
            piz_vcf_get_phase_data_line (vb, line_i);

        if (vb->has_haplotype_data) 
            piz_vcf_get_haplotype_data_line (vb, line_i, ht_columns_data);

        piz_vcf_merge_line (vb, line_i, false);
    }

    COPY_TIMER(vb->profile.piz_vcf_reconstruct_line_components);
}

void v1_piz_vcf_uncompress_all_sections (VBlockVCF *vb)
{
    unsigned *section_index = (unsigned *)vb->z_section_headers.data;

    // get the variant data - newline-seperated lines, each containing the first 8 (if no FORMAT field) or 9 fields (if FORMAT exists)
    zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[0], 
                              &vb->v1_variant_data_section_data, "v1_variant_data_section_data", SEC_VCF_VB_HEADER);
    
    v1_SectionHeaderVariantData *vardata_header = (v1_SectionHeaderVariantData *)(vb->z_data.data + section_index[0]);
    vb->first_line              = BGEN32 (vardata_header->first_line);
    vb->num_lines               = BGEN32 (vardata_header->num_lines);
    vb->phase_type              = (PhaseType)vardata_header->phase_type;
    vb->has_genotype_data       = vardata_header->has_genotype_data;
    vb->num_haplotypes_per_line = BGEN32 (vardata_header->num_haplotypes_per_line);
    vb->has_haplotype_data      = vb->num_haplotypes_per_line > 0;
    vb->num_sample_blocks       = BGEN32 (vardata_header->num_sample_blocks);
    vb->num_samples_per_block   = BGEN32 (vardata_header->num_samples_per_block);
    vb->ploidy                  = BGEN16 (vardata_header->ploidy);
    vb->num_dict_ids            = BGEN16 (vardata_header->num_dict_ids);
    // num_dictionary_sections is read in zfile_vcf_read_one_vb()
    vb->max_gt_line_len         = BGEN32 (vardata_header->max_gt_line_len);
    vb->vb_data_size            = BGEN32 (vardata_header->vb_data_size);
    
    ASSERT (global_vcf_num_samples == BGEN32 (vardata_header->num_samples), "Error: Expecting variant block to have %u samples, but it has %u", global_vcf_num_samples, BGEN32 (vardata_header->num_samples));

    // we allocate memory for the Buffer arrays only once the first time this VBlockVCF
    // is used. Subsequent blocks reusing the memory will have the same number of samples (by VCF spec)
    // BUG: this won't work if we're doing mutiple unrelated VCF on the command line
    if (vb->has_genotype_data && !vb->genotype_sections_data) 
        vb->genotype_sections_data  = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));

    if (vb->phase_type == PHASE_MIXED_PHASED && !vb->phase_sections_data) 
        vb->phase_sections_data     = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));
    
    if (vb->num_haplotypes_per_line && !vb->haplotype_sections_data) 
        vb->haplotype_sections_data = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));

    // unsqueeze permutation index - if this VCF has samples
    if (global_vcf_num_samples) {
        buf_alloc (vb, &vb->haplotype_permutation_index, vb->num_haplotypes_per_line * sizeof(unsigned), 0, 
                    "haplotype_permutation_index", vb->first_line);

        unsqueeze (vb,
                   (unsigned *)vb->haplotype_permutation_index.data, 
                   vardata_header->haplotype_index, 
                   BGEN16 (vardata_header->haplotype_index_checksum),
                   vb->num_haplotypes_per_line);
    }

  // get data for sample blocks - each block *may* have up to 3 file sections - genotype, phase and haplotype

    unsigned section_i=1;

    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sb = vb_vcf_num_samples_in_sb (vb, sb_i);

        // if genotype data exists, it appears first
        if (vb->has_genotype_data) 
            zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i++], 
                                      &vb->genotype_sections_data[sb_i], "genotype_sections_data", SEC_VCF_GT_DATA);
        
        // next, comes phase data
        if (vb->phase_type == PHASE_MIXED_PHASED) {
            
            zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i++], 
                                      &vb->phase_sections_data[sb_i], "phase_sections_data", SEC_VCF_PHASE_DATA);
            
            unsigned expected_size = vb->num_lines * num_samples_in_sb;
            ASSERT (vb->phase_sections_data[sb_i].len==expected_size, 
                    "Error: unexpected size of phase_sections_data[%u]: expecting %u but got %u", sb_i, expected_size, vb->phase_sections_data[sb_i].len)
        }

        // finally, comes haplotype data
        if (vb->has_haplotype_data) {
            
            zfile_uncompress_section ((VBlockP)vb, vb->z_data.data + section_index[section_i++], 
                                      &vb->haplotype_sections_data[sb_i], "haplotype_sections_data", SEC_VCF_HT_DATA );
            
            unsigned expected_size = vb->num_lines * num_samples_in_sb * vb->ploidy;
            ASSERT (vb->haplotype_sections_data[sb_i].len == expected_size, 
                    "Error: unexpected size of haplotype_sections_data[%u]: expecting %u but got %u", sb_i, expected_size, vb->haplotype_sections_data[sb_i].len)
        }
    }
}

#endif

#ifdef V1_VCF_HEADER

// returns true if there's a file or false if its an empty file
bool v1_header_genozip_to_vcf (Md5Hash *digest)
{
    // read vcf header if its not already read. it maybe already read, if we're in flag_split mode, for second component onwards
    // (it would have been read by previous component and stored for us)
    if (!buf_is_allocated (&z_file->v1_next_vcf_header) ||
        z_file->v1_next_vcf_header.len < sizeof(v1_SectionHeaderVCFHeader)) { // first VCF header, data moved here by zfile_read_genozip_header()
        
        int ret;
        ret = v1_zfile_read_section (evb, &z_file->v1_next_vcf_header, "z_file->v1_next_vcf_header", 
                                     sizeof(v1_SectionHeaderVCFHeader), SEC_TXT_HEADER, true);

        if (ret == EOF) {
            buf_free (&z_file->v1_next_vcf_header);
            return false; // empty file (or in case of split mode - no more components) - not an error
        }
    }

    // handle the GENOZIP header of the VCF header section
    v1_SectionHeaderVCFHeader *header = (v1_SectionHeaderVCFHeader *)z_file->v1_next_vcf_header.data;

    z_file->genozip_version = header->genozip_version;
    
    ASSERT0 (header->genozip_version == 1, "Error: v1_header_genozip_to_vcf() can only handle v1 VCF Header sections");

    ASSERT (BGEN32 (header->h.compressed_offset) == crypt_padded_len (sizeof(v1_SectionHeaderVCFHeader)), "Error: invalid VCF header's header size: header->h.compressed_offset=%u, expecting=%u", BGEN32 (header->h.compressed_offset), (unsigned)sizeof(v1_SectionHeaderVCFHeader));

    // in split mode - we open the output VCF file of the component
    if (flag_split) {
        ASSERT0 (!txt_file, "Error: not expecting txt_file to be open already in split mode");
        txt_file = file_open (header->txt_filename, WRITE, TXT_FILE, DATA_TYPE_VCF);
    }

    extern Buffer global_vcf_header_line; // defined in header_vcf.c
    bool first_vcf = !buf_is_allocated (&global_vcf_header_line);

    txt_file->max_lines_per_vb = 4096; // 4096 was the constant value in genozip version 1 - where header->max_lines_per_vb field is absent

    if (first_vcf || !flag_concat) 
        z_file->num_lines = BGEN64 (header->num_lines);

    *digest = flag_split ? header->md5_hash_single : header->md5_hash_concat;
        
    // now get the text of the VCF header itself
    static Buffer vcf_header_buf = EMPTY_BUFFER;
    zfile_uncompress_section (evb, header, &vcf_header_buf, "vcf_header_buf", SEC_TXT_HEADER);

    bool can_concatenate = header_vcf_set_globals (z_file->name, &vcf_header_buf);
    if (!can_concatenate) {
        buf_free (&z_file->v1_next_vcf_header);
        buf_free (&vcf_header_buf);
        return false;
    }

    // write vcf header if not in concat mode, or, in concat mode, we write the vcf header, only for the first genozip file
    if ((first_vcf || !flag_split) && !flag_no_header)
        txtfile_write_to_disk (&vcf_header_buf);
    
    // if we didn't write the header (bc 2nd+ file in concat mode) - just account for it in MD5 (this is normally done by txtfile_write_to_disk())
    else if (flag_md5)
        md5_update (&txt_file->md5_ctx_concat, vcf_header_buf.data, vcf_header_buf.len);

    buf_free (&z_file->v1_next_vcf_header);
    buf_free (&vcf_header_buf);

    return true;
}

// returns the the VCF header section's header from a GENOZIP file - used by main_list
bool v1_vcf_header_get_vcf_header (uint64_t *uncompressed_data_size,
                                   uint32_t *num_samples,
                                   uint64_t *num_items_concat,
                                   Md5Hash  *md5_hash_concat,
                                   char *created, unsigned created_len /* caller allocates space */)
{
    v1_SectionHeaderVCFHeader header;
    int bytes = fread ((char*)&header, 1, crypt_padded_len (sizeof(v1_SectionHeaderVCFHeader)), (FILE *)z_file->file);
    
    if (bytes < sizeof(v1_SectionHeaderVCFHeader)) return false;

    // case: encrypted, and we have a password
    if (BGEN32 (header.h.magic) != GENOZIP_MAGIC && // possibly encrypted
        crypt_have_password()) {
        v1_crypt_do (evb, (uint8_t *)&header, crypt_padded_len (sizeof (v1_SectionHeaderVCFHeader)), 0, -1);
    }

    // case: not encrypted, or encrypted and we successfully decrypted it
    if (BGEN32 (header.h.magic) == GENOZIP_MAGIC) {

        // this is actually just a bad v2+ file
        RETURNW (!(header.h.section_type == SEC_TXT_HEADER && BGEN32 (header.h.compressed_offset) == sizeof(SectionHeaderTxtHeader)), 
                 false,
                 "Error: failed to read file %s - it appears to be truncated or corrupted", z_name);
        
        *uncompressed_data_size = BGEN64 (header.txt_data_size);
        *num_samples            = BGEN32 (header.num_samples);
        *num_items_concat       = BGEN64 (header.num_lines);
        *md5_hash_concat        = header.md5_hash_concat;
        memcpy (created, header.created, MIN (created_len, v1_FILE_METADATA_LEN));
    }

    // case: we cannot read the header (we assume its encrypted) - we just return 0s
    else {
        *uncompressed_data_size = 0;
        *num_samples            = 0;
        *num_items_concat       = 0;
        md5_set_zero (md5_hash_concat);
        memset (created, 0, created_len);
    }

    return true; // all good
}

#endif // V1_VCF_HEADER

#ifdef V1_BASE250

#define BASE250_2_NUMERALS 253 // this number has 2 numerals, starting from numerals[1]
#define BASE250_3_NUMERALS 254 // this number has 3 numerals
#define BASE250_4_NUMERALS 255 // this number has 4 numerals

uint32_t v1_base250_decode (const uint8_t **str)
{
    if ((*str)[0] < 250) {
        (*str)++;
        return (*str)[-1]; // single numeral or control character
    }
    else if ((*str)[0] == BASE250_ONE_UP)     { (*str)++; return WORD_INDEX_ONE_UP; }
    else if ((*str)[0] == BASE250_EMPTY_SF)   { (*str)++; return WORD_INDEX_EMPTY_SF; }
    else if ((*str)[0] == BASE250_MISSING_SF) { (*str)++; return WORD_INDEX_MISSING_SF; }

    uint32_t result = 0;
    uint32_t factor = 1;
    unsigned num_numerals = (*str)[0] - BASE250_2_NUMERALS + 2; // 2, 3 or 4
    for (unsigned i=1; i <= num_numerals; i++) {
        result += (*str)[i] * factor;
        factor *= 250;
    }

    (*str) += num_numerals + 1;

    return result;
}

#endif // V1_BASE250

#ifdef V1_CRYPT

// 256 bit AES is a concatenation of 2 MD5 hashes of the password - each one of length 128 bit
// each hash is a hash of the password concatenated with a constant string
// we add data_len to the hash to give it a near-uniqueness for each section
static void v1_crypt_generate_aes_key (VBlock *vb,                
                                       uint32_t vb_i, int16_t sec_i, // used to generate an aes key unique to each block
                                       uint8_t *aes_key /* out */)
{
    const char *salt   = "frome";     
    const char *pepper = "vaughan";   
    static unsigned pw_len=0, salt_len=0, pepper_len=0;

    if (!pw_len) { // first call
        pw_len     = strlen (password);
        salt_len   = strlen (salt);
        pepper_len = strlen (pepper);
    }

    buf_alloc (vb, &vb->spiced_pw, pw_len + sizeof (uint32_t) + sizeof (int16_t) + salt_len + pepper_len, 1, "spiced_pw", 0);
    buf_add (&vb->spiced_pw, password, pw_len);

    Md5Hash salty_hash, peppered_hash;
    
    // add some salt to the password, mixed with vb_i and sec_i for uniqueness
    buf_add (&vb->spiced_pw, &vb_i, sizeof (uint32_t));
    buf_add (&vb->spiced_pw, &sec_i, sizeof (int16_t));
    buf_add (&vb->spiced_pw, salt, salt_len);
    md5_do (vb->spiced_pw.data, vb->spiced_pw.len, &salty_hash);

    // add some pepper
    buf_add (&vb->spiced_pw, pepper, pepper_len);
    md5_do (vb->spiced_pw.data, vb->spiced_pw.len, &peppered_hash);

    // get hash
    memcpy (aes_key, salty_hash.bytes, sizeof(Md5Hash)); // first half of key
    memcpy (aes_key + sizeof(Md5Hash), peppered_hash.bytes, sizeof(Md5Hash)); // 2nd half of key

    buf_free (&vb->spiced_pw);
}

// we generate a different key for each block by salting the password with vb_i and sec_i
// for sec_i we use (-1-section_i) for the section header and section_i for the section body
// the VCF header section: vb_i=0 (for all components) and sec_i=0 (i.e: 0 for the body, (-1 - 0)=-1 for header)
// the Variant Data section: vb_i={global consecutive number starting at 1}, sec_i=0 (body=0, header=-1)
// Other sections: vb_i same as Variant Data, sec_i consecutive running starting at 1 
void v1_crypt_do (VBlock *vb, uint8_t *data, unsigned data_len, uint32_t vb_i, int16_t sec_i) // used to generate an aes key unique to each block
{
    // generate an AES key just for this one section - combining the pasword with vb_i and sec_i
    uint8_t aes_key[AES_KEYLEN]; 
    v1_crypt_generate_aes_key (vb, vb_i, sec_i, aes_key);

    //fprintf (stderr, "command:%d id:%d vb_i=%d sec_i=%d data_len=%u key=%s\n", command, vb->id, vb_i, sec_i, data_len, aes_display_key (aes_key));

    aes_initialize (vb, aes_key);

    // encrypt in-place
    //fprintf (stderr, "BFRE: data len=%u: %s\n", data_len, aes_display_data (data, data_len));
    aes_xcrypt_buffer (vb, data, data_len);
    //fprintf (stderr, "AFTR: data len=%u: %s\n", data_len, aes_display_data (data, data_len));
}

#endif // V1_CRYPT
