// ------------------------------------------------------------------
//   v1.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// read section header - called from the I/O thread, but for a specific VB
// returns offset of header within data, EOF if end of file

// NOTE: this file contains all the functions that have a different version for v1 vs v2 genozip.
// it is #included at the end of the affected source files

#ifdef V1_ZFILE

// note: the first section (i.e. the first VCF header) is always read by v1_zfile_read_one_section() (the current version)
// if header.genozip_version is 1, then subsequent sections, including subsequent VCF headers, will be read by this function
int v1_zfile_read_one_section (VariantBlock *vb,
                               Buffer *data, const char *buf_name, /* buffer to append */
                               unsigned header_size, SectionType expected_sec_type,
                               bool allow_eof)
{
    bool is_encrypted = crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    unsigned header_offset = expected_sec_type == SEC_VCF_HEADER ? 0 : data->len;
    
    unsigned requested_bytes = header_size;
    if (expected_sec_type == SEC_VCF_HEADER) requested_bytes -= data->len;

    buf_alloc (vb, data, header_offset + header_size, 2, buf_name, 1);

    SectionHeader *header = zfile_read_from_disk (vb, data, requested_bytes, false); // note: header in file can be shorter than header_size if its an earlier version
    ASSERT (header || allow_eof, "Error: Failed to read header, section_type=%s", st_name(expected_sec_type));
    
    // if this is the first VCF header, part of it was already read by zfile_read_genozip_header() and placed in data (which is vb->z_file->v1_next_vcf_header)
    if (expected_sec_type == SEC_VCF_HEADER) header = (SectionHeader *)data->data;

    if (!header) return EOF; 

    // decrypt header
    if (is_encrypted) {
        ASSERT (BGEN32 (header->magic) != GENOZIP_MAGIC, 
                "Error: password provided, but file %s is not encrypted", file_printname (vb->z_file));

        crypt_do (vb, (uint8_t*)header, header_size, vb->variant_block_i, --vb->z_next_header_i); // negative section_i for a header
    }
    bool is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;
    if (!is_magical && !is_encrypted && expected_sec_type == SEC_VCF_HEADER) {

        // file appears to be encrypted but user hasn't provided a password    
        crypt_prompt_for_password();

        unsigned padding;
        is_encrypted = crypt_get_encrypted_len (&header_size, &padding); // update header size if encrypted
            
        if (padding) {
            char *header_extra_bytes = zfile_read_from_disk (vb, data, padding, false);
            ASSERT0 (header_extra_bytes, "Error: Failed to read header padding");
        }

        crypt_do (vb, (uint8_t*)header, header_size, vb->variant_block_i, --vb->z_next_header_i);
        is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;
    }

    if (!is_magical && is_encrypted && expected_sec_type == SEC_VCF_HEADER) {
        ABORT ("Error: password is wrong for file %s", file_printname (vb->z_file)); // mostly likely its because of a wrong password
    }

    // case: encryption failed because this is actually a SEC_VCF_HEADER of a new vcf component of a concatenated file
    bool new_vcf_component = false;
    unsigned new_header_size=0;

    if (!is_magical && is_encrypted && expected_sec_type == SEC_VB_HEADER) {
    
        // reverse failed decryption
        crypt_do (vb, (uint8_t*)header, header_size, vb->variant_block_i, vb->z_next_header_i);

        new_header_size = sizeof (SectionHeaderVCFHeader);
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
            crypt_do (vb, (uint8_t*)header, new_header_size, 0, -1); // vb_i=0 and sec_i=-1 always, for all VCFHeader section headers
            is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;

            vb->z_next_header_i++; // roll back
            new_vcf_component = true; 
        }
    } 

    ASSERT (is_magical, "Error: corrupt data (magic is wrong) when reading file %s", file_printname (vb->z_file));

    unsigned compressed_offset   = BGEN32 (header->compressed_offset);
    ASSERT (compressed_offset, "Error: header.compressed_offset is 0 when reading section_type=%s", st_name(expected_sec_type));

    unsigned data_compressed_len = BGEN32 (header->data_compressed_len);
    ASSERT (data_compressed_len, "Error: header.data_compressed_len is 0 when reading section_type=%s", st_name(expected_sec_type));

    unsigned data_encrypted_len  = BGEN32 (header->data_encrypted_len);

    unsigned data_len = MAX (data_compressed_len, data_encrypted_len);

    // We found a VCF header - possibly a result of concatenating files:
    // if regular mode (no split) - we will just read the data from disk to skip this header (unless its expected)
    // if in split mode - we will end processing this output file here and store the header for the next output vcf file
    bool found_a_vcf_header = header->section_type == SEC_VCF_HEADER;
    
    // check that we received the section type we expect, 
    ASSERT (header->section_type == expected_sec_type || (found_a_vcf_header && expected_sec_type == SEC_VB_HEADER),
            "Error: Unexpected section type: expecting %s, found %s", st_name(expected_sec_type), st_name(header->section_type));

    unsigned expected_header_size = new_vcf_component ? new_header_size : header_size;

    ASSERT (compressed_offset == crypt_padded_len (expected_header_size) || expected_sec_type == SEC_VB_HEADER, // for variant data, we also have the permutation index
            "Error: invalid header - expecting compressed_offset to be %u but found %u", expected_header_size, compressed_offset);

    // allocate more memory for the rest of the header + data (note: after this realloc, header pointer is no longer valid)
    buf_alloc (vb, data, header_offset + compressed_offset + data_len, 2, "v1_zfile_read_one_section", 2);
    header = (SectionHeader *)&data->data[header_offset]; // update after realloc

    // in case we're expecting SEC_VB_HEADER - read the rest of the header: 
    // if we found SEC_VB_HEADER, we need to read haplotype index, and if we found a SEC_VCF_HEADER of a concatenated
    // file componented - we need to read the read of the header as it is bigger than Variant Data header that was read
    if (expected_sec_type == SEC_VB_HEADER && !new_vcf_component) {

        int bytes_left_over = compressed_offset - header_size;
        ASSERT (bytes_left_over >= 0, "Error: expected bytes_left_over=%d to be >=0", bytes_left_over)

        if (bytes_left_over) { // there will be an Index only if this VCF has samples
            
            uint8_t *left_over_data = zfile_read_from_disk (vb, data, bytes_left_over, false);
            ASSERT (left_over_data, "Failed to read left over bytes. bytes_left_over=%u", bytes_left_over);

            if (is_encrypted) { // this is just for the ht index - we've already handle encrypted vcf header and set new_vcf_component=true 
                ASSERT (bytes_left_over == crypt_padded_len (bytes_left_over), "Error: bad length of bytes_left_over=%u", bytes_left_over); // we expected it to be aligned, bc the total encryption block is aligned and header_size is aligned 
                
                // for the haplotype index - it is part of the header so we just continue the encryption stream
                crypt_continue (vb, left_over_data, bytes_left_over);
            }
            if (header->section_type == SEC_VCF_HEADER) 
                new_vcf_component = true;
        }
    }

    // read section data
    if (data_len) {
        ASSERT (zfile_read_from_disk (vb, data, data_len, false), 
                "Error: failed to read section data, section_type=%s", st_name(header->section_type));
    }

    // deal with a VCF header that was encountered while expecting a SEC_VB_HEADER (i.e. 2nd+ component of a concatenated file)
    if (new_vcf_component) {

        if (flag_split) {
            // move the header from vb->z_data to z_file->v1_next_vcf_header
            buf_copy (vb, &vb->z_file->v1_next_vcf_header, data, 1, header_offset, data->len - header_offset, "z_file->v1_next_vcf_header", 0);
            data->len = header_offset; // shrink - remove the vcf header that doesn't belong to this component

            return EOF; // end of VCF component in flag_split mode
        }
        else {
            // since we're not in split mode, so we can skip this vcf header section - just return the next section instead
            return v1_zfile_read_one_section (vb, data, buf_name, header_size, expected_sec_type, allow_eof);
        }
    }

    return header_offset;
}


bool v1_zfile_read_one_vb (VariantBlock *vb)
{ 
    START_TIMER;

    int vardata_header_offset = v1_zfile_read_one_section (vb, &vb->z_data, "z_data",
                                                        sizeof(v1_SectionHeaderVariantData), SEC_VB_HEADER, true);
    if (vardata_header_offset == EOF) {

        // update size - in case they were not known (pipe, gzip etc)
        if (!vb->z_file->disk_size) 
            vb->z_file->disk_size = vb->z_file->disk_so_far;
            
        vb->z_file->eof = true;

        COPY_TIMER (vb->profile.zfile_read_one_vb);
        return false; // end of file
    }

    // note - copying values here z_data.data can get reallocated each call to v1_zfile_read_one_section
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

            unsigned start_i = vb->z_data.len; // vb->z_data.len is updated next, by v1_zfile_read_one_section()
            v1_zfile_read_one_section (vb, &vb->z_data, "z_data", sizeof(SectionHeaderDictionary), SEC_FRMT_SUBFIELD_DICT, false);    

            // update dictionaries in z_file->mtf_ctx with dictionary data from this VB
            mtf_integrate_dictionary_fragment (vb, &vb->z_data.data[start_i]);
        }

        vb->z_data.len = start_dictionary_sections; // shrink z_data back
    }

    // overlay all available dictionaries (not just those that have fragments in this variant block) to the vb
    mtf_overlay_dictionaries_to_vb (vb);

    // read the other sections

    buf_alloc (vb, &vb->z_section_headers, 1000 /* arbitrary initial value */ * sizeof(char*), 0, "z_section_headers", 1);
    
    ((unsigned *)vb->z_section_headers.data)[0] = vardata_header_offset; // variant data header is at index 0

    unsigned section_i=1;

    for (unsigned sb_i=0; sb_i < num_sample_blocks; sb_i++) {

        // make sure we have enough space for the section pointers
        buf_alloc (vb, &vb->z_section_headers, sizeof (unsigned) * (section_i + 3), 2, "z_section_headers", 2);

        if (has_genotype_data) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            v1_zfile_read_one_section (vb, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_GENOTYPE_DATA, false);
        }

        if (phase_type == PHASE_MIXED_PHASED) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            v1_zfile_read_one_section (vb, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_PHASE_DATA, false);
        }

        if (num_haplotypes_per_line) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            v1_zfile_read_one_section (vb, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_HAPLOTYPE_DATA, false);    
        }
    }
    
    COPY_TIMER (vb->profile.zfile_read_one_vb);

    return true; 
}

#endif

#ifdef V1_PIZ

#define V1_SAMPLES_PER_BLOCK 4096

// decode the delta-encoded value of the POS field
static inline void v1_piz_decode_pos (VariantBlock *vb, const char *str,
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

    COPY_TIMER(vb->profile.piz_decode_pos);
}

static void v1_piz_get_line_get_num_subfields (VariantBlock *vb, unsigned line_i, // line in vcf file
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

    COPY_TIMER (vb->profile.piz_get_format_info)
}

static void v1_piz_get_variant_data_line (VariantBlock *vb, 
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
                v1_piz_decode_pos (vb, next+1, pos_str, &delta_pos_start, &delta_pos_len, &add_len);
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
    COPY_TIMER(vb->profile.piz_get_variant_data_line);
    return;
}

static void v1_piz_get_line_subfields (VariantBlock *vb, unsigned line_i, // line in vcf file
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
        unsigned did_i=0 ; for (; did_i < vb->z_file->num_dict_ids; did_i++) 
            if (subfield.num == vb->z_file->mtf_ctx[did_i].dict_id.num) {
                // entry i corresponds to subfield i in FORMAT (excluding GT), and contains the index in mtf_ctx of this subfield
                line_subfields[i] = did_i;
                break;
            }
#ifdef DEBUG
        ASSERTW (did_i < vb->z_file->num_dict_ids, 
                 "Warning: subfield %.*s not found in dictionaries, line=%u. This can happen legitimately if the subfield is declared in FORMAT, but a value is never provided in any sample", DICT_ID_LEN, subfield.id, line_i);
#endif
        if (subfields_start[-1] == '\t' || subfields_start[-1] == '\n') break;
    } 

    COPY_TIMER (vb->profile.piz_get_line_subfields)
}

// convert genotype data from sample block format of indices in base-250 to line format
// of tab-separated genotypes
static void v1_piz_get_genotype_data_line (VariantBlock *vb, unsigned line_i, int *line_subfields)
{
    START_TIMER;

    SnipIterator *sample_iterator = (SnipIterator *)vb->sample_iterator.data; // for convenience

    char *next = vb->line_gt_data.data;
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned first_sample = sb_i*V1_SAMPLES_PER_BLOCK;
        unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);

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
                for (unsigned sf_i=0; sf_i < vb->num_subfields; sf_i++) {

                    if (line_subfields[sf_i] != NIL) {  // this line has this subfield (according to its FORMAT field)

                        // add a colon before, if needed
                        if (snip) *(next++) = ':'; // this works for empty "" snip too

                        unsigned snip_len;
                        mtf_get_next_snip (vb, &vb->mtf_ctx[line_subfields[sf_i]], // note: line_subfields[sf_i] maybe -2 (set in piz_get_line_subfields()), and this is an invalid value. this is ok, bc in this case sample_iterator[sample_i] will be a control character
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
                    "Error: line_gt_data buffer overflow. variant_block_i=%u vcf_line=%u sb_i=%u sample_i=%u",
                    vb->variant_block_i, line_i + vb->first_line, sb_i, sample_i);
        }
    }

    // change last terminator to a \n
    next[-1] = '\n';

    vb->line_gt_data.len = next - vb->line_gt_data.data;

    vb->data_lines[line_i].has_genotype_data = vb->line_gt_data.len > global_num_samples; // not all just \t

    COPY_TIMER(vb->profile.piz_get_genotype_data_line);
}

static void v1_piz_initialize_next_gt_in_sample (VariantBlock *vb, int *num_subfields)
{
    START_TIMER;
    
    buf_alloc (vb, &vb->sample_iterator, sizeof(SnipIterator) * global_num_samples, 1, "sample_iterator", 0);
    SnipIterator *sample_iterator = (SnipIterator *)vb->sample_iterator.data; 
    
    for (unsigned sb_i=0; sb_i < vb->num_sample_blocks; sb_i++) {

        unsigned num_samples_in_sb = vb_num_samples_in_sb (vb, sb_i);
        
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
                    next += base250_len (next, BASE250_ENCODING_8BIT); // v1 files are encoded in 8 bit
        }

        // sanity checks to see we read the correct amount of genotypes
        ASSERT (sample_i == sample_after, "Error: expected to find %u genotypes in sb_i=%u of variant_block_i=%u, but found only %u",
                vb->num_lines * num_samples_in_sb, sb_i, vb->variant_block_i, vb->num_lines * (sample_i - sb_i * V1_SAMPLES_PER_BLOCK));

        ASSERT (next == after, "Error: expected to find %u genotypes in sb_i=%u of variant_block_i=%u, but found more. ",
                vb->num_lines * num_samples_in_sb, sb_i, vb->variant_block_i);
    }

    COPY_TIMER (vb->profile.piz_initialize_sample_iterators)
}

// combine all the sections of a variant block to regenerate the variant_data, haplotype_data,
// genotype_data and phase_data for each row of the variant block
void v1_piz_reconstruct_line_components (VariantBlock *vb)
{
    START_TIMER;

    if (!vb->data_lines) 
        vb->data_lines = calloc (global_max_lines_per_vb, sizeof (DataLine));

    // initialize phase data if needed
    if (vb->phase_type == PHASE_MIXED_PHASED) 
        buf_alloc (vb, &vb->line_phase_data, global_num_samples, 1, "line_phase_data", 0);

    // initialize haplotype stuff
    const char **ht_columns_data=NULL;
    if (vb->has_haplotype_data) {

        //  memory - realloc for exact size, add 7 because depermuting_loop works on a word (32/64 bit) boundary
        buf_alloc (vb, &vb->line_ht_data, vb->num_haplotypes_per_line + 7, 1, "line_ht_data", 0);

        ht_columns_data = piz_get_ht_columns_data (vb);
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

    // initialize again - for piz_get_variant_data_line
    variant_data_next_line = vb->v1_variant_data_section_data.data;
    variant_data_length_remaining = vb->v1_variant_data_section_data.len;

    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {

        // de-permute variant data into vb->line_variant_data
        v1_piz_get_variant_data_line (vb, vb->first_line + line_i, &variant_data_length_remaining, &variant_data_next_line);

        // reset len for next line - no need to realloc as we have realloced what is needed already
        vb->line_ht_data.len = vb->line_gt_data.len = vb->line_phase_data.len = 0;

        // transform sample blocks (each block: n_lines x s_samples) into line components (each line: 1 line x ALL_samples)
        if (vb->has_genotype_data)  {
            int line_subfields[MAX_SUBFIELDS]; // entry i corresponds to subfield i in FORMAT (excluding GT), and contains the index in mtf_ctx of this subfield
            v1_piz_get_line_subfields (vb, vb->first_line + line_i,
                                       subfields_start[line_i], subfields_len[line_i], line_subfields);

            v1_piz_get_genotype_data_line (vb, line_i, line_subfields);
        }
        if (vb->phase_type == PHASE_MIXED_PHASED) 
            piz_get_phase_data_line (vb, line_i);

        if (vb->has_haplotype_data) 
            piz_get_haplotype_data_line (vb, line_i, ht_columns_data);

        piz_merge_line (vb, line_i);
    }

    COPY_TIMER(vb->profile.piz_reconstruct_line_components);
}

void v1_piz_uncompress_all_sections (VariantBlock *vb)
{
    unsigned *section_index = (unsigned *)vb->z_section_headers.data;

    // get the variant data - newline-seperated lines, each containing the first 8 (if no FORMAT field) or 9 fields (if FORMAT exists)
    zfile_uncompress_section (vb, vb->z_data.data + section_index[0], 
                              &vb->v1_variant_data_section_data, "v1_variant_data_section_data", SEC_VB_HEADER);
    
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
    // num_dictionary_sections is read in zfile_read_one_vb()
    vb->max_gt_line_len         = BGEN32 (vardata_header->max_gt_line_len);
    vb->vb_data_size            = BGEN32 (vardata_header->vb_data_size);
    
    ASSERT (global_num_samples == BGEN32 (vardata_header->num_samples), "Error: Expecting variant block to have %u samples, but it has %u", global_num_samples, BGEN32 (vardata_header->num_samples));

    // we allocate memory for the Buffer arrays only once the first time this VariantBlock
    // is used. Subsequent blocks reusing the memory will have the same number of samples (by VCF spec)
    // BUG: this won't work if we're doing mutiple unrelated VCF on the command line
    if (vb->has_genotype_data && !vb->genotype_sections_data) 
        vb->genotype_sections_data  = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));

    if (vb->phase_type == PHASE_MIXED_PHASED && !vb->phase_sections_data) 
        vb->phase_sections_data     = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));
    
    if (vb->num_haplotypes_per_line && !vb->haplotype_sections_data) 
        vb->haplotype_sections_data = (Buffer *)calloc (vb->num_sample_blocks, sizeof (Buffer));

    // unsqueeze permutation index - if this VCF has samples
    if (global_num_samples) {
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

        unsigned num_samples_in_sb = (sb_i == vb->num_sample_blocks-1 ? global_num_samples % vb->num_samples_per_block : vb->num_samples_per_block);

        // if genotype data exists, it appears first
        if (vb->has_genotype_data) 
            zfile_uncompress_section (vb, vb->z_data.data + section_index[section_i++], 
                                      &vb->genotype_sections_data[sb_i], "genotype_sections_data", SEC_GENOTYPE_DATA);
        
        // next, comes phase data
        if (vb->phase_type == PHASE_MIXED_PHASED) {
            
            zfile_uncompress_section (vb, vb->z_data.data + section_index[section_i++], 
                                      &vb->phase_sections_data[sb_i], "phase_sections_data", SEC_PHASE_DATA);
            
            unsigned expected_size = vb->num_lines * num_samples_in_sb;
            ASSERT (vb->phase_sections_data[sb_i].len==expected_size, 
                    "Error: unexpected size of phase_sections_data[%u]: expecting %u but got %u", sb_i, expected_size, vb->phase_sections_data[sb_i].len)
        }

        // finally, comes haplotype data
        if (vb->has_haplotype_data) {
            
            zfile_uncompress_section (vb, vb->z_data.data + section_index[section_i++], 
                                      &vb->haplotype_sections_data[sb_i], "haplotype_sections_data", SEC_HAPLOTYPE_DATA);
            
            unsigned expected_size = vb->num_lines * num_samples_in_sb * vb->ploidy;
            ASSERT (vb->haplotype_sections_data[sb_i].len == expected_size, 
                    "Error: unexpected size of haplotype_sections_data[%u]: expecting %u but got %u", sb_i, expected_size, vb->haplotype_sections_data[sb_i].len)
        }
    }

    // all dictionaries are 8bit in v1
    for (unsigned did_i=0; did_i < MAX_DICTS; did_i++)
        if (vb->mtf_ctx[did_i].b250_section_type == SEC_GENOTYPE_DATA)
            vb->mtf_ctx[did_i].encoding = BASE250_ENCODING_8BIT;

}

#endif

#ifdef V1_VCF_HEADER

// returns true if there's a file or false if its an empty file
bool v1_vcf_header_genozip_to_vcf (VariantBlock *vb, Md5Hash *digest)
{
    // read vcf header if its not already read. it maybe already read, if we're in flag_split mode, for second component onwards
    // (it would have been read by previous component and stored for us)
    if (!buf_is_allocated (&vb->z_file->v1_next_vcf_header) ||
        vb->z_file->v1_next_vcf_header.len < sizeof(v1_SectionHeaderVCFHeader)) { // first VCF header, data moved here by zfile_read_genozip_header()
        
        int ret;
        ret = v1_zfile_read_one_section (vb, &vb->z_file->v1_next_vcf_header, "z_file->v1_next_vcf_header", 
                                         sizeof(v1_SectionHeaderVCFHeader), SEC_VCF_HEADER, true);

        if (ret == EOF) {
            buf_free (&vb->z_file->v1_next_vcf_header);
            return false; // empty file (or in case of split mode - no more components) - not an error
        }
    }

    // handle the GENOZIP header of the VCF header section
    v1_SectionHeaderVCFHeader *header = (v1_SectionHeaderVCFHeader *)vb->z_file->v1_next_vcf_header.data;

    vb->z_file->genozip_version = header->genozip_version;
    
    ASSERT0 (header->genozip_version == 1, "Error: v1_vcf_header_genozip_to_vcf() can only handle v1 VCF Header sections");

    ASSERT (BGEN32 (header->h.compressed_offset) == crypt_padded_len (sizeof(v1_SectionHeaderVCFHeader)), "Error: invalid VCF header's header size: header->h.compressed_offset=%u, expecting=%u", BGEN32 (header->h.compressed_offset), (unsigned)sizeof(v1_SectionHeaderVCFHeader));

    // in split mode - we open the output VCF file of the component
    if (flag_split) {
        ASSERT0 (!vb->vcf_file, "Error: not expecting vb->vcf_file to be open already in split mode");
        vb->vcf_file = file_open (header->vcf_filename, WRITE, VCF);
    }

    extern Buffer global_vcf_header_line; // defined in vcf_header.c
    bool first_vcf = !buf_is_allocated (&global_vcf_header_line);

    global_max_lines_per_vb = 4096; // 4096 was the constant value in genozip version 1 - where header->max_lines_per_vb field is absent

    if (first_vcf || !flag_concat_mode) {
        vb->z_file->num_lines_concat     = vb->vcf_file->num_lines_concat     = BGEN64 (header->num_lines);
        vb->z_file->vcf_data_size_concat = vb->vcf_file->vcf_data_size_concat = BGEN64 (header->vcf_data_size);
    }

    *digest = flag_split ? header->md5_hash_single : header->md5_hash_concat;
    vb->vcf_file->has_md5 = !md5_is_zero (*digest); // has_md5 iff not all 0. note: a chance of 1 in about 10^38 that we will get all-0 by chance in which case will won't perform the md5 comparison
        
    // now get the text of the VCF header itself
    static Buffer vcf_header_buf = EMPTY_BUFFER;
    zfile_uncompress_section (vb, header, &vcf_header_buf, "vcf_header_buf", SEC_VCF_HEADER);

    bool can_concatenate = vcf_header_set_globals (vb, vb->z_file->name, &vcf_header_buf);
    if (!can_concatenate) {
        buf_free (&vb->z_file->v1_next_vcf_header);
        buf_free (&vcf_header_buf);
        return false;
    }

    // write vcf header if not in concat mode, or, in concat mode, we write the vcf header, only for the first genozip file
    if (first_vcf || !flag_concat_mode)
        vcffile_write_to_disk (vb->vcf_file, &vcf_header_buf);
    
    // if we didn't write the header (bc 2nd+ file in concat mode) - just account for it in MD5 if needed (this is normally done by vcffile_write_to_disk())
    else if (vb->vcf_file->has_md5)
        md5_update (&vb->vcf_file->md5_ctx_concat, vcf_header_buf.data, vcf_header_buf.len, true);

    buf_free (&vb->z_file->v1_next_vcf_header);
    buf_free (&vcf_header_buf);

    return true;
}

#endif