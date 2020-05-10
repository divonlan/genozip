
// ------------------------------------------------------------------
//   piz.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "txtfile.h"
#include "header.h"
#include "vblock.h"
#include "base250.h"
#include "dispatcher.h"
#include "move_to_front.h"
#include "file.h"
#include "endianness.h"
#include "piz.h"
#include "sections.h"
#include "random_access.h"
#include "regions.h"
#include "samples.h"
#include "dict_id.h"
#include "strings.h"
#include "seg.h"

static Buffer piz_iname_mapper_buf = EMPTY_BUFFER;

// Compute threads: decode the delta-encoded value of the POS field, and returns the new last_pos
int32_t piz_decode_pos (VBlock *vb, uint32_t txt_line_i,
                        int32_t last_pos, const char *delta_snip, unsigned delta_snip_len,
                        int32_t *last_delta, // optional in/out
                        char *pos_str, unsigned *pos_len) // out
{
    int32_t delta=0;

    // case: this is not a delta and snip is stored in dictionary (a result of a "nonsense number")
    if (delta_snip_len && delta_snip[0] == POS_NONSENSE) { 

        ASSERT0 (!last_delta, "Error: piz_decode_pos requires last_delta==NULL in delta fields that might be nonsense");
        
        memcpy (pos_str, delta_snip+1, delta_snip_len-1);
        *pos_len = delta_snip_len - 1;
        return last_pos; // unchanged
    }

    // case: this is not a delta and pos is stored in random_pos (a result of the delta being out of range)
    if (delta_snip[0] == POS_LOOKUP) { 
        ASSERT (delta_snip_len==1, "Expected POS_LOOKUP to be of length 1, but its length is %u", delta_snip_len);
        ASSERT (vb->next_random_pos < vb->random_pos_data.len, "Error reading sam_line=%u: unexpected end of RANDOM_POS data", txt_line_i);
        last_pos = BGEN32 (*ENT (uint32_t, vb->random_pos_data, vb->next_random_pos++));
        str_int (last_pos, pos_str, pos_len);
        return last_pos;
    }

    if (!delta_snip_len && last_delta)
        delta = -(*last_delta); // negated delta of last line - happens every other line in unsorted BAMs

    else {
        // we parse the string ourselves - this is hugely faster than sscanf.
        unsigned negative = *delta_snip == '-'; // 1 or 0

        const char *s; for (s=(delta_snip + negative); *s != '\t' ; s++)
            delta = delta*10 + (*s - '0');

        if (negative) delta = -delta;
    }

    last_pos += delta;    
    ASSERT (last_pos >= 0, "Error: last_pos=%d is negative", last_pos);

    if (last_delta) *last_delta = delta;

    // create number string without calling slow sprintf - start with creating the reverse string
    char reverse_pos_str[50];
    uint32_t n = last_pos;

    unsigned len=0; 
    if (n) {
        while (n) {
            reverse_pos_str[len++] = '0' + (n % 10);
            n /= 10;
        }

        // reverse it
        for (unsigned i=0; i < len; i++) pos_str[i] = reverse_pos_str[len-i-1];
        pos_str[len] = '\0';

        *pos_len = len;
    }
    else {  // n=0. 
        pos_str[0] = '0';
        pos_str[1] = '\0';
        *pos_len = 1;
    }

    return last_pos;
}

void piz_reconstruct_id (VBlock *vb, Buffer *id_buf, uint32_t *next_id, const char *id_snip, unsigned id_snip_len, 
                         bool *extra_bit, bool add_tab)
{
    bool my_extra_bit = (id_snip_len && id_snip[id_snip_len-1] == 2);
    if (my_extra_bit) id_snip_len--;
    if (extra_bit) *extra_bit = my_extra_bit;

    bool has_numeric = (id_snip_len && id_snip[id_snip_len-1] == 1);
    if (has_numeric) id_snip_len--;

    RECONSTRUCT (id_snip, id_snip_len);

    if (has_numeric) {
        uint32_t num = BGEN32 (*ENT (uint32_t, *id_buf, (*next_id)++));
        char num_str[20];
        unsigned num_str_len;
        str_int (num, num_str, &num_str_len);
        RECONSTRUCT (num_str, num_str_len);
    }
    
    if (add_tab) RECONSTRUCT1 ("\t");
}

void piz_uncompress_compound_field (VBlock *vb, SectionType field_b250_sec, SectionType sf_b250_sec, 
                                    SubfieldMapper *mapper, unsigned *section_i)
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);

    for (uint8_t sf_i=0; sf_i < mapper->num_subfields ; sf_i++) {
        
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[(*section_i)++]);
        
        if (zfile_is_skip_section (vb, field_b250_sec, DICT_ID_NONE)) continue; // don't create ctx if section is skipped

        MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, NULL, 
                                                  header->dict_id, sf_b250_sec-1);

        mapper->did_i[sf_i] = ctx->did_i;

        zfile_uncompress_section ((VBlockP)vb, header, &ctx->b250, "mtf_ctx.b250", sf_b250_sec);    
    }    
}

void piz_reconstruct_compound_field (VBlock *vb, SubfieldMapper *mapper, const char *separator, unsigned separator_len,
                                     const char *template, unsigned template_len, uint32_t txt_line_i)
{
    if (template_len==1 && template[0]=='*') template_len=0; // no separators, only one component

    for (unsigned i=0; i <= template_len; i++) { // for template_len, we have template_len+1 elements
        
        DECLARE_SNIP;
        mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[mapper->did_i[i]], NULL, &snip, &snip_len, txt_line_i);

        RECONSTRUCT (snip, snip_len);
        
        if (i < template_len)
            // add middle-field separator. note: in seg_compound_field we re-wrote \t at 1, now we re-write back
            RECONSTRUCT (template[i]==DICT_TAB_REWRITE_CHAR ? "\t" : &template[i], 1) // buf_add is a macro - no semicolon
        else if (separator_len)
            RECONSTRUCT (separator, separator_len); // add end-of-field separator if needed
    }
}

// Called from the I/O thread piz_*_read_one_vb - and immutable thereafter
// Constructs the iname mapper - mapping an iname word like "I1=I2=I3" to its subfields
// This uses dictionary data only, not b250 data, and hence can be done once after reading the dictionaries
void piz_map_iname_subfields (void)
{
    // terminology: we call a list of INFO subfield names, an "iname". An iname looks something like
    // this: "I1=I2=I3=". Each iname consists of info subfields. These fields are not unique to this
    // iname and can appear in other inames. The INFO field contains the iname, and values of the subfields.
    // iname_mapper maps these subfields. This function creates an iname_mapper for every unique iname.

    // in genozip genozip v2 and v3 we had a bug where when there was an INFO subfield with a value following one or more
    // subfields without a value, then the dictionary name of included the names of the valueless subfields
    // example: I2=a;N1;N2;I3=x - the name of the 2nd dictionary was "N1;N2;I3"
    bool v2v3_bug = !is_v4_or_above; // recover from a bug we had in v2/v3 VCF (we had only VCF in v2/3). See details in v2v3.c.

    const MtfContext *info_ctx = &z_file->mtf_ctx[DTFZ(info)];
    
    buf_free (&piz_iname_mapper_buf); // in case it was allocated by a previous file
    piz_iname_mapper_buf.len = info_ctx->word_list.len;
    buf_alloc (evb, &piz_iname_mapper_buf, sizeof (PizSubfieldMapper) * piz_iname_mapper_buf.len,
               1, "piz_iname_mapper_buf", 0);
    buf_zero (&piz_iname_mapper_buf);

    ARRAY (PizSubfieldMapper, all_iname_mappers, piz_iname_mapper_buf);
    ARRAY (const MtfWord, all_inames, info_ctx->word_list);

    for (unsigned iname_i=0; iname_i < piz_iname_mapper_buf.len; iname_i++) {

        // e.g. "I1=I2=I3=" - pointer into the INFO dictionary - immutable - global shared with all threads
        char *iname = ENT (char, info_ctx->dict, all_inames[iname_i].char_index); 
        unsigned iname_len = all_inames[iname_i].snip_len; 

        PizSubfieldMapper *iname_mapper = &all_iname_mappers[iname_i]; // iname_mapper of this specific set of names "I1=I2=I3="

        // get INFO subfield snips - which are the values of the INFO subfield, where the names are
        // in the INFO snip in the format "info1=info2=info3="). 
        iname_mapper->num_subfields = 0;

        // traverse the subfields of one iname. E.g. if the iname is "I1=I2=I3=" then we traverse I1, I2, I3
        for (unsigned i=0; i < iname_len; i++) {
            
            DictIdType *dict_id = &iname_mapper->dict_id[iname_mapper->num_subfields];

            // traverse the iname, and get the dict_id for each subfield name (using only the first 8 characers)
            dict_id->num = 0;
            unsigned j=0; 
            while (iname[i] != '=' && (v2v3_bug || iname[i] != ';') && iname[i] != '\t') { // value-less INFO names can be terminated by the end-of-word \t in the dictionary
                if (j < DICT_ID_LEN) 
                    dict_id->id[j] = iname[i]; // scan the whole name, but copy only the first 8 bytes to dict_id
                i++, j++;
            }
            
            if (is_v5_or_above) 
                *dict_id = dict_id_type_1 (dict_id_make (&iname[i-j], j));
            else {
                // up to version 4 the dict_id was the (up to) 8 first characters
                memcpy (dict_id->id, &iname[i-j], MIN (j, DICT_ID_LEN)); 
                *dict_id = dict_id_type_1 (*dict_id);
            }

            // case - INFO has a special added name "#" indicating that this VCF line has a Windows-style \r\n ending
            if (dict_id->num == dict_id_WindowsEOL) {
                iname_mapper->did_i[iname_mapper->num_subfields] = DID_I_HAS_13;
                
                // we now modify the global dictionary to match the true iname strings, so it can be used for reconstruction
                if (i>=2 && iname[i-2] == ';')
                    iname[iname_len-2] = '\t'; // chop off the ";#" at the end of the INFO word in the dictionary (happens if previous subfield is value-less)
                else     
                    iname[iname_len-1] = '\t'; // chop off the "#" at the end of the INFO word in the dictionary (happens if previous subfield has a value)
            }
            else 
                iname_mapper->did_i[iname_mapper->num_subfields] = mtf_get_existing_did_i_by_dict_id (*dict_id); // it will be NIL if this is an INFO name without values            
            
            iname_mapper->num_subfields++;
        }
    }
}

// used for VCF INFO and GFF3 ATTRIBUTES 
void piz_reconstruct_info (VBlock *vb, uint32_t iname_word_index, 
                           const char *iname_snip, unsigned iname_snip_len, 
                           PizReconstructSpecialInfoSubfields reconstruct_special_info_subfields,
                           uint32_t txt_line_i, bool *has_13)
{
    *has_13 = false; // false unless proven true

    // in genozip genozip v2 and v3 we had a bug where when there was an INFO subfield with a value following one or more
    // subfields without a value, then the dictionary name of included the names of the valueless subfields
    // example: I2=a;N1;N2;I3=x - the name of the 2nd dictionary was "N1;N2;I3"
    bool v2v3_bug = !is_v4_or_above; 

    ASSERT (iname_word_index < piz_iname_mapper_buf.len, "Error: expected iname_word_index=%d < piz_iname_mapper_buf.len=%u", 
            iname_word_index, (uint32_t)piz_iname_mapper_buf.len);

    const PizSubfieldMapper *iname_mapper = ENT (PizSubfieldMapper, piz_iname_mapper_buf, iname_word_index);

    for (unsigned sf_i = 0; sf_i < iname_mapper->num_subfields; sf_i++) {

        uint8_t did_i = iname_mapper->did_i[sf_i]; // DID_I_NONE for fields that have no ctx (either because they have no values or because the values are stored elsewhere)
        
        if (did_i == DID_I_HAS_13) {
            *has_13 = true; // line needs to end with a \r\n
            continue;
        }
        
        // get the name eg "AF="
        const char *start = iname_snip;
        for (; *iname_snip != '=' && (*iname_snip != ';' || v2v3_bug) && *iname_snip != '\t'; iname_snip++);
        bool has_value = (*iname_snip == '=');
        
        if (has_value) iname_snip++; // move past the '=', if one exists
        RECONSTRUCT (start, (unsigned)(iname_snip-start)); // name inc. '=' e.g. "Info1="

        // handle the value part of "name=value", if there is one
        if (has_value) {
        
            bool regular = reconstruct_special_info_subfields (vb, did_i, iname_mapper->dict_id[sf_i], txt_line_i);
            
            if (regular) { // no special treatment - we proceed with the regular treatment
                ASSERT (did_i != DID_I_NONE, "Error in piz_reconstruct_info: did_i=DID_I_NONE for dict_id=%s", err_dict_id (iname_mapper->dict_id[sf_i]));
                DECLARE_SNIP;
                LOAD_SNIP (did_i);
                RECONSTRUCT (snip, snip_len); // value e.g "value1"
            }
        }

        RECONSTRUCT1 (";"); // seperator between each two name=value pairs e.g "name1=value;name2=value2"
        if (*iname_snip == ';') iname_snip++;
    }

    // remove the semicolon for the last field. easier to remove now than calculate the cases involving has_13.
    vb->txt_data.len -= 1;
}

void piz_reconstruct_seq_qual (VBlock *vb, uint32_t seq_len, 
                               const Buffer *data, uint32_t *next,
                               SectionType sec, uint32_t txt_line_i, bool grepped_out)
{
    // seq and qual are expected to be either of length seq_len, or " " (unavailable)
    uint32_t len = (*next >= data->len || data->data[*next] == ' ') ? 1 : seq_len;
    ASSERT (*next + len <= data->len, "Error reading txt_line=%u: unexpected end of %s data", txt_line_i, st_name (sec));

    // rewrite unavailable back to "*"
    if (data->data[*next] == ' ') data->data[*next] = '*';

    if (!grepped_out && !zfile_is_skip_section (vb, sec, DICT_ID_NONE)) 
        RECONSTRUCT (&data->data[*next], len);
    
    *next += len;
}

// Called by PIZ I/O thread: read all the sections at the end of the file, before starting to process VBs
static DataType piz_read_global_area (Md5Hash *original_file_digest) // out
{
    DataType data_type = zfile_read_genozip_header (original_file_digest);
    
    dict_id_initialize(); // must run after zfile_read_genozip_header that sets z_file->data_type; needed by V1 too

    if (data_type == DT_VCF_V1 || data_type == DT_NONE) return data_type;
    
    // for FASTA and FASTQ we convert a "header_only" flag to "header_one" for consistency with other formats (header-only implies we don't show any data)
    if (flag_header_only && (data_type == DT_FASTA || data_type == DT_FASTQ)) {
        flag_header_only = false;
        flag_header_one = true;
    }

    // if the user wants to see only the header, we can skip the dictionaries, regions and random access
    if (!flag_header_only) {
        
        // read random access, but only if we are going to need it
        if (flag_regions || flag_show_index) {
            zfile_read_all_dictionaries (0, DICTREAD_CHROM_ONLY); // read all CHROM/RNAME dictionaries - needed for regions_make_chregs()

            // update chrom node indeces using the CHROM dictionary, for the user-specified regions (in case -r/-R were specified)
            regions_make_chregs (dt_fields[data_type].chrom);

            // if the regions are negative, transform them to the positive complement instead
            regions_transform_negative_to_positive_complement();

            SectionListEntry *ra_sl = sections_get_offset_first_section_of_type (SEC_RANDOM_ACCESS);
            zfile_read_section (evb, 0, NO_SB_I, &evb->z_data, "z_data", sizeof (SectionHeader), SEC_RANDOM_ACCESS, ra_sl);

            zfile_uncompress_section (evb, evb->z_data.data, &z_file->ra_buf, "z_file->ra_buf", SEC_RANDOM_ACCESS);

            z_file->ra_buf.len /= random_access_sizeof_entry();
            BGEN_random_access();

            if (flag_show_index) {
                random_access_show_index(false);
                if (exe_type == EXE_GENOCAT) exit(0); // in genocat --show-index, we only show the index, not the data
            }

            buf_free (&evb->z_data);
        }

        // get the last vb_i that included in the regions - returns -1 if no vb has the requested regions
        int32_t last_vb_i = flag_regions ? random_access_get_last_included_vb_i() : 0;

        // read dictionaries (this also seeks to the start of the dictionaries)
        if (last_vb_i >= 0)
            zfile_read_all_dictionaries (last_vb_i, flag_regions ? DICTREAD_EXCEPT_CHROM : DICTREAD_ALL);
    }
    
    file_seek (z_file, 0, SEEK_SET, false);

    return data_type;
}

static void enforce_v1_limitations (bool is_first_component)
{
    #define ENFORCE(flag,lflag) ASSERT (!(flag), "Error: %s option is not supported because %s was compressed with genozip version 1", (lflag), z_name);
    
    ENFORCE(flag_test, "--test");
    ENFORCE(flag_split, "--split");
    ENFORCE(flag_regions, "--regions");
    ENFORCE(flag_samples, "--samples");
    ENFORCE(flag_show_b250, "--show-b250");
    ENFORCE(flag_show_dict, "--show-dict");
    ENFORCE(dict_id_show_one_b250.num, "--show-one-b250");
    ENFORCE(dict_id_show_one_dict.num, "--show-one-dict");
    ENFORCE(dict_id_dump_one_b250.num, "--dump-one-b250");
    ENFORCE(flag_show_gheader, "--show-gheader");
    ENFORCE(flag_show_index, "--show-index");
    ENFORCE(flag_show_headers, "--show-headers");
    ENFORCE(flag_drop_genotypes, "--drop-genotypes");
    ENFORCE(flag_gt_only, "--flag_gt_only");
}

// returns true is successfully outputted a txt file
bool piz_dispatcher (const char *z_basename, unsigned max_threads, 
                     bool is_first_component, bool is_last_file)
{
    // static dispatcher - with flag_split, we use the same dispatcher when unzipping components
    static Dispatcher dispatcher = NULL;
    bool piz_successful = false;
    SectionListEntry *sl_ent = NULL;
    
    if (flag_split && !sections_has_more_components()) return false; // no more components

    if (!dispatcher) 
        dispatcher = dispatcher_init (max_threads, 0, flag_test, is_last_file, z_basename);
    
    // read genozip header
    Md5Hash original_file_digest;

    // read genozip header, dictionaries etc and set the data type when reading the first component of in case of --split, 
    static DataType data_type = DT_NONE; 
    if (is_first_component) {
        data_type = piz_read_global_area (&original_file_digest);

        if (data_type != DT_VCF_V1)  // genozip v2+ - move cursor past first txt header
            ASSERT (sections_get_next_header_type(&sl_ent, NULL, NULL) == SEC_TXT_HEADER, "Error: unable to find TXT Header data in %s", z_name);

        ASSERT (!flag_test || !md5_is_zero (original_file_digest), 
                "Error testing %s: --test cannot be used with this file, as it was not compressed with --md5 or --test", z_name);
    }

    // case: we couldn't open the file because we didn't know what type it is - open it now
    if (!flag_split && !txt_file->file) file_open_txt (txt_file);

    if (data_type == DT_NONE) goto finish;

    if (!is_v2_or_above) enforce_v1_limitations (is_first_component); // genozip_version will be 0 for v1, bc we haven't read the vcf header yet

    // read and write txt header. in split mode this also opens txt_file
    piz_successful = (data_type != DT_VCF_V1) ? header_genozip_to_txt (&original_file_digest)
                                              : v1_header_genozip_to_vcf (&original_file_digest);
    
    ASSERT (piz_successful || !is_first_component, "Error: failed to read %s header in %s", 
            dt_name (z_file->data_type), z_name);

    if (!piz_successful || flag_header_only) goto finish;

    if (flag_split) 
        dispatcher_resume (dispatcher); // accept more input 

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. In input is not exhausted, and a compute thread is available - read a variant block and compute it
    // 2. Wait for the first thread (by sequential order) to complete and write data

    bool header_only_file = true; // initialize
    do {
        // PRIORITY 1: In input is not exhausted, and a compute thread is available - read a variant block and compute it
        if (!dispatcher_is_input_exhausted (dispatcher) && dispatcher_has_free_thread (dispatcher)) {

            bool still_more_data = false, grepped_out = false;
            if (is_v2_or_above) {
                
                bool skipped_vb;
                static Buffer region_ra_intersection_matrix = EMPTY_BUFFER; // we will move the data to the VB when we get it
                SectionType header_type = sections_get_next_header_type (&sl_ent, &skipped_vb, &region_ra_intersection_matrix);
                switch (header_type) {
                    case SEC_VB_HEADER: {

                        // if we skipped VBs or we skipped the sample sections in the last vb (flag_drop_genotypes), we need to seek forward 
                        if (skipped_vb || flag_drop_genotypes) file_seek (z_file, sl_ent->offset, SEEK_SET, false); 

                        VBlock *next_vb = dispatcher_generate_next_vb (dispatcher, sl_ent->vblock_i);
                        
                        if (region_ra_intersection_matrix.data) {
                            buf_copy (next_vb, &next_vb->region_ra_intersection_matrix, &region_ra_intersection_matrix, 0,0,0, "region_ra_intersection_matrix", next_vb->vblock_i);
                            buf_free (&region_ra_intersection_matrix); // note: copy & free rather than move - so memory blocks are preserved for VB re-use
                        }
                        
                        // read one VB's genozip data
                        DTPZ(read_one_vb)(next_vb);

                        still_more_data = true; // not eof yet

                        // if grep found nothing next_vb->ready_to_dispatch is set to false
                        if (!next_vb->ready_to_dispatch) {
                            grepped_out = true;
                            dispatcher_abandon_next_vb (dispatcher); 
                        }

                        break;
                    }

                    case SEC_TXT_HEADER: // 2nd+ txt header of a concatenated file
                        if (!flag_split) {
                            header_genozip_to_txt (NULL); // skip 2nd+ txt header if concatenating
                            continue;
                        }
                        break; // eof if splitting
                    
                    case SEC_EOF: 
                        break; 
                    
                    default: ABORT ("Error in piz_dispatcher: unexpected section_type=%s", st_name (header_type));
                }
            }
            else still_more_data = v1_piz_vcf_read_one_vb ((VBlockVCF *)dispatcher_generate_next_vb (dispatcher, 0));  // genozip v1
            
            if (still_more_data) {
                if (!grepped_out) 
                    dispatcher_compute (dispatcher, DTPZ(uncompress));
                    
                header_only_file = false;                
            }
            else { // eof
                dispatcher_input_exhausted (dispatcher);

                if (header_only_file)
                    dispatcher_finalize_one_vb (dispatcher);
            }
        }

        // PRIORITY 2: Wait for the first thread (by sequential order) to complete and write data
        else { // if (dispatcher_has_processed_vb (dispatcher, NULL)) {
            VBlock *processed_vb = dispatcher_get_processed_vb (dispatcher, NULL); 

            txtfile_write_one_vblock (processed_vb);
            z_file->num_vbs++;
            
            z_file->txt_data_so_far_single += processed_vb->vb_data_size; 

            dispatcher_finalize_one_vb (dispatcher);
        }

    } while (!dispatcher_is_done (dispatcher));

    // verify file integrity, if the genounzip compress was run with --md5 or --test
    if (flag_md5) {
        Md5Hash decompressed_file_digest;
        md5_finalize (&txt_file->md5_ctx_concat, &decompressed_file_digest); // z_file might be a concatenation - this is the MD5 of the entire concatenation

        if (md5_is_zero (original_file_digest) && !flag_quiet) 
            fprintf (stderr, "MD5 = %s Note: unable to compare this to the original file as file was not originally compressed with --md5\n", md5_display (&decompressed_file_digest, false));
        
        else if (md5_is_equal (decompressed_file_digest, original_file_digest)) {

            if (flag_test && !flag_quiet) fprintf (stderr, "Success          \b\b\b\b\b\b\b\b\b\b\n");

            if (flag_md5 && !flag_quiet) 
                fprintf (stderr, "MD5 = %s verified as identical to the original %s\n", 
                         md5_display (&decompressed_file_digest, false), dt_name (txt_file->data_type));
        }
        else if (flag_test) 
            fprintf (stderr, "FAILED!!!          \b\b\b\b\b\b\b\b\b\b\nError: MD5 of original file=%s is different than decompressed file=%s\nPlease contact bugs@genozip.com to help fix this bug in genozip",
                    md5_display (&original_file_digest, false), md5_display (&decompressed_file_digest, false));
            
        else ASSERT (md5_is_zero (original_file_digest), // its ok if we decompressed only a partial file, or its a v1 files might be without md5
                     "File integrity error: MD5 of decompressed file %s is %s, but the original %s file's was %s", 
                     txt_file->name, md5_display (&decompressed_file_digest, false), dt_name (txt_file->data_type), 
                     md5_display (&original_file_digest, false));
    }

    if (flag_split) file_close (&txt_file, true); // close this component file

    if (!flag_test && !flag_quiet) 
        fprintf (stderr, "Done (%s)           \n", dispatcher_ellapsed_time (dispatcher, false));

finish:

    // in split mode - we continue with the same dispatcher in the next component. otherwise, we finish with it here
    if (!flag_split || !piz_successful) 
        dispatcher_finish (&dispatcher, NULL);
    else
        dispatcher_pause (dispatcher);

    return piz_successful;
}
