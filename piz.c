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

// Compute threads: decode the delta-encoded value of the POS field, and returns the new last_pos
int32_t piz_decode_pos (int32_t last_pos, const char *delta_snip, unsigned delta_snip_len,
                        char *pos_str, unsigned *pos_len) // out
{
    int32_t delta=0;

    // we parse the string ourselves - this is hugely faster than sscanf.
    unsigned negative = *delta_snip == '-'; // 1 or 0

    const char *s; for (s=(delta_snip + negative); *s != '\t' ; s++)
        delta = delta*10 + (*s - '0');

    if (negative) delta = -delta;

    last_pos += delta;
    
    ASSERT (last_pos >= 0, "Error: last_pos=%d is negative", last_pos);

    // create number string without calling slow sprintf

    // create reverse string
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

// Called by PIZ I/O thread: read all the sections at the end of the file, before starting to process VBs
static int16_t piz_read_global_area (Md5Hash *original_file_digest) // out
{
    int16_t data_type = zfile_read_genozip_header (original_file_digest);
    if (data_type == DATA_TYPE_VCF_V1 || data_type == EOF) return data_type;

    dict_id_initialize();
    
    // if the user wants to see only the header, we can skip the dictionaries, regions and random access
    if (!flag_header_only) {
        
        // read random access, but only if we are going to need it
        if (flag_regions || flag_show_index) {
            zfile_read_all_dictionaries (0, DICTREAD_CHROM_ONLY); // read all CHROM/RNAME dictionaries - needed for regions_make_chregs()

            // update chrom node indeces using the CHROM dictionary, for the user-specified regions (in case -r/-R were specified)
            regions_make_chregs (chrom_did_i_by_data_type[data_type]);

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

static void enforce_v1_limitations (bool is_first_vcf_component)
{
    #define ENFORCE(flag,lflag) ASSERT (!(flag), "Error: %s option is not supported because %s compressed with genozip version 1", (lflag), z_name);
    
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
    ENFORCE(flag_strip, "--flag_strip");
}

// returns true is successfully outputted a vcf file
bool piz_dispatcher (const char *z_basename, unsigned max_threads, 
                     bool is_first_vcf_component, bool is_last_file)
{
    // static dispatcher - with flag_split, we use the same dispatcher when unzipping components
    static Dispatcher dispatcher = NULL;
    bool piz_successful = false;
    SectionListEntry *sl_ent = NULL;
    
    if (flag_split && !sections_has_more_vcf_components()) return false; // no more components

    if (!dispatcher) 
        dispatcher = dispatcher_init (max_threads, 0, flag_test, is_last_file, z_basename);
    
    // read genozip header
    Md5Hash original_file_digest;

    // read genozip header, dictionaries etc and set the data type when reading the first vcf component of in case of --split, 
    static int16_t data_type = EOF; 
    if (is_first_vcf_component) {
        data_type = piz_read_global_area (&original_file_digest);

        if (data_type != DATA_TYPE_VCF_V1)  // genozip v2+ - move cursor past first txt header
            ASSERT (sections_get_next_header_type(&sl_ent, NULL, NULL) == SEC_TXT_HEADER, "Error: unable to find TXT Header data in %s", z_name);

        ASSERT (!flag_test || !md5_is_zero (original_file_digest), 
                "Error testing %s: --test cannot be used with this file, as it was not compressed with --md5", z_name);
    }

    // case: we couldn't open the file because we didn't know what type it is - open it now
    if (!txt_file->file) file_open_txt (txt_file);

    if (data_type == EOF) goto finish;

    if (z_file->genozip_version < 2) enforce_v1_limitations (is_first_vcf_component); // genozip_version will be 0 for v1, bc we haven't read the vcf header yet

    // read and write VCF header. in split mode this also opens txt_file
    piz_successful = (data_type != DATA_TYPE_VCF_V1) ? header_genozip_to_txt (&original_file_digest)
                                             : v1_header_genozip_to_vcf (&original_file_digest);
    
    ASSERT (piz_successful || !is_first_vcf_component, "Error: failed to read VCF header in %s", z_name);
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

            bool compute = false;
            if (z_file->genozip_version > 1) {
                
                bool skipped_vb;
                static Buffer region_ra_intersection_matrix = EMPTY_BUFFER; // we will move the data to the VB when we get it
                SectionType header_type = sections_get_next_header_type (&sl_ent, &skipped_vb, &region_ra_intersection_matrix);
                switch (header_type) {
                    case SEC_SAM_VB_HEADER: 
                    case SEC_VCF_VB_HEADER: {

                        // if we skipped VBs or we skipped the sample sections in the last vb (flag_drop_genotypes), we need to seek forward 
                        if (skipped_vb || flag_drop_genotypes) file_seek (z_file, sl_ent->offset, SEEK_SET, false); 

                        VBlock *next_vb = dispatcher_generate_next_vb (dispatcher, sl_ent->vblock_i);
                        
                        if (region_ra_intersection_matrix.data) {
                            buf_copy (next_vb, &next_vb->region_ra_intersection_matrix, &region_ra_intersection_matrix, 0,0,0, "region_ra_intersection_matrix", next_vb->vblock_i);
                            buf_free (&region_ra_intersection_matrix); // note: copy & free rather than move - so memory blocks are preserved for VB re-use
                        }
                        
                        if      (header_type == SEC_VCF_VB_HEADER) zfile_vcf_read_one_vb ((VBlockVCFP)next_vb); 
                        else if (header_type == SEC_SAM_VB_HEADER) zfile_sam_read_one_vb ((VBlockSAMP)next_vb); 

                        compute = true;
                        break;
                    }

                    case SEC_TXT_HEADER: // 2nd+ vcf header of a concatenated file
                        if (!flag_split) {
                            header_genozip_to_txt (NULL); // skip 2nd+ vcf header if concatenating
                            continue;
                        }
                        break; // eof if splitting
                    
                    case SEC_EOF: 
                        break; 
                    
                    default: ABORT ("Error in piz_dispatcher: unexpected section_type=%s", st_name (header_type));
                }
            }
            else compute = v1_zfile_vcf_read_one_vb ((VBlockVCF *)dispatcher_generate_next_vb (dispatcher, 0));  // genozip v1
            
            if (compute) {
                static DispatcherFuncType uncompress_funcs[NUM_DATATYPES] =  
                    { piz_vcf_uncompress_one_vb, piz_sam_uncompress_one_vb };

                dispatcher_compute (dispatcher, uncompress_funcs[z_file->data_type]);
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
    
            switch (z_file->data_type) {
                case DATA_TYPE_VCF: txtfile_write_one_vblock_vcf ((VBlockVCFP)processed_vb); break;
                case DATA_TYPE_SAM: txtfile_write_one_vblock_sam ((VBlockSAMP)processed_vb); break;
                default:            ABORT ("Error in piz_dispatcher: invalid data_type=%d", z_file->data_type);
            }
            
            z_file->txt_data_so_far_single += processed_vb->vb_data_size; 

            dispatcher_finalize_one_vb (dispatcher);
        }

    } while (!dispatcher_is_done (dispatcher));

    // verify file integrity, if the genounzip compress was run with --md5 or --test
    if (flag_md5) {
        Md5Hash decompressed_file_digest;
        md5_finalize (&txt_file->md5_ctx_concat, &decompressed_file_digest); // z_file might be a concatenation - this is the MD5 of the entire concatenation

        if (md5_is_zero (original_file_digest)) 
            fprintf (stderr, "MD5 = %s Note: unable to compare this to the original file as file was not originally compressed with --md5\n", md5_display (&decompressed_file_digest, false));
        
        else if (md5_is_equal (decompressed_file_digest, original_file_digest)) {

            if (flag_test && !flag_quiet) fprintf (stderr, "Success          \b\b\b\b\b\b\b\b\b\b\n");

            if (flag_md5) fprintf (stderr, "MD5 = %s verified as identical to the original %s\n", 
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
