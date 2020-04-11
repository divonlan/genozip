// ------------------------------------------------------------------
//   zip.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <math.h>
#include "genozip.h"
#include "profiler.h"
#include "vblock.h"
#include "buffer.h"
#include "file.h"
#include "zfile.h"
#include "txtfile.h"
#include "header.h"
#include "vblock.h"
#include "dispatcher.h"
#include "move_to_front.h"
#include "zip.h"
#include "base250.h"
#include "endianness.h"
#include "random_access_vcf.h"

static void zip_display_compression_ratio (Dispatcher dispatcher, bool is_last_file)
{
    const char *runtime = dispatcher_ellapsed_time (dispatcher, false);
    double z_bytes   = (double)z_file->disk_so_far;
    double vcf_bytes = (double)z_file->txt_data_so_far_concat;
    double ratio     = vcf_bytes / z_bytes;

    if (flag_concat) { // in concat, we don't show the compression ratio for files except for the last one

        static uint64_t txt_file_disk_size_concat = 0;
        static FileType source_file_type = UNKNOWN_FILE_TYPE;

        if (!txt_file_disk_size_concat) // first concat file
            source_file_type = txt_file->type;
        else if (source_file_type != txt_file->type) // heterogenous source file types
            source_file_type = UNKNOWN_FILE_TYPE;

        txt_file_disk_size_concat += txt_file->disk_size;

        fprintf (stderr, "Done (%s)                                     \n", runtime);

        if (is_last_file) {
            double ratio2 = (double)txt_file_disk_size_concat / z_bytes; // compression vs .gz/.bz2/.bcf/.xz... size

            if (txt_file->type == VCF || txt_file->type == SAM || ratio2 < 1)  // source file was plain txt or ratio2 is low (nothing to brag about)
                fprintf (stderr, "Time: %s, VCF compression ratio: %1.1f           \n", 
                            dispatcher_ellapsed_time (dispatcher, true), ratio);
            else
                fprintf (stderr, "Time: %s, VCF compression ratio: %1.1f - better than %s by a factor of %1.1f\n", 
                            dispatcher_ellapsed_time (dispatcher, true), ratio, 
                            source_file_type == UNKNOWN_FILE_TYPE ? "the input files" : file_exts[txt_file->type],
                            ratio2); // compression vs .gz/.bz2/.bcf/.xz... size
        }
    }
    else {
        double ratio2 = (double)txt_file->disk_size / z_bytes; // compression vs .gz/.bz2/.bcf/.xz... size
    
        if (txt_file->type == VCF || txt_file->type == SAM || ratio2 < 1)  // source file was plain txt or ratio2 is low (nothing to brag about)
            fprintf (stderr, "Done (%s, compression ratio: %1.1f)           \n", runtime, ratio);
        
        else // source was .vcf.gz or .vcf.bgz or .vcf.bz2
            fprintf (stderr, "Done (%s, VCF compression ratio: %1.1f - better than %s by a factor of %1.1f)\n", 
                     runtime, ratio, file_exts[txt_file->type], ratio2);
    }
}

// here we translate the mtf_i indeces creating during seg_* to their finally dictionary indeces in base-250.
// Note that the dictionary indeces have changed since segregate (which is why we needed this intermediate step)
// because: 1. the dictionary got integrated into the global one - some values might have already been in the global
// dictionary thanks to other threads working on other VBs ; 2. for the first VB, we sort the dictionary by frequency
void zip_generate_b250_section (VBlock *vb, MtfContext *ctx)
{
    buf_alloc (vb, &ctx->b250, ctx->mtf_i.len * MAX_BASE250_NUMERALS, // maximum length is if all entries are 5-numeral.
               1.1, "ctx->b250_buf", 0);

    bool show = flag_show_b250 || dict_id_printable (ctx->dict_id).num == dict_id_show_one_b250.num;

    if (show) 
        bufprintf (vb, &vb->show_b250_buf, "vb_i=%u %.*s: ", vb->vblock_i, DICT_ID_LEN, dict_id_printable(ctx->dict_id).id);

    int32_t prev = -1; 
    for (unsigned i=0; i < ctx->mtf_i.len; i++) {
        MtfNode *node = mtf_node (ctx, ((const uint32_t *)ctx->mtf_i.data)[i], NULL, NULL);

        uint32_t n            = node->word_index.n;
        unsigned num_numerals = base250_len (node->word_index.encoded.numerals);
        uint8_t *numerals     = node->word_index.encoded.numerals;
        
        bool one_up = (n == prev + 1) && (ctx->b250_section_type != SEC_VCF_GT_DATA) && (i > 0);

        if (one_up) // note: we can't do SEC_VCF_GT_DATA bc we can't PIZ it as many GT data types are in the same section 
            ((uint8_t *)ctx->b250.data)[ctx->b250.len++] = (uint8_t)BASE250_ONE_UP;

        else {
            memcpy (&ctx->b250.data[ctx->b250.len], numerals, num_numerals);
            ctx->b250.len += num_numerals;
        }

        if (show) {
            if (one_up) bufprintf (vb, &vb->show_b250_buf, "L%u:ONE_UP ", i)
            else        bufprintf (vb, &vb->show_b250_buf, "L%u:%u ", i, n)
        }

        prev = n;
    }
    if (show) {
        bufprintf (vb, &vb->show_b250_buf, "%s", "\n")
        fprintf (stderr, "%.*s", vb->show_b250_buf.len, vb->show_b250_buf.data);
        buf_free (&vb->show_b250_buf);
    }

    if (dict_id_printable (ctx->dict_id).num == dict_id_dump_one_b250.num) 
        fwrite (ctx->b250.data, 1, ctx->b250.len, stdout);
}


static void zip_output_processed_vb (VBlock *vb, Buffer *section_list_buf, bool update_vcf_file, bool is_z_data)
{
    START_TIMER;

    Buffer *data_buf = is_z_data ? &vb->z_data : &z_file->dict_data;

    if (section_list_buf) sections_list_concat (vb, section_list_buf);

    file_write (z_file, data_buf->data, data_buf->len);
    COPY_TIMER (vb->profile.write);

    z_file->disk_so_far += (int64_t)data_buf->len;
    
    if (is_z_data) {
        if (update_vcf_file) txt_file->num_lines += (int64_t)vb->num_lines; // lines in this VCF file
        z_file->num_lines                        += (int64_t)vb->num_lines; // lines in all concatenated VCF files in this z_file
        z_file->txt_data_so_far_single           += (int64_t)vb->vb_data_size;
        z_file->txt_data_so_far_concat           += (int64_t)vb->vb_data_size;
    }

    // update section stats
    for (unsigned sec_i=1; sec_i < NUM_SEC_TYPES; sec_i++) {
        if (update_vcf_file) txt_file->section_bytes[sec_i]  += vb->txt_section_bytes[sec_i];
        z_file->num_sections[sec_i]     += vb->z_num_sections[sec_i];
        z_file->section_bytes[sec_i]    += vb->z_section_bytes[sec_i];
        z_file->section_entries[sec_i]  += vb->z_section_entries[sec_i];
    }

    // this function holds the mutex and hence has a non-trival performance penalty. we call
    // it only if the user specifically requested --show-sections
    if (flag_show_sections) mtf_update_stats (vb);

    if (flag_show_headers && buf_is_allocated (&vb->show_headers_buf))
        buf_print (&vb->show_headers_buf, false);

    if (flag_show_threads) dispatcher_show_time ("Write genozip data done", -1, vb->vblock_i);
}

// write all the sections at the end of the file, after all VB stuff has been written
static void zip_write_global_area (const Md5Hash *single_component_md5)
{
    // output dictionaries to disk
    if (buf_is_allocated (&z_file->section_list_dict_buf)) // not allocated for vcf-header-only files
        zip_output_processed_vb (evb, &z_file->section_list_dict_buf, false, false);  
   
    // if this is VCF data, compress all random access records into evb->z_data
    if (z_file->data_type == DATA_TYPE_VCF) { 
        if (flag_show_index) random_access_show_index();
        
        BGEN_random_access(); // make ra_buf into big endian

        z_file->ra_buf.len *= random_access_sizeof_entry(); // change len to count bytes

        zfile_compress_section_data (evb, SEC_VCF_RANDOM_ACCESS, &z_file->ra_buf);
    }

    // compress genozip header (including its payload sectionlist and footer) into evb->z_data
    zfile_compress_genozip_header (single_component_md5);    

    // output to disk random access and genozip header sections to disk
    zip_output_processed_vb (evb, NULL, false, true);  
}

// this is the main dispatcher function. It first processes the VCF header, then proceeds to read 
// a variant block from the input file and send it off to a thread for computation. When the thread
// completes, this function proceeds to write the output to the output file. It can dispatch
// several threads in parallel.
void zip_dispatcher (const char *vcf_basename, unsigned max_threads, bool is_last_file)
{
    static unsigned last_vblock_i = 0; // used if we're concatenating files - the vblock_i will continue from one file to the next
    if (!flag_concat) last_vblock_i = 0; // reset if we're not concatenating

    // normally max_threads would be the number of cores available - we allow up to this number of compute threads, 
    // because the I/O thread is normally idling waiting for the disk, so not consuming a lot of CPU
    Dispatcher dispatcher = dispatcher_init (max_threads, last_vblock_i, false, is_last_file, vcf_basename);

    dict_id_initialize (z_file->data_type);

    unsigned txt_line_i = 1; // the next line to be read (first line = 1)
    
    // read the vcf header, assign the global variables, and write the compressed header to the GENOZIP file
    off64_t txt_header_header_pos = z_file->disk_so_far;
    bool success = header_txt_to_genozip (&txt_line_i);
    if (!success) goto finish;

    mtf_initialize_mutex();

    uint32_t max_lines_per_vb=0;

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. In there is a new variant block avaialble, and a compute thread available to take it - dispatch it
    // 2. If there is no new variant block available, but input is not exhausted yet - read one
    // 3. Wait for the first thread (by sequential order) to complete the compute and output the results
    VBlock *next_vb;
    do {
        next_vb = dispatcher_get_next_vb (dispatcher);
        bool has_vb_ready_to_compute = next_vb && next_vb->ready_to_dispatch;
        bool has_free_thread = dispatcher_has_free_thread (dispatcher);

        // PRIORITY 1: is there a block available and a compute thread available? in that case dispatch it
        if (has_vb_ready_to_compute && has_free_thread) {
            DispatcherFuncType compress_funcs[NUM_DATATYPES] = { zip_vcf_compress_one_vb, zip_sam_compress_one_vb };
            dispatcher_compute (dispatcher, compress_funcs[z_file->data_type]);
        }
        
        // PRIORITY 2: output completed vbs, so they can be released and re-used
        else if (dispatcher_has_processed_vb (dispatcher, NULL) ||  // case 1: there is a VB who's compute processing is completed
                 (has_vb_ready_to_compute && !has_free_thread)) {   // case 2: a VB ready to dispatch but all compute threads are occupied. wait here for one to complete
           
            VBlock *processed_vb = dispatcher_get_processed_vb (dispatcher, NULL); // this will block until one is available
            if (!processed_vb) continue; // no running compute threads 

            // update z_data in memory (its not written to disk yet)
            zfile_update_compressed_vb_header (processed_vb, txt_line_i);

            max_lines_per_vb = MAX (max_lines_per_vb, processed_vb->num_lines);
            txt_line_i += processed_vb->num_lines;

            zip_output_processed_vb (processed_vb, &processed_vb->section_list_buf, true, true);

            dispatcher_finalize_one_vb (dispatcher);
        }        
        
        // PRIORITY 3: If there is no variant block available to compute or to output, but input is not exhausted yet - read one
        else if (!next_vb && !dispatcher_is_input_exhausted (dispatcher)) {

            next_vb = dispatcher_generate_next_vb (dispatcher, 0);

            if (flag_show_threads) dispatcher_show_time ("Read input data", -1, next_vb->vblock_i);
            
            txtfile_read_variant_block (next_vb);

            if (flag_show_threads) dispatcher_show_time ("Read input data done", -1, next_vb->vblock_i);

            if (next_vb->txt_data.len)  // we found some data 
                next_vb->ready_to_dispatch = true;
            else {
                // this vb has no data
                dispatcher_input_exhausted (dispatcher);
                
                // update to the conclusive size. it might have been 0 (eg STDIN if HTTP) or an estimate (if compressed)
                txt_file->txt_data_size_single = txt_file->txt_data_so_far_single; 

                dispatcher_finalize_one_vb (dispatcher); 
            }
        }
    } while (!dispatcher_is_done (dispatcher));

    // go back and update some fields in the vcf header's section header and genozip header -
    // only if we can go back - i.e. is a normal file, not redirected
    Md5Hash single_component_md5;
    if (z_file && !z_file->redirected && txt_header_header_pos >= 0) 
        success = zfile_update_txt_header_section_header (txt_header_header_pos, max_lines_per_vb, &single_component_md5);

    // if this a non-concatenated file, or the last vcf component of a concatenated file - write the genozip header, random access and dictionaries
    if (is_last_file || !flag_concat) zip_write_global_area (&single_component_md5);

    if (!flag_quiet) zip_display_compression_ratio (dispatcher, is_last_file);

finish:
    z_file->disk_size              = z_file->disk_so_far;
    z_file->txt_data_so_far_single = 0;
    evb->z_data.len                = 0;
    evb->z_next_header_i           = 0;
    memset (&z_file->md5_ctx_single, 0, sizeof (Md5Context));

    dispatcher_finish (&dispatcher, &last_vblock_i);
}
