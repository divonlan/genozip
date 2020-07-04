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
#include "vblock.h"
#include "dispatcher.h"
#include "context.h"
#include "zip.h"
#include "seg.h"
#include "base250.h"
#include "endianness.h"
#include "random_access.h"
#include "dict_id.h"
#include "reference.h"

static void zip_display_compression_ratio (Dispatcher dispatcher, bool is_final_component)
{
    if (flag_quiet) return; 

    const char *runtime = dispatcher_ellapsed_time (dispatcher, false);
    double z_bytes   = (double)z_file->disk_so_far;
    double txt_bytes = (double)z_file->txt_data_so_far_concat;
    double ratio     = txt_bytes / z_bytes;

    // in concat, or when we store the refernce, we don't show the compression ratio for files except for the last one
    if (flag_concat || flag_reference == REF_INTERNAL || flag_reference == REF_EXT_STORE) { 

        static uint64_t txt_file_disk_size_concat = 0;
        static FileType source_file_type = UNKNOWN_FILE_TYPE;

        if (!txt_file_disk_size_concat) // first concat file
            source_file_type = txt_file->type;
        else if (source_file_type != txt_file->type) // heterogenous source file types
            source_file_type = UNKNOWN_FILE_TYPE;

        txt_file_disk_size_concat += txt_file->disk_size;

        if (!flag_quiet)
            fprintf (stderr, "Done (%s)                                     \n", runtime);

        if (is_final_component) {
            double ratio2 = (double)txt_file_disk_size_concat / z_bytes; // compression vs .gz/.bz2/.bcf/.xz... size

            if (txt_file->comp_alg == COMP_NONE || ratio2 < 1)  // source file was plain txt or ratio2 is low (nothing to brag about)
                fprintf (stderr, "Time: %s, %s compression ratio: %1.1f           \n", 
                            dispatcher_ellapsed_time (dispatcher, true), dt_name (z_file->data_type), ratio);
            else
                fprintf (stderr, "Time: %s, %s compression ratio: %1.1f - better than %s by a factor of %1.1f\n", 
                            dispatcher_ellapsed_time (dispatcher, true), dt_name (z_file->data_type), ratio, 
                            source_file_type == UNKNOWN_FILE_TYPE ? "the input files" : file_exts[txt_file->type],
                            ratio2); // compression vs .gz/.bz2/.bcf/.xz... size
        }
    }
    else {
        double ratio2 = (double)txt_file->disk_size / z_bytes; // compression vs .gz/.bz2/.bcf/.xz... size
    
        if (txt_file->comp_alg == COMP_NONE || ratio2 < 1)  // source file was plain txt or ratio2 is low (nothing to brag about)
            fprintf (stderr, "Done (%s, compression ratio: %1.1f)           \n", runtime, ratio);
        
        else // source was compressed
            fprintf (stderr, "Done (%s, %s compression ratio: %1.1f - better than %s by a factor of %1.1f)\n", 
                     runtime, dt_name (z_file->data_type), ratio, file_exts[txt_file->type], ratio2);
    }
}

// after segging - if any context appears to contain only singleton snips (eg a unique ID),
// we move it to local instead of needlessly cluttering the global dictionary
static void zip_handle_unique_words_ctxs (VBlock *vb)
{
    for (int did_i=0 ; did_i < vb->num_dict_ids ; did_i++) {
        Context *ctx = &vb->contexts[did_i];
    
        if (!ctx->mtf.len || ctx->mtf.len != ctx->mtf_i.len) continue; // check that all words are unique (and new to this vb)
        if (vb->data_type == DT_VCF && dict_id_is_vcf_format_sf (ctx->dict_id)) continue; // this doesn't work for FORMAT fields
        if (ctx->mtf.len < vb->lines.len / 5)   continue; // don't bother if this is a rare field less than 20% of the lines
        if (buf_is_allocated (&ctx->local))     continue; // skip if we are already using local to optimize in some other way

        // don't move to local if its on the list of special dict_ids that are always in dict (because local is used for something else - eg pos or id data)
        if ((ctx->flags & CTX_FL_NO_STONS) || ctx->ltype != CTX_LT_TEXT) continue; // NO_STONS is implicit if ctx isn't text

        buf_move (vb, &ctx->local, vb, &ctx->dict);
        buf_free (&ctx->mtf);
        buf_free (&ctx->mtf_i);
    }
}

// generate & write b250 data for all primary fields of this data type
void zip_generate_and_compress_ctxs (VBlock *vb)
{
    // generate & write b250 data for all primary fields
    for (int did_i=0 ; did_i < vb->num_dict_ids ; did_i++) {
        Context *ctx = &vb->contexts[did_i];

        if (ctx->mtf_i.len && 
            (vb->data_type != DT_VCF || !dict_id_is_vcf_format_sf (ctx->dict_id))) { // skip VCF FORMAT subfields, as they get compressed into SEC_GT_DATA instead
            
            zip_generate_b250_section (vb, ctx);
            zfile_compress_b250_data (vb, ctx, COMP_BZ2);
        }

        if (ctx->local.len || ctx->ltype == CTX_LT_SEQ_BITMAP) // bitmaps are always written, even if empty
            zfile_compress_local_data (vb, ctx);
    }
}

// here we translate the mtf_i indeces creating during seg_* to their finally dictionary indeces in base-250.
// Note that the dictionary indeces have changed since segregate (which is why we needed this intermediate step)
// because: 1. the dictionary got integrated into the global one - some values might have already been in the global
// dictionary thanks to other threads working on other VBs ; 2. for the first VB, we sort the dictionary by frequency
void zip_generate_b250_section (VBlock *vb, Context *ctx)
{
    ASSERT (ctx->b250.len==0, "Error in zip_generate_b250_section: ctx->mtf_i is not empty. Dict=%s", ctx->name);

    buf_alloc (vb, &ctx->b250, ctx->mtf_i.len * MAX_BASE250_NUMERALS, // maximum length is if all entries are 4-numeral.
               1.1, "ctx->b250_buf", 0);

    bool show = flag_show_b250 || dict_id_printable (ctx->dict_id).num == dict_id_show_one_b250.num;

    if (show) 
        bufprintf (vb, &vb->show_b250_buf, "vb_i=%u %s: ", vb->vblock_i, ctx->name);

    int32_t prev = -1; 
    for (unsigned i=0; i < ctx->mtf_i.len; i++) {

        uint32_t node_index = *ENT(uint32_t, ctx->mtf_i, i);

        if (node_index <= WORD_INDEX_MAX_INDEX) { // normal index

            MtfNode *node = mtf_node_vb (ctx, node_index, NULL, NULL);

            uint32_t n            = node->word_index.n;
            unsigned num_numerals = base250_len (node->word_index.encoded.numerals);
            uint8_t *numerals     = node->word_index.encoded.numerals;
            
            bool one_up = (n == prev + 1) && (ctx->dict_id.num != dict_id_fields[VCF_GT]) && (i > 0);

            if (one_up) { // note: we can't do SEC_VCF_GT_DATA bc we can't PIZ it as many GT data types are in the same section 
                NEXTENT(uint8_t, ctx->b250) = (uint8_t)BASE250_ONE_UP;
                if (show) bufprintf (vb, &vb->show_b250_buf, "L%u:ONE_UP ", i)
            }

            else {
                memcpy (AFTERENT (char, ctx->b250), numerals, num_numerals);
                ctx->b250.len += num_numerals;
                if (show) bufprintf (vb, &vb->show_b250_buf, "L%u:%u ", i, n)
            }
            prev = n;
        }

        else if (node_index == WORD_INDEX_MISSING_SF) {
            NEXTENT(uint8_t, ctx->b250) = (uint8_t)BASE250_MISSING_SF;
            prev = node_index;
            if (show) bufprintf (vb, &vb->show_b250_buf, "L%u:MISSING ", i)
        }
        
        else if (node_index == WORD_INDEX_EMPTY_SF) {
            NEXTENT(uint8_t, ctx->b250) = (uint8_t)BASE250_EMPTY_SF;
            prev = node_index;
            if (show) bufprintf (vb, &vb->show_b250_buf, "L%u:EMPTY ", i)
        }

        else ABORT ("Error in zip_generate_b250_section: invalid node_index=%u", node_index);        
    }

    if (show) {
        bufprintf (vb, &vb->show_b250_buf, "%s", "\n")
        fprintf (stderr, "%.*s", (uint32_t)vb->show_b250_buf.len, vb->show_b250_buf.data);
        buf_free (&vb->show_b250_buf);
    }

    if (dict_id_printable (ctx->dict_id).num == dict_id_dump_one_b250.num) 
        fwrite (ctx->b250.data, 1, ctx->b250.len, stdout);
}


void zip_output_processed_vb (VBlock *vb, Buffer *section_list_buf, bool update_txt_file, ProcessedDataType pd_type)
{
    START_TIMER;

    Buffer *data_buf = (pd_type == PD_DICT_DATA) ? &z_file->dict_data : &vb->z_data;

    if (section_list_buf) sections_list_concat (vb, section_list_buf);

    file_write (z_file, data_buf->data, data_buf->len);
    COPY_TIMER (vb->profile.write);

    z_file->disk_so_far += (int64_t)data_buf->len;
    
    if (pd_type == PD_VBLOCK_DATA) {
        if (update_txt_file) txt_file->num_lines += (int64_t)vb->lines.len; // lines in this txt file
        z_file->num_lines                        += (int64_t)vb->lines.len; // lines in all concatenated files in this z_file
        z_file->txt_data_so_far_single           += (int64_t)vb->vb_data_size;
        z_file->txt_data_so_far_concat           += (int64_t)vb->vb_data_size;
    }

    // this function holds the mutex and hence has a non-trival performance penalty. we call
    // it only if the user specifically requested --show-sections
    if (flag_show_sections) mtf_update_stats (vb);

    if (flag_show_headers && buf_is_allocated (&vb->show_headers_buf))
        buf_print (&vb->show_headers_buf, false);

    if (flag_show_threads) dispatcher_show_time ("Write genozip data done", -1, vb->vblock_i);
}

// write all the sections at the end of the file, after all VB stuff has been written
static void zip_write_global_area (Dispatcher dispatcher, const Md5Hash *single_component_md5)
{
    // output dictionaries (inc. aliases) to disk - they are in the "processed" data of evb
    if (buf_is_allocated (&z_file->section_list_dict_buf)) // not allocated for vcf-header-only files
        zip_output_processed_vb (evb, &z_file->section_list_dict_buf, false, PD_DICT_DATA);  

    // if we're making a reference, we need the RA data to populate the reference section chrome/first/last_pos ahead of ref_compress_ref
    if (DTPZ(has_random_access)) 
        random_access_finalize_entries();

    // output reference, if needed
    if (flag_reference == REF_INTERNAL || flag_reference == REF_EXT_STORE || flag_make_reference) {
        ref_compress_ref();
        zip_display_compression_ratio (dispatcher, true); // Done for reference + final compression ratio calculation
    }

    // add dict_id aliases list, if we have one
    Buffer *dict_id_aliases_buf = dict_id_create_aliases_buf();
    if (dict_id_aliases_buf->len) zfile_compress_section_data (evb, SEC_DICT_ID_ALIASES, dict_id_aliases_buf);

    // if this data has random access (i.e. it has chrom and pos), compress all random access records into evb->z_data
    if (DTPZ(has_random_access)) {

        // sort RA, update entries that don't yet have a chrom_index
        random_access_finalize_entries();

        if (flag_show_index) random_access_show_index(true);
        
        BGEN_random_access(); // make ra_buf into big endian

        z_file->ra_buf.len *= random_access_sizeof_entry(); // change len to count bytes

        zfile_compress_section_data_alg (evb, SEC_RANDOM_ACCESS, &z_file->ra_buf, 0,0, COMP_LZMA); // ra data compresses better with LZMA than BZLIB
    }

    // compress genozip header (including its payload sectionlist and footer) into evb->z_data
    zfile_compress_genozip_header (single_component_md5);    

    // output to disk random access and genozip header sections to disk
    zip_output_processed_vb (evb, NULL, false, PD_VBLOCK_DATA);  
}

// entry point for compute thread
static void zip_compress_one_vb (VBlock *vb)
{
    START_TIMER; 

    // if we're vb_i=1 lock, and unlock only when we're done merging. all other vbs need  
    // to wait for our merge. that is because our dictionaries are sorted 
    if (vb->vblock_i == 1) mtf_vb_1_lock(vb); 

    // allocate memory for the final compressed data of this vb. allocate 20% of the
    // vb size on the original file - this is normally enough. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, vb->vb_data_size / 5, 1.2, "z_data", 0);

    // clone global dictionaries while granted exclusive access to the global dictionaries
    mtf_clone_ctx (vb);

    // split each line in this variant block to its components
    seg_all_data_lines (vb);

    // identify dictionaries that contain almost only unique words (eg a unique id) and move the data from dict to local
    zip_handle_unique_words_ctxs (vb);

    // for the first vb only - sort dictionaries so that the most frequent entries get single digit
    // base-250 indices. This can be done only before any dictionary is written to disk, but likely
    // beneficial to all vbs as they are likely to more-or-less have the same frequent entries
    if (vb->vblock_i == 1) {
        mtf_sort_dictionaries_vb_1(vb);

        // for a file that's compressed (and hence we don't know the size of content a-priori) AND we know the
        // compressed file size (eg a local gz/bz2 file or a compressed file on a ftp server) - we now estimate the 
        // txt_data_size_single that will be used for the global_hash and the progress indicator
        txtfile_estimate_txt_data_size (vb);
    }

    if (vb->data_type == DT_VCF)
        vcf_zip_generate_ht_gt_compress_vb_header (vb);
    else
        zfile_compress_generic_vb_header (vb); // vblock header

    // merge new words added in this vb into the z_file.contexts, ahead of zip_generate_b250_section() and
    // vcf_zip_generate_genotype_one_section(). writing indices based on the merged dictionaries. dictionaries are compressed. 
    // all this is done while holding exclusive access to the z_file dictionaries.
    mtf_merge_in_vb_ctx(vb);

    // merge in random access - IF it is used
    if (DTP(has_random_access)) 
        random_access_merge_in_vb (vb);

    // optimize FORMAT/GL data in local (we have already optimized the data in dict elsewhere)
    if (vb->data_type == DT_VCF) {
        uint8_t gl_did_i =  mtf_get_existing_did_i (vb, (DictId)dict_id_FORMAT_GL);
        if (gl_did_i != DID_I_NONE) 
            gl_optimize_local (vb, &vb->contexts[gl_did_i].local);
    }

    // generate & compress b250 and local data for all ctxs (for reference files we don't output VBs)
    if (!flag_make_reference)
        zip_generate_and_compress_ctxs (vb);

    // compress data-type specific sections
    if (DTP(compress)) DTP(compress)(vb);

    // tell dispatcher this thread is done and can be joined.
    // thread safety: this isn't protected by a mutex as it will just be false until it at some point turns to true
    // this this operation needn't be atomic, but it likely is anyway
    vb->is_processed = true; 

    COPY_TIMER (vb->profile.compute);
}

// this is the main dispatcher function. It first processes the txt header, then proceeds to read 
// a variant block from the input file and send it off to a thread for computation. When the thread
// completes, this function proceeds to write the output to the output file. It can dispatch
// several threads in parallel.
void zip_dispatcher (const char *txt_basename, bool is_last_file)
{
    static DataType last_data_type = DT_NONE;
    static unsigned last_vblock_i = 0; // used if we're concatenating files - the vblock_i will continue from one file to the next
    if (!flag_concat) last_vblock_i = 0; // reset if we're not concatenating

    // we cannot concatenate files of different type
    ASSERT (!flag_concat || txt_file->data_type == last_data_type || last_data_type == DT_NONE, 
            "%s: cannot concatenate %s because it is a %s file, whereas the previous file was a %s",
             global_cmd, txt_name, dt_name (txt_file->data_type), dt_name (last_data_type));
    last_data_type =  txt_file->data_type;

    // normally global_max_threads would be the number of cores available - we allow up to this number of compute threads, 
    // because the I/O thread is normally idling waiting for the disk, so not consuming a lot of CPU
    Dispatcher dispatcher = dispatcher_init (global_max_threads, last_vblock_i, false, is_last_file, txt_basename);

    dict_id_initialize(z_file->data_type);

    if (DTPZ(zip_initialize)) DTPZ(zip_initialize)();

    unsigned txt_line_i = 1; // the next line to be read (first line = 1)
    
    // read the txt header, assign the global variables, and write the compressed header to the GENOZIP file
    off64_t txt_header_header_pos = z_file->disk_so_far;
    bool success = txtfile_header_to_genozip (&txt_line_i);
    if (!success) goto finish;

    mtf_initialize_for_zip();

    if (z_file->data_type == DT_VCF) vcf_zip_initialize();

    uint32_t max_lines_per_vb=0;

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. If there is a new variant block avaialble, and a compute thread available to take it - dispatch it
    // 2. If there is no new variant block available, but input is not exhausted yet - read one
    // 3. Wait for the first thread (by sequential order) to complete the compute and output the results
    VBlock *next_vb;
    do {
        next_vb = dispatcher_get_next_vb (dispatcher);
        bool has_vb_ready_to_compute = next_vb && next_vb->ready_to_dispatch;
        bool has_free_thread = dispatcher_has_free_thread (dispatcher);

        // PRIORITY 1: is there a block available and a compute thread available? in that case dispatch it
        if (has_vb_ready_to_compute && has_free_thread) 
            dispatcher_compute (dispatcher, zip_compress_one_vb);

        // PRIORITY 2: output completed vbs, so they can be released and re-used
        else if (dispatcher_has_processed_vb (dispatcher, NULL) ||  // case 1: there is a VB who's compute processing is completed
                 (has_vb_ready_to_compute && !has_free_thread)) {   // case 2: a VB ready to dispatch but all compute threads are occupied. wait here for one to complete
           
            VBlock *processed_vb = dispatcher_get_processed_vb (dispatcher, NULL); // this will block until one is available
            if (!processed_vb) continue; // no running compute threads 

            // update z_data in memory (its not written to disk yet)
            DTPZ(update_header)(processed_vb, txt_line_i); 

            max_lines_per_vb = MAX (max_lines_per_vb, processed_vb->lines.len);
            txt_line_i += (uint32_t)processed_vb->lines.len;

            if (!flag_make_reference)
                zip_output_processed_vb (processed_vb, &processed_vb->section_list_buf, true, PD_VBLOCK_DATA);
            else 
                z_file->txt_data_so_far_concat += (int64_t)processed_vb->vb_data_size;
                
            z_file->num_vbs++;
            
            dispatcher_finalize_one_vb (dispatcher);
        }        
        
        // PRIORITY 3: If there is no variant block available to compute or to output, but input is not exhausted yet - read one
        else if (!next_vb && !dispatcher_is_input_exhausted (dispatcher)) {

            next_vb = dispatcher_generate_next_vb (dispatcher, 0);
            if (flag_show_threads) dispatcher_show_time ("Read input data", -1, next_vb->vblock_i);
            
            txtfile_read_vblock (next_vb);

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
        else  // nothing for us to do right now, just wait
            usleep (100000); // 100ms

    } while (!dispatcher_is_done (dispatcher));

    // go back and update some fields in the txt header's section header and genozip header -
    // only if we can go back - i.e. is a normal file, not redirected
    Md5Hash single_component_md5;
    if (z_file && !z_file->redirected && txt_header_header_pos >= 0) 
        success = zfile_update_txt_header_section_header (txt_header_header_pos, max_lines_per_vb, &single_component_md5);

    // if this a non-concatenated file, or the last component of a concatenated file - write the genozip header, random access and dictionaries
    if (is_last_file || !flag_concat) {
        
        if (flag_reference == REF_INTERNAL || flag_reference == REF_EXT_STORE) 
            zip_display_compression_ratio (dispatcher, false); // Done for the file (if reference is to be written by zip_write_global_area)

        zip_write_global_area (dispatcher, &single_component_md5);

        if (flag_reference == REF_NONE || flag_reference == REF_EXTERNAL) 
            zip_display_compression_ratio (dispatcher, true); // Done for the file (if no stored reference - ratio includes global area in compression ratio)
    }
    else  // non-last file in concat mode
        zip_display_compression_ratio (dispatcher, false);

finish:
    z_file->disk_size              = z_file->disk_so_far;
    z_file->txt_data_so_far_single = 0;
    evb->z_data.len                = 0;
    evb->z_next_header_i           = 0;
    memset (&z_file->md5_ctx_single, 0, sizeof (Md5Context));

    dispatcher_finish (&dispatcher, &last_vblock_i);
}
