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
#include "refhash.h"
#include "progress.h"
#include "mutex.h"
#include "fastq.h"
#include "stats.h"
#include "codec.h"
#include "compressor.h"
#include "strings.h"
#include "bgzf.h"

static void zip_display_compression_ratio (Dispatcher dispatcher, Md5Hash md5, bool is_final_component)
{
    double z_bytes   = (double)z_file->disk_so_far;
    double txt_bytes = (double)z_file->txt_data_so_far_bind;
    double ratio     = txt_bytes / z_bytes;
    double ratio2    = -1;

    // in bind mode, we don't show compression ratio for files except for the last one
    if (flag.bind) { 

        static uint64_t txt_file_disk_size_bind = 0;
        static FileType source_file_type = UNKNOWN_FILE_TYPE;

        if (!txt_file_disk_size_bind) // first bound file
            source_file_type = txt_file->type;
        else if (source_file_type != txt_file->type) // heterogenous source file types
            source_file_type = UNKNOWN_FILE_TYPE;

        txt_file_disk_size_bind += txt_file->disk_size;

        if (is_final_component) 
            ratio2 = (double)txt_file_disk_size_bind / z_bytes; // compression vs .gz/.bz2/.bcf/.xz... size
        else 
            progress_finalize_component_time ("Done", md5);
    }
    else 
        ratio2 = (double)txt_file->disk_size / z_bytes; // compression vs .gz/.bz2/.bcf/.xz... size
    
    // when making a reference, we don't care about the compression
    if (flag.make_reference)
        progress_finalize_component_time ("Done", md5);

    // when compressing BAM report only ratio2 (compresstion against bgzipped BAM - we don't care about the underlying non-bgzipped BAM)
    else if (z_file->data_type == DT_BAM) 
            progress_finalize_component_time_ratio (dt_name (z_file->data_type), ratio2, md5);

    else if (ratio2 >= 0) {
        if (txt_file->codec == CODEC_NONE || ratio2 < 1.05)  // source file was plain txt or ratio2 is low (nothing to brag about)
            progress_finalize_component_time_ratio (dt_name (z_file->data_type), ratio, md5);
        
        else // source was compressed
            progress_finalize_component_time_ratio_better (dt_name (z_file->data_type), ratio, file_exts[txt_file->type], ratio2, md5);
    }
}

// here we translate the node_i indices creating during seg_* to their finally dictionary indices in base-250.
// Note that the dictionary indices have changed since segregate (which is why we needed this intermediate step)
// because: 1. the dictionary got integrated into the global one - some values might have already been in the global
// dictionary thanks to other threads working on other VBs ; 2. for the first VB, we sort the dictionary by frequency
// returns true if section should be dropped
static bool zip_generate_b250_section (VBlock *vb, Context *ctx, uint32_t sample_size)
{
    ASSERT (ctx->b250.len==0, "Error in zip_generate_b250_section: ctx->node_i is not empty. Dict=%s", ctx->name);

    buf_alloc (vb, &ctx->b250, ctx->node_i.len * MAX_BASE250_NUMERALS, // maximum length is if all entries are 4-numeral.
               1.1, "ctx->b250_buf");

    bool show = (flag.show_b250 || dict_id_printable (ctx->dict_id).num == flag.dict_id_show_one_b250.num) && !sample_size;

    if (show) 
        bufprintf (vb, &vb->show_b250_buf, "vb_i=%u %s: ", vb->vblock_i, ctx->name);

    // calculate number of node_i words to be generated - normally all of them, except if we're just sampling in zip_assign_best_codec
    uint32_t num_words = (uint32_t)ctx->node_i.len;
    if (sample_size && num_words > sample_size / sizeof (uint32_t))
        num_words = sample_size / sizeof (uint32_t);

    WordIndex first_node_index = *ENT (WordIndex, ctx->node_i, 0);
    bool all_the_same = true; // are all the node_index of this context the same in this VB

    WordIndex prev = WORD_INDEX_NONE; 
    for (uint32_t i=0; i < num_words; i++) {

        WordIndex node_index = *ENT(WordIndex, ctx->node_i, i);

        if (node_index != first_node_index) all_the_same = false;

        if (node_index >= 0) { // normal index

            CtxNode *node = ctx_node_vb (ctx, node_index, NULL, NULL);

            WordIndex n           = node->word_index.n;
            unsigned num_numerals = base250_len (node->word_index.encoded.numerals);
            uint8_t *numerals     = node->word_index.encoded.numerals;
            
            bool one_up = (n == prev + 1) && (i > 0);

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

    // if all the node_index of this context are the same in this VB, we store just one, and set a flag
    if (!sample_size) { // this logic is not to be applied if we're just sampling 
        if (all_the_same) {

            // if the entire section is word_index = 0 we can drop it, unless:
            // 1) it has local (bc if no-b250/dict-from-prev-vb/no-local piz can't distiguish between seg_id-with-b250-all-the-same-word-index-0 vs all-singleton-pushed-to-local)
            // 2) it has flags that need to be passed to piz (we can get rid of this limitation - bug 224), or
            // 3) the one snip is SELF_DELTA
            if (base250_decode ((const uint8_t **)&ctx->b250.data, false) == 0 && 
                !ctx->local.len &&
                !ctx->flags &&
                *FIRSTENT (char, (ctx->ol_dict.len ? ctx->ol_dict : ctx->dict)) != SNIP_SELF_DELTA) // word_index=0 is the first word in the dictionary
                return true; 

            ctx->b250.len = base250_len (ctx->b250.data);
            ctx->flags |= CTX_FL_ALL_THE_SAME;  // set flag
        }
        else
            ctx->flags &= ~CTX_FL_ALL_THE_SAME; // reset flag
    }

    if (show) {
        bufprintf (vb, &vb->show_b250_buf, "%s", "\n")
        fprintf (stderr, "%.*s", (uint32_t)vb->show_b250_buf.len, vb->show_b250_buf.data);
        buf_free (&vb->show_b250_buf);
    }

    return false; // don't drop this section
}

// Codecs may be assigned in 3 stages:
// 1. During Seg (a must for all complex codecs - eg HT, DOMQ, ACGT...)
// 2. At merge - inherit from z_file->context if not set in Seg
// 3. After merge before compress - if still not assigned - zip_assign_best_codec - which also commits back to z_file->context
//    (this is the only place we commit to z_file, therefore z_file will only contain simple codecs)
// Note: if vb=1 commits a codec, it will be during its lock, so that all subsequent VBs will inherit it. But for
// contexts not committed by vb=1 - multiple contexts running in parallel may commit their codecs overriding each other. that's ok.
static void zip_assign_best_codec (VBlock *vb)
{
    if (flag.show_codec && vb->vblock_i == 1)
        fprintf (stderr, "\n\nThe output of --show-codec-test: Testing a sample of up %u bytes on ctx.local of each context.\n"
                 "Results in the format [codec size clock] are in order of quality - the first was selected.\n", CODEC_ASSIGN_SAMPLE_SIZE);

    for (int did_i=0 ; did_i < vb->num_contexts ; did_i++) {
        Context *ctx = &vb->contexts[did_i];
       
        // local
        if (ctx->local.len) {
            ctx->local.len *= lt_desc[ctx->ltype].width;
            codec_assign_best_codec (vb, ctx, NULL, SEC_LOCAL);
            ctx->local.len /= lt_desc[ctx->ltype].width;
        }

        // b250
        if (ctx->node_i.len * sizeof (uint32_t) < MIN_LEN_FOR_COMPRESSION)
            continue; // case: size is too small even before shrinking during generation

        // generate a sample of ctx.b250 data from ctx.node_i data
        zip_generate_b250_section (vb, ctx, CODEC_ASSIGN_SAMPLE_SIZE);

        codec_assign_best_codec (vb, ctx, NULL, SEC_B250);

        ctx->b250.len = 0; // roll back
    }
}

// after segging - if any context appears to contain only singleton snips (eg a unique ID),
// we move it to local instead of needlessly cluttering the global dictionary
static void zip_handle_unique_words_ctxs (VBlock *vb)
{
    for (int did_i=0 ; did_i < vb->num_contexts ; did_i++) {
        Context *ctx = &vb->contexts[did_i];
    
        if (!ctx->nodes.len || ctx->nodes.len != ctx->node_i.len) continue; // check that all words are unique (and new to this vb)
        if (vb->data_type == DT_VCF && dict_id_is_vcf_format_sf (ctx->dict_id)) continue; // this doesn't work for FORMAT fields
        if (ctx->nodes.len < vb->lines.len / 5)   continue; // don't bother if this is a rare field less than 20% of the lines
        if (buf_is_allocated (&ctx->local))     continue; // skip if we are already using local to optimize in some other way

        // don't move to local if its on the list of special dict_ids that are always in dict (because local is used for something else - eg pos or id data)
        if ((ctx->inst & CTX_INST_NO_STONS) || ctx->ltype != LT_TEXT) continue; // NO_STONS is implicit if ctx isn't text

        buf_move (vb, &ctx->local, vb, &ctx->dict);
        buf_free (&ctx->nodes);
        buf_free (&ctx->node_i);
    }
}

// generate & write b250 data for all primary fields of this data type
static void zip_generate_and_compress_ctxs (VBlock *vb)
{
    START_TIMER;

    // generate & write b250 data for all primary fields
    for (int did_i=0 ; did_i < vb->num_contexts ; did_i++) {
        Context *ctx = &vb->contexts[did_i];

        if (ctx->node_i.len) {
            
            bool drop_section = zip_generate_b250_section (vb, ctx, 0);

            if (!drop_section) {
                if (dict_id_printable (ctx->dict_id).num == flag.dump_one_b250_dict_id.num) 
                    ctx_dump_binary (vb, ctx, false);

                if (flag.show_time) codec_show_time (vb, "B250", ctx->name, ctx->bcodec);
                zfile_compress_b250_data (vb, ctx);
            }
        }

        if (ctx->local.len || (ctx->inst & CTX_INST_LOCAL_ALWAYS)) { 

            if (dict_id_printable (ctx->dict_id).num == flag.dump_one_local_dict_id.num) 
                ctx_dump_binary (vb, ctx, true);

            if (ctx->ltype == LT_BITMAP) 
                LTEN_bit_array (buf_get_bitarray (&ctx->local));

            if (flag.show_time) codec_show_time (vb, "LOCAL", ctx->name, ctx->lcodec);
            zfile_compress_local_data (vb, ctx, 0);
        }
    }

    COPY_TIMER (zip_generate_and_compress_ctxs);
}

static void zip_update_txt_counters (VBlock *vb)
{
    // note: in case of an FASTQ with flag.optimize_DESC, we already updated this in fastq_txtfile_count_lines
    if (!(flag.optimize_DESC && vb->data_type == DT_FASTQ)) 
        txt_file->num_lines += (int64_t)vb->lines.len; // lines in this txt file
        
    z_file->num_lines                        += (int64_t)vb->lines.len; // lines in all bound files in this z_file
    z_file->txt_data_so_far_single           += (int64_t)vb->vb_data_size;
    z_file->txt_data_so_far_bind             += (int64_t)vb->vb_data_size;
}

// write all the sections at the end of the file, after all VB stuff has been written
static void zip_write_global_area (Dispatcher dispatcher, Md5Hash single_component_md5)
{
    // if we're making a reference, we need the RA data to populate the reference section chrome/first/last_pos ahead of ref_compress_ref
    if (DTPZ(has_random_access)) 
        random_access_finalize_entries (&z_file->ra_buf); // sort RA, update entries that don't yet have a chrom_index

    ctx_compress_dictionaries();
    
    // store a mapping of the file's chroms to the reference's contigs, if they are any different
    if (flag.reference == REF_EXT_STORE || flag.reference == REF_EXTERNAL)
        ref_alt_chroms_compress();
    
    // output reference, if needed
    bool store_ref = flag.reference == REF_INTERNAL || flag.reference == REF_EXT_STORE || flag.make_reference;
    if (store_ref) 
        ref_compress_ref();
        
    if (flag.make_reference)
        refhash_compress_refhash();

    // add dict_id aliases list, if we have one
    Buffer *dict_id_aliases_buf = dict_id_create_aliases_buf();
    if (dict_id_aliases_buf->len) zfile_compress_section_data (evb, SEC_DICT_ID_ALIASES, dict_id_aliases_buf);

    // if this data has random access (i.e. it has chrom and pos), compress all random access records into evb->z_data
    if (DTPZ(has_random_access)) 
        random_access_compress (&z_file->ra_buf, SEC_RANDOM_ACCESS, flag.show_index ? "Random-access index contents (result of --show-index)" : NULL);

    if (store_ref) 
        random_access_compress (&ref_stored_ra, SEC_REF_RAND_ACC, flag.show_ref_index ? "Reference random-access index contents (result of --show-ref-index)" : NULL);

    // flush to disk before stats
    //zfile_output_processed_vb (evb);  

    stats_compress();

    // compress genozip header (including its payload sectionlist and footer) into evb->z_data
    zfile_compress_genozip_header (single_component_md5);    

    // output to disk stats and genozip header sections to disk
    //zfile_output_processed_vb (evb);  
}

// entry point for compute thread
static void zip_compress_one_vb (VBlock *vb)
{
    START_TIMER; 

    // if the txt file is compressed with BGZF, we uncompress now, in the compute thread
    if (txt_file->codec == CODEC_BGZF) 
        bgzf_uncompress_vb (vb);    // some of the blocks might already have been decompressed while reading - we decompress the remaining

    // calculate the MD5 contribution of this VB to the single file and bound files, and the MD5 snapshot of this VB
    if (flag.md5 && !flag.make_reference) txtfile_md5_one_vb (vb); 

    // allocate memory for the final compressed data of this vb. allocate 33% of the
    // vb size on the original file - this is normally enough. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, vb->vb_data_size / 3, 1.2, "z_data");

    // clone global dictionaries while granted exclusive access to the global dictionaries
    if (flag.pair != PAIR_READ_2) // in case of PAIR_READ_2, we already cloned in zip_one_file
        ctx_clone_ctx (vb);

    // split each line in this variant block to its components
    seg_all_data_lines (vb);

    // identify dictionaries that contain almost only unique words (eg a unique id) and move the data from dict to local
    zip_handle_unique_words_ctxs (vb);

    // for the first vb only - sort dictionaries so that the most frequent entries get single digit
    // base-250 indices. This can be done only before any dictionary is written to disk, but likely
    // beneficial to all vbs as they are likely to more-or-less have the same frequent entries
    if (vb->vblock_i == 1) {
        ctx_sort_dictionaries_vb_1(vb);

        // for a file that's compressed (and hence we don't know the size of content a-priori) AND we know the
        // compressed file size (eg a local gz/bz2 file or a compressed file on a ftp server) - we now estimate the 
        // txt_data_size_single that will be used for the global_hash and the progress indicator
        txtfile_estimate_txt_data_size (vb);
    }

    zfile_compress_vb_header (vb); // vblock header

    // merge new words added in this vb into the z_file.contexts, ahead of zip_generate_b250_section().
    // writing indices based on the merged dictionaries. dictionaries are compressed. 
    // all this is done while holding exclusive access to the z_file dictionaries.
    // note: vb>=2 will block here, until vb=1 is completed
    ctx_merge_in_vb_ctx(vb);

    // for each context with CODEC_UNKNOWN - i.e. that was not assigned in Seg AND not inherited from
    // a previous VB during merge - find the best codecs for its local and b250 data. 
    // These codecs will be committed to zf_ctx so that subsequent VBs inherit it during their merge
    zip_assign_best_codec (vb);

    if (vb->vblock_i == 1) ctx_vb_1_unlock(vb); 

    // merge in random access - IF it is used
    if (DTP(has_random_access)) 
        random_access_merge_in_vb (vb);

    // generate & compress b250 and local data for all ctxs (for reference files we don't output VBs)
    if (!flag.make_reference && !flag.test_seg)
        zip_generate_and_compress_ctxs (vb);

    // compress data-type specific sections
    DT_FUNC (vb, compress)(vb);

    // tell dispatcher this thread is done and can be joined.
    // thread safety: this isn't protected by a mutex as it will just be false until it at some point turns to true
    // this this operation needn't be atomic, but it likely is anyway
    vb->is_processed = true; 

    COPY_TIMER (compute);
}

// this is the main dispatcher function. It first processes the txt header, then proceeds to read 
// a variant block from the input file and send it off to a thread for computation. When the thread
// completes, this function proceeds to write the output to the output file. It can dispatch
// several threads in parallel.
void zip_one_file (const char *txt_basename, bool is_last_file)
{
    static DataType last_data_type = DT_NONE;
    static uint32_t prev_file_first_vb_i=0, prev_file_last_vb_i=0; // used if we're binding files - the vblock_i will continue from one file to the next
    
    if (!flag.bind) prev_file_first_vb_i = prev_file_last_vb_i = 0; // reset if we're not binding

    // we cannot bind files of different type
    ASSERT (!flag.bind || txt_file->data_type == last_data_type || last_data_type == DT_NONE, 
            "%s: cannot bind %s because it is a %s file, whereas the previous file was a %s",
             global_cmd, txt_name, dt_name (txt_file->data_type), dt_name (last_data_type));
    last_data_type =  txt_file->data_type;

    // normally global_max_threads would be the number of cores available - we allow up to this number of compute threads, 
    // because the I/O thread is normally idling waiting for the disk, so not consuming a lot of CPU
    Dispatcher dispatcher = dispatcher_init (global_max_threads, prev_file_last_vb_i, false, is_last_file, txt_basename, PROGRESS_PERCENT, 0);

    uint32_t first_vb_i = prev_file_last_vb_i + 1;

    dict_id_initialize (z_file->data_type);

    uint32_t txt_line_i = 1; // the next line to be read (first line = 1)
    
    // read the txt header, assign the global variables, and write the compressed header to the GENOZIP file
    off64_t txt_header_header_pos = z_file->disk_so_far;
    bool success = txtfile_header_to_genozip (&txt_line_i);
    if (!success) goto finish; // 2nd+ VCF file cannot bind, because of different sample names

    DT_FUNC (txt_file, zip_initialize)();

    ctx_initialize_for_zip();

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
            zfile_update_compressed_vb_header (processed_vb, txt_line_i); 

            max_lines_per_vb = MAX (max_lines_per_vb, processed_vb->lines.len);
            txt_line_i += (uint32_t)processed_vb->lines.len;

            if (!flag.make_reference && !flag.test_seg)
                zfile_output_processed_vb (processed_vb);
            
            zip_update_txt_counters (processed_vb);

            z_file->num_vbs++;
            
            dispatcher_finalize_one_vb (dispatcher);
        }        
        
        // PRIORITY 3: If there is no variant block available to compute or to output, but input is not exhausted yet - read one
        else if (!next_vb && !dispatcher_is_input_exhausted (dispatcher)) {

            next_vb = dispatcher_generate_next_vb (dispatcher, 0);

            // if we're compressing the 2nd file in a fastq pair (with --pair) - look back at the z_file data
            // and copy the data we need for this vb. note: we need to do this before txtfile_read_vblock as
            // we need the num_lines of the pair file
            bool read_txt = true;
            if (flag.pair == PAIR_READ_2) {

                // normally we clone in the compute thread, because it might wait on mutex, but in this
                // case we need to clone (i.e. create all contexts before we can read the pair file data)
                ctx_clone_ctx (next_vb); 

                // returns false if their is no vb with vb_i in the previous file
                read_txt = fastq_read_pair_1_data (next_vb, prev_file_first_vb_i, prev_file_last_vb_i);
            }

            if (read_txt) {
                if (flag.show_threads) dispatcher_show_time ("Read input data", -1, next_vb->vblock_i);            
                txtfile_read_vblock (next_vb);
                if (flag.show_threads) dispatcher_show_time ("Read input data done", -1, next_vb->vblock_i);
            }

            if (next_vb->txt_data.len)   // we found some data 
                next_vb->ready_to_dispatch = true;
            
            else {
                // error if stdin is empty - can happen only when redirecting eg "cat empty-file|./genozip -" (we test for empty regular files in main_genozip)
                ASSERT0 (next_vb->vblock_i > 1 || txt_file->txt_data_so_far_single /* txt header data */, "Error: Cannot compress stdin data because its size is 0");

                // this vb has no data
                dispatcher_set_input_exhausted (dispatcher, true);
                
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
    Md5Hash single_component_md5 = MD5HASH_NONE;
    if (z_file && !flag.test_seg && !z_file->redirected && txt_header_header_pos >= 0) 
        success = zfile_update_txt_header_section_header (txt_header_header_pos, max_lines_per_vb, &single_component_md5);

    // write the BGZF section containing BGZF block sizes, if this txt file is compressed with BGZF
    bgzf_compress_bgzf_section();

    // if this a non-bound file, or the last component of a bound file - write the genozip header, random access and dictionaries
finish:
    if ((is_last_file || !flag.bind) && !flag.test_seg)
        zip_write_global_area (dispatcher, single_component_md5);

    zip_display_compression_ratio (dispatcher, flag.bind ? MD5HASH_NONE : single_component_md5, is_last_file || !flag.bind); // Done for reference + final compression ratio calculation
    
    if (flag.md5 && flag.bind && z_file->num_txt_components_so_far > 1 && is_last_file) 
        progress_concatenated_md5 (dt_name (z_file->data_type), md5_finalize (&z_file->md5_ctx_bound));

    z_file->disk_size              = z_file->disk_so_far;
    z_file->txt_data_so_far_single = 0;
    evb->z_data.len                = 0;
    evb->z_next_header_i           = 0;
    memset (&z_file->md5_ctx_single, 0, sizeof (Md5Context));

    prev_file_first_vb_i = first_vb_i;
    dispatcher_finish (&dispatcher, &prev_file_last_vb_i);

    DT_FUNC (txt_file, zip_finalize)();
}
