// ------------------------------------------------------------------
//   zip.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <math.h>
#include <errno.h>
#include <sys/types.h>
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
#include "random_access.h"
#include "dict_id.h"
#include "reference.h"
#include "refhash.h"
#include "ref_iupacs.h"
#include "progress.h"
#include "mutex.h"
#include "fastq.h"
#include "stats.h"
#include "codec.h"
#include "compressor.h"
#include "strings.h"
#include "bgzf.h"
#include "linesorter.h"
#include "txtheader.h"
#include "threads.h"
#include "endianness.h"
#include "segconf.h"
#include "contigs.h"
#include "chrom.h"
#include "biopsy.h"

static Mutex wait_for_vb_1_mutex = {};

static void zip_display_compression_ratio (Digest md5, bool is_final_component)
{
    float z_bytes        = MAX_((float)z_file->disk_so_far, 1.0); // at least one, to avoid division by zero in case of a z_bytes=0 issue
    float plain_bytes    = (float)z_file->txt_data_so_far_bind;
    float comp_bytes     = file_is_read_via_ext_decompressor (txt_file) 
                              ? (float)txt_file->disk_size    // 0 if via pipe or url, as we have no knowledge of file size
                              : (float)txt_file->disk_so_far; // unlike disk_size, works also for piped-in files (but not CRAM, BCF, XZ, ZIP)
    float ratio_vs_plain = plain_bytes / z_bytes;
    float ratio_vs_comp  = -1;

    if (flag.debug_progress) 
        iprintf ("Ratio calculation: ratio_vs_plain=%f = plain_bytes=%"PRIu64" / z_bytes=%"PRIu64"\n",
                    ratio_vs_plain, (uint64_t)plain_bytes, (uint64_t)z_bytes);

    // in bind mode, we don't show compression ratio for files except for the last one
    if (flag.bind) { 

        static float comp_bytes_bind = 0;
        static FileType source_file_type = UNKNOWN_FILE_TYPE;

        // reset for every set of bound files (we might have multiple sets if --pair)
        if (z_file->num_txt_components_so_far == 1) {
            comp_bytes_bind=0; 
            source_file_type = txt_file->type;
        }

        else if (source_file_type != txt_file->type) // heterogenous source file types
            source_file_type = UNKNOWN_FILE_TYPE;

        comp_bytes_bind += comp_bytes;

        if (is_final_component) { 
            ratio_vs_comp = comp_bytes_bind / z_bytes; // compression vs .gz/.bz2/.bcf/.xz... size
            if (flag.debug_progress) 
                iprintf ("Ratio calculation: ratio_vs_comp=%f = comp_bytes_bind=%"PRIu64" / z_bytes=%"PRIu64"\n",
                         ratio_vs_comp, (uint64_t)comp_bytes_bind, (uint64_t)z_bytes);
        }
        else 
            progress_finalize_component_time ("Done", md5);
    }
    else {
        ratio_vs_comp = comp_bytes / z_bytes; // compression vs .gz/.bz2/.bcf/.xz... size
        if (flag.debug_progress) 
            iprintf ("Ratio calculation: ratio_vs_comp=%f = comp_bytes=%"PRIu64" / z_bytes=%"PRIu64"\n",
                     ratio_vs_comp, (uint64_t)comp_bytes, (uint64_t)z_bytes);
    }

    // when making a reference, we don't care about the compression
    if (flag.make_reference)
        progress_finalize_component_time ("Done", md5);

    // when compressing BAM report only ratio_vs_comp (compare to BGZF-compress BAM - we don't care about the underlying plain BAM)
    else if (Z_DT(DT_BAM)) 
            progress_finalize_component_time_ratio (dt_name (z_file->data_type), ratio_vs_comp, md5);

    else if (ratio_vs_comp >= 0) {
        if (txt_file->codec == CODEC_NONE || ratio_vs_comp < 1.05) // disk_so_far doesn't give us the true txt file size 
            progress_finalize_component_time_ratio (dt_name (z_file->data_type), ratio_vs_plain, md5);
        
        else // source was compressed
            progress_finalize_component_time_ratio_better (dt_name (z_file->data_type), ratio_vs_plain, file_exts[txt_file->type], ratio_vs_comp, md5);
    }
}

// B250 generation means re-writing the b250 buffer with 2 modifications:
// 1) We write the word_index (i.e. the sequential number of the word in the dict in z_file) instead
// of the VB's node_index (that is private to this VB)
// 2) We optimize the representation of word_index giving privilage to 3 popular words (the most popular
//    words in VB=1) to be represented by a single byte
static inline void zip_generate_one_b250 (VBlockP vb, ContextP ctx, uint32_t word_i, WordIndex node_index,
                                          Buffer *b250_buf, 
                                          WordIndex *prev_word_index,  // in/out
                                          bool show)
{
    if (node_index >= 0) { // normal index

        Base250 b250 = ctx_node_vb (ctx, node_index, NULL, NULL)->word_index;

        WordIndex n           = b250.n;
        unsigned num_numerals = base250_len (b250.encoded.numerals);
        uint8_t *numerals     = b250.encoded.numerals;
        
        bool one_up = (n == *prev_word_index + 1) && (word_i > 0);

        if (one_up) { // note: we can't do SEC_VCF_GT_DATA bc we can't PIZ it as many GT data types are in the same section 
            NEXTENT(uint8_t, *b250_buf) = (uint8_t)BASE250_ONE_UP;
            if (show) bufprintf (vb, &vb->show_b250_buf, "L%u:ONE_UP ", word_i);
        }

        else {
            if (num_numerals == 4) { // assign word byte by byte, as it is not word-boundary aligned
                memcpy (AFTERENT (char, *b250_buf), numerals, 4); // hopefully the compiler will optimize this memcpy and it won't e a function call...
                b250_buf->len += 4;
            } 
            else
                NEXTENT (uint8_t, *b250_buf) = *numerals;

            if (show) bufprintf (vb, &vb->show_b250_buf, "L%u:%u ", word_i, n);
        }
        *prev_word_index = n;
    }

    else if (node_index == WORD_INDEX_MISSING) {
        NEXTENT(uint8_t, *b250_buf) = (uint8_t)BASE250_MISSING_SF;
        *prev_word_index = node_index;
        if (show) bufprintf (vb, &vb->show_b250_buf, "L%u:MISSING ", word_i);
    }
    
    else if (node_index == WORD_INDEX_EMPTY) {
        NEXTENT(uint8_t, *b250_buf) = (uint8_t)BASE250_EMPTY_SF;
        *prev_word_index = node_index;
        if (show) bufprintf (vb, &vb->show_b250_buf, "L%u:EMPTY ", word_i);
    }

    else ABORT ("Error in zip_generate_one_b250: invalid node_index=%u", node_index);        
}

// here we generate the b250 data: we convert the b250 buffer data from an index into context->nodes
// to an index into the to-be-generated-in-piz context->word_index, encoded in base250.
// Note that the word indices have changed since segmentation (which is why we needed this intermediate step)
// because: 1. the dictionary got integrated into the global one - some values might have already been in the global
// dictionary thanks to other threads working on other VBs  
// 2. for the first VB, we sort the dictionary by frequency returns true if section should be dropped
static void zip_generate_b250 (VBlock *vb, Context *ctx)
{
    START_TIMER;

    ASSERT (ctx->dict_id.num, "tag_name=%s did_i=%u: ctx->dict_id=0 despite ctx->b250 containing data", ctx->tag_name, (unsigned)(ctx - vb->contexts));

    bool show = flag.show_b250 || dict_id_typeless (ctx->dict_id).num == flag.dict_id_show_one_b250.num;
    if (show) bufprintf (vb, &vb->show_b250_buf, "vb_i=%u %s: ", vb->vblock_i, ctx->tag_name);

    // case: all-the-same b250 survived dropping - we just shorten it to one entry
    if (ctx->flags.all_the_same && ctx->b250.len > 1) { 
        if (flag.debug_generate) iprintf ("%s.b250 vb_i=%u is \"all_the_same\" - shortened b250 from len=%"PRIu64" to 1\n", ctx->tag_name, vb->vblock_i, ctx->b250.len);
        ctx->b250.len = 1;
    }
    else
        if (flag.debug_generate) iprintf ("%s.b250 vb_i=%u len=%"PRIu64" all_the_same=%u\n", ctx->tag_name, vb->vblock_i, ctx->b250.len, ctx->flags.all_the_same);

    ARRAY (WordIndex, b250, ctx->b250);

    ctx->b250.num_ctx_words = b250_len; // we move the number of words to param, as len will now contain the of bytes. used by ctx_update_stats()
    ctx->b250.len = 0; // we are going to overwrite b250 with the converted indices

    // convert b250 from an array of node_index (4 bytes each) to word_index (1-4 bytes each)    
    WordIndex prev = WORD_INDEX_NONE; 
    for (uint64_t word_i=0; word_i < b250_len; word_i++) 
        zip_generate_one_b250 (vb, ctx, word_i, b250[word_i], &ctx->b250, &prev, show);

    if (show) {
        bufprintf (vb, &vb->show_b250_buf, "%s", "\n");
        iprintf ("%.*s", (uint32_t)vb->show_b250_buf.len, vb->show_b250_buf.data);
        buf_free (&vb->show_b250_buf);
    }

    codec_assign_best_codec (vb, ctx, NULL, SEC_B250);

    COPY_TIMER (zip_generate_b250);
}

static void zip_resize_local (VBlock *vb, Context *ctx)
{
    ARRAY (int64_t, src, ctx->local);

    // search for the largest (stop if largest so far requires 32bit)
    int64_t largest = *FIRSTENT (int64_t, ctx->local);
    int64_t smallest = largest;

    for (uint64_t i=1; i < src_len; i++) 
        if (src[i] > largest) largest = src[i];
        else if (src[i] < smallest) smallest = src[i];

    // uint8
    if (smallest >= 0 && largest <= 0xffLL) {
        ctx->ltype = LT_UINT8;
        ARRAY (uint8_t, dst, ctx->local);
        for (uint64_t i=0; i < src_len; i++)
            dst[i] = (uint8_t)src[i];
    }

    // uint16
    else if (smallest >= 0 && largest <= 0xffffLL) {
        ctx->ltype = LT_UINT16;
        ARRAY (uint16_t, dst, ctx->local);
        for (uint64_t i=0; i < src_len; i++)
            dst[i] = BGEN16 ((uint16_t)src[i]);
    }

    // uint32
    else if (smallest >= 0 && largest <= 0xffffffffLL) {
        ctx->ltype = LT_UINT32;
        ARRAY (uint32_t, dst, ctx->local);
        for (uint64_t i=0; i < src_len; i++)
            dst[i] = BGEN32 ((uint32_t)src[i]);
    }

    // int8
    else if (smallest >= -0x80LL && largest <= 0x7fLL) {
        ctx->ltype = LT_INT8;
        ARRAY (uint8_t, dst, ctx->local);
        for (uint64_t i=0; i < src_len; i++)
            dst[i] = INTERLACE(int8_t, src[i]);
    }

    // int16
    else if (smallest >= -0x8000LL && largest <= 0x7fffLL) {
        ctx->ltype = LT_INT16;
        ARRAY (uint16_t, dst, ctx->local);
        for (uint64_t i=0; i < src_len; i++)
            dst[i] = BGEN16 (INTERLACE(int16_t, src[i]));
    }

    // int32
    else if (smallest >= -0x80000000LL && largest <= 0x7fffffffLL) {
        ctx->ltype = LT_INT32;
        ARRAY (uint32_t, dst, ctx->local);
        for (uint64_t i=0; i < src_len; i++)
            dst[i] = BGEN32 (INTERLACE(int32_t,src[i]));
    }

    // uint64
    else if (smallest >= 0) {
        ctx->ltype = LT_UINT64;
        for (uint64_t i=0; i < src_len; i++)
            src[i] = BGEN64 (src[i]);
    }

    // int64
    else {
        ctx->ltype = LT_INT64;
        for (uint64_t i=0; i < src_len; i++)
            src[i] = BGEN64 (INTERLACE(int64_t, src[i]));
    }
}

// selects the smallest size (8, 16, 32) for the data, transposes, and BGENs
static void zip_generate_transposed_local (VBlock *vb, Context *ctx)
{
    ARRAY (uint32_t, data, ctx->local);

    // get largest element (excluding 0xffffffff - representing a '.')
    uint32_t largest=0;
    for (uint64_t i=0; i < ctx->local.len; i++)
        if (data[i] > largest && data[i] != 0xffffffff) largest = data[i];

    if      (largest < 0xfe)   ctx->ltype = LT_UINT8_TR;  // -1 is reserved for "missing"
    else if (largest < 0xfffe) ctx->ltype = LT_UINT16_TR;
    
    buf_alloc (vb, &vb->compressed, 0, ctx->local.len * lt_desc[ctx->ltype].width, char, 1, "compressed");

    uint32_t cols = ctx->local.param;
    // we're restricted to 255 columns, because this number goes into uint8_t SectionHeaderCtx.param
    ASSERT (cols >= 0 && cols <= 255, "columns=%u out of range [1,255] in transposed matrix %s", cols, ctx->tag_name);

    if (!cols) cols = vcf_header_get_num_samples(); // 0 if not vcf/bcf (not restricted to 255)
    
    // case: matrix is not transposable - just BGEN it
    if (ctx->local.len % cols) {
        ctx->ltype = LT_UINT32; // not transposed

        for (unsigned i=0; i < ctx->local.len; i++)
            data[i] = BGEN32 (data[i]);
            
        goto done;
    }

    uint32_t rows = ctx->local.len / cols;

    ctx->flags.copy_local_param = true;
    
    for (uint32_t r=0; r < rows; r++) 
        for (uint32_t c=0; c < cols; c++) {

            uint32_t value = data[r * cols + c];

            switch (ctx->ltype) { // note: the casting aslo correctly converts 0xffffffff to eg 0xff
                case LT_UINT8_TR  : *ENT (uint8_t,  vb->compressed, c * rows + r) =         (uint8_t)value;   break;
                case LT_UINT16_TR : *ENT (uint16_t, vb->compressed, c * rows + r) = BGEN16 ((uint16_t)value); break;
                case LT_UINT32_TR : *ENT (uint32_t, vb->compressed, c * rows + r) = BGEN32 (value);           break;
                default: ABORT ("Error in zip_generate_transposed_local: Bad ltype=%s", lt_name (ctx->ltype));
            }
        }

    vb->compressed.len = ctx->local.len;
    buf_copy_do (vb, &ctx->local, &vb->compressed, lt_desc[ctx->ltype].width, 0, 0, __FUNCTION__,__LINE__, "contexts->local"); // copy and not move, so we can keep local's memory for next vb

done:
    buf_free (&vb->compressed);
}

// after segging - if any context appears to contain only singleton snips (eg a unique ID),
// we move it to local instead of needlessly cluttering the global dictionary
static void zip_handle_unique_words_ctxs (VBlock *vb)
{
    for (ContextP ctx=CTX(0); ctx < CTX(vb->num_contexts); ctx++) {
    
        if (!ctx->nodes.len || ctx->nodes.len != ctx->b250.len) continue; // check that all words are unique (and new to this vb)
        if (VB_DT(DT_VCF) && dict_id_is_vcf_format_sf (ctx->dict_id)) continue; // this doesn't work for FORMAT fields
        if (ctx->nodes.len < vb->lines.len / 5) continue; // don't bother if this is a rare field less than 20% of the lines
        if (buf_is_alloc (&ctx->local))     continue; // skip if we are already using local to optimize in some other way

        // don't move to local if its on the list of special dict_ids that are always in dict (because local is used for something else - eg pos or id data)
        if (!ctx_can_have_singletons (ctx) ||
            ctx->b250.len == 1) continue; // handle with all_the_same rather than singleton
         
        buf_move (vb, &ctx->local, vb, &ctx->dict);
        buf_free (&ctx->nodes);
        buf_free (&ctx->b250);
        ctx->flags.all_the_same = false;
    }
}

static void zip_generate_local (VBlockP vb, ContextP ctx)
{
    START_TIMER;

    if (!ctx->local.len) return;

    ASSERT (ctx->dict_id.num, "tag_name=%s did_i=%u: ctx->dict_id=0 despite ctx->local containing data", ctx->tag_name, (unsigned)(ctx - vb->contexts));

    if (ctx->dynamic_size_local) 
        zip_resize_local (vb, ctx);

    else if (ctx->ltype == LT_BITMAP) 
        LTEN_bit_array (buf_get_bitarray (&ctx->local));

    else if (ctx->ltype == LT_UINT32_TR)
        zip_generate_transposed_local (vb, ctx);

    else if (ctx->ltype == LT_FLOAT32)
        BGEN_u32_buf (&ctx->local, NULL);
        
    else if (ctx->ltype == LT_FLOAT64)
        BGEN_u64_buf (&ctx->local, NULL);
        
    codec_assign_best_codec (vb, ctx, NULL, SEC_LOCAL);

    if (flag.debug_generate) iprintf ("%s.local in vb_i=%u ltype=%s len=%"PRIu64" codec=%s\n", 
                                      ctx->tag_name, vb->vblock_i, lt_name (ctx->ltype), ctx->local.len, codec_name(ctx->lcodec));

    COPY_TIMER (zip_generate_local);
}

// generate & write b250 data for all contexts - do them in random order, to reduce the chance of multiple doing codec_assign_best_codec for the same context at the same time
// VBs doing codec_assign_best_codec at the same, so that they can benefit from pre-assiged codecs
void zip_compress_all_contexts_b250 (VBlock *vb)
{
    START_TIMER;
    threads_log_by_vb (vb, "zip", "START COMPRESSING B250", 0);
    
    // arrays of all contexts in this VB
    ContextP ctxs[vb->num_contexts];
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) ctxs[did_i] = CTX(did_i);

        // in each iteration, pick a context at random and remove it from the list 
    for (unsigned i=0; i < vb->num_contexts; i++) {
 
        int ctx_i = global_max_threads > 1 ? ((clock()+1) * (vb->vblock_i+1)) % (vb->num_contexts - i) : 0; // force predictability with single thread 
        
        ContextP ctx = ctxs[ctx_i];
        memmove (&ctxs[ctx_i], &ctxs[ctx_i+1], (vb->num_contexts - i - ctx_i - 1) * sizeof (ContextP));

        if (!ctx->b250.len || ctx->b250_compressed) continue;

        zip_generate_b250 (vb, ctx); // generate the final b250 buffers from their intermediate form

        if (dict_id_typeless (ctx->dict_id).num == flag.dump_one_b250_dict_id.num) 
            ctx_dump_binary (vb, ctx, false);

        if (flag.show_time) codec_show_time (vb, "B250", ctx->tag_name, ctx->bcodec);
        
        if (HAS_DEBUG_SEG(ctx)) iprintf ("zip_compress_all_contexts_b250: vb_i=%u %s: B250.len=%"PRIu64" NODES.len=%"PRIu64"\n", 
                                     vb->vblock_i, ctx->tag_name, ctx->b250.len, ctx->nodes.len);

        START_TIMER; // for compressor_time

        zfile_compress_b250_data (vb, ctx);

        if (flag.show_time) 
            ctx->compressor_time += CHECK_TIMER; // sum b250 and local

        ctx->b250_compressed = true;
    }

    COPY_TIMER (zip_compress_ctxs); // same profiler for b250 and local as we breakdown by ctx underneath it
}

// generate & write local data for all contexts - in random order, to reduce the chance of multiple doing codec_assign_best_codec for the same context at the same time
static void zip_compress_local (VBlock *vb)
{
    START_TIMER;
    threads_log_by_vb (vb, "zip", "START COMPRESSING LOCAL", 0);

    // first we handle local_dep=0 then local_dep=1 and finally local_dep=2
    for (int dep_level=DEP_L0 ; dep_level < NUM_LOCAL_DEPENDENCY_LEVELS; dep_level++) {

        // initialize list of contexts at this dependency level that need compression
        ContextP ctxs[vb->num_contexts];
        unsigned num_ctxs=0;
        for (ContextP ctx=CTX(0); ctx < CTX(vb->num_contexts); ctx++) 
            if ((ctx->local.len || ctx->local_always) && ctx->local_dep == dep_level && !ctx->local_compressed)
                ctxs[num_ctxs++] = ctx;

        while (num_ctxs) {
            // pick a context at "random" and remove it from the list (not random if single thread)
            int ctx_i = global_max_threads > 1 ? (65531 * (vb->vblock_i+1)) % num_ctxs : 0; 
            ContextP ctx = ctxs[ctx_i];
            memmove (&ctxs[ctx_i], &ctxs[ctx_i+1], (num_ctxs - (ctx_i+1)) * sizeof (ContextP));
            num_ctxs--;

            zip_generate_local (vb, ctx);

            if (dict_id_typeless (ctx->dict_id).num == flag.dump_one_local_dict_id.num) 
                ctx_dump_binary (vb, ctx, true);

            if (flag.show_time) codec_show_time (vb, "LOCAL", ctx->tag_name, ctx->lcodec);

            if (HAS_DEBUG_SEG(ctx)) iprintf ("zip_compress_local: vb_i=%u %s: LOCAL.len=%"PRIu64" LOCAL.param=%"PRIu64"\n", 
                                            vb->vblock_i, ctx->tag_name, ctx->local.len, ctx->local.param);

            START_TIMER; // for compressor_time

            zfile_compress_local_data (vb, ctx, 0);

            ctx->local_compressed = true; // so we don't compress it again
            ctx->no_stons = true; // since we had data on local, we don't allow ctx_commit_node to move singletons to local
                
            if (flag.show_time) 
                ctx->compressor_time += CHECK_TIMER; // sum b250 and local
        }
    }

    COPY_TIMER (zip_compress_ctxs); // same profiler for b250 and local as we breakdown by ctx underneath it
}

// called by main thread after VB has completed processing
static void zip_update_txt_counters (VBlock *vb)
{
    // note: in case of an FASTQ with flag.optimize_DESC or VCF with add_line_numbers, we already updated this in *_zip_read_one_vb
    if (!(flag.optimize_DESC && VB_DT(DT_FASTQ)) &&
        !(flag.add_line_numbers && (VB_DT(DT_VCF) || VB_DT(DT_BCF))))
        txt_file->num_lines += (int64_t)vb->lines.len; // lines in this txt file

    // counters of data AS IT APPEARS IN THE TXT FILE
    if (!flag.rejects_coord) {      
        z_file->num_lines                += (int64_t)vb->lines.len; // lines in all bound files in this z_file
        z_file->txt_data_so_far_single_0 += (int64_t)vb->txt_size;  // length of data before any modifications
        z_file->txt_data_so_far_bind_0   += (int64_t)vb->txt_size;
    }

    // counter of data FOR PROGRESS DISPLAY
    z_file->txt_data_so_far_single += (int64_t)vb->txt_size;   

    // counter of data in DEFAULT PRIMARY RECONSTRUCTION
    if (flag.rejects_coord == DC_NONE) z_file->txt_data_so_far_bind += vb->recon_size;
    if (flag.rejects_coord == DC_LUFT) z_file->txt_data_so_far_bind += vb->recon_size_luft; // luft reject VB shows in primary reconstruction

    // add up context compress time
    if (flag.show_time)
        ctx_add_compressor_time_to_zf_ctx (vb);
}

// write all the sections at the end of the file, after all VB stuff has been written
static void zip_write_global_area (Digest single_component_digest)
{
    // if we're making a reference, we need the RA data to populate the reference section chrome/first/last_pos ahead of ref_compress_ref
    random_access_finalize_entries (&z_file->ra_buf); // sort RA, update entries that don't yet have a chrom_index
    random_access_finalize_entries (&z_file->ra_buf_luft); 

    // if we used the aligner with REF_EXT_STORE, we make sure all the CHROMs referenced are in the CHROM context, so
    // as SEC_REF_CONTIGS refers to them. We do this by seeing which contigs have any bit set in is_set.
    // note: in REF_EXTERNAL we don't use is_set, so we populate all contigs in zip_initialize
    if (flag.aligner_available && flag.reference == REF_EXT_STORE)
        ref_contigs_populate_aligned_chroms();

    ctx_compress_dictionaries(); 

    ctx_compress_counts();
        
    // store a mapping of the file's chroms to the reference's contigs, if they are any different
    // note: not needed in REF_EXT_STORE, as we convert the stored ref_contigs to use chrom_index of the file's CHROM
    if (flag.reference == REF_EXTERNAL) 
        chrom_2ref_compress(gref);

    // output reference, if needed
    bool store_ref = (flag.reference & REF_STORED) || flag.make_reference;
    if (store_ref) 
        ref_compress_ref();
        
    if (flag.make_reference) {
        ref_iupacs_compress();
        refhash_compress_refhash();
    }

    // add dict_id aliases list, if we have one
    Buffer *dict_id_aliases_buf = dict_id_create_aliases_buf();
    if (dict_id_aliases_buf->len) zfile_compress_section_data (evb, SEC_DICT_ID_ALIASES, dict_id_aliases_buf);

    // if this data has random access (i.e. it has chrom and pos), compress all random access records into evb->z_data
    random_access_compress (&z_file->ra_buf,      SEC_RANDOM_ACCESS, DC_PRIMARY, flag.show_index ? RA_MSG_PRIM : NULL);
    random_access_compress (&z_file->ra_buf_luft, SEC_RANDOM_ACCESS, DC_LUFT,    flag.show_index ? RA_MSG_LUFT : NULL);

    if (store_ref) 
        random_access_compress (ref_get_stored_ra (gref), SEC_REF_RAND_ACC, 0, flag.show_ref_index ? RA_MSG_REF : NULL);

    stats_compress();

    // compress genozip header (including its payload sectionlist and footer) into evb->z_data
    zfile_compress_genozip_header (single_component_digest);    
}

// entry point for ZIP compute thread
static void zip_compress_one_vb (VBlock *vb)
{
    START_TIMER; 

    if (flag.biopsy) return; // we're just taking a biopsy of the txt data, so no need to actually compress

    // if the txt file is compressed with BGZF, we uncompress now, in the compute thread
    if (txt_file->codec == CODEC_BGZF && flag.pair != PAIR_READ_2) 
        bgzf_uncompress_vb (vb);    // some of the blocks might already have been decompressed while reading - we decompress the remaining

    // calculate the digest contribution of this VB to the single file and bound files, and the digest snapshot of this VB
    if (!flag.make_reference && !flag.data_modified && !flag.rejects_coord) 
        digest_one_vb (vb); 

    // allocate memory for the final compressed data of this vb. allocate 33% of the
    // vb size on the original file - this is normally enough. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, 0, vb->txt_size / 3, char, CTX_GROWTH, "z_data");

    // clone global dictionaries while granted exclusive access to the global dictionaries
    if (flag.pair != PAIR_READ_2) // in case of PAIR_READ_2, we already cloned in zip_one_file
        ctx_clone (vb);

    // split each line in this variant block to its components
    threads_log_by_vb (vb, "zip", "START SEG", 0);
    seg_all_data_lines (vb);

    // identify dictionaries that contain only singleton words (eg a unique id) and move the data from dict to local
    zip_handle_unique_words_ctxs (vb);

    // for the first vb only - sort dictionaries so that the most frequent entries get single digit
    // base-250 indices. This can be done only before any dictionary is written to disk, but likely
    // beneficial to all vbs as they are likely to more-or-less have the same frequent entries
    if (vb->vblock_i == 1) 
        ctx_sort_dictionaries_vb_1(vb);

    zfile_compress_vb_header (vb); // vblock header

    if (flag.show_codec) {
        DO_ONCE iprintf ("\n\nThe output of --show-codec-test: Testing a sample of up %u bytes on ctx.local of each context.\n"
                         "Results in the format [codec size clock] are in order of quality - the first was selected.\n", CODEC_ASSIGN_SAMPLE_SIZE);
    }

    if (vb->vblock_i != 1) {
        if (!flag.make_reference && !flag.seg_only) 
            // while vb_i=1 is busy merging, other VBs can handle local
            zip_compress_local (vb); // not yet locals that consist of singletons transferred from dict to local in ctx_merge_in_vb_ctx (these will have len=0 at this point)
        
        // let vb_i=1 merge first, as it has the sorted dictionaries, other vbs can go in arbitrary order. 
        START_TIMER;
        mutex_wait (wait_for_vb_1_mutex);
        COPY_TIMER(wait_for_vb_1_mutex);
    }

    // merge new words added in this vb into the z_file.contexts, ahead of zip_generate_b250().
    // writing indices based on the merged dictionaries. dictionaries are compressed. 
    // all this is done while holding exclusive access to the z_file dictionaries.
    // note: vb>=2 will block here, until vb=1 is completed
    threads_log_by_vb (vb, "zip", "START MERGE", 0);
    ctx_merge_in_vb_ctx(vb);

    if (vb->vblock_i == 1) 
        mutex_unlock (wait_for_vb_1_mutex); 

    if (!flag.make_reference && !flag.seg_only) {
        zip_compress_local (vb); // for vb=1 - all locals ; for vb>1 - locals which consist of singletons set in ctx_merge_in_vb_ctx (other locals were already compressed above)
        zip_compress_all_contexts_b250 (vb);
    }

    // case: we're sorting the file's lines - add line data to txt_file's line_info (VB order doesn't matter - we will sort them later)
    if (flag.sort) 
        linesorter_merge_vb (vb);

    // merge in random access - IF it is used
    random_access_merge_in_vb (vb, DC_PRIMARY);
    random_access_merge_in_vb (vb, DC_LUFT);
    
    // compress data-type specific sections
    DT_FUNC (vb, compress)(vb);

    // tell dispatcher this thread is done and can be joined.
    // thread safety: this isn't protected by a mutex as it will just be false until it at some point turns to true
    vb->is_processed = true; 

    COPY_TIMER (compute);
}

// data sent through dispatcher fan out functions - to do: make this an opaque struct
static uint32_t prev_file_first_vb_i=0, prev_file_last_vb_i=0; // used if we're binding files - the vblock_i will continue from one file to the next
static uint32_t max_lines_per_vb; // (resets in every file)

// returns true if successfully prepared a vb 
static void zip_prepare_one_vb_for_dispatching (VBlockP vb)
{
    // if we're compressing the 2nd file in a fastq pair (with --pair) - look back at the z_file data
    // and copy the data we need for this vb. note: we need to do this before txtfile_read_vblock as
    // we need the num_lines of the pair file
    bool read_txt = true;
    if (flag.pair == PAIR_READ_2) {

        // normally we clone in the compute thread, because it might wait on mutex, but in this
        // case we need to clone (i.e. create all contexts before we can read the pair file data)
        ctx_clone (vb); 
    
        uint32_t pair_vb_i = prev_file_first_vb_i + (vb->vblock_i - prev_file_last_vb_i - 1);
        
        read_txt = (pair_vb_i <= prev_file_last_vb_i)  // false if there is no vb with vb_i in the previous file
                 && fastq_read_pair_1_data (vb, pair_vb_i, false); // read here, decompressed in fastq_seg_initialize
    }

    if (read_txt) {
        txtfile_read_vblock (vb);

        // initializations after reading the first vb and before running any compute thread
        if (vb->vblock_i == 1) { 
            // we lock, so vb>=2 will block on merge, until vb=1 has completed its merge and unlocked it in zip_compress_one_vb()
            mutex_destroy (wait_for_vb_1_mutex); // destroy mutex of previous file (it will still be locked if that file had no VBs)
            mutex_initialize (wait_for_vb_1_mutex);
            mutex_lock (wait_for_vb_1_mutex); 
        }

        vb->vb_coords = !z_dual_coords ? DC_PRIMARY
                      : flag.rejects_coord == DC_NONE ? DC_BOTH
                      : flag.rejects_coord;

        vb->is_rejects_vb = (flag.rejects_coord != DC_NONE);
    }

    if (vb->txt_data.len)   // we found some data 
        vb->ready_to_dispatch = true;

    else {
        // error if stdin is empty - can happen only when redirecting eg "cat empty-file|./genozip -" (we test for empty regular files in main_genozip)
        ASSINP0 (vb->vblock_i > 1 || txt_file->txt_data_so_far_single /* txt header data */, 
                 "Error: Cannot compress stdin data because its size is 0");
    }
}

// called main thread, in order of VBs
static void zip_complete_processing_one_vb (VBlockP vb)
{
    if (DTP(zip_after_compute)) DTP(zip_after_compute)(vb);

    // update z_data in memory (its not written to disk yet)
    zfile_update_compressed_vb_header (vb); 

    max_lines_per_vb = MAX_(max_lines_per_vb, vb->lines.len);

    if (!flag.make_reference && !flag.seg_only)
        zfile_output_processed_vb (vb);
    
    zip_update_txt_counters (vb);

    z_file->num_vbs++;
    txt_file->num_vbs++;
}

// this is the main dispatcher function. It first processes the txt header, then proceeds to read 
// a variant block from the input file and send it off to a thread for computation. When the thread
// completes, this function proceeds to write the output to the output file. It can dispatch
// several threads in parallel.
void zip_one_file (const char *txt_basename, 
                   bool *is_last_file,     // in/out very last file in this execution
                   bool z_closes_after_me) // we will finalize this z_file after writing this component
{
    static DataType last_data_type = DT_NONE;
    Dispatcher dispatcher = 0;
    dispatcher_start_wallclock();

    z_file->txt_data_so_far_single = 0;
    evb->z_data.len                = 0;
    evb->z_next_header_i           = 0;
    memset (&z_file->digest_ctx_single, 0, sizeof (z_file->digest_ctx_single));

    if (!flag.bind) prev_file_first_vb_i = prev_file_last_vb_i = 0; // reset if we're not binding

    // we cannot bind files of different type
    ASSINP (!flag.bind || txt_file->data_type == last_data_type || last_data_type == DT_NONE, 
            "%s: cannot bind %s because it is a %s file, whereas the previous file was a %s",
            global_cmd, txt_name, dt_name (txt_file->data_type), dt_name (last_data_type));
    last_data_type =  txt_file->data_type;

    uint32_t first_vb_i = prev_file_last_vb_i + 1;

    if (!z_file->num_txt_components_so_far)  // first component of this z_file 
        ctx_initialize_predefined_ctxs (z_file->contexts, txt_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_contexts);

    segconf_initialize(); // before txtheader 

    // read the txt header, assign the global variables, and write the compressed header to the GENOZIP file
    uint64_t txt_header_header_pos = z_file->disk_so_far;
    uint64_t txt_header_size;
    bool success = txtheader_zip_read_and_compress (&txt_header_size); // also increments z_file->num_txt_components_so_far
    if (!success) goto finish; // eg 2nd+ VCF file cannot bind, because of different sample names

    if (z_dual_coords && !flag.rejects_coord) { 
        *is_last_file = false; // we're not the last file - at leasts, the liftover rejects file is after us
        z_closes_after_me = false;
    }

    DT_FUNC (txt_file, zip_initialize)();

    segconf_calculate();

    max_lines_per_vb=0;

    dispatcher = 
        dispatcher_fan_out_task ("zip", txt_basename, 
                                 txt_file->redirected ? PROGRESS_MESSAGE : PROGRESS_PERCENT,
                                 txt_file->redirected ? "Compressing..." : NULL, 
                                 false, *is_last_file, z_closes_after_me, 
                                 flag.xthreads, prev_file_last_vb_i, 300,
                                 zip_prepare_one_vb_for_dispatching, 
                                 zip_compress_one_vb, 
                                 zip_complete_processing_one_vb);

    // go back and update some fields in the txt header's section header and genozip header -
    // only if we can go back - i.e. is a normal file, not redirected
    Digest single_component_digest = DIGEST_NONE;
    if (z_file && !flag.seg_only && !z_file->redirected && txt_header_header_pos >= 0) 
        success = zfile_update_txt_header_section_header (txt_header_header_pos, max_lines_per_vb, &single_component_digest);

    // write the BGZF section containing BGZF block sizes, if this txt file is compressed with BGZF
    bgzf_compress_bgzf_section();

    // if this a non-bound file, or the last component of a bound file - write the genozip header, random access and dictionaries
finish:
    if (!flag.rejects_coord)
        z_file->txt_disk_so_far_bind += (int64_t)txt_file->disk_so_far + (txt_file->codec==CODEC_BGZF)*BGZF_EOF_LEN;

    // reconstruction plan for sorting the data (per txt file)
    if (flag.sort && !flag.seg_only)
        linesorter_compress_recon_plan();

    if (z_closes_after_me && !flag.seg_only) {
        DT_FUNC (txt_file, zip_after_vbs)();
    
        zip_write_global_area (single_component_digest);

        if (chain_is_loaded && !Z_DT(DT_CHAIN)) vcf_liftover_display_lift_report();
    }

    zip_display_compression_ratio (flag.bind ? DIGEST_NONE : single_component_digest, z_closes_after_me); // Done for reference + final compression ratio calculation
    
    if (flag.md5 && flag.bind && z_file->num_txt_components_so_far > 1 && z_closes_after_me) 
        progress_concatenated_md5 (dt_name (z_file->data_type), digest_finalize (&z_file->digest_ctx_bound, "file:digest_ctx_bound"));

    z_file->disk_size              = z_file->disk_so_far;

    if (z_closes_after_me) 
        prev_file_first_vb_i = prev_file_last_vb_i = 0; // reset statics
    
    prev_file_first_vb_i = first_vb_i;
    dispatcher_finish (&dispatcher, &prev_file_last_vb_i);

    DT_FUNC (txt_file, zip_finalize)();
}
