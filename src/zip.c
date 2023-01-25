// ------------------------------------------------------------------
//   zip.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

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
#include "txtheader.h"
#include "threads.h"
#include "endianness.h"
#include "segconf.h"
#include "contigs.h"
#include "chrom.h"
#include "biopsy.h"
#include "dict_io.h"
#include "gencomp.h"

static Mutex wait_for_vb_1_mutex = {};

static void zip_display_compression_ratio (Digest md5)
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

        // reset for every set of bound files (we might have multiple sets of --pair)
        if (z_file->num_txts_so_far == 1) {
            comp_bytes_bind=0; 
            source_file_type = txt_file->type;
        }

        else if (source_file_type != txt_file->type) // heterogenous source file types
            source_file_type = UNKNOWN_FILE_TYPE;

        comp_bytes_bind += comp_bytes;

        if (z_file->z_closes_after_me) { 
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

    // in bound files, for the non-last components, we already printed "Done" above
    else if (flag.bind && !z_file->z_closes_after_me) {}

    // when compressing BAM report only ratio_vs_comp (compare to BGZF-compress BAM - we don't care about the underlying plain BAM)
    // Likewise, doesn't have a compression extension (eg .gz), even though it may actually be compressed eg .tbi (which is actually BGZF)
    else if (Z_DT(BAM) || (txt_file && file_get_codec_by_txt_ft (txt_file->data_type, txt_file->type, false) == CODEC_NONE)) 
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
static inline void zip_generate_one_b250 (VBlockP vb, ContextP ctx, uint32_t word_i, WordIndex word_index,
                                          uint8_t **next, 
                                          WordIndex *prev_word_index,  // in/out
                                          bool show)
{
    if (word_index >= 0) { // normal index

        bool one_up = (word_index == *prev_word_index + 1) && (word_i > 0);

        if (one_up) {
            (*next)[0] = BASE250_ONE_UP;
            (*next)++;
        }

        else if (ctx->b250_size == B250_BYTES_1) {
            (*next)[0] = word_index;
            (*next)++;
        }

        else if (word_index <= 2) { // note: we don't use MOST_FREQ in case BYTES_1
            (*next)[0] = BASE250_MOST_FREQ0 + word_index;
            (*next)++;
        }

        else if (ctx->b250_size == B250_BYTES_2) {
            (*next)[0] = word_index >> 8; // big endian. note: data not aligned to word boundary
            (*next)[1] = word_index & 0xff;
            (*next) += 2;
        }
            
        else if (ctx->b250_size == B250_BYTES_3) {
            (*next)[0] = (word_index >> 16); 
            (*next)[1] = (word_index >> 8 ) & 0xff;
            (*next)[2] = (word_index >> 0 ) & 0xff;
            (*next) += 3;
        }
            
        else { // 4 bytes
            (*next)[0] = (word_index >> 24); 
            (*next)[1] = (word_index >> 16) & 0xff;
            (*next)[2] = (word_index >> 8 ) & 0xff;
            (*next)[3] = (word_index >> 0 ) & 0xff;
            (*next) += 4;
        }

        if (show) bufprintf (vb, &vb->show_b250_buf, one_up ? "L%u:ONE_UP " : "L%u:%u ", word_i, word_index);
    }

    else if (word_index == WORD_INDEX_MISSING) {
        (*next)[0] = BASE250_MISSING_SF;
        (*next)++;
        if (show) bufprintf (vb, &vb->show_b250_buf, "L%u:MISSING ", word_i);
    }
    
    else if (word_index == WORD_INDEX_EMPTY) {
        (*next)[0] = BASE250_EMPTY_SF;
        (*next)++;
        if (show) bufprintf (vb, &vb->show_b250_buf, "L%u:EMPTY ", word_i);
    }

    else ABORT ("invalid word_index=%u", word_index);        

    *prev_word_index = word_index;
}

// here we generate the b250 data: we convert the b250 buffer data from an index into context->nodes
// to an index into the to-be-generated-in-piz context->word_index, encoded in base250.
// Note that the word indices have changed since segmentation (which is why we needed this intermediate step)
// because: 1. the dictionary got integrated into the global one - some values might have already been in the global
// dictionary thanks to other threads working on other VBs  
// 2. for the first VB, we sort the dictionary by frequency returns true if section should be dropped
static void zip_generate_b250 (VBlockP vb, Context *ctx)
{
    START_TIMER;

    ASSERT (ctx->dict_id.num, "tag_name=%s did_i=%u: ctx->dict_id=0 despite ctx->b250 containing data", ctx->tag_name, (unsigned)(ctx - vb->contexts));

    bool show = flag.show_b250 || dict_id_typeless (ctx->dict_id).num == flag.dict_id_show_one_b250.num;
    if (show) bufprintf (vb, &vb->show_b250_buf, "%s %s: ", VB_NAME, ctx->tag_name);

    ctx->b250.count = ctx->b250.len; // we move the number of words to count, as len will now contain the of bytes. used by ctx_update_stats()

    // case: all-the-same b250 survived dropping (in ctx_drop_all_the_same) - we just shorten it to one entry
    if (ctx->flags.all_the_same && ctx->b250.len > 1) { 
        if (flag.debug_generate) iprintf ("%s: %s.b250 is \"all_the_same\" - shortened b250 from len=%"PRIu64" to 1\n", VB_NAME, ctx->tag_name, ctx->b250.len);
        ctx->b250.len = 1;
    }
    else
        if (flag.debug_generate) iprintf ("%s: %s.b250 len=%"PRIu64" all_the_same=%u\n", VB_NAME, ctx->tag_name, ctx->b250.len, ctx->flags.all_the_same);

    // replace node indices with word indices
    WordIndex largest_wi=0;
    for_buf (WordIndex, ent, ctx->b250) 
        if (*ent >= 0) 
            if ((*ent = ctx_node_vb (ctx, *ent, NULL, NULL)->word_index) > largest_wi)
                largest_wi = *ent;
    
    // determine size of word_index elements
    ctx->b250_size = (largest_wi <= B250_MAX_WI_1BYTE ) ? B250_BYTES_1
                   : (largest_wi <= B250_MAX_WI_2BYTES) ? B250_BYTES_2
                   : (largest_wi <= B250_MAX_WI_3BYTES) ? B250_BYTES_3
                   :                                      B250_BYTES_4;
    
    ARRAY (WordIndex, b250, ctx->b250);
    ctx->b250.len = 0; // we are going to overwrite b250 with the converted indices

    // convert b250 from an array of node_index (4 bytes each) to word_index (1-4 bytes each)    
    WordIndex prev_word_index = WORD_INDEX_NONE; 
    uint8_t *next = B1ST8 (ctx->b250);
    for (uint64_t i=0; i < b250_len; i++) 
        zip_generate_one_b250 (vb, ctx, i, b250[i], &next, &prev_word_index, show);
    ctx->b250.len = BNUM (ctx->b250, next);

    // in case we are pairing b250 - check if entire section is identical - and reduce it to all-the-same FASTQ_SPECIAL_mate_lookup if it is
    if (ctx->pair_b250 && !ctx->flags.all_the_same &&
        ctx->pair.len == ctx->b250.len && !memcmp (ctx->b250.data, ctx->pair.data, ctx->b250.len)) {
        ctx_convert_generated_b250_to_mate_lookup (vb, ctx);
        if (flag.debug_generate) iprintf ("%s: %s.b250 is all_the_same=1 mate_lookup=1\n", VB_NAME, ctx->tag_name);
    }
    else
        ctx->pair_b250 = false;
    
    if (show) {
        bufprintf (vb, &vb->show_b250_buf, "%s", "\n");
        iprintf ("%.*s", STRfb(vb->show_b250_buf));
        buf_free (vb->show_b250_buf);
    }

    COPY_TIMER (zip_generate_b250); // codec_assign measures its own time

    codec_assign_best_codec (vb, ctx, NULL, SEC_B250);
}

static void zip_resize_local (VBlockP vb, Context *ctx)
{
    if (!ctx->local.len) {
        ctx->ltype = LT_UINT8;
        return;
    }
    
    ARRAY (int64_t, src, ctx->local);

    // search for the largest (stop if largest so far requires 32bit)
    int64_t largest = *B1ST (int64_t, ctx->local);
    int64_t smallest = largest;

    for (uint64_t i=1; i < src_len; i++) 
        if (src[i] > largest) largest = src[i];
        else if (src[i] < smallest) smallest = src[i];

    // in VCF FORMAT fields, a max_int value is reconstructed as . (see reconstruct_from_local_int), which is not applicable for dynamic size -
    // therefore we increment largest by 1 to ensure that no actual value is max_int
    if ((VB_DT(VCF) || VB_DT(BCF)) && dict_id_is_vcf_format_sf (ctx->dict_id)) largest++;

    static const LocalType test_ltypes[3][8] = { { LT_UINT8, LT_UINT16, LT_UINT32, LT_INT8, LT_INT16, LT_INT32, LT_UINT64, LT_INT64 },  // LT_DYN_INT
                                                 { LT_hex8, LT_hex16, LT_hex32, LT_hex64 },   // LT_DYN_INT_h
                                                 { LT_HEX8, LT_HEX16, LT_HEX32, LT_HEX64 } }; // LT_DYN_INT_H

    const LocalType *my_test_ltypes = test_ltypes[ctx->ltype - LT_DYN_INT];
    for (LocalType lt_i=0; lt_i < (ctx->ltype == LT_DYN_INT ? 8 : 4); lt_i++)
        if (smallest >= lt_desc[my_test_ltypes[lt_i]].min_int && largest <= lt_desc[my_test_ltypes[lt_i]].max_int) {
            ctx->ltype = my_test_ltypes[lt_i];
            break; // found
        }

    switch (ctx->ltype) {
        case LT_UINT8: case LT_hex8: case LT_HEX8:  {
            ARRAY (uint8_t, dst, ctx->local);
            for (uint64_t i=0; i < src_len; i++)
                dst[i] = (uint8_t)src[i];
            break;
        }

        case LT_UINT16: case LT_hex16: case LT_HEX16: {
            ARRAY (uint16_t, dst, ctx->local);
            for (uint64_t i=0; i < src_len; i++)
                dst[i] = BGEN16 ((uint16_t)src[i]);
            break;
        }

        case LT_UINT32: case LT_hex32: case LT_HEX32: {
            ARRAY (uint32_t, dst, ctx->local);
            for (uint64_t i=0; i < src_len; i++)
                dst[i] = BGEN32 ((uint32_t)src[i]);
            break;
        }

        case LT_INT8: {
            ARRAY (uint8_t, dst, ctx->local);
            for (uint64_t i=0; i < src_len; i++)
                dst[i] = INTERLACE(int8_t, src[i]);
            break;
        }

        case LT_INT16: {
            ARRAY (uint16_t, dst, ctx->local);
            for (uint64_t i=0; i < src_len; i++)
                dst[i] = BGEN16 (INTERLACE(int16_t, src[i]));
            break;
        }

        case LT_INT32: {
            ARRAY (uint32_t, dst, ctx->local);
            for (uint64_t i=0; i < src_len; i++)
                dst[i] = BGEN32 (INTERLACE(int32_t,src[i]));
            break;
        }

        case LT_UINT64: case LT_hex64: case LT_HEX64: {
            for (uint64_t i=0; i < src_len; i++)
                src[i] = BGEN64 (src[i]);
            break;
        }

        case LT_INT64: {\
            for (uint64_t i=0; i < src_len; i++)
                src[i] = BGEN64 (INTERLACE(int64_t, src[i]));
            break;
        }

        default:
            ABORT ("%s: Cannot find ltype for ctx=%s: value_range=[%"PRId64",%"PRId64"] ltype=%s",
                   VB_NAME, ctx->tag_name, smallest, largest, lt_name (ctx->ltype));
    }
}

// selects the smallest size (8, 16, 32) for the data, transposes, and BGENs
// note: in PIZ, these are untransposed in eg BGEN_transpose_u32_buf
static void zip_generate_transposed_local (VBlockP vb, Context *ctx)
{
    ARRAY (uint32_t, data, ctx->local);

    // get largest element (excluding 0xffffffff - representing a '.')
    uint32_t largest=0;
    for (uint64_t i=0; i < ctx->local.len; i++)
        if (data[i] > largest && data[i] != 0xffffffff) largest = data[i];

    if      (largest < 0xfe)   ctx->ltype = LT_UINT8_TR;  // -1 is reserved for "missing"
    else if (largest < 0xfffe) ctx->ltype = LT_UINT16_TR;
    
    buf_alloc (vb, &vb->scratch, 0, ctx->local.len * lt_desc[ctx->ltype].width, char, 1, "scratch");

    uint32_t cols = ctx->local.count;
    // we're restricted to 255 columns, because this number goes into uint8_t SectionHeaderCtx.param
    ASSERT (cols >= 0 && cols <= 255, "columns=%u out of range [1,255] in transposed matrix %s", cols, ctx->tag_name);

    if (!cols) cols = vcf_header_get_num_samples(); // 0 if not vcf/bcf (not restricted to 255)
    
    // case: matrix is not transposable - just BGEN it
    if (ctx->local.len32 % cols) {
        ctx->ltype = LT_UINT32; // not transposed

        for (unsigned i=0; i < ctx->local.len32; i++)
            data[i] = BGEN32 (data[i]);
            
        goto done;
    }

    uint32_t rows = ctx->local.len32 / cols;

    ctx->local_param/*xxx flags.copy_local_param*/ = true;
    /* xxx I don't see where col goes into param??? need to test */
    
    for (uint32_t r=0; r < rows; r++) 
        for (uint32_t c=0; c < cols; c++) {

            uint32_t value = data[r * cols + c];

            switch (ctx->ltype) { // note: the casting aslo correctly converts 0xffffffff to eg 0xff
                case LT_UINT8_TR  : *B8 ( vb->scratch, c * rows + r) =         (uint8_t)value;   break;
                case LT_UINT16_TR : *B16 (vb->scratch, c * rows + r) = BGEN16 ((uint16_t)value); break;
                case LT_UINT32_TR : *B32 (vb->scratch, c * rows + r) = BGEN32 (value);           break;
                default: ABORT ("Bad ltype=%s", lt_name (ctx->ltype));
            }
        }

    vb->scratch.len = ctx->local.len;
    buf_copy_do (vb, &ctx->local, &vb->scratch, lt_desc[ctx->ltype].width, 0, 0, __FUNCTION__,__LINE__, "contexts->local"); // copy and not move, so we can keep local's memory for next vb

done:
    buf_free (vb->scratch);
}

// after segging - if any context appears to contain only singleton snips (eg a unique ID),
// we move it to local instead of needlessly cluttering the global dictionary
static void zip_handle_unique_words_ctxs (VBlockP vb)
{
    START_TIMER;

    for_ctx {
        if (!ctx->nodes.len || ctx->nodes.len != ctx->b250.len) continue; // check that all words are unique (and new to this vb)
        if ((VB_DT(VCF) || VB_DT(BCF)) && dict_id_is_vcf_format_sf (ctx->dict_id)) continue; // this doesn't work for FORMAT fields
        if (ctx->nodes.len < vb->lines.len / 5) continue; // don't bother if this is a rare field less than 20% of the lines
        if (buf_is_alloc (&ctx->local))     continue; // skip if we are already using local to optimize in some other way

        // don't move to local if its on the list of special dict_ids that are always in dict (because local is used for something else - eg pos or id data)
        if (!ctx_can_have_singletons (ctx) ||
            ctx->b250.len == 1) continue; // handle with all_the_same rather than singleton
         
        buf_move (vb, &ctx->local, vb, &ctx->dict);
        buf_free (ctx->nodes);
        buf_free (ctx->b250);
        ctx->flags.all_the_same = false;
    }

    COPY_TIMER (zip_handle_unique_words_ctxs);
}

static void zip_generate_local (VBlockP vb, ContextP ctx)
{
    START_TIMER;
    
    ASSERT (ctx->dict_id.num, "tag_name=%s did_i=%u: ctx->dict_id=0 despite ctx->local containing data", ctx->tag_name, (unsigned)(ctx - vb->contexts));

    bool resizeable = ctx->ltype == LT_DYN_INT || ctx->ltype == LT_DYN_INT_H || ctx->ltype == LT_DYN_INT_h; 
    if (resizeable) 
        zip_resize_local (vb, ctx);

    else {
        // case: local is LTEN (instead of native endianity) and machine is BGEN, so BGEN_*_buf ^ above did nothing.     
        bool need_lten = (ctx->local_is_lten && !flag.is_lten);
        
        switch (ctx->ltype) {
            case LT_BITMAP    : LTEN_bits ((BitsP)&ctx->local);    
                                ctx->local.prm8[0] = ((uint8_t)64 - (uint8_t)(ctx->local.nbits % 64)) % (uint8_t)64;
                                ctx->local_param   = true;
                                break;

            case LT_UINT32_TR : zip_generate_transposed_local (vb, ctx);            
                                break;

            case LT_FLOAT32   : 

            case LT_UINT32    : if (need_lten) LTEN_u32_buf (&ctx->local, NULL);   
                                else           BGEN_u32_buf (&ctx->local, NULL);    
                                break;

            case LT_UINT16    : if (need_lten) LTEN_u16_buf (&ctx->local, NULL);   
                                else           BGEN_u16_buf (&ctx->local, NULL);    
                                break;

            case LT_FLOAT64   :
            case LT_UINT64    : if (need_lten) LTEN_u64_buf (&ctx->local, NULL);   
                                else           BGEN_u64_buf (&ctx->local, NULL);    
                                break;

            case LT_INT8      : interlace_d8_buf  (&ctx->local, NULL);              
                                break;

            case LT_INT16     : if (need_lten) LTEN_interlace_d16_buf (&ctx->local, NULL);   
                                else           BGEN_interlace_d16_buf (&ctx->local, NULL);    
                                break;

            case LT_INT32     : if (need_lten) LTEN_interlace_d32_buf (&ctx->local, NULL);   
                                else           BGEN_interlace_d32_buf (&ctx->local, NULL);    
                                break;

            case LT_INT64     : if (need_lten) LTEN_interlace_d64_buf (&ctx->local, NULL);   
                                else           BGEN_interlace_d64_buf (&ctx->local, NULL);    
                                break;

            default           : break;        
        }
    }

    COPY_TIMER (zip_generate_local); // codec_assign measures its own time

    codec_assign_best_codec (vb, ctx, NULL, SEC_LOCAL);

    if (flag.debug_generate) 
        iprintf ("%s: %s.local ltype=%s%s len=%"PRIu64" codec=%s\n", VB_NAME, ctx->tag_name, 
                 lt_name (ctx->ltype), resizeable ? " (resized)" : "", ctx->local.len, codec_name(ctx->lcodec));
}

// generate & write b250 data for all contexts - do them in random order, to reduce the chance of multiple doing codec_assign_best_codec for the same context at the same time
// VBs doing codec_assign_best_codec at the same, so that they can benefit from pre-assiged codecs
void zip_compress_all_contexts_b250 (VBlockP vb)
{
    START_TIMER;
    threads_log_by_vb (vb, "zip", "START COMPRESSING B250", 0);
    
    // arrays of all contexts in this VB
    ContextP ctxs[vb->num_contexts];
    for (Did did_i=0; did_i < vb->num_contexts; did_i++) ctxs[did_i] = CTX(did_i);

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
        
        if (HAS_DEBUG_SEG(ctx)) iprintf ("zip_compress_all_contexts_b250: vb=%s %s: B250.len=%"PRIu64" NODES.len=%"PRIu64"\n", 
                                         VB_NAME, ctx->tag_name, ctx->b250.len, ctx->nodes.len);

        START_TIMER; // for compressor_time

        zfile_compress_b250_data (vb, ctx);

        if (flag.show_time) 
            ctx->compressor_time += CHECK_TIMER; // sum b250 and local

        ctx->b250_compressed = true;
    }

    COPY_TIMER (zip_compress_ctxs); // same profiler for b250 and local as we breakdown by ctx underneath it
}

// generate & write local data for all contexts - in random order, to reduce the chance of multiple doing codec_assign_best_codec for the same context at the same time
static void zip_compress_all_contexts_local (VBlockP vb)
{
    START_TIMER;
    threads_log_by_vb (vb, "zip", "START COMPRESSING LOCAL", 0);

    // first we handle local_dep=0 then local_dep=1 and finally local_dep=2
    for (int dep_level=DEP_L0 ; dep_level < NUM_LOCAL_DEPENDENCY_LEVELS; dep_level++) {

        // initialize list of contexts at this dependency level that need compression
        ContextP ctxs[vb->num_contexts];
        unsigned num_ctxs=0;
        for_ctx
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

            if (HAS_DEBUG_SEG(ctx)) iprintf ("%s: zip_compress_all_contexts_local: %s: LOCAL.len=%"PRIu64" LOCAL.param=%"PRIu64"\n", 
                                            VB_NAME, ctx->tag_name, ctx->local.len, ctx->local.param);

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

void zip_init_vb (VBlockP vb)
{
    vb->recon_size = vb->txt_data.len; // initial value. it may change if --optimize / --chain are used, or if dual coordintes - for the other coordinate
    vb->txt_size   = vb->txt_data.len; // this copy doesn't change with --optimize / --chain.

    if (DTPT(zip_init_vb)) DTPT(zip_init_vb)(vb); // data-type specific initialization of the VB    
}

// called by main thread after VB has completed processing
static void zip_update_txt_counters (VBlockP vb)
{
    // note: in case of an FASTQ with flag.optimize_DESC or VCF with add_line_numbers, we already updated this in *_zip_init_vb
    if (!(flag.optimize_DESC && VB_DT(FASTQ)) &&
        !(flag.add_line_numbers && (VB_DT(VCF) || VB_DT(BCF))))
        txt_file->num_lines += vb->lines.len; // lines in this txt file

    // counters of data AS IT APPEARS IN THE TXT FILE
    if (!(z_is_dvcf && vb->comp_i))  // in DVCF, the generated to components contain copied data so we don't need to count it again
        z_file->num_lines += vb->lines.len; // lines in all bound files in this z_file

    if (vb->comp_i==COMP_MAIN || !z_is_dvcf) { // DVCF, vb->txt_size of the main file, as read from disk, contains all the data
        z_file->txt_data_so_far_single_0 += (int64_t)vb->txt_size;  // length of data before any modifications
        z_file->txt_data_so_far_bind_0   += (int64_t)vb->txt_size;
    }

    // counter of data FOR PROGRESS DISPLAY
    z_file->txt_data_so_far_single += (int64_t)vb->txt_size;   

    // counter of data in DEFAULT RECONSTRUCTION (For DVCF, this is corrected in vcf_zip_update_txt_counters)
    z_file->txt_data_so_far_bind += vb->recon_size;

    // per-component data for stats
    z_file->txt_data_so_far_bind_0_comp[vb->comp_i] += (int64_t)vb->txt_size;
    
    // note: in case of SAM gencomp, MAIN, we add recon_size - assuming the discrepency vs txt_data.len
    // is only due to lines being deported to gencomp 
    z_file->txt_data_so_far_bind_comp[vb->comp_i] += 
        (z_sam_gencomp && vb->comp_i==SAM_COMP_MAIN) ? vb->recon_size : vb->txt_data.len;

    // add up context compress time
    if (flag.show_time)
        ctx_add_compressor_time_to_zf_ctx (vb);
}

// write all the sections at the end of the file, after all VB stuff has been written
static void zip_write_global_area (void)
{
    // if we're making a reference, we need the RA data to populate the reference section chrome/first/last_pos ahead of ref_compress_ref
    random_access_finalize_entries (&z_file->ra_buf); // sort RA, update entries that don't yet have a chrom_index
    random_access_finalize_entries (&z_file->ra_buf_luft); 

    // if we used the aligner with REF_EXT_STORE, we make sure all the CHROMs referenced are in the CHROM context, so
    // as SEC_REF_CONTIGS refers to them. We do this by seeing which contigs have any bit set in is_set.
    // note: in REF_EXTERNAL we don't use is_set, so we populate all contigs in zip_initialize
    if (flag.aligner_available && IS_REF_EXT_STORE)
        ref_contigs_populate_aligned_chroms();

    dict_io_compress_dictionaries(); 

    ctx_compress_counts();

    // store a mapping of the file's chroms to the reference's contigs, if they are any different
    // note: not needed in REF_EXT_STORE, as we convert the stored ref_contigs to use chrom_index of the file's CHROM
    if (IS_REF_EXTERNAL) 
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
    BufferP dict_id_aliases_buf = dict_id_create_aliases_buf();
    if (dict_id_aliases_buf->len) zfile_compress_section_data (evb, SEC_DICT_ID_ALIASES, dict_id_aliases_buf);

    // SAM/BAM: we don't compress RANDOM_ACCESS for non-sorted in --best (it can be very big, and we want to minimize file size in --best), 
    if (!flag.best || !(Z_DT(BAM) || Z_DT(SAM)) || segconf.is_sorted) {
        // if this data has random access (i.e. it has chrom and pos), compress all random access records into evb->z_data
        Codec codec = random_access_compress (&z_file->ra_buf, SEC_RANDOM_ACCESS, CODEC_UNKNOWN, 0, flag.show_index ? RA_MSG_PRIM : NULL);
    
        if (z_is_dvcf)
            random_access_compress (&z_file->ra_buf_luft, SEC_RANDOM_ACCESS, codec, 1, flag.show_index ? RA_MSG_LUFT : NULL);

        if (store_ref) 
            random_access_compress (ref_get_stored_ra (gref), SEC_REF_RAND_ACC, codec, 0, flag.show_ref_index ? RA_MSG_REF : NULL);
    }

    stats_generate();

    // compress genozip header (including its payload sectionlist and footer) into evb->z_data
    zfile_compress_genozip_header();    

    if (DTPZ(zip_free_end_of_z)) DTPZ(zip_free_end_of_z)();
}

// entry point for ZIP compute thread
static void zip_compress_one_vb (VBlockP vb)
{
    START_TIMER; 

    if (flag.biopsy) goto done; // we're just taking a biopsy of the txt data, so no need to actually compress

    // if the txt file is compressed with BGZF, we uncompress now, in the compute thread
    if (txt_file->codec == CODEC_BGZF && flag.pair != PAIR_READ_2) 
        bgzf_uncompress_vb (vb);    // some of the blocks might already have been decompressed while reading - we decompress the remaining

    // calculate the digest contribution of this VB, and the digest snapshot of this VB
    if (!flag.make_reference && !flag.data_modified) 
        digest_one_vb (vb, true, NULL); // serializes VBs in order

    // allocate memory for the final compressed data of this vb. allocate 33% of the
    // vb size on the original file - this is normally enough. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, 0, vb->txt_size / 3, char, CTX_GROWTH, "z_data");

    // clone global dictionaries while granted exclusive access to the global dictionaries
    if (flag.pair != PAIR_READ_2) // in case of PAIR_READ_2, we already cloned in zip_one_file
        ctx_clone (vb);

    // split each line in this VB to its components
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

    bool need_compress = !flag.make_reference && !flag.seg_only;

    if (vb->vblock_i != 1) {
        if (need_compress) 
            // while vb_i=1 is busy merging, other VBs can handle local
            zip_compress_all_contexts_local (vb); // not yet locals that consist of singletons transferred from dict to local in ctx_merge_in_vb_ctx (these will have len=0 at this point)
        
        // let vb_i=1 merge first, as it has the sorted dictionaries, other vbs can go in arbitrary order. 
        mutex_wait (wait_for_vb_1_mutex, true);
    }

    dispatcher_increment_progress ("compress1", vb->txt_size / 2); // 1/2 compression done

    // merge new words added in this vb into the z_file.contexts, ahead of zip_generate_b250().
    // writing indices based on the merged dictionaries. dictionaries are compressed. 
    // all this is done while holding exclusive access to the z_file dictionaries.
    // note: vb>=2 will block here, until vb=1 is completed
    threads_log_by_vb (vb, "zip", "START MERGE", 0);

    // for --make-reference we serialize merging by VB, so that contigs get their word_index in the order of the reference file
    if (flag.make_reference) serializer_lock (make_ref_merge_serializer, vb->vblock_i);

    ctx_merge_in_vb_ctx(vb);

    if (flag.make_reference) serializer_unlock (make_ref_merge_serializer);

    if (vb->vblock_i == 1) 
        mutex_unlock (wait_for_vb_1_mutex); 

    if (need_compress) {
        zip_compress_all_contexts_local (vb); // for vb=1 - all locals ; for vb>1 - locals which consist of singletons set in ctx_merge_in_vb_ctx (other locals were already compressed above)
        zip_compress_all_contexts_b250 (vb);
    }

    dispatcher_increment_progress ("compress2", 1 - vb->txt_size / 2); // 1/2 compression done

    // merge in random access - IF it is used
    random_access_merge_in_vb (vb, 0);
    random_access_merge_in_vb (vb, 1);
    
    // compress data-type specific sections
    DT_FUNC (vb, zip_after_compress)(vb);

    // tell dispatcher this thread is done and can be joined.
done:
    vb_set_is_processed (vb); 

    COPY_TIMER (compute);
}

// data sent through dispatcher fan out functions - to do: make this an opaque struct
static uint32_t prev_file_first_vb_i=0, prev_file_last_vb_i=0; // used if we're binding files - the vblock_i will continue from one file to the next
static uint32_t max_lines_per_vb; // (resets in every file)

// main thread: returns true if successfully prepared a vb 
static void zip_prepare_one_vb_for_dispatching (VBlockP vb)
{
    // if we're compressing the 2nd file in a fastq pair (with --pair) - look back at the z_file data
    // and copy the data we need for this vb. note: we need to do this before txtfile_read_vblock as
    // we need the num_lines of the pair file
    if (flag.pair == PAIR_READ_2) {
        // note: normally we clone in the compute thread, because it might wait on mutex, but in this
        // case we need to clone (i.e. create all contexts before we can read the pair file data)
        ctx_clone (vb); 
    
        uint32_t pair_vb_i = prev_file_first_vb_i + (vb->vblock_i - prev_file_last_vb_i - 1);
        
        if (pair_vb_i > prev_file_last_vb_i || // false if there is no vb with vb_i in the previous file
            !fastq_read_pair_1_data (vb, pair_vb_i, false)) { // read here, decompressed in fastq_seg_initialize
            vb->dispatch = DATA_EXHAUSTED;
            return;
        }
    }

    // case: we have out-of-band txt_data waiting (for generated components) - compress this data first,
    // before reading more data from the txt_file
    if (gencomp_get_txt_data(vb)) 
        goto dispatch;

    else {        
        txtfile_read_vblock (vb);

        // initializations after reading the first vb and before running any compute thread
        if (vb->vblock_i == 1) { 
            // we lock, so vb>=2 will block on merge, until vb=1 has completed its merge and unlocked it in zip_compress_one_vb()
            mutex_destroy (wait_for_vb_1_mutex); // destroy mutex of previous file (it will still be locked if that file had no VBs)
            mutex_initialize (wait_for_vb_1_mutex);
            mutex_lock (wait_for_vb_1_mutex); 
        }

        // --head advanced option in ZIP cuts a certain number of first lines from vb=1, and discards other VBs
        else if (flag.lines_last != NO_LINE) {
            vb->dispatch = DATA_EXHAUSTED;
            return;
        }

        if (vb->txt_data.len) 
            goto dispatch;

        else if (gencomp_am_i_expecting_more_txt_data()) // more data might be coming from MAIN VBs currently computing
            vb->dispatch = MORE_DATA_MIGHT_COME;

        else 
            vb->dispatch = DATA_EXHAUSTED;

        // error if stdin is empty - can happen only when redirecting eg "cat empty-file|./genozip -" (we test for empty regular files in main_genozip)
        ASSINP0 (vb->vblock_i > 1 || txt_file->txt_data_so_far_single /* txt header data */, 
                 "Error: Cannot compress stdin data because its size is 0");
        
        return;
    }

dispatch:
    vb->dispatch = READY_TO_COMPUTE;
    txt_file->num_vbs_dispatched++;

    if (vb->comp_i == COMP_MAIN) // note: we only update the MAIN comp from here, gen comps are updated
        gencomp_a_main_vb_has_been_dispatched();
}

// called main thread, as VBs complete (might be out-of-order)
static void zip_complete_processing_one_vb (VBlockP vb)
{
    DT_FUNC (vb, zip_after_compute)(vb);

    // update z_data in memory (its not written to disk yet)
    zfile_update_compressed_vb_header (vb); 

    max_lines_per_vb = MAX_(max_lines_per_vb, vb->lines.len);

    if (!flag.make_reference && !flag.seg_only)
        zfile_output_processed_vb (vb);
    
    zip_update_txt_counters (vb);

    dispatcher_increment_progress ("z_write", vb->txt_size); // writing done.

    z_file->num_vbs++;
    txt_file->num_vbs++;
}

// this is the main dispatcher function. It first processes the txt header, then proceeds to read 
// a VB from the input file and send it off to a thread for computation. When the thread
// completes, this function proceeds to write the output to the output file. It can dispatch
// several threads in parallel.
void zip_one_file (rom txt_basename, 
                   bool is_last_user_txt_file)  // the last user-specified txt file in this execution
{
    Dispatcher dispatcher = 0;
    dispatcher_start_wallclock();

    z_file->txt_data_so_far_single = 0;
    evb->z_data.len                = 0;
    evb->z_next_header_i           = 0;
    
    // we calculate digest for each component seperately, stored in SectionHeaderTxtHeader (always 0 for generated components, or if modified)
    if (gencomp_comp_eligible_for_digest(NULL)) // if generated component - keep digest to display in progress after the last component
        z_file->digest_ctx = DIGEST_CONTEXT_NONE;

    if (!flag.bind || flag.zip_comp_i == COMP_MAIN) {
        prev_file_first_vb_i = prev_file_last_vb_i = 0; // reset if we're not binding
        
        segconf_initialize(); // before txtheader 
    }

    uint32_t first_vb_i = prev_file_last_vb_i + 1;

    // initalize pre-defined ctxs after reading header 
    // note: in case of GENERIC, generic_is_header_done may change the data type and re-initialize the contexts
    if (z_file->num_txts_so_far == 0)  // first component of this z_file 
        ctx_initialize_predefined_ctxs (z_file->contexts, txt_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_contexts);

    // read the txt header, assign the global variables, and write the compressed header to the GENOZIP file
    int64_t txt_header_offset = -1;
    int64_t txt_header_len = txtheader_zip_read_and_compress (&txt_header_offset, flag.zip_comp_i); // also increments z_file->num_txts_so_far

    bool success = (txt_header_len >= -1);
    if (!success) goto finish; // eg 2nd+ VCF file cannot bind, because of different sample names

    DT_FUNC (txt_file, zip_initialize)();

    segconf_calculate();

    DT_FUNC (txt_file, zip_after_segconf)();

    max_lines_per_vb=0;

    static uint64_t target_progress=0;
    if (flag.pair != 2) { // note: if 2nd of a FASTQ file pair - we leave the target as it was in the first file as seggable_size is not calculated for the 2nd file
        int64_t est_seggable_size = txtfile_get_seggable_size(); 
    
        target_progress = est_seggable_size * 3 // read, seg, compress
                        + (!flag.make_reference && !flag.seg_only) * est_seggable_size; // write
    }

    dispatcher = 
        dispatcher_fan_out_task ("zip", txt_basename, 
                                 target_progress, // target progress: 1 for each read, compute, write
                                 target_progress ? NULL : "Compressing...",
                                 !flag.make_reference && !z_is_dvcf,   // allow callbacks to zip_complete_processing_one_vb not in order of VBs (not allowed for make-reference as contigs need to be in consistent order; not supported yet for DVCF)
                                 false,           // not test mode
                                 flag.xthreads, prev_file_last_vb_i, 5000,
                                 zip_prepare_one_vb_for_dispatching, 
                                 zip_compress_one_vb, 
                                 zip_complete_processing_one_vb);

    dispatcher_increment_progress ("txt_header", txt_header_len * 3); // txt_header was already read, computed and written

    // go back and update some fields in the txt header's section header and genozip header 
    if (txt_header_offset >= 0) // note: this will be -1 if we didn't write a SEC_TXT_HEADER section for any reason
        success = zfile_update_txt_header_section_header (txt_header_offset, max_lines_per_vb);

    // write the BGZF section containing BGZF block sizes, if this txt file is compressed with BGZF
    bgzf_compress_bgzf_section();

    // if this a non-bound file, or the last component of a bound file - write the genozip header, random access and dictionaries
finish:   
    if (!flag.zip_comp_i || Z_DT(FASTQ)) // in paired FASTQ we count found files
        z_file->txt_disk_so_far_bind += (int64_t)txt_file->disk_so_far + (txt_file->codec==CODEC_BGZF)*BGZF_EOF_LEN;

    // reconstruction plan (for VCF - for DVCF or --sort, for SAM - re-integrate supp/secondary alignments)
    if (!flag.seg_only && DTPZ(generate_recon_plan)) 
        DTPZ(generate_recon_plan)(); // should set z_file->z_closes_after_me if we need to close after this component after all

    if (z_file->z_closes_after_me && !flag.seg_only) { // note: for SAM, z_closes_after_me might be updated in sam_zip_generate_recon_plan
        DT_FUNC (txt_file, zip_after_vbs)();
    
        zip_write_global_area();

        if (chain_is_loaded && !Z_DT(CHAIN)) vcf_liftover_display_lift_report();
    }

    zip_display_compression_ratio (digest_snapshot (&z_file->digest_ctx, NULL)); // Done for reference + final compression ratio calculation
    
    if (flag.md5 && flag.bind && z_file->z_closes_after_me &&
        ((flag.bind == BIND_FQ_PAIR && z_file->num_txts_so_far == 2) ||
         ((flag.bind == BIND_DVCF || flag.bind == BIND_SAM) && z_file->num_txts_so_far == 3)))
        progress_concatenated_md5 (dt_name (z_file->data_type), digest_snapshot (&z_file->digest_ctx, "file"));

    z_file->disk_size = z_file->disk_so_far;

    prev_file_first_vb_i = first_vb_i;
    dispatcher_finish (&dispatcher, &prev_file_last_vb_i, 
                       z_file->z_closes_after_me && !is_last_user_txt_file,
                       flag.show_memory && z_file->z_closes_after_me && is_last_user_txt_file); // show memory

    if (!z_file->z_closes_after_me)
        ctx_reset_codec_commits(); //xxx

    DT_FUNC (txt_file, zip_finalize)(is_last_user_txt_file);
}
