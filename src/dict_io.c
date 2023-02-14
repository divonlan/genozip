// ------------------------------------------------------------------
//   dict_io.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include "genozip.h"
#include "sections.h"
#include "vblock.h"
#include "context.h"
#include "zfile.h"
#include "endianness.h"
#include "file.h"
#include "dict_id.h"
#include "strings.h"
#include "flags.h"
#include "buffer.h"
#include "dispatcher.h"
#include "piz.h"
#include "compressor.h"

// -------------------------------------
// ZIP: Assign codecs to dictionaries
// -------------------------------------

static ContextP next_ctx = 0;

static void dict_io_prepare_for_assign_codec (VBlockP vb)
{
    // get next non-too-small dict
    while (next_ctx->dcodec != CODEC_UNKNOWN && 
           next_ctx < ZCTX(z_file->num_contexts)) next_ctx++;

    if (next_ctx < ZCTX(z_file->num_contexts)) {
        vb->fragment_ctx = next_ctx;
        vb->dispatch = READY_TO_COMPUTE;
        next_ctx++;
    }
}

static void dict_io_assign_codec_one_dict (VBlockP vb)
{
    if (!vb->fragment_ctx->dcodec) // not small dict already assigned in dict_io_assign_codecs
        vb->fragment_ctx->dcodec = codec_assign_best_codec (vb, vb->fragment_ctx, NULL, SEC_DICT);
    
    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// called by main thread in zip_write_global_area
void dict_io_assign_codecs (void)
{
    START_TIMER;

    next_ctx = ZCTX(0);
    
    // handle some dictionaries here, so save on thread creation for trival dictionaries
    for_zctx {
        // assign CODEC_NONE to all the to-small-to-compress dictionaries, 
        if (zctx->dict.len < MIN_LEN_FOR_COMPRESSION) 
            zctx->dcodec = CODEC_NONE;

        // assign CODEC_ARTB to dictionaries under 1KB (unless --best)
        else if (!flag.best && zctx->dict.len < 1024) 
            zctx->dcodec = CODEC_ARITH8;
    }

    dispatcher_fan_out_task ("assign_dict_codecs", NULL, 0, "Writing dictionaries...", true, false, false, 0, 20000,
                             dict_io_prepare_for_assign_codec, 
                             dict_io_assign_codec_one_dict, 
                             NO_CALLBACK);

    COPY_TIMER_VB (evb, dict_io_assign_codecs);
}

// -------------------------------------
// ZIP: Compress and output dictionaries
// -------------------------------------

static Context *frag_ctx;
static const CtxNode *frag_next_node;
static unsigned frag_size = 0;

// compress the dictionary fragment - either an entire dict, or divide it to fragments if large to allow multi-threaded
// compression and decompression
static void dict_io_prepare_for_compress (VBlockP vb)
{
    while (frag_ctx < ZCTX(z_file->num_contexts)) {

        if (!frag_next_node) {
            if (!frag_ctx->nodes.len ||
                frag_ctx->please_remove_dict) { // this context belongs to a method that lost the test and is not used in this file
                frag_ctx++;
                continue; // unused context
            }

            frag_next_node = B1ST (const CtxNode, frag_ctx->nodes);

            ASSERT (frag_next_node->char_index + frag_next_node->snip_len <= frag_ctx->dict.len, 
                    "Corrupt nodes in ctx=%.8s did_i=%u", frag_ctx->tag_name, (int)(frag_ctx - z_file->contexts));

            ASSERT (frag_ctx->dict_id.num, "dict_id=0 for context \"%s\" did_i=%u", frag_ctx->tag_name, (unsigned)(frag_ctx - ZCTX(0)));

            // frag_size is set by default 1MB, but must be at least double the longest snip, or it will cause
            // mis-calculation of size_upper_bound in dict_io_read_one_vb
            frag_size = 1 MB;
            for (const CtxNode *node=B1ST (CtxNode, frag_ctx->nodes); node <= BLST (CtxNode, frag_ctx->nodes); node++)
                if (node->snip_len * 2 > frag_size)
                    frag_size = 2 * roundup2pow (node->snip_len); // must be power of 2 for dict_io_read_one_vb
        }

        vb->fragment_ctx   = frag_ctx;
        vb->fragment_start = Bc (frag_ctx->dict, frag_next_node->char_index);

        while (frag_next_node < BAFT (CtxNode, frag_ctx->nodes) && 
               vb->fragment_len + frag_next_node->snip_len + 1 < frag_size) {

            vb->fragment_len += frag_next_node->snip_len + 1;
            vb->fragment_num_words++;
            frag_next_node++;
        }

        if (frag_next_node == BAFT (CtxNode, frag_ctx->nodes)) {
            frag_ctx++;
            frag_next_node = NULL;
            frag_size = 0;
        }

        if (vb->fragment_len) {
            vb->dispatch = READY_TO_COMPUTE;
            break;
        }
    }
}

static void dict_io_compress_one_fragment (VBlockP vb)
{
    START_TIMER;

    SectionHeaderDictionary header = (SectionHeaderDictionary){ 
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .section_type          = SEC_DICT,
        .data_uncompressed_len = BGEN32 (vb->fragment_len),
        .compressed_offset     = BGEN32 (sizeof(SectionHeaderDictionary)),
        .codec                 = vb->fragment_ctx->dcodec,
        .vblock_i              = BGEN32 (vb->vblock_i),
        .num_snips             = BGEN32 (vb->fragment_num_words),
        .dict_id               = vb->fragment_ctx->dict_id
    };

    if (flag.show_dict) 
        iprintf ("%s (vb=%u, did=%u, num_snips=%u)\n", 
                 vb->fragment_ctx->tag_name, vb->vblock_i, vb->fragment_ctx->did_i, vb->fragment_num_words);
    
    if (dict_id_is_show (vb->fragment_ctx->dict_id))
        str_print_dict (info_stream, vb->fragment_start, vb->fragment_len, flag.show_dict, false);

    if (flag.list_chroms && vb->fragment_ctx->did_i == CHROM)
        str_print_dict (info_stream, vb->fragment_start, vb->fragment_len, false, VB_DT(SAM) || VB_DT(BAM));

    if (flag.show_time) codec_show_time (vb, st_name (SEC_DICT), vb->fragment_ctx->tag_name, vb->fragment_ctx->dcodec);

    comp_compress (vb, vb->fragment_ctx, &vb->z_data, (SectionHeader*)&header, vb->fragment_start, NO_CALLBACK, "SEC_DICT");

    COPY_TIMER (dict_io_compress_one_fragment)    

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// called by main thread in zip_write_global_area
void dict_io_compress_dictionaries (void)
{
    START_TIMER;

    frag_ctx = ZCTX(0);
    frag_next_node = NULL;

    dict_io_assign_codecs(); // assign codecs to all contexts' dicts

    dispatcher_fan_out_task ("compress_dicts", NULL, 0, "Writing dictionaries...", false, false, false, 0, 20000,
                             dict_io_prepare_for_compress, 
                             dict_io_compress_one_fragment, 
                             zfile_output_processed_vb);

    COPY_TIMER_VB (evb, dict_io_compress_dictionaries);
}

// -------------------------------------
// PIZ: Read and decompress dictionaries
// -------------------------------------
static Section dict_sec = NULL; 
static Context *dict_ctx;

static void dict_io_read_one_vb (VBlockP vb)
{
    vb->preprocessing = flag.preprocessing; // for debug_read_ctxs printing

    Section old_dict_sec = dict_sec;
    if (!sections_next_sec (&dict_sec, SEC_DICT)) 
        return; // we're done - no more SEC_DICT sections

    // while we could easily read non-consecutive sections, this is not expected to happen and may indicate multiple calls to zip_write_global_area
    ASSERT0 (!old_dict_sec || (old_dict_sec+1 == dict_sec), "Unexpectedly, not all SEC_DICT sections are consecutive in the Genozip file");

    // create context if if section is skipped, for containters to work (skipping a section should be mirror in 
    // a container filter)
    bool new_ctx = (!dict_ctx || dict_sec->dict_id.num != dict_ctx->dict_id.num);
    if (new_ctx)
        dict_ctx = ctx_get_ctx_do (z_file->contexts, z_file->data_type, z_file->d2d_map, &z_file->num_contexts, dict_sec->dict_id, 0, 0);

    if (piz_is_skip_section (SEC_DICT, COMP_NONE, dict_sec->dict_id, SKIP_PURPOSE_RECON)) {
        if (flag.debug_read_ctxs)
            iprintf ("%c Skipped loading DICT/%u %s.dict\n", sections_read_prefix, vb->vblock_i, dict_ctx->tag_name);
        goto done;
    }
    dict_ctx->is_loaded = true; // not skipped
    
    int32_t offset = zfile_read_section (z_file, vb, dict_sec->vblock_i, &vb->z_data, "z_data", SEC_DICT, dict_sec);    
    SectionHeaderDictionary *header = 
        (offset != SECTION_SKIPPED) ? (SectionHeaderDictionary *)vb->z_data.data : NULL;

    vb->fragment_len = header ? BGEN32 (header->data_uncompressed_len) : 0;

    ASSERT (!header || header->dict_id.num == dict_sec->dict_id.num, "Expecting dictionary fragment with DictId=%s but found one with DictId=%s",
            dis_dict_id (dict_sec->dict_id).s, dis_dict_id (header->dict_id).s);

    // new context
    // note: in v9+ same-dict fragments are consecutive in the file, and all but the last are frag_size or a bit less, allowing pre-allocation
    if (header && new_ctx && VER(9)) {
        unsigned num_fragments=0; 
        for (Section sec=dict_sec; sec->dict_id.num == dict_ctx->dict_id.num; sec++) num_fragments++;

        // get size: for multi-fragment dictionaries, first fragment will be at or less than a power of 2, but more than the previous power of two.
        // this allows us to calculate the frag_size with which this dictionary was compressed and hence an upper bound on the size
        uint64_t size_upper_bound = (num_fragments == 1) ? vb->fragment_len : ((uint64_t)roundup2pow (vb->fragment_len) * (uint64_t)num_fragments);
        
        buf_alloc (evb, &dict_ctx->dict, 0, size_upper_bound, char, 0, "contexts->dict");
        dict_ctx->dict.prm32[0] = num_fragments;     // for error reporting
        dict_ctx->dict.prm32[1] = vb->fragment_len;

        buf_set_overlayable (&dict_ctx->dict);
    }

    // when pizzing a v8 file, we run in single-thread since we need to do the following dictionary enlargement with which fragment
    if (!VER(9)) {
        buf_alloc (evb, &dict_ctx->dict, vb->fragment_len, 0, char, 0, "contexts->dict");
        buf_set_overlayable (&dict_ctx->dict);
    }

    if (header) {
        vb->fragment_ctx         = dict_ctx;
        vb->fragment_start       = Bc (dict_ctx->dict, dict_ctx->dict.len);
        dict_ctx->word_list.len += header ? BGEN32 (header->num_snips) : 0;
        dict_ctx->dict.len      += (uint64_t)vb->fragment_len;

        ASSERT (dict_ctx->dict.len <= dict_ctx->dict.size, "Dict %s len=%"PRIu64" exceeds allocated size=%"PRIu64" (num_fragments=%u fragment_1_len=%u)", 
                dict_ctx->tag_name, dict_ctx->dict.len, (uint64_t)dict_ctx->dict.size, dict_ctx->dict.prm32[0], dict_ctx->dict.prm32[1]);

        if (flag.debug_read_ctxs)
            sections_show_header ((SectionHeaderP)header, NULL, dict_sec->offset, sections_read_prefix);
    }

done: 
    // note: in cases we just "goto" here, no data is read, and a thread is needlessly created to decompress it
    // this is because the vb_i of the section needs to match the vb_i of the thread
    vb->dispatch = READY_TO_COMPUTE;
}

// entry point of compute thread of dictionary decompression
static void dict_io_uncompress_one_vb (VBlockP vb)
{
    if (!vb->fragment_ctx || flag.only_headers) goto done; // nothing to do in this thread
    SectionHeaderDictionary *header = (SectionHeaderDictionary *)vb->z_data.data;

    ASSERT (vb->fragment_start + BGEN32 (header->data_uncompressed_len) <= BAFTc (vb->fragment_ctx->dict), 
            "Buffer overflow when uncompressing dict=%s", vb->fragment_ctx->tag_name);

    // a hack for uncompressing to a location within the buffer - while multiple threads are uncompressing into 
    // non-overlappying regions in the same buffer in parallel
    Buffer copy = vb->fragment_ctx->dict;
    copy.data   = vb->fragment_start;
    zfile_uncompress_section (vb, header, &copy, NULL, 0, SEC_DICT); // NULL name prevents buf_alloc

done:
    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

static void dict_io_dict_build_word_list_one (ContextP zctx)
{
    if (!zctx->word_list.len || zctx->word_list.data) return; // skip if 1. no words, or 2. already built

    buf_alloc (evb, &zctx->word_list, 0, zctx->word_list.len, CtxWord, 0, "contexts->word_list");
    buf_set_overlayable (&zctx->word_list);

    rom word_start = zctx->dict.data;
    for (uint64_t snip_i=0; snip_i < zctx->word_list.len; snip_i++) {

        rom c=word_start; while (*c) c++;

        uint64_t index = BNUM64 (zctx->dict, word_start); 
        uint64_t len   = c - word_start;

        if (index > ZWORD_MAX_INDEX || len > ZWORD_MAX_LEN) {
            ASSERT (!VER(14), "A word was found in zctx=%s with index=%"PRIu64" amd len=%"PRIu64". This index/len is beyond current limits of Genozip. Use Genozip v13 to decompress this file.",
                    zctx->tag_name, (uint64_t)index, (uint64_t)len);

            ABORT ("A word was found in zctx=%s with index=%"PRIu64" amd len=%"PRIu64", which are beyond the limits",
                    zctx->tag_name, (uint64_t)index, (uint64_t)len);
        }       

        *B(CtxWord, zctx->word_list, snip_i) = (CtxWord){ .index = index, .len = len };

        word_start = c+1; // skip over the \0 separator
    }

    ASSERT (word_start == BAFTc(zctx->dict), "When building %s.word_list: expected to consume dict.len=%"PRIu64" bytes (in %"PRIu64" words), but consumed %"PRIu64" bytes",
            zctx->tag_name, zctx->dict.len, zctx->word_list.len, BNUM64(zctx->dict, word_start));
}

static void dict_io_build_word_lists (void)
{    
    START_TIMER;

    Context *last = ZCTX((flag.show_headers && is_genocat) ? CHROM : z_file->num_contexts);

    for (Context *zctx=z_file->contexts; zctx <= last; zctx++) 
        dict_io_dict_build_word_list_one (zctx);

    COPY_TIMER_VB (evb, dict_io_build_word_lists);
}

// PIZ main thread
void dict_io_read_all_dictionaries (void)
{
    START_TIMER;

    ctx_initialize_predefined_ctxs (z_file->contexts, z_file->data_type, z_file->d2d_map, &z_file->num_contexts);

    dict_sec = NULL;
    dict_ctx = NULL;

    dispatcher_fan_out_task (flag.reading_reference ? "read_dicts_ref" 
                            :flag.reading_chain     ? "read_dicts_chain" 
                            :flag.reading_kraken    ? "read_dicts_kraken" 
                            :                         "read_dicts",
                             NULL, 0, 0, 
                             true, flag.test,
                             !VER(9), // For v8 files, we read all fragments in the main thread as was the case in v8. This is because they are very small, and also we can't easily calculate the totel size of each dictionary.
                             0, 
                             10, // must be short - many dictionaries are just a few bytes
                             dict_io_read_one_vb, 
                             dict_io_uncompress_one_vb,
                             NO_CALLBACK);

    // build word lists in z_file->contexts with dictionary data 
    if (!flag.only_headers)
        dict_io_build_word_lists();

    // output the dictionaries if we're asked to
    if (flag.show_dict || flag.show_one_dict || flag.list_chroms) {
        for (uint32_t did_i=0; did_i < z_file->num_contexts; did_i++) {
            Context *ctx = ZCTX(did_i);
            if (!ctx->dict.len) continue;

            if (flag.list_chroms && ((!flag.luft && ctx->did_i == CHROM) || (flag.luft && ctx->did_i == VCF_oCHROM)))
                str_print_dict (info_stream, STRb(ctx->dict), true, Z_DT(SAM));
            
            if (flag.show_dict) 
                iprintf ("%s (did_i=%u, num_snips=%u, dict_size=%"PRIu64" bytes)\n", 
                         ctx->tag_name, did_i, ctx->word_list.len32, ctx->dict.len);

            if (dict_id_is_show (ctx->dict_id))
                str_print_dict (info_stream, STRb(ctx->dict), true, false);
        }
        iprint0 ("\n");

        if (is_genocat) exit_ok(); // if this is genocat - we're done
    }

    COPY_TIMER_VB (evb, dict_io_read_all_dictionaries);
}
