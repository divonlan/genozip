// ------------------------------------------------------------------
//   codec.c
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "codec.h"
#include "vblock.h"
#include "strings.h"
#include "dict_id.h"
#include "file.h"
#include "zfile.h"
#include "profiler.h"
#include "bgzf.h"

// --------------------------------------
// memory functions that serve the codecs
// --------------------------------------

// memory management for bzlib - tesing shows that compress allocates 4 times, and decompress 2 times. Allocations are the same set of sizes
// every call to compress/decompress with the same parameters, independent on the contents or size of the compressed/decompressed data.
void *codec_alloc_do (VBlockP vb, uint64_t size, float grow_at_least_factor, rom func, uint32_t code_line)
{
    rom names[NUM_CODEC_BUFS] = { "codec_bufs[0]", "codec_bufs[1]", "codec_bufs[2]", "codec_bufs[3]",
                                          "codec_bufs[4]", "codec_bufs[5]", "codec_bufs[6]" };

    // get the next buffer - allocations are always in the same order in bzlib and lzma -
    // so subsequent VBs will allocate roughly the same amount of memory for each buffer
    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        if (!buf_is_alloc (&vb->codec_bufs[i])) {
            buf_alloc_do (vb, &vb->codec_bufs[i], size, grow_at_least_factor, func, code_line, names[i]);
            //printf ("codec_alloc: %u bytes buf=%u\n", size, i);
            return vb->codec_bufs[i].data;
        }
    ABORT_R ("%s: called from %s:%u codec_alloc could not find a free buffer", VB_NAME, func, code_line);
}

void codec_free_do (void *vb_, void *addr, rom func, uint32_t code_line)
{
    VBlockP vb = (VBlockP)vb_;

    if (!addr) return; // already freed

    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        if (vb->codec_bufs[i].data == addr) {
            buf_free_do (&vb->codec_bufs[i], func, code_line);
            return;
        }

    ABORT ("Error: codec_free_do called from %s:%u: failed to find buffer to free. vb_i=%d codec=%s addr=%p", 
           func, code_line, vb->vblock_i, codec_name(vb->codec_using_codec_bufs), addr);
}

void codec_free_all (VBlockP vb)
{
    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        buf_free (vb->codec_bufs[i]);
}

void codec_verify_free_all (VBlockP vb, rom op, Codec codec)
{
    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        ASSERT (!buf_is_alloc (&vb->codec_bufs[i]), "About to call %s for %s, but codec_bufs[%u] is allocated, expecting it to be free: %s",
                op, codec_name (codec), i, buf_desc (&vb->codec_bufs[i]).s);
}

static bool codec_compress_error (VBlockP vb, SectionHeader *header, rom uncompressed, uint32_t *uncompressed_len, LocalGetLineCB callback,
                                  char *compressed, uint32_t *compressed_len, bool soft_fail, rom name) 
{
    ABORT_R ("compression of \"%s\": Unsupported codec: %s", name, codec_name (header->codec));
}


static void codec_uncompress_error (VBlockP vb, Codec codec, uint8_t param,
                                    rom compressed, uint32_t compressed_len,
                                    BufferP uncompressed_buf, uint64_t uncompressed_len,
                                    Codec sub_codec, rom name)
{
    ABORT ("Error in comp_uncompress: \"%s\": unsupported codec: %s. Please upgrade to the most recent version of Genozip.", name, codec_name (codec));
}

static void codec_reconstruct_error (VBlockP vb, Codec codec, ContextP ctx)
{
    ABORT ("Error in reconstruct_from_ctx_do: in ctx=%s - codec %s has no LT_CODEC reconstruction. Please upgrade to the most recent version of Genozip.", 
           dis_dict_id (ctx->dict_id).s, codec_name (codec));
}

static uint32_t codec_est_size_default (Codec codec, uint64_t uncompressed_len)
{
    return (uint32_t)MAX_(uncompressed_len / 2, 500);
}

// returns 4-character codec name
rom codec_name (Codec codec)
{
    return (codec >=0 && codec < NUM_CODECS) ? codec_args[codec].name : "BAD!";    
}

void codec_initialize (void)
{
    codec_bsc_initialize();
    bgzf_libdeflate_initialize();
}

// ------------------------------
// Automatic codec assignment
// ------------------------------

typedef struct {
    Codec codec;
    float size;
    float clock;
} CodecTest;

static SORTER (codec_assign_sorter)
{
    CodecTest *t1 = (CodecTest *)a;
    CodecTest *t2 = (CodecTest *)b;
    
    // in --best mode, we take the smallest size, regardless of speed
    if (flag.best) {
        if (t1->size != t2->size) return ASCENDING (CodecTest, size);
        else                      return ASCENDING (CodecTest, clock);
    }

    // in --fast mode - if one if significantly faster with a modest size hit, take it. Otherwise, take the best.
    if (flag.fast) {
        if (t1->clock < t2->clock * 0.80 && t1->size < t2->size * 1.3) return -1; // t1 has 20% or more better time with at most 30% size hit
        if (t2->clock < t1->clock * 0.80 && t2->size < t1->size * 1.3) return  1; 
    }

    // case: select for significant difference in size (more than 2%)
    if (t1->size  < t2->size  * 0.98) return -1; // t1 has significantly better size
    if (t2->size  < t1->size  * 0.98) return  1; // t2 has significantly better size

    // case: size is similar, select for significant difference in time (more than 50%)
    if (t1->clock < t2->clock * 0.50) return -1; // t1 has significantly better time
    if (t2->clock < t1->clock * 0.50) return  1; // t2 has significantly better time

    // case: size and time are quite similar, check 2nd level 

    // case: select for smaller difference in size (more than 1%)
    if (t1->size  < t2->size  * 0.99) return -1; // t1 has significantly better size
    if (t2->size  < t1->size  * 0.99) return  1; // t2 has significantly better size

    // case: select for smaller difference in time (more than 15%)
    if (t1->clock < t2->clock * 0.85) return -1; // t1 has significantly better time
    if (t2->clock < t1->clock * 0.85) return  1; // t2 has significantly better time

    // if size is exactly the same - choose the lower-index codex (so to pick non-packing version of RAN and ART)
    if (t1->size == t2->size)
        return ASCENDING(CodecTest, codec);

    // time and size are very similar (within %1 and 15% respectively) - select for smaller size
    return ASCENDING(CodecTest, size);
}

// this function tests each of our generic codecs on a 100KB sample of local or b250 data, and assigns the best one based on 
// compression ratio, or if the ratio is very similar, and the time is quite different, then based on time.
// the codec is then committed to zctx, so that future VBs that clone recieve it and needn't test again.
// This function is called from two places:
// 1. For contexts with generic codecs, left as CODEC_UNKNOWN by the segmenter, we are called from zip_assign_best_codec.
//    For vb=1, this is called while holding the vb=1 lock, so that for all such contexts that appear in vb=1, they are
//    guaranteed to be tested only once. For contexts that make a first appearance in a later VB, parallel VBs might test
//    in parallel. A bit wasteful, but no harm.
// 2. For "specific" codecs (DOMQUAL, LONGR...), subordinate contexts generated during compression of the primary
//    context (compression runs after zip_assign_best_codec is completed already) - those codecs explicitly call us to get the
//    codec for the subordinate context. Multiple of the early VBs may call in parallel, but future VBs will receive
//    the codec during cloning    
//
// Codecs for contexts may be assigned in 3 stages:
// 1. During Seg (a must for all complex codecs - eg HT, DOMQ, ACGT...)
// 2. At merge - inherit from z_file->context if not set in Seg
// 3. After merge before compress - if still not assigned - zip_assign_best_codec - which also commits back to z_file->context
//    (this is the only place we commit to z_file, therefore z_file will only contain simple codecs)
// Note: if vb=1 commits a codec, it will be during its lock, so that all subsequent VBs will inherit it. But for
// contexts not committed by vb=1 - multiple contexts running in parallel may commit their codecs overriding each other. that's ok.
// 
// Note: in --best mode, we lock-in a codec only after at least (BEST_LOCK_IN_THREASHOLD) VBs in a row selected it (counting restarts if a VB
// selected a different codec)

#define BEST_LOCK_IN_THREASHOLD 30 // this number has a significant impact on compression and time, it was chosen after trial & error for --best

Codec codec_assign_best_codec (VBlockP vb, 
                               ContextP ctx, /* for b250, local, dict */ 
                               BufferP data, /* non-context data, or context data that doesn't reside in its normal buffer */
                               SectionType st)
{
    START_TIMER;

    bool is_local = (st == SEC_LOCAL);
    bool is_b250  = (st == SEC_B250 );
    Codec non_ctx_codec = CODEC_UNKNOWN; // used for non-b250, non-local sections
    bool data_override = ctx && data;    // caller provided data for us to use *instead* of ctx's own local/b250/dict buffer
    
    ContextP zctx = ctx ? ctx_get_zctx_from_vctx (ctx) : NULL;

    Codec *selected_codec = is_local ? &ctx->lcodec  : 
                            is_b250  ? &ctx->bcodec  :
                                       &non_ctx_codec;

    Codec zselected_codec = !zctx    ? CODEC_UNKNOWN :
                            is_local ? zctx->lcodec  : 
                            is_b250  ? zctx->bcodec  :
                                       CODEC_UNKNOWN ;

    uint8_t zselected_codec_count = !zctx ? 0 : is_local ? zctx->lcodec_count : is_b250 ? zctx->bcodec_count : 0;

    // we don't change assigned non-simple codecs
    if (!codec_args[*selected_codec].is_simple) 
        return *selected_codec;

    // in --best mode, we accept a codec that's been selected by 5 previous VBs in a row
    else if (flag.best && zselected_codec != CODEC_UNKNOWN && 
             zselected_codec_count >= BEST_LOCK_IN_THREASHOLD)
        *selected_codec = zselected_codec;

    // in --best mode, we invalidate (=re-test) a codec previously seleceted if its not yet been selected by (BEST_LOCK_IN_THREASHOLD) VBs in row
    else if (flag.best && (*selected_codec != zselected_codec || zselected_codec_count < BEST_LOCK_IN_THREASHOLD))
        *selected_codec = CODEC_UNKNOWN;

    // in non-best mode, we inherit from z_file if its set by a previous VB
    else if (!flag.best && zselected_codec != CODEC_UNKNOWN) 
        *selected_codec = zselected_codec;

    if (*selected_codec != CODEC_UNKNOWN) 
        return *selected_codec; // if already assigned - no need to test

    CodecTest tests[] = { { CODEC_BZ2 }, { CODEC_NONE }, { CODEC_BSC }, { CODEC_LZMA }, 
                          { CODEC_RANS8 }, { CODEC_RANS32 }, { CODEC_RANS8_pack }, { CODEC_RANS32_pack }, 
                          { CODEC_ARITH8 }, { CODEC_ARITH32 }, { CODEC_ARITH8_pack }, { CODEC_ARITH32_pack } }; // not using for now due to bug in arith_dynamic
    const unsigned num_tests = 12;
    
    // set data
    if (!data)
        switch (st) {
            case SEC_DICT  : data = &ctx->dict  ; break;
            case SEC_B250  : data = &ctx->b250  ; break; 
            case SEC_LOCAL : data = &ctx->local ; break;
            default: ASSERT (data, "expecting non-NULL data for section=%s", st_name (st));
        }

    uint64_t save_data_len = data->len;
    uint64_t save_section_list = vb->section_list_buf.len; // save section list as comp_compress adds to it

    data->len = MIN_(data->len * (is_local ? lt_desc[ctx->ltype].width : 1), CODEC_ASSIGN_SAMPLE_SIZE);

    if (data->len < MIN_LEN_FOR_COMPRESSION) 
        goto done; // if too small - don't assign - compression will use the default BZ2 and the next VB can try to select

    vb->z_data_test.param = true; // testing

    // measure the compressed size and duration for a small sample of of the local data, for each codec
    for (unsigned t=0; t < num_tests; t++) {
        *selected_codec = tests[t].codec;

        if (flag.show_time) codec_show_time (vb, "Assign", ctx->tag_name, *selected_codec);

        clock_t start_time = clock();

        if (*selected_codec == CODEC_NONE) tests[t].size = data->len;
        else {
            LocalGetLineCB *callback = !data_override ? zfile_get_local_data_callback (vb->data_type, ctx) : NULL;

            zfile_compress_section_data_ex (vb, SEC_NONE, callback ? NULL : data, callback, data->len, *selected_codec, SECTION_FLAGS_NONE);
            tests[t].size = vb->z_data_test.len;
            vb->z_data_test.len = 0;
        }
                                                           
        tests[t].clock = (clock() - start_time);
    }
    
    // sort codec by our selection criteria
    qsort (tests, num_tests, sizeof (CodecTest), codec_assign_sorter);

    if (flag.show_codec) {
        iprintf ("vb_i=%-2u %-12s %-5s %6.1fX  *[%-4s %5d %4.1f] [%-4s %5d %4.1f] [%-4s %5d %4.1f] [%-4s %5d %4.1f]\n", 
                 vb->vblock_i, ctx ? ctx->tag_name : &st_name (st)[4], ctx ? &st_name (st)[4] : "SECT",
                 (float)data->len / tests[0].size,
                 codec_name (tests[0].codec), (int)tests[0].size, tests[0].clock,
                 codec_name (tests[1].codec), (int)tests[1].size, tests[1].clock,
                 codec_name (tests[2].codec), (int)tests[2].size, tests[2].clock,
                 codec_name (tests[3].codec), (int)tests[3].size, tests[3].clock);
        fflush (info_stream);
    }

    // if we already have an assigned codec (i.e. this is --best) and this codec came in 2nd but has the
    // compression as the first, then keep it - so we don't constantly oscillate between two very similar codecs (eg. ARTb and ARTB)
    if (zselected_codec != CODEC_UNKNOWN && zselected_codec == tests[1].codec && tests[0].size == tests[1].size)
        *selected_codec = zselected_codec;

    // assign the best codec - the first one in the sorted array - and commit it to zctx
    else
        *selected_codec = tests[0].codec;

    // save the assignment for future VBs, but not in --best, where each VB tests on its own.
    // note: for local (except in --fast), we don't commit for vb=1 bc less representative of data 
    // (ok for --best as we count (BEST_LOCK_IN_THREASHOLD) anyway))
    if ((is_b250 || (is_local && (flag.best || flag.fast || vb->vblock_i > 1))) && *selected_codec != CODEC_UNKNOWN && zctx) 
        ctx_commit_codec_to_zf_ctx (vb, ctx, is_local);

done:
    // roll back
    data->len                = save_data_len;
    vb->section_list_buf.len = save_section_list; 
    vb->z_data_test.param    = false; // no longer testing

    COPY_TIMER (codec_assign_best_codec);

    return *selected_codec;
}

void codec_assign_best_qual_codec (VBlockP vb, DidIType qual_did_i,  
                                   LocalGetLineCB callback, bool no_longr)
{
    ContextP ctx = CTX(qual_did_i);

    if (segconf.running) {
        if (!flag.fast && segconf.nontrivial_qual && segconf_is_long_reads()) // note: we can't use segconf.is_long_reads (variable) because it is not set yet 
            codec_longr_segconf_calculate_bins (vb, CTX(qual_did_i+1), callback);
        return;
    }

    if (!flag.fast && segconf.is_long_reads && !no_longr && segconf.nontrivial_qual)
        codec_longr_comp_init (vb, qual_did_i);

    else if (!flag.no_domqual && codec_domq_comp_init (vb, qual_did_i, callback)) 
        {} // domqual    
    
    else if (Z_DT(DT_SAM) || Z_DT (DT_BAM)) { 
        ctx->ltype  = LT_CODEC;
        ctx->lcodec = CODEC_NORMQ;
    }

    else { // FASTQ
        ctx->ltype = LT_SEQUENCE;  // codec to be assigned by codec_assign_best_codec
        ctx->flags.is_qual = true;
    }
    
    if (flag.show_codec) // aligned to the output of codec_assign_best_codec
        iprintf ("vb_i=%-2u %-12s %-5s          *[%s]\n", vb->vblock_i, "QUAL", "LOCAL", codec_name(CTX(qual_did_i)->lcodec));
}

// complex codec est size - result may be recompressed with RAN, ART
uint32_t codec_complex_est_size (Codec codec, uint64_t uncompressed_len)
{
    uint64_t preprocessed_len = uncompressed_len * 1.1; 
    return MAX_(codec_ARTB_est_size (0, preprocessed_len), codec_ARTB_est_size (0, preprocessed_len)) + 10000;
}

// print prefix to the string to be printed by COPY_TIMER in the codec compression functions in case
// of eg. --show-time=compressor_lzma
void codec_show_time (VBlockP vb, rom name, rom subname, Codec codec)
{
    if ((strcmp (flag.show_time, "compressor_lzma"  ) && codec==CODEC_LZMA) ||
        (strcmp (flag.show_time, "compressor_bsc"   ) && codec==CODEC_BSC ) || 
        (strcmp (flag.show_time, "compressor_acgt"  ) && codec==CODEC_ACGT) || 
        (strcmp (flag.show_time, "compressor_domq"  ) && codec==CODEC_DOMQ) || 
        (strcmp (flag.show_time, "compressor_pbwt"  ) && codec==CODEC_PBWT) || 
        (strcmp (flag.show_time, "compressor_longr" ) && codec==CODEC_LONGR) || 
        (strcmp (flag.show_time, "compressor_rans"  ) && (codec==CODEC_RANS32 || codec==CODEC_RANS32_pack || codec==CODEC_RANS8 || codec==CODEC_RANS32_pack)) || 
        (strcmp (flag.show_time, "compressor_arith" ) && (codec==CODEC_ARITH32 || codec==CODEC_ARITH32_pack || codec==CODEC_ARITH8 || codec==CODEC_ARITH32_pack)) || 
        (strcmp (flag.show_time, "compressor_bz2"   ) && codec==CODEC_BZ2 )) {

        vb->profile.next_name    = name;
        vb->profile.next_subname = subname;
    }
}

UNCOMPRESS (codec_gtshark_uncompress)
{
    ABORT0 ("Support for the gtshark codec has been discontinued. To decompress VCF files compressed with --gtshark, use Genozip v11.");
}

// needs to be after all the functions as it refers to them
CodecArgs codec_args[NUM_CODECS] = CODEC_ARGS;

