// ------------------------------------------------------------------
//   codec.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "codec.h"
#include "vblock.h"
#include "strings.h"
#include "dict_id.h"
#include "file.h"
#include "zfile.h"
#include "zip.h"
#include "profiler.h"
#include "mgzip.h"
#include "sorter.h"

#define BEST_LOCK_IN_THREASHOLD 4 // --best: number of VBs that agree on the codec, after which we lock it in --best
#define DATA_LEN_HTS_OK (1 KB)    // --best: if data is smaller than this, switching between RANs and ARITH codecs is still considered consistent
#define RETEST_VB_I 10            // (not --fast/--best): retest in this vblock_i, in case early test was not representative

// --------------------------------------
// memory functions that serve the codecs
// --------------------------------------

// memory management for codecs - tesing shows that compress allocates 4 times, and decompress 2 times. Allocations are the same set of sizes
// every call to compress/decompress with the same parameters, independent on the contents or size of the compressed/decompressed data.
void *codec_alloc_do (VBlockP vb, uint64_t size, float grow_at_least_factor, unsigned *buf_i, FUNCLINE)
{
    rom names[NUM_CODEC_BUFS] = { "codec_bufs[0]", "codec_bufs[1]", "codec_bufs[2]", "codec_bufs[3]",
                                  "codec_bufs[4]", "codec_bufs[5]", "codec_bufs[6]" };

    // get the next buffer - allocations are always in the same order in bzlib and lzma -
    // so subsequent VBs will allocate roughly the same amount of memory for each buffer
    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        if (!buf_is_alloc (&vb->codec_bufs[i])) {
            vb->codec_bufs[i].can_be_big = true; // LZMA, for example, can allocate ~4GB buffers in --best
            buf_alloc_(vb, &vb->codec_bufs[i], 0, size, 1, grow_at_least_factor, names[i], func, code_line);
            // printf ("%s:%u codec_alloc: buf_i=%u %"PRIu64" bytes\n", func, code_line, i, size);
            if (buf_i) *buf_i = i;
            return vb->codec_bufs[i].data;
        }
    ABORT ("%s: called from %s:%u codec_alloc could not find a free buffer", VB_NAME, func, code_line);
}

void codec_free_do (void *vb_, void *addr, FUNCLINE)
{
    VBlockP vb = (VBlockP)vb_;

    if (!addr) return; // already freed

    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        if (vb->codec_bufs[i].data == addr) {
            // printf ("%s:%u codec_free:buf=%u\n", func, code_line, i);
            buf_free_do (&vb->codec_bufs[i], func, code_line);
            return;
        }

    ABORT ("Error: codec_free_do called from %s:%u: failed to find buffer to free. vb_i=%d codec=%s addr=%p", 
           func, code_line, vb->vblock_i, codec_name(vb->codec_using_codec_bufs), addr);
}

void codec_destroy_all (VBlockP vb)
{
    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        buf_destroy (vb->codec_bufs[i]);
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

static COMPRESS (codec_compress_error) 
{
    ABORT ("compression of \"%s\": Unsupported codec: %s", name, codec_name (header->codec));
}


static UNCOMPRESS (codec_uncompress_error)
{
    ABORT ("in comp_uncompress: \"%s\": unsupported codec: %s. %s", name, codec_name (codec), genozip_update_msg());
}

static CODEC_RECONSTRUCT (codec_reconstruct_error)
{
    ABORT ("in reconstruct_from_ctx_do: in ctx=%s - codec %s has no LT_CODEC reconstruction. %s", 
           dis_dict_id (ctx->dict_id).s, codec_name (codec), genozip_update_msg());
}

static uint32_t codec_est_size_default (Codec codec, uint64_t uncompressed_len)
{
    return (uint32_t)MAX_(uncompressed_len / 2, 500);
}

// returns 4-character codec name
rom codec_name (Codec codec)
{
    return IN_RANGE(codec, 0, NUM_CODECS) ? codec_args[codec].name : "BAD!";    
}

void codec_initialize (void)
{
    codec_bsc_initialize();
    bgzf_libdeflate_1_7_initialize();
}

// ------------------------------
// Automatic codec assignment
// ------------------------------

typedef struct {
    Codec codec;
    float size;
    float clock; // POSIX clock ticks (CLOCKS_PER_SEC=1000000)
} CodecTest;

static SORTER (codec_assign_sorter)
{
    CodecTest *t1 = (CodecTest *)a;
    CodecTest *t2 = (CodecTest *)b;
    
    // in --fast mode - if one is significantly faster with a modest size hit, take it. Otherwise, take the best.
    if (flag.fast) {
        if (t1->clock < t2->clock * 0.80 && t1->size < t2->size * 1.3) return -1; // t1 has 20% or more better time with at most 30% size hit
        if (t2->clock < t1->clock * 0.80 && t2->size < t1->size * 1.3) return  1; 
    }

    // in --best mode, or in normal mode if speed is fast enough in both cases, we take the smallest size, regardless of speed
    if (flag.best || (t1->clock <= 5000 && t2->clock <= 5000)) {
        if (t1->size != t2->size) return ASCENDING (CodecTest, size);
        else                      return ASCENDING (CodecTest, clock);
    }

    // if both are tiny, take the fastest
    if (t1->size < 100 && t2->size < 100 && t1->clock != t2->clock)
        return ASCENDING (CodecTest, clock);

    static struct { float size, time; } threasholds[] = {
    //    size   time
        { 0.96,  0.20 }, 
        { 0.97,  0.33 }, 
        { 0.98,  0.50 }, 
        { 0.985, 0.67 },
        { 0.99,  0.85 } 
    };

    for (int level=0; level < ARRAY_LEN(threasholds); level++) {
        if (t1->size  < t2->size  * threasholds[level].size) return -1; // t1 has significantly better size
        if (t2->size  < t1->size  * threasholds[level].size) return  1; // t2 has significantly better size

        // case: size is similar, select for significant difference in time 
        if (t1->clock < t2->clock * threasholds[level].time) return -1; // t1 has significantly better time
        if (t2->clock < t1->clock * threasholds[level].time) return  1; // t2 has significantly better time
    }

    // if size is exactly the same - choose the lower-index codec (so to pick non-packing version of RAN and ART)
    if (t1->size == t2->size)
        return ASCENDING(CodecTest, codec);

    // time and size are very similar (within %1 and 15% respectively) - select for smaller size
    return ASCENDING(CodecTest, size);
}

// this function tests each of our generic codecs on a 100KB sample of local or b250 data, and assigns the best one based on 
// compression ratio, or if the ratio is very similar, and the time is quite different, then based on time.
// the codec is then committed to zctx, so that future VBs that clone recieve it and needn't test again.
// This function is called from two places:
// 1. For contexts with generic codecs, left as CODEC_UNKNOWN by the segmenter, we are called from codec_assign_best_codec.
//    For vb=1, this is called while holding the vb=1 lock, so that for all such contexts that appear in vb=1, they are
//    guaranteed to be tested only once. For contexts that make a first appearance in a later VB, parallel VBs might test
//    in parallel. A bit wasteful, but no harm.
// 2. For "specific" codecs (DOMQUAL, LONGR...), subordinate contexts generated during compression of the primary
//    context (compression runs after codec_assign_best_codec is completed already) - those codecs explicitly call us to get the
//    codec for the subordinate context. Multiple of the early VBs may call in parallel, but future VBs will receive
//    the codec during cloning    
//
// Codecs for contexts may be assigned in 3 stages:
// 1. During Seg (a must for all complex codecs - eg HT, DOMQ, ACGT...)
// 2. At merge - inherit from z_file->context if not set in Seg
// 3. After merge before compress - if still not assigned - codec_assign_best_codec - which also commits back to z_file->context
//    (this is the only place we commit to z_file, therefore z_file will only contain simple codecs)
// Note: if vb=1 commits a codec, it will be during its lock, so that all subsequent VBs will inherit it. But for
// contexts not committed by vb=1 - multiple contexts running in parallel may commit their codecs overriding each other. that's ok.
// 
// Note: in --best mode, we lock-in a codec only after at least (BEST_LOCK_IN_THREASHOLD) VBs in a row selected it (counting restarts if a VB
// selected a different codec)

static inline bool vb_1_not_representative (VBlockP vb, ContextP ctx, bool is_local, bool is_b250)
{
    if (!is_local && !is_b250) return false;

    uint8_t bits = DTP(vb_1_not_representative);

    return ((bits & 0b001) && dict_id_is_field  (ctx->dict_id)) || 
           ((bits & 0b010) && dict_id_is_type_1 (ctx->dict_id)) || 
           ((bits & 0b100) && dict_id_is_type_2 (ctx->dict_id));
}

// update zctx with codec as it is assigned - don't wait for merge, to increase the chance that subsequent
// VBs can get this codec and don't need to test for themselves.
// note: must be called with zctx->assign_codec_mutex[] locked
static void codec_commit_to_zf_ctx (uint32_t data_len, Codec v_codec, Codec *z_codec, uint8_t *z_codec_count)
{
    // not best: just commit the z_codec without discretion
    if (!flag.best)
        *z_codec = v_codec; 

    // case: --best: codec remains (at least within codec family) - increment counter
    else if (v_codec == *z_codec || 
        (is_rans (v_codec) && is_rans (*z_codec)) || 
        (is_arith (v_codec) && is_arith (*z_codec)) || 
        (data_len < DATA_LEN_HTS_OK && is_hts (v_codec) && is_hts (*z_codec))) // if data size for this field in this VB is tiny, accept as consistent even if flipping between rans and arith
        
        *z_codec_count = MIN_(255, 1+(int)*z_codec_count);    
   
    // codec family changed - update z_codec and reset counter
    else { 
        *z_codec_count = 0; 
        *z_codec = v_codec; 
    }
}

Codec codec_assign_best_codec (VBlockP vb, 
                               ContextP ctx, /* for b250, local, dict */ 
                               BufferP data, /* non-context data, or context data that doesn't reside in its normal buffer */
                               SectionType st)
{
    START_TIMER;

    bool is_local = ST(LOCAL);
    bool is_b250  = ST(B250);
    Codec non_ctx_codec = CODEC_UNKNOWN; // used for non-b250, non-local sections
    bool data_override = ctx && data;    // caller provided data for us to use *instead* of ctx's own local/b250/dict buffer
    
    ContextP zctx = ctx ? ctx_get_zctx_from_vctx (ctx, true, false) : NULL; // note: if is_dict, ctx==zctx
    
    // if another thread is busy testing, it will be faster to wait for their results than re-test ourselves
    if (is_local || is_b250)
        mutex_lock (zctx->assign_codec_mutex[is_local]);

    // THREAD SAFETY: 
    // - zctx->{lcodec,bcodec,lcodec_count,bcodec_count}: access protected by zctx->assign_codec_mutex[]
    // - zctx->lcodec_hard_coded immutable after segconf
    Codec *restrict selected_codec = is_local?&ctx->lcodec       : is_b250?&ctx->bcodec       : &non_ctx_codec;
    Codec zselected_codec          = is_local?zctx->lcodec       : is_b250?zctx->bcodec       : CODEC_UNKNOWN ;
    uint8_t zselected_codec_count  = is_local?zctx->lcodec_count : is_b250?zctx->bcodec_count : 0;

    // case: codec is hard_coded and should not be reassigned
    bool hard_coded = is_local && ((ctx ? ctx->lcodec_hard_coded : false) || (zctx ? zctx->lcodec_hard_coded : false));

    // we don't change assigned non-simple codecs
    if (!codec_args[*selected_codec].is_simple) 
        goto return_selected;

    // if --best, we accept a codec that's been selected by 5 previous VBs in a row
    else if (flag.best && zselected_codec != CODEC_UNKNOWN && zselected_codec_count >= BEST_LOCK_IN_THREASHOLD)
        *selected_codec = zselected_codec;

    // if --best, we invalidate (=re-test) a codec previously seleceted if its not yet been selected by (BEST_LOCK_IN_THREASHOLD) VBs in row
    else if (flag.best && !hard_coded && zselected_codec_count < BEST_LOCK_IN_THREASHOLD)
        *selected_codec = CODEC_UNKNOWN;

    // in normal mode, reasign in the 10th VB for fields for which the initial assignment might be unrepresentative
    else if (vb->vblock_i == RETEST_VB_I && !flag.best && !flag.fast && !hard_coded &&  
             vb_1_not_representative (vb, ctx, is_local, is_b250))
        *selected_codec = CODEC_UNKNOWN;

    // if not reassigning, we inherit from z_file if its set by a previous VB
    else if (!flag.best && zselected_codec != CODEC_UNKNOWN) 
        *selected_codec = zselected_codec;

    if (*selected_codec != CODEC_UNKNOWN) 
        goto return_selected; // if already assigned - no need to test

    CodecTest tests[] = { { CODEC_NONE }, 
                          { CODEC_RANB }, { CODEC_RANW }, { CODEC_RANb }, { CODEC_RANw }, 
                          { CODEC_ARTB }, { CODEC_ARTW }, { CODEC_ARTb }, { CODEC_ARTw },
                          { CODEC_BZ2  }, { CODEC_BSC  }, { CODEC_LZMA } }; 
    
    // don't allow LZMA or BSC in buffers being compressed in the main thread - too slow (unless --best)
    unsigned num_tests = ARRAY_LEN(tests) - (flag.best ? 0 : 2*(vb == evb)); 

    if (flag.no_lzma && num_tests == ARRAY_LEN(tests))
        num_tests--;

    // set data
    if (!data)
        switch (st) {
            case SEC_DICT  : data = &ctx->dict  ; break;
            case SEC_B250  : data = &ctx->b250  ; break; 
            case SEC_LOCAL : data = &ctx->local ; break;
            default: ASSERT (data, "expecting non-NULL data for section=%s", st_name (st));
        }

    uint64_t save_data_len = data->len;
    uint64_t save_section_list = vb->section_list.len; // save section list as comp_compress adds to it

    data->len = MIN_(data->len * (is_local ? lt_width(ctx) : 1), CODEC_ASSIGN_SAMPLE_SIZE);

    if (data->len < MIN_LEN_FOR_COMPRESSION) 
        goto done; // if too small - don't assign - compression will use CODEC_NONE and the next VB can try to select

    vb->z_data_test.param = true; // testing

    // measure the compressed size and duration for a small sample of of the local data, for each codec
    for (unsigned t=0; t < num_tests; t++) {
        *selected_codec = tests[t].codec;

        if (flag.show_time) codec_show_time (vb, "Assign", ctx->tag_name, *selected_codec);

        clock_t start_time = clock();

        if (*selected_codec == CODEC_NONE) tests[t].size = data->len;
        else {
            LocalGetLineCB *callback = (ST(LOCAL) && !data_override && !ctx->no_callback) ? zip_get_local_data_callback (vb->data_type, ctx) : NULL;

            zfile_compress_section_data_ex (vb, ctx, SEC_RANDOM_ACCESS/*a section with SectionType header*/, callback ? NULL : data, callback, data->len, *selected_codec, SECTION_FLAGS_NONE, 
                                            "codec_assign_best_codec");
            tests[t].size = vb->z_data_test.len;
            vb->z_data_test.len = 0;
        }
                                                           
        tests[t].clock = (clock() - start_time) * (1000000 / CLOCKS_PER_SEC); // note: POSIX requires CLOCKS_PER_SEC=1000000. In Windows it is 1000.
    }

    // sort codec by our selection criteria
    qsort (tests, num_tests, sizeof (CodecTest), codec_assign_sorter);

    // --best: if we already have an assigned codec and this codec came in 2nd but has the
    // compression as the first, then keep it - so we don't constantly oscillate between two very similar codecs (eg. ARTb and ARTB)
    if (flag.best && zselected_codec != CODEC_UNKNOWN && zselected_codec == tests[1].codec && tests[0].size == tests[1].size)
        *selected_codec = zselected_codec;

    // assign the best codec - the first one in the sorted array - and commit it to zctx
    else
        *selected_codec = tests[0].codec;

    // save the assignment for future VBs, but not in --best, where each VB tests on its own.
    // note: for local (except in --fast), we don't commit for vb=1 if vb_1_not_representative  
    // (ok for --best as we count (BEST_LOCK_IN_THREASHOLD) anyway))
    if (vb->txt_data.len > MIN_(4 MB, segconf.vb_size / 2) && // don't let tiny VBs (e.g. SAM PRIM/DEPN or final VB of the file) set the codec for everyone
        *selected_codec != CODEC_UNKNOWN && zctx) {

        if (is_b250)
            codec_commit_to_zf_ctx (save_data_len, ctx->bcodec, &zctx->bcodec, &zctx->bcodec_count);

        else if (is_local && (flag.best || 
                              flag.fast || 
                              vb->vblock_i > vb_1_not_representative (vb, ctx, is_local, is_b250) || 
                              vb->is_last_vb_in_txt_file))
            codec_commit_to_zf_ctx (save_data_len, ctx->lcodec, &zctx->lcodec, &zctx->lcodec_count);
    }

    if (flag_is_set (show_codec, zctx ? zctx->dict_id : DICT_ID_NONE)) {
        iprintf ("%-8s %-12s %-5s %6.1fX   *[%-4s %5d B %6d μs]  [%-4s %5d B %6d μs]  [%-4s %5d B %6d μs]  [%-4s %5d B %6d μs]%s\n", 
                 VB_NAME, ctx ? ctx->tag_name : &st_name (st)[4], ctx ? &st_name (st)[4] : "SECT",
                 (float)data->len / tests[0].size,
                 codec_name (tests[0].codec), (int)tests[0].size, (int)tests[0].clock,
                 codec_name (tests[1].codec), (int)tests[1].size, (int)tests[1].clock,
                 codec_name (tests[2].codec), (int)tests[2].size, (int)tests[2].clock,
                 codec_name (tests[3].codec), (int)tests[3].size, (int)tests[3].clock,
                 ((is_local && flag.best && zctx->lcodec_count >= BEST_LOCK_IN_THREASHOLD) || (is_b250 && flag.best && zctx->bcodec_count >= BEST_LOCK_IN_THREASHOLD)) ? " --best-LOCK-IN" : "");
        fflush (info_stream);
    }

done:
    // restore
    data->len             = save_data_len;
    vb->section_list.len  = save_section_list; 
    vb->z_data_test.param = false; // no longer testing

return_selected:    
    if (is_local || is_b250)
        mutex_unlock (zctx->assign_codec_mutex[is_local]);

    COPY_TIMER (codec_assign_best_codec);
    return *selected_codec;
}

void codec_assign_best_qual_codec (VBlockP vb, Did did_i,  
                                   LocalGetLineCB callback, 
                                   bool no_seq_dependency, // only consider codecs in which SEQ is not required to compress QUAL
                                   bool maybe_revcomped,
                                   bool *codec_requires_seq)
{
    decl_ctx (did_i);
    decl_zctx (did_i); // shortcut to zctx, because did_i is predefined
    ASSERT (did_i < DTF(num_fields), "%s is not predefined", ctx->tag_name); // because of ZCTX()

    mutex_lock (zctx->assign_codec_mutex[true]); 

    // case: a previous VB already determined that the did_i doesn't need one of the complex codec
    if (zctx->qual_codec == CODEC_NONE) {
        ctx->ltype = LT_BLOB;  
        goto done;
    }

    Codec forced_codec = zctx->qual_codec                             ? zctx->qual_codec 
                       : (did_i == SAM_QUAL && flag.force_qual_codec) ? flag.force_qual_codec // note: SAM_QUAL==FASTQ_QUAL
                       :                                                0;

    if (forced_codec)  
        switch (forced_codec) {
            case CODEC_LONGR : codec_longr_comp_init (vb, did_i, true);           break;
            case CODEC_SMUX  : codec_smux_comp_init  (vb, did_i, callback, true); break;
            case CODEC_PACB  : codec_pacb_comp_init  (vb, did_i, callback, true); break;
            case CODEC_HOMP  : codec_homp_comp_init  (vb, did_i, callback, true); break; 
            case CODEC_TMPL  : codec_tmpl_comp_init  (vb, did_i, true);           break; 
            case CODEC_DOMQ  : codec_domq_comp_init  (vb, did_i, callback, true); break;
            case CODEC_NORMQ : codec_normq_comp_init (vb, did_i, maybe_revcomped, true); break;
            case CODEC_ARTB  : break;
            case CODEC_BSC   : break;
            default          : ABORT ("Can't force codec %s", codec_name (forced_codec));
        }

    else if (!no_seq_dependency && codec_pacb_comp_init (vb, did_i, callback, false));
    
    else if (!no_seq_dependency && codec_longr_comp_init (vb, did_i, false));

    else if (!no_seq_dependency && codec_homp_comp_init (vb, did_i, callback, false)); // only if Ultima, it might succeed. takes precedence of DOMQ

    else if (!no_seq_dependency && codec_smux_comp_init (vb, did_i, callback, false));
        
    else if (codec_tmpl_comp_init (vb, did_i, false));

    else if (codec_domq_comp_init (vb, did_i, callback, false));

    else if (codec_normq_comp_init (vb, did_i, maybe_revcomped, false)); 

    else
        ctx->ltype = LT_BLOB;  // codec to be assigned by codec_assign_best_codec
    
    zctx->qual_codec = (ctx->ltype == LT_BLOB ? CODEC_NONE : ctx->lcodec);

    #define LCODEC(x) (ctx->lcodec == CODEC_##x)
    if (codec_requires_seq && (LCODEC(PACB) || LCODEC(LONGR) || LCODEC(HOMP) || LCODEC(SMUX))) 
        *codec_requires_seq = true;

    if (!forced_codec && (flag_is_set (show_codec, (DictId)_SAM_QUAL) || flag.show_qual) && ctx->lcodec) // printing aligned to the output of codec_assign_best_codec
        iprintf ("%-8s %-12s %-5s           *[%s]\n", VB_NAME, ctx->tag_name, "LOCAL", codec_name(CTX(did_i)->lcodec));

done:
    mutex_unlock (zctx->assign_codec_mutex[true]); 
}

// complex codec est size - result may be recompressed with RAN, ART
uint32_t codec_complex_est_size (Codec codec, uint64_t uncompressed_len)
{
    uint64_t preprocessed_len = uncompressed_len * 1.1; 
    return MAX_(codec_ARTB_est_size (0, preprocessed_len), codec_ARTB_est_size (0, preprocessed_len)) + 10000;
}

uint32_t codec_trivial_size (Codec codec, uint64_t uncompressed_len)
{
    return 1; // QUAL.compressed_len is one byte
}

// print prefix to the string to be printed by COPY_TIMER in the codec compression functions in case
// of eg. --show-time=compressor_lzma
void codec_show_time (VBlockP vb, rom name, rom subname, Codec codec)
{
    if (!flag.show_time[0] || // --show-time with no argument
        (strcmp (flag.show_time, "compressor_lzma"  ) && codec==CODEC_LZMA)  ||
        (strcmp (flag.show_time, "compressor_bsc"   ) && codec==CODEC_BSC )  || 
        (strcmp (flag.show_time, "compressor_acgt"  ) && codec==CODEC_ACGT)  || 
        (strcmp (flag.show_time, "compressor_domq"  ) && codec==CODEC_DOMQ)  || 
        (strcmp (flag.show_time, "compressor_ulti"  ) && codec==CODEC_HOMP)  || 
        (strcmp (flag.show_time, "compressor_pbwt"  ) && codec==CODEC_PBWT)  || 
        (strcmp (flag.show_time, "compressor_longr" ) && codec==CODEC_LONGR) || 
        (strcmp (flag.show_time, "compressor_smux"  ) && codec==CODEC_SMUX)  || 
        (strcmp (flag.show_time, "compressor_pacb"  ) && codec==CODEC_PACB)  || 
        (strcmp (flag.show_time, "compressor_t0"    ) && codec==CODEC_T0)    || 
        (strcmp (flag.show_time, "compressor_rans"  ) && (codec==CODEC_RANW || codec==CODEC_RANw || codec==CODEC_RANB || codec==CODEC_RANw)) || 
        (strcmp (flag.show_time, "compressor_arith" ) && (codec==CODEC_ARTW || codec==CODEC_ARTw || codec==CODEC_ARTB || codec==CODEC_ARTw)) || 
        (strcmp (flag.show_time, "compressor_bz2"   ) && codec==CODEC_BZ2 )) {

        vb->profile.next_name    = name;
        vb->profile.next_subname = subname;
    }
}

void codec_qual_show_stats (void)
{
    if (!z_file || !z_file->num_lines) return;

    uint32_t domq  = z_file->domq_lines;
    uint32_t divr  = z_file->divr_lines;
    uint32_t homp  = z_file->homp_lines;
    uint32_t pacb  = z_file->pacb_lines;
    uint32_t longr = z_file->longr_lines;
    uint32_t normq = z_file->normq_lines;

    int32_t other = z_file->num_lines - domq - divr - homp - pacb - longr - normq; // allow negative to catch bugs

    iprintf ("\nQUAL codec stats (# lines): DOMQ=%u (%.1f%%) DIVR=%u (%.1f%%) HOMP=%u (%.1f%%) PACB=%u (%.1f%%) LONGR=%u (%.1f%%) NORMQ=%u (%.1f%%) other=%d (%.1f%%)\n",
             domq,  percent (domq,  z_file->num_lines), 
             divr,  percent (divr,  z_file->num_lines), 
             homp,  percent (homp,  z_file->num_lines), 
             pacb,  percent (pacb,  z_file->num_lines), 
             longr, percent (longr, z_file->num_lines), 
             normq, percent (normq, z_file->num_lines), 
             other, percent (other, z_file->num_lines));

    ASSERT (other >= 0, "expected other=%d >= 0", other);
}

UNCOMPRESS (codec_hapmat_uncompress)
{
    ABORT0 ("Support for the hapmat codec has been discontinued. To unscompress this VCF file, use Genozip v14.");
}

UNCOMPRESS (codec_gtshark_uncompress)
{
    ABORT0 ("Support for the gtshark codec has been discontinued. To uncompress VCF files compressed with --gtshark, use Genozip v11.");
}

// needs to be after all the functions as it refers to them
CodecArgs codec_args[NUM_CODECS] = CODEC_ARGS;

