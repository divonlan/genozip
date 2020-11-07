// ------------------------------------------------------------------
//   codec.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "codec.h"
#include "vblock.h"
#include "strings.h"
#include "dict_id.h"
#include "file.h"
#include "zfile.h"

// --------------------------------------
// memory functions that serve the codecs
// --------------------------------------

// memory management for bzlib - tesing shows that compress allocates 4 times, and decompress 2 times. Allocations are the same set of sizes
// every call to compress/decompress with the same parameters, independent on the contents or size of the compressed/decompressed data.
void *codec_alloc (VBlock *vb, int size, double grow_at_least_factor)
{
    // get the next buffer - allocations are always in the same order in bzlib and lzma -
    // so subsequent VBs will allocate roughly the same amount of memory for each buffer
    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        if (!buf_is_allocated (&vb->codec_bufs[i])) {
            buf_alloc (vb, &vb->codec_bufs[i], size, grow_at_least_factor, "codec_bufs", i);
            //printf ("codec_alloc: %u bytes buf=%u\n", size, i);
            return vb->codec_bufs[i].data;
        }

    ABORT ("Error: codec_alloc could not find a free buffer. vb_i=%d", vb->vblock_i);
    return 0; // squash compiler warning
}

void codec_free (VBlock *vb, void *addr)
{
    if (!addr) return; // already freed

    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        if (vb->codec_bufs[i].data == addr) {
            buf_free (&vb->codec_bufs[i]);
            //printf ("codec_free: buf=%u\n", i);
            return;
        }

    char addr_str[POINTER_STR_LEN];
    ABORT ("Error: codec_free failed to find buffer to free. vb_i=%d addr=%s", 
           vb->vblock_i, str_pointer (addr, addr_str));
}

void codec_free_all (VBlock *vb)
{
    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        buf_free (&vb->codec_bufs[i]);
}

static bool codec_compress_error (VBlock *vb, SectionHeader *header, const char *uncompressed, uint32_t *uncompressed_len, LocalGetLineCB callback,
                                  char *compressed, uint32_t *compressed_len, bool soft_fail) 
{
    ABORT ("Error in comp_compress: Unsupported codec: %s", codec_name (header->codec));
    return false;
}


static void codec_uncompress_error (VBlock *vb, Codec codec, 
                                    const char *compressed, uint32_t compressed_len,
                                    Buffer *uncompressed_buf, uint64_t uncompressed_len,
                                    Codec sub_codec)
{
    ABORT ("Error in comp_uncompress: Unsupported codec: %s", codec_name (codec));
}

static void codec_reconstruct_error (VBlockP vb, Codec codec, ContextP ctx)
{
    ABORT ("Error in piz_reconstruct_from_ctx_do: in ctx=%s - codec %s has no LT_CODEC reconstruction", 
           err_dict_id (ctx->dict_id), codec_name (codec));
}

static uint32_t codec_est_size_default (Codec codec, uint64_t uncompressed_len)
{
    return (uint32_t)MAX (uncompressed_len / 2, 500);
}

// returns 4-character codec name
const char *codec_name (Codec codec)
{
    switch (codec) {
        case 0 ... NUM_CODECS : return codec_args[codec].name;
        default               : return "BAD!";    
    }
}

void codec_initialize (void)
{
    codec_bsc_initialize();
}

// ------------------------------
// Automatic codec selection
// ------------------------------

typedef struct {
    Codec codec;
    double size;
    double clock;
} CodecTest;

static int codec_assign_sorter (const CodecTest *t1, const CodecTest *t2)
{
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

    // time and size are very similar (within %1 and 15% respectively) - select for smaller size
    return t1->size - t2->size;
}

// this function tests each of our generic codecs on a 100KB sample of local or b250 data, and assigns the best one based on 
// compression ratio, or if the ratio is very similar, and the time is quite different, then based on time.
// the codec is then committed to zf_ctx, so that future VBs that clone recieve it and needn't test again.
// This function is called from two places:
// 1. For contexts with generic codecs, left as CODEC_UNKNOWN by the segmenter, we are called from zip_assign_best_codec.
//    For vb=1, this is called while holding the vb=1 lock, so that for all such contexts that appear in vb=1, they are
//    guaranteed to be tested only once. For contexts that make a first appearance in a later VB, parallel VBs might test
//    in parallel. A bit wasteful, but no harm.
// 2. For "specific" codecs (DOMQUAL, HAPMAT, GTSHARK...), subordinate contexts generated during compression of the primary
//    context (compression runs after zip_assign_best_codec is completed already) - those codecs explicitly call us to get the
//    codec for the subordinate context. Multiple of the early VBs may call in parallel, but future VBs will receive
//    the codec during cloning    
void codec_assign_best_codec (VBlockP vb, ContextP ctx, bool is_local, uint32_t len)
{
    RESET_FLAG (show_headers);
    uint64_t save_section_list = vb->section_list_buf.len; // save section list as comp_compress adds to it
    uint64_t save_z_data       = vb->z_data.len;

    CodecTest tests[] = { { CODEC_BZ2 }, { CODEC_NONE }, { CODEC_BSC }, { CODEC_LZMA } };
    const unsigned num_tests = flag.fast ? 2 : 4; // don't consider BSC or LZMA if --fast as they are slow

    Codec *selected_codec = is_local ? &ctx->lcodec : &ctx->bcodec;

    len = MIN (len, CODEC_ASSIGN_SAMPLE_SIZE);
    if (len < MIN_LEN_FOR_COMPRESSION ||  // if too small - don't assign - compression will use the default BZ2 and the next VB can try to select
        *selected_codec != CODEC_UNKNOWN) goto done; // if already selected - don't assign

    // last attempt to avoid double checking of the same context by parallel threads (as we're not locking, 
    // it doesn't prevent double testing 100% of time, but that's good enough) 
    Codec zf_codec = is_local ? z_file->contexts[ctx->did_i].lcodec :  // read without locking (1 byte)
                                z_file->contexts[ctx->did_i].bcodec;
    
    if (zf_codec != CODEC_UNKNOWN) {
        *selected_codec = zf_codec;
        goto done;
    }

    // measure the compressed size and duration for a small sample of of the local data, for each codec
    for (unsigned t=0; t < num_tests; t++) {
        *selected_codec = tests[t].codec;

        clock_t start_time = clock();
        tests[t].size  = (*selected_codec == CODEC_NONE) ? len : 
                                                           is_local ? zfile_compress_local_data (vb, ctx, len)
                                                                    : zfile_compress_b250_data (vb, ctx);
        tests[t].clock = (clock() - start_time);
    }

    // sort codec by our selection criteria
    qsort (tests, num_tests, sizeof (CodecTest), (int (*)(const void *, const void*))codec_assign_sorter);

    if (flag.show_codec_test)
        fprintf (stderr, "vb_i=%u %-8s %-5s [%-4s %5d %4.1f] [%-4s %5d %4.1f] [%-4s %5d %4.1f] [%-4s %5d %4.1f]\n", 
                vb->vblock_i, ctx->name, is_local ? "LOCAL" : "B250",
                codec_name (tests[0].codec), (int)tests[0].size, tests[0].clock,
                codec_name (tests[1].codec), (int)tests[1].size, tests[1].clock,
                codec_name (tests[2].codec), (int)tests[2].size, tests[2].clock,
                codec_name (tests[3].codec), (int)tests[3].size, tests[3].clock);

    // assign the best codec - the first one in the sorted array - and commit it to zf_ctx
    *selected_codec = tests[0].codec;
    ctx_commit_codec_to_zf_ctx (vb, ctx, is_local);

done:
    // roll back
    vb->z_data.len = save_z_data;
    vb->section_list_buf.len = save_section_list; 
    RESTORE_FLAG (show_headers);
}

// needs to be after all the functions as it refers to them
CodecArgs codec_args[NUM_CODECS] = CODEC_ARGS;

