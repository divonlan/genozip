// ------------------------------------------------------------------
//   codec_longr.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

// This codec is for quality scores of Nanopore and PacBio data, based on https://pubmed.ncbi.nlm.nih.gov/32470109/. This file contains 
// only original Genozip code, and the algorithm itself, derived from ENano source code, is located in codec_longr_alg.c

#include "genozip.h"
#include "codec.h"
#include "buffer.h"
#include "vblock.h"
#include "dict_id.h"
#include "reconstruct.h"
#include "strings.h"
#include "compressor.h"
#include "profiler.h"
#include "context.h"
#include "endianness.h"
#include "piz.h"
#include "strings.h"
#include "stats.h"
#include "sam_private.h"

#include "codec_longr_alg.c" // seperate source file for this, as it derived from external code with a different license

// similar structure to DOMQUAL 
void codec_longr_comp_init (VBlockP vb, DidIType qual_did_i)
{
    // lens_ctx contains the lens array - an array of uint32 - each entry is the length of the corresponding channel in values_ctx
    ContextP lens_ctx = CTX(qual_did_i);
    lens_ctx->ltype   = LT_CODEC;    // causes reconstruction to go to codec_longr_reconstruct
    lens_ctx->lcodec  = CODEC_LONGR;
    lens_ctx->counts_section = true; // we store the global value-to-bin mapper in a SEC_COUNTS

    // values_ctx contains the base quality data, sorted by channel
    ContextP values_ctx    = lens_ctx+1; // used
    values_ctx->ltype      = LT_UINT8;
    values_ctx->local_dep  = DEP_L1; 
    values_ctx->lcodec     = CODEC_ARITH8;
    
    stats_set_consolidation (vb, qual_did_i, 1, qual_did_i+1);
}

// upper bound of compressed size of length array
uint32_t codec_longr_est_size (Codec codec, uint64_t uncompressed_len) 
{ 
    return codec_complex_est_size (CODEC_LONGR, LONGR_NUM_CHANNELS * sizeof (uint32_t));
}

static inline void add_to_histogram (uint32_t *histogram, const uint8_t *values, uint32_t values_len)
{
    for (uint32_t i=0; i < values_len; i++) 
        histogram[values[i] - '!']++;
}

// ZIP, main thread, segconf.running
#define NUM_BINS (1 << QUAL_BITS)
void codec_longr_segconf_calculate_bins (VBlockP vb, ContextP ctx, 
                                         LocalGetLineCB callback) // option 2 - get one line
{
    ASSERT0 (segconf.running, "Expected segconf.running"); // must run in main thread as we are allocating from evb
    ContextP zctx = ZCTX(ctx->did_i);

    // create a histogram of values
    uint32_t histogram[256] = {};
    uint32_t num_values = 0;

    if (callback)
        for (uint64_t line_i=0; line_i < vb->lines.len; line_i++) {
            uint8_t *values; uint32_t values_len;
            callback (vb, line_i, (char **)pSTRa(values), vb->txt_data.len, NULL);
            add_to_histogram (histogram, STRa(values));

            num_values += values_len;
        }
    else {
        add_to_histogram (histogram, FIRSTENT (uint8_t, ctx->local), ctx->local.len);
        num_values = ctx->local.len;
    }

    // create value_to_bin mapper in zctx->con_cache (note: ctx but be predefinded, so zctx has the same did_i)
    buf_alloc (evb, &zctx->con_cache, 0, 256, uint8_t, 0, "value_to_bin");
    uint8_t *value_to_bin = FIRSTENT (uint8_t, zctx->con_cache);
    zctx->con_cache.len = NUM_BINS;
    
    uint8_t fixed = 11;

    // divide the values into as-equal-size-as-possible bins
    uint32_t next_val = 0;
    for (unsigned bin_i=0; bin_i < NUM_BINS; bin_i++) {

        uint32_t at_least = num_values / (NUM_BINS - bin_i); //  we need at least this number of values in this bin

        if (bin_i < fixed) {
            num_values -= histogram[next_val];
            value_to_bin[next_val] = bin_i;
            next_val++;

        }
        else {
            uint32_t bin_content = 0;
            while (bin_content < at_least && next_val < 256) {
                bin_content += histogram[next_val];
                num_values  -= histogram[next_val];
                value_to_bin[next_val] = bin_i;
                next_val++;
            }
        }
    }
    memset (&value_to_bin[next_val], NUM_BINS-1, 256 - next_val); // remaining high values (with 0 in their histogram) go to the top bin

    // store the value-to-bin map in the file as a global SEC_COUNTS section - convert to uint64 (luckily, it is small)
    buf_alloc (evb, &zctx->counts, 0, 256, uint64_t, 0, "contexts->ctx");
    for (unsigned i=0; i < 256; i++)
        NEXTENT (uint64_t, zctx->counts) = value_to_bin[i];
}

static void codec_longr_calc_channels (LongrState *state, STRp(seq), const uint8_t *qual, bool is_rev) 
{
    codec_longr_alg_init_read (state, STRa(seq), is_rev);

    uint8_t prev_q = 0;

    #define CALC_ONE(func_b) ({                                         \
        state->base_chan[state->next_base++] = state->chan.channel.n;   \
        state->chan_num_bases[state->chan.channel.n]++;                 \
        uint8_t b = (func_b);                                           \
        uint8_t q = (uint8_t)qual[i] - '!';                             \
        codec_longr_update_state (state, b, q, prev_q);                 \
        prev_q = q;                                                     \
    })

    if (!is_rev)
        for (uint32_t i=0; i < seq_len; i++) 
            CALC_ONE (acgt_encode[(uint8_t)codec_longr_next_base (STRa(seq), i)]);
    else
        for (int32_t i=seq_len-1; i >= 0; i--) 
            CALC_ONE (acgt_encode_comp[(uint8_t)codec_longr_next_base_rev (STRa(seq), i)]);
}

bool codec_longr_compress (VBlock *vb, 
                           SectionHeader *header,    
                           const char *uncompressed,     // option 1 - not supported
                           uint32_t *uncompressed_len, 
                           LocalGetLineCB qual_callback, // option 2 - get one line
                           char *compressed, uint32_t *compressed_len, // in/out 
                           bool soft_fail)               // soft fail not supported
{
    START_TIMER;
    ASSERTISNULL (uncompressed);

    ContextP lens_ctx     = ECTX (((SectionHeaderCtx *)header)->dict_id);
    ContextP values_ctx = lens_ctx + 1;
    
    LocalGetLineCB *seq_callback = (VB_DT(DT_FASTQ) ? fastq_zip_seq : sam_zip_seq);
    
    LongrState *state = codec_alloc (vb, sizeof (LongrState), 0);
    memset (state, 0, sizeof (LongrState));

    state->base_chan  = codec_alloc (vb, *uncompressed_len * sizeof (uint16_t), 0);
    state->value_to_bin = FIRSTENT (uint8_t, ZCTX(lens_ctx->did_i)->con_cache);
    codec_longr_alg_init (state);

    char *seq, *qual; 
    uint32_t seq_len, qual_len;
    bool is_rev;

    // calculate state->base_chan - the channel for each base. reads are treated in their
    // original (FASTQ) orientation.
    for (uint64_t line_i=0; line_i < vb->lines.len; line_i++) {

        seq_callback (vb, line_i,  pSTRa(seq), 0, &is_rev);
        qual_callback (vb, line_i, pSTRa(qual), vb->txt_data.len, NULL);
        
        codec_longr_calc_channels (state, STRa(seq), (uint8_t*)qual, is_rev);
    }
    
    // we now sort the quality data by channel (reversing qual of revcomp reads)
    buf_alloc (vb, &values_ctx->local, *uncompressed_len, 0, uint8_t, 0, "contexts->local");
    values_ctx->local.len = *uncompressed_len;

    uint8_t *sorted_qual = FIRSTENT (uint8_t, values_ctx->local);
    
    uint16_t *base_chan = state->base_chan;

    uint32_t next_of_chan[LONGR_NUM_CHANNELS];
    next_of_chan[0] = 0;
    for (uint32_t chan=1; chan < LONGR_NUM_CHANNELS; chan++) 
        next_of_chan[chan] = next_of_chan[chan-1] + state->chan_num_bases[chan-1];

    uint32_t base_i=0;
    for (uint64_t line_i=0; line_i < vb->lines.len; line_i++) {
        qual_callback (vb, line_i, pSTRa(qual), vb->txt_data.len, &is_rev);

        // we create a sorted_qual array, which contains LONGR_NUM_CHANNELS segments, one for each channel,
        // containing all the bases of that belong to that channel
        for (uint32_t i=0; i < qual_len; i++, base_i++) {
            uint16_t chan_i    = base_chan[base_i];
            uint32_t index     = next_of_chan[chan_i]++; // index into the segment belonging to chan_i in the sorted array
            sorted_qual[index] = qual[is_rev ? (qual_len-1-i) : i] - '!';
        }
    }

    // channel lengths 
    lens_ctx->local.len = 0; // overwrite previous QUAL->local.len
    buf_alloc (vb, &lens_ctx->local, LONGR_NUM_CHANNELS, 0, uint32_t, 0, "contexts->local");
    lens_ctx->local.len = LONGR_NUM_CHANNELS; 

    ARRAY (uint32_t, lens, lens_ctx->local);

    for (uint32_t chan=0; chan < LONGR_NUM_CHANNELS; chan++) 
        lens[chan] = BGEN32 (state->chan_num_bases[chan]); 

    codec_free_all (vb);

    // assigned codec to lens buffer (lens_ctx.local)
    lens_ctx->lcodec  = CODEC_UNKNOWN;
    lens_ctx->local.len *= sizeof (uint32_t);
    header->sub_codec = codec_assign_best_codec (vb, lens_ctx, &lens_ctx->local, SEC_LOCAL);
    if (header->sub_codec == CODEC_UNKNOWN) header->sub_codec = CODEC_NONE; // really small
    lens_ctx->lcodec  = CODEC_LONGR;

    CodecCompress *compress = codec_args[header->sub_codec].compress;
    *uncompressed_len = lens_ctx->local.len;

    // make sure we have enough memory - since this is a lengths array, our estimate is expected to always be sufficient, no need for soft_fail
    uint32_t min_required_compressed_len = codec_args[header->sub_codec].est_size (header->sub_codec, *uncompressed_len);
    if (*compressed_len < min_required_compressed_len) 
        ABORT ("Compressing %s in vb_i=%u with %s need %u bytes, but allocated only %u", lens_ctx->tag_name, vb->vblock_i, codec_name(header->sub_codec), min_required_compressed_len, *compressed_len);

    // length is a numeric context, and as such it must be written first, as piz_uncompress_all_ctxs needs to complete piz_adjust_one_local
    COPY_TIMER (compressor_longr);

    // actually compress the lens array
    return compress (vb, header, lens_ctx->local.data, uncompressed_len, NULL, compressed, compressed_len, false);
}

//--------------
// PIZ side
//--------------

static void codec_longr_recon_one_read (LongrState *state, STRp(seq), bool is_rev,
                                        const uint8_t *sorted_qual, uint32_t *next_of_chan, char *recon)
{
    codec_longr_alg_init_read (state, STRa(seq), is_rev);
    uint8_t prev_q = 0;

    #define RECON_ONE_QUAL                              \
        codec_longr_update_state (state, b, q, prev_q); \
        prev_q = q;                                     \
        recon[i] = q  + '!';

    if (!is_rev) // separate loops to save one "if" in the tight loop
        for (uint32_t i=0; i < seq_len; i++) {        
            uint8_t b = acgt_encode[(uint8_t)codec_longr_next_base (STRa(seq), i)];
            uint8_t q = sorted_qual[next_of_chan[state->chan.channel.n]++];
            RECON_ONE_QUAL;
        }
    else
        for (int32_t i=seq_len-1; i >= 0; i--) {        
            uint8_t b = acgt_encode_comp[(uint8_t)codec_longr_next_base_rev (STRa(seq), i)];
            uint8_t q = sorted_qual[next_of_chan[state->chan.channel.n]++];
            RECON_ONE_QUAL;
        }
}

// order of decompression: lens_ctx is decompressed, then baseq_ctx is decompressed with its codec, and then this function is called
// as a subcodec for baseq_ctx. 
// This function converts lens_ctx to be "next_of_chan" for each channel, and initializes LongrState
static void codec_longr_reconstruct_init (VBlock *vb, Context *lens_ctx, Context *values_ctx)
{
    // we adjust the buffer here, since it didn't get adjusted in piz_uncompress_all_ctxs because its ltype is LT_CODEC
    lens_ctx->local.len /= sizeof (uint32_t); 
    BGEN_u32_buf (&lens_ctx->local, NULL);

    ARRAY (uint32_t, next_of_chan, lens_ctx->local);

    // transform len array to next array
    uint32_t next=0;
    for (uint32_t chan=0; chan < next_of_chan_len; chan++) {
        uint32_t len = next_of_chan[chan];
        next_of_chan[chan] = next;
        next += len;
    }
    
    // retrieve the global value-to-bin mapper from SEC_COUNTS and store it in lens_ctx.con_cache
    buf_alloc (vb, &lens_ctx->con_cache, 0, 256, uint8_t, 0, "value_to_bin");
    lens_ctx->con_cache.len = 256;

    ARRAY (uint64_t, value_to_bin_src, ZCTX(lens_ctx->did_i)->counts); 
    ARRAY (uint8_t,  value_to_bin_dst, lens_ctx->con_cache);

    for (int i=0; i < 256; i++)
        value_to_bin_dst[i] = value_to_bin_src[i]; // uint64 -> uint8

    // initialize longr state - stored in values_ctx.con_cache
    buf_alloc_zero (vb, &values_ctx->con_cache, 1, 0, LongrState, 0, "contexts->local"); 
    codec_longr_alg_init (FIRSTENT (LongrState, values_ctx->con_cache));

    lens_ctx->is_initialized = true;
}

// When reconstructing a QUAL field on a specific line, piz calls the LT_CODEC reconstructor for CODEC_LONGR, 
// codec_domq_reconstruct, which combines data from the local buffers of lens and values to reconstruct the original QUAL field.
void codec_longr_reconstruct (VBlock *vb, Codec codec, Context *lens_ctx)
{
    ContextP values_ctx = lens_ctx + 1;
    ContextP seq_ctx    = CTX(VB_DT(DT_FASTQ) ? FASTQ_SQBITMAP : SAM_SQBITMAP);
    
    if (!lens_ctx->is_initialized) 
        codec_longr_reconstruct_init (vb, lens_ctx, values_ctx);

    bool is_rev = !VB_DT(DT_FASTQ) /* SAM or BAM */ && lens_ctx->dict_id.num == _SAM_QUAL && 
                  (CTX(SAM_FLAG)->last_value.i & SAM_FLAG_REV_COMP);

    ARRAY (uint8_t, sorted_qual, values_ctx->local);
    ARRAY (uint32_t, next_of_chan, lens_ctx->local);
    LongrState *state = FIRSTENT (LongrState, values_ctx->con_cache);
    state->value_to_bin = FIRSTENT (uint8_t, lens_ctx->con_cache);

    const char *seq = TXT_DT (DT_BAM) ? FIRSTENT (char, VB_SAM->textual_seq) : last_txtx (vb, seq_ctx); // note: textual_seq is prepared in sam_piz_sam2bam_SEQ

    codec_longr_recon_one_read (state, seq, vb->seq_len, is_rev, sorted_qual, next_of_chan, AFTERENT (char, vb->txt_data));
    vb->txt_data.len += vb->seq_len;
}

