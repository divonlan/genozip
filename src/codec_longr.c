// ------------------------------------------------------------------
//   codec_longr.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

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

#include "codec_longr_alg.c" // seperate source file for this, as it derived from external code with a different license

bool codec_longr_maybe_used (Did did_i)
{
    return (did_i == SAM_QUAL/*==FASTQ_QUAL*/ || did_i == SAM_CQUAL|| did_i == OPTION_OQ_Z) && 
           ((TECH(ONT) && segconf.nontrivial_qual && !flag.no_longr && !flag.fast) || flag.force_longr);
}

// similar structure to DOMQUAL 
void codec_longr_comp_init (VBlockP vb, Did qual_did_i)
{
    // lens_ctx contains the lens array - an array of uint32 - each entry is the length of the corresponding channel in values_ctx
    ContextP lens_ctx             = CTX(qual_did_i);
    lens_ctx->ltype               = LT_CODEC;    // causes reconstruction to go to codec_longr_reconstruct
    lens_ctx->lcodec              = CODEC_LONGR;

    // values_ctx contains the base quality data, sorted by channel
    ContextP values_ctx           = lens_ctx+1;  // used
    values_ctx->ltype             = LT_UINT8;
    values_ctx->local_dep         = DEP_L1; 
    values_ctx->lcodec            = CODEC_ARITH8;
    values_ctx->lcodec_hard_coded = true;
    values_ctx->counts_section    = true; // we store the global value-to-bin mapper in a SEC_COUNTS
    
    ctx_consolidate_stats (vb, qual_did_i, qual_did_i+1, DID_EOL);
}

// upper bound of compressed size of length array
uint32_t codec_longr_est_size (Codec codec, uint64_t uncompressed_len) 
{ 
    return codec_complex_est_size (CODEC_LONGR, LONGR_NUM_CHANNELS * sizeof (uint32_t));
}

static inline void add_to_histogram (uint32_t *histogram, bytes values, uint32_t values_len)
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

    ContextP zctx = ZCTX(ctx->did_i); // note: ctx is be predefinded, so zctx has the same did_i

    if (zctx->longr_bins_calculated) return; // bins already calculated - this happens in eg in Deep FASTQ, if bins were calculated by SAM before

    // create a histogram of values
    uint32_t histogram[256] = {};
    uint32_t num_values = 0;

    if (callback)
        for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {
            uint8_t *values; uint32_t values_len;
            callback (vb, ctx, line_i, (char **)pSTRa(values), Ltxt, NULL);
            add_to_histogram (histogram, STRa(values));

            num_values += values_len;
        }
    else {
        add_to_histogram (histogram, B1ST8 (ctx->local), ctx->local.len);
        num_values = ctx->local.len;
    }

    // create value_to_bin mapper in zctx->value_to_bin
    buf_alloc (evb, &zctx->value_to_bin, 0, 256, uint8_t, 0, "value_to_bin");
    uint8_t *value_to_bin = B1ST8 (zctx->value_to_bin);
    zctx->value_to_bin.len = NUM_BINS;
    
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
    buf_alloc (evb, &zctx->counts, 0, 256, uint64_t, 0, "zctx->counts");
    for (unsigned i=0; i < 256; i++)
        BNXT64 (zctx->counts) = value_to_bin[i];

    zctx->longr_bins_calculated = true;
}

static void codec_longr_calc_channels (LongrState *state, STRp(seq), bytes qual, bool is_rev) 
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

COMPRESS (codec_longr_compress)
{
    START_TIMER;
    ASSERTISNULL (uncompressed);
    ASSERT (soft_fail, "%s: second entry not expected, ctx=%s", VB_NAME, TAG_NAME);

    ContextP lens_ctx   = ctx; 
    ContextP values_ctx = ctx + 1;
    
    LocalGetLineCB *seq_callback = (VB_DT(FASTQ) ? fastq_zip_seq : sam_zip_seq);
    
    LongrState *state = codec_alloc (vb, sizeof (LongrState), 0);
    memset (state, 0, sizeof (LongrState));

    state->base_chan  = codec_alloc (vb, *uncompressed_len * sizeof (uint16_t), 0);
    state->value_to_bin = B1ST8 (ZCTX(values_ctx->did_i)->value_to_bin);
    codec_longr_alg_init (state);

    STRw(seq); STRw (qual); 
    bool is_rev;

    // calculate state->base_chan - the channel for each base. reads are treated in their
    // original (FASTQ) orientation.
    uint32_t total_len=0;
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {

        get_line_cb (vb, ctx, line_i, pSTRa(qual), Ltxt, NULL);
        if (!qual_len) continue; // this can happen, for example, if a SAM DEPN line is compressed against SA Group

        seq_callback (vb, ctx, line_i,  pSTRa(seq), CALLBACK_NO_SIZE_LIMIT, &is_rev);

        codec_longr_calc_channels (state, STRa(seq), (uint8_t*)qual, is_rev);

        ASSERT (seq_len == qual_len || str_is_1char(qual, ' '), "%s: \"%s\": Expecting seq_len=%u == qual_len=%u. ctx=%s", 
                LN_NAME, name, seq_len, qual_len, TAG_NAME);

        total_len += qual_len;
        ASSERT (total_len <= *uncompressed_len, "%s: \"%s\": Expecting total_len=%u <= total_len=%u. ctx=%s", 
                LN_NAME, name, total_len, *uncompressed_len, TAG_NAME);
    }
    
    // we now sort the quality data by channel (reversing qual of revcomp reads)
    buf_alloc (vb, &values_ctx->local, *uncompressed_len, 0, uint8_t, 0, CTX_TAG_LOCAL);
    values_ctx->local.len = *uncompressed_len;

    uint8_t *sorted_qual = B1ST8 (values_ctx->local);
    
    uint16_t *base_chan = state->base_chan;

    uint32_t next_of_chan[LONGR_NUM_CHANNELS];
    next_of_chan[0] = 0;
    for (uint32_t chan=1; chan < LONGR_NUM_CHANNELS; chan++) 
        next_of_chan[chan] = next_of_chan[chan-1] + state->chan_num_bases[chan-1];

    uint32_t base_i=0;
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {
        get_line_cb (vb, ctx, line_i, pSTRa(qual), Ltxt, &is_rev);

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
    buf_alloc (vb, &lens_ctx->local, LONGR_NUM_CHANNELS, 0, uint32_t, 0, CTX_TAG_LOCAL);
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
        ABORT ("%s: \"%s\": Compressing %s with %s need %u bytes, but allocated only %u", VB_NAME, name, lens_ctx->tag_name, codec_name(header->sub_codec), min_required_compressed_len, *compressed_len);

    // length is a numeric context, and as such it must be written first, as piz_uncompress_all_ctxs needs to complete piz_adjust_one_local
    COPY_TIMER_COMPRESS (compressor_longr);

    // actually compress the lens array
    return compress (vb, ctx, header, lens_ctx->local.data, uncompressed_len, NULL, compressed, compressed_len, false, name);
}

//--------------
// PIZ side
//--------------

static void codec_longr_recon_one_read (LongrState *state, STRp(seq), bool is_rev,
                                        bytes sorted_qual, uint32_t *next_of_chan, char *recon)
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
static void codec_longr_reconstruct_init (VBlockP vb, Context *lens_ctx, Context *values_ctx)
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
    
    // retrieve the global value-to-bin mapper from SEC_COUNTS and store it in values_ctx.value_to_bin
    buf_alloc (vb, &values_ctx->value_to_bin, 0, 256, uint8_t, 0, "value_to_bin");
    values_ctx->value_to_bin.len = 256;

    // until 13.0.11, stored in SEC_COUNTS of lens_ctx, and after in values_ctx
    BufferP value_to_bin = ZCTX(values_ctx->did_i)->counts.len ? &ZCTX(values_ctx->did_i)->counts 
                                                               : &ZCTX(lens_ctx->did_i)->counts;

    ARRAY (uint64_t, value_to_bin_src, *value_to_bin); 
    ARRAY (uint8_t,  value_to_bin_dst, values_ctx->value_to_bin);

    for (int i=0; i < 256; i++) 
        value_to_bin_dst[i] = value_to_bin_src[i]; // uint64 -> uint8

    // initialize longr state - stored in lens_ctx.longr_state
    buf_alloc_zero (vb, &lens_ctx->longr_state, 1, 0, LongrState, 0, CTX_TAG_LOCAL); 
    codec_longr_alg_init (B1ST (LongrState, lens_ctx->longr_state));

    lens_ctx->is_initialized = true;
}

// When reconstructing a QUAL field on a specific line, piz calls the LT_CODEC reconstructor for CODEC_LONGR, 
// codec_longr_reconstruct, which combines data from the local buffers of lens, values and SQBITMAP to reconstruct the original QUAL field.
CODEC_RECONSTRUCT (codec_longr_reconstruct)
{
    START_TIMER;

    ContextP lens_ctx   = ctx;
    ContextP values_ctx = ctx + 1;
    
    if (!lens_ctx->is_initialized) 
        codec_longr_reconstruct_init (vb, lens_ctx, values_ctx);

    bool is_rev = !VB_DT(FASTQ) /* SAM or BAM */ && sam_is_last_flags_rev_comp(vb);

    ARRAY (uint8_t, sorted_qual, values_ctx->local);
    ARRAY (uint32_t, next_of_chan, lens_ctx->local);
    LongrState *state = B1ST (LongrState, lens_ctx->longr_state);
    state->value_to_bin = B1ST8 (values_ctx->value_to_bin);

    rom seq = VB_DT(SAM) ? sam_piz_get_textual_seq(vb) : last_txtx (vb, CTX(FASTQ_SQBITMAP)); 

    // case: Deep, and len is only the trimmed suffix as the rest if copied from SAM (see fastq_special_deep_copy_QUAL)
    if (flag.deep && len < vb->seq_len)
        seq += (vb->seq_len - len); // advance seq to the trimmed part too
    else
        ASSPIZ (len == vb->seq_len, "expecting len=%u == vb->seq_len=%u", len, vb->seq_len);
    
    codec_longr_recon_one_read (state, seq, len, is_rev, sorted_qual, next_of_chan, BAFTtxt);
    Ltxt += len;

    COPY_TIMER(codec_longr_reconstruct);
}

