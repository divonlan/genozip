// ------------------------------------------------------------------
//   codec_pacb.c
//   Copyright (C) 2023-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "codec.h"
#include "reconstruct.h"
#include "profiler.h"
#include "compressor.h"
#include "sam.h"
#include "container.h"
#include "endianness.h"

// ZIP/PIZ: contribution of the SEQ environment (seq_i-1 to seq_i+2) to channel_i. This is part of codec / file format
#define NUM_Ks 7
static inline uint8_t QUAL_get_K_value (STRp(seq), uint32_t i)
{
    char base = seq[i];

    if (i > 0 && seq[i-1] == base)          return 6;                            // 2nd+ base in homopolymer
    if (i == seq_len-1 || seq[i+1] != base) return 4 + (base=='A' || base=='T'); // not a homopolymer 
    if (i == seq_len-2 || seq[i+2] != base) return 2 + (base=='A' || base=='T'); // first base in a homo-2-mer
    else                                    return 0 + (base=='A' || base=='T'); // first base in a homo-3-mer or higher
}

static void codec_pacb_init_ctxs (VBlockP vb, ContextP ctx, uint8_t n_channels, bool also_init_zctx, ContextP *subctxs/*alloced by caller*/)
{
    ContextP zctx = ZCTX(ctx->did_i);
    DictId *pb_dicts = B1ST (DictId, zctx->subdicts);

    for (uint8_t channel_i=0; channel_i < n_channels; channel_i++) {
        ContextP subctx = IS_ZIP ? ctx_get_ctx (vb, pb_dicts[channel_i]) : ECTX(pb_dicts[channel_i]);

        if (IS_ZIP) {
            subctx->local_dep = DEP_L1;
            subctx->st_did_i  = ctx->did_i;
        }

        if (subctxs) subctxs[channel_i] = subctx;

        if (also_init_zctx)
            ctx_get_zctx_from_vctx (subctx, true, false);
    }
}

//--------------
// ZIP side
//--------------

#define MAX_np 12 // ZIP: maximum value of np:i. Values beyond the maximum are treated as the maximum

static uint8_t get_max_np (void)
{
    return (TXT_DT(FASTQ) || !segconf.has[OPTION_np_i] /* CLR */) ? 1 : MAX_np;    
}

void codec_pacb_segconf_finalize (VBlockP vb)
{
    ContextP ctx  = CTX(SAM_QUAL);
    ContextP zctx = ZCTX(SAM_QUAL);
    uint8_t n_channels = get_max_np() * NUM_Ks;

    // store the dict_id array in the file as a global SEC_SUBDICTS section 
    ARRAY_alloc (DictId, subdicts, n_channels, false, zctx->subdicts, evb, "zctx->subdicts");

    // create dict_id for each channel  
    bytes id = dict_id_typeless (ctx->dict_id).id;

    for (uint8_t channel_i = 0; channel_i < n_channels; channel_i++)
        subdicts[channel_i] = dict_id_make((char[]){ id[0], '0'+channel_i/10, '0'+channel_i%10, '-',id[0], id[1], id[2], id[3]}, 8, dict_id_type (ctx->dict_id));

    // create zctx's for ctx_update_stats to work (since our subcontexts are not pre-defined) 
    codec_pacb_init_ctxs (vb, ctx, n_channels, true, NULL);

    zctx->subdicts_section = true; // we store the dict_id array in SEC_SUBDICTS (global section)
}

bool codec_pacb_comp_init (VBlockP vb, LocalGetLineCB get_line_cb)
{
    ContextP ctx = CTX(SAM_QUAL); // ==FASTQ_QUAL

    uint32_t (*get_seq_len)(VBlockP, uint32_t) = (VB_DT(FASTQ) ? fastq_zip_get_seq_len : sam_zip_get_seq_len);

    // make sure score is same length as SEQ, or 0 (i.e. compressed with another method), for all lines
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {   
        uint32_t score_len, seq_len;
        get_line_cb (vb, ctx, line_i, NULL, &score_len, 0, NULL);
        seq_len = get_seq_len (vb, line_i);

        if (score_len != seq_len && score_len != 0) 
            return false;
    }

    ctx->ltype  = LT_CODEC;
    ctx->lcodec = CODEC_PACB;
    
    return true;
}

// ZIP: called for QUAL-like dids
bool codec_pacb_maybe_used (Did did_i)
{
    return (TECH(PACBIO) && segconf.nontrivial_qual && !segconf.use_pacbio_iqsqdq && did_i == SAM_QUAL/*==FASTQ_QUAL*/);
}

// ZIP: calculate the channel_i for each score on the line based on its environment (SEQ)
static uint8_t *calc_channels_one_line (VBlockP vb, LineIType line_i, uint8_t np0, uint8_t *channel_i_p, uint32_t *lens)
{
    STRw(seq);
    (VB_DT(FASTQ) ? fastq_zip_seq : sam_zip_seq) (vb, NULL, line_i, pSTRa(seq), CALLBACK_NO_SIZE_LIMIT, NULL);

    for (uint32_t i=0; i < seq_len; i++) {
        uint8_t K = QUAL_get_K_value (STRa(seq), i);

        uint8_t channel_i = (np0 * NUM_Ks + K);
        *channel_i_p++ = channel_i; // channel_i where score i on line line_i will go
        lens[channel_i]++;
    }

    return channel_i_p;
}

static uint8_t *set_values_one_line (VBlockP vb, STRp(qual),
                                     uint8_t *channel_i_p, // channel into which each score should go
                                     char *ch_next[],      // array of pointers - "next" for each channel
                                     rom ch_after[])       // array of pointers - "after" for each channel
{
    rom after = qual + qual_len;
    while (qual < after) {
        uint8_t channel_i = *channel_i_p++; 

        ASSERT (ch_next[channel_i] < ch_after[channel_i], "%s: QUAL: channel_i=%u out of space", VB_NAME, channel_i);

        *ch_next[channel_i]++ = *qual++; 
    }

    return channel_i_p;
}

uint32_t codec_pacb_est_size (Codec codec, uint64_t uncompressed_len)
{
    return 1; // QUAL.compressed_len is one byte
}

COMPRESS (codec_pacb_compress)
{
    START_TIMER;
    
    uint8_t max_np = get_max_np();
    uint8_t n_channels = max_np * NUM_Ks;

    ContextP subctxs[n_channels];
    codec_pacb_init_ctxs (vb, ctx, n_channels, false, subctxs);

    uint32_t lens[n_channels];
    memset (lens, 0, n_channels * sizeof (uint32_t));
 
    buf_alloc (vb, &vb->codec_bufs[0], 0, ctx->local.len32, uint8_t, 0, "codec_bufs[0]");
    uint8_t *channel_i_p = B1ST8 (vb->codec_bufs[0]);

    // first pass - get the channel_i of each score into channels, and calculate lengths
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {   
        STRw(score);
        get_line_cb (vb, ctx, line_i, pSTRa (score), CALLBACK_NO_SIZE_LIMIT, NULL);

        if (!score_len) continue; // compressed by another method

        uint8_t np0 = (max_np > 1) ? MIN_(sam_zip_get_np (vb, line_i), max_np) - 1 : 0; // notes: np0 is 0-based np, i.e. (np-1) ; np0 is always 0 for FASTQ and CLR

        channel_i_p = calc_channels_one_line (vb, line_i, np0, channel_i_p, lens);
    }

    char *next[n_channels]; // pointer into next of each channel
    rom after[n_channels];

    // allocate local buffers
    for (uint8_t channel_i=0; channel_i < n_channels; channel_i++) {
        ContextP subctx = subctxs[channel_i];

        buf_alloc_exact (vb, subctx->local, lens[channel_i], char, CTX_TAG_LOCAL);

        subctx->txt_len += lens[channel_i];       // update accounting of the data to the channel context
        ctx->txt_len    -= lens[channel_i];

        next [channel_i] = B1STc (subctx->local); // initialize
        after[channel_i] = BAFTc (subctx->local);
    }

    // second pass: multiplex: copy each score to its designated channel
    channel_i_p = B1ST8 (vb->codec_bufs[0]);

    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {   
        STRw(score);
        get_line_cb (vb, ctx, line_i, pSTRa (score), CALLBACK_NO_SIZE_LIMIT, NULL);

        if (!score_len) continue; // compressed by another method

        channel_i_p = set_values_one_line (vb, STRa(score), channel_i_p, next, after);
    }

    buf_free (vb->codec_bufs[0]);

    // write one byte to the SEC_LOCAL section, just to trigger reconstruction
    header->sub_codec = CODEC_NONE;
    *compressed = 0;
    *uncompressed_len = *compressed_len = 1;
    
    COPY_TIMER_COMPRESS (compressor_pacb); 
    return true;
}

//--------------
// PIZ side
//--------------

static int32_t codec_pacb_piz_get_np (VBlockP vb)
{
    // case: AUX already recontructed (case 1: TOP2FQEX - reconstructed on AUX DESC line, case2: TOP2FQ reconstructed above)
    if (ctx_has_value_in_line_(vb, CTX(OPTION_np_i))) 
        return CTX(OPTION_np_i)->last_value.i;

    // case : AUX will be reconstructed later - just peek np:i now
    else 
        return reconstruct_peek (vb, CTX(OPTION_np_i), 0, 0).i;
}

CODEC_RECONSTRUCT (codec_pacb_reconstruct)
{
    START_TIMER;

    uint8_t max_np = ZCTX(ctx->did_i)->subdicts.len32 / NUM_Ks; 
    uint8_t n_channels = max_np * NUM_Ks;

    // case first recon of this field in this VB: initialize channels from dict_ids stored in zctz->subdicts
    if (!ctx->is_initialized) {
        buf_alloc_exact (vb, ctx->channel_data, n_channels * 2, rom, "channel_data"); // 2 arrays - "next" and "after"
        codec_pacb_init_ctxs (vb, ctx, n_channels, false, (void*)B1STc(ctx->channel_data));

        // set channel_data to hold the "next" and "after" for of each channel
        for (int i=0; i < n_channels; i++) {
            BufferP channel_buf = &(*B(ContextP, ctx->channel_data, i))->local;
            *B(rom, ctx->channel_data, i) = B1STc(*channel_buf);
            *B(rom, ctx->channel_data, n_channels + i) = BAFTc(*channel_buf); 
        }

        ctx->is_initialized = true;
    }

    int32_t np = (max_np > 1) ? codec_pacb_piz_get_np (vb) : 1;
    uint8_t np0 = MIN_(np, max_np) - 1;

    // arrays of the "next" score in each channel, and "after" for each channel 
    rom *next  = B1ST(rom, ctx->channel_data); // an array of rom
    const rom *after = B(rom, ctx->channel_data, n_channels);
    
    // reconstruct 
    ContextP seq_ctx = CTX(SAM_SQBITMAP); // same as FASTQ_SQBITMAP

    ConstBufferP textual_seq = VB_DT(SAM) ? sam_get_textual_seq(vb) : NULL; // note: textual_seq is prepared in sam_piz_sam2bam_SEQ sam_load_groups_add_seq

    rom seq = (textual_seq && textual_seq->len32) ? B1STc (*textual_seq) : last_txtx (vb, seq_ctx); 

    char *next_recon = BAFTtxt;
    for (uint32_t i=0; i < len; i++) {
        uint8_t K = QUAL_get_K_value (seq, len, i);

        uint8_t channel_i = NUM_Ks * np0 + K;
        ASSPIZ (next[channel_i] < after[channel_i], "out of data in channel_i=%u", channel_i);

        *next_recon++ = *next[channel_i]++;
    }

    if (reconstruct) Ltxt += len;

    COPY_TIMER(codec_pacb_reconstruct);
}
