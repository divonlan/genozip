// ------------------------------------------------------------------
//   sam_ultima.c
//   Copyright (C) 2022-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "reconstruct.h"
#include "codec.h"

sSTRl(copy_RNAME_snip, 32);
sSTRl(copy_QNAME_snip, 32);

// used in ZIP and PIZ - part of the file format
#define TP_NUM_BINS 7
static const uint8_t tp_bins[256] = { [0 ... 38]=0, [39 ... 41]=1, [42 ... 59]=2, [60 ...62]=3, [63 ... 64]=4, [65 ... 72]=5, [73 ... 255]=6 };

void sam_ultima_zip_initialize (void)
{
    seg_prepare_snip_other (SNIP_COPY, _SAM_RNAME, 0, 0, copy_RNAME_snip);
    seg_prepare_snip_other (SNIP_COPY, _SAM_QNAME, 0, 0, copy_QNAME_snip);
}

void sam_ultima_seg_initialize (VBlockSAMP vb)
{
    ContextP arr_ctx = CTX(OPTION_tp_B_ARR);
    arr_ctx->ltype = LT_INT8; // propagated to channels
    arr_ctx->st_did_i = OPTION_tp_B_c; 
    arr_ctx->flags.store = STORE_INT; // needed for BAM translation (propagated to channels)
    seg_mux_init (VB, arr_ctx, TP_NUM_BINS, SAM_SPECIAL_DEMUX_BY_QUAL, false, (MultiplexerP)&vb->mux_tp); 

    seg_by_ctx (VB,STRa(vb->mux_tp.snip), arr_ctx, 0);  // all-the-same (not ctx_create node bc we need the b250 to carry the flags)

    if (segconf.sam_has_ultima_t0) {
        if (IS_DEPN(vb))
            CTX(OPTION_t0_Z)->ltype = LT_SEQUENCE; // see bug 922
        else
            codec_t0_comp_init (VB); 
    } else
        CTX(OPTION_t0_Z)->no_callback = true; // seg t0 normally

    if (segconf.running) {
        segconf.sam_has_ultima_t0 = true; // optimistic
    }
}

void sam_ultima_finalize_segconf (VBlockSAMP vb)
{
    segconf.sam_has_ultima_t0 = segconf.has[OPTION_t0_Z] && codec_t0_data_is_a_fit_for_t0(VB);
}

// example: tp:B:c,1,1,1,1,1,1,-1,2,0,0,0,2,-1,1,1,-1,2,0,0,0,0,2,-1,1,1,0,0,
void sam_seg_ultima_tp (VBlockSAMP vb, ContextP arr_ctx, void *dl_, void *tp_, uint32_t tp_len)
{
    START_TIMER;

    char *tp = (char *)tp_; // this is OPTION_tp_B_ARR.local - we transfer it to the channels and free it
    ZipDataLineSAM *dl = (ZipDataLineSAM *)dl_;

    rom qual = Btxt (dl->QUAL.index);

    ContextP chan[TP_NUM_BINS];
    for (int b=0; b < TP_NUM_BINS; b++)     
        chan[b] = seg_mux_get_channel_ctx (VB, OPTION_tp_B_ARR, (MultiplexerP)&vb->mux_tp, b); 

    char *next[TP_NUM_BINS-1], *start[TP_NUM_BINS-1];
    for (int b=0; b < TP_NUM_BINS-1; b++) {
        buf_alloc (VB, &chan[b]->local, tp_len, 0, int8_t, CTX_GROWTH, "local");
        next[b]  = BAFTc(chan[b]->local);
        start[b] = next[b];
    }

    for (uint32_t i=0; i < tp_len; i++) {
        uint8_t b = tp_bins[(uint8_t)qual[i]];
        if (b == 6) { // seg as a snip, as its expected to be all-the-same '0'
            if (tp[i] == 0)     
                seg_by_ctx (VB, "0", 1, chan[6], 0); // shortcut
            else 
                seg_integer_as_snip (vb, chan[6]->did_i, tp[i], false);
        }
        else 
            *next[b]++ = tp[i]; 
    }

    // If BAM: divvy is txt_len between the channels. note: this is a lot more difficult to do for SAM
    bool is_bam = IS_BAM_ZIP;
    if (is_bam) {
        chan[6]->txt_len += tp_len;
        arr_ctx->txt_len -= tp_len;
    }

    for (int b=0; b < TP_NUM_BINS-1; b++) {
        chan[b]->local.len32 = BNUM (chan[b]->local, next[b]);

        if (is_bam) {
            chan[b]->txt_len += next[b] - start[b];
            chan[b]->txt_len -= next[b] - start[b];
        }
    }

    arr_ctx->local.len32 = 0; 

    COPY_TIMER (sam_seg_ultima_tp);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_BY_QUAL)
{
    if (!ctx_encountered_in_line (vb, OPTION_tp_B_ARR)) 
        CTX(OPTION_tp_B_c)->last_repeat = -1;

    int32_t qual_i = ++CTX(OPTION_tp_B_c)->last_repeat;
    
    int8_t qual = last_txtx (vb, CTX(IS_PRIM(vb) ? SAM_QUALSA : SAM_QUAL))[qual_i] + (OUT_DT(BAM) ? 33 : 0);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), tp_bins[qual], new_value, reconstruct);
}

static uint64_t sam_ultima_bi_prediction (VBlockP vb)
{
    STRlast (qname, SAM_QNAME);

    // note: we parse the qname rather than taking the last_value of Q3NAME or Q4NAME, bc if copying
    // QNAME from buddy or from SAG, we don't have the values of Q?NAME
    SAFE_NULT(qname);
    rom hyphen = strrchr (qname, '-');
    uint64_t value = hyphen ? strtoll (hyphen+1, NULL, 10) : 0; // discards any prefix of '0's, stops at the first non-digit 
    SAFE_RESTORE;

    return value;
}

void sam_seg_ultima_bi (VBlockSAMP vb, STRp(bi_str), unsigned add_bytes)
{
    int64_t bi;

    if (str_get_int (STRa(bi_str), &bi) && (bi == sam_ultima_bi_prediction (VB))) 
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_bi }, 2, OPTION_bi_Z, add_bytes);

    else
        seg_integer_or_not (VB, CTX(OPTION_bi_Z), STRa(bi_str), add_bytes); 
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_bi)
{
    new_value->i = sam_ultima_bi_prediction (vb);

    if (reconstruct)
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

static bool sam_seg_ultima_delta_POS (VBlockP vb, ContextP ctx, STRp(pos_str), uint32_t rep)
{
    seg_pos_field (vb, ctx->did_i, SAM_POS, 0, 0, STRa(pos_str), 0, pos_str_len);

    return true; // segged successfully
}

static bool sam_seg_ultima_XV_AS (VBlockP vb, ContextP ctx, STRp(as_str), uint32_t rep)
{
    int64_t as;
    if (!str_get_int (STRa(as_str), &as)) return false;

    int32_t delta = as - DATA_LINE(vb->line_i)->AS; 

    SNIP(32);
    seg_prepare_snip_other (SNIP_OTHER_DELTA, _OPTION_AS_i, true, delta, snip);

    seg_by_ctx (VB, STRa(snip), ctx, as_str_len);
    
    return true; // segged successfully
}

static bool sam_seg_ultima_XV_MAPQ (VBlockP vb, ContextP ctx, STRp(mapq_str), uint32_t rep)
{
    int64_t mapq;
    if (!str_get_int (STRa(mapq_str), &mapq)) return false;

    int32_t delta = mapq - DATA_LINE(vb->line_i)->MAPQ; 

    SNIP(32);
    seg_prepare_snip_other (SNIP_OTHER_DELTA, _SAM_MAPQ, true, delta, snip);

    seg_by_ctx (VB, STRa(snip), ctx, mapq_str_len);

    return true; // segged successfully
}

// example: XV:Z:chr1,25173438,306,60
void sam_seg_ultima_XV (VBlockSAMP vb, STRp(xv), unsigned add_bytes)
{
    static const MediumContainer container_XV = {
        .repeats      = 1, 
        .nitems_lo    = 4, 
        .items        = { { .dict_id.num = DICT_ID_MAKE2_L("X0V_RNAME"), .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_7("X1V_POS"),   .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_6("X2V_AS"),    .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_L("X3V_MAPQ")                    } } };

    SegCallback callbacks[4] = { 0, sam_seg_ultima_delta_POS, sam_seg_ultima_XV_AS, sam_seg_ultima_XV_MAPQ }; 

    seg_struct (VB, CTX(OPTION_XV_Z), container_XV, STRa(xv), callbacks, add_bytes, true);
}

// describes the first few mismatches - except for .repeats, it can be 100% determined by MD:Z and SEQ
// example: XW:Z:chr1,9206568,G,A;chr1,9206575,C,T;
void sam_seg_ultima_XW (VBlockSAMP vb, STRp(xw), unsigned add_bytes)
{
    static const MediumContainer container_XW = { // bug 892
        .repsep       = ";", 
        .nitems_lo    = 4, 
        .items        = { { .dict_id.num = DICT_ID_MAKE2_L("X0W_RNAME"), .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_7("X1W_POS"),   .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_6("X2V_REF"),   .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_L("X3V_ALT")                     } } };

    SegCallback callbacks[4] = { 0, sam_seg_ultima_delta_POS }; 

    seg_array_of_struct (VB, CTX(OPTION_XW_Z), container_XW, STRa(xw), callbacks, add_bytes);
}

// t0:Z : supplemental base quality information
// example: =77+**1119955))),,,..=I///;;;222888*****:;;>>>>AA<<<<IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII@@@@IIAAAA<<<<
void sam_seg_ultima_t0 (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(t0), unsigned add_bytes)    
{                                                          
    decl_ctx (OPTION_t0_Z);

    ASSINP (t0_len == dl->SEQ.len, "%s: Expecting t0.len=%u == SEQ.len=%u. %s=\"%.*s\"",       
            LN_NAME, t0_len, dl->SEQ.len, ctx->tag_name, STRf(t0));     

    dl->t0 = TXTWORD(t0);                                                           
    ctx->txt_len += add_bytes;                      
    ctx->local.len32 += t0_len;                  

    // next: compression by T0 codec
}              

// callback function for compress to get data of one line (used by T0 codec)
COMPRESSOR_CALLBACK (sam_zip_t0) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = MIN_(maximum_size, dl->t0.len);

    if (!line_data) return; // only lengths were requested

    *line_data = Btxt (dl->t0.index);
}

void sam_ultima_update_t0_len (VBlockP vb, uint32_t line_i, uint32_t new_len) 
{ 
    DATA_LINE(line_i)->t0.len = new_len; 
}

// MI:Z
void sam_seg_ultima_MI (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(mi), unsigned add_bytes)    
{
    // TODO: needs to be implemented with lookback as this condition fails in ~70% of cases
    if (str_issame_(STRa(mi), STRtxtw(dl->QNAME)))
        seg_by_did (VB, STRa(copy_QNAME_snip), OPTION_MI_Z, add_bytes);

    else
        seg_by_did (VB, STRa(mi), OPTION_MI_Z, add_bytes);
}

// used by SAM and FASTQ
void ultima_c_Q5NAME_cb (VBlockP vb, ContextP ctx, STRp(value))
{
    bool is_fq = VB_DT(FASTQ);
    Multiplexer2P mux = is_fq ? fastq_get_ultima_c_mux (vb) : &VB_SAM->mux_ultima_c;
    
    if (!ctx->is_initialized) {
        seg_mux_init (VB, ctx, 2, (is_fq ? FASTQ_SPECIAL_ULTIMA_C : SAM_SPECIAL_ULTIMA_C), false, (MultiplexerP)mux);

        seg_by_ctx (VB, STRa(mux->snip), ctx, 0);  // all-the-same (not ctx_create node bc we need the b250 to carry the flags)
        
        ctx->is_initialized = true;
        
        (ctx-1)->flags.store = STORE_INT;
    }

    int channel_i = ctx_has_value_in_line_(vb, ctx-1) && ((ctx-1)->last_value.i == 2);

    ContextP channel_ctx = seg_mux_get_channel_ctx (vb, ctx->did_i, (MultiplexerP)mux, channel_i);

    seg_integer_or_not (vb, channel_ctx, STRa(value), value_len);
}

SPECIAL_RECONSTRUCTOR (ultima_c_piz_special_DEMUX_BY_Q4NAME)
{
    int channel_i = ctx_has_value_in_line_(vb, ctx-1) && ((ctx-1)->last_value.i == 2);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);    
}
