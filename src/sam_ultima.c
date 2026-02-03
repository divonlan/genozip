// ------------------------------------------------------------------
//   sam_ultima.c
//   Copyright (C) 2022-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <math.h>
#include "sam_private.h"
#include "codec.h"

sSTRl(copy_RNAME_snip, 32);
sSTRl(copy_QNAME_snip, 32);

// used in ZIP and PIZ - part of the file format
#define TP_NUM_BINS 7
static const uint8_t tp_bins[256] = { [0 ... 38]=0, [39 ... 41]=1, [42 ... 59]=2, [60 ...62]=3, [63 ... 64]=4, [65 ... 72]=5, [73 ... 255]=6 };

#define TP_NON_CONDENSED 127

void sam_ultima_zip_initialize (void)
{
    DO_ONCE {
        seg_prepare_snip_other (SNIP_COPY, _SAM_RNAME, 0, 0, copy_RNAME_snip);
        seg_prepare_snip_other (SNIP_COPY, _SAM_QNAME, 0, 0, copy_QNAME_snip);
    }
}

void sam_ultima_seg_initialize (VBlockSAMP vb)
{
    ContextP arr_ctx = CTX(OPTION_tp_B_ARR);
    arr_ctx->ltype = LT_INT8; // propagated to channels
    arr_ctx->st_did_i = OPTION_tp_B_c; 
    arr_ctx->flags.store = STORE_INT; // needed for BAM translation (propagated to channels)
    seg_mux_init (vb, OPTION_tp_B_ARR, SAM_SPECIAL_ULTIMA_tp, false, tp); 

    seg_by_ctx (VB,STRa(vb->mux_tp.snip), arr_ctx, 0);  // all-the-same (not ctx_create node bc we need the b250 to carry the flags)

    if (segconf.sam_has_ultima_t0) {
        if (IS_DEPN(vb))
            CTX(OPTION_t0_Z)->ltype = LT_BLOB; // see bug 922
        else 
            codec_t0_comp_init (VB);
    } else
        CTX(OPTION_t0_Z)->no_callback = true; // seg t0 normally

    if (segconf_running) {
        segconf.sam_has_ultima_t0 = true; // optimistic
    }
}

void sam_ultima_finalize_segconf (VBlockSAMP vb)
{
    segconf.sam_has_ultima_t0 = segconf.has[OPTION_t0_Z] && codec_t0_data_is_a_fit_for_t0(VB);
}

static inline void seg_one_tp (VBlockSAMP vb, ContextP chan[TP_NUM_BINS], char *next[TP_NUM_BINS-1], char tp, char qual)
{
    uint8_t bin = tp_bins[(uint8_t)qual];
    if (bin == 6) { // seg as a snip, as its expected to be all-the-same '0'
        if (tp == 0)     
            seg_by_ctx (VB, "0", 1, chan[6], 0); // shortcut
        else 
            seg_integer_as_snip (vb, chan[6]->did_i, tp, false);
    }
    else 
        *next[bin]++ = tp; 
}

// example: tp:B:c,1,1,1,1,1,1,-1,2,0,0,0,2,-1,1,1,-1,2,0,0,0,0,2,-1,1,1,0,0,
ARRAY_ITEM_CALLBACK (sam_seg_ULTIMA_tp)
{
    START_TIMER;

    char *tp = (char *)array; // this is OPTION_tp_B_ARR.local - we transfer it to the channels and free it
    uint32_t tp_len = array_len;
    ZipDataLineSAMP dl = (ZipDataLineSAMP)cb_param;

    // note: flags.no_textual_seq is set to false, in sam_seg_finalize, as required by sam_piz_special_ULTIMA_tp
    ASSSEG (tp_len == dl->SEQ.len, "Expecting length of TP:B:c to be seq_len=%u but it is %u", dl->SEQ.len, tp_len);
    ASSSEG (tp_len < (1 << TP_LEN_BITS), "Genozip limitation: tp:B:c is supported up to length %u, but here tp_len=%u. %s", (1 << TP_LEN_BITS)-1, tp_len, report_support_if_unexpected());

    rom seq = vb->textual_seq_str;
    rom qual = Btxt (dl->QUAL.index);

    ContextP chan[TP_NUM_BINS];
    for (int bin=0; bin < TP_NUM_BINS; bin++)     
        chan[bin] = seg_mux_get_channel_ctx (VB, OPTION_tp_B_ARR, (MultiplexerP)&vb->mux_tp, bin); 

    char *next[TP_NUM_BINS-1];
    for (int bin=0; bin < TP_NUM_BINS-1; bin++) {
        buf_alloc (VB, &chan[bin]->local, tp_len, 0, int8_t, CTX_GROWTH, CTX_TAG_LOCAL);
        next[bin] = BAFTc(chan[bin]->local);
    }

    for (uint32_t i=0; i < dl->SEQ.len; i++) {
        unsigned hp_len = homopolymer_len (seq, dl->SEQ.len, i);  
        
        if (hp_len > 1) {
            ASSSEG (hp_len < (1 << HP_LEN_BITS), "Genozip limitation: homopolymers are supported up to length %u, but here tp_len=%u. %s", (1 << HP_LEN_BITS)-1, hp_len, report_support_if_unexpected());
        
            unsigned num_to_be_segged = (hp_len + 1) / 2; // optimistic

            // test condensability (rarely, homopolymer tp may be non-condensable)
            for (unsigned h=0; h < num_to_be_segged; h++) 
                if (tp[i + h] != tp[i + hp_len - 1 - h]) { // this tp doesn't equal its mirror -> not condensable
                    num_to_be_segged = hp_len; // unfortunately we have to seg the full length
                    seg_one_tp (vb, chan, next, TP_NON_CONDENSED, qual[i]); // before seg of first value - at offset i
                    break;
                }

            // actual seg - partial or full length, depending on condensablity
            for (unsigned h=0; h < num_to_be_segged; h++)
                seg_one_tp (vb, chan, next, tp[i + h], qual[i + h]);

            i += hp_len - 1; // skip to end of homopolymer
        }

        else
            seg_one_tp (vb, chan, next, tp[i], qual[i]);
    }

    for (int bin=0; bin < TP_NUM_BINS-1; bin++) 
        chan[bin]->local.len32 = BNUM (chan[bin]->local, next[bin]);

    ctx->local.len32 = 0; 

    COPY_TIMER (sam_seg_ULTIMA_tp);
}

// 15.0.10 to 15.0.27, see also bug 959
SPECIAL_RECONSTRUCTOR (sam_piz_special_ULTIMA_tp_old)
{
    int32_t qual_i = ctx->tp.repeat++;
    
    int8_t qual = last_txtx (vb, CTX(IS_PRIM(vb) ? SAM_QUALSA : SAM_QUAL))[qual_i] + (OUT_DT(SAM) ? 0 : 33);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), tp_bins[qual], new_value, reconstruct);
}

// since 15.0.28
SPECIAL_RECONSTRUCTOR (sam_piz_special_ULTIMA_tp)
{
    int32_t qual_i = ctx->tp.repeat++; // note: tp initialized in sam_reset_line
    
    if (!qual_i)
        buf_alloc_exact (vb, ctx->homopolymer, vb->seq_len, int8_t, "homopolymer");

    // case: homopolymer already ended
    if (ctx->tp.hp_len && qual_i - ctx->tp.hp_start >= ctx->tp.hp_len) {
        ctx->tp.hp_len = 0;
        ctx->tp.condensed = false;
    }

    // case: not in homopolymer, check if we should start a new one
    if (!ctx->tp.hp_len) {
        uint32_t hp_len = homopolymer_len (sam_piz_get_textual_seq (vb), vb->seq_len, qual_i);

        if (hp_len > 1) {
            ctx->tp.hp_start  = qual_i;
            ctx->tp.hp_len    = hp_len;
            ctx->tp.condensed = true; // optimisitc
        }
    }

    // case: we're in the second half of a condensed homopolymer - copy from mirror
    uint32_t hp_i, mid_i;
    if (ctx->tp.condensed && ((hp_i = (qual_i - ctx->tp.hp_start)) >= (mid_i = (ctx->tp.hp_len+1) / 2))) {
        int32_t half2_i = hp_i - mid_i; // index within second (copied) half of homopolymer
        int32_t last1_i = ctx->tp.hp_len / 2 - 1; // index of last first half, ignoring middle (non-copied) tp if odd length
        new_value->i = *B(int8_t, ctx->homopolymer, last1_i - half2_i);
        if (reconstruct) RECONSTRUCT_INT (new_value->i);
    }

    else {
        int8_t qual = last_txtx (vb, CTX(IS_PRIM(vb) ? SAM_QUALSA : SAM_QUAL))[qual_i] + (OUT_DT(SAM) ? 0 : 33);
        int channel_i = tp_bins[qual];

        uint32_t save_ltxt = Ltxt;
        reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
    
        if (new_value->i == TP_NON_CONDENSED) { // override homopolymer
            ctx->tp.condensed = false; // actually, its not condensed
            Ltxt = save_ltxt; // rollback
            reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
        }

        if (ctx->tp.condensed) // first half of condensed
            *B(int8_t, ctx->homopolymer, qual_i - ctx->tp.hp_start) = new_value->i;
    }

    return HAS_NEW_VALUE;
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
        seg_special0 (VB, SAM_SPECIAL_bi, CTX(OPTION_bi_Z), add_bytes);

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
        .items        = { { .dict_id.num = DICT_ID_MAKE2_8("X0V_RNAM"), .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_7("X1V_POS"),  .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_6("X2V_AS"),   .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_8("X3V_MAPQ")                   } } };

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
        .items        = { { .dict_id.num = DICT_ID_MAKE2_8("X0W_RNAM"), .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_7("X1W_POS"),  .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_7("X2W_REF"),  .separator = "," }, 
                          { .dict_id.num = DICT_ID_MAKE2_7("X3W_ALT")                     } } };

    SegCallback callbacks[4] = { 0, sam_seg_ultima_delta_POS }; 

    seg_array_of_struct (VB, CTX(OPTION_XW_Z), container_XW, STRa(xw), callbacks, 
                         segconf.sam_semcol_in_contig ? sam_seg_correct_for_semcol_in_contig : NULL,
                         add_bytes); 
}

// t0:Z : supplemental base quality information
// example: =77+**1119955))),,,..=I///;;;222888*****:;;>>>>AA<<<<IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII@@@@IIAAAA<<<<
void sam_seg_ultima_t0 (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(t0), unsigned add_bytes)    
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
    ZipDataLineSAMP dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = MIN_(maximum_size, dl->t0.len);

    if (!line_data) return; // only lengths were requested

    *line_data = Btxt (dl->t0.index);
}

void sam_ultima_update_t0_len (VBlockP vb, uint32_t line_i, uint32_t new_len) 
{ 
    DATA_LINE(line_i)->t0.len = new_len; 
}

// MI:Z - QNAME format: if its non-duplicate - same as QNAME. If its duplicate - MI equals 
// to the QNAME of the corresponding non-duplicate.
void sam_seg_ultima_MI (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(mi), unsigned add_bytes)    
{
    #define MI_HISTORY_LEN 64
    decl_ctx (OPTION_MI_Z);

    if (!ctx->mi_history.len32)
        buf_alloc_exact_255 (VB, ctx->mi_history, MI_HISTORY_LEN, TxtWord, "mi_history");

    ARRAY (TxtWord, mi_history, ctx->mi_history);

    if (!dl->FLAG.duplicate) {
        if (str_issame_(STRa(mi), STRqname(dl)))
            seg_special1 (VB, SAM_SPECIAL_ULTIMA_mi, '0', ctx, add_bytes);
        else
            goto fallback;
    }

    else {
        // search for a non-duplicaete or duplicate in the same MI-group that appears before
        int i=0; for (; i < MI_HISTORY_LEN; i++)
            if (str_issame_(STRa(mi), STRtxt(mi_history[i]))) {
                seg_special1 (VB, SAM_SPECIAL_ULTIMA_mi, '0' + i, ctx, add_bytes);
                break;
            }

        // no recently-previous read has the same MI (duplicate or not)
        if (i == MI_HISTORY_LEN) fallback: {
            SAFE_ASSIGNx(mi-3, SNIP_SPECIAL, 0);
            SAFE_ASSIGNx(mi-2, SAM_SPECIAL_ULTIMA_mi, 1);
            SAFE_ASSIGNx(mi-1, '/', 2);
            
            seg_by_ctx (VB, mi - 3, mi_len + 3, ctx, add_bytes);

            SAFE_RESTOREx(0); SAFE_RESTOREx(1); SAFE_RESTOREx(2);
        }
    }

    // shift history up, and add mi at entry 0 
    memmove (&mi_history[1], &mi_history[0], (MI_HISTORY_LEN-1) * sizeof(TxtWord));
    mi_history[0] = (TxtWord){ .index = BNUMtxt (mi), .len = mi_len }; 
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_ULTIMA_MI)
{
    if (!ctx->history.len32)
        buf_alloc_exact_255 (VB, ctx->history, MI_HISTORY_LEN, TxtWord, "history");

    ARRAY (TxtWord, mi_history, ctx->history);

    uint32_t recon_index = Ltxt;

    if (reconstruct) {
        if (!last_flags.duplicate && snip[0] == '0') {
            STRlast (qname, SAM_QNAME);
            RECONSTRUCT_str (qname);
        }

        else if (snip[0] != '/') {
            int back = snip[0] - '0';
            ASSPIZ (back < MI_HISTORY_LEN, "Expecting back=%u < MI_HISTORY_LEN=%u", back, MI_HISTORY_LEN);
            ASSPIZ (mi_history[back].index != 0xffffffff, "MI history at back=%u is not populated yet", back);

            RECONSTRUCT (Btxt(mi_history[back].index), mi_history[back].len);   
        }

        else
            RECONSTRUCT (snip+1, snip_len-1);
    }

    // shift history up, and add mi at entry 0 
    memmove (&mi_history[1], &mi_history[0], (MI_HISTORY_LEN-1) * sizeof(TxtWord));
    mi_history[0] = (TxtWord){ .index = recon_index, .len = Ltxt - recon_index }; 

    return NO_NEW_VALUE;
}

// a3:Z
void sam_seg_ultima_a3 (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t a3, unsigned add_bytes)    
{
    SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_delta_seq_len, (int64_t)dl->SEQ.len - a3);
    seg_by_did (VB, STRa(snip), OPTION_a3_i, add_bytes);
}

// used by SAM and FASTQ
void ultima_c_Q5NAME_cb (VBlockP vb, ContextP ctx, STRp(value))
{
    bool is_fq = VB_DT(FASTQ);
    Multiplexer2P mux = is_fq ? fastq_get_ultima_c_mux (vb) : &VB_SAM->mux_ultima_c;
    
    if (!ctx->is_initialized) {
        seg_mux_init_(VB, ctx->did_i, 2, (is_fq ? FASTQ_SPECIAL_ULTIMA_C : SAM_SPECIAL_ULTIMA_C), false, (MultiplexerP)mux);

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
