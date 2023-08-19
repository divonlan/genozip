// ------------------------------------------------------------------
//   vcf_samples.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "context.h"
#include "optimize.h"
#include "file.h"
#include "strings.h"
#include "dict_id.h"
#include "codec.h"
#include "reconstruct.h"
#include "base64.h"
#include "stats.h"
#include "piz.h"
#include "zfile.h"
#include "zip.h"
#include "lookback.h"

static SmallContainer con_AD={}, con_ADALL={}, con_ADF={}, con_ADR={}, con_SAC={}, con_F1R2={}, con_F2R1={}, 
    con_MB={}, con_SB={}, con_AF={};

static char sb_snips[2][32], mb_snips[2][32], f2r1_snips[VCF_MAX_ARRAY_ITEMS][32], adr_snips[VCF_MAX_ARRAY_ITEMS][32], adf_snips[VCF_MAX_ARRAY_ITEMS][32], 
    rdf_snip[32], rdr_snip[32], adf_snip[32], adr_snip[32], ad_varscan_snip[32], ab_snip[48], gq_by_pl[50], gq_by_gp[50],
    sac_snips[VCF_MAX_ARRAY_ITEMS/2][32], PL_to_PLn_redirect_snip[30], PL_to_PLy_redirect_snip[30];
STRl(af_snip,32);

static unsigned sb_snip_lens[2], mb_snip_lens[2], f2r1_snip_lens[VCF_MAX_ARRAY_ITEMS], adr_snip_lens[VCF_MAX_ARRAY_ITEMS], adf_snip_lens[VCF_MAX_ARRAY_ITEMS], 
    sac_snip_lens[VCF_MAX_ARRAY_ITEMS/2], rdf_snip_len, rdr_snip_len, adf_snip_len, adr_snip_len, ad_varscan_snip_len,
    ab_snip_len, gq_by_pl_len, gq_by_gp_len, PL_to_PLn_redirect_snip_len, PL_to_PLy_redirect_snip_len;

static DictId make_array_item_dict_id (uint64_t dict_id_num, unsigned item_i)
{
    bytes id = ((DictId)dict_id_num).id;
    char dict_id_str[8] = { id[0], base32(item_i), id[1], id[2], id[3], id[4], id[5], id[6] };
    
    return dict_id_make (dict_id_str, 8, DTYPE_2);
}

static SmallContainer vcf_samples_init_container_array (uint64_t dict_id_num)
{
    SmallContainer con = (SmallContainer){ .repeats = 1, .drop_final_item_sep = true };

    for (unsigned i=0; i < SMALL_CON_NITEMS; i++) {
        con.items[i].dict_id      = make_array_item_dict_id (dict_id_num, i);
        con.items[i].separator[0] = ',';
    }

    return con;
}

void vcf_samples_zip_initialize (void) 
{
    static bool done = false;
    if (done) return; // already initialized (in previous files)
    done = true;

    con_AF    = vcf_samples_init_container_array (_FORMAT_AF);    
    con_AD    = vcf_samples_init_container_array (_FORMAT_AD);    
    con_ADALL = vcf_samples_init_container_array (_FORMAT_ADALL);    
    con_ADF   = vcf_samples_init_container_array (_FORMAT_ADF);    
    con_ADR   = vcf_samples_init_container_array (_FORMAT_ADR);    
    con_F1R2  = vcf_samples_init_container_array (_FORMAT_F1R2);    
    con_F2R1  = vcf_samples_init_container_array (_FORMAT_F2R1);    
    con_MB    = vcf_samples_init_container_array (_FORMAT_MB);    
    con_SB    = vcf_samples_init_container_array (_FORMAT_SB);    
    con_SAC   = vcf_samples_init_container_array (_FORMAT_SAC);    

    // prepare special snips for the odd elements of SB and MB - (AD minus even item 0) 
    for (unsigned i=0; i < 2; i++) {
        seg_prepare_minus_snip_i (VCF, con_AD.items[i].dict_id, con_SB.items[i*2].dict_id, sb_snip, i);
        seg_prepare_minus_snip_i (VCF, con_AD.items[i].dict_id, con_MB.items[i*2].dict_id, mb_snip, i);
    }

    for (unsigned i=0; i < VCF_MAX_ARRAY_ITEMS/2; i++) 
        seg_prepare_minus_snip_i (VCF, con_AD.items[i].dict_id, con_SAC.items[i*2].dict_id, sac_snip, i);

    for (unsigned i=0; i < VCF_MAX_ARRAY_ITEMS; i++) {
        seg_prepare_minus_snip_i (VCF, con_AD.items[i].dict_id, con_F1R2.items[i].dict_id, f2r1_snip, i);
        seg_prepare_minus_snip_i (VCF, con_AD.items[i].dict_id, con_ADF.items[i].dict_id,  adr_snip,  i);
        seg_prepare_minus_snip_i (VCF, con_AD.items[i].dict_id, con_ADR.items[i].dict_id,  adf_snip,  i);
    }

    seg_prepare_minus_snip (VCF, _FORMAT_RD, _FORMAT_RDR, rdf_snip);
    seg_prepare_minus_snip (VCF, _FORMAT_RD, _FORMAT_RDF, rdr_snip);
    seg_prepare_minus_snip (VCF, _FORMAT_AD, _FORMAT_ADR, adf_snip);
    seg_prepare_minus_snip (VCF, _FORMAT_AD, _FORMAT_ADF, adr_snip);
    seg_prepare_minus_snip (VCF, _FORMAT_DP, _FORMAT_RD,  ad_varscan_snip);
        
    seg_prepare_snip_other (SNIP_COPY, _INFO_AF, 0, 0, af_snip);

    ab_snip_len = sizeof(ab_snip);
    DictId ad_dict_ids[3] = { make_array_item_dict_id(_FORMAT_AD, 0), 
                              make_array_item_dict_id(_FORMAT_AD, 1),
                              (DictId)_FORMAT_AB3 };
    seg_prepare_multi_dict_id_special_snip (VCF_SPECIAL_AB, 3, ad_dict_ids, ab_snip, &ab_snip_len);

    gq_by_gp_len = sizeof (gq_by_gp);
    seg_prepare_multi_dict_id_special_snip (VCF_SPECIAL_GQ, 1, (DictId[]){ (DictId)_FORMAT_GP }, gq_by_gp, &gq_by_gp_len);
    gq_by_gp[gq_by_gp_len++] = '\t';
    
    gq_by_pl_len = sizeof (gq_by_pl);
    seg_prepare_multi_dict_id_special_snip (VCF_SPECIAL_GQ, 1, (DictId[]){ (DictId)_FORMAT_PL }, gq_by_pl, &gq_by_pl_len);
    gq_by_pl[gq_by_pl_len++] = '\t';

    // two redirection snips that must be the same length - we might interchange them in-place in the dictionary
    PL_to_PLy_redirect_snip_len = PL_to_PLn_redirect_snip_len = sizeof (PL_to_PLy_redirect_snip);
    seg_prepare_snip_other (SNIP_REDIRECTION, _FORMAT_PLy, false, 0, PL_to_PLy_redirect_snip);
    seg_prepare_snip_other (SNIP_REDIRECTION, _FORMAT_PLn, false, 0, PL_to_PLn_redirect_snip);
    ASSERT (PL_to_PLy_redirect_snip_len == PL_to_PLn_redirect_snip_len, "Expecting PL_to_PLy_redirect_snip_len%u == PL_to_PLn_redirect_snip_len=%u", PL_to_PLy_redirect_snip_len, PL_to_PLn_redirect_snip_len);

    vcf_samples_zip_initialize_PS_PID();
    vcf_gwas_zip_initialize(); // no harm if not GWAS
}


void vcf_samples_seg_initialize (VBlockVCFP vb)
{
    ctx_set_store (VB, STORE_INT, FORMAT_ADALL, FORMAT_DP, FORMAT_AD, FORMAT_RD, FORMAT_RDR, FORMAT_RDF,
                   FORMAT_ADR, FORMAT_ADF, FORMAT_SDP, DID_EOL);

    ctx_set_ltype (VB, LT_DYN_INT, FORMAT_RD, FORMAT_GQ, FORMAT_RGQ, DID_EOL);
    
    CTX(FORMAT_GT)->no_stons  = true; // we store the GT matrix in local, so cannot accomodate singletons
    
    vb->ht_matrix_ctx = CTX(FORMAT_GT_HT); // different for different data types

    CTX(FORMAT_DP)->flags.same_line = true; // DP value, when delta's against another ctx, is always relative to the other value in the current sample (regardless of whether DP or the other value are reconstructed first)

    // create additional contexts as needed for compressing FORMAT/GT - must be done before merge
    if (vcf_num_samples) 
        codec_pbwt_seg_init (VB, CTX(FORMAT_PBWT_RUNS), CTX(FORMAT_PBWT_FGRC));

    ctx_consolidate_stats (VB, FORMAT_GT, FORMAT_GT_HT, FORMAT_GT_HT_INDEX, FORMAT_PBWT_RUNS, FORMAT_PBWT_FGRC, DID_EOL);
    ctx_consolidate_stats (VB, FORMAT_AB, FORMAT_AB3, DID_EOL);

    // determine which way to seg PL - Mux by dosage or Mux by dosageXDP, or test both options
    CTX(FORMAT_PL)->no_stons = true;
    vb->PL_mux_by_DP = (flag.best && !z_is_dvcf && !segconf.running && segconf.has_DP_before_PL) // only in --best, because it is very slow
        ? segconf.PL_mux_by_DP // set by a previous VB in vcf_FORMAT_PL_decide or still in its initial value of PL_mux_by_DP_TEST
        : PL_mux_by_DP_NO;

    // initialize dosage multiplexers
    #define init_mux_by_dosage(name) seg_mux_init ((VBlockP)vb, CTX(FORMAT_##name), 4, VCF_SPECIAL_MUX_BY_DOSAGE, CTX(FORMAT_##name)->no_stons, (MultiplexerP)&vb->mux_##name, "0123")
    init_mux_by_dosage(PRI);
    init_mux_by_dosage(GL);
    init_mux_by_dosage(DS);
    init_mux_by_dosage(PP);
    init_mux_by_dosage(GP);
    init_mux_by_dosage(PVAL);
    init_mux_by_dosage(FREQ);
    init_mux_by_dosage(RD);
    init_mux_by_dosage(PLn);

    seg_mux_init (VB, CTX(FORMAT_PLy), MUX_CAPACITY(vb->mux_PLy), VCF_SPECIAL_MUX_BY_DOSAGExDP, false, (MultiplexerP)&vb->mux_PLy, NULL);
    
    if (segconf.has[FORMAT_DP]) {
        seg_mux_init (VB, CTX(FORMAT_GQ),  MUX_CAPACITY(vb->mux_GQ),  VCF_SPECIAL_MUX_BY_DOSAGExDP, false, (MultiplexerP)&vb->mux_GQ, NULL);
        seg_mux_init (VB, CTX(FORMAT_RGQ), MUX_CAPACITY(vb->mux_RGQ), VCF_SPECIAL_RGQ, false, (MultiplexerP)&vb->mux_RGQ, NULL);
    }
    else
        init_mux_by_dosage(GQ);

    // flags to send to PIZ
    vb->flags.vcf.use_null_DP_method = segconf.use_null_DP_method;
}

void vcf_samples_seg_finalize (VBlockVCFP vb)
{
    if (segconf.running) {
        // In case of dependency DAG: DP->(sum)AD->(mux)GT we can't have GT->(null)DP
        if (segconf.FORMAT_DP_method == by_AD) segconf.use_null_DP_method = false;

        // percent of (samples x lines) that have a dosage value of 0,1 or 2 
        segconf.pc_has_dosage = (float)segconf.count_dosage[1] / (float)(segconf.count_dosage[0] + segconf.count_dosage[1]);

        // whether we should seg GQ as a function of GP or PL (the better of the two) - only if this works for at least 20% of the samples
        segconf.GQ_by_GP = (segconf.count_GQ_by_GP > vb->lines.len * vcf_num_samples / 5) && (segconf.count_GQ_by_GP >  segconf.count_GQ_by_PL);
        segconf.GQ_by_PL = (segconf.count_GQ_by_PL > vb->lines.len * vcf_num_samples / 5) && (segconf.count_GQ_by_PL >= segconf.count_GQ_by_GP);
    }
    else {
        vcf_samples_seg_finalize_PS_PID(vb);
    }
}

// returns true is the sample has a '.' value
static inline bool vcf_seg_sample_has_null_value (uint64_t dnum, ContextP *ctxs, STRps(sf))
{
    for (int i=1; i < n_sfs; i++) // start from 1 - we know 0 is GT
        if (ctxs[i]->dict_id.num == dnum)
            return sf_lens[i]==1 && sfs[i][0]=='.'; // null DP found on this line - return true if its not '.' or empty

    return false; // no DP at all on this line
}

//--------------------
// Multiplex by dosage
// -------------------

// returns: 0,1,2 correspond to 0/0 0/1 1/1, 3 means "no dosage" - multi-allelic or '.' or ploidy > 2
// -1 means caller should not use her SPECIAL
int vcf_seg_get_mux_channel_i (VBlockVCFP vb, bool fail_if_dvcf_refalt_switch)
{
    // fail if this is a DVCF ref<>alt switch, if caller requested so,
    // or if there is no valid GT in this sample
    if (fail_if_dvcf_refalt_switch && z_is_dvcf)// && LO_IS_OK_SWITCH (last_ostatus)) // TODO: no real need to fallback unless REF⇆ALT switch, but this doesn't work yet
        return -1;

    // fail if there is no GT in this variant
    if (!ctx_encountered_in_line (VB, FORMAT_GT))
        return -1;

    int64_t dosage = CTX(FORMAT_GT)->last_value.i; // dosage stored here by vcf_seg_FORMAT_GT
     
    return (dosage >= 0 && dosage <= 2) ? dosage : 3; // 3 happens if sample has ploidy > 2, or if one of the alleles is not 0 or 1
}

// if cell is NULL, leaves it up to the caller to seg to the channel 
ContextP vcf_seg_FORMAT_mux_by_dosage (VBlockVCFP vb, ContextP ctx, STRp(cell), const DosageMultiplexer *mux) 
{
    int channel_i = vcf_seg_get_mux_channel_i (vb, true);

    // we don't use the multiplexer if its a DVCF REF⇆ALT switch variant as GT changes
    if (channel_i == -1) {
        if (cell) seg_by_ctx (VB, STRa(cell), ctx, cell_len);
        return ctx;
    }

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, ctx->did_i, MUX, channel_i);

    if (cell) {
        if (channel_ctx->ltype == LT_DYN_INT)
            seg_integer_or_not (VB, channel_ctx, STRa(cell), cell_len);
        else
            seg_by_ctx (VB, STRa(cell), channel_ctx, cell_len);
    }

    // note: this is not necessarily all-the-same - there could be unmuxed snips due to REF⇆ALT switch, and/or WORD_INDEX_MISSING 
    seg_by_ctx (VB, STRa(mux->snip), ctx, 0);

    return channel_ctx;
}

// mirroring seg, we accept monoploid or diploid genotypes, with alleles 0 and 1.
// there return 0/0,0->0 ; 0/1,1/0,1->1 ; 1/1->2 ; others->3
int vcf_piz_get_mux_channel_i (VBlockP vb)
{
    STRlast (gt, FORMAT_GT);

    if (gt_len == 3) { // diploid
        if ((gt[0]!='0' && gt[0]!='1') || (gt[2]!='0' && gt[2]!='1')) return 3; 
        return (int)gt[0] + (int)gt[2] - 2*'0';
    }
    else if (gt_len == 1) { // monoploid
        if (gt[0]!='0' && gt[0]!='1') return 3; 
        return (int)gt[0]- '0';
    }
    else
        return 3;
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_DOSAGE)
{
    int channel_i = vcf_piz_get_mux_channel_i (vb);
    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

// if cell is NULL, leaves it up to the caller to seg to the channel 
static void vcf_seg_FORMAT_mux_by_dosagexDP (VBlockVCFP vb, ContextP ctx, STRp(cell), void *mux_p) 
{
    ConstMultiplexerP mux = (ConstMultiplexerP)mux_p;

    if (!ctx_encountered (VB, FORMAT_DP)) goto cannot_use_special; // no DP in the FORMAT of this line

    int64_t DP;
    if (!str_get_int (STRlst (FORMAT_DP), &DP)) // In some files, DP may be '.'
        DP=0;

    int channel_i = vcf_seg_get_mux_channel_i (vb, true); // we don't use the multiplexer if its a DVCF REF⇆ALT switch variant as GT changes
    if (channel_i == -1) goto cannot_use_special;

    unsigned num_dps = mux->num_channels / 3;
    DP = MAX_(0, MIN_(DP, num_dps-1));
    channel_i = (channel_i == 3) ? (num_dps * 3) : (DP*3 + channel_i);

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, ctx->did_i, MUX, channel_i);

    if (channel_ctx->ltype == LT_DYN_INT)
        seg_integer_or_not (VB, channel_ctx, STRa(cell), cell_len);
    else
        seg_by_ctx (VB, STRa(cell), channel_ctx, cell_len);

    // note: this is not necessarily all-the-same - there could be unmuxed snips due to REF⇆ALT switch, and/or WORD_INDEX_MISSING 
    seg_by_ctx (VB, MUX_SNIP(mux), MUX_SNIP_LEN(mux), ctx, 0);
    return;

cannot_use_special:
    seg_by_ctx (VB, STRa(cell), ctx, cell_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_DOSAGExDP)
{
    unsigned num_channels = ctx->ctx_cache.len32 ? ctx->ctx_cache.len32 : (1 + str_count_char (STRa(snip), '\t'));
    unsigned num_dps = num_channels / 3;

    rom DP_str;
    int64_t DP = reconstruct_peek (vb, CTX(FORMAT_DP), &DP_str, 0).i;
    DP = (*DP_str=='.') ? 0 : MAX_(0, MIN_(DP, num_dps-1));

    int channel_i = vcf_piz_get_mux_channel_i (vb); 
    channel_i = (channel_i == 3) ? (num_dps * 3) : (DP*3 + channel_i);

    ContextP channel_ctx = MCTX (channel_i, snip, snip_len);
    ASSPIZ (channel_ctx, "Cannot find channel context of channel_i=%d of multiplexed context %s", channel_i, ctx->tag_name);

    reconstruct_from_ctx (vb, channel_ctx->did_i, 0, reconstruct);

    if (ctx->flags.store == STORE_NONE) return NO_NEW_VALUE;

    // propagate last_value up
    new_value->i = channel_ctx->last_value.i; // note: last_value is a union, this copies the entire union
    return HAS_NEW_VALUE; 
}

// used when CTX is expected to be (BaseCtx-MinusCtx) - if it indeed is, we use a special snip
static WordIndex vcf_seg_FORMAT_minus (VBlockVCFP vb, ContextP ctx, 
                                       STRp(str), int64_t value, // option 1,2
                                       ContextP base_ctx, ContextP minus_ctx, STRp(minus_snip))
{
    // we can use the formula only if AD,F1R1 were encountered in this line, and that they have the number of items as us
    if (str && !str_get_int (STRa(str), &value)) goto fallback;

    ctx_set_last_value (VB, ctx, value);

    bool use_formula = ctx_has_value (VB, base_ctx->did_i) && ctx_has_value (VB, minus_ctx->did_i) &&
                       value == base_ctx->last_value.i - minus_ctx->last_value.i;

    // case: formula works - seg as minus
    if (use_formula)
        return seg_by_ctx (VB, STRa(minus_snip), ctx, str_len);

    // case: the formula doesn't work for this item - seg a normal snip
    else
        fallback:
        return seg_by_ctx (VB, STRa(str), ctx, str_len);
}

// used for DP, GQ, A0D and otheres - store in transposed matrix in local 
static inline void vcf_seg_FORMAT_transposed (VBlockVCFP vb, ContextP ctx, STRp(cell), unsigned add_bytes)
{
    ctx->ltype = LT_UINT32_TR;
    ctx->flags.store = STORE_INT;
    
    buf_alloc (vb, &ctx->local, 1, vb->lines.len * vcf_num_samples, uint32_t, 1, CTX_TAG_LOCAL);

    if (str_is_1char (cell, '.')) 
        BNXT32 (ctx->local) = 0xffffffff;
    
    else {
        ASSSEG (str_get_int (STRa(cell), &ctx->last_value.i) && ctx->last_value.i >= 0 && ctx->last_value.i <= 0xfffffffe, 
                "While compressing %s expecting an integer in the range [0, 0xfffffffe] or a '.', but found: %.*s", 
                ctx->tag_name, cell_len, cell);

        BNXT32 (ctx->local) = (uint32_t)ctx->last_value.i;
    }

    // add a LOOKUP to b250
    seg_by_ctx (VB, (char []){ SNIP_LOOKUP }, 1, ctx, add_bytes);
}

// a comma-separated array - each element goes into its own item context, single repeat
static WordIndex vcf_seg_FORMAT_A_R (VBlockVCFP vb, ContextP ctx, SmallContainer con /* by value */, STRp(value), StoreType item_store_type,
                                     void (*seg_item_cb)(VBlockVCFP, ContextP ctx, unsigned n_items, const char**, const uint32_t*, ContextP *item_ctxs, const int64_t*))
{   
    str_split (value, value_len, VCF_MAX_ARRAY_ITEMS, ',', item, false);
    
    if (!(con.nitems_lo = n_items)) 
        return seg_by_ctx (VB, value, value_len, ctx, value_len); // too many items - normal seg

    ContextP item_ctxs[con.nitems_lo];
    int64_t values[con.nitems_lo];

    for (unsigned i=0; i < con.nitems_lo; i++) {

        item_ctxs[i] = ctx_get_ctx (vb, con.items[i].dict_id);
        
        if (!item_ctxs[i]->is_initialized) {
            item_ctxs[i]->flags.store = item_store_type;
            item_ctxs[i]->is_initialized = true;
            ctx_consolidate_stats (VB, ctx->did_i, item_ctxs[i]->did_i, DID_EOL);
        }

        if (seg_item_cb) {
            if (str_get_int (STRi(item, i), &values[i])) 
                ctx_set_last_value (VB, item_ctxs[i], values[i]);
            else
                seg_item_cb = NULL; // can't use callback if not all items are int
        }
    }

    // case: seg items via callback
    if (seg_item_cb)
        seg_item_cb (vb, ctx, con.nitems_lo, items, item_lens, item_ctxs, values);

    // case: seg items as normal snips
    else 
        for (unsigned i=0; i < con.nitems_lo; i++) 
            seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);

    ctx->last_txt.len = con.nitems_lo; // seg only: for use by vcf_seg_*_items callbacks
    
    return container_seg (vb, ctx, (ContainerP)&con, 0, 0, con.nitems_lo-1); // account for the commas
}

//----------
// FORMAT/AD
// ---------

// <ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
static inline WordIndex vcf_seg_FORMAT_AD_varscan (VBlockVCFP vb, ContextP ctx, STRp(ad_str))
{
    // case: AD = DP-RD
    int64_t ad;
    if (ctx_has_value (VB, FORMAT_DP) &&
        ctx_has_value (VB, FORMAT_RD) && 
        str_get_int (STRa(ad_str), &ad) &&
        ad == CTX(FORMAT_DP)->last_value.i - CTX(FORMAT_RD)->last_value.i) 
    
        return vcf_seg_FORMAT_minus (vb, ctx, 0, ad_str_len, ad, CTX(FORMAT_DP), CTX(FORMAT_RD), STRa(ad_varscan_snip));

    // case: we have only one sample, and INFO/ADP - we expect FORMAT/AD and INFO/ADP to be related
    else if (ctx_has_value_in_line_(vb, CTX(INFO_ADP)) && vcf_num_samples==1)
        return seg_delta_vs_other (VB, ctx, CTX(INFO_ADP), STRa(ad_str));

    else
        return seg_by_ctx (VB, STRa(ad_str), ctx, ad_str_len);
}

// <ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
static void vcf_seg_AD_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values)
{       
    bool has_adall_this_sample = segconf.has[FORMAT_ADALL] && ctx_encountered (VB, FORMAT_ADALL); // note: we can't delta vs ADALL unless segconf says so, bc it can ruin other fields that rely on peeking AD, eg AB
    int64_t sum = 0; 

    for (unsigned i=0; i < n_items; i++) {

        // If we have ADALL in this file, we delta vs ADALL if we have it in this sample, or seg normally if not
        if (segconf.has[FORMAT_ADALL]) {
            // case: we had ADALL preceeding in this sample, seg as delta vs. ADALL 
            if (has_adall_this_sample)
                seg_delta_vs_other (VB, item_ctxs[i], ECTX (con_ADALL.items[i].dict_id), NULL, item_lens[i]);
            else 
                seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
        }

        // case: Since item 0 (depth of REF) is usually somewhat related to the overall sample depth,
        // and hence values within a sample are expected to be correlated - we store it transposed, and the other items - normally
        else {
            if (i==0 && vcf_num_samples > 1)
                vcf_seg_FORMAT_transposed (vb, item_ctxs[0], STRi(item, 0), item_lens[0]);
            
            else if (i==0 || i==1) {
                if (!vb->mux_AD[i].num_channels)
                    seg_mux_init (VB, item_ctxs[i], 4, VCF_SPECIAL_MUX_BY_DOSAGE, false, (MultiplexerP)&vb->mux_AD[i], "0123");
                
                vcf_seg_FORMAT_mux_by_dosage (vb, item_ctxs[i], STRi(item, i), &vb->mux_AD[i]);
            }
            else
                seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
        }
        sum += values[i];
    }

    // AD value is sum of its items
    ctx_set_last_value (VB, ctx, sum);

    memcpy (vb->ad_values, values, n_items * sizeof (values[0]));
}

//----------------------
// FORMAT/F2R1, ADF, ADR
//----------------------

// used when Vector is expected to be (AD-OtherVector) - if it indeed is, we use a special snip
static void vcf_seg_AD_complement_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values,
                                         DictId other_dict_id, const SmallContainer *other_con,
                                         char my_snips[][32], unsigned *my_snip_lens)
{
    // we can use the formula only if AD,F1R1 were encountered in this line, and that they have the number of items as us
    ContextP ad_ctx=CTX(FORMAT_AD), other_ctx;
    bool use_formula = ctx_encountered (VB, FORMAT_AD) &&
                       ctx_encountered_by_dict_id (VB, other_dict_id, &other_ctx) &&
                       ad_ctx->last_txt.len    == n_items &&  // last_txt_len is # of items stored by vcf_seg_FORMAT_A_R 
                       other_ctx->last_txt.len == n_items;

    for (unsigned i=0; i < n_items; i++) {

        // case: as expected, F1R2 + F2R1 = AD - seg as a F2R1 as a MINUS snip
        if (use_formula && vb->ad_values[i] == values[i] + ECTX (other_con->items[i].dict_id)->last_value.i) 
            seg_by_ctx (VB, my_snips[i], my_snip_lens[i], item_ctxs[i], item_lens[i]); 

        // case: the formula doesn't work for this item - seg a normal snip
        else
            seg_by_ctx (VB, STRi(item,i), item_ctxs[i], item_lens[i]);
    }
}

// F2R1 = AD - F1R2 (applied if AD and F1R2 are encountered before F2R1)
static void vcf_seg_F2R1_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values)
{
    vcf_seg_AD_complement_items (vb, ctx, STRas(item), item_ctxs, values, _FORMAT_F1R2, &con_F1R2, f2r1_snips, f2r1_snip_lens);
}

// ADF = AD - ADR (applied if AD and ADR are encountered before ADF)
static void vcf_seg_ADF_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values)
{
    vcf_seg_AD_complement_items (vb, ctx, STRas(item), item_ctxs, values, _FORMAT_ADR, &con_ADR, adf_snips, adf_snip_lens);
}

// ADR = AD - ADF (applied if AD and ADF are encountered before ADR)
static void vcf_seg_ADR_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values)
{
    vcf_seg_AD_complement_items (vb, ctx, STRas(item), item_ctxs, values, _FORMAT_ADF, &con_ADF, adr_snips, adr_snip_lens);
}

//----------
// FORMAT/SB
//----------

// For bi-allelic SNPs, sum every of two values is expected to equal the corresponding value in AD. Example: AD=59,28 SB=34,25,17,11. 
// seg the second of every pair as a MINUS snip
static void vcf_seg_SB_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values)
{
    // verify that AD was encountered in this line, and that it has exactly half the number of items as us
    ContextP ad_ctx=CTX(FORMAT_AD);
    bool use_formula = ctx_encountered (VB, FORMAT_AD) && ad_ctx->last_txt.len == 2 && n_items == 4; // note: last_txt_len = # of items stored by vcf_seg_FORMAT_A_R

    for (unsigned i=0; i < n_items; i++) {

        // seg odd-numbered element as AD - (even element), if the sum is correct
        if (use_formula && i%2 && vb->ad_values[i/2] == values[i-1] + values[i]) 
            seg_by_ctx (VB, sb_snips[i/2], sb_snip_lens[i/2], item_ctxs[i], item_lens[i]); 

        else {
            item_ctxs[i]->flags.store = STORE_INT; // consumed by the odd items ^
            seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
        }
    }
}

//-----------
// FORMAT/SAC
//-----------

// sum every of two values is expected to equal the corresponding value in AD. seg the second of every pair as a MINUS snip
static void vcf_seg_SAC_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values)
{
    // verify that AD was encountered in this line, and that it has exactly half the number of items as us
    bool use_formula = ctx_encountered (VB, FORMAT_AD) && 2 * CTX(FORMAT_AD)->last_txt.len == n_items; // note: last_txt_len = # of items stored by vcf_seg_FORMAT_A_R

    for (unsigned i=0; i < n_items; i++) {

        // seg odd-numbered element as AD - (even element), if the sum is correct
        if (use_formula && i%2 && vb->ad_values[i/2] == values[i-1] + values[i]) 
            seg_by_ctx (VB, sac_snips[i/2], sac_snip_lens[i/2], item_ctxs[i], item_lens[i]); 

        else {
            item_ctxs[i]->flags.store = STORE_INT; // consumed by the odd items ^
            seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
        }
    }
}

//--------------------------
// FORMAT/ICNT (DRAGEN gVCF)
//--------------------------
static void vcf_seg_ICNT (VBlockVCFP vb, ContextP ctx, STRp(ICNT))
{
    if (ctx_encountered (VB, FORMAT_AD)) {        
        str_split_ints (ICNT, ICNT_len, 0, ',', icnt, false);
        if (!n_icnts) goto fallback;

        uint8_t snip[2 + n_icnts];
        snip[0] = SNIP_SPECIAL;
        snip[1] = VCF_SPECIAL_ICNT;

        int64_t delta64 = vb->ad_values[0] - icnts[0];
        if (delta64 < -112 || delta64 > 111) goto fallback;    // map the range [-112,111]->[32,255]
        snip[2] = 144 + delta64;

        for (int i=1; i < n_icnts; i++) {
            if (icnts[i] < 0 || icnts[i] > 223) goto fallback; // map the range [0,223]->[32,255]
            snip[i+2] = icnts[i] + 32;
        }

        seg_by_ctx (VB, (char *)snip, 2 + n_icnts, ctx, ICNT_len);

        return;
    }

fallback:
    seg_by_ctx (VB, STRa(ICNT), ctx, ICNT_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_ICNT)
{
    rom ad0_str = last_txt (vb, FORMAT_AD);
    int64_t ad0 = atoll (ad0_str);

    int64_t delta = (int64_t)(uint8_t)snip[0] - 144;
    RECONSTRUCT_INT (ad0 - delta);

    for (int i=1; i < snip_len; i++) {
        RECONSTRUCT1 (',');
        RECONSTRUCT_INT ((int64_t)(uint8_t)snip[i] - 32);
    }

    return NO_NEW_VALUE;
}

// --------------------------
// FORMAT/SPL (DRAGEN gVCF)
// --------------------------
static void vcf_seg_SPL (VBlockVCFP vb, ContextP ctx, STRp(SPL))
{
    if (ctx_encountered (VB, FORMAT_PL)) {
        STRlast (PL, FORMAT_PL);

        str_split_ints (SPL, SPL_len, 0, ',', spl, false);
        str_split_ints (PL,  PL_len,  0, ',', pl , false);

        if (!n_pls || n_spls != n_pls) goto fallback;

        // verify prediction
        bool has_delta = false;
        uint8_t snip[2 + n_pls];
        snip[0] = SNIP_SPECIAL;
        snip[1] = VCF_SPECIAL_SPL;

        for (int i=0; i < n_pls; i++) {
            int delta64 = spls[i] - MIN_(pls[i], 255);
            if (delta64 < -112 || delta64 > 111) goto fallback; // doesn't fit 

            snip[i+2] = delta64 + 144; // map [-112,111] -> [32,255]
            if (delta64) has_delta = true;
        }

        seg_by_ctx (VB, (char *)snip, 2 + (has_delta ? n_pls : 0), ctx, SPL_len);
        return;
    }

fallback:
    seg_by_ctx (VB, STRa(SPL), ctx, SPL_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_SPL)
{
    STRlast (PL, FORMAT_PL);
    str_split_ints (PL, PL_len, snip_len/*0 or correct number*/, ',', pl,  false);

    for (int i=0; i < n_pls; i++) {
        int64_t delta = (snip_len ? ((int64_t)(uint8_t)snip[i] - 144) : 0);
        RECONSTRUCT_INT (MIN_(pls[i], 255) + delta);

        if (i < n_pls-1) 
            RECONSTRUCT1 (',');
    }

    return NO_NEW_VALUE;
}

//----------
// FORMAT/MB
//----------

// For bi-allelic SNPs: sum every of two items is expected to equal the corresponding value in AD. Example: AD=7,49 F2R1=3,28 MB=4,3,26,23 
// In addition, the even-numbered item is quite similar to the corresponding value in F2R1.
// Seg the even items as delta from F2R1 and odd items as a MINUS snip between AD and the preceding even item
static void vcf_seg_MB_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values)
{
    bool use_formula_even = ctx_encountered (VB, FORMAT_F2R1) && CTX(FORMAT_F2R1)->last_txt.len == 2 && n_items == 4;
    bool use_formula_odd  = ctx_encountered (VB, FORMAT_AD)   && CTX(FORMAT_AD)  ->last_txt.len == 2 && n_items == 4; // last_txt_len is # of items set by vcf_seg_FORMAT_A_R

    for (unsigned i=0; i < n_items; i++) {

        // if possible, seg even-numbered element delta vs the corresponding element in F2R1
        if (use_formula_even && !(i%2)) { 
            seg_delta_vs_other (VB, item_ctxs[i], ECTX (con_F2R1.items[i/2].dict_id), NULL, item_lens[i]);
            item_ctxs[i]->flags.store = STORE_INT; // consumed by the odd items (below)
        }

        // if possible, seg odd-numbered element as AD minus (even element), if the sum is correct
        else if (use_formula_odd && i%2 && vb->ad_values[i/2] == values[i-1] + values[i]) 
            seg_by_ctx (VB, mb_snips[i/2], mb_snip_lens[i/2], item_ctxs[i], item_lens[i]); 
        
        else { // fallback if formulas don't work
            seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
            item_ctxs[i]->flags.store = STORE_INT; // possibly consumed by the odd items (^)
        }
    }
}

//----------
// FORMAT/AF
// ---------

static inline WordIndex vcf_seg_FORMAT_AF (VBlockVCFP vb, ContextP ctx, STRp(cell))
{
    if (vcf_num_samples == 1 && // very little hope that INFO/AF is equal to FORMAT/AF if we have more than one sample
        !z_is_dvcf &&       // note: we can't use SNIP_COPY in dual coordinates, because when translating, it will translate the already-translated INFO/AF
        ctx_encountered_in_line (VB, INFO_AF) && 
        str_issame (cell, CTX(INFO_AF)->last_snip))
        return seg_by_ctx (VB, af_snip, af_snip_len, ctx, cell_len);
    else
        return vcf_seg_FORMAT_A_R (vb, ctx, con_AF, STRa(cell), STORE_NONE, NULL);
}

//----------
// FORMAT/GQ
// ---------

static SORTER (value_sorter)
{
    return DESCENDING_RAW (*(int64_t *)a, *(int64_t*)b); // sort in reverse order - usually faster as GP[0] / PL[0] are usually the biggest (corresponding to GT=0/0)
}

static int64_t vcf_predict_GQ (VBlockVCFP vb, Did src_did_i)
{
    bool is_gp = (src_did_i == FORMAT_GP);
    STR(src);

    if (IS_ZIP)
        SETlast (src, src_did_i);
    else
        reconstruct_peek (VB, CTX(src_did_i), pSTRa(src));

    str_split (src, src_len, 30, ',', item, false);

    int64_t values[n_items];
    unsigned n_values=0;

    for (int i=0; i < n_items; i++) {
        if (item_lens[i]==1 && items[i][0]=='.')
            values[n_values++] = 0; // consider a '.' to be 0

        else if (is_gp) {
            double f;
            if (str_get_float (STRi(item,i), &f, 0, 0)) 
                values[n_values++] = (int64_t)(f+0.5); // round to nearest integer
        }
        else
            n_values += str_get_int (STRi(item, i), &values[n_values]); // increment if successfully read an int
    }

    if (!n_values) return 0; // array to long (n_item=0) or none is an integer


    // now we have an array of integers that is the same or shorter that the GP/PL array. Now we sort it.
    qsort (values, n_values, sizeof(values[0]), value_sorter);

    int64_t mid_value = values[n_values/2];

    if (!is_gp) mid_value = MIN_(mid_value, 99); // for PL, it is capped by 99

    return mid_value;
}

static void vcf_seg_FORMAT_GQ_do_seg (VBlockVCFP vb, int64_t gq_value, int64_t prediction, STRp(template_snip), unsigned gq_len)
{
    STRl(snip, 50) = template_snip_len;
    memcpy (snip, template_snip, template_snip_len);
    snip_len += str_int (prediction - gq_value, &snip[snip_len]);

    seg_by_did (VB, STRa (snip), FORMAT_GQ, gq_len);
}

static inline void vcf_seg_FORMAT_GQ (VBlockVCFP vb)
{
    ContextP gq_ctx = CTX(FORMAT_GQ);
    STRlast (gq, FORMAT_GQ);

    // second best: mux by dosage, and possibly DP
    if (!segconf.running && !segconf.GQ_by_GP && !segconf.GQ_by_PL) {
        if (segconf.has[FORMAT_DP]) 
            vcf_seg_FORMAT_mux_by_dosagexDP (vb, gq_ctx, STRa(gq), &vb->mux_GQ);
        else                        
            vcf_seg_FORMAT_mux_by_dosage (vb, gq_ctx, STRa(gq), (DosageMultiplexer *)&vb->mux_GQ);
        return;
    }

    int64_t gq_value;
    if (!str_get_int (STRa(gq), &gq_value)) goto fallback;

    if (segconf.running) {
        if (ctx_encountered (VB, FORMAT_GP) && vcf_predict_GQ (vb, FORMAT_GP) == gq_value) segconf.count_GQ_by_GP++;
        if (ctx_encountered (VB, FORMAT_PL) && vcf_predict_GQ (vb, FORMAT_PL) == gq_value) segconf.count_GQ_by_PL++;
    }
    
    else {
        if (segconf.GQ_by_GP && ctx_encountered (VB, FORMAT_GP)) {
            int64_t prediction = vcf_predict_GQ (vb, FORMAT_GP);
            vcf_seg_FORMAT_GQ_do_seg (vb, gq_value, prediction, STRa(gq_by_gp), gq_len);
            return;
        }

        if (segconf.GQ_by_PL && ctx_encountered (VB, FORMAT_PL)) {
            int64_t prediction = vcf_predict_GQ (vb, FORMAT_PL);
            vcf_seg_FORMAT_GQ_do_seg (vb, gq_value, prediction, STRa(gq_by_pl), gq_len);
            return;
        }
    }

fallback:
    // fallback if not match to GP or PL, or segconf.running
    seg_by_ctx (VB, STRa(gq), gq_ctx, gq_len);    
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_FORMAT_GQ)
{
    rom tab = memchr (snip, '\t', snip_len);

    ContextP src_ctx = SCTX0 (snip);

    int64_t prediction = vcf_predict_GQ (VB_VCF, src_ctx->did_i);
    int64_t delta = atoi (tab+1);

    new_value->i = prediction - delta;
    RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

//-----------
// FORMAT/RGQ
// ----------

// <ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
// Appears in GVCF in lines which are no variants (i.e. no ALT)
static inline void vcf_seg_FORMAT_RGQ (VBlockVCFP vb, ContextP ctx, STRp(rgq), ContextP gt_ctx, STRp(gt))
{
    ConstMultiplexerP mux = (ConstMultiplexerP)&vb->mux_RGQ;
        
    // prediction: we have GT, and if GT[0]=. then RGQ=0. Fallback seg in case prediction fail
    if (gt_ctx->did_i != FORMAT_GT ||                    // prediction failed: first subfield isn't GT
        (gt[0] == '.' && (rgq_len != 1 || *rgq != '0'))) // prediction failed: GT[0]=. and yet RGQ!="0"
        goto fallback;

    // case: GT[0] is not '.' - seg the value of RGQ multiplexed by DP
    if (gt[0] != '.') {
        if (!segconf.has[FORMAT_DP]          ||    // segconf didn't detect FORMAT/DP so we didn't initialize the mux
            !ctx_encountered (VB, FORMAT_DP) ||    // no DP in the FORMAT of this line
            segconf.running) goto fallback;        // multiplexor not initalized yet 

        int64_t DP;
        if (!str_get_int (STRlst(FORMAT_DP), &DP)) // in some files, DP may be '.'
            DP=0;

        int channel_i = MAX_(0, MIN_(DP, mux->num_channels-1));
        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, FORMAT_RGQ, MUX, channel_i);

        seg_integer_or_not (VB, channel_ctx, STRa(rgq), rgq_len);
    }

    seg_by_ctx (VB, MUX_SNIP(mux), MUX_SNIP_LEN(mux), ctx, (gt[0] == '.' ? rgq_len : 0));
    return;

fallback:
    seg_integer_or_not (VB, ctx, STRa(rgq), rgq_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_RGQ)
{
    rom gt = last_txt(vb, FORMAT_GT);
    
    if (gt[0] == '.') {
        if (reconstruct) RECONSTRUCT1 ('0');
        new_value->i = 0;
    }
    
    // gt[0] != '.' - demulitplex by FORMAT_DP
    else {
        unsigned num_channels = ctx->ctx_cache.len32 ? ctx->ctx_cache.len32 : (1 + str_count_char (STRa(snip), '\t'));

        rom DP_str;
        int64_t DP = reconstruct_peek (vb, CTX(FORMAT_DP), &DP_str, 0).i;
        int channel_i = (*DP_str=='.') ? 0 : MAX_(0, MIN_(DP, num_channels-1));

        ContextP channel_ctx = MCTX (channel_i, snip, snip_len);
        ASSPIZ (channel_ctx, "Cannot find channel context of channel_i=%d of multiplexed context %s", channel_i, ctx->tag_name);

        reconstruct_from_ctx (vb, channel_ctx->did_i, 0, reconstruct);

        // propagate last_value up
        new_value->i = channel_ctx->last_value.i; // note: last_value is a union, this copies the entire union
    }

    return HAS_NEW_VALUE; 
}

//----------
// FORMAT/DP
// ---------

static inline void vcf_seg_FORMAT_DP (VBlockVCFP vb, ContextP ctx, STRp(cell))
{
    int64_t value;
    bool is_null = (cell_len==1 && *cell=='.');
    bool has_value = !is_null && str_get_int (STRa(cell), &value);

    if (segconf.running) {
        if (is_null) 
            segconf.use_null_DP_method = true;
        
        if (!segconf.FORMAT_DP_method) {
            if (ctx_has_value (VB, FORMAT_AD)) 
                segconf.FORMAT_DP_method = by_AD;

            else if (ctx_has_value (VB, FORMAT_SDP)) 
                segconf.FORMAT_DP_method = by_SDP;
        }
    }

    seg_set_last_txt (VB, ctx, STRa(cell));

    // case - we have FORMAT/AD - calculate delta vs the sum of AD components
    if (segconf.FORMAT_DP_method == by_AD && ctx_has_value (VB, FORMAT_AD))
        seg_delta_vs_other (VB, ctx, CTX(FORMAT_AD), STRa(cell));

    // case - we have FORMAT/SDP - calculate delta vs the sum of AD components
    else if (segconf.FORMAT_DP_method == by_SDP && ctx_has_value (VB, FORMAT_SDP))
        seg_delta_vs_other (VB, ctx, CTX(FORMAT_SDP), STRa(cell));
    
    // single-sample default: seg against INFO/DP
    else if (vcf_num_samples == 1) {
        bool info_dp_is_int = ctx_has_value_in_line_(VB, CTX(INFO_DP));

        // special means: if FORMAT/DP>=1, INFO/DP is an integer, if FORMAT/DP==0, INFO/DP has no integer value
        if (has_value && ((value > 0) == info_dp_is_int)) {
            SNIPi2 (SNIP_SPECIAL, VCF_SPECIAL_DP_by_DP_single, info_dp_is_int ? (value - CTX(INFO_DP)->last_value.i) : 0);
            seg_by_ctx (VB, STRa(snip), ctx, cell_len);
        }
        else // value is not an integer or '.', or INFO and FORMAT don't agree on existance of an integer DP
            seg_by_ctx (VB, STRa(cell), ctx, cell_len);
    }

    // multi-sample default: store in transposed matrix (or just LOOKUP from local if not transposable)
    else 
        vcf_seg_FORMAT_transposed (vb, ctx, STRa(cell), cell_len); // this handles DP that is an integer or '.'

    if (*cell != '.' && CTX(INFO_DP)->dp.by_format_dp) 
        CTX(INFO_DP)->dp.sum_format_dp += ctx->last_value.i; // we may delta INFO/DP against this sum

    // add up DP's in certain conditions, for consumption by INFO/QD predictor
    if (has_value)
        vcf_seg_sum_DP_for_QD (vb, value);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_DP_by_DP_single)
{
    int64_t info_dp = ctx_has_value_in_line_(vb, CTX(INFO_DP)) ? CTX(INFO_DP)->last_value.i : 0;

    new_value->i = atoi (snip) + info_dp;
    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    return HAS_NEW_VALUE;
}

//-----------
// FORMAT/PGT
// ----------

// by default, a ploidy-2 PGT would have 5 values: (0|0, 0|1, 1|0, 1|1, .) this method usually reduces the dictionary
// to 2 values: SPECIAL and 1|0.
static inline void vcf_seg_FORMAT_PGT (VBlockVCFP vb, ContextP ctx, STRp(pgt), ContextP *ctxs, STRps(sf))
{
    // case: PGT=. and PID=.
    if ((pgt_len==1 && *pgt=='.' && vcf_seg_sample_has_null_value (_FORMAT_PID, ctxs, STRas(sf)))

    // case: haplotypes of PGT are the same and in the same order as GT
    || (ctxs[0]->did_i == FORMAT_GT && sf_lens[0] == 3 && pgt_len==3 && 
        pgt[0]==sfs[0][0] && pgt[2]==sfs[0][2] && pgt[1] == '|'))

        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_PGT }), 2, ctx, pgt_len); 

    // other cases - normal seg: haplotypes appear in reverse order, . but not PID, ploidy!=2, non-standard PGT or GT format etc
    else
        seg_by_ctx (VB, STRa(pgt), ctx, pgt_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_PGT)
{
    if (reconstruct) {
        STR(pid);
        reconstruct_peek (vb, CTX(FORMAT_PID), pSTRa(pid)); // note: we can't use last_txt, because PS might be reconstructed before PID, as its peeked by GT
        
        // case: if this SPECIAL was used with PID='.'
        if (str_is_1char (pid, '.')) 
            RECONSTRUCT1 ('.');

        // case: ht0 and ht1 are the same as in GT
        else {
            rom gt = last_txt(vb, FORMAT_GT);
            RECONSTRUCT1 (gt[0]);
            RECONSTRUCT1 ('|');
            RECONSTRUCT1 (gt[2]);
        }
    }

    return NO_NEW_VALUE;
}

// ---------------------
// INFO/AF and FORMAT/AF
// ---------------------

// translate to (max_value - value).
static int32_t vcf_piz_luft_trans_complement_to_max_value (VBlockP vb, ContextP ctx, char *recon, int32_t recon_len, bool validate_only, double max_value)
{
    if (IS_TRIVAL_FORMAT_SUBFIELD) return true; // This is FORMAT field which is empty or "." - all good

    char format[20];
    double f;

    // if we're validating a FORMAT field with --chain (in vcf_seg_validate_luft_trans_one_sample, if REF⇆ALT) - accept a valid scientific notation
    // as it will be converted to normal notation in vcf_seg_one_sample
    if (validate_only && chain_is_loaded && dict_id_is_vcf_format_sf (ctx->dict_id) &&
        str_scientific_to_decimal (STRa(recon), NULL, NULL, &f) && f >= 0.0 && f <= max_value) return true; // scientific notation in the valid range

    // if item format is inconsistent with AF being a probability value - we won't translate it
    if (!str_get_float (STRa(recon), &f, format, NULL) || f < 0.0 || f > max_value) 
        return false;
    
    if (validate_only) return true; 

    Ltxt -= recon_len;
    char f_str[50];
    sprintf (f_str, format, max_value - f);
    RECONSTRUCT (f_str, strlen (f_str)); // careful not to use bufprintf as it adds a \0 and we are required to translate in-place for all FORMAT fields
    
    return true;
}


// Lift-over translator for INFO/AF, FORMAT/AF and similar fields, IF it is bi-allelic and we have a ALT<>REF switch.
// We change the probability value to 1-AF
// returns true if successful (return value used only if validate_only)
TRANSLATOR_FUNC (vcf_piz_luft_A_1)
{
    return vcf_piz_luft_trans_complement_to_max_value (vb, ctx, STRa(recon), validate_only, 1);
}

//----------
// FORMAT/GL
// ---------

// convert an array of probabilities to an array of integer phred scores capped at 60
static void vcf_convert_prob_to_phred (VBlockVCFP vb, rom flag_name, STRp(snip), char *optimized_snip, unsigned *optimized_snip_len)
{
    str_split_floats (snip, snip_len, 0, ',', prob, false);
    ASSVCF (n_probs, "cannot to apply %s to value \"%.*s\"", flag_name, STRf(snip)); // not an array of floats - abort, because we already changed the FORMAT field

    unsigned phred_len = 0;
    for (unsigned i=0; i < n_probs; i++) {
        
        int64_t phred = MIN_(60, (int64_t)(((-probs[i]) * 10)+0.5)); // round to the nearest int, capped at 60

        phred_len += str_int (phred, &optimized_snip[phred_len]);
        if (i < n_probs - 1)
            optimized_snip[phred_len++] = ',';
    }

    *optimized_snip_len = phred_len;
}

// converts an array of phred scores (possibly floats) to integers capped at 60
static bool vcf_phred_optimize (rom snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len /* in / out */)
{
    str_split_floats (snip, len, 0, ',', item, false);
    if (!n_items) return false; // not an array of floats

    unsigned out_len = 0;

    for (unsigned i=0; i < n_items; i++) {
        int64_t new_phred = MIN_(60, (int64_t)(items[i] + 0.5));
        out_len += str_int (new_phred, &optimized_snip[out_len]);
        if (i < n_items-1) optimized_snip[out_len++] = ',';
    }

    *optimized_snip_len = out_len;
    return true;
}

//----------
// FORMAT/AB
// ---------

// <ID=AB,Number=1,Type=Float,Description="Allele balance for each het genotype",RendAlg="NONE">
// Expecting: '.' if channel is 0 or 2, AD0/(AD0+AD1) if 1, and any value if 3 which is placed in AB3
// If expectation is met, SPECIAL is segged in AB. 
static inline void vcf_seg_FORMAT_AB (VBlockVCFP vb, ContextP ctx, STRp(ab))
{
    int channel_i = vcf_seg_get_mux_channel_i (vb, false);
    bool is_0_or_2 = (channel_i==0 || channel_i==2); // note: dos02 doesn't change in case of DVCF REF⇆ALT switch
    bool ab_missing = ab_len==1 && *ab=='.';

    if (channel_i==-1               || // GT didn't produce a mux channel
        (channel_i==1 && z_is_dvcf) || // we can't handle channel 1 in dual coordinates (TODO: limit to REF⇆ALT switch)
        segconf.has[FORMAT_ADALL]   || // we can't handle AD0/AD1 peeking in AD is a delta vs ADALL
        (is_0_or_2 && !ab_missing)) {  // if channel is 0 or 2, we expected a '.' 
    
        seg_by_ctx (VB, STRa(ab), ctx, ab_len);
        return;
    }

    // prepare rollback data: we will verify channel 1 in vcf_seg_FORMAT_AB_verify_channel1 and rollback if necessary
    if (channel_i==1) {
        seg_set_last_txt (VB, ctx, STRa(ab));
        ctx_set_last_value (VB, ctx, (ValueType){.i = 1}); // need verification

        seg_create_rollback_point (VB, NULL, 1, FORMAT_AB); 
    }

    if (channel_i==3) { 
        seg_by_ctx (VB, STRa(ab_snip), ctx, 0);
        seg_by_did (VB, STRa(ab), FORMAT_AB3, ab_len);
    }

    else // channel 0,1,2
        seg_by_ctx (VB, STRa(ab_snip), ctx, ab_len);
}

static inline void calculate_AB (double ad0, double ad1, char *snip)
{
    double ab = 0.0000001 + ad0 / (ad0 + ad1); // +epsilon to bring number over the 0.01 mark if it is almost almost there 
    sprintf (snip, "%.*g", ab < 0.1 ? 1 : 2, ab); 
}

// this is called if we segged AB channel 1 as a SPECIAL
static inline void vcf_seg_FORMAT_AB_verify_channel1 (VBlockVCFP vb)
{
    ContextP ab_ctx  = CTX(FORMAT_AB);
    ContextP ad0_ctx = ECTX (con_AD.items[0].dict_id);
    ContextP ad1_ctx = ECTX (con_AD.items[1].dict_id);

    STRlast (ab_str, FORMAT_AB);

    // rollback if we don't have AD0, AD1 this line, or if their value is not as expected by the formula
    if (!ad0_ctx || !ad1_ctx) goto rollback;

    if (!ctx_has_value (VB, ad0_ctx->did_i) || !ctx_has_value (VB, ad1_ctx->did_i)) goto rollback;

    double ad0 = ad0_ctx->last_value.i;
    double ad1 = ad1_ctx->last_value.i;
    if (ad0==0 && ad1==0) goto rollback; // formula would be division by zero 

    char recon_ab_str[32];
    calculate_AB (ad0, ad1, recon_ab_str);

    if (!str_issame_(STRa(ab_str), recon_ab_str, strlen (recon_ab_str))) goto rollback;
    
    ctx_unset_rollback (ab_ctx); // so we can ctx_set_rollback again in the next sample (note: we can't advance rback_id because it is set in vcf_seg_txt_line)
    return; // verified

rollback:
    seg_rollback (VB);
    seg_by_ctx (VB, STRa(ab_str), ab_ctx, ab_str_len); 
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_FORMAT_AB)
{
    int channel_i = vcf_piz_get_mux_channel_i (vb);

    if ((channel_i == 0 || channel_i == 2) && reconstruct) 
        RECONSTRUCT1 ('.');

    else if (channel_i == 1 && reconstruct) {
        ContextP ad0_ctx = MCTX (0, snip, snip_len);
        ContextP ad1_ctx = MCTX (1, snip, snip_len);

        double ad0 = reconstruct_peek (vb, ad0_ctx, 0, 0).i;
        double ad1 = reconstruct_peek (vb, ad1_ctx, 0, 0).i;
        ASSPIZ0 (ad0 || ad1, "Unexpectedly, ad0=ad1=0");

        char recon_ab_str[32];
        calculate_AB (ad0, ad1, recon_ab_str);
        RECONSTRUCT (recon_ab_str, strlen(recon_ab_str));
    }

    else if (channel_i == 3) {
        ContextP ab3_ctx = MCTX (2, snip, snip_len);
        reconstruct_from_ctx (vb, ab3_ctx->did_i, 0, reconstruct);
    }

    return NO_NEW_VALUE;
}

//----------
// FORMAT/PL
// ---------

// <ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">       
static inline void vcf_seg_FORMAT_PL (VBlockVCFP vb, ContextP ctx, STRp(PL))
{
    if (segconf.running && !segconf.has_DP_before_PL) 
        segconf.has_DP_before_PL = ctx_encountered (VB, FORMAT_DP);

    seg_set_last_txt (VB, ctx, STRa(PL)); // used by GQ and SPL (points into txt_data, before optimization)

    // attempt to optimize PL string, if requested
    unsigned optimized_len = PL_len*2 + 10;                 
    char optimized[PL_len]; // note: modifying functions need to make sure not to overflow this space

    if (flag.optimize_phred && vcf_phred_optimize (STRa(PL), qSTRa(optimized))) {
        int growth = (int)optimized_len - (int)PL_len;
        vb->recon_size      += growth;               
        vb->recon_size_luft += growth;
        STRset(PL, optimized);
    }
        
    // seg into either PLy or PLn, or into both if we're testing (vcf_FORMAT_PL_decide will drop one of them) 
    if (vb->PL_mux_by_DP == PL_mux_by_DP_YES || vb->PL_mux_by_DP == PL_mux_by_DP_TEST) 
        vcf_seg_FORMAT_mux_by_dosagexDP (vb, CTX(FORMAT_PLy), STRa(PL), &vb->mux_PLy);
    
    if (vb->PL_mux_by_DP == PL_mux_by_DP_NO  || vb->PL_mux_by_DP == PL_mux_by_DP_TEST) 
        vcf_seg_FORMAT_mux_by_dosage (vb, CTX(FORMAT_PLn), STRa(PL), &vb->mux_PLn);

    if (vb->PL_mux_by_DP == PL_mux_by_DP_TEST) 
        vb->recon_size += PL_len; // since we're segging twice, we need to pretent the recon_size is growing even thow it isn't (reversed in zfile_remove_ctx_group_from_z_data)

    // we seg "redirection to PLn" in all VBs regardless of PL_mux_by_DP, so that ZCTX(FORMAT_PL).dict has
    // only a single word. This word might be updated to PLy in vcf_FORMAT_PL_after_vbs.
    // note: this isn't always "all-the-same" - we can have WORD_INDEX_MISSING in b250 in case of missing PLs in samples
    seg_by_did (VB, STRa(PL_to_PLn_redirect_snip), FORMAT_PL, 0);
}

// in PL_mux_by_DP_TEST mode, this function is called to pick one of the two compressed options, and drop the others
void vcf_FORMAT_PL_decide (VBlockVCFP vb)
{
    Did keep_did_i, remove_did_i;

    mutex_lock (segconf.PL_mux_by_DP_mutex);

    // only first testing VB to lock this mutex gets here, to make the decision
    switch (segconf.PL_mux_by_DP) {
        case PL_mux_by_DP_TEST : { 
            // get combined compress size of all contexts associated with GL (mux by dosage)
            uint64_t PLy_size = ctx_get_ctx_group_z_len (VB, FORMAT_PLy);
            uint64_t PLn_size = ctx_get_ctx_group_z_len (VB, FORMAT_PLn);

            // note: this is not an accurate comparison bc it doesn't include dictionaries
            keep_did_i   = (PLy_size >= PLn_size) ? FORMAT_PLn : FORMAT_PLy;
            remove_did_i = (PLy_size >= PLn_size) ? FORMAT_PLy : FORMAT_PLn;

            // weigh in this VBs experience. New VBs will only test if the votes are still the same
            segconf.PL_mux_by_DP = (keep_did_i == FORMAT_PLy) ? PL_mux_by_DP_YES : PL_mux_by_DP_NO;
            
            if (flag.debug_generate) iprintf ("PL_mux_by_DP decision: %s\n", keep_did_i == FORMAT_PLy ? "YES" : "NO");
            
            break;
        }
        case PL_mux_by_DP_YES : keep_did_i = FORMAT_PLy; remove_did_i = FORMAT_PLn; break;

        default:
        case PL_mux_by_DP_NO  : keep_did_i = FORMAT_PLn; remove_did_i = FORMAT_PLy;
    }

    mutex_unlock (segconf.PL_mux_by_DP_mutex);

    zfile_remove_ctx_group_from_z_data (VB, remove_did_i); // removes data from VB and update z_file counters
}

// called after all VBs are compressed - before Global sections are compressed
void vcf_FORMAT_PL_after_vbs (void)
{
    if (!ZCTX(FORMAT_PL)->nodes.len) return; // no FORMAT/PL in this file

    if (ZCTX(FORMAT_PL)->nodes.len > 1) {
        str_print_dict (info_stream, STRb(ZCTX(FORMAT_PL)->dict), true, false);
        ABORT ("Expecting FORMAT_PL to have exactly one word in its dict, but it has %"PRIu64, ZCTX(FORMAT_PL)->nodes.len);
    }

    if (segconf.PL_mux_by_DP == PL_mux_by_DP_YES) { // Tested, and selected YES
        ctx_declare_winning_group (FORMAT_PLy, FORMAT_PLn, FORMAT_PL); // update groups for stats

        // if we select "Yes", we change the (only) snip in the PL dictionary in place
        memcpy (ZCTX(FORMAT_PL)->dict.data, PL_to_PLy_redirect_snip, PL_to_PLy_redirect_snip_len);
    }
    else // Tested, and selected NO, or not tested
        ctx_declare_winning_group (FORMAT_PLn, FORMAT_PLy, FORMAT_PL); 
}

//----------
// FORMAT/DS
// ---------

/* used by files compressed by Genozip up to 12.0.42. Kept here, commented out, for documentation 
// DS is a value [0,ploidy] which is a the sum of GP values that refer to allele!=0. Eg for bi-allelic diploid it is: GP[1] + 2*GP[2] (see: https://www.biostars.org/p/227346/)
// the DS (allele DoSage) value is usually close to or exactly the sum of '1' alleles in GT. we store it as a delta from that,
// along with the floating point format to allow exact reconstruction
static inline WordIndex vcf_seg_FORMAT_DS (VBlockVCFP vb, ContextP ctx, rom cell, unsigned cell_len)
{
    int64_t dosage = ctx_has_value (VB, FORMAT_GT) ? CTX(FORMAT_GT)->last_value.i : -1; // dosage stored here by vcf_seg_FORMAT_GT
    double ds_val;
    unsigned format_len;
    char snip[FLOAT_FORMAT_LEN + 20] = { SNIP_SPECIAL, VCF_SPECIAL_DS }; 

    if (dosage < 0 || !str_get_float (cell, cell_len, &ds_val, &snip[2], &format_len)) 
        return seg_by_ctx (VB, cell, cell_len, ctx, cell_len);

    unsigned snip_len = 2 + format_len;
    snip[snip_len++] = ' ';
    snip_len += str_int ((int64_t)((ds_val - dosage) * 1000000), &snip[snip_len]);

    return seg_by_ctx (VB, snip, snip_len, ctx, cell_len);
}
*/

// used for decompressing files compressed with version up to 12.0.42
SPECIAL_RECONSTRUCTOR (vcf_piz_special_FORMAT_DS_old)
{
    if (!reconstruct) goto done;

    char float_format[10];
    int32_t val;
    sscanf (snip, "%s %d", float_format, &val); // snip looks like eg: "%5.3f 50000"

    unsigned dosage = vcf_piz_get_mux_channel_i (vb); // seg guaranteed 0, 1 or 2 (or else this SPECIAL is not used)
    bufprintf (vb, &vb->txt_data, float_format, (double)val / 1000000 + dosage);

done:
    return NO_NEW_VALUE;
}

// Lift-over translator for FORMAT/DS, IF it is bi-allelic and we have a ALT<>REF switch.
// We change the value to (ploidy-value)
// returns true if successful (return value used only if validate_only)
TRANSLATOR_FUNC (vcf_piz_luft_PLOIDY)
{
    if (!ctx_encountered (VB, FORMAT_GT)) return false; // we can't translate unless this variant as GT

    // use gt_prev_ploidy: in Seg, set by vcf_seg_FORMAT_GT, in validate and piz set by vcf_piz_luft_GT 
    return vcf_piz_luft_trans_complement_to_max_value (vb, ctx, STRa(recon), validate_only, CTX(FORMAT_GT)->gt_actual_last_ploidy);
}

//------------------------------------------------
// FORMAT and INFO - subfields with G, R, R2 types
//------------------------------------------------

// Lift-over ALT<>REF switch translator for bi-allelic multi-value fields: 
// three cases: (1) R1,R2->R2,R1 (2) Ra1,Rb1,Ra2,Rb2->Ra2,Rb2,Ra1,Rb1 (3) G11,G12,G22->G22,G12,G11
// returns true if successful (return value used only if validate_only)
static bool vcf_piz_luft_switch_first_last (VBlockP vb, ContextP ctx, char *recon, int32_t recon_len, 
                                            unsigned num_items, char field_type, bool validate_only)
{
    if (IS_TRIVAL_FORMAT_SUBFIELD) return true; // This is FORMAT field which is empty or "." - all good

    char copy[recon_len];
    memcpy (copy, recon, recon_len);
    
    str_split (copy, recon_len, num_items, ',', item, true);
    if (!n_items) return false; // if item format is inconsistent with VCF header - we won't translate it

    if (validate_only) return true;
    
    Ltxt -= recon_len;
    
    if (num_items==2 || num_items == 3) {
        RECONSTRUCT_SEP (items[n_items-1], item_lens[n_items-1], ',');

        if (num_items==3)
            RECONSTRUCT_SEP (items[1], item_lens[1], ',');

        RECONSTRUCT (items[0], item_lens[0]);
    }
    else if (num_items == 4) { // Ra1,Rb1,Ra2,Rb2 -> Ra2,Rb2,Ra1,Rb1
        RECONSTRUCT_SEP (items[2], item_lens[2], ',');
        RECONSTRUCT_SEP (items[3], item_lens[3], ',');
        RECONSTRUCT_SEP (items[0], item_lens[0], ',');
        RECONSTRUCT     (items[1], item_lens[1]);
    }
    
    return true;
}

// Lift-over translator assigned to a Number=R item, IF it is bi-allelic and we have a ALT<>REF switch.
// 'R'   : We switch between the two comma-separated values.
// 'R2'  : We switch between the two PAIRS of comma-separated values.
// 'G'   : We have 3 values which represent the genotypes REF/REF,REF/ALT,ALT/ALT We switch between the 1st and 3rd value.
// 'NEG' : We negate a numeric value
// returns true if successful 
TRANSLATOR_FUNC (vcf_piz_luft_R)  { return vcf_piz_luft_switch_first_last (vb, ctx, STRa(recon), 2, 'R', validate_only); } // 2 bc we only handle bi-allelic
TRANSLATOR_FUNC (vcf_piz_luft_R2) { return vcf_piz_luft_switch_first_last (vb, ctx, STRa(recon), 4, '.', validate_only); } // 4 bc we only handle bi-allelic

TRANSLATOR_FUNC (vcf_piz_luft_G)  
{ 
    if (IS_TRIVAL_FORMAT_SUBFIELD) return true; // This is FORMAT field which is empty or "." - all good

    unsigned num_values = str_count_char (STRa(recon), ',')+1;
    if (num_values != 3 && num_values != 2) return false; // Genozip currently only support haploid (2 bi-allelic genotypes) and diploid (3 bi-allelic genotypes) 

    return vcf_piz_luft_switch_first_last (vb, ctx, STRa(recon), num_values, 'G', validate_only); 
}

//--------------------------------------------------------------------------------------------------------------
// LongRanger: <ID=BX,Number=.,Type=String,Description="Barcodes and Associated Qual-Scores Supporting Alleles">
// example: CCTAAAGGTATCGCCG-1_41;TTGTCCGTCGCTAGCG-1_55;TATCATCGTTGGAGGT-1_74
// See: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf
//--------------------------------------------------------------------------------------------------------------
static inline void vcf_seg_FORMAT_BX (VBlockVCFP vb, ContextP ctx, STRp(BX))
{
    static const MediumContainer con = {
        .nitems_lo   = 2, 
        .drop_final_repsep = true,
        .repsep      = {';'},
        .items       = { { .dict_id={ .id="BXbarcod" }, .separator = {'-'} },
                         { .dict_id={ .id="BXqual"   },                    } }
    };
    
    ctx_get_ctx (VB, con.items[0].dict_id)->no_stons = true;
    
    seg_array_of_array_of_struct (VB, CTX(FORMAT_BX), ',', con, STRa(BX), NULL);
}

//------------------------------------------------------------------------
// Validate that ALL subfields in ALL samples can luft-translate as needed
//------------------------------------------------------------------------

static rom error_format_field (unsigned n_items, ContextP *ctxs)
{
    static char format[256];
    unsigned len=0;
    for (unsigned i=0; i < n_items; i++) 
        len += strlen (ctxs[i]->tag_name) + 1;

    if (len > sizeof format-1) return "<FORMAT too long to display>";

    len=0;
    for (unsigned i=0; i < n_items; i++) {
        unsigned one_len = strlen (ctxs[i]->tag_name);
        memcpy (&format[len], ctxs[i]->tag_name, one_len); 
        len += one_len;
        format[len++] = ':';
    } 

    format[len-1] = 0;
    return format;
}

// if any context fails luft-translation, returns that context, or if all is good, returns NULL
static inline ContextP vcf_seg_validate_luft_trans_one_sample (VBlockVCFP vb, ContextP *ctxs, uint32_t num_items, char *sample, unsigned sample_len)
{
    str_split (sample, sample_len, num_items, ':', item, false);
    ASSVCF (n_items, "Sample %u has too many subfields - FORMAT field \"%s\" specifies only %u: \"%.*s\"", 
            vb->sample_i+1, error_format_field (num_items, ctxs), num_items, sample_len, sample);

    ContextP failed_ctx = NULL; // optimistic initialization - nothing failed

    uint32_t save_ploidy = CTX(FORMAT_GT)->gt_prev_ploidy; // ruined by vcf_piz_luft_GT 

    for (unsigned i=0; i < n_items; i++) {
        if (needs_translation (ctxs[i]) && item_lens[i]) {

            if ((vb->line_coords == DC_LUFT && !vcf_lo_seg_cross_render_to_primary (vb, ctxs[i], STRi(item,i), NULL, 0, true/*validate only*/)) ||
                (vb->line_coords == DC_PRIMARY && !(DT_FUNC(vb, translator)[ctxs[i]->luft_trans](VB, ctxs[i], (char *)STRi(item,i), 0, true)))) {
                failed_ctx = ctxs[i];  // failed translation
                break;
            }
        }
        ctx_set_encountered (VB, ctxs[i]); // might be needed for validation 
    }

    // restore values, in preparation for real Seg
    CTX(FORMAT_GT)->gt_prev_ploidy = save_ploidy;
    for (unsigned i=0; i < n_items; i++) 
        ctx_unset_encountered (VB, ctxs[i]); 

    return failed_ctx; 
}

// If ALL subfields in ALL samples can luft-translate as required and sets ctx->line_is_luft_trans for all contexts
// if NOT: ctx->line_is_luft_trans=false for all contexts, line is rejects (LO_FORMAT), and keeps samples in their original LUFT or PRIMARY coordinates.
static inline void vcf_seg_validate_luft_trans_all_samples (VBlockVCFP vb, uint32_t n_items, ContextP *ctxs, 
                                                            int32_t len, char *samples_start)
{
    rom field_start, next_field = samples_start;
    unsigned field_len=0;
    bool has_13;

    // initialize optimistically. we will roll back and set to false if ANY subfield in ANY sample fails to translate, and re-seg all samples
    for (unsigned sf_i=0; sf_i < n_items; sf_i++)
        ctxs[sf_i]->line_is_luft_trans = needs_translation (ctxs[sf_i]); 

    // 0 or more samples
    vb->sample_i=0;
    for (char separator=0 ; separator != '\n'; vb->sample_i++) {

        field_start = next_field;
        next_field = seg_get_next_item (VB, field_start, &len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, &has_13, "sample-subfield");
        ASSVCF (field_len, "unexpected tab character after sample # %u", vb->sample_i);

        ContextP failed_ctx = vcf_seg_validate_luft_trans_one_sample (vb, ctxs, n_items, (char *)field_start, field_len);
        if (failed_ctx) { // some context doesn't luft-translate as required
            REJECT_SUBFIELD (LO_FORMAT, failed_ctx, ".\tCannot cross-render sample due to field %s: \"%.*s\"", failed_ctx->tag_name, field_len, field_start);

            // make all contexts untranslateable in this line
            for (unsigned i=0; i < n_items; i++)  // iterate on the order as in the line
                ctxs[i]->line_is_luft_trans = false;
        }
    }
}

// ----------
// One sample
// ----------

// returns true is the sample has a non-'.' PS value
static inline bool vcf_seg_sample_has_PS (VBlockVCFP vb, ContextP *ctxs, STRps(sf))
{
    ContextP ps_ctx = CTX(FORMAT_PS);
    int16_t sf_i = ps_ctx->sf_i;

    if (sf_i == -1     ||                        // no SF in this line's FORMAT
        sf_i >= n_sfs  ||                        // SF field is missing in this sample despite being in FORMAT
        !sf_lens[sf_i] ||                        // SF is "" (i.e. empty)
        (sf_lens[sf_i]==1 && sfs[sf_i][0]=='.')) // SF is "."
        return false;

    // first sample in the VB is a proper PS - get its type 
    if (!ps_ctx->ps_type)
        vcf_samples_seg_initialize_PS_PID (vb, ps_ctx, STRi(sf, sf_i));

    return true;
}

// returns the number of colons in the sample
static inline unsigned vcf_seg_one_sample (VBlockVCFP vb, ZipDataLineVCF *dl, ContextP *ctxs, ContainerP samples, STRp(sample))
{
    #define COND0(condition, seg) if (condition) { seg; break; } else  
    #define COND(condition,  seg) if (condition) { seg; break; } else goto fallback; 

    str_split (sample, sample_len, con_nitems (*samples), ':', sf, false);

    ASSVCF (n_sfs, "Sample %u has too many subfields - FORMAT field \"%s\" specifies only %u: \"%.*s\"", 
            vb->sample_i+1, error_format_field (con_nitems (*samples), ctxs), con_nitems (*samples), STRf(sample));

    for (unsigned i=0; i < n_sfs; i++) { 

        DictId dict_id = samples->items[i].dict_id;
        ContextP ctx = ctxs[i];

        STRli(translated, sf_lens[i]*2 + 10); // for translations 
        STRli(optimized,  sf_lens[i]*2 + 10); // for optimizations - separate space - we can optimize a "translated" string

        #define SEG_OPTIMIZED_MUX_BY_DOSAGE(tag) ({                     \
            vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRa(optimized), &vb->mux_##tag);  \
            int32_t growth = (int)optimized_len - (int)sf_lens[i];      \
            vb->recon_size      += growth;                              \
            vb->recon_size_luft += growth; })

        // --chain: if this is RendAlg=A_1 and RendAlg=PLOIDY subfield, convert a eg 4.31e-03 to e.g. 0.00431. This is to
        // ensure primary->luft->primary is lossless (4.31e-03 cannot be converted losslessly as we can't preserve format info)
        // Warning: last_txt will be invalid for this context
        if (chain_is_loaded && 
            (ctx->luft_trans == VCF2VCF_A_1 || ctx->luft_trans == VCF2VCF_PLOIDY) && 
            str_scientific_to_decimal (STRi(sf, i), qSTRa(translated), NULL)) {          

            int32_t growth = (int32_t)translated_len - (int32_t)sf_lens[i];    
            if (growth) {                                               
                vb->recon_size      += growth;                          
                vb->recon_size_luft += growth;                          
            }                                                           
            sfs[i] = translated;                                                 
            sf_lens[i] = translated_len; 
        }

        if (vb->line_coords == DC_LUFT && needs_translation (ctx)) {

            // case: translate into "translated" if the translation might modify the size
            // Warning: last_txt will be invalid for this context
            if (ctx->luft_trans == VCF2VCF_NEG) {
                ASSERT (vcf_lo_seg_cross_render_to_primary (vb, ctx, STRi(sf,i), qSTRa(translated), false),
                        "Expecting translation of %s to succeed, because it was already validated", ctx->tag_name);
                
                sfs[i] = translated;   // note: recon_size was adjusted in vcf_lo_seg_cross_render_to_primary                                   
                sf_lens[i] = translated_len; 
            }

            else
                ASSERT (vcf_lo_seg_cross_render_to_primary (vb, ctx, STRi(sf,i), NULL, 0, false),
                        "Expecting translation of %s to succeed, because it was already validated", ctx->tag_name);
        }

        if (!sf_lens[i])
            seg_by_ctx (VB, "", 0, ctx, 0); // generates WORD_INDEX_EMPTY

        else switch (dict_id.num) {

        // <ID=GT,Number=1,Type=String,Description="Genotype">
        case _FORMAT_GT: {
            bool has_ps = vcf_seg_sample_has_PS (vb, ctxs, STRas(sf));
            
            bool has_null_dp = segconf.use_null_DP_method ? vcf_seg_sample_has_null_value (_FORMAT_DP, ctxs, STRas(sf)) 
                                                          : false;

            vcf_seg_FORMAT_GT (vb, ctx, dl, STRi(sf, i), has_ps, has_null_dp); 
            break;
        }

        // <ID=GL,Number=.,Type=Float,Description="Genotype Likelihoods">
        case _FORMAT_GL:
            // --GL-to-PL:  GL: 0.00,-0.60,-8.40 -> PL: 0,6,60
            // note: we changed the FORMAT field GL->PL in vcf_seg_format_field. data is still stored in the GL context.
            // <ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
            if (flag.GL_to_PL) {
                vcf_convert_prob_to_phred (vb, "--GL-to-PL", STRi(sf, i), qSTRa(optimized));
                SEG_OPTIMIZED_MUX_BY_DOSAGE(GL);
            }
            else // I tried muxing against DS instead of dosage (41 channels) - worse results than dosage in 1KGP-37 even with --best 
                vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_GL);
            break;

        // note: GP and PL - for non-optimized, I tested segging as A_R and seg_array - they are worse or not better than the default. likely because the values are correlated.
        case _FORMAT_GP:
            // convert GP (probabilities) to PP (phred values). PP was introduced in VCF v4.3.
            if (flag.GP_to_PP && vb->vcf_version >= VCF_v4_3) {
                vcf_convert_prob_to_phred (vb, "--GP-to-PP", STRi(sf, i), qSTRa(optimized));
                SEG_OPTIMIZED_MUX_BY_DOSAGE(GP);
            }
            else if (flag.optimize_phred && vb->vcf_version <= VCF_v4_2 &&
                     vcf_phred_optimize (STRi(sf, i), qSTRa(optimized))) 
                SEG_OPTIMIZED_MUX_BY_DOSAGE(GP);
            else
                vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_GP);

            seg_set_last_txt (VB, ctx, STRi(sf, i)); // used by GQ
            break;

        // <ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">       
        case _FORMAT_PL: vcf_seg_FORMAT_PL (vb, ctx, STRi (sf, i)); break;

        // <ID=PP,Number=G,Type=Integer,Description="Phred-scaled genotype posterior probabilities rounded to the closest integer">
        case _FORMAT_PP:
            if (segconf.vcf_is_pindel)
                goto fallback; // PP means something entirely different in Pindel
            else if (flag.optimize_phred && vcf_phred_optimize (STRi(sf, i), qSTRa(optimized))) 
                SEG_OPTIMIZED_MUX_BY_DOSAGE(PP);
            else
                vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_PP);
            break;
        
        // <ID=PRI,Number=G,Type=Float,Description="Phred-scaled prior probabilities for genotypes">
        case _FORMAT_PRI:
            if (flag.optimize_phred && vcf_phred_optimize (STRi(sf, i), qSTRa(optimized))) 
                SEG_OPTIMIZED_MUX_BY_DOSAGE(PRI);
            else
                vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_PRI);
            break;
        
        // <ID=DS,Number=1,Type=Float,Description="Genotype dosage from MaCH/Thunder"> (1000 Genome Project phase1 data)
        // See: https://genome.sph.umich.edu/wiki/Thunder
        case _FORMAT_DS   : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_DS)   ; break;
        
        // VarScan: <ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
        case _FORMAT_RD   : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_RD)   ; break;
        
        // VarScan: <ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
        case _FORMAT_PVAL : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_PVAL) ; break;   
        
        // VarScan: <ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
        case _FORMAT_FREQ : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_FREQ) ; break;
        
        #define ILLUM_GTYPING_MUX_BY_DOSAGE(tag) \
            if (segconf.vcf_illum_gtyping) { vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_##tag); break; }\
            else goto fallback;

        // Illumina Genotyping: <ID=BAF,Number=1,Type=Float,Description="B Allele Frequency">
        case _FORMAT_BAF  : vcf_seg_mux_by_adjusted_dosage (vb, ctx, STRi (sf, i), &vb->mux_BAF); break;

        // Illumina Genotyping: <ID=X,Number=1,Type=Integer,Description="Raw X intensity"> (same for Y) 
        case _FORMAT_X    : vcf_seg_mux_by_adjusted_dosage (vb, ctx, STRi (sf, i), &vb->mux_X); break;
        case _FORMAT_Y    : vcf_seg_mux_by_adjusted_dosage (vb, ctx, STRi (sf, i), &vb->mux_Y); break;
        
        case _FORMAT_PS   : 
        case _FORMAT_PID  : vcf_seg_FORMAT_PS_PID (vb, dl, ctx, STRi(sf, i)); break;

        // standard: <ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        // GIAB: <ID=GQ,Number=1,Type=Integer,Description="Net Genotype quality across all datasets, calculated from GQ scores of callsets supporting the consensus GT, using only one callset from each dataset">   
        case _FORMAT_GQ   : if (segconf.vcf_is_isaac) seg_set_last_txt_store_value (VB, ctx, STRi(sf, i), STORE_INT);
                            else                      seg_set_last_txt (VB, ctx, STRi(sf, i)); 
                            break; // postpone to later
            
        case _FORMAT_RGQ  : vcf_seg_FORMAT_RGQ (vb, ctx, STRi(sf, i), ctxs[0], STRi(sf,0)); break;

        // <ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        case _FORMAT_DP   : vcf_seg_FORMAT_DP (vb, ctx, STRi(sf, i)); break;
            
        // <ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
        // case: MIN_DP - it is slightly smaller and usually equal to DP - we store MIN_DP as the delta DP-MIN_DP
        // note: the delta is vs. the DP field that preceeds MIN_DP - we take the DP as 0 there is no DP that preceeds
        case _FORMAT_MIN_DP :
            COND (ctx_has_value (VB, FORMAT_DP), seg_delta_vs_other (VB, ctx, CTX(FORMAT_DP), STRi(sf, i)));

        case _FORMAT_SDP   :
            if (ctx_has_value_in_line_(VB, CTX(INFO_ADP)))
                seg_delta_vs_other (VB, ctx, CTX(INFO_ADP), STRi(sf, i));
            else goto fallback;
            break;
            
        case _FORMAT_AF    : vcf_seg_FORMAT_AF (vb, ctx, STRi(sf, i)); break;

        // standard: <ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">  
        // GIAB: <ID=AD,Number=R,Type=Integer,Description="Net allele depths across all unfiltered datasets with called genotype">
        case _FORMAT_AD    : if (segconf.vcf_is_varscan)
                                vcf_seg_FORMAT_AD_varscan (vb, ctx, STRi(sf, i));
                             else
                                vcf_seg_FORMAT_A_R (vb, ctx, con_AD, STRi(sf, i), STORE_INT, vcf_seg_AD_items); 
                             break;

        // <ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
        case _FORMAT_PGT   :  vcf_seg_FORMAT_PGT (vb, ctx, STRi(sf, i), ctxs, STRas(sf)); break;

        // VarScan: <ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">
        case _FORMAT_ADF   : 
            if (segconf.vcf_is_varscan)
                vcf_seg_FORMAT_minus (vb, ctx, STRi(sf, i), 0, CTX(FORMAT_AD), CTX(FORMAT_ADR), STRa(adf_snip)); 
            else
                vcf_seg_FORMAT_A_R (vb, ctx, con_ADF, STRi(sf, i), STORE_INT, vcf_seg_ADF_items); 
            break;

        // VarScan: <ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">                             
        case _FORMAT_ADR   : 
            if (segconf.vcf_is_varscan)
                vcf_seg_FORMAT_minus (vb, ctx, STRi(sf, i), 0, CTX(FORMAT_AD), CTX(FORMAT_ADF), STRa(adr_snip)); 
            else
                vcf_seg_FORMAT_A_R (vb, ctx, con_ADR, STRi(sf, i), STORE_INT, vcf_seg_ADR_items); 
            break;

        // <ID=F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting each allele">
        case _FORMAT_F1R2  : vcf_seg_FORMAT_A_R (vb, ctx, con_F1R2, STRi(sf, i), STORE_INT, NULL); break;

        // <ID=F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting each allele">
        case _FORMAT_F2R1  : vcf_seg_FORMAT_A_R (vb, ctx, con_F2R1, STRi(sf, i), STORE_INT, vcf_seg_F2R1_items); break;

        // <ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias">
        case _FORMAT_SB    : vcf_seg_FORMAT_A_R (vb, ctx, con_SB, STRi(sf, i), STORE_NONE, vcf_seg_SB_items); break;

        // <ID=MB,Number=4,Type=Integer,Description="Per-sample component statistics to detect mate bias">
        case _FORMAT_MB    : vcf_seg_FORMAT_A_R (vb, ctx, con_MB, STRi(sf, i), STORE_INT, vcf_seg_MB_items); break;

        // <ID=SAC,Number=.,Type=Integer,Description="Number of reads on the forward and reverse strand supporting each allele (including reference)">
        case _FORMAT_SAC   : vcf_seg_FORMAT_A_R (vb, ctx, con_SAC, STRi(sf, i), STORE_NONE, vcf_seg_SAC_items); break;

        // <ID=ICNT,Number=2,Type=Integer,Description="Counts of INDEL informative reads based on the reference confidence model">
        case _FORMAT_ICNT  : COND (segconf.vcf_is_gvcf, vcf_seg_ICNT (vb, ctx, STRi(sf, i))); 

        // <ID=SPL,Number=.,Type=Integer,Description="Normalized, Phred-scaled likelihoods for SNPs based on the reference confidence model">
        case _FORMAT_SPL  : COND (segconf.vcf_is_gvcf, vcf_seg_SPL (vb, ctx, STRi(sf, i))); 

        // VarScan: <ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
        case _FORMAT_RDF   : vcf_seg_FORMAT_minus (vb, ctx, STRi(sf, i), 0, CTX(FORMAT_RD), CTX(FORMAT_RDR), STRa(rdf_snip)); break;

        // VarScan: <ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
        case _FORMAT_RDR   : vcf_seg_FORMAT_minus (vb, ctx, STRi(sf, i), 0, CTX(FORMAT_RD), CTX(FORMAT_RDF), STRa(rdr_snip)); break;

        // <ID=AB,Number=1,Type=Float,Description="Allele balance for each het genotype">
        case _FORMAT_AB    : vcf_seg_FORMAT_AB (vb, ctx, STRi(sf, i)); break;

        // <ID=BX,Number=.,Type=String,Description="Barcodes and Associated Qual-Scores Supporting Alleles">
        // example: CCTAAAGGTATCGCCG-1_41;TTGTCCGTCGCTAGCG-1_55;TATCATCGTTGGAGGT-1_74
        case _FORMAT_BX    : vcf_seg_FORMAT_BX (vb, ctx, STRi(sf, i)); break;
        
        // GWAS-VCF fields
        #define IF_GWAS(f) ({ if (segconf.vcf_is_gwas) {f; break;} else goto fallback; })
        
        // VCF-GWAS: <ID=ID,Number=1,Type=String,Description="Study variant identifier">
        case _FORMAT_ID    : IF_GWAS(vcf_gwas_seg_FORMAT_ID (vb, ctx, STRi(sf, i)));

        // GIAB fields
        case _FORMAT_ADALL : vcf_seg_FORMAT_A_R (vb, ctx, con_ADALL, STRi(sf, i), STORE_INT, vcf_seg_ADALL_items); break;
        case _FORMAT_IGT   : COND (segconf.vcf_is_giab_trio, vcf_seg_FORMAT_IGT (vb, ctx, STRi(sf, i))); 
        case _FORMAT_IPS   : COND (segconf.vcf_is_giab_trio, vcf_seg_FORMAT_IPS (vb, dl, ctx, STRi(sf, i))); 

        // Illumina ISAAC fields
        case _FORMAT_GQX   : COND (segconf.vcf_is_isaac, vcf_seg_FORMAT_GQX (vb, ctx, STRi(sf, i)));

        // DRAGEN fields
        case _FORMAT_PE    : COND (segconf.vcf_is_dragen, seg_array (VB, ctx, ctx->did_i, STRi(sf, i), ',', 0, false, STORE_INT, DICT_ID_NONE, sf_lens[i]));
        case _FORMAT_BC    : 
        case _FORMAT_CN    : COND (segconf.vcf_is_dragen, seg_integer_or_not (VB, ctx, STRi(sf, i), sf_lens[i]));

        // manta fields
        case _FORMAT_SR    :
        case _FORMAT_PR    : COND (segconf.vcf_is_manta, seg_array (VB, ctx, ctx->did_i, STRi(sf, i), ',', 0, false, STORE_INT, DICT_ID_NONE, sf_lens[i]));

        default            :
        fallback           : seg_by_ctx (VB, STRi(sf, i), ctx, sf_lens[i]);
        }

        int64_t value;
        if (ctx->flags.store == STORE_INT && !ctx_has_value(VB, ctx->did_i) &&  // not already set
            str_get_int (STRi(sf, i), &value))
            ctx_set_last_value (VB, ctx, value);
        else        
            ctx_set_encountered (VB, ctx);
    }

    // missing subfields - defined in FORMAT but missing (not merely empty) in sample
    for (unsigned i=n_sfs; i < con_nitems (*samples); i++) {
        uint64_t dnum = samples->items[i].dict_id.num;

        // special handling for PS and PID
        if (dnum == _FORMAT_PS || dnum == _FORMAT_PID) 
            vcf_seg_FORMAT_PS_PID_missing_value (vb, ctxs[i], &sample[sample_len]);
        
        else
            seg_by_ctx (VB, NULL, 0, ctxs[i], 0); // generates WORD_INDEX_MISSING
    }
    
    // verify AB if its channel 1
    if (ctx_has_value (VB, FORMAT_AB))
        vcf_seg_FORMAT_AB_verify_channel1 (vb);

    // finally seg GQ if we have it (after GP, PL, DP)
    if (ctx_encountered (VB, FORMAT_GQ))
        vcf_seg_FORMAT_GQ (vb);

    return n_sfs - 1; // number of colons
}

//------------
// All samples
//------------

rom vcf_seg_samples (VBlockVCFP vb, ZipDataLineVCF *dl, int32_t *len, char *next_field, bool *has_13)
{
    // Container for samples - we have:
    // - repeats as the number of samples in the line (<= vcf_num_samples)
    // - n_items as the number of FORMAT subfields (inc. GT)

    Container samples = *B(Container, vb->format_mapper_buf, dl->format_node_i); // make a copy of the template
    ContextP *ctxs = B(ContextP, vb->format_contexts, dl->format_node_i * MAX_FIELDS);
    uint32_t n_items = con_nitems (samples);

    // check that all subfields in all samples can be luft-translated as required, or make this a LUFT-only / PRIMARY-only line.
    // Also, if the data is in LUFT coordinates and is indeed translatable, then this lifts-back the samples to PRIMARY coordinates
    if (z_is_dvcf && LO_IS_OK (last_ostatus))
        vcf_seg_validate_luft_trans_all_samples (vb, n_items, ctxs, *len, next_field);
        
    // set ctx->sf_i - for the line's FORMAT fields
    for (int sf_i=0; sf_i < con_nitems(samples); sf_i++) {
        ctxs[sf_i]->sf_i = sf_i;

        if (segconf.running) segconf.has[ctxs[sf_i]->did_i]++;
    }

    // initialize LOOKBACK if we have PS or PID
    if (!CTX(VCF_LOOKBACK)->is_initialized && (CTX(FORMAT_PID)->sf_i >= 0 || CTX(FORMAT_PS)->sf_i >= 0))
        vcf_samples_seg_initialize_LOOKBACK (vb);
    
    rom field_start;
    unsigned field_len=0, num_colons=0;

    // 0 or more samples
    for (char separator=0 ; separator != '\n'; samples.repeats++) {

        field_start = next_field;
        next_field = (char *)seg_get_next_item (VB, field_start, len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, has_13, "sample-subfield");

        ASSVCF (field_len, "Error: invalid VCF file - expecting sample data for sample # %u, but found a tab character", 
                samples.repeats+1);

        vb->sample_i = samples.repeats;
        num_colons += vcf_seg_one_sample (vb, dl, ctxs, &samples, (char *)field_start, field_len);

        ASSVCF (samples.repeats < vcf_num_samples || separator == '\n',
                "invalid VCF file - expecting a newline after the last sample (sample #%u)", vcf_num_samples);
    }
    vb->sample_i = 0;
    
    ASSVCF (samples.repeats <= vcf_num_samples, "according the VCF header, there should be %u sample%s per line, but this line has %u samples - that's too many",
            vcf_num_samples, vcf_num_samples==1 ? "" : "s", samples.repeats);

    // in some real-world files I encountered have too-short lines due to human errors. we pad them
    if (samples.repeats < vcf_num_samples) {
        WARN_ONCE ("FYI: the number of samples in variant CHROM=%.*s POS=%"PRId64" is %u, different than the VCF column header line which has %u samples",
                   vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS), samples.repeats, vcf_num_samples);

        if (dl->has_haplotype_data) {
            char *ht_data = Bc (CTX(FORMAT_GT_HT)->local, vb->line_i * vb->ploidy * vcf_num_samples + vb->ploidy * samples.repeats);
            unsigned num_missing = vb->ploidy * (vcf_num_samples - samples.repeats); 
            memset (ht_data, '*', num_missing);
        }
    }
    
    // assign all translators. note: we either have translators for all translatable items, or none at all.
    if (z_is_dvcf)
        for (uint32_t i=0; i < n_items; i++)
            if (ctxs[i]->line_is_luft_trans)
                samples.items[i].translator = ctxs[i]->luft_trans;

    container_seg (vb, CTX(VCF_SAMPLES), &samples, 0, 0, samples.repeats + num_colons); // account for : and \t \r \n separators

    ctx_set_last_value (VB, CTX(VCF_SAMPLES), (ValueType){ .i = samples.repeats });

    CTX(FORMAT_GT_HT)->local.len32 = (vb->line_i+1) * vb->ht_per_line;
 
    return next_field;
}

