// ------------------------------------------------------------------
//   vcf_samples.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

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

static SmallContainer con_AD={}, con_ADALL={}, con_ADF={}, con_ADR={}, con_SAC={}, con_F1R2={}, con_F2R1={}, 
    con_MB={}, con_SB={}, con_AF={};

static char sb_snips[2][32], mb_snips[2][32], f2r1_snips[VCF_MAX_ARRAY_ITEMS][32], adr_snips[VCF_MAX_ARRAY_ITEMS][32], adf_snips[VCF_MAX_ARRAY_ITEMS][32], 
    rdf_snip[32], rdr_snip[32], adf_snip[32], adr_snip[32], ad_varscan_snip[32], ab_snip[48], gq_by_pl[50], gq_by_gp[50],
    af_snip[32], sac_snips[VCF_MAX_ARRAY_ITEMS/2][32], PL_to_PLn_redirect_snip[30], PL_to_PLy_redirect_snip[30];

static unsigned sb_snip_lens[2], mb_snip_lens[2], f2r1_snip_lens[VCF_MAX_ARRAY_ITEMS], adr_snip_lens[VCF_MAX_ARRAY_ITEMS], adf_snip_lens[VCF_MAX_ARRAY_ITEMS], 
    af_snip_len, sac_snip_lens[VCF_MAX_ARRAY_ITEMS/2], rdf_snip_len, rdr_snip_len, adf_snip_len, adr_snip_len, ad_varscan_snip_len,
    ab_snip_len, gq_by_pl_len, gq_by_gp_len, PL_to_PLn_redirect_snip_len, PL_to_PLy_redirect_snip_len;

// prepare snip of A - B
static void vcf_seg_prepare_minus_snip (DictId dict_id_a, DictId dict_id_b, char *snip, unsigned *snip_len)
{
    snip[0] = SNIP_SPECIAL;
    snip[1] = VCF_SPECIAL_MINUS;
    
    DictId two_dicts[2] = { dict_id_a, dict_id_b };
    *snip_len = 2 + base64_encode ((uint8_t *)two_dicts, sizeof (two_dicts), &snip[2]);
}

static DictId make_array_item_dict_id (uint64_t dict_id_num, unsigned item_i)
{
    const uint8_t *id = ((DictId)dict_id_num).id;
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
        vcf_seg_prepare_minus_snip (con_AD.items[i].dict_id, con_SB.items[i*2].dict_id, sb_snips[i], &sb_snip_lens[i]);
        vcf_seg_prepare_minus_snip (con_AD.items[i].dict_id, con_MB.items[i*2].dict_id, mb_snips[i], &mb_snip_lens[i]);
    }

    for (unsigned i=0; i < VCF_MAX_ARRAY_ITEMS/2; i++) 
        vcf_seg_prepare_minus_snip (con_AD.items[i].dict_id, con_SAC.items[i*2].dict_id, sac_snips[i], &sac_snip_lens[i]);

    for (unsigned i=0; i < VCF_MAX_ARRAY_ITEMS; i++) {
        vcf_seg_prepare_minus_snip (con_AD.items[i].dict_id, con_F1R2.items[i].dict_id, f2r1_snips[i], &f2r1_snip_lens[i]);
        vcf_seg_prepare_minus_snip (con_AD.items[i].dict_id, con_ADF.items[i].dict_id,  adr_snips[i],  &adr_snip_lens[i]);
        vcf_seg_prepare_minus_snip (con_AD.items[i].dict_id, con_ADR.items[i].dict_id,  adf_snips[i],  &adf_snip_lens[i]);
    }

    vcf_seg_prepare_minus_snip ((DictId)_FORMAT_RD, (DictId)_FORMAT_RDR, rdf_snip, &rdf_snip_len);
    vcf_seg_prepare_minus_snip ((DictId)_FORMAT_RD, (DictId)_FORMAT_RDF, rdr_snip, &rdr_snip_len);
    vcf_seg_prepare_minus_snip ((DictId)_FORMAT_AD, (DictId)_FORMAT_ADR, adf_snip, &adf_snip_len);
    vcf_seg_prepare_minus_snip ((DictId)_FORMAT_AD, (DictId)_FORMAT_ADF, adr_snip, &adr_snip_len);
    vcf_seg_prepare_minus_snip ((DictId)_FORMAT_DP, (DictId)_FORMAT_RD,  ad_varscan_snip,  &ad_varscan_snip_len);

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
}

void vcf_samples_seg_initialize (VBlockVCFP vb)
{
    CTX(FORMAT_ADALL)->flags.store = STORE_INT;
    CTX(FORMAT_DP)->   flags.store = STORE_INT;
    CTX(FORMAT_AD)->   flags.store = STORE_INT;   // since v13
    CTX(FORMAT_RD)->   flags.store = STORE_INT;   // since v13
    CTX(FORMAT_RDR)->  flags.store = STORE_INT;   // since v13
    CTX(FORMAT_RDF)->  flags.store = STORE_INT;   // since v13
    CTX(FORMAT_ADR)->  flags.store = STORE_INT;   // since v13
    CTX(FORMAT_ADF)->  flags.store = STORE_INT;   // since v13
    CTX(FORMAT_SDP)->  flags.store = STORE_INT;   // since v13

    CTX(FORMAT_PS)->   flags.store = STORE_INT;   
    CTX(FORMAT_PS)->   no_stons = true;   

    CTX(FORMAT_GT)->no_stons  = true; // we store the GT matrix in local, so cannot accomodate singletons
    vb->ht_matrix_ctx = CTX(FORMAT_GT_HT); // different for different data types

    // create additional contexts as needed for compressing FORMAT/GT - must be done before merge
    if (vcf_num_samples) 
        codec_pbwt_seg_init (VB, CTX(FORMAT_PBWT_RUNS), CTX(FORMAT_PBWT_FGRC));

    stats_set_consolidation (VB, FORMAT_GT,  4, FORMAT_GT_HT, FORMAT_GT_HT_INDEX, FORMAT_PBWT_RUNS, FORMAT_PBWT_FGRC);
    stats_set_consolidation (VB, FORMAT_AB,  1, FORMAT_AB3);

    // determine which way to seg PL - Mux by dosage or Mux by dosageXDP, or test both options
    CTX(FORMAT_PL)->no_stons = true;
    vb->PL_mux_by_DP = (flag.best && !z_dual_coords && !segconf.running && segconf.has_DP_before_PL) 
        ? segconf.PL_mux_by_DP // set by a previous VB in vcf_FORMAT_PL_decide or still in its initial value of PL_mux_by_DP_TEST
        : PL_mux_by_DP_NO;

    // initialize dosage multiplexers
    #define init_mux(name,store_type) seg_mux_init ((VBlockP)vb, 4, VCF_SPECIAL_MUX_BY_DOSAGE, FORMAT_##name, FORMAT_##name, (store_type), (MultiplexerP)&vb->mux_##name, "0123")
    init_mux(PRI,  STORE_NONE);
    init_mux(GL,   STORE_NONE);
    init_mux(DS,   STORE_NONE);
    init_mux(PP,   STORE_NONE);
    init_mux(GP,   STORE_NONE);
    init_mux(GQ,   STORE_NONE);
    init_mux(PVAL, STORE_NONE);
    init_mux(FREQ, STORE_NONE);
    init_mux(RD,   STORE_INT);
    init_mux(PLn,  STORE_NONE);

    seg_mux_init (VB, 3*DOSAGExDP_NUM_DPs + 1, VCF_SPECIAL_MUX_BY_DOSAGExDP, FORMAT_PLy, FORMAT_PLy, STORE_NONE, (MultiplexerP)&vb->mux_PLy, NULL);
}

//--------------------
// Multiplex by dosage
// -------------------

// returns: 0,1,2 correspond to 0/0 0/1 1/1, 3 means "no dosage" - multi-allelic or '.' or ploidy > 2
// -1 means caller should not use her SPECIAL
static inline int vcf_seg_get_mux_channel_i (VBlockVCFP vb, bool fail_if_dvcf_refalt_switch)
{
    // fail if this is a DVCF ref<>alt switch, if caller requested so,
    // or if there is no valid GT in this sample
    if (fail_if_dvcf_refalt_switch && z_dual_coords)// && LO_IS_OK_SWITCH (last_ostatus)) // TODO: no real need to fallback unless REF<>ALT switch, but this doesn't work yet
        return -1;

    // fail if there is no GT in this variant
    if (!ctx_encountered_in_line (VB, FORMAT_GT))
        return -1;

    int64_t dosage = CTX(FORMAT_GT)->last_value.i; // dosage stored here by vcf_seg_FORMAT_GT
     
    return (dosage >= 0 && dosage <= 2) ? dosage : 3; // 3 happens if sample has ploidy > 2, or if one of the alleles is not 0 or 1
}

// if cell is NULL, leaves it up to the caller to seg to the channel 
static inline ContextP vcf_seg_FORMAT_mux_by_dosage (VBlockVCF *vb, Context *ctx, STRp(cell), const DosageMultiplexer *mux) 
{
    int channel_i = vcf_seg_get_mux_channel_i (vb, true);

    // we don't use the multiplexer if its a DVCF REF<>ALT switch variant as GT changes
    if (channel_i == -1) {
        if (cell) seg_by_ctx (VB, STRa(cell), ctx, cell_len);
        return ctx;
    }

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, MUX, channel_i);

    if (cell) seg_by_ctx (VB, STRa(cell), channel_ctx, cell_len);

    // note: this is not necessarily all-the-same - there could be unmuxed snips due to REF<>ALT switch, and/or WORD_INDEX_MISSING 
    seg_by_ctx (VB, STRa(mux->snip), ctx, 0);

    return channel_ctx;
}

// mirroring seg, we accept monoploid or diploid genotypes, with alleles 0 and 1.
// there return 0/0,0->0 ; 0/1,1/0,1->1 ; 1/1->2 ; others->3
static inline int vcf_piz_get_mux_channel_i (VBlockP vb)
{
    const char *gt = last_txt (vb, FORMAT_GT);
    unsigned gt_len = CTX(FORMAT_GT)->last_txt_len ;

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
    
    ContextP channel_ctx = MCTX (channel_i, snip, snip_len);
    ASSPIZ (channel_ctx, "Cannot find channel context of channel_i=%d of multiplexed context %s", channel_i, ctx->tag_name);

    reconstruct_from_ctx (vb, channel_ctx->did_i, 0, reconstruct);

    if (ctx->flags.store == STORE_NONE) return false;

    // propagate last_value up
    new_value->i = channel_ctx->last_value.i; // note: last_value is a union, this copies the entire union
    return true; 
}

// if cell is NULL, leaves it up to the caller to seg to the channel 
static inline void vcf_seg_FORMAT_mux_by_dosagexDP (VBlockVCF *vb, Context *ctx, STRp(cell), const DosageDPMultiplexer *mux) 
{
    if (!ctx_encountered (VB, FORMAT_DP)) goto cannot_use_special; // no DP in the FORMAT of this line

    int64_t DP;
    if (!str_get_int (last_txt(vb, FORMAT_DP), vb->last_txt_len(FORMAT_DP), &DP)) // In some files, DP may be '.'
        DP=0;

    int channel_i = vcf_seg_get_mux_channel_i (vb, true); // we don't use the multiplexer if its a DVCF REF<>ALT switch variant as GT changes
    if (channel_i == -1) goto cannot_use_special;

    DP = MAX_(0, MIN_(DP, DOSAGExDP_NUM_DPs-1));
    channel_i = (channel_i == 3) ? (DOSAGExDP_NUM_DPs * 3) : (DP*3 + channel_i);

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, MUX, channel_i);

    seg_by_ctx (VB, STRa(cell), channel_ctx, cell_len);

    // note: this is not necessarily all-the-same - there could be unmuxed snips due to REF<>ALT switch, and/or WORD_INDEX_MISSING 
    seg_by_ctx (VB, STRa(mux->snip), ctx, 0);
    return;

cannot_use_special:
    seg_by_ctx (VB, STRa(cell), ctx, cell_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_DOSAGExDP)
{
    unsigned num_channels = ctx->con_cache.len ? ctx->con_cache.len : (1 + str_count_char (STRa(snip), '\t'));
    unsigned num_dps = num_channels / 3;

    int64_t DP = ctx_has_value (vb, FORMAT_DP) ? CTX(FORMAT_DP)->last_value.i : 0; // possibly "."
    DP = MAX_(0, MIN_(DP, num_dps-1));

    int channel_i = vcf_piz_get_mux_channel_i (vb); 
    channel_i = (channel_i == 3) ? (num_dps * 3) : (DP*3 + channel_i);

    ContextP channel_ctx = MCTX (channel_i, snip, snip_len);
    ASSPIZ (channel_ctx, "Cannot find channel context of channel_i=%d of multiplexed context %s", channel_i, ctx->tag_name);

    reconstruct_from_ctx (vb, channel_ctx->did_i, 0, reconstruct);

    if (ctx->flags.store == STORE_NONE) return false;

    // propagate last_value up
    new_value->i = channel_ctx->last_value.i; // note: last_value is a union, this copies the entire union
    return true; 
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
static inline WordIndex vcf_seg_FORMAT_transposed (VBlockVCF *vb, Context *ctx, STRp(cell), unsigned add_bytes)
{
    ctx->ltype = LT_UINT32_TR;
    ctx->flags.store = STORE_INT;
    
    buf_alloc (vb, &ctx->local, 1, vb->lines.len * vcf_num_samples, uint32_t, 1, "contexts->local");

    if (cell_len == 1 && cell[0] == '.') {
        NEXTENT (uint32_t, ctx->local) = 0xffffffff;
    }
    else {
        ASSSEG (str_get_int (STRa(cell), &ctx->last_value.i) && ctx->last_value.i >= 0 && ctx->last_value.i <= 0xfffffffe, 
                cell, "While compressing %s expecting an integer in the range [0, 0xfffffffe] or a '.', but found: %.*s", 
                ctx->tag_name, cell_len, cell);

        NEXTENT (uint32_t, ctx->local) = (uint32_t)ctx->last_value.i;
    }

    // add a LOOKUP to b250
    seg_by_ctx (VB, (char []){ SNIP_LOOKUP }, 1, ctx, add_bytes);

    return 0;
}

// a comma-separated array - each element goes into its own item context, single repeat
static WordIndex vcf_seg_FORMAT_A_R (VBlockVCF *vb, Context *ctx, SmallContainer con /* by value */, STRp(value), StoreType item_store_type,
                                     void (*seg_item_cb)(VBlockVCFP, Context *ctx, unsigned num_items, ContextP *item_ctxs, 
                                                         const char**, const uint32_t*, const int64_t*))
{   
    str_split (value, value_len, VCF_MAX_ARRAY_ITEMS, ',', item, false);
    
    if (!(con.nitems_lo = n_items)) 
        return seg_by_ctx (VB, value, value_len, ctx, value_len); // too many items - normal seg

    Context *item_ctxs[con.nitems_lo];
    int64_t values[con.nitems_lo];

    for (unsigned i=0; i < con.nitems_lo; i++) {

        item_ctxs[i] = ctx_get_ctx (vb, con.items[i].dict_id);
        
        if (!item_ctxs[i]->seg_initialized) {
            item_ctxs[i]->flags.store = item_store_type;
            item_ctxs[i]->seg_initialized = true;
            stats_set_consolidation (VB, ctx->did_i, 1, item_ctxs[i]->did_i);
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
        seg_item_cb (vb, ctx, con.nitems_lo, item_ctxs, items, item_lens, values);

    // case: seg items as normal snips
    else 
        for (unsigned i=0; i < con.nitems_lo; i++) 
            seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);

    ctx->last_txt_len = con.nitems_lo; // seg only: for use by vcf_seg_*_items callbacks
    
    return container_seg (vb, ctx, (ContainerP)&con, 0, 0, con.nitems_lo-1); // account for the commas
}

//----------
// FORMAT/AD
// ---------

// <ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
static inline WordIndex vcf_seg_FORMAT_AD_varscan (VBlockVCF *vb, Context *ctx, STRp(ad_str))
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
static void vcf_seg_AD_items (VBlockVCFP vb, Context *ctx, unsigned num_items, ContextP *item_ctxs, 
                              STRps(item), const int64_t *values)
{       
    bool has_adall_this_sample = segconf.vcf_has_ADALL && ctx_encountered (VB, FORMAT_ADALL); // note: we can delta vs ADALL unless segconf says so, bc it can ruin other fields that rely on peeking AD, eg AB
    int64_t sum = 0; 

    for (unsigned i=0; i < num_items; i++) {

        // If we have ADALL in this file, we delta vs ADALL if we have it in this sample, or seg normally if not
        if (segconf.vcf_has_ADALL) {
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
                    seg_mux_init (VB, 4, VCF_SPECIAL_MUX_BY_DOSAGE, item_ctxs[i]->did_i, FORMAT_AD, STORE_INT, (MultiplexerP)&vb->mux_AD[i], "0123");
                
                vcf_seg_FORMAT_mux_by_dosage (vb, item_ctxs[i], STRi(item, i), &vb->mux_AD[i]);
            }
            else
                seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
        }
        sum += values[i];
    }

    // AD value is sum of its items
    ctx_set_last_value (VB, ctx, sum);

    memcpy (vb->ad_values, values, num_items * sizeof (values[0]));
}

//-------------
// FORMAT/ADALL
//-------------

// Sepcial treatment for item 0
static void vcf_seg_ADALL_items (VBlockVCFP vb, Context *ctx, unsigned num_items, ContextP *item_ctxs, 
                                 STRps(item), const int64_t *values)
{
    if (segconf.running) segconf.vcf_has_ADALL = true;

    for (unsigned i=0; i < num_items; i++) 
        if (i==0 || i==1) {
            if (!vb->mux_ADALL[i].num_channels)
                seg_mux_init (VB, 4, VCF_SPECIAL_MUX_BY_DOSAGE, item_ctxs[i]->did_i, FORMAT_ADALL, STORE_INT, (MultiplexerP)&vb->mux_ADALL[i], "0123");
            
            vcf_seg_FORMAT_mux_by_dosage (vb, item_ctxs[i], STRi(item, i), &vb->mux_ADALL[i]);
        }
        else 
            seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
}

//----------------------
// FORMAT/F2R1, ADF, ADR
//----------------------

// used when Vector is expected to be (AD-OtherVector) - if it indeed is, we use a special snip
static void vcf_seg_AD_complement_items (VBlockVCFP vb, Context *ctx, unsigned num_items, ContextP *item_ctxs, 
                                         STRps(item), const int64_t *values,
                                         DictId other_dict_id, const SmallContainer *other_con,
                                         char my_snips[][32], unsigned *my_snip_lens)
{
    // we can use the formula only if AD,F1R1 were encountered in this line, and that they have the number of items as us
    ContextP ad_ctx=CTX(FORMAT_AD), other_ctx;
    bool use_formula = ctx_encountered (VB, FORMAT_AD) &&
                       ctx_encountered_by_dict_id (VB, other_dict_id, &other_ctx) &&
                       ad_ctx->last_txt_len    == num_items &&  // last_txt_len is # of items stored by vcf_seg_FORMAT_A_R 
                       other_ctx->last_txt_len == num_items;

    for (unsigned i=0; i < num_items; i++) {

        // case: as expected, F1R2 + F2R1 = AD - seg as a F2R1 as a MINUS snip
        if (use_formula && vb->ad_values[i] == values[i] + ECTX (other_con->items[i].dict_id)->last_value.i) 
            seg_by_ctx (VB, my_snips[i], my_snip_lens[i], item_ctxs[i], item_lens[i]); 

        // case: the formula doesn't work for this item - seg a normal snip
        else
            seg_by_ctx (VB, STRi(item,i), item_ctxs[i], item_lens[i]);
    }
}

// F2R1 = AD - F1R2 (applied if AD and F1R2 are encountered before F2R1)
static void vcf_seg_F2R1_items (VBlockVCFP vb, Context *ctx, unsigned num_items, ContextP *item_ctxs, 
                                STRps(item), const int64_t *values)
{
    vcf_seg_AD_complement_items (vb, ctx, num_items, item_ctxs, items, item_lens, values, (DictId)_FORMAT_F1R2, &con_F1R2, f2r1_snips, f2r1_snip_lens);
}

// ADF = AD - ADR (applied if AD and ADR are encountered before ADF)
static void vcf_seg_ADF_items (VBlockVCFP vb, Context *ctx, unsigned num_items, ContextP *item_ctxs, 
                               STRps(item), const int64_t *values)
{
    vcf_seg_AD_complement_items (vb, ctx, num_items, item_ctxs, items, item_lens, values, (DictId)_FORMAT_ADR, &con_ADR, adf_snips, adf_snip_lens);
}

// ADR = AD - ADF (applied if AD and ADF are encountered before ADR)
static void vcf_seg_ADR_items (VBlockVCFP vb, Context *ctx, unsigned num_items, ContextP *item_ctxs, 
                               STRps(item), const int64_t *values)
{
    vcf_seg_AD_complement_items (vb, ctx, num_items, item_ctxs, items, item_lens, values, (DictId)_FORMAT_ADF, &con_ADF, adr_snips, adr_snip_lens);
}

//----------
// FORMAT/SB
//----------

// For bi-allelic SNPs, sum every of two values is expected to equal the corresponding value in AD. Example: AD=59,28 SB=34,25,17,11. 
// seg the second of every pair as a MINUS snip
static void vcf_seg_SB_items (VBlockVCFP vb, Context *ctx, unsigned num_items, ContextP *item_ctxs, 
                              STRps(item), const int64_t *values)
{
    // verify that AD was encountered in this line, and that it has exactly half the number of items as us
    ContextP ad_ctx=CTX(FORMAT_AD);
    bool use_formula = ctx_encountered (VB, FORMAT_AD) && ad_ctx->last_txt_len == 2 && num_items == 4; // note: last_txt_len = # of items stored by vcf_seg_FORMAT_A_R

    for (unsigned i=0; i < num_items; i++) {

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
static void vcf_seg_SAC_items (VBlockVCFP vb, Context *ctx, unsigned num_items, ContextP *item_ctxs, 
                               STRps(item), const int64_t *values)
{
    // verify that AD was encountered in this line, and that it has exactly half the number of items as us
    bool use_formula = ctx_encountered (VB, FORMAT_AD) && 2 * CTX(FORMAT_AD)->last_txt_len == num_items; // note: last_txt_len = # of items stored by vcf_seg_FORMAT_A_R

    for (unsigned i=0; i < num_items; i++) {

        // seg odd-numbered element as AD - (even element), if the sum is correct
        if (use_formula && i%2 && vb->ad_values[i/2] == values[i-1] + values[i]) 
            seg_by_ctx (VB, sac_snips[i/2], sac_snip_lens[i/2], item_ctxs[i], item_lens[i]); 

        else {
            item_ctxs[i]->flags.store = STORE_INT; // consumed by the odd items ^
            seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
        }
    }
}

//----------
// FORMAT/MB
//----------

// For bi-allelic SNPs: sum every of two items is expected to equal the corresponding value in AD. Example: AD=7,49 F2R1=3,28 MB=4,3,26,23 
// In addition, the even-numbered item is quite similar to the corresponding value in F2R1.
// Seg the even items as delta from F2R1 and odd items as a MINUS snip between AD and the preceding even item
static void vcf_seg_MB_items (VBlockVCFP vb, Context *ctx, unsigned num_items, ContextP *item_ctxs, 
                              STRps(item), const int64_t *values)
{
    bool use_formula_even = ctx_encountered (VB, FORMAT_F2R1) && CTX(FORMAT_F2R1)->last_txt_len == 2 && num_items == 4;
    bool use_formula_odd  = ctx_encountered (VB, FORMAT_AD)   && CTX(FORMAT_AD)  ->last_txt_len == 2 && num_items == 4; // last_txt_len is # of items set by vcf_seg_FORMAT_A_R

    for (unsigned i=0; i < num_items; i++) {

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

// parameter is two dict_id's (in base64). reconstructs dict1.last_value - dict2.last_value
SPECIAL_RECONSTRUCTOR (vcf_piz_special_MINUS)
{
    // decode and store the the contexts in the first call for ctx (only one MINUS snip allowed per ctx)
    if (!ctx->con_cache.len) {
        buf_alloc_zero (vb, &ctx->con_cache, 0, 2, ContextP, 1, "con_cache");

        DictId two_dicts[2];
        base64_decode (snip, &snip_len, (uint8_t *)two_dicts);

        *ENT(ContextP, ctx->con_cache, 0) = ECTX (two_dicts[0]);
        *ENT(ContextP, ctx->con_cache, 1) = ECTX (two_dicts[1]);
    }

    new_value->i = (*ENT(ContextP, ctx->con_cache, 0))->last_value.i - 
                   (*ENT(ContextP, ctx->con_cache, 1))->last_value.i;

    if (reconstruct)
        RECONSTRUCT_INT (new_value->i); 

    return true; // has new_value
}

//----------
// FORMAT/AF
// ---------

static inline WordIndex vcf_seg_FORMAT_AF (VBlockVCF *vb, Context *ctx, STRp(cell))
{
    if (vcf_num_samples == 1 && // very little hope that INFO/AF is equal to FORMAT/AF if we have more than one sample
        !z_dual_coords &&       // note: we can't use SNIP_COPY in dual coordinates, because when translating, it will translate the already-translated INFO/AF
        ctx_encountered_in_line (VB, INFO_AF) && 
        str_issame (cell, CTX(INFO_AF)->last_snip))
        return seg_by_ctx (VB, af_snip, af_snip_len, ctx, cell_len);
    else
        return vcf_seg_FORMAT_A_R (vb, ctx, con_AF, STRa(cell), STORE_NONE, NULL);
}

//----------
// FORMAT/PS
// ---------

static inline WordIndex vcf_seg_FORMAT_PS (VBlockVCF *vb, Context *ctx, STRp(cell))
{
    int64_t ps_value=0;
    if (str_get_int (STRa(cell), &ps_value) && ps_value == ctx->last_value.i) // same as previous line
        return seg_by_ctx (VB, ((char []){ SNIP_SELF_DELTA, '0' }), 2, ctx, cell_len);

    return seg_delta_vs_other_do (VB, ctx, CTX(VCF_POS), STRa(cell), 1000, cell_len);
}

//----------
// FORMAT/GQ
// ---------

static int value_sorter(const void *a, const void *b)  
{
    return *(int64_t *)b - *(int64_t*)a; // sort in reverse order - usually faster as GP[0] / PL[0] are usually the biggest (corresponding to GT=0/0)
}

static int64_t vcf_predict_GQ (VBlockVCFP vb, DidIType src_did_i)
{
    bool is_gp = (src_did_i == FORMAT_GP);
    STR(src);

    if (command == ZIP)
        CTXlast (src, CTX(src_did_i));
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

    seg_by_did_i (VB, STRa (snip), FORMAT_GQ, gq_len);
}

static inline void vcf_seg_FORMAT_GQ (VBlockVCF *vb)
{
    ContextP gq_ctx = CTX(FORMAT_GQ);
    STRlast (gq, gq_ctx);

    // second best: mux by dosage 
    if (!segconf.running && !segconf.GQ_by_GP && !segconf.GQ_by_PL) {
        vcf_seg_FORMAT_mux_by_dosage (vb, gq_ctx, STRa(gq), &vb->mux_GQ);
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
    const char *tab = memchr (snip, '\t', snip_len);

    ContextP src_ctx = MCTX (0, snip, tab - snip);

    int64_t prediction = vcf_predict_GQ (VB_VCF, src_ctx->did_i);
    int64_t delta = atoi (tab+1);

    RECONSTRUCT_INT (prediction - delta);

    return false; // no new value
}

//----------
// FORMAT/DP
// ---------

static inline WordIndex vcf_seg_FORMAT_DP (VBlockVCF *vb, Context *ctx, STRp(cell))
{
    seg_set_last_txt (VB, ctx, STRa(cell), STORE_NONE);

    // case - we have FORMAT/AD - calculate delta vs the sum of AD components
    if (ctx_has_value (VB, FORMAT_AD))
        return seg_delta_vs_other (VB, ctx, CTX(FORMAT_AD), STRa(cell));

    // case - we have FORMAT/SDP - calculate delta vs the sum of AD components
    else if (ctx_has_value (VB, FORMAT_SDP))
        return seg_delta_vs_other (VB, ctx, CTX(FORMAT_SDP), STRa(cell));
    
    // case: there is only one sample there is an INFO/DP too, we store a delta 
    else if (vcf_num_samples == 1 && ctx_has_value (VB, INFO_DP)) 
        return seg_delta_vs_other (VB, ctx, CTX(INFO_DP), STRa(cell));

    // case: no FORMAT/AD and no INFO/DP - store in transposed matrix
    else 
        return vcf_seg_FORMAT_transposed (vb, ctx, STRa(cell), cell_len); // this handles DP that is an integer or '.'
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

    // if we're validating a FORMAT field with --chain (in vcf_seg_validate_luft_trans_one_sample, if REF<>ALT) - accept a valid scientific notation
    // as it will be converted to normal notation in vcf_seg_one_sample
    if (validate_only && chain_is_loaded && dict_id_is_vcf_format_sf (ctx->dict_id) &&
        str_scientific_to_decimal (recon, recon_len, NULL, NULL, &f) && f >= 0.0 && f <= max_value) return true; // scientific notation in the valid range

    // if item format is inconsistent with AF being a probability value - we won't translate it
    if (!str_get_float (recon, recon_len, &f, format, NULL) || f < 0.0 || f > max_value) 
        return false;
    
    if (validate_only) return true; 

    vb->txt_data.len -= recon_len;
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
    return vcf_piz_luft_trans_complement_to_max_value (vb, ctx, recon, recon_len, validate_only, 1);
}

//----------
// FORMAT/GL
// ---------

// convert an array of probabilities to an array of integer phred scores capped at 60
static void vcf_convert_prob_to_phred (VBlockVCFP vb, const char *flag_name, STRp(snip), char *optimized_snip, unsigned *optimized_snip_len)
{
    str_split_floats (snip, snip_len, 0, ',', prob, false);
    ASSVCF (n_probs, "cannot to apply %s to value \"%.*s\"", flag_name, snip_len, snip); // not an array of floats - abort, because we already changed the FORMAT field

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
static bool vcf_phred_optimize (const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len /* in / out */)
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
static inline void vcf_seg_FORMAT_AB (VBlockVCF *vb, Context *ctx, STRp(ab))
{
    int channel_i = vcf_seg_get_mux_channel_i (vb, false);
    bool is_0_or_2 = (channel_i==0 || channel_i==2); // note: dos02 doesn't change in case of DVCF REF<>ALT switch
    bool ab_missing = ab_len==1 && *ab=='.';

    if (channel_i==-1                   || // GT didn't produce a mux channel
        (channel_i==1 && z_dual_coords) || // we can't handle channel 1 in dual coordinates (TODO: limit to REF<>ALT switch)
        segconf.vcf_has_ADALL           || // we can't handle AD0/AD1 peeking in AD is a delta vs ADALL
        (is_0_or_2 && !ab_missing)) {      // if channel is 0 or 2, we expected a '.' 
    
        seg_by_ctx (VB, STRa(ab), ctx, ab_len);
        return;
    }

    // prepare rollback data: we will verify channel 1 in vcf_seg_FORMAT_AB_verify_channel1 and rollback if necessary
    if (channel_i==1) {
        seg_set_last_txt (VB, ctx, STRa(ab), STORE_NONE);
        ctx_set_last_value (VB, ctx, (ValueType){.i = 1}); // need verification

        ctx_create_rollback_point (ctx); 
    }

    if (channel_i==3) { 
        seg_by_ctx (VB, STRa(ab_snip), ctx, 0);
        seg_by_did_i (VB, STRa(ab), FORMAT_AB3, ab_len);
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
static inline void vcf_seg_FORMAT_AB_verify_channel1 (VBlockVCF *vb)
{
    ContextP ab_ctx  = CTX(FORMAT_AB);
    ContextP ad0_ctx = ECTX (con_AD.items[0].dict_id);
    ContextP ad1_ctx = ECTX (con_AD.items[1].dict_id);

    const char *ab_str = last_txtx(vb, ab_ctx);
    unsigned ab_str_len = vb->last_txt_len(FORMAT_AB);

    // rollback if we don't have AD0, AD1 this line, or if their value is not as expected by the formula
    if (!ad0_ctx || !ad1_ctx) goto rollback;

    if (!ctx_has_value (VB, ad0_ctx->did_i) || !ctx_has_value (VB, ad1_ctx->did_i)) goto rollback;

    double ad0 = ad0_ctx->last_value.i;
    double ad1 = ad1_ctx->last_value.i;
    if (ad0==0 && ad1==0) goto rollback; // formula would be division by zero 

    char recon_ab_str[32];
    calculate_AB (ad0, ad1, recon_ab_str);

    if (strlen (recon_ab_str) != ab_str_len || memcmp (recon_ab_str, ab_str, ab_str_len)) goto rollback;
    
    return; // verified

rollback:
    ctx_rollback (VB, ab_ctx);
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

    return false; // no new value
}

//----------
// FORMAT/PL
// ---------

// <ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">       
static inline void vcf_seg_FORMAT_PL (VBlockVCF *vb, Context *ctx, STRp(PL))
{
    if (segconf.running && !segconf.has_DP_before_PL) 
        segconf.has_DP_before_PL = ctx_encountered (VB, FORMAT_DP);

    seg_set_last_txt (VB, ctx, STRa(PL), STORE_NONE); // used by GQ (points into txt_data, before optimization)

    // attempt to optimize PL string, if requested
    unsigned modified_len = PL_len*2 + 10;                 
    char modified[PL_len]; // note: modifying functions need to make sure not to overflow this space

    if (flag.optimize_phred && vcf_phred_optimize (STRa(PL), modified, &modified_len)) {
        int shrinkage = (int)PL_len - (int)modified_len;
        vb->recon_size      -= shrinkage;               
        vb->recon_size_luft -= shrinkage;
        PL = modified;
        PL_len = modified_len;
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
    seg_by_did_i (VB, STRa(PL_to_PLn_redirect_snip), FORMAT_PL, 0);
}

// in PL_mux_by_DP_TEST mode, this function is called to pick one of the two compressed options, and drop the others
void vcf_FORMAT_PL_decide (VBlockVCF *vb)
{
    DidIType keep_did_i, remove_did_i;

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
        str_print_dict (ZCTX(FORMAT_PL)->dict.data, ZCTX(FORMAT_PL)->dict.len, true, false);
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
static inline WordIndex vcf_seg_FORMAT_DS (VBlockVCF *vb, Context *ctx, const char *cell, unsigned cell_len)
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
    return false; // no new value
}

// Lift-over translator for FORMAT/DS, IF it is bi-allelic and we have a ALT<>REF switch.
// We change the value to (ploidy-value)
// returns true if successful (return value used only if validate_only)
TRANSLATOR_FUNC (vcf_piz_luft_PLOIDY)
{
    if (!ctx_encountered (VB, FORMAT_GT)) return false; // we can't translate unless this variant as GT

    // use gt_prev_ploidy: in Seg, set by vcf_seg_FORMAT_GT, in validate and piz set by vcf_piz_luft_GT 
    return vcf_piz_luft_trans_complement_to_max_value (vb, ctx, recon, recon_len, validate_only, VB_VCF->gt_prev_ploidy);
}

//----------
// FORMAT/GT
// ---------

// complete haplotypes of lines that don't have GT, if any line in the vblock does have GT.
// In this case, the haplotype matrix must include the lines without GT too
void vcf_seg_FORMAT_GT_complete_missing_lines (VBlockVCF *vb)
{
    buf_alloc (vb, &CTX(FORMAT_GT_HT)->local, 0, vb->lines.len * vb->ht_per_line, char, CTX_GROWTH, "contexts->local");

    for (vb->line_i=0; vb->line_i < (uint32_t)vb->lines.len; vb->line_i++) {

        if (CTX(FORMAT_GT_HT) && !DATA_LINE (vb->line_i)->has_haplotype_data) {
            char *ht_data = ENT (char, CTX(FORMAT_GT_HT)->local, vb->line_i * vb->ht_per_line);
            memset (ht_data, '*', vb->ht_per_line);

            // NOTE: we DONT set dl->has_haplotype_data to true bc downstream we still
            // count this row as having no GT field when analyzing gt data
        }
    }

    CTX(FORMAT_GT_HT)->local.len = vb->lines.len * vb->ht_per_line;
}

// increase ploidy of the previous lines, if higher ploidy was encountered
static void vcf_seg_FORMAT_GT_increase_ploidy (VBlockVCF *vb, unsigned new_ploidy, uint32_t max_new_size)
{
    // protect against highly unlikely case that we don't have enough consumed txt data to store increased-ploidy ht data 
    ASSVCF (new_ploidy * vb->line_i * vcf_num_samples <= max_new_size, 
            "haplotype data overflow due to increased ploidy on line %"PRIu64, vb->line_i);

    uint32_t num_samples = vb->line_i * vcf_num_samples + vb->sample_i; // all samples in previous lines + previous samples in current line
    char *ht_data = FIRSTENT (char, CTX(FORMAT_GT_HT)->local);

    // copy the haplotypes backwards (to avoid overlap), padding with '*' (which are NOT counted in .repeats of the GT container)
    for (int sam_i = num_samples-1; sam_i >= 0; sam_i--) {

        int ht_i=new_ploidy-1 ; for (; ht_i >= vb->ploidy; ht_i--) 
            ht_data[sam_i * new_ploidy + ht_i] = '*'; 

        for (; ht_i >= 0; ht_i--)
            ht_data[sam_i * new_ploidy + ht_i] = ht_data[sam_i * vb->ploidy + ht_i];
    }

    vb->ploidy = new_ploidy;
    vb->ht_per_line = vb->ploidy * vcf_num_samples;
}

static inline WordIndex vcf_seg_FORMAT_GT (VBlockVCF *vb, Context *ctx, ZipDataLineVCF *dl, STRp(cell))
{
    // the GT field is represented as a Container, with a single item repeating as required by poidy, and the separator 
    // determined by the phase
    MiniContainer gt = { .repeats = 1, 
                         .nitems_lo = 1, 
                         .drop_final_repeat_sep = true, 
                         .callback = (vb->use_special_sf == USE_SF_YES),
                         .items = { { .dict_id = (DictId)_FORMAT_GT_HT } },
                       };

    unsigned save_cell_len = cell_len;

    // update repeats according to ploidy, and separator according to phase
    for (unsigned i=1; i<cell_len-1; i++)
        if (cell[i] == '|' || cell[i] == '/') {
            gt.repeats++;
            gt.repsep[0] = cell[i];
        }

    ASSVCF (gt.repeats <= VCF_MAX_PLOIDY, "ploidy=%u exceeds the maximum of %u", gt.repeats, VCF_MAX_PLOIDY);
    
    // if the ploidy of this line is bigger than the ploidy of the data in this VB so far, then
    // we have to increase ploidy of all the haplotypes read in in this VB so far. This can happen for example in 
    // the X chromosome if initial samples are male with ploidy=1 and then a female sample with ploidy=2
    if (vb->ploidy && gt.repeats > vb->ploidy) 
        vcf_seg_FORMAT_GT_increase_ploidy (vb, gt.repeats, ENTNUM (vb->txt_data, cell));

    if (!vb->ploidy) {
        vb->ploidy = gt.repeats; // very first sample in the vb
        vb->ht_per_line = vb->ploidy * vcf_num_samples;
    }

    buf_alloc (vb, &CTX(FORMAT_GT_HT)->local, vb->ploidy, vb->ht_per_line * vb->lines.len, char, CTX_GROWTH, "contexts->local");

    // note - ploidy of this sample might be smaller than vb->ploidy (eg a male sample in an X chromosesome that was preceded by a female sample, or "." sample)
    Allele *ht_data = ENT (Allele, CTX(FORMAT_GT_HT)->local, vb->line_i * vb->ht_per_line + vb->ploidy * vb->sample_i);

    int64_t dosage=0; // sum of allele values
    for (unsigned ht_i=0; ht_i < gt.repeats; ht_i++) {

        Allele ht = *(cell++); 
        cell_len--;

        ASSVCF (IS_DIGIT(ht) || ht == '.', 
                "invalid VCF file - expecting an allele in a sample to be a number 0-99 or . , but seeing %c (ht_i=%u)", ht, ht_i);

        // single-digit allele numbers
        ht_data[ht_i] = ht;

        // calculate dosage contribution of this ht (to be used in vcf_seg_FORMAT_mux_by_dosage)
        // note: only set for ploidy=1 or 2, and only if GT has 0 or 1, so values are: 0/0->0 ; 0/1->1 ; 1/0->1 ; 1/1->2 ; other->-1
        if (dosage >= 0 && gt.repeats <= 2 && (ht == '0' || ht == '1'))
            dosage += ht - '0'; // dosage only works if alleles are 0 or 1
        else
            dosage = -1; // no dosage
    
        if (!cell_len) break;

        // handle 2-digit allele numbers
        if (ht != '.' && IS_DIGIT (*cell)) {
            unsigned allele = 10 * (ht-'0') + (*(cell++) - '0');
            cell_len--;

            // make sure there isn't a 3rd digit
            ASSVCF (!cell_len || !IS_DIGIT (*cell), "VCF file sample %u - genozip currently supports only alleles up to 99", vb->sample_i+1);

            ht_data[ht_i] = '0' + allele; // use ascii 48->147

            dosage = -1; // no dosage (since allele is not 0 or 1)
        }

        // read and verify phase
        if (gt.repeats > 1 && ht_i < gt.repeats-1) {
            
            char phase = *(cell++);
            cell_len--;

            ASSVCF (phase != ' ', "invalid VCF file - expecting a tab or newline after sample %u but seeing a space", vb->sample_i+1);
            ASSVCF (phase == gt.repsep[0], "invalid VCF file -  unable to parse sample %u: expecting a %c but seeing %c", vb->sample_i+1, gt.repsep[0], phase);
        }
    } // for characters in a sample

    // if the ploidy of the sample is lower than vb->ploidy, set missing ht as '-' (which will cause deletion of themselves and their separator)
    // and set the ploidy to vb->ploidy - to avoid increase in entroy of GT.b250
    if (gt.repeats != vb->ploidy) {
        
        for (unsigned ht_i=gt.repeats; ht_i < vb->ploidy; ht_i++) 
            ht_data[ht_i] = '-'; // unlike '*', we DO count '-' in .repeats (so that we can have the same number of repeats = lower entroy in GT.b250)

        gt.repeats = vb->ploidy;
        if (!gt.repsep[0]) gt.repsep[0] = vb->gt_prev_phase; // this happens in case if a 1-ploid sample
    }

    // if this sample is a "./." - replace it with "%|%" or "%/%" according to the previous sample's phase -  
    // so that the gt container is likely identical and we reduce GT.b250 entropy. Reason: many tools
    // (including bcftools merge) produce "./." for missing samples even if all other samples are phased
    if (ht_data[0]=='.' && gt.repeats==2 && ht_data[1]=='.' && gt.repsep[0]=='/') {
        gt.repsep[0] = vb->gt_prev_phase ? vb->gt_prev_phase : '|'; // '|' is arbitrary
        ht_data[0] = ht_data[1] = '%';
    }

    // in case we have INFO/SF, we verify that it is indeed the list of samples for which the first ht is not '.'
    if (vb->use_special_sf == USE_SF_YES && ht_data[0] != '.') 
        vcf_seg_INFO_SF_one_sample (vb);

    ctx_set_last_value (VB, ctx, dosage); // to be used in vcf_seg_FORMAT_mux_by_dosage
    
    if (segconf.running) 
        segconf.count_dosage[dosage >= 0 && dosage <= 2]++;
    
    ASSVCF (!cell_len, "Invalid GT data in sample_i=%u", vb->sample_i+1);

    // shortcut if we have the same ploidy and phase as previous GT (saves re-genetrating base64 in container_seg)
    if (gt.repeats == vb->gt_prev_ploidy && gt.repsep[0] == vb->gt_prev_phase) 
        return seg_duplicate_last (VB, ctx, save_cell_len);

    else {
        vb->gt_prev_ploidy = gt.repeats;
        vb->gt_prev_phase  = gt.repsep[0];
        return container_seg (vb, ctx, (ContainerP)&gt, 0, 0, save_cell_len); 
    }
}

// Lift-over translator assigned to a FORMAT/GT item, IF it is bi-allelic and we have a ALT<>REF switch. No limitations on ploidy.
// We switch 0<>1. If its unphased (only /) - we list the 0s first, then the 1s
TRANSLATOR_FUNC (vcf_piz_luft_GT)
{
    // validate. make sure this is a bi-allelic genotype (no ploidy limitation)
    for (uint32_t i=0; i < recon_len; i += 2)  
        if (recon[i] != '0' && recon[i] != '1' && recon[i] != '.') return false;

    for (uint32_t i=1; i < recon_len; i += 2)  
        if (recon[i] != '/' && recon[i] != '|') return false;

    VB_VCF->gt_prev_ploidy = (recon_len+1) / 2; // consumed by vcf_piz_luft_PLOIDY

    if (validate_only) return true;

    // exchange 0 <> 1
    for (uint32_t i=0; i < recon_len; i += 2)
        if      (recon[i] == '0') recon[i] = '1';
        else if (recon[i] == '1') recon[i] = '0';

    return true;    
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
    
    vb->txt_data.len -= recon_len;
    
    if (num_items==2 || num_items == 3) {
        RECONSTRUCT_SEP (items[num_items-1], item_lens[num_items-1], ',');

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
// 'R' : We switch between the two comma-separated values.
// 'R2': We switch between the two PAIRS of comma-separated values.
// 'G' : We have 3 values which represent the genotypes REF/REF,REF/ALT,ALT/ALT We switch between the 1st and 3rd value.
// returns true if successful 
TRANSLATOR_FUNC (vcf_piz_luft_R)  { return vcf_piz_luft_switch_first_last (vb, ctx, recon, recon_len, 2, 'R', validate_only); } // 2 bc we only handle bi-allelic
TRANSLATOR_FUNC (vcf_piz_luft_R2) { return vcf_piz_luft_switch_first_last (vb, ctx, recon, recon_len, 4, '.', validate_only); } // 4 bc we only handle bi-allelic

TRANSLATOR_FUNC (vcf_piz_luft_G)  
{ 
    if (IS_TRIVAL_FORMAT_SUBFIELD) return true; // This is FORMAT field which is empty or "." - all good

    unsigned num_values = str_count_char (recon, recon_len, ',')+1;
    if (num_values != 3 && num_values != 2) return false; // Genozip currently only support haploid (2 bi-allelic genotypes) and diploid (3 bi-allelic genotypes) 

    return vcf_piz_luft_switch_first_last (vb, ctx, recon, recon_len, num_values, 'G', validate_only); 
}

//------------------------------------------------------------------------
// Validate that ALL subfields in ALL samples can luft-translate as needed
//------------------------------------------------------------------------

static const char *error_format_field (unsigned n_items, ContextP *ctxs)
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
static inline Context *vcf_seg_validate_luft_trans_one_sample (VBlockVCF *vb, ContextP *ctxs, uint32_t num_items, char *sample, unsigned sample_len)
{
    str_split (sample, sample_len, num_items, ':', item, false);
    ASSVCF (n_items, "Sample %u has too many subfields - FORMAT field \"%s\" specifies only %u: \"%.*s\"", 
            vb->sample_i+1, error_format_field (num_items, ctxs), num_items, sample_len, sample);

    ContextP failed_ctx = NULL; // optimistic initialization - nothing failed

    uint32_t save_ploidy = vb->gt_prev_ploidy; // ruined by vcf_piz_luft_GT 
    
    for (unsigned i=0; i < n_items; i++) {
        if (needs_translation (ctxs[i]) && item_lens[i]) {
            if ((vb->line_coords == DC_LUFT && !vcf_lo_seg_cross_render_to_primary (vb, ctxs[i], STRi(item,i), NULL, NULL)) ||
                (vb->line_coords == DC_PRIMARY && !(DT_FUNC(vb, translator)[ctxs[i]->luft_trans](VB, ctxs[i], (char *)STRi(item,i), 0, true)))) {
                failed_ctx = ctxs[i];  // failed translation
                break;
            }
        }
        ctx_set_encountered (VB, ctxs[i]); // might be needed for validation 
    }

    // reset modified values, in preparation for real Seg
    vb->gt_prev_ploidy = save_ploidy;
    for (unsigned i=0; i < n_items; i++) 
        ctx_unset_encountered (VB, ctxs[i]); 

    return failed_ctx; 
}

// If ALL subfields in ALL samples can luft-translate as required: 1.sets ctx->line_is_luft_trans for all contexts 2.lifted-back if this is a LUFT lne
// if NOT: ctx->line_is_luft_trans=false for all contexts, line is rejects (LO_FORMAT), and keeps samples in their original LUFT or PRIMARY coordinates.
static inline void vcf_seg_validate_luft_trans_all_samples (VBlockVCF *vb, uint32_t num_items, ContextP *ctxs, 
                                                            int32_t len, char *samples_start,
                                                            const char *backup_luft_samples, uint32_t backup_luft_samples_len)
{
    const char *field_start, *next_field = samples_start;
    unsigned field_len=0;
    bool has_13;

    // initialize optimistically. we will roll back and set to false if ANY subfield in ANY sample fails to translate, and re-seg all samples
    for (unsigned sf_i=0; sf_i < num_items; sf_i++)
        ctxs[sf_i]->line_is_luft_trans = needs_translation (ctxs[sf_i]); 

    // 0 or more samples
    vb->sample_i=0;
    for (char separator=0 ; separator != '\n'; vb->sample_i++) {

        field_start = next_field;
        next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, &has_13, "sample-subfield");
        ASSVCF (field_len, "unexpected tab character after sample # %u", vb->sample_i);

        Context *failed_ctx = vcf_seg_validate_luft_trans_one_sample (vb, ctxs, num_items, (char *)field_start, field_len);
        if (failed_ctx) { // some context doesn't luft-translate as required
            REJECT_SUBFIELD (LO_FORMAT, failed_ctx, ".\tCannot cross-render sample due to field %s: \"%.*s\"", failed_ctx->tag_name, field_len, field_start);

            // make all contexts untranslateable in this line
            for (unsigned i=0; i < num_items; i++)  // iterate on the order as in the line
                ctxs[i]->line_is_luft_trans = false;

            // if this is an untranslatable LUFT-only line, recover the original LUFT-coordinates samples
            if (vb->line_coords == DC_LUFT) 
                memcpy (samples_start, backup_luft_samples, backup_luft_samples_len);
        }
    }
}

// ----------
// One sample
// ----------

// returns the number of colons in the sample
static inline unsigned vcf_seg_one_sample (VBlockVCF *vb, ZipDataLineVCF *dl, ContextP *ctxs, ContainerP samples, STRp(sample))
{
    str_split (sample, sample_len, con_nitems (*samples), ':', sf, false);

    ASSVCF (n_sfs, "Sample %u has too many subfields - FORMAT field \"%s\" specifies only %u: \"%.*s\"", 
            vb->sample_i+1, error_format_field (con_nitems (*samples), ctxs), con_nitems (*samples), sample_len, sample);

    for (unsigned i=0; i < n_sfs; i++) { 

        DictId dict_id = samples->items[i].dict_id;
        Context *ctx = ctxs[i];

        unsigned modified_len = sf_lens[i]*2 + 10;                 
        char modified[modified_len]; // note: modifying functions need to make sure not to overflow this space

        #define SEG_OPTIMIZED_MUX_BY_DOSAGE(tag) ({                     \
            vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRa(modified), &vb->mux_##tag);  \
            int32_t shrinkage = (int)sf_lens[i] - (int)modified_len;    \
            vb->recon_size      -= shrinkage;                           \
            vb->recon_size_luft -= shrinkage; })
               
        // --chain: if this is RendAlg=A_1 and RendAlg=PLOIDY subfield, convert a eg 4.31e-03 to e.g. 0.00431. This is to
        // ensure primary->luft->primary is lossless (4.31e-03 cannot be converted losslessly as we can't preserve format info)
        if (chain_is_loaded && (ctx->luft_trans == VCF2VCF_A_1 || ctx->luft_trans == VCF2VCF_PLOIDY) && 
            str_scientific_to_decimal (STRi(sf, i), modified, &modified_len, NULL)) {
            
            int32_t shrinkage = (int32_t)sf_lens[i] - (int32_t)modified_len; // possibly negative = growth
            vb->recon_size      -= shrinkage; 
            vb->recon_size_luft -= shrinkage; 
            sfs[i] = modified; 
            sf_lens[i] = modified_len; 
        }

        if (!sf_lens[i])
            seg_by_ctx (VB, "", 0, ctx, 0); // generates WORD_INDEX_EMPTY

        else switch (dict_id.num) {

        // <ID=GT,Number=1,Type=String,Description="Genotype">
        case _FORMAT_GT: vcf_seg_FORMAT_GT (vb, ctx, dl, STRi(sf, i)); break;

        // <ID=GL,Number=.,Type=Float,Description="Genotype Likelihoods">
        case _FORMAT_GL:
            // --GL-to-PL:  GL: 0.00,-0.60,-8.40 -> PL: 0,6,60
            // note: we changed the FORMAT field GL->PL in vcf_seg_format_field. data is still stored in the GL context.
            // <ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
            if (flag.GL_to_PL) {
                vcf_convert_prob_to_phred (vb, "--GL-to-PL", STRi(sf, i), modified, &modified_len);
                SEG_OPTIMIZED_MUX_BY_DOSAGE(GL);
            }
            else
                vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_GL);
            break;

        // note: GP and PL - for non-optimized, I tested segging as A_R and seg_array - they are worse or not better than the default. likely because the values are correlated.
        case _FORMAT_GP:
            // convert GP (probabilities) to PP (phred values). PP was introduced in VCF v4.3.
            if (flag.GP_to_PP && vb->vcf_version >= VCF_v4_3) {
                vcf_convert_prob_to_phred (vb, "--GP-to-PP", STRi(sf, i), modified, &modified_len);
                SEG_OPTIMIZED_MUX_BY_DOSAGE(GP);
            }
            else if (flag.optimize_phred && vb->vcf_version <= VCF_v4_2 &&
                     vcf_phred_optimize (STRi(sf, i), modified, &modified_len)) 
                SEG_OPTIMIZED_MUX_BY_DOSAGE(GP);
            else
                vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_GP);

            seg_set_last_txt (VB, ctx, STRi(sf, i), STORE_NONE); // used by GQ
            break;

        // <ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">       
        case _FORMAT_PL: vcf_seg_FORMAT_PL (vb, ctx, STRi (sf, i)); break;

        // <ID=PP,Number=G,Type=Integer,Description="Phred-scaled genotype posterior probabilities rounded to the closest integer">
        case _FORMAT_PP:
            if (flag.optimize_phred && vcf_phred_optimize (STRi(sf, i), modified, &modified_len)) 
                SEG_OPTIMIZED_MUX_BY_DOSAGE(PP);
            else
                vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_PP);
            break;
        
        // <ID=PRI,Number=G,Type=Float,Description="Phred-scaled prior probabilities for genotypes">
        case _FORMAT_PRI:
            if (flag.optimize_phred && vcf_phred_optimize (STRi(sf, i), modified, &modified_len)) 
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
        
        // case: PS ("Phase Set") - might be the same as POS (for example, if set by Whatshap: https://whatshap.readthedocs.io/en/latest/guide.html#features-and-limitations)
        // or might be the same as the previous line
        case _FORMAT_PS   : vcf_seg_FORMAT_PS (vb, ctx, STRi(sf, i)); break;

        // standard: <ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        // GIAB: <ID=GQ,Number=1,Type=Integer,Description="Net Genotype quality across all datasets, calculated from GQ scores of callsets supporting the consensus GT, using only one callset from each dataset">   
        case _FORMAT_GQ   : seg_set_last_txt (VB, ctx, STRi(sf, i), STORE_NONE); break; // postpone to later
            
        // <ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        case _FORMAT_DP   : vcf_seg_FORMAT_DP (vb, ctx, STRi(sf, i)); break;
            
        // <ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
        // case: MIN_DP - it is slightly smaller and usually equal to DP - we store MIN_DP as the delta DP-MIN_DP
        // note: the delta is vs. the DP field that preceeds MIN_DP - we take the DP as 0 there is no DP that preceeds
        case _FORMAT_MIN_DP :
            if (ctx_has_value (VB, FORMAT_DP)) 
                seg_delta_vs_other (VB, ctx, CTX(FORMAT_DP), STRi(sf, i));
            else goto fallback;
            break;

        case _FORMAT_SDP :
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

        // GIAB: <ID=ADALL,Number=R,Type=Integer,Description="Net allele depths across all datasets">
        case _FORMAT_ADALL : vcf_seg_FORMAT_A_R (vb, ctx, con_ADALL, STRi(sf, i), STORE_INT, vcf_seg_ADALL_items); break;

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
        case _FORMAT_MB    : vcf_seg_FORMAT_A_R (vb, ctx, con_MB, STRi(sf, i), STORE_NONE, vcf_seg_MB_items); break;

        // <ID=SAC,Number=.,Type=Integer,Description="Number of reads on the forward and reverse strand supporting each allele (including reference)">
        case _FORMAT_SAC   : vcf_seg_FORMAT_A_R (vb, ctx, con_SAC, STRi(sf, i), STORE_NONE, vcf_seg_SAC_items); break;

        // VarScan: <ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
        case _FORMAT_RDF   : vcf_seg_FORMAT_minus (vb, ctx, STRi(sf, i), 0, CTX(FORMAT_RD), CTX(FORMAT_RDR), STRa(rdf_snip)); break;

        // VarScan: <ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
        case _FORMAT_RDR   : vcf_seg_FORMAT_minus (vb, ctx, STRi(sf, i), 0, CTX(FORMAT_RD), CTX(FORMAT_RDF), STRa(rdr_snip)); break;

        // <ID=AB,Number=1,Type=Float,Description="Allele balance for each het genotype">
        case _FORMAT_AB    : vcf_seg_FORMAT_AB (vb, ctx, STRi(sf, i)); break;

        default            :
        fallback           : seg_by_ctx (VB, STRi(sf, i), ctx, sf_lens[i]);
        }

        int64_t value;
        if (ctx->flags.store == STORE_INT && !ctx_has_value(VB, ctx->did_i) &&  // not already set
            str_get_int (STRi(sf, i), &value))
            ctx_set_last_value(VB, ctx, value);
        else        
            ctx_set_encountered (VB, ctx);
    }

    // missing subfields - defined in FORMAT but missing (not merely empty) in sample
    for (unsigned i=n_sfs; i < con_nitems (*samples); i++)  
        seg_by_ctx (VB, NULL, 0, ctxs[i], 0); // generates WORD_INDEX_MISSING

    // verify AB if its channel 1
    if (ctx_has_value (VB, FORMAT_AB))
        vcf_seg_FORMAT_AB_verify_channel1 (vb);

    // finally seg GQ if we have it
    if (ctx_encountered (VB, FORMAT_GQ))
        vcf_seg_FORMAT_GQ (vb);

    return n_sfs - 1; // number of colons
}

//------------
// All samples
//------------

const char *vcf_seg_samples (VBlockVCF *vb, ZipDataLineVCF *dl, int32_t *len, char *next_field, bool *has_13,
                             const char *backup_luft_samples, uint32_t backup_luft_samples_len)
{
    // Container for samples - we have:
    // - repeats as the number of samples in the line (<= vcf_num_samples)
    // - num_items as the number of FORMAT subfields (inc. GT)

    Container samples = *ENT (Container, vb->format_mapper_buf, dl->format_node_i); // make a copy of the template
    ContextP *ctxs = ENT (ContextP, vb->format_contexts, dl->format_node_i * MAX_FIELDS);
    uint32_t num_items = con_nitems (samples);

    // check that all subfields in all samples can be luft-translated as required, or make this a LUFT-only / PRIMARY-only line.
    // Also, if the data is in LUFT coordinates and is indeed translatable, then this lifts-back the samples to PRIMARY coordinates
    if (z_dual_coords && LO_IS_OK (last_ostatus))
        vcf_seg_validate_luft_trans_all_samples (vb, num_items, ctxs, *len, next_field, backup_luft_samples, backup_luft_samples_len);

    const char *field_start;
    unsigned field_len=0, num_colons=0;

    // 0 or more samples
    for (char separator=0 ; separator != '\n'; samples.repeats++) {

        field_start = next_field;
        next_field = (char *)seg_get_next_item (vb, field_start, len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, has_13, "sample-subfield");

        ASSVCF (field_len, "Error: invalid VCF file - expecting sample data for sample # %u, but found a tab character", 
                samples.repeats+1);

        vb->sample_i = samples.repeats;
        num_colons += vcf_seg_one_sample (vb, dl, ctxs, &samples, (char *)field_start, field_len);

        ASSVCF (samples.repeats < vcf_num_samples || separator == '\n',
                "invalid VCF file - expecting a newline after the last sample (sample #%u)", vcf_num_samples);
    }

    ASSVCF (samples.repeats <= vcf_num_samples, "according the VCF header, there should be %u sample%s per line, but this line has %u samples - that's too many",
            vcf_num_samples, vcf_num_samples==1 ? "" : "s", samples.repeats);

    // in some real-world files I encountered have too-short lines due to human errors. we pad them
    if (samples.repeats < vcf_num_samples) {
        WARN_ONCE ("FYI: the number of samples in variant CHROM=%.*s POS=%"PRId64" is %u, different than the VCF column header line which has %u samples",
                   vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS), samples.repeats, vcf_num_samples);

        if (dl->has_haplotype_data) {
            char *ht_data = ENT (char, CTX(FORMAT_GT_HT)->local, vb->line_i * vb->ploidy * vcf_num_samples + vb->ploidy * samples.repeats);
            unsigned num_missing = vb->ploidy * (vcf_num_samples - samples.repeats); 
            memset (ht_data, '*', num_missing);
        }
    }
    
    // assign all translators. note: we either have translators for all translatable items, or none at all.
    if (z_dual_coords)
        for (uint32_t i=0; i < num_items; i++)
            if (ctxs[i]->line_is_luft_trans)
                samples.items[i].translator = ctxs[i]->luft_trans;

    container_seg (vb, CTX(VCF_SAMPLES), &samples, 0, 0, samples.repeats + num_colons); // account for : and \t \r \n separators

    CTX(FORMAT_GT_HT)->local.len = (vb->line_i+1) * vb->ht_per_line;
 
    return next_field;
}

