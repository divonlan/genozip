// ------------------------------------------------------------------
//   vcf_samples.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "vcf_private.h"
#include "codec.h"
#include "base64.h"
#include "stats.h"
#include "zfile.h"
#include "zip.h"
#include "lookback.h"
#include "tip.h"
#include "zip_dyn_int.h"

static SmallContainer con_AD={}, con_ADALL={}, con_ADF={}, con_ADR={}, con_SAC={}, con_F1R2={}, con_F2R1={}, 
    con_MB={}, con_SB={}, con_AF={};

sSTRl_ARRAY(sb_snip, 2, 32);
sSTRl_ARRAY(mb_snip, 2, 32);
sSTRl_ARRAY(f2r1_snip, VCF_MAX_ARRAY_ITEMS, 32);
sSTRl_ARRAY(adr_snip,  VCF_MAX_ARRAY_ITEMS, 32);
sSTRl_ARRAY(adf_snip,  VCF_MAX_ARRAY_ITEMS, 32);
sSTRl_ARRAY(sac_snip,  VCF_MAX_ARRAY_ITEMS/2, 32);
sSTRl(PL_to_PLn_redirect_snip, 30);
sSTRl(PL_to_PLy_redirect_snip, 30);
sSTRl(rdf_snip,32);
sSTRl(rdr_snip,32);
sSTRl(adf_snip,32);
sSTRl(adr_snip,32);
sSTRl(ad_varscan_snip,32);
sSTRl(gq_by_pl,32);
sSTRl(gq_by_gp,32);
STRl(snip_copy_af,32);

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
    DO_ONCE {
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
            
        seg_prepare_snip_other (SNIP_COPY, _INFO_AF, 0, 0, snip_copy_af);

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
}

void vcf_samples_seg_initialize (VBlockVCFP vb)
{
    ctx_set_store (VB, STORE_INT, FORMAT_ADALL, FORMAT_DP, FORMAT_AD, FORMAT_RD, FORMAT_RDR, FORMAT_RDF,
                   FORMAT_ADR, FORMAT_ADF, FORMAT_SDP, DID_EOL);

    ctx_set_dyn_int (VB, FORMAT_RD, FORMAT_GQ, FORMAT_RGQ, FORMAT_MIN_DP, FORMAT_SDP, 
                     T(segconf.vcf_is_GLIMPSE_phase, FORMAT_HS), 
                     DID_EOL);
            
    if (segconf.FMT_DP_method != FMT_DP_DEFAULT) {
        seg_mux_init (vb, FORMAT_DP, VCF_SPECIAL_MUX_FORMAT_DP, false, FORMAT_DP);
    
        seg_mux_get_channel_ctx (VB, FORMAT_DP, (MultiplexerP)&vb->mux_FORMAT_DP, 0)->dyn_transposed = true; // seg transposed for variants that don't have AD / SDP to seg against
    }

    else if (segconf.FMT_DP_method == FMT_DP_DEFAULT && vcf_num_samples > 1) 
        CTX(FORMAT_DP)->dyn_transposed = true;

    CTX(FORMAT_DP)->flags.same_line = true; // DP value, when delta's against another ctx, is always relative to the other value in the current sample (regardless of whether DP or the other value are reconstructed first)

    if (!segconf.vcf_is_varscan && vcf_num_samples > 1 && !segconf.has[FORMAT_ADALL])
        ctx_get_ctx (vb, con_AD.items[0].dict_id)->dyn_transposed = true; // first item = Σ(ADᵢ) is transposed
    
    // create additional contexts as needed for compressing FORMAT/GT - must be done before merge
    if (vcf_num_samples) 
        codec_pbwt_seg_init (VB);

    // determine which way to seg PL - Mux by dosage or Mux by dosageXDP, or test both options
    CTX(FORMAT_PL)->no_stons = true;
    vb->PL_mux_by_DP = (flag.best && !segconf.running && segconf.has_DP_before_PL) // only in --best, because it is very slow
        ? segconf.PL_mux_by_DP // set by a previous VB in vcf_FORMAT_PL_decide or still in its initial value of "unknown"
        : no;

    // initialize dosage multiplexers
    #define init_mux_by_dosage(name) seg_mux_init_(VB, FORMAT_##name, ZIP_NUM_DOSAGES_FOR_MUX, VCF_SPECIAL_MUX_BY_DOSAGE, CTX(FORMAT_##name)->no_stons, (MultiplexerP)&vb->mux_##name)
    init_mux_by_dosage(PRI);
    init_mux_by_dosage(GL);
    init_mux_by_dosage(DS);
    init_mux_by_dosage(PP);
    init_mux_by_dosage(GP);
    init_mux_by_dosage(VAF);
    init_mux_by_dosage(PVAL);
    init_mux_by_dosage(FREQ);
    init_mux_by_dosage(RD);
    init_mux_by_dosage(PLn);

    seg_mux_init (vb, FORMAT_PLy, VCF_SPECIAL_MUX_BY_DOSAGExDP, false, PLy);
    
    if (segconf.has[FORMAT_DP]) 
        seg_mux_init (vb, FORMAT_RGQ, VCF_SPECIAL_RGQ, false, RGQ);
    
    switch (segconf.FMT_GQ_method) {
        case MUX_DOSAGExDP : seg_mux_init (vb, FORMAT_GQ,  VCF_SPECIAL_MUX_BY_DOSAGExDP, false, GQ); break;
        case MUX_DOSAGE    : init_mux_by_dosage(GQ); break;
        default            : break;
    }

    // flags to send to PIZ
    vb->flags.vcf.use_null_DP_method = segconf.use_null_DP_method;
}

void vcf_samples_seg_finalize (VBlockVCFP vb)
{
    if (!segconf.running) 
        vcf_samples_seg_finalize_PS_PID(vb);
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

// returns: channel=dosage up to a maximum. -1 means caller should not use her SPECIAL
int vcf_seg_get_mux_channel_i (VBlockVCFP vb)
{
    // fail if there is no GT in this variant
    if (!ctx_encountered_in_line (VB, FORMAT_GT))
        return -1;

    int64_t dosage = CTX(FORMAT_GT)->last_value.i; // dosage stored here by vcf_seg_FORMAT_GT 
     
    return MIN_(dosage, ZIP_MAX_PLOIDY_FOR_MUX);
}

// if cell is NULL, leaves it up to the caller to seg to the channel 
ContextP vcf_seg_FORMAT_mux_by_dosage (VBlockVCFP vb, ContextP ctx, STRp(cell), const DosageMultiplexer *mux) 
{
    int channel_i = vcf_seg_get_mux_channel_i (vb);

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

static void vcf_seg_FORMAT_mux_by_dosage_int (VBlockVCFP vb, ContextP ctx, int64_t value, const DosageMultiplexer *mux, uint32_t add_bytes) 
{
    int channel_i = vcf_seg_get_mux_channel_i (vb);

    // we don't use the multiplexer if its a DVCF REF⇆ALT switch variant as GT changes
    if (channel_i == -1) {
        seg_integer_as_snip_do (VB, ctx, value, add_bytes);
        return;
    }

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, ctx->did_i, MUX, channel_i);
    
    seg_integer (VB, channel_ctx, value, true, add_bytes);

    // note: this is not necessarily all-the-same - there could be unmuxed snips due to REF⇆ALT switch, and/or WORD_INDEX_MISSING 
    seg_by_ctx (VB, STRa(mux->snip), ctx, 0);
}

// mirroring seg, we accept monoploid or diploid genotypes, with alleles 0 and 1.
int vcf_piz_get_mux_channel_i (VBlockP vb)
{
    // since 15.0.36 - count ht that are not '0' or '.'
    if (z_file->max_ploidy_for_mux)  
        return MIN_(vcf_piz_GT_get_last_dosage (vb), z_file->max_ploidy_for_mux);
    
    // up to 15.0.35 - return 0/0,0->0 ; 0/1,1/0,1->1 ; 1/1->2 ; others->3
    else { 
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
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_DOSAGE)
{
    int channel_i = vcf_piz_get_mux_channel_i (vb);
    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

// if cell is NULL, leaves it up to the caller to seg to the channel 
void vcf_seg_FORMAT_mux_by_dosagexDP (VBlockVCFP vb, ContextP ctx, STRp(cell), void *mux_p) 
{
    ConstMultiplexerP mux = (ConstMultiplexerP)mux_p;

    if (!ctx_encountered (VB, FORMAT_DP)) goto cannot_use_special; // no DP in the FORMAT of this line

    int64_t DP;
    if (!str_get_int (STRlst (FORMAT_DP), &DP)) // in some files, DP may be '.'
        DP=0;

    int channel_i = vcf_seg_get_mux_channel_i (vb); // we don't use the multiplexer if its a DVCF REF⇆ALT switch variant as GT changes
    if (channel_i == -1) goto cannot_use_special;

    unsigned num_dps = mux->num_channels / ZIP_NUM_DOSAGES_FOR_MUX;
    DP = MAX_(0, MIN_(DP, num_dps-1));
    channel_i = (DP * ZIP_NUM_DOSAGES_FOR_MUX + channel_i);

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, ctx->did_i, MUX, channel_i);

    if (channel_ctx->ltype == LT_DYN_INT)
        seg_integer_or_not (VB, channel_ctx, STRa(cell), cell_len);
    else
        seg_by_ctx (VB, STRa(cell), channel_ctx, cell_len);

    // note: this is not necessarily all-the-same - there could be unmuxed snips due to REF⇆ALT switch, and/or WORD_INDEX_MISSING 
    seg_by_ctx (VB, MUX_SNIP(mux), mux->snip_len, ctx, 0);
    return;

cannot_use_special:
    seg_by_ctx (VB, STRa(cell), ctx, cell_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_DOSAGExDP)
{
    bool new_method = (z_file->max_ploidy_for_mux > 0); // true iff version is 15.0.36 or newer
    
    unsigned num_channels = recon_multi_dict_id_get_num_dicts (ctx, STRa(snip));
    unsigned num_dosages = new_method ? (z_file->max_ploidy_for_mux + 1) : 3;
    unsigned num_dps = num_channels / num_dosages;

    rom DP_str;
    int64_t DP = reconstruct_peek (vb, CTX(FORMAT_DP), &DP_str, 0).i;
    DP = (*DP_str=='.') ? 0 : MAX_(0, MIN_(DP, num_dps-1));

    int channel_i = vcf_piz_get_mux_channel_i (vb); 

    channel_i = (channel_i == 3 && !new_method) ? (num_dps * 3) : (DP*num_dosages + channel_i);

    ContextP channel_ctx = MCTX (channel_i, snip, snip_len);
    ASSPIZ (channel_ctx, "Cannot find channel context of channel_i=%d of multiplexed context %s. snip=%s", 
            channel_i, ctx->tag_name, str_snip_ex (snip-2, snip_len+2, true).s);

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

// used for DP, D0P0, A0D - store in transposed matrix in local 
static void vcf_seg_FORMAT_transposed (VBlockVCFP vb, ContextP ctx, 
                                       STRp(cell),     // option 1
                                       uint32_t value, // option 2: used if cell=0
                                       STRp(lookup_snip),
                                       unsigned add_bytes)
{
#ifdef DEBUG
    ASSERT (ctx->dyn_transposed, "expecting dyn_transposed=true in %s", ctx->tag_name);
#endif

    if (IS_PERIOD (cell)) 
        dyn_int_append_nothing_char (VB, ctx, add_bytes);
    
    else {
        ASSSEG (!cell || str_get_uint32 (STRa(cell), &value), 
                "Expecting %s=\"%.*s\" to be an integer or '.'", ctx->tag_name, STRf(cell));

        ctx_set_last_value (VB, ctx, (int64_t)value);

        dyn_int_append (VB, ctx, value, add_bytes);
    }

    seg_by_ctx (VB, STRa(lookup_snip), ctx, 0); // note: not all-the-same, because some values might be missing with a b250 of WORD_INDEX_MISSING
}

// a comma-separated array - each element goes into its own item context, single repeat
static WordIndex vcf_seg_FORMAT_A_R (VBlockVCFP vb, ContextP ctx, SmallContainer con /* by value */, STRp(value), StoreType item_store_type,
                                     void (*seg_item_cb)(VBlockVCFP, ContextP ctx, unsigned n_items, const char**, const uint32_t*, ContextP *item_ctxs, const int64_t*))
{   
    if (IS_PERIOD(value)) 
        return seg_by_ctx (VB, ".", 1, ctx, value_len); // note: segging a '.' to ctx adds the same entropy as segging a 1-item container, so we refrain from adding entropy to the item_ctxs[0] too by not segging a container.

    str_split (value, value_len, VCF_MAX_ARRAY_ITEMS, ',', item, false);
    
    if (!(con.nitems_lo = n_items)) 
        return seg_by_ctx (VB, value, value_len, ctx, value_len); // too many items - normal seg

    ContextP item_ctxs[con.nitems_lo];
    int64_t values[con.nitems_lo];

    for (unsigned i=0; i < con.nitems_lo; i++) {

        item_ctxs[i] = ctx_get_ctx (vb, con.items[i].dict_id);
        
        if (!item_ctxs[i]->is_initialized) {
            item_ctxs[i]->flags.store = item_store_type;
            item_ctxs[i]->flags.ctx_specific_flag = ctx->flags.ctx_specific_flag; // inherit - needed for same_line
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

    // case: seg items here
    else {
        int64_t sum = 0; // since 15.0.37 we sum here too
        for (unsigned i=0; i < con.nitems_lo; i++) 
            if (item_store_type == STORE_INT) {
                if (seg_integer_or_not (VB, item_ctxs[i], STRi(item, i), item_lens[i]))
                    sum += item_ctxs[i]->last_value.i; //
            }

            else if (item_store_type == STORE_FLOAT) 
                vcf_seg_string (vb, item_ctxs[i], STRi(item, i));

            else
                seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
    
        if (ctx->flags.store == STORE_INT)
            ctx_set_last_value (VB, ctx, sum);
    }

    ctx->last_txt.len = con.nitems_lo; // seg only: for use by vcf_seg_*_items callbacks
    
    return container_seg (vb, ctx, (ContainerP)&con, 0, 0, con.nitems_lo-1); // account for the commas
}

//----------
// FORMAT/AD
// ---------

// <ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
static inline void vcf_seg_FORMAT_AD_varscan (VBlockVCFP vb, ContextP ctx, STRp(ad_str))
{
    // case: AD = DP-RD
    int64_t ad;
    if (ctx_has_value (VB, FORMAT_DP) &&
        ctx_has_value (VB, FORMAT_RD) && 
        str_get_int (STRa(ad_str), &ad) &&
        ad == CTX(FORMAT_DP)->last_value.i - CTX(FORMAT_RD)->last_value.i) 
    
        vcf_seg_FORMAT_minus (vb, ctx, 0, ad_str_len, ad, CTX(FORMAT_DP), CTX(FORMAT_RD), STRa(ad_varscan_snip));

    // case: we have only one sample, and INFO/ADP - we expect FORMAT/AD and INFO/ADP to be related
    else if (ctx_has_value_in_line_(vb, CTX(INFO_ADP)) && vcf_num_samples==1)
        seg_delta_vs_other_localS (VB, ctx, CTX(INFO_ADP), STRa(ad_str), -1);

    else
        seg_by_ctx (VB, STRa(ad_str), ctx, ad_str_len);
}

// <ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
static void vcf_seg_AD_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values)
{       
    // calculate sum of AD items
    int64_t sum_ad = 0; 
    for (unsigned i=0; i < n_items; i++) 
        sum_ad += values[i];

    ctx_set_last_value (VB, ctx, sum_ad); // AD value is sum of its items

    // if we have ADALL in this file, we delta vs ADALL if we have it in this sample, or seg normally if not
    if (segconf.has[FORMAT_ADALL]) { // note: we can't delta vs ADALL unless segconf says so, bc it can ruin other fields that rely on peeking AD, eg AB
        bool has_adall_this_sample = ctx_encountered (VB, FORMAT_ADALL); 
    
        for (unsigned i=0; i < n_items; i++) 
            // case: we had ADALL preceeding in this sample, seg as delta vs. ADALL 
            if (has_adall_this_sample)
                seg_delta_vs_other_localN (VB, item_ctxs[i], ECTX (con_ADALL.items[i].dict_id), values[i], -1, item_lens[i]);
            else 
                seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
    }

    else 
        for (unsigned i=0; i < n_items; i++) {
            if (i==0 && vcf_num_samples > 1)
                // we store the sum in the first item - it is expected to be roughly similar within each sample - hence transposing
                vcf_seg_FORMAT_transposed (vb, item_ctxs[0], 0, 0, sum_ad, (char []){ SNIP_SPECIAL, VCF_SPECIAL_FORMAT_AD0 }, 2, item_lens[0]);
            
            else if (i==0 || i==1) {
                if (!vb->mux_AD[i].num_channels) 
                    seg_mux_init (vb, item_ctxs[i]->did_i, VCF_SPECIAL_MUX_BY_DOSAGE, false, AD[i]);

                vcf_seg_FORMAT_mux_by_dosage_int (vb, item_ctxs[i], values[i], &vb->mux_AD[i], item_lens[i]);
            }
            else
                seg_integer_or_not (VB, item_ctxs[i], STRi(item, i), item_lens[i]);
        }

    memcpy (vb->ad_values, values, n_items * sizeof (values[0]));
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_FORMAT_AD0)
{
    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, (char[]){ SNIP_LOOKUP, 0 }, 1, false, __FUNCLINE); // note: nul-termianted as expected of a dictionary snip

    new_value->i = ctx->last_value.i; // calculated by ^, equals Σ(ADᵢ)

    // calculate: AD0 = Σ(ADᵢ) - AD₁ - AD₂...    
    for (int i=1; i < con_nitems (*current_con.con); i++)
        new_value->i -= reconstruct_peek_by_dict_id (vb, current_con.con->items[i].dict_id, 0, 0).i;

    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
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
            seg_by_ctx (VB, STRi(my_snip,i), item_ctxs[i], item_lens[i]); 

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

    for (int i=0; i < n_items; i++) {
        
        // seg odd-numbered element as AD - (even element), if the sum is correct
        if (use_formula && i%2 && vb->ad_values[i/2] == values[i-1] + values[i]) 
            seg_by_ctx (VB, sb_snips[i/2], sb_snip_lens[i/2], item_ctxs[i], item_lens[i]); 

        else {
            item_ctxs[i]->flags.store = STORE_INT; // consumed by the odd items ^
            seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
        }
    }

    if (segconf.AS_SB_TABLE_by_SB && n_items == 4) // our method only works for bi-allelic (i.e. SB.n_items==4)
        for (int i=0; i < 4; i++) 
            ctx->sum_sb[i] += values[i]; // consumed by vcf_seg_INFO_AS_SB_TABLE
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
            seg_delta_vs_other_localN (VB, item_ctxs[i], ECTX (con_F2R1.items[i/2].dict_id), values[i], -1, item_lens[i]);
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
        ctx_encountered_in_line (VB, INFO_AF) && 
        str_issame_(STRa(cell), STRlst(INFO_AF)))
        return seg_by_ctx (VB, STRa(snip_copy_af), ctx, cell_len);
    else
        return vcf_seg_FORMAT_A_R (vb, ctx, con_AF, STRa(cell), STORE_FLOAT, NULL);
}


//-----------
// FORMAT/RGQ
// ----------

// <ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
// Appears in GVCF in lines which are no variants (i.e. no ALT)
static inline void vcf_seg_FORMAT_RGQ (VBlockVCFP vb, ContextP ctx, STRp(rgq), ContextP gt_ctx, STRp(gt))
{
    ConstMultiplexerP mux = (ConstMultiplexerP)&vb->mux_RGQ;
        
    // prediction: we have GT, and if GT[0]=. then RGQ=0. Fallback seg in case prediction fails
    if (gt_ctx->did_i != FORMAT_GT ||               // prediction failed: first subfield isn't GT
        (gt[0] == '.' && !str_is_1char (rgq, '0'))) // prediction failed: GT[0]=. and yet RGQ!="0"
        goto fallback;

    // case: GT[0] is not '.' - seg the value of RGQ multiplexed by DP
    if (gt[0] != '.') {
        if (!segconf.has[FORMAT_DP]          ||    // segconf didn't detect FORMAT/DP so we didn't initialize the mux
            !ctx_encountered (VB, FORMAT_DP) ||    // no DP in the FORMAT of this line
            segconf.running) goto fallback;        // multiplexer not initalized yet 

        int64_t DP;
        if (!str_get_int (STRlst(FORMAT_DP), &DP)) // in some files, DP may be '.'
            DP=0;

        int channel_i = MAX_(0, MIN_(DP, mux->num_channels-1));
        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, FORMAT_RGQ, MUX, channel_i);

        seg_integer_or_not (VB, channel_ctx, STRa(rgq), rgq_len);
    }

    seg_by_ctx (VB, MUX_SNIP(mux), mux->snip_len, ctx, (gt[0] == '.' ? rgq_len : 0));
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
    
        return HAS_NEW_VALUE; 
    }
    
    // gt[0] != '.' - demulitplex by FORMAT_DP
    else {
        unsigned num_channels = recon_multi_dict_id_get_num_dicts (ctx, STRa(snip));

        rom DP_str;
        int64_t DP = reconstruct_peek (vb, CTX(FORMAT_DP), &DP_str, 0).i;
        int channel_i = (*DP_str=='.') ? 0 : MAX_(0, MIN_(DP, num_channels-1));

        return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
    }
}

//----------
// FORMAT/DP
// ---------

static inline void vcf_seg_FORMAT_DP (VBlockVCFP vb)
{
    decl_ctx(FORMAT_DP);
    STRlast (dp, FORMAT_DP);

    int64_t value;
    bool is_null = IS_PERIOD(dp);
    bool has_value = !is_null && str_get_int (STRa(dp), &value);

    if (segconf.running) {
        if (is_null) 
            segconf.use_null_DP_method = true;
        
        if (!segconf.FMT_DP_method) {
            if (ctx_has_value (VB, FORMAT_SDP)) // SDP has priority of AD (in files with SDP, AD means something else)
                segconf.FMT_DP_method = BY_SDP;
            
            else if (ctx_has_value (VB, FORMAT_AD)) 
                segconf.FMT_DP_method = BY_AD;

            else if (vcf_num_samples == 1 && ctx_has_value (VB, INFO_DP))
                segconf.FMT_DP_method = BY_INFO_DP;
        }

        goto fallback;
    }

    else if ((segconf.FMT_DP_method == BY_AD || segconf.FMT_DP_method == BY_SDP) && !segconf.running) {
        Did other_did_i = (segconf.FMT_DP_method == BY_AD ? FORMAT_AD : FORMAT_SDP);

        bool other_is_in_FORMAT = ctx_encountered (VB, other_did_i); // note: DP is always segged after AD and SDP, regardless of their order in FORMAT
        int channel_i = other_is_in_FORMAT;

        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, ctx->did_i, (MultiplexerP)&vb->mux_FORMAT_DP, channel_i);

        if (channel_i) 
            seg_delta_vs_other_localS (VB, channel_ctx, CTX(other_did_i), STRa(dp), -1);
        
        else
            vcf_seg_FORMAT_transposed (vb, channel_ctx, STRa(dp), 0, (char []){ SNIP_LOOKUP }, 1, dp_len); // this handles DP that is an integer or '.'
        
        seg_by_ctx (VB, STRa(vb->mux_FORMAT_DP.snip), ctx, 0); 

        ctx_set_last_value (VB, ctx, channel_ctx->last_value); // propagate up       
    }
    
    // case: seg against INFO/DP or not
    else if (segconf.FMT_DP_method == BY_INFO_DP) {
        bool info_dp_is_int = ctx_has_value_in_line_(VB, CTX(INFO_DP));

        // special means: if FORMAT/DP>=1, INFO/DP is an integer, if FORMAT/DP==0, INFO/DP has no integer value
        if (has_value && ((value > 0) == info_dp_is_int)) {
            SNIPi2 (SNIP_SPECIAL, VCF_SPECIAL_DP_by_DP_single, info_dp_is_int ? (value - CTX(INFO_DP)->last_value.i) : 0);
            seg_by_ctx (VB, STRa(snip), ctx, dp_len);
        }

        else fallback: // value is not an integer or '.', or INFO and FORMAT don't agree on existance of an integer DP
            seg_by_ctx (VB, STRa(dp), ctx, dp_len);
    }

    // multi-sample default: store in transposed matrix (or just LOOKUP from local if not transposable)
    else if (vcf_num_samples > 1) 
        vcf_seg_FORMAT_transposed (vb, ctx, STRa(dp), 0, (char []){ SNIP_LOOKUP }, 1, dp_len); // this handles DP that is an integer or '.'

    else
        seg_integer_or_not (VB, ctx, STRa(dp), dp_len);
        
    if (!IS_PERIOD(dp) && segconf.INFO_DP_method == BY_FORMAT_DP) 
        CTX(INFO_DP)->dp.sum_format_dp += ctx->last_value.i; // we may delta INFO/DP against this sum

    // add up DP's in certain conditions, for consumption by INFO/QD predictor
    if (has_value)
        vcf_seg_sum_DP_for_QD (vb, value);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_FORMAT_DP)
{
    bool other_is_in_FORMAT = curr_container_has (vb, (segconf.FMT_DP_method == BY_AD) ? _FORMAT_AD : _FORMAT_SDP);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), other_is_in_FORMAT, new_value, reconstruct);
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
        if (IS_PERIOD (pid)) 
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

//----------
// FORMAT/AB
// ---------

// <ID=AB,Number=1,Type=Float,Description="Allele balance for each het genotype",RendAlg="NONE">
// Expecting: '.' if channel is 0 or 2,  if 1. 
// If expectation is met, SPECIAL is segged in AB. 
static inline void vcf_seg_FORMAT_AB (VBlockVCFP vb, ContextP ctx, STRp(ab))
{
    int channel_i = vcf_seg_get_mux_channel_i (vb);

    if (channel_i == -1                     || // GT didn't produce a mux channel
        CTX(FORMAT_GT)->gt.prev_ploidy != 2 || // This method was tested only for ploidy == 2
        vb->n_alts != 1                     || // This method was tested only for n_alts == 1
        (channel_i != 1 && !IS_PERIOD(ab))) {  // prediction: '.' for channel 0,2. fail if prediction is wrong.
    
        seg_by_ctx (VB, STRa(ab), ctx, ab_len);
        return;
    }

    // prediction: AD0/(AD0+AD1) for channel 1. We will verify the prediction in 
    // vcf_seg_FORMAT_AB_verify_channel1 and rollback if prediction is wrong
    if (channel_i == 1) {
        seg_set_last_txt (VB, ctx, STRa(ab));
        ctx_set_last_value (VB, ctx, (ValueType){.i = 1}); // need verification

        seg_create_rollback_point (VB, NULL, 1, FORMAT_AB); 
    }

    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_AB }, 2, ctx, ab_len);

    // note: up to 15.0.35, we segged channel_i=3 (exceptions channel) into FORMAT_AB3. Since 15.0.36 we no longer have an exception channel.
}

static inline void calculate_AB (double ad0, double ad1, qSTRp(snip))
{
    double ab = 0.0000001 + ad0 / (ad0 + ad1); // +epsilon to bring number over the 0.01 mark if it is almost almost there 
    *snip_len = snprintf (snip, *snip_len, "%.*g", ab < 0.1 ? 1 : 2, ab); 
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

    STRlic(recon_ab_str,32);
    calculate_AB (ad0, ad1, qSTRa(recon_ab_str));

    if (!str_issame (ab_str, recon_ab_str)) goto rollback;
    
    ctx_unset_rollback (ab_ctx); // so we can ctx_set_rollback again in the next sample (note: we can't advance rback_id because it is set in vcf_seg_txt_line)
    return; // verified

rollback:
    seg_rollback (VB);
    seg_by_ctx (VB, STRa(ab_str), ab_ctx, ab_str_len); 
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_AB)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    int channel_i = vcf_piz_get_mux_channel_i (VB);

    if ((channel_i == 0 || channel_i == 2) && reconstruct) 
        RECONSTRUCT1 ('.');

    else if (channel_i == 1 && reconstruct) {
        STR(ad);
        reconstruct_peek (VB, CTX(FORMAT_AD), pSTRa(ad));

        // note: using MAX_ALLELES and not vb->n_alts, bc n_alts is not defined for files older than 14.0.12
        str_split_ints (ad, ad_len, MAX_ALLELES + 1, ',', ad, false); // note: up to 15.0.36 we could arrive here even with n_alts > 1 
        ASSPIZ (n_ads, "Failed to split AD=\"%.*s\" to %u integers", STRf(ad), VB_VCF->n_alts + 1);

        STRlic(recon_ab_str, 32);
        calculate_AB (ads[0], ads[1], qSTRa(recon_ab_str));
        RECONSTRUCT_str (recon_ab_str);
    }

    else if (channel_i == 3) { // back comp: only happens in files <= 15.0.35
        ContextP ab3_ctx = MCTX (2, snip, snip_len); // in old files, snip contained the dict_ids of AD0,AD1,AB3
        reconstruct_from_ctx (VB, ab3_ctx->did_i, 0, reconstruct);
    }

    return NO_NEW_VALUE;
}

//----------
// FORMAT/GP
// ---------

static inline void vcf_seg_FORMAT_GP (VBlockVCFP vb, ContextP ctx, STRp(gp))
{
    // according to the VCF spec, until VCF 4.2 GP contains phred values
    // and since 4.3, it contains probabilities. But tools don't adhere to this, and it might be either.
    // we test all GPs in the segconf data to identify their format and verify its consistency   
    if (segconf.running && segconf.FMT_GP_content == GP_unknown) {
        str_split_floats (gp, gp_len, 0, ',', gp, false, '.');
        if (!n_gps)
            segconf.FMT_GP_content = GP_other; // unrecognized format

        // case: probabilities format if numbers add up approximately 1
        else if (({ double sum=0 ; 
                    for (int i=0; i < n_gps; i++) sum += (isnan(gps[i]) ? 0 : gps[i]); 
                    (sum > 0.98 && sum < 1.02); }))
            segconf.FMT_GP_content = GP_probabilities; 

        // case: evidence of phred if all numbers are non-negative, but are not probabilities
        else if (({ bool all_non_negative=true; 
                    for (int i=0; i < n_gps; i++) if (!isnan(gps[i]) && gps[i] < 0) { all_non_negative=false; break; }; 
                    all_non_negative; }))
            segconf.FMT_GP_content = GP_phred; 

        // case: negative numbers - unknown format
        else
            segconf.FMT_GP_content = GP_other;
    }

    vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRa(gp), &vb->mux_GP);
    seg_set_last_txt (VB, ctx, STRa(gp)); // used by GQ
}

//----------
// FORMAT/PL
// ---------

// <ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">       
static inline void vcf_seg_FORMAT_PL (VBlockVCFP vb, ContextP ctx, STRp(PL))
{
    if (segconf.running && !segconf.has_DP_before_PL) 
        segconf.has_DP_before_PL = ctx_encountered (VB, FORMAT_DP);

    seg_set_last_txt (VB, ctx, STRa(PL)); // consumed by GQ and SPL 
       
    // seg into either PLy or PLn, or into both if we're testing (vcf_FORMAT_PL_decide will drop one of them) 
    if (vb->PL_mux_by_DP == yes || vb->PL_mux_by_DP == unknown) 
        vcf_seg_FORMAT_mux_by_dosagexDP (vb, CTX(FORMAT_PLy), STRa(PL), &vb->mux_PLy);
    
    if (vb->PL_mux_by_DP == no  || vb->PL_mux_by_DP == unknown) 
        vcf_seg_FORMAT_mux_by_dosage (vb, CTX(FORMAT_PLn), STRa(PL), &vb->mux_PLn);

    if (vb->PL_mux_by_DP == unknown) 
        vb->recon_size += PL_len; // since we're segging twice, we need to pretent the recon_size is growing even thow it isn't (reversed in zfile_remove_ctx_group_from_z_data)

    // we seg "redirection to PLn" in all VBs regardless of PL_mux_by_DP, so that ZCTX(FORMAT_PL).dict has
    // only a single word. This word might be updated to PLy in vcf_FORMAT_PL_after_vbs.
    // note: this isn't always "all-the-same" - we can have WORD_INDEX_MISSING in b250 in case of missing PLs in samples
    seg_by_did (VB, STRa(PL_to_PLn_redirect_snip), FORMAT_PL, 0);
}

// in unknown (=test) mode, this function is called to pick one of the two compressed options, and drop the others
void vcf_FORMAT_PL_decide (VBlockVCFP vb)
{
    Did keep_did_i, remove_did_i;

    mutex_lock (segconf.PL_mux_by_DP_mutex);

    // only first testing VB to lock this mutex gets here, to make the decision
    switch (segconf.PL_mux_by_DP) {
        case unknown : { 
            // get combined compress size of all contexts associated with GL (mux by dosage)
            uint64_t PLy_size = ctx_get_ctx_group_z_len (VB, FORMAT_PLy);
            uint64_t PLn_size = ctx_get_ctx_group_z_len (VB, FORMAT_PLn);

            // note: this is not an accurate comparison bc it doesn't include dictionaries
            keep_did_i   = (PLy_size >= PLn_size) ? FORMAT_PLn : FORMAT_PLy;
            remove_did_i = (PLy_size >= PLn_size) ? FORMAT_PLy : FORMAT_PLn;

            // weigh in this VBs experience. New VBs will only test if the votes are still the same
            segconf.PL_mux_by_DP = (keep_did_i == FORMAT_PLy) ? yes : no;
            
            if (flag.debug_generate) iprintf ("PL_mux_by_DP decision: %s\n", keep_did_i == FORMAT_PLy ? "YES" : "NO");
            
            break;
        }
        case yes : keep_did_i = FORMAT_PLy; remove_did_i = FORMAT_PLn; break;

        default:
        case no  : keep_did_i = FORMAT_PLn; remove_did_i = FORMAT_PLy;
    }

    mutex_unlock (segconf.PL_mux_by_DP_mutex);

    zfile_remove_ctx_group_from_z_data (VB, remove_did_i); // removes data from VB and update z_file counters
}

// called after all VBs are compressed - before Global sections are compressed
void vcf_FORMAT_PL_after_vbs (void)
{
    if (!ZCTX(FORMAT_PL)->nodes.len) return; // no FORMAT/PL in this file

    if (ZCTX(FORMAT_PL)->nodes.len > 1) {
        dict_io_print (info_stream, STRb(ZCTX(FORMAT_PL)->dict), true, false, true, false);
        ABORT ("Expecting FORMAT_PL to have exactly one word in its dict, but it has %"PRIu64, ZCTX(FORMAT_PL)->nodes.len);
    }

    if (segconf.PL_mux_by_DP == yes) { // Tested, and selected YES
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
static inline WordIndex vcf_seg_FORMAT_DS (VBlockVCFP vb, ContextP ctx, rom ds, unsigned ds_len)
{
    int64_t dosage = ctx_has_value (VB, FORMAT_GT) ? CTX(FORMAT_GT)->last_value.i : -1; // dosage stored here by vcf_seg_FORMAT_GT
    double ds_val;
    unsigned format_len;
    char snip[FLOAT_FORMAT_LEN + 20] = { SNIP_SPECIAL, VCF_SPECIAL_DS }; 

    if (dosage < 0 || !str_get_float (ds, ds_len, &ds_val, &snip[2], &format_len)) 
        return seg_by_ctx (VB, ds, ds_len, ctx, ds_len);

    unsigned snip_len = 2 + format_len;
    snip[snip_len++] = ' ';
    snip_len += str_int ((int64_t)((ds_val - dosage) * 1000000), &snip[snip_len]);

    return seg_by_ctx (VB, snip, snip_len, ctx, ds_len);
}
*/

// used for decompressing files compressed with version up to 12.0.42
SPECIAL_RECONSTRUCTOR (vcf_piz_special_DS_old)
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
static inline unsigned vcf_seg_one_sample (VBlockVCFP vb, ZipDataLineVCF *dl, ContextP *ctxs, ContainerP format, STRp(sample))
{
    #define COND0(condition, seg) if (condition) { seg; break; } else  
    #define COND(condition,  seg) if (condition) { seg; break; } else goto fallback; 

    str_split (sample, sample_len, con_nitems (*format), ':', sf, false);

    ASSVCF (n_sfs, "Sample %u has too many subfields - FORMAT field \"%s\" specifies only %u: \"%.*s\"", 
            vb->sample_i+1, error_format_field (con_nitems (*format), ctxs), con_nitems (*format), STRf(sample));

    for (unsigned i=0; i < n_sfs; i++) { 

        DictId dict_id = format->items[i].dict_id;
        ContextP ctx = ctxs[i];

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

        // note: GP and PL - for non-optimized, I tested segging as A_R and seg_array - they are worse or not better than the default, because the values are correlated.
        case _FORMAT_GP: vcf_seg_FORMAT_GP (vb, ctx, STRi (sf, i)); break;

        // <ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">       
        case _FORMAT_PL: vcf_seg_FORMAT_PL (vb, ctx, STRi (sf, i)); break;

        // <ID=PP,Number=G,Type=Integer,Description="Phred-scaled genotype posterior probabilities rounded to the closest integer">
        case _FORMAT_PP:
            if (segconf.vcf_is_pindel)
                goto fallback; // PP means something entirely different in Pindel
            else
                vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_PP);
            break;
        
        // <ID=GL,Number=.,Type=Float,Description="Genotype Likelihoods">
        // note: I tried muxing against DS instead of dosage (41 channels) - worse results than dosage in 1KGP-37 even with --best 
        case _FORMAT_GL   : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_GL); break;

        // <ID=PRI,Number=G,Type=Float,Description="Phred-scaled prior probabilities for genotypes">
        case _FORMAT_PRI  : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_PRI); break;
        
        case _FORMAT_CN   : seg_integer_or_not (VB, ctx, STRi(sf, i), sf_lens[i]); break;

        // <ID=DS,Number=1,Type=Float,Description="Genotype dosage from MaCH/Thunder"> (1000 Genome Project phase1 data)
        // See: https://genome.sph.umich.edu/wiki/Thunder
        case _FORMAT_DS   : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_DS)   ; break;
        
        // VarScan: <ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
        case _FORMAT_RD   : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_RD)   ; break;
        
        // VarScan: <ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
        case _FORMAT_PVAL : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_PVAL) ; break;   
        
        // VarScan: <ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
        case _FORMAT_FREQ : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_FREQ) ; break;
        
        // <ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">
        case _FORMAT_VAF  : vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRi (sf, i), &vb->mux_VAF) ; break;
        
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
        case _FORMAT_GQ   : if (segconf.has[FORMAT_GQ]) {
                                if (segconf.vcf_is_isaac) seg_set_last_txt_store_value (VB, ctx, STRi(sf, i), STORE_INT);
                                else seg_set_last_txt (VB, ctx, STRi(sf, i)); 
                                break; // postpone to later
                            }
                            else goto fallback;
            
        case _FORMAT_RGQ  : vcf_seg_FORMAT_RGQ (vb, ctx, STRi(sf, i), ctxs[0], STRi(sf,0)); break;

        case _FORMAT_DP   : seg_set_last_txt (VB, ctx, STRi(sf, i)); break; // defer to later
            
        // <ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
        // case: MIN_DP - it is slightly smaller and usually equal to DP - we store MIN_DP as the delta DP-MIN_DP
        // note: the delta is vs. the DP field that preceeds MIN_DP - we take the DP as 0 there is no DP that preceeds
        case _FORMAT_MIN_DP :
            COND (ctx_has_value (VB, FORMAT_DP), seg_delta_vs_other_localS (VB, ctx, CTX(FORMAT_DP), STRi(sf, i), -1));

        case _FORMAT_SDP   :
            if (ctx_has_value_in_line_(VB, CTX(INFO_ADP)))
                seg_delta_vs_other_localS (VB, ctx, CTX(INFO_ADP), STRi(sf, i), -1);
            else goto fallback;
            break;
            
        case _FORMAT_AF    : vcf_seg_FORMAT_AF (vb, ctx, STRi(sf, i)); break;

        // standard: <ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">  
        // GIAB:     <ID=AD,Number=R,Type=Integer,Description="Net allele depths across all unfiltered datasets with called genotype">
        case _FORMAT_AD    : if (segconf.vcf_is_varscan)
                                vcf_seg_FORMAT_AD_varscan (vb, ctx, STRi(sf, i));
                             else {
                                vcf_seg_FORMAT_A_R (vb, ctx, con_AD, STRi(sf, i), STORE_INT, vcf_seg_AD_items); 
                                if (segconf.has[INFO_BaseCounts]) seg_set_last_txt (VB, ctx, STRi(sf, i)); // consumed by vcf_seg_INFO_BaseCounts
                             }
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
        case _FORMAT_DPF   : COND (segconf.vcf_is_isaac, seg_integer_or_not (VB, ctx, STRi(sf, i), sf_lens[i]));

        // DRAGEN fields
        case _FORMAT_PE    : COND (segconf.vcf_is_dragen, seg_array (VB, ctx, ctx->did_i, STRi(sf, i), ',', 0, false, STORE_INT, DICT_ID_NONE, sf_lens[i]));
        case _FORMAT_BC    : COND (segconf.vcf_is_dragen, seg_integer_or_not (VB, ctx, STRi(sf, i), sf_lens[i]));

        // manta fields
        case _FORMAT_SR    :
        case _FORMAT_PR    : COND (segconf.vcf_is_manta, seg_array (VB, ctx, ctx->did_i, STRi(sf, i), ',', 0, false, STORE_INT, DICT_ID_NONE, sf_lens[i]));
        
        // Platypus fields
        case _FORMAT_GOF   : COND (segconf.vcf_is_platypus, vcf_seg_platypus_FORMAT_GOF (vb, ctx, STRi(sf, i)));

        // GLIMPSE_phase fields
        case _FORMAT_HS    : COND (segconf.vcf_is_GLIMPSE_phase, seg_integer_or_not (VB, ctx, STRi(sf, i), sf_lens[i]));

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
    for (unsigned i=n_sfs; i < con_nitems (*format); i++) {
        uint64_t dnum = format->items[i].dict_id.num;

        // special handling for PS and PID
        if (dnum == _FORMAT_PS || dnum == _FORMAT_PID) 
            vcf_seg_FORMAT_PS_PID_missing_value (vb, ctxs[i], &sample[sample_len]);
        
        else 
            seg_by_ctx (VB, NULL, 0, ctxs[i], 0); // generates WORD_INDEX_MISSING
    }
    
    // verify AB if its channel 1
    if (ctx_has_value (VB, FORMAT_AB))
        vcf_seg_FORMAT_AB_verify_channel1 (vb);

    // seg DP (must be after AD, SDP)
    if (ctx_encountered (VB, FORMAT_DP))
        vcf_seg_FORMAT_DP (vb); 

    // finally seg GQ if we have it (must be after after GP, PL, DP)
    if (segconf.has[FORMAT_GQ] && ctx_encountered (VB, FORMAT_GQ))
        vcf_seg_FORMAT_GQ (vb);

    return n_sfs - 1; // number of colons
}

//------------
// All samples
//------------

rom vcf_seg_samples (VBlockVCFP vb, ZipDataLineVCF *dl, int32_t len, char *next_field, bool *has_13)
{
    // Container for samples - we have:
    // - repeats as the number of samples in the line (<= vcf_num_samples)
    // - n_items as the number of FORMAT subfields (inc. GT)
    decl_ctx (VCF_SAMPLES);
    Container format = *B(Container, ctx->format_mapper_buf, dl->format_node_i); // make a copy of the template
    ContextP *ctxs = B(ContextP, ctx->format_contexts, dl->format_node_i * MAX_FIELDS);

    // structural variants: entire samples data is expected to be identical between BND mates
    if (segconf.vcf_is_svaba || segconf.vcf_is_manta) {
        SAFE_NUL (next_field + len);
        uint32_t samples_len = strcspn (next_field, "\n\r");
        SAFE_RESTORE;

        ContextP channel_ctx = vcf_seg_sv_SAMPLES (vb, next_field, len, ctxs, con_nitems(format));
        if (!channel_ctx)  
            return next_field + samples_len + (next_field[samples_len]=='\r') + 1/*\n*/; // segged as copy from mate
        else
            ctx = channel_ctx;
    }

    // set ctx->sf_i - for the line's FORMAT fields
    for (int sf_i=0; sf_i < con_nitems(format); sf_i++) {
        ctxs[sf_i]->sf_i = sf_i;

        if (segconf.running) segconf.has[ctxs[sf_i]->did_i]++;
    }

    // initialize LOOKBACK if we have PS or PID
    if (!CTX(VCF_LOOKBACK)->is_initialized && (CTX(FORMAT_PID)->sf_i >= 0 || CTX(FORMAT_PS)->sf_i >= 0))
        vcf_samples_seg_initialize_LOOKBACK (vb);
    
    rom field_start;
    unsigned field_len=0, num_colons=0;

    // 0 or more samples (note: we don't use str_split, because samples could be very numerous)
    for (char separator=0 ; separator != '\n'; format.repeats++) {

        field_start = next_field;
        next_field = (char *)seg_get_next_item (VB, field_start, &len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, has_13, "sample-subfield");

        ASSVCF (field_len, "Error: invalid VCF file - expecting sample data for sample # %u, but found a tab character", 
                format.repeats+1);

        vb->sample_i = format.repeats;
        num_colons += vcf_seg_one_sample (vb, dl, ctxs, &format, (char *)field_start, field_len);

        ASSVCF (format.repeats < vcf_num_samples || separator == '\n',
                "invalid VCF file - expecting a newline after the last sample (sample #%u)", vcf_num_samples);
    }
    vb->sample_i = 0;
    
    ASSVCF (format.repeats <= vcf_num_samples, "according the VCF header, there should be %u sample%s per line, but this line has %u samples - that's too many",
            STRfN(vcf_num_samples), format.repeats);

    // in some real-world files I encountered have too-short lines due to human errors. we pad them
    if (format.repeats < vcf_num_samples) 
        WARN_ONCE ("FYI: the number of samples in variant CHROM=%.*s POS=%"PRId64" is %u, different than the VCF column header line which has %u samples",
                   vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS), format.repeats, vcf_num_samples);

    vcf_seg_FORMAT_GT_finalize_line (vb, format.repeats);

    container_seg (vb, ctx, &format, 0, 0, format.repeats + num_colons); // account for : and \t \r \n separators

    ctx_set_last_value (VB, CTX(VCF_SAMPLES), (ValueType){ .i = format.repeats });
 
    return next_field;
}

