// ------------------------------------------------------------------
//   vcf_info.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "vcf_private.h"
#include "piz.h"
#include "optimize.h"
#include "file.h"
#include "dict_id.h"
#include "codec.h"
#include "reconstruct.h"
#include "stats.h"

static inline bool vcf_is_use_DP_by_DP (void); // forward
sSTRl(RAW_MQandDP_snip, 64 + 2 * 16);     

static SmallContainer RAW_MQandDP_con = {
    .nitems_lo = 2, 
    .repeats   = 1, 
    .items     = { { .dict_id={ _INFO_RAW_MQandDP_MQ }, .separator = {','} },
                   { .dict_id={ _INFO_RAW_MQandDP_DP }                     } }
};

// called after reading VCF header, before segconf
void vcf_info_zip_initialize (void) 
{
    container_prepare_snip ((ContainerP)&RAW_MQandDP_con, 0, 0, qSTRa(RAW_MQandDP_snip));

    vcf_dbsnp_zip_initialize(); // called even if not in VCF header, because can be discovered in segconf too
    if (segconf.vcf_is_vagrent)    vcf_vagrent_zip_initialize();
    if (segconf.vcf_is_mastermind) vcf_mastermind_zip_initialize();
    if (segconf.vcf_is_vep)        vcf_vep_zip_initialize();
    if (segconf.vcf_illum_gtyping) vcf_illum_gtyping_zip_initialize();
}

void vcf_info_seg_initialize (VBlockVCFP vb) 
{
    ctx_set_store (VB, STORE_INT, INFO_AN, INFO_AC, INFO_ADP, INFO_DP, INFO_MLEAC, 
                   INFO_DP4_RF, INFO_DP4_RR, INFO_DP4_AF, INFO_DP4_AR, 
                   INFO_AC_Hom, INFO_AC_Het, INFO_AC_Hemi,
                   DID_EOL);

    CTX(INFO_AF)->flags.store = STORE_FLOAT;
    // xxx (is this really needed for --indels-only?) CTX(INFO_SVTYPE)-> flags.store = STORE_INDEX; // since v13 - consumed by vcf_refalt_piz_is_variant_indel

    if (!vcf_is_use_DP_by_DP())
        CTX(INFO_DP)->ltype = LT_DYN_INT; 

    ctx_consolidate_stats (VB, INFO_RAW_MQandDP, INFO_RAW_MQandDP_MQ, INFO_RAW_MQandDP_DP, DID_EOL);
    
    if (segconf.has[INFO_CLNHGVS]) vcf_seg_hgvs_consolidate_stats (vb, INFO_CLNHGVS);
    if (segconf.has[INFO_HGVSG])   vcf_seg_hgvs_consolidate_stats (vb, INFO_HGVSG);
    if (segconf.has[INFO_ANN])     vcf_seg_hgvs_consolidate_stats (vb, INFO_ANN); // subfield HGVS_c

    CTX(INFO_SVLEN)->ltype = CTX(INFO_DP4_RF)->ltype = CTX(INFO_DP4_AF)->ltype = LT_DYN_INT;
}

//--------
// INFO/DP
// -------

static inline bool vcf_is_use_DP_by_DP (void)
{
    // note: this method causes genozip --drop-genotypes, --GT-only, --samples to show INFO/DP=-1
    // user can specify --secure-DP to avoid this
    return vcf_num_samples > 1 && !flag.secure_DP && segconf.has[FORMAT_DP];
}

// return true if caller still needs to seg 
static void vcf_seg_INFO_DP (VBlockVCFP vb, ContextP ctx, STRp(value_str))
{
    // used in: vcf_seg_one_sample (for 1-sample files), vcf_seg_INFO_DP_by_FORMAT_DP (multi sample files)
    int64_t value;
    bool has_value = str_get_int (STRa(value_str), &value);

    // also tried delta vs DP4, but it made it worse
    ContextP ctx_basecounts;
    if (ctx_has_value_in_line (vb, _INFO_BaseCounts, &ctx_basecounts)) 
        seg_delta_vs_other (VB, ctx, ctx_basecounts, STRa(value_str));

    // note: if we're doing multi-sample DP_by_DP, we will seg in vcf_seg_INFO_DP_by_FORMAT_DP
    else if (!has_value || !vcf_is_use_DP_by_DP()) 
        seg_integer_or_not (VB, ctx, STRa(value_str), value_str_len);

    // we're going to seg (in vcf_finalize_seg_info) INFO/DP against sum of FORMAT/DP (if it has a valid value)
    else
        ctx->dp.by_format_dp = has_value;

    if (has_value) 
        ctx_set_last_value (VB, ctx, value);
}

// used for multi-sample VCFs, IF FORMAT/DP is segged as a simple integer
static void vcf_seg_INFO_DP_by_FORMAT_DP (VBlockVCFP vb)
{
    decl_ctx (INFO_DP);

    int value_len = str_int_len (ctx->last_value.i);

    SNIP(32) = { SNIP_SPECIAL, VCF_SPECIAL_DP_by_DP };
    snip_len = 2 + str_int (value_len, &snip[2]);
    snip[snip_len++] = '\t';

    // note: INFO/DP >= sum(FORMAT/DP) as the per-sample value is filtered, see: https://gatk.broadinstitute.org/hc/en-us/articles/360036891012-DepthPerSampleHC
    snip_len += str_int (ctx->last_value.i - ctx->dp.sum_format_dp, &snip[snip_len]);

    seg_by_ctx (VB, STRa(snip), ctx, value_len); 
}

// initialize reconstructing INFO/DP by sum(FORMAT/DP) - save space in txt_data, and initialize delta
SPECIAL_RECONSTRUCTOR (vcf_piz_special_DP_by_DP)
{
    str_split (snip, snip_len, 2, '\t', item, 2);

    if (reconstruct) {
        if (!flag.drop_genotypes && !flag.gt_only && !flag.samples) {
            Ltxt += atoi (items[0]); // number of characters needed to reconstruct the INFO/DP integer
            ctx->dp.sum_format_dp = atoi (items[1]); // initialize with delta
            ctx->dp.by_format_dp = true;             // needs to be finalized
        }
        else
            RECONSTRUCT ("-1", 2); // bc we can't calculate INFO/DP in these cases bc we need FORMAT/DP of all samples
    }

    return NO_NEW_VALUE; 
}

// finalize reconstructing INFO/DP by sum(FORMAT/DP) - called after reconstructing all samples
void vcf_piz_finalize_DP_by_DP (VBlockVCFP vb)
{
    str_int_ex (CTX(INFO_DP)->dp.sum_format_dp, last_txt(VB, INFO_DP), false);
}

// used starting v13.0.5, replaced in v14 with a new vcf_piz_special_DP_by_DP
SPECIAL_RECONSTRUCTOR (vcf_piz_special_DP_by_DP_v13)
{
    str_split (snip, snip_len, 2, '\t', item, 2);

    int num_dps_this_line   = atoi (items[0]);
    int64_t value_minus_sum = atoi (items[1]);

    ContextP format_dp_ctx = CTX(FORMAT_DP);

    int64_t sum=0;

    ASSPIZ (format_dp_ctx->next_local + num_dps_this_line <= format_dp_ctx->local.len, "Not enough data in FORMAT/DP.local to reconstructed INFO/DP: next_local=%u local.len=%u but needed num_dps_this_line=%u",
            format_dp_ctx->next_local, format_dp_ctx->local.len32, num_dps_this_line);
            
    uint32_t invalid = lt_desc[format_dp_ctx->ltype].max_int; // represents '.'
    for (int i=0; i < num_dps_this_line; i++) {
        uint32_t format_dp = (format_dp_ctx->ltype == LT_UINT8)  ? (uint32_t)*B8 ( format_dp_ctx->local, format_dp_ctx->next_local + i)
                           : (format_dp_ctx->ltype == LT_UINT16) ? (uint32_t)*B16 (format_dp_ctx->local, format_dp_ctx->next_local + i)
                           : /* LT_UINT32 */                       (uint32_t)*B32 (format_dp_ctx->local, format_dp_ctx->next_local + i);

        if (format_dp != invalid) sum += format_dp; 
    }

    new_value->i = value_minus_sum + sum;

    RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

static bool vcf_seg_INFO_DP4_delta (VBlockP vb, ContextP ctx, STRp(value), uint32_t unused_rep)
{
    if (ctx_encountered_in_line (vb, (ctx-1)->did_i)) {
        seg_delta_vs_other (VB, ctx, ctx-1, STRa(value));
        return true; // segged successfully
    }
    else
        return false;
}

// <ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
// Expecting first two values to be roughly similar, as the two last bases roughly similar
static void vcf_seg_INFO_DP4 (VBlockVCFP vb, ContextP ctx, STRp(dp4))
{
    static const MediumContainer container_DP4 = {
        .repeats      = 1, 
        .nitems_lo    = 4, 
        .items        = { { .dict_id.num = _INFO_DP4_RF, .separator = "," }, 
                          { .dict_id.num = _INFO_DP4_RR, .separator = "," }, 
                          { .dict_id.num = _INFO_DP4_AF, .separator = "," }, 
                          { .dict_id.num = _INFO_DP4_AR                   } } };

    SegCallback callbacks[4] = { 0, vcf_seg_INFO_DP4_delta, 0, vcf_seg_INFO_DP4_delta }; 

    seg_struct (VB, ctx, container_DP4, STRa(dp4), callbacks, dp4_len, true);
}

// -------
// INFO/AA
// -------

// return 0 if the allele equals main REF, the alt number if it equals one of the ALTs, or -1 if none or -2 if '.'
static int vcf_INFO_ALLELE_get_allele (VBlockVCFP vb, STRp (value))
{
    // check for '.'
    if (IS_PERIOD (value)) return -2;

    // check if its equal main REF (which can by REF or oREF)
    if (str_issame (value, vb->main_ref)) return 0;

    // check if its equal one of the ALTs
    str_split (vb->main_alt, vb->main_alt_len, 0, ',', alt, false);

    for (int alt_i=0; alt_i < n_alts; alt_i++) 
        if (str_issame_(STRa(value), STRi(alt, alt_i)))
            return alt_i + 1;

    // case: not REF or any of the ALTs
    return -1;
}

// checks if value is identifcal to the REF or one of the ALT alleles, and if so segs a SPECIAL snip
// Used for INFO/AA, INFO/CSQ/Allele, INFO/ANN/Allele. Any field using this should have the VCF2VCF_ALLELE translator set in vcf_lo_luft_trans_id.
bool vcf_seg_INFO_allele (VBlockP vb_, ContextP ctx, STRp(value), uint32_t repeat)  
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    
    int allele = vcf_INFO_ALLELE_get_allele (vb, STRa(value));

    // case: this is one of the alleles in REF/ALT - our special alg will just copy from that allele
    if (allele >= 0) {
        char snip[] = { SNIP_SPECIAL, VCF_SPECIAL_ALLELE, '0' + vb->line_coords, '0' + allele /* ASCII 48...147 */ };
        seg_by_ctx (VB, snip, sizeof (snip), ctx, value_len);
    }

    // case: a unique allele and no xstrand - we just leave it as-is
    else
        seg_by_ctx (VB, STRa(value), ctx, value_len); 

    // validate that the primary value (as received from caller or lifted back) can be luft-translated 
    // note: for INFO/AA, but not for INFO/CSQ/Allele and INFO/ANN/Allele, this is done already in vcf_seg_info_one_subfield (no harm in redoing)
    if (vb->line_coords == DC_PRIMARY && needs_translation (ctx)) {
        if (allele != -1) 
            ctx->line_is_luft_trans = true; // assign translator to this item in the container, to be activated with --luft
        else 
            REJECT_SUBFIELD (LO_INFO, ctx, ".\tCannot cross-render INFO subfield %s: \"%.*s\"", ctx->tag_name, value_len, value);            
    }

    return true; // segged successfully
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_ALLELE)
{
    Coords seg_line_coord = snip[0] - '0';
    int allele = snip[1] - '0';
    LiftOverStatus ostatus = last_ostatus;

    if (LO_IS_OK_SWITCH (ostatus) && seg_line_coord != VB_VCF->vb_coords) {
        ASSPIZ (allele >= 0 && allele <= 1, "unexpected allele=%d with REF<>ALT switch", allele);
        allele = 1 - allele;
    }

    ContextP refalt_ctx = VB_VCF->vb_coords == DC_PRIMARY ? CTX (VCF_REFALT) : CTX (VCF_oREFALT); 

    STRlast (refalt, refalt_ctx->did_i);

    if (!refalt_len) goto done; // variant is single coordinate in the other coordinate

    char *tab = memchr (refalt, '\t', refalt_len);
    ASSPIZ (tab, "Invalid refalt: \"%.*s\"", MIN_(refalt_len, 100), refalt);

    // case: the allele is REF
    if (allele == 0)
        RECONSTRUCT (refalt, tab - refalt);
    
    // case: the allele is one of the alts
    else {
        str_split (tab+1, &refalt[refalt_len] - (tab+1), 0, ',', alt, false);
        RECONSTRUCT (alts[allele-1], alt_lens[allele-1]);
    }

done:
    return NO_NEW_VALUE;
}

// translator only validates - as vcf_piz_special_ALLELE copies verbatim (revcomp, if xstrand, is already done in REF/ALT)
TRANSLATOR_FUNC (vcf_piz_luft_ALLELE)
{
    VBlockVCFP vcf_vb = VB_VCF;

    // reject if LO_OK_REF_NEW_SNP and value is equal to REF
    if (validate_only && last_ostatus == LO_OK_REF_NEW_SNP && str_issame (recon, vcf_vb->main_ref)) 
        return false;

    // reject if the value is not equal to REF, any ALT or '.'
    if (validate_only && vcf_INFO_ALLELE_get_allele (vcf_vb, STRa(recon)) == -1) return false;

    return true;
}

// ---------------
// INFO/BaseCounts
// ---------------

// ##INFO=<ID=genozip BugP.vcf -ft ,Number=4,Type=Integer,Description="Counts of each base">
// Sorts BaseCounts vector with REF bases first followed by ALT bases, as they are expected to have the highest values
static bool vcf_seg_INFO_BaseCounts (VBlockVCFP vb, ContextP ctx_basecounts, STRp(value)) // returns true if caller still needs to seg 
{
    if (vb->main_ref_len != 1 || vb->main_alt_len != 1 || vb->line_coords == DC_LUFT) 
        return true; // not a bi-allelic SNP or line is a luft line without easy access to REFALT - caller should seg

    char *str = (char *)value;
    int64_t sum = 0;

    uint32_t counts[4], sorted_counts[4] = {}; // corresponds to A, C, G, T

    SAFE_NUL (&value[value_len]);
    for (unsigned i=0; i < 4; i++) {
        counts[i] = strtoul (str, &str, 10);
        str++; // skip comma separator
        sum += counts[i];
    }
    SAFE_RESTORE;

    if (str - value != value_len + 1 /* +1 due to final str++ */) return true; // invalid BaseCounts data - caller should seg

    unsigned ref_i = acgt_encode[(int)*vb->main_ref];
    unsigned alt_i = acgt_encode[(int)*vb->main_alt];

    bool used[4] = {};
    sorted_counts[0] = counts[ref_i]; // first - the count of the REF base
    sorted_counts[1] = counts[alt_i]; // second - the count of the ALT base
    used[ref_i] = used[alt_i] = true;

    // finally - the other two cases in their original order (usually these are 0)
    for (unsigned sc_i=2; sc_i <= 3; sc_i++)
        for (unsigned c_i=0; c_i <= 3; c_i++)
            if (!used[c_i]) { // found a non-zero count
                sorted_counts[sc_i] = counts[c_i];
                used[c_i] = true;
                break;
            }

    char snip[2 + value_len + 1]; // +1 for \0
    sprintf (snip, "%c%c%u,%u,%u,%u", SNIP_SPECIAL, VCF_SPECIAL_BaseCounts, 
             sorted_counts[0], sorted_counts[1], sorted_counts[2], sorted_counts[3]);

    seg_by_ctx (VB, snip, value_len+2, ctx_basecounts, value_len); 
    
    ctx_basecounts->flags.store = STORE_INT;
    ctx_set_last_value (VB, ctx_basecounts, sum);

    return false; // we already segged - caller needn't seg
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_BaseCounts)
{
    STR (refalt);
    reconstruct_peek (vb, CTX (VCF_REFALT), pSTRa(refalt));

    uint32_t counts[4], sorted_counts[4] = {}; // counts of A, C, G, T

    new_value->i = 0;
    char *str = (char *)snip;

    for (unsigned i=0; i < 4; i++) {
        sorted_counts[i] = strtoul (str, &str, 10);
        str++; // skip comma separator
        new_value->i += sorted_counts[i];
    }

    if (!reconstruct) goto done; // just return the new value

    ASSVCF (str - snip == snip_len + 1, "expecting (str-snip)=%d == (snip_len+1)=%u", (int)(str - snip), snip_len+1);

    unsigned ref_i = acgt_encode[(int)refalt[0]];
    unsigned alt_i = acgt_encode[(int)refalt[2]];
    
    counts[ref_i] = sorted_counts[0];
    counts[alt_i] = sorted_counts[1];
    
    unsigned sc_i=2;
    for (unsigned i=0; i <= 3; i++)
        if (ref_i != i && alt_i != i) counts[i] = sorted_counts[sc_i++];

    bufprintf (vb, &vb->txt_data, "%u,%u,%u,%u", counts[0], counts[1], counts[2], counts[3]);

done:
    return HAS_NEW_VALUE;
}

// currently used only for CountBases - reverses the vector in case of XSTRAND
TRANSLATOR_FUNC (vcf_piz_luft_XREV)
{
    if (validate_only) return true; // always possible

    if (IS_TRIVAL_FORMAT_SUBFIELD) return true; // This is FORMAT field which is empty or "." - all good

    char recon_copy[recon_len];
    memcpy (recon_copy, recon, recon_len);

    str_split_enforce (recon_copy, recon_len, 0, ',', item, true, "vcf_piz_luft_XREV");

    // re-reconstruct in reverse order
    Ltxt -= recon_len;
    for (int i=n_items-1; i >= 1; i--)
        RECONSTRUCT_SEP (items[i], item_lens[i], ',');
        
    RECONSTRUCT (items[0], item_lens[0]);

    return true;
}

// -------
// INFO/AC
// -------

static void vcf_seg_INFO_AC (VBlockVCFP vb, ContextP ac_ctx, STRp(ac_str))
{
    // case: AC = AN * AF (might not be, due to rounding errors, esp if AF is a very small fraction)
    if (ctx_has_value_in_line_(vb, CTX(INFO_AC)) && // AC exists and is a valid int (set in vcf_seg_info_subfields if AC is a single int)
        ctx_has_value_in_line_(vb, CTX(INFO_AN)) && // AN exists and is a valid int
        ctx_has_value_in_line_(vb, CTX(INFO_AF)) && // AF exists and is a valid float
        (int64_t)round (CTX(INFO_AF)->last_value.f * CTX(INFO_AN)->last_value.i) == ac_ctx->last_value.i) { // AF * AN == AC

        ac_ctx->no_stons = true;
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_AC }), 2, ac_ctx, ac_str_len);
    }

    // GIAB: AC = AC_Hom + AC_Het + AC_Hemo
    else if (ctx_has_value_in_line_(vb, CTX(INFO_AC_Hom)) &&
             ctx_has_value_in_line_(vb, CTX(INFO_AC_Het)) &&
             ctx_has_value_in_line_(vb, CTX(INFO_AC_Hemi)) &&
             ac_ctx->last_value.i == CTX(INFO_AC_Hom)->last_value.i + CTX(INFO_AC_Het)->last_value.i + CTX(INFO_AC_Hemi)->last_value.i) {
        ac_ctx->no_stons = true;
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_AC, '1' }), 3, ac_ctx, ac_str_len);
    }

    // case: AC is multi allelic, or an invalid value, or missing AN or AF or AF*AN != AC
    else 
        seg_by_ctx (VB, STRa(ac_str), ac_ctx, ac_str_len);
}

// reconstruct: AC = AN * AF
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_AC)
{
    // Backward compatability note: In files v6->11, snip has 2 bytes for AN, AF which mean: '0'=appears after AC, '1'=appears before AC. We ignore them.

    // note: update last_value too, so its available to vcf_piz_luft_A_AN, which is called becore last_value is updated
    if (!snip_len || !VER(15))
        ctx->last_value.i = new_value->i = (int64_t)round (reconstruct_peek (vb, CTX(INFO_AN), 0, 0).i * reconstruct_peek(vb, CTX(INFO_AF), 0, 0).f);

    else if (*snip == '1') 
        ctx->last_value.i = new_value->i = 
            reconstruct_peek (vb, CTX(INFO_AC_Hom),  0, 0).i +
            reconstruct_peek (vb, CTX(INFO_AC_Het),  0, 0).i +
            reconstruct_peek (vb, CTX(INFO_AC_Hemi), 0, 0).i;

    else 
        ABORT_PIZ ("unrecognized snip '%c'(%u). %s", *snip, (uint8_t)*snip, genozip_update_msg());

    if (reconstruct) RECONSTRUCT_INT (new_value->i); 

    return HAS_NEW_VALUE;
}

// Lift-over translator for INFO/AC fields, IF it is bi-allelic and we have a ALT<>REF switch AND we have an AN and AF field
// We change the AC to AN - AC and return true if successful 
TRANSLATOR_FUNC (vcf_piz_luft_A_AN)
{
    if (IS_TRIVAL_FORMAT_SUBFIELD) return true; // This is FORMAT field which is empty or "." - all good

    int64_t an, ac;

    if (IS_ZIP) {
        if (!ctx_has_value_in_line_(vb, CTX(INFO_AN)) ||
            !str_get_int_range64 (STRa(recon), 0, (an = CTX(INFO_AN)->last_value.i), &ac))
            return false; // in Seg, AC is always segged last, so this means there is no AN in the line for sure
        
        if (validate_only) return true;  // Yay! AC can be lifted - all it needs in AN in the line, which it has 
    }
    else {
        ac = ctx->last_value.i;
        an = reconstruct_peek (vb, CTX(INFO_AN), 0, 0).i;
    }

    // re-reconstruct: AN-AC
    Ltxt -= recon_len;
    RECONSTRUCT_INT (an - ac);

    return true;    
}

// ------------------------
// INFO/SVLEN & INFO/REFLEN
// ------------------------

static inline void vcf_seg_INFO_SVLEN (VBlockVCFP vb, ContextP ctx, STRp(svlen_str))
{
    int64_t svlen;
    if (!str_get_int (STRa(svlen_str), &svlen)) 
        seg_by_ctx (VB, STRa(svlen_str), ctx, svlen_str_len);

    // if SVLEN is negative, it is expected to be minus the delta between END and POS
    else if (-svlen == CTX(VCF_POS)->last_delta) // INFO_END is an alias of POS - so the last delta would be between END and POS
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_SVLEN }), 2, ctx, svlen_str_len);

    // for left-anchored deletions or insertions, SVLEN might be the length of the payload
    else if (svlen == MAX_(vb->main_alt_len, vb->main_ref_len) - 1)
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_SVLEN, '1' }), 3, ctx, svlen_str_len);

    else
        seg_integer_or_not (VB, ctx, STRa(svlen_str), svlen_str_len);
}

static inline void vcf_seg_INFO_REFLEN (VBlockVCFP vb, ContextP ctx, STRp(reflen_str)) // note: ctx is INFO/END *not* POS (despite being an alias)
{
    int64_t reflen;

    if (CTX(VCF_POS)->last_delta && str_get_int (STRa(reflen_str), &reflen) && reflen == CTX(VCF_POS)->last_delta)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_SVLEN, '2' }, 3, ctx, reflen_str_len);
    else
        seg_by_ctx (VB, STRa(reflen_str), ctx, reflen_str_len);
}


// the case where SVLEN is minus the delta between END and POS
SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_SVLEN)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    
    if (!snip_len) 
        new_value->i = -CTX(VCF_POS)->last_delta; // END is a alias of POS - they share the same data stream - so last_delta would be the delta between END and POS

    else if (*snip == '2') // introduced 15.0.13
        new_value->i = CTX(VCF_POS)->last_delta;

    else if (*snip == '1') // introduced 15.0.13
        new_value->i = MAX_(vb->main_alt_len, vb->main_ref_len) - 1;

    else
        ABORT_PIZ ("unrecognized snip '%c'(%u). %s", *snip, (uint8_t)*snip, genozip_update_msg());

    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// -----------
// INFO/SVTYPE
// -----------

static inline bool vcf_seg_SVTYPE (VBlockVCFP vb, ContextP ctx, STRp(svtype))
{
    uint32_t alt_len = vb->main_alt_len;
    uint32_t ref_len = vb->main_ref_len;
    rom alt = vb->main_alt;

    // TODO: need careful testing to see this is handled correctly in case of a REF/ALT switch
    if (z_is_dvcf) goto fallback;

    // prediction: ALT has a '[' or a ']', then SVTYPE is "BND"
    else if (memchr (alt, '[', alt_len) || memchr (alt, ']', alt_len)) 
        { if (memcmp (svtype, "BND", 3)) goto fallback; }

    // prediction: if ALT starts/ends with <>, then its the same as SVTYPE except <>
    else if (alt[0] == '<' && alt[alt_len-1] == '>') 
        { if (!str_issame_(STRa(svtype), alt+1, alt_len-2)) goto fallback; }

    // prediction: if ref_len>1 and alt_len=1, then SVTYPE is <DEL>
    else if (ref_len > 1 && alt_len == 1) 
        { if (memcmp (svtype, "DEL", 3)) goto fallback; }

    // prediction: if ref_len=1 and alt_len>1, then SVTYPE is <INS>
    else if (ref_len == 1 && alt_len > 1) 
        { if (memcmp (svtype, "INS", 3)) goto fallback; }

    else fallback:
        return true;

    // prediction succeeded
    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_SVTYPE }, 2, ctx, svtype_len);
    return false;  // no need for fallback
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_SVTYPE)
{    
    rom alt = VB_VCF->main_alt;
    uint32_t alt_len = VB_VCF->main_alt_len;
    uint32_t ref_len = VB_VCF->main_ref_len;

    // prediction: ALT has a '[' or a ']', then SVTYPE is "BND"
    if (memchr (alt, '[', alt_len) || memchr (alt, ']', alt_len)) 
        RECONSTRUCT ("BND", 3);
    
    // prediction: if ALT has 2 characters more than SVTYPE, starting/end with <>, then its the same as SVTYPE except <>
    else if (alt[0] == '<' && alt[alt_len-1] == '>') 
        RECONSTRUCT (alt+1, alt_len-2);

    // prediction: if ref_len>1 and alt_len=1, then SVTYPE is <DEL>
    else if (ref_len > 1 && alt_len == 1) 
        RECONSTRUCT ("DEL", 3);

    // prediction: if ref_len=1 and alt_len>1, then SVTYPE is <INS>
    else if (ref_len == 1 && alt_len > 1) 
        RECONSTRUCT ("INS", 3);

    else
        ABORT_PIZ ("failed to reconstruct SVTYPE: ref_len=%u alt=\"%.*s\"", ref_len, STRf(alt));

    return NO_NEW_VALUE;
}    

// ##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
// comma-seperated two numbers: RAW_MQandDP=720000,200: 1. sum of squared MQ values and 2. total reads over variant genotypes (note: INFO/MQ is sqrt(#1/#2))
static inline void vcf_seg_INFO_RAW_MQandDP (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    seg_by_container (VB, ctx, (ContainerP)&RAW_MQandDP_con, STRa(value), STRa(RAW_MQandDP_snip), NULL, true, value_len);
}

// --------------
// INFO container
// --------------

// for dual coordinate files (Primary, Luft and --chain) - add DVCF depending on ostatus (run after
// all INFO and FORMAT fields, so ostatus is final)
static void vcf_seg_info_add_DVCF_to_InfoItems (VBlockVCFP vb)
{
    // case: Dual coordinates file line has no PRIM, Lrej or Prej - this can happen if variants were added to the file,
    // for example, as a result of a "bcftools merge" with a non-DVCF file
    bool added_variant = false;
    if (!ctx_encountered_in_line (VB, INFO_LUFT) && // note: no need to check PRIM because LUFT and PRIM always appear together
        !ctx_encountered_in_line (VB, INFO_LREJ) &&
        !ctx_encountered_in_line (VB, INFO_PREJ)) {
        vcf_lo_seg_rollback_and_reject (vb, LO_ADDED_VARIANT, NULL); // note: we don't report this reject because it doesn't happen during --chain
        added_variant = true; // we added a REJX field in a variant that will be reconstructed in the current coordintes
    }

    // case: line originally had LIFTOVER or LIFTBACK. These can be fields from the txt files, or created by --chain
    bool has_luft    = ctx_encountered_in_line (VB, INFO_LUFT);
    bool has_prim    = ctx_encountered_in_line (VB, INFO_PRIM);
    bool has_lrej    = ctx_encountered_in_line (VB, INFO_LREJ);
    bool has_prej    = ctx_encountered_in_line (VB, INFO_PREJ);
    bool rolled_back = LO_IS_REJECTED (last_ostatus) && (has_luft || has_prim); // rejected in the Seg process
           
    // make sure we have either both LIFT/PRIM or both Lrej/Prej subfields in Primary and Luft
    ASSVCF ((has_luft && has_prim) || (has_lrej && has_prej), "%s", 
            vb->line_coords==DC_PRIMARY ? "Missing INFO/LUFT or INFO/Lrej subfield" : "Missing INFO/PRIM or INFO/Prej subfield");

    // case: --chain and INFO is '.' - remove the '.' as we are adding a DVCF field
    if (vb->info_items.len == 1 && B1ST (InfoItem, vb->info_items)->name_len == 1 && *B1ST (InfoItem, vb->info_items)->name == '.') {
        vb->info_items.len = 0;
        vb->recon_size--;
        vb->recon_size_luft--;
    }

    // dual coordinate line - we seg both options and vcf_piz_filter will decide which to render
    if (LO_IS_OK (last_ostatus)) {

        BNXT (InfoItem, vb->info_items) = (InfoItem) { 
            .name      = INFO_LUFT_NAME"=", 
            .name_len  = INFO_DVCF_LEN + 1, // +1 for the '='
            .ctx       = CTX (INFO_LUFT),
            .value     = "" // non-zero means value exists
        };  
        
        BNXT (InfoItem, vb->info_items) = (InfoItem) { 
            .name      = INFO_PRIM_NAME"=", 
            .name_len  = INFO_DVCF_LEN + 1, // +1 for the '='
            .ctx       = CTX (INFO_PRIM),
            .value     = "" // non-zero means value exists
        };  

        // case: --chain - we're adding ONE of these subfields to each of Primary and Luft reconstructions
        if (chain_is_loaded) {
            uint32_t growth = INFO_DVCF_LEN + 1 + (vb->info_items.len32 > 2); // +1 for '=', +1 for ';' if we already have item(s)
            vb->recon_size += growth;
            vb->recon_size_luft += growth;
        }
    }

    else { 
        BNXT (InfoItem, vb->info_items) = (InfoItem) { 
            .name      = INFO_LREJ_NAME"=", 
            .name_len  = INFO_DVCF_LEN + 1, 
            .ctx       = CTX (INFO_LREJ),
            .value     = "" // non-zero means value exists
        };

        BNXT (InfoItem, vb->info_items) = (InfoItem) { 
            .name      = INFO_PREJ_NAME"=", 
            .name_len  = INFO_DVCF_LEN + 1, 
            .ctx       = CTX (INFO_PREJ),
            .value     = "" // non-zero means value exists
        };

        // case: we added a REJX INFO field that wasn't in the TXT data: --chain or rolled back (see vcf_lo_seg_rollback_and_reject) or an added variant
        if (chain_is_loaded || rolled_back || added_variant) {
            uint32_t growth = INFO_DVCF_LEN + 1 + (vb->info_items.len32 > 2); // +1 for '=', +1 for ';' if we already have item(s) execpt for the DVCF items

            if (vb->line_coords == DC_PRIMARY) 
                vb->recon_size += growth;
            else 
                vb->recon_size_luft += growth;
        }
    }

    // add tags for the DVCF info items
    if (!vb->is_rejects_vb) {
        InfoItem *ii = BLST (InfoItem, vb->info_items) - 1;
        vcf_tags_add_tag (vb, ii[0].ctx, DTYPE_VCF_INFO, ii[0].ctx->tag_name, ii[0].name_len-1);
        vcf_tags_add_tag (vb, ii[1].ctx, DTYPE_VCF_INFO, ii[1].ctx->tag_name, ii[1].name_len-1);
    }
}

static void vcf_seg_info_one_subfield (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    unsigned modified_len = value_len + 20;
    char modified[modified_len]; // used for 1. fields that are optimized 2. fields translated luft->primary. A_1 transformed 4.321e-03->0.004321
        
    // note: since we use modified for both optimization and luft_back - we currently don't support
    // subfields having both translators and optimization. This can be fixed if needed.

    ctx->line_is_luft_trans = false; // initialize
    
    #define ADJUST_FOR_MODIFIED ({                                  \
        int32_t growth = (int32_t)modified_len - (int32_t)value_len;\
        if (growth) {                                               \
            vb->recon_size      += growth;                          \
            vb->recon_size_luft += growth;                          \
        }                                                           \
        STRset (value, modified); })                                       

    // --chain: if this is RendAlg=A_1 subfield in a REF⇆ALT variant, convert a eg 4.31e-03 to e.g. 0.00431. This is to
    // ensure primary->luft->primary is lossless (4.31e-03 cannot be converted losslessly as we can't preserve format info)
    if (chain_is_loaded && ctx->luft_trans == VCF2VCF_A_1 && LO_IS_OK_SWITCH (last_ostatus) && 
        str_scientific_to_decimal (STRa(value), qSTRa(modified), NULL)) {
        ADJUST_FOR_MODIFIED;
    }        

    // Translatable item on a Luft line: attempt to lift-back the value, so we can seg it as primary
    if (vb->line_coords == DC_LUFT && needs_translation (ctx)) {

        // If cross-rendering to Primary is successful - proceed to Seg this value in primary coords, and assign a translator for reconstructing if --luft
        if (vcf_lo_seg_cross_render_to_primary (vb, ctx, STRa(value), qSTRa(modified), false)) {
            STRset (value, modified); 
            ctx->line_is_luft_trans = true; // assign translator to this item in the container, to be activated with --luft
        } 

        // This item in Luft coordinates is not translatable to primary. It is therefore a luft-only line, and we seg the remainder of it items in
        // luft coords, and only ever reconstruct it in luft (as the line with have Coord=LUFT). Since this item is already in LUFT coords
        else 
            vcf_lo_seg_rollback_and_reject (vb, LO_INFO, ctx);
    }

    // validate that the primary value (as received from caller or lifted back) can be luft-translated 
    // note: looks at snips before optimization, we're counting on the optimization not changing the validation outcome
    if (vb->line_coords == DC_PRIMARY && needs_translation (ctx)) {

        if (DT_FUNC(vb, translator)[ctx->luft_trans](VB, ctx, (char *)value, value_len, 0, true)) 
            ctx->line_is_luft_trans = true; // assign translator to this item in the container, to be activated with --luft
        else 
            REJECT_SUBFIELD (LO_INFO, ctx, ".\tCannot cross-render INFO subfield %s: \"%.*s\"", ctx->tag_name, value_len, value);            
    }

    // ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
    if (z_is_dvcf && ctx->luft_trans == VCF2VCF_ALLELE && !(segconf.vcf_is_cosmic && ctx->dict_id.num == _INFO_AA)) // INFO/AA in COSMIC is something else
        vcf_seg_INFO_allele (VB, ctx, STRa(value), 0);

    // many fields, for example: ##INFO=<ID=AN_amr_male,Number=1,Type=Integer,Description="Total number of alleles in male samples of Latino ancestry">
    else if ((ctx->tag_name[0] == 'A' && (ctx->tag_name[1] == 'N' || ctx->tag_name[1] == 'C') && (ctx->tag_name[2] == '_' || ctx->tag_name[2] == '-')) ||
             !memcmp (ctx->tag_name, "nhomalt", 7)) {
        ctx->ltype = LT_DYN_INT;
        seg_integer_or_not (VB, ctx, STRa(value), value_len);
    }

    else switch (ctx->dict_id.num) {
        #define CALL(f) ({ (f); break; })
        #define CALL_IF(cond,f)  if (cond) { (f); break; } else goto standard_seg 
        #define CALL_WITH_FALLBACK(f) if (f(vb, ctx, STRa(value))) { seg_by_ctx (VB, STRa(value), ctx, value_len); } break
        #define STORE_AND_SEG(store_type) ({ seg_set_last_txt_store_value (VB, ctx, STRa(value), store_type); seg_by_ctx (VB, STRa(value), ctx, value_len); break; })

        // ##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds of being a true variant versus being false under the trained gaussian mixture model">
        // Optimize VQSLOD
        case _INFO_VQSLOD: 
            if (flag.optimize_VQSLOD && optimize_float_2_sig_dig (STRa(value), 0, modified, &modified_len)) 
                ADJUST_FOR_MODIFIED;
            goto standard_seg;
    
        // ##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
        // END is an alias of POS - they share the same delta stream - the next POS will be a delta vs this END)
        case _INFO_END:   // alias of POS
            CALL (vcf_seg_INFO_END (vb, ctx, STRa(value)));

        case _INFO_CIEND: // alias of INFO/CIPOS
            CALL (seg_by_did (VB, STRa(value), INFO_CIPOS, value_len));

        case _INFO_SVLEN:
            CALL (vcf_seg_INFO_SVLEN (vb, ctx, STRa(value)));

        case _INFO_SVTYPE:
            CALL_WITH_FALLBACK (vcf_seg_SVTYPE);

        case _INFO_REFLEN:
            CALL (vcf_seg_INFO_REFLEN (vb, ctx, STRa(value)));

        // ##INFO=<ID=BaseCounts,Number=4,Type=Integer,Description="Counts of each base">
        case _INFO_BaseCounts: 
            CALL_WITH_FALLBACK (vcf_seg_INFO_BaseCounts);
    
        // ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
        case _INFO_DP: 
            CALL (vcf_seg_INFO_DP (vb, ctx, STRa(value)));

        // Source File
        case _INFO_SF: 
            CALL_WITH_FALLBACK (vcf_seg_INFO_SF_init);

        //##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
        case _INFO_AN:
            STORE_AND_SEG (STORE_INT);

        // ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
        case _INFO_AF: 
            STORE_AND_SEG (STORE_FLOAT);

        // ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
        case _INFO_AC:
            CALL (vcf_seg_INFO_AC (vb, ctx, STRa(value))); 

        case _INFO_MLEAF:
            CALL_IF (!z_is_dvcf && // this doesn't work for DVCF (bug 680)
                     ctx_has_value_in_line_(vb, CTX(INFO_AF)) && str_issame_(STRa(value), STRtxtw(CTX(INFO_AF)->last_txt)), 
                     seg_by_ctx (VB, STRa(af_snip), ctx, value_len)); // copy AF

        case _INFO_MLEAC:
            CALL_IF (!z_is_dvcf && // this doesn't work for DVCF (bug 680)
                     ctx_has_value_in_line_(vb, CTX(INFO_AC)),
                     seg_delta_vs_other (VB, ctx, CTX(INFO_AC), STRa(value)));

        // ##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
        case _INFO_QD:
            seg_set_last_txt (VB, CTX(INFO_QD), STRa(value)); break;

        // ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele"> 
        case _INFO_AA: // But in COSMIC, INFO/AA is something entirely different
            CALL_IF (!segconf.vcf_is_cosmic, vcf_seg_INFO_allele (VB, ctx, STRa(value), 0));

        // ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral">
        // Originating from the VEP software
        case _INFO_CSQ:
        case _INFO_vep:
            CALL_IF (segconf.vcf_is_vep, vcf_seg_INFO_CSQ (vb, ctx, STRa(value)));
        
        // ##INFO=<ID=DP_HIST,Number=R,Type=String,Description="Histogram for DP; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">
        // ##INFO=<ID=GQ_HIST,Number=R,Type=String,Description="Histogram for GQ; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">
        // ##INFO=<ID=AGE_HISTOGRAM_HET,Number=A,Type=String,Description="Histogram of ages of allele carriers; Bins: <30|30|35|40|45|50|55|60|65|70|75|80+">
        // ##INFO=<ID=AGE_HISTOGRAM_HOM,Number=A,Type=String,Description="Histogram of ages of homozygous allele carriers; Bins: <30|30|35|40|45|50|55|60|65|70|75|80+">
        case _INFO_DP_HIST:
        case _INFO_GQ_HIST:
        case _INFO_AGE_HISTOGRAM_HET:
        case _INFO_AGE_HISTOGRAM_HOM: 
            CALL (seg_uint32_matrix (VB, ctx, ctx->did_i, STRa(value), ',', '|', false, value_len));

        case _INFO_DP4:
            CALL (vcf_seg_INFO_DP4 (vb, ctx, STRa(value)));

        // ##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
        case _INFO_CLNDN:
            CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), '|', 0, false, STORE_NONE, DICT_ID_NONE, value_len));

        // ##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
        case _INFO_CLNHGVS: // ClinVar & dbSNP
        case _INFO_HGVSG:   // COSMIC & Mastermind
            if (segconf.vcf_is_mastermind) CALL (vcf_seg_mastermind_HGVSG (vb, ctx, STRa(value)));
            else                           CALL (vcf_seg_INFO_HGVS (VB, ctx, STRa(value), 0)); 

        // ##INFO=<ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
        // example: CPIC:0b3ac4db1d8e6e08a87b6942|CPIC:647d4339d5c1ddb78daff52f|CPIC:9968ce1c4d35811e7175cd29|CPIC:PA166160951|CPIC:c6c73562e2b9e4ebceb0b8bc
        // I tried seg_array_of_struct - it is worse than simple seg

        // ##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
        case _INFO_ALLELEID:
            CALL (seg_integer_or_not (VB, ctx, STRa(value), value_len));

        case _INFO_RSID:
            CALL (seg_id_field (VB, ctx, STRa(value), false, value_len));

        // SnpEff
        case _INFO_ANN: CALL (vcf_seg_INFO_ANN (vb, ctx, STRa(value)));
        case _INFO_EFF: CALL (vcf_seg_INFO_EFF (vb, ctx, STRa(value)));

        // ##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
        case _INFO_RAW_MQandDP:
            CALL (vcf_seg_INFO_RAW_MQandDP (vb, ctx, STRa(value)));

        // ##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
        case _INFO_HaplotypeScore:
            CALL (seg_float_or_not (VB, ctx, STRa(value), value_len));

        // ICGC
        case _INFO_mutation:        CALL_IF (segconf.vcf_is_icgc, vcf_seg_INFO_mutation (vb, ctx, STRa(value)));
        case _INFO_CONSEQUENCE:     CALL_IF (segconf.vcf_is_icgc, seg_array (VB, ctx, INFO_CONSEQUENCE, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));
        case _INFO_OCCURRENCE:      CALL_IF (segconf.vcf_is_icgc, seg_array (VB, ctx, INFO_OCCURRENCE,  STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));

        // COSMIC
        case _INFO_LEGACY_ID:       CALL_IF (segconf.vcf_is_cosmic, vcf_seg_INFO_LEGACY_ID   (vb, ctx, STRa(value)));
        case _INFO_SO_TERM:         CALL_IF (segconf.vcf_is_cosmic, vcf_seg_INFO_SO_TERM     (vb, ctx, STRa(value)));

        // Mastermind
        case _INFO_MMID3:           CALL_IF (segconf.vcf_is_mastermind, vcf_seg_INFO_MMID3   (vb, ctx, STRa(value)));
        case _INFO_MMURI3:          CALL_IF (segconf.vcf_is_mastermind, vcf_seg_INFO_MMURI3  (vb, ctx, STRa(value)));
        case _INFO_MMURI:           CALL_IF (segconf.vcf_is_mastermind, seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));
        case _INFO_GENE:            CALL_IF (segconf.vcf_is_mastermind, STORE_AND_SEG (STORE_NONE)); // consumed by vcf_seg_INFO_MMID3

        // Illumina genotyping
        case _INFO_PROBE_A:         CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_PROBE_A      (vb, ctx, STRa(value)));
        case _INFO_PROBE_B:         CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_PROBE_B      (vb, ctx, STRa(value)));
        case _INFO_ALLELE_A:        CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ALLELE_A     (vb, ctx, STRa(value)));
        case _INFO_ALLELE_B:        CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ALLELE_B     (vb, ctx, STRa(value)));
        case _INFO_ILLUMINA_CHR:    CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ILLUMINA_CHR (vb, ctx, STRa(value)));
        case _INFO_ILLUMINA_POS:    CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ILLUMINA_POS (vb, ctx, STRa(value)));
        case _INFO_ILLUMINA_STRAND: CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ILLUMINA_STRAND (vb, ctx, STRa(value)));
        case _INFO_refSNP:          CALL_IF (segconf.vcf_illum_gtyping, seg_id_field         (VB, ctx, STRa(value), false, value_len));

        // dbSNP
        case _INFO_dbSNPBuildID:    CALL_IF (segconf.vcf_is_dbSNP, seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_RS:              CALL_IF (segconf.vcf_is_dbSNP, vcf_seg_INFO_RS (vb, ctx, STRa(value)));
        case _INFO_RSPOS:           CALL_IF (segconf.vcf_is_dbSNP, vcf_seg_INFO_RSPOS (vb, ctx, STRa(value)));
        case _INFO_GENEINFO:        CALL_IF (segconf.vcf_is_dbSNP, seg_array (VB, ctx, INFO_GENEINFO, STRa(value), '|', 0, false, STORE_NONE, DICT_ID_NONE, value_len));
        case _INFO_VC:              CALL_IF (segconf.vcf_is_dbSNP, vcf_seg_INFO_VC (vb, ctx, STRa(value)));
        case _INFO_FREQ:            CALL_IF (segconf.vcf_is_dbSNP, seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));
        // case _INFO_TOPMED: // better leave as simple snip as the items are allele frequencies which are correleted

        // dbNSFP
        case _INFO_Polyphen2_HDIV_score : 
        case _INFO_PUniprot_aapos       :
        case _INFO_SiPhy_29way_pi       :
            CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));

        case _INFO_VEST3_score  :
        case _INFO_FATHMM_score :
            CALL (seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));

        // gnomAD
        case _INFO_age_hist_het_bin_freq:
        //case _INFO_age_hist_hom_bin_freq: // same dict_id as _INFO_age_hist_het_bin_freq
        case _INFO_gq_hist_alt_bin_freq:
        //case _INFO_gq_hist_all_bin_freq:
        case _INFO_dp_hist_alt_bin_freq:
        //case _INFO_dp_hist_all_bin_freq:
        case _INFO_ab_hist_alt_bin_freq:
            CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), '|', 0, false, STORE_INT, DICT_ID_NONE, value_len));

        // VAGrENT
        case _INFO_VD:              CALL_IF (segconf.vcf_is_vagrent, vcf_seg_INFO_VD (vb, ctx, STRa(value)));
        case _INFO_VW:              CALL_IF (segconf.vcf_is_vagrent, vcf_seg_INFO_VW (vb, ctx, STRa(value)));

        // IsaacVariantCaller / starling
        case _INFO_RU:              CALL_IF (segconf.vcf_is_isaac, vcf_seg_INFO_RU (vb, ctx, STRa(value)));
        case _INFO_REFREP:          CALL_IF (segconf.vcf_is_isaac, seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_IDREP:           CALL_IF (segconf.vcf_is_isaac, vcf_seg_INFO_IDREP (vb, ctx, STRa(value)));
        case _INFO_CSQT:            CALL_IF (segconf.vcf_is_isaac, seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));
        case _INFO_cosmic:          CALL_IF (segconf.vcf_is_isaac, seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));

        // Ultima Genomics Deep Variant
        case _INFO_X_LM:
        case _INFO_X_RM:            CALL_IF (segconf.vcf_is_ultima_dv, vcf_seg_INFO_X_LM_RM (vb, ctx, STRa(value)));
        
        // manta
        // case _INFO_LEFT_SVINSSEQ: 
        // case _INFO_RIGHT_SVINSSEQ: // tried ACGT, better off without

        default: standard_seg:
            seg_by_ctx (VB, STRa(value), ctx, value_len);
            
            if (ctx->flags.store == STORE_INT) {
                int64_t val;
                if (str_get_int (STRa(value), &val))
                    ctx_set_last_value (VB, ctx, val);
            }
    }

    ctx_set_encountered (VB, ctx);
}

static SORTER (sort_by_subfield_name)
{ 
    InfoItem *ina = (InfoItem *)a;
    InfoItem *inb = (InfoItem *)b;
    
    return strncmp (ina->name, inb->name, MIN_(ina->name_len, inb->name_len));
}

void vcf_seg_info_subfields (VBlockVCFP vb, STRp(info))
{
    vb->info_items.len = 0; // reset from previous line

    // case: INFO field is '.' (empty) (but not in DVCF as we will need to deal with DVCF items)
    if (!z_is_dvcf && IS_PERIOD (info) && !segconf.vcf_is_isaac) { // note: in Isaac, it slightly better to mux the "."
        seg_by_did (VB, ".", 1, VCF_INFO, 2); // + 1 for \t or \n
        return;
    }

    // parse the info string
    str_split (info, info_len, MAX_FIELDS-2, ';', pair, false); // -2 - leave room for LUFT + PRIM
    ASSVCF (n_pairs, "Too many INFO subfields, Genozip supports up to %u", MAX_FIELDS-2);

    buf_alloc (vb, &vb->info_items, 0, n_pairs + 2, InfoItem, CTX_GROWTH, "info_items");

    int ac_i = -1; 
    InfoItem lift_ii = {}, rejt_ii = {};

    // pass 1: initialize info items + get indices of AC, and the DVCF items
    for (unsigned i=0; i < n_pairs; i++) {
        rom equal_sign = memchr (pairs[i], '=', pair_lens[i]);
        unsigned name_len = (unsigned)(equal_sign - pairs[i]); // nonsense if no equal sign
        unsigned tag_name_len = equal_sign ? name_len : pair_lens[i];

        InfoItem ii = { .name_len  = equal_sign ? name_len + 1 : pair_lens[i], // including the '=' if there is one
                        .value     = equal_sign ? equal_sign + 1 : NULL,
                        .value_len = equal_sign ? pair_lens[i] - name_len - 1 : 0  };
        memcpy (ii.name, pairs[i], ii.name_len); // note: we make a copy the name, because vcf_seg_FORMAT_GT might overwrite the INFO field
        
        // create context if it doesn't already exist (also verifies tag is not too long)
        DictId dict_id = dict_id_make (pairs[i], tag_name_len, DTYPE_1);
        ii.ctx = ctx_get_ctx_tag (vb, dict_id, pairs[i], tag_name_len); // create if it doesn't already exist
        
        if (z_is_dvcf && !vb->is_rejects_vb) vcf_tags_add_tag (vb, ii.ctx, DTYPE_VCF_INFO, pairs[i], tag_name_len);

        if (segconf.running) segconf.has[ii.ctx->did_i]++;

        ASSVCF (!z_is_dvcf || 
                  (((dict_id.num != _INFO_LUFT && dict_id.num != _INFO_LREJ) || vb->line_coords == DC_PRIMARY) && 
                   ((dict_id.num != _INFO_PRIM && dict_id.num != _INFO_PREJ) || vb->line_coords == DC_LUFT)),
                "Not expecting INFO/%.*s in a %s-coordinate line", tag_name_len, pairs[i], vcf_coords_name (vb->line_coords));

        if (dict_id.num == _INFO_AC) {
            int64_t ac;
            if (str_get_int (STRa(ii.value), &ac))
                ctx_set_last_value (VB, ii.ctx, ac); // needed for MLEAC        

            ac_i = vb->info_items.len;
        }

        else if (dict_id.num == _INFO_LUFT || dict_id.num == _INFO_PRIM) 
            { lift_ii = ii; continue; } // dont add LUFT and PRIM to Items yet

        else if (dict_id.num == _INFO_LREJ || dict_id.num == _INFO_PREJ) 
            { rejt_ii = ii; continue; } // dont add Lrej and Prej to Items yet

        BNXT (InfoItem, vb->info_items) = ii;
    }

    // case: we have a LUFT or PRIM item - Seg it now, but don't add it yet to InfoItems
    if (lift_ii.value) { 
        vcf_lo_seg_INFO_LUFT_and_PRIM (vb, lift_ii.ctx, lift_ii.value, lift_ii.value_len); 

        // case: we have both LIFT and REJT - could happen as a result of bcftools merge - discard the REJT for now, and let our Seg
        // decide if to reject it
        if (rejt_ii.value) {
            uint32_t shrinkage = rejt_ii.name_len + rejt_ii.value_len + 1; // unaccount for name, value and  
            vb->recon_size      -= shrinkage;
            vb->recon_size_luft -= shrinkage; // since its read from TXT, it is accounted for initialially in both recon_size and recon_size_luft
        }

        // case: line was reject - PRIM/LUFT changed to REJx (note: vcf_lo_seg_rollback_and_reject didn't decrement recon_size* in this case, bc PRIM/LUFT was not "encountered" yet)
        if (LO_IS_REJECTED (last_ostatus)) {
            vb->recon_size      -= lift_ii.value_len;
            vb->recon_size_luft -= lift_ii.value_len; // since its read from TXT, it is accounted for initialially in both recon_size and recon_size_luft
        }
    }
        
    // case: we have a *rej item - Seg it now, but don't add it yet to InfoItems
    else if (rejt_ii.value)
        vcf_lo_seg_INFO_REJX (vb, rejt_ii.ctx, rejt_ii.value, rejt_ii.value_len); 

    ARRAY (InfoItem, ii, vb->info_items);

    // pass 2: seg all subfields except AC (and PRIM/LUFT that weren't added)
    for (unsigned i=0; i < ii_len; i++) 
        if (ii[i].ctx->dict_id.num != _INFO_AC)
            vcf_seg_info_one_subfield (vb, ii[i].ctx, ii[i].value, ii[i].value_len);
    
    // last, seg AC (delayed, as we needed to seg AN and AF before)
    if (ac_i >= 0) 
        vcf_seg_info_one_subfield (vb, ii[ac_i].ctx, ii[ac_i].value, ii[ac_i].value_len);
}

// Adds the DVCF items according to ostatus, finalizes INFO/SF, INFO/QD and segs the INFO container
void vcf_finalize_seg_info (VBlockVCFP vb)
{
    if (!vb->info_items.len && !z_is_dvcf) return; // no INFO items on this line (except if dual-coords - we will add them in a sec)

    Container con = { .repeats             = 1, 
                      .drop_final_item_sep = true,
                      .filter_items        = z_is_dvcf,   // vcf_piz_filter chooses which (if any) DVCF item to show based on flag.luft and flag.single_coord
                      .callback            = z_is_dvcf }; // vcf_piz_container_cb appends oSTATUS to INFO if requested 
 
    // seg INFO/SF, if there is one
    if (vb->sf_txt.len) vcf_seg_INFO_SF_seg (vb);

    // seg INFO/DP, if against sum of FORMAT/DP
    if (CTX(INFO_DP)->dp.by_format_dp) 
        vcf_seg_INFO_DP_by_FORMAT_DP (vb);

    if (ctx_encountered_in_line (VB, INFO_QD))
        vcf_seg_INFO_QD (vb);

    // now that we segged all INFO and FORMAT subfields, we have the final ostatus and can add the DVCF items
    if (z_is_dvcf)
        vcf_seg_info_add_DVCF_to_InfoItems (vb);

    ARRAY (InfoItem, ii, vb->info_items);

    con_set_nitems (con, ii_len);

    // if requested, we will re-sort the info fields in alphabetical order. This will result less words in the dictionary
    // thereby both improving compression and improving --regions speed. 
    if (flag.optimize_sort && ii_len > 1) 
        qsort (ii, ii_len, sizeof(InfoItem), sort_by_subfield_name);

    char prefixes[CONTAINER_MAX_PREFIXES_LEN];  // these are the Container prefixes
    prefixes[0] = prefixes[1] = CON_PX_SEP; // initial CON_PX_SEP follow by separator of empty Container-wide prefix
    unsigned prefixes_len = 2;

    // Populate the Container 
    uint32_t total_names_len=0;
    for (unsigned i=0; i < ii_len; i++) {
        // Set the Container item and find (or create) a context for this name
        con.items[i] = (ContainerItem){ .dict_id   = !ii[i].value                        ? DICT_ID_NONE 
                                                   : ii[i].ctx->dict_id.num == _INFO_END ? (DictId)_VCF_POS
                                                   :                                       ii[i].ctx->dict_id,
                                        .separator = { ';' } }; 

        // if we're preparing a dual-coordinate VCF and this line needs translation to Luft - assign the liftover-translator for this item,
        if (ii[i].ctx && ii[i].ctx->line_is_luft_trans) { // item was segged in Primary coords and needs a luft translator to be reconstruced in --luft
            con.items[i].translator = ii[i].ctx->luft_trans;

            if (ii[i].ctx->luft_trans == VCF2VCF_A_AN)
                ii[i].ctx->flags.store = STORE_INT; // consumed by vcf_piz_luft_A_AN
        }
            
        // add to the prefixes
        ASSVCF (prefixes_len + ii[i].name_len + 1 <= CONTAINER_MAX_PREFIXES_LEN, 
                "INFO contains tag names that, combined (including the '='), exceed the maximum of %u characters", CONTAINER_MAX_PREFIXES_LEN);

        memcpy (&prefixes[prefixes_len], ii[i].name, ii[i].name_len);
        prefixes_len += ii[i].name_len;
        prefixes[prefixes_len++] = CON_PX_SEP;

        // don't include LIFTBACK or LIFTREJT because they are not reconstructed by default (genounzip) 
        // note: vcf_lo_seg_INFO_REJX / vcf_lo_seg_INFO_LUFT_and_PRIM already verified that this is a dual-coord file
        if (!ii[i].ctx || (ii[i].ctx->dict_id.num != _INFO_PRIM && ii[i].ctx->dict_id.num != _INFO_PREJ))
            total_names_len += ii[i].name_len + 1; // +1 for ; \t or \n separator
    }

    // --chain: if any tags need renaming we create a second, renames, prefixes string
    char ren_prefixes[con_nitems(con) * MAX_TAG_LEN]; 
    unsigned ren_prefixes_len = z_is_dvcf && !vb->is_rejects_vb ? vcf_tags_rename (vb, con_nitems(con), 0, 0, 0, B1ST (InfoItem, vb->info_items), ren_prefixes) : 0;

    // case GVCF: multiplex by has_RGQ or FILTER in Isaac
    if (!segconf.running && (segconf.has[FORMAT_RGQ] || segconf.vcf_is_isaac)) {
        ContextP channel_ctx = 
            seg_mux_get_channel_ctx (VB, VCF_INFO, (MultiplexerP)&vb->mux_INFO, (segconf.has[FORMAT_RGQ] ? vb->line_has_RGQ : vcf_isaac_info_channel_i (VB)));
        
        seg_by_did (VB, STRa(vb->mux_INFO.snip), VCF_INFO, 0);

        // if we're compressing a Luft rendition, swap the prefixes
        if (vb->line_coords == DC_LUFT && ren_prefixes_len) 
            container_seg_with_rename (vb, channel_ctx, &con, ren_prefixes, ren_prefixes_len, prefixes, prefixes_len, total_names_len /* names inc. = and separator */, NULL);
        else 
            container_seg_with_rename (vb, channel_ctx, &con, prefixes, prefixes_len, ren_prefixes, ren_prefixes_len, total_names_len /* names inc. = and separator */, NULL);
    }

    // case: not GVCF
    else {
        // if we're compressing a Luft rendition, swap the prefixes
        if (vb->line_coords == DC_LUFT && ren_prefixes_len) 
            container_seg_with_rename (vb, CTX(VCF_INFO), &con, ren_prefixes, ren_prefixes_len, prefixes, prefixes_len, total_names_len /* names inc. = and separator */, NULL);
        else 
            container_seg_with_rename (vb, CTX(VCF_INFO), &con, prefixes, prefixes_len, ren_prefixes, ren_prefixes_len, total_names_len /* names inc. = and separator */, NULL);
    }
}

