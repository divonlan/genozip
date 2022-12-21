// ------------------------------------------------------------------
//   vcf_info.c
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "optimize.h"
#include "file.h"
#include "strings.h"
#include "dict_id.h"
#include "codec.h"
#include "reconstruct.h"
#include "gff.h"
#include "stats.h"

static inline bool vcf_is_use_DP_by_DP (void); // forward
sSTRl(RAW_MQandDP_snip, 64 + 2 * 16);     
sSTRl(copy_ID_snip, 32);
sSTRl(copy_POS_snip, 32);

static SmallContainer RAW_MQandDP_con = {
    .nitems_lo = 2, 
    .repeats   = 1, 
    .items     = { { .dict_id={ _INFO_RAW_MQandDP_MQ }, .separator = {','} },
                   { .dict_id={ _INFO_RAW_MQandDP_DP }                     } }
};

void vcf_info_zip_initialize (void) 
{
    container_prepare_snip ((ContainerP)&RAW_MQandDP_con, 0, 0, qSTRa(RAW_MQandDP_snip));

    seg_prepare_snip_other (SNIP_OTHER_DELTA, _VCF_ID,  true, 0, copy_ID_snip);
    seg_prepare_snip_other (SNIP_OTHER_DELTA, _VCF_POS, true, 0, copy_POS_snip);
}

void vcf_info_seg_initialize (VBlockVCFP vb) 
{
    ctx_set_store (VB, STORE_INT, INFO_AN, INFO_AC, INFO_ADP, INFO_DP, INFO_MLEAC, DID_EOL);
    CTX(INFO_AF)->     flags.store = STORE_FLOAT;
    // xxx (is this really needed for --indels-only?) CTX(INFO_SVTYPE)-> flags.store = STORE_INDEX; // since v13 - consumed by vcf_refalt_piz_is_variant_indel

    if (!vcf_is_use_DP_by_DP())
        CTX(INFO_DP)->ltype = LT_DYN_INT; 

    seg_id_field_init (CTX(INFO_CSQ_Existing_variation));

    ctx_consolidate_stats (VB, INFO_RAW_MQandDP, INFO_RAW_MQandDP_MQ, INFO_RAW_MQandDP_DP, DID_EOL);
}

//--------
// INFO/DP
// -------

static inline bool vcf_is_use_DP_by_DP (void)
{
    // note: we restrict to best, as this method would cause genozip --drop-genotypes, --GT-only, --samples to show INFO/DP=-1
    return vcf_num_samples > 1 && flag.best;
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

    if (has_value)
        ctx_set_last_value (VB, ctx, value);
}

// used for multi-sample VCFs, IF FORMAT/DP is segged as a simple integer
static void vcf_seg_INFO_DP_by_FORMAT_DP (VBlockVCFP vb)
{
    ContextP ctx = CTX(INFO_DP);

    int value_len = str_int_len (ctx->last_value.i);

    SNIP(32) = { SNIP_SPECIAL, VCF_SPECIAL_DP_by_DP };
    snip_len = 2 + str_int (value_len, &snip[2]);
    snip[snip_len++] = '\t';

    // note: INFO/DP >= sum(FORMAT/DP) as the per-sample value is filtered, see: https://gatk.broadinstitute.org/hc/en-us/articles/360036891012-DepthPerSampleHC
    snip_len += str_int (ctx->last_value.i - ctx->sum_dp_this_line, &snip[snip_len]);

    seg_by_ctx (VB, STRa(snip), ctx, value_len); 
}

// initialize reconstructing INFO/DP by sum(FORMAT/DP) - save space in txt_data, and initialize delta
SPECIAL_RECONSTRUCTOR (vcf_piz_special_DP_by_DP)
{
    str_split (snip, snip_len, 2, '\t', item, 2);

    if (reconstruct) {
        if (!flag.drop_genotypes && !flag.gt_only && !flag.samples) {
            vb->txt_data.len32   += atoi (items[0]); // number of characters needed to reconstruct the INFO/DP integer
            ctx->sum_dp_this_line = atoi (items[1]); // initialize with delta
            ctx->is_initialized = true;              // needs to be finalized
        }
        else
            RECONSTRUCT ("-1", 2); // bc we can't calculate INFO/DP in these cases bc we need FORMAT/DP of all samples
    }

    return NO_NEW_VALUE; 
}

// finalize reconstructing INFO/DP by sum(FORMAT/DP) - called after reconstructing all samples
void vcf_piz_finalize_DP_by_DP (VBlockVCFP vb)
{
    str_int_ex (CTX(INFO_DP)->sum_dp_this_line, last_txt(VB, INFO_DP), false);
}

// used starting v13.0.5, replaced in v14 with a new vcf_piz_special_DP_by_DP
SPECIAL_RECONSTRUCTOR (vcf_piz_special_DP_by_DP_v13)
{
    str_split (snip, snip_len, 2, '\t', item, 2);

    int num_dps_this_line = atoi (items[0]);
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

// -------
// INFO/AA
// -------

// return 0 if the allele equals main REF, the alt number if it equals one of the ALTs, or -1 if none or -2 if '.'
static int vcf_INFO_ALLELE_get_allele (VBlockVCFP vb, STRp (value))
{
    // check for '.'
    if (value_len == 1 && *value == '.') return -2;

    // check if its equal main REF (which can by REF or oREF)
    if (value_len == vb->main_ref_len && !memcmp (value, vb->main_ref, value_len))
        return 0;

    // check if its equal one of the ALTs
    str_split (vb->main_alt, vb->main_alt_len, 0, ',', alt, false);

    for (int alt_i=0; alt_i < n_alts; alt_i++) 
        if (value_len == alt_lens[alt_i] && !memcmp (value, alts[alt_i], value_len)) 
            return alt_i + 1;

    // case: not REF or any of the ALTs
    return -1;
}

// checks if value is identifcal to the REF or one of the ALT alleles, and if so segs a SPECIAL snip
// Used for INFO/AA, INFO/CSQ/Allele, INFO/ANN/Allele. Any field using this should have the VCF2VCF_ALLELE translator set in vcf_lo_luft_trans_id.
static bool vcf_seg_INFO_allele (VBlockP vb_, ContextP ctx, STRp(value), uint32_t repeat) // returns true if caller still needs to seg 
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

    rom refalt = last_txtx (vb, refalt_ctx);
    unsigned refalt_len = refalt_ctx->last_txt.len;

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
    if (validate_only && last_ostatus == LO_OK_REF_NEW_SNP && 
        recon_len == vcf_vb->main_ref_len && !memcmp (recon, vcf_vb->main_ref, recon_len)) return false;

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
    vb->txt_data.len -= recon_len;
    for (int i=n_items-1; i >= 1; i--)
        RECONSTRUCT_SEP (items[i], item_lens[i], ',');
        
    RECONSTRUCT (items[0], item_lens[0]);

    return true;
}

// -------
// INFO/AC
// -------

static void vcf_seg_INFO_AC (VBlockVCFP vb, ContextP ac_ctx, STRp(field))
{
    // case: AC = AN * AF (might not be, due to rounding errors, esp if AF is a very small fraction)
    if (ctx_has_value_in_line_(vb, CTX(INFO_AC)) && // AC exists and is a valid int (set in vcf_seg_info_subfields if AC is a single int)
        ctx_has_value_in_line_(vb, CTX(INFO_AN)) && // AN exists and is a valid int
        ctx_has_value_in_line_(vb, CTX(INFO_AF)) && // AF exists and is a valid float
        (int64_t)round (CTX(INFO_AF)->last_value.f * CTX(INFO_AN)->last_value.i) == CTX(INFO_AC)->last_value.i) { // AF * AN == AC

        ac_ctx->no_stons = true;
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_AC }), 2, ac_ctx, field_len);
    }

    // case: AC is multi allelic, or an invalid value, or missing AN or AF or AF*AN != AC
    else 
        seg_by_ctx (VB, STRa(field), ac_ctx, field_len);
}

// reconstruct: AC = AN * AF
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_AC)
{
    if (!reconstruct) return false;

    // Backward compatability note: In files v6->11, snip has 2 bytes for AN, AF which mean: '0'=appears after AC, '1'=appears before AC. We ignore them.

    // note: update last_value too, so its available to vcf_piz_luft_A_AN, which is called becore last_value is updated
    ctx->last_value.i = new_value->i = (int64_t)round (reconstruct_peek (vb, CTX(INFO_AN), 0, 0).i * reconstruct_peek(vb, CTX(INFO_AF), 0, 0).f);
    RECONSTRUCT_INT (new_value->i); 

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
    vb->txt_data.len32 -= recon_len;
    RECONSTRUCT_INT (an - ac);

    return true;    
}

// --------
// INFO/END
// --------

static void vcf_seg_INFO_END (VBlockVCFP vb, ContextP end_ctx, rom end_str, unsigned end_len) // note: ctx is INFO/END *not* POS (despite being an alias)
{
    // END is an alias of POS
    seg_pos_field (VB, VCF_POS, VCF_POS, SPF_BAD_SNIPS_TOO | SPF_ZERO_IS_BAD | SPF_UNLIMITED_DELTA, 0, end_str, end_len, 0, end_len);

    // add end_delta to dl for sorting. it is used only in case chrom and pos are identical
    DATA_LINE (vb->line_i)->end_delta = vb->last_delta (VCF_POS);

    // case --chain: if we have lifted-over POS (as primary POS field or in INFO/LIFTBACK), 
    // check that lifting-over of END is delta-encoded and is lifted over to the same, non-xstrand, Chain alignment, and reject if not
    if (chain_is_loaded && LO_IS_OK (last_ostatus)) { 

        bool is_xstrand = (vb->last_index (VCF_oXSTRAND) > 0); // set in vcf_lo_seg_generate_INFO_DVCF
        PosType aln_last_pos = chain_get_aln_prim_last_pos (vb->pos_aln_i); 
        PosType end = vb->last_int (VCF_POS); 

        // case: we don't yet handle END translation in case of a reverse strand
        if (is_xstrand)            
            REJECT_SUBFIELD (LO_INFO, end_ctx, ".\tVariant with INFO/END and chain file alignment with a negative strand%s", "");

        // case: END goes beyond the end of the chain file alignment and its a <DEL>
        else if (vb->is_del_sv && end > aln_last_pos) {

            // case: END goes beyond end of alignment
            PosType gap_after = chain_get_aln_gap_after (vb->pos_aln_i);
            
            // case: END falls in the gap after - <DEL> is still valid but translated END needs to be closer to POS to avoid gap - 
            // we don't yet do this
            if (end <= aln_last_pos + gap_after)
                REJECT_SUBFIELD (LO_INFO, end_ctx, ".\t<DEL> variant: INFO/END=%.*s is in the gap after the end of the chain file alignment", end_len, end_str);
    
            // case: END falls on beyond the gap (next alignment or beyond) - this variant cannot be lifted
            else
                REJECT_SUBFIELD (LO_INFO, end_ctx, ".\t<DEL> variant: INFO/END=%.*s is beyond the end of the chain file alignment and also beyond the gap after the alignment", end_len, end_str);
        }

        // case: END goes beyond the end of the chain file alignment and its NOT a <DEL>
        else if (!vb->is_del_sv && end > aln_last_pos) 
            REJECT_SUBFIELD (LO_INFO, end_ctx, ".\tPOS and INFO/END=%.*s are not on the same chain file alignment", end_len, end_str);

        // case: invalid value. since we use SPF_UNLIMITED_DELTA, any integer value should succeed
        else if (!CTX(VCF_POS)->last_delta)
            REJECT_SUBFIELD (LO_INFO, end_ctx, ".\tINFO/END=%.*s has an invalid value", end_len, end_str);        
    }
}

// END data resides in POS (its an alias), but has a different translator as its a different container item. 
// For END, we didn't add an oPOS entry, because we can't consume it when showing Primary. Instead, we do delta arithmetic.
// returns true if successful (return value used only if validate_only)
TRANSLATOR_FUNC (vcf_piz_luft_END)
{
    // ZIP liftover validation: postpone to vcf_seg_INFO_END
    if (validate_only) return true; 

    PosType translated_end;
    ContextP pos_ctx  = CTX (VCF_POS);
    ContextP opos_ctx = CTX (VCF_oPOS);
    
    // ZIP liftback (POS is always before END, because we seg INFO/LIFTOVER first)
    if (IS_ZIP && VB_VCF->line_coords == DC_LUFT) { // liftback
        PosType oend;
        if (!str_get_int_range64 (STRa(recon), 0, MAX_POS, &oend))
            return false;

        translated_end = pos_ctx->last_value.i + (oend - opos_ctx->last_value.i);
    }

    // PIZ liftover: we have already reconstructed oPOS and POS (as an item in VCF_TOPLUFT)
    else {
        translated_end = opos_ctx->last_value.i + pos_ctx->last_delta; ; // delta for generated this END value (END - POS)

        CTX(INFO_END)->last_end_line_i = vb->line_i; // so vcf_piz_special_COPYPOS knows that END was reconstructed
    }

    // re-reconstruct END
    vb->txt_data.len -= recon_len; 
    RECONSTRUCT_INT (translated_end);
    
    return true;    
}

// Called to reconstruct the POS subfield of INFO/LIFTBACK, handling the possibility of INFO/END
// Also used for FORMAT/PS
SPECIAL_RECONSTRUCTOR (vcf_piz_special_COPYPOS)
{
    if (!reconstruct) return NO_NEW_VALUE;
    
    bool has_end = CTX(INFO_END)->last_end_line_i == vb->line_i; // true if INFO/END was encountered

    ContextP pos_ctx = CTX (VCF_POS);
    int64_t pos;

    if (has_end) {
        int64_t end   = pos_ctx->last_value.i;
        int64_t delta = pos_ctx->last_delta;
        pos = end - delta;
    }
    else
        pos = pos_ctx->last_value.i;

    if (snip_len)
        pos += atoi (snip); // add optional delta (since 13.0.5)

    RECONSTRUCT_INT (pos); // vcf_piz_luft_END makes sure it always contains the value of POS, not END
    return NO_NEW_VALUE;
}

// ----------
// INFO/SVLEN
// ----------

static inline bool vcf_seg_test_SVLEN (VBlockVCFP vb, STRp(svlen_str))
{
    int64_t svlen;
    if (!str_get_int (STRa(svlen_str), &svlen)) return false;

    int64_t last_delta = CTX (VCF_POS)->last_delta; // INFO_END is an alias of POS - so the last delta would be between END and POS
    return last_delta == -svlen;
}

// the case where SVLEN is minus the delta between END and POS
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_SVLEN)
{
    if (!reconstruct) goto done;

    int64_t value = -CTX (VCF_POS)->last_delta; // END is a alias of POS - they share the same data stream - so last_delta would be the delta between END and POS
    RECONSTRUCT_INT (value);

done:
    return NO_NEW_VALUE;
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

    return NO_NEW_VALUE;
}    

// ------------
// INFO/CLNHGVS
// ------------

// SNP case: "NC_000023.10:g.154507173T>G"
static bool vcf_seg_INFO_HGVS_snp (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    PosType pos = DATA_LINE (vb->line_i)->pos[0]; // data in variant
    char pos_str[30];
    unsigned pos_str_len = str_int (pos, pos_str);

    if (value_len < 3 + pos_str_len) return false;

    rom v = &value[value_len - 3];
    if (v[0] != vb->main_ref[0] || v[1] != '>' || v[2] != vb->main_alt[0]) return false; // REF/ALT differs

    v -= pos_str_len;
    if (memcmp (v, pos_str, pos_str_len)) return false; // POS differs

    SmallContainer con = { 
        .repeats   = 1,
        .nitems_lo = 2,
        .items = { { .dict_id = (DictId)_INFO_HGVS_snp_pos    },
                   { .dict_id = (DictId)_INFO_HGVS_snp_refalt } }
     }; 

    // temporarily surround prefix by separators, and seg container with prefix
    SAFE_ASSIGNx (&value[-1], CON_PX_SEP, 1);
    SAFE_ASSIGNx (v,          CON_PX_SEP, 2);

    container_seg (vb, ctx, (ContainerP)&con, &value[-1], v - value + 2, value_len - pos_str_len - 3);

    SAFE_RESTOREx(1);
    SAFE_RESTOREx(2);

    // seg special snips
    ctx_consolidate_stats (VB, ctx->did_i, INFO_HGVS_snp_pos, INFO_HGVS_snp_refalt, DID_EOL);

    seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_HGVS_SNP_POS    }), 2, CTX(INFO_HGVS_snp_pos),    pos_str_len);
    seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_HGVS_SNP_REFALT }), 2, CTX(INFO_HGVS_snp_refalt), 3);

    return true;
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_SNP_POS)
{
    ContextP pos_ctx = CTX (VCF_POS);
    
    if (VB_VCF->vb_coords == DC_PRIMARY)
        RECONSTRUCT (last_txtx (vb, pos_ctx), pos_ctx->last_txt.len); // faster than RECONSTRUCT_INT
    
    else  // if reconstructing Luft, VCF_POS just consumed and last_int set if in Luft coords (see top_luft container)
        RECONSTRUCT_INT (pos_ctx->last_value.i);
    
    return NO_NEW_VALUE;
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_SNP_REFALT)
{
    rom refalt;
    reconstruct_peek (vb, CTX (VCF_REFALT), &refalt, NULL); // this special works only on SNPs, so length is always 3

    RECONSTRUCT1 (refalt[0]); // this might overwrite the "peeked" data, but that's ok
    RECONSTRUCT1 ('>');
    RECONSTRUCT1 (refalt[2]);

    return NO_NEW_VALUE;
}

// Deletion case: "n.10571909_10571915delCCCGCCG" (POS=10571908 ; 10571909 are the bases except for the left-anchor, CCCGCCG is the payload of the indel)
//                "NC_000001.10:g.5993401_5993405del" REF=TAAAAC ALT=T POS=5993397
//                "NC_000001.10:g.6038478del" REF=AG ALT=A POS=6038476
// Insertion case: "n.10659939_10659940insTG" (POS=10659939 ; TG is the payload of the indel) 
// Delins: "NC_000001.10:g.5987727_5987729delinsCCACG" POS=5987727 REF=GTT ALT=CCACG
typedef enum { DEL, INS, DELINS } HgvsType;
static bool vcf_seg_INFO_HGVS_indel (VBlockVCFP vb, ContextP ctx, STRp(value), rom op, HgvsType t)
{
    rom payload = &op[3];
    unsigned payload_len = (unsigned)(&value[value_len] - payload);

    if (payload_len) switch (t) {
        case DEL    : // Payload is expected to be the same as the REF field, without the first, left-anchor, base
                      if (payload_len != vb->main_ref_len-1 || memcmp (payload, vb->main_ref + 1, payload_len)) return false;
                      break;
        case INS    : // Payload is expected to be the same as the ALT field, without the first, left-anchor, base
                      if (payload_len != vb->main_alt_len-1 || memcmp (payload, vb->main_alt + 1, payload_len)) return false;
                      break;
                      // Payload is expected to be the same as the entire ALT field
        case DELINS : if (payload_len != vb->main_alt_len || memcmp (payload, vb->main_alt, payload_len)) return false;
                      break;
    }

    // beginning of number is one after the '.' - scan backwards
    rom start_pos = op-1;
    while (start_pos[-1] != '.' && start_pos > value) start_pos--;
    if (start_pos[-1] != '.') return false;

    str_split (start_pos, op-start_pos, 2, '_', pos, false);

    PosType pos[2];
    if (!str_get_int (poss[0], pos_lens[0], &pos[0])) return false;
    
    if (n_poss == 2) {
        if (!str_get_int (poss[1], pos_lens[1], &pos[1])) return false;
        if (pos[0] + payload_len - 1 != pos[1]) return false;
    }
    
    static const DictId dict_id_start_pos[3] = { { _INFO_HGVS_del_start_pos }, { _INFO_HGVS_ins_start_pos }, { _INFO_HGVS_ins_start_pos  } }; // INS and DELINS use the same start_pos
    static const DictId dict_id_end_pos[3]   = { { _INFO_HGVS_del_end_pos   }, { _INFO_HGVS_ins_end_pos   }, { _INFO_HGVS_delins_end_pos } };
    static const DictId dict_id_payload[3]   = { { _INFO_HGVS_del_payload   }, { _INFO_HGVS_ins_payload   }, { _INFO_HGVS_delins_payload } };

    static const Did did_i_start_pos[3] = { INFO_HGVS_del_start_pos, INFO_HGVS_ins_start_pos, INFO_HGVS_ins_start_pos};
    static const Did did_i_end_pos[3]   = { INFO_HGVS_del_end_pos, INFO_HGVS_ins_end_pos, INFO_HGVS_delins_end_pos};
    static const Did did_i_payload[3]   = { INFO_HGVS_del_payload, INFO_HGVS_ins_payload, INFO_HGVS_delins_payload};

    static const uint8_t special_end_pos[3]  = { VCF_SPECIAL_HGVS_DEL_END_POS, VCF_SPECIAL_HGVS_INS_END_POS, VCF_SPECIAL_HGVS_DELINS_END_POS };
    static const uint8_t special_payload[3]  = { VCF_SPECIAL_HGVS_DEL_PAYLOAD, VCF_SPECIAL_HGVS_INS_PAYLOAD, VCF_SPECIAL_HGVS_DELINS_PAYLOAD };

    SmallContainer con = { 
        .repeats   = 1,
        .nitems_lo = 3,
        .items = { { .dict_id = dict_id_start_pos[t], .separator = "_" }, // separator deleted in container_reconstruct_do() if end_pos is missing
                   { .dict_id = dict_id_end_pos[t]                     },
                   { .dict_id = dict_id_payload[t]                     } } }; 

    // preper prefixes - a container-wide prefix, and the op is the prefix for item [2]
    int64_t header_len = start_pos - value;  
    int64_t op_len = (t==DELINS ? 6 : 3);
    int64_t prefixes_len = header_len + op_len + 5 /* separators */;

    // note: header_len, op_len and renamed_len prefixes_len to be int64_t to avoid -Wstringop-overflow warning in gcc 10
    char prefixes[prefixes_len];
    prefixes[0] = prefixes[header_len+1] = prefixes[header_len+2] = prefixes[header_len+3] = prefixes[header_len+4+op_len] = CON_PX_SEP;
    memcpy (&prefixes[1], value, header_len);
    memcpy (&prefixes[header_len+4], op, op_len);
    
    container_seg (vb, ctx, (ContainerP)&con, prefixes, prefixes_len, value_len - pos_lens[0] - (n_poss==2 ? pos_lens[1] : 0) - payload_len);

    // seg special snips
    ctx_consolidate_stats (VB, ctx->did_i, did_i_start_pos[t], did_i_end_pos[t], did_i_payload[t], DID_EOL);

    CTX(did_i_start_pos[t])->flags.store = STORE_INT; // consumed by vcf_piz_special_INFO_HGVS_DEL_END_POS

    seg_pos_field (VB, did_i_start_pos[t], VCF_POS, SPF_UNLIMITED_DELTA, 0, 0, 0, pos[0], pos_lens[0]);
    
    // We pos_lens[1] only if the payload is longer than 1
    if (n_poss == 2)
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, special_end_pos[t] }), 2, CTX(did_i_end_pos[t]), pos_lens[1]);
    else
        seg_by_ctx (VB, NULL, 0, CTX(did_i_end_pos[t]), 0); // becomes WORD_INDEX_MISSING - container_reconstruct_do will remove the preceding _

    // the del payload is optional - we may or may not have it
    if (payload_len)
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, special_payload[t] }), 2, CTX(did_i_payload[t]), payload_len);
    else
        seg_by_ctx (VB, NULL, 0, CTX(did_i_payload[t]), 0); 

    return true;
}

// <ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
static bool vcf_seg_INFO_HGVS (VBlockP vb_, ContextP ctx, STRp(value), uint32_t repeat)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    if (ctx_encountered_in_line (VB, INFO_END)) 
        goto fail; // we can't use this if there is an END before CLNHGVS, as it will change last_int(VCF_POS) during reconstruction

    // if main->refalt is different than REFALT then we can't use this function, as reconstruction is based on VCF_REFALT
    if (vb->line_coords == DC_LUFT && (LO_IS_OK_SWITCH (last_ostatus) || LO_IS_REJECTED (last_ostatus) || *CTX(VCF_oXSTRAND)->last_snip != '-'))
        goto fail;

    SAFE_NULT (value);

    bool success = false;
    rom op;
    if      (vb->main_ref_len == 1 && vb->main_alt_len == 1) success = vcf_seg_INFO_HGVS_snp   (vb, ctx, STRa(value));
    else if ((op = strstr (value, "delins")))                success = vcf_seg_INFO_HGVS_indel (vb, ctx, STRa(value), op, DELINS); // must be before del and ins
    else if ((op = strstr (value, "del")))                   success = vcf_seg_INFO_HGVS_indel (vb, ctx, STRa(value), op, DEL);
    else if ((op = strstr (value, "ins")))                   success = vcf_seg_INFO_HGVS_indel (vb, ctx, STRa(value), op, INS);
    
    SAFE_RESTORE;
    
    if (success) return true; // indeed segged

fail:
    seg_by_ctx (VB, STRa(value), ctx, value_len); 
    return true; // indeed segged
}

static void vcf_piz_special_INFO_HGVS_INDEL_END_POS (VBlockP vb, HgvsType t)
{
    STR (refalt);
    reconstruct_peek (vb, CTX (VCF_REFALT), pSTRa(refalt)); // this special works only on SNPs, so length is always 3

    SAFE_NULT (refalt);
    rom tab = strchr (refalt, '\t');
    ASSPIZ (tab, "invalid REF+ALT=\"%.*s\"", STRf(refalt));
    SAFE_RESTORE;

    // reconstruct the DEL payload - this is the REF except for the first (anchor) base. reconstruct with one (possibly overlapping) copy
    static uint64_t start_pos_dnum[3] = { _INFO_HGVS_del_start_pos, _INFO_HGVS_ins_start_pos, _INFO_HGVS_ins_start_pos }; // ins and delins share start_pos
    rom alt   = tab + 1;
    rom after = &refalt[refalt_len];

    PosType start_pos = ECTX (start_pos_dnum[t])->last_value.i;
    PosType end_pos = (t == DEL) ? (start_pos + tab - refalt - 2)
                    : (t == INS) ? (start_pos + after - alt - 2)
                    : /* DELINS */ (start_pos + after - alt - 1);

    RECONSTRUCT_INT (end_pos);
}

// three separate SPECIAL snips, so that they each because all_the_same
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_DEL_END_POS)    { vcf_piz_special_INFO_HGVS_INDEL_END_POS (vb, DEL);    return NO_NEW_VALUE; }
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_INS_END_POS)    { vcf_piz_special_INFO_HGVS_INDEL_END_POS (vb, INS);    return NO_NEW_VALUE; }
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_DELINS_END_POS) { vcf_piz_special_INFO_HGVS_INDEL_END_POS (vb, DELINS); return NO_NEW_VALUE; }

static void vcf_piz_special_INFO_HGVS_INDEL_PAYLOAD (VBlockP vb, HgvsType t)
{
    STR (refalt);
    reconstruct_peek (vb, CTX (VCF_REFALT), pSTRa(refalt)); // this special works only on SNPs, so length is always 3

    SAFE_NULT (refalt);
    rom tab = strchr (refalt, '\t');
    SAFE_RESTORE;

    ASSPIZ (tab, "invalid REF+ALT=\"%.*s\"", STRf(refalt));

    rom alt   = tab + 1;
    rom after = &refalt[refalt_len];

    rom payload = (t == DEL) ? (refalt + 1) // REF except for the anchor base
                        : (t == INS) ? (alt + 1)    // ALT except for the anchor base
                        : /* DELINS */ alt;         // the entire ALT

    PosType payload_len = (t == DEL) ? (tab - refalt - 1)
                        : (t == INS) ? (after - alt - 1)
                        : /* DELINS */ (after - alt);

    memmove (BAFTtxt, payload, payload_len);
    vb->txt_data.len += payload_len;
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_DEL_PAYLOAD)    { vcf_piz_special_INFO_HGVS_INDEL_PAYLOAD (vb, DEL);    return false; }
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_INS_PAYLOAD)    { vcf_piz_special_INFO_HGVS_INDEL_PAYLOAD (vb, INS);    return false; }
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_DELINS_PAYLOAD) { vcf_piz_special_INFO_HGVS_INDEL_PAYLOAD (vb, DELINS); return false; }

// --------
// INFO/CSQ
// --------

// ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral">
// Originating from the VEP software
// example: CSQ=-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000438504|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|YES||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|6/6||ENST00000450305.2:n.448_449delGC||448-449|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.734_735delGC||734-735|||||rs780379327|1||1||deletion|1|HGNC|37102|YES||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs780379327|1|917|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000515242|transcribed_unprocessed_pseudogene|3/3||ENST00000515242.2:n.727_728delGC||727-728|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000518655|transcribed_unprocessed_pseudogene|3/4||ENST00000518655.2:n.565_566delGC||565-566|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000538476|unprocessed_pseudogene||||||||||rs780379327|1|924|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000541675|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001576075|CTCF_binding_site||||||||||rs780379327|1||||deletion|1|||||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|
static inline void vcf_seg_INFO_CSQ (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    static const MediumContainer csq = {
        .nitems_lo   = 70, 
        .drop_final_repsep = true,
        .repsep      = {','},
        .items       = { { .dict_id={ _INFO_CSQ_Allele             }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_Consequence        }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_IMPACT             }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_SYMBOL             }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_Gene               }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_Feature            }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_Feature            }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_Feature            }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_EXON               }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_INTRON             }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_HGVSc              }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_HGVSp              }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_cDNA_position      }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_CDS_position       }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_Protein_position   }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_Amino_acids        }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_Codons             }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_Existing_variation }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_ALLELE_NUM         }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_DISTANCE           }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_STRAND             }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_FLAGS              }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_VARIANT_CLASS      }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_MINIMISED          }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_SYMBOL_SOURCE      }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_HGNC_ID            }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_CANONICAL          }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_TSL                }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_APPRIS             }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_CCDS               }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_ENSP               }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_SWISSPROT          }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_TREMBL             }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_UNIPARC            }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_GENE_PHENO         }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_SIFT               }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_PolyPhen           }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_DOMAINS            }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_HGVS_OFFSET        }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_AF                 }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_CLIN_SIG           }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_SOMATIC            }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_PHENO              }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_PUBMED             }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_MOTIF_NAME         }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_MOTIF_POS          }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_HIGH_INF_POS       }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_MOTIF_SCORE_CHANGE }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_LoF                }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_LoF_filter         }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_LoF_flags          }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_LoF_info           }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_context            }, .separator = {'|'} },
                         { .dict_id={ _INFO_CSQ_ancestral          }                     } } };

    SegCallback callbacks[70] = { [ 0] =  vcf_seg_INFO_allele,  // INFO_CSQ_Allele
                                  // [10] = vcf_seg_INFO_HGVS,  // INFO_CSQ_HGVSc - commented out bc compresses better without - bc POS with HGVS is relative rather than absolute
                                  [12] = seg_integer_or_not_cb, // INFO_CSQ_cDNA_position - usually "int", sometimes "int-int" where int is an integer or ? 
                                  [13] = seg_integer_or_not_cb, // INFO_CSQ_CDS_position
                                  [14] = seg_integer_or_not_cb, // INFO_CSQ_Protein_position
                                  [17] = seg_id_field_cb,       // INFO_CSQ_Existing_variation
                                  [19] = seg_integer_or_not_cb, // INFO_CSQ_DISTANCE 
                                };
    seg_array_of_struct (VB, ctx, csq, STRa(value), callbacks, value_len);
}

// --------
// INFO/ANN
// --------

// ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
// See: https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
// example: ANN=T|intergenic_region|MODIFIER|U2|ENSG00000277248|intergenic_region|ENSG00000277248|||n.10510103A>T||||||
static inline void vcf_seg_INFO_ANN (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    static const MediumContainer ann = {
        .nitems_lo   = 16, 
        .drop_final_repsep = true,
        .repsep      = {','},
        .items       = { { .dict_id={ _INFO_ANN_Allele             }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_Annotation         }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_Annotation_Impact  }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_Gene_Name          }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_Gene_ID            }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_Feature_Type       }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_Feature_ID         }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_Transcript_BioType }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_Rank               }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_HGVS_c             }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_HGVS_p             }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_cDNA               }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_CDS                }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_AA                 }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_Distance           }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_ANN_Errors             }                     } } };

    seg_array_of_struct (VB, ctx, ann, STRa(value), (SegCallback[]){vcf_seg_INFO_allele,0,0,0,0,0,0,0,0,vcf_seg_INFO_HGVS,0,0,0,0,0,0}, value_len);
}

// ##INFO=<ID=RS,Number=1,Type=Integer,Description="dbSNP ID (i.e. rs number)">
// *might* be the same as the numeric value of the ID
static inline void vcf_seg_INFO_RS (VBlockVCFP vb, ContextP ctx, STRp(rs))
{
    // case: eg ID=rs3844233 RS=3844233. We use a SNIP_DELTA_OTHER with delta=0 to copy last_value from ID.
    // (can't use SNIP_COPY bc it would copy the entire txt "rs3844233")
    int64_t rs_value;
    if (ctx_has_value_in_line_(VB, CTX(VCF_ID)) &&
        str_get_int (STRa(rs), &rs_value) &&
        rs_value == CTX(VCF_ID)->last_value.i) {

        seg_by_ctx (VB, STRa(copy_ID_snip), ctx, rs_len);
    }
    
    else 
        seg_integer_or_not (VB, ctx, STRa(rs), rs_len);
}

// ##INFO=<ID=RSPOS,Number=1,Type=Integer,Description="Chr position reported in dbSNP">
// *might* be the same as POS
static inline void vcf_seg_INFO_RSPOS (VBlockVCFP vb, ContextP ctx, STRp(rspos))
{
    // case: eg ID=rs3844233 RS=3844233. We use a SNIP_DELTA_OTHER with delta=0 to copy last_value from ID.
    // (can't use SNIP_COPY bc it would copy the entire txt "rs3844233")
    int64_t rspos_value;
    if (ctx_has_value_in_line_(VB, CTX(VCF_POS)) &&
        str_get_int (STRa(rspos), &rspos_value)) {

        if (rspos_value == CTX(VCF_POS)->last_value.i) // shortcut for most common case
            seg_by_ctx (VB, STRa(copy_POS_snip), ctx, rspos_len);

        else {
            STRl(snip, 32) = 32;
            seg_prepare_snip_other (SNIP_OTHER_DELTA, _VCF_POS, true, rspos_value - CTX(VCF_POS)->last_value.i, snip);
            seg_by_ctx (VB, STRa(snip), ctx, rspos_len);            
        }
    }
    
    else 
        seg_integer_or_not (VB, ctx, STRa(rspos), rspos_len);
}

// ##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
// comma-seperated two numbers: RAW_MQandDP=720000,200: 1. sum of squared MQ values and 2. total reads over variant genotypes (note: INFO/MQ is sqrt(#1/#2))
static inline void vcf_seg_INFO_RAW_MQandDP (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    seg_by_container (VB, ctx, (ContainerP)&RAW_MQandDP_con, STRa(value), STRa(RAW_MQandDP_snip), true, value_len);
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

            if (vb->line_coords == DC_PRIMARY) vb->recon_size += growth;
            else vb->recon_size_luft += growth;
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
    #define ADJUST_FOR_MODIFIED ({ \
        int32_t shrinkage = (int32_t)value_len - (int32_t)modified_len;\
        vb->recon_size -= shrinkage; \
        vb->recon_size_luft -= shrinkage; \
        value = modified; \
        value_len = modified_len; \
    })

    ctx->line_is_luft_trans = false; // initialize
    
    // --chain: if this is RendAlg=A_1 subfield in a REF<>ALT variant, convert a eg 4.31e-03 to e.g. 0.00431. This is to
    // ensure primary->luft->primary is lossless (4.31e-03 cannot be converted losslessly as we can't preserve format info)
    if (chain_is_loaded && ctx->luft_trans == VCF2VCF_A_1 && LO_IS_OK_SWITCH (last_ostatus) && 
        str_scientific_to_decimal (STRa(value), modified, &modified_len, NULL))
        ADJUST_FOR_MODIFIED;
        
    // Translatable item on a Luft line: attempt to lift-back the value, so we can seg it as primary
    if (vb->line_coords == DC_LUFT && needs_translation (ctx)) {

        // If cross-rendering to Primary is successful - proceed to Seg this value in primary coords, and assign a translator for reconstructing if --luft
        if (vcf_lo_seg_cross_render_to_primary (vb, ctx, STRa(value), modified, &modified_len)) {
            value = modified; 
            value_len = modified_len; 
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
    if (z_is_dvcf && ctx->luft_trans == VCF2VCF_ALLELE) 
        vcf_seg_INFO_allele (VB, ctx, STRa(value), 0);

    else switch (ctx->dict_id.num) {
        #define CALL(f) (f); break
        #define CALL_IF(cond,f) if (cond) { (f); break; } else goto standard_seg 
        #define CALL_WITH_FALLBACK(f) if (f(vb, ctx, STRa(value))) { seg_by_ctx (VB, STRa(value), ctx, value_len); } break
        #define STORE_AND_SEG(store_type) seg_set_last_txt_store_value (VB, ctx, STRa(value), store_type); seg_by_ctx (VB, STRa(value), ctx, value_len); break

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

        // if SVLEN is negative, it is expected to be minus the delta between END and POS
        case _INFO_SVLEN:
            CALL_IF (vcf_seg_test_SVLEN (vb, STRa(value)), 
                     seg_by_ctx (VB, ((char [2]){ SNIP_SPECIAL, VCF_SPECIAL_SVLEN }), 2, ctx, value_len));

        case _INFO_SVTYPE:
            CALL_WITH_FALLBACK (vcf_seg_SVTYPE);

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
        case _INFO_AA: 
            CALL (vcf_seg_INFO_allele (VB, ctx, STRa(value), 0));

        // ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral">
        // Originating from the VEP software
        // example: CSQ=-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000438504|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|YES||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|6/6||ENST00000450305.2:n.448_449delGC||448-449|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.734_735delGC||734-735|||||rs780379327|1||1||deletion|1|HGNC|37102|YES||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs780379327|1|917|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000515242|transcribed_unprocessed_pseudogene|3/3||ENST00000515242.2:n.727_728delGC||727-728|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000518655|transcribed_unprocessed_pseudogene|3/4||ENST00000518655.2:n.565_566delGC||565-566|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000538476|unprocessed_pseudogene||||||||||rs780379327|1|924|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000541675|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001576075|CTCF_binding_site||||||||||rs780379327|1||||deletion|1|||||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|
        case _INFO_CSQ:
            CALL (vcf_seg_INFO_CSQ (vb, ctx, STRa(value)));
        
        // ##INFO=<ID=DP_HIST,Number=R,Type=String,Description="Histogram for DP; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">
        // ##INFO=<ID=GQ_HIST,Number=R,Type=String,Description="Histogram for GQ; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">
        // ##INFO=<ID=AGE_HISTOGRAM_HET,Number=A,Type=String,Description="Histogram of ages of allele carriers; Bins: <30|30|35|40|45|50|55|60|65|70|75|80+">
        // ##INFO=<ID=AGE_HISTOGRAM_HOM,Number=A,Type=String,Description="Histogram of ages of homozygous allele carriers; Bins: <30|30|35|40|45|50|55|60|65|70|75|80+">
        case _INFO_vep:
        case _INFO_DP_HIST:
        case _INFO_GQ_HIST:
        case _INFO_AGE_HISTOGRAM_HET:
        case _INFO_AGE_HISTOGRAM_HOM: 
            CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), ',', '|', false, STORE_INT, DICT_ID_NONE, value_len));

        case _INFO_DP4:
            CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_INT, DICT_ID_NONE, value_len));

        // ##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
        case _INFO_CLNDN:
            CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), '|', 0, false, STORE_NONE, DICT_ID_NONE, value_len));

        // ##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
        case _INFO_CLNHGVS:
            CALL (vcf_seg_INFO_HGVS (VB, ctx, STRa(value), 0));

        // ##INFO=<ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
        // example: CPIC:0b3ac4db1d8e6e08a87b6942|CPIC:647d4339d5c1ddb78daff52f|CPIC:9968ce1c4d35811e7175cd29|CPIC:PA166160951|CPIC:c6c73562e2b9e4ebceb0b8bc
        // I tried seg_array_of_struct - it is worse than simple seg

        // ##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
        case _INFO_ALLELEID:
            CALL (seg_integer_or_not (VB, ctx, STRa(value), value_len));

        // ##INFO=<ID=dbSNPBuildID,Number=1,Type=Integer,Description="First dbSNP Build for RS">
        case _INFO_dbSNPBuildID:
            CALL (seg_integer_or_not (VB, ctx, STRa(value), value_len));

        // ##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
        case _INFO_RS:
            CALL (vcf_seg_INFO_RS (vb, ctx, STRa(value)));

        // ##INFO=<ID=RSPOS,Number=1,Type=Integer,Description="Chr position reported in dbSNP">
        case _INFO_RSPOS:
            CALL (vcf_seg_INFO_RSPOS (vb, ctx, STRa(value)));

        // ##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Pairs each of gene symbol:gene id.  The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
        case _INFO_GENEINFO:
            CALL (seg_array (VB, ctx, INFO_GENEINFO, STRa(value), '|', 0, false, STORE_NONE, DICT_ID_NONE, value_len));

        // case _INFO_TOPMED: // better leave as simple snip as the items are allele frequencies which are correleted

        // ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
        case _INFO_ANN: 
            CALL (vcf_seg_INFO_ANN (vb, ctx, STRa(value)));

        // ##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
        case _INFO_RAW_MQandDP:
            CALL (vcf_seg_INFO_RAW_MQandDP (vb, ctx, STRa(value)));

        case _INFO_PROBE_A:         CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_PROBE_A      (vb, ctx, STRa(value)));
        case _INFO_PROBE_B:         CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_PROBE_B      (vb, ctx, STRa(value)));
        case _INFO_ALLELE_A:        CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ALLELE_A     (vb, ctx, STRa(value)));
        case _INFO_ALLELE_B:        CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ALLELE_B     (vb, ctx, STRa(value)));
        case _INFO_ILLUMINA_CHR:    CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ILLUMINA_CHR (vb, ctx, STRa(value)));
        case _INFO_ILLUMINA_POS:    CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ILLUMINA_POS (vb, ctx, STRa(value)));
        case _INFO_ILLUMINA_STRAND: CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ILLUMINA_STRAND (vb, ctx, STRa(value)));
        case _INFO_refSNP:          CALL_IF (segconf.vcf_illum_gtyping, seg_id_field_do      (VB, ctx, STRa(value)));

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

void vcf_seg_info_subfields (VBlockVCFP vb, rom info_str, unsigned info_len)
{
    vb->info_items.len = 0; // reset from previous line

    // case: INFO field is '.' (empty) (but not in DVCF as we will need to deal with DVCF items)
    if (!z_is_dvcf && info_len == 1 && *info_str == '.') {
        seg_by_did (VB, ".", 1, VCF_INFO, 2); // + 1 for \t or \n
        return;
    }

    // parse the info string
    str_split (info_str, info_len, MAX_FIELDS-2, ';', pair, false); // -2 - leave room for LUFT + PRIM
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

        if (segconf.running) segconf.has[ii.ctx->did_i] = true;

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
    if (vcf_is_use_DP_by_DP() && ctx_has_value_in_line_(VB, CTX(INFO_DP))) 
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

    // case GVCF: multiplex by has_RGQ
    if (!segconf.running && segconf.has[FORMAT_RGQ]) {
        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, VCF_INFO, (MultiplexerP)&vb->mux_INFO, vb->line_has_RGQ);
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

