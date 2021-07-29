// ------------------------------------------------------------------
//   vcf_info.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

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
#include "gff3.h"
#include "coords.h"

typedef struct { char name[64];
                 const char *value; 
                 unsigned name_len, value_len; 
                 DictId dict_id;   } InfoItem;

//--------
// INFO/DP
// -------

// return true if caller still needs to seg 
static bool vcf_seg_INFO_DP (VBlockVCF *vb, ContextP ctx_dp, const char *value, int value_len)
{
    // also tried delta vs DP4, but it made it worse
    Context *ctx_basecounts;
    if (ctx_has_value_in_line (vb, dict_id_INFO_BaseCounts, &ctx_basecounts)) {
        seg_delta_vs_other (vb, ctx_dp, ctx_basecounts, value, value_len, -1);
        return false; // caller needn't seg
    }
    else {
        // store last_value of INFO/DP field in case we have FORMAT/DP as well (used in vcf_seg_one_sample)
        ctx_set_last_value (vb, ctx_dp, (int64_t)atoi (value));
        return true; // caller should seg
    }
}

//--------
// INFO/SF
// -------

#define adjustment vb->sf_ctx->last_delta
#define next param

// INFO/SF contains a comma-seperated list of the 0-based index of the samples that are NOT '.'
// example: "SF=0,1,15,22,40,51,59,78,88,89,90,112,124,140,147,155,156,164,168,183,189,197,211,215,216,217,222,239,244,256,269,270,277,281,290,291,299,323,338,340,348" 
// Algorithm: SF is segged either as an as-is string, or as a SPECIAL that includes the index of all the non-'.' samples. 
// if use_special_sf=YES, we use SNIP_SPECIAL and we validate the correctness during vcf_seg_FORMAT_GT -
// if it is wrong we set use_special_sf=NO. The assumption is that normally, it is either true for all lines or false.
static bool vcf_seg_INFO_SF_init (VBlockVCF *vb, Context *sf_ctx, const char *value, int value_len)
{
    switch (vb->use_special_sf) {

        case USE_SF_NO: 
            return true; // "special" is suppressed - caller should go ahead and seg normally

        case USE_SF_UNKNOWN: 
            vb->use_special_sf = USE_SF_YES; // first call to this function, after finding that we have an INFO/SF field - set and fall through
            vb->sf_ctx = sf_ctx;

        case USE_SF_YES: 
            // we store the SF value in a buffer, since seg_FORMAT_GT overlays the haplotype buffer onto txt_data and may override the SF field
            // we will need the SF data if the field fails verification in vcf_seg_INFO_SF_one_sample
            buf_alloc (vb, &vb->sf_txt, 0, value_len + 1, char, 2, "sf_txt"); // +1 for nul-terminator
            memcpy (FIRSTENT (char, vb->sf_txt), value, value_len);
            vb->sf_txt.len = value_len;
            *AFTERENT (char, vb->sf_txt) = 0; // nul-terminate
            vb->sf_txt.next = 0; 
            adjustment = 0;      
            
            // snip being contructed 
            buf_alloc (vb, &vb->sf_snip, 0, value_len + 20, char, 2, "sf_snip"); // initial value - we will increase if needed
            NEXTENT (char, vb->sf_snip) = SNIP_SPECIAL;
            NEXTENT (char, vb->sf_snip) = VCF_SPECIAL_SF;

            return false; // caller should not seg as we already did

        default:
            ABORT_R ("Error in vcf_seg_INFO_SF_init: invalid use_special_sf=%d", vb->use_special_sf);
    }
}

// verify next number on the list of the SF field is sample_i (called from vcf_seg_FORMAT_GT)
void vcf_seg_INFO_SF_one_sample (VBlockVCF *vb)
{
    // case: no more SF values left to compare - we ignore this sample
    while (vb->sf_txt.next < vb->sf_txt.len) {

        buf_alloc (vb, &vb->sf_snip, 10, 0, char, 2, "sf_snip");

        char *sf_one_value = ENT (char, vb->sf_txt, vb->sf_txt.next); 
        char *after;
        int32_t value = strtol (sf_one_value, &after, 10);

        int32_t adjusted_sample_i = (int32_t)(vb->sample_i + adjustment); // adjustment is the number of values in SF that are not in samples

        // case: badly formatted SF field
        if (*after != ',' && *after != 0) {
            vb->use_special_sf = USE_SF_NO; // failed - turn off for the remainder of this vb
            break;
        }
        
        // case: value exists in SF and samples
        else if (value == adjusted_sample_i) {
            NEXTENT (char, vb->sf_snip) = ',';
            vb->sf_txt.next = ENTNUM (vb->sf_txt, after) + 1; // +1 to skip comma
            break;
        }

        // case: value in SF file doesn't appear in samples - keep the value in the snip
        else if (value < adjusted_sample_i) {
            vb->sf_snip.len += str_int (value, AFTERENT (char, vb->sf_snip));
            NEXTENT (char, vb->sf_snip) = ',';
            adjustment++;
            vb->sf_txt.next = ENTNUM (vb->sf_txt, after) + 1; // +1 to skip comma
            // continue and read the next value
        }

        // case: value in SF is larger than current sample - don't advance iterator - perhaps future sample will cover it
        else { // value > adjusted_sample_i
            NEXTENT (char, vb->sf_snip) = '~'; // skipped sample
            break; 
        }
    }
}

static void vcf_seg_INFO_SF_seg (VBlockVCF *vb)
{   
    // case: SF data remains after all samples - copy it
    int32_t remaining_len = (uint32_t)(vb->sf_txt.len - vb->sf_txt.next); // -1 if all done, because we skipped a non-existing comma
    if (remaining_len > 0) {
        buf_add_more (vb, &vb->sf_snip, ENT (char, vb->sf_txt, vb->sf_txt.next), remaining_len, "sf_snip");
        NEXTENT (char, vb->sf_snip) = ','; // buf_add_more allocates one character extra
    }

    if (vb->use_special_sf == USE_SF_YES) 
        seg_by_ctx (vb, vb->sf_snip.data, vb->sf_snip.len, vb->sf_ctx, vb->sf_txt.len);
    
    else if (vb->use_special_sf == USE_SF_NO)
        seg_by_ctx (vb, vb->sf_txt.data, vb->sf_txt.len, vb->sf_ctx, vb->sf_txt.len);

    buf_free (&vb->sf_txt);
    buf_free (&vb->sf_snip);
}

#undef adjustment
#undef param

#define adjustment vcf_vb->sf_ctx->last_delta
#define sample_i   vcf_vb->sf_ctx->last_value.i
#define snip_i     vcf_vb->sf_snip.param

// leave space for reconstructing SF - actual reconstruction will be in vcf_piz_container_cb
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_SF)
{
    VBlockVCFP vcf_vb = (VBlockVCFP)vb;

    if (reconstruct) {
        vcf_vb->sf_ctx = ctx;
        adjustment    = 0;
        sample_i      = 0; 
        snip_i        = 0;

        // temporary place for SF
        buf_alloc_old (vb, &vcf_vb->sf_txt, 5 * vcf_header_get_num_samples(), 1, "sf_txt" ); // initial estimate, we may further grow it later
        vcf_vb->sf_txt.len = 0;

        // copy snip to sf_snip (note: the SNIP_SPECIAL+code are already removed)
        buf_alloc_old (vb, &vcf_vb->sf_snip, snip_len, 2, "sf_snip");
        vcf_vb->sf_snip.len = snip_len; 
        memcpy (vcf_vb->sf_snip.data, snip, snip_len); 
    }

    return false; // no new value
}

// While reconstructing the GT fields of the samples - calculate the INFO/SF field
void vcf_piz_GT_cb_calc_INFO_SF (VBlockVCFP vcf_vb, unsigned rep, char *recon, int32_t recon_len)
{
    if (rep != 0) return; // we only look at the first ht in a sample, and only if its not '.'/'%'

    if (*recon == '.' || *recon == '%') { // . can be written as % in vcf_seg_FORMAT_GT
        sample_i++;
        return;
    }

    ARRAY (const char, sf_snip, vcf_vb->sf_snip);

    while (snip_i < sf_snip_len) {

        buf_alloc (vcf_vb, &vcf_vb->sf_txt, 12, 0, char, 2, "sf_txt"); // sufficient for int32 + ','

        int32_t adjusted_sample_i = (int32_t)(sample_i + adjustment); // last_delta is the number of values in SF that are not in samples

        // case: we derive the value from sample_i
        if (sf_snip[snip_i] == ',') {
            snip_i++;

            buf_add_int ((VBlockP)vcf_vb, &vcf_vb->sf_txt, adjusted_sample_i);

            if (snip_i < sf_snip_len) // add comma if not done yet
                NEXTENT (char, vcf_vb->sf_txt) = ',';

            break;
        }

        // case: sample was not used to derive any value
        if (sf_snip[snip_i] == '~') {
            snip_i++;
            break;
        }

        // case: value quoted verbatim - output and continue searching for , or ~ for this sample
        else {
            char *comma = strchr (&sf_snip[snip_i], ',');
            unsigned len = comma - &sf_snip[snip_i];
            bool add_comma = (snip_i + len + 1 < sf_snip_len); // add comma if not last

            buf_add (&vcf_vb->sf_txt, &sf_snip[snip_i], len + add_comma);
            snip_i += len + 1;

            adjustment++;
        }
    }

    sample_i++;
}

// Upon completing the line - insert the calculated INFO/SF field to its place
void vcf_piz_TOPLEVEL_cb_insert_INFO_SF (VBlockVCFP vcf_vb)
{
    ARRAY (const char, sf_snip, vcf_vb->sf_snip);

    // if there are some items remaining in the snip (values that don't appear in samples) - copy them
    if (snip_i < sf_snip_len) {
        unsigned remaining_len = sf_snip_len - snip_i - 1; // all except the final comma
        buf_add_more (vcf_vb, &vcf_vb->sf_txt, &sf_snip[snip_i], remaining_len, "sf_txt"); 
    }

    // make room for the SF txt and copy it to its final location
    char *sf_txt = last_txt (vcf_vb, vcf_vb->sf_ctx->did_i);
    memmove (sf_txt + vcf_vb->sf_txt.len, sf_txt, AFTERENT (char, vcf_vb->txt_data) - sf_txt); // make room
    memcpy (sf_txt, vcf_vb->sf_txt.data, vcf_vb->sf_txt.len); // copy

    vcf_vb->txt_data.len += vcf_vb->sf_txt.len;

    buf_free (&vcf_vb->sf_snip);
    buf_free (&vcf_vb->sf_txt);
}

#undef sample_i
#undef adjustment
#undef snip_i

// -------
// INFO/AA
// -------

// checks if value is identifcal to the REF or one of the ALT alleles, and if so segs a SPECIAL snip
static bool vcf_seg_INFO_allele (VBlockVCF *vb, Context *ctx, const char *value, int value_len) // returns true if caller still needs to seg 
{
    // short circuit in common case of '.'
    if (value_len == 1 && *value == '.') return true; // caller should seg
    
    int allele = -1;
    
    // check if its equal main REF (which can by REF or oREF)
    if (value_len == vb->main_ref_len && !memcmp (value, vb->main_refalt, value_len))
        allele = 0;

    // check if its equal one of the ALTs
    else {
        str_split (&vb->main_refalt[vb->main_ref_len+1], vb->main_alt_len, 0, ',', alt, false);

        for (int alt_i=0; alt_i < n_alts; alt_i++) 
            if (value_len == alt_lens[alt_i] && !memcmp (value, alts[alt_i], value_len)) {
                allele = alt_i + 1;
                break;
            }
    }

    if (allele == -1) return true; // caller should seg

    char snip[] = { SNIP_SPECIAL, VCF_SPECIAL_ALLELE, '0' + vb->line_coords, '0' + allele /* ASCII 48...147 */ };
    seg_by_ctx (vb, snip, sizeof (snip), ctx, value_len);
    
    return false;
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_ALLELE)
{
    Coords seg_line_coord = snip[0] - '0';
    int allele = snip[1] - '0';
    LiftOverStatus ostatus = last_ostatus;

    if (LO_IS_OK_SWITCH (ostatus) && seg_line_coord != vb->vb_coords) {
        ASSPIZ (allele >= 0 && allele <= 1, "unexpected allele=%d with REF<>ALT switch", allele);
        allele = 1 - allele;
    }

    DidIType refalt_did = (vb->vb_coords == DC_PRIMARY ? VCF_REFALT : VCF_oREFALT); // we need the one that was already reconstructed in REF/ALT
    const char *refalt = last_txt (vb, refalt_did);
    unsigned refalt_len = vb->last_txt_len (refalt_did);

    if (!refalt_len) goto done; // variant is single coordinate in the other coordinate

    char *tab = memchr (refalt, '\t', refalt_len);
    ASSPIZ (tab, "Invalid refalt: \"%.*s\"", MIN (refalt_len, 100), refalt);

    // case: the allele is REF
    if (allele == 0)
        RECONSTRUCT (refalt, tab - refalt);
    
    // case: the allele is one of the alts
    else {
        str_split (tab+1, &refalt[refalt_len] - (tab+1), 0, ',', alt, false);
        RECONSTRUCT (alts[allele-1], alt_lens[allele-1]);
    }

done:
    return false; // no new value
}

// translator only validates - as vcf_piz_special_INFO_ALLELE copies verbatim
TRANSLATOR_FUNC (vcf_piz_luft_ALLELE)
{
    VBlockVCFP vcf_vb = ((VBlockVCFP)vb);

    // reject if LO_OK_REF_NEW_SNP and value is equal to REF
    if (validate_only && last_ostatus == LO_OK_REF_NEW_SNP && 
        recon_len == vcf_vb->main_ref_len && !memcmp (recon, vcf_vb->main_refalt, recon_len)) return false;

    return true;
}

// ---------------
// INFO/BaseCounts
// ---------------

// ##INFO=<ID=genozip BugP.vcf -ft ,Number=4,Type=Integer,Description="Counts of each base">
// Sorts BaseCounts vector with REF bases first followed by ALT bases, as they are expected to have the highest values
static bool vcf_seg_INFO_BaseCounts (VBlockVCF *vb, Context *ctx_basecounts, const char *value, int value_len) // returns true if caller still needs to seg 
{
    if (CTX(VCF_REFALT)->last_snip_len != 3 || vb->line_coords == DC_LUFT) 
        return true; // not a bi-allelic SNP or line is a luft line without easy access to REFALT - caller should seg

    char *str = (char *)value;
    int64_t sum = 0;

    uint32_t counts[4], sorted_counts[4] = {}; // corresponds to A, C, G, T

    SAFE_NUL (&value[value_len]);
    for (unsigned i=0; i < 4; i++) {
        counts[i] = strtoul (str, &str, 10);
        str++; // skip comma seperator
        sum += counts[i];
    }
    SAFE_RESTORE;

    if (str - value != value_len + 1 /* +1 due to final str++ */) return true; // invalid BaseCounts data - caller should seg

    const char *refalt = CTX(VCF_REFALT)->last_snip;

    unsigned ref_i = acgt_encode[(int)refalt[0]];
    unsigned alt_i = acgt_encode[(int)refalt[2]];

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

    seg_by_ctx (vb, snip, value_len+2, ctx_basecounts, value_len); 
    
    ctx_basecounts->flags.store = STORE_INT;
    ctx_set_last_value (vb, ctx_basecounts, sum);

    return false; // we already segged - caller needn't seg
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_BaseCounts)
{
    const char *refalt; unsigned refalt_len;
    reconstruct_peek (vb, CTX(VCF_REFALT), &refalt, &refalt_len);

    uint32_t counts[4], sorted_counts[4] = {}; // counts of A, C, G, T

    new_value->i = 0;
    char *str = (char *)snip;

    for (unsigned i=0; i < 4; i++) {
        sorted_counts[i] = strtoul (str, &str, 10);
        str++; // skip comma seperator
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
    return true; // has new value
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

// --------
// INFO/CSQ
// --------

// INFO fields with a format originating from the VEP software, eg
// vep=T|intergenic_variant|MODIFIER|||Intergenic||||||||||||1|||SNV||||||||||||||||||||||||
static inline void vcf_seg_INFO_CSQ (VBlockVCF *vb, Context *vep_ctx, const char *field, unsigned field_len)
{
    Container con = { .repsep = { ',' }, 
                      .drop_final_repeat_sep = true,
                      .keep_empty_item_sep   = true }; // don't delete the | before an empty item

    Context *sf_ctxs[MAX_FIELDS] = {};

    uint32_t item_i=0;
    const char *item_start = field;
    for (uint32_t i=0; i < field_len+1; i++) {
        if (i == field_len || field[i] == ',' || field[i] == '|') { // end of item
            if (item_i == con_nitems(con)) {
                ASSVCF (!con.repeats, 
                        "expecting all repeats of %s to have the same number of items, %u, as the first repeat, but repeat %u (0-based) has more: %.*s", 
                        vep_ctx->name, con_nitems(con), con.repeats, field_len, field);

                ASSVCF (item_i < MIN (126, MAX_FIELDS), "exceeded the max number of %s items=%u", 
                        vep_ctx->name, MIN (126, MAX_FIELDS)); // the 126 constraint is just the context naming scheme
                
                char name[8];
                sprintf (name, "%c%c_%.3s", item_i < 63 ? '_' : '`', '@' + (item_i % 63), vep_ctx->name);
                DictId dict_id = dict_id_make (name, 6, DTYPE_VCF_INFO);

                sf_ctxs[item_i] = ctx_get_ctx (vb, dict_id);
                sf_ctxs[item_i]->st_did_i = vep_ctx->did_i;

                con.items[item_i] = (ContainerItem){ .dict_id   = dict_id, 
                                                     .seperator = { field[i]=='|' ? '|' : 0 } };
                con_inc_nitems (con);                                     
            }

            unsigned item_len = &field[i] - item_start;
            seg_by_ctx (vb, item_start, item_len, sf_ctxs[item_i], item_len + (i != field_len));

            item_i++;
            item_start = &field[i+1];

            if (field[i] != '|') { // end of repeat
                ASSVCF (!con.repeats || item_i == con_nitems(con), 
                        "expecting all repeats of %s to have the same number of items, %u, as the first repeat, but repeat %u (0-based) has only %u items: %.*s", 
                        vep_ctx->name, con_nitems(con), con.repeats, item_i, field_len, field);
            
                con.repeats++;
                item_i=0;
            }
        }
    }

    container_seg_by_ctx (vb, vep_ctx, &con, 0, 0, 0);

    // TODO: recover and seg just the field as-is, rather than throw an error, if its not the CSQ format we expect
    // (need to roll back all changes to subfields)
}

// -------
// INFO/AC
// -------

static void vcf_seg_INFO_AC (VBlockVCF *vb, Context *ac_ctx, const char *field, unsigned field_len)
{
    Context *af_ctx, *an_ctx;
    int64_t ac;
    bool ac_has_value_in_line = str_get_int (field, field_len, &ac);
    bool an_has_value_in_line = ctx_has_value_in_line (vb, dict_id_INFO_AN, &an_ctx);

    // case: AC = AN * AF (might not be, due to rounding errors, esp if AF is a very small fraction)
    if (ac_has_value_in_line &&                                 // AC exists and is a valid int
        an_has_value_in_line &&                                 // AN exists and is a valid int
        ctx_has_value_in_line (vb, dict_id_INFO_AF, &af_ctx) && // AF exists and is a valid float
        (int64_t)round (af_ctx->last_value.f * an_ctx->last_value.i) == ac) { // AF * AN == AC

        ac_ctx->no_stons = true;
        seg_by_ctx (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_AC }), 2, ac_ctx, field_len);
    }

    // case: AC is multi allelic, or an invalid value, or missing AN or AF or AF*AN != AC
    else 
        seg_by_ctx (vb, field, field_len, ac_ctx, field_len);

    ac_ctx->flags.store = STORE_INT;
}

// reconstruct: AC = AN * AF
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_AC)
{
    if (!reconstruct) return false;

    // Backward compatability note: In files v6->11, snip has 2 bytes for AN, AF which mean: '0'=appears after AC, '1'=appears before AC. We ignore them.

    // note: update last_value too, so its available to vcf_piz_luft_A_AN, which is called becore last_value is updated
    ctx->last_value.i = new_value->i = (int64_t)round (reconstruct_peek_(vb, dict_id_INFO_AN, 0, 0).i * reconstruct_peek_(vb, dict_id_INFO_AF, 0, 0).f);
    RECONSTRUCT_INT (new_value->i); 

    return true;
}

// Lift-over translator for INFO/AC fields, IF it is bi-allelic and we have a ALT<>REF switch AND we have an AN and AF field
// We change the AC to AN - AC and return true if successful 
TRANSLATOR_FUNC (vcf_piz_luft_A_AN)
{
    if (IS_TRIVAL_FORMAT_SUBFIELD) return true; // This is FORMAT field which is empty or "." - all good

    int64_t an, ac;

    if (command == ZIP) {
        Context *an_ctx;
        if (!ctx_has_value_in_line (vb, dict_id_INFO_AN, &an_ctx) ||
            !str_get_int_range64 (recon, recon_len, 0, (an = an_ctx->last_value.i), &ac))
            return false; // in Seg, AC is always segged last, so this means there is no AN in the line for sure
        
        if (validate_only) return true;  // Yay! AC can be lifted - all it needs in AN in the line, which it has 
    }
    else {
        ac = ctx->last_value.i;
        an = reconstruct_peek_(vb, dict_id_INFO_AN, 0, 0).i;
    }

    // re-reconstruct: AN-AC
    vb->txt_data.len -= recon_len;
    RECONSTRUCT_INT (an - ac);

    return true;    
}

// --------
// INFO/END
// --------

static void vcf_seg_INFO_END (VBlockVCFP vb, Context *end_ctx, const char *end_str, unsigned end_len) // note: ctx is INFO/END *not* POS (despite being an alias)
{
    // END is an alias of POS
    seg_pos_field ((VBlockP)vb, VCF_POS, VCF_POS, SPF_BAD_SNIPS_TOO | SPF_ZERO_IS_BAD | SPF_UNLIMIED_DELTA, 0, end_str, end_len, 0, end_len);

    // add END to dl for sorting. it is used only in case chrom and pos are identical
    PosType end = vb->last_int (VCF_POS); 
    DATA_LINE (vb->line_i)->end = end;

    // case --chain: if we have lifted-over POS (as primary POS field or in INFO/LIFTBACK), 
    // check that lifting-over of END is delta-encoded and is lifted over to the same, non-xstrand, Chain alignment, and reject if not
    if (chain_is_loaded && LO_IS_OK (last_ostatus)) { 

        bool is_xstrand = vb->last_index (VCF_oXSTRAND); // set in vcf_lo_seg_generate_INFO_DVCF
        PosType aln_last_pos = chain_get_aln_last_pos (vb->pos_aln_i); 

        // case: we don't yet handle END translation in case of a reverse strand
        if (is_xstrand)            
            REJECT_SUBFIELD (LO_INFO, end_ctx, "Genozip limitation: variant with INFO/END and chain file alignment with a negative strand%s", "");

        // case: END goes beyond the end of the chain file alignment and its a <DEL>
        else if (vb->is_del_sv && end > aln_last_pos) {

            // case: END goes beyond end of alignment
            PosType gap_after = chain_get_aln_gap_after (vb->pos_aln_i);
            
            // case: END falls in the gap after - <DEL> is still valid but translated END needs to be closer to POS to avoid gap - 
            // we don't yet do this
            if (end <= aln_last_pos + gap_after)
                REJECT_SUBFIELD (LO_INFO, end_ctx, "Genozip limitation: <DEL> variant: INFO/END=%.*s is in the gap after the end of the chain file alignment", end_len, end_str);
    
            // case: END falls on beyond the gap (next alignment or beyond) - this variant cannot be lifted
            else
                REJECT_SUBFIELD (LO_INFO, end_ctx, "<DEL> variant: INFO/END=%.*s is beyond the end of the chain file alignment and also beyond the gap after the alignment", end_len, end_str);
        }

        // case: END goes beyond the end of the chain file alignment and its NOT a <DEL>
        else if (!vb->is_del_sv && end > aln_last_pos) 
            REJECT_SUBFIELD (LO_INFO, end_ctx, "POS and INFO/END=%.*s are not on the same chain file alignment", end_len, end_str);

        // case: invalid value. since we use SPF_UNLIMIED_DELTA, any integer value should succeed
        else if (!CTX(VCF_POS)->last_delta)
            REJECT_SUBFIELD (LO_INFO, end_ctx, "INFO/END=%.*s has an invalid value", end_len, end_str);        
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

    // ZIP liftback (POS is always before END, because we seg INFO/LIFTOVER first)
    if (command == ZIP && vb->line_coords == DC_LUFT) { // liftback
        PosType oend;
        if (!str_get_int_range64 (recon, recon_len, 0, MAX_POS, &oend))
            return false;

        translated_end = vb->last_int (VCF_POS) + (oend - vb->last_int (VCF_oPOS));
    }

    // PIZ liftover: we have already reconstructed oPOS and POS (as an item in VCF_TOPLUFT)
    else {
        translated_end = vb->last_int (VCF_oPOS) + vb->last_delta (VCF_POS); // delta for generated this END value (END - POS)

        ((VBlockVCFP)vb)->last_end_line_i = vb->line_i; // so vcf_piz_special_COPYPOS knows that END was reconstructed
    }

    // re-reconstruct END
    vb->txt_data.len -= recon_len; 
    RECONSTRUCT_INT (translated_end);
    
    return true;    
}

// Called to reconstruct the POS subfield of INFO/LIFTBACK, handling the possibility of INFO/END
SPECIAL_RECONSTRUCTOR (vcf_piz_special_COPYPOS)
{
    if (!reconstruct) return false; // no new value
    
    bool has_end = ((VBlockVCFP)vb)->last_end_line_i == vb->line_i; // true if INFO/END was encountered

    Context *pos_ctx = CTX(VCF_POS);
    int64_t pos;

    if (has_end) {
        int64_t end   = pos_ctx->last_value.i;
        int64_t delta = pos_ctx->last_delta;
        pos = end - delta;
    }
    else
        pos = pos_ctx->last_value.i;

    RECONSTRUCT_INT (pos); // vcf_piz_luft_END makes sure it always contains the value of POS, not END
    return false; // no new value
}

// ----------
// INFO/SVLEN
// ----------

static inline bool vcf_seg_test_SVLEN (VBlockVCF *vb, const char *svlen_str, unsigned svlen_str_len)
{
    int64_t svlen;
    if (!str_get_int (svlen_str, svlen_str_len, &svlen)) return false;

    int64_t last_delta = CTX(VCF_POS)->last_delta; // INFO_END is an alias of POS - so the last delta would be between END and POS
    return last_delta == -svlen;
}

// the case where SVLEN is minus the delta between END and POS
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_SVLEN)
{
    if (!reconstruct) goto done;

    int64_t value = -CTX(VCF_POS)->last_delta; // END is a alias of POS - they share the same data stream - so last_delta would be the delta between END and POS
    char str[30];
    unsigned str_len = str_int (value, str);
    RECONSTRUCT (str, str_len);

done:
    return false; // no new value
}

// --------------
// INFO container
// --------------

// for dual coordinate files (Primary, Luft and --chain) - add DVCF depending on ostatus (run after
// all INFO and FORMAT fields, so ostatus is final)
static void vcf_seg_info_add_DVCF_to_InfoItems (VBlockVCF *vb)
{
    // case: Dual coordinates file line has no PRIM, Lrej or Prej - this can happen if variants were added to the file,
    // for example, as a result of a "bcftools merge" with a non-DVCF file
    bool added_variant = false;
    if (!ctx_encountered_in_line (vb, dict_id_INFO_LUFT, NULL) && // note: no need to chech PRIM because LUFT and PRIM always appear together
        !ctx_encountered_in_line (vb, dict_id_INFO_LREJ, NULL) &&
        !ctx_encountered_in_line (vb, dict_id_INFO_PREJ, NULL)) {
        vcf_lo_seg_rollback_and_reject (vb, LO_ADDED_VARIANT, NULL); // note: we don't report this reject because it doesn't happen during --chain
        added_variant = true; // we added a REJX field in a variant that will be reconstructed in the current coordintes
    }

    // case: line originally had LIFTOVER or LIFTBACK. These can be fields from the txt files, or created by --chain
    bool has_luft    = ctx_encountered_in_line (vb, dict_id_INFO_LUFT, NULL);
    bool has_prim    = ctx_encountered_in_line (vb, dict_id_INFO_PRIM, NULL);
    bool has_lrej    = ctx_encountered_in_line (vb, dict_id_INFO_LREJ, NULL);
    bool has_prej    = ctx_encountered_in_line (vb, dict_id_INFO_PREJ, NULL);
    bool rolled_back = LO_IS_REJECTED (last_ostatus) && (has_luft || has_prim); // rejected in the Seg process
           
    // make sure we have either both LIFT/PRIM or both Lrej/Prej subfields in Primary and Luft
    ASSVCF ((has_luft && has_prim) || (has_lrej && has_prej), "%s", 
            vb->line_coords==DC_PRIMARY ? "Missing INFO/LUFT or INFO/Lrej subfield" : "Missing INFO/PRIM or INFO/Prej subfield");

    // case: --chain and INFO is '.' - remove the '.' as we are adding a DVCF field
    if (vb->info_items.len == 1 && FIRSTENT (InfoItem, vb->info_items)->name_len == 1 && *FIRSTENT (InfoItem, vb->info_items)->name == '.') {
        vb->info_items.len = 0;
        vb->recon_size--;
        vb->recon_size_luft--;
    }

    // dual coordinate line - we seg both options and vcf_piz_filter will decide which to render
    if (LO_IS_OK (last_ostatus)) {

        NEXTENT (InfoItem, vb->info_items) = (InfoItem) { 
            .name     = INFO_LUFT"=", 
            .name_len = INFO_DVCF_LEN + 1, // +1 for the '='
            .dict_id  = (DictId)dict_id_INFO_LUFT 
        };  
        
        NEXTENT (InfoItem, vb->info_items) = (InfoItem) { 
            .name     = INFO_PRIM"=", 
            .name_len = INFO_DVCF_LEN + 1, // +1 for the '='
            .dict_id  = (DictId)dict_id_INFO_PRIM 
        };  


        // case: --chain - we're adding ONE of these subfields to each of Primary and Luft reconstructions
        if (chain_is_loaded) {
            uint32_t growth = INFO_DVCF_LEN + 1 + (vb->info_items.len > 2); // +1 for '=', +1 for ';' if we already have item(s)
            vb->recon_size += growth;
            vb->recon_size_luft += growth;
        }
    }

    else { 
        NEXTENT (InfoItem, vb->info_items) = (InfoItem) { 
            .name     = INFO_LREJ"=", 
            .name_len = INFO_DVCF_LEN + 1, 
            .dict_id  = (DictId)dict_id_INFO_LREJ
        };

        NEXTENT (InfoItem, vb->info_items) = (InfoItem) { 
            .name     = INFO_PREJ"=", 
            .name_len = INFO_DVCF_LEN + 1, 
            .dict_id  = (DictId)dict_id_INFO_PREJ
        };

        // case: we added a REJX INFO field that wasn't in the TXT data: --chain or rolled back (see vcf_lo_seg_rollback_and_reject) or an added variant
        if (chain_is_loaded || rolled_back || added_variant) {
            uint32_t growth = INFO_DVCF_LEN + 1 + (vb->info_items.len > 2); // +1 for '=', +1 for ';' if we already have item(s) execpt for the DVCF items

            if (vb->line_coords == DC_PRIMARY) vb->recon_size += growth;
            else vb->recon_size_luft += growth;
        }
    }
}

static void vcf_seg_info_one_subfield (VBlockVCFP vb, DictId dict_id, const char *value, unsigned value_len)
{
    #define CALL(f) do { (f); not_yet_segged = false; } while(0) 
    #define CALL_WITH_FALLBACK(f) do { if (f (vb, ctx, value, value_len)) seg_by_ctx ((VBlockP)vb, value, value_len, ctx, value_len); \
                                       not_yet_segged = false; } while(0) 

    unsigned modified_len = value_len + 20;
    char modified[modified_len]; // used for 1. fields that are optimized 2. fields translated luft->primary. A_1 transformed 4.321e-03->0.004321
    
    Context *ctx = ctx_get_ctx (vb, dict_id), *other_ctx = NULL;
    bool not_yet_segged = true;
    
    // note: since we use modified for both optimization and luft_back - we currently don't support
    // subfields having both translators and optimization. This can be fixed if needed.
    #define ADJUST_FOR_MODIFIED do { \
        int32_t shrinkage = (int32_t)value_len - (int32_t)modified_len;\
        vb->recon_size -= shrinkage; \
        vb->recon_size_luft -= shrinkage; \
        value = modified; \
        value_len = modified_len; \
    } while(0)

    ctx->line_is_luft_trans = false; // initialize
    
    // --chain: if this is RendAlg=A_1 subfield in a REF<>ALT variant, convert a eg 4.31e-03 to e.g. 0.00431. This is to
    // ensure primary->luft->primary is lossless (4.31e-03 cannot be converted losslessly as we can't preserve format info)
    if (chain_is_loaded && ctx->luft_trans == VCF2VCF_A_1 && LO_IS_OK_SWITCH (last_ostatus) && 
        str_scientific_to_decimal (value, value_len, modified, &modified_len, NULL))
        ADJUST_FOR_MODIFIED;
        
    // Translatable item on a Luft line: attempt to lift-back the value, so we can seg it as primary
    if (vb->line_coords == DC_LUFT && needs_translation (ctx)) {

        // If cross-rendering to Primary is successful - proceed to Seg this value in primary coords, and assign a translator for reconstructing if --luft
        if (vcf_lo_seg_cross_render_to_primary (vb, ctx, value, value_len, modified, &modified_len)) {
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

        if (DT_FUNC(vb, translator)[ctx->luft_trans]((VBlockP)vb, ctx, (char *)value, value_len, true)) 
            ctx->line_is_luft_trans = true; // assign translator to this item in the container, to be activated with --luft
        else 
            REJECT_SUBFIELD (LO_INFO, ctx, "Cannot cross-render INFO subfield %s: \"%.*s\"", ctx->name, value_len, value);            
    }

    // ##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds of being a true variant versus being false under the trained gaussian mixture model">
    // Optimize VQSLOD
    if (dict_id.num == dict_id_INFO_VQSLOD && flag.optimize_VQSLOD && optimize_float_2_sig_dig (value, value_len, 0, modified, &modified_len)) 
        ADJUST_FOR_MODIFIED;
    
    // ##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
    // END is an alias of POS - they share the same delta stream - the next POS will be a delta vs this END)
    else if (dict_id.num == dict_id_INFO_END)
        CALL (vcf_seg_INFO_END (vb, ctx, value, value_len));

    // if SVLEN is negative, it is expected to be minus the delta between END and POS
    else if (dict_id.num == dict_id_INFO_SVLEN && vcf_seg_test_SVLEN (vb, value, value_len)) 
        CALL (seg_by_ctx (vb, ((char [2]){ SNIP_SPECIAL, VCF_SPECIAL_SVLEN }), 2, ctx, value_len));

    // ##INFO=<ID=BaseCounts,Number=4,Type=Integer,Description="Counts of each base">
    else if (dict_id.num == dict_id_INFO_BaseCounts) 
        CALL_WITH_FALLBACK (vcf_seg_INFO_BaseCounts);
    
    // ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    else if (dict_id.num == dict_id_INFO_DP) 
        CALL_WITH_FALLBACK (vcf_seg_INFO_DP);

    // Source File
    else if (dict_id.num == dict_id_INFO_SF) 
        CALL_WITH_FALLBACK (vcf_seg_INFO_SF_init);

    //##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    else if (dict_id.num == dict_id_INFO_AN) 
        seg_set_last_txt (vb, ctx, value, value_len, STORE_INT);

    // ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    else if (dict_id.num == dict_id_INFO_AF) {
        seg_set_last_txt (vb, ctx, value, value_len, STORE_FLOAT);
        ctx->keep_snip = true; // consumed by vcf_seg_FORMAT_AF
    }

    // ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    else if (dict_id.num == dict_id_INFO_AC)
        CALL (vcf_seg_INFO_AC (vb, ctx, value, value_len)); 

    else if (dict_id.num == dict_id_INFO_MLEAC && ctx_has_value_in_line (vb, dict_id_INFO_AC, &other_ctx)) 
        CALL (seg_delta_vs_other (vb, ctx, other_ctx, value, value_len, -1));

    // ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
    else if ((z_dual_coords  && ctx->luft_trans == VCF2VCF_ALLELE) || // if DVCF - apply to all fields (perhaps including AA) with RendAlg=ALLELE
             (!z_dual_coords && dict_id.num == dict_id_INFO_AA)) // apply to INFO/AA if not DVCF
        CALL_WITH_FALLBACK (vcf_seg_INFO_allele);

    // ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral">
    else if (dict_id.num == dict_id_INFO_CSQ) 
        CALL (vcf_seg_INFO_CSQ (vb, ctx, value, value_len));
    
    // ##INFO=<ID=DP_HIST,Number=R,Type=String,Description="Histogram for DP; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">
    // ##INFO=<ID=GQ_HIST,Number=R,Type=String,Description="Histogram for GQ; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">
    // ##INFO=<ID=AGE_HISTOGRAM_HET,Number=A,Type=String,Description="Histogram of ages of allele carriers; Bins: <30|30|35|40|45|50|55|60|65|70|75|80+">
    // ##INFO=<ID=AGE_HISTOGRAM_HOM,Number=A,Type=String,Description="Histogram of ages of homozygous allele carriers; Bins: <30|30|35|40|45|50|55|60|65|70|75|80+">
    else if (dict_id.num == dict_id_INFO_vep ||
             dict_id.num == dict_id_INFO_DP_HIST ||
             dict_id.num == dict_id_INFO_GQ_HIST ||
             dict_id.num == dict_id_INFO_AGE_HISTOGRAM_HET ||
             dict_id.num == dict_id_INFO_AGE_HISTOGRAM_HOM) 
        CALL (seg_array ((VBlockP)vb, ctx, ctx->did_i, value, value_len, ',', '|', false, true));

    else if (dict_id.num == dict_id_INFO_DP4) 
        CALL (seg_array ((VBlockP)vb, ctx, ctx->did_i, value, value_len, ',', 0, false, true));

    if (not_yet_segged) 
        seg_by_ctx (vb, value, value_len, ctx, value_len);
    
    ctx_set_encountered_in_line (ctx);
    
    #undef CALL
    #undef CALL_WITH_FALLBACK
}

static int sort_by_subfield_name (const void *a, const void *b)  
{ 
    InfoItem *ina = (InfoItem *)a;
    InfoItem *inb = (InfoItem *)b;
    
    return strncmp (ina->name, inb->name, MIN (ina->name_len, inb->name_len));
}

void vcf_seg_info_subfields (VBlockVCF *vb, const char *info_str, unsigned info_len)
{
    // parse the info string
    str_split (info_str, info_len, MAX_FIELDS-2, ';', pair, false); // -2 - leave room for LUFT + PRIM
    ASSVCF (n_pairs, "Too many INFO subfields, Genozip supports up to %u", MAX_FIELDS-2);

    buf_alloc (vb, &vb->info_items, 0, n_pairs + 2, InfoItem, CTX_GROWTH, "info_items");
    vb->info_items.len = 0; // reset from previous line

    int ac_i = -1; 
    InfoItem lift_ii = {}, rejt_ii = {};

    // pass 1: initialize info items + get indices of AC, and the DVCF items
    for (unsigned i=0; i < n_pairs; i++) {
        const char *equal_sign = memchr (pairs[i], '=', pair_lens[i]);
        unsigned name_len = (unsigned)(equal_sign - pairs[i]); // nonsense if no equal sign

        InfoItem ii = { .name_len  = equal_sign ? name_len + 1 : pair_lens[i], // including the '=' if there is one
                        .value     = equal_sign ? equal_sign + 1 : NULL,
                        .value_len = equal_sign ? pair_lens[i] - name_len - 1 : 0,
                        .dict_id   = equal_sign ? dict_id_make (pairs[i], name_len, DTYPE_1) : DICT_ID_NONE };
        
        // we make a copy the name, because vcf_seg_FORMAT_GT might overwrite the INFO field
        ASSVCF (ii.name_len <= sizeof (ii.name)-1, "INFO tag \"%s\" exceeds the maximum tag length supported by Genozip = %u", pairs[i], (int)sizeof (ii.name)-1);
        memcpy (ii.name, pairs[i], ii.name_len);

        ASSVCF (!z_dual_coords || 
                  (((ii.dict_id.num != dict_id_INFO_LUFT && ii.dict_id.num != dict_id_INFO_LREJ) || vb->line_coords == DC_PRIMARY) && 
                   ((ii.dict_id.num != dict_id_INFO_PRIM && ii.dict_id.num != dict_id_INFO_PREJ) || vb->line_coords == DC_LUFT)),
                "Not expecting %s in a %s-coordinate line", dis_dict_id_ex (ii.dict_id, true).s, coords_name (vb->line_coords));

        if (ii.dict_id.num == dict_id_INFO_AC) 
            ac_i = vb->info_items.len;

        else if (ii.dict_id.num == dict_id_INFO_LUFT || ii.dict_id.num == dict_id_INFO_PRIM) 
            { lift_ii = ii; continue; } // dont add LUFT and PRIM to Items yet

        else if (ii.dict_id.num == dict_id_INFO_LREJ || ii.dict_id.num == dict_id_INFO_PREJ) 
            { rejt_ii = ii; continue; } // dont add Lrej and Prej to Items yet

        NEXTENT (InfoItem, vb->info_items) = ii;
    }

    // case: we have a LUFT or PRIM item - Seg it now, but don't add it yet to InfoItems
    if (lift_ii.dict_id.num) { 
        vcf_lo_seg_INFO_LUFT_and_PRIM (vb, lift_ii.dict_id, lift_ii.value, lift_ii.value_len); 

        // case: we have both LIFT and REJT - could happen as a result of bcftools merge - discard the REJT for now, and let our Seg
        // decide if to reject it
        if (rejt_ii.dict_id.num) {
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
    else if (rejt_ii.dict_id.num)
        vcf_lo_seg_INFO_REJX (vb, rejt_ii.dict_id, rejt_ii.value, rejt_ii.value_len); 

    // pass 2: seg all subfields except AC (and PRIM/LUFT that weren't added)
    for (unsigned i=0; i < vb->info_items.len; i++) {
        InfoItem *ii = ENT (InfoItem, vb->info_items, i);
     
        if (ii->dict_id.num != dict_id_INFO_AC && ii->value)
            vcf_seg_info_one_subfield (vb, ii->dict_id, ii->value, ii->value_len);
    }
    
    // last, seg AC (delayed, as we needed to seg AN and AF before)
    if (ac_i >= 0) {
        InfoItem *ii = ENT (InfoItem, vb->info_items, ac_i);
        vcf_seg_info_one_subfield (vb, ii->dict_id, ii->value, ii->value_len);
    }
}

// Adds the DVCF items according to ostatus, finalizes INFO/SF and segs the INFO container
void vcf_finalize_seg_info (VBlockVCF *vb)
{
    Container con = { .repeats             = 1, 
                      .drop_final_item_sep = true,
                      .filter_items        = z_dual_coords,   // vcf_piz_filter chooses which DVCF item to show based on flag.luft 
                      .callback            = z_dual_coords }; // vcf_piz_container_cb appends oSTATUS to INFO if requested 

    // seg INFO/SF, if there is one
    if (vb->sf_txt.len) vcf_seg_INFO_SF_seg (vb);

    // now that we segged all INFO and FORMAT subfields, we have the final ostatus and can add the DVCF items
    if (z_dual_coords)
        vcf_seg_info_add_DVCF_to_InfoItems (vb);

    con_set_nitems (con, vb->info_items.len);

    // if requested, we will re-sort the info fields in alphabetical order. This will result less words in the dictionary
    // thereby both improving compression and improving --regions speed. 
    if (flag.optimize_sort && con_nitems(con) > 1) 
        qsort (vb->info_items.data, vb->info_items.len, sizeof(InfoItem), sort_by_subfield_name);

    char prefixes[CONTAINER_MAX_PREFIXES_LEN];  // these are the Container prefixes
    prefixes[0] = prefixes[1] = CON_PREFIX_SEP; // initial CON_PREFIX_SEP follow by seperator of empty Container-wide prefix
    unsigned prefixes_len = 2;

    // Populate the Container 
    uint32_t total_names_len=0;
    for (unsigned i=0; i < con_nitems(con); i++) {
        // Set the Container item and find (or create) a context for this name
        InfoItem *ii = ENT (InfoItem, vb->info_items, i);
        con.items[i] = (ContainerItem){ .dict_id   = (ii->dict_id.num == dict_id_INFO_END) ? (DictId)dict_id_fields[VCF_POS] : ii->dict_id,
                                        .seperator = { ';' } }; 

        // if we're preparing a dual-coordinate VCF and this line as a REF-ALT switch - assign the liftover-translator for this item,
        // which was calculated in vcf_header_set_translators
        Context *ii_ctx = ctx_get_existing_ctx (vb, ii->dict_id);
        if (ii_ctx && ii_ctx->line_is_luft_trans) { // item was segged in Primary coords and needs a luft translator to be reconstruced in --luft
            con.items[i].translator = ii_ctx->luft_trans;

            if (ii_ctx->luft_trans == VCF2VCF_A_AN)
                ii_ctx->flags.store = STORE_INT; // consumed by vcf_piz_luft_A_AN
        }
            
        // add to the prefixes
        ASSVCF (prefixes_len + ii->name_len + 1 <= CONTAINER_MAX_PREFIXES_LEN, 
                "INFO contains tag names that, combined (including the '='), exceed the maximum of %u characters", CONTAINER_MAX_PREFIXES_LEN);

        memcpy (&prefixes[prefixes_len], ii->name, ii->name_len);
        prefixes_len += ii->name_len;
        prefixes[prefixes_len++] = CON_PREFIX_SEP;

        // don't include LIFTBACK or LIFTREJT because they are not reconstructed by default (genounzip) 
        // note: vcf_lo_seg_INFO_REJX / vcf_lo_seg_INFO_LUFT_and_PRIM already verified that this is a dual-coord file
        if (ii->dict_id.num != dict_id_INFO_PRIM && ii->dict_id.num != dict_id_INFO_PREJ)
            total_names_len += ii->name_len + 1; // +1 for ; \t or \n separator
    }

    container_seg_by_ctx (vb, CTX(VCF_INFO), &con, prefixes, prefixes_len, total_names_len /* names inc. = and separator */);
}

