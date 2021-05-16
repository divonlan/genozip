// ------------------------------------------------------------------
//   vcf_info.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

typedef struct { const char *name, *value; 
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
void vcf_seg_INFO_SF_one_sample (VBlockVCF *vb, unsigned sample_i)
{
    // case: no more SF values left to compare - we ignore this sample
    while (vb->sf_txt.next < vb->sf_txt.len) {

        buf_alloc (vb, &vb->sf_snip, 10, 0, char, 2, "sf_snip");

        char *sf_one_value = ENT (char, vb->sf_txt, vb->sf_txt.next); 
        char *after;
        int32_t value = strtol (sf_one_value, &after, 10);

        int32_t adjusted_sample_i = (int32_t)(sample_i + adjustment); // adjustment is the number of values in SF that are not in samples

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

// ---------------
// INFO/BaseCounts
// ---------------

// ##INFO=<ID=BaseCounts,Number=4,Type=Integer,Description="Counts of each base">
// Sorts BaseCounts vector with REF bases first followed by ALT bases, as they are expected to have the highest values
static bool vcf_seg_INFO_BaseCounts (VBlockVCF *vb, Context *ctx_basecounts, const char *value, int value_len) // returns true if caller still needs to seg 
{
    // if XSTRAND, bases have changed and Genozip doesn't currently liftover BaseCounts - reject liftover
    if (chain_is_loaded && *vb->contexts[VCF_oXSTRAND].last_snip == 'X')
        vcf_lo_seg_rollback_and_reject (vb, LO_INFO, ctx_basecounts);

    if (vb->contexts[VCF_REFALT].last_snip_len != 3) return true; // not simple two bases - caller should seg

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

    const char *refalt = vb->contexts[VCF_REFALT].last_snip;

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
    reconstruct_peek (vb, &vb->contexts[VCF_REFALT], &refalt, &refalt_len);

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

    Context *sf_ctxs[MAX_SUBFIELDS] = {};

    uint32_t item_i=0;
    const char *item_start = field;
    for (uint32_t i=0; i < field_len+1; i++) {
        if (i == field_len || field[i] == ',' || field[i] == '|') { // end of item
            if (item_i == con_nitems(con)) {
                ASSVCF (!con.repeats, 
                        "expecting all repeats of %s to have the same number of items, %u, as the first repeat, but repeat %u (0-based) has more: %.*s", 
                        vep_ctx->name, con_nitems(con), con.repeats, field_len, field);

                ASSVCF (item_i < MIN (126, MAX_SUBFIELDS), "exceeded the max number of %s items=%u", 
                        vep_ctx->name, MIN (126, MAX_SUBFIELDS)); // the 126 constraint is just the context naming scheme
                
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
    ctx->last_value.i = new_value->i = reconstruct_peek_(vb, dict_id_INFO_AN, 0, 0).i * reconstruct_peek_(vb, dict_id_INFO_AF, 0, 0).f;
    RECONSTRUCT_INT (new_value->i); 

    return true;
}

// Lift-over translator for INFO/AC fields, IF it is bi-allelic and we have a ALT<>REF switch AND we have an AN and AF field
// We change the AC to AN - AC and return true if successful 
TRANSLATOR_FUNC (vcf_piz_luft_A_AN)
{
    int64_t an, ac;

    if (command == ZIP) {
        Context *an_ctx;
        if (!ctx_has_value_in_line (vb, dict_id_INFO_AN, &an_ctx) ||
            !str_get_int (recon, recon_len, &ac))
            return false; // in Seg, AC is always segged last, so this means there is no AN in the line for sure
        
        if (validate_only) return true;  // Yay! AC can be lifted - all it needs in AN in the line, which it has 

        an = an_ctx->last_value.i;
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

// ---------------------
// INFO/AF and FORMAT/AF
// ---------------------

// Lift-over translator for INFO/AF and FORMAT/AF fields, IF it is bi-allelic and we have a ALT<>REF switch.
// We change the probability value to 1-AF
// returns true if successful (return value used only if validate_only)
TRANSLATOR_FUNC (vcf_piz_luft_A_1)
{
    // if item format is inconsistent with AF being a probability value - we won't translate it
    double af = str_get_positive_float (recon, recon_len);
    
    if (af < 0 || af > 1) return false;
    
    if (validate_only) return true; 

    char format[20];
    str_get_float_format (recon, recon_len, format);

    vb->txt_data.len -= recon_len;
    char af_str[30];
    sprintf (af_str, format, 1 - af);
    RECONSTRUCT (af_str, strlen (af_str)); // careful not to use bufprintf as it adds a \0 and we are required to translate in-place for all FORMAT fields
    
    return true;
}

// --------
// INFO/END
// --------

static void vcf_seg_INFO_END (VBlockVCFP vb, Context *end_ctx, const char *end_str, unsigned end_len) // note: ctx is INFO/END *not* POS (despite being an alias)
{
    // END is an alias of POS
    seg_pos_field ((VBlockP)vb, VCF_POS, VCF_POS, true, true, 0, end_str, end_len, 0, end_len);
    
    // case --chain: if we have lifted-over POS (as primary POS field or in INFO/LIFTBACK), 
    // check that lifting-over of END is delta-encoded and is lifted over to the same, non-xstrand, Chain alignment, and reject if not
    if (chain_is_loaded && LO_IS_OK (last_ostatus)) { 

        bool bad_mapping = !vb->contexts[VCF_POS].last_delta; // reject if END was not stored as a delta - our translator won't work
        bool xstrand = false;
        PosType opos = 0;

        if (!bad_mapping) {
            uint32_t end_aln_i;
            LiftOverStatus end_ostatus = vcf_lo_get_liftover_coords (vb, 0, &opos, &xstrand, &end_aln_i); // uses vb->chrom_node_index and POS->last_int
            bad_mapping = LO_IS_REJECTED (end_ostatus) || end_aln_i != vb->pos_aln_i;
        }

        if (bad_mapping || xstrand) 
            vcf_lo_seg_rollback_and_reject (vb, LO_INFO, end_ctx);
    }
}

// END data resides in POS (its an alias), but has a different translator as its a different container item. 
// For END, we didn't add an oPOS entry, because we can't consume it when showing Primary. Instead, we do delta arithmetic.
// returns true if successful (return value used only if validate_only)
TRANSLATOR_FUNC (vcf_piz_luft_END)
{
    // ZIP liftover validation: postpone to vcf_seg_INFO_END
    if (validate_only) return true; 

  //  if (vb->rejects_vb || LO_IS_REJECTED (vb->last_index (VCF_oSTATUS))) return false; // not translated if rejected

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
        int64_t end   = ctx->last_value.i;
        int64_t delta = ctx->last_delta;
        int64_t pos   = end - delta;

        translated_end = vb->last_int (VCF_oPOS) + vb->last_delta (VCF_POS); // delta for generated this END value (END - POS)

        ((VBlockVCFP)vb)->last_end_line_i = vb->line_i; // so vcf_piz_special_COPY_POS knows that END was reconstructed
    }

    // re-reconstruct END
    vb->txt_data.len -= recon_len; 
    RECONSTRUCT_INT (translated_end);
    
    return true;    
}

// Called to reconstruct the POS subfield of INFO/LIFTBACK, handling the possibility of INFO/END
SPECIAL_RECONSTRUCTOR (vcf_piz_special_COPY_POS)
{
    bool has_end = ((VBlockVCFP)vb)->last_end_line_i == vb->line_i; // true if INFO/END was encountered

    Context *pos_ctx = &vb->contexts[VCF_POS];
    int64_t pos;

    if (has_end) {
        int64_t end   = pos_ctx->last_value.i;
        int64_t delta = pos_ctx->last_delta;
        pos = end - delta;
    }
    else
        pos = pos_ctx->last_value.i;;

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

    int64_t last_delta = vb->contexts[VCF_POS].last_delta; // INFO_END is an alias of POS - so the last delta would be between END and POS
    return last_delta == -svlen;
}

// the case where SVLEN is minus the delta between END and POS
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_SVLEN)
{
    if (!reconstruct) goto done;

    int64_t value = -vb->contexts[VCF_POS].last_delta; // END is a alias of POS - they share the same data stream - so last_delta would be the delta between END and POS
    char str[30];
    unsigned str_len = str_int (value, str);
    RECONSTRUCT (str, str_len);

done:
    return false; // no new value
}

// --------------
// INFO container
// --------------

// for dual coordinate files (Primary, Luft and --chain) - add LIFTXXXX depending on ostatus (run after
// all INFO and FORMAT fields, so ostatus is final)
static void vcf_seg_info_add_LIFTXXXX_items (VBlockVCF *vb)
{
    // case: line originally had LIFTOVER or LIFTBACK, but was rejected in the Seg process (perhaps new or modified
    // fields after created by --chain)
    bool has_lo      = ctx_encountered_in_line (vb, dict_id_INFO_LIFTOVER, NULL);
    bool has_lb      = ctx_encountered_in_line (vb, dict_id_INFO_LIFTBACK, NULL);
    bool has_ro      = ctx_encountered_in_line (vb, dict_id_INFO_REJTOVER, NULL);
    bool has_rb      = ctx_encountered_in_line (vb, dict_id_INFO_REJTBACK, NULL);
    bool rolled_back = LO_IS_REJECTED (last_ostatus) && (has_lo || has_lb);
   
    // make sure we have either both XXXXOVER or both XXXXREJT subfields in Primary and Luft
    ASSVCF0 ((has_lo && has_lb) || (has_ro && has_rb), "Missing INFO/LIFTXXXX or INFO/REJTXXXX subfield(s)");

    // case: --chain and INFO is '.' - remove the '.' as we are adding INFO/LIFTXXXX
    if (vb->info_items.len == 1 && FIRSTENT (InfoItem, vb->info_items)->name_len == 1 && *FIRSTENT (InfoItem, vb->info_items)->name == '.') {
        vb->info_items.len = 0;
        vb->vb_data_size--;
    }

    // dual coordinate line - we seg both options and vcf_piz_filter will decide which to render
    if (LO_IS_OK (last_ostatus)) {

        NEXTENT (InfoItem, vb->info_items) = (InfoItem) { 
            .name     = INFO_LIFTOVER"=", 
            .name_len = INFO_LIFTOVER_LEN + 1, // +1 for the '='
            .dict_id  = (DictId)dict_id_INFO_LIFTOVER 
        };  
        
        if (chain_is_loaded) // --chain - we're adding this subfield to the default reconstruction
            vb->vb_data_size += INFO_LIFTOVER_LEN + 1 + (vb->info_items.len > 1); // +1 for '=', +1 for ';' if we already have item(s).

        NEXTENT (InfoItem, vb->info_items) = (InfoItem) { 
            .name     = INFO_LIFTBACK"=", 
            .name_len = INFO_LIFTOVER_LEN + 1, // +1 for the '='
            .dict_id  = (DictId)dict_id_INFO_LIFTBACK 
        };  
        // note: we don't increase vb_data_size for LIFTBACK because its not displayed in the default reconstruction of the file
    }

    else { 
        NEXTENT (InfoItem, vb->info_items) = (InfoItem) { 
            .name     = INFO_REJTOVER"=", 
            .name_len = INFO_LIFTOVER_LEN + 1, 
            .dict_id  = (DictId)dict_id_INFO_REJTOVER
        };

        NEXTENT (InfoItem, vb->info_items) = (InfoItem) { 
            .name     = INFO_REJTBACK"=", 
            .name_len = INFO_LIFTOVER_LEN + 1, 
            .dict_id  = (DictId)dict_id_INFO_REJTBACK
        };

        // case: --chain or rolled back (see vcf_lo_seg_rollback_and_reject) - we're adding REJTOVER subfield to the default (genounzip) reconstruction
        if (chain_is_loaded || rolled_back) 
            vb->vb_data_size += INFO_LIFTOVER_LEN + 1 + (vb->info_items.len > 2); // +1 for '=', +1 for ';' if we already have item(s) execpt for the two LIFTXXXX or the two REJTXXXX
    }
}

static void vcf_seg_info_one_subfield (VBlockVCFP vb, DictId dict_id, const char *value, unsigned value_len)
{
    #define CALL(f) do { (f); not_yet_segged = false; } while(0) 
    #define CALL_WITH_FALLBACK(f) do { if (f) seg_integer_or_not ((VBlockP)vb, ctx, value, value_len, value_len); \
                                       not_yet_segged = false; } while(0) 

    char modified_snip[OPTIMIZE_MAX_SNIP_LEN]; // used for 1. fields that are optimized 2. fields translated luft->primary
    unsigned modified_snip_len;
    Context *ctx = dict_id.num ? ctx_get_ctx (vb, dict_id) : NULL, *other_ctx = NULL;
    bool not_yet_segged = true;
    
    // note: since we use modified_snip for both optimization and luft_back - we currently don't support
    // subfields having both translators and optimization. This can be fixed if needed.
    #define ADJUST_FOR_MODIFIED do { \
        vb->vb_data_size -= (int)value_len - (int)modified_snip_len; \
        value = modified_snip; \
        value_len = modified_snip_len; \
    } while(0)

    // lift-back the value, if we're segging a Luft file
    if (vb->line_coords == DC_LUFT && ctx->luft_trans &&
        (last_ostatus == LO_OK_REF_ALT_SWTCH || (dict_id.num == dict_id_INFO_END && LO_IS_OK (last_ostatus)))) {
        if (vcf_lo_seg_lift_back_to_primary (vb, ctx, value, value_len, modified_snip, &modified_snip_len)) {
            value = modified_snip; 
            value_len = modified_snip_len; 
        } 
        else {
            bool luft_only_line = true;
        }
    }

    // validate that the primary value (as received from caller or lifted back) can be luft-translated 
    // note: looks at snips before optimization, we're counting on the optimization not changing the validation outcome
    if (vb->line_coords == DC_PRIMARY && z_dual_coords && ctx->luft_trans && last_ostatus == LO_OK_REF_ALT_SWTCH &&
        !(DT_FUNC(vb, translator)[ctx->luft_trans]((VBlockP)vb, ctx, (char *)value, value_len, true)))
        vcf_lo_seg_rollback_and_reject (vb, LO_INFO, ctx);

    // ##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds of being a true variant versus being false under the trained gaussian mixture model">
    // Optimize VQSLOD
    if (flag.optimize_VQSLOD && (dict_id.num == dict_id_INFO_VQSLOD) &&
        optimize_float_2_sig_dig (value, value_len, 0, modified_snip, &modified_snip_len)) 
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
        CALL_WITH_FALLBACK (vcf_seg_INFO_BaseCounts (vb, ctx, value, value_len));
    
    // ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    else if (dict_id.num == dict_id_INFO_DP) 
        CALL_WITH_FALLBACK (vcf_seg_INFO_DP (vb, ctx, value, value_len));

    // Source File
    else if (dict_id.num == dict_id_INFO_SF) 
        CALL_WITH_FALLBACK (vcf_seg_INFO_SF_init (vb, ctx, value, value_len));

    //##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    else if (dict_id.num == dict_id_INFO_AN) 
        seg_set_last_txt (vb, ctx, value, value_len, STORE_INT);

    // ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    else if (dict_id.num == dict_id_INFO_AF) 
        seg_set_last_txt (vb, ctx, value, value_len, STORE_FLOAT);

    // ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    else if (dict_id.num == dict_id_INFO_AC)
        CALL (vcf_seg_INFO_AC (vb, ctx, value, value_len)); 

    else if (dict_id.num == dict_id_INFO_MLEAC && ctx_has_value_in_line (vb, dict_id_INFO_AC, &other_ctx)) 
        CALL (seg_delta_vs_other (vb, ctx, other_ctx, value, value_len, -1));

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
        CALL (seg_array ((VBlockP)vb, ctx, ctx->did_i, value, value_len, ',', '|', false));

    else if (dict_id.num == dict_id_INFO_DP4) 
        CALL (seg_array ((VBlockP)vb, ctx, ctx->did_i, value, value_len, ',', 0, false));

    if (not_yet_segged) 
        seg_integer_or_not ((VBlockP)vb, ctx, value, value_len, value_len);
    
    #undef CALL
}

static int sort_by_subfield_name (const void *a, const void *b)  
{ 
    InfoItem *ina = (InfoItem *)a;
    InfoItem *inb = (InfoItem *)b;
    
    return strncmp (ina->name, inb->name, MIN (ina->name_len, inb->name_len));
}

void vcf_seg_info_subfields (VBlockVCF *vb, const char *info_str, unsigned info_len)
{
    const char *pairs[MAX_SUBFIELDS];
    unsigned pair_lens[MAX_SUBFIELDS];

    // parse the info string
    unsigned num_items = str_split (info_str, info_len, MAX_SUBFIELDS-2, ';', pairs, pair_lens, false, 0); // -2 - leave room for LIFTXXXX
    ASSVCF (num_items, "Too many INFO subfields, Genozip supports up to %u", MAX_SUBFIELDS-2);

    buf_alloc (vb, &vb->info_items, 0, num_items + 2, InfoItem, 2, "info_items");
    vb->info_items.len = 0; // reset from previous line

    int ac_i = -1;

    // pass 1: initialize info items + set ac_i + seg LIFTXXXX (but don't add to info_items yet) 
    for (unsigned i=0; i < num_items; i++) {
        const char *equal_sign = memchr (pairs[i], '=', pair_lens[i]);
        unsigned name_len = (unsigned)(equal_sign - pairs[i]); // nonsense if no equal sign

        InfoItem ii = { .name      = pairs[i],
                        .name_len  = equal_sign ? name_len + 1 : pair_lens[i], // including the '=' if there is one
                        .value     = equal_sign ? equal_sign + 1 : NULL,
                        .value_len = equal_sign ? pair_lens[i] - name_len - 1 : 0,
                        .dict_id   = equal_sign ? dict_id_make (pairs[i], name_len, DTYPE_1) : DICT_ID_NONE };

        if (ii.dict_id.num == dict_id_INFO_AC) 
            ac_i = vb->info_items.len;
        
        else if (ii.dict_id.num == dict_id_INFO_LIFTOVER || ii.dict_id.num == dict_id_INFO_LIFTBACK) { 
            vcf_lo_seg_INFO_LIFTXXXX (vb, ii.dict_id, ii.value, ii.value_len); 
            continue; 
        }
        
        else if (ii.dict_id.num == dict_id_INFO_REJTOVER || ii.dict_id.num == dict_id_INFO_REJTBACK) { 
            vcf_lo_seg_INFO_REJTXXXX (vb, ii.dict_id, ii.value, ii.value_len); 
            continue; 
        }

        NEXTENT (InfoItem, vb->info_items) = ii;
    }

    // pass 2: seg all subfields except AC (and LIFTXXXX that weren't added)
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

// Adds INFO/LIFTXXXX and INFO/REJTXXXX according to ostatus, finalizes INFO/SF and segs the INFO container
void vcf_finalize_seg_info (VBlockVCF *vb)
{
    Container con = { .repeats             = 1, 
                      .drop_final_item_sep = true,
                      .filter_items        = z_dual_coords }; // vcf_piz_filter chooses which LIFTXXXX or REJTXXXX to show based on flag.luft 

    // seg INFO/SF, if there is one
    if (vb->sf_txt.len) vcf_seg_INFO_SF_seg (vb);

    // now that we segged all INFO and FORMAT subfields, we have the final ostatus and can add INFOXXXX
    if (z_dual_coords)
        vcf_seg_info_add_LIFTXXXX_items (vb);

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
        Context *ii_ctx;
        if ((last_ostatus == LO_OK_REF_ALT_SWTCH || (ii->dict_id.num == dict_id_INFO_END && LO_IS_OK (last_ostatus))) && 
            (ii_ctx = ctx_get_existing_ctx (vb, ii->dict_id))) // translator exists and validated for this line
            con.items[i].translator = ii_ctx->luft_trans;
            
        // add to the prefixes
        ASSVCF (prefixes_len + ii->name_len + 1 <= CONTAINER_MAX_PREFIXES_LEN, 
                "INFO contains tag names that, combined (including the '='), exceed the maximum of %u characters", CONTAINER_MAX_PREFIXES_LEN);

        memcpy (&prefixes[prefixes_len], ii->name, ii->name_len);
        prefixes_len += ii->name_len;
        prefixes[prefixes_len++] = CON_PREFIX_SEP;

        // don't include LIFTBACK or LIFTREJT because they are not reconstructed by default (genounzip) 
        // note: vcf_lo_seg_INFO_REJTXXXX / vcf_lo_seg_INFO_LIFTXXXX already verified that this is a dual-coord file
        if (ii->dict_id.num != dict_id_INFO_LIFTBACK && ii->dict_id.num != dict_id_INFO_REJTBACK)
            total_names_len += ii->name_len + 1; // +1 for ; \t or \n separator
    }

    container_seg_by_ctx (vb, &vb->contexts[VCF_INFO], &con, prefixes, prefixes_len, total_names_len /* names inc. = and separator */);
}

