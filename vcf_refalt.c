// ------------------------------------------------------------------
//   vcf_refalt.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <math.h>
#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "strings.h"
#include "codec.h"
#include "reconstruct.h"
#include "dict_id.h"
#include "file.h"

// ---------
// Seg stuff
// ---------

// optimize REF and ALT, for simple one-character REF/ALT (i.e. mostly a SNP or no-variant)
static void vcf_refalt_seg_ref_alt_snp (VBlockVCFP vb, char main_ref, char main_alt)
{
    char new_ref=0, new_alt=0;

    // if we have a reference, we use it (we treat the reference as PRIMARY)
    if ((flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) && vb->line_coords == DC_PRIMARY) {
        PosType pos = vb->last_int(VCF_POS);

        RefLock lock;
        Range *range = ref_seg_get_locked_range ((VBlockP)vb, vb->chrom_node_index, pos, 1, NULL, &lock);
        uint32_t index_within_range = pos - range->first_pos;

        ref_assert_nucleotide_available (range, pos);
        char ref = ref_get_nucleotide (range, index_within_range);

        if (main_ref == ref) new_ref = '-'; // this should always be the case...
        if (main_alt == ref) new_alt = '-'; 

        if (flag.reference == REF_EXT_STORE)
            bit_array_set (&range->is_set, index_within_range);

        ref_unlock (lock);
    }

    // replace the most common SNP with +
    // based on counting simple SNPs from from chr22 of 1000 genome project phase 3:
    // G to: A=239681 C=46244 T=44084
    // C to: T=238728 G=46508 A=43685
    // A to: G=111967 C=30006 T=26335
    // T to: C=111539 G=29504 A=25599

    if      (main_alt == 'A' && main_ref == 'G') new_alt = '+';
    else if (main_alt == 'C' && main_ref == 'T') new_alt = '+';
    else if (main_alt == 'G' && main_ref == 'A') new_alt = '+';
    else if (main_alt == 'T' && main_ref == 'C') new_alt = '+';

    // if anything was done, we create a "special" snip
    if (new_ref || new_alt) {
        char refalt_special[4] = { SNIP_SPECIAL, VCF_SPECIAL_main_REFALT };
        refalt_special[2] = new_ref ? new_ref : main_ref;
        refalt_special[3] = new_alt ? new_alt : main_alt;

        seg_by_did_i (vb, refalt_special, sizeof(refalt_special), SEL (VCF_REFALT, VCF_oREFALT), 0); // we do the account in vcf_refalt_seg_main_ref_alt
    }
    // if not - just the normal snip
    else {
        char refalt_normal[3] = { 0, '\t', 0 };
        refalt_normal[0] = main_ref;
        refalt_normal[2] = main_alt;

        seg_by_did_i (vb, refalt_normal, sizeof(refalt_normal), SEL (VCF_REFALT, VCF_oREFALT), 0); // we do the account in vcf_refalt_seg_main_ref_alt
    }
}

void vcf_refalt_seg_main_ref_alt (VBlockVCFP vb, const char *ref, unsigned ref_len, const char *alt, unsigned alt_len)
{
    // optimize ref/alt in the common case of single-character
    if (ref_len == 1 && alt_len == 1) 
        vcf_refalt_seg_ref_alt_snp (vb, *ref, *alt);

    else {
        char ref_alt[ref_len + 1 + alt_len];
        memcpy (ref_alt, ref, ref_len);
        ref_alt[ref_len] = '\t';
        memcpy (&ref_alt[ref_len+1], alt, alt_len);

        seg_by_did_i (vb, ref_alt, ref_len + alt_len + 1, SEL (VCF_REFALT, VCF_oREFALT), 0);
    }

    if (vb->line_coords == DC_PRIMARY) // in the default reconstruction, this REFALT is in the main fields along with 2 tabs
        vb->contexts[VCF_REFALT].txt_len += ref_len + alt_len + 2; \

    else // in the default reconstruction, oREF is in the INFO/LIFTOVER vector
        vb->contexts[VCF_LIFT_REF].txt_len += ref_len;

    // note: vb->recon_size doesn't change if txt_files.coords==LUFT, because currently Genozip only supports dual coords if the size
    // of REF and oREF are the same, and ALT and oALT.
}

// --------------
// Liftover stuff
// --------------

RefAltEquals vcf_refalt_ref_equals_ref2_or_alt (char ref, char ref2, const char *alt, unsigned alt_len, bool is_xstrand)
{
    if (ref == '.') return EQUALS_MISSING;

    char upper_ref  = UPPER_CASE (ref);
    char upper_ref2 = UPPER_CASE (ref2);
    char upper_alt  = alt_len == 1 ? UPPER_CASE (*alt) : 0;

    if ((!is_xstrand && upper_ref == upper_ref2) || 
        ( is_xstrand && REVCOMP[(int)upper_ref] == upper_ref2)) 
        return EQUALS_REF2;

    else if ((!is_xstrand && upper_alt == upper_ref2) || (is_xstrand && REVCOMP[(int)upper_alt] == upper_ref2)) 
        return EQUALS_ALT;

    else
        return EQUALS_NEITHER;
}

static inline char liftover_seg_get_oref (VBlockP vb, WordIndex ochrom, PosType opos)
{
    if (!opos) return '.'; // missing

    Range *range = ref_seg_get_locked_range (vb, ochrom, opos, 1, ENT (char, vb->txt_data, vb->line_start), NULL);
    ASSERT (range, "Failed to find range for ochrom=%d", ochrom);

    uint32_t index_within_range = opos - range->first_pos;

    ref_assert_nucleotide_available (range, opos);
    char ref = ref_get_nucleotide (range, index_within_range);

    return ref;
}

// Segging a NON-dual-coordinates file called when genozip --chain
// Seg: set oref and update ostatus
LiftOverStatus vcf_refalt_check_oref (VBlockVCFP vb, const ZipDataLineVCF *dl, bool is_xstrand, unsigned *oref_len)
{
    LiftOverStatus ostatus;

    if (vb->main_ref_len > 1) { *oref_len=0; return LO_REF_TOO_LONG; } // Genozip can only do --chain with REF of length 1 (even without change). It can liftover / liftback longer REFs in some conditions.

    char ref = *vb->main_refalt;
    const char *alt = &vb->main_refalt[vb->main_ref_len + 1]; // +1 for \t
    char oref = liftover_seg_get_oref ((VBlockP)vb, dl->chrom_index[1], dl->pos[1]);
    
    // oref preserves the case of ref (REF/ALT are case insensitive per VCF spec, but we want to ensure binary losslessness)
    if (IS_SLETTER (ref)) oref = LOWER_CASE(oref);
    else                  oref = UPPER_CASE(oref);

    switch (vcf_refalt_ref_equals_ref2_or_alt (ref, oref, alt, vb->main_alt_len, is_xstrand)) {
        case EQUALS_MISSING :
        case EQUALS_REF2    : ostatus = !is_xstrand || vb->main_alt_len==1 ? LO_OK_REF_SAME : LO_ALT_LONG_XSTRAND; break;
        case EQUALS_ALT     : ostatus = vb->main_alt_len==1 ? LO_OK_REF_ALT_SWTCH : LO_ALT_LONG_SWITCH; break; // Genozip can only switch REF<>ALT for an ALT of length 1
        case EQUALS_NEITHER : ostatus = LO_REF_CHNG_NOT_ALT; break; // Genozip can only liftover a changed REF if its switched with the ALT
        default             : ASSVCF0 (0, "should never reach here"); 
    }
    
    *oref_len = ostatus == LO_OK_REF_SAME      ? vb->main_ref_len
              : ostatus == LO_OK_REF_ALT_SWTCH ? vb->main_alt_len
              :                                  0;
    
    return ostatus;
}

// ---------
// PIZ stuff
// ---------

// Sometimes called to reconstruct the "main" refalt (main AT THE TIME OF SEGGING), from reference and/or common SNPs.
// The main refalt may also be a normal snip.
SPECIAL_RECONSTRUCTOR (vcf_piz_special_main_REFALT)
{
    if (!reconstruct) goto done;

    ASSPIZ (snip_len==2, "expecting snip_len=2 but seeing %u", snip_len);

    // snip is 3 characters - REF, \t, ALT
    char ref_alt[3] = { 0, '\t', 0 };
    char ref_value = 0;
    
    if (snip[0] == '-' || snip[1] == '-') { 
        PosType pos = vb->last_int (VCF_POS);

        const Range *range = ref_piz_get_range (vb, pos, 1);
        ASSPIZ (range, "failed to find range for chrom='%s' pos=%"PRId64, vb->chrom_name, pos);
        
        uint32_t idx = pos - range->first_pos;
        ASSPIZ (ref_is_nucleotide_set (range, idx), "reference is not set: chrom=%.*s pos=%"PRId64, range->chrom_name_len, range->chrom_name, pos);
        ref_value = ref_get_nucleotide (range, idx);
    }

    // recover ref
    if (snip[0] == '-') 
        ref_alt[0] = ref_value;
    else 
        ref_alt[0] = snip[0];

    // recover alt
    if (snip[1] == '+') { // the alt has the most common value for a SNP
        if      (ref_alt[0] == 'A') ref_alt[2] = 'G';
        else if (ref_alt[0] == 'C') ref_alt[2] = 'T';
        else if (ref_alt[0] == 'G') ref_alt[2] = 'A';
        else if (ref_alt[0] == 'T') ref_alt[2] = 'C';
    }
    else if (snip[1] == '-')  // the alt has the reference value
        ref_alt[2] = ref_value;

    else // the alt is specified verbatim
        ref_alt[2] = snip[1];

    RECONSTRUCT (ref_alt, sizeof (ref_alt));

done:
    return false; // no new value
}   

// Always called for reconstructing the "other" REFALT ("other" coordinates AT THE TIME OF SEGGING) (i.e. oREFALT if we segged a PRIMARY
// and REFALT if we segged a LUFT). This function may be called:
// 1. from a toplevel container as a main VCF field 
// 2. from vcf_piz_special_LIFT_REF when reconstructing the INFO/LIFTXXXX REF field
// Note that the function might be called before or after the main (at seg time) REFALT, depending on whether we are reconstructing in the same coordinates
// as we segged, or the reverse.
SPECIAL_RECONSTRUCTOR (vcf_piz_special_other_REFALT)
{
    if (!reconstruct) return false;

    Context *other_refalt_ctx = &vb->contexts[ctx->did_i == VCF_REFALT ? VCF_oREFALT : VCF_REFALT];

    const char *other_refalt; 
    unsigned other_refalt_len;
    reconstruct_peek (vb, other_refalt_ctx, &other_refalt, &other_refalt_len);
    
    const char *ref_alt[2]; unsigned ref_alt_lens[2];
    ASSPIZ (str_split (other_refalt, other_refalt_len, 2, '\t', ref_alt, ref_alt_lens, true, 0),
            "expecting one tab in the %s snip: \"%.*s\"", other_refalt_ctx->name, other_refalt_len, other_refalt);

    // make private copy as ref_alt is in txt_data and will be overwritten
    unsigned ref_len = ref_alt_lens[0]; char ref[ref_len]; memcpy (ref, ref_alt[0], ref_len); 
    unsigned alt_len = ref_alt_lens[1]; char alt[alt_len]; memcpy (alt, ref_alt[1], alt_len);

    const char *xstrand;
    reconstruct_peek (vb, &vb->contexts[VCF_oXSTRAND], &xstrand, 0);
    bool is_xstrand = (*xstrand == 'X');

    switch (last_ostatus) {
        case LO_OK_REF_SAME    : 
            if (is_xstrand) {
                // Note: REF and ALT are always of length 1 (or this line would have been rejected)
                RECONSTRUCT1 (REVCOMP[(int)ref[0]]);
                RECONSTRUCT1 ('\t');
                RECONSTRUCT1 (REVCOMP[(int)alt[0]]); 
            }
            else {            
                // Note: REF and ALT may be of any length
                RECONSTRUCT_SEP (ref, ref_len, '\t'); 
                RECONSTRUCT (alt, alt_len);  
            }
            break;

        case LO_OK_REF_ALT_SWTCH : 
            RECONSTRUCT1 (is_xstrand ? REVCOMP[(int)alt[0]] : alt[0]); 
            RECONSTRUCT1 ('\t');
            RECONSTRUCT1 (is_xstrand ? REVCOMP[(int)ref[0]] : ref[0]);
            break;
        
        default: 
            ASSPIZ (0, "oStatus=%s (coords=%u) unexpected for VCF_SPECIAL_other_REFALT", last_ostatus_name_piz, (unsigned)vb->last_index (VCF_COORDS));
    }

    return false; // no new value
}


// called for reconstructing the REF field of the INFO/LIFTOVER and INFO/LIFTBACK fields. It reconstructs the full REFALT and then discards the ALT.
SPECIAL_RECONSTRUCTOR (vcf_piz_special_LIFT_REF)
{
    DidIType refalt_did_i = flag.luft ? VCF_REFALT : VCF_oREFALT;

    snip = AFTERENT (char, vb->txt_data);
    reconstruct_from_ctx (vb, refalt_did_i, 0, true);
    snip_len = (unsigned)(AFTERENT (char, vb->txt_data) - snip);

    const char *after_ref = memchr (snip, '\t', snip_len);
    ASSPIZ (after_ref, "expected a \\t in the %s snip: \"%.*s\" ostatus=%s", 
            vb->contexts[refalt_did_i].name, snip_len, snip, last_ostatus_name_piz);

    vb->txt_data.len = ENTNUM (vb->txt_data, after_ref);

    return false; // no new value
}

