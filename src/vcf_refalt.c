// ------------------------------------------------------------------
//   vcf_refalt.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

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
#include "reference.h"
#include "ref_iupacs.h"
#include "dict_id_gen.h"

// ---------
// Seg stuff
// ---------

// optimize REF and ALT, for simple one-character REF/ALT (i.e. mostly a SNP or no-variant)
static inline void vcf_refalt_seg_ref_alt_snp (VBlockVCFP vb, char ref, char alt)
{
    char new_ref=0, new_alt=0;
    decl_acgt_decode;

    // if we have a reference, we use it 
    // except: if --match-chrom, we assume the user just wants to match, and we don't burden him with needing the reference to decompress
    if (((IS_REF_EXTERNAL && !flag.match_chrom_to_reference) || IS_REF_EXT_STORE)) {
        PosType32 pos = vb->last_int(VCF_POS);

        RefLock lock = REFLOCK_NONE;

        Range *range = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, 1, WORD_INDEX_NONE, 
                                          (IS_REF_EXT_STORE ? &lock : NULL));
        if (range) { // this chrom is in the reference
            uint32_t index_within_range = pos - range->first_pos;

            ref_assert_nucleotide_available (range, pos);
            
            // note: in GVCF, REF='N' for ~5% of human genome, but we can't identify as our reference doesn't store N
            if (ref == REF (index_within_range)) 
                new_ref = '-'; // normally, we expect our REF to match the reference...

            if (IS_REF_EXT_STORE)
                bits_set (&range->is_set, index_within_range);

            ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE
        }
    }

    switch (ref) {
        #define NEW_ALT(p0,p1,p2,p3) ((alt=='.'&&CTX(FORMAT_RGQ)->line_has_RGQ)||(alt==p0&&!CTX(FORMAT_RGQ)->line_has_RGQ))?'+' : alt==p0?'3' : alt==p1?'0' : alt==p2?'1' : alt==p3?'2' : 0
        //                           +/3  0   1   2
        case 'G' : new_alt = NEW_ALT('A','C','T','G'); break; // G to: A=239681 C=46244 T=44084 (based on counting simple SNPs from from chr22 of 1000 genome project phase 3)
        case 'C' : new_alt = NEW_ALT('T','G','A','C'); break; // C to: T=238728 G=46508 A=43685
        case 'A' : new_alt = NEW_ALT('G','C','T','A'); break; // A to: G=111967 C=30006 T=26335
        case 'T' : new_alt = NEW_ALT('C','G','A','T'); break; // T to: C=111539 G=29504 A=25599
        default  : break;
    }

    // if anything was done, we create a "special" snip
    if (new_ref || new_alt) {
        char refalt_special[4] = { SNIP_SPECIAL, VCF_SPECIAL_main_REFALT };
        refalt_special[2] = new_ref ? new_ref : ref;
        refalt_special[3] = new_alt ? new_alt : alt;

        seg_by_did (VB, refalt_special, sizeof(refalt_special), VCF_REFALT, 0); // we do the accounting in vcf_refalt_seg_main_ref_alt
    }
    // if not - just the normal snip
    else {
        char refalt_normal[3] = { 0, '\t', 0 };
        refalt_normal[0] = ref;
        refalt_normal[2] = alt;

        seg_by_did (VB, refalt_normal, sizeof(refalt_normal), VCF_REFALT, 0); // we do the account in vcf_refalt_seg_main_ref_alt
    }
}

// seg a deletion against the reference, if available
static inline bool vcf_refalt_seg_del_against_reference (VBlockVCFP vb, STRp(ref), char alt)
{
    decl_acgt_decode;

    // if we have a reference, we use it except: if --match-chrom, we assume the user just wants to match, and we don't burden him with needing the reference to decompress
    if ((!(IS_REF_EXTERNAL && !flag.match_chrom_to_reference) && !IS_REF_EXT_STORE) || 
        ref[0] != alt) return false;

    PosType32 pos = vb->last_int(VCF_POS);

    RefLock lock = REFLOCK_NONE;

    Range *range = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_len, WORD_INDEX_NONE, 
                                        (IS_REF_EXT_STORE ? &lock : NULL));
    
    if (!range || pos < range->first_pos || pos + ref_len - 1 > range->last_pos)
        return false; // region implied by REF doesn't fully exist in the reference
    
    uint32_t index_within_range = pos - range->first_pos;

    for (int i=0; i < ref_len; i++)
        if (ref[i] != REF (index_within_range + i)) {
            ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE
            return false; // REF doesn't match reference
        }

    SNIPi2 (SNIP_SPECIAL, VCF_SPECIAL_main_REFALT_DEL, ref_len);
    seg_by_did (VB, STRa(snip), VCF_REFALT, 0);

    if (IS_REF_EXT_STORE)
        bits_set_region (&range->is_set, index_within_range, ref_len);

    ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE

    return true;
}


#define IS_SINGLE_BASE_ALTS(n_alts, alt_len) ((n_alts)*2-1 == (alt_len))

void vcf_refalt_seg_main_ref_alt (VBlockVCFP vb, STRp(ref), STRp(alt))
{
    // optimize ref/alt in the common case of single-character
    if (vb->n_alts == 1 && vb->var_types[0] == VT_SNP) 
        vcf_refalt_seg_ref_alt_snp (vb, *ref, *alt);

    else if (vb->n_alts == 1 && vb->var_types[0] == VT_DEL && vcf_refalt_seg_del_against_reference (vb, STRa(ref), alt[0]))
        {}

    else {
        char ref_alt[ref_len + 1 + alt_len];
        memcpy (ref_alt, ref, ref_len);
        ref_alt[ref_len] = '\t';
        memcpy (&ref_alt[ref_len+1], alt, alt_len);

        seg_by_did (VB, ref_alt, ref_len + alt_len + 1, VCF_REFALT, 0);
    }
        
    CTX(VCF_REFALT)->txt_len += ref_len + alt_len + 2;
}

// ---------
// PIZ stuff
// ---------

// ZIP + PIZ
void vcf_parse_main_alt (VBlockVCFP vb)
{
    if (vb->n_alts != 0) return; // already parsed

    vb->n_alts = str_split_do (STRa(vb->main_alt), ARRAY_LEN(vb->alts), ',', vb->alts, vb->alt_lens, false, NULL);

    if (vb->n_alts == 0) vb->n_alts = -1; // failed

    else 
        for (int alt_i=0; alt_i < vb->n_alts; alt_i++) {       
            // missing due upstream deletion
            if (str_is_1chari (vb->alt, alt_i, '*'))
                vb->var_types[alt_i] = VT_MISSING;

            // symbolic alts e.g. "<DEL>"
            else if (!str_is_ACGTN (STRi(vb->alt, alt_i)))
                vb->var_types[alt_i] = VT_UNKNOWN; 
            
            // e.g. A C
            else if (vb->main_ref_len == 1 && vb->alt_lens[alt_i] == 1) 
                vb->var_types[alt_i] = VT_SNP;

            // e.g. ACG A
            else if (vb->main_ref_len > 1 && vb->alt_lens[alt_i] == 1 && *vb->main_ref == *vb->alts[alt_i])
                vb->var_types[alt_i] = VT_DEL;

            // e.g. A ACG
            else if (vb->main_ref_len == 1 && vb->alt_lens[alt_i] > 1)
                vb->var_types[alt_i] = VT_INS;

            // e.g. ACGCG ACG
            else if (vb->main_ref_len > vb->alt_lens[alt_i] && !memcmp (vb->main_ref, vb->alts[alt_i], vb->alt_lens[alt_i]))
                vb->var_types[alt_i] = VT_DEL_LONG;
                
            // e.g. ACG ACGCG
            else if (vb->main_ref_len < vb->alt_lens[alt_i] && !memcmp (vb->main_ref, vb->alts[alt_i], vb->main_ref_len))
                vb->var_types[alt_i] = VT_INS_LONG;

            // e.g. ACG TCG
            else if (vb->main_ref_len == vb->alt_lens[alt_i] && !memcmp (vb->main_ref+1, vb->alts[alt_i]+1, vb->main_ref_len-1))
                vb->var_types[alt_i] = VT_SNP_LONG;
            
            // e.g. ATATGTG ATG - left anchored, part of remaining bases after deletion are substituted 
            else if (vb->main_ref_len > vb->alt_lens[alt_i] && *vb->main_ref == *vb->alts[alt_i])
                vb->var_types[alt_i] = VT_SUBST_DEL;

            // e.g. ATG ATATGTG - left anchored, part of ref is substituted, and then insertion 
            else if (vb->main_ref_len < vb->alt_lens[alt_i] && *vb->main_ref == *vb->alts[alt_i])
                vb->var_types[alt_i] = VT_SUBST_INS;

            else  // To do - identify additional types 
                vb->var_types[alt_i] = VT_UNKNOWN; // need to reset as could be set by previous line
        }
}

// item callback of REFALT in TOPLEVEL, called with files compressed starting 14.0.12
// for files 14 and older, called from vcf_piz_filter
void vcf_piz_refalt_parse (VBlockVCFP vb)
{
    decl_ctx (VCF_REFALT);

    vb->main_ref     = last_txtx (vb, ctx);
    vb->main_alt     = memchr (vb->main_ref, '\t', ctx->last_txt.len) + 1;
    vb->main_ref_len = vb->main_alt - vb->main_ref - 1;
    vb->main_alt_len = ctx->last_txt.len - vb->main_ref_len - 1;

    vcf_parse_main_alt (vb);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_obsolete_dvcf)
{
    ABORT0 ("Error: use Genozip up to 15.0.41");
    return NO_NEW_VALUE;
}

// Sometimes called to reconstruct the "main" refalt (main AT THE TIME OF SEGGING), to reconstruct SNPs that were stored relative to a reference.
// This SPECIAL is only used for lines that are bi-allelic SNPs, and either REF or ALT match reference and user compressed with --reference. 
// Not used if compressed with --chain.
SPECIAL_RECONSTRUCTOR (vcf_piz_special_main_REFALT)
{
    if (!reconstruct) goto done;

    ASSPIZ (snip_len==2, "expecting snip_len=2 but seeing %u", snip_len);

    // snip is 3 characters - REF, \t, ALT
    char ref_alt[3] = { 0, '\t', 0 };
    char ref_value = 0;
    
    if (snip[0] == '-' || snip[1] == '-') { 
        PosType32 pos = CTX (VCF_POS)->last_value.i;

        ConstRangeP range = ref_piz_get_range (vb, gref, HARD_FAIL);
        
        uint32_t idx = pos - range->first_pos;

        ASSPIZ (!flag.debug || IS_REF_EXTERNAL || ref_is_nucleotide_set (range, idx), 
                "reference is not set: chrom=%.*s pos=%d", STRf(range->chrom_name), pos);

        decl_acgt_decode;
        ref_value = REF (idx);
    }

    // recover ref
    ref_alt[0] = (snip[0] == '-') ? ref_value : snip[0];

    // if --drop-genotypes, we don't normally conusme VCF_SAMPLES, so we do it here
    if (flag.drop_genotypes && segconf.has[FORMAT_RGQ]) {
        ASSPIZ0 (!vb->frozen_state.prm8[0], "Peek mode node supported - we're consuming");
        LOAD_SNIP (VCF_SAMPLES);
    }
    
    // recover ALT
    char r=ref_alt[0];
    switch (snip[1]) {
        case '+' : ref_alt[2] = vcf_piz_line_has_RGQ (VB_VCF)?'.' : r=='A'?'G' : r=='C'?'T' : r=='G'?'A' : 'C'; break; // the most common value
        case '0' : ref_alt[2] = (r=='C' || r=='T') ? 'G' : 'C';             break; // added v14
        case '1' : ref_alt[2] = (r=='C' || r=='T') ? 'A' : 'T';             break; // added v14
        case '3' : ref_alt[2] = r=='A'?'G' : r=='C'?'T' : r=='G'?'A' : 'C'; break; // !RGQ, added v14
        case '2' : ref_alt[2] = r;                                          break; // added v14
        case '-' : ref_alt[2] = ref_value;                                  break; // existed up to v13
        default  : ref_alt[2] = snip[1];                 
    }

    RECONSTRUCT (ref_alt, sizeof (ref_alt));

done:
    return NO_NEW_VALUE;
}   

// Sometimes called to reconstruct the "main" refalt (main AT THE TIME OF SEGGING), 
// to reconstruct deletions that were stored relative to a reference.
SPECIAL_RECONSTRUCTOR (vcf_piz_special_main_REFALT_DEL)
{
    if (!reconstruct) goto done;

    int ref_len = atoi (snip);
     
    PosType32 pos = CTX (VCF_POS)->last_value.i;

    ConstRangeP range = ref_piz_get_range (vb, gref, HARD_FAIL);
        
    uint32_t idx = pos - range->first_pos;

    decl_acgt_decode;
    char *next = BAFTtxt;
    for (int i=0; i < ref_len; i++)
        *next++ = REF (idx + i);

    *next++ = '\t';
    *next = REF (idx); // ALT is the anchor base

    Ltxt += ref_len + 2;
    
done:
    return NO_NEW_VALUE;
}   

// used by FORMAT/PS - reconstructs REF or ALT depending on the parameter - '0' REF or '1'..MAX_ALLELES-1 - ALT
SPECIAL_RECONSTRUCTOR (vcf_piz_special_COPY_REForALT)
{
    if (!reconstruct) return false;

    STRlast (refalt, VCF_REFALT);
    
    // shortcut in case of SNP
    if (refalt_len == 3)
        RECONSTRUCT1 (refalt[*snip=='0' ? 0 : 2]);

    else {
        str_split (refalt, refalt_len, 2, '\t', item, true);
        ASSPIZ0 (n_items, "failed splitting REFALT");

        if (*snip=='0') // REF
            RECONSTRUCT (items[0], item_lens[0]);

        else {
            str_split (items[1], item_lens[1], MAX_ALLELES-1, ',', alt, false);

            int allele = (snip_len == 1) ? (*snip - '0')
                                         : (snip[0] - '0') * 10 + (snip[1] - '0');
            
            RECONSTRUCT (alts[allele-1], alt_lens[allele-1]);
        }
    }

    return NO_NEW_VALUE;
}

// --snps-only implementation (called from vcf_piz_container_cb)
bool vcf_refalt_piz_is_variant_snp (VBlockVCFP vb)
{
    // true if ALL alts are SNPs
    for (int alt_i=0; alt_i < vb->n_alts; alt_i++)
        if (vb->var_types[alt_i] != VT_SNP) return false;

    return true;
}

// --indels-only implementation (called from vcf_piz_container_cb)
bool vcf_refalt_piz_is_variant_indel (VBlockVCFP vb)
{
    #define VT(x) (vb->var_types[alt_i] == VT_##x)

    // true is ANY alt is an indel
    for (int alt_i=0; alt_i < vb->n_alts; alt_i++)
        if (VT(DEL) || VT(INS) || VT(DEL_LONG) || VT(INS_LONG) || VT(SUBST_DEL) || VT(SUBST_INS)) 
            return true;

    return false;
}
