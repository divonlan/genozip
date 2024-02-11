// ------------------------------------------------------------------
//   vcf_ultima.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "reconstruct.h"
#include "reference.h"
#include "seg.h"
#include "context.h"

#define VT(x) (vb->var_types[alt_i] == VT_##x)

void vcf_ultima_seg_initialize (VBlockVCFP vb)
{
    ctx_set_store (VB, STORE_INT, INFO_ASSEMBLED_HAPS, INFO_FILTERED_HAPS, DID_EOL);

    ctx_set_dyn_int (VB, INFO_X_LM, INFO_X_RM, DID_EOL);
}

//----------------------
// INFO/X_LM & INFO/X_RM
//----------------------

// not 100% but close enough
static PosType64 X_LM_RM_predict_pos (VBlockVCFP vb, ConstRangeP range, PosType64 line_pos, bool is_rm, unsigned lm_seq_len, int alt_i)
{
    bool is_del = (VT(DEL) || VT(DEL_LONG) || VT(SUBST_DEL));
    bool is_ins = (VT(INS) || VT(INS_LONG) || VT(SUBST_INS)); 
    
    if (is_del || is_ins) {
        int n_reps = 0;
        if (is_rm) {
            rom payload     = (is_del ? vb->main_ref     : vb->alts[alt_i]    ) + 1;
            int payload_len = (is_del ? vb->main_ref_len : vb->alt_lens[alt_i]) - 1;

            PosType64 pos = line_pos + 1;
            int payload_i=0;
            decl_acgt_decode;
            while (REFp(pos) == payload[payload_i]) {
                pos++; 
                payload_i=(payload_i + 1) % payload_len;
                if (!payload_i) n_reps++; // counts full repeats
            }
            
            if ((is_del && !str_is_monochar (STRa(payload))) || (VT(INS) && !n_reps))
                pos -= payload_i; // exclude partial payload matches
            
            else if (is_ins && n_reps && payload_i) 
                pos = (pos - payload_i) + 1;

            return pos;
        }
        else
            return line_pos - (lm_seq_len - 1);
    }

    else if (VT(SNP)) {
        if (is_rm) 
            return line_pos + 1;
        else
            return line_pos - lm_seq_len;
    }

    return 0; // fail
} 

void vcf_seg_INFO_X_LM_RM (VBlockVCFP vb, ContextP ctx, STRp(seq))
{
    decl_acgt_decode;
    RefLock lock = REFLOCK_NONE;

    if (0) fallback: { // cannot be inside a variable-length array block
        ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE
        seg_by_ctx (VB, STRa(seq), ctx, seq_len);
        return;
    }

    if (!flag.reference) 
        goto fallback;

    str_split (seq, seq_len, 0, ',', seq, false);
    if (n_seqs != vb->n_alts) goto fallback;

    if (!str_is_ACGT (STRi(seq, 0), NULL)) goto fallback; // reference doesn't support N or IUPACs

    // all components are predicted to be identical (this is not true for multi-alt variants with for example, INS+SNP or DEL+SNP_LONG)
    for (int alt_i=1; alt_i < n_seqs; alt_i++)
        if (!str_issame_(seqs[0], seq_lens[0], seqs[alt_i], seq_lens[alt_i]))
            goto fallback;

    PosType64 line_pos = DATA_LINE(vb->line_i)->pos;

    RangeP range = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), line_pos, seq_lens[0], WORD_INDEX_NONE, 
                                      (IS_REF_EXT_STORE ? &lock : NULL));
    
    if (!range || line_pos < range->first_pos || line_pos + seq_lens[0] - 1 > range->last_pos)
        goto fallback;

    bool is_rm = (ctx->did_i == INFO_X_RM);
    PosType64 pos = X_LM_RM_predict_pos (vb, range, line_pos, is_rm, seq_lens[0], (vb->var_types[0] == VT_MISSING) ? 1 : 0);

    // check again, with pos
    if (!pos || pos < range->first_pos || pos + seq_lens[0] - 1 > range->last_pos)
        goto fallback;

    uint32_t index_within_range = pos - range->first_pos;

    // find delta
    int delta = 999;
    for (int delta_i=0; delta_i < 200; delta_i++) {
        int candidate_delta = DEINTERLACE(typeof(delta_i), delta_i); // 0, -1, 1, -2, 2...
        
        bool found = true; // optimistic
        for (int i=0; i < seq_lens[0]; i++)
            if (seq[i] != REFp (pos + i + candidate_delta)) {
                found = false;
                break;
            }
        
        if (found) {
            delta = candidate_delta;
            break;
        }
    }
    if (delta == 999) 
        goto fallback;
    
    SNIPi2_2 (SNIP_SPECIAL, VCF_SPECIAL_X_LM_RM, seq_lens[0], delta);
    seg_by_ctx (VB, STRa(snip), ctx, seq_len);

    if (IS_REF_EXT_STORE)
        bits_set_region (&range->is_set, index_within_range, seq_lens[0]);

    ref_unlock (gref, &lock); 
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_X_LM_RM)
{    
    if (!reconstruct) return NO_NEW_VALUE;

    str_split_ints (snip, snip_len, 2, ',', item, true);
    int seq_len = items[0];
    int delta   = items[1];

    ConstRangeP range = ref_piz_get_range (vb, gref, HARD_FAIL);

    PosType64 pos = delta + X_LM_RM_predict_pos (VB_VCF, range, CTX(VCF_POS)->last_value.i, (ctx->did_i == INFO_X_RM), seq_len, (VB_VCF->var_types[0] == VT_MISSING) ? 1 : 0);
    
    char *next = BAFTtxt;
    char *start = next;
    decl_acgt_decode;
    
    for (uint32_t i=0; i < seq_len; i++)
        *next++ = REFp (pos++);

    Ltxt += seq_len;

    // elements beyond the first
    for (int alt_i=1; alt_i < VB_VCF->n_alts; alt_i++) {
        RECONSTRUCT1 (',');
        RECONSTRUCT (start, seq_len);
    }

    return NO_NEW_VALUE;
}

//----------
// INFO/X_IC
//----------

static rom X_IC_prediction (VBlockVCFP vb, int alt_i)
{
    #define DELorNA(vt) ((vt)==VT_DEL || (vt)==VT_DEL_LONG) ? (vt) : VT_SNP;

    VariantType vt = vb->var_types[alt_i];
    if (vt == VT_MISSING) {
        if (alt_i == 0 && vb->n_alts > 1)   
            vt = DELorNA (vb->var_types[1]);
        if (alt_i > 0) 
            vt = DELorNA (vb->var_types[alt_i-1]);
    }
    
    switch (vt) {
        case VT_DEL_LONG  :
        case VT_SUBST_DEL :
        case VT_DEL       : return "del";
        case VT_SUBST_INS : 
        case VT_INS_LONG  :
        case VT_INS       : return "ins";
        default           : return "NA";
    }
}

void vcf_seg_INFO_X_IC (VBlockVCFP vb, ContextP ctx, STRp(ic_str))
{
    bool predicted = true; // optimistic

    str_split (ic_str, ic_str_len, 0, ',', ic, false);

    if (n_ics != vb->n_alts)
        predicted = false;

    else for (int alt_i=0; alt_i < n_ics; alt_i++) {
        rom prediction = X_IC_prediction (vb, alt_i);

        if (!str_issame_(ics[alt_i], ic_lens[alt_i], prediction, strlen (prediction))) {
            predicted = false;
            break;
        }
    }

    if (predicted)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_X_IC }, 2, ctx, ic_str_len);

    else
        seg_by_ctx (VB, STRa(ic_str), ctx, ic_str_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_X_IC)
{    
    if (reconstruct) 
        for (int alt_i=0; alt_i < VB_VCF->n_alts; alt_i++) {
            if (alt_i) RECONSTRUCT1 (',');

            rom prediction = X_IC_prediction (VB_VCF, alt_i);
            RECONSTRUCT (prediction, strlen (prediction));
        }

    return NO_NEW_VALUE;
}

//----------
// INFO/X_IL
//----------

// <ID=X_IL,Number=A,Type=Integer,Description="Flow: length of indel">
static int X_IL_prediction (VBlockVCFP vb, int alt_i)
{
    switch (vb->var_types[alt_i]) {
        case VT_DEL_LONG :
        case VT_SUBST_DEL:
        case VT_DEL      : return vb->main_ref_len - vb->alt_lens[alt_i];
        case VT_INS_LONG :
        case VT_SUBST_INS:
        case VT_INS      : return vb->alt_lens[alt_i] - vb->main_ref_len;
        default          : return 0;  // '.'
    }
}

void vcf_seg_INFO_X_IL (VBlockVCFP vb, ContextP ctx, STRp(il_str))
{
    bool predicted = true; // optimistic

    str_split (il_str, il_str_len, 0, ',', il, false);

    if (n_ils != vb->n_alts)
        predicted = false;

    else for (int alt_i=0; alt_i < n_ils; alt_i++) {
        int64_t il;

        if (str_is_1chari (il, alt_i, '.'))
            predicted &= (X_IL_prediction (vb, alt_i) == 0);

        else if (str_get_int (STRi(il, alt_i), &il))
            predicted &= (X_IL_prediction (vb, alt_i) == il);
    }

    if (predicted)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_X_IL }, 2, ctx, il_str_len);

    else
        seg_by_ctx (VB, STRa(il_str), ctx, il_str_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_X_IL)
{    
    if (reconstruct) 
        for (int alt_i=0; alt_i < VB_VCF->n_alts; alt_i++) {
            if (alt_i) RECONSTRUCT1 (',');

            int prediction = X_IL_prediction (VB_VCF, alt_i);

            if (prediction) RECONSTRUCT_INT (prediction);
            else            RECONSTRUCT1 ('.');
        }

    return NO_NEW_VALUE;
}

//-----------
// INFO/X_HIN
//-----------

// <ID=X_HIN,Number=A,Type=String,Description="Flow: nucleotide of the hmer indel, if so">
static char X_HIN_prediction (VBlockVCFP vb, int alt_i)
{
    STRlast (x_hil_str, INFO_X_HIL);
    str_split_ints (x_hil_str, x_hil_str_len, vb->n_alts, ',', x_hil, true);
    if (!n_x_hils) return 0; // invalid HIL -> not predictable
    
    if (!x_hils[alt_i])
        return '.';

    else switch (vb->var_types[alt_i]) {
        case VT_DEL_LONG :
        case VT_SUBST_DEL:
        case VT_DEL      : return vb->main_ref[vb->main_ref_len-1];
        case VT_INS_LONG :
        case VT_SUBST_INS:
        case VT_INS      : return vb->alts[alt_i][vb->alt_lens[alt_i]-1];
        default          : return 0;  
    }
}

void vcf_seg_INFO_X_HIN (VBlockVCFP vb, ContextP ctx, STRp(hin_str))
{
    str_split (hin_str, hin_str_len, 0, ',', hin, false);

    bool predicted = (n_hins == vb->n_alts) && ctx_encountered_in_line (VB, INFO_X_HIL); 

    if (predicted)
        for (int alt_i=0; alt_i < n_hins; alt_i++) 
            predicted &= hin_lens[alt_i] == 1 &&
                         *hins[alt_i] == X_HIN_prediction (vb, alt_i);

    if (predicted)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_X_HIN }, 2, ctx, hin_str_len);

    else
        seg_by_ctx (VB, STRa(hin_str), ctx, hin_str_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_X_HIN)
{    
    if (reconstruct) 
        for (int alt_i=0; alt_i < VB_VCF->n_alts; alt_i++) {
            if (alt_i) RECONSTRUCT1 (',');

            RECONSTRUCT1 (X_HIN_prediction (VB_VCF, alt_i));
        }

    return NO_NEW_VALUE;
}

//-----------
// INFO/X_HIL
//-----------

static int X_HIL_hmer_length (VBlockVCFP vb, char hmer)
{
    decl_acgt_decode;
    int hmer_len = 0;

    if (IS_ZIP) {    
        PosType64 pos = DATA_LINE(vb->line_i)->pos + 1; // hmer starts after anchor base

        RefLock lock = REFLOCK_NONE;

        #define MAX_HMER_LEN 64 // arbitrary, but cannot be changed due to backcomp

        RangeP range = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, MAX_HMER_LEN, WORD_INDEX_NONE, 
                                          (IS_REF_EXT_STORE ? &lock : NULL));

        if (range && pos >= range->first_pos && pos + MAX_HMER_LEN <= range->last_pos) { // range usable
            while (hmer_len < MAX_HMER_LEN && REFp(pos + hmer_len) == hmer)
                hmer_len++;

            if (IS_REF_EXT_STORE) 
                bits_set_region (&range->is_set, (pos - range->first_pos), MAX_HMER_LEN);
        }

        ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE
    }

    else {
        PosType64 pos = CTX(VCF_POS)->last_value.i + 1; // hmer starts after anchor base

        ConstRangeP range = ref_piz_get_range (VB, gref, HARD_FAIL);

        while (hmer_len < MAX_HMER_LEN && REFp(pos + hmer_len) == hmer)
            hmer_len++;
    }

    return hmer_len;
}

// <ID=X_HIL,Number=A,Type=Integer,Description="Flow: length of the hmer indel, if so">
static int X_HIL_prediction (VBlockVCFP vb, int alt_i, bool use_reference)
{
    bool is_non_h_indel = ctx_encountered_in_line (VB, INFO_VARIANT_TYPE) &&
                          str_issame_(STRlst(INFO_VARIANT_TYPE), _S("non-h-indel"));

    if (!VT(DEL) && !VT(INS)) return 0;

    rom payload     = (VT(DEL) ? vb->main_ref     : vb->alts[alt_i]    ) + 1;
    int payload_len = (VT(DEL) ? vb->main_ref_len : vb->alt_lens[alt_i]) - 1;

    if (is_non_h_indel) // non-h-indel can be determined by the payload
        return str_is_monochar (STRa(payload)) ? payload_len : 0;

    else if (!use_reference || !str_is_monochar (STRa(payload))) 
        return 0;

    else
        return X_HIL_hmer_length (vb, payload[0]) + (VT(INS) ? payload_len : 0);
}

void vcf_seg_INFO_X_HIL (VBlockVCFP vb, ContextP ctx, STRp(hil_str))
{
    str_split_ints (hil_str, hil_str_len, vb->n_alts, ',', hil, true);

    bool use_reference = !!flag.reference;
    bool predicted = (n_hils == vb->n_alts);

    if (predicted)
        for (int alt_i=0; alt_i < n_hils; alt_i++) 
            predicted &= (hils[alt_i] == X_HIL_prediction (vb, alt_i, use_reference));

    if (predicted)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_X_HIL, '0' + use_reference }, 3, ctx, hil_str_len);

    else 
        seg_by_ctx (VB, STRa(hil_str), ctx, hil_str_len);

    seg_set_last_txt (VB, ctx, STRa(hil_str));
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_X_HIL)
{    
    if (reconstruct) {
        bool use_reference = snip[0] - '0';

        for (int alt_i=0; alt_i < VB_VCF->n_alts; alt_i++) {
            if (alt_i) RECONSTRUCT1 (',');

            RECONSTRUCT_INT (X_HIL_prediction (VB_VCF, alt_i, use_reference));
        }
    }

    return NO_NEW_VALUE;
}

//------------------
// INFO/VARIANT_TYPE
//------------------

// if we have the reference, we can destinguish better between h-indel and non-h-indel
static bool VARIANT_TYPE_ref_confirms_hmer (VBlockVCFP vb, STRp(seq), bool use_reference)
{
    decl_acgt_decode;

    if (!use_reference) return true;

    if (IS_ZIP) {    
        bool confirm = true; // confirmed unless proven incorrect
        PosType64 pos = DATA_LINE(vb->line_i)->pos;

        RefLock lock = REFLOCK_NONE;

        RangeP range = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, seq_len, WORD_INDEX_NONE, 
                                        (IS_REF_EXT_STORE ? &lock : NULL));

        if (range && pos >= range->first_pos && pos + vb->main_ref_len <= range->last_pos && // range usable
            REFp(pos + vb->main_ref_len) != seq[seq_len-1]) // homopolymer does not continue one more base
            confirm = false;

        if (IS_REF_EXT_STORE)
            bits_set (&range->is_set, pos + vb->main_ref_len - range->first_pos);

        ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE
        return confirm; // actually, we can't use the reference
    }

    else {
        ConstRangeP range = ref_piz_get_range (VB, gref, HARD_FAIL);

        return REFp(CTX(VCF_POS)->last_value.i + vb->main_ref_len) == seq[seq_len-1];
    }
}

// <ID=VARIANT_TYPE,Number=1,Type=String,Description="Flow: type of variant: SNP/NON-H-INDEL/H-INDEL">
static void VARIANT_TYPE_prediction (VBlockVCFP vb, pSTRp(prediction), bool use_reference)
{
    int count_snps = 0;
    int count_hmers = 0;

    for (int alt_i=0; alt_i < vb->n_alts; alt_i++)    
        if (VT(SNP) || VT(MISSING)) 
            count_snps++;

        else if (VT(DEL) && str_is_monochar (vb->main_ref+1, vb->main_ref_len-1) &&
                 VARIANT_TYPE_ref_confirms_hmer (vb, STRa(vb->main_ref), use_reference))
            count_hmers++;
            
        else if (VT(INS) && str_is_monochar (vb->alts[alt_i]+1, vb->alt_lens[alt_i]-1) &&
                 VARIANT_TYPE_ref_confirms_hmer (vb, STRi(vb->alt, alt_i), use_reference))
            count_hmers++;
    
    if (count_snps == vb->n_alts)                    { *prediction = "snp";         *prediction_len = STRLEN("snp");         }
    else if (count_hmers + count_snps == vb->n_alts) { *prediction = "h-indel";     *prediction_len = STRLEN("h-indel");     }
    else                                             { *prediction = "non-h-indel"; *prediction_len = STRLEN("non-h-indel"); } 
}

void vcf_seg_INFO_VARIANT_TYPE (VBlockVCFP vb, ContextP ctx, STRp(vt))
{
    STR(prediction);
    VARIANT_TYPE_prediction (vb, pSTRa(prediction), flag.reference);

    if (str_issame (vt, prediction))
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_VARIANT_TYPE, '0' + !!flag.reference }, 3, ctx, vt_len);

    else
        seg_by_ctx (VB, STRa(vt), ctx, vt_len);

    seg_set_last_txt (VB, ctx, STRa(vt));
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_VARIANT_TYPE)
{    
    if (reconstruct) {
        bool use_reference = snip[0] - '0';

        STR(prediction);
        VARIANT_TYPE_prediction (VB_VCF, pSTRa(prediction), use_reference);

        RECONSTRUCT_str (prediction);
    }

    return NO_NEW_VALUE;
}

//-------------------
// INFO/FILTERED_HAPS
//-------------------

void vcf_seg_INFO_FILTERED_HAPS (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    if (ctx_has_value_in_line_(vb, CTX(INFO_ASSEMBLED_HAPS))) 
        seg_delta_vs_other (VB, ctx, CTX(INFO_ASSEMBLED_HAPS), STRa(value));

    else
        seg_integer_or_not (VB, ctx, STRa(value), value_len);
}
