// ------------------------------------------------------------------
//   vcf_format_GQ.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "base64.h"
#include "zip_dyn_int.h"

void vcf_segconf_finalize_GQ (VBlockVCFP vb)
{
    if ((segconf.count_GQ_by_GP > vb->lines.len * vcf_num_samples / 5) && (segconf.count_GQ_by_GP > segconf.count_GQ_by_PL))
        segconf.FMT_GQ_method = BY_GP;

    else if ((segconf.count_GQ_by_PL > vb->lines.len * vcf_num_samples / 5) && (segconf.count_GQ_by_PL >= segconf.count_GQ_by_GP))
        segconf.FMT_GQ_method = BY_PL;

    else if (segconf.has[FORMAT_DP])
        segconf.FMT_GQ_method = MUX_DOSAGExDP;

    else if (segconf.vcf_is_giggle)
        segconf.FMT_GQ_method = GQ_INTEGER;

    else
        segconf.FMT_GQ_method = MUX_DOSAGE;
}

static SORTER (value_sorter)
{
    return DESCENDING_RAW (*(int64_t *)a, *(int64_t*)b); // sort in reverse order - usually faster as GP[0] / PL[0] are usually the biggest (corresponding to GT=0/0)
}

static inline void get_value (VBlockVCFP vb, Did did_i, pSTRp(value))
{
    if (IS_ZIP)
        SETlast (*value, did_i);
    else
        reconstruct_peek (VB, CTX(did_i), STRa(value));
}

static int64_t vcf_predict_GQ_by_PL (VBlockVCFP vb)
{
    STR(pl);
    get_value (vb, FORMAT_PL, pSTRa(pl));

    if (IS_PERIOD(pl)) return -1;

    // get sorted array of PL values
    str_split_ints (pl, pl_len, 500, ',', val, false);
    if (!n_vals) return 0; // array too long or not all integers

    qsort (vals, n_vals, sizeof(int64_t), value_sorter);

    // usually, we take the second-lowest value
    int64_t prediction = n_vals == 1  ? vals[0]
                       : vals[0] == 0 ? -1             // all values are 0 - predicting a .
                       : n_vals >= 4 && vals[0] == 1 && vals[1] == 0 ? -1 // observed empircally
                       :                vals[n_vals-2]; // 2nd lowest

    return MIN_(prediction, 99); // cap at 99
}

static int64_t vcf_predict_GQ_by_GP (VBlockVCFP vb)
{
    STR(gp);
    get_value (vb, FORMAT_GP, pSTRa(gp));

    if (IS_PERIOD(gp)) return -1;

    // get array of GP values
    str_split_floats (gp, gp_len, 64, ',', val, false, 0);
    if (!n_vals) return 0; // array too long or not all float

    // round to integers (in-place)
    int64_t *val_ints = (int64_t *)vals; 
    for (int i=0; i < n_vals; i++) 
        val_ints[i] = (int64_t)(vals[i] + 0.5); 

    // take mid-value
    qsort (val_ints, n_vals, sizeof(int64_t), value_sorter);
    return val_ints[n_vals / 2]; 
}

void vcf_seg_FORMAT_GQ (VBlockVCFP vb)
{
    decl_ctx(FORMAT_GQ);
    STRlast (gq, FORMAT_GQ);

    int64_t gq_value;
    int64_t prediction;
    
    if (IS_PERIOD(gq)) gq_value = -1;

    else if (!str_get_int_range64 (STRa(gq), 0, 1000000, &gq_value)) fallback: {
        seg_by_ctx (VB, STRa(gq), ctx, gq_len); // simple snip and not seg_integer_or_not to not prevent singletons    
        return;
    }

    if (segconf_running) {
        if (ctx_encountered (VB, FORMAT_GP) && ABS (vcf_predict_GQ_by_GP (vb) - gq_value) <= 5) segconf.count_GQ_by_GP++;
        if (ctx_encountered (VB, FORMAT_PL) && ABS (vcf_predict_GQ_by_PL (vb) - gq_value) <= 5) segconf.count_GQ_by_PL++;
    }

    else switch (segconf.FMT_GQ_method) {
        case BY_GP: 
            if (!ctx_encountered (VB, FORMAT_GP)) goto fallback;
            prediction = vcf_predict_GQ_by_GP (vb);
            goto do_seg;

        case BY_PL: {
            if (!ctx_encountered (VB, FORMAT_PL)) goto fallback;
            prediction = vcf_predict_GQ_by_PL (vb);
        
        do_seg: {
            SNIPi2 (SNIP_SPECIAL, VCF_SPECIAL_GQ, prediction - gq_value);
            seg_by_ctx (VB, STRa(snip), ctx, gq_len);
            break;
        }}

        case MUX_DOSAGExDP:
            vcf_seg_FORMAT_mux_by_dosagexDP (vb, ctx, STRa(gq), &vb->mux_GQ);
            break;

        case MUX_DOSAGE:
            vcf_seg_FORMAT_mux_by_dosage (vb, ctx, STRa(gq), (DosageMultiplexer *)&vb->mux_GQ);
            break;

        case GQ_INTEGER:
            if (IS_PERIOD(gq)) 
                dyn_int_append_nothing_char (VB, ctx, gq_len);
            else
                dyn_int_append (VB, ctx, gq_value, gq_len);
            break;
            
        default:
            ABORT ("Invalid value segconf.FMT_GQ_method=%u", segconf.FMT_GQ_method);
    }
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_GQ_old); // forward

SPECIAL_RECONSTRUCTOR (vcf_piz_special_GQ)
{
    if (segconf.FMT_GQ_method == GQ_old) // up to 15.0.36
        return vcf_piz_special_GQ_old (vb, ctx, STRa(snip), new_value, reconstruct);

    ASSPIZ (segconf.FMT_GQ_method == BY_GP || segconf.FMT_GQ_method == BY_PL, "Invalid FMT_GQ_method=%u", segconf.FMT_GQ_method);

    int64_t prediction = segconf.FMT_GQ_method == BY_GP ? vcf_predict_GQ_by_GP (VB_VCF)
                                                        : vcf_predict_GQ_by_PL (VB_VCF);
    int64_t delta = atoi (snip);
    new_value->i = prediction - delta;

    if (new_value->i >= 0) {
        if (reconstruct) RECONSTRUCT_INT (new_value->i);
        return HAS_NEW_VALUE;
    }

    else {
        if (reconstruct) RECONSTRUCT1 ('.');
        return NO_NEW_VALUE;
    }    
}

// ----------------------------------
// old - used for files up to 15.0.36
// ----------------------------------

static int64_t vcf_piz_predict_GQ_old (VBlockVCFP vb, Did src_did_i)
{
    bool is_gp = (src_did_i == FORMAT_GP);
    STR(src);

    reconstruct_peek (VB, CTX(src_did_i), pSTRa(src));

    str_split (src, src_len, 30, ',', item, false);

    int64_t values[n_items];
    unsigned n_values=0;

    for (int i=0; i < n_items; i++) {
        if (item_lens[i]==1 && items[i][0]=='.')
            values[n_values++] = 0; // consider a '.' to be 0

        else if (is_gp) {
            double f;
            // until 15.0.59 str_get_float returned false for values with mantissa/exponent format, 
            // but now it returns return true, so we explicitly fail the condition for this format
            if (!memchr(items[i], 'e', item_lens[i]) &&
                !memchr(items[i], 'E', item_lens[i]) &&
                str_get_float (STRi(item,i), &f, 0, 0)) 
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

SPECIAL_RECONSTRUCTOR (vcf_piz_special_GQ_old)
{
    rom tab = memchr (snip, '\t', snip_len);

    ContextP src_ctx = SCTX0 (snip);

    int64_t prediction = vcf_piz_predict_GQ_old (VB_VCF, src_ctx->did_i);
    int64_t delta = atoi (tab+1);

    new_value->i = prediction - delta;
    RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}
