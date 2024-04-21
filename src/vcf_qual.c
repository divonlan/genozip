// ------------------------------------------------------------------
//   vcf_qual.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "piz.h"
#include "reconstruct.h"

void vcf_segconf_finalize_QUAL (VBlockVCFP vb)
{
    // gnomAD, dbSNP - store in local
    if (segconf.vcf_is_gnomad || segconf.vcf_is_dbSNP)
        segconf.vcf_QUAL_method = VCF_QUAL_local;

    // GATK GVCF - multiplex by has_RGQ
    else if (segconf.has[FORMAT_RGQ] || segconf.vcf_is_gatk_gvcf)
        segconf.vcf_QUAL_method = VCF_QUAL_by_RGQ;
    
    // Dragen - QUAL is predicted to equal for item of GP (TODO: may be applicable beyond dragen and single sample)
    else if (segconf.vcf_is_dragen && vcf_num_samples == 1) 
        segconf.vcf_QUAL_method = VCF_QUAL_by_GP;

    // SvABA - copy from mate, or store in local
    else if (segconf.vcf_is_svaba || segconf.vcf_is_manta) // note: in pbsv, QUAL="."
        segconf.vcf_QUAL_method = VCF_QUAL_mated;

    else
        segconf.vcf_QUAL_method = VCF_QUAL_DEFAULT;
}

static inline void str_float_remove_trailing_zeros (qSTRp(str))
{
    if (!memchr (str, '.', *str_len)) return; // cannot remove trailing zeros from an integer

    while (str[*str_len-1] == '0') (*str_len)--;
    if (str[*str_len-1] == '.') (*str_len)--;
}

// changes the number of decimals in a floating point number, where if number of decimals is decreased,
// number if rounded up or down to the nearest lesser-granularity float. 
// This function is MUCH faster than converting the string to double and then rounding.
static StrText str_round (STRp(in), 
                          bool boundary_round_down, // if excess decimals are 50000... round down if true, up if false
                          int out_decimals,  // number of digits after the period required 
                          bool truncate_trailing_zeros,
                          uint32_t *out_len)
{
    ASSERT (in_len < 20, "in_len=%u too big", in_len);
    ASSERT (out_decimals <= 7, "out_decimals=%u too big", out_decimals);
    
    StrText out = {};
    uint32_t in_decimals;

    if (str_is_simple_float (STRa(in), &in_decimals)) {

        // case: no rounding is needed, because out granularity is not less than in granularity
        if (in_decimals <= out_decimals) {
            memcpy (out.s, in, in_len);
            *out_len = in_len;

            // zero pad if needed
            if (in_decimals < out_decimals && !truncate_trailing_zeros) {
                if (!in_decimals) out.s[in_len] = '.';
                memset (&out.s[in_len + !in_decimals], '0', out_decimals - in_decimals); 
                *out_len += (out_decimals - in_decimals) + !in_decimals;
            }
        }

        // case: out needs less decimals than in - rounding might be needed
        else {
            // copy the integer and as many decimals as requested (if any) from in to out
            *out_len = (in_len - in_decimals - !out_decimals) + out_decimals;
            memcpy (out.s, in, *out_len);

            // "redundant" decimal digits 
            rom red = (in[*out_len] == '.') ? &in[*out_len + 1] : &in[*out_len]; // first redundant digit
            int red_len = &in[in_len] - red;

            // case where we need to round up
            if (*red >= '6' || // note: we know in[out_len] is a digit, because in's granualrity is larger that out's  
                (*red == '5' && (!boundary_round_down || (red_len > 1 && !str_is_monochar_(red+1, red_len-1, '0'))))) 
                
                for (int i=*out_len-1; i >= 0; i--) {
                    if (out.s[i] == '.') continue;
                    else if (out.s[i] <= '8') {
                        out.s[i]++;
                        break;
                    }
                    else { // '9'
                        out.s[i] = '0';
                        if (i == 0) {
                            memmove (out.s + 1, out.s, *out_len);
                            out.s[0] = '1';
                            (*out_len)++;
                        }
                    }
                }

            // remove trailing zeros
            if (out_decimals && truncate_trailing_zeros) 
                str_float_remove_trailing_zeros (out.s, out_len);
        }
    }

    // not a simple float: either an exponent/mantisa float, or not a float at all.
    else {
        double in_value;
        SAFE_NULT (in); // MUCH faster atof (at least in mingw)
        in_value = atof (in); // float string terminated by comma (or colon or newline) (add epsilon so that sprintf rounds up - eg .005 to .01)
        SAFE_RESTORE; 

        *out_len = snprintf (out.s, sizeof (out.s), "%.*f", out_decimals, in_value + (boundary_round_down ? -0.000000001 : 0.000000001));

        // remove trailing zeros
        if (out_decimals && truncate_trailing_zeros) 
            str_float_remove_trailing_zeros (out.s, out_len);
    }

    return out;
}

static inline void vcf_seg_QUAL_by_GP (VBlockVCFP vb, ContextP ctx, STRp(qual))
{
    // case: we have GP in this line
    if (ctx_encountered_in_line (VB, FORMAT_GP)) {
        STR(gp1_str);
        str_item_i (STRlst(FORMAT_GP), ',', 0, pSTRa (gp1_str)); // always returns true for item_i=0

        // test first with rounding boundary cases down
        uint32_t prediction_len;
        rom prediction = str_round (STRa(gp1_str), true, segconf.vcf_QUAL_decimals, segconf.vcf_QUAL_truncate_trailing_zeros, &prediction_len).s;

        if (str_issame (qual, prediction)) {
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_QUAL_BY_GP, 'D', '0' + segconf.vcf_QUAL_decimals, '0' + segconf.vcf_QUAL_truncate_trailing_zeros }, 5, ctx, qual_len+1);
            return;
        }

        // test first with rounding boundary cases up
        prediction = str_round (STRa(gp1_str), false, segconf.vcf_QUAL_decimals, segconf.vcf_QUAL_truncate_trailing_zeros, &prediction_len).s;

        if (str_issame (qual, prediction)) {
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_QUAL_BY_GP, 'U', '0' + segconf.vcf_QUAL_decimals, '0' + segconf.vcf_QUAL_truncate_trailing_zeros }, 5, ctx, qual_len+1);
            return;
        }
    }

    // case: no GP in this line
    else {
        if (str_is_1char (qual, '.')) {
            // note: decimals set as usual. They are not used for reconstructing '.' - just set in hope of "all the same"
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_QUAL_BY_GP, 'D', '0' + segconf.vcf_QUAL_decimals, '0' + segconf.vcf_QUAL_truncate_trailing_zeros }, 5, ctx, qual_len+1);
            return;
        }
    }

    // fallback
    seg_by_ctx (VB, STRa(qual), ctx, qual_len+1);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_QUAL_BY_GP)
{
    vcf_piz_defer_to_after_samples (QUAL);
    ctx->QUAL.by_GP        = true;
    ctx->QUAL.is_rounddown = (snip[0] == 'D');
    ctx->QUAL.decimals     = snip[1] - '0';
    ctx->QUAL.truncate_trailing_zeros = snip[2]- '0';

    return NO_NEW_VALUE;
}

// called from toplevel callback
void vcf_piz_insert_QUAL_by_GP (VBlockVCFP vb)
{
    decl_ctx (VCF_QUAL);

    if (ctx_encountered_in_line (VB, FORMAT_GP)) {
        STR(gp1_str);
        str_item_i (STRlst(FORMAT_GP), ',', 0, pSTRa (gp1_str)); // always returns true for item_i=0

        uint32_t qual_len;
        rom qual = str_round (STRa(gp1_str), ctx->QUAL.is_rounddown, ctx->QUAL.decimals, ctx->QUAL.truncate_trailing_zeros, &qual_len).s;

        if (IS_RECON_INSERTION(ctx)) 
            vcf_piz_insert_field (vb, ctx, STRa(qual), segconf.wid_QUAL.width);

        if (ctx->flags.store == STORE_FLOAT)
            // note: atof again, bc value might be slightly different than gp1 due to rounding in sprintf        
            // note: we explicitly build a ValueType union bc just passing a double doesn't work in clang
            ctx_set_last_value (VB, ctx, (ValueType){ .f = atof (qual) }); 
    }

    else {
        if (IS_RECON_INSERTION(ctx)) 
            vcf_piz_insert_field (vb, ctx, ".", 1, segconf.wid_QUAL.width);
    }
}

static inline void vcf_seg_QUAL_by_RGQ (VBlockVCFP vb, ContextP ctx, STRp(qual))
{
    bool has_rgq = CTX(FORMAT_RGQ)->line_has_RGQ;
    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, VCF_QUAL, (MultiplexerP)&vb->mux_QUAL, has_rgq);
    
    if (!has_rgq && vcf_num_samples > 1) // too many unique QUAL values when there are many samples
        seg_add_to_local_string (VB, channel_ctx, STRa(qual), LOOKUP_SIMPLE, qual_len+1);
    else
        seg_by_ctx (VB, STRa(qual), channel_ctx, qual_len+1);
    
    seg_by_ctx (VB, STRa(vb->mux_QUAL.snip), ctx, 0);
}

void vcf_seg_QUAL (VBlockVCFP vb, STRp(qual))
{
    START_TIMER;

    decl_ctx (VCF_QUAL);

    if (segconf.running) {
        SEGCONF_RECORD_WIDTH (QUAL, qual_len);

        if (!str_is_1char (qual, '.')) {
            rom period = memchr (qual, '.', qual_len); 
            uint8_t decimals = period ? qual_len - (period - qual)/*integer*/ - 1/*period*/ : 0;
            if (decimals < segconf.vcf_QUAL_decimals) segconf.vcf_QUAL_truncate_trailing_zeros = true;
            if (decimals > segconf.vcf_QUAL_decimals) segconf.vcf_QUAL_decimals = decimals;
        }
    }

    switch (segconf.vcf_QUAL_method) {
        case VCF_QUAL_local :
            seg_add_to_local_string (VB, ctx, STRa(qual), LOOKUP_NONE, qual_len+1);
            break;

        case VCF_QUAL_by_RGQ: 
            vcf_seg_QUAL_by_RGQ (vb, ctx, STRa(qual));
            break;

        case VCF_QUAL_mated: 
            vcf_seg_sv_copy_mate (vb, ctx, STRa(qual), TW_QUAL, TW_QUAL, false, qual_len+1);
            break;

        case VCF_QUAL_by_GP: 
            vcf_seg_QUAL_by_GP (vb, ctx, STRa(qual));
            break;

        case VCF_QUAL_DEFAULT: 
            seg_by_ctx (VB, STRa(qual), ctx, qual_len+1);
            break;

        default:
            ABORT ("invalid QUAL method %u", segconf.vcf_QUAL_method);
    }

    COPY_TIMER (vcf_seg_QUAL);
}

