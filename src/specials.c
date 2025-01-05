// ------------------------------------------------------------------
//   specials.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "seg.h"
#include "piz.h"
#include "zip_dyn_int.h"
#include "base64.h"

// decode and store the the contexts in the first call for ctx (only one MINUS snip allowed per ctx)
static void decode_dicts (VBlockP vb, ContextP ctx, pSTRp(snip), int num_dicts/*-1 if unlimited*/)
{
#ifdef DEBUG
    ASSPIZ (num_dicts <= MAX_DICTS, "ctx=%s: %u dicts requested, but MAX_DICTS=%u", ctx->tag_name, num_dicts, MAX_DICTS);
#endif

    buf_alloc_zero (vb, &ctx->ctx_cache, 0, MAX_SNIP_DICTS, ContextP, 1, "ctx_cache");

    DictId dicts[MAX_SNIP_DICTS];
    uint32_t b64_len = *snip_len; // might end up being shorter 
    ctx->ctx_cache.len32 = base64_decode (*snip, &b64_len, (uint8_t *)dicts, num_dicts * sizeof(DictId)) / sizeof (DictId);

    ASSPIZ (num_dicts == -1 || num_dicts == ctx->ctx_cache.len32, "ctx=%s: expecting %u dicts in snip, but found %u",
            ctx->tag_name, num_dicts, ctx->ctx_cache.len32);

    for (int i=0; i < ctx->ctx_cache.len32; i++)
        *B(ContextP, ctx->ctx_cache, i) = ECTX (dicts[i]);

    // advance snip to after b64
    STRinc (*snip, b64_len);
}

//-------------------------------------------------------------------
// LEN_OF - field value is expected to be the length of another field
//-------------------------------------------------------------------

void seg_LEN_OF (VBlockP vb, ContextP ctx, STRp(len_str), uint32_t other_str_len, STRp(special_snip))
{
    int64_t len;
    if (!str_get_int (STRa(len_str), &len)) goto fallback;
    if (other_str_len != len) goto fallback;

    seg_by_ctx (VB, STRa(special_snip), ctx, len_str_len);
    return;

fallback:
    seg_by_ctx (VB, STRa(len_str), ctx, len_str_len);
}

SPECIAL_RECONSTRUCTOR (piz_special_LEN_OF)
{    
    ContextP base_ctx = reconstruct_special_get_base_ctx (VB, ctx, pSTRa(snip));

    uint32_t len;
    reconstruct_peek (vb, base_ctx, NULL, &len);
    
    if (reconstruct) {
        if (snip_len) RECONSTRUCT_snip; // optional prefix

        RECONSTRUCT_INT (len);
    }

    new_value->i = len;
    return HAS_NEW_VALUE;
}

//------------------------------------------------------------------------------------------
// ARRAY_LEN_OF - field value is expected to be the number of items in another field's array
//------------------------------------------------------------------------------------------

// if our value is eg "3" and can be predicted by the length of an array "xx,yy,zz"
void seg_by_ARRAY_LEN_OF (VBlockP vb, ContextP ctx, STRp(value), STRp(other_array), STRp(snip))
{
    int64_t n_items;
    char sep = snip[snip_len-1]; // array separator

    if (str_get_int (STRa(value), &n_items) &&
        n_items == 1 + str_count_char (STRa(other_array), sep))

        seg_by_ctx (vb, STRa(snip), ctx, value_len);

    else
        seg_integer_or_not (vb, ctx, STRa(value), value_len);
}

SPECIAL_RECONSTRUCTOR (piz_special_ARRAY_LEN_OF)
{    
    ContextP base_ctx = reconstruct_special_get_base_ctx (VB, ctx, pSTRa(snip));

    STR(other_value);
    reconstruct_peek (vb, base_ctx, pSTRa(other_value));

    char sep = *snip;
    new_value->i = 1 + str_count_char (STRa(other_value), sep);

    if (reconstruct) 
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// parameter is two or more dict_id's (in base64) reconstructs sum(dict[i].last_value)
SPECIAL_RECONSTRUCTOR (piz_special_PLUS)
{
    if (!ctx->ctx_cache.len32) 
        decode_dicts (vb, ctx, pSTRa(snip), -1);

    new_value->i = 0;

    // case: same_line set: Seg guarantees that values are in the line (before or after current context)
    if (ctx->flags.same_line) 
        for (int i=0; i < ctx->ctx_cache.len32; i++)
            new_value->i += reconstruct_peek (vb, (*B(ContextP, ctx->ctx_cache, i)), 0, 0).i;

    // case: same_line not set: add values that have been already encountered in this line
    else  
        for (int i=0; i < ctx->ctx_cache.len32; i++)
            if (ctx_has_value_in_line_(vb, *B(ContextP, ctx->ctx_cache, i)))
                new_value->i += (*B(ContextP, ctx->ctx_cache, i))->last_value.i;

    if (reconstruct)
        RECONSTRUCT_INT (new_value->i); 

    return HAS_NEW_VALUE; 
}

// parameter is two dict_id's (in base64). reconstructs dict1.last_value - dict2.last_value
SPECIAL_RECONSTRUCTOR (piz_special_MINUS)
{
    if (!ctx->ctx_cache.len32) 
        decode_dicts (vb, ctx, pSTRa(snip), 2);

    if (ctx->flags.same_line) 
        new_value->i = reconstruct_peek (vb, (*B(ContextP, ctx->ctx_cache, 0)), 0, 0).i - 
                       reconstruct_peek (vb, (*B(ContextP, ctx->ctx_cache, 1)), 0, 0).i;
                       
    else
        new_value->i = (*B(ContextP, ctx->ctx_cache, 0))->last_value.i - 
                       (*B(ContextP, ctx->ctx_cache, 1))->last_value.i;

    if (str_is_1char (snip, 'A')) 
        new_value->i = ABS (new_value->i); // ABS option since 15.0.48

    if (reconstruct)
        RECONSTRUCT_INT (new_value->i); 

    return HAS_NEW_VALUE; 
}

// divides value of other integer field by an integer constant
SPECIAL_RECONSTRUCTOR (piz_special_DIVIDE_BY)
{
    if (!ctx->ctx_cache.len32) { 
        decode_dicts (vb, ctx, pSTRa(snip), 1);
        
        ASSPIZ0 (str_get_int (STRa(snip), &ctx->ctx_cache.param) && ctx->ctx_cache.param > 0, 
                "Failed to get positive integer value");
    }

    if (ctx->flags.same_line) 
        new_value->i = reconstruct_peek (vb, (*B(ContextP, ctx->ctx_cache, 0)), 0, 0).i / ctx->ctx_cache.param;
                       
    else
        new_value->i = (*B(ContextP, ctx->ctx_cache, 0))->last_value.i / ctx->ctx_cache.param;

    if (reconstruct)
        RECONSTRUCT_INT (new_value->i); 

    return HAS_NEW_VALUE; 
}

//-------------------------------------------------------------------
// TEXTUAL_FLOAT - field is a textual base-10 float (e.g. VCF, SAM) 
//-------------------------------------------------------------------

// seg using a separate contexts for mantissa, n_fraction_digits and sign
// note that 0 (m=0,f=0,s=+), 0.0 (m=0,f=1,s=+), 0.000 (m=0,f=3,s=+) and -0.000 (m=0,f=3,s=-) are all destinct legal values
// scientific notation is not supported
// NOTE: only rarely (depending on the field), this is better than segging as a string - with floats that are unique, non-scietific and contain many digits (e.g. 6 significant digits)
void seg_textual_float (VBlockP vb, ContextP ctx, STRp(f), unsigned add_bytes)
{    
    bool negative = (f[0] == '-');
    if (negative) STRinc (f, 1);

    // format check: cannot have eg 00.13 or 012
    if (f[0] == '0' && f_len > 1 && f[1] != '.') goto bad_format; 

    // find decimal point and verify max one point and everything else is a digit
    int frac_digits = 0; 
    int64_t mantissa=0;
    for (int i=0; i < f_len; i++)
        if (IS_DIGIT(f[i]))
            mantissa = (mantissa * 10) + (f[i] - '0');
        
        else if (f[i] == '.') {
            if (frac_digits || i == f_len-1) goto bad_format; // bad format if more than one decimal point, or decimal point is last character
            frac_digits = f_len - i - 1;
        }

        else goto bad_format; // not digit or decimal point

    if (frac_digits > 255) goto bad_format;
    uint8_t frac_digits8 = frac_digits;

    ContextP frac_ctx = ctx_get_ctx (vb, sub_dict_id (ctx->dict_id, '0'));
    ContextP sign_ctx = ctx_get_ctx (vb, sub_dict_id (ctx->dict_id, '1'));
    
    if (!ctx->is_initialized) {
        ctx_set_ltype (VB, LT_UINT8, frac_ctx->did_i, sign_ctx->did_i, DID_EOL);

        int num_per_line = (ctx->did_i < MAX_NUM_PREDEFINED) ? segconf.local_per_line[ctx->did_i] : 1;        
        buf_alloc (vb, &frac_ctx->local, 0, vb->lines.len32 * num_per_line, uint8_t, 0, CTX_TAG_LOCAL); // initial allocation
        buf_alloc (vb, &sign_ctx->local, 0, vb->lines.len32 * num_per_line, uint8_t, 0, CTX_TAG_LOCAL);
        
        ctx_consolidate_stats (VB, ctx->did_i, frac_ctx->did_i, sign_ctx->did_i, DID_EOL);
        ctx->is_initialized = true;
    }

    seg_special0 (VB, VCF_SPECIAL_TEXTUAL_FLOAT, ctx, add_bytes);

    dyn_int_append (vb, ctx, mantissa, 0);

    seg_add_to_local_fixed (VB, frac_ctx, &frac_digits8, 1, LOOKUP_NONE, 0);
    seg_add_to_local_fixed (VB, sign_ctx, &negative, 1, LOOKUP_NONE, 0);

    return;

bad_format:
    if (negative) STRdec (f, 1); // restore
    seg_by_ctx (vb, STRa(f), ctx, add_bytes);
}

SPECIAL_RECONSTRUCTOR (piz_special_TEXTUAL_FLOAT)
{
    ContextP frac_ctx = ECTX (sub_dict_id (ctx->dict_id, '0'));
    ContextP sign_ctx = ECTX (sub_dict_id (ctx->dict_id, '1'));

    bool negative = NEXTLOCAL(uint8_t, sign_ctx);


    int64_t mantissa = reconstruct_from_local_int (vb, ctx, 0, RECON_OFF);
    uint8_t mant_digits = str_int_len (mantissa); 
    uint8_t frac_digits = NEXTLOCAL(uint8_t, frac_ctx);
    uint8_t int_digits  = mant_digits - frac_digits;
    
    if (reconstruct) {    
        if (negative) RECONSTRUCT1 ('-');

        // number starts with 0.[0]*
        if (frac_digits >= mant_digits) {    
            RECONSTRUCT ("0.", 2);
            
            for (int i=0; i < frac_digits - mant_digits; i++)
                RECONSTRUCT1 ('0');
        
            RECONSTRUCT_INT (mantissa);
        }
         
        else {
            char *int_start = BAFTtxt;
            RECONSTRUCT_INT (mantissa);

            if (frac_digits) {
                memmove (int_start + int_digits + 1, int_start + int_digits, frac_digits); // move fraction digits one up
                Ltxt++;

                *(int_start + int_digits) = '.';
            }
        }
    }

    return NO_NEW_VALUE;
}
