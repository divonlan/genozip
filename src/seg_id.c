// ------------------------------------------------------------------
//   seg_id.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "context.h"
#include "seg.h"

static inline ContextP id_fallback_ctx (VBlockP vb, ContextP ctx)
{
    return ctx_get_ctx (vb, sub_dict_id_(ctx->dict_id, '0'));
}

static inline ContextP id_num1_ctx (VBlockP vb, ContextP ctx)
{
    return ctx_get_ctx (vb, sub_dict_id_(ctx->dict_id, '1'));
}

static inline ContextP id_num2_ctx (VBlockP vb, ContextP ctx)
{
    return ctx_get_ctx (vb, sub_dict_id_(ctx->dict_id, '2'));
}

static void seg_id_add_to_unknown (ContextP ctx, STRp(id))
{
    int id_i = -1;
    
    int i=0;  for (;i < NUM_COLLECTED_WORDS && segconf.unk_ids_tag_name[i][0]; i++)
        if (!memcmp (ctx->tag_name, segconf.unk_ids_tag_name[i], MAX_TAG_LEN)) {
            id_i = i;
            break;
        }

    if (i == NUM_COLLECTED_WORDS) return; // too many unknown ID contexts in file, no room to store

    if (id_i == -1) {
        id_i = i;
        memcpy (segconf.unk_ids_tag_name[i], ctx->tag_name, MAX_TAG_LEN);
    }

    for (i=0; i < NUM_UNK_ID_CTXS ; i++)
        if (!segconf.unk_ids[id_i][i][0]) { // we still have room
            memcpy (segconf.unk_ids[id_i][i], id, MIN_(id_len, UNK_ID_LEN));
            break;
        }
}

// Get type of ID (valid for this VB) and initialize contexts. Examples of types:
// "rs17030902" (SNP id) - IDT_ALPHA_NUM
// "ENSMUSE00000866652.1" (exon id) - IDT_ALPHA_NUM_DOT_INT
// "34324" - (integer) - treated as IDT_ALPHA_INT / IDT_ALPHA_NUM 
// "something" - treated as IDT_OTHER

static void seg_id_field_init (VBlockP vb, ContextP ctx, STRp(id), bool hint_zero_padded_fixed_len)
{
    SAFE_NULT(id);

    Did st_did_id = ctx->st_did_i == DID_NONE ? ctx->did_i : ctx->st_did_i;

    int prefix_len = strcspn (id, "0123456789");

    if (prefix_len == id_len) 
        goto idt_other;

    bool zero_padded = hint_zero_padded_fixed_len || (id[prefix_len] == '0');

    int digits_len = strspn (id + prefix_len, "0123456789");
    
    int an_len = prefix_len + digits_len;

    if (an_len == id_len) {
        ctx->id_type = zero_padded ? IDT_ALPHA_NUM : IDT_ALPHA_INT;
        ctx->ltype = LT_DYN_INT;

        ctx_consolidate_stats (vb, st_did_id, id_fallback_ctx(vb, ctx)->did_i, DID_EOL);
    }

    else if (id_len >= an_len + 2 && id[an_len] == '.' && str_is_int (&id[an_len+1], id_len - (an_len+1))) {
        ctx->id_type = zero_padded ? IDT_ALPHA_NUM_DOT_INT : IDT_ALPHA_INT_DOT_INT;

        ContextP ctx_num1 = id_num1_ctx (vb, ctx);
        ContextP ctx_num2 = id_num2_ctx (vb, ctx);

        ctx_set_ltype (vb, LT_DYN_INT, ctx_num1->did_i, ctx_num2->did_i, DID_EOL);
        ctx_consolidate_stats (vb, st_did_id, ctx_num1->did_i, ctx_num2->did_i, id_fallback_ctx(vb, ctx)->did_i, DID_EOL);
    }

    else idt_other: {
        ctx->id_type = IDT_OTHER;
        ctx->no_stons = true;
    }

    ctx->is_initialized = true;
    SAFE_RESTORE;
}

// index in string of number at end of sting (returns s_len if there is no such number)
static inline int terminal_number_index (STRp(s))
{
    int i, limit = MAX_(0, (int)s_len-18); // max 18 digits in int64_t 
    
    for (i=s_len - 1; i >= limit && IS_DIGIT (s[i]); i--); 

    return i+1;
}

static bool seg_id_split (STRp(id),  
                          bool has_num2,
                          bool allow_num1_leading_zeros, // if true, leading zeros are part of num. If flase - part of alpha
                          pSTRp(alpha), int64_t *num, unsigned *num_len, int64_t *num2, unsigned *num2_len) // out
{
    int num_i = terminal_number_index (STRa(id));

    // if num2 is requested: it must be an integer, no leading zeros, and prefixed with a '.'
    if (has_num2) {
        if (num_i < id_len && num_i >= 2 && id[num_i-1] == '.' && id[num_i] != '0') {
            *num2_len = id_len - num_i;
            str_get_int (&id[num_i], id_len - num_i, num2);

            id_len = num_i - 1;
            num_i = terminal_number_index (STRa(id)); // get num1
        }
        else 
            return false; // no second number
    }

    if (num_i == id_len) return false; // no number

    if (!allow_num1_leading_zeros && (num_i + 1 < id_len)/*not just a single zero*/) 
        while (id[num_i] == '0') num_i++; // exclude leading zeros - they will be part of alpha

    *num_len = id_len - num_i;
    if (!str_get_int_dec (&id[num_i], id_len - num_i, (uint64_t*)num))
        return false;

    *alpha = id;
    *alpha_len = num_i; // note: 0 is a valid length

    return true;
}                   

// this version compresses better if the numeric part is expected to be variable-width without leading zeros
void seg_id_field (VBlockP vb, ContextP ctx, STRp(id),
                   bool hint_zero_padded_fixed_len, // numeric part is expected to be zero padded (or not). if hint is correct, this will result in better compression
                   unsigned add_bytes)     
{
    #define T(x) (ctx->id_type == IDT_##x)

    if (!ctx->is_initialized) 
        seg_id_field_init (vb, ctx, STRa(id), hint_zero_padded_fixed_len);
        
    if (T(OTHER)) {
        seg_add_to_local_text (vb, ctx, STRa(id), LOOKUP_NONE, add_bytes);

        if (segconf.running && id_len > 1)
            seg_id_add_to_unknown (ctx, STRa(id));
    }

    else {
        STR(alpha); 
        int64_t num1, num2;
        unsigned num1_len, num2_len=0;
        ContextP ctx_num1 = ctx;

        if (!seg_id_split (STRa(id),  
                           T(ALPHA_INT_DOT_INT) || T(ALPHA_NUM_DOT_INT), // has_num2
                           T(ALPHA_NUM) || T(ALPHA_NUM_DOT_INT),        // allow_num1_leading_zeros
                           pSTRa(alpha), &num1, &num1_len, &num2, &num2_len))
            goto fallback;
            
        // case: eg "ENSMUSE00000866652.1": container goes into ctx, and id's components into other contexts
        if (T(ALPHA_INT_DOT_INT) || T(ALPHA_NUM_DOT_INT)) {
            ctx_num1 = id_num1_ctx (vb, ctx);
            ContextP ctx_num2 = id_num2_ctx (vb, ctx);
            
            SmallContainer con = {
                .repeats   = 1,
                .nitems_lo = 2,
                .items = { { .dict_id = ctx_num1->dict_id, .separator[0] = '.' },  // alpha and num1
                           { .dict_id = ctx_num2->dict_id                      } } // num2
            };

            container_seg (vb, ctx, (ContainerP)&con, 0, 0, 0);

            seg_add_to_local_resizable (vb, ctx_num2, num2, num2_len);
        }

        seg_add_to_local_resizable (vb, ctx_num1, num1, 0);
        if (ctx->flags.store == STORE_INT) ctx_set_last_value (vb, ctx, num1);

        // case: not expecting leading zeros - if there are any, they will be part of alpha
        if (T(ALPHA_INT) || T(ALPHA_INT_DOT_INT)) {
            SAFE_ASSIGN (&id[-1], SNIP_LOOKUP); // we assign it anyway bc of the macro convenience, but we included it only if num_digits>0
            seg_by_ctx (VB, id-1, 1 + alpha_len, ctx_num1, alpha_len + num1_len); // account for the entire length, and sometimes with \t
            SAFE_RESTORE;
        }

        // case: expecting leading zeros - seg as SNIP_NUMERIC
        else { 
            SAFE_ASSIGNx(&id[-3], SNIP_NUMERIC,      1);
            SAFE_ASSIGNx(&id[-2], '0'/*LT_DYN_INT*/, 2);
            SAFE_ASSIGNx(&id[-1], '0' + num1_len,    3);
            
            seg_by_ctx (vb, id-3, 3 + alpha_len, ctx_num1, alpha_len + num1_len); // SNIP_NUMERIC with prefix

            SAFE_RESTOREx(1); SAFE_RESTOREx(2); SAFE_RESTOREx(3);
        }

        ctx->txt_len += add_bytes - num1_len - num2_len - alpha_len;
    }

    return;

fallback: { 
    // in case of id type mismatch, add id to local, but of a different context, as ctx->local is DYN_INT
    ContextP ctx_fallback = id_fallback_ctx(vb, ctx);
    seg_by_ctx (vb, STRa(id), ctx_fallback, id_len);

    STRl(snip, 30);
    seg_prepare_snip_other (SNIP_REDIRECTION, ctx_fallback->dict_id, false, 0, snip);
    seg_by_ctx (vb, STRa(snip), ctx, add_bytes - id_len);
}
    #undef T
}

// optimized for variable-length integers without leading zeros (but its ok if they have leading zeros)
bool seg_id_field_varlen_int_cb (VBlockP vb, ContextP ctx, STRp(id), uint32_t repeat)
{
    seg_id_field (vb, ctx, STRa(id), false, id_len);
    return true; // segged successfully
}

// optimized for fixed-length zero-padded integers (but its ok if they are not fixed-width)
bool seg_id_field_fixed_int_cb (VBlockP vb, ContextP ctx, STRp(id), uint32_t repeat)
{
    seg_id_field (vb, ctx, STRa(id), true, id_len);
    return true; // segged successfully
}
