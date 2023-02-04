// ------------------------------------------------------------------
//   recon_history.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "reconstruct.h"
#include "vblock.h"
#include "file.h"
#include "piz.h"

//------------------------------------------
// Storing current line's value in "history"
//------------------------------------------

// upon dropping a line, we copy its textual history from txt_data to per_line
void reconstruct_copy_dropped_line_history (VBlockP vb)
{
    HistoryWord *hw;

    for_ctx 
        if ((ctx->flags.store_per_line || ctx->flags.spl_custom) && // note: we can't test history.len as the buffer might be used for other purposes if not storing history
            ctx->flags.store != STORE_INT && ctx->flags.store != STORE_INDEX &&  // store textual
            ctx->history.len32 && // 0 if context is skipped
            (hw = B(HistoryWord, ctx->history, vb->line_i))->lookup == LookupTxtData) {

            uint32_t per_line_start = ctx->per_line.len32;
            if (ctx->last_txt.len && hw->len) {
                buf_add_more (vb, &ctx->per_line, STRtxtw(*hw), "per_line");
                BNXTc (ctx->per_line) = 0; // nul-terminate
            }

            hw->lookup = LookupPerLine;
            hw->index  = per_line_start;
        }
}

// store last_value in context history - for copying by buddy line
void reconstruct_store_history (VBlockP vb, ContextP ctx)
{
    switch (ctx->flags.store) {
        case STORE_INT:   
            *B(int64_t, ctx->history, vb->line_i) = ctx->last_value.i; 
            break;

        case STORE_INDEX: 
            *B(WordIndex, ctx->history, vb->line_i) = ctx->last_value.i; 
            break;
        
        default:          
            // store as LookupTxtData, if not already stored by reconstruct_store_history_rollback_recon()
            if (B(HistoryWord, ctx->history, vb->line_i)->lookup != LookupPerLine) 
                *B(HistoryWord, ctx->history, vb->line_i) = (HistoryWord)
                    { .lookup       = LookupTxtData, 
                      .index        = ctx->last_txt.index, 
                      .len          = ctx->last_txt.len,
                      .ctx_specific = ctx->history.prm8[0] };
    }
}

// store in history (in per_line), and then rollback reconstruction
void reconstruct_store_history_rollback_recon (VBlockP vb, ContextP ctx, rom recon_start)
{
    HistoryWord *hw = B(HistoryWord, ctx->history, vb->line_i);
    *hw = (HistoryWord){ .index = ctx->per_line.len32, 
                         .len = (BAFTtxt - recon_start), 
                         .lookup = LookupPerLine };

    buf_add_more (VB, &ctx->per_line, recon_start, hw->len, "per_line");
    BNXTc (ctx->per_line) = 0;         // nul-terminate
    vb->txt_data.len32 = BNUMtxt(recon_start); // remove reconstructed text from txt_data
}

// consume and store in history the next value in the context without reconstructing it to txt_data
void reconstruct_to_history (VBlockP vb, ContextP ctx)
{
    STR(snip);
    HistoryWord *hw = B(HistoryWord, ctx->history, vb->line_i);

    // case: MC:Z is in dict (referred to from a b250)
    if (ctx->b250.len ||
        (!ctx->b250.len && !ctx->local.len && ctx->dict.len)) {  // all_the_same case - no b250 or local, but have dict      
        WordIndex wi = LOAD_SNIP(ctx->did_i); // note: if we have no b250, local but have dict, this will be word_index=0 (see ctx_get_next_snip)

        if (snip_len==1 && *snip == SNIP_LOOKUP) 
            goto snip_is_in_local;
        
        else if (!snip_len || IS_PRINTABLE(snip[0]))
            *hw = (HistoryWord){ .index = wi, .len = snip_len, .lookup = LookupDict };

        // not a textual snip (eg SNIP_SPECIAL) - reconstruct and then copy
        else { 
            rom txt = BAFTtxt;
            reconstruct_one_snip (vb, ctx, wi, STRa(snip), RECON_ON);
            
            *hw = (HistoryWord){ .index = ctx->per_line.len32, .len = (BAFTtxt - txt), .lookup = LookupPerLine };

            buf_add_more (vb, &ctx->per_line, txt, hw->len, "per_line");
            BNXTc (ctx->per_line) = 0;         // nul-terminate
            vb->txt_data.len32 = BNUMtxt(txt); // remove reconstructed text from txt_data
        } 
    }
    
    // case: snip is in local
    else snip_is_in_local: {
        uint32_t char_index = LOAD_SNIP_FROM_LOCAL (ctx);
        *hw = (HistoryWord){ .index = char_index, .len = snip_len, .lookup = LookupLocal }; 
    }
}

rom lookup_type_name (LookupType lookup)
{
    static rom names[] = LOOKUP_TYPE_NAMES;

    return (lookup < 0 || lookup >= ARRAY_LEN(names)) ? "Invalid" : names[lookup];
}
