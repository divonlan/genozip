// ------------------------------------------------------------------
//   recon_history.c
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

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
        case STORE_INT:   *B(int64_t,     ctx->history, vb->line_i) = ctx->last_value.i; break;

        case STORE_INDEX: *B(WordIndex,   ctx->history, vb->line_i) = ctx->last_value.i; break;
        
        default:          *B(HistoryWord, ctx->history, vb->line_i) = (HistoryWord){
                                .lookup       = LookupTxtData, 
                                .index        = ctx->last_txt.index, 
                                .len          = ctx->last_txt.len,
                                .ctx_specific = ctx->history.prm8[0] 
                          };
    }
}

//---------------------------------------------
// Consuming data stored by a "historical" line
//---------------------------------------------

void reconstruct_set_buddy (VBlockP vb)
{
    ASSPIZ0 (vb->buddy_line_i == NO_LINE, "Buddy line already set for the current line");

    ContextP buddy_ctx = ECTX (_SAM_BUDDY); // all data types using buddy are expected to have a dict_id identical to _SAM_BUDDY (but different did_i)

    int32_t num_lines_back = reconstruct_from_local_int (vb, buddy_ctx, 0, false);

    // a bug that existed 12.0.41-13.0.1 (see bug 367): we stored buddy in machine endianty instead of BGEN32.
    // When we reach this point pizzing in a buggy file, num_lines_back will be in BGEN since piz_adjust_one_local set it. 
    // We detect buggy files by local.prm8[0]=0 (see sam_seg_initialize) and convert it back to machine endianity.
    if (!buddy_ctx->local.prm8[0])    
        num_lines_back = BGEN32 ((uint32_t)num_lines_back);

    vb->buddy_line_i = vb->line_i - num_lines_back; // convert value passed (distance in lines to buddy) to 0-based buddy_line_i

    ASSPIZ (vb->buddy_line_i != NO_LINE, "Expecting vb->buddy_line_i=%d to be non-negative. num_lines_back=%d buddy_ctx->local.prm8[0]=%d", 
            vb->buddy_line_i, num_lines_back, buddy_ctx->local.prm8[0]);
}

void reconstruct_from_buddy_get_textual_snip (VBlockP vb, ContextP ctx, pSTRp(snip))
{
    ASSPIZ0 (vb->buddy_line_i != NO_LINE, "No buddy line is set for the current line");
    ASSPIZ (ctx->history.len32, "ctx->history not allocated for ctx=%s, perhaps seg_initialize did't set store_per_line?", ctx->tag_name);

    HistoryWord word = *B(HistoryWord, ctx->history, vb->buddy_line_i);
    BufferP buf=NULL; 
    CharIndex char_index = word.index;

    switch (word.lookup) {
        case LookupTxtData : buf = &vb->txt_data  ; break;
        case LookupDict    : buf = &ctx->dict; 
                             char_index = B(CtxWord, ctx->word_list, word.index)->char_index; 
                             break;
        case LookupLocal   : buf = &ctx->local    ; break;
        case LookupPerLine : buf = &ctx->per_line ; break;
        default : ASSPIZ (false, "Invalid value word.lookup=%d", word.lookup);
    }

    ASSPIZ (char_index < buf->len, "buddy word ctx=%s buddy_line_i=%d char_index=%"PRIu64" is out of range of buffer %s len=%"PRIu64, 
            ctx->tag_name, vb->buddy_line_i, char_index, buf->name, buf->len);

    *snip = Bc (*buf, char_index);
    *snip_len = word.len;
}

// Copy from buddy: buddy is data that appears on a specific "buddy line", in this context or another one. Not all lines need
// Note of difference vs. lookback: with buddy, not all lines need to have the data (eg MC:Z), so the line number is constant,
// but if we had have used lookback, the lookback value would have been different between different fields.  
HasNewValue reconstruct_from_buddy (VBlockP vb, ContextP ctx, STRp(snip), bool reconstruct, ValueType *new_value)
{
    ContextP base_ctx = ctx;

    // set buddy if needed, and not already set (in BAM it is already set in sam_piz_filter).
    if (vb->buddy_line_i == NO_LINE) 
        reconstruct_set_buddy (vb); 

    // optional: base context is different than ctx
    if (snip_len > 1) {
        snip--; snip_len++; // reconstruct_get_other_ctx_from_snip skips the first char
        base_ctx = reconstruct_get_other_ctx_from_snip (vb, ctx, pSTRa(snip));
    }

    ASSPIZ (vb->buddy_line_i < base_ctx->history.len32, "history not set for %s, perhaps seg forgot to set store_per_line? (vb->buddy_line_i=%d history.len=%u)", 
            base_ctx->tag_name, vb->buddy_line_i, base_ctx->history.len32);

    // case: numeric value 
    if (ctx->flags.store == STORE_INT) {
        new_value->i = *B(int64_t, base_ctx->history, vb->buddy_line_i);
        if (reconstruct) RECONSTRUCT_INT (new_value->i);
        return HAS_NEW_VALUE;
    }

    // case: word index (use if we always store word_index. to mix word_index with other method, use LookupDict)
    else if (ctx->flags.store == STORE_INDEX) {
        new_value->i = *B(WordIndex, base_ctx->history, vb->buddy_line_i);

        if (reconstruct) {
            ctx_get_snip_by_word_index (ctx, new_value->i, snip);
            RECONSTRUCT_snip;
        }
        return HAS_NEW_VALUE;
    }

    // case: textual value
    else {
        if (reconstruct) {
            reconstruct_from_buddy_get_textual_snip (vb, base_ctx, pSTRa(snip));
            RECONSTRUCT_snip;
        }

        return NO_NEW_VALUE; 
    }
}

