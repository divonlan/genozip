// ------------------------------------------------------------------
//   recon_peek.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "vblock.h"
#include "context.h"
#include "flags.h"
#include "piz.h"
#include "file.h"

// This is a sub-module of reconstruct, that provides the ability to peek contexts for their next item without changing
// their state, including if reconstructing the peeked contexts requires reconstructing or peeking other downstream contexts.
// By being able to do this, we gain the ability to seg one field based on another field that is still further down the line. 
//
// Example: Peek a context AA that is a container with one item BB and 2 repeats
// peek(AA)                    -> level=1
//   reconstruct(AA)
//     push (AA,level=1)      --> level=2, AA(1) on stack
//        reconstruct(BB)
//           push(BB,level=2) --> level=3, AA(1), BB(2) on stack
//           pop(BB)          --> level=2, no level=3 items to pop, so still AA(1), BB(2) on stack 
//                                and importantly, BB retains its new state (iterators etc)
//        reconstruct(BB)
//           push(BB,level=2) --> level=3, AA(1), BB(2), BB(2) on stack
//           pop(BB)          --> level=2, no level=3 items to pop, so still AA(1), BB(2), BB(2)
//     pop (AA)               --> level=1, pop level=2 items in reverse order (BB, BB), remaining stack: AA(1) 
//                            --> note: BB is now recovered to its state prior to the peek
// pop(done)                  --> level=0, pop level=1 items (AA), remaining stack: empty 
//                            --> note: AA is now recovered to its state prior to the peek 

#define stack_pointer (vb)->frozen_state.len32 // index of next item that will be pushed

#define RECON_STATE_SIZE 80 // needs adjustment if Context changes
typedef struct __attribute__ ((__packed__)) { 
    char ctx_p[sizeof (ContextP)]; // not pointer bc struct is not going to be word-aligned (except for first entry - assumed by ASSPIZ in reconstruct_peek)
    uint8_t level;
    char state[RECON_STATE_SIZE]; 
} Ice; // one item on the frozen stack

void recon_stack_initialize (void)
{
    ASSERT (RECON_STATE_SIZE == reconstruct_state_size_formula, "Expecting RECON_STATE_SIZE=%u == reconstruct_state_size_formula=%u - need to update RECON_STATE_SIZE in the code",
            RECON_STATE_SIZE, (int)reconstruct_state_size_formula);
}

static inline void increment_level (VBlockP vb, ContextP ctx)
{
    #define MAX_PEAK_STACK_LEVELS 50

    vb->peek_stack_level++; 

    ASSERT (vb->peek_stack_level <= MAX_PEAK_STACK_LEVELS, 
            "peek_stack_level exceeded MAX_PEAK_STACK_LEVELS=%u, while reconstructing %s. Dev note: Use --debug-peek to investigate infinite recursion.",
            MAX_PEAK_STACK_LEVELS, ctx->tag_name);
}

static inline void decrement_level (VBlockP vb, ContextP ctx)
{
    ASSERT (vb->peek_stack_level > 0, "Cannot decrement peek_stack_level because it is already 0 - while reconstructing %s", ctx->tag_name);

    vb->peek_stack_level--; 
}

// pops all downstream contexts used to reconstruct ctx
void recon_stack_pop (VBlockP vb, ContextP ctx, bool is_done_peek)
{
    ASSPIZ (stack_pointer, "stack is unexpectedly empty when popping ctx=%s. Dev tip: try --debug-peek", ctx->tag_name);  

    ARRAY (Ice, ice, vb->frozen_state);

    decrement_level (vb, ctx);

    // recover state of all downstream contexts - which have one level higher.
    // Note: Invariant: levels on the stack are always monotonously non-decreasing
    ContextP ice_ctx;
    while (stack_pointer && ice[stack_pointer-1].level > vb->peek_stack_level) {
        stack_pointer--;
        memcpy (&ice_ctx, ice[stack_pointer].ctx_p, sizeof (ContextP));
        memcpy (reconstruct_state_start(ice_ctx), ice[stack_pointer].state, RECON_STATE_SIZE);
    }

    if (flag.debug_peek) 
        iprintf ("%s: stack_pointer=%d level=%u %s %s\n", 
                 LN_NAME, stack_pointer, vb->peek_stack_level, is_done_peek ? "DONE" : "POP ", ctx->tag_name);
}

void recon_stack_push (VBlockP vb, ContextP ctx)
{
    // add this context to the frozen state
    buf_alloc (vb, &vb->frozen_state, 1, 10, Ice, 2, "frozen_state");

    Ice *ice = BAFT (Ice, vb->frozen_state);
    ice->level = vb->peek_stack_level;
    memcpy (ice->ctx_p, &ctx, sizeof (ContextP)); // copy pointer to context
    memcpy (ice->state, reconstruct_state_start(ctx), RECON_STATE_SIZE); // actual reconstruction state

    stack_pointer++;

    if (flag.debug_peek)
        iprintf ("%s: stack_pointer=%d level=%u PUSH %s\n", LN_NAME, stack_pointer, ice->level, ctx->tag_name);

    increment_level (vb, ctx);
}

// get reconstructed text without advancing the iterator or changing last_*. context may be already reconstructed or not.
// Note: txt points into txt_data (past reconstructed or BAFT) - caller should copy it elsewhere
ValueType reconstruct_peek (VBlockP vb, ContextP ctx, 
                            pSTRp(txt)) // optional in / out
{
    // case: already reconstructed in this line (or sample in the case of VCF/FORMAT)
    if (ctx_encountered (vb, ctx->did_i) && (ctx->last_encounter_was_reconstructed || (!txt && !txt_len))) {
        if (txt) *txt = last_txtx (vb, ctx);
        if (txt_len) *txt_len = ctx->last_txt.len;
        return ctx->last_value;
    }

    if (flag.debug_peek) 
        iprintf ("%s: stack_pointer=%d level=%u PEEK %s\n", LN_NAME, stack_pointer, vb->peek_stack_level, ctx->tag_name);

    // we "freeze" the reconstruction state of all contexts involved in this peek reconstruction (there could be more than
    // one context - eg items of a container, muxed contexts, "other" contexts (SNIP_OTHER etc). We freeze them by 
    // pushing their state on to a stack in vb->frozen_state when reconstruct_from_ctx is called for this ctx.
    
    uint32_t save_txt_data_len = Ltxt;
    uint32_t save_stack_pointer = stack_pointer;
    uint8_t save_level = vb->peek_stack_level;

    increment_level (vb, ctx); // the next ice that is pushed will be at this new level
    
    reconstruct_from_ctx (vb, ctx->did_i, 0, true);

    ValueType last_value = ctx->last_value;

    // after reconstruction, we still have all the first level ice, we pop it now
    recon_stack_pop (vb, ctx, true);

    // since we are reconstructing unaccounted for data, make sure we didn't go beyond the end of txt_data (this can happen if we are close to the end
    // of the VB, and reconstructed more than OVERFLOW_SIZE allocated in piz_reconstruct_one_vb)
    // note: leave txt_data.len 64bit to detect bugs
    ASSPIZ (vb->txt_data.len <= vb->txt_data.size, "txt_data overflow while peeking %s: len=%"PRIu64" size=%"PRIu64" last_txt_len=%u. vb->txt_data dumped to %s.gz", 
            ctx->tag_name, vb->txt_data.len, (uint64_t)vb->txt_data.size, ctx->last_txt.len, txtfile_dump_vb (vb, z_name));

    // sanity checks
    ASSPIZ (save_level == vb->peek_stack_level, "While peeking %s: expecting peek_stack_level=%u to be %u. Dev tip: try --debug-peek", 
            ctx->tag_name, vb->peek_stack_level, save_level);

    // after this final pop, the stack should be back to where we started
    ASSPIZ (save_stack_pointer == stack_pointer, "While peeking %s: expecting save_stack_pointer=%d to be %u. Dev tip: try --debug-peek", 
            ctx->tag_name, save_stack_pointer, stack_pointer);
    
    // capture data reconstructed for peeking, and "delete" it
    if (txt) *txt = Btxt (save_txt_data_len);
    if (txt_len) *txt_len = Ltxt - save_txt_data_len;

    Ltxt = save_txt_data_len; 

    return last_value;
}

ValueType reconstruct_peek_by_dict_id (VBlockP vb, DictId dict_id, pSTRp(txt)) 
{
    ContextP ctx = ECTX (dict_id); 
    ASSPIZ (ctx, "context doesn't exist for dict_id=%s", dis_dict_id (dict_id).s);

    return reconstruct_peek (vb, ctx, STRa(txt));
}

