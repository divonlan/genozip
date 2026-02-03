// ------------------------------------------------------------------
//   sam_minimap2.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <math.h>
#include "sam_private.h"
#include "chrom.h"

void sam_minimap2_seg_initialize (VBlockSAMP vb)
{
    ctx_set_store (VB, STORE_INT, OPTION_s1_i, OPTION_s2_i, OPTION_cm_i, DID_EOL);
}

// -------------------------------------------------------
// s1:i Chaining score
// -------------------------------------------------------

static inline int32_t s1_prediction (int32_t cm)
{
    return round ((double)cm * (double)segconf.s1_to_cm_32 / 32.0);
}

void sam_seg_s1_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t s1, unsigned add_bytes)
{
    ContextP ctx = CTX (OPTION_s1_i);
    
    int32_t cm = 0;
    sam_seg_peek_int_field (vb, OPTION_cm_i, vb->idx_cm_i, 1, 1000000, true/*needed for delta*/, &cm); // cm assigned only if successful

    if (segconf_running) {
        // calculate average (s1:i / cm:i) x32
        if (cm && s1 >= 1)
            segconf.s1_to_cm_32 += round (((double)s1 * 32.0) / (double)cm);
        goto fallback;
    }

    int32_t delta = (int32_t)s1 - s1_prediction (cm);

    if (segconf.s1_to_cm_32 && cm) {
        seg_integer (VB, ctx, delta, false, 0);
        seg_special0 (VB, SAM_SPECIAL_s1, ctx, add_bytes);
    }

    else fallback:
        seg_integer (VB, ctx, s1, true, add_bytes);

    ctx_set_last_value (VB, ctx, s1);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_s1)
{
    int32_t cm = reconstruct_peek (vb, CTX(OPTION_cm_i), 0, 0).i;

    int32_t delta = reconstruct_from_local_int (vb, ctx, 0, RECON_OFF);

    new_value->i = s1_prediction (cm) + delta;

    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// -------------------------------------------------------
// s2:i Chaining score of the best secondary chain
// -------------------------------------------------------

void sam_seg_s2_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t s2, unsigned add_bytes)
{
    if (s2 == 0)
        seg_by_did (VB, "0", 1, OPTION_s2_i, add_bytes);

    else if (!ctx_encountered_in_line (VB, OPTION_s1_i)) // just in case - unexpected 
        seg_integer_as_snip_do (VB, CTX(OPTION_s2_i), s2, add_bytes);

    else
        seg_delta_vs_other_localN (VB, CTX(OPTION_s2_i), CTX(OPTION_s1_i), s2, 127, add_bytes);
}

// -------------------------------------------------------
// cm:i Number of minimizers on the chain
// -------------------------------------------------------

void sam_seg_cm_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t cm, unsigned add_bytes)
{
    ContextP ctx = CTX (OPTION_cm_i);

    if (segconf_running) {
        // calculate average SEQ.len / cm:i
        segconf.seq_len_to_cm += (cm > 0) ? round ((double)dl->SEQ.len / (double)cm) : 0;
        goto fallback;
    }

    else if (segconf.seq_len_to_cm) {
        int32_t prediction = dl->SEQ.len / segconf.seq_len_to_cm;
        int32_t delta = (int32_t)cm - prediction;

        seg_integer (VB, ctx, delta, false, 0);
        seg_special0 (VB, SAM_SPECIAL_cm, ctx, add_bytes);
    }

    else fallback:
        seg_integer (VB, ctx, cm, true, add_bytes);

    ctx_set_last_value (VB, ctx, cm);    
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_cm)
{
    int32_t prediction = vb->seq_len / segconf.seq_len_to_cm;

    int32_t delta = snip_len ? atoi(snip) // 14.0.0 to 15.0.75
                             : reconstruct_from_local_int (vb, ctx, 0, RECON_OFF); // starting 15.0.76

    new_value->i = prediction + delta;

    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------
// ms:i DP score of the max scoring segment in the alignment
// ----------------------------------------------------------

// void sam_seg_ms_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t ms, unsigned add_bytes)
// {
//     // if (ms == vb->ref_and_seq_consumed)
// }
