// ------------------------------------------------------------------
//   sam_minimap2.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "chrom.h"
#include "codec.h"
#include "profiler.h"

// -------------------------------------------------------
// s1:i Chaining score
// -------------------------------------------------------

void sam_seg_s1_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t s1, unsigned add_bytes)
{
    if (has_AS &&
        (dl->AS || sam_seg_get_aux_int (vb, vb->idx_AS_i, &dl->AS, IS_BAM_ZIP, MIN_AS_i, MAX_AS_i, true))) { // AS can be before or after s1 - we have same_line=true

        CTX(OPTION_AS_i)->last_value.i = dl->AS;
        seg_delta_vs_other_do (VB, CTX(OPTION_s1_i), CTX(OPTION_AS_i), NULL, 0, s1, 255, add_bytes);
    }
    else
        seg_integer_as_text_do (VB, CTX(OPTION_s1_i), s1, add_bytes); // unlikely to be ever reached - seg as text to keep context as LT_TEXT
}

// -------------------------------------------------------
// cm:i Number of minimizers on the chain
// -------------------------------------------------------

void sam_seg_cm_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t cm, unsigned add_bytes)
{
    if (segconf.running) {
        // calculate average SEQ.len / cm:i
        segconf.seq_len_to_cm += (cm > 0) ? (int)((float)dl->SEQ.len / (float)cm + 0.5) : 0;
        goto fallback;
    }

    else if (segconf.seq_len_to_cm) {
        int32_t prediction = dl->SEQ.len / segconf.seq_len_to_cm;

        SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_cm, (int32_t)cm - prediction);
        seg_by_did (VB, STRa(snip), OPTION_cm_i, add_bytes); 
    }

    else fallback:
        seg_integer (VB, CTX(OPTION_cm_i), cm, true, add_bytes);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_cm)
{
    int32_t prediction = vb->seq_len / segconf.seq_len_to_cm;

    new_value->i = prediction + atoi(snip);

    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}
