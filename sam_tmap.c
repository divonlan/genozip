// ------------------------------------------------------------------
//   sam_gem3.c
//   Copyright (C) 2022-2022 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted,
//   under penalties specified in the license.

// Compresses auxilliary fields generated by tmap (IonXpress) mapper

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "reference.h"
#include "segconf.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"

void sam_seg_TMAP_XM_i (VBlockSAMP vb, ValueType XM, unsigned add_bytes)
{
    // XM is predicted to be ref_consumed
    if (XM.i == vb->ref_consumed)                              
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED }, 2, OPTION_XM_i, add_bytes);
        
    else
        seg_integer (VB, CTX(OPTION_XM_i), XM.i, true, add_bytes);
}

// optimization for Ion Torrent flow signal (ZM) - negative values become zero, positives are rounded to the nearest 10
void sam_optimize_TMAP_ZM (VBlockSAMP vb, ContextP ctx, void *cb_param, void *array_, uint32_t array_len)
{
    int16_t *array = (int16_t *)array_;

    for (uint32_t i=0; i < array_len; i++)
        if (array[i] >= 0) 
#ifdef __BIG_ENDIAN__
            array[i] = LTEN16 (((LTEN16(array[i]) + 5) / 10) * 10);
#else            
            array[i] = ((array[i] + 5) / 10) * 10;
#endif
        else array[i] = 0;
}

