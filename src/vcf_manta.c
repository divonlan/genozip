// ------------------------------------------------------------------
//   vcf_manta.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "reconstruct.h"

void vcf_manta_seg_initialize (VBlockVCFP vb)
{
    ContextP id_ctx[7];
    for (char c='0'; c < '7'; c++) 
        id_ctx[(int)c-'0'] = ctx_get_ctx (vb, dict_id_make ((char[]){'I',c,'D'}, 3, DTYPE_FIELD));

    id_ctx[1]->ltype = id_ctx[2]->ltype = id_ctx[4]->ltype = id_ctx[5]->ltype = id_ctx[6]->ltype = LT_DYN_INT;

    id_ctx[2]->flags.store = id_ctx[3]->flags.store = STORE_INT;
}

static bool vcf_seg_manta_ID_cb_3 (VBlockP vb, ContextP ctx, STRp(value), uint32_t unused_rep)
{
    seg_delta_vs_other (VB, ctx, ctx-1, STRa(value));
    return true; // segged successfully
}

void vcf_seg_manta_ID (VBlockVCFP vb, STRp(id))
{
    static MediumContainer con = {
        .nitems_lo = 7,
        .repeats   = 1,
        .items = { { .dict_id.num = DICT_ID_MAKEF_3("I0D"), .separator[0] = ':' },
                   { .dict_id.num = DICT_ID_MAKEF_3("I1D"), .separator[0] = ':' },
                   { .dict_id.num = DICT_ID_MAKEF_3("I2D"), .separator[0] = ':' },
                   { .dict_id.num = DICT_ID_MAKEF_3("I3D"), .separator[0] = ':' },
                   { .dict_id.num = DICT_ID_MAKEF_3("I4D"), .separator[0] = ':' },
                   { .dict_id.num = DICT_ID_MAKEF_3("I5D"), .separator[0] = ':' },
                   { .dict_id.num = DICT_ID_MAKEF_3("I6D"),                     } } };

    SegCallback callbacks[7] = { [3]=vcf_seg_manta_ID_cb_3 }; 

    seg_struct (VB, CTX(VCF_ID), con, STRa(id), callbacks, id_len + 1, true);
}
