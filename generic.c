// ------------------------------------------------------------------
//   generic.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "vblock.h"
#include "buffer.h"
#include "seg.h"
#include "generic.h"
#include "dict_id.h"

// all data is always consumed
int32_t generic_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    return 0;
}

void generic_seg_finalize (VBlockP vb)
{
    Context *data_ctx = CTX(GNRIC_DATA);
    data_ctx->ltype = LT_UINT8;
    buf_move (vb, &data_ctx->local, vb, &vb->txt_data);
    data_ctx->txt_len += data_ctx->local.len;

    Context *toplevel_ctx = CTX(GNRIC_TOPLEVEL);
    toplevel_ctx->no_stons = true; // keep in b250 so it can be eliminated as all_the_same
    
    static const char snip[2] = { SNIP_SPECIAL, GNRIC_SPECIAL_TOPLEVEL };
    seg_by_ctx (VB, snip, 2, toplevel_ctx, 0); 
}

bool generic_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return true; // contexts are expected to have small dictionaries
}

SPECIAL_RECONSTRUCTOR (generic_piz_TOPLEVEL)
{
    buf_destroy (&vb->txt_data);
    buf_move (vb, &vb->txt_data, vb, &CTX(GNRIC_DATA)->local);
    return false; // no new value
}
