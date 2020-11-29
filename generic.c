// ------------------------------------------------------------------
//   generic.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "vblock.h"
#include "buffer.h"
#include "seg.h"

void generic_seg_finalize (VBlockP vb)
{
    Context *ctx = &vb->contexts[GNRIC_TOPLEVEL];
    ctx->ltype = LT_UINT8;

    static const char snip[2] = { SNIP_SPECIAL, GNRIC_SPECIAL_TOPLEVEL };
    seg_by_ctx (vb, snip, 2, ctx, vb->txt_data.len, 0); 

    buf_move (vb, &vb->contexts[GNRIC_TOPLEVEL].local, vb, &vb->txt_data);
}

SPECIAL_RECONSTRUCTOR (generic_piz_TOPLEVEL)
{
    buf_destroy (&vb->txt_data);
    buf_move (vb, &vb->txt_data, vb, &ctx->local);
    return false; // no new value
}
