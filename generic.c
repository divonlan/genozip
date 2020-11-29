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
    Context *data_ctx = &vb->contexts[GNRIC_DATA];
    data_ctx->ltype = LT_UINT8;
    buf_move (vb, &data_ctx->local, vb, &vb->txt_data);
    data_ctx->txt_len += data_ctx->local.len;

    Context *toplevel_ctx = &vb->contexts[GNRIC_TOPLEVEL];
    toplevel_ctx->no_stons = true; // keep in b250 so it can be eliminated as all_the_same
    
    static const char snip[2] = { SNIP_SPECIAL, GNRIC_SPECIAL_TOPLEVEL };
    seg_by_ctx (vb, snip, 2, toplevel_ctx, 0, 0); 
}

SPECIAL_RECONSTRUCTOR (generic_piz_TOPLEVEL)
{
    buf_destroy (&vb->txt_data);
    buf_move (vb, &vb->txt_data, vb, &vb->contexts[GNRIC_DATA].local);
    return false; // no new value
}
