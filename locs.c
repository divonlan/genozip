// ------------------------------------------------------------------
//   locs.c
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "vblock.h"
#include "buffer.h"
#include "seg.h"
#include "locs.h"
#include "dict_id.h"
#include "reconstruct.h"
#include "piz.h"
#include "file.h"

//------------
// ZIP
//------------

// ZIP: called from txtfile_read_header
// returns header length if header read is complete, -1 not complete yet 
int32_t locs_is_header_done (bool is_eof)
{
    return evb->txt_data.len >= 12 ? 12 : -1;
}

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
int32_t locs_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i /* in/out */)
{
    return vb->txt_data.len % 8; // a line is an 8-byte cluster of {float x, y;}
}

bool locs_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return true; // all are small dictionaries
}

void locs_seg_initialize (VBlockP vb)
{
    CTX(LOCS_X)->ltype = LT_FLOAT32;
    CTX(LOCS_Y)->ltype = LT_FLOAT32;
    CTX(LOCS_X)->flags.store = STORE_FLOAT;
    CTX(LOCS_Y)->flags.store = STORE_FLOAT;

    // all-the-same contexts
    ctx_create_node (vb, LOCS_X, (char[]){SNIP_SPECIAL, LOCS_SPECIAL_DELTA_FLOAT }, 2);
    ctx_create_node (vb, LOCS_Y, (char[]){SNIP_SPECIAL, LOCS_SPECIAL_DELTA_FLOAT }, 2);
}

void locs_seg_finalize (VBlockP vb)
{
    SmallContainer top_level_locs = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .nitems_lo    = 2,
        .items        = { { .dict_id = { _LOCS_X  }, .separator = { CI0_TRANS_NOR }, LOCS2LOCS_LTEN_F32 }, 
                          { .dict_id = { _LOCS_Y  }, .separator = { CI0_TRANS_NOR }, LOCS2LOCS_LTEN_F32 } }
    };

    container_seg (vb, CTX(LOCS_TOPLEVEL), (ContainerP)&top_level_locs, 0, 0, 0);
}

rom locs_seg_txt_line (VBlockP vb, rom field_start_line, uint32_t remaining_txt_len, bool *has_special_eol)
{
    float *xy = (float *)field_start_line; 
    uint32_t n_clusters = remaining_txt_len / (sizeof (float) * 2);

    ContextP ctxs[2]={ CTX(LOCS_X), CTX(LOCS_Y)};

    buf_alloc (vb, &ctxs[0]->local, n_clusters, 0, float, 0, "contexts->local");
    buf_alloc (vb, &ctxs[1]->local, n_clusters, 0, float, 0, "contexts->local");

    // Note: volatile, otherwise the gcc optimizer screws up the "precisely reconstructable" test and all
    // clusters appear to be reconstructable even if they're not.
    volatile float last[2] = {}; 
    for (uint32_t i=0; i < n_clusters * 2; i++) {
        float X_or_Y = LTEN32F(xy[i]);
        float delta  = X_or_Y - last[i&1];
        
        // test if X_or_Y is precisely reconstructable from delta
        double test_correct = X_or_Y;
        double test_recon   = (double)last[i&1] + (double)delta;
        
        if ((float)test_correct == (float)test_recon) {
            BNXT (float, ctxs[i&1]->local) = delta; // BGEN in zip_generate_local
            seg_known_node_index (vb, ctxs[i&1], 0, sizeof (float));
        }
        
        // case: X_or_Y is not precisely reconstructable from delta due to dropping mantissa bits
        else {       
            BNXT (float, ctxs[i&1]->local) = X_or_Y; // store whole value, without delta
            seg_simple_lookup (vb, ctxs[i&1], sizeof (float));
        }
        last[i&1] = X_or_Y;
    }

    vb->lines.len = n_clusters;
    vb->line_i = n_clusters - 1;

    return field_start_line + n_clusters * (sizeof (float) * 2);
}

//------------ 
// PIZ
//------------

SPECIAL_RECONSTRUCTOR (locs_piz_special_DELTA_FLOAT)
{
    float delta = NEXTLOCAL (float, ctx);
    new_value->f = (float)ctx->last_value.f + delta;

    return true; // has new_value 
}
