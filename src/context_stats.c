// ------------------------------------------------------------------
//   context_stats.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "context.h"
#include "dyn_int.h"

void ctx_consolidate_stats_(VBlockP vb, ContextP parent_ctx, ContainerP con)
{
    uint32_t num_deps = con_nitems (con);

    // find the ultimate ancestor to be displayed in stats
    if (parent_ctx->st_did_i != DID_NONE) 
        parent_ctx = CTX(parent_ctx->st_did_i);

    for (uint32_t d=0; d < num_deps; d++) {
        if (!con.h->items[d].dict_id.num) continue;

        ContextP item_ctx = ctx_get_ctx (vb, con.h->items[d].dict_id);
        if (item_ctx->did_i != parent_ctx->did_i) {
            item_ctx->st_did_i = parent_ctx->did_i;
            item_ctx->header_info = parent_ctx->header_info;
        }
    }

    parent_ctx->is_stats_parent = true;
}

// consolidate a consecutive block of Dids
void ctx_consolidate_statsN (VBlockP vb, Did parent, Did first_dep, unsigned num_deps)
{
    // find the ultimate ancestor to be displayed in stats
    if (CTX(parent)->st_did_i != DID_NONE) 
        parent = CTX(parent)->st_did_i;

    for (ContextP ctx=CTX(first_dep); ctx < CTX(first_dep + num_deps); ctx++)
        if (ctx->did_i != parent) 
            ctx->st_did_i = parent;

    if (CTX(parent)->st_did_i == DID_NONE)
        CTX(parent)->is_stats_parent = true;
}

// consolidate an array of ContextP 
void ctx_consolidate_statsA (VBlockP vb, Did parent, ContextP ctxs[], unsigned num_deps)
{
    // find the ultimate ancestor to be displayed in stats
    if (CTX(parent)->st_did_i != DID_NONE) 
        parent = CTX(parent)->st_did_i;

    for (int i=0; i < num_deps; i++)
        if (ctxs[i]->did_i != parent) 
            ctxs[i]->st_did_i = parent;

    if (CTX(parent)->st_did_i == DID_NONE)
        CTX(parent)->is_stats_parent = true;
}

// simple malloc bc we don't know which thread this is
typedef struct {
    StrText tag_name;
    uint64_t size;
    uint64_t n_words;
    Pointeר buf_nameר;
    LocalType ltype;
    uint8_t dyn_lt_order;
} BigConsumers;

static DESCENDING_SORTER (ctx_big_consumers_sorter, BigConsumers, size);

// called by SIGUSR1, and can run in any thread (read-only access to z_file buffers)
void ctx_show_zctx_big_consumers (FILE *out)
{
    VBlockPoolP pool = vb_get_pool (POOL_MAIN, SOFT_FAIL);

    uint32_t n_ctxs = z_file->ca._num_contexts;  // snapshot lest it grows
    int n_bufs_per_ctx =        IS_ZIP ? 6 : 3; // zctx buffers
    if (pool) n_bufs_per_ctx += IS_ZIP ? 5 : 6; // vctx buffers

    BigConsumers *bc = MALLOC (n_ctxs * n_bufs_per_ctx/*# buffers for context*/ * sizeof (BigConsumers));
    BigConsumers *next = bc;

    ContextP vctx = NULL;

    for (ContextP zctx = ZCTX(0); zctx < ZCTX(n_ctxs); zctx++) { // can't use for zctx for the same reason
        StrText tag = ctx_tag_name_ex (zctx);
        
        *next++     = (BigConsumers){ .tag_name = tag, .buf_nameר = zctx->dict.nameר,        .size = zctx->dict.size };
        *next++     = (BigConsumers){ .tag_name = tag, .buf_nameר = zctx->counts.nameר,      .size = zctx->counts.size };

        if (IS_ZIP) {
            *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = zctx->global_hash.nameר, .size = zctx->global_hash.size };
            *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = zctx->nodes.nameר,       .size = zctx->nodes.size };
            *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = zctx->ston_hash.nameר,   .size = zctx->ston_hash.size };
            *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = zctx->ston_ents.nameר,   .size = zctx->ston_ents.size };
        }
        else {
            *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = zctx->word_list.nameר,   .size = zctx->word_list.size };
            *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = zctx->piz_word_list_hash.nameר, .size = zctx->piz_word_list_hash.size };
        }
        
        if (pool) {
            // for this zctx, get total size of these 4 buffers across all VBs
            // note: we are assuming did_i is equal in vctx and zctx - this will not be true in 
            // the first generation of VBs, so best to test mid-stream, otherwise mis-accounting will occur
            uint64_t dict=0, nodes=0, b250=0, local_hash=0, local=0, dyn_lt_order=0, 
                     dropped_txt=0, history=0, piz_ctx_specific_buf=0;
            LocalType ltype=0;
            uint32_t max_nodes = 0;

            for (VBIDIndex index=0; index < pool->num_vbs; index++) {
                if (!pool->vb[index] || !pool->vb[index]->in_use) continue;

                vctx = &pool->vb[index]->ca.contexts[zctx->did_i];

                if (IS_ZIP) {
                    dict       += vctx->dict.size;
                    nodes      += vctx->nodes.size;
                    b250       += vctx->b250.size;
                    local_hash += vctx->local_hash.size;
                    local      += vctx->local.size;

                    max_nodes = MAX_(vctx->nodes.len32 + vctx->ol_nodes.len32, max_nodes);

                    if (vctx->local.len) {
                        dyn_lt_order = MAX_(dyn_lt_order, vctx->dyn_lt_order);
                        ltype = vctx->ltype;
                    }
                }

                else { // PIZ
                    b250                 += vctx->b250.size;
                    local                += vctx->local.size;
                    dropped_txt          += vctx->dropped_txt.size;
                    history              += vctx->history.size;
                    piz_ctx_specific_buf += vctx->piz_ctx_specific_buf.size;
                }
                
            }

            *next++ =     (BigConsumers){ .tag_name = tag, .buf_nameר = ר(C_LOCAL),      .size = local, .ltype = ltype, .dyn_lt_order = dyn_lt_order };

            if (IS_ZIP) {
                *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = ר(C_DICT),       .size = dict };
                *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = ר(C_NODES),      .size = nodes };
                *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = ר(C_B250),       .size = b250, .n_words = max_nodes };
                *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = ר(C_LOCAL_HASH), .size = local_hash };
            }
            
            else {
                *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = ר(C_B250),           .size = b250, .n_words = zctx->word_list.len };
                *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = ר(C_DROPPED_TXT),    .size = dropped_txt };
                *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = ר(C_HISTORY),        .size = history };
                *next++ = (BigConsumers){ .tag_name = tag, .buf_nameר = ר("piz_ctx_specific_buf"), .size = piz_ctx_specific_buf };
            }
        }
    }   

    ASSERT0 (next - bc == n_bufs_per_ctx * n_ctxs, "bad number of nexts");   

    qsort (bc, n_ctxs * n_bufs_per_ctx, sizeof (BigConsumers), ctx_big_consumers_sorter);

    #define NUM_TO_PRINT 40
    fprintf (out, "%u largest buffers within zctx and (sum of all vctx's):\n", NUM_TO_PRINT);
    if (IS_ZIP) fprintf (out, "Note: ideally this should be run within the ZIP main loop, at the 2nd+ generation of contexts - so not too close to the start or end of the execution\n");

    for (int i=0; i < NUM_TO_PRINT; i++)
        fprintf (out, "%-15s: %-17s: %s%s%s%s\n", 
                 bc[i].tag_name.s, unר(bc[i].buf_nameר), str_size (bc[i].size).s,
                 cond_int (bc[i].n_words, " n_words=", bc[i].n_words),
                 cond_str (!strcmp (unר(bc[i].buf_nameר), C_LOCAL), " ltype=", lt_name (bc[i].ltype)),
                 cond_str (bc[i].dyn_lt_order, " dyn_ltype=", dyn_int_lt_order_name (bc[i].dyn_lt_order)));

    FREE (bc);
}