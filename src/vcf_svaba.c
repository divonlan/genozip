// ------------------------------------------------------------------
//   vcf_svaba.c
//   Copyright (C) 2024-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "zip_dyn_int.h"

static Did tw_dids[NUM_SVABA_TWs] = SVABA_TW_DIDS;

static SmallContainer id_con = {
    .nitems_lo = 2,
    .repeats   = 1,
    .items = { { .dict_id.num = DICT_ID_MAKEF_3("I0D"), .separator[0] = ':' },
               { .dict_id.num = DICT_ID_MAKEF_3("I1D"),                     } } };

sSTRl(con_id_snip, 64);

void vcf_svaba_zip_initialize (void)
{
    DO_ONCE {
        container_prepare_snip ((ContainerP)&id_con, 0, 0, qSTRa(con_id_snip));
    }

    // re-initialize for every file, as SV type might change
    vcf_sv_zip_initialize (tw_dids, NUM_SVABA_TWs);
}   

void vcf_svaba_seg_initialize (VBlockVCFP vb)
{
    vcf_sv_seg_initialize (vb, tw_dids, NUM_SVABA_TWs);

    ctx_consolidate_stats_(VB, CTX(VCF_ID), (ContainerP)&id_con);      

    segconf.MATEID_method = MATE_12;
}

void vcf_seg_svaba_ID (VBlockVCFP vb, STRp(id))
{
    decl_ctx(VCF_ID);

    str_split_ints (id, id_len, 2, ':', id, true);
    if (n_ids != 2) {
        seg_by_ctx (VB, STRa(id), ctx, id_len);
        return;
    }

    vcf_seg_BND_mate (vb, STRa(id), 0, 0, ids[0]);

    if (vcf_has_mate)
        seg_special0 (VB, VCF_SPECIAL_COPY_MATE, ctx, id_len + 1); // +1 for \t

    else {
        seg_integer (VB, ctx_get_ctx (vb, id_con.items[0].dict_id), ids[0], false, id_len - 2); // account without the eg :1
        seg_integer (VB, ctx_get_ctx (vb, id_con.items[1].dict_id), ids[1], false, 1); // account for 1 or 2
        seg_by_ctx (VB, STRa(con_id_snip), ctx, 2); // account for the colon and \t
    }
}

void vcf_seg_svaba_MATEID (VBlockVCFP vb, ContextP ctx, STRp(mate_id))
{
    STRlast (id, VCF_ID);

    if (id_len == mate_id_len && str_issame_(id, id_len-1, mate_id, mate_id_len-1) &&
        ((id[id_len-1]=='1' && mate_id[id_len-1]=='2') || (id[id_len-1]=='2' && mate_id[id_len-1]=='1')))

        seg_special0 (VB, VCF_SPECIAL_SVABA_MATEID, ctx, mate_id_len);
    else
        seg_by_ctx (VB, STRa(mate_id), ctx, mate_id_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_SVABA_MATEID)
{    
    if (reconstruct) {
        STRlast (id, VCF_ID);

        RECONSTRUCT (id, id_len-1);
        RECONSTRUCT1 (id[id_len-1]=='1' ? '2' : '1');
    }

    return NO_NEW_VALUE;
}

void vcf_seg_svaba_MAPQ (VBlockVCFP vb, ContextP ctx, STRp(mapq))
{
    ContextP channel_ctx = segconf.has[INFO_MATEMAPQ] ? vcf_seg_sv_copy_mate (vb, CTX(INFO_MAPQ), STRa(mapq), TW_MAPQ, TW_MATEMAPQ, true, mapq_len)
                         :                              ctx;

    // case: NOT segged as copy-mate  
    if (channel_ctx) { 
        bool has_disc = ctx_encountered_in_line (VB, INFO_DISC_MAPQ);

        // prediction: if DISC_MAPQ exists, MAPQ is the same, if not, it is 60
        if ((has_disc && str_issame_(STRa(mapq), STRlst(INFO_DISC_MAPQ))) ||
            (!has_disc && str_issame_(STRa(mapq), "60", 2)))

            seg_special0 (VB, VCF_SPECIAL_MAPQ, channel_ctx, mapq_len);
        else
            seg_by_ctx (VB, STRa(mapq), channel_ctx, mapq_len);
    }
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MAPQ)
{    

    if (reconstruct) {
        if (curr_container_has (vb, (DictId)_INFO_DISC_MAPQ)) {
            STR(disc_mapq);
            reconstruct_peek (vb, CTX(INFO_DISC_MAPQ), pSTRa(disc_mapq));
            RECONSTRUCT_str (disc_mapq);
        }

        else
            RECONSTRUCT ("60", 2);
    }

    return NO_NEW_VALUE;
}

void vcf_seg_svaba_SPAN (VBlockVCFP vb, ContextP ctx, STRp(span_str))
{
    int64_t span;
    bool is_bnd = VT0(BND); 
    bool same_chrom = is_bnd && str_issame (vb->chrom_name, vb->mate_chrom_name);
    PosType32 pos = DATA_LINE(vb->line_i)->pos;

    if (is_bnd && str_get_int (STRa(span_str), &span) &&
         ( (same_chrom  && span == ABS (pos - vb->mate_pos)) || // predicted SPAN:
           (!same_chrom && span == -1)))

        seg_special0 (VB, VCF_SPECIAL_SPAN, ctx, span_str_len);
    else
        seg_by_ctx (VB, STRa(span_str), ctx, span_str_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_SPAN)
{    
    STRlast (chrom_name, VCF_CHROM);
    bool same_chrom = str_issame (chrom_name, VB_VCF->mate_chrom_name);

    new_value->i = same_chrom ? ABS (CTX(VCF_POS)->last_value.i - VB_VCF->mate_pos) : -1;

    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}
