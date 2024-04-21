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
#include "libdeflate_1.19/libdeflate.h"

static MediumContainer id_con = {
    .nitems_lo = 7,
    .repeats   = 1,
    .items = { { .dict_id.num = DICT_ID_MAKEF_3("I0D"), .separator[0] = ':' },
               { .dict_id.num = DICT_ID_MAKEF_3("I1D"), .separator[0] = ':' },
               { .dict_id.num = DICT_ID_MAKEF_3("I2D"), .separator[0] = ':' },
               { .dict_id.num = DICT_ID_MAKEF_3("I3D"), .separator[0] = ':' },
               { .dict_id.num = DICT_ID_MAKEF_3("I4D"), .separator[0] = ':' },
               { .dict_id.num = DICT_ID_MAKEF_3("I5D"), .separator[0] = ':' },
               { .dict_id.num = DICT_ID_MAKEF_3("I6D"),                     } } };

static Did tw_dids[NUM_MANTA_TWs] = MANTA_TW_DIDS;

sSTRl(con_id_snip, 128);

void vcf_manta_zip_initialize (void)
{
    DO_ONCE
        container_prepare_snip ((ContainerP)&id_con, 0, 0, qSTRa(con_id_snip));

    // re-initialize for every file, as SV type might change
    segconf.MATEID_method = MATE_01;

    vcf_sv_zip_initialize (tw_dids, NUM_MANTA_TWs);
}   

void vcf_manta_seg_initialize (VBlockVCFP vb)
{
    vcf_sv_seg_initialize (vb, tw_dids, NUM_MANTA_TWs);

    ctx_consolidate_stats_(VB, CTX(VCF_ID), (ContainerP)&id_con);      

    ctx_set_dyn_int (VB, INFO_HOMLEN, INFO_BND_DEPTH, INFO_MATE_BND_DEPTH, DID_EOL);

    for_con2 (&id_con) {
        ContextP ctx = ctx_get_ctx (vb, item->dict_id);
        
        ctx_set_dyn_int (VB, ctx->did_i, DID_EOL);

        if (item_i==2 || item_i==3)
            ctx->flags.store = STORE_INT;
    }

    for (int tw=0; tw < NUM_MANTA_TWs; tw++) 
        seg_mux_init (VB, CTX(tw_dids[tw]), 2, VCF_SPECIAL_DEMUX_BY_MATE, false, (MultiplexerP)&vb->mate_mux[tw]);
}

static bool vcf_seg_manta_ID_cb_3 (VBlockP vb, ContextP ctx, STRp(value), uint32_t unused_rep)
{
    seg_delta_vs_other_localS (VB, ctx, ctx-1, STRa(value), -1);
    return true; // segged successfully
}

void vcf_seg_manta_ID (VBlockVCFP vb, STRp(id))
{
    SegCallback callbacks[7] = { [3]=vcf_seg_manta_ID_cb_3 }; 

    if (id_len > 15 && !memcmp (id, "MantaBND:", 9)) 
        vcf_seg_BND_mate (vb, STRa(id), 0, 0, crc32 (0, id, id_len-2)); // last digit is 0 or 1 - the mate

    if (vcf_has_mate)
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPY_MATE }, 2, VCF_ID, id_len + 1); // +1 for \t

    else
        seg_struct (VB, CTX(VCF_ID), id_con, STRa(id), callbacks, id_len + 1, true);
}

static void vcf_manta_predcited_CIGAR (VBlockVCFP vb, qSTRp (cigar))
{
    if (vb->n_alts == 1) {
        if      (ALT0(DEL)) 
            *cigar_len = snprintf (cigar, *cigar_len, "1M%uD", vb->REF_len-1);
        
        else if (ALT0(INS)) 
            *cigar_len = snprintf (cigar, *cigar_len, "1M%uI", vb->alt_lens[0]-1);
        
        else if (ALT0(SUBST) || ALT0(SUBST_INS) || ALT0(SUBST_DEL)) 
            *cigar_len = snprintf (cigar, *cigar_len, "1M%uI%uD", vb->alt_lens[0]-1, vb->REF_len-1);
    }
    else
        *cigar_len = 0;
} 

void vcf_seg_manta_CIGAR (VBlockVCFP vb, ContextP ctx, STRp(cigar))
{
    STRlic(predicted_cigar,32);
    vcf_manta_predcited_CIGAR (vb, qSTRa(predicted_cigar));

    if (str_issame (cigar, predicted_cigar)) 
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_MANTA_CIGAR }, 2, ctx, cigar_len);
    else
        seg_by_ctx (VB, STRa(cigar), ctx, cigar_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_manta_CIGAR)
{    
    if (reconstruct) {
        STRlic(cigar,32);
        vcf_manta_predcited_CIGAR (VB_VCF, qSTRa(cigar));
        RECONSTRUCT_str (cigar);
    }

    return NO_NEW_VALUE;
}

