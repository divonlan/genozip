// ------------------------------------------------------------------
//   vcf_gnomad.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "piz.h"
#include "reconstruct.h"

sSTRl(copy_QD_snip, 16);
sSTRl(copy_SOR_snip, 16);
sSTRl(copy_MQ_snip, 16);
sSTRl(copy_MQRankSum_snip, 16);
sSTRl(copy_FS_snip, 16);
sSTRl(copy_VarDP_snip, 16);
sSTRl(copy_QUALapprox_snip, 16);
sSTRl(copy_ReadPosRankSum_snip, 16);
sSTRl(con_VRS_End_snip, 48);
sSTRl(con_VRS_States_snip, 48);

#define _INFO_VRS_Ends_REF DICT_ID_MAKE1_8 ("V1RSEnds")
#define _INFO_VRS_Ends_ALT DICT_ID_MAKE1_8 ("V2RSEnds")
#define _INFO_VRS_States_REF DICT_ID_MAKE1_8 ("V3RSStat")
#define _INFO_VRS_States_ALT DICT_ID_MAKE1_8 ("V4RSStat")

void vcf_gnomad_zip_initialize (void)
{
    DO_ONCE {
        seg_prepare_snip_other (SNIP_COPY, _INFO_QD,             false, 0, copy_QD_snip);
        seg_prepare_snip_other (SNIP_COPY, _INFO_SOR,            false, 0, copy_SOR_snip);
        seg_prepare_snip_other (SNIP_COPY, _INFO_MQ,             false, 0, copy_MQ_snip);
        seg_prepare_snip_other (SNIP_COPY, _INFO_MQRankSum,      false, 0, copy_MQRankSum_snip);
        seg_prepare_snip_other (SNIP_COPY, _INFO_FS,             false, 0, copy_FS_snip);
        seg_prepare_snip_other (SNIP_COPY, _INFO_VarDP,          false, 0, copy_VarDP_snip);
        seg_prepare_snip_other (SNIP_COPY, _INFO_QUALapprox,     false, 0, copy_QUALapprox_snip);
        seg_prepare_snip_other (SNIP_COPY, _INFO_ReadPosRankSum, false, 0, copy_ReadPosRankSum_snip);

        SmallContainer con = { .nitems_lo = 2,
                            .repeats   = 1,
                            .items[0]  = { .dict_id.num = _INFO_VRS_Ends_REF, .separator[0] = ',' },
                            .items[1]  = { .dict_id.num = _INFO_VRS_Ends_ALT                      } };

        container_prepare_snip ((ContainerP)&con, 0, 0, qSTRa(con_VRS_End_snip));

        con = (SmallContainer){ .nitems_lo = 2,
                                .repeats   = 1,
                                .items[0]  = { .dict_id.num = _INFO_VRS_States_REF, .separator[0] = ',' },
                                .items[1]  = { .dict_id.num = _INFO_VRS_States_ALT                      } };


        container_prepare_snip ((ContainerP)&con, 0, 0, qSTRa(con_VRS_States_snip));
    }
}

void vcf_gnomad_seg_initialize (VBlockVCFP vb)
{
    // copy from same line (important if INFO fields have been sorted)
    ctx_set_same_line (VB, INFO_AS_QD, INFO_AS_SOR, INFO_AS_MQ, INFO_AS_MQRankSum, INFO_AS_FS,
                       INFO_AS_VarDP, INFO_AS_QUALapprox, INFO_AS_ReadPosRankSum, DID_EOL);

    ctx_set_store (VB, STORE_INT, ctx_get_ctx (VB, _INFO_VRS_Ends_REF)->did_i, ctx_get_ctx (VB, _INFO_VRS_Ends_ALT)->did_i, DID_EOL);

    ctx_consolidate_stats (VB, INFO_VRS_Ends, ctx_get_ctx (VB, _INFO_VRS_Ends_REF)->did_i, ctx_get_ctx (VB, _INFO_VRS_Ends_ALT)->did_i, DID_EOL);
    ctx_consolidate_stats (VB, INFO_VRS_States, ctx_get_ctx (VB, _INFO_VRS_States_REF)->did_i, ctx_get_ctx (VB, _INFO_VRS_States_ALT)->did_i, DID_EOL);
}

void vcf_seg_INFO_VRS_Allele_IDs (VBlockVCFP vb, ContextP ctx, STRp(ids))
{
    static MediumContainer VRS_Allele_IDs_con = {
        .repsep = { ',' },
        .nitems_lo = 1,
        .drop_final_repsep = true,
        .items[0] = { .dict_id.num = DICT_ID_MAKE1_8("V0RSAlID") }
    };

    #define VRS_ALL_ID_PREFIX "\4\4ga4gh:VA.\4" // \4 == CON_PX_SEP

    seg_array_of_struct_with_prefixes (VB, ctx, VRS_Allele_IDs_con, VRS_ALL_ID_PREFIX, STRLEN(VRS_ALL_ID_PREFIX), 
                                       STRa(ids), 
                                       (SegCallback[]){ seg_add_to_local_fixed_len_cb }, // expected to be all of length 32, so better seg as a fixed-len blob than a string
                                       NULL, ids_len);
}

#define VT(x) (vb->var_types[0] == VT_##x)

void vcf_seg_INFO_VRS_Starts (VBlockVCFP vb, ContextP ctx, STRp(arr))
{
    PosType32 pos = DATA_LINE(vb->line_i)->pos;

    str_split_ints (arr, arr_len, 2, ',', start, true);
    
    if (!n_starts || 
        (!VT(SNP) && !VT(DEL) && !VT(INS)) ||
        starts[0] != pos - 1 ||
        (VT(SNP) && starts[1] != pos - 1) ||
        (!VT(SNP) && starts[1] != pos))

        seg_by_ctx (VB, STRa(arr), ctx, arr_len); 

    else // value as predicted
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_VRS_Starts }, 2, ctx, arr_len); 
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_VRS_Starts)
{
    VBlockVCFP vb = (VBlockVCFP)vb_; 
    int64_t pos = CTX(VCF_POS)->last_value.i;

    RECONSTRUCT_INT (pos-1);
    RECONSTRUCT1 (',');
    RECONSTRUCT_INT (pos - VT(SNP));

    return NO_NEW_VALUE;
}

void vcf_seg_INFO_VRS_Ends (VBlockVCFP vb, ContextP ctx, STRp(arr))
{
    str_split_ints (arr, arr_len, 2, ',', end, true);
    
    if (!n_ends)
        seg_by_ctx (VB, STRa(arr), ctx, arr_len); 

    else { 
        seg_by_ctx (VB, STRa(con_VRS_End_snip), ctx, arr_len); 

        seg_delta_vs_other_localN (VB, ctx_get_ctx(vb, _INFO_VRS_Ends_REF), CTX(VCF_POS), ends[0], -1, 0);
        seg_delta_vs_other_localN (VB, ctx_get_ctx(vb, _INFO_VRS_Ends_ALT), CTX(VCF_POS), ends[1], -1, 0);
    }
}

void vcf_seg_INFO_VRS_States (VBlockVCFP vb, ContextP ctx, STRp(arr))
{
    str_split (arr, arr_len, 2, ',', state, true);
    
    if (!n_states)
        seg_by_ctx (VB, STRa(arr), ctx, arr_len); 

    else { 
        seg_by_ctx (VB, STRa(con_VRS_States_snip), ctx, 1); 

        if (str_issame_(STRi(state,0), STRa(vb->REF)))
            seg_by_dict_id (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_ALLELE, '*', '0' }, 4, _INFO_VRS_States_REF, state_lens[0]);
        else
            seg_by_dict_id (VB, STRi(state, 0), _INFO_VRS_States_REF, state_lens[0]);

        if (str_issame_(STRi(state,1), STRi(vb->alt, 0)))
            seg_by_dict_id (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_ALLELE, '*', '1' }, 4, _INFO_VRS_States_ALT, state_lens[1]);
        else
            seg_by_dict_id (VB, STRi(state, 1), _INFO_VRS_States_ALT, state_lens[1]);
    }
}

#define MAYBE_COPY(x) \
void vcf_seg_INFO_AS_##x (VBlockVCFP vb, ContextP ctx, STRp(value)) { seg_maybe_copy ((VBlockP)vb, ctx, INFO_##x, STRa(value), STRa(copy_##x##_snip)); }

MAYBE_COPY(QD)
MAYBE_COPY(SOR)
MAYBE_COPY(MQ)
MAYBE_COPY(MQRankSum)
MAYBE_COPY(FS)
MAYBE_COPY(VarDP)
MAYBE_COPY(QUALapprox)
MAYBE_COPY(ReadPosRankSum)
