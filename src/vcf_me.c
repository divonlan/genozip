// ------------------------------------------------------------------
//   vcf_me.c : Mobile Elements
//   Copyright (C) 2024-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"

sSTRl(END_minus_SVLEN_snip, 32);
sSTRl(START_plus_SVLEN_plus_DELTA_snip, 48);
sSTRl(copy_END_snip, 32);
sSTRl(con_meinfo_snip,96);

#define _MEINFO_NAME     DICT_ID_MAKE1_8("M0E_NAME")
#define _MEINFO_START    DICT_ID_MAKE1_8("M1E_STRT")
#define _MEINFO_DELTA    DICT_ID_MAKE1_8("M2E_DLTA")
#define _MEINFO_END      DICT_ID_MAKE1_7("M3E_END")
#define _MEINFO_POLARITY DICT_ID_MAKE1_8("M4E_POLR")

static SmallContainer meinfo_con = { 
    .repeats   = 1, 
    .nitems_lo = 5, 
    .items     = { { .dict_id.num = _MEINFO_NAME,  .separator = ","               }, 
                    { .dict_id.num = _MEINFO_START, .separator = ","               }, 
                    { .dict_id.num = _MEINFO_DELTA, .separator = { CI0_INVISIBLE } }, // MEINFO_END - (MEINFO_START + SVLEN) - calculated but not reconstructed
                    { .dict_id.num = _MEINFO_END,   .separator = ","               }, 
                    { .dict_id.num = _MEINFO_POLARITY                              } }
};

void vcf_me_zip_initialize (void)
{
    DO_ONCE {
        seg_prepare_minus_snip (VCF, _VCF_POS, _INFO_SVLEN, END_minus_SVLEN_snip);
        
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _VCF_POS, true, 0, copy_END_snip); // END is an alias of POS

        seg_prepare_plus_snip (VCF, 3, ((DictId[]){ {_MEINFO_START}, {_INFO_SVLEN}, {_MEINFO_DELTA} }), START_plus_SVLEN_plus_DELTA_snip);
        
        container_prepare_snip ((ContainerP)&meinfo_con, NULL, 0, qSTRa(con_meinfo_snip));
    }
}

void vcf_me_seg_initialize (VBlockVCFP vb)
{
    ctx_consolidate_stats_(VB, CTX(INFO_MEINFO), (ContainerP)&meinfo_con);

    CTX(INFO_ADJLEFT)->flags.same_line = true;
}

// <ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">
// Example: MEINFO=AluYb6_2,10,281,+
void vcf_seg_MEINFO (VBlockVCFP vb, ContextP ctx, STRp(meinfo))
{
    int64_t meinfo_start, meinfo_end, svlen;

    str_split (meinfo, meinfo_len, 4, ',', item, true);
    if (n_items != 4 || !has(SVLEN)                || 
        !str_get_int (STRi(item,1), &meinfo_start) || 
        !str_get_int (STRi(item,2), &meinfo_end)   ||
        !str_get_int (STRa(BII(SVLEN)->value), &svlen)) {

        seg_by_ctx (VB, STRa(meinfo), ctx, meinfo_len);
        return;
    }

    seg_by_dict_id (VB, STRi(item,0), _MEINFO_NAME,     item_lens[0]); 
    seg_by_dict_id (VB, STRi(item,3), _MEINFO_POLARITY, item_lens[3]); 
    
    seg_integer (VB, ctx_get_ctx (VB, _MEINFO_START), meinfo_start, false, item_lens[1]); 

    // MEINFO_END is expected to be close to MEINFO_START + SVLEN 
    int64_t delta = meinfo_end - (meinfo_start + svlen);
    seg_integer (VB, ctx_get_ctx (VB, _MEINFO_DELTA), delta, false, item_lens[2]); // piz sets last_value, but does not reconstruct this

    ContextP end_ctx = ctx_get_ctx (vb, _MEINFO_END);
    end_ctx->flags.same_line = true; // note: not done in vcf_me_seg_initialize to avoid unnecessarily creating the context

    seg_by_ctx (VB, STRa(START_plus_SVLEN_plus_DELTA_snip), end_ctx, 0); // MEINFO_END = MEINFO_START + SVLEN + delta

    seg_by_ctx (VB, STRa(con_meinfo_snip), ctx, 3); // account for 3 commas    
}

void vcf_seg_melt_ADJLEFT (VBlockVCFP vb, ContextP ctx, STRp(adjleft_str))
{
    int64_t end = CTX(VCF_POS)->last_value.i; // END is an alias of POS
    int64_t adjleft, svlen;

    if (has(SVLEN)                                    && 
        ctx_encountered_in_line (VB, INFO_END)        &&
        str_get_int (STRa(BII(SVLEN)->value), &svlen) &&
        str_get_int (STRa(adjleft_str), &adjleft)     &&
        adjleft == end - svlen)

        seg_by_ctx (VB, STRa(END_minus_SVLEN_snip), ctx, adjleft_str_len);
    
    else if (str_issame_(STRa(adjleft_str), "0", 1))
        seg_by_ctx (VB, "0", 1, ctx, 1); // seg as snip

    else
        seg_integer_or_not (VB, ctx, STRa(adjleft_str), adjleft_str_len);
}

void vcf_seg_melt_ADJRIGHT (VBlockVCFP vb, ContextP ctx, STRp(adjright_str))
{
    int64_t end = CTX(VCF_POS)->last_value.i; // END is an alias of POS
    int64_t adjright;

    if (ctx_encountered_in_line (VB, INFO_END)    &&
        str_get_int (STRa(adjright_str), &adjright) &&
        adjright == end)

        seg_by_ctx (VB, STRa(copy_END_snip), ctx, adjright_str_len);

    else if (str_issame_(STRa(adjright_str), "0", 1))
        seg_by_ctx (VB, "0", 1, ctx, 1); // seg as snip

    else
        seg_integer_or_not (VB, ctx, STRa(adjright_str), adjright_str_len);
}

// <ID=INTERNAL,Number=2,Type=String,Description="If insertion internal or close to a gene, listed here followed by a discriptor of the location in the gene (either INTRON, EXON_#, 5_UTR, 3_UTR, PROMOTER, or TERMINATOR)">
// Example: INTERNAL=NM_000384,PROMOTER
void vcf_seg_melt_INTERNAL (VBlockVCFP vb, ContextP ctx, STRp(internal))
{
    MediumContainer con = {
        .repeats   = 1,
        .nitems_lo = 2,
        .items     = { { .dict_id.num = DICT_ID_MAKE1_8("I0N_GENE"), .separator = "," },
                       { .dict_id.num = DICT_ID_MAKE1_8("I1N_DESC")                   } }
    };    

    seg_struct (VB, ctx, con, STRa(internal), NULL, internal_len, true);
}

static bool vcf_seg_melt_DIFF_arr (VBlockP vb, ContextP ctx, STRp(diff_arr), uint32_t repeat)
{
    // TO DO: Possible improvement, need to test: use callback to break the elements down further - t89c 3 item container: -> t c into the one local (lookup with len), 89 into another ctx.local
    seg_array (VB, ctx, INFO_DIFF, STRa(diff_arr), ',', 0, false, false, DICT_ID_NONE, diff_arr_len);
    return true; // segged successfully
}

// <ID=DIFF,Number=.,Type=String,Description="Coverage and Differences in relation to the ALU reference. Form is %2XCoverage:Differences, with differences delimited by ','">
// Example: DIFF=0.94:g73c,t89c,c96a,i127aaa,g145a,c174t,g237c
void vcf_seg_melt_DIFF (VBlockVCFP vb, ContextP ctx, STRp(diff))
{
    MediumContainer con = {
        .repeats   = 1,
        .nitems_lo = 2,
        .items     = { { .dict_id.num = DICT_ID_MAKE1_8("D0FF_VAL"), .separator = ":" },
                       { .dict_id.num = DICT_ID_MAKE1_8("D1FF_ARR")                   } }
    };    

    seg_struct (VB, ctx, con, STRa(diff), (SegCallback[]){ 0, vcf_seg_melt_DIFF_arr }, diff_len, true);
}
