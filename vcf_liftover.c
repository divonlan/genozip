// ------------------------------------------------------------------
//   vcf_liftover.c
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// Dual-coordinates genozip/genocat file flow:
//
// genozip -C source.vcf       --> dual-coord.genozip - component_i=0 contains all lines, each with LUFT or LIFTREJT 
//                                                                    LIFTREFD are also written to the rejects file
//                                                      component_i=1 is the rejects file containing LIFTREFD lines (again)
//
// genocat dual-crd.genozip    --> primary.vcf -        component_i=0 is reconstructed - lines have LUFT or LIFTREJT 
//                                                                    lines are IN ORDER
//                                                      component_i=1 is skipped.
//
// genozip primary.vcf         --> dual-coord.genozip - component_i=0 contains all lines, each with LUFT or LIFTREJT 
//                                                                    LIFTREFD are also written to the rejects file
//                                                      component_i=1 is the rejects file containing LIFTREFD lines (again)
//
// genocat -v dual-crd.genozip --> luft.vcf -           component_i=1 with LIFTREJT is reconstructed first and becomes part of the header
//                                                      component_i=0 is reconstructed - dropping LIFTREJT lines and lifting over LUFT->PRIM 
//                                                                    lines sorted by writer according the the RECON_PLAN section 
//
// genozip luft.vcf            --> dual-coord.genozip - LIFTREJT header lines are sent to vblock_i=1 via unconsumed_txt
//                                                      component_i=0 contains all lines first all LIFTREJT lines followed by all LUFT lines  
//                                                                    LIFTREFD lines are also written to the rejects file
//                                                      component_i=1 is the rejects file containing LIFTREFD lines (again)
//
// genocat dual-crd.genozip    --> primary.vcf          component_i=0 is reconstructed - first LIFTREJT then LUFT lines
//                                                                    lines are OUT OF ORDER (LIFTREJT first)
//                                                      component_i=1 is skipped.
//

#include <errno.h>
#include <math.h>
#include "vcf_private.h"
#include "context.h"
#include "flags.h"
#include "chain.h"
#include "data_types.h"
#include "file.h"
#include "buffer.h"
#include "vblock.h"
#include "container.h"
#include "seg.h"
#include "dict_id.h"
#include "codec.h"
#include "strings.h"
#include "random_access.h"
#include "reconstruct.h"
#include "coords.h"

const char *dvcf_status_names[] = DVCF_STATUS_NAMES;
const LuftTransLateProp ltrans_props[NUM_VCF_TRANS] = DVCF_TRANS_PROPS;

// constant snips initialized at the beginning of the first file zip
#define LO_SNIP_LEN 200
static char info_luft_snip[LO_SNIP_LEN]={}, info_prim_snip[LO_SNIP_LEN]={}, info_rejt_luft_snip[LO_SNIP_LEN]={}, info_rejt_prim_snip[LO_SNIP_LEN]={};
static unsigned info_luft_snip_len=LO_SNIP_LEN, info_prim_snip_len=LO_SNIP_LEN, info_rejt_luft_snip_len=LO_SNIP_LEN, info_rejt_prim_snip_len=LO_SNIP_LEN;
static bool lo_snip_initialized = false;

// ---------------
// ZIP & SEG stuff
// ---------------

// ZIP: called by the main thread from *_zip_initialize 
void vcf_lo_zip_initialize (void)
{
    if (lo_snip_initialized) return;
    lo_snip_initialized = true;

    // case: --chain : create SRCNAME context from liftover contigs and dict
    if (chain_is_loaded && z_file->num_txt_components_so_far == 1)  // copy only once
        chain_copy_contigs_to_z_file (DTFZ (ochrom));

    // prepare (constant) snips. note: we need these snips both for --chain and when zipping dual-coord files
    SmallContainer con = {
        .repeats     = 1,
        .nitems_lo   = NUM_IL_FIELDS,
        .drop_final_item_sep = true,
        .items       = { [IL_CHROM  ]={ .dict_id = (DictId)dict_id_fields[VCF_oCHROM],   .seperator = ","  },
                         [IL_POS    ]={ .dict_id = (DictId)dict_id_fields[VCF_oPOS],     .seperator = ","  },
                         [IL_REF    ]={ .dict_id = (DictId)dict_id_fields[VCF_LIFT_REF], .seperator = ","  },
                         [IL_XSTRAND]={ .dict_id = (DictId)dict_id_fields[VCF_oXSTRAND], .seperator = ","  } } };
    container_prepare_snip ((Container*)&con, 0, 0, info_luft_snip, &info_luft_snip_len);

    con = (SmallContainer){
        .repeats     = 1,
        .nitems_lo   = NUM_IL_FIELDS,
        .drop_final_item_sep = true,
        .items       = { [IL_CHROM  ]={ .dict_id = (DictId)dict_id_fields[VCF_CHROM],    .seperator = ","  },
                         [IL_POS    ]={ .dict_id = (DictId)dict_id_fields[VCF_COPYPOS],  .seperator = ","  },
                         [IL_REF    ]={ .dict_id = (DictId)dict_id_fields[VCF_LIFT_REF], .seperator = ","  },
                         [IL_XSTRAND]={ .dict_id = (DictId)dict_id_fields[VCF_oXSTRAND], .seperator = ","  } } };
    container_prepare_snip ((Container*)&con, 0, 0, info_prim_snip, &info_prim_snip_len);

    // for REJTOVER, appearing a line that has only PRIMARY coordinates, we include oCHROM, oPOS and oREFALT which will 
    // be segged as "" (empty snip) for this variant that cannot be lifted over. The reason is that when reconstructing 
    // with --luft, the reconstructor will reconstruct the entire line (consuming oCHROM, oPOS, oREFALT for the main VCF fields)
    // before finally being dropped by the container callback. Likewise for REJTBACK.
    con = (SmallContainer){
        .repeats     = 1,
        .nitems_lo   = 4,
        .items       = { { .dict_id = (DictId)dict_id_fields[VCF_COPYSTAT] },
                         { .dict_id = (DictId)dict_id_fields[VCF_oCHROM]   },
                         { .dict_id = (DictId)dict_id_fields[VCF_oPOS]     },
                         { .dict_id = (DictId)dict_id_fields[VCF_oREFALT]  } } };
    container_prepare_snip ((Container*)&con, 0, 0, info_rejt_luft_snip, &info_rejt_luft_snip_len);

    con = (SmallContainer){
        .repeats     = 1,
        .nitems_lo   = 4,
        .items       = { { .dict_id = (DictId)dict_id_fields[VCF_COPYSTAT] },
                         { .dict_id = (DictId)dict_id_fields[VCF_CHROM]    },
                         { .dict_id = (DictId)dict_id_fields[VCF_COPYPOS], .seperator = { CI_TRANS_NOR } }, // rather than segging "", we don't reconstruct. so we don't break the "all_the_same" of COPYPOS
                         { .dict_id = (DictId)dict_id_fields[VCF_REFALT]   } } };
    container_prepare_snip ((Container*)&con, 0, 0, info_rejt_prim_snip, &info_rejt_prim_snip_len);
}

// returns true if dict_id is AF_* or *_AF, excluding MAX_AF
static inline bool vcf_lo_is_INFO_AF_type (DictId dict_id)
{
    // case AF_*
    if (dict_id.id[0] == ('A' | 0xc0) && dict_id.id[1] == 'F' && dict_id.id[2] == '_') return true;

    // get dict_id string length (up to DICT_ID_LEN)
    unsigned len=DICT_ID_LEN;
    for (int i=DICT_ID_LEN-1; i >= 0; i--)
        if (!dict_id.id[i]) len--;
        else break;

    // case *_AF. Note this works on long names like "gnomAD_ASJ_AF" too due to how dict_id_make works
    // note: we exclude MAX_AF as it is the maximum of all population AFs - it will not be correct after flipping
    if (len > 3 && dict_id_is_vcf_info_sf (dict_id) &&  
        dict_id.id[len-3] == '_' && dict_id.id[len-2] == 'A' && dict_id.id[len-1] == 'F' && 
        dict_id.num != dict_id_INFO_MAX_AF)
        return true;

    return false;
}

// get liftover translator ID. For contexts for which there is an ##INFO or ##FORMAT field in the VCF header,
// these are entered to the context by vcf_header. Fields appearing in the file are queries by Seg. 
TranslatorId vcf_lo_luft_trans_id (DictId dict_id, char number)
{
    // note: in VCF 4.1 many files specified A,G and R fields as Number=. ; 'R' was introduced in 4.2. We set some common fields
    // regardless of their "Number" in the header.
    
    if (number == 'R' || dict_id.num == dict_id_FORMAT_AD || dict_id.num == dict_id_FORMAT_ADR || dict_id.num == dict_id_FORMAT_ADF || 
        dict_id.num == dict_id_FORMAT_ADALL || dict_id.num == dict_id_FORMAT_F2R1 || dict_id.num == dict_id_FORMAT_F1R2 ||
        dict_id.num == dict_id_INFO_DP_HIST || dict_id.num == dict_id_INFO_GQ_HIST) 
        return VCF2VCF_R;

    else if (number == 'G' || dict_id.num == dict_id_FORMAT_GL || dict_id.num == dict_id_FORMAT_PL || dict_id.num == dict_id_FORMAT_PRI ||
             dict_id.num == dict_id_FORMAT_GP || dict_id.num == dict_id_FORMAT_PP)
        return VCF2VCF_G;

    else if (dict_id.num == dict_id_FORMAT_GT)       return VCF2VCF_GT;
    else if (dict_id.num == dict_id_INFO_END )       return VCF2VCF_END;
    else if (dict_id.num == dict_id_INFO_AC
    ||       dict_id.num == dict_id_INFO_MLEAC)      return VCF2VCF_A_AN;
    else if (dict_id.num == dict_id_INFO_AF 
    ||       dict_id.num == dict_id_INFO_MLEAF  
    ||       dict_id.num == dict_id_INFO_LDAF  
    ||       dict_id.num == dict_id_FORMAT_AF  
    ||       vcf_lo_is_INFO_AF_type (dict_id))
          return VCF2VCF_A_1; // eg 0.150 -> 0.850 (upon REF<>ALT SWITCH)
    else if (dict_id.num == dict_id_FORMAT_DS)       return VCF2VCF_PLOIDY; // eg (if ploidy=2): 1.25 -> 0.75 
    else if (dict_id.num == dict_id_INFO_AA)         return VCF2VCF_ALLELE; 
    else if (dict_id.num == dict_id_FORMAT_SAC       // eg 25,1,6,2 --> 6,2,25,1 (upon REF<>ALT SWITCH)
    ||       dict_id.num == dict_id_FORMAT_SB        // SB and MB - I don't fully understand these, but at least for bi-allelic variants, empirically the sum of the first two values equals first value of AD and the sum of last two values adds up to the second value of AD
    ||       dict_id.num == dict_id_FORMAT_MB)       return VCF2VCF_R2;   
    else if (dict_id.num == dict_id_INFO_BaseCounts) return VCF2VCF_XREV; // eg 10,6,3,4 -> 4,3,6,10  (upon XSTRAND)
    
    return TRANS_ID_NONE;
}

// SEG: map coordinates primary->luft. In case of failure, sets dst_contig_index, dst_1pos, xstrand to 0.
LiftOverStatus vcf_lo_get_liftover_coords (VBlockVCFP vb, PosType pos, WordIndex *dst_contig_index, PosType *dst_1pos, bool *xstrand, uint32_t *aln_i) // out
{
    #define WORD_INDEX_MISSING -2
    #define WORD_INDEX_NOT_IN_PRIM_REF -3

    // extend vb->chrom_map_vcf_to_chain if needed - map between VCF node index (NOT word index) and chain primary contigs (NOT reference primary contigs)
    if (vb->chrom_node_index >= vb->chrom_map_vcf_to_chain.len) { 
        buf_alloc (vb, &vb->chrom_map_vcf_to_chain, 0, MAX (vb->chrom_node_index+1, chain_get_num_prim_contigs()), WordIndex, 2, "chrom_map_vcf_to_chain");

        // initialize new entries allocated to WORD_INDEX_MISSING
        uint32_t size = vb->chrom_map_vcf_to_chain.size / sizeof (WordIndex);
        for (uint32_t i=vb->chrom_map_vcf_to_chain.param; i < size; i++)
            *ENT (WordIndex, vb->chrom_map_vcf_to_chain, i) = WORD_INDEX_MISSING;

        vb->chrom_map_vcf_to_chain.param = size; // param holds the number of entries initialized >= len
        vb->chrom_map_vcf_to_chain.len   = vb->chrom_node_index + 1;
    }

    ARRAY (WordIndex, map, vb->chrom_map_vcf_to_chain);

    // if chrom is not yet mapped to src_contig, map it now:
    // map from chrom_node_index (not word index!) to entry in chain_index
    if (map[vb->chrom_node_index] == WORD_INDEX_MISSING) {
        map[vb->chrom_node_index] = chain_get_prim_contig_index_by_name (vb->chrom_name, vb->chrom_name_len);

        // check if the error is due to missing contig in primary reference, or only chain file
        if (map[vb->chrom_node_index] == WORD_INDEX_NONE &&
            ref_contigs_get_by_name (prim_ref, vb->chrom_name, vb->chrom_name_len, false, true) == WORD_INDEX_NONE)
            map[vb->chrom_node_index] = WORD_INDEX_NOT_IN_PRIM_REF;
    }

    if (map[vb->chrom_node_index] == WORD_INDEX_NONE || map[vb->chrom_node_index] == WORD_INDEX_NOT_IN_PRIM_REF) {
        if (dst_contig_index) *dst_contig_index = WORD_INDEX_NONE;
        if (dst_1pos) *dst_1pos = 0;
        if (xstrand) *xstrand = false;
        
        return map[vb->chrom_node_index] == WORD_INDEX_NONE ? LO_CHROM_NOT_IN_CHAIN : LO_CHROM_NOT_IN_PRIM_REF;
    }
    
    else {
        bool mapping_ok = chain_get_liftover_coords (map[vb->chrom_node_index], pos, dst_contig_index, dst_1pos, xstrand, aln_i); // if failure, sets output to 0.
        return mapping_ok ? LO_OK : LO_NO_MAPPING_IN_CHAIN;
    }
}                                

// Cross-render Luft values to Primary (when segging a Luft file)
// returns true if successfully translated, false if left as-is
// Notes: 1. all translators are symmetrical - running on primary results in luft and vice versa
// 2. AC and END translators might result in bigger or smaller text - we handle this
bool vcf_lo_seg_cross_render_to_primary (VBlockVCFP vb, ContextP ctx, 
                                         const char *this_value, unsigned this_value_len, 
                                         char *modified_snip, unsigned *modified_snip_len) // NULL means translate in-place (MUST be the same length)
{
    // save original value (and TXTFILE_READ_VB_PADDING extra bytes) and txt_data.len
    uint64_t save_txt_data_len = vb->txt_data.len;
    char save[this_value_len + TXTFILE_READ_VB_PADDING];
    
    if (modified_snip) 
        memcpy (save, this_value, this_value_len + TXTFILE_READ_VB_PADDING);

    // set txt_data as if we just reconsturcted value - len ending after value
    uint64_t len_before = vb->txt_data.len = ENTNUM (vb->txt_data, this_value + this_value_len);

    // convert Primary->Luft in place. Note: if untranslatable (because it was an untranslatable Primary value
    // to be begin with), this function will do nothing.
    if (DT_FUNC(vb, translator)[ctx->luft_trans]((VBlockP)vb, ctx, (char *)this_value, this_value_len, false)) {

        int len_diff = (int)vb->txt_data.len - (int)len_before;
        
        ASSERT (!len_diff || modified_snip, "len_diff=%d, but translation in place of %s is only allowed if len_diff=0", len_diff, ctx->name);

        // length only changes for AC and END - and since they are 32 bit ints in the VCF spec - they change by at most 9 
        // characters (1 character -> 10 characters)
        ASSERT (len_diff <= TXTFILE_READ_VB_PADDING, "When translating luft->primary of %s, length grew by %d - that's too much",
                ctx->name, len_diff);

        // copy result
        if (modified_snip) {
            *modified_snip_len = (int)this_value_len + len_diff;
            vb->recon_size -= (int)this_value_len - (int)*modified_snip_len; // Primary reconstruction is smaller

            memcpy (modified_snip, this_value, *modified_snip_len);

            // restore original
            memcpy ((char *)this_value, save, this_value_len + TXTFILE_READ_VB_PADDING);
        }

        vb->txt_data.len = save_txt_data_len;

        return true;
    }
    else { // failed lift over
        vcf_lo_seg_rollback_and_reject (vb, dict_id_is_vcf_format_sf (ctx->dict_id) ? LO_FORMAT : LO_INFO, ctx);
        vb->txt_data.len = save_txt_data_len;
        return false;
    }
}

static inline Context *vcf_lo_seg_lo_snip (VBlockVCFP vb, const char *snip, unsigned snip_len, uint64_t dict_id_num, unsigned add_bytes)
{
    Context *ctx = ctx_get_ctx (vb, dict_id_num);
    ctx_set_encountered_in_line (ctx);
    ctx->no_stons = true; // keep in b250 so it can be eliminated as all_the_same
    ctx->st_did_i = VCF_oSTATUS;
    seg_by_ctx (vb, snip, snip_len, ctx, add_bytes); 
    
    return ctx;
}

// Segs both INFO/REJTOVER and INFO/REJTBACK container and their direct items, but NOT ostatus
static void vcf_lo_seg_REJX_do (VBlockVCFP vb, unsigned add_bytes)
{
    ZipDataLineVCF *dl = DATA_LINE (vb->line_i);
    
    dl->chrom_index[SEL(1,0)] = WORD_INDEX_NONE;
    dl->pos[SEL(1,0)] = 0;

    // whether user views the file primary or luft coordinates, and whether this line is primary-only or luft-only:
    // each rejected line consumes all contexts exactly once, either in the main field, or in INFO/*rej
    // note: rejected lines are still reconstructed in their unsupported coordinates to consume all contexts, and then dropped.
    vcf_lo_seg_lo_snip (vb, info_rejt_luft_snip, info_rejt_luft_snip_len, dict_id_INFO_LREJ, 0); 
    vcf_lo_seg_lo_snip (vb, info_rejt_prim_snip, info_rejt_prim_snip_len, dict_id_INFO_PREJ, 0);
    
    // case: primary-only line: CHROM, POS and REFALT were segged in the main fields, trivial oCHROM, oPOS, oREFLT segged here
    // case: luft-only line: oCHROM, oPOS and oREFALT were segged in the main fields, trivial CHROM, POS, REFLT segged here
    seg_by_did_i (vb, "", 0, SEL (VCF_oCHROM,  VCF_CHROM),  0);  // 0 for all of these as these fields are not reconstructed 
    seg_by_did_i (vb, "", 0, SEL (VCF_oREFALT, VCF_REFALT), 0);
    seg_by_did_i (vb, "", 0, SEL (VCF_oPOS,    VCF_POS),    0);

    CTX(VCF_COPYSTAT)->txt_len += add_bytes; // we don't need to seg COPYSTAT because it is all_the_same. Just account for add_bytes.

    // if we're segging a Primary line, and rejecting the Luft rendition, then the reject will be only reconstructed in primary coordinates, and vice verse
    // we modified the text by adding this string (the "REJX" is accounted for in vcf_seg_info_add_DVCF_to_InfoItems)
    if (vb->line_coords == DC_PRIMARY) vb->recon_size      += add_bytes; 
    else                               vb->recon_size_luft += add_bytes;
}

#define NUM_ROLLBACK_CTXS 6
static const DidIType line_rollback_did_i[2][NUM_ROLLBACK_CTXS] = { 
    { VCF_oSTATUS, VCF_oCHROM, VCF_oPOS, VCF_oREFALT, VCF_LIFT_REF, VCF_oXSTRAND },
    { VCF_oSTATUS, VCF_CHROM,  VCF_POS,  VCF_REFALT,  VCF_LIFT_REF, VCF_oXSTRAND }
};

// Actually, reject this line (Primary or Luft line) from cross-rendering, possibly after already segging to o* fields in vcf_lo_seg_generate_INFO_DVCF
void vcf_lo_seg_rollback_and_reject (VBlockVCFP vb, LiftOverStatus ostatus, Context *ctx /* INFO or FORMAT context, or NULL */)
{
    if (LO_IS_REJECTED (last_ostatus)) return; // already rejected

    // unaccount for INFO/LUFT (in primary reconstruction) and INFO/PRIM (in Luft reconstruction)
    Context *luft_ctx; 
    if (ctx_encountered_in_line (vb, dict_id_INFO_LUFT, &luft_ctx)) { // we always have both LUFT and PRIM, or neither
        vb->recon_size      -= luft_ctx->last_txt_len;
        vb->recon_size_luft -= ctx_get_existing_ctx (vb, dict_id_INFO_PRIM)->last_txt_len;
    }

    for (unsigned i=0; i < NUM_ROLLBACK_CTXS; i++) 
        seg_rollback ((VBlockP)vb, &vb->contexts[line_rollback_did_i[SEL(0,1)][i]]);

    seg_rollback ((VBlockP)vb, ctx_get_ctx (vb, dict_id_INFO_LUFT));
    seg_rollback ((VBlockP)vb, ctx_get_ctx (vb, dict_id_INFO_PRIM));

    // seg Xrej and empty values for the fields of the rejected coord
    char str[MAX_VCF_ID_LEN + 20];
    const char *id = ctx ? vcf_header_get_VCF_ID_by_dict_id (ctx->dict_id, false) : NULL;
    sprintf (str, "%s%s%s", dvcf_status_names[ostatus], ctx ? "/" : "", 
             id ? id : ctx ? dis_dict_id (ctx->dict_id).s : "");
    unsigned str_len = strlen (str);

    vcf_lo_seg_REJX_do (vb, str_len); 

    // seg oSTATUS
    seg_by_did_i (vb, str, str_len, VCF_oSTATUS, 0);
    vcf_set_ostatus (ostatus); // ONLY this place can set a rejection ostatus, or else it will fail the if at the top of this function
}

void vcf_lo_set_rollback_point (VBlockVCFP vb)
{
    for (unsigned i=0; i < NUM_ROLLBACK_CTXS; i++) 
        seg_create_rollback_point (&vb->contexts[line_rollback_did_i[SEL(0,1)][i]]);

    seg_create_rollback_point (ctx_get_ctx (vb, dict_id_INFO_LUFT));
    seg_create_rollback_point (ctx_get_ctx (vb, dict_id_INFO_PRIM));
}

// --chain: Create:
// Dual coordinates line: LUFT and PRIM
// Primary and Luft only line: Lrej and Prej
void vcf_lo_seg_generate_INFO_DVCF (VBlockVCFP vb, ZipDataLineVCF *dl)
{
    bool is_xstrand;
    PosType pos = vb->last_int(VCF_POS);

    LiftOverStatus ostatus = vcf_lo_get_liftover_coords (vb, pos, &dl->chrom_index[1], &dl->pos[1], &is_xstrand, &vb->pos_aln_i);
    unsigned oref_len;

    if (ostatus == LO_NO_MAPPING_IN_CHAIN) {
        vcf_lo_seg_rollback_and_reject (vb, LO_NO_MAPPING_IN_CHAIN, NULL);
        REJECT_MAPPING ("No alignment in the chain file");
    }
    
    else if (ostatus == LO_CHROM_NOT_IN_CHAIN) {
        vcf_lo_seg_rollback_and_reject (vb, LO_CHROM_NOT_IN_CHAIN, NULL);
        REJECT_MAPPING ("CHROM missing in chain file");
    }

    else if (ostatus == LO_CHROM_NOT_IN_PRIM_REF) {
        vcf_lo_seg_rollback_and_reject (vb, LO_CHROM_NOT_IN_PRIM_REF, NULL);
        REJECT_MAPPING ("CHROM missing in Primary reference file");
    }

    // REF & ALT - update ostatus based the relationship between REF, ALT, oREF and is_xstrand.
    if (LO_IS_REJECTED (ostatus = vcf_refalt_lift (vb, dl, is_xstrand))) {
        vcf_lo_seg_rollback_and_reject (vb, ostatus, NULL);
        return;
    }

    oref_len = LO_IS_OK_SWITCH (ostatus) ? vb->main_alt_len : vb->main_ref_len;

    // case is_xstrand: oPOS now points to the old anchor base - the last base in oREF - change it to point to the new anchor - 
    // the base before oREF (vb->new_ref contains the new anchor base)
    if (is_xstrand && LO_IS_OK_INDEL(ostatus)) 
        dl->pos[1] -= oref_len;

    // Seg oCHROM
    Context *ochrom_ctx = &vb->contexts[VCF_oCHROM];
    const char *ochrom=""; unsigned ochrom_len=0;

    if (dl->chrom_index[1] != WORD_INDEX_NONE)
        ctx_get_snip_by_zf_node_index (&ochrom_ctx->ol_nodes, &ochrom_ctx->ol_dict, dl->chrom_index[1], &ochrom, &ochrom_len);

    seg_by_ctx (vb, ochrom, ochrom_len, ochrom_ctx, ochrom_len);
    random_access_update_chrom ((VBlockP)vb, DC_LUFT, dl->chrom_index[1], ochrom, ochrom_len);

    // we add oPOS, oCHROM, oREF and oXSTRAND-  only in case ostatus is LO_IS_OK - they will be consumed by the info_luft_snip container
    // note: no need to seg VCF_COPYPOS explicitly here, as it is all_the_same handled in vcf_seg_initialize
    unsigned opos_len = str_int_len (dl->pos[1]);
    seg_pos_field ((VBlockP)vb, VCF_oPOS, VCF_oPOS, SPF_ZERO_IS_BAD, 0, 0, 0, dl->pos[1], opos_len);
    random_access_update_pos ((VBlockP)vb, DC_LUFT, VCF_oPOS);

    // no need to seg VCF_LIFT_REF here as it is all_the_same, just account for txt_len
    CTX(VCF_LIFT_REF)->txt_len += oref_len; 

    // special snip for oREFALT - vcf_piz_special_oREF will reconstruct it based on REF, ALT, oStatus and xstrand
    vcf_refalt_seg_other_REFALT (vb, VCF_oREFALT, ostatus, is_xstrand, 0); // 0 as not reconstructed by genounzip

    vb->last_index(VCF_oXSTRAND) = seg_by_did_i (vb, is_xstrand ? "X" : "-", 1, VCF_oXSTRAND, 1); // index used by vcf_seg_INFO_END

    // Add LUFT container - we modified the txt by adding these 4 fields to INFO/LUFT. We account for them now, and we will account for the INFO name etc in vcf_seg_info_field
    Context *luft_ctx = vcf_lo_seg_lo_snip (vb, info_luft_snip, info_luft_snip_len, dict_id_INFO_LUFT, 3);
    luft_ctx->last_txt_len = opos_len + ochrom_len + oref_len + 4; // used for rollback: 4 = 3 commas + xstrand (the "PRIM="/"LUFT=" is accounted for in vcf_seg_info_add_DVCF_to_InfoItems)
    vb->recon_size += luft_ctx->last_txt_len; // We added INFO/LUFT to the Primary coordinates reconstruction

    // Do the same for PRIM. 
    Context *prim_ctx = vcf_lo_seg_lo_snip (vb, info_prim_snip, info_prim_snip_len, dict_id_INFO_PRIM, 0);
    prim_ctx->last_txt_len = 0; // 0 as not added to recon_size (txt_len is used only for Primary coordinates)
    vb->recon_size_luft += vb->last_txt_len (VCF_POS) + vb->chrom_name_len + vb->main_ref_len + 4; // 4=3 commas + xstrand. We added INFO/PRIM to the Luft coordinates reconstruction

    // ostatus is OK 
    seg_by_did_i (vb, dvcf_status_names[ostatus], strlen (dvcf_status_names[ostatus]), VCF_oSTATUS, 0);
    vcf_set_ostatus (ostatus);
}

// -----------------------------------------------------------------------
// Segging LUFT / PRIM / LIFTREJT records of dual-coordinates files
// Liftover record: CHROM,POS,STRAND,REF
// Liftback record: CHROM,POS,STRAND,REF
// Liftrejt snip:   REJECTION_REASON
// -----------------------------------------------------------------------

static inline bool vcf_lo_is_same_seq (const char *seq1, unsigned seq1_len, const char *seq2, unsigned seq2_len, bool is_snp, bool is_xstrand) 
{
    if (seq1_len != seq2_len) return false;

    if (!is_xstrand) {
        for (PosType i=0; i < seq1_len ; i++)
            if (UPPER_CASE (seq1[i]) != UPPER_CASE (seq2[i])) return false;
    }
    else { 
        // if xstrand indel - don't compare anchor bases. Note: a non-SNP is_xstrand is always an indel bc we don't support xstrand with structural variants
        if (!is_snp) {
            seq1++; seq1_len--; 
            seq2++; seq2_len--;
        }
        for (PosType i=0; i < seq1_len ; i++)
            if (UPPER_CASE (seq1[i]) != UPPER_COMPLEM[(int)seq2[seq1_len-i-1]]) return false;
    }

    return true;
}

// When segging INFO/LUFT or INFO/PRIM of a dual coordinate file, set the last_ostatus based on the observed relationship 
// between the REF and ALT in the main VCF fields, and the REF in the INFO/LIFT[OVER|BACK] field
// returns true if status OK and false if REJECTED
static LiftOverStatus vcf_lo_seg_ostatus_from_LUFT_or_PRIM (VBlockVCFP vb, const char *field_name, bool is_xstrand,
                                                            const char *main_ref, unsigned main_ref_len,
                                                            const char *info_ref, unsigned info_ref_len,
                                                            const char *main_alt, unsigned main_alt_len,
                                                            unsigned *info_alt_len /* out */)
{
    LiftOverStatus ostatus;

    // short-circuit the majority case of a non-xstrand, no-base-change, bi allelic SNP
    if (main_ref_len == 1 && main_alt_len == 1 && !is_xstrand && main_ref[0] == info_ref[0]) {
        ostatus = LO_OK_REF_SAME_SNP;
        goto finalize;
    }

    bool is_snp = main_ref_len == 1 && (main_alt_len==1 || str_count_char (main_alt, main_alt_len, ',')*2+1 == main_alt_len);

    // check for structural variant
    if (memchr (main_alt, '<', main_alt_len)) 
        ostatus = LO_OK_REF_SAME_SV;

    // check for no-base-change 
    else if (vcf_lo_is_same_seq (main_ref, main_ref_len, info_ref, info_ref_len, is_snp, is_xstrand))
        ostatus = is_snp ? LO_OK_REF_SAME_SNP : LO_OK_REF_SAME_INDEL;

    // check for REF<>ALT switch 
    else if (vcf_lo_is_same_seq (main_alt, main_alt_len, info_ref, info_ref_len, is_snp, is_xstrand))
        ostatus = is_snp ? LO_OK_REF_ALT_SWITCH_SNP : LO_OK_REF_ALT_SWITCH_INDEL;

    else if (is_snp && main_alt_len == 1)
        ostatus = LO_OK_REF_NEW_SNP;

    else { 
        // this can happen only when genozipping a DVCF file produced by a different tool, like a more advanced version of genozip...
        WARNVCF ("FYI: Unsupported REF/ALT configuration: REF=%.*s ALT=%.*s INFO/DVCF/REF=%.*s INFO/DVCF/XSTRAND=%c", 
                  main_ref_len, main_ref, main_alt_len, main_alt, info_ref_len, info_ref, is_xstrand ? 'X' : '-');
        ostatus = LO_UNSUPPORTED_REFALT;
    }

finalize:        
    if (LO_IS_OK (ostatus)) {
        seg_by_did_i (vb, dvcf_status_names[ostatus], strlen (dvcf_status_names[ostatus]), VCF_oSTATUS, 0); // 0 bc doesn't reconstruct by default
        vcf_set_ostatus (ostatus); // we can set if OK, but only vcf_lo_seg_rollback_and_reject can set if reject

        *info_alt_len = (LO_IS_OK_SWITCH (ostatus) ? main_ref_len : main_alt_len);
    }
    else 
        vcf_lo_seg_rollback_and_reject (vb, ostatus, NULL); // segs *rej and oSTATUS

    return ostatus;
}

// Segging a Primary line, INFO/LUFT field ; or Luft line, INFO/PRIM field:
// parse the of LUFT/PRIM record, seg the contained fields, and seg the LUFT and PRIM snips
void vcf_lo_seg_INFO_LUFT_and_PRIM (VBlockVCFP vb, DictId dict_id, const char *value, int value_len)
{
    ASSVCF (!chain_is_loaded, "%s: --chain cannot be used with this file because it is already a dual-coordinates VCF file - it contains INFO/"INFO_LUFT" or INFO/"INFO_PRIM" fields", txt_name);
    ASSVCF (txt_file->coords, "Found an %s subfield, but this is not a dual-coordinates VCF because it is missing \"" HK_DC_PRIMARY "\" or \"" HK_DC_LUFT "\" in the VCF header", dis_dict_id_ex (dict_id, true).s);

    Coords coord = (dict_id.num == dict_id_INFO_LUFT ? DC_PRIMARY : DC_LUFT);
    const char *field_name = dis_dict_id_ex (dict_id, true).s; // hopefully the returned-by-value structure stays in scope despite not assigning it...

    ASSVCF (vb->line_coords == coord, "Found an %s subfield, not expecting because this is a %s line", field_name, coords_name (vb->line_coords));

    ZipDataLineVCF *dl = DATA_LINE (vb->line_i);
        
    // parse and verify PRIM or LUFT record
    str_split (value, value_len, NUM_IL_FIELDS + 1, ',', str, false); // possibly oSTATUS in addition from --show-ostatus which we shall ignore
    ASSVCF (n_strs == 4, "Invalid %s field: \"%.*s\"", field_name, value_len, value);

    const char *info_chrom  = strs[IL_CHROM]; 
    unsigned info_chrom_len = str_lens[IL_CHROM];
    const char *info_ref    = strs[IL_REF];
    unsigned info_ref_len   = str_lens[IL_REF];
    unsigned info_alt_len;
    int64_t info_pos_value;
    unsigned info_pos_len   = str_lens[IL_POS];
    ASSVCF (str_get_int_range64 (strs[IL_POS], str_lens[IL_POS], 0, MAX_POS, &info_pos_value), "Invalid POS value in %s: \"%.*s\"", field_name, str_lens[IL_POS], strs[IL_POS]);
    ASSVCF ((strs[IL_XSTRAND][0] == 'X' || strs[IL_XSTRAND][0] == '-') && str_lens[IL_XSTRAND] ==1,
            "%s has an invalid XSTRAND=\"%.*s\" - expected \"-\" or \"X\"", field_name, str_lens[IL_XSTRAND], strs[IL_XSTRAND]);
    bool is_xstrand = (strs[IL_XSTRAND][0] == 'X');

    LiftOverStatus ostatus = 
        vcf_lo_seg_ostatus_from_LUFT_or_PRIM (vb, field_name, is_xstrand,
                                              vb->main_refalt, vb->main_ref_len, 
                                              info_ref, info_ref_len, 
                                              vb->main_refalt + vb->main_ref_len + 1, vb->main_alt_len,
                                              &info_alt_len); 
                                              
    if (LO_IS_REJECTED (ostatus)) return; // rolled back and segged *rej instead, due to rejection

    // Seg other (than vb->line_coords) coords' CHROM
    dl->chrom_index[SEL(1,0)] = seg_by_did_i (vb, info_chrom, info_chrom_len, SEL(VCF_oCHROM, VCF_CHROM), info_chrom_len);
    random_access_update_chrom ((VBlockP)vb, SEL(DC_LUFT, DC_PRIMARY), dl->chrom_index[SEL(1,0)], info_chrom, info_chrom_len);
    
    // Seg other coord's POS
    // note: no need to seg VCF_COPYPOS explicitly here, as it is all_the_same handled in vcf_seg_initialize
    dl->pos[SEL(1,0)] = seg_pos_field ((VBlockP)vb, SEL(VCF_oPOS, VCF_POS), SEL(VCF_oPOS, VCF_POS), SPF_ZERO_IS_BAD, 0, 0, 0, info_pos_value, info_pos_len);
    if (dl->pos[SEL(1,0)])random_access_update_pos ((VBlockP)vb, SEL(DC_LUFT, DC_PRIMARY), SEL (VCF_oPOS, VCF_POS));

    // Seg other coord's REFALT - this will be a SPECIAL that reconstructs from this coords REFALT + INFO_REF + XSTRAND
    if (ostatus == LO_OK_REF_NEW_SNP || (LO_IS_OK_INDEL (ostatus) && is_xstrand))
        vb->new_ref = info_ref[0]; // new anchor or new REF SNP (consumed by vcf_refalt_seg_other_REFALT)

    vcf_refalt_seg_other_REFALT (vb, SEL(VCF_oREFALT, VCF_REFALT), ostatus, is_xstrand, SEL (0, info_ref_len + 2 + info_alt_len)); // account for REFALT + 2 tabs (primary coords) if we are segging it now (but not for oREFALT)
    
    // "Seg" other coord's REF to appear in the this PRIM/LUFT - a SPECIAL that reconstructs the other REFALT, and discards the ALT
    // note: no need to actually seg here as it is all_the_same handled in vcf_seg_initialize, just account for txt_len
    CTX(VCF_LIFT_REF)->txt_len += SEL (info_ref_len, 0); // we account for oREF (to be shown in INFO/LUFT in default reconstruction). 

    // an indel with REF<>ALT switch can change the size of the oREF field in INFO/LUFT
    if (ostatus == LO_OK_REF_ALT_SWITCH_INDEL) 
        vb->recon_size += (int)vb->main_ref_len - (int)info_ref_len;

    // Seg the XSTRAND (same for both coordinates)
    seg_by_did_i (vb, (is_xstrand ? "X" : "-"), 1, VCF_oXSTRAND, 1);

    // LUFT and PRIM container snips (INFO container filter will determine which is reconstructed)
    Context *luft_ctx = vcf_lo_seg_lo_snip (vb, info_luft_snip, info_luft_snip_len, dict_id_INFO_LUFT, 3); // account for 3 commas
    Context *prim_ctx = vcf_lo_seg_lo_snip (vb, info_prim_snip, info_prim_snip_len, dict_id_INFO_PRIM, 0); // 0 as this container is not reconstructed in Primary coords
    
    if (dict_id.num == dict_id_INFO_LUFT) luft_ctx->last_txt_len = value_len;
    else                                  prim_ctx->last_txt_len = value_len;
}

// Segging a Primary or Luft dual-coordinates VCF file, INFO/Prej or Lrej field:
void vcf_lo_seg_INFO_REJX (VBlockVCFP vb, DictId dict_id, const char *value, int value_len)
{
    ASSVCF (!chain_is_loaded, "%s: --chain cannot be used with this file because it is already a dual-coordinates VCF file - it contains INFO/"INFO_LREJ" or INFO/"INFO_PREJ" fields", txt_name);
    ASSVCF (txt_file->coords, "Found an %s subfield, but this is not a dual-coordinates VCF because it is missing \"" HK_DC_PRIMARY "\" or \"" HK_DC_LUFT "\" in the VCF header", dis_dict_id_ex (dict_id, true).s);

    Coords coord = (dict_id.num == dict_id_INFO_LREJ ? DC_PRIMARY : DC_LUFT);
    ASSVCF (vb->line_coords == coord, "Found an %s subfield, not expecting because this is a %s line", dis_dict_id_ex (dict_id, true).s, coords_name (vb->line_coords));

    // value may be with or without a string eg "INFO/AC" ; "NO_MAPPING". We take the value up to the comma.
    const char *slash = memchr (value, '/', value_len);

    // get status from string
    SAFE_NUL (slash ? slash : &value[value_len]);

    unsigned i; for (i=LO_REJECTED; i < NUM_LO_STATUSES; i++)
        if (!strcmp (value, dvcf_status_names[i])) { // found
            vcf_set_ostatus (i);
            break;
        }
    
    if (i == NUM_LO_STATUSES)
        vcf_set_ostatus (LO_REJECTED); // fallback

    SAFE_RESTORE;

    seg_by_did_i (vb, value, value_len, VCF_oSTATUS, value_len); // note: word_index >= NUM_LO_STATUSES if unrecognized string, that's ok
    vcf_lo_seg_REJX_do (vb, 0);
}

// -------------------
// PIZ: reconstruction
// -------------------

SPECIAL_RECONSTRUCTOR (vcf_piz_special_COPYSTAT)
{
    if (!reconstruct) return false; // no new_value
    
    // note: we can't use last_txt because oStatus is not reconstructed
    ctx_get_snip_by_word_index (&vb->contexts[VCF_oSTATUS], vb->last_index(VCF_oSTATUS), &snip, &snip_len); 
    RECONSTRUCT (snip, snip_len); // works for unrecognized reject statuses too
    return false; // new_value was not set
}

// Callback at the end of reconstructing a VCF line of a dual-coordinates file, drops lines in some cases
void vcf_lo_piz_TOPLEVEL_cb_filter_line (VBlock *vb)
{
    Coords line_coord = vb->last_index (VCF_COORDS);

    // Case 1: lines can be PRIMARY, LUFT or BOTH. vb rendering can be PRIMARY or LUFT. Drop line if there is a mismatch.
    if (!(vb->vb_coords & line_coord))  // note: for sorted files, rejects are already excluded from the reconstruction plan in linesorter_merge_vb_do
        vb->drop_curr_line = "no_coords";

    // Case 2: drop ##primary_only / ##luft_only header lines in --no-header or --header-one 
    else if (vb->is_rejects_vb && (flag.no_header || flag.header_one)) 
        vb->drop_curr_line = "no_header_rejects";
}

// ----------------------------------------------------------------
// Rejects file stuff - file contains data in the native txt format
// ----------------------------------------------------------------

// ZIP: called when inspecting the txtheader to add header data, and after each VB to add rejected line
void vcf_lo_append_rejects_file (VBlockP vb, Coords coord)
{
    ASSERTNOTNULL (z_file);

    // create rejects file if not already open
    if (!z_file->rejects_file[coord-1]) {
        z_file->rejects_file_name[coord-1] = malloc (strlen (z_file->name) + 20);
        sprintf (z_file->rejects_file_name[coord-1], "%s.%s_ONLY%s", z_file->name, coords_name (coord), file_plain_ext_by_dt (z_file->data_type));
        
        z_file->rejects_file[coord-1] = fopen (z_file->rejects_file_name[coord-1], "wb+");
        ASSERT (z_file->rejects_file[coord-1], "fopen() failed to open %s: %s", z_file->rejects_file_name[coord-1], strerror (errno));
    }

    ASSERT0 (z_file->rejects_file[coord-1], "liftover rejects file is not open");

    fwrite (vb->lo_rejects[coord-1].data, 1, vb->lo_rejects[coord-1].len, z_file->rejects_file[coord-1]);
    z_file->rejects_disk_size[coord-1] += vb->lo_rejects[coord-1].len;

    buf_free (&vb->lo_rejects[coord-1]);
}

void vcf_liftover_display_lift_report (void)
{
    char report_fn[strlen (z_name) + 50];
    sprintf (report_fn, "%s.rejects.txt", z_name);

    if (z_file->rejects_report.len)
        buf_dump_to_file (report_fn, &z_file->rejects_report, 1, false, false, false);

    if (!flag.quiet)
        iprintf ("\nCongratulations! %s is a dual-coordinate VCF (DVCF). You can now:\n"
                "1. Render the data in Primary coordinates: genocat %s\n"
                "2. Render the data in Luft coordinates: genocat --luft %s\n"
                "3. See line by line variant status: genocat --show-dvcf %s (also works in combination with --luft)\n"
                "4. See the lift rejects report (only exists in case of rejects): see %s\n"
                "5. Learn more about DVCF: https://genozip.com/dvcf.html\n\n", 
                z_name, z_name, z_name, z_name, report_fn);
}
