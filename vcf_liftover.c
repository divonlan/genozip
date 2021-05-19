// ------------------------------------------------------------------
//   vcf_liftover.c
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// Dual-coordinates genozip/genocat file flow:
//
// genozip -C source.vcf       --> dual-coord.genozip - component_i=0 contains all lines, each with LIFTOVER or LIFTREJT 
//                                                                    LIFTREFD are also written to the rejects file
//                                                      component_i=1 is the rejects file containing LIFTREFD lines (again)
//
// genocat dual-crd.genozip    --> primary.vcf -        component_i=0 is reconstructed - lines have LIFTOVER or LIFTREJT 
//                                                                    lines are IN ORDER
//                                                      component_i=1 is skipped.
//
// genozip primary.vcf         --> dual-coord.genozip - component_i=0 contains all lines, each with LIFTOVER or LIFTREJT 
//                                                                    LIFTREFD are also written to the rejects file
//                                                      component_i=1 is the rejects file containing LIFTREFD lines (again)
//
// genocat -v dual-crd.genozip --> luft.vcf -           component_i=1 with LIFTREJT is reconstructed first and becomes part of the header
//                                                      component_i=0 is reconstructed - dropping LIFTREJT lines and lifting over LIFTOVER->LIFTBACK 
//                                                                    lines sorted by writer according the the RECON_PLAN section 
//
// genozip luft.vcf            --> dual-coord.genozip - LIFTREJT header lines are sent to vblock_i=1 via unconsumed_txt
//                                                      component_i=0 contains all lines first all LIFTREJT lines followed by all LIFTOVER lines  
//                                                                    LIFTREFD lines are also written to the rejects file
//                                                      component_i=1 is the rejects file containing LIFTREFD lines (again)
//
// genocat dual-crd.genozip    --> primary.vcf          component_i=0 is reconstructed - first LIFTREJT then LIFTOVER lines
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

const char *liftover_status_names[] = LIFTOVER_STATUS_NAMES;
const LuftTransLateProp ltrans_props[NUM_VCF_TRANS] = LUFT_TRANS_PROPS;

#define WORD_INDEX_MISSING -2 // a value in liftover_chrom2chainsrc

// constant snips initialized at the beginning of the first file zip
#define LO_SNIP_LEN 200
static char liftover_snip[LO_SNIP_LEN]={}, liftback_snip[LO_SNIP_LEN]={}, rejtover_snip[LO_SNIP_LEN]={}, rejtback_snip[LO_SNIP_LEN]={};
static unsigned liftover_snip_len=LO_SNIP_LEN, liftback_snip_len=LO_SNIP_LEN, rejtover_snip_len=LO_SNIP_LEN, rejtback_snip_len=LO_SNIP_LEN;
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
    if (chain_is_loaded && z_file->num_txt_components_so_far == 1) // copy only once
        chain_copy_dst_contigs_to_z_file (DTFZ (ochrom));

    // prepare (constant) snips. note: we need these snips both for --chain and when zipping dual-coord files
    SmallContainer con = {
        .repeats     = 1,
        .nitems_lo   = NUM_IL_FIELDS,
        .drop_final_item_sep = true,
        .items       = { [IL_CHROM  ]={ .dict_id = (DictId)dict_id_fields[VCF_oCHROM],   .seperator = ","  },
                         [IL_POS    ]={ .dict_id = (DictId)dict_id_fields[VCF_oPOS],     .seperator = ","  },
                         [IL_REF    ]={ .dict_id = (DictId)dict_id_fields[VCF_LIFT_REF], .seperator = ","  },
                         [IL_XSTRAND]={ .dict_id = (DictId)dict_id_fields[VCF_oXSTRAND], .seperator = ","  } } };
    container_prepare_snip ((Container*)&con, 0, 0, liftover_snip, &liftover_snip_len);

    con = (SmallContainer){
        .repeats     = 1,
        .nitems_lo   = NUM_IL_FIELDS,
        .drop_final_item_sep = true,
        .items       = { [IL_CHROM  ]={ .dict_id = (DictId)dict_id_fields[VCF_CHROM],    .seperator = ","  },
                         [IL_POS    ]={ .dict_id = (DictId)dict_id_fields[VCF_COPYPOS],  .seperator = ","  },
                         [IL_REF    ]={ .dict_id = (DictId)dict_id_fields[VCF_LIFT_REF], .seperator = ","  },
                         [IL_XSTRAND]={ .dict_id = (DictId)dict_id_fields[VCF_oXSTRAND], .seperator = ","  } } };
    container_prepare_snip ((Container*)&con, 0, 0, liftback_snip, &liftback_snip_len);

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
    container_prepare_snip ((Container*)&con, 0, 0, rejtover_snip, &rejtover_snip_len);

    con = (SmallContainer){
        .repeats     = 1,
        .nitems_lo   = 4,
        .items       = { { .dict_id = (DictId)dict_id_fields[VCF_COPYSTAT] },
                         { .dict_id = (DictId)dict_id_fields[VCF_CHROM]    },
                         { .dict_id = (DictId)dict_id_fields[VCF_COPYPOS], .seperator = { CI_TRANS_NOR } }, // rather than segging "", we don't reconstruct. so we don't break the "all_the_same" of COPYPOS
                         { .dict_id = (DictId)dict_id_fields[VCF_REFALT]   } } };
    container_prepare_snip ((Container*)&con, 0, 0, rejtback_snip, &rejtback_snip_len);
}

// get liftover translator ID. For contexts for which there is an ##INFO or ##FORMAT field in the VCF header,
// these are entered to the context by vcf_header. Fields appearing in the file are queries by Seg. 
TranslatorId vcf_lo_luft_trans_id (DictId dict_id, char number)
{
    // note: in VCF 4.1 many files specified A,G and R fields as Number=. ; 'R' was introduced in 4.2. We set some common fields
    // regardless of their "Number" in the header.
    
    if (number == 'R' || dict_id.num == dict_id_FORMAT_AD || dict_id.num == dict_id_FORMAT_F2R1 || dict_id.num == dict_id_FORMAT_F1R2 ||
        dict_id.num == dict_id_INFO_DP_HIST || dict_id.num == dict_id_INFO_GQ_HIST) 
        return VCF2VCF_R;

    else if (number == 'G' || dict_id.num == dict_id_FORMAT_GL || dict_id.num == dict_id_FORMAT_PL || dict_id.num == dict_id_FORMAT_PRI ||
             dict_id.num == dict_id_FORMAT_GP)
        return VCF2VCF_G;

    else if (dict_id.num == dict_id_FORMAT_GT)       return VCF2VCF_GT;
    else if (dict_id.num == dict_id_INFO_END )       return VCF2VCF_END;
    else if (dict_id.num == dict_id_INFO_AC  )       return VCF2VCF_A_AN;
    else if (dict_id.num == dict_id_INFO_AF 
    ||       dict_id.num == dict_id_FORMAT_AF  
    ||       (dict_id.id[0] == ('A' | 0xc0) && dict_id.id[1] == 'F' && dict_id.id[2] == '_') // INFO/AF_*
    ||       (dict_id_is_vcf_info_sf (dict_id) && dict_id.id[3] == '_' && dict_id.id[4] == 'A' && dict_id.id[5] == 'F' && !dict_id.id[6])) // INFO/???_AF
          return VCF2VCF_A_1; // eg 0.150 -> 0.850 (upon REF<>ALT SWITCH)
    else if (dict_id.num == dict_id_FORMAT_SAC)      return VCF2VCF_R2;   // eg 25,1,6,2 --> 6,2,25,1 (upon REF<>ALT SWITCH)
    else if (dict_id.num == dict_id_INFO_BaseCounts) return VCF2VCF_XREV; // eg 10,6,3,4 -> 4,3,6,10  (upon XSTRAND)
    
    return TRANS_ID_NONE;
}

// SEG: map coordinates primary->luft. In case of failure, sets dst_contig_index, dst_1pos, xstrand to 0.
LiftOverStatus vcf_lo_get_liftover_coords (VBlockVCFP vb, WordIndex *dst_contig_index, PosType *dst_1pos, bool *xstrand, uint32_t *aln_i) // out
{
    // extend vb->liftover if needed
    if (vb->chrom_node_index >= vb->liftover.len) { 
        buf_alloc (vb, &vb->liftover, 0, MAX (vb->chrom_node_index+1, 100), WordIndex, 0, "liftover");

        // initialize new entries allocated to WORD_INDEX_MISSING
        uint32_t size = vb->liftover.size / sizeof (WordIndex);
        for (uint32_t i=vb->liftover.param; i < size; i++)
            *ENT (WordIndex, vb->liftover, i) = WORD_INDEX_MISSING;

        vb->liftover.param = size; // param holds the number of entries initialized >= len
        vb->liftover.len   = vb->chrom_node_index + 1;
    }

    ARRAY (WordIndex, map, vb->liftover);

    // if chrom is not yet mapped to src_contig, map it now:
    // map from chrom_node_index (not word index!) to entry in chain_index
    if (map[vb->chrom_node_index] == WORD_INDEX_MISSING) 
        map[vb->chrom_node_index] = chain_get_src_contig_index (vb->chrom_name, vb->chrom_name_len);

    if (map[vb->chrom_node_index] == WORD_INDEX_NONE) {
        if (dst_contig_index) *dst_contig_index = WORD_INDEX_NONE;
        if (dst_1pos) *dst_1pos = 0;
        if (xstrand) *xstrand = false;
        return LO_NO_CHROM;
    }
    else {
        bool mapping_ok = chain_get_liftover_coords (map[vb->chrom_node_index], vb->last_int(VCF_POS), 
                                                     dst_contig_index, dst_1pos, xstrand, aln_i); // if failure, sets output to 0.
        return mapping_ok ? LO_OK : LO_NO_MAPPING;
    }
}                                

// lift-back Luft values to Primary (when segging a Luft file)
// returns true if successfully translated, false if left as-is
// Notes: 1. all translators are symetrical - running on primary results in luft and vice versa
// 2. AC and END translators might result in bigger or smaller text - we handle this
bool vcf_lo_seg_lift_back_to_primary (VBlockVCFP vb, ContextP ctx, 
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
            vb->vb_data_size -= (int)this_value_len - (int)*modified_snip_len;

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
static void vcf_lo_seg_REJTXXXX_do (VBlockVCFP vb, unsigned add_bytes)
{
    ZipDataLineVCF *dl = DATA_LINE (vb->line_i);
    
    dl->chrom_index[SEL(1,0)] = WORD_INDEX_NONE;
    dl->pos[SEL(1,0)] = 0;

    // whether user views the file primary or luft coordinates, and whether this line is primary-only or luft-only:
    // each rejected line consumes all contexts exactly once, either in the main field, or in INFO/REJTXXXX
    // note: rejected lines are still reconstructed in their unsupported coordinates to consume all contexts, and then dropped.
    vcf_lo_seg_lo_snip (vb, rejtover_snip, rejtover_snip_len, dict_id_INFO_REJTOVER, 0); 
    vcf_lo_seg_lo_snip (vb, rejtback_snip, rejtback_snip_len, dict_id_INFO_REJTBACK, 0);
    
    // case: primary-only line: CHROM, POS and REFALT were segged in the main fields, trivial oCHROM, oPOS, oREFLT segged here
    // case: luft-only line: oCHROM, oPOS and oREFALT were segged in the main fields, trivial CHROM, POS, REFLT segged here
    seg_by_did_i (vb, "", 0, SEL (VCF_oCHROM,  VCF_CHROM),  0);  // 0 for all of these as these fields are not reconstructed by default (genounzip) (oXXXX never reconstructed by default, CHROM, POS, REFALT not reconstruced when this variant has no primary coords)
    seg_by_did_i (vb, "", 0, SEL (VCF_oREFALT, VCF_REFALT), 0);
    seg_by_did_i (vb, "", 0, SEL (VCF_oPOS,    VCF_POS),    0);

    vb->contexts[VCF_COPYSTAT].txt_len += add_bytes; // we don't need to seg COPYSTAT because it is all_the_same. Just account for add_bytes.
    vb->vb_data_size += add_bytes; // we modified the vb in the default (genounzip) reconstruction, by adding this string (the "REJTXXXX=" is accounted for in vcf_seg_info_add_LIFTXXXX_items)
}

#define NUM_ROLLBACK_CTXS 6
static const DidIType line_rollback_did_i[2][NUM_ROLLBACK_CTXS] = { 
    { VCF_oSTATUS, VCF_oCHROM, VCF_oPOS, VCF_oREFALT, VCF_LIFT_REF, VCF_oXSTRAND },
    { VCF_oSTATUS, VCF_CHROM,  VCF_POS,  VCF_REFALT,  VCF_LIFT_REF, VCF_oXSTRAND }
};

// Actually, reject this line, after already segging to o* fields in vcf_lo_seg_generate_INFO_LIFTXXXX
void vcf_lo_seg_rollback_and_reject (VBlockVCFP vb, LiftOverStatus ostatus, Context *ctx /* INFO or FORMAT context, or NULL */)
{
    if (LO_IS_REJECTED (last_ostatus)) return; // already rejected

    // unaccount for INFO/LIFTOVER. No need to unaccount for LIFTBACK as it is not reconstructed by default so not included in vb_data_size. 
    Context *lo_ctx; 
    if (ctx_encountered_in_line (vb, dict_id_INFO_LIFTOVER, &lo_ctx))
        vb->vb_data_size -= lo_ctx->last_txt_len;

    for (unsigned i=0; i < NUM_ROLLBACK_CTXS; i++) 
        seg_rollback ((VBlockP)vb, &vb->contexts[line_rollback_did_i[SEL(0,1)][i]]);

    seg_rollback ((VBlockP)vb, ctx_get_ctx (vb, dict_id_INFO_LIFTOVER));
    seg_rollback ((VBlockP)vb, ctx_get_ctx (vb, dict_id_INFO_LIFTBACK));

    // seg REJT**** and empty values for the fields of the rejected coord
    char str[MAX_VCF_ID_LEN + 20];
    const char *id = ctx ? vcf_header_get_VCF_ID_by_dict_id (ctx->dict_id, false) : NULL;
    sprintf (str, "%s%s%s", liftover_status_names[ostatus], ctx ? "/" : "", 
             id ? id : ctx ? dis_dict_id (ctx->dict_id).s : "");
    unsigned str_len = strlen (str);

    vcf_lo_seg_REJTXXXX_do (vb, str_len); 

    // seg oSTATUS
    seg_by_did_i (vb, str, str_len, VCF_oSTATUS, 0);
    vcf_set_ostatus (ostatus); // ONLY this place can set a rejection ostatus, or else it will fail the if at the top of this function
}

void vcf_lo_set_rollback_point (VBlockVCFP vb)
{
    for (unsigned i=0; i < NUM_ROLLBACK_CTXS; i++) 
        seg_create_rollback_point (&vb->contexts[line_rollback_did_i[SEL(0,1)][i]]);

    seg_create_rollback_point (ctx_get_ctx (vb, dict_id_INFO_LIFTOVER));
    seg_create_rollback_point (ctx_get_ctx (vb, dict_id_INFO_LIFTBACK));
}

// --chain: Create:
// Dual coordinates line: LIFTOVER and LIFTBACK
// Primary only line:     REJTOVER
// Luft only line:        REJTBACK
void vcf_lo_seg_generate_INFO_LIFTXXXX (VBlockVCFP vb, ZipDataLineVCF *dl)
{
    bool xstrand;
    LiftOverStatus ostatus = vcf_lo_get_liftover_coords (vb, &dl->chrom_index[1], &dl->pos[1], &xstrand, &vb->pos_aln_i);
    unsigned oref_len;

    if (LO_IS_OK (ostatus)) // now update ostatus the relationship between REF, ALT, oREF and xstrand.
         ostatus = vcf_refalt_check_oref (vb, dl, xstrand, &oref_len);

    if (LO_IS_REJECTED (ostatus)) {
        vcf_lo_seg_rollback_and_reject (vb, ostatus, NULL);
        return;
    }

    // we add oPOS, oCHROM, oREF and oXSTRAND-  only in case ostatus is LO_IS_OK - they will be consumed by the liftover_snip container
    unsigned opos_len = str_int_len (dl->pos[1]);
    seg_pos_field ((VBlockP)vb, VCF_oPOS, VCF_oPOS, false, true, 0, 0, 0, dl->pos[1], opos_len);

    // special snip for POS reconstruction in INFO/LIFTBACK (handling the possibility of INFO/END)
    // note: no need to seg VCF_COPYPOS explicitly here, as it is all_the_same handled in vcf_seg_initialize

    // hacky shortcut to segging chrom as we know its index but not its name
    Context *chrom_ctx = &vb->contexts[VCF_oCHROM];
    buf_alloc (vb, &chrom_ctx->b250, 1, vb->lines.len, WordIndex, CTX_GROWTH, "contexts->b250");
    NEXTENT (WordIndex, chrom_ctx->b250) = dl->chrom_index[1];
    ctx_increment_count ((VBlockP)vb, chrom_ctx, dl->chrom_index[1]);
  
    unsigned ochrom_len = ENT (CtxNode, vb->contexts[VCF_oCHROM].ol_nodes, dl->chrom_index[1])->snip_len;
    vb->contexts[VCF_oCHROM].txt_len += ochrom_len;

    // special snip for oREFALT - vcf_piz_special_oREF will reconstruct it based on REF, ALT, oStatus and xstrand
    vb->contexts[VCF_LIFT_REF].txt_len += oref_len; // no need to seg VCF_LIFT_REF here as it is all_the_same, just account for txt_len
    seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_other_REFALT }), 2, VCF_oREFALT, 0); // 0 as not reconstructed by genounzip

    seg_by_did_i (vb, xstrand ? "X" : "-", 1, VCF_oXSTRAND, 1);

    // Add liftover container - we modified the txt by adding these 4 fields to INFO/LIFTOVER. We account for them now, and we will account for the INFO name etc in vcf_seg_info_field
    Context *lo_ctx = vcf_lo_seg_lo_snip (vb, liftover_snip, liftover_snip_len, dict_id_INFO_LIFTOVER, 3);
    lo_ctx->last_txt_len = opos_len + ochrom_len + 5; // 5 = 3 commas + oref + xstrand (the "LIFTXXXX=" is accounted for in vcf_seg_info_add_LIFTXXXX_items)
    vb->vb_data_size += lo_ctx->last_txt_len;         

    // Do the same for LIFTBACK. No change in vb_data_size because it is not reconstructed by default (ie in genounzip)
    Context *lb_ctx = vcf_lo_seg_lo_snip (vb, liftback_snip, liftback_snip_len, dict_id_INFO_LIFTBACK, 0);
    lb_ctx->last_txt_len = 0; // 0 as not added to vb_data_size

    // ostatus is OK 
    seg_by_did_i (vb, liftover_status_names[ostatus], strlen (liftover_status_names[ostatus]), VCF_oSTATUS, 0);
    vcf_set_ostatus (ostatus);
}

// -----------------------------------------------------------------------
// Segging LIFTOVER / LIFTBACK / LIFTREJT records of dual-coordinates files
// Liftover record: CHROM,POS,STRAND,REF
// Liftback record: CHROM,POS,STRAND,REF
// Liftrejt snip:   REJECTION_REASON
// -----------------------------------------------------------------------

// When segging INFO/LIFTOVER or INFO/LIFTBACK of a dual coordinate file, set the last_ostatus based on the observed relationship 
// between the REF and ALT in the main VCF fields, and the REF in the INFO/LIFT[OVER|BACK] field
// returns true if status OK and false if REJECTED
static bool vcf_lo_seg_ostatus_from_LIFTXXXX (VBlockVCFP vb, const char *field_name, bool is_xstrand,
                                              const char *main_ref, unsigned main_ref_len,
                                              const char *info_ref, unsigned info_ref_len,
                                              const char *main_alt, unsigned main_alt_len,
                                              unsigned *info_alt_len /* out */)
{
    LiftOverStatus ostatus;
    #define RESULT(o) do {ostatus=(o); goto done; } while(0)

    // if is_xstrand, expecting REF and oREF to be of length 1 (or else this line would have been rejected)
    if (!is_xstrand && main_ref_len > 1) RESULT (LO_REF_LONG_XSTRAND);
    if (!is_xstrand && info_ref_len > 1) RESULT (LO_ALT_LONG_XSTRAND);
    
    // if REF or oREF are longer than 1 (and !is_xstrand as verified above), then they must be identical
    if (main_ref_len > 1 || info_ref_len > 1) {

        if (main_ref_len != info_ref_len || memcmp (main_ref, info_ref, main_ref_len))
            RESULT (LO_REF_LONG_CHANGE);

        RESULT (LO_OK_REF_SAME);
    }

    // case: REF and oREF are both of length 1
    else switch (vcf_refalt_ref_equals_ref2_or_alt (main_ref[0], info_ref[0], main_alt, 1, is_xstrand)) {
        case EQUALS_ALT     : RESULT (LO_OK_REF_ALT_SWTCH);
        case EQUALS_REF2    : RESULT (LO_OK_REF_SAME);
        case EQUALS_MISSING : RESULT (LO_OK_REF_SAME); // both are '.'
        case EQUALS_NEITHER : RESULT (LO_REF_CHNG_NOT_ALT);
    }

    ABORT0_R ("Should never reach here");

done:
    if (LO_IS_OK (ostatus)) {
        seg_by_did_i (vb, liftover_status_names[ostatus], strlen (liftover_status_names[ostatus]), VCF_oSTATUS, 0); // 0 bc doesn't reconstruct by default
        vcf_set_ostatus (ostatus); // we can set if OK, but only vcf_lo_seg_rollback_and_reject can set if reject

        *info_alt_len = (ostatus == LO_OK_REF_ALT_SWTCH ? main_ref_len : main_alt_len);
        return true;
    }
    else {
        vcf_lo_seg_rollback_and_reject (vb, ostatus, NULL); // segs REJTXXXX and oSTATUS
        return false;
    }

    #undef RESULT
}

// Segging a Primary line, INFO/LIFTOVER field ; or Luft line, INFO/LIFTBACK field:
// parse the of LIFTOVER record, seg the contained fields, and add the liftover and liftback snips
void vcf_lo_seg_INFO_LIFTXXXX (VBlockVCFP vb, DictId dict_id, const char *value, int value_len)
{
    ASSVCF (!chain_is_loaded, "%s: --chain cannot be used with this file because it is already a dual-coordinates VCF file", txt_name);
    ASSVCF (txt_file->coords, "Found an %s subfield, but this is not a dual-coordinates VCF because it is missing \"" HK_DC_PRIMARY "\" or \"" HK_DC_LUFT "\" in the VCF header", dis_dict_id_ex (dict_id, true).s);

    Coords coord = (dict_id.num == dict_id_INFO_LIFTOVER ? DC_PRIMARY : DC_LUFT);
    const char *field_name = dis_dict_id_ex (dict_id, true).s; // hopefully the returned-by-value structure stays in scope despite not assigning it...

    ASSVCF (vb->line_coords == coord, "Found an %s subfield, not expecting because this is a %s line", field_name, coords_names[vb->line_coords]);

    ZipDataLineVCF *dl = DATA_LINE (vb->line_i);
        
    // parse and verify LIFTXXXX record
    const char *strs[NUM_IL_FIELDS]; unsigned str_lens[NUM_IL_FIELDS]; 
    ASSVCF (str_split (value, value_len, NUM_IL_FIELDS, ',', strs, str_lens, true, 0), "Invalid %s field: \"%.*s\"", field_name, value_len, value);
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

    if (!vcf_lo_seg_ostatus_from_LIFTXXXX (vb, field_name, is_xstrand,
                                           vb->main_refalt, vb->main_ref_len, 
                                           info_ref, info_ref_len, 
                                           vb->main_refalt + vb->main_ref_len + 1, vb->main_alt_len,
                                           &info_alt_len)) return; // rolled back and segged REJTXXXX instead, due to rejection

    // Seg other (than vb->line_coords) coords' CHROM
    dl->chrom_index[SEL(1,0)] = seg_by_did_i (vb, info_chrom, info_chrom_len, SEL(VCF_oCHROM, VCF_CHROM), info_chrom_len);
    random_access_update_chrom ((VBlockP)vb, SEL(DC_LUFT, DC_PRIMARY), dl->chrom_index[SEL(1,0)], info_chrom, info_chrom_len);
    
    // Seg other coord's POS
    dl->pos[SEL(1,0)] = seg_pos_field ((VBlockP)vb, SEL(VCF_oPOS, VCF_POS), SEL(VCF_oPOS, VCF_POS), false, true, 0, 0, 0, info_pos_value, info_pos_len);
    random_access_update_pos ((VBlockP)vb, SEL(DC_LUFT, DC_PRIMARY), SEL (VCF_oPOS, VCF_POS));

    // special snip for POS reconstruction in INFO/LIFTBACK (handling the possibility of INFO/END)
    // note: no need to seg VCF_COPYPOS explicitly here, as it is all_the_same handled in vcf_seg_initialize

    // Seg other coord's REFALT - this will be a SPECIAL that reconstructs from this coords REFALT + INFO_REF + XSTRAND
    seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_other_REFALT }), 2, SEL(VCF_oREFALT, VCF_REFALT), SEL (0, info_ref_len + 2 + info_alt_len)); // account for REFALT + 2 tabs (primary coords) if we are segging it now (but not for oREFALT)

    // "Seg" other coord's REF to appear in the this LIFTXXXX - a SPECIAL that reconstructs the other REFALT, and discards the ALT
    // note: no need to actually seg here as it is all_the_same handled in vcf_seg_initialize, just account for txt_len
    vb->contexts[VCF_LIFT_REF].txt_len += SEL (info_ref_len, 0); // we account for oREF (to be shown in INFO/LIFTOVER in default reconstruction). 

    // Seg the XSTRAND (same for both coordinates)
    seg_by_did_i (vb, (is_xstrand ? "X" : "-"), 1, VCF_oXSTRAND, 1);

    // LIFTOVER and LIFTBACK container snips (INFO container filter will determine which is reconstructed)
    Context *lo_ctx = vcf_lo_seg_lo_snip (vb, liftover_snip, liftover_snip_len, dict_id_INFO_LIFTOVER, 3); // account for 3 commas
    Context *lb_ctx = vcf_lo_seg_lo_snip (vb, liftback_snip, liftback_snip_len, dict_id_INFO_LIFTBACK, 0); // 0 as this container is not reconstructed by default
    
    if (dict_id.num == dict_id_INFO_LIFTOVER) lo_ctx->last_txt_len = value_len;
    else                                      lb_ctx->last_txt_len = value_len;
}

// Segging a Primary or Luft dual-coordinates VCF file, INFO/LIFTREJT field:
void vcf_lo_seg_INFO_REJTXXXX (VBlockVCFP vb, DictId dict_id, const char *value, int value_len)
{
    ASSVCF (!chain_is_loaded, "%s: --chain cannot be used with this file because it is already a dual-coordinates VCF file", txt_name);
    ASSVCF (txt_file->coords, "Found an %s subfield, but this is not a dual-coordinates VCF because it is missing \"" HK_DC_PRIMARY "\" or \"" HK_DC_LUFT "\" in the VCF header", dis_dict_id_ex (dict_id, true).s);

    Coords coord = (dict_id.num == dict_id_INFO_REJTOVER ? DC_PRIMARY : DC_LUFT);
    ASSVCF (vb->line_coords == coord, "Found an %s subfield, not expecting because this is a %s line", dis_dict_id_ex (dict_id, true).s, coords_names[vb->line_coords]);

    // value may be with or without a string eg "INFO/AC" ; "NO_MAPPING". We take the value up to the comma.
    const char *slash = memchr (value, '/', value_len);

    // get status from string
    SAFE_NUL (slash ? slash : &value[value_len]);

    unsigned i; for (i=LO_REJECTED; i < NUM_LO_STATUSES; i++)
        if (!strcmp (value, liftover_status_names[i])) { // found
            vcf_set_ostatus (i);
            break;
        }
    
    if (i == NUM_LO_STATUSES)
        vcf_set_ostatus (LO_REJECTED); // fallback

    SAFE_RESTORE;

    seg_by_did_i (vb, value, value_len, VCF_oSTATUS, value_len); // note: word_index >= NUM_LO_STATUSES if unrecognized string, that's ok
    vcf_lo_seg_REJTXXXX_do (vb, 0);
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
        sprintf (z_file->rejects_file_name[coord-1], "%s.%s_ONLY%s", z_file->name, coords_names[coord], file_plain_ext_by_dt (z_file->data_type));
        
        z_file->rejects_file[coord-1] = fopen (z_file->rejects_file_name[coord-1], "wb+");
        ASSERT (z_file->rejects_file[coord-1], "fopen() failed to open %s: %s", z_file->rejects_file_name[coord-1], strerror (errno));
    }

    ASSERT0 (z_file->rejects_file[coord-1], "liftover rejects file is not open");

    fwrite (vb->lo_rejects[coord-1].data, 1, vb->lo_rejects[coord-1].len, z_file->rejects_file[coord-1]);
    z_file->rejects_disk_size[coord-1] += vb->lo_rejects[coord-1].len;

    buf_free (&vb->lo_rejects[coord-1]);
}

