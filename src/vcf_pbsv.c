// ------------------------------------------------------------------
//   vcf_pbsv.c
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

#define _ID0  DICT_ID_MAKEF_8("I0D_VTYP")
#define _ID1  DICT_ID_MAKEF_7("I1D_SPC")
#define _ID_BND0 DICT_ID_MAKEF_7("I1chrom")
#define _ID_BND1 DICT_ID_MAKEF_5("I1pos")
#define _ID_BND2 DICT_ID_MAKEF_6("I1mate")

#define decl_pbsv_ID_ctxs                                                               \
    ContextP id0_ctx  __attribute__((unused)) = ctx_get_ctx (vb, _ID0);  \
    ContextP id1_ctx  __attribute__((unused)) = ctx_get_ctx (vb, _ID1);  \
    ContextP bnd0_ctx __attribute__((unused)) = ctx_get_ctx (vb, _ID_BND0); \
    ContextP bnd1_ctx __attribute__((unused)) = ctx_get_ctx (vb, _ID_BND1); \
    ContextP bnd2_ctx __attribute__((unused)) = ctx_get_ctx (vb, _ID_BND2); 

static SmallContainer id_bnd_con = { // goes into I1D
    .nitems_lo = 3,
    .repeats   = 1,
    .items = { { .dict_id.num = _ID_BND0, .separator[0] = ':' },     // copy of CHROM
               { .dict_id.num = _ID_BND1, .separator[0] = '-' },     // copy of POS
               { .dict_id.num = _ID_BND0,                     } } }; // copy of chrom:pos of mate as appears in ALT

static Did tw_dids[NUM_PBSV_TWs] = PBSV_TW_DIDS;

static int i0d_mux_channel[NUM_VTs] = { [VT_BND]=0, [VT_SYM_CNV]=1, [VT_DEL]=2, [VT_SYM_DEL]=2, [VT_INS]=3, [VT_SYM_INS]=3, [VT_SYM_INV]=4, [VT_SYM_DUP]=5 };
static int i1d_mux_channel[NUM_VTs] = { [VT_BND]=2, [VT_SYM_CNV]=1, [VT_DEL]=0, [VT_SYM_DEL]=0, [VT_INS]=0, [VT_SYM_INS]=0, [VT_SYM_INV]=0, [VT_SYM_DUP]=0 };

sSTRl(con_id_bnd_snip, 64);

void vcf_pbsv_zip_initialize (void)
{
    DO_ONCE {
        container_prepare_snip ((ContainerP)&id_bnd_con, 0, 0, qSTRa(con_id_bnd_snip));
    }
    
    // re-initialize for every file, as SV type might change
    segconf.MATEID_method = MATE_PBSV;

    vcf_sv_zip_initialize (tw_dids, NUM_PBSV_TWs);
}   

void vcf_pbsv_seg_initialize (VBlockVCFP vb)
{
    decl_pbsv_ID_ctxs;

    vcf_sv_seg_initialize (vb, tw_dids, NUM_PBSV_TWs);

    ctx_set_store_per_line (VB, VCF_ID, VCF_FILTER, INFO_MATEID, DID_EOL);
    
    ctx_consolidate_stats (VB, VCF_ID, id0_ctx->did_i, id1_ctx->did_i, DID_EOL);      
    ctx_consolidate_stats_(VB, CTX(VCF_ID), (ContainerP)&id_bnd_con);      

    ctx_set_dyn_int (VB, id1_ctx->did_i, bnd1_ctx->did_i, DID_EOL);

    seg_mux_init (vb, id0_ctx->did_i, VCF_SPECIAL_DEMUX_BY_VARTYPE, true,  pbsv_I0D);
    seg_mux_init (vb, id1_ctx->did_i, VCF_SPECIAL_DEMUX_BY_VARTYPE, false, pbsv_I1D);
}

StrTextLong vcf_pbsv_get_mate_id (STRp(id), uint32_t *len/*out*/)
{
    StrTextLong mate_id = {};
    *len = 0;

    if (id_len >= sizeof (mate_id.s) || id_len <= 9 || memcmp (id, "pbsv.BND.", 9))
        return mate_id; // failed

    str_split (id+9, id_len-9, 2, '-', id_item, true);
    if (!n_id_items) return mate_id; // failed   

    *len = snprintf (mate_id.s, sizeof(mate_id.s), "pbsv.BND.%.*s-%.*s", STRfi(id_item, 1), STRfi(id_item, 0));

    return mate_id;
}

// example: ID: pbsv.BND.hs37d5:35428322-hs37d5:24064194   MATEID: pbsv.BND.hs37d5:24064194-hs37d5:35428322
void vcf_seg_pbsv_MATEID (VBlockVCFP vb, ContextP ctx, STRp(mate_id))
{
    STRlast (id, VCF_ID);

    if (mate_id_len > 9 && !memcmp (mate_id, "pbsv.BND.", 9) &&
        id_len      > 9 && !memcmp (id,      "pbsv.BND.", 9)) {

        str_split (id+9,      id_len-9,      2, '-', id_item,      true);
        str_split (mate_id+9, mate_id_len-9, 2, '-', mate_id_item, true);

        if (n_id_items && n_mate_id_items && 
            str_issame_(STRi(id_item, 0), STRi(mate_id_item,1)) &&
            str_issame_(STRi(id_item, 1), STRi(mate_id_item,0))) {

            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_PBSV_MATEID }, 2, ctx, mate_id_len);
            return;
        }
    }

    // fallback
    seg_by_ctx (VB, STRa(mate_id), ctx, mate_id_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_PBSV_MATEID)
{    
    if (!reconstruct) return NO_NEW_VALUE;

    STRlast (id, VCF_ID);

    snip = vcf_pbsv_get_mate_id (STRa(id), &snip_len).s;
    RECONSTRUCT_snip;

    return NO_NEW_VALUE;
}

// returns id0_str1 or id0_str2 if matches, or NULL if not 
static rom match_vt (VariantType vt, pSTRp(id), STRp(id0_str1), STRp(id0_str2), uint32_t *vt_item_len)
{
    if (*id_len <= 5 || memcmp (*id, "pbsv.", 5)) 
        return NULL; // "pbsv." prefix doesn't match

    if (*id_len > 5 + id0_str1_len && !memcmp (*id + 5, id0_str1, id0_str1_len) && (*id)[5 + id0_str1_len]=='.') {
        *id     += 5 + id0_str1_len + 1; // after the '.' separator
        *id_len -= 5 + id0_str1_len + 1;
        *vt_item_len = id0_str1_len;
        return id0_str1;
    }

    if (id0_str2 && *id_len > 5 + id0_str2_len && !memcmp (*id + 5, id0_str2, id0_str2_len)  && (*id)[5 + id0_str2_len]=='.') {
        *id     += 5 + id0_str2_len + 1;
        *id_len -= 5 + id0_str2_len + 1;
        *vt_item_len = id0_str2_len;
        return id0_str2;
    }

    return NULL; // I0D string matches neither str1 nor str2        
} 

// seg "1:10002-11:176164" component of "pbsv.BND.1:10002-11:176164"
static bool vcf_seg_pbsv_ID_BND (VBlockVCFP vb, ContextP ctx, STRp(bnd))
{
    str_split (bnd, bnd_len, 2, '-', coord, true);
    if (!n_coords) return false;

    str_split (coords[0], coord_lens[0], 2, ':', mine, true);
    if (!n_mines) return false;

    str_split (coords[1], coord_lens[1], 2, ':', hers, true);
    if (!n_herss) return false;

    int64_t my_pos, her_pos;

    if (!str_issame_ (STRi(mine,0), STRa(vb->chrom_name)) ||                             // my chrom matches CHROM
        !str_get_int (STRi(mine,1), &my_pos) || my_pos != DATA_LINE(vb->line_i)->pos ||  // my pos matches POS
        !str_issame_ (STRi(hers,0), STRa(vb->mate_chrom_name)) ||                        // her chrom matches mate chrom from ALT
        !str_get_int (STRi(hers,1), &her_pos) || her_pos != vb->mate_pos)                // her pos matches mate pos from ALT
        return false;

    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_PBSV_ID_BND }, 2, ctx, bnd_len);

    return true;
} 

SPECIAL_RECONSTRUCTOR (vcf_piz_special_PBSV_ID_BND)
{    
    if (reconstruct) {
        RECONSTRUCT_str (vb->chrom_name);
        RECONSTRUCT1 (':');
        RECONSTRUCT_INT (CTX(VCF_POS)->last_value.i);
        RECONSTRUCT1 ('-');
        RECONSTRUCT_str (VB_VCF->mate_chrom_name);
        RECONSTRUCT1 (':');
        RECONSTRUCT_INT (VB_VCF->mate_pos);
    }

    return NO_NEW_VALUE;
}

void vcf_seg_pbsv_ID (VBlockVCFP vb, STRp(id))
{
    STR(mate_id);
    mate_id = vcf_pbsv_get_mate_id (STRa(id), &mate_id_len).s;

    if (mate_id_len)
        vcf_seg_BND_mate (vb, STRa(id), STRa(mate_id), crc32 (0, STRa(mate_id))); 
    
    if (vcf_has_mate)
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPY_MATE }, 2, VCF_ID, id_len + 1); // +1 for \t

    else {
        decl_pbsv_ID_ctxs;
    
        // seg I0D - the vartype component of e.g. pbsv.DEL.29378
        STR(vt_item);
        VariantType vt = vb->var_types[0];

        // note: special case: DUP is different than SVTYPE string, and has two options
        if ((ALT0(SYM_DUP) && (vt_item = match_vt (vt, pSTRa(id), "SPLIT.DUP", 9, "INS.DUP", 7, &vt_item_len))) ||
            (!ALT0(SYM_DUP) && svtype_by_vt[vt] && (vt_item = match_vt (vt, pSTRa(id), svtype_by_vt[vt], strlen (svtype_by_vt[vt]), 0, 0, &vt_item_len)))) {

            ContextP channel_ctx = 
                seg_mux_get_channel_ctx (VB, id0_ctx->did_i, (MultiplexerP)&vb->mux_pbsv_I0D, i0d_mux_channel[vt]);

            seg_by_ctx (VB, STRa(vt_item), channel_ctx, vt_item_len);
            seg_by_ctx (VB, STRa(vb->mux_pbsv_I0D.snip), id0_ctx, 0); // I0D de-multiplexer
        }

        else 
            goto fallback;

        // seg I1D - demux 3 types: 
        // channel 0 - INS, DEL, DUP, INV - close numbers
        // channel 1 - CNV - its own number range
        // channel 2 - BND - chrom/pos of this and mate
        ContextP channel_ctx = 
            seg_mux_get_channel_ctx (VB, id1_ctx->did_i, (MultiplexerP)&vb->mux_pbsv_I1D, i1d_mux_channel[vt]);

        if (ALT0(BND)) {
            if (!vcf_seg_pbsv_ID_BND (vb, channel_ctx, STRa(id)))
                goto ID1_fallback;
        }

        else {
            int64_t i1d;
            if (str_get_int (STRa(id), &i1d))
                seg_self_delta (VB, channel_ctx, i1d, 0, 0, id_len);
            
            else ID1_fallback:
                seg_by_ctx (VB, STRa(id), channel_ctx, id_len); // fallback
        }

        seg_by_ctx (VB, STRa(vb->mux_pbsv_I1D.snip), id1_ctx, 0); // I1D de-multiplexer

        seg_by_did (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_DEFER }, 2, VCF_ID, 7); // "pbsv" + 2 '.' separators + \t
    }

    return;
    
fallback:
    seg_by_did (VB, STRa(id), VCF_ID, id_len + 1);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_DEMUX_BY_VARTYPE)
{
    return reconstruct_demultiplex (vb, ctx, STRa(snip), (ctx->dict_id.num == _ID0 ? i0d_mux_channel : i1d_mux_channel)[VB_VCF->var_types[0]], new_value, reconstruct);
}

// called from vcf_piz_refalt_parse
void vcf_piz_insert_pbsv_ID (VBlockVCFP vb)
{
    decl_ctx (VCF_ID);
    
    if (IS_RECON_INSERTION(ctx)) {
        ContextP id0_ctx = ECTX (_ID0);
        ContextP id1_ctx = ECTX (_ID1);

        rom id_recon = BAFTtxt;

        RECONSTRUCT ("pbsv.", 5);
        reconstruct_from_ctx (VB, id0_ctx->did_i, '.', true);
        reconstruct_from_ctx (VB, id1_ctx->did_i, 0,   true);
        
        uint32_t id_len = BAFTtxt - id_recon;
        Ltxt -= id_len; 

        char id[id_len];
        memcpy (id, id_recon, id_len);
        
        vcf_piz_insert_field (vb, ctx, STRa(id), segconf.wid_ID.width);
    }
}
