// ------------------------------------------------------------------
//   sam_sag_piz.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "huffman.h"

// reconstruct PRIM or DEPN line from SA Group data: (rname, pos, strand, CIGAR, mapQ, NM ;)+
static void sam_sa_reconstruct_SA_from_SA_Group (VBlockSAMP vb, 
                                                 int16_t prim_aln_index_in_SA_Z) // DEPN: 1-based location of primary alignment withing DEPN. 0 up to 15.0.65, which is treated as "1"
{
    SAAln *vb_alns = B(SAAln, z_file->sag_alns, vb->sag->first_aln_i);
    if (!IS_DEPN(vb)) prim_aln_index_in_SA_Z = -1;
    else if (prim_aln_index_in_SA_Z >= 1) prim_aln_index_in_SA_Z--; // convert to 0-based

    for (uint32_t grp_aln_i=1, SA_aln_i=0; SA_aln_i < vb->sag->num_alns-1; SA_aln_i++) {   
        SAAln *a;
        bool cigar_abbreviated, nm_implied_by_cigar;

        // case: reconstruct the primary alignment as alignment number "prim_aln_index_in_SA_Z" in SA:Z 
        if (prim_aln_index_in_SA_Z == SA_aln_i) { 
            a = &vb_alns[0]; 

            // note: if SA_CIGARs are abbreviated (minimap2 etc), the SAG contains abbreviated cigars for all
            // alignments except for PRIM, for which we need to abbreviate the cigar during reconstruction         
            cigar_abbreviated   = segconf.SA_CIGAR_abbreviated;
            nm_implied_by_cigar = segconf.SA_NM_by_CIGAR_X;
        }
        
        // case: reconstruct a non-primary alignment in SA:Z according to order of group
        else {
            if (grp_aln_i == 0) grp_aln_i++; // don't reconstruct the primary alignment as its not its place (including: if PRIM, we don't reconstruct it at all with SA:Z)

            a = &vb_alns[grp_aln_i++];
            if (a == vb->sa_aln) a = &vb_alns[grp_aln_i++]; // skip my own alignment - already reconstructed in main SAM fields

            cigar_abbreviated = nm_implied_by_cigar = false;
        }

        // rname
        STR(rname);
        ctx_get_snip_by_word_index (CTX(OPTION_SA_RNAME), a->rname, rname);
        RECONSTRUCT_SEP (rname, rname_len, ',');

        // pos
        RECONSTRUCT_INT (a->pos);
        RECONSTRUCT1 (',');

        // strand
        RECONSTRUCT ((a->revcomp ? "-," : "+,"), 2);

        // cigar
        uint32_t X_bases = sam_reconstruct_SA_cigar_from_SA_Group (vb, a, cigar_abbreviated, nm_implied_by_cigar);

        // mapq
        RECONSTRUCT_INT (a->mapq);
        RECONSTRUCT1 (',');

        // nm
        RECONSTRUCT_INT (nm_implied_by_cigar ? X_bases : a->nm);
        RECONSTRUCT1 (';');
    }
}

// PIZ: reconstruction of a PRIM or DEPN VB: Called in reconstruction of lines with SPECIAL_pull_from_sag (=all lines if PRIM)
void sam_piz_set_sag (VBlockSAMP vb)
{
    if (SAM_PIZ_HAS_SAG) return; // vb->sag and sa_aln are already set for this line 

    // set vb->sag
    if (IS_DEPN(vb)) {
        // get the PRIM VB and the sag within it which are the SA Group of this DEPN line
        int64_t sa_grp_i = reconstruct_from_local_int (VB, CTX(SAM_SAG), 0, RECON_OFF);

        ASSPIZ (sa_grp_i < z_file->sag_grps.len, "sa_grp_i=%"PRId64" ∉ [0,%"PRId64"]", sa_grp_i, z_file->sag_grps.len-1);
        vb->sag = B(const Sag, z_file->sag_grps, sa_grp_i);
    }

    else  // PRIM
        vb->sag = vb->sag ? (vb->sag + 1) : sam_piz_get_prim_vb_first_sa_grp (vb);

    // set vb->sa_aln (only for SA:Z-type SA Groups)
    if (IS_SAG_SA) {
        // The primary alignment is the first alignment, then the alignments from the SA Groups
        int64_t sa_aln_i = vb->sag->first_aln_i + (IS_DEPN(vb) ? reconstruct_from_local_int (VB, CTX(SAM_SAALN), 0, RECON_OFF) : 0);

        ASSPIZ (sa_aln_i < z_file->sag_alns.len, "sa_aln_i=%"PRId64" ∉ [0,%"PRIu64"]", sa_aln_i, z_file->sag_alns.len-1);
        vb->sa_aln = B(const SAAln, z_file->sag_alns, sa_aln_i);
    }

    else if (IS_SAG_CC)
        vb->cc_aln = B(CCAln, z_file->sag_alns, ZGRP_I(vb->sag));

    else if (IS_SAG_SOLO)
        vb->solo_aln = B(SoloAln, z_file->sag_alns, ZGRP_I(vb->sag));

    // indicate that we have already updated sag and sa_aln for this line_i
    vb->sag_line_i = vb->line_i;
}

static void sam_reconstruct_solo_from_sag (VBlockSAMP vb, Did did_i, SoloTags solo)
{
    bytes comp = sam_solo_sag_data (vb, solo);
    uint32_t uncomp_len = vb->solo_aln->field_uncomp_len[solo];
            
    RECONSTRUCT_huffman_or_copy (VB, did_i, uncomp_len, comp);
}

// PIZ compute thread: called when reconstructing a PRIM or DEPN line - reconstruct pulling info from
// SA Groups loaded to z_file->sa_*
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_pull_from_sag)
{
    VBlockSAMP vb = (VBlockSAMP)vb_; 
    
    sam_piz_set_sag (vb);

    const Sag *g = vb->sag; ASSERTNOTNULL (g); // pointer into z_file->sag_grps
    const SAAln   *sa_aln   = NULL; 
    const CCAln   *cc_aln   = NULL;
    const SoloAln *solo_aln = NULL;

    ASSERTNOTZERO (segconf.sag_type);
    
    if    (IS_SAG_SOLO) ASSERTNOTNULL ((solo_aln = vb->solo_aln));  // pointer into z_file->sag_alns
    else if (IS_SAG_SA) ASSERTNOTNULL ((sa_aln = vb->sa_aln));  
    else if (IS_SAG_CC) ASSERTNOTNULL ((cc_aln = vb->cc_aln)); 

    switch (ctx->did_i) {
        case SAM_QNAME: 
            if (reconstruct) RECONSTRUCT_huffman_or_copy (VB, SAM_QNAME, g->qname_len, GRP_QNAME(g));
            return NO_NEW_VALUE;

        case SAM_RNAME: {
            STR(rname);
            ctx_get_snip_by_word_index (CTX(OPTION_SA_RNAME), sa_aln->rname, rname);
            if (reconstruct) RECONSTRUCT_str (rname);
            new_value->i = ctx_get_word_index_by_snip (VB, CTX(SAM_RNAME), STRa(rname)); // SAM_RNAME has STORE_INDEX. Note: SAM_RNAME has different word indices than OPTION_SA_RNAME
            return HAS_NEW_VALUE;
        }

        case SAM_POS:
            new_value->i = sa_aln->pos;
            if (reconstruct) RECONSTRUCT_INT (new_value->i);
            return HAS_NEW_VALUE;

        case SAM_MAPQ: 
            new_value->i = sa_aln->mapq;
            if (reconstruct) RECONSTRUCT_INT (new_value->i);
            return HAS_NEW_VALUE;

        case SAM_FLAG: case SAM_FLAG0: case SAM_FLAG1: { // only happens for SA-based, not NH-based, groups
            SamFlags sam_flags = { .value = atoi (snip) };
            if (IS_SAG_SA) sam_flags.rev_comp = sa_aln->revcomp; // revcomp is determined by alignment of SA:Z-type groups, and is in the snip for NH-tag groups
            sam_flags.multi_segs = g->multi_segs;
            sam_flags.is_first   = g->is_first;
            sam_flags.is_last    = g->is_last;

            new_value->i = sam_flags.value;
            if (reconstruct) RECONSTRUCT_INT (new_value->i);
            return HAS_NEW_VALUE;
        }

        case SAM_CIGAR:
            sam_reconstruct_main_cigar_from_sag (vb, segconf.SA_HtoS==yes && IS_DEPN(vb), reconstruct);
            return NO_NEW_VALUE;

        case OPTION_NM_i:
            new_value->i = sa_aln->nm;
            if (reconstruct) RECONSTRUCT_INT (new_value->i);
            return HAS_NEW_VALUE;

        case OPTION_SA_Z:
            if (reconstruct) sam_sa_reconstruct_SA_from_SA_Group (vb, atoi (snip));
            return NO_NEW_VALUE;

        case OPTION_NH_i:
            new_value->i = vb->sag->num_alns;
            if (reconstruct) RECONSTRUCT_INT (new_value->i);
            return HAS_NEW_VALUE;

        case OPTION_AS_i:
            new_value->i = vb->sag->as - atoi (snip);
            if (reconstruct) RECONSTRUCT_INT (new_value->i);
            return HAS_NEW_VALUE;

        case OPTION_CP_i:
            if (reconstruct) 
                RECONSTRUCT_INT (cc_aln->pos);
            
            new_value->i = cc_aln->pos;
            return HAS_NEW_VALUE;

        case SAM_SQBITMAP:
            ABORT_PIZ0 ("Expecting SAM_SQBITMAP to be handled by sam_piz_special_SEQ");
        
        case SAM_QUAL:
            ABORT_PIZ0 ("Expecting SAM_QUAL to be handled by sam_piz_special_QUAL");

        // solo tags
        case OPTION_BX_Z: if (reconstruct) sam_reconstruct_solo_from_sag (vb, ctx->did_i, SOLO_BX); return NO_NEW_VALUE;
        case OPTION_RX_Z: if (reconstruct) sam_reconstruct_solo_from_sag (vb, ctx->did_i, SOLO_RX); return NO_NEW_VALUE;
        case OPTION_CB_Z: if (reconstruct) sam_reconstruct_solo_from_sag (vb, ctx->did_i, SOLO_CB); return NO_NEW_VALUE;
        case OPTION_CR_Z: if (reconstruct) sam_reconstruct_solo_from_sag (vb, ctx->did_i, SOLO_CR); return NO_NEW_VALUE;
        case OPTION_BC_Z: if (reconstruct) sam_reconstruct_solo_from_sag (vb, ctx->did_i, SOLO_BC); return NO_NEW_VALUE;
        case OPTION_QX_Z: if (reconstruct) sam_reconstruct_solo_from_sag (vb, ctx->did_i, SOLO_QX); return NO_NEW_VALUE;
        case OPTION_CY_Z: if (reconstruct) sam_reconstruct_solo_from_sag (vb, ctx->did_i, SOLO_CY); return NO_NEW_VALUE;
        case OPTION_QT_Z: if (reconstruct) sam_reconstruct_solo_from_sag (vb, ctx->did_i, SOLO_QT); return NO_NEW_VALUE;

        default:
            ABORT_PIZ ("Unexpected ctx=%s(%u)", ctx->tag_name, ctx->did_i);
    }

    return NO_NEW_VALUE; // just to suppress compiler warning
}

// PIZ of a QNAME in a PRIM component (an all-the-same context) - copy QNAME from the in-memory SA Group
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_PRIM_QNAME)
{
    VBlockSAMP vb = (VBlockSAMP)vb_; 

    sam_piz_set_sag (vb);

    uint32_t last_txt_index = Ltxt;

    RECONSTRUCT_huffman_or_copy (VB, SAM_QNAME, vb->sag->qname_len, GRP_QNAME(vb->sag));

    CTX(SAM_QNAME)->last_txt = (TxtWord){ .index = last_txt_index,
                                          .len   = Ltxt - last_txt_index }; // 15.0.16 - needed by sam_ultima_bi_prediction

    // if seq_len is carried by a QNAME item, set the last value here - it was already reconstructed during sag loading.
    if (segconf.seq_len_dict_id.num)
        ctx_set_last_value (VB, ECTX(segconf.seq_len_dict_id), (int64_t)vb->sag->seq_len);

    return NO_NEW_VALUE;
}
