// ------------------------------------------------------------------
//   sam_sag_piz.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "sections.h"
#include "codec.h"
#include "piz.h"
#include "reconstruct.h"

// reconstruct PRIM or DEPN line from SA Group data: (rname, pos, strand, CIGAR, mapQ, NM ;)+
static void sam_sa_reconstruct_SA_from_SA_Group (VBlockSAMP vb)
{
    SAAln *vb_alns = B(SAAln, z_file->sag_alns, vb->sag->first_aln_i);

    for (uint32_t aln_i=0; aln_i < vb->sag->num_alns; aln_i++) {
        SAAln *a = &vb_alns[aln_i];
        if (a == vb->sa_aln) continue; // skip my own alignment - already reconstructed in main SAM fields

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
        sam_reconstruct_SA_cigar_from_SA_Group (vb, a);

        // mapq
        RECONSTRUCT_INT (a->mapq);
        RECONSTRUCT1 (',');

        // nm
        RECONSTRUCT_INT (a->nm);
        RECONSTRUCT1 (';');
    }
}

// PIZ: reconstruction of a PRIM or DEPN VB: Called in reconstruction of lines with SPECIAL_SAG (=all lines if PRIM)
void sam_piz_set_sag (VBlockSAMP vb)
{
    if (SAM_PIZ_HAS_SAG) return; // vb->sag and sa_aln are already set for this line 

    // set vb->sag
    if (IS_DEPN(vb)) {
        // get the PRIM VB and the sag within it which are the SA Group of this DEPN line
        int64_t sa_grp_i = reconstruct_from_local_int (VB, CTX(SAM_SAG), 0, RECON_OFF);

        ASSPIZ (sa_grp_i < z_file->sag_grps.len, "sa_grp_i=%"PRId64" is out of range [0,%"PRId64"]", sa_grp_i, z_file->sag_grps.len-1);
        vb->sag = B(const Sag, z_file->sag_grps, sa_grp_i);
    }

    else { // PRIM
        if (!vb->sag_line_i) 
            vb->sag = sam_piz_get_prim_vb_first_sa_grp (vb);
        else
            vb->sag++; // groups within a VB always have a consecutive grp_i
    }

    // set vb->sa_aln (only for SA:Z-type SA Groups)
    if (IS_SAG_SA) {
        // The primary alignment is the first alignment, then the alignments from the SA Groups
        int64_t sa_aln_i = vb->sag->first_aln_i + (IS_DEPN(vb) ? reconstruct_from_local_int (VB, CTX(SAM_SAALN), 0, RECON_OFF) : 0);

        ASSPIZ (sa_aln_i < z_file->sag_alns.len, "sa_aln_i=%"PRId64" is out of range [0,%"PRIu64"]", sa_aln_i, z_file->sag_alns.len-1);
        vb->sa_aln = B(const SAAln, z_file->sag_alns, sa_aln_i);
    }

    else if (IS_SAG_CC)
        vb->cc_aln = B(CCAln, z_file->sag_alns, ZGRP_I(vb->sag));

    else if (IS_SAG_SOLO)
        vb->solo_aln = B(SoloAln, z_file->sag_alns, ZGRP_I(vb->sag));

    // indicate that we have already updated sag and sa_aln for this line_i
    vb->sag_line_i = vb->line_i + 1; // +1 as sag_line_i==0 means "not set"
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

// printf ("g=%u g->first_aln=%u sa_aln=%u\n", ZGRP_I(g), g->first_aln_i, ZALN_I(sa_aln));
    switch (ctx->did_i) {
    
        case SAM_QNAME: {
            rom qname = GRP_QNAME(g);
            if (reconstruct) RECONSTRUCT (qname, g->qname_len);
            return NO_NEW_VALUE;
        }

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
            if (reconstruct) sam_sa_reconstruct_SA_from_SA_Group (vb);
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

        default:
            // solo tags
            for (SoloTags tag_i=0; tag_i < NUM_SOLO_TAGS; tag_i++)
                if (ctx->did_i == solo_props[tag_i].did_i) {
                    if (reconstruct) RECONSTRUCT (Bc(z_file->solo_data, vb->solo_aln->word[tag_i].index), vb->solo_aln->word[tag_i].len); 
                    return NO_NEW_VALUE;
                }

            // not a solo tag
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

    RECONSTRUCT (GRP_QNAME(vb->sag), vb->sag->qname_len);

    CTX(SAM_QNAME)->last_txt = (TxtWord){ .index = last_txt_index,
                                          .len   = Ltxt - last_txt_index }; // 15.0.16 - need by sam_ultima_bi_prediction

    // if seq_len is carried by a QNAME item, set the last value here - it was already reconstructed during sag loading.
    if (segconf.seq_len_dict_id.num)
        ctx_set_last_value (VB, ECTX(segconf.seq_len_dict_id), (int64_t)vb->sag->seq_len);

    return NO_NEW_VALUE;
}
