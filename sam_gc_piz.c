// ------------------------------------------------------------------
//   sam_gc_piz.c
//   Copyright (C) 2022-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "sam_private.h"
#include "sections.h"
#include "codec.h"
#include "piz.h"
#include "reconstruct.h"

// reconstruct PRIM or DEPN line from SA Group data: (rname, pos, strand, CIGAR, mapQ, NM ;)+
static void sam_sa_reconstruct_SA_from_SA_Group (VBlockSAMP vb)
{
    SAAlnType *vb_alns = B(SAAlnType, z_file->sa_alns, vb->sa_grp->first_aln_i);

    for (uint32_t aln_i=0; aln_i < vb->sa_grp->num_alns; aln_i++) {
        SAAlnType *a = &vb_alns[aln_i];
        if (a == vb->sa_aln) continue; // skip my own alignment - already reconstruct in main SAM fields

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

// PIZ: reconstruction of a PRIM or DEPN VB: Called in reconstruction of lines with SPECIAL_SAGROUP (=all lines if PRIM)
void sam_piz_set_sa_grp (VBlockSAMP vb)
{
    if (SAM_PIZ_HAS_SA_GROUP) return; // vb->sa_grp and sa_aln are already set for this line 

    // get sa_grp_i
    int64_t sa_grp_i, sa_aln_i;
    if (sam_is_depn_vb) {
        // get the PRIM VB and the sa_grp within it which are the SA Group of this DEPN line
        sa_grp_i = reconstruct_from_local_int (VB, CTX(SAM_SAGROUP), 0, false);

        ASSPIZ (sa_grp_i < z_file->sa_groups.len, "sa_grp_i=%"PRId64" is out of range [0,%"PRId64"]", sa_grp_i, z_file->sa_groups.len-1);
        vb->sa_grp = B(const SAGroupType, z_file->sa_groups, sa_grp_i);

        // get the alignment within the SA Group that is the alignment of this line
        sa_aln_i = vb->sa_grp->first_aln_i + reconstruct_from_local_int (VB, CTX(SAM_SAALN), 0, false); 
    }

    else { // PRIM
        if (!vb->sa_grp_line_i) 
            vb->sa_grp = sam_piz_get_prim_vb_first_sa_grp (vb);
        else
            vb->sa_grp++; // groups within a VB always have a consecutive grp_i

        // The primary alignment is the first alignment
        sa_aln_i = vb->sa_grp->first_aln_i;

        // set buddy if needed and not already set. Note: vb->sa_grp->prim_set_buddy is for PRIM lines. DEPN lines get
        // their buddy from SAM_QNAME (see sam_seg_QNAME)
        if (vb->sa_grp->prim_set_buddy && vb->buddy_line_i == NO_BUDDY)
            reconstruct_set_buddy (VB);
    }

    ASSPIZ (sa_aln_i < z_file->sa_alns.len, "sa_aln_i=%"PRId64" is out of range [0,%"PRIu64"]", sa_aln_i, z_file->sa_alns.len-1);
    vb->sa_aln = B(const SAAlnType, z_file->sa_alns, sa_aln_i);

    // indicate that we have already updated sa_grp and sa_aln for this line_i
    vb->sa_grp_line_i = vb->line_i + 1; // +1 as sa_grp_line_i==0 means "not set"
}

// PIZ compute thread: called when reconstructing a PRIM or DEPN line - reconstruct pulling info from
// SA Groups loaded to z_file->sa_*
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_pull_from_SAGROUP)
{
    VBlockSAMP vb = (VBlockSAMP)vb_; 
    
    sam_piz_set_sa_grp (vb);

    const SAGroupType *g = vb->sa_grp; ASSERTNOTNULL (g); // pointer into z_file->sa_groups
    const SAAlnType   *a = vb->sa_aln; ASSERTNOTNULL (a); // pointer into z_file->sa_alns
// printf ("g=%u g->first_aln=%u a=%u\n", ZGRP_I(g), g->first_aln_i, ZALN_I(a));
    switch (ctx->did_i) {
    
        case SAM_QNAME: {
            rom qname = GRP_QNAME(g);
            if (reconstruct) RECONSTRUCT (qname, g->qname_len);
            return false; // no new value
        }

        case SAM_RNAME: {
            STR(rname);
            ctx_get_snip_by_word_index (CTX(OPTION_SA_RNAME), a->rname, rname);
            if (reconstruct) RECONSTRUCT (rname, rname_len);
            new_value->i = ctx_get_word_index_by_snip (VB, CTX(SAM_RNAME), STRa(rname)); // SAM_RNAME has STORE_INDEX. Note: SAM_RNAME has different word indices than OPTION_SA_RNAME
            return true;
        }

        case SAM_POS:
            new_value->i = a->pos;
            if (reconstruct) RECONSTRUCT_INT (new_value->i);
            return true;

        case SAM_MAPQ: 
            new_value->i = a->mapq;
            if (reconstruct) RECONSTRUCT_INT (new_value->i);
            return true;

        case SAM_FLAG: {
            SamFlags sam_flags = { .value = atoi (snip) };
            sam_flags.bits.rev_comp        = a->revcomp;
            sam_flags.bits.multi_segments  = g->multi_segments;
            sam_flags.bits.is_first        = g->is_first;
            sam_flags.bits.is_last         = g->is_last;

            new_value->i = sam_flags.value;
            if (reconstruct) RECONSTRUCT_INT (new_value->i);
            return true;
        }

        case SAM_CIGAR:
            sam_reconstruct_main_cigar_from_SA_Group (vb, *snip == '0'/*replace S by H*/, reconstruct);
            return false;

        case OPTION_NM_i:
            new_value->i = a->nm;
            if (reconstruct) RECONSTRUCT_INT (new_value->i);
            return true;

        case OPTION_SA_Z:
            if (reconstruct) sam_sa_reconstruct_SA_from_SA_Group (vb);
            return false;

        case SAM_SQBITMAP:
            ASSPIZ0 (false, "Expecting SAM_SQBITMAP to be handled by sam_piz_special_DEPN_SEQ");
        
        case SAM_QUAL:
            ASSPIZ0 (false, "Expecting SAM_QUAL to be handled by sam_piz_special_QUAL");

        default:
            ASSPIZ ("Unexpected ctx=%s(%u)", ctx->tag_name, ctx->did_i);
    }

    return false; // just to suppress compiler warning
}

// PIZ of a QNAME in a PRIM component (an all-the-same context) - copy QNAME from the in-memory SA Group
SPECIAL_RECONSTRUCTOR (sam_piz_special_PRIM_QNAME)
{
    sam_piz_set_sa_grp (VB_SAM);

    RECONSTRUCT (GRP_QNAME(VB_SAM->sa_grp), VB_SAM->sa_grp->qname_len);
    
    return false; // no new value
}
