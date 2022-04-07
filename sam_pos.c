// ------------------------------------------------------------------
//   sam_pos.c
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

// a module for handling POS and PNEXT

#include "genozip.h"
#include "sam_private.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "random_access.h"
#include "segconf.h"

static void sam_seg_POS_segconf (VBlockSAMP vb, WordIndex prev_line_chrom, PosType pos, PosType prev_line_pos)
{
    // evidence of not being sorted: our RNAME is the same as the previous line, but POS has decreased
    if (segconf.sam_is_sorted && prev_line_chrom == vb->chrom_node_index && prev_line_pos > pos)
        segconf.sam_is_sorted = false;

    // evidence of not being sorted: our RNAME is different than previous line, but we encountered it before
    if (segconf.sam_is_sorted && (prev_line_chrom != NODE_INDEX_NONE) && (prev_line_chrom != vb->chrom_node_index) && 
        *B32 (CTX(SAM_RNAME)->counts, vb->chrom_node_index) > 1) // 1 if it has been segged on this line for the first time
        segconf.sam_is_sorted = false;
    
    // evidence of not being entirely unmapped: we have POS in at least one line
    if (pos)
        segconf.sam_is_unmapped = false; 
}

PosType sam_seg_POS (VBlockSAMP vb, ZipDataLineSAM *dl, WordIndex prev_line_chrom, unsigned add_bytes)
{
    ContextP ctx = CTX(SAM_POS);
    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1
    PosType pos = dl->POS;
    PosType prev_line_pos = vb->line_i ? (dl-1)->POS : 0;

    // case: DEPN or PRIM line.
    // Note: in DEPN, pos already verified in sam_sa_seg_depn_find_sagroup to be as in SA alignment
    if (sam_seg_has_SA_Group(vb)) {
        sam_seg_against_sa_group (vb, ctx, add_bytes);
        ctx_set_last_value (VB, ctx, pos);

        // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
        if (sam_is_prim_vb) {
            seg_pos_field (VB, OPTION_SA_POS, OPTION_SA_POS, 0, 0, 0, 0, dl->POS, 0);

            // count POS field contribution to OPTION_SA_POS, so sam_stats_reallocate can allocate the z_data between POS and SA:Z
            CTX(OPTION_SA_POS)->counts.count += add_bytes; 
        }
    }

    // in a collated file, we expect that in about half of the lines, the POS will be equal to the previous PNEXT.
    else if (segconf.sam_is_collated && pos != ctx->last_value.i && pos == CTX(SAM_PNEXT)->last_value.i) {

        seg_pos_field (VB, SAM_POS, SAM_PNEXT, 0, 0, 0, 0, pos, add_bytes);

        ctx->last_delta = pos - prev_line_pos; // always store a self delta, even if we delta'd against PNEXT. This mirrors flags.store_delta we set for PIZ
    }
    
    // seg against buddy
    else if (segconf.sam_is_sorted && vb->buddy_line_i != -1 && buddy_dl->PNEXT == pos) {
        seg_by_did_i (VB, STRa(POS_buddy_snip), SAM_POS, add_bytes); // copy POS from earlier-line buddy PNEXT
        ctx_set_last_value (VB, ctx, pos);
    }

    else  
        pos = seg_pos_field (VB, SAM_POS, SAM_POS, 0, 0, 0, 0, pos, add_bytes);
    
    random_access_update_pos (VB, 0, SAM_POS);

    dl->POS = pos;

    if (segconf.running) sam_seg_POS_segconf (vb, prev_line_chrom, pos, prev_line_pos);

    return pos;
}

void sam_seg_PNEXT (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(pnext_str)/* option 1 */, PosType pnext/* option 2 */, PosType prev_line_pos, unsigned add_bytes)
{
    PosType this_pnext=0;
    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // note: an invalid pointer if buddy_line_i is -1

    if (!pnext) str_get_int (STRa(pnext_str), &pnext);

    if (segconf.sam_is_collated && pnext) { // note: if sam_is_sorted, this will be handled by buddy

        // case: 2nd mate - PNEXT = previous line POS 
        if (pnext == prev_line_pos) {
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_PNEXT_IS_PREV_POS}, 2, CTX(SAM_PNEXT), add_bytes);
            ctx_set_last_value (VB, CTX(SAM_PNEXT), pnext);
        }

        // case: supplamentary alignment - PNEXT = previous line PNEXT
        else if (this_pnext == CTX(SAM_PNEXT)->last_value.i)
            seg_pos_field (VB, SAM_PNEXT, SAM_PNEXT, 0, 0, 0, 0, pnext, add_bytes);

        else
            seg_pos_field (VB, SAM_PNEXT, SAM_POS, 0, 0, 0, 0, pnext, add_bytes);
    }

    else if (segconf.sam_is_sorted && pnext && vb->buddy_line_i != -1 && buddy_dl->POS == pnext) {
        seg_by_did_i (VB, STRa(PNEXT_buddy_snip), SAM_PNEXT, add_bytes); // copy PNEXT from earlier-line buddy POS
        ctx_set_last_value (VB, CTX(SAM_PNEXT), pnext);
    }

    else
        pnext = seg_pos_field (VB, SAM_PNEXT, SAM_POS, 0, 0, 0, 0, pnext, add_bytes);

    dl->PNEXT = pnext;
}
