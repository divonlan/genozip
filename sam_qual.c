// ------------------------------------------------------------------
//   sam_qual.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "reference.h"
#include "segconf.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "optimize.h"
#include "codec.h"
#include "htscodecs/rANS_static4x16.h"

rom bam_qual_display (bytes qual, uint32_t l_seq) // caller should free memory
{
    char *str = MALLOC (l_seq + 2);

    for (uint32_t i=0; i < l_seq; i++) 
        str[i] = qual[i] + 33;

    str[l_seq] = 0;
    return str;
}

// get a score of the QUAL string - similar to that calculated by biobambam for ms:i:
// see here: https://github.com/gt1/libmaus/tree/master/src/libmaus/bambam/BamAlignmentDecoderBase.cpp getScore 
// The sum of phred values of the QUAL string, but only for phred values >= 15
static int64_t sam_get_QUAL_score (VBlockSAMP vb, STRp(qual))
{
    if (vb->qual_missing) return 0;

    int64_t score=0;
    for (uint32_t i=0; i < qual_len; i++) {
        char phred = qual[i] - 33;
        if (phred >= 15) score += phred;
    }
    
    return score;
}

// main thread (not thread-safe): called from sam_show_sa_one_grp for getting first few characters of alignment cigar
rom sam_display_qual_from_SA_Group (const SAGroupType *g)
{
    if (g->no_qual) return "*";

    static char qual[SA_QUAL_DISPLAY_LEN+1];
    memset (qual, 0, sizeof(qual));

    uint32_t uncomp_len = MIN_(SA_CIGAR_DISPLAY_LEN, g->seq_len); // possibly shorter than original cigar

    if (g->qual_comp_len) { // qual of this group is compressed
        uint32_t uncomp_len = MIN_(g->seq_len, SA_QUAL_DISPLAY_LEN);
        void *success = rans_uncompress_to_4x16 (evb, B8(z_file->sa_qual, g->qual), g->qual_comp_len,
                                                (uint8_t *)qual, &uncomp_len); 
        if (success && uncomp_len) qual[uncomp_len] = '\0';
    }

    else // not compressed
        memcpy (qual, B8(z_file->sa_qual, g->qual), uncomp_len);

    return qual;
}

//---------
// SEG
//---------

// callback function for compress to get data of one line
COMPRESSOR_CALLBACK (sam_zip_qual) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_data_len  = MIN_(maximum_size, dl->QUAL.len);

    if (!line_data) return; // only lengths were requested

    *line_data = Btxt(dl->QUAL.index);

    // if QUAL is just "*" (i.e. unavailable) replace it by " " because '*' is a legal PHRED quality value that will confuse PIZ
    if (dl->QUAL.len == 1 && (*line_data)[0] == '*') 
        *line_data = " "; // pointer to static string

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    else if (flag.optimize_QUAL) 
        optimize_phred_quality_string (STRa(*line_data));

    if (is_rev) *is_rev = dl->FLAG.bits.rev_comp;
}

void sam_seg_QUAL_initialize (VBlockSAMP vb)
{
    if (segconf.sam_ms_type == ms_BIOBAMBAM) {     // handle QUAL_scores for ms:i - added v13
        CTX(SAM_QUAL)->flags.store_per_line = true;
        CTX(SAM_QUAL)->no_stons = true; // since we're storing QUAL data in local
    }

    // in case of BIOBAMBAM or DEPN we go through the SPECIAL for extra functionality
    // otherwise, the QUAL string is just reconstructed directly from local
    if (segconf.sam_ms_type == ms_BIOBAMBAM || sam_is_depn_vb)
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_QUAL, '0'+(segconf.sam_ms_type == ms_BIOBAMBAM) }, 3, SAM_QUAL, 0);

    // if PRIM - TOPLEVEL reconstructs all-the-same SAM_QUALSA instead of SAM_QUAL. SAM_QUAL is consumed when loading SA Groups.
    // note: in DEPN QUALSA is used to store the diff vs primary
    if (sam_is_prim_vb)
        ctx_create_node (VB, SAM_QUALSA, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_QUAL, '0'+(segconf.sam_ms_type == ms_BIOBAMBAM) }, 3);
}

// ZIP/PIZ: decompresses grp qual of grp, into vb->scratch
static void sam_get_sa_grp_qual (VBlockSAMP vb)
{
    ASSERTNOTINUSE (vb->scratch);

    buf_alloc (vb, &vb->scratch, vb->sa_grp->seq_len, 0, char, 1, "scratch");

    if (vb->sa_grp->qual_comp_len) { // qual of this group is compressed
        uint32_t uncomp_len = vb->sa_grp->seq_len;
        void *success = rans_uncompress_to_4x16 (VB, B8(z_file->sa_qual, vb->sa_grp->qual), vb->sa_grp->qual_comp_len,
                                                 B1ST(uint8_t, vb->scratch), &uncomp_len); 

        ASSERTGOTO (success && uncomp_len == vb->sa_grp->seq_len, "rans_uncompress_to_4x16 failed to decompress an SA Group QUAL data: vb_i=%u line_i=%d grp_i=%u success=%u comp_len=%u uncomp_len=%u expected_uncomp_len=%u qual=%"PRIu64,
                    vb->vblock_i, vb->line_i, ZGRP_I(vb->sa_grp), !!success, vb->sa_grp->qual_comp_len, uncomp_len, vb->sa_grp->seq_len, vb->sa_grp->qual);

        vb->scratch.len = uncomp_len;
    }
    else // not compressed
        buf_add (&vb->scratch, B8(z_file->sa_qual, vb->sa_grp->qual), vb->sa_grp->seq_len);

    return;
    
error:
    sam_show_sa_one_grp (ZGRP_I(vb->sa_grp));
    exit_on_error(true);
}

static void sam_seg_QUAL_diff_vs_primary (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual), unsigned add_bytes)
{
    if (vb->qual_missing) return;

    buf_alloc (vb, &CTX(SAM_QUALSA)->local, qual_len, 0, int8_t, CTX_GROWTH, "contexts->local");
    int8_t *diff = BAFT (int8_t, CTX(SAM_QUALSA)->local);
    CTX(SAM_QUALSA)->local.len += qual_len;

    bool xstrand = (dl->FLAG.bits.rev_comp != vb->sa_grp->revcomp);

    sam_get_sa_grp_qual (vb); // decompress qual into vb->scratch

    if (!xstrand) {
        rom grp_qual = Bc (vb->scratch, vb->hard_clip[0]);
        for (uint32_t i=0; i < qual_len; i++) 
            diff[i] = qual[i] - grp_qual[i]; // a value [-93,93] as qual is [33,126]
    }
    else {
        rom grp_qual = Bc (vb->scratch, vb->sa_grp->seq_len - 1 - vb->hard_clip[0]);
        for (int32_t/*signed*/ i=0; i < qual_len; i++) 
            diff[i] = qual[i] - grp_qual[-i]; // a value [-93,93] as qual is [33,126]
    }

    buf_free (vb->scratch);
}

void sam_seg_QUAL (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual_data)/*always textual*/, unsigned add_bytes)
{

    Context *qual_ctx = CTX(SAM_QUAL);

    // DEPN case of SA Group
    if (vb->sa_grp && sam_is_depn_vb && !vb->qual_missing) {
        dl->QUAL.len = 0; // don't compress this line
        sam_seg_QUAL_diff_vs_primary (vb, dl, STRa(qual_data), add_bytes);
        CTX(SAM_QUALSA)->txt_len += add_bytes;
    }
    // MAIN & PRIM: compress. 
    // Note: in PRIM, QUAL is not reconstructed (as QUAL is not in TOPLEVEL container) - it is consumed when loading SA Groups
    //       Instead, all-the-same QUALSA is reconstructed (SPECIAL copying from the SA Group)
    else {
        qual_ctx->local.len += dl->QUAL.len;
        qual_ctx->txt_len += add_bytes;
    }   

    // get QUAL score, consumed by buddy ms:i
    if (!segconf.running && segconf.sam_ms_type == ms_BIOBAMBAM)
        dl->QUAL_score = qual_ctx->last_value.i = sam_get_QUAL_score (vb, STRa(qual_data));

    // get stats on qual scores
    if (segconf.running)
        segconf_update_qual (STRa (qual_data));
}

// ms:i: (output of bamsormadup and other biobambam tools - ms in small letters), created here: https://github.com/gt1/libmaus/tree/master/src/libmaus/bambam/BamAlignmentDecoderBase.cpp getScore 
// It is the sum of phred values of mate's QUAL, but only phred values >= 15
void sam_seg_ms_field (VBlockSAMP vb, ValueType ms, unsigned add_bytes)
{
    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1
    
    if (vb->buddy_line_i != -1 && buddy_dl->QUAL_score == ms.i)  // successful in ~97% of lines with buddy
        seg_by_did_i (VB, STRa(QUAL_buddy_snip), OPTION_ms_i, add_bytes); // copy ms from earlier-line buddy 

    else
        seg_integer (VB, CTX(OPTION_ms_i), ms.i, true, add_bytes);    
}

//---------
// PIZ
//---------

static void sam_piz_QUAL_undiff_vs_primary (VBlockSAMP vb)
{
    const SAGroupType *g = vb->sa_grp;

    ContextP flag_ctx   = LOADED_CTX(SAM_FLAG);
    ContextP qualsa_ctx = LOADED_CTX(SAM_QUALSA);

    SamFlags FLAG = { .value = flag_ctx->last_value.i };
    bool xstrand = (FLAG.bits.rev_comp != g->revcomp);

    uint32_t qual_len = g->seq_len - vb->hard_clip[0] - vb->hard_clip[1];
    rom diff = Bc (qualsa_ctx->local, qualsa_ctx->next_local);
    qualsa_ctx->next_local += qual_len;

    char *qual = BAFTtxt;
    vb->txt_data.len += qual_len;

    // uncompress PRIM qual to vb->scratch
    sam_get_sa_grp_qual (vb);

    if (!xstrand) {
        rom prim_qual = Bc (vb->scratch, vb->hard_clip[0]);
        for (uint32_t i=0; i < qual_len; i++) 
            qual[i] = prim_qual[i] + diff[i];
    }
    else {
        rom prim_qual = Bc(vb->scratch, g->seq_len - 1 - vb->hard_clip[0]);
        for (int32_t/*signed*/ i=0; i < qual_len; i++)
            qual[i] = prim_qual[-i] + diff[i];
    }

    buf_free (vb->scratch);
}

static void sam_piz_QUAL_primary (VBlockSAMP vb)
{
    sam_get_sa_grp_qual (vb); // uncompress PRIM qual to vb->scratch

    RECONSTRUCT (vb->scratch.data, vb->scratch.len);
    buf_free (vb->scratch);
}

void sam_reconstruct_missing_quality (VBlockP vb, char c, bool reconstruct)
{
    if (reconstruct) RECONSTRUCT1 ('*');

    VB_SAM->qual_missing = (c==' ') ? QUAL_MISSING_STANDARD // sam_zip_qual re-wrote a '*' marking 'unavailable' as ' ' to avoid confusing with '*' as a valid quality score
                         /*c==127*/ : QUAL_MISSING_PYSAM;   // 127 written by bam_seg_txt_line to support pysam BAM reconstruction (introduced V14)
}

// called in case file contains ms:i, to calculate the QUAL_score after reconstructing QUAL,
// or to calculate DEPN qual from PRIM qual
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_QUAL)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    char *qual = BAFTtxt;
    bool calc_score = !snip_len/*up to v13 we had no snip*/ || *snip=='1';
    const SAGroupType *g = vb->sa_grp;

    // case: reconstruct by copying from SA Group
    if (SAM_PIZ_HAS_SA_GROUP) {
        if (!reconstruct) {}

        else if (g->no_qual) {
            vb->qual_missing = g->no_qual; // needed, for correct translation to BAM
            RECONSTRUCT1('*');
        }

        else if (sam_is_depn_vb) 
            sam_piz_QUAL_undiff_vs_primary (vb);                
        
        else  // primary vb
            sam_piz_QUAL_primary (vb);
    }
        
    // case: reconstruct from data in local
    else switch (ctx->ltype) { // the relevant subset of ltypes from reconstruct_from_ctx_do
        case LT_CODEC:
            codec_args[ctx->lcodec].reconstruct (VB, ctx->lcodec, ctx); break;

        case LT_SEQUENCE: 
            reconstruct_from_local_sequence (VB, ctx, NULL, 0, reconstruct); break;

        default: ASSPIZ (false, "Invalid ltype=%s for QUAL", lt_name (ctx->ltype));
    }

    uint32_t qual_len = BAFTtxt - qual;

    if (calc_score)
        new_value->i = sam_get_QUAL_score (vb, STRa(qual));

    return calc_score; // has new value
}

// SAM-to-BAM translator: translate SAM ASCII (33-based) Phred values to BAM's 0-based
TRANSLATOR_FUNC (sam_piz_sam2bam_QUAL)
{
    // if QUAL is "*" there are two options:
    // 1. If l_seq is 0, the QUAL is empty
    // 2. If not (i.e. we have SEQ data but not QUAL) - it is a string of 0xff, length l_seq
    if (VB_SAM->qual_missing) {
        BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)Btxt(vb->line_start);
        uint32_t l_seq = LTEN32 (alignment->l_seq);

        if (!l_seq) // option 1
            vb->txt_data.len--;
        
        else if (VB_SAM->qual_missing == QUAL_MISSING_STANDARD) {  // option 2 - missing QUAL according to SAM spec
            memset (BLSTtxt, 0xff, l_seq); // override the '*' and l_seq-1 more
            vb->txt_data.len += l_seq - 1;
        }

        else /* QUAL_MISSING_PYSAM */ {  // option 3 - missing QUAL as created by pysam (non-compliant with SAM spec)
            *BLSTtxt = 0xff;
            memset (BAFTtxt, 0, l_seq-1); // filler is 0 instead of 0xff as in SAM SPEC
            vb->txt_data.len += l_seq - 1;
        }
    }
    
    else // we have QUAL - update Phred values
        for (uint32_t i=0; i < recon_len; i++)
            recon[i] -= 33; 

    return 0;
}

// SAM-to-FASTQ translator: reverse the sequence if needed, and drop if "*"
TRANSLATOR_FUNC (sam_piz_sam2fastq_QUAL)
{
    uint16_t sam_flag = (uint16_t)vb->last_int(SAM_FLAG);
    
    // case: QUAL is "*" - don't show this fastq record
    if (recon_len==1 && *recon == '*') 
        vb->drop_curr_line = "no_qual";

    // case: this sequence is reverse complemented - reverse the QUAL string
    else if (sam_flag & SAM_FLAG_REV_COMP) {

        // we move from the outside in, switching the left and right bases 
        for (unsigned i=0; i < recon_len / 2; i++) 
            SWAP (recon[i], recon[recon_len-1-i]);
    }

    return 0;
}

