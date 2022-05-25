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

static char ms_mate_snip[30];
static unsigned ms_mate_snip_len;

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
// QUAL SEG
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

    // note - we optimize just before compression - hopefully the string will remain in L1 cache
    else if (flag.optimize_QUAL) 
        optimize_phred_quality_string (STRa(*line_data));

    if (is_rev) *is_rev = dl->FLAG.bits.rev_comp;
}

void sam_zip_QUAL_initialize (void)
{
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _OPTION_ms_i, false, 0, ms_mate_snip);
}

void sam_seg_QUAL_initialize (VBlockSAMP vb)
{
    CTX(SAM_QUAL)->ltype = LT_UINT8;

    if (sam_is_main_vb) {
        CTX(SAM_QUAL)->flags.store_per_line = true; // for undiffing depn lines against prim (in same VB)
        CTX(SAM_FLAG)->flags.store_per_line = true; // also needed for undiffing
    }

    if (sam_is_depn_vb || sam_is_main_vb)
        CTX(SAM_QUALSA)->ltype = LT_INT8;           // diff of prim vs depn qual

    // if PRIM - TOPLEVEL reconstructs all-the-same SAM_QUALSA instead of SAM_QUAL. SAM_QUAL is consumed when loading SA Groups.
    // note: in DEPN and MAIN QUALSA is used to store the diff vs primary (note: we can't just ctx_create_node because we need to transfer flags to piz, so need b250)
    seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_QUAL }, 2, sam_is_prim_vb ? SAM_QUALSA : SAM_QUAL, 0);

    if (segconf.sam_ms_type == ms_BIOBAMBAM) {      // handle QUAL_scores for ms:i - added v13        
        CTX(OPTION_ms_i)->flags.spl_custom = true;  // custom store-per-line - SPECIAL will handle the storing
        CTX(OPTION_ms_i)->flags.store = STORE_INT;  // since v14 - store QUAL_score for mate ms:i (in v13 it was stored in QUAL)
    }
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

        ASSERTGOTO (success && uncomp_len == vb->sa_grp->seq_len, "%s: rans_uncompress_to_4x16 failed to decompress an SA Group QUAL data: grp_i=%u success=%u comp_len=%u uncomp_len=%u expected_uncomp_len=%u qual=%"PRIu64,
                    LN_NAME, ZGRP_I(vb->sa_grp), !!success, vb->sa_grp->qual_comp_len, uncomp_len, vb->sa_grp->seq_len, vb->sa_grp->qual);

        vb->scratch.len = uncomp_len;
    }
    else // not compressed
        buf_add (&vb->scratch, B8(z_file->sa_qual, vb->sa_grp->qual), vb->sa_grp->seq_len);

    return;
    
error:
    sam_show_sa_one_grp (ZGRP_I(vb->sa_grp));
    exit_on_error(true);
}

static void sam_seg_QUAL_diff_vs_primary (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual), 
                                          STRp (prim_qual), bool prim_revcomp, 
                                          unsigned add_bytes)
{
    if (vb->qual_missing) return;

    buf_alloc (vb, &CTX(SAM_QUALSA)->local, qual_len, 0, int8_t, CTX_GROWTH, "contexts->local");
    int8_t *diff = BAFT (int8_t, CTX(SAM_QUALSA)->local);
    CTX(SAM_QUALSA)->local.len += qual_len;

    bool xstrand = (dl->FLAG.bits.rev_comp != prim_revcomp);

    if (!xstrand) {
        rom grp_qual = &prim_qual[vb->hard_clip[0]];
        for (uint32_t i=0; i < qual_len; i++) 
            diff[i] = qual[i] - grp_qual[i]; // a value [-93,93] as qual is [33,126]
    }
    else {
        rom grp_qual = &prim_qual[prim_qual_len - 1 - vb->hard_clip[0]];
        for (int32_t/*signed*/ i=0; i < qual_len; i++) 
            diff[i] = qual[i] - grp_qual[-i]; // a value [-93,93] as qual is [33,126]
    }
}

void sam_seg_QUAL (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual_data)/*always textual*/, unsigned add_bytes)
{
    Context *qual_ctx = CTX(SAM_QUAL);
    ZipDataLineSAM *prim_dl;

    vb->has_qual |= !vb->qual_missing;

    // case: DEPN component, line has SA Group
    if (vb->sa_grp && sam_is_depn_vb && !vb->qual_missing) {
        dl->QUAL.len = 0; // don't compress this line

        sam_get_sa_grp_qual (vb); // decompress prim qual into vb->scratch

        sam_seg_QUAL_diff_vs_primary (vb, dl, STRa(qual_data), STRb(vb->scratch), vb->sa_grp->revcomp, add_bytes);
        CTX(SAM_QUALSA)->txt_len += add_bytes;

        buf_free (vb->scratch);
    }

    // case: MAIN component, depn line corresponding prim line in same VB (i.e. not sorted, or line failed the test to move to gc)
    else if (zip_has_prim && // note: prim_line_i is set by sam_seg_QNAME only for depn lines in MAIN component
             ({ prim_dl = DATA_LINE (vb->prim_line_i); prim_dl->QUAL.len; })) {
        dl->QUAL.len = 0; // don't compress this line

        sam_seg_QUAL_diff_vs_primary (vb, dl, STRa(qual_data), STRtxtw(prim_dl->QUAL), prim_dl->FLAG.bits.rev_comp, add_bytes);
        qual_ctx->txt_len += add_bytes;
    }

    // MAIN & PRIM: compress. 
    // Note: in PRIM, QUAL is not reconstructed (as QUAL is not in TOPLEVEL container) - it is consumed when loading SA Groups
    //       Instead, all-the-same QUALSA is reconstructed (SPECIAL copying from the SA Group)
    else {
        qual_ctx->local.len += dl->QUAL.len;
        qual_ctx->txt_len   += add_bytes;
    }   

    // get QUAL score, consumed by mate ms:i
    if (!segconf.running && segconf.sam_ms_type == ms_BIOBAMBAM)
        dl->QUAL_score = /*xxx qual_ctx->last_value.i =*/ sam_get_QUAL_score (vb, STRa(qual_data));
 
    // get stats on qual scores
    if (segconf.running)
        segconf_update_qual (STRa (qual_data));
}

//---------
// QUAL PIZ
//---------

static void sam_piz_QUAL_undiff_vs_primary (VBlockSAMP vb, STRp (prim_qual), bool prim_revcomp, bool prim_is_bam)
{
    ContextP qualsa_ctx = LOADED_CTX(SAM_QUALSA);

    bool xstrand = (last_flags.bits.rev_comp != prim_revcomp);

    uint32_t qual_len = prim_qual_len - vb->hard_clip[0] - vb->hard_clip[1];
    rom diff = Bc (qualsa_ctx->local, qualsa_ctx->next_local);
    qualsa_ctx->next_local += qual_len;

    char *qual = BAFTtxt;
    vb->txt_data.len32 += qual_len;
    int8_t bam_bump = (prim_is_bam ? 33 : 0);

    if (!xstrand) {
        prim_qual += vb->hard_clip[0];
        for (uint32_t i=0; i < qual_len; i++) 
            qual[i] = prim_qual[i] + diff[i] + bam_bump;
    }
    else {
        prim_qual += prim_qual_len - 1 - vb->hard_clip[0];
        for (int32_t/*signed*/ i=0; i < qual_len; i++)
            qual[i] = prim_qual[-i] + diff[i] + bam_bump;
    }
}

static void sam_piz_QUAL_primary (VBlockSAMP vb)
{
    sam_get_sa_grp_qual (vb); // uncompress PRIM qual to vb->scratch

    RECONSTRUCT (vb->scratch.data, vb->scratch.len);
    buf_free (vb->scratch);
}

void sam_reconstruct_missing_quality (VBlockP vb, char c, bool reconstruct)
{
    ASSERT (c==' ' || c==127, "%s: Expecting c=%u to 32 or 127", LN_NAME, c);

    if (reconstruct) RECONSTRUCT1 ('*');

    VB_SAM->qual_missing = (c==' ') ? QUAL_MISSING_STANDARD // sam_zip_qual re-wrote a '*' marking 'unavailable' as ' ' to avoid confusing with '*' as a valid quality score
                         /*c==127*/ : QUAL_MISSING_PYSAM;   // 127 written by bam_seg_txt_line to support pysam BAM reconstruction (introduced V14)
}

// Note: in PRIM, it is called with ctx=QUALSA, in MAIN and DEPN with ctx=QUAL
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_QUAL)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    char *qual = BAFTtxt;
    const SAGroupType *g = vb->sa_grp;

    // case: reconstruct by copying from SA Group 
    if (SAM_PIZ_HAS_SA_GROUP) {
        if (!reconstruct) {}

        else if (g->no_qual) 
            RECONSTRUCT1('*');

        else if (sam_is_depn_vb) {
            sam_get_sa_grp_qual (vb); // uncompress PRIM qual to vb->scratch
            
            sam_piz_QUAL_undiff_vs_primary (vb, STRb(vb->scratch), vb->sa_grp->revcomp, false);                
            buf_free (vb->scratch);
        }
        else  // primary vb
            sam_piz_QUAL_primary (vb);
    }

    // case: MAIN component, reconstruct depn line against prim line in this VB
    else if (sam_is_main_vb && 
             z_file->genozip_version >= 14 && // up to v13, we could have buddy lines for sup/sec which are not prim lines
             piz_has_buddy && (last_flags.bits.supplementary || last_flags.bits.secondary)) {

        HistoryWord word = *B(HistoryWord, ctx->history, vb->buddy_line_i); // QUAL is always stored as LookupTxtData or LookupPerLine
        SamFlags prim_flags = { .value = *B(int64_t, CTX(SAM_FLAG)->history, vb->buddy_line_i) };

        sam_piz_QUAL_undiff_vs_primary (vb, 
                                        (word.lookup == LookupTxtData) ? Btxt(word.index) : Bc(ctx->per_line, word.index),
                                        word.len, prim_flags.bits.rev_comp, flag.out_dt == DT_BAM);                
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

    // store quality score in ms:i history 
    if (z_file->genozip_version <= 13) // up to v13 this special was only when calculating score was needed
        new_value->i = sam_get_QUAL_score (vb, STRa(qual));
    
    else if (z_file->genozip_version >= 14 && segconf.sam_ms_type == ms_BIOBAMBAM) // note: segconf.sam_ms_type is populated in PIZ since v14
        *B(int64_t, CTX(OPTION_ms_i)->history, vb->line_i) = sam_get_QUAL_score (vb, STRa(qual));

    //xxx if (reconstruct) // note: At this, this is only used for translating to BAM, so not needed if not reconstructed
    //     vb->qual_missing = (BAFTtxt - qual == 1 && *qual == '*'); 

    return z_file->genozip_version <= 13; // has new value for files up to v13
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

//---------
// ms:i
//---------

// ms:i: (output of bamsormadup and other biobambam tools - ms in small letters), created here: https://github.com/gt1/libmaus/tree/master/src/libmaus/bambam/BamAlignmentDecoderBase.cpp getScore 
// It is the sum of phred values of mate's QUAL, but only phred values >= 15
void sam_seg_ms_field (VBlockSAMP vb, ValueType ms, unsigned add_bytes)
{
    ZipDataLineSAM *mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is NO_LINE

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, (MultiplexerP)&vb->mux_ms, zip_has_mate);

    if (zip_has_mate && mate_dl->QUAL_score == ms.i)     // successful in ~97% of lines with mate
        seg_by_ctx (VB, STRa(ms_mate_snip), channel_ctx, add_bytes); // note: prior to v14, we stored ms qual_score in QUAL history, not in ms:i history 
    else {
        channel_ctx->dynamic_size_local = true; 
        seg_integer (VB, channel_ctx, ms.i, true, add_bytes);    
    }

    seg_by_did_i (VB, STRa(vb->mux_ms.snip), OPTION_ms_i, 0); // de-multiplexor
}
