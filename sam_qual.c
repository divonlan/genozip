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
static uint32_t sam_get_QUAL_score (VBlockSAMP vb, STRp(qual))
{
    if (vb->qual_missing) return 0;

    uint32_t score=0;
    for (uint32_t i=0; i < qual_len; i++) {
        int8_t phred = qual[i] - 33;
        if (phred >= 15) score += phred;
    }
    
    return score;
}

// main thread (not thread-safe): called from sam_show_sag_one_grp for getting first few characters of alignment cigar
rom sam_display_qual_from_SA_Group (const Sag *g)
{
    if (g->no_qual) return "*";

    static char qual[SA_QUAL_DISPLAY_LEN+1];
    memset (qual, 0, sizeof(qual));

    uint32_t uncomp_len = MIN_(SA_CIGAR_DISPLAY_LEN, (uint32_t)g->seq_len); // possibly shorter than original cigar

    if (g->qual_comp_len) { // qual of this group is compressed
        uint32_t uncomp_len = MIN_((uint32_t)g->seq_len, SA_QUAL_DISPLAY_LEN);
        void *success = rans_uncompress_to_4x16 (evb, B8(z_file->sag_qual, g->qual), g->qual_comp_len,
                                                (uint8_t *)qual, &uncomp_len); 
        if (success && uncomp_len) qual[uncomp_len] = '\0';
    }

    else // not compressed
        memcpy (qual, B8(z_file->sag_qual, g->qual), uncomp_len);

    return qual;
}

//---------
// QUAL SEG
//---------

// callback function for compress to get data of one line
COMPRESSOR_CALLBACK (sam_zip_qual) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = dl->dont_compress_QUAL ? 0 : MIN_(maximum_size, dl->QUAL.len);

    if (!line_data) return; // only lengths were requested

    *line_data = Btxt(dl->QUAL.index);

    // if QUAL is just "*" (i.e. unavailable) replace it by " " because '*' is a legal PHRED quality value that will confuse PIZ
    if (dl->QUAL.len == 1 && (*line_data)[0] == '*') 
        *line_data = " "; // pointer to static string

    // note - we optimize just before compression - hopefully the string will remain in L1 cache
    else if (flag.optimize_QUAL) 
        optimize_phred_quality_string (STRa(*line_data));

    if (is_rev) *is_rev = dl->FLAG.rev_comp;
}

#define QUAL_ZIP_CALLBACK(tag)                          \
COMPRESSOR_CALLBACK (sam_zip_##tag)                     \
{                                                       \
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);         \
    *line_data_len  = dl->dont_compress_##tag ? 0 : MIN_(maximum_size, dl->tag.len); /* note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec */ \
    if (!line_data) return; /* only lengths were requested */ \
    *line_data = Btxt(dl->tag.index);                   \
    if (dl->OQ.len == 1 && (*line_data)[0] == '*') /* if OQ is just "*" (i.e. unavailable) replace it by " " (same as QUAL) */\
        *line_data = " "; /* pointer to static string */\
    if (is_rev) *is_rev = dl->FLAG.rev_comp;            \
}                               
QUAL_ZIP_CALLBACK(OQ)
QUAL_ZIP_CALLBACK(UY)

void sam_zip_QUAL_initialize (void)
{
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _OPTION_ms_i, false, 0, ms_mate_snip);
}

void sam_seg_QUAL_initialize (VBlockSAMP vb)
{
    CTX(SAM_QUAL)->ltype = LT_UINT8;

    if (sam_is_main_vb) 
        ctx_set_store_per_line (VB, 3, SAM_QUAL, SAM_FLAG, DID_EOL);

    // case: in DEPN and MAIN QUALSA is used to store the diff vs primary 
    if (sam_is_depn_vb || sam_is_main_vb)
        CTX(SAM_QUALSA)->ltype = LT_INT8; // diff of prim vs depn qual

    // case: PRIM - TOPLEVEL reconstructs all-the-same SAM_QUALSA instead of SAM_QUAL. SAM_QUAL is consumed when loading SA Groups.
    if (sam_is_prim_vb)
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_QUAL, '0' }, 3, SAM_QUALSA, 0); // note: we can't just ctx_create_node because we need to transfer flags to piz, so need b250

    if (segconf.sam_ms_type == ms_BIOBAMBAM) {      // handle QUAL_scores for ms:i - added v13        
        CTX(OPTION_ms_i)->flags.spl_custom = true;  // custom store-per-line - SPECIAL will handle the storing
        CTX(OPTION_ms_i)->flags.store = STORE_INT;  // since v14 - store QUAL_score for mate ms:i (in v13 it was stored in QUAL)
    }
}

// ZIP/PIZ: decompresses grp qual of grp, into vb->scratch
static void sam_get_sa_grp_qual (VBlockSAMP vb)
{
    ASSERTNOTINUSE (vb->scratch);

    buf_alloc (vb, &vb->scratch, vb->sag->seq_len, 0, char, 1, "scratch");

    if (vb->sag->qual_comp_len) { // qual of this group is compressed
        uint32_t uncomp_len = vb->sag->seq_len;
        void *success = rans_uncompress_to_4x16 (VB, B8(z_file->sag_qual, vb->sag->qual), vb->sag->qual_comp_len,
                                                 B1ST(uint8_t, vb->scratch), &uncomp_len); 

        ASSERTGOTO (success && uncomp_len == vb->sag->seq_len, "%s: rans_uncompress_to_4x16 failed to decompress an SA Group QUAL data: grp_i=%u success=%u comp_len=%u uncomp_len=%u expected_uncomp_len=%u qual=%"PRIu64,
                    LN_NAME, ZGRP_I(vb->sag), !!success, vb->sag->qual_comp_len, uncomp_len, vb->sag->seq_len, (uint64_t)vb->sag->qual);

        vb->scratch.len = uncomp_len;
    }
    else // not compressed
        buf_add (&vb->scratch, B8(z_file->sag_qual, vb->sag->qual), vb->sag->seq_len);

    return;
    
error:
    sam_show_sag_one_grp (ZGRP_I(vb->sag));
    exit_on_error(true);
}

static void sam_seg_QUAL_diff_vs_primary (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual), 
                                          STRp (prim_qual), bool prim_revcomp)
{
    buf_alloc (vb, &CTX(SAM_QUALSA)->local, qual_len, 0, int8_t, CTX_GROWTH, "contexts->local");
    int8_t *diff = BAFT (int8_t, CTX(SAM_QUALSA)->local);
    CTX(SAM_QUALSA)->local.len32 += qual_len;

    bool xstrand = (dl->FLAG.rev_comp != prim_revcomp);

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
    START_TIMER;

    Context *qual_ctx = CTX(SAM_QUAL);
    ZipDataLineSAM *prim_dl;
    bool prim_has_qual_but_i_dont = false; // will be set if this line has no QUAL, but its prim line does (very rare)

    vb->has_qual |= !vb->qual_missing;

    // note: if prim (of either type) has no QUAL, we don't attempt to diff - bc piz wouldn't be able to know whether 
    // the current line has QUAL or not

    // case: DEPN component, line has SA Group, and SA Group as qual
    if (vb->sag && sam_is_depn_vb && !vb->sag->no_qual) {
        if (vb->qual_missing)
            prim_has_qual_but_i_dont = true;

        // case: both this line and prim have QUAL - diff them
        else {
            sam_get_sa_grp_qual (vb); // decompress prim qual into vb->scratch

            sam_seg_QUAL_diff_vs_primary (vb, dl, STRa(qual_data), STRb(vb->scratch), vb->sag->revcomp);        
            buf_free (vb->scratch);
        }

        dl->dont_compress_QUAL = true; // don't compress this line
        CTX(SAM_QUALSA)->txt_len += add_bytes;
    }

    // case: MAIN component, depn line corresponding prim line in same VB, and prim line has qual
    else if (zip_has_prim && // note: prim_line_i is set by sam_seg_QNAME only for depn (with exceptions for STAR) lines in MAIN component
             ({ prim_dl = DATA_LINE (vb->prim_line_i); !prim_dl->no_qual; })) {

        if (vb->qual_missing) 
            prim_has_qual_but_i_dont = true;
        
        else 
            sam_seg_QUAL_diff_vs_primary (vb, dl, STRa(qual_data), STRtxtw(prim_dl->QUAL), prim_dl->FLAG.rev_comp);
        
        dl->dont_compress_QUAL = true; // don't compress this line
        qual_ctx->txt_len += add_bytes;
    }

    // case: standard
    // Note: in PRIM, QUAL is not reconstructed (as QUAL is not in TOPLEVEL container) - it is consumed when loading SA Groups
    //       Instead, all-the-same QUALSA is reconstructed (SPECIAL copying from the SA Group)
    else {
        qual_ctx->local.len32 += dl->QUAL.len;
        qual_ctx->txt_len     += add_bytes;
    }   

    // seg SPECIAL. note: prim this is all-the-same and segged in sam_seg_QUAL_initialize
    if (!sam_is_prim_vb)
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_QUAL, '0' + prim_has_qual_but_i_dont }, 3, SAM_QUAL, 0); 

    // get QUAL score, consumed by mate ms:i
    if (!segconf.running && segconf.sam_ms_type == ms_BIOBAMBAM)
        dl->QUAL_score = sam_get_QUAL_score (vb, STRa(qual_data));
 
    // get stats on qual scores
    if (segconf.running)
        segconf_update_qual (STRa (qual_data));

    COPY_TIMER (sam_seg_QUAL);
}

//-----------------------------------------------------------------------------------------------------
// OQ:Z - "Original QUAL" 
// "Original base quality, usually before recalibration. Same encoding as QUAL" - https://samtools.github.io/hts-specs/SAMtags.pdf
//-----------------------------------------------------------------------------------------------------
void sam_seg_OQ_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(oq), unsigned add_bytes)
{
    dl->OQ = (TxtWord){ .index = BNUMtxt (oq), .len = oq_len };
    CTX(OPTION_OQ_Z)->local.len32 += oq_len;
    CTX(OPTION_OQ_Z)->txt_len     += add_bytes;
}

//---------
// QUAL PIZ
//---------

static void sam_piz_QUAL_undiff_vs_primary (VBlockSAMP vb, STRp (prim_qual), bool prim_revcomp, bool prim_is_bam)
{
    ContextP qualsa_ctx = LOADED_CTX(SAM_QUALSA);

    bool xstrand = (last_flags.rev_comp != prim_revcomp);

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

void sam_reconstruct_missing_quality (VBlockP vb, bool reconstruct)
{
    if (reconstruct) 
        RECONSTRUCT1 ('*');

    VB_SAM->qual_missing = true;
}

// Note: in PRIM, it is called with ctx=QUALSA, in MAIN and DEPN with ctx=QUAL
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_QUAL)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    char *qual = BAFTtxt;
    const Sag *g = vb->sag;
    bool prim_has_qual_but_i_dont = snip[0] - '0';

    // case: reconstruct by copying from SA Group (except if we are depn and group has no qual)
    if (SAM_PIZ_HAS_SA_GROUP && (sam_is_prim_vb || (sam_is_depn_vb && !g->no_qual))) {
        if (!reconstruct) {}

        else if (g->no_qual || prim_has_qual_but_i_dont) 
            sam_reconstruct_missing_quality (VB, reconstruct);

        else if (sam_is_depn_vb) {
            sam_get_sa_grp_qual (vb); // uncompress PRIM qual to vb->scratch
            
            sam_piz_QUAL_undiff_vs_primary (vb, STRb(vb->scratch), vb->sag->revcomp, false);                
            buf_free (vb->scratch);
        }
        else  // primary vb
            sam_piz_QUAL_primary (vb);
    }

    // case: MAIN component, reconstruct depn line against prim line in this VB
    else if (sam_is_main_vb && 
             VER(14) && // up to v13, we could have buddy lines for sup/sec which are not prim lines
             piz_has_prim) {

        HistoryWord word = *B(HistoryWord, ctx->history, vb->buddy_line_i); // QUAL is always stored as LookupTxtData or LookupPerLine
        SamFlags prim_flags = { .value = history64 (SAM_FLAG, vb->buddy_line_i) };
        rom prim_qual = (word.lookup == LookupTxtData) ? Btxt(word.index) : Bc(ctx->per_line, word.index);

        if ((uint8_t)*prim_qual == (IS_RECON_BAM ? 0xff : '*'))
            goto no_diff; // in case prim has no qual - seg did not diff (the current line may or may not have qual)

        else if (prim_has_qual_but_i_dont)
            sam_reconstruct_missing_quality (VB, reconstruct);

        else
            sam_piz_QUAL_undiff_vs_primary (vb, prim_qual, word.len, prim_flags.rev_comp, flag.out_dt == DT_BAM);                
    } 
        
    // case: reconstruct from data in local
    else no_diff: 
        switch (ctx->ltype) { // the relevant subset of ltypes from reconstruct_from_ctx_do
            case LT_CODEC:
                codec_args[ctx->lcodec].reconstruct (VB, ctx->lcodec, ctx, NULL, 0); break;

            case LT_SEQUENCE: 
                reconstruct_from_local_sequence (VB, ctx, NULL, 0, reconstruct); break;

            default: ASSPIZ (false, "Invalid ltype=%s for QUAL", lt_name (ctx->ltype));
        }

    uint32_t qual_len = BAFTtxt - qual;

    // store quality score in ms:i history 
    if (!VER(14)) // up to v13 this special was only when calculating score was needed
        new_value->i = sam_get_QUAL_score (vb, STRa(qual));
    
    else if (VER(14) && segconf.sam_ms_type == ms_BIOBAMBAM) // note: segconf.sam_ms_type is populated in PIZ since v14
        *B(int64_t, CTX(OPTION_ms_i)->history, vb->line_i) = sam_get_QUAL_score (vb, STRa(qual));

    return !VER(14); // has new value for files up to v13
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
        
        else if (VB_SAM->qual_missing && !segconf.pysam_qual) {  // option 2 - missing QUAL according to SAM spec
            memset (BLSTtxt, 0xff, l_seq); // override the '*' and l_seq-1 more
            vb->txt_data.len += l_seq - 1;
        }

        else {  // option 3 - missing QUAL as created by pysam (non-compliant with SAM spec)
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
void sam_seg_ms_i (VBlockSAMP vb, ValueType ms, unsigned add_bytes)
{
    ZipDataLineSAM *mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is NO_LINE

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, (MultiplexerP)&vb->mux_ms, zip_has_mate);

    if (zip_has_mate && mate_dl->QUAL_score == ms.i)     // successful in ~97% of lines with mate
        seg_by_ctx (VB, STRa(ms_mate_snip), channel_ctx, add_bytes); // note: prior to v14, we stored ms qual_score in QUAL history, not in ms:i history 
    else {
        channel_ctx->ltype = LT_DYN_INT; 
        seg_integer (VB, channel_ctx, ms.i, true, add_bytes);    
    }

    seg_by_did (VB, STRa(vb->mux_ms.snip), OPTION_ms_i, 0); // de-multiplexor
}
