// ------------------------------------------------------------------
//   sam_qual.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted,
//   under penalties specified in the license.

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

#define QUAL_ZIP_CALLBACK(tag, f, may_be_revcomped)             \
COMPRESSOR_CALLBACK (sam_zip_##tag)                             \
{                                                               \
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);                 \
    *line_data_len = dl->dont_compress_##tag ? 0 : MIN_(maximum_size, dl->f.len); /* note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec */ \
    if (!line_data || ! *line_data_len) return; /* no dat, or only lengths were requested */   \
    *line_data = Btxt(dl->f.index);                             \
    if (is_rev) *is_rev = may_be_revcomped ? dl->FLAG.rev_comp : false;\
}                               

QUAL_ZIP_CALLBACK(OQ, OQ, true)
QUAL_ZIP_CALLBACK(TQ, TQ, true)
QUAL_ZIP_CALLBACK(GY, GY, false)
QUAL_ZIP_CALLBACK(2Y, _2Y, true)
QUAL_ZIP_CALLBACK(QX, solo_z_fields[SOLO_QX], false)
QUAL_ZIP_CALLBACK(CY, solo_z_fields[SOLO_CY], false)
QUAL_ZIP_CALLBACK(QT, solo_z_fields[SOLO_QT], false)

void sam_seg_QUAL_initialize (VBlockSAMP vb)
{
    if (sam_is_main_vb) 
        ctx_set_store_per_line (VB, SAM_QUAL, SAM_FLAG, DID_EOL);

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

// diff current QUAL against prim or saggy's QUAL
static void sam_seg_QUAL_diff (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual), 
                               STRp (prim_qual), bool prim_revcomp, uint32_t *prim_hard_clip, unsigned add_bytes)
{
    bool xstrand = (dl->FLAG.rev_comp != prim_revcomp);

    // flanks are the regions of qual that are not covered by prim_qual, due to prim_qual's hard clips.
    // example1: prim: HHHHQQQQQHH  (H=Hard clip, Q=some qual score)
    //           qual: HHQQQQQHHHH
    //                   ^^ flank[0]==2  flank[1]==0 (bc qual Q's don't go beyond prim's) overlap_len=3 (bc 3 Qs overlapping)
    //
    // example2: prim: HHHHHHHHQQH  (H=Hard clip, Q=some qual score)
    //           qual: HHQQQHHHHHH flank={3,0} overlap_len=0

    uint32_t flank[2] = { MIN_(qual_len, MAX_((int32_t)prim_hard_clip[xstrand]  - (int32_t)vb->hard_clip[0], 0)),   // left-flanking
                          MIN_(qual_len, MAX_((int32_t)prim_hard_clip[!xstrand] - (int32_t)vb->hard_clip[1], 0)) }; // right-flanking

    uint32_t overlap_len = qual_len - flank[0] - flank[1];

    CTX(SAM_QUAL)->txt_len += add_bytes - flank[0] - flank[1];
    CTX(SAM_QUAL_FLANK)->txt_len += flank[0] + flank[1];

    buf_alloc (vb, &CTX(SAM_QUALSA)->local, overlap_len, 0, int8_t, CTX_GROWTH, "contexts->local");
    int8_t *diff = BAFT (int8_t, CTX(SAM_QUALSA)->local);
    CTX(SAM_QUALSA)->local.len32 += overlap_len;

    rom overlap_this_qual = &qual[flank[0]];

    if (!xstrand) {
        uint32_t first_prim_overlap = flank[0] ? 0 : (vb->hard_clip[0] - prim_hard_clip[0]);
        rom overlap_prim_qual = &prim_qual[first_prim_overlap];

        for (uint32_t i=0; i < overlap_len; i++) 
            diff[i] = overlap_this_qual[i] - overlap_prim_qual[i]; // a value [-93,93] as qual is [33,126]
    }
    else {
        uint32_t last_prim_overlap = prim_qual_len - 1 - (flank[0] ? 0 : (vb->hard_clip[0] - prim_hard_clip[1]));
        rom overlap_prim_qual = &prim_qual[last_prim_overlap];

        for (int32_t/*signed*/ i=0; i < overlap_len; i++) 
            diff[i] = overlap_this_qual[i] - overlap_prim_qual[-i]; 
    }

    // add flanks (regions that could not be diffed) to SAM_QUAL_FLANK.local
    if (flank[0]) 
        seg_add_to_local_fixed (VB, CTX(SAM_QUAL_FLANK), qual, flank[0], LOOKUP_NONE, 0);

    if (flank[1])
        seg_add_to_local_fixed (VB, CTX(SAM_QUAL_FLANK), &qual[qual_len - flank[1]], flank[1], LOOKUP_NONE, 0);
}

void sam_seg_QUAL (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual_data)/*always textual*/, unsigned add_bytes)
{
    START_TIMER;

    Context *qual_ctx = CTX(SAM_QUAL);
    ZipDataLineSAM *saggy_dl;
    bool prim_has_qual_but_i_dont = false; // will be set if this line has no QUAL, but its prim line does (very rare)

    vb->has_qual |= !vb->qual_missing;

    // note: if prim (of either type) has no QUAL, we don't attempt to diff - bc piz wouldn't be able to know whether 
    // the current line has QUAL or not

    // case: DEPN component, line has a sag, and the saq has qual
    if (vb->sag && sam_is_depn_vb && !vb->sag->no_qual) {
        if (vb->qual_missing) {
            prim_has_qual_but_i_dont = true;
            CTX(SAM_QUALSA)->txt_len += add_bytes;
        }

        // case: both this line and prim have QUAL - diff them
        else {
            sam_get_sa_grp_qual (vb); // decompress prim qual into vb->scratch

            sam_seg_QUAL_diff (vb, dl, STRa(qual_data), STRb(vb->scratch), vb->sag->revcomp, (uint32_t[]){0, 0}, add_bytes);        
            buf_free (vb->scratch);
        }

        dl->dont_compress_QUAL = true; // don't compress this line
    }

    // case: MAIN component, the line has a corresponding saggy line in same VB, and saggy line has qual
    else if (sam_has_saggy && // note: saggy_line_i is set by sam_seg_saggy only for lines in MAIN component
             ({ saggy_dl = DATA_LINE (vb->saggy_line_i); 
                !saggy_dl->no_qual && 
                dl->hard_clip[0] + dl->hard_clip[1] + dl->SEQ.len == saggy_dl->hard_clip[0] + saggy_dl->hard_clip[1] + saggy_dl->SEQ.len; })) {

        if (vb->qual_missing) {
            prim_has_qual_but_i_dont = true;
            qual_ctx->txt_len += add_bytes;
        }    
        else 
            sam_seg_QUAL_diff (vb, dl, STRa(qual_data), STRtxtw(saggy_dl->QUAL), saggy_dl->FLAG.rev_comp, saggy_dl->hard_clip, add_bytes);
        
        dl->dont_compress_QUAL = true; // don't compress this line
    }

    // if we suspect entire file might be qual-less, seg as snip (without adding local) so that would become 
    // all-the-same in that. If we might have mixed qual/no-qual in the file, we add to local to not add entropy to b250
    else if (dl->no_qual && !segconf.nontrivial_qual) {
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_QUAL, '*' }, 3, SAM_QUAL, add_bytes); 
        goto done;
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

done:
    COPY_TIMER (sam_seg_QUAL);
}

//-----------------------------------------------------------------------------------------------------
// OQ:Z - Original QUAL (standard)
//        "Original base quality, usually before recalibration. Same encoding as QUAL" - https://samtools.github.io/hts-specs/SAMtags.pdf
// GY:Z - CellRanger
// 2Y:Z - CellRanger - R2 qual?
// TQ:Z - CellRanger & longranger: Quality values of the 7 trimmed bases following the barcode sequence at the start of R1. Can be used to reconstruct the original R1 quality values.
//-----------------------------------------------------------------------------------------------------
void sam_seg_other_qual (VBlockSAMP vb, TxtWord *dl_word, Did did_i, STRp(qual), bool len_is_seq_len, unsigned add_bytes)
{
    *dl_word = (TxtWord){ .index = BNUMtxt (qual), .len = qual_len };
    CTX(did_i)->local.len32 += qual_len;
    CTX(did_i)->txt_len     += add_bytes;

    if (!len_is_seq_len) 
        seg_lookup_with_length (VB, CTX(did_i), qual_len, 0);
}

//---------
// QUAL PIZ
//---------

static void sam_piz_QUAL_undiff_vs_primary (VBlockSAMP vb, STRp (prim_qual), bool prim_revcomp, const CigarAnalItem *saggy_anal, bool prim_is_bam, bool reconstruct)
{
    bool xstrand = (last_flags.rev_comp != prim_revcomp);

    uint32_t qual_len = prim_qual_len + (saggy_anal ? (saggy_anal->hard_clip[0] + saggy_anal->hard_clip[1]) : 0)
                      - (vb->hard_clip[0] + vb->hard_clip[1]);

    uint32_t flank[2] = { saggy_anal ? MIN_(qual_len, MAX_((int32_t)saggy_anal->hard_clip[xstrand]  - (int32_t)vb->hard_clip[0], 0)) : 0,   // left-flanking
                          saggy_anal ? MIN_(qual_len, MAX_((int32_t)saggy_anal->hard_clip[!xstrand] - (int32_t)vb->hard_clip[1], 0)) : 0 }; // right-flanking

    uint32_t overlap_len = qual_len - flank[0] - flank[1];

    rom diff = 0;
    if (overlap_len) {
        ContextP qualsa_ctx = LOADED_CTX(SAM_QUALSA);
        diff = Bc (qualsa_ctx->local, qualsa_ctx->next_local);
        qualsa_ctx->next_local += overlap_len;
    }

    char *qual = BAFTtxt;
    int8_t bam_bump = (prim_is_bam ? 33 : 0);

    uint32_t save_seq_len = vb->seq_len;

    if (flank[0]) {
        vb->seq_len = flank[0];
        reconstruct_from_local_sequence (VB, CTX(SAM_QUAL_FLANK), 0, 0, reconstruct); // reconstruct vb->seq_len quality scores
        qual += flank[0];
    }

    if (overlap_len && !xstrand && reconstruct) {
        uint32_t first_prim_overlap = flank[0] ? 0 : (vb->hard_clip[0] - (saggy_anal ? saggy_anal->hard_clip[0] : 0));
        prim_qual += first_prim_overlap;
        
        for (uint32_t i=0; i < overlap_len; i++) 
            qual[i] = prim_qual[i] + diff[i] + bam_bump;

        vb->txt_data.len32 += overlap_len;
    }

    else if (overlap_len && xstrand && reconstruct) {
        uint32_t last_prim_overlap = prim_qual_len - 1 - (flank[0] ? 0 : (vb->hard_clip[0] - (saggy_anal ? saggy_anal->hard_clip[1] : 0)));
        prim_qual += last_prim_overlap;

        for (int32_t/*signed*/ i=0; i < overlap_len; i++) {
            qual[i] = prim_qual[-i] + diff[i] + bam_bump;
}

        vb->txt_data.len32 += overlap_len;
    }

    if (flank[1]) {
        vb->seq_len = flank[1];
        reconstruct_from_local_sequence (VB, CTX(SAM_QUAL_FLANK), 0, 0, reconstruct);
    }

    vb->seq_len = save_seq_len;
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
    START_TIMER;

    VBlockSAMP vb = (VBlockSAMP)vb_;
    char *qual = BAFTtxt;
    const Sag *g = vb->sag;
    bool prim_has_qual_but_i_dont = (snip[0] == '1');
    const CigarAnalItem *saggy_anal;

    // case: reconstruct by copying from sag (except if we are depn and group has no qual)
    if (SAM_PIZ_HAS_SAG && (sam_is_prim_vb || (sam_is_depn_vb && !g->no_qual))) {
      
        if (!reconstruct) {}

        else if (g->no_qual || prim_has_qual_but_i_dont) 
            sam_reconstruct_missing_quality (VB, reconstruct);

        else if (sam_is_depn_vb) {
            sam_get_sa_grp_qual (vb); // uncompress PRIM qual to vb->scratch
            
            sam_piz_QUAL_undiff_vs_primary (vb, STRb(vb->scratch), vb->sag->revcomp, NULL, false, reconstruct);                
            buf_free (vb->scratch);
        }
        else  // primary vb
            sam_piz_QUAL_primary (vb);
    }

    // case: MAIN component, reconstruct depn line against prim line in this VB
    else if (sam_has_saggy && // only appears in the MAIN component and only since v14 (see sam_seg_saggy)
             ({ saggy_anal = B(CigarAnalItem, CTX(SAM_CIGAR)->cigar_anal_history, vb->saggy_line_i);
                // note: we've seen cases in the wild where a depn without hard clips is shorter than its prim (eg in test.NA12878.chr22.1x.bam), possibly due to GATK IndelRealigner
                vb->hard_clip[0] + vb->hard_clip[1] + vb->seq_len == saggy_anal->hard_clip[0] + saggy_anal->hard_clip[1] + saggy_anal->seq_len; })) {

        HistoryWord word = *B(HistoryWord, ctx->history, vb->saggy_line_i); // QUAL is always stored as LookupTxtData or LookupPerLine
        SamFlags saggy_flags = { .value = history64 (SAM_FLAG, vb->saggy_line_i) };
        rom saggy_qual = (word.lookup == LookupTxtData) ? Btxt(word.index) : Bc(ctx->per_line, word.index);

        if ((uint8_t)*saggy_qual == (IS_RECON_BAM ? 0xff : '*'))
            goto no_diff; // in case prim has no qual - seg did not diff (the current line may or may not have qual)

        else if (prim_has_qual_but_i_dont)
            sam_reconstruct_missing_quality (VB, reconstruct);

        else
            sam_piz_QUAL_undiff_vs_primary (vb, saggy_qual, word.len, saggy_flags.rev_comp, saggy_anal, flag.out_dt == DT_BAM, reconstruct);                
    } 
        
    // case: reconstruct from data in local
    else no_diff: 
        if (snip[0] == '*')  // case of seg suspecting entire file is qual-less
            sam_reconstruct_missing_quality (VB, reconstruct);

        else switch (ctx->ltype) { // the relevant subset of ltypes from reconstruct_from_ctx_do
            case LT_CODEC:
                codec_args[ctx->lcodec].reconstruct (VB, ctx->lcodec, ctx, NULL, 0); break;

            case LT_SEQUENCE: 
                reconstruct_from_local_sequence (VB, ctx, NULL, 0, reconstruct); break;

            default: ASSPIZ (false, "Invalid ltype=%s for %s", lt_name (ctx->ltype), ctx->tag_name);
        }

    uint32_t qual_len = BAFTtxt - qual;

#ifdef DEBUG
    for (uint32_t i=0; i < qual_len; i++)
        ASSPIZ (IS_NON_WS_PRINTABLE(qual[i]), "Invalid QUAL character reconstructed: i=%u char=%u\n", i, (uint8_t)qual[i]);
#endif

    // store quality score in ms:i history 
    if (!VER(14)) // up to v13 this special was only when calculating score was needed
        new_value->i = sam_get_QUAL_score (vb, STRa(qual));
    
    else if (VER(14) && segconf.sam_ms_type == ms_BIOBAMBAM) // note: segconf.sam_ms_type is populated in PIZ since v14
        *B(int64_t, CTX(OPTION_ms_i)->history, vb->line_i) = sam_get_QUAL_score (vb, STRa(qual));

    COPY_TIMER (sam_piz_special_QUAL);

    return !VER(14); // has new value for files up to v13
}

// SAM-to-BAM translator: translate SAM ASCII (33-based) Phred values to BAM's 0-based
TRANSLATOR_FUNC (sam_piz_sam2bam_QUAL)
{
    START_TIMER;

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

    COPY_TIMER (sam_piz_sam2bam_QUAL);
    return 0;
}

// SAM-to-FASTQ translator: reverse the sequence if needed, and drop if "*"
TRANSLATOR_FUNC (sam_piz_sam2fastq_QUAL)
{
    START_TIMER;

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

    COPY_TIMER (sam_piz_sam2fastq_QUAL);
    return 0;
}
