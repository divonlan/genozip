// ------------------------------------------------------------------
//   sam_qual.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
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

// get a score of the QUAL string - similar to that calculated by biobambam for ms:i:
// see here: https://github.com/gt1/libmaus/tree/master/src/libmaus/bambam/BamAlignmentDecoderBase.cpp getScore 
// The sum of phred values of the QUAL string, but only for phred values >= 15
static int64_t sam_get_QUAL_score (STRp(qual))
{
    int64_t score=0;
    for (unsigned i=0; i < qual_len; i++) {
        char phred = qual[i] - 33;
        if (phred >= 15) score += phred;
    }
    
    return score;
}

//---------
// SEG
//---------

// callback function for compress to get data of one line (called by codec_bz2_compress)
COMPRESSOR_CALLBACK (sam_zip_qual) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_data_len  = MIN_(maximum_size, dl->QUAL.snip_len);

    if (!line_data) return; // only lengths were requested

    *line_data = ENT (char, vb->txt_data, dl->QUAL.char_index);

    // if QUAL is just "*" (i.e. unavailable) replace it by " " because '*' is a legal PHRED quality value that will confuse PIZ
    if (dl->QUAL.snip_len == 1 && (*line_data)[0] == '*') 
        *line_data = " "; // pointer to static string

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    else if (flag.optimize_QUAL) 
        optimize_phred_quality_string (STRa(*line_data));
}

void sam_seg_QUAL_initialize (VBlockP vb)
{
    if (segconf.has_ms) {     // handle QUAL_scores for ms:i - added v13
        CTX(SAM_QUAL)->flags.store_per_line = true;
        CTX(SAM_QUAL)->no_stons = true; // since we're storing QUAL data in local

        seg_by_did_i (vb, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_QUAL }, 2, SAM_QUAL, 0); // all-the-same entry for calculating QUAL scores
    }
}

void sam_seg_QUAL (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(qual_data), unsigned add_bytes)
{
    dl->QUAL = WORD_IN_TXT_DATA (qual_data);

    Context *qual_ctx = CTX(SAM_QUAL);
    qual_ctx->local.len += dl->QUAL.snip_len;
    qual_ctx->txt_len   += add_bytes;

    // get QUAL score, consumed by buddy ms:i
    if (!segconf.running && segconf.has_ms)
        dl->QUAL_score = qual_ctx->last_value.i = sam_get_QUAL_score (STRa(qual_data));
}

// ms:i: (output of bamsormadup and other biobambam tools - ms in small letters), created here: https://github.com/gt1/libmaus/tree/master/src/libmaus/bambam/BamAlignmentDecoderBase.cpp getScore 
// It is the sum of phred values of mate's QUAL, but only phred values >= 15
void sam_seg_ms_field (VBlockSAM *vb, DictId dict_id, STRp(ms), unsigned add_bytes)
{
    if (segconf.running) segconf.has_ms = true;

    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1
    
    int64_t ms_value;
    if (!str_get_int (STRa(ms), &ms_value)) goto fallback;

    if (vb->buddy_line_i != -1 && buddy_dl->QUAL_score == ms_value)  // successful in ~97% of lines with buddy
        seg_by_did_i (VB, STRa(QUAL_buddy_snip), OPTION_ms_i, add_bytes); // copy MQ from earlier-line buddy 

    else
fallback:
        seg_by_did_i (VB, STRa(ms), OPTION_ms_i, add_bytes);    
}

//---------
// PIZ
//---------

// called in case file contains ms:i, to calculate the QUAL_score after reconstructing QUAL
SPECIAL_RECONSTRUCTOR (sam_piz_special_QUAL)
{
    const char *qual = AFTERENT (char, vb->txt_data);

    switch (ctx->ltype) { // the relevant subset of ltypes from reconstruct_from_ctx_do
        case LT_CODEC:
            codec_args[ctx->lcodec].reconstruct (vb, ctx->lcodec, ctx); break;

        case LT_SEQUENCE: 
            reconstruct_from_local_sequence (vb, ctx, NULL, 0); break;

        default: ASSPIZ (false, "Invalid ltype=%s for QUAL", lt_name (ctx->ltype));
    }

    new_value->i = sam_get_QUAL_score (qual, AFTERENT (char, vb->txt_data) - qual);
    return true; // new value
}

// SAM-to-BAM translator: translate SAM ASCII (33-based) Phread values to BAM's 0-based
TRANSLATOR_FUNC (sam_piz_sam2bam_QUAL)
{
    // if QUAL is "*" there are two options:
    // 1. If l_seq is 0, the QUAL is empty
    // 2. If not (i.e. we have SEQ data but not QUAL) - it is a string of 0xff, length l_seq
    if (recon_len==1 && *recon == '*') {
        BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)ENT (char, vb->txt_data, vb->line_start);
        uint32_t l_seq = LTEN32 (alignment->l_seq);

        if (!l_seq) // option 1
            vb->txt_data.len--;
        else {      // option 2
            memset (LASTENT (uint8_t, vb->txt_data), 0xff, l_seq); // override the '*' and l_seq-1 more
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
