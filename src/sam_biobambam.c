// ------------------------------------------------------------------
//   sam_biobambam.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "seg.h"

// get a score of the QUAL string - similar to that calculated by biobambam for ms:i and by samtools fixmate -m
// see here: https://github.com/gt1/libmaus/tree/master/src/libmaus/bambam/BamAlignmentDecoderBase.cpp getScore 
// and here: calc_mate_score in samtools source.
// The sum of phred values of the QUAL string, but only for phred values >= 15
uint32_t sam_get_QUAL_score (VBlockSAMP vb, STRp(qual))
{
    if (vb->qual_missing) return 0;

    uint32_t score=0;
    for (uint32_t i=0; i < qual_len; i++) {
        int8_t phred = qual[i] - 33;
        if (phred >= 15) score += phred; // same threadhold used by both samtools and biobambam
    }
    
    return score;
}

// ----------------------------------------------------------------------------------------------
// ms:i: (output of bamsormadup and other biobambam tools - ms in small letters), created here: https://github.com/gt1/libmaus/tree/master/src/libmaus/bambam/BamAlignmentDecoderBase.hpp getScore 
//       also output of samtools fixmate -m: See calc_mate_score in samtools source.
//       It is the sum of phred values of mate's QUAL, but only phred values >= 15
// ----------------------------------------------------------------------------------------------
void sam_seg_ms_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t ms, unsigned add_bytes)
{
    int32_t save = dl->QUAL_score;
    sam_seg_buddied_i_fields (vb, dl, OPTION_ms_i, ms, &dl->QUAL_score, (MultiplexerP)&vb->mux_ms, STRa(copy_mate_ms_snip), add_bytes);
    dl->QUAL_score = save;
}

// mc:i: (output of bamsormadup and other biobambam tools - mc in small letters) 
// appears to be a pos value usually close to PNEXT, but it is -1 is POS=PNEXT. Note: tried muxing by POS=PNEXT, but was worse off
// from bamsort manual: "adddupmarksupport=<0|1>: add information required for streaming duplicate marking in the aux fields MS and MC.
// Input is assumed to be collated by query name. This option is ignored unless fixmates=1. By default it is disabled."
// https://github.com/gt1/biobambam2/blob/master/src/programs/bamsort.cpp says: "biobambam used MC as a mate coordinate tag which now has a clash
// with the official SAM format spec.  New biobambam version uses mc."
// ms="MateBaseScore" - sum all characters in QUAL of the mate, where the value of each character is its ASCII minus 33 (i.e. the Phred score)
// mc="MateCoordinate"
void sam_seg_mc_i (VBlockSAMP vb, int64_t mc, unsigned add_bytes)
{
    // if snip is "-1", store as simple snip
    if (mc == -1)
        seg_by_did (VB, "-1", 2, OPTION_mc_i, add_bytes);
    
    // delta vs PNEXT
    else
        seg_pos_field (VB, OPTION_mc_i, SAM_PNEXT, SPF_BAD_SNIPS_TOO, 0, 0, 0, mc, add_bytes);
}
