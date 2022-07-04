// ------------------------------------------------------------------
//   sam_bowtie2.c
//   Copyright (C) 2022-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "chrom.h"
#include "codec.h"

// ----------------------------------------------------------------------------------------------
// YS:i mate alignment score (bowtie2 and hisat2)
// ----------------------------------------------------------------------------------------------

bool sam_has_bowtie2_YS_i (void)
{
    return MP(BOWTIE2) || MP(BSSEEKER2) || MP(HISAT2);
}

void sam_seg_bowtie2_YS_i (VBlockSAMP vb, ZipDataLineSAM *dl, ValueType YS, unsigned add_bytes)
{
    ASSERT (YS.i >= MIN_AS_i && YS.i <= MAX_AS_i, "%s: YS=%"PRId64" is out of range [%d,%d]", LN_NAME, YS.i, MIN_AS_i, MAX_AS_i);    

    ctx_set_last_value (VB, CTX (OPTION_YS_i), YS);
    dl->YS = YS.i;

    ZipDataLineSAM *mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is -1

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, (MultiplexerP)&vb->mux_YS, zip_has_mate);

    if (zip_has_mate && mate_dl->AS == YS.i) 
        seg_by_ctx (VB, STRa(copy_buddy_AS_snip), channel_ctx, add_bytes);

    else 
        seg_integer_as_text_do (VB, channel_ctx, YS.i, add_bytes);    

    seg_by_did (VB, STRa(vb->mux_YS.snip), OPTION_YS_i, 0); // de-multiplexor
}
