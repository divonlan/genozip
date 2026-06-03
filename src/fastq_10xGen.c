// ------------------------------------------------------------------
//   fastq_10xGen.c
//   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"
#include "codec.h"

#define BARCODE_LEN 16
#define UMI_LEN     12

void fastq_10xGen_segconf_finalize (VBlockP vb)
{
    // case: verified as 10xGenomics non-biological file: prepare container
    if (IS_NONBIO(10xGen) && segconf.std_seq_len == 28 && 
        percent (vb->num_aligned, vb->lines.len) <= 10) { // expected for non-biological file: aligner hardly worked
        
        static const Container(2) con = {
            .nitems_lo = 2,
            .repeats   = 1,
            .items = {
                { .dict_id.num = _FASTQ_NONBIO_BC0, .separator = { CI0_FIXED_0_PAD, BARCODE_LEN } }, // barcode
                { .dict_id.num = _FASTQ_NONBIO_UMI, .separator = { CI0_FIXED_0_PAD, UMI_LEN     } }, // UMI
            }
        };

        segconf.nonbio_con_snip_len = sizeof (segconf.nonbio_con_snip);
        container_prepare_snip ((ContainerP)&con, 0, 0, qSTRa(segconf.nonbio_con_snip));
    }

    // case: not 10xGenomics after all
    else 
        segconf.nonbio_type = NONBIO_NONE;
}

void fastq_10xGen_seg_initialize (VBlockFASTQP vb)
{
    ctx_consolidate_stats (VB, FASTQ_SQBITMAP, FASTQ_NONBIO, FASTQ_NONBIO_BC0, FASTQ_NONBIO_UMI, DID_EOL); 

    ctx_set_no_stons (VB, FASTQ_NONBIO_UMI, FASTQ_NONBIO_BC0, DID_EOL);
    
    ctx_set_ltype (VB, LT_BLOB, FASTQ_NONBIO_UMI, DID_EOL);

    CTX(FASTQ_NONBIO_UMI)->lcodec = CODEC_RANB;
    CTX(FASTQ_NONBIO_UMI)->lcodec_hard_coded = true;
} 

bool fastq_10xGen_seg_SEQ (VBlockFASTQP vb, STRp(seq))
{
    // if wrong length - seg verbatim
    if (seq_len != 28) 
        return false;

    seg_by_did (VB, seq, BARCODE_LEN, FASTQ_NONBIO_BC0, BARCODE_LEN);
    
    seg_add_to_local_fixed (VB, CTX(FASTQ_NONBIO_UMI), seq + BARCODE_LEN, UMI_LEN, LOOKUP_WITH_LENGTH, UMI_LEN);

    seg_by_did (VB, STRa(segconf.nonbio_con_snip), FASTQ_NONBIO, 0);

    return true;
}

