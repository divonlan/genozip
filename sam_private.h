// ------------------------------------------------------------------
//   sam_private.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SAM_PRIVATE_INCLUDED
#define SAM_PRIVATE_INCLUDED

#include "sam.h"
#include "vblock.h"

typedef struct {
    uint32_t qual_data_start, u2_data_start, bd_data_start, bi_data_start; // start within vb->txt_data
    uint32_t qual_data_len, u2_data_len, bd_data_len, bi_data_len;             // length within vb->txt_data
    uint32_t seq_len;        // actual sequence length determined from any or or of: CIGAR, SEQ, QUAL. If more than one contains the length, they must all agree
} ZipDataLineSAM;

typedef struct VBlockSAM {
    VBLOCK_COMMON_FIELDS
    SubfieldMapper qname_mapper;   // ZIP & PIZ
    Buffer optional_mapper_buf;    // PIZ: an array of type PizSubfieldMapper - one entry per entry in vb->contexts[SAM_OPTIONAL].mtf
    const char *last_cigar;        // ZIP/PIZ: last CIGAR
    uint32_t ref_consumed;         // ZIP/PIZ: how many bp of reference are consumed according to the last_cigar
    uint32_t ref_and_seq_consumed; // ZIP: how many bp in the last seq consumes both ref and seq, according to CIGAR
} VBlockSAM;

typedef VBlockSAM *VBlockSAMP;

#define DATA_LINE(i) ENT (ZipDataLineSAM, vb->lines, i)

#define CIGAR_DIGIT              1
#define CIGAR_CONSUMES_QUERY     2
#define CIGAR_CONSUMES_REFERENCE 4
#define CIGAR_INVALID            8
extern const uint8_t cigar_lookup[256];

extern void sam_analyze_cigar (const char *cigar, unsigned cigar_len, unsigned *seq_consumed, unsigned *ref_consumed, unsigned *seq_and_ref);

#endif
