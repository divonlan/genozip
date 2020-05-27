// ------------------------------------------------------------------
//   sam_private.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SAM_PRIVATE_INCLUDED
#define SAM_PRIVATE_INCLUDED

#include "sam.h"
#include "vblock.h"

typedef struct {
    uint32_t seq_data_start, qual_data_start, e2_data_start, u2_data_start, bd_data_start, bi_data_start; // start within vb->txt_data
    uint32_t seq_data_len, qual_data_len, e2_data_len, u2_data_len, bd_data_len, bi_data_len;             // length within vb->txt_data
    uint32_t seq_len;        // actual sequence length determined from any or or of: CIGAR, SEQ, QUAL. If more than one contains the length, they must all agree
} ZipDataLineSAM;

typedef struct VBlockSAM {
    VBLOCK_COMMON_FIELDS
    SubfieldMapper qname_mapper;         // ZIP & PIZ
    Buffer optional_mapper_buf;          // PIZ: an array of type PizSubfieldMapper - one entry per entry in vb->mtf_ctx[SAM_OPTIONAL].mtf
} VBlockSAM;

#define DATA_LINE(i) ENT (ZipDataLineSAM, vb->lines, i)

extern uint32_t sam_seq_len_from_cigar (const char *cigar, unsigned cigar_len);

#endif
