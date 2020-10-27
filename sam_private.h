// ------------------------------------------------------------------
//   sam_private.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SAM_PRIVATE_INCLUDED
#define SAM_PRIVATE_INCLUDED

#include "sam.h"
#include "vblock.h"

typedef struct {
    uint32_t qual_data_start, u2_data_start, bdbi_data_start[2]; // start within vb->txt_data
    uint32_t qual_data_len, u2_data_len; // length within vb->txt_data
    uint32_t seq_len;        // actual sequence length determined from any or or of: CIGAR, SEQ, QUAL. If more than one contains the length, they must all agree
} ZipDataLineSAM;

typedef struct VBlockSAM {
    VBLOCK_COMMON_FIELDS
    const char *last_cigar;        // ZIP/PIZ: last CIGAR
    Buffer textual_cigar;          // ZIP: Seg of BAM
    Buffer textual_seq;            // ZIP: Seg of BAM
    Buffer textual_opt;       // ZIP: Seg of BAM
    uint32_t ref_consumed;         // ZIP/PIZ: how many bp of reference are consumed according to the last_cigar
    uint32_t ref_and_seq_consumed; // ZIP: how many bp in the last seq consumes both ref and seq, according to CIGAR
    Buffer bd_bi_line;             // ZIP: interlaced BD and BI data for one line
} VBlockSAM;

typedef VBlockSAM *VBlockSAMP;

#define DATA_LINE(i) ENT (ZipDataLineSAM, vb->lines, i)

#define CIGAR_DIGIT              1
#define CIGAR_CONSUMES_QUERY     2
#define CIGAR_CONSUMES_REFERENCE 4
#define CIGAR_INVALID            8
extern const uint8_t cigar_lookup[256];

extern void sam_analyze_cigar (const char *cigar, unsigned cigar_len, unsigned *seq_consumed, unsigned *ref_consumed, unsigned *seq_and_ref);
extern void sam_seg_tlen_field (VBlockSAM *vb, const char *tlen, unsigned tlen_len, int64_t tlen_value, PosType pnext_pos_delta, int32_t cigar_seq_len);
extern void sam_seg_qual_field (VBlockSAM *vb, ZipDataLineSAM *dl, const char *qual, uint32_t qual_data_len, unsigned add_bytes);
extern void sam_seg_seq_field (VBlockSAM *vb, const char *seq, uint32_t seq_len, PosType pos, const char *cigar, unsigned recursion_level, uint32_t level_0_seq_len, const char *level_0_cigar, unsigned add_bytes);
extern const char *sam_seg_optional_all (VBlockSAM *vb, ZipDataLineSAM *dl, const char *next_field, int32_t len, bool *has_13, char separator, const char *after_field);
extern const char *bam_get_one_optional (VBlockSAM *vb, const char *next_field, const char **tag, char *type, const char **value, unsigned *value_len);

#endif
