// ------------------------------------------------------------------
//   sam_private.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// stuff needed by bamass shared between SAM and FASTQ

#pragma once
#include "genozip.h"

typedef enum { QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, AUX } SamFields __attribute__((unused)); // quick way to define constants

// fixed-field part of a BAM alignment, see https://samtools.github.io/hts-specs/SAMv1.pdf
typedef struct __attribute__((packed,aligned(1))) {
    // fixed-field
    uint32_t block_size;
    int32_t ref_id;
    PosType32 pos;
    uint8_t l_read_name;
    uint8_t mapq;
    uint16_t bin;
    uint16_t n_cigar_op;
    uint16_t flag;
    uint32_t l_seq;
    int32_t next_ref_id;
    PosType32 next_pos;  
    int32_t tlen;

    // variable-length fields (not included in sizeof(BAMAlignmentFixed))
    char read_name[/*l_read_name*/]; 
    // uint32_t cigar[n_cigar_op]
    // uint8_t seq[(l_seq+1)/2]
    // char qual[l_seq]
} BAMAlignmentFixed, *BAMAlignmentFixedP;

#define BAM_MAGIC  "BAM\1" // first 4 bytes of a BAM file
#define CRAM_MAGIC "CRAM"  // first 4 bytes of a CRAM file

// as defined in https://samtools.github.io/hts-specs/SAMv1.pdf section 1.4.2
#define SAM_FLAG_MULTI_SEG     ((uint16_t)0x0001) // 1     0000 0000 0001
#define SAM_FLAG_IS_ALIGNED    ((uint16_t)0x0002) // 2     0000 0000 0010
#define SAM_FLAG_UNMAPPED      ((uint16_t)0x0004) // 4     0000 0000 0100
#define SAM_FLAG_NEXT_UNMAPPED ((uint16_t)0x0008) // 8     0000 0000 1000
#define SAM_FLAG_REV_COMP      ((uint16_t)0x0010) // 16    0000 0001 0000
#define SAM_FLAG_NEXT_REV_COMP ((uint16_t)0x0020) // 32    0000 0010 0000
#define SAM_FLAG_IS_FIRST      ((uint16_t)0x0040) // 64    0000 0100 0000
#define SAM_FLAG_IS_LAST       ((uint16_t)0x0080) // 128   0000 1000 0000
#define SAM_FLAG_SECONDARY     ((uint16_t)0x0100) // 256   0001 0000 0000
#define SAM_FLAG_FILTERED      ((uint16_t)0x0200) // 512   0010 0000 0000
#define SAM_FLAG_DUPLICATE     ((uint16_t)0x0400) // 1024  0100 0000 0000
#define SAM_FLAG_SUPPLEMENTARY ((uint16_t)0x0800) // 2048  1000 0000 0000
#define SAM_MAX_FLAG           ((uint16_t)0x0FFF)

typedef union SamFlags {
    struct {
        uint8_t multi_segs    : 1;
        uint8_t is_aligned    : 1;
        uint8_t unmapped      : 1;
        uint8_t next_unmapped : 1;
        uint8_t rev_comp      : 1;
        uint8_t next_rev_comp : 1;
        uint8_t is_first      : 1;
        uint8_t is_last       : 1;
        uint8_t secondary     : 1;
        uint8_t filtered      : 1;
        uint8_t duplicate     : 1;
        uint8_t supplementary : 1;
        uint8_t unused        : 4;    
    };
    uint16_t value;
} SamFlags;

typedef union BamCigarOp {
    struct {
        uint32_t op : 4;  // BamCigarOpType values
        uint32_t n  : 28;
    };
    uint32_t value;
} BamCigarOp;

typedef packed_enum {
    BC_M=0, BC_I=1, BC_D=2, BC_N=3, BC_S=4, BC_H=5, BC_P=6, BC_E=7, BC_X=8, BC_NONE=15, BC_INVALID=255
} BamCigarOpType;

#define for_cigar(buf) for_buf (BamCigarOp, op, buf) switch (op->op) 

extern const char cigar_op_to_char[16]; // BAM to SAM cigar op

#define BAM_CIGAR_OP_NONE ((BamCigarOp){ .op=BC_NONE })

#define CIG(bin_cgar) (BamCigarOpP)STRb(bin_cgar)            

#define MAX_CIGAR_LEN_IN_DICT    7 // longer CIGARs are stored in local
extern void bam_seq_to_sam (VBlockP vb, bytes bam_seq, uint32_t seq_len, bool start_mid_byte, bool test_final_nibble, BufferP out, bool is_from_zip_cb);
extern bool sam_cigar_textual_to_binary (VBlockP vb, STRp(cigar), BufferP binary_cigar, rom buf_name);
extern void sam_cigar_binary_to_textual (VBlockP vb, ConstBamCigarOpP cigar, uint16_t n_cigar_op, bool reverse, BufferP textual_cigar);
extern void sam_seg_SEQ_initialize (VBlockP vb);
extern StrTextMegaLong dis_binary_cigar (VBlockP vb, const BamCigarOp *cigar, uint32_t cigar_len/*in ops*/, Buffer *working_buf); 
extern void sam_prepare_deep_cigar (VBlockP vb, ConstBamCigarOpP cigar, uint32_t cigar_len, bool reverse);
extern void sam_piz_produce_trivial_solo_huffmans (void);
extern uint64_t sam_deep_calc_hash_bits (void);

typedef enum { 
    SQUANK_BY_MAIN         = '0', // SAM:       from the seq_len+hard_clips implied by the primary CIGARs (used for SA_CIGAR/XA_CIGAR/OA_CIGAR and up to 15.0.68 also for MC:Z)
    SQUANK_BY_std_seq_len  = '1', // SAM/FASTQ: from the segconf.std_seq_len (up to 15.0.68: segconf.sam_seq_len which was the segconf-average seq len)
    SQUANK_BY_QNAME_length = '2'  // SAM/FASTQ: from the "length" component of QNAME (introduced 15.0.69) 
} SeqLenSource; // these values are part of the file format
#define SQUANK_NAMES { "BY_MAIN", "BY_std_seq_len", "BY_QNAME_length" }
extern bool squank_seg (VBlockP vb, ContextP ctx, STRp(cigar), uint32_t only_if_seq_len, SeqLenSource seq_len_source, uint32_t add_bytes);

// CIGAR snip opcodes - part of the file format
#define COPY_MATE_MC_Z             ((char)0x80)   // copy from mate's MC:Z
#define COPY_SAGGY_PRIM_SA_CIGAR   ((char)0x81)   // v14: copy from prim's SA_CIGAR
#define COPY_QNAME_LENGTH          ((char)0x82)   // v14: derive CIGAR from qname's length= component
#define SQUANK                     ((char)0x83)   // v14
#define COPY_QNAME_LENGTH_NO_CIGAR ((char)0x84)   // v15.0.26 - get seq_len from QNAME, and CIGAR is *
#define CIGAR_OP_NAMES { "COPY_MATE_MC_Z", "COPY_SAGGY_PRIM_SA_CIGAR", "COPY_QNAME_LENGTH", "SQUANK", "COPY_QNAME_LENGTH_NO_CIGAR" }