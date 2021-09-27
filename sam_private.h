// ------------------------------------------------------------------
//   sam_private.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef SAM_PRIVATE_INCLUDED
#define SAM_PRIVATE_INCLUDED

#include "sam.h"
#include "vblock.h"
#include "reference.h"
#include "contigs.h"
#include "file.h"

#define DTYPE_QNAME        DTYPE_1
#define DTYPE_SAM_OPTIONAL DTYPE_2

#define BAM_MAGIC "BAM\1" // first 4 characters of a BAM file

#define IS_BAM (command==ZIP ? txt_file->data_type==DT_BAM \
                             : (z_file->data_type==DT_SAM && z_file->z_flags.txt_is_bin))

// as defined in https://samtools.github.io/hts-specs/SAMv1.pdf 1.4.2
#define SAM_FLAG_MULTI_SEGMENTS 0x0001
#define SAM_FLAG_IS_ALIGNED     0x0002
#define SAM_FLAG_UNMAPPED       0x0004
#define SAM_FLAG_NEXT_UNMAPPED  0x0008
#define SAM_FLAG_REV_COMP       0x0010
#define SAM_FLAG_NEXT_REV_COMP  0x0020
#define SAM_FLAG_IS_FIRST       0x0040
#define SAM_FLAG_IS_LAST        0x0080
#define SAM_FLAG_SECONDARY      0x0100
#define SAM_FLAG_FAILED_FILTERS 0x0200
#define SAM_FLAG_DUPLICATE      0x0400
#define SAM_FLAG_SUPPLEMENTARY  0x0800

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
    Buffer textual_opt;            // ZIP: Seg of BAM

    // data set by sam_analyze_cigar
    uint32_t ref_consumed;         // how many bp of reference are consumed according to the last_cigar
    uint32_t ref_and_seq_consumed; // how many bp in the last seq consumes both ref and seq, according to CIGAR
    uint32_t mismatch_bases;       // mismatch bases according to NM:i definition in https://samtools.github.io/hts-specs/SAMtags.pdf. This includes 'I' and 'D' bases.
    uint32_t soft_clip;            // PIZnumber of bases that were soft-clipped in this line
    
    // data set by sam_seg_analyze_MD
    bool md_verified;              // Seg: MD is verified to be calculatable from CIGAR, SEQ and Reference
    Buffer md_M_is_ref;            // Seg: bitmap of length ref_and_seq_consumed (corresponding to M, X and = CIGAR ops): 1 for a base matching the reference, 0 for non-matching, according to MD:Z

    Buffer bd_bi_line;             // ZIP: interlaced BD and BI data for one line
    
    // data used in genocat --show-sex
    WordIndex x_index, y_index, a_index;    // word index of the X, Y and chr1 chromosomes
    uint64_t x_bases, y_bases, a_bases;     // counters of number chromosome X, Y and chr1 bases
} VBlockSAM;

// fixed-field part of a BAM alignment, see https://samtools.github.io/hts-specs/SAMv1.pdf
typedef struct __attribute__ ((__packed__)) {
    uint32_t block_size;
    int32_t ref_id;
    int32_t pos;
    uint8_t l_read_name;
    uint8_t mapq;
    uint16_t bin;
    uint16_t n_cigar_op;
    uint16_t flag;
    uint32_t l_seq;
    int32_t next_ref_id;
    int32_t next_pos;  
    int32_t tlen;
    char read_name[]; // char[l_read_name]
} BAMAlignmentFixed;

typedef VBlockSAM *VBlockSAMP;

#define DATA_LINE(i) ENT (ZipDataLineSAM, vb->lines, i)

#define MAX_POS_SAM ((PosType)0x7fffffff)

#define CIGAR_DIGIT              1
#define CIGAR_CONSUMES_QUERY     2
#define CIGAR_CONSUMES_REFERENCE 4
extern const uint8_t cigar_lookup_sam[256];
extern const uint8_t cigar_lookup_bam[16];

// loading a Little Endian uint32_t from an unaligned buffer
#define GET_UINT16(p) (((uint8_t*)(p))[0] | (((uint8_t*)(p))[1] << 8))
#define GET_UINT32(p) (((uint8_t*)(p))[0] | (((uint8_t*)(p))[1] << 8) | (((uint8_t*)(p))[2] << 16) | (((uint8_t*)(p))[3] << 24))

// getting integers from the BAM data
#define NEXT_UINT8  *((const uint8_t *)next_field++)
#define NEXT_UINT16 GET_UINT16 (next_field); next_field += sizeof (uint16_t);
#define NEXT_UINT32 GET_UINT32 (next_field); next_field += sizeof (uint32_t);

#define dict_id_is_sam_qname_sf dict_id_is_type_1
#define dict_id_sam_qname_sf dict_id_type_1

extern void sam_seg_qname_field (VBlockSAM *vb, STRp(qname), unsigned add_additional_bytes);
extern void sam_seg_rname_rnext (VBlockP vb, DidIType did_i, STRp (chrom), unsigned add_bytes);
extern void sam_analyze_cigar (VBlockSAMP vb, STRp(cigar), unsigned *seq_consumed);
extern void sam_seg_tlen_field (VBlockSAM *vb, STRp(tlen), int64_t tlen_value, PosType pnext_pos_delta, int32_t cigar_seq_len);
extern void sam_seg_qual_field (VBlockSAM *vb, ZipDataLineSAM *dl, const char *qual, uint32_t qual_data_len, unsigned add_bytes);
extern void sam_seg_seq_field (VBlockSAM *vb, DidIType bitmap_did, STRp(seq), PosType pos, const char *cigar, unsigned recursion_level, uint32_t level_0_seq_len, const char *level_0_cigar, unsigned add_bytes);
extern const char *sam_seg_optional_all (VBlockSAM *vb, ZipDataLineSAM *dl, const char *next_field, int32_t len, bool *has_13, char separator, const char *after_field);
extern const char *bam_get_one_optional (VBlockSAM *vb, const char *next_field, const char **tag, char *type, const char **value, unsigned *value_len);
extern uint16_t bam_reg2bin (int32_t first_pos, int32_t last_pos);
extern void bam_seg_bin (VBlockSAM *vb, uint16_t bin, uint16_t flag, PosType this_pos);
extern void sam_seg_verify_rname_pos (VBlockP vb, const char *p_into_txt, PosType this_pos);

// ----------
// MD:Z stuff
// -----=----
extern void sam_md_analyze (VBlockSAMP vb, STRp(md), PosType pos, const char *cigar);
extern void sam_md_seg (VBlockSAM *vb,  ZipDataLineSAM *dl, STRp(md), unsigned add_bytes);

extern const char optional_sep_by_type[2][256];

static inline TranslatorId optional_field_translator (char type)
{
    switch (type) {
        case 'c' : return SAM2BAM_I8;
        case 'C' : return SAM2BAM_U8;
        case 's' : return SAM2BAM_LTEN_I16;
        case 'S' : return SAM2BAM_LTEN_U16;
        case 'i' : return SAM2BAM_LTEN_I32;
        case 'I' : return SAM2BAM_LTEN_U32;
        case 'f' : return IS_BAM ? 0 : SAM2BAM_FLOAT; // when reconstucting BAM->BAM we use SPECIAL_FLOAT rather than a translator
        default  : return 0;
    }
}

static inline char sam_seg_bam_type_to_sam_type (char type)
{
    return (type=='c' || type=='C' || type=='s' || type=='S' || type=='I') ? 'i' : type;
}

extern void sam_aux_zip_initialize (void);
extern DictId sam_seg_optional_field (VBlockSAM *vb, ZipDataLineSAM *dl, bool is_bam, const char *tag, char bam_type, const char *value, unsigned value_len);

#endif
