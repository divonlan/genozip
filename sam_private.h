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
#define SAM_MAX_FLAG            0x0FFF

#pragma pack(1) 
typedef union {
    struct {
        uint8_t multi_segments : 1;
        uint8_t is_aligned     : 1;
        uint8_t unmapped       : 1;
        uint8_t next_unmapped  : 1;
        uint8_t rev_comp       : 1;
        uint8_t next_rev_comp  : 1;
        uint8_t is_first       : 1;
        uint8_t is_last        : 1;
        uint8_t secondary      : 1;
        uint8_t failed_filters : 1;
        uint8_t duplicate      : 1;
        uint8_t supplementary  : 1;
        uint8_t unused         : 4;    
    } bits;
    uint16_t value;
} SamFlags;
#pragma pack()

// Buddy & lookback parameters
#define BUDDY_MIN_RG_COUNT 3       // we only apply buddy to RG:Z if there are at least this many RGs
#define BUDDY_HASH_BITS 18    
#define COPY_BUDDY ((char)0x80)    // character entered in CIGAR and TLEN data to indicate a buddy copy (PART OF THE GENOZIP FILE FORMAT)

typedef struct {
    CtxWord QUAL, U2, BD_BI[2];    // coordinates in txt_data 
    CtxWord QNAME, RG, CIGAR, MC;  // coordinates in txt_data for buddy segging (except CIGAR in BAM - points instead into vb->buddy_textual_cigars)
    PosType POS, PNEXT, TLEN;
    SamFlags FLAG;
    uint32_t seq_len;              // actual sequence length determined from any or or of: CIGAR, SEQ, QUAL. If more than one contains the length, they must all agree
} ZipDataLineSAM;

typedef struct VBlockSAM {
    VBLOCK_COMMON_FIELDS
    const char *last_cigar;        // ZIP/PIZ: last CIGAR
    Buffer textual_cigar;          // ZIP: Seg of BAM, PIZ: store CIGAR in sam_cigar_analyze
    Buffer binary_cigar;           // PIZ: generate in sam_cigar_analyze, to reconstruct BAM
    Buffer textual_seq;            // ZIP: Seg of BAM
    Buffer textual_opt;            // ZIP: Seg of BAM

    // data set by sam_cigar_analyze
    uint32_t ref_consumed;         // how many bp of reference are consumed according to the last_cigar
    uint32_t ref_and_seq_consumed; // how many bp in the last seq consumes both ref and seq, according to CIGAR
    uint32_t mismatch_bases;       // mismatch bases according to NM:i definition in https://samtools.github.io/hts-specs/SAMtags.pdf. This includes 'I' and 'D' bases.
    uint32_t soft_clip;            // PIZ: number of bases that were soft-clipped in this line
    
    // data set by sam_seg_analyze_MD
    bool md_verified;              // Seg: MD is verified to be calculatable from CIGAR, SEQ and Reference
    Buffer md_M_is_ref;            // Seg: bitmap of length ref_and_seq_consumed (corresponding to M, X and = CIGAR ops): 1 for a base matching the reference, 0 for non-matching, according to MD:Z

    Buffer bd_bi_line;             // ZIP: interlaced BD and BI data for one line
    
    // data used in genocat --show-sex
    WordIndex x_index, y_index, a_index;    // word index of the X, Y and chr1 chromosomes
    uint64_t x_bases, y_bases, a_bases;     // counters of number chromosome X, Y and chr1 bases

    // buddied Seg
    Buffer qname_hash;             // Seg: each entry i contains a line number for which the hash(qname)=i (or -1)
    Buffer buddy_textual_cigars;   // Seg of BAM (not SAM): an array of textual CIGARs referred to from DataLine->CIGAR

} VBlockSAM;

#define VB_SAM ((VBlockSAMP)vb)

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

extern void sam_seg_QNAME (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qname), unsigned add_additional_bytes);
extern void sam_seg_FLAG (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(flag_str), unsigned add_bytes);
extern void sam_seg_RNAME_RNEXT (VBlockP vb, DidIType did_i, STRp (chrom), unsigned add_bytes);
extern PosType sam_seg_POS (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(pos_str)/* option 1 */, PosType pos/* option 2 */, WordIndex prev_line_chrom, PosType prev_line_pos, unsigned add_bytes);
extern void sam_seg_PNEXT (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(pnext_str)/* option 1 */, PosType pnext/* option 2 */, PosType prev_line_pos, unsigned add_bytes);
extern void sam_seg_TLEN (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(tlen), int64_t tlen_value, PosType pnext_pos_delta, int32_t cigar_seq_len);
extern void sam_seg_QUAL (VBlockSAM *vb, ZipDataLineSAM *dl, const char *qual, uint32_t qual_data_len, unsigned add_bytes);
extern void sam_seg_SEQ (VBlockSAM *vb, DidIType bitmap_did, STRp(seq), PosType pos, const char *cigar, uint32_t ref_consumed, uint32_t ref_and_seq_consumed, unsigned recursion_level, uint32_t level_0_seq_len, const char *level_0_cigar, unsigned add_bytes);
extern const char *sam_seg_optional_all (VBlockSAM *vb, ZipDataLineSAM *dl, const char *next_field, int32_t len, bool *has_13, char separator, const char *after_field);
extern const char *bam_get_one_optional (VBlockSAM *vb, const char *next_field, const char **tag, char *type, const char **value, unsigned *value_len);
extern uint16_t bam_reg2bin (int32_t first_pos, int32_t last_pos);
extern void bam_seg_BIN (VBlockSAM *vb, ZipDataLineSAM *dl, uint16_t bin, PosType this_pos);
extern void sam_seg_verify_RNAME_POS (VBlockP vb, const char *p_into_txt, PosType this_pos);

// ------------------
// CIGAR / MC:Z stuff
// -----=------------
extern void sam_cigar_analyze (VBlockSAMP vb, STRp(cigar), unsigned *seq_consumed);
extern void sam_cigar_binary_to_textual (VBlockSAM *vb, uint16_t n_cigar_op, const uint32_t *cigar, Buffer *textual_cigar);
extern void sam_cigar_seg_textual (VBlockSAM *vb, ZipDataLineSAM *dl, unsigned last_cigar_len, STRp(seq_data), STRp(qual_data));
extern void sam_cigar_seg_binary (VBlockSAM *vb, ZipDataLineSAM *dl, uint32_t l_seq, uint32_t n_cigar_op);
extern void sam_cigar_seg_MC (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(mc), unsigned add_bytes);
extern bool sam_cigar_reverse (char *dst, STRp(cigar));
extern bool sam_cigar_is_valid (STRp(cigar));

// ----------
// MD:Z stuff
// -----=----
extern void sam_md_analyze (VBlockSAMP vb, STRp(md), PosType pos, const char *cigar);
extern void sam_md_seg (VBlockSAM *vb,  ZipDataLineSAM *dl, STRp(md), unsigned add_bytes);

// ----------
// XA:Z stuff
// -----=----
extern void sam_piz_XA_field_insert_lookback (VBlockP vb);

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

extern DictId sam_seg_optional_field (VBlockSAM *vb, ZipDataLineSAM *dl, bool is_bam, const char *tag, char bam_type, const char *value, unsigned value_len);

extern char taxid_redirection_snip[100], xa_strand_pos_snip[100], XS_snip[30], XM_snip[30], MC_buddy_snip[30], xa_lookback_snip[30];
extern unsigned taxid_redirection_snip_len, xa_strand_pos_snip_len, XS_snip_len, XM_snip_len, MC_buddy_snip_len, xa_lookback_snip_len;

#endif
