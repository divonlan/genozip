// ------------------------------------------------------------------
//   sam_private.h
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "sam.h"
#include "vblock.h"
#include "reference.h"
#include "contigs.h"
#include "file.h"
#include "seg.h"
#include "gencomp.h"
#include "digest.h"

#define DTYPE_QNAME        DTYPE_1
#define DTYPE_SAM_AUX DTYPE_2

#define BAM_MAGIC "BAM\1" // first 4 characters of a BAM file

// as defined in https://samtools.github.io/hts-specs/SAMv1.pdf section 1.4.2
#define SAM_FLAG_MULTI_SEGMENTS ((uint16_t)0x0001) // 1     0000 0000 0001
#define SAM_FLAG_IS_ALIGNED     ((uint16_t)0x0002) // 2     0000 0000 0010
#define SAM_FLAG_UNMAPPED       ((uint16_t)0x0004) // 4     0000 0000 0100
#define SAM_FLAG_NEXT_UNMAPPED  ((uint16_t)0x0008) // 8     0000 0000 1000
#define SAM_FLAG_REV_COMP       ((uint16_t)0x0010) // 16    0000 0001 0000
#define SAM_FLAG_NEXT_REV_COMP  ((uint16_t)0x0020) // 32    0000 0010 0000
#define SAM_FLAG_IS_FIRST       ((uint16_t)0x0040) // 64    0000 0100 0000
#define SAM_FLAG_IS_LAST        ((uint16_t)0x0080) // 128   0000 1000 0000
#define SAM_FLAG_SECONDARY      ((uint16_t)0x0100) // 256   0001 0000 0000
#define SAM_FLAG_FILTERED       ((uint16_t)0x0200) // 512   0010 0000 0000
#define SAM_FLAG_DUPLICATE      ((uint16_t)0x0400) // 1024  0100 0000 0000
#define SAM_FLAG_SUPPLEMENTARY  ((uint16_t)0x0800) // 2048  1000 0000 0000
#define SAM_MAX_FLAG            ((uint16_t)0x0FFF)

#pragma pack(1) 
typedef union SamFlags {
    struct SamFlagsBits {
        uint8_t multi_segments : 1;
        uint8_t is_aligned     : 1;
        uint8_t unmapped       : 1;
        uint8_t next_unmapped  : 1;
        uint8_t rev_comp       : 1;
        uint8_t next_rev_comp  : 1;
        uint8_t is_first       : 1;
        uint8_t is_last        : 1;
        uint8_t secondary      : 1;
        uint8_t filtered       : 1;
        uint8_t duplicate      : 1;
        uint8_t supplementary  : 1;
        uint8_t unused         : 4;    
    } bits;
    uint16_t value;
} SamFlags;

typedef union {
    struct {
        uint32_t op : 4;  // BamCigarOpType values
        uint32_t n  : 28;
    };
    uint32_t value;
} BamCigarOp;

typedef enum __attribute__ ((__packed__)) {
    BC_M=0, BC_I=1, BC_D=2, BC_N=3, BC_S=4, BC_H=5, BC_P=6, BC_E=7, BC_X=8, BC_NONE=15 
} BamCigarOpType;

extern const char cigar_op_to_char[16]; // BAM to SAM cigar op

#pragma pack()

#define IS_BAM_ZIP TXT_DT(DT_BAM)
#define IS_BAM_PIZ (Z_DT(DT_SAM) && z_file->z_flags.txt_is_bin) // note: this means z_file is BAM, NOT that the reconstruction is BAM! 
#define IS_BAM_RECON (IS_BAM_PIZ && vb->translation.trans_containers)
#define IS_BAM (command==ZIP ? IS_BAM_ZIP : IS_BAM_PIZ)

// Buddy & lookback parameters
#define BUDDY_MIN_RG_COUNT 3       // we only apply buddy to RG:Z if there are at least this many RGs
#define BUDDY_HASH_BITS 18    
#define COPY_BUDDY ((char)0x80)    // character entered in CIGAR and TLEN data to indicate a buddy copy (PART OF THE GENOZIP FILE FORMAT)

typedef enum __attribute__ ((__packed__)) {
    QUAL_NOT_MISSING, 
    QUAL_MISSING_STANDARD,         // BAM as in SAM spec (all 0xff) or SAM 
    QUAL_MISSING_PYSAM             // BAM - non-spec-compliant 0xff followed by 0x00 as observed in a file created by BSeeker2 (which uses pysam).
} QualMissingType;

typedef struct {
    TxtWord QUAL, U2, BD_BI[2];    // coordinates in txt_data 
    TxtWord QNAME, RG, CIGAR, MC;  // coordinates in txt_data for buddy segging (except CIGAR in BAM - points instead into vb->buddy_textual_cigars)
    TxtWord SEQ;                   // coordinates in txt_data. Note: len is actual sequence length in bases (not bytes) determined from any or or of: CIGAR, SEQ, QUAL. If more than one contains the length, they must all agree
    PosType POS, PNEXT;
    int64_t QUAL_score, TLEN;
    int64_t AS, YS;                // used for bowtie2
    int64_t MQ;                    // MQ is expeceted to be [0,255] (same as MAPQ) but we don't enforce
    uint8_t MAPQ;                  // 8 bit by BAM specification
    SamFlags FLAG;
    uint32_t ref_consumed;
} ZipDataLineSAM;

typedef struct VBlockSAM {
    VBLOCK_COMMON_FIELDS

    Buffer textual_seq;            // ZIP/PIZ: BAM: contains the textual SEQ (PIZ: used only in some cases)

    bool seq_missing;              // ZIP/PIZ: current line: SAM "*" BAM: l_seq=0
    bool cigar_missing;            // ZIP/PIZ: current line: SAM "*" BAM: n cigar op=0
    QualMissingType qual_missing;  // ZIP/PIZ: current line: SAM "*" BAM: l_seq=0 or QUAL all 0xff (which is defined in the spec as equivalent to SAM "*")

    // data set by sam_seg_analyze_MD
    bool md_verified;              // Seg: MD:Z is verified to be calculatable from CIGAR, SEQ and Reference
    Buffer md_M_is_ref;            // Seg: bitmap of length ref_and_seq_consumed (corresponding to M, X and = CIGAR ops): 1 for a base matching the reference, 0 for non-matching, according to MD:Z

    Buffer XG;                     // Seg: bsseeker2 XG:Z field with the underscores removed. revcomped if FLAG.revcomp.
                                   // PIZ: bsseeker2 XG:Z field with the underscores removed, as reconstructed when peeking during XM:Z special recon
    XgIncSType XG_inc_S;           // Seg: whether to include soft_clip[0]

    Buffer bd_bi_line;             // ZIP: interlaced BD and BI data for one line

    // CIGAR stuff
    rom last_cigar;                // ZIP: last CIGAR (PIZ: use vb->textual_cigar instead)
    Buffer textual_cigar;          // ZIP: Seg of BAM, PIZ: store CIGAR in sam_cigar_analyze
    Buffer binary_cigar;           // ZIP/PIZ: BAM-style CIGAR representation, generated in sam_cigar_analyze
    uint32_t ref_consumed;         // how many bp of reference are consumed according to the last_cigar
    uint32_t ref_and_seq_consumed; // how many bp in the last seq consumes both ref and seq, according to CIGAR
    uint32_t mismatch_bases;       // mismatch bases according to NM:i definition in https://samtools.github.io/hts-specs/SAMtags.pdf
    uint32_t soft_clip[2];         // [0]=soft clip at start [1]=soft clip at end of this line
    uint32_t hard_clip[2];         // [0]=hard clip at start [1]=hard clip at end of this line

    // SA stuff
    uint32_t plsg_i;               // PIZ: prim_vb: index of this VB in the plsg array  
    LineIType sa_grp_line_i;       // PIZ: the vb->line_i for which sa_grp was set
    Buffer sa_groups;              // ZIP/PIZ: an SA group is a group of alignments, including the primary aligngment
    Buffer sa_alns;                // ZIP/PIZ: array of {RNAME, STRAND, POS, CIGAR, NM, MAPQ} of the alignment
    Buffer sa_prim_cigars;         // ZIP/PIZ: textual primary CIGARs of SA groups (SAM & BAM)
    const struct SAGroupType *sa_grp; // ZIP/PIZ of SAM_COMP_DEPN: SA Group of this line (pointer into ZIP: z_file->sa_groups PIZ:prim_vb->sa_groups), NULL if none
    const struct SAAlnType *sa_aln;   // ZIP/PIZ of SAM_COMP_DEPN: Alignment withing SA Group of this line 
    uint32_t comp_qual_len;        // ZIP PRIM: compressed length of QUAL as it appears in in-memory SA Group
    uint32_t comp_cigars_len;      // ZIP PRIM: compressed length of CIGARS as it appears in in-memory SA Group
    int64_t NM;                    // ZIP: Value of NM:i and its length
    uint8_t NM_len;
    bool check_for_gc;             // ZIP: true if Seg should check for gencomp lines
    SAGroup first_grp_i;           // ZIP PRIM: the index of first group of this PRIM VB, in z_file->sa_groups 
    
    // data used in genocat --show-sex
    WordIndex x_index, y_index, a_index; // word index of the X, Y and chr1 chromosomes
    uint64_t x_bases, y_bases, a_bases;  // counters of number chromosome X, Y and chr1 bases

    // buddied Seg
    Buffer qname_hash;             // Seg: each entry i contains a line number for which the hash(qname)=i (or -1)
    Buffer buddy_textual_cigars;   // Seg of BAM (not SAM): an array of textual CIGARs referred to from DataLine->CIGAR

    // QUAL stuff
    bool qual_codec_no_longr;      // true if we can compress qual with CODEC_LONGR
} VBlockSAM;

#define VB_SAM ((VBlockSAMP)vb)

#define dict_id_is_aux_sf dict_id_is_type_2
#define dict_id_aux_sf dict_id_type_2

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

// rname,pos,strand,CIGAR,mapQ,NM
typedef enum { SA_RNAME, SA_POS, SA_STRAND, SA_CIGAR, SA_MAPQ, SA_NM, NUM_SA_ITEMS } SAFields __attribute__((unused)); // quick way to define constants

// SA Group stuff

// Alignments containing SA are removed from the SAM_COMP_MAIN vbs and moved to PRIM / DEPN components.
// When segging the PRIM component (containing primary SA lines) we place the data in the data structures below.
// When segging the DEPN component, we seg referring to the data structures.
//
// Seg of PRIM: vb->{sa_groups,sa_alns} is segging a PRIM vb and are merged into 
// z_file->{sa_groups,sa_alns,sa_qnames,sa_seq,sa_qual} during merge. 


#define ALN_CIGAR_BITS          48 
#define ALN_CIGAR_LEN_BITS_HI   7 
#define ALN_CIGAR_LEN_BITS_LO   13
#define ALN_CIGAR_LEN_BITS (ALN_CIGAR_LEN_BITS_HI + ALN_CIGAR_LEN_BITS_LO)
#define ALN_CIGAR_COMP_LEN_BITS 19
#define ALN_NM_BITS             23

#define CIGAR_SIG_LEN 12
typedef struct __attribute__ ((__packed__)) { 
    uint8_t bytes[CIGAR_SIG_LEN]; 
} CigarSignature;

typedef struct __attribute__ ((__packed__)) SAAlnType {
    // 24 bytes
    union {
        struct {
            uint64_t is_word  : 1;                       // Piz: true is CIGAR is in OPTION_SA_CIGAR.dict, false if it is in OPTION_SA_CIGAR.local    
            uint64_t index    : ALN_CIGAR_BITS;          // Piz: cigar_is_word ? WordIndex in OPTION_SA_CIGAR : index into OPTION_SA_CIGAR.local
            uint64_t len_hi   : ALN_CIGAR_LEN_BITS_HI;   // Piz, only if cigar_is_word=false
            uint32_t len_lo   : ALN_CIGAR_LEN_BITS_LO;   // Piz, only if cigar_is_word=false
            uint32_t comp_len : ALN_CIGAR_COMP_LEN_BITS; // Piz, only if cigar_is_word=false
        } piz;
        #define ALN_CIGAR_LEN(a) ((a)->cigar.piz.len_lo | ((a)->cigar.piz.len_hi << ALN_CIGAR_LEN_BITS_LO))

        CigarSignature signature; // Zip: first 96 bits of MD5 for CIGARs 13 characters or longer, or the CIGAR itself if shorter
    } cigar;

    WordIndex rname;              // Zip: word_index into RNAME == SAM header contig number.
                                  // Piz: word_index into OPTION_SA_RNAME
    uint32_t pos;                 // 31 bits per SAM spec
    uint32_t mapq           : 8;  // 8 bit by SAM specification
    uint32_t revcomp        : 1;  // 1 for - and 0 for +
    uint32_t nm             : ALN_NM_BITS; 
} SAAlnType;

#define ALN_NUM_ALNS_BITS 6 

// note: fields ordered to packed and word-aligned. These structures are NOT written to the genozip file.
typedef struct __attribute__ ((__packed__)) SAGroupType {
    // 40 bytes
    uint64_t qname           : 59;// index into: vb: txt_data ; z_file: zfile->sa_qnames
    uint64_t prim_set_buddy  : 1; // When reconstructing primary QNAME, sam_piz_set_sa_grp needs to reconstruct_set_buddy. 
    uint64_t first_grp_in_vb : 1; // This group is the first group in its PRIM vb
    uint64_t multi_segments  : 1; // FLAG.multi_segnments is set for all alignments of this group
    uint64_t is_first        : 1; // FLAG.is_last is set for all alignments of this group
    uint64_t is_last         : 1; // FLAG.is_last is set for all alignments of this group
    uint64_t seq;                 // index into: vb: txt_data ; z_file: zfile->sa_seq
    uint64_t qual;                // index into: vb: txt_data ; z_file: zfile->sa_qual
    uint64_t first_aln_i     : 47;// index into z_file->sa_alns (up to 128 trillion)
    uint64_t qname_len       : 8; // limited to 8 by BAM format
    uint64_t no_qual         : 2; // (QualMissingType) this SA has no quality, all dependends are expected to not have quality either
    uint64_t revcomp         : 1; // FLAG.revcomp of primary is set (this is same as revcomp for the first alignment)
    uint64_t num_alns        : ALN_NUM_ALNS_BITS; // number of alignments in this SA group (including primary alignment)
    uint32_t seq_len;             // number of bases in the sequences (not bytes). if !no_qual, this is also the length of qual.
    uint32_t qual_comp_len;       // compressed length of QUAL
} SAGroupType;

#define ZGRP_I(g) ((SAGroup)BNUM (z_file->sa_groups, (g))) // group pointer to grp_i
#define ZALN_I(a) ((uint64_t)BNUM (z_file->sa_alns, (a)))  // aln pointer to aln_i
#define GRP_QNAME(g) Bc (z_file->sa_qnames, (g)->qname)

#define MAXB(x) ((uint32_t)((1ULL<<(x))-1))
#define MAX_SA_NUM_ALNS       MAXB(ALN_NUM_ALNS_BITS)  // our limit
#define MAX_SA_POS            MAXB(31)                 // BAM limit
#define MAX_SA_NM             MAXB(ALN_NM_BITS)        // our limit
#define MAX_SA_MAPQ           MAXB(8)                  // BAM limit
#define MAX_SA_QNAME_LEN      MAXB(8)                  // BAM limit
#define MAX_SA_SEQ_LEN        MAXB(32)                 // BAM limit
#define MAX_SA_CIGAR_LEN      MAXB(ALN_CIGAR_LEN_BITS) // our limit - 20 bit, BAM limit - 16 bit for the CIGAR field (may be extended with CG)
#define MAX_SA_CIGAR_COMP_LEN MAXB(ALN_CIGAR_COMP_LEN_BITS)
#define MAX_SA_CIGAR_INDEX    MAXB(ALN_CIGAR_BITS)

#define DATA_LINE(i) B(ZipDataLineSAM, vb->lines, i)

#define MAX_POS_SAM ((PosType)0x7fffffff)

extern void sam_seg_QNAME (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qname), unsigned add_additional_bytes);
extern void sam_seg_FLAG (VBlockSAMP vb, ZipDataLineSAM *dl, unsigned add_bytes);
extern void sam_seg_RNAME_RNEXT (VBlockSAMP vb, DidIType did_i, STRp (chrom), unsigned add_bytes);
extern void sam_seg_MAPQ (VBlockSAMP vb, ZipDataLineSAM *dl, unsigned add_bytes);
extern void sam_seg_TLEN (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(tlen), int64_t tlen_value, bool is_rname_rnext_same);
extern void sam_seg_aux_all (VBlockSAMP vb, ZipDataLineSAM *dl, STRps(aux));
extern rom bam_show_line (VBlockSAMP vb, rom alignment, uint32_t remaining_txt_len);
extern void bam_get_one_optional (VBlockSAMP vb, STRp(aux), rom *tag, char *type, char *array_subtype, pSTRp(value), ValueType *numeric);
extern uint16_t bam_reg2bin (int32_t first_pos, int32_t last_pos);
extern void bam_seg_BIN (VBlockSAMP vb, ZipDataLineSAM *dl, uint16_t bin, bool is_bam);
extern void sam_seg_verify_RNAME_POS (VBlockSAMP vb, rom p_into_txt, PosType this_pos);
extern rom sam_seg_get_aux_str (VBlockSAMP vb, rom name, STRps (aux), uint32_t *value_len, bool is_bam);
extern uint32_t sam_seg_get_aux_int (VBlockSAMP vb, rom name, STRps (aux), int64_t *number, bool is_bam);
extern void sam_piz_set_sa_grp (VBlockSAMP vb);

static inline bool sam_seg_has_SA_Group(VBlockSAMP vb) { return (sam_is_depn_vb && vb->sa_aln) || sam_is_prim_vb; }

// ------------------
// POS / PNEXT stuff
// ------------------
extern PosType sam_seg_POS (VBlockSAMP vb, ZipDataLineSAM *dl, WordIndex prev_line_chrom, unsigned add_bytes);
extern void sam_seg_PNEXT (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(pnext_str)/* option 1 */, PosType pnext/* option 2 */, PosType prev_line_pos, unsigned add_bytes);

// ------------------
// CIGAR / MC:Z stuff
// ------------------
#define CIGAR_DIGIT              1
#define CIGAR_CONSUMES_QUERY     2
#define CIGAR_CONSUMES_REFERENCE 4
extern const uint8_t cigar_lookup_sam[256];
extern const uint8_t cigar_lookup_bam[16];

#define MAX_CIGAR_LEN_IN_DICT    7 // longer CIGARs are stored in local

#define HAVANA_DID_I(cigar_did_i) ((cigar_did_i)==SAM_CIGAR ? SAM_CIGAROP0 : (cigar_did_i)+3) // OPTION_*A_CIGAROP0 are 3 after OPTION_*A_CIGAR

extern void sam_cigar_analyze (VBlockSAMP vb, STRp(cigar), bool cigar_is_in_textual_cigar, uint32_t *seq_consumed);
extern void bam_seg_cigar_analyze (VBlockSAMP vb, uint32_t *seq_consumed);
extern void sam_cigar_binary_to_textual (VBlockSAMP vb, uint16_t n_cigar_op, const uint32_t *cigar, Buffer *textual_cigar);
extern void sam_cigar_seg_textual (VBlockSAMP vb, ZipDataLineSAM *dl, uint32_t last_cigar_len, STRp(seq_data), STRp(qual_data));
extern void sam_cigar_seg_binary (VBlockSAMP vb, ZipDataLineSAM *dl, uint32_t l_seq, rom cigar, uint32_t n_cigar_op);
extern void sam_cigar_seg_MC (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(mc), uint32_t add_bytes);
extern bool sam_cigar_reverse (char *dst, STRp(cigar));
extern bool sam_cigar_is_valid (STRp(cigar));
extern void sam_cigar_S_to_H (STRc(cigar));
extern void sam_cigar_H_to_S (STRc(cigar));
extern bool sam_seg_0A_cigar_cb (VBlockP vb, ContextP ctx, STRp (cigar), uint32_t repeat);
extern uint32_t sam_cigar_get_seq_len (STRp(cigar));

extern uint32_t sam_cigar_get_MC_ref_consumed (STRp(mc));
extern void sam_reconstruct_main_cigar_from_SA_Group (VBlockSAMP vb, bool substitute_S, bool reconstruct);
extern void sam_reconstruct_SA_cigar_from_SA_Group (VBlockSAMP vb, SAAlnType *a);

extern CigarSignature cigar_sign (STRp(cigar));
extern bool cigar_is_same_signature (CigarSignature sig1, CigarSignature sig2) ;
typedef struct { char s[CIGAR_SIG_LEN*2 + 1]; } DisCigarSig;
extern DisCigarSig cigar_display_signature (CigarSignature sig);

#define SA_CIGAR_DISPLAY_LEN 12
extern rom sam_piz_display_aln_cigar (const SAAlnType *a);

// ----------
// SEQ stuff
// ----------
extern void sam_seg_SEQ_initialize (VBlockSAMP vb);
extern void sam_seg_SEQ (VBlockSAMP vb, DidIType bitmap_did, STRp(seq), PosType pos, rom cigar, bool is_revcomp, uint32_t ref_consumed, uint32_t ref_and_seq_consumed, unsigned recursion_level, uint32_t level_0_seq_len, rom level_0_cigar, unsigned add_bytes);
extern void sam_zip_prim_ingest_vb_pack_seq (VBlockSAMP vb, SAGroupType *vb_grps, uint32_t vb_grps_len, Buffer *underlying_buf, Buffer *packed_seq_buf, bool is_bam_format);
extern bool sam_sa_native_to_actg (VBlockSAMP vb, BitArray *packed, uint64_t next_bit, STRp(seq), bool bam_format);

// -------------------
// BAM sequence format
// -------------------
extern const char bam_base_codes[16];
extern void bam_seq_to_sam (VBlockSAMP vb, bytes bam_seq, uint32_t seq_len, bool start_mid_byte, bool test_final_nibble, Buffer *out);
extern void sam_seq_to_bam (STRp (seq_sam), Buffer *seq_bam_buf);
extern void bam_seq_copy (VBlockSAMP vb, bytes bam_seq, uint32_t seq_len, bool start_mid_byte, Buffer *out);
extern uint32_t bam_split_aux (VBlockSAMP vb, rom aux, rom after_aux, rom *auxs, uint32_t *aux_lens);
extern rom bam_seq_display (bytes seq, uint32_t seq_len);

// ----------
// QUAL stuff
// ----------
extern void sam_seg_QUAL_initialize (VBlockSAMP vb);
extern void sam_seg_QUAL (VBlockSAMP vb, ZipDataLineSAM *dl, rom qual, uint32_t qual_data_len, unsigned add_bytes);
extern void sam_seg_ms_field (VBlockSAMP vb, ValueType ms, unsigned add_bytes);
extern rom bam_qual_display (bytes qual, uint32_t l_seq); 

#define SA_QUAL_DISPLAY_LEN 12
extern rom sam_display_qual_from_SA_Group (const SAGroupType *g);

// ----------
// MD:Z stuff
// ----------
extern void sam_seg_MD_Z_analyze (VBlockSAMP vb, STRp(md), PosType pos, rom cigar);
extern void sam_MD_Z_seg (VBlockSAMP vb,  ZipDataLineSAM *dl, STRp(md), unsigned add_bytes);

// ----------------
// BS-Seeker2 stuff
// ----------------
extern void sam_seg_XO_Z_field (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(XO), unsigned add_bytes);
extern void sam_seg_XM_Z_field (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(XM), unsigned add_bytes);
extern void sam_seg_XG_Z_field (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(XG), unsigned add_bytes);
extern void sam_seg_XG_Z_analyze (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(XG), PosType line_pos);

// ----------
// XA:Z stuff
// ----------
extern void sam_piz_XA_field_insert_lookback (VBlockP vb);

// -----------------------------------
// supplementary/secondary group stuff
// -----------------------------------
#define QNAME_HASH(qname,qname_len,is_last) ((libdeflate_crc32 (0, (qname), (qname_len)) & 0xfffffffe) | is_last) // note: adler32 is not good for qnames - since qnames in a SAM tend to be similar, many can end of up with the same hash and adler32 will use only a tiny bit of its 32b space
extern void sam_gc_zip_init_vb (VBlockP vb);
extern bool sam_seg_is_gc_line (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(alignment), STRps(aux), bool is_bam);
extern bool sam_seg_prim_add_sa_group_SA (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(sa), int64_t this_nm, bool is_bam);
extern void sam_seg_sa_group_stuff (VBlockSAMP vb, ZipDataLineSAM *dl, STRps(aux), STRp(textual_cigar), rom textual_seq, bool is_bam);
extern void sam_zip_gc_after_compute_main (VBlockSAMP vb);
extern void sam_sa_prim_initialize_ingest(void);
extern void sam_zip_prim_ingest_vb (VBlockSAMP vb);
extern void sam_seg_against_sa_group (VBlockSAMP vb, ContextP ctx, uint32_t add_bytes);
extern void sam_seg_against_sa_group_int (VBlockSAMP vb, ContextP ctx, int64_t parameter, uint32_t add_bytes);
extern void sam_seg_against_sa_group_bool (VBlockSAMP vb, ContextP ctx, bool parameter, uint32_t add_bytes);
extern SAGroupType *sam_piz_get_prim_vb_first_sa_grp (VBlockSAMP vb);
extern void sam_show_sa (void);
extern bool sam_is_SA_Groups_loaded (void);
extern uint32_t sam_piz_get_plsg_i (VBIType vb_i);

typedef struct { char s[1024]; } ShowAln;
extern void sam_show_sa_one_grp (SAGroup grp_i);
extern ShowAln sam_show_sa_one_aln (const SAGroupType *g, const SAAlnType *a);

#define SAM_PIZ_HAS_SA_GROUP (((VBlockSAMP)vb)->sa_grp && ((VBlockSAMP)vb)->sa_grp_line_i == vb->line_i + 1)

extern const char aux_sep_by_type[2][256];

// Seg
static inline TranslatorId aux_field_translator (char type)
{
    switch (type) {
        case 'c' : return SAM2BAM_I8;
        case 'C' : return SAM2BAM_U8;
        case 's' : return SAM2BAM_LTEN_I16;
        case 'S' : return SAM2BAM_LTEN_U16;
        case 'i' : return SAM2BAM_LTEN_I32;
        case 'I' : return SAM2BAM_LTEN_U32;
        case 'f' : return IS_BAM_ZIP ? 0 : SAM2BAM_FLOAT; // when reconstucting BAM->BAM we use SPECIAL_FLOAT rather than a translator
        default  : return 0;
    }
}

static inline char sam_seg_bam_type_to_sam_type (char type)
{
    return (type=='c' || type=='C' || type=='s' || type=='S' || type=='I') ? 'i' : type;
}

static inline char sam_seg_sam_type_to_bam_type (char type, int64_t n)
{
    LocalType test[6] = { LT_UINT8, LT_INT8, LT_UINT16, LT_INT16, LT_UINT32, LT_INT32 }; // preference to UINT
    
    if (type != 'i') return type; // all SAM types except 'i' are the same in BAM

    // i converts to one of 6: C,c,S,s,I,i
    for (int i=0 ; i < 6; i++)
        if (n >= lt_desc[test[i]].min_int && n <= lt_desc[test[i]].max_int)
            return lt_desc[test[i]].sam_type;
    
    return 0; // number out of range
}

extern DictId sam_seg_aux_field (VBlockSAMP vb, ZipDataLineSAM *dl, bool is_bam, rom tag, char bam_type, char bam_array_subtype, STRp(value), ValueType numeric);

typedef struct { char s[200]; } DisFlagsStr;
extern DisFlagsStr sam_dis_flags (SamFlags flags);

// -------------------
// SAM-private globals
// -------------------

extern const uint8_t aux_width[256];

extern char taxid_redirection_snip[100], xa_strand_pos_snip[100], XS_snip[30], XM_snip[30], MC_buddy_snip[30], MQ_buddy_snip[30], AS_buddy_snip[30], YS_buddy_snip[30], QUAL_buddy_snip[30], XA_lookback_snip[30], POS_buddy_snip[100], PNEXT_buddy_snip[100];
extern unsigned taxid_redirection_snip_len, xa_strand_pos_snip_len, XS_snip_len, XM_snip_len, MC_buddy_snip_len, MQ_buddy_snip_len, AS_buddy_snip_len, YS_buddy_snip_len, QUAL_buddy_snip_len, XA_lookback_snip_len, POS_buddy_snip_len, PNEXT_buddy_snip_len;
