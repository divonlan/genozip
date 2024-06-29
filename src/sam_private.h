// ------------------------------------------------------------------
//   sam_private.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "vblock.h"
#include "contigs.h"
#include "file.h"
#include "seg.h"
#include "gencomp.h"
#include "digest.h"
#include "deep.h"
#include "piz.h"
#include "context.h"
#include "reconstruct.h"

#define DTYPE_QNAME   DTYPE_1
#define DTYPE_SAM_AUX DTYPE_2

#define BAM_MAGIC  "BAM\1" // first 4 characters of a BAM file
#define CRAM_MAGIC "CRAM"  // first 4 characters of a CRAM file

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

typedef union {
    struct {
        uint32_t op : 4;  // BamCigarOpType values
        uint32_t n  : 28;
    };
    uint32_t value;
} BamCigarOp;

typedef packed_enum {
    BC_M=0, BC_I=1, BC_D=2, BC_N=3, BC_S=4, BC_H=5, BC_P=6, BC_E=7, BC_X=8, BC_NONE=15, BC_INVALID=255
} BamCigarOpType;

#define BAM_CIGAR_OP_NONE ((BamCigarOp){ .op=BC_NONE })

extern const char cigar_op_to_char[16]; // BAM to SAM cigar op

typedef packed_enum {
    DEPN_CLIP_UNKNOWN, DEPN_CLIP_HARD, DEPN_CLIP_SOFT
} DepnClipping;

#define MAX_POS_SAM MAX_POS32 // according to SAM specification

#define MIN_TLEN -0x7fffffff // according to SAM specification
#define MAX_TLEN 0x7fffffff
typedef int32_t SamTlenType;

#define MIN_NM_i 0
#define MAX_NM_i 0x7fffffff  // at most the max length of a SAM contig
typedef int32_t SamNMType;

#define MIN_AS_i -0x7fffffff // note: in longranger and bowtie2, AS:i can be negative    
#define MAX_AS_i 0x7fffffff
typedef int32_t SamASType;

#define MAX_HI_NH 0x7fffffff

// Z fields which are expected to be the same as their prim - overlap . 
#define FIRST_SOLO_TAG MATED_BX
typedef enum         {                                                                                                                      SOLO_BX,     SOLO_RX,     SOLO_CB,     SOLO_CR,     SOLO_BC,/* solo-only-after-this*/ SOLO_QX, SOLO_CY, SOLO_QT, NUM_SOLO_TAGS } SoloTags;
                       /* only mated (not solo) tags first */
typedef enum         { MATED_RG,    MATED_PG,    MATED_PU,    MATED_LB,    MATED_OX,    MATED_MI,    MATED_ZA,    MATED_ZB,    MATED_YB,    MATED_BX,    MATED_RX,    MATED_CB,    MATED_CR,    MATED_BC,   NUM_MATED_Z_TAGS } MatedZFields;
#define MATED_Z_DIDs { OPTION_RG_Z, OPTION_PG_Z, OPTION_PU_Z, OPTION_LB_Z, OPTION_OX_Z, OPTION_MI_Z, OPTION_ZA_Z, OPTION_ZB_Z, OPTION_YB_Z, OPTION_BX_Z, OPTION_RX_Z, OPTION_CB_Z, OPTION_CR_Z, OPTION_BC_Z,                 }                
extern Did buddied_Z_dids[NUM_MATED_Z_TAGS];

#define SOLO_PROPS { /* in the order of SoloTags */ \
    { OPTION_BX_Z, false }, \
    { OPTION_RX_Z, true  }, \
    { OPTION_CB_Z, false }, \
    { OPTION_CR_Z, true  }, \
    { OPTION_BC_Z, false }, \
    { OPTION_QX_Z, false }, \
    { OPTION_CY_Z, false }, \
    { OPTION_QT_Z, false }, \
}

#define SOLO_CON_PEEK_ITEMS { /* in the order of SoloTags */ \
    { _OPTION_BX_Z, -1 }, \
    { _OPTION_RX_Z, -1 }, \
    { _OPTION_CB_Z, -1 }, \
    { _OPTION_CR_Z, -1 }, \
    { _OPTION_BC_Z, -1 }, \
    { _OPTION_QX_Z, -1 }, \
    { _OPTION_CY_Z, -1 }, \
    { _OPTION_QT_Z, -1 }, \
}

// Alignment used from SAG_BY_SOLO
typedef struct SoloAln {
    ZWord word[NUM_SOLO_TAGS]; // references into vb->solo_data (we can add more at any time, but not remove)
} SoloAln;

typedef struct {
    Did did_i;
    bool maybe_same_as_prev; // we might be able to save RAM due to the data in this tag often being identical to the previous tag in SoloTags
} SoloProp;

extern const SoloProp solo_props[NUM_SOLO_TAGS];

typedef struct {
    uint32_t hash;
    uint32_t count;
} QnameCount;

// number of alignments that are non-deepable for each of these reasons
typedef enum {      RSN_DEEPABLE, RSN_SECONDARY, RSN_SUPPLEMENTARY, RSN_NO_SEQ, RSN_MONOCHAR, RSN_CONSENSUS, NUM_DEEP_STATS_ZIP } DeepStatsZip; 
#define RSN_NAMES { "Deepable",   "Secondary",   "Supplementary",   "No_SEQ",   "Monochar",   "Consensus" }

// number of bytes consumed by QNAME/SEQ/QUAL in z_file->deep_ents
typedef enum {           QNAME_BYTES, SEQ_BYTES, QUAL_BYTES, NUM_DEEP_STATS_PIZ } DeepStatsPiz;
#define DEEP_ENT_NAMES { "QNAME",     "SEQ",     "QUAL" }

#if (NUM_DEEP_STATS_ZIP > NUM_DEEP_STATS) || (NUM_DEEP_STATS_PIZ > NUM_DEEP_STATS)
#error NUM_DEEP_STATS is too small
#endif

typedef struct {
    // we use this construction so that the same fields (eg UR) are in same memory location and can be 
    // accessed through either mated_z_fields or solo_z_fields
    union {
        struct { // mated 
            TxtWord mated_z_fields[NUM_MATED_Z_TAGS];
            TxtWord unused1[FIRST_SOLO_TAG + NUM_SOLO_TAGS - NUM_MATED_Z_TAGS];
        };
        struct { // solo 
            TxtWord unused2[FIRST_SOLO_TAG];
            TxtWord solo_z_fields[NUM_SOLO_TAGS];
        };
    };

    DeepHash deep_hash;            // Set for primary (i.e. not supplementary or secondary) alignments
    TxtWord QUAL, OQ, TQ, _2Y, GY, U2, BD_BI[2], SA, QNAME, MC, YB;// coordinates in txt_data 
    union {
    struct {                       // PacBio
    TxtWord dq, iq, sq;            // coordinates in txt_data
    int32_t np;
    };
    TxtWord t0;                    // Ultima - coordinates in txt_data
    int32_t nM;                    // STAR - used in paired-end STAR files
    };
    TxtWord rb, mb;                // NanoSeq - coordinates in txt_data
    TxtWord CIGAR;                 // SAM: coordinates in txt_data (always); BAM: coordinates in vb->line_textual_cigars
    TxtWord SEQ;                   // coordinates in txt_data. Note: len is actual sequence length in bases (not bytes) determined from any or or of: CIGAR, SEQ, QUAL. If more than one contains the length, they must all agree
    TxtWord BQ;
    int32_t QUAL_score;            // used by ms:i
    SamASType AS;
    SamASType YS;                  // used for bowtie2 and some other mappers
    WordIndex RNAME, RNEXT;
    PosType32 POS, PNEXT;
    PosType32 GP, MP;              // used for CellRanger-DNA         
    uint32_t seq_consumed;
    uint32_t ref_consumed;
    uint32_t hard_clip[2];        
    int32_t NH;                    // used by sam_seg_HI_i
    uint16_t PQ;
    SamNMType NM;                  // value of NM:i and its length
    SamFlags FLAG; 
    uint8_t NM_len;                
    uint8_t MAPQ, MQ;              // MAPQ is 8 bit by BAM specification, and MQ is mate MAPQ by SAMtags specification
    bool no_seq             : 1;   // SEQ is missing for this line
    bool is_deepable        : 1;   // True if this line is deepable: Not SEC/DUPP and SEQ exists and is not mono-char.
    bool is_consensus       : 1;   // This alignment is a consensus read
    bool no_qual            : 1;   // QUAL is missing this line
    bool dont_compress_QUAL : 1;   // don't compress QUAL in this line
    bool dont_compress_OQ   : 1;   // don't compress OQ in this line
    bool dont_compress_TQ   : 1;   // 
    bool dont_compress_QX   : 1;   // 
    bool dont_compress_CY   : 1;   // 
    bool dont_compress_QT   : 1;   // 
    bool dont_compress_GY   : 1;   // 
    bool dont_compress_2Y   : 1;   // 
} ZipDataLineSAM, *ZipDataLineSAMP;

// note: making VBlockSAM packed will ruin the math in sam_seg_buddied_i_fields
typedef struct VBlockSAM {
    VBLOCK_COMMON_FIELDS
    
    uint32_t arith_compress_bound_longest_seq_len; // PIZ

    // current line
    Buffer textual_seq;            // ZIP/PIZ: BAM: contains the textual SEQ (PIZ: used unless flags.no_textual_seq)
    rom textual_seq_str;           // ZIP/PIZ: BAM: points into textual_seq SAM: points into txt_data
    bool seq_missing;              // ZIP/PIZ: current line: SAM "*" BAM: l_seq=0
    bool seq_is_monochar;          // PIZ only
    bool cigar_missing;            // ZIP/PIZ: current line: SAM "*" BAM: n cigar op=0
    bool qual_missing;             // ZIP/PIZ: current line: SAM "*" BAM: l_seq=0 or QUAL all 0xff (of if segconf.pysam_qual: 0xff followed by 0s
    bool RNEXT_is_equal;           // RNEXT represents the same contig as RNAME (though the word index may differ in SAM, or RNEXT might be "=")
    bool line_not_deepable;        // PIZ only
    PizDeepSeq deep_seq;           // PIZ Deep: reconstruction instructions of SEQ, used only if cigar is perfect (eg 151M) and there is at most one mismatch vs reference

    uint32_t n_auxs;               // ZIP: AUX data of this line
    rom *auxs;                     
    uint32_t *aux_lens;

    ConstContainerP aux_con;       // AUX container being reconstructed

    // REF_INTERNAL and REF_EXT_STORE: the current length of a consecutive range of is_set
    WordIndex consec_is_set_chrom;
    PosType32 consec_is_set_pos;
    uint32_t consec_is_set_len;

    // Seg: 0-based index into AUX fields, -1 means field is not present in this line
    // When adding, also update sam_seg_idx_aux
    #define first_idx idx_NM_i
    int16_t idx_NM_i, idx_MD_Z, idx_SA_Z, idx_XG_Z, idx_NH_i, idx_HI_i, idx_IH_i,
            idx_X0_i, idx_X1_i, idx_XA_Z, idx_AS_i, 
            idx_CC_Z, idx_CP_i, idx_ms_i, idx_SM_i,
            idx_UB_Z, idx_BX_Z, idx_CB_Z, idx_GX_Z, idx_CR_Z, idx_CY_Z,
            idx_XO_Z, idx_YS_Z, idx_XB_A, idx_XM_Z, idx_XB_Z,
            idx_dq_Z, idx_iq_Z, idx_sq_Z, idx_ZA_Z, idx_ZB_Z,
            idx_pr_i, idx_qs_i, idx_ws_i, idx_ZM_B, idx_xq_i;
    #define has(f)   (vb->idx_##f  != -1)
    #define has_MD   (has(MD_Z) && segconf.has[OPTION_MD_Z])

    // data set by sam_seg_analyze_MD
    #define after_idx md_verified
    bool md_verified;              // Seg: MD:Z is verified to be calculatable from CIGAR, SEQ and Reference
    Buffer md_M_is_ref;            // Seg: bitmap of length ref_and_seq_consumed (corresponding to M, X and = CIGAR ops): 1 for a base matching the reference, 0 for non-matching, according to MD:Z

    Buffer deep_index;             // PIZ: entry per prim_line - uint32_t index into deep_ents
 
    #define first_mux mux_XS
    Multiplexer4 mux_XS;
    Multiplexer4 mux_PNEXT, mux_MAPQ;
    Multiplexer3 mux_POS;// ZIP: DEMUX_BY_MATE_PRIM multiplexers
    Multiplexer3 mux_GP, mux_MP;   // ZIP: channels: is_first, is_last, is_mated
    Multiplexer2 mux_FLAG, mux_MQ, mux_MC, mux_ms, mux_AS, mux_YS, mux_nM, // ZIP: DEMUX_BY_MATE or DEMUX_BY_BUDDY multiplexers
                 mux_mated_z_fields[NUM_MATED_Z_TAGS], mux_ultima_c, mux_dragen_sd, mux_YY, mux_XO, mux_PQ,
                 mux_sn, mux_rb, mux_mb;
    Multiplexer3 mux_NH;           // ZIP: DEMUX_BY_BUDDY_MAP
    Multiplexer7 mux_tp;           // ZIP: ULTIMA_tp (number of channels matches TP_NUM_BINS)

    // CIGAR stuff
    #define after_mux last_cigar
    rom last_cigar;                // ZIP: last CIGAR (PIZ: use vb->textual_cigar instead)
    Buffer textual_cigar;          // ZIP: Seg of BAM, PIZ: store CIGAR in sam_cigar_analyze
    Buffer binary_cigar;           // ZIP/PIZ: BAM-style CIGAR representation, generated in sam_cigar_analyze. binary_cigar.next is used by sam_seg_SEQ
    uint32_t ref_consumed;         // how many bp of reference are consumed according to the last_cigar
    uint32_t ref_and_seq_consumed; // how many bp in the last seq consumes both ref and seq, according to CIGAR
    int32_t mismatch_bases_by_SEQ; // ZIP/PIZ: mismatch bases according to SEQ - M/X/= bases that differ between SEQ and Reference
    int32_t mismatch_bases_by_MD;  // ZIP: mismatch bases according to MD - M/X/= bases that MD reports to be different between SEQ and Reference (in PIZ, this is stored in MD.last_value)
    uint32_t soft_clip[2];         // [0]=soft clip at start [1]=soft clip at end of this line
    uint32_t hard_clip[2];         // [0]=hard clip at start [1]=hard clip at end of this line
    uint32_t deletions;            // number of deleted bases according to CIGAR
    uint32_t insertions;           // number of inserted bases according to CIGAR
    uint128_t introns;             // number of introns ('N' CIGAR ops) - number of ops, not number of bases

    union {      
    // Deep (continued)
    Buffer deep_ents;              // PIZ: QNAME(plain), SEQ(packed) and QUAL(compressed) for each reconstructed list

    // Bisulfite conversion stuff
    Buffer unconverted_bitmap;     // ZIP: mismatches vs unconverted reference - used for MD:Z prediction
    };
    Buffer meth_call;              // ZIP/PIZ: prediction of methylation call in Bismark format: z/Z: unmethylated/methylated C in CpG ; x/X in CHG h/H in CHH ; u/U unknown context ; . not C. prm8[0] is bisulfite_strand ('C' or 'G') 
    
    // sag stuff
    uint32_t plsg_i;               // PIZ: prim_vb: index of this VB in the plsg array  
    LineIType sag_line_i;          // PIZ: the vb->line_i for which sag was set
    Buffer sag_grps;               // ZIP/PIZ: an SA group is a group of alignments, including the primary aligngment
    Buffer sag_alns;               // ZIP/PIZ: array of {RNAME, STRAND, POS, CIGAR, NM, MAPQ} of the alignment
    Buffer sa_prim_cigars;         // ZIP/PIZ: textual primary CIGARs of SAG_BY_SA (SAM & BAM)
    Buffer qname_count;            // ZIP: count the number of each qname_hash in this VB (if needed)

    const struct Sag *sag;         // ZIP/PIZ of SAM_COMP_DEPN: SA Group of this line (pointer into ZIP: z_file->sag_grps PIZ:prim_vb->sag_grps), NULL if none
    union {
    const struct SAAln *sa_aln;    // ZIP/PIZ of SAM_COMP_DEPN: SAG_BY_SA: Alignment withing SA Group of this line 
    const struct CCAln *cc_aln;    // ZIP/PIZ of SAM_COMP_DEPN: SAG_BY_CC: Prim alignment
    const struct SoloAln *solo_aln;// ZIP/PIZ of SAM_COMP_DEPN: SAG_BY_CC: Prim alignment
    };
    uint32_t comp_qual_len;        // ZIP PRIM: compressed length of QUAL as it appears in in-memory SA Group
    union {
    uint32_t comp_cigars_len;      // ZIP PRIM SAG_BY_SA: compressed length of CIGARS as it appears in in-memory SA Group
    uint32_t solo_data_len;        // ZIP PRIM SAG_BY_SOLO: size of solo_data
    };
    bool check_for_gc;             // ZIP: true if Seg should check for gencomp lines
    DepnClipping depn_clipping_type; // ZIP: For depn lines that have a clipping, what type of clipping do they have
    SAGroup first_grp_i;           // ZIP PRIM: the index of first group of this PRIM VB, in z_file->sag_grps 
    bool seg_found_prim_line;      // Seg MAIN: a prim line was found in the VB
    bool seg_found_depn_line;      // Seg MAIN: a prim line was found in the VB

    // data used in genocat --show-sex
    WordIndex x_index, y_index, a_index; // word index of the X, Y and chr1 chromosomes
    uint64_t x_bases, y_bases, a_bases;  // counters of number chromosome X, Y and chr1 bases

    // buddied Seg
    Buffer line_textual_cigars;    // Seg of BAM (not SAM): an array of textual CIGARs referred to from DataLine->CIGAR
    uint32_t saggy_near_count, mate_line_count, prim_far_count; // for stats
    LineIType mate_line_i;         // Seg/PIZ: the mate of this line. Goes into vb->buddy_line_i in Piz.
    LineIType saggy_line_i;        // Seg/PIZ: without gencomp: a previous line with the same sag as the current line.
    bool saggy_is_prim;            // Seg/PIZ: The line in saggy_line_i is 1. not secondary 2. not supplementary 3. has no hard clips

    // QUAL stuff
    bool has_qual;                 // Seg: This VB has at least one line with non-missing qual
    bool codec_requires_seq;       // Seg: one or more fields uses a codec that requires SEQ for reconstruction

    // stats
    uint32_t deep_stats[NUM_DEEP_STATS]; // ZIP/PIZ: stats collection regarding Deep - one entry for each in DeepStatsZip/DeepStatsPiz
} VBlockSAM;

#define VB_SAM ((VBlockSAMP)vb)

#define line_textual_cigars_used (segconf.has[OPTION_MC_Z] || /* sam_cigar_seg_MC_Z might be called  */ \
                                  IS_MAIN(vb))             /* sam_seg_SA_field_is_depn_from_prim might be called */
#define line_cigar(dl) Bc (*(IS_BAM_ZIP ? &vb->line_textual_cigars : &vb->txt_data), (dl)->CIGAR.index)

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
} BAMAlignmentFixed;

#define bam_seg_get_aux_B_template(array_subtype, element_type) \
typedef struct __attribute__((packed,aligned(1))) { char tag[2]; char B; char subtype; uint32_t array_len; element_type elem[]; } BamArray_##array_subtype;
bam_seg_get_aux_B_template (s, int16_t); // defines BamArray_s
bam_seg_get_aux_B_template (f, float);   // defines BamArray_f

typedef enum { QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, AUX } SamFields __attribute__((unused)); // quick way to define constants

typedef VBlockSAM *VBlockSAMP;

// rname,pos,strand,CIGAR,mapQ,NM
typedef enum { SA_RNAME, SA_POS, SA_STRAND, SA_CIGAR, SA_MAPQ, SA_NM, NUM_SA_ITEMS } SAFields __attribute__((unused)); // quick way to define constants

// SA Group stuff

// Alignments containing SA are removed from the SAM_COMP_MAIN vbs and moved to PRIM / DEPN components.
// When segging the PRIM component (containing primary SA lines) we place the data in the data structures below.
// When segging the DEPN component, we seg referring to the data structures.
//
// Seg of PRIM: vb->{sag_grps,sag_alns} is segging a PRIM vb and are merged into 
// z_file->{sag_grps,sag_alns,sag_qnames,sag_seq,sag_qual} during merge. 


#define ALN_CIGAR_BITS          48 
#define ALN_CIGAR_LEN_BITS_HI   7 
#define ALN_CIGAR_LEN_BITS_LO   13
#define ALN_CIGAR_LEN_BITS (ALN_CIGAR_LEN_BITS_HI + ALN_CIGAR_LEN_BITS_LO)
#define ALN_CIGAR_COMP_LEN_BITS 19
#define ALN_NM_BITS             23

#define CIGAR_SIG_LEN 12
typedef struct { 
    uint8_t bytes[CIGAR_SIG_LEN]; 
} CigarSignature;

// Alignment used from SAG_BY_SA
typedef struct SAAln {
    // 24 bytes
    union {
        struct __attribute__ ((packed,aligned(4))) {
            uint64_t is_word  : 1;                       // Piz: true is CIGAR is in OPTION_SA_CIGAR.dict, false if it is in OPTION_SA_CIGAR.local    
            uint64_t index    : ALN_CIGAR_BITS;          // Piz: cigar_is_word ? WordIndex in OPTION_SA_CIGAR : index into OPTION_SA_CIGAR.local
            uint64_t len_hi   : ALN_CIGAR_LEN_BITS_HI;   // Piz: only if cigar_is_word=false
            uint32_t len_lo   : ALN_CIGAR_LEN_BITS_LO;   // Piz: only if cigar_is_word=false
            uint32_t comp_len : ALN_CIGAR_COMP_LEN_BITS; // Piz: only if cigar_is_word=false
        } piz;
        #define ALN_CIGAR_LEN(a) ((a)->cigar.piz.len_lo | ((a)->cigar.piz.len_hi << ALN_CIGAR_LEN_BITS_LO))

        CigarSignature signature; // Zip: first 96 bits of MD5 for CIGARs 13 characters or longer, or the CIGAR itself if shorter
    } cigar;

    WordIndex rname;       // Zip: word_index into RNAME == SAM header contig number.
                           // Piz: word_index into OPTION_SA_RNAME
    PosType32 pos;         // 31 bits per SAM spec
    uint32_t mapq    : 8;  // 8 bit by SAM specification
    uint32_t revcomp : 1;  // 1 for - and 0 for +
    uint32_t nm      : ALN_NM_BITS; 
} SAAln;

// Alignment used from SAG_BY_CC
typedef struct CCAln {
    WordIndex rname;  // RNAME of prim alignment (alignment with NH>2 but no CC,CP) ; sometimes WORD_INDEX_NONE if rname is not in sam header / ref file
    PosType32 pos;   // POS of prim alignment
} CCAln;

// PIZ: history of cigar analysis - one item per line
typedef struct CigarAnalItem {
    uint32_t seq_len      : 31;  // equivalent to dl->SEQ.len in ZIP - set if ANY of SEQ, QUAL, CIGAR-implied-seq-consumed have it
    uint32_t qual_missing : 1;   // this alignment has no QUAL (or '*')
    uint32_t ref_consumed;
    uint32_t hard_clip[2];
} CigarAnalItem;

#define ALN_NUM_ALNS_BITS 12
#define GRP_SEQ_LEN_BITS  32
#define GRP_QUAL_BITS     48
#define GRP_AS_BITS       16

// note: fields ordered so to be packed and word-aligned. These structures are NOT written to the genozip file.
typedef struct Sag {
    // 40 bytes
    uint64_t qname           : 60;                // index into: vb: txt_data ; z_file: zfile->sag_qnames
    uint64_t first_grp_in_vb : 1;                 // This group is the first group in its PRIM vb
    uint64_t multi_segs      : 1;                 // FLAG.multi_segnments is set for all alignments of this group
    uint64_t is_first        : 1;                 // FLAG.is_last is set for all alignments of this group
    uint64_t is_last         : 1;                 // FLAG.is_last is set for all alignments of this group
    uint64_t seq;                                 // index into: vb: txt_data ; z_file: zfile->sag_seq
    uint64_t qual            : GRP_QUAL_BITS;     // index into: vb: txt_data ; z_file: zfile->sag_qual (256 TB)
    uint64_t as              : GRP_AS_BITS;       // AS:i ([0,65535] - capped at 0 and 65535)
    uint64_t first_aln_i     : 41;                // index into z_file->sag_alns (up to 2 trillion)
    uint64_t qname_len       : 8;                 // limited to 8 by BAM format
    uint64_t no_qual         : 2;                 // (QualMissingType) this SA has no quality, all dependends are expected to not have quality either
    uint64_t revcomp         : 1;                 // FLAG.revcomp of primary is set (in SA:Z-defined groups, this is same as revcomp for the first alignment)
    uint64_t num_alns        : ALN_NUM_ALNS_BITS; // number of alignments in this SA group (including primary alignment)
    uint64_t seq_len         : GRP_SEQ_LEN_BITS;  // number of bases in the sequences (not bytes). if !no_qual, this is also the length of qual.
    uint64_t qual_comp_len   : GRP_SEQ_LEN_BITS;  // compressed length of QUAL
} Sag; // Originally named "SA Group" representing the group of alignments in an SA:Z field, but now extended to other types of SAGs

#define ZGRP_I(g) ((SAGroup)BNUM (z_file->sag_grps, (g)))   // group pointer to grp_i
#define ZALN_I(a) (BNUM64 (z_file->sag_alns, (a)))  // aln pointer to aln_i
#define GRP_QNAME(g) Bc (z_file->sag_qnames, (g)->qname)

#define MAX_SA_NUM_ALNS       MAXB(ALN_NUM_ALNS_BITS)  // our limit
#define MAX_SA_POS            MAXB(31)                 // BAM limit
#define MAX_SA_NM             MAXB(ALN_NM_BITS)        // our limit
#define MAX_SA_MAPQ           MAXB(8)                  // BAM limit
#define MAX_SA_SEQ_LEN        MAXB(GRP_SEQ_LEN_BITS)   // our limit
#define MAX_SA_CIGAR_LEN      MAXB(ALN_CIGAR_LEN_BITS) // our limit - 20 bit, BAM limit - 16 bit for the CIGAR field (may be extended with CG)
#define MAX_SA_CIGAR_COMP_LEN MAXB(ALN_CIGAR_COMP_LEN_BITS)
#define MAX_SA_CIGAR_INDEX    MAXB64(ALN_CIGAR_BITS)
#define CAP_SA_AS(as)         MAX_(0, MIN_((as), (SamASType)MAXB(GRP_AS_BITS)))  // our limit - [0,65535] - capped at 0 and 65535

#define DATA_LINE(i) ((i) >= 0 ? B(ZipDataLineSAM, vb->lines, (i)) : NULL)

#define IS_MAIN(x) ((x)->comp_i == SAM_COMP_MAIN)
#define IS_PRIM(x) ((x)->comp_i == SAM_COMP_PRIM)
#define IS_DEPN(x) ((x)->comp_i == SAM_COMP_DEPN)

#define DICT_ID_ARRAY(dict_id) (DictId){ .id = { (dict_id).id[0], (dict_id).id[1], type,  '_','A','R','R' } } // DTYPE_2

extern void      sam_seg_QNAME (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(qname), unsigned add_additional_bytes);
extern WordIndex sam_seg_RNAME (VBlockSAMP vb, ZipDataLineSAMP dl, STRp (chrom), bool against_sa_group_ok, unsigned add_bytes);
extern WordIndex sam_seg_RNEXT (VBlockSAMP vb, ZipDataLineSAMP dl, STRp (chrom), unsigned add_bytes);
extern void      sam_seg_MAPQ  (VBlockSAMP vb, ZipDataLineSAMP dl, unsigned add_bytes);
extern void      sam_seg_TLEN  (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(tlen), SamTlenType tlen_value, bool is_rname_rnext_same);
extern void      bam_seg_BIN   (VBlockSAMP vb, ZipDataLineSAMP dl, uint16_t bin, bool is_bam);

// mismatch counters
extern void sam_seg_NM_i (VBlockSAMP vb, ZipDataLineSAMP dl, SamNMType NM, unsigned add_bytes);
extern void sam_seg_XM_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t xm, int16_t idx, unsigned add_bytes);

extern uint16_t bam_reg2bin (int32_t first_pos, int32_t last_pos);
extern void sam_seg_verify_RNAME (VBlockSAMP vb);
extern void sam_piz_set_sag (VBlockSAMP vb);
extern bool sam_seg_test_biopsy_line (VBlockP vb, STRp(line));
extern void sam_seg_float_as_snip (VBlockSAMP vb, ContextP ctx, STRp(sam_value), ValueType bam_value, unsigned add_bytes);

static inline bool sam_is_depn (SamFlags f) { return f.supplementary || f.secondary; }

static bool inline sam_line_is_depn (ZipDataLineSAMP dl) { return  sam_is_depn (dl->FLAG); }
static bool inline sam_line_is_prim (ZipDataLineSAMP dl) { return !sam_is_depn (dl->FLAG); }

#define sam_has_mate  (VB_SAM->mate_line_i  != NO_LINE)
#define sam_has_saggy (VB_SAM->saggy_line_i != NO_LINE)
#define sam_has_prim  (VB_SAM->saggy_line_i != NO_LINE && VB_SAM->saggy_is_prim) 

#define zip_has_real_prim (has_saggy && VB_SAM->saggy_is_prim)

#define piz_has_real_prim (piz_has_buddy && ((segconf.is_paired && sam_is_depn(last_flags)) || \
                                             (!segconf.is_paired && !sam_is_depn ((SamFlags){ .value = history64(SAM_FLAG, VB_SAM->buddy_line_i)}))))

// segconf stuff
extern void sam_segconf_set_by_MP (void);
extern void sam_segconf_finalize_optimizations (void);

// Header stuff
extern void sam_header_zip_inspect_SQ_lines_in_cram (rom cram_filename);

// BUDDY stuff
extern void sam_piz_set_buddy_v13 (VBlockP vb);
extern void sam_reconstruct_from_buddy_get_textual_snip (VBlockSAMP vb, ContextP ctx, BuddyType bt, pSTRp(snip));
extern rom buddy_type_name (BuddyType bt);

// FLAG stuff
extern StrTextLong sam_dis_FLAG (SamFlags f);
extern void sam_seg_FLAG (VBlockSAMP vb, ZipDataLineSAMP dl, unsigned add_bytes);
#define last_flags ((SamFlags){ .value = CTX(SAM_FLAG)->last_value.i })

// AUX stuff
extern void sam_seg_aux_all (VBlockSAMP vb, ZipDataLineSAMP dl);
extern rom bam_show_line (VBlockSAMP vb, rom alignment, uint32_t remaining_txt_len);
extern void bam_get_one_aux (VBlockSAMP vb, int16_t idx, rom *tag, char *type, char *array_subtype, pSTRp(value), ValueType *numeric);
extern void sam_seg_idx_aux (VBlockSAMP vb);
extern bool sam_piz_line_has_aux_field (VBlockSAMP vb, DictId dict_id);

extern uint32_t sam_seg_get_aux_int (VBlockSAMP vb, int16_t idx, int32_t *number, bool is_bam, int32_t min_value, int32_t max_value, FailType soft_fail);
#define sam_seg_get_aux_int_(vb,tag) ({ int32_t value; sam_seg_get_aux_int (vb, vb->idx_##tag, &value, IS_BAM_ZIP, -0x80000000, 0x7fffffff, HARD_FAIL); value; })

extern uint32_t sam_seg_get_aux_float (VBlockSAMP vb, int16_t idx, double *number, bool is_bam, FailType soft_fail);

extern void sam_seg_get_aux_Z (VBlockSAMP vb, int16_t idx, pSTRp (snip), bool is_bam);
extern char sam_seg_get_aux_A (VBlockSAMP vb, int16_t idx, bool is_bam);

extern bool sam_seg_peek_int_field (VBlockSAMP vb, Did did_i, int16_t idx, int32_t min_value, int32_t max_value, bool set_last_value, int32_t *value);

extern uint32_t bam_split_aux (VBlockSAMP vb, rom alignment, rom aux, rom after_aux, rom *auxs, uint32_t *aux_lens);

typedef void (*SegBuddiedCallback)(VBlockSAMP, ContextP, STRp(value), unsigned add_bytes);
extern void sam_seg_buddied_Z_fields (VBlockSAMP vb, ZipDataLineSAMP dl, MatedZFields f, STRp(value), SegBuddiedCallback seg_cb, unsigned add_bytes);
extern void sam_seg_buddied_i_fields (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, int64_t my_value, int32_t *mate_value, MultiplexerP mux, STRp(copy_snip), unsigned add_bytes);

#define STRauxZ(name,is_bam) (vb->auxs[vb->idx_##name]+((is_bam) ? 3 : 5)), (vb->aux_lens[vb->idx_##name]-((is_bam) ? 4 : 5))

// POS / PNEXT stuff
extern PosType32 sam_seg_POS (VBlockSAMP vb, ZipDataLineSAMP dl, WordIndex prev_line_chrom, unsigned add_bytes);
extern void sam_seg_PNEXT (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(pnext_str)/* option 1 */, PosType32 pnext/* option 2 */, unsigned add_bytes);

// CIGAR / MC:Z stuff
#define CIGAR_DIGIT              1
#define CIGAR_CONSUMES_QUERY     2
#define CIGAR_CONSUMES_REFERENCE 4
extern const uint8_t cigar_lookup_sam[256];
extern const uint8_t cigar_lookup_bam[16];

#define MAX_CIGAR_LEN_IN_DICT    7 // longer CIGARs are stored in local

#define HAVANA_DID_I(cigar_did_i) ((cigar_did_i)==SAM_CIGAR ? SAM_CIGAROP0 : (cigar_did_i)+3) // OPTION_*A_CIGAROP0 are 3 after OPTION_*A_CIGAR

extern void sam_seg_cigar_initialize (VBlockSAMP vb);
extern void sam_cigar_analyze (VBlockSAMP vb, STRp(cigar), bool cigar_is_in_textual_cigar, uint32_t *seq_consumed);
extern void bam_seg_cigar_analyze (VBlockSAMP vb, ZipDataLineSAMP dl, uint32_t *seq_consumed);
extern void sam_cigar_binary_to_textual (VBlockSAMP vb, uint16_t n_cigar_op, const BamCigarOp *cigar, BufferP textual_cigar);
extern bool sam_cigar_textual_to_binary (VBlockSAMP vb, STRp(cigar), BufferP binary_cigar);
extern bool sam_is_cigar (STRp(cigar), bool allow_empty);
extern void sam_seg_CIGAR (VBlockSAMP vb, ZipDataLineSAMP dl, uint32_t last_cigar_len, STRp(seq_data), STRp(qual_data), uint32_t add_bytes);
extern void sam_cigar_seg_binary (VBlockSAMP vb, ZipDataLineSAMP dl, uint32_t l_seq, BamCigarOp *cigar, uint32_t n_cigar_op);
extern void sam_cigar_seg_MC_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(mc), uint32_t add_bytes);
extern bool sam_cigar_reverse (char *dst, STRp(cigar));
extern bool sam_cigar_is_valid (STRp(cigar));
extern void sam_seg_other_CIGAR (VBlockSAMP vb, ContextP ctx, STRp (cigar), bool squanking_allowed, unsigned add_bytes);
extern rom display_binary_cigar (VBlockSAMP vb);

extern bool sam_seg_0A_rname_cb (VBlockP vb, ContextP ctx, STRp(oa_rname), uint32_t repeat);
extern bool sam_seg_0A_pos_cb (VBlockP vb, ContextP ctx, STRp(oa_pos), uint32_t repeat);
extern bool sam_seg_0A_cigar_cb (VBlockP vb, ContextP ctx, STRp (cigar), uint32_t repeat);

typedef struct { char s[10000]; } CigarStr;
extern CigarStr dis_binary_cigar (VBlockSAMP vb, const BamCigarOp *cigar, uint32_t cigar_len/*in ops*/, Buffer *working_buf); 

typedef struct { char *left, *right; bool is_binary; } HtoS, StoH;
extern HtoS sam_cigar_H_to_S (VBlockSAMP vb, STRc(cigar), bool is_binary);
extern StoH sam_cigar_S_to_H (VBlockSAMP vb, STRc(cigar), bool is_binary);
static inline void sam_cigar_restore_H (HtoS htos) { if (htos.is_binary) { if (htos.left) ((BamCigarOp*)htos.left)->op = BC_H; if (htos.right) ((BamCigarOp*)htos.right)->op = BC_H; } \
                                                     else                { if (htos.left)  *htos.left  = 'H';                  if (htos.right) *htos.right = 'H'; } }
static inline void sam_cigar_restore_S (StoH stoh) { if (stoh.is_binary) { if (stoh.left) ((BamCigarOp*)stoh.left)->op = BC_S; if (stoh.right) ((BamCigarOp*)stoh.right)->op = BC_S; } \
                                                     else                { if (stoh.left)  *stoh.left  = 'S';  if (stoh.right) *stoh.right = 'S'; } }
extern bool sam_cigar_has_H (STRp(cigar)); // textual

extern uint32_t sam_cigar_get_MC_ref_consumed (STRp(mc));
extern void sam_reconstruct_main_cigar_from_sag (VBlockSAMP vb, bool do_htos, ReconType reconstruct);
extern void sam_reconstruct_SA_cigar_from_SA_Group (VBlockSAMP vb, SAAln *a);

extern CigarSignature cigar_sign (STRp(cigar));
extern bool cigar_is_same_signature (CigarSignature sig1, CigarSignature sig2) ;
extern StrText cigar_display_signature (CigarSignature sig);

#define SA_CIGAR_DISPLAY_LEN 12
extern rom sam_piz_display_aln_cigar (const SAAln *a);

// SEQ stuff
extern void sam_seg_SEQ_initialize (VBlockSAMP vb);
extern void sam_seg_SEQ (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(seq), unsigned add_bytes);
extern void sam_zip_prim_ingest_vb_pack_seq (VBlockSAMP vb, Sag *vb_grps, uint32_t vb_grps_len, BufferP underlying_buf, BufferP packed_seq_buf, bool is_bam_format);
extern bool sam_seq_pack (VBlockSAMP vb, Bits *packed, uint64_t next_bit, STRp(seq), bool bam_format, bool revcomp, FailType soft_fail);
extern rom sam_seg_analyze_set_one_ref_base (VBlockSAMP vb, bool is_depn, PosType32 pos, char base, uint32_t ref_consumed, RangeP *range_p, RefLock *lock);
extern void sam_zip_report_monochar_inserts (void);

// BAM sequence format
extern const char bam_base_codes[16];
extern void bam_seq_to_sam (VBlockSAMP vb, bytes bam_seq, uint32_t seq_len, bool start_mid_byte, bool test_final_nibble, BufferP out);
extern void sam_seq_to_bam (STRp (seq_sam), Buffer *seq_bam_buf);
extern rom bam_seq_display (bytes seq, uint32_t seq_len);
extern uint32_t sam_seq_copy (char *dst, rom src, uint32_t src_start_base, uint32_t n_bases, bool revcomp, bool is_bam_format);

// QUAL stuff
#define QUAL_ZIP_CALLBACK(tag, f, may_be_revcomped)             \
COMPRESSOR_CALLBACK (sam_zip_##tag)                             \
{                                                               \
    ZipDataLineSAMP dl = DATA_LINE (vb_line_i);                 \
    *line_data_len = dl->dont_compress_##tag ? 0 : MIN_(maximum_size, dl->f.len); /* note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec */ \
    if (!line_data || ! *line_data_len) return; /* no data, or only lengths were requested */   \
    *line_data = Btxt (dl->f.index);                             \
    if (is_rev) *is_rev = may_be_revcomped ? dl->FLAG.rev_comp : false;\
}                               

extern void sam_seg_QUAL_initialize (VBlockSAMP vb);
extern void sam_seg_QUAL (VBlockSAMP vb, ZipDataLineSAMP dl, rom qual, uint32_t qual_data_len, unsigned add_bytes);
extern rom bam_qual_display (bytes qual, uint32_t l_seq); 
extern void sam_seg_other_qual (VBlockSAMP vb, ZipDataLineSAMP dl, TxtWord *dl_word, Did did_i, STRp(qual), bool len_is_seq_len, unsigned add_bytes);

// --deep stuff
// Stuff that happens during SAM seg
extern void sam_deep_zip_finalize (void);
extern void sam_deep_set_QNAME_hash (VBlockSAMP vb, ZipDataLineSAMP dl, QType q, STRp(qname));
extern void sam_deep_set_SEQ_hash (VBlockSAMP vb,ZipDataLineSAMP dl, STRp(textual_seq));
extern void sam_deep_set_QUAL_hash (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(qual));
extern void sam_piz_deep_initialize (void);
extern void sam_piz_deep_init_vb (VBlockSAMP vb, ConstSectionHeaderVbHeaderP header);
extern void sam_piz_deep_grab_deep_ents (VBlockSAMP vb);

#define SA_QUAL_DISPLAY_LEN 20
extern rom sam_display_qual_from_SA_Group (const Sag *g);

// MD:Z stuff
extern void sam_seg_MD_Z_analyze (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(md), PosType32 pos);
extern void sam_seg_MD_Z (VBlockSAMP vb,  ZipDataLineSAMP dl, STRp(md), unsigned add_bytes);
extern void sam_MD_Z_verify_due_to_seq (VBlockSAMP vb, STRp(seq), PosType32 pos, BitsP sqbitmap, uint64_t sqbitmap_start);

// SA:Z stuff
extern void sam_seg_SA_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(sa), unsigned add_bytes);
extern bool sam_seg_0A_mapq_cb (VBlockP vb, ContextP ctx, STRp (mapq_str), uint32_t repeat);
extern bool sam_seg_SA_get_prim_item (VBlockSAMP vb, int sa_item, pSTRp(out));
extern void sam_piz_SA_get_prim_item (VBlockSAMP vb, int sa_item, pSTRp(out));
extern bool sam_seg_is_item_predicted_by_prim_SA (VBlockSAMP vb, int sa_item_i, int64_t value);

// MM:Z stuff
extern void sam_MM_zip_initialize (void);

// BWA stuff
extern void sam_seg_BWA_X1_i (VBlockSAMP vb, int64_t X1, unsigned add_bytes);
extern void sam_seg_BWA_XT_A (VBlockSAMP vb, char XT, unsigned add_bytes);
extern void sam_seg_BWA_XC_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t XC, unsigned add_bytes);
extern void sam_seg_BWA_XS_i (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, ValueType XS, unsigned add_bytes);
extern void sam_seg_BWA_XA_Z (VBlockSAMP vb, STRp(xa), unsigned add_bytes);
extern void sam_piz_XA_field_insert_lookback_v13 (VBlockP vb);
extern void sam_seg_correct_for_semcol_in_contig (uint32_t *n_repeats, rom *repeats, uint32_t *repeat_lens);

// bowtie2 stuff
extern void sam_seg_bowtie2_YS_i (VBlockSAMP vb, ZipDataLineSAMP dl, ValueType YS, unsigned add_bytes);

// minimap2 stuff
extern void sam_seg_s1_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t s1, unsigned add_bytes);
extern void sam_seg_cm_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t cm, unsigned add_bytes);

// PacBio stuff
extern void sam_pacbio_seg_initialize (VBlockSAMP vb);
extern void sam_seg_pacbio_xq (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, TxtWord *dl_word, STRp(value), unsigned add_bytes);
extern void sam_seg_pacbio_np (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t np, unsigned add_bytes);
extern void sam_seg_pacbio_qs (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t qs, unsigned add_bytes);
extern void sam_seg_pacbio_qe (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t qe, unsigned add_bytes);
extern void sam_seg_pacbio_we (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t we, unsigned add_bytes);
extern void sam_seg_pacbio_zm (VBlockSAMP vb, int64_t zm, unsigned add_bytes);
extern bool sam_seg_pacbio_qual (VBlockSAMP vb, STRp(qual), unsigned add_bytes);
extern void sam_seg_pacbio_sn (VBlockSAMP vb, ZipDataLineSAMP dl, rom sn, int sn_len);
extern void sam_recon_pacbio_qual (VBlockSAMP vb, ContextP ctx, bool reconstruct);
extern void sam_recon_skip_pacbio_qual (VBlockSAMP vb);

// Ultima stuff
extern void sam_ultima_zip_initialize (void);
extern void sam_ultima_seg_initialize (VBlockSAMP vb);
extern void sam_ultima_finalize_segconf (VBlockSAMP vb);
extern void sam_seg_ULTIMA_tp (VBlockSAMP vb, ContextP ctx, void *cb_param, void *tp_, uint32_t tp_len);
extern void sam_seg_ultima_bi (VBlockSAMP vb, STRp(bi_str), unsigned add_bytes);
extern void sam_seg_ultima_XV (VBlockSAMP vb, STRp(xv), unsigned add_bytes);
extern void sam_seg_ultima_XW (VBlockSAMP vb, STRp(xw), unsigned add_bytes);
extern void sam_seg_ultima_t0 (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(t0), unsigned add_bytes);
extern void sam_seg_ultima_MI (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(mi), unsigned add_bytes);
extern void sam_seg_ultima_a3 (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t a3, unsigned add_bytes);
extern void sam_seg_ultima_pt (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t pt, unsigned add_bytes);

// STAR stuff
extern void sam_star_zip_initialize (void);
extern void sam_star_seg_initialize (VBlockSAMP vb);
extern void sam_seg_STAR_nM (VBlockSAMP vb, ZipDataLineSAMP dl, SamNMType nm, unsigned add_bytes);
extern void sam_seg_STAR_jI (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(raw), bool is_bam);

// Dragen stuff
extern void sam_dragen_seg_initialize (VBlockSAMP vb);
extern void sam_dragen_seg_sd_f (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(sd), ValueType numeric, unsigned add_bytes);

// hisat2 stuff
extern void sam_seg_HISAT2_Zs_Z (VBlockSAMP vb, STRp(zs), unsigned add_bytes);

// tmap (IonXpress) stuff
extern void sam_seg_TMAP_XM_i (VBlockSAMP vb, ValueType XM, unsigned add_bytes);

extern void sam_seg_init_bisulfite (VBlockSAMP vb, ZipDataLineSAMP dl);

// Bismark stuff
#define bisulfite_strand meth_call.prm8[0]
extern void sam_seg_bismark_XM_Z_analyze (VBlockSAMP vb, ZipDataLineSAMP dl);
extern void sam_seg_bismark_XM_Z (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, int special_code, STRp(xm), unsigned add_bytes);
extern void sam_seg_bismark_XG_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(xg), unsigned add_bytes);
extern void sam_bismark_zip_update_meth_call (VBlockSAMP vb, RangeP range, uint32_t range_len, int32_t idx, bool methylated, uint32_t M_i, uint32_t Mseg_len, uint32_t i);
extern void sam_bismark_piz_update_meth_call (VBlockSAMP vb, bytes ref, int32_t idx, uint32_t seq_i, uint32_t M_i, uint32_t Mseg_len, const BamCigarOp *cigar, char bisulfite, bool methylated);

// BS-Seeker2 stuff
extern void sam_seg_bsseeker2_XO_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(XO), unsigned add_bytes);
extern void sam_seg_bsseeker2_XM_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(XM), unsigned add_bytes);
extern void sam_seg_bsseeker2_XG_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(XG), unsigned add_bytes);
extern void sam_seg_bsseeker2_XG_Z_analyze (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(XG), PosType32 line_pos);

// BSBolt stuff
extern void sam_seg_bsbolt_XB_Z_analyze (VBlockSAMP vb, ZipDataLineSAMP dl);
extern void sam_seg_bsbolt_XB (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(XB), unsigned add_bytes);
extern void sam_seg_bsbolt_YS_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(YS), unsigned add_bytes);

// Gem3 mapper stuff
extern void sam_seg_gem3_XB_A (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(xb), unsigned add_bytes);
extern bool sam_seg_gem3_XA_strand_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep);

// BLASR stuff
extern void sam_seg_blasr_FI_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t fi, unsigned add_bytes);

// scRNA-seq stuff (STARsolo and cellranger)
extern void sam_segconf_retag_UBURUY (void);
extern void sam_10xGen_seg_initialize (VBlockSAMP vb);
extern void sam_seg_TX_AN_Z (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, STRp(value), unsigned add_bytes);

extern void sam_seg_CR_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(cr), unsigned add_bytes);
extern void sam_seg_CB_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(cb), unsigned add_bytes);
extern void sam_seg_RX_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(ur), unsigned add_bytes);
extern void sam_seg_BX_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(ub), unsigned add_bytes);
extern void sam_seg_QX_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(uy), unsigned add_bytes);
extern void sam_seg_fx_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(fx), unsigned add_bytes);
extern void sam_seg_BC_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(bc), unsigned add_bytes);
extern void sam_seg_other_seq (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, STRp(seq),unsigned add_bytes);
extern void sam_seg_GR_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(gr), unsigned add_bytes);
extern void sam_seg_GY_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(gy), unsigned add_bytes);
extern void sam_seg_GP_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t value, unsigned add_bytes);
extern void sam_seg_MP_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t value, unsigned add_bytes);

extern bool sam_seg_barcode_qual (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, SoloTags solo, uint8_t n_seps, STRp(qual), qSTRp (con_snip), MiniContainerP con, unsigned add_bytes);
extern void sam_seg_gene_name_id (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, STRp(value), unsigned add_bytes);

// biobambam2 stuff
extern uint32_t sam_get_QUAL_score (VBlockSAMP vb, STRp(qual));
extern void sam_seg_ms_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t ms, unsigned add_bytes);
extern void sam_seg_mc_i (VBlockSAMP vb, int64_t mc, unsigned add_bytes);

// Agilent stuff
extern void agilent_seg_initialize (VBlockP vb);
extern void agilent_seg_RX (VBlockP vb, ContextP ctx, STRp(rx), unsigned add_bytes); // RX and QX are also in fastq_private.h.
extern void agilent_seg_QX (VBlockP vb, ContextP ctx, STRp(qx), unsigned add_bytes); 

// xcons stuff
extern void sam_seg_xcons_YY_i (VBlockSAMP vb, int64_t yy, unsigned add_bytes);
extern void sam_seg_xcons_XC_i (VBlockSAMP vb, int64_t xc, unsigned add_bytes);
extern void sam_seg_xcons_XO_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t xo, unsigned add_bytes);

// NanoSeq stuff
extern void sam_seg_rb (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(rb), unsigned add_bytes);
extern void sam_seg_mb (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(mb), unsigned add_bytes);

// Abra2 stuff
extern void sam_abra2_zip_initialize (void);
extern void sam_abra2_seg_initialize (VBlockSAMP vb);
extern void sam_seg_ABRA2_YA_Z (VBlockSAMP vb, STRp(field), unsigned add_bytes);
extern void sam_seg_ABRA2_YO_Z (VBlockSAMP vb, STRp(field), unsigned add_bytes);

// SAG stuff
static inline bool sam_seg_has_sag_by_SA (VBlockSAMP vb)    { return  IS_SAG_SA && ((IS_DEPN(vb) && vb->sa_aln) || IS_PRIM(vb)); }
static inline bool sam_seg_has_sag_by_nonSA (VBlockSAMP vb) { return !IS_SAG_SA && ((IS_DEPN(vb) && vb->sag)    || IS_PRIM(vb)); }

extern void sam_sag_zip_init_vb (VBlockSAMP vb);
extern void sam_seg_gc_initialize (VBlockSAMP vb);
extern bool sam_seg_is_gc_line (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(alignment), bool is_bam);
extern bool sam_seg_prim_add_sag (VBlockSAMP vb, ZipDataLineSAMP dl, uint16_t num_alns/* inc primary aln*/, bool is_bam);
extern int32_t sam_seg_prim_add_sag_SA (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(sa), int64_t this_nm, bool is_bam);
extern void sam_seg_prim_add_sag_NH (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t nh);
extern void sam_seg_prim_add_sag_CC (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t nh);
extern void sam_seg_prim_add_sag_SOLO (VBlockSAMP vb, ZipDataLineSAMP dl);
extern void sam_seg_sag_stuff (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(textual_cigar), rom textual_seq, bool is_bam);
extern void sam_zip_gc_after_compute_main (VBlockSAMP vb);
extern void sam_sa_prim_initialize_ingest(void);
extern void sam_zip_prim_ingest_vb (VBlockSAMP vb);
extern void sam_seg_against_sa_group (VBlockSAMP vb, ContextP ctx, uint32_t add_bytes);
extern void sam_seg_against_sa_group_int (VBlockSAMP vb, ContextP ctx, int64_t parameter, uint32_t add_bytes);
extern Sag *sam_piz_get_prim_vb_first_sa_grp (VBlockSAMP vb);
extern void sam_show_sag (void);
extern bool sam_is_sag_loaded (void);
extern void sam_sag_by_flag_scan_for_depn (void);
extern bool sam_might_have_saggies_in_other_VBs (VBlockSAMP vb, ZipDataLineSAMP dl, int32_t n_alns/*0 if unknown*/);
extern void scan_index_qnames_seg (VBlockSAMP vb);
extern uint32_t sam_piz_get_plsg_i (VBIType vb_i);

typedef struct { char s[1024]; } ShowAln;
extern void sam_show_sag_one_grp (SAGroup grp_i);
extern ShowAln sam_show_sag_one_aln (const Sag *g, const SAAln *a);

typedef void (*ArrayItemCallback) (VBlockSAMP vb, ContextP ctx, void *cb_param, void *array, uint32_t array_len);
extern void sam_seg_array_one_ctx (VBlockSAMP vb, ZipDataLineSAMP dl, DictId dict_id, uint8_t type, rom array, int array_len, ArrayItemCallback callback, void *cb_param, PizSpecialReconstructor length_predictor);

#define SAM_PIZ_HAS_SAG (((VBlockSAMP)vb)->sag && ((VBlockSAMP)vb)->sag_line_i == vb->line_i + 1)

extern const char aux_sep_by_type[2][256];

typedef enum        { HD_SO_UNKNOWN, HD_SO_UNSORTED, HD_SO_QUERYNAME, HD_SO_COORDINATE } HdSoType;
#define HD_SO_NAMES { "unknown",     "unsorted",     "queryname",     "coordinate" }
extern HdSoType sam_hd_so;

typedef enum        { HD_GO_UNKNOWN, HD_GO_NONE, HD_GO_QUERY, HD_GO_REFERENCE } HdGoType;
#define HD_GO_NAMES { "unknown",     "none",     "query",     "reference"     }
extern HdGoType sam_hd_go;

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
        if (IN_RANGE (n, lt_min (test[i]), lt_max (test[i])))
            return lt_desc[test[i]].sam_type;
    
    return 0; // number out of range
}

extern DictId sam_seg_aux_field (VBlockSAMP vb, ZipDataLineSAMP dl, bool is_bam, rom tag, char bam_type, char bam_array_subtype, STRp(value), ValueType numeric, int16_t idx);

typedef struct { char s[200]; } DisFlagsStr;
extern DisFlagsStr sam_dis_flags (SamFlags flags);

// -------------------
// SAM-private globals
// -------------------

extern const uint8_t aux_width[256];

extern rom ERR_ANALYZE_RANGE_NOT_AVAILABLE, ERR_ANALYZE_DEPN_NOT_IN_REF, ERR_ANALYZE_INCORRECT_REF_BASE;

eSTRl(copy_GX_snip);
eSTRl(copy_sn_snip);
eSTRl(copy_POS_snip);
eSTRl(copy_Q1NAME_int);
eSTRl(copy_Q2NAME_int);
eSTRl(copy_Q3NAME_int);
eSTRl(copy_mate_CIGAR_snip);
eSTRl(copy_mate_MAPQ_snip);
eSTRl(copy_mate_MQ_snip);
eSTRl(copy_mate_PQ_snip);
eSTRl(XA_lookback_snip);
eSTRl(TX_lookback_snip);
eSTRl(AN_lookback_snip);
eSTRl(copy_mate_YS_snip);
eSTRl(copy_mate_AS_snip);
eSTRl(copy_mate_PNEXT_snip);
eSTRl(copy_saggy_PNEXT_snip);
eSTRl(copy_mate_POS_snip);
eSTRl(copy_mate_ms_snip);
eSTRl(copy_mate_nM_snip);
eSTRl(copy_buddy_NH_snip);
eSTRl(copy_mate_rb_snip);
eSTRl(copy_mate_mb_snip);
eSTRl(copy_mate_GP_snip);
eSTRl(copy_mate_MP_snip);
eSTRl(redirect_to_CR_X_snip);
eSTRl(redirect_to_GR_X_snip);
eSTRl(redirect_to_GY_X_snip);
eSTRl(redirect_to_RX_X_snip);
eSTRl(XC_snip);

eSTRl_ARRAY(copy_buddy_Z_snip, NUM_MATED_Z_TAGS, 30);
