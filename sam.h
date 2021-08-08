// ------------------------------------------------------------------
//   sam.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef SAM_INCLUDED
#define SAM_INCLUDED

#include "genozip.h"
#include "sections.h"

#define DTYPE_QNAME        DTYPE_1
#define DTYPE_SAM_OPTIONAL DTYPE_2

// Fields
#define _SAM_RNAME       DICT_ID_MAKEF_5 ("RNAME")
#define _SAM_QNAME       DICT_ID_MAKEF_5 ("QNAME")
#define _SAM_FLAG        DICT_ID_MAKEF_4 ("FLAG")
#define _SAM_POS         DICT_ID_MAKEF_3 ("POS")
#define _SAM_MAPQ        DICT_ID_MAKEF_4 ("MAPQ")
#define _SAM_CIGAR       DICT_ID_MAKEF_5 ("CIGAR")
#define _SAM_RNEXT       DICT_ID_MAKEF_5 ("RNEXT")
#define _SAM_PNEXT       DICT_ID_MAKEF_5 ("PNEXT")
#define _SAM_TLEN        DICT_ID_MAKEF_4 ("TLEN")
#define _SAM_OPTIONAL    DICT_ID_MAKEF_L ("OPTIONAL")
#define _SAM_SQBITMAP    DICT_ID_MAKEF_L ("SQBITMAP")
#define _SAM_NONREF      DICT_ID_MAKEF_6 ("NONREF")
#define _SAM_NONREF_X    DICT_ID_MAKEF_L ("NONREF_X")
#define _SAM_GPOS        DICT_ID_MAKEF_4 ("GPOS")
#define _SAM_STRAND      DICT_ID_MAKEF_6 ("STRAND")
#define _SAM_QUAL        DICT_ID_MAKEF_4 ("QUAL") 
#define _SAM_DOMQRUNS    DICT_ID_MAKEF_L ("DOMQRUNS")
#define _SAM_EOL         DICT_ID_MAKEF_3 ("EOL")
#define _SAM_BAM_BIN     DICT_ID_MAKEF_7 ("BAM_BIN")
#define _SAM_TOPLEVEL    DICT_ID_MAKEF_L (TOPLEVEL)
#define _SAM_TOP2BAM     DICT_ID_MAKEF_7 ("TOP2BAM")
#define _SAM_TOP2FQ      DICT_ID_MAKEF_6 ("TOP2FQ")
#define _SAM_E2_Z        DICT_ID_MAKEF_4 ("E2:Z") 
#define _SAM_2NONREF     DICT_ID_MAKEF_7 ("2NONREF")
#define _SAM_N2ONREFX    DICT_ID_MAKEF_L ("N2ONREFX")
#define _SAM_2GPOS       DICT_ID_MAKEF_5 ("2GPOS")
#define _SAM_S2TRAND     DICT_ID_MAKEF_7 ("S2TRAND")
#define _SAM_U2_Z        DICT_ID_MAKEF_4 ("U2:Z") 
#define _SAM_D2OMQRUN    DICT_ID_MAKEF_L ("D2OMQRUN")
#define _SAM_TAXID       DICT_ID_MAKEF_5 ("TAXID")

#define _OPTION_AM       DICT_ID_MAKE2_4 ("AM:i")
#define _OPTION_AS       DICT_ID_MAKE2_4 ("AS:i")
#define _OPTION_CC       DICT_ID_MAKE2_4 ("CC:Z")
#define _OPTION_BD       DICT_ID_MAKE2_4 ("BD:Z")
#define _OPTION_BI       DICT_ID_MAKE2_4 ("BI:Z")
#define _OPTION_BD_BI    DICT_ID_MAKE2_5 ("BD_BI")
#define _OPTION_CM       DICT_ID_MAKE2_4 ("CM:i")
#define _OPTION_E2       DICT_ID_MAKE2_4 ("E2:Z")
#define _OPTION_FI       DICT_ID_MAKE2_4 ("FI:i")
#define _OPTION_H0       DICT_ID_MAKE2_4 ("H0:i")
#define _OPTION_H1       DICT_ID_MAKE2_4 ("H1:i")
#define _OPTION_H2       DICT_ID_MAKE2_4 ("H2:i")
#define _OPTION_LB       DICT_ID_MAKE2_4 ("LB:Z")
#define _OPTION_MC       DICT_ID_MAKE2_4 ("MC:Z")
#define _OPTION_MD       DICT_ID_MAKE2_4 ("MD:Z")
#define _OPTION_MQ       DICT_ID_MAKE2_4 ("MQ:i")
#define _OPTION_NH       DICT_ID_MAKE2_4 ("NH:i")
#define _OPTION_NM       DICT_ID_MAKE2_4 ("NM:i")
#define _OPTION_OA       DICT_ID_MAKE2_4 ("OA:Z")
#define _OPTION_OC       DICT_ID_MAKE2_4 ("OC:Z")
#define _OPTION_PG       DICT_ID_MAKE2_4 ("PG:Z")
#define _OPTION_PQ       DICT_ID_MAKE2_4 ("PQ:i")
#define _OPTION_PU       DICT_ID_MAKE2_4 ("PU:Z")
#define _OPTION_RG       DICT_ID_MAKE2_4 ("RG:Z")
#define _OPTION_SA       DICT_ID_MAKE2_4 ("SA:Z")
#define _OPTION_SM       DICT_ID_MAKE2_4 ("SM:i")
#define _OPTION_TC       DICT_ID_MAKE2_4 ("TC:i")
#define _OPTION_UQ       DICT_ID_MAKE2_4 ("UQ:i")
#define _OPTION_U2       DICT_ID_MAKE2_4 ("U2:Z")
                
// Ion Torrent flow signal array
#define _OPTION_ZM       DICT_ID_MAKE2_4 ("ZM:B")

// bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
#define _OPTION_X0       DICT_ID_MAKE2_4 ("X0:i") 
#define _OPTION_X1       DICT_ID_MAKE2_4 ("X1:i") 
#define _OPTION_XA       DICT_ID_MAKE2_4 ("XA:Z") 
#define _OPTION_XA_RNAME DICT_ID_MAKE2_L ("X0ARNAME") 
#define _OPTION_XN       DICT_ID_MAKE2_4 ("XN:i") 
#define _OPTION_XM       DICT_ID_MAKE2_4 ("XM:i") 
#define _OPTION_XO       DICT_ID_MAKE2_4 ("XO:i")
#define _OPTION_XG       DICT_ID_MAKE2_4 ("XG:i") 
#define _OPTION_XS       DICT_ID_MAKE2_4 ("XS:i") 
#define _OPTION_XE       DICT_ID_MAKE2_4 ("XE:i")

// biobambam tags
#define _OPTION_mc       DICT_ID_MAKE2_4 ("mc:i")
#define _OPTION_ms       DICT_ID_MAKE2_4 ("ms:i")

// added by GATK's BQSR (Base Quality Score Recalibration)
#define _OPTION_BD       DICT_ID_MAKE2_4 ("BD:Z") // not used in newer versions of GATK
#define _OPTION_BI       DICT_ID_MAKE2_4 ("BI:Z") // not used in newer versions of GATK

// our private dictionary for + or 0 strands
#define _OPTION_STRAND   DICT_ID_MAKE2_7 ("@STRAND")
#define _OPTION_RNAME    DICT_ID_MAKE2_6 ("@RNAME")
#define _OPTION_POS      DICT_ID_MAKE2_4 ("@POS")
#define _OPTION_CIGAR    DICT_ID_MAKE2_6 ("@CIGAR")
#define _OPTION_MAPQ     DICT_ID_MAKE2_5 ("@MAPQ")
#define _OPTION_TX       DICT_ID_MAKE2_4 ("TX:i")

// did_i to dict_i mapping - only needed for did_i's referred to explicitly
// the CHROM field MUST be the first field (because of ctx_build_zf_ctx_from_contigs)
typedef enum { SAM_RNAME, SAM_QNAME, SAM_FLAG, SAM_POS, SAM_MAPQ, SAM_CIGAR, SAM_RNEXT, SAM_PNEXT, SAM_TLEN, SAM_OPTIONAL, 
               SAM_SQBITMAP, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND, 
               SAM_QUAL, SAM_DOMQRUNS, SAM_EOL, SAM_BAM_BIN,
               SAM_TOPLEVEL, SAM_TOP2BAM, SAM_TOP2FQ, 
               SAM_E2_Z, SAM_2NONREF, SAM_N2ONREFX, SAM_2GPOS, SAM_S2TRAND, // This fields need to be sequential. E2 and U2 are in primary fields instead of FORMAT due to historical reasons.
               SAM_U2_Z, SAM_D2OMQRUN,                                      // This fields need to be sequential
               SAM_TAXID,
               OPTION_BD_BI, OPTION_XA, OPTION_XA_RNAME,
               NUM_SAM_FIELDS } SamFields;

#define SAM_MAPPING { V(SAM_RNAME), V(SAM_QNAME), V(SAM_FLAG), V(SAM_POS), V(SAM_MAPQ), V(SAM_CIGAR), V(SAM_RNEXT), V(SAM_PNEXT), V(SAM_TLEN), \
                      V(SAM_OPTIONAL), V(SAM_SQBITMAP), V(SAM_NONREF), V(SAM_NONREF_X), V(SAM_GPOS), V(SAM_STRAND), V(SAM_QUAL), V(SAM_DOMQRUNS), V(SAM_EOL), \
                      V(SAM_BAM_BIN) ,V(SAM_TOPLEVEL), V(SAM_TOP2BAM), V(SAM_TOP2FQ), \
                      V(SAM_E2_Z), V(SAM_2NONREF), V(SAM_N2ONREFX), V(SAM_2GPOS), V(SAM_S2TRAND), V(SAM_U2_Z), V(SAM_D2OMQRUN), V(SAM_TAXID), \
                      V(OPTION_BD_BI), V(OPTION_XA), V(OPTION_XA_RNAME) }

// ZIP Stuff
COMPRESSOR_CALLBACK(sam_zip_qual)
COMPRESSOR_CALLBACK(sam_zip_u2)
COMPRESSOR_CALLBACK(sam_zip_bd_bi)
extern void sam_zip_initialize (void);
extern bool sam_zip_is_unaligned_line (const char *line, int len);
extern bool sam_zip_dts_flag (void);

// HEADER stuff
extern bool sam_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);

// SEG Stuff
extern void sam_zip_initialize (void);
extern void sam_seg_initialize (VBlockP vb);
extern void sam_seg_finalize (VBlockP vb);
extern bool sam_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern const char *sam_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern uint32_t sam_seg_seq_len_from_cigar (const char *cigar, unsigned cigar_len);
extern uint32_t sam_seg_get_seq_len_by_MD_field (const char *md_str, unsigned md_str_len);

// PIZ Stuff
extern bool sam_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
extern void sam_reconstruct_vb ();
extern void sam_reconstruct_seq (VBlockP vb, ContextP ctx, const char *unused, unsigned unused2);
extern void sam_set_FLAG_filter (const char *optarg);
extern void sam_set_MAPQ_filter (const char *optarg);

// BAM Stuff
extern void bam_seg_initialize (VBlockP vb);
extern int32_t bam_is_header_done (bool is_eof);
extern int32_t bam_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern void bam_read_vblock (VBlockP vb);
extern void bam_seg_initialize (VBlockP vb);
extern const char *bam_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);

// SAM-to-FASTQ stuff
CONTAINER_CALLBACK (sam_piz_container_cb);

// VB stuff
extern void sam_vb_release_vb();
extern void sam_vb_destroy_vb();
extern unsigned sam_vb_size (DataType dt);
extern unsigned sam_vb_zip_dl_size (void);

// Special - used for SAM & BAM
#define SAM_SPECIAL { sam_piz_special_CIGAR, sam_piz_special_TLEN, sam_piz_special_BD_BI, sam_piz_special_AS, \
                      sam_piz_special_MD, bam_piz_special_FLOAT, bam_piz_special_BIN, sam_piz_special_XA_POS }
SPECIAL (SAM, 0, CIGAR, sam_piz_special_CIGAR);
SPECIAL (SAM, 1, TLEN,  sam_piz_special_TLEN);
SPECIAL (SAM, 2, BDBI,  sam_piz_special_BD_BI);
SPECIAL (SAM, 3, AS,    sam_piz_special_AS);
SPECIAL (SAM, 4, MD,    sam_piz_special_MD);
SPECIAL (SAM, 5, FLOAT, bam_piz_special_FLOAT); // used in BAM to represent float optional values
SPECIAL (SAM, 6, BIN,   bam_piz_special_BIN);   
SPECIAL (SAM, 7, XA_POS, sam_piz_special_XA_POS);   
#define NUM_SAM_SPECIAL 8

// note: we can't alias RNEXT to RNAME, because we can't alias to CHROM - see comment in reconstruct_from_ctx_do 
#define SAM_DICT_ID_ALIASES \
    /*         alias                       maps to this ctx            */  \
    { DT_SAM,  _OPTION_MC,          _OPTION_CIGAR         }, \
    { DT_SAM,  _OPTION_OC,          _OPTION_CIGAR         }, \
    { DT_SAM,  _OPTION_E2,          _SAM_E2_Z             }, \
    { DT_SAM,  _OPTION_U2,          _SAM_U2_Z             },

#define SAM_LOCAL_GET_LINE_CALLBACKS  \
    { DT_SAM,   _OPTION_BD_BI,      sam_zip_bd_bi                }, \
    { DT_SAM,   _SAM_QUAL,          sam_zip_qual                 }, 

#define BAM_DICT_ID_ALIASES \
    /*         alias                       maps to this ctx            */  \
    { DT_BAM,  _OPTION_MC,          _OPTION_CIGAR         }, \
    { DT_BAM,  _OPTION_OC,          _OPTION_CIGAR         }, \
    { DT_BAM,  _OPTION_E2,          _SAM_E2_Z             }, \
    { DT_BAM,  _OPTION_U2,          _SAM_U2_Z             },

#define BAM_LOCAL_GET_LINE_CALLBACKS  \
    { DT_BAM,  _OPTION_BD_BI,       sam_zip_bd_bi                }, \
    { DT_BAM,  _SAM_QUAL,           sam_zip_qual                 }, 


// Important: Numbers (and order) of translators cannot be changed, as they are part of the file format
// (they are included in the TOP2BAM container)
// translator numbers must start from 1 - 0 is reserved for "none"
TRANSLATOR (SAM, BAM,   1,  I8,         container_translate_I8)   // reconstruct binary little endian functions
TRANSLATOR (SAM, BAM,   2,  U8,         container_translate_U8)   // 
TRANSLATOR (SAM, BAM,   3,  LTEN_I16,   container_translate_LTEN_I16) 
TRANSLATOR (SAM, BAM,   4,  LTEN_U16,   container_translate_LTEN_U16) 
TRANSLATOR (SAM, BAM,   5,  LTEN_I32,   container_translate_LTEN_I32) 
TRANSLATOR (SAM, BAM,   6,  LTEN_U32,   container_translate_LTEN_U32) 
TRANSLATOR (SAM, BAM,   7,  FLOAT,      sam_piz_sam2bam_FLOAT)      // reconstructs SAM-stored textual floating point as little endian 32bit float
TRANSLATOR (SAM, BAM,   8,  ARRAY_SELF, sam_piz_sam2bam_ARRAY_SELF) // remove the comma from the prefix that contains the type, eg "i,"->"i"
TRANSLATOR (SAM, BAM,   9,  RNAME,      sam_piz_sam2bam_RNAME)      // reconstructs the b250 index or -1 if "*"
TRANSLATOR (SAM, BAM,   10, POS,        sam_piz_sam2bam_POS)        // reconstructs Little Endian U32 0-based POS. 
TRANSLATOR (SAM, BAM,   11, SEQ,        sam_piz_sam2bam_SEQ)        // textual SEQ to BAM-format SEQ 
TRANSLATOR (SAM, BAM,   12, QUAL,       sam_piz_sam2bam_QUAL)       // textual QUAL to BAM-format QUAL 
TRANSLATOR (SAM, BAM,   13, TLEN,       sam_piz_sam2bam_TLEN)       // place TLEN last_value in BAM alignment 
TRANSLATOR (SAM, BAM,   14, OPTIONAL,   sam_piz_sam2bam_OPTIONAL)   // used up to v11, kept for for backward compatability as old files expect it
TRANSLATOR (SAM, BAM,   15, OPTIONAL_SELF, sam_piz_sam2bam_OPTIONAL_SELF) // transform prefixes in Optional Container from SAM to BAM format 
TRANSLATOR (SAM, FASTQ, 16, SEQ,        sam_piz_sam2fastq_SEQ)      // reverse-complement the sequence if needed, and drop if "*"
TRANSLATOR (SAM, FASTQ, 17, QUAL,       sam_piz_sam2fastq_QUAL)     // reverse the QUAL if reverse-complemented and drop fastq records with QUAL="*"
TRANSLATOR (SAM, FASTQ, 18, FLAG,       sam_piz_sam2fastq_FLAG)     // emit 1 if (FLAGS & 0x40) or 2 of (FLAGS & 0x80)

#define NUM_SAM_TRANS   19 // including "none"
#define SAM_TRANSLATORS { NULL /* none */, container_translate_I8, container_translate_U8, container_translate_LTEN_I16, \
                          container_translate_LTEN_U16, container_translate_LTEN_I32, container_translate_LTEN_U32, \
                          sam_piz_sam2bam_FLOAT, sam_piz_sam2bam_ARRAY_SELF, sam_piz_sam2bam_RNAME, sam_piz_sam2bam_POS, sam_piz_sam2bam_SEQ, \
                          sam_piz_sam2bam_QUAL, sam_piz_sam2bam_TLEN, sam_piz_sam2bam_OPTIONAL, sam_piz_sam2bam_OPTIONAL_SELF, \
                          sam_piz_sam2fastq_SEQ, sam_piz_sam2fastq_QUAL, sam_piz_sam2fastq_FLAG }

TXTHEADER_TRANSLATOR (txtheader_bam2sam);
TXTHEADER_TRANSLATOR (txtheader_sam2bam);
TXTHEADER_TRANSLATOR (txtheader_sam2fq);

#define SAM_CONTIG_FMT "@SQ	SN:%.*s	LN:%"PRId64

#endif
