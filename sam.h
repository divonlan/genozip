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

#pragma GENDICT_PREFIX SAM

// Fields
#pragma GENDICT SAM_RNAME=DTYPE_FIELD=RNAME // RNAME must be first
#pragma GENDICT SAM_QNAME=DTYPE_FIELD=QNAME
#pragma GENDICT SAM_FLAG=DTYPE_FIELD=FLAG
#pragma GENDICT SAM_POS=DTYPE_FIELD=POS
#pragma GENDICT SAM_MAPQ=DTYPE_FIELD=MAPQ
#pragma GENDICT SAM_CIGAR=DTYPE_FIELD=CIGAR
#pragma GENDICT SAM_RNEXT=DTYPE_FIELD=RNEXT
#pragma GENDICT SAM_PNEXT=DTYPE_FIELD=PNEXT
#pragma GENDICT SAM_TLEN=DTYPE_FIELD=TLEN
#pragma GENDICT SAM_OPTIONAL=DTYPE_FIELD=OPTIONAL
#pragma GENDICT SAM_SQBITMAP=DTYPE_FIELD=SQBITMAP
#pragma GENDICT SAM_NONREF=DTYPE_FIELD=NONREF    // these 4 fields must be in this order, right after SAM_SQBITMAP
#pragma GENDICT SAM_NONREF_X=DTYPE_FIELD=NONREF_X
#pragma GENDICT SAM_GPOS=DTYPE_FIELD=GPOS
#pragma GENDICT SAM_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT SAM_QUAL=DTYPE_FIELD=QUAL 
#pragma GENDICT SAM_DOMQRUNS=DTYPE_FIELD=DOMQRUNS // must be right after SAM_QUAL
#pragma GENDICT SAM_EOL=DTYPE_FIELD=EOL
#pragma GENDICT SAM_BAM_BIN=DTYPE_FIELD=BAM_BIN
#pragma GENDICT SAM_TOPLEVEL=DTYPE_FIELD=TOPLEVEL // must be called TOPLEVEL
#pragma GENDICT SAM_TOP2BAM=DTYPE_FIELD=TOP2BAM
#pragma GENDICT SAM_TOP2FQ=DTYPE_FIELD=TOP2FQ
#pragma GENDICT SAM_TAXID=DTYPE_FIELD=TAXID

#pragma GENDICT OPTION_AM=DTYPE_2=AM:i
#pragma GENDICT OPTION_AS=DTYPE_2=AS:i
#pragma GENDICT OPTION_CC=DTYPE_2=CC:Z
#pragma GENDICT OPTION_CM=DTYPE_2=CM:i
#pragma GENDICT OPTION_E2=DTYPE_2=E2:Z
#pragma GENDICT OPTION_2NONREF=DTYPE_2=N2ONREF // these 4 fields must be in this order, right after OPTION_E2
#pragma GENDICT OPTION_N2ONREFX=DTYPE_2=n2ONREFX
#pragma GENDICT OPTION_2GPOS=DTYPE_FIELD=G2POS
#pragma GENDICT OPTION_S2TRAND=DTYPE_2=S2TRAND
#pragma GENDICT OPTION_FI=DTYPE_2=FI:i
#pragma GENDICT OPTION_H0=DTYPE_2=H0:i
#pragma GENDICT OPTION_H1=DTYPE_2=H1:i
#pragma GENDICT OPTION_H2=DTYPE_2=H2:i
#pragma GENDICT OPTION_LB=DTYPE_2=LB:Z
#pragma GENDICT OPTION_MC=DTYPE_2=MC:Z
#pragma GENDICT OPTION_MD=DTYPE_2=MD:Z
#pragma GENDICT OPTION_MQ=DTYPE_2=MQ:i
#pragma GENDICT OPTION_NH=DTYPE_2=NH:i
#pragma GENDICT OPTION_NM=DTYPE_2=NM:i
#pragma GENDICT OPTION_OA=DTYPE_2=OA:Z
#pragma GENDICT OPTION_OC=DTYPE_2=OC:Z
#pragma GENDICT OPTION_PG=DTYPE_2=PG:Z
#pragma GENDICT OPTION_PQ=DTYPE_2=PQ:i
#pragma GENDICT OPTION_PU=DTYPE_2=PU:Z
#pragma GENDICT OPTION_RG=DTYPE_2=RG:Z
#pragma GENDICT OPTION_SA=DTYPE_2=SA:Z
#pragma GENDICT OPTION_SM=DTYPE_2=SM:i
#pragma GENDICT OPTION_TC=DTYPE_2=TC:i
#pragma GENDICT OPTION_UQ=DTYPE_2=UQ:i
#pragma GENDICT OPTION_U2=DTYPE_2=U2:Z
#pragma GENDICT OPTION_D2OMQRUN=DTYPE_2=D2OMQRUN // must be right after OPTION_U2
                
// Ion Torrent flow signal array
#pragma GENDICT OPTION_ZM=DTYPE_2=ZM:B

// bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
#pragma GENDICT OPTION_X0=DTYPE_2=X0:i 
#pragma GENDICT OPTION_X1=DTYPE_2=X1:i 
#pragma GENDICT OPTION_XA=DTYPE_2=XA:Z 
#pragma GENDICT OPTION_XA_RNAME=DTYPE_2=X0ARNAME 
#pragma GENDICT OPTION_XN=DTYPE_2=XN:i 
#pragma GENDICT OPTION_XM=DTYPE_2=XM:i 
#pragma GENDICT OPTION_XO=DTYPE_2=XO:i
#pragma GENDICT OPTION_XG=DTYPE_2=XG:i 
#pragma GENDICT OPTION_XS=DTYPE_2=XS:i 
#pragma GENDICT OPTION_XE=DTYPE_2=XE:i

// biobambam tags
#pragma GENDICT OPTION_mc=DTYPE_2=mc:i
#pragma GENDICT OPTION_ms=DTYPE_2=ms:i

// added by GATK's BQSR (Base Quality Score Recalibration)
#pragma GENDICT OPTION_BD=DTYPE_2=BD:Z // not used in newer versions of GATK
#pragma GENDICT OPTION_BI=DTYPE_2=BI:Z // not used in newer versions of GATK
#pragma GENDICT OPTION_BD_BI=DTYPE_2=BD_BI

// our private dictionary for + or 0 strands
#pragma GENDICT OPTION_STRAND=DTYPE_2=@STRAND
#pragma GENDICT OPTION_RNAME=DTYPE_2=@RNAME
#pragma GENDICT OPTION_POS=DTYPE_2=@POS
#pragma GENDICT OPTION_CIGAR=DTYPE_2=@CIGAR
#pragma GENDICT OPTION_MAPQ=DTYPE_2=@MAPQ
#pragma GENDICT OPTION_TX=DTYPE_2=TX:i

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
#define SAM_DICT_ID_ALIASES                               \
    /*         alias                 maps to this ctx */  \
    { DT_SAM,  _OPTION_MC,          _OPTION_CIGAR      }, \
    { DT_SAM,  _OPTION_OC,          _OPTION_CIGAR      }, \

#define BAM_DICT_ID_ALIASES                               \
    { DT_BAM,  _OPTION_MC,          _OPTION_CIGAR      }, \
    { DT_BAM,  _OPTION_OC,          _OPTION_CIGAR      }, \

#define SAM_LOCAL_GET_LINE_CALLBACKS                      \
    { DT_SAM,  _OPTION_BD_BI,       sam_zip_bd_bi      }, \
    { DT_SAM,  _SAM_QUAL,           sam_zip_qual       }, \
    { DT_SAM,  _OPTION_U2,          sam_zip_u2         }, 

#define BAM_LOCAL_GET_LINE_CALLBACKS                      \
    { DT_BAM,  _OPTION_BD_BI,       sam_zip_bd_bi      }, \
    { DT_BAM,  _SAM_QUAL,           sam_zip_qual       }, \
    { DT_BAM,  _OPTION_U2,          sam_zip_u2         }, 


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
