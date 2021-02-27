// ------------------------------------------------------------------
//   sam.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SAM_INCLUDED
#define SAM_INCLUDED

#include "genozip.h"
#include "sections.h"

// ZIP Stuff
COMPRESSOR_CALLBACK(sam_zip_qual)
COMPRESSOR_CALLBACK(sam_zip_u2)
COMPRESSOR_CALLBACK(sam_zip_bd_bi)
extern void sam_zip_initialize (void);
extern bool sam_zip_is_unaligned_line (const char *line, int len);
extern bool sam_zip_dts_flag (void);

// HEADER stuff
extern bool sam_header_inspect (BufferP txt_header);
extern void sam_header_get_contigs (ConstBufferP *contigs_dict, ConstBufferP *contigs);

// SEG Stuff
extern void sam_seg_initialize (VBlockP vb);
extern void sam_header_finalize (void);
extern void sam_seg_finalize (VBlockP vb);
extern const char *sam_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern uint32_t sam_seg_seq_len_from_cigar (const char *cigar, unsigned cigar_len);
extern uint32_t sam_seg_get_seq_len_by_MD_field (const char *md_str, unsigned md_str_len);

// PIZ Stuff
extern bool sam_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
extern void sam_reconstruct_vb ();
extern void sam_reconstruct_seq (VBlockP vb, ContextP ctx, const char *unused, unsigned unused2);
extern void sam_piz_show_sex (void);
extern void sam_piz_show_sex_count_one_vb (VBlockP vb);

// BAM Stuff
extern void bam_seg_initialize (VBlockP vb);
extern int32_t bam_is_header_done (void);
extern int32_t bam_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern void txtheader_bam2sam (BufferP txt);
extern void bam_read_vblock (VBlockP vb);
extern void bam_seg_initialize (VBlockP vb);
extern const char *bam_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);

// SAM-to-FASTQ stuff
CONTAINER_FILTER_FUNC (sam_piz_sam2fq_filter);

// VB stuff
extern void sam_vb_release_vb();
extern void sam_vb_destroy_vb();
extern unsigned sam_vb_size (void);
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

// SAM field types 
#define DTYPE_QNAME    DTYPE_1
#define DTYPE_SAM_OPTIONAL DTYPE_2

// note: we can't alias RNEXT to RNAME, because we can't alias to CHROM - see comment in reconstruct_from_ctx_do 
#define SAM_DICT_ID_ALIASES \
    /*         alias                        maps to this ctx          */     \
    { DT_SAM,  &dict_id_OPTION_MC,          &dict_id_OPTION_CIGAR         }, \
    { DT_SAM,  &dict_id_OPTION_OC,          &dict_id_OPTION_CIGAR         }, \
    { DT_SAM,  &dict_id_OPTION_E2,          &dict_id_fields[SAM_E2_Z]     }, \
    { DT_SAM,  &dict_id_OPTION_U2,          &dict_id_fields[SAM_U2_Z]     },

#define SAM_LOCAL_GET_LINE_CALLBACKS  \
    { DT_SAM,   &dict_id_OPTION_BD_BI,      sam_zip_bd_bi }, \
    { DT_SAM,   &dict_id_fields[SAM_QUAL],  sam_zip_qual  }, 

#define BAM_DICT_ID_ALIASES \
    /*         alias                        maps to this ctx              */ \
    { DT_BAM,  &dict_id_OPTION_MC,          &dict_id_OPTION_CIGAR         }, \
    { DT_BAM,  &dict_id_OPTION_OC,          &dict_id_OPTION_CIGAR         }, \
    { DT_BAM,  &dict_id_OPTION_E2,          &dict_id_fields[SAM_E2_Z]     }, \
    { DT_BAM,  &dict_id_OPTION_U2,          &dict_id_fields[SAM_U2_Z]     },

#define BAM_LOCAL_GET_LINE_CALLBACKS  \
    { DT_BAM,   &dict_id_OPTION_BD_BI,      sam_zip_bd_bi }, \
    { DT_BAM,   &dict_id_fields[SAM_QUAL],  sam_zip_qual  }, 


// Important: Numbers (and order) of translators cannot be changed, as they are part of the file format
// (included in the TOP2BAM container)
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
TRANSLATOR (SAM, BAM,   14, OPTIONAL,   sam_piz_sam2bam_OPTIONAL)   // set block_size after Optional reconstruction
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

#endif
