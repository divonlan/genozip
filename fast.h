// ------------------------------------------------------------------
//   fast.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FAST_INCLUDED
#define FAST_INCLUDED

#include "genozip.h"

// ZIP Stuff
COMPRESSOR_CALLBACK(fast_zip_get_start_len_line_i_seq)
COMPRESSOR_CALLBACK(fastq_zip_get_start_len_line_i_qual)

// SEG Stuff
extern void fasta_seg_initialize();
extern void fastq_seg_initialize();
extern const char *fastq_seg_txt_line();
extern const char *fasta_seg_txt_line();

// PIZ Stuff
extern bool fast_piz_read_one_vb (VBlockP vb, SectionListEntryP sl);
extern void fasta_piz_reconstruct_vb(); // no parameter - implicit casting of VBlockP
extern void fastq_piz_reconstruct_vb();

// VBlock stuff
extern void fast_vb_release_vb();
extern unsigned fast_vb_size (void);
extern unsigned fast_vb_zip_dl_size (void);

#define FASTA_SPECIAL { fasta_piz_special_SEQ, fasta_piz_special_COMMENT, fasta_piz_special_DESC }
SPECIAL (FASTA, 0, SEQ, fasta_piz_special_SEQ);
SPECIAL (FASTA, 1, COMMENT, fasta_piz_special_COMMENT);
SPECIAL (FASTA, 2, DESC, fasta_piz_special_DESC);
#define NUM_FASTA_SPECIAL 3

#define FASTQ_DICT_ID_ALIASES \
    /*          alias                           maps to this ctx          */  \
    { DT_FASTQ, &dict_id_fields[FASTQ_E2L],     &dict_id_fields[FASTQ_E1L]  }, /* note: the lowest did_i must be the non-alias */ \
    { DT_FASTQ, &dict_id_fields[FASTQ_E3L],     &dict_id_fields[FASTQ_E1L]  }, \
    { DT_FASTQ, &dict_id_fields[FASTQ_E4L],     &dict_id_fields[FASTQ_E1L]  }, 

#define FASTA_DICT_ID_ALIASES 

#define FASTQ_LOCAL_COMPRESSOR_CALLBACKS  \
    { DT_FASTQ, &dict_id_fields[FASTQ_QUAL], fastq_zip_get_start_len_line_i_qual }, \
    { DT_FASTQ, &dict_id_fields[FASTQ_SEQ],  fast_zip_get_start_len_line_i_seq  }, 

#define FASTA_LOCAL_COMPRESSOR_CALLBACKS  \
    { DT_FASTA, &dict_id_FASTA_SEQ,          fast_zip_get_start_len_line_i_seq  }, 

#define dict_id_is_fast_desc_sf dict_id_is_type_1
#define dict_id_fast_desc_sf dict_id_type_1

#endif
