// ------------------------------------------------------------------
//   fastq.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FASTQ_INCLUDED
#define FASTQ_INCLUDED

#include "genozip.h"

// ZIP Stuff
COMPRESSOR_CALLBACK(fastq_zip_get_start_len_line_i_qual)

// SEG Stuff
extern void fastq_seg_initialize();
extern void fastq_seg_finalize();
extern const char *fastq_seg_txt_line();

// PIZ Stuff
extern bool fastq_piz_read_one_vb (VBlockP vb, SectionListEntryP sl);
extern void fastq_piz_reconstruct_vb();
extern bool fastq_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
extern void fastq_piz_reconstruct_seq (VBlockP vb, ContextP bitmap_ctx, const char *seq_len_str, unsigned seq_len_str_len);

// VBlock stuff
extern void fast_vb_release_vb();
extern unsigned fast_vb_size (void);
extern unsigned fast_vb_zip_dl_size (void);

// file pairing (--pair) stuff
extern bool fastq_read_pair_1_data (VBlockP vb, uint32_t first_vb_i_of_pair_1, uint32_t last_vb_i_of_pair_1);
extern bool fastq_txtfile_have_enough_lines (VBlockP vb, uint32_t *unconsumed_len);
extern uint32_t fastq_get_pair_vb_i (VBlockP vb);

#define FASTQ_DICT_ID_ALIASES \
    /*          alias                       maps to this ctx          */  \
    { DT_FASTQ, &dict_id_fields[FASTQ_E2L], &dict_id_fields[FASTQ_E1L] }, /* note: the lowest did_i must be the non-alias */ \
    { DT_FASTQ, &dict_id_fields[FASTQ_E3L], &dict_id_fields[FASTQ_E1L] }, \
    { DT_FASTQ, &dict_id_fields[FASTQ_E4L], &dict_id_fields[FASTQ_E1L] }, 

#define FASTQ_LOCAL_GET_LINE_CALLBACKS  \
    { DT_FASTQ, &dict_id_fields[FASTQ_QUAL], fastq_zip_get_start_len_line_i_qual },

#endif
