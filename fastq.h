// ------------------------------------------------------------------
//   fastq.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef FASTQ_INCLUDED
#define FASTQ_INCLUDED

#include "genozip.h"
#include "sections.h"

#define _FASTQ_CONTIG   DICT_ID_MAKEF_6 ("CONTIG")
#define _FASTQ_DESC     DICT_ID_MAKEF_4 ("DESC")
#define _FASTQ_E1L      DICT_ID_MAKEF_3 ("E1L")
#define _FASTQ_SQBITMAP DICT_ID_MAKEF_L ("SQBITMAP")
#define _FASTQ_NONREF   DICT_ID_MAKEF_6 ("NONREF")
#define _FASTQ_NONREF_X DICT_ID_MAKEF_L ("NONREF_X")
#define _FASTQ_GPOS     DICT_ID_MAKEF_4 ("GPOS")
#define _FASTQ_STRAND   DICT_ID_MAKEF_6 ("STRAND")
#define _FASTQ_E2L      DICT_ID_MAKEF_3 ("E2L")
#define _FASTQ_QUAL     DICT_ID_MAKEF_4 ("QUAL") 
#define _FASTQ_DOMQRUNS DICT_ID_MAKEF_L ("DOMQRUNS")
#define _FASTQ_TOPLEVEL DICT_ID_MAKEF_L (TOPLEVEL)
#define _FASTQ_TAXID    DICT_ID_MAKEF_5 ("TAXID")
#define _FASTQ_LINE3    DICT_ID_MAKEF_5 ("LINE3")

typedef enum { FASTQ_CONTIG /* copied from reference */, FASTQ_DESC, FASTQ_E1L, FASTQ_SQBITMAP, FASTQ_NONREF, FASTQ_NONREF_X, FASTQ_GPOS, FASTQ_STRAND, FASTQ_E2L, FASTQ_QUAL, FASTQ_DOMQRUNS, FASTQ_TOPLEVEL, FASTQ_TAXID, FASTQ_LINE3, NUM_FASTQ_FIELDS } FastqFields;

#define FASTQ_MAPPING { V(FASTQ_CONTIG), V(FASTQ_DESC), V(FASTQ_E1L), V(FASTQ_SQBITMAP), V(FASTQ_NONREF), V(FASTQ_NONREF_X), V(FASTQ_GPOS), V(FASTQ_STRAND), V(FASTQ_E2L), V(FASTQ_QUAL), V(FASTQ_DOMQRUNS), V(FASTQ_TOPLEVEL), V(FASTQ_TAXID), V(FASTQ_LINE3) }

// Txtfile stuff
extern int32_t fastq_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern bool fastq_txtfile_have_enough_lines (VBlockP vb, uint32_t *unconsumed_len, uint32_t *my_lines, uint32_t *her_lines);

// ZIP Stuff
extern void fastq_zip_initialize (void);
extern void fastq_zip_read_one_vb (VBlockP vb);
extern bool fastq_zip_dts_flag (void);
COMPRESSOR_CALLBACK(fastq_zip_qual)

// SEG Stuff
extern void fastq_seg_initialize();
extern void fastq_seg_finalize();
extern bool fastq_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern const char *fastq_seg_txt_line();

// PIZ Stuff
extern bool fastq_piz_is_paired (void);
extern bool fastq_piz_read_one_vb (VBlockP vb, Section sl);
CONTAINER_FILTER_FUNC (fastq_piz_filter);
extern bool fastq_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
extern void fastq_reconstruct_seq (VBlockP vb, ContextP bitmap_ctx, const char *seq_len_str, unsigned seq_len_str_len);
extern CONTAINER_CALLBACK (fastq_piz_container_cb);

// VBlock stuff
extern void fastq_vb_release_vb();
extern void fastq_vb_destroy_vb();
extern unsigned fastq_vb_size (DataType dt);
extern unsigned fastq_vb_zip_dl_size (void);

// file pairing (--pair) stuff
extern bool fastq_read_pair_1_data (VBlockP vb, uint32_t pair_vb_i, bool must_have);
extern uint32_t fastq_get_pair_vb_i (VBlockP vb);

#define FASTQ_DICT_ID_ALIASES 

#define FASTQ_LOCAL_GET_LINE_CALLBACKS  \
    { DT_FASTQ, _FASTQ_QUAL, fastq_zip_qual },

#endif
