// ------------------------------------------------------------------
//   fastq.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "sections.h"

#pragma GENDICT_PREFIX FASTQ
#pragma GENDICT FASTQ_CONTIG=DTYPE_FIELD=CONTIG // must be first - did_i=0=CHROM
#pragma GENDICT FASTQ_DESC=DTYPE_FIELD=DESC     // Q?NAME must immediately follow
#pragma GENDICT FASTQ_Q0NAME=DTYPE_1=Q0NAME     // MAX_QNAME_ITEMS fixed qname items must have a did_i directly after container's (MUST be the same dict_id as in sam.h)
#pragma GENDICT FASTQ_Q1NAME=DTYPE_1=Q1NAME 
#pragma GENDICT FASTQ_Q2NAME=DTYPE_1=Q2NAME
#pragma GENDICT FASTQ_Q3NAME=DTYPE_1=Q3NAME
#pragma GENDICT FASTQ_Q4NAME=DTYPE_1=Q4NAME
#pragma GENDICT FASTQ_Q5NAME=DTYPE_1=Q5NAME
#pragma GENDICT FASTQ_Q6NAME=DTYPE_1=Q6NAME 
#pragma GENDICT FASTQ_Q7NAME=DTYPE_1=Q7NAME 
#pragma GENDICT FASTQ_QmatNAME=DTYPE_1=QmatNAME // QmatNAME reserved for mate number (always the last dict_id in the container)
#pragma GENDICT FASTQ_E1L=DTYPE_FIELD=E1L
#pragma GENDICT FASTQ_SQBITMAP=DTYPE_FIELD=SQBITMAP
#pragma GENDICT FASTQ_NONREF=DTYPE_FIELD=NONREF
#pragma GENDICT FASTQ_NONREF_X=DTYPE_FIELD=NONREF_X
#pragma GENDICT FASTQ_GPOS=DTYPE_FIELD=GPOS
#pragma GENDICT FASTQ_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT FASTQ_E2L=DTYPE_FIELD=E2L
#pragma GENDICT FASTQ_QUAL=DTYPE_FIELD=QUAL 
#pragma GENDICT FASTQ_DOMQRUNS=DTYPE_FIELD=DOMQRUNS
#pragma GENDICT FASTQ_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT FASTQ_TAXID=DTYPE_FIELD=TAXID
#pragma GENDICT FASTQ_LINE3=DTYPE_FIELD=LINE3

// Txtfile stuff
extern int32_t fastq_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern bool fastq_txtfile_have_enough_lines (VBlockP vb, uint32_t *unconsumed_len, uint32_t *my_lines, uint32_t *her_lines);

// ZIP Stuff
extern void fastq_zip_initialize (void);
extern void fastq_zip_read_one_vb (VBlockP vb);
extern bool fastq_zip_dts_flag (void);
COMPRESSOR_CALLBACK(fastq_zip_qual);

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
extern void fastq_reconstruct_seq (VBlockP vb, ContextP bitmap_ctx, STRp(seq_len_str));
extern CONTAINER_CALLBACK (fastq_piz_container_cb);

// VBlock stuff
extern void fastq_vb_release_vb();
extern void fastq_vb_destroy_vb();
extern unsigned fastq_vb_size (DataType dt);
extern unsigned fastq_vb_zip_dl_size (void);

// file pairing (--pair) stuff
extern bool fastq_read_pair_1_data (VBlockP vb, uint32_t pair_vb_i, bool must_have);
extern uint32_t fastq_get_pair_vb_i (VBlockP vb);

#define FASTQ_LOCAL_GET_LINE_CALLBACKS  \
    { DT_FASTQ, _FASTQ_QUAL, fastq_zip_qual },
