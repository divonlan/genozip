// ------------------------------------------------------------------
//   fastq.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
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
#pragma GENDICT FASTQ_Q8NAME=DTYPE_1=Q8NAME 
#pragma GENDICT FASTQ_Q9NAME=DTYPE_1=Q9NAME 
#pragma GENDICT FASTQ_QmNAME=DTYPE_1=QmNAME     // QmNAME reserved for mate number (always the last dict_id in the container)
#pragma GENDICT FASTQ_QNAME2=DTYPE_1=QNAME2     // QNAME embedded in DESC (QNAME2 items immediately follow)
#pragma GENDICT FASTQ_Q0NAME2=DTYPE_1=q0NAME    
#pragma GENDICT FASTQ_Q1NAME2=DTYPE_1=q1NAME 
#pragma GENDICT FASTQ_Q2NAME2=DTYPE_1=q2NAME
#pragma GENDICT FASTQ_Q3NAME2=DTYPE_1=q3NAME
#pragma GENDICT FASTQ_Q4NAME2=DTYPE_1=q4NAME
#pragma GENDICT FASTQ_Q5NAME2=DTYPE_1=q5NAME
#pragma GENDICT FASTQ_Q6NAME2=DTYPE_1=q6NAME 
#pragma GENDICT FASTQ_Q7NAME2=DTYPE_1=q7NAME 
#pragma GENDICT FASTQ_Q8NAME2=DTYPE_1=q8NAME 
#pragma GENDICT FASTQ_Q9NAME2=DTYPE_1=q9NAME 
#pragma GENDICT FASTQ_E1L=DTYPE_FIELD=E1L
#pragma GENDICT FASTQ_SQBITMAP=DTYPE_FIELD=SQBITMAP
#pragma GENDICT FASTQ_NONREF=DTYPE_FIELD=NONREF
#pragma GENDICT FASTQ_NONREF_X=DTYPE_FIELD=NONREF_X
#pragma GENDICT FASTQ_GPOS=DTYPE_FIELD=GPOS
#pragma GENDICT FASTQ_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT FASTQ_E2L=DTYPE_FIELD=E2L
#pragma GENDICT FASTQ_QUAL=DTYPE_FIELD=QUAL 
#pragma GENDICT FASTQ_DOMQRUNS=DTYPE_FIELD=DOMQRUNS // must be 1 after QUAL. Also used by LONGR. For backwards compatability, we can never change its name.
#pragma GENDICT FASTQ_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT FASTQ_TAXID=DTYPE_FIELD=TAXID

#pragma GENDICT FASTQ_LINE3=DTYPE_FIELD=LINE3
#define MAX_LINE3_ITEMS 9
#pragma GENDICT FASTQ_T0HIRD=DTYPE_1=t0NAME    // must be directly after FASTQ_LINE3
#pragma GENDICT FASTQ_T1HIRD=DTYPE_1=t1NAME 
#pragma GENDICT FASTQ_T2HIRD=DTYPE_1=t2NAME
#pragma GENDICT FASTQ_T3HIRD=DTYPE_1=t3NAME
#pragma GENDICT FASTQ_T4HIRD=DTYPE_1=t4NAME
#pragma GENDICT FASTQ_T5HIRD=DTYPE_1=t5NAME
#pragma GENDICT FASTQ_T6HIRD=DTYPE_1=t6NAME 
#pragma GENDICT FASTQ_T7HIRD=DTYPE_1=t7NAME 
#pragma GENDICT FASTQ_COPY_Q=DTYPE_1=tQcopy 

// Txtfile stuff
extern int32_t fastq_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern bool fastq_txtfile_have_enough_lines (VBlockP vb, uint32_t *unconsumed_len, uint32_t *my_lines, uint32_t *her_lines);

// ZIP Stuff
extern void fastq_zip_initialize (void);
extern void fastq_zip_read_one_vb (VBlockP vb);
extern bool fastq_zip_dts_flag (void);
COMPRESSOR_CALLBACK (fastq_zip_seq);
COMPRESSOR_CALLBACK(fastq_zip_qual); // used by codec_longr_compress

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

