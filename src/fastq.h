// ------------------------------------------------------------------
//   fastq.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "sections.h"

#pragma GENDICT_PREFIX FASTQ

// -----------------------------------------------------------------------------------------------------------
// Common contexts of FASTQ and SAM - these MUST be first in same order exactly in SAM/FASTQ for Deep to work.
// -----------------------------------------------------------------------------------------------------------
#pragma GENDICT FASTQ_CONTIG=DTYPE_FIELD=CONTIG     // must be first - did_i=0=CHROM. Note: not used in --deep, so nevermind that it conflicts with SAM's RNAME

#pragma GENDICT FASTQ_QNAME=DTYPE_FIELD=QNAME       // Q?NAME must immediately follow (up to v14 the dict_id was DESC)
#pragma GENDICT FASTQ_Q0NAME=DTYPE_1=Q0NAME         // MAX_QNAME_ITEMS fixed qname items must have a did_i directly after container's (MUST be the same dict_id as in sam.h)
#pragma GENDICT FASTQ_Q1NAME=DTYPE_1=Q1NAME 
#pragma GENDICT FASTQ_Q2NAME=DTYPE_1=Q2NAME
#pragma GENDICT FASTQ_Q3NAME=DTYPE_1=Q3NAME
#pragma GENDICT FASTQ_Q4NAME=DTYPE_1=Q4NAME
#pragma GENDICT FASTQ_Q5NAME=DTYPE_1=Q5NAME
#pragma GENDICT FASTQ_Q6NAME=DTYPE_1=Q6NAME 
#pragma GENDICT FASTQ_Q7NAME=DTYPE_1=Q7NAME 
#pragma GENDICT FASTQ_Q8NAME=DTYPE_1=Q8NAME 
#pragma GENDICT FASTQ_Q9NAME=DTYPE_1=Q9NAME 
#pragma GENDICT FASTQ_QANAME=DTYPE_1=QANAME     
#pragma GENDICT FASTQ_QBNAME=DTYPE_1=QBNAME         // if adding more Q*NAMEs - add to kraken.h and sam.h and update MAX_QNAME_ITEMS
#pragma GENDICT FASTQ_QmNAME=DTYPE_1=QmNAME         // QmNAME reserved for mate number (always the last dict_id in the container)

#pragma GENDICT FASTQ_QNAME2=DTYPE_FIELD=QNAME2     // QNAME2 is embedded in QNAME (QNAME2 items immediately follow)
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
#pragma GENDICT FASTQ_QANAME2=DTYPE_1=qANAME 
#pragma GENDICT FASTQ_QBNAME2=DTYPE_1=qBNAME 
#pragma GENDICT FASTQ_QmNAME2=DTYPE_1=qmNAME

#pragma GENDICT FASTQ_EXTRA=DTYPE_1=DESC            // until v14, the entire line1 was segged in FASTQ_DESC, since v15 is it used just for the "extra info" field in line1. dict_id remains "DESC" for backward compatibility

#pragma GENDICT FASTQ_AUX=DTYPE_FIELD=AUX

#pragma GENDICT FASTQ_SQBITMAP=DTYPE_FIELD=SQBITMAP
#pragma GENDICT FASTQ_NONREF=DTYPE_FIELD=NONREF
#pragma GENDICT FASTQ_NONREF_X=DTYPE_FIELD=NONREF_X
#pragma GENDICT FASTQ_GPOS=DTYPE_FIELD=GPOS
#pragma GENDICT FASTQ_GPOS_DELTA=DTYPE_FIELD=G0POS  // v15
#pragma GENDICT FASTQ_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT FASTQ_SEQMIS_A=DTYPE_FIELD=SEQMIS_A // v14: mismatch bases vs the reference, when ref=A
#pragma GENDICT FASTQ_SEQMIS_C=DTYPE_FIELD=SEQMIS_C
#pragma GENDICT FASTQ_SEQMIS_G=DTYPE_FIELD=SEQMIS_G
#pragma GENDICT FASTQ_SEQMIS_T=DTYPE_FIELD=SEQMIS_T

#pragma GENDICT FASTQ_QUAL=DTYPE_FIELD=QUAL 
#pragma GENDICT FASTQ_DOMQRUNS=DTYPE_FIELD=DOMQRUNS // these 3 must be right after FASTQ_QUAL. DOMQRUNS is also used by LONGR. For backwards compatability, we can never change its name.
#pragma GENDICT FASTQ_QUALMPLX=DTYPE_FIELD=QUALMPLX // v14.0.0. DOMQUAL alg: dom multiplexer 
#pragma GENDICT FASTQ_DIVRQUAL=DTYPE_FIELD=DIVRQUAL // v14.0.0. DOMQUAL alg: lines that don't have enough dom. NORMQ codec: lines that are of length other than segconf.sam_seq_len 

#pragma GENDICT FASTQ_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT FASTQ_BUDDY=DTYPE_FIELD=BUDDY       
#pragma GENDICT FASTQ_TAXID=DTYPE_FIELD=TAXID
#pragma GENDICT FASTQ_DEBUG_LINES=DTYPE_FIELD=DBGLINES  // used by --debug-lines

#pragma GENDICT FASTQ_DEEP=DTYPE_FIELD=DEEP         // v15: deep: SAM vb/line corresponding to this read  
#pragma GENDICT FASTQ_DEEP_DELTA=DTYPE_FIELD=D0EEP  // v15: in pair-2, if there is no delta vs pair-1

#pragma GENDICT FASTQ_E1L=DTYPE_FIELD=E1L
#pragma GENDICT FASTQ_E2L=DTYPE_FIELD=E2L

#pragma GENDICT FASTQ_LINE3=DTYPE_FIELD=LINE3  
#pragma GENDICT FASTQ_T0HIRD=DTYPE_1=t0NAME         // must be directly after FASTQ_LINE3
#pragma GENDICT FASTQ_T1HIRD=DTYPE_1=t1NAME 
#pragma GENDICT FASTQ_T2HIRD=DTYPE_1=t2NAME
#pragma GENDICT FASTQ_T3HIRD=DTYPE_1=t3NAME
#pragma GENDICT FASTQ_T4HIRD=DTYPE_1=t4NAME
#pragma GENDICT FASTQ_T5HIRD=DTYPE_1=t5NAME
#pragma GENDICT FASTQ_T6HIRD=DTYPE_1=t6NAME 
#pragma GENDICT FASTQ_T7HIRD=DTYPE_1=t7NAME 
#pragma GENDICT FASTQ_T8HIRD=DTYPE_1=t8NAME 
#pragma GENDICT FASTQ_T9HIRD=DTYPE_1=t9NAME 
#pragma GENDICT FASTQ_TAHIRD=DTYPE_1=tANAME 
#pragma GENDICT FASTQ_TBHIRD=DTYPE_1=tBNAME 
#pragma GENDICT FASTQ_TmHIRD=DTYPE_1=tmNAME 

#pragma GENDICT FASTQ_AUX_LENGTH=DTYPE_2=length 

// -----------------------------------------------------------------------------------------------------------
// End of common contexts of FASTQ and SAM
// -----------------------------------------------------------------------------------------------------------

// Txtfile stuff
extern int32_t fastq_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern bool fastq_txtfile_have_enough_lines (VBlockP vb, uint32_t *unconsumed_len, uint32_t *my_lines, VBIType *pair_vb_i, uint32_t *pair_lines, uint32_t *pair_txt_data_len);
extern bool is_fastq (STRp(header), bool *need_more);
extern bool is_fastq_pair_2 (VBlockP vb);
extern void fastq_zip_set_txt_header_flags (struct FlagsTxtHeader *f);

// ZIP Stuff
extern void fastq_zip_initialize (void);
extern void fastq_segconf_set_r1_or_r2 (void);
extern void fastq_zip_finalize (bool is_last_user_txt_file);
extern void fastq_zip_init_vb (VBlockP vb);
extern void fastq_zip_after_compute (VBlockP vb);
extern bool fastq_zip_use_pair_assisted (DictId dict_id, SectionType st);
extern bool fastq_zip_use_pair_identical (DictId dict_id);

COMPRESSOR_CALLBACK (fastq_zip_seq);
COMPRESSOR_CALLBACK(fastq_zip_qual); // used by codec_longr_compress

// SEG Stuff
extern void fastq_seg_initialize();
extern void fastq_seg_finalize();
extern bool fastq_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern rom fastq_seg_txt_line();
extern void fastq_seg_pair2_gpos (VBlockP vb, PosType64 pair1_pos, PosType64 pair2_gpos);

// PIZ Stuff
extern void fastq_piz_process_recon (VBlockP vb);
extern void fastq_piz_before_read (VBlockP vb);
extern bool fastq_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header, uint32_t *txt_data_so_far_single_0_increment);
CONTAINER_FILTER_FUNC (fastq_piz_filter);
extern IS_SKIP (fastq_piz_is_skip_section);
extern void fastq_recon_aligned_SEQ (VBlockP vb, STRp(seq_len_str), ReconType reconstruct);
extern CONTAINER_CALLBACK (fastq_piz_container_cb);

// VBlock stuff
extern void fastq_vb_release_vb();
extern unsigned fastq_vb_size (DataType dt);
extern unsigned fastq_vb_zip_dl_size (void);
extern void fastq_reset_line (VBlockP vb);

// file pairing (--pair) stuff
extern void fastq_read_pair_1_data (VBlockP vb, VBIType pair_vb_i);

// FASTQ-specific fields in genozip header
extern void fastq_zip_genozip_header (SectionHeaderGenozipHeaderP header);
extern void fastq_piz_genozip_header (ConstSectionHeaderGenozipHeaderP header);

#define FASTQ_LOCAL_GET_LINE_CALLBACKS           \
    { DT_FASTQ, _FASTQ_QUAL,   fastq_zip_qual }, 
 // { DT_FASTQ, _FASTQ_NONREF, fastq_zip_seq  },  this callback is called directly from codec_longr, it is not added to the list, since it should NOT be called in other codecs

#define FASTQ_DICT_ID_ALIASES                               \
    /*           type        alias        maps to      */   \
    { DT_FASTQ,  ALIAS_DICT, _FASTQ_E1L, _FASTQ_E2L     },  /* each of E1L and E2L *may* have different b250s if lines differ in their EOL, but normally they are both all-the-same, so only the dict of E2L survives, and no b250 sections */

typedef enum { FQ_COMP_R1, FQ_COMP_R2 } FastqComponentType;
#define FASTQ_COMP_NAMES { "FQR1", "FQR2" }
 
#define FASTQ_SPECIAL { fastq_special_unaligned_SEQ, fastq_special_PAIR2_GPOS, fastq_special_mate_lookup, \
                        fastq_special_set_deep, fastq_special_deep_copy_QNAME, fastq_special_deep_copy_SEQ, fastq_special_deep_copy_QUAL,\
                        fastq_special_backspace, fastq_special_copy_line1 }

SPECIAL (FASTQ, 0,  unaligned_SEQ,   fastq_special_unaligned_SEQ);    // v14
SPECIAL (FASTQ, 1,  PAIR2_GPOS,      fastq_special_PAIR2_GPOS);       // v14
SPECIAL (FASTQ, 2,  mate_lookup,     fastq_special_mate_lookup);      // v14
SPECIAL (FASTQ, 3,  set_deep,        fastq_special_set_deep);         // v15
SPECIAL (FASTQ, 4,  deep_copy_QNAME, fastq_special_deep_copy_QNAME);  // v15
SPECIAL (FASTQ, 5,  deep_copy_SEQ,   fastq_special_deep_copy_SEQ);    // v15
SPECIAL (FASTQ, 6,  deep_copy_QUAL,  fastq_special_deep_copy_QUAL);   // v15
SPECIAL (FASTQ, 7,  backspace,       fastq_special_backspace);        // v15
SPECIAL (FASTQ, 8,  copy_line1,      fastq_special_copy_line1);       // v15
#define NUM_FASTQ_SPECIAL 9
