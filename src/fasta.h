// ------------------------------------------------------------------
//   fasta.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "sections.h"

#pragma GENDICT_PREFIX FASTA

// -----------------------------------------------------------------------------------------------------------
// Common contexts of FASTA and GFF - these MUST be first; in same order; same dict_id
// -----------------------------------------------------------------------------------------------------------

#pragma GENDICT FASTA_CONTIG=DTYPE_FIELD=CONTIG

#pragma GENDICT FASTA_QNAME=DTYPE_FIELD=QNAME // MAX_QNAME_ITEMS - same did as SAM/FASTQ/GFF - used for testing if file is actually a QUAL-less FASTQ
#pragma GENDICT FASTA_Q0NAME=DTYPE_1=Q0NAME // must have a did_i directly after QNAME
#pragma GENDICT FASTA_Q1NAME=DTYPE_1=Q1NAME 
#pragma GENDICT FASTA_Q2NAME=DTYPE_1=Q2NAME
#pragma GENDICT FASTA_Q3NAME=DTYPE_1=Q3NAME
#pragma GENDICT FASTA_Q4NAME=DTYPE_1=Q4NAME
#pragma GENDICT FASTA_Q5NAME=DTYPE_1=Q5NAME
#pragma GENDICT FASTA_Q6NAME=DTYPE_1=Q6NAME 
#pragma GENDICT FASTA_Q7NAME=DTYPE_1=Q7NAME 
#pragma GENDICT FASTA_Q8NAME=DTYPE_1=Q8NAME 
#pragma GENDICT FASTA_Q9NAME=DTYPE_1=Q9NAME 
#pragma GENDICT FASTA_QANAME=DTYPE_1=QANAME 
#pragma GENDICT FASTA_QBNAME=DTYPE_1=QBNAME 
#pragma GENDICT FASTA_QCNAME=DTYPE_1=QCNAME 
#pragma GENDICT FASTA_QDNAME=DTYPE_1=QDNAME 
#pragma GENDICT FASTA_QENAME=DTYPE_1=QENAME // if adding more Q*NAMEs - add to fastq.h too, and update MAX_QNAME_ITEMS
#pragma GENDICT FASTA_QmNAME=DTYPE_1=QmNAME // QmNAME reserved for mate number (always the last dict_id in the container)

#pragma GENDICT FASTA_LINEMETA=DTYPE_FIELD=LINEMETA
#pragma GENDICT FASTA_EOL=DTYPE_FIELD=EOL
#pragma GENDICT FASTA_DESC=DTYPE_FIELD=DESC
#pragma GENDICT FASTA_COMMENT=DTYPE_FIELD=COMMENT
#pragma GENDICT FASTA_NONREF=DTYPE_FIELD=NONREF 
#pragma GENDICT FASTA_NONREF_X=DTYPE_FIELD=NONREF_X
#pragma GENDICT FASTA_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT FASTA_TAXID=DTYPE_FIELD=TAXID
#pragma GENDICT FASTA_DEBUG_LINES=DTYPE_FIELD=DBGLINES      // used by --debug-lines

// -----------------------------------------------------------------------------------------------------------
// End of common contexts of FASTA and GFF
// -----------------------------------------------------------------------------------------------------------

// Txtfile stuff
extern int32_t fasta_unconsumed (VBlockP vb, uint32_t first_i);
extern bool is_fasta (STRp(header), bool *need_more);

// ZIP Stuff
extern void fasta_zip_initialize (void);
COMPRESSOR_CALLBACK(fasta_zip_seq);
extern void fasta_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeaderP vb_header);
extern void fasta_zip_after_compute (VBlockP vb);

// SEG Stuff
extern void fasta_seg_initialize(VBlockP vb);
extern void fasta_segconf_finalize (VBlockP vb);
extern void fasta_seg_finalize (VBlockP vb);
extern bool fasta_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern rom fasta_seg_txt_line (VBlockP vb, rom line_start, uint32_t remaining, bool *has_13);
extern bool fasta_seg_is_big (ConstVBlockP vb, DictId dict_id, DictId st_dict_id);

// PIZ Stuff
extern bool fasta_piz_initialize (CompIType comp_i);
extern void fasta_piz_process_recon (VBlockP vb);
extern bool fasta_piz_is_vb_needed (VBIType vb_i);
extern bool fasta_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header);
extern void fasta_reconstruct_vb(); // no parameter - implicit casting of VBlockP
extern IS_SKIP (fasta_piz_is_skip_section);
extern bool fastq_piz_get_r2_is_forward (VBlockP vb);
extern bool fastq_piz_get_interleaved_r2_is_forward (VBlockP vb);

// VBlock stuff
extern unsigned fasta_vb_size (DataType dt);
extern unsigned fasta_vb_zip_dl_size (void);
extern void fasta_get_data_line (VBlockP vb_, uint32_t line_i, uint32_t *seq_data_start, uint32_t *seq_len);

// SPECIALs
SPECIAL (FASTA, 0, SEQ,     fasta_piz_special_SEQ);
SPECIAL (FASTA, 1, COMMENT, fasta_piz_special_COMMENT);
SPECIAL (FASTA, 2, DESC,    fasta_piz_special_DESC);

#define FASTA_LOCAL_GET_LINE_CALLBACKS  \
    { DT_FASTA, _FASTA_NONREF, fasta_zip_seq }, 

CONTAINER_FILTER_FUNC (fasta_piz_filter);
