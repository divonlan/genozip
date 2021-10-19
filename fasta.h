// ------------------------------------------------------------------
//   fasta.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef FASTA_INCLUDED
#define FASTA_INCLUDED

#include "genozip.h"
#include "sections.h"

#pragma GENDICT_PREFIX FASTA
#pragma GENDICT FASTA_CONTIG=DTYPE_FIELD=CONTIG
#pragma GENDICT FASTA_LINEMETA=DTYPE_FIELD=LINEMETA
#pragma GENDICT FASTA_EOL=DTYPE_FIELD=EOL
#pragma GENDICT FASTA_DESC=DTYPE_FIELD=DESC
#pragma GENDICT FASTA_COMMENT=DTYPE_FIELD=COMMENT
#pragma GENDICT FASTA_NONREF=DTYPE_FIELD=NONREF 
#pragma GENDICT FASTA_NONREF_X=DTYPE_FIELD=NONREF_X
#pragma GENDICT FASTA_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT FASTA_TAXID=DTYPE_FIELD=TAXID

// Txtfile stuff
extern int32_t fasta_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);

// ZIP Stuff
extern void fasta_zip_initialize (void);
COMPRESSOR_CALLBACK(fasta_zip_seq)
extern void vcf_zip_after_compute (VBlockP vb);

// SEG Stuff
extern void fasta_seg_initialize(VBlockP vb);
extern void fasta_seg_finalize (VBlockP vb);
extern bool fasta_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern const char *fasta_seg_txt_line();

// PIZ Stuff
extern bool fasta_piz_initialize (void);
extern bool fasta_piz_read_one_vb (VBlockP vb, Section sl);
extern void fasta_reconstruct_vb(); // no parameter - implicit casting of VBlockP
extern bool fasta_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
extern bool fasta_piz_initialize_contig_grepped_out (VBlockP vb, bool does_vb_have_any_desc, bool last_desc_in_this_vb_matches_grep);

// VBlock stuff
extern void fasta_vb_release_vb();
extern void fasta_vb_destroy_vb();
extern unsigned fasta_vb_size (DataType dt);
extern unsigned fasta_vb_zip_dl_size (void);
extern void fasta_get_data_line (VBlockP vb_, uint32_t line_i, uint32_t *seq_data_start, uint32_t *seq_len);

#define FASTA_SPECIAL { fasta_piz_special_SEQ, fasta_piz_special_COMMENT, fasta_piz_special_DESC }
SPECIAL (FASTA, 0, SEQ, fasta_piz_special_SEQ);
SPECIAL (FASTA, 1, COMMENT, fasta_piz_special_COMMENT);
SPECIAL (FASTA, 2, DESC, fasta_piz_special_DESC);
#define NUM_FASTA_SPECIAL 3

#define FASTA_LOCAL_GET_LINE_CALLBACKS  \
    { DT_FASTA, _FASTA_NONREF, fasta_zip_seq }, 

// Multifasta to Phylip translation stuff

// Important: Numbers (and order) of translators cannot be changed, as they are part of the file format
// (they are included in the TOPLEVEL container)
// translator numbers must start from 1 - 0 is reserved for "none"
TRANSLATOR (FASTA, PHYLIP, 1, EOL, fasta_piz_fa2phy_EOL)   // drop EOL after DESC

#define NUM_FASTA_TRANS 2 // including "none"
#define FASTA_TRANSLATORS { NULL /* none */, fasta_piz_fa2phy_EOL }

TXTHEADER_TRANSLATOR (txtheader_fa2phy);
CONTAINER_FILTER_FUNC (fasta_piz_filter);

#endif
