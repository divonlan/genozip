// ------------------------------------------------------------------
//   fasta.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FASTA_INCLUDED
#define FASTA_INCLUDED

#include "genozip.h"
#include "sections.h"

// ZIP Stuff
COMPRESSOR_CALLBACK(fasta_zip_get_start_len_line_i_seq)
extern void ref_make_create_range (VBlockP vb);

// SEG Stuff
extern void fasta_seg_initialize();
extern const char *fasta_seg_txt_line();

// PIZ Stuff
extern void fasta_piz_initialize (void);
extern bool fasta_piz_read_one_vb (VBlockP vb, SectionListEntryP sl);
extern void fasta_piz_reconstruct_vb(); // no parameter - implicit casting of VBlockP
extern bool fasta_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);

// VBlock stuff
extern void fast_vb_release_vb();
extern unsigned fast_vb_size (void);
extern unsigned fast_vb_zip_dl_size (void);


#define FASTA_SPECIAL { fasta_piz_special_SEQ, fasta_piz_special_COMMENT, fasta_piz_special_DESC }
SPECIAL (FASTA, 0, SEQ, fasta_piz_special_SEQ);
SPECIAL (FASTA, 1, COMMENT, fasta_piz_special_COMMENT);
SPECIAL (FASTA, 2, DESC, fasta_piz_special_DESC);
#define NUM_FASTA_SPECIAL 3

#define FASTA_DICT_ID_ALIASES 

#define FASTA_LOCAL_GET_LINE_CALLBACKS  \
    { DT_FASTA, &dict_id_fields[FASTA_SEQ],          fasta_zip_get_start_len_line_i_seq  }, 

#endif
