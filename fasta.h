// ------------------------------------------------------------------
//   fasta.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FASTA_INCLUDED
#define FASTA_INCLUDED

#include "genozip.h"
#include "sections.h"

// Txtfile stuff
extern int32_t fasta_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);

// ZIP Stuff
COMPRESSOR_CALLBACK(fasta_zip_seq)
extern void ref_make_create_range (VBlockP vb);

// SEG Stuff
extern void fasta_seg_initialize();
extern void fasta_seg_finalize (VBlockP vb);
extern bool fasta_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern const char *fasta_seg_txt_line();

// PIZ Stuff
extern void fasta_piz_initialize (void);
extern bool fasta_piz_read_one_vb (VBlockP vb, ConstSectionListEntryP sl);
extern void fasta_reconstruct_vb(); // no parameter - implicit casting of VBlockP
extern bool fasta_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
extern bool fasta_piz_initialize_contig_grepped_out (VBlockP vb, bool does_vb_have_any_desc, bool last_desc_in_this_vb_matches_grep);
extern bool fasta_piz_is_grepped_out_due_to_regions (VBlockP vb, const char *line_start);

CONTAINER_FILTER_FUNC (fasta_piz_filter);

// VBlock stuff
extern void fasta_vb_release_vb();
extern void fasta_vb_destroy_vb();
extern unsigned fasta_vb_size (void);
extern unsigned fasta_vb_zip_dl_size (void);
extern void fasta_get_data_line (VBlockP vb_, uint32_t line_i, uint32_t *seq_data_start, uint32_t *seq_len);

#define FASTA_SPECIAL { fasta_piz_special_SEQ, fasta_piz_special_COMMENT, fasta_piz_special_DESC }
SPECIAL (FASTA, 0, SEQ, fasta_piz_special_SEQ);
SPECIAL (FASTA, 1, COMMENT, fasta_piz_special_COMMENT);
SPECIAL (FASTA, 2, DESC, fasta_piz_special_DESC);
#define NUM_FASTA_SPECIAL 3

#define FASTA_DICT_ID_ALIASES 

#define FASTA_LOCAL_GET_LINE_CALLBACKS  \
    { DT_FASTA, &dict_id_fields[FASTA_NONREF], fasta_zip_seq }, 

TXTHEADER_TRANSLATOR (txtheader_fa2phy);

#endif
