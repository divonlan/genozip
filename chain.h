// ------------------------------------------------------------------
//   chain.h
//   Copyright (C) 2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef CHAIN_INCLUDED
#define CHAIN_INCLUDED

#include "genozip.h"

#define _CHAIN_NAMELUFT   DICT_ID_MAKEF_L ("NaMELUFT")
#define _CHAIN_STRNDLUFT  DICT_ID_MAKEF_L ("SrNDLUFT")
#define _CHAIN_STARTLUFT  DICT_ID_MAKEF_L ("StRTLUFT")
#define _CHAIN_ENDLUFT    DICT_ID_MAKEF_7 ("EnDLUFT")
#define _CHAIN_SIZELUFT   DICT_ID_MAKEF_L ("SiZELUFT")
#define _CHAIN_NAMEPRIM   DICT_ID_MAKEF_L ("NAMEPRIM")
#define _CHAIN_STRNDPRIM  DICT_ID_MAKEF_L ("SRNDPRIM")
#define _CHAIN_STARTPRIM  DICT_ID_MAKEF_L ("STRTPRIM")
#define _CHAIN_ENDPRIM    DICT_ID_MAKEF_7 ("ENDPRIM")
#define _CHAIN_SIZEPRIM   DICT_ID_MAKEF_L ("SIZEPRIM")
#define _CHAIN_CHAIN      DICT_ID_MAKEF_5 ("CHAIN")
#define _CHAIN_SCORE      DICT_ID_MAKEF_5 ("SCORE")
#define _CHAIN_ID         DICT_ID_MAKEF_2 ("ID")
#define _CHAIN_VERIFIED   DICT_ID_MAKEF_L ("VERIFIED")
#define _CHAIN_SET        DICT_ID_MAKEF_3 ("SET")
#define _CHAIN_SIZE       DICT_ID_MAKEF_4 ("SIZE")
#define _CHAIN_GAPS       DICT_ID_MAKEF_4 ("GAPS")
#define _CHAIN_EOL        DICT_ID_MAKEF_3 ("EOL")
#define _CHAIN_TOPLEVEL   DICT_ID_MAKEF_L (TOPLEVEL)
#define _CHAIN_SEP        DICT_ID_MAKEF_3 ("SEP")

typedef enum { CHAIN_NAMELUFT, CHAIN_STRNDLUFT, CHAIN_STARTLUFT, CHAIN_ENDLUFT, CHAIN_SIZELUFT,
               CHAIN_NAMEPRIM, CHAIN_STRNDPRIM, CHAIN_STARTPRIM, CHAIN_ENDPRIM, CHAIN_SIZEPRIM,
               CHAIN_CHAIN, CHAIN_SCORE, CHAIN_ID, CHAIN_VERIFIED, CHAIN_SET, 
               CHAIN_SIZE, CHAIN_GAPS, CHAIN_EOL, CHAIN_TOPLEVEL, CHAIN_SEP, NUM_CHAIN_FIELDS } ChainFields;

#define CHAIN_MAPPING { V(CHAIN_NAMELUFT), V(CHAIN_STRNDLUFT), V(CHAIN_STARTLUFT), V(CHAIN_ENDLUFT), V(CHAIN_SIZELUFT), V(CHAIN_NAMEPRIM), V(CHAIN_STRNDPRIM), V(CHAIN_STARTPRIM), V(CHAIN_ENDPRIM), V(CHAIN_SIZEPRIM), V(CHAIN_CHAIN), V(CHAIN_SCORE), V(CHAIN_ID), V(CHAIN_VERIFIED), V(CHAIN_SET), V(CHAIN_SIZE), V(CHAIN_GAPS), V(CHAIN_EOL), V(CHAIN_TOPLEVEL), V(CHAIN_SEP) }

// zip of a chain file - txtfile
extern void chain_zip_initialize (void);
extern int32_t chain_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i /* in/out */);

// seg of a chain file
extern void chain_seg_initialize (VBlockP vb);
extern void chain_seg_finalize (VBlockP vb);
extern const char *chain_seg_txt_line (VBlockP vb, const char *line, uint32_t remaining_txt_len, bool *has_13);
extern bool chain_seg_is_small (ConstVBlockP vb, DictId dict_id);

// piz of a chain file
extern bool chain_piz_initialize (void);
extern CONTAINER_FILTER_FUNC (chain_piz_filter);

// using the chain data in genozip --chain
extern void chain_load (void);
extern void chain_destroy (void);
extern const char *chain_get_luft_contig (uint32_t contig_i, PosType *length);
extern void chain_copy_contigs_to_z_file (DidIType luft_contig_did_i);
extern WordIndex chain_get_prim_contig_index_by_name (const char *contig, unsigned contig_len, bool recursive);
extern uint64_t chain_get_num_prim_contigs (void);
extern void chain_append_all_luft_contig_index (const char *prim_contig_name, unsigned prim_contig_name_len, Buffer *luft_contigs);
extern bool chain_get_liftover_coords (WordIndex prim_contig_index,  PosType prim_1pos, 
                                       WordIndex *luft_contig_index, PosType *luft_1pos, bool *is_xstrand, uint32_t *aln_i);
extern PosType chain_get_aln_last_pos (uint32_t aln_i);
extern PosType chain_get_aln_gap_after (uint32_t aln_i);

#define CHAIN_SPECIAL { chain_piz_special_BACKSPACE, chain_piz_special_ENDLUFT, chain_piz_special_SIZE }
SPECIAL (CHAIN, 0, BACKSPACE, chain_piz_special_BACKSPACE);
SPECIAL (CHAIN, 1, ENDLUFT,    chain_piz_special_ENDLUFT);
SPECIAL (CHAIN, 2, SIZE,      chain_piz_special_SIZE);
#define NUM_CHAIN_SPECIAL 3

extern char *chain_filename; // global
#define chain_is_loaded ((bool)chain_filename)

#endif
