// ------------------------------------------------------------------
//   chain.h
//   Copyright (C) 2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef CHAIN_INCLUDED
#define CHAIN_INCLUDED

#include "genozip.h"

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
