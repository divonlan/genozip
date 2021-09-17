// ------------------------------------------------------------------
//   chain.h
//   Copyright (C) 2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef CHAIN_INCLUDED
#define CHAIN_INCLUDED

#include "genozip.h"

#pragma GENDICT_PREFIX CHAIN
#pragma GENDICT CHAIN_NAMELUFT=DTYPE_FIELD=NaMELUFT
#pragma GENDICT CHAIN_STRNDLUFT=DTYPE_FIELD=SrNDLUFT
#pragma GENDICT CHAIN_STARTLUFT=DTYPE_FIELD=StRTLUFT
#pragma GENDICT CHAIN_ENDLUFT=DTYPE_FIELD=EnDLUFT
#pragma GENDICT CHAIN_SIZELUFT=DTYPE_FIELD=SiZELUFT
#pragma GENDICT CHAIN_NAMEPRIM=DTYPE_FIELD=NAMEPRIM
#pragma GENDICT CHAIN_STRNDPRIM=DTYPE_FIELD=SRNDPRIM
#pragma GENDICT CHAIN_STARTPRIM=DTYPE_FIELD=STRTPRIM
#pragma GENDICT CHAIN_ENDPRIM=DTYPE_FIELD=ENDPRIM
#pragma GENDICT CHAIN_SIZEPRIM=DTYPE_FIELD=SIZEPRIM
#pragma GENDICT CHAIN_CHAIN=DTYPE_FIELD=CHAIN
#pragma GENDICT CHAIN_SCORE=DTYPE_FIELD=SCORE
#pragma GENDICT CHAIN_ID=DTYPE_FIELD=ID
#pragma GENDICT CHAIN_VERIFIED=DTYPE_FIELD=VERIFIED
#pragma GENDICT CHAIN_SET=DTYPE_FIELD=SET
#pragma GENDICT CHAIN_SIZE=DTYPE_FIELD=SIZE
#pragma GENDICT CHAIN_GAPS=DTYPE_FIELD=GAPS
#pragma GENDICT CHAIN_EOL=DTYPE_FIELD=EOL
#pragma GENDICT CHAIN_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT CHAIN_SEP=DTYPE_FIELD=SEP

// vblock stuff
extern unsigned chain_vb_size (DataType dt);
extern void chain_vb_release_vb();

// zip of a chain file - txtfile
extern void chain_zip_initialize (void);
extern int32_t chain_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i /* in/out */);

// seg of a chain file
extern void chain_seg_initialize (VBlockP vb);
extern void chain_seg_finalize (VBlockP vb);
extern const char *chain_seg_txt_line (VBlockP vb, const char *line, uint32_t remaining_txt_len, bool *has_13);
extern bool chain_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern bool chain_zip_dts_flag (void);

// piz of a chain file
extern bool chain_piz_initialize (void);
extern CONTAINER_FILTER_FUNC (chain_piz_filter);

// using the chain data in genozip --chain
extern void chain_load (void);
extern void chain_destroy (void);
extern const char *chain_get_luft_contig (uint32_t contig_i, PosType *length);
extern void chain_copy_contigs_to_z_file (DidIType luft_contig_did_i);
extern WordIndex chain_get_prim_contig_index (const char *contig, unsigned contig_len, PosType LN);
extern uint64_t chain_get_num_prim_contigs (void);
extern void chain_append_all_luft_contig_index (const char *prim_contig_name, unsigned prim_contig_name_len, PosType LN, Buffer *luft_contigs);
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
#define chain_is_loaded ((bool)chain_filename || flag.reference == REF_MAKE_CHAIN)

#endif
