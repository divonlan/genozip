// ------------------------------------------------------------------
//   random_access.h
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef RANDOM_ACCESS_INCLUDED
#define RANDOM_ACCESS_INCLUDED

#include <stdint.h>
#include "genozip.h"
#include "sections.h"

#define RA_MISSING_RA_MIN (MAX_POS + 1)
#define RA_MISSING_RA_MAX (-1)

extern void random_access_initialize(void);
extern void random_access_update_chrom (VBlockP vb, WordIndex chrom_node_index, const char *chrom_name, unsigned chrom_name_len);
extern void random_access_update_pos (VBlockP vb, DidIType did_i_pos);
extern void random_access_increment_last_pos (VBlockP vb, PosType increment);
extern void random_access_update_last_pos (VBlockP vb, PosType last_pos);
extern void random_access_update_to_entire_chrom (VBlockP vb, PosType first_pos_of_chrom, PosType last_pos_of_chrom);
extern void random_access_merge_in_vb (VBlockP vb);
extern void random_access_finalize_entries (BufferP ra_buf);
extern void BGEN_random_access (Buffer *ra_buf);
extern void random_access_show_index (ConstBufferP ra_buf, bool from_zip, const char *msg);
extern bool random_access_is_vb_included (uint32_t vb_i, BufferP region_ra_intersection_matrix);
extern int32_t random_access_get_last_included_vb_i (void);
extern void random_access_get_first_chrom_of_vb (VBlockP vb, PosType *first_pos, PosType *last_pos);
extern bool random_access_does_last_chrom_continue_in_next_vb (uint32_t vb_i);
extern uint32_t random_access_num_chroms_start_in_this_vb (uint32_t vb_i);
extern void random_access_alloc_ra_buf (VBlockP vb, WordIndex chrom_node_index);
extern void random_access_pos_of_chrom (WordIndex chrom_word_index, PosType *min_pos, PosType *max_pos);
extern void random_access_get_ra_info (uint32_t vblock_i, WordIndex *chrom_index, PosType *min_pos, PosType *max_pos);
extern void random_access_load_ra_section (SectionType section_type, BufferP ra_buf, const char *buf_name, const char *show_index_msg);
extern void random_access_compress (Buffer *ra_buf, SectionType sec_type, const char *msg);
extern uint32_t random_access_verify_all_contigs_same_length (void);

#endif