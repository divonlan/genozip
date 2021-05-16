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

// ZIP
extern void random_access_initialize(void);
extern void random_access_update_chrom (VBlockP vb, Coords dc, WordIndex chrom_node_index, const char *chrom_name, unsigned chrom_name_len);
extern void random_access_update_pos (VBlockP vb, Coords dc, DidIType did_i_pos);
extern void random_access_increment_last_pos (VBlockP vb, Coords dc, PosType increment);
extern void random_access_update_last_pos (VBlockP vb, Coords dc, PosType last_pos);
extern void random_access_update_to_entire_chrom (VBlockP vb, Coords dc, PosType first_pos_of_chrom, PosType last_pos_of_chrom);
extern void random_access_merge_in_vb (VBlockP vb, Coords dc);
extern void random_access_finalize_entries (BufferP ra_buf);
extern void random_access_compress (Buffer *ra_buf, SectionType sec_type, Coords dc, const char *msg);

// PIZ
extern bool random_access_is_vb_included (uint32_t vb_i);
extern int32_t random_access_get_last_included_vb_i (void);
extern bool random_access_does_last_chrom_continue_in_next_vb (uint32_t vb_i);
extern uint32_t random_access_num_chroms_start_in_this_vb (uint32_t vb_i);
extern void random_access_alloc_ra_buf (VBlockP vb, Coords dc, WordIndex chrom_node_index);
extern void random_access_get_ra_info (uint32_t vblock_i, WordIndex *chrom_index, PosType *min_pos, PosType *max_pos);
extern void random_access_load_ra_section (SectionType section_type, BufferP ra_buf, const char *buf_name, const char *show_index_msg);
extern uint32_t random_access_verify_all_contigs_same_length (void);

// misc
extern void BGEN_random_access (Buffer *ra_buf);
extern void random_access_show_index (ConstBufferP ra_buf, bool from_zip, const char *msg);

#endif