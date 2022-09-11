// ------------------------------------------------------------------
//   random_access.h
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#pragma once

#include <stdint.h>
#include "genozip.h"
#include "sections.h"

#define RA_MISSING_RA_MIN (MAX_POS + 1)
#define RA_MISSING_RA_MAX (-1)

// ZIP
extern void random_access_initialize(void);
extern void random_access_update_chrom (VBlockP vb, int ra_i, WordIndex chrom_node_index, STRp(chrom_name));
extern void random_access_update_pos (VBlockP vb, int ra_i, Did did_i_pos);
extern void random_access_increment_last_pos (VBlockP vb, int ra_i, PosType increment);
extern void random_access_update_last_pos (VBlockP vb, int ra_i, PosType last_pos);
extern void random_access_update_first_last_pos (VBlockP vb, int ra_i, WordIndex chrom_node_index, STRp (first_pos), STRp (last_pos));
extern void random_access_update_to_entire_chrom (VBlockP vb, int ra_i, PosType first_pos_of_chrom, PosType last_pos_of_chrom);
extern void random_access_merge_in_vb (VBlockP vb, int ra_i);
extern void random_access_finalize_entries (BufferP ra_buf);
extern Codec random_access_compress (ConstBufferP ra_buf, SectionType sec_type, Codec codec, int ra_i, rom msg);
extern void random_access_get_ra_info (VBIType vblock_i, WordIndex *chrom_index, PosType *min_pos, PosType *max_pos);

// PIZ
extern bool random_access_has_filter (void);
extern bool random_access_is_vb_included (VBIType vb_i);
extern uint32_t random_access_num_chroms_start_in_this_vb (VBIType vb_i);
extern void random_access_alloc_ra_buf (VBlockP vb, int ra_i, WordIndex chrom_node_index);
extern void random_access_load_ra_section (SectionType section_type, Did chrom_did_i, BufferP ra_buf, rom buf_name, rom show_index_msg);

// PIZ - FASTA
extern bool random_access_does_last_chrom_continue_in_next_vb (VBIType vb_i);
extern uint32_t random_access_verify_all_contigs_same_length (void);

#define RA_MSG_PRIM "Random-access index contents (result of --show-index)"
#define RA_MSG_LUFT "Luft random-access index contents (result of --show-index)"
#define RA_MSG_REF  "Reference random-access index contents (result of --show-ref-index)"
