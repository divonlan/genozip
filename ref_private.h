// ------------------------------------------------------------------
//   ref_private.h
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REF_PRIVATE_INCLUDED
#define REF_PRIVATE_INCLUDED

#include "genozip.h"
#include "reference.h"
#include "sections.h"

extern void ref_make_prepare_range_for_compress (VBlockP vb);

// lock stuff
extern void ref_lock_initialize_loaded_genome (void);
extern void ref_lock_initialize_denovo_genome (void);
extern void ref_lock_free (void);
extern void ref_lock_destroy (void);

// contigs stuff
extern uint32_t ref_contigs_num_contigs (void);
extern void ref_contigs_free (void);
extern void ref_contigs_destroy (void);
extern const char *ref_contigs_get_chrom_snip (WordIndex chrom_index, const char **snip, uint32_t *snip_len);
extern const RefContig *ref_contigs_get_contig (WordIndex chrom_index, bool soft_fail);
extern PosType ref_contigs_get_genome_nbases (void);
extern void ref_contigs_generate_data_if_denovo (void);
extern WordIndex ref_seg_get_alt_chrom (VBlockP vb);
extern void ref_contigs_compress (void);
extern WordIndex ref_contigs_get_by_accession_number (const char *ac, unsigned ac_len);
extern WordIndex ref_alt_chroms_zip_get_alt_index (const char *chrom, unsigned chrom_len, GetWordIndexType where_is_alt, WordIndex fallback_index);

extern void sam_seg_verify_pos (VBlockP vb, PosType this_pos);

extern Buffer ranges; // param is RangesType

#define IS_REF_INTERNAL(f) (((f)->data_type == DT_SAM || (f)->data_type == DT_BAM) && (f)->z_flags.dts_ref_internal)

#define ROUNDUP32(x) (((x) + 32) & ~(typeof(x))0x1f) // round up to the nearest 32
#define ROUNDDOWN32(x) ((x)      & ~(typeof(x))0x1f) // round down to the nearest 32

#define ROUNDUP64(x) (((x) + 63) & ~(typeof(x))0x3f) // round up to the nearest 64
#define ROUNDDOWN64(x) ((x)      & ~(typeof(x))0x3f) // round down to the nearest 64

#endif
