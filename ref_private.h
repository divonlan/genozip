// ------------------------------------------------------------------
//   ref_private.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
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

// contigs stuff
extern uint32_t ref_contigs_num_contigs (void);
typedef enum { WI_REF_CONTIG, WI_ZFILE_CHROM } GetWordIndexType;
extern WordIndex ref_contigs_get_word_index (const char *chrom_name, unsigned chrom_name_len, GetWordIndexType wi_type, bool soft_fail);
extern void ref_contigs_free (void);
extern const char *ref_contigs_get_chrom_snip (WordIndex chrom_index, const char **snip, uint32_t *snip_len);
extern WordIndex ref_alt_chroms_zip_get_alt_index (const char *chrom, unsigned chrom_len, GetWordIndexType where_is_alt, WordIndex fallback_index);
extern const RefContig *ref_contigs_get_contig (WordIndex chrom_index, bool soft_fail);
extern PosType ref_contigs_get_genome_size (void);
extern void ref_contigs_generate_data_if_denovo (void);
extern WordIndex ref_seg_get_alt_chrom (VBlockP vb);
extern void ref_contigs_compress (void);
extern WordIndex ref_contigs_get_by_accession_number (const char *ac, unsigned ac_len);

extern Buffer ranges; // param is RangesType

// used when compressing a reference and contigs buffer. note: ref_prepare_range_for_compress sets is_set.num_of_bits=0 
// for unused ranges. also, in make_ref mode, is_set is always 0, and ref is 0 for unused ranges. 
#define ref_is_range_used(r) ((r)->ref.num_of_bits && ((r)->is_set.num_of_bits || flag_make_reference))

#define ROUNDUP64(x) (((x) + 63) & ~(typeof(x))0x3f) // round up to the nearest 64
#define ROUNDDOWN64(x) ((x)      & ~(typeof(x))0x3f) // round down to the nearest 64
#endif
