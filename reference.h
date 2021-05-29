// ------------------------------------------------------------------
//   reference.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REFERENCE_INCLUDED
#define REFERENCE_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "digest.h"
#include "bit_array.h"
#include "flags.h"

// reference sequences - one per range of 1MB. ranges (chrom, pos) are mapped here with a hash function. In the rare case two unrelated ranges
// are mapped to the same entry - only the first range will have a reference, and the other ones will not. this will hurt the compression ratio,
// but not the correctness.
// Thread safety: "ref" is set atomically as the last set in initialization. If its set, then it is correct and will never change. If it appears to be
// not set yet, it is necesarry to lock the mutex and test again, and initialize if still not set.

typedef struct Range {
    uint32_t range_id;           // index of range within Buffer ranges
    BitArray ref;                // actual reference data - 2-bit array
    BitArray is_set;             // a 1-bit array - SEG: a pos is set if seg set this reference PIZ: is set if SEC_REF_IS_SET said so
    const char *chrom_name;
    unsigned chrom_name_len;
    WordIndex chrom;             // index to the chrom of the external reference
    uint32_t range_i;            // range ordinal number within contig
    PosType first_pos, last_pos; // the range that includes all locii (note: in ZIP-INTERNAL it might include unset locii too)
    PosType gpos;                // position of this range in the "global position" 
    uint32_t copied_first_index, copied_len; // ZIP with REF_EXT_STORE: the subset of this range that was copied directly from the fasta file and doesn't need to be compressed
} Range;

#define ref_size(r) ((r) ? ((r)->last_pos - (r)->first_pos + 1) : 0)

#define REF_NUM_DENOVO_RANGES (1 << 20)
#define REF_NUM_DENOVO_SITES_PER_RANGE (1 << 20) // 1 Mbp

// ZIP ONLY: access range_i and index within range, for ranges configured for ZIP
#define range_i2pos(range_i) ((PosType)(range_i) * REF_NUM_DENOVO_SITES_PER_RANGE)
#define pos2range_i(pos)   ((uint32_t)((pos) / REF_NUM_DENOVO_SITES_PER_RANGE))
#define pos2range_idx(pos) ((uint32_t)((pos) % REF_NUM_DENOVO_SITES_PER_RANGE))

// locks
typedef struct { int32_t first_mutex, last_mutex; } RefLock;
#define REFLOCK_NONE ((RefLock){-1,-1})

extern RefLock ref_lock (PosType gpos_start, uint32_t seq_len);
extern RefLock ref_unlock (RefLock lock);
extern RefLock ref_lock_range (int32_t range_id);

typedef enum { RT_NONE,      // value of ranges.param if ranges is unallocated
               RT_MAKE_REF,  // used in --make-ref one range per vb of fasta reference file - ranges in order of the fasta file
               RT_DENOVO,    // used in SAM with REF_INTERNAL - an large array of Range's, hashed by the chrom and pos
               RT_LOADED,    // one Range per chrom (contig), overlayed on genome
               RT_CACHED     // same as RT_LOADED, but data is mmap'ed to a cache file
             } RangesType;
#define ranges_type ranges.param
              
extern void ref_initialize_ranges (RangesType type);
extern void ref_compress_ref (void);
extern void ref_load_external_reference (bool display);
extern void ref_load_stored_reference (void);
extern bool ref_is_reference_loaded (void);
extern void ref_set_reference (const char *filename, ReferenceType ref_type, bool is_explicit);
extern void ref_set_ref_file_info (Digest md5, const char *fasta_name);
extern void ref_unload_reference (void);
extern void ref_destroy_reference (bool destroy_only_if_not_mmap);
extern MemStats ref_memory_consumption (void);
extern const Range *ref_piz_get_range (VBlockP vb, PosType first_pos_needed, uint32_t num_nucleotides_needed);
extern void ref_consume_ref_fasta_global_area (void);
extern Range *ref_seg_get_locked_range (VBlockP vb, WordIndex chrom, PosType pos, uint32_t seq_len, const char *field /* used for ASSSEG */, RefLock *lock);
extern const char *ref_get_cram_ref (void);
extern void ref_make_ref_init (void);
extern void ref_generate_reverse_complement_genome (void);
extern Range *ref_get_range_by_chrom (WordIndex chrom, const char **chrom_name);

// cache stuff
extern bool ref_mmap_cached_reference (void);
extern void ref_create_cache_in_background (void);
extern void ref_create_cache_join (void);
extern void ref_remove_cache (void);

// contigs stuff
typedef enum { WI_REF_CONTIG, WI_ZFILE_CHROM } GetWordIndexType;
extern WordIndex ref_contigs_get_word_index (const char *chrom_name, unsigned chrom_name_len, GetWordIndexType wi_type, bool soft_fail);

extern void ref_contigs_get (ConstBufferP *out_contig_dict, ConstBufferP *out_contigs);
extern uint32_t ref_num_loaded_contigs (void);
extern PosType ref_contigs_get_contig_length (WordIndex chrom_index, const char *chrom_name, unsigned chrom_name_len, bool enforce);
extern WordIndex ref_contigs_ref_chrom_from_header_chrom (const char *chrom_name, unsigned chrom_name_len, PosType *last_pos, WordIndex header_chrom);
extern void ref_contigs_sort_chroms (void);
extern void ref_contigs_load_contigs (void);

typedef void (*RefContigsIteratorCallback)(const char *chrom_name, unsigned chrom_name_len, PosType last_pos, void *callback_param);
extern void ref_contigs_iterate (RefContigsIteratorCallback callback, void *callback_param);
extern WordIndex ref_contig_get_by_gpos (PosType gpos);
extern WordIndex ref_contig_get_by_chrom (VBlockP vb, WordIndex txt_chrom_index);

// alt chroms stuff
extern void ref_alt_chroms_load (void);
extern void ref_alt_chroms_compress (void);
extern WordIndex ref_alt_chroms_zip_get_alt_index (const char *chrom, unsigned chrom_len, GetWordIndexType where_is_alt, WordIndex fallback_index);

// gets the index of the matching chrom in the reference - either its the chrom itself, or one with an alternative name
// eg 'chr22' instead of '22'
#define ref_alt_get_final_index(chrom_index) \
    (buf_is_alloc (&z_file->alt_chrom_map) ? *ENT (WordIndex, z_file->alt_chrom_map, (chrom_index)) : (chrom_index))

// note that the following work on idx and not pos! (idx is the index within the range)
#define ref_set_nucleotide(range,idx,value) { bit_array_assign (&(range)->ref, (idx) * 2,      acgt_encode[(uint8_t)value] & 1)       ;  \
                                              bit_array_assign (&(range)->ref, (idx) * 2 + 1, (acgt_encode[(uint8_t)value] & 2) >> 1) ; }

#define ref_is_nucleotide_set(range,idx) ((bool)bit_array_get (&(range)->is_set, (idx)))

#define ref_is_idx_in_range(range,idx) ((idx) < (range)->ref.nbits / 2)

#define ref_get_nucleotide(range,idx)  acgt_decode[(bit_array_get (&(range)->ref, (idx) * 2 + 1) << 1) | \
                                                    bit_array_get (&(range)->ref, (idx) * 2)]

static inline void ref_assert_nucleotide_available (const Range *range, PosType pos) {
    bool available;
    switch (flag.reference) {
        case REF_STORED    : available = ref_is_nucleotide_set (range, pos2range_idx (pos)); break;
        case REF_INTERNAL  : available = true; break;
        default            : available = (pos >= range->first_pos && pos <= range->last_pos); break;
    }
    ASSERT (available, "reference is not set: chrom=%.*s pos=%"PRId64, (range)->chrom_name_len, (range)->chrom_name, (pos));
}

// display
typedef struct { char s[300]; } RangeStr;
extern RangeStr ref_display_range (const Range *r);
extern void ref_print_bases_region (FILE *file, ConstBitArrayP bitarr, ConstBitArrayP is_set, PosType first_pos, uint64_t start_base, uint64_t num_of_bases, bool is_forward);
extern void ref_print_subrange (const char *msg, const Range *r, PosType start_pos, PosType end_pos, FILE *file);
extern void ref_print_is_set (const Range *r, PosType around_pos, FILE *file);

// globals
extern const char REVCOMP[];
extern const char *ref_filename;
extern Digest ref_file_md5;
extern Buffer ref_stored_ra;
extern Buffer loaded_contigs;
extern BitArrayP genome, emoneg, genome_is_set;
extern PosType genome_nbases;

#endif
