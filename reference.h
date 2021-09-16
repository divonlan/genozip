// ------------------------------------------------------------------
//   reference.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef REFERENCE_INCLUDED
#define REFERENCE_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "digest.h"
#include "bit_array.h"
#include "flags.h"
#include "sections.h"

// reference sequences - one per range of 1MB. ranges (chrom, pos) are mapped here with a hash function. In the rare case two unrelated ranges
// are mapped to the same entry - only the first range will have a reference, and the other ones will not. this will hurt the compression ratio,
// but not the correctness.
// Thread safety: "ref" is set atomically as the last set in initialization. If its set, then it is correct and will never change. If it appears to be
// not set yet, it is necesarry to lock the mutex and test again, and initialize if still not set.

typedef struct Range {
    uint32_t range_id;           // index of range within Buffer ranges
    BitArray ref;                // actual reference data - 2-bit array
    BitArray is_set;             // a 1-bit array - SEG: a pos is set if seg set this reference PIZ: is set if SEC_REF_IS_SET said so
    int64_t num_set;             // used by ref_prepare_range_for_compress: number of set bits in in_set
    const char *chrom_name;
    unsigned chrom_name_len;
    WordIndex chrom;             // index to the chrom of the external reference
    uint32_t range_i;            // range ordinal number within contig
    PosType first_pos, last_pos; // the range that includes all locii (note: in ZIP-INTERNAL it might include unset locii too)
    PosType gpos;                // position of this range in the "global position" 
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

extern RefLock ref_lock (Reference ref, PosType gpos_start, uint32_t seq_len);
extern RefLock ref_unlock (Reference ref, RefLock lock);
extern RefLock ref_lock_range (Reference ref, int32_t range_id);

typedef enum { RT_NONE,      // value of ranges.param if ranges is unallocated
               RT_MAKE_REF,  // used in --make-ref one range per vb of fasta reference file - ranges in order of the fasta file
               RT_DENOVO,    // used in SAM with REF_INTERNAL - an large array of Range's, hashed by the chrom and pos
               RT_LOADED,    // one Range per chrom (contig), overlayed on genome
               RT_CACHED     // same as RT_LOADED, but data is mmap'ed to a cache file
             } RangesType;
              
extern void ref_initialize_ranges (Reference ref, RangesType type);
extern void ref_load_external_reference (Reference ref, ContextP chrom_ctx);
extern void ref_load_stored_reference (Reference ref);
extern bool ref_is_loaded (const Reference ref);
extern bool ref_is_external_loaded (const Reference ref);
extern void ref_display_ref (const Reference ref);
extern void ref_set_reference (Reference ref, const char *filename, ReferenceType ref_type, bool is_explicit);
extern void ref_set_ref_file_info (Reference ref, Digest md5, const char *fasta_name, uint8_t genozip_version);
extern void ref_unload_reference (Reference ref);
extern void ref_destroy_reference (Reference ref, bool destroy_only_if_not_mmap);
extern MemStats ref_memory_consumption (const Reference ref);
extern const Range *ref_piz_get_range (VBlockP vb, Reference ref, PosType first_pos_needed, uint32_t num_nucleotides_needed);
extern Range *ref_get_range_by_ref_index (VBlockP vb, Reference ref, WordIndex ref_chrom_index);
extern Range *ref_seg_get_locked_range (VBlockP vb, Reference ref, WordIndex chrom, const char *chrom_name, unsigned chrom_name_len, PosType pos, uint32_t seq_len, const char *field /* used for ASSSEG */, RefLock *lock);
extern const char *ref_get_cram_ref (const Reference ref);
extern void ref_generate_reverse_complement_genome (Reference ref);

extern const char *ref_get_filename (const Reference ref);
extern uint8_t ref_get_genozip_version (const Reference ref);
extern BufferP ref_get_stored_ra (Reference ref);
extern Digest ref_get_file_md5 (const Reference ref);
extern void ref_get_genome (Reference ref, const BitArray **genome, const BitArray **emoneg, PosType *genome_nbases);
extern void ref_set_genome_is_used (Reference ref, PosType gpos, uint32_t len);


// ZIPping a reference
extern void ref_compress_ref (void);

// make-reference stuff
extern void ref_make_ref_init (void);
extern void ref_consume_ref_fasta_global_area (void);
extern void ref_make_create_range (VBlockP vb);
extern void ref_make_after_compute (VBlockP vb);
extern ConstBufferP ref_make_get_contig_metadata (void);
extern void ref_make_finalize (void);

// cache stuff
extern bool ref_mmap_cached_reference (Reference ref);
extern void ref_create_cache_in_background (Reference ref);
extern void ref_create_cache_join (Reference ref);
extern void ref_remove_cache (Reference ref);

// contigs stuff
typedef struct { char s[REFCONTIG_MD_LEN]; } ContigMetadata;

extern WordIndex ref_contigs_get_by_name (Reference ref, const char *chrom_name, unsigned chrom_name_len, bool alt_ok, bool soft_fail);
extern WordIndex ref_contigs_get_by_name_or_alt (Reference ref, const char *chrom_name, unsigned chrom_name_len, bool soft_fail);

extern ConstBufferP ref_get_contigs (const Reference ref);
extern void ref_contigs_get (const Reference ref, ConstBufferP *out_contig_dict, ConstBufferP *out_contigs);
extern uint32_t ref_num_loaded_contigs (const Reference ref);
extern PosType ref_contigs_get_contig_length (const Reference ref, WordIndex chrom_index, const char *chrom_name, unsigned chrom_name_len, bool enforce);
extern WordIndex ref_contigs_ref_chrom_from_header_chrom (const Reference ref, const char *chrom_name, unsigned chrom_name_len, PosType *last_pos, WordIndex header_chrom);
extern void ref_contigs_sort_chroms (void);
extern void ref_contigs_load_contigs (Reference ref);

typedef void (*RefContigsIteratorCallback)(const char *chrom_name, unsigned chrom_name_len, PosType last_pos, void *callback_param);
extern void ref_contigs_iterate (const Reference ref, RefContigsIteratorCallback callback, void *callback_param);
extern WordIndex ref_contig_get_by_gpos (const Reference ref, PosType gpos, PosType *pos);
extern WordIndex ref_contig_get_by_chrom (ConstVBlockP vb, const Reference ref, WordIndex txt_chrom_index, const char *chrom_name, unsigned chrom_name_len, PosType *max_pos);

// alt chroms stuff
extern void ref_alt_chroms_load (Reference ref);
extern void ref_alt_chroms_compress (Reference ref);
extern WordIndex ref_alt_chroms_get_alt_index (Reference ref, const char *chrom, unsigned chrom_len, PosType chrom_LN, WordIndex fallback_index);

// gets the index of the matching chrom in the reference - either its the chrom itself, or one with an alternative name
// eg 'chr22' instead of '22'
#define ref_alt_get_final_index(chrom_index) \
    (buf_is_alloc (&z_file->alt_chrom_map) ? *ENT (WordIndex, z_file->alt_chrom_map, (chrom_index)) : (chrom_index))

// note that the following work on idx and not pos! (idx is the index within the range)
#define ref_set_nucleotide(range,idx,value) { bit_array_assign (&(range)->ref, (idx) * 2,      acgt_encode[(uint8_t)value] & 1)       ;  \
                                              bit_array_assign (&(range)->ref, (idx) * 2 + 1, (acgt_encode[(uint8_t)value] & 2) >> 1) ; }

#define ref_is_nucleotide_set(range,idx) ((bool)bit_array_get (&(range)->is_set, (idx)))

#define ref_is_idx_in_range(range,idx) ((idx) < (range)->ref.nbits / 2)

#define ref_base_by_idx(range,idx)  acgt_decode[(bit_array_get (&(range)->ref, (idx) * 2 + 1) << 1) | \
                                                 bit_array_get (&(range)->ref, (idx) * 2)]
#define ref_base_by_pos(range,pos) ref_base_by_idx((range), (pos)-(range)->first_pos)

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
extern char *ref_dis_subrange (Reference ref, const Range *r, PosType start_pos, PosType len, char *seq, bool revcomp);

// globals
extern Reference gref, prim_ref;

#endif
