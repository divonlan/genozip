// ------------------------------------------------------------------
//   reference.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REFERENCE_INCLUDED
#define REFERENCE_INCLUDED

#include <pthread.h>
#include "genozip.h"
#include "buffer.h"
#include "md5.h"
#include "bit_array.h"
#include "compressor.h"

// reference sequences - one per range of 1MB. ranges (chrom, pos) are mapped here with a hash function. In the rare case two unrelated ranges
// are mapped to the same entry - only the first range will have a reference, and the other ones will not. this will hurt the compression ratio,
// but not the correctness.
// Thread safety: "ref" is set atomically as the last set in initialization. If its set, then it is correct and will never change. If it appears to be
// not set yet, it is necesarry to lock the mutex and test again, and initialize if still not set.

#define REF_NUM_RANGES (1 << 20)
#define REF_NUM_SITES_PER_RANGE (1 << 20) // 1 Mbp

typedef struct Range {
    BitArray ref;                // actual reference data - 2-bit array
    BitArray is_set;             // a 1-bit array - SEG: a pos is set if seg set this reference PIZ: is set if SEC_REF_IS_SET said so
    const char *chrom_name;
    unsigned chrom_name_len;
    int32_t chrom;               // index to the chrom of the external reference
    uint32_t range_i;
    int64_t first_pos, last_pos; // the range that includes all locii (note: in ZIP-INTERNAL it might include unset locii too)
    int64_t gpos;                // position of this range in the "global position" 
    uint32_t copied_first_index, copied_len; // ZIP with REF_EXT_STORE: the subset of this range that was copied directly from the fasta file and doesn't need to be compressed
    pthread_mutex_t mutex;
} Range;

#define ref_size(r) ((r) ? ((r)->last_pos - (r)->first_pos + 1) : 0)

typedef enum { RT_SMALL_RANGES, RT_RANGE_PER_CONTIG, RT_WHOLE_GENOME } RangesType;
extern void ref_initialize_ranges (RangesType ranges_type);
extern void ref_compress_ref (void);
extern void ref_load_stored_reference (void);
extern void ref_set_reference (const char *filename);
extern void ref_set_md5 (Md5Hash md5);
extern void ref_load_external_reference (bool display);
extern void ref_unload_reference (bool force_clean_all);
extern MemStats ref_memory_consumption (void);
extern const Range *ref_piz_get_range (VBlockP vb, int64_t first_pos_needed, uint32_t num_nucleotides_needed);
extern void ref_consume_ref_fasta_global_area (void);
extern void ref_contig_word_index_from_name (const char *chrom_name, unsigned chrom_name_len, const char **snip, int32_t *word_index);
extern void ref_get_contigs (ConstBufferP *out_contig_dict, ConstBufferP *out_contig_words);
extern int64_t ref_min_max_of_chrom (int32_t chrom, bool get_max);

extern Range *ref_zip_get_locked_range (VBlockP vb, int64_t pos);

extern void ref_print_subrange (const char *msg, const Range *r, int64_t start_pos, int64_t end_pos);
extern void ref_print_is_set (const Range *r);

// Make-reference stuff (ZIP of FASTA with --make-reference)
extern void ref_make_ref_init (void);
extern Range *ref_make_ref_get_range (uint32_t vblock_i);
extern void ref_output_vb (VBlockP vb);

// ZIP ONLY: access range_i and index within range, for ranges configured for ZIP
#define ridx2pos(range_i,idx) (((int64_t)(range_i) * REF_NUM_SITES_PER_RANGE) | (idx))
#define range_i2pos(range_i) ((int64_t)(range_i) * REF_NUM_SITES_PER_RANGE)
#define pos2range_i(pos)   ((uint32_t)((pos) / REF_NUM_SITES_PER_RANGE))
#define pos2range_idx(pos) ((uint32_t)((pos) % REF_NUM_SITES_PER_RANGE))
#define pos2startrange(pos) ((pos) & ~(uint64_t)(REF_NUM_SITES_PER_RANGE-1)) // zeros the LSbs

#define ref_assert_nucleotide_available(range,pos) \
    ASSERT (/* piz w stored ref */ (flag_reference == REF_STORED && ref_is_nucleotide_set ((range), pos2range_idx(pos))) ||  \
            /* zip w ext ref    */ ((flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE) && ((pos) >= (range)->first_pos && (pos) <= (range)->last_pos)) || \
            /* zip internal ref */ flag_reference == REF_INTERNAL, \
        "Error in %s:%u: reference is not set: chrom=%.*s pos=%"PRId64, __FUNCTION__, __LINE__, (range)->chrom_name_len, (range)->chrom_name, (pos))

// note that the following work on idx and not pos! (idx is the index within the range)
#define ref_set_nucleotide(range,idx,value) { bit_array_assign (&(range)->ref, (idx) * 2,      acgt_encode[(uint8_t)value] & 1)       ;  \
                                              bit_array_assign (&(range)->ref, (idx) * 2 + 1, (acgt_encode[(uint8_t)value] & 2) >> 1) ; }

#define ref_is_nucleotide_set(range,idx) ((bool)bit_array_get (&(range)->is_set, (idx)))

#define ref_get_nucleotide(range,idx)   acgt_decode[(bit_array_get (&(range)->ref, (idx) * 2 + 1) << 1) | \
                                                     bit_array_get (&(range)->ref, (idx) * 2)]

// globals
extern const char *ref_filename;
extern Md5Hash ref_md5;
extern Buffer ref_stored_ra;

extern Range *genome, *genome_rev;
extern int64_t genome_size;

#define GENOME_BASES_PER_MUTEX (1 << 16) // 2^16 = 64K
#define GPOS2MUTEX(gpos) ((gpos) / GENOME_BASES_PER_MUTEX)
extern pthread_mutex_t *genome_muteces; // one spinlock per 16K bases - protects genome->is_set
extern uint32_t genome_num_muteces;

#endif
