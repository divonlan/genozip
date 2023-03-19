// ------------------------------------------------------------------
//   reference.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "buffer.h"
#include "digest.h"
#include "bits.h"
#include "flags.h"
#include "sections.h"
#include "mutex.h"

#pragma GENDICT_PREFIX REF
#pragma GENDICT REF_CONTIG=DTYPE_FIELD=CONTIG 

// reference sequences - 
// Thread safety: "ref" is set atomically as the last set in initialization. If its set, then it is correct and will never change. If it appears to be
// not set yet, it is necesarry to lock the mutex and test again, and initialize if still not set.

typedef struct Range {
    Bits ref;                    // actual reference data - 2-bit array
    Bits is_set;                 // a 1-bit array - SEG: a pos is set if seg set this reference PIZ: is set if SEC_REF_IS_SET said so
    int64_t num_set;             // used by ref_prepare_range_for_compress: number of set bits in in_set
    STR (chrom_name);
    WordIndex chrom;             // index to the contig of the in the CHROM of the file from which this reference was loaded.
    uint32_t range_id;           // index of range within Buffer ranges
    uint32_t range_i;            // range ordinal number within contig
    PosType64 first_pos, last_pos; // the range that includes all loci (note: in ZIP-INTERNAL it might include unset loci too)
    PosType64 gpos;                // position of this range in the "global position" 
} Range;

#define ref_size(r) ((r) ? ((r)->last_pos - (r)->first_pos + 1) : 0)

// locks
typedef struct { int32_t first_mutex, last_mutex, first_mutex_rr, last_mutex_rr; } RefLock;
#define REFLOCK_NONE ((RefLock){-1,-1,-1,-1})

extern RefLock ref_lock (Reference ref, PosType64 gpos_start, uint32_t seq_len);
extern RangeP ref_seg_get_range (VBlockP vb, Reference ref, WordIndex chrom, STRp(chrom_name), PosType64 pos, uint32_t ref_consumed, WordIndex ref_index, rom field /* used for ASSSEG */, RefLock *lock);
extern void ref_unlock (Reference ref, RefLock *lock);

// replace range if POS has moved to next rage
extern RangeP ref_seg_renew_locked_range_do (VBlockP vb, Reference ref, RangeP range, PosType64 pos, PosType64 seq_len, RefLock *lock);
static inline RangeP ref_seg_renew_locked_range (VBlockP vb, Reference ref, RangeP range, PosType64 pos, PosType64 seq_len, RefLock *lock)
{
    return (!range || range->last_pos < pos) ? ref_seg_renew_locked_range_do (vb, ref, range, pos, seq_len, lock) : range/*same range*/;
}

typedef enum { RT_NONE,     // value of ranges.param if ranges is unallocated
               RT_MAKE_REF, // used in --make-ref one range per vb of fasta reference file - ranges in order of the fasta file
               RT_DENOVO,   // used in SAM with REF_INTERNAL 
               RT_LOADED    // one Range per chrom (contig), overlayed on genome
             } RangesType;
              
extern void ref_initialize_ranges (Reference ref, RangesType type);
extern void ref_finalize (void);
extern void ref_load_external_reference (Reference ref, ContextP chrom_ctx);
extern void ref_load_stored_reference (Reference ref);
extern bool ref_is_loaded (const Reference ref);
extern bool ref_is_external_loaded (const Reference ref);
extern void ref_diff_ref (void);
extern void ref_set_reference (Reference ref, rom filename, ReferenceType ref_type, bool is_explicit);
extern void ref_set_ref_file_info (Reference ref, Digest genome_digest, bool is_adler, rom fasta_name, uint8_t genozip_version);
extern void ref_unload_reference (Reference ref);
extern void ref_destroy_reference (Reference ref);
extern ConstRangeP ref_piz_get_range (VBlockP vb, Reference ref, FailType soft_fail);
extern RangeP ref_get_range_by_ref_index (VBlockP vb, Reference ref, WordIndex ref_contig_index);
extern rom ref_get_cram_ref (const Reference ref);
extern void ref_generate_reverse_complement_genome (Reference ref);
extern rom ref_get_filename (const Reference ref);
extern uint8_t ref_get_genozip_version (const Reference ref);
extern BufferP ref_get_stored_ra (Reference ref);
extern Digest ref_get_genome_digest (const Reference ref);
extern rom ref_get_digest_name (const Reference ref);
extern void ref_get_genome (Reference ref, const Bits **genome, const Bits **emoneg, PosType64 *genome_nbases);
extern void ref_set_genome_is_used (Reference ref, PosType64 gpos, uint32_t len);
extern bool ref_is_digest_adler (const Reference ref);

// ZIPping a reference
extern void ref_compress_ref (void);

// make-reference stuff
extern bool is_ref (STRp(data), bool *need_more);
extern void ref_make_ref_init (void);
extern void ref_make_seg_initialize (VBlockP vb);
extern void ref_consume_ref_fasta_global_area (void);
extern void ref_make_create_range (VBlockP vb);
extern void ref_make_after_compute (VBlockP vb);
extern ConstBufferP ref_make_get_contig_metadata (void);
extern void ref_make_genozip_header (SectionHeaderGenozipHeaderP header);
extern void ref_make_finalize (bool unused);
extern void ref_fasta_to_ref (FileP file);

// contigs stuff
extern void ref_contigs_populate_aligned_chroms (void);
extern WordIndex ref_contigs_get_by_name (const Reference ref, STRp(chrom_name), bool alt_ok, FailType soft_fail);
extern WordIndex ref_contigs_get_matching (const Reference ref, PosType64 LN, STRp(txt_chrom), STRp(*ref_contig), bool strictly_alt, bool *is_alt, int32_t *chrom_name_growth);
extern rom ref_contigs_get_name (const Reference ref, WordIndex ref_index, unsigned *contig_name_len);
extern ContigPkgP ref_get_ctgs (Reference ref);
extern uint32_t ref_num_contigs (const Reference ref);
extern PosType64 ref_contigs_get_contig_length (const Reference ref, WordIndex ref_contig_index, STRp(chrom_name), bool enforce);
extern WordIndex ref_contigs_ref_chrom_from_header_chrom (const Reference ref, STRp(chrom_name), PosType64 *hdr_LN);
extern void ref_contigs_load_contigs (Reference ref);

extern uint32_t ref_contigs_get_num_contigs (Reference ref);
extern PosType64 ref_contigs_get_genome_nbases (Reference ref);

extern WordIndex ref_contig_get_by_gpos (const Reference ref, PosType64 gpos, int32_t seq_len, PosType32 *pos);

extern const uint8_t acgt_encode[256];
extern const uint8_t acgt_encode_comp[256];

// note that the following work on idx and not pos! (idx is the index within the range)
static inline void ref_set_nucleotide (RangeP range, uint32_t idx, uint8_t value) 
    { bits_assign2 (&range->ref, idx*2, acgt_encode[value]); }

static inline bool ref_is_nucleotide_set (ConstRangeP range, uint32_t idx) { return (bool)bits_get (&range->is_set, idx); }

static inline bool ref_is_idx_in_range (ConstRangeP range, uint32_t idx) { return idx < range->ref.nbits / 2; }

static inline char acgt_decode (uint8_t b2) { switch (b2) { case 0:return'A' ; case 1:return'C' ; case 2:return'G' ; default:return'T'; }}
static inline char base_by_idx (ConstBitsP bitarr, uint64_t idx) { return acgt_decode (bits_get2 (bitarr, idx * 2)); } 
static inline char ref_base_by_idx (ConstRangeP range, uint64_t idx) { return base_by_idx (&range->ref, idx); }
static inline char ref_base_by_pos (ConstRangeP range, PosType64 pos)  { return ref_base_by_idx (range, pos - range->first_pos); }

static inline void ref_assert_nucleotide_available (ConstRangeP range, PosType64 pos) {
    bool available;
    switch (flag.reference) {
        case REF_STORED   : available = ref_is_nucleotide_set (range, pos); break;
        default           : available = (pos >= range->first_pos && pos <= range->last_pos); break;
    }
    ASSERT (available, "reference is not set: chrom=%.*s pos=%"PRId64, (range)->chrom_name_len, (range)->chrom_name, (pos));
}

#define REF(idx) ref_base_by_idx (range, (idx))
#define REFp(pos) ref_base_by_idx (range, (pos)-range->first_pos)

// round-robin around range
#define RR_IDX(idx) ( ((idx) > range_len) ? ((idx) - range_len) \
                    : ((idx) < 0)         ? ((idx) + range_len) \
                    :                        (idx)) // faster than mod (hopefully)

// display
typedef struct { char s[300]; } RangeStr;
extern RangeStr ref_display_range (ConstRangeP r);
extern void ref_display_all_ranges (Reference ref);
extern void ref_print_bases_region (FILE *file, ConstBitsP bitarr, ConstBitsP is_set, PosType64 first_pos, uint64_t start_base, uint64_t num_of_bases, bool is_forward);
extern void ref_print_subrange (rom msg, ConstRangeP r, PosType64 start_pos, PosType64 end_pos, FILE *file);
extern void ref_print_is_set (ConstRangeP r, PosType64 around_pos, FILE *file);
extern char *ref_dis_subrange (Reference ref, ConstRangeP r, PosType64 start_pos, PosType64 len, char *seq, bool revcomp);
extern void ref_display_ref (const Reference ref);

// globals
extern Reference gref, prim_ref;
extern Serializer make_ref_merge_serializer;
 