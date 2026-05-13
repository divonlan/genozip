// ------------------------------------------------------------------
//   refhash_make.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "compressor.h" // codec.h included in compressor.h must be included before reference.h
#include "file.h"
#include "dispatcher.h"
#include "zfile.h"
#include "refhash_friend.h"
#include "threads.h"
#include "reference.h"
#include "fasta_friend.h"

// ref_hash logic:
// we use the 28 bits (14 nucleotides) following a "G" hook, as the hash value, to index into a hash table with multiple
// layers - each layer being of size of half of the size of the previous layer.
// the full hash value is used to index into the first hash table layer, the 27 LSb of the hash value into the 2nd layer etc.
//
// The size of the first hash layer is 256M entries (28 bits) * sizeof(uint32_t) = 1GB, and all layers < 2 GB. 
//
// With the simplified assumption that G appears 1/4 of all bases, then a human genome will have ~ 3B / 4 =~ 750M hash values 
// competing to enter these 256M (one layer) - 512M slots (28 layers).
//
// The expected number of times any particular 29 bit ("G" + 28) string will appear in a 3B genome is 3B / 750M = 4.
// On a typical 150 bases short read, we expected 150/4 =~ 37.5 'G' hooks - hopefully one of them will point to the correct
// reference gpos
//

#define BASES_PER_HASH 18 // Can be up to 30 (not 32 bc of NO_KMER). was 14 up to 15.0.80 (see internal-docs/internal-docs/2026-03-14-refhash-stats)

#define VBLOCK_MEMORY_REFHASH (10 MB) // amount of refhash data each each vb compresses.

// make-ref: array of ref_hash_len. During setup, used to count how many instances of each hash value appear in the reference. 
// Thereafter, used to select which of the instances gets to have their gpos occupy in the hash table
static Buffer hash_hits_by_entry = {}; 

static uint64_t hooks_in_genome     = 0; // number of hooks in genome (not counting hooks on the boundary between make-ref 1 Mbps sranges - 0.0017% of hooks)
static uint64_t hooks_in_hash_table = 0; // number of hooks whose GPOS is in the hash table (the rest were dropped due to hash collisions - more than one hook wishes to occupy the same hash entry)

static bool refhash_calc_bits_per_hash_out (void)
{
    #define LARGEST_BITS 32 // our largest refhash. in "default" mode: 2^32 x gpos_bytes=5 = 20 GB RAM)

    uint64_t genome_nbases = ROUNDUP64 (ref_contigs_get_genome_nbases()) + 64; // note: gref.genome_nbases is not set yet
    ASSERTNOTZERO (genome_nbases);

    uint64_t threshhold = 1 MB; // 17 bits. note: for human (< 4GB threshold) it will be 29 bits.
    for (bits_per_hash_out=17; bits_per_hash_out < LARGEST_BITS; bits_per_hash_out++, threshhold *= 2)
        if (genome_nbases <= threshhold) 
            break; // note: set to LARGEST_BITS if loop completes with no threshold condition met

    // now modify based on command line option
    switch (flag.make_reference) {
        case MAKE_REF_TINY   : bits_per_hash_out -= 2; break;
        case MAKE_REF_MEDIUM :                         break;
        case MAKE_REF_SMALL  : bits_per_hash_out -= 1; break;
        case MAKE_REF_LARGE  : bits_per_hash_out += 1; break;
        
        default: ABORT ("Invalid make_reference=%u", flag.make_reference);
    }

    hash_shift = (64 - bits_per_hash_out);

    return true;
}

void refhash_make_initialize (void)
{
    uint64_t max_gpos = ref_get_max_gpos();

    // we cannot use the reference data beyond the first MAX_ALIGNER_GPOS for creating the refhash
    if (max_gpos > MAX_ALIGNER_GPOS)
        WARN_ONCE ("%s contains %s bases. When compressing a FASTQ or unaligned (i.e. missing RNAME, POS) "
                   "SAM/BAM file using the reference being generated, only the first %s bases of the reference will be used, "
                   "possibly affecting the compression ratio. This limitation doesn't apply to aligned SAM/BAM files and VCF files. "
                   "If you need to use reference files larger than %s, let us know at " EMAIL_SUPPORT ".", 
                   txt_name, str_bases (max_gpos).s, str_bases (MAX_ALIGNER_GPOS+1).s, str_bases (MAX_ALIGNER_GPOS+1).s);

    gpos_bytes = max_gpos < NO_HASH_ENT32 ? 4 : 5; // either up to 4 Gps or up to 1 Tbp

    refhash_calc_bits_per_hash_out();

    buf_alloc_exact_zero (evb, hash_hits_by_entry, ref_hash_len, uint8_t, "hash_hits_by_entry");
}

// compute thread: pass two

// main thread: the "read" part of pass two
static void refhash_prepare_for_decide_occupier (VBlockP vb_)
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;

    // since there is no "read" time in this dispatcher, each core (roughly) get to do one iteration
    uint32_t ents_per_vb = (ref_hash_len / global_max_threads) + 1;

    if (vb->vblock_i > global_max_threads) return; // we're done

    vb->first_ent = (vb->vblock_i - 1) * ents_per_vb;
    vb->num_ents  = MIN_(ents_per_vb, ref_hash_len - vb->first_ent);
    vb->dispatch  = READY_TO_COMPUTE;
}

// compute thread of pass two
static void refhash_decide_occupier_one_vb (VBlockP vb_)
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;
    START_TIMER;

    ARRAY (uint8_t, counts, hash_hits_by_entry);

    // decide which of the counts[i] k-mers (1 to counts[i]) of this hash in the reference data gets 
    // to have its gpos occupy the hash table, and reduce counts[i] to that number.
    // note: the number is not sequential in the reference, but by in some random order in which compute threads 
    //       consider hash[i] refhash_calc_one_range.
    // note: for extremely common k-mers - we ignore those beyond the first 255. while it might seem as a bias
    // towards the reference beginning, in fact this is a useless entry anyway, as the chances of the aligner hitting the right one are low 
    if (!flag.show_ref_hash) { // fast loop if no stats collection
        for (uint32_t i=vb->first_ent; i < vb->first_ent + vb->num_ents; i++)
            if (counts[i] > 1) 
                counts[i] = (i % counts[i]) + 1; // ∈ [1,counts[i]] 
    }

    else
        for (uint32_t i=vb->first_ent; i < vb->first_ent + vb->num_ents; i++)
            if (counts[i] > 0) {
                vb->hooks_in_hash_table++;
                vb->hooks_in_genome += counts[i];

                if (counts[i] > 1) 
                    counts[i] = (i % counts[i]) + 1; // ∈ [1,counts[i]] 
            }

    vb_set_is_processed (VB); // tell dispatcher this thread is done and can be joined.

    COPY_TIMER (refhash_p2_decide_occupier);
}

static void refhash_aggregate_stats (VBlockP vb_)
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;

    hooks_in_genome     += vb->hooks_in_genome;
    hooks_in_hash_table += vb->hooks_in_hash_table;
}

// increment a value atomically, but never past 255
static inline void increment_relaxed_ceiling_255 (uint8_t *addr)
{
    uint8_t value = load_relaxed (*addr);
    while (value < 255) // don't try to increment if already 255
        // increment, but if another thread beat me to it, loop back and try again (value is updated if failed)
        if (__atomic_compare_exchange_n (addr, &value, value+1, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) 
            break;
}

// decrement a value atomically, but never beneath 0. return true if we decremented to 0
static inline bool decrement_relaxed_floor_0 (uint8_t *addr)
{
    uint8_t value = load_relaxed (*addr);
    while (value > 0) // loop until we successfully decrement value, or another thread decrements it down to zero
        // decrement, but if another thread beat me to it, loop back and try again (value is updated, but only if failed)
        if (__atomic_compare_exchange_n (addr, &value, value-1, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) 
            return (value == 1); // true if it is we that decremented the value to 0 (i.e. the expected value was 1), false if we decremented but not to zero

    return false; // another thread decremented the value to zero, not us
}

// make-reference compute thread: called by two tasks: pass one and pass three of refhash generation
void refhash_calc_one_range (VBlockP vb_, RefhashCalcType calc_type) // VB of reference or refhash dispatcher
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;

    START_TIMER;

    ConstRangeP r = vb->range;
    PosType64 num_bases = ref_size (r);
    
    ASSERT (num_bases * 2 == r->ref.nbits, "mismatch between num_bases=%"PRId64" (x2 = %"PRId64") and r->ref.nbits=%"PRIu64". Expecting the latter to be exactly double the former. chrom=%s r->first_pos=%"PRId64" r->last_pos=%"PRId64" r->range_id=%u", 
            num_bases, num_bases*2, r->ref.nbits, Bc (ZCTX(0)->dict, B(CtxNode, ZCTX(0)->nodes, r->chrom)->char_index), r->first_pos, r->last_pos, r->range_id);

    if (r->gpos + (num_bases - BASES_PER_HASH) > MAX_ALIGNER_GPOS)
        return; // ignore the tail end of the reference that exceeds the maximum GPOS

    ARRAY (uint8_t, counts, hash_hits_by_entry);
    ARRAY (uint32_t, refhash32, refhash_buf); // add restrict?
    ARRAY (uint40_t, refhash40, refhash_buf); // add restrict?

    #define encoded_ref_base(idx) bits_get2 (&r->ref, (idx) * 2)

    // note: we lose a handful of hooks which might be in the final BASES_PER_HASH of each range - negligible vs the number of hooks we lose due to hash collisions
    for (PosType64 base_i=0; base_i < num_bases - BASES_PER_HASH; base_i++)

        // take only the final hook in a homopolymer of hooks (i.e. the last G in a e.g. GGGGG)
        if (encoded_ref_base (base_i) == encoded_HOOK && encoded_ref_base (base_i+1) != encoded_HOOK) {
            uint64_t kmer = bits_get_wordn (&r->ref, (base_i+1) * 2, BASES_PER_HASH * 2); // starting from the base after the hook
            uint64_t hash = fibonacci_hash (kmer);

            if (calc_type == REFHASH_COUNT_INSTANCES) // just count (pass 1)
                increment_relaxed_ceiling_255 (&counts[hash]); // increment count, but not past 255

            else { // (calc_type == REFHASH_OCCUPY) - place gpos values in the hash table (pass 3)
                // logic: when a kmer in reference that maps to this hash entry arrives here, we decrement count, 
                // and it is the kmer the is lucky one to reduce count to 0, that gets to occupy hash with its GPOS.
                if (decrement_relaxed_floor_0 (&counts[hash])) {
                    if (gpos_bytes == 4) refhash32[hash] = LTEN32 (r->gpos + base_i);
                    else /* 5 */         refhash40[hash] = LTEN40 (r->gpos + base_i);
                }
            }
        }

    if (calc_type == REFHASH_OCCUPY) COPY_TIMER (refhash_p3_occupy); else COPY_TIMER (refhash_p1_count);
}

static void refhash_occupy_one_vb (VBlockP vb)
{
    refhash_calc_one_range (vb, REFHASH_OCCUPY);

    dispatcher_increment_progress ("refhash_occupy_one_vb", vb->range->ref.nbits / 2); 

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// compress the reference - one section at the time, using Dispatcher to do them in parallel (make_ref version)
static void refhash_prepare_for_compress (VBlockP vb_)
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;

    static uint64_t next_first_ent = 0;
    if (vb->vblock_i == 1) next_first_ent = 0; // initialize

    if (next_first_ent == ref_hash_len) return; // we're done
        
    // tell this vb what to do
    vb->first_ent = next_first_ent;
    vb->num_ents  = MIN_(VBLOCK_MEMORY_REFHASH / gpos_bytes, ref_hash_len - next_first_ent); 
    next_first_ent += vb->num_ents;

    dispatcher_increment_progress ("refhash_prepare_for_compress", vb->num_ents); 

    vb->dispatch = READY_TO_COMPUTE;
}

// part of --make-reference - compute thread for compressing part of the hash
static void refhash_compress_one_vb (VBlockP vb_)
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;

    START_TIMER;

    SectionHeaderRefHash header = { .magic                 = BGEN32 (GENOZIP_MAGIC),
                                    .section_type          = SEC_REF_HASH, 
                                    .codec                 = CODEC_RANW, // Much!! faster than LZMA (compress and uncompress), 8% worse compression of human refs, MUCH better on small refs (sparse hash)
                                    .data_uncompressed_len = BGEN32 (vb->num_ents * gpos_bytes),
                                    .vblock_i              = BGEN32 (vb->vblock_i),
                                    .first_ent             = vb->first_ent }; // no BGEN - bit field

    comp_compress (VB, NULL, &vb->z_data, &header, Bc(refhash_buf, gpos_bytes * vb->first_ent), NO_CALLBACK, "SEC_REF_HASH");

    dispatcher_increment_progress ("refhash_compress_one_vb", vb->num_ents); 

    vb_set_is_processed (VB); // tell dispatcher this thread is done and can be joined.

    COPY_TIMER(refhash_p4_compress);
}

static void refhash_compress_write (VBlockP vb)
{
    zfile_output_processed_vb (vb);

    dispatcher_increment_progress ("refhash_compress_write", VB_FASTA->num_ents); 
}

// make-reference: called by main thread in zip_write_global_area
void refhash_make_refhash (void)
{
    START_TIMER;

    // note: time estimate breakdown: 50% zip, 17% write reference, 16% occupy, 17% write refhash

    // pass one: (dispatch ranges) calls refhash_calc_one_range() with REFHASH_COUNT_INSTANCES to count 
    // how many instances each hash appears in the reference data (called from ref_compress_one_range)

    // pass two: (dispatch refhash blocks) decide which of the instances is going to occupy each hash table entry
    dispatcher_fan_out_task (TASK_MRH_DECIDE_OCCUPIER, NULL, 0, JOIN_OUT_OF_ORDER, false, 0, 100,
                             refhash_prepare_for_decide_occupier, 
                             refhash_decide_occupier_one_vb, 
                             refhash_aggregate_stats);

    // now we will need refhash_buf, so its time to allocate it
    refhash_buf.can_be_big = hash_hits_by_entry.can_be_big = true; // suppress warning
    if (gpos_bytes == 4) buf_alloc_exact_255 (evb, refhash_buf, ref_hash_len, uint32_t, "refhash_buf"); 
    else                 buf_alloc_exact_255 (evb, refhash_buf, ref_hash_len, uint40_t, "refhash_buf"); 

    // pass three: (dispatch ranges) occupy the refhash table
    dispatcher_fan_out_task (TASK_MRH_OCCUPY, NULL,
                             /*target_progress=*/ref_contigs_get_genome_nbases(),
                             JOIN_OUT_OF_ORDER, false, 0, 100,
                             ref_make_prepare_one_range_for_dispatch, 
                             refhash_occupy_one_vb,
                             NO_CALLBACK);

    // we can free a bunch memory now
    buf_destroy (hash_hits_by_entry);
    ref_destroy_genome();
    return_freed_memory_to_kernel();

    // pass four: (dispatch refhash blocks) compress and write the hash table to the reference file being created
    dispatcher_fan_out_task (TASK_MRH_COMPRESS, NULL, 
                             /*target_progress=*/ref_hash_len * 3, 
                             JOIN_OUT_OF_ORDER, // SEC_REF_HASH sections may be written out-of-order
                             false, 0, 5000,
                             refhash_prepare_for_compress, 
                             refhash_compress_one_vb, 
                             refhash_compress_write);

    // pass five: (main thread) calculate a digest of the ref table
    { START_TIMER;
    // note: by design, the refhash table and its digest are expected to be different in each run, because the gpos that wins contention depends on the random order of threads
    refhash_digest = digest_do (B1STc(refhash_buf), refhash_buf.len * gpos_bytes, flag.md5 ? DIGEST_MD5 : DIGEST_XXH3, "refhash"); // gets written to genozip_headers
    COPY_TIMER_EVB (refhash_p5_digest); }

    if (flag.show_ref_hash) {
        iprintf ("\n\nReference hash summary for %s:\n", z_name);
        iprintf ("Algorithm parameters:\n\trefhash_size = %s\n\tbases_per_hash_in = %u\n\tbits_per_hash_out = %u\n\tn_hash_ents = %s\n\n", 
                 make_ref_size_name (flag.make_reference), BASES_PER_HASH, bits_per_hash_out, str_int_commas (ref_hash_len).s);
        iprintf ("FASTA:\n\tgenome_nbases = %s\n\tgpos_bytes = %u\n\thooks_in_genome = %s (%1.1f%% of bases)\n\n", // note: genome_nbases include padding on every range, so a bit more than the number of biological bases
                 str_int_commas (*p_genome_nbases).s, (int)gpos_bytes, str_int_commas (hooks_in_genome).s, percent (hooks_in_genome, *p_genome_nbases)); 
        iprintf ("Performance:\n\thash_ents_occupied = %1.1f%%\n\thooks_in_hash_table = %s (%1.1f%% of hooks in genome)\n",
                 percent (hooks_in_hash_table, ref_hash_len), str_int_commas (hooks_in_hash_table).s, percent (hooks_in_hash_table, hooks_in_genome));
    }

    COPY_TIMER_EVB (refhash_make_refhash);
}

// callback from zfile_compress_genozip_header (for reference file GENOZIP_HEADER)
void ref_make_genozip_header (SectionHeaderGenozipHeaderP header)
{
    header->genome_digest         = z_file->digest;
    header->ref.bits_per_hash_out = bits_per_hash_out;
    header->ref.gpos_bytes        = gpos_bytes;
    header->ref.bases_per_hash    = BASES_PER_HASH;
    header->ref.make_ref_size     = flag.make_reference;
}
