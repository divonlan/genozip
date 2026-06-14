// ------------------------------------------------------------------
//   aligner.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "seg.h"
#include "refhash_friend.h"
#include "file.h"
#include "codec.h"
#include "reconstruct.h"
#include "piz.h"
#include "aligner.h"

// splicing
#define MIN_MATCH_PC       47  // fail alignment if less than this number of bases match
#define SPLICE_MIN_PC      32  // min and max match percent in which we search for splicing
#define SPLICE_MAX_PC      90  
#define SPLICE_MIN_SEQ_LEN 32  // no splicing for sequences shorter than this - not worth it
#define MAX_SPLICE_GAP     (32 KB - 1)  // max allowed gap beween gpos and gpos2
#define MIN_SPLICE_MATCH_CONTRIBTION 10 // for splicing to compress better, the spliced alignment needs to have this many more matching bases 

// pairing
#define PAIR_MAX_DISTANCE  500 // if pair2's gpos' gap from pair1's is larger than this, it will get a penalty
#define NON_PAIR_PENALTY   13  // penalty reduction in matching bases, for a pair2 gpos that is too far from pair1 gpos

// good enough result - stop searching
#define NEAR_PERFECT_NORM  3   // if we find a genome sequence with this many mismatches - we stop searching and declare "near perfect"
#define NEAR_PERFECT_FAST  20
#define NEAR_PERFECT_BEST  0

// low level
#define CPU_CACHE_LINE_BYTES 64
#define CPU_CACHE_LINE(addr) (~bitmask64_(6) & (uintptr_t)(addr)) // integer value of pointer to first byte in cache line

typedef enum { NOT_FOUND=-1, REVERSE=0, FORWARD=1 } Direction;

typedef struct {
    PosType64 gpos;
    int32_t n_matching; // matching bits per hamming
    int32_t score;        // matching_bits after applying penalty if needed (might be negative)
    bool is_forward; 
    bool is_perfect;
} BestAlignment;

#include "aligner_layered.c" // backward compatabilities with reference files up to 15.0.80

static inline uint32_t nuke_encode_dir (char c, bool is_forward) 
{ 
    return is_forward ? nuke_encode(c) : nuke_encode_comp(c); 
}

static inline int percent_match (uint32_t n_matching, uint32_t seq_len)
{
    return n_matching * 100 / seq_len; 
}

// Foward example: If seq is: G-AGGGCT  (G is the hook)  -- matches reference AGGGCT       - function returns 110110101000 (A=00 is the LSb)
// Reverse       : If seq is: CGCCCT-C  (C is the hook)  -- also matches reference AGGGCT  - function returns 110110101000 - the same
// calculates a refhash word from 14 nucleotides following a 'G' (only last G in a sequenece of GGGG...)
static inline uint64_t aligner_get_kmer (Bits𐤐 seq, uint64_t bit_i, bool is_forward)
{
    // note on 'N' in sequence: aligner_seq_to_bitmap() encoded non-ACGT as 00 ('A').
    if (!is_forward)
        bit_i = bit_i + 2 - bits_per_hash_in;

    uint64_t word_index = bitset64_wrd(bit_i);
    int word_offset     = bitset64_idx(bit_i);

    uint64_t kmer = seq->words[word_index] >> word_offset; // all the bits of the word, starting from bit_i (so the most significant bits)

    int bits_taken = WORD_SIZE - word_offset;
    int bits_still_needed = bits_per_hash_in - bits_taken; 
    
    // case: too many bits 
    if (bits_still_needed < 0) 
        kmer &= bitmask64_(bits_per_hash_in); // note: even though bits_per_hash_in is immuteable, calculating the bitmask by this macro is faster than looking up a global variable

    // case: not enough bits: get more from next word (caller must ensure sufficient data is available)    
    else if (bits_still_needed > 0) 
        kmer |= (seq->words[word_index+1] & bitmask64_(bits_still_needed)) << (WORD_SIZE - word_offset);

    if (!is_forward) 
        kmer = bits_revcomp_word (kmer) >> (64 - bits_per_hash_in);

    return kmer;
}

static __attribute__((always_inline)) inline 
uint64_t aligner_prefetch_gpos (VBlock𐤐 vb, Bits𐤐 seq, uint64_t base_i, bool is_forward)                                              
{   
    uint64_t kmer = aligner_get_kmer (seq, base_i*2, is_forward); 
    uint64_t hash = fibonacci_hash (kmer);

    // Locality=1 (only L3 caching) is best (tested). Reason: after reading the GPOS
    // from refhash, we probably will never access that location again, so caching just causes cache pollution.
    // However, we cache just in L3, in case the Line Fill Buffer (LFB) was already evicted by the time the code was ready to pick it up. 
    if (gpos_bytes == 4)
        __builtin_prefetch (B32(refhash_buf, hash), 0/*read-only*/, 1/*locality: cache just in L3*/);

    if (gpos_bytes == 5) {
        char *addr = (char *)B40(refhash_buf, hash);
        __builtin_prefetch (addr, 0, 1);
 
        // case: refhash entry spans two cache lines - prefetch the second one too (applicable to 1/16 of 5-byte words)
        if (CPU_CACHE_LINE (addr) != CPU_CACHE_LINE (addr + 4)) 
            __builtin_prefetch (addr + 4, 0, 1);
    }

    return hash;
}

// converts a string sequence to a 2-bit bitmap
static Bits aligner_seq_to_bitmap (VBlock𐤐 vb, rom𐤐 seq, uint64_t seq_len, 
                                   uint64_t *restrict bitmap_words, // allocated by caller
                                   bool𐤐 seq_is_all_acgt)  
{
    START_TIMER;

    // convert seq to 2-bit array
    Bits seq_bits = { .nbits  = seq_len * 2, 
                      .nwords = roundup_bits2words64(seq_len * 2), 
                      .words  = bitmap_words,
                      .type   = BUF_REGULAR };

    *seq_is_all_acgt = true; // starting optimistically

    for (uint64_t base_i=0; base_i < seq_len; base_i++) {
        uint8_t encoding = nuke_encode (seq[base_i]);
    
        if (encoding == 4) { // not A, C, G or T - usually N
            *seq_is_all_acgt = false;
            encoding = 0;    // arbitrarily convert 4 (any non-ACGT is 4) to 0 ('A')
        }
    
        bits_assign2 (&seq_bits, (base_i << 1), encoding);
    }

    bits_clear_excess_bits_in_top_word (&seq_bits, false); // bc bitmap_words is uninitialized

    COPY_TIMER (aligner_seq_to_bitmap);

    return seq_bits;
}

// return number of mismatching bases between seq_bis and the candidate genomic region
// note: we already verified (in aligner_collect_gpos) that seq_bits doesn't go before or after the edge of the genome
static inline uint32_t aligner_count_mismatches (ConstBits𐤐 seq_bits, // the entire bit array (restrict-ed!)
                                                 ConstBits𐤐 genome, 
                                                 PosType64 gpos,
                                                 bool is_forward)
{
    uint64_t genome_bit_i = (gpos * 2) + (is_forward ? 0 : seq_bits->nbits);  // index of bit in genome that matches first bit in seq_bits 

    const uint64_t *restrict words_1 = seq_bits->words;
    const uint64_t *restrict words_2 = &genome->words[genome_bit_i >> 6];
    const uint64_t *after_1 = words_1 + seq_bits->nwords;
    uint64_t diff=0;
    int shift_2 = genome_bit_i & 0b111111;
    uint32_t n_mismatches=0; // number of non-matching bits

    #define calc_n_mismatches(w1,w2) ({   \
        uint64_t bit_diff = (w1) ^ (w2);  \
        diff = ((bit_diff >> 1) | bit_diff) & 0x5555555555555555ULL; /* OR the two bits of every base - if either is 1 then count as a diff. */ \
        __builtin_popcountll (diff); }) /* note: expected to be _mm_popcnt_u64 with SSE4.2 */

    uint64_t next_w2 = *words_2; 
    if (is_forward)
        for (const uint64_t *w1=words_1, *w2=words_2; w1 < after_1; w1++, w2++) {
            uint64_t this_w2 = next_w2;
            next_w2 = *(w2+1);
            // note: if last seq_bits word is partial, we still calculate it in its entirety and correct later. However, in case we reach 
            // the end of genome and have shift_2, we will access the word beyond the end of the genome->words. counting on words being a Buffer to not cause a segfault.
            uint64_t w2_value = _bits_combined_word (this_w2, next_w2, shift_2);
            n_mismatches += calc_n_mismatches (*w1, w2_value);
        }

    else 
        for (const uint64_t *w1=words_1, *w2=words_2; w1 < after_1; w1++, w2--) {
            uint64_t this_w2 = next_w2;
            next_w2 = *(w2-1);
            uint64_t w2_revcomp = bits_revcomp_word (_bits_combined_word (next_w2, this_w2, shift_2)); 
            n_mismatches += calc_n_mismatches (*w1, w2_revcomp);
        }

    // decrease n_mismatches due to the unused part of the last seq_bits word
    if (seq_bits->nbits & 0b111111)
        n_mismatches -= __builtin_popcountll (diff & ~bitmask64 (seq_bits->nbits & 0b111111));

    return n_mismatches; 
    #undef calc_n_mismatches
}

static inline bool aligner_update_best (VBlock𐤐 vb, PosType64 gpos, 
                                        PosType64 gpos_R1, // if this is for R2: gpos of R1
                                        PosType64 gpos1,   // if splicing: main gpos of sequence
                                        Bits𐤐 seq_bits, uint32_t seq_len, bool fwd, 
                                        ConstBits𐤐 genome, PosType64 genome_nbases, uint32_t near_perfect_max_mismatches,
                                        BestAlignment *restrict best) // in/out
{
    // START_TIMER; // this has a small performance impact as it is called in a tight loop - uncomment when needed

    // note on 'N' in seq: since 'N' was encoded as 'A' in aligner_seq_to_bitmap(), we might be
    // slightly underestimating the mismatches in cases were the genome also happens to have an 'A' there. That's ok - piz sees A in the genome too. 
    // note: if rev, the genome sequence compared is still [gpos, gpos+seq_len), just revcomped
    int32_t n_matching = seq_len - aligner_count_mismatches (seq_bits, genome, gpos, fwd);
        
    // penalty for remote GPOS in 2nd pair
    uint32_t penalty = (gpos_R1 != NO_GPOS // this is R2
                     && gpos1 == NO_GPOS   // this is not the 2nd segment in a spliced alignment
                     && ABS(gpos-gpos_R1) > PAIR_MAX_DISTANCE) // the GPOS is too far from R1 GPOS
        ? NON_PAIR_PENALTY : 0;
    
    uint32_t score = n_matching - penalty; 
              
    // 32 bit arithmetic
    if (score > best->score) { 
        *best = (BestAlignment){ .gpos       = gpos,
                                 .n_matching = n_matching,
                                 .score      = score,
                                 .is_forward = fwd,
                                 .is_perfect = (n_matching == seq_len) };
        
        // note: we allow "near_perfect_max_mismatches" mismatches and we still consider the match good enough and stop looking further    
        // compared to stopping only if n_matching==seq_len, this adds about 1% to the file size, but is significantly faster 
        if (score >= seq_len - near_perfect_max_mismatches) { // we found (almost) the best possible match              
            // COPY_TIMER (aligner_update_best); 
            return true; // perfect or near perfect                                                                                                     
        }                                                                                                                       
    }        

    // COPY_TIMER (aligner_update_best); 
    return false; // not perfect (possibly because of penalty) - consider other GPOSs before deciding
}                                                                                                                               

static __attribute__((always_inline)) inline 
PosType64 aligner_collect_gpos (VBlock𐤐 vb, 
                                uint64_t s1_hash, Direction s1_is_fwd, PosType64 s1_gpos_decrement,
                                PosType64 gpos1, uint32_t seq_len, 
                                ConstBits𐤐 genome, PosType64 genome_nbases,
                                BestAlignment *restrict best)
{
    PosType64 gpos;

    // pick up gpos we previously scheduled for pre-fetching
    if (gpos_bytes == 4) {
        gpos = LTEN32 (*B32(refhash_buf, s1_hash));
        if (gpos == NO_HASH_ENT32) return NO_GPOS;
    }

    else { // gpos_bytes == 5
        gpos = U40to64_LTEN (*B40(refhash_buf, s1_hash));
        if (gpos == NO_HASH_ENT40) return NO_GPOS;
    }

    gpos -= s1_gpos_decrement;
    
    // ignore this gpos if the seq wouldn't fall completely within reference genome  
    if (__builtin_expect (gpos < 0 || gpos + seq_len >= genome_nbases, false))  
        return NO_GPOS;

    // ignore it if it has already been evaluated via previous hook.
    if (gpos == best->gpos) 
        return NO_GPOS;

    // consider splicing, only if gpos is close enough to (but not the same as) the main gpos
    if (gpos1 != NO_GPOS && 
        (gpos == gpos1 || ABS(gpos - gpos1) > MAX_SPLICE_GAP))
        return NO_GPOS;

    // prefetch max 2 cache lines (each cache line = 64 bytes = 256 bases)
    char *first_word = (char *)&genome->words[gpos / 32/*32 bases per word*/];
    char *last_word  = (char *)&genome->words[(gpos + seq_len - 1) / 32];
    
    // Locality=3 (fully cache) is best for seg_all_data_lines time (tested). Reason: if this genome 
    // sequence is the best, we will access it again shortly after to calculate mismatches

    // note: in case of REVERSE, aligner_count_mismatches accesses the last word first, so we fetch its cache line first
    __builtin_prefetch (s1_is_fwd ? first_word : last_word, 0/*read-only*/, 3/*fully cached*/);
    
    // case: beginning and end of genome sequence are on different cache lines prefetch one more 64B cache line
    if (CPU_CACHE_LINE (first_word) != CPU_CACHE_LINE (last_word))
        __builtin_prefetch (s1_is_fwd ? (first_word + CPU_CACHE_LINE_BYTES) : (last_word - CPU_CACHE_LINE_BYTES), 0 , 3);

    return gpos;
}

static void aligner_evaluate_hooks (
    VBlock𐤐 vb, 
    Bits𐤐 seq_bits, rom𐤐 seq, PosType64 seq_len, // sequence to be aligned 
    ConstBits𐤐 genome, PosType64 genome_nbases,  // genome aligned against 
    PosType64 gpos_R1,  // we are R2, this is the gpos of our counterpart in R1
    PosType64 gpos1,    // we are looking for the 2nd segment of a splice alignment, this is the gpos of the first segment
    BestAlignment *restrict best) // out
{
    START_TIMER;
    
    // in case of --fast, we check only 1/3 of the bases, and we are content with a match (not searching any further) if it 
    // has at most 10 SNPs. On our test file, this reduced the number of calls to aligner_update_best by about 4X, 
    // at the cost of the compressed file being about 11% larger
    uint32_t near_perfect_max_mismatches = (flag.fast?NEAR_PERFECT_FAST : flag.best?NEAR_PERFECT_BEST : NEAR_PERFECT_NORM);
    
    // pipeline stage 1: fetching gpos
    #define NOT_FETCHING UINT64_MAX
    uint64_t s1_hash = NOT_FETCHING; 
    bool s1_is_fwd = 0;
    PosType64 s1_gpos_decrement = 0;

    // pipeline stage 2: fetching genome words
    PosType64 s2_gpos = NO_GPOS;
    bool s2_is_fwd = 0;

    char seq_prev=0, seq_curr=0, seq_next=seq[0];
    uint8_t flushing = 0;
    uint64_t i; // declare before the loop to allow "goto" back into the loop for flushing

    // We identify locations that might contain our sequence in the genome based on kmers appearing after a hook. 
    // We do so in a 3-stage pipeline to parallelize random access RAM reads (from refhash and then genome) with processing:
    // Pipeline stage 1: Initiate prefetch of gpos from refhash[hash]. hash is calculated from the kmer starting at the base after the hook.
    // Pipeline stage 2: Initiate prefetch of the cache lines containing the genomic sequence from genome[gpos]
    // Pipeline stage 3: Compare hamming distance of the genome[gpos] sequence (already fetched) and our sequence.
    for (i=0; i < seq_len; i++) {          
        seq_prev = seq_curr;
        seq_curr = seq_next;
        seq_next = seq[i+1]; // note: no problem i+1 overflowing - this Buffer memory: in FASTQ/SAM this is vb->txt_data and BAM vb->textual_seq
        
        if ((seq_curr == HOOK   // 75% failure: first condition to fail fast 
          && seq_next != HOOK // 25% failure 
          && i < seq_len - bases_per_hash // room for the hash word: ~10% failure, but pure register test
          && (gpos1 == NO_GPOS || best->is_forward == true)) // 100% success in the non-splice aligning
     ||
            (seq_curr == HOOK_REV
          && seq_prev != HOOK_REV 
          && i >= bases_per_hash
          && (gpos1 == NO_GPOS || best->is_forward == false))) {

        flush_pipeline:
            // pipeline stage 3: calculate outcome: use genome sequence prefetched in stage 2 to calculate hamming distance and evaluate the result
            if (s2_gpos != NO_GPOS &&
                aligner_update_best (vb, s2_gpos, gpos_R1, gpos1, seq_bits, seq_len, s2_is_fwd, 
                                     genome, genome_nbases, near_perfect_max_mismatches, best)) // true if near-perfect match
                goto done; // near-perfect match, search no longer

            // pipeline stage 2: prefetch genome: use gpos fetched in stage 1, to prefetch the genomic sequence at genome[gpos]
            if (s1_hash != NOT_FETCHING) {
                s2_gpos = aligner_collect_gpos (vb, s1_hash, s1_is_fwd, s1_gpos_decrement, gpos1, seq_len, 
                                                genome, genome_nbases, best);
                s2_is_fwd = s1_is_fwd;
            }
            else
                s2_gpos = NO_GPOS;

            // pipeline stage 1: prefetch gpos: calculate hash by kmer, and start prefetching gpos from refhash[hash]
            if (!flushing) {
                uint64_t base_i;
                if (seq_curr == HOOK) {
                    s1_is_fwd = FORWARD;
                    s1_gpos_decrement = i; // after decrementing, gpos is the first base on the reference, that aligns to the first base of seq
                    base_i = i+1; // this is faster than ?: another conditional in the function call
                }
                else { 
                    s1_is_fwd = REVERSE;
                    s1_gpos_decrement = seq_len-1 - i; // after decrementing, gpos is the first base of the reference, that aligns wit the LAST base of seq
                    base_i = i-1;
                }
                s1_hash = aligner_prefetch_gpos (vb, seq_bits, base_i, s1_is_fwd);  
            }
            else
                s1_hash = NOT_FETCHING;
        }
    }

    // flush pipeline twice, to process what is currently in stage 1 and 2
    if (++flushing < 3) goto flush_pipeline; // like this, to avoid expanding inline functions again

done:
    if (gpos1 == NO_GPOS)
        COPY_TIMER (aligner_evaluate_hooks);    
    else
        COPY_TIMER (aligner_evaluate_hooks2);    
}

// returns gpos aligned with seq with M (as in CIGAR) length, containing the longest match to the reference. 
// returns false if no match found.
// note: matches that imply a negative GPOS (i.e. their beginning is aligned to before the start of the genome), aren't consisdered
static inline PosType64 aligner_best_match (
    VBlock𐤐 vb, STR𐤐(seq), PosType64 gpos_R1,
    ConstBits𐤐 genome, PosType64 genome_nbases,
    bool𐤐 is_forward, bool𐤐 is_perfect, PosType64𐤐 gpos2, bool𐤐 can_match_as_nonspliced) // out
{
    START_TIMER;

    // convert seq to a bitmap
    uint64_t seq_bits_words[roundup_bits2words64(seq_len * 2)];
    bool seq_is_all_acgt;
    Bits seq_bits = aligner_seq_to_bitmap (vb, STRa(seq), seq_bits_words, &seq_is_all_acgt);
    
    //ref_print_bases (&seq_bits, "\nseq_bits fwd", true);
    //ref_print_bases (&seq_bits, "seq_bits rev", false);
       
    BestAlignment best = { .gpos=NO_GPOS }, best2 = { .gpos=NO_GPOS }/*must be initialized even if not splicing*/;
    aligner_evaluate_hooks (vb, &seq_bits, STRa(seq), genome, genome_nbases, gpos_R1, NO_GPOS, &best);

    int match_percent = percent_match (best.n_matching, seq_len); // note: this is matching bits, not matching bases

    *is_perfect = best.is_perfect && seq_is_all_acgt; // if !seq_is_all_acgt (=seq has an 'N'), it is never a perfect match even n_matching says so, because 'N's are converted to 'A's in the seq bitmap, and might therefore inflate n_matching by matching an 'A' in the genome
    *is_forward = best.is_forward;
    *can_match_as_nonspliced = (match_percent >= MIN_MATCH_PC);                                         

    // case: spliced alignment allowed, and match is relatively short - check for another segment (might be a biological exon, a biological indel variant, or a sequencer indel error)
    if (!flag.no_splicing 
     && IN_RANGX (match_percent, SPLICE_MIN_PC, SPLICE_MAX_PC) 
     && seq_len >= SPLICE_MIN_SEQ_LEN) {
        
        best2.is_forward = best.is_forward; // spliced segments must be in the same orientation
        aligner_evaluate_hooks (vb, &seq_bits, STRa(seq), genome, genome_nbases, NO_GPOS, best.gpos, &best2);
        *gpos2 = best2.gpos;
    }

    COPY_TIMER (aligner_best_match);

    // minimum critiera to be worth segging the alignment. note for spliced: this is an initial filtering, verifying that 
    // sum of matching bits (even double-counting overlaps) matches the criteria. aligner_seg_mismatches() will re-test more tightly after junction is known.
    if (percent_match (best.n_matching + best2.n_matching, seq_len) < MIN_MATCH_PC)
        return NO_GPOS;
    else
        return best.gpos;
}


void aligner_seg_gpos_and_fwd (VBlock𐤐 vb, 
                               uint32_t seq_len, // note: less that vb->seq_len when aligning the excess of a nonbio read
                               PosType64 gpos, PosType64 gpos2, bool is_forward, 
                               uint32_t junction, // only if spliced, i.e. gpos2 is set
                               bool am_i_R2, // R2 file, or R2 read in an interleaved file (whether or not we have gpos_R1)
                               PosType64 gpos_R1, bool is_forward_R1, // used iff am_i_R2
                               PosType64𐤐 G1) // out
{
    declare_seq_contexts;

    ContextP this_strand_ctx = (am_i_R2 && segconf.is_interleaved) ? strand_r2_ctx : strand_ctx;

    if (!this_strand_ctx->local.data || (this_strand_ctx->local.size & ~3ULL) * 8 == this_strand_ctx->local.nbits)
        buf_alloc_do (vb, &this_strand_ctx->local, this_strand_ctx->local.size + sizeof(uint64_t), CTX_GROWTH, NULL, __FUNCLINE);

    buf_add_bit (&this_strand_ctx->local, gpos_R1 == NO_GPOS ? is_forward // we are not r2, or we are r2 and r1 is unaligned - just store the strand
                                                             : (is_forward == is_forward_R1)); // pair 1 is aligned - store equality, expected to he 1 in most cases
    
    // note: this schema is copied also in aligner_recon_get_gpos_and_fwd                                                               
    //
    // **FORWARD**    ————————SEG————————       ———————————PIZ———————————
    // seq:           0·····j···········L       0·····j    0···········L2   (g=gpos j=junction L=(seq_len-1) L2=L-j)
    // 1st segment ⟶ g···············g+L       G1·····    G2············   SEG: GPOS ⇐ g ; JUNCTION ⇐ j ; GAP ⇐ g2-g1
    //                ‾‾‾‾‾‾_____________                                   PIZ: G1 ⇐ GPOS ; L2 ⇐ seq_len-JUNCTION   
    // 2st segment ⟶ g2·············g2+L                                        G2 ⇐ G1+GAP+JUNCTION
    //
    // **REVERSE**    ————————SEG————————       ———————————PIZ———————————   SEG: GPOS ⇐ g+(seq_len-j)  
    // seq:           0·····j···········L       0·····j    0···········L2        JUNCTION ⇐ j 
    // 1st segment ⟵ g+L···············g       ·····G1    ············G2        GAP ⇐ g2-g1
    //                ‾‾‾‾‾‾_____________                                   PIZ: G1 ⇐ GPOS ; L2 ⇐ seq_len-JUNCTION
    // 2st segment ⟵ g2+L·············g2                                        G2 ⇐ G1+GAP-L2 = G1+GAP+JUNCTION-seq_len
    
    // revcomp case of splicing -  GPOS ⇐ g+(seq_len-j) (see above)
    *G1 = (gpos2 != NO_GPOS && !is_forward) ? (gpos + seq_len - junction) 
                                            : gpos;
    if (am_i_R2) 
        fastq_seg_gpos_R2 (vb, gpos_R1, *G1, is_forward);

    // store the GPOS in local if its not a 2nd pair, or if it is, but the delta is not small enough
    else 
        seg_integer (VB, gpos_ctx, *G1, false, 0);
    
    // if there is a spliced alignment, store its gpos2 as a delta
    if (gpos2 != NO_GPOS)           
        // note: (gpos2 - gpos) is the gap: the number of bases in the "intron" between the two spliced sub-sequences. can be negative.
        seg_integer (vb, gpos_gap_ctx, gpos2 - gpos, false, 0);

    // note regarding interleaved: we always store last_value in gpos/strand (never in gpos_r2/strand_r2)
    gpos_ctx->last_value.i = gpos; // note: gpos, not G1
    strand_ctx->last_value.i = is_forward;
}

// to find out whether 1 or ref2 are first, we compare the first few bases - at least 75% match
static thool who_is_first (STR𐤐(seq), rom ref1, rom ref2, bool test_left)
{
    #define WHO_IS_FIRST_N_SEQ 16

    int matches1=0, matches2=0;
    int first_i = test_left ? 0 : seq_len - WHO_IS_FIRST_N_SEQ;
    
    for (int i=first_i; i < first_i + WHO_IS_FIRST_N_SEQ; i++) {
        if (seq[i] == ref1[i]) matches1++;
        if (seq[i] == ref2[i]) matches2++;
    }

    #define MATCH_THREASH (3*WHO_IS_FIRST_N_SEQ / 4)
    thool result = (matches1 > MATCH_THREASH && matches2 < MATCH_THREASH) ? yes
                 : (matches2 > MATCH_THREASH && matches1 < MATCH_THREASH) ? no 
                 :                                                          unknown;

    // if we're testing the right end, then the result is reversed
    if (!test_left && result != unknown)
        result = !result;

    // if we're testing the left end, and can't reach a conclusion - test the right end (recursively)
    if (test_left && result == unknown)
        return who_is_first (STRa(seq), ref1, ref2, false);
    else
        return result;
}

// assuming the beginning of seq aligns with ref1, and the end with ref2, find the ideal junction index
static void aligner_get_junction_do (STR𐤐(seq), rom ref1, rom ref2, int *restrict max_matches, int *restrict junction)
{
    int my_max_matches   = 0;
    int my_junction = seq_len;

    // count matches against ref1 alone (i.e. junction point at end of seq - 0 bases from ref2)
    for (int i=0; i < seq_len; i++) 
        if (seq[i] == ref1[i]) my_max_matches++;

    // now replace one base at the of ref1 with a base from ref2 and see if that improves things
    int matches = my_max_matches; // initialize - matches when it is all ref1 
    for (int i=seq_len-1; i >= 0; i--) {
        if (seq[i] == ref1[i]) matches--; // remove match due to ref1 in base i
        if (seq[i] == ref2[i]) matches++; // add a match if it matches ref 2 in base i
        
        // now matches has matches vs ref1 until i-1 PLUS matches vs ref2 for bases starting i
        if (matches > my_max_matches) {
            my_max_matches = matches;
            my_junction = i;
        }
    }

    *max_matches = my_max_matches;
    *junction    = my_junction;
}

// finds index in seq which minimizes mismatches vs ref and ref2
static bool aligner_get_junction (VBlock𐤐 vb, STR𐤐(seq), rom ref1, rom ref2, // note: ref1,2 might be aliased pointers
                                  uint32_t *junction, bool *ref1_then_ref2)    // out
{          
    START_TIMER;

    _Static_assert (SPLICE_MIN_SEQ_LEN >= WHO_IS_FIRST_N_SEQ, "WHO_IS_FIRST_N_SEQ too long");

    // get matches vs ref1 (the ref with more bit matches), if splicing is not applied
    int matches_without_splicing=0; 
        for (int i=0; i<seq_len; i++) 
            if (ref1[i] == seq[i]) 
                matches_without_splicing++; 

    // ref1 is a better match that ref2, but it is not necessary appears first in spicing.
    thool ref1_is_first = who_is_first (STRa(seq), ref1, ref2, true); // see if we can determine from the first 16 bases
    
    int max_matches_ref1_first=0, max_matches_ref2_first=0, junction_if_ref1_first, junction_if_ref2_first=0;

    // ref1 *might* be first - get max_matches assuming it is
    if (ref1_is_first != no) {
        aligner_get_junction_do (STRa(seq), ref1, ref2, &max_matches_ref1_first, &junction_if_ref1_first);

        // if we get at least 95% matching bases when ref1 is first, then our doubts are over
        if (ref1_is_first == unknown && percent (max_matches_ref1_first, seq_len) >= 95)
            ref1_is_first = yes;
    }

    // ref2 *might* be first - get max_matches assuming it is (possibly after already assuming ref1 is)
    if (ref1_is_first != yes) { 
        SWAP (ref1, ref2);
        aligner_get_junction_do (STRa(seq), ref1, ref2, &max_matches_ref2_first, &junction_if_ref2_first);
    }
 
    // check if splicing helps enough to be worth while
    if (MAX_(max_matches_ref1_first, max_matches_ref2_first) - matches_without_splicing < MIN_SPLICE_MATCH_CONTRIBTION) {
        COPY_TIMER (aligner_get_junction);
        return false; // no splicing
    }

    *ref1_then_ref2 = (max_matches_ref1_first >= max_matches_ref2_first);
    *junction =  *ref1_then_ref2 ? junction_if_ref1_first : junction_if_ref2_first;
    
    COPY_TIMER (aligner_get_junction);
    return true; // splice it away!
}

static char *alloc_ref (VBlock𐤐 vb, STRc(ref_space), uint32_t ref_len)
{
    ASSERTNOTINUSE (vb->scratch);

    if (ref_len <= ref_space_len - 8) // leave 4 bytes before and after for overwrites by ref_get_textual_seq
        return ref_space + 4;

    else {
        buf_alloc_exact (vb, vb->scratch, ref_len, char, "scratch");
        return B1STc(vb->scratch);
    }
}

// get textual reference starting at the lower of gpos, gpos2, and containing seq_len starting for each gpos and gpos2
static inline void aligner_get_textual_ref (VBlock𐤐 vb, 
                                            PosType64 gpos, PosType64 gpos2, uint32_t seq_len, bool is_forward,
                                            STRc (ref_space), // use if possible, or scratch if not
                                            rom *ref, rom *ref2)
{    
    if (gpos2 == NO_GPOS) gpos2 = gpos;

    uint32_t abs_gpos_gap = ABS(gpos2 - gpos);
    
    // case: single sequence
    if (abs_gpos_gap <= seq_len) {
        uint32_t ref_len = seq_len + abs_gpos_gap;
        char *ref_data = alloc_ref (vb, STRa(ref_space), ref_len);

        PosType64 min_gpos = MIN_(gpos, gpos2);

        // note: if reverse, seq generated is still between [min_gpos, (min_gpos + ref_len - 1)], just revcomped 
        ref_get_textual_seq (vb, min_gpos, ref_data, ref_len, !is_forward);

        *ref  = (gpos  == min_gpos) ? ref_data : (ref_data + abs_gpos_gap);
        *ref2 = (gpos2 == min_gpos) ? ref_data : (ref_data + abs_gpos_gap);

        // note: if reverse, entire ref_data is revcomped, so ref1 last base is gpos
        if (!is_forward) SWAP (*ref, *ref2);
    } 

    // two disjoint sequences, each
    else {
        char *ref_data = alloc_ref (vb, STRa(ref_space), seq_len * 2);

        // note: if reverse, ref and ref2 are revcomped in-place independently
        ref_get_textual_seq (vb, gpos, ref_data, seq_len, !is_forward);
        *ref = ref_data;

        ref_get_textual_seq (vb, gpos2, &ref_data[seq_len], seq_len, !is_forward);
        *ref2 = ref_data + seq_len;
    }
}

// segs SQBITMAP(local), SEQMIS_*, JUNCTION. Returns true is is_perfect (no mismatches)
static bool aligner_seg_mismatches (VBlock𐤐 vb, STR𐤐(seq), 
                                    PosType64𐤐 gpos, PosType64𐤐 gpos2, // might be swapped
                                    bool is_forward, bool can_match_as_nonspliced,
                                    uint32_t *restrict junction) // out
{
    START_TIMER;

    declare_seq_contexts;

    Bits𐤐 bitmap = buf_alloc_bits (vb, &bitmap_ctx->local, seq_len, vb->lines.len32 / 16 * segconf.std_seq_len, 
                                   SET/*initialize to "no mismatches"*/, CTX_GROWTH, C_LOCAL); 

    for (int i=0; i < 4; i++)
        buf_alloc (vb, &seqmis_ctx[i].local, seq_len, 64 KB, char, CTX_GROWTH, NULL); 

    // get textual reference segment (possibly rev-comped)
    char ref_space[4 KB]; // not too big, so "scratch" code path gets milage
    rom ref, ref2; // note: may be aliased
    aligner_get_textual_ref (vb, *gpos, *gpos2, seq_len, is_forward, ref_space, sizeof(ref_space), &ref, &ref2);

#ifdef DEBUG // DEBUG only due to performance
    if (flag.debug_aligner)
        iprintf ("%s: mismatches muxed by these ref bases: gpos=%-10"PRIu64"%s %s seq_len=%-3u: ", 
                 LN_NAME, *gpos, cond_int (*gpos2 != NO_GPOS, " gpos2=", *gpos2), 
                 is_forward ? "FWD" : "REV", seq_len);
#endif

    *junction = seq_len;
    if (*gpos2 != NO_GPOS) {
        bool ref1_then_ref2;

        if (aligner_get_junction (vb, STRa(seq), ref, ref2, junction, &ref1_then_ref2)) {
            seg_integer (vb, junction_ctx, *junction, false, 0);

            if (!ref1_then_ref2) {  // ref2 is first - swap
                SWAP (ref, ref2); 
                SWAP (*gpos, *gpos2);
            }
        }

        // case: turns out that ref2 doesn't contribute much to matches. discard it.
        else {
            *gpos2 = NO_GPOS;

            // if, absent spliced alignment, the alignment now falls outside the success criteria, discard it
            if (!can_match_as_nonspliced) {
                *gpos = NO_GPOS; // we're better off segging as verbatim
                return false; 
            }
        }
    }

    // handle mismatches
    uint64_t bit_i = bitmap_ctx->next_local;  // automatic variable for performance

    rom r = ref;
    bool is_perfect = true; // note: relevant for spliced - we identify perfection here (for non-spliced, this function is not called if is_perfect)
    for (uint32_t i=0; i < seq_len; i++, r++, seq++, bit_i++) {
        if (i == *junction) {
            r = &ref2[i]; // starting at the junction, compare to ref2 instead
#ifdef DEBUG 
            if (flag.debug_aligner) iprint0 ("+"); // separator between two spliced segments
#endif
        }
        if (*seq != *r) {
#ifdef DEBUG 
            if (flag.debug_aligner) iprintf ("%c", is_forward ? *r : COMPLEM[(uint8_t)*r]); 
#endif
            bits_clear (bitmap, bit_i); 
            BNXTc (seqmis_ctx[nuke_encode_dir (*r, is_forward)].local) = *seq;
            
            is_perfect = false;
        }
    }

#ifdef DEBUG 
    if (flag.debug_aligner) iprint_newline();
#endif

    if (!is_perfect)
        bitmap_ctx->next_local += seq_len;

    buf_free (vb->scratch);

    COPY_TIMER (aligner_seg_mismatches);
    return is_perfect;
}

// seg with --REFERENCE: set ref_is_set to store the required subset of the reference that needs storing
static void aligner_set_is_set (PosType64 gpos, PosType64 gpos2/*if spliced*/, uint32_t seq_len, uint32_t junction,/*if set*/ bool is_forward)
{
    // update is_set    
    if (IS_REF_EXT_STORE) {
        if (gpos2 == NO_GPOS) {
            RefLock lock = ref_lock (gpos, seq_len); // lock region of reference to protect is_set
            ref_set_genome_is_used (gpos, seq_len);  // this region of the reference is used (in case we want to store it with REF_EXT_STORE)
            ref_unlock (&lock);
        }

        else { // spliced
            // lock all of ref, ref2 and the gap between them (if there is a gap)
            RefLock lock = ref_lock (MIN_(gpos, gpos2), seq_len + ABS(gpos - gpos2)); 

            if (is_forward) {
                // seq:           0·····j···········L  (g=gpos j=junction L=(seq_len-1))
                // 1st segment ⟶ g···············g+L       
                //                ‾‾‾‾‾‾_____________
                // 2st segment ⟶ g2·············g2+L
                ref_set_genome_is_used (gpos, junction); 
                ref_set_genome_is_used (gpos2 + junction, seq_len - junction); 
            }            
            else {
                // seq:           0·····j···········L
                // 1st segment ⟵ g+L···············g      
                //                ‾‾‾‾‾‾_____________
                // 2st segment ⟵ g2+L·············g2
                ref_set_genome_is_used (gpos + (seq_len - junction), junction); 
                ref_set_genome_is_used (gpos2, seq_len - junction); 
            }

            ref_unlock (&lock);
        }
    }
}

MappingType aligner_seg_seq (VBlock𐤐 vb, STR𐤐(seq),  
                             bool am_i_R2, // R2 file, or R2 read in an interleaved file (whether or not we have gpos_R1)
                             PosType64 gpos_R1, bool is_forward_R1)
{
    START_TIMER;
    declare_seq_contexts;
    
    ConstBitsP genome;
    PosType64 genome_nbases;
    
    ref_get_genome (&genome, &genome_nbases);

    bool is_forward=false, is_perfect=false, is_perfect_spliced=false, can_match_as_nonspliced=true;

    PosType64 gpos2 = NO_GPOS; 
    PosType64 gpos =  
        (seq_len > MAX_SHORT_READ_LEN) ? NO_GPOS
      : refhash_is_flat                ? aligner_best_match    (VB, STRa(seq), gpos_R1, genome, genome_nbases, &is_forward, &is_perfect, &gpos2, &can_match_as_nonspliced) 
      :                                  aligner_best_match_LAYERED (VB, STRa(seq), gpos_R1, genome, genome_nbases, &is_forward, &is_perfect);
    uint32_t junction = 0;

    if (gpos == NO_GPOS || gpos > genome_nbases - seq_len || gpos > MAX_ALIGNER_GPOS - seq_len || gpos < 0/*never happens*/) {
        gpos_ctx->last_value.i = NO_GPOS;

no_mapping:
        COPY_TIMER (aligner_seg_seq);
        return MAPPING_NO_MAPPING;
    }

    if (gpos2 != NO_GPOS && (gpos2 > genome_nbases - seq_len || gpos2 > MAX_ALIGNER_GPOS - seq_len || gpos2 < 0))
         gpos2 = NO_GPOS; // if 2nd segment in a spliced alignment has bad gpos2, cancel it

    // skip checking mismatches if perfect (note: this works only for non-spliced)         
    if (!is_perfect) {
        // segs SQBITMAP(local), SEQMIS_*, JUNCTION (may set gpos / gpos2 to NO_GPOS if failed)
        // NOTE ON GPOS: IN: gpos is of segment with most matches, gpos2 is second most matches (gpos2 exists only if splicing)
        //              OUT: gpos is the segment matching the beginning of seq, and gpos2 matches the end of seq
        //                   i.e., if the ref starting at gpos2 turns out to be first, we switch gpos ⇔ gpos2
        is_perfect_spliced = aligner_seg_mismatches (vb, STRa(seq), &gpos, &gpos2, is_forward, can_match_as_nonspliced, &junction);
                
        if (gpos == NO_GPOS) goto no_mapping;

        if (is_perfect_spliced) 
            vb->num_perfect_spliced++;
    }

    if (is_perfect || is_perfect_spliced) {
        vb->num_aligned_perfect++; // for stats: no mismatches
        
        if (am_i_R2)
            vb->num_perfect_r2++; // use this to verify the numer of perfect in R2 makes sense
    }

    // segs STRAND, GPOS, GPOS_DELTA, GPOS_R2, GPOS_GAP
    PosType64 G1;
    aligner_seg_gpos_and_fwd (vb, seq_len, gpos, gpos2, is_forward, junction, am_i_R2, gpos_R1, is_forward_R1, &G1);

    if (IS_REF_EXT_STORE)
        aligner_set_is_set (gpos, gpos2, seq_len, junction, is_forward);

    vb->num_aligned++; // for stats
    if (gpos2 != NO_GPOS) vb->num_aligned_spliced++;

    if (flag.show_aligner)
        iprintf ("%s: gpos=%-10"PRId64"%s%s%s%s%s%s %s %s\n", LN_NAME, 
                 gpos, 
                 cond_int (ctx_has_value_in_line_(vb, gpos_Δ_ctx), " gpos_Δ=", gpos_Δ_ctx->last_value.i), // set in fastq_seg_gpos_R2
                 cond_int (gpos2 != NO_GPOS, " gpos2=", gpos2), 
                 cond_int (gpos2 != NO_GPOS, " gpos_gap=", gpos2 - gpos), 
                 cond_int (gpos2 != NO_GPOS, " G1=", G1), 
                 cond_int (gpos2 != NO_GPOS, " G2=", gpos2 + (is_forward ? junction : 0)), 
                 cond_int (gpos2 != NO_GPOS, " junc=", junction),
                 is_forward ? "FWD" : "REV", is_perfect ? "PERFECT" : "");

    COPY_TIMER (aligner_seg_seq);

    return is_perfect_spliced ? MAPPING_PERFECT_SPLICED
         : (gpos2 != NO_GPOS) ? MAPPING_SPLICED
         : is_perfect         ? MAPPING_PERFECT
         :                      MAPPING_ALIGNED; 
}

void aligner_recon_get_gpos_and_fwd (VBlock𐤐 vb, bool am_i_R2, 
                                     bool spliced_2nd_segment, 
                                     PosType64𐤐 gpos, bool𐤐 is_forward) // out
{    
    declare_seq_contexts;
    bool is_interleaved_r2 = (segconf.is_interleaved && VB_DT(FASTQ) && vb->line_i % 2); // FASTQ file or FASTQ component of a Deep file

    // spliced alignment: (this schema is copied also in aligner_seg_gpos_and_fwd)
    //
    // **FORWARD**    ————————SEG————————       ———————————PIZ———————————
    // seq:           0·····j···········L       0·····j    0···········L2   (g=gpos j=junction L=(seq_len-1) L2=L-j)
    // 1st segment ⟶ g···············g+L       G1·····    G2············   SEG: GPOS ⇐ g ; JUNCTION ⇐ j ; GAP ⇐ g2-g1
    //                ‾‾‾‾‾‾_____________                                   PIZ: G1 ⇐ GPOS ; L2 ⇐ seq_len-JUNCTION   
    // 2st segment ⟶ g2·············g2+L                                        G2 ⇐ G1+GAP+JUNCTION
    //
    // **REVERSE**    ————————SEG————————       ———————————PIZ———————————   SEG: GPOS ⇐ g+(seq_len-j)  
    // seq:           0·····j···········L       0·····j    0···········L2        JUNCTION ⇐ j 
    // 1st segment ⟵ g+L···············g       ·····G1    ············G2        GAP ⇐ g2-g1
    //                ‾‾‾‾‾‾_____________                                   PIZ: G1 ⇐ GPOS ; L2 ⇐ seq_len-JUNCTION
    // 2st segment ⟵ g2+L·············g2                                        G2 ⇐ G1+GAP-L2 = G1+GAP+JUNCTION-seq_len

    // second segment in a spliced alignment: gpos is delta vs first segment's gpos, and is_forward is the same
    if (spliced_2nd_segment) {
        *is_forward  = strand_ctx->last_value.i;   // from first segment (this function)
        int64_t G1   = gpos_ctx->last_value.i;     // from first segment (this function)
        int64_t jnct = junction_ctx->last_value.i; // from fastq_recon_aligned_SEQ 
        int64_t gap  = reconstruct_from_local_int (VB, gpos_gap_ctx, 0, false); 
        *gpos        = G1 + gap + jnct - (*is_forward ? 0 : vb->seq_len);

        ctx_set_last_value (vb, gpos_gap_ctx, gap);
    }

    // first file of a pair (R1) or a non-pair fastq or sam
    else if (!am_i_R2 && !is_interleaved_r2) {
        *gpos = reconstruct_from_local_int (vb, gpos_ctx, 0, RECON_OFF);        
        *is_forward = NEXTLOCALBIT (strand_ctx);  
        ctx_set_last_value (vb, strand_ctx, (int64_t)*is_forward);

        if (segconf.is_interleaved) bitmap_ctx->r1_is_aligned = PAIR1_ALIGNED;      
    }

    // 2nd file of a pair (R2) 
    else if (am_i_R2) {
        *is_forward = fastq_piz_get_r2_is_forward (vb); // MUST be called before gpos reconstruction as it inquires GPOS.localR1.next
        ctx_set_last_value (vb, strand_ctx, (int64_t)*is_forward);

        reconstruct_from_ctx (vb, gpos_ctx->did_i, 0, RECON_OFF); // calls fastq_special_PAIR2_GPOS
        *gpos = gpos_ctx->last_value.i;
    }
    
    // R2 in a single-file interleaved FASTQ (v15)
    else {
        *is_forward = fastq_piz_get_interleaved_r2_is_forward (vb);
        ctx_set_last_value (vb, strand_ctx, (int64_t)*is_forward);

        reconstruct_from_ctx (vb, gpos_r2_ctx->did_i, 0, RECON_OFF); // calls fastq_special_PAIR2_GPOS
        *gpos = gpos_r2_ctx->last_value.i;
    }

    // backcomp: early versions stored NO_GPOS which use to be 0xfffffffff in the GPOS context
    if (*gpos == NO_HASH_ENT32 && gpos_bytes == 4) 
        *gpos = NO_GPOS; 

    ctx_set_last_value (vb, gpos_ctx, *gpos); 
}

// PIZ: SEQ reconstruction - only for reads compressed with the aligner
void aligner_reconstruct_seq (VBlockP vb, 
                              uint32_t seq_len, bool am_i_R2, bool is_spliced_seg2, bool is_perfect_alignment, 
                              ReconType reconstruct,
                              int max_deep_mismatches,        // length of mismatch_base and mismatch_offset arrays  
                              char *mismatch_base,       // optional out
                              uint32_t *mismatch_offset, // optional out
                              uint32_t *num_mismatches)  // optional out
{
    START_TIMER;
    declare_seq_contexts;

    if (!bitmap_ctx->is_loaded) return; // if case we need to skip the SEQ field (for the entire file)
    
    if (is_perfect_alignment || buf_is_alloc (&bitmap_ctx->local)) { // not all non-ref

        bool is_forward;
        PosType64 gpos;
        aligner_recon_get_gpos_and_fwd (vb, am_i_R2, is_spliced_seg2, &gpos, &is_forward);

        if (flag.show_aligner)
            iprintf ("%s: gpos=%-10"PRId64" %s %s\n", 
                     LN_NAME, gpos, is_forward ? "FWD" : "REV", is_perfect_alignment ? "PERFECT" : "");

        if (gpos == NO_GPOS) { 
            bitmap_ctx->next_local += seq_len; // up to v13 ; not sure if this still can happen since v14
            goto all_nonref;
        }

        // reconstruct reference sequence
        ref_get_textual_seq (vb, gpos, BAFTtxt, seq_len, !is_forward);

        // non-perfect alignment: update the mismatches
        if (!is_perfect_alignment) {
            
#ifdef DEBUG // DEBUG only due to performance
            decl_acgt_decode;
            if (flag.debug_aligner)
                iprintf ("%s: mismatches muxed by these ref bases: gpos=%-10"PRIu64" %s seq_len=%-3u: %s", 
                        LN_NAME, gpos, is_forward ? "FWD" : "REV", vb->seq_len,
                        is_spliced_seg2 ? "+" : "");
#endif
            char *recon = BAFTtxt;
            for (uint32_t i=0; i < seq_len; i++) {
                uint32_t matches = bits_get_run ((BitsP)&bitmap_ctx->local, bitmap_ctx->next_local, seq_len - i);

                bool has_mismatch = (i + matches < seq_len);
                recon += matches;
                bitmap_ctx->next_local += matches + has_mismatch;
                i += matches - 1 + has_mismatch;

                // case: we have a mismatch: overwrite base copied from reference with mismatch
                if (has_mismatch) {
                    uint8_t base_2bit = nuke_encode_dir (*recon, is_forward); // if is_forward, *recon (the reconstructed ref base) is revcomped - we complement it back to its original

#ifdef DEBUG 
                    if (flag.debug_aligner) iprintf ("%c", acgt_decoder[base_2bit]);
#endif
                    
                    ContextP mis_ctx = VER(14) ? &seqmis_ctx[base_2bit] : nonref_ctx;
                    ASSPIZ (mis_ctx->next_local < mis_ctx->local.len32, NEXT_ERRFMT". "_TIP"genozip-debug --debug-aligner, compare seg and piz output", 
                            __FUNCTION__, mis_ctx->tag_name, mis_ctx->next_local, 1, mis_ctx->local.len32); 

                    *recon++ = *Bc(mis_ctx->local, mis_ctx->next_local++);

                    if (num_mismatches) {
                        if (*num_mismatches <= max_deep_mismatches) {
                            mismatch_base  [*num_mismatches] = *(recon-1);
                            mismatch_offset[*num_mismatches] = i;
                        }
                        (*num_mismatches)++;
                    }
                }
            }
#ifdef DEBUG 
            if (flag.debug_aligner) iprint_newline();
#endif
        }

        if (reconstruct) Ltxt += seq_len;
    }

    else all_nonref: { 
        if (flag.show_aligner)
            iprintf ("%s: verbatim\n", LN_NAME);

        ASSPIZ (nonref_ctx->next_local + seq_len <= nonref_ctx->local.len32, "NONREF exhausted: next_local=%u + seq_len=%u > local.len=%u",
                nonref_ctx->next_local, seq_len, nonref_ctx->local.len32);

        if (reconstruct) RECONSTRUCT (Bc (nonref_ctx->local, nonref_ctx->next_local), seq_len);
        nonref_ctx->next_local += seq_len;
    }

    ASSPIZ (nonref_ctx->next_local <= nonref_ctx->local.len, "nonref_ctx->next_local=%u is out of range: nonref_ctx->local.len=%u",
            nonref_ctx->next_local, nonref_ctx->local.len32);

    COPY_TIMER (aligner_reconstruct_seq);
}
