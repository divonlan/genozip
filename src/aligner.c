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
#define MIN_MATCH_PC       73  // fail alignment if less than 73% of bits match
#define SPLICE_MIN_PC      50  // min and max bit match percent in which we search for splicing
#define SPLICE_MAX_PC      90  
#define SPLICE_MIN_SEQ_LEN 32  // no splicing for sequences shorter than this - not worth it
#define MAX_SPLICE_GAP     (32 KB - 1)  // max allowed gap beween gpos and gpos2
#define MIN_SPLICE_MATCH_CONTRIBTION 10 // for splicing to compress better, the spliced alignment needs to have this many more matching bases 

// pairing
#define PAIR_MAX_DISTANCE  500 // if pair2 is more distant from pair1 than this, it will get a penalty
#define NON_PAIR_PENALTY   17  // penalty reduction in matching bits, for a pair2 gpos that is too distant from pair1 gpos

typedef enum { NOT_FOUND=-1, REVERSE=0, FORWARD=1 } Direction;

static inline uint32_t nuke_encode_dir (char c, bool is_forward) 
{ 
    return is_forward ? nuke_encode(c) : nuke_encode_comp(c); 
}

static inline int percent_match (uint32_t match_n_bits, uint32_t seq_len)
{
    return match_n_bits * 50 / seq_len; 
}

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

// Foward example: If seq is: G-AGGGCT  (G is the hook)  -- matches reference AGGGCT       - function returns 110110101000 (A=00 is the LSb)
// Reverse       : If seq is: CGCCCT-C  (C is the hook)  -- also matches reference AGGGCT  - function returns 110110101000 - the same
// calculates a refhash word from 14 nucleotides following a 'G' (only last G in a sequenece of GGGG...)
static inline PosType64 aligner_get_gpos_FLAT (VBlock𐤐 vb, Bits𐤐 seq, uint64_t base_i, bool is_forward)                                              
{   
    uint64_t kmer = aligner_get_kmer (seq, base_i*2, is_forward); 
    uint64_t hash = fibonacci_hash (kmer);

    if (gpos_bytes == 4) {
        PosType64 gpos = LTEN32 (*B32(refhash_buf, hash));
        return (gpos == NO_HASH_ENT32) ? NO_GPOS : gpos;
    }

    else { // gpos_bytes == 5
        PosType64 gpos = U40to64_LTEN (*B40(refhash_buf, hash));
        return (gpos == NO_HASH_ENT40) ? NO_GPOS : gpos;
    }
}

// used with reference files up to 15.0.80
#define layer_bitmask ((uint64_t[NUM_LAYERS]){ 0x0fffffffULL, 0x07ffffffULL, 0x03ffffffULL, 0x01ffffffULL })  // LSb set for 28,27,26,25 bits
static inline PosType64 aligner_get_gpos_LAYERED (VBlockP vb, BitsP seq, uint64_t base_i, uint64_t *kmer, bool is_forward)                                              
{   
    *kmer = aligner_get_kmer (seq, base_i*2, is_forward); // note: LAYERED is 28 bits per kmer, so uses 32bit variables

    // prefetch cache for all layers of refhash (this improves aligner performance by ~12% and genozip performance on a simple FASTQ compression by ~8%)
    for (unsigned layer_i=0; layer_i < NUM_LAYERS; layer_i++) 
        __builtin_prefetch (B32(refhash_buf, LAYER_START[layer_i] + (*kmer & layer_bitmask[layer_i])), 0/*read-only*/, 1/*best option, empirically*/);
        
    // Performance note: 50% of the aligner time is taken by memory latency - looking up from refhash_buf - ~0.28 lookups every base in the sequence.
    PosType64 gpos = (PosType64)BGEN32(*B32(refhash_buf, (*kmer & layer_bitmask[0])));
    return (gpos == NO_HASH_ENT32) ? NO_GPOS : gpos;
}

// converts a string sequence to a 2-bit bitmap
static inline Bits aligner_seq_to_bitmap (VBlock𐤐 vb, rom𐤐 seq, uint64_t seq_len, 
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
    
        if (__builtin_expect (encoding == 4, false)) { // not A, C, G or T - usually N
            *seq_is_all_acgt = false;
            encoding = 0;    // arbitrarily convert 4 (any non-ACGT is 4) to 0 ('A')
        }
    
        bits_assign2 (&seq_bits, (base_i << 1), encoding);
    }

    bits_clear_excess_bits_in_top_word (&seq_bits, false); // bc bitmap_words is uninitialized

    COPY_TIMER (aligner_seq_to_bitmap);

    return seq_bits;
}

typedef struct {
    PosType64 gpos;
    int32_t match_n_bits; // matching bits per hamming
    int32_t score;        // matching_bits after applying penalty if needed (might be negative)
    bool is_forward; 
    bool is_perfect;
} BestAlignment;

static inline bool aligner_update_best (VBlock𐤐 vb, PosType64 gpos, 
                                        PosType64 gpos_R1, // if this is for R2: gpos of R1
                                        PosType64 gpos1,   // if splicing: main gpos of sequence
                                        Bits𐤐 seq_bits, uint32_t seq_len, bool fwd, 
                                        ConstBits𐤐 genome, PosType64 genome_nbases, uint32_t max_snps_for_near_perfect,
                                        BestAlignment *restrict best) // in/out
{
    // START_TIMER; // this has a small performance impact as it is called in a tight loop - uncomment when needed

    if (gpos == best->gpos) goto no_update; // a previous hook yielded the same gpos (which is not "near perfect") - no need to test this gpos again

    // consider splicing, only if gpos is close enough to (but not the same as) the main gpos
    if (gpos1 != NO_GPOS && 
        (gpos == gpos1 || ABS(gpos - gpos1) > MAX_SPLICE_GAP))
        goto no_update; 

    // note on 'N' in seq: since 'N' was encoded as 'A' in aligner_seq_to_bitmap(), we might be
    // slightly underestimating the distance in cases were the genome also happens to have an 'A' there. That's ok. 
    // note: if rev, the genome sequence compared is still [gpos, gpos+seq_len), just revcomped
    int32_t distance = bits_hamming_distance (seq_bits, genome, 2 * gpos, fwd);
    int32_t match_n_bits = (uint32_t)seq_bits->nbits - distance;     
        
    // penalty for remote GPOS in 2nd pair
    uint32_t penalty = (gpos_R1 != NO_GPOS // this is R2
                     && gpos1 == NO_GPOS   // this is not the 2nd segment in a spliced alignment
                     && ABS(gpos-gpos_R1) > PAIR_MAX_DISTANCE) // the GPOS is too far from R1 GPOS
        ? NON_PAIR_PENALTY : 0;
    
    uint32_t score = match_n_bits - penalty; 
              
    // 32 bit arithmetic
    if (score > best->score) { 
        *best = (BestAlignment){ .gpos         = gpos,
                                 .match_n_bits = match_n_bits,
                                 .score        = score,
                                 .is_forward   = fwd,
                                 .is_perfect   = (match_n_bits == seq_len * 2) };
        
        // note: we allow "max_snps_for_near_perfect" snps and we still consider the match good enough and stop looking further    
        // compared to stopping only if match_n_bits==seq_len, this adds about 1% to the file size, but is significantly faster 
        if (score >= (seq_len - max_snps_for_near_perfect) * 2) { // we found (almost) the best possible match              
            // COPY_TIMER (aligner_update_best); 
            return true; // perfect or near perfect                                                                                                     
        }                                                                                                                       
    }        

no_update:
    // COPY_TIMER (aligner_update_best); 
    return false; // not perfect (possibly because of penalty) - consider other GPOSs before deciding
}                                                                                                                               

static void aligner_best_match_FLAT_search (
    VBlock𐤐 vb, 
    Bits𐤐 seq_bits, rom𐤐 seq, PosType64 seq_len, // sequence to be aligned 
    ConstBits𐤐 genome, PosType64 genome_nbases,  // genome aligned against 
    PosType64 gpos_R1,  // we are R2, this is the gpos of our counterpart in R1
    PosType64 gpos1,    // we are looking for the 2nd segment of a splice alignment, this is the gpos of the first segment
    BestAlignment *restrict best) // out
{
    START_TIMER;
    
    PosType64 gpos = NO_GPOS; // match not found yet

    // in case of --fast, we check only 1/3 of the bases, and we are content with a match (not searching any further) if it 
    // has at most 10 SNPs. On our test file, this reduced the number of calls to aligner_update_best by about 4X, 
    // at the cost of the compressed file being about 11% larger
    uint64_t density = (flag.fast ? 3 : 1);
    uint32_t max_snps_for_near_perfect = (flag.fast ? 10 : 2); // note: this is only approximately SNPs as we actually measure mismatched bits, not bases

    // we search - checking both forward hooks and reverse hooks, we check only the first layer for now
    for (uint64_t i=0; i < seq_len; i += density) {          
        Direction found = NOT_FOUND;

        if ((gpos1 == NO_GPOS || best->is_forward == true) 
         && __builtin_expect (i < seq_len - bases_per_hash, true/*probability 93%*/) // room for the hash word
         && __builtin_expect (seq[i] == HOOK,  false/*probability 75%*/) 
         && __builtin_expect (seq[i+1] != HOOK, true/*probability 75%*/) // take the G - if there is a homopolymer GGGG... take the last one
         && (gpos = aligner_get_gpos_FLAT (vb, seq_bits, i+1, FORWARD)) != NO_GPOS) { 

            gpos -= i; // gpos is the first base on the reference, that aligns to the first base of seq
            found = FORWARD;
        }

        // note: if the previous condition is true, then seq[i]==HOOK and no need to check this condition. hence the "else"
        else if ((gpos1 == NO_GPOS || best->is_forward == false)
         && __builtin_expect (i >= bases_per_hash, true) 
         && __builtin_expect (seq[i] == HOOK_REV,  false) 
         && __builtin_expect (seq[i-1] != HOOK_REV, true)  
         && (gpos = aligner_get_gpos_FLAT (vb, seq_bits, i-1, REVERSE)) != NO_GPOS) {

            gpos -= seq_len-1 - i; // gpos is the first base of the reference, that aligns wit the LAST base of seq
            found = REVERSE;
        }

        if (found != NOT_FOUND && 
            __builtin_expect (gpos >= 0 && gpos + seq_len < genome_nbases, true)) { // ignore this gpos if the seq wouldn't fall completely within reference genome
            
            if (aligner_update_best (vb, gpos, gpos_R1, gpos1, seq_bits, seq_len, found, 
                                     genome, genome_nbases, max_snps_for_near_perfect, best)) {
                return; // near-perfect match, search no longer
            }
        }
    }

    if (gpos1 == NO_GPOS)
        COPY_TIMER (aligner_best_match_FLAT_search);    
    else
        COPY_TIMER (aligner_best_match_FLAT_search2);    
}

// returns gpos aligned with seq with M (as in CIGAR) length, containing the longest match to the reference. 
// returns false if no match found.
// note: matches that imply a negative GPOS (i.e. their beginning is aligned to before the start of the genome), aren't consisdered
static inline PosType64 aligner_best_match_FLAT (
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
    aligner_best_match_FLAT_search (vb, &seq_bits, STRa(seq), genome, genome_nbases, gpos_R1, NO_GPOS,
                                    &best);

    int match_percent = percent_match (best.match_n_bits, seq_len); // note: this is matching bits, not matching bases

    *is_perfect = best.is_perfect && seq_is_all_acgt; // if !seq_is_all_acgt (=seq has an 'N'), it is never a perfect match even match_n_bits says so, because 'N's are converted to 'A's in the seq bitmap, and might therefore inflate match_n_bits by matching an 'A' in the genome
    *is_forward = best.is_forward;
    *can_match_as_nonspliced = (match_percent >= MIN_MATCH_PC);                                         

    // case: spliced alignment allowed, and match is relatively short - check for another segment (might be a biological exon, a biological indel variant, or a sequencer indel error)
    if (!flag.no_splicing 
     && IN_RANGX (match_percent, SPLICE_MIN_PC, SPLICE_MAX_PC) 
     && seq_len >= SPLICE_MIN_SEQ_LEN) {
        
        best2.is_forward = best.is_forward; // spliced segments must be in the same orientation
        aligner_best_match_FLAT_search (vb, &seq_bits, STRa(seq), genome, genome_nbases, NO_GPOS, best.gpos, &best2);
        *gpos2 = best2.gpos;
    }

    COPY_TIMER (aligner_best_match);

    // minimum critiera to be worth segging the alignment. note for spliced: this is an initial filtering, verifying that 
    // sum of matching bits (even double-counting overlaps) matches the criteria. aligner_seg_mismatches() will re-test more tightly after junction is known.
    if (percent_match (best.match_n_bits + best2.match_n_bits, seq_len) < MIN_MATCH_PC)
        return NO_GPOS;
    else
        return best.gpos;
}

static inline PosType64 aligner_best_match_LAYERED (VBlockP vb, STRp(seq), PosType64 gpos_R1,
                                                    ConstBitsP genome, PosType64 genome_nbases,
                                                    bool *is_forward, bool *is_perfect) // out
{
    uint64_t kmer;
    const PosType64 seq_len_64 = (PosType64)seq_len; // 64 bit version of seq_len
    PosType64 gpos = NO_GPOS;
    BestAlignment best = { .gpos = NO_GPOS };
    
    // convert seq to a bitmap
    bool seq_is_all_acgt;
    uint64_t seq_bits_words[roundup_bits2words64(seq_len * 2)];
    Bits seq_bits = aligner_seq_to_bitmap (vb, STRa(seq), seq_bits_words, &seq_is_all_acgt);
    
    *is_perfect = false;

    typedef enum { NOT_FOUND=-1, REVERSE=0, FORWARD=1 } Direction;

    // each "find" corresponds to a place in seq that has the hook base (or if a homopolymer of the hook - the last base in the homopolymer)
    struct Finds { 
        uint64_t kmer;   // the sequence (in 2bit) directly before (if fwd) or after (if rev) the hook. Used as a key into refhash. 
        uint32_t i;      // index of hook in seq
        Direction found; // orientation of find
    } finds[seq_len]; 
    
    uint32_t num_finds = 0;
    
    uint32_t density = (flag.fast ? 3 : 1);
    uint32_t max_snps_for_near_perfect = (flag.fast ? 10 : 2); // note: this is only approximately SNPs as we actually measure mismatched bits, not bases

    // we search - checking both forward hooks and reverse hooks, we check only the first layer for now
    for (int i=0; i < seq_len; i += density) {          
        Direction found = NOT_FOUND;

        if (__builtin_expect (i < seq_len - bases_per_hash, true/*probability 93%*/) && // room for the hash word
            __builtin_expect (seq[i] == HOOK,  false/*probability 75%*/) && 
            __builtin_expect (seq[i+1] != HOOK, true/*probability 75%*/) &&  // take the G - if there is a homopolymer GGGG... take the last one
            (gpos = aligner_get_gpos_LAYERED (vb, &seq_bits, i+1, &kmer, FORWARD)) != NO_GPOS) { 

            gpos -= i; // gpos is the first base on the reference, that aligns to the first base of seq
            found = FORWARD;
        }

        // note: if the previous condition is true, then seq[i]==HOOK and no need to check this condition. hence the "else"
        else if (__builtin_expect (i >= bases_per_hash, true/*probability 93%*/) && // room for the hash word
            __builtin_expect (seq[i] == HOOK_REV,  false/*probability 75%*/) && 
            __builtin_expect (seq[i-1] != HOOK_REV, true/*probability 75%*/) &&  // take the G - if there is a polymer GGGG... take the last one
            (gpos = aligner_get_gpos_LAYERED (vb, &seq_bits, i-1, &kmer, REVERSE)) != NO_GPOS) { 

            gpos -= seq_len_64-1 - i; // gpos is the first base of the reference, that aligns wit the LAST base of seq
            found = REVERSE;
        }

        if (found != NOT_FOUND && 
            __builtin_expect (gpos >= 0 && gpos + seq_len_64 < genome_nbases, true)) { // ignore this gpos if the seq wouldn't fall completely within reference genome

            finds[num_finds++] = (struct Finds){ .kmer = kmer, .i = i, .found = found };
            
            if (aligner_update_best (vb, gpos, gpos_R1, NO_GPOS, &seq_bits, seq_len, found, 
                                     genome, genome_nbases, max_snps_for_near_perfect, &best)) {
                goto done; // near-perfect match, search no longer
            }
        }
    }

    // if still no near-perfect matches found, search the additional layers 
    // (each kmer can have up to NUM_LAYERS corresponding GPOSes - we already tested one, now we test the rest)    
    for (uint32_t find_i=0; find_i < num_finds; find_i++) {
        // amount to adjust gpos, which was found for base #i of seq, to correspond to the beginning of seq
        int64_t gpos_shift = (finds[find_i].found == FORWARD ? (int64_t)finds[find_i].i : (seq_len_64-1 - (int64_t)finds[find_i].i));

        for (unsigned layer_i=1; layer_i < NUM_LAYERS; layer_i++) {
            gpos = (PosType64)BGEN32 (*B32(refhash_buf, LAYER_START[layer_i] + (finds[find_i].kmer & layer_bitmask[layer_i]))); 

            // case: no more GPOSes for this kmer - this layer and all subsequent are NO_GPOS
            if (gpos == NO_GPOS) break;

            gpos -= gpos_shift;

            // check if this gpos results in a better match than the one we already have
            if (__builtin_expect ((gpos >= 0) && (gpos + seq_len_64 < genome_nbases), true) &&
                aligner_update_best (vb, gpos, gpos_R1, NO_GPOS, &seq_bits, seq_len, finds[find_i].found, 
                                     genome, genome_nbases, max_snps_for_near_perfect, &best)) {
                goto done; // near-perfect match, search no longer       
            }
        }
    }
    
done:
    *is_forward = best.is_forward;
    *is_perfect = seq_is_all_acgt && best.is_perfect;

    // check that the best match is above a threshold, where mapping against the reference actually improves compression
    if (percent_match (best.match_n_bits, seq_len) >= MIN_MATCH_PC) 
        return best.gpos;
    else
        return NO_GPOS;
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
    // **FORWARD**    ------- SEG -------       ---------- PIZ ----------
    // seq:           0·····j···········L       0·····j    0···········L2   (g=gpos j=junction L=(seq_len-1) L2=L-j)
    // 1st segment ⟶ g···············g+L       G1         G2               SEG: GPOS ⇐ g ; JUNCTION ⇐ j ; GAP ⇐ g2-g1
    //                ‾‾‾‾‾‾_____________                                   PIZ: G1 ⇐ GPOS ; L2 ⇐ seq_len-JUNCTION   
    // 2st segment ⟶ g2·············g2+L                                        G2 ⇐ G1+GAP+JUNCTION
    //
    // **REVERSE**    ------- SEG -------       ---------- PIZ ----------   SEG: GPOS ⇐ g+(seq_len-j)  
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
        ref_get_textual_seq (min_gpos, ref_data, ref_len, !is_forward);

        *ref  = (gpos  == min_gpos) ? ref_data : (ref_data + abs_gpos_gap);
        *ref2 = (gpos2 == min_gpos) ? ref_data : (ref_data + abs_gpos_gap);

        // note: if reverse, entire ref_data is revcomped, so ref1 last base is gpos
        if (!is_forward) SWAP (*ref, *ref2);
    } 

    // two disjoint sequences, each
    else {
        char *ref_data = alloc_ref (vb, STRa(ref_space), seq_len * 2);

        // note: if reverse, ref and ref2 are revcomped in-place independently
        ref_get_textual_seq (gpos, ref_data, seq_len, !is_forward);
        *ref = ref_data;

        ref_get_textual_seq (gpos2, &ref_data[seq_len], seq_len, !is_forward);
        *ref2 = ref_data + seq_len;
    }
}

// segs SQBITMAP(local), SEQMIS_*, JUNCTION. Returns true is is_perfect (no mismatches)
static bool aligner_seg_mismatches (VBlock𐤐 vb, STR𐤐(seq), 
                                    PosType64𐤐 gpos, PosType64𐤐 gpos2, // might be swapped
                                    bool is_forward, bool can_match_as_nonspliced,
                                    uint32_t *restrict junction) // out
{
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
            if (flag.debug_aligner) iprintf ("%c", is_forward ? *r : COMPLEM[(int)*r]); 
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
      : refhash_is_flat                ? aligner_best_match_FLAT    (VB, STRa(seq), gpos_R1, genome, genome_nbases, &is_forward, &is_perfect, &gpos2, &can_match_as_nonspliced) 
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
    // **FORWARD**    ------- SEG -------       ---------- PIZ ----------
    // seq:           0·····j···········L       0·····j    0···········L2   (g=gpos j=junction L=(seq_len-1) L2=L-j)
    // 1st segment ⟶ g···············g+L       G1         G2               SEG: GPOS ⇐ g ; JUNCTION ⇐ j ; GAP ⇐ g2-g1
    //                ‾‾‾‾‾‾_____________                                   PIZ: G1 ⇐ GPOS ; L2 ⇐ seq_len-JUNCTION   
    // 2st segment ⟶ g2·············g2+L                                        G2 ⇐ G1+GAP+JUNCTION
    //
    // **REVERSE**    ------- SEG -------       ---------- PIZ ----------   SEG: GPOS ⇐ g+(seq_len-j)  
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
        ref_get_textual_seq (gpos, BAFTtxt, seq_len, !is_forward);

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
