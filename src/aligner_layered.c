// ------------------------------------------------------------------
//   aligner_layers.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// backward compatability with reference files generated with genozip up to 15.0.80 (layered refhash)

#include "seg.h"
#include "refhash_friend.h"
#include "file.h"
#include "codec.h"
#include "reconstruct.h"
#include "piz.h"
#include "aligner.h"

// functions in aligner.c
static inline uint64_t aligner_get_kmer (Bits𐤐 seq, uint64_t bit_i, bool is_forward);
static Bits aligner_seq_to_bitmap (VBlockP vb, rom𐤐 seq, uint64_t seq_len, uint64_t *restrict bitmap_words, bool𐤐 seq_is_all_acgt);
static inline int percent_match (uint32_t n_matching, uint32_t seq_len);
static inline bool aligner_update_best (VBlockP vb, PosType64 gpos, PosType64 gpos_R1, PosType64 gpos1, Bits𐤐 seq_bits, uint32_t seq_len, bool fwd, ConstBits𐤐 genome, PosType64 genome_nbases, uint32_t near_perfect_max_mismatches, BestAlignment *restrict best);

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
    
    uint32_t near_perfect_max_mismatches = (flag.fast?NEAR_PERFECT_FAST : flag.best?NEAR_PERFECT_BEST : NEAR_PERFECT_NORM);

    // we search - checking both forward hooks and reverse hooks, we check only the first layer for now
    for (int i=0; i < seq_len; i++) {          
        Direction found = NOT_FOUND;

        if (i < seq_len - bases_per_hash && // room for the hash word
            seq[i] == HOOK && 
            seq[i+1] != HOOK &&  // take the G - if there is a homopolymer GGGG... take the last one
            (gpos = aligner_get_gpos_LAYERED (vb, &seq_bits, i+1, &kmer, FORWARD)) != NO_GPOS) { 

            gpos -= i; // gpos is the first base on the reference, that aligns to the first base of seq
            found = FORWARD;
        }

        // note: if the previous condition is true, then seq[i]==HOOK and no need to check this condition. hence the "else"
        else if (i >= bases_per_hash && // room for the hash word
            seq[i] == HOOK_REV && 
            seq[i-1] != HOOK_REV &&  // take the G - if there is a polymer GGGG... take the last one
            (gpos = aligner_get_gpos_LAYERED (vb, &seq_bits, i-1, &kmer, REVERSE)) != NO_GPOS) { 

            gpos -= seq_len_64-1 - i; // gpos is the first base of the reference, that aligns wit the LAST base of seq
            found = REVERSE;
        }

        if (found != NOT_FOUND && 
            __builtin_expect (gpos >= 0 && gpos + seq_len_64 < genome_nbases, true)) { // ignore this gpos if the seq wouldn't fall completely within reference genome

            finds[num_finds++] = (struct Finds){ .kmer = kmer, .i = i, .found = found };
            
            if (gpos != best.gpos &&
                aligner_update_best (vb, gpos, gpos_R1, NO_GPOS, &seq_bits, seq_len, found, 
                                     genome, genome_nbases, near_perfect_max_mismatches, &best)) {
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
                gpos != best.gpos &&
                aligner_update_best (vb, gpos, gpos_R1, NO_GPOS, &seq_bits, seq_len, finds[find_i].found, 
                                     genome, genome_nbases, near_perfect_max_mismatches, &best)) {
                goto done; // near-perfect match, search no longer       
            }
        }
    }
    
done:
    *is_forward = best.is_forward;
    *is_perfect = seq_is_all_acgt && best.is_perfect;

    // check that the best match is above a threshold, where mapping against the reference actually improves compression
    if (percent_match (best.n_matching, seq_len) >= MIN_MATCH_PC) 
        return best.gpos;
    else
        return NO_GPOS;
}
