// ------------------------------------------------------------------
//   aligner.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "context.h"
#include "buffer.h"
#include "vblock.h"
#include "seg.h"
#include "refhash.h"
#include "endianness.h"
#include "file.h"
#include "strings.h"
#include "codec.h"
#include "reconstruct.h"
#include "piz.h"
#include "segconf.h"
#include "fastq.h"
#include "aligner.h"
#include "profiler.h"

#define COMPLIMENT(b) (3-(b))

// strict encoding of A,C,G,T - everything else in non-encodable (a 4 here)
static const uint32_t nuke_encode[256] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 4
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 16
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 32
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 48
                                           4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,   // 64  A(65)->0 C(67)->1 G(71)->2
                                           4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 84  T(84)->3
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 96  
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 112 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 128
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

// Foward example: If seq is: G-AGGGCT  (G is the hook)  -- matches reference AGGGCT       - function returns 110110101000 (A=00 is the LSb)
// Reverse       : If seq is: CGCCCT-C  (C is the hook)  -- also matches reference AGGGCT  - function returns 110110101000 - the same
// calculates a refhash word from 14 nucleotides following a 'G' (only last G in a sequenece of GGGG...)
static inline bool aligner_get_word_from_seq (VBlockP vb, rom seq, uint32_t *refhash_word, int direction /* 1 forward, -1 reverse */)
                                              
{   
    START_TIMER;

    *refhash_word = 0;

    for (int i=0; direction * i < nukes_per_hash; i += direction) {   
        uint32_t base = nuke_encode[(uint8_t)seq[i]];
        if (base == 4) {
            COPY_TIMER (aligner_get_word_from_seq);
            return false; // not a A,C,G,T
        }

        if (direction == -1) base = COMPLIMENT (base);

        *refhash_word |= (base << (direction * i * 2)); // 2-LSb of word is the first base
    }

    // remove MSb in case of an odd number of base layer bits
    if (bits_per_hash_is_odd) 
        *refhash_word &= layer_bitmask[0]; 

    COPY_TIMER (aligner_get_word_from_seq);
    return true;
}

// converts a string sequence to a 2-bit bitmap
Bits aligner_seq_to_bitmap (rom seq, uint64_t seq_len, 
                                uint64_t *bitmap_words,  // allocated by caller
                                bool *seq_is_all_actg) // optional out
{
    // covert seq to 2-bit array
    Bits seq_bits = { .nbits  = seq_len * 2, 
                      .nwords = roundup_bits2words64(seq_len * 2), 
                      .words  = bitmap_words,
                      .type   = BITARR_REGULAR };

    if (seq_is_all_actg) *seq_is_all_actg = true; // starting optimistic

    for (uint64_t base_i=0; base_i < seq_len; base_i++) {
        uint8_t encoding = nuke_encode[(uint8_t)seq[base_i]];
    
        if (encoding == 4) { // not A, C, G or T - usually N
            if (seq_is_all_actg) *seq_is_all_actg = false;
            encoding = 0; // arbitrarily convert 4 (any non-actg is 4) to 0 ('A')
        }
    
        bits_assign2 (&seq_bits, (base_i << 1), encoding);
    }

    bits_clear_excess_bits_in_top_word (&seq_bits);

    return seq_bits;
}

// returns gpos aligned with seq with M (as in CIGAR) length, containing the longest match to the reference. 
// returns false if no match found.
// note: matches that imply a negative GPOS (i.e. their beginning is aligned to before the start of the genome), aren't consisdered
static inline PosType aligner_best_match (VBlockP vb, STRp(seq), PosType pair_gpos,
                                          ConstBitsP genome, ConstBitsP emoneg, PosType genome_nbases,
                                          bool *is_forward, bool *is_all_ref) // out
{
    START_TIMER;

    int32_t/*signed*/ longest_len=0; // longest number of bits (not bases!) that match
    uint32_t refhash_word;
    const PosType seq_len_64 = (PosType)seq_len; // 64 bit version of seq_len
    
    PosType gpos = NO_GPOS, best_gpos = NO_GPOS; // match not found yet
    bool best_is_forward = false;
    
    // convert seq to a bitmap
    bool maybe_perfect_match;
    uint64_t seq_bits_words[roundup_bits2words64(seq_len * 2)];
    Bits seq_bits = aligner_seq_to_bitmap (seq, seq_len, seq_bits_words, &maybe_perfect_match);
    
    //ref_print_bases (&seq_bits, "\nseq_bits fwd", true);
    //ref_print_bases (&seq_bits, "seq_bits rev", false);

    *is_all_ref = false;

    typedef enum { NOT_FOUND=-1, REVERSE=0, FORWARD=1 } Direction;

    struct Finds { uint32_t refhash_word; uint32_t i; Direction found; } finds[seq_len];
    uint32_t num_finds = 0;
    
    // in case of --fast, we check only 1/5 of the bases, and we are content with a match (not searching any further) if it 
    // has at most 10 SNPs. On our test file, this reduced the number of calls to aligner_get_match_len by about 4X, 
    // at the cost of the compressed file being about 11% larger
    uint32_t density = (flag.fast ? 5 : 1);
    uint32_t max_snps_for_perfection = (flag.fast ? 10 : 2);

    // we search - checking both forward hooks and reverse hooks, we check only the first layer for now
    for (uint32_t i=0; i < seq_len; i += density) {    
        
        Direction found = NOT_FOUND;

        if (i < seq_len - nukes_per_hash && // room for the hash word
            seq[i] == HOOK && seq[i+1] != HOOK &&  // take the G - if there is a polymer GGGG... take the last one
            aligner_get_word_from_seq (vb, &seq[i+1], &refhash_word, 1)) { 

            gpos = (PosType)BGEN32 (refhashs[0][refhash_word & layer_bitmask[0]]); // position of the start of the G... sequence in the genome
            
            if (gpos != NO_GPOS) {
                gpos -= i; // gpos is the first base on the reference, that aligns to the first base of seq
                found = FORWARD;
            }
        }

        else if (i >= nukes_per_hash && // room for the hash word
            seq[i] == HOOK_REV && seq[i-1] != HOOK_REV &&  // take the G - if there is a polymer GGGG... take the last one
            aligner_get_word_from_seq (vb, &seq[i-1], &refhash_word, -1)) { 

            gpos = (PosType)BGEN32 (refhashs[0][refhash_word & layer_bitmask[0]]); // position of the start of the G... sequence in the FORWARD genome
            
            if (gpos != NO_GPOS) {
                gpos -= seq_len_64-1 - i; // gpos is the first base of the reference, that aligns wit the LAST base of seq
                found = REVERSE;
            }
        }

#       define UPDATE_BEST(fwd)  ({              \
            if (gpos != best_gpos) {             \
                int32_t match_len = (uint32_t)seq_bits.nbits - \
                    bits_hamming_distance ((fwd) ? genome : emoneg, \
                                               ((fwd) ? gpos : genome_nbases-1 - (gpos + seq_bits.nbits/2 -1)) * 2, \
                                               &seq_bits, 0, \
                                               seq_bits.nbits); \
                if (pair_gpos != NO_GPOS && ABS(gpos-pair_gpos) > 500) match_len -= 17; /* penalty for remote GPOS in 2nd pair */ \
                if (match_len > longest_len) {   \
                    longest_len     = match_len; \
                    best_gpos       = gpos;      \
                    best_is_forward = (fwd);     \
                    /* note: we allow 2 snps and we still consider the match good enough and stop looking further */\
                    /* compared to stopping only if match_len==seq_len, this adds about 1% to the file size, but is significantly faster */\
                    if (match_len >= (seq_len - max_snps_for_perfection) * 2) { /* we found (almost) the best possible match */ \
                        *is_all_ref = maybe_perfect_match && (match_len == seq_len*2); /* perfect match */ \
                        goto done;               \
                    }                            \
                }                                \
            }                                    \
        })

        if (found != NOT_FOUND && (gpos >= 0) && (gpos != NO_GPOS) && (gpos + seq_len_64 < genome_nbases)) { // ignore this gpos if the seq wouldn't fall completely within reference genome
            finds[num_finds++] = (struct Finds){ .refhash_word = refhash_word, .i = i, .found = found };
            UPDATE_BEST (found);
        }
    }

    // if still no near-perfect matches found, search the additional layers
    for (unsigned layer_i=1; layer_i < num_layers; layer_i++) {
        for (uint32_t find_i=0; find_i < num_finds; find_i++) {

            if (finds[find_i].found == NOT_FOUND) continue;

            gpos = (PosType)BGEN32 (refhashs[layer_i][finds[find_i].refhash_word & layer_bitmask[layer_i]]); 

            if (gpos == NO_GPOS) {
                finds[find_i].found = NOT_FOUND; // if we can't find it in this layer, we won't find it in the next layers either
                continue;
            }

            gpos -= (finds[find_i].found == FORWARD ? finds[find_i].i : seq_len_64-1 - finds[find_i].i);

            if ((gpos >= 0) && (gpos + seq_len_64 < genome_nbases)) 
                UPDATE_BEST (finds[find_i].found);
        }
    }

done:
    *is_forward = best_is_forward;

    // check that the best match is above a threshold, where mapping against the reference actually improves compression
    if (longest_len * 50 / seq_len < 73) best_gpos = NO_GPOS; // note that longest_len==seq_len*2 is perfect match (longest_len is in bits)

    COPY_TIMER (aligner_best_match);
    return best_gpos;
}

MappingType aligner_seg_seq (VBlockP vb, ContextP bitmap_ctx, STRp(seq), bool no_bitmap_if_perfect,
                             bool is_pair_2, PosType pair_gpos, bool pair_is_forward)
{
    START_TIMER;

    ConstBitsP genome, emoneg;
    PosType genome_nbases;
    ref_get_genome (gref, &genome, &emoneg, &genome_nbases);

    // these 4 contexts are consecutive and in the same order for all relevant data_types in data_types.h
    ContextP gpos_ctx   = bitmap_ctx + 3; // GPOS
    ContextP strand_ctx = bitmap_ctx + 4; // STRAND
    ContextP seqmis_ctx = bitmap_ctx + 5; // SEQMIS_A/C/G/T (4 contexts)

    bool is_forward=false, is_all_ref=false;

    // our aligner algorithm only works for short reads - long reads tend to have many Indel differences (mostly errors) vs the reference
    #define MAX_SHORT_READ_LEN 2500
    PosType gpos = (seq_len <= MAX_SHORT_READ_LEN && !segconf.running) 
                 ? aligner_best_match (VB, STRa(seq), pair_gpos, genome, emoneg, genome_nbases, &is_forward, &is_all_ref) : NO_GPOS; 

    if (gpos == NO_GPOS || gpos > genome_nbases - seq_len || gpos > MAX_ALIGNER_GPOS - seq_len || gpos < 0/*never happens*/) {
        COPY_TIMER (aligner_seg_seq);
        return MAPPING_NO_MAPPING;
    }

    if (flag.show_aligner)
        iprintf ("%s: gpos=%"PRId64" forward=%s\n", LN_NAME, gpos, TF(is_forward));

    Bits *bitmap = (BitsP)&bitmap_ctx->local;

    // allocate bitmaps - don't provide name to avoid re-writing param which would overwrite nbits that overlays it 
    if (!no_bitmap_if_perfect || !is_all_ref) {
        int64_t missing_bits = (int64_t)seq_len - ((int64_t)bitmap_ctx->local.nbits - (int64_t)bitmap_ctx->next_local);
        if (missing_bits > 0) {
            buf_alloc_do (vb, &bitmap_ctx->local,roundup_bits2bytes64 (bitmap_ctx->local.nbits + seq_len), CTX_GROWTH, __FUNCLINE, NULL); 
            buf_extend_bits (&bitmap_ctx->local, missing_bits);
        }
    }

    if ((strand_ctx->local.size & ~3ULL) * 8 == strand_ctx->local.nbits)
        buf_alloc_do (vb, &strand_ctx->local, strand_ctx->local.size + sizeof(uint64_t), CTX_GROWTH, __FUNCLINE, NULL);
        
    buf_alloc (vb, &gpos_ctx->local, 1, 0, uint32_t, CTX_GROWTH, NULL); 

    buf_add_bit (&strand_ctx->local, pair_gpos == NO_GPOS ? is_forward // pair 1 is unaligned - just store the strand
                                   :                        (is_forward == pair_is_forward)); // pair 1 is aligned aligned - store equality, expected to he 1 in most cases
            
    if (is_pair_2) 
        fastq_seg_pair2_gpos (vb, pair_gpos, gpos);
    
    // store the GPOS in local if its not a 2nd pair, or if it is, but the delta is not small enough
    else
        BNXT32 (gpos_ctx->local) = (uint32_t)gpos;

    // lock region of reference to protect is_set
    RefLock lock = (IS_REF_EXT_STORE) ? ref_lock (gref, gpos, seq_len) : REFLOCK_NONE;

    if (IS_REF_EXT_STORE) 
        ref_set_genome_is_used (gref, gpos, seq_len); // this region of the reference is used (in case we want to store it with REF_EXT_STORE)

    // shortcut if we have a full reference match
    if (is_all_ref) {
        if (!no_bitmap_if_perfect) {
            bits_set_region (bitmap, bitmap_ctx->next_local, seq_len); // all bases match the reference
            bitmap_ctx->next_local += seq_len;
        }

        goto done;
    }

    for (int i=0; i < 4; i++)
        buf_alloc (vb, &seqmis_ctx[i].local, seq_len, 0, char, CTX_GROWTH, NULL); 

    uint64_t next_bit = bitmap_ctx->next_local; // copy to automatic variable (optimized to a register) for performace
    for (uint32_t i=0; i < seq_len; i++) {
                
        // case our seq is identical to the reference at this site
        char seq_base = is_forward ? seq[i] : complement[(uint8_t)seq[i]];
        
        PosType ref_i = gpos + (is_forward ? i : seq_len-1-i);
        
        uint8_t ref_base_2bit = bits_get2 (genome, ref_i * 2);
        char ref_base = acgt_decode(ref_base_2bit);

        if (seq_base == ref_base) 
            bits_set (bitmap, next_bit++); 

        // case: we can't use the reference (different value than base or we have passed the end of the reference)
        else {
            bits_clear (bitmap, next_bit++); 
            BNXTc (seqmis_ctx[ref_base_2bit].local) = seq[i];
        }
    }
    bitmap_ctx->next_local = next_bit;

done:
    ref_unlock (gref, &lock);

    COPY_TIMER (aligner_seg_seq);
    return is_all_ref ? MAPPING_PERFECT : MAPPING_ALIGNED; // successful
}

// PIZ: SEQ reconstruction - only for reads compressed with the aligner
void aligner_reconstruct_seq (VBlockP vb, ContextP bitmap_ctx, uint32_t seq_len, bool is_pair_2, bool is_perfect_alignment, bool reconstruct)
{
    START_TIMER;

    if (!bitmap_ctx->is_loaded) return; // if case we need to skip the SEQ field (for the entire file)

    ContextP nonref_ctx = bitmap_ctx + 1;
    ContextP gpos_ctx   = bitmap_ctx + 3;
    ContextP strand_ctx = bitmap_ctx + 4;
    ContextP seqmis_ctx = bitmap_ctx + 5; // SEQMIS_A/C/G/T (4 contexts)
    
    if (is_perfect_alignment || buf_is_alloc (&bitmap_ctx->local)) { // not all non-ref

        bool is_forward;
        PosType gpos;

        // first file of a pair ("pair 1") or a non-pair fastq or sam
        if (!is_pair_2) {
            gpos = NEXTLOCAL (uint32_t, gpos_ctx);
            is_forward = NEXTLOCALBIT (strand_ctx);        
        }

        // 2nd file of a pair ("pair 2")
        else {
            is_forward = fastq_piz_get_pair2_is_forward (vb); // MUST be called before gpos reconstruction as it inquires GPOS.pair.next

            // gpos: don't reconstruct just get last_value
            reconstruct_from_ctx (vb, gpos_ctx->did_i, 0, false); // calls fastq_special_PAIR2_GPOS
            gpos = gpos_ctx->last_value.i;
        }

        if (flag.show_aligner)
            iprintf ("%s: gpos=%"PRId64" forward=%s\n", LN_NAME, gpos, TF(is_forward));

        const Bits *genome=NULL;
        PosType genome_nbases;
        
        if (gpos != NO_GPOS) { 
            ref_get_genome (gref, &genome, NULL, &genome_nbases);

            // sanity check - the sequence is supposed to fit in the 
            ASSPIZ (gpos + seq_len <= genome->nbits / 2, "gpos=%"PRId64" is out of range: seq_len=%u and genome_nbases=%"PRIu64,
                    gpos, seq_len, genome->nbits / 2);
        }
        else {
            bitmap_ctx->next_local += seq_len; // up to v13 ; not sure if this still can happen since v14
            goto all_nonref;
        }

        BASE_ITER_INIT (genome, gpos, seq_len, is_forward); // careful: segfault if no genome

        // faster loop in the common case of a perfect (= no mismatches) alignment (is_perfect_alignment introduced in v14)
        if (is_perfect_alignment) {
            if (reconstruct) {
                char *next = BAFTtxt;

                for (uint32_t i=0; i < seq_len; i++) 
                    *next++ = acgt_decode(BASE_NEXT);                    
                
                vb->txt_data.len32 += seq_len;
            }
        }

        else {
            for (uint32_t i=0; i < seq_len; i++) {
                uint8_t base = BASE_NEXT; // advance iterator whether v14 or not
                if (NEXTLOCALBIT (bitmap_ctx)) { // get base from reference
                    if (reconstruct) RECONSTRUCT1 (acgt_decode(base));
                }
                else  { // get base from seqmis_ctx
                    ContextP ctx = VER(14) ? (seqmis_ctx + (is_forward ? base : (3 - base))) // if is_forward, ref is revcomped - we complement ref[i] back to its original
                                           : nonref_ctx;
                    RECONSTRUCT_NEXT (ctx, 1); 
                }
            }
        }
    }

    else all_nonref: { 
        ASSPIZ (nonref_ctx->next_local + seq_len <= nonref_ctx->local.len32, "NONREF exhausted: next_local=%u + seq_len=%u > local.len=%u",
                nonref_ctx->next_local, seq_len, nonref_ctx->local.len32);

        if (reconstruct) RECONSTRUCT (Bc (nonref_ctx->local, nonref_ctx->next_local), seq_len);
        nonref_ctx->next_local += seq_len;
    }

    ASSPIZ (nonref_ctx->next_local <= nonref_ctx->local.len, "nonref_ctx->next_local=%u is out of range: nonref_ctx->local.len=%u",
            nonref_ctx->next_local, nonref_ctx->local.len32);

    COPY_TIMER (aligner_reconstruct_seq);
}
