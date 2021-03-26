// ------------------------------------------------------------------
//   aligner.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

#define MAX_GPOS_DELTA 1000 // paired reads are usually with a delta less than 300 - so this is more than enough

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
static inline bool aligner_get_word_from_seq (VBlock *vb, const char *seq, uint32_t *refhash_word, int direction /* 1 forward, -1 reverse */)
                                              
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
BitArray aligner_seq_to_bitmap (const char *seq, word_t seq_len, 
                                word_t *bitmap_words,  // allocated by caller
                                bool *seq_is_all_actg) // optional out
{
    // covert seq to 2-bit array
    BitArray seq_bits = { .nbits  = seq_len * 2, 
                          .nwords = roundup_bits2words64(seq_len * 2), 
                          .words        = bitmap_words,
                          .type         = BITARR_REGULAR };

    if (seq_is_all_actg) *seq_is_all_actg = true; // starting optimistic

    for (word_t base_i=0; base_i < seq_len; base_i++) {
        uint8_t encoding = nuke_encode[(uint8_t)seq[base_i]];
    
        if (encoding == 4) { // not A, C, G or T - usually N
            if (seq_is_all_actg) *seq_is_all_actg = false;
            encoding = 0; // arbitrarily convert 4 (any non-actg is 4) to 0 ('A')
        }
    
        bit_array_assign2 (&seq_bits, (base_i << 1), encoding);
    }

    bit_array_clear_excess_bits_in_top_word (&seq_bits);

    return seq_bits;
}

// returns gpos aligned with seq with M (as in CIGAR) length, containing the longest match to the reference. 
// returns false if no match found.
// note: matches that imply a negative GPOS (i.e. their beginning is aligned to before the start of the genome), aren't consisdered
static inline PosType aligner_best_match (VBlock *vb, const char *seq, const uint32_t seq_len,
                                          bool *is_forward, bool *is_all_ref) // out
{
    START_TIMER;

    uint32_t longest_len=0; // longest number of bits (not bases!) that match
    uint32_t refhash_word;
    const PosType seq_len_64 = (PosType)seq_len; // 64 bit version of seq_len
    
    PosType gpos = NO_GPOS, best_gpos = NO_GPOS; // match not found yet
    bool best_is_forward = false;

    // convert seq to a bitmap
    bool maybe_perfect_match;
    word_t seq_bits_words[roundup_bits2words64(seq_len * 2)];
    BitArray seq_bits = aligner_seq_to_bitmap (seq, seq_len, seq_bits_words, &maybe_perfect_match);
    
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

#       define UPDATE_BEST(fwd)  {               \
            if (gpos != best_gpos) {             \
                uint32_t match_len = (uint32_t)seq_bits.nbits - \
                    bit_array_manhattan_distance ((fwd) ? genome : emoneg, \
                                                  ((fwd) ? gpos : genome_nbases-1 - (gpos + seq_bits.nbits/2 -1)) * 2, \
                                                  &seq_bits, 0, \
                                                  seq_bits.nbits); \
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
        }

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

            if ((gpos >= 0) && (gpos + seq_len_64 < genome_nbases)) {
                UPDATE_BEST (finds[find_i].found);
            }
        }
    }

done:
    *is_forward = best_is_forward;

    COPY_TIMER (aligner_best_match);
    return best_gpos;
}

void aligner_seg_seq (VBlockP vb, ContextP bitmap_ctx, const char *seq, uint32_t seq_len)
{
    // these 4 contexts are consecutive and in the same order for all relevant data_types in data_types.h
    Context *nonref_ctx = bitmap_ctx + 1; // NONREF
    Context *gpos_ctx   = bitmap_ctx + 3; // GPOS
    Context *strand_ctx = bitmap_ctx + 4; // STRAND

    BitArray *bitmap = buf_get_bitarray (&bitmap_ctx->local);

    // allocate bitmaps - provide name only if buffer is not allocated, to avoid re-writing param which would overwrite nbits that overlays it + param must be 0
    buf_alloc_old (vb, &bitmap_ctx->local, MAX (bitmap_ctx->local.len + roundup_bits2bytes64 (seq_len), vb->lines.len * (seq_len+5) / 8), CTX_GROWTH, 
               buf_is_allocated (&bitmap_ctx->local) ? NULL : "contexts->local"); 

    buf_alloc_old (vb, &strand_ctx->local, MAX (nonref_ctx->local.len + sizeof (int64_t), roundup_bits2bytes64 (vb->lines.len)), CTX_GROWTH, 
               buf_is_allocated (&strand_ctx->local) ? NULL : "contexts->local"); 

    buf_alloc_old (vb, &nonref_ctx->local, MAX (nonref_ctx->local.len + seq_len + 3, vb->lines.len * seq_len / 4), CTX_GROWTH, "contexts->local"); 
    buf_alloc_old (vb, &gpos_ctx->local,   MAX (nonref_ctx->local.len + sizeof (uint32_t), vb->lines.len * sizeof (uint32_t)), CTX_GROWTH, "contexts->local"); 

    bool is_forward, is_all_ref;
    PosType gpos = aligner_best_match ((VBlockP)vb, seq, seq_len, &is_forward, &is_all_ref);

    // case: we're the 2nd of the pair - the bit represents whether this strand is equal to the pair's strand (expecting
    // it to be 1 in most cases - making the bitmap highly compressible)
    if (gpos_ctx->pair_local) {
        const BitArray *pair_strand = buf_get_bitarray (&strand_ctx->pair);
        
        ASSERTE (vb->line_i < pair_strand->nbits, "vb=%u cannot get pair-1 STRAND bit for line_i=%u because pair-1 strand bitarray has only %u bits",
                 vb->vblock_i, vb->line_i, (unsigned)pair_strand->nbits);

        bool pair_is_forward = bit_array_get (pair_strand, vb->line_i); // same location, in the pair's local
        buf_add_bit (&strand_ctx->local, is_forward == pair_is_forward);
    }
    // case: not 2nd in a pair - just store the strange
    else 
        buf_add_bit (&strand_ctx->local, is_forward);
    
    buf_extend_bits (&bitmap_ctx->local, seq_len);

    ASSSEG ((gpos >= 0 && gpos <= MAX_GPOS) || gpos == NO_GPOS, seq, "gpos=%"PRId64" is out of range [0,%"PRId64"]", gpos, MAX_GPOS);
    
    // case: we're the 2nd of the pair - store a delta if its small enough, or a lookup from local if not
    bool store_local = true;
    if (gpos_ctx->pair_local) {

        ASSERTE (vb->line_i < gpos_ctx->pair.len, "vb=%u cannot get pair-1 GPOS for line_i=%u because pair-1 GPOS.len=%"PRIu64,
                 vb->vblock_i, vb->line_i, gpos_ctx->pair.len);

        PosType pair_gpos = (PosType)*ENT (uint32_t, gpos_ctx->pair, vb->line_i); // same location, in the pair's local
        PosType gpos_delta = gpos - pair_gpos; 

        if (gpos != NO_GPOS && gpos_delta <= MAX_GPOS_DELTA && gpos_delta >= -MAX_GPOS_DELTA) {
            store_local = false;      

            char delta_snip[30] = { SNIP_PAIR_DELTA };
            unsigned delta_str_len = str_int (gpos_delta, &delta_snip[1]);
            seg_by_ctx (vb, delta_snip, delta_str_len + 1, gpos_ctx, 0);
        }
        else {
            static const char lookup[1] = { SNIP_LOOKUP }; // lookup from local
            seg_by_ctx (vb, lookup, 1, gpos_ctx, 0);
        }
    }
    
    // store the GPOS in local if its not a 2nd pair, or if it is, but the delta is not small enough
    if (store_local)
        NEXTENT (uint32_t, gpos_ctx->local) = BGEN32 ((uint32_t)gpos);

    // shortcut if there's no reference match
    if (gpos == NO_GPOS) {
        bit_array_clear_region (bitmap, bitmap_ctx->next_local, seq_len); // no bases match the reference
        bitmap_ctx->next_local += seq_len;
        buf_add (&nonref_ctx->local, seq, seq_len);
        return;
    }

    // lock region of reference to protect is_set
    RefLock lock = (flag.reference == REF_EXT_STORE) ? ref_lock (gpos, seq_len) : REFLOCK_NONE;

    // shortcut if we have a full reference match
    if (is_all_ref) {
        bit_array_set_region (bitmap, bitmap_ctx->next_local, seq_len); // all bases match the reference
        bitmap_ctx->next_local += seq_len;
        
        if (flag.reference == REF_EXT_STORE) 
            bit_array_set_region (genome_is_set, gpos, seq_len); // this region of the reference is used (in case we want to store it with REF_EXT_STORE)

        goto done;
    }

    PosType room_fwd = genome_nbases - gpos; // how much reference forward might contain a match

    bit_index_t next_bit = bitmap_ctx->next_local; // copy to automatic variable (optimized to a register) for performace
    for (uint32_t i=0; i < seq_len; i++) {
                
        bool use_reference = false;

        // case our seq is identical to the reference at this site
        if (i < room_fwd) {
            char seq_base = is_forward ? seq[i] : complement[(uint8_t)seq[i]];
            
            PosType ref_i = gpos + (is_forward ? i : seq_len-1-i);
            char ref_base = ACGT_DECODE (genome, ref_i);

            if (seq_base == ref_base) {
                
                // TO DO: replace this with bit_array_or_with (dst, start, len, src, start) (dst=is_set, src=bitmap) (bug 174)
                if (flag.reference == REF_EXT_STORE) 
                    bit_array_set (genome_is_set, ref_i); // we will need this ref to reconstruct

                use_reference = true;
            }
        }

        // case: we can't use the reference (different value than base or we have passed the end of the reference)
        if (!use_reference) 
            NEXTENT (char, nonref_ctx->local) = seq[i];

        bit_array_assign (bitmap, next_bit, use_reference);
        next_bit++; // can't increment inside macro
    }
    bitmap_ctx->next_local = next_bit;

done:
    ref_unlock (lock);
}

// PIZ: SEQ reconstruction 
void aligner_reconstruct_seq (VBlockP vb, ContextP bitmap_ctx, uint32_t seq_len, bool is_pair_2) // true if this is a second file of a pair
{
    if (piz_is_skip_section (vb, SEC_LOCAL, bitmap_ctx->dict_id)) return; // if case we need to skip the SEQ field (for the entire file)

    Context *nonref_ctx = bitmap_ctx + 1;
    Context *gpos_ctx   = bitmap_ctx + 3;
    Context *strand_ctx = bitmap_ctx + 4;
    
    if (buf_is_allocated (&bitmap_ctx->local)) { // not all non-ref

        bool is_forward;
        PosType gpos;

        // first file of a pair ("pair 1") or a non-pair fastq or sam
        if (!is_pair_2) {
            gpos = NEXTLOCAL (uint32_t, gpos_ctx);
            is_forward = NEXTLOCALBIT (strand_ctx);        
        }

        // 2nd file of a pair ("pair 2")
        else {
            // the strand bit is a 1 iff the strand is the same as the pair
            bool is_forward_pair_1 = PAIRBIT (strand_ctx);
            is_forward = NEXTLOCALBIT (strand_ctx) ? is_forward_pair_1 : !is_forward_pair_1;

            // gpos: reconstruct, then cancel the reconstruction and just use last_value
            int32_t reconstructed_len = reconstruct_from_ctx (vb, gpos_ctx->did_i, 0, true);
            vb->txt_data.len -= reconstructed_len; // roll back reconstruction
            gpos = gpos_ctx->last_value.i;
        }

        // sanity check - the sequence is supposed to fit in the 
        ASSERTE (gpos == NO_GPOS || gpos + seq_len <= genome->nbits / 2, "gpos=%"PRId64" is out of range: seq_len=%u and genome_nbases=%"PRIu64,
                 gpos, seq_len, genome->nbits / 2);

        if (is_forward)  // normal (note: this condition test is outside of the tight loop)
            for (uint32_t i=0; i < seq_len; i++)
                if (NEXTLOCALBIT (bitmap_ctx))  // get base from reference
                    RECONSTRUCT1 (ACGT_DECODE (genome, gpos + i));
                else  // get base from nonref
                    RECONSTRUCT1 (NEXTLOCAL (char, nonref_ctx));

        else // reverse complement
            for (uint32_t i=0; i < seq_len; i++) 
                if (NEXTLOCALBIT (bitmap_ctx))  // case: get base from reference
                    RECONSTRUCT1 (complement [(uint8_t)ACGT_DECODE (genome, gpos + seq_len-1 - i)]);
                else  // case: get base from nonref
                    RECONSTRUCT1 (NEXTLOCAL (char, nonref_ctx));
    }
    else {
        RECONSTRUCT (ENT (char, nonref_ctx->local, nonref_ctx->next_local), seq_len);
        nonref_ctx->next_local += seq_len;
    }

    ASSERTE (nonref_ctx->next_local <= nonref_ctx->local.len, "nonref_ctx->next_local=%u is out of range: nonref_ctx->local.len=%u",
             nonref_ctx->next_local, (uint32_t)nonref_ctx->local.len);
}
