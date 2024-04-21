// ------------------------------------------------------------------
//   aligner.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
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
#include "contigs.h"

// Foward example: If seq is: G-AGGGCT  (G is the hook)  -- matches reference AGGGCT       - function returns 110110101000 (A=00 is the LSb)
// Reverse       : If seq is: CGCCCT-C  (C is the hook)  -- also matches reference AGGGCT  - function returns 110110101000 - the same
// calculates a refhash word from 14 nucleotides following a 'G' (only last G in a sequenece of GGGG...)
static inline PosType64 aligner_get_word_from_seq (VBlockP vb, rom seq, uint32_t *refhash_word, int direction /* 1 forward, -1 reverse */)                                              
{   
    // START_TIMER; // this has a small performance impact as it is called in a tight loop - uncomment when needed

    *refhash_word = 0;

    for (int i=0; direction * i < nukes_per_hash; i += direction) {   
        uint32_t base = (direction == -1) ? nuke_encode_comp (seq[i]) : nuke_encode (seq[i]);
        if (__builtin_expect (base == 4, false)) {
            // COPY_TIMER (aligner_get_word_from_seq);
            return NO_GPOS; // not a A,C,G,T
        }

        *refhash_word |= (base << ((direction==-1 ? -i : i) * 2)); // 2-LSb of word is the first base
    }

    // remove MSb in case of an odd number of base layer bits
    if (bits_per_hash_is_odd) 
        *refhash_word &= layer_bitmask[0]; 

    // Performance note: 50% of the aligner time is taken by memory latency - looking up from 
    // refhashs[] - ~0.28 lookups every base in the sequence.
    PosType64 gpos = (PosType64)BGEN32 (refhashs[0][*refhash_word & layer_bitmask[0]]);
    // COPY_TIMER (aligner_get_word_from_seq);
    return gpos;
}

// converts a string sequence to a 2-bit bitmap
static inline Bits aligner_seq_to_bitmap (VBlockP vb, rom seq, uint64_t seq_len, 
                                          uint64_t *bitmap_words, // allocated by caller
                                          bool *seq_is_all_acgt)  
{
    START_TIMER;

    // covert seq to 2-bit array
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

    bits_clear_excess_bits_in_top_word (&seq_bits); // bc bitmap_words is uninitialized

    COPY_TIMER (aligner_seq_to_bitmap);

    return seq_bits;
}

static inline bool aligner_has_failed (uint32_t match_len, uint32_t seq_len)
{
    return match_len * 50 / seq_len < 73;
}

static inline bool aligner_update_best (VBlockP vb, PosType64 gpos, PosType64 pair_gpos, 
                                        BitsP seq_bits, uint32_t seq_len, bool fwd, bool maybe_perfect_match,
                                        ConstBitsP genome, ConstBitsP emoneg, PosType64 genome_nbases, uint32_t max_snps_for_perfection,
                                        PosType64 *best_gpos, int32_t *longest_len, bool *best_is_forward, bool *is_all_ref) // in/out
{
    // START_TIMER; // this has a small performance impact as it is called in a tight loop - uncomment when needed

    if (gpos == *best_gpos) goto not_near_perfect;
    
    uint64_t genome_start_bit = (fwd ? gpos : genome_nbases-1 - (gpos + seq_bits->nbits/2 -1)) * 2;                            
    
    int32_t distance = bits_hamming_distance (fwd ? genome : emoneg, genome_start_bit, seq_bits, 0, seq_bits->nbits);
    int32_t match_len = (uint32_t)seq_bits->nbits - distance;     
                                                                                                                                                                                                                                                 
    if (pair_gpos != NO_GPOS && ABS(gpos-pair_gpos) > 500) match_len -= 17; // penalty for remote GPOS in 2nd pair
                                                                                                                                
    if (match_len > *longest_len) {                                                                                              
        *longest_len     = match_len;                                                                                            
        *best_gpos       = gpos;                                                                                                 
        *best_is_forward = fwd;                                                                                                
        
        // note: we allow 2 snps and we still consider the match good enough and stop looking further    
        // compared to stopping only if match_len==seq_len, this adds about 1% to the file size, but is significantly faster 
        if (match_len >= (seq_len - max_snps_for_perfection) * 2) { // we found (almost) the best possible match             
            if (maybe_perfect_match && (match_len == seq_len*2)) *is_all_ref = true; // perfect match                 
            // COPY_TIMER (aligner_update_best); 
            return true; // we're done                                                                                                          
        }                                                                                                                       
    }        

not_near_perfect:
    // COPY_TIMER (aligner_update_best); 
    return false; // get other GPOSes matching this refhash_word - in the additional layers 
}                                                                                                                               

// returns gpos aligned with seq with M (as in CIGAR) length, containing the longest match to the reference. 
// returns false if no match found.
// note: matches that imply a negative GPOS (i.e. their beginning is aligned to before the start of the genome), aren't consisdered
static inline PosType64 aligner_best_match (VBlockP vb, STRp(seq), PosType64 pair_gpos,
                                            ConstBitsP genome, ConstBitsP emoneg, PosType64 genome_nbases,
                                            bool *is_forward, bool *is_all_ref) // out
{
    START_TIMER;

    int32_t/*signed*/ longest_len=0; // longest number of bits (not bases!) that match
    uint32_t refhash_word;
    const PosType64 seq_len_64 = (PosType64)seq_len; // 64 bit version of seq_len
    
    PosType64 gpos = NO_GPOS, best_gpos = NO_GPOS; // match not found yet
    bool best_is_forward = false;
    
    // convert seq to a bitmap
    bool maybe_perfect_match;
    uint64_t seq_bits_words[roundup_bits2words64(seq_len * 2)];
    Bits seq_bits = aligner_seq_to_bitmap (vb, seq, seq_len, seq_bits_words, &maybe_perfect_match);
    
    //ref_print_bases (&seq_bits, "\nseq_bits fwd", true);
    //ref_print_bases (&seq_bits, "seq_bits rev", false);

    *is_all_ref = false;

    typedef enum { NOT_FOUND=-1, REVERSE=0, FORWARD=1 } Direction;

    // each "find" corresponds to a place in seq that has the hook base (or if a homopolymer of the hook - the last base in the homopolymer)
    struct Finds { 
        uint32_t refhash_word;  // the sequence (in 2bit) directly before (if fwd) or after (if rev) the hook. Used as a key into refhash. 
        uint32_t i;             // index of hook in seq
        Direction found;        // orientation of find
    } finds[seq_len]; 
    
    uint32_t num_finds = 0;
    
    // in case of --fast, we check only 1/3 of the bases, and we are content with a match (not searching any further) if it 
    // has at most 10 SNPs. On our test file, this reduced the number of calls to aligner_update_best by about 4X, 
    // at the cost of the compressed file being about 11% larger
    uint32_t density = (flag.fast ? 3 : 1);
    uint32_t max_snps_for_perfection = (flag.fast ? 10 : 2); // note: this is only approximately SNPs as we actually measure mismatched bits, not bases

    {START_TIMER;
    // we search - checking both forward hooks and reverse hooks, we check only the first layer for now
    for (uint32_t i=0; i < seq_len; i += density) {          
        Direction found = NOT_FOUND;

        if (i < seq_len - nukes_per_hash && // room for the hash word
            seq[i] == HOOK && seq[i+1] != HOOK &&  // take the G - if there is a homopolymer GGGG... take the last one
            (gpos = aligner_get_word_from_seq (vb, &seq[i+1], &refhash_word, 1)) != NO_GPOS) { 

            gpos -= i; // gpos is the first base on the reference, that aligns to the first base of seq
            found = FORWARD;
        }

        // note: if the previous condition is true, then seq[i]==HOOK and no need to check this condition. hence the "else"
        else if (i >= nukes_per_hash && // room for the hash word
            seq[i] == HOOK_REV && seq[i-1] != HOOK_REV &&  // take the G - if there is a polymer GGGG... take the last one
            (gpos = aligner_get_word_from_seq (vb, &seq[i-1], &refhash_word, -1)) != NO_GPOS) { 

            gpos -= seq_len_64-1 - i; // gpos is the first base of the reference, that aligns wit the LAST base of seq
            found = REVERSE;
        }

        if (found != NOT_FOUND && (gpos >= 0) && (gpos + seq_len_64 < genome_nbases)) { // ignore this gpos if the seq wouldn't fall completely within reference genome
            finds[num_finds++] = (struct Finds){ .refhash_word = refhash_word, .i = i, .found = found };
            
            if (aligner_update_best (vb, gpos, pair_gpos, &seq_bits, seq_len, found, maybe_perfect_match, 
                                     genome, emoneg, genome_nbases, max_snps_for_perfection, 
                                     &best_gpos, &longest_len, &best_is_forward, is_all_ref)) {
                COPY_TIMER(aligner_first_layer);
                goto done; // near-perfect match, search no longer
            }
        }
    }
    COPY_TIMER(aligner_first_layer);}

    // if still no near-perfect matches found, search the additional layers 
    // (each refhash_word can have up to num_layers corresponding GPOSes - we already tested one, now we test the rest)
    {START_TIMER;
    for (uint32_t find_i=0; find_i < num_finds; find_i++) {
        // amount to adjust gpos, which was found for base #i of seq, to correspond to the beginning of seq
        int64_t gpos_shift = (finds[find_i].found == FORWARD ? finds[find_i].i : seq_len_64-1 - finds[find_i].i);

        for (unsigned layer_i=1; layer_i < num_layers; layer_i++) {
            gpos = (PosType64)BGEN32 (refhashs[layer_i][finds[find_i].refhash_word & layer_bitmask[layer_i]]); 

            // case: no more GPOSes for this refhash_word - this layer and all subsequent are NO_GPOS
            if (gpos == NO_GPOS) break;

            gpos -= gpos_shift;

            // check if this gpos results in a better match than the one we already have
            if ((gpos >= 0) && (gpos + seq_len_64 < genome_nbases) &&
                aligner_update_best (vb, gpos, pair_gpos, &seq_bits, seq_len, finds[find_i].found, maybe_perfect_match, 
                                     genome, emoneg, genome_nbases, max_snps_for_perfection, 
                                     &best_gpos, &longest_len, &best_is_forward, is_all_ref)) {
                COPY_TIMER(aligner_additional_layers);
                goto done; // near-perfect match, search no longer       
            }
        }
    }
    COPY_TIMER(aligner_additional_layers);}
    
done:
    *is_forward = best_is_forward;

    // check that the best match is above a threshold, where mapping against the reference actually improves compression
    if (aligner_has_failed (longest_len, seq_len)) 
        best_gpos = NO_GPOS; // note that longest_len==seq_len*2 is perfect match (longest_len is in bits)

    COPY_TIMER (aligner_best_match);
    return best_gpos;
}

MappingType aligner_seg_seq (VBlockP vb, STRp(seq), bool is_pair_2, PosType64 pair_gpos, bool pair_is_forward)
{
    START_TIMER;
    
    declare_seq_contexts;
    ConstBitsP genome, emoneg;
    PosType64 genome_nbases;
    ref_get_genome (gref, &genome, &emoneg, &genome_nbases);

    bool is_forward=false, is_all_ref=false;

    // our aligner algorithm only works for short reads - long reads tend to have many Indel differences (mostly errors) vs the reference
    PosType64 gpos = (seq_len <= MAX_SHORT_READ_LEN) 
        ? aligner_best_match (VB, STRa(seq), pair_gpos, genome, emoneg, genome_nbases, &is_forward, &is_all_ref) : NO_GPOS; 

    if (gpos == NO_GPOS || gpos > genome_nbases - seq_len || gpos > MAX_ALIGNER_GPOS - seq_len || gpos < 0/*never happens*/) {
        COPY_TIMER (aligner_seg_seq);
        return MAPPING_NO_MAPPING;
    }

    if (flag.show_aligner)
        iprintf ("%s: gpos=%"PRId64" forward=%s\n", LN_NAME, gpos, TF(is_forward));

    if ((strand_ctx->local.size & ~3ULL) * 8 == strand_ctx->local.nbits)
        buf_alloc_do (vb, &strand_ctx->local, strand_ctx->local.size + sizeof(uint64_t), CTX_GROWTH, NULL, __FUNCLINE);
        
    buf_alloc (vb, &gpos_ctx->local, 1, 0, uint32_t, CTX_GROWTH, NULL); 

    buf_add_bit (&strand_ctx->local, pair_gpos == NO_GPOS ? is_forward // pair 1 is unaligned - just store the strand
                                                          : (is_forward == pair_is_forward)); // pair 1 is aligned - store equality, expected to he 1 in most cases
            
    if (is_pair_2) 
        fastq_seg_pair2_gpos (vb, pair_gpos, gpos);
    
    // store the GPOS in local if its not a 2nd pair, or if it is, but the delta is not small enough
    else
        BNXT32 (gpos_ctx->local) = (uint32_t)gpos;

    // lock region of reference to protect is_set
    if (IS_REF_EXT_STORE) {
        RefLock lock = ref_lock (gref, gpos, seq_len);
     
        ref_set_genome_is_used (gref, gpos, seq_len); // this region of the reference is used (in case we want to store it with REF_EXT_STORE)
     
        ref_unlock (gref, &lock);
    }

    // shortcut if we have a full reference match
    if (is_all_ref) {
        vb->num_perfect_matches++; // for stats
        
        goto done;
    }

    BitsP bitmap = (BitsP)&bitmap_ctx->local;

    // allocate bitmaps - don't provide name to avoid re-writing param which would overwrite nbits that overlays it 
    int64_t missing_bits = (int64_t)seq_len - ((int64_t)bitmap_ctx->local.nbits - (int64_t)bitmap_ctx->next_local);
    if (missing_bits > 0) {
        buf_alloc (vb, &bitmap_ctx->local, roundup_bits2bytes64 (missing_bits + 64), 
                   seq_len * vb->lines.len32 / 10, // 8 for bits-to-bytes, but up to 10 due to perfect matches 
                   char, CTX_GROWTH, CTX_TAG_LOCAL); 
        buf_extend_bits (&bitmap_ctx->local, missing_bits);
    }

    for (int i=0; i < 4; i++)
        buf_alloc (vb, &seqmis_ctx[i].local, seq_len, 0, char, CTX_GROWTH, NULL); 

    uint64_t next_bit = bitmap_ctx->next_local; // copy to automatic variable (optimized to a register) for performace
    decl_acgt_decode;
    for (uint32_t i=0; i < seq_len; i++) {

        char seq_base = is_forward?seq[i] : seq[i]=='A'?'T' : seq[i]=='T'?'A' : seq[i]=='C'?'G' : seq[i]=='G'?'C' : 0;  
        
        PosType64 ref_i = gpos + (is_forward ? i : seq_len-1-i);
        
        uint8_t ref_base_2bit = bits_get2 (genome, ref_i * 2);
        char ref_base = acgt_decode(ref_base_2bit);

        // case our seq is identical to the reference at this site
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
    vb->num_aligned++; // for stats

    COPY_TIMER (aligner_seg_seq);

    return is_all_ref ? MAPPING_PERFECT : MAPPING_ALIGNED; // successful
}

// PIZ: SEQ reconstruction - only for reads compressed with the aligner
void aligner_reconstruct_seq (VBlockP vb, uint32_t seq_len, bool is_pair_2, bool is_perfect_alignment, ReconType reconstruct,
                              char *first_mismatch_base,       // optional out: caller should initialize to 0
                              uint32_t *first_mismatch_offset, // optional out
                              uint32_t *num_mismatches)        // optional out: caller should initialize to 0
{
    START_TIMER;
    declare_seq_contexts;

    if (!bitmap_ctx->is_loaded) return; // if case we need to skip the SEQ field (for the entire file)
    
    if (is_perfect_alignment || buf_is_alloc (&bitmap_ctx->local)) { // not all non-ref

        bool is_forward;
        PosType64 gpos;

        // first file of a pair (R1) or a non-pair fastq or sam
        if (!is_pair_2) {
            gpos = gpos_ctx->last_value.i = NEXTLOCAL (uint32_t, gpos_ctx);
            is_forward = NEXTLOCALBIT (strand_ctx);        
        }

        // 2nd file of a pair (R2)
        else {
            is_forward = fastq_piz_get_pair2_is_forward (vb); // MUST be called before gpos reconstruction as it inquires GPOS.localR1.next

            // gpos: don't reconstruct just get last_value
            reconstruct_from_ctx (vb, gpos_ctx->did_i, 0, false); // calls fastq_special_PAIR2_GPOS
            gpos = gpos_ctx->last_value.i;
        }

        if (flag.deep) ctx_set_last_value (vb, strand_ctx, (int64_t)is_forward); // consumed by sam_piz_deep_add_seq

        if (flag.show_aligner)
            iprintf ("%s: gpos=%"PRId64" forward=%s\n", LN_NAME, gpos, TF(is_forward));

        const Bits *genome=NULL;
        PosType64 genome_nbases;
        
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

        decl_acgt_decode;

        // faster loop in the common case of a perfect (= no mismatches) alignment (is_perfect_alignment introduced in v14)
        if (is_perfect_alignment) {
            if (reconstruct) {
                char *next = BAFTtxt;

                for (uint32_t i=0; i < seq_len; i++) 
                    *next++ = acgt_decode(BASE_NEXT);                    
                
                Ltxt += seq_len;
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

                    if (first_mismatch_base) {
                        if (! *first_mismatch_base) {
                            *first_mismatch_base   = is_forward ? *BLSTtxt : COMPLEM[(int)*BLSTtxt];
                            *first_mismatch_offset = is_forward ? i : seq_len - i -1;
                        }
                        (*num_mismatches)++;
                    }
                }
            }
        }
    }

    else all_nonref: { 
        if (flag.show_aligner)
            iprintf ("%s: all nonref\n", LN_NAME);

        ASSPIZ (nonref_ctx->next_local + seq_len <= nonref_ctx->local.len32, "NONREF exhausted: next_local=%u + seq_len=%u > local.len=%u",
                nonref_ctx->next_local, seq_len, nonref_ctx->local.len32);

        if (reconstruct) RECONSTRUCT (Bc (nonref_ctx->local, nonref_ctx->next_local), seq_len);
        nonref_ctx->next_local += seq_len;
    }

    ASSPIZ (nonref_ctx->next_local <= nonref_ctx->local.len, "nonref_ctx->next_local=%u is out of range: nonref_ctx->local.len=%u",
            nonref_ctx->next_local, nonref_ctx->local.len32);

    COPY_TIMER (aligner_reconstruct_seq);
}
