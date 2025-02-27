// ------------------------------------------------------------------
//   aligner.c
//   Copyright (C) 2020-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "seg.h"
#include "refhash.h"
#include "file.h"
#include "codec.h"
#include "reconstruct.h"
#include "piz.h"
#include "aligner.h"

// Foward example: If seq is: G-AGGGCT  (G is the hook)  -- matches reference AGGGCT       - function returns 110110101000 (A=00 is the LSb)
// Reverse       : If seq is: CGCCCT-C  (C is the hook)  -- also matches reference AGGGCT  - function returns 110110101000 - the same
// calculates a refhash word from 14 nucleotides following a 'G' (only last G in a sequenece of GGGG...)
static inline PosType64 aligner_get_word_from_seq (VBlockP vb, rom seq, uint32_t *refhash_word, int direction /* 1 forward, -1 reverse */)                                              
{   
    // START_TIMER; // this has a small performance impact as it is called in a tight loop - uncomment when needed

    *refhash_word = 0;

    for (int i=0; direction * i < nukes_per_hash; i += direction) {   
        uint32_t base = nuke_encode_dir (seq[i], direction == 1);
        if (__builtin_expect (base == 4, false)) {
            // COPY_TIMER (aligner_get_word_from_seq);
            return NO_GPOS; // not a A,C,G,T
        }

        *refhash_word |= (base << ((direction==-1 ? -i : i) * 2)); // 2-LSb of word is the first base
    }

    // remove MSb in case of an odd number of base layer bits
    if (bits_per_hash_is_odd) 
        *refhash_word &= layer_bitmask[0]; 

    // prefetch cache for all layers of refhash (this improves aligner performance by ~12% and genozip performance on a simple FASTQ compression by ~8%)
    for (unsigned layer_i=0; layer_i < num_layers; layer_i++) 
        __builtin_prefetch (&refhashs[layer_i][*refhash_word & layer_bitmask[layer_i]], 0/*read-only*/, 1/*best option, empirically*/);
        
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

    bits_clear_excess_bits_in_top_word (&seq_bits, false); // bc bitmap_words is uninitialized

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

        if (__builtin_expect (i < seq_len - nukes_per_hash, true/*probability 93%*/) && // room for the hash word
            __builtin_expect (seq[i] == HOOK,  false/*probability 75%*/) && 
            __builtin_expect (seq[i+1] != HOOK, true/*probability 75%*/) &&  // take the G - if there is a homopolymer GGGG... take the last one
            (gpos = aligner_get_word_from_seq (vb, &seq[i+1], &refhash_word, 1)) != NO_GPOS) { 

            gpos -= i; // gpos is the first base on the reference, that aligns to the first base of seq
            found = FORWARD;
        }

        // note: if the previous condition is true, then seq[i]==HOOK and no need to check this condition. hence the "else"
        else if (__builtin_expect (i >= nukes_per_hash, true/*probability 93%*/) && // room for the hash word
            __builtin_expect (seq[i] == HOOK_REV,  false/*probability 75%*/) && 
            __builtin_expect (seq[i-1] != HOOK_REV, true/*probability 75%*/) &&  // take the G - if there is a polymer GGGG... take the last one
            (gpos = aligner_get_word_from_seq (vb, &seq[i-1], &refhash_word, -1)) != NO_GPOS) { 

            gpos -= seq_len_64-1 - i; // gpos is the first base of the reference, that aligns wit the LAST base of seq
            found = REVERSE;
        }

        if (found != NOT_FOUND && 
            __builtin_expect (gpos >= 0 && gpos + seq_len_64 < genome_nbases, true)) { // ignore this gpos if the seq wouldn't fall completely within reference genome
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
        int64_t gpos_shift = (finds[find_i].found == FORWARD ? (int64_t)finds[find_i].i : (seq_len_64-1 - (int64_t)finds[find_i].i));

        for (unsigned layer_i=1; layer_i < num_layers; layer_i++) {
            gpos = (PosType64)BGEN32 (refhashs[layer_i][finds[find_i].refhash_word & layer_bitmask[layer_i]]); 

            // case: no more GPOSes for this refhash_word - this layer and all subsequent are NO_GPOS
            if (gpos == NO_GPOS) break;

            gpos -= gpos_shift;

            // check if this gpos results in a better match than the one we already have
            if (__builtin_expect ((gpos >= 0) && (gpos + seq_len_64 < genome_nbases), true) &&
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

void aligner_seg_gpos_and_fwd (VBlockP vb, PosType64 gpos, bool is_forward, 
                               bool is_pair_2, PosType64 pair_gpos, bool pair_is_forward)
{
    declare_seq_contexts;

    bool is_interleaved_r2 = (segconf.is_interleaved && vb->line_i % 2);
    ContextP this_strand_ctx = is_interleaved_r2 ? strand_r2_ctx : strand_ctx;

    if (!this_strand_ctx->local.data || (this_strand_ctx->local.size & ~3ULL) * 8 == this_strand_ctx->local.nbits)
        buf_alloc_do (vb, &this_strand_ctx->local, this_strand_ctx->local.size + sizeof(uint64_t), CTX_GROWTH, NULL, __FUNCLINE);

    buf_add_bit (&this_strand_ctx->local, pair_gpos == NO_GPOS ? is_forward // we are not r2, or we are r2 and r1 is unaligned - just store the strand
                                                               : (is_forward == pair_is_forward)); // pair 1 is aligned - store equality, expected to he 1 in most cases

    if (is_pair_2 || is_interleaved_r2) 
        fastq_seg_r2_gpos (vb, pair_gpos, gpos);
    
    // store the GPOS in local if its not a 2nd pair, or if it is, but the delta is not small enough
    else 
        seg_integer (VB, gpos_ctx, gpos, false, 0);
    
    // note regarding interleaved: we always store last_value in gpos/strand (never in gpos_r2/strand_r2)
    gpos_ctx->last_value.i = gpos;
    strand_ctx->last_value.i = is_forward;
}

MappingType aligner_seg_seq (VBlockP vb, STRp(seq), bool is_pair_2, PosType64 pair_gpos, bool pair_is_forward)
{
    START_TIMER;
    
    declare_seq_contexts;
    ConstBitsP genome, emoneg;
    PosType64 genome_nbases;
    ref_get_genome (&genome, &emoneg, &genome_nbases);

    bool is_forward=false, is_all_ref=false;

    // our aligner algorithm only works for short reads - long reads tend to have many Indel differences (mostly errors) vs the reference
    PosType64 gpos = (seq_len <= MAX_SHORT_READ_LEN) 
        ? aligner_best_match (VB, STRa(seq), pair_gpos, genome, emoneg, genome_nbases, &is_forward, &is_all_ref) : NO_GPOS; 

    if (gpos == NO_GPOS || gpos > genome_nbases - seq_len || gpos > MAX_ALIGNER_GPOS - seq_len || gpos < 0/*never happens*/) {
        gpos_ctx->last_value.i = NO_GPOS;

        COPY_TIMER (aligner_seg_seq);
        return MAPPING_NO_MAPPING;
    }

    if (flag.show_aligner)
        iprintf ("%s: gpos=%-10"PRId64" %s %s\n", LN_NAME, gpos, is_forward ? "FWD" : "REV", is_all_ref ? "PERFECT" : "");

    aligner_seg_gpos_and_fwd (vb, gpos, is_forward, is_pair_2, pair_gpos, pair_is_forward);

    // lock region of reference to protect is_set
    if (IS_REF_EXT_STORE) {
        RefLock lock = ref_lock (gpos, seq_len);
     
        ref_set_genome_is_used (gpos, seq_len); // this region of the reference is used (in case we want to store it with REF_EXT_STORE)
     
        ref_unlock (&lock);
    }

    // shortcut if we have a full reference match
    if (is_all_ref) {
        vb->num_perfect_matches++; // for stats
        
        goto done;
    }

    BitsP bitmap = (BitsP)&bitmap_ctx->local;

    buf_alloc_bits (vb, &bitmap_ctx->local, seq_len, vb->lines.len32 / 16 * segconf.std_seq_len, 
                    SET/*initialize to "no mismatches"*/, CTX_GROWTH, CTX_TAG_LOCAL); 

    for (int i=0; i < 4; i++)
        buf_alloc (vb, &seqmis_ctx[i].local, seq_len, 64 KB, char, CTX_GROWTH, NULL); 

    // get textual reference segment (possibly rev-comped)
    char ref_data[4 KB + 8], *ref; // not too big, so "scratch" code path gets milage
    if (seq_len <= 4 KB)
        ref = ref_data + 4; // +4 to avoid underflow in ref_get_textual_seq: save_before assignment (likewise we +4 at the end - total +8)
    else {
        ASSERTNOTINUSE (vb->scratch);
        buf_alloc (vb, &vb->scratch, 0, seq_len, char, 0, "scratch");
        ref = B1STc(vb->scratch);
    }

    ref_get_textual_seq (gpos, ref, seq_len, !is_forward);

    // handle mismatches
    uint64_t bit_i = bitmap_ctx->next_local;  // automatic variable for performace
    // bits_set_region (bitmap, bit_i, seq_len); // initialize to "no mismatches"

    for (uint32_t i=0; i < seq_len; i++, ref++, seq++, bit_i++) 
        if (*seq != *ref) {
            bits_clear (bitmap, bit_i); 
            BNXTc (seqmis_ctx[nuke_encode_dir (*ref, is_forward)].local) = *seq;
        }
    
    bitmap_ctx->next_local += seq_len;

    if (seq_len > sizeof(ref_data)) buf_free (vb->scratch);

done:
    vb->num_aligned++; // for stats

    COPY_TIMER (aligner_seg_seq);

    return is_all_ref ? MAPPING_PERFECT : MAPPING_ALIGNED; // successful
}

void aligner_recon_get_gpos_and_fwd (VBlockP vb, bool is_pair_2, 
                                     PosType64 *gpos, bool *is_forward) // out
{    
    declare_seq_contexts;
    bool is_interleaved_r2 = (segconf.is_interleaved && VB_DT(FASTQ) && vb->line_i % 2); // FASTQ file or FASTQ component of a Deep file

    // first file of a pair (R1) or a non-pair fastq or sam
    if (!is_pair_2 && !is_interleaved_r2) {
        *gpos = gpos_ctx->last_value.i = reconstruct_from_local_int (vb, gpos_ctx, 0, false);
        *is_forward = NEXTLOCALBIT (strand_ctx);  

        if (segconf.is_interleaved) bitmap_ctx->r1_is_aligned = PAIR1_ALIGNED;      
    }

    // 2nd file of a pair (R2) 
    else if (is_pair_2) {
        *is_forward = fastq_piz_get_r2_is_forward (vb); // MUST be called before gpos reconstruction as it inquires GPOS.localR1.next
        reconstruct_from_ctx (vb, gpos_ctx->did_i, 0, false); // calls fastq_special_PAIR2_GPOS
        *gpos = gpos_ctx->last_value.i;
    }
    
    // R2 in a single-file interleaved FASTQ
    else {
        *is_forward = fastq_piz_get_interleaved_r2_is_forward (vb);
        reconstruct_from_ctx (vb, gpos_r2_ctx->did_i, 0, false); // calls fastq_special_PAIR2_GPOS
        *gpos = gpos_r2_ctx->last_value.i;
    }

    ctx_set_last_value (vb, strand_ctx, (int64_t)*is_forward); // consumed by sam_piz_deep_add_seq and next line's fastq_piz_get_interleaved_r2_is_forward
}

// PIZ: SEQ reconstruction - only for reads compressed with the aligner
void aligner_reconstruct_seq (VBlockP vb, uint32_t seq_len, bool is_pair_2, bool is_perfect_alignment, ReconType reconstruct,
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
        aligner_recon_get_gpos_and_fwd (vb, is_pair_2, &gpos, &is_forward);

        if (flag.show_aligner)
            iprintf ("%s: gpos=%"PRId64" forward=%s perfect=%s\n", LN_NAME, gpos, TF(is_forward), TF(is_perfect_alignment));

        if (gpos == NO_GPOS) { 
            bitmap_ctx->next_local += seq_len; // up to v13 ; not sure if this still can happen since v14
            goto all_nonref;
        }

        // reconstruct reference sequence
        ref_get_textual_seq (gpos, BAFTtxt, seq_len, !is_forward);

        // non-perfect alignment: update the mismatches
        if (!is_perfect_alignment) {
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
                    
                    ContextP mis_ctx = VER(14) ? &seqmis_ctx[base_2bit] : nonref_ctx;
                    ASSPIZ (mis_ctx->next_local < mis_ctx->local.len32, NEXT_ERRFMT, __FUNCTION__, mis_ctx->tag_name, mis_ctx->next_local, 1, mis_ctx->local.len32); 

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
        }

        if (reconstruct) Ltxt += seq_len;
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
