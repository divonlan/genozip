// ------------------------------------------------------------------
//   sam_seq.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "sam_private.h"
#include "random_access.h"
#include "aligner.h"
#include "refhash.h"

//---------------
// Shared ZIP/PIZ
//---------------

// called when SEQ in a prim or supp/sec line is segged against ref against prim: 
// Sets vb->mismatch_bases_by_SEQ. ZIP: also sets vb->md_verified and. PIZ: also returns sqbitmap of this line in line_sqbitmap 
static bool sam_analyze_copied_SEQ (VBlockSAMP vb, STRp(seq), const PosType32 pos, bool is_revcomp,
                                    uint32_t ref_consumed, uint32_t ref_and_seq_consumed,
                                    BufferP line_sqbitmap, 
                                    bool has_deep, 
                                    PosType64 *range_last_pos,       // optional out
                                    PosType64 *range_gpos)           // optional out

{
    START_TIMER;

    declare_seq_contexts;
    RefLock lock = REFLOCK_NONE;

    BitsP bitmap = (BitsP)line_sqbitmap;
    uint32_t bit_i=0;

    if (vb->cigar_missing) goto fail;
    
    // we initialize all the bits to "set", and clear as needed.   
    ASSERT (!line_sqbitmap->len32, "%s: line_sqbitmap is in use", LN_NAME);
    buf_alloc_bits_exact (VB, line_sqbitmap, ref_and_seq_consumed, SET, 0, line_sqbitmap->name ? line_sqbitmap->name : "line_sqbitmap"); 

    // if we're going to store seq for Deep, we need to record the NONREF (i.e. S, I) bases
    char *deep_nonref = NULL; 
    if (has_deep) {
        buf_alloc_exact (vb, nonref_ctx->deep_nonref, vb->seq_len - ref_and_seq_consumed, char, "contexts->deep_nonref");
        deep_nonref = is_revcomp ? BLSTc(nonref_ctx->deep_nonref) : B1STc(nonref_ctx->deep_nonref);
    }

    ConstRangeP range = IS_ZIP ? ref_seg_get_range (VB, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_consumed, WORD_INDEX_NONE,  
                                                    IS_REF_EXTERNAL ? NULL : &lock)
                               : ref_piz_get_range (VB, SOFT_FAIL);
    if (!range) goto fail; // can happen in ZIP/REF_INTERNAL due to contig cache contention ; in ZIP/PIZ with an external reference - contig is not in the reference

    if (range_last_pos) *range_last_pos = range->last_pos;
    if (range_gpos)     *range_gpos     = range->gpos;

    uint32_t pos_index = pos - range->first_pos;
    uint32_t next_ref  = pos_index;

    uint32_t i=0, n=0;
    BamCigarOpType op = BC_NONE;

    uint32_t range_len = (range->last_pos - range->first_pos + 1);
    decl_acgt_decode;

    while (i < seq_len || next_ref < pos_index + ref_consumed) {

        ASSERT0 (i <= seq_len && next_ref <= pos_index + ref_consumed, "i or next_ref are out of range");

        if (vb->binary_cigar.next < vb->binary_cigar.len) {
            BamCigarOpP next_op = B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next++);
            op = next_op->op;
            n  = next_op->n;
        }

        switch (op) {

            case BC_M: case BC_E: case BC_X: // alignment match or sequence match or mismatch

                ASSERT (n > 0 && n <= (seq_len - i), "%s: CIGAR implies seq_len longer than actual seq_len=%u (line=%s n=%u). CIGAR=\"%.*s\"", 
                        LN_NAME, seq_len, line_name(VB).s, n, STRfb(vb->textual_cigar));

                for (; n && next_ref < pos_index + ref_consumed; n--, next_ref++, i++) {

                    // circle around to beginning of chrom if out of range (can only happen with external reference, expected only with circular chromosomes) 
                    uint32_t idx = (next_ref >= range_len) ? next_ref - range_len : next_ref; 

                    // case: we don't have the value for the REF site, therefore we can't seg MD against REF (its not worth it to add the REF site just for MD)
                    if (IS_ZIP && (flag.reference & REF_STORED) && !ref_is_nucleotide_set (range, idx)) 
                        goto fail;

                    // case our seq is identical to the reference at this site
                    else if (IS_ACGT (seq[i]) && seq[i] == REF(idx)) 
                        bit_i++; // bit remains set. 
 
                    // case: ref is set to a different value - update the bitmap
                    else {
                        if (has_deep) {
                            if (vb->num_deep_mismatches <= MAX_DEEP_SEQ_MISMATCHES) {
                                vb->deep_mismatch_base[vb->num_deep_mismatches] = seq[i];
                                vb->deep_mismatch_offset[vb->num_deep_mismatches] = i;
                            }
                            vb->num_deep_mismatches++;
                        }

                        bits_clear (bitmap, bit_i++);
                        vb->mismatch_bases_by_SEQ++;
                    } 
                }
                break; // end if 'M', '=', 'X'

            // for Insertion or Soft clipping - this SEQ segment doesn't align with the reference - we leave it as is 
            case BC_I: case BC_S: 
                ASSERT (n > 0 && n <= (seq_len - i), "%s: CIGAR \"%s\" implies seq_len longer than actual seq_len=%u", LN_NAME, 
                        dis_binary_cigar (VB, CIG(vb->binary_cigar), &vb->scratch).s, seq_len);
                
                if (has_deep) 
                    for (uint32_t j=i; j < i + n; j++) {
                        if (is_revcomp) *deep_nonref-- = COMPLEM[(uint8_t)seq[j]];
                        else            *deep_nonref++ = seq[j];
                    }

                i += n;
                n = 0;
                break;

            // for Deletion or Skipping - we move the next_ref ahead
            case BC_D: case BC_N: {
                unsigned ref_consumed_skip = (IS_REF_INTERNAL ? MIN_(n, range_len - next_ref) : n);
                next_ref += ref_consumed_skip;
                n        -= ref_consumed_skip;
                break;
            }

            // Hard clippping (H) or padding (P) - nothing much to do
            case BC_H: case BC_P: 
                n = 0;
                break;

            case BC_NONE:
                ABORT ("%s: End of CIGAR reached but we still have %u reference and %u sequence bases to consume"
                        "(cigar=%s pos=%d seq_len=%u) (ref_consumed=%d next_ref=%u pos_index=%u op.n=%u range=[%.*s %"PRId64"-%"PRId64"])",
                        LN_NAME, pos_index + ref_consumed - next_ref, seq_len-i, vb->last_cigar ? vb->last_cigar : "", pos, seq_len,
                        ref_consumed, next_ref, pos_index, n, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);        

            default:
                ABORT ("%s: Invalid CIGAR op=%u", LN_NAME, op);        
        }

        // case: we're at the end of the reference AND we want more of it
        if (next_ref == pos_index + ref_consumed && n) break;
    }

    ref_unlock (&lock); // does nothing if REFLOCK_NONE

    // an error here can indicate that CIGAR is inconsistent with sequence
    ASSERT (i == seq_len, "%s: expecting i(%u) == seq_len(%u) pos=%d range=[%.*s %"PRId64"-%"PRId64"] (possibly reason: inconsistency between seq_len and CIGAR=\"%s\")", 
            LN_NAME, i, seq_len, pos, STRf(range->chrom_name), range->first_pos, range->last_pos, (vb->last_cigar ? vb->last_cigar : ""));

    // ZIP: final verification step - does MD:Z correctly reflect matches and mismatches of M/X/=
    if (IS_ZIP) {
        sam_MD_Z_verify_due_to_seq (vb, STRa(seq), pos, bitmap, 0); // sets vb->md_verified
        buf_free (*line_sqbitmap);
    }
    
    COPY_TIMER (sam_analyze_copied_SEQ);
    return true;

fail:
    // case ZIP: if we cannot verify against the reference, the MD:Z is not verified, and we don't use the SPECIAL for segging
    ref_unlock (&lock); // does nothing if REFLOCK_NONE
    buf_free (*line_sqbitmap);
    vb->md_verified = false;
    vb->mismatch_bases_by_SEQ = -1; // therefore NM cannot seg against mismatch_bases_by_SEQ
    
    COPY_TIMER (sam_analyze_copied_SEQ);
    return false;
}

//---------
// ZIP
//---------

void sam_seg_SEQ_initialize (VBlockP vb)
{
    declare_seq_contexts;

    bitmap_ctx->ltype = LT_BITMAP;
    strand_ctx->ltype = LT_BITMAP;

    ctx_set_ltype (VB, LT_UINT8, SAM_SEQINS_A, SAM_SEQINS_C, SAM_SEQINS_G, SAM_SEQINS_T,
                   SAM_SEQMIS_A, SAM_SEQMIS_C, SAM_SEQMIS_G, SAM_SEQMIS_T, DID_EOL);

    ctx_set_dyn_int (VB, SAM_GPOS, DID_EOL); // also sets STORE_INT
    
    // MAIN: we may seg depn lines against in-VB prim lines
    if (IS_MAIN(vb) && !flag.bam_assist)
        bitmap_ctx->flags.store_per_line = true; // 14.0.0

    // initial allocations, these might grow during segging if needed
    int factor = segconf.sam_is_unmapped ? 1 : 32; // if mapped, we allocate assuming 1 out of 32 lines is unmapped
    buf_alloc (vb, &bitmap_ctx->local, 1, Ltxt / (4 * factor), uint8_t, 0, CTX_TAG_LOCAL); 
    buf_alloc (vb, &strand_ctx->local, 1, roundup_bits2bytes64 (vb->lines.len / factor), uint8_t, 0, CTX_TAG_LOCAL); 
    buf_alloc (vb, &gpos_ctx->local, 1, vb->lines.len / factor, uint32_t, CTX_GROWTH, CTX_TAG_LOCAL); 

    if (!segconf.sam_is_unmapped || flag.aligner_available) {
        buf_alloc (vb, &nonref_ctx->local, 0, Ltxt / 64, char, 0, CTX_TAG_LOCAL);

        for (int i=0; i < 4; i++)
            buf_alloc (vb, &seqmis_ctx[i].local, 1, Ltxt / 128, char, 0, CTX_TAG_LOCAL); 

        if (segconf.use_insertion_ctxs)    
            for (int i=0; i < 4; i++) 
                buf_alloc (vb, &seqins_ctx[i].local, Ltxt / 256, 0, char, 0, CTX_TAG_LOCAL); 
    }
    else // we store seq vertabim - no reference and no CIGAR
        buf_alloc (vb, &nonref_ctx->local, 0, Ltxt / 3, char, 0, CTX_TAG_LOCAL);
}

// align nonref_ctx->local to a 4-character boundary. this is because CODEC_ACGT squeezes every 4 characters into a byte,
// before compressing it with LZMA. In sorted SAM, we want subsequent identical sequences to have the same byte alignment
// so that LZMA can catch their identicality. (note: as of v13, we assign a codec rather than hard-coding LZMA, but it is usually LZMA anyway)
static void sam_seg_SEQ_pad_nonref (VBlockP vb)
{
    decl_ctx (SAM_NONREF);
    uint32_t add_chars = (4 - (ctx->local.len32 & 3)) & 3;
    if (add_chars) buf_add (&ctx->local, "AAA", add_chars); // add 1 to 3 As
}

// SEG: called from sam_seg_SEQ_vs_ref to handle an M segment of a bisulfite-treated SEQ
static void sam_seg_bisulfite_M (VBlockSAMP vb, 
                                 RangeP range, int32_t range_len, int32_t idx,
                                 STRp(Mseg), uint32_t M_i, // M segment of a SEQ string
                                 BitsP bitmap, uint32_t first_bit,
                                 Context *seqmis_ctx) // 4 contexts
{
    START_TIMER;

    char bisulfite = vb->bisulfite_strand;
    bool MD_NM_by_unconverted = segconf.MD_NM_by_unconverted;
    idx = RR_IDX (idx);
    decl_acgt_decode;

    // iterating on an M segmenet of a SEQ
    for (uint32_t i=0; i < Mseg_len; i++) {
    
        char base = Mseg[i];
        char unconverted_ref_base = REF(idx);
        
        // C->T or G->A conversion
        char converted_ref_base = (bisulfite == unconverted_ref_base) ? ((bisulfite == 'C') ? 'T' : 'A') 
                                                                      : unconverted_ref_base;

        // if needed update unconverted_bitmap - bitmap vs unconverted reference
        if (MD_NM_by_unconverted) {
            vb->unconverted_bitmap.nbits++;
        
            // mismatch vs unconverted
            if (base != unconverted_ref_base) { 
                bits_clear ((BitsP)&vb->unconverted_bitmap, vb->unconverted_bitmap.nbits - 1);
                vb->mismatch_bases_by_SEQ++; // used by NM:i and MD:Z
            }
        }
        
        // case: seq base mismatches converted reference base
        if (base != converted_ref_base) {
            uint8_t ref_base_2bit = acgt_encode[(uint8_t)unconverted_ref_base];

            BNXTc (seqmis_ctx[ref_base_2bit].local) = base;
            bits_clear (bitmap, first_bit + i);

            if (!MD_NM_by_unconverted)
                vb->mismatch_bases_by_SEQ++; // used by NM:i and MD:Z - in case we are calculate MD and NM vs converted reference
        }

        // set methylation call. note: bisulfite is not supported in REF_INTERNAL, so we know range is the entire contig
        if (segconf.sam_predict_meth_call &&
            (bisulfite == unconverted_ref_base && (base == converted_ref_base || base == unconverted_ref_base)))

            sam_bismark_zip_update_meth_call (vb, range, range_len, idx, (bisulfite == base), M_i, Mseg_len, i);

        idx = RR_IDX(idx+1); // round robin around contig (useful in case of circular contigs)
    }

    COPY_TIMER(sam_seg_bisulfite_M);
}

static uint64_t count_non_monochar_ins=0, count_S=0;
static uint64_t count_ins[4] = {};           // monochar insertions by base
static uint64_t count_ins_after [4][4] = {}; // monochar insertions: first index: base before insertion (A,C,G,T) second index: base inserted (A,C,G,T)
static uint64_t count_ins_before[4][4] = {}; // monochar insertions: first index: base after insertion  (A,C,G,T) second index: base inserted (A,C,G,T)

static void sam_seg_analyze_monochar_inserts (VBlockSAMP vb, STRp(seq))
{
    if (vb->cigar_missing) return;

    uint32_t seq_i=0;
    for_cigar (vb->binary_cigar) {        
        case BC_I:
            if (IS_ACGT(seq[seq_i]) && str_is_monochar (&seq[seq_i], op->n)) {
                if (seq_i > 0 && IS_ACGT(seq[seq_i-1]))
                    count_ins_after[acgt_encode[(int8_t)seq[seq_i-1]]][acgt_encode[(int8_t)seq[seq_i]]] += op->n;

                if (seq_i + op->n < seq_len && IS_ACGT(seq[seq_i+op->n]))
                    count_ins_before[acgt_encode[(int8_t)seq[seq_i+op->n]]][acgt_encode[(int8_t)seq[seq_i]]] += op->n;

                count_ins[acgt_encode[(int8_t)seq[seq_i]]] += op->n;
            }

            else
                count_non_monochar_ins += op->n;
            break;

        case BC_S: // the other contributor to NONREF
            count_S += op->n;
            break;

        case BC_M: case BC_E: case BC_X: case BC_D: case BC_N: 
            seq_i += op->n;
            break;

        default: {}
    }
}

// report for --analyze-insertions
void sam_zip_report_monochar_inserts (void)
{
    decl_acgt_decode;
    uint64_t count_monochar_ins = count_ins[0] + count_ins[1] + count_ins[2] + count_ins[3];

    if (!count_monochar_ins) {
        iprintf ("\n%s does not contain any monochar insertions.\n", txt_name);
        return;
    }

    iprintf ("\nSequencing technology: %s\nAligner: %s\n", tech_name(segconf.tech), segconf_sam_mapper_name());

    double total_nonref = count_monochar_ins + count_non_monochar_ins + count_S;
    iprint0 ("\nBreakdown of NONREF data:\n");
    iprintf ("Monochar insertions: %f%% (%"PRIu64")\n", 100.0 * (double)count_monochar_ins / total_nonref, count_monochar_ins);
    iprintf ("Non-monochar insertions: %f%% (%"PRIu64")\n", 100.0 * count_non_monochar_ins / total_nonref, count_non_monochar_ins);
    iprintf ("Soft clips: %f%% (%"PRIu64")\n", 100.0 * count_S / total_nonref, count_S);

             
    iprint0 ("\nBreakdown of monochar insertions by base: (monochar = single base or homopolymer)\n");
    for (int b=0; b < 4; b++)
        iprintf ("%c: %f%% (%"PRIu64"); of them: same_as_base_before=%f%% same_as_base_after=%f%%\n", 
                 acgt_decode(b), 100.0 * (double)count_ins[b] / (double)count_monochar_ins, count_ins[b],
                 count_ins[b] ? 100.0 *(double)count_ins_after[b][b] / (double)count_ins[b] : 0.0,
                 count_ins[b] ? 100.0 *(double)count_ins_before[b][b] / (double)count_ins[b] : 0.0);

    iprint0 ("\nBreakdown of monochar insertions by base and PREVIOUS base in sequence:\n");
    for (int prev_b=0; prev_b < 4; prev_b++)
        for (int b=0; b < 4; b++)
            iprintf ("prev=%c insertion=%c: %f%% (%"PRIu64")\n", 
                     acgt_decode(prev_b), acgt_decode(b), 100.0 * (double)count_ins_after[prev_b][b] / (double)count_monochar_ins, count_ins_after[prev_b][b]);

    iprint0 ("\nBreakdown of monochar insertions by base and NEXT base in sequence:\n");
    for (int next_b=0; next_b < 4; next_b++)
        for (int b=0; b < 4; b++)
            iprintf ("next=%c insertion=%c: %f%% (%"PRIu64")\n", 
                     acgt_decode(next_b), acgt_decode(b), 100.0 * (double)count_ins_before[next_b][b] / (double)count_monochar_ins, count_ins_before[next_b][b]);

    // reset for next file
    memset (count_ins, 0, sizeof (count_ins));
    memset (count_ins_after, 0, sizeof (count_ins_after));
    memset (count_ins_before, 0, sizeof (count_ins_before));   
    count_non_monochar_ins = count_S = 0; 
}

// true if next seq-consuming op is I
static inline bool next_op_is_I (VBlockSAMP vb, uint32_t this_op_i)
{
    for (uint32_t op_i=this_op_i+1; op_i < vb->binary_cigar.len32; op_i++) {
        BamCigarOpType op = B(BamCigarOp, vb->binary_cigar, op_i)->op;
        if (op == BC_I) 
            return true; // next seq-consuming op is indeed an I
        else if (op == BC_M || op == BC_S || op == BC_E || op == BC_X)
            return false; // next seq-consuming op is not an I
    }
    
    return false; // there is no further seq-consuming op
}

// Creates a bitmap from seq data - exactly one bit per base that consumes reference (i.e. not bases with CIGAR I/S)
// - Normal SEQ: tracking CIGAR, we compare the sequence to the reference, indicating in a SAM_SQBITMAP whether this
//   base in the same as the reference or not. In case of REF_INTERNAL, if the base is not already in the reference, we add it.
//   bases that differ from the reference are stored in SAM_NONREF
// - Edge case: no POS (i.e. unaligned read) - we just store the sequence in SAM_NONREF (unless we can use an aligner)
// - Edge case: no CIGAR (it is "*") - we just store the sequence in SAM_NONREF
// - Edge case: no SEQ (it is "*") - we '*' in SAM_NONREF and indicate "different from reference" in the bitmap. We store a
//   single entry, regardless of the number of entries indicated by CIGAR
static MappingType sam_seg_SEQ_vs_ref (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(seq), const PosType32 pos, bool is_revcomp,
                                       uint32_t ref_consumed, uint32_t ref_and_seq_consumed)
{
    declare_seq_contexts;   
    START_TIMER;

    if (flag.analyze_ins) 
        sam_seg_analyze_monochar_inserts (vb, STRa(seq));

    BitsP bitmap = (BitsP)&bitmap_ctx->local;
    uint32_t bitmap_start = bitmap_ctx->next_local;        

    ASSERTW (seq_len < 100000 || segconf_running || segconf.is_long_reads, 
             "%s: Warning: sam_seg_SEQ: seq_len=%u is suspiciously high and might indicate a bug", LN_NAME, seq_len);

    // we don't need to lock if the entire ref_consumed of this read was already is_set by previous reads (speed optimization)
    bool no_lock = IS_REF_EXTERNAL || ((vb->chrom_node_index == vb->consec_is_set_chrom) && (pos >= vb->consec_is_set_pos) && (pos + ref_consumed <= vb->consec_is_set_pos + vb->consec_is_set_len));

    RefLock lock = REFLOCK_NONE;
    RangeP range = vb->cigar_missing ? NULL : ref_seg_get_range (VB, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_consumed, WORD_INDEX_NONE, (no_lock ? NULL : &lock));

    // Cases where we don't consider the refernce and just copy the seq as-is
    // 1. (denovo:) this contig defined in @SQ went beyond the maximum genome size of 4B and is thus ignored
    // 2. (loaded:) case contig doesn't exist in the reference, or POS is out of range of contig (observed in the wild with chrM)    
    if (!range) {
        COPY_TIMER (sam_seg_SEQ_vs_ref);
        return MAPPING_NO_MAPPING; // seg as a verbatim copy or use aligner
    }

    buf_alloc_bits (vb, &bitmap_ctx->local, ref_and_seq_consumed, vb->lines.len32 / 16 * segconf.std_seq_len, 
                    SET/*initialize to "no mismatches"*/, CTX_GROWTH, CTX_TAG_LOCAL); 
    
    if (vb->bisulfite_strand) {
        if (segconf.sam_predict_meth_call) {
            buf_alloc_exact (vb, vb->meth_call, seq_len, char, "meth_call");
            memset (vb->meth_call.data, '.', vb->meth_call.len32);
        }

        if (segconf.MD_NM_by_unconverted) {
            // we initialize all the bits to "set", and clear as needed.
            buf_alloc_bits_exact (vb, &vb->unconverted_bitmap, ref_and_seq_consumed, SET, 0, "unconverted_bitmap");
            vb->unconverted_bitmap.nbits = 0; // we will grow it back to ref_and_seq_consumed as we go along
        }
    }

    for (int i=0; i < 4; i++) 
        buf_alloc (vb, &seqmis_ctx[i].local, ref_and_seq_consumed, 0, char, CTX_GROWTH, CTX_TAG_LOCAL); 

    if (segconf.use_insertion_ctxs)    
        for (int i=0; i < 4; i++) 
            buf_alloc (vb, &seqins_ctx[i].local, vb->insertions,   0, char, CTX_GROWTH, CTX_TAG_LOCAL); 

    buf_alloc (vb, &nonref_ctx->local, seq_len + 3, 0, uint8_t, CTX_GROWTH, CTX_TAG_LOCAL); 

    bitmap_ctx->local_num_words++;

    uint32_t pos_index = pos - range->first_pos;
    uint32_t next_ref  = pos_index;

    uint32_t i=0, n=0;
    BamCigarOpType op = BC_NONE;

    uint32_t range_len = (range->last_pos - range->first_pos + 1);
    
    if (IS_REF_EXT_STORE && !no_lock) 
        bits_set_region (&range->is_set, pos_index, ref_consumed); // we will need this ref to reconstruct

    bool has_D_N = false;
    decl_acgt_decode;

    bool is_ref_internal = IS_REF_INTERNAL;

    while (i < seq_len || next_ref < pos_index + ref_consumed) {

        ASSERT0 (i <= seq_len && next_ref <= pos_index + ref_consumed, "i or next_ref are out of range");

        if (vb->binary_cigar.next < vb->binary_cigar.len) {
            BamCigarOp *next_op = B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next++);
            op = next_op->op;
            n  = next_op->n;
        }

        switch (op) {

            case BC_M: case BC_E: case BC_X: // alignment match or sequence match or mismatch

                ASSERT (n > 0 && n <= (seq_len - i), 
                        "%s: CIGAR implies seq_len longer than actual seq_len=%u (n=%u ). CIGAR=\"%s\"", 
                        LN_NAME, seq_len, n, vb->last_cigar ? vb->last_cigar : "");

                uint32_t bit_i = bitmap_ctx->next_local; // copy to automatic variable for performance
                uint32_t save_n = n;

                if (vb->bisulfite_strand) {
                    sam_seg_bisulfite_M (vb, range, range_len, next_ref, &seq[i], n, i, bitmap, bit_i, seqmis_ctx);
                    next_ref += n;
                    i += n;
                    n = 0;
                }
                
                else while (n && next_ref < pos_index + ref_consumed) {

                    bool normal_base = IS_ACGT (seq[i]);

                    // circle around to beginning of chrom if out of range (can only happen with external reference, expected only with circular chromosomes) 
                    uint32_t idx = (next_ref >= range_len) ? next_ref - range_len : next_ref; 

                    // case: we have not yet set a value for this site - we set it now. note: in ZIP, is_set means that the site
                    // will be needed for pizzing. With REF_INTERNAL, this is equivalent to saying we have set the ref value for the site
                    if (is_ref_internal && !no_lock && !ref_is_nucleotide_set (range, idx)) { 

                        // note: in case this is a non-normal base (eg N), set the reference to an arbitrarily to 'A' as we 
                        // we will store this non-normal base in seqmis_ctx multiplexed by the reference base (i.e. in seqmis_ctx['A']).
                        ref_set_nucleotide (range, idx, normal_base ? seq[i] : 'A');
                        
                        bits_set (&range->is_set, idx); // we will need this ref to reconstruct

                        if (normal_base) 
                            bit_i++; 
                        else
                            goto mismatch;
                    }

                    // case: our seq is identical to the reference at this site
                    else if (normal_base && seq[i] == REF(idx)) 
                        bit_i++; // bit remains set. 

                    // case: ref is set to a different value - we store our value in seqmis_ctx
                    else mismatch: {
                        uint8_t ref_base_2bit = bits_get2 (&range->ref, idx * 2);
                        BNXTc (seqmis_ctx[ref_base_2bit].local) = seq[i];
                        bits_clear (bitmap, bit_i++);
                        vb->mismatch_bases_by_SEQ++; // used by NM:i and MD:Z
                    } 
  
                    n--;
                    next_ref++;
                    i++;
                }
                
                bitmap_ctx->next_local += save_n - n;
                break; // end if 'M', '=', 'X'

            // for Insertion or Soft clipping - this SEQ segment doesn't align with the reference - we leave it as is 
            case BC_I: case BC_S:
                ASSSEG (n > 0 && n <= (seq_len - i), 
                        "CIGAR %s implies seq_len longer than actual seq_len=%u", vb->last_cigar ? vb->last_cigar : "", seq_len);

                if (op == BC_I && segconf.use_insertion_ctxs) {
                    // mux by base after insertion in SEQ. 
                    // note: if i+n is the end of SEQ, then the byte after is 0 in BAM (see bam_seq_to_sam) and \t in SAM - both mapped to 0.
                    // note: if the next seq-consuming op is also I (e.g. in RNA: 621I2611N310I), map to 0 as reconstruct won't yet have seq[i+n] (since 15.0.61 - see defect 2024-06-16)
                    int ins_ctx_i = next_op_is_I (vb, vb->binary_cigar.next - 1) ? 0 : acgt_encode[(int8_t)seq[i + n]];

                    buf_add (&(seqins_ctx + ins_ctx_i)->local, &seq[i], n);
                }
                else 
                    buf_add (&nonref_ctx->local, &seq[i], n);

                i += n;
                n = 0;
                break;

            // for Deletion or Skipping - we move the next_ref ahead
            case BC_D: case BC_N: {
                unsigned ref_consumed_skip = (is_ref_internal ? MIN_(n, range_len - next_ref) : n);
                next_ref += ref_consumed_skip;
                n        -= ref_consumed_skip;
                has_D_N  =  true;
                break;
            }

            // Hard clippping (H) or padding (P) - nothing much to do
            case BC_H: case BC_P: 
                n = 0;
                break;

            case BC_NONE:
                ABOSEG ("End of CIGAR reached but we still have %u reference and %u sequence bases to consume"
                        "(cigar=%s pos=%d seq_len=%u) (ref_consumed=%d next_ref=%u pos_index=%u op.n=%u range=[%.*s %"PRId64"-%"PRId64"])",
                        pos_index + ref_consumed - next_ref, seq_len-i, vb->last_cigar, pos, seq_len,
                        ref_consumed, next_ref, pos_index, n, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);        

            default:
                ABOSEG ("Invalid CIGAR op=%u", op);        
        }

        // case: we're at the end of the reference AND we want more of it
        if (next_ref == pos_index + ref_consumed && n) break;
    }

    ref_unlock (&lock); // does nothing if REFLOCK_NONE

    // an error here can indicate that CIGAR is inconsistent with sequence
    ASSERT (i == seq_len, "%s: expecting i(%u) == seq_len(%u) pos=%d range=[%.*s %"PRId64"-%"PRId64"] (possibly reason: inconsistency between seq_len and CIGAR=\"%s\")", 
            LN_NAME, i, seq_len, pos, STRf(range->chrom_name), range->first_pos, range->last_pos, (vb->last_cigar ? vb->last_cigar : ""));

    // if we set the entire consecutive reference range covered by this read - extend consec_is_set 
    if (segconf.is_sorted && !no_lock && (IS_REF_EXT_STORE || !has_D_N/*REF_INTERNAL*/)) {
        // case: current region is not consecutive - start a new region
        if (vb->consec_is_set_chrom != vb->chrom_node_index || vb->consec_is_set_pos + vb->consec_is_set_len < pos) {
            vb->consec_is_set_chrom = vb->chrom_node_index;
            vb->consec_is_set_pos   = pos;
            vb->consec_is_set_len   = ref_consumed;
        }

        // case: extend current region
        else if (vb->consec_is_set_pos + vb->consec_is_set_len < pos + ref_consumed)
            vb->consec_is_set_len = (pos + ref_consumed) - vb->consec_is_set_pos;
    }
    
    // verify - does MD:Z correctly reflect matches and mismatches of M/X/=
    bool use_un = segconf.MD_NM_by_unconverted && vb->bisulfite_strand;
    sam_MD_Z_verify_due_to_seq (vb, STRa(seq), pos, 
                                use_un ? (BitsP)&vb->unconverted_bitmap : bitmap, 
                                use_un ? 0                              : bitmap_start);

    // note: change of logic in v14: up to v13, we use to align every recursion level, it seems that this was a bug
    sam_seg_SEQ_pad_nonref (VB);

    vb->num_seq_by_aln++; // for stats

    COPY_TIMER (sam_seg_SEQ_vs_ref);

    if (vb->mismatch_bases_by_SEQ != 0 || vb->bisulfite_strand) // we don't support MAPPING_PERFECT with bisulfite (bug 649)
        return MAPPING_ALIGNED;
    
    else {
        bitmap_ctx->next_local -= vb->ref_and_seq_consumed; // note: we truncate bitmap, if needed, in sam_seg_finalize
        return MAPPING_PERFECT; // no mismatches
    } 
}

// verify that the depn line is identical (modulo rev_comp and hard_clips) to the prim line
static bool sam_seg_verify_saggy_line_SEQ (VBlockSAMP vb, ZipDataLineSAMP my_dl, ZipDataLineSAMP prim_dl, rom my_seq/*textual*/)
{
    START_TIMER;

    bool xstrand = (my_dl->FLAG.rev_comp != prim_dl->FLAG.rev_comp); // primary and dependent are on opposite strands
    uint32_t seq_len = my_dl->SEQ.len;

    rom prim_seq;
     // TO DO: compare BAM native-to-native without converting to textual
    if (IS_BAM_ZIP) {
        ASSERTNOTINUSE (vb->scratch);
        buf_alloc (vb, &vb->scratch, 0, prim_dl->SEQ.len + 2/* +1 for last half-byte and \0 */, char, 0, "scratch");
        bam_seq_to_sam (VB, (bytes)STRtxt (prim_dl->SEQ), false, false, &vb->scratch, false);
        prim_seq = Bc(vb->scratch, vb->hard_clip[xstrand]);
    }
    else
        prim_seq = Btxt (prim_dl->SEQ.index) + vb->hard_clip[xstrand];

    bool success = true;
    if (!xstrand)
        success = !memcmp (my_seq, prim_seq, seq_len);

    else  {
        rom my_next = my_seq;
        rom prim_next = prim_seq + seq_len - 1;
        for (uint32_t i=0; i < seq_len; i++)
            if (*my_next++ != COMPLEM[(uint8_t)*prim_next--]) {
                success = false;
                break;
            }
    }

    if (IS_BAM_ZIP)
        buf_free (vb->scratch);

    COPY_TIMER(sam_seg_verify_saggy_line_SEQ);
    return success;
}

void sam_seg_SEQ (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(textual_seq), unsigned add_bytes)
{
    START_TIMER;

    declare_seq_contexts;   
    bool force_sequence = false; // true if reconstructor would normally diff, but in this case it should reconstruct from SQBITMAP.local
    bool force_verbatim = false; // true if reconstructor would normally be vs_ref, but it verbatim
    bool aligner_used = false;
    bool perfect = false;
    bool vs_prim = false;

    ZipDataLineSAMP saggy_dl;
    bool unmapped = dl->FLAG.unmapped || vb->cigar_missing || !dl->POS || str_issame_(STRa(vb->chrom_name), "*", 1);

    ASSSEG (textual_seq[0] != '*' || textual_seq_len == 1, "Invalid SEQ - starts with '*' but length=%u", textual_seq_len);
    
    // case segconf: we created the contexts for segconf_set_vb_size accounting. that's enough - actually segging will mark is_set and break sam_seg_MD_Z_analyze.
    if (segconf_running) {
        if ((vb->bisulfite_strand=='G') != dl->FLAG.rev_comp)
            segconf.bs_strand_not_by_rev_comp = true; // note: in files that bisulfite_strand is determined by rev_comp, it is true for all lines

        if (!vb->line_i)
            segconf_mark_as_used (VB, 4 + !segconf.sam_is_unmapped, SAM_SQBITMAP, SAM_NONREF_X, SAM_GPOS, SAM_STRAND, SAM_NONREF);

        if (segconf.sam_is_unmapped)
            goto add_seq_verbatim; // used to test multiseq in segconf_finalize

        return; 
    }

    if (textual_seq[0] == '*')
        vb->seq_missing = dl->no_seq = true;

    if (flag.deep || flag.show_deep == SHOW_DEEP_ONE_HASH)
        sam_deep_set_SEQ_hash (vb, dl, STRa(textual_seq));

    bool aligner_ok = flag.aligner_available && !segconf_running;

    // case: unmapped line and we have refhash: align to reference
    if (unmapped && aligner_ok) {
        use_aligner:
        switch (aligner_seg_seq (VB, STRa(textual_seq), false, NO_GPOS, false)) {
            case MAPPING_NO_MAPPING : force_verbatim = true; goto add_seq_verbatim; 
            case MAPPING_PERFECT    : perfect        = true; // fallthrough
            case MAPPING_ALIGNED    : aligner_used   = true; break;
            default                 : ABORT0 ("bad value");
        }

        buf_alloc (vb, &nonref_ctx->local, 3, 0, uint8_t, CTX_GROWTH, CTX_TAG_LOCAL); 
        sam_seg_SEQ_pad_nonref (VB);

        vb->md_verified = false;    
        vb->mismatch_bases_by_SEQ = -1; // we can't seg NM:i with special
    }

    // case: unmapped line, no refhash: just store the sequence in nonref without an indication in the bitmap
    else if (unmapped && !aligner_ok) { 
        add_seq_verbatim:        
        buf_alloc (vb, &nonref_ctx->local, textual_seq_len + 3, 0, uint8_t, CTX_GROWTH, CTX_TAG_LOCAL); 
        memcpy (BAFTc(nonref_ctx->local), textual_seq, textual_seq_len);
        nonref_ctx->local.len32 += textual_seq_len;

        sam_seg_SEQ_pad_nonref (VB);

        vb->md_verified = false;    
        vb->mismatch_bases_by_SEQ = -1; 
        force_verbatim = true;    
        vb->num_verbatim++; // for stats
    }

    // case: no SEQ. we already add '-' to CIGAR - no data added here
    else if (textual_seq[0] == '*') {
        vb->md_verified = false;    
        vb->mismatch_bases_by_SEQ = -1; 
    }

    // in case of DEPN line with SEQ confirmed by sam_seg_depn_is_subseq_of_prim to be a subset of the sag sequence, 
    else if (IS_DEPN(vb) && vb->sag) { // depn 
        vs_prim = true;
        vb->num_vs_prim++;
    }
    
    // case: DEPN line (or STAR line - see sam_seg_buddy) vs same-VB prim
    else if (({ saggy_dl = DATA_LINE (vb->saggy_line_i) ; sam_has_saggy && !saggy_dl->hard_clip[0] && !saggy_dl->hard_clip[1]; })) {
        
        // note: these conditions must generate "force_sequence" because PIZ cannot check for them
        if (saggy_dl->no_seq || 
            (textual_seq_len + vb->hard_clip[0] + vb->hard_clip[1] != saggy_dl->SEQ.len) ||
            !sam_seg_verify_saggy_line_SEQ (vb, dl, saggy_dl, textual_seq)) {

            force_sequence = true;    
            goto vs_sequence;
        }
        
        vs_prim = true;
        vb->num_vs_prim++;
    }
    
    else vs_sequence: 
        switch (sam_seg_SEQ_vs_ref (vb, dl, STRa(textual_seq), dl->POS, dl->FLAG.rev_comp, vb->ref_consumed, vb->ref_and_seq_consumed)) {
            case MAPPING_NO_MAPPING : if (flag.aligner_available) goto use_aligner;
                                      else { force_verbatim = true; goto add_seq_verbatim; };

            case MAPPING_PERFECT    : perfect = true; break;

            default                 : break;
        }

    // set vb->md_verified and vb->mismatch_bases_by_SEQ needed by MD:Z and NM:i
    bool force_no_analyze_depn_SEQ = false;
    if (vs_prim && segconf.has_MD_or_NM && !vb->cigar_missing && !vb->seq_missing &&
        !sam_analyze_copied_SEQ (vb, STRa(textual_seq), dl->POS, dl->FLAG.rev_comp, vb->ref_consumed, vb->ref_and_seq_consumed, &vb->codec_bufs[0], false, NULL, NULL)) 
            force_no_analyze_depn_SEQ = true;


    // case: PRIM - SQBITMAP context is reconstructed when loading SA Groups (preprocessing) - invoking SPECIAL_SEQ and reconstructing the sequence
    //       and SQBITMAP is reconstructed again during recon - but this SPECIAL_SEQ copies from SA Groups
    // case: MAIN/DEPN - SPECIAL_SEQ is invoked by SQBITMAP - it either reconstructs sequence or diffs vs prim.
    //       (prim SEQ may be taken from SA Groups or from earlier in the current VB (the latter only in MAIN VBs)
    #define SEQ_SNIP_FORCE_SEQ           0
    #define SEQ_SNIP_ALIGNER_USER        1
    #define SEQ_SNIP_PERFECT             2
    #define SEQ_SNIP_NO_ANALYZE_DEPN_SEQ 3
    #define SEQ_SNIP_BISULFATE           4
    #define SEQ_FORCE_VERBATIM           5

    seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SEQ, '0'+force_sequence, '0'+aligner_used, '0'+perfect, '0'+force_no_analyze_depn_SEQ, 
                                (!vb->bisulfite_strand) ? '0'
                              : (vb->bisulfite_strand && !segconf.bs_strand_not_by_rev_comp && (vb->bisulfite_strand=='G') == dl->FLAG.rev_comp) ? '*' // bisulfite_strand is predicted by rev_comp
                              :                           vb->bisulfite_strand,
                              '0'+force_verbatim },
                8, SAM_SQBITMAP, add_bytes); 

    COPY_TIMER (sam_seg_SEQ);
}

// converts native SAM/BAM format to 2bit ACGT - if soft_fail, returns false if any base is not A,C,G or T 
bool sam_seq_pack (VBlockSAMP vb, Bits *packed, uint64_t next_bit, STRp(seq), bool bam_format, bool revcomp, FailType soft_fail)
{
    if (bam_format) {
        if (!revcomp)
            for (uint32_t i=0; i < seq_len; i++, next_bit += 2) {
                uint8_t b = (!(i&1)) ? (((uint8_t*)seq)[i>>1] >> 4) : (((uint8_t*)seq)[i>>1] & 0xf);
                switch (b) {
                    case 0b0001 : bits_assign2 (packed, next_bit, 0); break;
                    case 0b0010 : bits_assign2 (packed, next_bit, 1); break;
                    case 0b0100 : bits_assign2 (packed, next_bit, 2); break;
                    case 0b1000 : bits_assign2 (packed, next_bit, 3); break;
                    default     : if (soft_fail) return false;
                                  ABORT ("%s: Unexpected base: '%c' i=%u seq_len=%u", LN_NAME, bam_base_codes[b], i, seq_len);
                }
            }    
        else
            for (int32_t i=seq_len-1; i >= 0; i--, next_bit += 2) {
                uint8_t b = (!(i&1)) ? (((uint8_t*)seq)[i>>1] >> 4) : (((uint8_t*)seq)[i>>1] & 0xf);
                switch (b) {
                    case 0b0001 : bits_assign2 (packed, next_bit, 3); break;
                    case 0b0010 : bits_assign2 (packed, next_bit, 2); break;
                    case 0b0100 : bits_assign2 (packed, next_bit, 1); break;
                    case 0b1000 : bits_assign2 (packed, next_bit, 0); break;
                    default     : if (soft_fail) return false;
                                  ABORT ("%s: Unexpected base: '%c' i=%u seq_len=%u", LN_NAME, bam_base_codes[b], i, seq_len);
                }
            }    
    }

    else { // SAM
        if (!revcomp)
            for (uint32_t i=0; i < seq_len; i++, next_bit += 2) 
                switch (seq[i]) {
                    case 'A' : bits_assign2 (packed, next_bit, 0); break;
                    case 'C' : bits_assign2 (packed, next_bit, 1); break;
                    case 'G' : bits_assign2 (packed, next_bit, 2); break;
                    case 'T' : bits_assign2 (packed, next_bit, 3); break;
                    default  : if (soft_fail) return false;
                               ABORT ("%s: Unexpected base: '%s'(ASCII %u) i=%u seq_len=%u", LN_NAME, char_to_printable(seq[i]).s, (uint8_t)seq[i], i, seq_len);
                }
        else
            for (int32_t i=seq_len-1; i >= 0; i--, next_bit += 2) 
                switch (seq[i]) {
                    case 'A' : bits_assign2 (packed, next_bit, 3); break;
                    case 'C' : bits_assign2 (packed, next_bit, 2); break;
                    case 'G' : bits_assign2 (packed, next_bit, 1); break;
                    case 'T' : bits_assign2 (packed, next_bit, 0); break;
                    default  : if (soft_fail) return false;
                               ABORT ("%s: Unexpected base: '%s'(ASCII %u) i=%u seq_len=%u", LN_NAME, char_to_printable(seq[i]).s, (uint8_t)seq[i], i, seq_len);
                }
    }

    bits_clear_excess_bits_in_top_word (packed, true);

    return true;
}

// setting ref bases by analyze_* functions
rom ERR_ANALYZE_RANGE_NOT_AVAILABLE = "Range not available";
rom ERR_ANALYZE_DEPN_NOT_IN_REF     = "Depn alignment base not is_set in reference";
rom ERR_ANALYZE_INCORRECT_REF_BASE  = "incorrect reference base";
rom sam_seg_analyze_set_one_ref_base (VBlockSAMP vb, bool is_depn, PosType32 pos, char base, 
                                      uint32_t ref_consumed, // remaining ref_consumed starting at pos                                      
                                      RangeP *range_p, RefLock *lock)
{
    // case: pos is beyond the existing range
    if ((*range_p) && (*range_p)->last_pos < pos) {
        ref_unlock (lock);
        *range_p = NULL;
    }

    // get range (and lock it if needed)
    if (! *range_p) {
        *range_p = ref_seg_get_range (VB, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_consumed, WORD_INDEX_NONE,  
                                      (IS_REF_EXTERNAL || is_depn) ? NULL : lock);
        if (! *range_p) return ERR_ANALYZE_RANGE_NOT_AVAILABLE; // cannot access this range in the reference
    }

    uint32_t pos_index = pos - (*range_p)->first_pos; // index within range

    // case: depn line, but we haven't set the reference data for it. since we don't need the reference data for reconstructing
    // SEQ (it is copied from prim), we don't store it just for MD:Z - so we won't use SPECIAL for this MD:Z
    if (is_depn && (flag.reference & REF_STORED) && !ref_is_nucleotide_set (*range_p, pos_index))
        return ERR_ANALYZE_DEPN_NOT_IN_REF;

    bool internal_pos_is_populated = IS_REF_INTERNAL && ref_is_nucleotide_set (*range_p, pos_index);

    // case: reference already contains a base - but unfortunately it is not "base" - so MD is not reconstractable from reference
    decl_acgt_decode;
    if (((flag.reference & REF_ZIP_LOADED) || internal_pos_is_populated) && (base != ref_base_by_pos (*range_p, pos))) 
        return ERR_ANALYZE_INCORRECT_REF_BASE; // encountered in the wild when the reference base is a IUPAC

    // case: reference is not set yet - set it now (extra careful never to set anything in depn, as range is not locked)
    if (!is_depn && IS_REF_INTERNAL && !internal_pos_is_populated)
        ref_set_nucleotide (*range_p, pos_index, base);

    // set is_set - we will need this base in the reference to reconstruct MD
    if (!is_depn && flag.reference & REF_STORED)
        bits_set (&(*range_p)->is_set, pos_index); // we will need this ref to reconstruct

    return NULL; // success
}

uint32_t sam_zip_get_seq_len (VBlockP vb, uint32_t line_i) 
{ 
    return DATA_LINE (line_i)->SEQ.len;
}

// used by codec_longr, codec_homp, codec_pacb and codec_smux
COMPRESSOR_CALLBACK_DT (sam_zip_seq) 
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    ZipDataLineSAMP dl = DATA_LINE (vb_line_i);
        
    *line_data_len = dl->SEQ.len;

    if (__builtin_expect (VB_DT(SAM), false)) 
        *line_data = Btxt (dl->SEQ.index);
    
    else { // BAM
        vb->textual_seq.len = 0; // reset
        buf_alloc (vb, &vb->textual_seq, 0, dl->SEQ.len+2 /* +1 for last half-byte and \0 */, char, 1.5, "textual_seq");

        if (__builtin_expect (dl->no_seq, false)) 
            bam_seq_to_sam (VB, 0, 0, false, false, &vb->textual_seq, true);
        else
            bam_seq_to_sam (VB, (uint8_t *)STRtxt(dl->SEQ), false, false, &vb->textual_seq, true);
        
        *line_data = B1STc (VB_SAM->textual_seq);
    }

    if (is_rev) *is_rev = dl->FLAG.rev_comp;

    if ((*line_data)[0] == '*')
        *line_data_len = 1; // note: if we have a CIGAR but SEQ=*, SEQ.len is determined by the CIGAR
}

//---------
// PIZ
//---------

// get a bytemap of ref_consumed values. returns NULL if no range exists. get 2 bases before and after
// in case needed for methylation calling.
static inline rom sam_reconstruct_SEQ_get_textual_ref (VBlockSAMP vb, bool v14, PosType32 pos, bool predict_meth_call, bool to_txt_data)
{
    START_TIMER;

    declare_seq_contexts;

    // note: in an edge case, when all is_set bits are zero, the range might not even be written to the file
    bool uses_ref_data = vb->ref_consumed && 
                            (v14 || // v14: we always have is_set for all the ref data (mismatches use it for SEQMIS)
                             bits_num_set_bits_region ((BitsP)&bitmap_ctx->local, bitmap_ctx->next_local, vb->ref_and_seq_consumed) > 0);

    vb->range = uses_ref_data ? (RangeP)ref_piz_get_range (VB, SOFT_FAIL) : NULL; 
    if (!vb->range) return NULL; 

    PosType32 range_len = vb->range ? (vb->range->last_pos - vb->range->first_pos + 1) : 0;
    uint32_t num_bases = vb->ref_consumed + (predict_meth_call ? 4 : 0); // 2 bases before and after in case needed for methylation

    char *ref;
    if (to_txt_data) 
        ref = BAFTtxt;

    else {
        ASSERTNOTINUSE (vb->scratch);
        buf_alloc_exact (vb, vb->scratch, num_bases, char, "scratch");
        ref = B1STc(vb->scratch);
    }

#ifdef DEBUG
    if (!IS_REF_EXTERNAL) {
        ASSERTNOTINUSE (bitmap_ctx->piz_is_set);
        buf_alloc_exact (vb, bitmap_ctx->piz_is_set, num_bases, uint8_t, "piz_is_set");
    }
#endif

    int32_t idx = RR_IDX ((int32_t)(pos - vb->range->first_pos) - (predict_meth_call ? 2 : 0)); // two bases before pos in case needed for methylation 

    uint32_t num_ref_bases = MIN_(num_bases, range_len - idx);
    ref_get_textual_seq (vb->range->gpos + idx, ref, num_ref_bases, false);

#ifdef DEBUG
    if (!IS_REF_EXTERNAL)
        bits_bit_to_byte (B1ST8(bitmap_ctx->piz_is_set), &vb->range->is_set, idx, num_ref_bases); 
#endif

    // if ref_consumed goes beyond end of range, take the rest from the beginning of range (i.e. circling around)
    if (num_ref_bases < num_bases) {
        ref_get_textual_seq (vb->range->gpos, &ref[num_ref_bases], num_bases - num_ref_bases, false);

#ifdef DEBUG
        if (!IS_REF_EXTERNAL)
            bits_bit_to_byte (B8(bitmap_ctx->piz_is_set, num_ref_bases), &vb->range->is_set, 0, num_bases - num_ref_bases); 
#endif
    }
    
    if (to_txt_data)
        Ltxt += num_bases;
        
    COPY_TIMER (sam_reconstruct_SEQ_get_textual_ref);
    return ref + (predict_meth_call ? 2 : 0); // start of refconsumed
}

// PIZ: SEQ reconstruction 
void sam_reconstruct_SEQ_vs_ref (VBlockP vb_, STRp(snip), ReconType reconstruct)
{
    START_TIMER;

    VBlockSAMP vb = (VBlockSAMP)vb_;
    declare_seq_contexts;

    #define adjusted(base_i) ((base_i) + (predict_meth_call ? 2 : 0))
#ifdef DEBUG // this function is a performance hotspot, therefore this code is only a compilation, not runtime, option
    #define verify_is_set(base_i) ASSPIZ (!flag.debug || IS_REF_EXTERNAL || (adjusted(base_i) < bitmap_ctx->piz_is_set.len32 && *Bc(bitmap_ctx->piz_is_set, adjusted(base_i)) == 1), \
                                          "Expecting POS=%u + base_i=%u to have is_set=1", (PosType32)vb->contexts[SAM_POS].last_value.i, (base_i))
    #define verify_is_set_n(ref, n) ({ for (uint32_t i=0; i < (n); i++) verify_is_set ((ref) + i); })
#else
    #define verify_is_set(base_i)
    #define verify_is_set_n(ref, n)
#endif

    if (!bitmap_ctx->is_loaded) goto done; // if case we need to skip the SEQ field (for the entire VB)

    bool v14            = VER(14);
    char *nonref        = Bc(nonref_ctx->local, nonref_ctx->next_local); // possibly, this VB has no nonref (i.e. everything is ref), in which case nonref would be an invalid pointer. That's ok, as it will not be accessed.
    rom start_nonref    = nonref;
    const PosType32 pos = vb->last_int(SAM_POS);

    // note: prim alignments (loaded in preprocessing) are always mapped - enforced by sam_seg_is_gc_line
    bool unmapped = v14 ? (!vb->preprocessing && (last_flags.unmapped || vb->cigar_missing || !pos || IS_ASTERISK (vb->chrom_name)))
                        : (!pos || IS_ASTERISK (vb->chrom_name)); // the criterion seg used in up to v13
                                                  
    bool aligner_used = v14 ? (snip && snip[SEQ_SNIP_ALIGNER_USER] == '1') // starting v14, seg tells us explicitly
                            : (unmapped && z_file->z_flags.aligner); 

    bool is_perfect = v14 && snip && snip[SEQ_SNIP_PERFECT] == '1';

    bool force_verbatim = v14 && snip && snip[SEQ_FORCE_VERBATIM] == '1';

    char bisulfite = vb->bisulfite_strand; // automatic var for efficiency in tight loop
    bool predict_meth_call = bisulfite && segconf.sam_predict_meth_call;
    
    bool MD_NM_by_unconverted = v14 && bisulfite && segconf.MD_NM_by_unconverted; // we calculate MD:Z and NM:i vs unconverted reference

    bool deep_seq_by_ref = flag.deep && !flag.deep_sam_only && !last_flags.secondary && !last_flags.supplementary;

    if (deep_seq_by_ref && bisulfite) {
        sam_piz_set_deep_seq (vb, false, 0);
        deep_seq_by_ref = false;
    }
    
    if (predict_meth_call && vb->ref_and_seq_consumed) {
        buf_alloc_exact (vb, vb->meth_call, vb->seq_len, char, "meth_call");
        memset (vb->meth_call.data, '.', vb->meth_call.len32);
    }
    
    // case: unmapped, segged against reference using our aligner
    if (aligner_used) {
        aligner_reconstruct_seq (VB, vb->seq_len, false, is_perfect, reconstruct,
                                 MAX_DEEP_SEQ_MISMATCHES, 
                                 deep_seq_by_ref ? vb->deep_mismatch_base   : NULL, 
                                 deep_seq_by_ref ? vb->deep_mismatch_offset : NULL, 
                                 deep_seq_by_ref ? &vb->num_deep_mismatches : NULL);
        
        nonref_ctx->next_local = ROUNDUP4 (nonref_ctx->next_local);
        PosType64 gpos = gpos_ctx->last_value.i; // set in aligner_reconstruct_seq

        if (deep_seq_by_ref) 
            sam_piz_set_deep_seq (vb, true, gpos); 

        return; // aligner accounts for time separately, so we don't account for the time here
    }

    // case: unmapped, and copied verbatim
    if (unmapped || force_verbatim) unmapped: {
        ASSPIZ (nonref_ctx->next_local + vb->seq_len <= nonref_ctx->local.len32, "nonref_ctx.local overflow: next_local=%u seq_len=%u local.len=%u",
                nonref_ctx->next_local, vb->seq_len, nonref_ctx->local.len32);

        uint32_t recon_len = (*nonref == '*') ? 1 : vb->seq_len;
        if (reconstruct) RECONSTRUCT (nonref, recon_len); 
        nonref_ctx->next_local += ROUNDUP4 (recon_len);

        if (deep_seq_by_ref) 
            sam_piz_set_deep_seq (vb, true, NO_GPOS64); 

        goto done;
    }

    // case: missing sequence - sequence is '*' - just reconstruct the '*' (but only if it is the primary SEQ field,
    // which we can tell by not being encountered earlier on the line - bc in v8 (at least) E2 was also stored in the same SAM_SQBITMAP)
    if (vb->seq_missing && !ctx_encountered_in_line (VB, bitmap_ctx->did_i)) {
        if (reconstruct) RECONSTRUCT1 ('*');
        goto done;
    }

    // shortcut in case of a simple cigar (e.g. "151M") with no mismatches or bisulfite
    if (is_perfect && vb->binary_cigar.len32 == 1 && B1ST(BamCigarOp, vb->binary_cigar)->op == BC_M && !bisulfite) {
        sam_reconstruct_SEQ_get_textual_ref (vb, v14, pos, false, true);

        if (deep_seq_by_ref) 
            sam_piz_set_deep_seq (vb, pos + vb->seq_len - 1 <= vb->range->last_pos,  // doesn't round robin at edge of the chromosome
                                  vb->range->gpos + (pos-1)); // gpos
        goto done;
    }

    char *deep_nonref = NULL;
    if (deep_seq_by_ref) {
        buf_alloc_exact (vb, nonref_ctx->deep_nonref, vb->insertions + vb->soft_clip[0] + vb->soft_clip[1], char, "contexts->deep_nonref");
        deep_nonref = B1STc (nonref_ctx->deep_nonref);
    }

    ASSERT0 (is_perfect || bitmap_ctx->local.len32 || nonref_ctx->local.len32, "No SEQ data, perhaps sections were skipped?");

    bitmap_ctx->last_value.i = bitmap_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence

    rom ref = sam_reconstruct_SEQ_get_textual_ref (vb, v14, pos, predict_meth_call, false);
    if (v14 && !ref) goto unmapped; // starting v14, a missing range is another case of unmapped

    // up to v13, alignments with CIGAR=* but with a sequence where segged as nonref
    if (!vb->binary_cigar.len32) {
        buf_alloc_exact (vb, vb->binary_cigar, 1, BamCigarOp, "binary_cigar");
        *B1ST(BamCigarOp, vb->binary_cigar) = (BamCigarOp){ .n = vb->seq_len, .op = BC_S };
    }

    char *recon = BAFTtxt;
    uint32_t seq_consumed=0, ref_consumed=0;

    for_cigar (vb->binary_cigar) {
        case BC_M: case BC_E: case BC_X:
            if (ref) memcpy (recon, &ref[ref_consumed], op->n); // note: up to v13 ref is NULL if all ref_consumed data are mismatches

            if (v14) verify_is_set_n (ref_consumed, op->n); // prior to v14, we only set matching bits (not mismatches)

            if (is_perfect)  //  no mismatches
                recon += op->n;

            else {     
                for (int i=0; i < op->n; i++) {                        
                    // use faster method to find consecutive matches (possible if !bisulfite)
                    if (!bisulfite) {
                        uint32_t matches = bits_get_run ((BitsP)&bitmap_ctx->local, bitmap_ctx->next_local, op->n - i);
                        if (matches) {
                            bool has_mismatch = (i + matches < op->n);
                            recon += matches;
                            bitmap_ctx->next_local += matches + has_mismatch;
                            i += matches - 1 + has_mismatch;

                            if (has_mismatch) goto mismatch; // there is a mismatch after this sequence of matches (i.e. this sequence of matches is not the end of the op)
                        }

                        else {
                            bitmap_ctx->next_local++;
                            goto mismatch;
                        }
                    }

                    // bisulfite: copy converted reference                    
                    else if (NEXTLOCALBIT (bitmap_ctx)) { 
                        // case: bisulfite data, we segged against the converted reference (C->T or G->A) - converted bases represent unmethylated bases (methylated bases remain unconverted)
                        if (bisulfite == ref[ref_consumed + i]) {
                            *recon = (bisulfite == 'C') ? 'T' : 'A'; // convert reference base C->T or G->A
                            
                            if (MD_NM_by_unconverted) {
                                vb->mismatch_bases_by_SEQ++; // this is actually a mismatch vs the unconverted reference
                                bits_clear ((BitsP)&bitmap_ctx->local, bitmap_ctx->next_local-1); // adjust sqbitmap for use for reconstructing MD:Z
                            }

                            if (predict_meth_call)
                                sam_bismark_piz_update_meth_call (vb, ref, ref_consumed+i, seq_consumed+i, i, op->n, op+1, bisulfite, false); // converted therefore unmethylated
                        }

                        recon++;
                    }

                    else mismatch: {
                        char base = v14 ? ({ ContextP mis_ctx = &seqmis_ctx[nuke_encode (ref[ref_consumed+i])] ; // note: a "not enough data" error here is often an indication that pos is wrong
                                                *Bc(mis_ctx->local, mis_ctx->next_local++); })
                                        : *nonref++;

                        *recon++ = base;

                        if (deep_seq_by_ref && vb->mismatch_bases_by_SEQ <= MAX_DEEP_SEQ_MISMATCHES) {
                            vb->deep_mismatch_offset[vb->mismatch_bases_by_SEQ] = seq_consumed + i;
                            vb->deep_mismatch_base[vb->mismatch_bases_by_SEQ] = base;
                        }

                        vb->mismatch_bases_by_SEQ++;

                        if (bisulfite && base == ref[ref_consumed+i]) { // uncoverted base = metyhlated
                            if (segconf.sam_predict_meth_call)
                                sam_bismark_piz_update_meth_call (vb, ref, ref_consumed+i, seq_consumed+i, i, op->n, op+1, bisulfite, true); 

                            if (MD_NM_by_unconverted) {
                                vb->mismatch_bases_by_SEQ--; // this is actually a match vs the unconverted reference, undo the mismatch counter
                                bits_set ((BitsP)&bitmap_ctx->local, bitmap_ctx->next_local-1); // adjust sqbitmap for use for reconstructing MD:Z                        
                            }
                        }
                    }
                }
            }

            seq_consumed += op->n;
            ref_consumed += op->n;
            break;

        case BC_I:
            if (segconf.use_insertion_ctxs) { // skip for now, we will fill in the gaps later
                recon += op->n; 
                seq_consumed += op->n;
                if (deep_seq_by_ref) deep_nonref += op->n; // skip for now
                break;
            }
            // fallthrough

        case BC_S:
            if (op->n == 1) *recon++ = *nonref;  // optimization for common case (helps mostly in pacbio/nanoport that have lots of I)
            else            recon = mempcpy (recon, nonref, op->n);
                
            if (deep_seq_by_ref) {
                if (op->n == 1) *deep_nonref++ = *nonref;
                else            deep_nonref = mempcpy (deep_nonref, nonref, op->n);
            }

            nonref       += op->n;
            seq_consumed += op->n;
            break;

        case BC_D: case BC_N:
            ref_consumed += op->n;
            break;

        default: {} // BC_P, BC_H
    }

    // case: insertions were muxed by the base after - we can fill them in now (since 15.0.30)
    if (segconf.use_insertion_ctxs) {
        SAFE_ASSIGN (recon, 0); // in case last insertion goes to the end of SEQ - it will use this "base" for muxing
        recon -= vb->seq_len;   // go back to the start of the SEQ reconstruction
        deep_nonref = B1STc(nonref_ctx->deep_nonref); // rewind

        for_buf2 (BamCigarOp, op, op_i, vb->binary_cigar) {        
            if (op->op == BC_I) {
                int ins_ctx_i = next_op_is_I (vb, op_i) ? 0 : acgt_encode[(int8_t)recon[op->n]];
                ContextP ins_ctx = seqins_ctx + ins_ctx_i; // mux by base after insertion in SEQ. 
                
                ASSPIZ (ins_ctx->next_local + op->n <= ins_ctx->local.len32, NEXT_ERRFMT " op_i=%u cigar=\"%s\"%s", 
                        __FUNCTION__, ins_ctx->tag_name, ins_ctx->next_local, op->n, ins_ctx->local.len32, op_i, display_binary_cigar (vb),
                        VER2(15,61) ? "" : ". See defect 2024-06-16."); 
       
                recon = mempcpy (recon, Bc(ins_ctx->local, ins_ctx->next_local), op->n);
                
                if (deep_seq_by_ref) 
                    deep_nonref = mempcpy (deep_nonref, Bc(ins_ctx->local, ins_ctx->next_local), op->n);

                ins_ctx->next_local += op->n;
            }

            else if (op->op == BC_M || op->op == BC_S || op->op == BC_E || op->op == BC_X) {
                recon += op->n;

                if (deep_seq_by_ref && op->op == BC_S)
                    deep_nonref += op->n;
            }
        }

        SAFE_RESTORE;
    }
    
    // handle deep seqs passed from SAM reconstruction to FASTQ reconstruction
    if (deep_seq_by_ref) {
        vb->num_deep_mismatches = vb->mismatch_bases_by_SEQ;

        if (nonref_ctx->deep_nonref.len32 && last_flags.rev_comp) 
            str_revcomp_in_place (STRb(nonref_ctx->deep_nonref));

        sam_piz_set_deep_seq (vb, pos + vb->seq_len - 1 <= vb->range->last_pos,  // doesn't round robin at edge of the chromosome
                              vb->range->gpos + (pos-1)); // gpos
    }

    if (reconstruct)
        Ltxt += vb->seq_len;

    ASSPIZ (seq_consumed == vb->seq_len,      "expecting seq_consumed(%u) == vb->seq_len(%u)", seq_consumed, vb->seq_len);
    ASSPIZ (ref_consumed == vb->ref_consumed, "expecting ref_consumed(%u) == vb->ref_consumed(%u)", ref_consumed, vb->ref_consumed);

    nonref_ctx->next_local += ROUNDUP4 (nonref - start_nonref); // I and S data, except if use_insertion_ctxs in which case it is only S

    buf_free (vb->scratch); // allocated by sam_reconstruct_SEQ_get_textual_ref

done:
#ifdef DEBUG
    buf_free (bitmap_ctx->piz_is_set); 
#endif
    COPY_TIMER (sam_reconstruct_SEQ_vs_ref);
}


// reconstruct from a 2bit array - start_base and seq_len are in bases (not in bits)
static void reconstruct_SEQ_acgt (VBlockSAMP vb, BitsP seq_2bits, uint64_t start_base, uint32_t seq_len, bool revcomp)
{
    // Reconstruct as SAM. If BAM, we later translate to BAM. To do: reconstruct directly as BAM if needed, bug 530
    char *next = BAFTtxt;

    if (!revcomp)
        for (uint64_t i=0; i < seq_len; i++) {
            uint8_t b = bits_get2 (seq_2bits, (start_base + i)*2);
            *next++ = b==0 ? 'A' : b==1 ? 'C' : b==2 ? 'G' : 'T';
        }
    
    else
        for (int64_t i=seq_len-1; i >= 0; i--) {
            uint8_t b = bits_get2 (seq_2bits, (start_base + i)*2);
            *next++ = b==3 ? 'A' : b==2 ? 'C' : b==1 ? 'G' : 'T';
        }

    Ltxt += seq_len;
}

// PRIM or DEPN VB
static void reconstruct_SEQ_copy_sag_prim (VBlockSAMP vb, ContextP ctx, 
                                           BitsP prim, uint64_t prim_start_base, uint32_t prim_seq_len, 
                                           bool xstrand, ReconType reconstruct, bool force_no_analyze_depn_SEQ)
{
    START_TIMER;

    if (!reconstruct || !ctx->is_loaded) return;

    rom seq = BAFTtxt;

    uint32_t recon_seq_len = prim_seq_len;
    uint64_t recon_start_base = prim_start_base;

    if (IS_DEPN(vb)) {
        recon_seq_len -= vb->hard_clip[0] + vb->hard_clip[1]; // depn sequence is a sub-sequence of the prim sequence
        recon_start_base += xstrand ? (prim_seq_len - recon_seq_len - vb->hard_clip[0]) : vb->hard_clip[0];
    }

    reconstruct_SEQ_acgt (vb, prim, recon_start_base, recon_seq_len, xstrand);

    bool has_deep = (flag.deep && !flag.deep_sam_only && !last_flags.secondary && !last_flags.supplementary && !vb->bisulfite_strand);
    
    // if needed, get vb->mismatch_bases_by_SEQ and line_sqbitmap for use in reconstructing MD:Z and NM:i, as well as Deep
    // NOTE: if SEQ is not loaded or !reconstruct - it is expected that NM and MD also don't reconstruct
    if ((segconf.has_MD_or_NM || has_deep) && !vb->cigar_missing && !vb->seq_missing && !force_no_analyze_depn_SEQ) {
        PosType64 range_last_pos=0, range_gpos=0;
        PosType32 pos = vb->last_int(SAM_POS);

        sam_analyze_copied_SEQ (vb, seq, recon_seq_len, pos, last_flags.rev_comp,
                                vb->ref_consumed, vb->ref_and_seq_consumed, &ctx->line_sqbitmap, 
                                has_deep, &range_last_pos, &range_gpos);        

        // case: reconstructing a PRIM line - calculate vb->deep_seq if it has at most one mismatch vs reference
        if (has_deep) 
            sam_piz_set_deep_seq (vb, pos + recon_seq_len - 1 <= range_last_pos,  // doesn't round robin at edge of the chromosome
                                  range_gpos + (pos-1));
    }

    COPY_TIMER (reconstruct_SEQ_copy_sag_prim);
}

static void reconstruct_SEQ_copy_saggy (VBlockSAMP vb, ContextP ctx, ReconType reconstruct, bool force_no_analyze_depn_SEQ)
{
    START_TIMER;
    if (!reconstruct || !ctx->is_loaded) return;

    SamFlags saggy_flags = (SamFlags){ .value = *B(int64_t, CTX(SAM_FLAG)->history, vb->saggy_line_i) };
    bool xstrand = saggy_flags.rev_comp != last_flags.rev_comp;

    HistoryWord word = *B(HistoryWord, ctx->history, vb->saggy_line_i); // SEQ is always stored as LookupTxtData or LookupPerLine
    uint32_t seq_len = word.len - vb->hard_clip[0] - vb->hard_clip[1]; // our sequence is a sub-sequence of the saggh sequence

    char *seq = BAFTtxt;

    rom saggy_seq = ((word.lookup == LookupTxtData) ? Btxt (word.index) : Bc(ctx->per_line, word.index));

    if (OUT_DT(SAM)) {
        if (xstrand)
            str_revcomp (seq, saggy_seq + vb->hard_clip[1], seq_len);
        else
            memcpy (seq, saggy_seq + vb->hard_clip[0], seq_len);
    
        Ltxt += seq_len;
    }

    else { // bam or cram
        uint32_t clip = vb->hard_clip[xstrand];
        bam_seq_to_sam (VB, (uint8_t*)saggy_seq + clip/2, seq_len, clip%2, false, &vb->txt_data, false);

        if (xstrand)
            str_revcomp_in_place (STRa(seq));
    }
    
    bool has_deep = flag.deep && !flag.deep_sam_only && !last_flags.secondary && !last_flags.supplementary && !vb->bisulfite_strand;

    // in needed, get vb->mismatch_bases_by_SEQ and line_sqbitmap for use in reconstructing MD:Z and NM:i
    // NOTE: if SEQ is not loaded or !reconstruct - it is expected that NM and MD also don't reconstruct
    if ((segconf.has_MD_or_NM && !vb->cigar_missing && !vb->seq_missing && !force_no_analyze_depn_SEQ) || 
        has_deep) {
        
        PosType64 range_last_pos=0, range_gpos=0;
        PosType32 pos = vb->last_int(SAM_POS);

        sam_analyze_copied_SEQ (vb, STRa(seq), pos, last_flags.rev_comp,
                                vb->ref_consumed, vb->ref_and_seq_consumed, &ctx->line_sqbitmap, 
                                has_deep, &range_last_pos, &range_gpos);        

        if (has_deep)
            sam_piz_set_deep_seq (vb, pos + seq_len - 1 <= range_last_pos, // doesn't round robin at edge of the chromosome  
                                  range_gpos + (pos-1)); // gpos 
    }

    COPY_TIMER (reconstruct_SEQ_copy_saggy);
}

// PIZ of a SEQ in a MAIN and DEPN components - this is an all-the-same context - SEQ all lines in a DEPN component
// come here. We multiplex between SAGroup-based reconstruction vs normal SEQ reconstruction
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_SEQ)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    START_TIMER;

    if (!CTX(SAM_SQBITMAP)->flags.no_textual_seq)
        reconstruct = true; // SEQ should always be reconstructed because it is consumed by other fields

    vb->bisulfite_strand = (!VER(14) || snip[SEQ_SNIP_BISULFATE]=='0') ? 0
                         : snip[SEQ_SNIP_BISULFATE] == '*'             ? "CG"[last_flags.rev_comp]
                         :                                               snip[SEQ_SNIP_BISULFATE];
    // case: reconstructor would normally diff, but in this case it should reconstruct from SQBITMAP.local    
    if (snip[SEQ_SNIP_FORCE_SEQ] == '1') 
        goto vs_seq__aligner__verbatim;

    // case: PRIM component - copy from SA Group 
    else if (IS_PRIM(vb) && !vb->preprocessing) 
        reconstruct_SEQ_copy_sag_prim (vb, ctx, (BitsP)&z_file->sag_seq, STRa(vb->sag->seq), false, reconstruct, false);

    else if (snip[SEQ_SNIP_ALIGNER_USER] == '1' || snip[SEQ_FORCE_VERBATIM] == '1')
        goto vs_seq__aligner__verbatim;

    // case: DEPN component - line with sag - copy from sag loaded from prim VB
    else if (IS_DEPN(vb) && SAM_PIZ_HAS_SAG) {
        const Sag *grp = vb->sag;
        ASSPIZ (grp->seq_len > vb->hard_clip[0] + vb->hard_clip[1], "grp_i=%u aln_i=%"PRIu64" grp->seq_len=%u <= hard_clip[LEFT]=%u + hard_clip[RIGHT]=%u", 
                ZGRP_I(grp), ZALN_I(vb->sa_aln), grp->seq_len, vb->hard_clip[0], vb->hard_clip[1]);

        reconstruct_SEQ_copy_sag_prim (vb, ctx, (BitsP)&z_file->sag_seq, STRa(grp->seq), 
                                       (last_flags.rev_comp != grp->revcomp), reconstruct, snip[SEQ_SNIP_NO_ANALYZE_DEPN_SEQ]-'0');
    }

    // case: copy from saggy line in this VB
    else if (sam_has_saggy && !vb->seq_missing &&
             !B(CigarAnalItem, CTX(SAM_CIGAR)->cigar_anal_history, vb->saggy_line_i)->hard_clip[0] &&
             !B(CigarAnalItem, CTX(SAM_CIGAR)->cigar_anal_history, vb->saggy_line_i)->hard_clip[1])
        reconstruct_SEQ_copy_saggy (vb, ctx, reconstruct, snip[SEQ_SNIP_NO_ANALYZE_DEPN_SEQ]-'0');

    // case: reconstruct sequence directly
    else 
        vs_seq__aligner__verbatim:
        sam_reconstruct_SEQ_vs_ref (VB, STRa(snip), reconstruct);

    COPY_TIMER (sam_piz_special_SEQ);
    return NO_NEW_VALUE;
}

// SAM-to-BAM translator: translate SAM ASCII sequence characters to BAM's 4-bit characters:
TRANSLATOR_FUNC (sam_piz_sam2bam_SEQ)
{
    START_TIMER;

    declare_seq_contexts;

    // before translating - add to Deep if needed
    if (flag.deep) 
        sam_piz_deep_SEQ_cb (VB_SAM, STRa(recon));
    
    BAMAlignmentFixedP alignment = (BAMAlignmentFixedP)Btxt (vb->line_start);
    uint32_t l_seq = alignment->l_seq;

    // backward compatability note: prior to v14 the condition was:
    // if (CTX(SAM_QUAL)->lcodec == CODEC_LONGR) ...
    // Since v14, it is determined by a flag. Since this flag is 0 in V<=13, earlier files will always

    // case downstream contexts need access to the textual_seq: copy it.
    // Examples: LONGR, HOMP, t0 codecs; sam_piz_special_BSSEEKER2_XM ; sam_piz_special_ULTIMA_tp
    if (!bitmap_ctx->flags.no_textual_seq) {
        buf_alloc (vb, &VB_SAM->textual_seq, 0, recon_len, char, 0, "textual_seq");
        memcpy (B1STc(VB_SAM->textual_seq), recon, recon_len);
        VB_SAM->textual_seq.len = recon_len;
    }

    if (vb->drop_curr_line) goto done; // sequence was not reconstructed - nothing to translate
    
    // if l_seq=0, just remove the '*'
    if (!l_seq) {
        Ltxt--;
        goto done;
    }

    // if l_seq is odd, 0 the next byte that will be half of our last result byte
    if (l_seq % 2) *BAFTtxt = 0; 

    // the characters "=ACMGRSVTWYHKDBN" are mapped to BAM 0->15, in this matrix we add 0x80 as a validity bit. All other characters are 0x00 - invalid
    static const uint8_t sam2bam_seq_map[256] = { ['=']=0x0, ['A']=0x1, ['C']=0x2, ['M']=0x3, ['G']=0x4, ['R']=0x5, ['S']=0x6, ['V']=0x7, 
                                                  ['T']=0x8, ['W']=0x9, ['Y']=0xa, ['H']=0xb, ['K']=0xc, ['D']=0xd, ['B']=0xe, ['N']=0xf };

    uint8_t *in=(uint8_t *)recon, *out=in; 
    uint32_t out_len = (l_seq + 1) / 2;
    
    for (uint32_t i=0; i < out_len; i++, out++, in += 2) 
        *out = (sam2bam_seq_map[in[0]] << 4) | sam2bam_seq_map[in[1]]; // note: invalid characters are encoded as 0 (piz doesn't verify)
    
    Ltxt = Ltxt - l_seq + out_len;

    done: 
    COPY_TIMER (sam_piz_sam2bam_SEQ);
    return 0;
}

// PIZ
rom sam_piz_get_textual_seq (VBlockP vb_)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    if (!vb->textual_seq_str) { // first call
        vb->textual_seq_str = OUT_DT(SAM)           ? last_txt (VB, SAM_SQBITMAP)
                            : vb->textual_seq.len32 ? B1STc(vb->textual_seq) // generated by sam_piz_sam2bam_SEQ
                            :                         NULL;
        
        ASSPIZ (vb->textual_seq_str, "textual_seq_str=NULL (ctx->flags.no_textual_seq=%s)", TF(CTX(SAM_SQBITMAP)->flags.no_textual_seq));  
    }

    return vb->textual_seq_str; // note: textual_seq_str is prepared in sam_piz_sam2bam_SEQ sam_load_groups_add_seq
}

