// ------------------------------------------------------------------
//   sam_md.c
//   Copyright (C) 2021-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// ----------------------------------------------------------------------------------------
// MD:Z - "Mismatch & Deleted bases" - see https://samtools.github.io/hts-specs/SAMtags.pdf
// ----------------------------------------------------------------------------------------

#include "sam_private.h"

//---------
// SEG
//---------

// called when segging SEQ: if MD is already verified as consistent with CIGAR,
// we now verify that the mismatches are consistent with SEQ.
void sam_MD_Z_verify_due_to_seq (VBlockSAMP vb, STRp(seq), PosType32 pos, BitsP sqbitmap, uint64_t sqbitmap_start)
{
    BitsP M_is_ref = (BitsP)&vb->md_M_is_ref;

    bool bitmap_matches_MD = vb->md_verified && !bits_hamming_distance (M_is_ref, 0, sqbitmap, sqbitmap_start, M_is_ref->nbits);

    if (flag.show_wrong_md && vb->md_verified && !bitmap_matches_MD) {

        iprintf ("%s RNAME=%.*s POS=%d CIGAR=%s MD=%.*s SEQ=%.*s\n", 
                LN_NAME, STRf(vb->chrom_name), pos, vb->last_cigar, vb->last_txt_len(OPTION_MD_Z), last_txt(VB, OPTION_MD_Z), STRf(seq));
        if (sqbitmap) 
            bits_print_substr ("SEQ match to ref", sqbitmap, sqbitmap_start, M_is_ref->nbits, info_stream);
        else
            iprint0 ("SEQ: no bitmap\n");
        bits_print_substr ("MD implied match", M_is_ref, 0, M_is_ref->nbits, info_stream); 
    }

    vb->md_verified = bitmap_matches_MD;
}

static inline rom sam_md_consume_D (VBlockSAMP vb, bool is_depn, char **md_in_out, uint32_t *M_D_bases, PosType32 *pos, int D_bases, 
                                    RangeP *range_p, RefLock *lock, bool *critical_error)
{
    char *md = *md_in_out;

    if (! *md || *md != '^' || !IS_ACGT(md[1])) {
        *critical_error = true;
        return "Malformed MD while parsing D - expecting '^'"; // expecting a deletion (must include at least one base).
    }

    if (!IS_DIGIT (md[-1])) {
        *critical_error = true;
        return "Malformed MD while parsing D - expecting preceding character to be a digit"; 
    }
    
    md++;

    rom error=NULL;
    while (IS_ACGT(*md) && D_bases) {
        if (!error)
            error = sam_seg_analyze_set_one_ref_base (vb, is_depn, *pos, *md, *M_D_bases, range_p, lock); 
            
        D_bases--;
        (*M_D_bases)--;
        (*pos)++;
        md++;
    }

    if (IS_ACGT(*md) || D_bases) {
        *critical_error = true;
        return "CIGAR D length and MD number of deleted bases don't match"; 
    }

    *md_in_out = md;
    return error; // NULL if success
}

// verifies that the reference matches as required, and updates reference bases if missing
static inline rom sam_md_consume_M (VBlockSAMP vb, bool is_depn, char **md_in_out, uint32_t *M_D_bases, PosType32 *pos, int M_bases,
                                    Bits *M_is_ref, uint64_t *M_is_ref_i,
                                    RangeP *range_p, RefLock *lock, bool *critical_error)
{
    char *md = *md_in_out;
    rom error = NULL;

    // expecting a series of <number><base> where number can be 0 and the base of the last pair can be missing. eg: 0T12A4
    while (M_bases) {
        if (!IS_DIGIT(*md)) {
            *critical_error = true; 
            return "Malformed MD while parsing M - expecting digit"; 
        }

        // matching bases
        int match_len = strtod (md, &md); // get number and advance past number

        // case: MD number is bigger than needed by current CIGAR op (perhaps partially covering the next CIGAR op) - update MD in-place
        if (match_len > M_bases || (match_len == M_bases && *md && *md != '^')) {
            int int_len = str_int_ex (match_len - M_bases, *md_in_out, false);
            memmove (*md_in_out + int_len, md, strlen(md) + 1/*\0*/); // note: memmove and not sprintf, bc moving to overlapping memory
            md = *md_in_out;
            match_len = M_bases;
        }

        M_bases     -= match_len;
        *M_D_bases  -= match_len;
        *pos        += match_len;
        *M_is_ref_i += match_len;              
        
        ASSERT (M_bases >= 0, "Expecting M_bases=%d >= 0", M_bases);

        // if we still need more M_bases, the next one should be a mismatch nucleotide. Note that the SAM standard permits IUPAC
        // "bases" (eg N), but we apply the MD special alg only for ACGT.
        if (M_bases) {
            if (*md == '\0') {
                *critical_error = true; 
                return "MD string implies a shorter SEQ than CIGAR does"; 
            }

            else if (!IS_ACGT (*md))
                error = (*md=='N' ? "Encountered 'N' base while parsing M" : "Not A,C,G,T,N while parsing M"); // Genozip reference supports only A,C,G,T, but this "base" in the MD string is not one of them

            else { // set base (if A,C,G,T) even if previous bases had an error
                rom result = sam_seg_analyze_set_one_ref_base (vb, is_depn, *pos, *md, *M_D_bases, range_p, lock); // continue counting mismatch_bases_by_MD despite error
                if (result && !error) error = result;
            }

            bits_clear (M_is_ref, *M_is_ref_i); // base in SEQ is expected to be NOT equal to the reference base
            (*M_is_ref_i)++;              

            M_bases--;
            (*M_D_bases)--;
            (*pos)++;
            md++;
            vb->mismatch_bases_by_MD++;
        }
    }

    *md_in_out = md;
    return error; // NULL if success
}

// called for analyzing MD before segging SEQ and later MD
// - Verifies that MD is consistent with CIGAR
// - Verifies that the mismatched and deleted bases are the same as the reference, or if not 
//   in the reference yet (in REF_INTERNAL), adds them
// - Sets md_M_is_ref, of length ref_and_seq_consumed (corresponding to M, X and = CIGAR ops): 1 for a base matching the reference, 0 for not.
//   sam_seg_SEQ will conduct the final verification step of comparing this bitmap to the one calculated from the SEQ data.
// Note: a correct MD:Z may still appear sometimes as "unverified" in REF_INTERNAL, if the locus in the internal reference was populated
// by an earlier read with a base different than the base in the reference with which this SAM file was generated.
void sam_seg_MD_Z_analyze (VBlockSAMP vb, ZipDataLineSAM *dl, rom md_orig, uint32_t md_len, PosType32 pos)
{
    START_TIMER;

    if (segconf_running || dl->FLAG.unmapped) return; // sometimes unmapped reads have a bogus CIGAR and MD:Z, however, if we try to verify sam_seg_SEQ_vs_ref (which we no longer do), verification fails.

    rom reason=NULL;
    #define not_verified(s) { reason=s ; goto not_verified; }

    RangeP range = NULL;
    RefLock lock = REFLOCK_NONE;
    
    bool is_depn = (IS_DEPN(vb) && vb->sag) || sam_has_saggy; 

    if (flag.show_wrong_md)
        seg_set_last_txt (VB, CTX(OPTION_MD_Z), md_orig, md_len); // consumed in sam_seg_SEQ

    // copy of MD as we are going to modify it (but we still need the original intact for sam_seg_MD_Z)
    char md_data[md_len+2];
    char *md = memcpy (md_data+1, md_orig, md_len);
    md[-1] = md[md_len] = 0; // nul-terminator on both ends (so we can test preceding character with md[-1])

    if (!pos) not_verified ("No POS")
    if (IS_ASTERISK (vb->chrom_name)) not_verified ("No RNAME");
    if (vb->cigar_missing) not_verified ("No CIGAR");

    // According to the specification (https://samtools.github.io/hts-specs/SAMtags.pdf), an MD string may start or end with a mismatch or D sequence.
    // however, in actual BAM files in the wild MD always starts and ends with a digit (possibly 0). 
    // Therefore, Genozip only activates the special MD alg if its starts and ends with digit.
    if (!md_len || !IS_DIGIT(md[0]) || !IS_DIGIT(md[md_len-1])) not_verified ("Malformed MD:Z field");

    // start by marking all as matching, and clear the SNPs later
    buf_alloc_bits_exact (vb, &vb->md_M_is_ref, vb->ref_and_seq_consumed, SET, CTX_GROWTH, "md_M_is_ref"); 
    BitsP M_is_ref = (BitsP)&vb->md_M_is_ref;
    
    uint64_t M_is_ref_i=0;
    
    uint32_t M_D_bases = vb->ref_consumed; // M/=/X and D
    
    bool critical_error=false;
    rom error=NULL;
    for_cigar (vb->binary_cigar) {
        case BC_M: case BC_E: case BC_X:
            if ((error = sam_md_consume_M (vb, is_depn, &md, &M_D_bases, &pos, op->n, M_is_ref, &M_is_ref_i, &range, &lock, &critical_error))
                && critical_error) // break loop now if critical error, else continue to count mismatch_bases_by_MD despite error
                not_verified (error);
            if (!reason) reason = error; // non-critical error - continue
            break;

        case BC_D: 
            if ((error = sam_md_consume_D (vb, is_depn, &md, &M_D_bases, &pos, op->n, &range, &lock, &critical_error)) 
                && critical_error)
                not_verified (error);
            if (!reason) reason = error; // non-critical error - continue
            break;

        case BC_N: 
            pos += op->n;
            M_D_bases -= op->n;
            break;

        default: {}
    }

    // case: we didn't consume the entire MD (except an allowed trailing '0')
    if (*md && !(*md=='0' && !md[1])) not_verified ("Failed to consume entire MD"); 

    if (reason) goto not_verified_buf_mismatch_count_ok; // now we can handle the non-critical error raised in sam_md_consume_M

    if (range) ref_unlock (&lock);

    vb->md_verified = true;    
    goto done; // verified

not_verified:
    vb->mismatch_bases_by_MD = -1; // fallthrough

not_verified_buf_mismatch_count_ok:    
    if (range) ref_unlock (&lock);

    vb->md_verified = false;

    if (flag.show_wrong_md)
        iprintf ("%s\tRNAME=%.*s\tPOS=%d\tCIGAR=%s\tMD=%.*s\tReason=%s\n", 
                 LN_NAME, STRf(vb->chrom_name), pos, vb->last_cigar, vb->last_txt_len(OPTION_MD_Z), last_txt(VB, OPTION_MD_Z), reason);
    // fallthrough
    
done:
    COPY_TIMER (sam_seg_MD_Z_analyze);
}

// MD's logical length is normally the same as seq_len, we use this to optimize it.
// In the common case that it is just a number equal the seq_len, we replace it with an empty string.
// if MD value can be derived from the seq_len, we don't need to store - store just an empty string
void sam_seg_MD_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(md), unsigned add_bytes)
{
    decl_ctx (OPTION_MD_Z);
    uint32_t md_value;

    ctx_set_encountered (VB, CTX(OPTION_MD_Z));
    
    if (vb->md_verified && !vb->cigar_missing && !vb->seq_missing) // note: md_verified might have been reset in sam_seg_SEQ 
        seg_special0 (VB, SAM_SPECIAL_MD, ctx, add_bytes);

    // in case not verified (eg bc it is depn) but (common case) it is equal to (seq_len - soft_clips)
    else if (!vb->md_verified && !vb->seq_missing && str_get_uint32 (STRa(md), &md_value) && 
             md_value == dl->SEQ.len - vb->soft_clip[0] - vb->soft_clip[1]) 
        seg_special4 (VB, SAM_SPECIAL_SEQ_LEN, '+', '0', '-', '-', ctx, add_bytes); // also used in sam_seg_SEQ_END

    // store in local. note: this is tested to be much better than splitting by length and storing short ones in dict
    else  
        seg_add_to_local_string (VB, ctx, STRa(md), LOOKUP_SIMPLE, add_bytes);
}

//---------
// PIZ
//---------

SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_MD)
{
    START_TIMER;

    VBlockSAMP vb = (VBlockSAMP)vb_;
    bool perfect_match = VER(14) && !vb->mismatch_bases_by_SEQ; // "perfect" encoding introduced v14

    PosType32 pos = LOADED_CTX(SAM_POS)->last_value.i;

    ASSERTNOTEMPTY (vb->binary_cigar);

    ContextP sqbitmap_ctx = LOADED_CTX(SAM_SQBITMAP);
    
    BitsP line_sqbitmap = (BitsP)&sqbitmap_ctx->line_sqbitmap;
    bool is_depn = !!line_sqbitmap->nwords;
    uint32_t line_sqbitmap_next = 0;

    uint32_t save_next_local = sqbitmap_ctx->next_local;
    
    if (!perfect_match && !is_depn)
        sqbitmap_ctx->next_local = sqbitmap_ctx->last_value.i; // rewind back to beginning of the bits of this line (value stored by sam_reconstruct_SEQ_vs_ref)

    uint32_t count_match=0;
    decl_acgt_decode;
    
    for_cigar (vb->binary_cigar) {
        case BC_M : case BC_E : case BC_X : {
            uint32_t n = op->n;      

            if (perfect_match) {
                count_match += n;
                pos += n;
            }
            
            // reconstruct a series of <number><base> where number can be 0 and the base of the last pair can be missing. eg: 0T12A4
            else while (n) {
                uint32_t run = is_depn ? bits_get_run (line_sqbitmap, line_sqbitmap_next, n)
                                       : bits_get_run ((BitsP)&sqbitmap_ctx->local, sqbitmap_ctx->next_local, n);
                count_match += run;
                n           -= run;
                pos         += run;

                if (is_depn) line_sqbitmap_next       += run + (n > 0);
                else         sqbitmap_ctx->next_local += run + (n > 0);

                if (n) {
                    // usually vb->range set in sam_reconstruct_SEQ_vs_ref, except if no base in SEQ matched the reference. In that case, we set it here
                    if (!vb->range) 
                        // note: if this throws an error, it is possibly related to private/defects/SAM-MD-external-ref-contig-with-no-ref.txt
                        vb->range = (RangeP)ref_piz_get_range (VB, HARD_FAIL);

                    if (reconstruct) RECONSTRUCT_INT (count_match); // flush matches before reconstructing mismatch
                    count_match=0;
                    
                    if (reconstruct) RECONSTRUCT1 (ref_base_by_pos (vb->range, pos)); 
                    n--;
                    pos++;
                }
            }
            break;
        }

        case BC_D :
            // flush matches before reconstructing deletion (but not if deletion is first)
            if (BNUM(vb->binary_cigar, op)) { // not first op
                if (reconstruct) RECONSTRUCT_INT (count_match);
                count_match=0;
            }

            if (reconstruct) RECONSTRUCT1 ('^');

            if (!vb->range) 
                vb->range = (RangeP)ref_piz_get_range (VB, HARD_FAIL);

            for (uint32_t i=0; i < op->n; i++) { 
                if (reconstruct) RECONSTRUCT1 (ref_base_by_pos (vb->range, pos)); 
                pos++;
            }
            break;

        case BC_N :  // skipping without deletion in MD
            pos += op->n;
            break;

        default : {}
    }

    if (count_match || !IS_DIGIT (*BLSTtxt)) 
        if (reconstruct) RECONSTRUCT_INT (count_match); // flush matches if any unflushed yet

    ASSPIZ (perfect_match || save_next_local == sqbitmap_ctx->next_local, "expecting save_next_local=%u == sqbitmap_ctx->next_local=%u (RNAME=%s POS=%"PRId64" CIGAR=%s)", 
            save_next_local, sqbitmap_ctx->next_local, vb->chrom_name, CTX(SAM_POS)->last_value.i, VB_SAM->textual_cigar.data);

    COPY_TIMER (sam_piz_special_MD);
    return NO_NEW_VALUE;
}

// Used in files compressed with Genozip up to 12.0.36.
// logic: snip is eg "119C" (possibly also "") - we reconstruct the original, eg "119C31" 
// by concating a number which is (seq_len - partial_seq_len_by_md_field)
SPECIAL_RECONSTRUCTOR (sam_piz_special_MD_old)
{
    if (!reconstruct) return false;
    
    if (snip_len) RECONSTRUCT_snip;

    unsigned partial_seq_len_by_md_field=0, curr_num=0;

    for (unsigned i=0; i < snip_len; i++) {   
        if (IS_DIGIT (snip[i])) 
            curr_num = curr_num * 10 + (snip[i] - '0');

        else {
            partial_seq_len_by_md_field += curr_num + 1; // number terminates here + one character
            curr_num = 0;
        }
    }

    partial_seq_len_by_md_field += curr_num; // in case the string ends with a number

    RECONSTRUCT_INT (vb->seq_len - partial_seq_len_by_md_field);

    return NO_NEW_VALUE;
}
