// ------------------------------------------------------------------
//   sam_md.c
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

// ----------------------------------------------------------------------------------------
// MD:Z - "Mismatch & Deleted bases" - see https://samtools.github.io/hts-specs/SAMtags.pdf
// ----------------------------------------------------------------------------------------

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "reference.h"
#include "segconf.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"

//---------
// SEG
//---------

// called when segging SEQ
void sam_MD_Z_verify_due_to_seq (VBlockSAMP vb, STRp(seq), SamPosType pos, BitArrayP sqbitmap, uint64_t sqbitmap_start)
{
    BitArrayP M_is_ref = (BitArrayP)&vb->md_M_is_ref;

    bool bitmap_matches_MD = vb->md_verified && !bit_array_hamming_distance (M_is_ref, 0, sqbitmap, sqbitmap_start, M_is_ref->nbits);

    if (flag.show_wrong_md && vb->md_verified && !bitmap_matches_MD) {

        iprintf ("%s RNAME=%.*s POS=%d CIGAR=%s MD=%.*s SEQ=%.*s\n", 
                LN_NAME, STRf(vb->chrom_name), pos, vb->last_cigar, vb->last_txt_len(OPTION_MD_Z), last_txt(VB, OPTION_MD_Z), STRf(seq));
        if (sqbitmap) 
            bit_array_print_substr ("SEQ match to ref", sqbitmap, sqbitmap_start, M_is_ref->nbits, info_stream);
        else
            iprint0 ("SEQ: no bitmap\n");
        bit_array_print_substr ("MD implied match", M_is_ref, 0, M_is_ref->nbits, info_stream); 
    }

    vb->md_verified = bitmap_matches_MD;
}

static rom sam_md_set_one_ref_base (VBlockSAMP vb, bool is_depn, SamPosType pos, char base, uint32_t M_D_bases,                                      
                                    RangeP *range_p, RefLock *lock)
{
    // case: pos is beyond the existing range
    if ((*range_p) && (*range_p)->last_pos < pos && !is_depn) {
        ref_unlock (gref, *lock);
        *range_p = NULL;
    }

    // get range (and lock it if needed)
    if (! *range_p) {
        *range_p = ref_seg_get_locked_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, M_D_bases, WORD_INDEX_NONE, NULL, 
                                             (flag.reference == REF_EXTERNAL || is_depn) ? NULL : lock);
        if (! *range_p) return "Range not available"; // cannot access this range in the reference
    }

    uint32_t pos_index = pos - (*range_p)->first_pos; // index within range

    // case: depn line, but we don't haven't set the reference data for it. since we don't need the reference data for reconstructing
    // SEQ (it is copied from prim), we don't store it just for MD:Z - so we won't use SPECIAL for this MD:Z
    if (is_depn && (flag.reference & REF_STORED) && !ref_is_nucleotide_set (*range_p, pos_index))
        return "Depn alignment base not is_set in reference";

    bool internal_pos_is_populated = (flag.reference == REF_INTERNAL) && ref_is_nucleotide_set (*range_p, pos_index);

    // case: reference already contains a base - but unfortunately it is not "base" - so MD is not reconstractable from reference
    if (((flag.reference & REF_ZIP_LOADED) || internal_pos_is_populated) && (base != ref_base_by_pos (*range_p, pos))) 
        return "MD has incorrect reference base"; // encountered in the wild when the reference base is a IUPAC
    
    // case: reference is not set yet - set it now
    if (flag.reference == REF_INTERNAL && !internal_pos_is_populated)
        ref_set_nucleotide (*range_p, pos_index, base);

    // set is_set - we will need this base in the reference to reconstruct MD
    if (flag.reference & REF_STORED)
        bit_array_set (&(*range_p)->is_set, pos_index); // we will need this ref to reconstruct

    return NULL; // success
}

static inline rom sam_md_consume_D (VBlockSAMP vb, bool is_depn, char **md_in_out, uint32_t *M_D_bases, SamPosType *pos, int D_bases, 
                                    RangeP *range_p, RefLock *lock, bool *critical_error)
{
    char *md = *md_in_out;

    if (! *md || *md != '^' || !IS_NUCLEOTIDE(md[1])) {
        *critical_error = true;
        return "Malformed MD while parsing D - expecting '^'"; // expecting a deletion (must include at least one base)
    }

    md++;

    rom error=NULL;
    while (IS_NUCLEOTIDE(*md) && D_bases) {
        if (!error)
            error = sam_md_set_one_ref_base (vb, is_depn, *pos, *md, *M_D_bases, range_p, lock); // // continue counting mismatch_bases_by_MD despite error
            
        D_bases--;
        (*M_D_bases)--;
        (*pos)++;
        md++;
    }

    if (IS_NUCLEOTIDE(*md) || D_bases) {
        *critical_error = true;
        return "CIGAR D length and MD number of deleted bases don't match"; 
    }

    *md_in_out = md;
    return error; // NULL if success
}

// verifies that the reference matches as required, and updates reference bases if missing
static inline rom sam_md_consume_M (VBlockSAMP vb, bool is_depn, char **md_in_out, uint32_t *M_D_bases, SamPosType *pos, int M_bases,
                                    BitArray *M_is_ref, uint64_t *M_is_ref_i,
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
            sprintf (*md_in_out, "%u%s", match_len - M_bases, md);
            md = *md_in_out;
            match_len = M_bases;
        }

        M_bases     -= match_len;
        *M_D_bases  -= match_len;
        *pos        += match_len;
        *M_is_ref_i += match_len;              
        
        // if we still need more M_bases, the next one should be a mismatch nucleotide. Note that the SAM standard permits IUPAC
        // "bases" (eg N), but we apply the MD special alg only for ACGT.
        if (M_bases) {

            if (!IS_NUCLEOTIDE (*md))
                error = (*md=='N' ? "Encountered 'N' base while parsing M" : "Not A,C,G,T,N while parsing M"); // Genozip reference supports only A,C,G,T, but this "base" in the MD string is not one of them

            if (!error) 
                error = sam_md_set_one_ref_base (vb, is_depn, *pos, *md, *M_D_bases, range_p, lock); // continue counting mismatch_bases_by_MD despite error

            bit_array_clear (M_is_ref, *M_is_ref_i); // base in SEQ is expected to be NOT equal to the reference base
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
void sam_seg_MD_Z_analyze (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(md), SamPosType pos)
{
    if (dl->FLAG.bits.unmapped) return; // sometimes unmapped reads have a bogus CIGAR and MD:Z, however, if we try to verify sam_seg_SEQ_vs_ref (which we no longer do), verification fails.

    rom reason=NULL;
    #define not_verified(s) { reason=s ; goto not_verified; }

    RangeP range = NULL;
    RefLock lock = REFLOCK_NONE;
    
    bool is_depn = (sam_is_depn_vb && vb->sa_grp) || zip_has_prim;

    if (flag.show_wrong_md)
        seg_set_last_txt (VB, CTX(OPTION_MD_Z), STRa(md)); // consumed in sam_seg_SEQ

    // copy of MD as we are going to modify it (but we still need the original intact for sam_MD_Z_seg)
    char md_data[md_len+1];
    memcpy (md_data, md, md_len);
    md_data[md_len] = 0;
    md = md_data;

    if (!pos) not_verified ("No POS")
    if (vb->chrom_name_len==1 && vb->chrom_name[0]=='*') not_verified ("No RNAME");
    if (vb->cigar_missing) not_verified ("No CIGAR");

    // According to the specification (https://samtools.github.io/hts-specs/SAMtags.pdf), an MD string may start or end with a mismatch or D sequence.
    // however, in actual BAM files in the wild MD always starts and ends with a digit (possibly 0). 
    // Therefore, Genozip only activates the special MD alg if its starts and ends with digit.
    if (!md_len || !IS_DIGIT(md[0]) || !IS_DIGIT(md[md_len-1])) not_verified ("Malformed MD:Z field");

    buf_alloc_bitarr (vb, &vb->md_M_is_ref, vb->ref_and_seq_consumed, "md_M_is_ref"); 
    BitArrayP M_is_ref = (BitArrayP)&vb->md_M_is_ref;
    bit_array_set_all (M_is_ref); // start by marking all as matching, and clear the SNPs later
    uint64_t M_is_ref_i=0;
    
    uint32_t M_D_bases = vb->ref_consumed; // M/=/X and D
    
    bool critical_error=false;
    for_buf (BamCigarOp, op, vb->binary_cigar) {

        rom error=NULL;
        if (op->op==BC_M || op->op==BC_E || op->op==BC_X) { 
            if ((error = sam_md_consume_M (vb, is_depn, (char**)&md, &M_D_bases, &pos, op->n, M_is_ref, &M_is_ref_i, &range, &lock, &critical_error))
                && critical_error) // break loop now if critical error, else continue to count mismatch_bases_by_MD despite error
                goto not_verified;
        }

        else if (op->op==BC_D) {
            if (!IS_DIGIT (md[-1])) not_verified ("No digit before D");

            if ((error = sam_md_consume_D (vb, is_depn, (char**)&md, &M_D_bases, &pos, op->n, &range, &lock, &critical_error)) 
                && critical_error)
                goto not_verified;
        }

        else if (op->op==BC_N) { 
            pos += op->n;
            M_D_bases -= op->n;
        }

        if (!reason) reason = error;
    }

    if (*md && !(*md=='0' && !md[1])) not_verified ("Failed to consume entire MD"); // case: we didn't consume the entire MD (except an allowed trailing '0')

    if (reason) goto not_verified_buf_mismatch_count_ok; // now we can handle the non-critical error raised in sam_md_consume_M

    if (range) ref_unlock (gref, lock);

    vb->md_verified = true;    
    return; // verified

not_verified:
    vb->mismatch_bases_by_MD = -1; // fallthrough

not_verified_buf_mismatch_count_ok:    
    if (range) ref_unlock (gref, lock);
    vb->md_verified = false;

    if (flag.show_wrong_md)
        iprintf ("%s\tRNAME=%.*s\tPOS=%d\tCIGAR=%s\tMD=%.*s\tReason=%s\n", 
                 LN_NAME, STRf(vb->chrom_name), pos, vb->last_cigar, vb->last_txt_len(OPTION_MD_Z), last_txt(VB, OPTION_MD_Z), reason);
}

// MD's logical length is normally the same as seq_len, we use this to optimize it.
// In the common case that it is just a number equal the seq_len, we replace it with an empty string.
// if MD value can be derived from the seq_len, we don't need to store - store just an empty string
void sam_MD_Z_seg (VBlockSAMP vb,  ZipDataLineSAM *dl, STRp(md), unsigned add_bytes)
{
    segconf_set_has (OPTION_MD_Z);

    if (vb->md_verified && !vb->cigar_missing && !vb->seq_missing) // note: md_verified might have been reset in sam_seg_SEQ 
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_MD }, 2, OPTION_MD_Z, add_bytes);

    else // store in local. note: this is tested to be much better than splitting by length and storing short ones in dict 
        seg_add_to_local_text (VB, CTX(OPTION_MD_Z), STRa(md), true, add_bytes);
}

//---------
// PIZ
//---------

SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_MD)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    bool perfect_match = z_file->genozip_version >= 14 && !vb->mismatch_bases_by_SEQ; // "perfect" encoding introduced v14

    SamPosType pos = LOADED_CTX(SAM_POS)->last_value.i;

    ASSERTNOTEMPTY (vb->binary_cigar);
    ARRAY (BamCigarOp, cigar, vb->binary_cigar);

    ContextP sqbitmap_ctx = LOADED_CTX(SAM_SQBITMAP);
    
    BitArrayP line_sqbitmap = (BitArrayP)&sqbitmap_ctx->line_sqbitmap;
    bool is_depn = !!line_sqbitmap->nwords;
    uint32_t line_sqbitmap_next = 0;

    uint32_t save_next_local = sqbitmap_ctx->next_local;
    
    if (!perfect_match && !is_depn)
        sqbitmap_ctx->next_local = sqbitmap_ctx->last_value.i; // rewind back to beginning of the bits of this line (value stored by sam_reconstruct_SEQ)

    uint32_t count_match=0;

    for (uint32_t op_i=0; op_i < cigar_len; op_i++) { 
        uint32_t n = cigar[op_i].n;

        switch (cigar[op_i].op) {
            case BC_M : case BC_E : case BC_X :               
                if (perfect_match) {
                    count_match += n;
                    pos += n;
                }
                
                // reconstruct a series of <number><base> where number can be 0 and the base of the last pair can be missing. eg: 0T12A4
                else while (n) {

                    if (!is_depn)
                        while (n && NEXTLOCALBIT (sqbitmap_ctx)) 
                            { count_match++; n--; pos++; }
                    else 
                        while (n && bit_array_get (line_sqbitmap, line_sqbitmap_next++))
                            { count_match++; n--; pos++; }

                    if (n) {
                        // usually vb->range set in sam_reconstruct_SEQ, except if no base in SEQ matched the reference. In that case, we set it here
                        if (!vb->range) 
                            // note: if this throws an error, it is possibly related to private/defects/SAM-MD-external-ref-contig-with-no-ref.txt
                            vb->range = (RangeP)ref_piz_get_range (VB, gref, false);

                        if (reconstruct) RECONSTRUCT_INT (count_match); // flush matches before reconstructing mismatch
                        count_match=0;
                        
                        if (reconstruct) RECONSTRUCT1 (ref_base_by_pos (vb->range, pos)); 
                        n--;
                        pos++;
                    }
                }
                break;

            case BC_D :
                // flush matches before reconstructing deletion (but not if deletion is first)
                if (op_i) { 
                    if (reconstruct) RECONSTRUCT_INT (count_match);
                    count_match=0;
                }

                if (reconstruct) RECONSTRUCT1 ('^');

                if (!vb->range) 
                    vb->range = (RangeP)ref_piz_get_range (VB, gref, false);

                while (n) { 
                    if (reconstruct) RECONSTRUCT1 (ref_base_by_pos (vb->range, pos)); 
                    pos++;
                    n--;
                }
                break;

            case BC_N :  // skipping without deletion in MD
                pos += n;
                break;

            default : {}
        }
    }

    if (count_match || !IS_DIGIT (*BLSTtxt)) 
        if (reconstruct) RECONSTRUCT_INT (count_match); // flush matches if any unflushed yet

    ASSPIZ (perfect_match || save_next_local == sqbitmap_ctx->next_local, "expecting save_next_local=%u == sqbitmap_ctx->next_local=%u (RNAME=%s POS=%"PRId64" CIGAR=%s)", 
            save_next_local, sqbitmap_ctx->next_local, vb->chrom_name, CTX(SAM_POS)->last_value.i, VB_SAM->textual_cigar.data);

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
