// ------------------------------------------------------------------
//   sam_md.c
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
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

static inline bool sam_md_set_one_ref_base (VBlockSAMP vb, PosType pos, char base, uint32_t M_D_bases,                                      
                                            RangeP *range_p, RefLock *lock)
{
    // case: pos is beyond the existing range
    if ((*range_p) && (*range_p)->last_pos < pos) {
        ref_unlock (gref, *lock);
        *range_p = NULL;
    }

    // lock range
    if (! *range_p) {
        *range_p = ref_seg_get_locked_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, M_D_bases, WORD_INDEX_NONE, NULL, lock);
        if (! *range_p) return false; // cannot access this range in the reference
    }

    uint32_t pos_index = pos - (*range_p)->first_pos; // index within range

    bool internal_pos_is_populated = flag.reference == REF_INTERNAL && ref_is_nucleotide_set (*range_p, pos_index);

    // case: reference already contains a base - but unfortunately it is not "base" - so MD is not reconstractable from reference
    if (((flag.reference & REF_ZIP_LOADED) || internal_pos_is_populated) && (base != ref_base_by_pos (*range_p, pos))) return false;
    
    // case: reference is not set yet - set it now
    if (flag.reference == REF_INTERNAL && !internal_pos_is_populated)
        ref_set_nucleotide (*range_p, pos_index, base);

    // set is_set - we will need this base in the reference to reconstruct MD
    if (flag.reference & REF_STORED)
        bit_array_set (&(*range_p)->is_set, pos_index); // we will need this ref to reconstruct

    return true;
}

static inline bool sam_md_consume_D (VBlockSAMP vb, char **md_in_out, uint32_t *M_D_bases, PosType *pos, int D_bases, 
                                     RangeP *range_p, RefLock *lock)
{
    char *md = *md_in_out;

    if (! *md || *md != '^' || !IS_NUCLEOTIDE(md[1])) 
        return false; // expecting a deletion (must include at least one base)

    md++;

    while (IS_NUCLEOTIDE(*md) && D_bases) {
        if (!sam_md_set_one_ref_base (vb, *pos, *md, *M_D_bases, range_p, lock))
            return false; // base doesn't match base already in reference
        D_bases--;
        (*M_D_bases)--;
        (*pos)++;
        md++;
    }

    if (IS_NUCLEOTIDE(*md) || D_bases) return false; // CIGAR D length and MD number of deleted bases don't match

    *md_in_out = md;
    return true;
}

// verifies that the reference matches as required, and updates reference bases if missing
static inline bool sam_md_consume_M (VBlockSAMP vb, char **md_in_out, uint32_t *M_D_bases, PosType *pos, int M_bases,
                                     BitArray *M_is_ref, uint64_t *M_is_ref_i,
                                     RangeP *range_p, RefLock *lock)
{
    char *md = *md_in_out;

    // expecting a series of <number><base> where number can be 0 and the base of the last pair can be missing. eg: 0T12A4
    while (M_bases) {
        if (!IS_DIGIT(*md)) 
            return false; 

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
        // "bases" (eg N), but we apply the MD special alg only for ACTG.
        if (M_bases) {

            if (!IS_NUCLEOTIDE (*md))
                return false; // Genozip reference supports only A,C,T,G

            if (!sam_md_set_one_ref_base (vb, *pos, *md, *M_D_bases, range_p, lock))
                return false; // base doesn't match base already in reference

            bit_array_clear (M_is_ref, *M_is_ref_i); // base in SEQ is expected to be NOT equal to the reference base
            (*M_is_ref_i)++;              

            M_bases--;
            (*M_D_bases)--;
            (*pos)++;
            md++;
        }
    }

    *md_in_out = md;
    return true;
}

// called after analyzing CIGAR but before segging SEQ and later MD
// - Verifies that MD is consistent with CIGAR
// - Verifies that the mismatched and deleted bases are the same as the reference, or if not 
//   in the reference yet (in REF_INTERNAL), adds them
// - Sets md_M_is_ref, of length ref_and_seq_consumed (corresponding to M, X and = CIGAR ops): 1 for a base matching the reference, 0 for not.
//   sam_seg_seq will conduct the final verification step of comparing this bitmap to the one calculated from the SEQ data.
void sam_md_analyze (VBlockSAMP vb, STRp(md), PosType pos, const char *cigar)
{
    RangeP range = NULL;
    RefLock lock;

    if (flag.show_wrong_md)
        seg_set_last_txt (VB, CTX(OPTION_MD_Z), STRa(md)); // consumed in sam_seg_SEQ

    vb->md_verified = true; // initialize optimistically
    
    // copy of MD as we are going to modify it (but we still need the original intact for sam_md_seg)
    char md_data[md_len+1];
    memcpy (md_data, md, md_len);
    md_data[md_len] = 0;
    md = md_data;

    if (!pos || 
        (vb->chrom_name_len==1 && vb->chrom_name[0]=='*') ||
        (cigar[0] == '*' && cigar[1] == 0)) goto not_verified;

    // According to the specification (https://samtools.github.io/hts-specs/SAMtags.pdf), an MD string may start or end with a mismatch or D sequence.
    // however, in actual BAM files in the wild MD always starts and ends with a digit (possibly 0). 
    // Therefore, Genozip only activates the special MD alg if its starts and ends with digit.
    if (!md_len || !IS_DIGIT(md[0]) || !IS_DIGIT(md[md_len-1])) goto not_verified;

    buf_alloc_bitarr (vb, &vb->md_M_is_ref, vb->ref_and_seq_consumed, "md_M_is_ref"); 
    BitArray *M_is_ref = buf_get_bitarray (&vb->md_M_is_ref);
    bit_array_set_all (M_is_ref); // start by marking all as matching, and clear the SNPs later
    uint64_t M_is_ref_i=0;
    
    uint32_t M_D_bases = vb->ref_consumed; // M/=/X and D
    
    while (*cigar) {

        int subcigar_len = strtod (cigar, (char **)&cigar); // get number and advance next_cigar
        char cigar_op = *(cigar++);

        if (cigar_op=='M' || cigar_op=='=' || cigar_op=='X') {
            
            if (!sam_md_consume_M (vb, (char**)&md, &M_D_bases, &pos, subcigar_len, M_is_ref, &M_is_ref_i, &range, &lock))
                goto not_verified;
        }

        else if (cigar_op=='D') {
            if (!IS_DIGIT (md[-1])) goto not_verified; 

            if (!sam_md_consume_D (vb, (char**)&md, &M_D_bases, &pos, subcigar_len, &range, &lock))
                goto not_verified;
        }

        else if (cigar_op=='N') { 
            pos += subcigar_len;
            M_D_bases -= subcigar_len;
        }
    }

    if (*md && !(*md=='0' && !md[1])) goto not_verified; // case: we didn't consume the entire MD (except an allowed trailing '0')

    if (range) ref_unlock (gref, lock);
    return;

not_verified:
    if (range) ref_unlock (gref, lock);
    vb->md_verified = false;

    if (flag.show_wrong_md)
        iprintf ("vb=%u line=%"PRIu64" RNAME=%.*s POS=%"PRId64" CIGAR=%s MD=%.*s Special MD algorithm doesn't work on this MD (no harm)\n", 
                vb->vblock_i, vb->line_i, STRf(vb->chrom_name), pos, vb->last_cigar, vb->last_txt_len(OPTION_MD_Z), last_txt(vb, OPTION_MD_Z));
}

// MD's logical length is normally the same as seq_len, we use this to optimize it.
// In the common case that it is just a number equal the seq_len, we replace it with an empty string.
// if MD value can be derived from the seq_len, we don't need to store - store just an empty string
void sam_md_seg (VBlockSAM *vb,  ZipDataLineSAM *dl, STRp(md), unsigned add_bytes)
{
    segconf_set_has (OPTION_MD_Z);

    if (vb->md_verified) 
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_MD}, 2, OPTION_MD_Z, add_bytes);
    else
        seg_by_did_i (VB, STRa(md), OPTION_MD_Z, add_bytes);
}

//---------
// PIZ
//---------

SPECIAL_RECONSTRUCTOR (sam_piz_special_MD)
{
    VBlockSAMP vb_sam = VB_SAM;

    PosType pos = CTX(SAM_POS)->last_value.i;

    ASSERTISALLOCED (vb_sam->textual_cigar);
    const char *cigar = FIRSTENT (char, vb_sam->textual_cigar);  // note: we can't use last_txt as it is already translated to binary in BAM (snips in dict/local are textual CIGARs)

    ContextP sqbitmap_ctx = CTX(SAM_SQBITMAP);
    uint32_t save_next_local = sqbitmap_ctx->next_local;
    sqbitmap_ctx->next_local = sqbitmap_ctx->last_value.i; // rewind back to beginning of the bits of this line (value stored by sam_reconstruct_SEQ)

    uint32_t count_match=0;
    for (uint32_t op_i=0; *cigar && *cigar != '\t' && *cigar != '\n'; op_i++) { 

        int subcigar_len = strtod (cigar, (char **)&cigar); // get number and advance next_cigar
        char cigar_op = *(cigar++);

        if (cigar_op=='M' || cigar_op=='=' || cigar_op=='X') {
            
            // reconstruct a series of <number><base> where number can be 0 and the base of the last pair can be missing. eg: 0T12A4
            while (subcigar_len) {

                while (subcigar_len && NEXTLOCALBIT (sqbitmap_ctx)) {
                    count_match++;
                    subcigar_len--;
                    pos++;
                }

                if (subcigar_len) {
                    // case: no base in SEQ matched the reference - sam_reconstruct_SEQ set vb->range to NULL, so we set it here
                    if (!vb->range) {
                        vb->range = (RangeP)ref_piz_get_range (vb, gref, vb->last_int(SAM_POS), VB_SAM->ref_consumed);
                        ASSERTNOTNULL (vb->range);
                    }

                    RECONSTRUCT_INT (count_match); // flush matches before reconstructing mismatch
                    count_match=0;
                    
                    RECONSTRUCT1 (ref_base_by_pos (vb->range, pos)); // vb->range set in sam_reconstruct_SEQ
                    subcigar_len--;
                    pos++;
                }
            }
        }

        else if (cigar_op=='D') {
            
            // flush matches before reconstructing deletion (but not if deletion is first)
            if (op_i) { 
                RECONSTRUCT_INT (count_match);
                count_match=0;
            }

            RECONSTRUCT1 ('^');
            while (subcigar_len) { 
                RECONSTRUCT1 (ref_base_by_pos (vb->range, pos)); 
                pos++;
                subcigar_len--;
            }
        }

        else if (cigar_op=='N')  // skipping without deletion in MD
            pos += subcigar_len;

    }

    if (count_match || !IS_DIGIT (*LASTENT (char, vb->txt_data))) 
        RECONSTRUCT_INT (count_match); // flush matches if any unflushed yet

    ASSPIZ (save_next_local == sqbitmap_ctx->next_local, "expecting save_next_local=%u == sqbitmap_ctx->next_local=%u (buddy_line_i=%d RNAME=%s POS=%"PRId64" CIGAR=%s)", 
            save_next_local, sqbitmap_ctx->next_local, vb->buddy_line_i, vb->chrom_name, CTX(SAM_POS)->last_value.i, VB_SAM->textual_cigar.data);

    return false; // no new value
}

// Used in files compressed with Genozip up to 12.0.36 - 
// logic: snip is eg "119C" (possibly also "") - we reconstruct the original, eg "119C31" 
// by concating a number which is (seq_len - partial_seq_len_by_md_field)
SPECIAL_RECONSTRUCTOR (sam_piz_special_MD_old)
{
    if (!reconstruct) return false;
    
    if (snip_len) RECONSTRUCT (snip, snip_len);

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

    return false; // no new value
}
