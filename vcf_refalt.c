// ------------------------------------------------------------------
//   vcf_refalt.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <math.h>
#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "strings.h"
#include "codec.h"
#include "reconstruct.h"
#include "dict_id.h"
#include "file.h"
#include "reference.h"
#include "chain.h"
#include "coords.h"
#include "ref_iupacs.h"

// ---------
// Seg stuff
// ---------

// optimize REF and ALT, for simple one-character REF/ALT (i.e. mostly a SNP or no-variant)
static void vcf_refalt_seg_ref_alt_snp (VBlockVCFP vb, char main_ref, char main_alt)
{
    char new_ref=0, new_alt=0;

    // if we have a reference, we use it (we treat the reference as PRIMARY)
    // except: if --match-chrom, we assume the user just wants to match, and we don't burden him with needing the reference to decompress
    if (((flag.reference == REF_EXTERNAL && !flag.match_chrom_to_reference) || flag.reference == REF_EXT_STORE) && vb->line_coords == DC_PRIMARY) {
        PosType pos = vb->last_int(VCF_POS);

        RefLock lock;

        Range *range = ref_seg_get_locked_range (VB, gref, vb->chrom_node_index, vb->chrom_name, vb->chrom_name_len, pos, 1, WORD_INDEX_NONE, NULL, &lock);
        if (range) { // this chrom is in the reference
            uint32_t index_within_range = pos - range->first_pos;

            ref_assert_nucleotide_available (range, pos);
            char ref = ref_base_by_idx (range, index_within_range);

            if (main_ref == ref) new_ref = '-'; // this should always be the case...
            if (main_alt == ref) new_alt = '-'; 

            if (flag.reference == REF_EXT_STORE)
                bit_array_set (&range->is_set, index_within_range);

            ref_unlock (gref, lock);
        }
    }

    // replace the most common SNP with +
    // based on counting simple SNPs from from chr22 of 1000 genome project phase 3:
    // G to: A=239681 C=46244 T=44084
    // C to: T=238728 G=46508 A=43685
    // A to: G=111967 C=30006 T=26335
    // T to: C=111539 G=29504 A=25599

    if      (main_alt == 'A' && main_ref == 'G') new_alt = '+';
    else if (main_alt == 'C' && main_ref == 'T') new_alt = '+';
    else if (main_alt == 'G' && main_ref == 'A') new_alt = '+';
    else if (main_alt == 'T' && main_ref == 'C') new_alt = '+';

    // if anything was done, we create a "special" snip
    if (new_ref || new_alt) {
        char refalt_special[4] = { SNIP_SPECIAL, VCF_SPECIAL_main_REFALT };
        refalt_special[2] = new_ref ? new_ref : main_ref;
        refalt_special[3] = new_alt ? new_alt : main_alt;

        seg_by_did_i (VB, refalt_special, sizeof(refalt_special), SEL (VCF_REFALT, VCF_oREFALT), 0); // we do the accounting in vcf_refalt_seg_main_ref_alt
    }
    // if not - just the normal snip
    else {
        char refalt_normal[3] = { 0, '\t', 0 };
        refalt_normal[0] = main_ref;
        refalt_normal[2] = main_alt;

        seg_by_did_i (VB, refalt_normal, sizeof(refalt_normal), SEL (VCF_REFALT, VCF_oREFALT), 0); // we do the account in vcf_refalt_seg_main_ref_alt
    }
}

void vcf_refalt_seg_main_ref_alt (VBlockVCFP vb, STRp(ref), STRp(alt))
{
    // optimize ref/alt in the common case of single-character
    if (ref_len == 1 && alt_len == 1) 
        vcf_refalt_seg_ref_alt_snp (vb, *ref, *alt);

    else {
        char ref_alt[ref_len + 1 + alt_len];
        memcpy (ref_alt, ref, ref_len);
        ref_alt[ref_len] = '\t';
        memcpy (&ref_alt[ref_len+1], alt, alt_len);

        seg_by_did_i (VB, ref_alt, ref_len + alt_len + 1, SEL (VCF_REFALT, VCF_oREFALT), 0);
    }
        
    if (vb->line_coords == DC_PRIMARY) // in the default reconstruction, this REFALT is in the main fields along with 2 tabs
        CTX(VCF_REFALT)->txt_len += ref_len + alt_len + 2;

    else if (vb->vb_coords == DC_BOTH) // LUFT line in dual coord VB - in the default reconstruction of dual-coord line, oREF is in the INFO/LUFT vector
        CTX(VCF_LIFT_REF)->txt_len += ref_len;

    else // ##luft-only VB - we account for the main fields and 2 tabs
        CTX(VCF_oREFALT)->txt_len += ref_len + alt_len + 2;

    // note: vb->recon_size doesn't change if txt_files.coords==LUFT, because currently Genozip only supports dual coords if the size
    // of REF and oREF are the same, and ALT and oALT.
}

// ----------------------------
// Lifting with genozip --chain
// ----------------------------

// we accept an INDEL or complex variant for REF⇄ALT switch if (1) it is short enough (2) flanking regions are the same
#define FLANKING_SEQ_LEN 4
#define MAX_LEN_REF_ALT_SWITCH_INDEL_BASED_ON_FLANKING 16

#define MAX_BASES_REJECTS_FILE 256
#define ECLIPSE(len) ((len) > MAX_BASES_REJECTS_FILE ? "..." : "")
#define ECLIPSED(seq,len) MIN_(MAX_BASES_REJECTS_FILE,(len)), seq, ECLIPSE(len)
#define XSTRANDF "%s"
#define XSTRAND (is_xstrand ? "xstrand=true " : "")

#define IS_SINGLE_BASE_ALTS(n_alts, alt_len) ((n_alts)*2-1 == (alt_len))

#define ALTF "%.*s%s\t"
#define ALT ECLIPSED (alt, alt_len)

#define DEF_PRIM(len) char prim_str[MIN_((MAX_BASES_REJECTS_FILE+1), (len)+1)]; unsigned prim_str_len = (unsigned)(len)
#define PRIMF "PRIM=%s%s "
#define PRIM ref_dis_subrange (prim_ref, prim_range, pos, sizeof (prim_str), prim_str, false), ECLIPSE (prim_str_len)
#define PRIMFLANKING ref_dis_subrange (prim_ref, (prim_range), (pos)-FLANKING_SEQ_LEN, sizeof (prim_str), (prim_str), false), ECLIPSE (prim_str_len)

#define DEF_LUFT(len) char luft_str[MIN_((MAX_BASES_REJECTS_FILE+1), (len)+1)]; unsigned luft_str_len = (unsigned)(len)
#define LUFTF "LUFT(%"PRId64")=%s%s "
#define LUFT (opos), ref_dis_subrange (gref, (luft_range), (opos), sizeof (luft_str), (luft_str), (is_xstrand)), ECLIPSE (luft_str_len)
#define LUFTFLANKING (opos)-FLANKING_SEQ_LEN, ref_dis_subrange (gref, (luft_range), (opos)-FLANKING_SEQ_LEN, sizeof (luft_str), (luft_str), (is_xstrand)), ECLIPSE (luft_str_len)

#define SAME_AS_LUFT_REF(vcf_base) is_same_base (vb, luft_range, opos, (vcf_base), is_xstrand)
static inline bool is_same_base (VBlockVCFP vb, const Range *luft_range, PosType opos, char vcf_base, bool is_xstrand)
{
    if (is_xstrand) vcf_base = COMPLEM[(int)vcf_base];

    char ref_base = ref_base_by_pos (luft_range, opos); // counting on compiler optimizer to avoid re-calculating if already calculated

    return vcf_base == ref_base ||
           ref_iupacs_is_included (gref, VB, luft_range, opos, vcf_base);
}

// false is not the same or if either goes beyond the end of the range 
static inline bool is_same_seq (Reference ref, VBlockVCFP vb, const Range *range, PosType pos, const char *seq, PosType seq_len, bool is_xstrand) 
{
    if (!is_xstrand) {
        if (pos + seq_len - 1 > range->last_pos) return false; // seq goes beyond the end of range

        for (PosType i=0; i < seq_len ; i++) {
            char vcf_base = UPPER_CASE (seq[i]);
            char ref_base = ref_base_by_pos (range, pos + i);
            if (ref_base != vcf_base && !ref_iupacs_is_included (ref, VB, range, pos, vcf_base)) return false;
        }
    }
    else { // xstrand
        if (pos < seq_len) return false; // seq goes beyond the start of range

        for (PosType i=0; i < seq_len ; i++) {
            char vcf_base = UPPER_COMPLEM[(int)seq[i]];
            char ref_base = ref_base_by_pos (range, pos - i);
            if (ref_base != vcf_base && !ref_iupacs_is_included (ref, VB, range, pos, vcf_base)) return false;
        }
    }

    return true;
}

// false is not the same or if either goes beyond the end of the range 
static inline bool is_same_refs (const Range *prim_range, PosType prim_pos, 
                                 const Range *luft_range, PosType luft_pos, 
                                 unsigned seq_len, bool is_xstrand) 
{
    if (prim_pos + seq_len - 1 > prim_range->last_pos) return false; // seq goes beyond the end of range

    if (!is_xstrand) {
        if (luft_pos + seq_len - 1 > luft_range->last_pos) return false; // seq goes beyond the end of range

        for (PosType i=0; i < seq_len ; i++) 
            if (ref_base_by_pos (prim_range, prim_pos + i) != ref_base_by_pos (luft_range, luft_pos + i)) return false;
    }
    else { // xstrand
        if (luft_pos < seq_len) return false; // seq goes beyond the start of range

        for (PosType i=0; i < seq_len ; i++)
            if (ref_base_by_pos (prim_range, prim_pos + i) != COMPLEM[(int)ref_base_by_pos (luft_range, luft_pos - i)]) return false;
    }

    return true;
}

// counts repeated bases (forwards or backwards depending on direction)
static unsigned vcf_refalt_lift_get_repeats (const Range *range, const char *payload, PosType payload_len, PosType pos, bool revcomp)
{
    unsigned num_reps = 0;
    if (!revcomp) {
        for (PosType p=pos; p <= range->last_pos; p++, num_reps++) 
            if (ref_base_by_pos (range, p) != payload[num_reps % payload_len]) break;
    }
    else {
        for (PosType p=pos; p >= 1; p--, num_reps++) 
            if (ref_base_by_pos (range, p) != COMPLEM[(int)payload[num_reps % payload_len]]) break;
    }

    return num_reps;
}

// we accept an INDEL for REF⇄ALT switch if (1) it is short enough (2) flanking regions are the same
static bool vcf_refalt_lift_same_flanking_regions (bool is_xstrand,
                                                   const Range *prim_range, PosType pos, PosType ref_len,
                                                   const Range *luft_range, PosType opos, PosType alt_len)
{
    if (ref_len > MAX_LEN_REF_ALT_SWITCH_INDEL_BASED_ON_FLANKING || alt_len > MAX_LEN_REF_ALT_SWITCH_INDEL_BASED_ON_FLANKING)
        return false; // INDEL too long to be considered for REF⇄ALT switch on basis of flanking regions
        
    if (!is_xstrand) 
        return is_same_refs (prim_range, pos + ref_len, luft_range, opos + alt_len, FLANKING_SEQ_LEN, false) &&
               is_same_refs (prim_range, pos - FLANKING_SEQ_LEN, luft_range, opos - FLANKING_SEQ_LEN, FLANKING_SEQ_LEN, false);

    else
        return is_same_refs (prim_range, pos + ref_len, luft_range, opos - 1, FLANKING_SEQ_LEN, true) &&
               is_same_refs (prim_range, pos - FLANKING_SEQ_LEN, luft_range, opos + alt_len + FLANKING_SEQ_LEN, FLANKING_SEQ_LEN, true);
}

// if is_xstrand, change anchor base to the right of the payload (so it ends up on the left after revcomp due to xstrand)
// note: we only change the anchor to be to the left of the rev-comped oREF, we don't left align as that would
// risk data corruption, see Supplamentary Information for an explanation why
static inline LiftOverStatus 
    vcf_refalt_xstrand_flip_anchor (VBlockVCFP vb, 
                                    const Range *prim_range, PosType pos, unsigned ref_len,
                                    const Range *luft_range, PosType opos)
{
    PosType last_aln_pos = chain_get_aln_last_pos (vb->pos_aln_i);

    REJECTIF (last_aln_pos < pos + ref_len, LO_NO_MAPPING_IN_CHAIN, ".\tNew right-anchor base is beyond chain file alignment: new_achor_pos=%"PRId64" last_aln_pos=%"PRId64,
              pos + ref_len, last_aln_pos);

    // new anchor - the base to the left of the oREF
    vb->new_ref = ref_base_by_pos (luft_range, opos - ref_len); // ... and the luft reference
    
    return LO_OK;                
}
#define IF_XSTRAND_FLIP_ANCHOR \
    if (is_xstrand) {\
        LiftOverStatus ostatus = vcf_refalt_xstrand_flip_anchor (vb, prim_range, pos, ref_len, luft_range, opos);\
        if (ostatus != LO_OK) return ostatus; /* rejected by vcf_refalt_xstrand_flip_anchor */\
    }

static inline LiftOverStatus vcf_refalt_lift_report_indel_outcome (VBlockVCFP vb, STRp(alt),
                                                                   const Range *prim_range, PosType pos, const Range *luft_range, PosType opos,
                                                                   unsigned num_alts, bool is_xstrand,
                                                                   int rep_len_prim, int rep_len_luft, bool all_same, int switch_alt_i, int maybe_switch_alt_i,
                                                                   PosType ref_len, PosType switch_alt_len, PosType maybe_switch_alt_len)
{    
    // we can lift a Insertion IF: either - all ALTs agree that it is REF_SAME ; or - single ALT and REF⇄ALT switch
    if (all_same) {
        IF_XSTRAND_FLIP_ANCHOR;
        LIFTOK (LO_OK_REF_SAME_INDEL, ".\tINDEL REF unchanged " XSTRANDF, XSTRAND);
    }

    DEF_PRIM (rep_len_prim + 2);
    DEF_LUFT (rep_len_luft + 2);

    if (switch_alt_i==0 && num_alts == 1) {
        IF_XSTRAND_FLIP_ANCHOR;
        LIFTOK (LO_OK_REF_ALT_SWITCH_INDEL, ALTF "INDEL REF<>ALT switch with repeats: " XSTRANDF PRIMF LUFTF, ALT, XSTRAND, PRIM, LUFT); 
    }

    // encountered a switch ALT in a multi allelic DEL variant
    REJECTIF (switch_alt_i >= 0, LO_REF_MULTIALT_SWITCH_INDEL, ALTF "Genozip limitation: " PRIMF " and ALT[%d] matches " LUFTF ", but can't switch in multi-allelic INDEL variant", 
              ALT, PRIM, switch_alt_i+1, LUFT);

    // case: possibly a REF⇄ALT switch, call based on flanking regions
    if (maybe_switch_alt_i >= 0) {
        bool same_flanking_regions = vcf_refalt_lift_same_flanking_regions (is_xstrand, prim_range, pos, ref_len, luft_range, opos, maybe_switch_alt_len);

        DEF_PRIM (rep_len_prim + FLANKING_SEQ_LEN*2 + 2*(ref_len > 1)); 
        DEF_LUFT (rep_len_luft + FLANKING_SEQ_LEN*2 + 2*(ref_len > 1));

        REJECTIF (!same_flanking_regions, LO_NEW_ALLELE_INDEL, ALTF "Genozip limitation: REF switched with ALT[%d], but flanking regions differ. " PRIMF LUFTF "(flanking %d bases on either side)", 
                ALT, maybe_switch_alt_i+1, PRIMFLANKING, LUFTFLANKING, FLANKING_SEQ_LEN); 

        REJECTIF (num_alts > 1, LO_REF_MULTIALT_SWITCH_INDEL, ALTF "Genozip limitation: " PRIMF " and ALT[%d] matches " LUFTF "(flanking %d bases on either side), but can't switch in multi-allelic INDEL variant", 
                ALT, PRIMFLANKING, maybe_switch_alt_i+1, LUFTFLANKING, FLANKING_SEQ_LEN);

        IF_XSTRAND_FLIP_ANCHOR;

        LIFTOK (LO_OK_REF_ALT_SWITCH_INDEL, ALTF "INDEL REF<>ALT switch with matching flanking regions. " XSTRANDF PRIMF LUFTF "(flanking %d bases on either side)", 
                ALT, XSTRAND, PRIMFLANKING, LUFTFLANKING, FLANKING_SEQ_LEN); 
    }

    // LUFT is a new allele
    REJECTIF (true, LO_NEW_ALLELE_INDEL, ALTF "Genozip limitation: new INDEL allele: " PRIMF LUFTF XSTRANDF "num_alts=%u rep_len_prim=%u rep_len_luft=%u switch_alt_i=%d maybe_switch_alt_i=%d", 
              ALT, PRIM, LUFT, XSTRAND, num_alts, rep_len_prim, rep_len_luft, switch_alt_i, maybe_switch_alt_i); 

    return 0; // never reaches here
}

// returns true is REF is the "same" taking into account that is is_xstrand the anchor base shifts sides, and the payload rev-comps
static bool vcf_refalt_is_REF_same_in_luft (VBlockVCFP vb, STRp(ref), bool is_xstrand,
                                            const Range *prim_range, PosType pos,
                                            const Range *luft_range, PosType opos)
{
    // if LUFT and PRIMARY is not the same for ref_len bases: example: "GAC G" Primary: "GAC" 
    // Luft=GTT: REF⇄ALT switch: "G GAC"
    char anchor = is_xstrand ? ref_base_by_pos (prim_range, pos + ref_len) : ref[0]; // if is_xstrand, we re-anchor REF to be right-anchored
    PosType anchor_opos = is_xstrand ? opos - ref_len : opos; 

    bool same_anchor  = is_same_base (vb, luft_range, anchor_opos, anchor, is_xstrand); 
    bool same_payload = ref_len > 1 ? is_same_seq (gref, vb, luft_range, opos + (is_xstrand ? -1 : 1), ref+1, ref_len-1, is_xstrand) : true;

    return same_anchor && same_payload;
}

// when --chain: handle a Deletion
static LiftOverStatus vcf_refalt_lift_deletion (VBlockVCFP vb, const char *ref, PosType ref_len,
                                                unsigned num_alts, const char **alts, STRp(alt), 
                                                bool is_xstrand,
                                                const Range *prim_range, PosType pos, 
                                                const Range *luft_range, PosType opos)                                        
{
    if (!vcf_refalt_is_REF_same_in_luft (vb, ref, ref_len, is_xstrand, prim_range, pos, luft_range, opos)) {
        REJECTIF (num_alts > 1, LO_NEW_ALLELE_INDEL, ALTF "Genozip limitation: REF changed in a Deletion variant that has more than one ALT", ALT);

        { DEF_PRIM (ref_len); DEF_LUFT (ref_len);
          REJECTIF (!is_same_seq (gref, vb, luft_range, opos, alts[0], 1, is_xstrand), LO_NEW_ALLELE_INDEL, 
                    ALTF "Genozip limitation: REF changed in Deletion variant, but not to ALT. " PRIMF LUFTF, ALT, PRIM, LUFT); }

        // if the REF is not in the Luft range, we check if the flanking regions of the ALT in Luft are the same as REF's to decide
        // whether this is a REF⇄ALT switch
        bool same_flanking_regions = vcf_refalt_lift_same_flanking_regions (is_xstrand, prim_range, pos, ref_len, luft_range, opos, 1);

        DEF_PRIM (ref_len + FLANKING_SEQ_LEN*2); 
        DEF_LUFT (1 + FLANKING_SEQ_LEN*2);
        REJECTIF (!same_flanking_regions, LO_NEW_ALLELE_INDEL, 
                  ALTF "Genozip limitation: REF changed in Deletion variant to ALT, but this is not a REF<>ALT switch because flanking regions differ " PRIMF LUFTF "(shown with flanking %d on either side)",
                  ALT, PRIMFLANKING, LUFTFLANKING, FLANKING_SEQ_LEN);

        LIFTOK (LO_OK_REF_ALT_SWITCH_INDEL, ALTF "REF<>ALT INDEL switch with matching flanking region: " PRIMF LUFTF "(shown with flanking %d on either side)", 
                ALT, PRIMFLANKING, LUFTFLANKING, FLANKING_SEQ_LEN); // simple REF⇄ALT switch - deletion already included in Luft
    }

    // now we know that REF is stringwise the same as the reference, but we need to test for repeats

    const char *payload = ref + 1;
    PosType payload_len = ref_len-1;
    
    unsigned rep_len_prim=0, rep_len_luft=0;
    bool all_same=true;
    int switch_alt_i=-1, maybe_switch_alt_i=-1;

    for (unsigned alt_i=0; alt_i < num_alts; alt_i++) {

        rep_len_prim = rep_len_luft = payload_len; // we already verified that LUFT and PRIMARY are the same - those count as repeats

        // left-anchored eg "GAC G" - count payload repeats ('G' being the "achor base" and 'AC' being the payload)
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACACT   (4 repeats)  -> REF/ALT Unchanged
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACT     (3 repeats)  -> REF <> ALT switch. Luft: G GAC
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACACACT (5+ repeats) -> new allele: GACACAC G,GAC (old REF was a 4-repeats: now GAC   ; old ALT was a 3-repeats: now G)
        // PRIM: GACACACACT (4 repeats) LUFT: GACACT       (2- repeats) -> new allele: G GAC,GACAC   (old REF was a 4-repeats: now GACAC ; old ALT was a 3-repeats: now GAC)
        
        rep_len_prim +=               vcf_refalt_lift_get_repeats (prim_range, payload, payload_len, pos  + ref_len, false); // count EXCLUDING this the present repeat (already set to 1)
        rep_len_luft += !is_xstrand ? vcf_refalt_lift_get_repeats (luft_range, payload, payload_len, opos + ref_len, false)
                                    : vcf_refalt_lift_get_repeats (luft_range, payload, payload_len, opos - ref_len, true);

        if (rep_len_luft != rep_len_prim) all_same = false; // at least one of the ALTs is not REF_SAME
        
        if (rep_len_luft == rep_len_prim - payload_len) {
            if (rep_len_luft >= payload_len)  // this is a region with repeats 
                switch_alt_i = alt_i; // an ALT that is definitely an REF⇄ALT switch
            else
                maybe_switch_alt_i = alt_i; // possible REF⇄ALT switch without repeats - need more evidence by comparing flanking regions
        } 
    }

    return vcf_refalt_lift_report_indel_outcome (vb, alt, alt_len, prim_range, pos, luft_range, opos, 
                                                 num_alts, is_xstrand, rep_len_prim, rep_len_luft, all_same, 
                                                 switch_alt_i, maybe_switch_alt_i, ref_len, 1, 1);
}

// when --chain: handle an Insertion
static LiftOverStatus vcf_refalt_lift_insertion (VBlockVCFP vb, STRp(ref), /* must be 1 - needed for macros */ 
                                                 unsigned num_alts, const char **alts, const unsigned *alt_lens, STRp(alt), 
                                                 bool is_xstrand,
                                                 const Range *prim_range, PosType pos, 
                                                 const Range *luft_range, PosType opos)
                                                
{    
    unsigned rep_len_prim=0, rep_len_luft=0; 

    bool all_same=true;
    int switch_alt_i=-1, maybe_switch_alt_i=-1; // the ALT number that is a REF⇄ALT switch
    for (unsigned alt_i=0; alt_i < num_alts; alt_i++) {

        rep_len_prim = rep_len_luft = 0; 

        // check if anchor remains - eg "A AGC" (left anchored) or "A GCA" (right anchored, ie is_xstrand - anchor ia A)
        if (!vcf_refalt_is_REF_same_in_luft (vb, ref, 1, is_xstrand, prim_range, pos, luft_range, opos)) {
            all_same = false;
            continue;
        }

        // left-anchored eg "GAC G" - count payload repeats ('G' being the "achor base" and 'AC' being the payload)
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACACT   (4 repeats)  -> REF/ALT Unchanged
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACT     (3 repeats)  -> REF <> ALT switch. Luft: G GAC
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACACACT (5+ repeats) -> new allele: GACACAC G,GAC (old REF was a 4-repeats: now GAC   ; old ALT was a 3-repeats: now G)
        // PRIM: GACACACACT (4 repeats) LUFT: GACACT       (2- repeats) -> new allele: G GAC,GACAC   (old REF was a 4-repeats: now GACAC ; old ALT was a 3-repeats: now GAC)

        const char *payload = alts[alt_i] + 1;
        PosType payload_len = alt_lens[alt_i]-1;

        rep_len_prim +=               vcf_refalt_lift_get_repeats (prim_range, payload, payload_len, pos  + 1, false); 
        rep_len_luft += !is_xstrand ? vcf_refalt_lift_get_repeats (luft_range, payload, payload_len, opos + 1, false)
                                    : vcf_refalt_lift_get_repeats (luft_range, payload, payload_len, opos - 1, true);
    
        if (rep_len_luft != rep_len_prim) all_same = false; // at least one of the ALTs is not REF_SAME
        
        if (rep_len_luft == rep_len_prim + payload_len) {
            if (rep_len_prim >= payload_len) // this is a region with repeats            
                switch_alt_i = alt_i; // found an ALT that is definitely a REF⇄ALT switch
            else
                maybe_switch_alt_i = alt_i; // possible REF⇄ALT switch without repeats - need more evidence by comparing flanking regions
        } 
    }

    return vcf_refalt_lift_report_indel_outcome (vb, alt, alt_len, prim_range, pos, luft_range, opos, 
                                                 num_alts, is_xstrand, rep_len_prim, rep_len_luft, all_same, 
                                                 switch_alt_i, maybe_switch_alt_i, 1, 
                                                 switch_alt_i >= 0 ? alt_lens[switch_alt_i] : 0, 
                                                 maybe_switch_alt_i >= 0 ? alt_lens[maybe_switch_alt_i] : 0);
}

// when --chain: handle an structural variant with an ALT containing a <***> ID  (eg <DEL>)
static LiftOverStatus vcf_refalt_lift_with_sym_allele (VBlockVCFP vb, STRp(ref), 
                                          unsigned num_alts, const char **alts, const unsigned *alt_lens, STRp(alt),
                                          bool is_xstrand,
                                          const Range *prim_range, PosType pos, 
                                          const Range *luft_range, PosType opos)
{
    REJECTIF (is_xstrand, LO_XSTRAND_SV, ALTF "Genozip limitation: Variant with a symbolic allele is mapped to the reverse strand" XSTRANDF, ALT, XSTRAND);

    // case: structural variant - but REF unchanged 
    if (is_same_seq (gref, vb, luft_range, opos, ref, ref_len, is_xstrand))
        LIFTOK0 (LO_OK_REF_SAME_SV, ".\tREF unchanged - structural variant");

    else
        REJECTIF (true, LO_NEW_ALLELE_SV, ALTF "Genozip limitation: REF changed in variant with a symbolic allele", ALT);

    return 0; // never reaches here
}

// when --chain: handle a complex variant (anything either than SNPs, left-anchored indels, complex rearrangements or symbolic ALTs)
static LiftOverStatus vcf_refalt_lift_complex (VBlockVCFP vb, STRp(ref), 
                                               unsigned num_alts, const char **alts, const unsigned *alt_lens, STRp(alt),
                                               bool is_xstrand, bool is_left_anchored,
                                               const Range *prim_range, PosType pos, 
                                               const Range *luft_range, PosType opos)
{
    // rejecting a non-left-anchored complex with xstrand, because vcf_lo_seg_ostatus_from_LUFT_or_PRIM can't tell if its a
    // left-anchored indel vs a non-left-anchored complex
    //
    { DEF_LUFT (ref_len);
      REJECTIF (is_xstrand && !is_left_anchored, LO_XSTRAND_NLA, 
                ALTF "Genozip limitation: " LUFTF "strand reversal not supported for non-left-anchored indels", ALT, LUFT); }

    // case: complex variant - but REF is unchanged (note: we don't check flanking regions, it is debatable if we should)
    if (vcf_refalt_is_REF_same_in_luft (vb, ref, ref_len, is_xstrand, prim_range, pos, luft_range, opos)) {
        if (is_left_anchored) {
            IF_XSTRAND_FLIP_ANCHOR;
            LIFTOK0 (LO_OK_REF_SAME_INDEL, ".\tREF unchanged - left-anchored complex INDEL");
        }
        else 
            LIFTOK0 (LO_OK_REF_SAME_NLA, ".\tREF unchanged - non-left-anchored COMPLEX variant");
    }

    // case: complex variant - REF⇄ALT switch - supported only if bi-allelic and not SV
    else if (num_alts==1) {
        if (is_same_seq (gref, vb, luft_range, opos, alts[0], alt_lens[0], is_xstrand) &&
            vcf_refalt_lift_same_flanking_regions (is_xstrand, prim_range, pos, ref_len, luft_range, opos, alt_lens[0])) {
            if (is_left_anchored) {
                IF_XSTRAND_FLIP_ANCHOR;
                LIFTOK (LO_OK_REF_ALT_SWITCH_INDEL, ALTF "REF<>ALT switch, left-anchored complex INDEL", ALT);
            }
            else
                LIFTOK (LO_OK_REF_ALT_SWITCH_NLA, ALTF "REF<>ALT switch, non-left-anchored COMPLEX variant", ALT);
        }
    }

    // case: complex variant - multiallelic REF⇄ALT switch - detect, but unsupported 
    else
        for (unsigned alt_i=0; alt_i < num_alts; alt_i++)
            if (is_same_seq (gref, vb, luft_range, opos, alts[alt_i], alt_lens[alt_i], is_xstrand) &&
                vcf_refalt_lift_same_flanking_regions (is_xstrand, prim_range, pos, ref_len, luft_range, opos, alt_lens[alt_i])) {
                DEF_LUFT (alt_lens[alt_i]);
                REJECTIF (true, LO_REF_MULTIALT_SWITCH_INDEL, ALTF "Genozip limitation: " LUFTF "== ALT[%d], but can't switch in multi-allelic variant", 
                          ALT, LUFT, alt_i+1);
            }

    DEF_LUFT (ref_len);
    REJECTIF (true, is_left_anchored ? LO_NEW_ALLELE_INDEL : LO_NEW_ALLELE_NLA, ALTF "Genozip limitation: " LUFTF "is a new allele", ALT, LUFT);

    return 0; // never reaches here
}

// an ineffecient way of parsing a VCF line to extract AF - used in exceptional cases. Returns AF if its a valid single number, or -1 otherwise
static double vcf_refalt_get_INFO (VBlockVCFP vb, const char *tag)
{
    // split the VCF line into fields
    const char *line = ENT (char, vb->txt_data, vb->line_start);
    str_split (line, strcspn (line, "\n\r"), 8 + (vcf_num_samples>0) + vcf_num_samples, '\t', field, false);
    if (n_fields < 8) return -1; // no INFO field (we will deal with this error in vcf_seg_txt_line)

    // split INFO fields into subfields
    str_split (fields[7], field_lens[7], 0, ';', subfield, true);

    // look for tag
    for (unsigned sf=0; sf < n_subfields; sf++) 
        if (!memcmp (subfields[sf], tag, 3)) {
            SAFE_NUL (&subfields[sf][subfield_lens[sf]]);
            char *after;
            double value = strtod (&subfields[sf][3], &after);
            SAFE_RESTORE;

            return (subfield_lens[sf] && after == subfields[sf] + subfield_lens[sf]) ? value : -1; // not a valid single-number AF field
        }

    return -1; // AF not found in INFO
}

// when --chain: handle a SNP (one or more ALTs)
static LiftOverStatus vcf_refalt_lift_snp (VBlockVCFP vb, STRp(ref)/* must be 1 */, 
                                           unsigned num_alts, const char **alts, const unsigned *alt_lens, STRp(alt),
                                           bool is_xstrand,
                                           const Range *prim_range, PosType pos, 
                                           const Range *luft_range, PosType opos)
{
    char oref = ref_base_by_pos (luft_range, opos); // range is 0-based

    if (SAME_AS_LUFT_REF (ref[0])) 
        LIFTOK0 (LO_OK_REF_SAME_SNP, ".\tSNP: REF unchanged");

    if (num_alts == 1 && SAME_AS_LUFT_REF (alts[0][0])) 
        LIFTOK (LO_OK_REF_ALT_SWITCH_SNP, ALTF "SNP: REF<>ALT switch", ALT);

    for (unsigned alt_i=0; alt_i < num_alts; alt_i++)
        if (SAME_AS_LUFT_REF (alts[alt_i][0]))
            REJECTIF (true, LO_REF_MULTIALT_SWITCH_SNP, ALTF "Genozip limitation: LUFT(%"PRId64")=%c == ALT[%d], can't switch in multi-allelic SNP", 
                      ALT, opos, oref, alt_i+1);

    // case: oREF is neither REF nor ALT, but sum AF=1 OR AC=AN (i.e. no REF in the samples) - we allow the new REF as it doesn't affect any annotation
    // note: it is sufficient that either the AC or AF indicate that there is no REF in the samples (eg bcftools view --samples only updates AC, not AF)
    // TO DO: extend this to num_alts>1     
    double an;
    if ((num_alts == 1)
    &&    ((vcf_refalt_get_INFO (vb, "AF=") == 1)
       ||  ((an = vcf_refalt_get_INFO (vb, "AN=")) >= 2 && vcf_refalt_get_INFO (vb, "AC=") == an))) {
        vb->new_ref = oref;
        LIFTOK (LO_OK_REF_NEW_SNP, ".\tSNP: REF changed to %c, but lifted because AF=1 or AC=AN", oref);
    }

    REJECTIF (true, LO_NEW_ALLELE_SNP, ALTF "Genozip limitation: LUFT(%"PRId64")=%c is a new SNP allele", ALT, opos, oref);

    return 0; // never reaches here
}

static inline bool vcf_refalt_is_left_anchored (const char *ref, const char **alts, unsigned n_alts)
{
    for (unsigned alt_i=0; alt_i < n_alts; alt_i++)
        if (alts[alt_i][0] != ref[0]) return false; // first base of ALT[alt_i] differs from first base of REF

    return true;
}

// Segging a NON-dual-coordinates file called when genozip --chain
// --chain: set oref and update ostatus
// Analyzes the variant's REF, ALT relative to the Primary and Luft references and the chain file data
// and returns ostatus, possibly also updating the rejects report. It doesn't seg.
LiftOverStatus vcf_refalt_lift (VBlockVCFP vb, const ZipDataLineVCF *dl, bool is_xstrand, WordIndex luft_ref_index,
                                bool *is_left_anchored) // out - only relevant is return value is LO_IS_OK_INDEL()
{
    ASSERT0 (luft_ref_index != WORD_INDEX_NONE, "not expecting luft_ref_index=WORD_INDEX_NONE");

    PosType pos  = dl->pos[0];
    PosType opos = dl->pos[1];
    const unsigned ref_len = vb->main_ref_len;
    const unsigned alt_len = vb->main_alt_len;
    char ref[ref_len], alt[alt_len];

    if (!opos) return LO_OK_REF_SAME_SNP; // POS==oPOS==0 and REF==oREF=='.'

    // note: as we're using external references, the range is always the entire contig, and range->first_pos is always 1.
    const Range *prim_range = ref_seg_get_locked_range (VB, prim_ref, dl->chrom[0], vb->chrom_name, vb->chrom_name_len, pos, 1, WORD_INDEX_NONE, ENT (char, vb->txt_data, vb->line_start), NULL); // doesn't lock as lock=NULL
    ASSVCF (prim_range, "Failed to find PRIM range for chrom=\"%.*s\"", vb->chrom_name_len, vb->chrom_name);

    const Range *luft_range = ref_seg_get_locked_range (VB, gref, dl->chrom[1], NULL, 0, opos, 1, luft_ref_index, ENT (char, vb->txt_data, vb->line_start), NULL);
    ASSVCF (luft_range, "Failed to find LUFT range for chrom=%d", dl->chrom[1]);

    str_toupper_(vb->main_refalt, ref, ref_len);
    str_toupper_(&vb->main_refalt[ref_len + 1], alt, alt_len); // +1 for \t
    
    // split ALT 
    str_split (alt, alt_len, alt_len == 1 ? 1 : 0, ',', alt, false); // short circuit if alt_len=1
    ASSVCF (n_alts, "Invalid ALT=\"%.*s\"", alt_len, alt);
    vb->is_del_sv = false; // initialize

    bool is_single_base_alts = IS_SINGLE_BASE_ALTS (n_alts, alt_len); // All ALTs are a single character - SNPs or deletions
    bool is_snp = (ref_len == 1 && is_single_base_alts);

    // check that the whole REF is mapped in LUFT to the same alignment
    if (ref_len > 1) {
        uint32_t last_ref_aln_i;
        PosType last_pos_in_ref = vb->last_int(VCF_POS) + ref_len - 1;
        LiftOverStatus last_ref_ostatus = vcf_lo_get_liftover_coords (vb, last_pos_in_ref, NULL, NULL, NULL, &last_ref_aln_i); // uses vb->chrom_node_index and POS->last_int

        // Exception: a Deletion that has the anchor base on the alignment, and the rest of REF entirely in the gap after the alignment - 
        // is a REF⇄ALT switch (currently handled only for the non-xstrand case)
        if (last_ref_ostatus == LO_NO_MAPPING_IN_CHAIN && 
            n_alts == 1 && alt_lens[0] == 1 && // a simple Deletion
            SAME_AS_LUFT_REF (ref[0]) && !is_xstrand &&
            chain_get_aln_last_pos (vb->pos_aln_i) == pos && 
            chain_get_aln_gap_after (vb->pos_aln_i) >= (ref_len-1)) 
                return LO_OK_REF_ALT_SWITCH_INDEL;   
            // TO DO: 1. right-anchor on second aligntmnet 2. xstrand - anchor base on opposite alignments

        REJECTIF0 (LO_IS_REJECTED (last_ref_ostatus) || last_ref_aln_i != vb->pos_aln_i,
                   LO_NO_MAPPING_IN_CHAIN, "REF is not fully within a single alignment in the chain file");
    }

    // verify that the REF is consistent between the VCF file and prim_range (if not - there's an error in the VCF or the wrong reference file is used)
    if (!is_same_seq (prim_ref, vb, prim_range, pos, ref, ref_len, false)) {
        char seq[1 + MIN_(ref_len, MAX_BASES_REJECTS_FILE)]; // cap size of automatic variable
        REJECT (LO_REF_MISMATCHES_REFERENCE, ".\tReferenceFile=%s. Possibly this VCF is based on reads that were mapped to a reference other than %s, or the VCF has been erroneously modified",
                ref_dis_subrange (prim_ref, prim_range, pos, sizeof (seq), seq, false), ref_get_filename (prim_ref));
    }
    
    // case: SNP
    if (is_snp) 
        return vcf_refalt_lift_snp (vb, ref, 1, n_alts, alts, alt_lens, STRa(alt), is_xstrand, prim_range, pos, luft_range, opos);

    // case: complex rearrangements (i.e. with [ or ])) (see VCF spec) - Genozip cannot lift these
    REJECTIF (memchr (alt, '[', alt_len) || memchr (alt, ']', alt_len), 
              LO_COMPLEX_REARRANGEMENTS, ALTF "Genozip limitation: Rearrangements with breakends", ALT);

    // case: variant with a symbolic ALT allele (eg <DEL>, <INS> etc)
    if (memchr (alt, '<', alt_len)) {
        vb->is_del_sv = alt_len == 5 && !memcmp (alt, "<DEL>", 5); // used by vcf_seg_INFO_END

        return vcf_refalt_lift_with_sym_allele (vb, STRa(ref), n_alts, alts, alt_lens, STRa(alt), is_xstrand, prim_range, pos, luft_range, opos);
    }

    *is_left_anchored = vcf_refalt_is_left_anchored (ref, alts, n_alts);

    // case: deletion - defined as single-base ALTs - eg GAC G 
    if (ref_len > 1 && is_single_base_alts && *is_left_anchored) 
        return vcf_refalt_lift_deletion (vb, STRa(ref), n_alts, alts, STRa(alt), is_xstrand, prim_range, pos, luft_range, opos); 

    // case insertion - defined as single-base REF, and ALTs that are not all SNPs, eg. "A ACAG"  "A AC,ACA"
    if (ref_len == 1 && !is_single_base_alts && *is_left_anchored)
        return vcf_refalt_lift_insertion (vb, ref, 1, n_alts, alts, alt_lens, STRa(alt), is_xstrand, prim_range, pos, luft_range, opos);
    
    // case: complex variant - anything else 
    return vcf_refalt_lift_complex (vb, STRa(ref), n_alts, alts, alt_lens, STRa(alt), is_xstrand, *is_left_anchored, prim_range, pos, luft_range, opos);
}

// other_REFALT special (other to main REF/ALT fields when segging --chain/primary/luft) snip consists of an integer followed by an optional nucleotide sequence:
// -1 means "new REF" - sequence is new REF
// 0  means - same REF - followed by an optional sequence of reanchoring / left-alignment (left alignment is implemented not implemented yet)
// 1-(MAX_ALLELES-1) - REF⇄ALT switch with these allele. likewise followed by an optional sequence of reanchoring / left-alignment
void vcf_refalt_seg_other_REFALT (VBlockVCFP vb, DidIType did_i, LiftOverStatus ostatus, bool is_xstrand, unsigned add_bytes)
{
    char snip[5] = { SNIP_SPECIAL, VCF_SPECIAL_other_REFALT };
    int snip_len = 2;

    // case: New REF
    if (ostatus == LO_OK_REF_NEW_SNP) {
        snip[2] = '-';  snip[3] = '1'; snip[4] = vb->new_ref;
        snip_len = 5;
    }
    
    // case: Same or Switch
    else {        
        snip[snip_len++] = LO_IS_OK_SWITCH (ostatus) ? '1' : '0'; // currently, --chain only supports switch with the first ALT, but vcf_piz_special_other_REFALT has no such limit

        // add re-anchoring for xstrand INDELS 
        if (is_xstrand && LO_IS_OK_INDEL(ostatus)) 
            snip[snip_len++] = vb->new_ref;
    }

    seg_by_did_i (VB, snip, snip_len, did_i, add_bytes); 
}

// ---------
// PIZ stuff
// ---------

// single-base re-anchoring, eg REF=ACTG ALT=G anchor="T" --> REF=TACT ALT=T  (assuming reference is TACTG)
// left aligning: REF=ACTG ALT=G anchor="TCT" --> REF=TCTA ALT=T (assuming reference is TCTACTG)
static inline void vcf_refalt_shift_right (char *seq, unsigned seq_len, STRp(anchor))
{
    // move what remains of sequence after anchor is inserted
    if (seq_len > anchor_len)
        memcpy (seq + anchor_len, seq, seq_len - anchor_len);
    
    memcpy (seq, anchor, MIN_(seq_len, anchor_len));
}

// actuall reconstruct other REFALT
static void vcf_reconstruct_other_REFALT_do (VBlockP vb, STRp(snip), bool is_xstrand,
                                             const char *other_refalt, unsigned refalt_len, const char *other_refalt_ctx_name)
{
    // short-circuit the most common case of a "same" bi-allelic SNP
    if (refalt_len == 3 && *snip == '0' && snip_len == 1 && !is_xstrand) {
        RECONSTRUCT (other_refalt, refalt_len);
        return;
    }

    char refalt[refalt_len]; // make private copy as ref_alt is in txt_data and will be overwritten
    memcpy (refalt, other_refalt, refalt_len);

    str_split (refalt, refalt_len, 2, '\t', ref_alt, true);
    ASSPIZ (n_ref_alts, "expecting one tab in the %s snip: \"%.*s\"", other_refalt_ctx_name, refalt_len, other_refalt);

    const char *ref = (char *)ref_alts[0]; unsigned ref_len = ref_alt_lens[0];  // assignments for code readability
    const char *alt = (char *)ref_alts[1]; unsigned alt_len = ref_alt_lens[1]; 

    str_split (alt, alt_len, alt_len == 1 ? 1 : 0, ',', alt, false); 

    // case: we need to rev-comp
    if (is_xstrand) {
        str_revcomp ((char*)ref, ref_len);
        for (unsigned alt_i = 0; alt_i < n_alts; alt_i++)
            str_revcomp ((char*)alts[alt_i], alt_lens[alt_i]);
    }

    // parse snip - an integer followed by an optional nucleotide sequence
    char *seq;
    long allele = strtol (snip, &seq, 10);
    unsigned seq_len = &snip[snip_len] - seq;

    // case: we need to reanchor - right-shift and insert new left-anchor for all REF, ALTs 
    if (seq_len && allele >= 0) {
        vcf_refalt_shift_right ((char*)ref, ref_len, seq, seq_len); // snip contains new anchor
        for (unsigned alt_i = 0; alt_i < n_alts; alt_i++)
            vcf_refalt_shift_right ((char *)alts[alt_i], alt_lens[alt_i], seq, seq_len); // alts are pointers into alt
    }

    // case: REF⇄ALT switch. 
    if (allele >= 1) {
        SWAP (ref, alts[allele-1]);
        SWAP (ref_len, alt_lens[allele-1]);
    }

    // case: new REF
    else if (allele == -1) {
        ref = seq;
        ref_len = seq_len;
    }

    RECONSTRUCT_TABBED (ref, ref_len);
    
    for (unsigned alt_i = 0; alt_i < n_alts-1; alt_i++)
        RECONSTRUCT_SEP (alts[alt_i], alt_lens[alt_i], ',');
    
    RECONSTRUCT (alts[n_alts-1], alt_lens[n_alts-1]); // last ALT
}

// Always called for reconstructing the "other" REFALT ("other" coordinates AT THE TIME OF SEGGING) (i.e. oREFALT if we segged a PRIMARY
// and REFALT if we segged a LUFT). This function may be called:
// 1. from a toplevel container as a main VCF field 
// 2. from vcf_piz_special_LIFT_REF when reconstructing the INFO/LIFTXXXX REF field
// Note that the function might be called before or after the main (at seg time) REFALT, depending on whether we are reconstructing in the same coordinates
// as we segged, or the reverse.
SPECIAL_RECONSTRUCTOR (vcf_piz_special_other_REFALT)
{
    if (!reconstruct) return false;

    LiftOverStatus ostatus = last_ostatus;
    ASSPIZ (LO_IS_OK (ostatus), "oStatus=%s unexpected for VCF_SPECIAL_other_REFALT coords=%s ctx=%s", 
            last_ostatus_name_piz, coords_name (CTX (VCF_COORDS)->last_value.i), ctx->tag_name);

    const char *xstrand;
    reconstruct_peek (vb, CTX(VCF_oXSTRAND), &xstrand, 0);
    bool is_xstrand = (*xstrand == 'X'); // process xstrand here, because it will be overwritten in the following reconstruct_peek

    Context *other_refalt_ctx = CTX(ctx->did_i == VCF_REFALT ? VCF_oREFALT : VCF_REFALT);

    const char *other_refalt; 
    unsigned refalt_len;
    reconstruct_peek (vb, other_refalt_ctx, &other_refalt, &refalt_len);
    
    vcf_reconstruct_other_REFALT_do (vb, STRa(snip), is_xstrand, other_refalt, refalt_len, other_refalt_ctx->tag_name);

    return false; // no new value
}

// When segging Luft line: convert REF\tALT from Luft to Primary, in-place, based on ostatus -
// so we can calculate the tie-breaker. Called from vcf_seg_txt_line
void vcf_refalt_seg_convert_to_primary (VBlockVCFP vb, LiftOverStatus ostatus)
{
    // "reconstruct" to vb->txt_data, overwriting current REF\tALT
    uint64_t save_txt_len = vb->txt_data.len;
    vb->txt_data.len = ENTNUM (vb->txt_data, vb->main_refalt);

    vcf_reconstruct_other_REFALT_do (VB, CTX(VCF_REFALT)->last_snip + 2, CTX(VCF_REFALT)->last_snip_len - 2, 
                                     *CTX(VCF_oXSTRAND)->last_snip == 'X', 
                                     vb->main_refalt, vb->main_ref_len+1+vb->main_alt_len, CTX(VCF_oREFALT)->tag_name);

    vb->txt_data.len = save_txt_len; // restore
}

// called for reconstructing the REF field of the INFO/LIFTOVER and INFO/LIFTBACK fields. It reconstructs the full REFALT and then discards the ALT.
SPECIAL_RECONSTRUCTOR (vcf_piz_special_LIFT_REF)
{
    ContextP refalt_ctx = (vb->vb_coords == DC_LUFT ? CTX (VCF_REFALT) : CTX (VCF_oREFALT));

    snip = AFTERENT (char, vb->txt_data);
    reconstruct_from_ctx (vb, refalt_ctx->did_i, 0, true);
    snip_len = (unsigned)(AFTERENT (char, vb->txt_data) - snip);

    const char *after_ref = memchr (snip, '\t', snip_len);
    ASSPIZ (after_ref, "expected a \\t in the %s snip: \"%.*s\" ostatus=%s", 
            refalt_ctx->tag_name, snip_len, snip, last_ostatus_name_piz);

    vb->txt_data.len = ENTNUM (vb->txt_data, reconstruct ? after_ref : snip);

    return false; // no new value
}

// Sometimes called to reconstruct the "main" refalt (main AT THE TIME OF SEGGING), to reconstruct SNPs that were stored relative to a reference.
// This SPECIAL is only used for lines that are bi-allelic SNPs, and either REF or ALT match reference and user compressed with --reference. 
// Not used if compressed with --chain.
SPECIAL_RECONSTRUCTOR (vcf_piz_special_main_REFALT)
{
    if (!reconstruct) goto done;

    ASSPIZ (snip_len==2, "expecting snip_len=2 but seeing %u", snip_len);

    // snip is 3 characters - REF, \t, ALT
    char ref_alt[3] = { 0, '\t', 0 };
    char ref_value = 0;
    
    if (snip[0] == '-' || snip[1] == '-') { 
        PosType pos = CTX (VCF_POS)->last_value.i;

        const Range *range = ref_piz_get_range (vb, gref, pos, 1);
        ASSPIZ (range, "failed to find range for chrom='%s' pos=%"PRId64" in %s", vb->chrom_name, pos, ref_get_filename(gref));
        
        uint32_t idx = pos - range->first_pos;
        ASSPIZ (ref_is_nucleotide_set (range, idx), "reference is not set: chrom=%.*s pos=%"PRId64, range->chrom_name_len, range->chrom_name, pos);
        ref_value = ref_base_by_idx (range, idx);
    }

    // recover ref
    if (snip[0] == '-') 
        ref_alt[0] = ref_value;
    else 
        ref_alt[0] = snip[0];

    // recover alt
    if (snip[1] == '+') { // the alt has the most common value for a SNP
        if      (ref_alt[0] == 'A') ref_alt[2] = 'G';
        else if (ref_alt[0] == 'C') ref_alt[2] = 'T';
        else if (ref_alt[0] == 'G') ref_alt[2] = 'A';
        else if (ref_alt[0] == 'T') ref_alt[2] = 'C';
    }
    else if (snip[1] == '-')  // the alt has the reference value
        ref_alt[2] = ref_value;

    else // the alt is specified verbatim
        ref_alt[2] = snip[1];

    RECONSTRUCT (ref_alt, sizeof (ref_alt));

done:
    return false; // no new value
}   

// --snps-only implementation (called from vcf_piz_container_cb)
bool vcf_refalt_piz_is_variant_snp (VBlockP vb)
{
    DidIType refalt = vb->vb_coords == DC_PRIMARY ? VCF_REFALT : VCF_oREFALT;
    unsigned txt_len = vb->last_txt_len (refalt);
    if (txt_len <= 3) return true; // short circuit most common case  of a bi-allelic SNP (<3 can never happen, here to avoid issues)

    const char *txt = last_txt (vb, refalt);
    return txt[1] == '\t' && str_count_char (&txt[2], txt_len-2, ',') * 2 == txt_len-3; // true if multi-allelic SNP
}

// --snps-only implementation (called from vcf_piz_container_cb)
bool vcf_refalt_piz_is_variant_indel (VBlockP vb)
{
    DidIType refalt = vb->vb_coords == DC_PRIMARY ? VCF_REFALT : VCF_oREFALT;
    unsigned txt_len = vb->last_txt_len (refalt);
    if (txt_len == 3 || vcf_refalt_piz_is_variant_snp (vb)) return false; // Not an INDEL: its a SNP (short circuit most common case of a bi-allelic SNP)

    const char *txt = last_txt (vb, refalt);
    if (memchr (txt, '<', txt_len) || memchr (txt, '[', txt_len) || memchr (txt, ']', txt_len)) return false; // Not an INDEL: its an SV

    return true;
}
