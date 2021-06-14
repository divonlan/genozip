// ------------------------------------------------------------------
//   vcf_refalt.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

// ---------
// Seg stuff
// ---------

// optimize REF and ALT, for simple one-character REF/ALT (i.e. mostly a SNP or no-variant)
static void vcf_refalt_seg_ref_alt_snp (VBlockVCFP vb, char main_ref, char main_alt)
{
    char new_ref=0, new_alt=0;

    // if we have a reference, we use it (we treat the reference as PRIMARY)
    if ((flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) && vb->line_coords == DC_PRIMARY) {
        PosType pos = vb->last_int(VCF_POS);

        RefLock lock;
        Range *range = ref_seg_get_locked_range ((VBlockP)vb, gref, vb->chrom_node_index, vb->chrom_name, vb->chrom_name_len, pos, 1, NULL, &lock);
        uint32_t index_within_range = pos - range->first_pos;

        ref_assert_nucleotide_available (range, pos);
        char ref = ref_base_by_idx (range, index_within_range);

        if (main_ref == ref) new_ref = '-'; // this should always be the case...
        if (main_alt == ref) new_alt = '-'; 

        if (flag.reference == REF_EXT_STORE)
            bit_array_set (&range->is_set, index_within_range);

        ref_unlock (gref, lock);
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

        seg_by_did_i (vb, refalt_special, sizeof(refalt_special), SEL (VCF_REFALT, VCF_oREFALT), 0); // we do the account in vcf_refalt_seg_main_ref_alt
    }
    // if not - just the normal snip
    else {
        char refalt_normal[3] = { 0, '\t', 0 };
        refalt_normal[0] = main_ref;
        refalt_normal[2] = main_alt;

        seg_by_did_i (vb, refalt_normal, sizeof(refalt_normal), SEL (VCF_REFALT, VCF_oREFALT), 0); // we do the account in vcf_refalt_seg_main_ref_alt
    }
}

void vcf_refalt_seg_main_ref_alt (VBlockVCFP vb, const char *ref, unsigned ref_len, const char *alt, unsigned alt_len)
{
    // optimize ref/alt in the common case of single-character
    if (ref_len == 1 && alt_len == 1) 
        vcf_refalt_seg_ref_alt_snp (vb, *ref, *alt);

    else {
        char ref_alt[ref_len + 1 + alt_len];
        memcpy (ref_alt, ref, ref_len);
        ref_alt[ref_len] = '\t';
        memcpy (&ref_alt[ref_len+1], alt, alt_len);

        seg_by_did_i (vb, ref_alt, ref_len + alt_len + 1, SEL (VCF_REFALT, VCF_oREFALT), 0);
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

#define ECLIPSE(len) ((len) > 100 ? "..." : "")
#define ECLIPSED(seq,len) MIN (100,(len)), seq, ECLIPSE(len)
#define XSTRAND (is_xstrand ? ". XSTRAND=X" : "")

#define ALTF "ALT=%.*s%s "
#define ALT ECLIPSED (alt, alt_len)

#define DEF_PRIM(len) char prim_str[MIN (101, (len)+1)]; unsigned prim_str_len = (unsigned)(len)
#define PRIMF "PRIM=%s%s "
#define PRIM ref_dis_subrange (prim_range, pos, sizeof (prim_str), prim_str, false), ECLIPSE (prim_str_len)

#define DEF_LUFT(len) char luft_str[MIN (101, (len)+1)]; unsigned luft_str_len = (unsigned)(len)
#define LUFTF "LUFT(%"PRId64")=%s%s "
#define LUFT opos, ref_dis_subrange (luft_range, opos, sizeof (luft_str), luft_str, is_xstrand), ECLIPSE (luft_str_len)

#define SAME(base1,base2) (is_xstrand ? REVCOMP[(int)base1]==(base2) : (base1)==(base2))

// false is not the same or if either goes beyond the end of the range 
static inline bool vcf_refalt_lift_is_same (const Range *range, PosType pos, const char *seq, PosType seq_len, bool is_xstrand) 
{
    if (!is_xstrand) {
        if (pos + seq_len - 1 > range->last_pos) return false; // seq goes beyond the end of range

        for (PosType i=0; i < seq_len ; i++) 
            if (ref_base_by_pos (range, pos + i) != UPPER_CASE (seq[i])) return false;
    }
    else { // xstrand
        if (pos < seq_len) return false; // seq goes beyond the start of range

        for (PosType i=0; i < seq_len ; i++)
            if (ref_base_by_pos (range, pos - i) != UPPER_REVCOMP[(int)seq[i]]) return false;
    }

    return true;
}

// counts repeated bases (forwards or backwards depending on direction)
typedef enum { REP_NO_REVCOMP, REP_REVCOMP} RepIsRevComp;
typedef enum { REP_RIGHT_ALIGNED, REP_LEFT_ALIGNED } RepAlignment;
static unsigned vcf_refalt_lift_get_repeats (const Range *range, const char *rep, PosType rep_len, PosType pos, 
                                             RepIsRevComp revcomp, RepAlignment alignment)
{
    unsigned num_reps = 0;
    if (alignment == REP_LEFT_ALIGNED) {
        if (revcomp == REP_NO_REVCOMP) {
            for (PosType p=pos; p <= range->last_pos; p++, num_reps++) 
                if (ref_base_by_pos (range, p) != rep[num_reps % rep_len]) break;
        }
        else {
            for (PosType p=pos; p >= 1; p--, num_reps++) 
                if (ref_base_by_pos (range, p) != REVCOMP[(int)rep[num_reps % rep_len]]) break;
        }
    }
    else {
        char reverse_rep[rep_len];
        str_reverse (rep, reverse_rep, rep_len);

        if (revcomp == REP_NO_REVCOMP) {
            for (PosType p=pos; p >= 1; p--, num_reps++) 
                if (ref_base_by_pos (range, p) != reverse_rep[num_reps % rep_len]) break;
        }
        else {
            for (PosType p=pos; p <= range->last_pos; p++, num_reps++) 
                if (ref_base_by_pos (range, p) != REVCOMP[(int)reverse_rep[num_reps % rep_len]]) break;
        }
    }    
    return num_reps;
}

// when --chain: handle a Deletion
static LiftOverStatus vcf_refalt_lift_deletion (VBlockVCFP vb, const char *ref, PosType ref_len,
                                                unsigned num_alts, const char **alts, const char *alt, unsigned alt_len, 
                                                bool is_xstrand,
                                                const Range *prim_range, PosType pos, 
                                                const Range *luft_range, PosType opos)
                                                
{
    PosType rep_len = ref_len-1;

    // if LUFT and PRIMARY is not the same for ref_len bases: example: "GAC G" Primary: "GAC" 
    // Luft=GTT: REF<>ALT switch: "G GAC"
    // Luft=TAC or TTT: too complicated change. 
    // note: this works also for deletions that are not left or right aligned eg "GAC A" or even "GC A"
    if (!vcf_refalt_lift_is_same (luft_range, opos, ref, ref_len, is_xstrand)) {
        REJECTIF0 (num_alts > 1, LO_NEW_ALLELE_INDEL, "Genozip limitation: REF changed in a DEL variant that has more than one ALT");
        REJECTIF0 (!vcf_refalt_lift_is_same (luft_range, opos, alts[0], 1, is_xstrand), LO_NEW_ALLELE_INDEL, "Genozip limitation: REF changed in DEL variant, but not to ALT");
        return LO_OK_REF_ALT_SWITCH_INDEL; // simple REF<>ALT switch - deletion already included in Luft
    }

    // now we know that REF is stringwise the same as the reference, but we need to test for repeats

    unsigned num_rep_prim, num_rep_luft;
    bool all_same=true;
    int switch_alt_i=-1;

    for (unsigned alt_i=0; alt_i < num_alts; alt_i++) {

        bool left_aligned  = ref[0] == alts[alt_i][0];
        bool right_aligned = ref[ref_len-1] == alts[alt_i][0];

        num_rep_prim = num_rep_luft = rep_len; // we already verified that LUFT and PRIMARY are the same - those count as repeats

        // case: none-aligned, but REF==oREF, eg "GAC A" or even "GC A" (no repeats issue in this case) - we already tested that LUFT contains REF
        if (!left_aligned && !right_aligned) continue;

        // case: left-aligned eg "GAC G" - count right repeats ('G' being the "achor base" and 'AC' is the repeat)
        // case: right-aligned eg "ACT T" - count right repeats:
        // note: this logic also covers a case of "both-aligned": GAG G
        // note: 
        // in both cases:
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACACT   (4 repeats)  -> REF/ALT Unchanged
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACT     (3 repeats)  -> REF <> ALT switch. Luft: G GAC
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACACACT (5+ repeats) -> new allele: GACACAC G,GAC (old REF was a 4-repeats: now GAC   ; old ALT was a 3-repeats: now G)
        // PRIM: GACACACACT (4 repeats) LUFT: GACACT       (2- repeats) -> new allele: G GAC,GACAC   (old REF was a 4-repeats: now GACAC ; old ALT was a 3-repeats: now GAC)
        
        if (left_aligned) {
            const char *rep = ref + 1;
            num_rep_prim +=               vcf_refalt_lift_get_repeats (prim_range, rep, rep_len, pos  + ref_len, REP_NO_REVCOMP, REP_LEFT_ALIGNED); // count EXCLUDING this the present repeat (already set to 1)
            num_rep_luft += !is_xstrand ? vcf_refalt_lift_get_repeats (luft_range, rep, rep_len, opos + ref_len, REP_NO_REVCOMP, REP_LEFT_ALIGNED)
                                        : vcf_refalt_lift_get_repeats (luft_range, rep, rep_len, opos - ref_len, REP_REVCOMP,    REP_LEFT_ALIGNED);
        }

        if (right_aligned) {
            const char *rep = ref;
            num_rep_prim +=               vcf_refalt_lift_get_repeats (prim_range, rep, rep_len, pos  - rep_len, REP_NO_REVCOMP, REP_RIGHT_ALIGNED);
            num_rep_luft += !is_xstrand ? vcf_refalt_lift_get_repeats (luft_range, rep, rep_len, opos - rep_len, REP_NO_REVCOMP, REP_RIGHT_ALIGNED)
                                        : vcf_refalt_lift_get_repeats (luft_range, rep, rep_len, opos + ref_len, REP_REVCOMP, REP_RIGHT_ALIGNED);
        }

        if (num_rep_luft != num_rep_prim) all_same = false; // at least one of the ALTs is not REF_SAME
        if (num_rep_luft == num_rep_prim - rep_len) switch_alt_i = alt_i; // found an ALT that is a REF<>ALT switch
    }


    // we can lift a Insertion IF: either - all ALTs agree that it is REF_SAME ; or - single ALT and REF<>ALT switch
    if (all_same) return LO_OK_REF_SAME_INDEL;
    if (switch_alt_i==0 && num_alts == 1) return LO_OK_REF_ALT_SWITCH_INDEL;

    DEF_PRIM (num_rep_prim + 2);
    DEF_LUFT (num_rep_luft + 2);
    
    // encountered a switch ALT in a multi allelic DEL variant
    REJECTIF (switch_alt_i >= 0, LO_REF_MULTIALT_SWITCH_INDEL, ALTF "Genozip limitation: " PRIMF " and ALT[%d] matches " LUFTF ", but can't switch in multi-allelic DELETION variant", 
              ALT, PRIM, switch_alt_i+1, LUFT);
    
    // LUFT is a new allele
    REJECTIF (true, LO_NEW_ALLELE_INDEL, ALTF "Genozip limitation: new INDEL allele: " PRIMF LUFTF, ALT, PRIM, LUFT); 

    return 0; // never reaches here
}

// when --chain: handle an Insertion
static LiftOverStatus vcf_refalt_lift_insertion (VBlockVCFP vb, const char *ref, unsigned ref_len, /* must be 1 */ 
                                                 unsigned num_alts, const char **alts, const unsigned *alt_lens, const char *alt, unsigned alt_len, 
                                                 bool is_xstrand,
                                                 const Range *prim_range, PosType pos, 
                                                 const Range *luft_range, PosType opos)
                                                
{    
    unsigned num_rep_prim=0, num_rep_luft=0; 

    bool all_same=true;
    int switch_alt_i=-1;
    for (unsigned alt_i=0; alt_i < num_alts; alt_i++) {

        bool left_aligned  = (ref[0] == alts[alt_i][0]);
        bool right_aligned = (ref[0] == alts[alt_i][alt_lens[alt_i]-1]);

        num_rep_prim = num_rep_luft = 0; 

        // case: none-aligned, but REF==oREF, eg "A GCC"
        if (!left_aligned && !right_aligned) {
            if (vcf_refalt_lift_is_same (luft_range, opos, alts[alt_i], alt_lens[alt_i], is_xstrand))
                switch_alt_i = alt_i;
            continue;
        }

        // check if anchor remains - eg "A AGC" (left aligned) or "A GCA" (right aligned) (anchor ia A)
        if (!SAME (ref_base_by_pos (luft_range, opos), ref[0])) {
            all_same = false;
            continue;
        }

        // case: left-aligned eg "GAC G" - count right repeats ('G' being the "achor base" and 'AC' is the repeat)
        // case: right-aligned eg "ACT T" - count right repeats:
        // note: this logic also covers a case of "both-aligned": GAG G
        // note: 
        // in both cases:
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACACT   (4 repeats)  -> REF/ALT Unchanged
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACT     (3 repeats)  -> REF <> ALT switch. Luft: G GAC
        // PRIM: GACACACACT (4 repeats) LUFT: GACACACACACT (5+ repeats) -> new allele: GACACAC G,GAC (old REF was a 4-repeats: now GAC   ; old ALT was a 3-repeats: now G)
        // PRIM: GACACACACT (4 repeats) LUFT: GACACT       (2- repeats) -> new allele: G GAC,GACAC   (old REF was a 4-repeats: now GACAC ; old ALT was a 3-repeats: now GAC)

        PosType rep_len = alt_lens[alt_i]-1;

        if (left_aligned) {
            const char *rep = alts[alt_i] + 1;
            num_rep_prim +=               vcf_refalt_lift_get_repeats (prim_range, rep, rep_len, pos  + 1, REP_NO_REVCOMP, REP_LEFT_ALIGNED); 
            num_rep_luft += !is_xstrand ? vcf_refalt_lift_get_repeats (luft_range, rep, rep_len, opos + 1, REP_NO_REVCOMP, REP_LEFT_ALIGNED)
                                        : vcf_refalt_lift_get_repeats (luft_range, rep, rep_len, opos - 1, REP_REVCOMP,    REP_LEFT_ALIGNED);
        }

        if (right_aligned) {
            const char *rep = alts[alt_i];
            num_rep_prim +=               vcf_refalt_lift_get_repeats (prim_range, rep, rep_len, pos  - 1, REP_NO_REVCOMP, REP_RIGHT_ALIGNED);
            num_rep_luft += !is_xstrand ? vcf_refalt_lift_get_repeats (luft_range, rep, rep_len, opos - 1, REP_NO_REVCOMP, REP_RIGHT_ALIGNED)
                                        : vcf_refalt_lift_get_repeats (luft_range, rep, rep_len, opos + 1, REP_REVCOMP,    REP_RIGHT_ALIGNED);
        }
    
        if (num_rep_luft != num_rep_prim) all_same = false; // at least one of the ALTs is not REF_SAME
        if (num_rep_luft == num_rep_prim + rep_len) switch_alt_i = alt_i; // found an ALT that is a REF<>ALT switch
    }

    // we can lift a Insertion IF: either - all ALTs agree that it is REF_SAME ; or - single ALT and REF<>ALT switch
    if (all_same) return LO_OK_REF_SAME_INDEL;
    if (switch_alt_i==0 && num_alts == 1) return LO_OK_REF_ALT_SWITCH_INDEL;

    DEF_PRIM (num_rep_prim + 2); // +2 for the bases before and after the repeats
    DEF_LUFT (num_rep_luft + 2);

    // encountered a switch ALT in a multi allelic INS variant
    REJECTIF (switch_alt_i >= 0, LO_REF_MULTIALT_SWITCH_INDEL, ALTF "Genozip limitation: " PRIMF " and ALT[%d] matches " LUFTF ", but can't switch in multi-allelic INSERTION variant", 
              ALT, PRIM, switch_alt_i+1, LUFT);
    
    // LUFT is a new allele
    REJECTIF (true, LO_NEW_ALLELE_INDEL, ALTF "Genozip limitation: new INDEL allele: " PRIMF LUFTF, ALT, PRIM, LUFT);

    return 0; // never reaches here
}

// when --chain: handle an structural variant with an ALT containing a <***> ID  (eg <DEL>)
static LiftOverStatus vcf_refalt_lift_sv (VBlockVCFP vb, const char *ref, unsigned ref_len, 
                                          unsigned num_alts, const char **alts, const unsigned *alt_lens, const char *alt, unsigned alt_len,
                                          bool is_xstrand,
                                          const Range *prim_range, PosType pos, 
                                          const Range *luft_range, PosType opos)
{
    // case: structural variant - but REF unchanged 
    if (vcf_refalt_lift_is_same (luft_range, opos, ref, ref_len, is_xstrand))
        return LO_OK_REF_SAME_SV;

    else
        REJECTIF (true, LO_COMPLEX_REFALT, ALTF "Genozip limitation: Structural variant and REF changed", ALT);

    return 0; // never reaches here
}

// when --chain: handle a complex variant - REF and at least one of the ALTs is longer than 1 base
static LiftOverStatus vcf_refalt_lift_complex (VBlockVCFP vb, const char *ref, unsigned ref_len, 
                                               unsigned num_alts, const char **alts, const unsigned *alt_lens, const char *alt, unsigned alt_len,
                                               bool is_xstrand,
                                               const Range *prim_range, PosType pos, 
                                               const Range *luft_range, PosType opos)
{
    // case: complex variant - but REF unchanged 
    if (vcf_refalt_lift_is_same (luft_range, opos, ref, ref_len, is_xstrand))
        return LO_OK_REF_SAME_INDEL;

    // case: complex variant - REF<>ALT switch - supported only if bi-allelic and not SV
    if (num_alts==1 && vcf_refalt_lift_is_same (luft_range, opos, alts[0], alt_lens[0], is_xstrand)) 
        return LO_OK_REF_ALT_SWITCH_INDEL;

    for (unsigned alt_i=0; alt_i < num_alts; alt_i++)
        if (vcf_refalt_lift_is_same (luft_range, opos, alts[alt_i], alt_lens[alt_i], is_xstrand)) {
            DEF_LUFT (alt_lens[alt_i]);
            REJECTIF (true, LO_REF_MULTIALT_SWITCH_INDEL, ALTF "Genozip limitation: " LUFTF "== ALT[%d], but can't switch in multi-allelic variant", 
                      ALT, LUFT, alt_i+1);
        }

    DEF_LUFT (ref_len);
    REJECTIF (true, LO_NEW_ALLELE_INDEL, ALTF "Genozip limitation: " LUFTF "is a new allele", ALT, LUFT);

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
static LiftOverStatus vcf_refalt_lift_snp (VBlockVCFP vb, const char *ref, unsigned ref_len/* must be 1 */, 
                                           unsigned num_alts, const char **alts, const unsigned *alt_lens, const char *alt, unsigned alt_len,
                                           bool is_xstrand,
                                           const Range *prim_range, PosType pos, 
                                           const Range *luft_range, PosType opos)
{
    char oref = ref_base_by_pos (luft_range, opos); // range is 0-based

    if (SAME (oref, ref[0])) 
        return LO_OK_REF_SAME_SNP;

    if (num_alts == 1 && SAME (oref, alts[0][0])) 
        return LO_OK_REF_ALT_SWITCH_SNP;

    for (unsigned alt_i=0; alt_i < num_alts; alt_i++)
        if (SAME (oref, alts[alt_i][0]))
            REJECTIF (true, LO_REF_MULTIALT_SWITCH_SNP, ALTF "Genozip limitation: LUFT(%"PRId64")=%c == ALT[%d], can't switch in multi-allelic SNP", 
                      ALT, opos, oref, alt_i+1);

    // case: oREF is neither REF nor ALT, but sum AF=1 OR AC=AN (i.e. no REF in the samples) - we allow the new REF as it doesn't affect any annotation
    // note: it is sufficient that either the AC or AF indicate that there is no REF in the samples (eg bcftools view --samples only updates AC, not AF)
    // TO DO: extend this to num_alts>1     
    double an;
    if ((num_alts == 1)
    &&    ((vcf_refalt_get_INFO (vb, "AF=") == 1)
       ||  ((an = vcf_refalt_get_INFO (vb, "AN=")) >= 2 && vcf_refalt_get_INFO (vb, "AC=") == an))) {
        vb->new_ref = ref[0];
        return LO_OK_REF_NEW_SNP;
    }

    REJECTIF (true, LO_NEW_ALLELE_SNP, ALTF "Genozip limitation: LUFT(%"PRId64")=%c is a new SNP allele", ALT, opos, oref);

    return 0; // never reaches here
}

// Segging a NON-dual-coordinates file called when genozip --chain
// --chain: set oref and update ostatus
// Analyzes the variant's REF, ALT relative to the Primary and Luft references and the chain file data
// and returns ostatus, possibly also updating the rejects report. It doesn't seg.
LiftOverStatus vcf_refalt_lift (VBlockVCFP vb, const ZipDataLineVCF *dl, bool is_xstrand)
{
    const PosType pos      = dl->pos[0];
    const PosType opos     = dl->pos[1];
    const unsigned ref_len = vb->main_ref_len;
    const unsigned alt_len = vb->main_alt_len;
    char ref[ref_len], alt[alt_len];
    
    if (!opos) return LO_OK_REF_SAME_SNP; // POS==oPOS==0 and REF==oREF=='.'

    // note: as we're using external references, the range is always the entire contig, and range->first_pos is always 1.
    const Range *prim_range = ref_seg_get_locked_range ((VBlockP)vb, prim_ref, dl->chrom_index[0], vb->chrom_name, vb->chrom_name_len, pos, 1, ENT (char, vb->txt_data, vb->line_start), NULL); // doesn't lock as lock=NULL
    ASSVCF (prim_range, "Failed to find PRIM range for chrom=\"%.*s\"", vb->chrom_name_len, vb->chrom_name);

    const Range *luft_range = ref_seg_get_locked_range ((VBlockP)vb, gref, dl->chrom_index[1], NULL, 0, opos, 1, ENT (char, vb->txt_data, vb->line_start), NULL);
    ASSVCF (luft_range, "Failed to find LUFT range for chrom=%d", dl->chrom_index[1]);

    str_toupper_(vb->main_refalt, ref, ref_len);
    str_toupper_(&vb->main_refalt[ref_len + 1], alt, alt_len); // +1 for \t
    
    // split ALT 
    str_split (alt, alt_len, alt_len == 1 ? 1 : 0, ',', alt, false); // short circuit if alt_len=1
    ASSVCF (n_alts, "Invalid ALT=\"%.*s\"", alt_len, alt);
    vb->is_del_sv = false; // initialize

    bool is_single_base_alts = (n_alts*2 -1 == alt_len); // All ALTs are a single character - SNPs or deletions

    // check that the whole REF is mapped in LUFT to the same alignment
    if (ref_len > 1) {
        uint32_t last_ref_aln_i;
        PosType last_pos_in_ref = vb->last_int(VCF_POS) + ref_len - 1;
        LiftOverStatus last_ref_ostatus = vcf_lo_get_liftover_coords (vb, last_pos_in_ref, NULL, NULL, NULL, &last_ref_aln_i); // uses vb->chrom_node_index and POS->last_int

        // Exception: a Deletion that is entirely in the gap after the alignment - is a REF<>ALT switch
        if (last_ref_ostatus == LO_NO_MAPPING_IN_CHAIN && 
            n_alts == 1 && alt_lens[0] == 1 && // a simple Deletion
            SAME (ref_base_by_pos (luft_range, opos), ref[0]) &&
            chain_get_aln_last_pos (vb->pos_aln_i) == pos && 
            chain_get_aln_gap_after (vb->pos_aln_i) >= (ref_len-1)) 
                return LO_OK_REF_ALT_SWITCH_INDEL;

        REJECTIF0 (LO_IS_REJECTED (last_ref_ostatus) || last_ref_aln_i != vb->pos_aln_i,
                   LO_NO_MAPPING_IN_CHAIN, "REF is not fully within a single alignment in the chain file");
    }

    // verify that the REF is consistent between the VCF file and prim_range (if not - there's an error in the VCF or the wrong reference file is used)
    if (!vcf_refalt_lift_is_same (prim_range, pos, ref, ref_len, false)) {
        char seq[MIN (ref_len, 100)]; // cap size of automatic variable
        REJECT (LO_REF_MISMATCHES_REFERENCE, "doesn't match the reference file=\"%s\" - either error in the VCF or wrong reference file",
                ref_dis_subrange (prim_range, pos, MIN (ref_len, 100), seq, false));
    }
    
    // case: SNP
    if (ref_len == 1 && is_single_base_alts) 
        return vcf_refalt_lift_snp (vb, ref, 1, n_alts, alts, alt_lens, alt, alt_len, is_xstrand, prim_range, pos, luft_range, opos);

    // case: complex re-arrangements
    else if (memchr (alt, '[', alt_len) || memchr (alt, ']', alt_len)) // VCF 4.3 complex rearrangements (i.e. with [ or ]))
        REJECTIF0 (true, LO_COMPLEX_REARRANGEMENTS, "Genozip limitation: Rearrangements with breakends");

    // case: structural variant - at leaset one of the ALTs eith eg <DEL>
    else if (memchr (alt, '<', alt_len)) {
        vb->is_del_sv = alt_len == 5 && !memcmp (alt, "<DEL>", 5); // used by vcf_seg_INFO_END

        return vcf_refalt_lift_sv (vb, ref, ref_len, n_alts, alts, alt_lens, alt, alt_len, is_xstrand, prim_range, pos, luft_range, opos);
    }

    // case: complex variant - REF *and* at least one of the ALTs are multi-base
    else if (ref_len > 1 && !is_single_base_alts)
        return vcf_refalt_lift_complex (vb, ref, ref_len, n_alts, alts, alt_lens, alt, alt_len, is_xstrand, prim_range, pos, luft_range, opos);

    // case: deletion - defined as single-base ALTs - eg GAC G / GAC C / GAC A / GAC T / GAC A,G
    else if (ref_len > 1 && is_single_base_alts) 
        return vcf_refalt_lift_deletion (vb, ref, ref_len, n_alts, alts, alt, alt_len, is_xstrand, prim_range, pos, luft_range, opos); 

    // case insertion - defined as single-base REF, and ALTs that are not all SNPs, eg. "A ACAG"  "A AC,ACA" "A GT,TC"
    else if (ref_len == 1 && !is_single_base_alts)
        return vcf_refalt_lift_insertion (vb, ref, 1, n_alts, alts, alt_lens, alt, alt_len, is_xstrand, prim_range, pos, luft_range, opos);
    
    ASSVCF (false, "unidentified variant type REF=%.*s ALT=%.*s", ref_len, ref, alt_len, alt);

    return 0; // never reaches here
}

// ---------
// PIZ stuff
// ---------

// Always called for reconstructing the "other" REFALT ("other" coordinates AT THE TIME OF SEGGING) (i.e. oREFALT if we segged a PRIMARY
// and REFALT if we segged a LUFT). This function may be called:
// 1. from a toplevel container as a main VCF field 
// 2. from vcf_piz_special_LIFT_REF when reconstructing the INFO/LIFTXXXX REF field
// Note that the function might be called before or after the main (at seg time) REFALT, depending on whether we are reconstructing in the same coordinates
// as we segged, or the reverse.
SPECIAL_RECONSTRUCTOR (vcf_piz_special_other_REFALT)
{
    if (!reconstruct) return false;

    Context *other_refalt_ctx = CTX(ctx->did_i == VCF_REFALT ? VCF_oREFALT : VCF_REFALT);

    const char *other_refalt; 
    unsigned other_refalt_len;
    reconstruct_peek (vb, other_refalt_ctx, &other_refalt, &other_refalt_len);
    
    str_split (other_refalt, other_refalt_len, 2, '\t', ref_alt, true);
    ASSPIZ (n_ref_alts, "expecting one tab in the %s snip: \"%.*s\"", other_refalt_ctx->name, other_refalt_len, other_refalt);

    // make private copy as ref_alt is in txt_data and will be overwritten
    unsigned ref_len = ref_alt_lens[0]; char ref[ref_len]; memcpy (ref, ref_alts[0], ref_len); 
    unsigned alt_len = ref_alt_lens[1]; char alt[alt_len]; memcpy (alt, ref_alts[1], alt_len);

    const char *xstrand;
    reconstruct_peek (vb, CTX(VCF_oXSTRAND), &xstrand, 0);
    bool is_xstrand = (*xstrand == 'X');    
    char ref_revcomp[ref_len], alt_revcomp[alt_len];

    // case: REF<>ALT switch. Genozip only supports this for bi-allelic SNPs. REF and ALT can be of any length, possibly with xstrand.
    LiftOverStatus ostatus = last_ostatus;
    if (ostatus == LO_OK_REF_ALT_SWITCH_SNP || ostatus == LO_OK_REF_ALT_SWITCH_INDEL) { 
        RECONSTRUCT_SEP ((is_xstrand ? str_to_revcomp (alt, alt_revcomp, alt_len) : alt), alt_len, '\t');
        RECONSTRUCT ((is_xstrand ? str_to_revcomp (ref, ref_revcomp, ref_len) : ref), ref_len);
    }

    else if (ostatus == LO_OK_REF_NEW_SNP) {
        RECONSTRUCT_SEP (snip, snip_len, '\t'); // REF was provided explicitly
        RECONSTRUCT ((is_xstrand ? str_to_revcomp (alt, alt_revcomp, alt_len) : alt), alt_len);
    }
    
    else if (LO_IS_OK (ostatus)) {
        RECONSTRUCT_SEP ((is_xstrand ? str_to_revcomp (ref, ref_revcomp, ref_len) : ref), ref_len, '\t');
        RECONSTRUCT ((is_xstrand ? str_to_revcomp (alt, alt_revcomp, alt_len) : alt), alt_len);
    }
    
    else
        ASSPIZ (0, "oStatus=%s (coords=%u) unexpected for VCF_SPECIAL_other_REFALT", last_ostatus_name_piz, (unsigned)vb->last_index (VCF_COORDS));

    return false; // no new value
}

// called for reconstructing the REF field of the INFO/LIFTOVER and INFO/LIFTBACK fields. It reconstructs the full REFALT and then discards the ALT.
SPECIAL_RECONSTRUCTOR (vcf_piz_special_LIFT_REF)
{
    DidIType refalt_did_i = flag.luft ? VCF_REFALT : VCF_oREFALT;

    snip = AFTERENT (char, vb->txt_data);
    reconstruct_from_ctx (vb, refalt_did_i, 0, true);
    snip_len = (unsigned)(AFTERENT (char, vb->txt_data) - snip);

    const char *after_ref = memchr (snip, '\t', snip_len);
    ASSPIZ (after_ref, "expected a \\t in the %s snip: \"%.*s\" ostatus=%s", 
            CTX(refalt_did_i)->name, snip_len, snip, last_ostatus_name_piz);

    vb->txt_data.len = ENTNUM (vb->txt_data, after_ref);

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
        PosType pos = vb->last_int (VCF_POS);

        const Range *range = ref_piz_get_range (vb, gref, pos, 1);
        ASSPIZ (range, "failed to find range for chrom='%s' pos=%"PRId64, vb->chrom_name, pos);
        
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
