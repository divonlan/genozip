// ------------------------------------------------------------------
//   ref_iupacs.c
//   Copyright (C) 2021-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "data_types.h"
#include "reference.h"
#include "vblock.h"
#include "ref_private.h"
#include "context.h"
#include "file.h"
#include "endianness.h"
#include "zfile.h"
#include "buffer.h"

//----------------------
// make-reference side
//----------------------

typedef struct {
    uint64_t idx;
    VBIType vblock_i;
    char iupac;
} MakeIupac;

// --make-reference: called by compute thread, when iupac is found
void ref_iupacs_add_do (VBlockP vb, uint64_t idx, char iupac)
{
    buf_alloc (vb, &CTX(FASTA_NONREF_X)->local, 1, 100, MakeIupac, 2, "contexts->local");
    BNXT (MakeIupac, CTX(FASTA_NONREF_X)->local) = (MakeIupac){ .vblock_i = vb->vblock_i, .idx = idx, .iupac = iupac };
}

// --make-reference: called by main thread, after each FASTA vb, in order of VBs which is also the order of ranges
void ref_iupacs_after_compute (VBlockP vb)
{
    if (CTX(FASTA_NONREF_X)->local.len) {
        buf_add_buf (evb, &gref->iupacs_buf, &CTX(FASTA_NONREF_X)->local, MakeIupac, "iupacs_buf");
        buf_free (CTX(FASTA_NONREF_X)->local);
    }
}

// --make-reference: create a REF_IUPACS section if there are any IUPACs in this reference. The section consists of pairs (GPOS, IUPAC) sorted by GPOS
void ref_iupacs_compress (void)
{
    ARRAY (MakeIupac, make_iupacs, gref->iupacs_buf);
    if (!make_iupacs_len) return;

    ASSERTNOTINUSE (evb->scratch);
    buf_alloc (evb, &evb->scratch, 0, make_iupacs_len, Iupac, 0, "scratch");
    ARRAY (Iupac, iupacs, evb->scratch);

    if (flag.show_ref_iupacs) iprintf ("\nIUPACs found in %s:\n", z_name);

    PosType64 last_gpos = 0;
    for (uint64_t i=0; i < make_iupacs_len; i++) {
        const Range *r = B(Range, gref->ranges, make_iupacs[i].vblock_i-1);

        PosType64 gpos = r->gpos + make_iupacs[i].idx;

        if (flag.show_ref_iupacs)   
            iprintf ("IUPAC=%c\nCHROM=%s\tPOS=%"PRIu64"\tGPOS=%"PRIu64"\n", 
                     make_iupacs[i].iupac, ctx_snip_from_zf_nodes (ZCTX(FASTA_CONTIG), r->chrom, 0, 0), r->first_pos + make_iupacs[i].idx, gpos);

        iupacs[i] = (Iupac) { .gpos  = BGEN64 (gpos - last_gpos), // store delta
                              .iupac = make_iupacs[i].iupac };
        last_gpos = gpos;
    }

    evb->scratch.len = sizeof (Iupac) * make_iupacs_len;
    zfile_compress_section_data_ex (evb, NULL, SEC_REF_IUPACS, &evb->scratch, 0,0, CODEC_BSC, SECTION_FLAGS_NONE, NULL); 

    buf_free (evb->scratch);
}


//--------------------------------------------------------------------------------------------
// using a reference side (iupacs are only used for lifting with --chain, not for compressing)
//--------------------------------------------------------------------------------------------

void ref_iupacs_load (Reference ref)
{
    Section sec = sections_last_sec (SEC_REF_IUPACS, SOFT_FAIL);
    if (!sec) {
        if (flag.show_ref_iupacs) 
            iprintf ("\nThere are no IUPACs in %s\n", z_name);

        goto done; // we don't have iupacs
    }

    zfile_get_global_section (SectionHeader, sec, &ref->iupacs_buf, "iupacs_buf");

    if (flag.show_ref_iupacs) 
        iprintf ("\nIUPACs found in %s:\n", z_name);

    ref->iupacs_buf.len /= sizeof (Iupac);
    ARRAY (Iupac, iupacs, ref->iupacs_buf);
    
    for (uint64_t i=0; i < iupacs_len; i++) { // first is 0
        iupacs[i].gpos = (i ? iupacs[i-1].gpos : 0) + BGEN64 (iupacs[i].gpos);

        if (flag.show_ref_iupacs)  {
            PosType32 pos;
            WordIndex chrom_index = ref_contig_get_by_gpos (ref, iupacs[i].gpos, 0, &pos);
            iprintf ("IUPAC=%c\tCHROM=%s\tPOS=%u\tGPOS=%"PRIu64"\n", 
                     iupacs[i].iupac, ctx_get_snip_by_word_index0 (ZCTX(FASTA_CONTIG), chrom_index), pos, iupacs[i].gpos);
        }
    }

done:
    if (is_genocat && flag.show_ref_iupacs) exit_ok();
}

static const Iupac *ref_iupacs_find (Iupac *iupacs, int64_t first, int64_t last, PosType64 gpos)
{
    if (first > last) return &iupacs[first]; // gpos not found in iupacs - return one after (possibly beyond end of array)

    int64_t mid = (first + last) / 2;
    int64_t delta = gpos - iupacs[mid].gpos;

    if (!delta) return &iupacs[mid];
    else return ref_iupacs_find (iupacs, delta < 0 ? first : mid+1, delta < 0 ? mid-1 : last, gpos);
}

// note: we don't really need the lower case as both reference IUPACs and the bases tested are upper-cased, but no harm specifying for future-proofing.
static const char IUPAC_IS_INCLUDED[128][128] = { ['R']={ ['A']=1, ['G']=1, ['a']=1, ['g']=1 },
                                                  ['r']={ ['A']=1, ['G']=1, ['a']=1, ['g']=1 },
                                                  ['Y']={ ['C']=1, ['T']=1, ['c']=1, ['t']=1 },
                                                  ['y']={ ['C']=1, ['T']=1, ['c']=1, ['t']=1 },
                                                  ['S']={ ['C']=1, ['G']=1, ['c']=1, ['g']=1 },
                                                  ['s']={ ['C']=1, ['G']=1, ['c']=1, ['g']=1 },
                                                  ['W']={ ['A']=1, ['T']=1, ['a']=1, ['t']=1 },
                                                  ['w']={ ['A']=1, ['T']=1, ['a']=1, ['t']=1 },
                                                  ['K']={ ['T']=1, ['G']=1, ['t']=1, ['g']=1 },
                                                  ['k']={ ['T']=1, ['G']=1, ['t']=1, ['g']=1 },
                                                  ['M']={ ['A']=1, ['C']=1, ['a']=1, ['c']=1 },
                                                  ['m']={ ['A']=1, ['C']=1, ['a']=1, ['c']=1 },
                                                  ['B']={ ['C']=1, ['G']=1, ['T']=1, ['c']=1, ['g']=1, ['t']=1 },
                                                  ['b']={ ['C']=1, ['G']=1, ['T']=1, ['c']=1, ['g']=1, ['t']=1 },
                                                  ['D']={ ['A']=1, ['G']=1, ['T']=1, ['a']=1, ['g']=1, ['t']=1 },
                                                  ['d']={ ['A']=1, ['G']=1, ['T']=1, ['a']=1, ['g']=1, ['t']=1 },
                                                  ['H']={ ['A']=1, ['C']=1, ['T']=1, ['a']=1, ['c']=1, ['t']=1 },
                                                  ['h']={ ['A']=1, ['C']=1, ['T']=1, ['a']=1, ['c']=1, ['t']=1 },
                                                  ['V']={ ['A']=1, ['C']=1, ['G']=1, ['a']=1, ['c']=1, ['g']=1 },
                                                  ['v']={ ['A']=1, ['C']=1, ['G']=1, ['a']=1, ['c']=1, ['g']=1 } };

bool ref_iupacs_is_included_do (Reference ref, VBlockP vb, const Range *r, PosType64 pos, char vcf_base)
{
    PosType64 gpos = r->gpos + (pos - r->first_pos);
    unsigned ref_i = (ref==gref); // [0]=prim_ref [1]=gref

    ARRAY (Iupac, iupacs, ref->iupacs_buf);
    const Iupac *iupac_ent = ref_iupacs_find (iupacs, 0, (int64_t)iupacs_len-1, gpos);
    bool is_after = !ref->iupacs_buf.data || BISAFT (ref->iupacs_buf, iupac_ent);

    vb->iupacs_last_range[ref_i] = r;

    // case: iupac is found at GPOS - now test if it includes vcf_base
    if (!is_after && iupac_ent->gpos == gpos) {
        vb->iupacs_last_pos[ref_i] = pos;
        vb->iupacs_next_pos[ref_i] = BISLST (ref->iupacs_buf, iupac_ent) ? 1000000000000000LL : ((iupac_ent+1)->gpos - gpos) + pos;
        bool is_included = IUPAC_IS_INCLUDED[(int)iupac_ent->iupac][(int)vcf_base];

        if (flag.show_ref_iupacs) 
            iprintf (ref_i ? "Luft: oCHROM=%s oPOS=%"PRIu64" REF=%c reference=%c is_considered_same_REF=%s\n" 
                           : "Primary: CHROM=%s POS=%"PRIu64" REF=%c reference=%c is_considered_same_REF=%s\n",
                     r->chrom_name, pos, vcf_base, iupac_ent->iupac, is_included ? "YES" : "NO");

        return is_included;
    }
    
    // case: the reference base at GPOS is not a IUPAC
    else {
        bool is_before = !ref->iupacs_buf.data || BISBEFORE(ref->iupacs_buf, iupac_ent-1);
        vb->iupacs_last_pos[ref_i] = is_before ? 0 : ((iupac_ent-1)->gpos - gpos) + pos;
        vb->iupacs_next_pos[ref_i] = is_after ? 1000000000000000LL : (iupac_ent->gpos - gpos) + pos;
        return false;
    }   
}

// iterator to check if there is a IUPAC at a position. Can be called in a sequence of calls - next call can be delayed until next_pos,
// this way this function is called scarcely. returns iupac if found (or 0). 
char ref_iupacs_get (Reference ref, const Range *r, PosType64 pos, bool reverse,
                     PosType64 *next_pos) // out 
{
    #define IUPAC2POS(iup) (r->first_pos + (iup)->gpos - r->gpos)

    // case: no iupacs in this reference
    if (!ref->iupacs_buf.len) {
        *next_pos = MAX_POS;
        return 0;
    }

    const Iupac *last = BLST(const Iupac, ref->iupacs_buf);
    PosType64 gpos = r->gpos + (pos - r->first_pos);
    
    const Iupac *iupac = ref_iupacs_find (B1ST(Iupac, ref->iupacs_buf), 0, ref->iupacs_buf.len-1, gpos);

    if (!reverse) {
        *next_pos = (iupac > last)       ? MAX_POS            // our position is beyond the last iupac
                  : (gpos < iupac->gpos) ? IUPAC2POS(iupac)   // our position is smaller than iupac - next time search when we reach this iupac
                  : (iupac < last)       ? IUPAC2POS(iupac+1) // iupac found - next, check the next iupac
                  :                        MAX_POS;           // iupac found - and this is the last iupac
            
        return (iupac <= last && iupac->gpos == gpos) ? iupac->iupac : 0;
    }
    else { // reverse
        const Iupac *first = B1ST(const Iupac, ref->iupacs_buf);

        *next_pos = (iupac > last)                         ? IUPAC2POS (last)   // our position is beyond the last iupac, next try when we get back to the last iupac
                  : (gpos < iupac->gpos && iupac == first) ? MAX_POS            // our position is smaller than the first iupac (and it will only get smaller as we are in reverse)
                  : (gpos < iupac->gpos)                   ? IUPAC2POS(iupac-1) // our position is smaller than iupac - next time search when we reach back the previous iupac
                  : (iupac > first)                        ? IUPAC2POS(iupac-1) // iupac found - next, check the previous iupac
                  :                                          MAX_POS;           // iupac found - and this is the first iupac
            
        return (iupac <= last && iupac->gpos == gpos) ? iupac->iupac : 0;    
    }
    #undef IUPAC2POS
}
