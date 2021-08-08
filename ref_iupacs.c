// ------------------------------------------------------------------
//   ref_iupacs.c
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

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
    uint32_t vblock_i;
    char iupac;
} MakeIupac;

// --make-reference: called by compute thread, when iupac is found
void ref_iupacs_add_do (VBlockP vb, uint64_t idx, char iupac)
{
    buf_alloc (vb, &CTX(FASTA_NONREF_X)->local, 1, 100, MakeIupac, 2, "contexts->local");
    NEXTENT (MakeIupac, CTX(FASTA_NONREF_X)->local) = (MakeIupac){ .vblock_i = vb->vblock_i, .idx = idx, .iupac = iupac };
}

// --make-reference: called by main thread, after each FASTA vb, in order of VBs which is also the order of ranges
void ref_iupacs_after_compute (VBlockP vb)
{
    if (CTX(FASTA_NONREF_X)->local.len) {
        buf_add_buf (evb, &gref->iupacs_buf, &CTX(FASTA_NONREF_X)->local, MakeIupac, "iupacs_buf");
        buf_free (&CTX(FASTA_NONREF_X)->local);
    }
}

// --make-reference: create a REF_IUPACS section if there are any IUPACs in this reference. The section consists of pairs (GPOS, IUPAC) sorted by GPOS
void ref_iupacs_compress (void)
{
    ARRAY (MakeIupac, make_iupacs, gref->iupacs_buf);
    if (!make_iupacs_len) return;

    ASSERTNOTINUSE (evb->compressed);
    buf_alloc (evb, &evb->compressed, 0, make_iupacs_len, Iupac, 0, "compressed");
    ARRAY (Iupac, iupacs, evb->compressed);

    if (flag.show_ref_iupacs) iprintf ("\nIUPACs found in %s:\n", z_name);

    PosType last_gpos = 0;
    for (uint64_t i=0; i < make_iupacs_len; i++) {
        const Range *r = ENT (Range, gref->ranges, make_iupacs[i].vblock_i-1);

        PosType gpos = r->gpos + make_iupacs[i].idx;

        if (flag.show_ref_iupacs)   
            iprintf ("CHROM=%s POS=%"PRIu64" GPOS=%"PRIu64" IUPAC=%c\n", 
                     ctx_get_zf_nodes_snip (ZCTX(FASTA_CONTIG), r->chrom), r->first_pos + make_iupacs[i].idx, gpos, make_iupacs[i].iupac);

        iupacs[i] = (Iupac) { .gpos  = BGEN64 (gpos - last_gpos), // store delta
                              .iupac = make_iupacs[i].iupac };
        last_gpos = gpos;
    }

    evb->compressed.len = sizeof (Iupac) * make_iupacs_len;
    zfile_compress_section_data_ex (evb, SEC_REF_IUPACS, &evb->compressed, 0,0, CODEC_BSC, SECTION_FLAGS_NONE); 

    buf_free (&evb->compressed);
}


//--------------------------------------------------------------------------------------------
// using a reference side (iupacs are only used for lifting with --chain, not for compressing)
//--------------------------------------------------------------------------------------------

void ref_iupacs_load (Reference ref)
{
    Section sec = sections_last_sec (SEC_REF_IUPACS, true);
    if (!sec) return; // we don't have iupacs

    zfile_get_global_section (SectionHeader, SEC_REF_IUPACS, sec, &ref->iupacs_buf, "iupacs_buf");

    if (flag.show_ref_iupacs) 
        iprintf ("\nIUPACs found in %s:\n", z_name);

    ref->iupacs_buf.len /= sizeof (Iupac);
    ARRAY (Iupac, iupacs, ref->iupacs_buf);
    
    for (uint64_t i=0; i < iupacs_len; i++) { // first is 0
        iupacs[i].gpos = (i ? iupacs[i-1].gpos : 0) + BGEN64 (iupacs[i].gpos);

        if (flag.show_ref_iupacs) 
            iprintf ("GPOS=%"PRIu64" IUPAC=%c\n", iupacs[i].gpos, iupacs[i].iupac);
    }

    if (exe_type == EXE_GENOCAT && flag.show_ref_iupacs) exit_ok();
}

static const Iupac *ref_iupacs_find (Iupac *iupacs, int64_t first, int64_t last, PosType gpos)
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

bool ref_iupacs_is_included_do (VBlockP vb, const Range *luft_range, PosType opos, char vcf_base)
{
    PosType gpos = luft_range->gpos + (opos - luft_range->first_pos);

    ARRAY (Iupac, iupacs, gref->iupacs_buf);
    const Iupac *iupac_ent = ref_iupacs_find (iupacs, 0, (int64_t)iupacs_len-1, gpos);
    bool is_after = ISAFTERENT (gref->iupacs_buf, iupac_ent);

    vb->iupacs_last_range = luft_range;

    // case: iupac is found at GPOS - now test if it includes vcf_base
    if (!is_after && iupac_ent->gpos == gpos) {
        vb->iupacs_last_opos = opos;
        vb->iupacs_next_opos = ISLASTENT (gref->iupacs_buf, iupac_ent) ? 1000000000000000LL : ((iupac_ent+1)->gpos - gpos) + opos;
        bool is_included = IUPAC_IS_INCLUDED[(int)iupac_ent->iupac][(int)vcf_base];

        if (flag.show_ref_iupacs) 
            iprintf ("oCHROM=%s oPOS=%"PRIu64" REF=%c reference=%c is_considered_same_REF=%s\n",
                     luft_range->chrom_name, opos, vcf_base, iupac_ent->iupac, is_included ? "YES" : "NO");

        return is_included;
    }
    
    // case: the reference base at GPOS is not a IUPAC
    else {
        bool is_before = ISBEFOREENT(gref->iupacs_buf, iupac_ent-1);
        vb->iupacs_last_opos = is_before ? 0 : ((iupac_ent-1)->gpos - gpos) + opos;
        vb->iupacs_next_opos = is_after ? 1000000000000000LL : (iupac_ent->gpos - gpos) + opos;
        return false;
    }   
}
