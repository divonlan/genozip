// ------------------------------------------------------------------
//   sam_sa.c
//   Copyright (C) 2021-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

// ---------------------------------------------------------
// SA:Z "Other canonical alignments in a chimeric alignment"
//
// SA format is: (rname, pos, strand, CIGAR, mapQ, NM ;)+ 
// Example SA:Z:chr13,52863337,-,56S25M70S,0,0;chr6,145915118,+,97S24M30S,0,0;chr18,64524943,-,13S22M116S,0,0;chr7,56198174,-,20M131S,0,0;chr7,87594501,+,34S20M97S,0,0;chr4,12193416,+,58S19M74S,0,0;
// See: https://samtools.github.io/hts-specs/SAMtags.pdf
// ---------------------------------------------------------

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "chrom.h"


// Alignments containing SA are removed from the SAM_COMP_MAIN vbs and moved to PRIM / DEPN components.
// When segging the PRIM component (containing primary SA lines) we place the data in the data structures below.
// When segging the DEPN component, we seg referring to the data structures.
//
// Seg of PRIM: vb->{sa_groups,sa_alns} is segging a PRIM vb and are merged into 
// z_file->{sa_groups,sa_alns,sa_qnames,sa_seq,sa_qual} during merge. 

// note: fields ordered to packed and word-aligned
typedef struct __attribute__ ((__packed__)) {
    uint64_t qname;       // index into: vb: txt_data ; z_file: zfile->sa_qnames
    uint64_t seq;         // index into: vb: txt_data ; z_file: zfile->sa_seq
    uint64_t qual;        // index into: vb: txt_data ; z_file: zfile->sa_qual
    uint16_t qname_len;
    uint16_t num_alns;    // number of alignments in this SA group (including primary alignment)
    uint32_t first_aln;   // index into sa_alns
    uint32_t seq;         // index into: vb: txt_data ; z_file: zfile->sa_seq. seq_len is (qual_len+1)/2
    uint32_t qual;        // index into: vb: txt_data ; z_file: zfile->sa_qual
    uint32_t qual_len;
} SAGroupType;

typedef struct __attribute__ ((__packed__)) {
    uint64_t cigar;       // index into: vb : txt_data ; z_file: zfile->sa_cigars. Note: SA cigar is textual form even in BAM
                          // BAM vb primary alignment: index into vb->sa_cigars (because BAM primary alignment CIGAR is binary)
    uint32_t cigar_len;
    WordIndex rname;      // vb: node_index z_file: word_index
    uint32_t pos;         // 31 bits per SAM spec
    uint32_t mapq   : 8;  // 8 bit by SAM specification
    uint32_t strand : 1;  // 0 for - and 1 for +
    uint32_t nm     : 23; // plenty, NM is very rarely above 1000
} SAAlnType;

static void sam_sa_add_sa_group (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(field))
{
    uint64_t save_sa_groups_len      = vb->sa_groups.len;
    uint64_t save_sa_alns_len        = vb->sa_alns.len;
    uint64_t save_sa_prim_cigars_len = vb->sa_prim_cigars.len;
    
    // check that bits values are within limits defined in SAGroupType and SAAlnType - otherwise, we cannot add this SA
    if ((dl->QNAME.len > 65535) || (dl->POS >= (1<<32)) || (dl->MAPQ >= (1<<8)) || (vb->NM >= (1<<23)))
        return;

    str_split (field, field_len, 0, ';', aln, false);
    if (!n_alns || n_alns > 65536 /*-- in a sec*/) return; // we cannot add this SA - either invalid or over 65535 alignments (SAGroupType.n_alns in 16 bit)

    uint32_t num_depn_alns = n_alns - 1; // drop extra item due to terminal ';'

    buf_alloc (vb, &vb->sa_groups, 1, 32, SAGroupType, CTX_GROWTH, "sa_groups");
    buf_alloc (vb, &vb->sa_alns, num_depn_alns+1 /*+1 for primary aln*/, 64, SAAlnType, CTX_GROWTH, "sa_alns");

    // add SA group
    NEXTENT (SAGroupType, vb->sa_groups) = (SAGroupType){
        .first_aln = vb->sa_alns.len,
        .num_alns  = num_depn_alns + 1,
        .qname     = dl->QNAME.index,
        .qname_len = dl->QNAME.len,
        .qual      = dl->QUAL.index,
        .qual_len  = dl->QUAL.len,
        .seq       = dl->seq_data_start
    };

    // add primary alignment as first alignment in SA group
    NEXTENT (SAAlnType, vb->sa_alns) = (SAAlnType){
        .rname     = vb->chrom_node_index,
        .pos       = dl->POS,
        .cigar     = IS_BAM ? vb->sa_prim_cigars.len : dl->CIGAR.index,
        .cigar_len = IS_BAM ? vb->textual_cigar.len  : dl->CIGAR.len,
        .mapq      = dl->MAPQ,
        .strand    = !dl->FLAG.bits.rev_comp,
        .nm        = vb->NM
    };

    // in BAM, add textual CIGAR to vb->sa_cigars, as txt_data CIGAR is binary
    if (IS_BAM)
        buf_add_buf (vb, &vb->sa_prim_cigars, &vb->textual_cigar, char, "sa_prim_cigars");

    // add all dependent alignments
    for (uint32_t i=0; i < num_depn_alns; i++) {
        
        // get items - rname, pos, strand, CIGAR, mapQ, NM
        str_split (alns[i], aln_lens[i], 6, ',', item, true);
        if (n_items != 6) goto rollback;

        // rname
        WordIndex rname = ctx_get_existing_node_index (vb, STRi(item,0));
        if (rname == WORD_INDEX_NONE) goto rollback;

        // pos, mapq, nm - get integers and verify limits
        int64_t pos, mapq, nm;
        if (!str_get_int_range64 (STRi(item,1), 0, 0xffffffff, &pos )) goto rollback;
        if (!str_get_int_range64 (STRi(item,4), 0, 255,        &mapq)) goto rollback;
        if (!str_get_int_range64 (STRi(item,5), 0, (1<<23)-1,  &nm  )) goto rollback;

        // strand
        bool positive = *items[2] == '+';
        if (item_lens[2] != 1 || (!positive && *items[2] != '-')) goto rollback;

        NEXTENT (SAAlnType, vb->sa_alns) = (SAAlnType)
            { .rname = rname, .pos = pos, .strand = positive, .cigar = items[3], .cigar_len = item_lens[3], .nm = nm };
    }

    return;

rollback:
    vb->sa_groups.len      = save_sa_groups_len;
    vb->sa_alns.len        = save_sa_alns_len;
    vb->sa_prim_cigars.len = save_sa_prim_cigars_len;
}       

void sam_seg_SA_field (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(field))
{
    bool is_prim = (flag.gencomp_num == SAM_COMP_PRIM);
    ASSERT (is_prim || flag.gencomp_num == SAM_COMP_DEPN, "Unexpected gencomp_num=%u", flag.gencomp_num);
    
    static const MediumContainer container_SA = { .nitems_lo = 6,      
                                                  .repsep    = { ';' }, // including on last repeat    
                                                  .items     = { { .dict_id = { _OPTION_SA_RNAME  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_POS    }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_STRAND }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_CIGAR  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_MAPQ   }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_NM     },                  } } };

    if (is_prim) {

        // seg normally
        SegCallback callbacks[6] = { [0]=chrom_seg_cb, [1]=seg_pos_field_cb, [3]=sam_seg_0A_cigar_cb };            
        seg_array_of_struct (VB, CTX(OPTION_SA_Z), container_SA, STRa(field), callbacks);

        // build SA data structure in memory
        if (!segconf.running) 
            sam_sa_add_sa_group (vb, dl, STRa (field));
    }
    
    // supplamentary / secondary alignment: seg against SA data structure
    else {

    }

    CTX(OPTION_SA_Z)->txt_len++; // 1 for \t in SAM and \0 in BAM 
}

void sam_sa_set_NM (VBlockSAMP vb, STRps(aux), bool is_bam)
{
    STR(NM);
    NM = sam_seg_get_aux ("NM:i", STRas(aux), &NM_len, is_bam);
    
    if (!NM || !str_get_int (STRa(NM), &vb->NM)) 
        vb->NM = 0xfffffffffffffULL; // (value > 2^23) = no value NM value on this line
}
