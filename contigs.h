// ------------------------------------------------------------------
//   contigs.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef CONTIGS_INCLUDED
#define CONTIGS_INCLUDED

#include "genozip.h"
#include "buffer.h"

//--------------------------------------------------------------------------------------------------------------
// NOTE: these structs are part of the Genozip file format - section SEC_REF_CONTIGS contains an array of Contig
//--------------------------------------------------------------------------------------------------------------

#pragma pack(1) // structures that are part of the genozip format are packed.

// Accession Number: https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
typedef struct {
    char version; // version character is version , eg '1'. If version is not specified, then it is set to '1'. 
    #define ACCESSION_LEN 15 // max by according to ^^^ is 12
    char AC[ACCESSION_LEN];  // upper-case letters followed by numerals ; zero padded
} AccessionNumber; 

typedef union {
    #define REFCONTIG_MD_LEN 96
    char str[REFCONTIG_MD_LEN];  // as appears in a reference file on disk: nul-termianted/padded string: Properties, as they appear in the DESC line of the reference FASTA. note: prior to v12 this space was occupied by specific fields, which were never utilized
    struct {
        int64_t count;
        AccessionNumber ac;      // as appears in a loaded reference file in memory
    } parsed;
} ContigMetadata;

typedef struct Contig {
    CharIndex char_index;        // char index in CHROM dictionary of this contig
    uint32_t snip_len;
    WordIndex ref_index;         // index in reference
    PosType min_pos, max_pos;    // POS field value of smallest and largest POS value of this contig
    PosType gpos;                // The GPOS in genome matching min_pos in contig.
    ContigMetadata metadata;
} Contig; 

#pragma pack()

typedef enum { SORT_BY_NONE=0, SORT_BY_NAME=1, SORT_BY_AC=2, SORT_BY_REF_INDEX=4 } SortBy;

typedef struct ContigPkg {
    #define cp_next        contigs.param // iterator
    #define cp_num_contigs contigs.len   
    const char *name;
    bool has_counts;
    Buffer contigs;
    Buffer dict;
    Buffer by_name, by_LN, by_AC, by_ref_index; // sorters
    SortBy sorted_by;
} ContigPkg;

//-------------------------------------------------------------------------

// initialization & finalization
extern void contigs_build_contig_pkg_from_ctx (ContigPkg *ctgs, ConstContextP ctx, SortBy sort_by);
extern void contigs_create_index (ContigPkg *ctgs, SortBy sort_by);
extern void contigs_free (ContigPkg *ctg);
extern void contigs_destroy (ContigPkg *ctg);

// finding
#define WORD_INDEX_NOT_UNIQUE (-2)
extern WordIndex contigs_get_by_name (ConstContigPkgP ctgs, STRp(contig_name));
static inline ContigP contigs_get_by_index (ConstContigPkgP ctgs, WordIndex index) { return ENT (Contig, ctgs->contigs, index); }
extern WordIndex contigs_get_matching (ConstContigPkgP ctgs, STRp(name), PosType LN /* optional */, bool strictly_alt, bool *is_alt);
extern const char *contigs_get_name (ConstContigPkgP ctgs, WordIndex index, unsigned *contig_name_len /* optional */);
extern WordIndex contigs_get_by_ref_index (ConstContigPkgP ctgs,WordIndex ref_index);
static inline PosType contigs_get_LN (ConstContigPkgP ctgs, WordIndex index) { return ENT (Contig, ctgs->contigs, index)->max_pos; }

// iterator
typedef void (*ContigsIteratorCallback)(STRp(contig_name), PosType last_pos, void *callback_param);
extern void foreach_contig (ConstContigPkgP ctgs, ContigsIteratorCallback callback, void *callback_param);

// accession numbers
typedef struct { char s[ACCESSION_LEN+5]; } AccNumText;
extern AccNumText display_acc_num (const AccessionNumber *ac);

#define CONTIG(ctg_pkg,ctg_i) ENT (Contig, ((ctg_pkg).contigs), (ctg_i))
#endif