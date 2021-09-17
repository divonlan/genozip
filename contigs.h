// ------------------------------------------------------------------
//   contigs.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef CONTIGS_INCLUDED
#define CONTIGS_INCLUDED

#include "genozip.h"
#include "buffer.h"

//-------------------------------------------------------------------------
// these structs are part of the Genozip file format
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
    WordIndex chrom_index;
    PosType min_pos, max_pos;    // POS field value of smallest and largest POS value of this contig
    PosType gpos;                // The GPOS in genome matching min_pos in contig.
    ContigMetadata metadata;
} Contig; 

#pragma pack()

typedef enum { SORT_BY_NONE=0, SORT_BY_NAME=1, SORT_BY_AC=2, SORT_BY_ALL=3 } SortBy;

typedef struct ContigPkg {
    const char *name;
    bool has_counts;
    Buffer contigs;
    Buffer dict;
    Buffer by_name, by_LN, by_AC; // sorters
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
extern WordIndex contigs_get_matching (ConstContigPkgP ctgs, STRp(name), PosType LN /* optional */, bool *is_alt);

// accession numbers
typedef struct { char s[ACCESSION_LEN+5]; } AccNumText;
extern AccNumText display_acc_num (const AccessionNumber *ac);

#define CONTIG(ctg_pkg,ctg_i) ENT (Contig, ((ctg_pkg).contigs), (ctg_i))
#endif
