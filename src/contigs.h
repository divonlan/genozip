// ------------------------------------------------------------------
//   contigs.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"
#include "buffer.h"

//--------------------------------------------------------------------------------------------------------------
// NOTE: these structs are part of the Genozip file format - section SEC_REF_CONTIGS contains an array of Contig
//--------------------------------------------------------------------------------------------------------------

#pragma pack(1) // structures that are part of the genozip format are packed.

// Accession Number: https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
typedef struct {
    char version;                // version character is version , eg '1'. If version is not specified, then it is set to '1'. 
    #define ACCESSION_LEN 15     // max bytes according to ^^^ is 12
    char AC[ACCESSION_LEN];      // upper-case letters followed by numerals ; zero padded
    char version2;               // second digit of version (since v13)
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
    PosType64 min_pos, max_pos;  // POS field value of smallest and largest POS value of this contig
    PosType64 gpos;              // The GPOS in genome matching min_pos in contig.
    ContigMetadata metadata;
} Contig; 

// used for the SEC_REF_CONTIG section in non-refernce files, starting 14.0.10
typedef struct CompactContig {
    WordIndex ref_index;         // index in reference. zero if sequential_ref_index=true
    uint32_t min_pos, max_pos;   // POS field value of smallest and largest POS value of this contig
    PosType64 gpos;
} CompactContig; 

#pragma pack()

typedef enum { SORT_BY_NONE=0, SORT_BY_NAME=1, SORT_BY_AC=2, SORT_BY_REF_INDEX=4, SORT_BY_LN=8 } SortBy;

typedef struct ContigPkg {
    #define cp_next        contigs.param // iterator
    #define cp_num_contigs contigs.len   
    uint64_t unique_id; // id of this specific instance, unique across the entire execution of possibley multiple files
    rom name;
    bool has_counts;
    Buffer contigs;
    Buffer dict;
    Buffer by_name, by_LN, by_AC, by_ref_index; // sorters
    SortBy sorted_by;
} ContigPkg;

//-------------------------------------------------------------------------

// initialization & finalization
extern void contigs_create_index (ContigPkg *ctgs, SortBy sort_by);
extern void contigs_free (ContigPkg *ctg);
extern void contigs_destroy (ContigPkg *ctg);

// finding
#define WORD_INDEX_NOT_UNIQUE (-2)
extern WordIndex contigs_get_by_name (ConstContigPkgP ctgs, STRp(contig_name));
static inline ContigP contigs_get_by_index (ConstContigPkgP ctgs, WordIndex index) { return B(Contig, ctgs->contigs, index); }
extern WordIndex contigs_get_matching (ConstContigPkgP ctgs, STRp(name), PosType64 LN /* optional */, bool strictly_alt, bool *is_alt);
extern rom contigs_get_name (ConstContigPkgP ctgs, WordIndex index, unsigned *contig_name_len /* optional */);
extern WordIndex contigs_get_by_ref_index (ConstContigPkgP ctgs,WordIndex ref_index);
static inline PosType64 contigs_get_LN (ConstContigPkgP ctgs, WordIndex index) { return index < ctgs->contigs.len32 ? B(Contig, ctgs->contigs, index)->max_pos : 0; }
static inline PosType64 contigs_get_gpos (ConstContigPkgP ctgs, WordIndex index) { return B(Contig, ctgs->contigs, index)->gpos; }
extern uint64_t contigs_get_nbases (ConstContigPkgP ctgs);

// iterator
typedef void (*ContigsIteratorCallback)(STRp(contig_name), PosType64 last_pos, void *callback_param);
extern void foreach_contig (ConstContigPkgP ctgs, ContigsIteratorCallback callback, void *callback_param);

// accession numbers
typedef struct { char s[ACCESSION_LEN+20]; } AccNumText;
extern AccNumText display_acc_num (const AccessionNumber *ac);

#define CONTIG(ctg_pkg,ctg_i) B(Contig, ((ctg_pkg).contigs), (ctg_i))

