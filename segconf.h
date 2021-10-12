// ------------------------------------------------------------------
//   segconf.c
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef SEGCONF_INCLUDED
#define SEGCONF_INCLUDED

#include "genozip.h"

#define MIN_VBLOCK_MEMORY  1    // in MB
#define MAX_VBLOCK_MEMORY  2048 

// seg configuration set prior to starting to seg a file during segconfig_calculate or txtheader_zip_read_and_compress
typedef struct {

    // Seg parameters - general
    uint64_t vb_size;
    bool running;               // currently in segconf_calculate()
    uint32_t line_len;          // approx line len

    // SAM/BAM stuff
    bool sam_use_aligner;       // use of aligner is possible if its flag.aligner_available and there are no header contigs
    bool NM_is_integer;         // true if NM is integer, false if it binary
    bool has_MC, has_MD;        // MC, MD field was detected in the data
    enum { XA_NONE, XA_BWA, XA_IONTORRENT, XA_UNKNOWN } has_XA; // IonTorret and BWA have different XA:Z
    bool sam_is_collated;       // Every QNAME appears in two or more consecutive lines
    bool sam_is_sorted;         // every two consecutive lines that have the same RNAME, have non-decreasing POS
    bool sam_buddy_RG;          // attempt to use the same mate for RG:Z as QNAME
    uint64_t sam_cigar_len;     // approx average CIGAR len (during running==true - total len)
    uint32_t XA_reps, SA_reps, OA_reps; // approx average repeats per line of XA, SA, OA elements

    // FASTA stuff
    bool fasta_has_contigs;     // the sequences in this FASTA represent contigs (as opposed to reads) - in which case we have a FASTA_CONTIG dictionary and RANDOM_ACCESS

    // Chain stuff
    bool chain_mismatches_ref;  // Some contigs mismatch the reference files, so this chain file cannot be used with --chain

    // read name characteristics (SAM/BAM and FASTQ)
    bool is_bgi_E9L1C3R10;      // Read names have a fixed-length format that looks like "E100020409L1C001R0030000801". 
    bool is_illumina_7;         // Read names look like: "A00488:61:HMLGNDSXX:4:1101:4345:1000"
} SegConf;

extern SegConf segconf;

extern void segconf_initialize (void);
extern void segconf_calculate (void);

#endif
