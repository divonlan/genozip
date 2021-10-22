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
    float b250_per_line[MAX_DICTS]; // b250.len / num_lines
    #define AT_LEAST(did_i) ((uint64_t)(10.0 + (segconf.b250_per_line[did_i] * (float)(vb->lines.len))))

    // SAM/BAM stuff
    bool sam_use_aligner;       // use of aligner is possible if its flag.aligner_available and there are no header contigs
    bool NM_is_integer;         // true if NM is integer, false if it binary
    bool has_MC, has_MD, has_MQ;// MC, MD, MQ field was detected in the data
    enum { XA_NONE, XA_BWA, XA_IONTORRENT, XA_UNKNOWN } has_XA; // IonTorret and BWA have different XA:Z
    bool sam_is_collated;       // Every QNAME appears in two or more consecutive lines
    bool sam_is_sorted;         // every two consecutive lines that have the same RNAME, have non-decreasing POS
    bool sam_buddy_RG;          // attempt to use the same mate for RG:Z as QNAME
    uint64_t sam_cigar_len;     // approx average CIGAR len (during running==true - total len)
    int64_t MAPQ_value;         // used during segconf.running to calculate sam_mapq_has_single_value
    bool MAPQ_has_single_value; // all non-0 MAPQ have the same value

    // FASTA stuff
    bool fasta_has_contigs;     // the sequences in this FASTA represent contigs (as opposed to reads) - in which case we have a FASTA_CONTIG dictionary and RANDOM_ACCESS

    // Chain stuff
    bool chain_mismatches_ref;  // Some contigs mismatch the reference files, so this chain file cannot be used with --chain

    // read name characteristics (SAM/BAM, KRAKEN and FASTQ)
    unsigned qname_flavor;
} SegConf;

extern SegConf segconf;

extern void segconf_initialize (void);
extern void segconf_calculate (void);

#endif
