// ------------------------------------------------------------------
//   segconf.c
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

#define MIN_VBLOCK_MEMORY  1    // in MB
#define MAX_VBLOCK_MEMORY  2048 

typedef enum { TECH_UNKNOWN, TECH_ILLUM_7, TECH_ILLUM_5, TECH_PACBIO, TECH_ONP, TECH_454, TECH_BGI, TECH_IONTORR, TECH_HELICOS, TECH_NCBI } SeqTech;

typedef enum { SQT_UNKNOWN, SQT_NUKE, SQT_AMINO, SQT_NUKE_OR_AMINO } SeqType;

typedef enum { PL_mux_by_DP_TEST, PL_mux_by_DP_NO, PL_mux_by_DP_YES } PLMuxByDP;

typedef enum { PS_NONE, PS_POS, PS_POS_REF_ALT, PS_UNKNOWN } PSType;

typedef enum { DP_DEFAULT, by_AD, by_SDP, by_INFO_DP } FormatDPMethod;

// seg configuration set prior to starting to seg a file during segconfig_calculate or txtheader_zip_read_and_compress
typedef struct {

    // Seg parameters - general
    uint64_t vb_size;
    bool running;               // currently in segconf_calculate()
    uint32_t line_len;          // approx line len
    float b250_per_line[MAX_DICTS]; // b250.len / num_lines
    #define AT_LEAST(did_i) ((uint64_t)(10.0 + (segconf.b250_per_line[did_i] * (float)(vb->lines.len))))

    bool has[MAX_DICTS];        // for select did_i's, states whether this field was encountered during segconf.running

    // SAM/BAM stuff
    bool sam_use_aligner;       // use of aligner is possible if its flag.aligner_available and there are no header contigs
    bool NM_is_integer;         // true if NM is integer, false if it binary
    bool has_TLEN_non_zero;
    bool has_DP_before_PL;
    enum { XA_NONE, XA_BWA, XA_IONTORRENT, XA_UNKNOWN } has_XA; // IonTorret and BWA have different XA:Z
    bool sam_is_collated;       // Every QNAME appears in two or more consecutive lines
    bool sam_is_sorted;         // every two consecutive lines that have the same RNAME, have non-decreasing POS
    bool sam_buddy_RG;          // attempt to use the same mate for RG:Z as QNAME
    uint64_t sam_cigar_len;     // approx average CIGAR len (during running==true - total len)
    uint64_t sam_seq_len;       // approx average CIGAR len (during running==true - total len)
    int64_t MAPQ_value;         // used during segconf.running to calculate sam_mapq_has_single_value
    bool MAPQ_has_single_value; // all non-0 MAPQ have the same value
    
    // VCF stuff
    bool vcf_is_varscan;        // this VCF file was produced by VarScan
    uint64_t count_dosage[2];   // used to calculate pc_has_dosage
    float pc_has_dosage;        // % of the samples x lines that have a valid (0-2) dosage value [0.0,1.0]
    PSType ps_type;
    bool use_null_DP_method;    // A method for predicting GT=./. by DP=.
    bool INFO_DP_by_FORMAT_DP;  // INFO_DP = Delta vs SUM(FORMAT_DP) 
    FormatDPMethod FORMAT_DP_method;

    PLMuxByDP PL_mux_by_DP;
    Mutex PL_mux_by_DP_mutex;

    uint64_t count_GQ_by_PL, count_GQ_by_GP; // used tp calculate GQ_by_PL, GQ_by_GP
    bool GQ_by_PL, GQ_by_GP;
    
    // FASTA stuff
    bool fasta_has_contigs;     // the sequences in this FASTA represent contigs (as opposed to reads) - in which case we have a FASTA_CONTIG dictionary and RANDOM_ACCESS
    
    SeqType seq_type;           // nucleotide or protein
    unsigned seq_type_counter;  // used for calculating seq_type 

    // Chain stuff
    bool chain_mismatches_ref;  // Some contigs mismatch the reference files, so this chain file cannot be used with --chain

    // read name characteristics (SAM/BAM, KRAKEN and FASTQ)
    unsigned qname_flavor, qname_flavor2;  
    SeqTech tech;

} SegConf;

extern SegConf segconf;

extern void segconf_initialize (void);
extern void segconf_calculate (void);
static inline bool segconf_is_long_reads(void) { return segconf.tech == TECH_PACBIO || segconf.tech == TECH_ONP; }
