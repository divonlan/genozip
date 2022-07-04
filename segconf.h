// ------------------------------------------------------------------
//   segconf.c
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

// Documented range for users
#define MIN_VBLOCK_MEMORY  1    // in MB 
#define MAX_VBLOCK_MEMORY  2048 

// Range in developer options eg --vblock 10000B 
#define ABSOLUTE_MIN_VBLOCK_MEMORY ((uint64_t)1000) // in Bytes
#define ABSOLUTE_MAX_VBLOCK_MEMORY ((uint64_t)MAX_VBLOCK_MEMORY<<20)

typedef enum __attribute__ ((__packed__)) { TECH_UNKNOWN, TECH_ILLUM_7, TECH_ILLUM_5, TECH_PACBIO, TECH_ONP, TECH_454, TECH_BGI, TECH_IONTORR, TECH_HELICOS } SeqTech;

typedef enum __attribute__ ((__packed__)) { SQT_UNKNOWN, SQT_NUKE, SQT_AMINO, SQT_NUKE_OR_AMINO } SeqType;

typedef enum __attribute__ ((__packed__)) { PL_mux_by_DP_TEST, PL_mux_by_DP_NO, PL_mux_by_DP_YES } PLMuxByDP;

typedef enum __attribute__ ((__packed__)) { ms_NONE, ms_BIOBAMBAM, ms_MINIMAP2 } msType; // type of SAM ms:i field 

typedef enum __attribute__ ((__packed__)) { DP_DEFAULT, by_AD, by_SDP } FormatDPMethod;

typedef enum __attribute__ ((__packed__)) { L3_UNKNOWN, L3_EMPTY, L3_COPY_DESC, L3_QF, NUM_L3s } FastqLine3Type;

typedef enum __attribute__ ((__packed__)) { XG_S_UNKNOWN, XG_WITHOUT_S, XG_WITH_S } XgIncSType;

typedef enum __attribute__ ((__packed__)) { MP_UNKNOWN, MP_BSBOLT,             MP_bwa,   MP_BWA,   MP_MINIMAP2,   MP_STAR,   MP_BOWTIE2,   MP_DRAGEN,    MP_GEM3,         MP_BISMARK,   MP_BSSEEKER2,     MP_WINNOWMAP,   MP_BAZ2BAM,    MP_BBMAP,   MP_TMAP,   MP_HISAT2,   MP_BOWTIE,   MP_NOVOALIGN, MP_RAZER3,    MP_BLASR,   MP_NGMLR,           MP_DELVE,   MP_TOPHAT,   NUM_MAPPERS    } SamMapperType;
#define SAM_MAPPER_NAME             { "Unknown_mapper", "bsbolt",              "bwa",    "BWA",    "minimap2",    "STAR",    "bowtie2",    "dragen",     "gem3",          "bismark",    "bsseeker2",      "Winnowmap",    "baz2bam",     "BBMap",    "tmap",    "hisat2",    "Bowtie",    "NovoAlign",  "razers3",    "blasr",    "ngmlr",            "Delve",    "TopHat",    }
#define SAM_MAPPER_SIGNATURE        { "Unknown_mapper", "PN:bwa	VN:BSB"/*\t*/, "PN:bwa", "PN:BWA", "PN:minimap2", "PN:STAR", "PN:bowtie2", "ID: DRAGEN", "PN:gem-mapper", "ID:Bismark", "PN:BS Seeker 2", "PN:Winnowmap", "PN:baz2bam",  "PN:BBMap", "ID:tmap", "PN:hisat2", "ID:Bowtie", "xxx",        "PN:razers3", "ID:BLASR", "PN:nextgenmap-lr", "ID:Delve", "ID:TopHat", }

// seg configuration set prior to starting to seg a file during segconfig_calculate or txtheader_zip_read_and_compress
typedef struct {

    // Seg parameters - general
    uint64_t vb_size;           // ZIP/PIZ: compression VBlock size in bytes (PIZ: passed in SectionHeaderGenozipHeader.vb_size)
    bool running;               // currently in segconf_calculate()
    bool has[MAX_DICTS];        // for select did_i's, states whether this field was encountered during segconf.running
    uint32_t line_len;          // approx line len
    float b250_per_line[MAX_DICTS]; // b250.len / num_lines
    #define AT_LEAST(did_i) ((uint64_t)(10.0 + (segconf.b250_per_line[did_i] * (float)(vb->lines.len32))))

    // read characteristics (SAM/BAM, KRAKEN and FASTQ)
    QnameFlavor qname_flavor, qname_flavor2;  
    SeqTech tech;

    // SAM/BAM and FASTQ
    uint32_t longest_seq_len;   // length of the longest seq_len in the segconf data 
    DictId qname_seq_len_dict_id; // dict_id of one of the Q?NAME contexts, which is expected to hold the seq_len for this read. 0 if there is no such item.
    
    // SAM/BAM stuff
    STRl (std_cigar, 16);       // first CIGAR in the file - used in case all CIGARs in the file are the same
    bool sam_is_unmapped;       // all POS fields in the segconf block were 0
    SamMapperType sam_mapper;           
    bool NM_is_integer;         // true if NM is integer, false if it binary
    bool has_TLEN_non_zero;
    bool has_DP_before_PL;
    bool is_collated;       // Every QNAME appears in two or more consecutive lines
    bool evidence_of_collated;  // during segconf: at least a one pair of consecutive lines has the same QNAME
    bool is_sorted;         // every two consecutive lines that have the same RNAME, have non-decreasing POS
    bool is_paired;         // file has a least one read that is marked as "last" in FLAG
    bool evidence_of_sorted;    // during segconf: at least a one pair of consecutive lines has the same RNAME and increasing POS
    bool sam_multi_RG;          // evidence that file has more than one type of RG
    bool sam_bisulfite;         // this BAM file is of reads that have been treated with bisulfite to detect methylation
    bool has_MD_or_NM;          // ZIP/PIZ: call sam_analyze_copied_SEQ for SEQ copied from prim unless no cigar, no seq or (PIZ) explicitly told not to.
    bool NM_after_MD;           // in all segconf lines that had both NM and MD, NM appeared after MD
    bool pysam_qual;            // ZIP/PIZ: BAM missing QUAL is generated by old versions of pysam: first byte is 0xff, followed by 0s (SAM spec requires all bytes to be 0xff)
    bool sam_has_depn;          // ZIP: true if any line in the segconf sample had SAM_FLAG_SECONDARY or SAM_FLAG_SUPPLEMENTARY
    bool star_solo;             // ZIP: using STARsolo or cellranger
    uint32_t AS_is_2ref_consumed; // ZIP/PIZ: AS value tends to be double ref_consumed (counter during segconf, bool during seg/piz)
    uint8_t MAPQ_value;         // used during segconf.running to calculate sam_mapq_has_single_value
    bool MAPQ_has_single_value; // all non-0 MAPQ have the same value
    msType sam_ms_type;         // ZIP/PIZ: type of ms:i 
    XgIncSType sam_XG_inc_S;    // Does XG include soft_clip[0]
    bool is_long_reads;
    SagType sag_type;           // Type of SA Groups
    bool sag_has_AS;       // SA Groups store the AS:i values of prim lines that have them. Set if its beneficial to seg AS:i in depn lines against prim
    uint32_t sam_cigar_len;     // approx average CIGAR len rounded up (during running==true - total len)
    uint32_t sam_seq_len;       // ZIP/PIZ: approx average (SEQ.len+hard-clips) rounded to the nearest (during running - total len)
    uint32_t seq_len_to_cm;     // ZIP/PIZ: approx average of (seq_len/cm:i) (during running - cumulative)
    uint8_t n_CR_CB_seps;
    char CR_CB_seperator;       // ZIP: seperator within CR:Z and CB:Z fields
    MiniContainer CY_con;
    STRl(CY_con_snip, 64);

    // SAM/BAM and FASTQ
    bool nontrivial_qual;       // true if we know that not all QUAL values are the same (as they are in newer PacBio files)

    // VCF stuff
    bool vcf_is_varscan;        // this VCF file was produced by VarScan
    uint64_t count_dosage[2];   // used to calculate pc_has_dosage
    float pc_has_dosage;        // % of the samples x lines that have a valid (0-2) dosage value [0.0,1.0]
    bool use_null_DP_method;    // A method for predicting GT=./. by DP=.
    FormatDPMethod FORMAT_DP_method;
    PLMuxByDP PL_mux_by_DP;
    Mutex PL_mux_by_DP_mutex;
    uint64_t count_GQ_by_PL, count_GQ_by_GP; // used tp calculate GQ_by_PL, GQ_by_GP
    bool GQ_by_PL, GQ_by_GP;
    
    // FASTQ
    FastqLine3Type line3;       // format of line3
    QnameFlavor line3_flavor;   // in case of L3_QF 

    // FASTA stuff
    bool fasta_has_contigs;     // the sequences in this FASTA represent contigs (as opposed to reads) - in which case we have a FASTA_CONTIG dictionary and RANDOM_ACCESS
    SeqType seq_type;           // nucleotide or protein
    unsigned seq_type_counter;  // used for calculating seq_type 

    // Chain stuff
    bool chain_mismatches_ref;  // Some contigs mismatch the reference files, so this chain file cannot be used with --chain
} SegConf;

extern SegConf segconf; // ZIP: set based on segging a sample of a few first lines of the file
                        // PIZ: select fields are transferred through SectionHeaderGenozipHeader

extern void segconf_initialize (void);
extern void segconf_calculate (void);
extern void segconf_update_qual (STRp (qual));
extern bool segconf_is_long_reads(void);
extern void segconf_mark_as_used (VBlockP vb, unsigned num_ctxs, ...);

static inline void segconf_set_has (Did did_i)
{
    if (segconf.running) segconf.has[did_i] = true;
}
