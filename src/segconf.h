// ------------------------------------------------------------------
//   segconf.c
//   Copyright (C) 2021-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

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
#define ms_type_NAME { "None", "biobambam", "minimap2"}

typedef enum __attribute__ ((__packed__)) { DP_DEFAULT, by_AD, by_SDP } FormatDPMethod;

typedef enum __attribute__ ((__packed__)) { L3_UNKNOWN, L3_EMPTY, L3_COPY_DESC, L3_QF, NUM_L3s } FastqLine3Type;

typedef enum __attribute__ ((__packed__)) { XG_S_UNKNOWN, XG_WITHOUT_S, XG_WITH_S } XgIncSType;
#define XG_INC_S_NAME { "Unknown", "without-S", "with-S"}

// SamMapperType is part of the file format and values should not be changed (new ones can be added)
typedef enum __attribute__ ((__packed__)) { MP_UNKNOWN, MP_BSBOLT,             MP_bwa,   MP_BWA,   MP_MINIMAP2,   MP_STAR,   MP_BOWTIE2,   MP_DRAGEN,    MP_GEM3,         MP_GEM2SAM,     MP_BISMARK,   MP_BSSEEKER2,     MP_WINNOWMAP,   MP_BAZ2BAM,    MP_BBMAP,   MP_TMAP,   MP_HISAT2,   MP_BOWTIE,   MP_NOVOALIGN,   MP_RAZER3,    MP_BLASR,   MP_NGMLR,           MP_DELVE,   MP_TOPHAT,   MP_CPU,   MP_LONGRANGER,          MP_CLC,             NUM_MAPPERS } SamMapperType;
#define SAM_MAPPER_NAME             { "Unknown_mapper", "bsbolt",              "bwa",    "BWA",    "minimap2",    "STAR",    "bowtie2",    "dragen",     "gem3",          "gem2sam",      "bismark",    "bsseeker2",      "Winnowmap",    "baz2bam",     "BBMap",    "tmap",    "hisat2",    "Bowtie",    "NovoAlign",    "razers3",    "blasr",    "ngmlr",            "Delve",    "TopHat",    "cpu",    "longranger",           "CLCGenomicsWB",                }
#define SAM_MAPPER_SIGNATURE        { "Unknown_mapper", "PN:bwa	VN:BSB"/*\t*/, "PN:bwa", "PN:BWA", "PN:minimap2", "PN:STAR", "PN:bowtie2", "ID: DRAGEN", "PN:gem-mapper", "PN:gem-2-sam", "ID:Bismark", "PN:BS Seeker 2", "PN:Winnowmap", "PN:baz2bam",  "PN:BBMap", "ID:tmap", "PN:hisat2", "ID:Bowtie", "PN:novoalign", "PN:razers3", "ID:BLASR", "PN:nextgenmap-lr", "ID:Delve", "ID:TopHat", "PN:cpu", "PN:longranger.lariat", "PN:clcgenomicswb",             }

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
    #define UNK_QNANE_LEN 127
    #define NUM_COLLECTED_WORDS 6
    char unknown_flavor_qnames[NUM_COLLECTED_WORDS][UNK_QNANE_LEN+1];

    // SAM/BAM stuff
    STRl (std_cigar, 16);       // first CIGAR in the file - used in case all CIGARs in the file are the same
    int num_mapped;             // number of segconf reads that are mapped - defined as having (!FLAG.unmapped, RNAME, POS, CIGAR)
    bool sam_is_unmapped;       // false if there is at least one read in segconf with (!FLAG.unmapped, RNAME, POS, CIGAR)
    SamMapperType sam_mapper;       
    bool is_biobambam2_sort;    // PG records indicate running biobambam tag-creating programs    
    bool has_bqsr;              // PG records indicate running GATK ApplyBQSR
    bool NM_is_integer;         // true if NM is integer, false if it binary
    bool has_TLEN_non_zero;
    bool has_DP_before_PL;
    bool is_collated;           // Every QNAME appears in two or more consecutive lines
    bool evidence_of_collated;  // during segconf: at least a one pair of consecutive lines has the same QNAME
    bool is_sorted;             // every two consecutive lines that have the same RNAME, have non-decreasing POS
    bool is_paired;             // file has a least one read that is marked as "last" in FLAG
    bool evidence_of_sorted;    // during segconf: at least a one pair of consecutive lines has the same RNAME and increasing POS
    bool sam_multi_RG;          // evidence that file has more than one type of RG
    bool sam_bisulfite;         // this BAM file is of reads that have been treated with bisulfite to detect methylation
    bool sam_predict_meth_call; // ZIP/PIZ: true if segging SEQ should also predict the methylation call in vb->meth_call
    bool bs_strand_not_by_rev_comp; // true if vb->bisulfite_strand cannot be predicted by FLAG.rev_comp
    bool MD_NM_by_unconverted;  // in bisulfate data, we still calculate MD:Z and NM:i vs unconverted reference
    bool has_MD_or_NM;          // ZIP/PIZ: call sam_analyze_copied_SEQ for SEQ copied from prim unless no cigar, no seq or (PIZ) explicitly told not to.
    bool NM_after_MD;           // in all segconf lines that had both NM and MD, NM appeared after MD
    bool nM_after_MD;           // same, for nM:i
    bool pysam_qual;            // ZIP/PIZ: BAM missing QUAL is generated by old versions of pysam: first byte is 0xff, followed by 0s (SAM spec requires all bytes to be 0xff)
    bool sam_has_depn;          // ZIP: true if any line in the segconf sample had SAM_FLAG_SECONDARY or SAM_FLAG_SUPPLEMENTARY
    bool has_barcodes;          // ZIP: file uses barcodes
    bool star_solo;             // ZIP: using STARsolo or cellranger
    uint32_t AS_is_2ref_consumed; // ZIP/PIZ: AS value tends to be double ref_consumed (counter during segconf, bool during seg/piz)
    uint8_t MAPQ_value;         // used during segconf.running to calculate sam_mapq_has_single_value
    bool MAPQ_has_single_value; // all non-0 MAPQ have the same value
    msType sam_ms_type;         // ZIP/PIZ: type of ms:i 
    XgIncSType sam_XG_inc_S;    // Does XG include soft_clip[0]
    bool is_long_reads;
    SagType sag_type;           // Type of sag
    bool depn_CIGAR_can_have_H; // some DEPN CIGARs (of alignments with SA:Z) have H (set while segging MAIN)
    bool SA_CIGAR_can_have_H;   // some SA_CIGARs have H (set while segging MAIN)
    thool SA_HtoS;              // when a DEPN CIGAR has H, the corresponding SA_CIGAR has S
    bool sag_has_AS;            // sag store the AS:i values of prim lines that have them. Set if its beneficial to seg AS:i in depn lines against prim
    uint32_t sam_cigar_len;     // approx average CIGAR len rounded up (during running==true - total len)
    uint32_t sam_seq_len;       // ZIP/PIZ: approx average (SEQ.len+hard-clips) rounded to the nearest (during running - total len)
    uint32_t seq_len_to_cm;     // ZIP/PIZ: approx average of (seq_len/cm:i) (during running - cumulative)
    char CR_CB_seperator;       // ZIP: seperator within CR:Z and CB:Z fields
    bool abort_gencomp;         // ZIP: vb=1 found out that the file actually has no depn or no prim, so we stop sending lines to prim/depn
    bool has_cellranger;        // ZIP/PIZ: if TX:Z and/or AN:Z fields are present, they were generated by cellranger
    uint8_t n_CR_CB_CY_seps, n_BC_QT_seps, n_RX_seps;
    char BC_sep, RX_sep;
    MiniContainer CY_con, QT_con;
    SmallContainer CB_con;

    STRl(CY_con_snip, 64);
    STRl(QT_con_snip, 64);
    STRl(CB_con_snip, 64 + SMALL_CON_NITEMS * 16);     

    // SAM/BAM and FASTQ
    bool nontrivial_qual;       // true if we know that not all QUAL values are the same (as they are in newer PacBio files)

    // VCF stuff
    bool vcf_is_varscan;        // this VCF file was produced by VarScan
    bool vcf_is_gvcf;
    bool vcf_is_beagle;
    bool vcf_is_gwas;           // GWAS-VCF format: https://github.com/MRCIEU/gwas-vcf-specification
    bool vcf_illum_gtyping;     // tags from Illumina GenCall genotyping software
    bool vcf_infinium;
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
    int r1_or_r2;               // in case compression is WITHOUT --pair: our guess of whether this file is R1 or R2
    
    // FASTA stuff
    bool fasta_has_contigs;     // the sequences in this FASTA represent contigs (as opposed to reads) - in which case we have a FASTA_CONTIG dictionary and RANDOM_ACCESS
    SeqType seq_type;           // nucleotide or protein
    unsigned seq_type_counter;  // used for calculating seq_type 

    // GFF stuff
    int gff_version;
    bool has_embdedded_fasta;

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
extern rom sam_mapper_name (SamMapperType mp);

