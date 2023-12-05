// ------------------------------------------------------------------
//   segconf.c
//   Copyright (C) 2021-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "container.h"
#include "mutex.h"

// Documented range for users
#define MIN_VBLOCK_MEMORY  1    // in MB 
#define MAX_VBLOCK_MEMORY  1024 // up to v14 this was 2048. brought down bc when reconstructing translated, could exceed 2048 which causes all kinds of problems. also, there's no compression benefit beyond 512MB. 

// Range in developer options eg --vblock 10000B 
#define ABSOLUTE_MIN_VBLOCK_MEMORY ((uint64_t)1000) // in Bytes
#define ABSOLUTE_MAX_VBLOCK_MEMORY ((uint64_t)MAX_VBLOCK_MEMORY MB)

typedef enum __attribute__ ((__packed__)) { TECH_NONE=-1, TECH_ANY=-2, TECH_CONS=-3, TECH_UNKNOWN=0,   TECH_ILLUM, TECH_PACBIO, TECH_ONT,          TECH_454, TECH_MGI,   TECH_IONTORR, TECH_HELICOS, TECH_NCBI, TECH_ULTIMA, TECH_SINGLR, TECH_ELEMENT, TECH_ONSO, NUM_TECHS } SeqTech;
#define TECH_NAME                         {                                          "Unknown_tech",   "Illumina", "PacBio",    "Oxford_Nanopore", "454",    "MGI_Tech", "IonTorrent", "Helicos",    "NCBI",    "Ultima",    "Singular",  "Element",    "Onso"               }
#define TECH(x) (segconf.tech == TECH_##x)

typedef enum __attribute__ ((__packed__)) { SQT_UNKNOWN, SQT_NUKE, SQT_AMINO, SQT_NUKE_OR_AMINO } SeqType;

typedef enum __attribute__ ((__packed__)) { PL_mux_by_DP_TEST, PL_mux_by_DP_NO, PL_mux_by_DP_YES } PLMuxByDP;

typedef enum __attribute__ ((__packed__)) { ms_NONE, ms_BIOBAMBAM, ms_MINIMAP2 } msType; // type of SAM ms:i field 
#define ms_type_NAME { "None", "biobambam", "minimap2"}

typedef enum __attribute__ ((__packed__)) { DP_DEFAULT, by_AD, by_SDP } FormatDPMethod;

typedef enum __attribute__ ((__packed__)) { L3_UNKNOWN, L3_EMPTY, L3_COPY_LINE1, L3_NCBI, NUM_L3s } FastqLine3Type;

// SamMapperType is part of the file format and values should not be changed (new ones can be added)
typedef enum  {                       MP_UNKNOWN,       MP_BSBOLT,             MP_bwa,   MP_BWA,   MP_MINIMAP2,   MP_STAR,   MP_BOWTIE2,   MP_DRAGEN,    MP_GEM3,         MP_GEM2SAM,     MP_BISMARK,   MP_BSSEEKER2,     MP_WINNOWMAP,   MP_BAZ2BAM,    MP_BBMAP,   MP_TMAP,   MP_HISAT2,   MP_BOWTIE,   MP_NOVOALIGN,   MP_RAZER3,    MP_BLASR,   MP_NGMLR,           MP_DELVE,   MP_TOPHAT,   MP_CPU,   MP_LONGRANGER,          MP_CLC,              MP_PBMM2,   MP_CCS,  MP_SNAP,   MP_BWA_MEM2,   MP_PARABRICKS,     MP_ISAAC,   MP_ULTIMA, MP_TORRENT_BC,   MP_BIONANO,      NUM_MAPPERS } SamMapperType;
#define SAM_MAPPER_NAME             { "Unknown_mapper", "bsbolt",              "bwa",    "BWA",    "minimap2",    "STAR",    "bowtie2",    "dragen",     "gem3",          "gem2sam",      "bismark",    "bsseeker2",      "Winnowmap",    "baz2bam",     "BBMap",    "tmap",    "hisat2",    "Bowtie1",   "NovoAlign",    "razers3",    "blasr",    "ngmlr",            "Delve",    "TopHat",    "cpu",    "longranger",           "CLCGenomicsWB",    "pbmm2",    "ccs",    "snap",    "bwa-mem2",    "parabricks",      "iSAAC",    "Ultima",  "Torrent_BC",    "Bionano",                   }
#define SAM_MAPPER_SIGNATURE        { "Unknown_mapper", "PN:bwa	VN:BSB"/*\t*/, "PN:bwa", "PN:BWA", "PN:minimap2", "PN:STAR", "PN:bowtie2", "ID: DRAGEN", "PN:gem-mapper", "PN:gem-2-sam", "ID:Bismark", "PN:BS Seeker 2", "PN:Winnowmap", "PN:baz2bam",  "PN:BBMap", "ID:tmap", "PN:hisat2", "ID:Bowtie", "PN:novoalign", "PN:razers3", "ID:BLASR", "PN:nextgenmap-lr", "ID:Delve", "ID:TopHat", "PN:cpu", "PN:longranger.lariat", "PN:clcgenomicswb", "PN:pbmm2", "PN:ccs", "PN:SNAP", "PN:bwa-mem2", "PN:pbrun fq2bam", "PN:iSAAC", "ID:UA-",  "PN:BaseCaller", "ID:xmap_to_bam",            }   
#define MP(x) (segconf.sam_mapper == MP_##x)

#define MAX_SHORT_READ_LEN 2500

// also defined in sam.h
#define SAM_MAX_QNAME_LEN 255/*exc. \0*/     // In initial SAM specification verions, the max QNAME length was 255, and reduced to 254 in Aug 2015. We support 255 to support old SAM/BAM files too. BAM specifies 255 including \0 (so 254).

typedef struct { char s[SAM_MAX_QNAME_LEN+1]; } QnameStr;

typedef enum { QHT_QUAL, QHT_CONSENSUS, QHT_OQ, NUM_QHT } QualHistType;

// seg configuration set prior to starting to seg a file during segconfig_calculate or txtheader_zip_read_and_compress
typedef struct {

    // Seg parameters - general
    uint64_t vb_size;            // ZIP/PIZ: compression VBlock size in bytes (PIZ: passed in SectionHeaderGenozipHeader.vb_size)
    bool running;                // currently in segconf_calculate()
    int has[MAX_DICTS];          // for select did_i's, counts the numner of times this field was encountered during segconf.running
    uint32_t line_len;           // approx line len
    float b250_per_line[MAX_DICTS]; // b250.len / num_lines
    #define AT_LEAST(did_i) ((uint64_t)(10.0 + (segconf.b250_per_line[did_i] * (float)(vb->lines.len32))))
    bool disable_random_acccess; // random_access section is not to be outputted

    // qname characteristics (SAM/BAM, KRAKEN and FASTQ)
    QnameFlavor qname_flavor[NUM_QTYPES+1];     // 0-QNAME 1-QNAME2(FASTQ) 1=secondary flavor(SAM) 2=NCBI LINE3 (FASTQ)
    QnameFlavor deep_sam_qname_flavor[2];       // save for stats, in --deep (SAM's QNAME and QNAME2)
    QnameFlavorProp flav_prop[NUM_QTYPES];      // flavor properties
    bool sorted_by_qname[NUM_QTYPES];           // qnames appear in the file in a sorted order
    bool qname_flavor_rediscovered[NUM_QTYPES]; // true if flavor has been modified (used only during segconf.running)
    QnameStr qname_line0[NUM_QTYPES];           // qname of line_i=0 (by which flavor is determined) (nul-terminated)
    SeqTech tech;
    
    // FASTQ and FASTA
    bool multiseq;              // sequences in file are variants of each others 

    // SAM/BAM and FASTQ
    uint32_t longest_seq_len;   // length of the longest seq_len in the segconf data 
    DictId seq_len_dict_id;     // dict_id of one of the QNAME/QNAME2/LINE3/FASTQ_AUX contexts, which is expected to hold the seq_len for this read. 0 if there is no such item.
    #define UNK_QNANE_LEN 191
    #define NUM_COLLECTED_WORDS 6
    char unk_flav_qnames[NUM_QTYPES][NUM_COLLECTED_WORDS][UNK_QNANE_LEN+1]; // first 6 qnames if flavor is unknown
    uint8_t n_unk_flav_qnames[NUM_QTYPES];  // number of entries populated in unk_flav_qnames
    #define NUM_UNK_ID_CTXS 10
    #define UNK_ID_LEN 32
    char unk_ids_tag_name[NUM_UNK_ID_CTXS][MAX_TAG_LEN];
    char unk_ids[NUM_UNK_ID_CTXS][NUM_COLLECTED_WORDS][UNK_ID_LEN+1];
    
    #define NUM_QUAL_SCORES 94
    uint32_t qual_histo[NUM_QHT][NUM_QUAL_SCORES]; // histogram of qual values in segconf data for each QHT
    uint32_t qual_histo_acgt[4][NUM_QUAL_SCORES];  // histogram of qual values corresponding to A,C,G,T in SEQ (QHT_QUAL only)
    bool nontrivial_qual;       // true if we know that not all QUAL values are the same (as they are in newer PacBio files)
    bool has_agent_trimmer;     // has fields generated by Agilent AGeNT Trimmer
    bool is_pacbio_ccs;
    
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
    bool is_bwa;                // aligner used is based on bwa
    bool is_minimap2;           // aligner used is based on minimap2
    bool is_bowtie2;            // aligner used is based on bowtie2
    bool has_bwa_meth;           // file was treated with bwa-meth before and after alignment with bwa
    bool pacbio_subreads;       // this is a pacbio subreads file
    bool sam_has_SA_Z;
    thool sam_has_BWA_XA_Z;
    bool sam_has_BWA_XS_i;
    bool sam_has_XM_i_is_mismatches;
    bool sam_has_BWA_XT_A ;
    bool sam_has_BWA_XC_i;
    bool sam_has_BWA_X01_i;
    bool sam_has_bowtie2_YS_i;
    bool sam_has_bismark_XM_XG_XR;
    bool sam_has_ultima_t0;
    bool sam_has_zm_by_Q1NAME;
    bool sam_has_xcons;
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
    uint32_t AS_is_ref_consumed;// ZIP/PIZ: AS value tends to be double ref_consumed (counter during segconf, bool during seg/piz)
    uint8_t MAPQ_value;         // used during segconf.running to calculate sam_mapq_has_single_value
    bool MAPQ_has_single_value; // all non-0 MAPQ have the same value
    msType sam_ms_type;         // ZIP/PIZ: type of ms:i 
    thool sam_XG_inc_S;         // Does XG include soft_clip[0]
    bool is_long_reads;
    SagType sag_type;           // Type of sag
    bool depn_CIGAR_can_have_H; // some DEPN CIGARs (of alignments with SA:Z) have H (set while segging MAIN)
    bool SA_CIGAR_can_have_H;   // some SA_CIGARs have H (set while segging MAIN)
    thool SA_HtoS;              // when a DEPN CIGAR has H, the corresponding SA_CIGAR has S
    bool sag_has_AS;            // sag store the AS:i values of prim lines that have them. Set if its beneficial to seg AS:i in depn lines against prim
    uint32_t sam_cigar_len;     // approx average CIGAR len rounded up (during segconf.running==true - total len)
    uint32_t sam_seq_len;       // ZIP/PIZ: approx average (SEQ.len+hard-clips) rounded to the nearest (during segconf.running - total len)
    uint32_t seq_len_to_cm;     // ZIP/PIZ: approx average of (seq_len/cm:i) (during segconf.running - cumulative)
    uint32_t sam_cropped_at;    // ZIP: it appears that the FASTQ reads were sam_is_cropped to an identical length before aligning. If --deep, this field (as calculated in SAM) is also used when segging FASTQ.
    bool use_pacbio_iqsqdq;     // ZIP: if iq:Z sq:Z and dq:Z are present in all lines, we can compress them, as well as QUAL, better. 
    char CR_CB_seperator;       // ZIP: seperator within CR:Z and CB:Z fields
    bool abort_gencomp;         // ZIP: vb=1 found out that the file actually has no depn or no prim, so we stop sending lines to prim/depn
    bool has_cellranger;        // ZIP/PIZ: if TX:Z and/or AN:Z fields are present, they were generated by cellranger
    bool has_RSEM;              // ZIP: RSEM is used (https://github.com/bli25/RSEM_tutorial)
    uint8_t n_CR_CB_CY_seps, n_BC_QT_seps, n_RX_seps;
    char BC_sep, RX_sep;
    MiniContainer CY_con, QT_con;
    SmallContainer CB_con, MM_con;
    char sam_deep_filename[1024]; // name of the SAM filename in --deep 
    uint32_t est_segconf_sam_size;// BAM only: est size of segconf vblock translated to SAM
    #define SAM_FACTOR_MULT 32  // multiplication for storing est_sam_factor as an int in SectionHeaderGenozipHeader
    double est_sam_factor;      // BAM only: est size factor when translating BAM to SAM based on segconf data

    STRl(CY_con_snip, 64);
    STRl(QT_con_snip, 64);
    STRl(CB_con_snip, 64 + SMALL_CON_NITEMS * 16);     
    STRl(MM_con_snip, 64 + SMALL_CON_NITEMS * 16);     

    // VCF stuff
    bool vcf_is_varscan;        // this VCF file was produced by VarScan
    bool vcf_is_gvcf;
    bool vcf_is_beagle;
    bool vcf_is_dragen;
    bool vcf_is_manta;
    bool vcf_is_cosmic;
    bool vcf_is_clinvar;    
    bool vcf_is_pindel;
    bool vcf_is_caveman;
    bool vcf_is_vagrent;
    bool vcf_is_gwas;           // GWAS-VCF format: https://github.com/MRCIEU/gwas-vcf-specification
    bool vcf_illum_gtyping;     // tags from Illumina GenCall genotyping software
    bool vcf_is_infinium;
    bool vcf_is_dbSNP;
    bool vcf_is_giab_trio;
    bool vcf_is_vep;
    char *vcf_vep_spec;
    bool vcf_is_gnomad;
    bool vcf_is_icgc;
    bool vcf_is_exac;
    bool vcf_is_mastermind;
    bool vcf_is_isaac;          // IsaacVariantCaller / starling
    bool vcf_is_deep_variant;   // Google Deep Variant
    bool vcf_is_ultima_dv;      // Ultima Genomics version of Deep Variant
    uint64_t count_dosage[2];   // used to calculate pc_has_dosage
    float pc_has_dosage;        // % of the samples x lines that have a valid (0-2) dosage value [0.0,1.0]
    bool use_null_DP_method;    // A method for predicting GT=./. by DP=.
    FormatDPMethod FORMAT_DP_method;
    PLMuxByDP PL_mux_by_DP;
    Mutex PL_mux_by_DP_mutex;
    uint64_t count_GQ_by_PL, count_GQ_by_GP; // used tp calculate GQ_by_PL, GQ_by_GP
    bool GQ_by_PL, GQ_by_GP;
    
    // FASTQ
    union {
        struct {            
            uint8_t has_desc     : 4; // non-zero if any of the 4 bitfields are set
            uint8_t unused       : 4;
        };
        struct {
            uint8_t has_qname2   : 1;               
            uint8_t has_extra    : 1;
            uint8_t has_aux      : 1; // aux data in the format "length=7 pooptiz=2"
            uint8_t has_saux     : 1; // aux data in SAM format "BC:Z:TATTCATA+TCCAAGCG        ZX:Z:TTAA"
            uint8_t saux_tab_sep : 1; // has_saux AND seperator before SAUX is '\t'
            uint8_t desc_is_l3   : 1; // either L1 or L3 can have these properties, but not both
            uint8_t unused2      : 2;
        };
    };

    FastqLine3Type line3;       // format of line3
    int r1_or_r2;               // in case compression is WITHOUT --pair: our guess of whether this file is R1 or R2
    QType deep_qtype;           // Deep: QNAME1 or QNAME2 if for most segconf lines for which have Deep, SAM qname matches FASTQ's QNAME1 or QNAME2, or QNONE if not
    bool deep_no_qual;          // Deep: true if for most segconf lines which have Deep, qual doesn't match (eg, bc of undocumented BQSR) 
    bool deep_has_trimmed;      // Deep: some FASTQ reads in segconf appear in SAM trimmed (beyond cropping)
    bool deep_has_trimmed_left; // Deep: some FASTQ reads in segconf are trimmed on the left too (not just the right)
    char deep_1st_desc[256];    // Deep: DESC line of first FASTQ read of first FASTQ file
    char deep_1st_qname[256];   // Deep: QNAME of first SAM alignment in file
    unsigned n_full_mch[2];     // Deep: count segconf lines where hash matches with at least one SAM line - (QNAME1 or QNAME2), SEQ, QUAL
    unsigned n_seq_qual_mch;    // Deep: count segconf lines where hash matches with at least one SAM line - SEQ and QUAL
    unsigned n_seq_qname_mch[2];// Deep: count segconf lines where FASTQ (QNAME1 or QNAME2) hash matches with at least one SAM line - SEQ and QNAME 
    unsigned n_seq_mch;         // Deep: count segconf lines where hash matches with at least one SAM line - SEQ
    unsigned n_no_mch;          // Deep: count segconf lines that don't match any SAM line (perhaps because SAM is filtered)

    // FASTA stuff (including FASTA embedded in GFF)
    bool fasta_has_contigs;     // the sequences in this FASTA represent contigs (as opposed to reads) - in which case we have a FASTA_CONTIG dictionary and RANDOM_ACCESS
    SeqType seq_type;           // nucleotide or protein
    unsigned seq_type_counter;  // used for calculating seq_type 

    // GFF stuff
    int gff_version;
    bool has_embedded_fasta;

    // Chain stuff
    bool chain_mismatches_ref;  // Some contigs mismatch the reference files, so this chain file cannot be used with --chain

    // BED stuff
    int bed_num_flds;
} SegConf;

extern SegConf segconf; // ZIP: set based on segging a sample of a few first lines of the file
                        // PIZ: select fields are transferred through SectionHeaderGenozipHeader

extern void segconf_initialize (void);
extern void segconf_calculate (void);
extern bool segconf_is_long_reads(void);
extern void segconf_mark_as_used (VBlockP vb, unsigned num_ctxs, ...);
extern rom segconf_sam_mapper_name (void);
extern rom segconf_tech_name (void);
extern rom segconf_deep_trimming_name (void);
extern void segconf_test_sorted (VBlockP vb, WordIndex prev_line_chrom, PosType32 pos, PosType32 prev_line_pos);
extern StrText segconf_get_qual_histo (QualHistType qht);
