// ------------------------------------------------------------------
//   segconf.h
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "container.h"
#include "mutex.h"
#include "data_types.h"
#include "dict_id_gen.h"

// Documented range for users
#define MIN_VBLOCK_MEMORY  1    // in MB 
#define MAX_VBLOCK_MEMORY  1024 // up to v14 this was 2048. brought down bc when reconstructing translated, could exceed 2048 which causes all kinds of problems. also, there's no compression benefit beyond 512MB. 

// Range in developer options eg --vblock 10000B 
#define ABSOLUTE_MIN_VBLOCK_MEMORY ((uint64_t)1000) // in Bytes
#define ABSOLUTE_MAX_VBLOCK_MEMORY ((uint64_t)MAX_VBLOCK_MEMORY MB)

#define MAX_SEGCONF_LINES 1000 // max lines tested in segconf (even if VB is large e.g. due to reading full MGZIP block)

typedef packed_enum { TECH_NONE=-1, TECH_ANY=-2, TECH_CONS=-3, TECH_NCBI=-4, TECH_UNKNOWN=0,   TECH_ILLUMINA, TECH_PACBIO, TECH_NANOPORE,     TECH_LS454, TECH_MGI,   TECH_IONTORR, TECH_HELICOS, TECH_ULTIMA, TECH_SINGULAR, TECH_ELEMENT, TECH_ONSO, TECH_CAPILLARY, TECH_SOLID, TECH_SIKUN, NUM_TECHS } SeqTech;
#define TECH_NAME   {                                                            "Unknown_tech",   "Illumina",    "PacBio",    "Oxford_Nanopore", "LS454",    "MGI_Tech", "IonTorrent", "Helicos",    "Ultima",    "Singular",    "Element",    "Onso",    "Capillary",    "SOLiD",    "Sikun"               }
#define TECH(x) (segconf.tech == TECH_##x)

typedef packed_enum { SQT_UNKNOWN, SQT_NUKE, SQT_AMINO } SeqType;

typedef packed_enum { ms_NONE, ms_BIOBAMBAM, ms_MINIMAP2 } msType; // type of SAM ms:i field 
#define ms_type_NAME { "None", "biobambam", "minimap2"}

typedef packed_enum { FMT_DP_DEFAULT=0, BY_AD=1, BY_SDP=2, BY_INFO_DP=3 } FormatDPMethod; // part of file format: values go into SectionHeaderGenozipHeader.segconf_FMT_DP_method

typedef packed_enum { BY_FORMAT_DP=0, INFO_DP_DEFAULT=1, BY_BaseCounts=2 } InfoDPMethod; // part of file format: values go into SectionHeaderGenozipHeader.segconf_INF_DP_method since 15.0.52

typedef packed_enum { VCF_QUAL_DEFAULT, VCF_QUAL_local, VCF_QUAL_by_RGQ, VCF_QUAL_mated, VCF_QUAL_by_GP } VcfQualMethod;

typedef packed_enum { VCF_INFO_DEFAULT, VCF_INFO_by_RGQ, VCF_INFO_by_FILTER } VcfInfoMethod;

typedef packed_enum { L3_UNKNOWN, L3_EMPTY, L3_COPY_LINE1, L3_NCBI } FastqLine3Type;

typedef packed_enum { INFO_VT_UNKNOWN, INFO_VT_VAGrENT, INFO_VT_1KG, INFO_VT_CALLMOM } InfoVTType; // part of the file format: values go into the snip of VCF_SPECIAL_VT

typedef packed_enum { RG_DEFAULT, RG_BY_ILLUM_QNAME, RG_BY_NCBI_QNAME } RGMethod; // part of the file format: values go into the snip of SAM_SPECIAL_RG_by_QNAME

// SamMapperType is part of the file format and values should not be changed (new ones can be added)
typedef packed_enum  {                 MP_UNKNOWN,       MP_BSBOLT,             MP_bwa,   MP_BWA,   MP_MINIMAP2,   MP_STAR,   MP_BOWTIE2,   MP_DRAGEN,    MP_GEM3,         MP_GEM2SAM,     MP_BISMARK,   MP_BSSEEKER2,     MP_WINNOWMAP,   MP_BAZ2BAM,    MP_BBMAP,   MP_TMAP,   MP_HISAT2,   MP_BOWTIE,   MP_NOVOALIGN,   MP_RAZER3,    MP_BLASR,   MP_NGMLR,           MP_DELVE,   MP_TOPHAT,   MP_CPU,   MP_LONGRANGER,          MP_CLC,              MP_PBMM2,   MP_CCS,  MP_SNAP,   MP_BWA_MEM2,   MP_PARABRICKS,     MP_ISAAC,   MP_ULTIMA, MP_TORRENT_BC,   MP_BIONANO,       MP_CRDNA,        MP_VG,   MP_CRATAC,            MP_CELLRANGER, NUM_MAPPERS } SamMapperType;
#define SAM_MAPPER_NAME             { "Unknown_mapper", "bsbolt",              "bwa",    "BWA",    "minimap2",    "STAR",    "bowtie2",    "dragen",     "gem3",          "gem2sam",      "bismark",    "bsseeker2",      "Winnowmap",    "baz2bam",     "BBMap",    "tmap",    "hisat2",    "Bowtie1",   "NovoAlign",    "razers3",    "blasr",    "ngmlr",            "Delve",    "TopHat",    "cpu",    "longranger",           "CLCGenomicsWB",    "pbmm2",    "ccs",    "snap",    "bwa-mem2",    "parabricks",      "iSAAC",    "Ultima",  "Torrent_BC",    "Bionano",        "CellRangerDNA", "vg",    "CellRangerATAC",     "CellRanger",              }
#define SAM_MAPPER_SIGNATURE        { "Unknown_mapper", "PN:bwa	VN:BSB"/*\t*/, "PN:bwa", "PN:BWA", "PN:minimap2", "PN:STAR", "PN:bowtie2", "ID: DRAGEN", "PN:gem-mapper", "PN:gem-2-sam", "ID:Bismark", "PN:BS Seeker 2", "PN:Winnowmap", "PN:baz2bam",  "PN:BBMap", "ID:tmap", "PN:hisat2", "ID:Bowtie", "PN:novoalign", "PN:razers3", "ID:BLASR", "PN:nextgenmap-lr", "ID:Delve", "ID:TopHat", "PN:cpu", "PN:longranger.lariat", "PN:clcgenomicswb", "PN:pbmm2", "PN:ccs", "PN:SNAP", "PN:bwa-mem2", "PN:pbrun fq2bam", "PN:iSAAC", "ID:UA-",  "PN:BaseCaller", "ID:xmap_to_bam", "ID:crdna",      "PN:vg", "VN:cellranger-atac", "ID:cellranger",           }   
#define MP(x) (segconf.sam_mapper == MP_##x)

#define MAX_SHORT_READ_LEN 1500

// also defined in sam.h
#define SAM_MAX_QNAME_LEN 255/*exc. \0*/     // In initial SAM specification verions, the max QNAME length was 255, and reduced to 254 in Aug 2015. We support 255 to support old SAM/BAM files too. BAM specifies 255 including \0 (so 254).

typedef struct { char s[SAM_MAX_QNAME_LEN+1]; } QnameStr;

typedef enum { QHT_QUAL, QHT_OQ, NUM_QHT } QualHistType;
extern QualHistType did_i_to_qht (Did did_i);

typedef struct { uint8_t q; int count; } QualHisto;

typedef packed_enum { GQ_old=0/*up to 15.0.36*/, BY_PL=1, BY_GP=2, MUX_DOSAGExDP, MUX_DOSAGE, GQ_INTEGER } GQMethodType; // values go into SectionHeaderGenozipHeader.segconf_GQ_method (only 0,1,2 are used in PIZ)

typedef packed_enum { RO_AO_none, RO_AO_by_AD, RO_AO_by_DP } ROAOMethodType;

typedef packed_enum { GP_unknown, GP_probabilities, GP_phred, GP_other } GPContentType;

typedef packed_enum { TRIM_IS_M, TRIM_IS_S, TRIM_IS_MS/*extends existing M or S*/ } BamAssTrimCigarTreatment;
#define TRIM_IS_NAMES { "M",     "S",       "MS" }

typedef packed_enum         { MATE_NONE, MATE_01, MATE_12, MATE_PBSV } MateIDMethodType;
#define MATEID_METHOD_NAMES { "NONE",    "01",    "12",    "PBSV"    }

// seg configuration set prior to starting to seg a file during segconfig_calculate or txtheader_zip_read_and_compress
// Notes:
// - gcc_struct causes mingw struct to be the same as Linux. Specifically, bit-fields with a different type DO NOT start a new aligned field.
// - segconf is NOT ((packed)), but we take care to order members so that it is as packed as possible
#ifdef __clang__
typedef struct { 
#else //gcc
typedef struct __attribute__((gcc_struct)) { 
#endif
    // ________________________________________________________________________________________
    // For each field, there exist files where field is accessed in every line of NORMAL 
    // (no special flags) SEG or RECON. Keep close together + use bits for L1 cache effeciency
    // ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾ 
    // Seg parameters - general (frequent)
    struct { uint16_t b250, local; } per_line[MAX_NUM_PREDEFINED]; // b250.len (or local.len) / num_lines

    // qname characteristics (SAM/BAM and FASTQ) (frequent)
    QnameFlavor qname_flavor[NUM_QTYPES+1];     // 0-QNAME 1-QNAME2(FASTQ) 1=secondary flavor(SAM) 2=NCBI LINE3 (FASTQ)
    QnameFlavorProp flav_prop[NUM_QTYPES];      // ZIP: flavor properties (in PIZ: this is in z_file->flav_prop)
    bool qname_flavor_rediscovered[NUM_QTYPES]; // true if flavor has been modified (used only during segconf.running)
    bool sorted_by_qname[NUM_QTYPES];    // qnames appear in the file in a sorted order
    
    // SAM/BAM and FASTQ (frequent)
    uint32_t std_seq_len;                // length of the longest seq_len in the segconf data. For FASTQ: if paired, applies to R1
    uint32_t std_seq_lR2;                // FASTQ R2: std_seq_len
    DictId seq_len_dict_id;              // dict_id of one of the QNAME/QNAME2/LINE3/FASTQ_AUX contexts, which is expected to hold the seq_len for this read. 0 if there is no such item.

    // SAM/BAM stuff (frequent)
    STRl (std_cigar, 16);                // first CIGAR in the file - used in case all CIGARs in the file are the same
    SamMapperType sam_mapper;    

    // fields accessed atomically, so can't be bits (frequent)
    thool sam_has_BWA_XA_Z;              // not a bit field because we access it atomically
    bool pysam_qual;                     // ZIP/PIZ: BAM missing QUAL is generated by old versions of pysam: first byte is 0xff, followed by 0s (SAM spec requires all bytes to be 0xff) (not bitfield because accessed atomicallys)
    thool SA_HtoS;                       // when a DEPN CIGAR has H, the corresponding SA_CIGAR has S

    // bit fields applicable to multiple data types (frequent)
    #define segconf_has(did_i) bitset_get (segconf.has_bits, (did_i))
    const uint64_t has_bits[(MAX_DICTS+63)/64]; // bit array: true if did_i was encountered during segconf

    uint64_t running                : 1; // currently in segconf_calculate() (first bit: no shift)
    SeqTech tech                    : 6; // sequencer technology
    uint64_t disable_random_access  : 1; // random_access section is not to be outputted (SAM, VCF etc)
    uint64_t nontrivial_qual        : 1; // true if we know that not all QUAL values are the same (as they are in newer PacBio files)
    uint64_t has_agent_trimmer      : 1; // has fields generated by Agilent AGeNT Trimmer
    uint64_t is_pacbio_ccs          : 1; // not yet used

    // SAM/BAM bit fields (frequent)
    uint64_t sam_is_unmapped        : 1; // false if there is at least one read in segconf with (!FLAG.unmapped, RNAME, POS, CIGAR)
    SagType sag_type                : 3; // Type of sag
    uint64_t is_biobambam2_sort     : 1; // PG records indicate running biobambam tag-creating programs    
    uint64_t has_bqsr               : 1; // PG records indicate running GATK ApplyBQSR
    uint64_t NM_is_integer          : 1; // true if NM is integer, false if it binary
    uint64_t has_TLEN_non_zero      : 1;
    uint64_t has_DP_before_PL       : 1;
    uint64_t is_collated            : 1; // Every QNAME appears in two or more consecutive lines
    uint64_t is_sorted              : 1; // every two consecutive lines that have the same RNAME, have non-decreasing POS
    uint64_t is_paired              : 1; // file has a least one read that is marked as "last" in FLAG
    uint64_t sam_multi_RG           : 1; // file probably has more than one type of RG
    uint64_t is_bwa                 : 1; // aligner used is based on bwa
    uint64_t is_minimap2            : 1; // aligner used is based on minimap2
    uint64_t is_bowtie2             : 1; // aligner used is based on bowtie2
    uint64_t has_bwa_meth           : 1; // file was treated with bwa-meth before and after alignment with bwa
    uint64_t pacbio_subreads        : 1; // this is a pacbio subreads file
    uint64_t use_insertion_ctxs     : 1; // use separate contexts for SEQ insertions
    uint64_t sam_use_sn_mux         : 1;
    uint64_t sam_semcol_in_contig   : 1; // some contig names contain a semicolon, eg "ScxkALA_1850;HRSCAF=2697" (see: https://hgdownload.soe.ucsc.edu/hubs/GCF/005/870/125/GCF_005870125.1/GCF_005870125.1.chromAlias.txtwcs)
    uint64_t sam_has_SA_Z           : 1; // file is expected to have SA:Z based on the mapper (not necessary discovered by segconf)
    uint64_t sam_has_BWA_XS_i       : 1;
    uint64_t sam_has_XM_i_is_mismatches : 1;
    uint64_t sam_has_BWA_XT_A       : 1;
    uint64_t sam_has_BWA_XC_i       : 1;
    uint64_t sam_has_BWA_X01_i      : 1;
    uint64_t sam_has_bowtie2_YS_i   : 1;
    uint64_t sam_has_bismark_XM_XG_XR : 1;
    uint64_t sam_has_ultima_t0      : 1;
    uint64_t sam_has_zm_by_Q1NAME   : 1;
    uint64_t sam_is_nanoseq         : 1;
    uint64_t sam_has_abra2          : 1;
    uint64_t sam_diverse_qs         : 1; // true if not all qs:i values in segconf are equal
    uint64_t sam_bisulfite          : 1; // this BAM file is of reads that have been treated with bisulfite to detect methylation
    uint64_t sam_predict_meth_call  : 1; // ZIP/PIZ: true if segging SEQ should also predict the methylation call in vb->meth_call
    uint64_t bs_strand_not_by_rev_comp : 1; // true if vb->bisulfite_strand cannot be predicted by FLAG.rev_comp
    uint64_t MD_NM_by_unconverted   : 1; // in bisulfate data, we still calculate MD:Z and NM:i vs unconverted reference
    uint64_t has_MD_or_NM           : 1; // ZIP/PIZ: call sam_analyze_copied_SEQ for SEQ copied from prim unless no cigar, no seq or (PIZ) explicitly told not to.
    uint64_t NM_after_MD            : 1; // in all segconf lines that had both NM and MD, NM appeared after MD
    uint64_t nM_after_MD            : 1; // same, for nM:i
    uint64_t sam_has_depn           : 1; // ZIP: true if any line in the segconf sample had SAM_FLAG_SECONDARY or SAM_FLAG_SUPPLEMENTARY
    uint64_t has_barcodes           : 1; // ZIP: file uses barcodes
    uint64_t star_solo              : 1; // ZIP: using STARsolo or cellranger
    uint64_t AS_is_2ref_consumed    : 1; // ZIP/PIZ: AS value tends to be double ref_consumed (counter during segconf, uint64_t during seg/piz)
    uint64_t AS_is_ref_consumed     : 1; // ZIP/PIZ: AS value tends to be double ref_consumed (counter during segconf, uint64_t during seg/piz)
    uint64_t MAPQ_has_single_value  : 1; // all non-0 MAPQ have the same value
    uint64_t MAPQ_use_xq            : 1; // ZIP: DRAGEN: use the exists(xq:i) for segging MAPQ     
    msType sam_ms_type              : 2; // ZIP/PIZ: type of ms:i 
    thool sam_XG_inc_S              : 2; // Does XG include soft_clip[0]
    uint64_t is_long_reads          : 1;
    uint64_t HI_has_two_plus        : 1; // a HI:i value of >= 2 was detected in segconf 
    uint64_t qual_in_mem_use_rans   : 1; // Valid only if huffman for QUAL is missing: If true, use RANS for in-memory compression of QUAL (gencomp/deep) ; if false use ARITH
    uint64_t CIGAR_has_eqx          : 1; // segconf detected lines with '=' and/or 'X' in their CIGAR (always in pbmm2, with --eqx in minimap2)
    uint64_t SA_NM_by_CIGAR_X       : 1; // NM:i is not used, instead we get SA_NM it from the number of X bases in CIGAR
    uint64_t depn_CIGAR_can_have_H  : 1; // some DEPN CIGARs (of alignments with SA:Z) have H (set while segging MAIN)
    uint64_t SA_CIGAR_can_have_H    : 1; // some SA_CIGARs have H (set while segging MAIN)
    thool SA_CIGAR_abbreviated      : 2; // CIGAR strings in SA:Z are abbreviated as in minimap2, see https://github.com/lh3/minimap2/blob/master/format.c : mm_write_sam3
    uint64_t sag_has_AS             : 1; // sag store the AS:i values of prim lines that have them. Set if its beneficial to seg AS:i in depn lines against prim
    uint64_t use_pacbio_iqsqdq      : 1; // ZIP: if iq:Z sq:Z and dq:Z are present in all lines, we can compress them, as well as QUAL, better. 
    uint64_t no_gc_checking         : 1; // ZIP: if true: vb=1 had no depn lines, so we heuristically decided not to check for gc in future VBs (note that some VBs running in parallel to vb=1 might have had depn lines)
    uint64_t has_10xGen             : 1; // SAM: ZIP/PIZ: has 10xGenomics tags 
    uint64_t has_Parse              : 1; // ZIP: has Parse Biosciences fields
    uint64_t has_TR_TQ              : 1; // ZIP: use cellrangerATAC-style TR:Z / TQ:Z methods             
    uint64_t has_RSEM               : 1; // ZIP: RSEM is used (https://github.com/bli25/RSEM_tutorial)
    RGMethod RG_method              : 2;
    uint64_t sam_has_xcons          : 1;
    uint64_t xcons_std_seq_len      : 9; // Seg / PIZ: the most common standard xcons qual length ("standard": alignments with XC:i but without XO:i) (values: 100 to 355)

    #define FAF segconf.fasta_as_fastq
    uint64_t fasta_as_fastq         : 1; // ZIP/PIZ: Segging a FASTA file as a QUAL-less FASTQ (also in SAM with Deep)
    
    uint16_t s1_to_cm_32;                // ZIP/PIZ: approx average of (s1:i/cm:i) x 32 (during segconf.running - cumulative)
    uint8_t seq_len_to_cm;               // ZIP/PIZ: approx average of (seq_len/cm:i) (during segconf.running - cumulative)

    uint8_t n_CR_CB_CY_seps, n_BC_QT_seps, n_RX_seps;
    char BC_sep, RX_sep, CR_CB_seperator;
    
    // each container and its snip are accessed together - store them in proximity (frequent)
    #define MAX_CB_ITEMS 3
    Container(MAX_CB_ITEMS) CB_con; STRl(CB_con_snip, con_snip_sizeof(MAX_CB_ITEMS));     
    Container(1) CY_con;            STRl(CY_con_snip, con_snip_sizeof(1)); 
    Container(1) QT_con;            STRl(QT_con_snip, con_snip_sizeof(1));
    Container(2) MM_con;            STRl(MM_con_snip, con_snip_sizeof(2));     

    // FASTQ / FASTA (frequent) (packed)
    char aux_sep;                        // separator between name and value in aux fields (either '=' or ':')
    char desc_char;                      // FASTQ/FASTA: FASTA: '>' by default, but can also be '@' ; FASTQ: normally '@', but can be '>' if fasta_as_fastq
    
    #define DC segconf.desc_char
/*⸨*/uint8_t deep_paired_qname      : 1; // Deep: QNAME hash includes deep_is_last
    thool deep_is_last              : 2; // Deep: whether this FASTQ file corresponds to is_first or is_last alignments in BAM, or unknown if !segconf.deep_paired_qname 
    uint8_t deep_no_qual            : 1; // Deep: true if for most segconf lines which have Deep, qual doesn't match (eg, bc of undocumented BQSR) 
    uint8_t deep_has_trimmed        : 1; // Deep: some FASTQ reads in segconf appear in SAM trimmed (beyond cropping)
    uint8_t deep_has_trimmed_left   : 1; // Deep: some FASTQ reads in segconf are trimmed on the left too (not just the right)
    int8_t deep_qtype               : 2; // Deep ZIP/PIZ: QNONE, QNAME1 or QNAME2 if for most segconf lines for which have Deep, SAM qname matches FASTQ's QNAME1 or QNAME2. QNONE means no deep, or in PIZ of Deep files up to 15.0.66, it means matching was by SEQ and QUAL only (not QNAME) 
    char deep_N_sam_score;               // Deep: ZIP only: Base qualities of 'N' bases in the SAM are this value, regardless of their value in FASTQ
    char deep_N_fq_score;                // Deep: ZIP/PIZ:  Base qualities of 'N' bases in the FASTQ are this value
/*⸩*/
/*⸨*/FastqLine3Type line3           : 2; // format of line3
    thool is_interleaved            : 2; // whether FASTQ file is identified as interleaved

    // FASTA stuff (including FASTA embedded in GFF)
    SeqType fasta_seq_type          : 2; // nucleotide or protein
    uint8_t fasta_has_contigs       : 1; // the sequences in this FASTA represent contigs (as opposed to reads) - in which case we have a FASTA_CONTIG dictionary and RANDOM_ACCESS
/*⸩*/
    enum { NONBIO_NONE, NONBIO_Parse/*R2 is non-biological*/, NONBIO_10xGen/*R1 is non-biological*/ } nonbio_type : 2;
    #define IS_NONBIO(x) (segconf.nonbio_type == NONBIO_##x)
    
    union { // one byte
        struct {            
            uint8_t has_desc        : 4; // non-zero if any of the 4 bitfields are set
            uint8_t unused          : 4;
        };
        struct {
            uint8_t has_qname2      : 1;               
            uint8_t has_extra       : 1;
            uint8_t has_aux         : 1; // aux data in the format "length=7 pooptiz=2"
            uint8_t has_saux        : 1; // aux data in SAM format "BC:Z:TATTCATA+TCCAAGCG        ZX:Z:TTAA"
            uint8_t saux_tab_sep    : 1; // has_saux AND seperator before SAUX is '\t'
            uint8_t desc_is_l3      : 1; // either L1 or L3 can have these properties, but not both
            uint8_t unused2         : 2;
        };
    };

    // FASTQ - non-biological file (if paired - this refers to R2) (frequent)
    struct split_seq {
        Did umi_did_i;                   // optional: qname item that is expected to be identical to the UMI (first 10 bases of sequence).
        uint8_t non_bio_len, linker_index[2], linker_len[2], barcode_index[3]; 
    } split_seq;
    
    STRl(nonbio_con_snip, con_snip_sizeof(7));     
    STRl(copy_qname_umi_snip, 32);

    // VCF (frequent)
/*⸨*/uint64_t vcf_is_varscan        : 1; // this VCF file was produced by VarScan
    uint64_t vcf_is_gvcf            : 1;
    uint64_t vcf_is_gatk_gvcf       : 1;
    uint64_t vcf_is_beagle          : 1;
    uint64_t vcf_is_dragen          : 1;
    uint64_t vcf_is_hail            : 1;
    uint64_t vcf_is_manta           : 1;
    uint64_t vcf_is_cosmic          : 1;
    uint64_t vcf_is_clinvar         : 1;    
    uint64_t vcf_is_pindel          : 1;
    uint64_t vcf_is_caveman         : 1;
    uint64_t vcf_is_vagrent         : 1;
    uint64_t vcf_is_platypus        : 1;
    uint64_t vcf_is_gwas            : 1; // GWAS-VCF format: https://github.com/MRCIEU/gwas-vcf-specification
    uint64_t vcf_illum_gtyping      : 1; // tags from Illumina GenCall genotyping software
    uint64_t vcf_is_infinium        : 1;
    uint64_t vcf_is_dbSNP           : 1;
    uint64_t vcf_is_giab            : 1;
    uint64_t vcf_is_giab_trio       : 1;
    uint64_t vcf_is_vep             : 1;
    uint64_t vcf_is_gnomad          : 1;
    uint64_t vcf_is_icgc            : 1;
    uint64_t vcf_is_exac            : 1;
    uint64_t vcf_is_mastermind      : 1;
    uint64_t vcf_is_isaac           : 1; // IsaacVariantCaller / starling
    uint64_t vcf_is_deep_variant    : 1; // Google Deep Variant
    uint64_t vcf_is_ultima          : 1; // Ultima Genomics version of Deep Variant
    uint64_t vcf_is_svaba           : 1;
    uint64_t vcf_is_pbsv            : 1;
    uint64_t vcf_is_sv              : 1;
    uint64_t vcf_is_callmom         : 1;
    uint64_t vcf_is_melt            : 1;
    uint64_t vcf_is_GLIMPSE_phase   : 1;
    uint64_t vcf_is_gencove         : 1;
    uint64_t vcf_is_freebayes       : 1; // generated by https://github.com/freebayes/freebayes
    uint64_t vcf_is_giggle          : 1;
    uint64_t use_null_DP_method     : 1; // A method for predicting GT=./. by DP=.
    uint64_t vcf_del_svlen_is_neg   : 1;
    uint64_t vcf_local_alleles      : 1;
    uint64_t vcf_sample_copy        : 1; // possibly copy this sample on previous line
    GPContentType FMT_GP_content    : 2; 
    FormatDPMethod FMT_DP_method    : 2;
    InfoDPMethod INFO_DP_method     : 2;
    thool PL_mux_by_DP              : 2;
    uint64_t FI_by_DP               : 1;
    uint64_t AS_SB_TABLE_by_SB      : 1;
    InfoVTType INFO_VT_type         : 2;
    GQMethodType FMT_GQ_method      : 3; // values go into SectionHeaderGenozipHeader.segconf_GQ_method (only 0,1,2 are used in PIZ)
    MateIDMethodType MATEID_method  : 2; // method to convert between VCF_ID and the BND mate's VCF_ID
    ROAOMethodType FMT_RO_AO_method : 2;
    VcfInfoMethod vcf_INFO_method   : 2;
/*⸩*/     
    float Q_to_O;                        // freebayes: average ratio of QR/RO and QA/AO rounded (running: sum of)
/*⸨*/uint8_t vcf_QUAL_truncate_trailing_zeros : 1;
    VcfQualMethod vcf_QUAL_method   : 3;
    #define MAX_VCF_QUAL_DECIMALS 15
    uint8_t vcf_QUAL_decimals       : 4; // decimals in QUAL value
/*⸩*/ 
    uint8_t vcf_max_MAPQ;                // maximum MAPQ of BAM alignments that contributed to this variant, as derived from RAW_MQandDP, but not more than 223
    char vcf_ID_is_variant;              // ID format is eg 1_2704352_AT_A or 1:2704352:AT:A which is CHROM_POS_REF_ALT. value is '_' or ':' or 0. 
        
    // GFF (frequent)
    uint8_t gff_version;

    // BED (frequent)
    uint8_t bed_num_flds;
    
    // ________________________________________________________________________________________
    // Used for every VB, or very rarely in seg, but not used in every line of NORMAL (no special flags) SEG or RECON
    // ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾ 
    // General (every VB)
    QnameStr qname_line0[NUM_QTYPES];    // qname of line_i=0 (by which flavor is determined) (nul-terminated). used for qname re-discovery.
    uint64_t vb_size;                    // ZIP/PIZ: compression VBlock size in bytes (PIZ: passed in SectionHeaderGenozipHeader.vb_size)
    double gz_comp_ratio;                // GZ compression ratio of segconf data
    uint32_t gz_comp_size;               // size of segconf data in source gz compression, used if discover_during_segconf 
    uint32_t line_len;                   // approx line len
    bool zip_txt_modified;               // ZIP/PIZ: txt data is/was modified during Seg (e.g. by --optimize, --add-line-numbers). Before segconf: true if data *might* be modifed. After segconf: true iff data is modified.

    // SAM/BAM and FASTQ: not used during normal segging or recon of every line (every VB)
    #define NUM_QUAL_SCORES 94
    QualHisto qual_histo[NUM_QHT][NUM_QUAL_SCORES]; // histogram of qual values in segconf data for each QHT. Also accessed when in initiazing HOMP and reported in Stats

    #define segconf_optimize(did_i) bitset_get (segconf.needs_optimize, (did_i))
    const uint64_t needs_optimize[(MAX_DICTS+63)/64]; // true if --optimize indicates that this field should be optimized 
    double smux_max_stdv;                // for a particular quality score q - stdv across the 5 values (corresponding to b=A,C,T,G,other) which are the % of qual data which has score q, looking only at the quality values of associated with b. smux_max_stdv is maximal of the 94 stdv's asscoiated with all q's.
    #define SAM_FACTOR_MULT 32           // multiplication for storing est_sam_factor as an int in SectionHeaderGenozipHeader
    double est_sam_factor;               // BAM only: calculated in Segconf and used in PIZ: est size factor when translating BAM to SAM based on segconf data
    QnameFlavor deep_sam_qname_flavor[2];// save for stats, in --deep (SAM's QNAME and QNAME2)

    char smux_max_stdv_q;                // the quality score for which results in smux_max_stdv 
    uint8_t sam_cigar_len;               // Seg: approx average CIGAR len rounded up (max 255) - used for initial allocation

    // FASTQ (every VB)
    uint32_t optimized_qname_len;
    StrText optimized_qname;             // --optimize: prefix of optimized qname
    char interleaved_r1, interleaved_r2; // valid if is_interleaved: character representing R1/R2: usually '0' or '1' for R1 and '1' or '2' for R2
    bool multiseq;                       // FASTQ/FASTA: sequences in file are variants of each others 

    // VCF
    Mutex PL_mux_by_DP_mutex;

    // ________________________________________________________________________________________
    // Used ONLY during segconf and possibly right at the end - 
    // stats, writing global sections, error messages...
    // These fields can get evicted out of all cache levels during the bulk of the zip/piz time
    // ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾ 
    // General (only beginning or end of execution)
    bool qname_barcode2_is_alias[NUM_QTYPES];   // qname has two barcodes which are aliases

    // SAM/BAM (only beginning or end of execution)
    char sam_deep_filename[1024];        // name of the SAM filename in --deep 
    QnameStr sam_qname_line0;            // copy of qname_line0[QNAME1], but survives to fastq in deep/bammass (nul-terminated)
    uint32_t AS_is_2ref_consumed_evidence; // Segconf: AS value tends to be double ref_consumed (counter during segconf, bool during seg/piz)
    uint32_t AS_is_ref_consumed_evidence;  // Segconf: AS value tends to be double ref_consumed (counter during segconf, bool during seg/piz)
    uint32_t num_mapped;                 // Segconf: number of segconf reads that are mapped - defined as having (!FLAG.unmapped, RNAME, POS, CIGAR)
    uint32_t cummul_sam_cigars_len;      // Segconf: cumulative length of all cigars
    uint32_t cummul_seq_len_to_cm;       // ZIP/PIZ: approx average of (seq_len/cm:i) (during segconf.running - cumulative)
    uint32_t cummul_s1_to_cm_32;         // ZIP/PIZ: approx average of (s1:i/cm:i) x 32 (during segconf.running - cumulative)
    uint32_t commul_translated_sam_size; // BAM only: est size of segconf vblock translated to SAM

    #define MIN_XCONS_STD_QUAL_LEN 100
    #define MAX_XCONS_STD_QUAL_LEN (MIN_XCONS_STD_QUAL_LEN + 255) // SectionHeaderGenozipHeader.xcons_std_seq_len_M100 is uint8_t
    uint32_t xcons_std_line_histogram[MAX_XCONS_STD_QUAL_LEN - MIN_XCONS_STD_QUAL_LEN + 1]; // Segconf: used to find the most common standard xcons qual length ("standard": alignments with XC:i but without XO:i)

    uint8_t MAPQ_value;                  // Segconf: calculate sam_mapq_has_single_value
    SeqTech tech_by_RG;                  // tech by first @RG PL field
    bool evidence_of_collated;           // Segconf: at least a one pair of consecutive lines has the same QNAME
    bool evidence_of_sorted;             // Segconf: at least a one pair of consecutive lines has the same RNAME and increasing POS

    // FASTQ (only beginning or end of execution)
    #define NUM_INSTS 6
    unsigned n_full_mch[NUM_INSTS];      // Deep: count segconf lines where hash matches with at least one SAM line - (QNAME1 or QNAME2), SEQ, QUAL
    unsigned n_full_mch_trimmed;         // Deep: the size of the subset of n_full_mch which are trimmed     
    unsigned n_seq_qname_mch[NUM_INSTS]; // Deep: count segconf lines where FASTQ (QNAME1 or QNAME2) hash matches with at least one SAM line - SEQ and QNAME 
    unsigned n_no_mch;                   // Deep: count segconf lines that don't match any SAM line (perhaps because SAM is filtered)
    uint32_t total_usable_len;           // used by bamass_segconf to calculate line_len
    char deep_1st_desc[256];             // Deep: DESC line of first FASTQ read of first FASTQ file

    BamAssTrimCigarTreatment bamass_trims; // bamass: how to treat a trim in the CIGAR
    PairType r1_or_r2;                   // in case compression is WITHOUT --pair: our guess of whether this file is R1 or R2

    // VCF (only beginning or end of execution)
    uint32_t count_GQ_by_PL, count_GQ_by_GP; // Segconf: used tp calculate GQ_by_PL, GQ_by_GP
    uint32_t n_Q_to_O;                   // Segconf: number of values added up in Q_to_O
    bool vcf_evidence_not_gvcf;          // Segconf

    // GFF
    bool has_embedded_fasta;             // GFF: embedded FASTA encountered - from when set until the end the file 

} SegConf;

extern SegConf segconf; // ZIP: set based on segging a sample of a few first lines of the file
                        // PIZ: select fields are transferred through SectionHeaderGenozipHeader

static inline void segconf_set_has (Did did_i)
{
    bitset_set ((uint64_t *)segconf.has_bits, did_i); // only place in the code where we override the "const" of has_bits
}

// PIZ: copy segconf value sent from ZIP via GenozipHeader
static inline void segconf_load_has (Did did_i, bool value)
{
    bitset_cpy ((uint64_t *)segconf.has_bits, did_i, value); // only place in the code where we override the "const" of has_bits
}

static inline void segconf_set_optimize (Did did_i, bool value)
{
    bitset_cpy ((uint64_t *)segconf.needs_optimize, did_i, value); // only place in the code where we override the "const" of has_bits
}

extern void segconf_zip_initialize (void);
extern void segconf_free (void);
extern void segconf_calculate (void);
extern void segconf_set_vb_size (VBlockP vb, uint64_t curr_vb_size);
extern bool segconf_is_long_reads(void);
extern void segconf_set_use_insertion_ctxs (void);
extern void segconf_mark_as_used (VBlockP vb, unsigned num_ctxs, ...);
extern rom segconf_sam_mapper_name (void);
extern rom tech_name (SeqTech tech);
extern rom segconf_deep_trimming_name (void);
extern rom VCF_QUAL_method_name (VcfQualMethod method);
extern rom VCF_INFO_method_name (VcfInfoMethod method);
extern rom FMT_GQ_method_name (GQMethodType method);
extern rom FMT_ROAO_method_name (ROAOMethodType method);
extern rom FMT_DP_method_name (FormatDPMethod method);
extern rom INFO_DP_method_name (InfoDPMethod method);
extern rom FMT_GP_content_name (GPContentType method);
extern rom RG_method_name (RGMethod method);
extern void segconf_test_sorted (VBlockP vb, WordIndex prev_line_chrom, PosType32 pos, PosType32 prev_line_pos);
extern void segconf_test_multiseq (VBlockP vb, Did nonref);
extern StrText segconf_get_qual_histo (QualHistType qht);
extern unsigned segconf_get_num_qual_scores (QualHistType qht);
extern StrText1K segconf_get_optimizations (void);
extern void segconf_piz_initialize (void);

#define segconf_running __builtin_expect (segconf.running, false)

static bool inline tech_is_unknown (void) { return TECH(UNKNOWN) || TECH(CONS) || TECH(NCBI); }
