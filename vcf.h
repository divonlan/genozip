// ------------------------------------------------------------------
//   vcf.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef VCF_INCLUDED
#define VCF_INCLUDED

#include "genozip.h"
#include "digest.h"
#include "sections.h"

#define DTYPE_VCF_INFO   DTYPE_1
#define DTYPE_VCF_FORMAT DTYPE_2

#pragma GENDICT_PREFIX VCF

// Fields
#pragma GENDICT VCF_CHROM=DTYPE_FIELD=CHROM // CHROM must be first
#pragma GENDICT VCF_POS=DTYPE_FIELD=POS    
#pragma GENDICT VCF_ID=DTYPE_FIELD=ID
#pragma GENDICT VCF_REFALT=DTYPE_FIELD=REF+ALT
#pragma GENDICT VCF_QUAL=DTYPE_FIELD=QUAL
#pragma GENDICT VCF_FILTER=DTYPE_FIELD=FILTER
#pragma GENDICT VCF_INFO=DTYPE_FIELD=INFO
#pragma GENDICT VCF_FORMAT=DTYPE_FIELD=FORMAT
#pragma GENDICT VCF_SAMPLES=DTYPE_FIELD=SAMPLES
#pragma GENDICT VCF_EOL=DTYPE_FIELD=EOL
#pragma GENDICT VCF_TOPLEVEL=DTYPE_FIELD=TOPLEVEL // must be called TOPLEVEL
#pragma GENDICT VCF_oCHROM=DTYPE_FIELD=oCHROM     
#pragma GENDICT VCF_oPOS=DTYPE_FIELD=oPOS
#pragma GENDICT VCF_oREFALT=DTYPE_FIELD=oREFALT
#pragma GENDICT VCF_oXSTRAND=DTYPE_FIELD=oXSTRAND
#pragma GENDICT VCF_COORDS=DTYPE_FIELD=COORDS
#pragma GENDICT VCF_oSTATUS=DTYPE_FIELD=o$TATUS
#pragma GENDICT VCF_COPYPOS=DTYPE_FIELD=C0PYPOS
#pragma GENDICT VCF_LIFT_REF=DTYPE_FIELD=LIFT_REF
#pragma GENDICT VCF_COPYSTAT=DTYPE_FIELD=CoPYSTAT
#pragma GENDICT VCF_TOPLUFT=DTYPE_FIELD=ToPLUFT
#pragma GENDICT VCF_LINE_NUM=DTYPE_FIELD=LINE_NUM

// FORMAT fields
#pragma GENDICT FORMAT_AD=DTYPE_2=AD       // <ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
#pragma GENDICT FORMAT_ADF=DTYPE_2=ADF     // <ID=ADF,Number=R,Type=Float,Description="Allele dosage on fwd strand">
#pragma GENDICT FORMAT_ADR=DTYPE_2=ADR     // <ID=ADR,Number=R,Type=Float,Description="Allele dosage on rev strand">
#pragma GENDICT FORMAT_ADALL=DTYPE_2=ADALL // from GIAB: <ID=ADALL,Number=R,Type=Integer,Description="Net allele depths across all datasets">
#pragma GENDICT FORMAT_AF=DTYPE_2=AF       // <ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed">
#pragma GENDICT FORMAT_DP=DTYPE_2=DP       // <ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
#pragma GENDICT FORMAT_DS=DTYPE_2=DS
#pragma GENDICT FORMAT_GL=DTYPE_2=GL
#pragma GENDICT FORMAT_GP=DTYPE_2=GP       // <ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification">
#pragma GENDICT FORMAT_GQ=DTYPE_2=GQ       // <ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#pragma GENDICT FORMAT_GT=DTYPE_2=GT       // <ID=GT,Number=1,Type=String,Description="Genotype">
#pragma GENDICT FORMAT_PL=DTYPE_2=PL       // <ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
#pragma GENDICT FORMAT_PRI=DTYPE_2=PRI     // <ID=PRI,Number=G,Type=Float,Description="Phred-scaled prior probabilities for genotypes">
#pragma GENDICT FORMAT_F1R2=DTYPE_2=F1R2   // <ID=F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting each allele"> see: https://github.com/broadinstitute/gatk/blob/master/docs/mutect/mutect.pdf
#pragma GENDICT FORMAT_F2R1=DTYPE_2=F2R1   // <ID=F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting each allele">
#pragma GENDICT FORMAT_MB=DTYPE_2=MB       // <ID=MB,Number=4,Type=Integer,Description="Per-sample component statistics to detect mate bias">
#pragma GENDICT FORMAT_PP=DTYPE_2=PP       // <ID=PP,Number=G,Type=Integer,Description="Phred-scaled genotype posterior probabilities rounded to the closest integer">
#pragma GENDICT FORMAT_SAC=DTYPE_2=SAC     // <ID=SAC,Number=.,Type=Integer,Description="Number of reads on the forward and reverse strand supporting each allele (including reference)">
#pragma GENDICT FORMAT_SB=DTYPE_2=SB       // <ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias">
#pragma GENDICT FORMAT_PS=DTYPE_2=PS 

// PBWT fields - same dict_id for all data types using PBWT, as codec_pbwt_uncompress relies on it
#pragma GENDICT FORMAT_GT_HT=DTYPE_2=@HT    
#pragma GENDICT FORMAT_GT_HT_INDEX=DTYPE_2=@INDEXHT // different first 2 letters
#pragma GENDICT FORMAT_PBWT_RUNS=DTYPE_2=@1BWTRUN // PBWT runs - MUST have a did_i higher that FORMAT_GT_HT's
#pragma GENDICT FORMAT_PBWT_FGRC=DTYPE_2=@2BWTFGR // PBWT foreground run count - MUST be right after FORMAT_PBWT_RUNS

#ifndef _PBWT_RUNS
#define _PBWT_RUNS        _FORMAT_PBWT_RUNS       
#define _PBWT_FGRC        _FORMAT_PBWT_FGRC       
#define _PBWT_HT_MATRIX   _FORMAT_GT_HT           
#define _PBWT_GT_HT_INDEX _FORMAT_GT_HT_INDEX     
#endif

// INFO fields
#pragma GENDICT INFO_AC=DTYPE_1=AC         // <ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_AF=DTYPE_1=AF         // <ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_AN=DTYPE_1=AN         // <ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#pragma GENDICT INFO_AA=DTYPE_1=AA         // <ID=AA,Number=1,Type=String,Description="Ancestral Allele"> - defined in the VCF specification
#pragma GENDICT INFO_BaseCounts=DTYPE_1=BaseCounts
#pragma GENDICT INFO_DP=DTYPE_1=DP         // <ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
#pragma GENDICT INFO_DP4=DTYPE_1=DP4   
#pragma GENDICT INFO_SF=DTYPE_1=SF
#pragma GENDICT INFO_VQSLOD=DTYPE_1=VQSLOD // <ID=VQSLOD,Number=1,Type=Float,Description="Log odds of being a true variant versus being false under the trained Gaussian mixture model">

// Ann field defined: https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
#pragma GENDICT INFO_ANN=DTYPE_1=ANN       // <ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
#pragma GENDICT INFO_ANN_Allele=DTYPE_1=@ANN_Allele // note: nice looking tag, as might be displayed in INFO/Lref
#pragma GENDICT INFO_ANN_Annotation=DTYPE_1=A1Annotation
#pragma GENDICT INFO_ANN_Annotation_Impact=DTYPE_1=A2Annotation_Impact
#pragma GENDICT INFO_ANN_Gene_Name=DTYPE_1=A3Gene_Name
#pragma GENDICT INFO_ANN_Gene_ID=DTYPE_1=A4Gene_ID
#pragma GENDICT INFO_ANN_Feature_Type=DTYPE_1=A5Feature_Type
#pragma GENDICT INFO_ANN_Feature_ID=DTYPE_1=A6Feature_ID
#pragma GENDICT INFO_ANN_Transcript_BioType=DTYPE_1=A7Transcript_BioType
#pragma GENDICT INFO_ANN_Rank=DTYPE_1=A8Rank
#pragma GENDICT INFO_ANN_HGVS_c=DTYPE_1=A9HGVS_c
#pragma GENDICT INFO_ANN_HGVS_p=DTYPE_1=AaHGVS_p
#pragma GENDICT INFO_ANN_cDNA=DTYPE_1=AbcDNA
#pragma GENDICT INFO_ANN_CDS=DTYPE_1=AcCDS
#pragma GENDICT INFO_ANN_AA=DTYPE_1=AdAA
#pragma GENDICT INFO_ANN_Distance=DTYPE_1=AeDistance
#pragma GENDICT INFO_ANN_Errors=DTYPE_1=AfErrors

// Added by GATK HaplotypeCaller in a gVCF: https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
#pragma GENDICT INFO_END=DTYPE_1=END           // <ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
#pragma GENDICT INFO_MLEAC=DTYPE_1=MLEAC       // <ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_MLEAF=DTYPE_1=MLEAF       // <ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_LDAF=DTYPE_1=LDAF         //  MLE Allele Frequency Accounting for LD
#pragma GENDICT INFO_MQ0=DTYPE_1=MQ0    
#pragma GENDICT FORMAT_MIN_DP=DTYPE_1=MIN_DP   // <ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">

// Ensembl VEP (Variant Effect Predictor) fields: https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html
#pragma GENDICT INFO_vep=DTYPE_1=vep
#pragma GENDICT INFO_AGE_HISTOGRAM_HET=DTYPE_1=AGE_HISTOGRAM_HET 
#pragma GENDICT INFO_AGE_HISTOGRAM_HOM=DTYPE_1=AGE_HISTOGRAM_HOM
#pragma GENDICT INFO_MAX_AF=DTYPE_1=MAX_AF     // highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD

#pragma GENDICT INFO_CSQ=DTYPE_1=CSQ           // <ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral">
#pragma GENDICT INFO_CSQ_Allele=DTYPE_1=@CSQ_Allele              // 0 (note: nice looking tag, as might be displayed in INFO/Lref)
#pragma GENDICT INFO_CSQ_Consequence=DTYPE_1=c1SQ_Consequence
#pragma GENDICT INFO_CSQ_IMPACT=DTYPE_1=c2SQ_IMPACT
#pragma GENDICT INFO_CSQ_SYMBOL=DTYPE_1=c3SQ_SYMBOL
#pragma GENDICT INFO_CSQ_Gene=DTYPE_1=c4SQ_Gene
#pragma GENDICT INFO_CSQ_Feature=DTYPE_1=c6SQ_Feature            // 5-7 together as Feature implies Feature_Type and Transcript_BioType
#pragma GENDICT INFO_CSQ_EXON=DTYPE_1=c8SQ_EXON
#pragma GENDICT INFO_CSQ_INTRON=DTYPE_1=c9SQ_INTRON
#pragma GENDICT INFO_CSQ_HGVSc=DTYPE_1=caSQ_HGVSc                // 10
#pragma GENDICT INFO_CSQ_HGVSp=DTYPE_1=cbSQ_HGVSp
#pragma GENDICT INFO_CSQ_cDNA_position=DTYPE_1=ccSQ_cDNA_position
#pragma GENDICT INFO_CSQ_CDS_position=DTYPE_1=cdSQ_CDS_position
#pragma GENDICT INFO_CSQ_Protein_position=DTYPE_1=ceSQ_Protein_position
#pragma GENDICT INFO_CSQ_Amino_acids=DTYPE_1=cfSQ_Amino_acids
#pragma GENDICT INFO_CSQ_Codons=DTYPE_1=cgSQ_Codons
#pragma GENDICT INFO_CSQ_Existing_variation=DTYPE_1=chSQ_Existing_variation
#pragma GENDICT INFO_CSQ_ALLELE_NUM=DTYPE_1=ciSQ_ALLELE_NUM
#pragma GENDICT INFO_CSQ_DISTANCE=DTYPE_1=cjSQ_DISTANCE
#pragma GENDICT INFO_CSQ_STRAND=DTYPE_1=ckSQ_STRAND              // 20
#pragma GENDICT INFO_CSQ_FLAGS=DTYPE_1=clSQ_FLAGS
#pragma GENDICT INFO_CSQ_VARIANT_CLASS=DTYPE_1=cmSQ_VARIANT_CLASS
#pragma GENDICT INFO_CSQ_MINIMISED=DTYPE_1=cnSQ_MINIMISED
#pragma GENDICT INFO_CSQ_SYMBOL_SOURCE=DTYPE_1=coSQ_SYMBOL_SOURCE
#pragma GENDICT INFO_CSQ_HGNC_ID=DTYPE_1=cpSQ_HGNC_ID
#pragma GENDICT INFO_CSQ_CANONICAL=DTYPE_1=cqSQ_CANONICAL
#pragma GENDICT INFO_CSQ_TSL=DTYPE_1=crSQ_TSL
#pragma GENDICT INFO_CSQ_APPRIS=DTYPE_1=csSQ_APPRIS
#pragma GENDICT INFO_CSQ_CCDS=DTYPE_1=ctSQ_CCDS
#pragma GENDICT INFO_CSQ_ENSP=DTYPE_1=cuSQ_ENSP                  // 30
#pragma GENDICT INFO_CSQ_SWISSPROT=DTYPE_1=cvSQ_SWISSPROT
#pragma GENDICT INFO_CSQ_TREMBL=DTYPE_1=cwSQ_TREMBL
#pragma GENDICT INFO_CSQ_UNIPARC=DTYPE_1=czSQ_UNIPARC
#pragma GENDICT INFO_CSQ_GENE_PHENO=DTYPE_1=cySQ_GENE_PHENO
#pragma GENDICT INFO_CSQ_SIFT=DTYPE_1=czSQ_SIFT
#pragma GENDICT INFO_CSQ_PolyPhen=DTYPE_1=cASQ_PolyPhen
#pragma GENDICT INFO_CSQ_DOMAINS=DTYPE_1=cBSQ_DOMAINS
#pragma GENDICT INFO_CSQ_HGVS_OFFSET=DTYPE_1=cCSQ_HGVS_OFFSET
#pragma GENDICT INFO_CSQ_AF=DTYPE_1=cDSQ_AF                      // 39-55 all 17 AF fields together, as their sequence in b250 is repetative
#pragma GENDICT INFO_CSQ_CLIN_SIG=DTYPE_1=cUSQ_CLIN_SIG
#pragma GENDICT INFO_CSQ_SOMATIC=DTYPE_1=cVSQ_SOMATIC
#pragma GENDICT INFO_CSQ_PHENO=DTYPE_1=cQSQ_PHENO
#pragma GENDICT INFO_CSQ_PUBMED=DTYPE_1=cXSQ_PUBMED
#pragma GENDICT INFO_CSQ_MOTIF_NAME=DTYPE_1=cYSQ_MOTIF_NAME      // 60
#pragma GENDICT INFO_CSQ_MOTIF_POS=DTYPE_1=cZSQ_MOTIF_POS
#pragma GENDICT INFO_CSQ_HIGH_INF_POS=DTYPE_1=c@SQ_HIGH_INF_POS
#pragma GENDICT INFO_CSQ_MOTIF_SCORE_CHANGE=DTYPE_1=c$SQ_MOTIF_SCORE_CHANGE
#pragma GENDICT INFO_CSQ_LoF=DTYPE_1=c%SQ_LoF
#pragma GENDICT INFO_CSQ_LoF_filter=DTYPE_1=c^SQ_LoF_filter
#pragma GENDICT INFO_CSQ_LoF_flags=DTYPE_1=c-SQ_LoF_flags
#pragma GENDICT INFO_CSQ_LoF_info=DTYPE_1=c*SQ_LoF_info
#pragma GENDICT INFO_CSQ_context=DTYPE_1=c_SQ_context
#pragma GENDICT INFO_CSQ_ancestral=DTYPE_1=c+SQ_ancestral

// clinvar
#pragma GENDICT INFO_ALLELEID=DTYPE_1=ALLELEID // <ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
#pragma GENDICT INFO_CLNDN=DTYPE_1=CLNDN       // <ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
#pragma GENDICT INFO_RS=DTYPE_1=RS             // <ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
#pragma GENDICT INFO_CLNHGVS=DTYPE_1=CLNHGVS   // <ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
#pragma GENDICT INFO_CLNVI=DTYPE_1=CLNVI       // <ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">

// HGVS - used for CLNHGVS and and also INFO/ANN/HGVS
#pragma GENDICT INFO_HGVS_snp_pos=DTYPE_1=H0GVS_snp_pos   
#pragma GENDICT INFO_HGVS_snp_refalt=DTYPE_1=H1GVS_snp_refalt
#pragma GENDICT INFO_HGVS_del_start_pos=DTYPE_1=H2GVS_startpos_del
#pragma GENDICT INFO_HGVS_del_end_pos=DTYPE_1=H3GVS_endpos_del
#pragma GENDICT INFO_HGVS_del_payload=DTYPE_1=H4GVS_payload_del
#pragma GENDICT INFO_HGVS_ins_start_pos=DTYPE_1=H5GVS_startpos_ins // used for ins and delins
#pragma GENDICT INFO_HGVS_ins_end_pos=DTYPE_1=H6GVS_endpos_ins     
#pragma GENDICT INFO_HGVS_ins_payload=DTYPE_1=H7GVS_payload_ins
#pragma GENDICT INFO_HGVS_delins_end_pos=DTYPE_1=H8GVS_endpos_delins     
#pragma GENDICT INFO_HGVS_delins_payload=DTYPE_1=H9GVS_payload_delins

// ExAC fields
#pragma GENDICT INFO_DP_HIST=DTYPE_1=DP_HIST  // from ExAC: Depth (DP) histogram in 20 equal intervals between 0-100 : See https://www.biorxiv.org/content/biorxiv/suppl/2015/10/30/030338.DC1/030338-1.pdf
#pragma GENDICT INFO_GQ_HIST=DTYPE_1=GQ_HIST  // from ExAC: Genotype Quality (GQ) histogram in 20 equal intervals between 0-100

// Structural variants (also uses INFO/END): https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/VCF%20(Variant%20Call%20Format)%20version%204.0/encoding-structural-variants/
#pragma GENDICT INFO_SVLEN=DTYPE_1=SVLEN

// genozip INFO fields
#pragma GENDICT INFO_LUFT=DTYPE_1=LUFT
#pragma GENDICT INFO_PRIM=DTYPE_1=PRIM
#pragma GENDICT INFO_LREJ=DTYPE_1=Lrej
#pragma GENDICT INFO_PREJ=DTYPE_1=Prej

#define VCF_MAX_PLOIDY 100  // set to a reasonable 100 to avoid memory allocation explosion in case of an error in the VCF file
#if VCF_MAX_PLOIDY > 65535
#error "VCF_MAX_PLOIDY cannot go beyond 65535 VBlockVCF.ploidy are uint16_t"
#endif

#define MAX_ALLELES 100 // REF (allele #0) + 99 ALTs (alleles # 1-99)
typedef uint8_t Allele; // elements of ht_matrix: values 48->147 for allele 0 to 99, '*' for unused, '%', '-'

// ZIP stuff
extern void vcf_zip_initialize (void);
extern void vcf_zip_read_one_vb (VBlockP vb);
extern void vcf_liftover_display_lift_report (void);

// SEG stuff
extern const char *vcf_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void vcf_seg_initialize (VBlockP vb_);
extern void vcf_zip_after_compute (VBlockP vb);
extern void vcf_seg_finalize (VBlockP vb_);
extern bool vcf_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern TranslatorId vcf_lo_luft_trans_id (DictId dict_id, char number);

// PIZ stuff
extern bool vcf_piz_read_one_vb (VBlockP vb, Section sl);
extern bool vcf_vb_is_luft (VBlockP vb);
extern bool vcf_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
CONTAINER_FILTER_FUNC (vcf_piz_filter);
CONTAINER_CALLBACK (vcf_piz_container_cb);

// VCF Header stuff
extern void vcf_header_piz_init (void);
extern bool vcf_inspect_txt_header (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);
extern uint32_t vcf_header_get_num_samples (void);
extern bool vcf_header_get_has_fileformat (void);

// VBlock stuff
extern void vcf_vb_release_vb();
extern void vcf_vb_destroy_vb();
extern void vcf_vb_cleanup_memory();
extern unsigned vcf_vb_size (DataType dt);
extern unsigned vcf_vb_zip_dl_size (void);
extern bool vcf_vb_has_haplotype_data (VBlockP vb);

// Liftover - INFO fields
#define INFO_LUFT_NAME  "LUFT"
#define INFO_PRIM_NAME  "PRIM"
#define INFO_LREJ_NAME  "Lrej"
#define INFO_PREJ_NAME  "Prej" // lower case so Prej doesn't have the same first 2 chars as PRIM (to not conflict in dict_id_to_did_i_map)
#define INFO_DVCF_LEN (sizeof INFO_LUFT_NAME - 1) // these 4 must be the same length

#define VCF_CONTIG_FMT "##contig=<ID=%.*s,length=%"PRId64">"

// DVCF tag renaming stuff
extern void vcf_tags_cmdline_drop_option(void);
extern void vcf_tags_cmdline_rename_option(void);

// Samples stuff
extern void vcf_samples_add  (const char *samples_str);

#define VCF_SPECIAL { vcf_piz_special_main_REFALT, vcf_piz_special_FORMAT, vcf_piz_special_INFO_AC, vcf_piz_special_INFO_SVLEN, \
                      vcf_piz_special_FORMAT_DS, vcf_piz_special_INFO_BaseCounts, vcf_piz_special_INFO_SF, vcf_piz_special_MINUS,  \
                      vcf_piz_special_LIFT_REF, vcf_piz_special_COPYSTAT, vcf_piz_special_other_REFALT, vcf_piz_special_COPYPOS, vcf_piz_special_ALLELE, \
                      vcf_piz_special_INFO_HGVS_SNP_POS, vcf_piz_special_INFO_HGVS_SNP_REFALT, \
                      vcf_piz_special_INFO_HGVS_DEL_END_POS, vcf_piz_special_INFO_HGVS_DEL_PAYLOAD, \
                      vcf_piz_special_INFO_HGVS_INS_END_POS, vcf_piz_special_INFO_HGVS_INS_PAYLOAD, \
                      vcf_piz_special_INFO_HGVS_DELINS_END_POS, vcf_piz_special_INFO_HGVS_DELINS_PAYLOAD }

SPECIAL (VCF, 0,  main_REFALT,         vcf_piz_special_main_REFALT);
SPECIAL (VCF, 1,  FORMAT,              vcf_piz_special_FORMAT)
SPECIAL (VCF, 2,  AC,                  vcf_piz_special_INFO_AC);
SPECIAL (VCF, 3,  SVLEN,               vcf_piz_special_INFO_SVLEN);
SPECIAL (VCF, 4,  DS,                  vcf_piz_special_FORMAT_DS);
SPECIAL (VCF, 5,  BaseCounts,          vcf_piz_special_INFO_BaseCounts);
SPECIAL (VCF, 6,  SF,                  vcf_piz_special_INFO_SF);
SPECIAL (VCF, 7,  MINUS,               vcf_piz_special_MINUS);                    // added v12.0.0
SPECIAL (VCF, 8,  LIFT_REF,            vcf_piz_special_LIFT_REF);                 // added v12.0.0
SPECIAL (VCF, 9,  COPYSTAT,            vcf_piz_special_COPYSTAT);                 // added v12.0.0
SPECIAL (VCF, 10, other_REFALT,        vcf_piz_special_other_REFALT);             // added v12.0.0
SPECIAL (VCF, 11, COPYPOS,             vcf_piz_special_COPYPOS);                  // added v12.0.0
SPECIAL (VCF, 12, ALLELE,              vcf_piz_special_ALLELE);                   // added v12.0.0
SPECIAL (VCF, 13, HGVS_SNP_POS,        vcf_piz_special_INFO_HGVS_SNP_POS);        // added v12.0.15
SPECIAL (VCF, 14, HGVS_SNP_REFALT,     vcf_piz_special_INFO_HGVS_SNP_REFALT);     // added v12.0.15
SPECIAL (VCF, 15, HGVS_DEL_END_POS,    vcf_piz_special_INFO_HGVS_DEL_END_POS);    // added v12.0.34
SPECIAL (VCF, 16, HGVS_DEL_PAYLOAD,    vcf_piz_special_INFO_HGVS_DEL_PAYLOAD);    // added v12.0.34
SPECIAL (VCF, 17, HGVS_INS_END_POS,    vcf_piz_special_INFO_HGVS_INS_END_POS);    // added v12.0.34
SPECIAL (VCF, 18, HGVS_INS_PAYLOAD,    vcf_piz_special_INFO_HGVS_INS_PAYLOAD);    // added v12.0.34
SPECIAL (VCF, 19, HGVS_DELINS_END_POS, vcf_piz_special_INFO_HGVS_DELINS_END_POS); // added v12.0.34
SPECIAL (VCF, 20, HGVS_DELINS_PAYLOAD, vcf_piz_special_INFO_HGVS_DELINS_PAYLOAD); // added v12.0.34
#define NUM_VCF_SPECIAL 21

// Translators for Luft (=secondary coordinates)
TRANSLATOR (VCF, VCF,   1,  G,      vcf_piz_luft_G)       // same order as LiftOverStatus starting LO_CANT_G
TRANSLATOR (VCF, VCF,   2,  R,      vcf_piz_luft_R)
TRANSLATOR (VCF, VCF,   3,  R2,     vcf_piz_luft_R2)
TRANSLATOR (VCF, VCF,   4,  A_AN,   vcf_piz_luft_A_AN)
TRANSLATOR (VCF, VCF,   5,  A_1,    vcf_piz_luft_A_1)
TRANSLATOR (VCF, VCF,   6,  PLOIDY, vcf_piz_luft_PLOIDY)
TRANSLATOR (VCF, VCF,   7,  GT,     vcf_piz_luft_GT)      
TRANSLATOR (VCF, VCF,   8,  END,    vcf_piz_luft_END)      
TRANSLATOR (VCF, VCF,   9,  XREV,   vcf_piz_luft_XREV)      
TRANSLATOR (VCF, VCF,   10, ALLELE, vcf_piz_luft_ALLELE)      

#define NUM_VCF_TRANS   11 // including "none"
#define VCF_TRANSLATORS { NULL /* none */, vcf_piz_luft_G, vcf_piz_luft_R, vcf_piz_luft_R2, vcf_piz_luft_A_AN, \
                          vcf_piz_luft_A_1, vcf_piz_luft_PLOIDY, vcf_piz_luft_GT, vcf_piz_luft_END, vcf_piz_luft_XREV, vcf_piz_luft_ALLELE }

typedef struct {
    const char *alg_name;
    enum { TW_NEVER, TW_ALWAYS, TW_REF_ALT_SWITCH, TW_XSTRAND } upon;
} LuftTransLateProp;

// names of INFO / FORMAT algorithms, goes into VCF header's ##INFO / ##FORMAT "RendAlg" attribute
                           /* Algorithm   Trigger          */
#define DVCF_TRANS_PROPS { { "NONE",      TW_NEVER          },   /* never translate */\
                           { "G",         TW_REF_ALT_SWITCH },   /* reshuffle a  'G' vector (one element per genotype) if REF<>ALT changed */\
                           { "R",         TW_REF_ALT_SWITCH },   /* reshuffle an 'R' vector (one element per ref/alt allele) if REF<>ALT changed */\
                           { "R2",        TW_REF_ALT_SWITCH },   /* reshuffle a vector with 2 elements per ref/alt allele, if REF<>ALT changed */\
                           { "A_AN",      TW_REF_ALT_SWITCH },   /* recalculate an 'A' vector (one element per ALT allele) if REF<>ALT changed, who's elements, including a missing element for REF, add up to AN (example: AC). */ \
                           { "A_1",       TW_REF_ALT_SWITCH },   /* recalculate an 'A' vector (one element per ALT allele) if REF<>ALT changed, who's elements, including a missing element for REF, add up to 1 (example: AF). */ \
                           { "PLOIDY",    TW_REF_ALT_SWITCH },   /* recalculate a float to (ploidy-value) */ \
                           { "GT",        TW_REF_ALT_SWITCH },   /* recalculate the allele numbers FORMAT/GT if REF<>ALT changed */ \
                           { "END",       TW_ALWAYS         },   /* recalculate INFO/END */\
                           { "XREV",      TW_XSTRAND        },   /* reverse the elements of a vector if XSTRAND. Example: INFO/BaseCounts */\
                           { "ALLELE",    TW_ALWAYS         } }  /* copy an allele verbatim including if changes or changes order. example: INFO/AA */ 
                           
extern const LuftTransLateProp ltrans_props[NUM_VCF_TRANS];

#define needs_translation(ctx)  (z_dual_coords && (ctx)->luft_trans && \
    ((ltrans_props[(ctx)->luft_trans].upon == TW_REF_ALT_SWITCH && LO_IS_OK_SWITCH (last_ostatus)) || \
     (ltrans_props[(ctx)->luft_trans].upon == TW_ALWAYS         && LO_IS_OK (last_ostatus))        || \
     (ltrans_props[(ctx)->luft_trans].upon == TW_XSTRAND        && LO_IS_OK (last_ostatus) && *CTX(VCF_oXSTRAND)->last_snip != '-')))

#define VCF_DICT_ID_ALIASES \
    /*         alias             maps to this ctx  */  \
    { DT_VCF,  _INFO_END, _VCF_POS    }, \

#define VCF_LOCAL_GET_LINE_CALLBACKS

#define dict_id_is_vcf_info_sf   dict_id_is_type_1
#define dict_id_is_vcf_format_sf dict_id_is_type_2

#endif
