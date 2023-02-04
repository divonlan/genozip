// ------------------------------------------------------------------
//   vcf.h
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "digest.h"
#include "sections.h"

#define DTYPE_VCF_INFO   DTYPE_1
#define DTYPE_VCF_FORMAT DTYPE_2

#pragma GENDICT_PREFIX VCF

// Fields
#pragma GENDICT VCF_CHROM=DTYPE_FIELD=CHROM         // CHROM must be first
#pragma GENDICT VCF_POS=DTYPE_FIELD=POS    
#pragma GENDICT VCF_ID=DTYPE_FIELD=ID
#pragma GENDICT VCF_REFALT=DTYPE_FIELD=REF+ALT
#pragma GENDICT VCF_QUAL=DTYPE_FIELD=QUAL
#pragma GENDICT VCF_FILTER=DTYPE_FIELD=FILTER
#pragma GENDICT VCF_INFO=DTYPE_FIELD=INFO
#pragma GENDICT VCF_FORMAT=DTYPE_FIELD=FORMAT
#pragma GENDICT VCF_SAMPLES=DTYPE_FIELD=SAMPLES
#pragma GENDICT VCF_LOOKBACK=DTYPE_FIELD=LOOKBACK   // samples lookback
#pragma GENDICT VCF_EOL=DTYPE_FIELD=EOL
#pragma GENDICT VCF_TOPLEVEL=DTYPE_FIELD=TOPLEVEL   // must be called TOPLEVEL
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
#pragma GENDICT VCF_DEBUG_LINES=DTYPE_FIELD=DBGLINES// used by --debug-lines

// FORMAT fields
#pragma GENDICT FORMAT_AD=DTYPE_2=AD                // <ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
                                                    // BUT in VarScan conflicting: <ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
#pragma GENDICT FORMAT_ADALL=DTYPE_2=ADALL          // from GIAB: <ID=ADALL,Number=R,Type=Integer,Description="Net allele depths across all datasets">
#pragma GENDICT FORMAT_ADF=DTYPE_2=ADF              // <ID=ADF,Number=R,Type=Float,Description="Allele dosage on fwd strand">
                                                    // Conflicting VarScan: <ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">
#pragma GENDICT FORMAT_ADR=DTYPE_2=ADR              // <ID=ADR,Number=R,Type=Float,Description="Allele dosage on rev strand">
                                                    // Conflicting VarScan: <ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">
#pragma GENDICT FORMAT_AF=DTYPE_2=AF                // <ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed">
#pragma GENDICT FORMAT_DP=DTYPE_2=DP                // <ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">. See also: https://gatk.broadinstitute.org/hc/en-us/articles/360036891012-DepthPerSampleHC
#pragma GENDICT FORMAT_DS=DTYPE_2=DS                // <ID=DS,Number=1,Type=Float,Description="Genotype dosage from MaCH/Thunder">. See: https://genome.sph.umich.edu/wiki/Thunder Beagle: "estimated ALT dose [P(RA) + P(AA)]"
#pragma GENDICT FORMAT_GL=DTYPE_2=GL                // <ID=GL,Number=.,Type=Float,Description="Genotype Likelihoods">
#pragma GENDICT FORMAT_GP=DTYPE_2=GP                // <ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification">
#pragma GENDICT FORMAT_GQ=DTYPE_2=GQ                // <ID=GQ,Number=1,Type=Integer,Description="Genotype Quality"> VCF spec: "conditional genotype quality, encoded as a phred quality −10log10 p(genotype call is wrong, conditioned on the site’s being variant) (Integer)"
#pragma GENDICT FORMAT_GT=DTYPE_2=GT                // <ID=GT,Number=1,Type=String,Description="Genotype">
#pragma GENDICT FORMAT_PL=DTYPE_2=PL                // <ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
#pragma GENDICT FORMAT_PLy=DTYPE_2=PLy              // Alternative PL context while testing PL_mux_by_DP (YES Mux by DP)
#pragma GENDICT FORMAT_PLn=DTYPE_2=PLn              // Alternative PL context while testing PL_mux_by_DP (NO don't Mux by DP)
#pragma GENDICT FORMAT_PRI=DTYPE_2=PRI              // <ID=PRI,Number=G,Type=Float,Description="Phred-scaled prior probabilities for genotypes">
#pragma GENDICT FORMAT_F1R2=DTYPE_2=F1R2            // <ID=F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting each allele"> see: https://github.com/broadinstitute/gatk/blob/master/docs/mutect/mutect.pdf
#pragma GENDICT FORMAT_F2R1=DTYPE_2=F2R1            // <ID=F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting each allele">
#pragma GENDICT FORMAT_MB=DTYPE_2=MB                // <ID=MB,Number=4,Type=Integer,Description="Per-sample component statistics to detect mate bias">
#pragma GENDICT FORMAT_PP=DTYPE_2=PP                // <ID=PP,Number=G,Type=Integer,Description="Phred-scaled genotype posterior probabilities rounded to the closest integer">
#pragma GENDICT FORMAT_SAC=DTYPE_2=SAC              // <ID=SAC,Number=.,Type=Integer,Description="Number of reads on the forward and reverse strand supporting each allele (including reference)">
#pragma GENDICT FORMAT_SB=DTYPE_2=SB                // <ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias">
#pragma GENDICT FORMAT_PS=DTYPE_2=PS                // seen 1: <ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
                                                    // seen 2: <ID=PS,Number=1,Type=Integer,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
#pragma GENDICT FORMAT_PSpos=DTYPE_2=PSpos          // 3 contexts used to seg PS and PID
#pragma GENDICT FORMAT_PSalt=DTYPE_2=PSalt
#pragma GENDICT FORMAT_PSref=DTYPE_2=PSref
#pragma GENDICT FORMAT_PID=DTYPE_2=PID              // <ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
#pragma GENDICT FORMAT_PGT=DTYPE_2=PGT              // <ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
#pragma GENDICT FORMAT_FL=DTYPE_2=FL                // Seen in Reich's ancient DNA datasets: <ID=FL,Number=1,Type=Character,Description="filter level in range 0-9 or no value (non-integer: N,?) with zero being least reliable; to threshold at FL=n, use all levels n-9">

#pragma GENDICT FORMAT_AB=DTYPE_2=AB                // <ID=AB,Number=1,Type=Float,Description="Allele balance for each het genotype",RendAlg="NONE">
#pragma GENDICT FORMAT_AB3=DTYPE_2=AB3              // must be the next context after AB

// PBWT fields - same dict_id for all data types using PBWT, as codec_pbwt_uncompress relies on it
#pragma GENDICT FORMAT_GT_HT=DTYPE_2=@HT    
#pragma GENDICT FORMAT_GT_HT_INDEX=DTYPE_2=@INDEXHT // different first 2 letters
#pragma GENDICT FORMAT_PBWT_RUNS=DTYPE_2=@1BWTRUN   // PBWT runs - MUST have a did_i higher that FORMAT_GT_HT's
#pragma GENDICT FORMAT_PBWT_FGRC=DTYPE_2=@2BWTFGR   // PBWT foreground run count - MUST be right after FORMAT_PBWT_RUNS

#ifndef _PBWT_RUNS
#define _PBWT_RUNS        _FORMAT_PBWT_RUNS       
#define _PBWT_FGRC        _FORMAT_PBWT_FGRC       
#define _PBWT_HT_MATRIX   _FORMAT_GT_HT           
#define _PBWT_GT_HT_INDEX _FORMAT_GT_HT_INDEX     
#endif

// INFO fields
#pragma GENDICT INFO_AC=DTYPE_1=AC                  // <ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_AF=DTYPE_1=AF                  // <ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_AN=DTYPE_1=AN                  // <ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#pragma GENDICT INFO_AA=DTYPE_1=AA                  // <ID=AA,Number=1,Type=String,Description="Ancestral Allele"> - defined in the VCF specification
#pragma GENDICT INFO_BaseCounts=DTYPE_1=BaseCounts
#pragma GENDICT INFO_DP=DTYPE_1=DP                  // <ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
#pragma GENDICT INFO_DP4=DTYPE_1=DP4                // <ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
#pragma GENDICT INFO_SF=DTYPE_1=SF                  // <ID=SF,Number=.,Type=String,Description="Source File (index to sourceFiles, f when filtered)">
#pragma GENDICT INFO_VQSLOD=DTYPE_1=VQSLOD          // <ID=VQSLOD,Number=1,Type=Float,Description="Log odds of being a true variant versus being false under the trained Gaussian mixture model">
#pragma GENDICT INFO_MQ=DTYPE_1=MQ                  // <ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">    
#pragma GENDICT INFO_MQ0=DTYPE_1=MQ0                // <ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads"> (VCF spec: "Number of MAPQ == 0 reads covering this record")
#pragma GENDICT INFO_CIPOS=DTYPE_1=CIPOS            // <ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
#pragma GENDICT INFO_CIEND=DTYPE_1=CIEND            // <ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
#pragma GENDICT INFO_NS=DTYPE_1=NS                  // <ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">

#pragma GENDICT INFO_LDAF=DTYPE_1=LDAF              // <ID=LDAF,Number=1,Type=Float,Description="MLE Allele Frequency Accounting for LD">
#pragma GENDICT INFO_AVGPOST=DTYPE_1=AVGPOST        // <ID=AVGPOST,Number=1,Type=Float,Description="Average posterior probability from MaCH/Thunder">
#pragma GENDICT INFO_RSQ=DTYPE_1=RSQ                // <ID=RSQ,Number=1,Type=Float,Description="Genotype imputation quality from MaCH/Thunder">
#pragma GENDICT INFO_ERATE=DTYPE_1=ERATE            // <ID=ERATE,Number=1,Type=Float,Description="Per-marker Mutation rate from MaCH/Thunder">
#pragma GENDICT INFO_THETA=DTYPE_1=THETA            // <ID=THETA,Number=1,Type=Float,Description="Per-marker Transition rate from MaCH/Thunder">

// Ann field defined: https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
#pragma GENDICT INFO_ANN=DTYPE_1=ANN                // <ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
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
#pragma GENDICT INFO_END=DTYPE_1=END                // <ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
#pragma GENDICT INFO_MLEAC=DTYPE_1=MLEAC            // <ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_MLEAF=DTYPE_1=MLEAF            // <ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_SOR=DTYPE_1=SOR                // <ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">. See: https://gatk.broadinstitute.org/hc/en-us/articles/360036361772-StrandOddsRatio
#pragma GENDICT INFO_QD=DTYPE_1=QD                  // <ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">. See: https://gatk.broadinstitute.org/hc/en-us/articles/360041414572-QualByDepth

#pragma GENDICT INFO_RAW_MQandDP=DTYPE_1=RAW_MQandDP// <ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
#pragma GENDICT INFO_RAW_MQandDP_MQ=DTYPE_1=R0AW_MQandDP
#pragma GENDICT INFO_RAW_MQandDP_DP=DTYPE_1=R1AW_MQandDP

#pragma GENDICT FORMAT_RGQ=DTYPE_2=RGQ              // <ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
#pragma GENDICT FORMAT_MIN_DP=DTYPE_2=MIN_DP        // <ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">

// 10xGenomics: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf
#pragma GENDICT FORMAT_BX=DTYPE_2=BX                // <ID=BX,Number=.,Type=String,Description="Barcodes and Associated Qual-Scores Supporting Alleles">
#pragma GENDICT FORMAT_PQ=DTYPE_2=PQ                // <ID=PQ,Number=1,Type=Integer,Description="Phred QV indicating probability at this variant is incorrectly phased">
#pragma GENDICT FORMAT_JQ=DTYPE_2=JQ                // <ID=JQ,Number=1,Type=Integer,Description="Phred QV indicating probability of a phasing switch error in gap prior to this variant">

// Ensembl VEP (Variant Effect Predictor) fields: https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html
#pragma GENDICT INFO_vep=DTYPE_1=vep
#pragma GENDICT INFO_AGE_HISTOGRAM_HET=DTYPE_1=AGE_HISTOGRAM_HET 
#pragma GENDICT INFO_AGE_HISTOGRAM_HOM=DTYPE_1=AGE_HISTOGRAM_HOM
#pragma GENDICT INFO_MAX_AF=DTYPE_1=MAX_AF          // highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD

#pragma GENDICT INFO_CSQ=DTYPE_1=CSQ                // <ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral">
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

// clinvar (also in dbSNP)
#pragma GENDICT INFO_ALLELEID=DTYPE_1=ALLELEID      // <ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
#pragma GENDICT INFO_CLNDN=DTYPE_1=CLNDN            // <ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
#pragma GENDICT INFO_CLNHGVS=DTYPE_1=CLNHGVS        // <ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
#pragma GENDICT INFO_CLNVI=DTYPE_1=CLNVI            // <ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
#pragma GENDICT INFO_CLNORIGIN=DTYPE_1=CLNORIGIN    // <ID=CLNORIGIN,Number=.,Type=String,Description="Allele Origin. One or more of the following values may be summed: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other">
#pragma GENDICT INFO_CLNSIG=DTYPE_1=CLNSIG          // <ID=CLNSIG,Number=.,Type=String,Description="Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 8 - confers sensitivity, 9 - risk-factor, 10 - association, 11 - protective, 12 - conflict, 13 - affects, 255 - other">
#pragma GENDICT INFO_CLNDISDB=DTYPE_1=CLNDISDB      // <ID=CLNDISDB,Number=.,Type=String,Description="Variant disease database name and ID, separated by colon (:)">
#pragma GENDICT INFO_CLNREVSTAT=DTYPE_1=CLNREVSTAT  // <ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar Review Status: no_assertion - No asserition provided by submitter, no_criteria - No assertion criteria provided by submitter, single - Classified by single submitter, mult - Classified by multiple submitters, conf - Criteria provided conflicting interpretations, exp - Reviewed by expert panel, guideline - Practice guideline">
#pragma GENDICT INFO_CLNACC=DTYPE_1=CLNACC          // <ID=CLNACC,Number=.,Type=String,Description="For each allele (comma delimited), this is a pipe-delimited list of the Clinvar RCV phenotype accession.version strings associated with that allele.">

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
#pragma GENDICT INFO_DP_HIST=DTYPE_1=DP_HIST        // from ExAC: Depth (DP) histogram in 20 equal intervals between 0-100 : See https://www.biorxiv.org/content/biorxiv/suppl/2015/10/30/030338.DC1/030338-1.pdf
#pragma GENDICT INFO_GQ_HIST=DTYPE_1=GQ_HIST        // from ExAC: Genotype Quality (GQ) histogram in 20 equal intervals between 0-100

// Structural variants (also uses INFO/END): https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/VCF%20(Variant%20Call%20Format)%20version%204.0/encoding-structural-variants/
#pragma GENDICT INFO_SVLEN=DTYPE_1=SVLEN
#pragma GENDICT INFO_SVTYPE=DTYPE_1=SVTYPE

// VarScan FORMAT and INFO fields: http://varscan.sourceforge.net/using-varscan.html
#pragma GENDICT FORMAT_RDF=DTYPE_2=RDF              // <ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
#pragma GENDICT FORMAT_RDR=DTYPE_2=RDR              // <ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
#pragma GENDICT FORMAT_SDP=DTYPE_2=SDP              // <ID=SDP,Number=1,Type=Integer,Description="Raw Read Depth as reported by SAMtools">
#pragma GENDICT FORMAT_RD=DTYPE_2=RD                // 

#pragma GENDICT FORMAT_FREQ=DTYPE_2=FREQ            // <ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
#pragma GENDICT FORMAT_PVAL=DTYPE_2=PVAL            // <ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
#pragma GENDICT FORMAT_RBQ=DTYPE_2=RBQ              // <ID=RBQ,Number=1,Type=Integer,Description="Average quality of reference-supporting bases (qual1)">
#pragma GENDICT FORMAT_ABQ=DTYPE_2=ABQ              // <ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">
#pragma GENDICT INFO_ADP=DTYPE_1=ADP                // <ID=ADP,Number=1,Type=Integer,Description="Average per-sample depth of bases with Phred score >= 0">
#pragma GENDICT INFO_WT=DTYPE_1=WT                  // <ID=WT,Number=1,Type=Integer,Description="Number of samples called reference (wild-type)">
#pragma GENDICT INFO_HET=DTYPE_1=HET                // <ID=HET,Number=1,Type=Integer,Description="Number of samples called heterozygous-variant">
#pragma GENDICT INFO_HOM=DTYPE_1=HOM                // <ID=HOM,Number=1,Type=Integer,Description="Number of samples called homozygous-variant">
#pragma GENDICT INFO_NC=DTYPE_1=NC                  // <ID=NC,Number=1,Type=Integer,Description="Number of samples not called">

// dbSNP. See: https://github.com/ncbi/dbsnp/blob/master/tutorials/Variation%20Services/Jupyter_Notebook/Data/test_vcf.vcf
// also: https://vatlab.github.io/vat-docs/documentation/keyconcepts/
#pragma GENDICT INFO_RS=DTYPE_1=RS                  // <ID=RS,Number=1,Type=Integer,Description="dbSNP ID (i.e. rs number)">
#pragma GENDICT INFO_RSPOS=DTYPE_1=RSPOS            // <ID=RSPOS,Number=1,Type=Integer,Description="Chr position reported in dbSNP">
#pragma GENDICT INFO_TOPMED=DTYPE_1=TOPMED          // <ID=TOPMED,Number=.,Type=String,Description="An ordered, comma delimited list of allele frequencies based on TOPMed, starting with the reference allele followed by alternate alleles as ordered in the ALT column. The TOPMed minor allele is the second largest value in the list.">
#pragma GENDICT INFO_GENEINFO=DTYPE_1=GENEINFO      // <ID=GENEINFO,Number=1,Type=String,Description="Pairs each of gene symbol:gene id.  The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
#pragma GENDICT INFO_dbSNPBuildID=DTYPE_1=dbSNPBuildID     // <ID=dbSNPBuildID,Number=1,Type=Integer,Description="First dbSNP Build for RS">
#pragma GENDICT INFO_PSEUDOGENEINFO=DTYPE_1=PSEUDOGENEINFO // <ID=PSEUDOGENEINFO,Number=1,Type=String,Description="Pairs each of pseudogene symbol:gene id.  The pseudogene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
#pragma GENDICT INFO_SAO=DTYPE_1=SAO                // <ID=SAO,Number=1,Type=Integer,Description="Variant Allele Origin: 0 - unspecified, 1 - Germline, 2 - Somatic, 3 - Both">
#pragma GENDICT INFO_SSR=DTYPE_1=SSR                // <ID=SSR,Number=1,Type=Integer,Description="Variant Suspect Reason Codes (may be more than one value added together) 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other">
#pragma GENDICT INFO_VC=DTYPE_1=VC                  // <ID=VC,Number=1,Type=String,Description="Variation Class">
#pragma GENDICT INFO_PM=DTYPE_1=PM                  // <ID=PM,Number=0,Type=Flag,Description="Variant has associated publication">
#pragma GENDICT INFO_NSF=DTYPE_1=NSF                // <ID=NSF,Number=0,Type=Flag,Description="Has non-synonymous frameshift A coding region variation where one allele in the set changes all downstream amino acids. FxnClass = 44">
#pragma GENDICT INFO_NSM=DTYPE_1=NSM                // <ID=NSM,Number=0,Type=Flag,Description="Has non-synonymous missense A coding region variation where one allele in the set changes protein peptide. FxnClass = 42">
#pragma GENDICT INFO_NSN=DTYPE_1=NSN                // <ID=NSN,Number=0,Type=Flag,Description="Has non-synonymous nonsense A coding region variation where one allele in the set changes to STOP codon (TER). FxnClass = 41">
#pragma GENDICT INFO_SYN=DTYPE_1=SYN                // <ID=SYN,Number=0,Type=Flag,Description="Has synonymous A coding region variation where one allele in the set does not change the encoded amino acid. FxnCode = 3">
#pragma GENDICT INFO_U3=DTYPE_1=U3                  // <ID=U3,Number=0,Type=Flag,Description="In 3' UTR Location is in an untranslated region (UTR). FxnCode = 53">
#pragma GENDICT INFO_U5=DTYPE_1=U5                  // <ID=U5,Number=0,Type=Flag,Description="In 5' UTR Location is in an untranslated region (UTR). FxnCode = 55">
#pragma GENDICT INFO_ASS=DTYPE_1=ASS                // <ID=ASS,Number=0,Type=Flag,Description="In acceptor splice site FxnCode = 73">
#pragma GENDICT INFO_DSS=DTYPE_1=DSS                // <ID=DSS,Number=0,Type=Flag,Description="In donor splice-site FxnCode = 75">
#pragma GENDICT INFO_INT=DTYPE_1=INT                // <ID=INT,Number=0,Type=Flag,Description="In Intron FxnCode = 6">
#pragma GENDICT INFO_R3=DTYPE_1=R3                  // <ID=R3,Number=0,Type=Flag,Description="In 3' gene region FxnCode = 13">
#pragma GENDICT INFO_R5=DTYPE_1=R5                  // <ID=R5,Number=0,Type=Flag,Description="In 5' gene region FxnCode = 15">
#pragma GENDICT INFO_GNO=DTYPE_1=GNO                // <ID=GNO,Number=0,Type=Flag,Description="Genotypes available.">
#pragma GENDICT INFO_PUB=DTYPE_1=PUB                // <ID=PUB,Number=0,Type=Flag,Description="RefSNP or associated SubSNP is mentioned in a publication">
#pragma GENDICT INFO_FREQ=DTYPE_1=FREQ              // <ID=FREQ,Number=.,Type=String,Description="An ordered list of allele frequencies as reported by various genomic studies, starting with the reference allele followed by alternate alleles as ordered in the ALT column. When not already in the dbSNP allele set, alleles from the studies are added to the ALT column.  The minor allele, which was previuosly reported in VCF as the GMAF, is the second largest value in the list.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter">
#pragma GENDICT INFO_COMMON=DTYPE_1=COMMON          // <ID=COMMON,Number=0,Type=Flag,Description="RS is a common SNP.  A common SNP is one that has at least one 1000Genomes population with a minor allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.">
#pragma GENDICT INFO_VP=DTYPE_1=VP                  // <ID=VP,Number=1,Type=String,Description="Variation Property.  Documentation is at ftp://ftp.ncbi.nlm.nih.gov/snp/specs/dbSNP_BitField_latest.pdf">
#pragma GENDICT INFO_CAF=DTYPE_1=CAF                // <ID=CAF,Number=.,Type=String,Description="An ordered, comma delimited list of allele frequencies based on 1000Genomes, starting with the reference allele followed by alternate alleles as ordered in the ALT column. Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previuosly reported in VCF as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter"> 
#pragma GENDICT INFO_G5A=DTYPE_1=G5A                // >5% minor allele frequency in each and all populations
#pragma GENDICT INFO_G5=DTYPE_1=G5                  // >5% minor allele frequency in 1+ populations

// Illumina Genotyping. See explanation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4751548/
#pragma GENDICT INFO_PROBE_A=DTYPE_1=PROBE_A        // <ID=PROBE_A,Number=1,Type=String,Description="Probe base pair sequence">
#pragma GENDICT INFO_PROBE_B=DTYPE_1=PROBE_B        // <ID=PROBE_B,Number=1,Type=String,Description="Probe base pair sequence; not missing for strand-ambiguous SNPs">
#pragma GENDICT INFO_ALLELE_A=DTYPE_1=ALLELE_A      // <ID=ALLELE_A,Number=1,Type=String,Description="A allele">
#pragma GENDICT INFO_ALLELE_B=DTYPE_1=ALLELE_B      // <ID=ALLELE_B,Number=1,Type=String,Description="B allele">
#pragma GENDICT INFO_refSNP=DTYPE_1=refSNP          // <ID=refSNP,Number=1,Type=String,Description="dbSNP rsID">
#pragma GENDICT INFO_ILLUMINA_CHR=DTYPE_1=ILLUMINA_CHR      // <ID=ILLUMINA_CHR,Number=1,Type=String,Description="Chromosome in Illumina manifest">
#pragma GENDICT INFO_ILLUMINA_POS=DTYPE_1=ILLUMINA_POS      // <ID=ILLUMINA_POS,Number=1,Type=Integer,Description="Position in Illumina manifest">
#pragma GENDICT INFO_ILLUMINA_STRAND=DTYPE_1=ILLUMINA_STRAND// <ID=ILLUMINA_STRAND,Number=1,Type=String,Description="Probe strand">
#pragma GENDICT FORMAT_BAF=DTYPE_2=BAF              // <ID=BAF,Number=1,Type=Float,Description="B Allele Frequency">
#pragma GENDICT FORMAT_X=DTYPE_2=X                  // <ID=X,Number=1,Type=Integer,Description="Raw X intensity">
#pragma GENDICT FORMAT_Y=DTYPE_2=Y                  // <ID=Y,Number=1,Type=Integer,Description="Raw Y intensity">

// infiniumFinalReportConverter
#pragma GENDICT INFO_INFINIUM_CR=DTYPE_1=CR         // <ID=CR,Number=.,Type=Float,Description="SNP Callrate"
#pragma GENDICT INFO_INFINIUM_GentrainScore=DTYPE_1=GentrainScore // <ID=GentrainScore,Number=.,Type=Float,Description="Gentrain Score">
#pragma GENDICT INFO_INFINIUM_HW=DTYPE_1=HW         // <ID=HW,Number=.,Type=Float,Description="Hardy-Weinberg Equilibrium">

// Beagle. See: https://faculty.washington.edu/browning/beagle/beagle_5.0_26May18.pdf
#pragma GENDICT INFO_AR2=DTYPE_1=AR2                // <ID=AR2,Number=1,Type=Float,Description="Allelic R-Squared: estimated squared correlation between most probable REF dose and true REF dose">
#pragma GENDICT INFO_DR2=DTYPE_1=DR2                // <ID=DR2,Number=1,Type=Float,Description="Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose">
#pragma GENDICT INFO_IMP=DTYPE_1=IMP                // <ID=IMP,Number=0,Type=Flag,Description="Imputed marker">

// Genocove loimpute fields: https://docs.gencove.com/main/faq/
#pragma GENDICT FORMAT_RC=DTYPE_2=RC                // <ID=RC,Number=1,Type=Integer,Description="Count of reads with REF allele">
#pragma GENDICT FORMAT_AC=DTYPE_2=AC                // <ID=AC,Number=1,Type=Integer,Description="Count of reads with ALT allele">

// Gwas2VCF: https://github.com/MRCIEU/gwas-vcf-specification
#pragma GENDICT INFO_RSID=DTYPE_1=RSID              // <ID=RSID,Number=1,Type=String,Description="dbSNP identifier">
#pragma GENDICT FORMAT_NS=DTYPE_2=NS                // <ID=NS,Number=A,Type=Float,Description="Variant-specific number of samples/individuals with called genotypes used to test association with specified trait">
#pragma GENDICT FORMAT_EZ=DTYPE_2=EZ                // <ID=EZ,Number=A,Type=Float,Description="Z-score provided if it was used to derive the ES and SE fields">
#pragma GENDICT FORMAT_SI=DTYPE_2=SI                // <ID=SI,Number=A,Type=Float,Description="Accuracy score of summary association statistics imputation">
#pragma GENDICT FORMAT_NC=DTYPE_2=NC                // <ID=NC,Number=A,Type=Float,Description="Variant-specific number of cases used to estimate genetic effect (binary traits only)">
#pragma GENDICT FORMAT_ES=DTYPE_2=ES                // <ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">
#pragma GENDICT FORMAT_SE=DTYPE_2=SE                // <ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">
#pragma GENDICT FORMAT_LP=DTYPE_2=LP                // <ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">
#pragma GENDICT FORMAT_ID=DTYPE_2=ID                // (v1 of GWAS-VCF)  <ID=ID,Number=1,Type=String,Description="Study variant identifier">
//#pragma GENDICT FORMAT_AF=DTYPE_2=AF              // (overlap) <ID=AF,Number=A,Type=Float,Description="Alternative allele frequency in trait subset">
//#pragma GENDICT FORMAT_AC=DTYPE_2=AC              // (overlap) <ID=AC,Number=A,Type=Float,Description="Alternative allele count in the trait subset">

// Genozip INFO fields
#pragma GENDICT INFO_LUFT=DTYPE_1=LUFT
#pragma GENDICT INFO_PRIM=DTYPE_1=PRIM
#pragma GENDICT INFO_LREJ=DTYPE_1=Lrej
#pragma GENDICT INFO_PREJ=DTYPE_1=Prej

#define VCF_MAX_PLOIDY 100  // set to a reasonable 100 to avoid memory allocation explosion in case of an error in the VCF file
#if VCF_MAX_PLOIDY > 255
#error "VCF_MAX_PLOIDY cannot go beyond 255 because Ploidy is uint8_t"
#endif

#define MAX_ALLELES 100 // REF (allele #0) + 99 ALTs (alleles # 1-99)
typedef uint8_t Allele; // elements of ht_matrix: values 48->147 for allele 0 to 99, '*' for unused, '%', '-'

// ZIP stuff
extern void vcf_zip_initialize (void);
extern void vcf_zip_finalize (bool is_last_user_txt_file);
extern void vcf_zip_genozip_header (SectionHeaderGenozipHeader *header);
extern void vcf_zip_init_vb (VBlockP vb);
extern void vcf_liftover_display_lift_report (void);
extern void vcf_zip_after_compress (VBlockP vb);
extern void vcf_zip_after_vbs (void);
extern void vcf_zip_set_txt_header_specific (SectionHeaderTxtHeader *txt_header);
extern void vcf_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeader *vb_header);
extern bool vcf_zip_vb_has_count (VBlockP vb);
extern void vcf_zip_generate_recon_plan (void);
extern void vcf_zip_update_txt_counters (VBlockP vb);
extern bool is_vcf (STRp(header), bool *need_more);

// SEG stuff
extern rom vcf_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void vcf_seg_initialize (VBlockP vb_);
extern void vcf_zip_after_compute (VBlockP vb);
extern void vcf_seg_finalize (VBlockP vb_);
extern bool vcf_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern bool vcf_seg_is_big (ConstVBlockP vb, DictId dict_id);
extern TranslatorId vcf_lo_luft_trans_id (DictId dict_id, char number);
extern uint32_t vcf_seg_get_vb_recon_size (VBlockP vb);

// PIZ stuff
extern void vcf_piz_genozip_header (const SectionHeaderGenozipHeader *header);
extern bool vcf_piz_maybe_reorder_lines (void);
extern bool vcf_piz_init_vb (VBlockP vb, const SectionHeaderVbHeader *header, uint32_t *txt_data_so_far_single_0_increment);
extern void vcf_piz_recon_init (VBlockP vb);
extern IS_SKIP (vcf_piz_is_skip_section);
extern CONTAINER_FILTER_FUNC (vcf_piz_filter);
extern CONTAINER_CALLBACK (vcf_piz_container_cb);
extern CONTAINER_ITEM_CALLBACK (vcf_piz_con_item_cb);

// VCF Header stuff
extern void vcf_header_piz_init (void);
extern bool vcf_inspect_txt_header (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);
extern uint32_t vcf_header_get_num_samples (void);
extern bool vcf_header_get_has_fileformat (void);
extern void vcf_piz_finalize (void);

// VBlock stuff
extern void vcf_vb_release_vb();
extern void vcf_vb_destroy_vb();
extern void vcf_vb_cleanup_memory();
extern unsigned vcf_vb_size (DataType dt);
extern unsigned vcf_vb_zip_dl_size (void);
extern void vcf_reset_line (VBlockP vb);
extern bool vcf_vb_has_haplotype_data (VBlockP vb);
extern bool vcf_vb_is_primary (VBlockP vb);
extern bool vcf_vb_is_luft (VBlockP vb);
extern int32_t vcf_vb_get_reject_bytes (VBlockP vb);
extern rom vcf_coords_name (int coord);

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
extern void vcf_samples_add  (rom samples_str);

#define VCF_SPECIAL { vcf_piz_special_main_REFALT, vcf_piz_special_FORMAT, vcf_piz_special_INFO_AC, vcf_piz_special_INFO_SVLEN, \
                      vcf_piz_special_FORMAT_DS_old, vcf_piz_special_INFO_BaseCounts, vcf_piz_special_INFO_SF, piz_special_MINUS,  \
                      vcf_piz_special_LIFT_REF, vcf_piz_special_COPYSTAT, vcf_piz_special_other_REFALT, vcf_piz_special_COPYPOS, vcf_piz_special_ALLELE, \
                      vcf_piz_special_INFO_HGVS_SNP_POS, vcf_piz_special_INFO_HGVS_SNP_REFALT, \
                      vcf_piz_special_INFO_HGVS_DEL_END_POS, vcf_piz_special_INFO_HGVS_DEL_PAYLOAD, \
                      vcf_piz_special_INFO_HGVS_INS_END_POS, vcf_piz_special_INFO_HGVS_INS_PAYLOAD, \
                      vcf_piz_special_INFO_HGVS_DELINS_END_POS, vcf_piz_special_INFO_HGVS_DELINS_PAYLOAD,\
                      vcf_piz_special_MUX_BY_DOSAGE, vcf_piz_special_FORMAT_AB, vcf_piz_special_FORMAT_GQ, \
                      vcf_piz_special_MUX_BY_DOSAGExDP, vcf_piz_special_COPY_REForALT, vcf_piz_special_DP_by_DP_v13, \
                      vcf_piz_special_PS_by_PID, vcf_piz_special_PGT, vcf_piz_special_DP_by_DP, vcf_piz_special_DP_by_DP_single,\
                      vcf_piz_special_RGQ, vcf_piz_special_MUX_BY_HAS_RGQ, vcf_piz_special_SVTYPE, \
                      vcf_piz_special_ALLELE_A, vcf_piz_special_ALLELE_B, vcf_piz_special_MUX_BY_ADJ_DOSAGE,\
                      vcf_piz_special_PROBE_A, vcf_piz_special_PROBE_B, vcf_piz_special_QD, \
                      vcf_piz_special_MUX_BY_VARTYPE }

SPECIAL (VCF, 0,  main_REFALT,         vcf_piz_special_main_REFALT);
SPECIAL (VCF, 1,  FORMAT,              vcf_piz_special_FORMAT)
SPECIAL (VCF, 2,  AC,                  vcf_piz_special_INFO_AC);
SPECIAL (VCF, 3,  SVLEN,               vcf_piz_special_INFO_SVLEN);
SPECIAL (VCF, 4,  DS_old,              vcf_piz_special_FORMAT_DS_old);            // used in files up to 12.0.42
SPECIAL (VCF, 5,  BaseCounts,          vcf_piz_special_INFO_BaseCounts);
SPECIAL (VCF, 6,  SF,                  vcf_piz_special_INFO_SF);
SPECIAL (VCF, 7,  MINUS,               piz_special_MINUS);                        // added v12.0.0 
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
SPECIAL (VCF, 21, MUX_BY_DOSAGE,       vcf_piz_special_MUX_BY_DOSAGE);            // added v13.0.0
SPECIAL (VCF, 22, AB,                  vcf_piz_special_FORMAT_AB);                // added v13.0.0
SPECIAL (VCF, 23, GQ,                  vcf_piz_special_FORMAT_GQ);                // added v13.0.0
SPECIAL (VCF, 24, MUX_BY_DOSAGExDP,    vcf_piz_special_MUX_BY_DOSAGExDP);         // added v13.0.3
SPECIAL (VCF, 25, COPY_REForALT,       vcf_piz_special_COPY_REForALT);            // added v13.0.5
SPECIAL (VCF, 26, DP_by_DP_v13,        vcf_piz_special_DP_by_DP_v13);             // added v13.0.5, removed in v14
SPECIAL (VCF, 27, PS_BY_PID,           vcf_piz_special_PS_by_PID);                // added v13.0.11
SPECIAL (VCF, 28, PGT,                 vcf_piz_special_PGT);                      // added v14.0.0
SPECIAL (VCF, 29, DP_by_DP,            vcf_piz_special_DP_by_DP);                 // added v14.0.0 - multiple samples: INFO/DP by sum(FORMAT/DP)
SPECIAL (VCF, 30, DP_by_DP_single,     vcf_piz_special_DP_by_DP_single);          // added v14.0.0 - single sample: FORMAT/DP by INFO/DP
SPECIAL (VCF, 31, RGQ,                 vcf_piz_special_RGQ);                      // added v14.0.0
SPECIAL (VCF, 32, MUX_BY_HAS_RGQ,      vcf_piz_special_MUX_BY_HAS_RGQ);           // added v14.0.0
SPECIAL (VCF, 33, SVTYPE,              vcf_piz_special_SVTYPE);                   // added v14.0.12
SPECIAL (VCF, 34, ALLELE_A,            vcf_piz_special_ALLELE_A);                 // added v14.0.12
SPECIAL (VCF, 35, ALLELE_B,            vcf_piz_special_ALLELE_B);                 // added v14.0.12
SPECIAL (VCF, 36, MUX_BY_ADJ_DOSAGE,   vcf_piz_special_MUX_BY_ADJ_DOSAGE);        // added v14.0.12
SPECIAL (VCF, 37, PROBE_A,             vcf_piz_special_PROBE_A);                  // added v14.0.12
SPECIAL (VCF, 38, PROBE_B,             vcf_piz_special_PROBE_B);                  // added v14.0.12
SPECIAL (VCF, 39, QD,                  vcf_piz_special_QD);                       // added v14.0.17
SPECIAL (VCF, 40, MUX_BY_VARTYPE,      vcf_piz_special_MUX_BY_VARTYPE);           // added v15.0.0
#define NUM_VCF_SPECIAL 41

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
TRANSLATOR (VCF, VCF,   11, NEG,    vcf_piz_luft_NEG)      

#define NUM_VCF_TRANS   12 // including "none"
#define VCF_TRANSLATORS { NULL /* none */, vcf_piz_luft_G, vcf_piz_luft_R, vcf_piz_luft_R2, vcf_piz_luft_A_AN,         \
                          vcf_piz_luft_A_1, vcf_piz_luft_PLOIDY, vcf_piz_luft_GT, vcf_piz_luft_END, vcf_piz_luft_XREV, \
                          vcf_piz_luft_ALLELE, vcf_piz_luft_NEG }

typedef struct {
    rom alg_name;
    enum { TW_NEVER, TW_ALWAYS, TW_REF_ALT_SWITCH, TW_XSTRAND } upon;
} LuftTransLateProp;

// names of INFO / FORMAT algorithms, goes into VCF header's ##INFO / ##FORMAT "RendAlg" attribute
                           /* Algorithm   Trigger          */
#define DVCF_TRANS_PROPS { { "NONE",      TW_NEVER          },   /* never translate */\
                           { "G",         TW_REF_ALT_SWITCH },   /* reshuffle a  'G' vector (one element per genotype) if REF⇆ALT changed */\
                           { "R",         TW_REF_ALT_SWITCH },   /* reshuffle an 'R' vector (one element per ref/alt allele) if REF⇆ALT changed */\
                           { "R2",        TW_REF_ALT_SWITCH },   /* reshuffle a vector with 2 elements per ref/alt allele, if REF⇆ALT changed */\
                           { "A_AN",      TW_REF_ALT_SWITCH },   /* recalculate an 'A' vector (one element per ALT allele) if REF⇆ALT changed, who's elements, including a missing element for REF, add up to AN (example: AC). */ \
                           { "A_1",       TW_REF_ALT_SWITCH },   /* recalculate an 'A' vector (one element per ALT allele) if REF⇆ALT changed, who's elements, including a missing element for REF, add up to 1 (example: AF). */ \
                           { "PLOIDY",    TW_REF_ALT_SWITCH },   /* recalculate a float to (ploidy-value) */ \
                           { "GT",        TW_REF_ALT_SWITCH },   /* recalculate the allele numbers FORMAT/GT if REF⇆ALT changed */ \
                           { "END",       TW_ALWAYS         },   /* recalculate INFO/END */\
                           { "XREV",      TW_XSTRAND        },   /* reverse the elements of a vector if XSTRAND. Example: INFO/BaseCounts */\
                           { "ALLELE",    TW_ALWAYS         },   /* copy an allele verbatim including if changes or changes order. example: INFO/AA */ \
                           { "NEG",       TW_REF_ALT_SWITCH } }  /* negate numeric value if REF⇆ALT */
extern const LuftTransLateProp ltrans_props[NUM_VCF_TRANS];

#define needs_translation(ctx)  (z_is_dvcf && (ctx)->luft_trans && \
    ((ltrans_props[(ctx)->luft_trans].upon == TW_REF_ALT_SWITCH && LO_IS_OK_SWITCH (last_ostatus)) || \
     (ltrans_props[(ctx)->luft_trans].upon == TW_ALWAYS         && LO_IS_OK (last_ostatus))        || \
     (ltrans_props[(ctx)->luft_trans].upon == TW_XSTRAND        && LO_IS_OK (last_ostatus) && *CTX(VCF_oXSTRAND)->last_snip != '-')))

#define VCF_DICT_ID_ALIASES                 \
    /*         alias        maps to     */  \
    { DT_VCF,  _INFO_END,   _VCF_POS    },  \
    { DT_VCF,  _INFO_CIEND, _INFO_CIPOS },  \

#define dict_id_is_vcf_info_sf   dict_id_is_type_1
#define dict_id_is_vcf_format_sf dict_id_is_type_2

typedef enum { VCF_COMP_MAIN, VCF_COMP_PRIM_ONLY, VCF_COMP_LUFT_ONLY } VcfComponentType;
#define VCF_COMP_NAMES { "MAIN", "PRIM", "LUFT" }
