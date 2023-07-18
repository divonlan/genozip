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
#pragma GENDICT FORMAT_PSpos=DTYPE_2=PSpos          // 3 contexts used to seg PS and PID (must be in this order, consecutive)
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
#pragma GENDICT INFO_NS=DTYPE_1=NS                  // <ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">

#pragma GENDICT INFO_LDAF=DTYPE_1=LDAF              // <ID=LDAF,Number=1,Type=Float,Description="MLE Allele Frequency Accounting for LD">
#pragma GENDICT INFO_AVGPOST=DTYPE_1=AVGPOST        // <ID=AVGPOST,Number=1,Type=Float,Description="Average posterior probability from MaCH/Thunder">
#pragma GENDICT INFO_RSQ=DTYPE_1=RSQ                // <ID=RSQ,Number=1,Type=Float,Description="Genotype imputation quality from MaCH/Thunder">
#pragma GENDICT INFO_ERATE=DTYPE_1=ERATE            // <ID=ERATE,Number=1,Type=Float,Description="Per-marker Mutation rate from MaCH/Thunder">
#pragma GENDICT INFO_THETA=DTYPE_1=THETA            // <ID=THETA,Number=1,Type=Float,Description="Per-marker Transition rate from MaCH/Thunder">

// SnpEFF
#pragma GENDICT INFO_ANN=DTYPE_1=ANN                // <ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">. See: https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
#pragma GENDICT INFO_ANN_Allele=DTYPE_1=@ANN_Allele // note: nice looking tag, as might be displayed in INFO/Lref

#pragma GENDICT INFO_EFF=DTYPE_1=EFF                // See: https://pcingola.github.io/SnpEff/se_inputoutput/#eff-field-vcf-output-files
 	 
// Added by GATK HaplotypeCaller in a gVCF: https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
#pragma GENDICT INFO_END=DTYPE_1=END                // <ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
#pragma GENDICT INFO_MLEAC=DTYPE_1=MLEAC            // <ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_MLEAF=DTYPE_1=MLEAF            // <ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">

// GATK Hard-filtering germline short variants : https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
#pragma GENDICT INFO_QD=DTYPE_1=QD                  // <ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">. See: https://gatk.broadinstitute.org/hc/en-us/articles/360041414572-QualByDepth
#pragma GENDICT INFO_FS=DTYPE_1=FS                  // Strand bias estimated using Fisher's exact test : https://gatk.broadinstitute.org/hc/en-us/articles/360037592371-FisherStrand
#pragma GENDICT INFO_SOR=DTYPE_1=SOR                // <ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">. See: https://gatk.broadinstitute.org/hc/en-us/articles/360036361772-StrandOddsRatio
#pragma GENDICT INFO_MQRankSum=DTYPE_1=MQRankSum    // Rank sum test for mapping qualities of REF versus ALT reads : https://gatk.broadinstitute.org/hc/en-us/articles/360037426091-MappingQualityRankSumTest
#pragma GENDICT INFO_ReadPosRankSum=DTYPE_1=ReadPosRankSum         // Rank sum test for relative positioning of REF versus ALT alleles within read : https://gatk.broadinstitute.org/hc/en-us/articles/360036367052-ReadPosRankSumTest

#pragma GENDICT INFO_HaplotypeScore=DTYPE_1=HaplotypeScore // <ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
#pragma GENDICT INFO_ExcessHet=DTYPE_1=ExcessHet    // Phred-scaled p-value for exact test of excess heterozygosity : https://gatk.broadinstitute.org/hc/en-us/articles/360037261311-ExcessHet
#pragma GENDICT INFO_AS_QD=DTYPE_1=AS_QD            // Allele-specific call confidence normalized by depth of sample reads supporting the allele : https://gatk.broadinstitute.org/hc/en-us/articles/360036485832-AS-QualByDepth
#pragma GENDICT INFO_BaseQRankSum=DTYPE_1=BaseQRankSum             // Rank sum test of REF versus ALT base quality scores : https://gatk.broadinstitute.org/hc/en-us/articles/360036863231-BaseQualityRankSumTest
#pragma GENDICT INFO_InbreedingCoeff=DTYPE_1=InbreedingCoeff       // Likelihood-based test for the consanguinity among samples : https://gatk.broadinstitute.org/hc/en-us/articles/360036351032-InbreedingCoeff
#pragma GENDICT INFO_AS_InbreedingCoeff=DTYPE_1=AS_InbreedingCoeff // Allele-specific likelihood-based test for the consanguinity among samples : https://gatk.broadinstitute.org/hc/en-us/articles/360036827291-AS-InbreedingCoeff
#pragma GENDICT INFO_RAW_MQ=DTYPE_1=RAW_MQ          // <ID=RAW_MQ,Number=1,Type=Float,Description="Raw data for RMS Mapping Quality">

#pragma GENDICT INFO_RAW_MQandDP=DTYPE_1=RAW_MQandDP// <ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
#pragma GENDICT INFO_RAW_MQandDP_MQ=DTYPE_1=R0AW_MQandDP
#pragma GENDICT INFO_RAW_MQandDP_DP=DTYPE_1=R1AW_MQandDP

#pragma GENDICT FORMAT_RGQ=DTYPE_2=RGQ              // <ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
#pragma GENDICT FORMAT_MIN_DP=DTYPE_2=MIN_DP        // <ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">

// DRAGEN gVCF, see: https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/Homeref_Blocks_Format_fDG.htm
#pragma GENDICT FORMAT_SPL=DTYPE_2=SPL              // <ID=SPL,Number=.,Type=Integer,Description="Normalized, Phred-scaled likelihoods for SNPs based on the reference confidence model">
#pragma GENDICT FORMAT_ICNT=DTYPE_2=ICNT            // <ID=ICNT,Number=2,Type=Integer,Description="Counts of INDEL informative reads based on the reference confidence model">

// Illumina IsaacVariantCaller (discontinued) : https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/basespace/isaac-wgs-user-guide-15050954b.pdf
// Also: https://github.com/sequencing/isaac_variant_caller
#pragma GENDICT INFO_SNVSB=DTYPE_1=SNVSB            // SNV site strand bias
#pragma GENDICT INFO_SNVHPOL=DTYPE_1=SNVHPOL        // SNV contextual homopolymer length
#pragma GENDICT INFO_CIGAR=DTYPE_1=CIGAR            // CIGAR alignment for each alternate indel allele
#pragma GENDICT INFO_RU=DTYPE_1=RU                  // Smallest repeating sequence unit extended or contracted in the indel allele relative to the reference. RUs longer than 20 bases are not reported.
#pragma GENDICT INFO_REFREP=DTYPE_1=REFREP          // Number of times RU is repeated in reference
#pragma GENDICT INFO_IDREP=DTYPE_1=IDREP            // Number of times RU is repeated in indel allele.
#pragma GENDICT INFO_BLOCKAVG_min30p3a=DTYPE_1=BLOCKAVG_min30p3a // Non-variant site block. All sites in a block are constrained to be non-variant, have the same filter value, and have all sample values in range [x,y], y <= max(x+3,(x*1.3)). All printed site block sample values are the minimum observed in the region spanned by the block
#pragma GENDICT FORMAT_GQX=DTYPE_2=GQX              // <ID=GQX,Number=1,Type=Integer,Description="Empirically calibrated variant quality score for variant sites, otherwise Minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}">
#pragma GENDICT FORMAT_DPF=DTYPE_2=DPF              // <ID=DPF,Number=1,Type=Integer,Description="Basecalls filtered from input prior to site genotyping">
#pragma GENDICT FORMAT_DPI=DTYPE_2=DPI              // <ID=DPI,Number=1,Type=Integer,Description="Read depth associated with indel, taken from the site preceding the indel.">

// Illumina starling (looks like evolved IsaacVariantCaller): https://support.illumina.com/help/BS_App_TS_Amplicon_OLH_15055858/Content/Source/Informatics/Apps/IsaacVariantCaller_appENR.htm
// (also SNVSB, SNVHPOL, CIGAR, RU, REFREP, IDREP, BLOCKAVG_min30p3a, GQX, DPF, DPI as in IsaacVariantCaller and standard tags)
#pragma GENDICT INFO_cosmic=DTYPE_1=cosmic          // <ID=cosmic,Number=.,Type=String,Description="The numeric identifier for the variant in the Catalogue of Somatic Mutations in Cancer (COSMIC) database. Format: GenotypeIndex|Significance">
#pragma GENDICT INFO_phyloP=DTYPE_1=phyloP          // <ID=phyloP,Number=A,Type=Float,Description="PhyloP conservation score. Denotes how conserved the reference sequence is between species throughout evolution">
#pragma GENDICT INFO_AF1000G=DTYPE_1=AF1000G        // <ID=AF1000G,Number=A,Type=Float,Description="The allele frequency from all populations of 1000 genomes data">
//#pragma GENDICT INFO_AA=DTYPE_1=AA                // (dup) <ID=AA,Number=A,Type=String,Description="The inferred allele ancestral (if determined) to the chimpanzee/human lineage.">
#pragma GENDICT INFO_GMAF=DTYPE_1=GMAF              // <ID=GMAF,Number=A,Type=String,Description="Global minor allele frequency (GMAF); technically, the frequency of the second most frequent allele.  Format: GlobalMinorAllele|AlleleFreqGlobalMinor">
#pragma GENDICT INFO_clinvar=DTYPE_1=clinvar        // <ID=clinvar,Number=.,Type=String,Description="Clinical significance. Format: GenotypeIndex|Significance">
#pragma GENDICT INFO_EVS=DTYPE_1=EVS                // <ID=EVS,Number=A,Type=String,Description="Allele frequency, coverage and sample count taken from the Exome Variant Server (EVS). Format: AlleleFreqEVS|EVSCoverage|EVSSamples.">
#pragma GENDICT INFO_RefMinor=DTYPE_1=RefMinor      // <ID=RefMinor,Number=0,Type=Flag,Description="Denotes positions where the reference base is a minor allele and is annotated as though it were a variant">
#pragma GENDICT INFO_CSQT=DTYPE_1=CSQT              // <ID=CSQT,Number=.,Type=String,Description="Consequence type as predicted by IAE. Format: GenotypeIndex|HGNC|Transcript ID|Consequence">
#pragma GENDICT INFO_CSQR=DTYPE_1=CSQR              // <ID=CSQR,Number=.,Type=String,Description="Predicted regulatory consequence type. Format: GenotypeIndex|RegulatoryID|Consequence">
#pragma GENDICT INFO_Unphased=DTYPE_1=Unphased      // <ID=Unphased,Number=0,Type=Flag,Description="Indicates a record that is within the specified phasing window of another variant but could not be phased due to lack of minimum read support.">
#pragma GENDICT FORMAT_VF=DTYPE_2=VF                // <ID=VF,Number=1,Type=Float,Description="Variant frequency">
#
// 10xGenomics: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf
#pragma GENDICT FORMAT_BX=DTYPE_2=BX                // <ID=BX,Number=.,Type=String,Description="Barcodes and Associated Qual-Scores Supporting Alleles">
#pragma GENDICT FORMAT_PQ=DTYPE_2=PQ                // <ID=PQ,Number=1,Type=Integer,Description="Phred QV indicating probability at this variant is incorrectly phased">
#pragma GENDICT FORMAT_JQ=DTYPE_2=JQ                // <ID=JQ,Number=1,Type=Integer,Description="Phred QV indicating probability of a phasing switch error in gap prior to this variant">

// Ensembl VEP (Variant Effect Predictor) fields: https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html https://www.ensembl.org/info/docs/tools/vep/index.html
#pragma GENDICT INFO_AGE_HISTOGRAM_HET=DTYPE_1=AGE_HISTOGRAM_HET 
#pragma GENDICT INFO_AGE_HISTOGRAM_HOM=DTYPE_1=AGE_HISTOGRAM_HOM
#pragma GENDICT INFO_MAX_AF=DTYPE_1=MAX_AF          // highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD

#pragma GENDICT INFO_CSQ=DTYPE_1=CSQ                // <ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral">
#pragma GENDICT INFO_vep=DTYPE_1=vep/*in gnomAD*/   // <ID=vep,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info">

// Mastermind (also as VEP extension) : https://www.genomenon.com/wp-content/uploads/2019/09/MM-Integration-Technical-Documentation-2019-08-23.pdf
// INFO fields: GENE,HGVSG,MMCNT1,MMCNT2,MMCNT3,MMID3,MMURI3
//#pragma GENDICT INFO_GENE=DTYPE_1=GENE            // (dup) <ID=GENE,Number=.,Type=String,Description="Genes for this variant">
//#pragma GENDICT INFO_HGVSG=DTYPE_1=HGVSG          // (dup) <ID=HGVSG,Number=1,Type=String,Description="HGVS genomic notation for this variant">
#pragma GENDICT INFO_MMCNT=DTYPE_1=MMCNT           
#pragma GENDICT INFO_MMCNT1=DTYPE_1=MMCNT1          // <ID=MMCNT1,Number=1,Type=Integer,Description="Count of Mastermind articles with cDNA matches for this specific variant">
#pragma GENDICT INFO_MMCNT2=DTYPE_1=MMCNT2          // <ID=MMCNT2,Number=1,Type=Integer,Description="Count of Mastermind articles with variants either explicitly matching at the cDNA level or given only at protein level">
#pragma GENDICT INFO_MMCNT3=DTYPE_1=MMCNT3          // <ID=MMCNT3,Number=1,Type=Integer,Description="Count of Mastermind articles including other DNA-level variants resulting in the same amino acid change">
#pragma GENDICT INFO_MMID3=DTYPE_1=MMID3            // <ID=MMID3,Number=.,Type=String,Description="Mastermind variant identifiers, as gene:key, for MMCNT3">
#pragma GENDICT INFO_MMURI3=DTYPE_1=MMURI3          // <ID=MMURI3,Number=1,Type=String,Description="Mastermind search URI for articles including other DNA-level variants resulting in the same amino acid change">
#pragma GENDICT INFO_MMURI=DTYPE_1=MMURI

// ClinVar (also in dbSNP)
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
#pragma GENDICT INFO_HGVS_dup_end_pos=DTYPE_1=H6GVS_endpos_dup     
#pragma GENDICT INFO_HGVS_no_payload=DTYPE_1=H9GVS_payload_dup

// ICGC see: https://icgc-data-parser.readthedocs.io/en/master/icgc-ssm-file.html and https://docs.icgc.org/submission/guide/icgc-simple-somatic-mutation-format/
#pragma GENDICT INFO_CONSEQUENCE=DTYPE_1=CONSEQUENCE// <ID=CONSEQUENCE,Number=.,Type=String,Description="Mutation consequence predictions annotated by SnpEff (subfields: gene_symbol|gene_affected|gene_strand|transcript_name|transcript_affected|protein_affected|consequence_type|cds_mutation|aa_mutation)">
#pragma GENDICT INFO_OCCURRENCE=DTYPE_1=OCCURRENCE  // <ID=OCCURRENCE,Number=.,Type=String,Description="Mutation occurrence counts broken down by project (subfields: project_code|affected_donors|tested_donors|frequency)">
#pragma GENDICT INFO_mutation=DTYPE_1=mutation      // <ID=mutation,Number=1,Type=String,Description="Somatic mutation definition">
#pragma GENDICT INFO_studies=DTYPE_1=studies        // <ID=studies,Number=1,Type=Integer,Description="Produced from study.">
#pragma GENDICT INFO_affected_donors=DTYPE_1=affected_donors // <ID=affected_donors,Number=1,Type=Integer,Description="Number of donors with the current mutation">
#pragma GENDICT INFO_project_count=DTYPE_1=project_count     // <ID=project_count,Number=1,Type=Integer,Description="Number of projects with the current mutation">
#pragma GENDICT INFO_tested_donors=DTYPE_1=tested_donors     // <ID=tested_donors,Number=1,Type=Integer,Description="Total number of donors with SSM data available">

// ExAC (Exome Aggregation Consortium) fields: https://gnomad.broadinstitute.org/downloads#exac-variants
#pragma GENDICT INFO_DP_HIST=DTYPE_1=DP_HIST        // from ExAC: Depth (DP) histogram in 20 equal intervals between 0-100 : See https://www.biorxiv.org/content/biorxiv/suppl/2015/10/30/030338.DC1/030338-1.pdf
#pragma GENDICT INFO_GQ_HIST=DTYPE_1=GQ_HIST        // from ExAC: Genotype Quality (GQ) histogram in 20 equal intervals between 0-100

// gnomAD: https://gnomad.broadinstitute.org/downloads 
// note: pairs of fields, eg INFO_age_hist_het_bin_freq and INFO_age_hist_hom_bin_freq, map to the same dict_id
#pragma GENDICT INFO_age_hist_het_bin_freq=DTYPE_1=age_hist_het_bin_freq     // <ID=age_hist_het_bin_freq,Number=A,Type=String,Description="Histogram of ages of heterozygous individuals; bin edges are: 30.0|35.0|40.0|45.0|50.0|55.0|60.0|65.0|70.0|75.0|80.0; total number of individuals of any genotype bin: bin_edges|bin_freq|n_smaller|n_larger">
//#pragma GENDICT INFO_age_hist_hom_bin_freq=DTYPE_1=age_hist_hom_bin_freq   // <ID=age_hist_hom_bin_freq,Number=A,Type=String,Description="Histogram of ages of homozygous alternate individuals; bin edges are: 30.0|35.0|40.0|45.0|50.0|55.0|60.0|65.0|70.0|75.0|80.0; total number of individuals of any genotype bin: bin_edges|bin_freq|n_smaller|n_larger">
#pragma GENDICT INFO_gq_hist_alt_bin_freq=DTYPE_1=gq_hist_alt_bin_freq       // <ID=gq_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for GQ in heterozygous individuals calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100">
//#pragma GENDICT INFO_gq_hist_all_bin_freq=DTYPE_1=gq_hist_all_bin_freq     // <ID=dp_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for DP in heterozygous individuals calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100">
#pragma GENDICT INFO_dp_hist_alt_bin_freq=DTYPE_1=dp_hist_alt_bin_freq       // <ID=dp_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for DP in heterozygous individuals calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100">
//#pragma GENDICT INFO_dp_hist_all_bin_freq=DTYPE_1=dp_hist_all_bin_freq     // <ID=dp_hist_all_bin_freq,Number=A,Type=String,Description="Histogram for DP calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100">
#pragma GENDICT INFO_ab_hist_alt_bin_freq=DTYPE_1=ab_hist_alt_bin_freq       // <ID=ab_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for AB in heterozygous individuals calculated on high quality genotypes; bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00">

// Structural variants (also uses INFO/END): https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/VCF%20(Variant%20Call%20Format)%20version%204.0/encoding-structural-variants/
#pragma GENDICT INFO_SVLEN=DTYPE_1=SVLEN
#pragma GENDICT INFO_SVTYPE=DTYPE_1=SVTYPE
#pragma GENDICT INFO_CIPOS=DTYPE_1=CIPOS            // <ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
#pragma GENDICT INFO_CIEND=DTYPE_1=CIEND            // <ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">

// PacBio pbsv
#pragma GENDICT INFO_SVANN=DTYPE_1=SVANN            // <ID=SVANN,Number=.,Type=String,Description="Repeat annotation of structural variant">
#pragma GENDICT INFO_MATEID=DTYPE_1=MATEID          // <ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
#pragma GENDICT INFO_MATEDIST=DTYPE_1=MATEDIST      // <ID=MATEDIST,Number=1,Type=Integer,Description="Distance to the mate breakend for mates on the same contig">
#pragma GENDICT INFO_IMPRECISE=DTYPE_1=IMPRECISE    // <ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
#pragma GENDICT INFO_SHADOWED=DTYPE_1=SHADOWED      // <ID=SHADOWED,Number=0,Type=Flag,Description="CNV overlaps with or is encapsulated by deletion">
#pragma GENDICT FORMAT_CN=DTYPE_2=CN                // <ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">

// bcftools call
#pragma GENDICT INFO_PV4=DTYPE_1=PV4                // <ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
#pragma GENDICT INFO_RPB=DTYPE_1=RPB                // <ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
#pragma GENDICT INFO_MQB=DTYPE_1=MQB                // <ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
#pragma GENDICT INFO_BQB=DTYPE_1=BQB                // <ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
#pragma GENDICT INFO_MQSB=DTYPE_1=MQSB              // <ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)”>

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
#pragma GENDICT FORMAT_FRQ=DTYPE_2=FRQ              // <ID=FRQ,Number=.,Type=Float,Description="Frequency of each alternate allele." : https://www.ncbi.nlm.nih.gov/projects/SNP/docs/dbSNP_VCF_Submission.pdf

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

//  COSMIC: https://cancer.sanger.ac.uk/cosmic/download
// #pragma GENDICT INFO_AA=DTYPE_1=AA               // (dup) <ID=AA,Number=1,Type=String,Description="Peptide annotation">
#pragma GENDICT INFO_CDS=DTYPE_1=CDS                // <ID=CDS,Number=1,Type=String,Description="CDS annotation">
#pragma GENDICT INFO_GENE=DTYPE_1=GENE              // <ID=GENE,Number=1,Type=String,Description="Gene name">
#pragma GENDICT INFO_HGVSC=DTYPE_1=HGVSC            // <ID=HGVSC,Number=1,Type=String,Description="HGVS cds syntax">
#pragma GENDICT INFO_HGVSG=DTYPE_1=HGVSG            // <ID=HGVSG,Number=1,Type=String,Description="HGVS genomic syntax">
#pragma GENDICT INFO_HGVSP=DTYPE_1=HGVSP            // <ID=HGVSP,Number=1,Type=String,Description="HGVS peptide syntax">
#pragma GENDICT INFO_LEGACY_ID=DTYPE_1=LEGACY_ID    // <ID=LEGACY_ID,Number=1,Type=String,Description="Legacy Mutation ID">
#pragma GENDICT INFO_SO_TERM=DTYPE_1=SO_TERM        // <ID=SO_TERM,Number=1,Type=String,Description="SO term for this mutation">
#pragma GENDICT INFO_STRAND=DTYPE_1=STRAND          // <ID=STRAND,Number=1,Type=String,Description="Gene strand">
#pragma GENDICT INFO_TIER=DTYPE_1=TIER              // <ID=TIER,Number=1,Type=String,Description="Indicates to which tier of the Cancer Gene Census the gene belongs (1/2)">
#pragma GENDICT INFO_TRANSCRIPT=DTYPE_1=TRANSCRIPT  // <ID=TRANSCRIPT,Number=1,Type=String,Description="Transcript accession">
#pragma GENDICT INFO_CNT=DTYPE_1=CNT                // <ID=CNT,Number=1,Type=Integer,Description="How many samples have this mutation">
#pragma GENDICT INFO_IS_CANONICAL=DTYPE_1=IS_CANONICAL // <ID=IS_CANONICAL,Number=1,Type=String,Description="The Ensembl Canonical transcript is a single, representative transcript identified at every locus. For details see: https://www.ensembl.org/info/genome/genebuild/canonical.html">
#pragma GENDICT INFO_OLD_VARIANT=DTYPE_1=OLD_VARIANT   // <ID=OLD_VARIANT,Number=.,Type=String,Description="Original chr:pos:ref:alt encoding">
#pragma GENDICT INFO_SAMPLE_COUNT=DTYPE_1=SAMPLE_COUNT // <ID=SAMPLE_COUNT,Number=1,Type=Integer,Description="How many samples have this mutation">

// CaVEMan: https://github.com/cancerit/CaVEMan
#pragma GENDICT INFO_MP=DTYPE_1=MP                  // <ID=MP,Number=1,Type=Float,Description="Sum of CaVEMan somatic genotype probabilities">
#pragma GENDICT INFO_GP=DTYPE_1=GP                  // <ID=GP,Number=1,Type=Float,Description="Sum of CaVEMan germline genotype probabilities">
#pragma GENDICT INFO_TG=DTYPE_1=TG                  // <ID=TG,Number=1,Type=String,Description="Most probable genotype as called by CaVEMan">
#pragma GENDICT INFO_TP=DTYPE_1=TP                  // <ID=TP,Number=1,Type=Float,Description="Probability of most probable genotype as called by CaVEMan">
#pragma GENDICT INFO_SG=DTYPE_1=SG                  // <ID=SG,Number=1,Type=String,Description="2nd most probable genotype as called by CaVEMan">
#pragma GENDICT INFO_SP=DTYPE_1=SP                  // <ID=SP,Number=1,Type=Float,Description="Probability of 2nd most probable genotype as called by CaVEMan">
#pragma GENDICT INFO_DS=DTYPE_1=DS                  // <ID=DS,Number=.,Type=String,Description="DBSnp ID of known SNP">
#pragma GENDICT INFO_CA=DTYPE_1=CA                  // <ID=CA,Number=0,Type=Flag,Description="Position could not be annotated to a coding region of a transcript.">
#pragma GENDICT INFO_SNP=DTYPE_1=SNP                // <ID=SNP,Number=0,Type=Flag,Description="Position matches a dbSNP entry.">
#pragma GENDICT FORMAT_AA=DTYPE_2=AA                // <ID=AA,Number=1,Type=Integer,Description="Reads presenting a A for this position">
#pragma GENDICT FORMAT_CA=DTYPE_2=CA                // <ID=CA,Number=1,Type=Integer,Description="Reads presenting a C for this position">
#pragma GENDICT FORMAT_GA=DTYPE_2=GA                // <ID=GA,Number=1,Type=Integer,Description="Reads presenting a G for this position">
#pragma GENDICT FORMAT_TA=DTYPE_2=TA                // <ID=TA,Number=1,Type=Integer,Description="Reads presenting a T for this position">
#pragma GENDICT FORMAT_PM=DTYPE_2=PM                // <ID=PM,Number=1,Type=Float,Description="Proportion of mut allele">

// pindel: https://gmt.genome.wustl.edu/packages/pindel/user-manual.html 
#pragma GENDICT INFO_PC=DTYPE_1=PC                  // <ID=PC,Number=1,Type=String,Description="Pindel call">
//#pragma GENDICT INFO_RS=DTYPE_1=RS                // (dup) <ID=RS,Number=1,Type=Integer,Description="Range start">
#pragma GENDICT INFO_RE=DTYPE_1=RE                  // <ID=RE,Number=1,Type=Integer,Description="Range end">
#pragma GENDICT INFO_LEN=DTYPE_1=LEN                // <ID=LEN,Number=1,Type=Integer,Description="Length">
#pragma GENDICT INFO_S1=DTYPE_1=S1                  // <ID=S1,Number=1,Type=Integer,Description="S1">
#pragma GENDICT INFO_S2=DTYPE_1=S2                  // <ID=S2,Number=1,Type=Float,Description="S2">
#pragma GENDICT INFO_PA=DTYPE_1=PA                  // <ID=PA,Number=1,Type=String,Description="Annotation on the positive strand">
#pragma GENDICT INFO_NA=DTYPE_1=NA                  // <ID=NA,Number=1,Type=String,Description="Annotation on the negative strand">
#pragma GENDICT INFO_REP=DTYPE_1=REP                // <ID=REP,Number=1,Type=Integer,Description="Change repeat count with range">
#pragma GENDICT INFO_PRV=DTYPE_1=PRV                // <ID=PRV,Number=1,Type=Float,Description="Fraction of reads reporting the most prevalent change">
#pragma GENDICT INFO_F017=DTYPE_1=F017              // <ID=F017,Number=0,Type=Flag,Description="Variant must not overlap with a simple repeat">
//#pragma GENDICT FORMAT_PP=DTYPE_2=PP              // (dup) <ID=PP,Number=1,Type=Integer,Description="Pindel calls on the positive strand">
#pragma GENDICT FORMAT_NP=DTYPE_2=NP                // <ID=NP,Number=1,Type=Integer,Description="Pindel calls on the negative strand">
#pragma GENDICT FORMAT_PB=DTYPE_2=PB                // <ID=PB,Number=1,Type=Integer,Description="BWA calls on the positive strand">
#pragma GENDICT FORMAT_NB=DTYPE_2=NB                // <ID=NB,Number=1,Type=Integer,Description="BWA calls on the negative strand">
#pragma GENDICT FORMAT_PD=DTYPE_2=PD                // <ID=PD,Number=1,Type=Integer,Description="BWA mapped reads on the positive strand">
#pragma GENDICT FORMAT_ND=DTYPE_2=ND                // <ID=ND,Number=1,Type=Integer,Description="BWA mapped reads on the negative strand">
#pragma GENDICT FORMAT_PR=DTYPE_2=PR                // <ID=PR,Number=1,Type=Integer,Description="Total mapped reads on the positive strand">
#pragma GENDICT FORMAT_NR=DTYPE_2=NR                // <ID=NR,Number=1,Type=Integer,Description="Total mapped reads on the negative strand">
#pragma GENDICT FORMAT_PU=DTYPE_2=PU                // <ID=PU,Number=1,Type=Integer,Description="Unique calls on the positive strand">
#pragma GENDICT FORMAT_NU=DTYPE_2=NU                // <ID=NU,Number=1,Type=Integer,Description="Unique calls on the negative strand">

// VAGrENT - Variation Annotation Generator: 
// https://github.com/cancerit/VAGrENT   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6101192/
// https://www.researchgate.net/publication/309117351_VAGrENT_Variation_Annotation_Generator/figures?lo=1
#pragma GENDICT INFO_VD=DTYPE_1=VD                  // <ID=VD,Number=1,Type=String,Description="Vagrent Default Annotation">
#pragma GENDICT INFO_VW=DTYPE_1=VW                  // <ID=VW,Number=1,Type=String,Description="Vagrent Most Deleterious Annotation">
#pragma GENDICT INFO_VDVW_ARR=DTYPE_1=VDVW_ARR
//#pragma GENDICT INFO_VC=DTYPE_1=VC                // (dup) <ID=VW,Number=1,Type=String,Description="Variant consequence based on Vagrent default annotation">
#pragma GENDICT INFO_VT=DTYPE_1=VT                  // <ID=VW,Number=1,Type=String,Description="Variant type based on Vagrent default annotation">

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

// GIAB
#pragma GENDICT FORMAT_ADALL=DTYPE_2=ADALL          // <ID=ADALL,Number=R,Type=Integer,Description="Net allele depths across all datasets">
#pragma GENDICT FORMAT_IGT=DTYPE_2=IGT              // <ID=IGT,Number=1,Type=String,Description="Original input genotype">
#pragma GENDICT FORMAT_IPS=DTYPE_2=IPS              // <ID=IPS,Number=1,Type=String,Description="Phase set for IGT">

// dbNSFP: https://gist.github.com/sahilseth/78721eada1f0007c7afd and also https://hzhou.scholar.harvard.edu/blog/dbnsfp
#pragma GENDICT INFO_Polyphen2_HDIV_score=DTYPE_1=Polyphen2_HDIV_score // Polyphen2 score based on HumDiv, i.e. hdiv_prob.
#pragma GENDICT INFO_PUniprot_aapos=DTYPE_1=Uniprot_aapos              // amino acid position as to Uniprot_acc_Polyphen2.
#pragma GENDICT INFO_VEST3_score=DTYPE_1=VEST3_score                   // VEST 3.0 score. Score ranges from 0 to 1. The larger the score the more likely the mutation may cause functional change. 
#pragma GENDICT INFO_FATHMM_score=DTYPE_1=FATHMM_score                 // FATHMM default score (weighted for human inherited-disease mutations with Disease Ontology) (FATHMMori). Scores range from -16.13 to 10.64. The smaller the score the more likely the SNP has damaging effect.
#pragma GENDICT INFO_SiPhy_29way_pi=DTYPE_1=SiPhy_29way_pi             // The estimated stationary distribution of A, C, G and T at the site, using SiPhy algorithm based on 29 mammals genomes. 

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
extern void vcf_zip_genozip_header (SectionHeaderGenozipHeaderP header);
extern void vcf_zip_init_vb (VBlockP vb);
extern void vcf_liftover_display_lift_report (void);
extern void vcf_zip_after_compress (VBlockP vb);
extern void vcf_zip_after_vbs (void);
extern void vcf_zip_set_txt_header_flags (struct FlagsTxtHeader *f);
extern void vcf_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeaderP vb_header);
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
extern bool vcf_seg_is_big (ConstVBlockP vb, DictId dict_id, DictId st_dict_id);
extern uint32_t vcf_seg_get_vb_recon_size (VBlockP vb);

// PIZ stuff
extern void vcf_piz_genozip_header (ConstSectionHeaderGenozipHeaderP header);
extern bool vcf_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header, uint32_t *txt_data_so_far_single_0_increment);
extern void vcf_piz_recon_init (VBlockP vb);
extern IS_SKIP (vcf_piz_is_skip_section);
extern CONTAINER_FILTER_FUNC (vcf_piz_filter);
extern CONTAINER_CALLBACK (vcf_piz_container_cb);
extern CONTAINER_ITEM_CALLBACK (vcf_piz_con_item_cb);

// VCF Header stuff
extern void vcf_piz_header_init (void);
extern bool vcf_inspect_txt_header (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);
extern uint32_t vcf_header_get_num_samples (void);
extern bool vcf_header_get_has_fileformat (void);
extern void vcf_piz_finalize (void);

// VBlock stuff
extern void vcf_header_finalize(void);
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
#define INFO_PREJ_NAME  "Prej" // lower case so Prej doesn't have the same first 2 chars as PRIM (to not conflict in d2d_map)
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
                      vcf_piz_special_MUX_BY_VARTYPE, vcf_piz_special_ICNT, vcf_piz_special_SPL,\
                      vcf_piz_special_MUX_BY_SAMPLE_I, vcf_piz_special_IGT, \
                      vcf_piz_special_MUX_BY_IGT_PHASE, vcf_piz_special_main_REFALT_DEL, vcf_piz_special_mutation, \
                      vcf_piz_special_SO_TERM, vcf_piz_special_MMURI, vcf_piz_special_DEMUX_GQX }

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
SPECIAL (VCF, 41, ICNT,                vcf_piz_special_ICNT);                     // added v15.0.0
SPECIAL (VCF, 42, SPL,                 vcf_piz_special_SPL);                      // added v15.0.0
SPECIAL (VCF, 43, MUX_BY_SAMPLE_I,     vcf_piz_special_MUX_BY_SAMPLE_I);          // added v15.0.0
SPECIAL (VCF, 44, IGT,                 vcf_piz_special_IGT);                      // added v15.0.0
SPECIAL (VCF, 45, MUX_BY_IGT_PHASE,    vcf_piz_special_MUX_BY_IGT_PHASE);         // added v15.0.0
SPECIAL (VCF, 46, main_REFALT_DEL,     vcf_piz_special_main_REFALT_DEL);          // added v15.0.0
SPECIAL (VCF, 47, mutation,            vcf_piz_special_mutation);                 // added v15.0.8
SPECIAL (VCF, 48, SO_TERM,             vcf_piz_special_SO_TERM);                  // added v15.0.8
SPECIAL (VCF, 49, MMURI,               vcf_piz_special_MMURI);                    // added v15.0.8
SPECIAL (VCF, 50, DEMUX_GQX,           vcf_piz_special_DEMUX_GQX);                // added v15.0.11
#define NUM_VCF_SPECIAL 51

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

#define VCF_DICT_ID_ALIASES                           \
    /*        type       alias        maps to     */  \
    { DT_VCF, ALIAS_CTX, _INFO_END,   _VCF_POS    },  \
    { DT_VCF, ALIAS_CTX, _INFO_CIEND, _INFO_CIPOS },  \

#define dict_id_is_vcf_info_sf   dict_id_is_type_1
#define dict_id_is_vcf_format_sf dict_id_is_type_2

typedef enum { VCF_COMP_MAIN, VCF_COMP_PRIM_ONLY, VCF_COMP_LUFT_ONLY } VcfComponentType;
#define VCF_COMP_NAMES { "MAIN", "PRIM", "LUFT" }
