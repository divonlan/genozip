// ------------------------------------------------------------------
//   vcf.h
//   Copyright (C) 2020-2025 Genozip Limited. Patent Pending.
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
#pragma GENDICT VCF_MATE_POS=DTYPE_FIELD=MATEPOS    // POS that appears in a BND ALT (mate's POS)
#pragma GENDICT VCF_ID=DTYPE_FIELD=ID
#pragma GENDICT VCF_REFALT=DTYPE_FIELD=REF+ALT
#pragma GENDICT VCF_MATE_CHROM=DTYPE_FIELD=MATECROM // CHROM that appears in a BND ALT (mate's CHROM)
#pragma GENDICT VCF_MATE_CHROM0=DTYPE_FIELD=M0ATECRO// define channel=0 of mux so we can use in alias
#pragma GENDICT VCF_QUAL=DTYPE_FIELD=QUAL
#pragma GENDICT VCF_FILTER=DTYPE_FIELD=FILTER
#pragma GENDICT VCF_INFO=DTYPE_FIELD=INFO
#pragma GENDICT VCF_FORMAT=DTYPE_FIELD=FORMAT
#pragma GENDICT VCF_SAMPLES=DTYPE_FIELD=SAMPLES
#pragma GENDICT VCF_SAMPLES_0=DTYPE_FIELD=S0AMPLES  // channel_i=0 ("no mate" channel) of SAMPLES multiplexing by mate
#pragma GENDICT VCF_COPY_SAMPLE=DTYPE_FIELD=COPY_SMP// used to copy same sample in previous line
#pragma GENDICT VCF_LOOKBACK=DTYPE_FIELD=LOOKBACK   // samples lookback
#pragma GENDICT VCF_EOL=DTYPE_FIELD=EOL
#pragma GENDICT VCF_TOPLEVEL=DTYPE_FIELD=TOPLEVEL   // must be called TOPLEVEL
#pragma GENDICT VCF_COORDS=DTYPE_FIELD=COORDS       // exist in files 12.0.0-15.0.41 even if not DVCF, and needed by vcf_piz_filter to reconstruct them
#pragma GENDICT VCF_oSTATUS=DTYPE_FIELD=o$TATUS     // - " -
#pragma GENDICT VCF_LINE_NUM=DTYPE_FIELD=LINE_NUM
#pragma GENDICT VCF_MATE=DTYPE_FIELD=MATE           // mate of this variant (used by svaba) 
#pragma GENDICT VCF_DEBUG_LINES=DTYPE_FIELD=DBGLINES// used by --debug-lines

// FORMAT fields
#define VCF_FIRST_OPTIONAL_DID FORMAT_AD
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
#pragma GENDICT FORMAT_GP=DTYPE_2=GP                // (Phred)         <ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification">
                                                    // (probabilities) <ID=GP,Number=3,Type=Float,Description="Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 ">
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
#pragma GENDICT FORMAT_VAF=DTYPE_2=VAF              // <ID=VAF,Number=A,Type=Float,Description="Variant allele fractions."> see: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0664-4#:~:text=The%20variant%20allele%20frequency%20(VAF,approximately%2050%25%20or%20100%25.
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
#pragma GENDICT FORMAT_AB3=DTYPE_2=AB3              // AB exceptions channel

#pragma GENDICT FORMAT_RNC=DTYPE_2=RNC              // <ID=RNC,Number=2,Type=Character,Description="Reason for No Call in GT: . = n/a, M = Missing data, P = Partial data, I = gVCF input site is non-called, D = insufficient Depth of coverage, - = unrepresentable overlapping deletion, L = Lost/unrepresentable allele (other than deletion), U = multiple Unphased variants present, O = multiple Overlapping variants present, 1 = site is Monoallelic, no assertion about presence of REF or ALT allele">

#pragma GENDICT FORMAT_FI=DTYPE_2=FI                // <ID=FI,Number=1,Type=Integer,Description="High confidence (1) or low confidence (0) based on soft filtering values">

// PBWT fields 
#pragma GENDICT FORMAT_GT_HT=DTYPE_2=@HT    
#pragma GENDICT FORMAT_PBWT_RUNS=DTYPE_2=@1BWTRUN   // PBWT runs - MUST have a did_i higher that FORMAT_GT_HT's
#pragma GENDICT FORMAT_PBWT_FGRC=DTYPE_2=@2BWTFGR   // PBWT foreground run count - MUST be right after FORMAT_PBWT_RUNS
#pragma GENDICT FORMAT_GT_HT_BIG=DTYPE_2=@3HT_BIG   // Values of Alleles greater than NUM_SMALL_ALLELES

// INFO fields
#pragma GENDICT INFO_AC=DTYPE_1=AC                  // <ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_AF=DTYPE_1=AF                  // <ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_AN=DTYPE_1=AN                  // <ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#pragma GENDICT INFO_AA=DTYPE_1=AA                  // <ID=AA,Number=1,Type=String,Description="Ancestral Allele"> - defined in the VCF specification
#pragma GENDICT INFO_BaseCounts=DTYPE_1=BaseCounts  // <ID=BaseCounts,Number=4,Type=Integer,Description="Counts of each base">
#pragma GENDICT INFO_DP=DTYPE_1=DP                  // <ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
#pragma GENDICT INFO_SF=DTYPE_1=SF                  // <ID=SF,Number=.,Type=String,Description="Source File (index to sourceFiles, f when filtered)">
#pragma GENDICT INFO_MQ=DTYPE_1=MQ                  // <ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">    
#pragma GENDICT INFO_MQ0=DTYPE_1=MQ0                // <ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads"> (VCF spec: "Number of MAPQ == 0 reads covering this record")
#pragma GENDICT INFO_NS=DTYPE_1=NS                  // <ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">

#pragma GENDICT INFO_DP4=DTYPE_1=DP4                // <ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
#pragma GENDICT INFO_DP4_RF=DTYPE_1=DrfP4                
#pragma GENDICT INFO_DP4_RR=DTYPE_1=DrrP4                
#pragma GENDICT INFO_DP4_AF=DTYPE_1=DafP4                
#pragma GENDICT INFO_DP4_AR=DTYPE_1=DarP4          

#pragma GENDICT INFO_LDAF=DTYPE_1=LDAF              // <ID=LDAF,Number=1,Type=Float,Description="MLE Allele Frequency Accounting for LD">
#pragma GENDICT INFO_AVGPOST=DTYPE_1=AVGPOST        // <ID=AVGPOST,Number=1,Type=Float,Description="Average posterior probability from MaCH/Thunder">
#pragma GENDICT INFO_RSQ=DTYPE_1=RSQ                // <ID=RSQ,Number=1,Type=Float,Description="Genotype imputation quality from MaCH/Thunder">
#pragma GENDICT INFO_ERATE=DTYPE_1=ERATE            // <ID=ERATE,Number=1,Type=Float,Description="Per-marker Mutation rate from MaCH/Thunder">
#pragma GENDICT INFO_THETA=DTYPE_1=THETA            // <ID=THETA,Number=1,Type=Float,Description="Per-marker Transition rate from MaCH/Thunder">

// SnpEFF
#pragma GENDICT INFO_ANN=DTYPE_1=ANN                // <ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">. See: https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
#pragma GENDICT INFO_ANN_Allele=DTYPE_1=@ANN_Allele // note: nice looking tag, as might be displayed in INFO/Lref

#pragma GENDICT INFO_EFF=DTYPE_1=EFF                // See: https://pcingola.github.io/SnpEff/se_inputoutput/#eff-field-vcf-output-files

// 1000 Genomes stuff 
#pragma GENDICT INFO_ID=DTYPE_1=ID                  // <ID=ID,Number=A,Type=String,Description="Variant IDs per ALT allele.">
#pragma GENDICT INFO_MAF=DTYPE_1=MAF                // <ID=MAF,Number=1,Type=Float,Description="Frequency of the second most common allele">
#pragma GENDICT INFO_HWE=DTYPE_1=HWE                // <ID=HWE,Number=A,Type=Float,Description="HWE test (PMID:15789306); 1=good, 0=bad">
#pragma GENDICT INFO_ExcHet=DTYPE_1=ExcHet          // <ID=ExcHet,Number=A,Type=Float,Description="Test excess heterozygosity; 1=good, 0=bad">
// (dup) #pragma GENDICT FORMAT_VAF=DTYPE_2=VAF     // <ID=VAF,Number=A,Type=Float,Description="The fraction of reads with alternate allele (nALT/nSumAll)">
#pragma GENDICT FORMAT_VAF1=DTYPE_2=VAF1            // <ID=VAF1,Number=1,Type=Float,Description="The fraction of reads with alternate alleles (nSumALT/nSumAll)">
#
// GATK fields 
// Added by GATK HaplotypeCaller in a gVCF: https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
#pragma GENDICT INFO_END=DTYPE_1=END                // <ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
#pragma GENDICT INFO_MLEAC=DTYPE_1=MLEAC            // <ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
#pragma GENDICT INFO_MLEAF=DTYPE_1=MLEAF            // <ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">

#pragma GENDICT INFO_VQSLOD=DTYPE_1=VQSLOD          // <ID=VQSLOD,Number=1,Type=Float,Description="Log odds of being a true variant versus being false under the trained Gaussian mixture model">
#pragma GENDICT INFO_AS_FilterStatus=DTYPE_1=AS_FilterStatus  // <ID=AS_FilterStatus,Number=A,Type=String,Description="Filter status for each allele, as assessed by ApplyRecalibration. Note that the VCF filter field will reflect the most lenient/sensitive status across all alleles.">
#pragma GENDICT INFO_AS_SB_TABLE=DTYPE_1=AS_SB_TABLE// <ID=AS_SB_TABLE,Number=1,Type=String,Description="Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.">
#pragma GENDICT INFO_AS_UNIQ_ALT_READ_COUNT=DTYPE_1=AS_UNIQ_ALT_READ_COUNT // <ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=Integer,Description="Number of reads with unique start and mate end positions for each alt at a variant site">
#pragma GENDICT INFO_CONTQ=DTYPE_1=CONTQ            // <ID=CONTQ,Number=1,Type=Float,Description="Phred-scaled qualities that alt allele are not due to contamination">
#pragma GENDICT INFO_ECNT=DTYPE_1=ECNT              // <ID=ECNT,Number=1,Type=Integer,Description="Number of events in this haplotype">
#pragma GENDICT INFO_GERMQ=DTYPE_1=GERMQ            // <ID=GERMQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not germline variants">
#pragma GENDICT INFO_MBQ=DTYPE_1=MBQ                // <ID=MBQ,Number=R,Type=Integer,Description="median base quality">
#pragma GENDICT INFO_MFRL=DTYPE_1=MFRL              // <ID=MFRL,Number=R,Type=Integer,Description="median fragment length">
#pragma GENDICT INFO_MMQ=DTYPE_1=MMQ                // <ID=MMQ,Number=R,Type=Integer,Description="median mapping quality">
#pragma GENDICT INFO_MPOS=DTYPE_1=MPOS              // <ID=MPOS,Number=A,Type=Integer,Description="median distance from end of read">
#pragma GENDICT INFO_NALOD=DTYPE_1=NALOD            // <ID=NALOD,Number=A,Type=Float,Description="Negative log 10 odds of artifact in normal with same allele fraction as tumor">
#pragma GENDICT INFO_NCount=DTYPE_1=NCount          // <ID=NCount,Number=1,Type=Integer,Description="Count of N bases in the pileup">
#pragma GENDICT INFO_NLOD=DTYPE_1=NLOD              // <ID=NLOD,Number=A,Type=Float,Description="Normal log 10 likelihood ratio of diploid het or hom alt genotypes">
#pragma GENDICT INFO_OCM=DTYPE_1=OCM                // <ID=OCM,Number=1,Type=Integer,Description="Number of alt reads whose original alignment doesn't match the current contig.">
#pragma GENDICT INFO_PON=DTYPE_1=PON                // <ID=PON,Number=0,Type=Flag,Description="site found in panel of normals">
#pragma GENDICT INFO_POPAF=DTYPE_1=POPAF            // <ID=POPAF,Number=A,Type=Float,Description="negative log 10 population allele frequencies of alt alleles">
#pragma GENDICT INFO_ROQ=DTYPE_1=ROQ                // <ID=ROQ,Number=1,Type=Float,Description="Phred-scaled qualities that alt allele are not due to read orientation artifact">
#pragma GENDICT INFO_RPA=DTYPE_1=RPA                // <ID=RPA,Number=R,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
// (dup) #pragma GENDICT INFO_RU=DTYPE_1=RU // <ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
#pragma GENDICT INFO_SEQQ=DTYPE_1=SEQQ              // <ID=SEQQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not sequencing errors">
#pragma GENDICT INFO_STR=DTYPE_1=STR                // <ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
#pragma GENDICT INFO_STRANDQ=DTYPE_1=STRANDQ        // <ID=STRANDQ,Number=1,Type=Integer,Description="Phred-scaled quality of strand bias artifact">
#pragma GENDICT INFO_STRQ=DTYPE_1=STRQ              // <ID=STRQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors">
#pragma GENDICT INFO_TLOD=DTYPE_1=TLOD              // <ID=TLOD,Number=A,Type=Float,Description="Log 10 likelihood ratio score of variant existing versus not existing">
#pragma GENDICT INFO_R2_5P_bias=DTYPE_1=R2_5P_bias  // <ID=R2_5P_bias,Number=1,Type=Float,Description="Score based on mate bias and distance from 5 prime end">

// GATK Hard-filtering germline short variants : https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
#pragma GENDICT INFO_QD=DTYPE_1=QD                  // <ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">. See: https://gatk.broadinstitute.org/hc/en-us/articles/360041414572-QualByDepth
#pragma GENDICT INFO_FS=DTYPE_1=FS                  // Strand bias estimated using Fisher's exact test : https://gatk.broadinstitute.org/hc/en-us/articles/360037592371-FisherStrand
#pragma GENDICT INFO_SOR=DTYPE_1=SOR                // <ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">. See: https://gatk.broadinstitute.org/hc/en-us/articles/360036361772-StrandOddsRatio
#pragma GENDICT INFO_MQRankSum=DTYPE_1=MQRankSum    // <ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities"> : https://gatk.broadinstitute.org/hc/en-us/articles/360037426091-MappingQualityRankSumTest
#pragma GENDICT INFO_ReadPosRankSum=DTYPE_1=ReadPosRankSum         // <ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias"> : https://gatk.broadinstitute.org/hc/en-us/articles/360036367052-ReadPosRankSumTest
#pragma GENDICT INFO_BaseQRankSum=DTYPE_1=BaseQRankSum             // <ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities"> : https://gatk.broadinstitute.org/hc/en-us/articles/360036863231-BaseQualityRankSumTest
#pragma GENDICT INFO_ClippingRankSum=DTYPE_1=ClippingRankSum       // <ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">

#pragma GENDICT INFO_HaplotypeScore=DTYPE_1=HaplotypeScore         // <ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
#pragma GENDICT INFO_InbreedingCoeff=DTYPE_1=InbreedingCoeff       // Likelihood-based test for the consanguinity among samples : https://gatk.broadinstitute.org/hc/en-us/articles/360036351032-InbreedingCoeff
#pragma GENDICT INFO_AS_InbreedingCoeff=DTYPE_1=AS_InbreedingCoeff // Allele-specific likelihood-based test for the consanguinity among samples : https://gatk.broadinstitute.org/hc/en-us/articles/360036827291-AS-InbreedingCoeff
#pragma GENDICT INFO_ExcessHet=DTYPE_1=ExcessHet    // Phred-scaled p-value for exact test of excess heterozygosity : https://gatk.broadinstitute.org/hc/en-us/articles/360037261311-ExcessHet
#pragma GENDICT INFO_RAW_MQ=DTYPE_1=RAW_MQ          // <ID=RAW_MQ,Number=1,Type=Float,Description="Raw data for RMS Mapping Quality">
#pragma GENDICT INFO_RAW_MQandDP=DTYPE_1=RAW_MQandDP// <ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">

#pragma GENDICT INFO_QUALapprox=DTYPE_1=QUALapprox  // <ID=QUALapprox,Number=1,Type=Integer,Description="Sum of PL[0] values; used to approximate the QUAL score">
#pragma GENDICT INFO_VarDP=DTYPE_1=VarDP            // <ID=VarDP,Number=1,Type=Integer,Description="Depth over variant genotypes (does not include depth of reference samples)">

#pragma GENDICT INFO_AS_QD=DTYPE_1=AS_QD            // Allele-specific call confidence normalized by depth of sample reads supporting the allele : https://gatk.broadinstitute.org/hc/en-us/articles/360036485832-AS-QualByDepth
#pragma GENDICT INFO_AS_SOR=DTYPE_1=AS_SOR          // <ID=AS_SOR,Number=A,Type=Float,Description="Allele-specific strand bias estimated by the symmetric odds ratio test">
#pragma GENDICT INFO_AS_MQ=DTYPE_1=AS_MQ            // <ID=AS_MQ,Number=A,Type=Float,Description="Allele-specific root mean square of the mapping quality of reads across all samples">
#pragma GENDICT INFO_AS_MQRankSum=DTYPE_1=AS_MQRankSum    // <ID=AS_MQRankSum,Number=A,Type=Float,Description="Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities">
#pragma GENDICT INFO_AS_FS=DTYPE_1=AS_FS            // <ID=AS_FS,Number=A,Type=Float,Description="Allele-specific phred-scaled p-value of Fisher's exact test for strand bias">
#pragma GENDICT INFO_AS_QUALapprox=DTYPE_1=AS_QUALapprox  // <ID=AS_QUALapprox,Number=A,Type=Integer,Description="Allele-specific sum of PL[0] values; used to approximate the QUAL score">
#pragma GENDICT INFO_AS_ReadPosRankSum=DTYPE_1=AS_ReadPosRankSum  // <ID=AS_ReadPosRankSum,Number=A,Type=Float,Description="Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read position bias">
#pragma GENDICT INFO_AS_VarDP=DTYPE_1=AS_VarDP     // <ID=AS_VarDP,Number=A,Type=Integer,Description="Allele-specific depth over variant genotypes (does not include depth of reference samples)">

#pragma GENDICT FORMAT_RGQ=DTYPE_2=RGQ              // <ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
#pragma GENDICT FORMAT_MIN_DP=DTYPE_2=MIN_DP        // <ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">

// DRAGEN gVCF, see: https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/Homeref_Blocks_Format_fDG.htm
#pragma GENDICT FORMAT_SPL=DTYPE_2=SPL              // <ID=SPL,Number=.,Type=Integer,Description="Normalized, Phred-scaled likelihoods for SNPs based on the reference confidence model">
#pragma GENDICT FORMAT_ICNT=DTYPE_2=ICNT            // <ID=ICNT,Number=2,Type=Integer,Description="Counts of INDEL informative reads based on the reference confidence model">

// DRAGEN CNV files
// also: FORMAT/CN
#pragma GENDICT INFO_REFLEN=DTYPE_1=REFLEN          // <ID=REFLEN,Number=1,Type=Integer,Description="Number of REF positions included in this record"> (observed in DRAGEN)
#pragma GENDICT FORMAT_PE=DTYPE_2=PE                // <ID=PE,Number=2,Type=Integer,Description="Number of improperly paired end reads at start and stop breakpoints">
#pragma GENDICT FORMAT_BC=DTYPE_2=BC                // <ID=BC,Number=1,Type=Integer,Description="Number of bins in the region">

// DRAGEN trio fields
#pragma GENDICT FORMAT_DN=DTYPE_2=DN                // <ID=DN,Number=1,Type=String,Description="Possible values are 'Inherited', 'DeNovo' or 'LowDQ'. Threshold for passing de novo call: SNPs: 0.05, INDELs: 0.02">
#pragma GENDICT FORMAT_DPL=DTYPE_2=DPL              // <ID=DPL,Number=.,Type=Integer,Description="Normalized, Phred-scaled likelihoods used for DQ calculation">
#pragma GENDICT FORMAT_DQ=DTYPE_2=DQ                // <ID=DQ,Number=1,Type=Float,Description="De novo quality">

// Illumina IsaacVariantCaller (discontinued) : https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/basespace/isaac-wgs-user-guide-15050954b.pdf
// Also: https://github.com/sequencing/isaac_variant_caller
#pragma GENDICT INFO_SNVSB=DTYPE_1=SNVSB            // SNV site strand bias
#pragma GENDICT INFO_SNVHPOL=DTYPE_1=SNVHPOL        // SNV contextual homopolymer length
#pragma GENDICT INFO_CIGAR=DTYPE_1=CIGAR            // CIGAR alignment for each alternate indel allele
#pragma GENDICT INFO_RU=DTYPE_1=RU                  // Smallest repeating sequence unit extended or contracted in the indel allele relative to the reference. RUs longer than 20 bases are not reported.
#pragma GENDICT INFO_REFREP=DTYPE_1=REFREP          // Number of times RU is repeated in reference
#pragma GENDICT INFO_IDREP=DTYPE_1=IDREP            // Number of times RU is repeated in indel allele.
#pragma GENDICT INFO_TI=DTYPE_1=TI                  // <ID=TI,Number=.,Type=String,Description="Transcript ID">
#pragma GENDICT INFO_GI=DTYPE_1=GI                  // <ID=GI,Number=.,Type=String,Description="Gene ID">
#pragma GENDICT INFO_FC=DTYPE_1=FC                  // <ID=FC,Number=.,Type=String,Description="Functional Consequence">
#pragma GENDICT INFO_RefMinor=DTYPE_1=RefMinor      // <ID=RefMinor,Number=0,Type=Flag,Description="Denotes positions where the reference base is a minor allele and is annotated as though it were a variant">
#pragma GENDICT FORMAT_GQX=DTYPE_2=GQX              // <ID=GQX,Number=1,Type=Integer,Description="Empirically calibrated variant quality score for variant sites, otherwise Minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}">
#pragma GENDICT FORMAT_DPF=DTYPE_2=DPF              // <ID=DPF,Number=1,Type=Integer,Description="Basecalls filtered from input prior to site genotyping">
#pragma GENDICT FORMAT_DPI=DTYPE_2=DPI              // <ID=DPI,Number=1,Type=Integer,Description="Read depth associated with indel, taken from the site preceding the indel.">
//#pragma GENDICT INFO_=DTYPE_1=EXON                // (flag) <ID=EXON,Number=0,Type=Flag,Description="Exon Region">
//#pragma GENDICT INFO_BLOCKAVG_min30p3a=DTYPE_1=BLOCKAVG_min30p3a // (flag) <ID=BLOCKAVG_min30p3a,Number=0,Type=Flag,Description="Non-variant site block. All sites in a block are constrained to be non-variant, have the same filter value, and have all sample values in range [x,y], y <= max(x+3,(x*1.3)). All printed site block sample values are the minimum observed in the region spanned by the block">

// Illumina starling (looks like evolved IsaacVariantCaller): https://support.illumina.com/help/BS_App_TS_Amplicon_OLH_15055858/Content/Source/Informatics/Apps/IsaacVariantCaller_appENR.htm
// (also SNVSB, SNVHPOL, CIGAR, RU, REFREP, IDREP, BLOCKAVG_min30p3a, GQX, DPF, DPI as in IsaacVariantCaller and standard tags)
#pragma GENDICT INFO_cosmic=DTYPE_1=cosmic          // <ID=cosmic,Number=.,Type=String,Description="The numeric identifier for the variant in the Catalogue of Somatic Mutations in Cancer (COSMIC) database. Format: GenotypeIndex|Significance">
#pragma GENDICT INFO_phyloP=DTYPE_1=phyloP          // <ID=phyloP,Number=A,Type=Float,Description="PhyloP conservation score. Denotes how conserved the reference sequence is between species throughout evolution">
#pragma GENDICT INFO_AF1000G=DTYPE_1=AF1000G        // <ID=AF1000G,Number=A,Type=Float,Description="The allele frequency from all populations of 1000 genomes data">
#pragma GENDICT INFO_GMAF=DTYPE_1=GMAF              // <ID=GMAF,Number=A,Type=String,Description="Global minor allele frequency (GMAF); technically, the frequency of the second most frequent allele.  Format: GlobalMinorAllele|AlleleFreqGlobalMinor">
#pragma GENDICT INFO_clinvar=DTYPE_1=clinvar        // <ID=clinvar,Number=.,Type=String,Description="Clinical significance. Format: GenotypeIndex|Significance">
#pragma GENDICT INFO_EVS=DTYPE_1=EVS                // <ID=EVS,Number=A,Type=String,Description="Allele frequency, coverage and sample count taken from the Exome Variant Server (EVS). Format: AlleleFreqEVS|EVSCoverage|EVSSamples.">
#pragma GENDICT INFO_CSQT=DTYPE_1=CSQT              // <ID=CSQT,Number=.,Type=String,Description="Consequence type as predicted by IAE. Format: GenotypeIndex|HGNC|Transcript ID|Consequence">
#pragma GENDICT INFO_CSQR=DTYPE_1=CSQR              // <ID=CSQR,Number=.,Type=String,Description="Predicted regulatory consequence type. Format: GenotypeIndex|RegulatoryID|Consequence">
#pragma GENDICT FORMAT_VF=DTYPE_2=VF                // <ID=VF,Number=1,Type=Float,Description="Variant frequency">
//#pragma GENDICT INFO_AA=DTYPE_1=AA                // (dup) <ID=AA,Number=A,Type=String,Description="The inferred allele ancestral (if determined) to the chimpanzee/human lineage.">
//#pragma GENDICT INFO_RefMinor=DTYPE_1=RefMinor    // (flag) <ID=RefMinor,Number=0,Type=Flag,Description="Denotes positions where the reference base is a minor allele and is annotated as though it were a variant">
//#pragma GENDICT INFO_Unphased=DTYPE_1=Unphased    // (flag) <ID=Unphased,Number=0,Type=Flag,Description="Indicates a record that is within the specified phasing window of another variant but could not be phased due to lack of minimum read support.">

// 10xGenomics: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf
#pragma GENDICT FORMAT_BX=DTYPE_2=BX                // <ID=BX,Number=.,Type=String,Description="Barcodes and Associated Qual-Scores Supporting Alleles">
#pragma GENDICT FORMAT_PQ=DTYPE_2=PQ                // <ID=PQ,Number=1,Type=Integer,Description="Phred QV indicating probability at this variant is incorrectly phased">
#pragma GENDICT FORMAT_JQ=DTYPE_2=JQ                // <ID=JQ,Number=1,Type=Integer,Description="Phred QV indicating probability of a phasing switch error in gap prior to this variant">

// Ensembl VEP (Variant Effect Predictor) fields: https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html https://www.ensembl.org/info/docs/tools/vep/index.html
#pragma GENDICT INFO_AGE_HISTOGRAM_HET=DTYPE_1=AGE_HISTOGRAM_HET 
#pragma GENDICT INFO_AGE_HISTOGRAM_HOM=DTYPE_1=AGE_HISTOGRAM_HOM
#pragma GENDICT INFO_MAX_AF=DTYPE_1=MAX_AF          // highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD
#pragma GENDICT INFO_NCC=DTYPE_1=NCC                // <ID=NCC,Number=1,Type=Integer,Description="Number of no-called samples">

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

// ClinVar (also in dbSNP): https://cloud.tiledb.com/arrays/details/TileDB-Inc/clinvar-annotations-data/overview
#pragma GENDICT INFO_ALLELEID=DTYPE_1=ALLELEID      // <ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
#pragma GENDICT INFO_CLNID=DTYPE_1=CLNID            // <ID=CLNID,Number=1,Type=Integer,Description="the ClinVar Variation ID">
#pragma GENDICT INFO_CLNDN=DTYPE_1=CLNDN            // <ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
#pragma GENDICT INFO_CLNHGVS=DTYPE_1=CLNHGVS        // <ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
#pragma GENDICT INFO_CLNVI=DTYPE_1=CLNVI            // <ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
#pragma GENDICT INFO_CLNORIGIN=DTYPE_1=CLNORIGIN    // <ID=CLNORIGIN,Number=.,Type=String,Description="Allele Origin. One or more of the following values may be summed: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other">
#pragma GENDICT INFO_CLNSIG=DTYPE_1=CLNSIG          // <ID=CLNSIG,Number=.,Type=String,Description="Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 8 - confers sensitivity, 9 - risk-factor, 10 - association, 11 - protective, 12 - conflict, 13 - affects, 255 - other">
#pragma GENDICT INFO_CLNDISDB=DTYPE_1=CLNDISDB      // <ID=CLNDISDB,Number=.,Type=String,Description="Variant disease database name and ID, separated by colon (:)">
#pragma GENDICT INFO_CLNREVSTAT=DTYPE_1=CLNREVSTAT  // <ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar Review Status: no_assertion - No asserition provided by submitter, no_criteria - No assertion criteria provided by submitter, single - Classified by single submitter, mult - Classified by multiple submitters, conf - Criteria provided conflicting interpretations, exp - Reviewed by expert panel, guideline - Practice guideline">
#pragma GENDICT INFO_CLNACC=DTYPE_1=CLNACC          // <ID=CLNACC,Number=.,Type=String,Description="For each allele (comma delimited), this is a pipe-delimited list of the Clinvar RCV phenotype accession.version strings associated with that allele.">
#pragma GENDICT INFO_MC=DTYPE_1=MC                  // <ID=MC,Number=.,Type=String,Description="comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence">

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
#pragma GENDICT INFO_VRS_Allele_IDs=DTYPE_1=VRS_Allele_IDs                   // <ID=VRS_Allele_IDs,Number=R,Type=String,Description="The computed identifiers for the GA4GH VRS Alleles corresponding to the values in the REF and ALT fields">
#pragma GENDICT INFO_VRS_Starts=DTYPE_1=VRS_Starts                           // <ID=VRS_Starts,Number=R,Type=Integer,Description="Interresidue coordinates used as the location starts for the GA4GH VRS Alleles corresponding to the values in the REF and ALT fields">
#pragma GENDICT INFO_VRS_Ends=DTYPE_1=VRS_Ends                               // <ID=VRS_Ends,Number=R,Type=Integer,Description="Interresidue coordinates used as the location ends for the GA4GH VRS Alleles corresponding to the values in the REF and ALT fields">
#pragma GENDICT INFO_VRS_States=DTYPE_1=VRS_States                           // <ID=VRS_States,Number=.,Type=String,Description="The literal sequence states used for the GA4GH VRS Alleles corresponding to the values in the REF and ALT fields">
#pragma GENDICT INFO_Genes=DTYPE_1=Genes            // <ID=Genes,Number=1,Type=String,Description="Genes predicted to be impacted by varian">

// Structural variants - 1000 Genome Project conventions : https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/VCF%20(Variant%20Call%20Format)%20version%204.0/encoding-structural-variants/
#pragma GENDICT INFO_SVLEN=DTYPE_1=SVLEN            // <ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles"> 
#pragma GENDICT INFO_SVTYPE=DTYPE_1=SVTYPE          // <ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant"> 
// (dup) #pragma GENDICT INFO_END=DTYPE_1=END       // <ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record"> 
#pragma GENDICT INFO_CIPOS=DTYPE_1=CIPOS            // <ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
#pragma GENDICT INFO_CIEND=DTYPE_1=CIEND            // <ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
#pragma GENDICT INFO_HOMSEQ=DTYPE_1=HOMSEQ          // <ID=HOMSEQ,Number=1,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints. Plus strand sequence displayed.">
#pragma GENDICT INFO_HOMLEN=DTYPE_1=HOMLEN          // <ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
#pragma GENDICT INFO_BKPTID=DTYPE_1=BKPTID          // <ID=BKPTID,Number=-1,Type=String,Description="ID of the assembled alternate allele in the assembly file">
#pragma GENDICT INFO_MEINFO=DTYPE_1=MEINFO          // <ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">
#pragma GENDICT INFO_METRANS=DTYPE_1=METRANS        // <ID=METRANS,Number=4,Type=String,Description="Mobile element transduction info of the form CHR,START,END,POLARITY">
#pragma GENDICT INFO_DGVID=DTYPE_1=DGVID            // <ID=DGVID,Number=1,Type=String,Description="ID of this element in Database of Genomic Variation">
#pragma GENDICT INFO_DBVARID=DTYPE_1=DBVARID        // <ID=DBVARID,Number=1,Type=String,Description="ID of this element in DBVARID of this element in DBVAR">
#pragma GENDICT INFO_DBRIPID=DTYPE_1=DBRIPID        // <ID=DBRIPID,Number=1,Type=String,Description="ID of this element in DBRIP">
#pragma GENDICT INFO_IMPRECISE=DTYPE_1=IMPRECISE    // <ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">

// Copy Number standard fields - https://samtools.github.io/hts-specs/VCFv4.3.pdf
#pragma GENDICT FORMAT_CN=DTYPE_2=CN                // <ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
#pragma GENDICT FORMAT_CICN=DTYPE_2=CICN            // <ID=CICN,Number=2,Type=Float,Description="Confidence interval around copy number">
#pragma GENDICT FORMAT_CNQ=DTYPE_2=CNQ              // <ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
#pragma GENDICT FORMAT_CNL=DTYPE_2=CNL              // <ID=CNL,Number=.,Type=Float,Description=\"Log10-scaled copy-number likelihoods\">
                                                    // alternative: "Copy number genotype likelihood for imprecise events"
#pragma GENDICT FORMAT_CNP=DTYPE_2=CNP              // <ID=CNP,Number=G,Type=Float,Description="Copy number posterior probabilities">

// PacBio pbsv
#pragma GENDICT INFO_SVANN=DTYPE_1=SVANN            // <ID=SVANN,Number=.,Type=String,Description="Repeat annotation of structural variant">
#pragma GENDICT INFO_MATEID=DTYPE_1=MATEID          // <ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
#pragma GENDICT INFO_MATEDIST=DTYPE_1=MATEDIST      // <ID=MATEDIST,Number=1,Type=Integer,Description="Distance to the mate breakend for mates on the same contig">
#pragma GENDICT INFO_SHADOWED=DTYPE_1=SHADOWED      // <ID=SHADOWED,Number=0,Type=Flag,Description="CNV overlaps with or is encapsulated by deletion">

// SvABA https://github.com/walaj/svaba
// Also dup with similar semantics: INFO/BX INFO/SVTYPE INFO/IMPRECISE INFO/MATEID
#pragma GENDICT INFO_REPSEQ=DTYPE_1=REPSEQ          // <ID=REPSEQ,Number=1,Type=String,Description="Repeat sequence near the event">
#pragma GENDICT INFO_READNAMES=DTYPE_1=READNAMES    // <ID=READNAMES,Number=.,Type=String,Description="IDs of ALT reads">
#pragma GENDICT INFO_NM=DTYPE_1=NM                  // <ID=NM,Number=1,Type=Integer,Description="Number of mismatches of this alignment fragment to reference">
#pragma GENDICT INFO_MATENM=DTYPE_1=MATENM          // <ID=MATENM,Number=1,Type=Integer,Description="Number of mismatches of partner alignment fragment to reference">
#pragma GENDICT INFO_SECONDARY=DTYPE_1=SECONDARY    // <ID=SECONDARY,Number=0,Type=Flag,Description="SV calls comes from a secondary alignment">
#pragma GENDICT INFO_MAPQ=DTYPE_1=MAPQ              // <ID=MAPQ,Number=1,Type=Integer,Description="Mapping quality (BWA-MEM) of this fragement of the contig (-1 if discordant only)">
#pragma GENDICT INFO_MATEMAPQ=DTYPE_1=MATEMAPQ      // <ID=MATEMAPQ,Number=1,Type=Integer,Description="Mapping quality of the partner fragment of the contig">
#pragma GENDICT INFO_SUBN=DTYPE_1=SUBN              // <ID=SUBN,Number=1,Type=Integer,Description="Number of secondary alignments associated with this contig fragment">
#pragma GENDICT INFO_NUMPARTS=DTYPE_1=NUMPARTS      // <ID=NUMPARTS,Number=1,Type=Integer,Description=">
#pragma GENDICT INFO_EVDNC=DTYPE_1=EVDNC            // <ID=EVDNC,Number=1,Type=String,Description="Evidence for variant. ASSMB assembly only, ASDIS assembly+discordant. DSCRD discordant only, TSI_L templated-sequence insertion (local, e.g. AB or BC of an ABC), TSI_G global (e.g. AC of ABC)">
#pragma GENDICT INFO_SCTG=DTYPE_1=SCTG              // <ID=SCTG,Number=1,Type=String,Description="Identifier for the contig assembled by svaba to make the SV call"> (or indel instead of SV)
#pragma GENDICT INFO_INSERTION=DTYPE_1=INSERTION    // <ID=INSERTION,Number=1,Type=String,Description="Sequence insertion at the breakpoint.">
#pragma GENDICT INFO_SPAN=DTYPE_1=SPAN              // <ID=SPAN,Number=1,Type=Integer,Description="Distance between the breakpoints. -1 for interchromosomal">
#pragma GENDICT INFO_DISC_MAPQ=DTYPE_1=DISC_MAPQ    // <ID=DISC_MAPQ,Number=1.,Type=Integer,Description="Mean mapping quality of discordant reads mapped here">
// (dup) #pragma GENDICT FORMAT_SR=DTYPE_2=SR       // <ID=SR,Number=1,Type=Integer,Description="Number of spanning reads for this variant">
#pragma GENDICT FORMAT_CR=DTYPE_2=CR                // <ID=CR,Number=1,Type=Integer,Description="Number of cigar-supported reads for this variant">
#pragma GENDICT FORMAT_LR=DTYPE_2=LR                // <ID=LR,Number=1,Type=Float,Description="Log-odds that this variant is AF=0 vs AF>=0.5">
#pragma GENDICT FORMAT_LO=DTYPE_2=LO                // <ID=LO,Number=1,Type=Float,Description="Log-odds that this variant is real vs artifact">
#pragma GENDICT FORMAT_SL=DTYPE_2=SL                // <ID=SL,Number=1,Type=Float,Description="Alignment-quality Scaled log-odds, where LO is LO * (MAPQ - 2*NM)/60">
// (dup) #pragma GENDICT INFO_PON=DTYPE_1=PON       // <ID=PON,Number=1,Type=Integer,Description="Number of normal samples that have this indel present">
#pragma GENDICT INFO_DBSNP=DTYPE_1=DBSNP            // <ID=DBSNP,Number=0,Type=Flag,Description="Variant found in dbSNP">
#pragma GENDICT INFO_LOD=DTYPE_1=LOD                // <ID=LOD,Number=1,Type=Float,Description="Log-odds that this variant is real vs artifact">

// DRAGEN manta : https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/MantaVCFINFOFields_fDG.htm
// (also: IMPRECISE, SVTYPE, SVLEN, END, CIPOS, CIEND, CIGAR, MATEID, HOMSEQ, HOMLEN)
#pragma GENDICT FORMAT_FT=DTYPE_2=FT                // <ID=FT,Number=1,Type=String,Description="Sample filter, 'PASS' indicates that all filters have passed for this sample">
#pragma GENDICT FORMAT_SR=DTYPE_2=SR                // <ID=SR,Number=.,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999">
// (dup) #pragma GENDICT FORMAT_PR=DTYPE_2=PR       // <ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
#pragma GENDICT INFO_EVENT=DTYPE_1=EVENT            // <ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
#pragma GENDICT INFO_SVINSLEN=DTYPE_1=SVINSLEN      // <ID=SVINSLEN,Number=.,Type=Integer,Description="Length of insertion">
#pragma GENDICT INFO_SVINSSEQ=DTYPE_1=SVINSSEQ      // <ID=SVINSSEQ,Number=.,Type=String,Description="Sequence of insertion">
#pragma GENDICT INFO_DUPSVINSLEN=DTYPE_1=DUPSVINSLEN// <ID=DUPSVINSLEN,Number=.,Type=Integer,Description="Length of inserted sequence after duplicated reference sequence">
#pragma GENDICT INFO_DUPSVINSSEQ=DTYPE_1=DUPSVINSSEQ// <ID=DUPSVINSSEQ,Number=.,Type=String,Description="Inserted sequence after duplicated reference sequence">
#pragma GENDICT INFO_DUPHOMLEN=DTYPE_1=DUPHOMLEN    // <ID=DUPHOMLEN,Number=.,Type=Integer,Description="Length of base pair identical homology at event breakpoints excluding duplicated reference sequence">
#pragma GENDICT INFO_DUPHOMSEQ=DTYPE_1=DUPHOMSEQ    // <ID=DUPHOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical homology at event breakpoints excluding duplicated reference sequence">
#pragma GENDICT INFO_BND_DEPTH=DTYPE_1=BND_DEPTH    // <ID=BND_DEPTH,Number=1,Type=Integer,Description="Read depth at local translocation breakend">
#pragma GENDICT INFO_MATE_BND_DEPTH=DTYPE_1=MATE_BND_DEPTH // <ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description="Read depth at remote translocation mate breakend">
#pragma GENDICT INFO_LEFT_SVINSSEQ=DTYPE_1=LEFT_SVINSSEQ   // <ID=LEFT_SVINSSEQ,Number=.,Type=String,Description="Known left side of insertion for an insertion of unknown length">
#pragma GENDICT INFO_RIGHT_SVINSSEQ=DTYPE_1=RIGHT_SVINSSEQ // <ID=RIGHT_SVINSSEQ,Number=.,Type=String,Description="Known right side of insertion for an insertion of unknown length">
#pragma GENDICT INFO_JUNCTION_QUAL=DTYPE_1=JUNCTION_QUAL   // <ID=JUNCTION_QUAL,Number=1,Type=Integer,Description="If the SV junction is part of an EVENT (ie. a multi-adjacency variant), this field provides the QUAL value for the adjacency in question only">

// Delly: https://github.com/dellytools/delly
// Also: INFO: CIEND, CIPOS, END, MAPQ, SVLEN, IMPRECISE, SVTYPE, HOMLEN, MAPQ. FORMAT: GT, CN, GQ, FT. 
#pragma GENDICT INFO_SOMATIC=DTYPE_1=SOMATIC        // <ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic copy-number variant.\">
#pragma GENDICT INFO_PGERM=DTYPE_1=PGERM            // <ID=PGERM,Number=1,Type=Float,Description=\"Probability of being germline.\">
#pragma GENDICT INFO_CNDIFF=DTYPE_1=CNDIFF          // <ID=CNDIFF,Number=1,Type=Float,Description=\"Absolute tumor-normal CN difference.\">
#pragma GENDICT INFO_CNSHIFT=DTYPE_1=CNSHIFT        // <ID=CNSHIFT,Number=1,Type=Float,Description=\"Estimated CN shift.\">
#pragma GENDICT INFO_CNSD=DTYPE_1=CNSD              // <ID=CNSD,Number=1,Type=Float,Description=\"Estimated CN standard deviation.\">
// (dup) #pragma GENDICT INFO_MP=DTYPE_1=MP         // <ID=MP,Number=1,Type=Float,Description=\"Mappable fraction of CNV\">
#pragma GENDICT INFO_SVMETHOD=DTYPE_1=SVMETHOD      // <ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect CNV\">
#pragma GENDICT INFO_LINKID=DTYPE_1=LINKID          // <ID=LINKID,Number=2,Type=String,Description=\"Linked paired-end IDs.\">
#pragma GENDICT INFO_REGION=DTYPE_1=REGION          // <ID=REGION,Number=3,Type=String,Description=\"Entire spanning region of the complex SV as chr, start, end.\">
#pragma GENDICT INFO_REGION1=DTYPE_1=REGION1        // <ID=REGION1,Number=3,Type=String,Description=\"Sub-Region1 of the complex SV as chr, start, end.\">
#pragma GENDICT INFO_REGION2=DTYPE_1=REGION2        // <ID=REGION1,Number=3,Type=String,Description=\"Sub-Region2 of the complex SV as chr, start, end.\">
#pragma GENDICT INFO_REGION3=DTYPE_1=REGION3        // <ID=REGION1,Number=3,Type=String,Description=\"Sub-Region3 of the complex SV as chr, start, end.\">
#pragma GENDICT INFO_CARCONC=DTYPE_1=CARCONC        // <ID=CARCONC,Number=1,Type=Float,Description=\"Carrier concordance of the linked paired-end calls.\">
#pragma GENDICT INFO_RDRATIO=DTYPE_1=RDRATIO        // <ID=RDRATIO,Number=1,Type=Float,Description=\"Read-depth ratio of tumor vs. normal.\">" or: "of carrier vs. non-carrier." or: "of SV"
#pragma GENDICT INFO_CHR2=DTYPE_1=CHR2              // <ID=CHR2,Number=1,Type=String,Description=\"Chromosome for POS2 coordinate in case of an inter-chromosomal translocation\">
#pragma GENDICT INFO_POS2=DTYPE_1=POS2              // <ID=POS2,Number=1,Type=Integer,Description=\"Genomic position for CHR2 in case of an inter-chromosomal translocation\">
#pragma GENDICT INFO_PE=DTYPE_1=PE                  // <ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">
#pragma GENDICT INFO_SRMAPQ=DTYPE_1=SRMAPQ          // <ID=SRMAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of split-reads\">
#pragma GENDICT INFO_SR=DTYPE_1=SR                  // <ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">
#pragma GENDICT INFO_SRQ=DTYPE_1=SRQ                // <ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">
#pragma GENDICT INFO_CONSENSUS=DTYPE_1=CONSENSUS    // <ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">
#pragma GENDICT INFO_CONSBP=DTYPE_1=CONSBP          // <ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">
#pragma GENDICT INFO_CE=DTYPE_1=CE                  // <ID=CE,Number=1,Type=Float,Description=\"Consensus sequence entropy\">
#pragma GENDICT INFO_CT=DTYPE_1=CT                  // <ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">
#pragma GENDICT INFO_PRECISE=DTYPE_1=PRECISE        // <ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">
#pragma GENDICT INFO_INSLEN=DTYPE_1=INSLEN          // <ID=INSLEN,Number=1,Type=Integer,Description=\"Predicted length of the insertion\">
#pragma GENDICT FORMAT_RDCN=DTYPE_2=RDCN            // <ID=RDCN,Number=1,Type=Float,Description=\"Read-depth based copy-number estimate\"> or: "... for autosomal sites"
#pragma GENDICT FORMAT_RDSD=DTYPE_2=RDSD            // <ID=RDSD,Number=1,Type=Float,Description=\"Read-depth standard deviation\">
// (dup) #pragma GENDICT FORMAT_RC=DTYPE_2=RC       // <ID=RC,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the SV\">
#pragma GENDICT FORMAT_RCL=DTYPE_2=RCL              // <ID=RCL,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the left control region\">
#pragma GENDICT FORMAT_RCR=DTYPE_2=RCR              // <ID=RCR,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the right control region\">
#pragma GENDICT FORMAT_DR=DTYPE_2=DR                // <ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">
// (dup) #pragma GENDICT FORMAT_DV=DTYPE_2=DV       // <ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">
#pragma GENDICT FORMAT_RR=DTYPE_2=RR                // <ID=,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">
#pragma GENDICT FORMAT_RV=DTYPE_2=RV                // <ID=,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">

// MELT - Mobile Element Locator Tool: https://melt.igs.umaryland.edu/manual.php
// Also: MEINFO, SVTYPE, SVLEN
#pragma GENDICT INFO_ASSESS=DTYPE_1=ASSESS          // <ID=ASSESS,Number=1,Type=Integer,Description="Provides information on evidence availible to decide insertion site.0 = No overlapping reads at site;1 = Imprecise breakpoint due to greater than expected distance between evidence;2 = discordant pair evidence only -- No split read information;3 = left side TSD evidence only;4 = right side TSD evidence only;5 = TSD decided with split reads, highest possible quality.">
#pragma GENDICT INFO_TSD=DTYPE_1=TSD                // <ID=TSD,Number=1,Type=String,Description="Precise Target Site Duplication for bases, if unknown, value will be "null"">
#pragma GENDICT INFO_INTERNAL=DTYPE_1=INTERNAL      // <ID=INTERNAL,Number=2,Type=String,Description="If insertion internal or close to a gene, listed here followed by a discriptor of the location in the gene (either INTRON, EXON_#, 5_UTR, 3_UTR, PROMOTER, or TERMINATOR)">
#pragma GENDICT INFO_DIFF=DTYPE_1=DIFF              // <ID=DIFF,Number=.,Type=String,Description="Coverage and Differences in relation to the ALU reference. Form is %2XCoverage:Differences, with differences delimited by ','">
#pragma GENDICT INFO_LP=DTYPE_1=LP                  // <ID=LP,Number=1,Type=Integer,Description="Total number of discordant pairs supporting the left side of the breakpont">
#pragma GENDICT INFO_RP=DTYPE_1=RP                  // <ID=RP,Number=1,Type=Integer,Description="Total number of discordant pairs supporting the right side of the breakpont">
#pragma GENDICT INFO_RA=DTYPE_1=RA                  // <ID=RA,Number=1,Type=Float,Description="Ratio between LP and RP, reported as log2(LP / RP)">
#pragma GENDICT INFO_PRIOR=DTYPE_1=PRIOR            // <ID=PRIOR,Number=1,Type=String,Description="True if this site was not discovered in this dataset, but was included on a provided priors list.">
// (dup) #pragma GENDICT INFO_SR=DTYPE_1=SR         // <ID=SR,Number=1,Type=Integer,Description="Total number of SRs at the estimated breakpoint for this site. Recomended to filter sites with <= 2 SRs">
#pragma GENDICT INFO_ADJLEFT=DTYPE_1=ADJLEFT        // <ID=ADJLEFT,Number=1,Type=Integer,Description=""Left bound of deletion if different from repeatmasker, else 0"">"											
#pragma GENDICT INFO_ADJRIGHT=DTYPE_1=ADJRIGHT      // <ID=ADJRIGHT,Number=1,Type=Integer,Description=""Right bound of deletion if different from repeatmasker, else 0"">"	

// Ultima Genomics
#pragma GENDICT INFO_VARIANT_TYPE=DTYPE_1=VARIANT_TYPE // <ID=VARIANT_TYPE,Number=1,Type=String,Description="Flow: type of variant: SNP/NON-H-INDEL/H-INDEL">
#pragma GENDICT INFO_SUSP_NOISY_ADJACENT_TP_VARIANT=DTYPE_1=SUSP_NOISY_ADJACENT_TP_VARIANT  // <ID=SUSP_NOISY_ADJACENT_TP_VARIANT,Number=0,Type=Flag,Description="Indicates a locus where false positive allele might be affecting a true positive allele">
#pragma GENDICT INFO_UG_HCR=DTYPE_1=UG_HCR          // <ID=UG_HCR,Number=1,Type=String,Description="Genomic Region Annotation: High confidence regions as described in Almogy et al., 2022 (gs://concordanz/hg38/UG-High-Confidence-Regions/v1.3/ug_hcr.bed)">
#pragma GENDICT INFO_XC=DTYPE_1=XC                  // <ID=XC,Number=1,Type=Integer,Description="Indicates longer hmer collapsing took place (this is a flow-based specific tag)">
#pragma GENDICT INFO_X_CSS=DTYPE_1=X_CSS            // <ID=X_CSS,Number=A,Type=String,Description="Flow: cycle skip status: cycle-skip, possible-cycle-skip, non-skip">
#pragma GENDICT INFO_X_GCC=DTYPE_1=X_GCC            // <ID=X_GCC,Number=1,Type=Float,Description="Flow: percentage of G or C in the window around hmer">
#pragma GENDICT INFO_X_HIL=DTYPE_1=X_HIL            // <ID=X_HIL,Number=A,Type=Integer,Description="Flow: length of the hmer indel, if so">
#pragma GENDICT INFO_X_HIN=DTYPE_1=X_HIN            // <ID=X_HIN,Number=A,Type=String,Description="Flow: nucleotide of the hmer indel, if so">
#pragma GENDICT INFO_X_IC=DTYPE_1=X_IC              // <ID=X_IC,Number=A,Type=String,Description="Flow: indel class: ins, del, NA">
#pragma GENDICT INFO_X_IL=DTYPE_1=X_IL              // <ID=X_IL,Number=A,Type=Integer,Description="Flow: length of indel">
#pragma GENDICT INFO_X_LM=DTYPE_1=X_LM              // <ID=X_LM,Number=A,Type=String,Description="Flow: motif to the left of the indel">
#pragma GENDICT INFO_X_RM=DTYPE_1=X_RM              // <ID=X_RM,Number=A,Type=String,Description="Flow: motif to the right of the indel">
#pragma GENDICT INFO_HPOL_RUN=DTYPE_1=HPOL_RUN      // <ID=HPOL_RUN,Number=1,Type=Flag,Description="In or close to homopolymer run">
#pragma GENDICT INFO_BLACKLST=DTYPE_1=BLACKLST      // <ID=BLACKLST,Number=.,Type=String,Description="blacklist">
#pragma GENDICT INFO_TREE_SCORE=DTYPE_1=TREE_SCORE  // <ID=TREE_SCORE,Number=1,Type=Float,Description="Filtering score">
#pragma GENDICT INFO_ASSEMBLED_HAPS=DTYPE_1=ASSEMBLED_HAPS  // <ID=ASSEMBLED_HAPS,Number=1,Type=Integer,Description="Haplotypes detected by the assembly region before haplotype filtering is applied">
#pragma GENDICT INFO_FILTERED_HAPS=DTYPE_1=FILTERED_HAPS    // <ID=FILTERED_HAPS,Number=1,Type=Integer,Description="Haplotypes filtered out by the haplotype filtering code">
#pragma GENDICT INFO_HAPDOM=DTYPE_1=HAPDOM          // <ID=HAPDOM,Number=A,Type=Float,Description="For each alt allele, fraction of read support that best fits the most-supported haplotype containing the allele">
#pragma GENDICT INFO_HAPCOMP=DTYPE_1=HAPCOMP        // <ID=HAPCOMP,Number=A,Type=Integer,Description="Edit distances of each alt allele's most common supporting haplotype from closest germline haplotype, excluding differences at the site in question.">

// population AF (always single value)
#pragma GENDICT INFO_GNOMAD_AF=DTYPE_1=GNOMAD_AF    // <ID=GNOMAD_AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed (from /cromwell_root/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz)">
#pragma GENDICT INFO_AFR_AF=DTYPE_1=AFR_AF          // <ID=AMR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AFR based on AC/AN">
#pragma GENDICT INFO_AMR_AF=DTYPE_1=AMR_AF          // <ID=AMR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AMR based on AC/AN">
#pragma GENDICT INFO_EUR_AF=DTYPE_1=EUR_AF          // <ID=AMR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from EUR based on AC/AN">
#pragma GENDICT INFO_ASN_AF=DTYPE_1=ASN_AF          // <ID=AMR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from ASN based on AC/AN">
#pragma GENDICT INFO_SAS_AF=DTYPE_1=SAS_AF          
#pragma GENDICT INFO_EAS_AF=DTYPE_1=EAS_AF          

// bcftools csq
#pragma GENDICT INFO_BCSQ=DTYPE_1=BCSQ              // <ID=BCSQ,Number=.,Type=String,Description="Local consequence annotation from BCFtools/csq, see http://samtools.github.io/bcftools/howtos/csq-calling.html for details. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change">

// bcftools call
#pragma GENDICT INFO_PV4=DTYPE_1=PV4                // <ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
#pragma GENDICT INFO_RPB=DTYPE_1=RPB                // <ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
#pragma GENDICT INFO_MQB=DTYPE_1=MQB                // <ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
#pragma GENDICT INFO_BQB=DTYPE_1=BQB                // <ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
#pragma GENDICT INFO_MQSB=DTYPE_1=MQSB              // <ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)”>

// samtools mpileup / bcftools call
#pragma GENDICT INFO_INDEL=DTYPE_1=INDEL            // <ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
#pragma GENDICT INFO_IDV=DTYPE_1=IDV                // <ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
#pragma GENDICT INFO_IMF=DTYPE_1=IMF                // <ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">
#pragma GENDICT INFO_VDB=DTYPE_1=VDB                // <ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
#pragma GENDICT INFO_RPB2=DTYPE_1=RPB2              // <ID=RPB2,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias [CDF] (bigger is better)">
#pragma GENDICT INFO_MQB2=DTYPE_1=MQB2              // <ID=MQB2,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias [CDF] (bigger is better)">
#pragma GENDICT INFO_BQB2=DTYPE_1=BQB2              // <ID=BQB2,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias [CDF] (bigger is better)">
#pragma GENDICT INFO_MQSB2=DTYPE_1=MQSB2            // <ID=MQSB2,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias [CDF] (bigger is better)">
#pragma GENDICT INFO_SGB=DTYPE_1=SGB                // <ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
#pragma GENDICT INFO_MQ0F=DTYPE_1=MQ0F              // <ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
#pragma GENDICT INFO_I16=DTYPE_1=I16                // <ID=I16,Number=16,Type=Float,Description="Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h">
#pragma GENDICT INFO_QS=DTYPE_1=QS                  // <ID=QS,Number=R,Type=Float,Description="Auxiliary tag used for calling">
#pragma GENDICT INFO_DPR=DTYPE_1=DPR                // <ID=DPR,Number=R,Type=Integer,Description="Number of high-quality bases observed for each allele">
#pragma GENDICT INFO_AD=DTYPE_1=AD                  // <ID=AD,Number=R,Type=Integer,Description="Total allelic depths">
#pragma GENDICT INFO_ADF=DTYPE_1=ADF                // <ID=ADF,Number=R,Type=Integer,Description="Total allelic depths on the forward strand">
#pragma GENDICT INFO_ADR=DTYPE_1=ADR                // <ID=ADR,Number=R,Type=Integer,Description="Total allelic depths on the reverse strand">
#pragma GENDICT FORMAT_SP=DTYPE_2=SP                // <ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
#pragma GENDICT FORMAT_DV=DTYPE_2=DV                // <ID=DV,Number=1,Type=Integer,Description="Number of high-quality non-reference bases">
#pragma GENDICT FORMAT_DPR=DTYPE_2=DPR              // <ID=DPR,Number=R,Type=Integer,Description="Number of high-quality bases observed for each allele">

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
#pragma GENDICT INFO_VT=DTYPE_1=VT                  // <ID=VT,Number=1,Type=String,Description="Variant type based on Vagrent default annotation">
                        /* callMom version */       // <ID=VT,Number=.,Type=String,Description="Alternate allele type. S=SNP, M=MNP, I=Indel"> 
                        /* 1KGP version    */       // <ID=VT,Number=1,Type=String,Description="indicates what type of variant the line represents">
                        
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
#pragma GENDICT INFO_AC_Hom=DTYPE_1=AC_Hom          // <ID=AC_Hom,Number=A,Type=Integer,Description="Allele counts in homozygous genotypes">
#pragma GENDICT INFO_AC_Het=DTYPE_1=AC_Het          // <ID=AC_Het,Number=A,Type=Integer,Description="Allele counts in heterozygous genotypes">
#pragma GENDICT INFO_AC_Hemi=DTYPE_1=AC_Hemi        // <ID=AC_Hemi,Number=A,Type=Integer,Description="Allele counts in hemizygous genotypes">
#pragma GENDICT INFO_platforms=DTYPE_1=platforms    // <ID=platforms,Number=1,Type=Integer,Description="Number of different platforms for which at least one callset called this genotype, whether filtered or not">
#pragma GENDICT INFO_datasets=DTYPE_1=datasets      // <ID=datasets,Number=1,Type=Integer,Description="Number of different datasets for which at least one callset called this genotype, whether filtered or not">
#pragma GENDICT INFO_callsets=DTYPE_1=callsets      // <ID=callsets,Number=1,Type=Integer,Description="Number of different callsets that called this genotype, whether filtered or not">
#pragma GENDICT INFO_platformnames=DTYPE_1=platformnames // <ID=platformnames,Number=.,Type=String,Description="Names of platforms for which at least one callset called this genotype, whether filtered or not">
#pragma GENDICT INFO_datasetnames=DTYPE_1=datasetnames   // <ID=datasetnames,Number=.,Type=String,Description="Names of datasets for which at least one callset called this genotype, whether filtered or not">
#pragma GENDICT INFO_callsetnames=DTYPE_1=callsetnames   // <ID=callsetnames,Number=.,Type=String,Description="Names of callsets that called this genotype, whether filtered or not">

// dbNSFP: https://gist.github.com/sahilseth/78721eada1f0007c7afd and also https://hzhou.scholar.harvard.edu/blog/dbnsfp
#pragma GENDICT INFO_Polyphen2_HDIV_score=DTYPE_1=Polyphen2_HDIV_score // Polyphen2 score based on HumDiv, i.e. hdiv_prob.
#pragma GENDICT INFO_PUniprot_aapos=DTYPE_1=Uniprot_aapos              // amino acid position as to Uniprot_acc_Polyphen2.
#pragma GENDICT INFO_VEST3_score=DTYPE_1=VEST3_score                   // VEST 3.0 score. Score ranges from 0 to 1. The larger the score the more likely the mutation may cause functional change. 
#pragma GENDICT INFO_FATHMM_score=DTYPE_1=FATHMM_score                 // FATHMM default score (weighted for human inherited-disease mutations with Disease Ontology) (FATHMMori). Scores range from -16.13 to 10.64. The smaller the score the more likely the SNP has damaging effect.
#pragma GENDICT INFO_SiPhy_29way_pi=DTYPE_1=SiPhy_29way_pi             // The estimated stationary distribution of A, C, G and T at the site, using SiPhy algorithm based on 29 mammals genomes. 

// Platypus: https://github.com/andyrimmer/Platypus
#pragma GENDICT INFO_BE=DTYPE_1=BE                  // <ID=BE,Number=.,Type=Integer,Description="End position of reference call block">
#pragma GENDICT INFO_FR=DTYPE_1=FR                  // <ID=FR,Number=.,Type=Float,Description="Estimated population frequency of variant">
#pragma GENDICT INFO_MMLQ=DTYPE_1=MMLQ              // <ID=MMLQ,Number=1,Type=Float,Description="Median minimum base quality for bases around variant">
#pragma GENDICT INFO_TC=DTYPE_1=TC                  // <ID=TC,Number=1,Type=Integer,Description="Total coverage at this locus">
#pragma GENDICT INFO_TCR=DTYPE_1=TCR                // <ID=TCR,Number=1,Type=Integer,Description="Total reverse strand coverage at this locus">
#pragma GENDICT INFO_TCF=DTYPE_1=TCF                // <ID=TCF,Number=1,Type=Integer,Description="Total forward strand coverage at this locus">
#pragma GENDICT INFO_HP=DTYPE_1=HP                  // <ID=HP,Number=1,Type=Integer,Description="Homopolymer run length around variant locus">
#pragma GENDICT INFO_WS=DTYPE_1=WS                  // <ID=WS,Number=1,Type=Integer,Description="Starting position of calling window">
#pragma GENDICT INFO_WE=DTYPE_1=WE                  // <ID=WE,Number=1,Type=Integer,Description="End position of calling window">
#pragma GENDICT INFO_Source=DTYPE_1=Source          // <ID=Source,Number=.,Type=String,Description="Was this variant suggested by Playtypus, Assembler, or from a VCF?">
#pragma GENDICT INFO_BS=DTYPE_1=BS                  // <ID=BS,Number=.,Type=Integer,Description="Start position of reference call block">
#pragma GENDICT INFO_TR=DTYPE_1=TR                  // <ID=TR,Number=.,Type=Integer,Description="Total number of reads containing this variant">
#pragma GENDICT INFO_NF=DTYPE_1=NF                  // <ID=NF,Number=.,Type=Integer,Description="Total number of forward reads containing this variant">
#pragma GENDICT INFO_NR=DTYPE_1=NR                  // <ID=NR,Number=.,Type=Integer,Description="Total number of reverse reads containing this variant">
#pragma GENDICT INFO_MGOF=DTYPE_1=MGOF              // <ID=MGOF,Number=.,Type=Integer,Description="Worst goodness-of-fit value reported across all samples">
#pragma GENDICT INFO_SbPval=DTYPE_1=SbPval          // <ID=SbPval,Number=.,Type=Float,Description="Binomial P-value for strand bias test">
#pragma GENDICT INFO_SC=DTYPE_1=SC                  // <ID=SC,Number=1,Type=String,Description="Genomic sequence 10 bases either side of variant position">
#pragma GENDICT INFO_PP=DTYPE_1=PP                  // <ID=PP,Number=.,Type=Float,Description="Posterior probability (phred scaled) that this variant segregates">
//#pragma GENDICT INFO_FS=DTYPE_1=FS                // (dup) <ID=FS,Number=.,Type=Float,Description="Fisher's exact test for strand bias (Phred scale)">
//#pragma GENDICT INFO_ReadPosRankSum=DTYPE_1=ReadPosRankSum // (dup) <ID=ReadPosRankSum,Number=.,Type=Float,Description="Mann-Whitney Rank sum test for difference between in positions of variants in reads from ref and alt">
//#pragma GENDICT INFO_MQ=DTYPE_1=MQ                // (dup) <ID=MQ,Number=.,Type=Float,Description="Root mean square of mapping qualities of reads at the variant position">
//#pragma GENDICT INFO_QD=DTYPE_1=QD                // (dup) <ID=QD,Number=1,Type=Float,Description="Variant-quality/read-depth for this variant">
#pragma GENDICT INFO_BRF=DTYPE_1=BRF                // <ID=BRF,Number=1,Type=Float,Description="Fraction of reads around this variant that failed filters">
#pragma GENDICT INFO_HapScore=DTYPE_1=HapScore      // <ID=HapScore,Number=.,Type=Integer,Description="Haplotype score measuring the number of haplotypes the variant is segregating into in a window">
#pragma GENDICT FORMAT_GOF=DTYPE_2=GOF              // <ID=GOF,Number=.,Type=Float,Description="Goodness of fit value">
#pragma GENDICT FORMAT_NV=DTYPE_2=NV                // <ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">
//#pragma GENDICT FORMAT_NR=DTYPE_2=NR              // (dup) <ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
//#pragma GENDICT FORMAT_GL=DTYPE_2=GL              // (dup) <ID=GL,Number=.,Type=Float,Description="Genotype log10-likelihoods for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites">

// GLIMPSE_phase
#pragma GENDICT FORMAT_HS=DTYPE_2=HS

// Gencove
#pragma GENDICT INFO_RAF=DTYPE_1=RAF                // <ID=RAF,Number=A,Type=Float,Description="ALT allele frequency in the reference panel">

// freebayes: https://github.com/freebayes/freebayes
// also: GT, GQ, GL, AD, DP 
#pragma GENDICT FORMAT_RO=DTYPE_2=RO                // <ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
#pragma GENDICT FORMAT_QR=DTYPE_2=QR                // <ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
#pragma GENDICT FORMAT_AO=DTYPE_2=AO                // <ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
#pragma GENDICT FORMAT_QA=DTYPE_2=QA                // <ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
#pragma GENDICT INFO_DPB=DTYPE_1=DPB                // <ID=DPB,Number=1,Type=Float,Description="Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype">

// local allele annotations (e.g. DRAGEN msVCF files)
// see: https://help.dragen.illumina.com/product-guides/dragen-v4.3/dragen-dna-pipeline/small-variant-calling/joint-analysis
// see: https://samtools.github.io/bcftools/howtos/scaling.html
#pragma GENDICT FORMAT_LAD=DTYPE_2=LAD              // <ID=LAD,Number=.,Type=Integer,Description="Localized field: Allelic Depths">
#pragma GENDICT FORMAT_LPL=DTYPE_2=LPL              // <ID=LPL,Number=.,Type=Integer,Description="Local normalized, Phred-scaled likelihoods for genotypes as in original gVCF (without allele reordering)">
#pragma GENDICT FORMAT_LAA=DTYPE_2=LAA              // <ID=LAA,Number=.,Type=String,Description="Mapping of alt allele index from original gVCF to msVCF, comma-separated, 1-based (each value is the allele index in the msVCF)">
#pragma GENDICT FORMAT_LAF=DTYPE_2=LAF              // <ID=LAF,Number=A,Type=Float,Description="Local allele fractions for alt alleles in the order listed">
#pragma GENDICT FORMAT_QL=DTYPE_2=QL                // <ID=QL,Number=1,Type=Float,Description="Phred-scaled probability that the site has no variant in this sample (original gVCF QUAL)">

// not clear which software generates these (see https://ngdc.cncb.ac.cn/gvm/download)
#pragma GENDICT INFO_MA=DTYPE_1=MA                  // <ID=MA,Number=A,Type=String,Description="Minor Allele">
// (dup) #pragma GENDICT INFO_MAF=DTYPE_1=MAF       // <ID=MAF,Number=1,Type=Float,Description="Minor Allele Frequency">

#define VCF_MAX_PLOIDY 100  // set to a reasonable 100 to avoid memory allocation explosion in case of an error in the VCF file
#if VCF_MAX_PLOIDY > 255
#error "VCF_MAX_PLOIDY cannot go beyond 255 because Ploidy is uint8_t"
#endif

#define NUM_SMALL_ALLELES 245 // REF (allele #0) + 244 ALTs (alleles # 1-244) (note: these map to ASCII 48->255,0->36 where ASCII 37-47 are reserved for special values). Important: this value is part of the file format (see codec_pbwt.c)
typedef uint8_t Allele; // elements of ht_matrix: values 48->255,0->36 for allele 0 to 244, '*' for unused, '%', '-'

// ZIP stuff
extern void vcf_zip_initialize (void);
extern void vcf_zip_finalize (bool is_last_user_txt_file);
extern void vcf_zip_genozip_header (SectionHeaderGenozipHeaderP header);
extern void vcf_zip_init_vb (VBlockP vb);
extern void vcf_zip_after_compress (VBlockP vb);
extern void vcf_zip_after_vbs (void);
extern int32_t vcf_is_header_done (bool is_eof);
extern void vcf_zip_set_txt_header_flags (struct FlagsTxtHeader *f);
extern void vcf_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeaderP vb_header);
extern bool is_vcf (STRp(header), bool *need_more);
extern bool is_bcf (STRp(header), bool *need_more);
extern char *optimize_float_3_sig_dig (VBlockP vb, ContextP ctx, STRp(snip), char *out);

// SEG stuff
extern rom vcf_zip_modify (VBlockP vb_, rom line_start, uint32_t remaining);
extern rom vcf_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void vcf_seg_initialize (VBlockP vb_);
extern void vcf_zip_after_compute (VBlockP vb);
extern void vcf_segconf_finalize (VBlockP vb);
extern void vcf_seg_finalize (VBlockP vb_);
extern bool vcf_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern bool vcf_seg_is_big (ConstVBlockP vb, DictId dict_id, DictId st_dict_id);

// PIZ stuff
extern void vcf_piz_genozip_header (ConstSectionHeaderGenozipHeaderP header);
extern bool vcf_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header);
extern void vcf_piz_vb_recon_init (VBlockP vb);
extern IS_SKIP (vcf_piz_is_skip_section);
extern CONTAINER_FILTER_FUNC (vcf_piz_filter);
extern CONTAINER_CALLBACK (vcf_piz_container_cb);
extern CONTAINER_ITEM_CALLBACK (vcf_piz_con_item_cb);

// VCF Header stuff
extern void vcf_piz_header_init (CompIType comp_i);
extern bool vcf_inspect_txt_header (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);
extern uint32_t vcf_header_get_num_samples (void);
extern uint32_t vcf_header_get_num_contigs (void);
extern uint64_t vcf_header_get_nbases (void);
extern bool vcf_header_get_has_fileformat (void);
extern void vcf_piz_finalize (bool is_last_z_file);

// VBlock stuff
extern void vcf_header_finalize(void);
extern unsigned vcf_vb_size (DataType dt);
extern unsigned vcf_vb_zip_dl_size (void);
extern void vcf_reset_line (VBlockP vb);

#define VCF_CONTIG_FMT "##contig=<ID=%.*s,length=%"PRId64">"

// Samples stuff
extern void vcf_samples_add  (rom samples_str);

// SPECIALs
SPECIAL (VCF, 0,  REFALT,              vcf_piz_special_REFALT);
SPECIAL (VCF, 1,  FORMAT,              vcf_piz_special_FORMAT)
SPECIAL (VCF, 2,  INFO_AC,             vcf_piz_special_INFO_AC);
SPECIAL (VCF, 3,  SVLEN,               vcf_piz_special_SVLEN);
SPECIAL (VCF, 4,  DS_old,              vcf_piz_special_DS_old);                   // used in files up to 12.0.42
SPECIAL (VCF, 5,  BaseCounts,          vcf_piz_special_INFO_BaseCounts);
SPECIAL (VCF, 6,  SF,                  vcf_piz_special_INFO_SF);
SPECIAL (VCF, 7,  MINUS,               piz_special_MINUS);                        // added v12.0.0 
SPECIAL (VCF, 8,  LIFT_REF,            vcf_piz_special_obsolete_dvcf);            // added v12.0.0 up to 15.0.41
SPECIAL (VCF, 9,  COPYSTAT,            vcf_piz_special_obsolete_dvcf);            // added v12.0.0 up to 15.0.41
SPECIAL (VCF, 10, other_REFALT,        vcf_piz_special_obsolete_dvcf);            // added v12.0.0 up to 15.0.41
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
SPECIAL (VCF, 22, AB,                  vcf_piz_special_AB);                       // added v13.0.0
SPECIAL (VCF, 23, GQ,                  vcf_piz_special_GQ);                       // added v13.0.0
SPECIAL (VCF, 24, MUX_BY_DOSAGExDP,    vcf_piz_special_MUX_BY_DOSAGExDP);         // added v13.0.3
SPECIAL (VCF, 25, COPY_REForALT,       vcf_piz_special_COPY_REForALT);            // added v13.0.5
SPECIAL (VCF, 26, DP_by_DP_v13,        vcf_piz_special_DP_by_DP_v13);             // added v13.0.5, removed in v14
SPECIAL (VCF, 27, PS_BY_PID,           vcf_piz_special_PS_by_PID);                // added v13.0.11
SPECIAL (VCF, 28, PGT,                 vcf_piz_special_PGT);                      // added v14.0.0
SPECIAL (VCF, 29, deferred_DP,         vcf_piz_special_deferred_DP);              // added v14.0.0 - multiple samples: INFO/DP by sum(FORMAT/DP) or sum(BaseCounts). Called DP_by_DP until 15.0.51.
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
SPECIAL (VCF, 43, MUX_BY_IS_SAMPLE_0,  vcf_piz_special_MUX_BY_IS_SAMPLE_0);       // added v15.0.0
SPECIAL (VCF, 44, IGT,                 vcf_piz_special_IGT);                      // added v15.0.0
SPECIAL (VCF, 45, MUX_BY_IGT_PHASE,    vcf_piz_special_MUX_BY_IGT_PHASE);         // added v15.0.0
SPECIAL (VCF, 46, REFALT_DEL,          vcf_piz_special_REFALT_DEL);               // added v15.0.0
SPECIAL (VCF, 47, mutation,            vcf_piz_special_mutation);                 // added v15.0.8
SPECIAL (VCF, 48, SO_TERM,             vcf_piz_special_SO_TERM);                  // added v15.0.8
SPECIAL (VCF, 49, MMURI,               vcf_piz_special_MMURI);                    // added v15.0.8
SPECIAL (VCF, 50, MUX_GQX,             vcf_piz_special_MUX_GQX);                  // added v15.0.11
SPECIAL (VCF, 51, RU,                  vcf_piz_special_RU);                       // added v15.0.13
SPECIAL (VCF, 52, IDREP,               vcf_piz_special_IDREP);                    // added v15.0.13
SPECIAL (VCF, 53, next_ALT,            vcf_piz_special_next_ALT);                 // added v15.0.25
SPECIAL (VCF, 54, MUX_BY_END,          vcf_piz_special_MUX_BY_END);               // added v15.0.26
SPECIAL (VCF, 55, MUX_BY_ISAAC_FILTER, vcf_piz_special_MUX_BY_ISAAC_FILTER);      // added v15.0.26
SPECIAL (VCF, 56, X_LM_RM,             vcf_piz_special_X_LM_RM);                  // added v15.0.28
SPECIAL (VCF, 57, X_IL,                vcf_piz_special_X_IL);                     // added v15.0.30
SPECIAL (VCF, 58, X_IC,                vcf_piz_special_X_IC);                     // added v15.0.30
SPECIAL (VCF, 59, X_HIN,               vcf_piz_special_X_HIN);                    // added v15.0.30
SPECIAL (VCF, 60, X_HIL,               vcf_piz_special_X_HIL);                    // added v15.0.30
SPECIAL (VCF, 61, VARIANT_TYPE,        vcf_piz_special_VARIANT_TYPE);             // added v15.0.30
SPECIAL (VCF, 62, PLATYPUS_SC,         vcf_piz_special_PLATYPUS_SC);              // added v15.0.30
SPECIAL (VCF, 63, PLATYPUS_HP,         vcf_piz_special_PLATYPUS_HP);              // added v15.0.30
SPECIAL (VCF, 64, INFO_MLEAF,          vcf_piz_special_INFO_MLEAF);               // added v15.0.36
SPECIAL (VCF, 65, FORMAT_AD0,          vcf_piz_special_FORMAT_AD0);               // added v15.0.37
SPECIAL (VCF, 66, MUX_FORMAT_DP,       vcf_piz_special_MUX_FORMAT_DP);            // added v15.0.37
//SPECIAL (VCF, 67, AN,                vcf_piz_special_INFO_AN);                  // In the code from v15.0.37-60, but VCF_SPECIAL_AN was never segged, so we reused its number for QR_QA
SPECIAL (VCF, 67, QR_QA,               vcf_piz_special_QR_QA);                    // added v15.0.61
SPECIAL (VCF, 68, DEFER,               vcf_piz_special_DEFER);                    // added v15.0.41
SPECIAL (VCF, 69, RPA,                 vcf_piz_special_RPA);                      // added v15.0.41
SPECIAL (VCF, 70, SVABA_MATEID,        vcf_piz_special_SVABA_MATEID);             // added v15.0.48
SPECIAL (VCF, 71, MAPQ,                vcf_piz_special_MAPQ);                     // added v15.0.48
SPECIAL (VCF, 72, SPAN,                vcf_piz_special_SPAN);                     // added v15.0.48
SPECIAL (VCF, 73, COPY_MATE,           vcf_piz_special_COPY_MATE);                // added v15.0.48
SPECIAL (VCF, 74, DEMUX_BY_MATE,       vcf_piz_special_DEMUX_BY_MATE);            // added v15.0.48
SPECIAL (VCF, 75, PBSV_MATEID,         vcf_piz_special_PBSV_MATEID);              // added v15.0.48
SPECIAL (VCF, 76, DEMUX_BY_VARTYPE,    vcf_piz_special_DEMUX_BY_VARTYPE);         // added v15.0.48
SPECIAL (VCF, 77, PBSV_ID_BND,         vcf_piz_special_PBSV_ID_BND);              // added v15.0.48
SPECIAL (VCF, 78, MANTA_CIGAR,         vcf_piz_special_manta_CIGAR);              // added v15.0.48
SPECIAL (VCF, 79, LEN_OF,              piz_special_LEN_OF);                       // added v15.0.48
SPECIAL (VCF, 80, HOMSEQ,              vcf_piz_special_HOMSEQ);                   // added v15.0.48
SPECIAL (VCF, 81, RAW_MQandDP_MQ,      vcf_piz_special_RAW_MQandDP_MQ);           // added v15.0.49
SPECIAL (VCF, 82, VT,                  vcf_piz_special_VT);                       // added v15.0.51
SPECIAL (VCF, 83, VRS_Starts,          vcf_piz_special_VRS_Starts);               // added v15.0.51
SPECIAL (VCF, 84, QUAL_BY_GP,          vcf_piz_special_QUAL_BY_GP);               // added v15.0.51
SPECIAL (VCF, 85, N_ALTS,              vcf_piz_special_N_ALTS);                   // added v15.0.51
SPECIAL (VCF, 86, N_ALLELES,           vcf_piz_special_N_ALLELES);                // added v15.0.51
SPECIAL (VCF, 87, GMAF_allele,         vcf_piz_special_GMAF_allele);              // added v15.0.57
SPECIAL (VCF, 88, PLUS,                piz_special_PLUS);                         // added v15.0.58
SPECIAL (VCF, 89, ARRAY_LEN_OF,        piz_special_ARRAY_LEN_OF);                 // added v15.0.58
SPECIAL (VCF, 90, DIVIDE_BY,           piz_special_DIVIDE_BY);                    // added v15.0.58
SPECIAL (VCF, 91, GMAF_AF,             vcf_piz_special_GMAF_AF);                  // added v15.0.58
SPECIAL (VCF, 92, COPY_SAMPLE,         vcf_piz_special_COPY_SAMPLE);              // added v15.0.69
SPECIAL (VCF, 93, LAA,                 vcf_piz_special_LAA);                      // added v15.0.69
SPECIAL (VCF, 94, MUX_BY_PREV_COPIED,  vcf_piz_special_MUX_BY_PREV_COPIED);       // added v15.0.69
SPECIAL (VCF, 95, SNVHPOL,             vcf_piz_special_SNVHPOL)                   // added v15.0.71
SPECIAL (VCF, 96, TEXTUAL_FLOAT,       piz_special_TEXTUAL_FLOAT)                 // added v15.0.71
SPECIAL (VCF, 97, DEMUX_BY_DP_CUTOFF,  vcf_piz_special_DEMUX_BY_DP_CUTOFF)        // added v15.0.71
SPECIAL (VCF, 98, DEMUX_BY_COMMON,     vcf_piz_special_DEMUX_BY_COMMON);          // added v15.0.72

#define VCF_DICT_ID_ALIASES                                                 \
    /*        type        alias                   maps to               */  \
    { DT_VCF, ALIAS_CTX,  _INFO_END,              _VCF_POS              },  \
    { DT_VCF, ALIAS_DICT, _INFO_CIEND,            _INFO_CIPOS           },  \
    { DT_VCF, ALIAS_DICT, _VCF_MATE_CHROM0,       _VCF_CHROM            },  \
    { DT_VCF, ALIAS_DICT, _INFO_MQRankSum,        _INFO_ReadPosRankSum  },  \
    { DT_VCF, ALIAS_DICT, _INFO_ClippingRankSum,  _INFO_ReadPosRankSum  },  \
    { DT_VCF, ALIAS_DICT, _INFO_BaseQRankSum,     _INFO_ReadPosRankSum  },  \
 
#define dict_id_is_vcf_info_sf   dict_id_is_type_1
#define dict_id_is_vcf_format_sf dict_id_is_type_2
