// ------------------------------------------------------------------
//   sam.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "sections.h"

#pragma GENDICT_PREFIX SAM

// -----------------------------------------------------------------------------------------------------------
// Common contexts of FASTQ and SAM - these MUST be first; in same order; same dict_id; for Deep to work.
// -----------------------------------------------------------------------------------------------------------
#pragma GENDICT SAM_RNAME=DTYPE_FIELD=RNAME // RNAME must be first

#pragma GENDICT SAM_QNAME=DTYPE_FIELD=QNAME // MAX_QNAME_ITEMS
#pragma GENDICT SAM_Q0NAME=DTYPE_1=Q0NAME // must have a did_i directly after QNAME
#pragma GENDICT SAM_Q1NAME=DTYPE_1=Q1NAME 
#pragma GENDICT SAM_Q2NAME=DTYPE_1=Q2NAME
#pragma GENDICT SAM_Q3NAME=DTYPE_1=Q3NAME
#pragma GENDICT SAM_Q4NAME=DTYPE_1=Q4NAME
#pragma GENDICT SAM_Q5NAME=DTYPE_1=Q5NAME
#pragma GENDICT SAM_Q6NAME=DTYPE_1=Q6NAME 
#pragma GENDICT SAM_Q7NAME=DTYPE_1=Q7NAME 
#pragma GENDICT SAM_Q8NAME=DTYPE_1=Q8NAME 
#pragma GENDICT SAM_Q9NAME=DTYPE_1=Q9NAME 
#pragma GENDICT SAM_QANAME=DTYPE_1=QANAME 
#pragma GENDICT SAM_QBNAME=DTYPE_1=QBNAME // if adding more Q*NAMEs - add to fastq.h and kraken.h too, and update MAX_QNAME_SEGMENTS
#pragma GENDICT SAM_QmNAME=DTYPE_1=QmNAME // QmNAME reserved for mate number (always the last dict_id in the container)

#pragma GENDICT SAM_SQBITMAP=DTYPE_FIELD=SQBITMAP
#pragma GENDICT SAM_NONREF=DTYPE_FIELD=NONREF     // these 4 fields must be in this order, right after SAM_SQBITMAP
#pragma GENDICT SAM_NONREF_X=DTYPE_FIELD=NONREF_X
#pragma GENDICT SAM_GPOS=DTYPE_FIELD=GPOS
#pragma GENDICT SAM_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT SAM_SEQMIS_A=DTYPE_FIELD=SEQMIS_A // v14: mismatch bases vs the reference, when ref=A
#pragma GENDICT SAM_SEQMIS_C=DTYPE_FIELD=SEQMIS_C // v14: mismatch bases vs the reference, when ref=A
#pragma GENDICT SAM_SEQMIS_G=DTYPE_FIELD=SEQMIS_G // v14: mismatch bases vs the reference, when ref=A
#pragma GENDICT SAM_SEQMIS_T=DTYPE_FIELD=SEQMIS_T // v14: mismatch bases vs the reference, when ref=A

#pragma GENDICT SAM_QUAL=DTYPE_FIELD=QUAL 
#pragma GENDICT SAM_DOMQRUNS=DTYPE_FIELD=DOMQRUNS // these 3 must be right after SAM_QUAL. DOMQRUNS is also used by LONGR. For backwards compatability, we can never change its name.
#pragma GENDICT SAM_QUALMPLX=DTYPE_FIELD=QUALMPLX // v14.0.0. DOMQUAL codec: dom multiplexer 
#pragma GENDICT SAM_DIVRQUAL=DTYPE_FIELD=DIVRQUAL // v14.0.0. DOMQUAL codec: lines that don't have enough dom.

#pragma GENDICT SAM_TOPLEVEL=DTYPE_FIELD=TOPLEVEL // must be called TOPLEVEL
#pragma GENDICT SAM_TAXID=DTYPE_FIELD=TAXID
#pragma GENDICT SAM_DEBUG_LINES=DTYPE_FIELD=DBGLINES      // used by --debug-lines

// contexts that exist in FASTQ and not SAM - we put them here to reserve the Did so its not 
// occupied by another SAM contexts causing dictionaries to become mingled in Deep

#pragma GENDICT UNUSED_FASTQ_DEEP_VB=DTYPE_FIELD=DEEP

#pragma GENDICT UNUSED_FASTQ_QNAME2=DTYPE_1=QNAME2     // QNAME2 is embedded in QNAME (QNAME2 items immediately follow)
#pragma GENDICT UNUSED_FASTQ_Q0NAME2=DTYPE_1=q0NAME    
#pragma GENDICT UNUSED_FASTQ_Q1NAME2=DTYPE_1=q1NAME 
#pragma GENDICT UNUSED_FASTQ_Q2NAME2=DTYPE_1=q2NAME
#pragma GENDICT UNUSED_FASTQ_Q3NAME2=DTYPE_1=q3NAME
#pragma GENDICT UNUSED_FASTQ_Q4NAME2=DTYPE_1=q4NAME
#pragma GENDICT UNUSED_FASTQ_Q5NAME2=DTYPE_1=q5NAME
#pragma GENDICT UNUSED_FASTQ_Q6NAME2=DTYPE_1=q6NAME 
#pragma GENDICT UNUSED_FASTQ_Q7NAME2=DTYPE_1=q7NAME 
#pragma GENDICT UNUSED_FASTQ_Q8NAME2=DTYPE_1=q8NAME 
#pragma GENDICT UNUSED_FASTQ_Q9NAME2=DTYPE_1=q9NAME 
#pragma GENDICT UNUSED_FASTQ_QANAME2=DTYPE_1=qANAME 

#pragma GENDICT UNUSED_FASTQ_E1L=DTYPE_FIELD=E1L
#pragma GENDICT UNUSED_FASTQ_E2L=DTYPE_FIELD=E2L

#pragma GENDICT UNUSED_FASTQ_LINE3=DTYPE_FIELD=LINE3
#pragma GENDICT UNUSED_FASTQ_T0HIRD=DTYPE_1=t0NAME    // must be directly after FASTQ_LINE3
#pragma GENDICT UNUSED_FASTQ_T1HIRD=DTYPE_1=t1NAME 
#pragma GENDICT UNUSED_FASTQ_T2HIRD=DTYPE_1=t2NAME
#pragma GENDICT UNUSED_FASTQ_T3HIRD=DTYPE_1=t3NAME
#pragma GENDICT UNUSED_FASTQ_T4HIRD=DTYPE_1=t4NAME
#pragma GENDICT UNUSED_FASTQ_T5HIRD=DTYPE_1=t5NAME
#pragma GENDICT UNUSED_FASTQ_T6HIRD=DTYPE_1=t6NAME 
#pragma GENDICT UNUSED_FASTQ_T7HIRD=DTYPE_1=t7NAME 
#pragma GENDICT UNUSED_FASTQ_T8HIRD=DTYPE_1=t8NAME 
#pragma GENDICT UNUSED_FASTQ_T9HIRD=DTYPE_1=t9NAME 
#pragma GENDICT UNUSED_FASTQ_TAHIRD=DTYPE_1=tANAME 
#pragma GENDICT UNUSED_FASTQ_TBHIRD=DTYPE_1=tBNAME 
#pragma GENDICT UNUSED_FASTQ_COPY_Q=DTYPE_1=tQcopy 

// -----------------------------------------------------------------------------------------------------------
// End of common contexts of FASTQ and SAM
// -----------------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------------
// Contexts that are unique to SAM - not in FASTQ - these must start after the highest possible FASTQ Did
// -----------------------------------------------------------------------------------------------------------

#pragma GENDICT SAM_FLAG=DTYPE_FIELD=FLAG
#pragma GENDICT SAM_FLAG0=DTYPE_FIELD=F0LAG0
#pragma GENDICT SAM_FLAG1=DTYPE_FIELD=F1LAG1
#pragma GENDICT SAM_POS=DTYPE_FIELD=POS
#pragma GENDICT SAM_MAPQ=DTYPE_FIELD=MAPQ
#pragma GENDICT SAM_CIGAR=DTYPE_FIELD=CIGAR
#pragma GENDICT SAM_RNEXT=DTYPE_FIELD=RNEXT
#pragma GENDICT SAM_PNEXT=DTYPE_FIELD=PNEXT
#pragma GENDICT SAM_P0NEXT=DTYPE_FIELD=P0NEXT // multiplexed values
#pragma GENDICT SAM_P1NEXT=DTYPE_FIELD=P1NEXT
#pragma GENDICT SAM_P2NEXT=DTYPE_FIELD=P2NEXT
#pragma GENDICT SAM_P3NEXT=DTYPE_FIELD=P3NEXT
#pragma GENDICT SAM_TLEN=DTYPE_FIELD=TLEN
#pragma GENDICT SAM_AUX=DTYPE_FIELD=AUX
#pragma GENDICT SAM_QNAMESA=DTYPE_FIELD=QNAMESA   // v14.0.0. PRIM: copy from SA Group
#pragma GENDICT SAM_QUALSA=DTYPE_FIELD=QUALSA     // v14.0.0. DEPN: Qual diff vs SA Group PRIM: copy from SA Group 
#pragma GENDICT SAM_QUAL_FLANK=DTYPE_FIELD=QUALFLNK  // v14.0.0. DEPN: Flanking regions of QUAL that are not diff'ed
#pragma GENDICT SAM_QUAL_FLANK_DOMQRUNS=DTYPE_2=Q0F_DOMQ // these 3 must be right after SAM_QUAL_FLANK. 
#pragma GENDICT SAM_QUAL_FLANK_QUALMPLX=DTYPE_2=Q1F_MPLX 
#pragma GENDICT SAM_QUAL_FLANK_DIVRQUAL=DTYPE_2=Q2F_DIVR 
#pragma GENDICT SAM_EOL=DTYPE_FIELD=EOL
#pragma GENDICT SAM_BAM_BIN=DTYPE_FIELD=BAM_BIN
#pragma GENDICT SAM_TOP2BAM=DTYPE_FIELD=TOP2BAM
#pragma GENDICT SAM_TOP2FQ=DTYPE_FIELD=TOP2FQ
#pragma GENDICT SAM_TOP2FQEX=DTYPE_FIELD=TOP2FQEX
#pragma GENDICT SAM_BUDDY=DTYPE_FIELD=BUDDY  // must be called BUDDY, expected by sam_reconstruct_from_buddy 
#pragma GENDICT SAM_SAG=DTYPE_FIELD=SAG      // PRIM and DEPN: the sag from which to copy data
#pragma GENDICT SAM_SAALN=DTYPE_FIELD=SAALN  // DEPN: sags: the alignment within sag which is this line (not needed for PRIM, as the aln_i is always 0)
#pragma GENDICT SAM_FQ_AUX=DTYPE_FIELD=MC_Z  // used for consuming some AUX fields in case of translation to FASTQ (name is "MC_Z" for as up to 14.0.25 it was called SAM_MC_Z)

// Standard AUX fields - section 1.1 here: https://samtools.github.io/hts-specs/SAMtags.pdf
#pragma GENDICT OPTION_AM_i=DTYPE_2=AM:i     // The smallest template-independent mapping quality in the template
#pragma GENDICT OPTION_AS_i=DTYPE_2=AS:i     // SAM: Alignment score generated by aligner ; STAR: the local alignment score (paired for paired-end reads) ; bowtie2: "Alignment score. Can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode)."
#pragma GENDICT OPTION_CC_Z=DTYPE_2=CC:Z     // Reference name of the next hit
#pragma GENDICT OPTION_CP_i=DTYPE_2=CP:i     // Leftmost coordinate of the next hit
#pragma GENDICT OPTION_CM_i=DTYPE_2=CM:i     // Edit distance between the color sequence and the color reference (see also NM)
#pragma GENDICT OPTION_FI_i=DTYPE_2=FI:i     // The index of segment in the template
#pragma GENDICT OPTION_H0_i=DTYPE_2=H0:i     // Number of perfect hits
#pragma GENDICT OPTION_H1_i=DTYPE_2=H1:i     // Number of 1-difference hits (see also NM)
#pragma GENDICT OPTION_H2_i=DTYPE_2=H2:i     // Number of 2-difference hits
#pragma GENDICT OPTION_MC_Z=DTYPE_2=MC:Z     // CIGAR string for mate/next segment
#pragma GENDICT OPTION_MD_Z=DTYPE_2=MD:Z     // String encoding mismatched and deleted reference bases
#pragma GENDICT OPTION_MQ_i=DTYPE_2=MQ:i     // Mapping quality of the mate/next segment
#pragma GENDICT OPTION_NH_i=DTYPE_2=NH:i     // Number of reported alignments that contain the query in the current record
#pragma GENDICT OPTION_IH_i=DTYPE_2=IH:i     // Query hit total count. Novoalign: Number of stored alignments in SAM that contains the query in the current record. Only present if there is more than one alignment reported for the read (i.e. IH <= NH)
#pragma GENDICT OPTION_HI_i=DTYPE_2=HI:i     // Query hit index (a number [1,NH])
#pragma GENDICT OPTION_NM_i=DTYPE_2=NM:i     // Edit distance to the reference
#pragma GENDICT OPTION_PQ_i=DTYPE_2=PQ:i     // Phred likelihood of the template. Novoalign: "Phred likelihood of the template, conditional on both the mapping being correct. Only for multi-template reads."
#pragma GENDICT OPTION_SM_i=DTYPE_2=SM:i     // Template-independent mapping quality
#pragma GENDICT OPTION_TC_i=DTYPE_2=TC:i     // The number of segments in the template
#pragma GENDICT OPTION_UQ_i=DTYPE_2=UQ:i     // Phred likelihood of the segment, conditional on the mapping being correct

#pragma GENDICT OPTION_U2_Z=DTYPE_2=U2:Z     // Phred probability of the 2nd call being wrong conditional on the best being wrong
#pragma GENDICT OPTION_U2_DOMQRUNS=DTYPE_2=DOMQRUNS // these 3 must be right after SAM_QUAL. DOMQRUNS is also used by LONGR. For backwards compatability, we can never change its name.
#pragma GENDICT OPTION_U2_QUALMPLX=DTYPE_2=QUALMPLX // v14.0.0. DOMQUAL alg: dom multiplexer 
#pragma GENDICT OPTION_U2_DIVRQUAL=DTYPE_2=DIVRQUAL // v14.0.0. DOMQUAL alg: lines that don't have enough dom. NORMQ codec: lines that are of length other than segconf.sam_seq_len 

#pragma GENDICT OPTION_E2_Z=DTYPE_2=E2:Z     // The 2nd most likely base calls
#pragma GENDICT OPTION_2NONREF=DTYPE_2=N2ONREF // these 4 fields must be in this order, right after OPTION_E2_Z
#pragma GENDICT OPTION_N2ONREFX=DTYPE_2=n2ONREFX
#pragma GENDICT OPTION_2GPOS=DTYPE_FIELD=G2POS
#pragma GENDICT OPTION_S2TRAND=DTYPE_2=S2TRAND

#pragma GENDICT OPTION_SA_Z=DTYPE_2=SA:Z     // Other canonical alignments in a chimeric alignment
#pragma GENDICT OPTION_SA_RNAME=DTYPE_2=S0A_RNAME 
#pragma GENDICT OPTION_SA_STRAND=DTYPE_2=S1A_STRAND 
#pragma GENDICT OPTION_SA_POS=DTYPE_2=S2A_POS 
#pragma GENDICT OPTION_SA_CIGAR=DTYPE_2=S3A_CIGAR 
#pragma GENDICT OPTION_SA_NM=DTYPE_2=S4A_NM 
#pragma GENDICT OPTION_SA_MAPQ=DTYPE_2=S5A_MAPQ 
#pragma GENDICT OPTION_SA_MAIN=DTYPE_2=S6A_MAIN 

// Standard "metadata" fields - section 1.2 here: https://samtools.github.io/hts-specs/SAMtags.pdf
#pragma GENDICT OPTION_PG_Z=DTYPE_2=PG:Z     // Program. Value matches the header PG-ID tag if @PG is present.
#pragma GENDICT OPTION_PU_Z=DTYPE_2=PU:Z     // The platform unit in which the read was sequenced. If @RG headers are present, then platformunit must match the RG-PU field of one of the headers.
#pragma GENDICT OPTION_RG_Z=DTYPE_2=RG:Z     // The read group to which the read belongs. If @RG headers are present, then readgroup must match the RG-ID field of one of the headers.
#pragma GENDICT OPTION_LB_Z=DTYPE_2=LB:Z     // The library from which the read has been sequenced. If @RG headers are present, then library must match the RG-LB field of one of the headers.

// Standard "Barcode data" fields - section 1.3 here: https://samtools.github.io/hts-specs/SAMtags.pdf
#pragma GENDICT OPTION_BC_Z=DTYPE_2=BC:Z     // Barcode sequence (Identifying the sample/library)
#pragma GENDICT OPTION_BC_ARR=DTYPE_2=BC_ARR // array items (must be one after)
#pragma GENDICT OPTION_QT_Z=DTYPE_2=QT:Z     // Phred quality of the sample barcode sequence in the BC:Z tag
#pragma GENDICT OPTION_QT_ARR=DTYPE_2=QT_ARR // array items (must be one after)
#pragma GENDICT OPTION_QT_DOMQRUNS=DTYPE_FIELD=Q0T_DOMQ // these 3 must be right after OPTION_QT_Z (similar to SAM_QUAL).
#pragma GENDICT OPTION_QT_QUALMPLX=DTYPE_FIELD=Q1T_MPLX 
#pragma GENDICT OPTION_QT_DIVRQUAL=DTYPE_FIELD=Q2T_DEVQ 

#pragma GENDICT OPTION_CR_Z=DTYPE_2=CR:Z     // Cellular barcode. The uncorrected sequence bases of the cellular barcode as reported by the sequencing machine
#pragma GENDICT OPTION_CR_Z_X=DTYPE_2=C0R_X  
#pragma GENDICT OPTION_CB_Z=DTYPE_2=CB:Z     // Cell identifier, consisting of the optionally-corrected cellular barcode sequence and an optional suffix.
#pragma GENDICT OPTION_CB_ARR=DTYPE_2=C0B_ARR 
#pragma GENDICT OPTION_CB_SUFFIX=DTYPE_2=C1B_SUFF 

#pragma GENDICT OPTION_CY_Z=DTYPE_2=CY:Z     // Phred quality of the cellular barcode sequence in the CR tag
#pragma GENDICT OPTION_CY_ARR=DTYPE_2=CY_ARR // array items (must be one after)
#pragma GENDICT OPTION_CY_DOMQRUNS=DTYPE_FIELD=C0Y_DOMQ // these 3 must be right after OPTION_CY_Z (similar to SAM_QUAL).
#pragma GENDICT OPTION_CY_QUALMPLX=DTYPE_FIELD=C1Y_MPLX 
#pragma GENDICT OPTION_CY_DIVRQUAL=DTYPE_FIELD=C2Y_DEVQ 

#pragma GENDICT OPTION_OX_Z=DTYPE_2=OX:Z     // Original unique molecular barcode bases
#pragma GENDICT OPTION_MI_Z=DTYPE_2=MI:Z     // Molecular identifier; a string that uniquely identifies the molecule from which the record was derived

// Standard "Original data" fields - section 1.4 here: https://samtools.github.io/hts-specs/SAMtags.pdf
#pragma GENDICT OPTION_OQ_Z=DTYPE_2=OQ:Z     // Original Quality - "Original base quality, usually before recalibration. Same encoding as QUAL" 
#pragma GENDICT OPTION_OQ_DOMQRUNS=DTYPE_FIELD=O0Q_DOMQ // these 3 must be right after OPTION_OQ_Z (similar to SAM_QUAL).
#pragma GENDICT OPTION_OQ_QUALMPLX=DTYPE_FIELD=O1Q_MPLX 
#pragma GENDICT OPTION_OQ_DIVRQUAL=DTYPE_FIELD=O2Q_DEVQ 

#pragma GENDICT OPTION_OA_Z=DTYPE_2=OA:Z     // Original Alignment - "The original alignment information of the record prior to realignment or unalignment by a subsequent tool" 
#pragma GENDICT OPTION_OA_RNAME=DTYPE_2=O0A_RNAME  
#pragma GENDICT OPTION_OA_STRAND=DTYPE_2=O1A_STRAND 
#pragma GENDICT OPTION_OA_POS=DTYPE_2=O2A_POS 
#pragma GENDICT OPTION_OA_CIGAR=DTYPE_2=O3A_CIGAR 
#pragma GENDICT OPTION_OA_NM=DTYPE_2=O4A_NM 
#pragma GENDICT OPTION_OA_MAPQ=DTYPE_2=O5A_MAPQ 

#pragma GENDICT OPTION_OC_Z=DTYPE_2=OC:Z     // Original CIGAR - "usually before realignment. Deprecated in favour of the more general OA."
#pragma GENDICT OPTION_OP_i=DTYPE_2=OP:i     // Original POS -   "usually before realignment. Deprecated in favour of the more general OA."
                
// bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
#pragma GENDICT OPTION_X0_i=DTYPE_2=X0:i     // Number of best hits
#pragma GENDICT OPTION_X1_i=DTYPE_2=X1:i     // Number of suboptimal hits found by BWA
#pragma GENDICT OPTION_XC_i=DTYPE_2=XC:i     // Undocumented: usually seq_len minus the final soft-clip (right if forward and left if rev-comp) 
#pragma GENDICT OPTION_XN_i=DTYPE_2=XN:i     // Number of ambiguous bases in the reference (also bowtie2, bsbolt, hisat2)
#pragma GENDICT OPTION_XM_i=DTYPE_2=XM:i     // Number of mismatches in the alignment (also bowtie2, bsbolt, hisat2, tophat)
#pragma GENDICT OPTION_XO_i=DTYPE_2=XO:i     // Number of gap opens (also bowtie2, bsbolt, hisat2, tophat)
#pragma GENDICT OPTION_XG_i=DTYPE_2=XG:i     // Number of gap extentions (also bowtie2, bsbolt, hisat2, tophat)
#pragma GENDICT OPTION_XT_A=DTYPE_2=XT:A     // U=Unique alignment R=Repeat N=Not mapped M=Mate-sw (Read is fixed due to paired end rescue)
#pragma GENDICT OPTION_XS_i=DTYPE_2=XS:i     // Suboptimal alignment score. Also: Bowtie, BSSeeker2, BSBolt, tmap, gem
#pragma GENDICT OPTION_XE_i=DTYPE_2=XE:i     // Number of supporting seeds
// type not known #pragma GENDICT OPTION_XF_?=DTYPE_2=XF:?  // Support from forward/reverse alignment  

#pragma GENDICT OPTION_XA_Z=DTYPE_2=XA:Z     // (OVERLAP WITH Ion Torrent XA:Z) Alternative hits; format: (chr,pos,CIGAR,NM;)*. Also produced by gem-mapper.
#pragma GENDICT OPTION_XA_LOOKBACK=DTYPE_2=X^A_LOOKBACK
#pragma GENDICT OPTION_XA_RNAME=DTYPE_2=X0A_RNAME
#pragma GENDICT OPTION_XA_STRAND=DTYPE_2=X1A_STRAND 
#pragma GENDICT OPTION_XA_POS=DTYPE_2=X2A_POS 
#pragma GENDICT OPTION_XA_CIGAR=DTYPE_2=X3A_CIGAR 
#pragma GENDICT OPTION_XA_NM=DTYPE_2=X4A_NM 
#pragma GENDICT OPTION_XA_STRAND_POS=DTYPE_2=X5A_STRAND_POS 

#pragma GENDICT OPTION_TS_A=DTYPE_2=TS:A     // Same as as STAR's XS:A, minimap2's ts:A - defined in the SAM spec since 2020

// bowtie2 fields - http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml ; https://github.com/BenLangmead/bowtie2/sam.h
// Standard: AS:i NM:i MD:Z
// Inherited from bwa-sw: XS:i, Xs:i YN:i Yn:i XN:i 
// Inherited from bwa: X0:i X1:i XM:i, XO:i, XG:i
#pragma GENDICT OPTION_YF_Z=DTYPE_2=YF:Z     // String indicating reason why the read was filtered out
#pragma GENDICT OPTION_YS_i=DTYPE_2=YS:i     // Alignment score for opposite mate in the paired-end alignment (also Hisat2)
#pragma GENDICT OPTION_YT_Z=DTYPE_2=YT:Z     // Value of UU indicates the read was not part of a pair. Value of CP indicates the read was part of a pair and the pair aligned concordantly. Value of DP indicates the read was part of a pair and the pair aligned discordantly. Value of UP indicates the read was part of a pair but the pair failed to aligned either concordantly or discordantly (tag overlaps)
// many more, see bowtie2's sam.h

// bowtie fields - http://bowtie-bio.sourceforge.net/manual.shtml
// #pragma GENDICT OPTION_XM_i=DTYPE_2=XM:i  // (OVERLAP WITH BWA XM:i) Number of alignments (up to cap+1)
#pragma GENDICT OPTION_XA_i=DTYPE_2=XA:i     // Aligned read belongs to stratum <N>

// Hisat2 fields - http://daehwankimlab.github.io/hisat2/manual/
// YS:i, YF:Z. YT:Z - same as bowtie2
// XN:i, XM:i, XO:i, XG:i - same as bwa, bowtie2
// XS:A - same as STAR, TopHat, same as standard TS:A, minimap2's ts:A
// many more undocumented tags - see hisat2's sam.h: YN:i, Yn:i, YP:i, YM:i, XT:i, XD:i, Xd:i, XU:i, Xu:i, YE:i, Ye:i, YL:i, Yl:i, YU:i, Yu:i, XP:B:I, YR:i, ZB:i, ZF:i, Zf:i, ZM:Z, ZI:i
#pragma GENDICT OPTION_ZS_i=DTYPE_2=ZS:i     // Alignment score for the best-scoring alignment found other than the alignment reported. Can be negative. (but code says: "Pseudo-random seed for read")
#pragma GENDICT OPTION_Zs_Z=DTYPE_2=Zs:Z     // When the alignment of a read involves SNPs that are in the index, this option is used to indicate where exactly the read involves the SNPs
#pragma GENDICT OPTION_Zs_POS=DTYPE_2=Z0s_POS
#pragma GENDICT OPTION_Zs_TYPE=DTYPE_2=Z1s_TYPE
#pragma GENDICT OPTION_Zs_RS=DTYPE_2=Z2s_RS

// Ion Torrent Base Caller tags (source: "Torrent Suite Software 5.12 Help": http://192.167.218.6/ion-docs/GUID-965C5ED4-20C8-45D5-AF07-8B0008AF74AD.html)
#pragma GENDICT OPTION_ZA_i=DTYPE_2=ZA:i     // Number of library insert bases, where the library insert is defined as the sequence after the key and barcode adapter, and before the 3' adapter. (Only present if a 3' adapter was found.)
#pragma GENDICT OPTION_ZB_i=DTYPE_2=ZB:i     // Number of overlapping adapter bases. (Only present if a 3' adapter was found.)
#pragma GENDICT OPTION_ZC_B_i=DTYPE_2=ZC:B:i // B:i A vector of the following four values (only present if a 3' adapter was found): 1. The zero-based flow during which the first base of the adapter was incorporated (same as ZG). 2. The zero-based flow corresponding to the last insert base 3. Length of the last insert homopolymer 4. Zero-based index of adapter type found.
#pragma GENDICT OPTION_ZF_i=DTYPE_2=ZF:i     // The zero-indexed flow position corresponding to the first template base after 5' trimmed region
#pragma GENDICT OPTION_ZG_i=DTYPE_2=ZG:i     // The zero-based flow during which the first base of the adapter was incorporated. (Present only if a 3' adapter was found.)
#pragma GENDICT OPTION_ZM_B_s=DTYPE_2=ZM:B:s // B:s Normalized signals, which include phasing effects. Stored as floor(256*value).
#pragma GENDICT OPTION_ZP_B_f=DTYPE_2=ZP:B:f // B:f The estimated phase parameters for the read. The values are stored in the order CF (carry forward), IE (incomplete extension), and DR (droop).
#pragma GENDICT OPTION_ZT_Z=DTYPE_2=ZT:Z     // The trimmed 5’ unique molecular tag sequence. Written only if a tag was trimmed.
//#pragma GENDICT OPTION_YT_Z=DTYPE_2=YT:Z   // (overlap) The trimmed 3’ unique molecular tag sequence. Written only if a tag was trimmed.
#pragma GENDICT OPTION_ZE_Z=DTYPE_2=ZE:Z     // The 5’ trimmed sequence removed by the extra-trim-left command. Written only if a sequence was trimmed.
#pragma GENDICT OPTION_YE_Z=DTYPE_2=YE:Z     // The 3’ trimmed sequence removed by the extra-trim-right command. Written only if a sequence was trimmed.
#pragma GENDICT OPTION_ZK_Z=DTYPE_2=ZK:Z     // The trimmed 3' portion of read group specific identifiers that can vary within a read group. Written only if a tag was trimmed.
#pragma GENDICT OPTION_YK_Z=DTYPE_2=YK:Z     // The trimmed 3' portion of read group specific identifiers that can vary within a read group. Written only if a sequence was trimmed.

// Ion Torrent TMAP tags - https://github.com/iontorrent/TMAP/blob/master/doc/tmap-book.pdf
// XS:i - same as BWA, bowtie2 etc
// #pragma GENDICT OPTION_XA_Z=DTYPE_2=XA:Z  // (overlap) The algorithm that produced this mapping and from what stage. The format is the algorithm name and the zero-based stage (separated by a dash).
// #pragma GENDICT OPTION_XM_i=DTYPE_2=XM:i  // (overlap) The target length, that is, the number of reference bases spanned by the alignment.

// STAR aligner tags. Source: https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf
// full list according to the source code Parameters.h:  NH,HI,AS,NM,MD,nM,jM,jI,RG,XS,rB,vG,vA,vW,ha,ch,MC,CR,CY,UR,UY,CB,UB,GX,GN,gx,gn,sM,sS,sQ,cN
#pragma GENDICT OPTION_nM_i=DTYPE_2=nM:i     // the number of mismatches per (paired) alignment, not to be confused with NM, which is the number of mismatches in each mate.
#pragma GENDICT OPTION_jM_B_c=DTYPE_2=jM:B:c // jM:B:c,M1,M2,... intron motifs for all junctions (i.e. N in CIGAR): 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT. If splice junctions database is used, and a junction is annotated, 20 is added to its motif value.
#pragma GENDICT OPTION_jI_B_I=DTYPE_2=jI:B:I // jI:B:I,Start1,End1,Start2,End2,... Start and End of introns for all junctions (1-based).
#pragma GENDICT OPTION_XS_A=DTYPE_2=XS:A     // Transcript strand (also hisat2), same as standard TS:A, minimap2's ts:A
#pragma GENDICT OPTION_uT_A=DTYPE_2=uT:A     // Unmapped Type: 0='no acceptable seed/windows', 1='best alignment shorter than min allowed mapped length', 2='best alignment has more mismatches than max allowed number of mismatches', 3='read maps to more loci than the max number of multimappng loci', 4='unmapped mate of a mapped paired-end read'
#pragma GENDICT OPTION_vA_i=DTYPE_2=vA:i     // Variant Allele: 1 or 2 match one of the genotype alleles, 3 - no match to genotype
#pragma GENDICT OPTION_vG_Z=DTYPE_2=vG:Z     // Variant Genomic coordinate, of the variant overlapped by the read
#pragma GENDICT OPTION_vW_i=DTYPE_2=vW:i     // WASP filtering tag: 1='passed WASP filtering', 2='multi-mapping read' 3='variant base in the read is N (non-ACGT)' 4='remapped read did not map' 5='remapped read multi-maps' 6='remapped read maps to a different locus' 7='read overlaps too many variants'

// scRNA-seq fields:
// STARsolo: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
// 10xgenomics cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
// also outputs: CB:Z, CR:Z, CY:Z standard fields
#pragma GENDICT OPTION_UR_Z=DTYPE_2=UR:Z     // (alias of RX:Z) Chromium molecular barcode sequence as reported by the sequencer.
#pragma GENDICT OPTION_UB_Z=DTYPE_2=UB:Z     // (alias of BX:Z) Chromium molecular barcode sequence that is error-corrected among other molecular barcodes with the same cellular barcode and gene alignment.
#pragma GENDICT OPTION_UY_Z=DTYPE_2=UY:Z     // (alias of QX:Z) Chromium molecular barcode read quality. Phred scores as reported by sequencer.

#pragma GENDICT OPTION_GN_Z=DTYPE_2=GN:Z     // STARsolo: Gene name for unique-gene reads. CellRanger: (;-seperated list) 
#pragma GENDICT OPTION_GX_Z=DTYPE_2=GX:Z     // STARsolo: Gene ID for for unique-gene reads. CellRanger: (;-seperated list)
#pragma GENDICT OPTION_gn_Z=DTYPE_2=gn:Z     // Gene names for unique- and multi-gene reads (;-seperated list)
#pragma GENDICT OPTION_gx_Z=DTYPE_2=gx:Z     // Gene IDs for unique- and multi-gene reads (;-seperated list)
#pragma GENDICT OPTION_sS_Z=DTYPE_2=sS:Z     // Sequence of the entire barcode (CB,UMI,adapter...)
#pragma GENDICT OPTION_sQ_Z=DTYPE_2=sQ:Z     // Quality of the entire barcode
#pragma GENDICT OPTION_sM_Z=DTYPE_2=sM:Z     // Assessment of CB and UMI
#pragma GENDICT OPTION_fx_Z=DTYPE_2=fx:Z     // CellRanger: Feature identifier matched to this Feature Barcode read. Specified in the id column of the feature reference.
#pragma GENDICT OPTION_TX_Z=DTYPE_2=TX:Z     // Transcript list
#pragma GENDICT OPTION_TX_LOOKBACK=DTYPE_2=T^X_LOOKBACK
#pragma GENDICT OPTION_TX_NEGATIVE=DTYPE_2=T1X_NEG 
#pragma GENDICT OPTION_TX_GENE=DTYPE_2=T2X_GENE
#pragma GENDICT OPTION_TX_STRAND=DTYPE_2=T3X_STRAND 
#pragma GENDICT OPTION_TX_POS=DTYPE_2=T4X_POS 
#pragma GENDICT OPTION_TX_CIGAR=DTYPE_2=T5X_CIGAR 
#pragma GENDICT OPTION_TX_SAM_POS=DTYPE_2=T6X_SPOS 

#pragma GENDICT OPTION_AN_Z=DTYPE_2=AN:Z     // Antisense transcript list
#pragma GENDICT OPTION_AN_LOOKBACK=DTYPE_2=A^N_LOOKBACK
#pragma GENDICT OPTION_AN_NEGATIVE=DTYPE_2=A1N_NEG
#pragma GENDICT OPTION_AN_GENE=DTYPE_2=A2N_GENE
#pragma GENDICT OPTION_AN_STRAND=DTYPE_2=A3N_STRAND 
#pragma GENDICT OPTION_AN_POS=DTYPE_2=A4N_POS 
#pragma GENDICT OPTION_AN_CIGAR=DTYPE_2=A5N_CIGAR 
#pragma GENDICT OPTION_AN_SAM_POS=DTYPE_2=A6N_SPOS 

#pragma GENDICT OPTION_GR_Z=DTYPE_2=GR:Z     // CellRanger
#pragma GENDICT OPTION_GR_Z_X=DTYPE_2=G0R_X  

#pragma GENDICT OPTION_GY_Z=DTYPE_2=GY:Z     // CellRanger
#pragma GENDICT OPTION_GY_Z_X=DTYPE_2=G0Y_X  

#pragma GENDICT OPTION_2R_Z=DTYPE_2=2R:Z     // Sequence - perhaps the R2?
#pragma GENDICT OPTION_2Y_Z=DTYPE_2=2Y:Z     // Quality related to the sequence 2R - perhaps the R2 QUAL?
#pragma GENDICT OPTION_2Y_DOMQRUNS=DTYPE_FIELD=20Y_DOMQ // these 3 must be right after OPTION_2Y_Z (similar to SAM_QUAL).
#pragma GENDICT OPTION_2Y_QUALMPLX=DTYPE_FIELD=21Y_MPLX 
#pragma GENDICT OPTION_2Y_DIVRQUAL=DTYPE_FIELD=22Y_DEVQ 

// cellranger has many more

// 10x Genomics longranger lariat aligner: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam

#pragma GENDICT OPTION_RX_Z=DTYPE_2=RX:Z     // Standard: Sequence bases of the (possibly corrected) unique molecular identifier
#pragma GENDICT OPTION_RX_Z_X=DTYPE_2=R0X_X  // longranger: Raw Chromium barcode sequence. This read is subject to sequencing errors.
  
#pragma GENDICT OPTION_BX_Z=DTYPE_2=BX:Z     // longranger: Chromium barcode sequence that is error-corrected and confirmed against a list of known-good barcode sequences. 

#pragma GENDICT OPTION_QX_Z=DTYPE_2=QX:Z     // longranger: Raw Chromium barcode read quality. Phred scores as reported by sequencer.
#pragma GENDICT OPTION_QX_DOMQRUNS=DTYPE_FIELD=Q0X_DOMQ // these 3 must be right after OPTION_QX_Z (similar to SAM_QUAL).
#pragma GENDICT OPTION_QX_QUALMPLX=DTYPE_FIELD=Q1X_MPLX 
#pragma GENDICT OPTION_QX_DIVRQUAL=DTYPE_FIELD=Q2X_DEVQ 

#pragma GENDICT OPTION_TR_Z=DTYPE_2=TR:Z     // longranger & cellranger: Sequence of the 7 trimmed bases following the barcode sequence at the start of R1. Can be used to reconstruct the original R1 sequence.

#pragma GENDICT OPTION_TQ_Z=DTYPE_2=TQ:Z     // longranger: Quality values of the 7 trimmed bases following the barcode sequence at the start of R1. Can be used to reconstruct the original R1 quality values.
#pragma GENDICT OPTION_TQ_DOMQRUNS=DTYPE_FIELD=T0Q_DOMQ // these 3 must be right after OPTION_TQ_Z (similar to SAM_QUAL).
#pragma GENDICT OPTION_TQ_QUALMPLX=DTYPE_FIELD=T1Q_MPLX 
#pragma GENDICT OPTION_TQ_DIVRQUAL=DTYPE_FIELD=T2Q_DEVQ 

#pragma GENDICT OPTION_PC_i=DTYPE_2=PC:i     // longranger: Phred-scaled confidence that this read was phased correctly.
#pragma GENDICT OPTION_PS_i=DTYPE_2=PS:i     // longranger: Phase set containing this read. This corresponds to the phase set (PS) field in the VCF file. The value is the position of the first SNP in the phase block.
#pragma GENDICT OPTION_HP_i=DTYPE_2=HP:i     // longranger: Haplotype of the molecule that generated the read.
#pragma GENDICT OPTION_MI_i=DTYPE_2=MI:i     // longranger: Global molecule identifier for molecule that generated this read.

#pragma GENDICT OPTION_AM_A=DTYPE_2=AM:A     // longranger: 1 if this alignment in a long molecule, 0 otherwise. Alignments in long molecules will have their MAPQ boosted above alternative alignments not in molecules.
#pragma GENDICT OPTION_XM_A=DTYPE_2=XM:A     // longranger: 1 if second best alignment is in a long molecule, 0 otherwise.
//#pragma GENDICT OPTION_XT_i=DTYPE_2=XT:i   // (overlap) longranger: Indicate if there is tandem duplication affecting this alignment. 1 if second best alignment is in the same molecule as the best alignment, 0 otherwise.
#pragma GENDICT OPTION_DM_Z=DTYPE_2=DM:Z     // longranger: Mean NM among all reads of this molecule

// cpu (ChIA-PET Utilities): https://github.com/cheehongsg/CPU/wiki
// background: https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-S12-S11

// BLASR aligner tags: Source: https://github.com/jcombs1/blasr/blob/master/common/datastructures/alignmentset/SAMAlignment.h
//#pragma GENDICT OPTION_NM_i=DTYPE_2=NM:i   // (overlap) # of subreads
//#pragma GENDICT OPTION_FI_i=DTYPE_2=FI:i   // (overlap) read alignment start position (1 based)
//#pragma GENDICT OPTION_XS_i=DTYPE_2=XS:i   // (overlap) read alignment start position without counting previous soft clips (1 based)
//#pragma GENDICT OPTION_XE_i=DTYPE_2=XE:i   // (overlap) read alignment end position without counting previous soft clips (1 based)
#pragma GENDICT OPTION_XL_i=DTYPE_2=XL:i     // aligned read length
#pragma GENDICT OPTION_XQ_i=DTYPE_2=XQ:i     // query read length
#pragma GENDICT OPTION_XT_i=DTYPE_2=XT:i     // # of continues reads, always 1 for blasr
    
// PacBio tags. Source: https://pacbiofileformats.readthedocs.io/en/10.0/BAM.html
#pragma GENDICT OPTION_cx_i=DTYPE_2=cx:i     // per-read: Subread local context Flags: enum LocalContextFlags { ADAPTER_BEFORE = 1, ADAPTER_AFTER = 2, BARCODE_BEFORE = 4, BARCODE_AFTER = 8, FORWARD_PASS = 16, REVERSE_PASS = 32 }
#pragma GENDICT OPTION_qs_i=DTYPE_2=qs:i     // per-read: 0-based start of query in the ZMW read (absent in CCS)
#pragma GENDICT OPTION_qe_i=DTYPE_2=qe:i     // per-read: 0-based end of query in the ZMW read (absent in CCS)
#pragma GENDICT OPTION_ws_i=DTYPE_2=ws:i     // per-read: Start of first base of the query (‘qs’) in approximate raw frame count since start of movie. For a CCS read, the start of the first base of the first incorporated subread.
#pragma GENDICT OPTION_we_i=DTYPE_2=we:i     // per-read: Start of last base of the query (‘qe - 1’) in approximate raw frame count since start of movie. For a CCS read, the start of the last base of the last incorporated subread.
#pragma GENDICT OPTION_zm_i=DTYPE_2=zm:i     // per-read: ZMW hole number
#pragma GENDICT OPTION_np_i=DTYPE_2=np:i     // per-read: NumPasses (1 for subreads, variable for CCS—encodes number of complete passes of the insert)
#pragma GENDICT OPTION_ec_f=DTYPE_2=ec:f     // per-read: Effective coverage for CCS reads, the average subread coverage across all windows (only present in CCS reads)
#pragma GENDICT OPTION_rq_f=DTYPE_2=rq:f     // per-read: Float in [0, 1] encoding expected accuracy
#pragma GENDICT OPTION_sn_B_f=DTYPE_2=sn:B:f // per-read: 4 floats for the average signal-to-noise ratio of A, C, G, and T (in that order) over the HQRegion
#pragma GENDICT OPTION_dq_Z=DTYPE_2=dq:Z     // per-base: DeletionQV
#pragma GENDICT OPTION_dt_Z=DTYPE_2=dt:Z     // per-base: DeletionTag
#pragma GENDICT OPTION_iq_Z=DTYPE_2=iq:Z     // per-base: InsertionQV
#pragma GENDICT OPTION_mq_Z=DTYPE_2=mq:Z     // per-base: MergeQV
#pragma GENDICT OPTION_sq_Z=DTYPE_2=sq:Z     // per-base: SubstitutionQV
#pragma GENDICT OPTION_st_Z=DTYPE_2=st:Z     // per-base: SubstitutionTag
#pragma GENDICT OPTION_ip_B_C=DTYPE_2=ip:B:C // per-base: Interpulse duration (IPD) measured in frames (raw frames or codec V1)
#pragma GENDICT OPTION_pw_B_C=DTYPE_2=pw:B:C // per-base: PulseWidth measured in frames (raw frames or codec V1)
#pragma GENDICT OPTION_fi_B_C=DTYPE_2=fi:B:C // per-base (Hi-Fi kinetic info): Forward IPD (codec V1)
#pragma GENDICT OPTION_ri_B_C=DTYPE_2=ri:B:C // per-base (Hi-Fi kinetic info): Reverse IPD (codec V1)
#pragma GENDICT OPTION_fp_B_C=DTYPE_2=fp:B:C // per-base (Hi-Fi kinetic info): Forward PulseWidth (codec V1)
#pragma GENDICT OPTION_rp_B_C=DTYPE_2=rp:B:C // per-base (Hi-Fi kinetic info): Reverse PulseWidth (codec V1)
#pragma GENDICT OPTION_fn_i=DTYPE_2=fn:i     // per-base (Hi-Fi kinetic info): Forward number of complete passes (zero or more)
#pragma GENDICT OPTION_rn_i=DTYPE_2=rn:i     // per-base (Hi-Fi kinetic info): Reverse number of complete passes (zero or more)
#pragma GENDICT OPTION_sz_A=DTYPE_2=sz_A     // scrap read: ZMW classification annotation, one of N:=Normal, C:=Control, M:=Malformed, or S:=Sentinel
#pragma GENDICT OPTION_sc_A=DTYPE_2=sc_A     // scrap read: Scrap region-type annotation, one of A:=Adapter, B:=Barcode, L:=LQRegion, or F:=Filtered

// PacBio Lima tags: https://lima.how/output/bam.html
#pragma GENDICT OPTION_bc_B_S=DTYPE_2=bc:B:S // Barcode pair indices, integer codes represent 0-based position in the FASTA file of barcodes
#pragma GENDICT OPTION_bq_i=DTYPE_2=bq:i     // Barcode score / quality, normalized between 0 and 100
#pragma GENDICT OPTION_bl_Z=DTYPE_2=bl:Z     // Barcode sequence clipped from leading end
#pragma GENDICT OPTION_bt_Z=DTYPE_2=bt:Z     // Barcode sequence clipped from trailing end
#pragma GENDICT OPTION_ql_Z=DTYPE_2=ql:Z     // Qualities of barcode bases clipped from leading end, stored as a FASTQ string
#pragma GENDICT OPTION_qt_Z=DTYPE_2=qt:Z     // Qualities of barcode bases clipped from trailing end, stored as a FASTQ string
#pragma GENDICT OPTION_bx_B_i=DTYPE_2=bx_B_i // Pair of clipped barcode sequence lengths

// samtools non-standard tags
#pragma GENDICT OPTION_ms_i=DTYPE_2=ms:i     // mate score. Produced by samtools fixmate -m. Sum of those QUAL scores of mate that are >= 15. See calc_mate_score in samtools source.

// biobambam2 tags (when updating this list here, also update biobambam_programs)
#pragma GENDICT OPTION_mc_i=DTYPE_2=mc:i     // MateCoordinate.  Produced by bamsort / bamsormadup
//#pragma GENDICT OPTION_ms_i=DTYPE_2=ms:i   // MateBaseScore.   Produced by bamsort / bamsormadup. Identical to samtools' ms:i. Sum of those QUAL scores of mate that are >= 15. See https://github.com/gt1/libmaus/tree/master/src/libmaus/bambam/BamAlignmentDecoderBase.hpp getScore 
#pragma GENDICT OPTION_cq_Z=DTYPE_2=cq:Z     // Clipped Quality. Produced by bamclipXT, see: https://gitlab.com/german.tischler/biobambam2/-/blob/master/src/programs/bamclipXT.1

// Novocraft tags: http://www.novocraft.com/documentation/novosort-2/
#pragma GENDICT OPTION_Z5_i=DTYPE_2=Z5:i     // During the input phase for unsorted or name sorted alignments, Novosort calculates the signature of each read and adds this as a SAM tag (Z5:i:) to other segments of the template. Later, during the processing of the sorted alignments, we can determine the signature of a read and, from the Z5 tag, the signature of it’s pair. Reads are then grouped according to the two signatures, strand & library and duplicates detected within a group.
#pragma GENDICT OPTION_Zq_i=DTYPE_2=Zq:i

// Novoalign tags
#pragma GENDICT OPTION_YH_Z=DTYPE_2=YH:Z     // Hard Clipped bases. 5’ & 3’ bases are separated by a ‘|’. Not present if there is no hard clipping.
#pragma GENDICT OPTION_YQ_Z=DTYPE_2=YQ:Z     // Qualities for hard clipped bases.

// ngmlr tags: https://github.com/philres/ngmlr/blob/master/src/SAMWriter.cpp
//#pragma GENDICT OPTION_XS_i=DTYPE_2=XS:i   // (overlap) Unclear what this is
//#pragma GENDICT OPTION_XE_i=DTYPE_2=XE:i   // (overlap) Unclear what this is
#pragma GENDICT OPTION_XR_i=DTYPE_2=XR:i     // Query length minus soft clips
#pragma GENDICT OPTION_QS_i=DTYPE_2=QS:i     // Query start
#pragma GENDICT OPTION_QE_i=DTYPE_2=QE:i     // Query end
#pragma GENDICT OPTION_XI_f=DTYPE_2=XI:f     // Identity
#pragma GENDICT OPTION_CV_f=DTYPE_2=CV:f     // Covered

// Devle aligner tags: https://fantom.gsc.riken.jp/5/sstar/images/9/95/Delve_User_Manual.pdf
// XA:Z same as bwa
#pragma GENDICT OPTION_XP_Z=DTYPE_2=XP:Z     // posterior probability of each nucleotide in the read to be aligned to the genome

// TopHat: http://ccb.jhu.edu/software/tophat/manual.shtml https://github.com/DaehwanKimLab/tophat
// XM:i, XO:i, XG:i - same as bwa, bowtie2, HiSat ? XF:Z XP:Z MD:Z
// XS:A - same as in STAR
// Produces standard tags: NM:i, AS:i, NH:i, HI:i CP:i, CC:Z, RG:Z  

// BS-Seeker2 tags: https://github.com/BSSeeker/BSseeker2#output-files
// Also: all bowtie2 fields
#pragma GENDICT OPTION_XO_Z=DTYPE_2=XO:Z     // Orientation, from forward/reverted
//#pragma GENDICT OPTION_XS_i=DTYPE_2=XS:i   // (overlap) 1 when read is recognized as not fully converted by bisulfite treatment, or else 0
#pragma GENDICT OPTION_XM_Z=DTYPE_2=XM:Z     // number of sites for mismatch: X=methylated CG x=un-methylated CG Y=methylated CHG y=un-methylated CHG Z=methylated CHH z=un-methylated CHH
#pragma GENDICT OPTION_XG_Z=DTYPE_2=XG:Z     // genome sequences, with 2bp extended on both ends, from 5' to 3'

// BSBolt tags: https://bsbolt.readthedocs.io/en/latest/alignment_output/
// also: all BWA tags
#pragma GENDICT OPTION_XB_Z=DTYPE_2=XB:Z     // Read bisulfite conversion position and context
#pragma GENDICT OPTION_YS_Z=DTYPE_2=YS:Z     // Mapping strand (C=Crick, W=Watson) and alignment conversion pattern (C2T or G2A)

// Bismark tags: https://github.com/FelixKrueger/Bismark/tree/master/Docs#bismark-bamsam-output-default and https://github.com/FelixKrueger/Bismark/bismark
// Dragen (compatible XR,XG,XM with Bismark): https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/MPipelineMeth_fDG.htm
//#pragma GENDICT OPTION_XM_Z=DTYPE_2=XM:Z   // (overlap) (also Dragen) methylation call string
#pragma GENDICT OPTION_XR_Z=DTYPE_2=XR:Z     // read conversion state for the alignment
//#pragma GENDICT OPTION_XR_Z=DTYPE_2=XR:Z   // (overlap) (also Dragen) read conversion state for the alignment
//#pragma GENDICT OPTION_XA_Z=DTYPE_2=XA:Z   // (overlap) the number of non-bisulfite mismatches in read-1 (note: this is a numeric field despite the :Z)
//#pragma GENDICT OPTION_XB_Z=DTYPE_2=XB:Z   // (overlap) the number of non-bisulfite mismatches in read-2 (note: this is a numeric field despite the :Z)
//#pragma GENDICT OPTION_YS_Z=DTYPE_2=YS:Z   // (overlap) strand identity, see https://github.com/FelixKrueger/Bismark/issues/455
//#pragma GENDICT OPTION_XG_Z=DTYPE_2=XG:Z   // (overlap) (also Dragen) genome conversion state for the alignment

// gem3-mapper: https://github.com/smarco/gem3-mapper/blob/master/src/io/output_sam.c
// also: XS:i - same as BWA
#pragma GENDICT OPTION_XB_A=DTYPE_2=XB:A     // Bisulfite conversion: the version of the reference to which the read is mapped (either CT or GA): http://dcc.blueprint-epigenome.eu/#/md/bs_seq_grch38
#pragma GENDICT OPTION_X4_Z=DTYPE_2=X4:Z     // Column 4 of GEM output (MAP) => MAP Counters
#pragma GENDICT OPTION_X5_Z=DTYPE_2=X5:Z     // Column 5 of GEM output (MAP) => MAPs

// gem-2-sam
// also: XA:Z, XT:A - same as BWA
#pragma GENDICT OPTION_md_Z=DTYPE_2=md:Z     // eg "md:Z:27G19>1-1>1-4>3-18"

// added by GATK's BQSR (Base Quality Score Recalibration)
#pragma GENDICT OPTION_BD_Z=DTYPE_2=BD:Z     // Deletion base quality  (not used in newer versions of GATK)
#pragma GENDICT OPTION_BI_Z=DTYPE_2=BI:Z     // Insertion base quality (not used in newer versions of GATK)
#pragma GENDICT OPTION_BD_BI=DTYPE_2=BD_BI

// GATK MergeBamAlignment: https://gatk.broadinstitute.org/hc/en-us/articles/360047217211-MergeBamAlignment-Picard-
// #pragma GENDICT OPTION_XB_Z=DTYPE_2=XB:Z  // (overlap) Hard cliped bases
#pragma GENDICT OPTION_XQ_Z=DTYPE_2=XQ:Z     // Hard cliped base qualities

// minimap2 tags: https://lh3.github.io/minimap2/minimap2.html#10
#pragma GENDICT OPTION_tp_A=DTYPE_2=tp:A     // Type of aln: P/primary, S/secondary and I,i/inversion
#pragma GENDICT OPTION_cm_i=DTYPE_2=cm:i     // Number of minimizers on the chain
#pragma GENDICT OPTION_s1_i=DTYPE_2=s1:i     // Chaining score
#pragma GENDICT OPTION_s2_i=DTYPE_2=s2:i     // Chaining score of the best secondary chain
//#pragma GENDICT OPTION_ms_i=DTYPE_2=ms:i   // Conflict with biobambam and samtools: DP score of the max scoring segment in the alignment
#pragma GENDICT OPTION_nn_i=DTYPE_2=nn:i     // Number of ambiguous bases in the alignment
#pragma GENDICT OPTION_ts_A=DTYPE_2=ts:A     // Transcript strand (splice mode only)
#pragma GENDICT OPTION_cs_Z=DTYPE_2=cs:Z     // Difference string, see: https://github.com/lh3/minimap2#cs and https://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag
#pragma GENDICT OPTION_dv_f=DTYPE_2=dv:f     // Approximate per-base sequence divergence
#pragma GENDICT OPTION_de_f=DTYPE_2=de:f     // Gap-compressed per-base sequence divergence
#pragma GENDICT OPTION_rl_i=DTYPE_2=rl:i     // Length of query regions harboring repetitive seeds

#pragma GENDICT OPTION_tx_i=DTYPE_2=tx:i     // Genozip tag for taxonomy ID

// backward compatability for decompressing files compressed with older versions that had aliases to these destinations (used by ctx_initialize_predefined_ctxs)
#pragma GENDICT OPTION_CIGAR=DTYPE_2=@CIGAR  // For files compressed with 12.0.37 or older which had aliases MC:Z, OC:Z -> @CIGAR
#pragma GENDICT SAM_E2_Z=DTYPE_FIELD=E2:Z    // This used to be the destination alias from OPTION_E2_Z 
#pragma GENDICT SAM_U2_Z=DTYPE_FIELD=U2:Z    // This used to be the destination alias from OPTION_U2_Z 

// Note: 32 bit to save memory (we haven't seen a scenario in which there are more than 4B primary alignments), but can trivially extend to 64b if needed.
// if we ever want to extend to 64 bits, consider than number of PRIM groups is also limited to 1B unique CIGARs - see sam_seg_prim_add_sag_SA
typedef uint32_t SAGroup; 

typedef enum { BUDDY_NONE=0, BUDDY_MATE=1, BUDDY_SAGGY=2, BUDDY_EITHER=3 } BuddyType; // part of the file format

extern rom sag_type_name (SagType sagt);

// ZIP Stuff
COMPRESSOR_CALLBACK(sam_zip_seq); // used by longr codec
COMPRESSOR_CALLBACK(sam_zip_qual);
COMPRESSOR_CALLBACK(sam_zip_OQ);
COMPRESSOR_CALLBACK(sam_zip_TQ);
COMPRESSOR_CALLBACK(sam_zip_QX);
COMPRESSOR_CALLBACK(sam_zip_2Y);
COMPRESSOR_CALLBACK(sam_zip_U2);
COMPRESSOR_CALLBACK(sam_zip_BD_BI);
extern void sam_zip_initialize (void);
extern void sam_set_sag_type (void);
extern void sam_zip_finalize (bool is_last_user_txt_file);
extern bool sam_zip_dts_flag (int dts);
extern void sam_zip_after_compute (VBlockP vb);
extern void sam_zip_after_vbs (void);
extern void sam_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeader *vb_header);
extern void sam_sa_prim_finalize_ingest (void);

// HEADER stuff
extern bool sam_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);
extern void sam_header_finalize (void);
extern void sam_zip_free_end_of_z (void);
extern bool is_sam (STRp(header), bool *need_more);
extern bool is_bam (STRp(header), bool *need_more);

extern ContigPkgP sam_hdr_contigs;
extern uint32_t sam_num_header_contigs (void);

// ZIP/SEG Stuff
extern void sam_seg_initialize (VBlockP vb);
extern void sam_seg_finalize (VBlockP vb);
extern bool sam_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern bool sam_seg_is_big (ConstVBlockP vb, DictId dict_id);
extern rom sam_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern uint32_t sam_seg_seq_len_from_cigar (STRp(cigar));
extern uint32_t sam_seg_get_seq_len_by_MD_field (STRp(md_str));
extern void sam_zip_generate_recon_plan (void);
extern void sam_zip_init_vb (VBlockP vb);
extern void sam_zip_after_compress (VBlockP vb);
extern void sam_stats_reallocate (void);
extern void sam_zip_genozip_header (SectionHeaderGenozipHeader *header);
extern void sam_deep_merge (VBlockP vb);

// PIZ Stuff
extern void sam_piz_genozip_header (const SectionHeaderGenozipHeader *header);
extern void sam_piz_finalize (void);
extern IS_SKIP (sam_piz_is_skip_section);
extern bool sam_piz_maybe_reorder_lines (void);
extern bool sam_piz_init_vb (VBlockP vb, const SectionHeaderVbHeader *header, uint32_t *txt_data_so_far_single_0_increment);
extern void sam_piz_recon_init (VBlockP vb);
extern void sam_piz_after_recon (VBlockP vb);
extern void sam_piz_process_recon (VBlockP vb);
extern CONTAINER_FILTER_FUNC (sam_piz_filter);
extern void sam_reconstruct_SEQ_vs_ref (VBlockP vb, ContextP ctx, rom unused, unsigned unused2, ReconType reconstruct);
extern void sam_set_FLAG_filter (rom optarg);
extern void sam_set_MAPQ_filter (rom optarg);
extern void sam_piz_load_sags (void);
extern bool sam_piz_dispatch_one_load_sag_vb (Dispatcher dispatcher);
extern void sam_reconstruct_missing_quality (VBlockP vb, ReconType reconstruct);
extern void sam_piz_xtra_line_data (VBlockP vb);
extern void sam_piz_after_preproc (VBlockP vb);
extern CONTAINER_ITEM_CALLBACK (sam_piz_con_item_cb);

// BAM Stuff
extern void bam_seg_initialize (VBlockP vb);
extern int32_t bam_is_header_done (bool is_eof);
extern int32_t bam_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern void bam_read_vblock (VBlockP vb);
extern void bam_seg_initialize (VBlockP vb);
extern rom bam_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);

// SAM-to-FASTQ stuff
CONTAINER_CALLBACK (sam_piz_container_cb);

// VB stuff
extern void sam_vb_release_vb();
extern void sam_vb_destroy_vb();
extern unsigned sam_vb_size (DataType dt);
extern unsigned sam_vb_zip_dl_size (void);
extern void sam_reset_line (VBlockP vb);

// Special - used for SAM & BAM
#define SAM_SPECIAL { sam_cigar_special_CIGAR, sam_piz_special_TLEN_old, sam_piz_special_BD_BI, sam_piz_special_AS_old, \
                      sam_piz_special_MD_old, bam_piz_special_FLOAT, bam_piz_special_BIN, sam_piz_special_NM,   \
                      sam_piz_special_MD, sam_piz_special_REF_CONSUMED, \
                      sam_piz_special_PNEXT_IS_PREV_POS_old, sam_piz_special_COPY_MATE_FLAG, sam_piz_special_COPY_MATE_TLEN_old, \
                      sam_piz_special_COPY_BUDDY_CIGAR, sam_piz_special_FASTQ_CONSUME_AUX, sam_piz_special_TLEN, sam_piz_special_QUAL,\
                      sam_piz_special_pull_from_sag, sam_piz_special_SEQ, \
                      sam_piz_special_PRIM_QNAME, sam_piz_special_SQUANK, \
                      sam_piz_special_BSSEEKER2_XO, sam_piz_special_BSSEEKER2_XG, sam_piz_special_BSSEEKER2_XM,\
                      sam_piz_special_SA_main, sam_piz_special_COPY_PRIM, \
                      sam_piz_special_BWA_XC, sam_piz_special_BWA_XT, sam_piz_special_BWA_X1, sam_piz_special_BWA_XS, \
                      sam_piz_special_SM, sam_piz_special_AM, sam_piz_special_PNEXT, \
                      sam_piz_special_DEMUX_BY_MATE, sam_piz_special_DEMUX_BY_MATE_PRIM, sam_piz_special_DEMUX_BY_BUDDY, \
                      sam_piz_special_GEM3_XB, sam_piz_special_BSBOLT_YS, sam_piz_special_OA_RNAME, \
                      sam_piz_special_BISMARK_XG, sam_piz_special_HI, sam_piz_special_DEMUX_BY_BUDDY_MAP, sam_piz_special_SEQ_LEN,\
                      sam_piz_special_FI, sam_piz_special_cm, sam_piz_special_COPY_BUDDY, sam_piz_special_SET_BUDDY, \
                      sam_piz_special_TX_AN_POS, sam_piz_special_COPY_TEXTUAL_CIGAR, sam_piz_special_BISMARK_XM, \
                      sam_piz_special_BSBOLT_XB, sam_piz_special_UQ \
                    }
SPECIAL (SAM, 0,  CIGAR,                 sam_cigar_special_CIGAR);
SPECIAL (SAM, 1,  TLEN_old,              sam_piz_special_TLEN_old);            // used up to 12.0.42
SPECIAL (SAM, 2,  BDBI,                  sam_piz_special_BD_BI);
SPECIAL (SAM, 3,  AS_old,                sam_piz_special_AS_old);              // Reconstructs seq_len minus delta. Note: called "AS" until 12.0.37 and used to reconstruct AS:i
SPECIAL (SAM, 4,  MD_old,                sam_piz_special_MD_old);              // used in files compressed with Genozip up to 12.0.36
SPECIAL (SAM, 5,  FLOAT,                 bam_piz_special_FLOAT);               // used in BAM to represent float optional values
SPECIAL (SAM, 6,  BIN,                   bam_piz_special_BIN);   
SPECIAL (SAM, 7,  NM,                    sam_piz_special_NM);                  // introduced 12.0.37
SPECIAL (SAM, 8,  MD,                    sam_piz_special_MD);                  // introduced 12.0.37
SPECIAL (SAM, 9,  REF_CONSUMED,          sam_piz_special_REF_CONSUMED);        // introduced 12.0.41: Reconstructs based on ref_consumed
SPECIAL (SAM, 10, PNEXT_IS_PREV_POS_old, sam_piz_special_PNEXT_IS_PREV_POS_old); // introduced 12.0.41 up to v13
SPECIAL (SAM, 11, COPY_MATE_FLAG,        sam_piz_special_COPY_MATE_FLAG);      // introduced 12.0.41
SPECIAL (SAM, 12, COPY_MATE_TLEN_old,    sam_piz_special_COPY_MATE_TLEN_old);  // Used in 12.0.41 and 12.0.42
SPECIAL (SAM, 13, COPY_BUDDY_CIGAR,      sam_piz_special_COPY_BUDDY_CIGAR);    // introduced 12.0.41
SPECIAL (SAM, 14, FASTQ_CONSUME_AUX,     sam_piz_special_FASTQ_CONSUME_AUX);   // introduced 12.0.41 (up to 14.0.25 called CONSUME_MC_Z)
SPECIAL (SAM, 15, TLEN,                  sam_piz_special_TLEN);                // introduced 13.0.1
SPECIAL (SAM, 16, QUAL,                  sam_piz_special_QUAL);                // introduced 13.0.1
SPECIAL (SAM, 17, SAG,                   sam_piz_special_pull_from_sag);       // introduced 14.0.0
SPECIAL (SAM, 18, SEQ,                   sam_piz_special_SEQ);                 // introduced 14.0.0
SPECIAL (SAM, 19, PRIM_QNAME,            sam_piz_special_PRIM_QNAME);          // introduced 14.0.0
SPECIAL (SAM, 20, SQUANK,                sam_piz_special_SQUANK);              // introduced 14.0.0
SPECIAL (SAM, 21, BSSEEKER2_XO,          sam_piz_special_BSSEEKER2_XO);        // introduced 14.0.0
SPECIAL (SAM, 22, BSSEEKER2_XG,          sam_piz_special_BSSEEKER2_XG);        // introduced 14.0.0
SPECIAL (SAM, 23, BSSEEKER2_XM,          sam_piz_special_BSSEEKER2_XM);        // introduced 14.0.0
SPECIAL (SAM, 24, SA_main,               sam_piz_special_SA_main);             // introduced 14.0.0
SPECIAL (SAM, 25, COPY_PRIM,             sam_piz_special_COPY_PRIM);           // introduced 14.0.0
SPECIAL (SAM, 26, BWA_XC,                sam_piz_special_BWA_XC);              // introduced 14.0.0
SPECIAL (SAM, 27, BWA_XT,                sam_piz_special_BWA_XT);              // introduced 14.0.0
SPECIAL (SAM, 28, BWA_X1,                sam_piz_special_BWA_X1);              // introduced 14.0.0
SPECIAL (SAM, 29, BWA_XS,                sam_piz_special_BWA_XS);              // introduced 14.0.0
SPECIAL (SAM, 30, SM,                    sam_piz_special_SM);                  // introduced 14.0.0
SPECIAL (SAM, 31, AM,                    sam_piz_special_AM);                  // introduced 14.0.0
SPECIAL (SAM, 32, PNEXT,                 sam_piz_special_PNEXT);               // introduced 14.0.0
SPECIAL (SAM, 33, DEMUX_BY_MATE,         sam_piz_special_DEMUX_BY_MATE);       // introduced 14.0.0
SPECIAL (SAM, 34, DEMUX_BY_MATE_PRIM,    sam_piz_special_DEMUX_BY_MATE_PRIM);  // introduced 14.0.0
SPECIAL (SAM, 35, DEMUX_BY_BUDDY,        sam_piz_special_DEMUX_BY_BUDDY);      // introduced 14.0.0
SPECIAL (SAM, 36, GEM3_XB,               sam_piz_special_GEM3_XB);             // introduced 14.0.0
SPECIAL (SAM, 37, BSBOLT_YS,             sam_piz_special_BSBOLT_YS);           // introduced 14.0.0
SPECIAL (SAM, 38, OA_RNAME,              sam_piz_special_OA_RNAME);            // introduced 14.0.0
SPECIAL (SAM, 39, BISMARK_XG,            sam_piz_special_BISMARK_XG);          // introduced 14.0.0
SPECIAL (SAM, 40, HI,                    sam_piz_special_HI);                  // introduced 14.0.0
SPECIAL (SAM, 41, DEMUX_BY_BUDDY_MAP,    sam_piz_special_DEMUX_BY_BUDDY_MAP);  // introduced 14.0.0
SPECIAL (SAM, 42, SEQ_LEN,               sam_piz_special_SEQ_LEN);             // introduced 14.0.0
SPECIAL (SAM, 43, FI,                    sam_piz_special_FI);                  // introduced 14.0.0
SPECIAL (SAM, 44, cm,                    sam_piz_special_cm);                  // introduced 14.0.0
SPECIAL (SAM, 45, COPY_BUDDY,            sam_piz_special_COPY_BUDDY);          // introduced 14.0.0
SPECIAL (SAM, 46, SET_BUDDY,             sam_piz_special_SET_BUDDY);           // introduced 14.0.0
SPECIAL (SAM, 47, TX_AN_POS,             sam_piz_special_TX_AN_POS);           // introduced 14.0.0
SPECIAL (SAM, 48, COPY_TEXTUAL_CIGAR,    sam_piz_special_COPY_TEXTUAL_CIGAR);  // introduced 14.0.0
SPECIAL (SAM, 49, BISMARK_XM,            sam_piz_special_BISMARK_XM);          // introduced 14.0.0
SPECIAL (SAM, 50, BSBOLT_XB,             sam_piz_special_BSBOLT_XB);           // introduced 14.0.0
SPECIAL (SAM, 51, UQ,                    sam_piz_special_UQ);                  // introduced 14.0.10
#define NUM_SAM_SPECIAL 52
 
#define SAM_LOCAL_GET_LINE_CALLBACKS                          \
    { DT_SAM, _OPTION_BD_BI,       sam_zip_BD_BI           }, \
    { DT_SAM, _SAM_QUAL,           sam_zip_qual            }, \
    { DT_SAM, _OPTION_U2_Z,        sam_zip_U2              }, \
    { DT_SAM, _OPTION_OQ_Z,        sam_zip_OQ              }, \
    { DT_SAM, _OPTION_TQ_Z,        sam_zip_TQ              }, \
    { DT_SAM, _OPTION_QX_Z,        sam_zip_QX              }, \
    { DT_SAM, _OPTION_2Y_Z,        sam_zip_2Y              }, 

#define BAM_LOCAL_GET_LINE_CALLBACKS                          \
    { DT_BAM, _OPTION_BD_BI,       sam_zip_BD_BI           }, \
    { DT_BAM, _SAM_QUAL,           sam_zip_qual            }, \
    { DT_BAM, _OPTION_U2_Z,        sam_zip_U2              }, \
    { DT_BAM, _OPTION_OQ_Z,        sam_zip_OQ              }, \
    { DT_BAM, _OPTION_TQ_Z,        sam_zip_TQ              }, \
    { DT_BAM, _OPTION_QX_Z,        sam_zip_QX              }, \
    { DT_BAM, _OPTION_2Y_Z,        sam_zip_2Y              }, 
    
// Important: Numbers (and order) of translators cannot be changed, as they are part of the file format
// (they are included in the TOP2BAM container)
// translator numbers must start from 1 - 0 is reserved for "none"
TRANSLATOR (SAM, BAM,   1,  I8,         container_translate_I8)     // reconstruct binary little endian functions
TRANSLATOR (SAM, BAM,   2,  U8,         container_translate_U8)     // 
TRANSLATOR (SAM, BAM,   3,  LTEN_I16,   container_translate_LTEN_I16) 
TRANSLATOR (SAM, BAM,   4,  LTEN_U16,   container_translate_LTEN_U16) 
TRANSLATOR (SAM, BAM,   5,  LTEN_I32,   container_translate_LTEN_I32) 
TRANSLATOR (SAM, BAM,   6,  LTEN_U32,   container_translate_LTEN_U32) 
TRANSLATOR (SAM, BAM,   7,  FLOAT,      sam_piz_sam2bam_FLOAT)      // reconstructs SAM-stored textual floating point as little endian 32bit float
TRANSLATOR (SAM, BAM,   8,  ARRAY_SELF, sam_piz_sam2bam_ARRAY_SELF) // remove the comma from the prefix that contains the type, eg "i,"->"i"
TRANSLATOR (SAM, BAM,   9,  RNAME,      sam_piz_sam2bam_RNAME)      // reconstructs the b250 index or -1 if "*"
TRANSLATOR (SAM, BAM,   10, POS,        sam_piz_sam2bam_POS)        // reconstructs Little Endian U32 0-based POS. 
TRANSLATOR (SAM, BAM,   11, SEQ,        sam_piz_sam2bam_SEQ)        // textual SEQ to BAM-format SEQ 
TRANSLATOR (SAM, BAM,   12, QUAL,       sam_piz_sam2bam_QUAL)       // textual QUAL to BAM-format QUAL 
TRANSLATOR (SAM, BAM,   13, TLEN,       sam_piz_sam2bam_TLEN)       // place TLEN last_value in BAM alignment 
TRANSLATOR (SAM, BAM,   14, AUX,        sam_piz_sam2bam_AUX)        // used up to v11, kept for for backward compatability as old files expect it
TRANSLATOR (SAM, BAM,   15, AUX_SELF,   sam_piz_sam2bam_AUX_SELF)   // transform prefixes in Aux Container from SAM to BAM format 
TRANSLATOR (SAM, FASTQ, 16, SEQ,        sam_piz_sam2fastq_SEQ)      // reverse-complement the sequence if needed, and drop if "*"
TRANSLATOR (SAM, FASTQ, 17, QUAL,       sam_piz_sam2fastq_QUAL)     // reverse the QUAL if reverse-complemented and drop fastq records with QUAL="*"
TRANSLATOR (SAM, FASTQ, 18, FLAG,       sam_piz_sam2fastq_FLAG)     // emit 1 if (FLAGS & 0x40) or 2 of (FLAGS & 0x80)
TRANSLATOR (SAM, BAM,   19, QNAME,      sam_piz_sam2bam_QNAME)      // v15: if --deep: add QNAME to vb->deep_ents in case of BAM
#define NUM_SAM_TRANS   20 // including "none"

#define SAM_TRANSLATORS { NULL /* none */, container_translate_I8, container_translate_U8, container_translate_LTEN_I16, \
                          container_translate_LTEN_U16, container_translate_LTEN_I32, container_translate_LTEN_U32, \
                          sam_piz_sam2bam_FLOAT, sam_piz_sam2bam_ARRAY_SELF, sam_piz_sam2bam_RNAME, sam_piz_sam2bam_POS, sam_piz_sam2bam_SEQ, \
                          sam_piz_sam2bam_QUAL, sam_piz_sam2bam_TLEN, sam_piz_sam2bam_AUX, sam_piz_sam2bam_AUX_SELF, \
                          sam_piz_sam2fastq_SEQ, sam_piz_sam2fastq_QUAL, sam_piz_sam2fastq_FLAG, sam_piz_sam2bam_QNAME }

TXTHEADER_TRANSLATOR (sam_header_bam2sam);
TXTHEADER_TRANSLATOR (sam_header_sam2bam);
TXTHEADER_TRANSLATOR (txtheader_sam2fq);

#define SAM_DICT_ID_ALIASES                    \
    /*         alias         maps to       */  \
    { DT_SAM,  _OPTION_UB_Z, _OPTION_BX_Z },   /* 10xGenomics tags - CellRanger: UB/UR/UY = LongRanger BX/RX/QX */ \
    { DT_SAM,  _OPTION_UR_Z, _OPTION_RX_Z },   \
    { DT_SAM,  _OPTION_UY_Z, _OPTION_QX_Z },  

#define BAM_DICT_ID_ALIASES                    \
    /*         alias         maps to      */   \
    { DT_BAM,  _OPTION_UB_Z, _OPTION_BX_Z },   \
    { DT_BAM,  _OPTION_UR_Z, _OPTION_RX_Z },   \
    { DT_BAM,  _OPTION_UY_Z, _OPTION_QX_Z },  

#define SAM_CONTIG_FMT "@SQ	SN:%.*s	LN:%"PRId64

typedef enum { SAM_COMP_NONE=255, SAM_COMP_MAIN=0, SAM_COMP_PRIM=1, SAM_COMP_DEPN=2, SAM_COMP_FQ00=3, SAM_COMP_FQ01=4 } SamComponentType;
#define SAM_COMP_NAMES { "MAIN", "PRIM", "DEPN", "FQ00", "FQ01" }

#define sam_is_main_vb (vb->comp_i == SAM_COMP_MAIN)
#define sam_is_prim_vb (vb->comp_i == SAM_COMP_PRIM)
#define sam_is_depn_vb (vb->comp_i == SAM_COMP_DEPN)

#define dict_id_is_aux_sf dict_id_is_type_2
#define dict_id_aux_sf dict_id_type_2

#define IS_BAM_ZIP TXT_DT(BAM)
#define IS_SRC_BAM_PIZ (Z_DT(SAM) && z_file->z_flags.txt_is_bin)       // source file was BAM, NOT that the reconstruction is BAM! 
#define IS_RECON_BAM (flag.out_dt == DT_BAM)                              // reconstructed file is BAM, regardless of source file type 
#define IS_SRC_and_RECON_BAM_PIZ (IS_SRC_BAM_PIZ && vb->translation.trans_containers) 
#define IS_SRC_BAM (command==ZIP ? IS_BAM_ZIP : IS_SRC_BAM_PIZ)
