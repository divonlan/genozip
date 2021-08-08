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

// Fields
#define _VCF_POS                DICT_ID_MAKEF_4 ("POS")    
#define _VCF_CHROM              DICT_ID_MAKEF_5 ("CHROM")
#define _VCF_ID                 DICT_ID_MAKEF_2 ("ID")
#define _VCF_REFALT             DICT_ID_MAKEF_7 ("REF+ALT")
#define _VCF_QUAL               DICT_ID_MAKEF_4 ("QUAL")
#define _VCF_FILTER             DICT_ID_MAKEF_6 ("FILTER")
#define _VCF_INFO               DICT_ID_MAKEF_4 ("INFO")
#define _VCF_FORMAT             DICT_ID_MAKEF_6 ("FORMAT")
#define _VCF_SAMPLES            DICT_ID_MAKEF_7 ("SAMPLES")
#define _VCF_EOL                DICT_ID_MAKEF_3 ("EOL")
#define _VCF_TOPLEVEL           DICT_ID_MAKEF_L (TOPLEVEL)
#define _VCF_oCHROM             DICT_ID_MAKEF_6 ("oCHROM")
#define _VCF_oPOS               DICT_ID_MAKEF_4 ("oPOS")
#define _VCF_oREFALT            DICT_ID_MAKEF_7 ("oREFALT")
#define _VCF_oXSTRAND           DICT_ID_MAKEF_L ("oXSTRAND")
#define _VCF_COORDS             DICT_ID_MAKEF_6 ("COORDS")
#define _VCF_oSTATUS            DICT_ID_MAKEF_7 ("o$TATUS")
#define _VCF_COPYPOS            DICT_ID_MAKEF_7 ("C0PYPOS")
#define _VCF_LIFT_REF           DICT_ID_MAKEF_L ("LIFT_REF")
#define _VCF_COPYSTAT           DICT_ID_MAKEF_L ("CoPYSTAT")
#define _VCF_TOPLUFT            DICT_ID_MAKEF_7 ("ToPLUFT")
#define _VCF_LINE_NUM           DICT_ID_MAKEF_L ("LINE_NUM")

// FORMAT fields
#define _FORMAT_AD              DICT_ID_MAKE2_2 ("AD")     // <ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
#define _FORMAT_ADF             DICT_ID_MAKE2_3 ("ADF")    // <ID=ADF,Number=R,Type=Float,Description="Allele dosage on fwd strand">
#define _FORMAT_ADR             DICT_ID_MAKE2_3 ("ADR")    // <ID=ADR,Number=R,Type=Float,Description="Allele dosage on rev strand">
#define _FORMAT_ADALL           DICT_ID_MAKE2_4 ("ADALL")  // from GIAB: <ID=ADALL,Number=R,Type=Integer,Description="Net allele depths across all datasets">
#define _FORMAT_AF              DICT_ID_MAKE2_2 ("AF")     // <ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed">
#define _FORMAT_DP              DICT_ID_MAKE2_2 ("DP")     // <ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
#define _FORMAT_DS              DICT_ID_MAKE2_2 ("DS")
#define _FORMAT_GL              DICT_ID_MAKE2_2 ("GL")
#define _FORMAT_GP              DICT_ID_MAKE2_2 ("GP")     // <ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification">
#define _FORMAT_GQ              DICT_ID_MAKE2_2 ("GQ")     // <ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#define _FORMAT_GT              DICT_ID_MAKE2_2 ("GT")     // <ID=GT,Number=1,Type=String,Description="Genotype">
#define _FORMAT_PL              DICT_ID_MAKE2_2 ("PL")     // <ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
#define _FORMAT_PRI             DICT_ID_MAKE2_3 ("PRI")    // <ID=PRI,Number=G,Type=Float,Description="Phred-scaled prior probabilities for genotypes">
#define _FORMAT_F1R2            DICT_ID_MAKE2_4 ("F1R2")   // <ID=F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting each allele"> see: https://github.com/broadinstitute/gatk/blob/master/docs/mutect/mutect.pdf
#define _FORMAT_F2R1            DICT_ID_MAKE2_4 ("F2R1")   // <ID=F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting each allele">
#define _FORMAT_MB              DICT_ID_MAKE2_2 ("MB")     // <ID=MB,Number=4,Type=Integer,Description="Per-sample component statistics to detect mate bias">
#define _FORMAT_PP              DICT_ID_MAKE2_2 ("PP")     // <ID=PP,Number=G,Type=Integer,Description="Phred-scaled genotype posterior probabilities rounded to the closest integer">
#define _FORMAT_SAC             DICT_ID_MAKE2_3 ("SAC")    // <ID=SAC,Number=.,Type=Integer,Description="Number of reads on the forward and reverse strand supporting each allele (including reference)">
#define _FORMAT_SB              DICT_ID_MAKE2_2 ("SB")     // <ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias">
#define _FORMAT_PS              DICT_ID_MAKE2_2 ("PS") 

// PBWT fields - same dict_id for all data types using PBWT, as codec_pbwt_uncompress relies on it
#ifndef _PBWT_RUNS
#define _PBWT_RUNS              DICT_ID_MAKE2_L ("@1BWTRUN") // PBWT runs
#define _PBWT_FGRC              DICT_ID_MAKE2_L ("@2BWTFGR") // PBWT foreground run count
#define _PBWT_HT_MATRIX         DICT_ID_MAKE2_3 ("@HT")    
#define _PBWT_GT_HT_INDEX       DICT_ID_MAKE2_L ("@INDEXHT") // different first 2 letters
#endif

#define _FORMAT_PBWT_RUNS       _PBWT_RUNS
#define _FORMAT_PBWT_FGRC       _PBWT_FGRC
#define _FORMAT_GT_HT           _PBWT_HT_MATRIX
#define _FORMAT_GT_HT_INDEX     _PBWT_GT_HT_INDEX

// INFO fields
#define _INFO_AC                DICT_ID_MAKE1_2 ("AC")     // <ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
#define _INFO_AF                DICT_ID_MAKE1_2 ("AF")     // <ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
#define _INFO_AN                DICT_ID_MAKE1_2 ("AN")     // <ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#define _INFO_AA                DICT_ID_MAKE1_2 ("AA")     // <ID=AA,Number=1,Type=String,Description="Ancestral Allele"> - defined in the VCF specification
#define _INFO_BaseCounts        DICT_ID_MAKE1_L ("BaseCounts")
#define _INFO_DP                DICT_ID_MAKE1_2 ("DP")     // <ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
#define _INFO_DP4               DICT_ID_MAKE1_3 ("DP4")   
#define _INFO_SF                DICT_ID_MAKE1_2 ("SF")
#define _INFO_VQSLOD            DICT_ID_MAKE1_6 ("VQSLOD") // <ID=VQSLOD,Number=1,Type=Float,Description="Log odds of being a true variant versus being false under the trained Gaussian mixture model">

// Added by GATK HaplotypeCaller in a gVCF: https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
#define _INFO_END               DICT_ID_MAKE1_3 ("END")    // <ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
#define _INFO_MLEAC             DICT_ID_MAKE1_5 ("MLEAC")  // <ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
#define _INFO_MLEAF             DICT_ID_MAKE1_5 ("MLEAF")  // <ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
#define _INFO_LDAF              DICT_ID_MAKE1_4 ("LDAF")   //  MLE Allele Frequency Accounting for LD
#define _INFO_MQ0               DICT_ID_MAKE1_3 ("MQ0")    
#define _FORMAT_MIN_DP          DICT_ID_MAKE1_6 ("MIN_DP") // <ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">

// Ensembl VEP (Variant Effect Predictor) fields: https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html
#define _INFO_vep               DICT_ID_MAKE1_3 ("vep")
#define _INFO_AGE_HISTOGRAM_HET DICT_ID_MAKE1_L ("AGE_HISTOGRAM_HET") 
#define _INFO_AGE_HISTOGRAM_HOM DICT_ID_MAKE1_L ("AGE_HISTOGRAM_HOM")
#define _INFO_CSQ               DICT_ID_MAKE1_3 ("CSQ")    // "Consequences" 
#define _INFO_MAX_AF            DICT_ID_MAKE1_6 ("MAX_AF") // highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD
        
// clinvar
#define _INFO_ALLELEID          DICT_ID_MAKE1_L ("ALLELEID")// <ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
#define _INFO_CLNDN             DICT_ID_MAKE1_5 ("CLNDN")  // <ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
#define _INFO_RS                DICT_ID_MAKE1_2 ("RS")     // <ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
#define _INFO_CLNHGVS           DICT_ID_MAKE1_7 ("CLNHGVS")// <ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
#define _INFO_CLNHGVS_pos       DICT_ID_MAKE1_L ("CpLNHGVS")
#define _INFO_CLNHGVS_refalt    DICT_ID_MAKE1_L ("CrLNHGVS")

// ExAC fields
#define _INFO_DP_HIST           DICT_ID_MAKE1_7 ("DP_HIST")// from ExAC: Depth (DP) histogram in 20 equal intervals between 0-100 : See https://www.biorxiv.org/content/biorxiv/suppl/2015/10/30/030338.DC1/030338-1.pdf
#define _INFO_GQ_HIST           DICT_ID_MAKE1_7 ("GQ_HIST")// from ExAC: Genotype Quality (GQ) histogram in 20 equal intervals between 0-100

// Structural variants (also uses INFO/END): https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/VCF%20(Variant%20Call%20Format)%20version%204.0/encoding-structural-variants/
#define _INFO_SVLEN             DICT_ID_MAKE1_5 ("SVLEN")

// genozip INFO fields
#define _INFO_LUFT              DICT_ID_MAKE1_4 (INFO_LUFT_NAME)
#define _INFO_PRIM              DICT_ID_MAKE1_4 (INFO_PRIM_NAME)
#define _INFO_LREJ              DICT_ID_MAKE1_4 (INFO_LREJ_NAME)
#define _INFO_PREJ              DICT_ID_MAKE1_4 (INFO_PREJ_NAME)

// did_i to dict_i mapping - only needed for did_i's referred to explicitly
// the CHROM field MUST be the first field (because of ctx_build_zf_ctx_from_contigs)
typedef enum { VCF_CHROM, VCF_POS, VCF_ID, VCF_REFALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT, VCF_SAMPLES, VCF_EOL, VCF_TOPLEVEL, 
               VCF_oCHROM, VCF_oPOS, VCF_oREFALT, VCF_oXSTRAND, VCF_COORDS, VCF_oSTATUS, VCF_COPYPOS, VCF_LIFT_REF, VCF_COPYSTAT, VCF_TOPLUFT, VCF_LINE_NUM, // Liftover data - must appear in same order in any data type that has it
               INFO_LUFT, INFO_PRIM, INFO_LREJ, INFO_PREJ, INFO_SF, INFO_CLNHGVS_pos, INFO_CLNHGVS_refalt, INFO_AF, INFO_DP,
               FORMAT_GT, FORMAT_GT_HT, FORMAT_GT_HT_INDEX, FORMAT_ADALL, FORMAT_AD, FORMAT_F2R1, 
               FORMAT_PBWT_RUNS, FORMAT_PBWT_FGRC, // RUNS must be before FGRC so it is emitted in the file in this order
               NUM_VCF_FIELDS } VcfFields;

#define VCF_MAPPING { V(VCF_CHROM), V(VCF_POS), V(VCF_ID), V(VCF_REFALT), V(VCF_QUAL), V(VCF_FILTER), V(VCF_INFO), V(VCF_FORMAT), V(VCF_SAMPLES), \
                      V(VCF_EOL), V(VCF_TOPLEVEL), V(VCF_oCHROM), V(VCF_oPOS), V(VCF_oREFALT), V(VCF_oXSTRAND), V(VCF_COORDS), V(VCF_oSTATUS), V(VCF_COPYPOS), \
                      V(VCF_LIFT_REF), V(VCF_COPYSTAT), V(VCF_TOPLUFT), V(VCF_LINE_NUM),\
                      V(INFO_LUFT), V(INFO_PRIM), V(INFO_LREJ), V(INFO_PREJ), V(INFO_SF), V(INFO_CLNHGVS_pos), V(INFO_CLNHGVS_refalt), V(INFO_AF), V(INFO_DP), \
                      V(FORMAT_GT), V(FORMAT_GT_HT), V(FORMAT_GT_HT_INDEX), V(FORMAT_ADALL), V(FORMAT_AD), V(FORMAT_F2R1), V(FORMAT_PBWT_RUNS), V(FORMAT_PBWT_FGRC), }

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
                      vcf_piz_special_INFO_HGVS_POS, vcf_piz_special_INFO_HGVS_REFALT }
SPECIAL (VCF, 0,  main_REFALT,  vcf_piz_special_main_REFALT);
SPECIAL (VCF, 1,  FORMAT,       vcf_piz_special_FORMAT)
SPECIAL (VCF, 2,  AC,           vcf_piz_special_INFO_AC);
SPECIAL (VCF, 3,  SVLEN,        vcf_piz_special_INFO_SVLEN);
SPECIAL (VCF, 4,  DS,           vcf_piz_special_FORMAT_DS);
SPECIAL (VCF, 5,  BaseCounts,   vcf_piz_special_INFO_BaseCounts);
SPECIAL (VCF, 6,  SF,           vcf_piz_special_INFO_SF);
SPECIAL (VCF, 7,  MINUS,        vcf_piz_special_MINUS);            // added v12.0.0
SPECIAL (VCF, 8,  LIFT_REF,     vcf_piz_special_LIFT_REF);         // added v12.0.0
SPECIAL (VCF, 9,  COPYSTAT,     vcf_piz_special_COPYSTAT);         // added v12.0.0
SPECIAL (VCF, 10, other_REFALT, vcf_piz_special_other_REFALT);     // added v12.0.0
SPECIAL (VCF, 11, COPYPOS,      vcf_piz_special_COPYPOS);          // added v12.0.0
SPECIAL (VCF, 12, ALLELE,       vcf_piz_special_ALLELE);           // added v12.0.0
SPECIAL (VCF, 13, HGVS_POS,     vcf_piz_special_INFO_HGVS_POS);    // added v12.0.15
SPECIAL (VCF, 14, HGVS_REFALT,  vcf_piz_special_INFO_HGVS_REFALT); // added v12.0.15
#define NUM_VCF_SPECIAL 15

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
     (ltrans_props[(ctx)->luft_trans].upon == TW_XSTRAND        && LO_IS_OK (last_ostatus) && *CTX(VCF_oXSTRAND)->last_snip == 'X')))

#define VCF_DICT_ID_ALIASES \
    /*         alias             maps to this ctx  */  \
    { DT_VCF,  _INFO_END, _VCF_POS    }, \

#define VCF_LOCAL_GET_LINE_CALLBACKS

#define dict_id_is_vcf_info_sf   dict_id_is_type_1
#define dict_id_is_vcf_format_sf dict_id_is_type_2

enum { oCHROM, oPOS, oREF, oXSTRAND, oSTATUS }; // order of fields as defined in data_types.h
#define ODID(offset) (DTFZ(ochrom)+(DidIType)(offset))

#endif
