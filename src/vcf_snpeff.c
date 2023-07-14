// ------------------------------------------------------------------
//   vcf_snpeff.c
//   Copyright (C) 2022-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "reconstruct.h"

// --------
// INFO/ANN
// --------

// ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
// See: https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
// example: ANN=T|intergenic_region|MODIFIER|U2|ENSG00000277248|intergenic_region|ENSG00000277248|||n.10510103A>T||||||
void vcf_seg_INFO_ANN (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    static const MediumContainer ann = {
        .nitems_lo   = 16, 
        .drop_final_repsep = true,
        .repsep      = {','},
        .items       = { { .dict_id.num=_INFO_ANN_Allele,                       .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_L("Annotation"),          .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_L("Annotation_Impact"),   .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_L("Gene_Name"),           .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_7("Gene_ID"),             .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_L("Feature_Type"),        .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_L("Feature_ID"),          .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_L("Transcript_BioType"),  .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_4("Rank"),                .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_6("HGVS_c"),              .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_6("HGVS_p"),              .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_4("cDNA"),                .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_3("CDS"),                 .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_2("AA"),                  .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_L("Distance"),            .separator = {'|'} }, 
                         { .dict_id.num=DICT_ID_MAKE1_6("Errors")                                  } } };

    seg_array_of_struct (VB, ctx, ann, STRa(value), (SegCallback[]){vcf_seg_INFO_allele,0,0,0,0,0,0,0,0,vcf_seg_INFO_HGVS,0,0,0,0,0,0}, value_len);
}

// See: https://pcingola.github.io/SnpEff/se_inputoutput/#eff-field-vcf-output-files
// Example: INTERGENIC(MODIFIER||||||||||1),INTERGENIC(MODIFIER||||||||||2),UPSTREAM(MODIFIER|||||Gene_Qt01t0103210-03|||Os01t013410-23||1),UPSTREAM(MODIFIER|||||Gene_Qt01t0100231-33|||Os01t0102314-31||2)
void vcf_seg_INFO_EFF (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len);
    // seg_add_to_local_text (VB, ctx, STRa(value), LOOKUP_NONE, value_len);
    // return;
#if 0    
    // TODO: move to this after we improve to exploit correlations between fields
    bool is_xstream = value_len > 10 && (!memcmp (value, "UPSTREAM", 8) || !memcmp (value, "DOWNSTREAM", 10));

    MediumContainer eff = {
        .nitems_lo   = 12, 
        .repsep      = {','},
        .drop_final_repsep = true,
        .items       = { { .dict_id={ _INFO_EFF_Effect             }, .separator = {'('} }, 
                         { .dict_id={ _INFO_EFF_Effect_impact      }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_EFF_Functional_Class   }, .separator = {'|'} }, 
                         { .dict_id={ is_xstream ? _INFO_EFF_Distance : _INFO_EFF_Codon_Change  }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_EFF_Amino_Acid_Change  }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_EFF_Amino_Acid_Length  }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_EFF_Gene_Name          }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_EFF_Transcript_BioType }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_EFF_Gene_Coding        }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_EFF_Transcript_ID      }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_EFF_Exon_Intron_Rank   }, .separator = {'|'} }, 
                         { .dict_id={ _INFO_EFF_Genotype_Number    }, .separator = {')'} }, 
                         
                        // TODO: handle case where 1 or more of the array have an extra "Warnings_Errors" field
                        // { .dict_id={ _INFO_EFF_Warnings_Errors       }, .separator = {')'} }, 
                       } };

    seg_array_of_struct (VB, ctx, eff, STRa(value), NULL, value_len);
#endif
}
