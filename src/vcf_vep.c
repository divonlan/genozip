// ------------------------------------------------------------------
//   vcf_vep.c
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

// Handle VEP fields: https://www.ensembl.org/info/docs/tools/vep/index.html

static MediumContainer csq_con;
static SegCallback csq_cbs[MEDIUM_CON_NITEMS];

uint64_t _INFO_Allele; // needs to be translated in DVCF
static uint64_t _INFO_Existing_variation, _INFO_cDNA_position, _INFO_CDS_position, _INFO_Protein_position, _INFO_DISTANCE;

// called from vcf_inspect_txt_header_zip, NOT vcf_info_zip_initialize
void vcf_vep_zip_initialize (rom spec, rom buf_1st) // nul-terminated string containing list of fields. list of fields expected to be terminated by a '"'
{
    _INFO_Allele             = dict_id_make ("Allele",             STRLEN("Allele"),             DTYPE_VCF_INFO).num;
    _INFO_Existing_variation = dict_id_make ("Existing_variation", STRLEN("Existing_variation"), DTYPE_VCF_INFO).num;
    _INFO_cDNA_position      = dict_id_make ("cDNA_position",      STRLEN("cDNA_position"),      DTYPE_VCF_INFO).num;
    _INFO_CDS_position       = dict_id_make ("CDS_position",       STRLEN("CDS_position"),       DTYPE_VCF_INFO).num;
    _INFO_Protein_position   = dict_id_make ("Protein_position",   STRLEN("Protein_position"),   DTYPE_VCF_INFO).num;
    _INFO_DISTANCE           = dict_id_make ("CSQ_DISTANCE",       STRLEN("CSQ_DISTANCE"),       DTYPE_VCF_INFO).num;

    rom after = strchr (spec, '"');
    ASSERT (after, "VCF header errror: INFO/CSQ Description does not end with a quote:\n%.2000s", spec);
    int spec_len = after - spec;

    str_split (spec, spec_len, 0, '|', name, false);

    if (n_names > MEDIUM_CON_NITEMS) {
        WARN ("FYI: VEP CSQ field has %u annotations, but compression will be sub-optimal as it has more than %u fields. - please report to " EMAIL_SUPPORT ":\n%.*s", n_names, CONTAINER_MAX_DICTS, spec_len, spec);    
        segconf.vcf_is_vep = false;
        return;
    }

    memset (&csq_con, 0, sizeof (csq_con));
    memset (csq_cbs,  0, sizeof (csq_cbs));
    csq_con.nitems_lo = n_names;
    csq_con.drop_final_repsep = true;
    csq_con.repsep[0] = ',';

    // check if field is CSQ or vep
    bool is_vep = false;
    for (rom c=spec-1; c > buf_1st; c--)
        if (*c == '#') {
            SAFE_NUL(spec);
            is_vep = !!strstr (c, "ID=vep");
            SAFE_RESTORE;
            break;
        }

    for (unsigned i=0; i < n_names; i++) {
        csq_con.items[i] = (ContainerItem){ .dict_id = dict_id_make (names[i], name_lens[i], DTYPE_VCF_INFO), .separator[0] = ((i < n_names-1) ? '|'  : 0) };

        ContextP zctx = ctx_add_new_zf_ctx_from_txtheader (names[i], name_lens[i], csq_con.items[i].dict_id, TRANS_ID_NONE); // zctx, but will be the same did for vb contexts because pre-created
        if (!zctx) // already exist - predefined
            zctx = ctx_get_existing_zctx (csq_con.items[i].dict_id);

        ASSERT (zctx, "Failed to get zctx for INFO/CTX subfield %.*s", name_lens[i], names[i]);
        zctx->st_did_i = is_vep ? INFO_vep : INFO_CSQ;

        #define CB(dnum,cb) if (zctx->dict_id.num == (dnum)) csq_cbs[i] = (cb); else
        CB (_INFO_Allele,             vcf_seg_INFO_allele  )
        CB (_INFO_Existing_variation, seg_id_field_cb      )
        CB (_INFO_cDNA_position,      seg_integer_or_not_cb)
        CB (_INFO_CDS_position,       seg_integer_or_not_cb)
        CB (_INFO_Protein_position,   seg_integer_or_not_cb)
        CB (_INFO_DISTANCE,           seg_integer_or_not_cb)
        /* else */                        {};
        #undef CB
    }

// 5-7 together as Feature implies Feature_Type and Transcript_BioType
// 39-55 all 17 AF fields together, as their sequence in b250 is repetative
}

void vcf_vep_seg_initialize (VBlockVCFP vb)
{
    Did did_i;
    #define IF_HAS(dnum,action) if ((did_i = ctx_get_existing_did_i (vb, (DictId)(dnum))) != DID_NONE) action
    IF_HAS (_INFO_cDNA_position     , CTX(did_i)->ltype = LT_DYN_INT); 
    IF_HAS (_INFO_CDS_position      , CTX(did_i)->ltype = LT_DYN_INT); 
    IF_HAS (_INFO_Protein_position  , CTX(did_i)->ltype = LT_DYN_INT); 
    IF_HAS (_INFO_DISTANCE          , CTX(did_i)->ltype = LT_DYN_INT); 
    IF_HAS (_INFO_Existing_variation, seg_id_field_init (CTX(did_i)));
    #undef IF_HAS
}

// CSQ subfields vary according to VEP parameters. Two examples:
// ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral">
// ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|AA_MAF|EA_MAF|ALLELE_NUM|EXON|INTRON|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|DISTANCE|STRAND|CLIN_SIG|CANONICAL|SYMBOL|SYMBOL_SOURCE|SIFT|PolyPhen|GMAF|BIOTYPE|ENSP|DOMAINS|CCDS|HGVSc|HGVSp|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|PUBMED|LoF_info|LoF_flags|LoF_filter|LoF">

// example: CSQ=-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000438504|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|YES||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|6/6||ENST00000450305.2:n.448_449delGC||448-449|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.734_735delGC||734-735|||||rs780379327|1||1||deletion|1|HGNC|37102|YES||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs780379327|1|917|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000515242|transcribed_unprocessed_pseudogene|3/3||ENST00000515242.2:n.727_728delGC||727-728|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000518655|transcribed_unprocessed_pseudogene|3/4||ENST00000518655.2:n.565_566delGC||565-566|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000538476|unprocessed_pseudogene||||||||||rs780379327|1|924|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000541675|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001576075|CTCF_binding_site||||||||||rs780379327|1||||deletion|1|||||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|
void vcf_seg_INFO_CSQ (VBlockVCFP vb, ContextP ctx, STRp(csq))
{
    seg_array_of_struct (VB, ctx, csq_con, STRa(csq), csq_cbs, csq_len);
}
