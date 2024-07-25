// ------------------------------------------------------------------
//   vcf_vep.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"

// Handle VEP fields: https://www.ensembl.org/info/docs/tools/vep/index.html

static MediumContainer csq_con;
static SegCallback csq_cbs[MEDIUM_CON_NITEMS];

static uint64_t _INFO_Allele, _INFO_Existing_variation, _INFO_cDNA_position, _INFO_CDS_position, _INFO_Protein_position, _INFO_DISTANCE, 
                _INFO_AFR_MAF, _INFO_AMR_MAF, _INFO_EAS_MAF, _INFO_EUR_MAF, _INFO_SAS_MAF, _INFO_AA_MAF, _INFO_EA_MAF, _INFO_ExAC_MAF;

static DictId maf_dict_id[2];
sSTRl(maf_container_snip,100);
           
static bool vcf_vep_Existing_var_cb (VBlockP vb, ContextP ctx, STRp(ev), uint32_t repeat)
{
    seg_array_by_callback (VB, ctx, STRa(ev), '&', seg_id_field_varlen_int_cb, 0, 0, ev_len);
    return true;
}

static bool vcf_vep_af_item (VBlockP vb_, ContextP ctx, STRp(af_item), uint32_t repeat)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    str_split (af_item, af_item_len, 2, ':', item, true);
    
    bool is_snp = vb->REF_len == 1 && vb->alt_lens[repeat] == 1;
    
    if (n_items == 2 && item_lens[0] == 1 &&
          ((is_snp  && *af_item == *vb->alts[repeat]) || 
           (!is_snp && *af_item == '-' ))) {

        seg_by_dict_id (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_next_ALT }, 2, maf_dict_id[0], 0);
        
        seg_by_dict_id (VB, STRi(item,1), maf_dict_id[1], 0);
        
        seg_by_ctx (VB, STRa(maf_container_snip), ctx, af_item_len);
    }

    else 
        seg_by_ctx (VB, STRa(af_item), ctx, af_item_len);

    return true;
}

static bool vcf_vep_af_field (VBlockP vb, ContextP ctx, STRp(af_arr), uint32_t repeat)
{
    seg_array_by_callback (VB, ctx, STRa(af_arr), '&', vcf_vep_af_item, VCF_SPECIAL_N_ALTS, VB_VCF->n_alts, af_arr_len);
    return true;
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_next_ALT)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    int8_t alt_i = vb->con_stack[vb->con_stack_len-2].repeat;

    if (reconstruct) {
        if (vb->REF_len == 1 && vb->alt_lens[alt_i] == 1) // SNP
            RECONSTRUCT1 (*vb->alts[alt_i]);
        else
            RECONSTRUCT1('-');
    }

    return NO_NEW_VALUE;
}

// called from vcf_inspect_txt_header_zip, NOT vcf_info_zip_initialize
void vcf_vep_zip_initialize (void) // nul-terminated string containing list of fields. list of fields expected to be terminated by a '"'
{
    ASSERTNOTNULL (segconf.vcf_vep_spec); // allocated in vcf_inspect_txt_header_zip
    
    _INFO_Allele             = dict_id_make (_S("Allele"),             DTYPE_VCF_INFO).num;
    _INFO_Existing_variation = dict_id_make (_S("Existing_variation"), DTYPE_VCF_INFO).num;
    _INFO_cDNA_position      = dict_id_make (_S("cDNA_position"),      DTYPE_VCF_INFO).num;
    _INFO_CDS_position       = dict_id_make (_S("CDS_position"),       DTYPE_VCF_INFO).num;
    _INFO_Protein_position   = dict_id_make (_S("Protein_position"),   DTYPE_VCF_INFO).num;
    _INFO_DISTANCE           = dict_id_make (_S("DISTANCE"),           DTYPE_VCF_INFO).num;
    _INFO_AFR_MAF            = dict_id_make (_S("AFR_MAF"),            DTYPE_VCF_INFO).num;
    _INFO_AMR_MAF            = dict_id_make (_S("AMR_MAF"),            DTYPE_VCF_INFO).num;
    _INFO_EAS_MAF            = dict_id_make (_S("EAS_MAF"),            DTYPE_VCF_INFO).num;
    _INFO_EUR_MAF            = dict_id_make (_S("EUR_MAF"),            DTYPE_VCF_INFO).num;
    _INFO_SAS_MAF            = dict_id_make (_S("SAS_MAF"),            DTYPE_VCF_INFO).num;
    _INFO_AA_MAF             = dict_id_make (_S("AA_MAF"),             DTYPE_VCF_INFO).num;
    _INFO_EA_MAF             = dict_id_make (_S("EA_MAF"),             DTYPE_VCF_INFO).num;
    _INFO_ExAC_MAF           = dict_id_make (_S("ExAC_MAF"),           DTYPE_VCF_INFO).num; // note: this includes all ExAC_*_MAF fields that get mapped to the same context due to dict_id_make logic

    rom spec = segconf.vcf_vep_spec;
    int spec_len = strlen (segconf.vcf_vep_spec);

    str_split (spec, spec_len, 0, '|', name, false);

    if (n_names > MEDIUM_CON_NITEMS) { // if this every happens, we should switch to a larger container
        WARN ("FYI: VEP CSQ field has %u annotations, but compression will be sub-optimal as it has more than %u fields. - please report to " EMAIL_SUPPORT ":\n%.*s", n_names, CONTAINER_MAX_DICTS, STRf(spec));    
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
    for (rom c=spec + spec_len -1; c >= spec; c--)
        if (*c == '#') {
            is_vep = !!strstr (c, "ID=vep");
            break;
        }

    for (unsigned i=0; i < n_names; i++) {
        // special case: we store AF in VEP_AF, to not conflict with INFO_AF (we don't reconstruct VEP field names, this will only affect STATS output)
        // MAX_AF is another overlap item between vcf.h and the VEP field list https://ensembl.org/info/docs/tools/vep/vep_formats.html but currently there is no special method for it, so no conflict worries
        if (str_issame_(STRi(name,i), "AF", 2)) {
            names[i] = "VEP_AF"; name_lens[i]=6;
        }

        csq_con.items[i] = (ContainerItem){ .dict_id = dict_id_make (names[i], name_lens[i], DTYPE_VCF_INFO), .separator[0] = ((i < n_names-1) ? '|'  : 0) };

        ContextP zctx = ctx_add_new_zf_ctx_at_init (names[i], name_lens[i], csq_con.items[i].dict_id); // zctx, but will be the same did for vb contexts because pre-created
        if (!zctx) // already exist - predefined
            zctx = ctx_get_existing_zctx (csq_con.items[i].dict_id);

        ASSERT (zctx, "Failed to get zctx for %s subfield \"%.*s\"", 
                ZCTX(is_vep ? INFO_vep : INFO_CSQ)->tag_name, name_lens[i], names[i]);

        zctx->st_did_i = is_vep ? INFO_vep : INFO_CSQ;

        #define CB(dnum,cb) if (zctx->dict_id.num == (dnum)) csq_cbs[i] = (cb); else

        // if and not switch, because DICT_ID_MAKE1_* do not produce an integer constant
        CB (_INFO_Allele,             vcf_seg_INFO_allele    )
        CB (_INFO_Existing_variation, vcf_vep_Existing_var_cb)
        CB (_INFO_cDNA_position,      seg_integer_or_not_cb  )
        CB (_INFO_CDS_position,       seg_integer_or_not_cb  )
        CB (_INFO_Protein_position,   seg_integer_or_not_cb  )
        CB (_INFO_DISTANCE,           seg_integer_or_not_cb  )

        CB (_INFO_AFR_MAF,  vcf_vep_af_field)
        CB (_INFO_AMR_MAF,  vcf_vep_af_field)
        CB (_INFO_EAS_MAF,  vcf_vep_af_field)
        CB (_INFO_EUR_MAF,  vcf_vep_af_field)
        CB (_INFO_SAS_MAF,  vcf_vep_af_field)
        CB (_INFO_AA_MAF,   vcf_vep_af_field)
        CB (_INFO_EA_MAF,   vcf_vep_af_field)
        CB (_INFO_ExAC_MAF, vcf_vep_af_field)
        /* else */          {};

        #undef CB
    }

    maf_dict_id[0] = (DictId)DICT_ID_MAKE1_8("VEP_MAF0");
    maf_dict_id[1] = (DictId)DICT_ID_MAKE1_8("VEP_MAF0");
    
    SmallContainer con = { .repeats   = 1, 
                           .nitems_lo = 2,
                           .items[0]  = { .dict_id = maf_dict_id[0], .separator = ":" },
                           .items[1]  = { .dict_id = maf_dict_id[1]                   } };

    maf_container_snip_len = sizeof (maf_container_snip); // re-initialize for every file
    container_prepare_snip ((ContainerP)&con, 0, 0, qSTRa (maf_container_snip));

    FREE (segconf.vcf_vep_spec); // allocated in vcf_inspect_txt_header_zip
}

void vcf_vep_seg_initialize (VBlockVCFP vb)
{
    Did did_i;
    #define IF_HAS(dnum,action) if ((did_i = ctx_get_existing_did_i (vb, (DictId)(dnum))) != DID_NONE) action
    IF_HAS (_INFO_cDNA_position   , ctx_set_dyn_int (VB, did_i, DID_EOL)); 
    IF_HAS (_INFO_CDS_position    , ctx_set_dyn_int (VB, did_i, DID_EOL)); 
    IF_HAS (_INFO_Protein_position, ctx_set_dyn_int (VB, did_i, DID_EOL)); 
    IF_HAS (_INFO_DISTANCE        , ctx_set_dyn_int (VB, did_i, DID_EOL)); 
    #undef IF_HAS

    ctx_get_ctx (vb, maf_dict_id[0])->st_did_i = ctx_get_ctx (vb, maf_dict_id[1])->st_did_i = INFO_CSQ;
}

// CSQ subfields vary according to VEP parameters. Two examples:
// ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info|context|ancestral">
// ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|AA_MAF|EA_MAF|ALLELE_NUM|EXON|INTRON|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|DISTANCE|STRAND|CLIN_SIG|CANONICAL|SYMBOL|SYMBOL_SOURCE|SIFT|PolyPhen|GMAF|BIOTYPE|ENSP|DOMAINS|CCDS|HGVSc|HGVSp|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|PUBMED|LoF_info|LoF_flags|LoF_filter|LoF">

// example: CSQ=-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000423562|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000438504|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|YES||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene|6/6||ENST00000450305.2:n.448_449delGC||448-449|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.734_735delGC||734-735|||||rs780379327|1||1||deletion|1|HGNC|37102|YES||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs780379327|1|917|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000515242|transcribed_unprocessed_pseudogene|3/3||ENST00000515242.2:n.727_728delGC||727-728|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|non_coding_transcript_exon_variant&non_coding_transcript_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000518655|transcribed_unprocessed_pseudogene|3/4||ENST00000518655.2:n.565_566delGC||565-566|||||rs780379327|1||1||deletion|1|HGNC|37102|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000538476|unprocessed_pseudogene||||||||||rs780379327|1|924|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000541675|unprocessed_pseudogene||||||||||rs780379327|1|876|-1||deletion|1|HGNC|38034|||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001576075|CTCF_binding_site||||||||||rs780379327|1||||deletion|1|||||||||||||||||-:0||||||||-:0|-:1.128e-05|-:0|-:0|-:0|-:0|-:0|-:0|||||||||||||AGCT|
void vcf_seg_INFO_CSQ (VBlockVCFP vb, ContextP ctx, STRp(csq))
{
    seg_array_of_struct (VB, ctx, csq_con, STRa(csq), csq_cbs, NULL, csq_len);
}
