// ------------------------------------------------------------------
//   vcf_seg.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "vcf_private.h"
#include "random_access.h"
#include "zip.h"
#include "chrom.h"
#include "stats.h"
#include "tip.h"
#include "arch.h"

// called by main thread after reading the header
void vcf_zip_initialize (void)
{
    vcf_refalt_zip_initialize();
    vcf_samples_zip_initialize();
    vcf_info_zip_initialize();
    vcf_dbsnp_zip_initialize(); // called even if not in VCF header, because can be discovered in segconf too
    vcf_gatk_zip_initialize();
    vcf_gnomad_zip_initialize();
    vcf_samples_zip_initialize_PS_PID();
    vcf_gwas_zip_initialize(); 
    vcf_me_zip_initialize();
    vcf_giab_zip_initialize();
    vcf_1000G_zip_initialize();

    if (segconf.vcf_is_vagrent)    vcf_vagrent_zip_initialize();
    if (segconf.vcf_is_mastermind) vcf_mastermind_zip_initialize();
    if (segconf.vcf_is_vep)        vcf_vep_zip_initialize();
    if (segconf.vcf_illum_gtyping) vcf_illum_gtyping_zip_initialize();
    if (segconf.vcf_is_platypus)   vcf_platypus_zip_initialize();
    if (segconf.vcf_is_freebayes)  vcf_freebayes_zip_initialize();
    if (segconf.vcf_is_svaba)      vcf_svaba_zip_initialize();
    if (segconf.vcf_is_manta)      vcf_manta_zip_initialize();
    if (segconf.vcf_is_pbsv)       vcf_pbsv_zip_initialize();
    if (!segconf.vcf_is_sv)        vcf_sv_zip_initialize (0, 0); // if not any known SV caller (manta, pbsv, svaba), we still initialize for the non-specific stuff

    // set header_info for QUAL, consumed by get_max_items
    ZCTX(VCF_QUAL)->header_info.vcf.Type   = VCF_Float;
    ZCTX(VCF_QUAL)->header_info.vcf.Number = 1;
}

// called after each file
void vcf_zip_finalize (bool is_last_user_txt_file)
{
    // if REFALT takes more than 10% of z_file, advise on using --reference (note: if some cases, we already advised in vcf_segconf_finalize)
    decl_zctx(VCF_REFALT);
    int refalt_z_pc = flag.zip_no_z_file ? 0 : (100 * (zctx->dict.count + zctx->b250.count + zctx->local.count) / z_file->disk_size);

    if (!flag.reference && refalt_z_pc > 10 && !flag.seg_only)
        TIP ("Compressing a this %s file using a reference file can reduce the compressed file's size by %d%%-%d%%.\n"
             "Use: \"%s --reference 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 %s\". 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 may be a FASTA file or a .ref.genozip file.\n",
             z_dt_name(), refalt_z_pc / 3, (int)((float)refalt_z_pc / 1.5), arch_get_argv0(), txt_file->name);

    if (!flag.let_OS_cleanup_on_exit) {
        if (IS_REF_EXT_STORE)
            ref_destroy_reference();
    }
}

// detect if a generic file is actually a vcf
bool is_vcf (STRp(header), bool *need_more)
{
    if (header_len < STRLEN(VCF_MAGIC)) {
        *need_more = true;
        return false;
    }

    return str_isprefix_(STRa(header), _S(VCF_MAGIC));
}

// detect if a generic file is actually a bcf
bool is_bcf (STRp(header), bool *need_more)
{
    if (header_len < STRLEN(BCF_MAGIC)) {
        if (need_more) *need_more = true;
        return false;
    }

    return str_isprefix_(STRa(header), _S(BCF_MAGIC));
}


// main thread: writing data-type specific fields to genozip header
void vcf_zip_genozip_header (SectionHeaderGenozipHeaderP header)
{
    header->vcf.segconf_has_RGQ          = (segconf.has[FORMAT_RGQ] > 0); // introduced in v14
    header->vcf.segconf_GQ_method        = segconf.FMT_GQ_method;         // since 15.0.37
    header->vcf.segconf_INF_DP_method    = segconf.INFO_DP_method;        // since 15.0.52
    header->vcf.segconf_FMT_DP_method    = segconf.FMT_DP_method;         // since 15.0.37
    header->vcf.max_ploidy_for_mux       = ZIP_MAX_PLOIDY_FOR_MUX;        // since 15.0.36
    header->vcf.segconf_MATEID_method    = segconf.MATEID_method;         // since 15.0.48
    header->vcf.segconf_del_svlen_is_neg = segconf.vcf_del_svlen_is_neg;  // since 15.0.48
    header->vcf.segconf_sample_copy      = segconf.vcf_sample_copy;       // since 15.0.69
    header->vcf.segconf_Q_to_O           = BGEN32F (segconf.Q_to_O);      // since 15.0.61
    header->vcf.width.MLEAC              = segconf.wid[INFO_MLEAC].width; // since 15.0.37
    header->vcf.width.AC                 = segconf.wid[INFO_AC].width;    // since 15.0.37
    header->vcf.width.AF                 = segconf.wid[INFO_AF].width;    // since 15.0.37
    header->vcf.width.AN                 = segconf.wid[INFO_AN].width;    // since 15.0.37
    header->vcf.width.DP                 = segconf.wid[INFO_DP].width;    // since 15.0.37
    header->vcf.width.SF                 = segconf.wid[INFO_SF].width;    // since 15.0.37
    header->vcf.width.QD                 = segconf.wid[INFO_QD].width;    // since 15.0.37
    header->vcf.width.AS_SB_TABLE        = segconf.wid[INFO_AS_SB_TABLE].width; // since 15.0.41
    header->vcf.width.ID                 = segconf.wid[VCF_ID].width;     // since 15.0.48
    header->vcf.width.QUAL               = segconf.wid[VCF_QUAL].width;   // since 15.0.51
    header->vcf.width.BaseCounts         = segconf.wid[INFO_BaseCounts].width;  // since 15.0.52
    header->vcf.width.DPB                = segconf.wid[INFO_DPB].width;   // since 15.0.61
}

// ZIP main thread
void vcf_zip_init_vb (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    // set vcf_version in VB, since the the vcf_version in vcf_header might change as we might be reading the next txt file
    vb->vcf_version = vcf_header_get_version();

    // in case we're replacing ID with the line number
    if (flag.zip_lines_counted_at_init_vb) 
        vcf_zip_add_line_numbers_init_vb (vb);
}

// called by main thread, as VBs complete (might be out-of-order)
void vcf_zip_after_compute (VBlockP vb)
{
    z_file->max_ploidy = MAX_(z_file->max_ploidy, VB_VCF->ploidy);

    z_file->vcf_num_samples_copied += CTX(VCF_COPY_SAMPLE)->num_samples_copied;
}   


void vcf_zip_set_txt_header_flags (struct FlagsTxtHeader *f)
{
}

void vcf_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeaderP vb_header)
{
    // note: we write 0xffffffff instead of 0, to differentiate from the case where we didn't write at all prior to 15.0.48
    decl_ctx (FORMAT_GT_HT); 
    vb_header->vcf_HT_n_lines = 
        ctx->HT_n_lines ? BGEN32 (ctx->HT_n_lines) : 0xffffffff; // since 15.0.48
}

// called by Compute threadfrom seg_all_data_lines
void vcf_seg_initialize (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    ctx_consolidate_stats (VB, VCF_ID, VCF_MATE, DID_EOL);
    ctx_consolidate_stats (VB, FORMAT_GT, FORMAT_GT_HT, FORMAT_GT_HT_BIG, FORMAT_PBWT_RUNS, FORMAT_PBWT_FGRC, DID_EOL);

    ctx_set_no_stons (VB, VCF_CHROM, VCF_FORMAT, VCF_INFO, VCF_TOPLEVEL, 
                      VCF_POS, VCF_LINE_NUM, INFO_HGVS_del_start_pos, INFO_HGVS_ins_start_pos, INFO_HGVS_ins_start_pos, // as required by seg_pos_field
                      DID_EOL);

    ctx_set_store (VB, STORE_INDEX, VCF_CHROM, VCF_FORMAT, FORMAT_GT/*=stores ploidy*/, DID_EOL);

    ctx_set_store (VB, STORE_INT, VCF_POS, VCF_LINE_NUM, 
                   T(!segconf.vcf_is_sv, VCF_ID),    // svaba/manta needs to store history as text and segs ID differently
                   T(segconf.vcf_is_sv, VCF_QUAL),
                   T(segconf.vcf_sample_copy, VCF_COPY_SAMPLE),
                   FORMAT_DP, FORMAT_MIN_DP, 
                   INFO_DP,
                   DID_EOL); 

    ctx_set_store (VB, STORE_FLOAT, T(segconf.has[INFO_QD], VCF_QUAL), DID_EOL); // consumed by vcf_piz_insert_INFO_QD

    ctx_set_ltype (VB, LT_STRING, 
                   T(segconf.vcf_QUAL_method == VCF_QUAL_local, VCF_QUAL), 
                   T(segconf.vcf_is_mastermind, INFO_MMURI), 
                   T(segconf.vcf_is_dbSNP, INFO_FREQ),
                   INFO_FATHMM_score, INFO_VEST3_score, INFO_MEINFO, DID_EOL);

    ctx_set_dyn_int (VB, T(segconf.vcf_is_dbSNP, INFO_dbSNPBuildID), 
                     T(segconf.INFO_DP_method == BY_BaseCounts, INFO_DP), // we store deltas in local. also, we can't have stons bc piz relies on ctx->last_wi.
                     DID_EOL);
    
    CTX(FORMAT_DP)->flags.same_line = true; // delta against AD or SDP regardless if before or after on line
    CTX(FORMAT_AD)->flags.same_line = true; // delta against ADALL regardless if before or after on line
    
    // counts sections
    CTX(VCF_CHROM)-> counts_section = true;

    buf_alloc_exact (vb, CTX(VCF_FORMAT)->sf_i, MAX_DICTS, uint16_t, "contexts->sf_i");
    
    // room for already existing FORMATs from previous VBs
    ContextP samples_ctx = CTX(VCF_SAMPLES);
    uint32_t n_fmts = CTX(VCF_FORMAT)->ol_nodes.len;
    buf_alloc_exact_zero (vb, samples_ctx->format_mapper_buf, n_fmts, Container, "contexts->format_mapper_buf");
    buf_alloc_exact_zero (vb, samples_ctx->format_contexts, n_fmts, ContextPBlock, "contexts->format_contexts");

    if (segconf.vcf_QUAL_method == VCF_QUAL_by_RGQ) {
        seg_mux_init (vb, VCF_QUAL, VCF_SPECIAL_MUX_BY_HAS_RGQ, false, QUAL);

        if (vcf_num_samples > 1) // too many unique QUAL values when there are many samples
            ctx_get_ctx (vb, vb->mux_QUAL.dict_ids[0])->ltype = LT_STRING;
    }

    if (segconf.vcf_INFO_method == VCF_INFO_by_RGQ) 
        seg_mux_init (vb, VCF_INFO, VCF_SPECIAL_MUX_BY_HAS_RGQ, false, INFO);        
    
    else if (segconf.vcf_INFO_method == VCF_INFO_by_FILTER)
        seg_mux_init (vb, VCF_INFO, VCF_SPECIAL_MUX_BY_ISAAC_FILTER, false, INFO);        

    if (segconf.vcf_is_gvcf)
        seg_mux_init (vb, VCF_POS, VCF_SPECIAL_MUX_BY_END, false, POS);

    vcf_refalt_seg_initialize (vb);
    vcf_info_seg_initialize(vb);
    vcf_samples_seg_initialize(vb);
    vcf_gatk_seg_initialize (vb);
    vcf_gnomad_seg_initialize (vb);
    vcf_me_seg_initialize (vb);
    vcf_1000G_seg_initialize (vb);
    vcf_copy_sample_seg_initialize (vb);
    if (flag.add_line_numbers)      vcf_add_line_numbers_seg_initialize (vb);
    if (segconf.vcf_illum_gtyping)  vcf_illum_gtyping_seg_initialize (vb);
    if (segconf.vcf_is_gwas)        vcf_gwas_seg_initialize (vb);
    if (segconf.vcf_is_cosmic)      vcf_cosmic_seg_initialize (vb);
    if (segconf.vcf_is_vep)         vcf_vep_seg_initialize (vb);
    if (segconf.vcf_is_mastermind)  vcf_mastermind_seg_initialize (vb);
    if (segconf.vcf_is_dbSNP)       vcf_dbsnp_seg_initialize (vb);
    if (segconf.vcf_is_giab_trio || segconf.vcf_is_giab) vcf_giab_seg_initialize (vb);
    if (segconf.vcf_is_isaac)       vcf_isaac_seg_initialize (vb);
    if (segconf.vcf_is_freebayes)   vcf_freebayes_seg_initialize (vb);
    if (segconf.vcf_is_ultima)      vcf_ultima_seg_initialize (vb);
    if (segconf.vcf_is_platypus)    vcf_platypus_seg_initialize (vb);
    if (segconf.vcf_is_svaba)       vcf_svaba_seg_initialize (vb);
    if (segconf.vcf_is_manta)       vcf_manta_seg_initialize (vb);
    if (segconf.vcf_is_pbsv)        vcf_pbsv_seg_initialize (vb);
    if (!segconf.vcf_is_sv)         vcf_sv_seg_initialize (vb, 0, 0); // not known SV caller - initialize generic stuff
}             

void vcf_segconf_finalize (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    if (!segconf.vcf_evidence_not_gvcf || // all POS values in segconf are consecutive
        (segconf.has[FORMAT_ICNT] && segconf.has[FORMAT_SPL])) // DRAGEN GVCF - variants are consolidated using INFO/END so POS are not consecutive
        segconf.vcf_is_gvcf = true;

    vcf_segconf_finalize_QUAL (vb);
    vcf_copy_samples_segconf_finalize (vb);

    if (segconf.has[FORMAT_RGQ] || segconf.vcf_is_gatk_gvcf)
        segconf.vcf_INFO_method = VCF_INFO_by_RGQ;

    else if (segconf.vcf_is_isaac) 
        segconf.vcf_INFO_method = VCF_INFO_by_FILTER;
            
    // GATK GVCF: set fields as if they were encountered, as often they are encountered starting deep in the file
    if (segconf.vcf_is_gatk_gvcf) {
        Did gvcf_dids[] = { FORMAT_DP, FORMAT_RGQ, FORMAT_GT, FORMAT_PL, FORMAT_AD, FORMAT_GQ };
        for (int i=0; i < ARRAY_LEN(gvcf_dids); i++)
            segconf.has[gvcf_dids[i]] = true;

        segconf.FMT_GQ_method = BY_PL;
        segconf.FMT_DP_method = BY_AD; // override potential setting to BY_AD in segconf, bc most lines don't have AD
        segconf.has_DP_before_PL = true;
        if (flag.best) segconf.PL_mux_by_DP = yes;
    }

    if (segconf.has[INFO_AS_SB_TABLE] && segconf.has[FORMAT_SB]) 
        segconf.AS_SB_TABLE_by_SB = true;
    
    if (segconf.has[INFO_SNVHPOL] || (segconf.has[FORMAT_DPF] && segconf.has[FORMAT_GQX]))
        segconf.vcf_is_isaac = true; // other Illumina tools that generate isaac-like annotations
        
    if (segconf.has[INFO_DP]) {
        if (segconf.has[INFO_BaseCounts] && !flag.secure_DP)
            segconf.INFO_DP_method = BY_BaseCounts;
            
        else if (segconf.has[FORMAT_DP] && segconf.FMT_DP_method != BY_INFO_DP && !flag.secure_DP)
            segconf.INFO_DP_method = BY_FORMAT_DP;

        else 
            segconf.INFO_DP_method = INFO_DP_DEFAULT; // note: INFO_DP_DEFAULT≠0, so set explicitly 
    }

    if (segconf.has[FORMAT_DP] && segconf.has[FORMAT_FI])
        segconf.FI_by_DP = true;

    if (segconf.has[FORMAT_IGT] && segconf.has[FORMAT_IPS] && segconf.has[FORMAT_ADALL])
        segconf.vcf_is_giab_trio = true;
    
    if (segconf.has[INFO_callsetnames] && segconf.has[FORMAT_ADALL])
        segconf.vcf_is_giab = true;

    // an alternative way to discover dbSNP if ##source is missing
    if (!!segconf.has[INFO_RS] + !!segconf.has[INFO_RSPOS] + !!segconf.has[INFO_SAO] + !!segconf.has[INFO_SSR] +
        !!segconf.has[INFO_dbSNPBuildID] + !!segconf.has[INFO_VC] + !!segconf.has[INFO_NSM] + !!segconf.has[INFO_U3] >= 3)
        segconf.vcf_is_dbSNP = true;

    if (segconf.has[INFO_X_LM] && segconf.has[INFO_X_RM]) {
        segconf.vcf_is_ultima = true;

        if (!flag.reference && !flag.seg_only)
            TIP ("Compressing a Ultima Genomics %s file using a reference file can reduce the size by 12%%.\n"
                 "Use: \"%s --reference 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 %s\". 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 may be a FASTA file or a .ref.genozip file.\n",
                 z_dt_name(), arch_get_argv0(), txt_file->name);
    }

    if (segconf.has[INFO_X_HIL] && segconf.has[INFO_X_HIN]) 
        segconf.vcf_is_ultima = true; // another way to identify Ultima

    if (segconf.has[FORMAT_LAA] && segconf.has[FORMAT_LAD] && !segconf.has[FORMAT_AD] && !segconf.has[FORMAT_PL])
        segconf.vcf_local_alleles = true;

    // in gnomAD, we have a huge number of INFO fields in various permutations - generating a huge INFO dictionary, but which compresses very very well
    if (segconf.vcf_is_gnomad)
        ZCTX(VCF_INFO)->dict_len_excessive = true; // don't warn if excessive
        
    if (!flag.reference && segconf.vcf_is_gvcf && !flag.seg_only && vcf_num_samples < 5)
        TIP ("Compressing a GVCF file using a reference file can reduce the compressed file's size by 10%%-30%%.\n"
             "Use: \"%s --reference 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 %s\". 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 may be a FASTA file or a .ref.genozip file.\n",
             arch_get_argv0(), txt_file->name);

    else if (!flag.reference && segconf.vcf_is_platypus && (segconf.has[INFO_SC] || segconf.has[INFO_HP]) && !flag.seg_only)
        TIP ("Compressing a Platypus %s file using a reference file can reduce the compressed file's size by 30%%.\n"
             "Use: \"%s --reference 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 %s\". 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 may be a FASTA file or a .ref.genozip file.\n",
             z_dt_name(), arch_get_argv0(), txt_file->name);

    else if (!flag.reference && segconf.vcf_is_sv && !flag.seg_only)
        TIP ("Compressing a structrual variants %s file using a reference file can reduce the compressed file's size by 20%%-60%%.\n"
            "Use: \"%s --reference 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 %s\". 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 may be a FASTA file or a .ref.genozip file.\n", 
            z_dt_name(), arch_get_argv0(), txt_file->name);

    else if (segconf.has[INFO_SNVHPOL]) // isaac field
        TIP ("Compressing this particular %s using --reference could result in better compression", z_dt_name());

    // In case of dependency DAG: DP->(sum)AD->(mux)GT we can't have GT->(null)DP
    if (segconf.FMT_DP_method == BY_AD) segconf.use_null_DP_method = false;

    // whether we should seg GQ as a function of GP or PL (the better of the two) - only if this works for at least 20% of the samples
    if (segconf.has[FORMAT_GQ] && !segconf.FMT_GQ_method) 
        vcf_segconf_finalize_GQ (vb);

    if (segconf.has_DP_before_PL && !flag.best)
        TIP ("Compressing this particular %s with --best could result in significantly better compression", z_dt_name());

    if (segconf.vcf_is_freebayes) {
        if      (segconf.has[FORMAT_AD]) segconf.FMT_RO_AO_method = RO_AO_by_AD;
        else if (segconf.has[FORMAT_DP]) segconf.FMT_RO_AO_method = RO_AO_by_DP;

        segconf.Q_to_O = segconf.Q_to_O / (float)segconf.n_Q_to_O; 
    }

    // IMPORTANT: number of bits in each field MUST be the same as in SectionHeaderGenozipHeader
    struct { Did did_i; uint8_t width; } wid[] = {
        {VCF_ID,6}, {VCF_QUAL,4}, {INFO_AC,3}, {INFO_MLEAC,3}, {INFO_AN,3}, {INFO_AF,3}, {INFO_SF,3},
        {INFO_QD,3}, {INFO_DP,3}, {INFO_AS_SB_TABLE, 4}, {INFO_BaseCounts,5}, {INFO_DPB,3} 
    };

    for (int i=0; i < ARRAY_LEN(wid); i++)
        segconf_set_width (&segconf.wid[wid[i].did_i], wid[i].width);   

    if (flag.optimize) 
        vcf_segconf_finalize_optimizations (vb);

    // decide which Floats should be segged to local (note: child ctxs inherit the header_info from their parent during consolidate)
    for_ctx_that (ctx->header_info.vcf.Type == VCF_Float && !ctx->is_stats_parent) 
        if (ctx->nodes.len32 > ctx->b250.count / 2) {
            ContextP parent_ctx = ctx->st_did_i == DID_NONE ? ctx : CTX(ctx->st_did_i);
            if (parent_ctx->did_i < z_file->num_contexts)
                ZCTX(parent_ctx->did_i)->seg_to_local = STORE_FLOAT; // note: zctx expected to exist, as it was created from the header
#ifdef DEBUG
            else 
                WARN ("parent_ctx=%s doesn't exist in zctx", parent_ctx->tag_name);
#endif
        }
}
void vcf_seg_finalize (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
        
    // Toplevel snip for reconstructing this VBp
    SmallContainer top_level = { 
        .repeats      = vb->lines.len32,
        .is_toplevel  = true,
        .callback     = (CTX(INFO_SF)->sf.SF_by_GT == yes),   // cases where we need a callback
        .filter_items = true,
        .nitems_lo    = 10,                                                                 
        .items        = { { .dict_id = { _VCF_CHROM },   .separator = "\t" },
                          { .dict_id = { _VCF_POS },     .separator = "\t" },
                          { .dict_id = { _VCF_ID },      .separator = "\t" },
                          { .dict_id = { _VCF_REFALT },  .separator = { '\t', CI1_ITEM_CB } }, // piz calls vcf_piz_refalt_parse
                          { .dict_id = { _VCF_QUAL },    .separator = "\t" },
                          { .dict_id = { _VCF_FILTER },  .separator = "\t" },
                          { .dict_id = { _VCF_INFO },    .separator = "\t" }, 
                          { .dict_id = { _VCF_FORMAT },  .separator = "\t" },
                          { .dict_id = { _VCF_SAMPLES }, .separator = ""   },
                          { .dict_id = { _VCF_EOL },     .separator = ""   } },
    };

    ContextP ctx = CTX(VCF_TOPLEVEL);

    container_seg (vb_, ctx, (ContainerP)&top_level, 0, 0, 0); 
    
    vcf_samples_seg_finalize (vb);

    if (segconf.vcf_sample_copy) vcf_copy_sample_seg_finalize (vb);
    
    __atomic_add_fetch (&z_file->mate_line_count, (uint64_t)vb->mate_line_count,  __ATOMIC_RELAXED);
}

// after each VB is compressed and merge (VB order is arbitrary)
void vcf_zip_after_compress (VBlockP vb)
{
    if (VB_VCF->PL_mux_by_DP == unknown) 
        vcf_FORMAT_PL_decide (VB_VCF);
}

// called after all VBs are compressed - before Global sections are compressed
void vcf_zip_after_vbs (void)
{
    vcf_FORMAT_PL_after_vbs (segconf.vcf_local_alleles ? FORMAT_LPL : FORMAT_PL);
}

bool vcf_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    switch (dict_id.num) {
        case _VCF_TOPLEVEL: case _VCF_CHROM: case _VCF_FORMAT: case _VCF_FILTER: case _VCF_EOL:
        case _VCF_SAMPLES: case _INFO_DP: case _INFO_LDAF: case _INFO_MQ0:
            return true;

        case _INFO_AC: case _INFO_AF: case _INFO_AN: case _INFO_MLEAC: case _INFO_MLEAF:
            return vcf_num_samples >= 1 && vcf_num_samples < 10;

        case _INFO_AA:
            return !segconf.vcf_is_cosmic;

        case _FORMAT_RO: case _FORMAT_AO:
            return segconf.vcf_is_freebayes;

        default: return false;
    }
}

bool vcf_seg_is_big (ConstVBlockP vb, DictId dict_id, DictId st_dict_id/*dict_id of st_did_i*/)
{
    switch (dict_id.num) {
        case _VCF_REFALT: case _INFO_RAW_MQ:
        case _VCF_QUAL:   // QUAL is sometimes in upper tens of thousands, but hash calculated to 64K
            return true;

        case _INFO_AC: case _INFO_AF: case _INFO_MLEAC: case _INFO_MLEAF:
            return vcf_num_samples > 100;

        default: switch (st_dict_id.num) {
            case _FORMAT_PLn: case _FORMAT_PLy: case _FORMAT_GP:
                return true;
            
            default:
                return false;
        }
    }
}

// ---------------------------------------------------------------
// Seg array of integers or floats, expected to be of length n_alt
// ---------------------------------------------------------------

void vcf_seg_array_of_N_ALTS_numbers (VBlockVCFP vb, ContextP ctx, STRp(value), StoreType type)
{
    str_split (value, value_len, N_ALTS, ',', number, true);

    if (!n_numbers) {
        seg_by_ctx (VB, STRa(value), ctx, value_len); // fallback
        return;
    }

    MiniContainer con = { .repeats           = CON_REPEATS_IS_SPECIAL,
                          .repsep            = ",",
                          .drop_final_repsep = true,
                          .nitems_lo         = 1, 
                          .items             = { { .dict_id = sub_dict_id (ctx->dict_id, '0') } } };

    if (!ctx->is_initialized) { 
        ctx_consolidate_stats_(VB, ctx, (ContainerP)&con);
        ctx->con_rep_special = VCF_SPECIAL_N_ALTS;
        ctx->is_initialized  = true;
        ctx->nothing_char    = '.';
    }

    ContextP item_ctx = ctx_get_ctx (vb, con.items[0].dict_id);

    for (uint32_t i=0; i < n_numbers; i++) 
        if (type == STORE_INT)
            seg_integer_or_not (VB, item_ctx, STRi(number,i), 0); // an integer or '.'
        else // float
            seg_add_to_local_string (VB, item_ctx, STRi(number,i), LOOKUP_SIMPLE, 0);

    container_seg (VB, ctx, (ContainerP)&con, 0, 0, value_len);
}

// used by container_reconstruct to retrieve the number of repeats
SPECIAL_RECONSTRUCTOR (vcf_piz_special_N_ALTS)
{
    new_value->i = N_ALTS;
    return HAS_NEW_VALUE;
}

// used by container_reconstruct to retrieve the number of repeats
SPECIAL_RECONSTRUCTOR (vcf_piz_special_N_ALLELES)
{
    new_value->i = N_ALTS + 1;
    return HAS_NEW_VALUE;
}

static inline bool vcf_refalt_seg_ref_alt_line_has_RGQ (rom str)
{
    int count_tabs = 0;
    while (*str != '\n' && (*str != '\t' || (++count_tabs < 4))) str++;

    if (*str == '\n' || str[2]=='\t' || str[2]=='\n' || str[3]=='\t' || str[3]=='\n') return false; // a line without FORMAT, or FORMAT not long enough for RGQ

    for (str++; str[2] != '\t' && str[2] != '\n'; str++) 
        if (str[0]=='R' && str[1]=='G' && str[2]=='Q' && (str[3]=='\t' || str[3]==':'))
            return true;

    return false; // no RGQ in this FORMAT
}

static void vcf_seg_FILTER (VBlockVCFP vb, STRp(filter))
{
    decl_ctx (VCF_FILTER);

    if (segconf.vcf_is_manta || segconf.vcf_is_pbsv)
        vcf_seg_sv_copy_mate (vb, ctx, STRa(filter), TW_FILTER, TW_FILTER, false, filter_len+1);
    
    else
        // note: better than ;-separated array, bc when there are multiple filters, they are often correlated
        seg_by_ctx (VB, STRa(filter), ctx, filter_len+1);

    set_last_txt (VCF_FILTER, filter);
}

// check if filter contains a specific string
bool vcf_FILTER_has (VBlockVCFP vb, rom test_for/*nul-terminated*/)
{
    STRlast (filter, VCF_FILTER);
    
    SAFE_NULT (filter);
    bool has = strstr (filter, test_for);
    SAFE_RESTORE;

    return has;
}

// true is ID is of the format 1_2704352_AT_A which is CHROM_POS_REF_ALT
static inline bool vcf_seg_ID_is_variant (VBlockVCFP vb, ZipDataLineVCF *dl, STRp(id))
{
    str_split (id, id_len, 4, segconf.vcf_ID_is_variant, item, true);
    
    return n_items == 4 &&
           str_issame_(STRi(item,0), STRa(vb->chrom_name))   &&
           str_issame_(STRi(item,1), STRtxt(dl->tw[TW_POS])) &&
           str_issame_(STRi(item,2), STRa(vb->REF))          &&
           str_issame_(STRi(item,3), STRa(vb->ALT));
}

static inline void vcf_seg_ID (VBlockVCFP vb, ZipDataLineVCF *dl, STRp(id))
{
    decl_ctx(VCF_ID);
    
    if (segconf_running && !flag.add_line_numbers) {
        // segconf: set vcf_ID_is_variant in first line
        if (vb->line_i == 0 && !segconf.vcf_is_pbsv) {
            if (!({ segconf.vcf_ID_is_variant = '_' ; vcf_seg_ID_is_variant (vb, dl, STRa(id)); }) &&
                !({ segconf.vcf_ID_is_variant = ':' ; vcf_seg_ID_is_variant (vb, dl, STRa(id)); })) 
                segconf.vcf_ID_is_variant = 0;
        }

        // cases where we want PIZ to insert ID after REF,ALT
        if (segconf.vcf_is_pbsv || segconf.vcf_ID_is_variant)
            SEGCONF_RECORD_WIDTH (VCF_ID, id_len);
    }

    if (flag.add_line_numbers) 
        vcf_seg_line_number_ID (vb, STRa(id));

    else if (segconf.vcf_is_manta)
        vcf_seg_manta_ID (vb, STRa(id));
    
    else if (segconf.vcf_is_svaba)
        vcf_seg_svaba_ID (vb, STRa(id));

    else if (segconf.vcf_is_pbsv)
        vcf_seg_pbsv_ID (vb, STRa(id));

    else if (segconf.vcf_ID_is_variant && vcf_seg_ID_is_variant (vb, dl, STRa(id)))
        // defer reconstruction of ID to after reconstruction of REF/ALT - called from vcf_piz_refalt_parse
        seg_special1 (VB, VCF_SPECIAL_DEFER, segconf.vcf_ID_is_variant/*'_' or ':'*/, ctx, id_len+1);  // often - all the same

    else if (IS_PERIOD(id))
        seg_by_ctx (VB, STRa(id), ctx, id_len+1); // often - all the same

    else 
        seg_id_field (VB, ctx, STRa(id), false, id_len+1);
    
    seg_set_last_txt (VB, ctx, STRa(id));
}

// segment a VCF line into its fields
rom vcf_seg_txt_line (VBlockP vb_, rom line_start, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    ZipDataLineVCF *dl = DATA_LINE (vb->line_i);

    rom next_field=line_start, field_start;
    int32_t len = remaining_txt_len;
    unsigned field_len=0;
    char separator;

    vcf_reset_line (VB);

    GET_NEXT_ITEM (VCF_CHROM);
    dl->chrom = chrom_seg_by_did_i (VB, VCF_CHROM, STRd(VCF_CHROM),VCF_CHROM_len+1);

    GET_NEXT_ITEM (VCF_POS);
    vcf_seg_pos (vb, dl, STRd(VCF_POS));

    // set the tw here (for mating BND lines) (normally set in vcf_seg_sv_copy_mate)
    DATA_LINE (vb->line_i)->tw[TW_CHROM] = TXTWORD (vb->chrom_name); 
    DATA_LINE (vb->line_i)->tw[TW_POS]   = ((TxtWord){ .index = BNUMtxt (VCF_POS_str), .len = VCF_POS_len });

    GET_NEXT_ITEM (VCF_ID);

    // REF + ALT 
    GET_NEXT_ITEM (VCF_REF);
    GET_NEXT_ITEM (VCF_ALT);

    // save REF and ALT (in primary or luft coordinates) to be used for INFO fields
    vb->REF     = VCF_REF_str;
    vb->ALT     = VCF_ALT_str;
    vb->REF_len = VCF_REF_len;
    vb->ALT_len = VCF_ALT_len;
    vb_parse_ALT (vb);

    // ID: seg after REFALT is parsed as in PBSV is predicted based on ALT
    vcf_seg_ID (vb, dl, STRd(VCF_ID));

    CTX(FORMAT_RGQ)->line_has_RGQ = !segconf_running && segconf.has[FORMAT_RGQ] && vcf_refalt_seg_ref_alt_line_has_RGQ (VCF_ALT_str);

    // note: we treat REF+\t+ALT as a single field because REF and ALT are highly corrected, in the case of SNPs:
    // e.g. GG has a probability of 0 and GC has a higher probability than GA.
    vcf_refalt_seg_REF_ALT (vb, STRd(VCF_REF), STRd(VCF_ALT));
    
    GET_NEXT_ITEM (VCF_QUAL);
    seg_set_last_txt (VB, CTX(VCF_QUAL), field_start, field_len); // consumed by vcf_seg_INFO_QD

    GET_NEXT_ITEM (VCF_FILTER);
    seg_set_last_txt (VB, CTX(VCF_FILTER), field_start, field_len);
        
    // INFO
    if (vcf_num_samples) {
        GET_NEXT_ITEM (VCF_INFO); // to do: currently, we require FORMAT (and hence a \t after INFO) in the case the file has samples but this line doesn't. VCF spec doesn't require FORMAT in this case.
    } else {
        GET_MAYBE_LAST_ITEM (VCF_INFO); // may or may not have a FORMAT field
    }

    vcf_seg_info_subfields (vb, field_start, field_len);

    vcf_seg_FILTER (vb, STRd (VCF_FILTER)); // seg after INFO as it might depends on DP

    bool has_samples = false;
    if (separator != '\n') { // has a FORMAT field

        // FORMAT
        if (vcf_num_samples) {
            GET_MAYBE_LAST_ITEM (VCF_FORMAT); // possibly no samples this line
        } else {
            GET_LAST_ITEM (VCF_FORMAT);
        }

        vcf_seg_FORMAT (vb, dl, field_start, field_len);

        if ((has_samples = (separator != '\n')))          
            next_field = vcf_seg_samples (vb, dl, len, (char*)next_field, has_13); 
        else 
            seg_by_did (VB, NULL, 0, VCF_SAMPLES, 0); // case no samples: WORD_INDEX_MISSING
    }

    // case no format or samples
    else {
        seg_by_did (VB, NULL, 0, VCF_FORMAT, 0); 
        seg_by_did (VB, NULL, 0, VCF_SAMPLES, 0);
    }

    vcf_seg_QUAL (vb, STRd (VCF_QUAL)); // seg after samples as it might depends on FORMAT/GP or INFO/DP

    // segs deferred INFO fields and the INFO container
    vcf_seg_finalize_INFO_fields (vb);

    if (!has_samples && vcf_num_samples)
        WARN_ONCE ("FYI: variant CHROM=%.*s POS=%"PRId64" has no samples", vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS));

    SEG_EOL (VCF_EOL, false);

    return next_field;
}
