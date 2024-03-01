// ------------------------------------------------------------------
//   vcf_seg.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "random_access.h"
#include "zip.h"
#include "chrom.h"
#include "stats.h"
#include "tip.h"
#include "arch.h"

static MiniContainer line_number_container = {};

// called by main thread after reading the header
void vcf_zip_initialize (void)
{
    vcf_refalt_zip_initialize();
    vcf_samples_zip_initialize ();
    vcf_info_zip_initialize ();

    vcf_dbsnp_zip_initialize(); // called even if not in VCF header, because can be discovered in segconf too
    vcf_gatk_zip_initialize();
    if (segconf.vcf_is_vagrent)    vcf_vagrent_zip_initialize();
    if (segconf.vcf_is_mastermind) vcf_mastermind_zip_initialize();
    if (segconf.vcf_is_vep)        vcf_vep_zip_initialize();
    if (segconf.vcf_illum_gtyping) vcf_illum_gtyping_zip_initialize();
    if (segconf.vcf_is_platypus)   vcf_platypus_zip_initialize();
    if (segconf.vcf_is_svaba)      vcf_svaba_zip_initialize();
    if (segconf.vcf_is_manta)      vcf_manta_zip_initialize();
    if (segconf.vcf_is_pbsv)       vcf_pbsv_zip_initialize();
    if (!segconf.vcf_is_sv)        vcf_sv_zip_initialize (0, 0); // if not any known SV caller (manta, pbsv, svaba), we still initialize for the non-specific stuff

    // container just for adding a prefix to the delta-encoded line number (the container is all_the_same)
    if (flag.add_line_numbers && !line_number_container.repeats) { // possibly already initialized by previous files
        line_number_container = (MiniContainer) {
            .repeats   = 1,
            .nitems_lo = 1,
            .items     = { { .dict_id = { _VCF_LINE_NUM } } }
        };
    }
}

// called after each file
void vcf_zip_finalize (bool is_last_user_txt_file)
{
    // if REFALT takes more than 10% of z_file, advise on using --reference (note: if some cases, we already advised in vcf_seg_finalize_segconf)
    decl_zctx(VCF_REFALT);
    int refalt_z_pc = flag.zip_no_z_file ? 0 : (100 * (zctx->dict.count + zctx->b250.count + zctx->local.count) / z_file->disk_size);

    if (!flag.reference && refalt_z_pc > 10 && !flag.seg_only)
        TIP ("Compressing a this %s file using a reference file can reduce the compressed file's size by %d%%-%d%%.\n"
             "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n",
             z_dt_name(), refalt_z_pc / 3, (int)((float)refalt_z_pc / 1.5), arch_get_argv0(), txt_file->name);
}

// detect if a generic file is actually a vcf
bool is_vcf (STRp(header), bool *need_more)
{
    return header_len >= 16 && !memcmp (header, "##fileformat=VCF", 16);
}

// main thread: writing data-type specific fields to genozip header
void vcf_zip_genozip_header (SectionHeaderGenozipHeaderP header)
{
    header->vcf.segconf_has_RGQ       = (segconf.has[FORMAT_RGQ] > 0); // introduced in v14
    header->vcf.segconf_GQ_method     = segconf.GQ_method;             // since 15.0.37
    header->vcf.segconf_FMT_DP_method = segconf.FMT_DP_method;         // since 15.0.37
    header->vcf.max_ploidy_for_mux    = ZIP_MAX_PLOIDY_FOR_MUX;        // since 15.0.36
    header->vcf.segconf_MATEID_method = segconf.MATEID_method;         // since 15.0.47
    header->vcf.segconf_del_svlen_is_neg = segconf.vcf_del_svlen_is_neg; // since 15.0.47
    header->vcf.width.MLEAC           = segconf.wid_MLEAC.width;       // since 15.0.37
    header->vcf.width.AC              = segconf.wid_AC.width;          // since 15.0.37
    header->vcf.width.AF              = segconf.wid_AF.width;          // since 15.0.37
    header->vcf.width.AN              = segconf.wid_AN.width;          // since 15.0.37
    header->vcf.width.DP              = segconf.wid_DP.width;          // since 15.0.37
    header->vcf.width.SF              = segconf.wid_SF.width;          // since 15.0.37
    header->vcf.width.QD              = segconf.wid_QD.width;          // since 15.0.37
    header->vcf.width.AS_SB_TABLE     = segconf.wid_AS_SB_TABLE.width; // since 15.0.41
    header->vcf.width.ID              = segconf.wid_ID.width;          // since 15.0.47
}

void vcf_zip_init_vb (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    // set vcf_version in VB, since the the vcf_version in vcf_header might change as we might be reading the next txt file
    vb->vcf_version = vcf_header_get_version();

    // in case we're replacing ID with the line number
    if (flag.add_line_numbers) {
        vb->first_line = txt_file->num_lines + 1; 
        txt_file->num_lines += str_count_char (STRb(vb->txt_data), '\n');  // update here instead of in zip_update_txt_counters;
    }
}

// called by main thread, as VBs complete (might be out-of-order)
void vcf_zip_after_compute (VBlockP vb)
{
    z_file->max_ploidy = MAX_(z_file->max_ploidy, VB_VCF->ploidy);
}   


void vcf_zip_set_txt_header_flags (struct FlagsTxtHeader *f)
{
}

void vcf_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeaderP vb_header)
{
    // note: we write 0xffffffff instead of 0, to differentiate from the case where we didn't write at all prior to 15.0.47
    decl_ctx (FORMAT_GT_HT); 
    vb_header->vcf_HT_n_lines = 
        ctx->HT_n_lines ? BGEN32 (ctx->HT_n_lines) : 0xffffffff; // since 15.0.47
}

// called by Compute threadfrom seg_all_data_lines
void vcf_seg_initialize (VBlockP vb_)
{
    #define T(cond, did_i) ((cond) ? (did_i) : DID_NONE)

    VBlockVCFP vb = (VBlockVCFP)vb_;

    ctx_consolidate_stats (VB, VCF_ID, VCF_MATE, DID_EOL);

    ctx_set_no_stons (VB, VCF_CHROM, VCF_FORMAT, VCF_INFO, VCF_TOPLEVEL, 
                      VCF_POS, VCF_LINE_NUM, INFO_HGVS_del_start_pos, INFO_HGVS_ins_start_pos, INFO_HGVS_ins_start_pos, // as required by seg_pos_field
                      DID_EOL);

    ctx_set_store (VB, STORE_INDEX, VCF_CHROM, DID_EOL);

    ctx_set_store (VB, STORE_INT, VCF_POS, VCF_LINE_NUM, 
                   T(!segconf.vcf_is_sv, VCF_ID),    // svaba/manta needs to store history as text and segs ID differently
                   T(segconf.vcf_is_sv, VCF_QUAL),
                   FORMAT_DP, FORMAT_MIN_DP, 
                   INFO_DP, 
                   DID_EOL); 

    ctx_set_store (VB, STORE_FLOAT, T(segconf.has[INFO_QD], VCF_QUAL), DID_EOL); // consumed by vcf_piz_insert_INFO_QD

    ctx_set_ltype (VB, LT_STRING, 
                   T(segconf.vcf_QUAL_method == VCF_QUAL_local, VCF_QUAL), 
                   T(segconf.vcf_is_mastermind, INFO_MMURI), 
                   T(segconf.vcf_is_dbSNP, INFO_FREQ),
                   INFO_FATHMM_score, INFO_VEST3_score, DID_EOL);

    ctx_set_dyn_int (VB, INFO_RAW_MQandDP_MQ, INFO_RAW_MQandDP_DP, 
                     T(segconf.vcf_is_dbSNP, INFO_dbSNPBuildID), 
                     DID_EOL);
    
    CTX(FORMAT_DP)->flags.same_line = true; // delta against AD or SDP regardless if before or after on line
    CTX(FORMAT_AD)->flags.same_line = true; // delta against ADALL regardless if before or after on line
    
    // counts sections
    CTX(VCF_CHROM)-> counts_section = true;

    // room for already existing FORMATs from previous VBs
    ContextP samples_ctx = CTX(VCF_SAMPLES);
    samples_ctx->format_mapper_buf.len = samples_ctx->format_contexts.len = CTX(VCF_FORMAT)->ol_nodes.len;
    buf_alloc_zero (vb, &samples_ctx->format_mapper_buf, 0, samples_ctx->format_mapper_buf.len, Container, CTX_GROWTH, "contexts->format_mapper_buf");
    buf_alloc_zero (vb, &samples_ctx->format_contexts, 0, samples_ctx->format_contexts.len, ContextPBlock, CTX_GROWTH, "contexts->format_contexts");
    
    if (flag.add_line_numbers) {
        // create a b250 and dict entry for VCF_LINE_NUM, VCF_ID - these become "all_the_same" so no need to seg them explicitly hereinafter        
        #define LN_PREFIX_LEN 3 // "LN="
        static const char prefix[] = { CON_PX_SEP, // has prefix 
                                       CON_PX_SEP, // end of (empty) container-wide prefix
                                       'L', 'N', '=', CON_PX_SEP };  // NOTE: if changing prefix, update LN_PREFIX_LEN

        container_seg (vb, CTX(VCF_ID), (Container *)&line_number_container, prefix, sizeof prefix, 0); 
        ctx_decrement_count (vb_, CTX(VCF_ID), 0);
        CTX(VCF_ID)->no_stons = true;
    }
        
    if (segconf.vcf_QUAL_method == VCF_QUAL_by_RGQ) {
        seg_mux_init (VB, CTX(VCF_QUAL), 2, VCF_SPECIAL_MUX_BY_HAS_RGQ, false, (MultiplexerP)&vb->mux_QUAL);

        if (vcf_num_samples > 1) // too many unique QUAL values when there are many samples
            ctx_get_ctx (vb, vb->mux_QUAL.dict_ids[0])->ltype = LT_STRING;
    }

    if (segconf.vcf_INFO_method == VCF_INFO_by_RGQ) 
        seg_mux_init (VB, CTX(VCF_INFO), 2, VCF_SPECIAL_MUX_BY_HAS_RGQ, false, (MultiplexerP)&vb->mux_INFO);        
    
    else if (segconf.vcf_INFO_method == VCF_INFO_by_FILTER)
        seg_mux_init (VB, CTX(VCF_INFO), 2, VCF_SPECIAL_MUX_BY_ISAAC_FILTER, false, (MultiplexerP)&vb->mux_INFO);        

    if (segconf.vcf_is_gvcf)
        seg_mux_init (VB, CTX(VCF_POS), 2, VCF_SPECIAL_MUX_BY_END, false, (MultiplexerP)&vb->mux_POS);

    vcf_refalt_seg_initialize (vb);
    vcf_info_seg_initialize(vb);
    vcf_samples_seg_initialize(vb);

    vcf_gatk_seg_initialize (vb);
    if (segconf.vcf_illum_gtyping)  vcf_illum_gtyping_seg_initialize (vb);
    if (segconf.vcf_is_gwas)        vcf_gwas_seg_initialize (vb);
    if (segconf.vcf_is_cosmic)      vcf_cosmic_seg_initialize (vb);
    if (segconf.vcf_is_vep)         vcf_vep_seg_initialize (vb);
    if (segconf.vcf_is_mastermind)  vcf_mastermind_seg_initialize (vb);
    if (segconf.vcf_is_dbSNP)       vcf_dbsnp_seg_initialize (vb);
    if (segconf.vcf_is_giab_trio)   vcf_giab_seg_initialize (vb);
    if (segconf.vcf_is_isaac)       vcf_isaac_seg_initialize (vb);
    if (segconf.vcf_is_ultima)      vcf_ultima_seg_initialize (vb);
    if (segconf.vcf_is_platypus)    vcf_platypus_seg_initialize (vb);
    if (segconf.vcf_is_svaba)       vcf_svaba_seg_initialize (vb);
    if (segconf.vcf_is_manta)       vcf_manta_seg_initialize (vb);
    if (segconf.vcf_is_pbsv)        vcf_pbsv_seg_initialize (vb);
    if (!segconf.vcf_is_sv)         vcf_sv_seg_initialize (vb, 0, 0); // not known SV caller - initialize generic stuff
    #undef T
}             

static void vcf_segconf_finalize_QUAL (VBlockVCFP vb); // forward

static void vcf_seg_finalize_segconf (VBlockVCFP vb)
{
    vcf_segconf_finalize_QUAL (vb);

    if (segconf.has[FORMAT_RGQ] || segconf.vcf_is_gatk_gvcf)
        segconf.vcf_INFO_method = VCF_INFO_by_RGQ;

    else if (segconf.vcf_is_isaac) 
        segconf.vcf_INFO_method = VCF_INFO_by_FILTER;

    // identify DRAGEN and Isaac GVCF. GATK's is identified in vcf_inspect_txt_header_zip()
    if ((segconf.has[FORMAT_ICNT] && segconf.has[FORMAT_SPL]) || // DRAGEN GVCF
         segconf.vcf_is_isaac) // Isaac is always GVCF
        segconf.vcf_is_gvcf = true;

    // GATK GVCF: set fields as if they were encountered, as often they are encountered starting deep in the file
    if (segconf.vcf_is_gatk_gvcf) {
        Did gvcf_dids[] = { FORMAT_DP, FORMAT_RGQ, FORMAT_GT, FORMAT_PL, FORMAT_AD, FORMAT_GQ };
        for (int i=0; i < ARRAY_LEN(gvcf_dids); i++)
            segconf.has[gvcf_dids[i]] = true;

        segconf.GQ_method = BY_PL;
        segconf.FMT_DP_method = BY_AD; // override potential setting to BY_AD in segconf, bc most lines don't have AD
        segconf.has_DP_before_PL = true;
        if (flag.best) segconf.PL_mux_by_DP = yes;
    }

    if (segconf.has[INFO_AS_SB_TABLE] && segconf.has[FORMAT_SB]) 
        segconf.AS_SB_TABLE_by_SB = true;
    
    if (vcf_num_samples > 1 && !flag.secure_DP && segconf.has[INFO_DP] && segconf.has[FORMAT_DP])
        segconf.INFO_DP_method = BY_FORMAT_DP;
        
    if (segconf.has[FORMAT_IGT] && segconf.has[FORMAT_IPS] && segconf.has[FORMAT_ADALL])
        segconf.vcf_is_giab_trio = true;
        
    // an alternative way to discover dbSNP if ##source is missing
    if (!!segconf.has[INFO_RS] + !!segconf.has[INFO_RSPOS] + !!segconf.has[INFO_SAO] + !!segconf.has[INFO_SSR] +
        !!segconf.has[INFO_dbSNPBuildID] + !!segconf.has[INFO_VC] + !!segconf.has[INFO_NSM] + !!segconf.has[INFO_U3] >= 3)
        segconf.vcf_is_dbSNP = true;

    if (segconf.has[INFO_X_LM] && segconf.has[INFO_X_RM]) {
        segconf.vcf_is_ultima = true;

        if (!flag.reference && !flag.seg_only)
            TIP ("Compressing a Ultima Genomics VCF file using a reference file can reduce the size by 12%%.\n"
                 "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n",
                 arch_get_argv0(), txt_file->name);
    }

    if (segconf.has[INFO_X_HIL] && segconf.has[INFO_X_HIN]) 
        segconf.vcf_is_ultima = true; // another way to identify Ultima

    // in gnomAD, we have a huge number of INFO fields in various permutations - generating a huge INFO dictionary, but which compresses very very well
    if (segconf.vcf_is_gnomad)
        ZCTX(VCF_INFO)->dict_len_excessive = true; // don't warn if excessive
        
    if (!flag.reference && segconf.vcf_is_gvcf && !flag.seg_only)
        TIP ("Compressing a GVCF file using a reference file can reduce the compressed file's size by 10%%-30%%.\n"
             "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n",
             arch_get_argv0(), txt_file->name);

    if (!flag.reference && segconf.vcf_is_platypus && (segconf.has[INFO_SC] || segconf.has[INFO_HP]) && !flag.seg_only)
        TIP ("Compressing a Platypus VCF file using a reference file can reduce the compressed file's size by 30%%.\n"
             "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n",
             arch_get_argv0(), txt_file->name);

    if (!flag.reference && segconf.vcf_is_sv && !flag.seg_only)
        TIP ("Compressing a structrual variants VCF file using a reference file can reduce the compressed file's size by 20%%-60%%.\n"
            "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n", 
            arch_get_argv0(), txt_file->name);

    // In case of dependency DAG: DP->(sum)AD->(mux)GT we can't have GT->(null)DP
    if (segconf.FMT_DP_method == BY_AD) segconf.use_null_DP_method = false;

    // whether we should seg GQ as a function of GP or PL (the better of the two) - only if this works for at least 20% of the samples
    if (segconf.has[FORMAT_GQ] && !segconf.GQ_method) 
        vcf_segconf_finalize_GQ (vb);

    if (segconf.has_DP_before_PL && !flag.best)
        TIP0 ("Compressing this particular VCF with --best could result in significantly better compression");

    segconf_set_width (&segconf.wid_ID, 6); // non-zero only in PBSV
    segconf_set_width (&segconf.wid_AF, 3);
    segconf_set_width (&segconf.wid_AC, 3);
    segconf_set_width (&segconf.wid_AN, 3);
    segconf_set_width (&segconf.wid_DP, 3);
    segconf_set_width (&segconf.wid_QD, 3);
    segconf_set_width (&segconf.wid_SF, 3);
    segconf_set_width (&segconf.wid_MLEAC, 3);
    segconf_set_width (&segconf.wid_AS_SB_TABLE, 4);
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
                          { .dict_id = { _VCF_INFO },    .separator = "\t" }, // in dual-coordinates, contains INFO/LIFTOVER or INFO/REJTOVER that reconstructs oCHROM, oPOS, oREF, oXSTRAND
                          { .dict_id = { _VCF_FORMAT },  .separator = "\t" },
                          { .dict_id = { _VCF_SAMPLES }, .separator = ""   },
                          { .dict_id = { _VCF_EOL },     .separator = ""   } },
    };

    ContextP ctx = CTX(VCF_TOPLEVEL);

    container_seg (vb_, ctx, (ContainerP)&top_level, 0, 0, 0); 
    
    vcf_samples_seg_finalize (vb);

    if (segconf.running)
        vcf_seg_finalize_segconf (vb);

    else
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
    vcf_FORMAT_PL_after_vbs ();
}

bool vcf_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return 
        dict_id.num  == _VCF_TOPLEVEL ||
        dict_id.num  == _VCF_TOPLUFT  ||
        dict_id.num  == _VCF_CHROM    ||
        dict_id.num  == _VCF_oCHROM   ||
        dict_id.num  == _VCF_FORMAT   || // note: NOT including _VCF_INFO - there are cases where it is not small
        dict_id.num  == _VCF_FILTER   || 
        dict_id.num  == _VCF_EOL      ||
        dict_id.num  == _VCF_SAMPLES  ||
        dict_id.num  == _VCF_oCHROM   ||
        dict_id.num  == _VCF_oXSTRAND ||
        dict_id.num  == _VCF_oSTATUS  ||
        dict_id.num  == _VCF_COORDS   ||
        dict_id.num  == _VCF_LIFT_REF ||
        (dict_id.num == _INFO_AC    && vcf_num_samples >= 1 && vcf_num_samples < 10) || // note: _INFO_AN, _INFO_AC, _INFO_AF, _INFO_MLEAC, _INFO_MLEAF - can be big (e.g. big in GWAS VCF, ExAC, gnomad...)
        (dict_id.num == _INFO_AF    && vcf_num_samples >= 1 && vcf_num_samples < 10) ||
        (dict_id.num == _INFO_AN    && vcf_num_samples >= 1 && vcf_num_samples < 10) ||
        (dict_id.num == _INFO_MLEAC && vcf_num_samples >= 1 && vcf_num_samples < 10) ||
        (dict_id.num == _INFO_MLEAF && vcf_num_samples >= 1 && vcf_num_samples < 10) ||
        dict_id.num  == _INFO_DP      ||
        (dict_id.num == _INFO_AA && !segconf.vcf_is_cosmic) || // stored as a SPECIAL snip
        dict_id.num  == _INFO_LDAF    ||
        dict_id.num  == _INFO_MQ0     ||
        dict_id.num  == _INFO_LUFT    ||
        dict_id.num  == _INFO_PRIM    ||
        dict_id.num  == _INFO_LREJ    ||
        dict_id.num  == _INFO_PREJ;       
}

bool vcf_seg_is_big (ConstVBlockP vb, DictId dict_id, DictId st_dict_id/*dict_id of st_did_i*/)
{
    return 
        (dict_id.num   == _INFO_AC    && vcf_num_samples > 100) ||
        (dict_id.num   == _INFO_AF    && vcf_num_samples > 100) ||
        (dict_id.num   == _INFO_MLEAC && vcf_num_samples > 100) ||
        (dict_id.num   == _INFO_MLEAF && vcf_num_samples > 100) ||
        dict_id.num    == _VCF_REFALT  ||
        dict_id.num    == _VCF_oREFALT ||
        dict_id.num    == _VCF_QUAL    || // QUAL is sometimes in upper tens of thousands, but hash calculated to 64K
        dict_id.num    == _INFO_RAW_MQ ||
        st_dict_id.num == _FORMAT_PLn  || // P1Ln1 etc are often big
        st_dict_id.num == _FORMAT_PLy  ||
        st_dict_id.num == _FORMAT_GP    ; // G1P1 etc are often big
}

static void vcf_seg_add_line_number (VBlockVCFP vb, unsigned VCF_ID_len)
{
    char line_num[20];
    unsigned line_num_len = str_int (vb->first_line + vb->line_i, line_num);

    seg_pos_field (VB, VCF_LINE_NUM, VCF_LINE_NUM, SPF_ZERO_IS_BAD, 0, line_num, line_num_len, 0, line_num_len + 1 + LN_PREFIX_LEN); // +1 for tab +3 for "LN="
    
    int shrinkage = (int)VCF_ID_len - line_num_len - LN_PREFIX_LEN;
    vb->recon_size -= shrinkage;
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

static void vcf_segconf_finalize_QUAL (VBlockVCFP vb)
{
    // gnomAD, dbSNP - store in local
    if (segconf.vcf_is_gnomad || segconf.vcf_is_dbSNP)
        segconf.vcf_QUAL_method = VCF_QUAL_local;

    // GVCF - multiplex by has_RGQ
    else if (segconf.has[FORMAT_RGQ] || segconf.vcf_is_gatk_gvcf)
        segconf.vcf_QUAL_method = VCF_QUAL_by_RGQ;
    
    // SvABA - copy from mate, or store in local
    else if (segconf.vcf_is_svaba || segconf.vcf_is_manta) // note: in pbsv, QUAL="."
        segconf.vcf_QUAL_method = VCF_QUAL_mated;

    else
        segconf.vcf_QUAL_method = VCF_QUAL_DEFAULT;
}

static void vcf_seg_QUAL (VBlockVCFP vb, STRp(qual))
{
    decl_ctx (VCF_QUAL);

    switch (segconf.vcf_QUAL_method) {
        case VCF_QUAL_local :
            seg_add_to_local_string (VB, ctx, STRa(qual), LOOKUP_NONE, qual_len+1);
            break;

        case VCF_QUAL_by_RGQ: {
            bool has_rgq = CTX(FORMAT_RGQ)->line_has_RGQ;
            ContextP channel_ctx = seg_mux_get_channel_ctx (VB, VCF_QUAL, (MultiplexerP)&vb->mux_QUAL, has_rgq);
            
            if (!has_rgq && vcf_num_samples > 1) // too many unique QUAL values when there are many samples
                seg_add_to_local_string (VB, channel_ctx, STRa(qual), LOOKUP_SIMPLE, qual_len+1);
            else
                seg_by_ctx (VB, STRa(qual), channel_ctx, qual_len+1);
            
            seg_by_ctx (VB, STRa(vb->mux_QUAL.snip), ctx, 0);
            break;
        }

        case VCF_QUAL_mated: 
            vcf_seg_sv_copy_mate (vb, ctx, STRa(qual), TW_QUAL, TW_QUAL, false, qual_len+1);
            break;

        case VCF_QUAL_DEFAULT:
            seg_by_ctx (VB, STRa(qual), ctx, qual_len+1);
            break;

        default:
            ABORT ("invalid QUAL method %u", segconf.vcf_QUAL_method);
    }
}

static void vcf_seg_FILTER (VBlockVCFP vb, STRp(filter))
{
    decl_ctx (VCF_FILTER);

    if (segconf.vcf_is_manta || segconf.vcf_is_pbsv)
        vcf_seg_sv_copy_mate (vb, ctx, STRa(filter), TW_FILTER, TW_FILTER, false, filter_len+1);
    
    else
        seg_by_ctx (VB, STRa(filter), ctx, filter_len+1);
}

static inline void vcf_seg_ID (VBlockVCFP vb, STRp(id))
{
    decl_ctx(VCF_ID);
    
    if (flag.add_line_numbers) 
        vcf_seg_add_line_number (vb, id_len);

    else {
        if (segconf.vcf_is_manta)
            vcf_seg_manta_ID (vb, STRa(id));
        
        else if (segconf.vcf_is_svaba)
            vcf_seg_svaba_ID (vb, STRa(id));

        else if (segconf.vcf_is_pbsv)
            vcf_seg_pbsv_ID (vb, STRa(id));

        else if (IS_PERIOD(id))
            seg_by_ctx (VB, STRa(id), ctx, id_len+1); // often - all the same

        else
            seg_id_field (VB, ctx, STRa(id), false, id_len+1);
        
        seg_set_last_txt (VB, ctx, STRa(id));
    }
}

/* segment a VCF line into its fields:
   fields CHROM->FORMAT are "field" contexts
   all samples go into the SAMPLES context, which is a Container
   Haplotype and phase data are stored in a separate buffers + a SNIP_SPECIAL in the GT context  */
rom vcf_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    ZipDataLineVCF *dl = DATA_LINE (vb->line_i);

    rom next_field=field_start_line, field_start;
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
    vcf_seg_ID (vb, STRd(VCF_ID));

    CTX(FORMAT_RGQ)->line_has_RGQ = !segconf.running && segconf.has[FORMAT_RGQ] && vcf_refalt_seg_ref_alt_line_has_RGQ (VCF_ALT_str);

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

    vcf_seg_QUAL (vb, STRd (VCF_QUAL)); // seg after INFO as it might depends on DP

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

    // Adds DVCF items according to ostatus, finalizes INFO/SF and segs the INFO container
    vcf_seg_finalize_INFO_fields (vb);

    if (!has_samples && vcf_num_samples)
        WARN_ONCE ("FYI: variant CHROM=%.*s POS=%"PRId64" has no samples", vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS));

    SEG_EOL (VCF_EOL, false);

    return next_field;
}
