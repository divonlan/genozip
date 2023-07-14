// ------------------------------------------------------------------
//   vcf_seg.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "context.h"
#include "random_access.h"
#include "file.h"
#include "zip.h"
#include "dict_id.h"
#include "codec.h"
#include "strings.h"
#include "chrom.h"
#include "libdeflate/libdeflate.h"
#include "stats.h"
#include "segconf.h"
#include "gencomp.h"
#include "tip.h"
#include "arch.h"

static MiniContainer line_number_container = {};

// called by main thread after reading the header
void vcf_zip_initialize (void)
{
    vcf_samples_zip_initialize ();
    vcf_info_zip_initialize ();

    if (z_is_dvcf) 
        vcf_lo_zip_initialize ();

    // container just for adding a prefix to the delta-encoded line number (the container is all_the_same)
    if (flag.add_line_numbers && !line_number_container.repeats) { // possibly already initialized by previous files
        line_number_container = (MiniContainer) {
            .repeats   = 1,
            .nitems_lo = 1,
            .items     = { { .dict_id = { _VCF_LINE_NUM } } }
        };
    }

    // if sorting - we have to pre-populate header (in vcf_header_consume_contig) & reference contigs so that vb->is_unsorted is calculated correctly 
    // in vcf_seg_evidence_of_unsorted() (if VBs have their chrom nodes, they might be in reverse order vs the eventual z_file chroms)

    // NOTE: unlikely edge cases in which some variants sorted inconsistently:
    // 1. Two variants with the same chrom, start_pos, end_pos and tie_breaker(=Adler32 of REF+ALT)
    // 2. If two (or more) contigs appears in the variants but not in the VCF header nor in the reference file. If these first appear in
    //    reverse order, in two VBs running in parallel, then the "wrong" VB (compared to the eventual zctx) will be mis-sorted.
    if (vcf_is_sorting (VCF_COMP_MAIN)) {
        if (chain_is_loaded) {
            ctx_populate_zf_ctx_from_contigs (prim_ref, VCF_CHROM,  ref_get_ctgs (prim_ref)); 
            ctx_populate_zf_ctx_from_contigs (gref,     VCF_oCHROM, ref_get_ctgs (gref));
        }
        else
            ctx_populate_zf_ctx_from_contigs (gref, VCF_CHROM, ref_get_ctgs (gref)); 
    }
}

// called after each file
void vcf_zip_finalize (bool is_last_user_txt_file)
{
    if (flag.let_OS_cleanup_on_exit) return; // no need to waste time freeing if this is the last file - the process will die momentarily

    if (z_is_dvcf) gencomp_destroy();
}

// detect if a generic file is actually a vcf
bool is_vcf (STRp(header), bool *need_more)
{
    return header_len >= 16 && !memcmp (header, "##fileformat=VCF", 16);
}

// main thread: writing data-type specific fields to genozip header
void vcf_zip_genozip_header (SectionHeaderGenozipHeaderP header)
{
    header->vcf.segconf_has_RGQ = (segconf.has[FORMAT_RGQ] > 0); // introduced in v14
}

void vcf_zip_init_vb (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    // ZIP of a dual-coordinates file: calculate how much of the VB is rejected lines originating from ##primary_only/##luft_only
    vb->reject_bytes = MIN_(vb->recon_size, txt_file->reject_bytes);
    txt_file->reject_bytes -= vb->reject_bytes;

    vb->recon_size_luft = Ltxt; // initial value. it may change if --optimize / --chain are used, or if dual coordintes - for the other coordinate

    // set vcf_version in VB, since the the vcf_version in vcf_header might change as we might be reading the next txt file
    vb->vcf_version = vcf_header_get_version();

    vb->vb_coords = !z_is_dvcf                         ? DC_PRIMARY
                    : vb->comp_i == VCF_COMP_MAIN      ? DC_BOTH
                    : vb->comp_i == VCF_COMP_PRIM_ONLY ? DC_PRIMARY
                    : vb->comp_i == VCF_COMP_LUFT_ONLY ? DC_LUFT
                    :                                    -1; // invalid value

    vb->is_rejects_vb = z_is_dvcf && (vb->comp_i != VCF_COMP_MAIN);

    vb->sort = vcf_is_sorting (vb->comp_i);

    // in case we're replacing ID with the line number
    if (flag.add_line_numbers) {
        vb->first_line = txt_file->num_lines + 1; 
        txt_file->num_lines += str_count_char (STRb(vb->txt_data), '\n');  // update here instead of in zip_update_txt_counters;
    }
}

bool vcf_zip_vb_has_count (VBlockP vb)
{
    return !VB_VCF->is_rejects_vb;  // don't count DVCF rejects VB - these are duplicate lines counted in the normal VBs.
}

// called by main thread, as VBs complete (might be out-of-order)
void vcf_zip_after_compute (VBlockP vb)
{
    // note: VBs are out of order, impacting the neatness of the report. To solve, move to end of compute thread with serializer.
    if (chain_is_loaded && VB_VCF->rejects_report.len)
        buf_add_buf (evb, &z_file->rejects_report, &VB_VCF->rejects_report, char, "rejects_report");
}   

// ZIP main thread, called by zip_update_txt_counters after that vb has finished processing
void vcf_zip_update_txt_counters (VBlockP vb)
{
    // if we're compressing the primary-only rejects, they are not reconstructed in the default (primary) reconstruction
    if (vb->comp_i == VCF_COMP_PRIM_ONLY) 
        z_file->txt_data_so_far_bind -= vb->recon_size; // cancel increment by recon_size done by zip_update_txt_counters

    // if we're compressing the luft-only rejects, the default (primary) reconstruction will show the these lines (in their luft reconstruction, as they are luft-only)
    else if (vb->comp_i == VCF_COMP_LUFT_ONLY) 
        z_file->txt_data_so_far_bind += VB_VCF->recon_size_luft - vb->recon_size; // add recon_size_luft instead of (already added) recon_size
}

void vcf_zip_set_txt_header_flags (struct FlagsTxtHeader *f)
{
    f->is_txt_luft = (txt_file->coords == DC_LUFT);
}

void vcf_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeaderP vb_header)
{
    vb_header->dvcf_recon_size_luft = BGEN32 (VB_VCF->recon_size_luft);
}

uint32_t vcf_seg_get_vb_recon_size (VBlockP vb)
{
    ASSERT (VB_VCF->recon_size_luft >= 0, "recon_size_luft=%d is negative for vb_i=%u, coord=%s", VB_VCF->recon_size_luft, vb->vblock_i, vcf_coords_name(VB_VCF->vb_coords));

    return VB_VCF->vb_coords == DC_LUFT ? VB_VCF->recon_size_luft : vb->recon_size; // in primary reconstruction, ##luft_only VB is reconstructed in luft coords
}

// called by Compute threadfrom seg_all_data_lines
void vcf_seg_initialize (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    ctx_set_no_stons (VB, VCF_CHROM, VCF_oCHROM, VCF_FORMAT, VCF_INFO, VCF_oSTATUS, VCF_COORDS, 
                      VCF_TOPLEVEL, VCF_TOPLUFT, VCF_LIFT_REF, VCF_COPYPOS, VCF_oXSTRAND, 
                      VCF_POS, VCF_oPOS, VCF_LINE_NUM, INFO_HGVS_del_start_pos, INFO_HGVS_ins_start_pos, INFO_HGVS_ins_start_pos, // as required by seg_pos_field
                      DID_EOL);

    ctx_set_store (VB, STORE_INDEX, VCF_oSTATUS, VCF_COORDS, VCF_oXSTRAND, VCF_CHROM, VCF_oCHROM, DID_EOL);

    ctx_set_store (VB, STORE_INT, VCF_POS, VCF_oPOS, VCF_ID, VCF_LINE_NUM, FORMAT_DP, FORMAT_MIN_DP, INFO_DP, DID_EOL); 

    ctx_set_store (VB, STORE_FLOAT, VCF_QUAL, DID_EOL); // consumed by vcf_piz_special_QD

    ctx_set_ltype (VB, LT_DYN_INT, INFO_RAW_MQandDP_MQ, INFO_RAW_MQandDP_DP, INFO_dbSNPBuildID, DID_EOL);
    
    CTX(VCF_oCHROM)->  no_vb1_sort = true; // indices need to remain as in the Chain file
    CTX(VCF_oSTATUS)-> no_vb1_sort = true; // indices need to remaining matching to LiftOverStatus
    CTX(VCF_COORDS)->  no_vb1_sort = true; // indices need to remaining matching to Coords
    CTX(VCF_oXSTRAND)->no_vb1_sort = true; // indices need to order of ctx_create_node

    seg_id_field_init (CTX(VCF_ID));

    // counts sections
    CTX(VCF_CHROM)-> counts_section = true;
    CTX(VCF_oCHROM)->counts_section = true;

    if (z_is_dvcf) {
        CTX(VCF_oSTATUS)->counts_section = true;
        CTX(VCF_COORDS)-> counts_section = true;
    }

    // consolidate stats
    ctx_consolidate_stats (VB, VCF_REFALT, VCF_oREFALT, VCF_LIFT_REF, DID_EOL);
    ctx_consolidate_stats (VB, VCF_POS,    VCF_oPOS, VCF_COPYPOS, DID_EOL);
    ctx_consolidate_stats (VB, VCF_CHROM,  VCF_oCHROM, DID_EOL);
    ctx_consolidate_stats (VB, VCF_COORDS, INFO_PRIM, INFO_PREJ, INFO_LUFT, INFO_LREJ, VCF_oSTATUS, VCF_COPYSTAT, VCF_oXSTRAND, DID_EOL);

    // room for already existing FORMATs from previous VBs
    vb->format_mapper_buf.len = vb->format_contexts.len = CTX(VCF_FORMAT)->ol_nodes.len;
    buf_alloc_zero (vb, &vb->format_mapper_buf, 0, vb->format_mapper_buf.len, Container, CTX_GROWTH, "format_mapper_buf");
    buf_alloc_zero (vb, &vb->format_contexts, 0, vb->format_contexts.len, ContextPBlock, CTX_GROWTH, "format_contexts");
    
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

    // evaluate oSTATUS and COORDS, XSTRAND snips in order, as we rely on their indices being identical to the order of these arrays
    for (int i=0; i < NUM_LO_STATUSES; i++) 
        ctx_create_node (VB, VCF_oSTATUS, dvcf_status_names[i], strlen (dvcf_status_names[i]));

    for (int i=0; i < NUM_COORDS; i++) 
        ctx_create_node (VB, VCF_COORDS, vcf_coords_name(i), strlen (vcf_coords_name(i)));

    ctx_create_node (VB, VCF_oXSTRAND, cSTR("-")); // is_xstrand=false
    ctx_create_node (VB, VCF_oXSTRAND, cSTR("0")); // is_xstrand=true - REF and ALTs where rev-comped in place
    ctx_create_node (VB, VCF_oXSTRAND, cSTR("1")); // is_xstrand=true - REF and ALTs were rotated one base to the left due to re-left-anchoring

    // create a nodes and dict entry for LIFT_REF, COPYSTAT and COPYPOS - these become "all_the_same" so no need to seg them explicitly hereinafter
    ctx_create_node (VB, VCF_LIFT_REF, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_LIFT_REF }, 2);
    ctx_create_node (VB, VCF_COPYSTAT, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPYSTAT }, 2);
    ctx_create_node (VB, VCF_COPYPOS,  (char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPYPOS  }, 2);

    if (segconf.has[FORMAT_RGQ]) {
        seg_mux_init (VB, CTX(VCF_QUAL), 2, VCF_SPECIAL_MUX_BY_HAS_RGQ, false, (MultiplexerP)&vb->mux_QUAL, "01");
        seg_mux_init (VB, CTX(VCF_INFO), 2, VCF_SPECIAL_MUX_BY_HAS_RGQ, false, (MultiplexerP)&vb->mux_INFO, "01");        
    }

    vcf_info_seg_initialize(vb);
    vcf_samples_seg_initialize(vb);

    if (segconf.vcf_illum_gtyping)  vcf_illum_gtyping_initialize (vb);
    if (segconf.vcf_is_gwas)        vcf_gwas_seg_initialize (vb);
    if (segconf.vcf_is_cosmic)      vcf_cosmic_seg_initialize (vb);
    if (segconf.vcf_is_vep)         vcf_vep_seg_initialize (vb);
    if (segconf.vcf_is_mastermind)  vcf_mastermind_seg_initialize (vb);
    if (segconf.vcf_is_dbSNP)       vcf_dbsnp_seg_initialize (vb);
    if (segconf.vcf_is_giab_trio)   vcf_giab_seg_initialize (vb);
    if (segconf.vcf_is_gnomad)      CTX(VCF_QUAL)->no_stons = true;
}             

static void vcf_seg_finalize_segconf (VBlockVCFP vb)
{
    // identify DRAGEN GVCF. GATK's is identified in vcf_inspect_txt_header_zip()
    if (segconf.has[FORMAT_ICNT] && segconf.has[FORMAT_SPL]) 
        segconf.vcf_is_gvcf = true;

    if (segconf.has[FORMAT_IGT] && segconf.has[FORMAT_IPS] && segconf.has[FORMAT_ADALL])
        segconf.vcf_is_giab_trio = true;
        
    // an alternative way to discover dbSNP if ##source is missing
    if (!!segconf.has[INFO_RS] + !!segconf.has[INFO_RSPOS] + !!segconf.has[INFO_SAO] + !!segconf.has[INFO_SSR] +
        !!segconf.has[INFO_dbSNPBuildID] + !!segconf.has[INFO_VC] + !!segconf.has[INFO_NSM] + !!segconf.has[INFO_U3] >= 3)
        segconf.vcf_is_dbSNP = true;

    // in gnomAD, we have a huge number of INFO fields in various permutations - generating a huge INFO dictionary, but which compresses very very well
    if (segconf.vcf_is_gnomad)
        ZCTX(VCF_INFO)->dict_len_excessive = true; // don't warn if excessive
        
    if (!flag.reference && segconf.vcf_is_gvcf)
        TIP ("Compressing a GVCF file using a reference file can reduce the compressed file's size by 10%%-30%%.\n"
             "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n",
             arch_get_argv0(), txt_file->name);
}

void vcf_seg_finalize (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    
    if (vb->ht_matrix_ctx) 
        vcf_seg_FORMAT_GT_complete_missing_lines (vb);

    // for a dual-coordinates VCF, we offer 2 ways to reconstruct it: normally, it is reconstructed in the
    // primary coordinates. --luft invokes top_level_luft in translated mode, which reconstructs in luft coordintes.
    
    // Toplevel snip for reconstructing this VB a PRIMARY
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .callback     = (vb->use_special_sf == USE_SF_YES) || z_is_dvcf, // cases where we need a callback
        .filter_items = true,
        .nitems_lo    = 12,                                                                 
        .items        = { { .dict_id = { _VCF_COORDS },  .separator = "\t" }, // suppressed by vcf_piz_filter unless --show-dvcf                                   
                          { .dict_id = { _VCF_oSTATUS }, .separator = "\t" }, // suppressed by vcf_piz_filter unless --show-dvcf                                   
                          { .dict_id = { _VCF_CHROM },   .separator = "\t" },
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

    Context *ctx = CTX(VCF_TOPLEVEL);

    if (vb->vb_coords == DC_BOTH || !z_is_dvcf)
        container_seg (vb_, ctx, (ContainerP)&top_level, 0, 0, 0); 

    // when processing the rejects file containing variants that are primary-only, we add a "##primary_only=" prefix to 
    // first item of each line, so that it reconstructs as part of the VCF header 
    else if (vb->vb_coords == DC_PRIMARY) { // primary-only variants 
        static const char primary_only_prefix[] = CON_PX_SEP_ CON_PX_SEP_ HK_PRIM_ONLY CON_PX_SEP_;
        container_seg (vb_, ctx, (ContainerP)&top_level, primary_only_prefix, strlen (primary_only_prefix), 0);
        vb->recon_size += (sizeof HK_PRIM_ONLY - 1) * vb->lines.len; // when reconstructing primary-only rejects, we also reconstruct a prefix for each line
        ctx->txt_len += (sizeof HK_PRIM_ONLY - 1) * vb->lines.len;
    }
    
    // Toplevel snip for reconstructing this VB a LUFT
    SmallContainer top_luft = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .callback     = (vb->use_special_sf == USE_SF_YES) || z_is_dvcf, // cases where we need a callback
        .filter_items = true,
        .nitems_lo    = 13,                                                                 
        .items        = { { .dict_id = { _VCF_COORDS },  .separator = "\t" }, // suppressed by vcf_piz_filter unless --show-dvcf                                   
                          { .dict_id = { _VCF_oSTATUS }, .separator = "\t" }, // suppressed by vcf_piz_filter unless --show-dvcf                                   
                          { .dict_id = { _VCF_oCHROM },  .separator = "\t" },
                          { .dict_id = { _VCF_oPOS },    .separator = "\t" },
                          { .dict_id = { _VCF_ID },      .separator = "\t" },
                          { .dict_id = { _VCF_oREFALT }, .separator = "\t" },
                          { .dict_id = { _VCF_QUAL },    .separator = "\t" },
                          { .dict_id = { _VCF_FILTER },  .separator = "\t" },
                          { .dict_id = { _VCF_POS },     .separator = { CI0_TRANS_NOR } }, // consume POS before INFO, in case we have INFO/END
                          { .dict_id = { _VCF_INFO },    .separator = "\t" }, // in dual-coordinates, contains INFO/LIFTOVER or INFO/REJTOVER that reconstructs oCHROM, oPOS, oREF, oXSTRAND
                          { .dict_id = { _VCF_FORMAT },  .separator = "\t" },
                          { .dict_id = { _VCF_SAMPLES }, .separator = ""   },
                          { .dict_id = { _VCF_EOL },     .separator = ""   } }
    };

    if (vb->vb_coords == DC_BOTH)
        container_seg (vb_, CTX(VCF_TOPLUFT), (ContainerP)&top_luft, 0, 0, 0);

    // similarly, when processing the rejects file containing variants that are luft-only, we add a "##luft_only=" prefix
    else if (vb->vb_coords == DC_LUFT) { // luft-only variants 
        static const char luft_only_prefix[] = CON_PX_SEP_ CON_PX_SEP_ HK_LUFT_ONLY CON_PX_SEP_;
        container_seg (vb_, CTX(VCF_TOPLUFT), (ContainerP)&top_luft, luft_only_prefix, strlen (luft_only_prefix), (sizeof HK_LUFT_ONLY - 1) * vb->lines.len);
        vb->recon_size_luft += (sizeof HK_LUFT_ONLY - 1) * vb->lines.len; // when reconstructing luft-only rejects, we also reconstruct a prefix for each line
        // note: there is no equivalent of ctx->txt_len for Luft coordinates
    }

    vb->flags.vcf.coords = vb->vb_coords;

    vcf_samples_seg_finalize (vb);

    if (segconf.running)
        vcf_seg_finalize_segconf (vb);
}

// after each VB is compressed and merge (VB order is arbitrary)
void vcf_zip_after_compress (VBlockP vb)
{
    // case: we're sorting the file's lines - add line data to txt_file's line_info (VB order doesn't matter - we will sort them later). 
    if (VB_VCF->sort || z_is_dvcf) 
        vcf_linesort_merge_vb (vb);

    if (VB_VCF->PL_mux_by_DP == PL_mux_by_DP_TEST) 
        vcf_FORMAT_PL_decide (VB_VCF);

    // Only the MAIN component produces gencomp lines, however we are processing VBs in order, so out-of-band VBs
    // need to be sent too, just to advance the serializing mutex
    if (z_is_dvcf) 
        gencomp_absorb_vb_gencomp_lines (vb);
}

// called after all VBs are compressed - before Global sections are compressed
void vcf_zip_after_vbs (void)
{
    vcf_FORMAT_PL_after_vbs ();
}

bool vcf_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return 
        dict_id.num == _VCF_TOPLEVEL ||
        dict_id.num == _VCF_TOPLUFT  ||
        dict_id.num == _VCF_CHROM    ||
        dict_id.num == _VCF_oCHROM   ||
        dict_id.num == _VCF_FORMAT   || // note: NOT including _VCF_INFO - there are cases where it is not small
        dict_id.num == _VCF_FILTER   || 
        dict_id.num == _VCF_EOL      ||
        dict_id.num == _VCF_SAMPLES  ||
        dict_id.num == _VCF_oCHROM   ||
        dict_id.num == _VCF_oXSTRAND ||
        dict_id.num == _VCF_oSTATUS  ||
        dict_id.num == _VCF_COORDS   ||
        dict_id.num == _VCF_LIFT_REF ||
        dict_id.num == _INFO_DP      ||
        (dict_id.num == _INFO_AA && !segconf.vcf_is_cosmic) || // stored as a SPECIAL snip
        dict_id.num == _INFO_LDAF    ||
        dict_id.num == _INFO_MQ0     ||
        dict_id.num == _INFO_LUFT    ||
        dict_id.num == _INFO_PRIM    ||
        dict_id.num == _INFO_LREJ    ||
        dict_id.num == _INFO_PREJ;
        // note: _INFO_AN, _INFO_AC, _INFO_AF, _INFO_MLEAC, _INFO_MLEAF - can be big (e.g. big in GWAS VCF, ExAC, gnomad...)
}

bool vcf_seg_is_big (ConstVBlockP vb, DictId dict_id, DictId st_dict_id/*dict_id of st_did_i*/)
{
    return 
        dict_id.num    == _VCF_REFALT  ||
        dict_id.num    == _VCF_oREFALT ||
        dict_id.num    == _INFO_RAW_MQ ||
        st_dict_id.num == _FORMAT_PLn  || // P1Ln1 etc are often big
        st_dict_id.num == _FORMAT_PLy  ||
        st_dict_id.num == _FORMAT_GP    ; // G1P1 etc are often big
}

// returns length of gencomp before the copying
static uint32_t vcf_seg_get_line_len (VBlockVCFP vb, rom field_start_line, uint32_t remaining_txt_len)
{
    rom last = memchr (field_start_line, '\n', remaining_txt_len);
    ASSERT (last, "Line has no newline: %.*s", remaining_txt_len, field_start_line);

    return last - field_start_line + 1;
}

static inline LineCmpInfo vcf_seg_make_lci (ZipDataLineVCF *dl, bool is_luft)
{
    return (LineCmpInfo){ .chrom_wi    = dl->chrom[is_luft],
                          .start_pos   = dl->pos[is_luft],
                          .end_pos     = dl->pos[is_luft] + MAX_(0, dl->end_delta),
                          .tie_breaker = dl->tie_breaker };
}  

static inline void vcf_seg_evidence_of_unsorted (VBlockVCFP vb, ZipDataLineVCF *dl, Did chrom_did_i)
{
    bool is_luft = !!chrom_did_i; // 0 for primary, 1 for last

    if (!vb->is_unsorted[is_luft] && vb->line_i > 0)     
        vb->is_unsorted[is_luft] = vcf_linesort_cmp (vcf_seg_make_lci (dl-1, is_luft), vcf_seg_make_lci (dl, is_luft)) > 0;
}

static void vcf_seg_add_line_number (VBlockVCFP vb, unsigned VCF_ID_len)
{
    char line_num[20];
    unsigned line_num_len = str_int (vb->first_line + vb->line_i, line_num);

    seg_pos_field (VB, VCF_LINE_NUM, VCF_LINE_NUM, SPF_ZERO_IS_BAD, 0, line_num, line_num_len, 0, line_num_len + 1 + LN_PREFIX_LEN); // +1 for tab +3 for "LN="
    
    int shrinkage = (int)VCF_ID_len - line_num_len - LN_PREFIX_LEN;
    vb->recon_size -= shrinkage;
    vb->recon_size_luft -= shrinkage;
}

// calculate tie-breaker for sorting 
static void vcf_seg_assign_tie_breaker (VBlockVCFP vb, ZipDataLineVCF *dl)
{
    // if this is a dual-coord line, and it is a LUFT line, convert REF\tALT to PRIMARY
    LiftOverStatus ostatus = last_ostatus;
    if (vb->line_coords == DC_LUFT && LO_IS_OK (ostatus))
        vcf_refalt_seg_convert_to_primary (vb, ostatus);

    // tie breaker is crc32 of the Primary REF\tALT, except for Luft-only lines, where it is the Adler32 of the Luft REF\tALT
    dl->tie_breaker = crc32 (1, vb->main_ref, vb->main_ref_len + 1 + vb->main_alt_len);
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

static inline void vcf_seg_QUAL (VBlockVCFP vb, STRp(qual))
{
    set_last_txt (VCF_QUAL, qual);

    // gnomAD - store in local
    if (segconf.vcf_is_gnomad || segconf.vcf_is_dbSNP)
        seg_add_to_local_text (VB, CTX(VCF_QUAL), STRa(qual), LOOKUP_NONE, qual_len+1);

    // case: GVCF - multiplex by has_RGQ
    else if (!segconf.running && segconf.has[FORMAT_RGQ]) {
        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, VCF_QUAL, (MultiplexerP)&vb->mux_QUAL, vb->line_has_RGQ);
        seg_by_ctx (VB, STRa(qual), channel_ctx, qual_len+1);
        seg_by_did (VB, STRa(vb->mux_QUAL.snip), VCF_QUAL, 0);
    }
    
    // case: not GVCF
    else
        seg_by_did (VB, STRa(qual), VCF_QUAL, qual_len+1);
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

    vb->line_coords = (vb->comp_i == VCF_COMP_PRIM_ONLY) ? DC_PRIMARY // this is the Primary-only rejects component
                    : (vb->comp_i == VCF_COMP_LUFT_ONLY) ? DC_LUFT    // this is the Luft-only rejects component
                    : (txt_file->coords == DC_LUFT)      ? DC_LUFT    // The main component, when compressing a DVCF with ##dual_coordinates=LUFT
                    :                                      DC_PRIMARY;

    ASSERT (!vb->reject_bytes || (vb->reject_bytes > 0 && vb->comp_i == VCF_COMP_MAIN), 
            "Expecting reject_bytes=%d >= 0 and they can only appear in comp_i=%s == MAIN", vb->reject_bytes, comp_name (vb->comp_i));

    if (vb->reject_bytes) 
        vb->line_coords = OTHER_COORDS (vb->line_coords); // dual coordinates file - this line originates from ##primary_only/##luft_only as is in the opposite coordinates

    // allow un-segging of o* and PRIM/LUFT in case due to INFO or FORMAT ostatus changes from OK to REJECT
    int32_t save_recon_size=0; uint32_t line_len=0; 
    unsigned save_txt_len_len = vb->num_contexts;
    uint64_t save_txt_len[save_txt_len_len];
    if (z_is_dvcf) {
        line_len = vcf_seg_get_line_len (vb, field_start_line, remaining_txt_len);
        vcf_lo_set_rollback_point (vb); // rollback point for rolling back seg of the OTHER coord if it turns out this line cannot be luft-translated

        // LUFT-only lines are not reconstructed by default. Save the txt_len of all contexts to be able to undo their increment if it turns
        // out this LUFT-coordinate line is in fact LUFT-only.
        if (vb->line_coords == DC_LUFT) {
            save_recon_size = vb->recon_size;
            for (unsigned i=0; i < save_txt_len_len; i++)
                save_txt_len[i] = CTX(i)->txt_len; // txt_len is for PRIMARY reconstruction ; we don't calculate txt_len for LUFT
        }
        
        // make sure we don't change recon_size_luft if segging a primary-only line
        else if (vb->line_coords == DC_PRIMARY) 
            save_recon_size = vb->recon_size_luft;
    }

    GET_NEXT_ITEM (VCF_CHROM);
    if (vb->line_coords == DC_PRIMARY) 
        dl->chrom[0] = chrom_seg_by_did_i (VB, VCF_CHROM, STRd(VCF_CHROM),VCF_CHROM_len+1);
    
    else { // LUFT
        dl->chrom[1] = chrom_seg_by_did_i (VB, VCF_oCHROM, STRd(VCF_CHROM), VCF_CHROM_len);
        CTX(vb->vb_coords==DC_LUFT ? VCF_oCHROM : VCF_CHROM)->txt_len++; // account for the tab - in oCHROM in the ##luft_only VB and in CHROM (on behalf on the primary CHROM) if this is a Dual-coord line (we will rollback accounting later if its not)
    }

    GET_NEXT_ITEM (VCF_POS);
    set_last_txt_(VCF_POS, VCF_POS_str, VCF_POS_len); // consumed by vcf_seg_FORMAT_PS, vcf_seg_ILLUMINA_POS

    if (vb->line_coords == DC_PRIMARY) {
        PosType64 pos = dl->pos[0] = seg_pos_field (vb_, VCF_POS, VCF_POS, 0, '.', STRd(VCF_POS), 0, VCF_POS_len+1);

        if (pos == 0 && !(*VCF_POS_str == '.' && VCF_POS_len == 1)) // POS == 0 - invalid value return from seg_pos_field
            WARN_ONCE ("FYI: invalid POS=%"PRId64" value in chrom=%.*s vb_i=%u vb_line_i=%d: line will be compressed, but not indexed", 
                       pos, vb->chrom_name_len, vb->chrom_name, vb->vblock_i, vb->line_i);
                
        if (pos) random_access_update_pos (vb_, 0, VCF_POS);
    }

    else { // LUFT
        dl->pos[1] = seg_pos_field (vb_, VCF_oPOS, VCF_oPOS, 0, '.', STRd(VCF_POS), 0, VCF_POS_len);
        if (dl->pos[1]) random_access_update_pos (vb_, 1, VCF_oPOS);
        CTX(vb->vb_coords==DC_LUFT ? VCF_oPOS : VCF_POS)->txt_len++; // account for the tab - in oPOS in the ##luft_only VB and in POS (on behalf on the primary POS) if this is a Dual-coord line (we will rollback accounting later if its not)
    }

    GET_NEXT_ITEM (VCF_ID);

    if (flag.add_line_numbers) 
        vcf_seg_add_line_number (vb, VCF_ID_len);
    else {
        seg_id_field (VB, CTX(VCF_ID), VCF_ID_str, VCF_ID_len, true);
        seg_set_last_txt (VB, CTX(VCF_ID), STRd(VCF_ID));
    }

    // REF + ALT 
    GET_NEXT_ITEM (VCF_REF);
    GET_NEXT_ITEM (VCF_ALT);

    // save REF and ALT (in primary or luft coordinates) to be used for INFO fields
    vb->main_ref     = VCF_REF_str;
    vb->main_alt     = VCF_ALT_str;
    vb->main_ref_len = VCF_REF_len;
    vb->main_alt_len = VCF_ALT_len;
 
    vb->line_has_RGQ = !segconf.running && segconf.has[FORMAT_RGQ] && vcf_refalt_seg_ref_alt_line_has_RGQ (VCF_ALT_str);

    // note: we treat REF+\t+ALT as a single field because REF and ALT are highly corrected, in the case of SNPs:
    // e.g. GG has a probability of 0 and GC has a higher probability than GA.
    vcf_refalt_seg_main_ref_alt (vb, STRd(VCF_REF), STRd(VCF_ALT));
    
    GET_NEXT_ITEM (VCF_QUAL);

    vcf_seg_QUAL (vb, STRd (VCF_QUAL));

    SEG_NEXT_ITEM (VCF_FILTER);
    
    // if --chain, seg dual coordinate record - lift over CHROM, POS and REFALT to luft coordinates
    if (chain_is_loaded)
        vcf_lo_seg_generate_INFO_DVCF (vb, dl);

    // INFO
    if (vcf_num_samples) {
        GET_NEXT_ITEM (VCF_INFO); // to do: currently, we require FORMAT (and hence a \t after INFO) in the case the file has samples but this line doesn't. VCF spec doesn't require FORMAT in this case.
    } else {
        GET_MAYBE_LAST_ITEM (VCF_INFO); // may or may not have a FORMAT field
    }

    // set oSTATUS except --chain (already set)
    if (!chain_is_loaded) 
        vcf_set_ostatus (vb->reject_bytes ? LO_REJECTED : // reject lines in Luft are all rejected (happens only in txt_file->coords==DC_LUFT)
                         txt_file->coords ? LO_UNKNOWN  : // we don't know yet, we will test for existance of INFO/*rej in vcf_seg_info_field_correct_for_dual_coordinates
                                            LO_NA)      ; // this z_file is not a dual-coordinates file

    vcf_seg_info_subfields (vb, field_start, field_len);

    bool has_samples = false;
    if (separator != '\n') { // has a FORMAT field

        // FORMAT
        if (vcf_num_samples) {
            GET_MAYBE_LAST_ITEM (VCF_FORMAT); // possibly no samples this line
        } else {
            GET_LAST_ITEM (VCF_FORMAT);
        }

        vcf_seg_FORMAT (vb, dl, field_start, field_len);

        if ((has_samples = (separator != '\n'))) { 

            ASSVCF0 (dl->has_genotype_data || dl->has_haplotype_data, "expecting line to end as it has no sample data, but it has not");
            
            // seg all samples. note that this is destructive: Samples are first lift-back to PRIMARY if needed 
            next_field = vcf_seg_samples (vb, dl, &len, (char*)next_field, has_13); 
        }
        else 
            seg_by_did (VB, NULL, 0, VCF_SAMPLES, 0); // case no samples: WORD_INDEX_MISSING
    }

    // case no format or samples
    else {
        seg_by_did (VB, NULL, 0, VCF_FORMAT, 0); 
        seg_by_did (VB, NULL, 0, VCF_SAMPLES, 0);
    }

    // Adds DVCF items according to ostatus, finalizes INFO/SF and segs the INFO container
    vcf_finalize_seg_info (vb);

    // calculate tie-breaker for sorting - do after INFO before segging the samples as GT data can overwrite REF ALT
    if (vb->sort) vcf_seg_assign_tie_breaker (vb, dl);

    if (!has_samples && vcf_num_samples)
        WARN_ONCE ("FYI: variant CHROM=%.*s POS=%"PRId64" has no samples", vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS));

    if (z_has_gencomp) { // DVCF
        // note: we don't seg Coord for non-DC files, despite being in TOPLEVEL. That's ok, it will be treated as
        // as "all_the_same" and have word_index=0 == NONE for all lines
        bool lo_ok = LO_IS_OK (last_ostatus);

        Coords reconstructable_coords = lo_ok ? DC_BOTH : vb->line_coords;
        rom name = vcf_coords_name (reconstructable_coords); 
        seg_by_did (VB, name, strlen (name), VCF_COORDS, 0); // 0 as its not in the txt data
                
        // case: line was rejected
        if (!lo_ok) {
            gencomp_seg_add_line (VB, vb->line_coords, field_start_line, line_len);

            // in a VCF_COMP_MAIN VB, single-coordinate lines won't be displayed in the opposite reconstruction
            if (vb->line_coords == DC_PRIMARY) vb->recon_size_luft = save_recon_size - line_len;
            else  /* DC_LUFT  */               vb->recon_size      = save_recon_size - line_len;
        
            // unaccount for this line, if its a Luft-only line in a dual-coordinate variant (i.e. main VCF data line) as it won't appear in the default reconstruction
            if (vb->vb_coords == DC_BOTH && vb->line_coords == DC_LUFT)
                for (unsigned i=0; i < vb->num_contexts; i++) 
                    CTX(i)->txt_len = (i < save_txt_len_len ? save_txt_len[i] : 0);
        }
        
        // in case of a reject line - update reject_bytes to indicate its consumption (note: rejects lines, if any, are at the beginning of the VB)
        if (vb->reject_bytes)
            vb->reject_bytes -= next_field - field_start_line;
    } 

    SEG_EOL (VCF_EOL, false);

    // test if still sorted
    if (vb->sort) {
        vcf_seg_evidence_of_unsorted (vb, dl, VCF_CHROM);
        if (z_is_dvcf) vcf_seg_evidence_of_unsorted (vb, dl, VCF_oCHROM);
    }

    return next_field;
}
