// ------------------------------------------------------------------
//   vcf_seg.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "vcf_private.h"
#include "seg.h"
#include "context.h"
#include "random_access.h"
#include "file.h"
#include "zip.h"
#include "dict_id.h"
#include "codec.h"
#include "strings.h"
#include "coords.h"
#include "chrom.h"
#include "linesorter.h"
#include "libdeflate/libdeflate.h"
#include "stats.h"
#include "segconf.h"

static MiniContainer line_number_container = {};

// called by main thread after reading the header
void vcf_zip_initialize (void)
{
    vcf_lo_zip_initialize ();
    vcf_samples_zip_initialize ();
    vcf_info_zip_initialize ();

    // container just for adding a prefix to the delta-encoded line number (the container is all_the_same)
    if (flag.add_line_numbers && !line_number_container.repeats) { // possibly already initialized by previous files
        line_number_container = (MiniContainer) {
            .repeats   = 1,
            .nitems_lo = 1,
            .items     = { { .dict_id = { _VCF_LINE_NUM } } }
        };
    }

    // if flag.sort - we have to pre-populate header (in vcf_header_consume_contig) & reference contigs so that vb->is_unsorted is calculated correctly 
    // in vcf_seg_evidence_of_unsorted() (if VBs have their chrom nodes, they might be in reverse order vs the eventual z_file chroms)

    // NOTE: unlikely edge cases in which some variants sorted consistently:
    // 1. Two variants with the same chrom, start_pos, end_pos and tie_breaker(=Adler32 of REF+ALT)
    // 2. If two (or more) contigs appears in the variants but not in the VCF header nor in the reference file. If these first appear in
    //    reverse order, in two VBs running in parallel, then the "wrong" VB (compared to the eventual zctx) will be mis-sorted.
    if (flag.sort && z_file->num_txt_components_so_far == 1) {
        if (chain_is_loaded) {
            ctx_populate_zf_ctx_from_contigs (prim_ref, VCF_CHROM,  ref_get_ctgs (prim_ref)); 
            ctx_populate_zf_ctx_from_contigs (gref,     VCF_oCHROM, ref_get_ctgs (gref));
        }
        else
            ctx_populate_zf_ctx_from_contigs (gref, VCF_CHROM, ref_get_ctgs (gref)); 
    }
}

void vcf_zip_read_one_vb (VBlockP vb)
{
    // set vcf_version in VB, since the the vcf_version in vcf_header might change as we might be reading the next txt file
    VB_VCF->vcf_version = vcf_header_get_version();

    // in case we're replacing ID with the line number
    if (flag.add_line_numbers) {
        vb->first_line = txt_file->num_lines + 1; // this is normally not used in ZIP
        txt_file->num_lines += str_count_char (STRb(vb->txt_data), '\n');  // update here instead of in zip_update_txt_counters;
    }
}

void vcf_zip_after_compute (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    if (z_dual_coords && !flag.rejects_coord) { // normal file, not zipping rejects files
        if (vb->lo_rejects[DC_PRIMARY-1].len) vcf_lo_append_rejects_file (vb_, DC_PRIMARY);
        if (vb->lo_rejects[DC_LUFT-1].len)    vcf_lo_append_rejects_file (vb_, DC_LUFT);
    }

    if (chain_is_loaded && vb->rejects_report.len) 
        buf_add_buf (evb, &z_file->rejects_report, &vb->rejects_report, char, "rejects_report");
}   

// called by Compute threadfrom seg_all_data_lines
void vcf_seg_initialize (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    CTX(VCF_CHROM)->   no_stons    = true; // needs b250 node_index for random access and reconstrution plan
    CTX(VCF_oCHROM)->  no_stons    = true; // same
    CTX(VCF_FORMAT)->  no_stons    = true;
    CTX(VCF_INFO)->    no_stons    = true;
    CTX(VCF_oSTATUS)-> no_stons    = true;
    CTX(VCF_COORDS)->  no_stons    = true;
    CTX(VCF_TOPLEVEL)->no_stons    = true; 
    CTX(VCF_TOPLUFT)-> no_stons    = true; 
    CTX(VCF_LIFT_REF)->no_stons    = true; 
    CTX(VCF_COPYPOS)-> no_stons    = true; 
    CTX(VCF_oXSTRAND)->no_stons    = true; // keep in b250 so it can be eliminated as all_the_same
    CTX(VCF_oCHROM)->  no_vb1_sort = true; // indices need to remain as in the Chain file
    CTX(VCF_oSTATUS)-> no_vb1_sort = true; // indices need to remaining matching to LiftOverStatus
    CTX(VCF_COORDS)->  no_vb1_sort = true; // indices need to remaining matching to Coords
    CTX(VCF_oXSTRAND)->no_vb1_sort = true; // indices need to order of ctx_create_node

    CTX(VCF_oSTATUS)-> flags.store = STORE_INDEX;
    CTX(VCF_COORDS)->  flags.store = STORE_INDEX;
    CTX(VCF_oXSTRAND)->flags.store = STORE_INDEX;
    CTX(VCF_CHROM)->   flags.store = STORE_INDEX; // since v12
    CTX(VCF_oCHROM)->  flags.store = STORE_INDEX; // used by regions_is_site_included
    CTX(VCF_POS)->     flags.store = STORE_INT;   // since v12
    CTX(VCF_oPOS)->    flags.store = STORE_INT;   // used by vcf_piz_luft_END

    seg_id_field_init (CTX(VCF_ID));

    // counts sections
    CTX(VCF_CHROM)->   counts_section = true;
    CTX(VCF_oCHROM)->  counts_section = true;

    if (z_dual_coords) {
        CTX(VCF_oSTATUS)->counts_section = true;
        CTX(VCF_COORDS)-> counts_section = true;
    }

    // consolidate stats
    stats_set_consolidation (VB, VCF_REFALT, 2, VCF_oREFALT, VCF_LIFT_REF);
    stats_set_consolidation (VB, VCF_POS,    2, VCF_oPOS, VCF_COPYPOS);
    stats_set_consolidation (VB, VCF_CHROM,  1, VCF_oCHROM);
    stats_set_consolidation (VB, VCF_COORDS, 7, INFO_PRIM, INFO_PREJ, INFO_LUFT, INFO_LREJ, VCF_oSTATUS, VCF_COPYSTAT, VCF_oXSTRAND);

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
        ctx_create_node (VB, VCF_COORDS, coords_name(i), strlen (coords_name(i)));

    ctx_create_node (VB, VCF_oXSTRAND, cSTR("-")); // is_xstrand=false
    ctx_create_node (VB, VCF_oXSTRAND, cSTR("0")); // is_xstrand=true - REF and ALTs where rev-comped in place
    ctx_create_node (VB, VCF_oXSTRAND, cSTR("1")); // is_xstrand=true - REF and ALTs were rotated one base to the left due to re-left-anchoring

    // create a nodes and dict entry for LIFT_REF, COPYSTAT and COPYPOS - these become "all_the_same" so no need to seg them explicitly hereinafter
    seg_by_did_i (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_LIFT_REF }), 2, VCF_LIFT_REF, 0); 
    ctx_decrement_count (vb_, CTX(VCF_LIFT_REF), 0);

    seg_by_did_i (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPYSTAT }), 2, VCF_COPYSTAT, 0); 
    ctx_decrement_count (vb_, CTX(VCF_COPYSTAT), 0);

    seg_by_did_i (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPYPOS }), 2, VCF_COPYPOS, 0);
    ctx_decrement_count (vb_, CTX(VCF_COPYPOS), 0);

    vcf_info_seg_initialize(vb);
    vcf_samples_seg_initialize(vb);

    if (segconf.running) {
        if (vcf_num_samples > 1) segconf.INFO_DP_by_FORMAT_DP = true; // unless proven otherwise
    }
}             

void vcf_seg_finalize (VBlockP vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    
    if (vb->ht_matrix_ctx) 
        vcf_seg_FORMAT_GT_complete_missing_lines (vb);

    // for a dual-coordinates VCF, we offer 2 ways to reconstruct it: normally, it is reconstructed in the
    // primary coordinates. --luft invokes top_level_luft in translated mode, which reconstructs in luft coordintes.
    
    // Toplevel snip for reconstructing this VB a PRIMARY
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .callback     = (vb->use_special_sf == USE_SF_YES) || z_dual_coords, // cases where we need a callback
        .filter_items = true,
        .nitems_lo    = 12,                                                                 
        .items        = { { .dict_id = { _VCF_COORDS },  .separator = "\t" }, // suppressed by vcf_piz_filter unless --show-dvcf                                   
                          { .dict_id = { _VCF_oSTATUS }, .separator = "\t" }, // suppressed by vcf_piz_filter unless --show-dvcf                                   
                          { .dict_id = { _VCF_CHROM },   .separator = "\t" },
                          { .dict_id = { _VCF_POS },     .separator = "\t" },
                          { .dict_id = { _VCF_ID },      .separator = "\t" },
                          { .dict_id = { _VCF_REFALT },  .separator = "\t" },
                          { .dict_id = { _VCF_QUAL },    .separator = "\t" },
                          { .dict_id = { _VCF_FILTER },  .separator = "\t" },
                          { .dict_id = { _VCF_INFO },    .separator = "\t" }, // in dual-coordinates, contains INFO/LIFTOVER or INFO/REJTOVER that reconstructs oCHROM, oPOS, oREF, oXSTRAND
                          { .dict_id = { _VCF_FORMAT },  .separator = "\t" },
                          { .dict_id = { _VCF_SAMPLES }, .separator = ""   },
                          { .dict_id = { _VCF_EOL },     .separator = ""   } },
    };

    Context *ctx = CTX(VCF_TOPLEVEL);

    if (vb->vb_coords == DC_BOTH || !z_dual_coords)
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
        .callback     = (vb->use_special_sf == USE_SF_YES) || z_dual_coords, // cases where we need a callback
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
}

// after each VB is compressed
void vcf_zip_after_compress (VBlockP vb)
{
    if (VB_VCF->PL_mux_by_DP == PL_mux_by_DP_TEST) 
        vcf_FORMAT_PL_decide (VB_VCF);
}

// called after all VBs are compressed - before Global sections are compressed
void zip_vcf_after_vbs (void)
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
        dict_id.num == _VCF_FORMAT   ||
        dict_id.num == _VCF_INFO     ||
        dict_id.num == _VCF_REFALT   ||
        dict_id.num == _VCF_oREFALT  ||
        dict_id.num == _VCF_FILTER   ||
        dict_id.num == _VCF_EOL      ||
        dict_id.num == _VCF_SAMPLES  ||
        dict_id.num == _VCF_oCHROM   ||
        dict_id.num == _VCF_oXSTRAND ||
        dict_id.num == _VCF_oSTATUS  ||
        dict_id.num == _VCF_COORDS   ||
        dict_id.num == _VCF_LIFT_REF ||
        dict_id.num == _INFO_AC              ||
        dict_id.num == _INFO_AF              ||
        dict_id.num == _INFO_AN              ||
        dict_id.num == _INFO_DP              ||
        dict_id.num == _INFO_AA              || // stored as a SPECIAL snip
        dict_id.num == _INFO_MLEAC           ||
        dict_id.num == _INFO_MLEAF           ||
        dict_id.num == _INFO_LDAF            ||
        dict_id.num == _INFO_MQ0             ||
        dict_id.num == _INFO_LUFT            ||
        dict_id.num == _INFO_PRIM            ||
        dict_id.num == _INFO_LREJ            ||
        dict_id.num == _INFO_PREJ            ||

        // INFO/ AC_* AN_* AF_* and ???_AF are small
        ((dict_id.id[0] == ('A' | 0xc0)) && (dict_id.id[1] == 'C' || dict_id.id[1] == 'F' || dict_id.id[1] == 'N') && dict_id.id[2] == '_') ||
        (dict_id_is_vcf_info_sf (dict_id) && dict_id.id[3] == '_' && dict_id.id[4] == 'A' && dict_id.id[5] == 'F' && !dict_id.id[6]);
}


// converts a FORMAT field name to a dict_id
static DictId vcf_seg_get_format_sf_dict_id (STRp (sf_name)) 
{
    DictId dict_id;
    // case: normal field - starts with a letter or another character in the range
    if (sf_name[0] >= 64 && sf_name[0] <= 127) 
        dict_id = dict_id_make (sf_name, sf_name_len, DTYPE_VCF_FORMAT);
    
    // case: unusual field - starts with an out-of-range character, eg a digit - prefix with @ so its a legal FORMAT dict_id
    else {
        SAFE_ASSIGN (sf_name - 1, '@');
        dict_id = dict_id_make (sf_name-1, sf_name_len+1, DTYPE_VCF_FORMAT);
        SAFE_RESTORE;
    }

    return dict_id; 
}

static void vcf_seg_FORMAT (VBlockVCF *vb, ZipDataLineVCF *dl, STRp(fmt))
{
    ASSVCF0 (fmt_len >= 1, "missing or invalid FORMAT field");

    if (!vcf_num_samples) {
        seg_by_did_i_ex (VB, fmt, fmt_len, VCF_FORMAT, fmt_len + 1 /* \n */, NULL);
        return; // if we're not expecting any samples, no need to analyze the FORMAT field
    }

    dl->has_haplotype_data = (fmt[0] == 'G' && fmt[1] == 'T' && (fmt[2] == ':' || fmt_len==2)); // GT tag in FORMAT field - must always appear first per VCF spec (if it appears)
    dl->has_genotype_data  = (fmt_len > 2 || !dl->has_haplotype_data);

    bool possibly_rename = z_dual_coords && !vb->is_rejects_vb;
    bool possible_conditional_renaming = possibly_rename && (vb->last_index (VCF_oXSTRAND) || LO_IS_OK_SWITCH (last_ostatus));

    // case: FORMAT is the same as previous line - just use the same node_index, but only if there is no chance of conditional renaming
    if (fmt_len == vb->last_format.len && !memcmp (vb->last_format.data, fmt, fmt_len) && !possible_conditional_renaming) {
        seg_duplicate_last (VB, CTX(VCF_FORMAT), fmt_len + 1 /* \t or \n */);
        dl->format_node_i = (dl-1)->format_node_i;
        return;
    }   

    // save FORMAT field for potential duplication on the next line - only if we are certain there is no conditional renaming
    vb->last_format.len = 0; // reset 
    if (!possible_conditional_renaming)
        buf_add_moreS (vb, &vb->last_format, fmt, "last_format");

    str_split (fmt, fmt_len, 0, ':', sf_name, false);

    Container format_mapper = (Container){ 
        .drop_final_repeat_sep = true,
        .drop_final_item_sep   = true,
        .callback              = true,
        .filter_items          = true,
        .filter_repeats        = true,
        .repsep                = "\t"
    };

    ASSVCF (n_sf_names < MAX_FIELDS,
            "FORMAT field has too many subfields, the maximum allowed is %u: \"%.*s\"", MAX_FIELDS, STRf(fmt));

    con_set_nitems (format_mapper, n_sf_names);

    ContextPBlock ctxs;

    for (unsigned i=0; i < n_sf_names; i++) {

        DictId dict_id = vcf_seg_get_format_sf_dict_id (sf_names[i], sf_name_lens[i]);
        ctxs[i] = ctx_get_ctx_tag (vb, dict_id, sf_names[i], sf_name_lens[i]);

        if (possibly_rename) vcf_tags_add_tag (vb, ctxs[i], DTYPE_VCF_FORMAT, sf_names[i], sf_name_lens[i]);

        format_mapper.items[i] = (ContainerItem) { .dict_id = dict_id, .separator = {':'} };

        if (dict_id.num == _FORMAT_PS)
            format_mapper.items[i].separator[1] = CI1_ITEM_CB;

        // case: GL_to_PL:  FORMAT field snip is changed here to GL. Note: dict_id remains _FORMAT_GL.
        // so that vcf_seg_one_sample treats it as GL, and converts it to PL.
        if (dict_id.num == _FORMAT_GL && flag.GL_to_PL)
            *((char *)sf_names[i]) = 'P'; // change GL to GP (note: FORMAT changes and field changes, but still stored in dict_id=GL)

        // case: GP_to_PP - only relevant to VCF 4.3 where GP is probabilities and PP is Phred values (up to 4.2 GP was Phred values)
        if (dict_id.num == _FORMAT_GP && flag.GP_to_PP && vb->vcf_version >= VCF_v4_3)
            *((char *)sf_names[i]) = 'P'; // change GP to PP (note: FORMAT changes and field changes, but still stored in dict_id=GP)
    } 
    
    char renamed_data[n_sf_names * MAX_TAG_LEN];
    const char *renamed = renamed_data;
    int64_t renamed_len = possibly_rename ? vcf_tags_rename (vb, n_sf_names, ctxs, sf_names, sf_name_lens, NULL, renamed_data) : 0;

    // note: fmt_len and renamed_len need to be int64_t to avoid -Wstringop-overflow warning in gcc 10
    SNIP(6 + fmt_len + renamed_len); // maximum length in case of dual-snip with renaming

    if (renamed_len) {

        if (vb->line_coords == DC_LUFT) {
            vb->recon_size += renamed_len - fmt_len;
            SWAP (fmt, renamed);
            SWAP (fmt_len, renamed_len);
        }

        ASSERT0 (6 + fmt_len + renamed_len == snip_len, "String index out of range");

        snip[0] = snip[3+fmt_len] = SNIP_DUAL;
        snip[1] = snip[4+fmt_len] = SNIP_SPECIAL;
        snip[2] = snip[5+fmt_len] = VCF_SPECIAL_FORMAT;
        memcpy (&snip[3], fmt, fmt_len);
        memcpy (&snip[6+fmt_len], renamed, renamed_len);
    }

    // case: not tag renaming
    else {
        snip[0] = SNIP_SPECIAL;
        snip[1] = VCF_SPECIAL_FORMAT;
        memcpy (&snip[2], fmt, fmt_len);
        snip_len = 2 + fmt_len;
    }

    bool is_new;
    uint32_t node_index = seg_by_did_i_ex (VB, STRa(snip), VCF_FORMAT, fmt_len + 1 /* \t or \n */, &is_new);
    
    dl->format_node_i = node_index;

    if (is_new) {
        ASSVCF (node_index == vb->format_mapper_buf.len, 
                "node_index=%u different than vb->format_mapper_buf.len=%u", node_index, (uint32_t)vb->format_mapper_buf.len);

        buf_alloc (vb, &vb->format_mapper_buf, 1, 0, Container, CTX_GROWTH, "format_mapper_buf");
        buf_alloc (vb, &vb->format_contexts, 1, 0, ContextPBlock, CTX_GROWTH, "format_contexts");
        vb->format_contexts.len =vb->format_mapper_buf.len = vb->format_mapper_buf.len + 1;
    }    

    ContainerP con = ENT (Container, vb->format_mapper_buf, node_index);
    if (is_new || !con_nitems (*con)) { // assign if not already assigned (con->nitem=0 if format_mapper was already segged in another VB (so it is in ol_nodes), but not yet this VB) 
        *con = format_mapper; 
        memcpy (ENT (ContextPBlock, vb->format_contexts, node_index), ctxs, sizeof (ctxs));
    }
}

// returns length of lo_rejects before the copying
static uint32_t vcf_seg_copy_line_to_reject (VBlockVCF *vb, const char *field_start_line, uint32_t remaining_txt_len)
{
    const char *last = memchr (field_start_line, '\n', remaining_txt_len);
    ASSERT (last, "Line has no newline: %.*s", remaining_txt_len, field_start_line);

    uint32_t line_len = last - field_start_line + 1;
    buf_add_more (VB, &vb->lo_rejects[vb->line_coords-1], field_start_line, line_len, "lo_rejects");

    return line_len;
}

static inline LineCmpInfo vcf_seg_make_lci (ZipDataLineVCF *dl, bool is_luft)
{
    return (LineCmpInfo){ .chrom_wi    = dl->chrom[is_luft],
                          .start_pos   = dl->pos[is_luft],
                          .end_pos     = dl->pos[is_luft] + MAX_(0, dl->end_delta),
                          .tie_breaker = dl->tie_breaker };
}  

static inline void vcf_seg_evidence_of_unsorted (VBlockVCF *vb, ZipDataLineVCF *dl, DidIType chrom_did_i)
{
    bool is_luft = !!chrom_did_i; // 0 for primary, 1 for last

    if (!vb->is_unsorted[is_luft] && vb->line_i > 0)     
        vb->is_unsorted[is_luft] = linesorter_cmp (vcf_seg_make_lci (dl-1, is_luft), vcf_seg_make_lci (dl, is_luft)) > 0;
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

    // tie breaker is Adler32 of the Primary REF\tALT, except for Luft-only lines, where it is the Adler32 of the Luft REF\tALT
    dl->tie_breaker = libdeflate_adler32 (1, vb->main_refalt, vb->main_ref_len + 1 + vb->main_alt_len);
}

/* segment a VCF line into its fields:
   fields CHROM->FORMAT are "field" contexts
   all samples go into the SAMPLES context, which is a Container
   Haplotype and phase data are stored in a separate buffers + a SNIP_SPECIAL in the GT context  */
const char *vcf_seg_txt_line (VBlock *vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    ZipDataLineVCF *dl = DATA_LINE (vb->line_i);
    vb->sum_dp_this_line  = 0;
    vb->num_dps_this_line = 0;

    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    vb->line_coords = (txt_file->coords == DC_NONE) ? DC_PRIMARY : txt_file->coords;
    if (vb->reject_bytes) 
        vb->line_coords = OTHER_COORDS (vb->line_coords); // dual coordinates file - this line originates from ##primary_only/##luft_only as is in the opposite coordinates

    // allow un-segging of o* and PRIM/LUFT in case due to INFO or FORMAT ostatus changes from OK to REJECT
    int32_t save_recon_size=0; uint32_t line_len=0, save_lo_rejects_len=0; 
    unsigned save_txt_len_len = vb->num_contexts;
    uint64_t save_txt_len[save_txt_len_len];
    if (z_dual_coords) {
        save_lo_rejects_len = vb->lo_rejects[vb->line_coords-1].len; // copy line to lo_rejects now, because we are going to destroy it. We remove it later if its not a reject.
        line_len = vcf_seg_copy_line_to_reject (vb, field_start_line, remaining_txt_len);
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
    CTX(VCF_POS)->last_txt_index = ENTNUM (vb->txt_data, VCF_POS_str); // consumed by vcf_seg_FORMAT_PS
    vb->last_txt_len (VCF_POS) = VCF_POS_len;

    if (vb->line_coords == DC_PRIMARY) {
        PosType pos = dl->pos[0] = seg_pos_field (vb_, VCF_POS, VCF_POS, 0, '.', VCF_POS_str, VCF_POS_len, 0, VCF_POS_len+1);

        if (pos == 0 && !(*VCF_POS_str == '.' && VCF_POS_len == 1)) // POS == 0 - invalid value return from seg_pos_field
            WARN_ONCE ("FYI: invalid POS=%"PRId64" value in chrom=%.*s vb_i=%u vb_line_i=%"PRIu64": line will be compressed, but not indexed", 
                       pos, vb->chrom_name_len, vb->chrom_name, vb->vblock_i, vb->line_i);
                
        if (pos) random_access_update_pos (vb_, DC_PRIMARY, VCF_POS);
    }

    else { // LUFT
        dl->pos[1] = seg_pos_field (vb_, VCF_oPOS, VCF_oPOS, 0, '.', VCF_POS_str, VCF_POS_len, 0, VCF_POS_len);
        if (dl->pos[1]) random_access_update_pos (vb_, DC_LUFT, VCF_oPOS);
        CTX(vb->vb_coords==DC_LUFT ? VCF_oPOS : VCF_POS)->txt_len++; // account for the tab - in oPOS in the ##luft_only VB and in POS (on behalf on the primary POS) if this is a Dual-coord line (we will rollback accounting later if its not)
    }

    GET_NEXT_ITEM (VCF_ID);

    if (flag.add_line_numbers) 
        vcf_seg_add_line_number (vb, VCF_ID_len);
    else
        seg_id_field (vb_, CTX(VCF_ID), VCF_ID_str, VCF_ID_len, true);

    // REF + ALT 
    GET_NEXT_ITEM (VCF_REF);
    GET_NEXT_ITEM (VCF_ALT);

    // save REF and ALT (in primary or luft coordinates) to be used for INFO fields
    vb->main_refalt  = VCF_REF_str;
    vb->main_ref_len = VCF_REF_len;
    vb->main_alt_len = VCF_ALT_len;
 
    // note: we treat REF+\t+ALT as a single field because REF and ALT are highly corrected, in the case of SNPs:
    // e.g. GG has a probability of 0 and GC has a higher probability than GA.
    vcf_refalt_seg_main_ref_alt (vb, VCF_REF_str, VCF_REF_len, VCF_ALT_str, VCF_ALT_len);
    
    SEG_NEXT_ITEM (VCF_QUAL);
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
            
            const char *backup_luft_samples = NULL; uint32_t backup_luft_samples_len=0;
            if (vb->line_coords == DC_LUFT) {
                 backup_luft_samples = ENT (char, vb->lo_rejects[DC_LUFT-1], save_lo_rejects_len + (next_field - field_start_line));
                 backup_luft_samples_len = line_len - (next_field - field_start_line);
            }

            // seg all samples. note that this is destructive: 1. Samples are first lift-back to PRIMARY if needed ; 2. FORMAT/GT overwrite txt_data 
            next_field = vcf_seg_samples (vb, dl, &len, (char*)next_field, has_13, backup_luft_samples, backup_luft_samples_len); 
        }
        else 
            seg_by_did_i (VB, NULL, 0, VCF_SAMPLES, 0); // case no samples: WORD_INDEX_MISSING
    }

    // case no format or samples
    else {
        seg_by_did_i (VB, NULL, 0, VCF_FORMAT, 0); 
        seg_by_did_i (VB, NULL, 0, VCF_SAMPLES, 0);
    }

    // Adds DVCF items according to ostatus, finalizes INFO/SF and segs the INFO container
    vcf_finalize_seg_info (vb);

    // calculate tie-breaker for sorting - do after INFO before segging the samples as GT data can overwrite REF ALT
    if (flag.sort) vcf_seg_assign_tie_breaker (vb, dl);

    if (!has_samples && vcf_num_samples)
        WARN_ONCE ("FYI: variant CHROM=%.*s POS=%"PRId64" has no samples", vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS));

    if (z_dual_coords) {
        // note: we don't seg Coord for non-DC files, despite being in TOPLEVEL. That's ok, it will be treated as
        // as "all_the_same" and have word_index=0 == NONE for all lines
        bool lo_ok = LO_IS_OK (last_ostatus);

        Coords reconstructable_coords = lo_ok ? DC_BOTH : vb->line_coords;
        const char *name = coords_name (reconstructable_coords); 
        seg_by_did_i (VB, name, strlen (name), VCF_COORDS, 0); // 0 as its not in the txt data
        
        if (reconstructable_coords & DC_PRIMARY) vb->recon_num_lines++;   // PRIMARY or BOTH
        if (reconstructable_coords & DC_LUFT) vb->recon_num_lines_luft++; // LUFT or BOTH
        
        // if line was NOT rejected (the default, if not dual coordinates), we can delete the text from lo_rejects
        if (lo_ok)
            vb->lo_rejects[vb->line_coords-1].len = save_lo_rejects_len;

        // in a dual-coordinate VB (i.e. the main VCF data lines), single-coordinate lines won't be displayed in the opposite reconstruction
        else {
            if (vb->line_coords == DC_PRIMARY) vb->recon_size_luft = save_recon_size - line_len;
            else  /* DC_LUFT  */               vb->recon_size      = save_recon_size - line_len;
        
            // unaccount for this line, if its a Luft-only line in a dual-coordinate variant (i.e. main VCF data line) as it won't appear in the default reconstruction
            if (vb->vb_coords == DC_BOTH && vb->line_coords == DC_LUFT)
                for (unsigned i=0; i < vb->num_contexts; i++) 
                    CTX(i)->txt_len = (i < save_txt_len_len ? save_txt_len[i] : 0);
        }
        
        // in case of a reject line - update reject_bytes to indicate its consumption
        if (txt_file->coords && txt_file->coords != vb->line_coords)
            vb->reject_bytes -= next_field - field_start_line;
    } 

    SEG_EOL (VCF_EOL, false);

    // test if still sorted
    if (flag.sort) {
        vcf_seg_evidence_of_unsorted (vb, dl, VCF_CHROM);
        if (z_dual_coords) vcf_seg_evidence_of_unsorted (vb, dl, VCF_oCHROM);
    }

    return next_field;
}
