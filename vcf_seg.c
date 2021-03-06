// ------------------------------------------------------------------
//   vcf_seg.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

static MiniContainer line_number_container = {};

// called by main thread after reading the header
void vcf_zip_initialize (void)
{
    vcf_lo_zip_initialize ();
    vcf_seg_samples_initialize ();

    // container just for adding a prefix to the delta-encoded line number (the container is all_the_same)
    if (flag.add_line_numbers && !line_number_container.repeats) { // possibly already initialized by previous files
        line_number_container = (MiniContainer) {
            .repeats   = 1,
            .nitems_lo = 1,
            .items     = { { .dict_id = (DictId)dict_id_fields[VCF_LINE_NUM] } }
        };
    }
}

void vcf_zip_read_one_vb (VBlockP vb)
{
    // set vcf_version in VB, since the the vcf_version in vcf_header might change as we might be reading the next txt file
    ((VBlockVCFP)vb)->vcf_version = vcf_header_get_version();

    // in case we're replacing ID with the line number
    if (flag.add_line_numbers) {
        vb->first_line = txt_file->num_lines + 1; // this is normally not used in ZIP
        txt_file->num_lines += str_count_char (vb->txt_data.data, vb->txt_data.len, '\n');  // update here instead of in zip_update_txt_counters;
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
    CTX(VCF_REFALT)->  keep_snip   = true; // set ctx->last_snip after evaluating
    CTX(VCF_oXSTRAND)->keep_snip   = true;

    CTX(VCF_oSTATUS)-> flags.store = STORE_INDEX;
    CTX(VCF_COORDS)->  flags.store = STORE_INDEX;
    CTX(VCF_oXSTRAND)->flags.store = STORE_INDEX;
    CTX(VCF_CHROM)->   flags.store = STORE_INDEX; // since v12
    CTX(VCF_oCHROM)->  flags.store = STORE_INDEX; // used by regions_is_site_included
    CTX(VCF_POS)->     flags.store = STORE_INT;   // since v12
    CTX(VCF_oPOS)->    flags.store = STORE_INT;   // used by vcf_piz_luft_END

    // counts sections
    CTX(VCF_CHROM)->   counts_section = true;
    CTX(VCF_oCHROM)->  counts_section = true;

    if (z_dual_coords) {
        CTX(VCF_oSTATUS)->counts_section = true;
        CTX(VCF_COORDS)-> counts_section = true;
    }

    // consolidate stats
    CTX(VCF_oREFALT)->st_did_i = CTX(VCF_LIFT_REF)->st_did_i = VCF_REFALT;
    CTX(VCF_oPOS)->st_did_i = CTX(VCF_COPYPOS)->st_did_i = VCF_POS;
    CTX(VCF_oCHROM)->st_did_i = VCF_CHROM;
    CTX(VCF_COPYSTAT)->st_did_i = CTX(VCF_oSTATUS)->st_did_i = VCF_COORDS; // dual-coordinate non-reconstruted stuff goes here
    
    Context *gt_gtx   = ctx_get_ctx (vb, dict_id_FORMAT_GT);
    gt_gtx->no_stons  = true; // we store the GT matrix in local, so cannot accomodate singletons
    vb->ht_matrix_ctx = ctx_get_ctx (vb, dict_id_FORMAT_GT_HT);

    // room for already existing FORMATs from previous VBs
    vb->format_mapper_buf.len = CTX(VCF_FORMAT)->ol_nodes.len;
    buf_alloc_zero (vb, &vb->format_mapper_buf, 0, vb->format_mapper_buf.len, Container, 1.2, "format_mapper_buf");
    
    // create additional contexts as needed for compressing FORMAT/GT - must be done before merge
    if (vcf_num_samples) {
        Context *runs_ctx = ctx_get_ctx (vb, dict_id_PBWT_RUNS); // must be created before FGRC so it is emitted in the file in this order
        Context *fgrc_ctx = ctx_get_ctx (vb, dict_id_PBWT_FGRC);
        codec_pbwt_seg_init (vb_, runs_ctx, fgrc_ctx, gt_gtx->did_i);
    }

    if (flag.add_line_numbers) {
        // create a b250 and dict entry for VCF_LINE_NUM, VCF_ID - these become "all_the_same" so no need to seg them explicitly hereinafter        
        #define LN_PREFIX_LEN 3 // "LN="
        static const char prefix[] = { CON_PREFIX_SEP, // has prefix 
                                       CON_PREFIX_SEP, // end of (empty) container-wide prefix
                                       'L', 'N', '=', CON_PREFIX_SEP };  // NOTE: if changing prefix, update LN_PREFIX_LEN

        container_seg_by_ctx (vb, CTX(VCF_ID), (Container *)&line_number_container, prefix, sizeof prefix, 0); 
        ctx_decrement_count (vb_, CTX(VCF_ID), 0);
        CTX(VCF_ID)->no_stons = true;
    }

    // evaluate oSTATUS and COORDS, XSTRAND snips in order, as we rely on their indices being identical to the order of these arrays
    for (int i=0; i < NUM_LO_STATUSES; i++) {
        WordIndex node_index = ctx_evaluate_snip_seg (vb_, CTX(VCF_oSTATUS), dvcf_status_names[i], strlen (dvcf_status_names[i]), NULL);
        ctx_decrement_count (vb_, CTX(VCF_oSTATUS), node_index);
    }

    for (int i=0; i < NUM_COORDS; i++) {
        WordIndex node_index = ctx_evaluate_snip_seg (vb_, CTX(VCF_COORDS), coords_name(i), strlen (coords_name(i)), NULL);
        ctx_decrement_count (vb_, CTX(VCF_COORDS), node_index);
    }

    ctx_evaluate_snip_seg (vb_, CTX(VCF_oXSTRAND), "-", 1, NULL);
    ctx_decrement_count (vb_, CTX(VCF_oXSTRAND), 0);

    ctx_evaluate_snip_seg (vb_, CTX(VCF_oXSTRAND), "X", 1, NULL);
    ctx_decrement_count (vb_, CTX(VCF_oXSTRAND), 1);

    // create a b250 and dict entry for LIFT_REF, COPYSTAT and COPYPOS - these become "all_the_same" so no need to seg them explicitly hereinafter
    seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_LIFT_REF }), 2, VCF_LIFT_REF, 0); 
    ctx_decrement_count (vb_, CTX(VCF_LIFT_REF), 0);

    seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPYSTAT }), 2, VCF_COPYSTAT, 0); 
    ctx_decrement_count (vb_, CTX(VCF_COPYSTAT), 0);

    seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPYPOS }), 2, VCF_COPYPOS, 0);
    ctx_decrement_count (vb_, CTX(VCF_COPYPOS), 0);
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
        .items        = { { .dict_id = (DictId)dict_id_fields[VCF_COORDS],  .seperator = "\t" }, // suppressed by vcf_piz_filter unless --show-ostatus                                   
                          { .dict_id = (DictId)dict_id_fields[VCF_oSTATUS], .seperator = "\t" }, // suppressed by vcf_piz_filter unless --show-ostatus                                   
                          { .dict_id = (DictId)dict_id_fields[VCF_CHROM],   .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_POS],     .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_ID],      .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_REFALT],  .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_QUAL],    .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_FILTER],  .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_INFO],    .seperator = "\t" }, // in dual-coordinates, contains INFO/LIFTOVER or INFO/REJTOVER that reconstructs oCHROM, oPOS, oREF, oXSTRAND
                          { .dict_id = (DictId)dict_id_fields[VCF_FORMAT],  .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_SAMPLES], .seperator = ""   },
                          { .dict_id = (DictId)dict_id_fields[VCF_EOL],     .seperator = ""   } },
    };

    Context *ctx = CTX(VCF_TOPLEVEL);

    if (vb->vb_coords == DC_BOTH || !z_dual_coords)
        container_seg_by_ctx (vb_, ctx, (ContainerP)&top_level, 0, 0, 0); 

    // when processing the rejects file containing variants that are primary-only, we add a "##primary_only=" prefix to 
    // first item of each line, so that it reconstructs as part of the VCF header 
    else if (vb->vb_coords == DC_PRIMARY) { // primary-only variants 
        static const char primary_only_prefix[] = CON_PREFIX_SEP_ CON_PREFIX_SEP_ HK_PRIM_ONLY CON_PREFIX_SEP_;
        container_seg_by_ctx (vb_, ctx, (ContainerP)&top_level, primary_only_prefix, strlen (primary_only_prefix), 0);
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
        .items        = { { .dict_id = (DictId)dict_id_fields[VCF_COORDS],  .seperator = "\t" }, // suppressed by vcf_piz_filter unless --show-ostatus                                   
                          { .dict_id = (DictId)dict_id_fields[VCF_oSTATUS], .seperator = "\t" }, // suppressed by vcf_piz_filter unless --show-ostatus                                   
                          { .dict_id = (DictId)dict_id_fields[VCF_oCHROM],  .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_oPOS],    .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_ID],      .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_oREFALT], .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_QUAL],    .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_FILTER],  .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_POS],     .seperator = { CI_TRANS_NOR } }, // consume POS before INFO, in case we have INFO/END
                          { .dict_id = (DictId)dict_id_fields[VCF_INFO],    .seperator = "\t" }, // in dual-coordinates, contains INFO/LIFTOVER or INFO/REJTOVER that reconstructs oCHROM, oPOS, oREF, oXSTRAND
                          { .dict_id = (DictId)dict_id_fields[VCF_FORMAT],  .seperator = "\t" },
                          { .dict_id = (DictId)dict_id_fields[VCF_SAMPLES], .seperator = ""   },
                          { .dict_id = (DictId)dict_id_fields[VCF_EOL],     .seperator = ""   } }
    };

    if (vb->vb_coords == DC_BOTH)
        container_seg_by_ctx (vb_, CTX(VCF_TOPLUFT), (ContainerP)&top_luft, 0, 0, 0);

    // similarly, when processing the rejects file containing variants that are luft-only, we add a "##luft_only=" prefix
    else if (vb->vb_coords == DC_LUFT) { // luft-only variants 
        static const char luft_only_prefix[] = CON_PREFIX_SEP_ CON_PREFIX_SEP_ HK_LUFT_ONLY CON_PREFIX_SEP_;
        container_seg_by_ctx (vb_, CTX(VCF_TOPLUFT), (ContainerP)&top_luft, luft_only_prefix, strlen (luft_only_prefix), (sizeof HK_LUFT_ONLY - 1) * vb->lines.len);
        vb->recon_size_luft += (sizeof HK_LUFT_ONLY - 1) * vb->lines.len; // when reconstructing luft-only rejects, we also reconstruct a prefix for each line
        // note: there is no equivalent of ctx->txt_len for Luft coordinates
    }

    // consolidate DVCF info fields to VCF_COORDS (just the container, not the values)
    if ((ctx = ctx_get_existing_ctx (vb, dict_id_INFO_PRIM))) ctx->st_did_i = VCF_COORDS;
    if ((ctx = ctx_get_existing_ctx (vb, dict_id_INFO_PREJ))) ctx->st_did_i = VCF_COORDS;
    if ((ctx = ctx_get_existing_ctx (vb, dict_id_INFO_LUFT))) ctx->st_did_i = VCF_COORDS;
    if ((ctx = ctx_get_existing_ctx (vb, dict_id_INFO_LREJ))) ctx->st_did_i = VCF_COORDS;
}

bool vcf_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return 
        dict_id.num == dict_id_fields[VCF_TOPLEVEL] ||
        dict_id.num == dict_id_fields[VCF_TOPLUFT]  ||
        dict_id.num == dict_id_fields[VCF_CHROM]    ||
        dict_id.num == dict_id_fields[VCF_oCHROM]   ||
        dict_id.num == dict_id_fields[VCF_FORMAT]   ||
        dict_id.num == dict_id_fields[VCF_INFO]     ||
        dict_id.num == dict_id_fields[VCF_REFALT]   ||
        dict_id.num == dict_id_fields[VCF_oREFALT]  ||
        dict_id.num == dict_id_fields[VCF_FILTER]   ||
        dict_id.num == dict_id_fields[VCF_EOL]      ||
        dict_id.num == dict_id_fields[VCF_SAMPLES]  ||
        dict_id.num == dict_id_fields[VCF_oCHROM]   ||
        dict_id.num == dict_id_fields[VCF_oXSTRAND] ||
        dict_id.num == dict_id_fields[VCF_oSTATUS]  ||
        dict_id.num == dict_id_fields[VCF_COORDS]   ||
        dict_id.num == dict_id_fields[VCF_LIFT_REF] ||
        dict_id.num == dict_id_INFO_AC              ||
        dict_id.num == dict_id_INFO_AF              ||
        dict_id.num == dict_id_INFO_AN              ||
        dict_id.num == dict_id_INFO_DP              ||
        dict_id.num == dict_id_INFO_AA              || // stored as a SPECIAL snip
        dict_id.num == dict_id_INFO_MLEAC           ||
        dict_id.num == dict_id_INFO_MLEAF           ||
        dict_id.num == dict_id_INFO_LDAF            ||
        dict_id.num == dict_id_INFO_MQ0             ||
        dict_id.num == dict_id_INFO_LUFT            ||
        dict_id.num == dict_id_INFO_PRIM            ||
        dict_id.num == dict_id_INFO_LREJ            ||
        dict_id.num == dict_id_INFO_PREJ            ||

        // INFO/ AC_* AN_* AF_* and ???_AF are small
        ((dict_id.id[0] == ('A' | 0xc0)) && (dict_id.id[1] == 'C' || dict_id.id[1] == 'F' || dict_id.id[1] == 'N') && dict_id.id[2] == '_') ||
        (dict_id_is_vcf_info_sf (dict_id) && dict_id.id[3] == '_' && dict_id.id[4] == 'A' && dict_id.id[5] == 'F' && !dict_id.id[6]);
}


// traverses the FORMAT field, gets ID of subfield, and moves to the next subfield
static DictId vcf_seg_get_format_subfield (const char **str, uint32_t *len) // remaining length of line 
{
    unsigned i=0; for (; i < *len && (*str)[i] != ':' && (*str)[i] != '\t' && (*str)[i] != '\n'; i++);
    
    DictId dict_id;
    // case: normal field - starts with a letter or another character in the range
    if ((*str)[0] >= 64 && (*str)[0] <= 127) 
        dict_id = dict_id_make (*str, i, DTYPE_VCF_FORMAT);
    
    // case: unusual field - starts with an out-of-range character, eg a digit - prefix with @ so its a legal FORMAT dict_id
    else {
        SAFE_ASSIGN (*str - 1, '@');
        dict_id = dict_id_make (*str-1, i+1, DTYPE_VCF_FORMAT);
        SAFE_RESTORE;
    }

    *str += i+1;
    *len -= i+1;
    return dict_id; 
}


// assign contexts to a format mapper, if not already assigned
static void vcf_seg_format_get_contexts (VBlockVCF *vb, ZipDataLineVCF *dl)
{
    ContainerP format_mapper = ENT (Container, vb->format_mapper_buf, dl->format_node_i);

    // an array of blocks of MAX_FIELDS of ContextP. The number of such blocks is format_contexts.len
    buf_alloc_zero (vb, &vb->format_contexts, 0, vb->format_mapper_buf.len * MAX_FIELDS, ContextP, CTX_GROWTH, "format_contexts");
    ContextP *ctxs = ENT (ContextP, vb->format_contexts, dl->format_node_i * MAX_FIELDS);

    for (unsigned i=0; i < con_nitems (*format_mapper); i++) 
        if (!ctxs[i])
            ctxs[i] = ctx_get_ctx (vb, format_mapper->items[i].dict_id);
}

static void vcf_seg_format_field (VBlockVCF *vb, ZipDataLineVCF *dl, const char *field_start, int field_len)
{
    const char *str = field_start;
    int len = field_len;

    if (!vcf_num_samples) {
        seg_by_did_i_ex (vb, field_start, field_len, VCF_FORMAT, field_len + 1 /* \n */, NULL);
        return; // if we're not expecting any samples, no need to analyze the FORMAT field
    }

    ASSVCF0 (field_len >= 2, "missing or invalid FORMAT field");

    Container format_mapper = (Container){ 
        .drop_final_repeat_sep = true,
        .drop_final_item_sep   = true,
        .callback              = true,
        .filter_items          = true,
        .filter_repeats        = true,
        .repsep                = "\t"
    };

    dl->has_haplotype_data = (str[0] == 'G' && str[1] == 'T' && (str[2] == ':' || field_len==2)); // GT tag in FORMAT field - must always appear first per VCF spec (if it appears)
    dl->has_genotype_data  = (field_len > 2 || (!dl->has_haplotype_data && field_len > 0));

    bool last_item = false;
    do {
        ASSVCF (con_nitems (format_mapper) < MAX_FIELDS,
                "FORMAT field has too many subfields, the maximum allowed is %u: \"%.*s\"",  
                MAX_FIELDS, field_len, field_start);

        DictId dict_id = vcf_seg_get_format_subfield (&str, (unsigned *)&len);
        last_item = (str[-1] == '\t' || str[-1] == '\n');

        format_mapper.items[con_nitems (format_mapper)] = (ContainerItem) {
            .dict_id   = dict_id,
            .seperator = {':'}
        };
        con_inc_nitems (format_mapper);

        ASSVCF (dict_id_is_vcf_format_sf (dict_id),
                "string %.*s in the FORMAT field \"%.*s\" is not a legal subfield", 
                DICT_ID_LEN, dict_id.id, field_len, field_start);

        // case: GL_to_PL:  FORMAT field snip is changed here to GL. Note: dict_id remains dict_id_FORMAT_GL.
        // so that vcf_seg_one_sample treats it as GL, and converts it to PL.
        if (dict_id.num == dict_id_FORMAT_GL && flag.GL_to_PL)
            ((char *)str)[-3] = 'P'; // change GL to GP (note: FORMAT changes and field changes, but still stored in dict_id=GL)

        // case: optimize_GP - only relevant to VCF 4.3 where GP is probabilities and PP is Phred values (up to 4.2 GP was Phred values)
        if (dict_id.num == dict_id_FORMAT_GP && flag.GP_to_PP && vb->vcf_version >= VCF_v4_3)
            ((char *)str)[-3] = 'P'; // change GP to PP (note: FORMAT changes and field changes, but still stored in dict_id=GP)
    } 

    while (!last_item && len > 0);
    
    bool is_new;
    char snip[field_len+2];
    snip[0] = SNIP_SPECIAL;
    snip[1] = VCF_SPECIAL_FORMAT;
    memcpy (&snip[2], field_start, field_len);

    uint32_t node_index = seg_by_did_i_ex (vb, snip, field_len+2, VCF_FORMAT, field_len + 1 /* \t or \n */, &is_new);

    dl->format_node_i = node_index;

    if (is_new) {
        ASSVCF (node_index == vb->format_mapper_buf.len, 
                "node_index=%u different than vb->format_mapper_buf.len=%u", node_index, (uint32_t)vb->format_mapper_buf.len);

        vb->format_mapper_buf.len++;
        buf_alloc (vb, &vb->format_mapper_buf, 0, vb->format_mapper_buf.len, Container, 2, "format_mapper_buf");
    }    

    ContainerP con = ENT (Container, vb->format_mapper_buf, node_index);
    if (is_new || !con_nitems (*con)) // assign if not already assigned. 
        *con = format_mapper; 
}

// returns length of lo_rejects before the copying
static uint32_t vcf_seg_copy_line_to_reject (VBlockVCF *vb, const char *field_start_line, uint32_t remaining_txt_len)
{
    const char *last = memchr (field_start_line, '\n', remaining_txt_len);
    ASSERT (last, "Line has no newline: %.*s", remaining_txt_len, field_start_line);

    uint32_t line_len = last - field_start_line + 1;
    buf_add_more (vb, &vb->lo_rejects[vb->line_coords-1], field_start_line, line_len, "lo_rejects");

    return line_len;
}

static inline void vcf_seg_assign_dl_sorted (VBlockVCF *vb, ZipDataLineVCF *dl, DidIType chrom_did_i)
{
    bool is_luft = !!chrom_did_i; // 0 for primary, 1 for last

    if (!vb->is_unsorted[is_luft] && vb->line_i > 0) {
        // cases where we have evidence this VB is NOT sorted:

        // 1. pos is out of order
        if( (dl->chrom_index[is_luft] == (dl-1)->chrom_index[is_luft] && dl->pos[is_luft] < (dl-1)->pos[is_luft])

        // 2. rejected lift-over 
        ||  (is_luft && dl->chrom_index[is_luft] == WORD_INDEX_NONE)

        // 3. chrom os out of order (if chroms are in header, then this is by the header, otherwise by first appearance in file)
        ||  (dl->chrom_index[is_luft] < (dl-1)->chrom_index[is_luft]))
           
            vb->is_unsorted[is_luft] = true;
    }
}

static void vcf_seg_add_line_number (VBlockVCFP vb, unsigned VCF_ID_len)
{
    char line_num[20];
    unsigned line_num_len = str_int (vb->first_line + vb->line_i, line_num);

    seg_pos_field ((VBlockP)vb, VCF_LINE_NUM, VCF_LINE_NUM, SPF_ZERO_IS_BAD, 0, line_num, line_num_len, 0, line_num_len + 1 + LN_PREFIX_LEN); // +1 for tab +3 for "LN="
    
    int shrinkage = (int)VCF_ID_len - line_num_len - LN_PREFIX_LEN;
    vb->recon_size -= shrinkage;
    vb->recon_size_luft -= shrinkage;
}

/* segment a VCF line into its fields:
   fields CHROM->FORMAT are "field" contexts
   all samples go into the SAMPLES context, which is a Container
   Haplotype and phase data are stored in a separate buffers + a SNIP_SPECIAL in the GT context  */
const char *vcf_seg_txt_line (VBlock *vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    ZipDataLineVCF *dl = DATA_LINE (vb->line_i);

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
        dl->chrom_index[0] = seg_chrom_field (vb_, VCF_CHROM_str, VCF_CHROM_len);
    
    else { // LUFT
        dl->chrom_index[1] = seg_by_did_i (vb_, VCF_CHROM_str, VCF_CHROM_len, VCF_oCHROM, VCF_CHROM_len); // we will add_bytes of the CHROM field (not oCHROM) as genounzip reconstructs the PRIMARY
        random_access_update_chrom (vb_, DC_LUFT, dl->chrom_index[1], VCF_CHROM_str, VCF_CHROM_len); // also sets vb->chrom_name
        CTX(vb->vb_coords==DC_LUFT ? VCF_oCHROM : VCF_CHROM)->txt_len++; // account for the tab - in oCHROM in the ##luft_only VB and in CHROM (on behalf on the primary CHROM) if this is a Dual-coord line (we will rollback accounting later if its not)
    }

    GET_NEXT_ITEM (VCF_POS);
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
        seg_id_field (vb_, dict_id_fields[VCF_ID], VCF_ID_str, VCF_ID_len, true);

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

    bool has_samples = true;
    if (separator != '\n') { // has a FORMAT field

        // FORMAT
        if (vcf_num_samples) {
            GET_MAYBE_LAST_ITEM (VCF_FORMAT); // possibly no samples this line
        } else {
            GET_LAST_ITEM (VCF_FORMAT);
        }

        vcf_seg_format_field (vb, dl, field_start, field_len);

        if (separator != '\n') { // has samples

            ASSVCF0 (dl->has_genotype_data || dl->has_haplotype_data, "expecting line to end as it has no sample data, but it has not");
            
            vcf_seg_format_get_contexts (vb, dl);

            const char *backup_luft_samples = NULL; uint32_t backup_luft_samples_len=0;
            if (vb->line_coords == DC_LUFT) {
                 backup_luft_samples = ENT (char, vb->lo_rejects[DC_LUFT-1], save_lo_rejects_len + (next_field - field_start_line));
                 backup_luft_samples_len = line_len - (next_field - field_start_line);
            }

            // seg all samples. note that this is destructive: 1. Samples are first lift-back to PRIMARY if needed ; 2. FORMAT/GT overwrite txt_data 
            next_field = vcf_seg_samples (vb, dl, &len, (char*)next_field, has_13, backup_luft_samples, backup_luft_samples_len); 
        }
        else {
            has_samples = false;
            seg_by_did_i (vb, NULL, 0, VCF_SAMPLES, 0); // case no samples: WORD_INDEX_MISSING_SF
        }
    }

    // case no format or samples
    else {
        has_samples = false;
        seg_by_did_i (vb, NULL, 0, VCF_FORMAT, 0); 
        seg_by_did_i (vb, NULL, 0, VCF_SAMPLES, 0);
    }

    // Adds DVCF items according to ostatus, finalizes INFO/SF and segs the INFO container
    vcf_finalize_seg_info (vb);

    if (!has_samples && vcf_num_samples)
        WARN_ONCE ("FYI: variant CHROM=%.*s POS=%"PRId64" has no samples", vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS));

    if (z_dual_coords) {
        // note: we don't seg Coord for non-DC files, despite being in TOPLEVEL. That's ok, it will be treated as
        // as "all_the_same" and have word_index=0 == NONE for all lines
        bool lo_ok = LO_IS_OK (last_ostatus);

        Coords reconstructable_coords = lo_ok ? DC_BOTH : vb->line_coords;
        const char *name = coords_name (reconstructable_coords); 
        seg_by_did_i (vb, name, strlen (name), VCF_COORDS, 0); // 0 as its not in the txt data
        
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
        vcf_seg_assign_dl_sorted (vb, dl, VCF_CHROM);
        if (z_dual_coords) vcf_seg_assign_dl_sorted (vb, dl, VCF_oCHROM);
    }

    return next_field;
}
