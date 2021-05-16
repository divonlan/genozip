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

const char *coords_names[4] = COORDS_NAMES;

static void vcf_seg_complete_missing_lines (VBlockVCF *vb);

// called by main thread after reading the header
void vcf_zip_initialize (void)
{
    vcf_lo_zip_initialize ();
    vcf_format_zip_initialize ();
}

void vcf_zip_after_compute (VBlockP vb)
{
    if (z_dual_coords && !flag.rejects_coord) { // normal file, not zipping rejects files
        if (vb->lo_rejects[DC_PRIMARY-1].len) vcf_lo_append_rejects_file (vb, DC_PRIMARY);
        if (vb->lo_rejects[DC_LUFT-1].len)    vcf_lo_append_rejects_file (vb, DC_LUFT);
    }
}   

// called by Compute threadfrom seg_all_data_lines
void vcf_seg_initialize (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    vb->contexts[VCF_CHROM]   .no_stons    = true; // needs b250 node_index for random access and reconstrution plan
    vb->contexts[VCF_oCHROM]  .no_stons    = true; // same
    vb->contexts[VCF_FORMAT]  .no_stons    = true;
    vb->contexts[VCF_INFO]    .no_stons    = true;
    vb->contexts[VCF_oSTATUS] .no_stons    = true;
    vb->contexts[VCF_COORDS]  .no_stons    = true;
    vb->contexts[VCF_TOPLEVEL].no_stons    = true; 
    vb->contexts[VCF_TOPLUFT] .no_stons    = true; 
    vb->contexts[VCF_LIFT_REF].no_stons    = true; 
    vb->contexts[VCF_COPY_POS].no_stons    = true; 
    vb->contexts[VCF_oXSTRAND].no_stons    = true; // keep in b250 so it can be eliminated as all_the_same
    vb->contexts[VCF_oCHROM]  .no_vb1_sort = true; // indices need to remain as in the Chain file
    vb->contexts[VCF_oSTATUS] .no_vb1_sort = true; // indices need to remaining matching to LiftOverStatus
    vb->contexts[VCF_REFALT]  .keep_snip   = true; // set ctx->last_snip after evaluating
    vb->contexts[VCF_oXSTRAND].keep_snip   = true;

    vb->contexts[VCF_oSTATUS] .flags.store = STORE_INDEX;
    vb->contexts[VCF_COORDS]  .flags.store = STORE_INDEX;
    vb->contexts[VCF_CHROM]   .flags.store = STORE_INDEX; // since v12
    vb->contexts[VCF_oCHROM]  .flags.store = STORE_INDEX; // used by regions_is_site_included
    vb->contexts[VCF_POS]     .flags.store = STORE_INT;   // since v12
    vb->contexts[VCF_oPOS]    .flags.store = STORE_INT;   // used by vcf_piz_luft_END

    // consolidate stats
    vb->contexts[VCF_oREFALT].st_did_i = vb->contexts[VCF_LIFT_REF].st_did_i = VCF_REFALT;
    vb->contexts[VCF_oPOS].st_did_i = VCF_POS;
    vb->contexts[VCF_oCHROM].st_did_i = VCF_CHROM;
    vb->contexts[VCF_COPYSTAT].st_did_i = VCF_oSTATUS;
    
    Context *gt_gtx   = ctx_get_ctx (vb, dict_id_FORMAT_GT);
    gt_gtx->no_stons  = true; // we store the GT matrix in local, so cannot accomodate singletons
    vb->ht_matrix_ctx = ctx_get_ctx (vb, dict_id_FORMAT_GT_HT);

    // room for already existing FORMATs from previous VBs
    vb->format_mapper_buf.len = vb->contexts[VCF_FORMAT].ol_nodes.len;
    buf_alloc (vb, &vb->format_mapper_buf, 0, vb->format_mapper_buf.len, Container, 1.2, "format_mapper_buf");
    buf_zero (&vb->format_mapper_buf);

    // create additional contexts as needed for compressing FORMAT/GT - must be done before merge
    if (vcf_num_samples) {
        Context *runs_ctx = ctx_get_ctx (vb, dict_id_PBWT_RUNS); // must be created before FGRC so it is emitted in the file in this order
        Context *fgrc_ctx = ctx_get_ctx (vb, dict_id_PBWT_FGRC);
        codec_pbwt_seg_init (vb_, runs_ctx, fgrc_ctx, gt_gtx->did_i);
    }

    // evaluate oSTATUS and COORDS snips in order, as we rely on their indices being identical to the order of these arrays
    for (int i=0; i < NUM_LO_STATUSES; i++) {
        WordIndex node_index = ctx_evaluate_snip_seg (vb_, &vb->contexts[VCF_oSTATUS], liftover_status_names[i], strlen (liftover_status_names[i]), NULL);
        ctx_decrement_count (vb_, &vb->contexts[VCF_oSTATUS], node_index);
    }

    for (int i=0; i < NUM_COORDS; i++) {
        WordIndex node_index = ctx_evaluate_snip_seg (vb_, &vb->contexts[VCF_COORDS], coords_names[i], strlen (coords_names[i]), NULL);
        ctx_decrement_count (vb_, &vb->contexts[VCF_COORDS], node_index);
    }
    
    // create dict entry for LIFT_REF, COPYSTAT and COPY_POS - these become "all_the_same" so no need to seg them explicitly
    seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_LIFT_REF }), 2, VCF_LIFT_REF, 0); 
    ctx_decrement_count (vb_, &vb->contexts[VCF_LIFT_REF], 0);

    seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPYSTAT }), 2, VCF_COPYSTAT, 0); 
    ctx_decrement_count (vb_, &vb->contexts[VCF_COPYSTAT], 0);

    seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPY_POS }), 2, VCF_COPY_POS, 0);
    ctx_decrement_count (vb_, &vb->contexts[VCF_COPY_POS], 0);

}             

void vcf_seg_finalize (VBlockP vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    
    if (vb->ht_matrix_ctx) 
        vcf_seg_complete_missing_lines (vb);

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

    if (vb->vb_coords == DC_BOTH || !z_dual_coords)
        container_seg_by_ctx (vb_, &vb->contexts[VCF_TOPLEVEL], (ContainerP)&top_level, 0, 0, 0); 

    // when processing the rejects file containing variants that are primary-only, we add a "##primary_only=" prefix to 
    // first item of each line, so that it reconstructs as part of the VCF header 
    else if (vb->vb_coords == DC_PRIMARY) { // primary-only variants 
        static const char primary_only_prefix[] = CON_PREFIX_SEP_ CON_PREFIX_SEP_ HK_PRIM_ONLY CON_PREFIX_SEP_;
        container_seg_by_ctx (vb_, &vb->contexts[VCF_TOPLEVEL], (ContainerP)&top_level, primary_only_prefix, strlen (primary_only_prefix), 0);
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
        container_seg_by_ctx (vb_, &vb->contexts[VCF_TOPLUFT], (ContainerP)&top_luft, 0, 0, 0);

    // similarly, when processing the rejects file containing variants that are luft-only, we add a "##luft_only=" prefix
    else if (vb->vb_coords == DC_LUFT) { // luft-only variants 
        static const char luft_only_prefix[] = CON_PREFIX_SEP_ CON_PREFIX_SEP_ HK_LUFT_ONLY CON_PREFIX_SEP_;
        container_seg_by_ctx (vb_, &vb->contexts[VCF_TOPLUFT], (ContainerP)&top_luft, luft_only_prefix, strlen (luft_only_prefix), 0);
    }
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
        dict_id.num == dict_id_INFO_MLEAC           ||
        dict_id.num == dict_id_INFO_MLEAF           ||
        dict_id.num == dict_id_INFO_MQ0             ||
        dict_id.num == dict_id_INFO_LIFTOVER        ||
        dict_id.num == dict_id_INFO_LIFTBACK        ||
        dict_id.num == dict_id_INFO_REJTOVER        ||
        dict_id.num == dict_id_INFO_REJTBACK        ||

        // AC_* AN_* AF_* are small
        ((dict_id.id[0] == ('A' | 0xc0)) && (dict_id.id[1] == 'C' || dict_id.id[1] == 'F' || dict_id.id[1] == 'N') && dict_id.id[2] == '_');
}


// traverses the FORMAT field, gets ID of subfield, and moves to the next subfield
static DictId vcf_seg_get_format_subfield (const char **str, uint32_t *len) // remaining length of line 
{
    unsigned i=0; for (; i < *len && (*str)[i] != ':' && (*str)[i] != '\t' && (*str)[i] != '\n'; i++);
    
    DictId dict_id;
    // case: normal field - starts with a letter or another character in the range
    if ((*str)[0] >= 64 && (*str)[0] <= 127) 
        dict_id = dict_id_make (*str, i, DTYPE_VCF_FORMAT);
    
    // case: unusual field - starts with an out-range character, eg a digit - prefix with @ so its a legal FORMAT dict_id
    else {
        SAFE_ASSIGN (*str - 1, '@');
        dict_id = dict_id_make (*str-1, i+1, DTYPE_VCF_FORMAT);
        SAFE_RESTORE;
    }

    *str += i+1;
    *len -= i+1;
    return dict_id; 
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

    dl->has_haplotype_data = (str[0] == 'G' && str[1] == 'T' && (str[2] == ':' || field_len==2)); // GT field in FORMAT columns - must always appear first per VCF spec (if at appears)
    dl->has_genotype_data  = (field_len > 2 || (!dl->has_haplotype_data && field_len > 0));

    bool last_item = false;
    do {
        ASSVCF (con_nitems (format_mapper) < MAX_SUBFIELDS,
                "FORMAT field has too many subfields, the maximum allowed is %u: \"%.*s\"",  
                MAX_SUBFIELDS, field_len, field_start);

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


static const char *vcf_seg_samples (VBlockVCF *vb, ZipDataLineVCF *dl, int32_t *len, const char *next_field, 
                                    bool *has_13)
{
    // Container for samples - we have:
    // - repeats as the number of samples in the line (<= vcf_num_samples)
    // - num_items as the number of FORMAT subfields (inc. GT)

    Container samples = *ENT (Container, vb->format_mapper_buf, dl->format_node_i); // make a copy of the template

    const char *field_start;
    unsigned field_len=0, num_colons=0;

    bool is_ref_alt_switch          = last_ostatus == LO_OK_REF_ALT_SWTCH;
    bool needs_luft_validation      = is_ref_alt_switch && z_dual_coords && vb->line_coords == DC_PRIMARY;
    bool needs_lift_back_to_primary = is_ref_alt_switch && vb->line_coords == DC_LUFT;

    // 0 or more samples
    for (char separator=0 ; separator != '\n'; samples.repeats++) {

        field_start = next_field;
        next_field = seg_get_next_item (vb, field_start, len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, has_13, "sample-subfield");

        ASSVCF (field_len, "Error: invalid VCF file - expecting sample data for sample # %u, but found a tab character", 
                samples.repeats+1);

        vcf_seg_one_sample (vb, dl, &samples, samples.repeats, (char *)field_start, field_len, true, 
                            needs_luft_validation, needs_lift_back_to_primary,
                            &num_colons, has_13);

        ASSVCF (samples.repeats < vcf_num_samples || separator == '\n',
                "invalid VCF file - expecting a newline after the last sample (sample #%u)", vcf_num_samples);
    }

    ASSVCF (samples.repeats <= vcf_num_samples, "according the VCF header, there should be %u sample%s per line, but this line has %u samples - that's too many",
            vcf_num_samples, vcf_num_samples==1 ? "" : "s", samples.repeats);

    // in some real-world files I encountered have too-short lines due to human errors. we pad them
    if (samples.repeats < vcf_num_samples) {
        WARN_ONCE ("FYI: the number of samples in vb->line_i=%u is %u, different than the VCF column header line which has %u samples",
                   vb->line_i, samples.repeats, vcf_num_samples);

        if (dl->has_haplotype_data) {
            char *ht_data = ENT (char, vb->ht_matrix_ctx->local, vb->line_i * vb->ploidy * vcf_num_samples + vb->ploidy * samples.repeats);
            unsigned num_missing = vb->ploidy * (vcf_num_samples - samples.repeats); 
            memset (ht_data, '*', num_missing);
        }
    }
    
    container_seg_by_ctx (vb, &vb->contexts[VCF_SAMPLES], &samples, 0, 0, samples.repeats + num_colons); // account for : and \t \r \n separators

    if (vb->ht_matrix_ctx)
        vb->ht_matrix_ctx->local.len = (vb->line_i+1) * vb->ht_per_line;
 
    return next_field;
}

// complete haplotypes of lines that don't have GT, if any line in the vblock does have GT.
// In this case, the haplotype matrix must include the lines without GT too
static void vcf_seg_complete_missing_lines (VBlockVCF *vb)
{
    for (vb->line_i=0; vb->line_i < (uint32_t)vb->lines.len; vb->line_i++) {

        if (vb->ht_matrix_ctx && !DATA_LINE (vb->line_i)->has_haplotype_data) {
            char *ht_data = ENT (char, vb->ht_matrix_ctx->local, vb->line_i * vb->ht_per_line);
            memset (ht_data, '*', vb->ht_per_line);

            // NOTE: we DONT set dl->has_haplotype_data to true bc downstream we still
            // count this row as having no GT field when analyzing gt data
        }
    }

    vb->ht_matrix_ctx->local.len = vb->lines.len * vb->ht_per_line;
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

    // allow un-segging of o* and LIFTXXXX in case due to INFO or FORMAT ostatus changes from OK to REJECT
    if (z_dual_coords) vcf_lo_set_rollback_point (vb); 

    // copy line to lo_rejects now, because we are going to destroy it. We remove it later if its not a reject.
    uint64_t save_lo_rejects_len = vb->lo_rejects[vb->line_coords-1].len;
    uint32_t line_len;
    if (z_dual_coords)
        line_len = vcf_seg_copy_line_to_reject (vb, field_start_line, remaining_txt_len);
        
    GET_NEXT_ITEM (VCF_CHROM);
    if (vb->line_coords == DC_PRIMARY) 
        dl->chrom_index[0] = seg_chrom_field (vb_, VCF_CHROM_str, VCF_CHROM_len);
    
    else { // LUFT
        dl->chrom_index[1] = seg_by_did_i (vb_, VCF_CHROM_str, VCF_CHROM_len, VCF_oCHROM, VCF_CHROM_len); // we will add_bytes of the CHROM field (not oCHROM) as genounzip reconstructs the PRIMARY
        random_access_update_chrom (vb_, DC_LUFT, dl->chrom_index[1], VCF_CHROM_str, VCF_CHROM_len);
        vb->contexts[VCF_CHROM].txt_len++; // account for the tab
    }

    GET_NEXT_ITEM (VCF_POS);

    if (vb->line_coords == DC_PRIMARY) {
        PosType pos = dl->pos[0] = seg_pos_field (vb_, VCF_POS, VCF_POS, false, false, '.', VCF_POS_str, VCF_POS_len, 0, VCF_POS_len+1);

        // POS <= 0 not expected in a VCF file
        if (pos <= 0 && !(*VCF_POS_str == '.' && VCF_POS_len == 1))
            WARN_ONCE ("FYI: invalid POS=%"PRId64" value in vb_i=%u vb_line_i=%u: line will be compressed, but not indexed", 
                       pos, vb->vblock_i, vb->line_i);
                
        random_access_update_pos (vb_, DC_PRIMARY, VCF_POS);
    }

    else { // LUFT
        dl->pos[1] = seg_pos_field (vb_, VCF_oPOS, VCF_oPOS, false, false, '.', VCF_POS_str, VCF_POS_len, 0, VCF_POS_len);
        random_access_update_pos (vb_, DC_LUFT, VCF_oPOS);
        vb->contexts[VCF_POS].txt_len++; // account for the tab
    }

    GET_NEXT_ITEM (VCF_ID);
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
        vcf_lo_seg_generate_INFO_LIFTXXXX (vb, dl);

    // INFO
    if (vcf_num_samples) {
        GET_NEXT_ITEM (VCF_INFO); // to do: currently, we require FORMAT (and hence a \t after INFO) in the case the file has samples but this line doesn't. VCF spec doesn't require FORMAT in this case.
    } else {
        GET_MAYBE_LAST_ITEM (VCF_INFO); // may or may not have a FORMAT field
    }

    // set oSTATUS except --chain (already set)
    if (!chain_is_loaded) 
        vcf_set_ostatus (vb->reject_bytes ? LO_REJECTED : // reject lines in Luft are all rejected (happens only in txt_file->coords==DC_LUFT)
                         txt_file->coords ? LO_UNKNOWN  : // we don't know yet, we will test for existance of INFO/REJTXXXX in vcf_seg_info_field_correct_for_dual_coordinates
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

        ASSVCF0 (separator == '\n' || dl->has_genotype_data || dl->has_haplotype_data,
                 "expecting line to end as it has no genotype or haplotype data, but it is not");

        if (separator != '\n') // has samples
            next_field = vcf_seg_samples (vb, dl, &len, next_field, has_13); // All sample columns
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

    // Adds INFO/LIFTXXXX according to ostatus, finalizes INFO/SF and segs the INFO container
    vcf_finalize_seg_info (vb);

    if (!has_samples && vcf_num_samples)
        WARN_ONCE ("FYI: vb->line_i=%u has no samples", vb->line_i);

    // seg coords - native or dual
    bool lo_ok = true; // default if not dual coords
    if (z_dual_coords) {
        const char *name = coords_names[LO_IS_OK (last_ostatus) ? DC_BOTH : vb->line_coords]; 
        seg_by_did_i (vb, name, strlen (name), VCF_COORDS, 0); // 0 as its not in the txt data

        lo_ok = LO_IS_OK (last_ostatus);
    } 
    
    SEG_EOL (VCF_EOL, false);

    // if line was NOT rejected (the default, if not dual coordinates), we can delete the text from lo_rejects
    if (lo_ok)
        vb->lo_rejects[vb->line_coords-1].len = save_lo_rejects_len;

    // if this line is a LUFT-only line, it won't be reconstructed by default
    if (!lo_ok && !vb->reject_bytes && vb->line_coords == DC_LUFT)
        vb->vb_data_size -= line_len; 
        
    // in case of a reject line - update reject_bytes to indicate its consumption
    if (txt_file->coords && txt_file->coords != vb->line_coords)
        vb->reject_bytes -= next_field - field_start_line;

// TODO: rollback txt_len of all contexts...currently Genozip only supports oREFALT same size as REFALT 
// TODO: in rejects ##luft-only, we account for LUFT as if it is primary?

    // test if still sorted
    if (flag.sort) {
        vcf_seg_assign_dl_sorted (vb, dl, VCF_CHROM);
        if (z_dual_coords) vcf_seg_assign_dl_sorted (vb, dl, VCF_oCHROM);
    }

    return next_field;
}
