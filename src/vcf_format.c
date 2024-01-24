// ------------------------------------------------------------------
//   vcf_format.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "reconstruct.h"
#include "file.h"
#include "codec.h"

// -----------
// Seg
// -----------

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

void vcf_seg_FORMAT (VBlockVCFP vb, ZipDataLineVCF *dl, STRp(fmt))
{
    ASSVCF0 (fmt_len >= 1, "missing or invalid FORMAT field");

    if (!vcf_num_samples) {
        seg_by_did_i_ex (VB, fmt, fmt_len, VCF_FORMAT, fmt_len + 1 /* \n */, NULL);
        return; // if we're not expecting any samples, no need to analyze the FORMAT field
    }

    dl->has_haplotype_data = (fmt[0] == 'G' && fmt[1] == 'T' && (fmt[2] == ':' || fmt_len==2)); // GT tag in FORMAT field - must always appear first per VCF spec (if it appears)
    dl->has_genotype_data  = (fmt_len > 2 || !dl->has_haplotype_data);

    bool possibly_rename = z_is_dvcf && !vb->is_rejects_vb;
    bool possible_conditional_renaming = possibly_rename && (vb->last_index (VCF_oXSTRAND) || LO_IS_OK_SWITCH (last_ostatus));

    // case: FORMAT is the same as previous line - just use the same node_index, but only if there is no chance of conditional renaming
    if (str_issame_(STRa(fmt), STRb(vb->last_format)) && !possible_conditional_renaming) {
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
        .drop_final_repsep   = true,
        .drop_final_item_sep = true,
        .callback            = true,
        .filter_items        = true,
        .filter_repeats      = true,
        .repsep              = "\t"
    };

    ASSVCF (n_sf_names < MAX_FIELDS,
            "FORMAT field has too many subfields, the maximum allowed is %u: \"%.*s\"", MAX_FIELDS, STRf(fmt));

    con_set_nitems (format_mapper, n_sf_names);

    ContextPBlock ctxs;

    for (unsigned i=0; i < n_sf_names; i++) {

        DictId dict_id = vcf_seg_get_format_sf_dict_id (sf_names[i], sf_name_lens[i]);
        ctxs[i] = ctx_get_ctx_tag (vb, dict_id, sf_names[i], sf_name_lens[i]);

        if (possibly_rename) vcf_tags_add_tag (vb, ctxs[i], DTYPE_VCF_FORMAT, sf_names[i], sf_name_lens[i]);

        format_mapper.items[i] = (ContainerItem){ .dict_id = dict_id, .separator = {':'} };

        if (dict_id.num == _FORMAT_PS || dict_id.num == _FORMAT_PID || dict_id.num == _FORMAT_DP)
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
    rom renamed = renamed_data;
    uint32_t renamed_len = possibly_rename ? vcf_tags_rename (vb, n_sf_names, ctxs, sf_names, sf_name_lens, NULL, renamed_data) : 0;

    // note: fmt_len and renamed_len need to be int64_t to avoid -Wstringop-overflow warning in gcc 10
    SNIP(6 + fmt_len + renamed_len); // maximum length in case of dual-snip with renaming

    if (renamed_len) {

        if (vb->line_coords == DC_LUFT) {
            vb->recon_size += (int)renamed_len - (int)fmt_len;
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
        ASSVCF (node_index == vb->format_mapper_buf.len32, 
                "node_index=%d different than vb->format_mapper_buf.len=%u", node_index, vb->format_mapper_buf.len32);

        buf_alloc (vb, &vb->format_mapper_buf, 1, 0, Container, CTX_GROWTH, "format_mapper_buf");
        buf_alloc (vb, &vb->format_contexts, 1, 0, ContextPBlock, CTX_GROWTH, "format_contexts");
        vb->format_contexts.len =vb->format_mapper_buf.len = vb->format_mapper_buf.len + 1;
    }    

    ContainerP con = B(Container, vb->format_mapper_buf, node_index);
    if (is_new || !con_nitems (*con)) { // assign if not already assigned (con->nitem=0 if format_mapper was already segged in another VB (so it is in ol_nodes), but not yet this VB) 
        *con = format_mapper; 
        memcpy (B(ContextPBlock, vb->format_contexts, node_index), ctxs, sizeof (ctxs));
    }
}

// -----------
// PIZ
// -----------

// FORMAT - obey GT-only ; load haplotype line (note: if drop_genotypes we don't reach here, as FORMAT is eliminated by vcf_piz_is_skip_section)
SPECIAL_RECONSTRUCTOR (vcf_piz_special_FORMAT)
{
    bool has_GT = (snip_len>=2 && snip[0]=='G' && snip[1] == 'T' && (snip_len==2 || snip[2] == ':'));
    
    if (reconstruct) {
        if (flag.gt_only) {
            if (has_GT)
                RECONSTRUCT ("GT", 2);
        }
        else 
            RECONSTRUCT_snip;
    }

    // initialize haplotype stuff
    if (!VER(14) && has_GT && !vb->ht_matrix_ctx) 
        codec_hapmat_piz_calculate_columns (vb);

    return NO_NEW_VALUE;
}
