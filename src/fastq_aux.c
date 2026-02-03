// ------------------------------------------------------------------
//   fastq_aux.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"

//--------------------------------------------
// AUX: eg "length=1234" 
//--------------------------------------------

static inline void fastq_seg_one_aux (VBlockFASTQP vb, STRp(tag_name), STRp(value))
{
    DictId dict_id = dict_id_make (STRa(tag_name), DTYPE_2);
    ContextP ctx = ctx_get_ctx_tag (vb, dict_id, tag_name, tag_name_len); // create if it doesn't already exist
    
    switch (dict_id.num) {
        #define COND0(condition, seg) if (condition) { seg; break; } else  
        #define COND(condition,  seg) if (condition) { seg; break; } else goto fallback

        case _FASTQ_AUX_length         : // "length=1234" or "l:1234" 
            // seq_len_dict_id must be detected in the first line, and if it appears, it must appear in all lines
            if (segconf_running && !vb->line_i && str_is_int (STRa(value))) 
                segconf.seq_len_dict_id = ctx->dict_id; // possibly overriding the value set in qname_segconf_discover_flavor
            goto fallback;
        
        case _FASTQ_AUX_parent_read_id : 
            COND(TECH(NANOPORE), seg_maybe_copy (VB, ctx, FASTQ_QNAME, STRa(value), STRa(copy_qname_snip)));
        
        case _FASTQ_AUX_start_time     : 
            COND(TECH(NANOPORE), seg_diff (VB, ctx, NULL, STRa(value), false, value_len));
        
        default : fallback             : 
            seg_integer_or_not (VB, ctx, STRa(value), value_len); // also sets last_value    
    }

    seg_set_last_txt (VB, ctx, STRa(value));
}

static inline void fastq_seg_aux_container (VBlockFASTQP vb, STRps(tag), uint32_t total_tag_len) // tags including '='/':'
{
    Container con = { .repeats = 1, .drop_final_item_sep = true };
    con_set_nitems (con, n_tags);

    char prefixes[n_tags + total_tag_len + 3];        // each name is tag= followed CON_PX_SEP ; +3 for the initial CON_PX_SEP.
    prefixes[0] = prefixes[1] = CON_PX_SEP; // initial CON_PX_SEP follow by separator of empty Container-wide prefix followed by separator for empty prefix for translator-only item[0]
    unsigned prefixes_len = 2;
    
    // add AUX fields to the container
    for (int16_t tag_i = 0; tag_i < n_tags; tag_i++) {

        con.items[tag_i] = (ContainerItem){ .dict_id = dict_id_make (STRi(tag, tag_i), DTYPE_2), .separator  = " " };
    
        mempcpy (&prefixes[prefixes_len], tags[tag_i], tag_lens[tag_i]); 
        prefixes_len += tag_lens[tag_i]; // including '='/':'
        prefixes[prefixes_len++] = segconf.aux_sep;
        prefixes[prefixes_len++] = CON_PX_SEP;
    }

    container_seg (vb, CTX(FASTQ_AUX), &con, prefixes, prefixes_len, 
                   total_tag_len + n_tags/*leading and internal ' '*/); 
}

int fastq_seg_aux (VBlockFASTQP vb, STRps(item))
{
    int n_auxes=0, total_tag_len=0;
    int i; for (i=n_items-1; i >= 0; i--, n_auxes++) {
        str_split (items[i], item_lens[i], 2, segconf.aux_sep, side, true); // an AUX field is a name=value pair, eg "length=151"
        if (n_sides != 2) break;

        fastq_seg_one_aux (vb, STRi(side,0), STRi(side,1));

        ((uint32_t *)item_lens)[i] = side_lens[0]; // shorten to the tag_name only - will be used to build the container
        total_tag_len += item_lens[i] + 1; // +1 for '='
    }
    
    if (n_auxes) 
        fastq_seg_aux_container (vb, n_auxes, &items[i+1], &item_lens[i+1], total_tag_len);
    
    // case: we expected to have aux - but we don't. delete the ' ' added by the toplevel container prefix
    else 
        seg_special0 (VB, FASTQ_SPECIAL_backspace, CTX(FASTQ_AUX), 0);

    return n_auxes;
}