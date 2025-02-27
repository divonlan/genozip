// ------------------------------------------------------------------
//   fastq_desc.c
//   Copyright (C) 2020-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"
#include "qname.h"
#include "tokenizer.h"

void fastq_seg_QNAME (VBlockFASTQP vb, STRp(qname), uint32_t line1_len, bool deep, uint32_t uncanonical_suffix_len)
{
    if (deep)
        fastq_deep_seg_QNAME (vb, FASTQ_QNAME, STRa(qname), uncanonical_suffix_len, qname_len + 1); // +1 for '@'

    else
        qname_seg (VB, QNAME1, STRa(qname), 1); // account for the '@' (segged as a toplevel container prefix)

    set_last_txt (FASTQ_QNAME, qname);
}

bool fastq_is_line3_copy_of_line1 (STRp(qname), STRp(line3), uint32_t desc_len)
{
    return !segconf.desc_is_l3 && 
           (line3_len == qname_len + (segconf.has_desc ? 1 + desc_len : 0)) &&
           !memcmp (line3, qname, line3_len);
}

void fastq_seg_LINE3 (VBlockFASTQP vb, STRp(qline3), STRp(qline1), STRp(desc))
{
    switch (segconf.line3) {
        case L3_EMPTY:  // no segging - we will drop the line from top_level
            ASSSEG (!qline3_len || segconf.optimize[FASTQ_QNAME], "Invalid FASTQ file format (#1): expecting middle line to be a \"+\", but it is \"+%.*s\"", STRf(qline3));
            CTX(FASTQ_LINE3)->txt_len++;      // account for the '+' (it is segged in the toplevel container)
            break;

        case L3_COPY_LINE1:
            ASSSEG (fastq_is_line3_copy_of_line1 (STRa(qline1), STRa(qline3), desc_len),
                    "Invalid FASTQ file format (#2): expecting middle line to be a \"+\" followed by a copy of the description line, but it is \"%.*s\"", STRf(qline3)); 
            seg_special0 (VB, FASTQ_SPECIAL_copy_line1, CTX(FASTQ_LINE3), 1 + qline3_len); // +1 - account for the '+' (segged as a toplevel container prefix)
            break;

        case L3_NCBI: 
            // note: we only need to seg QLINE3, DESC (if exists) was already segged in fastq_seg_txt_line
            qname_seg (VB, QLINE3, STRa(qline3), 1);
            break;

        default:
            ABOSEG ("Invalid segconf.line3=%d. Perhaps segconf failed to run?", segconf.line3);
    }
}

//--------------------------------------------
// AUX: eg "length=1234" or "l:1234" 
//--------------------------------------------

static void fastq_seg_one_aux (VBlockFASTQP vb, STRp(tag_name), STRp(value))
{
    DictId dict_id = dict_id_make (STRa(tag_name), DTYPE_2);
    ContextP ctx = ctx_get_ctx_tag (vb, dict_id, tag_name, tag_name_len); // create if it doesn't already exist
    
    seg_integer_or_not (VB, ctx, STRa(value), value_len); // also sets last_value

    // seq_len_dict_id must be detected in the first line, and if it appears, it must appear in all lines
    if (segconf_running && dict_id.num == _FASTQ_AUX_LENGTH && !vb->line_i && str_is_int (STRa(value))) 
        segconf.seq_len_dict_id = ctx->dict_id; // possibly overriding the value set in qname_segconf_discover_flavor
}

static void fastq_seg_aux_container (VBlockFASTQP vb, STRps(tag), uint32_t total_tag_len) // tags including '='/':'
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

// called for segconf line_i=0. DESC = QNAME2 + EXTRA + AUX <or> DESC = SAUX
void fastq_segconf_analyze_DESC (VBlockFASTQP vb, STRp(desc))
{
    ASSERTNOTZERO (segconf_running);

    str_split (desc, desc_len, MAX_DESC_FIELDS - 100, ' ', item, false);
    if (desc_len && !n_items) {
        segconf.has_extra = true;
        return; // too many fields - entire desc will be treated as "extra" and added to DESC.local
    }

    int n_auxes=0;

    segconf.aux_sep = (str_count_char (STRi(item, n_items-1), ':') == 1) ? ':' : '=';

    // count auxes (an aux is eg "length=151")
    for (int i=n_items-1; i >= 0; i--, n_auxes++) {
        str_split (items[i], item_lens[i], 2, segconf.aux_sep, side, true); // an AUX field is a name=value pair, eg "length=151"
        if (n_sides != 2) break;

        // set segconf.has[]. 
        DictId dict_id = dict_id_make (STRi(side,0), DTYPE_2);
        ContextP ctx = ctx_get_ctx_tag (vb, dict_id, sides[0], side_lens[0]); // create if it doesn't already exist
        segconf.has[ctx->did_i] = true;
    }

    if (n_items - n_auxes > 0)
        qname_segconf_discover_flavor (VB, QNAME2, STRi(item,0));   // note: also discovers the original TECH, if file has NCBI qnames or when optimize_DESC

    segconf.has_qname2 |= (segconf.qname_flavor[QNAME2] != 0);      // first item, apart from the auxes is the QNAME2 - but only if it has a flavor
    segconf.has_extra  |= (n_items - n_auxes > segconf.has_qname2); // extra info that is not the QNAME2 and not AUX
    segconf.has_aux    |= (n_auxes > 0);    
}

void fastq_seg_DESC (VBlockFASTQP vb, STRp(desc), bool deep_qname2, uint32_t uncanonical_suffix_len)
{
    START_TIMER;

    str_split (desc, desc_len, 0, ' ', item, false);
    
    // edge case: too many fields - just make it 1
    if (n_items > MAX_DESC_FIELDS - 100) { 
        n_items = 1;
        items[0] = desc;
        item_lens[0] = desc_len;
    }

    // case: we are allowed to seg auxes
    if (segconf.has_aux) { // we are allowed to seg AUX
        int n_auxes=0, total_tag_len=0;
        int i; for (i=n_items-1; i >= 0; i--, n_auxes++) {
            str_split (items[i], item_lens[i], 2, segconf.aux_sep, side, true); // an AUX field is a name=value pair, eg "length=151"
            if (n_sides != 2) break;

            fastq_seg_one_aux (vb, STRi(side,0), STRi(side,1));

            item_lens[i] = side_lens[0]; // shorten to the tag_name only - will be used to build the container
            total_tag_len += item_lens[i] + 1; // +1 for '='
        }
        
        if (n_auxes) 
            fastq_seg_aux_container (vb, n_auxes, &items[i+1], &item_lens[i+1], total_tag_len);
        
        // case: we expected to have aux - but we don't. delete the ' ' added by the toplevel container prefix
        else 
            seg_special0 (VB, FASTQ_SPECIAL_backspace, CTX(FASTQ_AUX), 0);
    
        n_items -= n_auxes;
    }

    // edge case: segconf determined that desc has aux fields only, and we have now encountered one field has non-aux data.
    if (n_items && !segconf.has_qname2 && !vb->has_extra) {
        // set has_extra, and seg a backspace for all previous lines
        vb->has_extra = true;
        for (int32_t line_i=0; line_i < (int32_t)vb->line_i - 1; line_i++) // signed, so it works with line_i==0 too
            seg_special0 (VB, FASTQ_SPECIAL_backspace, CTX(FASTQ_EXTRA), 0);
    }
    
    if (vb->has_extra) {
        // case: we have fields beyond a possible QNAME2
        if (n_items > segconf.has_qname2) {
            int32_t extra_len = -1; // -1 bc last item doesn't have an internal separator
            for (int i=segconf.has_qname2; i < n_items; i++)
                extra_len += item_lens[i] + 1; // +1 for internal separator

            seg_by_did (VB, items[segconf.has_qname2], extra_len, FASTQ_EXTRA, extra_len + 1); // +1 - account for ' ' before extra
        }

        // case: we are expected to have extra - but we don't. delete the ' ' added by the toplevel container prefix
        else 
            seg_special0 (VB, FASTQ_SPECIAL_backspace, CTX(FASTQ_EXTRA), 0);

        n_items = MIN_((int)segconf.has_qname2, n_items); // 0 or 1
    }

    // seg QNAME2, if there is one and we are allowed to seg it
    if (segconf.has_qname2) {
        // case: we are expected to have qname2 - but we don't. delete the ' ' added by the toplevel container prefix
        if (n_items == 0)
            seg_special0 (VB, FASTQ_SPECIAL_backspace, CTX(FASTQ_QNAME2), 0);

        else if (n_items == 1 && deep_qname2) 
            fastq_deep_seg_QNAME (vb, FASTQ_QNAME2, STRi(item,0), uncanonical_suffix_len, item_lens[0] + 1); // +1 for the ' '

        else if (n_items == 1) 
            qname_seg (VB, QNAME2, STRi(item,0), 1); // account for the ' ' (segged as a toplevel container prefix)

        // case: n_items > 1, this happens when we have extra but vb->has_extra=false. we tokenize.
        else {
            int32_t qname2_len = -1; // -1 bc last item doesn't have an internal separator
            for (int i=0; i < n_items; i++)
                qname2_len += item_lens[i] + 1; // +1 for internal separator

            tokenizer_seg (VB, CTX(FASTQ_QNAME2), items[0], qname2_len, sep_with_space, 1/*+1 for leading ' '*/);
        }        
    }

    COPY_TIMER (fastq_seg_DESC);
}

Multiplexer2P fastq_get_ultima_c_mux (VBlockP vb)
{
    return &VB_FASTQ->mux_ultima_c;
}

//-----------
// PIZ
//-----------

SPECIAL_RECONSTRUCTOR (fastq_special_backspace)
{
    ASSERTNOTZERO (Ltxt);
    Ltxt--;

    return NO_NEW_VALUE;
}

SPECIAL_RECONSTRUCTOR (fastq_special_copy_line1)
{
    if (!reconstruct) return NO_NEW_VALUE;

    SETlast (snip, FASTQ_QNAME);
    RECONSTRUCT_snip;

    static Did dids[] = { FASTQ_QNAME2, FASTQ_EXTRA, FASTQ_AUX };
    for (int i=0; i < ARRAY_LEN (dids); i++)
        if (ctx_encountered_in_line (vb, dids[i])) {
            SETlast (snip, dids[i]);
            RECONSTRUCT1 (' ');
            RECONSTRUCT_snip;
        };

    return NO_NEW_VALUE;
}
