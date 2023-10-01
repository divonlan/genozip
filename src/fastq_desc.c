// ------------------------------------------------------------------
//   fastq_desc.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"
#include "seg.h"
#include "piz.h"
#include "deep.h"
#include "aligner.h"
#include "refhash.h"
#include "coverage.h"
#include "qname.h"
#include "tokenizer.h"

#define MAX_DESC_FIELDS (MAX_FIELDS-100)

void fastq_seg_QNAME (VBlockFASTQP vb, STRp(qname), uint32_t line1_len, bool deep, uint32_t uncanonical_suffix_len)
{
    // case: --optimize_DESC: we replace the description with "filename.line_i" (optimized string stored in vb->optimized_desc)
    unsigned optimized_len = 0; 
    if (flag.optimize_DESC) {
        optimized_len  = vb->optimized_desc_len + str_int (vb->first_line + vb->line_i, &vb->optimized_desc[vb->optimized_desc_len]);   
        vb->recon_size -= line1_len - optimized_len;

        qname = vb->optimized_desc;
        qname_len = optimized_len;
    }

    if (deep)
        fastq_deep_seg_QNAME (vb, FASTQ_QNAME, STRa(qname), uncanonical_suffix_len, qname_len + 1); // +1 for '@'

    else
        qname_seg (VB, QNAME1, STRa(qname), 1); // account for the '@' (segged as a toplevel container prefix)
}

bool fastq_is_line3_copy_of_line1 (STRp(qname), STRp(line3), uint32_t desc_len)
{
    return !segconf.desc_is_l3 && 
           (line3_len == qname_len + (segconf.has_desc ? 1 + desc_len : 0)) &&
           !memcmp (line3, qname, line3_len);
}

void fastq_seg_LINE3 (VBlockFASTQP vb, STRp(qline3), STRp(qline1), STRp(desc))
{
    if (flag.optimize_DESC) {
        vb->recon_size -= qline3_len;    // no segging - we will drop the line from top_level
        CTX(FASTQ_LINE3)->txt_len++;     // account for the '+' (it is segged in the toplevel container)
    }
    
    else if (segconf.line3 == L3_EMPTY) { // no segging - we will drop the line from top_level
        ASSSEG (!qline3_len, "Invalid FASTQ file format: expecting middle line to be a \"+\", but it is \"%.*s\"", STRf(qline3));
        CTX(FASTQ_LINE3)->txt_len++;      // account for the '+' (it is segged in the toplevel container)
    }

    else if (segconf.line3 == L3_COPY_LINE1) {
        ASSSEG (fastq_is_line3_copy_of_line1 (STRa(qline1), STRa(qline3), desc_len),
                "Invalid FASTQ file format: expecting middle line to be a \"+\" followed by a copy of the description line, but it is \"%.*s\"", STRf(qline3)); 
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_copy_line1 }, 2, FASTQ_LINE3, 1 + qline3_len); // +1 - account for the '+' (segged as a toplevel container prefix)
    }

    else if (segconf.line3 == L3_NCBI) 
        // note: we only need to seg QLINE3, DESC (if exists) was already segged in fastq_seg_txt_line
        qname_seg (VB, QLINE3, STRa(qline3), 1);

    else 
        ABOSEG ("Invalid FASTQ file format: expecting middle line to be a \"+\" with or without a copy of the description, but it is \"%.*s\"",
                STRf(qline3));
}

//-----------------
// AUX: eg len=1234
//-----------------

static void fastq_seg_one_aux (VBlockFASTQP vb, STRp(tag_name), STRp(value))
{
    DictId dict_id = dict_id_make (STRa(tag_name), DTYPE_2);
    ContextP ctx = ctx_get_ctx_tag (vb, dict_id, tag_name, tag_name_len); // create if it doesn't already exist
    
    switch (dict_id.num) {
        
        case _FASTQ_AUX_LENGTH: 
            seg_integer_or_not (VB, ctx, STRa(value), value_len); // also sets last_value

            // seq_len_dict_id must be detected in the first line, and if it appears, it must appear in all lines
            if (segconf.running && !vb->line_i && str_is_int (STRa(value))) 
                segconf.seq_len_dict_id = ctx->dict_id; // possibly overriding the value set in qname_segconf_discover_flavor

            break;

        default: 
            seg_by_ctx (VB, STRa(value), ctx, value_len); // note: the tag and '=' are accounted for by the AUX container
    }
}

static void fastq_seg_aux_container (VBlockFASTQP vb, STRps(tag), uint32_t total_tag_len) // tags including '='
{
    Container con = { .repeats = 1, .drop_final_item_sep = true };
    con_set_nitems (con, n_tags + kraken_is_loaded);

    char prefixes[n_tags + total_tag_len + 3 + 7];        // each name is tag= followed CON_PX_SEP ; +3 for the initial CON_PX_SEP. +7 for kraken
    prefixes[0] = prefixes[1] = CON_PX_SEP; // initial CON_PX_SEP follow by separator of empty Container-wide prefix followed by separator for empty prefix for translator-only item[0]
    unsigned prefixes_len = 2;
    
    // add AUX fields to the container
    for (int16_t tag_i = 0; tag_i < n_tags; tag_i++) {

        con.items[tag_i] = (ContainerItem){ .dict_id = dict_id_make (STRi(tag, tag_i), DTYPE_2), .separator  = " " };
    
        mempcpy (&prefixes[prefixes_len], tags[tag_i], tag_lens[tag_i]); 
        prefixes_len += tag_lens[tag_i]; // including '='
        prefixes[prefixes_len++] = '=';
        prefixes[prefixes_len++] = CON_PX_SEP;
    }

    // add taxid= tag if needed
    if (kraken_is_loaded) {
        con.items[n_tags] = (ContainerItem){ .dict_id = { _FASTQ_TAXID } };  

        memcpy(&prefixes[prefixes_len], ((char[]){ 't','a','x','i','d','=',CON_PX_SEP }), 7);
        prefixes_len += 7;

        vb->recon_size += 6; // "taxid="
    }

    container_seg (vb, CTX(FASTQ_AUX), &con, prefixes, prefixes_len, 
                   total_tag_len + (kraken_is_loaded ? 6 : 0) + n_tags/*leading and internal ' '*/); 
}

//-----------------------------------------------------------------------------------
// SAUX: SAM style AUX. eg BC:Z:CTTA, separated by an arbitrary number of ' ' or '\t'
//-----------------------------------------------------------------------------------

static inline bool fastq_get_saux_field_len (rom f, uint32_t *len) // false if format is invalid
{
    *len = strcspn (f, " \t");

    if (*len < 5 || (f[3] == 'B' && *len < 7))
        return false;

    if (f[2] != ':' || f[4] != ':' || 
        !IS_LETTER(f[0]) || !IS_ALPHANUMERIC(f[1]) ||
        (f[3] != 'i' && f[3] != 'Z' && f[3] != 'B' && f[3] != 'A' && f[3] != 'H' && f[3] != 'f'))
        return false;

    if (f[3] == 'B' && 
        (f[6] != ':' || (f[5] != 'c' && f[5] != 'C' && f[5] != 's' && f[5] != 'S' && f[5] != 'i' && f[5] != 'I' && f[5] != 'f')))
        return false;

    return true;
}

// check if DESC consists of SAM-style fields, separated by one or more spaces or tabs
static bool fastq_segconf_is_saux (STRp(desc)) 
{
    SAFE_NULT (desc);
    uint32_t n_fields = 0, n_spaces, field_len;
    rom c = desc;

    while ((n_spaces = strspn (c, " \t"))) { // skip to c field
        c += n_spaces;

        if (!fastq_get_saux_field_len (c, &field_len)) {            
            n_fields = 0; // not saux format
            break;
        }

        n_fields++;
        c += field_len;
    }

    SAFE_RESTORE;
    return n_fields > 0;
}

static DictId fastq_seg_one_saux (VBlockFASTQP vb, rom tag, char sam_type, char array_subtype, STRp (value), unsigned add_bytes)
{
    char dict_name[6] = { tag[0], tag[1], ':', sam_type, ':', array_subtype }; // last 2 are ignored if not array
    DictId dict_id = dict_id_make (dict_name, (sam_type=='B' ? 6 : 4), DTYPE_FASTQ_AUX); // match dict_id as declared in #pragma GENDICT

    ValueType numeric = {};
    if (sam_type == 'i') {
        ASSSEG(str_get_int (STRa (value), &numeric.i), "%s: Expecting integer value for auxiliary field %c%c but found \"%.*s\"",
                LN_NAME, tag[0], tag[1], STRf (value));
        value = 0;
    }
        
    sam_seg_aux_field_fallback (VB, NULL, dict_id, sam_type, array_subtype, STRa(value), numeric, add_bytes);

    return dict_id;
}

static void fastq_seg_saux (VBlockFASTQP vb, STRp(desc))
{
    SAFE_NULT (desc);
    uint32_t n_fields = 0, n_spaces, field_len;
    rom c = desc;

    ASSSEG (desc_len < 64 KB, "desc too long:\n%.*s", STRf(desc));

    STRl(prefixes, desc_len) = 2;
    prefixes[0] = CON_PX_SEP;
    prefixes[1] = CON_PX_SEP;

    Container con = { .repeats = 1 };

    while ((n_spaces = strspn (c, " \t"))) { // skip to c field
        memcpy (&prefixes[prefixes_len], c, n_spaces);
        prefixes_len += n_spaces;

        c += n_spaces;
        
        ASSSEG (fastq_get_saux_field_len (c, &field_len), "Invalid SAM tags format in field %u: \"%.*s\"", n_fields, STRf(desc));

        uint32_t tag_name_len = (c[3] == 'B') ? 7 : 5;

        memcpy (&prefixes[prefixes_len], c, tag_name_len);
        prefixes[prefixes_len + tag_name_len] = CON_PX_SEP;
        prefixes_len += tag_name_len + 1;

        con.items[n_fields++].dict_id = fastq_seg_one_saux (vb, c, c[3], (c[3] == 'B') ? c[5] : 0, c + tag_name_len, field_len - tag_name_len, field_len - tag_name_len);

        c += field_len;
    }

    con_set_nitems (con, n_fields);

    container_seg (VB, CTX(FASTQ_AUX), &con, prefixes, prefixes_len, 
                   prefixes_len - n_fields - 1  ); // account for tags and whitespace
    SAFE_RESTORE;
}

// DESC = QNAME2 + EXTRA + AUX <or> DESC = SAUX
void fastq_segconf_analyze_DESC (VBlockFASTQP vb, STRp(desc))
{
    ASSERTNOTZERO (segconf.running);

    if (fastq_segconf_is_saux (STRa(desc))) {
        segconf.has_saux = true; // aux data in SAM format
        return;
    }

    str_split (desc, desc_len, MAX_DESC_FIELDS - 100, ' ', item, false);
    if (desc_len && !n_items) {
        segconf.has_extra = true;
        return; // too many fields - entire desc will be treated as "extra" and added to DESC.local
    }

    // count auxes
    int n_auxes=0;
    for (int i=n_items-1; i >= 0; i--, n_auxes++) {
        str_split (items[i], item_lens[i], 2, '=', side, true); // an AUX field is a name=value pair, eg "length=151"
        if (n_sides != 2) break;
    }

    if (n_items - n_auxes > 0)
        qname_segconf_discover_flavor (VB, QNAME2, STRi(item,0));   // note: also discovers the original TECH, if file has NCBI qnames or when optimize_DESC

    segconf.has_qname2 |= (segconf.qname_flavor[QNAME2] != 0);      // first item, apart from the auxes is the QNAME2 - but only if it has a flavor
    segconf.has_extra  |= (n_items - n_auxes > segconf.has_qname2); // extra info that is not the QNAME2 and not AUX
    segconf.has_aux    |= (n_auxes > 0) || kraken_is_loaded;        // with --kraken, we add an aux field taxid=
}

void fastq_seg_DESC (VBlockFASTQP vb, STRp(desc), bool deep_qname2, uint32_t uncanonical_suffix_len)
{
    START_TIMER;

    if (segconf.has_saux) {
        fastq_seg_saux (vb, STRa(desc));
        return;
    }

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
            str_split (items[i], item_lens[i], 2, '=', side, true); // an AUX field is a name=value pair, eg "length=151"
            if (n_sides != 2) break;

            fastq_seg_one_aux (vb, STRi(side,0), STRi(side,1));

            item_lens[i] = side_lens[0]; // shorten to the tag_name only - will be used to build the container
            total_tag_len += item_lens[i] + 1; // +1 for '='
        }
        
        if (n_auxes || kraken_is_loaded) 
            fastq_seg_aux_container (vb, n_auxes, &items[i+1], &item_lens[i+1], total_tag_len);
        
        // case: we expected to have aux - but we don't. delete the ' ' added by the toplevel container prefix
        else 
            seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_backspace }, 2, FASTQ_AUX, 0);
    
        n_items -= n_auxes;
    }

    // edge case: segconf determined that desc has aux fields only, and we have now encountered one field has non-aux data.
    if (n_items && !segconf.has_qname2 && !vb->has_extra) {
        // set has_extra, and seg a backspace for all previous lines
        vb->has_extra = true;
        for (int32_t line_i=0; line_i < (int32_t)vb->line_i - 1; line_i++) // signed, so it works with line_i==0 too
            seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_backspace }, 2, FASTQ_EXTRA, 0);
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
            seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_backspace }, 2, FASTQ_EXTRA, 0);

        n_items = MIN_((int)segconf.has_qname2, n_items); // 0 or 1
    }

    // seg QNAME2, if there is one and we are allowed to seg it
    if (segconf.has_qname2) {
        // case: we are expected to have qname2 - but we don't. delete the ' ' added by the toplevel container prefix
        if (n_items == 0)
            seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_backspace }, 2, FASTQ_QNAME2, 0);

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
