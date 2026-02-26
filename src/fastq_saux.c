// ------------------------------------------------------------------
//   fastq_saux.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"
#include "stats.h"

//-----------------------------------------------------------------------------------
// SAUX: SAM style AUX. eg BC:Z:CTTA, separated by an arbitrary number of ' ' or '\t'
//-----------------------------------------------------------------------------------

static inline DictId fastq_get_SAM_AUX_dict_id (STRp(f)) // DICT_ID_NONE if format is invalid
{
    if (f_len < 5 || (f[3] == 'B' && f_len < 6))
        return DICT_ID_NONE;

    if (f[2] != ':' || f[4] != ':' || 
        !IS_LETTER(f[0]) || !IS_ALPHANUMERIC(f[1]) ||
        (f[3] != 'i' && f[3] != 'Z' && f[3] != 'B' && f[3] != 'A' && f[3] != 'H' && f[3] != 'f'))
        return DICT_ID_NONE;

    if (f[3] == 'B' && 
        ((f[6] != ',' && f[6] != '\t' && f[6] != '\n' && f[6] != '\r') || 
         (f[5] != 'c' && f[5] != 'C' && f[5] != 's' && f[5] != 'S' && f[5] != 'i' && f[5] != 'I' && f[5] != 'f')))
        return DICT_ID_NONE;


    char dict_name[6] = { f[0], f[1], ':', f[3], ':', f[5] }; // last 2 are ignored if not array
    DictId dict_id = dict_id_make (dict_name, (f[3]=='B' ? 6 : 4), DTYPE_FASTQ_AUX); // match dict_id as declared in #pragma GENDICT

    return dict_id;
}

// called for segconf line_i=0. check if DESC consists of SAM-style fields, separated by tabs (i.e. parseable by eg bwa mem -C)
bool fastq_segconf_analyze_saux (VBlockFASTQP vb, STRp(saux)) 
{
    ASSERTNOTZERO (segconf_running);
    
    str_split (saux, saux_len, 0, '\t', field, false);

    if (!n_fields || n_fields > MAX_DESC_FIELDS) return false; // not SAUX

    for (int i=0; i < n_fields; i++) {
        DictId dict_id = fastq_get_SAM_AUX_dict_id (STRi(field, i));
        
        if (!dict_id.num) 
            return false; // invalid format - not SAUX
        
        segconf.has[ctx_get_ctx (VB, dict_id)->did_i] = true;
    }

    // note: same logic as in sam_seg_finalize_segconf
    if (segconf.has[OPTION_ZA_Z] && segconf.has[OPTION_ZB_Z] && segconf.has[OPTION_RX_Z] && segconf.has[OPTION_QX_Z] && segconf.has[OPTION_BC_Z]) {
        segconf.has_agent_trimmer = true;
        stats_add_one_program (_S("AGeNT_Trimmer"));
    }

    segconf.has_saux = true;
    segconf.saux_tab_sep = (saux[-1] == '\t');

    return true;
}

static void fastq_seg_one_saux (VBlockFASTQP vb, DictId dict_id, unsigned value_offset, STRp(field))
{
    char sam_type = field[3];
    char array_subtype = (field[3] == 'B') ? field[5] : 0;
    
    rom value = field + value_offset;
    unsigned value_len = field_len - value_offset;
    unsigned add_bytes = value_len + (array_subtype ? (1 + (field[6] == ',')) : 0) + 1/*preceding \t*/;

    ValueType numeric = {};
    if (sam_type == 'i') {
        ASSSEG (str_get_int (STRa (value), &numeric.i), "%s: Expecting integer value for SAM-format AUX field: \"%.*s\"",
                LN_NAME, STRf (field));
        value = 0;
    }
    
    ContextP ctx = ctx_get_ctx (vb, dict_id);

    #define COND0(condition, seg) if (condition) { seg; break; } else  
    #define COND(condition,  seg) if (condition) { seg; break; } else goto fallback

    if (segconf_running)
        segconf.has[ctx_get_ctx (VB, dict_id)->did_i]++;

    switch (dict_id.num) {
        case _OPTION_RX_Z : COND (segconf.has_agent_trimmer, agilent_seg_RX (VB, ctx, STRa(value), add_bytes)); // AGeNT Trimmer eg RX:Z:GGC-CCA
        case _OPTION_QX_Z : COND (segconf.has_agent_trimmer, agilent_seg_QX (VB, ctx, STRa(value), add_bytes)); // AGeNT Trimmer eg QX:Z:FFF FFD

        // note: for now sam_seg_MM_Z is quite useless, as its benefit is predicting length by ML:B
        case _OPTION_MM_Z : sam_seg_MM_Z (VB, STRa(value), add_bytes); break;

        default: fallback:
            sam_seg_aux_field_fallback (VB, NULL, dict_id, sam_type, array_subtype, STRa(value), numeric, add_bytes);
    }

    ctx_set_encountered (VB, ctx);
    set_last_txtC (ctx, value, value_len);
}

void fastq_seg_saux (VBlockFASTQP vb, STRp(saux))
{
    START_TIMER;

    str_split (saux, saux_len, 0, '\t', field, false);

    ASSSEG (n_fields <= MAX_DESC_FIELDS, "SAM-format AUX has %u fields, beyond the maximum of %u: \"%.*s\"", 
            n_fields, MAX_DESC_FIELDS, STRf(saux));

    STRl(prefixes, MIN_(saux_len, 8 * n_fields + 8)) = 2;
    prefixes[0] = CON_PX_SEP;
    prefixes[1] = CON_PX_SEP;

    Container con = { .repeats = 1, .drop_final_item_sep = true };
    con_set_nitems (con, n_fields);

    for (int f=0; f < n_fields; f++) {
        con.items[f] = (ContainerItem){ .dict_id   = fastq_get_SAM_AUX_dict_id (STRi(field, f)),
                                        .separator = "\t" };

        ASSSEG (con.items[f].dict_id.num, "Malformated SAM-format field: \"%.*s\"", STRfi(field, f));

        uint32_t value_offset = (fields[f][3] == 'B') ? (6 + (fields[f][6] == ',')) : 5;

        memcpy (&prefixes[prefixes_len], fields[f], 5); // e.g. ML:B:
        prefixes[prefixes_len + 5] = CON_PX_SEP;
        prefixes_len += 6;

        fastq_seg_one_saux (vb, con.items[f].dict_id, value_offset, STRi(field,f));
    }

    container_seg (VB, CTX(FASTQ_AUX), &con, prefixes, prefixes_len, n_fields * 5); // account for tags eg 'ML:B:'

    COPY_TIMER (fastq_seg_saux);
}
