// ------------------------------------------------------------------
//   sam_parse.c
//   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "context.h"

typedef enum { PCQ_pN='0', PCQ_pB='1', PCQ_CR='2', PCQ_CB='3', PCQ_Q0NAME='4' } ParseCopyQnameType; // part of the file format!

void sam_parse_seg_initialize (VBlockSAMP vb)
{
    seg_mux_init (vb, OPTION_RE_Z, SAM_SPECIAL_DEMUX_by_empty_field, true, RE);
}

// Parse Biosciences barcode indices - also in QNAME
void sam_seg_pB_Z (VBlockSAMP vb, ZipDataLineSAM𐤐 dl, STRp(pb), unsigned add_bytes)
{
    STRlast (qname, SAM_QNAME);

    if (segconf_qf_id (QNAME1) == QF_PARSE_ILLUM && 
        str_count_char (STRa(pb), '_') && // exactly to underscores (e.g. 97_2_81)
        qname_len >= 13 + pb_len && !memcmp (pb, qname + 13, pb_len) && // included in qname starting character 13
        qname[13 + pb_len] == '_') // following qname character is '_'

        seg_special1 (VB, SAM_SPECIAL_Parse_copy_QNAME, PCQ_pB, CTX(OPTION_pB_Z), add_bytes);
    else 
        seg_by_did (VB, STRa(pb), OPTION_pB_Z, add_bytes);
}

// Parse Biosciences polyN (=UMI)
void sam_seg_pN_Z (VBlockSAMP vb, ZipDataLineSAM𐤐 dl, STRp(pn), unsigned add_bytes)
{
    STRlast (qname, SAM_QNAME); // note: using QNAME and not Q4NAME, as its available in PRIM/DEPN too
    STR(pn_in_qname);

    if (segconf_qf_id (QNAME1) == QF_PARSE_ILLUM && 
        str_item_i (STRa(qname), '_', 14, pSTRa(pn_in_qname)) &&
        str_issame (pn, pn_in_qname))

        seg_special1 (VB, SAM_SPECIAL_Parse_copy_QNAME, PCQ_pN, CTX(OPTION_pN_Z), add_bytes);
        
    else 
        seg_by_did (VB, STRa(pn), OPTION_pN_Z, add_bytes);
}

void sam_seg_CB_Z_Parse (VBlockSAMP vb, ZipDataLineSAM𐤐 dl, STRp(cb), unsigned add_bytes)
{
    decl_ctx(OPTION_CB_Z);
    STRlast(qname, SAM_QNAME);

    // Parse - we expect CB to be the first 8 characters of QNAME
    if (segconf_qf_id (QNAME1) == QF_PARSE_ILLUM &&
        cb_len == 8 && qname_len > 8 && !memcmp (cb, qname, 8))

        seg_special1 (VB, SAM_SPECIAL_Parse_copy_QNAME, PCQ_CB, ctx, add_bytes);
        
    else
        seg_by_ctx (VB, STRa(cb), ctx, add_bytes);
}

// e.g. CATCATCC_AACGTGAT_AAACATCG (length=26)
void sam_seg_CR_Z_Parse (VBlockSAMP vb, ZipDataLineSAM𐤐 dl, STRp(cr), unsigned add_bytes)
{
    decl_ctx(OPTION_CR_Z);    
    STRlast(qname, SAM_QNAME);
    STR(cr_in_qname);
    
    if (segconf_qf_id (QNAME1) == QF_PARSE_ILLUM && 
        cr_len == 26 &&
        str_item_i (STRa(qname), '_', 10, pSTRa(cr_in_qname)) &&
        cr_in_qname + 26 <= qname + qname_len &&
        !memcmp (cr, cr_in_qname, 26)) 

        seg_special1 (VB, SAM_SPECIAL_Parse_copy_QNAME, PCQ_CR, CTX(OPTION_CR_Z), add_bytes);
        
    else 
        seg_by_ctx (VB, STRa(cr), ctx, add_bytes);
}

// typically, the same Q0NAME appears in several consecutive rows, and nowhere else in the file.
// we store the first occurance in local, and the subsequent ones copy from previous QNAME.
// note: copying from QNAME and not Q0NAME allows it to work in PRIM/DEPN too
void seg_qname_parse_QNAME0_cb (VBlockP vb, ContextP q0name_ctx, STRp(q0name))
{
    STRlast (prev_qname, SAM_QNAME);

    if (!segconf.is_collated) // not expected
        seg_by_ctx (vb, STRa(q0name), q0name_ctx, q0name_len);

    else if (prev_qname_len > q0name_len && 
        !memcmp (q0name, prev_qname, q0name_len) &&
        q0name[q0name_len-1] == '_' && str_count_char (STRa(q0name), '_') == 10) // q0name has 10 '_', and ends after the 10th.
        
        seg_special1 (VB, SAM_SPECIAL_Parse_copy_QNAME, PCQ_Q0NAME, q0name_ctx, q0name_len);
        
    else
        seg_add_to_local_string (vb, q0name_ctx, STRa(q0name), LOOKUP_SIMPLE, q0name_len);
}

// Copy substrings of QNAME into pN, pB, CR and CB
SPECIAL_RECONSTRUCTOR (sam_piz_special_Parse_copy_QNAME)
{
    STRlast(qname, SAM_QNAME);
    STR0(sub_qname);

    if (!reconstruct) goto done;

    switch (*snip) {
        case PCQ_CB:
            RECONSTRUCT (qname, 8);
            break;

        case PCQ_CR:
            ASSPIZ0 (str_item_i (STRa(qname), '_', 10, pSTRa(sub_qname)), "failed to find CR in QNAME");
            RECONSTRUCT (sub_qname, 26);
            break;

        case PCQ_pB:
            for (int i=13, underscores=0; qname[i] != '_' || ++underscores < 3 ; i++) // stop at the third '_' starting at [13]
                sub_qname_len++;
            RECONSTRUCT (&qname[13], sub_qname_len);
            break;

        case PCQ_pN:
            ASSPIZ0 (str_item_i (STRa(qname), '_', 14, pSTRa(sub_qname)), "failed to find pN in QNAME");
            RECONSTRUCT_str (sub_qname);
            break;

        case PCQ_Q0NAME: // note: qname is previous lines' qname
            for (int i=0, underscores=0; qname[i] != '_' || ++underscores < 10 ; i++) // stop at the 10th '_' 
                sub_qname_len++;
            RECONSTRUCT (qname, sub_qname_len + 1);
            break;

        default:
            ASSPIZ (false, "Unexpected opcode '%c'(%u)", *snip, *snip);
    }

done:    
    return NO_NEW_VALUE;
}

// REgion (I=Intron, E=Exon, N=Inter-gene) (note: same field in cellranger is RE:A)
void sam_seg_RE_Z (VBlockSAMP vb, ZipDataLineSAM𐤐 dl, STRp(re), unsigned add_bytes)
{
    int channel_i = (ctx_encountered_in_line (VB, OPTION_GX_Z) && CTX(OPTION_GX_Z)->last_txt.len > 0);
    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, OPTION_RE_Z, (MultiplexerP)&vb->mux_RE, channel_i);

    seg_by_ctx (VB, STRa(re), channel_ctx, add_bytes);
    seg_by_did (VB, STRa(vb->mux_RE.snip), OPTION_RE_Z, 0);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_by_empty_field)
{
    int channel_i = (ctx_encountered_in_line (VB, OPTION_GX_Z) && CTX(OPTION_GX_Z)->last_txt.len > 0);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}
