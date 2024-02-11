// ------------------------------------------------------------------
//   vcf_clinvar.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "piz.h"
#include "reconstruct.h"

void vcf_seg_hgvs_consolidate_stats (VBlockVCFP vb, Did parent)
{
    ctx_consolidate_stats (VB, parent, INFO_HGVS_snp_pos, INFO_HGVS_snp_refalt, INFO_HGVS_del_start_pos, INFO_HGVS_ins_start_pos,
                           INFO_HGVS_del_end_pos,   INFO_HGVS_ins_end_pos, INFO_HGVS_delins_end_pos, INFO_HGVS_dup_end_pos,
                           INFO_HGVS_del_payload,   INFO_HGVS_ins_payload, INFO_HGVS_delins_payload, INFO_HGVS_no_payload, DID_EOL);
}

// ------------
// INFO/CLNHGVS
// ------------

// SNP case: "NC_000023.10:g.154507173T>G"
static bool vcf_seg_INFO_HGVS_snp (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    PosType64 pos = DATA_LINE (vb->line_i)->pos; // data in variant
    char pos_str[30];
    unsigned pos_str_len = str_int (pos, pos_str);

    if (value_len < 3 + pos_str_len) return false;

    rom v = &value[value_len - 3];
    if (v[0] != vb->main_ref[0] || v[1] != '>' || v[2] != vb->main_alt[0]) return false; // REF/ALT differs

    v -= pos_str_len;
    if (memcmp (v, pos_str, pos_str_len)) return false; // POS differs

    SmallContainer con = { 
        .repeats   = 1,
        .nitems_lo = 2,
        .items = { { .dict_id.num = _INFO_HGVS_snp_pos    },
                   { .dict_id.num = _INFO_HGVS_snp_refalt } }
     }; 

    // temporarily surround prefix by separators, and seg container with prefix
    SAFE_ASSIGNx (&value[-1], CON_PX_SEP, 1);
    SAFE_ASSIGNx (v,          CON_PX_SEP, 2);

    container_seg (vb, ctx, (ContainerP)&con, &value[-1], v - value + 2, value_len - pos_str_len - 3);

    SAFE_RESTOREx(1);
    SAFE_RESTOREx(2);

    // seg special snips
    seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_HGVS_SNP_POS    }), 2, CTX(INFO_HGVS_snp_pos),    pos_str_len);
    seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_HGVS_SNP_REFALT }), 2, CTX(INFO_HGVS_snp_refalt), 3);

    return true;
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_SNP_POS)
{
    ContextP pos_ctx = CTX (VCF_POS);
    
    RECONSTRUCT_LAST_TXT (pos_ctx); // faster than RECONSTRUCT_INT    
    return NO_NEW_VALUE;
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_SNP_REFALT)
{
    rom refalt;
    reconstruct_peek (vb, CTX (VCF_REFALT), &refalt, NULL); // this special works only on SNPs, so length is always 3

    RECONSTRUCT1 (refalt[0]); // this might overwrite the "peeked" data, but that's ok
    RECONSTRUCT1 ('>');
    RECONSTRUCT1 (refalt[2]);

    return NO_NEW_VALUE;
}

// Deletion case:  "n.10571909_10571915delCCCGCCG" (POS=10571908 ; 10571909 are the bases except for the left-anchor, CCCGCCG is the payload of the indel)
//                 "NC_000001.10:g.5993401_5993405del" REF=TAAAAC ALT=T POS=5993397
//                 "NC_000001.10:g.6038478del" REF=AG ALT=A POS=6038476
// Duplication:    "NC_000001.10:g.977452dup" REF=A ALT=AG POS=977450
//                 "NC_000001.10:g.981867_981868dup" REF=G ALT=GCC POS=981860
// Insertion case: "n.10659939_10659940insTG" (POS=10659939 ; TG is the payload of the indel) 
// Delins:         "NC_000001.10:g.5987727_5987729delinsCCACG" POS=5987727 REF=GTT ALT=CCACG
typedef enum { DEL, INS, DELINS, DUP, INV, NUM_HGVS_TYPES } HgvsType;
static bool vcf_seg_INFO_HGVS_indel (VBlockVCFP vb, ContextP ctx, STRp(value), rom op, HgvsType t)
{
    rom payload = &op[t==DELINS ? 6 : 3];
    unsigned payload_len = (unsigned)(&value[value_len] - payload);

    if (payload_len) switch (t) {
        case DEL    : // Payload is expected to be the same as the REF field, without the first, left-anchor, base
                      if (!str_issame_(STRa(payload), vb->main_ref+1, vb->main_ref_len-1)) return false;
                      break;
        case INS    : // Payload is expected to be the same as the ALT field, without the first, left-anchor, base
                      if (!str_issame_(STRa(payload), vb->main_alt+1, vb->main_alt_len-1)) return false;
                      break;
                      // Payload is expected to be the same as the entire ALT field
        case DELINS : if (!str_issame (payload, vb->main_alt)) return false;
                      break;
        case DUP    : 
        case INV    : return false; // dup and inv are not expected to have a payload
        case NUM_HGVS_TYPES : {} // avoid compiler warning
    }

    // beginning of number is one after the '.' - scan backwards
    rom start_pos = op-1;
    while (start_pos[-1] != '.' && start_pos > value) start_pos--;
    if (start_pos[-1] != '.') return false;

    str_split (start_pos, op-start_pos, 2, '_', pos, false);

    PosType64 pos[2];
    if (!str_get_int (poss[0], pos_lens[0], &pos[0])) return false;
    
    if (n_poss == 2) {
        if (!str_get_int (poss[1], pos_lens[1], &pos[1])) return false;
        if (pos[0] + payload_len - 1 != pos[1] && t != INV) return false;
    }
    
    static const DictId dict_id_start_pos[NUM_HGVS_TYPES] = { { _INFO_HGVS_del_start_pos }, { _INFO_HGVS_ins_start_pos }, { _INFO_HGVS_ins_start_pos  }, { _INFO_HGVS_del_start_pos }, { _INFO_HGVS_ins_start_pos  } }; // INS, INV and DELINS use the same start_pos
    static const DictId dict_id_end_pos[NUM_HGVS_TYPES]   = { { _INFO_HGVS_del_end_pos   }, { _INFO_HGVS_ins_end_pos   }, { _INFO_HGVS_delins_end_pos }, { _INFO_HGVS_dup_end_pos   }, { _INFO_HGVS_delins_end_pos } };
    static const DictId dict_id_payload[NUM_HGVS_TYPES]   = { { _INFO_HGVS_del_payload   }, { _INFO_HGVS_ins_payload   }, { _INFO_HGVS_delins_payload }, { _INFO_HGVS_no_payload    }, { _INFO_HGVS_no_payload     } };                              

    static const Did did_i_start_pos[NUM_HGVS_TYPES] = { INFO_HGVS_del_start_pos, INFO_HGVS_ins_start_pos, INFO_HGVS_ins_start_pos,  INFO_HGVS_del_start_pos, INFO_HGVS_ins_start_pos  };
    static const Did did_i_end_pos[NUM_HGVS_TYPES]   = { INFO_HGVS_del_end_pos,   INFO_HGVS_ins_end_pos,   INFO_HGVS_delins_end_pos, INFO_HGVS_dup_end_pos,   INFO_HGVS_delins_end_pos };
    static const Did did_i_payload[NUM_HGVS_TYPES]   = { INFO_HGVS_del_payload,   INFO_HGVS_ins_payload,   INFO_HGVS_delins_payload, INFO_HGVS_no_payload,    INFO_HGVS_no_payload     };

    static const uint8_t special_end_pos[NUM_HGVS_TYPES]  = { VCF_SPECIAL_HGVS_DEL_END_POS, VCF_SPECIAL_HGVS_INS_END_POS, VCF_SPECIAL_HGVS_DELINS_END_POS, VCF_SPECIAL_HGVS_INS_END_POS, VCF_SPECIAL_HGVS_DELINS_END_POS };
    static const uint8_t special_payload[NUM_HGVS_TYPES]  = { VCF_SPECIAL_HGVS_DEL_PAYLOAD, VCF_SPECIAL_HGVS_INS_PAYLOAD, VCF_SPECIAL_HGVS_DELINS_PAYLOAD, 0                           , 0 };

    SmallContainer con = { 
        .repeats   = 1,
        .nitems_lo = 3,
        .items = { { .dict_id = dict_id_start_pos[t], .separator = "_" }, // separator deleted in container_reconstruct() if end_pos is missing
                   { .dict_id = dict_id_end_pos[t]                     },
                   { .dict_id = dict_id_payload[t]                     } } }; 

    // preper prefixes - a container-wide prefix, and the op is the prefix for item [2]
    // note: header_len, op_len and renamed_len prefixes_len to be int64_t to avoid -Wstringop-overflow warning in gcc 10
    int64_t header_len = start_pos - value;  
    int64_t op_len = (t==DELINS ? 6 : 3);
    int64_t prefixes_len = header_len + op_len + 5 /* separators */;

    char prefixes[prefixes_len];
    prefixes[0] = prefixes[header_len+1] = prefixes[header_len+2] = prefixes[header_len+3] = prefixes[header_len+4+op_len] = CON_PX_SEP;
    memcpy (&prefixes[1], value, header_len);
    memcpy (&prefixes[header_len+4], op, op_len);
    
    container_seg (vb, ctx, (ContainerP)&con, prefixes, prefixes_len, value_len - pos_lens[0] - (n_poss==2 ? pos_lens[1] : 0) - payload_len);

    CTX(did_i_start_pos[t])->flags.store = STORE_INT; // consumed by vcf_piz_special_INFO_HGVS_DEL_END_POS

    seg_pos_field (VB, did_i_start_pos[t], VCF_POS, 0, 0, 0, 0, pos[0], pos_lens[0]);
    
    // We pos_lens[1] only if the payload is longer than 1
    if (n_poss == 2)
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, special_end_pos[t] }), 2, CTX(did_i_end_pos[t]), pos_lens[1]);
    else
        seg_by_ctx (VB, NULL, 0, CTX(did_i_end_pos[t]), 0); // becomes WORD_INDEX_MISSING - container_reconstruct will remove the preceding _

    // the del payload is optional - we may or may not have it ; dup never has payload
    if (payload_len)
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, special_payload[t] }), 2, CTX(did_i_payload[t]), payload_len);
    else
        seg_by_ctx (VB, NULL, 0, CTX(did_i_payload[t]), 0); 

    return true;
}

// ClinVar: <ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
// COSMIC:  <ID=HGVSG,Number=1,Type=String,Description="HGVS genomic syntax">
bool vcf_seg_INFO_HGVS (VBlockP vb_, ContextP ctx, STRp(value), uint32_t repeat)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    if (ctx_encountered_in_line (VB, INFO_END)) 
        goto fail; // we can't use this if there is an END before CLNHGVS, as it will change last_int(VCF_POS) during reconstruction

    SAFE_NULT (value);

    bool success = false;
    rom op;
    if      (vb->main_ref_len == 1 && vb->main_alt_len == 1) success = vcf_seg_INFO_HGVS_snp   (vb, ctx, STRa(value));
    else if ((op = strstr (value, "delins")))                success = vcf_seg_INFO_HGVS_indel (vb, ctx, STRa(value), op, DELINS); // must be before del and ins
    else if ((op = strstr (value, "del")))                   success = vcf_seg_INFO_HGVS_indel (vb, ctx, STRa(value), op, DEL);
    else if ((op = strstr (value, "ins")))                   success = vcf_seg_INFO_HGVS_indel (vb, ctx, STRa(value), op, INS);
    else if ((op = strstr (value, "dup")))                   success = vcf_seg_INFO_HGVS_indel (vb, ctx, STRa(value), op, DUP);
    else if ((op = strstr (value, "inv")))                   success = vcf_seg_INFO_HGVS_indel (vb, ctx, STRa(value), op, INV);
    
    SAFE_RESTORE;
    
    if (success) return true; // indeed segged

fail:
    seg_by_ctx (VB, STRa(value), ctx, value_len); 
    return true; // indeed segged
}

static void vcf_piz_special_INFO_HGVS_INDEL_END_POS (VBlockP vb, HgvsType t)
{
    STR (refalt);
    reconstruct_peek (vb, CTX (VCF_REFALT), pSTRa(refalt)); // this special works only on SNPs, so length is always 3

    SAFE_NULT (refalt);
    rom tab = strchr (refalt, '\t');
    ASSPIZ (tab, "invalid REF+ALT=\"%.*s\"", STRf(refalt));
    SAFE_RESTORE;

    // reconstruct the DEL payload - this is the REF except for the first (anchor) base. reconstruct with one (possibly overlapping) copy
    static uint64_t start_pos_dnum[NUM_HGVS_TYPES] = { _INFO_HGVS_del_start_pos, _INFO_HGVS_ins_start_pos, _INFO_HGVS_ins_start_pos, _INFO_HGVS_del_start_pos }; // ins and delins share start_pos
    rom alt   = tab + 1;
    rom after = &refalt[refalt_len];

    PosType64 start_pos = ECTX (start_pos_dnum[t])->last_value.i;
    PosType64 end_pos = (t == DEL) ? (start_pos + tab - refalt - 2)
                      : (t == INS) ? (start_pos + after - alt - 2)
                      : /* DELINS */ (start_pos + after - alt - 1);

    RECONSTRUCT_INT (end_pos);
}

// three separate SPECIAL snips, so that they each because all_the_same
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_DEL_END_POS)    { vcf_piz_special_INFO_HGVS_INDEL_END_POS (vb, DEL);    return NO_NEW_VALUE; }
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_INS_END_POS)    { vcf_piz_special_INFO_HGVS_INDEL_END_POS (vb, INS);    return NO_NEW_VALUE; }
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_DELINS_END_POS) { vcf_piz_special_INFO_HGVS_INDEL_END_POS (vb, DELINS); return NO_NEW_VALUE; }

static void vcf_piz_special_INFO_HGVS_INDEL_PAYLOAD (VBlockP vb, HgvsType t)
{
    STR (refalt);
    reconstruct_peek (vb, CTX (VCF_REFALT), pSTRa(refalt)); // this special works only on SNPs, so length is always 3

    SAFE_NULT (refalt);
    rom tab = strchr (refalt, '\t');
    SAFE_RESTORE;

    ASSPIZ (tab, "invalid REF+ALT=\"%.*s\"", STRf(refalt));

    rom alt   = tab + 1;
    rom after = &refalt[refalt_len];

    rom payload = (t == DEL) ? (refalt + 1) // REF except for the anchor base
                : (t == INS) ? (alt + 1)    // ALT except for the anchor base
                : /* DELINS */ alt;         // the entire ALT

    PosType64 payload_len = (t == DEL) ? (tab - refalt - 1)
                          : (t == INS) ? (after - alt - 1)
                          : /* DELINS */ (after - alt);

    memmove (BAFTtxt, payload, payload_len);
    Ltxt += payload_len;
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_DEL_PAYLOAD)    { vcf_piz_special_INFO_HGVS_INDEL_PAYLOAD (vb, DEL);    return false; }
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_INS_PAYLOAD)    { vcf_piz_special_INFO_HGVS_INDEL_PAYLOAD (vb, INS);    return false; }
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_HGVS_DELINS_PAYLOAD) { vcf_piz_special_INFO_HGVS_INDEL_PAYLOAD (vb, DELINS); return false; }
