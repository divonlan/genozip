// ------------------------------------------------------------------
//   sam_star.c
//   Copyright (C) 2024-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"

void sam_star_zip_initialize (void)
{
    ZCTX(OPTION_jM_B_c)->con_rep_special = SAM_SPECIAL_jM_length;
}

void sam_star_seg_initialize (VBlockSAMP vb)
{
    ctx_set_store (VB, STORE_INT, OPTION_jM_B_c, DID_EOL); // required for con_rep_special
}

// jI:B:i,Start1,End1,Start2,End2,... Start and End of introns for all junctions (1-based).
void sam_seg_STAR_jI (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(raw), bool is_bam)
{
    decl_ctx (OPTION_jI_B_i);
    
    uint32_t n_vals = is_bam ? raw_len : (str_count_char (STRa(raw), ',') + 1); 
    int32_t vals[n_vals];

    if (is_bam) 
        memcpy (vals, raw, raw_len * sizeof (uint32_t)); // works because little endian
#ifdef __BIG_ENDIAN__
#error BIG ENDIAN not implemented here yet
#endif

    else {
        str_split_ints (raw, raw_len, n_vals, ',', val64, true);
        if (!n_val64s) goto fallback;

        for (int i=0; i < n_vals; i++)
            vals[i] = val64s[i]; // int64-> uint32_t
    }

    ARRAY (BamCigarOp, cigar, vb->binary_cigar);

    // case: no intron: expecting vals to have a single -1
    if (!vb->introns) {
        if (n_vals != 1 || vals[0] != -1) 
            goto fallback;
    }

    // case we have at least one intron: expecting vals to be pair of {first,last} POS of each intron
    else { 
        PosType32 pos = dl->POS;
        int op_i = 0;

        for (int i=0; i < n_vals; i += 2) {
            int M_len = 0;
            while (op_i < cigar_len && cigar[op_i].op != BC_N) {
                if (cigar[op_i].op == BC_M || cigar[op_i].op == BC_D || cigar[op_i].op == BC_E || cigar[op_i].op == BC_X)
                    M_len += cigar[op_i].n;
                op_i++;
            }

            if (op_i >= cigar_len ||  
                pos + M_len != vals[i] || // POS of first base of intron
                pos + M_len + cigar[op_i].n - 1 != vals[i+1]) // POS of last base of intron
                goto fallback;

            pos += M_len + cigar[op_i++].n;
        }
    }

    unsigned add_bytes = is_bam ? (4/*count*/ + 1/*type*/ + n_vals * sizeof (uint32_t)) 
                                : (2/*type - eg "i,"*/ + 1/*\t or \n*/ + raw_len);

    seg_special0 (VB, SAM_SPECIAL_jI, ctx, add_bytes);

    return;

fallback:
    sam_seg_array_one_ctx (vb, dl, ctx->dict_id, 'i', STRa(raw), NULL, NULL, NULL);
}

SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_jI)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    PosType32 pos = CTX(SAM_POS)->last_value.i;
    uint32_t count_N = 0;
    char *bam_array_len_p = 0; // pointer to char and not int32 bc not word-aligned
    char *next = 0;
    bool is_bam = OUT_DT(BAM) || OUT_DT(CRAM);

    RECONSTRUCT1 ('i'); // array type

    if (is_bam) {
        // leave room for count_N
        bam_array_len_p = BAFTtxt;
        next = BAFTtxt + 4;
    }

    for_cigar (vb->binary_cigar) {
        case BC_N:
            if (is_bam) {
                PUT_UINT32 (next, pos);
                next += sizeof (uint32_t);
                PUT_UINT32 (next, pos + op->n - 1);
                next += sizeof (uint32_t);
            }

            else { // SAM
                RECONSTRUCT1 (',');
                RECONSTRUCT_INT (pos); // intron start
                RECONSTRUCT1 (',');
                RECONSTRUCT_INT (pos + op->n - 1); // intron end
            }

            count_N++;
            // fallthrough

        case BC_M: case BC_D: case BC_E: case BC_X:
            pos += op->n;

        default: {} // BC_I, BC_S, BC_H, BC_P
    }

    if (!count_N) { // no intron
        if (is_bam) {
            PUT_UINT32 (next, (uint32_t)-1);
            next += sizeof (uint32_t);
        }

        else 
            RECONSTRUCT (",-1", 3);
    }

    if (is_bam) {
        PUT_UINT32 (bam_array_len_p, count_N ? count_N * 2 : 1);

        vb->txt_data.len32 = BNUM (vb->txt_data, next);
    }

    return NO_NEW_VALUE; 
}

// ZIP / PIZ
SPECIAL_RECONSTRUCTOR (sam_piz_special_jM_length)
{
    new_value->i = MAX_(1, VB_SAM->introns);

    // note: array length item reconstructor never actually reconstructs
    return HAS_NEW_VALUE;
}
