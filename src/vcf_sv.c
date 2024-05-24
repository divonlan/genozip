// ------------------------------------------------------------------
//   vcf_sv.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "reconstruct.h"
#include "zip_dyn_int.h"
#include "reference.h"

static char copy_mate_snips[NUM_TWs][30];
uint32_t copy_mate_snip_lens[NUM_TWs];
sSTRl(snip_copy_cipos,16);
STRl(homlen_snip, 16);
STRl(duphomlen_snip, 16);
STRl(svinslen_snip, 16);
STRl(dupsvinslen_snip, 16);
sSTRl(cipos_snip,18);

void vcf_sv_zip_initialize (Did *tw_dids, int num_tw_dids)
{
    for (int tw=0; tw < num_tw_dids; tw++)
        if (tw_dids[tw] != DID_NONE)
            seg_prepare_snip_special_otheri (VCF_SPECIAL_COPY_MATE, ZCTX(tw_dids[tw])->dict_id, copy_mate_snip, tw, 0);

    seg_prepare_snip_special_other (VCF_SPECIAL_LEN_OF, _INFO_HOMSEQ, homlen_snip, 0);
    seg_prepare_snip_special_other (VCF_SPECIAL_LEN_OF, _INFO_DUPHOMSEQ, duphomlen_snip, 0);
    seg_prepare_snip_special_other (VCF_SPECIAL_LEN_OF, _INFO_SVINSSEQ, svinslen_snip, 0);
    seg_prepare_snip_special_other (VCF_SPECIAL_LEN_OF, _INFO_DUPSVINSSEQ, dupsvinslen_snip, 0);
    seg_prepare_snip_other (SNIP_COPY, _INFO_CIPOS, 0, 0, snip_copy_cipos);

    cipos_snip_len = snprintf (cipos_snip, sizeof (cipos_snip), "%.*s%s", STRf(homlen_snip), "0,"); // LEN_OF with "0," prefix
}   

void vcf_sv_seg_initialize (VBlockVCFP vb, Did *tw_dids, int num_tw_dids)
{
    if (tw_dids) { // a known SV caller
        // follow the logic of SAM_QNAME and qname_hash
        decl_ctx (VCF_ID);
        ctx->id_hash.prm8[0] = MIN_(20, MAX_(14, 32 - __builtin_clz (vb->lines.len32 * 5))); // between 14 and 20 bits - tested - no additional compression benefit beyond 20 bits
        buf_alloc_255(vb, &ctx->id_hash, 0, (1ULL << ctx->id_hash.prm8[0]), int32_t, 1, "contexts->id_hash");

        ctx_set_store_per_line (VB, VCF_ID, DID_EOL); // note: VCF_ID is not in the TW list because we don't mux it (we can't mux it as we don't know the mate yet)

        ctx_set_dyn_int (VB, VCF_MATE_POS, DID_EOL);

        for (int tw=0; tw < num_tw_dids; tw++) 
            if (tw_dids[tw] != DID_NONE) {
                ctx_set_store_per_line (VB, tw_dids[tw], DID_EOL);
                seg_mux_init (vb, tw_dids[tw], VCF_SPECIAL_DEMUX_BY_MATE, false, mate[tw]);
            }
    }

    CTX(INFO_CIEND)->flags.same_line = true; // copy CIPOS from same line, regardless if before or after
}

//-------------------
// copying from mate
//-------------------

// in case a field is predicted to be the same as mate's
ContextP vcf_seg_sv_copy_mate (VBlockVCFP vb, ContextP ctx, STRp(value), int my_tw, int her_tw, bool seg_only_if_mated, unsigned add_bytes)
{
    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, ctx->did_i, (MultiplexerP)&vb->mux_mate[my_tw], vcf_has_mate);
    ContextP ret = NULL; // default: caller shouldn't seg 

    if (vcf_has_mate) {
        txtSTR (mate_value, DATA_LINE(vb->mate_line_i)->tw[her_tw]);

        if (str_issame (value, mate_value)) 
            seg_by_ctx (VB, STRi(copy_mate_snip, her_tw), channel_ctx, add_bytes); 
        else 
            goto fallback;
    }

    else fallback: {
        DATA_LINE (vb->line_i)->tw[my_tw] = TXTWORD (value); // we haven't found mate - perhaps our mate will find us...
    
        if (seg_only_if_mated) 
            ret = channel_ctx; // caller will seg - NOT copying mate

        else if (ctx->flags.store == STORE_INT)
            seg_integer_or_not (VB, channel_ctx, STRa(value), add_bytes);
        
        else
            seg_by_ctx (VB, STRa(value), channel_ctx, add_bytes); 
    }

    seg_by_ctx (VB, STRa(vb->mux_mate[my_tw].snip), ctx, 0); // de-multiplexer

    return ret; 
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_DEMUX_BY_MATE)
{
    return reconstruct_demultiplex (vb, ctx, STRa(snip), vcf_has_mate, new_value, reconstruct);
}

static void vcf_piz_sv_recon_ID_from_mate (VBlockVCFP vb, ContextP ctx, bool reconstruct)
{
    vb->mate_line_i = vb->line_i - reconstruct_from_local_int (VB, CTX(VCF_MATE), 0, RECON_OFF);

    if (reconstruct) {
        STR(mate_id);
        recon_history_get_historical_snip (VB, ctx, vb->mate_line_i, pSTRa(mate_id));

        if (ctx->did_i == VCF_ID) {
            switch (segconf.MATEID_method) {
                case MATE_01: 
                    RECONSTRUCT_str (mate_id);
                    *BLSTtxt = (*BLSTtxt == '0') ? '1' : '0';
                    break;

                case MATE_12: 
                    RECONSTRUCT_str (mate_id);
                    *BLSTtxt = (*BLSTtxt == '1') ? '2' : '1';
                    break;

                case MATE_PBSV: {
                    extern StrTextLong vcf_pbsv_get_mate_id (STRp(id), uint32_t *len/*out*/);
                    STR(id);
                    id = vcf_pbsv_get_mate_id (STRa(mate_id), &id_len).s;
                    RECONSTRUCT_str (id);
                    break;
                }
                default:
                    ABORT_PIZ ("Invalid MATEID_method=%u", segconf.MATEID_method);
            }
        }
    }
}

// reconstructs ID from mate, and also sets vb->mate_line_i
SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_COPY_MATE)
{    
    VBlockVCFP vb = (VBlockVCFP)vb_;
    ContextP base_ctx = ctx;

    // optional: base context is different than ctx
    if (snip_len > 1) 
        base_ctx = reconstruct_special_get_base_ctx (VB, ctx, pSTRa(snip));

    // only VCF_ID gets mate from its context, other fields just use vb->mate_line_i
    if (base_ctx->did_i == VCF_ID) {
        if (reconstruct) 
            vcf_piz_sv_recon_ID_from_mate (vb, base_ctx, reconstruct);
        
        return NO_NEW_VALUE;
    }

    // case: string 
    else if (base_ctx->flags.store == STORE_NONE) {
        if (reconstruct) {
            recon_history_get_historical_snip (VB, base_ctx, vb->mate_line_i, pSTRa(snip));
            RECONSTRUCT_snip;
        }
        return NO_NEW_VALUE;
    }
    
    // case: numeric value 
    else if (base_ctx->flags.store == STORE_INT) {
        new_value->i = (vb->mate_line_i < base_ctx->history.len32) ? *B64(base_ctx->history, vb->mate_line_i) : 0; // note: a non-existant STORE_INT mate is taken as 0
        
        if (reconstruct) RECONSTRUCT_INT (new_value->i);
        
        return HAS_NEW_VALUE;
    }

    else if (base_ctx->flags.store == STORE_INDEX) {
        if (reconstruct) {
            WordIndex mate_index = *B(WordIndex, base_ctx->history, vb->mate_line_i);
            ctx_get_snip_by_word_index (base_ctx, mate_index, snip);
            RECONSTRUCT_snip;
        }
            
        return NO_NEW_VALUE;
    }

    else
        ABORT_PIZ ("Invalid flags.store=%s", store_type_name (ctx->flags.store));
}

// ------------------------
// INFO/SVLEN & INFO/REFLEN
// ------------------------

static int64_t vcf_INFO_SVLEN_prediction (VBlockVCFP vb)
{
    if (ALT0(INS) || ALT0(SUBST_INS)) 
        return vb->alt_lens[0] - 1;
    
    if ((ALT0(SYM_DUP) || ALT0(SYM_CNV)) && (IS_PIZ || vb->idx_SVLEN > vb->idx_END)) 
        return CTX(VCF_POS)->last_delta; // note: can predict only if END is before SVLEN

    if (ALT0(SYM_DEL) && (IS_PIZ || vb->idx_SVLEN > vb->idx_END)) 
        return CTX(VCF_POS)->last_delta * (segconf.vcf_del_svlen_is_neg ? -1 : 1); // same comment ^

    if ((ALT0(DEL) || ALT0 (SUBST_DEL) || ALT0(SUBST))) 
        return (int64_t)(vb->REF_len - 1) * (segconf.vcf_del_svlen_is_neg ? -1 : 1);

    return 0x7fffffffffffffff; // unknown
}

void vcf_seg_INFO_SVLEN (VBlockVCFP vb, ContextP ctx, STRp(svlen_str))
{
    int64_t svlen;
    if (!str_get_int (STRa(svlen_str), &svlen)) {
        seg_by_ctx (VB, STRa(svlen_str), ctx, svlen_str_len);
        return;
    }

    int64_t predicted_sv_len = vcf_INFO_SVLEN_prediction (vb);

    if (segconf.running && !segconf.vcf_del_svlen_is_neg && svlen == -predicted_sv_len)
        segconf.vcf_del_svlen_is_neg = true;

    // prediction based on variant type (since 15.0.48)
    if (svlen == predicted_sv_len) {
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, VCF_SPECIAL_SVLEN, 'P' }), 3, ctx, svlen_str_len);
        ctx_set_last_value (VB, ctx, svlen);
    }

    else
        seg_integer_or_not (VB, ctx, STRa(svlen_str), svlen_str_len);
}

void vcf_seg_INFO_REFLEN (VBlockVCFP vb, ContextP ctx, STRp(reflen_str)) // note: ctx is INFO/END *not* POS (despite being an alias)
{
    int64_t reflen;

    if (CTX(VCF_POS)->last_delta && str_get_int (STRa(reflen_str), &reflen) && reflen == CTX(VCF_POS)->last_delta)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_SVLEN, '2' }, 3, ctx, reflen_str_len);
    else
        seg_by_ctx (VB, STRa(reflen_str), ctx, reflen_str_len);
}


// the case where SVLEN is minus the delta between END and POS
SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_SVLEN)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    
    if (!snip_len) 
        new_value->i = -CTX(VCF_POS)->last_delta; // END is a alias of POS - they share the same data stream - so last_delta would be the delta between END and POS

    else if (*snip == 'P')  // introduced 15.0.48
        new_value->i = vcf_INFO_SVLEN_prediction (vb);

    else if (*snip == '2') // introduced 15.0.13 (used for REFLEN ; used for SVLEN up to 15.0.46)
        new_value->i = CTX(VCF_POS)->last_delta;

    else if (*snip == '1') // introduced 15.0.13 ()
        new_value->i = MAX_(vb->ALT_len, vb->REF_len) - 1;

    else
        ABORT_PIZ ("unrecognized snip '%c'(%u). %s", *snip, (uint8_t)*snip, genozip_update_msg());

    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------
// INFO/CIPOS & CIEND
// ----------------------

void vcf_seg_INFO_CIPOS (VBlockVCFP vb, ContextP ctx, STRp(cipos))
{
    str_split_ints (cipos, cipos_len, 2, ',', item, true);

    // sometimes CIPOS is "0,n" where n is the length of HOMSEQ
    if (has(HOMSEQ) && n_items == 2 && items[0] == 0 && items[1] == BII(HOMSEQ)->value_len) 
        seg_by_ctx (VB, STRa(cipos_snip), ctx, cipos_len); // LEN_OF with "0," prefix

    else
        seg_by_ctx (VB, STRa(cipos), ctx, cipos_len);
}

void vcf_seg_INFO_CIEND (VBlockVCFP vb, ContextP ctx, STRp(ciend))
{
    if (has(CIPOS) && str_issame (BII(CIPOS)->value, ciend))
        seg_by_ctx (VB, STRa(snip_copy_cipos), ctx, ciend_len);

    else
        seg_by_ctx (VB, STRa(ciend), ctx, ciend_len);
}

// -----------
// INFO/SVTYPE
// -----------

rom svtype_by_vt[] = SVTYPE_BY_VT; // NULL if vt's SVTYPE is not known (usually bc this vt is not expected to have an SVTYPE)

void vcf_seg_SVTYPE (VBlockVCFP vb, ContextP ctx, STRp(svtype))
{
    VariantType vt = vb->var_types[0];
    bool lowercase_cnv = segconf.vcf_is_pbsv;

    // note in pbsv the SVTYPE for <CNV> is "cnv"
    rom predicted = (vt == VT_SYM_CNV && lowercase_cnv) ? "cnv" : svtype_by_vt[vt];
    uint32_t predicted_len = predicted ? strlen (predicted) : 0;

    // case: we don't expect this VT to have a SVTYPE field
    if (!predicted) fallback:
        seg_by_ctx (VB, STRa(svtype), ctx, svtype_len);

    // case: prediction succeeded
    else if (str_issame (svtype, predicted)) 
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_SVTYPE, segconf.vcf_is_pbsv ? '1' : '0' }, 3, ctx, svtype_len);
    
    // case: SVTYPE field is different than expected
    else
        goto fallback;
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_SVTYPE)
{    
    VBlockVCFP vb = (VBlockVCFP)vb_;

    VariantType vt = vb->var_types[0];
    bool lowercase_cnv = snip_len && snip[0] == '1';

    rom predicted = (vt == VT_SYM_CNV && lowercase_cnv) ? "cnv" : svtype_by_vt[vt];
    uint32_t predicted_len = strlen (predicted);
    
    RECONSTRUCT_str (predicted);

    return NO_NEW_VALUE;
}    

// -----------
// INFO/HOMSEQ
// -----------

// false is not the same or if either goes beyond the end of the range 
static bool is_same_seq (STRp(homseq), 
                         bool mate,            // true - compare mate coordinates, fales - compare our coordinates 
                         bool starting_at_pos, // true compare to reference string starting at POS, false - compare to reference string ending in POS
                         bool revcomp)         // true - compare revcomp of HOMSEQ
{
    return false;

    // xxx TO DO

    // decl_acgt_decode;

    //     RangeP range = IS_REF_EXTERNAL ? ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, homseq_len + 1, WORD_INDEX_NONE, NULL) : NULL;
    //     if (!range) goto fallback;

    // if (!revcomp) {
    //     if (pos + seq_len - 1 > range->last_pos) return false; // seq goes beyond the end of range

    //     for (PosType32 i=0; i < seq_len ; i++) 
    //         if (REFp(pos + i) != seq[i]) return false;
    // }

    // else { 
    //     if (pos < seq_len) return false; // seq goes beyond the start of range

    //     for (PosType32 i=0; i < seq_len ; i++) {
    //         char vcf_base = UPPER_COMPLEM[(int)seq[i]];
    //         char ref_base = REFp (pos - i);
    //         if (ref_base != vcf_base) return false;
    //     }
    // }

    // return true;
}

void vcf_seg_HOMSEQ (VBlockVCFP vb, ContextP ctx, STRp(homseq))
{
    decl_acgt_decode;

    PosType32 pos = DATA_LINE(vb->line_i)->pos;
    char method = 0;

    // method 1: with reference
    if (IS_REF_EXTERNAL && !segconf.vcf_is_svaba) {
        method = '1';

        RangeP range = IS_REF_EXTERNAL ? ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos-homseq_len, homseq_len*2, WORD_INDEX_NONE, NULL) : NULL;
        if (!range) goto fallback;

        for (uint32_t i=0; i < homseq_len; i++) 
            if (REFp (pos + i + 1) != homseq[i]) 
                goto fallback;
    }

    // method 2: Svaba (BND variants only)
    else if (segconf.vcf_is_svaba) {
        ctx = vcf_seg_sv_copy_mate (vb, ctx, STRa(homseq), TW_HOMSEQ, TW_HOMSEQ, true, homseq_len); // set ctx to channel_ctx
        if (!ctx) return; // segged as copy from mate

        if (!IS_REF_EXTERNAL || !ALT0(BND)) goto fallback; // cannot apply method 2
        
        method = '2';

        if ((vb->BND_type == 0 && !is_same_seq (STRa(homseq), true,  true,  false)) ||
            (vb->BND_type == 1 && !is_same_seq (STRa(homseq), false, false, false)) ||
            (vb->BND_type == 2 && !is_same_seq (STRa(homseq), false, true,  false)) || 
            (vb->BND_type == 3 && !is_same_seq (STRa(homseq), true,  false, true )))
            goto fallback;
    }

    // method 3: based on other fields
    else if ((ALT0(DEL) && !str_isprefix_(vb->REF+1, vb->REF_len-1, STRa(homseq))) ||
            (ALT0(INS) && !str_isprefix_(vb->alts[0]+1, vb->alt_lens[0]-1, STRa(homseq))) ||
            (ALT0(SYM_INS) && (!has(LEFT_SVINSSEQ) || !str_isprefix (BII(LEFT_SVINSSEQ)->value, homseq))) ||
            (!ALT0(DEL) && !ALT0(INS) && !ALT0(SYM_INS)))
        goto fallback;
    else
        method = '3';
    
    seg_integer (VB, ctx, homseq_len, false, 0); // store here rather than relying on HOMSEQ_LEN, as HOMSEQ_LEN is not always present
    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_HOMSEQ, method }, 3, ctx, homseq_len);
    return;

fallback:
    seg_by_ctx (VB, STRa(homseq), ctx, homseq_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_HOMSEQ)
{    
    decl_acgt_decode;
    uint32_t homseq_len = reconstruct_from_local_int (vb, ctx, 0, false);
    char method = snip[0];

    // method 1: with reference
    if (method == '1') {
        PosType32 pos = CTX(VCF_POS)->last_value.i;
        ConstRangeP range = ref_piz_get_range (VB, gref, false);

        char *next = BAFTtxt;
        for (uint32_t i=0; i < homseq_len; i++) 
            *next++ = REFp (pos + i + 1);

        Ltxt += homseq_len;
    }

    // method 2: SvABA BNDs based on reference
    else if (method == '2') {
//xxx TODO
    }

    // method 3: based on other fields
    else if (method == '3') {
        if (ALT0(DEL)) 
            RECONSTRUCT (VB_VCF->REF+1, homseq_len);

        else if (ALT0(INS))
            RECONSTRUCT (VB_VCF->alts[0]+1, homseq_len);

        else if (ALT0(SYM_INS)) {
            rom left_svinsseq;
            reconstruct_peek (vb, CTX(INFO_LEFT_SVINSSEQ), &left_svinsseq, NULL);
            RECONSTRUCT (left_svinsseq, homseq_len);
        }
        else
            ABORT_PIZ ("unsupported vartype=%u", VB_VCF->var_types[0]);
    }

    else
        ABORT_PIZ ("Invalid method %u", method);

    return NO_NEW_VALUE;
}    

//----------------
// BND mates
//----------------

#define LINE_BY_HASH(hash) *B(LineIType, CTX(VCF_ID)->id_hash, (hash))

// seg mate as buddy and return true if this line has one 
void vcf_seg_BND_mate (VBlockVCFP vb, STRp(id), 
                       STRp(mate_id), // optional: used if given, if not, test one character difference based on segconf.vcf_mate_id_chars
                       uint64_t hash)
{
    hash &= MAXB(CTX(VCF_ID)->id_hash.prm8[0]);
    LineIType candidate = LINE_BY_HASH(hash);
    
    if (candidate == NO_LINE) goto not_found;

    txtSTR (cand_id, DATA_LINE(candidate)->BND_id);

    // case: test based on the final character changing 0<>1 or 1<>2
    if (!mate_id) {
        char me  = id[id_len-1];
        char her = cand_id[cand_id_len-1];

        rom chars = segconf.MATEID_method == MATE_01 ? "01" : "12";

        // case: mate is found
        if (str_issame_(id, id_len-1, cand_id, cand_id_len-1) &&
            ((me == chars[0] && her == chars[1]) || (me == chars[1] && her == chars[0]))) {
            
            found:
            vb->mate_line_i = candidate;
            vb->mate_line_count++; // for stats
            dyn_int_append (VB, CTX(VCF_MATE), vb->line_i - candidate, 0); // add buddy (delta) >= 1 .
        }
        
        // case: we haven't found our mate - store this line in the hash table so our mate can find us (don't store
        // unneccessarily to reduce hash contention)
        else not_found: {
            DATA_LINE(vb->line_i)->BND_id = TXTWORD(id); 
            LINE_BY_HASH(hash) = vb->line_i;
        }
    }

    // case: test based on given mate_id 
    else {
        if (str_issame (mate_id, cand_id))
            goto found;
        else
            goto not_found;
    }
}

//------------------------------------------------
// copying entire SAMPLES line from mate
//------------------------------------------------

static void vcf_seg_sv_SAMPLES_account (VBlockVCFP vb, STRp(samples), ContextP *ctxs, uint32_t n_ctxs)
{
    str_split (samples, samples_len, vcf_num_samples, '\t', smp, false);

    for (int smp_i=0; smp_i < n_smps; smp_i++) {
        str_split (smps[smp_i], smp_lens[smp_i], n_ctxs, ':', sf, false);

        for (uint32_t sf_i=0; sf_i < n_sfs; sf_i++) 
            ctxs[sf_i]->txt_len += sf_lens[sf_i];
    }
}

ContextP vcf_seg_sv_SAMPLES (VBlockVCFP vb, rom samples, uint32_t remaining_txt_len, ContextP *ctxs, uint32_t n_ctxs)
{
    // case INFO fields depend on sample fields: we can't copy mate
    if (segconf.INFO_DP_method == BY_FORMAT_DP || // INFO/DP depends on FORMAT/DP
        segconf.has[INFO_QD]                   || // INFO/QD sometimes depends on FORMAT/DP
        segconf.AS_SB_TABLE_by_SB              || // INFO/AS_SB_TABLE depends on FORMAT/SB
        CTX(INFO_SF)->sf.SF_by_GT == yes)         // INFO/SF depends on FORMAT/GT
        return CTX(VCF_SAMPLES);

    SAFE_NUL (samples + remaining_txt_len);
    uint32_t samples_len = strcspn (samples, "\n\r");
    SAFE_RESTORE;

    ASSSEG0 (samples_len < remaining_txt_len, "cannot find end of line");

    bool do_account = (vcf_num_samples < 1000); // small enough for str_split
    uint32_t add_bytes = 1/*\n*/ + (do_account ? (str_count_char (STRa(samples), '\t') + str_count_char (STRa(samples), ':')) : samples_len);

    ContextP channel_ctx = vcf_seg_sv_copy_mate (vb, CTX(VCF_SAMPLES), STRa(samples), TW_SAMPLES, TW_SAMPLES, true, add_bytes); 

    if (!channel_ctx) { // segged as copy-mate
        CTX(FORMAT_GT_HT)->use_HT_matrix = false; // we are not going to place GT data in the HT matrix

        if (do_account)
            vcf_seg_sv_SAMPLES_account (vb, STRa(samples), ctxs, n_ctxs);
    }
    
    return channel_ctx;
}

