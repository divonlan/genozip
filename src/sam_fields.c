// ------------------------------------------------------------------
//   sam_fields.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "sam_private.h"
#include "chrom.h"
#include "zip_dyn_int.h"

static const StoreType aux_field_store_flag[256] = {
    ['c']=STORE_INT, ['C']=STORE_INT, 
    ['s']=STORE_INT, ['S']=STORE_INT,
    ['i']=STORE_INT, ['I']=STORE_INT,
    ['f']=STORE_FLOAT
};

static const LocalType aux_field_to_ltype[256] = {
    ['c']=LT_INT8,   ['C']=LT_UINT8, 
    ['s']=LT_INT16,  ['S']=LT_UINT16,
    ['i']=LT_INT32,  ['I']=LT_UINT32,
    ['f']=LT_FLOAT32
};

const char aux_sep_by_type[2][256] = { { // compressing from SAM
    ['c']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['C']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // reconstruct number and \t separator is SAM, and don't reconstruct anything if BAM (reconstruction will be done by translator)
    ['s']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['S']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // -"-
    ['i']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['I']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // -"-
    ['f']=CI0_NATIVE_NEXT | CI0_TRANS_NOR,                                        // compressing SAM - a float is stored as text, and when piz with translate to BAM - is not reconstructed, instead - translated
    ['Z']=CI0_NATIVE_NEXT | CI0_TRANS_NUL, ['H']=CI0_NATIVE_NEXT | CI0_TRANS_NUL, // reconstruct text and then \t separator if SAM and \0 if BAM 
    ['A']=CI0_NATIVE_NEXT,                                                        // reconstruct character and then \t separator if SAM and no separator for BAM
    ['B']=CI0_NATIVE_NEXT                                                         // reconstruct array and then \t separator if SAM and no separator for BAM
}, 
{ // compressing from BAM
    ['c']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['C']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // reconstruct number and \t separator is SAM, and don't reconstruct anything if BAM (reconstruction will be done by translator)
    ['s']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['S']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // -"-
    ['i']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['I']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // -"-
    ['f']=CI0_NATIVE_NEXT,                                                        // compressing SAM - a float is stored as a SPECIAL, and the special reconstructor handles the SAM and BAM reconstructing
    ['Z']=CI0_NATIVE_NEXT | CI0_TRANS_NUL, ['H']=CI0_NATIVE_NEXT | CI0_TRANS_NUL, // reconstruct text and then \t separator if SAM and \0 if BAM 
    ['A']=CI0_NATIVE_NEXT,                                                        // reconstruct character and then \t separator if SAM and no separator for BAM
    ['B']=CI0_NATIVE_NEXT                                                         // reconstruct array and then \t separator if SAM and no separator for BAM
} };

// gets integer field value - may appear before or after this field. 
// return false if not found or not integer or out of [min_value,max_value]
bool sam_seg_peek_int_field (VBlockSAMP vb, Did did_i, int16_t idx, int32_t min_value, int32_t max_value, 
                             bool set_last_value, // set last value if not already set (only if returning true)
                             int32_t *value)      // optional out
{
    decl_ctx (did_i);
    int32_t my_value;

    if (idx == -1)
        return false; // this tag doesn't appear in this alignment

    else if (ctx_has_value_in_line_(VB, ctx)) {
        if (ctx->last_value.i < min_value || ctx->last_value.i > max_value) 
            return false; // note: we test 64 bit variable before truncating to 32 bit

        else {
            if (value) *value = ctx->last_value.i;
            return true;
        }
    }

    else if (idx != -1 && sam_seg_get_aux_int (vb, idx, value ? value : &my_value, IS_BAM_ZIP, min_value, max_value, SOFT_FAIL)) {
        if (set_last_value)
            ctx_set_last_value (VB, ctx, (int64_t)(value ? *value : my_value));
        return true;
    }

    else
        return false;
}

//--------------
// FLAG
//--------------

StrTextLong sam_dis_FLAG (SamFlags f)
{
    StrTextLong result;
    int result_len = 0;

    if (f.multi_segs)     SNPRINTF0 (result, "MultiSeg,");
    if (f.is_aligned)     SNPRINTF0 (result, "Aligned,");
    if (f.unmapped)       SNPRINTF0 (result, "Unmapped,");
    if (f.next_unmapped)  SNPRINTF0 (result, "NextUnmapped,");
    if (f.rev_comp)       SNPRINTF0 (result, "Revcomp,");
    if (f.next_rev_comp)  SNPRINTF0 (result, "NextRevcomp,");
    if (f.is_first)       SNPRINTF0 (result, "First,");
    if (f.is_last)        SNPRINTF0 (result, "Last,");
    if (f.secondary)      SNPRINTF0 (result, "Secondary,");
    if (f.filtered)       SNPRINTF0 (result, "Filtered,");
    if (f.duplicate)      SNPRINTF0 (result, "Duplicate,");
    if (f.supplementary)  SNPRINTF0 (result, "Supplementary,");
    
    if (result_len) result.s[result_len-1] = 0; // removal comma separator from final item

    return result;
}

void sam_seg_FLAG (VBlockSAMP vb, ZipDataLineSAMP dl, unsigned add_bytes)
{    
    // note: for simplicity, we only mux in MAIN.
    // note: we can only mux FLAG by buddy, not mate, bc sam_piz_special_DEMUX_BY_MATE needs FLAG to determine mate
    // note: in collated files, we're better off without this.
    bool do_mux = (IS_MAIN(vb) && segconf.is_paired && !segconf.is_collated);
    ContextP channel_ctx = do_mux ? seg_mux_get_channel_ctx (VB, SAM_FLAG, (MultiplexerP)&vb->mux_FLAG, sam_has_mate || sam_has_saggy) 
                                  : CTX(SAM_FLAG);

    // case: PRIM line: 
    if (IS_PRIM(vb)) {
        if (!IS_SAG_SA) goto normal; // NH-based SA Groups: FLAG will be used in preprocessing to set revcomp, and then also reconstructed normally

        // case: SAG_SA: we store FLAG.rev_comp in OPTION_SA_STRAND instead of FLAG
        sam_seg_against_sa_group_int (vb, CTX(SAM_FLAG), dl->FLAG.value & ~SAM_FLAG_REV_COMP, add_bytes);

        ContextP sa_strand_ctx = CTX(OPTION_SA_STRAND);
        seg_by_ctx (VB, dl->FLAG.rev_comp ? "-" : "+", 1, sa_strand_ctx, 0); 

        // count FLAG field contribution to OPTION_SA_CIGAR, so sam_stats_reallocate can allocate the z_data between CIGAR and SA:Z
        sa_strand_ctx->counts.count++; // contributed z_data due to a single-byte + or -
    }

    // case: DEPN line with SA Group: we know some flags from SA Groups, so we store FLAG without them (reduces FLAG entropy)
    else if (IS_DEPN(vb) && sam_seg_has_sag_by_nonSA(vb)) // First, Last, Multi are the same for the group (see sam_sa_seg_depn_find_sagroup_noSA)
        sam_seg_against_sa_group_int (vb, CTX(SAM_FLAG), dl->FLAG.value & ~(SAM_FLAG_MULTI_SEG | SAM_FLAG_IS_FIRST | SAM_FLAG_IS_LAST), add_bytes);

    else if (IS_DEPN(vb) && sam_seg_has_sag_by_SA(vb)) // in SA, we can get the revcomp from the SA alignment (see sam_sa_seg_depn_find_sagroup_SAtag)
        sam_seg_against_sa_group_int (vb, CTX(SAM_FLAG), dl->FLAG.value & ~(SAM_FLAG_MULTI_SEG | SAM_FLAG_IS_FIRST | SAM_FLAG_IS_LAST | SAM_FLAG_REV_COMP), add_bytes);
    
    // case: depn line in main VB
    // tested this - if snip is not identical to mate's case, the the compression worsens, if snip
    // is identical, we will need to figure out if this is depn from another source - such of a H in CIGAR
    // not doing - too much effort for too small benefit

    // case: we can retrieve the FLAG from this line's mate
    #define SAME_AS_MATE (SAM_FLAG_MULTI_SEG | SAM_FLAG_IS_ALIGNED | SAM_FLAG_SECONDARY | SAM_FLAG_FILTERED | \
                          SAM_FLAG_DUPLICATE | SAM_FLAG_SUPPLEMENTARY)
    else if (sam_has_mate && !IS_PRIM(vb) && // TO DO: improve sam_load_groups_add_flags to allow FLAG mating in prim VB (bug 620)
        ({ ZipDataLineSAMP mate_dl = DATA_LINE (vb->mate_line_i); 
           (dl->FLAG.value & SAME_AS_MATE) == (mate_dl->FLAG.value & SAME_AS_MATE) &&
            dl->FLAG.unmapped              == mate_dl->FLAG.next_unmapped          &&
            dl->FLAG.next_unmapped         == mate_dl->FLAG.unmapped               &&
            dl->FLAG.rev_comp              == mate_dl->FLAG.next_rev_comp          &&
            dl->FLAG.next_rev_comp         == mate_dl->FLAG.rev_comp               &&
            dl->FLAG.is_first              == mate_dl->FLAG.is_last                &&
            dl->FLAG.is_last               == mate_dl->FLAG.is_first; }))
        seg_special0 (VB, SAM_SPECIAL_COPY_MATE_FLAG, channel_ctx, add_bytes); // added 12.0.41

    // case: normal snip
    else normal:
        seg_integer_as_snip_do (VB, channel_ctx, dl->FLAG.value, add_bytes);

    if (do_mux)
        seg_by_did (VB, STRa(vb->mux_FLAG.snip), SAM_FLAG, 0); // de-multiplexer

    // first pairing test (another test is in sam_seg_finalize_segconf)
    if (segconf_running) {
        if (dl->FLAG.is_last && !dl->FLAG.is_first)
            segconf.is_paired = true;
        
        if (sam_is_depn (dl->FLAG))
            segconf.sam_has_depn = true;
    }

    if (dl->FLAG.secondary)     vb->secondary_count++;
    if (dl->FLAG.supplementary) vb->supplementary_count++; 
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_COPY_MATE_FLAG)
{
    if (!VER(14)) sam_piz_set_buddy_v13 (vb); 

    SamFlags flag = { .value = history64(SAM_FLAG, VB_SAM->mate_line_i) }; 

    // flip fields to opposite mate
    SWAPbits (flag.unmapped, flag.next_unmapped);
    SWAPbits (flag.rev_comp, flag.next_rev_comp);
    SWAPbits (flag.is_first, flag.is_last);

    new_value->i = flag.value;
    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE; 
}

// ---------
// BD and BI
// ---------

static void sam_seg_BD_BI_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(field), DictId dict_id, unsigned add_bytes)
{
    bool is_bi = (dict_id.num == _OPTION_BI_Z);
    ContextP ctx  = CTX (is_bi ? OPTION_BI_Z : OPTION_BD_Z);

    // TO DO: we should not enforce BD/BI format correctness, bug 924 (challenge: if either BD or BI are misformatted, then *both* need to be segged normally and dl->BD_BI=(0,0))
    ASSSEG (field_len == dl->SEQ.len, "Expecting %s.len=%u == SEQ.len==%u. %s=\"%.*s\"",
            ctx->tag_name, field_len, dl->SEQ.len, ctx->tag_name, STRf(field));

    dl->BD_BI[is_bi] = BNUMtxt (field);

    CTX(OPTION_BD_BI)->txt_len += add_bytes; 
    CTX(OPTION_BD_BI)->local.len32 += field_len;

    seg_special0 (VB, SAM_SPECIAL_BDBI, ctx, 0);
}

// callback function for compress to get BD_BI data of one line: this is an
// interlaced line containing a character from BD followed by a character from BI - since these two fields are correlated
COMPRESSOR_CALLBACK (sam_zip_BD_BI)
{
    ZipDataLineSAMP dl = DATA_LINE (vb_line_i);

    uint32_t bd_len = dl->BD_BI[0] ? dl->SEQ.len : 0;
    uint32_t bi_len = dl->BD_BI[1] ? dl->SEQ.len : 0;
     
    rom bd = dl->BD_BI[0] ? Btxt (dl->BD_BI[0]) : "";
    rom bi = dl->BD_BI[1] ? Btxt (dl->BD_BI[1]) : "";
    
    if (!bd_len && !bi_len) {
        *line_data_len = 0;
        return;
    }
    
    // both BD and BI must exist. bug 924 
    ASSERT (bd_len && bi_len, "%s/%u Expecting both BD and BI to exist if either exists but BD=\"%.*s\" BI=\"%.*s\"",
            VB_NAME, dl->SEQ.len, STRf(bd), STRf(bi));
        
    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = MIN_(maximum_size, bd_len * 2);
    if (is_rev) *is_rev = dl->FLAG.rev_comp;

    if (!line_data) return; // only length was requested

    ctx = CTX(OPTION_BD_BI);
    buf_alloc_exact (vb, ctx->interlaced, bd_len * 2, uint8_t, "interlaced");

    uint8_t *next = B1ST8 (ctx->interlaced);
    rom after = bd + bd_len;
    while (bd < after) {
        *next++ = *bd;
        *next++ = *bi++ - *bd++;
    }

    *line_data = B1STc (ctx->interlaced);
}   

// BD and BI - reconstruct from BD_BI context which contains interlaced BD and BI data. 
SPECIAL_RECONSTRUCTOR (sam_piz_special_BD_BI)
{
    if (!vb->seq_len || !reconstruct) goto done;

    ContextP bdbi_ctx = CTX(OPTION_BD_BI);

    // note: bd and bi use their own next_local to retrieve data from bdbi_ctx. the actual index
    // in bdbi_ctx.local is calculated given the interlacing
    ASSPIZ (ctx->next_local + vb->seq_len * 2 <= bdbi_ctx->local.len, "Error reading: unexpected end of %s data. Expecting ctx->next_local=%u + vb->seq_len=%u * 2 <= bdbi_ctx->local.len=%"PRIu64, 
            dis_dict_id (bdbi_ctx->dict_id).s, ctx->next_local, vb->seq_len, bdbi_ctx->local.len);

    char *dst   = BAFTtxt;
    char *after = dst + vb->seq_len;
    rom src     = Bc (bdbi_ctx->local, ctx->next_local * 2);

    if (ctx->dict_id.num == _OPTION_BD_Z)
        for (; dst < after; src+=2, dst++) *dst = *src;
    else
        for (; dst < after; src+=2, dst++) *dst = *src + *(src+1);
    
    Ltxt += vb->seq_len;    
    ctx->next_local += vb->seq_len;

done:
    return NO_NEW_VALUE;
}

// ---------------------------------------------
// BQ:Z "Offset to base alignment quality (BAQ)"
// ---------------------------------------------

static void sam_seg_BQ (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(bq), unsigned add_bytes)
{
    ContextP ctx = CTX(OPTION_BQ_Z);

    dl->BQ = TXTWORD (bq);
    ctx->local.len32 += bq_len;

    if (bq_len == dl->SEQ.len) seg_simple_lookup (VB, ctx, add_bytes);
    else                       seg_lookup_with_length (VB, ctx, bq_len, add_bytes);
}

COMPRESSOR_CALLBACK (sam_zip_BQ)
{
    ZipDataLineSAMP dl = DATA_LINE (vb_line_i);
    
    *line_data_len = dl->BQ.len;
    *line_data = Btxt (dl->BQ.index);
}

// ----------------------------------------------------------------------------------------------------------
// MM:Z: Base modifications / methylation
// Example: "C+m,5,12,0;C+h,5,12,0;"
// ----------------------------------------------------------------------------------------------------------

void sam_MM_zip_initialize (void)
{
    segconf.MM_con = (SmallContainer){ 
        .nitems_lo = 2, 
        .repeats   = 1,
        .items[0]  = { .dict_id = sub_dict_id (_OPTION_MM_Z, 'T'), .separator = "," }, // not '0' as it is used by seg_array_by_callback
        .items[1]  = { .dict_id = sub_dict_id (_OPTION_MM_Z, 'N')                   } 
    };

    segconf.MM_con_snip_len = sizeof (segconf.MM_con_snip);
    container_prepare_snip ((ContainerP)&segconf.MM_con, 0, 0, qSTRa(segconf.MM_con_snip));
}

static bool sam_seg_MM_Z_item (VBlockP vb, ContextP ctx, 
                               STRp(mm_item),  // e.g. "C+m,5,12,0" 
                               uint32_t repeat)
{
    rom comma = memchr (mm_item, ',', mm_item_len);

    if (comma) {
        seg_by_dict_id (vb, mm_item, comma - mm_item, segconf.MM_con.items[0].dict_id, comma - mm_item);
        
        rom arr = comma + 1;
        uint32_t arr_len = mm_item_len + mm_item - arr;

        // note: we seg MNM:Z repeats as "special" and copy them from ML:Z if confirmed to be the same
        seg_array_(vb, ctx_get_ctx (vb, segconf.MM_con.items[1].dict_id), ctx->st_did_i, STRa(arr), ',', 0, false, true, sub_dict_id (_OPTION_MM_Z, 'A'), 
                   SAM_SPECIAL_ML_REPEATS, ctx_encountered_in_line (VB, OPTION_ML_B_C) ? CTX(OPTION_ML_B_C)->last_value.i/*ML's n_repeats*/ : -1, // note: will only use special if MM:Z has a single item (otherwise ML:B repeats are distributed between all MM:Z items)
                   arr_len);
    }

    // trivial but legal MM:Z field: "MM:Z:C+m?;"
    else {
        seg_by_dict_id (vb, mm_item, mm_item_len, segconf.MM_con.items[0].dict_id, mm_item_len);
        seg_by_dict_id (vb, NULL, 0,              segconf.MM_con.items[1].dict_id, 0); // Note: NULL rather than "" causes the preceding ',' separator to be deleted by container_reconstruct, allowing us to keep the container intact
    }

    seg_by_ctx (vb, STRa(segconf.MM_con_snip), ctx, !!comma); // account for comma is there is one

    return true;
}

static void sam_seg_MM_Z (VBlockSAMP vb, STRp(mm), unsigned add_bytes)
{
    seg_array_by_callback (VB, CTX(OPTION_MM_Z), STRa(mm), ';', sam_seg_MM_Z_item, 0, 0, add_bytes);
}

// used by container_reconstruct to retrieve the number of repeats of ML:B, and use that for MM:Z subcontext MNM:Z
SPECIAL_RECONSTRUCTOR (sam_piz_special_ML_REPEATS)
{
    new_value->i = CTX(OPTION_ML_B_C)->last_value.i; // number of repeats because ML.store=STORE_INDEX
    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// SM:i: Template-independent mapping quality
// ----------------------------------------------------------------------------------------------------------
static void sam_seg_SM_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t SM, unsigned add_bytes)
{
    decl_ctx (OPTION_SM_i);

    if (IN_RANGX (SM, 0, 255) && 
        SM != 254 &&           // note: 254 is a valid, but highly improbable value - we use 254 for "copy from MAPQ" so a actual 254 is segged as an exception
        !(SM && !dl->MAPQ)) {  // we're expecting SM=0 if MAPQ=0
        
        // store value in local (254 means "copy from MAPQ"), except if MAPQ=0 - we know SM=0 so no need to store
        if (dl->MAPQ) {
            uint8_t SM8 = (SM == dl->MAPQ) ? 254 : SM;
            seg_integer_fixed (VB, ctx, &SM8, false, 0);
        }

        seg_special0 (VB, SAM_SPECIAL_SM, ctx, add_bytes); // this usually is an all-the-same
    }

    else 
        seg_integer_as_snip_do (VB, ctx, SM, add_bytes); 

    ctx_set_last_value (VB, ctx, SM);
}    

SPECIAL_RECONSTRUCTOR (sam_piz_special_SM)
{
    uint8_t MAPQ = CTX(SAM_MAPQ)->last_value.i;

    if (!MAPQ)
        new_value->i = 0;

    else {
        uint8_t value = NEXTLOCAL(uint8_t, ctx);
        new_value->i = (value==254) ? MAPQ : value;
    }

    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// AM:i The smallest template-independent mapping quality in the template
// ----------------------------------------------------------------------------------------------------------
static void sam_seg_AM_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t AM, unsigned add_bytes)
{
    decl_ctx (OPTION_AM_i);

    // note: currently we only support for this algorithm AM appearing after SM. Easily fixable if ever needed.
    // AM is often one of 3 options: 0, =SM =MAPQ-SM. If SM=0 then AM is expected to be 0.
    if (has(SM_i) && 
        IN_RANGX (AM, 0, 255) &&  // valid value
        AM != 253 && AM != 254) { // note: 253,254 are valid, but highly improbable values

        int32_t SM;
        if (!sam_seg_peek_int_field (vb, OPTION_SM_i, vb->idx_SM_i, 0, 255, false, &SM)) goto fallback;

        if (!SM && AM) goto fallback; // usually, AM=0 if SM=0
        
        if (SM) { // no need to store if SM=0, as we know AM=0
            uint8_t AM8 = (AM == SM)                  ? 253 // copy SM (note: AM is not 0, since SM is not 0)
                        : (AM && AM == dl->MAPQ - SM) ? 254 // note subtelty: substraction in uint8_t arithmetics. we are careful to do the same in sam_piz_special_AM.
                        :                               AM;

            seg_integer_fixed (VB, ctx, &AM8, false, 0);
        }
        
        seg_special0 (VB, SAM_SPECIAL_AM, ctx, add_bytes);
    }

    else fallback:
        seg_integer_as_snip_do (VB, ctx, AM, add_bytes);    
}    

SPECIAL_RECONSTRUCTOR (sam_piz_special_AM)
{
    uint8_t MAPQ = CTX(SAM_MAPQ)->last_value.i;

    uint8_t SM = reconstruct_peek (vb, CTX(OPTION_SM_i), 0, 0).i;

    if (!SM)
        new_value->i = 0;

    else {
        uint8_t value = NEXTLOCAL(uint8_t, ctx);
        new_value->i = (value==253) ? SM 
                     : (value==254) ? MAPQ - SM
                     :                value;
    }

    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// UQ:i Phred likelihood of the segment, conditional on the mapping being correct
// ----------------------------------------------------------------------------------------------------------
static void sam_seg_UQ_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t UQ, unsigned add_bytes)
{
    decl_ctx (OPTION_UQ_i);
    int32_t other; // either AS or NM

    // in Novoalign, usually UQ==AS
    if (MP(NOVOALIGN) && sam_seg_peek_int_field (vb, OPTION_AS_i, vb->idx_AS_i, 0, 1000, true, &other) && other == UQ)
        seg_special1 (VB, SAM_SPECIAL_UQ, '1', ctx, add_bytes);

    // In GATK produced data, in many cases (~half) UQ==2*NM
    else if (!MP(NOVOALIGN) && sam_seg_peek_int_field (vb, OPTION_NM_i, vb->idx_NM_i, 0, 0x3fffffff, 
             false, // bc seg_NM_isam_seg_NM_i inteprets would interpret a set value as "segged already"
             &other) && other*2 == UQ)
        seg_special1 (VB, SAM_SPECIAL_UQ, '2', ctx, add_bytes);

    else 
        seg_integer (VB, ctx, UQ, true, add_bytes);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_UQ)
{
    ASSPIZ (*snip=='1' || *snip=='2', "Unrecognized snip. %s", genozip_update_msg());

    int64_t other = reconstruct_peek (vb, CTX(*snip=='1' ? OPTION_AS_i : OPTION_NM_i), 0, 0).i; 

    new_value->i = other * (*snip - '0');

    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// NH:i Number of reported alignments that contain the query in the current record
// ----------------------------------------------------------------------------------------------------------

static inline void sam_seg_NH_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t nh, unsigned add_bytes)
{
    decl_ctx (OPTION_NH_i);
        
    dl->NH = nh; 
    ctx_set_encountered (VB, ctx);  // needed by sam_seg_HI_i

    // build SAG structure in VB, to be later ingested into z_file->sag_*
    if (IS_DEPN(vb) && vb->sag && vb->sag->num_alns == nh) {
        sam_seg_against_sa_group (vb, ctx, add_bytes);
        add_bytes = 0; // if this is prim, we seg again below, but we already added the bytes
    }

    else
        sam_seg_buddied_i_fields (vb, dl, OPTION_NH_i, nh, &dl->NH, (MultiplexerP)&vb->mux_NH, 
                                  STRa (copy_buddy_NH_snip), add_bytes);
}

// ----------------------------------------------------------------------------------------------------------
// HI:i Query hit index (a number [1,NH])
// ----------------------------------------------------------------------------------------------------------

SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_BY_BUDDY_MAP)
{
    int channel_i = (last_flags.unmapped?2 : (sam_has_mate || sam_has_saggy)?1 : 0);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

static inline void sam_seg_HI_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t hi, unsigned add_bytes)
{
    decl_ctx (OPTION_HI_i);

    if (hi >= 2 && segconf.running)
        segconf.HI_has_two_plus = true;

    if (segconf.is_collated) {
        int64_t last = ctx->last_value.i;
        int prediction = (ctx_encountered_in_line (VB, OPTION_NH_i) && dl->NH==0) ? 0         // NH=0 -> HI=0
                       : (!vb->line_i || last == (dl-1)->NH)                      ? 1         // first in new SA Group 
                       :                                                            last + 1; // next HI in current SA Group

        if (prediction == hi)
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_HI }, 2, ctx, add_bytes);

        else  
            seg_integer (VB, ctx, hi, true, add_bytes);
    }
    else  // not collated
        seg_integer (VB, ctx, hi, false, add_bytes);

    ctx_set_last_value (VB, ctx, hi);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_HI)
{
    ContextP nh_ctx = CTX(OPTION_NH_i);
    int64_t last_nh = nh_ctx->last_value.i;
    int64_t last_hi = ctx->last_value.i;

    new_value->i = (ctx_encountered_in_line (VB, OPTION_NH_i) && last_nh == 0)            ? 0 // NH=0 -> HI=0
                 : (!vb->line_i || last_hi == *B(int64_t, nh_ctx->history, vb->line_i-1)) ? 1 // first in new SA Group
                 :                                                                          last_hi + 1; // next HI in current SA Group

    if (reconstruct)
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------
// CP:i Leftmost coordinate of the next hit
// ----------------------------------------

static inline void sam_seg_CP_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t cp, unsigned add_bytes)
{
    if (IS_DEPN(vb) && vb->sag && IS_SAG_CC && vb->cc_aln->pos == cp) 
        sam_seg_against_sa_group (vb, CTX(OPTION_CP_i), add_bytes);

    // if there's no gencomp - CP tends to compress well by delta, as subsequent lines also tend to have subsequent secondaries
    else if (IS_MAIN(vb) && segconf.is_sorted) 
        seg_self_delta (VB, CTX(OPTION_CP_i), cp, 0, 0, add_bytes); 

    else {
        uint32_t cp32 = cp;
        seg_integer (VB, CTX(OPTION_CP_i), cp32, true, add_bytes);
    }
}

// -------------------------
// OA:Z "Original alignment"
// -------------------------

bool sam_seg_0A_rname_cb (VBlockP vb, ContextP ctx, STRp(oa_rname), uint32_t repeat)
{
    if (str_issame(oa_rname, vb->chrom_name)) // note: 0A are merely ALIAS_DICT of RNAME - they have their own b250 so this is beneficial
        seg_special0 (vb, SAM_SPECIAL_COPY_RNAME, ctx, oa_rname_len);

    else
        chrom_seg_ex (vb, ctx->did_i, STRa(oa_rname), 0, oa_rname_len, NULL);
    
    return true; // segged successfully
}

// copy from RNAME (note: we always recon textual, unlike RNAME that is reconstructed binary in BAM)
SPECIAL_RECONSTRUCTOR (sam_piz_special_COPY_RNAME)
{
    if (!reconstruct) goto done;

    ContextP rname_ctx = CTX(SAM_RNAME);
    ASSPIZ0 (ctx_has_value_in_line_(vb, rname_ctx), "No value for RNAME in line");

    ctx_get_snip_by_word_index (rname_ctx, rname_ctx->last_value.i, snip);
    
    RECONSTRUCT_snip;

done:
    return NO_NEW_VALUE;
}

// also called for YA:Z and YO:Z
bool sam_seg_0A_pos_cb (VBlockP vb, ContextP ctx, STRp(oa_pos), uint32_t repeat)
{
    seg_pos_field (vb, ctx->did_i, SAM_POS, 0, 0, STRa(oa_pos), 0, oa_pos_len);
    return true; // segged successfully
}

// OA is: (rname, pos, strand, CIGAR, mapQ, NM ;)+ . NM is optional (but its , is not)
// Example OA:Z:chr13,52863337,-,56S25M70S,0,0;chr6,145915118,+,97S24M30S,0,0;chr18,64524943,-,13S22M116S,0,0;chr7,56198174,-,20M131S,0,0;chr7,87594501,+,34S20M97S,0,0;chr4,12193416,+,58S19M74S,0,0;
// See: https://samtools.github.io/hts-specs/SAMtags.pdf
static void sam_seg_OA_Z (VBlockSAMP vb, STRp(field), unsigned add_bytes)
{
    static const MediumContainer container_OA = { .nitems_lo = 6,          
                                                  .repsep    = { ';' }, // including on last repeat    
                                                  .items     = { { .dict_id = { _OPTION_OA_RNAME  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_POS    }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_STRAND }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_CIGAR  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_MAPQ   }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_NM     },                    } } };

    SegCallback callbacks[6] = { [SA_RNAME]=sam_seg_0A_rname_cb, [SA_POS]=sam_seg_0A_pos_cb, [SA_CIGAR]=sam_seg_0A_cigar_cb, [SA_MAPQ]=sam_seg_0A_mapq_cb };
     
    seg_array_of_struct (VB, CTX(OPTION_OA_Z), container_OA, STRa(field), callbacks, 
                         segconf.sam_semcol_in_contig ? sam_seg_correct_for_semcol_in_contig : NULL,
                         add_bytes);
}

// ----------------------------------------------------------------------------------------------
// AS:i "Alignment score generated by aligner" (https://samtools.github.io/hts-specs/SAMtags.pdf)
// ----------------------------------------------------------------------------------------------

// AS has a value set (at least as set by BWA and IonTorrent TMAP) of at most vb->ref_consumed, and often equal to it. we modify
// it to be new_value=(value-ref_consumed) 
static inline void sam_seg_AS_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t as, unsigned add_bytes)
{
    START_TIMER;

    if (segconf_running) {
        if (ABS((int32_t)vb->ref_consumed     - as) < 20) segconf.AS_is_ref_consumed++;
        if (ABS((int32_t)vb->ref_consumed * 2 - as) < 10) segconf.AS_is_2ref_consumed++;
    }

    // note: dl->AS was already set in sam_seg_txt_line/bam_seg_txt_line
    ctx_set_last_value (VB, CTX (OPTION_AS_i), as);

    // depn VB - seg against prim line AS (sag_has_AS determines if its beneficial to do so)
    if (IS_DEPN(vb) && segconf.sag_has_AS && vb->sag) 
        sam_seg_against_sa_group_int (vb, CTX(OPTION_AS_i), (int64_t)vb->sag->as - as, add_bytes);

    // in bowtie2-like data, we might be able to copy from mate
    else if (segconf.is_bowtie2) {
        ASSERT (IN_RANGX (as, MIN_AS_i, MAX_AS_i), "%s: AS=%"PRId64" is ∉ [%d,%d]", LN_NAME, as, MIN_AS_i, MAX_AS_i);    
        
        ZipDataLineSAMP mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is -1

        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, OPTION_AS_i, (MultiplexerP)&vb->mux_AS, sam_has_mate);

        if (sam_has_mate && segconf.has[OPTION_YS_i] && mate_dl->YS == as) 
            seg_by_ctx (VB, STRa(copy_mate_YS_snip), channel_ctx, add_bytes);
        
        // case: AS:i tends to be close to 2 X ref_consumed
        else if (segconf.AS_is_2ref_consumed && !segconf_running && ABS((int32_t)2*vb->ref_consumed - as) < 10) {
            SNIPi3 (SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED, 'x', (int32_t)2*vb->ref_consumed - as);
            seg_by_ctx (VB, STRa(snip), channel_ctx, add_bytes); 
        }

        else
            // TODO: AS prediction, see bug 520
            //seg_integer_as_snip_do (VB, channel_ctx, as, add_bytes);    
            seg_integer (VB, channel_ctx, as, true, add_bytes);

        seg_by_did (VB, STRa(vb->mux_AS.snip), OPTION_AS_i, 0); // de-multiplexer
    }

    // case: in STAR paired files, AS expected to be the same value as its mate's
    else if (MP(STAR) && segconf.is_paired && !IS_DEPN(vb)) 
        sam_seg_buddied_i_fields (vb, dl, OPTION_AS_i, as, &dl->AS, (MultiplexerP)&vb->mux_AS, STRa(copy_mate_AS_snip), add_bytes);

    // not bowtie2: store a special snip with delta from ref_consumed
    else if (segconf.AS_is_ref_consumed) {
        SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED, (int32_t)vb->ref_consumed-as);
        seg_by_did (VB, STRa(snip), OPTION_AS_i, add_bytes); 
    }

    // if we have minimap2-produced ms:i, AS is close to it
    else if (segconf.sam_ms_type == ms_MINIMAP2 && sam_seg_peek_int_field (vb, OPTION_ms_i, vb->idx_ms_i, -0x8000000, 0x7fffffff, true/*needed for delta*/, NULL)) 
        seg_delta_vs_other_localN (VB, CTX(OPTION_AS_i), CTX(OPTION_ms_i), as, -1, add_bytes);

    else
        seg_integer (VB, CTX (OPTION_AS_i), as, true, add_bytes);

    COPY_TIMER(sam_seg_AS_i);
}

// ----------------------------------------------------------------------------------------------
// delta vs seq_len. Used for Ultima a3:i, minimap2 ms:i and until 12.0.37 also AS:i 
// ----------------------------------------------------------------------------------------------

static void sam_seg_delta_seq_len (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, int64_t value, bool in_local, unsigned add_bytes)    
{
    decl_ctx (did_i);

    if (in_local) { // since 15.0.61
        seg_special1 (VB, SAM_SPECIAL_delta_seq_len, '$'/*=in_local*/, ctx, add_bytes);
        dyn_int_append (VB, ctx, (int64_t)dl->SEQ.len - value, 0); 
    }

    else if (!value)
        seg_special0 (VB, SAM_SPECIAL_delta_seq_len, ctx, add_bytes);

    else {
        SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_delta_seq_len, (int64_t)dl->SEQ.len - value);
        seg_by_did (VB, STRa(snip), did_i, add_bytes);
    }
}

// reconstruct seq_len or (seq_len-snip)
SPECIAL_RECONSTRUCTOR (sam_piz_special_delta_seq_len)
{
    int64_t delta;
    if (str_is_1char (snip, '$')) // since 15.0.61
        delta = reconstruct_from_local_int (vb, ctx, 0, false);    
    else
        delta = atoi (snip); // 0 if snip==""

    new_value->i = (int32_t)vb->seq_len - delta; 
    
    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// reconstruct ref_consumed or (ref_consumed-snip)
SPECIAL_RECONSTRUCTOR (sam_piz_special_REF_CONSUMED)
{
    if (snip_len && *snip == 'x')
        new_value->i = 2 * (int32_t)VB_SAM->ref_consumed - atoi (snip+1); // snip="" is 0
    else
        new_value->i =     (int32_t)VB_SAM->ref_consumed - atoi (snip); 

    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------
// Prediction based on SEQ.len and soft_clip. Use for:
// qs:i (pacbio) : per-read: 0-based start of query in the ZMW read (absent in CCS)
// qe:i (pacbio) : per-read: 0-based end of query in the ZMW read (absent in CCS)
// XS:i (blasr)  : read alignment start position without counting previous soft clips (1 based)
// XE:i (blasr)  : read alignment end position without counting previous soft clips (1 based)
// XL:i (blasr)  : aligned read length
// XQ:i (blasr)  : query read length
// QS:i (ngmlr)  : query start
// QE:i (ngmlr)  : query end
// XR:i (ngmlr)  : query length minus soft clips
// ----------------------------------------------------------------------------------------------

static inline int32_t sam_SEQ_END_prediction (rom op, int32_t seq_len, int32_t soft0, int32_t soft1, bool revcomp) // all signed
{
    return (op[0]=='+'             ? seq_len : 0)
         + (op[1]=='+'             ? 1       : 0)
         + (op[2]=='+'             ? soft0   : 0)    
         + (op[2]=='F' && !revcomp ? soft0   : 0)  // like +, but executed only if !revcomp (F and R introduced 15.0.69)
         + (op[2]=='R' &&  revcomp ? soft0   : 0)    
         - (op[2]=='-'             ? soft0   : 0)    
         + (op[3]=='+'             ? soft1   : 0)    
         + (op[3]=='F' && !revcomp ? soft1   : 0)    
         + (op[3]=='R' &&  revcomp ? soft1   : 0)    
         - (op[3]=='-'             ? soft1   : 0);
} 

static inline void sam_seg_SEQ_END (VBlockSAMP vb, ZipDataLineSAMP dl, ContextP ctx, int64_t value, rom op, unsigned add_bytes)
{
    int32_t prediction = sam_SEQ_END_prediction (op, dl->SEQ.len, vb->soft_clip[0], vb->soft_clip[1], dl->FLAG.rev_comp);

    if (value == prediction)
        seg_special4 (VB, SAM_SPECIAL_SEQ_LEN, op[0], op[1], op[2], op[3], ctx, add_bytes); // note: this SPECIAL is also used in sam_seg_MD_Z

    else 
        seg_integer (VB, ctx, value, true, add_bytes);
}

// reconstruct a modified dl->SEQ.len based on a 4-parameter snip
SPECIAL_RECONSTRUCTOR (sam_piz_special_SEQ_LEN)
{
    new_value->i = sam_SEQ_END_prediction (snip, CTX(SAM_SQBITMAP)->last_txt.len, VB_SAM->soft_clip[0], VB_SAM->soft_clip[1], last_flags.rev_comp);

    if (reconstruct) 
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------
// MAPQ and MQ:i (mate MAPQ)
// ----------------------------------------------------------------------------------------------

// We seg against a previous mate line's MQ if one exists, but not if this is a single-MAPQ-value file
void sam_seg_MAPQ (VBlockSAMP vb, ZipDataLineSAMP dl, unsigned add_bytes)
{
    if (segconf_running && dl->MAPQ) {
        if (!segconf.MAPQ_value) 
            segconf.MAPQ_value = dl->MAPQ;
        else if (segconf.MAPQ_value != dl->MAPQ) 
            segconf.MAPQ_has_single_value = false;
    }

    ctx_set_last_value (VB, CTX(SAM_MAPQ), (int64_t)dl->MAPQ);

    bool do_mux = IS_MAIN(vb) && !segconf.sam_is_unmapped; // MAIN-only for simplicity. To do: also for prim/depn components
    int channel_i = (segconf.MAPQ_use_xq && has(xq_i))         ? 3  // DRAGEN: if xq:i exists in the alignment, MAPQ is expected to be all-the-same
                  : (segconf.MAPQ_use_xq && has(XQ_i))         ? 3  // DRAGEN: if XQ:i exists in the alignment, MAPQ is expected to be all-the-same
                  : (segconf.has[OPTION_MQ_i] && sam_has_mate) ? 1  // MAPQ is expected to be == mate's MQ 
                  : sam_has_prim                               ? 2  // seg saggy depn to separate channel as expected to be lower MAPQ
                  :                                              0;

    ContextP channel_ctx = do_mux ? seg_mux_get_channel_ctx (VB, SAM_MAPQ, (MultiplexerP)&vb->mux_MAPQ, channel_i) 
                                  : CTX(SAM_MAPQ);

    ZipDataLineSAMP mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is -1

    if (segconf.sam_is_unmapped) 
        seg_integer_as_snip_do (VB, channel_ctx, dl->MAPQ, add_bytes); // expecting all-the-same
    
    else if (sam_seg_has_sag_by_SA (vb)) {
        sam_seg_against_sa_group (vb, channel_ctx, add_bytes);

        // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
        if (IS_PRIM(vb)) {
            ContextP sa_mapq_ctx = CTX(OPTION_SA_MAPQ);
            seg_integer_fixed (VB, sa_mapq_ctx, &dl->MAPQ, false, 0);

            // count MAPQ field contribution to OPTION_SA_MAPQ, so sam_stats_reallocate can allocate the z_data between MAPQ and SA:Z
            sa_mapq_ctx->counts.count += add_bytes; 
        }
    }

    else if (channel_i == 3) 
        seg_integer_as_snip_do (VB, channel_ctx, dl->MAPQ, add_bytes); // expecting all-the-same

    // case: seg against mate
    else if (!segconf.MAPQ_has_single_value && segconf.has[OPTION_MQ_i] &&
             sam_has_mate && dl->MAPQ == mate_dl->MQ)
        seg_by_ctx (VB, STRa(copy_mate_MQ_snip), channel_ctx, add_bytes); // copy MQ from earlier-line mate 

    // // case: seg against predicted alignment in prim line SA:Z 
    // else if (!segconf.MAPQ_has_single_value && zip_has_prim && 
    //          dl->MAPQ && dl->MAPQ != (dl-1)->MAPQ &&
    //          sam_seg_is_item_predicted_by_prim_SA (vb, SA_MAPQ, dl->MAPQ)) 
    //     seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_COPY_PRIM, '0'+SA_MAPQ }, 3, channel_ctx, add_bytes); 

    else {
        channel_ctx->ltype = LT_UINT8;
        seg_integer_fixed (VB, channel_ctx, &dl->MAPQ, true, add_bytes);
    }

    if (do_mux)
        seg_by_did (VB, STRa(vb->mux_MAPQ.snip), SAM_MAPQ, 0); // de-multiplexer
}

// since 15.0.61. until then MAPQ used DEMUX_BY_MATE_PRIM
SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_MAPQ)
{
    int channel_i;
    if (segconf.MAPQ_use_xq && (
              container_peek_has_item (vb, CTX(SAM_AUX), _OPTION_xq_i, false) || // since 15.0.61
              (VER2 (15,69) && container_peek_has_item (vb, CTX(SAM_AUX), _OPTION_XQ_i, false))))  // since 15.0.69
        channel_i = 3; // since 15.0.61
    
    else {
        // note: when reconstructing MAPQ in BAM, FLAG is not known yet, so we peek it here
        if (!ctx_has_value_in_line_(vb, CTX(SAM_FLAG)))
            ctx_set_last_value (vb, CTX(SAM_FLAG), reconstruct_peek (vb, CTX(SAM_FLAG), 0, 0));

        channel_i = (segconf.has[OPTION_MQ_i] && sam_has_mate)?1 : sam_has_prim?2 : 0;
    }

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

// MQ:i Mapping quality of the mate/next segment
// Seg against mate if we have one, or else against MAPQ as it is often very similar
static inline void sam_seg_MQ_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t mq, unsigned add_bytes)
{
    ASSERT (IN_RANGX (mq, 0, 255), "%s: Invalid MQ:i=%"PRId64": expecting an integer [0,255]", LN_NAME, mq);
    dl->MQ = mq; 
    
    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, OPTION_MQ_i, (MultiplexerP)&vb->mux_MQ, sam_has_mate);

    if (sam_has_mate && 
        dl->MQ == DATA_LINE (vb->mate_line_i)->MAPQ) 
        seg_by_ctx (VB, STRa(copy_mate_MAPQ_snip), channel_ctx, add_bytes); // copy MAPQ from earlier-line mate 

    else 
        seg_delta_vs_other_localN (VB, channel_ctx, CTX(SAM_MAPQ), mq, -1, add_bytes);

    seg_by_did (VB, STRa(vb->mux_MQ.snip), OPTION_MQ_i, 0); // de-multiplexer
}

// PQ:i Phred likelihood of the template, conditional on the mapping locations of both/all segments being correct.
static inline void sam_seg_PQ_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t pq, unsigned add_bytes)
{
    if (IN_RANGX (pq, 0, 65534)) // dl->PQ is uint16_t
        dl->PQ = pq + 1; // +1, so that if pq is out of this range, leave dl as 0, which will mean "no valid PQ"
    
    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, OPTION_PQ_i, (MultiplexerP)&vb->mux_PQ, sam_has_mate);

    uint16_t mate_pq; 
    if (sam_has_mate && (mate_pq = DATA_LINE (vb->mate_line_i)->PQ) && pq+1 == mate_pq) // this implies 0<=pq<=65534
        seg_by_ctx (VB, STRa(copy_mate_PQ_snip), channel_ctx, add_bytes); // copy MAPQ from earlier-line mate 

    else if (has(AS_i)) {
        CTX(OPTION_AS_i)->last_value.i = dl->AS;
        seg_delta_vs_other_localN (VB, channel_ctx, CTX(OPTION_AS_i), pq, -1, add_bytes);
    }

    else
        seg_integer (VB, channel_ctx, pq, true, add_bytes);

    seg_by_did (VB, STRa(vb->mux_PQ.snip), OPTION_PQ_i, 0); // de-multiplexer
}

// RG:Z, PG:Z, PU:Z, LB:Z, RX:Z and others: we predict that the value will be the same as the buddy
void sam_seg_buddied_Z_fields (VBlockSAMP vb, ZipDataLineSAMP dl, MatedZFields f, STRp(value), 
                               SegBuddiedCallback seg_cb, // optional
                               unsigned add_bytes)
{
    ContextP ctx = CTX(buddied_Z_dids[f]);

    ZipDataLineSAMP buddy_dl = sam_has_mate  ? DATA_LINE (vb->mate_line_i)  // mate before saggy, consistent with sam_piz_special_COPY_BUDDY
                             : sam_has_saggy ? DATA_LINE (vb->saggy_line_i) // if collated, it will likely be the same as prev line, so better not use special 
                             : NULL;

    bool do_mux = segconf.is_paired && !segconf.is_collated; // muxing just on prim if not paired is counter-beneficial

    ContextP channel_ctx = do_mux ? seg_mux_get_channel_ctx (VB, buddied_Z_dids[f], (MultiplexerP)&vb->mux_mated_z_fields[f], !!buddy_dl)
                                  : ctx;

    if (do_mux && buddy_dl && str_issame_(STRa(value), STRline (buddy_dl, mated_z_fields[f]))) 
        seg_by_ctx (VB, STRi(copy_buddy_Z_snip,f), channel_ctx, add_bytes); 

    else if (seg_cb && !segconf_running) 
        seg_cb (vb, channel_ctx, STRa(value), add_bytes);

    else 
        seg_by_ctx (VB, STRa(value), channel_ctx, add_bytes); 

    if (do_mux)
        seg_by_did (VB, STRa(vb->mux_mated_z_fields[f].snip), buddied_Z_dids[f], 0); // de-multiplexer

    set_LineWord_str (dl, mated_z_fields[f], value);
}

void sam_seg_buddied_i_fields (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, 
                               int64_t my_value, 
                               int32_t *dl_value, // pointer to my dl value (value will be set)
                               MultiplexerP mux,
                               STRp(copy_snip),   // buddy_type embedded in snip needs to be consistent with mux->special_code
                               unsigned add_bytes)
{
    decl_ctx (did_i);

    ASSERT (ctx->flags.store_per_line || ctx->flags.spl_custom || segconf_running,  
            "%s: expecting ctx=%s to have store_per_line=true or spl_custom=true", LN_NAME, ctx->tag_name);

    // BAM spec permits values up to 0xffffffff, and SAM is unlimited, however for code covenience we limit
    // values segged with this method to int32_t. If this is ever an issue, it can be solved.
    ASSERT (IN_RANGX (my_value, -0x80000000LL, 0x7fffffffLL), "%s: Value of %s is %"PRId64", outside the supported range by Genozip of [%d,%d]",
            LN_NAME, ctx->tag_name, my_value, -0x80000000, 0x7fffffff);

    #define by_mate      (mux->special_code == SAM_SPECIAL_DEMUX_BY_MATE)
    #define by_buddy     (mux->special_code == SAM_SPECIAL_DEMUX_BY_BUDDY)
    #define by_buddy_map (mux->special_code == SAM_SPECIAL_DEMUX_BY_BUDDY_MAP)
    #define by_mate_prim (mux->special_code == SAM_SPECIAL_DEMUX_BY_MATE_PRIM)

    if (!segconf_running && segconf.is_paired && !segconf.is_collated) {
        int channel_i = by_mate      ? sam_has_mate
                      : by_buddy     ? (sam_has_mate || sam_has_saggy)
                      : by_buddy_map ? (dl->FLAG.unmapped?2 : (sam_has_mate || sam_has_saggy)?1 : 0)
                      : by_mate_prim ? (sam_has_mate?1 : sam_has_prim?2 : 0)
                      : -1; // invalid

        ASSERT (channel_i >= 0, "mux of %s has unsupported special code %u", CTX(did_i)->tag_name, mux->special_code);

        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, did_i, mux, channel_i);

        // calculating mate value's pointer. Note: this only works if ZipDataLineSAM is NOT packed
        int32_t *mate_value  = (vb->mate_line_i  >= 0) ? &dl_value[(vb->mate_line_i  - vb->line_i) * sizeof (ZipDataLineSAM) / sizeof (int32_t)] : NULL;
        int32_t *saggy_value = (vb->saggy_line_i >= 0) ? &dl_value[(vb->saggy_line_i - vb->line_i) * sizeof (ZipDataLineSAM) / sizeof (int32_t)] : NULL;

        // note: we don't check if our mate has this field. the math still works, because if it doesn't have 
        // this field, the dl value in Seg will be 0, as will be the History value in recon.

        if (by_buddy_map && dl->FLAG.unmapped)
            seg_integer_as_snip_do (VB, channel_ctx, my_value, add_bytes); // seg as snip, as this will likely become all-the-same

        else if (sam_has_mate && my_value == *mate_value)    
            seg_by_ctx (VB, STRa(copy_snip), channel_ctx, add_bytes); // note: prior to v14, we stored ms qual_score in QUAL history, not in ms:i history 
        
        else if (((by_buddy && sam_has_saggy) || (by_mate_prim && sam_has_prim)) && my_value == *saggy_value)
            seg_by_ctx (VB, STRa(copy_snip), channel_ctx, add_bytes); 

        else 
            seg_integer (VB, channel_ctx, my_value, true, add_bytes);    

        seg_by_ctx (VB, MUX_SNIP(mux), mux->snip_len, ctx, 0); // de-multiplexer
    }
    else
        seg_integer (VB, ctx, my_value, true, add_bytes);        
}

// seg Z field which is expected to be equal to a different field on mate (eg rb:Z <> mb:Z)
static void sam_seg_cross_mated_Z_fields (VBlockSAMP vb, Did did_i, ZipDataLineSAMP dl, STRp(value), 
                                          LineWordL *my_lw, // to be set
                                          const LineWordL *mate_lw, // pointer to CURRENT LINE value of other field in dl
                                          const Multiplexer2 *mux,
                                          STRp(copy_snip),  
                                          unsigned add_bytes)
{
    *my_lw = LINEWORDL (dl, value);

    ZipDataLineSAMP mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is -1

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, did_i, (MultiplexerP)mux, sam_has_mate);
    if (!sam_has_mate) goto no_mate;

    mate_lw = (LineWordL *)((rom)mate_lw - (rom)dl + (rom)mate_dl); // update to actual line
    rom mate_value = Btxt (mate_dl->line_start + mate_lw->index);

    if (str_issame_(mate_value, mate_lw->len, STRa(value)))
        seg_by_ctx (VB, STRa(copy_snip), channel_ctx, add_bytes);

    else no_mate:
        seg_by_ctx (VB, STRa(value), channel_ctx, add_bytes);

    seg_by_did (VB, STRa(mux->snip), did_i, 0); // de-multiplexer
}

// E2 - SEQ data. Currently broken. To do: fix (bug 403)
/*static void sam_seg_E2_field (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(field), unsigned add_bytes)
{
    ASSSEG0 (dl->SEQ.len, "E2 tag without a SEQ"); 
    ASSINP (field_len == dl->SEQ.len, 
            "Error in %s: Expecting E2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
            txt_name, dl->SEQ.len, field_len, field_len, field);

    PosType32 this_pos = vb->last_int(SAM_POS);

    sam_seg_SEQ (vb, OPTION_E2_Z, (char *)STRa(field), this_pos, vb->last_cigar, vb->ref_consumed, vb->ref_and_seq_consumed, 0, field_len, // remove const bc SEQ data is actually going to be modified
                        vb->last_cigar, add_bytes); 
}*/

// U2 - QUAL data
static void sam_seg_U2_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(field), unsigned add_bytes)
{
    ASSSEG0 (dl->SEQ.len, "U2 tag without a SEQ"); 
    ASSINP (field_len == dl->SEQ.len, 
            "%s: Expecting U2 data to be of length %u as indicated by CIGAR, but it is %u. U2=%.*s",
            LN_NAME, dl->SEQ.len, field_len, STRf(field));

    dl->U2 = BNUMtxt(field);
    CTX(OPTION_U2_Z)->txt_len   += add_bytes;
    CTX(OPTION_U2_Z)->local.len += field_len;
}

static void sam_seg_RG_Z (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(rg), unsigned add_bytes)
{
    decl_ctx (OPTION_RG_Z);
    
    // this pattern was observed in CellRanger files, but we don't limit it to only CellRanger 
    if (segconf.RG_method == RG_BY_ILLUM_QNAME ||
        (segconf_running && segconf.tech == TECH_ILLUM && segconf.sam_multi_RG)) {
        STRlast (qname, SAM_QNAME);
        int64_t wi_plus_1;
        if (!str_item_i_int (qname, qname_len, ':', 3, &wi_plus_1)) // note: we use str_item_i and not Q3NAME.last_value because QNAME might be segged by copy buddy, and different Illumina flavors have the RG in different items
            goto fallback;

        if (wi_plus_1 < 1 || wi_plus_1 > ctx->ol_nodes.len)
            goto fallback;
            
        STR0(snip);
        ctx_get_vb_snip_ex (ctx, wi_plus_1-1, pSTRa(snip));

        if (!str_issame (snip, rg)) goto fallback;

        seg_special1 (VB, SAM_SPECIAL_RG_by_QNAME, '0' + segconf.RG_method, ctx, add_bytes);
    }
    
    else fallback: {
        if (segconf.sam_multi_RG)
            sam_seg_buddied_Z_fields (vb, dl, MATED_RG, STRa(rg), 0, add_bytes);

        else
            seg_by_ctx (VB, STRa(rg), ctx, add_bytes);
    }

}

SPECIAL_RECONSTRUCTOR (sam_piz_special_RG_by_QNAME)
{
    if (reconstruct) {
        switch (snip[0] - '0') { // RG_method
            case RG_BY_ILLUM_QNAME: {
                STRlast (qname, SAM_QNAME);

                int64_t wi_plus_1;
                str_item_i_int (STRa(qname), ':', 3, &wi_plus_1); // note: we use str_item_i and not Q3NAME.last_value because QNAME might be segged by copy buddy, and different Illumina flavors have the RG in different items

                STR0(snip);

                ctx_get_snip_by_word_index (ctx, wi_plus_1-1, snip);
                RECONSTRUCT_snip;
                break;
            }

            default: ABORT_PIZ ("Invalid RG_method=%u", snip[0] - '0');
        }
    }

    return NO_NEW_VALUE; 
}

static inline unsigned sam_seg_aux_add_bytes (char type, unsigned value_len, bool is_bam)
{
    if (!is_bam || type=='Z' || type=='H')
        return value_len + 1; // +1 for \0 (BAM Z/H) or \t (SAM)

    else
        return aux_width[(uint8_t)type]; // BAM

    // note: this will return 0 for 'B'
}

static void sam_seg_initialize_for_float (VBlockSAMP vb, ContextP ctx)
{
    ctx->flags.store = STORE_FLOAT; // needs this be reconstructable as BAM
    ctx->ltype = IS_BAM_ZIP ? LT_FLOAT32 : LT_STRING;
    ctx->is_initialized = true;

    // case BAM: add all-the-same special (note: since we have local, we must seg one b250, its not enough to just add it to dict)
    if (IS_BAM_ZIP)
        seg_special0 (VB, SAM_SPECIAL_FLOAT, ctx, 0); 

    buf_alloc (VB, &ctx->local, 0, 16384, float, 0, CTX_TAG_LOCAL); // initial allocation, so we can use buf_append_one
} 

static inline SmallContainerP sam_seg_array_one_ctx_get_con (VBlockSAMP vb, ContextP con_ctx, uint8_t type, bool is_bam,
                                                             ContextP *elem_ctx) // out
{
    // case: cached with correct type
    if (con_ctx->con_cache.param == type) { // already initialized
        SmallContainerP con = B1ST (SmallContainer, con_ctx->con_cache);
        *elem_ctx = ctx_get_ctx (vb, con->items[1].dict_id);
        return con; // note: we return a pointer into con_cache- sam_seg_array_one_ctx may modify the number of repeats. that's ok.
    }

    buf_alloc (vb, &con_ctx->con_cache, 0, 1, SmallContainer, 0, "con_cache");

    SmallContainerP con = B1ST (SmallContainer, con_ctx->con_cache);
    *con = (SmallContainer){ .nitems_lo = 2, 
                             .drop_final_item_sep_of_final_repeat = true, // TODO - get rid of this flag and move to making the seperators to be repeat seperators as they should have been, using drop_final_repsep and obsoleting this flag 
                             .repsep    = {0,0}, 
                             .items     = { { .translator = SAM2BAM_ARRAY_SELF_1  },  // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM 
                                            { .separator  = {0, ','}            } } // item[1] is actual array item
                           };            
    
    // prepare context where array elements will go in
    con->items[1].dict_id      = DICT_ID_ARRAY (con_ctx->dict_id);
    
    con->items[1].translator   = aux_field_translator (type); // instructions on how to transform array items if reconstructing as BAM (array[0] is the subtype of the array)
    con->items[1].separator[0] = aux_sep_by_type[IS_BAM_ZIP][type];
    
    *elem_ctx = ctx_get_ctx (vb, con->items[1].dict_id);
    
    StoreType store_type = aux_field_store_flag[type];
    ASSERT (store_type, "%s: Invalid type \"%c\" in array of %s", LN_NAME, type, con_ctx->tag_name);

    (*elem_ctx)->flags.store = store_type;

    ctx_consolidate_stats_(VB, con_ctx, (ContainerP)con);

    if (store_type == STORE_INT || is_bam) {
        (*elem_ctx)->ltype = aux_field_to_ltype[type];
        (*elem_ctx)->local_is_lten = true; // we store in local in LTEN (as in BAM) and *not* in machine endianity
    }

    if (store_type == STORE_FLOAT)
        sam_seg_initialize_for_float (vb, *elem_ctx);

    con_ctx->con_cache.param = type; // this also used as "is_initialized"

    return con;
}

// an array - all elements go into a single item context, multiple repeats. items are segged as dynamic integers or floats, or a callback is called to seg them.
void sam_seg_array_one_ctx (VBlockSAMP vb, ZipDataLineSAMP dl, DictId dict_id, uint8_t type, 
                            rom array, int/*signed*/ array_len, // SAM: comma separated array ; BAM : arrays original width and machine endianity
                            ArrayItemCallback callback, void *cb_param, // optional - call back with array
                            PizSpecialReconstructor length_predictor)   // optional - SPECIAL function for predicting repeats
{   
    bool is_bam = IS_BAM_ZIP;

    uint32_t repeats = (is_bam || !array_len) ? array_len : (1 + str_count_char (STRa(array), ','));

    // prepare array container - a single item, with number of repeats of array element. array type is stored as a prefix
    ContextP con_ctx = ctx_get_ctx (vb, dict_id), elem_ctx;
    SmallContainerP con = sam_seg_array_one_ctx_get_con (vb, con_ctx, type, is_bam, &elem_ctx);
    ASSERTNOTNULL (elem_ctx);

    int width = aux_width[type];

    int array_bytes = is_bam ? (width * array_len) : array_len;
    elem_ctx->txt_len += array_bytes;

    // edge case note: if repeats==SEQ.len==0, we don't use CON_REPEATS_IS_SEQ_LEN because we have special handling of empty arrays in sam_piz_sam2bam_ARRAY_SELF_1
    con->repeats = ((dl && repeats == dl->SEQ.len && repeats) ? CON_REPEATS_IS_SEQ_LEN 
                 : length_predictor && repeats == ({ ValueType predicted={}; length_predictor (VB, con_ctx, 0, 0, &predicted, false); predicted.i; }) ? CON_REPEATS_IS_SPECIAL
                 : repeats);

    ASSERT (repeats < CONTAINER_MAX_REPEATS, "%s: array has too many elements, more than %u", LN_NAME, CONTAINER_MAX_REPEATS);

    bool is_int = (elem_ctx->flags.store == STORE_INT);
    if (is_int || is_bam) {
        elem_ctx->local.len *= width; // len will be calculated in bytes in this function
        buf_alloc (vb, &elem_ctx->local, repeats * width, repeats * width * 256, char, CTX_GROWTH, CTX_TAG_LOCAL); // careful not * line.len - we can get OOM
    }
     
    ASSERT (is_int || (type=='f' && !callback), "%s: Type not supported for SAM/BAM arrays '%c'(%u) in ctx=%s",
            LN_NAME, type, type, con_ctx->tag_name);

    char *local_start = BAFTc (elem_ctx->local);

    int64_t sum = 0;

    if (is_bam) {
        buf_add (&elem_ctx->local, array, array_bytes); // LTEN, not machine endianity (local_is_lten is set)
        
        if (con_ctx->flags.store == STORE_INT) {
            switch (type) {
                #define SUM_TYPE(type) for (type *n = (type *)local_start; n < (type *)(local_start + array_bytes); n++) sum += *n 
                case 'c':  SUM_TYPE(int8_t);   break; 
                case 'C':  SUM_TYPE(uint8_t);  break; 
                case 's':  SUM_TYPE(int16_t);  break; 
                case 'S':  SUM_TYPE(uint16_t); break; 
                case 'i':  SUM_TYPE(int32_t);  break; 
                case 'I':  SUM_TYPE(uint32_t); break;                 
                default : ABOSEG ("not implemented for type '%c'", type);
            }
            ctx_set_last_value (VB, con_ctx, sum);
        }
    }
    else { // SAM
        // note: we're not using str_split on array, because the number of elements can be very large (eg one per base in PacBio ip:B) - possibly stack overflow
        for (uint32_t i=0; i < repeats; i++) { // str_len will be -1 after last number

            rom snip = array;
            for (; array_len && *array != ','; array++, array_len--) {};

            unsigned snip_len = (unsigned)(array - snip);

            if (is_int) { 
                int64_t value;
                ASSERT (str_get_int (STRa(snip), &value), "%s: Invalid array: \"%.*s\", expecting an integer in an array element of %s", 
                        LN_NAME, STRf(snip), con_ctx->tag_name);
                value = LTEN64 (value); // consistent with BAM and avoid a condition before the memcpy below 
                memcpy (&local_start[i*width], &value, width);

                if (con_ctx->flags.store == STORE_INT) 
                    sum += value;
            }
            else // float
                seg_add_to_local_string (VB, elem_ctx, STRa(snip), LOOKUP_NONE, 0);
            
            array_len--; // skip comma
            array++;
        }

        if (is_int) elem_ctx->local.len += repeats * width;

        if (con_ctx->flags.store == STORE_INT)
            ctx_set_last_value (VB, con_ctx, sum);
    }

    if (con_ctx->flags.store == STORE_INDEX)
        ctx_set_last_value (VB, con_ctx, (int64_t)repeats);

    if (callback)
        callback (vb, elem_ctx, cb_param, local_start, repeats);

    if (is_int || is_bam)
       elem_ctx->local.len /= width; // return len back to counting in units of ltype
 
    if (repeats) {
        unsigned container_add_bytes = is_bam ? (4/*count*/ + 1/*type*/) : (2/*type - eg "i,"*/ + 1/*\t or \n*/);
        container_seg (vb, con_ctx, (ContainerP)con, ((char[]){ CON_PX_SEP, type, ',', CON_PX_SEP }), 4, container_add_bytes);
    }
    else {
        unsigned container_add_bytes = is_bam ? (4/*count*/ + 1/*type*/) : (1/*type - eg "i"*/ + 1/*\t or \n*/);
        container_seg (vb, con_ctx, (ContainerP)con, ((char[]){ CON_PX_SEP, type, CON_PX_SEP }), 3, container_add_bytes);
    }
}


static inline SmallContainerP sam_seg_array_multi_ctx_get_con (VBlockSAMP vb, ContextP con_ctx, uint32_t n_items, uint8_t type, bool is_bam)
{
    // case: cached with correct type
    if (con_ctx->con_cache.param == -(int8_t)type) { // already initialized
        SmallContainerP con = B1ST (SmallContainer, con_ctx->con_cache);

        if (con->nitems_lo == 1 + n_items) return con; // this is the container we need
    }

    ASSERT (n_items+1 <= SMALL_CON_NITEMS, "n_times=%u is too manyt", n_items);

    buf_alloc (vb, &con_ctx->con_cache, 0, 1, SmallContainer, 0, "con_cache");

    SmallContainerP con = B1ST (SmallContainer, con_ctx->con_cache);
    *con = (SmallContainer){ .nitems_lo           = 1 + n_items, 
                             .repeats             = 1,
                             .drop_final_item_sep = true,
                             .items[0]            = { .translator = SAM2BAM_ARRAY_SELF_M  } }; // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM 

    StoreType store_type = aux_field_store_flag[type];
    ASSERT (store_type, "%s: Invalid type \"%c\" in array of %s", LN_NAME, type, con_ctx->tag_name);

    for (uint32_t i=0; i < n_items; i++) {
        con->items[i+1] = (ContainerItem){ .dict_id    = sub_dict_id (con_ctx->dict_id, '0'+i),
                                           .separator  = { [0]=aux_sep_by_type[IS_BAM_ZIP][type], [1]=',' },
                                           .translator = aux_field_translator (type) }; // instructions on how to transform array items if reconstructing as BAM (array[0] is the subtype of the array)
        
        ContextP item_ctx = ctx_get_ctx (vb, con->items[i+1].dict_id);
        item_ctx->flags.store = store_type;

        if (store_type == STORE_INT || is_bam) {
            item_ctx->ltype = aux_field_to_ltype[type];
            item_ctx->local_is_lten = true; // we store in local in LTEN (as in BAM) and *not* in machine endianity
        }

        if (store_type == STORE_FLOAT)
            sam_seg_initialize_for_float (vb, item_ctx);
    }

    ctx_consolidate_stats_(VB, con_ctx, (ContainerP)con);

    con_ctx->con_cache.param = -(int8_t)type; // this also used as "is_initialized". negative to indicate multi-context (positive is single-context)

    return con;
}

// an array - each element go into a its own context, multiple repeats. items are segged as dynamic integers or floats, or a callback is called to seg them.
void sam_seg_array_multi_ctx (VBlockSAMP vb, ZipDataLineSAMP dl, ContextP con_ctx, uint8_t type, uint32_t expected_n_subctxs,
                              rom array, int/*signed*/ array_len) // SAM: comma separated array ; BAM : arrays original width and machine endianity
{   
    ASSERT (expected_n_subctxs >= 2, "expecting n_subctxs=%u >= 2", expected_n_subctxs);

    bool is_bam = IS_BAM_ZIP;

    uint32_t n_subctxs = (is_bam || !array_len) ? array_len : (1 + str_count_char (STRa(array), ','));
    
    // if array does not have the expected length, seg as single-context array 
    if (expected_n_subctxs != n_subctxs) {
        sam_seg_array_one_ctx (vb, dl, con_ctx->dict_id, type, STRa(array), NULL, 0, NULL);
        return;
    }
    
    // prepare array container - a one context per element. array type is stored as a prefix
    SmallContainerP con = sam_seg_array_multi_ctx_get_con (vb, con_ctx, n_subctxs, type, is_bam);

    int width = aux_width[type];

    if (is_bam)
        for (uint32_t i=0; i < n_subctxs; i++) {
            ContextP item_ctx = ctx_get_ctx (vb, con->items[i+1].dict_id);
            
            buf_insert_do (VB, &item_ctx->local, width, item_ctx->local.len, &array[i * width], 1, CTX_TAG_LOCAL, __FUNCLINE);
            
            item_ctx->txt_len += width;
        }

    else if (aux_field_store_flag[type] == STORE_INT) { // SAM - integers
        str_split_ints (array, array_len, n_subctxs, ',', item, true);

        for (uint32_t i=0; i < n_subctxs; i++) {
            ContextP item_ctx = ctx_get_ctx (vb, con->items[i+1].dict_id);
            
            // note: &items[i] contains the data because in little endian LSB comes first. TO DO: verify that value fits in width
            buf_insert_do (VB, &item_ctx->local, width, item_ctx->local.len, &items[i], 1, CTX_TAG_LOCAL, __FUNCLINE);
        }

        con_ctx->txt_len += array_len; // not ideal, but saves us the need to measure lengths of the items
    }

    else { // SAM - floats
        str_split (array, array_len, n_subctxs, ',', item, true);

        for (uint32_t i=0; i < n_subctxs; i++) {
            ContextP item_ctx = ctx_get_ctx (vb, con->items[i+1].dict_id);
            seg_add_to_local_string (VB, item_ctx, STRi(item, i), LOOKUP_NONE, item_lens[i]);
        }

        con_ctx->txt_len += n_subctxs-1; // commas
    }

    unsigned container_add_bytes = is_bam ? (4/*count*/ + 1/*type*/) : (2/*type - eg "i,"*/ + 1/*\t or \n*/);
    container_seg (vb, con_ctx, (ContainerP)con, ((char[]){ CON_PX_SEP, type, ',', CON_PX_SEP }), 4, container_add_bytes);
}

static void sam_seg_set_last_value_f_from_aux (VBlockSAMP vb, Did did_i, bool is_bam,
                                               STRp(value),       // used if SAM 
                                               ValueType numeric) // .f32 needs to be set if BAM
{
    ContextP ctx = CTX(did_i);

    if (is_bam) 
        ctx_set_last_value (VB, ctx, (ValueType){ .f = numeric.f32.f }); // note: just casting f32.f to double works in gcc but not clang

    else {
        SAFE_NULT(value);
        char *after;
        numeric.f = strtod (value, &after); // same as in reconstruct_one_snip
        SAFE_RESTORE;

        ASSSEG (after == value + value_len, "%s: Expecting float value for auxiliary field %s but found \"%.*s\"",
                LN_NAME, ctx->tag_name, STRf (value));

        ctx_set_last_value (VB, ctx, numeric);
    }
}

void sam_seg_float_as_snip (VBlockSAMP vb, ContextP ctx, STRp(sam_value), ValueType bam_value, unsigned add_bytes)
{
    if (IS_BAM_ZIP) {
        // in case of BAM, numeric contains the float value - we seg its binary representation
        SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_FLOAT, bam_value.i);
        seg_by_dict_id (VB, STRa(snip), ctx->dict_id, add_bytes); 
    }
    
    else 
        // in case of SAM, we seg the string so we maintain the precise textual representation. 
        seg_by_ctx (VB, STRa(sam_value), ctx, add_bytes);

    if (ctx->flags.store == STORE_FLOAT)
        sam_seg_set_last_value_f_from_aux (vb, ctx->did_i, IS_BAM_ZIP, STRa(sam_value), bam_value);
}

// note: also called from fastq_seg_one_saux()
void sam_seg_aux_field_fallback (VBlockP vb, void *dl, DictId dict_id, char sam_type, char array_subtype,
                                 STRp(value), ValueType numeric, unsigned add_bytes)
{
    // all types of integer
    if (sam_type == 'i') {
        ContextP ctx = ctx_get_ctx (vb, dict_id);
        seg_integer (vb, ctx, numeric.i, true, add_bytes);
    }

    else if (sam_type == 'f') {
        ContextP ctx = ctx_get_ctx (vb, dict_id);

        if (!ctx->is_initialized) 
            sam_seg_initialize_for_float (VB_SAM, ctx);

        // note: for some fields, sam_seg_float_as_snip() is better.
        if (IS_BAM_ZIP) {
            uint32_t f = numeric.i; // bam_get_one_aux stores float values as binary-identical 32bit integer
            seg_integer_fixed (VB, ctx, &f, LOOKUP_NONE, add_bytes);
        }
        
        else // sam
            seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, add_bytes);
    }

    // Numeric array
    else if (sam_type == 'B') 
        sam_seg_array_one_ctx (VB_SAM, (ZipDataLineSAMP )dl, dict_id, array_subtype, STRa(value), NULL, NULL, NULL);

    // Z,H,A - normal snips in their own dictionary
    else 
        seg_by_dict_id (VB, STRa(value), dict_id, add_bytes);     
}

// process an optional subfield, that looks something like MX:Z:abcdefg. We use "MX" for the field name, and
// the data is abcdefg. The full name "MX:Z:" is stored as part of the AUX dictionary entry
DictId sam_seg_aux_field (VBlockSAMP vb, ZipDataLineSAMP dl, bool is_bam, 
                          rom tag, char bam_type, char array_subtype, 
                          STRp(value), ValueType numeric, // two options 
                          int16_t idx) // 0-based index of this field within AUX
{
    char sam_type = sam_seg_bam_type_to_sam_type (bam_type);
    char dict_name[6] = { tag[0], tag[1], ':', sam_type, ':', array_subtype }; // last 2 are ignored if not array
    DictId dict_id = dict_id_make (dict_name, (sam_type=='B' ? 6 : 4), DTYPE_SAM_AUX); // match dict_id as declared in #pragma GENDICT

    unsigned add_bytes = sam_seg_aux_add_bytes (bam_type, value_len, is_bam);

    #define COND0(condition, seg) if (condition) { seg; break; } else  
    #define COND(condition,  seg) if (condition) { seg; break; } else goto fallback

    if (segconf_running)
        segconf.has[ctx_get_ctx (VB, dict_id)->did_i]++;

    switch (dict_id.num) {

        // ---------------------
        // Standard fields
        // ---------------------
        case _OPTION_SA_Z: sam_seg_SA_Z (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_OA_Z: sam_seg_OA_Z (vb, STRa(value), add_bytes); break;

        case _OPTION_OC_Z: sam_seg_other_CIGAR (vb, CTX(OPTION_OC_Z), STRa(value), true, add_bytes); break;

        case _OPTION_OP_i: seg_pos_field (VB, OPTION_OP_i, SAM_POS, 0, 0, 0, 0, numeric.i, add_bytes); break;

        case _OPTION_OQ_Z: sam_seg_other_qual (vb, dl, &dl->OQ, NULL, OPTION_OQ_Z, STRa(value), add_bytes); break;

        // case _OPTION_CC_Z: // see bug 609

        case _OPTION_MC_Z: sam_cigar_seg_MC_Z (vb, dl, STRa(value), add_bytes); break;

        // note: we cannot seg MD using SPECIAL if we're segging the line against SA Group, because
        // the MD alg requires the SQBITMAP.
        case _OPTION_MD_Z: sam_seg_MD_Z (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_NM_i: COND (!MP(BLASR), sam_seg_NM_i (vb, dl, (SamNMType)numeric.i, add_bytes)); // blasr uses NM:i in a non-standard way

        case _OPTION_BQ_Z: sam_seg_BQ (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_BD_Z:
        case _OPTION_BI_Z: sam_seg_BD_BI_Z (vb, dl, STRa(value), dict_id, add_bytes); break;
        
        case _OPTION_RG_Z: sam_seg_RG_Z (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_PG_Z: sam_seg_buddied_Z_fields (vb, dl, MATED_PG, STRa(value), 0, add_bytes); break;
        case _OPTION_PU_Z: sam_seg_buddied_Z_fields (vb, dl, MATED_PU, STRa(value), 0, add_bytes); break;
        case _OPTION_LB_Z: sam_seg_buddied_Z_fields (vb, dl, MATED_LB, STRa(value), 0, add_bytes); break;
        case _OPTION_OX_Z: sam_seg_buddied_Z_fields (vb, dl, MATED_OX, STRa(value), 0, add_bytes); break;

        case _OPTION_YB_Z: COND (segconf.sam_has_xcons, sam_seg_buddied_Z_fields (vb, dl, MATED_YB, STRa(value), 0, add_bytes));
        
        case _OPTION_MI_Z: COND0 (!MP(ULTIMA), sam_seg_buddied_Z_fields (vb, dl, MATED_MI, STRa(value), 0, add_bytes))
                           COND (  MP(ULTIMA), sam_seg_ultima_MI (vb, dl, STRa(value), add_bytes));
        case _OPTION_CY_Z: dl->dont_compress_CY = sam_seg_barcode_qual (vb, dl, OPTION_CY_Z, SOLO_CY, segconf.n_CR_CB_CY_seps, STRa(value), qSTRa(segconf.CY_con_snip), &segconf.CY_con, add_bytes); break;
        case _OPTION_QT_Z: dl->dont_compress_QT = sam_seg_barcode_qual (vb, dl, OPTION_QT_Z, SOLO_QT, segconf.n_BC_QT_seps,    STRa(value), qSTRa(segconf.QT_con_snip), &segconf.QT_con, add_bytes); break;

        case _OPTION_CR_Z: sam_seg_CR_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_CB_Z: sam_seg_CB_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_UR_Z: // aliases of RX (cellranger)
        case _OPTION_RX_Z: sam_seg_RX_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_UB_Z: // aliases of BX (cellranger)
        case _OPTION_BX_Z: sam_seg_BX_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_UY_Z: // aliases of QX (cellranger)
        case _OPTION_QX_Z: sam_seg_QX_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_BC_Z: sam_seg_BC_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_MM_Z: sam_seg_MM_Z (vb, STRa(value), add_bytes); break;

        case _OPTION_GP_i: COND (MP(CRDNA), sam_seg_GP_i (vb, dl, numeric.i, add_bytes));
        case _OPTION_MP_i: COND (MP(CRDNA), sam_seg_MP_i (vb, dl, numeric.i, add_bytes));

        case _OPTION_fx_Z: COND (segconf.has_10xGen, sam_seg_fx_Z (vb, dl, STRa(value), add_bytes));
        case _OPTION_xf_i: COND (segconf.has_10xGen, sam_seg_xf_i (vb, dl, numeric.i, add_bytes));
        case _OPTION_mm_i: COND (segconf.has_10xGen, seg_integer_as_snip_do (VB, CTX(OPTION_mm_i), numeric.i, add_bytes));
        case _OPTION_MM_i: COND (segconf.has_10xGen, seg_integer_as_snip_do (VB, CTX(OPTION_MM_i), numeric.i, add_bytes));
        case _OPTION_pa_i: COND (segconf.has_10xGen, sam_seg_SEQ_END (vb, dl, CTX(OPTION_pa_i), numeric.i, "00RF", add_bytes));
        case _OPTION_ts_i: COND (segconf.has_10xGen, sam_seg_SEQ_END (vb, dl, CTX(OPTION_ts_i), numeric.i, "00FR", add_bytes));
        case _OPTION_2Y_Z: COND (segconf.has_10xGen, sam_seg_other_qual (vb, dl, NULL, &dl->_2Y, OPTION_2Y_Z, STRa(value), add_bytes));        
        case _OPTION_GR_Z: COND (segconf.has_10xGen, sam_seg_GR_Z (vb, dl, STRa(value), add_bytes));
        case _OPTION_GY_Z: COND (segconf.has_10xGen, sam_seg_GY_Z (vb, dl, STRa(value), add_bytes));
        case _OPTION_2R_Z: COND (segconf.has_10xGen, sam_seg_other_seq (vb, dl, OPTION_2R_Z, STRa(value), add_bytes));
        case _OPTION_TR_Z: COND (segconf.has_TR_TQ,  sam_seg_other_seq (vb, dl, OPTION_TR_Z, STRa(value), add_bytes));
        case _OPTION_TQ_Z: COND (segconf.has_TR_TQ,  sam_seg_other_qual (vb, dl, NULL, &dl->TQ, OPTION_TQ_Z, STRa(value), add_bytes));
        case _OPTION_TX_Z: COND (segconf.has_10xGen, sam_seg_TX_AN_Z (vb, dl, OPTION_TX_Z, STRa(value), add_bytes));
        case _OPTION_AN_Z: COND (segconf.has_10xGen, sam_seg_TX_AN_Z (vb, dl, OPTION_AN_Z, STRa(value), add_bytes));
        case _OPTION_CQ_Z: COND (segconf.has_10xGen, sam_seg_other_qual (vb, dl, NULL, &dl->CQ, OPTION_CQ_Z, STRa(value), add_bytes));

        case _OPTION_GN_Z: sam_seg_gene_name_id (vb, dl, OPTION_GN_Z, STRa(value), add_bytes); break;
        case _OPTION_GX_Z: sam_seg_gene_name_id (vb, dl, OPTION_GX_Z, STRa(value), add_bytes); break;
        case _OPTION_gn_Z: sam_seg_gene_name_id (vb, dl, OPTION_gn_Z, STRa(value), add_bytes); break;
        case _OPTION_gx_Z: sam_seg_gene_name_id (vb, dl, OPTION_gx_Z, STRa(value), add_bytes); break;

        
        //case _OPTION_E2: sam_seg_E2_field (vb, dl, STRa(value), add_bytes); // BROKEN. To do: fix.

        case _OPTION_U2_Z: sam_seg_U2_Z (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_AS_i: sam_seg_AS_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_MQ_i: sam_seg_MQ_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_PQ_i: sam_seg_PQ_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_SM_i: sam_seg_SM_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_AM_i: sam_seg_AM_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_HI_i: sam_seg_HI_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_NH_i: sam_seg_NH_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_CP_i: sam_seg_CP_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_UQ_i: sam_seg_UQ_i (vb, dl, numeric.i, add_bytes); break;

        // ---------------------
        // Non-standard fields
        // ---------------------
        case _OPTION_XA_Z: COND (segconf.sam_has_BWA_XA_Z != no, sam_seg_BWA_XA_Z (vb, STRa(value), add_bytes)); // yes or unknown
        
        case _OPTION_XS_i: COND0 (segconf.sam_has_BWA_XS_i, sam_seg_BWA_XS_i (vb, dl, OPTION_XS_i, numeric.i, add_bytes))
                           COND (MP(BLASR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_XS_i), numeric.i, "0+00", add_bytes)); 

        case _OPTION_XC_i: COND0 (segconf.sam_has_BWA_XC_i, sam_seg_BWA_XC_i (vb, dl, numeric.i, add_bytes))
                           COND (segconf.sam_has_xcons, sam_seg_xcons_XC_i (vb, numeric.i, add_bytes));

        case _OPTION_YY_i: COND (segconf.sam_has_xcons, sam_seg_xcons_YY_i (vb, numeric.i, add_bytes));
        
        case _OPTION_XO_i: COND (segconf.sam_has_xcons, sam_seg_xcons_XO_i (vb, dl, numeric.i, add_bytes));

        case _OPTION_XT_A: COND (segconf.sam_has_BWA_XT_A, sam_seg_BWA_XT_A (vb, value[0], add_bytes));

        case _OPTION_nM_i: COND (MP(STAR), sam_seg_STAR_nM (vb, dl, (SamNMType)numeric.i, add_bytes)); 
        
        case _OPTION_jI_B_i: COND (MP(STAR), sam_seg_STAR_jI (vb, dl, STRa(value), is_bam)); 

        case _OPTION_jM_B_c: COND (MP(STAR), sam_seg_array_one_ctx (VB_SAM, dl, dict_id, array_subtype, STRa(value), NULL, NULL, sam_piz_special_jM_length));

        case _OPTION_XM_i: COND0 (segconf.sam_has_XM_i_is_mismatches, sam_seg_XM_i (vb, dl, numeric.i, idx, add_bytes))
                           COND (MP(TMAP), sam_seg_tmap_XM_i (vb, numeric.i, add_bytes));

        case _OPTION_XT_i: COND (MP(TMAP), sam_seg_tmap_XT_i (vb, dl, numeric.i, add_bytes)); 
        
        case _OPTION_X1_i: COND (segconf.sam_has_BWA_X01_i, sam_seg_BWA_X1_i (vb, numeric.i, add_bytes)); 

        case _OPTION_ZS_i: COND (MP(HISAT2), sam_seg_BWA_XS_i (vb, dl, OPTION_ZS_i, numeric.i, add_bytes)); 

        case _OPTION_YS_i: COND (segconf.sam_has_bowtie2_YS_i, sam_seg_bowtie2_YS_i (vb, dl, numeric, add_bytes));

        // case _OPTION_XR_Z: COND (segconf.sam_has_bismark_XM_XG_XR, sam_seg_bismark_XR_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_XO_Z: COND (MP(BSSEEKER2), sam_seg_bsseeker2_XO_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_XG_Z: COND0 (segconf.sam_has_bismark_XM_XG_XR, sam_seg_bismark_XG_Z (vb, dl, STRa(value), add_bytes))
                           COND (MP(BSSEEKER2), sam_seg_bsseeker2_XG_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_XM_Z: COND0 (segconf.sam_has_bismark_XM_XG_XR, sam_seg_bismark_XM_Z (vb, dl, OPTION_XM_Z, SAM_SPECIAL_BISMARK_XM, STRa(value), add_bytes))
                           COND (MP(BSSEEKER2), sam_seg_bsseeker2_XM_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_XB_A: COND (MP(GEM3), sam_seg_gem3_XB_A (vb, dl, STRa(value), add_bytes));

        case _OPTION_XB_Z: COND (MP(BSBOLT), sam_seg_bsbolt_XB (vb, dl, STRa(value), add_bytes));

        case _OPTION_YS_Z: COND (MP(BSBOLT), sam_seg_bsbolt_YS_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_mc_i: COND (segconf.is_biobambam2_sort, sam_seg_mc_i (vb, numeric.i, add_bytes));

        case _OPTION_ms_i: COND0 (segconf.sam_ms_type == ms_BIOBAMBAM && !segconf.optimize[SAM_QUAL], sam_seg_ms_i (vb, dl, numeric.i, add_bytes)) // ms:i produced by biobambam or samtools
                           COND  (segconf.sam_ms_type == ms_MINIMAP2, sam_seg_delta_seq_len (vb, dl, OPTION_ms_i, numeric.i, true, add_bytes));    // ms:i produced by minimap2

        case _OPTION_s1_i: COND (segconf.is_minimap2, sam_seg_s1_i (vb, dl, numeric.i, add_bytes));

        case _OPTION_cm_i: COND (segconf.is_minimap2, sam_seg_cm_i (vb, dl, numeric.i, add_bytes));

        case _OPTION_Z5_i: seg_pos_field (VB, OPTION_Z5_i, SAM_PNEXT, 0, 0, 0, 0, numeric.i, add_bytes); break;

        case _OPTION_Zs_Z: COND (MP(HISAT2), sam_seg_HISAT2_Zs_Z (vb, STRa(value), add_bytes));

        case _OPTION_YH_Z: COND (MP(NOVOALIGN), seg_add_to_local_string (VB, CTX(OPTION_YH_Z), STRa(value), LOOKUP_NONE, add_bytes)); break;
        case _OPTION_YQ_Z: COND (MP(NOVOALIGN), seg_add_to_local_string (VB, CTX(OPTION_YQ_Z), STRa(value), LOOKUP_NONE, add_bytes)); break;

        case _OPTION_QS_i: COND (MP(NGMLR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_QS_i), numeric.i, "00+0", add_bytes));
        case _OPTION_QE_i: COND (MP(NGMLR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_QE_i), numeric.i, "+00-", add_bytes));
        case _OPTION_XR_i: COND (MP(NGMLR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_XR_i), numeric.i, "+0--", add_bytes));
        case _OPTION_XE_i: COND (MP(BLASR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_XE_i), numeric.i, "++00", add_bytes));
        case _OPTION_XQ_i: COND (MP(BLASR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_XQ_i), numeric.i, "+000", add_bytes));
        case _OPTION_XL_i: COND (MP(BLASR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_XL_i), numeric.i, "+0--", add_bytes));

        case _OPTION_FI_i: COND (MP(BLASR), sam_seg_blasr_FI_i (vb, dl, numeric.i, add_bytes));

        // PacBio fields
        case _OPTION_dq_Z: COND0 (segconf.use_pacbio_iqsqdq, sam_seg_pacbio_xq (vb, dl, OPTION_dq_Z, STRa(value), add_bytes))
                           COND (TECH(PACBIO), seg_add_to_local_string (VB, CTX(OPTION_dq_Z), STRa(value), LOOKUP_SIMPLE, add_bytes)); 
        
        case _OPTION_iq_Z: COND0 (segconf.use_pacbio_iqsqdq, sam_seg_pacbio_xq (vb, dl, OPTION_iq_Z, STRa(value), add_bytes))
                           COND (TECH(PACBIO), seg_add_to_local_string (VB, CTX(OPTION_iq_Z), STRa(value), LOOKUP_SIMPLE, add_bytes)); 
        
        case _OPTION_sq_Z: COND0 (segconf.use_pacbio_iqsqdq, sam_seg_pacbio_xq (vb, dl, OPTION_sq_Z, STRa(value), add_bytes))
                           COND (TECH(PACBIO), seg_add_to_local_string (VB, CTX(OPTION_sq_Z), STRa(value), LOOKUP_SIMPLE, add_bytes)); 

        case _OPTION_qs_i: COND (TECH(PACBIO), sam_seg_pacbio_qs (vb, dl, numeric.i, add_bytes));
        case _OPTION_qe_i: COND (TECH(PACBIO), sam_seg_pacbio_qe (vb, dl, numeric.i, add_bytes));
        case _OPTION_we_i: COND (TECH(PACBIO), sam_seg_pacbio_we (vb, dl, numeric.i, add_bytes));
        case _OPTION_np_i: COND (TECH(PACBIO), sam_seg_pacbio_np (vb, dl, numeric.i, add_bytes));

        case _OPTION_sn_B_f : COND (TECH(PACBIO), sam_seg_pacbio_sn (vb, dl, STRa(value)));

        case _OPTION_ec_f: COND (TECH(PACBIO) && segconf.has[OPTION_ec_f], ({ sam_seg_set_last_value_f_from_aux (vb, OPTION_ec_f, is_bam, STRa(value), numeric); goto fallback; })); // consumed by np:i

        case _OPTION_dt_Z: COND (TECH(PACBIO), seg_add_to_local_string (VB, CTX(OPTION_dt_Z), STRa(value), LOOKUP_SIMPLE, add_bytes));
        case _OPTION_mq_Z: COND (TECH(PACBIO), seg_add_to_local_string (VB, CTX(OPTION_mq_Z), STRa(value), LOOKUP_SIMPLE, add_bytes));
        case _OPTION_st_Z: COND (TECH(PACBIO), seg_add_to_local_string (VB, CTX(OPTION_st_Z), STRa(value), LOOKUP_SIMPLE, add_bytes));
        
        case _OPTION_zm_i: COND(segconf.sam_has_zm_by_Q1NAME, sam_seg_pacbio_zm (vb, numeric.i, add_bytes));

        // case _OPTION_ql_Z: _OPTION_qt_Z: // better off as snips that QUAL-like context

        // Illumina DRAGEN fields
        case _OPTION_sd_f: COND (MP(DRAGEN), sam_dragen_seg_sd_f (vb, dl, STRa(value), numeric, add_bytes));

        // Ultima Genomics fields
        case _OPTION_tp_B_c: COND (segconf.has[OPTION_tp_B_c] && !dl->no_qual && segconf.has[OPTION_tp_B_c], sam_seg_array_one_ctx (vb, dl, _OPTION_tp_B_c, array_subtype, STRa(value), sam_seg_ULTIMA_tp, dl, NULL));
        case _OPTION_bi_Z: COND (MP(ULTIMA), sam_seg_ultima_bi (vb, STRa(value), add_bytes));
        case _OPTION_a3_i: COND (MP(ULTIMA), sam_seg_delta_seq_len (vb, dl, OPTION_a3_i, numeric.i, false, add_bytes));
        case _OPTION_XV_Z: COND (MP(ULTIMA), sam_seg_ultima_XV (vb, STRa(value), add_bytes));
        case _OPTION_XW_Z: COND (MP(ULTIMA), sam_seg_ultima_XW (vb, STRa(value), add_bytes));
        case _OPTION_rq_f: COND (MP(ULTIMA) || segconf.pacbio_subreads, sam_seg_float_as_snip (vb, CTX(OPTION_rq_f), STRa(value), numeric, add_bytes)); // expecting all-the-same for subreads

        case _OPTION_t0_Z: COND (segconf.sam_has_ultima_t0, sam_seg_ultima_t0 (vb, dl, STRa(value), add_bytes));
        
        // cpu
        case _OPTION_Y0_i: COND (MP(CPU) && sam_seg_peek_int_field (vb, OPTION_AS_i, vb->idx_AS_i, 0, 10000, true, NULL), seg_delta_vs_other_dictN (VB, CTX(OPTION_Y0_i), CTX(OPTION_AS_i), numeric.i, 10, add_bytes));
        case _OPTION_Y1_i: COND (MP(CPU) && sam_seg_peek_int_field (vb, OPTION_XS_i, vb->idx_XS_i, 0, 10000, true, NULL), seg_delta_vs_other_localN (VB, CTX(OPTION_Y1_i), CTX(OPTION_XS_i), numeric.i, 1000, add_bytes));
        case _OPTION_XL_Z: COND (MP(CPU), sam_seg_CPU_XL_Z (vb, STRa(value), add_bytes));
        
        // Agilent
        case _OPTION_ZA_Z: COND (segconf.has_agent_trimmer, sam_seg_buddied_Z_fields (vb, dl, MATED_ZA, STRa(value), 0, add_bytes));
        case _OPTION_ZB_Z: COND (segconf.has_agent_trimmer, sam_seg_buddied_Z_fields (vb, dl, MATED_ZB, STRa(value), 0, add_bytes));

        // NanoSeq
        case _OPTION_rb_Z: COND (segconf.sam_is_nanoseq, sam_seg_cross_mated_Z_fields (vb, OPTION_rb_Z, dl, STRa(value), &dl->rb, &dl->mb, &vb->mux_rb, STRa(copy_mate_mb_snip), add_bytes));
        case _OPTION_mb_Z: COND (segconf.sam_is_nanoseq, sam_seg_cross_mated_Z_fields (vb, OPTION_mb_Z, dl, STRa(value), &dl->mb, &dl->rb, &vb->mux_mb, STRa(copy_mate_rb_snip), add_bytes));
        
        // Abra2
        case _OPTION_YA_Z: COND (segconf.sam_has_abra2, sam_seg_ABRA2_YA_Z (vb, STRa(value), add_bytes));
        case _OPTION_YO_Z: COND (segconf.sam_has_abra2, sam_seg_ABRA2_YO_Z (vb, STRa(value), add_bytes));

        default: fallback:
            sam_seg_aux_field_fallback (VB, dl, dict_id, sam_type, array_subtype, STRa(value), numeric, add_bytes);
    }
    
    return dict_id;
}
