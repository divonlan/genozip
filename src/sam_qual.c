// ------------------------------------------------------------------
//   sam_qual.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "piz.h"
#include "reconstruct.h"
#include "optimize.h"
#include "codec.h"
#include "htscodecs/rANS_static4x16.h"

rom bam_qual_display (bytes qual, uint32_t l_seq) // caller should free memory
{
    bool valid_qual = true;
    for (uint32_t i=0; i < l_seq; i++) 
        if (qual[i] > 93) {
            valid_qual = false;
            break;
        }

    if (valid_qual) {
        char *str = MALLOC (l_seq + 2);

        for (uint32_t i=0; i < l_seq; i++) 
            str[i] = qual[i] + 33;

        str[l_seq] = 0;
        return str;
    }

    else {
        char *str = MALLOC (l_seq*3 + 2);
        return str_to_hex (qual, l_seq, str, true);
    }
}

// main thread (not thread-safe): called from sam_show_sag_one_grp for getting first few characters of alignment cigar
rom sam_display_qual_from_SA_Group (const Sag *g)
{
    if (g->no_qual) return "*";

    static char qual[SA_QUAL_DISPLAY_LEN+1];
    memset (qual, 0, sizeof(qual));

    uint32_t uncomp_len = MIN_(SA_CIGAR_DISPLAY_LEN, (uint32_t)g->seq_len); // possibly shorter than original cigar

    if (g->qual_comp_len) { // qual of this group is compressed
        uint32_t uncomp_len = MIN_((uint32_t)g->seq_len, SA_QUAL_DISPLAY_LEN);
        void *success = rans_uncompress_to_4x16 (evb, B8(z_file->sag_qual, g->qual), g->qual_comp_len,
                                                (uint8_t *)qual, &uncomp_len); 
        if (success && uncomp_len) qual[uncomp_len] = '\0';
    }

    else // not compressed
        memcpy (qual, B8(z_file->sag_qual, g->qual), uncomp_len);

    return qual;
}

//---------
// QUAL SEG
//---------

// callback function for compress to get data of one line
COMPRESSOR_CALLBACK (sam_zip_qual) 
{
    ZipDataLineSAMP dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = (dl->dont_compress_QUAL || dl->is_consensus) ? 0 : MIN_(maximum_size, dl->QUAL.len);

    if (__builtin_expect (!line_data, false)) return; // only lengths were requested

    *line_data = Btxt (dl->QUAL.index);

    // if QUAL is just "*" (i.e. unavailable) replace it by " " because '*' is a legal PHRED quality value that will confuse PIZ
    if (__builtin_expect (dl->QUAL.len == 1 && (*line_data)[0] == '*', false))
        *line_data = " "; // pointer to static string

    if (is_rev) *is_rev = dl->FLAG.rev_comp;
}

COMPRESSOR_CALLBACK (sam_zip_cqual) 
{
    ZipDataLineSAMP dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = (dl->dont_compress_QUAL || !dl->is_consensus) ? 0 : MIN_(maximum_size, dl->QUAL.len);

    if (!line_data) return; // only lengths were requested

    *line_data = Btxt (dl->QUAL.index);

    if (is_rev) *is_rev = dl->FLAG.rev_comp;
}

#define QUAL_ZIP_CALLBACK(tag, f, may_be_revcomped)             \
COMPRESSOR_CALLBACK (sam_zip_##tag)                             \
{                                                               \
    ZipDataLineSAMP dl = DATA_LINE (vb_line_i);                 \
    *line_data_len = dl->dont_compress_##tag ? 0 : MIN_(maximum_size, dl->f.len); /* note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec */ \
    if (!line_data || ! *line_data_len) return; /* no dat, or only lengths were requested */   \
    *line_data = Btxt (dl->f.index);                             \
    if (is_rev) *is_rev = may_be_revcomped ? dl->FLAG.rev_comp : false;\
}                               

QUAL_ZIP_CALLBACK(OQ, OQ, true)
QUAL_ZIP_CALLBACK(TQ, TQ, true)
QUAL_ZIP_CALLBACK(GY, GY, false)
QUAL_ZIP_CALLBACK(2Y, _2Y, true)
QUAL_ZIP_CALLBACK(QX, solo_z_fields[SOLO_QX], false)
QUAL_ZIP_CALLBACK(CY, solo_z_fields[SOLO_CY], false)
QUAL_ZIP_CALLBACK(QT, solo_z_fields[SOLO_QT], false)

void sam_seg_QUAL_initialize (VBlockSAMP vb)
{
    if (IS_MAIN(vb)) 
        ctx_set_store_per_line (VB, SAM_QUAL, SAM_FLAG, DID_EOL);

    // case: in DEPN and MAIN QUALSA is used to store the diff vs primary 
    if (IS_DEPN(vb) || IS_MAIN(vb))
        CTX(SAM_QUALSA)->ltype = LT_INT8; // diff of prim vs depn qual

    // case: PRIM - TOPLEVEL reconstructs all-the-same SAM_QUALSA instead of SAM_QUAL. SAM_QUAL is consumed when loading SA Groups.
    if (IS_PRIM(vb)) 
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_QUAL, '0' }, 3, SAM_QUALSA, 0); // note: we can't just ctx_create_node because we need to transfer flags to piz, so need b250
}

// ZIP/PIZ: decompresses grp qual of grp, into vb->scratch
static void sam_get_sa_grp_qual (VBlockSAMP vb)
{
    ASSERTNOTINUSE (vb->scratch);

    buf_alloc (vb, &vb->scratch, vb->sag->seq_len, 0, char, 1, "scratch");

    if (vb->sag->qual_comp_len) { // qual of this group is compressed
        uint32_t uncomp_len = vb->sag->seq_len;
        void *success = rans_uncompress_to_4x16 (VB, B8(z_file->sag_qual, vb->sag->qual), vb->sag->qual_comp_len,
                                                 B1ST8(vb->scratch), &uncomp_len); 

        ASSGOTO (success && uncomp_len == vb->sag->seq_len, "%s: rans_uncompress_to_4x16 failed to decompress an SA Group QUAL data: grp_i=%u success=%u comp_len=%u uncomp_len=%u expected_uncomp_len=%u qual=%"PRIu64,
                    LN_NAME, ZGRP_I(vb->sag), !!success, vb->sag->qual_comp_len, uncomp_len, vb->sag->seq_len, (uint64_t)vb->sag->qual);

        vb->scratch.len = uncomp_len;
    }
    else // not compressed
        buf_add (&vb->scratch, Bc(z_file->sag_qual, vb->sag->qual), vb->sag->seq_len);

    return;
    
error:
    sam_show_sag_one_grp (ZGRP_I(vb->sag));
    exit_on_error(true);
}

static bool sam_seg_QUAL_diff_do (VBlockSAMP vb, rom my_qual, rom other_qual, uint32_t qual_len, bool reverse_other, 
                                  int8_t *diff)  // out
{
    uint32_t num_diff = 0; // number of base scores for which diff=0
    memset (diff, 0, qual_len);

    if (!reverse_other) {
        for (uint32_t i=0; i < qual_len; i++) 
            if (my_qual[i] != other_qual[i]) {
                diff[i] = my_qual[i] - other_qual[i]; // a value [-93,93] as qual is [33,126]
                num_diff++;
            }
    }
    else {
        // note: in case of reverse_other, other_qual should point to the LAST character
        for (int32_t/*signed*/ i=0; i < qual_len; i++) 
            if (my_qual[i] != other_qual[-i]) {
                diff[i] = my_qual[i] - other_qual[-i]; 
                num_diff++;
            }
    }

    return num_diff < qual_len/2; // low enough number of changes to make it worth diffing
}

// diff current QUAL against prim or saggy's QUAL. return false if diff was aborted
typedef enum { QDT_DEFAULT=0, QDT_ABORTED=1, QDT_REVERSED=2/*v15*/ } QualDiffType; // this values are part of the file format
static QualDiffType sam_seg_QUAL_diff (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(my_qual), 
                               STRp (other_qual), bool other_revcomp, uint32_t *other_hard_clip, unsigned add_bytes)
{
    QualDiffType diff_type = QDT_DEFAULT; // optimistic
    bool xstrand = (dl->FLAG.rev_comp != other_revcomp);

    // flanks are the regions of my_qual that are not covered by other_qual, due to other_qual's hard clips.
    // example1: prim: HHHHQQQQQHH  (H=Hard clip, Q=some qual score)
    //           qual: HHQQQQQHHHH
    //                   ^^ flank[0]==2  flank[1]==0 (bc qual Q's don't go beyond prim's) overlap_len=3 (bc 3 Qs overlapping)
    //
    // example2: prim: HHHHHHHHQQH  (H=Hard clip, Q=some qual score)
    //           qual: HHQQQHHHHHH flank={3,0} overlap_len=0

    uint32_t flank[2] = { MIN_(my_qual_len, MAX_((int32_t)other_hard_clip[xstrand]  - (int32_t)vb->hard_clip[0], 0)),   // left-flanking
                          MIN_(my_qual_len, MAX_((int32_t)other_hard_clip[!xstrand] - (int32_t)vb->hard_clip[1], 0)) }; // right-flanking

    uint32_t overlap_len = my_qual_len - flank[0] - flank[1];

    buf_alloc (vb, &CTX(SAM_QUALSA)->local, overlap_len, 0, int8_t, CTX_GROWTH, CTX_TAG_LOCAL);
    int8_t *diff = BAFT (int8_t, CTX(SAM_QUALSA)->local);

    rom overlap_my_qual = &my_qual[flank[0]];

    // update other_qual if partial
    if (!xstrand) {
        uint32_t first_other_overlap = flank[0] ? 0 : (vb->hard_clip[0] - other_hard_clip[0]);
        other_qual = &other_qual[first_other_overlap]; // FIRST character of other qual
    }
    else {
        uint32_t last_other_overlap = other_qual_len - 1 - (flank[0] ? 0 : (vb->hard_clip[0] - other_hard_clip[1]));
        other_qual = &other_qual[last_other_overlap]; // LAST character of other qual
    }

    // abort diff if we have too many different base scores - we're better off not diffing
    if (!sam_seg_QUAL_diff_do (vb, overlap_my_qual, other_qual, overlap_len, xstrand, diff)) {  
        // rescue QUAL that suffer from an NGMLR bug - where QUAL strings of supplementary alignments sometimes have the wrong orientation
        if (segconf.sam_mapper == MP_NGMLR && 
            !vb->hard_clip[0] && !vb->hard_clip[1] && !other_hard_clip[0] && !other_hard_clip[1]) {

            if (xstrand) other_qual -= overlap_len-1; // previously other_qual was the last character, change it to the first character
            else         other_qual += overlap_len-1; // previously other_qual was the first character, change it to the last character

            if (!sam_seg_QUAL_diff_do (vb, overlap_my_qual, other_qual, overlap_len, !xstrand, diff))
                return QDT_ABORTED;
            else  
                diff_type = QDT_REVERSED;
        }
        else 
            return QDT_ABORTED;
    }

    // add flanks (regions that could not be diffed) to SAM_QUAL_FLANK.local
    if (flank[0]) 
        seg_add_to_local_fixed (VB, CTX(SAM_QUAL_FLANK), my_qual, flank[0], LOOKUP_NONE, 0);

    if (flank[1])
        seg_add_to_local_fixed (VB, CTX(SAM_QUAL_FLANK), &my_qual[my_qual_len - flank[1]], flank[1], LOOKUP_NONE, 0);

    CTX(SAM_QUALSA)->local.len32 += overlap_len;
    CTX(SAM_QUAL)->txt_len += add_bytes - flank[0] - flank[1];
    CTX(SAM_QUAL_FLANK)->txt_len += flank[0] + flank[1];

    return diff_type; // all good
}

static void sam_seg_QUAL_segconf (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(qual)/*always textual*/, bool monochar)
{
    if (!monochar && !vb->qual_missing) segconf.nontrivial_qual = true;

    for (uint32_t i=0; i < qual_len; i++)
        if (IS_QUAL_SCORE(qual[i]))
            segconf.qual_histo[dl->is_consensus ? QHT_CONSENSUS : QHT_QUAL][qual[i]-33].count++;
}

void sam_seg_QUAL (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(qual)/*always textual*/, unsigned add_bytes)
{
    START_TIMER;

    ContextP qual_ctx  = CTX(SAM_QUAL);
    ContextP cqual_ctx = CTX(SAM_CQUAL);
    ZipDataLineSAMP saggy_dl;
    bool prim_has_qual_but_i_dont = false; // will be set if this line has no QUAL, but its prim line does (very rare)
    QualDiffType diff_type = QDT_DEFAULT;  // quality scores are too different, we're better off not diffing (added 14.0.10)
    bool pacbio_diff = false;
    char monochar = 0;
    
    if (!vb->qual_missing) vb->has_qual = true;
    
    // case --optimize_QUAL: optimize in place, must be done before sam_deep_set_QUAL_hash calculates a hash
    if (flag.optimize_QUAL && !vb->qual_missing) 
        optimize_phred_quality_string ((char*)STRa(qual));

    if ((flag.deep || flag.show_deep == 2) && !segconf.running)
        sam_deep_set_QUAL_hash (vb, dl, STRa(qual));

    // case: monochar 
    if (!vb->qual_missing && IS_MAIN(vb) && str_is_monochar (STRa(qual))) { // note: not for prim vb - complicates loading + no need for depn - getting from prim
        monochar = qual[0];
        qual_ctx->txt_len += add_bytes;
        dl->dont_compress_QUAL = true; // don't compress this line
    }
        
    // note: if prim (of either type) has no QUAL, we don't attempt to diff - bc piz wouldn't be able to know whether 
    // the current line has QUAL or not

    // case: DEPN component, line has a sag, and the sag has qual
    else if (vb->sag && IS_DEPN(vb) && !vb->sag->no_qual) {
        if (vb->qual_missing) {
            prim_has_qual_but_i_dont = true;
            CTX(SAM_QUALSA)->txt_len += add_bytes;
        }

        // case: both this line and prim have QUAL - diff them
        else {
            sam_get_sa_grp_qual (vb); // decompress prim qual into vb->scratch

            diff_type = sam_seg_QUAL_diff (vb, dl, STRa(qual), STRb(vb->scratch), vb->sag->revcomp, (uint32_t[]){0, 0}, add_bytes);        
            buf_free (vb->scratch);

            if (diff_type == QDT_ABORTED) goto standard;
        }

        dl->dont_compress_QUAL = true; // don't compress this line
    }

    // case: MAIN component (note: vb->saggy_line_i is only set in MAIN), the line has a corresponding 
    // saggy line in same VB, and saggy line has qual
    else if (sam_has_saggy && // note: saggy_line_i is set by sam_seg_saggy only for lines in MAIN component
             ({ saggy_dl = DATA_LINE (vb->saggy_line_i); 
                !saggy_dl->no_qual && 
                // note: SEQ.len is set if ANY of QUAL, SEQ, CIGAR-implied-seq_consumed have it
                dl->hard_clip[0] + dl->hard_clip[1] + dl->SEQ.len == saggy_dl->hard_clip[0] + saggy_dl->hard_clip[1] + saggy_dl->SEQ.len; })) {

        if (vb->qual_missing) {
            prim_has_qual_but_i_dont = true;
            qual_ctx->txt_len += add_bytes;
        }    
        else {
            diff_type = sam_seg_QUAL_diff (vb, dl, STRa(qual), STRtxt(saggy_dl->QUAL), saggy_dl->FLAG.rev_comp, saggy_dl->hard_clip, add_bytes);
            if (diff_type == QDT_ABORTED) goto standard;
        }

        dl->dont_compress_QUAL = true; // don't compress this line
    }

    // 1. if we suspect entire file might be qual-less, seg as snip (without adding local) so that would become 
    // all-the-same in that. If we might have mixed qual/no-qual in the file, we add to local to not add entropy to b250
    // 2. LONGR codec can't handle missing QUAL
    else if (!IS_PRIM(vb) && dl->no_qual && (!segconf.nontrivial_qual || codec_longr_maybe_used (VB, SAM_QUAL))) {
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_QUAL, '*' }, 3, SAM_QUAL, add_bytes); 
        dl->dont_compress_QUAL = true; // don't compress this line
        goto done;
    }

    // case: predict QUAL from dq, iq, sq
    else if (!dl->no_qual && segconf.use_pacbio_iqsqdq &&
             sam_seg_pacbio_qual (vb, STRa (qual), add_bytes)) {
        
        pacbio_diff = true;
        dl->dont_compress_QUAL = true; // don't compress this line with QUAL codec
    }

    // case: standard
    // Note: in PRIM, QUAL is not reconstructed (as QUAL is not in TOPLEVEL container) - it is consumed when loading SA Groups
    //       Instead, all-the-same QUALSA is reconstructed (SPECIAL copying from the SA Group)
    else standard: {
        ContextP ctx = dl->is_consensus ? cqual_ctx : qual_ctx;
        ctx->local.len32 += dl->QUAL.len;
        ctx->txt_len     += add_bytes;
    }   

    // seg SPECIAL. note: in prim this is all-the-same and segged in sam_seg_QUAL_initialize
    seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_QUAL, '0' + prim_has_qual_but_i_dont, '0' + diff_type, '0' + pacbio_diff, monochar }, 5 + !!monochar, 
                IS_PRIM(vb) ? SAM_QUALSA : SAM_QUAL, 0); 
   
    // get QUAL score, consumed by mate ms:i
    if (!segconf.running && segconf.sam_ms_type == ms_BIOBAMBAM && !flag.optimize_QUAL)
        dl->QUAL_score = sam_get_QUAL_score (vb, STRa(qual));
 
    if (segconf.running) 
        sam_seg_QUAL_segconf (vb, dl, STRa(qual), monochar);

done:
    COPY_TIMER (sam_seg_QUAL);
}

//-----------------------------------------------------------------------------------------------------
// OQ:Z - Original QUAL (standard)
//        "Original base quality, usually before recalibration. Same encoding as QUAL" - https://samtools.github.io/hts-specs/SAMtags.pdf
// GY:Z - CellRanger
// 2Y:Z - CellRanger - R2 qual?
// TQ:Z - CellRanger & longranger: Quality values of the 7 trimmed bases following the barcode sequence at the start of R1. Can be used to reconstruct the original R1 quality values.
//-----------------------------------------------------------------------------------------------------
void sam_seg_other_qual (VBlockSAMP vb, ZipDataLineSAMP dl, TxtWord *dl_word, Did did_i, STRp(qual), bool len_is_seq_len, unsigned add_bytes)
{
    decl_ctx (did_i);

    *dl_word = TXTWORD(qual);
    ctx->local.len32 += qual_len;
    ctx->txt_len     += add_bytes;

    if (!len_is_seq_len) 
        seg_lookup_with_length (VB, ctx, qual_len, 0);
    else {
        ASSSEG (qual_len == dl->SEQ.len, "Expecting %s to have length == SEQ.len=%u but its length is %u. %s=\"%.*s\"",
                ctx->tag_name, dl->SEQ.len, qual_len, ctx->tag_name, STRf(qual));
    }

    if (segconf.running && did_i == OPTION_OQ_Z) 
        for (uint32_t i=0; i < qual_len; i++)
            if (IS_QUAL_SCORE(qual[i]))
                segconf.qual_histo[QHT_OQ][qual[i]-33].count++;
}

void sam_update_qual_len (VBlockP vb, uint32_t line_i, uint32_t new_len) 
{ 
    DATA_LINE(line_i)->QUAL.len = new_len; 
}

//---------
// QUAL PIZ
//---------

// undiff vs primary or saggy
static void sam_piz_QUAL_undiff_vs_other (VBlockSAMP vb, STRp (other_qual), QualDiffType diff_type, bool other_revcomp, const CigarAnalItem *saggy_anal, bool other_is_bam, ReconType reconstruct)
{
    bool xstrand = (last_flags.rev_comp != other_revcomp);

    if (diff_type == QDT_REVERSED) xstrand = !xstrand; // QDT_REVERSED introduced in v15

    int32_t qual_len = other_qual_len + (saggy_anal ? (saggy_anal->hard_clip[0] + saggy_anal->hard_clip[1]) : 0)
                     - (vb->hard_clip[0] + vb->hard_clip[1]);

    ASSERT (qual_len >= 0, "qual_len=%d: other_qual_len=%u saggy_anal->hard_clip=[%u,%u] vb->hard_clip=[%u,%u]", 
            qual_len, other_qual_len, saggy_anal ? saggy_anal->hard_clip[0] : 0, saggy_anal ? saggy_anal->hard_clip[1] : 0, vb->hard_clip[0], vb->hard_clip[1]);

    uint32_t flank[2] = { saggy_anal ? MIN_(qual_len, MAX_((int32_t)saggy_anal->hard_clip[xstrand]  - (int32_t)vb->hard_clip[0], 0)) : 0,   // left-flanking
                          saggy_anal ? MIN_(qual_len, MAX_((int32_t)saggy_anal->hard_clip[!xstrand] - (int32_t)vb->hard_clip[1], 0)) : 0 }; // right-flanking

    int32_t overlap_len = qual_len - flank[0] - flank[1];

    ASSPIZ (overlap_len >= 0, "Expecting qual_len=%u >= flank=[%u,%u] xstrand=%u saggy_anal->hard_clip=[%u,%u] vb->hard_clip=[%u,%u] CIGAR=\"%.*s\" vb->saggy_line_i=%u is_depn_vb=%s has_saggy=%s",
            qual_len, flank[0], flank[1], xstrand, saggy_anal ? saggy_anal->hard_clip[0] : 0, saggy_anal ? saggy_anal->hard_clip[1] : 0, 
            vb->hard_clip[0], vb->hard_clip[1], STRfb(vb->textual_cigar), vb->saggy_line_i, TF(IS_DEPN(vb)), TF(sam_has_saggy));

    rom diff = 0;
    if (overlap_len) {
        ContextP qualsa_ctx = LOADED_CTX(SAM_QUALSA);
        diff = Bc (qualsa_ctx->local, qualsa_ctx->next_local);
        qualsa_ctx->next_local += overlap_len;

        ASSPIZ (qualsa_ctx->next_local <= qualsa_ctx->local.len32, "Out of diff-vs-primary data in QUALSA.local: overlap_len=%u qual_len=%u flank=[%u,%u] QUALSA.local.len=%u is_depn_vb=%s has_saggy=%s", 
                overlap_len, qual_len, flank[0], flank[1], qualsa_ctx->local.len32, TF(IS_DEPN(vb)), TF(sam_has_saggy)); 
    }

    char *qual = BAFTtxt;
    int8_t bam_bump = (other_is_bam ? 33 : 0);

    uint32_t save_seq_len = vb->seq_len;

    if (flank[0]) {
        vb->seq_len = flank[0];
        reconstruct_from_local_sequence (VB, CTX(SAM_QUAL_FLANK), vb->seq_len, reconstruct); // reconstruct vb->seq_len quality scores
        qual += flank[0];
    }

    if (overlap_len && !xstrand && reconstruct) {
        uint32_t first_other_overlap = flank[0] ? 0 : (vb->hard_clip[0] - (saggy_anal ? saggy_anal->hard_clip[0] : 0));
        ASSPIZ (first_other_overlap + overlap_len <= other_qual_len, "other_qual overflow: expecting first_other_overlap=%u + overlap_len=%u <= other_qual_len=%u. flank[0]=%u vb->hard_clip[0]=%u saggy_anal->seq_len=%u saggy_anal->hard_clip[0]=%u is_depn_vb=%s has_saggy=%s bam_bump=%d",
                first_other_overlap, overlap_len, other_qual_len, flank[0], vb->hard_clip[0], saggy_anal->seq_len, saggy_anal->hard_clip[0], TF(IS_DEPN(vb)), TF(sam_has_saggy), bam_bump);

        other_qual += first_other_overlap;
        
        for (uint32_t i=0; i < overlap_len; i++) 
            qual[i] = other_qual[i] + diff[i] + bam_bump;

        Ltxt += overlap_len;
    }

    else if (overlap_len && xstrand && reconstruct) {
        uint32_t last_other_overlap = other_qual_len - 1 - (flank[0] ? 0 : (vb->hard_clip[0] - (saggy_anal ? saggy_anal->hard_clip[1] : 0)));

        ASSPIZ (last_other_overlap >= overlap_len-1, "other_qual underflow: expecting last_other_overlap=%u >= overlap_len=%u-1. flank[0]=%u vb->hard_clip[0]=%u saggy_anal->seq_len=%d saggy_anal->hard_clip[1]=%d is_depn_vb=%s has_saggy=%s bam_bump=%d",
                last_other_overlap, overlap_len, flank[0], vb->hard_clip[0], saggy_anal ? saggy_anal->seq_len : -1, saggy_anal ? saggy_anal->hard_clip[1] : -1, TF(IS_DEPN(vb)), TF(sam_has_saggy), bam_bump);

        other_qual += last_other_overlap;

        for (int32_t/*signed*/ i=0; i < overlap_len; i++) 
            qual[i] = other_qual[-i] + diff[i] + bam_bump;

        Ltxt += overlap_len;
    }

    if (flank[1]) {
        vb->seq_len = flank[1];
        reconstruct_from_local_sequence (VB, CTX(SAM_QUAL_FLANK), vb->seq_len, reconstruct);
    }

    vb->seq_len = save_seq_len;
}

static void sam_piz_QUAL_primary (VBlockSAMP vb)
{
    sam_get_sa_grp_qual (vb); // uncompress PRIM qual to vb->scratch

    RECONSTRUCT_BUF (vb->scratch);
    buf_free (vb->scratch);
}

void sam_reconstruct_missing_quality (VBlockP vb, ReconType reconstruct)
{
    if (reconstruct) 
        RECONSTRUCT1 ('*');

    VB_SAM->qual_missing = true;

    if (!vb->preprocessing) 
        B(CigarAnalItem, CTX(SAM_CIGAR)->cigar_anal_history, vb->line_i)->qual_missing = true;
}

// Note: in PRIM, it is called with ctx=QUALSA, in MAIN and DEPN with ctx=QUAL
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_QUAL)
{
    START_TIMER;

    VBlockSAMP vb = (VBlockSAMP)vb_;
    char *qual = BAFTtxt;
    const Sag *g = vb->sag;
    bool prim_has_qual_but_i_dont = (snip[0] == '1');
    QualDiffType diff_type        = (snip[1] - '0');
    bool pacbio_diff              = (snip_len >= 3 && snip[2] == '1'); // v15
    char monochar                 = (snip_len >= 4 ? snip[3] : 0);     // v15
    char is_missing               = (snip[0] == '*');

    const CigarAnalItem *saggy_anal;

    if (monochar) {
        char *c = BAFTtxt;
        char *aft = c + vb->seq_len;
        while (c < aft) *c++ = monochar;
        Ltxt += vb->seq_len;
    }
    
    // case: diff against prediction based on PacBio dq, iq, sq fields
    else if (pacbio_diff)
        sam_recon_pacbio_qual (vb, ctx, reconstruct);

    // case: reconstruct by copying from sag (except if we are depn and group has no qual)
    else if (SAM_PIZ_HAS_SAG && (IS_PRIM(vb) || (IS_DEPN(vb) && !g->no_qual)) && !is_missing) {
      
        if (!reconstruct) {}

        else if (g->no_qual || prim_has_qual_but_i_dont) // note: g->no_qual is set if the line has no QUAL, and also if QUAL is not loaded due to flags
            sam_reconstruct_missing_quality (VB, reconstruct);

        else if (IS_DEPN(vb)) {
            if (diff_type == QDT_ABORTED) goto no_diff;

            sam_get_sa_grp_qual (vb); // uncompress PRIM qual to vb->scratch
            
            sam_piz_QUAL_undiff_vs_other (vb, STRb(vb->scratch), diff_type, vb->sag->revcomp, NULL, false, reconstruct);                
            buf_free (vb->scratch);
        }
        else  // primary vb
            sam_piz_QUAL_primary (vb);
    }

    // case: MAIN component, reconstruct depn line against prim line in this VB
    else if (sam_has_saggy && // only appears in the MAIN component and only since v14 (see sam_seg_saggy)
             ({ saggy_anal = B(CigarAnalItem, CTX(SAM_CIGAR)->cigar_anal_history, vb->saggy_line_i);
                !saggy_anal->qual_missing &&
                // note: we've seen cases in the wild where a depn without hard clips is shorter than its prim (eg in test.NA12878.chr22.1x.bam), possibly due to GATK IndelRealigner
                vb->hard_clip[0] + vb->hard_clip[1] + vb->seq_len == saggy_anal->hard_clip[0] + saggy_anal->hard_clip[1] + saggy_anal->seq_len; })) {

        HistoryWord word = *B(HistoryWord, ctx->history, vb->saggy_line_i); // QUAL is always stored as LookupTxtData or LookupPerLine
        SamFlags saggy_flags = { .value = history64 (SAM_FLAG, vb->saggy_line_i) };
        rom saggy_qual = (word.lookup == LookupTxtData) ? Btxt (word.index) : Bc(ctx->per_line, word.index);

        if (word.len == 1 && (uint8_t)*saggy_qual == (OUT_DT(BAM) ? 0xff : '*'))
            goto no_diff; // in case prim has no qual - seg did not diff (the current line may or may not have qual)

        else if (prim_has_qual_but_i_dont)
            sam_reconstruct_missing_quality (VB, reconstruct);

        else if (diff_type == QDT_ABORTED) 
            goto no_diff;

        else
            sam_piz_QUAL_undiff_vs_other (vb, saggy_qual, word.len, diff_type, saggy_flags.rev_comp, saggy_anal, OUT_DT(BAM), reconstruct);                
    } 
        
    // case: reconstruct from data in local
    else no_diff: 
        if (is_missing)  // case of seg suspecting entire file is qual-less
            sam_reconstruct_missing_quality (VB, reconstruct);

        else {
            bool is_consensus = (!IS_PRIM(vb)/*bug 949*/ && ctx_encountered_in_line (VB, SAM_QNAME2) && segconf.flav_prop[QNAME2].is_consensus);
            if (is_consensus) ctx = CTX(SAM_CQUAL);
            
            switch (ctx->ltype) { // the relevant subset of ltypes from reconstruct_from_ctx_do
                case LT_CODEC:
                    codec_args[ctx->lcodec].reconstruct (VB, ctx->lcodec, ctx, vb->seq_len, reconstruct); break;

                case LT_BLOB: 
                    reconstruct_from_local_sequence (VB, ctx, vb->seq_len, reconstruct); break;

                default: ABORT_PIZ ("Invalid ltype=%s for %s. local.len=%u", lt_name (ctx->ltype), ctx->tag_name, ctx->local.len32);
            }
        }

    uint32_t qual_len = BAFTtxt - qual;

    // if we're using pacbio diff in this VB, but skipped it this line (eg due to monochar), update next_local
    if (CTX(OPTION_iq_sq_dq)->is_loaded && !pacbio_diff)
        sam_recon_skip_pacbio_qual (vb);
        
#ifdef DEBUG
    for (uint32_t i=0; i < qual_len; i++)
        ASSPIZ (IS_QUAL_SCORE(qual[i]), "Invalid QUAL character reconstructed: i=%u char='%c'(%u)\n", 
                i, qual[i], (uint8_t)qual[i]);
#endif

    // store quality score in ms:i history 
    if (!VER(14)) // up to v13 this special was only when calculating score was needed
        new_value->i = sam_get_QUAL_score (vb, STRa(qual));
    
    // calculate QUAL_score to be copied into our future mate's ms:i
    else if (CTX(OPTION_ms_i)->flags.spl_custom) // since v14 for ms_BIOBAMBAM
        *B(int64_t, CTX(OPTION_ms_i)->history, vb->line_i) = sam_get_QUAL_score (vb, STRa(qual));

    COPY_TIMER (sam_piz_special_QUAL);

    return !VER(14); // has new value for files up to v13 (QUAL history is used for TxtWord since v14)
}

// SAM-to-BAM translator: translate SAM ASCII (33-based) Phred values to BAM's 0-based
TRANSLATOR_FUNC (sam_piz_sam2bam_QUAL)
{
    START_TIMER;

    // before translating - add to Deep if needed
    if (flag.deep) 
        sam_piz_con_item_cb (vb, &(ContainerItem){ .dict_id = ctx->dict_id }, STRa(recon));

    // if QUAL is "*" there are two options:
    // 1. If l_seq is 0, the QUAL is empty
    // 2. If not (i.e. we have SEQ data but not QUAL) - it is a string of 0xff, length l_seq
    if (VB_SAM->qual_missing) {
        BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)Btxt (vb->line_start);
        uint32_t l_seq = GET_UINT32_(alignment, l_seq);
        
        if (!l_seq) // option 1
            Ltxt--;
        
        else if (VB_SAM->qual_missing && !segconf.pysam_qual) {  // option 2 - missing QUAL according to SAM spec
            memset (BLSTtxt, 0xff, l_seq); // override the '*' and l_seq-1 more
            Ltxt += l_seq - 1;
        }

        else {  // option 3 - missing QUAL as created by pysam (non-compliant with SAM spec)
            *BLSTtxt = 0xff;
            memset (BAFTtxt, 0, l_seq-1); // filler is 0 instead of 0xff as in SAM SPEC
            Ltxt += l_seq - 1;
        }
    }
    
    else // we have QUAL - update Phred values
        for (uint32_t i=0; i < recon_len; i++)
            recon[i] -= 33; 

    COPY_TIMER (sam_piz_sam2bam_QUAL);
    return 0;
}

