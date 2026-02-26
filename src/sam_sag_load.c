// ------------------------------------------------------------------
//   sam_sag_load.c
//   Copyright (C) 2022-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "sam_private.h"
#include "compressor.h"
#include "dispatcher.h"
#include "zfile.h"
#include "qname.h"
#include "writer.h"
#include "huffman.h"
#include "htscodecs/arith_dynamic.h"

typedef struct { // one per PRIM VB
    VBIType vblock_i;
    SAGroup first_grp_i; 
    uint32_t num_grps, num_alns;
    uint64_t first_aln_i;
    uint64_t comp_seq_start;    // in bases
    uint32_t comp_seq_len;           // in bases
    uint32_t comp_qual_len;     // size as of compressed in-memory in Sag: huffman or arith compress as in ingest for files since 15.0.68, and arith-compressed for older files (length not know a-priori as compressed ingest using another codec) 
    uint64_t comp_qual_start;
    uint64_t comp_solo_data_start; // used by SAG_BY_SOLO 
    uint32_t comp_solo_data_len;// used by SAG_BY_SOLO: length of data for this Plsg in z_file->solo_data (each solo field is huffman-compressed)
    uint32_t comp_cigars_len;   // used by SAG_BY_SA: length of data for this Plsg in z_file->sag_cigars (each cigar is huffman-compressed for files since 15.0.68 and uncompressed for older files)
    uint64_t comp_cigars_start; // used by SAG_BY_SA
    uint64_t comp_qnames_start;
    uint32_t comp_qnames_len;   // length of data for this Plsg in z_file->sag_qnames (each qname is huffman-compressed for files since 15.0.65 and uncompressed for older files)
} PlsgVbInfo;

#define plsg_info z_file->vb_info[1]
static uint32_t next_plsg_i = 0; // iterator for dispatching compute threads
static VBIType num_prim_vbs_loaded = 0; // number of PRIM VBs whose loading is complete

static Mutex copy_qnames_mutex={}, copy_qual_mutex={}, copy_cigars_mutex={}, copy_solo_data_mutex={};

// buffers for in-memory compressed data in case total length consumed by this data is not
// known in advance because we changed the compression alogirhtm since zip (i.e. not PRECISE)
#define vb_cigars_buf    vb->codec_bufs[0]
#define vb_qnames_buf    vb->codec_bufs[1]
#define vb_solo_data_buf vb->codec_bufs[2]
#define vb_qual_buf      vb->z_data
static Flags save_flag; 

static bool is_precise_qual, is_precise_cigars, is_precise_qnames, is_precise_solo;

//-------------------------------------------------------------------------
// PIZ: SA Group loader
// Method: decompress directly into memory, parallelizing with dispatcher and 3 threads per PRIM VB: SEQ, QUAL and Alignments
//-------------------------------------------------------------------------

ShowAln sam_show_sag_one_SA_aln (VBlockP vb, const Sag *g, const SAAln *a)
{
    ShowAln s;
    
    if (IS_PIZ) {
        char cigar_info[64];
        if (a->cigar.piz.is_word) 
            snprintf (cigar_info, sizeof (cigar_info), "word=%u", (WordIndex)a->cigar.piz.index);
        else
            snprintf (cigar_info, sizeof (cigar_info), "in=%"PRIu64" l=%u", (uint64_t)a->cigar.piz.index, a->cigar.piz.len);

        char rname_str[64]; // should be big enough...
        snprintf (rname_str, sizeof (rname_str), "\"%.48s\"(%u)", ctx_get_snip_by_word_index0 (ZCTX(VER(15) ? SAM_RNAME/*alias since v15*/ : OPTION_SA_RNAME), a->rname), a->rname);

        snprintf (s.s, sizeof (s.s), "aln_i=%u.%u: sa_rname=%-13s pos=%-10u mapq=%-3u strand=%c nm=%-3u cigar[%u]=\"%s\"(%s)",
                  ZGRP_I(g), (unsigned)(ZALN_I(a) - g->first_aln_i), rname_str,
                  a->pos, a->mapq, "+-"[a->revcomp], a->nm,
                  SA_CIGAR_DISPLAY_LEN, sam_piz_display_aln_cigar(vb, a).s, cigar_info);
    }
    else 
        snprintf (s.s, sizeof (s.s), "aln_i=%u.%u: cigar_sig=%016"PRIx64" rname=%-4u pos=%-10u mapq=%-3u strand=%c  nm=%-3u",
                ZGRP_I(g), (unsigned)(ZALN_I(a) - g->first_aln_i), a->cigar.signature, a->rname, a->pos, a->mapq, "+-"[a->revcomp], a->nm);
    
    return s;
}

void sam_show_sag_one_grp (VBlockP vb, SAGroup grp_i)
{
    const Sag *g = B(Sag, z_file->sag_grps, grp_i);

    StrTextMegaLong extra = {};
    unsigned extra_len=0;

    if (IS_SAG_SA)
        SNPRINTF (extra, " first_aln=%"PRIu64, (uint64_t)g->first_aln_i);

    else if (IS_SAG_CC) {
        CCAln *aln = B(CCAln, z_file->sag_alns, grp_i);
        SNPRINTF (extra, " rname=%s pos=%u", 
                  (command==ZIP) ? ctx_get_z_snip_ex (ZCTX(SAM_RNAME), aln->rname, 0, 0) 
                                 : ctx_get_snip_by_word_index0 (ZCTX(SAM_RNAME), aln->rname), 
                  aln->pos);
    }

    else if (IS_SAG_SOLO) {
        SoloAlnP aln = B(SoloAln, z_file->sag_alns, grp_i);
        SNPRINTF (extra, " aln_index=%"PRIu64, U40to64(aln->index));
        for (SoloTags solo=0; solo < NUM_SOLO_TAGS; solo++) 
            if (aln->field_uncomp_len[solo]) {
                char uncomp[4 KB];
                huffman_uncompress (solo_props[solo].did_i, sam_solo_sag_data (VB_SAM, aln, solo), uncomp, aln->field_uncomp_len[solo]);
                SNPRINTF (extra, " %s=\"%.*s\"", ZCTX(solo_props[solo].did_i)->tag_name, aln->field_uncomp_len[solo], uncomp);
            }
    }

    char qname[g->qname_len];
    huffman_uncompress (SAM_QNAME, GRP_QNAME(g), qname, g->qname_len);

    iprintf ("grp_i=%u: qname(i=%"PRIu64",l=%u,hash=%u)=\"%.*s\" seq=(i=%"PRIu64",l=%u) qual=(i=%"PRId64",comp_l=%u)[%u]=\"%s\" AS=%u strand=%c mul/fst/lst=%u,%u,%u %s=%u%s%s\n",
             grp_i, (uint64_t)g->qname, g->qname_len, 
             (uint32_t)qname_calc_hash (QNAME1, COMP_NONE, qname, g->qname_len, g->is_last, false, CRC32, NULL), // the hash is part of the index in ZIP, and not part of the SAGroup struct. Its provided here for its usefulness in debugging
             g->qname_len, qname, (uint64_t)g->seq, g->seq_len, 
             g->qual == SA_QUAL_ERROR ? -1 : (int64_t)g->qual, g->qual_comp_len, 
             SA_QUAL_DISPLAY_LEN, g->qual == SA_QUAL_ERROR ? "Error": sam_display_qual_from_sag (vb, g).s,  
             (int)g->as, "+-"[g->revcomp], g->multi_segs, g->is_first, g->is_last,
             ((IS_SAG_NH || IS_SAG_CC || IS_SAG_SOLO) ? "NH" : "num_alns"), g->num_alns, 
             extra.s, g->first_grp_in_vb ? " FIRST_GRP_IN_VB" : "");

    if (IS_SAG_SA)
        for (uint64_t aln_i=g->first_aln_i; aln_i < g->first_aln_i + g->num_alns; aln_i++) 
            iprintf ("    %s\n", sam_show_sag_one_SA_aln (vb, g, B(SAAln, z_file->sag_alns, aln_i)).s);
}

void sam_show_sag (void)
{
    if (flag.show_sag >= 1) { // show a specific group
        ASSINP (flag.show_sag <= z_file->sag_grps.len, "--show-sag argument, for this file, should be between 0 and %"PRIu64, 
                z_file->sag_grps.len-1);
        
        sam_show_sag_one_grp (evb, flag.show_sag - 1);
    }

    else {
        iprintf ("SAG_TYPE=%s NUM_SA_GROUPS=%"PRIu64" NUM_ALIGNMENTS=%"PRIu64"\n", 
                 sag_type_name (segconf.sag_type), z_file->sag_grps.len, z_file->sag_alns.count);

        for (SAGroup grp_i=0; grp_i < z_file->sag_grps.len; grp_i++) 
            sam_show_sag_one_grp (evb, grp_i);
    }
}

static inline void sam_load_groups_add_qnames (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps, Sag *g,
                                               bytes start, bytes after, uint8_t **next) 
{
    START_TIMER;
    decl_ctx (SAM_QNAME);
    
    if (!ctx->is_loaded) return;

    ContextP seq_len_ctx = segconf.seq_len_dict_id.num ? ECTX(segconf.seq_len_dict_id) : NULL;
        
    rom uncomp_qname = BAFTtxt;
    ctx->mate_copied_exactly = false; 
    reconstruct_from_ctx (vb, SAM_QNAME, 0, RECON_ON);  // reconstructs into vb->txt_data, sets vb->buddy_line_i if SNIP_COPY_BUDDY
    uint32_t uncomp_qname_len = BAFTtxt - uncomp_qname;
    
    ASSPIZ (uncomp_qname_len <= SAM_MAX_QNAME_LEN, "unexpectedly, uncomp_qname_len=%u > SAM_MAX_QNAME_LEN=%u", uncomp_qname_len, SAM_MAX_QNAME_LEN);

    g->qname_len = uncomp_qname_len;

    if (ctx->mate_copied_exactly)
        g->qname = vb_grps[ctx->last_value.i].qname; // last_value is buddy_line_i     

    else if (start) { // precise - compress directly into z_file
        g->qname = *next - start;
    
        uint32_t comp_len = after - *next; // theoretical maximum 
        *next += huffman_compress (VB, SAM_QNAME, STRa(uncomp_qname), *next, &comp_len); // using huffman sent via SEC_HUFFMAN

        ASSERT (after >= *next, "PRIM/%u: plsg->comp_qnames_len=%u exceeded: (after-next)=%"PRId64, 
                vb->plsg_i, plsg->comp_qnames_len, (int64_t)(after - *next));
    }

    else {
        uint32_t comp_len = huffman_get_theoretical_max_comp_len (SAM_QNAME, uncomp_qname_len);
        buf_alloc (vb, &vb_qnames_buf, comp_len, 0, uint8_t, 1, NULL); // likely already allocated in sam_load_groups_add_grps

        huffman_compress (VB, SAM_QNAME, STRa(uncomp_qname), BAFT8(vb_qnames_buf), &comp_len);

        g->qname = vb_qnames_buf.len; // relative to the VB. we will update later to be relative to z_file->sag_qnames.
        vb_qnames_buf.len += comp_len;
    }

    // if seq_len is carried by a QNAME item, set the last value here (needed for sam_cigar_special_CIGAR)
    if (seq_len_ctx) g->seq_len = seq_len_ctx->last_value.i;

    COPY_TIMER (sam_load_groups_add_qnames);
}

// FLAGS: get multi_segs, is_first, is_last from SAM_FLAGS. 
static inline void sam_load_groups_add_flags (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g) 
{
    START_TIMER;

    STR0(snip);
    LOAD_SNIP(SAM_FLAG); 

    // in SAG_SA, FLAG is segged as a SPECIAL
    if (IS_SAG_SA) {
        ASSERT0 (snip[0]==SNIP_SPECIAL && snip[1]==SAM_SPECIAL_pull_from_sag, "expecting FLAG to be SAM_SPECIAL_pull_from_sag");
        snip += 2;
    }

    SamFlags sam_flags = { .value = atoi (snip) };
    g->multi_segs = sam_flags.multi_segs;
    g->is_first   = sam_flags.is_first;
    g->is_last    = sam_flags.is_last;
    
    // case: non-SAG_SA - revcomp is segged in SAM_FLAG. Note: For SA-defined, it is always 0, and we retrieve it later from OPTION_SA_STRAND of the primary alignment
    if (!IS_SAG_SA)
        g->revcomp = sam_flags.rev_comp; 

    COPY_TIMER (sam_load_groups_add_flags);
}

static inline void sam_load_groups_add_seq (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g, bytes start, bytes after, uint8_t **next)
{
    START_TIMER;

    if (!CTX(SAM_SQBITMAP)->is_loaded) return; // sequence is skipped

    ASSERT (vb->seq_len <= MAX_SA_SEQ_LEN, "PRIM/%u: while loading SAGs: g->seq_len=%u > MAX_SA_SEQ_LEN=%u", 
            plsg->vblock_i, vb->seq_len, MAX_SA_SEQ_LEN);

    // since are start/after/next are in bases vs start, we will ignore the pointer value and only consider the distances
    g->seq     = *next - start; // in bases
    g->seq_len = vb->seq_len;   // in bases

    // reconstruct SEQ to vb->textual_seq, needed by sam_piz_prim_add_QUAL
    buf_alloc (vb, &vb->textual_seq, 0, vb->seq_len, char, 0, "textual_seq");

    SWAP (vb->textual_seq, vb->txt_data); // careful as this doesn't modify buffer_list 
    reconstruct_from_ctx (vb, SAM_SQBITMAP, 0, RECON_ON);
    buf_verify (vb->txt_data, "SWAP-txt_data-textual_seq");
    SWAP (vb->textual_seq, vb->txt_data);

    vb->textual_seq_str = B1STc (vb->textual_seq);

    ASSERT (vb->textual_seq.len32 == vb->seq_len, "Expecting textual_seq.len=%u == seq_len=%u", vb->textual_seq.len32, vb->seq_len);

    // pack SEQ data into z_file->sag_seq
    BitsP z_sa_seq = (BitsP)&z_file->sag_seq;
    { START_TIMER; 
    sam_seq_pack (vb, z_sa_seq, (*next - start) * 2/*2 bits per base*/, B1STc(vb->textual_seq), vb->seq_len, false, false, HARD_FAIL); 
    COPY_TIMER (sam_load_groups_add_seq_pack); }

    *next += vb->seq_len; // ignoring the pointer value, we only consider (*next - start) to be the name of bases (not bytes) in the buffer

    ASSERT (after >= *next, "PRIM/%u: plsg->comp_seq_len=%u exceeded: (after-next)=%"PRId64, 
            vb->plsg_i, plsg->comp_qual_len, (int64_t)(after - *next));

    vb->textual_seq.len = 0;

    COPY_TIMER (sam_load_groups_add_seq);
}

// PIZ: loads QUAL data of PRIM, and compresses it for in-memory storage, using fast huffman or arith compression.
static inline void sam_load_groups_add_qual (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g,
                                             bytes start, bytes after, uint8_t **next)
{
    START_TIMER;
   
    if (!CTX(SAM_QUAL)->is_loaded) { 
        g->no_qual = true; // prevent attempting to reconstruct QUAL if its not loaded
        return;
    }

    // reconstructs into vb->txt_data. note that in PRIM we segged the SPECIAL into QUALSA, not QUAL.
    // so this reconstructs as a LOOKUP from LT_CODEC / LT_BLOB, without going through sam_piz_special_QUAL.
    rom qual = BAFTtxt;
    reconstruct_from_ctx (vb, SAM_QUAL, 0, RECON_ON); 
    uint32_t qual_len = BAFTtxt - qual;
    uint32_t comp_len = 0;
    
    g->no_qual = IS_ASTERISK(qual);
    if (g->no_qual) goto done; // no quality data for this group - we're done

    ASSERT (qual_len == vb->seq_len, "Expecting qual_len=%u == vb->seq_len=%u", qual_len, vb->seq_len);
    
    // back comp: old file with no SEC_HUFFMAN section for SAM_QUAL - create it based on first QUAL reconstructed
    if (!huffman_exists (SAM_QUAL)) // note: non-atomic: if it says it exists, then it exists. If it says it doesn't, it might or might not exist.
        huffman_piz_backcomp_produce_qual (STRa(qual)); // also handles thread-safety

    // case: precise: compress directly into z_file->sag_qual
    if (start) {
        comp_len = after - *next;
        huffman_compress (VB, SAM_QUAL, STRa(qual), *next, &comp_len);

        ASSERT (after >= *next + comp_len, "PRIM/%u: plsg->comp_qual_len=%u exceeded: (after-next)=%"PRId64, 
                vb->plsg_i, plsg->comp_qual_len, (int64_t)(after - (*next + comp_len)));

        g->qual = *next - start;
        g->qual_comp_len = comp_len; 
        *next += comp_len;  
    }
    
    // case: we don't have precise: arith-compressed qual in vb_qual_buf - to be copied to z_file->sag_qual later
    // note: huffman_piz_read_all doen't produce a default huffman for QUAL bc we fallback on arith
    else {
        comp_len = huffman_get_theoretical_max_comp_len (SAM_QUAL, qual_len);
        buf_alloc (vb, &vb_qual_buf, comp_len, 0, uint8_t, 1, NULL); // likely already allocated in sam_load_groups_add_grps

        huffman_compress (VB, SAM_QUAL, STRa(qual), BAFT8(vb_qual_buf), &comp_len);

        g->qual = vb_qual_buf.len; // relative to the VB. we will update later to be relative to z_file->sag_qual.
        g->qual_comp_len = comp_len;  
        vb_qual_buf.len += comp_len;
    }

done:
    COPY_TIMER (sam_load_groups_add_qual);
}

// populates the cigar data of one alignment
static void sam_load_groups_add_aln_cigar (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g, SAAln *a, bool is_first_aln,
                                           bool is_all_the_same_LOOKUP, bool is_all_the_same_SQUANK,
                                           bytes start, bytes after, uint8_t **next, 
                                           pSTRp(out_cigar)) // optional out
{
    START_TIMER;

    ContextP ctx = CTX (OPTION_SA_CIGAR);

    STR0(snip);
    WordIndex word_index = WORD_INDEX_NONE;
    
    // case: word_list.len==0 if all snips are in local, so dictionary was an all-the-same SNIP_LOOKUP, therefore dropped by ctx_drop_all_the_same, possibly dropping the SEC_DICT as well.
    if (ctx->word_list.len32)
        word_index = ctx_get_next_snip (VB, ctx, false, pSTRa(snip)); 

    bool is_lookup = is_all_the_same_LOOKUP || (snip_len && *snip == SNIP_LOOKUP);
    bool is_squank = !is_lookup && (is_all_the_same_SQUANK || (snip_len > 2 && snip[0] == SNIP_SPECIAL && snip[1] == SAM_SPECIAL_SQUANK));

    if (!start)
        buf_alloc (vb, &vb_cigars_buf, 0, 128 KB, char, 0, "vb_cigars_buf"); // initial alloction
        
    // case: the CIGAR is in local of this vb we need to copy it (compressed) as the "SA Group loader" VB will be soon released. 
    if (is_lookup || is_squank) {

        if (is_lookup)
            ctx_get_next_snip_from_local (VB, ctx, pSTRa(snip)); 

        else { // squank - temporarily reconstruct into txt_data - just for compressing. note: prim cigar is never segged as squank
            cigar_special_SQUANK (VB, ctx, &snip[2], snip_len-2, NULL/*reconstruct to vb->scratch*/, true);
            snip     = vb->scratch.data;
            snip_len = vb->scratch.len;
        }

        // case: we know the exact size of huffman for this VB, so we write directly to z_file 
        if (start) {
            a->cigar.piz.index = *next - start;
    
            uint32_t comp_len = MIN_(after - *next, 0xffffffffULL); // careful that it doesn't go beyond 32b

            *next += nico_compress_textual_cigar (VB, OPTION_SA_CIGAR, STRa(snip), *next, comp_len);

            ASSERT (after >= *next, "PRIM/%u: plsg->comp_cigars_len=%u exceeded: (after-next)=%"PRId64, 
                    vb->plsg_i, plsg->comp_qnames_len, (int64_t)(after - *next));
        }

        // case: ZIP used a different version of huffman: we write to VB buffer to be added to z_file later
        else {
            a->cigar.piz.index = vb_cigars_buf.len;  
            buf_alloc (vb, &vb_cigars_buf, snip_len * 4/*plenty*/, 0, char, CTX_GROWTH, "vb_cigars_buf");

            vb_cigars_buf.len32 += nico_compress_textual_cigar (VB, OPTION_SA_CIGAR, STRa(snip), BAFT8(vb_cigars_buf), BFREE8(vb_cigars_buf));
        }

        a->cigar.piz.is_word = false; // cigar is in z_file->sag_cigars (long, possibly compressed)
        a->cigar.piz.len     = snip_len;
    }

    // case: the cigar is in dict - we therefore don't store - reconstruction will copy it from the dictionary
    else {
        ASSERT (IN_RANGE(word_index, 0, ctx->word_list.len32), "word_index=%d of %s ∉ [0,%d) snip_len=%u snip[0]=%u snip=\"%.*s\"", 
                word_index, ctx->tag_name, ctx->word_list.len32, snip_len, snip_len ? (uint8_t)snip[0] : 0, STRf(snip));

        a->cigar.piz.is_word = true; // cigar is in ZCTX(OPTION_SA_CIGAR).dict (short cigar)
        a->cigar.piz.index   = word_index;
    }

    if (out_cigar) STRset (*out_cigar, snip);
    if (is_squank) buf_free (vb->scratch);
}

static void sam_load_groups_add_cigars (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g, bytes start, bytes after, uint8_t **next) // out
{
    START_TIMER;
    decl_ctx (OPTION_SA_CIGAR);

    if (!ctx->is_loaded) return; // skipped

    bool is_all_the_same_LOOKUP = !ctx->word_list.len && ((ctx->dict.len && *B1STc(ctx->dict)==SNIP_LOOKUP) || !ctx->dict.len);
    bool is_all_the_same_SQUANK = !ctx->word_list.len && ctx->dict.len >= 3 && 
                                  *B1STc(ctx->dict)==SNIP_SPECIAL && *Bc(ctx->dict, 1)==SAM_SPECIAL_SQUANK;

    SAAln *prim_aln = B(SAAln, z_file->sag_alns, g->first_aln_i); // first Aln of group

    // primary alignment CIGAR 
    STR(prim_cigar);
    sam_load_groups_add_aln_cigar (vb, plsg, g, prim_aln, ZGRP_I(g)==0, is_all_the_same_LOOKUP, is_all_the_same_SQUANK, start, after, next, pSTRa(prim_cigar));

    // calculate seq_len from primary CIGAR, and copy cigar to vb->textual_cigar
    sam_cigar_analyze (vb, STRa(prim_cigar), false, &vb->seq_len);

    // populate SAAln.cigar with the CIGAR word index, for the non-primary alignments of this group
    for (uint32_t aln_i=1; aln_i < g->num_alns; aln_i++) 
        sam_load_groups_add_aln_cigar (vb, plsg, g, prim_aln + aln_i, false, is_all_the_same_LOOKUP, is_all_the_same_SQUANK, start, after, next, 0, 0);

    COPY_TIMER (sam_load_groups_add_cigars);
}

static inline bool reconstruct_to_solo_aln (VBlockSAMP vb, Did did_i, pSTRp(prev), bool check_copy, pSTRp(recon))
{
    decl_ctx (did_i);

    if (!ctx->is_loaded) return false;

    *recon = BAFTtxt;
    reconstruct_from_ctx (VB, did_i, 0, RECON_ON); 
    *recon_len = BAFTtxt - *recon;

    // if identical, UR and UB, CR and CB, point to the same data in solo_data
    if (check_copy && str_issame (*recon, *prev)) {
        *recon = NULL;
        return true; // note: no need to update prev, as its identical
    }

    STRset (*prev, *recon);
    return true;
}

static inline void sam_load_groups_add_solo_data (VBlockSAMP vb, PlsgVbInfo *plsg, SoloAln *solo_aln, bytes start, bytes first, bytes after, uint8_t **next)
{
    START_TIMER;

    // get the index of the tag in this alignment's AUX container. we're really interested if they exist or not.
    ContainerPeekItem idxs[NUM_SOLO_TAGS] = SOLO_CON_PEEK_ITEMS;
    container_peek_get_idxs (VB, CTX(SAM_AUX), ARRAY_LEN(idxs), idxs, &vb->aux_con, true);

    STR0(prev); // pointer into txt_data of previous tag reconstructed        

    // set index
    if (start) { 
        solo_aln->index = U64to40 (*next - start);
        ASSERT0 (*next - start <= MAX_SA_SOLO_DATA_INDEX, "sag_alns for Solo is too big"); // this should never happen in PIZ as it was already tested in ZIP
    }

    else
        solo_aln->index = U64to40 (vb_solo_data_buf.len);

    // compress each solo tag
    for (SoloTags solo=0; solo < NUM_SOLO_TAGS; solo++)
        if (idxs[solo].idx != -1) {
            STR(recon);
            if (!reconstruct_to_solo_aln (vb, solo_props[solo].did_i, pSTRa(prev), solo_props[solo].maybe_same_as_prev, pSTRa(recon)))
                continue; // not loaded

            ASSERTNOTZERO (recon_len);
            solo_aln->field_uncomp_len[solo] = recon_len;

            // case: field is a copy of the previous field
            if (!recon) {      
                solo_aln->field_comp_len[solo] = solo_aln->field_comp_len[solo-1];
                solo_aln->field_comp_len[solo-1] = 0;     
            }

            // case: precise - write directly to z_file
            else if (start) { 
                uint32_t comp_len = after - *next;
                *next += huffman_compress (VB, solo_props[solo].did_i, STRa(recon), *next, &comp_len);
                solo_aln->field_comp_len[solo] = comp_len;

                ASSERT (after >= *next, "PRIM/%u: plsg->comp_solo_data_len=%u exceeded: (after-next)=%"PRId64, 
                        vb->plsg_i, plsg->comp_solo_data_len, (int64_t)(after - *next));
            }

            // case: write to vb buffer, to be copied to z_file later
            else { 
                uint32_t comp_len = huffman_get_theoretical_max_comp_len (solo_props[solo].did_i, recon_len);
                buf_alloc (vb, &vb_solo_data_buf, comp_len, 0, uint8_t, 1, NULL); // likely already allocated in sam_load_groups_add_grps

                vb_solo_data_buf.len += huffman_compress (VB, solo_props[solo].did_i, STRa(recon), BAFT8(vb_solo_data_buf), &comp_len);
                solo_aln->field_comp_len[solo] = comp_len;
            }
        }

    COPY_TIMER (sam_load_groups_add_solo_data);
}

static inline WordIndex v14_chrom_node_index (VBlockSAMP vb, SAAln *prim_aln)
{
    ctx_get_snip_by_word_index (CTX(OPTION_SA_RNAME), prim_aln->rname, vb->chrom_name);
    
    WordIndex wi = ctx_get_word_index_by_snip (SAM_RNAME, STRa(vb->chrom_name)); // convert OPTION_SA_RNAME word_index to RNAME word_index
    ASSERT (wi != WORD_INDEX_NONE, "Cannot find rname=%u in context RNAME", prim_aln->rname);

    return wi;
}

// Alns - populate RNAME, POS, MAPQ, STRAND and NM (CIGAR is added sam_load_groups_add_grps)
static inline void sam_load_groups_add_SA_alns (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps) 
{
    START_TIMER;

    ContextP sa_ctx    = LOADED_CTX(OPTION_SA_Z);
    ContextP rname_ctx = CTX(OPTION_SA_RNAME); // possibly not loaded

    ASSERT (sa_ctx->local.len == plsg->num_grps, "PRIM/%u: Mismatch in number of PRIM lines: expecting OPTION_SA_Z.local.len=%"PRIu64" == num_lines=%u OPTION_SA_Z.is_loaded=%s",
            plsg->vblock_i, sa_ctx->local.len, plsg->num_grps, TF(sa_ctx->is_loaded));
    
    ASSERT (sa_ctx->ltype == LT_UINT8 || sa_ctx->ltype == LT_UINT16, "Unexpected OPTION_SA_Z.ltype=%s", lt_name (sa_ctx->ltype));

    for (uint32_t grp_i=0; grp_i < plsg->num_grps; grp_i++) {
        uint64_t first_aln_i       = grp_i ? vb_grps[grp_i-1].first_aln_i + vb_grps[grp_i-1].num_alns : plsg->first_aln_i;
        ASSERT (first_aln_i <= MAX_SA_GRP_ALNS, "PRIM/%u: while loading SAGs vb_grp_i=%u: first_aln_i%"PRIu64" > MAX_SA_GRP_ALNS=%"PRIu64, 
                plsg->vblock_i, grp_i, first_aln_i, MAX_SA_GRP_ALNS);

        vb_grps[grp_i].first_aln_i = first_aln_i;

        uint32_t grp_num_alns = (sa_ctx->ltype == LT_UINT8) ? *B8(sa_ctx->local, grp_i) : *B16(sa_ctx->local, grp_i);

        ASSERT (grp_num_alns <= (uint32_t)MAX_SA_NUM_ALNS, "PRIM/%u: while loading SAGs vb_grp_i=%u: grp_num_alns=%d > MAX_SA_NUM_ALNS=%u",
                plsg->vblock_i, grp_i, grp_num_alns, MAX_SA_NUM_ALNS);

        vb_grps[grp_i].num_alns = grp_num_alns;
    }

    int32_t vb_num_alns = (vb_grps[plsg->num_grps-1].first_aln_i - plsg->first_aln_i) + vb_grps[plsg->num_grps-1].num_alns;

    ASSERT (vb_num_alns == plsg->num_alns, "PRIM/%u: Alignment count mismatch: expecting vb_num_alns=%d == plsg->num_alns=%u",
            plsg->vblock_i, vb_num_alns, plsg->num_alns);

    uint32_t sa_rname_len = rname_ctx->word_list.len32;

    ARRAY (SAAln, alns, z_file->sag_alns);

    for (uint32_t grp_i=0; grp_i < plsg->num_grps; grp_i++) 
        for (uint32_t aln_i=vb_grps[grp_i].first_aln_i; aln_i < vb_grps[grp_i].first_aln_i + vb_grps[grp_i].num_alns; aln_i++) {
            SAAln *a = &alns[aln_i];

            // RNAME
            if (rname_ctx->is_loaded) { // RNAME not skipped
                reconstruct_from_ctx (vb, OPTION_SA_RNAME, 0, RECON_OFF);  // don't reconstruct, only set last_value to word_index
                a->rname = rname_ctx->last_value.i;         // WordIndex into OPTION_SA_RNAME, NOT RNAME!

                ASSERT (a->rname >= 0 && a->rname < sa_rname_len, 
                        "While loading aln_i=%u: rname=%u out of range: OPTION_SA_RNAME.word_len.len=%u", aln_i, a->rname, sa_rname_len);
            }

            // POS
            if (CTX(OPTION_SA_POS)->is_loaded) {
                reconstruct_from_ctx (vb, OPTION_SA_POS, 0, RECON_OFF);    // don't reconstruct, only set last_value to pos
                a->pos = CTX(OPTION_SA_POS)->last_value.i;
            }

            // STRAND - nodes initialized in sam_seg_0X_initialize to 0="-" 1="+"
            if (CTX(OPTION_SA_STRAND)->is_loaded) {
                reconstruct_from_ctx (vb, OPTION_SA_STRAND, 0, RECON_OFF); // don't reconstruct, only set last_value to word_index
                a->revcomp = !CTX(OPTION_SA_STRAND)->last_value.i;
            }
            
            // MAPQ
            if (CTX(OPTION_SA_MAPQ)->is_loaded) { // MAPQ is not skipped
                reconstruct_from_ctx (vb, OPTION_SA_MAPQ, 0, RECON_OFF);   // don't reconstruct, only set last_value to pos
                a->mapq = CTX(OPTION_SA_MAPQ)->last_value.i;
            }

            // NM
            if (CTX(OPTION_SA_NM)->is_loaded &&   // NM is not skipped
                !(segconf.SA_NM_by_CIGAR_X && aln_i == vb_grps[grp_i].first_aln_i)) { // in case of SA_NM_by_CIGAR_X, we get the NM of the primary alignment from the its CIGAR 
                reconstruct_from_ctx (vb, OPTION_SA_NM, 0, RECON_OFF);     // don't reconstruct, only set last_value to pos
                a->nm = CTX(OPTION_SA_NM)->last_value.i;
            }
        }

    // group revcomp is the revcomp of its primary alignment
    for (uint32_t grp_i=0; grp_i < plsg->num_grps; grp_i++) 
        vb_grps[grp_i].revcomp = alns[vb_grps[grp_i].first_aln_i].revcomp;

    COPY_TIMER (sam_load_groups_add_SA_alns);
}

// add AS:i and MC:Z
static inline void sam_load_groups_add_aux (VBlockSAMP vb, Sag *g) 
{
    // get the index of the tag in this alignment's AUX container. we're really interested if they exist or not.
    ContainerPeekItem idxs[2] = { { _OPTION_AS_i, -1 }, { _OPTION_MC_Z, -1 } };
    container_peek_get_idxs (VB, CTX(SAM_AUX), ARRAY_LEN(idxs), idxs, &vb->aux_con, true);

    // set AS.last_value, unless AS:i is skipped entirely
    if (segconf.sag_has_AS          && // depn lines' AS:i was segged against prim 
        CTX(OPTION_AS_i)->is_loaded && // we didn't skip AS:i due to subsetting command line options
        idxs[0].idx != -1) {           // this line has AS:i

        reconstruct_from_ctx (VB, OPTION_AS_i, 0, RECON_OFF);
        g->as = CAP_SA_AS (vb->last_int(OPTION_AS_i));  // capped at 0 and 65535
    }

    // if line contains MC:Z - reconstruct it into txt_data, as our mate's CIGAR might copy it
    // note: a better way to do this would be using reconstruct_to_history
    if (!IS_SAG_SA && CTX(OPTION_MC_Z)->is_loaded && idxs[1].idx != -1/*this line has MC:Z*/)  
        reconstruct_from_ctx (VB, OPTION_MC_Z, 0, RECON_ON);
}

// sets SAM_POS.last_value, SAM_FLAG.last_value, vb->chrom_node_index, vb->chrom_name, vb->seq_len, vb->ref_consumed
static inline void sam_load_groups_set_CHROM_POS_FLAG (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g, SAGroup grp_i, ContextP seq_len_ctx)
{
    // if this is an SA:Z-defined SAG, load primary alignment
    if (IS_SAG_SA) {
        ASSERT (IN_RANGE(g->first_aln_i, plsg->first_aln_i, z_file->sag_alns.len),
                "grp_i=%u z_grp_i=%u has first_aln_i=%"PRIu64" ∉ [0,%"PRId64")", 
                grp_i, ZGRP_I(g), (uint64_t)g->first_aln_i, z_file->sag_alns.len);

        SAAln *prim_aln = B(SAAln, z_file->sag_alns, g->first_aln_i); // first Aln of group
                
        // populate vb with alignment details, needed for QUAL reconstruction for several codecs
        vb->last_int(SAM_POS) = prim_aln->pos;
        CTX(SAM_FLAG)->last_value.i = prim_aln->revcomp ? SAM_FLAG_REV_COMP : 0; 

        ASSERT (IN_RANGE(prim_aln->rname, 0, CTX(OPTION_SA_RNAME)->word_list.len32), 
                "%s: rname=%d out of range: OPTION_SA_RNAME.word_len.len=%u. grp_i=%u",
                VB_NAME, prim_aln->rname, CTX(OPTION_SA_RNAME)->word_list.len32, grp_i);

        // note: since v15, OPTION_SA_RNAME is an dict alias of SAM_RNAME, and v14 they were independent
        vb->chrom_node_index = VER(15) ? prim_aln->rname : v14_chrom_node_index (vb, prim_aln); 
    }

    // non-SAG_SA 
    else {
        // needed to reconstruct SEQ, needed for QUAL reconstruction for several codecs - g->revcomp was populated in sam_load_groups_add_flags
        CTX(SAM_FLAG)->last_value.i = g->revcomp ? SAM_FLAG_REV_COMP : 0;

        // set vb->chrom_* (needed to reconstruct SEQ)
        reconstruct_from_ctx (VB, SAM_RNAME, 0, RECON_OFF);
        vb->chrom_node_index = vb->last_int (SAM_RNAME);

        // set POS.last_value (which also requires RNEXT and PNEXT)
        reconstruct_from_ctx (VB, SAM_POS,   0, RECON_OFF);
        reconstruct_from_ctx (VB, SAM_RNEXT, 0, RECON_OFF); // needed by PNEXT...
        reconstruct_from_ctx (VB, SAM_PNEXT, 0, RECON_OFF); // store in history, in case POS of future-line mates needs it (if copy_buddy(PNEXT))
        
        // need to set before reconstructing CIGAR: if seq_len is carried by a QNAME item, set the last value here
        if (seq_len_ctx) ctx_set_last_value (VB, seq_len_ctx, (int64_t)g->seq_len);

        // analyze CIGAR (setting vb->seq_len, vb->ref_consumed etc)
        reconstruct_from_ctx (VB, SAM_CIGAR, 0, RECON_OFF);
    }

    ctx_get_snip_by_word_index (CTX(SAM_RNAME), vb->chrom_node_index, vb->chrom_name);
}

// SEQ - ACGT-pack uncompressed sequence directly to z_file->sag_seq
// Alns - populate CIGAR the cigar field of the alignments
// Grp  - populate seq, seq_len, num_alns
static inline void sam_load_groups_add_grps (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps, void *vb_alns) 
{
    START_TIMER;

    // txt_data is used to retrieve qnames, solo data and MC:Z. Note we don't reset between lines, because some fields rely on previous lines
    buf_alloc (vb, &vb->txt_data, 0, vb->recon_size, char, 0, "txt_data"); 

    // if MC:Z attempts to copy an earlier CIGAR, it will get an empty snip bc CIGAR is not in txt_data. 
    // That's ok, bc this means this MC:Z will not be used to reconstruct its mate's CIGAR as its mate already appeared (and we don't use MC:Z for anything else during pre-processing)
    CTX(SAM_CIGAR)->empty_lookup_ok = true; 

    ContextP seq_len_ctx = segconf.seq_len_dict_id.num ? ECTX(segconf.seq_len_dict_id) : NULL;

    #define INIT(x, cond, precise_cond, imprecise_buf)                                                \
        bytes start_##x=0, first_##x=0, after_##x=0; uint8_t *next_##x=0;                       \
        if (cond && precise_cond) {                                                             \
            start_##x = B1ST8(z_file->sag_##x);                                                 \
            first_##x = B8(z_file->sag_##x, plsg->comp_##x##_start); /* first x for this plsg */\
            after_##x = first_##x + plsg->comp_##x##_len; /* after x data of this plsg */       \
            next_##x  = (uint8_t *)first_##x;                                                   \
        }                                                                                       \
        else if (cond) {                                                                        \
            ASSERTNOTINUSE (imprecise_buf);                                                     \
            uint32_t imprecise_allocation = z_file->sag_##x.len / plsg_info.len + 1 MB/*arith*/;\
            buf_alloc (vb, &imprecise_buf, imprecise_allocation, 0, char, 1.1, #imprecise_buf); \
        }
                    /* cond for writing to z_file | buffer if writing locally */
    INIT(qnames,    true,        is_precise_qnames, vb_qnames_buf);
    INIT(cigars,    IS_SAG_SA,   is_precise_cigars, vb_cigars_buf);
    INIT(seq,       true,        true,              vb->scratch/*dummy - code never reached*/); // note: all *_seq variables are in terms of bases, not bytes
    INIT(qual,      true,        is_precise_qual,   vb_qual_buf);
    INIT(solo_data, IS_SAG_SOLO, is_precise_solo,   vb_solo_data_buf);
    
    #undef INIT

    if (IS_SAG_SA)
        sam_load_groups_add_SA_alns (vb, plsg, vb_grps);

    // every alignment in the Primary VB represents a grp
    for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) {
        Sag *g = &vb_grps[grp_i];

        vb->line_i = grp_i;
        sam_reset_line (VB);

        reconstruct_from_ctx (VB, SAM_BUDDY, 0, RECON_OFF/*this means: don't consume QNAME (see sam_piz_special_COPY_BUDDY)*/); // set buddy

        sam_load_groups_add_qnames (vb, plsg, vb_grps, g, start_qnames, after_qnames, &next_qnames); 
        
        sam_load_groups_add_flags(vb, plsg, g);

        sam_load_groups_set_CHROM_POS_FLAG (vb, plsg, g, grp_i, seq_len_ctx); // needed for QUAL reconstruction in several codecs

        switch (segconf.sag_type) {
            case SAG_BY_SA: 
                sam_load_groups_add_cigars (vb, plsg, g, start_cigars, after_cigars, &next_cigars);
                break;

            case SAG_BY_CC : 
                ((CCAlnP)vb_alns)[grp_i] = (CCAln){ .rname = vb->chrom_node_index, .pos = vb->last_int(SAM_POS) };
                goto NH;

            case SAG_BY_SOLO :
                sam_load_groups_add_solo_data (vb, plsg, &((SoloAlnP)vb_alns)[grp_i], start_solo_data, first_solo_data, after_solo_data, &next_solo_data);
                goto NH;

            case SAG_BY_NH : NH :
                reconstruct_from_ctx (VB, OPTION_NH_i, 0, RECON_OFF);
                g->num_alns = vb->last_int (OPTION_NH_i);
                break;

            case SAG_BY_FLAG : break;

            default : ABORT ("invalid sag_type=%d", segconf.sag_type);   
        }

        sam_load_groups_add_seq (vb, plsg, g, start_seq, after_seq, &next_seq);

        sam_load_groups_add_qual (vb, plsg, g, start_qual, after_qual, &next_qual);

        sam_load_groups_add_aux (vb, g); // MC:Z, AS:i + consumes AUX container
    }

    #define FINALIZE(cond, x, imprecise_buf)                                                \
        /* sanity checks in case if IS_PRECISE */                                           \
        ASSERT (after_##x == next_##x || !(cond),                                           \
                "PRIM/%u: Expecting compressed %s for this prim vb to have length=%u, but its length=%"PRIu64, \
                vb->plsg_i, #x, plsg->comp_##x##_len, (uint64_t)(next_##x - first_##x));    \
        if (!start_##x) plsg->comp_##x##_len = imprecise_buf.len32; // update to compressed lengths actual if not precise

    FINALIZE (CTX(SAM_QNAME)->is_loaded,    qnames,    vb_qnames_buf);
    FINALIZE (IS_SAG_SA,                    cigars,    vb_cigars_buf);
    FINALIZE (CTX(SAM_SQBITMAP)->is_loaded, seq,       vb->scratch/*dummy, code not reached*/)
    FINALIZE (CTX(SAM_QUAL)->is_loaded,     qual,      vb_qual_buf);
    FINALIZE (sam_is_solo_loaded (vb),      solo_data, vb_solo_data_buf)
    #undef FINALIZE
    
    COPY_TIMER (sam_load_groups_add_grps);
}

// Grp loader compute thread: copy to z_file->sag_*
static inline void sam_load_groups_move_comp_to_zfile (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps, void *vb_alns) 
{
    START_TIMER;
    
    // copy buffers in arbitrary order, based on mutex availability. no need to copy if precise length was known and hence data was written directly.
    bool qual_done      = is_precise_qual;  
    bool cigars_done    = is_precise_cigars || !IS_SAG_SA; 
    bool solo_data_done = is_precise_solo   || !IS_SAG_SOLO;
    bool qnames_done    = is_precise_qnames; 

    uint64_t start_qnames=0, start_qual=0, start_cigars=0, start_solo_data=0;

    while (!qual_done || !cigars_done || !solo_data_done || !qnames_done) {
        bool achieved_something = false;
        
        // note: we lock as we might realloc (therefore other threads can't access this Buffer concurrently)
        #define APPEND(x, vb_buf, max_index)                                    \
            if (!x##_done && mutex_trylock (copy_##x##_mutex)) {                \
                start_##x = z_file->sag_##x.len;                                \
                buf_append_buf (evb, &z_file->sag_##x, &vb_buf, uint8_t, NULL); \
                x##_done = achieved_something = true;                           \
                ASSERT (z_file->sag_##x.len <= max_index, "PRIM/%u: while loading SAGs: z_file->sag_"#x".len=%"PRIu64" > "#max_index"=%"PRIu64, \
                        plsg->vblock_i, z_file->sag_##x.len, max_index);        \
                mutex_unlock (copy_##x##_mutex);                                \
            }

        APPEND (qual,      vb_qual_buf,      MAX_SA_QUAL_INDEX);
        APPEND (cigars,    vb_cigars_buf,    MAX_SA_CIGAR_INDEX);
        APPEND (solo_data, vb_solo_data_buf, MAX_SA_SOLO_DATA_INDEX);
        APPEND (qnames,    vb_qnames_buf,    MAX_SA_QNAME_INDEX);

        // case: all muteces were locked by other threads - wait a bit
        if (!achieved_something) { 
            START_TIMER;
            usleep (1000); // 1 ms
            COPY_TIMER (sam_load_groups_move_comp_to_zfile_idle);
        }
    }

    // update indices outside of the muteces

    if (!is_precise_qnames) {
        buf_free (vb_qnames_buf);
        for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) 
            vb_grps[grp_i].qname += start_qnames;
    }

    if (!is_precise_qual) {
        buf_free (vb_qual_buf);
        for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) 
            vb_grps[grp_i].qual += start_qual;
    }

    if (!is_precise_solo && IS_SAG_SOLO) {
        buf_free (vb_solo_data_buf);
        SoloAln *solo_alns = (SoloAln *)vb_alns;
        for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) 
            solo_alns[grp_i].index = U64to40(start_solo_data + solo_alns[grp_i].index.lo); // note: index.hi not used within a VB
    }

    if (!is_precise_cigars && IS_SAG_SA) {
        buf_free (vb_cigars_buf);
        SAAln *sa_alns = (SAAln *)vb_alns;
        for (uint32_t aln_i=0; aln_i < plsg->num_alns ; aln_i++) 
            if (!sa_alns[aln_i].cigar.piz.is_word)
                sa_alns[aln_i].cigar.piz.index += start_cigars; 
    }

    COPY_TIMER (sam_load_groups_move_comp_to_zfile);
}

// entry point of compute thread of sag loading
static void sam_load_groups_add_one_prim_vb (VBlockP vb_)
{
    START_TIMER;
    VBlockSAMP vb = (VBlockSAMP)vb_;

    piz_uncompress_all_ctxs (VB, PUR_SAM_LOAD_SAG);

    buf_free (vb->z_data); // we will use z_data for vb_qual_buf to avoid allocating a separate buffer

    PlsgVbInfo *plsg = B(PlsgVbInfo, plsg_info, vb->plsg_i);
    Sag *vb_grps     = B(Sag, z_file->sag_grps, plsg->first_grp_i);
    void *vb_alns    = IS_SAG_SA   ? (void *)B(SAAln,   z_file->sag_alns, plsg->first_aln_i) // multiple alignments per group
                     : IS_SAG_CC   ? (void *)B(CCAln,   z_file->sag_alns, plsg->first_grp_i) // single "alignment" data per group
                     : IS_SAG_SOLO ? (void *)B(SoloAln, z_file->sag_alns, plsg->first_grp_i) // single "alignment" data per group
                     :               NULL;

    vb_grps[0].first_grp_in_vb = true;

    sam_load_groups_add_grps (vb, plsg, vb_grps, vb_alns);

    // move compressed QNAME, QUAL, CIGAR to from vb to z_file 
    sam_load_groups_move_comp_to_zfile (vb, plsg, vb_grps, vb_alns); 

    __atomic_fetch_add (&num_prim_vbs_loaded, (int)1, __ATOMIC_ACQ_REL);

    vb_set_is_processed (VB); // tell dispatcher this thread is done and can be joined.

    COPY_TIMER_EVB (sam_load_groups_add_one_prim_vb);                          
}

// PIZ main thread - dispatches a compute thread do read SA Groups from one PRIM VB. returns true if dispatched
bool sam_piz_dispatch_one_load_sag_vb (Dispatcher dispatcher)
{
    // check if we're done
    if (next_plsg_i == plsg_info.len) goto done_loading;

    PlsgVbInfo *plsg = B(PlsgVbInfo, plsg_info, next_plsg_i++);

    VBlockP vb = dispatcher_generate_next_vb (dispatcher, plsg->vblock_i, SAM_COMP_PRIM);
    vb->comp_i          = SAM_COMP_PRIM;
    vb->preprocessing   = true; // we're done dispatching preprocessing VBs, however, preprocessing VB compute threads may still be running.
    vb->show_containers = false; 

    piz_read_one_vb (vb, false); // only the needed contexts - filtered by sam_piz_is_skip_section
    dispatcher_compute (dispatcher, sam_load_groups_add_one_prim_vb);

    return true;

done_loading: // note: at this point, gencomp load VBs are still running in compute threads
    return false;
}

// PIZ main thread: after joining a preprocessing VB
void sam_piz_after_preproc_vb (VBlockP vb)
{
    z_file->num_preproc_vbs_joined++;

    // case: done loading - all groups are loaded
    if (sections_get_num_vbs (SAM_COMP_PRIM) == z_file->num_preproc_vbs_joined) {
        if (flag.show_sag) {
            sam_show_sag();
            if (is_genocat) exit_ok;
        }

        sam_gencomp_trim_memory();
    }

    if (flag_is_show_vblocks (PREPROCESSING_TASK_NAME) || flag_is_show_vblocks (PIZ_TASK_NAME))
        iprintf ("LOADED_SA(id=%d) vb=%s\n", vb->id, VB_NAME);
}

// PIZ main thread: after joining call preprocessing VBs
void sam_piz_preproc_finalize (Dispatcher dispatcher)
{
    // restore globals
    flag = save_flag;  // also resets flag.preprocessing
    dispatcher_set_task_name (dispatcher, PIZ_TASK_NAME);
}

static void sam_sag_load_alloc_z (BufferP buf, bool is_precise, uint64_t precise_len, uint64_t estimated_len, MutexP mutex, rom buf_name)
{
    // if total length is known precisely, we allocate it, and each thread puts the data directly into its location
    if (is_precise) 
        buf_alloc_exact (evb, *buf, precise_len, uint8_t, buf_name);    

    // if total length is only an estimate, each thread will put the data into a VB buffer, and it will be added to the z_file buffer under mutex lock
    else { 
        buf_alloc (evb, buf, 0, estimated_len, uint8_t, 1, buf_name); 
        mutex_initialize (*mutex); 
    }
}

static void sam_sag_load_set_is_precise (void)
{
    is_precise_qnames = ZCTX(SAM_QNAME)->huffman.param == HUFF_PRODUCED_BY_ZIP; // huffman exists for files since 15.0.65
    is_precise_qual   = ZCTX(SAM_QUAL) ->huffman.param == HUFF_PRODUCED_BY_ZIP; // huffman exists for files since 15.0.69
    is_precise_cigars = false; // nico is chewed during segging of the first prim VB, so we can't calculate length in seg (bug 1147). 15.0.68 files had huffman instead of nico, and older files had no compression

    // for solo - we declare precise is *all* solo fields used by this file are precise (note: solo huffmans exist since 15.0.68 for the solo fields defined at that point)
    if (IS_SAG_SOLO) {
        is_precise_solo = true; // optimistic

        for (SoloTags solo=0; solo < NUM_SOLO_TAGS; solo++) {
            ContextP zctx = ZCTX(solo_props[solo].did_i); 

            if (zctx->z_data_exists && zctx->huffman.param == HUFF_PRODUCED_BY_PIZ) {
                is_precise_solo = false;
                break;
            }
        }
    }

    if (flag.show_sag)
        #define PA(x) (x ? "PRECISE" : "APPEND_VB")
        iprintf ("memory allocation: qual=%s%s qnames=%s%s\n", 
                 PA(is_precise_qual),   cond_str (IS_SAG_SA, " cigars=", PA(is_precise_cigars)), 
                 PA(is_precise_qnames), cond_str (IS_SAG_SOLO, " solo=", PA(is_precise_solo)));
        #undef PA
}

// PIZ main thread: a callback of piz_after_global_area 
void sam_piz_load_sags (void)
{
    next_plsg_i = num_prim_vbs_loaded = 0; // reset for new z_file

    // piz of a file with gencomp needs at least 2 threads
    if (z_has_gencomp)
        MAXIMIZE (global_max_threads, 2); 

    // cases in which there is no need to load sags
    if (!VER(14) || flag.genocat_no_reconstruct || flag.header_only
        || (flag.one_vb && IS_MAIN (sections_vb_header (flag.one_vb)))) // --one-vb of a MAIN VB - no need for PRIM/DEPN
        return; // no need to load sags

    sam_sag_load_set_is_precise();

    uint32_t num_prim_vbs = sections_get_num_vbs (SAM_COMP_PRIM);
    if (!num_prim_vbs) return; // no actual PRIM lines (either not gencomp file, or only DEPN lines)

    ARRAY_alloc (PlsgVbInfo, plsg, num_prim_vbs, false, plsg_info, evb, "z_file->plsg");

    uint64_t total_seq_len=0, total_qual_len=0, total_cigars_len=0, total_solo_data_len=0, total_qnames_len=0, total_alns=0; // total across all PRIM VBs of the file
    SAGroup total_grps=0;

    // surveys all PRIM VBs in the order of vb_i (i.e. order of creation) (note: this is not the order they appear in the file)
    Section vb_header_sec = NULL;
    for (uint32_t i=0; i < num_prim_vbs; i++) {
        
        sections_get_next_vb_header_sec (SAM_COMP_PRIM, &vb_header_sec);

        ASSERT0 (vb_header_sec->st==SEC_VB_HEADER && IS_PRIM(vb_header_sec), "expecting a PRIM VB Header");

        SectionHeaderVbHeader header = zfile_read_section_header (evb, vb_header_sec, SEC_VB_HEADER).vb_header;

        plsg[i] = (PlsgVbInfo){ .vblock_i             = vb_header_sec->vblock_i,
                                .first_grp_i          = BGEN32 (header.sam_prim_first_grp_i),
                                .first_aln_i          = total_alns,
                                .comp_seq_start       = total_seq_len, // in bases
                                .comp_seq_len         = BGEN32 (header.sam_prim_seq_len), // in bases
                                .num_grps             = vb_header_sec->num_lines,
                                .num_alns             = BGEN32 (header.sam_prim_num_sag_alns),
                                .comp_qual_start      = is_precise_qual   ? total_qual_len : 0, 
                                .comp_qual_len        = is_precise_qual   ? BGEN32 (header.sam_prim_comp_qual_len) : 0, 
                                .comp_qnames_start    = is_precise_qnames ? total_qnames_len : 0,
                                .comp_qnames_len      = is_precise_qnames ? BGEN32 (header.sam_prim_comp_qname_len) : 0,
                                .comp_cigars_start    = (IS_SAG_SA && is_precise_cigars) ? total_cigars_len : 0,
                                .comp_cigars_len      = (IS_SAG_SA && is_precise_cigars) ? BGEN32 (header.sam_prim_comp_cigars_len) : 0,
                                .comp_solo_data_start = (IS_SAG_SOLO && is_precise_solo) ? total_solo_data_len : 0,
                                .comp_solo_data_len   = (IS_SAG_SOLO && is_precise_solo) ? BGEN32 (header.sam_prim_solo_data_len) : 0 };  

        // sum up the precise length is known (word-aligning bit fields for multi-threaded writing), or an estimate if not
        total_grps          += plsg[i].num_grps;
        total_alns          += plsg[i].num_alns;
        total_seq_len       += ROUNDUP32 (plsg[i].comp_seq_len);  // 32 bases = 8 bytes

        total_qnames_len    += is_precise_qnames ? ROUNDUP8 (plsg[i].comp_qnames_len) 
                             : VER2(15,65)       ? BGEN32 (header.sam_prim_comp_qname_len)  // future proof
                             :                     (BGEN32 (header.sam_prim_comp_qname_len) / 3); // prior to 15.0.65, sam_prim_comp_qname_len contains uncompressed length
        
        total_qual_len      += is_precise_qual  ? ROUNDUP8 (plsg[i].comp_qual_len)  // since 15.0.69, comp_qual_len is an exact length of huffman compression (in 68, it was either arith or huffman)
                             : VER(15)          ? (BGEN32 (header.longest_seq_len) * plsg->num_grps / 3) // v15: estimate length: we use longest_seq_len which is available since v15
                             :                    (151 * plsg->num_grps / 3);                            // v14: estimate length: we arbitrarily take seq_len=151

        total_cigars_len    += (!IS_SAG_SA)      ? 0
                             : is_precise_cigars ? ROUNDUP8 (plsg[i].comp_cigars_len)
                             :                     BGEN32 (header.sam_prim_comp_cigars_len); // v14 to 15.0.68 - a variety of older, non-nico compressed size. use as estimate.

        total_solo_data_len += (!IS_SAG_SOLO)    ? 0
                             : is_precise_solo   ? ROUNDUP8 (plsg[i].comp_solo_data_len)
                             : VER2(15,68)       ? BGEN32 (header.sam_prim_solo_data_len) // future proof
                             :                     BGEN32 (header.sam_prim_solo_data_len) / 3;    // prior to 15.0.68, sam_prim_solo_data_len contains uncompressed length        
    }

    z_file->sag_grps.can_be_big = z_file->sag_qnames.can_be_big = z_file->sag_alns.can_be_big = 
    z_file->sag_solo_data.can_be_big = z_file->sag_seq.can_be_big = z_file->sag_qual.can_be_big = 
    z_file->sag_cigars.can_be_big = true; // suppress warning at alloc

    z_file->sag_alns.count = total_alns;

    buf_alloc_exact_zero (evb, z_file->sag_grps, total_grps, Sag, "z_file->sag_grps"); // also sets .len ; zero in case some fields not loaded
    
    // SEQ: we know the precise length     
    buf_alloc_bits_exact (evb, &z_file->sag_seq, total_seq_len * 2, NOINIT, 0, "z_file->sag_seq");   // 2 bits per base

    // QUAL: we know the precise compressed length since 15.0.68
    sam_sag_load_alloc_z (&z_file->sag_qual, is_precise_qual, total_qual_len, total_qual_len, &copy_qual_mutex, "z_file->sag_qual");

    // QNAMES: we know the precise huffman-compressed length since 15.0.65
    sam_sag_load_alloc_z (&z_file->sag_qnames, is_precise_qnames, total_qnames_len, 
                          VER2(15,65) ? total_qnames_len : (total_qnames_len / 3), // estimated length: prior to 15.0.65, sam_prim_comp_qname_len was uncompressed length
                          &copy_qnames_mutex, "z_file->sag_qnames");

    if (IS_SAG_SA) {
        buf_alloc_exact_zero (evb, z_file->sag_alns, total_alns, SAAln, "z_file->sag_alns");       
        
        // CIGARS: currently, we never know cigars precisely
        sam_sag_load_alloc_z (&z_file->sag_cigars, is_precise_cigars, total_cigars_len, 
                              MAX_(total_grps * 3, 1 MB), // very imprecise estimate
                              &copy_cigars_mutex, "z_file->sag_cigars");
    }

    else if (IS_SAG_CC)
        buf_alloc_exact_zero (evb, z_file->sag_alns, total_grps, CCAln, "z_file->sag_alns");       

    else if (IS_SAG_SOLO) { // we know the exact length we are going to store in memory: huffman-compressed since 15.0.68, not compressed for earlier files
        buf_alloc_exact_zero (evb, z_file->sag_alns, total_grps, SoloAln, "z_file->sag_alns");       

        // SOLO DATA: we know the precise huffman-compressed length with 15.0.68     
        sam_sag_load_alloc_z (&z_file->sag_solo_data, is_precise_solo, total_solo_data_len, 
                              VER2(15,68) ? total_solo_data_len : (total_solo_data_len/3), // estimated length: prior to 15.0.68, comp_solo_data_len was uncompressed length
                              &copy_solo_data_mutex, "z_file->solo_data"); 
    }

    save_flag = flag; // save in global variable

    flag.preprocessing = PREPROC_RUNNING; // we're currently dispatching compute threads for preprocessing (= loading SA Groups)
}

uint32_t sam_piz_get_plsg_i (VBIType vb_i)
{
    ASSERTNOTEMPTY (plsg_info);

    ARRAY (PlsgVbInfo, plsg, plsg_info);
    for (int i=0; i < plsg_len; i++)
        if (plsg[i].vblock_i == vb_i)
            return i;

    ABORT ("vb_i=%u is not a PRIM VB", vb_i);
}

// PIZ reconstruction compute thread
Sag *sam_piz_get_prim_vb_first_sa_grp (VBlockSAMP vb)
{
    ASSERTNOTEMPTY (plsg_info);
    PlsgVbInfo *plsg = B(PlsgVbInfo, plsg_info, vb->plsg_i);
    return B(Sag, z_file->sag_grps, plsg->first_grp_i);
}

bool sam_is_sag_loaded (void)
{
    return load_acquire (num_prim_vbs_loaded) == plsg_info.len;
}
