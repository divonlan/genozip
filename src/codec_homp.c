// ------------------------------------------------------------------
//   codec_homp.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// WARNING: THIS FILE CONTAINS A METHOD THAT IS PATENT PENDING.

#include "reconstruct.h"
#include "profiler.h"
#include "compressor.h"

#define TOP_QUAL 'I'

//--------------
// ZIP side
//--------------

static bool is_domq_better (Did did_i)
{
    QualHistType qht = did_i_to_qht (did_i);
    if (qht < 0) return false;

    int used_scores=0, total_count=0, highest_score=0;
    for (int i=0; i < 94; i++) {
        int score_i = segconf.qual_histo[qht][i].count;
        total_count += score_i;

        if (score_i) used_scores++;

        MAXIMIZE (highest_score, score_i);
    }

    return used_scores <= 4 || // Illumina binned
           percent (highest_score, total_count) > 80; // over 80% of scores are the dom
}

static bool codec_homp_qual_data_is_a_fit_for_homp (VBlockP vb, ContextP qual_ctx, LocalGetLineCB qual_callback)
{
    // case: if this is illumina (or illumina-style) binned, it may look like HOMP bc of long dom runs,
    // but we're better off with DOMQ
    if (is_domq_better (qual_ctx->did_i)) 
        return false;

    LocalGetLineCB *seq_callback = (VB_DT(FASTQ) ? fastq_zip_seq : sam_zip_seq);

    #define NUM_HPS_IN_SAMPLE 100  // test at least this number of HPs (or as many as are available)

    uint32_t n_hp_examined = 0;

    for (LineIType line_i=0; line_i < vb->lines.len32 && n_hp_examined < NUM_HPS_IN_SAMPLE; line_i++) {   
        STRw(qual);
        qual_callback (vb, qual_ctx, line_i, pSTRa (qual), CALLBACK_NO_SIZE_LIMIT, NULL);
        
        if (!qual_len) continue;

        STRw(seq);
        seq_callback (vb, NULL, line_i,  pSTRa(seq), CALLBACK_NO_SIZE_LIMIT, NULL);

        for (unsigned i=0; i < seq_len; ) {
            unsigned hp_len = homopolymer_len (STRa(seq), i);

            if (hp_len > 1 && 
                i && i + hp_len < seq_len) { // don't examine homopolymers at the edge of the sequence, as they might have been trimmed

                char prev_qual = 0;

                for (unsigned h=0; h < (hp_len + 1) / 2; h++) {
                    char this_qual   = qual[i + h];
                    char mirror_qual = qual[i + hp_len - 1 - h];

                    if (this_qual != mirror_qual || (prev_qual == TOP_QUAL && this_qual != TOP_QUAL)) {
                        printf ("HOMP: %s: found non-compliant '%c' homopolymer: line_i=%u seq_i=%u hp_len=%u\nSEQ =%.*s\nQUAL=%.*s\n", 
                                VB_NAME, seq[i], line_i, i, hp_len, STRf(seq), seq_len, qual);
                        return false; // found proof of non-homp
                    }

                    prev_qual = this_qual;
                }
                n_hp_examined++;
            }

            i += hp_len; // skip to after homopolymer
        }
    }
    
    if (flag.show_qual) {
        if (n_hp_examined) printf ("HOMP: %s: examined %u homopolymers and all are compliant with HOMP\n", VB_NAME, n_hp_examined);
        else               printf ("HOMP: %s: no homopolymers found\n", VB_NAME);
    }

    return n_hp_examined > 0; // fail if no HPs
}

bool codec_homp_comp_init (VBlockP vb, Did qual_did_i, LocalGetLineCB get_line_cb, bool force)
{
    ContextP qual_ctx = CTX(qual_did_i);
    bool is_ultima = TECH(ULTIMA) || MP(ULTIMA);
    bool is_fit = false;

    if (force || flag.force_qual_codec == CODEC_HOMP) {}

    // All Ultimate generate homp qual
    else if (flag.no_homp                     
         || (!is_ultima && !TECH(UNKNOWN))   // either the TECH or the mapper are an indication of potentially Ultima data 
         || qual_did_i != FASTQ_QUAL/*=SAM_QUAL*/ 
         || (!is_ultima && !(is_fit = codec_homp_qual_data_is_a_fit_for_homp (vb, qual_ctx, get_line_cb)))) // note: we don't test if knowt to be Ultima - HOMP is always better 
        return false;
    
    qual_ctx->ltype     = LT_CODEC;
    qual_ctx->lcodec    = CODEC_HOMP;
    qual_ctx->local_dep = DEP_L1; // yield to other codecs (eg CODEC_OQ) that need to query QUAL before we destroy it
    
    if (TECH(UNKNOWN)) segconf.tech = TECH_ULTIMA; // a HOMP is definitive signature of Ultima

    // show_qual: if is_fit was not run (because is_ultima=true) - run it just to display stats
    if (flag.show_qual && !is_fit)
        codec_homp_qual_data_is_a_fit_for_homp (vb, qual_ctx, get_line_cb);

    return true;
}

COMPRESS (codec_homp_compress)
{
    START_TIMER;

    // case: this is our second entry, after soft-failing. Just continue from where we stopped
    if (!soft_fail) goto do_compress;

    LocalGetLineCB *seq_callback = (VB_DT(FASTQ) ? fastq_zip_seq : sam_zip_seq);
    void (*update_qual_len)(VBlockP, uint32_t, uint32_t) = (VB_DT(FASTQ) ? fastq_update_qual_len : sam_update_qual_len);

    if (ctx->did_i == SAM_QUAL/*==FASTQ_QUAL*/) 
        add_relaxed (z_file->homp_lines, vb->lines.len);

    // first pass - condese qual scores of homopolymers by removing redundant scores
    for_line {   
        STRw(qual);
        get_line_cb (vb, ctx, line_i, pSTRa (qual), CALLBACK_NO_SIZE_LIMIT, NULL);

        if (qual_len <= 1) continue; // including no qual and qual=" "

        char *next_condensed = qual; // copy in-place;

        STRw(seq);
        seq_callback (vb, NULL, line_i,  pSTRa(seq), CALLBACK_NO_SIZE_LIMIT, NULL);

        for (uint32_t i=0; i < seq_len; i++) {
            unsigned hp_len = homopolymer_len (STRa(seq), i);  
            
            if (hp_len > 1) {
                bool condensable = true; // optimistic

                // test condensability (rarely, homopolymer qual may be non-condensable)
                char prev_qual = 0;
                for (unsigned h=0; h < (hp_len + 1) / 2; h++) {
                    char this_qual   = qual[i + h];
                    char mirror_qual = qual[i + hp_len - 1 - h];
                    
                    if ((this_qual != mirror_qual) || (prev_qual == TOP_QUAL && this_qual != TOP_QUAL)) {
                        condensable = false;
                        break;
                    } 

                    prev_qual = this_qual;
                }

                if (condensable) {
                    for (unsigned h=0; h < (hp_len + 1) / 2; h++) 
                        if ((*next_condensed++ = qual[i + h]) == TOP_QUAL) break;
                }
                else {
                    *next_condensed++ = qual[i] | 0x80; // "not-condensed" flag

                    for (unsigned h=1; h < hp_len; h++)
                        *next_condensed++ = qual[i + h];
                }

                i += hp_len - 1; // skip to end of homopolymer
            }
            
            else  
                *next_condensed++ = qual[i];
        }

        uint32_t condensed_len = next_condensed - qual;

        if (condensed_len != qual_len) {
            // shrink qual_len because sub_codec uses the callback to get the condensed qual
            update_qual_len (vb, line_i, condensed_len); 
            ctx->local.len32 -= qual_len - condensed_len;
        }
    }

    // second pass - compress qual data with redundant homopolymer qual scores removed
    ctx->lcodec = CODEC_UNKNOWN;
    header->sub_codec = codec_assign_best_codec (vb, ctx, NULL, SEC_LOCAL); 
    if (header->sub_codec == CODEC_UNKNOWN) header->sub_codec = CODEC_NONE; // really small

do_compress: ({});
    CodecCompress *compress = codec_args[header->sub_codec].compress;
    *uncompressed_len = ctx->local.len32;

    // make sure we have enough memory
    uint32_t min_required_compressed_len = codec_args[header->sub_codec].est_size (header->sub_codec, ctx->local.len);
    if (*compressed_len < min_required_compressed_len) {
        if (soft_fail) return false; // call me again with more memory
        ABORT ("%s: Compressing %s with %s need %u bytes, but allocated only %u", 
               VB_NAME, ctx->tag_name, codec_name(header->sub_codec), min_required_compressed_len, *compressed_len);
    }

    COPY_TIMER_COMPRESS (compressor_homp); // don't account for sub-codec compressor, it accounts for itself

    return compress (vb, ctx, header, 0, uncompressed_len, get_line_cb, STRa(compressed), false, name);
}

//--------------
// PIZ side
//--------------

CODEC_RECONSTRUCT (codec_homp_reconstruct)
{
    START_TIMER;
 
    rom seq = VB_DT(SAM) ? sam_piz_get_textual_seq(vb) : last_txtx (vb, CTX(FASTQ_SQBITMAP)); 

    // case: translating to FASTQ, and alignment is revcomped in SAM - use reversed SEQ as in SAM file
    bool out_to_fq_and_revcomp = (OUT_DT(FASTQ) && VB_DT(SAM) && sam_is_last_flags_rev_comp (vb));

    if (out_to_fq_and_revcomp) {
        buf_alloc (vb, &vb->scratch, 0, len, char, 0, "scratch");
        str_reverse (B1STc(vb->scratch), seq, len);
        seq = B1STc(vb->scratch);
    }

    // case: Deep, and len is only the trimmed suffix as the rest if copied from SAM (see fastq_special_deep_copy_QUAL)
    if (flag.deep && len < vb->seq_len)
        seq += (vb->seq_len - len); // advance seq to the trimmed part too
    else
        ASSPIZ (len == vb->seq_len, "expecting len=%u == vb->seq_len=%u", len, vb->seq_len);

    // reconsutrct by expanding homopolymers 
    char *next_recon = BAFTtxt;
    char *next_condensed = Bc(ctx->local, ctx->next_local);

    if (*next_condensed == ' ') { // SAM missing quality (expressed as a ' ')
        sam_reconstruct_missing_quality (vb, reconstruct);
        ctx->next_local++;
    }
    
    else {
        for (uint32_t i=0; i < len; i++) {
            unsigned hp_len = homopolymer_len (seq, len, i);  

            if (hp_len > 1) {
                // case: non-condensable
                if (*next_condensed & 0x80) {
                    *next_recon++ = *next_condensed++ & 0x7f;
                    for (unsigned h=1; h < hp_len; h++)
                        *next_recon++ = *next_condensed++;
                }

                else { 
                    char prev_qual = 0;
                    for (unsigned h=0; h < (hp_len + 1) / 2; h++) 
                        *next_recon++ = prev_qual = (prev_qual == TOP_QUAL) ? TOP_QUAL : (*next_condensed++);

                    rom c = next_recon - 1 - (hp_len&1);
                    for (unsigned h=0; h < hp_len / 2; h++)
                        *next_recon++ = *c--; // copy mirror
                }

                i += hp_len - 1; // skip to end of homopolymer
            }
            
            else
                *next_recon++ = *next_condensed++;
        }

        if (reconstruct) Ltxt += len;

        ctx->next_local = BNUM (ctx->local, next_condensed);
    }

    COPY_TIMER(codec_homp_reconstruct);
}
