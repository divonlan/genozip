// ------------------------------------------------------------------
//   codec_homp.c
//   Copyright (C) 2020-2024 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "codec.h"
#include "reconstruct.h"
#include "profiler.h"
#include "compressor.h"

#define TOP_QUAL 'I'

//--------------
// ZIP side
//--------------

static bool codec_homp_qual_data_is_a_fit_for_homp (VBlockP vb, ContextP qual_ctx, LocalGetLineCB qual_callback)
{
    LocalGetLineCB *seq_callback = (VB_DT(FASTQ) ? fastq_zip_seq : sam_zip_seq);

    #define NUM_HPS_IN_SAMPLE 100  // test at least this number of HPs (or as many as are available)
    #define SUCCESS_CRITERION 0.85 // at least 85% of tested homopolymers are condensable

    uint32_t count[2] = {}; // [0]=count non-condensable homopolymers [1]=condensable

    for (LineIType line_i=0; line_i < vb->lines.len32 && (count[false] + count[true] < NUM_HPS_IN_SAMPLE); line_i++) {   
        STRw(qual);
        qual_callback (vb, qual_ctx, line_i, pSTRa (qual), CALLBACK_NO_SIZE_LIMIT, NULL);
        
        if (!qual_len) continue;

        STRw(seq);
        seq_callback (vb, NULL, line_i,  pSTRa(seq), CALLBACK_NO_SIZE_LIMIT, NULL);

        for (unsigned i=0; i < seq_len; i++) {
            unsigned hp_len = homopolymer_len (STRa(seq), i);

            if (hp_len > 1) {
                char prev_qual = 0;
                bool condensable = true; // optimistic

                for (unsigned h=0; h < (hp_len + 1) / 2; h++) {
                    char this_qual   = qual[i + h];
                    char mirror_qual = qual[i + hp_len - 1 - h];

                    if (this_qual != mirror_qual || (prev_qual == TOP_QUAL && this_qual != TOP_QUAL)) {
                        condensable = false;
                        break;
                    }

                    prev_qual = this_qual;
                }

                i += hp_len - 1; // skip to end of homopolymer
                count[condensable]++;
            }
        }
    }

    bool success = (count[true] || count[false]) && 
                   ((double)count[true] / (double)(count[true] + count[false])) > SUCCESS_CRITERION; 
    
    if (flag.show_qual) 
        printf ("HOMP: breakdown of first %u homopolymers in %10s: [condensable]=%u [non-condensable]=%u success=%s\n",
                 count[true] + count[false], VB_NAME, count[true], count[false], TF(success));

    return success;
}

bool codec_homp_comp_init (VBlockP vb, Did qual_did_i, LocalGetLineCB get_line_cb)
{
    ContextP qual_ctx = CTX(qual_did_i);

    if ((!TECH(ULTIMA) && !TECH(UNKNOWN) && !MP(ULTIMA)) || // either the TECH or the mapper are an indication of potentially Ultima data 
        qual_did_i != FASTQ_QUAL /* =SAM_QUAL */                         ||
        !codec_homp_qual_data_is_a_fit_for_homp (vb, qual_ctx, get_line_cb))
        return false;

    qual_ctx->ltype  = LT_CODEC;
    qual_ctx->lcodec = CODEC_HOMP;

    if (segconf.running && TECH(UNKNOWN)) 
        segconf.tech = TECH_ULTIMA; // if tech is unknown, given HOMP compatability, it is likely Ultima
        
    return true;
}

COMPRESS (codec_homp_compress)
{
    START_TIMER;

    // case: this is our second entry, after soft-failing. Just continue from where we stopped
    if (!soft_fail) goto do_compress;

    LocalGetLineCB *seq_callback = (VB_DT(FASTQ) ? fastq_zip_seq : sam_zip_seq);
    void (*update_qual_len)(VBlockP, uint32_t, uint32_t) = (VB_DT(FASTQ) ? fastq_update_qual_len : sam_update_qual_len);

    // first pass - condese qual scores of homopolymers by removing redundant scores
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {   
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
