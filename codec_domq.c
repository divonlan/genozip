// ------------------------------------------------------------------
//   codec_domq.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// compression algorithm for QUAL value that has a dominant value ("dom") (roughly over 50%) - as typically does binned Illumina
// data with a dominant 'F'. We use two contexts' local buffer:
// QUAL - will contain an array of all values EXCEPT for doms. It is assumed that before any non-dom value,
// including at the beginning of the data, there is a run of doms. If there isn't a run of doms at that point, 
// we insert a NO_DOMS. In addition, if there is a run of doms at the end of the data, there will be a terminating NO_DOMS in qual
// DOMQRUNS - for each dom run we have a value 0-254 that is interpreted as a dom run of 1-255 doms,
// or a value of 255 means a run of 255 and then continue the run with the next dom value, 
// thereby allowing runs of over 255 (e.g. a run "255-255-5" would be a run of 255+255+5=515 doms)

/* feasibility test on Illumina binned data:
-rw-r--r-- 1 USER 197610 14236798 Jul 25 16:49 qual
-rw-r--r-- 1 USER 197610   924542 Jul 26 18:17 qual.bz2 

-rw-r--r-- 1 USER 197610   295298 Jul 26 17:07 qual.chars.bz2
-rw-r--r-- 1 USER 197610   600946 Jul 26 17:07 qdomruns.bz2
                           896244 3% better and 10% faster

-rw-r--r-- 1 USER 197610   274560 Jul 26 17:11 qual.chars.xz
-rw-r--r-- 1 USER 197610   545944 Jul 26 17:11 qdomruns.xz
                           820504 12% better and 3.5X more time
*/

#include "vblock.h"
#include "data_types.h"
#include "piz.h"
#include "profiler.h"
#include "codec.h"

#define NO_DOMS '\x1'

//--------------
// ZIP side
//--------------

// Sample a few lines, and check that at least 50% of the Phred scores are a single character. 
// This is typically with Illumina binning and "normal" samples where most scores are F
// but might apply with other technologies too, including in combination with our optimize-QUAL
// Returns the character that appears more than 50% of the sample lines tested, or -1 if there isn't one.
bool codec_domq_comp_init (VBlock *vb, DidIType qual_did_i, LocalGetLineCB callback)
{
#   define DOMQUAL_THREADSHOLD_DOM_OF_TOTAL 0.5 // minimum % - doms of of total to trigger domqual
#   define DOMQUAL_THREADSHOLD_NUM_CHARS 5      // not worth it if less than this (and will fail in SAM with 1)
#   define DOMQUAL_LINE_SAMPLE_LEN 500          // we don't need more than this to find the dom (in case of long reads of 10s of thousands)
#   define NUM_LINES_IN_SAMPLE 5

    Context *qual_ctx = &vb->contexts[qual_did_i];
    qual_ctx->lcodec = CODEC_UNKNOWN; // cancel possible inheritence from previous VB

    uint32_t char_counter[256] = { 0 };
    uint32_t total_len = 0;
    for (uint32_t line_i=0; line_i < MIN (NUM_LINES_IN_SAMPLE, vb->lines.len); line_i++) {   
        char *qual_data;
        uint32_t qual_data_len;
        callback (vb, line_i, &qual_data, &qual_data_len, CALLBACK_NO_SIZE_LIMIT);
    
        if (qual_data_len > DOMQUAL_LINE_SAMPLE_LEN) qual_data_len = DOMQUAL_LINE_SAMPLE_LEN; 
    
        total_len += qual_data_len;

        for (unsigned j=0; j < qual_data_len; j++)
            char_counter[(uint8_t)qual_data[j]]++;
    }

    unsigned threshold = MAX ((unsigned)((double)total_len * DOMQUAL_THREADSHOLD_DOM_OF_TOTAL), DOMQUAL_THREADSHOLD_NUM_CHARS);

    for (unsigned c=33; c <= 126; c++)  // legal Phred scores only
        if (char_counter[c] > threshold) {
            qual_ctx->local.param = c;
            qual_ctx->local_param = true;
            qual_ctx->ltype   = LT_CODEC;
            qual_ctx->lcodec  = CODEC_DOMQ;

            Context *domqruns_ctx  = qual_ctx + 1;
            domqruns_ctx->ltype    = LT_UINT8;
            domqruns_ctx->st_did_i = qual_ctx->did_i;
            return true;
        }

    return false; // no value is dominant
}

static inline void codec_domq_add_runs (Buffer *qdomruns_buf, uint32_t runlen)
{
    // add one more bytes to represent the run
    while (runlen) {
        uint8_t subrun_len = (uint8_t)MIN (runlen, 254);

        NEXTENT (uint8_t, *qdomruns_buf) = (runlen <= 254 ? subrun_len : 255);
        runlen -= subrun_len;
    }
}

bool codec_domq_compress (VBlock *vb, 
                          SectionHeader *header,    // out
                          const char *uncompressed, // option 1 - not supported
                          uint32_t *uncompressed_len, 
                          LocalGetLineCB callback,  // option 2 - callback to fetch one line of qual data
                          char *compressed, uint32_t *compressed_len /* in/out */, 
                          bool soft_fail)
{
    START_TIMER;

    ASSERT0 (!uncompressed && callback, "Error in codec_domq_compress: only callback option is supported");

    SectionHeaderCtx *local_header = (SectionHeaderCtx *)header;
    Context *qual_ctx = ctx_get_existing_ctx (vb, local_header->dict_id);
    Context *qdomruns_ctx = qual_ctx + 1;

    const char dom = qual_ctx->local.param;
    ASSERT0 (dom, "Error in codec_domq_compress: dom is not set");

    Buffer *qual_buf     = &qual_ctx->local;
    Buffer *qdomruns_buf = &qdomruns_ctx->local;

    // this is usually enough, but might not be in some edge cases
    // note: qual_buf->len is the total length of all qual lines
    buf_alloc (vb, qual_buf, qual_buf->len / 5, 1, "context->local"); 
    qual_buf->param = dom; // dom goes into param, and eventually into SectionHeaderCtx.local_param

    buf_alloc (vb, qdomruns_buf, qual_buf->len / 10, 1, "context->local");

    qual_buf->len = 0; 
    uint32_t runlen = 0;
    
    for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) {   
        char *qual = 0;
        uint32_t qual_len = 0;
        callback (vb, line_i, &qual, &qual_len, CALLBACK_NO_SIZE_LIMIT);

        // grow if needed
        buf_alloc_more (vb, qual_buf, 2 * qual_len, 0, char, 1.5, 0); // theoretical worst case is 2 characters (added NO_DOMS) per each original character
        buf_alloc_more (vb, qdomruns_buf, qual_len, 0, uint8_t, 1.5, 0);

        if (!qual) continue;

        for (uint32_t i=0; i < qual_len; i++) {    
            if (qual[i] == dom) 
                runlen++;
            
            else {
                // this non-dom value terminates a run of doms
                if (runlen) {
                    codec_domq_add_runs (qdomruns_buf, runlen);
                    runlen = 0;
                }

                // this non-dom does not terminate a run of doms - add NO_DOMs to indicate the missing dom run
                else 
                    NEXTENT (char, *qual_buf) = NO_DOMS;

                // add the non-dom character
                NEXTENT (char, *qual_buf) = qual[i];
            }
        }
    }

    // case: we have a final dom run. note: we mark the terminating run, for example to avoid a situation
    // where QUAL is empty if qual is just one run. We use NO_DOMS rather than another marker, to avoid introducing
    // another letter into the compressed alphabet
    if (runlen) {
        codec_domq_add_runs (qdomruns_buf, runlen); // add final dom runs
        NEXTENT (char, *qual_buf) = NO_DOMS;
    }

    Codec sub_codec = codec_args[CODEC_DOMQ].sub_codec;
    CodecCompress *compress = codec_args[sub_codec].compress;
    uint32_t min_required_compressed_len = codec_args[sub_codec].est_size (sub_codec, qual_buf->len);

    COPY_TIMER (compressor_domq); // don't account for sub-codec compressor, it accounts for itself

    // case: all good - compress the QUAL context; the DOMQRUNS will be compressed after us, as its the subsequent context.
    if (*compressed_len >= min_required_compressed_len) {

        // since codecs were already assigned to contexts before compression of all contexts begun, but
        // we just created this context now, we assign a codec manually
        codec_assign_best_codec (vb, qdomruns_ctx, NULL, SEC_LOCAL);

        *uncompressed_len = (uint32_t)qual_buf->len;
        return compress (vb, header, qual_buf->data, uncompressed_len, NULL, compressed, compressed_len, soft_fail);
    }

    // case: our uncompressed length is too long vs the allocation of compressed (in a rare case that domqual enlengthen QUAL)
    // fallback on compressing the QUAL data using sub_codec directly (by which compressed_len was alloceted in comp_compress)
    else {
        ((SectionHeaderCtx *)header)->ltype = LT_SEQUENCE; // not LD_CODEC any more
        header->codec     = sub_codec;
        header->sub_codec = CODEC_UNKNOWN;
        header->flags.ctx.copy_param = 0; // cancel flag
        buf_free (qual_buf);
        buf_free (qdomruns_buf);
        return compress (vb, header, NULL, uncompressed_len, callback, compressed, compressed_len, soft_fail);
    }
}

//--------------
// PIZ side
//--------------

// shorten a run, including handling multi-bytes run - preparing the run length for the next line, 
// by deducting the amount that was consumed by this line
static inline uint32_t shorten_run (uint8_t *run, uint32_t old_num_bytes, uint32_t old_runlen, uint32_t dec)
{
    int32_t new_runlen = old_runlen - dec;
    ASSERT (new_runlen >= 0, "Error in shorten_run: new_runlen=%d is out of range", new_runlen);

    uint32_t new_num_bytes = MAX (1, ((uint32_t)new_runlen + 253) / 254); // roundup (if runlen=0, we still need 1 byte)

    // update run
    uint32_t increment = old_num_bytes - new_num_bytes;
    for (uint32_t i=increment; i < old_num_bytes; i++) { // strat the run bytes pushed forward (by increment), if we need less bytes
        run[i] = (new_runlen > 254 ? 255 : new_runlen);
        new_runlen -= 254;
    }

    return increment;
}

// reconstructed a run of the dominant character
static inline uint32_t codec_domq_reconstruct_dom_run (VBlockP vb, Context *domqruns_ctx, char dom, uint32_t max_len)
{
    ASSERT (domqruns_ctx->next_local < domqruns_ctx->local.len, "Error in codec_domq_reconstruct_dom_run: unexpectedly reached the end of vb->domqruns_ctx in vb_i=%u (first_line=%u len=%u)", 
            vb->vblock_i, vb->first_line, (uint32_t)vb->lines.len);

    // read the entire runlength (even bytes that are in excess of max_len)
    uint32_t runlen=0, num_bytes;
    uint8_t this_byte=255;
    uint32_t start_next_local = domqruns_ctx->next_local;
    for (num_bytes=0; this_byte == 255 && domqruns_ctx->next_local < domqruns_ctx->local.len; num_bytes++) {
        this_byte = NEXTLOCAL (uint8_t, domqruns_ctx); // might be 0 in beginning of line if previous line consumed entire run
        runlen += (this_byte < 255 ? this_byte : 254); // 0-254 is the length, 255 means length 254 continued in next byte
    }

    // case: a run spans multiple lines - take only what we need, and leave the rest for the next line
    // note: if we use max_len exactly, then we still leave a run of 0 length, so next line can start with a "run" as usual
    if (runlen >= max_len) { 
        uint32_t increment = shorten_run (ENT (uint8_t, domqruns_ctx->local, start_next_local), num_bytes, runlen, max_len);
        
        domqruns_ctx->next_local = start_next_local + increment; // unconsume this run as we will consume it again in the next line (but shorter)
        runlen = max_len;
    }

    memset (AFTERENT (char, vb->txt_data), dom, runlen);
    vb->txt_data.len += runlen;

    return runlen;
}

// Explanation of the reconstruction process of QUAL data compressed with the DOMQ codec:
// 1) The QUAL and DOMQ sections are decompressed normally using their sub_codecs
// 2) When reconstructing a QUAL field on a specific line, piz calls the LT_CODEC reconstructor for CODEC_DOMQUAL, 
//    codec_domq_reconstruct,which combines data from the local buffers of QUAL and DOMQRUNS to reconstruct the original QUAL field.
void codec_domq_reconstruct (VBlockP vb, Codec codec, ContextP qual_ctx)
{
    bool reconstruct = !piz_is_skip_section (vb, SEC_LOCAL, qual_ctx->dict_id);

    Context *domqruns_ctx = qual_ctx + 1;   // the qdomruns context is always one after qual context
    char dom = (char)qual_ctx->local.param; // passed from SectionHeaderCtx.local_param

    uint32_t qual_len=0;
    while (qual_len < vb->seq_len) {

        char c = NEXTLOCAL (char, qual_ctx);
        if (c != NO_DOMS) {
            qual_len += codec_domq_reconstruct_dom_run (vb, domqruns_ctx, dom, vb->seq_len - qual_len);

            // case: we're at an end of a line that ended with a run
            if (qual_len == vb->seq_len) {
                qual_ctx->next_local--; // unconsume c
                break;
            }
        }

        else if ((uint32_t)qual_ctx->local.len == qual_ctx->next_local) { // this is an final-run indicator
            qual_len += codec_domq_reconstruct_dom_run (vb, domqruns_ctx, dom, vb->seq_len - qual_len);
            qual_ctx->next_local--; // leave it unconsumed as it might be needed by the next lines
            break;
        }
        else
            c = NEXTLOCAL (char, qual_ctx);

        if (reconstruct)
            RECONSTRUCT1 (c==' ' ? '*' : c); // in SAM, sam_zip_qual re-wrote a '*' marking 'unavailable' as ' ' to avoid confusing with '*' as a valid quality score

        qual_len++;
    }

    ASSERT (qual_len == vb->seq_len, "Error in codec_domq_reconstruct: expecting qual_len(%u) == vb->seq_len(%u) in vb_i=%u (last_line=%u, num_lines=%u) line_i=%u", 
            qual_len, vb->seq_len, vb->vblock_i, vb->first_line + (uint32_t)vb->lines.len-1, (uint32_t)vb->lines.len, vb->line_i);   
}
