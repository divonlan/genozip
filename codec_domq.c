// ------------------------------------------------------------------
//   codec_domq.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

// compression algorithm for QUAL value that has a dominant value ("dom") (roughly over 50%) - as typically does binned Illumina
// data with a dominant 'F'. We use two contexts' local buffer:
// QUAL - will contain an array of all values EXCEPT for doms. It is assumed that before any non-dom value,
// including at the beginning of the data, there is a run of doms. If there isn't a run of doms at that point, 
// we insert a NO_DOMS. In addition, if there is a run of doms at the end of the data, there will be a terminating NO_DOMS in qual
// DOMQRUNS - for each dom run we have a value 0-254 that is interpreted as a dom run of 1-255 doms,
// or a value of 255 means a run of 255 and then continue the run with the next dom value, 
// thereby allowing runs of over 255 (e.g. a run "255-255-5" would be a run of 255+255+5=515 doms)

#include "vblock.h"
#include "data_types.h"
#include "piz.h"
#include "reconstruct.h"
#include "profiler.h"
#include "codec.h"
#include "context.h"
#include "stats.h"

#define NO_DOMS '\x1'

//--------------
// ZIP side
//--------------

// Sample a few lines, and check that at least 50% of the Phred scores are a single character. 
// This is typically with Illumina binning and "normal" samples where most scores are F
// but might apply with other technologies too, including in combination with our optimize-QUAL
// Returns the character that appears more than 50% of the sample lines tested, or -1 if there isn't one.
bool codec_domq_comp_init (VBlockP vb, DidIType qual_did_i, LocalGetLineCB callback)
{
#   define DOMQUAL_THREADSHOLD_DOM_OF_TOTAL 0.45 // minimum % - doms of of total to trigger domqual
#   define DOMQUAL_THREADSHOLD_NUM_CHARS 5       // not worth it if less than this (and will fail in SAM with 1)
#   define DOMQUAL_SAMPLE_LEN 2500               // we don't need more than this to find the dom 
#   define NUM_LINES_IN_SAMPLE 5

    Context *qual_ctx = CTX(qual_did_i);
    qual_ctx->lcodec = CODEC_UNKNOWN; // cancel possible inheritence from previous VB

    uint32_t char_counter[256] = { 0 };
    uint32_t total_len = 0;
    uint32_t num_sampled_lines = MIN_(NUM_LINES_IN_SAMPLE, vb->lines.len);
    uint32_t sampled_one_line = DOMQUAL_SAMPLE_LEN / MAX_(1,num_sampled_lines);

    for (uint32_t line_i=0; line_i < num_sampled_lines; line_i++) {   
        STRw(qual_data);
        callback (vb, line_i, pSTRa (qual_data), CALLBACK_NO_SIZE_LIMIT, NULL);
    
        if (qual_data_len > sampled_one_line) qual_data_len = sampled_one_line; 
    
        total_len += qual_data_len;

        for (unsigned j=0; j < qual_data_len; j++)
            char_counter[(uint8_t)qual_data[j]]++;
    }

    unsigned threshold = MAX_((unsigned)((float)total_len * DOMQUAL_THREADSHOLD_DOM_OF_TOTAL), DOMQUAL_THREADSHOLD_NUM_CHARS);

    for (unsigned c=33; c <= 126; c++) { // legal Phred scores only
        
        //if (char_counter[c]) printf ("'%c' : %u\n", c, char_counter[c]);
        
        if (char_counter[c] > threshold) {
            qual_ctx->local.param   = c;
            qual_ctx->local_param   = true;
            qual_ctx->ltype         = LT_CODEC;
            qual_ctx->lcodec        = CODEC_DOMQ;

            Context *domqruns_ctx   = qual_ctx + 1;
            domqruns_ctx->ltype     = LT_UINT8;
            domqruns_ctx->local_dep = DEP_L1;  // DOMQRUNS.local is created with QUAL.local is compressed
            stats_set_consolidation (vb, qual_ctx->did_i, 1, domqruns_ctx->did_i);
            return true;
        }
    }

    return false; // no value is dominant
}

static inline void codec_domq_add_runs (Buffer *qdomruns_buf, uint32_t runlen)
{
    // add one more bytes to represent the run
    while (runlen) {
        uint8_t subrun_len = (uint8_t)MIN_(runlen, 254);

        BNXT8 (*qdomruns_buf) = (runlen <= 254 ? subrun_len : 255);
        runlen -= subrun_len;
    }
}

bool codec_domq_compress (VBlockP vb, 
                          SectionHeader *header,    // out
                          rom uncompressed,         // option 1 - not supported
                          uint32_t *uncompressed_len, 
                          LocalGetLineCB callback,  // option 2 - callback to fetch one line of qual data
                          STRe (compressed),        // in/out 
                          bool soft_fail, rom name)
{
    START_TIMER;

    ASSERT0 (!uncompressed && callback, "only callback option is supported");

    SectionHeaderCtx *local_header = (SectionHeaderCtx *)header;
    Context *qual_ctx = ECTX (local_header->dict_id);
    Context *qdomruns_ctx = qual_ctx + 1;

    const char dom = qual_ctx->local.param;
    ASSERT (dom, "dom is not set for %s in vb=%u \"%s\"", qual_ctx->tag_name, vb->vblock_i, name);

    Buffer *qual_buf     = &qual_ctx->local;
    Buffer *qdomruns_buf = &qdomruns_ctx->local;

    // case: this is our second entry, after soft-failing. Just continue from where we stopped
    if (!soft_fail) goto do_compress;

    // this is usually enough, but might not be in some edge cases
    // note: qual_buf->len is the total length of all qual lines
    buf_alloc (vb, qual_buf, 0, qual_buf->len / 5, char, 1, "contexts->local"); 
    qual_buf->param = dom; // dom goes into param, and eventually into SectionHeaderCtx.local_param

    buf_alloc (vb, qdomruns_buf, 0, qual_buf->len / 10, char, 1, "contexts->local");

    qual_buf->len = 0; 
    uint32_t runlen = 0;
    
    for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) {   
        STRw0 (qual);
        callback (vb, line_i, pSTRa(qual), CALLBACK_NO_SIZE_LIMIT, NULL);

        // grow if needed
        buf_alloc (vb, qual_buf, 2 * qual_len, 0, char, 1.5, 0); // theoretical worst case is 2 characters (added NO_DOMS) per each original character
        buf_alloc (vb, qdomruns_buf, qual_len, 0, uint8_t, 1.5, 0);

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
                    BNXTc (*qual_buf) = NO_DOMS;

                // add the non-dom character
                BNXTc (*qual_buf) = qual[i];
            }
        }
    }

    // case: we have a final dom run. note: we mark the terminating run, for example to avoid a situation
    // where QUAL is empty if qual is just one run. We use NO_DOMS rather than another marker, to avoid introducing
    // another letter into the compressed alphabet
    if (runlen) {
        codec_domq_add_runs (qdomruns_buf, runlen); // add final dom runs
        BNXTc (*qual_buf) = NO_DOMS;
    }

    qual_ctx->lcodec = CODEC_UNKNOWN;
    header->sub_codec = codec_assign_best_codec (vb, qual_ctx, &qual_ctx->local, SEC_LOCAL); // provide BufferP to override callback
    if (header->sub_codec == CODEC_UNKNOWN) header->sub_codec = CODEC_NONE; // really small

do_compress: ({});
    CodecCompress *compress = codec_args[header->sub_codec].compress;
    *uncompressed_len = (uint32_t)qual_buf->len;

    // make sure we have enough memory
    uint32_t min_required_compressed_len = codec_args[header->sub_codec].est_size (header->sub_codec, qual_buf->len);
    if (*compressed_len < min_required_compressed_len) {
        if (soft_fail) return false; // call me again with more memory
        ABORT ("Compressing %s in vb_i=%u with %s need %u bytes, but allocated only %u", qual_ctx->tag_name, vb->vblock_i, codec_name(header->sub_codec), min_required_compressed_len, *compressed_len);
    }

    COPY_TIMER_COMPRESS (compressor_domq); // don't account for sub-codec compressor, it accounts for itself

    return compress (vb, header, qual_buf->data, uncompressed_len, NULL, compressed, compressed_len, false, name);
}

//--------------
// PIZ side
//--------------

// shorten a run, including handling multi-bytes run - preparing the run length for the next line, 
// by deducting the amount that was consumed by this line
static inline uint32_t shorten_run (uint8_t *run, uint32_t old_num_bytes, uint32_t old_runlen, uint32_t dec)
{
    int32_t new_runlen = old_runlen - dec;
    ASSERT (new_runlen >= 0, "new_runlen=%d is out of range", new_runlen);

    uint32_t new_num_bytes = MAX_(1, ((uint32_t)new_runlen + 253) / 254); // roundup (if runlen=0, we still need 1 byte)

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
    ASSERT (domqruns_ctx->next_local < domqruns_ctx->local.len, "unexpectedly reached the end of vb->domqruns_ctx in vb_i=%u (first_line=%"PRIu64" len=%u)", 
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
        uint32_t increment = shorten_run (B8 (domqruns_ctx->local, start_next_local), num_bytes, runlen, max_len);
        
        domqruns_ctx->next_local = start_next_local + increment; // unconsume this run as we will consume it again in the next line (but shorter)
        runlen = max_len;
    }

    memset (BAFTc (vb->txt_data), dom, runlen);
    vb->txt_data.len += runlen;

    return runlen;
}

// Explanation of the reconstruction process of QUAL data compressed with the DOMQ codec:
// 1) The QUAL and DOMQ sections are decompressed normally using their sub_codecs
// 2) When reconstructing a QUAL field on a specific line, piz calls the LT_CODEC reconstructor for CODEC_DOMQUAL, 
//    codec_domq_reconstruct,which combines data from the local buffers of QUAL and DOMQRUNS to reconstruct the original QUAL field.
void codec_domq_reconstruct (VBlockP vb, Codec codec, ContextP qual_ctx)
{
    ASSERTNOTEMPTY (qual_ctx->local);
    
    if (!qual_ctx->is_loaded) return;
    bool reconstruct = true;
    
    Context *domqruns_ctx = qual_ctx + 1;   // the qdomruns context is always one after qual context
    char dom = (char)qual_ctx->local.param; // passed from SectionHeaderCtx.local_param

    uint32_t qual_len=0;
    uint32_t expected_qual_len = vb->seq_len;

    while (qual_len < expected_qual_len) {

        char c = NEXTLOCAL (char, qual_ctx);
        if (c != NO_DOMS) {
            qual_len += codec_domq_reconstruct_dom_run (vb, domqruns_ctx, dom, expected_qual_len - qual_len);

            // case: we're at an end of a line that ended with a run
            if (qual_len == expected_qual_len) {
                qual_ctx->next_local--; // unconsume c
                break;
            }
        }

        else if ((uint32_t)qual_ctx->local.len == qual_ctx->next_local) { // this is an final-run indicator
            qual_len += codec_domq_reconstruct_dom_run (vb, domqruns_ctx, dom, expected_qual_len - qual_len);
            qual_ctx->next_local--; // leave it unconsumed as it might be needed by the next lines
            break;
        }
        else
            c = NEXTLOCAL (char, qual_ctx);

        // case: handle SAM missing quality (may be expressed as a ' ' or ASCII 127)
        if (c == ' ' || c == 127) {
            expected_qual_len = 1;
            sam_reconstruct_missing_quality (vb, c, reconstruct);
        }
        
        else if (reconstruct) 
            RECONSTRUCT1 (c); 

        qual_len++;
    }

    ASSERT (qual_len == expected_qual_len, "expecting qual_len(%u) == expected_qual_len(%u) in vb_i=%u (last_line=%"PRIu64", num_lines=%"PRIu64") line_i=%"PRIu64"", 
            qual_len, expected_qual_len, vb->vblock_i, vb->first_line + vb->lines.len-1, vb->lines.len, vb->line_i);   
}
