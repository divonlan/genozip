// ------------------------------------------------------------------
//   codec_domq.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// compression algorithm for QUAL value that has a dominant value ("dom") (roughly over 50%) - as typically does binned Illumina
// data with a dominant 'F'. We use two contexts' local buffer:
// QUAL - will contain an array of all values EXCEPT for doms. It is assumed that before any non-dom value,
// including at the beginning of the data, there is a run of doms. If there isn't a run of doms at that point, 
// we insert a no_doms. In addition, if there is a run of doms at the end of the data, there will be a terminating no_doms in qual
// DOMQRUNS - for each dom run we have a value 0-254 that is interpreted as a dom run of 1-255 doms,
// or a value of 255 means a run of 255 and then continue the run with the next dom value, 
// thereby allowing runs of over 255 (e.g. a run "255-255-5" would be a run of 255+255+5=515 doms)

#include "vblock.h"
#include "reconstruct.h"
#include "codec.h"
#include "base64.h"
#include "seg.h"

typedef struct {
    uint8_t *qual;
    uint32_t qual_len; 
    uint8_t dom;     // dom of this line
    bool is_diverse; // not enough % of of this line is dom
    bool is_rev;
} QualLine;

#define ql_buf domqruns_ctx->qual_line // Seg: an entry of QualLine per line. We can use this buffer this is called from seg_finalize and local_hash is no longer needed
#define normalize_buf qual_ctx->normalize_buf

#define NO_DOMS_v13 '\x1'

#define FIRST_Q 32  // first valid quality (in SAM terms) (' '=32 if QUAL="*")
#define LAST_Q  126 // last printable ascii
#define NUM_Qs (LAST_Q-FIRST_Q+1)

#define declare_domq_contexts(ctx)                                \
             qual_ctx     __attribute__((unused)) = (ctx);        \
    ContextP domqruns_ctx __attribute__((unused)) = qual_ctx + 1; \
    ContextP qualmplx_ctx __attribute__((unused)) = qual_ctx + 2; \
    ContextP divrqual_ctx __attribute__((unused)) = qual_ctx + 3
 
static void show_denormalize (VBlockP vb, bytes denormalize, const uint32_t lines_with_dom[NUM_Qs], rom dom_to_ascii, uint8_t width, uint8_t num_doms)
{
    iprintf ("%s: Doms present and the normalized QUAL values in each dom:\n", VB_NAME);

    for (uint8_t dom=0; dom < num_doms; dom++) {
        if (dom_to_ascii)
            iprintf ("dom='%c' (lines=%u): ", dom_to_ascii[dom], lines_with_dom[dom_to_ascii[dom] - FIRST_Q]);
        else 
            iprintf ("dom_i=%u: ", dom);

        char d;
        for (uint8_t q=0; q < width && (d=denormalize[dom * width + q]); q++) 
            iprintf ("'%c' ", d);

        iprint0 ("\n");
    }
}

void codec_domq_update_qual_len (VBlockP vb, ContextP ctx, uint32_t line_i, uint32_t new_len) 
{ 
    ContextP declare_domq_contexts (ctx);
    B(QualLine, ql_buf, line_i)->qual_len = new_len; 
}

//--------------
// ZIP side
//--------------

// get line hisogram
#define line_histogram(qual, qual_len) \
    uint32_t line_ascii_histogram[256] = {}; /* 256 and not NUM_Qs so we can validate */ \
    for (uint32_t base_i=0; base_i < (qual_len); base_i++) \
        line_ascii_histogram[((uint8_t*)qual)[base_i]]++

static bool codec_domq_qual_data_is_a_fit_for_domq (VBlockP vb, ContextP qual_ctx, LocalGetLineCB get_line_cb)
{
    ASSERT (!qual_ctx->local.len32 || qual_ctx->local.data || get_line_cb, "%s: ctx=%s: since len=%u but data=NULL, expecting a callback, but there is none", 
            VB_NAME, qual_ctx->tag_name, qual_ctx->local.len32);

#   define DOMQUAL_SAMPLE_LEN 2500   // we don't need more than this to find out
#   define NUM_LINES_IN_SAMPLE 10
#   define MINIMUM_PERCENT_DOM_PER_LINE   50 // a lower bar for testing - we just need to see that this file is of the time of doms, even if these few reads would be "diversity"
#   define MINIMUM_PERCENT_LINES_WITH_DOM 50

    uint32_t num_sampled_lines = get_line_cb ? MIN_(NUM_LINES_IN_SAMPLE, vb->lines.len32) : 1;
    uint32_t sampled_one_line  = DOMQUAL_SAMPLE_LEN / MAX_(1,num_sampled_lines);
    uint32_t num_tested_lines=0, num_lines_with_dom=0;

    for (LineIType line_i=0; line_i < num_sampled_lines; line_i++) {   
        STRw(qual);
        
        if (get_line_cb)
            get_line_cb (vb, qual_ctx, line_i, pSTRa (qual), CALLBACK_NO_SIZE_LIMIT, NULL);
        else {
            qual = B1STc (qual_ctx->local);
            qual_len = qual_ctx->local.len32;
        }

        qual_len = MIN_(qual_len, sampled_one_line);

        if (!qual_len) { 
            if (get_line_cb && num_sampled_lines < vb->lines.len32) { // happens when eg depn line qual is copied from prim
                num_sampled_lines++;
                continue;
            }
            else
                break;
        }

        // check if this line as a score value that appears over 50%
        line_histogram (qual, qual_len);
        for (int ascii_i=FIRST_Q; ascii_i <= LAST_Q; ascii_i++)
            if (line_ascii_histogram[ascii_i] * 100 / qual_len > MINIMUM_PERCENT_DOM_PER_LINE) { 
                num_lines_with_dom++;
                break;
            }
        
        num_tested_lines++;
    }
    
    if (flag.show_codec)
        printf ("%s fit for DOMQ: num_lines_with_dom=%u (i.e. min %u%% dom) out of of num_tested_lines=%u\n", 
                qual_ctx->tag_name, num_lines_with_dom, MINIMUM_PERCENT_DOM_PER_LINE, num_tested_lines);

    return num_tested_lines && (num_lines_with_dom * 100 / num_tested_lines > MINIMUM_PERCENT_LINES_WITH_DOM);
}

typedef struct { uint8_t q; uint32_t count; } Mapping;
static DESCENDING_SORTER (mapping_sorter, Mapping, count)

static void codec_domq_calc_histogram (VBlockP vb, ContextP qual_ctx, ContextP domqruns_ctx, uint32_t histogram[NUM_Qs][NUM_Qs], uint32_t lines_with_dom[NUM_Qs], uint8_t *has_diverse)
{
#   define DOMQ_THREADHOLD 85 // minimum % of scores in a line required for domq
    
    *has_diverse = false;
    uint32_t count_dom_lines = 0;

    for_buf2 (QualLine, ql, line_i, ql_buf) {
        if (!ql->qual_len) continue;

        // get line hisogram
        line_histogram (ql->qual, ql->qual_len);

        // validate bases and get line dom
        uint32_t max_score_count=0;
        for (int ascii_i=0; ascii_i < 256; ascii_i++) {
            ASSERT ((ascii_i >= FIRST_Q && ascii_i <= LAST_Q) || !line_ascii_histogram[ascii_i],
                    "%s/%u: QUAL value=%u is out of range [%u, %u] for %s", 
                    VB_NAME, line_i, ascii_i, FIRST_Q, LAST_Q, qual_ctx->tag_name);

            if (line_ascii_histogram[ascii_i] >= max_score_count) { // if equal, the higher ascii_i is the dom
                max_score_count = line_ascii_histogram[ascii_i];
                ql->dom = ascii_i - FIRST_Q; // dom of the line
            }
        }

        uint32_t percent_dom = 100 * line_ascii_histogram[ql->dom + FIRST_Q] / ql->qual_len;
        
        if (!(ql->is_diverse = (percent_dom < DOMQ_THREADHOLD))) {
            lines_with_dom[ql->dom]++;
    
            // add line_ascii_histogram to histogram according to its dom
            for (uint8_t q=0; q < NUM_Qs; q++)
                histogram[ql->dom][q] += line_ascii_histogram[q + FIRST_Q];
        }
        else
            *has_diverse = true;

        if (flag.show_qual) 
            count_dom_lines += !ql->is_diverse;
    }

    if (flag.show_qual) 
        iprintf ("\n%s: %u/%u dominant reads (%u%%)\n", VB_NAME, count_dom_lines, vb->lines.len32, 100*count_dom_lines / vb->lines.len32);
}

static uint8_t inline codec_domq_compact_histogram (uint32_t histogram[NUM_Qs][NUM_Qs], uint32_t lines_with_dom[NUM_Qs], 
                                                    uint8_t score_i_to_dom_i[NUM_Qs], char dom_to_ascii[NUM_Qs])
{
    uint8_t num_doms=0;

    for (uint8_t q=0; q < NUM_Qs; q++)
        if (lines_with_dom[q]) {
            score_i_to_dom_i[q] = num_doms;
            dom_to_ascii[num_doms] = q + FIRST_Q;

            if (num_doms != q)
                memcpy (histogram[num_doms], histogram[q], NUM_Qs * sizeof (uint32_t));
            
            num_doms++;
        }

    return num_doms;
}

static uint8_t codec_domq_calc_norm_table (VBlockP vb, ContextP qual_ctx, ContextP domqruns_ctx,
                                           uint32_t histogram[NUM_Qs][NUM_Qs], uint8_t num_doms,
                                           uint8_t denormalize[][NUM_Qs])
{

    // initialize
    buf_free (normalize_buf); // we can free global_hash as this function is called from seg_finalize 
    ARRAY_alloc (uint8_t, normalize, num_doms * NUM_Qs, true, normalize_buf, vb, "contexts->normalize_buf");
    memset (denormalize, 0, num_doms * NUM_Qs);

    uint8_t num_norm_qs=0; // highest q_norm of any dom
    
    for (uint8_t dom=0; dom < num_doms; dom++) {

        Mapping mapping[NUM_Qs];        
        for (uint8_t q=0; q < NUM_Qs; q++) // initialize
            mapping[q] = (Mapping){ .q = q, .count = histogram[dom][q] };

        qsort (mapping, NUM_Qs, sizeof(Mapping), mapping_sorter); // sort by count

        // prepare normalization and de-normalization tables:
        // the most common q for this dom is normalized to q_norm=0, the second most common to 1 etc
        uint8_t q_norm = 0; 
        for (; q_norm < NUM_Qs && mapping[q_norm].count; q_norm++) {
            normalize[dom * NUM_Qs + mapping[q_norm].q] = q_norm;
            denormalize[dom][q_norm] = mapping[q_norm].q + FIRST_Q;
        }

        num_norm_qs = MAX_(num_norm_qs, q_norm);
    }

    // compact and seg de-normalizing table (to be transferred to PIZ in QUALNORM.local)
    // table has a line for each dom, a column for each normalized value, and cell contents are the ascii q
    uint8_t denorm[num_doms * num_norm_qs];
    qual_ctx->local_param = true; 
    qual_ctx->local.prm8[0] = num_norm_qs | 0x80; // send this to PIZ. note: 0x80 means a diverse read is marked by dom_i==255. In 14.0.0-14.0.4 this was 0, and diverse read marker was num_norm_qs

    for (uint8_t dom_i=0; dom_i < num_doms; dom_i++)
        for (uint8_t q_norm=0; q_norm < num_norm_qs; q_norm++)
            denorm[dom_i * num_norm_qs + q_norm] = denormalize[dom_i][q_norm];

    // seg the denormalization table into DOMQRUNS - it is normally the same for many VBs
    char denorm_snip[base64_size ((uint32_t)num_doms * (uint32_t)num_norm_qs)];
    unsigned denorm_snip_len = base64_encode (denorm, num_doms * num_norm_qs, denorm_snip);

    seg_by_ctx (vb, STRa(denorm_snip), domqruns_ctx, 0);

    return num_norm_qs;
}

// modifies QUAL data in place, to make it compress better
static uint8_t codec_domq_prepare_normalize (VBlockP vb, ContextP ctx, LocalGetLineCB get_line_cb, ContextP qual_ctx, ContextP domqruns_ctx)
{
    uint32_t histogram[NUM_Qs][NUM_Qs] = {}; // a histogram summarizing quality scores - multiplexed by dom
    uint32_t lines_with_dom[NUM_Qs] = {}; // entry q is true if there exists any line in the VB for which this q is dom

    // get quality scores 
    buf_free (ql_buf); // we can free local_hash as this function is called from seg_finalize 
    buf_alloc_exact (vb, ql_buf, get_line_cb ? vb->lines.len : 1, QualLine, "contexts->qual_line");

    if (get_line_cb)
        for_buf2 (QualLine, ql, line_i, ql_buf) 
            get_line_cb (vb, ctx, line_i, (char**)pSTRa(ql->qual), CALLBACK_NO_SIZE_LIMIT, &ql->is_rev);
    else
        *B1ST (QualLine, ql_buf) = (QualLine){ .qual = B1ST8(qual_ctx->local), .qual_len = qual_ctx->local.len32 };

    // verify validy of QUAL, get dom of each line, and populate histogram, lines_with_dom
    codec_domq_calc_histogram (vb, qual_ctx, domqruns_ctx, histogram, lines_with_dom, &qual_ctx->local.prm8[1]); // we store has_diverse is prm8[1] for use by codec_domq_compress

    // shrink histogram to only qual scores that are a dom of any line
    uint8_t score_i_to_dom_i[NUM_Qs] = {};
    char dom_to_ascii[NUM_Qs] = {};
    uint8_t num_doms = codec_domq_compact_histogram (histogram, lines_with_dom, score_i_to_dom_i, dom_to_ascii);

    // calculate normalization table from histogram
    uint8_t denormalize[num_doms][NUM_Qs]; // maps q to its normalized value - a mapping for each dom
    uint8_t num_norm_qs = codec_domq_calc_norm_table (vb, qual_ctx, domqruns_ctx, histogram, num_doms, denormalize);

    // normalize the dom of each line
    for_buf (QualLine, ql, ql_buf) 
        if (ql->qual_len && !ql->is_diverse)
            ql->dom = score_i_to_dom_i[ql->dom]; 

    if (flag.show_qual) 
        show_denormalize (vb, (uint8_t*)denormalize, lines_with_dom, dom_to_ascii, NUM_Qs, num_doms);

    return num_norm_qs; 
}

// Sample a few lines, and check that at least 50% of the Phred scores are a single character. 
// This is typically with Illumina binning and "normal" samples where most scores are F
// but might apply with other technologies too, including in combination with our optimize-QUAL
// Returns the character that appears more than 50% of the sample lines tested, or -1 if there isn't one.
bool codec_domq_comp_init (VBlockP vb, Did qual_did_i, LocalGetLineCB get_line_cb)
{
    ContextP declare_domq_contexts (CTX(qual_did_i));
 
    if (codec_domq_qual_data_is_a_fit_for_domq (vb, qual_ctx, get_line_cb)) {
        qual_ctx->ltype         = LT_CODEC;
        qual_ctx->lcodec        = CODEC_DOMQ;

        domqruns_ctx->ltype     = qualmplx_ctx->ltype     = divrqual_ctx->ltype     = LT_UINT8;
        domqruns_ctx->local_dep = qualmplx_ctx->local_dep = divrqual_ctx->local_dep = DEP_L1;  
        domqruns_ctx->no_stons  = true; // no singletons as we use local for the data

        // normalize quality scores in preparation for compression. note: we do it here in seg, as we need to seg the denormalization table
        codec_domq_prepare_normalize (vb, qual_ctx, get_line_cb, qual_ctx, domqruns_ctx);

        return true;
    }
    else {
        qual_ctx->lcodec = CODEC_UNKNOWN; // cancel possible inheritence from previous VB
        return false; // sampled VB qual scores not a good fit for domqual
    }
}

// normalize qual in lines, so that lines can have a different dom characters - doms are always normalized to 0
static inline void codec_domq_normalize_qual (ContextP qual_ctx, ContextP domqruns_ctx)
{
    ARRAY (uint8_t, normalize, normalize_buf);

    // PASS 2: normalize qual scores of the lines that have dom
    for_buf2 (QualLine, ql, line_i, ql_buf) {
        if (!ql->qual_len || ql->is_diverse) continue;

        // in case of QUAL="*" (missing qual) it is a static pointer in memory which we cannot change in place. We replace to a static pointer to '\0'.
        // (since qual_len=1, the single value definitely got mapped to 0)
        if (ql->qual_len==1 && ql->qual[0]==' ')
            ql->qual = (uint8_t *)"";

        else 
            for (uint32_t base_i=0; base_i < ql->qual_len; base_i++) 
                ql->qual[base_i] = normalize[ql->dom * NUM_Qs + ql->qual[base_i]-FIRST_Q];
    }
}

static inline void codec_domq_add_runs (BufferP qdomruns_buf, uint32_t runlen)
{
    // add one more bytes to represent the run
    while (runlen) {
        uint8_t subrun_len = (uint8_t)MIN_(runlen, 254);

        BNXT8 (*qdomruns_buf) = (runlen <= 254 ? subrun_len : 255);
        runlen -= subrun_len;
    }
}

COMPRESS (codec_domq_compress)
{
    START_TIMER;
    
    ContextP declare_domq_contexts (ctx);

    codec_domq_normalize_qual (qual_ctx, domqruns_ctx);

    BufferP qual_buf     = &qual_ctx->local;
    BufferP qdomruns_buf = &domqruns_ctx->local;

    // case: this is our second entry, after soft-failing. Just continue from where we stopped
    if (!soft_fail) goto do_compress;

    uint8_t no_doms = qual_ctx->local.prm8[0] & 0x7f; // set by codec_domq_prepare_normalize

    uint32_t vb_qual_len = qual_buf->len32;

    BufferP non_dom_buf;
    if (get_line_cb) {
        non_dom_buf = qual_buf;
        non_dom_buf->len = 0;
    }

    // case: qual data is in qual buf - we will construct non-dom in scratch and copy later
    else {
        ASSERTNOTINUSE (vb->scratch);
        non_dom_buf = &vb->scratch;
    }

    // this is usually enough, but might not be in some edge cases
    // note: vb_qual_len is the total length of all qual lines
    buf_alloc (vb, non_dom_buf,          0, 1 + vb_qual_len / 5,  char, 1, CTX_TAG_LOCAL); 
    buf_alloc (vb, qdomruns_buf,         0, 1 + vb_qual_len / 10, char, 1, CTX_TAG_LOCAL);
    buf_alloc (vb, &qualmplx_ctx->local, 0, 1 + vb->lines.len,      char, 0, CTX_TAG_LOCAL);
    
    if (qual_ctx->local.prm8[1]) // has_diverse
        buf_alloc (vb, &divrqual_ctx->local, 0, 1 + vb_qual_len / 5, char, 1, CTX_TAG_LOCAL);

    uint32_t runlen = 0;
    
    for_buf2 (QualLine, ql, line_i, ql_buf) {

        // case: qual line should not be compressed. Might happen eg in a SAM DEPN component - line segged against SA Group
        if (!ql->qual_len) continue; 

        // case: diverse read (note: criteria may be changed without impacting file format)
        if (ql->is_diverse) {
            buf_alloc (vb, &divrqual_ctx->local, ql->qual_len, 0, char, 1.5, CTX_TAG_LOCAL); // larger alloc than buf_add_more

            memcpy (BAFTc(divrqual_ctx->local), ql->qual, ql->qual_len);
            divrqual_ctx->local.len32 += ql->qual_len;

            BNXT8 (qualmplx_ctx->local) = 255; 
        }

        // case: quality string dominated by one particular quality score
        else {
            buf_alloc (vb, non_dom_buf, 2 * ql->qual_len + 1, 0, char, 1.5, 0); // theoretical worst case is 2 characters (added no_doms) per each original character + 1 for final dom run
            buf_alloc (vb, qdomruns_buf, 0, qdomruns_buf->len + ql->qual_len + runlen / 254 + 1, uint8_t, 1.5, 0);

            BNXTc (qualmplx_ctx->local) = ql->dom; // this is dom_i - i.e. the line in the denormalization table

            for (uint32_t i=0; i < ql->qual_len; i++) {    
                if (ql->qual[i] == 0) // dom
                    runlen++;
                
                else {
                    // this non-dom value terminates a run of doms
                    if (runlen) {
                        codec_domq_add_runs (qdomruns_buf, runlen);
                        runlen = 0;
                    }

                    // this non-dom does not terminate a run of doms - add NO_DOMs to indicate the missing dom run
                    else 
                        BNXTc (*non_dom_buf) = no_doms;

                    // add the non-dom character
                    BNXTc (*non_dom_buf) = ql->qual[i];
                }
            }
        }
    }

    // case: we have a final dom run. note: we mark the terminating run, for example to avoid a situation
    // where QUAL is empty if qual is just one run. We use no_doms rather than another marker, to avoid introducing
    // another letter into the compressed alphabet
    if (runlen &&
        // note: if qdomruns_buf->len32=0 - this final run is the only run - i.e. the qual data of this VB 
        // consists of initial non-doms followed by a single final run. We can therefore refrain from having a 
        // qdomruns section and figure it out in recon.
        (qdomruns_buf->len32 || runlen < BLST(QualLine, ql_buf)->qual_len)) {
        
        buf_alloc (vb, qdomruns_buf, runlen / 254 + 1, 0, uint8_t, 0, 0);
        codec_domq_add_runs (qdomruns_buf, runlen); // add final dom runs
        BNXTc (*non_dom_buf) = no_doms;
    }

    if (!get_line_cb) {
        buf_copy (vb, qual_buf, non_dom_buf, char, 0, 0, CTX_TAG_LOCAL);
        buf_free (vb->scratch);
    }

    qual_ctx->lcodec = CODEC_UNKNOWN;
    
    // all diverse - compress 1 byte in local anyway, just so codec_domq_reconstruct gets called
    if (!qual_ctx->local.len32) {
        BNXTc(qual_ctx->local) = 'X';
        header->sub_codec = CODEC_NONE;
    }

    else {
        header->sub_codec = codec_assign_best_codec (vb, qual_ctx, &qual_ctx->local, SEC_LOCAL); // provide BufferP to override callback
        if (header->sub_codec == CODEC_UNKNOWN) header->sub_codec = CODEC_NONE; // really small
    }

do_compress: ({});
    CodecCompress *compress = codec_args[header->sub_codec].compress;
    *uncompressed_len = qual_buf->len32;

    // make sure we have enough memory
    uint32_t min_required_compressed_len = codec_args[header->sub_codec].est_size (header->sub_codec, qual_buf->len);
    if (*compressed_len < min_required_compressed_len) {
        if (soft_fail) return false; // call me again with more memory
        ABORT ("%s: Compressing %s with %s need %u bytes, but allocated only %u", VB_NAME, qual_ctx->tag_name, codec_name(header->sub_codec), min_required_compressed_len, *compressed_len);
    }

    COPY_TIMER_COMPRESS (compressor_domq); // don't account for sub-codec compressor, it accounts for itself

    return compress (vb, qual_ctx, header, B1STc(*qual_buf), uncompressed_len, NULL, STRa(compressed), false, name);
}

//--------------
// PIZ side
//--------------

// shorten a run, including handling multi-bytes run - preparing the run length for the next line, 
// by deducting the amount that was consumed by this line
uint32_t codec_domq_shorten_run (uint8_t *run, uint32_t full_num_bytes, uint32_t full_runlen, uint32_t this_runlen)
{
    uint32_t next_runlen = full_runlen - this_runlen;

    uint32_t new_num_bytes = MAX_(1, (next_runlen + 253) / 254); // roundup (if runlen=0, we still need 1 byte)

    // update run - first bytes remain all 255 (each counting for 254 dom qual scores)
    uint32_t increment = full_num_bytes - new_num_bytes;
    
    if (next_runlen) {
        // last byte is the remainder 1-254 (example: runlen=254*3=762 -> 255,255,254)
        uint8_t mod = next_runlen % 254;
        run[full_num_bytes-1] = mod ? mod : 254;
    }
    else
        // next line will start with a run of length 0
        run[full_num_bytes-1] = 0;

    return increment;
}

// reconstructed a run of the dominant character
static inline uint32_t codec_domq_reconstruct_dom_run (VBlockP vb, ContextP domqruns_ctx, char dom, uint32_t max_len, ReconType reconstruct)
{
    START_TIMER;

    // note: absent this test, in case it would have failed, we would be in an infinite loop
    ASSPIZ0 (domqruns_ctx->next_local < domqruns_ctx->local.len32, "unexpectedly reached the end of vb->domqruns_ctx");

    // read the entire runlength (even bytes that are in excess of max_len)
    uint8_t *runs  = B8 (domqruns_ctx->local, domqruns_ctx->next_local);
    uint8_t *start = runs;

    while (*runs++ == 255); // advance runs to after the first non-255 (note: a run is always terminated by a non-255)

    // sanity - if overflowing, runs will terminate in the overflow fence of the buffer (since its not 255)
    ASSPIZ (runs <= BAFT8(domqruns_ctx->local), "%s.local exhausted: len=%u", domqruns_ctx->tag_name, domqruns_ctx->local.len32); 

    uint32_t num_bytes = runs - start;

    uint32_t runlen = (num_bytes - 1) * 254 + *(runs-1);

    // case: a run spans multiple lines - take only what we need, and leave the rest for the next line
    // note: if we use max_len exactly, then we still leave a run of 0 length, so next line can start with a "run" as usual
    if (runlen >= max_len) { 
        domqruns_ctx->next_local += codec_domq_shorten_run (start, num_bytes, runlen, max_len); // unconsume this run as we will consume it again in the next line (but shorter)
        runlen = max_len;
    }
    else
        domqruns_ctx->next_local += num_bytes;

    if (reconstruct) {
        memset (BAFTtxt, dom, runlen);
        Ltxt += runlen;
    }

    COPY_TIMER (codec_domq_reconstruct_dom_run);

    return runlen;
}

static inline void codec_domq_reconstruct_do_v13 (VBlockP vb, ContextP qual_ctx, ReconType reconstruct)
{
    ContextP domqruns_ctx = qual_ctx + 1;   // the qdomruns context is always one after qual context

    char dom = (char)qual_ctx->local.prm8[0]; // up to v13: passed from SectionHeaderCtx.local_param. starting v14: always 0

    uint32_t qual_len=0;
    uint32_t expected_qual_len = vb->seq_len;

    while (qual_len < expected_qual_len) {
        ASSERT (qual_ctx->next_local < qual_ctx->local.len32, "%s: QUAL.local exhausted prematurely for ctx=%s: next=%u len=%u", 
                VB_NAME, qual_ctx->tag_name, qual_ctx->next_local, qual_ctx->local.len32);

        char c = NEXTLOCAL (char, qual_ctx);
        if (c != NO_DOMS_v13) {
            qual_len += codec_domq_reconstruct_dom_run (vb, domqruns_ctx, dom, expected_qual_len - qual_len, reconstruct);

            // case: we're at an end of a line that ended with a run
            if (qual_len == expected_qual_len) {
                qual_ctx->next_local--; // unconsume c
                break;
            }
        }

        else if (qual_ctx->local.len32 == qual_ctx->next_local) { // this is an final-run indicator
            qual_len += codec_domq_reconstruct_dom_run (vb, domqruns_ctx, dom, expected_qual_len - qual_len, reconstruct);
            qual_ctx->next_local--; // leave it unconsumed as it might be needed by the next lines
            break;
        }
        
        else
            c = NEXTLOCAL (char, qual_ctx);

        // case: handle SAM missing quality (expressed as a ' ')
        if (c == ' ') {
            expected_qual_len = 1;
            sam_reconstruct_missing_quality (vb, reconstruct);
        }
        
        else if (reconstruct) 
            RECONSTRUCT1 (c); 

        qual_len++;
    }

    ASSPIZ (qual_len == expected_qual_len, "expecting qual_len(%u) == expected_qual_len(%u)", qual_len, expected_qual_len);   
}

// PIZ: returns de-normalization vector for dom_i (i.e. for current line)
bytes codec_domq_piz_get_denorm (VBlockP vb, ContextP domqruns_ctx, uint8_t dom_i, uint8_t num_norm_qs)
{
    // initialize, if not already initialized
    if (!domqruns_ctx->domq_denorm.len) {
        STR(snip);
        PEEK_SNIP (domqruns_ctx->did_i); // only one snip per VB - used for all lines

        buf_alloc_exact (vb, domqruns_ctx->domq_denorm, snip_len, uint8_t, "domq_denorm");
        domqruns_ctx->domq_denorm.len32 = base64_decode (snip, &snip_len, B1ST8(domqruns_ctx->domq_denorm));

        if (flag.show_qual)
            show_denormalize (vb, B1ST8(domqruns_ctx->domq_denorm), NULL, NULL, num_norm_qs, domqruns_ctx->domq_denorm.len32 / (uint32_t)num_norm_qs);
    }

    return B8(domqruns_ctx->domq_denorm, (uint32_t)dom_i * (uint32_t)num_norm_qs);    
}

// Explanation of the reconstruction process of QUAL data compressed with the DOMQ codec:
// 1) The QUAL and DOMQ sections are decompressed normally using their sub_codecs
// 2) When reconstructing a QUAL field on a specific line, piz calls the LT_CODEC reconstructor for CODEC_DOMQUAL, 
//    codec_domq_reconstruct, which combines data from the local buffers of QUAL and DOMQRUNS to reconstruct the original QUAL field.
// 3) Starting v14, the resulting data is de-normalized
static inline void codec_domq_reconstruct_do (VBlockP vb, ContextP qual_ctx, ContextP domqruns_ctx, 
                                              uint32_t len, uint8_t dom_i, uint8_t no_dom, ReconType reconstruct)
{
    bytes denormalize = codec_domq_piz_get_denorm (vb, domqruns_ctx, dom_i, no_dom);

    uint8_t dom = denormalize[0];

    bool missing_qual = (dom==' '); // this is QUAL="*"

    uint32_t qual_len=0;
    uint32_t expected_qual_len = missing_qual ? 1 : len; 

    // case: all non-diverse lines contain only dom (after possibly a few initial characters - see defect 2023-04-27). 
    // we identify this by domqruns_ctx.local being empty.
    if (!domqruns_ctx->local.len32) {
        if (reconstruct) {            
            // possibly an initial string of non-dom values, followed by only dom values (in the domq lines) until the end of the VB
            while (qual_len < expected_qual_len && qual_ctx->next_local < qual_ctx->local.len32 - 1) { // -1 is important, bc a VB with no non-doms, would have a single 'X' in qual_ctx
                uint8_t expecting_no_dom = *B8(qual_ctx->local, qual_ctx->next_local);
                uint8_t q_norm           = *B8(qual_ctx->local, qual_ctx->next_local + 1);
                
                ASSPIZ (expecting_no_dom == no_dom, "expecting non-dom qual_len=%u but found %u, qual_ctx->local.len=%u", 
                        qual_len, expecting_no_dom, qual_ctx->local.len32);
                RECONSTRUCT1 (denormalize[q_norm]);
                qual_ctx->next_local += 2;
                qual_len++;
            } 

            // remainer is dom
            memset (BAFTtxt, dom, expected_qual_len - qual_len);
            Ltxt += expected_qual_len - qual_len;
        }

        qual_len = expected_qual_len;
    }

    // normal case: reconstruct runs
    else while (qual_len < expected_qual_len) {
        uint8_t q_norm = NEXTLOCAL (uint8_t, qual_ctx);

        if (q_norm != no_dom) { 
            qual_len += codec_domq_reconstruct_dom_run (vb, domqruns_ctx, dom, expected_qual_len - qual_len, reconstruct);

            // case: we're at an end of a line that ended with a run
            if (qual_len == expected_qual_len) {
                qual_ctx->next_local--; // unconsume q_norm
                break;
            }
        }

        else if (qual_ctx->local.len32 == qual_ctx->next_local) { // this is an final-run indicator
            qual_len += codec_domq_reconstruct_dom_run (vb, domqruns_ctx, dom, expected_qual_len - qual_len, reconstruct);
            qual_ctx->next_local--; // leave it unconsumed as it might be needed by the next lines
            break;
        }
        
        else 
            q_norm = NEXTLOCAL (uint8_t, qual_ctx);

        if (reconstruct) 
            RECONSTRUCT1 (denormalize[q_norm]); 
        
        qual_len++;
    }

    if (missing_qual) { 
        if (reconstruct) Ltxt--; // undo
        sam_reconstruct_missing_quality (vb, reconstruct);
    }
    
    ASSPIZ (qual_len == expected_qual_len, "expecting qual_len=%u == expected_qual_len=%u in ctx=%s", 
            qual_len, expected_qual_len, qual_ctx->tag_name);   
}

CODEC_RECONSTRUCT (codec_domq_reconstruct)
{   
    START_TIMER;

    if (!ctx->is_loaded && !(ctx+1)->is_loaded && !(ctx+2)->is_loaded && !(ctx+3)->is_loaded) return;

    ContextP declare_domq_contexts (ctx);

    ReconType reconstruct = RECON_ON;

    // case: up to v13, all reads were compressed with the same dom (no multiplexing)
    if (!VER(14)) 
        codec_domq_reconstruct_do_v13 (vb, qual_ctx, reconstruct);
    
    else { // v14 and above
        uint8_t dom_i = (qualmplx_ctx->local.len32 == 1) ? *B1ST8 (qualmplx_ctx->local) // context was not comprssed line-by-line (single buffer) (or had just one line - same effect)
                                                         : NEXTLOCAL(uint8_t, qualmplx_ctx); 
        uint8_t num_norm_qs = qual_ctx->local.prm8[0] & 0x7f;

        // note: in v14.0.0-14.0.4 a diverse read was represented by dom_i=num_norm_qs - this was a flawed design, bc dom_i could legitimately be num_norm_qs or larger.
        bool is_diverse = (dom_i == ((qual_ctx->local.prm8[0] & 0x80/*true starting v14.0.5*/) ? 255 : num_norm_qs)); 

        // case: dom read
        if (!is_diverse)
            codec_domq_reconstruct_do (vb, qual_ctx, domqruns_ctx, len, dom_i, num_norm_qs, reconstruct);

        // case: diverse (non-dom) read.
        else if (reconstruct) 
            RECONSTRUCT_NEXT (divrqual_ctx, len);
    }

    COPY_TIMER(codec_domq_reconstruct);
}