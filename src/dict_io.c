// ------------------------------------------------------------------
//   dict_io.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include "zfile.h"
#include "dispatcher.h"
#include "piz.h"
#include "compressor.h"
#include "dict_io.h"
#include "base64.h"
#include "sam_friend.h"

// -------------------------------------
// ZIP: Assign codecs to dictionaries
// -------------------------------------

static ContextP next_ctx = 0;

static void dict_io_prepare_for_assign_codec (VBlockP vb)
{
    // get next non-too-small dict
    while (next_ctx->dcodec != CODEC_UNKNOWN && 
           next_ctx < ZCTX(z_file->num_contexts)) next_ctx++;

    if (next_ctx < ZCTX(z_file->num_contexts)) {
        vb->fragment_ctx = next_ctx;
        vb->dispatch = READY_TO_COMPUTE;
        next_ctx++;
    }
}

static void dict_io_assign_codec_one_dict (VBlockP vb)
{
    if (!vb->fragment_ctx->dcodec) // not small dict already assigned in dict_io_assign_codecs
        vb->fragment_ctx->dcodec = codec_assign_best_codec (vb, vb->fragment_ctx, NULL, SEC_DICT);
    
    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// called by main thread in zip_write_global_area
void dict_io_assign_codecs (void)
{
    START_TIMER;

    next_ctx = ZCTX(0);
    
    // handle some dictionaries here, so save on thread creation for trival dictionaries
    for_zctx {
        // assign CODEC_NONE to all the to-small-to-compress dictionaries, 
        if (zctx->dict.len < MIN_LEN_FOR_COMPRESSION) 
            zctx->dcodec = CODEC_NONE;

        // assign CODEC_ARTB to dictionaries under 1KB (unless --best)
        else if (!flag.best && zctx->dict.len < 1 KB) 
            zctx->dcodec = CODEC_ARTB;
    }

    dispatcher_fan_out_task ("assign_dict_codecs", NULL, 0, "Writing dictionaries...", true, false, false, 0, 20000, true,
                             dict_io_prepare_for_assign_codec, 
                             dict_io_assign_codec_one_dict, 
                             NO_CALLBACK);

    COPY_TIMER_EVB (dict_io_assign_codecs);
}

// -------------------------------------
// ZIP: Compress and output dictionaries
// -------------------------------------

static Context *frag_ctx;
static const CtxNode *frag_next_node;
static unsigned frag_size;

// compress the dictionary fragment - either an entire dict, or divide it to fragments if large to allow multi-threaded
// compression and decompression
static void dict_io_prepare_for_compress (VBlockP vb)
{
    while (frag_ctx < ZCTX(z_file->num_contexts)) {

        if (!frag_next_node) {
            if (!frag_ctx->nodes.len ||
                frag_ctx->please_remove_dict ||  // this context belongs to a method that lost the test and is not used in this file
                (frag_ctx->rm_dict_all_the_same && !frag_ctx->override_rm_dict_ats)) { // if only word in dict is a SNIP_LOOKUP, we don't need the dict, as absent a dict, lookup will happen (unless any VB objects) 
                frag_ctx++;
                continue; // unused context
            }

            frag_next_node = B1ST (const CtxNode, frag_ctx->nodes);

            ASSERT (frag_next_node->char_index + frag_next_node->snip_len <= frag_ctx->dict.len, 
                    "Corrupt nodes in ctx=%.8s did_i=%u", frag_ctx->tag_name, (int)(frag_ctx - z_file->contexts));

            ASSERT (frag_ctx->dict_id.num, "dict_id=0 for context \"%s\" did_i=%u", frag_ctx->tag_name, (unsigned)(frag_ctx - ZCTX(0)));

            // frag_size is set by default 1MB, but must be at least double the longest snip, or it will cause
            // mis-calculation of size_upper_bound in dict_io_read_one_vb
            frag_size = 1 MB;
            for (const CtxNode *node=B1ST (CtxNode, frag_ctx->nodes); node <= BLST (CtxNode, frag_ctx->nodes); node++)
                if (node->snip_len * 2 > frag_size)
                    frag_size = 2 * roundup2pow ((uint64_t)node->snip_len); // must be power of 2 for dict_io_read_one_vb

            if (flag.show_compress)
                iprintf ("DICT:  %s: len=%"PRIu64" words=%"PRIu64"\n", ctx_tag_name_ex (frag_ctx).s, frag_ctx->dict.len, frag_ctx->nodes.len);
        }

        vb->fragment_ctx   = frag_ctx;
        vb->fragment_start = Bc (frag_ctx->dict, frag_next_node->char_index);

        ctx_zip_z_data_exist (frag_ctx);

        while (frag_next_node < BAFT (CtxNode, frag_ctx->nodes) && 
               vb->fragment_len + frag_next_node->snip_len + 1 < frag_size) {

            vb->fragment_len += frag_next_node->snip_len + 1;
            vb->fragment_num_words++;
            frag_next_node++;
        }

        if (frag_next_node == BAFT (CtxNode, frag_ctx->nodes)) {
            frag_ctx++;
            frag_next_node = NULL;
            frag_size = 0;
        }

        if (vb->fragment_len) {
            vb->dispatch = READY_TO_COMPUTE;
            break;
        }
    }
}

static void dict_io_compress_one_fragment (VBlockP vb)
{
    START_TIMER;

    SectionHeaderDictionary header = (SectionHeaderDictionary){ 
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .section_type          = SEC_DICT,
        .dict_helper           = vb->fragment_ctx->dict_helper, // 15.0.42
        .data_uncompressed_len = BGEN32 (vb->fragment_len),
        .codec                 = vb->fragment_ctx->dcodec,
        .vblock_i              = BGEN32 (vb->vblock_i),
        .num_snips             = BGEN32 (vb->fragment_num_words),
        .dict_id               = vb->fragment_ctx->dict_id,
        .flags                 = { .dictionary = vb->fragment_ctx->dict_flags } // v15
    };

    if (flag.show_dict || dict_id_is_show (vb->fragment_ctx->dict_id)) 
        iprintf ("\nShowing dicts of %s (did=%u fragment_num_snips=%u fragment_vb=%u)\n", 
                 ctx_tag_name_ex (vb->fragment_ctx).s, vb->fragment_ctx->did_i, vb->fragment_num_words, vb->vblock_i);
    
    if (dict_id_is_show (vb->fragment_ctx->dict_id))
        dict_io_print (info_stream, vb->fragment_ctx->dict_id, vb->fragment_start, vb->fragment_len, true, true, true, false);

    if (flag.show_contigs && vb->fragment_ctx->did_i == CHROM)
        dict_io_print (info_stream, vb->fragment_ctx->dict_id, vb->fragment_start, vb->fragment_len, true, true, true, VB_DT(SAM) || VB_DT(BAM));

    if (flag.show_time) codec_show_time (vb, st_name (SEC_DICT), vb->fragment_ctx->tag_name, vb->fragment_ctx->dcodec);

    comp_compress (vb, vb->fragment_ctx, &vb->z_data, &header, vb->fragment_start, NO_CALLBACK, "SEC_DICT");

    COPY_TIMER (dict_io_compress_one_fragment)    

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// called by main thread in zip_write_global_area
void dict_io_compress_dictionaries (void)
{
    START_TIMER;

    frag_ctx = ZCTX(0);
    frag_next_node = NULL;
    frag_size = 0;

    dict_io_assign_codecs(); // assign codecs to all contexts' dicts

    dispatcher_fan_out_task ("compress_dicts", NULL, 0, "Writing dictionaries...", 
                             false, false, 
                             flag.show_dict || flag.show_one_dict, // force single thread if displaying 
                             0, 20000, true,
                             dict_io_prepare_for_compress, 
                             dict_io_compress_one_fragment, 
                             zfile_output_processed_vb);

    COPY_TIMER_EVB (dict_io_compress_dictionaries);
}

// -------------------------------------
// PIZ: Read and decompress dictionaries
// -------------------------------------
static Section dict_sec; 
static ContextP dict_ctx;

// V8 files: calculate dict size by adding up the size of all fragments. 
// Note: We need this for supporting v8 reference files.
static uint64_t v8_get_dict_size (VBlockP vb, Section dict_sec)
{
    uint64_t size = 0;

    for (Section sec=dict_sec; sec->dict_id.num == dict_sec->dict_id.num ; sec++) 
        size += BGEN32 (zfile_read_section_header (vb, sec, SEC_DICT).common.data_uncompressed_len);

    return size;
}

static void dict_io_read_one_vb (VBlockP vb) // one vb = one fragment
{
    Section old_dict_sec = dict_sec;
    if (!sections_next_sec (&dict_sec, SEC_DICT)) 
        return; // we're done - no more SEC_DICT sections

    // while we could easily read non-consecutive sections, this is not expected to happen and may indicate multiple calls to zip_write_global_area
    ASSERT0 (!old_dict_sec || (old_dict_sec+1 == dict_sec), "Unexpectedly, not all SEC_DICT sections are consecutive in the Genozip file");

    bool new_ctx = (!dict_ctx || dict_sec->dict_id.num != dict_ctx->dict_id.num);
    if (new_ctx)
        dict_ctx = ctx_get_ctx_do (z_file->contexts, z_file->data_type, z_file->d2d_map, &z_file->num_contexts, dict_sec->dict_id, 0, 0);

    // note: skipping a section should be mirrored in a container filter
    if (piz_is_skip_section (z_file->data_type, SEC_DICT, COMP_NONE, dict_sec->dict_id, dict_sec->flags.flags, SKIP_PURPOSE_RECON)) {
        if (flag.debug_read_ctxs)
            iprintf ("%c Skipped loading DICT/%u %s.dict\n", sections_read_prefix (flag.preprocessing), vb->vblock_i, dict_ctx->tag_name);
        goto done;
    }

    dict_ctx->dict_flags  = dict_sec->flags.dictionary; // v15
    
    int32_t offset = zfile_read_section (z_file, vb, dict_sec->vblock_i, &vb->z_data, "z_data", SEC_DICT, dict_sec);    
    SectionHeaderDictionaryP header = 
        (offset != SECTION_SKIPPED) ? (SectionHeaderDictionaryP)vb->z_data.data : NULL;

    dict_ctx->is_loaded = (offset != SECTION_SKIPPED);

    vb->fragment_len = header ? BGEN32 (header->data_uncompressed_len) : 0;

    ASSERT (!header || header->dict_id.num == dict_sec->dict_id.num, "Expecting dictionary fragment with DictId=%s but found one with DictId=%s",
            dis_dict_id (dict_sec->dict_id).s, dis_dict_id (header->dict_id).s);

    // new context
    // same-dict fragments are consecutive in the file, and all but the last are frag_size or a bit less, allowing pre-allocation
    if (header && new_ctx) {
        unsigned num_fragments=0; 
        for (Section sec=dict_sec; sec->dict_id.num == dict_ctx->dict_id.num; sec++) num_fragments++;

        // get size: for multi-fragment dictionaries, first fragment will be at or less than a power of 2, but more than the previous power of two.
        // this allows us to calculate the frag_size with which this dictionary was compressed and hence an upper bound on the size
        uint64_t size_upper_bound = (num_fragments == 1) ? vb->fragment_len : ((uint64_t)roundup2pow (vb->fragment_len) * (uint64_t)num_fragments);
        
        buf_alloc (evb, &dict_ctx->dict, 0, VER(9) ? size_upper_bound : v8_get_dict_size (vb, dict_sec), char, 0, "zctx->dict");
        dict_ctx->dict.prm32[0] = num_fragments;     // for error reporting
        dict_ctx->dict.prm32[1] = vb->fragment_len;

        buf_set_shared (&dict_ctx->dict);
    }

    if (header) {
        vb->fragment_ctx         = dict_ctx;
        vb->fragment_start       = Bc(dict_ctx->dict, dict_ctx->dict.len);
        dict_ctx->word_list.len += header ? BGEN32 (header->num_snips) : 0;
        dict_ctx->dict.len      += (uint64_t)vb->fragment_len;
        dict_ctx->dict_helper    = header->dict_helper; // 15.0.42 (union with con_rep_special)

        ASSERT (dict_ctx->dict.len <= dict_ctx->dict.size, "%s: Dict %s vb=%u len=%"PRIu64" exceeds allocated size=%"PRIu64" (num_fragments=%u fragment_1_len=%u)", 
                z_name, dict_ctx->tag_name, vb->vblock_i, dict_ctx->dict.len, (uint64_t)dict_ctx->dict.size, dict_ctx->dict.prm32[0], dict_ctx->dict.prm32[1]);

        if (flag.debug_read_ctxs)
            sections_show_header ((SectionHeaderP)header, NULL, COMP_NONE, dict_sec->offset, sections_read_prefix (flag.preprocessing));
    }

done: 
    // note: in cases we just "goto" here, no data is read, and a thread is needlessly created to decompress it
    // this is because the vb_i of the section needs to match the vb_i of the thread
    vb->dispatch = READY_TO_COMPUTE;
}

// entry point of compute thread of dictionary decompression of one fragment
static void dict_io_uncompress_one_vb (VBlockP vb)
{
    if (!vb->fragment_ctx || flag.only_headers) goto done; // nothing to do in this thread
    SectionHeaderDictionaryP header = (SectionHeaderDictionaryP)vb->z_data.data;

    ASSERT (vb->fragment_start + BGEN32 (header->data_uncompressed_len) <= BAFTc (vb->fragment_ctx->dict), 
            "Buffer overflow when uncompressing dict=%s", vb->fragment_ctx->tag_name);

    // uncompress to a location within the dict buffer - while multiple threads are uncompressing into 
    // non-overlappying regions in the same buffer in parallel
    buf_overlay_partial (vb, &vb->scratch, &vb->fragment_ctx->dict, BNUM64(vb->fragment_ctx->dict, vb->fragment_start), "scratch");
    zfile_uncompress_section (vb, header, &vb->scratch, NULL, 0, SEC_DICT); // NULL name prevents buf_alloc
    buf_destroy (vb->scratch);

done:
    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// PIZ
static void dict_io_dict_build_word_list_one (ContextP zctx)
{
    if (!zctx->word_list.len || zctx->word_list.data) return; // skip if 1. no words, or 2. already built

    buf_alloc (evb, &zctx->word_list, 0, zctx->word_list.len, CtxWord, 0, "zctx->word_list");
    buf_set_shared (&zctx->word_list);

    rom word_start = zctx->dict.data;
    for (uint64_t snip_i=0; snip_i < zctx->word_list.len; snip_i++) {

        rom c=word_start; while (*c) c++;

        uint64_t index = BNUM64 (zctx->dict, word_start); 
        uint64_t len   = c - word_start;

        if (index > CTX_MAX_DICT_LEN || len > CTX_MAX_SNIP_LEN) {
            ASSERT (!VER(14), "A word was found in zctx=%s with index=%"PRIu64" amd len=%"PRIu64". This index/len is beyond current limits of Genozip. Use Genozip v13 to decompress this file.",
                    zctx->tag_name, (uint64_t)index, (uint64_t)len);

            ABORT ("A word was found in zctx=%s with index=%"PRIu64" amd len=%"PRIu64", which are beyond the limits",
                    zctx->tag_name, (uint64_t)index, (uint64_t)len);
        }       

        *B(CtxWord, zctx->word_list, snip_i) = (CtxWord){ .char_index = index, .snip_len = len };

        word_start = c+1; // skip over the \0 separator
    }

    ASSERT (word_start == BAFTc(zctx->dict), "When building %s.word_list: expected to consume dict.len=%"PRIu64" bytes (in %"PRIu64" words), but consumed %"PRIu64" bytes",
            zctx->tag_name, zctx->dict.len, zctx->word_list.len, BNUM64(zctx->dict, word_start));
}

static void dict_io_build_word_lists (void)
{    
    START_TIMER;

    if (flag.show_headers && is_genocat)
        dict_io_dict_build_word_list_one (ZCTX(CHROM)); // only CHROM is needed
    
    else for_zctx
        dict_io_dict_build_word_list_one (zctx);

    COPY_TIMER_EVB (dict_io_build_word_lists);
}

// PIZ main thread
void dict_io_read_all_dictionaries (void)
{
    START_TIMER;

    dict_sec = NULL;
    dict_ctx = NULL;

    dispatcher_fan_out_task (flag.reading_reference ? "read_dicts_ref" : "read_dicts",
                             NULL, 0, 0, 
                             true, flag.test,
                             false,
                             0, 
                             10, // must be short - many dictionaries are just a few bytes
                             true,
                             dict_io_read_one_vb, 
                             dict_io_uncompress_one_vb,
                             NO_CALLBACK);

    // build word lists in z_file->contexts with dictionary data 
    if (!flag.only_headers)
        dict_io_build_word_lists();

    // output the dictionaries if we're asked to
    if (flag.show_dict || flag.show_one_dict || flag.show_contigs) {
        for (uint32_t did_i=0; did_i < z_file->num_contexts; did_i++) {
            ContextP ctx = ZCTX(did_i);
            if (!ctx->dict.len) continue;

            if (flag.show_contigs && ctx->did_i == CHROM)
                dict_io_print (info_stream, ctx->dict_id, STRb(ctx->dict), true, false, true, Z_DT(SAM));
            
            if (flag.show_dict || (flag.show_one_dict && dict_id_is_show (ctx->dict_id))) 
                iprintf ("%-8s\tdid_i=%-3u\tnum_snips=%u\tdict_size=%s\tall_the_same_wi=%u\tdeep_sam=%s\tdeep_fastq=%s\n", 
                         ctx_tag_name_ex (ctx).s, did_i, ctx->word_list.len32, str_size (ctx->dict.len).s,
                         ctx->dict_flags.all_the_same_wi, TF(ctx->dict_flags.deep_sam), TF(ctx->dict_flags.deep_fastq));

            if (dict_id_is_show (ctx->dict_id))
                dict_io_print (info_stream, ctx->dict_id, STRb(ctx->dict), true, true, true, false);
        }
        
        if (flag.show_dict) iprint_newline(); // but not for show_contigs or show_one_dict, as often wc is used to count them 

        if (is_genocat) exit_ok; // if this is genocat - we're done
    }

    COPY_TIMER_EVB (dict_io_read_all_dictionaries);
}

StrTextMegaLong str_snip_ex (DataType dt, STRp(snip), bool add_quote)
{
    StrTextMegaLong s;
    int s_len=0;

    if (dt == DT_NONE)
        dt = txt_file ? txt_file->data_type : z_file->data_type;

    char op = (snip_len && snip[0] > 0 && snip[0] < 32) ? snip[0] : 0;
    int i=1;

    static rom special_names[NUM_DATATYPES][MAX_NUM_SPECIAL] = 
        { [DT_VCF]=VCF_SPECIAL_NAMES,       [DT_SAM]=SAM_SPECIAL_NAMES,     [DT_BAM]=SAM_SPECIAL_NAMES,
          [DT_FASTQ]=FASTQ_SPECIAL_NAMES,   [DT_FASTA]=FASTA_SPECIAL_NAMES, [DT_GFF]=GFF_SPECIAL_NAMES, 
          [DT_GNRIC]=GENERIC_SPECIAL_NAMES, [DT_LOCS]=LOCS_SPECIAL_NAMES,   [DT_BED]=BED_SPECIAL_NAMES };

    switch (op) {
        case 0                         : i--;                           break;
        case SNIP_LOOKUP               : SNPRINTF0(s, "[LOOKUP]");      break; 
        case SNIP_OTHER_LOOKUP         : SNPRINTF0(s, "[OLOOKUP]");     break;
        case v13_SNIP_MATE_LOOKUP      : SNPRINTF0(s, "[MLOOKUP]");     break;
        case SNIP_CONTAINER            : SNPRINTF0(s, "[CONTAINER]");   break;
        case SNIP_SELF_DELTA           : SNPRINTF0(s, "[DELTA]");       break;
        case SNIP_OTHER_DELTA          : SNPRINTF0(s, "[ODELTA]");      break; 
        case v13_SNIP_FASTQ_PAIR2_GPOS : SNPRINTF0(s, "[PAIR2GPOS]");   break; 
        case SNIP_REDIRECTION          : SNPRINTF0(s, "[REDIRECTION]"); break;
        case SNIP_DONT_STORE           : SNPRINTF0(s, "[DONT_STORE]");  break; 
        case SNIP_COPY                 : SNPRINTF0(s, "[COPY]");        break; 
        case dvcf_SNIP_DUAL            : SNPRINTF0(s, "[DUAL]");        break; 
        case SNIP_LOOKBACK             : SNPRINTF0(s, "[LOOKBACK]");    break;
        case v13_SNIP_COPY_BUDDY       : SNPRINTF0(s, "[BCOPY]");       break; 
        case SNIP_DIFF                 : SNPRINTF0(s, "[DIFF]");        break; 
        case SNIP_NUMERIC              : SNPRINTF0(s, "[NUMERIC]");     break; 
        case SNIP_SPECIAL              : if (z_file && special_names[dt][snip[1]-32]) 
                                            SNPRINTF (s, "[%s_SPECIAL_%s]", dt_name(dt), special_names[dt][snip[1]-32]);
                                         else 
                                            SNPRINTF (s, "[SPECIAL-%u]", snip[1]-32); 
                                         i++; 
                                         break;
        default                        : SNPRINTF (s, "\\x%x", (uint8_t)op);
    }

    #define X_SAM(sp)   (op == SNIP_SPECIAL && (dt==DT_SAM || dt==DT_BAM) && (snip[1] == SAM_SPECIAL_##sp))
    #define X_VCF(sp)   (op == SNIP_SPECIAL && dt==DT_VCF   && (snip[1] == VCF_SPECIAL_##sp))
    #define X_GFF(sp)   (op == SNIP_SPECIAL && dt==DT_GFF   && (snip[1] == GFF_SPECIAL_##sp))
    #define X_FASTQ(sp) (op == SNIP_SPECIAL && dt==DT_FASTQ && (snip[1] == FASTQ_SPECIAL_##sp))

    if (op == SNIP_OTHER_LOOKUP || op == SNIP_OTHER_DELTA || op == SNIP_COPY || op == SNIP_REDIRECTION ||
        X_VCF(LEN_OF) || X_VCF(ARRAY_LEN_OF) || X_VCF(COPY_MATE) || (X_VCF(GQ) && snip_len > 8) ||
        (X_SAM(COPY_BUDDY) && snip_len > 3))
    {
        unsigned b64_len = base64_sizeof (DictId);
        DictId dict_id;
        base64_decode (&snip[1 + (op==SNIP_SPECIAL)], &b64_len, dict_id.id, sizeof(DictId));
        SNPRINTF (s, "(%s)", dis_dict_id (dict_id).s);
        i += b64_len;
    }

    else if (op == SNIP_CONTAINER) {
        // decode
        Container con;
        unsigned b64_len = snip_len - 1; // maximum length of b64 
        base64_decode (&snip[1], &b64_len, (uint8_t*)&con, -1);
        i += b64_len;

        con.repeats = BGEN24 (con.repeats);
        unsigned prefixes_len = snip_len - i;
        SNPRINTF (s, "%.*s", (int)sizeof(s)-20, container_to_json (&con, &snip[i], prefixes_len).s);
        i += prefixes_len;
    }

    // case: MINUS
    else if (X_VCF(MINUS) || X_GFF(MINUS)) {
        DictId dicts[2]; 
        unsigned b64_len = snip_len - 2; 
        base64_decode (&snip[2], &b64_len, (uint8_t *)dicts, 2 * sizeof(DictId));

        i += b64_len;

        SNPRINTF (s, " (%s - %s)", dis_dict_id (dicts[0]).s, dis_dict_id (dicts[1]).s); 
    }
        
    // case: PLUS
    else if (X_VCF(PLUS) || X_SAM(PLUS)) {
        DictId dicts[MAX_SNIP_DICTS]; 
        unsigned b64_len = snip_len - 2; 
        int num_dicts = base64_decode (&snip[2], &b64_len, (uint8_t *)dicts, -1) / sizeof (DictId);

        SNPRINTF0 (s, " (");

        for (int i=0; i < num_dicts; i++) {
            SNPRINTF (s, "%s", dis_dict_id (dicts[i]).s);
            SNPRINTF0 (s, (i < num_dicts-1) ? " + " : ")");
        }

        i += b64_len;
    }

    else if (X_VCF(DEFER) && snip_len > 2) {
        SNPRINTF (s, "%.*s", (int)sizeof (s)-50, str_snip_ex (dt, snip+2, snip_len-2, false).s); // recursive call
        return s;
    }  

    else if ((X_SAM(SQUANK)) && snip_len > 2) {
        SNPRINTF (s, "[%s]", (rom[])SQUANK_NAMES[snip[2] - '0']);
        return s;
    }  

    else if (X_SAM(CIGAR) && snip[2] < 0) {
        SNPRINTF (s, "[%s]%s", (rom[])CIGAR_OP_NAMES[(uint8_t)snip[2] - 0x80], 
                  snip[2] == SQUANK ? (rom[])SQUANK_NAMES[snip[3] - '0'] : &snip[3]);
        return s;
    }  

    else if (X_SAM(SEQ)) {
        SNPRINTF (s, "{%s%s%s%s%.12s%s }", 
                  (snip_len >= 3 && snip[2]=='1') ? " FORCE_SEQ"    : "",
                  (snip_len >= 4 && snip[3]=='1') ? " ALIGNER"      : "",
                  (snip_len >= 5 && snip[4]=='1') ? " PERFECT"      : "",
                  (snip_len >= 6 && snip[5]=='1') ? " NO_ANAL_DEPN" : "",
                  cond_str (snip_len >= 7 && snip[6]!='0', " BISULFITE=", &snip[6]),
                  (snip_len >= 8 && snip[7]=='1') ? " FRC_VERBATIM" : "");
        return s;
    }
    else if (X_FASTQ(SEQ_by_bamass)) {
        SNPRINTF (s, "[%s]", snip[2]=='0' ? "PERFECT" : "BITMAP");
        return s;
    }

    // case: SPECIAL with base64-encoded dictionaries: 
    // last item in a \t-separated array should be eg "MUX0"SNIP_SPECIAL
    // SNIP_SPECIAL tell us to print the dictionaries. 0 is the offset from the beginning.
    else if (op == SNIP_SPECIAL && snip[snip_len-1] == SNIP_SPECIAL) {
        unsigned start_at = i + (snip[snip_len-2] - '0');
        str_split (snip + start_at, snip_len - start_at, 0, '\t', item, false);

        DictId dict_id_1st, dict_id_lst;
        unsigned b64_len = item_lens[0]; 
        base64_decode (items[0], &b64_len, dict_id_1st.id, sizeof(DictId));

        if (n_items > 2) {
            b64_len = item_lens[1]; 
            base64_decode (items[n_items-2], &b64_len, dict_id_lst.id, sizeof(DictId));
        }

        SNPRINTF (s, " %.*s∈{%s%s}", STRfi (item, n_items-1), dis_dict_id (dict_id_1st).s, 
                  cond_str (n_items > 2, " → ", dis_dict_id (dict_id_lst).s));

        snip_len = start_at; // print the bytes before the dictionaries, if any
    }

    if (add_quote && !op) SNPRINTF0 (s, "\"");    

    int len = MAX_(0, MIN_(snip_len - i, sizeof(s) - add_quote*2 - 1));
    s_len += str_to_printable (snip+i, len, &s.s[s_len], sizeof (s.s) - s_len - add_quote - 1); // also nul-terminates

    if (add_quote && !op) SNPRINTF0 (s, "\"");

    return s;
}

// decide for a SPECIAL snip in a Deep file, whether it is a FASTQ or SAM special (ugly code)
static DataType dict_io_print_deep_dt_by_special (char special, DictId dict_id)
{
    #if NUM_FASTQ_SPECIAL != 16
    #error need to update dict_io_print_deep_dt_by_special()
    #endif

    #define SP(sam_name, fq_name) ({ \
        ASSERT ((int)SAM_SPECIAL_##sam_name == (int)FASTQ_SPECIAL_##fq_name, "expecting SAM_SPECIAL_%s == FASTQ_SPECIAL_%s", #sam_name, #fq_name); /* optimized away if ok */ \
        special == SAM_SPECIAL_##sam_name; \
    })

    //     SAM                    FASTQ               // special index as appears in sam.h and fastq.h
    //     ---------------------  --------------      ------------------------------------------------
    
    // SAM specials that don't appear in v15 (when deep was introduced), so these specials are definitely FASTQ
    if (SP(TLEN_old,              PAIR2_GPOS)      || // 1
        SP(MD_old,                deep_copy_QNAME) || // 4
        SP(PNEXT_IS_PREV_POS_old, ULTIMA_C)        || // 10
        SP(COPY_MATE_TLEN_old,    AGENT_QX))          // 12
        return DT_FASTQ;

    // FASTQ specials that can't appear in deep files - so these are definitely SAM
    if (SP(TLEN,                  SEQ_by_bamass))     // 15
        return DT_SAM;
        
    // specials that apply to OPTION in FASTQ but FIELD in SAMs
    if (SP(COPY_MATE_FLAG,        AGENT_RX))          // 11
        return dict_id_is_field (dict_id) ? DT_SAM : DT_FASTQ;

    // specials that apply to OPTION in SAM but FIELD in FASTQ
    if (SP(BDBI,                  mate_lookup)     || // 2
        SP(delta_seq_len,         set_deep)        || // 3
        SP(FLOAT,                 deep_copy_SEQ)   || // 5
        SP(NM,                    backspace)       || // 7
        SP(MD,                    copy_line1)      || // 8
        SP(REF_CONSUMED,          monochar_QUAL)   || // 9
        SP(COPY_BUDDY_CIGAR,      qname_rng2seq_len)) // 13
        return dict_id_is_field (dict_id) ? DT_FASTQ : DT_SAM;

    // specials that in SAM apply to only one specific field
    if (SP(BIN,                   deep_copy_QUAL))    // 6
        return dict_id.num == _SAM_BAM_BIN    ? DT_SAM : DT_FASTQ;

    if (SP(FASTQ_CONSUME_AUX,     DEMUX_by_R))        // 14
        return dict_id.num == _SAM_FQ_AUX     ? DT_SAM : DT_FASTQ;

    // specials that in FASTQ apply to only one specific field
    if (SP(CIGAR,                 unaligned_SEQ))     // 0
        return dict_id.num == _FASTQ_SQBITMAP ? DT_FASTQ : DT_SAM;
    
    return DT_SAM; // larger than the max FASTQ special

    #undef SP
}

// print one or more words in Context.dict
void dict_io_print (FILE *fp, DictId dict_id, STRp(data), bool with_word_index, bool add_quotation_marks, bool add_newline, bool remove_non_contigs)
{
    rom word = data, after = data + data_len;
    bool is_deep = flag.deep;

    for (WordIndex wi=0; word < after; wi++) {
        int word_len = strlen (word);

        DataType dt = is_deep && word_len > 2 && word[0] == SNIP_SPECIAL 
                    ? dict_io_print_deep_dt_by_special (word[1], dict_id) : DT_NONE;

        StrTextMegaLong snip = str_snip_ex (dt, STRa(word), add_quotation_marks);

        // in case we are showing chrom data in --list-chroms in SAM - don't show * and =
        bool remove_one = remove_non_contigs && 
                          (IS_ASTERISK(word) || IS_EQUAL_SIGN(word) || snip.s[0] == '[');

        if (!remove_one) {
            if (add_newline) fprintf (fp, "[%u] ", wi);
            fprintf (fp, "%s", snip.s);
        }

        word += word_len + 1;

        if (!remove_one)
            fputc ((add_newline || word == after) ? '\n' : ' ', fp);
    }

    fflush (fp);
}

void dict_io_show_singletons (VBlockP vb, ContextP ctx)
{
    if (ctx->ltype == LT_SINGLETON) {
        iprintf ("%s: %s.local contains singletons:\n", VB_NAME, ctx->tag_name);
        dict_io_print (info_stream, ctx->dict_id, STRb(ctx->local), true, true, true, false);
        iprint_newline();
    }

    else
        iprintf ("%s: %s.ltype=%s, it does not contain singletons\n", VB_NAME, ctx->tag_name, lt_name (ctx->ltype));
}