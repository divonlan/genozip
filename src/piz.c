// ------------------------------------------------------------------
//   piz.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "dispatcher.h"
#include "piz.h"
#include "random_access.h"
#include "regions.h"
#include "ref_iupacs.h"
#include "refhash.h"
#include "progress.h"
#include "profiler.h"
#include "stats.h"
#include "reconstruct.h"
#include "coverage.h"
#include "writer.h"
#include "threads.h"
#include "endianness.h"
#include "chrom.h"
#include "txtheader.h"
#include "base64.h"
#include "dict_io.h"
#include "user_message.h"

TRANSLATOR_FUNC (piz_obsolete_translator)
{
    return 0;
}

// output coordinates of current line (for error printing) - very carefully as we are in an error condition - we can't assume anything
PizDisCoords piz_dis_coords (VBlockP vb)
{
    PizDisCoords out = {};
    if (DTF(chrom) == DID_NONE || !ctx_has_value (vb, DTF(chrom))) return out;
    
    ContextP chrom_ctx = CTX(DTF(chrom));
    WordIndex chrom = chrom_ctx->last_value.i;
    if (chrom < 0 || chrom >= chrom_ctx->word_list.len) return out; // not a valid chrom value

    STR(chrom_str);
    ctx_get_snip_by_word_index (chrom_ctx, chrom, chrom_str);
    if (strlen (chrom_str) > sizeof(out.s)-20) return out;

    int out_len = snprintf (out.s, sizeof (out.s), " CHROM=\"%.64s\"(%d)", str_to_printable_(STRa(chrom_str)).s, chrom); // with leading space

    if (DTF(pos) == DID_NONE || !ctx_has_value (vb, DTF(pos))) return out;
    
    snprintf (&out.s[out_len], sizeof (out.s)-out_len, " POS=%"PRId64, CTX(DTF(pos))->last_value.i);
    return out;
}

// output a data-type-specific id of the line (for ASSPIZ) - very carefully as we are in an error condition - we can't assume anything
PizDisQname piz_dis_qname (VBlockP vb) 
{
    PizDisQname out = {};

    if (DTF(qname) != DID_NONE && ctx_encountered_in_line (vb, DTF(qname)) && !vb->preprocessing) {
        ContextP ctx = CTX(DTF(qname));
        snprintf (out.s, sizeof (out.s), " %.10s=\"%.*s\"", ctx->tag_name, MIN_(80, ctx->last_txt.len), last_txtx(vb, ctx));
    }

    return out;
}

void asspiz_text (VBlockP vb, FUNCLINE)
{
    StrTextSuperLong stack;

    char *next = stack.s;
    for (int i=0; i < vb->con_stack_len; i++)
        next += sprintf (next, "%s[%u]->", CTX(vb->con_stack[i].did_i)->tag_name, vb->con_stack[i].repeat);

    sprintf (next, "%s", (vb->curr_item != DID_NONE ? CTX(vb->curr_item)->tag_name : "N/A"));

    progress_newline(); 
    fprintf (stderr, "%s %s: Error in %s:%u line_in_file(1-based)=%"PRId64"%s %s%s stack=%s %s: ", 
             str_time().s, LN_NAME, func, code_line, 
             writer_get_txt_line_i ((VBlockP)(vb), vb->line_i), 
             cond_int (Z_DT(VCF) || Z_DT(BCF), " sample_i=", vb->sample_i), 
             piz_dis_coords((VBlockP)(vb)).s, piz_dis_qname((VBlockP)(vb)).s, stack.s, version_str().s); 
}

bool piz_grep_match (rom start, rom after)
{
    bool found = false;
    SAFE_NUL (after);

    if (!flag.grepw) {
        found = !!strstr (start, flag.grep);
        goto done;
    }

    // case: --grepw - grep whole word
    rom s = start;
    while (s <= after - flag.grep_len) {
        if (!(s = strstr (s, flag.grep))) break;

        char before = (s == start ? ' ' : s[-1]);
        char after  = s[flag.grep_len];
    
        if (!IS_ALPHANUMERIC(before) && before != '_' &&
            !IS_ALPHANUMERIC(after)  && after  != '_') {
            
            found = true;
            break;
        }

        s += flag.grep_len;
    }

done:                
    SAFE_RESTORE;
    return found;
}

bool piz_default_skip_section (SectionType st, DictId dict_id)
{
    // --show-dict=DICT - read only the one dictionary
    if (ST(DICT) && flag.show_one_dict && is_genocat && !dict_id_is_show (dict_id)) return true; // skip

    // B250, LOCAL, COUNT sections
    bool skip = is_genocat && dict_id.num 
                && dict_id.num != DTFZ(predefined)[CHROM].dict_id.num 
                && (
    
    // sometimes we don't need dictionaries. but we always load CHROM.
        (flag.genocat_no_dicts && dict_id_typeless (dict_id).num != flag.show_one_counts.num)

    // if show_counts - we only need the requested section and CHROM (note: not true for dump_one_b250_dict_id,
    // as we need to reconstruct to dump it)
    ||  (flag.show_one_counts.num && dict_id_typeless (dict_id).num != flag.show_one_counts.num)

    // if --counts, we filter here - TOPLEVEL only - unless there's a skip_section function which will do the filtering
    ||  (flag.count && !DTPZ(is_skip_section) && dict_id.num != DTFZ(toplevel).num) 
    );

    skip |= flag.dont_load_ref_file && (ST(REFERENCE) || st == SEC_REF_HASH || ST(REF_IS_SET));

    if (skip && is_genocat && dict_id.num && (dict_id.num == flag.show_singletons_dict_id.num || dict_id.num == flag.dump_one_local_dict_id.num))
        skip = false;
        
    return skip;
}

static inline void piz_adjust_one_local (ContextP ctx, BufferP local_buf, LocalType *ltype, uint8_t param, bool uncompress_to_pair) 
{ 
    const LocalTypeDesc *ltd = &lt_desc[*ltype];

    ASSERT (local_buf->len % ltd->width == 0, "%s.local has %u bytes - but expecting the number of bytes to be a multiple of %u since ltype=%s",
            ctx->tag_name, local_buf->len32, ltd->width, ltd->name);

    local_buf->len /= ltd->width; 

    // note: in ZIP, if we're loading R1.local to localR1 just to verify identicality, then we should 
    // keep it in "file format" for zip_generate_local(), not "native fomrat"
    if (!uncompress_to_pair || IS_PIZ || !fastq_zip_use_pair_identical (ctx->dict_id)) {
        
        if (*ltype == LT_BITMAP) { 
            local_buf->nbits = local_buf->len * 64 - param ; 
            LTEN_bits ((BitsP)local_buf); 
        }

        else if (*ltype >= LT_UINT8_TR && *ltype <= LT_UINT64_TR)
            local_buf->n_cols = param; // 0 means vcf_num_samples
        
        if (ltd->file_to_native)   
            ltd->file_to_native (local_buf, ltype); // BGEN, transpose etc - updates ltype in case of Transpose, after untransposing
    }
}

// PIZ compute thread: decompress all contexts (in pair-2 of paired FASTQ: z_data contains contexts of both pairs)
// ZIP compute thread in FASTQ: decompress pair_1 contexts when compressing pair_2
void piz_uncompress_all_ctxs (VBlockP vb)
{
    bool vb_is_pair_2 = is_fastq_pair_2(vb); // is this VB a pair-2 FASTQ VB (either in a FASTQ or SAM z_file)

    for_buf (uint32_t, header_offset, vb->z_section_headers) {
        SectionHeaderCtxP header = (SectionHeaderCtxP)Bc(vb->z_data, *header_offset);

        bool is_local = HEADER_IS(LOCAL);
        bool is_b250  = HEADER_IS(B250);
        if (!is_b250 && !is_local) continue;

        ContextP ctx = ctx_get_ctx (vb, header->dict_id); // gets the context (creating it if it doesn't already exist)
        
        // back comp: bug observed with E2:Z in v11.0.10: OPTION_E2_Z has LOCAL section despite being an alias to SAM_E2_Z
        if (ctx->is_ctx_alias && !VER(12)) { 
            ctx = CTX(ctx->did_i); 
            header->dict_id = ctx->dict_id; 
        }

        bool is_pair_section = vb_is_pair_2 && (BGEN32 (header->vblock_i) != vb->vblock_i); // is this a section of R1 read into an R2 vb 
        bool uncompress_to_pair = is_pair_section && (!header->flags.ctx.paired/*not pair-identical*/ || IS_ZIP); // ZIP: always; PIZ: if pair-assisted

        ASSERT (is_b250 || header->ltype < NUM_LTYPES, "in vb=%u ctx=%s.%s: ltype=%u >= NUM_LTYPES=%u. This can possibly be solved by upgrading Genozip to the latest version", 
                vb->vblock_i, ctx->tag_name, is_local ? "local" : "b250", header->ltype, NUM_LTYPES);

        // PIZ only: load normal section, or a pair-identical section of from the R1 VB
        if (!uncompress_to_pair) {
            
            // case: buffer has already been decompressed  - either during pre-processing, or 
            // R2 section was decompressed so no need for the pair-identical R1 section
            if ((is_local && ctx->local_uncompressed) || (is_b250 && ctx->b250_uncompressed))
                continue;
                        
            if (is_local) {
                ctx->lcodec = header->codec;
                ctx->ltype  = header->ltype; 
                ctx->nothing_char = !lt_max(ctx->ltype)          ? 0  // nothing char is only relevant for integer ltypes   
                                  : header->nothing_char == 0xff ? 0 // no nothing char
                                  : header->nothing_char == 0    ? 1 // use hard-coded logic up to (0 always, and only, appears in files up to 15.0.37)
                                  :                                header->nothing_char;
            }

            else { // b250
                // old logic - not clear if/why it is needed but not removed to avoid break back comp:
                if (!VER(15) && !ctx->ltype) 
                    ctx->ltype = header->ltype; 

                ctx->iterator    = (SnipIterator){ .next_b250 = B1ST8 (ctx->b250), .prev_word_index = WORD_INDEX_NONE };
                ctx->b250_size   = header->b250_size; // note: for files<=v13, this was always 0, ie B250_BYTES_4
            }
        }

        // A pair section (but only pair-assisted in PIZ)
        else {
            if (is_b250) {
                ctx->pair_b250_iter = (SnipIterator){ .next_b250 = B1ST8 (ctx->b250R1), .prev_word_index = WORD_INDEX_NONE };
                ctx->pair_b250_size = header->b250_size;
            }
            else 
                ctx->pair_ltype = header->ltype;
        }

        BufferP target_buf = uncompress_to_pair ? (is_local ? &ctx->localR1 : &ctx->b250R1)
                                                : (is_local ? &ctx->local   : &ctx->b250);
        
        rom target_buf_name = uncompress_to_pair ? (is_local ? "contexts->localR1" : "contexts->b250R1")
                                                 : (is_local ? CTX_TAG_LOCAL   : CTX_TAG_B250  );

        START_TIMER;

        zfile_uncompress_section (vb, header, target_buf, target_buf_name, BGEN32 (header->vblock_i), header->section_type); 

        if (is_local && dict_id_typeless (ctx->dict_id).num == flag.show_singletons_dict_id.num && !is_pair_section) 
            dict_io_show_singletons (vb, ctx);
            
        if (is_local && dict_id_typeless (ctx->dict_id).num == flag.dump_one_local_dict_id.num && !is_pair_section) 
            ctx_dump_binary (vb, ctx, true);

        if (!is_local && dict_id_typeless (ctx->dict_id).num == flag.dump_one_b250_dict_id.num && !is_pair_section) 
            ctx_dump_binary (vb, ctx, false);

        // BGEN32, transpose, fix len
        if (is_local && uncompress_to_pair) 
            piz_adjust_one_local (ctx, &ctx->localR1, &ctx->pair_ltype, header->param, true);
        
        else if (is_local && !uncompress_to_pair)  
            piz_adjust_one_local (ctx, &ctx->local, &ctx->ltype, header->param, false);

        if (is_local && !is_pair_section) 
            ctx->local_uncompressed = true;

        else if (!is_local && !is_pair_section) // b250
            ctx->b250_uncompressed = true;

        if (!uncompress_to_pair/*added this condition in v15*/ && 
            ((VER(14) && ctx->ltype != LT_BITMAP) ||                // starting v14: assign to all except LT_BITMAP (in which param is used to determine nbits)
             (!VER(14) && header->flags.ctx.v13_copy_local_param))) // up to v13: copy if v13_copy_local_param is set
            target_buf->prm8[0] = header->param;

        if (flag.debug_read_ctxs) 
            iprintf ("%c Uncompressed %s: %s[%u].len=%u into %s\n", sections_read_prefix (is_pair_section || vb->preprocessing), 
                     VB_NAME, ctx->tag_name, ctx->did_i, target_buf->len32, target_buf_name);
    }

    if (IS_PIZ) {
        // initialize history buffer (eg for SAM buddy)
        for_ctx_that (ctx->flags.store_per_line || ctx->flags.spl_custom) 
            switch (ctx->flags.store) {
                // we zero the history, bc when seg compares to a dl value for a field that didn't exist,
                // it sees 0. It might seg against that 0. So we need history to be 0 too.
                case STORE_INT   : buf_alloc_exact_zero (vb, ctx->history, vb->lines.len, int64_t,     "history"); break;
                case STORE_FLOAT : buf_alloc_exact_zero (vb, ctx->history, vb->lines.len, double,      "history"); break;
                case STORE_INDEX : buf_alloc_exact_zero (vb, ctx->history, vb->lines.len, WordIndex,   "history"); break;
                default          : buf_alloc_exact_zero (vb, ctx->history, vb->lines.len, HistoryWord, "history"); break;
            }

        // prepare context index
        for_ctx
            vb->ctx_index[ctx->did_i] = (ContextIndex){ .did_i = ctx->did_i, .dict_id = ctx->dict_id };
    
        qsort (vb->ctx_index, vb->num_contexts, sizeof (ContextIndex), sort_by_dict_id);

        vb->has_ctx_index = true;
    }

    if (flag.debug_or_test) buflist_test_overflows(vb, __FUNCTION__); 
}

// PIZ compute thread entry point
static void piz_reconstruct_one_vb (VBlockP vb)
{
    START_TIMER;

    ASSERTNOTNULL (vb);
    ASSERT (vb->vblock_i, "vb->vblock_i is 0: vb->compute_thread_id=%d pthread=%"PRIu64, 
            vb->compute_thread_id, (uint64_t)pthread_self());

    ASSERT (!flag.reference || ref_is_loaded (gref) || flag.dont_load_ref_file,
            "%s: reference is not loaded correctly", VB_NAME);

    ASSERT (vb->recon_size >= 0, "Invalid vb->recon_size=%d", vb->recon_size);

    // note: txt_data is fully allocated in advance and cannot be extended mid-reconstruction (container_reconstruct and possibly others rely on this)
    #define OVERFLOW_SIZE (1 MB) // allow some overflow space as sometimes we reconstruct unaccounted for data: 1. container templates 2. reconstruct_peek and others
    
    buf_alloc (vb, &vb->txt_data, 0, 
               vb->recon_size * vb->translation.factor/*see TRANSLATIONS*/ + OVERFLOW_SIZE, 
               char, 1.1, "txt_data"); 

    piz_uncompress_all_ctxs (vb);

    DT_FUNC (vb, piz_recon_init)(vb);

    // reconstruct from top level snip
    Did top_level_did_i = ctx_get_existing_did_i (vb, vb->translation.toplevel); 
    reconstruct_from_ctx (vb, top_level_did_i, 0, true);
    
    ASSERT (!vb->con_stack_len, "%s: Expecting container stack to be empty, but con_stack_len=%u", VB_NAME, vb->con_stack_len);

    // calculate the digest contribution of this VB, and the digest snapshot of this VB
    // note: if we have generated components from which lines might be inserted into the VB - we verify in writer instead 
    // note: for Deep with gencomp - the SAM components are verified in writer, while the FASTQ components are verified here.
    if (piz_need_digest && (!z_has_gencomp || VB_DT(FASTQ)) && !(flag.deep_fq_only && !VB_DT(FASTQ)))
        digest_one_vb (vb, true, NULL); // LOOKING FOR A DEADLOCK BUG? CHECK HERE

    if (DTP(piz_after_recon)) DTP(piz_after_recon)(vb);

    vb_set_is_processed (vb); /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ 

    if (flag.debug_or_test) buflist_test_overflows(vb, __FUNCTION__); 

    COPY_TIMER (compute);
}

static void piz_initialize_ctx_flags_from_vb_1 (VBlockP vb)
{
    // ctx.flags defaults to vb_i=1 flags, overridden if a b250 or local section is read. this will not be overridden if all_the_same, i.e. no b250/local sections.
    // note: we use section_list_save and not section_list_buf, because the latter might not contain vb=1, if removed by writer_create_plan
    Section vb_1_first_sec = B(SectionEnt, z_file->section_list_save, z_file->section_list_save.prm32[0]);
    Section vb_1_last_sec  = B(SectionEnt, z_file->section_list_save, z_file->section_list_save.prm32[1]);

    for (Section sec = vb_1_first_sec+1; sec <= vb_1_last_sec; sec++) {
        ContextP ctx = ECTX (sec->dict_id); // will exist if it has a dict (all_the_same sections always have a dict)
        if (ctx) {
            ctx->flags = sec->flags.ctx;
            ctx->flags.paired = false; // flags.paired is VB-specific and is not inherited
        }
    }
}

void piz_read_all_ctxs (VBlockP vb, Section *sec/* VB_HEADER section */, bool is_pair_data) 
{
    START_TIMER;

    for ((*sec)++; (*sec)->st == SEC_B250 || (*sec)->st == SEC_LOCAL; (*sec)++) {
        ASSERT (is_pair_data || vb->vblock_i == (*sec)->vblock_i, "expecting vb->vblock_i=%u == sec->vblock_i=%u", 
                vb->vblock_i, (*sec)->vblock_i); // sanity

        // create a context even if section is skipped, for containers to work (skipping a section should be mirrored in a container filter)
        ContextP zctx = ctx_get_existing_zctx ((*sec)->dict_id); 
        ContextP vctx = CTX(zctx->did_i); // in PIZ z and vb contexts always have same did_i. This is also true for ZIP of R2, bc context was created by R1 and overlayed on this R2 VB.
        bool is_local = ((*sec)->st == SEC_LOCAL);
        
        // don't assert for <=v11 due to bug (see comment in piz_uncompress_all_ctxs)
        ASSERT (!zctx->is_ctx_alias || !VER(12), "Found a %s section of %s, this is unexpected because %s is an alias (of %s)",
                st_name((*sec)->st), zctx->tag_name, zctx->tag_name, ZCTX(zctx->did_i)->tag_name);

        // if we're a FASTQ R2 VB loading R1 data, decide if we need to load this section
        bool pair_assisted=false, pair_identical=false, skip_R1=false;
        if (is_pair_data) {
            // note: pair_assisted available since early versions. when loading old versions, flags are fixed in sections_list_file_to_memory_format to be consistent with current version. 
            pair_assisted  = (IS_ZIP && fastq_zip_use_pair_assisted ((*sec)->dict_id, (*sec)->st)) ||
                             (IS_PIZ && is_local  && vctx->pair_assist_type == SEC_LOCAL) ||  // R2 section indicated that recon requires assistence of R1 data
                             (IS_PIZ && !is_local && vctx->pair_assist_type == SEC_B250);

            // note: pair_identical was introduced in v15
            pair_identical = (IS_ZIP && fastq_zip_use_pair_identical ((*sec)->dict_id)) ||
                             (IS_PIZ && is_local && !vctx->local_in_z && ((*sec)->flags.ctx.paired)) || // R2 section is missing and R1 is willing to take its stead
                             (IS_PIZ && !is_local && !vctx->b250_in_z && ((*sec)->flags.ctx.paired));

            if (!pair_assisted && !pair_identical) skip_R1 = true; // this R1 section is not needed by R2

            ASSERT (!pair_assisted || !pair_identical, "%s: %s.%s is invalidly both pair_assisted and pair_identical",
                    VB_NAME, vctx->tag_name, is_local ? "local" : "b250");
        }     

        uint32_t section_start = vb->z_data.len32;
        int32_t offset = skip_R1 ? SECTION_SKIPPED 
                                 : zfile_read_section (z_file, vb, (*sec)->vblock_i, &vb->z_data, "z_data", (*sec)->st, *sec); // returns 0 if section is skipped
        
        bool section_read = (offset != SECTION_SKIPPED); // section could be skipped either bc of skip_R1 or piz_is_skip_section() called from zfile_read_section 

        if (section_read) {
            BNXT32 (vb->z_section_headers) = section_start; 

            if  (!is_pair_data) {
                if (is_local) vctx->local_in_z = true;
                else          vctx->b250_in_z  = true;

                // note: |= (instead of =) to overcome bug in ZIP --pair in some versions <= 13 (lost track of which): 
                // local sections with junk data created in addition to the expected b250: if both b250 and local indicate pair_assist, we take the b250.
                if ((*sec)->flags.ctx.paired && is_fastq_pair_2(vb))
                    vctx->pair_assist_type = (*sec)->st; // set when reading R2 sections, consumed when reading pair R1 sections 
            }

            if (pair_assisted || (pair_identical && IS_ZIP)) 
                vctx->pair_flags = (*sec)->flags.ctx;
            else 
                vctx->flags = (*sec)->flags.ctx; // override flags inherited from vb=1 and possibly the other B250/LOCAL section
        } 

        // note: vctx->is_loaded possibly already true if it has a dictionary - set in ctx_overlay_dictionaries_to_vb - and now sets to false if section is skipped
        if (IS_PIZ && !pair_assisted && !skip_R1) 
            vctx->is_loaded = section_read; 

        if (flag.debug_read_ctxs) {
            if (section_read)
                sections_show_header ((SectionHeaderP)Bc (vb->z_data, section_start), NULL, (*sec)->offset, sections_read_prefix (is_pair_data || vb->preprocessing));
            else
                iprintf ("%c Skipped loading %s/%u %s.%s\n", sections_read_prefix (is_pair_data || vb->preprocessing), 
                         comp_name((*sec)->comp_i), vb->vblock_i, zctx->tag_name, st_name ((*sec)->st));
        }
    }

    if (flag.debug_or_test) buflist_test_overflows(vb, __FUNCTION__); 

    COPY_TIMER (piz_read_all_ctxs);
}

// Called by PIZ main thread: read all the sections at the end of the file, before starting to process VBs
DataType piz_read_global_area (Reference ref)
{
    START_TIMER;

    bool success = zfile_read_genozip_header (0, SOFT_FAIL); // already read if normal file, but not if auxilliary file

    if (flag.show_stats) {
        stats_read_and_display();
        if (is_genocat) return DT_NONE;
    }

    user_message_display();

    if (!success) return DT_NONE;

    if (flags_writer_counts()) goto done;

    // check if the genozip file includes a reference
    bool has_ref_sections = !!sections_last_sec (SEC_REFERENCE, SOFT_FAIL);

    ASSERTW (!has_ref_sections || !IS_REF_EXTERNAL || flag.reading_reference, 
             "FYI: ignoring reference file %s because %s was not compressed with --reference", ref_get_filename (ref), z_name);

    if (!flag.reading_reference && has_ref_sections) {
        ref_destroy_reference (ref);  // destroy an old reference, if one is loaded
        flag.reference = REF_STORED;  // possibly override REF_EXTERNAL (it will be restored for the next file in )
    }

    // read all dictionaries - CHROM/RNAME is needed for regions_make_chregs(). 
    // Note: some dictionaries are skipped based on skip() and all flag logic should implemented there
    dict_io_read_all_dictionaries(); 

    if (!flag.header_only) {
        // mapping of the file's chroms to the reference chroms (for files originally compressed with REF_EXTERNAL/EXT_STORE and have alternative chroms)
        chrom_2ref_load (ref); 

        ref_contigs_load_contigs (ref); // note: in case of REF_EXTERNAL, reference is already pre-loaded
    }

    // if the user wants to see only the header, we can skip regions and random access
    if (!flag.header_only) {

        ctx_read_all_counts();   // read all SEC_COUNTS sections

        ctx_read_all_subdicts(); // read all SEC_SUBDICTS sections

        // update chrom node indices using the CHROM dictionary, for the user-specified regions (in case -r/-R were specified)
        if (flag.regions && ZCTX(DTFZ(chrom))->word_list.len) // FASTQ compressed without reference, has no CHROMs 
            regions_make_chregs (ZCTX(DTFZ(chrom)));

        // if the regions are negative, transform them to the positive complement instead
        regions_transform_negative_to_positive_complement();

        // if this is a stored reference we load the reference random access that will determined which reference sections
        // should be read & uncompressed in case of --regions.
        // note: in case of a data file with stored reference - SEC_REF_RAND_ACC will contain the random access of the reference
        // and SEC_RANDOM_ACCESS will contain the random access of the data. In case of a .ref.genozip file, both sections exist 
        // and are identical. It made the coding easier and their size is negligible.
        random_access_load_ra_section (SEC_RANDOM_ACCESS, DTFZ(chrom), &z_file->ra_buf, "z_file->ra_buf", 
                                       !flag.show_index ? NULL : RA_MSG_PRIM);

        random_access_load_ra_section (SEC_REF_RAND_ACC, CHROM, ref_get_stored_ra (ref), "ref_stored_ra", 
                                       flag.show_ref_index && !flag.reading_reference ? RA_MSG_REF : NULL);

        if (IS_REF_CHROM2REF && !flag.reading_reference && !flag.genocat_no_reconstruct)
            // xxx is this actually used? 
            chrom_index_by_name (CHROM); // create alphabetically sorted index for user file (not reference) chrom word list

        // case: reading reference file
        if (flag.reading_reference) {

            // when reading the reference for genocat --coverage/idxstats, don't need the actual REF sections 
            if (is_genocat && (flag.show_coverage || flag.idxstats)) 
                goto done;  

            bool ref_loaded_from_disk = !flag.dont_load_ref_file && ref_load_stored_reference (ref);

            // load the IUPACs list of the reference (rare non-ACGT "bases")
            ref_iupacs_load (ref);

            // load the refhash, if we are compressing FASTA or FASTQ, or if user requested to see it
            if (  (primary_command == ZIP && flag.aligner_available) ||
                  (flag.show_ref_hash && is_genocat) ||
                  ref_cache_is_populating (ref))
                refhash_load (ref);

            // exit now if all we wanted was just to see the reference (we've already shown it)
            if ((flag.show_reference || flag.show_is_set || flag.show_ref_hash) && is_genocat) exit_ok;

            if (ref_loaded_from_disk) progress_finalize_component ("Done");
        }

        // case: non-reference file has stored reference sections
        else if (has_ref_sections) {
            if (!flag.dont_load_ref_file) { 
                ref_load_stored_reference (gref);

                // exit now if all we wanted was just to see the reference (we've already shown it)
                if ((flag.show_reference || flag.show_is_set || flag.show_ref_hash) && is_genocat) exit_ok;
            }
        }
    }
    
done:
    COPY_TIMER_EVB (piz_read_global_area);

    return z_file->data_type;
}

// main thread
bool piz_read_one_vb (VBlockP vb, bool for_reconstruction)
{
    START_TIMER; 
   
    Section sec = sections_vb_header (vb->vblock_i); 
    
    int32_t vb_header_offset = zfile_read_section (z_file, vb, vb->vblock_i, &vb->z_data, "z_data", SEC_VB_HEADER, sec); 
    ASSERT0 (vb_header_offset >= 0, "Unexpectedly VB_HEADER section was skipped");

    SectionHeaderVbHeader header = *(SectionHeaderVbHeaderP)Bc (vb->z_data, vb_header_offset); // copy of header as it will be overwritten in piz_read_all_ctxs

    // any of these might be overridden by callback
    vb->flags            = header.flags.vb_header;
    vb->recon_size       = BGEN32 (header.recon_size);   
    vb->longest_line_len = BGEN32 (header.longest_line_len);
    vb->longest_seq_len  = VER(15) ? BGEN32 (header.longest_seq_len) : 0;
    vb->expected_digest  = header.digest;
    vb->chrom_node_index = WORD_INDEX_NONE;
    vb->lines.len        = VER(14) ? sec->num_lines : BGEN32 (header.v13_top_level_repeats);
    vb->comp_i           = sec->comp_i; 
    vb->show_containers  = (flag.show_containers == SHOW_CONTAINERS_ALL_VBs || flag.show_containers == vb->vblock_i); // a per-VB value bc in SAM Load-Prim VBs =false vs normal VBs have the flag value (set in sam_piz_dispatch_one_load_sag_vb)

    if (txt_file) { // sometimes we don't have a txtfile, eg when genocat is used with some flags that emit other data, no the file
        vb->vb_position_txt_file = txt_file->txt_data_so_far_single_0; // position in original txt file (before any ZIP or PIZ modifications)
        txt_file->num_lines += vb->lines.len; // source file lines
    }

    // in case of unbind, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every txt component
    if (flag.unbind) vb->vblock_i = BGEN32 (header.vblock_i);

    if (flag_is_show_vblocks (PIZ_TASK_NAME)) 
        iprintf ("READING(id=%d) vb=%s num_lines=%u recon_size=%u genozip_size=%u longest_line_len=%u\n",
                 vb->id, VB_NAME, vb->lines.len32, vb->recon_size, BGEN32 (header.z_data_bytes), vb->longest_line_len);

    ctx_overlay_dictionaries_to_vb (VB); // overlay all dictionaries to the vb 

    buf_alloc (vb, &vb->z_section_headers, MAX_DICTS * 2, 0, uint32_t, 0, "z_section_headers"); // room for section headers  

    BNXT32 (vb->z_section_headers) = vb_header_offset; // vb_header_offset is always 0 for VB header

    piz_initialize_ctx_flags_from_vb_1 (vb);

    DT_FUNC (vb, piz_before_read)(vb);

    // read all b250 and local of all fields and subfields
    piz_read_all_ctxs (vb, &sec, false);

    bool ok_to_compute = DT_FUNC_OPTIONAL (vb, piz_init_vb, true)(vb, &header);

    vb->translation = dt_get_translation (vb); // must be after piz_init_vb, as in VCF we set vb->vb_chords there, needed for dt_get_translation

    if (txt_file) 
        txt_file->txt_data_so_far_single_0 += BGEN32 (header.recon_size); // cumulative expected recon size without piz-side modifications

    if (ok_to_compute && for_reconstruction && flag.collect_coverage) 
        coverage_initialize (vb);

    COPY_TIMER (piz_read_one_vb); 

    return ok_to_compute;
}

static void piz_handover_or_discard_vb (Dispatcher dispatcher, VBlockP *vb)
{
    bool is_handed_over = false;

    if ((*vb)->preprocessing) 
        DT_FUNC (z_file, piz_after_preproc)(*vb);
    
    else if (!flag.no_writer_thread)  // note: in SAM with gencomp - writer does the digest calculation
        is_handed_over = writer_handover_data (vb);

    if (!is_handed_over && !(*vb)->preprocessing &&
        (!flag.one_component || writer_does_vb_need_write ((*vb)->vblock_i)))
        txt_file->txt_data_so_far_single += (*vb)->txt_data.len; // note: if writing (or SAM with gencomp), this is done in writer_flush_vb, caputring the processing in writer too

    dispatcher_recycle_vbs (dispatcher, !is_handed_over); // don't release VB if handed over - it will be released in writer_release_vb when writing is completed
}

// returns false if VB was dispatched, and true if vb was skipped
static void piz_dispatch_one_vb (Dispatcher dispatcher, Section sec)
{
    VBlockP next_vb = dispatcher_generate_next_vb (dispatcher, sec->vblock_i, sec->comp_i);

    // read one VB's data from z_file
    ReconType reconstruct = piz_read_one_vb (next_vb, true) && // read even if no_reconstruct
                            !flag.genocat_no_reconstruct; 

    if (reconstruct) {
        if (flag_is_show_vblocks (PIZ_TASK_NAME)) 
            iprintf ("BEFORE_COMPUTE(id=%d) vb=%s/%u num_running_compute_threads(before)=%u\n", 
                     next_vb->id, comp_name(next_vb->comp_i), next_vb->vblock_i, dispatcher_get_num_running_compute_threads(dispatcher));

        dispatcher_compute (dispatcher, piz_reconstruct_one_vb);
        dispatcher_increment_progress ("read", 1); // done reading
    }

    // case: we won't proceed to uncompressing, reconstructing we're done reading - just handover
    // an empty VB as it appears in the recon plan, and writer might be already blocking on waiting for it
    else {
        dispatcher_increment_progress ("all_no_reconstruct", 3); // done reading, skipped reconstructing and writing
        dispatcher_abandon_next_vb (dispatcher); // just moves the to processed_vb so dispatcher_recycle_vbs can recycle it
        piz_handover_or_discard_vb (dispatcher, &next_vb);
    }
}

// main thread: usually called in order of VBs, but out-of-order if --test with no writer
static void piz_handle_reconstructed_vb (Dispatcher dispatcher, VBlockP vb, uint64_t *num_nondrop_lines)
{
    // verify that files are the same size, unless we intended to modify the data
    ASSERTW (Ltxt == vb->recon_size || flag.piz_txt_modified || (flag.deep_fq_only && !VB_DT(FASTQ)),
            "Warning: vblock_i=%s/%u (num_lines=%u) had %s bytes in the original %s file but %s bytes in the reconstructed file (diff=%d)", 
            comp_name (vb->comp_i), vb->vblock_i, vb->lines.len32, str_int_commas (vb->recon_size).s, dt_name (txt_file->data_type), 
            str_int_commas (Ltxt).s, 
            (int32_t)Ltxt - (int32_t)vb->recon_size);

    *num_nondrop_lines += vb->num_nondrop_lines;
    if (flag.count == CNT_VBs)
        iprintf ("vb=%s lines=%u nondropped_lines=%u txt_data.len=%u\n", 
                 VB_NAME, vb->lines.len32, vb->num_nondrop_lines, Ltxt);

    DT_FUNC (vb, piz_process_recon)(vb);
    
    z_file->txt_data_so_far_single += vb->recon_size; 

    piz_handover_or_discard_vb (dispatcher, &vb);
}

// dispatcher of PIZ of main (i.e. not auxilliary) files
static Dispatcher main_dispatcher = NULL; 
void piz_set_main_dispatcher (Dispatcher dispatcher)
{
    main_dispatcher = dispatcher;
}

// allow out of order joining of VBs (to reverse non-allowing set in dispatcher_init)
void piz_allow_out_of_order (void)
{
    ASSERTNOTNULL (main_dispatcher);
    dispatcher_allow_out_of_order (main_dispatcher);
}

static uint64_t piz_target_progress (CompIType comp_i)
{
    if (comp_i == COMP_MAIN && Z_DT(SAM))      
        return 3 * sections_get_num_vbs_(SAM_COMP_MAIN, SAM_COMP_DEPN) + sections_get_num_vbs (SAM_COMP_PRIM); // VBs pre-processed
        
    else if (Z_DT(FASTQ) && flag.interleaved) 
        return 3 * sections_get_num_vbs_(FQ_COMP_R1, FQ_COMP_R2);
    
    else if (Z_DT(SAM) && flag.deep && flag.interleaved)
        return 3 * sections_get_num_vbs_(SAM_COMP_FQ00, SAM_COMP_FQ01);

    else {
        if (comp_i == COMP_NONE) comp_i = COMP_MAIN;
        return 3 * sections_get_num_vbs_(comp_i, comp_i);
    }
}

Dispatcher piz_z_file_initialize (void)
{
    // read all non-VB non-TxtHeader sections
    DataType data_type = piz_read_global_area (gref);
    if (data_type == DT_NONE || flag.reading_reference) 
        return NULL; // no components in this file (as is always the case with reference files) 

    if (flag.genocat_global_area_only) return NULL;

    if (!flag_loading_auxiliary && DTPZ(piz_after_global_area)) // must be before writer_create_plan messes up the section list
        DTPZ(piz_after_global_area)();

    if (!writer_create_plan())
        return NULL; // --count, and it was reported already

    if (flag.test || flag.md5) 
        ASSINP0 (dt_get_translation(NULL).is_src_dt, "Error: --test or --md5 cannot be used when converting a file to another format"); 

    // note: if --unbind, we will recalculate the target progress in dispatcher_resume()
    Dispatcher dispatcher = dispatcher_init (flag.reading_reference ? "piz-ref" : PIZ_TASK_NAME, // also referred to in dispatcher_recycle_vbs()
                                             PREPROCESSING_TASK_NAME, 
                                             POOL_MAIN,
                                             flag.xthreads ? 1 : global_max_threads, 0, 
                                             flag.test && flag.no_writer_thread, // out-of-order if --test with no writer thread (note: SAM gencomp always have writer thread to do digest). 
                                             flag.test,
                                             flag.out_filename ? flag.out_filename : txtheader_get_txt_filename_from_section().s, 
                                             piz_target_progress (COMP_MAIN), 
                                             0);

    return dispatcher;
}

// main thread: called once per txt_file created: i.e. once, except if unbinding a paired FASTQ, or a Deep file.
// returns true if piz completed, false if piz aborted by piz_initialize
bool piz_one_txt_file (Dispatcher dispatcher, bool is_first_z_file, bool is_last_z_file,
                       CompIType first_comp_i, CompIType last_comp_i, // COMP_NONE unless flag.unbind
                       bool allow_skip_cleaning)
{
    dispatcher_start_wallclock();
    
    recon_stack_initialize();
    
    if (DTPZ(piz_initialize) && !DTPZ(piz_initialize)(first_comp_i))
        return false; // abort PIZ if piz_initialize says so
      
    bool header_only_file = true; // initialize - true until we encounter a VB header
    uint64_t num_nondrop_lines = 0;

    // note: may be NULL if txt_header was removed by writer, eg when loading auxillary files
    Section txt_header_sec = (first_comp_i != COMP_NONE) ? sections_get_comp_txt_header_sec (first_comp_i) : NULL;

    Section sec = sections_one_before (txt_header_sec);

    // traverse section list as re-arranged by writer_create_plan
    while (!dispatcher_is_done (dispatcher)) {

        bool achieved_something = false;
        
        // we're pre-processing data (SAM: loading sag)
        if (flag.preprocessing && dispatcher_has_free_thread (dispatcher) && !vb_pool_is_full (POOL_MAIN))
            achieved_something = DTPZ(piz_preprocess)(dispatcher);
    
        // In input is not exhausted, and a compute thread is available - read a vblock and dispatch it
        else if (!dispatcher_is_input_exhausted (dispatcher) && dispatcher_has_free_thread (dispatcher) && !vb_pool_is_full (POOL_MAIN)) {
            achieved_something = true;

            bool found_header = sections_next_sec2 (&sec, SEC_TXT_HEADER, SEC_VB_HEADER);
            bool is_sec_in_comp = (first_comp_i==COMP_NONE || (sec->comp_i >= first_comp_i && sec->comp_i <= last_comp_i));

            // case SEC_TXT_HEADER
            if (found_header && sec->st == SEC_TXT_HEADER && is_sec_in_comp) { 
                if (sec->vblock_i >= 2) continue; // fragments >= 2 were already handled together with the first fragment
                
                // note: also starts writer, and if unbinding, also opens the txt file and hands data over to the writer
                txtheader_piz_read_and_reconstruct (sec); 

                // case --unbind: unpausing after previous txt_file pause (requires txt file to be open)
                uint64_t target_progress = piz_target_progress (first_comp_i ? first_comp_i : COMP_MAIN);
                dispatcher_resume (dispatcher, target_progress);
            }

            // case SEC_VB_HEADER
            else if (found_header && sec->st == SEC_VB_HEADER && is_sec_in_comp) {
                
                if (!writer_does_vb_need_recon (sec->vblock_i)) {
                    dispatcher_increment_progress ("vb_no_recon", 3); // skipped reading, reconstructing, writing

                    if (flag_is_show_vblocks (PIZ_TASK_NAME)) 
                        iprintf ("SKIPPED vb=%s/%u\n", comp_name (sec->comp_i), sec->vblock_i); 
                    continue;
                }

                piz_dispatch_one_vb (dispatcher, sec);
                header_only_file = false;
            }

            // case: we're done with this txt_file (either no header bc EOF, or TXT_HEADER belongs to the next txt_file when unbinding)
            else {
                if (flag_is_show_vblocks (PIZ_TASK_NAME)) 
                    iprintf ("INPUT EXHAUSTED - no more SEC_VB_HEADER or SEC_TXT_HEADER for txt_file_i=%u\n", z_file->num_txts_so_far);                

                dispatcher_set_no_data_available (dispatcher, false, DATA_EXHAUSTED);

                if (header_only_file)
                    dispatcher_recycle_vbs (dispatcher, true); // note: this is normally done in piz_handover_or_discard_vb
            }
        }

        // if the next thread (usually in order or VBs in recon_plan, but out-of-order if --test with no writer) is ready, handle the reconstructed VB
        VBlockP vb = dispatcher_get_processed_vb (dispatcher, NULL, false);  // non-blocking
        if (vb) {    
            if (flag_is_show_vblocks (PIZ_TASK_NAME)) 
                iprintf ("AFTER_COMPUTE(task=piz id=%d) vb=%s/%u num_running_compute_threads(after)=%u\n", 
                         vb->id, comp_name(vb->comp_i), vb->vblock_i, dispatcher_get_num_running_compute_threads(dispatcher));
            
            dispatcher_increment_progress ("preproc_or_recon", 1 + (!vb->preprocessing && flag.no_writer_thread)); // done preprocessing or reconstructing (+1 if skipping writing)
            if (!vb->preprocessing) 
                piz_handle_reconstructed_vb (dispatcher, vb, &num_nondrop_lines);
            else
                piz_handover_or_discard_vb (dispatcher, &vb);
        }

        if (!achieved_something) {
            START_TIMER;
            usleep (30000); // nothing for us to do right now - wait 30ms
            COPY_TIMER_EVB (piz_main_loop_idle);
        }
    }

    // make sure memory writes by compute threads are visible to the main thread
    __atomic_thread_fence (__ATOMIC_ACQUIRE); 

    if (flag_is_show_vblocks (PIZ_TASK_NAME)) 
        iprintf ("DISPATCHER is done for txt_file_i=%u: %s\n", z_file->num_txts_so_far, txt_file ? txt_file->name : "(no filename)");                

    dispatcher_calc_avg_compute_vbs (dispatcher);

    z_file->num_txts_so_far++;

    // finish writing the txt_file (note: the writer thread also calculates digest in SAM/BAM with PRIM/DEPN)
    writer_finish_writing (z_file->num_txts_so_far == z_file->num_txt_files || is_genocat);

    // verify amount of data written (if writing) or reconstructed (if --test) sizes adds up as expected 
    ASSINP (txt_file->txt_data_so_far_single/*accumulated when reconstructing/writing*/ == 
            txt_file->txt_data_so_far_single_0/*accummulated from section headers    */ || flag.piz_txt_modified,
            "Data integrity error: Size of original file (without source compression) was %"PRIu64", but reconstructed file is %"PRIu64,
            txt_file->txt_data_so_far_single_0, txt_file->txt_data_so_far_single);

    // verifies reconstructed file against MD5 or Adler2 and/or codec_args (if bgzf)
    if (piz_need_digest)
        digest_piz_verify_one_txt_file (z_file->num_txts_so_far - 1);

    progress_finalize_component_time ("Done", DIGEST_NONE);

    // --sex and --coverage - output results
    if (txt_file && !flag_loading_auxiliary) {
        if (flag.show_coverage) coverage_show_coverage();
        if (flag.idxstats) coverage_show_idxstats();
        if (flag.count == CNT_TOTAL) iprintf ("%"PRIu64"\n", num_nondrop_lines);
    }

    if (is_genocat || (z_file->num_txts_so_far == z_file->num_txt_files)) // genocat always produces exactly one txt file 
        dispatcher_finish (&dispatcher, NULL, !is_last_z_file || flag.test,
                           flag.show_memory && is_last_z_file);
    else 
        dispatcher_pause (dispatcher); // we're unbinding and still have more txt_files
    
    if (txt_file)
        DT_FUNC (txt_file, piz_finalize)();

    if (flag_is_show_vblocks (PIZ_TASK_NAME)) 
        iprintf ("Finished PIZ of %s\n", txt_file ? txt_file->name : "(no filename)");                

    // if we're loading an aux file for ZIP - destroy VBs as contexts are unions of ZIP and PIZ
    if (primary_command == ZIP && flag_loading_auxiliary)
        vb_destroy_pool_vbs (POOL_MAIN, true);
        
    file_close (&txt_file); 

    if (flag.show_time && ((flag.show_time_comp_i >= first_comp_i && flag.show_time_comp_i <= last_comp_i) || 
                           (first_comp_i == COMP_NONE && flag.show_time_comp_i != COMP_ALL)))
        profiler_add_evb_and_print_report();

    return true;
}
