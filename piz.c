// ------------------------------------------------------------------
//   piz.c
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted,
//   under penalties specified in the license.

#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "vblock.h"
#include "dispatcher.h"
#include "context.h"
#include "file.h"
#include "piz.h"
#include "sections.h"
#include "random_access.h"
#include "regions.h"
#include "dict_id.h"
#include "reference.h"
#include "ref_iupacs.h"
#include "refhash.h"
#include "progress.h"
#include "profiler.h"
#include "stats.h"
#include "codec.h"
#include "bgzf.h"
#include "flags.h"
#include "reconstruct.h"
#include "coverage.h"
#include "writer.h"
#include "strings.h"
#include "threads.h"
#include "endianness.h"
#include "website.h"
#include "chrom.h"
#include "txtheader.h"
#include "base64.h"
#include "dict_io.h"

bool piz_digest_failed = false;

// output coordinates of current line (for error printing) - very carefully as we are in an error condition - we can't assume anything
PizDisCoords piz_dis_coords (VBlockP vb)
{
    PizDisCoords out = {};
    if (DTF(prim_chrom) == DID_NONE || !ctx_has_value (vb, DTF(prim_chrom))) return out;
    
    ContextP chrom_ctx = CTX(DTF(prim_chrom));
    WordIndex chrom = chrom_ctx->last_value.i;
    if (chrom < 0 || chrom >= chrom_ctx->word_list.len) return out; // not a valid chrom value

    STR(chrom_str);
    ctx_get_snip_by_word_index (chrom_ctx, chrom, chrom_str);
    if (strlen (chrom_str) > sizeof(out.s)-20) return out;

    char printable_chrom[2*chrom_str_len+1];
    str_to_printable (STRa(chrom_str), printable_chrom);

    sprintf (out.s, " CHROM=\"%s\"(%d)", printable_chrom, chrom); // with leading space

    if (DTF(pos) == DID_NONE || !ctx_has_value (vb, DTF(pos))) return out;
    
    sprintf (&out.s[strlen(out.s)], " POS=%"PRId64, CTX(DTF(pos))->last_value.i);
    return out;
}

// output a data-type-specific id of the line (for ASSPIZ) - very carefully as we are in an error condition - we can't assume anything
PizDisQname piz_dis_qname (VBlockP vb) 
{
    PizDisQname out = {};

    if (DTF(qname) != DID_NONE && ctx_encountered_in_line (vb, DTF(qname)) && !vb->preprocessing) {
        ContextP ctx = CTX(DTF(qname));
        sprintf (out.s, " %.10s=\"%.*s\"", ctx->tag_name, MIN_(80, ctx->last_txt.len), last_txtx(vb, ctx));
    }

    return out;
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
    
        if (!IS_LETTER(before) && !IS_DIGIT(before) && before != '_' &&
            !IS_LETTER(after) && !IS_DIGIT(after) && after != '_') {
            
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
    if (st == SEC_DICT && flag.show_one_dict && is_genocat && !dict_id_is_show (dict_id)) return true; // skip

    // B250, LOCAL, COUNT sections
    bool skip = is_genocat && dict_id.num 
                && dict_id.num != DTFZ(predefined)[CHROM].dict_id.num 
                && (!flag.luft || dict_id.num != DTFZ(predefined)[DTFZ(luft_chrom)].dict_id.num) && (
    
    // sometimes we don't need dictionaries. but we always load CHROM.
        (flag.genocat_no_dicts && dict_id_typeless (dict_id).num != flag.show_one_counts.num)

    // if show_counts - we only need the requested section and CHROM (note: not true for dump_one_b250_dict_id,
    // as we need to reconstruct to dump it)
    ||  (flag.show_one_counts.num && dict_id_typeless (dict_id).num != flag.show_one_counts.num)

    // if --counts, we filter here - TOPLEVEL only - unless there's a skip_section function which will do the filtering
    ||  (flag.count && !DTPZ(is_skip_section) && dict_id.num != DTFZ(toplevel).num) 
    );

    skip |= flag.genocat_no_ref_file && (st == SEC_REFERENCE || st == SEC_REF_HASH || st == SEC_REF_IS_SET);

    if (skip && (is_genocat) && dict_id.num && dict_id.num == flag.dump_one_local_dict_id.num)
        skip = false;
        
    return skip;
}

static inline void piz_adjust_one_local (ContextP ctx, BufferP local_buf, uint8_t num_bits) 
{ 
    const LocalTypeDesc *ltd = &lt_desc[ctx->ltype];

    ASSERT (local_buf->len % ltd->width == 0, "%s.local has %u bytes - but expecting the number of bytes to be a multiple of %u since ltype=%s",
            ctx->tag_name, local_buf->len32, ltd->width, ltd->name);

    local_buf->len /= ltd->width; 

    if (ctx->ltype == LT_BITMAP) { 
        local_buf->nbits = local_buf->len * 64 - num_bits ; 
        LTEN_bits ((BitsP)local_buf); 
    } 
    else if (ltd->file_to_native)   
        ltd->file_to_native (local_buf, &ctx->ltype); // BGEN, transpose etc - updates ltype in case of Transpose, after untransposing
}

// PIZ compute thread: decompress all contexts (in pair-2 of paired FASTQ: z_data contains contexts of both pairs)
// ZIP compute thread in FASTQ: decompress pair_1 contexts when compressing pair_2
uint32_t piz_uncompress_all_ctxs (VBlockP vb)
{
    ARRAY (const uint32_t, section_index, vb->z_section_headers);
    
    uint32_t section_i = (IS_ZIP) ? 0 : 1; // normally, we skip the VB header (starting from 1), but when uncompressing paired sections there is no VB header
    for ( ; section_i < vb->z_section_headers.len32; section_i++) {

        SectionHeaderCtx *header = (SectionHeaderCtx *)Bc (vb->z_data, section_index[section_i]);

        bool is_local = (header->h.section_type == SEC_LOCAL);
        bool is_b250  = (header->h.section_type == SEC_B250);
        if (!is_b250 && !is_local) break;

        ContextP ctx = ctx_get_ctx (vb, header->dict_id); // gets the context (creating it if it doesn't already exist)

        bool is_pair_section = (BGEN32 (header->h.vblock_i) != vb->vblock_i); // is this a section of "pair 1" 

        if (!is_pair_section) {
            
            // case: buffer has already been decompressed during pre-processing, no need to decompress again
            if ((is_local && ctx->local_uncompressed) || (is_b250 && ctx->b250_uncompressed))
                continue;

            ctx->flags = header->h.flags.ctx; // overrides default inherited from vb_i=1 (assigned in piz_read_all_ctxs)
            
            if (is_local || !ctx->ltype) // a b250 can set ltype if its not already set by an earlier SEC_LOCAL section  
                ctx->ltype = header->ltype; // in case of b250
            
            if (is_local) 
                ctx->lcodec = header->h.codec;

            else { // b250
                ctx->iterator  = (SnipIterator){ .next_b250 = B1ST8 (ctx->b250), .prev_word_index = WORD_INDEX_NONE };
                ctx->b250_size = header->b250_size; // note: for files<=v13, this was always 0, ie B250_BYTES_4
            }
        }

        // initialize pair stuff (happens in FASTQ only)
        else {
            // overcome bug in ZIP --pair (in some versions <= 13, lost track of which) - local sections with junk data 
            // created in addition to the expected b250. we ignore these local sections (or allow b250 to overwrite them)    
            if (is_local && ctx->pair_b250) continue;

            ctx->pair_flags = header->h.flags.ctx;            
            ctx->pair_b250  = is_b250;                
            ctx->pair_local = is_local;

            if (is_b250) {
                ctx->pair_b250_iter = (SnipIterator){ .next_b250 = B1ST8 (ctx->pair), .prev_word_index = WORD_INDEX_NONE };
                ctx->pair_b250_size = header->b250_size;
            }
        }

        BufferP target_buf = is_local ? &ctx->local : &ctx->b250;

        START_TIMER;

        if (is_pair_section) vb->preprocessing = true; // inform fastq_piz_is_skip_section
        zfile_uncompress_section (vb, header, 
                                  is_pair_section ? &ctx->pair      : target_buf, 
                                  is_pair_section ? "context->pair" : is_local ? "contexts->local" : "contexts->b250", 
                                  BGEN32 (header->h.vblock_i),
                                  header->h.section_type); 

        if (flag.show_time && exe_type != EXE_GENOZIP)
            ctx->compressor_time = CHECK_TIMER;

        ctx_add_compressor_time_to_zf_ctx (vb);

        if (is_local && dict_id_typeless (ctx->dict_id).num == flag.dump_one_local_dict_id.num && !is_pair_section) 
            ctx_dump_binary (vb, ctx, true);

        if (!is_local && dict_id_typeless (ctx->dict_id).num == flag.dump_one_b250_dict_id.num && !is_pair_section) 
            ctx_dump_binary (vb, ctx, false);

        // BGEN32, transpose, fix len
        if (is_local) {
            if (is_pair_section) piz_adjust_one_local (ctx, &ctx->pair, header->param);
            else {
                piz_adjust_one_local (ctx, &ctx->local, header->param);
                ctx->local_uncompressed = true;
            }
        }
        else if (!is_pair_section) // b250
            ctx->b250_uncompressed = true;

        if ((VER(14) && ctx->ltype != LT_BITMAP) ||                // starting v14: assign to all except LT_BITMAP (in which param is used to determine nbits)
            (!VER(14) && header->h.flags.ctx.v13_copy_local_param)) // up to v13: copy if v13_copy_local_param is set
            target_buf->prm8[0] = header->param;

        if (flag.debug_read_ctxs)
            iprintf ("%c Uncompressed %s: %s[%u].%s.len=%"PRIu64"\n", sections_read_prefix, 
                     VB_NAME, ctx->tag_name, ctx->did_i, is_local ? "local" : "b250", 
                     (is_pair_section ? &ctx->pair : target_buf)->len);

        if (is_pair_section) vb->preprocessing = false;
    }

    // initialize history buffer (eg for SAM buddy)
    for_ctx
        if (ctx->flags.store_per_line || ctx->flags.spl_custom) 
            switch (ctx->flags.store) {
                // we zero the history, bc when seg compares to a dl value for a field that didn't existed,
                // it sees 0. It might seg against that 0. So we need history to be 0 too.
                case STORE_INT   : buf_alloc_exact_zero (vb, ctx->history, vb->lines.len, int64_t,     "history"); break;
                case STORE_INDEX : buf_alloc_exact_zero (vb, ctx->history, vb->lines.len, WordIndex,   "history"); break;
                default          : buf_alloc_exact_zero (vb, ctx->history, vb->lines.len, HistoryWord, "history"); break;
            }

    // prepare context index
    if (IS_PIZ) {
        for (Did did_i=0; did_i < vb->num_contexts; did_i++)
            vb->ctx_index[did_i] = (ContextIndex){ .did_i = did_i, .dict_id = CTX(did_i)->dict_id };
    
        qsort (vb->ctx_index, vb->num_contexts, sizeof (ContextIndex), sort_by_dict_id);

        vb->has_ctx_index = true;
    }

    return section_i;
}

// PIZ compute thread entry point
static void piz_reconstruct_one_vb (VBlockP vb)
{
    START_TIMER;

    ASSERTNOTNULL (vb);
    ASSERT (vb->vblock_i, "vb->vblock_i is 0: vb->compute_thread_id=%d pthread=%"PRIu64, 
            vb->compute_thread_id, (uint64_t)pthread_self());

    ASSERT (!flag.reference || ref_is_loaded (gref) || flag.genocat_no_ref_file,
            "%s: reference is not loaded correctly", VB_NAME);

    ASSERT (vb->recon_size >= 0, "Invalid vb->recon_size=%d", vb->recon_size);

    // note: txt_data is fully allocated in advance and cannot be extended mid-reconstruction (container_reconstruct_do and possibly others rely on this)
    #define OVERFLOW_SIZE (1 << 20) // allow some overflow space as sometimes we reconstruct unaccounted for data: 1. container templates 2. reconstruct_peek and others

    buf_alloc (vb, &vb->txt_data, 0, vb->recon_size * vb->translation.factor + OVERFLOW_SIZE, char, 1.1, "txt_data"); 
    
    piz_uncompress_all_ctxs (vb);

    DT_FUNC (vb, piz_recon_init)(vb);

    // reconstruct from top level snip
    Did top_level_did_i = ctx_get_existing_did_i (vb, vb->translation.toplevel); 
    reconstruct_from_ctx (vb, top_level_did_i, 0, true);

    // calculate the digest contribution of this VB, and the digest snapshot of this VB
    // note: if we have generated components from which lines might be inserted into the VB - we verify in writer instead 
    if (piz_need_digest && !z_has_gencomp) 
        digest_one_vb (vb, true, NULL); // LOOKING FOR A DEADLOCK BUG? CHECK HERE

    if (DTPZ(piz_after_recon)) DTPZ(piz_after_recon)(vb);

    vb_set_is_processed (vb); /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ 

    COPY_TIMER (compute);
}

void piz_read_all_ctxs (VBlockP vb, Section *sec/*first VB section after VB_HEADER */, bool is_pair_data) 
{
    if (!is_pair_data) {
        // ctxs that have dictionaries are already initialized, but others (eg local data only) are not
        ctx_initialize_predefined_ctxs (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_contexts);

        // ctx.flags defaults to vb_i=1 flags, overridden if a b250 or local section is read. this will not be overridden if all_the_same, i.e. no b250/local sections.
        // note: we use section_list_vb1 and not section_list_buf, because the latter might not contain vb=1, if removed by writer_create_plan
        for (Section sec = B1ST (SectionEnt, z_file->section_list_vb1)+1; sec < BAFT (SectionEnt, z_file->section_list_vb1); sec++) {
            ContextP ctx = ECTX (sec->dict_id); // will exist if it has a dict (all_the_same sections always have a dict)
            if (ctx) ctx->flags = sec->flags.ctx;
        }
    }
    else 
        vb->preprocessing = true; // inform fastq_piz_is_skip_section

    while ((*sec)->st == SEC_B250 || (*sec)->st == SEC_LOCAL) {
        uint32_t section_start = vb->z_data.len32;
        *B32 (vb->z_section_headers, vb->z_section_headers.len) = section_start; 

        ASSERT (is_pair_data || vb->vblock_i == (*sec)->vblock_i, "expecting vb->vblock_i=%u == sec->vblock_i=%u", vb->vblock_i, (*sec)->vblock_i); // sanity

        // create a context even if section is skipped, for containers to work (skipping a section should be mirrored in a container filter)
        ContextP zctx = ctx_get_ctx_do (z_file->contexts, z_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_contexts, (*sec)->dict_id, 0, 0);
        ContextP vctx = ctx_get_ctx (vb, zctx->dict_id);

        int32_t offset = zfile_read_section (z_file, vb, (*sec)->vblock_i, &vb->z_data, "z_data", (*sec)->st, *sec); // returns 0 if section is skipped
        
        if (offset != SECTION_SKIPPED) {
            vb->z_section_headers.len++;
            
            // mark as not skipped 
            if (!is_pair_data)
                vctx->is_loaded = true; // not skipped. note: possibly already true if it has a dictionary - set in ctx_overlay_dictionaries_to_vb

            if (flag.debug_read_ctxs)
                sections_show_header ((SectionHeaderP)Bc (vb->z_data, section_start), NULL, (*sec)->offset, sections_read_prefix);
        }
        else {
            if (!is_pair_data)
                vctx->is_loaded = false; // skipped 

            if (flag.debug_read_ctxs) 
                iprintf ("%c Skipped loading %s/%u %s.%s\n", sections_read_prefix, 
                         comp_name((*sec)->comp_i), vb->vblock_i, zctx->tag_name, (*sec)->st==SEC_LOCAL ? "local" : "b250");
        }
        
        (*sec)++;                             
    }

    if (is_pair_data) vb->preprocessing = false;
}

// Called by PIZ main thread: read all the sections at the end of the file, before starting to process VBs
DataType piz_read_global_area (Reference ref)
{
    START_TIMER;

    bool success = zfile_read_genozip_header (0);

    if (flag.show_stats) stats_read_and_display();

    if (!success) return DT_NONE;

    // check if the genozip file includes a reference
    bool has_ref_sections = !!sections_last_sec (SEC_REFERENCE, true);

    ASSERTW (!has_ref_sections || !IS_REF_EXTERNAL || flag.reading_reference, 
             "%s: ignoring reference file %s because it was not compressed with --reference", z_name, ref_get_filename (ref));

    if (!flag.reading_reference && has_ref_sections) {
        ref_destroy_reference (ref, false);  // destroy an old reference, if one is loaded
        flag.reference = REF_STORED;         // possibly override REF_EXTERNAL (it will be restored for the next file in )
    }

    // read all dictionaries - CHROM/RNAME is needed for regions_make_chregs(). 
    // Note: some dictionaries are skipped based on skip() and all flag logic should implemented there
    dict_io_read_all_dictionaries(); 

    if (!flag.header_only || z_is_dvcf) { // dual coordinates need this stuff of for the rejects part of the header

        // mapping of the file's chroms to the reference chroms (for files originally compressed with REF_EXTERNAL/EXT_STORE and have alternative chroms)
        chrom_2ref_load (ref); 

        ref_contigs_load_contigs (ref); // note: in case of REF_EXTERNAL, reference is already pre-loaded

        // read dict_id aliases, if there are any
        dict_id_read_aliases();
    }

    // if the user wants to see only the header, we can skip regions and random access
    if (!flag.header_only) {

        ctx_read_all_counts(); // read all counts section

        // update chrom node indices using the CHROM dictionary, for the user-specified regions (in case -r/-R were specified)
        if (flag.regions && ZCTX(flag.luft ? DTFZ(luft_chrom) : DTFZ(prim_chrom))->word_list.len) // FASTQ compressed without reference, has no CHROMs 
            regions_make_chregs (ZCTX(flag.luft ? DTFZ(luft_chrom) : DTFZ(prim_chrom)));

        // if the regions are negative, transform them to the positive complement instead
        regions_transform_negative_to_positive_complement();

        // if this is a stored reference we load the reference random access that will determined which reference sections
        // should be read & uncompressed in case of --regions.
        // note: in case of a data file with stored reference - SEC_REF_RAND_ACC will contain the random access of the reference
        // and SEC_RANDOM_ACCESS will contain the random access of the data. In case of a .ref.genozip file, both sections exist 
        // and are identical. It made the coding easier and their size is negligible.
        random_access_load_ra_section (SEC_RANDOM_ACCESS, flag.luft ? DTFZ(luft_chrom) : DTFZ(prim_chrom), &z_file->ra_buf, "z_file->ra_buf", 
                                       !flag.show_index ? NULL : flag.luft ? RA_MSG_LUFT : RA_MSG_PRIM);

        random_access_load_ra_section (SEC_REF_RAND_ACC, CHROM, ref_get_stored_ra (ref), "ref_stored_ra", 
                                       flag.show_ref_index && !flag.reading_reference ? RA_MSG_REF : NULL);

        if ((flag.reference & REF_ZIP_CHROM2REF) && !flag.reading_reference && !flag.genocat_no_reconstruct)
            // xxx is this actually used? 
            chrom_index_by_name (CHROM); // create alphabetically sorted index for user file (not reference) chrom word list

        // case: reading reference file
        if (flag.reading_reference) {

            // when reading the reference for genocat --sex/coverage/idxstats, don't need the actual REF sections 
            if (is_genocat && (flag.show_sex || flag.show_coverage || flag.idxstats)) 
                goto done;  

            bool dispatcher_invoked = false;

            // attempt to mmap a cached reference, and if one doesn't exist, uncompress the reference file and cache it
            if (!flag.genocat_no_ref_file) {
                if (!ref_mmap_cached_reference (ref)) {
                    ref_load_stored_reference (ref);
                    ref_generate_reverse_complement_genome (ref);

                    // start creating the genome cache now in a background thread, but only if we loaded the entire reference
                    if (!flag.regions) ref_create_cache_in_background (ref); 

                    dispatcher_invoked = true;
                }
            }

            // load the IUPACs list of the reference (rare non-ACGT "bases")
            ref_iupacs_load (ref);

            // load the refhash, if we are compressing FASTA or FASTQ, or if user requested to see it
            if (  (primary_command == ZIP && flag.aligner_available) ||
                  (flag.show_ref_hash && is_genocat)) {
                
                refhash_initialize (&dispatcher_invoked);
            }

            // exit now if all we wanted was just to see the reference (we've already shown it)
            if ((flag.show_reference || flag.show_is_set || flag.show_ref_hash) && is_genocat) exit_ok();

            if (dispatcher_invoked) progress_finalize_component ("Done");
        }

        // case: non-reference file has stored reference sections
        else if (has_ref_sections) {
            if (!flag.genocat_no_ref_file) { 
                ref_load_stored_reference (gref);

                // exit now if all we wanted was just to see the reference (we've already shown it)
                if ((flag.show_reference || flag.show_is_set || flag.show_ref_hash) && is_genocat) exit_ok();
            }
        }
    }
    
done:
    file_seek (z_file, 0, SEEK_SET, false);

    COPY_TIMER_VB (evb, piz_read_global_area);

    return z_file->data_type;
}

// main thread
bool piz_read_one_vb (VBlockP vb, bool for_reconstruction)
{
    START_TIMER; 
   
    Section sec = sections_vb_header (vb->vblock_i, false); 
    
    int32_t vb_header_offset = zfile_read_section (z_file, vb, vb->vblock_i, &vb->z_data, "z_data", SEC_VB_HEADER, sec++); 
    ASSERT0 (vb_header_offset >= 0, "Unexpectedly VB_HEADER section was skipped");

    SectionHeaderVbHeader header = *(SectionHeaderVbHeader *)Bc (vb->z_data, vb_header_offset); // copy of header as it will be overwritten in piz_read_all_ctxs

    // any of these might be overridden by callback
    vb->flags            = header.h.flags.vb_header;
    vb->recon_size       = BGEN32 (header.recon_size_prim);   // might be modified by callback (if DVCF Luft)
    vb->longest_line_len = BGEN32 (header.longest_line_len);
    vb->expected_digest  = header.digest;
    vb->chrom_node_index = WORD_INDEX_NONE;
    vb->lines.len        = VER(14) ? (sec-1)->num_lines : BGEN32 (header.v13_top_level_repeats);
    vb->comp_i           = (sec-1)->comp_i; 
    vb->show_containers  = (flag.show_containers == SHOW_CONTAINERS_ALL_VBs || flag.show_containers == vb->vblock_i); // a per-VB value bc in SAM Load-Prim VBs =false vs normal VBs have the flag value (set in sam_piz_dispatch_one_load_sag_vb)

    if (txt_file) { // sometimes we don't have a txtfile, eg when genocat is used with some flags that emit other data, no the file
        vb->vb_position_txt_file = txt_file->txt_data_so_far_single_0; // position in original txt file (before any ZIP or PIZ modifications)
        txt_file->num_lines += vb->lines.len; // source file lines
    }

    uint32_t txt_data_so_far_single_0_increment = BGEN32 (header.recon_size_prim); // might be modified by callback

    // in case of unbind, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every txt component
    if (flag.unbind) vb->vblock_i = BGEN32 (header.h.vblock_i);

    if (flag.show_vblocks) 
        iprintf ("READING(id=%d) vb=%s num_lines=%u recon_size=%u genozip_size=%u longest_line_len=%u\n",
                 vb->id, VB_NAME, vb->lines.len32, vb->recon_size, BGEN32 (header.z_data_bytes), vb->longest_line_len);

    ctx_overlay_dictionaries_to_vb (VB); /* overlay all dictionaries (not just those that have fragments in this vblock) to the vb */ 

    buf_alloc (vb, &vb->z_section_headers, 0, MAX_DICTS * 2 + 50, uint32_t, 0, "z_section_headers"); // room for section headers  

    BNXT32 (vb->z_section_headers) = vb_header_offset; // vb_header_offset is always 0 for VB header

    // read all b250 and local of all fields and subfields
    piz_read_all_ctxs (vb, &sec, false);

    bool ok_to_compute = DTPZ(piz_init_vb) ? DTPZ(piz_init_vb)(vb, &header, &txt_data_so_far_single_0_increment) : true;

    vb->translation = dt_get_translation (vb); // must be after piz_init_vb, as in VCF we set vb->vb_chords there, needed for dt_get_translation

    if (txt_file) 
        txt_file->txt_data_so_far_single_0 += txt_data_so_far_single_0_increment;

    if (ok_to_compute && for_reconstruction && flag.collect_coverage) 
        coverage_initialize (vb);

    COPY_TIMER (piz_read_one_vb); 

    return ok_to_compute;
}

static Digest piz_one_verify_digest (void)
{
    char s[200];  

    // since v14, if Alder32, we verify each VB, but we don't create a cumulative digest for the entire file. 
    // if we reached here, txt header and all VBs are verified.
    if (VER(14) && z_file->z_flags.adler) {
        if (flag.test) { 
            sprintf (s, "verified as identical to the original %s", dt_name (txt_file->data_type));
            progress_finalize_component (s); 
        }

        return DIGEST_NONE; // we can't calculate the digest for some reason
    }

    if (v8_digest_is_zero (z_file->digest) || digest_is_zero (z_file->digest) || 
        flag.genocat_no_reconstruct || flag.data_modified || flag_loading_auxiliary) 
        return DIGEST_NONE; // we can't calculate the digest for some reason

    Digest decompressed_file_digest = digest_snapshot (&z_file->digest_ctx, "file"); 

    if (digest_recon_is_equal (decompressed_file_digest, z_file->digest)) {
        if (flag.test) { 
            sprintf (s, "verified as identical to the original %s (%s=%s)", 
                     dt_name (txt_file->data_type), digest_name(), digest_display (decompressed_file_digest).s);
            progress_finalize_component (s); 
        }
    }
    
    else if (flag.test) {
        progress_finalize_component ("FAILED!");
        ABORT ("Error: %s of original file=%s is different than decompressed file=%s\n",
               digest_name(), digest_display (z_file->digest).s, digest_display (decompressed_file_digest).s);
    }

    // if decompressed incorrectly - warn, but still give user access to the decompressed file
    else if (!digest_is_zero (z_file->digest)) { // its ok if we decompressed only a partial file
        piz_digest_failed = true; // inspected by main_genounzip
        WARN ("File integrity error: %s of decompressed file %s is %s, but %s of the original %s file was %s", 
              digest_name(), txt_file->name, digest_display (decompressed_file_digest).s, digest_name(), 
              dt_name (txt_file->data_type), digest_display (z_file->digest).s);
    }

    return decompressed_file_digest;
}

static void piz_handover_or_discard_vb (Dispatcher dispatcher, VBlockP *vb)
{
    // add up context decompress time
    if (flag.show_time)
        ctx_add_compressor_time_to_zf_ctx (*vb);

    if (!flag.no_writer_thread && !(*vb)->preprocessing) { // note: in SAM with gencomp - writer does the digest calculation
        writer_handover_data (vb);
        dispatcher_recycle_vbs (dispatcher, false); // don't release VB- it will be released in writer_release_vb when writing is completed
    }
    else {
        if ((*vb)->preprocessing) 
            DT_FUNC (z_file, piz_after_preproc)(*vb);

        dispatcher_recycle_vbs (dispatcher, true);  // also release VB
    }
}

// returns false if VB was dispatched, and true if vb was skipped
static void piz_dispatch_one_vb (Dispatcher dispatcher, Section sec)
{
    VBlockP next_vb = dispatcher_generate_next_vb (dispatcher, sec->vblock_i, sec->comp_i);

    // read one VB's data from z_file
    bool reconstruct = piz_read_one_vb (next_vb, true)  // read even if no_reconstruct
                    && !flag.genocat_no_reconstruct; 

    if (reconstruct) {
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

// main thread: called in order of VBs
static void piz_handle_reconstructed_vb (Dispatcher dispatcher, VBlockP vb, uint64_t *num_nondrop_lines)
{
    ASSERTW (vb->txt_data.len == vb->recon_size || flag.data_modified, // files are the same size, unless we intended to modify the data
            "Warning: vblock_i=%s/%u (num_lines=%u) had %s bytes in the original %s file but %s bytes in the reconstructed file (diff=%d)", 
            comp_name (vb->comp_i), vb->vblock_i, vb->lines.len32, str_int_commas (vb->recon_size).s, dt_name (txt_file->data_type), 
            str_int_commas (vb->txt_data.len).s, 
            (int32_t)vb->txt_data.len - (int32_t)vb->recon_size);

    *num_nondrop_lines += vb->num_nondrop_lines;
    if (flag.count == COUNT_VBs)
        iprintf ("vb=%s lines=%u nondropped_lines=%u txt_data.len=%u\n", 
                  VB_NAME, vb->lines.len32, vb->num_nondrop_lines, vb->txt_data.len32);

    if (flag.collect_coverage)    
        coverage_add_one_vb (vb);
    
    else if (flag.reading_kraken)
        kraken_piz_handover_data (vb);

    z_file->txt_data_so_far_single += vb->recon_size; 

    piz_handover_or_discard_vb (dispatcher, &vb);
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

    writer_create_plan();

    if (flag.test || flag.md5) 
        ASSINP0 (dt_get_translation(NULL).is_src_dt, "Error: --test or --md5 cannot be used when converting a file to another format"); 
    
    
    uint64_t target_progress = (Z_DT(SAM) ? sections_get_num_vbs (SAM_COMP_PRIM) : 0) // VBs pre-processed
                             + 3 * z_file->num_vbs; // VBs read, reconstructed, written

    Dispatcher dispatcher = dispatcher_init (flag.reading_chain     ? "piz-chain"
                                            :flag.reading_reference ? "piz-ref"
                                            :flag.reading_kraken    ? "piz-kraken"
                                            :flag.preprocessing     ? "preprocessing"
                                            :                         PIZ_TASK_NAME, // also referred to in dispatcher_recycle_vbs()
                                             POOL_MAIN,
                                             flag.xthreads ? 1 : global_max_threads, 0, 
                                             flag.test && flag.no_writer_thread, // out-of-order if --test with no writer thread (note: SAM gencomp always have writer thread to do digest). 
                                             flag.test,
                                             z_file->basename, target_progress, 0);

    return dispatcher;
}

// main thread: called once per txt_file created: i.e. once, except if unbinding a paired FASTQ.
// returns true if piz completed, false if piz aborted by piz_initialize
bool piz_one_txt_file (Dispatcher dispatcher, bool is_first_z_file, bool is_last_z_file,
                       CompIType unbind_comp_i) // COMP_NONE unless flag.unbind
{
    dispatcher_start_wallclock();
    
    recon_stack_initialize();
    
    if (DTPZ(piz_initialize) && !DTPZ(piz_initialize)())
        return false; // abort PIZ if piz_initialize says so
      
    bool header_only_file = true; // initialize - true until we encounter a VB header
    uint64_t num_nondrop_lines = 0;

    Section sec = unbind_comp_i != COMP_NONE ? sections_one_before (sections_get_comp_txt_header_sec (unbind_comp_i)) : NULL;

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

            // case SEC_TXT_HEADER
            if (found_header && sec->st == SEC_TXT_HEADER && (unbind_comp_i==COMP_NONE || unbind_comp_i==sec->comp_i)) { 

                if (!writer_does_txtheader_need_recon (sec)) continue;

                if (sec->vblock_i >= 2) continue; // fragments >= 2 where already handled together with the first fragment
                
                // note: also starts writer, and if unbinding, also opens the txt file and hands data over to the writer
                txtheader_piz_read_and_reconstruct (sec); 

                // case --unbind: unpausing after previous txt_file pause (requires txt file to be open)
                if (flag.unbind) dispatcher_resume (dispatcher);  
            }

            // case SEC_VB_HEADER
            else if (found_header && sec->st == SEC_VB_HEADER && (unbind_comp_i==COMP_NONE || unbind_comp_i==sec->comp_i)) {
                
                if (!writer_does_vb_need_recon (sec->vblock_i)) {
                    dispatcher_increment_progress ("vb_no_recon", 3); // skipped reading, reconstructing, writing

                    if (flag.show_vblocks) 
                        iprintf ("SKIPPED vb=%s/%u\n", comp_name (sec->comp_i), sec->vblock_i); 
                    continue;
                }

                piz_dispatch_one_vb (dispatcher, sec);
                header_only_file = false;
            }

            // case: we're done with this txt_file (either no header bc EOF, or TXT_HEADER belongs to the next txt_file when unbinding)
            else {
                if (flag.show_vblocks) 
                    iprintf ("INPUT EXHAUSTED - no more SEC_VB_HEADER or SEC_TXT_HEADER for txt_file_i=%u\n", z_file->num_txts_so_far);                

                dispatcher_set_no_data_available (dispatcher, false, DATA_EXHAUSTED);

                if (header_only_file)
                    dispatcher_recycle_vbs (dispatcher, true); // note: this is normally done in piz_handover_or_discard_vb
            }
        }

        // if the next thread (by sequential order) is ready, handle the reconstructed VB
        VBlockP recon_vb = dispatcher_get_processed_vb (dispatcher, NULL, false);  // non-blocking
        if (recon_vb) {    
            if (flag.show_vblocks) 
                iprintf ("END_OF_COMPUTE(id=%d) vb=%s/%u num_running_compute_threads(after)=%u\n", 
                         recon_vb->id, comp_name(recon_vb->comp_i), recon_vb->vblock_i, dispatcher_get_num_running_compute_threads(dispatcher));
            
            dispatcher_increment_progress ("preproc_or_recon", 1 + (!recon_vb->preprocessing && flag.no_writer_thread)); // done preprocessing or reconstructing (+1 if skipping writing)
            if (!recon_vb->preprocessing) 
                piz_handle_reconstructed_vb (dispatcher, recon_vb, &num_nondrop_lines);
            else
                piz_handover_or_discard_vb (dispatcher, &recon_vb);
        }

        if (!achieved_something) usleep (30000); // nothing for us to do right now - wait 30ms
    }

    if (flag.show_vblocks) 
        iprintf ("DISPATCHER is done for txt_file_i=%u: %s\n", z_file->num_txts_so_far, txt_file ? txt_file->name : "(no filename)");                

    z_file->num_txts_so_far++;

    // finish writing the txt_file (note: the writer thread also calculates digest in SAM/BAM with PRIM/DEPN)
    writer_finish_writing (z_file->num_txts_so_far == (flag.unbind ? 2 : 1));

    // verifies reconstructed file against MD5 (if compressed with --md5 or --test) or Adler2 and/or codec_args (if bgzf)
    Digest decompressed_file_digest = piz_one_verify_digest();

    if (!flag.test) 
        progress_finalize_component_time ("Done", decompressed_file_digest);

    // --sex and --coverage - output results
    if (txt_file && !flag_loading_auxiliary) {
        if (flag.show_coverage) coverage_show_coverage();
        if (flag.show_sex) coverage_sex_classifier (is_first_z_file);
        if (flag.idxstats) coverage_show_idxstats();
        if (flag.count == CNT_TOTAL) iprintf ("%"PRIu64"\n", num_nondrop_lines);
    }

    if (z_file->num_txts_so_far == (flag.unbind ? 2 : 1)) 
        dispatcher_finish (&dispatcher, NULL, !is_last_z_file || flag.test,
                           flag.show_memory && is_last_z_file);
    else                  
        dispatcher_pause (dispatcher); // we're unbinding and still have more txt_files
     
    DT_FUNC (z_file, piz_finalize)();

    if (flag.show_vblocks) 
        iprintf ("Finished PIZ of %s\n", txt_file ? txt_file->name : "(no filename)");                

    return true;
}
