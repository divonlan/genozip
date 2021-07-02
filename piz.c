// ------------------------------------------------------------------
//   piz.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "txtfile.h"
#include "vblock.h"
#include "base250.h"
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

bool piz_grep_match (const char *start, const char *after)
{
    bool found = false;
    SAFE_NUL (after);

    if (!flag.grepw) {
        found = !!strstr (start, flag.grep);
        goto done;
    }

    // case: --grepw - grep whole word
    const char *s = start;
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

// called by main thread in FASTA and FASTQ, in case of --grep, to decompress and reconstruct the desc line, to 
// see if this vb is included. 
bool piz_test_grep (VBlock *vb)
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
    vb->first_line       = 1 + txt_file->num_lines; // doesn't count a dropped txtheader
    vb->recon_num_lines  = BGEN32 (header->num_lines_prim);
    vb->recon_size       = BGEN32 (header->recon_size_prim);
    vb->longest_line_len = BGEN32 (header->longest_line_len);

    // in case of unbind, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every sam component
    if (flag.unbind) vb->vblock_i = BGEN32 (header->h.vblock_i);
    
    // we only need room for one line for now 
    buf_alloc (vb, &vb->txt_data, 0, vb->longest_line_len, char, 1.1, "txt_data");

    // uncompress & map desc field (filtered by piz_is_skip_section)
    vb->grep_stages = GS_TEST; // tell piz_is_skip_section to skip decompressing sections not needed for determining the grep
    piz_uncompress_all_ctxs ((VBlockP)vb, 0);
    vb->grep_stages = GS_UNCOMPRESS; // during uncompress in the compute thread, uncompress only what was not already uncompressed here

    // reconstruct each description line and check for string matching with flag.grep
    bool found = false, match = false;

    Context *desc_ctx =  CTX(vb->data_type == DT_FASTQ ? FASTQ_DESC : FASTA_DESC);
    desc_ctx->iterator.next_b250 = FIRSTENT (uint8_t, desc_ctx->b250); 

    uint32_t num_descs; // number of desciption lines in this VB
    if (vb->data_type == DT_FASTQ) {
        vb->line_i = 4 * vb->first_line;
        num_descs = (uint32_t)vb->lines.len; // every read has a description line
    }
    else {
        vb->line_i = vb->first_line;
        num_descs = random_access_num_chroms_start_in_this_vb (vb->vblock_i);
    }

    // iterate on all DESCs of this VB
    for (uint32_t desc_i=0; desc_i < num_descs; desc_i++) { 

        reconstruct_from_ctx (vb, desc_ctx->did_i, 0, true);

        match = flag.grep && piz_grep_match (FIRSTENT (char, vb->txt_data), AFTERENT (char, vb->txt_data));

        vb->txt_data.len = 0; // reset

        if (match) { 
            found = true; // we've found a match to the grepped string
            if (vb->data_type == DT_FASTQ) break; // for FASTA, we need to go until the last line, for FASTQ, we can break here
        }

        if (vb->data_type == DT_FASTQ) vb->line_i += 4; // note: for FASTA we have no idea what txt line we're on, because we're only tracking DESC lines
    }

    // last FASTA - carry over whether its grepped to the next VB - in case next VB starts not from the description line
    // similarly, note whether the previous VB ended with a grepped sequence. If previous VB didn't have any description
    // i.e the entire VB was a sequence that started in an earlier VB - the grep status of the easier VB is carried forward
    if (vb->data_type == DT_FASTA) 
        // if the last contig of the previous vb was grepped in - then include this VB anyway
        found = fasta_piz_initialize_contig_grepped_out (vb, desc_ctx->b250.len > 0, match) || found;

    // reset iterators - piz_fast*_reconstruct_vb will use them again 
    ctx_init_iterator (desc_ctx);
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) {
        
        Context *ctx = CTX(did_i);
        if (dict_id_is_type_1 (ctx->dict_id)) {
            ctx_init_iterator (ctx);
            ctx->last_delta = ctx->last_value.f = 0;
        }
    }

    return found; 
}

bool piz_default_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (!vb) return false; // we don't skip reading any SEC_DICT / SEC_COUNTS sections

    // B250, LOCAL, COUNT sections
    bool skip = exe_type == EXE_GENOCAT && dict_id.num && dict_id.num != dict_id_fields[CHROM] && (!flag.luft || dict_id.num != dict_id_fields[ODID(oCHROM)]) && (
    
    // sometimes we don't need dictionaries. but we always load CHROM.
        (flag.genocat_no_dicts && dict_id_typeless (dict_id).num != flag.show_one_counts.num)

    // if show_counts - we only need the requested section and CHROM (note: not true for dump_one_b250_dict_id,
    // as we need to reconstruct to dump it)
    ||  (flag.show_one_counts.num && dict_id_typeless (dict_id).num != flag.show_one_counts.num)

    // if --counts, we filter here - TOPLEVEL only - unless there's a skip_section function which will do the filtering
    ||  (flag.count && !DTPZ(is_skip_section) && dict_id.num != dict_id_fields[DTFZ(toplevel)]) 
    );

    skip |= flag.genocat_no_ref_file && (st == SEC_REFERENCE || st == SEC_REF_HASH || st == SEC_REF_IS_SET);

    return skip;
}

// PIZ compute thread: decompress all contexts
// ZIP compute thread in FASTQ: decompress pair_1 contexts when compressing pair_2
uint32_t piz_uncompress_all_ctxs (VBlock *vb, 
                                  uint32_t pair_vb_i) // used in ZIP when uncompressing previous file's paired sections
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);
    
    uint64_t section_i = pair_vb_i ? 0 : 1; // normally, we skip the VB header (starting from 1), but when uncompressing paired sections there is no VB header
    for ( ; section_i < vb->z_section_headers.len; section_i++) {

        SectionHeaderCtx *header = (SectionHeaderCtx *)ENT (char, vb->z_data, section_index[section_i]);

        bool is_local = header->h.section_type == SEC_LOCAL;
        bool is_b250  = header->h.section_type == SEC_B250;
        if (!is_b250 && !is_local) break;

        Context *ctx = ctx_get_ctx (vb, header->dict_id); // gets the context (creating it if it doesn't already exist)

        // case: in ZIP & we are PAIR_2, reading PAIR_1 info - save flags, ltype, lcodec to restore later
        //struct FlagsCtx save_flags={}; LocalType save_ltype=0; Codec save_lcodec=0;
        //if (pair_vb_i && command == ZIP) { save_flags=ctx->flags; save_ltype=ctx->ltype; save_lcodec=ctx->lcodec; }
        
        bool is_pair_section = (BGEN32 (header->h.vblock_i) == pair_vb_i); // is this a section of "pair 1" 

        if (!is_pair_section) {
            ctx->flags = header->h.flags.ctx; // overrides default inherited from vb_i=1 (assigned in piz_read_all_ctxs)
            ctx->ltype = header->ltype;
            if (is_local) ctx->lcodec = header->h.codec;
        }
        else {
            // overcome bug in ZIP --pair - local sections with junk data created in addition to the expected b250.
            // we ignore these local sections (or allow b250 to overwrite them)
            if (is_local && ctx->pair_b250) continue;

            ctx->pair_flags = header->h.flags.ctx;            
            ctx->pair_b250  = is_b250;                
            ctx->pair_local = is_local;
        }

        // case: in PIZ: flags.paired appears on the sections the "pair 2" VB (that come first in section_index)
        if (ctx->flags.paired && !pair_vb_i) 
            pair_vb_i = fastq_get_pair_vb_i (vb);

        // initialize b250 iterator
        if (is_b250) {
            if (is_pair_section) ctx->pair_b250_iter = (SnipIterator){ .next_b250 = FIRSTENT (uint8_t, ctx->pair), .prev_word_index = WORD_INDEX_NONE };
            else                 ctx->iterator       = (SnipIterator){ .next_b250 = FIRSTENT (uint8_t, ctx->b250), .prev_word_index = WORD_INDEX_NONE };
        }

        Buffer *target_buf = is_local ? &ctx->local : &ctx->b250;

        zfile_uncompress_section (vb, header, 
                                  is_pair_section ? &ctx->pair      : target_buf, 
                                  is_pair_section ? "context->pair" : is_local ? "contexts->local" : "contexts->b250", 
                                  is_pair_section ? pair_vb_i : vb->vblock_i,
                                  header->h.section_type); 

        if (!is_pair_section && is_local && dict_id_typeless (ctx->dict_id).num == flag.dump_one_local_dict_id.num) 
            ctx_dump_binary (vb, ctx, true);

        if (!is_pair_section && !is_local && dict_id_typeless (ctx->dict_id).num == flag.dump_one_b250_dict_id.num) 
            ctx_dump_binary (vb, ctx, false);

#           define adjust_lens(buf) { \
            buf.len /= lt_desc[ctx->ltype].width; \
            if (ctx->ltype == LT_BITMAP) { \
                buf.param = buf.len * 64 - header->param ; /* number of bits */ \
                LTEN_bit_array (buf_get_bitarray (&buf)); \
            } \
            else if (lt_desc[ctx->ltype].file_to_native)   \
                lt_desc[ctx->ltype].file_to_native (&buf, &ctx->ltype); \
        }

        if      (is_pair_section) adjust_lens (ctx->pair)
        else if (is_local)        adjust_lens (ctx->local)

        if (header->h.flags.ctx.copy_param)
            target_buf->param = header->param;

//skip_uncompress_due_to_old_bug:
        // restore
//        if (pair_vb_i && command == ZIP) { ctx->flags=save_flags; ctx->ltype=save_ltype; ctx->lcodec=save_lcodec; }
    }

    // initialize pair iterators (pairs only exist in fastq)
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) {
        Context *ctx = CTX(did_i);
        if (buf_is_alloc (&ctx->pair))
            ctx->pair_b250_iter = (SnipIterator){ .next_b250 = FIRSTENT (uint8_t, ctx->pair),
                                                  .prev_word_index = -1 };
    }

    ctx_map_aliases (vb);

    return section_i;
}

// PIZ compute thread entry point
static void piz_reconstruct_one_vb (VBlock *vb)
{
    START_TIMER;

    ASSERTNOTNULL (vb);
    ASSERT (vb->vblock_i, "vb->vblock_i is 0: vb->compute_thread_id=%d pthread=%"PRIu64, 
            vb->compute_thread_id, (uint64_t)pthread_self());

    ASSERT (!flag.reference || ref_is_loaded (gref) || flag.genocat_no_ref_file,
            "reference is not loaded correctly (vb=%u)", vb->vblock_i);

    // note: txt_data is fully allocated in advance and cannot be extended mid-reconstruction (container_reconstruct_do and possibly others rely on this)
    #define OVERFLOW_SIZE 65536 // allow some overflow space as sometimes we reconstruct unaccounted for data: 1. container templates 2. reconstruct_peek and others
    buf_alloc (vb, &vb->txt_data, 0, vb->recon_size * vb->translation.factor + OVERFLOW_SIZE, char, 1.1, "txt_data"); 
    
    piz_uncompress_all_ctxs (vb, 0);

    // reconstruct from top level snip
    reconstruct_from_ctx (vb, vb->translation.toplevel, 0, true);

    // compress txt_data into BGZF blocks (in vb->compressed) if applicable
    if (txt_file && txt_file->codec == CODEC_BGZF && !flag.no_writer &&
        !flag.maybe_vb_modified_by_writer)  // if --downsample, --interleave or sorting - writer will BGZF-compress
        bgzf_compress_vb (vb);

    // calculate the digest contribution of this VB to the single file and bound files, and the digest snapshot of this VB
    if (!v8_digest_is_zero (vb->digest_so_far) && !flag.data_modified && !flag.reading_chain) 
        digest_one_vb (vb); // LOOKING FOR A DEADLOCK BUG? CHECK HERE

    vb->is_processed = true; /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ 
    COPY_TIMER (compute);
}

static void piz_read_all_ctxs (VBlock *vb, Section *next_sl)
{
    // ctxs that have dictionaries are already initialized, but others (eg local data only) are not
    ctx_initialize_primary_field_ctxs (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_contexts);

    // ctx.flags defaults to vb_i=1 flags, overridden if a b250 or local section is read. this will not be overridden if all_the_same, i.e. no b250/local sections.
    for (Section sec = sections_vb_first (1, false) + 1; sec->st == SEC_B250 || sec->st == SEC_LOCAL; sec++) {
        ContextP ctx = ctx_get_existing_ctx (vb, sec->dict_id); // will exist if it has a dict (all_the_same sections always have a dict)
        if (ctx) ctx->flags = sec->flags.ctx;
    }

    while ((*next_sl)->st == SEC_B250 || (*next_sl)->st == SEC_LOCAL) {
        uint32_t section_start = vb->z_data.len;
        *ENT (uint32_t, vb->z_section_headers, vb->z_section_headers.len) = section_start; 

        // create a context even if section is skipped, for containers to work (skipping a section should be mirrored in a container filter)
        ctx_get_ctx_do (z_file->contexts, z_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_contexts, (*next_sl)->dict_id);
        int32_t offset = zfile_read_section (z_file, vb, vb->vblock_i, &vb->z_data, "z_data", (*next_sl)->st, *next_sl); // returns 0 if section is skipped
        if (offset != SECTION_SKIPPED) vb->z_section_headers.len++;
        
        (*next_sl)++;                             
    }
}

// Called by PIZ main thread: read all the sections at the end of the file, before starting to process VBs
DataType piz_read_global_area (Reference ref)
{
    bool success = zfile_read_genozip_header (0);

    if (flag.show_stats) stats_read_and_display();

    if (!success) return DT_NONE;

    dict_id_initialize (z_file->data_type); // must run after zfile_read_genozip_header that sets z_file->data_type
    
    // check if the genozip file includes a reference
    bool has_ref_sections = !!sections_last_sec (SEC_REFERENCE, true);

    ASSERTW (!has_ref_sections || flag.reference != REF_EXTERNAL || flag.reading_reference, 
             "%s: ignoring reference file %s because it was not compressed with --reference", z_name, ref_get_filename (ref));

    if (!flag.reading_reference && has_ref_sections) {
        ref_destroy_reference (ref, false);     // destroy an old reference, if one is loaded
        flag.reference = REF_STORED; // possibly override REF_EXTERNAL (it will be restored for the next file in )
    }

    // read all dictionaries - CHROM/RNAME is needed for regions_make_chregs(). 
    // Note: some dictionaries are skipped based on skip() and all flag logic should implemented there
    ctx_read_all_dictionaries(); 

//    if (!flag.header_only && !z_dual_coords) { // dual coordinates need this stuff of for the rejects part of the header
    if (!flag.header_only || z_dual_coords) { // dual coordinates need this stuff of for the rejects part of the header

        // mapping of the file's chroms to the reference chroms (for files originally compressed with REF_EXTERNAL/EXT_STORE and have alternative chroms)
        ref_alt_chroms_load (ref); 

        ref_contigs_load_contigs (ref); // note: in case of REF_EXTERNAL, reference is already pre-loaded

        // read dict_id aliases, if there are any
        dict_id_read_aliases();
    }

    // reading external references for lifting with --chain - load the IUPACs list of the reference (rare non-ACTG "bases")
    if ((flag.reading_chain && flag.reading_reference) || flag.show_ref_iupacs)
        ref_iupacs_load(ref);

    // if the user wants to see only the header, we can skip regions and random access
    if (!flag.header_only) {

        ctx_read_all_counts(); // read all counts section

        // update chrom node indices using the CHROM dictionary, for the user-specified regions (in case -r/-R were specified)
        if (flag.regions && ZCTX(flag.luft ? ODID(oCHROM) : CHROM)->word_list.len) // FASTQ compressed without reference, has no CHROMs 
            regions_make_chregs (&z_file->contexts[flag.luft ? ODID(oCHROM) : CHROM]);

        // if the regions are negative, transform them to the positive complement instead
        regions_transform_negative_to_positive_complement();

        // if this is a stored reference we load the reference random access that will determined which reference sections
        // should be read & uncompressed in case of --regions.
        // note: in case of a data file with stored reference - SEC_REF_RAND_ACC will contain the random access of the reference
        // and SEC_RANDOM_ACCESS will contain the random access of the data. In case of a .ref.genozip file, both sections exist 
        // and are identical. It made the coding easier and their size is negligible.
        random_access_load_ra_section (SEC_RANDOM_ACCESS, flag.luft ? ODID(oCHROM) : CHROM, &z_file->ra_buf, "z_file->ra_buf", 
                                       !flag.show_index ? NULL : flag.luft ? RA_MSG_LUFT : RA_MSG_PRIM);

        random_access_load_ra_section (SEC_REF_RAND_ACC, CHROM, ref_get_stored_ra (ref), "ref_stored_ra", 
                                       flag.show_ref_index && !flag.reading_reference ? RA_MSG_REF : NULL);

        if ((flag.reference == REF_STORED || flag.reference == REF_EXTERNAL || flag.reference == REF_LIFTOVER) && 
            !flag.reading_reference && !flag.genocat_no_reconstruct)
            ref_contigs_sort_chroms(); // create alphabetically sorted index for user file (not reference) chrom word list

        // case: reading reference file
        if (flag.reading_reference) {

            // when reading the reference for genocat --show-sex/coverage/idxstats, don't need the actual REF sections 
            if (exe_type == EXE_GENOCAT && (flag.show_sex || flag.show_coverage || flag.idxstats)) 
                goto done;  

            bool dispatcher_invoked = false;

            // attempt to mmap a cached reference, and if one doesn't exist, uncompress the reference file and cache it
            if (!ref_mmap_cached_reference (ref)) {
                ref_load_stored_reference (ref);
                ref_generate_reverse_complement_genome (ref);

                // start creating the genome cache now in a background thread, but only if we loaded the entire reference
                if (!flag.regions) ref_create_cache_in_background (ref); 

                dispatcher_invoked = true;
            }

            // load the refhash, if we are compressing FASTA or FASTQ, or if user requested to see it
            if (  (primary_command == ZIP && flag.ref_use_aligner) ||
                  (flag.show_ref_hash && exe_type == EXE_GENOCAT)) {
                
                refhash_initialize (&dispatcher_invoked);
            }

            // exit now if all we wanted was just to see the reference (we've already shown it)
            if ((flag.show_reference || flag.show_is_set || flag.show_ref_hash) && exe_type == EXE_GENOCAT) exit_ok;

            if (dispatcher_invoked) progress_finalize_component ("Done");
        }

        // case: non-reference file has stored reference sections
        else if (has_ref_sections) {
            if (!flag.genocat_no_ref_file) { 
                ref_load_stored_reference (gref);

                // exit now if all we wanted was just to see the reference (we've already shown it)
                if ((flag.show_reference || flag.show_is_set || flag.show_ref_hash) && exe_type == EXE_GENOCAT) exit_ok;
            }
            else
                z_file->disk_size_minus_skips -= sections_get_ref_size();    
        }
    }
    
done:
    file_seek (z_file, 0, SEEK_SET, false);

    return z_file->data_type;
}

static bool piz_read_one_vb (VBlock *vb)
{
    START_TIMER; 

    if (flag.lines_last >= 0 && txt_file->num_lines > flag.lines_last) // lines_last is 0-based
        return false; // we don't need this VB as we have read all the data needed according to --lines

    Section sl = sections_vb_first (vb->vblock_i, false); 

    vb->vb_position_txt_file = txt_file->txt_data_so_far_single_0; // position in original txt file (before any ZIP or PIZ modifications)
    
    int32_t vb_header_offset = zfile_read_section (z_file, vb, vb->vblock_i, &vb->z_data, "z_data", SEC_VB_HEADER, sl++); 

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)ENT (char, vb->z_data, vb_header_offset);

    // calculate the coordinates in which this VB will be rendered - PRIMARY or LUFT
    vb->vb_coords        = !z_dual_coords ? DC_PRIMARY // non dual-coordinates file - always PRIMARY
                         : header->h.flags.vb_header.coords == DC_PRIMARY ? DC_PRIMARY // reject component ##primary_only
                         : header->h.flags.vb_header.coords == DC_LUFT    ? DC_LUFT    // reject component ##luft_only
                         : flag.luft ? DC_LUFT // dual component - render as LUFT
                         : DC_PRIMARY;         // dual component - render as PRIMARY

    vb->recon_size       = BGEN32 (vb->vb_coords==DC_PRIMARY ? header->recon_size_prim : header->recon_size_luft); 
    vb->first_line       = 1 + txt_file->num_lines; // doesn't count a dropped txtheader
    vb->lines.len        = BGEN32 (header->top_level_repeats);   
    vb->recon_num_lines  = z_file->genozip_version >= 12 ? BGEN32 (vb->vb_coords==DC_PRIMARY ? header->num_lines_prim : header->num_lines_luft) : vb->lines.len;
    vb->longest_line_len = BGEN32 (header->longest_line_len);
    vb->digest_so_far    = header->digest_so_far;
    vb->chrom_node_index = WORD_INDEX_NONE;

    vb->is_rejects_vb    = z_dual_coords && (header->h.flags.vb_header.coords != DC_BOTH);

    vb->translation      = dt_get_translation (vb); // vb->vb_chords needs to be set first

    txt_file->num_lines += vb->lines.len; // source file lines

    // accounting for data as in the original source file 
    txt_file->txt_data_so_far_single_0 += BGEN32 (txt_file->txt_flags.is_txt_luft ? header->recon_size_luft : header->recon_size_prim); 

    // in case of unbind, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every txt component
    if (flag.unbind) vb->vblock_i = BGEN32 (header->h.vblock_i);

    if (flag.show_vblocks) 
        iprintf ("READING(id=%d) vb_i=%u first_line=%"PRIu64" num_lines=%u recon_size=%u genozip_size=%u longest_line_len=%u\n",
                 vb->id, vb->vblock_i, vb->first_line, vb->recon_num_lines, vb->recon_size, BGEN32 (header->z_data_bytes), vb->longest_line_len);

    ctx_overlay_dictionaries_to_vb ((VBlockP)vb); /* overlay all dictionaries (not just those that have fragments in this vblock) to the vb */ 

    buf_alloc (vb, &vb->z_section_headers, 0, MAX_DICTS * 2 + 50, uint32_t, 0, "z_section_headers"); // room for section headers  

    NEXTENT (uint32_t, vb->z_section_headers) = vb_header_offset; // vb_header_offset is always 0 for VB header

    // read all b250 and local of all fields and subfields
    if (txt_file->num_lines > flag.lines_first) // read only if it might be ok_to_compute
        piz_read_all_ctxs (vb, &sl);

    // check some flags (--grep, --lines...)
    bool ok_to_compute = (txt_file->num_lines > flag.lines_first) // --lines: we've reached the start line (note: if we passed lines_end, VB is dropped before calling this function)
                      && (DTPZ(piz_read_one_vb) ? DTPZ(piz_read_one_vb)(vb, sl) : true); // logic specific to this data type (--grep for FASTQ, --grep,--regions for FASTA)

    // calculate the BGZF blocks from SEC_BGZF that the compute thread is expected to re-create,
    // unless isizes was not loaded (bc flag.data_modified).
    if (flag.bgzf == FLAG_BGZF_BY_ZFILE && txt_file->codec == CODEC_BGZF)     
        bgzf_calculate_blocks_one_vb (vb, vb->recon_size); // does nothing if isizes is not loaded

    // initialize coverage counters
    if (flag.collect_coverage)
        coverage_initialize (vb);
            
    COPY_TIMER (piz_read_one_vb); 

    return ok_to_compute;
}

static Digest piz_one_verify_digest (void)
{
    Digest original_digest = flag.unbind ? txt_file->digest /* digest_single */ : z_file->digest /* digest_bound */;

    if (v8_digest_is_zero (original_digest) || digest_is_zero (original_digest) || 
        flag.genocat_no_reconstruct || flag.data_modified || flag_loading_auxiliary) 
        return DIGEST_NONE; // we can't calculate the digest for some reason

    // Note: in piz, we compare txt_file->digest_ctx_bound to original bound or single, depending on flag.unbind
    Digest decompressed_file_digest = digest_finalize (&txt_file->digest_ctx_bound, "file:digest_ctx_bound"); 
    char s[200]; 

    if (digest_recon_is_equal (decompressed_file_digest, original_digest)) {
        if (flag.test) { 
            sprintf (s, "%s = %s verified as identical to the original %s", 
                     digest_name(), digest_display (decompressed_file_digest).s, dt_name (txt_file->data_type));
            progress_finalize_component (s); 
        }
    }
    
    else if (flag.test) {
        progress_finalize_component ("FAILED!!!");
        ABORT ("Error: %s of original file=%s is different than decompressed file=%s\nPlease contact bugs@genozip.com to help fix this bug in genozip\n",
               digest_name(), digest_display (original_digest).s, digest_display (decompressed_file_digest).s);
    }

    // if compressed incorrectly - warn, but still give user access to the decompressed file
    else ASSERTW (digest_is_zero (original_digest), // its ok if we decompressed only a partial file
                  "File integrity error: %s of decompressed file %s is %s, but %s of the original %s file was %s", 
                  digest_name(), txt_file->name, digest_display (decompressed_file_digest).s, digest_name(), 
                  dt_name (txt_file->data_type), digest_display (original_digest).s);

    return decompressed_file_digest;
}

static void piz_handover_or_discard_vb (Dispatcher dispatcher, VBlockP *vb)
{
    if (!flag.no_writer) {
        writer_handover_data (vb);
        dispatcher_recycle_vbs (dispatcher, false); // don't release VB- it will be released in writer_release_vb when writing is completed
    }
    else
        dispatcher_recycle_vbs (dispatcher, true); // also release VB
}

// returns false if VB was dispatched, and true if vb was skipped
static bool piz_dispatch_one_vb (Dispatcher dispatcher, Section sl_ent)
{
    if (writer_is_vb_no_read (sl_ent->vblock_i)) 
        return true; // skip this VB - we don't need it

    VBlock *next_vb = dispatcher_generate_next_vb (dispatcher, sl_ent->vblock_i);
    next_vb->component_i = z_file->num_txt_components_so_far;
    
    // read one VB's genozip data
    bool reconstruct = piz_read_one_vb (next_vb)  // read even if no_reconstruct
                    && !flag.genocat_no_reconstruct; 

    if (reconstruct) 
        dispatcher_compute (dispatcher, piz_reconstruct_one_vb);
    
    // case: we won't proceed to uncompressing, reconstructing we're done reading - just handover
    // an empty VB as it appears in the recon plan, and writer might be already blocking on waiting for it
    else {
        dispatcher_abandon_next_vb (dispatcher); // just moves the to processed_vb so dispatcher_recycle_vbs can recycle it
        piz_handover_or_discard_vb (dispatcher, &next_vb);
    }

    return false;
}

static void piz_handle_reconstructed_vb (Dispatcher dispatcher, VBlock *vb, uint64_t *num_nondrop_lines)
{
    ASSERTW (vb->txt_data.len == vb->recon_size || flag.data_modified, // files are the same size, unless we intended to modify the data
            "Warning: vblock_i=%u (num_lines=%u vb_start_line_in_file=%"PRIu64") had %s bytes in the original %s file but %s bytes in the reconstructed file (diff=%d)", 
            vb->vblock_i, (unsigned)vb->lines.len, vb->first_line, str_uint_commas (vb->recon_size).s, dt_name (txt_file->data_type), 
            str_uint_commas (vb->txt_data.len).s, 
            (int32_t)vb->txt_data.len - (int32_t)vb->recon_size);

    *num_nondrop_lines += vb->num_nondrop_lines;
    if (flag.count == COUNT_VBs)
        iprintf ("vb_i=%u lines=%u nondropped_lines=%u txt_data.len=%u\n", 
                  vb->vblock_i, (unsigned)vb->lines.len, vb->num_nondrop_lines, (unsigned)vb->txt_data.len);

    if (flag.collect_coverage)    
        coverage_add_one_vb (vb);
    
    else if (flag.reading_kraken)
        kraken_piz_handover_data (vb);

    z_file->txt_data_so_far_single += vb->recon_size; 

    piz_handover_or_discard_vb (dispatcher, &vb);
}

Dispatcher piz_z_file_initialize (bool is_last_z_file)
{
    digest_initialize();

    // read genozip header
    DataType data_type = piz_read_global_area (gref);
    if (data_type == DT_NONE || flag.reading_reference) 
        return NULL; // no components in this file (as is always the case with reference files) 

    if (flag.genocat_global_area_only) return NULL;

    writer_create_plan();

    ASSINP (!flag.test || !digest_is_zero (z_file->digest), 
            "Error testing %s: --test cannot be used with this file, as it was not compressed with digest information. See " WEBSITE_DIGEST, z_name);

    if (flag.test || flag.md5) 
        ASSINP0 (dt_get_translation(NULL).is_src_dt, "Error: --test or --md5 cannot be used when converting a file to another format"); 

    Dispatcher dispatcher = dispatcher_init (flag.reading_chain     ? "piz-chain"
                                            :flag.reading_reference ? "piz-ref"
                                            :flag.reading_kraken    ? "piz-kraken"
                                            :                         "piz",
                                             flag.xthreads ? 1 : global_max_threads, 0, flag.test, 
                                             is_last_z_file, true, z_file->basename, PROGRESS_PERCENT, 0);
    return dispatcher;
}

// called once per txt_file created: i.e. if concatenating - a single call, if unbinding there will be multiple calls to this function
// returns true if piz completed, false if piz aborted by piz_initialize
bool piz_one_txt_file (Dispatcher dispatcher, bool is_first_z_file)
{
    bool is_last_txt_file = (z_file->num_txt_components_so_far == z_file->txt_file_info.len-1);
 
    if (DTPZ(piz_initialize) && !DTPZ(piz_initialize)())
        return false; // abort PIZ if piz_initialize says so
      
    bool header_only_file = true; // initialize - true until we encounter a VB header
    uint32_t first_comp_this_txt, num_comps_this_txt;
    Section sl;
    uint64_t num_nondrop_lines = 0;

    writer_get_txt_file_info (&first_comp_this_txt, &num_comps_this_txt, &sl);

    while (!dispatcher_is_done (dispatcher)) {

        bool achieved_something = false;
        
        // In input is not exhausted, and a compute thread is available - read a vblock and dispatch it
        if (!dispatcher_is_input_exhausted (dispatcher) && dispatcher_has_free_thread (dispatcher) && vb_has_free_vb()) {
            achieved_something = true;

            bool found_header = sections_next_sec2 (&sl, SEC_TXT_HEADER, SEC_VB_HEADER);

            // case SEC_TXT_HEADER
            if (found_header && sl->st == SEC_TXT_HEADER && z_file->num_txt_components_so_far < first_comp_this_txt + num_comps_this_txt) { 

                // case: skip entire component 
                if (writer_is_component_no_read (z_file->num_txt_components_so_far)) 
                    sl = sections_component_last (sl);

                else {
                    // note: also starts writer, and if unbinding, also opens the txt file and hands data over to the writer
                    Coords rejects_coord = txtheader_piz_read_and_reconstruct (z_file->num_txt_components_so_far, sl); 

                    // case --unbind: unpausing after previous txt_file pause (requires txt file to be open)
                    if (flag.unbind) dispatcher_resume (dispatcher);  

                    if (flag.header_only && !rejects_coord)
                        sl = sections_component_last (sl); // skip all VBs
                }
                z_file->num_txt_components_so_far++;
            }

            // case SEC_VB_HEADER
            else if (found_header && sl->st == SEC_VB_HEADER) 
                header_only_file &= piz_dispatch_one_vb (dispatcher, sl);  // function returns true if VB was skipped

            // case: we're done with this txt_file (either no header bc EOF, or TXT_HEADER belongs to the next txt_file when unbinding)
            else {
                dispatcher_set_input_exhausted (dispatcher, true);

                if (header_only_file)
                    dispatcher_recycle_vbs (dispatcher, true);
            }
        }

        // if the next thread (by sequential order) is ready, handle the reconstructed VB
        VBlock *recon_vb = dispatcher_get_processed_vb (dispatcher, NULL, false);  // non-blocking
        if (recon_vb) piz_handle_reconstructed_vb (dispatcher, recon_vb, &num_nondrop_lines);

        if (!achieved_something) usleep (30000); // nothing for us to do right now - wait 30ms
    }

    // verifies reconstructed file against MD5 (if compressed with --md5 or --test) or Adler2 and/or codec_args (if bgzf)
    Digest decompressed_file_digest = piz_one_verify_digest();

    if (!flag.test) progress_finalize_component_time ("Done", decompressed_file_digest);

    // finish writing the txt_file
    writer_finish_writing (is_last_txt_file);

    // --show-sex and --show-coverage - output results
    if (txt_file && !flag_loading_auxiliary) {
        if (flag.show_coverage) coverage_show_coverage();
        if (flag.show_sex) coverage_sex_classifier (is_first_z_file);
        if (flag.idxstats) coverage_show_idxstats();
        if (flag.count == CNT_TOTAL) iprintf ("%"PRIu64"\n", num_nondrop_lines);
    }

    if (is_last_txt_file) dispatcher_finish (&dispatcher, NULL);
    else                  dispatcher_pause (dispatcher); // we're unbinding and still have more txt_files
     
    DT_FUNC (z_file, piz_finalize)();

    return true;
}
