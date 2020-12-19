// ------------------------------------------------------------------
//   piz.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
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
#include "endianness.h"
#include "piz.h"
#include "sections.h"
#include "random_access.h"
#include "regions.h"
#include "dict_id.h"
#include "reference.h"
#include "refhash.h"
#include "progress.h"
#include "profiler.h"
#include "stats.h"
#include "codec.h"
#include "bgzf.h"
#include "flags.h"
#include "reconstruct.h"

// PIZ compute thread: decompress all contexts
// ZIP compute thread in FASTQ: decompress pair_1 contexts when compressing pair_2
uint32_t piz_uncompress_all_ctxs (VBlock *vb, 
                                  uint32_t pair_vb_i) // used in ZIP when uncompressing previous file's paired sections
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);
    
    uint32_t section_i = pair_vb_i ? 0 : 1; // normally, we skip the VB header, but when uncompressing paired sections there is no VB header
    while (section_i < vb->z_section_headers.len) {

        if (section_i == vb->z_section_headers.len) break; // no more sections left

        SectionHeaderCtx *header = (SectionHeaderCtx *)ENT (char, vb->z_data, section_index[section_i]);

        bool is_local = header->h.section_type == SEC_LOCAL;
        if (header->h.section_type == SEC_B250 || is_local) {
            Context *ctx = ctx_get_ctx (vb, header->dict_id); // creates the context

            // case: in ZIP & we are PAIR_2, reading PAIR_1 info - save flags, ltype, lcodec to restore later
            struct FlagsCtx save_flags={}; LocalType save_ltype=0; Codec save_lcodec=0;
            if (pair_vb_i && command == ZIP) { save_flags=ctx->flags; save_ltype=ctx->ltype; save_lcodec=ctx->lcodec; }
            
            *(uint8_t*)&ctx->flags |= header->h.flags.flags; // ??? is |= really needed? or is assignment sufficient?

            ctx->ltype  = header->ltype;
            if (is_local) ctx->lcodec = header->h.codec;

            // case: in PIZ: flags.paired appears on the sections the "pair 2" VB (that come first in section_index)
            if (ctx->flags.paired && !pair_vb_i) 
                pair_vb_i = fastq_get_pair_vb_i (vb);

            bool is_pair_section = (BGEN32 (header->h.vblock_i) == pair_vb_i); // is this a section of "pair 1" 
            
            if (is_pair_section) ctx->pair_flags = header->h.flags.ctx;            

            Buffer *target_buf = is_local ? &ctx->local : &ctx->b250;

            zfile_uncompress_section (vb, header, 
                                      is_pair_section ? &ctx->pair      : target_buf, 
                                      is_pair_section ? "contexts.pair" : is_local ? "contexts.local" : "contexts.b250", 
                                      is_pair_section ? pair_vb_i : vb->vblock_i,
                                      header->h.section_type); 

            if (!is_pair_section && is_local && dict_id_printable (ctx->dict_id).num == flag.dump_one_local_dict_id.num) 
                ctx_dump_binary (vb, ctx, true);

            if (!is_pair_section && !is_local && dict_id_printable (ctx->dict_id).num == flag.dump_one_b250_dict_id.num) 
                ctx_dump_binary (vb, ctx, false);

#           define adjust_lens(buf) { \
                buf.len /= lt_desc[ctx->ltype].width; \
                if (ctx->ltype == LT_BITMAP) { \
                    buf.param = buf.len * 64 - header->param ; /* number of bits */ \
                    LTEN_bit_array (buf_get_bitarray (&buf)); \
                } \
                else if (ctx->ltype >= LT_INT8 && ctx->ltype <= LT_UINT64)    \
                    lt_desc[ctx->ltype].file_to_native (&buf); \
            }

            if      (is_pair_section) adjust_lens (ctx->pair)
            else if (is_local)        adjust_lens (ctx->local)

            if (header->h.flags.ctx.copy_param)
                target_buf->param = header->param;

            // restore
            if (pair_vb_i && command == ZIP) { ctx->flags=save_flags; ctx->ltype=save_ltype; ctx->lcodec=save_lcodec; }

            section_i++;
        }    

        else break;
    }

    // initialize pair iterators (pairs only exist in fastq)
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) {
        Context *ctx = &vb->contexts[did_i];
        if (buf_is_allocated (&ctx->pair))
            ctx->pair_b250_iter = (SnipIterator){ .next_b250 = FIRSTENT (uint8_t, ctx->pair),
                                                  .prev_word_index = -1 };
    }

    ctx_map_aliases (vb);

    return section_i;
}

static void piz_uncompress_one_vb (VBlock *vb)
{
    START_TIMER;

    ASSERT0 (!flag.reference || (genome && genome->nbits), "Error in piz_uncompress_one_vb: reference is not loaded correctly");

    DtTranslation trans = dt_get_translation(); // in case we're translating from one data type to another

    // note: txt_data is fully allocated in advance and cannot be extended mid-reconstruction (container_reconstruct_do and possibly others rely on this)
    buf_alloc (vb, &vb->txt_data, vb->vb_data_size * trans.factor + 10000, 1.1, "txt_data"); // +10000 as sometimes we pre-read control data (eg container templates) and then roll back

    piz_uncompress_all_ctxs (vb, 0);

    // genocat flags that with which we needn't reconstruct
    if (exe_type == EXE_GENOCAT && flag.genocat_info_only) goto done;

    // reconstruct from top level snip
    reconstruct_from_ctx (vb, trans.toplevel, 0, true);

    // compress txt_data into BGZF blocks (in vb->compressed) if applicable
    if (flag.bgzf) bgzf_compress_vb (vb);

    // calculate the digest contribution of this VB to the single file and bound files, and the digest snapshot of this VB
    if (!v8_digest_is_zero (vb->digest_so_far) && !flag.data_modified) 
        digest_one_vb (vb); 

done:
    vb->is_processed = true; /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ 
    COPY_TIMER (compute);
}

static void piz_read_all_ctxs (VBlock *vb, ConstSectionListEntryP *next_sl)
{
    // ctxs that have dictionaries are already initialized, but others (eg local data only) are not
    ctx_initialize_primary_field_ctxs (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_contexts);

    while ((*next_sl)->section_type == SEC_B250 || (*next_sl)->section_type == SEC_LOCAL) {
        uint32_t section_start = vb->z_data.len;
        *ENT (uint32_t, vb->z_section_headers, vb->z_section_headers.len) = section_start; 

        int32_t ret = zfile_read_section (z_file, vb, vb->vblock_i, &vb->z_data, "z_data", (*next_sl)->section_type, *next_sl); // returns 0 if section is skipped

        if (ret) vb->z_section_headers.len++;
        (*next_sl)++;                             
    }
}

// Called by PIZ I/O thread: read all the sections at the end of the file, before starting to process VBs
static DataType piz_read_global_area (Digest *original_file_digest) // out
{
    bool success = zfile_read_genozip_header (original_file_digest, 0, 0, 0);
    
    if (flag.show_stats) stats_read_and_display();

    if (!success) return DT_NONE;

    dict_id_initialize (z_file->data_type); // must run after zfile_read_genozip_header that sets z_file->data_type
    
    // check if the genozip file includes a reference
    bool has_ref_sections = !!sections_get_first_section_of_type (SEC_REFERENCE, true);

    ASSERT (!has_ref_sections || flag.reference != REF_EXTERNAL || flag.reading_reference, 
            "Error: cannot use --reference with %s because it was not compressed with --reference", z_name);

    if (!flag.reading_reference && has_ref_sections) 
        flag.reference = REF_STORED; // possibly override REF_EXTERNAL (it will be restored for the next file in )

    // if the user wants to see only the header, we can skip the dictionaries, regions and random access
    if (!flag.header_only) {
        
        ctx_read_all_dictionaries (DICTREAD_ALL); // read all CHROM/RNAME dictionaries - needed for regions_make_chregs()

        // update chrom node indices using the CHROM dictionary, for the user-specified regions (in case -r/-R were specified)
        regions_make_chregs();

        // if the regions are negative, transform them to the positive complement instead
        regions_transform_negative_to_positive_complement();

        // if this is a stored reference we load the reference random access that will determined which reference sections
        // should be read & uncompressed in case of --regions.
        // note: in case of a data file with stored reference - SEC_REF_RAND_ACC will contain the random access of the reference
        // and SEC_RANDOM_ACCESS will contain the random access of the data. In case of a .ref.genozip file, both sections exist 
        // and are identical. It made the coding easier and their size is negligible.
        random_access_load_ra_section (SEC_RANDOM_ACCESS, &z_file->ra_buf, "z_file->ra_buf", 
                                       flag.show_index ? "Random-access index contents (result of --show-index)" : NULL);

        random_access_load_ra_section (SEC_REF_RAND_ACC, &ref_stored_ra, "ref_stored_ra", 
                                       flag.show_ref_index && !flag.reading_reference ? "Reference random-access index contents (result of --show-index)" : NULL);

        if ((flag.reference == REF_STORED || flag.reference == REF_EXTERNAL) && 
            !flag.reading_reference &&
            !(flag.show_headers && exe_type == EXE_GENOCAT))
            ref_contigs_sort_chroms(); // create alphabetically sorted index for user file chrom word list

        ref_contigs_load_contigs(); // note: in case of REF_EXTERNAL, reference is already pre-loaded

        // mapping of the file's chroms to the reference chroms (for files originally compressed with REF_EXTERNAL/EXT_STORE and have alternative chroms)
        ref_alt_chroms_load();

        // case: reading reference file
        if (flag.reading_reference) {
            bool dispatcher_invoked = false;

            // attempt to mmap a cached reference, and if one doesn't exist, uncompress the reference file and cache it
            if (!ref_mmap_cached_reference()) {
                ref_load_stored_reference();
                ref_generate_reverse_complement_genome();

                // start creating the genome cache now in a background thread, but only if we loaded the entire reference
                if (!flag.regions) ref_create_cache_in_background(); 

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
            ref_load_stored_reference();

            // exit now if all we wanted was just to see the reference (we've already shown it)
            if ((flag.show_reference || flag.show_is_set || flag.show_ref_hash) && exe_type == EXE_GENOCAT) exit_ok;
        }

        // read dict_id aliases, if there are any
        dict_id_read_aliases();
    }
    
    file_seek (z_file, 0, SEEK_SET, false);

    return z_file->data_type;
}

static bool piz_read_one_vb (VBlock *vb)
{
    START_TIMER; 

    ConstSectionListEntryP sl = sections_vb_first (vb->vblock_i, false); 

    vb->vb_position_txt_file = txt_file->txt_data_so_far_single;
    
    int32_t vb_header_offset = zfile_read_section (z_file, vb, vb->vblock_i, &vb->z_data, "z_data", SEC_VB_HEADER, sl++); 

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)ENT (char, vb->z_data, vb_header_offset);
    vb->vb_data_size     = BGEN32 (header->vb_data_size); 
    vb->first_line       = BGEN32 (header->first_line);      
    vb->lines.len        = BGEN32 (header->num_lines);       
    vb->longest_line_len = BGEN32 (header->longest_line_len);
    vb->digest_so_far    = header->digest_so_far;

    // in case of --unbind, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every txt component
    if (flag.unbind) vb->vblock_i = BGEN32 (header->h.vblock_i);

    if (flag.show_vblocks) 
        fprintf (info_stream, "vb_i=%u first_line=%u num_lines=%u txt_size=%u genozip_size=%u longest_line_len=%u\n",
                 vb->vblock_i, vb->first_line, (uint32_t)vb->lines.len, vb->vb_data_size, BGEN32 (header->z_data_bytes), vb->longest_line_len);

    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);
    ctx_overlay_dictionaries_to_vb ((VBlockP)vb); /* overlay all dictionaries (not just those that have fragments in this vblock) to the vb */ 

    buf_alloc (vb, &vb->z_section_headers, (MAX_DICTS * 2 + 50) * sizeof(uint32_t), 0, "z_section_headers"); // room for section headers  

    NEXTENT (uint32_t, vb->z_section_headers) = vb_header_offset; // vb_header_offset is always 0 for VB header

    // read all b250 and local of all fields and subfields
    piz_read_all_ctxs (vb, &sl);

    // read additional sections and other logic specific to this data type
    bool ok_to_compute = DTPZ(piz_read_one_vb) ? DTPZ(piz_read_one_vb)(vb, sl) : true; // true if we should go forward with computing this VB (otherwise skip it)

    // calculate the BGZF blocks that the compute thread is expected to compress
    if (flag.bgzf) bgzf_calculate_blocks_one_vb (vb, vb->vb_data_size);

    COPY_TIMER (piz_read_one_vb); 

    return ok_to_compute;
}

static Digest piz_one_file_verify_digest (Digest original_file_digest)
{
    if (v8_digest_is_zero (original_file_digest) || flag.genocat_info_only || flag.data_modified) return DIGEST_NONE; // we can't calculate the digest for some reason

    Digest decompressed_file_digest = digest_finalize (&txt_file->digest_ctx_bound, "file:digest_ctx_bound"); // z_file might be a bound file - this is the MD5 of the entire bound file
    char s[200]; 

    if (digest_is_zero (original_file_digest)) { 
        sprintf (s, "%s = %s", digest_name(), digest_display (decompressed_file_digest).s);
        progress_finalize_component (s); 
    }

    else if (digest_is_equal (decompressed_file_digest, original_file_digest)) {

        if (flag.test) { 
            sprintf (s, "%s = %s verified as identical to the original %s", 
                     digest_name(), digest_display (decompressed_file_digest).s, dt_name (txt_file->data_type));
            progress_finalize_component (s); 
        }
    }
    
    else if (flag.test) {
        progress_finalize_component ("FAILED!!!");
        ABORT ("Error: %s of original file=%s is different than decompressed file=%s\nPlease contact bugs@genozip.com to help fix this bug in genozip\n",
               digest_name(), digest_display (original_file_digest).s, digest_display (decompressed_file_digest).s);
    }

    // if compressed incorrectly - warn, but still give user access to the decompressed file
    else ASSERTW (digest_is_zero (original_file_digest), // its ok if we decompressed only a partial file
                  "File integrity error: %s of decompressed file %s is %s, but the original %s file's was %s", 
                  txt_file->name, digest_display (decompressed_file_digest).s, dt_name (txt_file->data_type), 
                  digest_name(), digest_display (original_file_digest).s);

    return decompressed_file_digest;
}

// returns false if VB was dispatched, and true if vb was skipped
bool piz_dispatch_one_vb (Dispatcher dispatcher, ConstSectionListEntryP sl_ent, uint32_t component_i)
{
    static Buffer region_ra_intersection_matrix = EMPTY_BUFFER; // we will move the data to the VB when we get it
    if (!random_access_is_vb_included (sl_ent->vblock_i, &region_ra_intersection_matrix)) return true; // skip this VB if not included

    VBlock *next_vb = dispatcher_generate_next_vb (dispatcher, sl_ent->vblock_i);
    next_vb->component_i = component_i;

    if (region_ra_intersection_matrix.data) {
        buf_copy (next_vb, &next_vb->region_ra_intersection_matrix, &region_ra_intersection_matrix, 0,0,0, "region_ra_intersection_matrix");
        buf_free (&region_ra_intersection_matrix); // note: copy & free rather than move - so memory blocks are preserved for VB re-use
    }
    
    // read one VB's genozip data
    bool grepped_out = !piz_read_one_vb (next_vb);

    if (grepped_out                                    || // this VB was filtered out by grep
        (flag.show_headers && exe_type == EXE_GENOCAT))   // we're not reconstructing VBs at all - only showing headers
        dispatcher_abandon_next_vb (dispatcher); 
    else
        dispatcher_compute (dispatcher, piz_uncompress_one_vb);

    return false;
}

// called once per components reconstructed txt_file: i.e.. if unbinding there will be multiple calls.
void piz_one_file (uint32_t component_i /* 0 if not unbinding */, bool is_last_z_file)
{
    static Dispatcher dispatcher = NULL; // static dispatcher - with flag.unbind, we use the same dispatcher when pizzing components
    static ConstSectionListEntryP sl_ent = NULL, sl_ent_leaf_2=NULL; // preserve for unbinding multiple files

    // read genozip header
    Digest original_file_digest = DIGEST_NONE;

    // read genozip header, dictionaries etc and set the data type when reading the first component of in case of --unbind, 
    static DataType data_type = DT_NONE; 
    bool no_more_headers = false;

    if (component_i == 0) {

        data_type = piz_read_global_area (&original_file_digest);
        if (data_type == DT_NONE || flag.reading_reference) {
            no_more_headers = true; // reference file has no VBs
            goto finish; 
        }

        sl_ent = sl_ent_leaf_2 = NULL; // reset

        ASSERT (!flag.test || !digest_is_zero (original_file_digest), 
                "Error testing %s: --test cannot be used with this file, as it was not compressed (in genozip v8) with --md5 or --test", z_name);

        if (flag.test || flag.md5) 
            ASSERT0 (dt_get_translation().is_src_dt, "Error: --test or --md5 cannot be used when converting a file to another format"); 

        dispatcher = dispatcher_init ("piz", flag.xthreads ? 1 : global_max_threads, 
                                      0, flag.test, is_last_z_file, true, z_file->basename, PROGRESS_PERCENT, 0);
    }
    
    else if (flag.unbind) {
        if (!sections_get_next_section_of_type (&sl_ent, SEC_TXT_HEADER, false, true)) return; // unbinding - no more components
        sl_ent--; // rewind;
    
        dispatcher_set_input_exhausted (dispatcher, false); // accept more input 
    }
  
    if (DTPZ(piz_initialize)) DTPZ(piz_initialize)();

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. In input is not exhausted, and a compute thread is available - read a variant block and compute it
    // 2. Wait for the first thread (by sequential order) to complete and write data

    bool do_interleave = (flag.interleave && !flag.reading_reference);
    bool header_only_file = true; // initialize - true until we encounter a VB header
    bool first_component_this_txtfile = true;
    bool is_leaf_2 = false; // used when interleaving - true if we are reading a VB from the 2nd leaf

    while (!dispatcher_is_done (dispatcher)) {

        // PRIORITY 1: In input is not exhausted, and a compute thread is available - read a vblock and compute it
        if (!dispatcher_is_input_exhausted (dispatcher) && dispatcher_has_free_thread (dispatcher)) {

            ConstSectionListEntryP *sl_p = is_leaf_2 ? &sl_ent_leaf_2 : &sl_ent;

            // note when interleaving: either leaf 1 or 2 may encounter a VB_HEADER, but only leaf_1 will encounter a TXT_HEADER bc we skipped it for leaf 2
            bool another_header = sections_get_next_section_of_type2 (sl_p, SEC_TXT_HEADER, SEC_VB_HEADER, false, true);

            // if we're interleaving and leaf 1 encountered a TXT_HEADER - this is the header of leaf_2. We skip to the next one
            if (  another_header && sl_ent->section_type == SEC_TXT_HEADER && // this is leaf 1 and we encoutered a TXT_HEADER
                  do_interleave && sl_ent_leaf_2) { // leaf 2 was already read (we set sl_ent_leaf_2 when we skipped its TXT_HEADER)
                sl_ent_leaf_2 = NULL;
                another_header = sections_get_next_section_of_type (&sl_ent, SEC_TXT_HEADER, false, false); // move to the next TXT_HEADER beyond leaf 2
            }

            // case: SEC_VB_HEADER
            if (another_header && sl_ent->section_type == SEC_VB_HEADER) {
                
                if (flag.one_vb && flag.one_vb != (*sl_p)->vblock_i) // we want only one VB, but not this one
                    { (*sl_p)++; continue; }

                header_only_file &= piz_dispatch_one_vb (dispatcher, (*sl_p), component_i + is_leaf_2);  // function returns true if VB was skipped

                // if we're interleaving - next VB will be from the opposite leaf
                if (do_interleave) is_leaf_2 = !is_leaf_2;
            } 

            // case SEC_TXT_HEADER when concatenating or first TXT_HEADER when unbinding: proceed with this component 
            // note: this never happens in the 2nd leaf when interleaving, because we skipped the header
            else if (another_header && (!flag.unbind || first_component_this_txtfile)) {

                txtfile_genozip_to_txt_header (sl_ent,  
                                               first_component_this_txtfile ? &original_file_digest : NULL); // NULL means skip txt header (2nd+ component if concatenating)
                
                if (!first_component_this_txtfile) 
                    component_i += flag.interleave ? 2 : 1; 

                first_component_this_txtfile = false;
                
                dispatcher_resume (dispatcher);  // in case it was paused by previous component when unbinding

                if (flag.header_only) goto finish;

                // if interleaving, set the start of leaf_2 to after its TXT_HEADER
                if (do_interleave) { 
                    sl_ent_leaf_2 = sl_ent;
                    ASSERT (sections_get_next_section_of_type (&sl_ent_leaf_2, SEC_TXT_HEADER, false, false), // next sections_get_next_section_of_type2 will read after the TXT_HEADER
                            "Error in piz_one_file: cannot find SEC_TXT_HEADER of 2nd leaf when trying to interleave in %s", z_name);
                }
            }

            // case: we're done (concatenating: no more VBs in the entire file ; unbinding: no more VBs in our component)
            else { 
                no_more_headers = !another_header;

                sl_ent--; // re-read in next call to this function if unbinding
                dispatcher_set_input_exhausted (dispatcher, true);

                if (header_only_file)
                    dispatcher_recycle_vbs (dispatcher);
            }
        }

        // PRIORITY 2 (interleaving): Wait for two VBs, one from each leaf (by sequential order) and interleave them
        else if (do_interleave) { 
            VBlock *processed_vb_1 = dispatcher_get_processed_vb (dispatcher, NULL); 
            VBlock *processed_vb_2 = dispatcher_get_processed_vb (dispatcher, NULL); 

            // write the two VBs - interleaving their lines
            fastq_txtfile_write_one_vblock_interleave (processed_vb_1, processed_vb_2);

            z_file->num_vbs += 2;
            z_file->txt_data_so_far_single += processed_vb_1->vb_data_size + processed_vb_2->vb_data_size; 

            dispatcher_recycle_vbs (dispatcher);
        }

        // PRIORITY 2 (not interleaving): Wait for the first thread (by sequential order) to complete and write data
        else {
            VBlock *processed_vb = dispatcher_get_processed_vb (dispatcher, NULL); 

            // read of a normal file - output uncompressed block (unless we're reading a reference - we don't need to output it)
            if (!flag.reading_reference) txtfile_write_one_vblock (processed_vb);

            z_file->num_vbs++;
            z_file->txt_data_so_far_single += processed_vb->vb_data_size; 

            dispatcher_recycle_vbs (dispatcher);
        }
    }

    // verifies reconstructed file against MD5 (if compressed with --md5 or --test) or Adler2 and/or codec_args (if bgzf)
    Digest decompressed_file_digest = piz_one_file_verify_digest (original_file_digest);

    if (!flag.test) progress_finalize_component_time ("Done", decompressed_file_digest);

finish:
    // case: we're unbinding and still have more components - we continue with the same dispatcher in the next component.
    if (!no_more_headers) 
        dispatcher_pause (dispatcher);

    // case: we're done with reconstructing this z_file - either concatenated, or this was the last component unbound
    else dispatcher_finish (&dispatcher, NULL);    

    DT_FUNC (z_file, piz_finalize)();
}
