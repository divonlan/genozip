// ------------------------------------------------------------------
//   zip.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "dispatcher.h"
#include "zip.h"
#include "seg.h"
#include "random_access.h"
#include "reference.h"
#include "refhash.h"
#include "ref_iupacs.h"
#include "progress.h"
#include "fastq.h"
#include "stats.h"
#include "compressor.h"
#include "bgzf.h"
#include "txtheader.h"
#include "threads.h"
#include "contigs.h"
#include "chrom.h"
#include "biopsy.h"
#include "dict_io.h"
#include "gencomp.h"
#include "aliases.h"
#include "arch.h"
#include "user_message.h"
#include "zriter.h"
#include "b250.h"
#include "zip_dyn_int.h"

static void zip_display_compression_ratio (Digest md5)
{
    float z_bytes        = MAX_((float)z_file->disk_so_far, 1.0); // at least one, to avoid division by zero in case of a z_bytes=0 issue
    float plain_bytes    = (float)z_file->txt_data_so_far_bind;
    float comp_bytes     = file_is_read_via_ext_decompressor (txt_file) 
                              ? (float)txt_file->disk_size    // 0 if via pipe or url, as we have no knowledge of file size
                              : (float)txt_file->disk_so_far; // unlike disk_size, works also for piped-in files (but not CRAM, BCF, XZ, ZIP)
    float ratio_vs_plain = plain_bytes / z_bytes;
    float ratio_vs_comp  = -1;

    if (flag.debug_progress) 
        iprintf ("Ratio calculation: ratio_vs_plain=%f = plain_bytes=%"PRIu64" / z_bytes=%"PRIu64"\n",
                    ratio_vs_plain, (uint64_t)plain_bytes, (uint64_t)z_bytes);

    // in bind mode, we don't show compression ratio for files except for the last one
    if (flag.bind) { 

        static float comp_bytes_bind = 0;
        static FileType source_file_type = UNKNOWN_FILE_TYPE;

        // reset for every set of bound files (we might have multiple sets of --pair)
        if (z_file->num_txts_so_far == 1) {
            comp_bytes_bind=0; 
            source_file_type = txt_file->type;
        }

        else if (source_file_type != txt_file->type) // heterogenous source file types
            source_file_type = UNKNOWN_FILE_TYPE;

        comp_bytes_bind += comp_bytes;

        if (z_file->z_closes_after_me) { 
            ratio_vs_comp = comp_bytes_bind / z_bytes; // compression vs .gz/.bz2/.bcf/.xz... size
            if (flag.debug_progress) 
                iprintf ("Ratio calculation: ratio_vs_comp=%f = comp_bytes_bind=%"PRIu64" / z_bytes=%"PRIu64"\n",
                         ratio_vs_comp, (uint64_t)comp_bytes_bind, (uint64_t)z_bytes);
        }
        else 
            progress_finalize_component_time ("Done", md5);
    }
    else {
        ratio_vs_comp = comp_bytes / z_bytes; // compression vs .gz/.bz2/.bcf/.xz... size
        if (flag.debug_progress) 
            iprintf ("Ratio calculation: ratio_vs_comp=%f = comp_bytes=%"PRIu64" / z_bytes=%"PRIu64"\n",
                     ratio_vs_comp, (uint64_t)comp_bytes, (uint64_t)z_bytes);
    }

    // in bound files, for the non-last components, we already printed "Done" above
    if (flag.bind && !z_file->z_closes_after_me) {}

    // when making a reference, we don't care about the compression
    else if (flag.make_reference || flag.zip_no_z_file)
        progress_finalize_component_time ("Done", md5);

    // Deep - to complicated to communicate compression vs FASTQ and BAM/SAM source files - show only vs  compression
    else if (flag.deep) 
        progress_finalize_component_time_ratio ("Deep", ratio_vs_comp, md5);

    // when compressing BAM report only ratio_vs_comp (compare to BGZF-compress BAM - we don't care about the underlying plain BAM)
    // Likewise, doesn't have a compression extension (eg .gz), even though it may actually be compressed eg .tbi (which is actually BGZF)
    else if (Z_DT(BAM) || (txt_file && file_get_codec_by_txt_ft (txt_file->data_type, txt_file->type, false) == CODEC_NONE)) 
        progress_finalize_component_time_ratio (txt_file->source_codec == CODEC_CRAM ? "CRAM" : z_dt_name(), ratio_vs_comp, md5);

    else if (ratio_vs_comp >= 0) {
        if (txt_file->codec == CODEC_NONE || ratio_vs_comp < 1.05) // disk_so_far doesn't give us the true txt file size 
            progress_finalize_component_time_ratio (z_dt_name(), ratio_vs_plain, md5);
        
        else // source was compressed
            progress_finalize_component_time_ratio_better (z_dt_name(), ratio_vs_plain, file_exts[txt_file->type], ratio_vs_comp, md5);
    }
}

static struct { DataType dt; const uint64_t dict_id_num; LocalGetLineCB *func; } callbacks[] = LOCAL_GET_LINE_CALLBACKS;

LocalGetLineCB *zip_get_local_data_callback (DataType dt, ContextP ctx)
{
    if (ctx && !ctx->no_callback)
        for (unsigned i=0; i < ARRAY_LEN(callbacks); i++)
            if (callbacks[i].dt == dt && callbacks[i].dict_id_num == ctx->dict_id.num) 
                return callbacks[i].func;

    return NULL;
}

void zip_set_no_stons_if_callback (VBlockP vb)
{
    for (unsigned i=0; i < ARRAY_LEN(callbacks); i++)
        if (callbacks[i].dt == vb->data_type) {
            ContextP ctx = ctx_get_ctx (vb, callbacks[i].dict_id_num);
            
            if (!ctx->no_callback) ctx->no_stons = true;
        }
}

// after segging - if any context appears to contain only singleton snips (eg a unique ID),
// we move it to local instead of needlessly cluttering the global dictionary
static void zip_handle_unique_words_ctxs (VBlockP vb)
{
    START_TIMER;

    for_ctx {
        if (ctx->local.len ||  // local is not free to accept our singletons
            ctx->no_stons  ||  // don't change to LT_SINGLETON if we were explicitly forbidden having singletons
            ctx->ltype == LT_SUPP || // local data might be created by codec (later)
            ((VB_DT(VCF) || VB_DT(BCF)) && dict_id_is_vcf_format_sf (ctx->dict_id))) // this doesn't work for FORMAT fields
            continue;
            
        // reset ltype to LT_SINGLETON so that we can use local for singletons (either here or in ctx_commit_node - subject to conditions).
        // note: ltype was possibly assigned a different value in *_seg_initialize, but then local not utilized.
        ctx->ltype = LT_SINGLETON; 

        if (!ctx->nodes.len                    || // no new words in this VB 
            ctx->nodes.len != ctx->b250.count  || // not all new words in this VB are singletons
            ctx->nodes.len < vb->lines.len / 5 || // don't bother if this is a rare field less than 20% of the lines
            !ctx_can_have_singletons (ctx)     || // this context is not allowed to have singletons
            ctx->b250.count == 1)                 // only one word - better to handle with all_the_same rather than singleton
            continue;
        
        buf_free (ctx->local); // possibly local was allocated, but then not utilized
        buf_move (vb, ctx->local, CTX_TAG_LOCAL, ctx->dict);
        buf_free (ctx->nodes);
        buf_free (ctx->b250);
    }

    COPY_TIMER (zip_handle_unique_words_ctxs);
}

static bool zip_generate_local (VBlockP vb, ContextP ctx)
{
    START_TIMER;
    
    ASSERT (ctx->dict_id.num, "tag_name=%s did_i=%u: ctx->dict_id=0 despite ctx->local containing data", ctx->tag_name, (unsigned)(ctx - vb->contexts));

    ctx->ltype = dyn_int_get_ltype (ctx);

    // case: local is LTEN (instead of native endianity) and machine is BGEN, so BGEN_*_buf ^ above did nothing.     
    bool need_lten = (ctx->local_is_lten && !flag.is_lten);
    
    switch (ctx->ltype) {
        case LT_BITMAP    : 
            LTEN_bits ((BitsP)&ctx->local);    
            ctx->local.prm8[0] = ((uint8_t)64 - (uint8_t)(ctx->local.nbits % 64)) % (uint8_t)64;
            ctx->local_param   = true;
            break;

        case LT_UINT32 : case LT_hex32 : case LT_HEX32 : case LT_FLOAT32 : 
            if (need_lten) LTEN_u32_buf (&ctx->local, NULL);   
            else           BGEN_u32_buf (&ctx->local, NULL);    
            break;

        case LT_UINT16 : case LT_hex16 : case LT_HEX16 : 
            if (need_lten) LTEN_u16_buf (&ctx->local, NULL);   
            else           BGEN_u16_buf (&ctx->local, NULL);    
            break;

        case LT_UINT64 : case LT_hex64 : case LT_HEX64 : case LT_FLOAT64 :
            if (need_lten) LTEN_u64_buf (&ctx->local, NULL);   
            else           BGEN_u64_buf (&ctx->local, NULL);    
            break;

        case LT_INT8      : interlace_d8_buf  (&ctx->local, NULL);              
                            break;

        case LT_INT16     : if (need_lten) LTEN_interlace_d16_buf (&ctx->local, NULL);   
                            else           BGEN_interlace_d16_buf (&ctx->local, NULL);    
                            break;

        case LT_INT32     : if (need_lten) LTEN_interlace_d32_buf (&ctx->local, NULL);   
                            else           BGEN_interlace_d32_buf (&ctx->local, NULL);    
                            break;

        case LT_INT64     : if (need_lten) LTEN_interlace_d64_buf (&ctx->local, NULL);   
                            else           BGEN_interlace_d64_buf (&ctx->local, NULL);    
                            break;

        default           : break;        
    }

    // transpose if needed AND local is rectangular
    if (ctx->dyn_transposed)
        dyn_int_transpose (vb, ctx);            

    COPY_TIMER (zip_generate_local); // codec_assign measures its own time

    // in case we are using "pair identical", drop this section if it is an R2 section identical to its R1 counterpart
    if (is_fastq_pair_2 (vb) && fastq_zip_use_pair_identical (ctx->dict_id) && 
        buf_issame (&ctx->local, &ctx->localR1, lt_width(ctx))) {
        
        if (flag.debug_generate) 
            iprintf ("%s: %s[%u].local dropped because it is an R2 section which is identical to its R1 counterpart\n", 
                     VB_NAME, ctx->tag_name, ctx->did_i);

        // note: careful not to set local.len=0 bc this is called before merge, and ctx_drop_all_the_same (that 
        // is part of merge), relies on local.len to decide if it can drop b250 of pair-identical R2.b250 sections    
        return false;
    }

    codec_assign_best_codec (vb, ctx, NULL, SEC_LOCAL);

    if (flag.debug_generate) 
        iprintf ("%s: %s[%u].local ltype=%s len=%"PRIu64" codec=%s\n", VB_NAME, ctx->tag_name, ctx->did_i,
                 lt_name (ctx->ltype), ctx->local.len, codec_name(ctx->lcodec));

    return true;
}

// generate & write b250 data for all contexts - do them in random order, to reduce the chance of multiple doing codec_assign_best_codec for the same context at the same time
// VBs doing codec_assign_best_codec at the same, so that they can benefit from pre-assiged codecs
void zip_compress_all_contexts_b250 (VBlockP vb)
{
    START_TIMER;
    threads_log_by_vb (vb, "zip", "START COMPRESSING B250", 0);
    
    // arrays of all contexts in this VB
    ContextP ctxs[vb->num_contexts];
    for (Did did_i=0; did_i < vb->num_contexts; did_i++) ctxs[did_i] = CTX(did_i);

        // in each iteration, pick a context at random and remove it from the list 
    for (unsigned i=0; i < vb->num_contexts; i++) {
 
        int ctx_i = global_max_threads > 1 ? ((clock()+1) * (vb->vblock_i+1)) % (vb->num_contexts - i) : 0; // force predictability with single thread 
        
        ContextP ctx = ctxs[ctx_i];
        memmove (&ctxs[ctx_i], &ctxs[ctx_i+1], (vb->num_contexts - i - ctx_i - 1) * sizeof (ContextP));

        if (!ctx->b250.len || ctx->b250_compressed) continue;

        if (!b250_zip_generate (vb, ctx))  // generate the final b250 buffers from their intermediate form
            continue; // dropped

        if (dict_id_typeless (ctx->dict_id).num == flag.dump_one_b250_dict_id.num) 
            ctx_dump_binary (vb, ctx, false);

        if (flag.show_time) codec_show_time (vb, "B250", ctx->tag_name, ctx->bcodec);
        
        if (HAS_DEBUG_SEG(ctx)) iprintf ("zip_compress_all_contexts_b250: vb=%s %s: B250.len=%"PRIu64" NODES.len=%"PRIu64"\n", 
                                         VB_NAME, ctx->tag_name, ctx->b250.len, ctx->nodes.len);

        START_TIMER; 
        zfile_compress_b250_data (vb, ctx);
        COPY_TIMER(fields[ctx->did_i]);    

        ctx->b250_compressed = true;
    }

    COPY_TIMER (zip_compress_ctxs); // same profiler for b250 and local as we breakdown by ctx underneath it
}

// generate & write local data for all contexts - in random order, to reduce the chance of multiple doing codec_assign_best_codec for the same context at the same time
static void zip_compress_all_contexts_local (VBlockP vb)
{
    START_TIMER;
    threads_log_by_vb (vb, "zip", "START COMPRESSING LOCAL", 0);

    // first we handle local_dep=0 then local_dep=1 and finally local_dep=2
    for (int dep_level=DEP_L0 ; dep_level < NUM_LOCAL_DEPENDENCY_LEVELS; dep_level++) {

        // initialize list of contexts at this dependency level that need compression
        ContextP ctxs[vb->num_contexts];
        unsigned num_ctxs=0;
        for_ctx_that ((ctx->local.len || ctx->local_always) && ctx->local_dep == dep_level && !ctx->local_compressed)
            ctxs[num_ctxs++] = ctx;

        while (num_ctxs) {
            // pick a context at "random" and remove it from the list (not random if single thread)
            int ctx_i = global_max_threads > 1 ? (65531 * (vb->vblock_i+1)) % num_ctxs : 0; 
            ContextP ctx = ctxs[ctx_i];
            memmove (&ctxs[ctx_i], &ctxs[ctx_i+1], (num_ctxs - (ctx_i+1)) * sizeof (ContextP));
            num_ctxs--;

            ctx->local_compressed = true; // so we don't compress it again

            if (!zip_generate_local (vb, ctx))
                continue; // section dropped

            if (dict_id_typeless (ctx->dict_id).num == flag.show_singletons_dict_id.num) 
                dict_io_show_singletons (vb, ctx);

            if (dict_id_typeless (ctx->dict_id).num == flag.dump_one_local_dict_id.num) 
                ctx_dump_binary (vb, ctx, true);

            if (flag.show_time) codec_show_time (vb, "LOCAL", ctx->tag_name, ctx->lcodec);

            if (HAS_DEBUG_SEG(ctx)) iprintf ("%s: zip_compress_all_contexts_local: %s: LOCAL.len=%"PRIu64" LOCAL.param=%"PRIu64"\n", 
                                            VB_NAME, ctx->tag_name, ctx->local.len, ctx->local.param);

            START_TIMER; 
            zfile_compress_local_data (vb, ctx, 0);
            COPY_TIMER(fields[ctx->did_i]);    

            if (!ctx->dict_merged) // note: if dict_merged, we are in the second call to this function, and local consists of singletons
                ctx->no_stons = true; // since we had data in local, we don't allow ctx_commit_node to move singletons to local

        }
    }

    COPY_TIMER (zip_compress_ctxs); // same profiler for b250 and local as we breakdown by ctx underneath it
}

void zip_init_vb (VBlockP vb)
{
    vb->recon_size = Ltxt; // initial value. it may change if --optimize / --chain are used, or if dual coordintes - for the other coordinate
    vb->txt_size   = Ltxt; // this copy doesn't change with --optimize / --chain.

    if (DTPT(zip_init_vb)) DTPT(zip_init_vb)(vb); // data-type specific initialization of the VB    
}

// called by main thread after VB has completed processing
static void zip_update_txt_counters (VBlockP vb)
{
    // note: in case of an FASTQ with flag.optimize_DESC or VCF with add_line_numbers, we already updated this in *_zip_init_vb
    if (!(flag.optimize_DESC && VB_DT(FASTQ)) &&
        !(flag.add_line_numbers && (VB_DT(VCF) || VB_DT(BCF))))
        txt_file->num_lines += vb->lines.len; // lines in this txt file

    // counters of data AS IT APPEARS IN THE TXT FILE
    z_file->num_lines += vb->lines.len; // lines in all bound files in this z_file

    z_file->comp_num_lines[vb->comp_i] += vb->lines.len; // set also for DVCF rejects

    z_file->txt_data_so_far_single_0 += (int64_t)vb->txt_size;  // length of data before any modifications
    z_file->txt_data_so_far_bind_0   += (int64_t)vb->txt_size;

    // counter of data FOR PROGRESS DISPLAY
    z_file->txt_data_so_far_single += (int64_t)vb->txt_size;   

    // counter of data in DEFAULT RECONSTRUCTION 
    z_file->txt_data_so_far_bind += vb->recon_size;

    // per-component data for stats
    z_file->txt_data_so_far_bind_0_comp[vb->comp_i] += (int64_t)vb->txt_size;
    
    // note: in case of SAM gencomp, MAIN, we add recon_size - assuming the discrepency vs txt_data.len
    // is only due to lines being deported to gencomp 
    z_file->txt_data_so_far_bind_comp[vb->comp_i] += 
        (z_sam_gencomp && vb->comp_i==SAM_COMP_MAIN) ? vb->recon_size : Ltxt;

    z_file->num_components = MAX_(z_file->num_components, vb->comp_i+1);

    // Note: no data-type-specific code here, instead, put in *_zip_after_compute
}

// ZIP: free z_file buffers no longer needed before global sections
static void zip_free_undeeded_zctx_bufs_after_seg (void)
{
    START_TIMER;

    for_zctx {
        buf_destroy (zctx->ston_hash);
        buf_destroy (zctx->ston_ents);
        buf_destroy (zctx->global_hash);
    }

    buf_destroy (z_file->sag_grps);
    buf_destroy (z_file->sag_grps_index);
    buf_destroy (z_file->sag_alns);
    buf_destroy (z_file->sag_qnames);
    buf_destroy (z_file->sag_depn_index);
    buf_destroy (z_file->sag_cigars);
    buf_destroy (z_file->sag_seq);
    buf_destroy (z_file->sag_qual);
    buf_destroy (z_file->vb_start_deep_line);
    buf_destroy (z_file->deep_ents);

    buf_low_level_release_memory_back_to_kernel();

    COPY_TIMER_EVB (zip_free_undeeded_zctx_bufs_after_seg);
}

// write all the sections at the end of the file, after all VB stuff has been written
static void zip_write_global_area (void)
{
    START_TIMER;

    #define THREAD_DEBUG(x) threads_log_by_vb (evb, "main_thread:global_area", #x, 0);

    if (!flag.show_memory) // in show-mem, keep these, so we can report them.
        zip_free_undeeded_zctx_bufs_after_seg();

    codec_qual_show_stats();

    // if we're making a reference, we need the RA data to populate the reference section chrome/first/last_pos ahead of ref_compress_ref
    THREAD_DEBUG (finalize_random_access);
    random_access_finalize_entries (&z_file->ra_buf); // sort RA, update entries that don't yet have a chrom_index

    THREAD_DEBUG (compress_dictionaries);
    dict_io_compress_dictionaries(); 

    THREAD_DEBUG (compress_counts);
    ctx_compress_counts();

    THREAD_DEBUG (compress_subdicts);
    ctx_compress_subdicts();

    // store a mapping of the file's chroms to the reference's contigs, if they are any different
    // note: not needed in REF_EXT_STORE, as we convert the stored ref_contigs to use chrom_index of the file's CHROM
    if (IS_REF_EXTERNAL && DTFZ(chrom) != DID_NONE) {
        THREAD_DEBUG (compress_chrom_2ref);
        chrom_2ref_compress(gref);
    }

    // output reference, if needed
    bool store_ref = (flag.reference & REF_STORED) || flag.make_reference;
    if (store_ref) {
        THREAD_DEBUG (compress_ref);
        ref_compress_ref();
    }

    if (flag.make_reference) {
        THREAD_DEBUG (compress_iupacs);
        ref_iupacs_compress();

        THREAD_DEBUG (compress_refhash);
        refhash_compress_refhash();
    }

    // compress alias list, if this data_type has any aliases defined
    THREAD_DEBUG (compress_aliases);
    aliases_compress();

    if (!segconf.disable_random_acccess) {
        THREAD_DEBUG (compress_random_access);
        // if this data has random access (i.e. it has chrom and pos), compress all random access records into evb->z_data
        Codec codec = random_access_compress (&z_file->ra_buf, SEC_RANDOM_ACCESS, CODEC_UNKNOWN, flag.show_index ? RA_MSG_PRIM : NULL);
    
        if (store_ref) 
            random_access_compress (ref_get_stored_ra (gref), SEC_REF_RAND_ACC, codec, flag.show_ref_index ? RA_MSG_REF : NULL);
    }

    THREAD_DEBUG (user_message);
    user_message_compress();

    THREAD_DEBUG (stats);
    stats_generate();

    // compress genozip header (including its payload sectionlist and footer) into evb->z_data
    zfile_compress_genozip_header();    

    stats_finalize();
    
    if (DTPZ(zip_free_end_of_z)) DTPZ(zip_free_end_of_z)();

    COPY_TIMER_EVB (zip_write_global_area);
}

// entry point for ZIP compute thread
static void zip_compress_one_vb (VBlockP vb)
{
    START_TIMER; 

    // we're just taking a biopsy of the txt data, so no need to actually compress. 
    if (flag.biopsy && 
        !(segconf.sag_type && vb->comp_i == SAM_COMP_MAIN)) // except in MAIN of SAM/BAM gencomp - need to generate PRIM and DEPN VBs 
        goto after_compress; 

    // if the txt file is compressed with BGZF, we uncompress now, in the compute thread
    if (txt_file->codec == CODEC_BGZF && flag.pair != PAIR_R2) 
        bgzf_uncompress_vb (vb);    // some of the blocks might already have been decompressed while reading - we decompress the remaining

    // calculate the digest contribution of this VB, and the digest snapshot of this VB
    if (zip_need_digest) 
        digest_one_vb (vb, true, NULL); // serializes VBs in order if MD5

    // allocate memory for the final compressed data of this vb. allocate 1/8 of the
    // vb size on the (uncompressed) txt file - this is normally plenty. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, 0, vb->txt_size / 8, char, CTX_GROWTH, "z_data");

    // clone global dictionaries while granted exclusive access to the global dictionaries
    ctx_clone (vb);

    // split each line in this VB to its components
    threads_log_by_vb (vb, "zip", "START SEG", 0);

    seg_all_data_lines (vb);

    if (flag.biopsy) goto after_compress; // in case of MAIN VB of SAM/BAM gencomp: we end our biopsy journey here

    // identify dictionaries that contain only singleton words (eg a unique id) and move the data from dict to local
    zip_handle_unique_words_ctxs (vb);

    zfile_compress_vb_header (vb); // vblock header

    if (flag.show_codec) {
        DO_ONCE iprintf ("\n\nThe output of --show-codec-test: Testing a sample of up %u bytes on ctx.local of each context.\n"
                         "Results in the format [codec bytes μsec] are in order of quality - the first was selected.\n", CODEC_ASSIGN_SAMPLE_SIZE);
    }

    bool need_compress = !flag.make_reference && !flag.seg_only;

    // while vb_i=1 is busy merging, other VBs can handle local
    if (vb->vblock_i != 1 && need_compress) 
        zip_compress_all_contexts_local (vb); // not yet locals that consist of singletons transferred from dict to local in ctx_merge_in_vb_ctx (these will have len=0 at this point)
    
    dispatcher_increment_progress ("compress1", PROGRESS_UNIT/2); // 1/2 compression done

    threads_log_by_vb (vb, "zip", "START MERGE", 0);

    // for --make-reference we serialize merging by VB, so that contigs get their word_index in the order of the reference file
    if (flag.make_reference) serializer_lock (make_ref_merge_serializer, vb->vblock_i);

    // merge new words added in this vb into the z_file.contexts (zctx), ahead of b250_zip_generate().
    // writing indices based on the merged dictionaries. all this is done while locking a mutex for each zctx.
    // note: vb>=2 will block here, until vb=1 is completed
    ctx_merge_in_vb_ctx (vb);

    if (flag.make_reference) serializer_unlock (make_ref_merge_serializer);

    if (need_compress) {
        zip_compress_all_contexts_local (vb); // for vb=1 - all locals ; for vb>1 - locals which consist of singletons set in ctx_merge_in_vb_ctx (other locals were already compressed above)
        zip_compress_all_contexts_b250 (vb);
    }

    dispatcher_increment_progress ("compress2", PROGRESS_UNIT-(PROGRESS_UNIT/2)); // 1/2 compression done

    // merge in random access - IF it is used
    if (!segconf.disable_random_acccess) 
        random_access_merge_in_vb (vb);
    
after_compress:
    // examples: compress data-type specific sections ; absorb gencomp lines
    DT_FUNC (vb, zip_after_compress)(vb);

    // tell dispatcher this thread is done and can be joined.
    vb_set_is_processed (vb); 

    COPY_TIMER (compute);
}

// data sent through dispatcher fan out functions - to do: make this an opaque struct
static VBIType prev_file_first_vb_i=0, prev_file_last_vb_i=0; // used if we're binding files - the vblock_i will continue from one file to the next

// main thread: returns true if successfully prepared a vb 
static void zip_prepare_one_vb_for_dispatching (VBlockP vb)
{
    // if we're compressing the 2nd file in a fastq pair (with --pair) - look back at the z_file data
    // and copy the data we need for this vb. note: we need to do this before txtfile_read_vblock as
    // we need the num_lines of the pair VB
    bool R1_data_exhausted = false;
    if (flag.pair == PAIR_R2) {
        uint32_t pair_vb_i = prev_file_first_vb_i + (vb->vblock_i-1 - prev_file_last_vb_i);
        
        if (pair_vb_i <= prev_file_last_vb_i)
            fastq_read_pair_1_data (vb, pair_vb_i); // add the R1 sections z_data after the R2 sections 
        else
            R1_data_exhausted = true; // R1 data is already exhausted. This is ok if R2 data is exhausted too.
    }

    // case: we have out-of-band txt_data waiting (for generated components) - compress this data first,
    // before reading more data from the txt_file
    if (gencomp_get_txt_data(vb)) 
        goto dispatch;

    else {        
        txtfile_read_vblock (vb);

        // --head diagnostic option in ZIP cuts a certain number of first lines from vb=1, and discards other VBs
        if (vb->vblock_i != 1 && flag.lines_last != NO_LINE) {
            vb->dispatch = DATA_EXHAUSTED;
            return;
        }

        if (Ltxt) 
            goto dispatch;

        else if (gencomp_am_i_expecting_more_txt_data()) // more data might be coming from MAIN VBs currently computing
            vb->dispatch = MORE_DATA_MIGHT_COME;

        else 
            vb->dispatch = DATA_EXHAUSTED;

        // note: the opposite case where R2 has less reads than R1 is caught in txtfile_read_vblock
        ASSINP (!R1_data_exhausted || vb->dispatch == DATA_EXHAUSTED, // we are expecting that if our pair R1 data is exhausted, then our R2 data is exhausted too
                "Error: File %s has more FASTQ reads than its R1 mate (vb=%s)", txt_name, VB_NAME);

        // error if stdin is empty - can happen only when redirecting eg "cat empty-file|./genozip -" (we test for empty regular files in main_genozip)
        ASSINP0 (vb->vblock_i > 1 || txt_file->txt_data_so_far_single /* txt header data */, 
                 "Error: Cannot compress stdin data because its size is 0");

        if (flag.biopsy && vb->dispatch == DATA_EXHAUSTED)
            biopsy_data_is_exhausted();

        if (flag.debug_or_test) buflist_test_overflows(vb, __FUNCTION__); 

        return;
    }

dispatch:
    vb->dispatch = READY_TO_COMPUTE;
    txt_file->num_vbs_dispatched++;

    if (vb->comp_i == COMP_MAIN) // note: we only update the MAIN comp from here, gen comps are updated
        gencomp_a_main_vb_has_been_dispatched();
}

// called main thread, as VBs complete (might be out-of-order)
static void zip_complete_processing_one_vb (VBlockP vb)
{
    DT_FUNC (vb, zip_after_compute)(vb);

    // update z_data in memory (its not written to disk yet)
    zfile_update_compressed_vb_header (vb); 
        
    txt_file->max_lines_per_vb = MAX_(txt_file->max_lines_per_vb, vb->lines.len);

    if (!flag.make_reference && !flag.seg_only)
        zfile_output_processed_vb_ext (vb, true);
    
    zip_update_txt_counters (vb);

    // destroy some buffers of "first generation" contexts (those that didn't clone any nodes)  
    if (vb->vblock_i < 100) // don't bother checking for high vb_i 
        for_ctx_that (ctx->nodes.len32 && !ctx->ol_nodes.len32) {
            buf_destroy (ctx->b250);       // 1st generation likely to have excessive length due to being all-new 4B nodes
            buf_destroy (ctx->local_hash); // 1st generation allocated based on wild guess
            buf_destroy (ctx->nodes);      // 1st generation likely to have a lot more new nodes (+dict) that subsequent generations
            buf_destroy (ctx->dict);       
        }

    dispatcher_increment_progress ("z_write", PROGRESS_UNIT); // writing done.

    z_file->num_vbs = MAX_(z_file->num_vbs, vb->vblock_i); // note: VBs are written out of order, so this can increase by 0, 1, or more than 1
    txt_file->num_vbs++;
}

// this is the main dispatcher function. It first processes the txt header, then proceeds to read 
// a VB from the input file and send it off to a thread for computation. When the thread
// completes, this function proceeds to write the output to the output file. It can dispatch
// several threads in parallel.
void zip_one_file (rom txt_basename, 
                   bool is_last_user_txt_file)  // the last user-specified txt file in this execution
{
    Dispatcher dispatcher = 0;
    dispatcher_start_wallclock();
    if (flag.show_time_comp_i == flag.zip_comp_i) profiler_initialize(); // re-start wallclock

    z_file->txt_data_so_far_single = 0;
    z_file->num_components         = MAX_(z_file->num_components, flag.zip_comp_i+1); // may increase further with generated components (in zip_update_txt_counters())
    evb->z_data.len                = 0;
    evb->z_next_header_i           = 0;
    
    // we calculate digest for each component seperately, stored in SectionHeaderTxtHeader (always 0 for generated components, or if modified)
    if (gencomp_comp_eligible_for_digest(NULL)) // if generated component - keep digest to display in progress after the last component
        z_file->digest_ctx = DIGEST_CONTEXT_NONE;

    if (!flag.bind || flag.zip_comp_i == COMP_MAIN) 
        prev_file_first_vb_i = prev_file_last_vb_i = 0; // reset if we're not binding

    // initalize pre-defined ctxs before segconf
    // note: generic_is_header_done as well as segconf may change the data type and re-initialize the contexts
    if (z_file->num_txts_so_far == 0)  // first component of this z_file 
        ctx_initialize_predefined_ctxs (z_file->contexts, txt_file->data_type, z_file->d2d_map, &z_file->num_contexts);

    segconf_initialize(); // before txtheader 

    uint32_t first_vb_i = prev_file_last_vb_i + 1;

    // read the txt header, assign the global variables, and write the compressed header to the GENOZIP file
    int64_t txt_header_offset = -1;
    int64_t txt_header_len = txtheader_zip_read_and_compress (&txt_header_offset, flag.zip_comp_i); // also increments z_file->num_txts_so_far

    bool success = (txt_header_len >= -1);
    if (!success) goto finish; // eg 2nd+ VCF file cannot bind, because of different sample names

    DT_FUNC (txt_file, zip_initialize)();

    segconf_calculate();

    DT_FUNC (txt_file, zip_after_segconf)();

    static uint64_t target_progress=0;
    if ((Z_DT(FASTQ) && flag.pair != PAIR_R2) ||   // note: if 2nd of a FASTQ file pair - we leave the target as it was in the first file as seggable_size is not calculated for the 2nd file
        (flag.deep && flag.zip_comp_i <= SAM_COMP_FQ00) ||
        (!flag.deep && !Z_DT(FASTQ))) {

        int64_t progress_unit = txt_file->est_num_lines ? txt_file->est_num_lines : txtfile_get_seggable_size(); 

        target_progress = progress_unit * 3 // read, seg, compress
                        + (!flag.make_reference && !flag.seg_only && !flag.biopsy) * progress_unit; // write
    }

    if (flag.debug_progress)
        iprintf ("zip_comp_i=%u : target_progress=%"PRIu64"\n", flag.zip_comp_i, target_progress);

    dispatcher = 
        dispatcher_fan_out_task (ZIP_TASK_NAME, txt_basename, 
                                 target_progress, // target progress: 1 for each read, compute, write
                                 target_progress ? NULL : txt_file->is_remote ? "Downloading & compressing..." : "Compressing...",
                                 !flag.make_reference,   // allow callbacks to zip_complete_processing_one_vb not in order of VBs (not allowed for make-reference as contigs need to be in consistent order)
                                 false,           // not test mode
                                 flag.xthreads, prev_file_last_vb_i, 5000, false,
                                 zip_prepare_one_vb_for_dispatching, 
                                 zip_compress_one_vb, 
                                 zip_complete_processing_one_vb);

    // verify that entire file was read (if file size is known)
    ASSERT (txt_file->disk_so_far == txt_file->disk_size || !txt_file->disk_size || file_is_read_via_ext_decompressor (txt_file),
            "Failed to compress entire file: file size is %"PRIu64", but only %"PRIu64" bytes were compressed",
            txt_file->disk_size, txt_file->disk_so_far);

    bgzf_finalize_discovery();

    zriter_wait_for_bg_writing(); // complete writing VBs before moving on

    dispatcher_calc_avg_compute_vbs (dispatcher);

    dispatcher_increment_progress ("txt_header", txt_file->est_num_lines ? 3 : (txt_header_len * 3)); //  account for txt_header read, computed and written

    // go back and update some fields in the txt header's section header and genozip header 
    if (txt_header_offset >= 0) // note: this will be -1 if we didn't write a SEC_TXT_HEADER section for any reason (e.g. SAM PRIM/DEPN, --make-reference...)
        zfile_update_txt_header_section_header (txt_header_offset);

    ASSERT0 (!flag.biopsy || biopsy_is_done(), "Biopsy request not complete - some VBs missing");

    // write the BGZF section containing BGZF block sizes, if this txt file is compressed with BGZF
    if (txt_file->codec == CODEC_BGZF)
        bgzf_compress_bgzf_section();

    // if this a non-bound file, or the last component of a bound file - write the genozip header, random access and dictionaries
finish:   
    z_file->txt_file_disk_sizes[flag.zip_comp_i] = txt_file->disk_size ? txt_file->disk_size // actual file size on disk, if we know it (we don't if its a remote or stdin file)
                                                                       : (int64_t)txt_file->disk_so_far + (txt_file->codec==CODEC_BGZF ? BGZF_EOF_LEN : 0); // data (plain, BGZF, GZ or BZ2) read from the file descriptor (we won't have correct src data here if reading through an external decompressor - but luckily txt_file->disk_size will capture that case)
    z_file->txt_file_disk_sizes_sum += z_file->txt_file_disk_sizes[flag.zip_comp_i];

    z_file->comp_codec[flag.zip_comp_i]        = txt_file->codec;
    z_file->comp_source_codec[flag.zip_comp_i] = txt_file->source_codec;

    // (re-)index sections after adding this txt_file 
    sections_create_index(); 

    // reconstruction plan (for VCF - for DVCF or --sort, for SAM - re-integrate supp/secondary alignments)
    if (!flag.seg_only && DTPZ(generate_recon_plan)) 
        DTPZ(generate_recon_plan)(); // should set z_file->z_closes_after_me if we need to close after this component after all

    if (z_file->z_closes_after_me && !flag.seg_only) { // note: for SAM, z_closes_after_me might be updated in sam_zip_generate_recon_plan
        // if we used the aligner with REF_EXT_STORE, we make sure all the CHROMs referenced are in the CHROM context, so
        // as SEC_REF_CONTIGS refers to them. We do this by seeing which contigs have any bit set in is_set.
        // note: in REF_EXTERNAL we don't use is_set, so we populate all contigs in zip_initialize
        // note: must be before zip_after_vbs() bc sam_zip_after_vbs() removes unused dict words (they are marked as used in ref_contigs_populate_aligned_chroms)
        if (flag.aligner_available && IS_REF_EXT_STORE) {
            THREAD_DEBUG (populate_aligned_chroms);
            ref_contigs_populate_aligned_chroms();
        }

        DT_FUNC (txt_file, zip_after_vbs)();
    
        zip_write_global_area();
    }

    zip_display_compression_ratio (digest_snapshot (&z_file->digest_ctx, NULL)); // Done for reference + final compression ratio calculation
    
    if (flag.md5 && flag.bind && z_file->z_closes_after_me &&
        ((flag.bind == BIND_FQ_PAIR && z_file->num_txts_so_far == 2) ||
         (flag.bind == BIND_SAM && z_file->num_txts_so_far == 3)))
        progress_concatenated_md5 (z_dt_name(), digest_snapshot (&z_file->digest_ctx, "file"));

    z_file->disk_size = z_file->disk_so_far;

    prev_file_first_vb_i = first_vb_i;
    dispatcher_finish (&dispatcher, &prev_file_last_vb_i, 
                       z_file->z_closes_after_me && !is_last_user_txt_file,
                       flag.show_memory && z_file->z_closes_after_me && is_last_user_txt_file); // show memory

    if (!z_file->z_closes_after_me) 
        ctx_reset_codec_commits(); 

    // no need to waste time freeing memory of the last file, the process termination will do that
    flag.let_OS_cleanup_on_exit = is_last_user_txt_file && z_file->z_closes_after_me && !arch_is_valgrind(); 

    DT_FUNC (txt_file, zip_finalize)(is_last_user_txt_file);

    if (flag.show_time_comp_i == flag.zip_comp_i) 
        profiler_add_evb_and_print_report();
}
