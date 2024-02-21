// ------------------------------------------------------------------
//   profiler.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "vblock.h"
#include "profiler.h"
#include "flags.h"
#include "file.h"
#include "segconf.h"
#include "mutex.h"

static ProfilerRec profile = {};    // data for this z_file 
static Mutex profile_mutex = {};
static TimeSpecType profiler_timer; // wallclock

void profiler_initialize (void)
{
    mutex_bottleneck_analysis_init();
    mutex_initialize (profile_mutex);
}

void profiler_new_z_file (void)
{
    memset (&profile, 0, sizeof (profile));
    clock_gettime (CLOCK_REALTIME, &profiler_timer); // initialze wallclock
}

void profiler_set_avg_compute_vbs (float avg_compute_vbs)
{
    ASSERT0 (profile.num_txt_files >= 0 && profile.num_txt_files < MAX_NUM_TXT_FILES_IN_ZFILE, "too many txt files");

    profile.avg_compute_vbs[profile.num_txt_files++] = avg_compute_vbs;
}


StrTextSuperLong profiler_get_avg_compute_vbs (char sep)
{
    StrTextSuperLong s = {};
    for (int i=0; i < profile.num_txt_files; i++)
        sprintf (&s.s[strlen(s.s)], "%.1f%c", profile.avg_compute_vbs[i], sep);

    if (s.s[0]) s.s[strlen(s.s)-1] = 0; // remove final separator

    return s;
}

void profiler_add (ConstVBlockP vb)
{
    mutex_lock (profile_mutex);

    if (Ltxt) {
        profile.num_vbs++;
        profile.max_vb_size_mb = MAX_(profile.max_vb_size_mb, segconf.vb_size >> 20);
    }

    int num_profiled = sizeof (profile.nanosecs) / sizeof (uint64_t) - MAX_DICTS + 
                       (IS_PIZ ? vb->num_contexts : 0); 
    
    for (int i=0; i < num_profiled; i++) 
        if (((uint64_t *)&vb->profile.count)[i]) {
            ((uint64_t *)&profile.nanosecs)[i] += ((uint64_t *)&vb->profile.nanosecs)[i];
            ((uint64_t *)&profile.count)[i]    += ((uint64_t *)&vb->profile.count)[i];
        }

    // ZIP: add compressor data by zctx, while collected by vctx
    if (IS_ZIP) 
        for (int v=num_profiled; v < num_profiled + vb->num_contexts; v++) 
            if (((uint64_t *)&vb->profile.count)[v]) {
                ContextP zctx = ctx_get_zctx_from_vctx (CTX(v-num_profiled), false, true);
                if (!zctx) continue; // should never happen

                int z = zctx->did_i + num_profiled;
                ((uint64_t *)&profile.nanosecs)[z] += ((uint64_t *)&vb->profile.nanosecs)[v];
                ((uint64_t *)&profile.count)[z]    += ((uint64_t *)&vb->profile.count)[v];
            }

    mutex_unlock (profile_mutex);
}

static inline uint32_t ms(uint64_t ns) { return (uint32_t)(ns / 1000000);}

rom profiler_print_short (const ProfilerRec *p)
{
    static char str[300]; // not thread safe
    sprintf (str, "read: %s compute:%s write: %s", str_int_commas (ms(p->nanosecs.read)).s, str_int_commas (ms(p->nanosecs.compute)).s, str_int_commas (ms(p->nanosecs.write)).s);
    return str;
}

void profiler_add_evb_and_print_report (void)
{
    profiler_add (evb);

    static rom space = "                                                   ";
#   define PRINT_(x, label, level) if (profile.nanosecs.x)          \
        iprintf ("%.*s %s: %s (N=%s)\n", (level)*3, space, (label), \
                 str_int_commas (ms(profile.nanosecs.x)).s,         \
                 str_int_commas (profile.count.x).s);

#   define PRINT(x, level) PRINT_(x, #x, (level)) 
    
    rom os = flag.is_windows ? "Windows"
           : flag.is_mac     ? "MacOS"
           : flag.is_linux   ? "Linux"
           :                   "Unknown OS";

    iprintf ("\n%s PROFILER:\n", IS_ZIP ? "ZIP" : "PIZ");
    iprintf ("OS=%s\n", os);
    iprintf ("Build=%s\n", flag.debug ? "Debug" : "Optimized");

    iprintf ("Wallclock: %s milliseconds\n", str_int_commas (ms (CHECK_TIMER)).s);

    if (command == PIZ) { // this is a uncompress operation

        iprint0 ("GENOUNZIP main thread (piz_one_txt_file):\n");
        PRINT (ref_load_stored_reference, 1);
        PRINT (ref_initialize_ranges, 2);
        PRINT (ref_read_one_range, 2);
        PRINT (ref_uncompress_one_range, 2);
        PRINT (piz_read_global_area, 1); // sometimes also includes ref_load_stored_reference, but usually not
        PRINT (dict_io_read_all_dictionaries, 2);
        PRINT (dict_io_build_word_lists, 3);
        PRINT (txtheader_piz_read_and_reconstruct, 1);
        PRINT (digest_txt_header, 2);
        PRINT (vb_get_vb, 1);
        PRINT (piz_read_one_vb, 1);
        PRINT (read, 2);
        PRINT (bgzf_io_thread, 1);
        PRINT (bgzf_writer_thread, 1);
        PRINT (write, 1);
        PRINT (sam_sa_prim_finalize_ingest, 1);
        PRINT (sam_piz_deep_grab_deep_ents, 1);
        PRINT (sam_piz_deep_finalize_ents, 2);
        PRINT (piz_main_loop_idle, 1);

        iprintf ("GENOUNZIP compute threads: %s\n", str_int_commas (ms(profile.nanosecs.compute)).s);
        PRINT (zfile_uncompress_section, 1);
        PRINT (compressor_bz2,   2);
        PRINT (compressor_lzma,  2);
        PRINT (compressor_bsc,   2);
        PRINT (compressor_rans,  2);
        PRINT (compressor_arith, 2);
        PRINT (compressor_domq,  2);
        PRINT (compressor_normq, 2);
        PRINT (compressor_actg,  2);
        PRINT (compressor_pbwt,  2);
        PRINT (compressor_longr, 2);
        PRINT (compressor_homp,  2);
        PRINT (compressor_t0,    2);
        PRINT (compressor_pacb,  2);
        PRINT (compressor_smux,  2);
        
        PRINT (reconstruct_vb, 1);
        for (Did did_i=0; did_i < z_file->num_contexts; did_i++) 
            PRINT_(fields[did_i], ZCTX(did_i)->tag_name, 2);

        PRINT (sam_reconstruct_SEQ_vs_ref, 2);
        PRINT (reconstruct_SEQ_copy_sag_prim, 3);
        PRINT (sam_analyze_copied_SEQ, 3);
        PRINT (sam_bismark_piz_update_meth_call, 3);
        PRINT (aligner_reconstruct_seq, 2);
        PRINT (sam_piz_special_QUAL, 2);
        if (Z_DT(SAM) || Z_DT(BAM)) {
            PRINT (codec_longr_reconstruct,3);
            PRINT (codec_homp_reconstruct, 3);
            PRINT (codec_t0_reconstruct,   3);
            PRINT (codec_smux_reconstruct, 3);
            PRINT (codec_pacb_reconstruct, 3);
            PRINT (codec_domq_reconstruct, 3);
            PRINT (codec_domq_reconstruct_dom_run, 4);
        }
        PRINT (fastq_special_monochar_QUAL, 2);        
        PRINT (sam_piz_sam2fastq_QUAL, 2); 
        PRINT (sam_piz_sam2bam_QUAL, 2);        
        PRINT (sam_cigar_special_CIGAR, 2);

        PRINT (sam_piz_con_item_cb, 2); 
        PRINT (sam_piz_deep_add_qname, 3); 
        PRINT (sam_piz_deep_add_seq, 3); 
        PRINT (sam_piz_deep_add_qual, 3);
        PRINT (sam_piz_deep_compress, 4); // mostly under qual, a tiny bit of seq

        PRINT (fastq_special_set_deep, 2); 
        PRINT (fastq_special_deep_copy_QNAME, 2);
        PRINT (fastq_special_deep_copy_SEQ, 2);
        PRINT (fastq_special_deep_copy_QUAL, 2);
         
        PRINT (sam_zip_prim_ingest_vb, 1);
        PRINT (digest, 1); // note: in SAM/BAM digest is done in the writer thread, otherwise its done in the compute thread. TODO: change level to 0 in case of SAM/BAM
        PRINT (piz_get_line_subfields, 2);
        PRINT (codec_hapmat_piz_get_one_line, 2);
        PRINT (sam_load_groups_add_one_prim_vb, 1);
        if (Z_DT(FASTQ)) {
            PRINT (codec_longr_reconstruct, 3);
            PRINT (codec_domq_reconstruct, 2);
            PRINT (codec_domq_reconstruct_dom_run, 3);
            PRINT (codec_homp_reconstruct, 2);
            PRINT (codec_smux_reconstruct, 2);
            PRINT (codec_pacb_reconstruct, 2);
        }
        
        if (profile.nanosecs.bgzf_compute_thread) {
            iprintf ("GENOUNZIP BGZF threads: %s\n", str_int_commas (ms(profile.nanosecs.bgzf_compute_thread)).s);
            PRINT (bgzf_compute_thread, 1);
            PRINT (bgzf_compress_one_block, 2);
        }
    }

    else { // compress
        iprint0 ("GENOZIP main thread (zip_one_file):\n");
        PRINT (ref_load_stored_reference, 1);
        PRINT (ref_read_one_range, 2);
        PRINT (ref_uncompress_one_range, 2);
        PRINT (ref_load_digest, 2); 
        PRINT (refhash_load, 1);
        PRINT (refhash_load_digest, 2);
        PRINT (refhash_read_one_vb, 2);
        PRINT (refhash_compress_digest, 2); // make-ref
        PRINT (refhash_uncompress_one_vb, 2);
        PRINT (txtheader_zip_read_and_compress, 1);
        PRINT (txtfile_read_header, 2);
        PRINT (sam_header_add_contig, 2); 
        PRINT (contigs_create_index, 2);
        PRINT (sam_header_zip_inspect_PG_lines, 2); 
        PRINT (sam_header_zip_inspect_HD_line, 2);
        PRINT (ref_initialize_ranges, 2);
        PRINT (txtheader_compress, 2);
        PRINT (txtheader_compress_one_fragment, 3); 
        PRINT (digest_txt_header, 2);
        PRINT (vb_get_vb, 1);
        PRINT (fastq_read_pair_1_data, 1);
        PRINT (piz_read_all_ctxs, 2);
        PRINT (txtfile_read_vblock, 1);
        PRINT (read, 2);
        PRINT (txtfile_read_block_zlib, 3);
        PRINT (txtfile_read_block_gz, 3);
        PRINT (txtfile_read_block_bz2, 3);
        PRINT (txtfile_read_block_bgzf, 3);
        PRINT (bgzf_read_block, 4);
        PRINT (txtfile_read_block_bgzf_uncompress, 4);
        PRINT (fastq_txtfile_have_enough_lines, 2);
        PRINT (txtfile_get_unconsumed_to_pass_to_next_vb, 2);
        PRINT (bgzf_copy_unconsumed_blocks, 2);
        PRINT (zriter_write, 1);
        PRINT (write_fg, 2);
        PRINT (write_bg, 2);
        PRINT (bgzf_io_thread, 1);
        PRINT (sam_sa_prim_finalize_ingest, 1);
        PRINT (zip_main_loop_idle, 1);
        PRINT (zip_free_undeeded_zctx_bufs_after_seg, 1);
        PRINT (generate_recon_plan, 1);
        PRINT (sam_zip_recon_plan_add_gc_lines, 2);
        PRINT (sam_zip_recon_plan_count_writers, 3);
        PRINT (recon_plan_compress, 2);
        PRINT (recon_plan_deltify, 3);
        PRINT (recon_plan_compress_one_fragment, 3);
        PRINT (zip_write_global_area, 1);
        PRINT (dict_io_compress_dictionaries, 2); 
        PRINT (dict_io_assign_codecs, 3); 
        PRINT (dict_io_compress_one_fragment, 3);
        PRINT (ref_compress_ref, 2);
        PRINT (ref_compress_one_range, 3);
        PRINT (refhash_calc_one_range, 4);
        PRINT (refhash_compress_refhash, 2);
        PRINT (refhash_compress_one_vb, 3);
        PRINT (refhash_compress_digest, 3);
        PRINT (ref_make_calculate_digest, 3);
        PRINT (ref_contigs_compress, 3);
        PRINT (ref_copy_compressed_sections_from_reference_file, 3);
        PRINT (random_access_finalize_entries, 2);
        PRINT (random_access_compress, 2);
        PRINT (ctx_compress_counts, 2);
        PRINT (zfile_compress_genozip_header, 2);

        iprintf ("GENOZIP compute threads %s\n", str_int_commas (ms(profile.nanosecs.compute)).s);
        PRINT (bgzf_uncompress_vb, 1);
        PRINT (ctx_clone, 1);
        PRINT (scan_index_qnames_preprocessing, 1);
        PRINT (seg_all_data_lines, 1);
        PRINT (seg_initialize, 2);
        PRINT (qname_seg, 2);
        PRINT (sam_cigar_seg, 2);
        PRINT (squank_seg, 3);
        PRINT (fastq_seg_get_lines, 2);
        PRINT (seg_get_next_line, 3);
        PRINT (seg_get_next_item, 3);
        PRINT (fastq_seg_deep, 2);
        PRINT (fastq_deep_seg_find_subseq, 3);
        PRINT (fastq_seg_deep_consume_unique_matching_ent, 3);
        PRINT (fastq_seg_SEQ, 2);
        PRINT (sam_seg_SEQ, 2);
        PRINT (sam_seg_SEQ_vs_ref, 3);
        PRINT (sam_seg_bisulfite_M, 4);
        PRINT (sam_seg_verify_saggy_line_SEQ, 3); 
        PRINT (sam_analyze_copied_SEQ, 3);
        PRINT (aligner_seg_seq, 3);
        PRINT (aligner_best_match, 4);
        PRINT (aligner_additional_layers, 5);
        PRINT (aligner_update_best, 5);
        PRINT (aligner_get_word_from_seq, 5);
        PRINT (aligner_seq_to_bitmap, 5);
        PRINT (fastq_seg_DESC, 2);
        PRINT (fastq_seg_saux, 2);
        PRINT (bam_seq_to_sam, 2);
        PRINT (fastq_seg_QUAL, 2);
        PRINT (sam_seg_QUAL, 2);
        PRINT (sam_seg_is_gc_line, 2);
        PRINT (sam_seg_MD_Z_analyze, 2);
        PRINT (sam_cigar_binary_to_textual, 2);
        PRINT (sam_seg_bsseeker2_XG_Z_analyze, 2);
        PRINT (sam_seg_sag_stuff, 2);
        PRINT (sam_seg_aux_all, 2);
        PRINT (sam_seg_SA_Z, 3);
        PRINT (sam_seg_AS_i, 3);
        PRINT (sam_seg_NM_i, 3);                       
        PRINT (sam_seg_BWA_XA_Z, 3);
        PRINT (sam_seg_BWA_XA_pos, 4);
        PRINT (sam_seg_BWA_XS_i, 3);
        PRINT (sam_seg_bismark_XM_Z, 3);
        PRINT (sam_seg_bsbolt_XB, 3);
        PRINT (sam_seg_TX_AN_Z, 3);
        PRINT (sam_seg_barcode_qual, 3);
        PRINT (sam_seg_CB_Z, 3);
        PRINT (sam_seg_CR_Z, 3);
        PRINT (sam_seg_RX_Z, 3); 
        PRINT (sam_seg_BX_Z, 3);
        PRINT (sam_seg_QX_Z, 3);
        PRINT (sam_seg_BC_Z, 3);
        PRINT (sam_seg_gene_name_id, 3);
        PRINT (sam_seg_fx_Z, 3);
        PRINT (sam_seg_other_seq, 3);
        PRINT (sam_seg_GR_Z, 3);
        PRINT (sam_seg_GY_Z, 3);
        PRINT (sam_seg_ULTIMA_tp, 3);
        PRINT (vcf_seg_PROBE_A, 2);
        PRINT (random_access_merge_in_vb, 1); 
        PRINT (gencomp_absorb_add_to_queue, 1);
        PRINT (gencomp_flush, 2);
        PRINT (gencomp_offload_DEPN_to_disk, 3);
        PRINT (gencomp_reread_lines_as_prescribed, 1);
        PRINT (bgzf_uncompress_one_prescribed_block, 2);
        PRINT (ctx_merge_in_vb_ctx, 1);
        PRINT (wait_for_merge, 2);
        PRINT (sam_deep_merge, 2);
        PRINT (zip_compress_ctxs, 1);
        PRINT (b250_zip_generate, 2);
        PRINT (zip_generate_local, 2);
        PRINT (codec_assign_best_codec, 2);
        
        PRINT (compressor_bz2,   2);
        PRINT (compressor_lzma,  2);
        PRINT (compressor_bsc,   2);
        PRINT (compressor_rans,  2);
        PRINT (compressor_arith, 2);
        PRINT (compressor_domq,  2);
        PRINT (compressor_normq, 2);
        PRINT (compressor_actg,  2);
        PRINT (compressor_pbwt,  2);
        PRINT (compressor_longr, 2);
        PRINT (compressor_homp,  2);
        PRINT (compressor_t0,    2);
        PRINT (compressor_pacb,  2);
        PRINT (compressor_smux,  2);

        for_zctx 
            PRINT_(fields[zctx->did_i], zctx->tag_name, 2);

        PRINT (sam_zip_prim_ingest_vb, 1);
        PRINT (digest, 1);
    }    
    
    // ZIP and PIZ compute thread
    iprint0 ("COMPUTE THREADS ADMIN\n");
    PRINT (buf_alloc_compute, 1);
    PRINT (buflist_add_buf, 2);
    PRINT (buf_destroy_do_do_compute, 1);
    PRINT (buf_free_compute, 1);
    PRINT (buflist_remove_buf, 2)
    PRINT (buf_overlay_do, 1);

    PRINT (file_open_z, 0);
    PRINT (file_close, 0);
    PRINT (buf_alloc_main, 0);
    PRINT (buf_destroy_do_do_main, 0);
    PRINT (buflist_test_overflows_do, 0);
    PRINT (buflist_sort, 0);
    PRINT (sections_create_index, 0);
    PRINT (buflist_find_buf, 0);
    PRINT (buf_low_level_free, 0);
    PRINT (vb_release_vb_do, 0);
    PRINT (vb_destroy_vb, 0);
    PRINT (dispatcher_recycle_vbs, 0);
    PRINT (refhash_generate_emoneg, 0);
    PRINT (tmp1, 0); 
    PRINT (tmp2, 0); 
    PRINT (tmp3, 0); 
    PRINT (tmp4, 0); 
    PRINT (tmp5, 0); 

    if (profile.num_vbs) {
        iprint0 ("\nVblock stats:\n");
        iprintf ("  Vblocks: %u\n", profile.num_vbs);
        iprintf ("  Maximum vblock size: %u MB\n", profile.max_vb_size_mb);
        iprintf ("  Average number of VBs in compute (per txt file): %s\n", profiler_get_avg_compute_vbs(',').s); // average during the lifetime of the ZIP/PIZ dispatcher, i.e. excluding global_area time etc
        iprintf ("  Average read time: %u ms\n", ms(profile.nanosecs.read) / profile.num_vbs);
        iprintf ("  Average compute time: %u ms\n", ms(profile.nanosecs.compute) / profile.num_vbs);
        iprintf ("  Average write time: %u ms\n", ms(profile.nanosecs.write) / profile.num_vbs);
    }
    
    iprint0 ("\n\n");

    mutex_show_bottleneck_analsyis();
}
