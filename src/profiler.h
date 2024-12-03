// ------------------------------------------------------------------
//   profiler.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#ifdef __APPLE__ 
#include "mac_compat.h"
#endif
#include <time.h>

#define profiled \
        file_open_z, file_close, buf_low_level_free, buflist_find_buf, buflist_sort, buflist_test_overflows_do,\
        read, compute, compressor_bz2, compressor_lzma, compressor_bsc, \
        write, write_fg, write_bg, zriter_write, piz_read_one_vb, vb_get_vb,\
        compressor_domq, compressor_actg, mgzip_uncompress_during_read, igzip_uncompress_during_read, \
        piz_get_line_subfields, b250_zip_generate, zip_generate_local, zip_compress_ctxs, ctx_merge_in_vb_ctx, wait_for_merge,\
        zfile_uncompress_section, codec_assign_best_codec, compressor_pbwt, compressor_longr, compressor_homp, compressor_t0, \
        compressor_rans, compressor_arith, compressor_normq, compressor_pacb, compressor_smux, compressor_oq, \
        codec_domq_reconstruct, codec_domq_reconstruct_dom_run, codec_longr_reconstruct, codec_homp_reconstruct, \
        codec_t0_reconstruct, codec_pacb_reconstruct, codec_smux_reconstruct, codec_oq_reconstruct, \
        reconstruct_vb, buf_alloc_main, buf_alloc_compute, buf_destroy_do_do_main, buf_destroy_do_do_compute, buf_overlay_do, buf_trim_do, \
        buf_free_main, buf_free_compute, buflist_add_buf, buflist_remove_buf, \
        dispatcher_recycle_vbs, sections_create_index, \
        txtfile_discover_specific_gz, txtfile_read_header, txtfile_read_vblock, txtfile_get_unconsumed_callback, fastq_txtfile_sync_to_R1_by_num_lines, \
        txtfile_read_block_mgzip, txtfile_read_block_zlib, txtfile_read_block_igzip, txtfile_read_block_bz2, \
        bgzf_io_thread, bgzf_compute_thread, bgzf_writer_thread, mgzip_uncompress_vb, mgzip_copy_unconsumed_blocks, mgzip_read_block_with_bsize, \
        bgzf_compress_one_block, bgzf_uncompress_one_prescribed_block, \
        mgzip_read_block_no_bsize, \
        zip_modify, vcf_zip_modify, vcf_optimize_samples, vcf_optimize_QUAL, vcf_optimize_INFO, vcf_convert_probabilites_to_phred, \
        vcf_convert_likelihoods_to_phred, vcf_phred_optimize, optimize_float_3_sig_dig, \
        vcf_seg_samples, vcf_seg_copy_one_sample, vcf_seg_one_sample, vcf_seg_info_subfields, vcf_seg_finalize_INFO_fields, vcf_seg_analyze_copied_GT, \
        seg_all_data_lines, seg_get_next_line, seg_get_next_item, seg_initialize,\
        ctx_clone, qname_seg, sam_cigar_seg, sam_seg_BWA_XA_Z, sam_seg_BWA_XA_pos, sam_sa_prim_finalize_ingest, \
        sam_zip_prim_ingest_vb, sam_zip_prim_ingest_idle, sam_zip_prim_ingest_wait_for_seq_mutex, sam_zip_prim_ingest_wait_for_qual_mutex, sam_zip_prim_ingest_wait_for_qname_mutex, sam_zip_prim_ingest_wait_for_aln_mutex, \
        sam_zip_prim_ingest_vb_pack_seq, sam_zip_prim_ingest_vb_compress_qual, sam_zip_prim_ingest_vb_compress_qnames, sam_zip_prim_ingest_solo_data, sam_zip_prim_ingest_vb_create_index, \
        sam_seg_SEQ, sam_seg_verify_saggy_line_SEQ, sam_seg_SEQ_vs_ref, sam_seg_bisulfite_M, reconstruct_SEQ_copy_sag_prim, \
        sam_analyze_copied_SEQ, sam_cigar_special_CIGAR, sam_piz_special_QUAL, \
        sam_sag_by_flag_scan_for_depn, sam_sag_by_flag_scan_sag_depn_index, sam_sag_by_flag_scan_sort_qname_index, scan_remove_single_vb_depns, \
        sam_seg_QUAL, sam_seg_is_gc_line, sam_seg_aux_all, sam_seg_MD_Z_analyze, sam_seg_bsseeker2_XG_Z_analyze,\
        sam_seg_bismark_XM_Z, sam_seg_bsbolt_XB, sam_seg_AS_i, sam_seg_NM_i, sam_seg_SA_Z, sam_seg_BWA_XS_i,\
        sam_seg_TX_AN_Z, sam_seg_barcode_qual, sam_seg_CB_Z, sam_seg_CR_Z, sam_seg_RX_Z, sam_seg_BX_Z, sam_seg_ULTIMA_tp,\
        sam_seg_QX_Z, sam_seg_BC_Z, sam_seg_gene_name_id, sam_seg_fx_Z, sam_seg_other_seq, sam_seg_GR_Z, sam_seg_GY_Z,\
        sam_seg_sag_stuff, sam_cigar_binary_to_textual, squank_seg, bam_seq_to_sam, aligner_seg_seq, sam_header_inspect,\
        sam_header_add_contig, contigs_create_index, sam_header_zip_inspect_PG_lines, sam_header_zip_inspect_RG_lines, sam_header_zip_inspect_HD_line, \
        sam_header_zip_inspect_SQ_lines, cram_inspect_file, \
        sam_deep_zip_merge, sam_piz_con_item_cb, sam_piz_deep_compress, sam_piz_deep_add_qname, sam_piz_deep_add_seq, sam_piz_deep_add_qual,\
        sam_piz_deep_finalize_ents, sam_piz_deep_grab_deep_ents, fastq_seg_find_deep, \
        scan_index_qnames_preprocessing, sam_piz_sam2fastq_QUAL, sam_piz_sam2bam_QUAL,\
        fastq_read_R1_data, piz_read_all_ctxs, fastq_seg_get_lines, fastq_seg_SEQ, fastq_seg_QUAL, \
        fastq_seg_deep, fastq_deep_seg_find_subseq, fastq_seg_DESC, fastq_seg_saux, fastq_seg_deep_consume_unique_matching_ent,\
        fastq_bamass_populate, bamass_read_one_vb, bamass_append_z_ents, bamass_link_entries, bamass_generate, bamass_link,\
        bamass_generate_bamass_ents, fastq_bamass_seg_SEQ, fastq_special_SEQ_by_bamass, fastq_seg_find_bamass, bamass_after_link_entries,\
        bamass_get_one_bam_aln, bamass_get_one_sam_aln, bamass_prepare_cigar, bamass_generate_bamass_ents_hash, fastq_bamass_seg_CIGAR, bamass_segconf, \
        fastq_bamass_retrieve_ent, \
        ref_initialize_ranges,\
        fastq_special_set_deep, fastq_special_deep_copy_QNAME, fastq_special_deep_copy_SEQ, fastq_special_deep_copy_SEQ_by_ref, \
        fastq_special_deep_copy_QUAL, fastq_special_monochar_QUAL, \
        refhash_calc_one_range, refhash_compress_one_vb, refhash_compress_refhash, refhash_load, refhash_uncompress_one_vb, refhash_read_one_vb,\
        txtheader_zip_read_and_compress, txtheader_compress, txtheader_compress_one_fragment, txtheader_piz_read_and_reconstruct,\
        digest, digest_txt_header, ref_make_calculate_digest, refhash_load_digest, ref_load_digest, refhash_compress_digest, \
        dict_io_compress_dictionaries, dict_io_assign_codecs, dict_io_compress_one_fragment, \
        aligner_best_match, aligner_get_word_from_seq, aligner_update_best, aligner_seq_to_bitmap, aligner_first_layer, aligner_additional_layers, \
        refhash_generate_emoneg, ref_contigs_compress,\
        zip_write_global_area, zip_finalize, \
        piz_read_global_area, ref_load_stored_reference, dict_io_read_all_dictionaries, dict_io_build_word_lists, \
        ref_read_one_range, ref_uncompress_one_range, vb_release_vb_do, vb_destroy_vb,\
        sam_load_groups_add_one_prim_vb, sam_load_groups_move_comp_to_zfile, sam_load_groups_move_comp_to_zfile_idle, \
        sam_load_groups_add_qnames, sam_load_groups_add_flags, sam_load_groups_add_seq, sam_load_groups_add_seq_pack, sam_load_groups_add_qual, sam_load_groups_add_cigars, \
        sam_load_groups_add_SA_alns, sam_load_groups_add_solo_data, sam_load_groups_add_grps, \
        sam_zip_calculate_max_conc_writing_vbs, sam_reconstruct_SEQ_vs_ref, sam_reconstruct_SEQ_get_textual_ref, aligner_reconstruct_seq, sam_piz_sam2bam_SEQ,\
        sam_piz_special_SEQ, reconstruct_SEQ_copy_saggy, sam_piz_special_MD, \
        writer_main_loop, writer_create_plan, gencomp_piz_initialize_vb_info, gencomp_piz_update_reading_list, gencomp_piz_vb_to_plan, \
        sam_bismark_piz_update_meth_call, sam_gencomp_trim_memory, \
        zip_handle_unique_words_ctxs, random_access_merge_in_vb, \
        vcf_seg_PROBE_A, vcf_seg_QUAL, \
        random_access_finalize_entries, random_access_compress, ctx_compress_counts, zfile_compress_genozip_header,\
        ref_compress_ref, ref_compress_one_range, ref_copy_compressed_sections_from_reference_file,\
        piz_main_loop_idle, zip_main_loop_idle, zip_free_undeeded_zctx_bufs_after_seg, \
        piz_uncompress_all_ctxs__recon, piz_uncompress_all_ctxs__fasta_writer_init, piz_uncompress_all_ctxs__fastq_read_r1, piz_uncompress_all_ctxs__sam_load_sag,\
        gencomp_absorb_vb_gencomp_lines, gencomp_flush, gencomp_offload_DEPN_to_disk, gencomp_reread_lines_as_prescribed, gencomp_do_offload_write, \
        compress_depn_buf,  \
        tmp1, tmp2, tmp3, tmp4, tmp5, \
        fields[MAX_DICTS]/* ZIP: compression time (all ctxs); PIZ: recon time (fields only). must be last for profiler_add. */ \

typedef struct {
        struct { int64_t profiled; } nanosecs; 
        struct { int64_t profiled; } count;
        rom next_name, next_subname;
        unsigned num_vbs, max_vb_size_mb, num_txt_files;
        float avg_compute_vbs[MAX_NUM_TXT_FILES_IN_ZFILE];  // ZIP/PIZ: average number of compute threads active at any given time during the lifetime of the ZIP/PIZ dispatcher
} ProfilerRec;

typedef struct timespec TimeSpecType;

#define HAS_SHOW_TIME(vb) (__builtin_expect(flag.show_time_comp_i != COMP_NONE/*fail fast*/, false) && /* __builtin_expect to reduce overhead in normal execution without --show-time */ \
                           (flag.show_time_comp_i == COMP_ALL || flag.show_time_comp_i == (vb)->comp_i))

#define START_TIMER TimeSpecType profiler_timer; \
                    if (__builtin_expect(flag.show_time_comp_i != COMP_NONE, false)) clock_gettime(CLOCK_REALTIME, &profiler_timer); 

#define START_TIMER_ALWAYS TimeSpecType profiler_timer; clock_gettime(CLOCK_REALTIME, &profiler_timer); 

#define CHECK_TIMER ({ TimeSpecType tb; \
                       clock_gettime(CLOCK_REALTIME, &tb); \
                       ((uint64_t)((tb).tv_sec-(profiler_timer).tv_sec))*1000000000ULL + ((int64_t)(tb).tv_nsec-(int64_t)(profiler_timer).tv_nsec); })

#define COPY_TIMER_FULL(vb,res,atomic) { /* str - print in case of specific show-time=<res> */ \
    if (HAS_SHOW_TIME(vb)) { \
        uint64_t delta = CHECK_TIMER; \
        if (flag.show_time[0] && strstr (#res, flag.show_time)) { \
            iprintf ("%s %s%s%s: %"PRIu64" μsec\n", #res, \
                     ((vb)->profile.next_name    ? (vb)->profile.next_name : ""),\
                     ((vb)->profile.next_subname ? "." : ""),\
                     ((vb)->profile.next_subname ? (vb)->profile.next_subname : ""),\
                     delta/1000);\
            (vb)->profile.next_name = (vb)->profile.next_subname = NULL;\
        }\
        if (atomic) { \
            __atomic_fetch_add (&(vb)->profile.nanosecs.res, delta, __ATOMIC_RELAXED); \
            increment_relaxed ((vb)->profile.count.res); \
        } \
        else { \
            (vb)->profile.nanosecs.res += delta; \
            (vb)->profile.count.res++; \
        } \
    } \
}

#define COPY_TIMER(res) COPY_TIMER_FULL(vb, res, false)
#define COPY_TIMER_EVB(res) ({ if (evb) COPY_TIMER_FULL(evb, res, true); })
#define COPY_TIMER_COMPRESS(res)  ({ if (!in_assign_codec) COPY_TIMER(res); }) // account only if not running from codec_assign_best_codec, because it accounts for itself

#define PAUSE_TIMER(vb) \
    TimeSpecType on_hold_timer; \
    rom save_name=0, save_subname=0; \
    if (HAS_SHOW_TIME(vb)) { \
        clock_gettime(CLOCK_REALTIME, &on_hold_timer); \
        save_name = vb->profile.next_name; \
        save_subname = vb->profile.next_subname; \
    }

#define RESUME_TIMER(vb,res)\
    if (HAS_SHOW_TIME(vb)) { \
        TimeSpecType tb; \
        clock_gettime(CLOCK_REALTIME, &tb); \
        vb->profile.next_name = save_name; \
        vb->profile.next_subname = save_subname; \
        (vb)->profile.nanosecs.res -= (tb.tv_sec-on_hold_timer.tv_sec)*1000000000ULL + ((int64_t)tb.tv_nsec-(int64_t)on_hold_timer.tv_nsec); \
    } 

#define PRINT_TIMER(str) { TimeSpecType tb; \
                           clock_gettime(CLOCK_REALTIME, &tb); \
                           iprintf ("%u.%06u: %s\n", (uint32_t)tb.tv_sec, (uint32_t)(tb.tv_nsec/1000), (str)); }

extern void profiler_new_z_file (void);
extern void profiler_initialize (void);
extern void profiler_add (ConstVBlockP vb);
extern rom profiler_print_short (const ProfilerRec *p);
extern void profiler_add_evb_and_print_report (void);

extern void profiler_set_avg_compute_vbs (float avg_compute_vbs);
extern StrTextSuperLong profiler_get_avg_compute_vbs (char sep);

