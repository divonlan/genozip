// ------------------------------------------------------------------
//   profiler.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include <time.h>
#ifdef __APPLE__ 
#include "compatibility/mac_gettime.h"
#endif

typedef struct {
    int64_t file_open, file_close, buf_low_level_free, buf_remove_from_buffer_list,
        read, compute, compressor_bz2, compressor_lzma, compressor_bsc, 
        write, piz_read_one_vb, codec_hapmat_piz_get_one_line, vb_get_vb, buf_mmap_do,
        compressor_domq, compressor_actg, bgzf_io_thread, bgzf_compute_thread, bgzf_writer_thread,
        piz_get_line_subfields, zip_generate_b250, zip_generate_local, zip_compress_ctxs, ctx_merge_in_vb_ctx,
        zfile_uncompress_section, codec_assign_best_codec, compressor_pbwt, compressor_longr, 
        compressor_rans, compressor_arith, compressor_normq,
        codec_domq_reconstruct, codec_domq_reconstruct_dom_run, codec_longr_reconstruct,
        reconstruct_vb, buf_alloc, dispatcher_recycle_vbs, txtfile_read_header, txtfile_read_vblock,
        seg_all_data_lines, codec_hapmat_count_alt_alleles, seg_initialize,
        ctx_clone, qname_seg, sam_cigar_seg, sam_seg_BWA_XA_Z, sam_seg_BWA_XA_pos, sam_sa_prim_finalize_ingest, sam_zip_prim_ingest_vb,
        sam_seg_SEQ, sam_seg_verify_saggy_line_SEQ, sam_seg_SEQ_vs_ref, sam_seg_bisulfite_M, reconstruct_SEQ_copy_sag_prim, 
        sam_analyze_copied_SEQ, sam_cigar_special_CIGAR, sam_piz_special_QUAL,
        sam_seg_QUAL, sam_seg_is_gc_line, sam_seg_aux_all, sam_seg_MD_Z_analyze, sam_seg_bsseeker2_XG_Z_analyze,
        sam_seg_bismark_XM_Z, sam_seg_bsbolt_XB, sam_seg_AS_i, sam_seg_NM_i, sam_seg_SA_Z, sam_seg_BWA_XS_i,
        sam_seg_TX_AN_Z, sam_seg_barcode_qual, sam_seg_CB_Z, sam_seg_CR_Z, sam_seg_RX_Z, sam_seg_BX_Z,
        sam_seg_QX_Z, sam_seg_BC_Z, sam_seg_gene_name_id, sam_seg_fx_Z, sam_seg_other_seq, sam_seg_GR_Z, sam_seg_GY_Z,
        scan_index_qnames_preprocessing, sam_piz_sam2fastq_QUAL, sam_piz_sam2bam_QUAL,
        fastq_seg_SEQ, fastq_seg_QUAL,
        sam_seg_sag_stuff, sam_cigar_binary_to_textual, squank_seg, bam_seq_to_sam, aligner_seg_seq,
        sam_header_add_contig, contigs_create_index, sam_header_zip_inspect_PG_lines, sam_header_zip_inspect_HD_line, ref_initialize_ranges,
        sam_deep_merge,
        txtheader_zip_read_and_compress, txtheader_compress, txtheader_compress_one_fragment, txtheader_piz_read_and_reconstruct,
        digest, digest_txt_header,
        dict_io_compress_dictionaries, dict_io_assign_codecs, dict_io_compress_one_fragment, 
        aligner_best_match, aligner_get_word_from_seq, aligner_get_match_len, 
        generate_rev_complement_genome, ref_contigs_compress,
        vcf_linesort_compress_qsort, generate_recon_plan, 
        piz_read_global_area, ref_load_stored_reference, dict_io_read_all_dictionaries, dict_io_build_word_lists, 
        ref_read_one_range, ref_uncompress_one_range, vb_release_vb_do, vb_destroy_vb,
        sam_load_groups_add_one_prim_vb, recon_plan_compress, recon_plan_compress_one_fragment,
        sam_zip_recon_plan_add_gc_lines, sam_zip_gc_calc_depn_vb_info, sam_reconstruct_SEQ_vs_ref, aligner_reconstruct_seq,
        sam_bismark_piz_update_meth_call,
        zip_handle_unique_words_ctxs, ctx_sort_dictionaries_vb_1, random_access_merge_in_vb, gencomp_absorb_vb_gencomp_lines,
        vcf_linesort_merge_vb, vcf_seg_PROBE_A,
        random_access_finalize_entries, random_access_compress, ctx_compress_counts, zfile_compress_genozip_header,
        ref_compress_ref, ref_compress_one_range, ref_copy_compressed_sections_from_reference_file,
        tmp1, tmp2, tmp3, tmp4, tmp5;
        rom next_name, next_subname;
        unsigned num_vbs, max_vb_size_mb;
} ProfilerRec;

typedef struct timespec TimeSpecType;

#define START_TIMER TimeSpecType profiler_timer; \
                    if (flag.show_time) clock_gettime(CLOCK_REALTIME, &profiler_timer); 

#define START_TIMER_ALWAYS TimeSpecType profiler_timer; clock_gettime(CLOCK_REALTIME, &profiler_timer); 

#define CHECK_TIMER ({ TimeSpecType tb; \
                       clock_gettime(CLOCK_REALTIME, &tb); \
                       ((uint64_t)((tb).tv_sec-(profiler_timer).tv_sec))*1000000000ULL + ((int64_t)(tb).tv_nsec-(int64_t)(profiler_timer).tv_nsec); })

#define COPY_TIMER_FULL(vb,res) { /* str - print in case of specific show-time=<res> */ \
    if (flag.show_time) { \
        uint64_t delta = CHECK_TIMER; \
        if (flag.show_time[0] && strstr (#res, flag.show_time)) { \
            iprintf ("%s %s%s%s: %"PRIu64" microsec\n", #res, \
                     ((vb)->profile.next_name    ? (vb)->profile.next_name : ""),\
                     ((vb)->profile.next_subname ? "." : ""),\
                     ((vb)->profile.next_subname ? (vb)->profile.next_subname : ""),\
                     delta/1000);\
            (vb)->profile.next_name = (vb)->profile.next_subname = NULL;\
        }\
        (vb)->profile.res += delta; \
    } \
}

#define COPY_TIMER(res)       COPY_TIMER_FULL(vb, res)
#define COPY_TIMER_VB(vb,res) COPY_TIMER_FULL((vb), res)
#define COPY_TIMER_COMPRESS(res)  ({ if (!in_assign_codec) COPY_TIMER(res); }) // account only if not running from codec_assign_best_codec, because it accounts for itself

#define PAUSE_TIMER(vb) \
    TimeSpecType on_hold_timer; \
    rom save_name=0, save_subname=0; \
    if (flag.show_time) { \
        clock_gettime(CLOCK_REALTIME, &on_hold_timer); \
        save_name = vb->profile.next_name; \
        save_subname = vb->profile.next_subname; \
    }

#define RESUME_TIMER(vb,res)\
    if (flag.show_time) { \
        TimeSpecType tb; \
        clock_gettime(CLOCK_REALTIME, &tb); \
        vb->profile.next_name = save_name; \
        vb->profile.next_subname = save_subname; \
        (vb)->profile.res -= (tb.tv_sec-on_hold_timer.tv_sec)*1000000000ULL + ((int64_t)tb.tv_nsec-(int64_t)on_hold_timer.tv_nsec); \
    } 

#define PRINT_TIMER(str) { TimeSpecType tb; \
                           clock_gettime(CLOCK_REALTIME, &tb); \
                           iprintf ("%u.%06u: %s\n", (uint32_t)tb.tv_sec, (uint32_t)(tb.tv_nsec/1000), (str)); }

extern void profiler_initialize (void);
extern void profiler_add (ConstVBlockP vb);
extern rom profiler_print_short (const ProfilerRec *p);
extern void profiler_print_report (void);

