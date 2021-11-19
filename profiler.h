// ------------------------------------------------------------------
//   profiler.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

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
        sam_seg_SEQ, compressor_domq, compressor_actg, bgzf_io_thread, bgzf_compute_thread,
        piz_get_line_subfields, zip_generate_b250, zip_generate_local, zip_compress_ctxs, ctx_merge_in_vb_ctx,
        zfile_uncompress_section, codec_assign_best_codec, compressor_pbwt, compressor_enano, compressor_rans, compressor_arith,
        reconstruct_vb, buf_alloc, dispatcher_recycle_vbs, txtfile_read_header, txtfile_read_vblock,
        seg_all_data_lines, codec_hapmat_count_alt_alleles, seg_initialize,
        ctx_clone, qname_seg, sam_cigar_seg, sam_seg_XA_pos, 
        md5,ctx_compress_one_dict_fragment, aligner_best_match, aligner_get_word_from_seq,
        aligner_get_match_len, generate_rev_complement_genome, ref_contigs_compress,
        linesorter_compress_qsort, linesorter_compress_recon_plan, 
        piz_read_global_area, ref_load_stored_reference, ctx_read_all_dictionaries, ctx_dict_build_word_lists, 
        ref_read_one_range, ref_uncompress_one_range, vb_release_vb_do, vb_destroy_vb,
        wait_for_vb_1_mutex,
        tmp1, tmp2, tmp3, tmp4, tmp5;
        const char *next_name, *next_subname;
        unsigned num_vbs, max_vb_size_mb;
} ProfilerRec;

typedef struct timespec TimeSpecType;

#define START_TIMER TimeSpecType profiler_timer; \
                    if (flag.show_time) clock_gettime(CLOCK_REALTIME, &profiler_timer); 

#define START_TIMER_ALWAYS TimeSpecType profiler_timer; clock_gettime(CLOCK_REALTIME, &profiler_timer); 

#define CHECK_TIMER ({ TimeSpecType tb; \
                       clock_gettime(CLOCK_REALTIME, &tb); \
                       ((uint64_t)((tb).tv_sec-(profiler_timer).tv_sec))*1000000000ULL + ((tb).tv_nsec-(profiler_timer).tv_nsec); })

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

#define PAUSE_TIMER(vb) \
    TimeSpecType on_hold_timer; \
    const char *save_name=0, *save_subname=0; \
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
        (vb)->profile.res -= (tb.tv_sec-on_hold_timer.tv_sec)*1000000000ULL + (tb.tv_nsec-on_hold_timer.tv_nsec); \
    } 

#define PRINT_TIMER(str) { TimeSpecType tb; \
                           clock_gettime(CLOCK_REALTIME, &tb); \
                           iprintf ("%u.%06u: %s\n", (uint32_t)tb.tv_sec, (uint32_t)(tb.tv_nsec/1000), (str)); }

extern void profiler_initialize (void);
extern void profiler_add (ConstVBlockP vb);
extern const char *profiler_print_short (const ProfilerRec *p);
extern void profiler_print_report (void);

