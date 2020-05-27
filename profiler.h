// ------------------------------------------------------------------
//   profiler.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PROFILER_INCLUDED
#define PROFILER_INCLUDED

#include "genozip.h"
#include <time.h>
#if defined _MSC_VER // Microsoft compiler
#include "compatibility/visual_c_gettime.h"
#elif defined __APPLE__ 
#include "compatibility/mac_gettime.h"
#endif

typedef struct {
    int64_t wallclock, read, compute, compressor, write, piz_read_one_vb, vcf_piz_get_variant_data_line, 
        vcf_piz_get_haplotype_data_line,
        vcf_piz_initialize_sample_iterators, piz_get_line_subfields, vcf_piz_reconstruct_samples, 
        vcf_piz_get_phase_data_line, vcf_piz_reconstruct_genotype_data_line, zfile_uncompress_section,
        piz_reconstruct_vb, buf_alloc, txtfile_read_header, txtfile_read_vblock,
        seg_all_data_lines, vcf_zip_generate_haplotype_sections, sample_haplotype_data, count_alt_alleles,
        zip_generate_genotype_sections, vcf_zip_generate_phase_sections, zip_generate_variant_data_section,
        mtf_integrate_dictionary_fragment, mtf_clone_ctx, mtf_merge_in_vb_ctx_one_dict_id,
        md5,zfile_compress_dictionary_data,
        lock_mutex_compress_dict, lock_mutex_zf_ctx,
        tmp1, tmp2, tmp3, tmp4, tmp5;
} ProfilerRec;

#ifdef _MSC_VER
typedef struct my_timespec TimeSpecType;
#else
typedef struct timespec TimeSpecType;
#endif

#define START_TIMER     TimeSpecType profiler_timer; \
                        if (flag_show_time) clock_gettime(CLOCK_REALTIME, &profiler_timer); 

#define COPY_TIMER(res) if (flag_show_time) { \
                            TimeSpecType tb; \
                            clock_gettime(CLOCK_REALTIME, &tb); \
                            res += (tb.tv_sec-profiler_timer.tv_sec)*1000000000ULL + (tb.tv_nsec-profiler_timer.tv_nsec); \
                        }

#define PRINT_TIMER(str) { TimeSpecType tb; \
                           clock_gettime(CLOCK_REALTIME, &tb); \
                           printf ("%u.%06u: %s\n", (uint32_t)tb.tv_sec, (uint32_t)(tb.tv_nsec/1000), (str)); }

extern void profiler_add (ProfilerRec *dst, const ProfilerRec *src);
extern const char *profiler_print_short (const ProfilerRec *p);
extern void profiler_print_report (const ProfilerRec *p, unsigned max_threads, unsigned used_threads, const char *filename, unsigned num_vbs);

#endif
