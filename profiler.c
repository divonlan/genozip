// ------------------------------------------------------------------
//   profiler.c
//   Copyright (C) 2019 Divon Lan <genozip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <stdlib.h>
#include <pthread.h> 

#include "genozip.h"

void profiler_add (ProfilerRec *dst, const ProfilerRec *src)
{
    dst->compressor                        += src->compressor;
    dst->piz_uncompress_variant_block      += src->piz_uncompress_variant_block;
    dst->read                              += src->read;
    dst->write                             += src->write;
    dst->piz_reconstruct_line_components   += src->piz_reconstruct_line_components;
    dst->piz_get_variant_data_line         += src->piz_get_variant_data_line;
    dst->piz_decode_pos                    += src->piz_decode_pos;
    dst->piz_get_haplotype_data_line       += src->piz_get_haplotype_data_line;
    dst->piz_get_phase_data_line           += src->piz_get_phase_data_line;
    dst->piz_get_genotype_data_line        += src->piz_get_genotype_data_line;
    dst->zfile_uncompress_section          += src->zfile_uncompress_section;
    dst->squeeze                           += src->squeeze;
    dst->buffer                            += src->buffer;
    dst->piz_get_line_get_num_subfields    += src->piz_get_line_get_num_subfields;
    dst->piz_get_genotype_sample_starts    += src->piz_get_genotype_sample_starts;
    dst->piz_get_line_subfields            += src->piz_get_line_subfields;
    dst->piz_merge_line                    += src->piz_merge_line;

    dst->zip_compress_variant_block        += src->zip_compress_variant_block;
    dst->seg_all_data_lines          += src->seg_all_data_lines;
    dst->zip_generate_haplotype_sections   += src->zip_generate_haplotype_sections;
    dst->count_alt_alleles                 += src->count_alt_alleles;
    dst->sample_haplotype_data             += src->sample_haplotype_data;
    dst->zip_generate_genotype_sections    += src->zip_generate_genotype_sections;
    dst->zip_generate_phase_sections       += src->zip_generate_phase_sections;
    dst->zip_generate_variant_data_section += src->zip_generate_variant_data_section;

    dst->mtf_merge_in_vb_ctx               += src->mtf_merge_in_vb_ctx;
    dst->mtf_clone_ctx                     += src->mtf_clone_ctx;
    dst->mtf_integrate_dictionary_fragment += src->mtf_integrate_dictionary_fragment;

    dst->tmp1                              += src->tmp1;
    dst->tmp2                              += src->tmp2;
    dst->tmp3                              += src->tmp3;
    dst->tmp4                              += src->tmp4;
    dst->tmp5                              += src->tmp5;
}

static inline unsigned ms(long long ns) { return (unsigned)(ns / 1000000);}

const char *profiler_print_short (const ProfilerRec *p)
{
    static char str[2000]; // not thread safe

    sprintf (str, "read: %u compute:%u (compressor: %u) write: %u", 
             ms(p->read), ms(p->piz_uncompress_variant_block + p->zip_compress_variant_block), ms(p->compressor), ms(p->write));

    return str;
}

const char *profiler_print_report (const ProfilerRec *p, unsigned max_threads, const char *filename)
{
    static char str[2000]; // not thread safe

#if __WIN32__
    static const char *os ="Windows";
#elif defined __APPLE__
    static const char *os ="MacOS";
#else
    static const char *os ="Linux";
#endif

    if (p->piz_uncompress_variant_block)
        sprintf (str, "OS=%s\n"
#ifdef DEBUG
                      "Build=DEBUG\n"
#else
                      "Build=OPTIMIZED\n"
#endif
                      "max_threads=%u\n"
                      "file=%s\n\n"

                      "PIZ dispatcher thread:\n"
                      "   read: %u\n"
                      "   mtf_integrate_dictionary_fragment: %u\n"
                      "   write: %u\n"
                      "piz_uncompress_variant_block: %u\n"
                      "   zfile_uncompress_section: %u\n"
                      "   piz_reconstruct_line_components: %u\n"
                      "      piz_get_variant_data_line: %u\n"
                      "          piz_decode_pos: %u\n"
                      "      piz_get_line_get_num_subfields: %u\n"
                      "      piz_get_line_subfields: %u\n"
                      "      piz_get_haplotype_data_line: %u\n"
                      "      piz_get_genotype_sample_starts: %u\n"
                      "      piz_get_genotype_data_line: %u\n"
                      "      piz_get_phase_data_line: %u\n"
                      "      piz_merge_line: %u\n"
                      "   squeeze: %u\n"
                      "tmp1: %u tmp2: %u tmp3: %u tmp4: %u tmp5: %u\n",
                 os, max_threads, filename ? filename : "(not file)",
                 ms(p->read), ms(p->mtf_integrate_dictionary_fragment), ms(p->write), 
                 ms(p->piz_uncompress_variant_block), ms(p->zfile_uncompress_section),
                 ms(p->piz_reconstruct_line_components), ms(p->piz_get_variant_data_line), ms(p->piz_decode_pos),
                 ms(p->piz_get_line_get_num_subfields), ms(p->piz_get_line_subfields),
                 ms(p->piz_get_haplotype_data_line), 
                 ms(p->piz_get_genotype_sample_starts), ms(p->piz_get_genotype_data_line), 
                 ms(p->piz_get_phase_data_line), ms(p->piz_merge_line),
                 ms(p->squeeze), 
                 ms(p->tmp1), ms(p->tmp2), ms(p->tmp3), ms(p->tmp4), ms(p->tmp5));
                 
    else // compress
        sprintf (str, "OS=%s\n"
#ifdef DEBUG
                      "Build=DEBUG\n"
#else
                      "Build=OPTIMIZED\n"
#endif
                      "max_threads=%u\n"
                      "file=%s\n\n"

                      "ZIP dispatcher thread:\n"
                      "   read: %u\n"
                      "   write: %u\n"
                      "zip_compress_variant_block: %u\n"
                      "   compressor: %u\n"
                      "   seg_all_data_lines: %u\n"
                      "   zip_generate_haplotype_sections: %u\n"
                      "      count_alt_alleles: %u\n"
                      "      sample_haplotype_data: %u\n"
                      "   zip_generate_genotype_sections: %u\n"
                      "   zip_generate_phase_sections: %u\n"
                      "   zip_generate_variant_data_section: %u\n"
                      "   mtf_clone_ctx: %u\n"
                      "   mtf_merge_in_vb_ctx: %u\n"
                      "tmp1: %u tmp2: %u tmp3: %u tmp4: %u tmp5: %u\n",
                 os, max_threads, filename ? filename : "(not file)",
                 ms(p->read), ms(p->write), ms(p->zip_compress_variant_block), ms(p->compressor), 
                 ms(p->seg_all_data_lines),
                 ms(p->zip_generate_haplotype_sections), ms(p->count_alt_alleles), ms(p->sample_haplotype_data), 
                 ms(p->zip_generate_genotype_sections),
                 ms(p->zip_generate_phase_sections), ms(p->zip_generate_variant_data_section),
                 ms(p->mtf_clone_ctx), ms(p->mtf_merge_in_vb_ctx), 
                 ms(p->tmp1), ms(p->tmp2), ms(p->tmp3), ms(p->tmp4), ms(p->tmp5));
    return str;
}
