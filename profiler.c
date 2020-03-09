// ------------------------------------------------------------------
//   profiler.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"

void profiler_add (ProfilerRec *dst, const ProfilerRec *src)
{
    dst->read                              += src->read;
    dst->compute                           += src->compute;
    dst->write                             += src->write;
    dst->compressor                        += src->compressor;
    dst->piz_reconstruct_line_components   += src->piz_reconstruct_line_components;
    dst->piz_get_variant_data_line         += src->piz_get_variant_data_line;
    dst->piz_decode_pos                    += src->piz_decode_pos;
    dst->piz_get_haplotype_data_line       += src->piz_get_haplotype_data_line;
    dst->piz_get_phase_data_line           += src->piz_get_phase_data_line;
    dst->piz_get_genotype_data_line        += src->piz_get_genotype_data_line;
    dst->zfile_uncompress_section          += src->zfile_uncompress_section;
    dst->squeeze                           += src->squeeze;
    dst->buf_alloc                         += src->buf_alloc;
    dst->piz_get_format_info               += src->piz_get_format_info;
    dst->piz_initialize_sample_iterators    += src->piz_initialize_sample_iterators;
    dst->piz_get_line_subfields            += src->piz_get_line_subfields;
    dst->piz_merge_line                    += src->piz_merge_line;
    dst->zfile_read_one_vb                 += src->zfile_read_one_vb;

    dst->vcffile_read_variant_block        += src->vcffile_read_variant_block;
    dst->vcffile_read_vcf_header           += src->vcffile_read_vcf_header;
    dst->seg_all_data_lines                += src->seg_all_data_lines;
    dst->zip_generate_haplotype_sections   += src->zip_generate_haplotype_sections;
    dst->count_alt_alleles                 += src->count_alt_alleles;
    dst->sample_haplotype_data             += src->sample_haplotype_data;
    dst->zip_generate_genotype_sections    += src->zip_generate_genotype_sections;
    dst->zip_generate_phase_sections       += src->zip_generate_phase_sections;
    dst->zip_generate_variant_data_section += src->zip_generate_variant_data_section;
    dst->md5                               += src->md5;
    
    dst->mtf_merge_in_vb_ctx               += src->mtf_merge_in_vb_ctx;
    dst->gl_optimize_dictionary            += src->gl_optimize_dictionary;
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
    static char str[200]; // not thread safe
    sprintf (str, "read: %u compute:%u write: %u", ms(p->read), ms(p->compute), ms(p->write));
    return str;
}

void profiler_print_report (const ProfilerRec *p, unsigned max_threads, unsigned used_threads, const char *filename, unsigned num_vbs)
{
#if defined _WIN32 && defined _WIN64
    static const char *os ="Windows (64 bit)";
#elif defined _WIN32 
    static const char *os ="Windows (32 bit)";
#elif defined __APPLE__
    static const char *os ="MacOS";
#else
    static const char *os ="Linux";
#endif

    fprintf (stderr, "\n%s PROFILER:\n", command == ZIP ? "ZIP" : "UNZIP");
    fprintf (stderr, "OS=%s\n", os);
#ifdef DEBUG
    fprintf (stderr, "Build=Debug\n");
#else
    fprintf (stderr, "Build=Optimized\n");
#endif
    fprintf (stderr, "Compute threads: max_permitted=%u actually_used=%u\n", max_threads, used_threads);
    fprintf (stderr, "file=%s\n\n", filename ? filename : "(not file)");
    
    fprintf (stderr, "Wallclock: %u milliseconds\n", ms (p->wallclock));

    if (command != ZIP) { // this is a uncompress operation

        fprintf (stderr, "GENOUNZIP I/O thread (piz_dispatcher):\n");
        fprintf (stderr, "   zfile_read_one_vb: %u\n", ms(p->zfile_read_one_vb));
        fprintf (stderr, "      read: %u\n", ms(p->read));
        fprintf (stderr, "      mtf_integrate_dictionary_fragment: %u\n", ms(p->mtf_integrate_dictionary_fragment));
        fprintf (stderr, "   write: %u\n", ms(p->write));
        fprintf (stderr, "GENOUNZIP compute threads (piz_uncompress_variant_block): %u\n", ms(p->compute));
        fprintf (stderr, "   zfile_uncompress_section: %u\n", ms(p->zfile_uncompress_section));
        fprintf (stderr, "   piz_reconstruct_line_components: %u\n", ms(p->piz_reconstruct_line_components));
        fprintf (stderr, "      piz_get_variant_data_line: %u\n", ms(p->piz_get_variant_data_line));
        fprintf (stderr, "          piz_decode_pos: %u\n", ms(p->piz_decode_pos) );
        fprintf (stderr, "      piz_get_format_info: %u\n", ms(p->piz_get_format_info));
        fprintf (stderr, "      piz_get_line_subfields: %u\n", ms(p->piz_get_line_subfields));
        fprintf (stderr, "      piz_get_haplotype_data_line: %u\n", ms(p->piz_get_haplotype_data_line));
        fprintf (stderr, "      piz_initialize_sample_iterators: %u\n", ms(p->piz_initialize_sample_iterators));
        fprintf (stderr, "      piz_get_genotype_data_line: %u\n", ms(p->piz_get_genotype_data_line));
        fprintf (stderr, "      piz_get_phase_data_line: %u\n", ms(p->piz_get_phase_data_line));
        fprintf (stderr, "      piz_merge_line: %u\n", ms(p->piz_merge_line));
        fprintf (stderr, "   squeeze: %u\n", ms(p->squeeze));
    }
    else { // compress
        fprintf (stderr, "GENOZIP I/O thread (zip_dispatcher):\n");
        fprintf (stderr, "   vcffile_read_vcf_header: %u\n", ms(p->vcffile_read_vcf_header));
        fprintf (stderr, "   vcffile_read_variant_block: %u\n", ms(p->vcffile_read_variant_block));
        fprintf (stderr, "      read: %u\n", ms(p->read));
        fprintf (stderr, "      md5: %u\n", ms(p->md5));
        fprintf (stderr, "   write: %u\n", ms(p->write));
        fprintf (stderr, "GENOZIP compute threads (zip_compress_one_vb): %u\n", ms(p->compute));
        fprintf (stderr, "   compressor: %u\n", ms(p->compressor));
        fprintf (stderr, "   seg_all_data_lines: %u\n", ms(p->seg_all_data_lines));
        fprintf (stderr, "   zip_generate_haplotype_sections: %u\n", ms(p->zip_generate_haplotype_sections));
        fprintf (stderr, "      count_alt_alleles: %u\n", ms(p->count_alt_alleles));
        fprintf (stderr, "      sample_haplotype_data: %u\n", ms(p->sample_haplotype_data));
        fprintf (stderr, "   zip_generate_genotype_sections: %u\n", ms(p->zip_generate_genotype_sections));
        fprintf (stderr, "   zip_generate_phase_sections: %u\n", ms(p->zip_generate_phase_sections));
        fprintf (stderr, "   zip_generate_variant_data_section: %u\n", ms(p->zip_generate_variant_data_section));
        fprintf (stderr, "   mtf_clone_ctx: %u\n", ms(p->mtf_clone_ctx));
        fprintf (stderr, "   mtf_merge_in_vb_ctx: %u\n", ms(p->mtf_merge_in_vb_ctx));
        fprintf (stderr, "      gl_optimize_dictionary: %u\n", ms(p->gl_optimize_dictionary));
    }    
    fprintf (stderr, "buf_alloc: %u\n", ms(p->buf_alloc));
    fprintf (stderr, "tmp1: %u tmp2: %u tmp3: %u tmp4: %u tmp5: %u\n\n", ms(p->tmp1), ms(p->tmp2), ms(p->tmp3), ms(p->tmp4), ms(p->tmp5));

    fprintf (stderr, "\nVariant Block stats:\n");
    fprintf (stderr, "  Variant blocks: %u\n", num_vbs);
    fprintf (stderr, "  Maximum variant block size: %u MB\n", global_max_memory_per_vb / (1024 * 1024));
    fprintf (stderr, "  Average wallclock: %u\n", ms(p->wallclock) / num_vbs);
    fprintf (stderr, "  Average read time: %u\n", ms(p->read) / num_vbs);
    fprintf (stderr, "  Average compute time: %u\n", ms(p->compute) / num_vbs);
    fprintf (stderr, "  Average write time: %u\n", ms(p->write) / num_vbs);
    fprintf (stderr, "\n");
}
