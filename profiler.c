// ------------------------------------------------------------------
//   profiler.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"

void profiler_add (ProfilerRec *dst, const ProfilerRec *src)
{
#   define ADD(x) dst->x += src->x

    ADD(read);
    ADD(compute);
    ADD(write);
    ADD(compressor_bz2);
    ADD(compressor_lzma);
    ADD(piz_reconstruct_vb);
    ADD(vcf_piz_get_variant_data_line);
    ADD(vcf_piz_get_haplotype_data_line);
    ADD(vcf_piz_get_phase_data_line);
    ADD(vcf_piz_reconstruct_genotype_data_line);
    ADD(zfile_uncompress_section);
    ADD(buf_alloc);
    ADD(vcf_piz_initialize_sample_iterators);
    ADD(piz_get_line_subfields);
    ADD(vcf_piz_reconstruct_samples);
    ADD(piz_read_one_vb);
    ADD(zfile_compress_dictionary_data);
    ADD(txtfile_read_vblock);
    ADD(txtfile_read_header);
    ADD(seg_all_data_lines);
    ADD(vcf_zip_generate_haplotype_sections);
    ADD(count_alt_alleles);
    ADD(sample_haplotype_data);
    ADD(zip_generate_genotype_sections);
    ADD(vcf_zip_generate_phase_sections);
    ADD(zip_generate_variant_data_section);
    ADD(md5);
    ADD(lock_mutex_compress_dict);
    ADD(lock_mutex_zf_ctx);
    ADD(mtf_merge_in_vb_ctx_one_dict_id);
    ADD(mtf_clone_ctx);
    ADD(mtf_integrate_dictionary_fragment);
    ADD(refhash_best_match);
    ADD(refhash_get_match_len);
    ADD(refhash_get_word_from_seq);
    ADD(generate_rev_complement_genome);
    ADD(tmp1);
    ADD(tmp2);
    ADD(tmp3);
    ADD(tmp4);
    ADD(tmp5);

#undef ADD
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
    static const char *space = "                                                   ";
#   define PRINT(x, level) if (p->x) fprintf (stderr, "%.*s" #x ": %u\n", level*3, space, ms(p->x));
    
#if defined _WIN32
    static const char *os ="Windows";
#elif defined __APPLE__
    static const char *os ="MacOS";
#elif defined __linux__
    static const char *os ="Linux";
#else
    static const char *os ="Unknown OS";
#endif

    fprintf (stderr, "\n%s PROFILER:\n", command == ZIP ? "ZIP" : "PIZ");
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
        PRINT (piz_read_one_vb, 1);
        PRINT (read, 2);
        PRINT (mtf_integrate_dictionary_fragment, 2);
        PRINT (write, 1);
        fprintf (stderr, "GENOUNZIP compute threads (vcf_piz_uncompress_vb): %u\n", ms(p->compute));
        PRINT (zfile_uncompress_section, 1);
        PRINT (piz_reconstruct_vb, 1);
        PRINT (vcf_piz_get_variant_data_line, 2);
        PRINT (piz_get_line_subfields, 2);
        PRINT (vcf_piz_get_haplotype_data_line, 2);
        PRINT (vcf_piz_initialize_sample_iterators, 2);
        PRINT (vcf_piz_reconstruct_genotype_data_line, 2);
        PRINT (vcf_piz_get_phase_data_line, 2);
        PRINT (vcf_piz_reconstruct_samples, 2);
    }
    else { // compress
        fprintf (stderr, "GENOZIP I/O thread (zip_dispatcher):\n");
        PRINT (txtfile_read_header, 1);
        PRINT (txtfile_read_vblock, 1);
        PRINT (read, 2);
        PRINT (md5, 2);
        PRINT (write, 1);
        fprintf (stderr, "GENOZIP compute threads (vcf_zip_compress_one_vb): %u\n", ms(p->compute));
        PRINT (compressor_bz2, 1);
        PRINT (compressor_lzma, 1);
        PRINT (seg_all_data_lines, 1);
        PRINT (refhash_best_match, 2);
        PRINT (refhash_get_match_len, 3);
        PRINT (refhash_get_word_from_seq, 3);
        PRINT (vcf_zip_generate_haplotype_sections, 1);
        PRINT (count_alt_alleles, 2);
        PRINT (sample_haplotype_data, 2);
        PRINT (zip_generate_genotype_sections, 1);
        PRINT (vcf_zip_generate_phase_sections, 1);
        PRINT (zip_generate_variant_data_section, 1);
        PRINT (mtf_clone_ctx, 1);
        PRINT (lock_mutex_zf_ctx, 1);
        PRINT (mtf_merge_in_vb_ctx_one_dict_id, 2);
        PRINT (lock_mutex_compress_dict, 2);
        PRINT (zfile_compress_dictionary_data, 3);
    }    

    PRINT (buf_alloc, 0);
    PRINT (generate_rev_complement_genome, 0);
    
    fprintf (stderr, "tmp1: %u tmp2: %u tmp3: %u tmp4: %u tmp5: %u\n\n", ms(p->tmp1), ms(p->tmp2), ms(p->tmp3), ms(p->tmp4), ms(p->tmp5));

    fprintf (stderr, "\nVblock stats:\n");
    fprintf (stderr, "  Vblocks: %u\n", num_vbs);
    fprintf (stderr, "  Maximum vblock size: %u MB\n", global_max_memory_per_vb / (1024 * 1024));
    fprintf (stderr, "  Average wallclock: %u\n", ms(p->wallclock) / num_vbs);
    fprintf (stderr, "  Average read time: %u\n", ms(p->read) / num_vbs);
    fprintf (stderr, "  Average compute time: %u\n", ms(p->compute) / num_vbs);
    fprintf (stderr, "  Average write time: %u\n", ms(p->write) / num_vbs);
    fprintf (stderr, "\n");
}
