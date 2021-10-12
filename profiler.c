// ------------------------------------------------------------------
//   profiler.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "profiler.h"
#include "flags.h"
#include "vblock.h"
#include "file.h"
#include "segconf.h"

static Mutex evb_profile_mutex = {};

void profiler_initialize (void)
{
    mutex_initialize (evb_profile_mutex);
}

void profiler_add (ConstVBlockP vb)
{
#   define ADD(x) evb->profile.x += vb->profile.x
    
    mutex_lock (evb_profile_mutex);

    ADD(read);
    ADD(compute);
    ADD(write);
    ADD(vb_get_vb);
    ADD(compressor_bz2);
    ADD(compressor_lzma);
    ADD(compressor_bsc);
    ADD(compressor_domq);
    ADD(compressor_actg);
    ADD(compressor_pbwt);
    ADD(compressor_hapmat);
    ADD(zip_generate_ctxs);
    ADD(zip_compress_ctxs);
    ADD(codec_assign_best_codec);
    ADD(bgzf_io_thread);
    ADD(bgzf_compute_thread);
    ADD(reconstruct_vb);
    ADD(piz_get_line_subfields);
    ADD(piz_read_one_vb);
    ADD(codec_hapmat_piz_get_one_line);
    ADD(sam_seg_SEQ);
    ADD(sam_cigar_seg);
    ADD(ctx_compress_one_dict_fragment);
    ADD(zfile_uncompress_section);
    ADD(buf_alloc);
    ADD(dispatcher_recycle_vbs);
    ADD(txtfile_read_vblock);
    ADD(txtfile_read_header);
    ADD(seg_all_data_lines);
    ADD(seg_initialize);
    ADD(sam_seg_XA_pos);
    ADD(compound_seg);
    ADD(ctx_merge_in_vb_ctx);
    ADD(codec_hapmat_count_alt_alleles);
    ADD(md5);
    ADD(lock_mutex_zf_ctx);
    ADD(ctx_merge_in_vb_ctx_one_dict_id);
    ADD(ctx_clone);
    ADD(ctx_dict_build_word_lists);
    ADD(aligner_best_match);
    ADD(aligner_get_match_len);
    ADD(aligner_get_word_from_seq);
    ADD(generate_rev_complement_genome);
    ADD(ctx_read_all_dictionaries);
    ADD(ref_contigs_compress);
    ADD(linesorter_compress_recon_plan);
    ADD(linesorter_compress_qsort);
    ADD(ref_load_stored_reference);
    ADD(ref_read_one_range);
    ADD(ref_uncompress_one_range);
    ADD(tmp1);
    ADD(tmp2);
    ADD(tmp3);
    ADD(tmp4);
    ADD(tmp5);

    mutex_unlock (evb_profile_mutex);

#undef ADD
}

static inline uint32_t ms(uint64_t ns) { return (uint32_t)(ns / 1000000);}

const char *profiler_print_short (const ProfilerRec *p)
{
    static char str[200]; // not thread safe
    sprintf (str, "read: %u compute:%u write: %u", ms(p->read), ms(p->compute), ms(p->write));
    return str;
}

static void print_ctx_compressor_times (void)
{
    for (DidIType did_i=0; did_i < z_file->num_contexts; did_i++) {
        ContextP ctx = ZCTX(did_i);
        if (ctx->compressor_time)
            iprintf ("      %s: %u\n", ctx->tag_name, ms(ctx->compressor_time));
    }
}

void profiler_print_report (const ProfilerRec *p, unsigned max_threads, unsigned used_threads, const char *filename, unsigned num_vbs)
{
    static const char *space = "                                                   ";
#   define PRINT(x, level) if (p->x) iprintf ("%.*s" #x ": %u\n", level*3, space, ms(p->x));
    
    const char *os = flag.is_windows ? "Windows"
                   : flag.is_mac     ? "MacOS"
                   : flag.is_linux   ? "Linux"
                   :                   "Unknown OS";

    iprintf ("\n%s PROFILER:\n", command == ZIP ? "ZIP" : "PIZ");
    iprintf ("OS=%s\n", os);
    iprintf ("Build=%s\n", flag.debug ? "Debug" : "Optimized");
    iprintf ("Compute threads: max_permitted=%u actually_used=%u\n", max_threads, used_threads);
    iprintf ("file=%s\n\n", filename ? filename : "(not file)");
    
    iprintf ("Wallclock: %u milliseconds\n", ms (p->wallclock));

    if (command != ZIP) { // this is a uncompress operation

        iprint0 ("GENOUNZIP main thread (piz_one_txt_file):\n");
        PRINT (piz_read_global_area, 1);
        PRINT (ref_load_stored_reference, 2);
        PRINT (ref_read_one_range, 3);
        PRINT (ref_uncompress_one_range, 3);
        PRINT (ctx_read_all_dictionaries, 2);
        PRINT (ctx_dict_build_word_lists, 3);
        PRINT (vb_get_vb, 1);
        PRINT (piz_read_one_vb, 1);
        PRINT (read, 2);
        PRINT (bgzf_io_thread, 1);
        PRINT (write, 1);
        iprintf ("GENOUNZIP compute threads: %u\n", ms(p->compute));
        PRINT (zfile_uncompress_section, 1);
        PRINT (compressor_bz2,  2);
        PRINT (compressor_lzma, 2);
        PRINT (compressor_bsc,  2);
        PRINT (compressor_domq, 2);
        PRINT (compressor_actg, 2);
        PRINT (compressor_hapmat, 2);
        PRINT (compressor_pbwt, 2);
        print_ctx_compressor_times();
        PRINT (reconstruct_vb, 1);
        PRINT (md5, 1);
        PRINT (bgzf_compute_thread, 1);
        PRINT (piz_get_line_subfields, 2);
        PRINT (codec_hapmat_piz_get_one_line, 2);
    }
    else { // compress
        iprint0 ("GENOZIP main thread (zip_one_file):\n");
        PRINT (txtfile_read_header, 1);
        PRINT (vb_get_vb, 1);
        PRINT (txtfile_read_vblock, 1);
        PRINT (read, 2);
        PRINT (write, 1);
        PRINT (bgzf_io_thread, 1);
        PRINT (ref_contigs_compress, 1);
        PRINT (linesorter_compress_recon_plan, 1);
        PRINT (linesorter_compress_qsort, 2);
        iprintf ("GENOZIP compute threads %u\n", ms(p->compute));
        PRINT (ctx_clone, 1);
        PRINT (seg_all_data_lines, 1);
        PRINT (sam_seg_XA_pos, 2);
        PRINT (aligner_best_match, 2);
        PRINT (aligner_get_match_len, 3);
        PRINT (aligner_get_word_from_seq, 3);
        PRINT (seg_initialize, 2);
        PRINT (compound_seg,2);
        PRINT (sam_cigar_seg,2);
        PRINT (sam_seg_SEQ,2);
        PRINT (ctx_merge_in_vb_ctx, 1);
        PRINT (lock_mutex_zf_ctx, 2);
        PRINT (ctx_merge_in_vb_ctx_one_dict_id, 2);
        PRINT (zip_generate_ctxs, 1);
        PRINT (zip_compress_ctxs, 1);
        print_ctx_compressor_times();

        PRINT (ctx_compress_one_dict_fragment, 1);
        PRINT (codec_assign_best_codec, 1);
        PRINT (md5, 1);
        PRINT (bgzf_compute_thread, 1);
        PRINT (compressor_bz2,  1);
        PRINT (compressor_lzma, 1);
        PRINT (compressor_bsc,  1);
        PRINT (compressor_domq, 1);
        PRINT (compressor_actg, 1);
        PRINT (compressor_hapmat, 1);
        PRINT (compressor_pbwt, 1);
        PRINT (codec_hapmat_count_alt_alleles, 2);
    }    

    PRINT (buf_alloc, 0);
    PRINT (buf_mmap_do, 0);
    PRINT (vb_release_vb_do, 0);
    PRINT (vb_destroy_vb, 0);
    PRINT (dispatcher_recycle_vbs, 0);
    PRINT (generate_rev_complement_genome, 0);
    
    iprintf ("tmp1: %u tmp2: %u tmp3: %u tmp4: %u tmp5: %u\n\n", ms(p->tmp1), ms(p->tmp2), ms(p->tmp3), ms(p->tmp4), ms(p->tmp5));

    iprint0 ("\nVblock stats:\n");
    iprintf ("  Vblocks: %u\n", num_vbs);
    iprintf ("  Maximum vblock size: %"PRIu64" MB\n", segconf.vb_size >> 20);
    iprintf ("  Average wallclock: %u\n", ms(p->wallclock) / num_vbs);
    iprintf ("  Average read time: %u\n", ms(p->read) / num_vbs);
    iprintf ("  Average compute time: %u\n", ms(p->compute) / num_vbs);
    iprintf ("  Average write time: %u\n", ms(p->write) / num_vbs);
    iprint0 ("\n");
}
