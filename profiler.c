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

static ProfilerRec profile = {};    // data for entire execution 
static Mutex profile_mutex = {};
static TimeSpecType profiler_timer; // wallclock

void profiler_initialize (void)
{
    mutex_initialize (profile_mutex);
    clock_gettime(CLOCK_REALTIME, &profiler_timer); // initialze wallclock
}

void profiler_add (ConstVBlockP vb)
{
#   define ADD(x) profile.x += vb->profile.x
    
    mutex_lock (profile_mutex);

    if (vb->txt_data.len) {
        profile.num_vbs++;
        profile.max_vb_size_mb = MAX_(profile.max_vb_size_mb, segconf.vb_size >> 20);
    }

    ADD(file_open);
    ADD(file_close);
    ADD(buf_low_level_free);
    ADD(buf_remove_from_buffer_list);
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
    ADD(compressor_longr);
    ADD(compressor_rans);
    ADD(compressor_arith);
    ADD(zip_generate_b250);
    ADD(zip_generate_local);
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
    ADD(qname_seg);
    ADD(ctx_merge_in_vb_ctx);
    ADD(codec_hapmat_count_alt_alleles);
    ADD(md5);
    ADD(wait_for_vb_1_mutex);
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

    mutex_unlock (profile_mutex);

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
// this prints the compressor time for each context. temporarily disabled bc we moved profile to be global rather than file by file
// TO DO: revive this, when compressing a single file
/*    for (DidIType did_i=0; did_i < z_file->num_contexts; did_i++) {
        ContextP ctx = ZCTX(did_i);
        if (ms(ctx->compressor_time))
            iprintf ("      %s: %u\n", ctx->tag_name, ms(ctx->compressor_time));
    }
*/
}

void profiler_print_report (void)
{
    static const char *space = "                                                   ";
#   define PRINT(x, level) if (ms(profile.x)) iprintf ("%.*s" #x ": %u\n", level*3, space, ms(profile.x));
    
    const char *os = flag.is_windows ? "Windows"
                   : flag.is_mac     ? "MacOS"
                   : flag.is_linux   ? "Linux"
                   :                   "Unknown OS";

    iprintf ("\n%s PROFILER:\n", command == ZIP ? "ZIP" : "PIZ");
    iprintf ("OS=%s\n", os);
    iprintf ("Build=%s\n", flag.debug ? "Debug" : "Optimized");

    iprintf ("Wallclock: %s milliseconds\n", str_uint_commas (ms (CHECK_TIMER)).s);

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
        iprintf ("GENOUNZIP compute threads: %u\n", ms(profile.compute));
        PRINT (zfile_uncompress_section, 1);
        PRINT (compressor_bz2,  2);
        PRINT (compressor_lzma, 2);
        PRINT (compressor_bsc,  2);
        PRINT (compressor_rans, 2);
        PRINT (compressor_arith, 2);
        PRINT (compressor_domq, 2);
        PRINT (compressor_actg, 2);
        PRINT (compressor_pbwt, 2);
        PRINT (compressor_longr, 2);
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
        iprintf ("GENOZIP compute threads %u\n", ms(profile.compute));
        PRINT (ctx_clone, 1);
        PRINT (seg_all_data_lines, 1);
        PRINT (sam_seg_XA_pos, 2);
        PRINT (aligner_best_match, 2);
        PRINT (aligner_get_match_len, 3);
        PRINT (aligner_get_word_from_seq, 3);
        PRINT (seg_initialize, 2);
        PRINT (qname_seg,2);
        PRINT (sam_cigar_seg,2);
        PRINT (sam_seg_SEQ,2);
        PRINT (wait_for_vb_1_mutex, 1);
        PRINT (ctx_merge_in_vb_ctx, 1);
        PRINT (zip_compress_ctxs, 1);
        PRINT (zip_generate_b250, 2);
        PRINT (zip_generate_local, 2);
        print_ctx_compressor_times();

        PRINT (ctx_compress_one_dict_fragment, 1);
        PRINT (codec_assign_best_codec, 1);
        PRINT (md5, 1);
        PRINT (bgzf_compute_thread, 1);
        PRINT (compressor_bz2,  1);
        PRINT (compressor_lzma, 1);
        PRINT (compressor_bsc,  1);
        PRINT (compressor_rans, 1);
        PRINT (compressor_arith, 1);
        PRINT (compressor_domq, 1);
        PRINT (compressor_actg, 1);
        PRINT (compressor_pbwt, 1);
        PRINT (compressor_longr, 1);
        PRINT (codec_hapmat_count_alt_alleles, 2);
    }    

    PRINT (file_open, 0);
    PRINT (file_close, 0);
    PRINT (buf_alloc, 0);
    PRINT (buf_remove_from_buffer_list, 0);
    PRINT (buf_low_level_free, 0);
    PRINT (buf_mmap_do, 0);
    PRINT (vb_release_vb_do, 0);
    PRINT (vb_destroy_vb, 0);
    PRINT (dispatcher_recycle_vbs, 0);
    PRINT (generate_rev_complement_genome, 0);
    
    iprintf ("tmp1: %u tmp2: %u tmp3: %u tmp4: %u tmp5: %u\n\n", ms(profile.tmp1), ms(profile.tmp2), ms(profile.tmp3), ms(profile.tmp4), ms(profile.tmp5));

    if (profile.num_vbs) {
        iprint0 ("\nVblock stats:\n");
        iprintf ("  Vblocks: %u\n", profile.num_vbs);
        iprintf ("  Maximum vblock size: %u MB\n", profile.max_vb_size_mb);
        iprintf ("  Average wallclock: %u\n", (unsigned)(ms(CHECK_TIMER) / profile.num_vbs));
        iprintf ("  Average read time: %u\n", ms(profile.read) / profile.num_vbs);
        iprintf ("  Average compute time: %u\n", ms(profile.compute) / profile.num_vbs);
        iprintf ("  Average write time: %u\n", ms(profile.write) / profile.num_vbs);
    }
    
    iprint0 ("\n");
}
