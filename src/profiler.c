// ------------------------------------------------------------------
//   profiler.c
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "profiler.h"
#include "flags.h"
#include "vblock.h"
#include "file.h"
#include "segconf.h"
#include "mutex.h"

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

    ADD (file_open);
    ADD (file_close);
    ADD (buf_low_level_free);
    ADD (buf_remove_from_buffer_list);
    ADD (read);
    ADD (compute);
    ADD (write);
    ADD (vb_get_vb);
    ADD (compressor_bz2);
    ADD (compressor_lzma);
    ADD (compressor_bsc);
    ADD (compressor_domq);
    ADD (compressor_actg);
    ADD (compressor_pbwt);
    ADD (compressor_longr);
    ADD (compressor_rans);
    ADD (compressor_arith);
    ADD (compressor_normq);
    ADD (zip_generate_b250);
    ADD (zip_generate_local);
    ADD (zip_compress_ctxs);
    ADD (codec_assign_best_codec);
    ADD (bgzf_io_thread);
    ADD (bgzf_writer_thread);
    ADD (bgzf_compute_thread);
    ADD (reconstruct_vb);
    ADD (piz_get_line_subfields);
    ADD (piz_read_one_vb);
    ADD (codec_hapmat_piz_get_one_line);
    ADD (codec_longr_reconstruct);
    ADD (codec_domq_reconstruct);
    ADD (codec_domq_reconstruct_dom_run);    
    ADD (sam_seg_SEQ);
    ADD (sam_seg_QUAL);
    ADD (sam_seg_is_gc_line);
    ADD (sam_seg_aux_all);
    ADD (sam_cigar_binary_to_textual);
    ADD (squank_seg);
    ADD (sam_seg_MD_Z_analyze);
    ADD (bam_seq_to_sam);
    ADD (sam_seg_bsseeker2_XG_Z_analyze);
    ADD (sam_seg_sag_stuff);
    ADD (sam_cigar_seg);
    ADD (sam_seg_bismark_XM_Z);
    ADD (sam_seg_bsbolt_XB);
    ADD (sam_piz_sam2fastq_QUAL); 
    ADD (sam_piz_sam2bam_QUAL);        
    ADD (sam_cigar_special_CIGAR);
    ADD (sam_bismark_piz_update_meth_call);
    ADD (sam_piz_special_QUAL);
    ADD (sam_seg_SEQ_vs_ref);
    ADD (sam_seg_bisulfite_M);
    ADD (sam_reconstruct_SEQ_vs_ref);
    ADD (sam_header_add_contig); 
    ADD (contigs_create_index);
    ADD (sam_header_zip_inspect_PG_lines); 
    ADD (sam_header_zip_inspect_HD_line);
    ADD (ref_initialize_ranges);
    ADD (txtheader_zip_read_and_compress);
    ADD (txtheader_compress);
    ADD (txtheader_compress_one_fragment);
    ADD (txtheader_piz_read_and_reconstruct);
    ADD (digest_txt_header);
    ADD (scan_index_qnames_preprocessing);
    ADD (aligner_seg_seq);
    ADD (aligner_reconstruct_seq);    
    ADD (dict_io_compress_one_fragment);
    ADD (dict_io_compress_dictionaries); 
    ADD (dict_io_assign_codecs); 
    ADD (recon_plan_compress);
    ADD (recon_plan_compress_one_fragment);
    ADD (zfile_uncompress_section);
    ADD (buf_alloc);
    ADD (dispatcher_recycle_vbs);
    ADD (txtfile_read_vblock);
    ADD (txtfile_read_header);
    ADD (seg_all_data_lines);
    ADD (seg_initialize);
    ADD (sam_seg_BWA_XA_pos);
    ADD (sam_seg_BWA_XA_Z);
    ADD (sam_seg_BWA_XS_i);
    ADD (sam_seg_AS_i);
    ADD (sam_seg_NM_i);
    ADD (sam_seg_SA_Z);
    ADD (sam_seg_TX_AN_Z);
    ADD (sam_seg_barcode_qual);
    ADD (sam_seg_CB_Z);
    ADD (sam_seg_CR_Z);
    ADD (sam_seg_RX_Z); 
    ADD (sam_seg_BX_Z);
    ADD (sam_seg_QX_Z);
    ADD (sam_seg_BC_Z);
    ADD (sam_seg_gene_name_id);
    ADD (sam_seg_fx_Z);
    ADD (sam_seg_other_seq);
    ADD (sam_seg_GR_Z);
    ADD (sam_seg_GY_Z);
    ADD (qname_seg);
    ADD (ctx_merge_in_vb_ctx);
    ADD (codec_hapmat_count_alt_alleles);
    ADD (digest);
    ADD (ctx_clone);
    ADD (dict_io_build_word_lists);
    ADD (aligner_best_match);
    ADD (aligner_get_match_len);
    ADD (aligner_get_word_from_seq);
    ADD (sam_seg_verify_saggy_line_SEQ); 
    ADD (reconstruct_SEQ_copy_sag_prim);
    ADD (sam_analyze_copied_SEQ);
    ADD (generate_rev_complement_genome);
    ADD (dict_io_read_all_dictionaries);
    ADD (ref_contigs_compress);
    ADD (generate_recon_plan);
    ADD (vcf_linesort_compress_qsort);
    ADD (vcf_linesort_merge_vb);
    ADD (vcf_seg_PROBE_A);
    ADD (ref_load_stored_reference);
    ADD (ref_read_one_range);
    ADD (ref_uncompress_one_range);
    ADD (ref_compress_ref);
    ADD (ref_compress_one_range);
    ADD (ref_copy_compressed_sections_from_reference_file);
    ADD (sam_sa_prim_finalize_ingest);
    ADD (sam_zip_prim_ingest_vb);
    ADD (sam_zip_recon_plan_add_gc_lines);
    ADD (sam_zip_gc_calc_depn_vb_info);
    ADD (sam_load_groups_add_one_prim_vb);
    ADD (zip_handle_unique_words_ctxs);
    ADD (ctx_sort_dictionaries_vb_1);
    ADD (random_access_merge_in_vb);
    ADD (gencomp_absorb_vb_gencomp_lines);
    ADD (zip_handle_unique_words_ctxs);
    ADD (random_access_finalize_entries);
    ADD (random_access_compress);
    ADD (ctx_compress_counts);
    ADD (zfile_compress_genozip_header);
    ADD (tmp1);
    ADD (tmp2);
    ADD (tmp3);
    ADD (tmp4);
    ADD (tmp5);

    mutex_unlock (profile_mutex);

#undef ADD
}

static inline uint32_t ms(uint64_t ns) { return (uint32_t)(ns / 1000000);}

rom profiler_print_short (const ProfilerRec *p)
{
    static char str[300]; // not thread safe
    sprintf (str, "read: %s compute:%s write: %s", str_int_commas (ms(p->read)).s, str_int_commas (ms(p->compute)).s, str_int_commas (ms(p->write)).s);
    return str;
}

static void print_ctx_compressor_times (void)
{
// this prints the compressor time for each context. temporarily disabled bc we moved profile to be global rather than file by file
// TO DO: revive this, when compressing a single file
/*    for (Did did_i=0; did_i < z_file->num_contexts; did_i++) {
        ContextP ctx = ZCTX(did_i);
        if (ms(ctx->compressor_time))
            iprintf ("      %s: %u\n", ctx->tag_name, ms(ctx->compressor_time));
    }
*/
}

void profiler_print_report (void)
{
    static rom space = "                                                   ";
#   define PRINT(x, level) if (ms(profile.x)) iprintf ("%.*s" #x ": %s\n", level*3, space, str_int_commas (ms(profile.x)).s);
    
    rom os = flag.is_windows ? "Windows"
           : flag.is_mac     ? "MacOS"
           : flag.is_linux   ? "Linux"
           :                   "Unknown OS";

    iprintf ("\n%s PROFILER:\n", IS_ZIP ? "ZIP" : "PIZ");
    iprintf ("OS=%s\n", os);
    iprintf ("Build=%s\n", flag.debug ? "Debug" : "Optimized");

    iprintf ("Wallclock: %s milliseconds\n", str_int_commas (ms (CHECK_TIMER)).s);

    if (command != ZIP) { // this is a uncompress operation

        iprint0 ("GENOUNZIP main thread (piz_one_txt_file):\n");
        PRINT (piz_read_global_area, 1);
        PRINT (ref_load_stored_reference, 2);
        PRINT (ref_initialize_ranges, 3);
        PRINT (ref_read_one_range, 3);
        PRINT (ref_uncompress_one_range, 3);
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
        iprintf ("GENOUNZIP compute threads: %s\n", str_int_commas (ms(profile.compute)).s);
        PRINT (zfile_uncompress_section, 1);
        PRINT (compressor_bz2,  2);
        PRINT (compressor_lzma, 2);
        PRINT (compressor_bsc,  2);
        PRINT (compressor_rans, 2);
        PRINT (compressor_arith, 2);
        PRINT (compressor_domq, 2);
        PRINT (compressor_normq, 2);
        PRINT (compressor_actg, 2);
        PRINT (compressor_pbwt, 2);
        PRINT (compressor_longr, 2);
        print_ctx_compressor_times();
        PRINT (reconstruct_vb, 1);
        PRINT (sam_reconstruct_SEQ_vs_ref, 2);
        PRINT (reconstruct_SEQ_copy_sag_prim, 3);
        PRINT (sam_analyze_copied_SEQ, 3);
        PRINT (sam_bismark_piz_update_meth_call, 3);
        PRINT (aligner_reconstruct_seq, 2);
        PRINT (sam_piz_special_QUAL, 2);
        if (last_z_dt == DT_SAM || last_z_dt == DT_BAM) {
            PRINT (codec_longr_reconstruct, 3);
            PRINT (codec_domq_reconstruct, 3);
            PRINT (codec_domq_reconstruct_dom_run, 4);
        }
        PRINT (sam_piz_sam2fastq_QUAL, 2); 
        PRINT (sam_piz_sam2bam_QUAL, 2);        
        PRINT (sam_cigar_special_CIGAR, 2);
        PRINT (sam_zip_prim_ingest_vb, 1);
        PRINT (digest, 1); // note: in SAM/BAM digest is done in the writer thread, otherwise its done in the compute thread. TODO: change level to 0 in case of SAM/BAM
        PRINT (bgzf_compute_thread, 1);
        PRINT (piz_get_line_subfields, 2);
        PRINT (codec_hapmat_piz_get_one_line, 2);
        PRINT (sam_load_groups_add_one_prim_vb, 1);
        if (last_z_dt == DT_FASTQ) {
            PRINT (codec_longr_reconstruct, 3);
            PRINT (codec_domq_reconstruct, 2);
            PRINT (codec_domq_reconstruct_dom_run, 3);
        }
        
        if (profile.bgzf_compute_thread)
            iprintf ("GENOUNZIP bgzf threads: %s\n", str_int_commas (ms(profile.bgzf_compute_thread)).s);
    }
    else { // compress
        iprint0 ("GENOZIP main thread (zip_one_file):\n");
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
        PRINT (txtfile_read_vblock, 1);
        PRINT (read, 2);
        PRINT (write, 1);
        PRINT (bgzf_io_thread, 1);
        PRINT (ref_contigs_compress, 1);
        PRINT (generate_recon_plan, 1);
        PRINT (vcf_linesort_compress_qsort, 2);
        PRINT (sam_zip_gc_calc_depn_vb_info, 1);
        PRINT (sam_zip_recon_plan_add_gc_lines, 1);
        PRINT (sam_sa_prim_finalize_ingest, 1);
        iprintf ("GENOZIP compute threads %s\n", str_int_commas (ms(profile.compute)).s);
        PRINT (ctx_clone, 1);
        PRINT (scan_index_qnames_preprocessing, 1);
        PRINT (seg_all_data_lines, 1);
        PRINT (seg_initialize, 2);
        PRINT (qname_seg, 2);
        PRINT (sam_cigar_seg, 2);
        PRINT (squank_seg, 3);
        PRINT (sam_seg_SEQ, 2);
        PRINT (sam_seg_SEQ_vs_ref, 3);
        PRINT (sam_seg_bisulfite_M, 4);
        PRINT (sam_seg_verify_saggy_line_SEQ, 3); 
        PRINT (sam_analyze_copied_SEQ, 3);
        PRINT (aligner_seg_seq, last_z_dt == DT_FASTQ ? 2 : 3);
        PRINT (aligner_best_match, last_z_dt == DT_FASTQ ? 3 : 4);
        PRINT (aligner_get_match_len, last_z_dt == DT_FASTQ ? 3 : 4);
        PRINT (aligner_get_word_from_seq, last_z_dt == DT_FASTQ ? 3 : 4);
        PRINT (bam_seq_to_sam, 2);
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
        PRINT (vcf_seg_PROBE_A, 2);
        PRINT (random_access_merge_in_vb, 1); 
        PRINT (gencomp_absorb_vb_gencomp_lines, 1);
        PRINT (vcf_linesort_merge_vb, 1);
        PRINT (ctx_sort_dictionaries_vb_1, 1); 
        PRINT (ctx_merge_in_vb_ctx, 1);
        PRINT (zip_compress_ctxs, 1);
        PRINT (zip_generate_b250, 2);
        PRINT (zip_generate_local, 2);
        print_ctx_compressor_times();
        PRINT (codec_assign_best_codec, 2);
        PRINT (compressor_bz2,  2);
        PRINT (compressor_lzma, 2);
        PRINT (compressor_bsc,  2);
        PRINT (compressor_rans, 2);
        PRINT (compressor_arith, 2);
        PRINT (compressor_domq, 2);
        PRINT (compressor_normq, 2);
        PRINT (compressor_actg, 2);
        PRINT (compressor_pbwt, 2);
        PRINT (compressor_longr, 2);
        PRINT (codec_hapmat_count_alt_alleles, 2);
        PRINT (dict_io_compress_dictionaries, 1); 
        PRINT (dict_io_assign_codecs, 2); 
        PRINT (dict_io_compress_one_fragment, 2);
        PRINT (recon_plan_compress, 1);
        PRINT (recon_plan_compress_one_fragment, 2);
        PRINT (ref_compress_ref, 1);
        PRINT (ref_compress_one_range, 2);
        PRINT (ref_copy_compressed_sections_from_reference_file, 2);
        PRINT (random_access_finalize_entries, 1);
        PRINT (random_access_compress, 1);
        PRINT (ctx_compress_counts, 1);
        PRINT (zfile_compress_genozip_header, 1);
        PRINT (sam_zip_prim_ingest_vb, 1);
        PRINT (digest, 1);
        PRINT (bgzf_compute_thread, 1);
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
        iprintf ("  Average read time: %u ms\n", ms(profile.read) / profile.num_vbs);
        iprintf ("  Average compute time: %u ms\n", ms(profile.compute) / profile.num_vbs);
        iprintf ("  Average write time: %u ms\n", ms(profile.write) / profile.num_vbs);
    }
    
    iprint0 ("\n\n");

    mutex_show_bottleneck_analsyis();
}