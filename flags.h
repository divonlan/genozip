// ------------------------------------------------------------------
//   flags.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

typedef enum { 
    REF_NONE          = 0,  // ZIP (except SAM) and PIZ when user didn't specify an external reference
    REF_INTERNAL      = 1,  // ZIP SAM only: use did not specify an external reference - reference is calculated from file(s) data
    REF_EXTERNAL      = 2,  // ZIP & PIZ: user specified --reference
    REF_EXT_STORE     = 4,  // ZIP: user specified --REFERENCE
    REF_MAKE_CHAIN    = 8,  // ZIP of a chain file
    REF_LIFTOVER      = 16, // ZIP: user specified --chain which uses a reference

    REF_STORED        = REF_INTERNAL | REF_EXT_STORE,                                 // ZIP/PIZ: file contains REFERENCE sections (==REF_STORED is only in PIZ; in ZIP only one of the bits is set)
    REF_ZIP_LOADED    = REF_EXTERNAL | REF_EXT_STORE | REF_MAKE_CHAIN | REF_LIFTOVER, // ZIP: external reference is loaded
    REF_ZIP_CHROM2REF = REF_EXTERNAL | REF_EXT_STORE | REF_LIFTOVER,                  // ZIP: chrom2ref mapping is stored
} ReferenceType;
extern rom ref_type_name(void);

#define IS_REF_INTERNAL   (flag.reference == REF_INTERNAL)
#define IS_REF_EXTERNAL   (flag.reference == REF_EXTERNAL)
#define IS_REF_EXT_STORE  (flag.reference == REF_EXT_STORE)
#define IS_REF_MAKE_CHAIN (flag.reference == REF_MAKE_CHAIN)
#define IS_REF_LIFTOVER   (flag.reference == REF_LIFTOVER)

#define IS_REF_STORED_PIZ (flag.reference == REF_STORED) // 2 bits set
#define IS_REF_INTERNAL_PIZ ((Z_DT(DT_SAM) || Z_DT(DT_BAM)) && z_file->z_flags.dts_ref_internal)

typedef enum { STATS_NONE=0, STATS_SHORT=1, STATS_LONG=2, STATS_SHORT_GREP=-1, STATS_LONG_GREP=-2 } StatsType;

typedef struct {
    
    // genozip options that affect the compressed file
    int fast, best, make_reference, multiseq, md5;
    rom vblock;
    
    // ZIP: data modifying options
    int optimize, optimize_sort, optimize_phred, GL_to_PL, GP_to_PP, optimize_VQSLOD,  // optimize flags
        optimize_QUAL, optimize_Vf, optimize_ZM, optimize_DESC,
        allow_ambiguous, add_line_numbers, match_chrom_to_reference;
    
    char *dvcf_rename, *dvcf_drop;
    
    #define NOT_PAIRED_END 0 
    #define PAIR_READ_1    1
    #define PAIR_READ_2    2
    int pair; // unfortunately we can't rely on enum being sizeof(int)
    int undocumented_dts_paired; // undocumented flag to uncompress paired files older than 9.0.13 when genozip_header.dts_paired was introduced

    // genounzip options
    int32_t bgzf;   // can be set by --bgzf, or by various other conditions. values 0-12 indicate the level of libdeflate, BGZF_BY_ZFILE means use SEC_BGZF or default level if it is absent
    
    int out_dt; // used to indicate the desired dt of the output txt - consumed by file_open, and thereafter equal to txt_file->data_type
    #define out_dt_is_binary dt_props[flag.out_dt].is_binary

    // PIZ: data-modifying genocat options for showing only a subset of the file, or otherwise modify the file 
    int header_one, header_only_fast, no_header, header_only, // how to handle the txt header
        seq_only, qual_only, single_coord, 
        regions, gpos, samples, 
        drop_genotypes, gt_only, luft, sort, unsorted, snps_only, indels_only, // VCF options
        sequential, no_pg, extended_translation,  one_component,
        kraken_taxid;
    rom regions_file;
    int64_t lines_first, lines_last, tail;  // set by --head, --tail, --lines 
    rom grep; int grepw; unsigned grep_len; // set by --grep and --grep-w
    uint32_t one_vb, downsample, shard ;
    enum { SAM_FLAG_INCLUDE_IF_ALL=1, SAM_FLAG_INCLUDE_IF_NONE, SAM_FLAG_EXCLUDE_IF_ALL } sam_flag_filter;
    enum { SAM_MAPQ_INCLUDE_IF_AT_LEAST=1, SAM_MAPQ_EXCLUDE_IF_AT_LEAST } sam_mapq_filter;
    enum { INTERLEAVE_NONE, INTERLEAVE_EITHER, INTERLEAVE_BOTH } interleaved;
    uint16_t FLAG; // the value for sam_flag_filter
    uint8_t MAPQ;  // the value for sam_mapq_filter
    enum { IUP_NONE, IUP_POSITIVE, IUP_NEGATIVE } bases;

    // genols options
    int bytes;

    // options affecting the software interaction (but not the file contents)
    int force, quiet, show_filename,
        to_stdout,   // set implicitly if genocat without --output
        replace, 
        lic_width,   // width of license output, 0=dynamic (undocumented parameter of --license)
        test, no_test,       
        index_txt,   // create an index
        subdirs,     // recursively traversing subdirectories
        list;        // a genols option
    rom threads_str, out_filename, out_dirname, files_from, do_register;

    ReferenceType reference;

    // stats / metadata flags for end users
    int show_dvcf, show_ostatus, show_lift, show_wrong_md, show_wrong_xg, show_wrong_xm, show_wrong_xb; 
    enum { VLD_NONE, VLD_REPORT_INVALID, VLD_REPORT_VALID, VLD_INVALID_FOUND } validate; // genocat: tests if this is a valid genozip file (z_file opens correctly)
    StatsType show_stats;
    CompIType show_stats_comp_i;

    // analysis
    int list_chroms, show_sex, idxstats;
    enum { CNT_NONE, CNT_TOTAL, COUNT_VBs } count; 
    enum { COV_NONE, COV_ALL, COV_CHROM, COV_ONE } show_coverage;
    enum { KRK_NONE, KRK_ALL, KRK_INCLUDED, KRK_EXCLUDED } show_kraken;

    // stats / debug useful mostly for developers
    int debug, show_memory, show_sag, show_depn, show_dict, show_b250, show_aliases, show_digest, log_digest, show_recon_plan,
        show_index, show_gheader, show_ref_contigs, show_chain_contigs, show_ref_seq, show_ref_diff,
        show_reference, show_ref_hash, show_ref_index, show_chrom2ref, show_ref_iupacs, show_chain, show_ranges,
        show_codec, 
        show_alleles, show_bgzf, show_txt_contigs, show_lines,
        show_vblocks, show_threads, show_uncompress, biopsy,
        debug_progress, show_hash, debug_memory, debug_threads, debug_stats, debug_generate, debug_recon_size, debug_seg,
        debug_LONG, show_qual, debug_qname, debug_read_ctxs, debug_sag, debug_gencomp, debug_lines, debug_latest,
        debug_peek,
        no_gencomp, force_gencomp, no_domqual, verify_codec, seg_only, show_bam, xthreads, show_flags, show_rename_tags,
        #define SHOW_CONTAINERS_ALL_VBs (-1)
        show_containers, show_aligner, show_buddy,
        echo,         // show the command line in case of an error
        show_headers; // (1 + SectionType to display) or 0=flag off or -1=all sections
    rom help, dump_section, show_is_set, show_time, show_mutex;
    struct biopsy_line { VBIType vb_i; int32_t line_i/*within vb*/; } biopsy_line; // argument of --biopsy-line (line_i=-1 means: not used)
    
    DictId dict_id_show_one_b250,   // argument of --show-b250-one
           show_one_counts,
           dump_one_b250_dict_id,   // argument of --dump-b250-one
           dump_one_local_dict_id,  // argument of --dump-local-one
           dict_id_show_containers, // argument of --show-containers
           dict_id_debug_seg;       // argument of --debug-seg
    rom show_one_dict;      // argument of --show-dict-one

    #define HAS_DEBUG_SEG(ctx) (flag.debug_seg && (!flag.dict_id_debug_seg.num || dict_id_typeless ((ctx)->dict_id).num == flag.dict_id_debug_seg.num))

    // internal flags set by the system, not the command line
    CompIType zip_comp_i;    // ZIP only: currently zipping component (copied into VBlock.comp_i, FlagsTxtHeader.comp_i, SectionEntFileFormat.comp_i)
    bool debug_top,
         is_windows, is_mac, is_linux, // set according to OS
         is_lten,            // set according to endianness   
         explicit_out_dt,    // genocat - out txt file data type set explicitly from command line
         aligner_available,  // ZIP: compression requires using the aligner
         genocat_no_ref_file,// PIZ (genocat): we don't need to load the reference data
         genocat_no_dicts,   // PIZ (genocat): we don't need to read the dicts
         genocat_global_area_only, // PIZ (genocat): we quit after processing the global area
         genocat_no_reconstruct,   // PIZ: User requested to genocat with only metadata to be shown, not file contents
         no_writer,          // PIZ: User requested to genocat with only metadata to be shown, not file contents (but we still might do reconstruction without output)
         no_writer_thread,   // PIZ: Don't use a Writer thread. Sometimes when no_writer, we still need a writer thread (eg to calculate digest with SAM generated components)
         zip_no_z_file,      // ZIP: a mode where no z_file is created (eg biopsy)
         preprocessing,      // PIZ: we're currently dispatching compute threads for preprocessing (= loading SA Groups)
         multiple_files,     // Command line includes multiple files
         reconstruct_as_src, // the reconstructed data type is the same as the source data type
         maybe_txt_header_modified,
         maybe_lines_dropped_by_reconstructor,
         maybe_lines_dropped_by_writer,
         maybe_vb_modified_by_reconstructor,
         maybe_lines_out_of_order,
         maybe_vb_dropped_by_writer,
         maybe_vb_dropped_after_read,
         missing_contexts_allowed, // PIZ: its not an error if contexts are missing - just reconstruct as an empty string
         data_modified,      // PIZ: output is NOT precisely identical to the compressed source, and hence we cannot use its BZGF blocks
                             // ZIP: txt data is modified during Seg
         explicit_ref,       // ref->filename was set by --reference or --REFERENCE (as opposed to being read from the genozip header)
         collect_coverage;   // PIZ: collect coverage data for show_sex/show_coverage/idxstats

    int only_headers,        // genocat --show_headers (not genounzip) show only headers
        check_latest;        // PIZ: run with "genozip --decompress --test": ZIP passes this to PIZ upon testing of the last file

    Reference reading_reference;  // system is currently reading a reference  as a result of --chain (not normal PIZ of a .chain.genozip)

#define flag_loading_auxiliary (flag.reading_reference || flag.reading_chain || flag.reading_kraken) // PIZ: currently reading auxiliary file (reference, chain etc)

    rom reading_chain; // system is currently loading a chain file by this name
    rom reading_kraken;// system is currently loading a kraken file by this name
    rom unbind;
    rom log_filename;  // output to info_stream goes here

    enum { BIND_NONE, BIND_FQ_PAIR, BIND_DVCF, BIND_SAM } bind; // ZIP: cases where we have more than one txt_file bound in a z_file
    uint64_t stdin_size;
    unsigned longest_filename; // length of longest filename of the txt/z files on the command line
    
} Flags;

extern Flags flag;

#define NO_LINE (-1) // used for lines_first, lines_last, biopsy_lines, VB_SAM->mate_line_i, VB_SAM->saggy_line_i

#define SAVE_FLAGS Flags save_flag = flag
#define RESTORE_FLAGS flag = save_flag

// save and reset flags that are intended to operate on the compressed file rather than the reference file
#define SAVE_FLAGS_AUX Flags save_flag = flag; \
    flag.test = flag.md5 = flag.show_memory = flag.show_stats = flag.no_header = flag.show_bgzf =\
    flag.header_one = flag.header_only = flag.regions = flag.show_index = flag.show_dict =  \
    flag.show_b250 = flag.show_ref_contigs = flag.list_chroms = flag.count = \
    flag.downsample = flag.shard = flag.one_vb = flag.one_component = flag.xthreads = \
    flag.show_sex = flag.show_coverage = flag.idxstats = flag.collect_coverage = 0; /* int */ \
    flag.bases = IUP_NONE; \
    flag.interleaved = INTERLEAVE_NONE; \
    flag.grep = flag.show_time = flag.unbind = flag.show_one_dict = flag.out_filename = NULL; /* char* */ \
    flag.dict_id_show_one_b250 = flag.dump_one_b250_dict_id = flag.dump_one_local_dict_id = DICT_ID_NONE; /* DictId */ 
    
#define SAVE_FLAG(f) typeof(flag.f) save_##f = flag.f 
#define TEMP_FLAG(f,value) SAVE_FLAG(f) ; flag.f=(typeof(flag.f))(uint64_t)value
#define SET_FLAG(f) SAVE_FLAG(f) ; flag.f=(typeof(flag.f))(uint64_t)1
#define CLEAR_FLAG(f) SAVE_FLAG(f) ; flag.f=(typeof(flag.f))(uint64_t)0
#define RESTORE_FLAG(f) flag.f = save_##f

// check for incompatabilities between flags
extern bool option_is_short[];
#define OT(l,s) option_is_short[(int)s[0]] ? "-"s : "--"l

extern void flags_init_from_command_line (int argc, char **argv);
extern void flags_update (unsigned num_files, rom *filenames);
extern void flags_update_zip_one_file (void);
extern void flags_update_piz_one_file (int z_file_i);

extern void flags_store_command_line (int argc, char **argv);
extern rom flags_command_line (void);
extern rom flags_pipe_in_process_name (void);
extern unsigned flags_pipe_in_pid (void);
extern bool flags_pipe_in_process_died (void);
extern bool flags_is_genocat_global_area_only (void);
