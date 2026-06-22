// ------------------------------------------------------------------
//   flags.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "file_types.h"
#include "endianness.h"
#include "dict_id.h"

typedef packed_enum { 
    REF_NONE          = 0,  // ZIP (except SAM) and PIZ when user didn't specify an external reference
    REF_INTERNAL      = 1,  // ZIP SAM only: use did not specify an external reference - reference is calculated from file(s) data
    REF_EXTERNAL      = 2,  // ZIP & PIZ: user specified --reference
    REF_EXT_STORE     = 4,  // ZIP: user specified --REFERENCE

    REF_STORED        = REF_INTERNAL | REF_EXT_STORE, // ZIP/PIZ: file contains REFERENCE sections (==REF_STORED is only in PIZ; in ZIP only one of the bits is set)
    REF_ZIP_LOADED    = REF_EXTERNAL | REF_EXT_STORE, // ZIP: external reference is loaded
    REF_ZIP_CHROM2REF = REF_EXTERNAL | REF_EXT_STORE, // ZIP: chrom2ref mapping is stored
} ReferenceType;
extern rom ref_type_name(void);

#define IS_REF_INTERNAL   (flag.reference == REF_INTERNAL)
#define IS_REF_EXTERNAL   (flag.reference == REF_EXTERNAL)
#define IS_REF_EXT_STORE  (flag.reference == REF_EXT_STORE)

#define IS_REF_LOADED_ZIP (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE)
#define IS_REF_CHROM2REF  (flag.reference & REF_ZIP_CHROM2REF)

#define IS_REF_STORED_PIZ (flag.reference == REF_STORED) // 2 bits set
#define IS_REF_INTERNAL_PIZ ((Z_DT(BAM) || Z_DT(SAM)) && z_file->z_flags.is_ref_internal)

typedef enum { STATS_NONE=0, STATS_SHORT=1, STATS_LONG=2, STATS_SHORT_GREP=-1, STATS_LONG_GREP=-2 } StatsType;

typedef packed_enum { NOT_PAIRED,       // ZIP and PIZ
                      PAIR_R1, PAIR_R2, // ZIP: currently compressing R1 or R2 
                      PAIRED,           // PIZ: z_file is paired ; ZIP: --pair or --deep with paired FASTQ
                    } PairType; 
#define PAIR_TYPE_NAMES { "NOT_PAIRED", "PAIR_R1", "PAIR_R2", "PAIRED" }
#define IS_R1 (flag.pair == PAIR_R1 && !flag.preprocessing) // ZIP only
#define IS_R2 (flag.pair == PAIR_R2 && !flag.preprocessing) // ZIP only


// Notes:
// - gcc_struct causes mingw struct to be the same as Linux. Specifically, bit-fields with a different type DO NOT start a new aligned field.
// - Flag is NOT ((packed)), but we take care to order members so that it is as packed as possible
#ifdef __clang__
typedef struct { 
#else //gcc
typedef struct __attribute__((gcc_struct)) { 
#endif
    // _____________________________________________________________________________________________________________________
    // HIGH FREQUENCY FLAGS: accessed during seg/recon AT LEAST ONCE PER LINE during NORMAL (no special flags) seg/recon 
    // Fit everything in a single 64B L1 cache line that will become pinned during seg/recon

    // recon - used in every line
    uint64_t maybe_lines_dropped_by_reconstructor : 1;
    uint64_t maybe_lines_dropped_by_writer : 1;
    uint64_t missing_contexts_allowed      : 1; // PIZ: its not an error if contexts are missing - just reconstruct as an empty string
    uint64_t genocat_no_reconstruct : 1;   // PIZ: User requested to genocat with only metadata to be shown, not file contents
    packed_enum { IUP_NONE, IUP_POSITIVE, IUP_NEGATIVE } bases : 2;
    uint64_t collect_coverage   : 1;        // PIZ: collect coverage data for show_coverage/idxstats
    uint64_t show_lines         : 1;
    uint64_t show_stack         : 1;
    uint64_t show_snips         : 1;
    uint64_t seq_only           : 1;
    uint64_t qual_only          : 1;
    uint64_t show_aligner       : 1;
    uint64_t header_only_fast   : 1;
    int8_t qname_filter         : 2; // possible values: -1, 0, 1
    DataType out_dt             : 6; // (within on byte) used to indicate the desired dt of the output txt - consumed by file_open_z, and thereafter equal to txt_file->data_type

    // seg & recon - used in both seg and recon of every line 
    PairType pair               : 2; 
    CommandType command         : 8; // (on byte boundary) command running now (eg genozip -d and genounzip are both PIZ)
    CompIType show_time_comp_i  : 8; // (on byte boundary) comp_i for which to show time (possibly COMP_NONE or COMP_ALL) - used by COPY_TIMER    
    ReferenceType reference     : 3; // (within byte) 
    uint64_t deep               : 1; // deep is set with --deep in ZIP and from SectionHeaderGenozipHeader.flags.genozip_header.is_deep in PIZ

    // seg - used in every line
    uint64_t aligner_available  : 1; // ZIP: compression requires using the aligner
    packed_enum { SHOW_DEEP_SUMMARY=1, SHOW_DEEP_ONE_HASH=2, SHOW_DEEP_ALL=3 } show_deep : 2;
    uint64_t fast               : 1; 
    uint64_t best               : 1; 
    uint64_t low_memory         : 1;
    
    // developer options that are nonetheless tested in normal seg/recon
    uint64_t debug_seg          : 1;
    uint64_t show_depn          : 1;
    uint64_t debug_split        : 1;
    uint64_t debug_peek         : 1;
    uint64_t show_wrong_md      : 1;
    uint64_t show_wrong_xg      : 1;
    uint64_t show_wrong_xm      : 1;
    uint64_t show_wrong_xb      : 1;
    uint64_t show_tlen_pred     : 1;
    uint64_t debug_dyn_int      : 1;
    uint64_t debug_aligner      : 1;
    uint64_t show_buddy         : 1;
    uint64_t debug_sag          : 1;
    uint64_t debug_lines        : 1;
    // ——— 64 bits up to here ———

    #define flag_after_bits bam_assist // first field after 64-bits
    rom bam_assist, show_time, show_mutex; 
    struct biopsy_line { VBIType vb_i; int32_t line_i/*within vb*/; } biopsy_line; // argument of --biopsy-line (line_i=-1 means: not used)
    
    // other high frequency flags
    uint64_t no_BDBI            : 1; // ZIP SAM: disable the BD/BI method 

    // _____________________________________________________________________________________________________________________

    
    // Other bits

    // OTHER BIT FLAGS: not accessed in every line in NORMAL seg/recon
    uint64_t is_windows         : 1; 
    uint64_t is_mac             : 1; 
    uint64_t is_linux           : 1; 
    uint64_t is_wsl             : 1; 
    uint64_t is_lten            : 1; // set according to endianness  
    uint64_t is_sanitize_thread : 1; // build includes -fsanitize=thread 
    uint64_t is_valgrind        : 1; // running under valgrind
    uint64_t is_docker          : 1; // running in a docker container
    uint64_t explicit_out_dt    : 1; // genocat - out txt file data type set explicitly from command line
    packed_enum { NO_PREPROC, PREPROC_RUNNING, PREPROC_FINALIZING } preprocessing : 2; // we're currently dispatching compute threads for preprocessing (PIZ: loading SA Groups, ZIP: loading bamass ents)
    uint64_t dont_load_ref_file : 1; // PIZ (genocat): we don't need to load the reference data
    uint64_t genocat_no_dicts   : 1; // PIZ (genocat): we don't need to read the dicts
    uint64_t genocat_global_area_only : 1; // PIZ (genocat): we quit after processing the global area
    uint64_t no_writer          : 1; // PIZ: User requested to genocat with only metadata to be shown, not file contents (but we still might do reconstruction without output)
    uint64_t no_writer_thread   : 1; // PIZ: Don't use a Writer thread. Sometimes when no_writer, we still need a writer thread (eg to calculate digest with SAM generated components)
    uint64_t zip_no_z_file      : 1; // ZIP: a mode where no z_file is created (eg biopsy)
    uint64_t multiple_files     : 1; // Command line includes multiple files
    uint64_t reconstruct_as_src : 1; // the reconstructed data type is the same as the source data type
    uint64_t has_reconstructor_filter : 1;
    uint64_t maybe_vb_modified_by_reconstructor : 1;
    uint64_t maybe_lines_out_of_order : 1;
    uint64_t piz_txt_modified   : 1; // PIZ: output is NOT precisely identical to the compressed source, and hence we cannot use its BZGF blocks or verify digest
    uint64_t zip_lines_counted_at_init_vb : 1; // ZIP: VB lines need to be counted at zip_init_vb instead of zip_update_txt_counters, requiring BGZF-uncompression of a VB by the main thread
    uint64_t explicit_ref       : 1; // ref->filename was set by --reference or --REFERENCE (as opposed to being read from the genozip header)
    uint64_t deep_fq_only       : 1; // PIZ: SAM data is reconstructed by not written; only FASTQ data is written
    uint64_t deep_sam_only      : 1; // PIZ: only SAM data is needed, no need to create deep_ents
    uint64_t deep_no_qual       : 1; // ZIP: don't deep QUAL data
    uint64_t removing_cache     : 1; // genocat: running --no-cache with only -e 𝑟𝑒𝑓-𝑓𝑖𝑙𝑒 to remove cache
    uint64_t let_OS_cleanup_on_exit : 1; // don't release resources as we are about to exit - the OS does it faster
    uint64_t reading_reference  : 1; // system is currently reading a reference
    uint64_t md5                : 1;
    uint64_t secure_DP          : 1; 
    uint64_t not_paired         : 1;
    uint64_t default_make_ref   : 1;

    // File modifying options
    enum { NO_OPTIMIZE, FULL_OPTIMIZE, OPTIMIZE_POS_LIST, OPTIMIZE_NEG_LIST } optimize : 2;
    uint64_t add_line_numbers   : 1; 
    uint64_t add_seq            : 1;        
    uint64_t truncate           : 1; // allow truncated file - compress only available full lines. note: we don't consider this option data modifying as its used for debugging - digest is calculated only after truncation

    // ZIP options
    uint64_t explicitly_generic : 1; // user explicitly set the type to generic
    enum { BIND_NONE, BIND_FQ_PAIR, BIND_SAM, BIND_DEEP } bind : 2; // ZIP: cases where we have more than one txt_file bound in a z_file
    
    // PIZ: data-modifying genocat options for showing only a subset of the file, or otherwise modify the file 
    uint64_t header_one         : 1; // how to handle the txt header
    uint64_t no_header          : 2; // 0, 1 or 2
    uint64_t header_only        : 1;
    uint64_t qname_only         : 1;
    uint64_t regions            : 1;
    uint64_t gpos               : 1;
    uint64_t samples            : 1;
    int8_t seq_filter           : 2; // 1 positive, -1 negative filter
    uint64_t drop_genotypes     : 1;
    uint64_t gt_only            : 1;
    uint64_t snps_only          : 1; // VCF options
    uint64_t indels_only        : 1;
    uint64_t sequential         : 1;
    uint64_t no_pg              : 1;
    uint64_t explicit_head      : 1;
    enum { SAM_FLAG_INCLUDE_IF_ALL=1, SAM_FLAG_INCLUDE_IF_NONE, SAM_FLAG_EXCLUDE_IF_ALL } sam_flag_filter : 2;
    enum { SAM_MAPQ_INCLUDE_IF_AT_LEAST=1, SAM_MAPQ_EXCLUDE_IF_AT_LEAST } sam_mapq_filter : 2;
    enum { INTERLEAVE_NONE, INTERLEAVE_EITHER, INTERLEAVE_BOTH } interleaved : 2;

    uint64_t check_latest       : 1; // PIZ: run with "genozip --decompress --test": ZIP passes this to PIZ upon testing of the last file

    #define MAKE_REF_DEFAULT MAKE_REF_MEDIUM
    MakeRefSize make_reference  : 3; // note: values are part of the file format: GenozipHeader.ref.make_ref_size

    // genols options
    uint64_t ls_bytes           : 1; 
    uint64_t ls_cache           : 1;
    uint64_t list               : 1; // a genols option

    // analysis
    uint64_t show_contigs       : 1;
    uint64_t idxstats           : 1;
    enum { CNT_NONE, CNT_TOTAL, CNT_VBs } count : 2; 
    enum { COV_NONE, COV_ALL, COV_CHROM, COV_ONE } show_coverage : 2;
    enum { TELEMETRY_OFF, TELEMETRY_SEND, TELEMETRY_FILE } telemetry : 2;

    // diagnostics
    enum { VLD_NONE, VLD_REPORT_INVALID, VLD_REPORT_VALID, VLD_INVALID_FOUND, VLD_NO_REPORT } validate : 3; // genocat: tests if this is a valid genozip file (z_file opens correctly)
    StatsType show_stats        : 3;
    uint64_t show_gheader       : 2; // =1 show gheader as in file, =2 show shows section list after possible modiciation by writer_create_plan 
    uint64_t debug              : 1;
    uint64_t debug_or_test      : 1;
    uint64_t show_dict          : 1;
    uint64_t show_b250          : 1;
    uint64_t show_aliases       : 1;
    uint64_t show_digest        : 1;
    uint64_t log_digest         : 1;
    uint64_t show_recon_plan    : 1;
    uint64_t show_index         : 1;
    uint64_t show_reading_list  : 1;
    uint64_t show_ref_contigs   : 1;
    uint64_t show_ref_seq       : 1;
    uint64_t show_reference     : 1;
    uint64_t show_ref_hash      : 1;
    uint64_t show_ref_index     : 1;
    uint64_t show_chrom2ref     : 1;
    uint64_t show_ref_iupacs    : 1;
    uint64_t show_ranges        : 1;
    uint64_t show_codec         : 1;
    uint64_t show_cache         : 1;
    uint64_t show_memory        : 1;
    uint64_t show_regions       : 1;
    uint64_t show_seg_summary   : 1;
    uint64_t show_alleles       : 1;
    uint64_t show_bgzf          : 1;
    uint64_t show_isizes        : 1;
    uint64_t show_is_exactable  : 1;
    uint64_t show_txt_contigs   : 1;
    uint64_t show_gz_uncomp     : 1;
    uint64_t show_threads       : 1;
    uint64_t show_uncompress    : 1;
    uint64_t biopsy             : 1;
    uint64_t skip_segconf       : 1;
    uint64_t show_data_type     : 1;
    uint64_t show_tasks         : 1;
    uint64_t show_hash          : 1;
    uint64_t debug_threads      : 1;
    uint64_t debug_stats        : 1;
    uint64_t debug_generate     : 1;
    uint64_t debug_recon_size   : 1;
    uint64_t debug_LONG         : 1;
    uint64_t show_qual          : 1;
    uint64_t debug_qname        : 1; // accessed during qname discovery
    uint64_t debug_read_ctxs    : 1;
    uint64_t debug_gencomp      : 1;
    uint64_t debug_latest       : 1;
    uint64_t show_segconf_has   : 1;
    uint64_t debug_upgrade      : 1;
    uint64_t debug_expiration   : 1;
    uint64_t debug_debug        : 1; // a flag with no functionality - used for ad-hoc debugging  
    uint64_t debug_valgrind     : 1;
    uint64_t debug_tar          : 1;
    uint64_t debug_bai          : 1;
    uint64_t show_compress      : 1;
    uint64_t show_sec_gencomp   : 1;
    uint64_t show_scan          : 1;
    uint64_t no_gencomp         : 1; // accessed when setting Sag type (beginning of compression)
    uint64_t force_gencomp      : 1;
    uint64_t force_reread       : 1;
    uint64_t force_deep         : 1;
    uint64_t force_PLy          : 1; // VCF: accessed during segconf initialization
    uint64_t no_domqual         : 1;
    uint64_t no_pacb            : 1;
    uint64_t no_longr           : 1;
    uint64_t no_homp            : 1;
    uint64_t no_smux            : 1;
    uint64_t no_tmpl            : 1;
    uint64_t no_lzma            : 1;
    uint64_t no_bgzf            : 1;
    uint64_t no_faf             : 1;
    uint64_t no_interleaved     : 1;
    uint64_t no_splicing        : 1;
    uint64_t verify_codec       : 1;
    uint64_t seg_only           : 1;
    uint64_t show_bam           : 1;
    uint64_t no_index           : 1;
    uint64_t skip_index         : 1;
    uint64_t xthreads           : 1;
    uint64_t echo               : 2; // 0,1 or 2. show the command line in case of an error (including echo and its optional argument)
    uint64_t show_headers       : 1;
    uint64_t only_headers       : 1; // genocat --show_headers (not genounzip) show only headers

    uint64_t show_txt_offsets   : 1; 
    enum { SHOW_BAI_NONE, SHOW_BAI_UNSORTED, SHOW_BAI_SORT, SHOW_BAI_CHUNKS, SHOW_BAI_RAW, SHOW_BAI_LINEAR } show_bai : 3;

    // options affecting the software interaction (but not the file contents)
    uint64_t force              : 1; // accessed in segconf but not in actual segging
    uint64_t quiet              : 1;
    uint64_t explicit_quiet     : 1;
    uint64_t noisy              : 1;
    uint64_t no_tip             : 1;
    uint64_t show_filename      : 1;
    uint64_t analyze_ins        : 1;
    uint64_t to_stdout          : 1; // set implicitly if genocat without --output
    uint64_t replace            : 1;
    uint64_t test               : 1;
    uint64_t no_test            : 1;
    uint64_t explicit_test      : 1;
    uint64_t ext_indexing       : 1; // create an index file using an external indexer (samtools index etc)
    uint64_t subdirs            : 1; // recursively traversing subdirectories
    uint64_t no_cache           : 1; // don't load cache or delete cache
    uint64_t no_upgrade         : 1; // disable upgrade checks
    uint64_t no_eval            : 1; // don't allow features on eval basis (used for testing permissions)
    uint64_t from_url           : 1; // used for stats
    uint64_t do_activate        : 1; // activate license
    uint64_t restarted          : 1; // this process was run as a restart after an error (undocumented command line options)

    // ——— zip options (ordered from 64-bit types and down) ———
    rom vblock;
    int64_t sendto;
    uint64_t stdin_size;

    DictId *optimize_dict_ids;

// we track 64-bit words with ⸨ ... ⸩
/*⸨*/uint16_t optimize_dict_ids_len; // up to MAX_DICTS
    FileType stdin_type; // (1 byte) set by the --input command line option
    CompIType deep_num_fastqs;

    // ——— piz options (ordered from 64-bit types and down) ———
    uint32_t one_vb; 
/*⸩*/ 
    int64_t lines_last; // force member to be aligned
    rom regions_file, qnames_file, qnames_opt;
    int64_t lines_first, tail;  // set by --head, --tail, --lines 
    rom grep; int grepw; unsigned grep_len; // set by --grep and --grep-w
    int64_t t_offset, t_size;   // PIZ: offset and size of a z_file within a tar file - for genounzip --test executed from genozip
    rom unbind;         // unbinding prefix or ""
    rom log_filename;   // output to info_stream goes here

/*⸨*/uint32_t downsample, shard;
/*⸩*/ 
/*⸨*/uint16_t FLAG; // the value for sam_flag_filter
    
    #define MAX_FLAG_BGZF 5
    MgzipLevel bgzf;   // (1 byte) PIZ: can be set by --bgzf, or by various other conditions. values 0-MAX_FLAG_BGZF indicate a compress level out of bgzf_recompression_levels
        
    CompIType one_component; // (1 byte) 1-based ; 0=option unset (i.e. comp_i = one_component-1)

    CompIType one_vb_comp_i; // PIZ: COMP_NONE is --one-vb is
    uint8_t MAPQ;  // the value for sam_mapq_filter

    // ——— internal flags set by the system, not the command line ———
    CompIType zip_comp_i;    // ZIP only: currently zipping component (copied into VBlock.comp_i, FlagsTxtHeader.comp_i, SectionEntFileFormat.comp_i)
    volatile bool make_bai;  // PIZ BAM: write .bai file (volatile, because any thread can change its value at any time)
    /* 2 unused bits here */
/*⸩*/ 
    // ——— options affecting the software interaction (but not the file contents) ———
    rom test_i;      // test of test.sh currently running (undocumented)
    rom threads_str, out_filename, out_dirname, files_from;
    rom lic_param;   // format: width,type - invoked by Makefile
    rom license_filename;
    rom help;

    // ——— diagnostics ———
    rom dump_section, show_is_set, show_header_dict_name, show_flavor;

    struct biopsy_bytes { uint64_t start, len; } biopsy_bytes;
    
    DictId show_b250_δ,      // argument of --show-b250
           show_counts_δ,
           show_codec_δ,
           show_huffman_δ,
           debug_huffman_δ,
           show_singletons_δ,// argument of --show-singletons
           dump_b250_δ,      // argument of --dump-b250-one
           dump_local_δ,     // argument of --dump-local-one
           debug_seg_δ,      // argument of --debug-seg
           show_containers_δ;// argument of --show-containers
    rom show_one_dict;       // argument of --show-dict-one

    DeepHash debug_deep_hash;// qname, seq, qual hashes

/*⸨*/int show_sag; // -1=show all, >=1 - show grp_i=show_sag-1. accessed during ingest (ZIP) and load (PIZ)
    uint32_t debug_memory; // if > 1 show allocations at least this number of bytes
/*⸩*/ 
/*⸨*/int32_t dump_section_i;
    int32_t show_header_section_i; 
/*⸩*/ 
/*⸨*/int32_t dump_gz_block;

    #define SHOW_SEC_HEADER_ALL ((uint32_t)-1)
    uint32_t show_sec_headers; // bit array - bit per section type
/*⸩*/ 
/*⸨*/#define SHOW_CONTAINERS_ALL_VBs ((VBIType)-1)
    VBIType show_containers; // ALL_VBs, or vblock_i, or (ALL_VBs and show_containers_δ)
    
    CompIType show_stats_comp_i;
    
    Codec force_qual_codec;

    #define SHOW_VBLOCKS_ALL_TASKS ((Task)-1)
    Task show_vblocks;
    /* ⸨room for uint8⸩ */
/*⸩*/ 
} Flags;
extern Flags flag;

#define flag_loading_auxiliary (flag.reading_reference) // PIZ: currently reading auxiliary file (reference etc)

// use __builtin_expect for show/debug flags that are tested throughout execution (i.e. not just during initialization or finalization)
#define flag_show_deep          __builtin_expect (flag.show_deep > 0,   false)
#define flag_show_threads       __builtin_expect (flag.show_threads,    false)
#define flag_show_bgzf          __builtin_expect (flag.show_bgzf,       false)
#define flag_show_snips         __builtin_expect (flag.show_snips,      false)
#define flag_show_stack         __builtin_expect (flag.show_stack,      false)
#define flag_show_containers    __builtin_expect (flag.show_containers, false)
#define flag_show_memory        __builtin_expect (flag.show_memory,     false)
#define flag_debug_peek         __builtin_expect (flag.debug_peek,      false)
#define flag_debug_gencomp      __builtin_expect (flag.debug_gencomp,   false)
#define flag_has_biopsy_line    __builtin_expect (flag.biopsy_line.line_i != NO_LINE, false) // ZIP: --biopsy-line is used
#define flag_no_biopsy_line     __builtin_expect (flag.biopsy_line.line_i == NO_LINE, true)  // ZIP: --biopsy-line is not used

#define flag_show_sec_headers(st) ((flag.show_sec_headers >> (st)) & 1)

#define flag_has_head   (flag.lines_last != NO_LINE)  // ZIP: --head PIZ: --head or --lines is used
#define flag_has_tail   (flag.lines_last != NO_LINE || flag.tail)  // PIZ: --tail or --lines is used    
#define zip_is_biopsy   __builtin_expect (flag.biopsy || flag_has_biopsy_line, false) // ZIP: either --biopsy or --biopsy-line is used

#define deep_or_bamass (flag.deep ? "deep" : "bamass")

#define SAVE_FLAGS Flags save_flag = flag
#define RESTORE_FLAGS flags_restore (&save_flag)
extern void flags_restore (Flags *save_flag);

// save and reset flags that are intended to operate on the compressed file rather than the reference file
#define SAVE_FLAGS_AUX(aux_name) Flags save_flag = flag;                                                                        \
    if (flag.xthreads) { /* we assume user is interested only in PIZ/ZIP traffic */                                             \
        if ((aux_name) && (flag.show_threads || flag.debug_memory || flag.show_vblocks))                                        \
            WARN_ONCE (_FYI "xthreads: Turning off show_vblocks, show_threads, debug_memory while loading %s\n", (aux_name) ? (aux_name) : "");   \
        flag.show_threads = flag.debug_memory = false;                                                                          \
        flag.show_vblocks = TASK_NONE;                                                                                               \
    }                                                                                                                           \
    flag.test = flag.md5 = flag.show_memory = flag.show_stats = flag.no_header = flag.show_bgzf =                               \
    flag.header_one = flag.header_only = flag.regions = flag.show_index = flag.show_dict =                                      \
    flag.show_b250 = flag.show_ref_contigs = flag.show_contigs = flag.count =                                                   \
    flag.downsample = flag.shard = flag.one_vb = flag.one_component = flag.xthreads =                                           \
    flag.show_coverage = flag.idxstats = flag.collect_coverage = 0; /* int */                                   \
    flag.bases = IUP_NONE;                                                                                                      \
    flag.interleaved = INTERLEAVE_NONE;                                                                                         \
    flag.grep = flag.unbind = flag.show_one_dict = flag.out_filename = NULL; /* char* */                                        \
    flag.show_b250_δ = flag.dump_b250_δ = flag.dump_local_δ = flag.show_singletons_δ = DICT_ID_NONE; /* DictId */ 
    
#define SAVE_FLAG(f) int64_t save_##f = flag.f 
#define TEMP_FLAG(f,value) SAVE_FLAG(f) ; flag.f=value
#define SAVE_FLAG_(f) typeof(flag.f) save_##f = flag.f 
#define TEMP_FLAG_(f,value) SAVE_FLAG_(f) ; flag.f=(typeof(flag.f))(uint64_t)value
#define SET_FLAG(f)   SAVE_FLAG(f) ; flag.f=true
#define CLEAR_FLAG(f) SAVE_FLAG(f) ; flag.f=false
#define RESTORE_FLAG(f) flag.f = save_##f

// check for incompatabilities between flags
extern bool option_is_short[];
#define OT(l,s) (option_is_short[(int)s[0]] ? "-"s : "--"l)

extern void flags_init_from_command_line (int argc, char **argv);
extern void flags_finalize (void);
extern void flags_update (unsigned num_files, rom *filenames);
extern void flags_update_zip_one_file (void);
extern void flags_update_piz_one_z_file (int z_file_i);
extern void flags_update_piz_no_ref_file (void);
extern void flags_store_command_line (int argc, char **argv);
extern rom flags_command_line (void);
extern rom flags_pipe_in_process_name (void);
extern unsigned flags_pipe_in_pid (void);
extern bool flags_pipe_in_process_died (void);
extern bool flags_is_genocat_global_area_only (void);
extern rom pair_type_name (PairType p);
extern bool flags_writer_counts (void);

static inline bool flag_is_show_vblocks (Task task)
{
    return flag.show_vblocks && (!task || flag.show_vblocks == SHOW_VBLOCKS_ALL_TASKS || flag.show_vblocks == task);
}

static inline bool flag_is_set_(int flag_without_dict, DictId flag_with_dict, DictId dict/*may be DICT_ID_NONE*/)
{
    return flag_without_dict || (dict.num && flag_with_dict.num == dict_id_typeless(dict).num);    
}
#define flag_is_set(f, dict) flag_is_set_(flag.f, flag.f##_δ, (dict)) // true if either the "all" flag is set, or the specific dict_id (the *_δ flag)
#define flag_is_δ(f, dict)   flag_is_set_(0, flag.f##_δ, (dict))      // checks only the dict_id (not the "all" option)

#define HAS_DEBUG_SEG(ctx) flag_is_set (debug_seg, ctx->dict_id)
