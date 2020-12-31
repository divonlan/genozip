// ------------------------------------------------------------------
//   flags.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FLAGS_INCLUDED
#define FLAGS_INCLUDED

#include "genozip.h"

typedef struct {
    
    // genozip options that affect the compressed file
    int gtshark, fast, make_reference, md5;
    char *vblock;
    
    // ZIP: data modifying options
    int optimize, optimize_sort, optimize_PL, optimize_GL, optimize_GP, optimize_VQSLOD,  // optimize flags
        optimize_QUAL, optimize_Vf, optimize_ZM, optimize_DESC;
    
    #define NOT_PAIRED_END 0 
    #define PAIR_READ_1    1
    #define PAIR_READ_2    2
    int pair; // unfortunately we can't rely on enum being sizeof(int)

    // genounzip options
    #define FLAG_BGZF_BY_ZFILE -1
    int bgzf;   // can be set by --bgzf, or by various other conditions. values 0-12 indicate the level of libdeflate, FLAG_BGZF_BY_ZFILE means use SEC_BGZF or default level if it is absent
    
    int out_dt; // used to indicate the desired dt of the output txt - consumed by file_open, and thereafter equal to txt_file->data_type
    char *unbind;

    // PIZ: data-modifying genocat options for showing only a subset of the file 
    int header_one, header_only_fast, no_header, header_only, // how to handle the txt header
        regions, samples, drop_genotypes, gt_only, sequential, no_pg, interleave;
    char *grep;
    uint32_t one_vb, downsample;

    // genols options
    int bytes;

    // options affecting the software interaction (but not the file contents)
    int force, quiet, 
        to_stdout,   // redirect txt output upon decompression to stdout 
        replace, 
        do_register,
        test,        // implies md5
        index_txt;   // create an index
    char *threads_str, *out_filename;

    enum { REF_NONE,      // ZIP (except SAM) and PIZ when user didn't specify an external reference
           REF_INTERNAL,  // ZIP SAM only: use did not specify an external reference - reference is calculated from file(s) data
           REF_EXTERNAL,  // ZIP & PIZ: user specified --reference
           REF_EXT_STORE, // ZIP: user specified --REFERENCE
           REF_STORED     // PIZ: file contains REFERENCE sections (user cannot specify --reference)
    } reference;

    // undocumented options for internal use
    char *genobwa; // --genobwa=<contig-name> is used by the genobwa script in genozip / genocat to filter a fastq to a superset that includes all the reads that *might* be mapped to a chromosome 

    // stats / metadata flags for end users
    int list_chroms, show_stats; 
    
    // stats / debug useful mostly for developers
    int show_memory, show_dict, show_b250, show_aliases, show_digest,
        show_index, show_gheader, show_ref_contigs, show_ref_seq,
        show_reference, show_ref_hash, show_ref_index, show_ref_alts,
        show_codec, show_containers, show_alleles, show_bgzf, show_txt_contigs,
        debug_progress, show_hash, debug_memory, show_vblocks, show_threads,
        seg_only, xthreads,
        show_headers; // (1 + SectionType to display) or 0=flag off or -1=all sections
    char *help, *dump_section, *show_is_set, *show_time, *show_mutex;

    DictId dict_id_show_one_b250,   // argument of --show-b250-one
           dict_id_show_one_dict,   // argument of --show-dict-one
           dump_one_b250_dict_id,   // argument of --dump-b250-one
           dump_one_local_dict_id;  // argument of --dump-local-one

    // internal flags set by the system, not the command line
    bool ref_use_aligner,    // ZIP: compression requires using the aligner
         const_chroms,       // ZIP: chroms dictionary created from reference or file header and no more chroms can be added
         reading_reference,  // system is currently reading a reference file
         trans_containers,   // PIZ: decompression invokes container translators
         genocat_info_only,  // User requested to genocat with only metadata to be shown, not file contents
         multiple_files,     // Command line includes multiple files
         reconstruct_as_src, // the reconstructed data type is the same as the source data type
         data_modified;      // PIZ: output is NOT precisely identical to the compressed source, and hence we cannot use its BZGF blocks

    enum { BIND_NONE, BIND_ALL, BIND_PAIRS } bind; // ZIP: user used --output to bind all files or --pair without --output to bind every 2
    uint64_t stdin_size;

    // default max amount of txt data in each variant block. this is user-configurable with --vblock
    #define MAX_VBLOCK_MEMORY      2048       // in MB
    #define VBLOCK_MEMORY_MIN_DYN  (16   << 20) // VB memory - min/max when set in zip_dynamically_set_max_memory
    #define VBLOCK_MEMORY_MAX_DYN  (512  << 20) 
    #define VBLOCK_MEMORY_FAST     (16   << 20) // VB memory with --fast
    #define VBLOCK_MEMORY_MAKE_REF (1    << 20) // VB memory with --make-reference - reference data 
    #define VBLOCK_MEMORY_REFHASH  (16   << 20) // VB memory with --make-reference - refhash data (overridable with --vblock)
    #define VBLOCK_MEMORY_GENERIC  (16   << 20) // VB memory for the generic data type
    uint64_t vblock_memory;
} Flags;

extern Flags flag;
extern FILE *info_stream;

#define SAVE_FLAGS Flags save_flag = flag
#define RESTORE_FLAGS flag = save_flag

#define SAVE_FLAG(f) typeof(flag.f) save_##f = flag.f 
#define SET_FLAG(f) SAVE_FLAG(f) ; flag.f=(typeof(flag.f))(uint64_t)1
#define CLEAR_FLAG(f) SAVE_FLAG(f) ; flag.f=(typeof(flag.f))(uint64_t)0
#define RESTORE_FLAG(f) flag.f = save_##f

// check for incompatabilities between flags
extern bool option_is_short[];
#define OT(l,s) option_is_short[(int)s[0]] ? "-"s : "--"l

extern void flags_init_from_command_line (int argc, char **argv);
extern void flags_update (unsigned num_files, const char **filenames);
extern void flags_update_zip_one_file (void);
extern void flags_update_piz_one_file (void);

extern void flags_store_command_line (int argc, char **argv);
const BufferP flags_command_line (void);
const char *flags_pipe_in_process_name (void);
unsigned flags_pipe_in_pid (void);
bool flags_pipe_in_process_died (void);

#endif