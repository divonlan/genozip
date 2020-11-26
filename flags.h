// ------------------------------------------------------------------
//   flags.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FLAGS_INCLUDED
#define FLAGS_INCLUDED

#include "genozip.h"

typedef struct {
    
    // genozip options that affect the compressed file
    int vblock, gtshark, fast, make_reference, md5;

    // ZIP: data modifying optinos
    int optimize, optimize_sort, optimize_PL, optimize_GL, optimize_GP, optimize_VQSLOD,  // optimize flags
        optimize_QUAL, optimize_Vf, optimize_ZM, optimize_DESC;
    
    #define NOT_PAIRED_END 0 
    #define PAIR_READ_1    1
    #define PAIR_READ_2    2
    int pair; // unfortunately we can't rely on enum being sizeof(int)

    // genounzip options
    int plain, bgzf, out_dt;
    char *unbind;

    // PIZ: data-modifying genocat options for showing only a subset of the file 
    int header_one, no_header, header_only, // how to handle the txt header
        regions, samples, drop_genotypes, gt_only, sequential, no_pg;
    char *grep;

    // genols options
    int bytes;

    // options affecting the software interaction (but not the file contents)
    int force, quiet, to_stdout, replace, do_register,
        test; // implies md5
    char *threads_str, *out_filename;

    enum { REF_NONE,      // ZIP (except SAM) and PIZ when user didn't specify an external reference
           REF_INTERNAL,  // ZIP SAM only: use did not specify an external reference - reference is calculated from file(s) data
           REF_EXTERNAL,  // ZIP & PIZ: user specified --reference
           REF_EXT_STORE, // ZIP: user specified --REFERENCE
           REF_STORED     // PIZ: file contains REFERENCE sections (user cannot specify --reference)
    } reference;

    // stats / metadata flags for end users
    int list_chroms, show_stats; 
    
    // stats / debug useful mostly for developers
    int show_memory, show_dict, show_b250, show_headers, show_aliases, show_digest,
        show_index, show_gheader, show_ref_contigs, show_ref_seq,
        show_reference, show_ref_hash, show_ref_index, show_ref_alts,
        show_codec, show_containers, show_alleles, show_bgzf, show_txt_contigs,
        debug_progress, show_hash, debug_memory, show_vblocks, show_threads,
        test_seg;
    char *dump_section, *show_is_set, *show_time;

    DictId dict_id_show_one_b250,   // argument of --show-b250-one
           dict_id_show_one_dict,   // argument of --show-dict-one
           dump_one_b250_dict_id,   // argument of --dump-b250-one
           dump_one_local_dict_id;  // argument of --dump-local-one

    // internal flags set by the system, not the command line
    bool bind,               // ZIP: user used --output to bind 2+ files
         ref_use_aligner,    // ZIP: compression requires using the aligner
         const_chroms,       // ZIP: chroms dictionary created from reference or file header and no more chroms can be added
         reading_reference,  // system is currently reading a reference file
         do_translate,       // PIZ: decompression requires translation to another data type
         genocat_info_only,  // User requested to genocat with only metadata to be shown, not file contents
         multiple_files,     // Command line includes multiple files
         reconstruct_as_src, // the reconstructed data type is the same as the source data type
         data_modified;      // PIZ: output is NOT precisely identical to the compressed source, and hence we cannot use its BZGF blocks

    uint64_t stdin_size;
} Flags;

extern Flags flag;

#define SAVE_FLAGS Flags save_flag = flag
#define RESTORE_FLAGS flag = save_flag

#define SAVE_FLAG(f) typeof(flag.f) save_##f = flag.f 
#define RESET_FLAG(f) SAVE_FLAG(f) ; flag.f=(typeof(flag.f))(uint64_t)0
#define RESTORE_FLAG(f) flag.f = save_##f

// check for incompatabilities between flags
#define OT(l,s) is_short[(int)s[0]] ? "-"s : "--"l

extern void flags_init_from_command_line (int argc, char **argv, bool *is_short);
extern void flags_update (unsigned num_files, char **filenames, const bool *is_short);
extern void flags_update_zip_one_file (void);
extern void flags_update_piz_one_file (void);

#endif