// ------------------------------------------------------------------
//   genozip.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef GENOZIP_INCLUDED
#define GENOZIP_INCLUDED

#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#include <unistd.h> 
#else
#include "compatibility/visual_c_stdint.h"
#include "compatibility/visual_c_unistd.h"
#include "compatibility/visual_c_misc_funcs.h"
#endif

// -----------------
// system parameters
// -----------------
#define GENOZIP_EXT ".genozip"

#define MAX_POS ((PosType)0xffffffff) // maximum allowed value for POS (constraint: fit into uint32 ctx.local)

// default max amount of txt data in each variant block. this is user-configurable with --vblock
#define TXT_DATA_PER_VB_DEFAULT "16" // MB in default mode
#define TXT_DATA_PER_VB_FAST    "16" // MB with --fast
extern void vb_set_global_max_memory_per_vb (const char *mem_size_mb_str);

#define MAX_SUBFIELDS 100   // maximum number of VCF_FORMAT subfield types (except for GT), VCF_INFO, SAM_QNAME, SAM_OPTIONAL, GFF3_ATTRS subfield types that is supported in one line.
#define MAX_DICTS     253   // 254 is for future use and 255 is DID_I_NONE

#define DEFAULT_MAX_THREADS 8 // used if num_cores is not discoverable and the user didn't specifiy --threads

// ------------------------------------------------------------------------------------------------------------------------
// pointers used in header files - so we don't need to include the whole .h (and avoid cyclicity and save compilation time)
// ------------------------------------------------------------------------------------------------------------------------
typedef struct VBlock *VBlockP;
typedef const struct VBlock *ConstVBlockP;
typedef struct SubfieldMapper *SubfieldMapperP; 
typedef const struct SubfieldMapper *ConstSubfieldMapperP; 
typedef struct File *FileP;
typedef const struct File *ConstFileP;
typedef struct Buffer *BufferP;
typedef const struct Buffer *ConstBufferP;
typedef struct Structured *StructuredP;
typedef const struct Structured *ConstStructuredP;
typedef struct Context *ContextP;
typedef const struct Context *ConstContextP;
typedef struct MtfNode *MtfNodeP;
typedef const struct MtfNode *ConstMtfNodeP;
typedef struct SectionHeader *SectionHeaderP;
typedef struct SectionListEntry *SectionListEntryP;
typedef const struct SectionListEntry *ConstSectionListEntryP;
typedef struct Range *RangeP;
typedef struct BitArray *BitArrayP;
typedef const struct BitArray *ConstBitArrayP;
typedef struct RAEntry *RAEntryP;
typedef const struct RAEntry *ConstRAEntryP;

typedef enum { EXE_GENOZIP, EXE_GENOUNZIP, EXE_GENOLS, EXE_GENOCAT } ExeType;

#pragma pack(1) // structures that are part of the genozip format are packed.
#define DICT_ID_LEN    ((int)sizeof(uint64_t))    // VCF/GFF3 specs don't limit the field name (tag) length, we limit it to 8 chars. zero-padded. (note: if two fields have the same 8-char prefix - they will just share the same dictionary)
typedef union DictId {
    uint64_t num;            // num is just for easy comparisons - it doesn't have a numeric value and endianity should not be changed
    uint8_t id[DICT_ID_LEN]; // \0-padded IDs 
    uint16_t map_key;        // we use the first two bytes as they key into vb/z_file->dict_id_mapper
} DictId;
#pragma pack()

typedef uint8_t DidIType;   // index of a context in vb->contexts or z_file->contexts / a counter of contexts
typedef uint64_t CharIndex; // index within dictionary
typedef int32_t WordIndex;  // used for word and node indices
typedef int64_t PosType;    // used for position coordinate within a genome

// global parameters - set before any thread is created, and never change
extern uint32_t global_max_threads, global_max_memory_per_vb;
extern const char *global_cmd;            // set once in main()
extern ExeType exe_type;

// PIZ / ZIP inspired by "We don't sell Duff. We sell Fudd"
typedef enum { NO_COMMAND=-1, ZIP='z', PIZ='d' /* this is unzip */, LIST='l', LICENSE='L', VERSION='V', HELP='h', TEST_AFTER_ZIP } CommandType;
extern CommandType command, primary_command;

#define SAVE_FLAG(flag) typeof(flag) save_##flag = flag 
#define TEMP_FLAG(flag,temp) typeof(flag) save_##flag = flag ; flag = (temp)
#define RESET_FLAG(flag) SAVE_FLAG(flag) ; flag=(typeof(flag))(uint64_t)0
#define RESTORE_FLAG(flag) flag = save_##flag

// flags set by user's command line options
extern int flag_force, flag_quiet, flag_bind, flag_md5, flag_unbind, flag_show_alleles, flag_show_time, flag_bgzip, flag_bam, flag_bcf,
           flag_show_memory, flag_show_dict, flag_show_gt_nodes, flag_show_b250, flag_show_stats, flag_show_headers, flag_show_aliases,
           flag_show_index, flag_show_gheader, flag_show_ref_contigs, flag_stdout, flag_replace, flag_test, flag_regions,  
           flag_samples, flag_drop_genotypes, flag_no_header, flag_header_only, flag_show_threads, flag_list_chroms, 
           flag_show_vblocks, flag_gtshark, flag_sblock, flag_vblock, flag_gt_only, 
           flag_header_one, flag_fast, flag_multiple_files, flag_fasta_sequential, flag_register,
           flag_show_reference, flag_show_ref_hash, flag_show_ref_index, flag_show_ref_alts, flag_pair, flag_genocat_info_only, 
           flag_debug_progress, flag_show_hash, flag_debug_memory, flag_debug_no_singletons, flag_make_reference, flag_reading_reference,

           flag_optimize, flag_optimize_sort, flag_optimize_PL, flag_optimize_GL, flag_optimize_GP, flag_optimize_VQSLOD, 
           flag_optimize_QUAL, flag_optimize_Vf, flag_optimize_ZM, flag_optimize_DESC, flag_optimize_SEQ,

// flags set in code, that impact reference
            flag_ref_use_aligner, flag_ref_originates_from_internal;
           
           ;

// values of flag_reference
typedef enum { REF_NONE,      // ZIP (except SAM) and PIZ when user didn't specify an external reference
               REF_INTERNAL,  // ZIP SAM only: use did not specify an external reference - reference is calculated from file(s) data
               REF_EXTERNAL,  // ZIP & PIZ: user specified --reference
               REF_EXT_STORE, // ZIP: user specified --REFERENCE
               REF_STORED     // PIZ: file contains REFERENCE sections (user cannot specify --reference)
} ReferenceType;
extern ReferenceType flag_reference;           

// values of flag_pair
#define NOT_PAIRED_END 0
#define PAIR_READ_1    1
#define PAIR_READ_2    2

extern char *flag_grep, *flag_show_is_set;
extern uint64_t flag_stdin_size;

// external vb - used when an operation is needed outside of the context of a specific variant block;
extern VBlockP evb;

// macros
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b) )
#define MAX(a, b) (((a) > (b)) ? (a) : (b) )
#endif

// we defined these ourselves (normally defined in stdbool.h), as not always available on all platforms (namely issues with Docker Hub)
typedef _Bool bool;
#define true 1
#define false 0

#define SPECIAL(dt,num,name,func) \
    extern void func (VBlockP vb, ContextP ctx, const char *snip, unsigned snip_len); \
    static const int dt##_SPECIAL_##name = (num + 32); /* +32 to make it printable ASCII that can go into a snip */

// IMPORTANT: This is part of the genozip file format. 
// If making any changes, update arrays in 1. comp_compress 2. file_viewer 3. txtfile_estimate_txt_data_size
typedef enum __attribute__ ((__packed__)) { // 1 byte
    COMP_UNKNOWN=-1, 
    COMP_NONE=0, COMP_GZ=1, COMP_BZ2=2, COMP_LZMA=3, // internal compressors
    
    // novel codecs
    // compress a sequence of A,C,G,T nucleotides - first squeeze into 2 bits and then LZMA. It's about 25X faster and 
    // slightly better compression ratio than LZMA. Any characters that are not ACGT are stored in a complementary 
    // COMP_NON_ACGT compression - which is \0 for ACGT locations, \1 for acgt (smaller letters) locations and verbatim for other characters
    COMP_ACGT=10, COMP_NON_ACGT=11, 

    COMP_BGZ=20, COMP_XZ=21, COMP_BCF=22, COMP_BAM=23, COMP_CRAM=24, COMP_ZIP=25,  // external compressors

    NUM_COMPRESSION_ALGS
} CompressionAlg; 

// aligned with CompressionAlg ; used in --show-header
#define COMP_ALG_NAMES {"NONE", "GZ",   "BZ2",  "LZMA", "FFU4", "FFU5", "FFU6", "FFU7", "FFU8", "FFU9", \
                        "ACGT", "~CGT", "FF12", "FF13", "FF14", "FF15", "FF16", "FF17", "FF18", "FF19", \
                        "BGZ",  "XZ",   "BCF",  "BAM" , "CRAM", "ZIP" }

#define COMPRESSOR_CALLBACK(func) \
extern void func (VBlockP vb, uint32_t vb_line_i, \
                  char **line_data_1, uint32_t *line_data_len_1,\
                  char **line_data_2, uint32_t *line_data_len_2);
                  
// sanity checks
extern void main_exit (bool show_stack, bool is_error);
#define exit_on_error(show_stack) main_exit (show_stack, true)
#define exit_ok main_exit (false, false)

// check for a user error
#define ASSINP(condition, format, ...)       { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(false); }}
#define ASSINP0(condition, string)           { if (!(condition)) { fprintf (stderr, "\n%s\n", string); exit_on_error(false); }}

// check for a bug - prints stack
#define ASSERT(condition, format, ...)       { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); }}
#define ASSERT0(condition, string)           { if (!(condition)) { fprintf (stderr, "\n%s\n", string); exit_on_error(true); }}
#define ASSERTW(condition, format, ...)      { if (!(condition) && !flag_quiet) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); }}
#define ASSERTW0(condition, string)          { if (!(condition) && !flag_quiet) { fprintf (stderr, "\n%s\n", string); } }
#define RETURNW(condition, ret, format, ...) { if (!(condition)) { if (!flag_quiet) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); } return ret; }}
#define RETURNW0(condition, ret, string)     { if (!(condition)) { if (!flag_quiet) { fprintf (stderr, "\n%s\n", string); } return ret; } }
#define ABORT(format, ...)                   { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true);}
#define ABORT0(string)                       { fprintf (stderr, "\n%s\n", string); exit_on_error(true);}
#define WARN(format, ...)                    { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); }
#define WARN0(string)                        { fprintf (stderr, "\n%s\n", string); }
#define ASSERTGOTO(condition, format, ...)   { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); goto error; }}

#endif