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

// default max amount of VCF data in each variant block. this is user-configurable with --vblock
#define TXT_DATA_PER_VB_DEFAULT "16" // MB in default mode
#define TXT_DATA_PER_VB_FAST    "16" // MB with --fast

#define MAX_SUBFIELDS 100   // maximum number of VCF_FORMAT subfield types (except for GT), VCF_INFO, SAM_QNAME, SAM_OPTIONAL, GFF3_ATTRS subfield types that is supported in one line.
#define MAX_DICTS     253   // 254 is for future use and 255 is DID_I_NONE

#define DEFAULT_MAX_THREADS 8 // used if num_cores is not discoverable and the user didn't specifiy --threads

// ------------------------
// VCF stuff
// ------------------------

// default max number of samples in each sample block within a variant block. user configurable with --sblock
#define VCF_SAMPLES_PER_VBLOCK "4096" 

#define VCF_MAX_PLOIDY     100  // set to a reasonable 100 to avoid memory allocation explosion in case of an error in the VCF file
#if VCF_MAX_PLOIDY > 65535
#error "VCF_MAX_PLOIDY cannot go beyond 65535 as SectionHeaderVbHeaderVCF.ploidy and VBlockVCF.ploidy are uint16_t"
#endif

#define VCF_MAX_ALLELE_VALUE   99 // the code currently allows for 2-digit alleles.

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
typedef struct MtfContext *MtfContextP;
typedef struct MtfNode *MtfNodeP;
typedef const struct MtfNode *ConstMtfNodeP;
typedef struct SectionHeader *SectionHeaderP;
typedef struct SectionListEntry *SectionListEntryP;
typedef const struct SectionListEntry *ConstSectionListEntryP;

typedef enum { EXE_GENOZIP, EXE_GENOUNZIP, EXE_GENOLS, EXE_GENOCAT } ExeType;

#pragma pack(push, 1) // structures that are part of the genozip format are packed.
#define DICT_ID_LEN    ((int)sizeof(uint64_t))    // VCF/GFF3 specs don't limit the field name (tag) length, we limit it to 8 chars. zero-padded. (note: if two fields have the same 8-char prefix - they will just share the same dictionary)
typedef union DictIdType {
    uint64_t num;            // num is just for easy comparisons - it doesn't have a numeric value and endianity should not be changed
    uint8_t id[DICT_ID_LEN]; // \0-padded IDs 
    uint16_t map_key;        // we use the first two bytes as they key into vb/z_file->dict_id_mapper
} DictIdType;
#pragma pack(pop)

// global parameters - set before any thread is created, and never change
extern uint32_t global_max_threads, global_max_memory_per_vb;
extern const char *global_cmd;            // set once in main()
extern ExeType exe_type;

#define ZIP        'z'
#define UNZIP      'd'
#define LIST       'l'
#define LICENSE    'L'
#define VERSION    'V'
#define HELP       'h'
extern int command;

// flags set by user's command line options
extern int flag_force, flag_quiet, flag_concat, flag_md5, flag_split, flag_show_alleles, flag_show_time, flag_bgzip, flag_bam, flag_bcf,
           flag_show_memory, flag_show_dict, flag_show_gt_nodes, flag_show_b250, flag_show_sections, flag_show_headers,
           flag_show_index, flag_show_gheader, flag_stdout, flag_replace, flag_test, flag_regions,
           flag_samples, flag_drop_genotypes, flag_no_header, flag_header_only, flag_show_threads,
           flag_show_vblocks, flag_optimize, flag_gtshark, flag_sblock, flag_vblock, flag_gt_only,
           flag_header_one, flag_fast, flag_multiple_files, flag_fasta_sequential, flag_register,
           flag_debug_progress, flag_show_hash, flag_debug_memory, flag_debug_no_singletons,

           flag_optimize_sort, flag_optimize_PL, flag_optimize_GL, flag_optimize_GP, flag_optimize_VQSLOD, 
           flag_optimize_QUAL, flag_optimize_Vf, flag_optimize_ZM;
           
extern char *flag_grep;
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
    extern void func (VBlockP vb, MtfContextP ctx, const char *snip, unsigned snip_len); \
    static const int dt##_SPECIAL_##name = (num + 32); /* +32 to make it printable ASCII that can go into a snip */

#define COMPRESSOR_CALLBACK(func) \
extern void func (VBlockP vb, uint32_t vb_line_i, \
                  char **line_data_1, uint32_t *line_data_len_1,\
                  char **line_data_2, uint32_t *line_data_len_2);
                  
// sanity checks
extern void exit_on_error(void);

#define ASSERT(condition, format, ...)       { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(); }}
#define ASSERT0(condition, string)           { if (!(condition)) { fprintf (stderr, "\n%s\n", string); exit_on_error(); }}
#define ASSERTW(condition, format, ...)      { if (!(condition) && !flag_quiet) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); }}
#define ASSERTW0(condition, string)          { if (!(condition) && !flag_quiet) { fprintf (stderr, "\n%s\n", string); } }
#define RETURNW(condition, ret, format, ...) { if (!(condition)) { if (!flag_quiet) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); } return ret; }}
#define RETURNW0(condition, ret, string)     { if (!(condition)) { if (!flag_quiet) { fprintf (stderr, "\n%s\n", string); } return ret; } }
#define ABORT(format, ...)                   { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error();}
#define ABORT0(string)                       { fprintf (stderr, "\n%s\n", string); exit_on_error();}
#define ASSERTGOTO(condition, format, ...)   { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); goto error; }}

#endif