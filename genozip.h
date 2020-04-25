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
#include <stdbool.h>
#include <unistd.h>
#else
#include "compatibility/visual_c_stdint.h"
#include "compatibility/visual_c_stdbool.h"
#include "compatibility/visual_c_unistd.h"
#include "compatibility/visual_c_misc_funcs.h"
#endif

// -----------------
// system parameters
// -----------------
#define GENOZIP_EXT ".genozip"

// default max amount of VCF data in each variant block. this is user-configurable with --vblock
#define TXT_DATA_PER_VB_DEFAULT "128" // MB in default mode
#define TXT_DATA_PER_VB_FAST    "16"  // MB with --fast

#define MAX_SUBFIELDS      63   // maximum number of VCF_FORMAT subfield types (except for GT), VCF_INFO, SAM_QNAME and SAM_OPTIONAL subfield types that is supported in one GENOZIP file.

#define MAX_DICTS          (MAX_SUBFIELDS + MAX_SUBFIELDS + 8)   // dictionaries of subfields, infos and the 9 first fields (tabs) of the VCF file (+8 because REF and ALT are combined). 
#if MAX_DICTS > 253 // 254 and 255 are used for special purpose
#error "MAX_DICTS cannot go beyond 255 as SubfieldMapper and SubfieldInfoMapperPiz represent did_i as uint8_t, and NIL=255"
#endif

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
typedef struct VBlockVCF *VBlockVCFP;
typedef const struct VBlockVCF *ConstVBlockVCFP;
typedef struct VBlockSAM *VBlockSAMP;
typedef const struct VBlockSAM *ConstVBlockSAMP;
typedef struct VBlockFAST *VBlockFASTP;
typedef const struct VBlockFAST *ConstVBlockFASTP;
typedef struct VBlockME23 *VBlockME23P;
typedef const struct VBlockME23 *ConstVBlockME23P;

typedef struct SubfieldMapper *SubfieldMapperP; 
typedef const struct SubfieldMapper *ConstSubfieldMapperP; 
typedef struct File *FileP;
typedef const struct File *ConstFileP;
typedef struct Buffer *BufferP;
typedef const struct Buffer *ConstBufferP;
typedef struct MtfContext *MtfContextP;
typedef struct MtfNode *MtfNodeP;
typedef struct SectionHeader *SectionHeaderP;
typedef struct SectionListEntry *SectionListEntryP;
typedef const struct SectionListEntry *ConstSectionListEntryP;

typedef enum { EXE_GENOZIP, EXE_GENOUNZIP, EXE_GENOLS, EXE_GENOCAT } ExeType;

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
           flag_show_index, flag_show_gheader, flag_stdout, flag_replace, flag_show_content, flag_test, flag_regions,
           flag_samples, flag_drop_genotypes, flag_no_header, flag_header_only, flag_show_threads, flag_debug_memory,
           flag_show_vblocks, flag_optimize, flag_gtshark, flag_sblock, flag_vblock, flag_strip, flag_gt_only,
           flag_header_one, flag_fast, flag_multiple_files;

extern uint64_t flag_stdin_size;

// external vb - used when an operation is needed outside of the context of a specific variant block;
extern VBlockP evb;

// macros
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b) )
#define MAX(a, b) (((a) > (b)) ? (a) : (b) )
#endif

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