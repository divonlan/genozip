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
#include "compatability/visual_c_stdint.h"
#include "compatability/visual_c_stdbool.h"
#include "compatability/visual_c_unistd.h"
#include "compatability/visual_c_misc_funcs.h"
#endif

#define GENOZIP_EXT ".genozip"

// this was carefully picked as the optimal number based on testing with 1000-genomes chromosome 22 - 1024 samples:
// there is a tradeoff: we sort haplotypes by number of 1s for this entire variant block, and maintain and Index that is 
// #haplotypes x roundup(log2(#haplotypes)) - and we the number of indices we have in the final file depends on the number
// of variant blocks - the less the better. On the other hand - smaller variant blocks result in more linkage disequilibrium between
// variants in the block (assuming the VCF is sorted by POS), resulting in our sorting by number of 1s more likely to result
// in similar haplotypes grouped together - improving the compression.

#define SAMPLES_PER_BLOCK  4096 // tradeoff: larger is better compression, but in some cases might be slower retrieval speed
#define MAX_PLOIDY         100  // this can be any number up to 65535, it is set to 100 to avoid memory allocation
                                // explosion in case of an error in the VCF file

#define MAX_SUBFIELDS      64   // maximum number of FORMAT subfield types (except for GT) and INFO subfield types that is supported in one GENOZIP file. This constant can be increased if needed.
#define MAX_DICTS          (MAX_SUBFIELDS + MAX_SUBFIELDS + 8)   // dictionaries of subfields, infos and the 9 first fields (tabs) of the VCF file (+8 because REF and ALT are combined)

#define MAX_CHROM_LEN      64   // maximum length of chromosome (contig) name

// Note: the algorithm will use as many cores as it can - but there's no speed penalty for a higher MAX_COMPUTE_THREADS
// that number of cores - it will not be faster, but also not slower.
// However, each thread adds memory consumption approximately linearly

// the optimal number of compute threads is determined by the ratio between CPU time and I/O time. 
// For uncompress, 3 or more threads result in similar performance for HD and SSD, with 7 seaming to be about optimal (but not a big difference than 3). 
// Adding threads doesn't help. 2 or less threads result in significantly slower execution time. 
// Memory consumption is linear with the number of threads (each allocated a VB)

#define DEFAULT_MAX_THREADS 8 // maximum compute threads created - one I/O thread, and the rest of compute threads. This can be changed with command line option -@
 
#define MAX_32BIT_WINDOWS_MEMORY (1.7*1024*1024*1024) // 1.7GB - so Windows 32bit code doesn't explode at 2GB. TO DO - make this platform specific or get ulimits

// note: the numbering of the sections cannot be modified, for backward compatability
typedef enum {
    // data sections - statring in v1
    SEC_VCF_HEADER         = 0,  SEC_VB_HEADER           = 1, 
    SEC_GENOTYPE_DICT      = 2,  SEC_GENOTYPE_DATA       = 3, 
    SEC_PHASE_DATA         = 4,  SEC_HAPLOTYPE_DATA      = 5,

    // data sections added in v2
    SEC_GENOZIP_HEADER     = 6,

    SEC_CHROM_DICT         = 7,   SEC_CHROM_B250         = 8,
    SEC_POS_DICT           = 9,   SEC_POS_B250           = 10,
    SEC_ID_DICT            = 11,  SEC_ID_B250            = 12,  
    SEC_REFALT_DICT        = 13,  SEC_REFALT_B250        = 14,  
    SEC_QUAL_DICT          = 15,  SEC_QUAL_B250          = 16, 
    SEC_FILTER_DICT        = 17,  SEC_FILTER_B250        = 18, 
    SEC_INFO_DICT          = 19,  SEC_INFO_B250          = 20, 
    SEC_FORMAT_DICT        = 21,  SEC_FORMAT_B250        = 22,
    SEC_INFO_SUBFIELD_DICT = 23,  SEC_INFO_SUBFIELD_B250 = 24,
    
    // These sections are not real sections - they don't appear in the genozip file - just for stats. They can be changed if needed.
    SEC_STATS_HT_SEPERATOR, 
} SectionType;

// we put the names here in a #define so we can eyeball their identicality to SectionType
#define SECTIONTYPE_NAMES { \
    "SEC_VCF_HEADER"        ,  "SEC_VB_HEADER",\
    "SEC_GENOTYPE_DICT"     ,  "SEC_GENOTYPE_DATA",\
    "SEC_PHASE_DATA"        ,  "SEC_HAPLOTYPE_DATA",\
    \
    "SEC_GENOZIP_HEADER"    ,\
    \
    "SEC_CHROM_DICT"        ,  "SEC_CHROM_B250",\
    "SEC_POS_DICT"          ,  "SEC_POS_B250",\
    "SEC_ID_DICT"           ,  "SEC_ID_B250",\
    "SEC_REFALT_DICT"       ,  "SEC_REFALT_B250",\
    "SEC_QUAL_DICT"         ,  "SEC_QUAL_B250",\
    "SEC_FILTER_DICT"       ,  "SEC_FILTER_B250",\
    "SEC_INFO_DICT"         ,  "SEC_INFO_B250",\
    "SEC_FORMAT_DICT"       ,  "SEC_FORMAT_B250",\
    "SEC_INFO_SUBFIELD_DICT",  "SEC_INFO_SUBFIELD_B250",\
    \
    "SEC_STATS_HT_SEPERATOR" \
}

#define NUM_SEC_TYPES (SEC_STATS_HT_SEPERATOR+1) // put this here and not in sections.h as its used in vb.h that is widely used

#define section_type_is_dictionary(s) (((s) >= SEC_CHROM_DICT && (s) <= SEC_FORMAT_DICT && (s % 2 == 1)) ||       \
                                        (s) == SEC_INFO_SUBFIELD_DICT || (s) == SEC_GENOTYPE_DICT)

#define section_type_is_b250(s)       (((s) >= SEC_CHROM_B250 && (s) <= SEC_FORMAT_B250 && (s % 2 == 0)) ||       \
                                        (s) == SEC_INFO_SUBFIELD_B250)

// pointers used in header files - so we don't need to include the whole .h (and avoid cyclicity)
typedef struct variant_block_ *VariantBlockP;
typedef const struct variant_block_ *ConstVariantBlockP;
typedef struct file_ *FileP;
typedef const struct file_ *ConstFileP;

// global parameters - set before any thread is created, and never change
extern unsigned    global_num_samples, global_max_lines_per_vb;
extern const char *global_cmd;            // set once in main()

// flags set by user's command line options
extern int flag_force, flag_quiet, flag_concat_mode, flag_md5, flag_split, flag_show_alleles, flag_show_time, 
           flag_show_memory, flag_show_dict, flag_show_gt_nodes, flag_show_b250, flag_show_sections;

// macros
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b) )
#define MAX(a, b) (((a) > (b)) ? (a) : (b) )
#endif

// sanity checks
static inline void my_exit() { exit(1); }// an exit function so we can put a debugging break point when ASSERT exits
#define ASSERT(condition, format, ...)  { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); my_exit(); }}
#define ASSERT0(condition, string)      { if (!(condition)) { fprintf (stderr, "\n%s\n", string); my_exit(); }}
#define ASSERTW(condition, format, ...) { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); }}
#define ASSERTW0(condition, string)     { if (!(condition)) { fprintf (stderr, "\n%s\n", string); } }
#define ABORT(format, ...)              { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); my_exit();}
#define ABORT0(string)                  { fprintf (stderr, "\n%s\n", string); my_exit();}

#endif