// ------------------------------------------------------------------
//   sections.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SECTIONS_INCLUDED
#define SECTIONS_INCLUDED

#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#include <stdbool.h>
#else
#include "compatibility/visual_c_stdint.h"
#include "compatibility/visual_c_stdbool.h"
#endif

#include "md5.h"
#include "dict_id.h"

// note: the numbering of the sections cannot be modified, for backward compatibility
typedef enum {
    SEC_EOF                = -1, // doesn't appear in the file - just a value to indicate there are no more sections

    // data sections - statring in v1
    SEC_TXT_HEADER         = 0,                                // used by VCF, SAM, ME23
    SEC_VB_HEADER          = 1,                                // used by all data types

    SEC_VCF_FRMT_SF_DICT   = 2,  
    SEC_VCF_GT_DATA        = 3,  SEC_VCF_PHASE_DATA     = 4,
    SEC_HT_DATA            = 5,                                // used by VCF, ME23

    // data sections added in v2
    SEC_GENOZIP_HEADER     = 6,   SEC_RANDOM_ACCESS      = 7,

    SEC_CHROM_DICT         = 8,   SEC_CHROM_B250         = 9,  // used by VCF, ME23
    SEC_POS_DICT           = 10,  SEC_POS_B250           = 11, // used by VCF, ME23
    SEC_ID_DICT            = 12,  SEC_ID_B250            = 13, // used by VCF until v4 (VCF's ID moved to SEC_NUMERIC_ID_DATA in v5)
    SEC_VCF_REFALT_DICT    = 14,  SEC_VCF_REFALT_B250    = 15,  
    SEC_VCF_QUAL_DICT      = 16,  SEC_VCF_QUAL_B250      = 17, 
    SEC_VCF_FILTER_DICT    = 18,  SEC_VCF_FILTER_B250    = 19, 
    SEC_VCF_INFO_DICT      = 20,  SEC_VCF_INFO_B250      = 21, 
    SEC_VCF_FORMAT_DICT    = 22,  SEC_VCF_FORMAT_B250    = 23,
    SEC_VCF_INFO_SF_DICT   = 24,  SEC_VCF_INFO_SF_B250   = 25,
    
    // added gtshark support in file version 3
    SEC_HT_GTSHARK_DB_DB   = 26,  SEC_HT_GTSHARK_DB_GT   = 27,
    SEC_HT_GTSHARK_X_LINE  = 28,  SEC_HT_GTSHARK_X_HTI   = 29,
    SEC_HT_GTSHARK_X_ALLELE= 30,

    // added SAM, FASTQ, FASTA, 23andMe support in file version 5
    SEC_SAM_RAND_POS_DATA  = 31,  SEC_SAM_MD_DATA        = 32, 
    SEC_SAM_BD_DATA        = 33,  SEC_SAM_BI_DATA        = 34,
    SEC_SEQ_DATA           = 35,                               // used by SAM, FASTA, FASTQ 
    SEC_QUAL_DATA          = 36,                               // used by SAM, FASTQ
    SEC_NUMERIC_ID_DATA    = 37,                               // used by VCF (starting v5) and ME23
    SEC_SAM_QNAME_SF_DICT  = 38,  SEC_SAM_QNAME_SF_B250  = 39,
    SEC_SAM_OPTNL_SF_DICT  = 40,  SEC_SAM_OPTNL_SF_B250  = 41,
    SEC_SAM_QNAME_DICT     = 42,  SEC_SAM_QNAME_B250     = 43,
    SEC_SAM_FLAG_DICT      = 44,  SEC_SAM_FLAG_B250      = 45,
    SEC_SAM_RNAME_DICT     = 46,  SEC_SAM_RNAME_B250     = 47, 
    SEC_SAM_POS_DICT       = 48,  SEC_SAM_POS_B250       = 49, 
    SEC_SAM_MAPQ_DICT      = 50,  SEC_SAM_MAPQ_B250      = 51, 
    SEC_SAM_CIGAR_DICT     = 52,  SEC_SAM_CIGAR_B250     = 53, 
    SEC_SAM_PNEXT_DICT     = 54,  SEC_SAM_PNEXT_B250     = 55, 
    SEC_SAM_TLEN_DICT      = 56,  SEC_SAM_TLEN_B250      = 57, 
    SEC_SAM_OPTIONAL_DICT  = 58,  SEC_SAM_OPTIONAL_B250  = 59, 

    SEC_FAST_DESC_SF_DICT  = 60,  SEC_FAST_DESC_SF_B250  = 61, // used by FASTQ & FASTA
    SEC_FAST_DESC_DICT     = 62,  SEC_FAST_DESC_B250     = 63, // used by FASTQ & FASTA
    SEC_FAST_LINEMETA_DICT = 64,  SEC_FAST_LINEMETA_B250 = 65, // used by FASTQ & FASTA
    SEC_FASTA_COMMENT_DATA = 66,

    SEC_GFF3_SEQID_DICT    = 67,  SEC_GFF3_SEQID_B250    = 68, 
    SEC_GFF3_SOURCE_DICT   = 69,  SEC_GFF3_SOURCE_B250   = 70,
    SEC_GFF3_TYPE_DICT     = 71,  SEC_GFF3_TYPE_B250     = 72, 
    SEC_GFF3_START_DICT    = 73,  SEC_GFF3_START_B250    = 74,
    SEC_GFF3_END_DICT      = 75,  SEC_GFF3_END_B250      = 76, 
    SEC_GFF3_SCORE_DICT    = 77,  SEC_GFF3_SCORE_B250    = 78,
    SEC_GFF3_STRAND_DICT   = 79,  SEC_GFF3_STRAND_B250   = 80, 
    SEC_GFF3_PHASE_DICT    = 81,  SEC_GFF3_PHASE_B250    = 82,
    SEC_GFF3_ATTRS_DICT    = 83,  SEC_GFF3_ATTRS_B250    = 84, 
    SEC_GFF3_ATTRS_SF_DICT = 85,  SEC_GFF3_ATTRS_SF_B250 = 86, 

    SEC_ENST_DATA          = 87,
    // This sections is not a real section - it doesn't appear in the genozip file. It can be changed if needed.
    SEC_STATS_HT_SEPERATOR
} SectionType;

// this data must be perfectly aligned with SectionType. it contains: 1. name 2. 1 if the section is stripped out in --strip 
#define SECTIONTYPE_ABOUT { \
    {"SEC_TXT_HEADER",          0},  {"SEC_VB_HEADER",          0},\
    \
    {"SEC_VCF_FRMT_SF_DICT",    1},  {"SEC_VCF_GT_DATA",        1},\
    {"SEC_VCF_PHASE_DATA",      0},  {"SEC_HT_DATA",            0},\
    \
    {"SEC_GENOZIP_HEADER",      0},  {"SEC_RANDOM_ACCESS",      0},\
    \
    {"SEC_CHROM_DICT",          0},  {"SEC_CHROM_B250",         0},\
    {"SEC_POS_DICT",            0},  {"SEC_POS_B250",           0},\
    {"SEC_ID_DICT",             1},  {"SEC_ID_B250",            1},\
    {"SEC_VCF_REFALT_DICT",     0},  {"SEC_VCF_REFALT_B250",    0},\
    {"SEC_VCF_QUAL_DICT",       1},  {"SEC_VCF_QUAL_B250",      1},\
    {"SEC_VCF_FILTER_DICT",     1},  {"SEC_VCF_FILTER_B250",    1},\
    {"SEC_VCF_INFO_DICT",       1},  {"SEC_VCF_INFO_B250",      1},\
    {"SEC_VCF_FORMAT_DICT",     1},  {"SEC_VCF_FORMAT_B250",    1},\
    {"SEC_VCF_INFO_SF_DICT",    1},  {"SEC_VCF_INFO_SF_B250",   1},\
    \
    {"SEC_HT_GTSHARK_DB_DB",    0},  {"SEC_HT_GTSHARK_DB_GT",   0},\
    {"SEC_HT_GTSHARK_X_LINE",   0},  {"SEC_HT_GTSHARK_X_HTI",   0},\
    {"SEC_HT_GTSHARK_X_ALLELE", 0},\
    \
    {"SEC_SAM_RAND_POS_DATA",   0},  {"SEC_SAM_MD_DATA",        1},\
    {"SEC_SAM_BD_DATA",         1},  {"SEC_SAM_BI_DATA",        1},\
    {"SEC_SEQ_DATA",            0},\
    {"SEC_QUAL_DATA",           1},\
    {"SEC_NUMERIC_ID_DATA",     1},\
    {"SEC_SAM_QNAME_SF_DICT",   1},  {"SEC_SAM_QNAME_SF_B250",  1},\
    {"SEC_SAM_OPTNL_SF_DICT",   0},  {"SEC_SAM_OPTNL_SF_B250",  0},\
    {"SEC_SAM_QNAME_DICT",      1},  {"SEC_SAM_QNAME_B250",     1},\
    {"SEC_SAM_FLAG_DICT",       1},  {"SEC_SAM_FLAG_B250",      1},\
    {"SEC_SAM_RNAME_DICT",      0},  {"SEC_SAM_RNAME_B250",     0},\
    {"SEC_SAM_POS_DICT",        0},  {"SEC_SAM_POS_B250",       0},\
    {"SEC_SAM_MAPQ_DICT",       1},  {"SEC_SAM_MAPQ_B250",      1},\
    {"SEC_SAM_CIGAR_DICT",      0},  {"SEC_SAM_CIGAR_B250",     0},\
    {"SEC_SAM_PNEXT_DICT",      1},  {"SEC_SAM_PNEXT_B250",     1},\
    {"SEC_SAM_TLEN_DICT",       1},  {"SEC_SAM_TLEN_B250",      1},\
    {"SEC_SAM_OPTIONAL_DICT",   0},  {"SEC_SAM_OPTIONAL_B250",  0},\
    \
    {"SEC_FAST_DESC_SF_DICT",   1},  {"SEC_FAST_DESC_SF_B250",  1},\
    {"SEC_FAST_DESC_DICT",      1},  {"SEC_FAST_DESC_B250",     1},\
    {"SEC_FAST_LINEMETA_DICT",  0},  {"SEC_FAST_LINEMETA_B250", 0},\
    {"SEC_FASTA_COMMENT_DATA",  1},\
    \
    {"SEC_GFF3_SEQID_DICT",     0}, {"SEC_GFF3_SEQID_B250",     0},\
    {"SEC_GFF3_SOURCE_DICT",    0}, {"SEC_GFF3_SOURCE_B250",    0},\
    {"SEC_GFF3_TYPE_DICT",      0}, {"SEC_GFF3_TYPE_B250",      0},\
    {"SEC_GFF3_START_DICT",     0}, {"SEC_GFF3_START_B250",     0},\
    {"SEC_GFF3_END_DICT",       0}, {"SEC_GFF3_END_B250",       0},\
    {"SEC_GFF3_SCORE_DICT",     0}, {"SEC_GFF3_SCORE_B250",     0},\
    {"SEC_GFF3_STRAND_DICT",    0}, {"SEC_GFF3_STRAND_B250",    0},\
    {"SEC_GFF3_PHASE_DICT",     0}, {"SEC_GFF3_PHASE_B250",     0},\
    {"SEC_GFF3_ATTRS_DICT",     0}, {"SEC_GFF3_ATTRS_B250",     0},\
    {"SEC_GFF3_ATTRS_SF_DICT",  0}, {"SEC_GFF3_ATTRS_SF_B250",  0},\
    \
    {"SEC_ENST_DATA",           0},\
    \
    {"SEC_STATS_HT_SEPERATOR",  0} \
}

#define NUM_SEC_TYPES (SEC_STATS_HT_SEPERATOR+1) 

#define section_type_is_dictionary(s) (((s) >= SEC_CHROM_DICT        && (s) <= SEC_VCF_INFO_SF_DICT   && ((s)%2) == (SEC_CHROM_DICT % 2)) ||       \
                                        (s) == SEC_VCF_FRMT_SF_DICT || \
                                       ((s) >= SEC_SAM_QNAME_SF_DICT && (s) <= SEC_FAST_LINEMETA_DICT && ((s)%2) == (SEC_SAM_QNAME_SF_DICT % 2)) ||\
                                       ((s) >= SEC_GFF3_SEQID_DICT   && (s) <= SEC_GFF3_ATTRS_SF_DICT && ((s)%2) == (SEC_GFF3_SEQID_DICT % 2)))

#define section_type_is_b250(s)       (section_type_is_dictionary((s)-1) && (s) != SEC_VCF_GT_DATA /* one exception */)

extern const SectionType first_field_dict_section[NUM_DATATYPES];

#define FIELD_TO_DICT_SECTION(dt,f) (first_field_dict_section[dt] + (f)*2)
#define FIELD_TO_B250_SECTION(dt,f) (FIELD_TO_DICT_SECTION((dt), (f)) + 1)

// Section headers - big endian

#define GENOZIP_MAGIC 0x27052012

// encryption types
#define NUM_ENCRYPTION_TYPES   2
#define ENCRYPTION_TYPE_NONE   0
#define ENCRYPTION_TYPE_AES256 1
#define ENCRYPTION_TYPE_NAMES { "No encryption", "AES 256 bit" }

#pragma pack(push, 1) // structures that are part of the genozip format are packed.

// section headers are encoded in Big Endian (see https://en.wikipedia.org/wiki/Endianness)
// the reason for selecting big endian is that I am developing on little endian CPU (Intel) so
// endianity bugs will be discovered more readily this way

typedef struct SectionHeader {
    uint32_t magic; 
    uint32_t compressed_offset;      // number of bytes from the start of the header that is the start of compressed data (sizeof header + header encryption padding)
    uint32_t data_encrypted_len;     // = data_compressed_len + padding if encrypted, 0 if not
    uint32_t data_compressed_len;
    uint32_t data_uncompressed_len;
    uint32_t vblock_i;               // VB with in file starting from 1 ; 0 for Txt Header
    uint16_t section_i;              // section within VB - 0 for Variant Data
    uint8_t  section_type;          
    uint8_t  sec_compression_alg : 4; // one of CompressionAlg. introduced in genozip v5 (before that it this field was unused)
    uint8_t  unused              : 4;
} SectionHeader; 

typedef struct {
    SectionHeader h;
    uint8_t  genozip_version;
    uint8_t  encryption_type;         // one of ENC_TYPE_*
    uint16_t data_type;               // one of DATA_TYPE_*
    uint32_t num_samples;             // number of samples. "samples" is data_type-dependent. 
    uint64_t uncompressed_data_size;  // data size of uncompressed file, if uncompressed as a single file
    uint64_t num_items_concat;        // number of items in a concatenated file. "item" is data_type-dependent. For VCF, it is lines.
    uint32_t num_sections;            // number sections in this file (including this one)
    uint32_t num_components;          // number of vcf concatenated components in this file (1 if no concatenation)

    Md5Hash  md5_hash_concat;         // md5 of original VCF file, or 0s if no hash was calculated. if this is a concatenation - this is the md5 of the entire concatenation.

    uint8_t  password_test[16];       // short encrypted block - used to test the validy of a password
#define FILE_METADATA_LEN 72
    char created[FILE_METADATA_LEN];    
    Md5Hash  license_hash;            // added in genozip v5. MD5(license_num)
} SectionHeaderGenozipHeader;

// this footer appears AFTER the genozip header data, facilitating reading the genozip header in reverse from the end of the file
typedef struct {
    uint64_t genozip_header_offset;
    uint32_t magic;
} SectionFooterGenozipHeader;

// The text file header section appears once in the file (or multiple times in case of concatenation), and includes the VCF file header 
typedef struct {
    SectionHeader h;
    uint64_t txt_data_size;    // number of bytes in the original VCF file. 
#define NUM_LINES_UNKNOWN ((uint64_t)-1) 
    uint64_t num_lines;        // number of data (non-header) lines in the original txt file. Concat mode: entire file for first SectionHeaderTxtHeader, and only for that txt if not first
    uint32_t num_samples;      // VCF only: number of samples in the original VCF file
    uint32_t max_lines_per_vb; // upper bound on how many data lines a VB can have in this file
    uint8_t  compression_type; // compression type of original file, one of CompressionAlg 
    Md5Hash  md5_hash_single;  // non-0 only if this genozip file is a result of concatenatation with --md5. md5 of original single txt file.

#define TXT_FILENAME_LEN 256
    char txt_filename[TXT_FILENAME_LEN]; // filename of this single component. without path, 0-terminated. always a .vcf or .sam, even if the original was eg .vcf.gz or .bam

} SectionHeaderTxtHeader; 

// A generic VB header - it will suffice for all data types but VCF. introduced in v5.
typedef struct {
    SectionHeader h;
    uint32_t first_line;               // line (starting from 1) of this vblock in the single VCF file
                                       // new in v2: if this value is 0, then this is the terminating section of the file. after it is either EOF or a VCF Header section of the next concatenated file
    uint32_t num_lines;                // number of records in this vblock
    
    // features of the data
    uint32_t vb_data_size;             // size of vblock as it appears in the source file
    uint32_t z_data_bytes;             // total bytes of this vblock in the genozip file including all sections and their headers
    uint32_t longest_line_len;         // length of the longest line in this vblock

    uint32_t ffu;                      // for future use
} SectionHeaderVbHeader; 

typedef struct {
    SectionHeader h;                   // in v1, the section_type was SEC_DICTIONARY, in v2 is is the specific SEC_*_DICT
    uint32_t num_snips;                // number of items in dictionary
    DictIdType dict_id;           
} SectionHeaderDictionary; 

// b250 data encoded according to a dictionary
typedef struct {
    SectionHeader h;
    uint32_t num_b250_items;           // number of items in b250 items
    DictIdType dict_id;           
} SectionHeaderBase250;         

// the data of SEC_SECTION_LIST is an array of the following type, as is the z_file->section_list_buf
typedef struct SectionListEntry {
    uint64_t offset;                   // offset of this section in the file
    DictIdType dict_id;                // used if this section is a DICT or a B250 section
    uint32_t vblock_i;
    uint8_t section_type;
    uint8_t unused[3];                 // padding
} SectionListEntry;


//----------------------------------------
// VCF Sections
//----------------------------------------

// The is the first section in every VCF VBlock
typedef struct {
    SectionHeader h;
    uint32_t first_line;               // line (starting from 1) of this variant block in the single VCF file
                                       // new in v2: if this value is 0, then this is the terminating section of the file. after it is either EOF or a VCF Header section of the next concatenated file
    uint32_t num_lines;                // number of variants in this block
    uint8_t phase_type;
    
    // flags
    uint8_t has_genotype_data : 1;     // 1 if there is at least one variant in the block that has FORMAT with have anything except for GT 
    uint8_t is_gtshark        : 1;     // 1 if the haplotype sections are compressed with gtshark instead of bzip2
    uint8_t for_future_use    : 6;
    uint16_t ploidy;

    // features of the data
    uint32_t num_samples;
    uint32_t num_haplotypes_per_line;  // 0 if no haplotypes
    uint32_t num_sample_blocks;
    uint32_t num_samples_per_block;
    uint32_t max_gt_line_len;

    uint32_t vb_data_size;             // size of variant block as it appears in the source file
    uint32_t z_data_bytes;             // total bytes of this variant block in the genozip file including all sections and their headers
    uint16_t haplotype_index_checksum;
    uint16_t unused3;                  // new in v2: padding / ffu
} SectionHeaderVbHeaderVCF; 

// the data of SEC_RANDOM_ACCESS is an array of the following type, as is the z_file->ra_buf and vb->ra_buf
// we maintain one RA entry per vb per every chrom in the the VB
typedef struct {
    uint32_t vblock_i;                 // the vb_i in which this range appears
    uint32_t chrom_index;              // before merge: node index into chrom context mtf, after merge - word index in CHROM dictionary
    uint32_t min_pos, max_pos;         // POS field value of smallest and largest POS value of this chrom in this VB (regardless of whether the VB is sorted)
} RAEntry; 

// ------------------------------------------------------------------------------------------------------
// GENOZIP_FILE_FORMAT_VERSION==1 historical version (VCF only) - we support uncomrpessing old version files

typedef struct {
    SectionHeader h;
    uint8_t  genozip_version;
    uint32_t num_samples;    // number of samples in the original VCF file
    uint64_t txt_data_size;  // number of bytes in the original VCF file. Concat mode: entire file for first SectionHeaderTxtHeader, and only for that VCF if not first
    uint64_t num_lines;      // number of variants (data lines) in the original VCF file. Concat mode: entire file for first SectionHeaderTxtHeader, and only for that VCF if not first
    Md5Hash md5_hash_concat; // md5 of original VCF file, or 0s if no hash was calculated. if this is a concatenation - this is the md5 of the entire concatenation.
    Md5Hash md5_hash_single; // non-0 only if this genozip file is a result of concatenatation with --md5. md5 of original single VCF file.

#define v1_VCF_FILENAME_LEN 256
    char txt_filename[v1_VCF_FILENAME_LEN];    // filename of this single component. without path, 0-terminated.
    
#define v1_FILE_METADATA_LEN 72
    char created[v1_FILE_METADATA_LEN];    
} v1_SectionHeaderVCFHeader; 

typedef struct {
    SectionHeader h;
    uint32_t first_line;               // line (starting from 1) of this variant block in the VCF file
    uint32_t num_lines;                // number of variants in this block
    uint8_t phase_type;
    
    // flags
    uint8_t has_genotype_data : 1;     // 1 if there is at least one variant in the block that has FORMAT with have anything except for GT 
    uint8_t is_sorted_by_pos  : 1;     // 1 if variant block is sorted by POS
    uint8_t for_future_use    : 6;

    // features of the data
    uint32_t num_samples;
    uint32_t num_haplotypes_per_line;  // 0 if no haplotypes
    uint32_t num_sample_blocks;
    uint32_t num_samples_per_block;
    uint16_t ploidy;
    uint16_t num_dict_ids;            // number of subfields used in this VB. the actual fields are deduced from the FORMAT column on each line
    uint16_t num_dictionary_sections;  // we have a dictionary for each dict_id, in which new words were introduced in this VB that weren't already in the dictionary. so num_dictionary_sections <= num_dict_ids
    uint32_t max_gt_line_len;
#define v1_MAX_CHROM_LEN      64       // maximum length of chromosome (contig) name
    char chrom[v1_MAX_CHROM_LEN];      // a null-terminated ID of the chromosome
    int64_t min_pos, max_pos;          // minimum and maximum POS values in this VB. -1 if unknown. Note: our format support 64bit POS, but VCF specification as well as the POS dictionary supports 32bit (values 0 to 2M-1)
    uint32_t vb_data_size;             // size of variant block as it appears in the source file
    uint32_t z_data_bytes;             // total bytes of this variant block in the genozip file including all sections and their headers
    uint16_t haplotype_index_checksum;
    uint8_t haplotype_index[];         // length is num_haplotypes. e.g. the first entry shows for the first haplotype in the original file, its index into the permuted block. # of bits per entry is roundup(log2(num_samples*ploidy)).
} v1_SectionHeaderVariantData; 

typedef struct {
    SectionHeader h;
    uint8_t  genozip_version;
    uint8_t  encryption_type;         // one of ENC_TYPE_*
    uint16_t data_type;               // one of DATA_TYPE_*
    uint32_t num_samples;             // number of samples. "samples" is data_type-dependent. 
    uint64_t uncompressed_data_size;  // data size of uncompressed file, if uncompressed as a single file
    uint64_t num_items_concat;        // number of items in a concatenated file. "item" is data_type-dependent. For VCF, it is lines.
    uint32_t num_sections;            // number sections in this file (including this one)
    uint32_t num_components;          // number of vcf concatenated components in this file (1 if no concatenation)
    Md5Hash  md5_hash_concat;         // md5 of original VCF file, or 0s if no hash was calculated. if this is a concatenation - this is the md5 of the entire concatenation.
    uint8_t  password_test[16];       // short encrypted block - used to test the validy of a password
#define v2v3v4_FILE_METADATA_LEN 72
    char created[v2v3v4_FILE_METADATA_LEN];    
} v2v3v4_SectionHeaderGenozipHeader;

#pragma pack(pop)

// zip stuff
extern void sections_add_to_list (VBlockP vb, const SectionHeader *header);
extern void sections_list_concat (VBlockP vb, BufferP section_list_buf);

// piz stuff
extern uint8_t sections_count_sec_type (unsigned vb_i, SectionType sec);
extern SectionType sections_get_next_header_type(SectionListEntry **sl_ent, bool *skipped_vb, BufferP region_ra_intersection_matrix);
extern bool sections_get_next_dictionary(SectionListEntry **sl_ent);
extern bool sections_has_more_components(void);
extern SectionListEntry *sections_get_offset_first_section_of_type (SectionType st);
extern SectionListEntry *sections_vb_first (uint32_t vb_i);
//extern uint64_t sections_vb_next (SectionListEntry **sl_ent /* in / out */);

extern void BGEN_sections_list(void);
extern const char *st_name (SectionType sec_type);
extern void sections_show_gheader (SectionHeaderGenozipHeader *header);
extern void sections_get_sizes (DictIdType dict_id, uint32_t *dict_compressed_size, uint32_t *b250_compressed_size);

#endif
