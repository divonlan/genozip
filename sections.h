// ------------------------------------------------------------------
//   sections.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SECTIONS_INCLUDED
#define SECTIONS_INCLUDED

#include "genozip.h"
#include "md5.h"

// note: the numbering of the sections cannot be modified, for backward compatibility
typedef enum __attribute__ ((__packed__)) { // 1 byte
    SEC_NONE            = -1, // doesn't appear in the file 

    SEC_RANDOM_ACCESS   = 0,
    SEC_DICT_ID_ALIASES = 1,
    SEC_REFERENCE       = 2,
    SEC_REF_IS_SET      = 3,
    SEC_REF_HASH        = 4,
    SEC_REF_RANDOM_ACC  = 5,
    SEC_GENOZIP_HEADER  = 6, // SEC_GENOZIP_HEADER remains 6 as in v2-v5, to be able to read old versions' genozip header
    SEC_TXT_HEADER      = 7, 
    SEC_VB_HEADER       = 8,
    SEC_DICT            = 9, 
    SEC_B250            = 10, 
    SEC_LOCAL           = 11, 

    // vcf specific    
    SEC_VCF_GT_DATA     = 20,  
    SEC_VCF_PHASE_DATA  = 21,
    SEC_VCF_HT_DATA     = 22,                               
    SEC_VCF_HT_GTSHARK  = 23,

    NUM_SEC_TYPES // fake section for counting
} SectionType;

// this data must be perfectly aligned with SectionType.
#define SECTIONTYPE_ABOUT {  \
    {"SEC_RANDOM_ACCESS"},   \
    {"SEC_DICT_ID_ALIASES"}, \
    {"SEC_REFERENCE"},       \
    {"SEC_REF_IS_SET"},      \
    {"SEC_REF_HASH"},        \
    {"SEC_REF_RANDOM_ACC"},  \
    {"SEC_GENOZIP_HEADER"},  \
    {"SEC_TXT_HEADER"},      \
    {"SEC_VB_HEADER"},       \
    {"SEC_DICT"},            \
    {"SEC_B250"},            \
    {"SEC_LOCAL"},           \
    {}, {}, {}, {}, {}, {}, {}, {},  \
    {"SEC_VCF_GT_DATA"},     \
    {"SEC_VCF_PHASE_DATA"},  \
    {"SEC_VCF_HT_DATA"},     \
    {"SEC_VCF_HT_GTSHARK"},  \
}

// Section headers - big endian

#define GENOZIP_MAGIC 0x27052012

#pragma pack(1) // structures that are part of the genozip format are packed.

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
    SectionType section_type;        // 1 byte
    CompressionAlg sec_compression_alg; // 1 byte
    uint8_t  flags;                  // CTX_FL_*
    uint8_t  ffu;
} SectionHeader; 

typedef struct {
    SectionHeader h;
    uint8_t  genozip_version;
    uint8_t  encryption_type;         // one of ENC_TYPE_*
    uint16_t data_type;               // one of DATA_TYPE_*
    uint32_t num_samples;             // number of samples. "samples" is data_type-dependent. 
    uint64_t uncompressed_data_size;  // data size of uncompressed` file, if uncompressed as a single file
    uint64_t num_items_bound;         // number of items in a bound file. "item" is data_type-dependent. For VCF, it is lines.
    uint32_t num_sections;            // number sections in this file (including this one)
    uint32_t num_components;          // number of txt bound components in this file (1 if no binding)

    Md5Hash  md5_hash_bound;          // md5 of original txt file, or 0s if no hash was calculated. if this is a binding - this is the md5 of the entire bound file.

    uint8_t  password_test[16];       // short encrypted block - used to test the validy of a password
#define FILE_METADATA_LEN 72
    char created[FILE_METADATA_LEN];    
    Md5Hash  license_hash;            // MD5(license_num)
#define REF_FILENAME_LEN 256
    char ref_filename[REF_FILENAME_LEN]; // external reference filename, null-terimated. ref_filename[0]=0 if there is no external reference.
    Md5Hash ref_file_md5;             // SectionHeaderGenozipHeader.md5_hash_bound of the reference FASTA genozip file
} SectionHeaderGenozipHeader;

// this footer appears AFTER the genozip header data, facilitating reading the genozip header in reverse from the end of the file
typedef struct {
    uint64_t genozip_header_offset;
    uint32_t magic;
} SectionFooterGenozipHeader;

// The text file header section appears once in the file (or multiple times in case of bound file), and includes the VCF file header 
typedef struct {
    SectionHeader h;
    uint64_t txt_data_size;    // number of bytes in the original txt file. 
#define NUM_LINES_UNKNOWN ((uint64_t)-1) 
    uint64_t num_lines;        // number of data (non-header) lines in the original txt file. Concat mode: entire file for first SectionHeaderTxtHeader, and only for that txt if not first
    uint32_t num_samples;      // VCF only: number of samples in the original VCF file
    uint32_t max_lines_per_vb; // upper bound on how many data lines a VB can have in this file
    CompressionAlg compression_type; // compression type of original file
    Md5Hash  md5_hash_single;  // non-0 only if this genozip file is a result of binding with --md5. md5 of original single txt file.

#define TXT_FILENAME_LEN 256
    char txt_filename[TXT_FILENAME_LEN]; // filename of this single component. without path, 0-terminated. always a .vcf or .sam, even if the original was eg .vcf.gz or .bam

} SectionHeaderTxtHeader; 

#define VB_HEADER_COMMON_FIELDS \
    SectionHeader h;            \
    uint32_t first_line;       /* line (starting from 1) of this vblock in the single VCF file */ \
                               /* if this value is 0, then this is the terminating section of the file. after it is either EOF or a VCF Header section of the next bound file */ \
    uint32_t num_lines;        /* number of records in this vblock */ \
                                \
    /* features of the data */  \
    uint32_t vb_data_size;     /* size of vblock as it appears in the source file */ \
    uint32_t z_data_bytes;     /* total bytes of this vblock in the genozip file including all sections and their headers */ \
    uint32_t longest_line_len; /* length of the longest line in this vblock */ \
                                \
    Md5Hash md5_hash_so_far;   /* partial calculation of MD5 up to and including this VB */

// A generic VB header - used for all data types but VCF
typedef struct {
    VB_HEADER_COMMON_FIELDS
} SectionHeaderVbHeader; 

// VCF VB header
typedef struct {
    VB_HEADER_COMMON_FIELDS

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

    uint16_t haplotype_index_checksum;
    uint16_t unused;
} SectionHeaderVbHeaderVCF; 

typedef struct {
    SectionHeader h;           // we use h.flags to store the first 4 bits of ctx->flags, also stored in SectionHeaderCtx. This is because VCF FORMAT contexts have no SectionHeaderCtx.
    uint32_t num_snips;        // number of items in dictionary
    DictId dict_id;           
} SectionHeaderDictionary; 

// LT_* values are consistent with BAM optional 'B' types (and extend them)
// IMPORTANT - if adding or chaging LTs, all lt_* arrays in sections.c need to be updated
typedef enum __attribute__ ((__packed__)) { // 1 byte
    LT_TEXT     = 0,
    LT_INT8     = 1,    
    LT_UINT8    = 2,
    LT_INT16    = 3,
    LT_UINT16   = 4,
    LT_INT32    = 5,
    LT_UINT32   = 6,
    LT_INT64    = 7,   // ffu
    LT_UINT64   = 8,   // ffu
    LT_FLOAT32  = 9,   // ffu
    LT_FLOAT64  = 10,  // ffu
    LT_SEQUENCE = 11,  // length of data extracted is determined by vb->seq_len
    LT_BITMAP   = 12,  // a bitmap
    LT_DOMQUAL  = 13,  // quality scores in domqual format
    NUM_LOCAL_TYPES
} LocalType;

// 3 letter names for --show-headers
#define LOCALTYPE_NAMES { "TXT", "I8 ", "U8 ", "I16", "U16", "I32", "U32", "I64", "U64", "F32", "F64", "SEQ", "BMP", "DOM" }
extern const char *lt_names[NUM_LOCAL_TYPES];
extern const char lt_to_sam_map[NUM_LOCAL_TYPES];
extern const int lt_sizeof_one[NUM_LOCAL_TYPES];
extern const bool lt_is_signed[NUM_LOCAL_TYPES];
extern const int64_t lt_min[NUM_LOCAL_TYPES], lt_max[NUM_LOCAL_TYPES];

typedef struct {
    SectionHeader h;
    LocalType ltype; // used by SEC_LOCAL
    uint8_t param;   // goes into ctx.b250/local.param if flags contains CTX_FL_COPY_PARAM
    uint8_t ffu[2];
    DictId dict_id;           
} SectionHeaderCtx;         

// two ways of storing a range:
// uncompacted - we will have one section, SEC_REFERENCE, containing the data, and first/last_pos containing the coordinates of this range
// compacting we will have 2 sections:
// - SEC_REF_IS_SET - containing a bitmap (1 bit per base), and chrom,first/last_pos containing the coordinates of this range
// - SEC_REFERENCE - containing the compacted range (i.e. excluding bases that have a "0" in the bitmap), 
//   with chrom,first_pos same as SEC_REF_IS_SET and last_pos set so that (last_pos-first_pos+1) = number of '1' bits in SEC_REF_IS_SET
// SEC_REFERENCE (in both cases) contains 2 bits per base, and SEC_REF_IS_SET contains 1 bit per location.
typedef struct {
    SectionHeader h;
    uint64_t pos;              // first pos within chrom (1-based) of this range         
    uint64_t gpos;             // first pos within genome (0-based) of this range
    uint32_t num_bases;        // number of bases (nucleotides) in this range
    uint32_t chrom_word_index; // index in context->word_list of the chrom of this reference range    
} SectionHeaderReference;

typedef struct {
    SectionHeader h;
    uint8_t num_layers;        // total number of layers
    uint8_t layer_i;           // layer number of this section (0=base layer, with the most bits)
    uint8_t layer_bits;        // number of bits in layer
    uint8_t ffu;
    uint32_t start_in_layer;   // start index within layer
} SectionHeaderRefHash;

// the data of SEC_SECTION_LIST is an array of the following type, as is the z_file->section_list_buf
typedef struct SectionListEntry {
    uint64_t offset;           // offset of this section in the file
    DictId dict_id;            // used if this section is a DICT, LOCAL or a B250 section
    uint32_t vblock_i;
    SectionType section_type;  // 1 byte
    uint8_t unused[3];         
} SectionListEntry;

// the data of SEC_RANDOM_ACCESS is an array of the following type, as is the z_file->ra_buf and vb->ra_buf
// we maintain one RA entry per vb per every chrom in the the VB
typedef struct RAEntry {
    uint32_t vblock_i;            // the vb_i in which this range appears
    WordIndex chrom_index;        // before merge: node index into chrom context mtf, after merge - word index in CHROM dictionary
    uint64_t min_pos, max_pos;    // POS field value of smallest and largest POS value of this chrom in this VB (regardless of whether the VB is sorted)
} RAEntry; 

#pragma pack()

// zip stuff
extern void sections_add_to_list (VBlockP vb, const SectionHeader *header);
extern void sections_list_concat (VBlockP vb, BufferP section_list_buf);

// piz stuff
extern SectionType sections_get_next_header_type (SectionListEntry **sl_ent, bool *skipped_vb, BufferP region_ra_intersection_matrix);
extern bool sections_get_next_section_of_type (SectionListEntry **sl_ent, uint32_t *cursor, SectionType st1, SectionType st2);
extern SectionType sections_peek (uint32_t cursor);
extern uint32_t sections_count_sections (SectionType st);

extern bool sections_has_more_components(void);
extern SectionListEntry *sections_get_offset_first_section_of_type (SectionType st);
extern SectionListEntry *sections_vb_first (uint32_t vb_i, bool soft_fail);
extern bool sections_seek_to (SectionType st, bool first);
extern void sections_get_prev_file_vb_i (const SectionListEntry *sl, uint32_t *prev_file_first_vb_i, uint32_t *prev_file_last_vb_i);

extern void BGEN_sections_list(void);
extern const char *st_name (SectionType sec_type);
extern void sections_show_gheader (SectionHeaderGenozipHeader *header);

extern bool sections_has_reference(void);
extern void sections_get_refhash_details (uint32_t *num_layers, uint32_t *base_layer_bits);
extern int64_t sections_get_genome_size (void);

#endif
