// ------------------------------------------------------------------
//   sections.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SECTIONS_INCLUDED
#define SECTIONS_INCLUDED

#include "genozip.h"
#include "section_types.h"
#include "md5.h"

// Section headers - big endian

#define GENOZIP_MAGIC 0x27052012

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
    uint8_t  sec_compression_alg : 4; // one of CompressionAlg
    uint8_t  flags               : 4; // section-type specific flags SEC_FLAG_*
} SectionHeader; 

#define SEC_FLAG_GENOZIP_HEADER_IS_REFERENCE 0x01
typedef struct {
    SectionHeader h;
    uint8_t  genozip_version;
    uint8_t  encryption_type;         // one of ENC_TYPE_*
    uint16_t data_type;               // one of DATA_TYPE_*
    uint32_t num_samples;             // number of samples. "samples" is data_type-dependent. 
    uint64_t uncompressed_data_size;  // data size of uncompressed file, if uncompressed as a single file
    uint64_t num_items_concat;        // number of items in a concatenated file. "item" is data_type-dependent. For VCF, it is lines.
    uint32_t num_sections;            // number sections in this file (including this one)
    uint32_t num_components;          // number of txt concatenated components in this file (1 if no concatenation)

    Md5Hash  md5_hash_concat;         // md5 of original txt file, or 0s if no hash was calculated. if this is a concatenation - this is the md5 of the entire concatenation.

    uint8_t  password_test[16];       // short encrypted block - used to test the validy of a password
#define FILE_METADATA_LEN 72
    char created[FILE_METADATA_LEN];    
    Md5Hash  license_hash;            // MD5(license_num)
#define REF_FILENAME_LEN 255
    char ref_filename[REF_FILENAME_LEN]; // external reference filename, null-terimated. ref_filename[0]=0 if there is no external reference.
    Md5Hash ref_file_md5;             // SectionHeaderGenozipHeader.md5_hash_concat of the reference FASTA genozip file
} SectionHeaderGenozipHeader;

// this footer appears AFTER the genozip header data, facilitating reading the genozip header in reverse from the end of the file
typedef struct {
    uint64_t genozip_header_offset;
    uint32_t magic;
} SectionFooterGenozipHeader;

// The text file header section appears once in the file (or multiple times in case of concatenation), and includes the VCF file header 
typedef struct {
    SectionHeader h;
    uint64_t txt_data_size;    // number of bytes in the original txt file. 
#define NUM_LINES_UNKNOWN ((uint64_t)-1) 
    uint64_t num_lines;        // number of data (non-header) lines in the original txt file. Concat mode: entire file for first SectionHeaderTxtHeader, and only for that txt if not first
    uint32_t num_samples;      // VCF only: number of samples in the original VCF file
    uint32_t max_lines_per_vb; // upper bound on how many data lines a VB can have in this file
    uint8_t  compression_type; // compression type of original file, one of CompressionAlg 
    Md5Hash  md5_hash_single;  // non-0 only if this genozip file is a result of concatenatation with --md5. md5 of original single txt file.

#define TXT_FILENAME_LEN 256
    char txt_filename[TXT_FILENAME_LEN]; // filename of this single component. without path, 0-terminated. always a .vcf or .sam, even if the original was eg .vcf.gz or .bam

} SectionHeaderTxtHeader; 

// A generic VB header - used for all data types but VCF
typedef struct {
    SectionHeader h;
    uint32_t first_line;       // line (starting from 1) of this vblock in the single VCF file
                               // if this value is 0, then this is the terminating section of the file. after it is either EOF or a VCF Header section of the next concatenated file
    uint32_t num_lines;        // number of records in this vblock
    
    // features of the data
    uint32_t vb_data_size;     // size of vblock as it appears in the source file
    uint32_t z_data_bytes;     // total bytes of this vblock in the genozip file including all sections and their headers
    uint32_t longest_line_len; // length of the longest line in this vblock

    uint32_t ffu;              // for future use
} SectionHeaderVbHeader; 

// VCF VB header
typedef struct {
    SectionHeader h;
    uint32_t first_line;               // line (starting from 1) of this variant block in the single VCF file
                                       // if this value is 0, then this is the terminating section of the file. after it is either EOF or a VCF Header section of the next concatenated file
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
    uint16_t unused3;                  
} SectionHeaderVbHeaderVCF; 

typedef struct {
    SectionHeader h;           // we use h.flags to store the first 4 bits of ctx->flags, also stored in SectionHeaderCtx. This is because VCF FORMAT contexts have no SectionHeaderCtx.
    uint32_t num_snips;        // number of items in dictionary
    DictId dict_id;           
} SectionHeaderDictionary; 

typedef struct {
    SectionHeader h;
    uint8_t ltype;             // CTX_*
    uint8_t flags;             // CTX_FL_*
    uint16_t ffu;
    DictId dict_id;           
} SectionHeaderCtx;         

// the data of SEC_SECTION_LIST is an array of the following type, as is the z_file->section_list_buf
typedef struct SectionListEntry {
    uint64_t offset;           // offset of this section in the file
    DictId dict_id;            // used if this section is a DICT, LOCAL or a B250 section
    uint32_t vblock_i;
    uint8_t section_type;
    uint8_t unused[3];         
} SectionListEntry;

typedef struct {
    SectionHeader h;
    int64_t first_pos, last_pos; // first and last pos within chrom of this range         
    uint32_t chrom_word_index;   // index in context->word_list of the chrom of this reference range    
} SectionHeaderReference;

// the data of SEC_RANDOM_ACCESS is an array of the following type, as is the z_file->ra_buf and vb->ra_buf
// we maintain one RA entry per vb per every chrom in the the VB
typedef struct {
    uint32_t vblock_i;           // the vb_i in which this range appears
    uint32_t chrom_index;        // before merge: node index into chrom context mtf, after merge - word index in CHROM dictionary
    int64_t min_pos, max_pos;    // POS field value of smallest and largest POS value of this chrom in this VB (regardless of whether the VB is sorted)
} RAEntry; 

#pragma pack(pop)

// zip stuff
extern void sections_add_to_list (VBlockP vb, const SectionHeader *header);
extern void sections_list_concat (VBlockP vb, BufferP section_list_buf);

// piz stuff
extern SectionType sections_get_next_header_type (SectionListEntry **sl_ent, bool *skipped_vb, BufferP region_ra_intersection_matrix);

typedef bool (*IsSectionTypeFunc)(SectionType);
extern bool sections_get_next_section_of_type (SectionListEntry **sl_ent, uint32_t *cursor, SectionType st);

extern bool sections_has_more_components(void);
extern SectionListEntry *sections_get_offset_first_section_of_type (SectionType st);
extern SectionListEntry *sections_vb_first (uint32_t vb_i);
extern bool sections_seek_to (SectionType st);

extern void BGEN_sections_list(void);
extern const char *st_name (SectionType sec_type);
extern void sections_show_gheader (SectionHeaderGenozipHeader *header);

extern bool sections_has_reference(void);

#endif
