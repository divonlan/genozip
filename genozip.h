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
#include <math.h>
#include <time.h>
#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#include <stdbool.h>
#include <unistd.h>
#include <pthread.h>
#else
#include "compatability/visual_c_stdint.h"
#include "compatability/visual_c_stdbool.h"
#include "compatability/visual_c_unistd.h"
#include "compatability/visual_c_gettime.h"
#include "compatability/visual_c_misc_funcs.h"
#include "compatability/visual_c_pthread.h"
#endif

#if defined __APPLE__ 
#include "compatability/mac_gettime.h"
#endif
#include "version.h" // automatically incremented by the make when we create a new distribution

#define GENOZIP_EXT ".genozip"

// this was carefully picked as the optimal number based on testing with 1000-genomes chromosome 22 - 1024 samples:
// there is a tradeoff: we sort haplotypes by number of 1s for this entire variant block, and maintain and Index that is 
// #haplotypes x roundup(log2(#haplotypes)) - and we the number of indices we have in the final file depends on the number
// of variant blocks - the less the better. On the other hand - smaller variant blocks result in more linkage disequilibrium between
// variants in the block (assuming the VCF is sorted by POS), resulting in our sorting by number of 1s more likely to result
// in similar haplotypes grouped together - improving the compression.

#define VARIANTS_PER_BLOCK 4096 // Max legal value 65535. tradeoff: larger is better compression, but in some cases might be slower retrieval speed
#define SAMPLES_PER_BLOCK  1024 // tradeoff: larger is better compression, but in some cases might be slower retrieval speed
#define MAX_PLOIDY         100  // this can be any number up to 65535, it is set to 100 to avoid memory allocation
                                // explosion in case of an error in the VCF file

#define MAX_SUBFIELDS      32   // maximum number of subfield types (except for GT) that is supported in one GENOZIP file. This constant can be increased if needed.
#define SUBFIELD_ID_LEN    ((int)sizeof(uint64_t))    // VCF spec doesn't limit the ID length, we limit it to 8 chars. zero-padded.
typedef union {
    char id[SUBFIELD_ID_LEN];   // \0-padded IDs 
    uint64_t num;
} SubfieldIdType;
#define EMPTY_SUBFIELD_ID { {0,0,0,0,0,0,0,0} }

#define AES_BLOCKLEN       16   // aes encrypts blocks of 128 bits
#define AES_KEYLEN         32   // 256 bit

// Note: the algorithm will use as many cores as it can - but there's no speed penalty for a higher MAX_COMPUTE_THREADS
// that number of cores - it will not be faster, but also not slower.
// However, each thread adds memory consumption approximately linearly

// the optimal number of compute threads is determined by the ratio between CPU time and I/O time. 
// For uncompress, 3 or more threads result in similar performance for HD and SSD, with 7 seaming to be about optimal (but not a big difference than 3). 
// Adding threads doesn't help. 2 or less threads result in significantly slower execution time. 
// Memory consumption is linear with the number of threads (each allocated a VB)

#define DEFAULT_MAX_THREADS 8 // maximum compute threads created - one I/O thread, and the rest of compute threads. This can be changed with command line option -@
 
#define MAX_32BIT_WINDOWS_MEMORY (1.7*1024*1024*1024) // 1.7GB - so Windows 32bit code doesn't explode at 2GB. TO DO - make this platform specific or get ulimits

typedef enum { PHASE_UNKNOWN      = '-',
               PHASE_HAPLO        = '1',
               PHASE_PHASED       = '|',
               PHASE_NOT_PHASED   = '/',
               PHASE_MIXED_PHASED = '+'    } PhaseType;

typedef enum {
    SEC_VCF_HEADER, SEC_VARIANT_DATA, SEC_DICTIONARY, SEC_GENOTYPE_DATA, SEC_PHASE_DATA, SEC_HAPLOTYPE_DATA,
    SEC_STATS_HT_SEPERATOR // not a real section, just for stats
} SectionType;
#define NUM_SEC_TYPES (SEC_STATS_HT_SEPERATOR+1)

// Section headers - big endian

#define GENOZIP_MAGIC 0x27052012

#pragma pack(push, 1) // structures that are part of the genozip format are packed.

typedef union { 
    uint8_t  bytes[16]; 
    uint32_t words[4];
    uint64_t ulls[2];
} Md5Hash;

typedef struct {
    uint32_t magic; 
    uint32_t compressed_offset;     // number of bytes from the start of the header that is the start of compressed data (sizeof header + header encryption padding)
    uint32_t data_encrypted_len;    // = data_compressed_len + padding if encrypted, 0 if not
    uint32_t data_compressed_len;
    uint32_t data_uncompressed_len;
    uint32_t variant_block_i;       // VB with in file starting from 1 ; 0 for VCF Header
    uint16_t section_i;             // section within VB - 0 for Variant Data
    uint8_t  section_type;          
    uint8_t  unused;                // padding + for future use
} SectionHeader; 

// The VCF header section appears once in the file (or multiple times in case of concatenation), and includes the VCF file header 
typedef struct {
    SectionHeader h;
    uint8_t  genozip_version;
    uint32_t num_samples;    // number of samples in the original VCF file

    uint64_t vcf_data_size;  // number of bytes in the original VCF file. Concat mode: entire file for first SectionHeaderVCFHeader, and only for that VCF if not first
#define NUM_LINES_UNKNOWN ((uint64_t)-1) 
    uint64_t num_lines;      // number of variants (data lines) in the original VCF file. Concat mode: entire file for first SectionHeaderVCFHeader, and only for that VCF if not first
    Md5Hash md5_hash_concat; // md5 of original VCF file, or 0s if no hash was calculated. if this is a concatenation - this is the md5 of the entire concatenation.
    Md5Hash md5_hash_single; // non-0 only if this genozip file is a result of concatenatation with --md5. md5 of original single VCF file.

#define VCF_FILENAME_LEN 256
    char vcf_filename[VCF_FILENAME_LEN];    // filename of this single component. without path, 0-terminated.

#define FILE_METADATA_LEN 72
    char created[FILE_METADATA_LEN];    
} SectionHeaderVCFHeader; 

// The variant data section appears for each variant block

// Note: the reason we index the haplotypes across all samples rather than for each sample group, despite more bits in haplotype_index entry
// is that the better global sorting results in overall smaller file. 
// 
// Tested with the first 1024 variant block of chr22 of 1000-genomes (5096 haplotypes):
//   - Sorting all haplotypes together:                                        62.8KB Index 13bit x 5096 = 8.3K Total: 71.1K <--- better sort together despite bigger index
//   - Sorting each 1024 haplotypes (5 sample groups) and sorting separately - 68.1KB Index 10bit x 5096 = 6.4K Total: 74.5K
//
// Note: this doesn't affect retrieval time of a specific sample, because the sample block is just looked up via the index
 
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
    uint16_t num_subfields;            // number of subfields used in this VB. the actual fields are deduced from the FORMAT column on each line
    uint16_t num_dictionary_sections;  // we have a dictionary for each subfield, in which new words were introduced in this VB that weren't already in the dictionary. so num_dictionary_sections <= num_subfields
    uint32_t max_gt_line_len;

#define MAX_CHROM_LEN 64
    char chrom[MAX_CHROM_LEN];         // a null-terminated ID of the chromosome
    int64_t min_pos, max_pos;          // minimum and maximum POS values in this VB. -1 if unknown

    uint32_t vb_data_size;             // size of variant block as it appears in the source file
    uint32_t z_data_bytes;             // total bytes of this variant block in the genozip file including all sections and their headers
    uint16_t haplotype_index_checksum;
    uint8_t haplotype_index[];         // length is num_haplotypes. e.g. the first entry shows for the first haplotype in the original file, its index into the permuted block. # of bits per entry is roundup(log2(num_samples*ploidy)).
} SectionHeaderVariantData; 

typedef struct {
    SectionHeader h;
    uint32_t num_snips;                // number of items in dictionary
    char subfield_id[SUBFIELD_ID_LEN]; // \0-padded id
} SectionHeaderDictionary; 

#pragma pack(pop)

typedef enum {BUF_UNALLOCATED=0, BUF_REGULAR, BUF_FULL_OVERLAY, BUF_PARTIAL_OVERLAY} BufferType; // BUF_UNALLOCATED must be 0

typedef struct {
    BufferType type;
    bool overlayable; // this buffer may be fully overlaid by one or more overlay buffers
    const char *name; // name of allocator - used for memory debugging & statistics
    unsigned param;   // parameter provided by allocator - used for memory debugging & statistics
    unsigned size;    // number of bytes allocated to memory
    unsigned len;     // used by the buffer user according to its internal logic. not modified by malloc/realloc, zeroed by buf_free
    char *data;       // =memory+2*sizeof(long long) if buffer is allocated or NULL if not
    char *memory;     // memory allocated to this buffer - amount is: size + 2*sizeof(longlong) to allow for OVERFLOW and UNDERFLOW)
} Buffer;
#define EMPTY_BUFFER {BUF_UNALLOCATED,false,NULL,0,0,0,NULL,NULL}

#define NUM_VB_POOLS 2 
typedef enum { POOL_ID_UNIT_TEST=-1, POOL_ID_ZIP=0, POOL_ID_UNZIP=1 } VariantBlockPoolID;

typedef struct {
    long long wallclock, read, compute, compressor, write, zfile_read_one_vb, piz_get_variant_data_line, 
        piz_get_haplotype_data_line, piz_get_line_get_num_subfields,
        piz_get_genotype_sample_starts, piz_get_line_subfields, piz_merge_line, 
        piz_get_phase_data_line, piz_get_genotype_data_line, zfile_uncompress_section,
        piz_reconstruct_line_components, squeeze, piz_decode_pos, buf_alloc,
        seg_all_data_lines, zip_generate_haplotype_sections, sample_haplotype_data, count_alt_alleles,
        zip_generate_genotype_sections, zip_generate_phase_sections, zip_generate_variant_data_section,
        mtf_integrate_dictionary_fragment, mtf_clone_ctx, mtf_merge_in_vb_ctx, gl_optimize_dictionary,
        tmp1, tmp2, tmp3, tmp4, tmp5;
} ProfilerRec;

#pragma pack(push, 1) // packing Base250, HashEnt, MtfNode - as these might consume the majority of RAM when zipping in really big files

// values 0 to 249 are used as numerals in base-250. 
// The remaining 6 values below are control characters, and can only appear in numerals[0].
#define BASE250_EMPTY_SF   250 // subfield declared in FORMAT is empty, terminating : present
#define BASE250_MISSING_SF 251 // subfield declared in FORMAT is missing at end of cell, no :
#define BASE250_2_NUMERALS 253 // this number has 2 numerals, starting from numerals[1]
#define BASE250_3_NUMERALS 254 // this number has 3 numerals
#define BASE250_4_NUMERALS 255 // this number has 4 numerals
typedef struct {
    uint8_t numerals[5];       // number in base-250 (i.e. 0-249), digit[0] is least significant digit. If num_numerals=2,3 or 4 then numerals[1] is 250,251 or 252
    uint8_t num_numerals;      // legal values - 1,2,3,4
    uint8_t unused[2];         // padding to 64 bit
} Base250;

extern Base250 base250_encode (uint32_t n);
extern uint32_t base250_decode (const uint8_t *str);

// fake mtf index values that go into genotype_data after segregation if subfields are missing
#define SEG_MAX_INDEX  0xfffffffdUL // the number just smaller than all the special values below
#define SEG_EMPTY_SF   0xfffffffeUL // subfield is missing, terminating : present
#define SEG_MISSING_SF 0xffffffffUL // subfield is missing at end of cell, no :

#define NIL -1
typedef struct {
    uint32_t char_index;       // character index into dictionary array
    uint32_t snip_len;         // not including \t terminator present in dictionary array
    Base250 word_index;        // word index into dictionary to be written to gt section
} MtfNode;

typedef struct {        
    int32_t mtf_i;             // index into MtfContext.ol_mtf (if < ol_mtf.len) or MtfContext.mtf or NIL
    int32_t next;              // linked list - index into MtfContext.hash or NIL
} HashEnt;

// used by variant_block_i=1 to collect frequency statistics
typedef struct {
    int32_t mtf_i;             // index into MtfContext.mtf
    uint32_t count;            // number of times this snip has been encoutered so far
} SorterEnt;

#pragma pack(pop)

typedef struct {
    SubfieldIdType subfield;   // which subfield is this MTF dealing with

    Buffer ol_dict;            // tab-delimited list of all unique snips - overlayed all previous VB dictionaries
    Buffer ol_mtf;             // MTF nodes - overlayed all previous VB dictionaries. char/word indeces are into ol_dict.
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    Buffer mtf;                // array of MtfNode - in this VB that don't exist in ol_mtf. char/word indeces are into dict.
    Buffer hash;               // hash table of entries HashEnt - initialized as a copy of all previous VBs. For entries are
                               // obtained by hash function hash(snip) and the rest of linked to them by linked list
    Buffer sorter;             // used by the first VB only of ZIP to sort the dictionary - entries of SorterEnt
    Buffer word_list;          // word list - used for PIZ only
} MtfContext;

typedef struct {
    uint32_t     lo, hi;
    uint32_t     a, b,c, d;
    union {
        uint8_t  bytes[64];
        uint32_t words[16];
    } buffer;
    uint32_t     block[16];
} Md5Context;

typedef enum {UNKNOWN, VCF, VCF_GZ, VCF_BZ2, GENOZIP, GENOZIP_TEST, PIPE, STDIN, STDOUT} FileType;

typedef struct {
    void *file;
    char *name;                        // allocated by file_open(), freed by file_close()
    FileType type;
    // these relate to actual bytes on the disk
    uint64_t disk_size;                // 0 if not known (eg stdin)
    uint64_t disk_so_far;              // data read/write to/from "disk" (using fread/fwrite)

    // this relate to the VCF data represented. In case of READ - only data that was picked up from the read buffer.
    uint64_t vcf_data_size_single;     // VCF: size of the VCF data (if known)
                                       // GENOZIP: GENOZIP: size of original VCF data in the VCF file currently being processed
    uint64_t vcf_data_size_concat;     // concatenated vcf_data_size of all files compressed
    uint64_t vcf_data_so_far;          // VCF: data sent to/from the caller (after coming off the read buffer and/or decompression)
                                       // GENOZIP: VCF data so far of original VCF file currently being processed

    // Used for READING VCF/VCF_GZ/VCF_BZ2 files: stats used to optimize memory allocation
    double avg_header_line_len, avg_data_line_len;   // average length of data line so far. 
    uint32_t header_lines_so_far, data_lines_so_far; // number of lines read so far

    // Used for READING & WRITING VCF files - but stored in the z_file structure for zip to support concatenation (and in the vcf_file structure for piz)
    bool has_md5;
    Md5Context md5_ctx_concat;         // md5 context of vcf file. in concat mode - of the resulting concatenated vcf file
    Md5Context md5_ctx_single;         // used only in concat mode - md5 of the single vcf component

    // Used for READING GENOZIP files
    Buffer next_vcf_header;            // next VCF header - used when reading in --split mode
    
    // Used for WRITING GENOZIP files
    uint64_t disk_at_beginning_of_this_vcf_file;     // the value of disk_size when starting to read this vcf file
    uint64_t num_lines_concat;                // number of lines in concatenated vcf file
    uint64_t num_lines_single;                // number of lines in single vcf file
    
    SectionHeaderVCFHeader vcf_header_first;  // store the first VCF header - we might need to update it at the very end;
    uint8_t vcf_header_enc_padding[AES_BLOCKLEN-1]; // just so we can overwrite vcf_header with encryption padding

    SectionHeaderVCFHeader vcf_header_single; // store the VCF header of single component in concat mode
    uint8_t vcf_header_enc_padding2[AES_BLOCKLEN-1]; // same

    // dictionary information used for writing GENOZIP files - can be accessed only when holding mutex
    pthread_mutex_t mutex;
    bool mutex_initialized;
    unsigned next_variant_i_to_merge;  // merging vb's dictionaries into mtf_ctx needs to be in the variant_block_i order
    unsigned num_subfields;            // length of populated subfield_ids and mtx_ctx;
    MtfContext mtf_ctx[MAX_SUBFIELDS]; // a merge of dictionaries of all VBs

    // Information content stats - how many bytes does this file have in each section type
    uint64_t section_bytes[NUM_SEC_TYPES];   

    // USED FOR READING ALL FILES
#   define READ_BUFFER_SIZE (1<<19)    // 512KB
    uint32_t next_read, last_read;     // indices into read_buffer
    bool eof;                          // we reached EOF
    char read_buffer[];                // only allocated for mode=READ files   
} File;

// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {

    uint32_t line_i;

    // initially, data from vcf file line, later segregated to components "stored" in the overlay buffers below
    Buffer line;             

    // the following 4 buffers are overlay buffers onto line, so that they don't consume memory
    Buffer variant_data;     // string terminated by a newline. len includes the newline.
    Buffer genotype_data;    // \t separated genotype data for the line. last one has \t too. no \0. exists if any variant in the variant blck has FORMAT other that "GT"
    Buffer haplotype_data;   // length=ploidy*num_samples. exists if the GT subfield exists in any variant in the variant block

    Buffer phase_data;       // used only if phase is mixed. length=num_samples. exists if haplotype data exists and ploidy>=2
    PhaseType phase_type;    // phase type of this line
    
    bool has_haplotype_data; // FORMAT field contains GT
    bool has_genotype_data;  // FORMAT field contains subfields other than GT

    unsigned num_subfields;
    unsigned sf_i[MAX_SUBFIELDS]; // array in the order it appears in FORMAT - each an index into vb->mtf_ctx[]
} DataLine;

// IMPORTANT: if changing fields in VariantBlock, also update vb_release_vb
typedef struct variant_block_ {

    unsigned id;               // id of vb within the vb pool
    VariantBlockPoolID pool_id; 

    File *vcf_file, *z_file;  // pointers to objects that span multiple VBs

    // memory management
    Buffer buffer_list;        // a buffer containing an array of pointers to all buffers allocated for this VB (either by the I/O thread or its compute thread)

    bool ready_to_dispatch;    // line data is read, and dispatcher can dispatch this VB to a compute thread
    bool is_processed;         // thread completed processing this VB - it is ready for outputting
    bool in_use;               // this vb is in use
        
    DataLine data_lines[VARIANTS_PER_BLOCK];
    uint32_t num_lines;        // number of lines in this variant block
    uint32_t first_line;       // line number in VCF file (counting from 1), of this variant block
    uint32_t variant_block_i;  // number of variant block within VCF file

    // tracking execution
    uint32_t vb_data_size;     // size of variant block as it appears in the source file
    uint32_t max_gt_line_len;  // length of longest gt line in this vb after segregation

    ProfilerRec profile;

    // charactaristics of the data
    uint16_t ploidy;
    uint32_t num_sample_blocks;
    uint32_t num_samples_per_block; // except last sample block that may have less
    uint32_t num_haplotypes_per_line;
    bool has_genotype_data;    // if any variant has genotype data, then the block is considered to have it
    bool has_haplotype_data;   // ditto for haplotype data
    PhaseType phase_type;      // phase type of this variant block

    // chrom and pos
    char chrom[MAX_CHROM_LEN]; // a null-terminated ID of the chromosome
    int64_t min_pos, max_pos;  // minimum and maximum POS values in this VB. -1 if unknown
    uint64_t last_pos;         // value of POS field of the previous line, to do delta encoding
    bool is_sorted_by_pos;     // true if it this variant block is sorted by POS    

    // working memory for segregate - we segregate a line components into these buffers, and when done
    // we copy it back to DataLine - the buffers overlaying the line field
    Buffer line_variant_data;  // string terminated by a newline. len includes the newline.
    Buffer line_gt_data;       // \t separated genotype data for the line. last one has \t too. no \0. exists if any variant in the variant blck has FORMAT other that "GT"
    Buffer line_ht_data;       // length=ploidy*num_samples. exists if the GT subfield exists in any variant in the variant block
    Buffer line_phase_data;    // used only if phase is mixed. length=num_samples. exists if haplotype data exists and ploidy>=2

    // crypto stuff
    Buffer spiced_pw;  // used by crypt_generate_aes_key()
    uint8_t aes_round_key[240];// for 256 bit aes
    uint8_t aes_iv[AES_BLOCKLEN];
    int bi;

    // section data - ready to compress
    Buffer variant_data_section_data;    // all fields until FORMAT, newline-separated, \0-termianted. .len includes the terminating \0
    Buffer haplotype_permutation_index;
    Buffer optimized_gl_dict;  // GL dictionary data after extra optimization

    // these are Buffer arrays of size vb->num_sample_blocks allocated once when used for the first time and never freed. 
    // Subsequent variant blocks that re-use the memory have the same number of samples by VCF spec. Each buffer is a string of the data as written to the GENOZIP file

    // Note: The sample blocks for phase and genotype data are subsequent blocks of num_samples_per_block samples.
    // in contrast the sample blocks of haplotypes are num_samples_per_block*ploidy haplotypes as permuted. 
    // e.g. genotypes and phase in sample_block_i=1 don't necessary correspond to haplotypes in sample_block_i=1.

    Buffer *haplotype_sections_data;  // this is the haplotype character for each haplotype in the transposed sample group
    Buffer *phase_sections_data;      // this is the phase character for each genotype in the sample group
    Buffer *genotype_sections_data;   // this is for piz - each entry is a sample block, scanned columns first, each cell containing num_subfields indices (in base250 - 1 to 5 bytes each) into the subfield dictionaries
    Buffer genotype_one_section_data; // for zip we need only one section data

    // compresssed file data 
    Buffer z_data;                    // all headers and section data as read from disk

    int16_t z_next_header_i;          // next header of this VB to be encrypted or decrypted

    Buffer z_section_headers;         // (used by piz) an array of unsigned offsets of section headers within z_data

    Buffer gt_sb_line_starts_buf,     // used by zip_get_genotype_vb_start_len 
           gt_sb_line_lengths_buf,
           genotype_section_lens_buf; 

    Buffer helper_index_buf;          // used by zip_do_haplotypes

    Buffer vardata_header_buf;        // used by zfile_compress_variant_data

    Buffer compressed;                // used by various zfile functions
 
    Buffer ht_columns_data;           // used by piz_get_ht_permutation_lookups

    Buffer next_gt_in_sample;         // used for reconstructing genotype data by piz

    // subfields stuff 
    unsigned num_subfields;
    MtfContext mtf_ctx[MAX_SUBFIELDS];

    // Information content stats - how many bytes does this section have more than the corresponding part of the vcf file    
    int add_bytes[NUM_SEC_TYPES];                
    uint32_t vcf_section_bytes[NUM_SEC_TYPES];  // how many bytes did each section have in the original vcf file - should add up to the file size
    uint32_t z_section_bytes[NUM_SEC_TYPES];    // how many bytes does each section type have (including headers) in the genozip file - should add up to the file size

#   define NUM_COMPRESS_BUFS 4                  // bzlib2 compress requires 4 and decompress requires 2
    Buffer compress_bufs[NUM_COMPRESS_BUFS];    // memory allocation for compressor so it doesn't do its own malloc/free

} VariantBlock;

typedef enum {READ, WRITE} FileMode;
extern File *file_open (const char *filename, FileMode mode, FileType expected_type);
extern File *file_fdopen (int fd, FileMode mode, FileType type, bool initialize_mutex);
extern void file_close (File **vcf_file_p);
extern size_t file_write (File *file, const void *data, unsigned len);
extern void file_remove (const char *filename);
extern bool file_has_ext (const char *filename, const char *extension);
extern const char *file_basename (const char *filename, bool remove_exe, const char *default_basename,
                                  char *basename /* optional pre-allocated memory */, unsigned basename_size /* basename bytes */);

extern bool vcffile_get_line(VariantBlock *vb, unsigned line_i_in_file, bool skip_md5_vcf_header, Buffer *line, const char *buf_name);
extern void vcffile_write_one_variant_block (File *vcf_file, VariantBlock *vb);
extern unsigned vcffile_write_to_disk(File *vcf_file, const Buffer *buf);
extern void vcffile_compare_pipe_to_file (FILE *from_pipe, File *vcf_file);

// reads VCF header and writes its compressed form to the GENOZIP file. returns num_samples.
extern bool vcf_header_vcf_to_genozip (VariantBlock *vb, unsigned *line_i, Buffer **first_data_line);
extern bool vcf_header_genozip_to_vcf (VariantBlock *vb, Md5Hash *digest /* out */);
extern bool vcf_header_get_vcf_header (File *z_file, SectionHeaderVCFHeader *vcf_header_header, bool *encrypted);

extern void vb_cleanup_memory (VariantBlockPoolID pool_id);
extern VariantBlock *vb_get_vb (VariantBlockPoolID pool_id, File *vcf_file, File *z_file, unsigned variant_block_i);
extern unsigned vb_num_samples_in_sb (const VariantBlock *vb, unsigned sb_i);
extern unsigned vb_num_sections(VariantBlock *vb);
extern void vb_release_vb (VariantBlock **vb_p);

typedef struct {
    unsigned num_vbs;
    VariantBlock vb[]; // variable length
} VariantBlockPool;
extern void vb_create_pool (VariantBlockPoolID pool_id, unsigned num_vbs);
extern VariantBlockPool *vb_get_pool (VariantBlockPoolID pool_id);


extern void seg_all_data_lines (VariantBlock *vb, Buffer *lines_orig /* for testing */);
extern SubfieldIdType seg_get_subfield (const char **data, unsigned len, unsigned line_i);

extern int32_t mtf_evaluate_snip (VariantBlock *vb, MtfContext *ctx, const char *snip, uint32_t snip_len, bool overlayable, MtfNode **node /* out */);
extern void mtf_get_snip_by_word_index (VariantBlock *vb, MtfContext *ctx, const uint8_t *word_index_base250, char **snip, uint32_t *snip_len);
extern void mtf_clone_ctx (VariantBlock *vb);
extern MtfNode *mtf_node (const MtfContext *ctx, uint32_t mtf_i, const char **snip_in_dict /* optional out */);
extern unsigned mtf_merge_in_vb_ctx (VariantBlock *vb);
extern unsigned mtf_get_sf_i_by_subfield (MtfContext *mtf_ctx, unsigned *num_subfields, SubfieldIdType subfield);
extern void mtf_integrate_dictionary_fragment (VariantBlock *vb, char *data);
extern void mtf_overlay_dictionaries_to_vb (VariantBlock *vb);
extern void mtf_sort_dictionaries_vb_1(VariantBlock *vb);
extern void mtf_initialize_mutex (File *z_file, unsigned next_variant_i_to_merge);

extern void mtf_free_context (MtfContext *ctx);
#ifdef DEBUG
extern void mtf_tree_test (const MtfContext *ctx);
#endif

extern const char *gl_optimize_dictionary (VariantBlock *vb, Buffer *dict, MtfNode *nodes, unsigned dict_start_char, unsigned num_words);
extern void gl_deoptimize_dictionary (char *data, int len);

extern void zip_dispatcher (const char *vcf_basename, File *vcf_file, 
                            File *z_file, bool test_mode, unsigned max_threads, bool is_last_file);

// returns true is successfully outputted a vcf file
extern bool piz_dispatcher (const char *z_basename, File *z_file, File *vcf_file, 
                            bool test_mode, unsigned max_threads);

extern void piz_reconstruct_line_components (VariantBlock *vb);
extern void piz_merge_all_lines (VariantBlock *vb);

typedef void *Dispatcher;
extern Dispatcher dispatcher_init (unsigned max_threads, VariantBlockPoolID pool_id, unsigned previous_vb_i, File *vcf_file, File *z_file,
                                   bool test_mode, bool show_progress, const char *filename);
extern void dispatcher_pause (Dispatcher dispatcher);
extern void dispatcher_resume (Dispatcher dispatcher, File *vcf_file);
extern void dispatcher_finish (Dispatcher *dispatcher, unsigned *last_vb_i);

typedef void (*DispatcherFuncType)(VariantBlock *);
extern void dispatcher_compute (Dispatcher dispatcher, DispatcherFuncType func);
extern VariantBlock *dispatcher_generate_next_vb (Dispatcher dispatcher);       
extern bool dispatcher_has_processed_vb (Dispatcher dispatcher, bool *is_final);                                  
extern VariantBlock *dispatcher_get_processed_vb (Dispatcher dispatcher, bool *is_final);
extern bool dispatcher_has_free_thread (Dispatcher dispatcher);
extern VariantBlock *dispatcher_get_pseudo_vb (Dispatcher dispatcher);
extern VariantBlock *dispatcher_get_next_vb (Dispatcher dispatcher);
extern void dispatcher_finalize_one_vb (Dispatcher dispatcher, const File *file, long long vcf_data_written_so_far, uint64_t bytes_compressed);
extern void dispatcher_input_exhausted (Dispatcher dispatcher);
extern bool dispatcher_is_done (Dispatcher dispatcher);
extern bool dispatcher_is_input_exhausted (Dispatcher dispatcher);

extern void zfile_write_vcf_header (VariantBlock *vb, Buffer *vcf_header_text, bool is_first_vcf);
extern void zfile_compress_variant_data (VariantBlock *vb);
extern void zfile_update_compressed_variant_data_header (VariantBlock *vb, unsigned pos, unsigned num_dictionary_sections);
extern void zfile_compress_section_data (VariantBlock *vb, SectionType section_type, Buffer *section_data);
extern void zfile_compress_dictionary_data (VariantBlock *vb, SubfieldIdType subfield, 
                                            uint32_t num_words, const char *data, uint32_t num_chars);

extern bool zfile_read_one_vb (VariantBlock *vb);

// returns offset of header within data, EOF if end of file (or end of VCF component in the case of flag_split)
extern int zfile_read_one_section (VariantBlock *vb, 
                                   Buffer *data /* buffer to append */, const char *buf_name,
                                   unsigned header_size, SectionType expected_sec_type,
                                   bool allow_eof);

extern void zfile_uncompress_section(VariantBlock *vb, void *section_header, Buffer *uncompressed_data, SectionType expected_section_type);

#ifdef __APPLE__
#define off64_t __int64_t // needed for conda mac - otherwise zlib.h throws compilation errors
#endif
extern void zfile_update_vcf_header_section_header (VariantBlock *vb, off64_t vcf_header_header_pos_single, bool final_for_concat);

extern void crypt_set_password (char *new_password);
extern bool crypt_have_password ();
extern bool crypt_prompt_for_password();
extern unsigned crypt_padded_len (unsigned len);
extern bool crypt_get_encrypted_len (unsigned *data_encrypted_len /* in/out */, unsigned *padding_len /* out */);
extern void crypt_do (VariantBlock *vb, uint8_t *data, unsigned data_len, uint32_t vb_i, int16_t sec_i);
extern void crypt_continue (VariantBlock *vb, uint8_t *data, unsigned data_len);
extern void crypt_pad (uint8_t *data, unsigned data_len, unsigned padding_len);
extern unsigned crypt_max_padding_len();

extern void aes_initialize (VariantBlock *vb, const uint8_t *key);
extern void aes_xcrypt_buffer (VariantBlock *vb, uint8_t *data, uint32_t length);
extern char *aes_display_key (const uint8_t* key);
extern char *aes_display_data (const uint8_t *data, unsigned data_len);

extern void md5_do (const void *data, unsigned len, Md5Hash *digest);
extern void md5_update (Md5Context *ctx, const void *data, unsigned len, bool initialize);
extern void md5_finalize (Md5Context *ctx, Md5Hash *digest);
const char *md5_display (const Md5Hash *digest, bool prefix_space);
extern void md5_display_ctx (const Md5Context *x);

extern void squeeze (VariantBlock *vb,
                     uint8_t *dst, // memory should be pre-allocated by caller
                     uint16_t *squeezed_checksum,
                     const unsigned *src, 
                     unsigned src_len);

extern void unsqueeze (VariantBlock *vb,
                       unsigned *normal, // memory should be pre-allocated by caller
                       const uint8_t *squeezed, 
                       uint16_t squeezed_checksum,
                       unsigned normal_len);

extern unsigned squeeze_len(unsigned int len);

extern void buf_initialize();
extern unsigned buf_alloc (VariantBlock *vb,
                           Buffer *buf, 
                           unsigned requested_size, // whether contents of memory should be zeroed
                           float grow_at_least_factor, // grow more than new_size    
                           const char *name, unsigned param); // for debugging
static inline void buf_set_overlayable (Buffer *buf) { buf->overlayable = true;}
extern void buf_overlay (Buffer *overlaid_buf, Buffer *regular_buf, const Buffer *copy_from, unsigned *regular_buf_offset, const char *name, unsigned param);
extern void buf_free (Buffer *buf); // free buffer - without freeing memory. A future buf_alloc of this buffer will reuse the memory if possible.
extern void buf_destroy (VariantBlock *vb, Buffer *buf);

static inline bool buf_is_allocated (const Buffer *buf) {return buf->data != NULL && (buf->type == BUF_REGULAR || buf->type == BUF_FULL_OVERLAY || buf->type == BUF_PARTIAL_OVERLAY);}

extern void buf_copy (VariantBlock *vb, Buffer *dst, const Buffer *src, unsigned bytes_per_entry,
                      unsigned start_entry, unsigned max_entries, // if 0 copies the entire buffer
                      const char *name, unsigned param);

extern void buf_move (VariantBlock *vb, Buffer *dst, Buffer *src);

static inline void buf_add (Buffer *buf, const void *data, unsigned len) { memcpy (&buf->data[buf->len], data, len);  buf->len += len; }
#define buf_add_string(buf,str) buf_add (buf, str, strlen (str));

extern void buf_test_overflows(const VariantBlock *vb);

extern long long buf_vb_memory_consumption (const VariantBlock *vb);
extern void buf_display_memory_usage (VariantBlockPoolID pool_id, bool memory_full);

extern char *buf_human_readable_size (uint64_t size, char *str /* out */);

// global parameters - set before any thread is created, and never change
extern unsigned    global_num_samples;
extern const char *global_cmd;            // set once in main()
extern bool        global_little_endian;  // set in main()

// flags set by user's command line options
extern int flag_force, flag_quiet, flag_concat_mode, flag_md5, flag_split, flag_show_alleles, flag_show_time, flag_show_memory;

// macros
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b) )
#define MAX(a, b) (((a) > (b)) ? (a) : (b) )
#endif

// encode section headers in Big Endian (see https://en.wikipedia.org/wiki/Endianness)
// the reason for selecting big endian is that I am developing on little endian CPU (Intel) so
// endianity bugs will be discovered more readily this way
#ifdef _MSC_VER
#define __builtin_bswap16 _byteswap_ushort
#define __builtin_bswap32 _byteswap_ulong
#define __builtin_bswap64 _byteswap_uint64
#endif

#define BGEN16(x) (global_little_endian ? __builtin_bswap16(x) : (x))
#define BGEN32(x) (global_little_endian ? __builtin_bswap32(x) : (x))
#define BGEN64(x) (global_little_endian ? __builtin_bswap64(x) : (x))

#define LTEN16(x) (global_little_endian ? (x) : __builtin_bswap16(x))
#define LTEN32(x) (global_little_endian ? (x) : __builtin_bswap32(x))
#define LTEN64(x) (global_little_endian ? (x) : __builtin_bswap64(x))

// sanity checks
static inline void my_exit() { exit(1); }// an exit function so we can put a debugging break point when ASSERT exits
#define ASSERT(condition, format, ...)  { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); my_exit(); }}
#define ASSERT0(condition, string)      { if (!(condition)) { fprintf (stderr, "\n%s\n", string); my_exit(); }}
#define ASSERTW(condition, format, ...) { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); }}
#define ASSERTW0(condition, string)     { if (!(condition)) { fprintf (stderr, "\n%s\n", string); } }
#define ABORT(format, ...)              { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); my_exit();}
#define ABORT0(string)                  { fprintf (stderr, "\n%s\n", string); my_exit();}

#ifdef _MSC_VER
typedef struct my_timespec TimeSpecType;
#else
typedef struct timespec TimeSpecType;
#endif

#define START_TIMER     TimeSpecType profiler_timer; \
                        if (flag_show_time) clock_gettime(CLOCK_REALTIME, &profiler_timer); 

#define COPY_TIMER(res) if (flag_show_time) { \
                            TimeSpecType tb; \
                            clock_gettime(CLOCK_REALTIME, &tb); \
                            res += (tb.tv_sec-profiler_timer.tv_sec)*1000000000ULL + (tb.tv_nsec-profiler_timer.tv_nsec); \
                        }

extern void profiler_add (ProfilerRec *dst, const ProfilerRec *src);
extern const char *profiler_print_short (const ProfilerRec *p);
extern void profiler_print_report (const ProfilerRec *p, unsigned max_threads, unsigned used_threads, const char *filename, unsigned num_vbs);

// a hacky addition to bzip2
extern unsigned long long BZ2_bzoffset (void* b);

// Windows compatibility stuff
#ifdef _WIN32
#define stat64  _stat64
#define fstat64 _fstat64
#else // this needs more work - there are more cases, depending if gcc is 32 or 64
#define stat64  stat
#define fstat64 fstat
#endif

#endif


