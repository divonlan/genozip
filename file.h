// ------------------------------------------------------------------
//   file.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FILE_INCLUDED
#define FILE_INCLUDED

#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#include <stdbool.h>
#include <pthread.h>
#else
#include "compatability/visual_c_stdint.h"
#include "compatability/visual_c_stdbool.h"
#include "compatability/visual_c_pthread.h"
#endif

#include "buffer.h"
#include "sections.h"
#include "move_to_front.h"
#include "aes.h"

typedef enum {UNKNOWN, VCF, VCF_GZ, VCF_BZ2, GENOZIP, GENOZIP_TEST, PIPE, STDIN, STDOUT} FileType;

typedef struct file_ {
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

typedef enum {READ, WRITE} FileMode;
extern File *file_open (const char *filename, FileMode mode, FileType expected_type);
extern File *file_fdopen (int fd, FileMode mode, FileType type, bool initialize_mutex);
extern void file_close (FileP *vcf_file_p);
extern size_t file_write (FileP file, const void *data, unsigned len);
extern void file_remove (const char *filename);
extern bool file_has_ext (const char *filename, const char *extension);
extern const char *file_basename (const char *filename, bool remove_exe, const char *default_basename,
                                  char *basename /* optional pre-allocated memory */, unsigned basename_size /* basename bytes */);

// a hacky addition to bzip2
extern unsigned long long BZ2_bzoffset (void* b);
extern const char *BZ2_errstr (int err);

// Windows compatibility stuff
#ifdef _WIN32
#define stat64  _stat64
#define fstat64 _fstat64
#else // this needs more work - there are more cases, depending if gcc is 32 or 64
#define stat64  stat
#define fstat64 fstat
#endif

#endif
