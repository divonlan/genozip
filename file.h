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
#include "compatibility/visual_c_stdint.h"
#include "compatibility/visual_c_stdbool.h"
#include "compatibility/visual_c_pthread.h"
#endif

#include "buffer.h"
#include "sections.h"
#include "move_to_front.h"
#include "aes.h"

#define VCF_         ".vcf"
#define VCF_GZ_      ".vcf.gz"
#define VCF_BGZ_     ".vcf.bgz"
#define VCF_BZ2_     ".vcf.bz2"
#define VCF_XZ_      ".vcf.xz"
#define BCF_         ".bcf"
#define BCF_GZ_      ".bcf.gz"
#define BCF_BGZ_     ".bcf.bgz"
#define VCF_GENOZIP_ ".vcf" GENOZIP_EXT

typedef enum      {UNKNOWN_FILE_TYPE, VCF,  VCF_GZ,  VCF_BGZ,  VCF_BZ2,  VCF_XZ,  BCF,  BCF_GZ,  BCF_BGZ,  VCF_GENOZIP,  STDIN,   STDOUT} FileType;
#define FILE_EXTS {"Unknown",   VCF_, VCF_GZ_, VCF_BGZ_, VCF_BZ2_, VCF_XZ_, BCF_, BCF_GZ_, BCF_BGZ_, VCF_GENOZIP_, "stdin", "stdout" }
extern char *file_exts[];

#define VCF_EXTENSIONS VCF_ " " VCF_GZ_ " " VCF_BGZ_ " " VCF_BZ2_ " " VCF_XZ_ " " BCF_ " " BCF_GZ_ " " BCF_BGZ_

typedef const char *FileMode;
extern FileMode READ, WRITE; // this are pointers to static strings - so they can be compared eg "if (mode==READ)"

#define file_is_zip_read(file) ((file)->mode == READ && command == ZIP)

#define file_is_via_ext_decompressor(file) ((file)->type == VCF_XZ || \
                                            (file)->type == BCF || (file)->type == BCF_GZ  || (file)->type == BCF_BGZ)

// files that are read by ZIP as plain VCF - they might have been decompressed by an external decompressor
#define file_is_plain_vcf(file) (((file)->type == VCF || (file)->type == STDIN) || file_is_via_ext_decompressor(file))

#define file_is_vcf(file) (file_is_plain_vcf(file) || (file)->type == VCF_GZ || (file)->type == VCF_BGZ || (file)->type == VCF_BZ2)

typedef struct file_ {
    void *file;
    char *name;                        // allocated by file_open(), freed by file_close()
    FileMode mode;
    FileType type;
    
    // these relate to actual bytes on the disk
    int64_t disk_size;                 // 0 if not known (eg stdin or http stream). 
                                       // note: this is different from vcf_data_size_single as disk_size might be compressed (gz, bz2 etc)
    int64_t disk_so_far;               // data read/write to/from "disk" (using fread/fwrite)

    // this relate to the VCF data represented. In case of READ - only data that was picked up from the read buffer.
    int64_t vcf_data_size_single;      // vcf_file: size of the VCF data. ZIP: if its a plain VCF file, then its the disk_size. If not, we initially do our best to estimate the size, and update it when it becomes known.
    int64_t vcf_data_so_far_single;    // vcf_file: data read (ZIP) or written (PIZ) to/from vcf file so far
                                       // z_file: VCF data represented in the GENOZIP data written (ZIP) or read (PIZ) to/from the genozip file so far for the current VCF
    int64_t vcf_data_so_far_concat;    // z_file & ZIP only: VCF data represented in the GENOZIP data written so far for all VCFs
    int64_t num_lines;                 // z_file: number of lines in all vcf files concatenated into this z_file
                                       // vcf_file: number of lines in single vcf file

    // Used for READING & WRITING VCF files - but stored in the z_file structure for zip to support concatenation (and in the vcf_file structure for piz)
    Md5Context md5_ctx_concat;         // md5 context of vcf file. in concat mode - of the resulting concatenated vcf file
    Md5Context md5_ctx_single;         // used only in concat mode - md5 of the single vcf component
    uint32_t max_lines_per_vb;         // ZIP & PIZ - in ZIP, discovered while segmenting, in PIZ - given by SectionHeaderVCFHeader

    // Used for READING GENOZIP files
    Buffer v1_next_vcf_header;         // genozip v1 only: next VCF header - used when reading in --split mode
    uint8_t genozip_version;           // GENOZIP_FILE_FORMAT_VERSION of the genozip file being read
    uint32_t num_vcf_components;       // set from genozip header

    // Used for WRITING GENOZIP files
    uint64_t disk_at_beginning_of_this_vcf_file;     // z_file: the value of disk_size when starting to read this vcf file
    
    SectionHeaderVCFHeader vcf_header_first;         // store the first VCF header - we might need to update it at the very end;
    uint8_t vcf_header_enc_padding[AES_BLOCKLEN-1];  // just so we can overwrite vcf_header with encryption padding

    SectionHeaderVCFHeader vcf_header_single;        // store the VCF header of single component in concat mode
    uint8_t vcf_header_enc_padding2[AES_BLOCKLEN-1]; // same

    // dictionary information used for writing GENOZIP files - can be accessed only when holding mutex
    pthread_mutex_t dicts_mutex;
    bool dicts_mutex_initialized;      // this mutex protects mtf_ctx and num_dict_ids from concurrent adding of a new dictionary
    unsigned num_dict_ids;             // length of populated subfield_ids and mtx_ctx;
    
    MtfContext mtf_ctx[MAX_DICTS];     // a merge of dictionaries of all VBs
    Buffer ra_buf;                     // RAEntry records - in a format ready to write to disk (Big Endian etc)
    Buffer dict_data;                  // Dictionary data accumulated from all VBs and written near the end of the file

    // section list - used for READING and WRITING genozip files
    Buffer section_list_buf;           // section list to be written as the payload of the genotype header section
    Buffer section_list_dict_buf;      // ZIP: a subset of section_list_buf - dictionaries are added here by VBs as they are being constructed
    uint32_t sl_cursor, sl_dir_cursor; // PIZ: next index into section_list for searching for sections
    uint32_t num_vcf_components_so_far;

    // Information content stats - how many bytes and how many sections does this file have in each section type
    int64_t section_bytes[NUM_SEC_TYPES];   
    int64_t section_entries[NUM_SEC_TYPES]; // used only for Z files - number of entries of this type (dictionary entries or base250 entries)
    uint32_t num_sections[NUM_SEC_TYPES];    // used only for Z files - number of sections of this type

#   define READ_BUFFER_SIZE (1<<19)    // 512KB
    // Used for reading vcf files
    Buffer vcf_unconsumed_data;         // excess data read from the vcf file - moved to the next VB
    
    // Used for reading genozip files
    uint32_t next_read, last_read;     // indices into read_buffer
    char read_buffer[];                // only allocated for mode=READ files   
} File;

// globals
extern File *z_file, *vcf_file; 

// methods
extern File *file_open (const char *filename, FileMode mode, FileType expected_type);
extern File *file_fdopen (int fd, FileMode mode, FileType type, bool initialize_mutex);
extern void file_close (FileP *file_p, bool cleanup_memory /* optional */);
extern size_t file_write (FileP file, const void *data, unsigned len);
extern bool file_seek (File *file, int64_t offset, int whence, bool soft_fail); // SEEK_SET, SEEK_CUR or SEEK_END
extern uint64_t file_tell (File *file);
extern uint64_t file_get_size (const char *filename);
extern bool file_is_dir (const char *filename);
extern void file_remove (const char *filename, bool fail_quietly);
extern bool file_has_ext (const char *filename, const char *extension);
extern const char *file_basename (const char *filename, bool remove_exe, const char *default_basename,
                                  char *basename /* optional pre-allocated memory */, unsigned basename_size /* basename bytes */);
extern void file_get_file (VariantBlockP vb, const char *filename, Buffer *buf, const char *buf_name, unsigned buf_param, bool add_string_terminator);
extern void file_assert_ext_decompressor (void);
extern void file_kill_external_compressors (void);

#define file_printname(file) ((file)->name ? (file)->name : ((file)->mode==READ ? "(stdin)" : "(stdout)"))

#define CLOSE(fd,name)  { ASSERTW (!close (fd),  "Warning in %s:%u: Failed to close %s: %s",  __FUNCTION__, __LINE__, (name), strerror(errno));}
#define FCLOSE(fp,name) { if (fp) { ASSERTW (!fclose (fp), "Warning in %s:%u: Failed to fclose %s: %s", __FUNCTION__, __LINE__, (name), strerror(errno)); fp = NULL; } }
 
// a hacky addition to bzip2
extern unsigned long long BZ2_consumed (void* b);
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
