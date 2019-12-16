// ------------------------------------------------------------------
//   vczip.h
//   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VCZIP_INCLUDED
#define VCZIP_INCLUDED

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#ifdef __APPLE__
#include "mac/mach_gettime.h"
#endif

#define VCZIP_VERSION 1 // legal value 0-255. this needs to be incremented when the dv file format changes

// this was carefully picked as the optimal number based on testing with 1000-genomes chromosome 22 - 1024 samples:
// there is a tradeoff: we sort haplotypes by number of 1s for this entire variant block, and maintain and Index that is 
// #haplotypes x roundup(log2(#haplotypes)) - and we the number of indices we have in the final file depends on the number
// of variant blocks - the less the better. On the other hand - smaller variant blocks result in more linkage disequilibrium between
// variants in the block (assuming the VCF is sorted by POS), resulting in our sorting by number of 1s more likely to result
// in similar haplotypes grouped together - improving the compression.

#define VARIANTS_PER_BLOCK 4096 // Max legal value 65535. tradeoff: larger is better compression, but in some cases might be slower retrieval speed
#define SAMPLES_PER_BLOCK  1024 // tradeoff: larger is better compression, but in some cases might be slower retrieval speed

// Note: the algorithm will use as many cores as it can - but there's no speed penalty for a higher MAX_COMPUTE_THREADS
// that number of cores - it will not be faster, but also not slower.
// However, each thread adds memory consumption approximately linearly

// the optimal number of compute threads is determined by the ratio between CPU time and I/O time. 
// For uncompress, 3 or more threads result in similar performance for HD and SSD, with 7 seaming to be about optimal (but not a big difference than 3). 
// Adding threads doesn't help. 2 or less threads result in significantly slower execution time. 
// Memory consumption is linear with the number of threads (each allocated a VB)

#define DEFAULT_MAX_THREADS 8 // maximum compute threads created - one I/O thread, and the rest of compute threads. This can be changed with command line option -@
 
#define MAX_MEMORY (1.7*1024*1024*1024) // 1.7GB - so Windows 32bit code doesn't explode at 2GB. TO DO - make this platform specific or get ulimits

//#define PROFILER 
typedef enum { PHASE_UNKNOWN      = '-',
               PHASE_HAPLO        = '1',
               PHASE_PHASED       = '|',
               PHASE_NOT_PHASED   = '/',
               PHASE_MIXED_PHASED = '+'    } PhaseType;

typedef enum {
    SEC_VCF_HEADER, SEC_VARIANT_DATA, SEC_GENOTYPE_DATA, SEC_PHASE_DATA, SEC_HAPLOTYPE_DATA,
    SEC_STATS_HT_SEPERATOR // not a real section, just for stats
} SectionType;
#define NUM_SEC_TYPES 6

// Section headers - unsigned is 32 bit, unsigned short is 16 bit, little endian

typedef struct  __attribute__((__packed__)) {
    unsigned char section_type;
    unsigned compressed_offset; // number of bytes from the start of the header that is the start of compressed data
    unsigned data_compressed_len;
    unsigned data_uncompressed_len;
} SectionHeader; 

// The VCF header section appears once in the file, and includes the VCF file header 
#define FILE_METADATA_LEN 72

#define VCZIP_MAGIC 0x27052012

#define COMPRESSION_ALG_BZLIB 0

typedef struct  __attribute__((__packed__)) {
    SectionHeader h;
    unsigned magic; 
    unsigned char vczip_version;
    unsigned char compression_alg;
    long long vcf_data_size; // number of bytes in the original VCF file
    long long num_lines;     // number of variants (data lines) in the original vCF file
    unsigned num_samples;             // number of samples in the original VCF file
    char created [FILE_METADATA_LEN];
    char modified[FILE_METADATA_LEN];
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
 
typedef struct  __attribute__((__packed__)) {
    SectionHeader h;
    unsigned first_line;                  // line (starting from 1) of this variant block in the VCF file
    unsigned short num_lines;             // number of variants in this block
    unsigned char phase_type;
    
    // flags
    unsigned char has_genotype_data : 1;  // 1 if there is at least one variant in the block that has FORMAT with have anything except for GT 
    unsigned char for_future_use    : 7;

    unsigned num_samples;
    unsigned num_haplotypes_per_line;     // 0 if no haplotypes
    unsigned num_sample_blocks;
    unsigned num_samples_per_block;
    unsigned short ploidy;
    unsigned vcf_data_size;               // size of variant block as it appears in the source file
    unsigned z_data_bytes;                // total bytes of this variant block in the dv file including all sections and their headers
    unsigned short haplotype_index_checksum;
    unsigned char haplotype_index[]    ;  // length is num_haplotypes. e.g. the first entry shows for the first haplotype in the original file, its index into the permuted block. # of bits per entry is roundup(log2(num_samples*ploidy)).
} SectionHeaderVariantData; 

typedef struct {
    enum {BUF_UNALLOCATED, BUF_REGULAR, BUF_OVERLAYED} type;
    const char *name; // name of allocator - used for memory debugging & statistics
    unsigned param;   // parameter provided by allocator - used for memory debugging & statistics
    unsigned size;    // number of bytes allocated to memory
    unsigned len;     // used by the buffer user according to its internal logic. not modified by malloc/realloc, zeroed by buf_free
    char *data;       // =memory+2*sizeof(long long) if buffer is allocated or NULL if not
    char *memory;     // memory allocated to this buffer - amount is: size + 2*sizeof(longlong) to allow for OVERFLOW and UNDERFLOW)
} Buffer;
#define EMPTY_BUFFER {0,0,0,0,0,0,0}

typedef struct {
#ifdef PROFILER
    long long read, piz_uncompress_variant_block, compressor, write, piz_get_variant_data_line, 
        piz_get_haplotype_data_line, piz_get_ht_permutation_lookups,
        piz_get_phase_data_line, piz_get_genotype_data_line, zfile_uncompress_section,
        piz_reconstruct_line_components, squeeze, piz_decode_pos, buffer,
        zip_compress_variant_block, segregate_all_data_lines, zip_generate_haplotype_sections, sample_haplotype_data, count_alt_alleles,
        zip_generate_genotype_sections, zip_generate_phase_sections, zip_generate_variant_data_section, gl_optimize_do, gl_optimize_undo,
        tmp1, tmp2, tmp3, tmp4, tmp5;

#else
    int dummy;
#endif
} ProfilerRec;

typedef enum {UNKNOWN, VCF, VCF_GZ, VCZ, VCZ_TEST, PIPE, STDIN, STDOUT} FileType;

typedef struct {
    void *file;
    const char *name;
    FileType type;
    // these relate to actual bytes on the disk
    long long disk_size; // 0 if not known (eg stdin)
    long long disk_so_far;  // data read/write to/from "disk" (using fread/fwrite)

    // this relate to the VCF data represented. In case of READ - only data that was picked up from the read buffer.
    long long vcf_data_size;    // VCF: size of the VCF data (if known)
                                         // VCZ: VCZ: size of original VCF data in the VCF file currently being processed
    long long vcf_data_so_far;  // VCF: data sent to/from the caller (after coming off the read buffer and/or decompression)
                                         // VCZ: VCF data so far of original VCF file currently being processed

    // Used for READING VCF/VCF_GZ files: stats used to optimize memory allocation
    double avg_header_line_len, avg_data_line_len;// average length of data line so far. 
    unsigned header_lines_so_far, data_lines_so_far; // number of lines read so far

    // Used for WRITING VCZ files
    long long disk_at_beginning_of_this_vcf_file; // the value of disk_size when starting to read this vcf file
    long long vcf_concat_data_size; // concatenated vcf_data_size of all files compressed
    long long num_lines;            // number of lines in concatenated (or single) vcf file

    // Information content stats - how many bytes does this file have in each section type
    long long section_bytes[6];   

    // USED FOR READING ALL FILES
#   define READ_BUFFER_SIZE (1<<19) // 512KB
    unsigned next_read, last_read; // indices into read_buffer
    bool eof; // we reached EOF
    char read_buffer[]; // only allocated for mode=READ files
    
} File;

// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {

    unsigned line_i;

    // initially, data from vcf file line, later segregated to components "stored" in the overlay buffers below
    Buffer line;             

    // the following 4 buffers are overlay buffers onto line, so that they don't consume memory
    Buffer variant_data;     // string terminated by a newline. len includes the newline.
    Buffer genotype_data;    // \t separated genotype data for the line. last one has \t too. no \0. exists if any variant in the variant blck has FORMAT other that "GT"
    Buffer haplotype_data;   // length=ploidy*num_samples. exists if the GT subfield exists in any variant in the variant block

    Buffer phase_data;       // used only if phase is mixed. length=num_samples. exists if haplotype data exists and ploidy>=2
    PhaseType phase_type;         // phase type of this line
    
    bool has_haplotype_data; // FORMAT field contains GT
    bool has_genotype_data;  // FORMAT field contains subfields other than GT

    unsigned gl_subfield;    // the index (1-based) of the subfield with the gt data, for GL optimizaiton. 0 means no GL subfield in the this line, or no GL optimization is needed.
} DataLine;


// IMPORTANT: if changing fields in VariantBlock, also update vb_release_vb
typedef struct variant_block_ {

    unsigned id;               // id of vb within the vb pool. compress vbs start with 100 and decompress 200

    File *vcf_file, *z_file;  // pointers to objects that span multiple VBs

    // memory management
    Buffer buffer_list;        // a buffer containing an array of pointers to all buffers allocated for this VB (either by the I/O thread or its compute thread)

    bool ready_to_dispatch;    // line data is read, and dispatcher can dispatch this VB to a compute thread
    bool in_use;               // this vb is in use
        
    DataLine data_lines[VARIANTS_PER_BLOCK];
    unsigned num_lines;        // number of lines in this variant block
    unsigned first_line;       // line number in VCF file (counting from 1), of this variant block
    unsigned variant_block_i;  // number of variant block within VCF file
    unsigned longest_line_genotype_data; // used for predictive memory allocation for gt data

    long long last_pos; // value of POS field of the previous line, to do delta encoding

    // tracking execution
    unsigned vcf_data_size;    // size of variant block as it appears in the source file

    ProfilerRec profile;

    // charactaristics of the data
    unsigned ploidy;
    unsigned num_sample_blocks;
    unsigned num_samples_per_block; // except last sample block that may have less
    unsigned num_haplotypes_per_line;
    bool has_genotype_data;    // if any variant has genotype data, then the block is considered to have it
    bool has_haplotype_data;   // ditto for haplotype data
    bool optimize_gl;          // do we apply the GL subfield optimization to this VB
    PhaseType phase_type;      // phase type of this variant block

    // working memory for segregate - we segregate a line components into these buffers, and when done
    // we copy it back to DataLine - the buffers overlaying the line field
    Buffer line_variant_data;     // string terminated by a newline. len includes the newline.
    Buffer line_genotype_data;    // \t separated genotype data for the line. last one has \t too. no \0. exists if any variant in the variant blck has FORMAT other that "GT"
    Buffer line_haplotype_data;   // length=ploidy*num_samples. exists if the GT subfield exists in any variant in the variant block
    Buffer line_phase_data;       // used only if phase is mixed. length=num_samples. exists if haplotype data exists and ploidy>=2

    // section data - ready to compress
    Buffer variant_data_section_data;    // all fields until FORMAT, newline-separated, \0-termianted. .len includes the terminating \0

    Buffer haplotype_permutation_index;
    
    // these are Buffer arrays of size vb->num_sample_blocks allocated once when used for the first time and never freed. 
    // Subsequent variant blocks that re-use the memory have the same number of samples by VCF spec. Each buffer is a string of the data as written to the VCZ file

    // Note: The sample blocks for phase and genotype data are subsequent blocks of num_samples_per_block samples.
    // in contrast the sample blocks of haplotypes are num_samples_per_block*ploidy haplotypes as permuted. 
    // e.g. genotypes and phase in sample_block_i=1 don't necessary correspond to haplotypes in sample_block_i=1.

    Buffer *haplotype_sections_data; // this is the haplotype character for each haplotype in the transposed sample group
    Buffer *phase_sections_data;     // this is the phase character for each genotype in the sample group
    Buffer *genotype_sections_data;  // this is for piz - this is \t-seperated concatenation of all genotype data for the sample block. .len is strlen() 
    Buffer genotype_one_section_data; // for zip we need only one section data

    // compresssed file data 
    Buffer z_data;                  // all headers and section data as read from disk
    Buffer z_section_headers;       // (used by piz) an array of unsigned offsets of section headers within z_data

    Buffer genotype_sample_block_line_starts_buf,   // used by zip_get_genotype_vb_start_len 
           genotype_sample_block_line_lengths_buf,
           genotype_section_lens_buf; 

    Buffer helper_index_buf;         // used by zip_do_haplotypes

    Buffer vardata_header_buf;       // used by zfile_compress_variant_data

    Buffer compressed;               // used by various zfile functions

    Buffer ht_columns_data;          // used by piz_get_ht_permutation_lookups

    Buffer next_gt_in_sample;        // used for depermuting genotype data
    Buffer gt_line_lengths;          // line length of each line of genotypes - used for depermuting genotype data

// Information content stats - how many bytes does this section have more than the corresponding part of the vcf file    
    int add_bytes[6];                
    unsigned vcf_section_bytes[6];   // how many bytes did each section have in the original vcf file - should add up to the file size
    unsigned z_section_bytes[6];        // how many bytes does each section type have (including headers) in the vcz file - should add up to the file size

#define NUM_COMPRESS_BUFS 4  // bzlib2 compress requires 4 and decompress requires 2
    Buffer compress_bufs[NUM_COMPRESS_BUFS];    // memory allocation for compressor so it doesn't do its own malloc/free
    
    // use for testing
    struct variant_block_ *test_vb;
} VariantBlock;

typedef enum {READ, WRITE} FileMode;
extern File *file_open (const char *filename, FileMode mode, FileType expected_type);
extern File *file_fdopen (int fd, FileMode mode, FileType type);
extern void file_close (File **vcf_file_p);
extern void file_remove (const char *filename);
extern bool file_has_ext (const char *filename, const char *extension);

extern bool vcffile_get_line(VariantBlock *vb, unsigned line_i_in_file, Buffer *line);
extern void vcffile_write_one_variant_block (File *vcf_file, VariantBlock *vb);
extern unsigned vcffile_write_to_disk(File *vcf_file, const Buffer *buf);
extern void vcffile_compare_pipe_to_file (FILE *from_pipe, File *vcf_file);

// reads VCF header and writes its compressed form to the VCZ file. returns num_samples.
extern bool vcf_header_vcf_to_vcz (VariantBlock *vb, bool concat_mode, unsigned *line_i, Buffer **first_data_line);
extern bool vcf_header_vcz_to_vcf (VariantBlock *vb, bool concat_mode);
extern bool vcf_header_get_vcf_header (File *z_file, SectionHeaderVCFHeader *vcf_header_header);

typedef struct {
    unsigned num_vbs;
    unsigned vb_id_prefix;
    VariantBlock vb[];
} VariantBlockPool;
extern VariantBlockPool *vb_construct_pool (unsigned num_vbs, unsigned vb_id_prefix);

extern VariantBlock *vb_get_vb(VariantBlockPool *pool, File *vcf_file, File *z_file, unsigned variant_block_i);

extern unsigned vb_num_samples_in_sb (const VariantBlock *vb, unsigned sb_i);
extern unsigned vb_num_sections(VariantBlock *vb);

extern void vb_release_vb (VariantBlock **vb_p);

extern void vb_show_progress(double *last_percent, const File *file, long long vcf_data_written_so_far,                           
                             long long bytes_compressed, // may be 0
                             unsigned seconds_so_far, bool done, bool test_mode);

extern void vb_adjust_max_threads_by_vb_size(const VariantBlock *vb, unsigned *max_threads, bool test_mode);

extern void segregate_all_data_lines (VariantBlock *vb, Buffer *lines_orig /* for testing */);

extern void zip_dispatcher (char *vcf_basename, File *vcf_file, 
                            File *z_file, bool concat_mode, bool test_mode, unsigned max_threads);

extern void piz_dispatcher (char *z_basename, File *z_file, File *vcf_file, 
                            bool concat_mode, bool test_mode, unsigned max_threads);

extern void piz_reconstruct_line_components (VariantBlock *vb);
extern void piz_merge_all_lines (VariantBlock *vb);

extern void gl_optimize_do(VariantBlock *vb, char *data, unsigned data_len, unsigned gl_subfield_no_gt);
extern void gl_optimize_undo (VariantBlock *vb, char *data, unsigned len, unsigned gl_subfield_no_gt);
extern int gl_optimize_get_gl_subfield_index(const char *data);

extern void zfile_write_vcf_header (VariantBlock *vb, Buffer *vcf_header_text);
extern void zfile_compress_variant_data (VariantBlock *vb);
extern void zfile_compress_section_data (VariantBlock *vb, SectionType section_type, Buffer *section_data);

extern bool zfile_read_one_vb (VariantBlock *vb);

// returns offset of header within data, -1 if EOF
extern int zfile_read_one_section (VariantBlock *vb, 
                                   Buffer *data /* buffer to append */, 
                                   unsigned header_size, SectionType section_type,
                                   bool allow_eof);

extern void zfile_uncompress_section(VariantBlock *vb, const void *section_header, Buffer *uncompressed_data, SectionType expected_section_type);
extern void zfile_update_vcf_header_section_header (File *z_file, long long vcf_data_size, long long vcf_num_lines);

extern void squeeze (VariantBlock *vb,
                     unsigned char *dst, // memory should be pre-allocated by caller
                     unsigned short *squeezed_checksum,
                     const unsigned *src, 
                     unsigned src_len);

extern void unsqueeze (VariantBlock *vb,
                       unsigned *normal, // memory should be pre-allocated by caller
                       const unsigned char *squeezed, 
                       unsigned short squeezed_checksum,
                       unsigned normal_len);

extern unsigned squeeze_len(unsigned int len);

extern unsigned buf_alloc (VariantBlock *vb,
                           Buffer *buf, 
                           unsigned requested_size, // whether contents of memory should be zeroed
                           float grow_at_least_factor, // grow more than new_size    
                           const char *name, unsigned param); // for debugging

extern void buf_overlay (Buffer *overlaid_buf, Buffer *regular_buf, const Buffer *copy_from,
                         unsigned *regular_buf_offset, const char *name, unsigned param);

extern void buf_free(Buffer *buf); // free buffer - without freeing memory. A future buf_malloc of this buffer will reuse the memory if possible.

static inline bool buf_is_allocated(Buffer *buf) {return buf->data != NULL && buf->type == BUF_REGULAR;}

extern void buf_copy (VariantBlock *vb, Buffer *dst, Buffer *src, 
                      unsigned start, unsigned max_size /* if 0 copies the entire buffer */);

extern void buf_test_overflows(const VariantBlock *vb);
extern bool buf_has_overflowed(const Buffer *buf);
extern bool buf_has_underflowed(const Buffer *buf);

extern long long buf_vb_memory_consumption (const VariantBlock *vb);
extern void buf_display_memory_usage(bool memory_full);

extern char *buf_human_readable_size (long long size, char *str /* out */);

// global parameters - set before any thread is created, and never change
extern unsigned    global_num_samples;
extern const char *global_cmd;         // set once in main()
extern unsigned    global_max_threads; // set in main()
extern bool        global_little_endian;  // set in main()
extern int flag_stdout, flag_force, flag_replace, flag_quiet, flag_gzip, flag_concat_mode, flag_show_content;

// unit tests
extern void segregate_data_line_unit_test();
extern void squeeze_unit_test();
extern void zip_compress_fp_unit_test();

// macros
    
#define MIN(a, b) ((a < b) ? a : b )
#define MAX(a, b) ((a > b) ? a : b )

// encode section headers in Big Endian (see https://en.wikipedia.org/wiki/Endianness)
// the reason for selecting big endian is that I am developing on little endian CPU (Intel) so
// endianity bugs will be discovered more readily this way
#define ENDN16(x) (global_little_endian ? __builtin_bswap16(x) : (x))
#define ENDN32(x) (global_little_endian ? __builtin_bswap32(x) : (x))
#define ENDN64(x) (global_little_endian ? __builtin_bswap64(x) : (x))

// sanity checks
static inline void exit_assert() { exit(1); }// an exit function so we can put a debugging break point when ASSERT exits
#define ASSERT(condition, format, ...)  { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); exit_assert(); }}
#define ASSERT0(condition, string)      { if (!(condition)) { fprintf (stderr, "\n%s\n", string); exit_assert(); }}
#define ASSERTW(condition, format, ...) { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); }}
#define ASSERTW0(condition, string)     { if (!(condition)) { fprintf (stderr, "\n%s\n", string); }
#define ABORT(format, ...)              { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); exit_assert();}
#define ABORT0(string)                  { fprintf (stderr, "\n%s\n", string); exit_assert();}

#define START_WALLCLOCK       struct timespec wallclock; clock_gettime(CLOCK_REALTIME, &wallclock); 
#define COPY_WALLCLOCK(res)   struct timespec tb; clock_gettime(CLOCK_REALTIME, &tb); \
                              res += (tb.tv_sec-wallclock.tv_sec)*1000 + (tb.tv_nsec-wallclock.tv_nsec) / 1000000; 

static inline int strcpy_tab (char *dst, const char *src)
{
    const char *start = src;

    do {
#ifdef DEBUG
        ASSERT (*src, "Error: a string has no \t terminator. src=%s", start);
#endif
        *(dst++) = *(src++);
    } while (*(src-1) != '\t');

    return *(src-1) == '\t' ? src - start : -1; // return length, or -1 if not enough dst space to complete the copy;
}

#ifdef PROFILER

#define START_TIMER     struct timespec profiler_timer; clock_gettime(CLOCK_REALTIME, &profiler_timer); 
#define COPY_TIMER(res) { struct timespec tb; \
                          clock_gettime(CLOCK_REALTIME, &tb); \
                          res += (tb.tv_sec-profiler_timer.tv_sec)*1000000000ULL + (tb.tv_nsec-profiler_timer.tv_nsec); \
                        }

extern void profiler_add (ProfilerRec *dst, const ProfilerRec *src);
extern const char *profiler_print_short (const ProfilerRec *p);
extern const char *profiler_print_report (const ProfilerRec *p);

#else
#define START_TIMER
#define COPY_TIMER(res) {}
#endif


// Windows compatibility stuff
#if __WIN32__
#define stat64  _stat64
#define fstat64 _fstat64
#else // this needs more work - there are more cases, depending if gcc is 32 or 64
#define stat64  stat
#define fstat64 fstat
#endif

#endif


