// ------------------------------------------------------------------
//   file.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "context.h"
#include "file_types.h"
#include "aes.h"
#include "mutex.h"
#include "buffer.h"

#define MAX_GEN_COMP 2

typedef rom FileMode;
extern FileMode READ, WRITE, WRITEREAD;// this are pointers to static strings - so they can be compared eg "if (mode==READ)"

// number of alignments that are deepable or non-deepable for each of these reasons
typedef enum {          NDP_FQ_READS, NDP_DEEPABLE, NDP_DEEPABLE_TRIM, NDP_NO_ENTS, NDP_MONOSEQ, NDP_MONOQUAL, NDP_MULTI_MATCH, NDP_MULTI_TRIMMED, NDP_NO_MATCH , NUM_DEEP_STATS } DeepStatsFastq; 
#define NO_DEEP_NAMES { "FQ_reads",   "Deepable",   "Dpable_trim",     "No_ents",   "Monoseq",   "Monoqual",   "Multi_match",   "MultM_trim",      "No_match" }

typedef struct File {
    void *file;
    char *name;                        // allocated by file_open_z(), freed by file_close()
    rom basename;                      // basename of name
    FileMode mode;
    FileSupertype supertype;            
    FileType type;
    bool is_remote;                    // true if file is downloaded from a url
    bool redirected;                   // TXT_FILE: true if this file is redirected from stdin/stdout or a pipe
    bool is_eof;                       // we've read the entire file
    bool is_in_tar;                    // z_file: file is embedded in tar file
    bool is_scanned;                   // TXT_FILE: sam_sag_by_flag_scan_for_depn has been performed for this file
    DataType data_type;
    Codec source_codec;                // TXT_FILE ZIP: codec of txt file before redirection (eg CRAM, XZ, ZIP...). Note: CODEC_BAM if BAM (with or without internal bgzf compression)
                                       // Z_FILE PIZ: set to CODEC_BCF or CODEC_CRAM iff GenozipHeader.data_type is DT_BCF/DT_CRAM
    Codec codec;                       // TXT_FILE ZIP: generic codec used with this file. If redirected - as read by txtfile (eg for cram files this is BGZF)

    // these relate to actual bytes on the disk
    int64_t disk_size;                 // ZIP: size of actual file on disk. 0 if not known (eg stdin or http stream). 
    int64_t disk_so_far;               // ZIP: Z/TXT_FILE: data actually read/write to/from "disk" (using fread/fwrite), (TXT_FILE: possibley bgzf/gz/bz2 compressed ; 0 if external compressor is used for reading).
    int64_t est_seggable_size;         // TXT_FILE ZIP, access via txtfile_get_seggable_size(). Estimated size of txt_data in file, i.e. excluding the header. It is exact for plain files, or based on test_vb if the file has source compression
    int64_t est_num_lines;             // TXT_FILE ZIP, an alternative for progress bar - by lines instead of bytes (used for CRAM)
    
    // this relate to the textual data represented. In case of READ - only data that was picked up from the read buffer.
    int64_t txt_data_size;             // TXT_FILE: PIZ: value copied from SectionHeaderTxtHeader.txt_data_size. At the end of piz, expected to be equal to txt_data_so_far_single if no PIZ modifications (eg filters) 
    int64_t txt_data_so_far_single;    // TXT_FILE: data read (ZIP) or written (PIZ with writer) to/from txt file so far or reconsturcted (PIZ with --test)
                                       // Z_FILE: txt data represented in the GENOZIP data written (ZIP) or read (PIZ) to/from the genozip file so far for the current txt file
    int64_t header_size;               // TXT_FILE ZIP: size of txt header z_file ZIP: size of MAIN txt header
    int64_t header_size_bgzf;          // TXT_FILE ZIP: compressed size of header - or 0 if not whole BGZF blocks
    
    int64_t txt_data_so_far_bind;      // Z_FILE only: uncompressed txt data represented in the GENOZIP data written so far for all bound files
                                       // note: txt_data_so_far_single/bind accounts for txt modifications due to e.g. --optimize 
    int64_t txt_data_so_far_single_0;  // TXT_FILE PIZ: accounting for expected recon size (so far) as reported by TxtHeader and VbHeader sections, without any piz-side modifications
                                       // Z_FILE & ZIP only: same as txt_data_so_far_single/bind, but original sizes without modifications due to e.g. --optimize
    int64_t txt_data_so_far_bind_0;    // Z_FILE & ZIP only: similar to txt_data_so_far_single_0, but for bound
    int64_t num_lines;                 // Z_FILE: number of lines in all txt files bound into this z_file (sum of comp_num_lines)
                                       // TXT_FILE: number of lines, in source file terms, (read so far) in single txt file

    // Digest stuff - stored in z_file (ZIP & PIZ)
    Serializer digest_serializer;      // ZIP/PIZ: used for serializing VBs so they are MD5ed in order (not used for Adler32)
    DigestContext digest_ctx;          // ZIP/PIZ: Z_FILE: digest context of txt file being compressed / reconstructed (used for MD5 and, in v9-13, for Adler32, starting v15 also for make-reference)
    DigestContext v13_commulative_digest_ctx; // PIZ: z_file: used for multi-component up-to-v13 files - VB digests (adler and md5) are commulative since the beginning of the data, while txt file digest are commulative only with in the component.
    Digest digest;                     // ZIP: Z_FILE: digest of txt data read from input file (make-ref since v15: digest of in-memory genome)  PIZ: z_file: as read from TxtHeader section (used for MD5 and, in v9-13, for Adler32)
    bool vb_digest_failed;             // PIZ: TXT_FILE: At least one VB has an unexpected digest when decompressing
    
    // PIZ: reference file name and digest as appears in the z_file header 
    char ref_filename_used_in_zip[REF_FILENAME_LEN]; // PIZ: Z_FILE: ref filename as appears in the z_file header - this is not necessarily the file used - it could be overridden with --reference or $GENOZIP_REFERENCE
    Digest ref_genome_digest;          // PIZ: Z_FILE: digest of REF_EXTERNAL file with which this file was compressed
    
    // Used for READING & WRITING txt files - but stored in the z_file structure for zip to support bindenation (and in the txt_file structure for piz)
    uint32_t max_lines_per_vb;         // ZIP & PIZ - in ZIP, discovered while segmenting, in PIZ - given by SectionHeaderTxtHeader
    bool piz_header_init_has_run;      // PIZ: true if we called piz_header_init to initialize (only once per outputted txt_file, even if concatenated)

    // Used for READING GENOZIP files
    uint8_t genozip_version;           // major version of the genozip file being created / read
    uint16_t genozip_minor_ver;
            
    CompIType num_txt_files;           // PIZ Z_FILE: set from genozip header (ZIP: see num_txts_so_far)
    CompIType num_txts_so_far;         // ZIP Z_FILE: number of txt files compressed into this z_file - each becomes one or more components
                                       // PIZ Z_FILE: number of txt files written from this z_file - each generated from one or more components 
    char txt_filename[TXT_FILENAME_LEN];// PIZ Z_FILE: name of txt filename (same length as in SectionHeaderTxtHeader)

    union {
    struct FlagsGenozipHeader z_flags; // Z_FILE PIZ: genozip file flags as read from SectionHeaderGenozipHeader.flags
    struct FlagsTxtHeader txt_flags;   // TXT_FILE PIZ: genozip file flags as read from SectionHeaderTxtHeader.flags
    };                                 

    Buffer vb_sections_index;          // PIZ Z_FILE: an index into VB sections
    Buffer comp_sections_index;        // PIZ Z_FILE: an index into Txt sections

    // Used for WRITING GENOZIP files
    Mutex zriter_mutex;                // Z_FILE ZIP: forces one writer at a time
    Buffer zriter_threads;             // Z_FILE ZIP: list of currently writing threads and their buffers
    bool zriter_last_was_fg;

    #define Z_READ_FP(f) (((f)->mode == READ || flag.no_zriter) ? (FILE *)(f)->file : (f)->z_reread_file) 
    #define GET_FP(f,md) (((f)->mode == READ || (md) != READ || (f)->supertype != Z_FILE || flag.no_zriter) ? (FILE *)(f)->file : (f)->z_reread_file)

    FILE *z_reread_file;               // Z_FILE ZIP: file pointer for re-reading files as they are being writting (used for pair reading)

    SectionHeaderTxtHeader txt_header_hdr; // Z_FILE ZIP: store the txt header of single component in bound mode
    bool z_closes_after_me;            // Z_FILE ZIP: this z_file will close after completing the current txt_file
    
    // dictionary information used for writing GENOZIP files - can be accessed only when holding mutex
    Mutex dicts_mutex;                 // this mutex protects contexts and num_contexts from concurrent adding of a new dictionary
    Did num_contexts;                  // length of populated subfield_ids and mtx_ctx;
    
    DictIdtoDidMap d2d_map; // map for quick look up of did_i from dict_id : 64K for key_map, 64K for alt_map 
    ContextArray contexts;             // Z_FILE ZIP/PIZ: a merge of dictionaries of all VBs
    Buffer ra_buf;                     // ZIP/PIZ:  RAEntry records
    
    // section list - used for READING and WRITING genozip files
    Buffer section_list_buf;           // section list to be written as the payload of the genotype header section
    Buffer section_list_save;          // a copy of the section_list in case it is modified due to recon plan.

    // TXT file: reading
    Buffer unconsumed_txt;             // ZIP: excess uncompressed data read from the txt file - moved to the next VB

    // TXT file: BGZF stuff reading and writing compressed txt files 
    Buffer unconsumed_bgzf_blocks;  // ZIP TXT BGZF: unconsumed or partially consumed bgzf blocks - moved to the next VB
    Buffer bgzf_isizes;             // ZIP/PIZ: uncompressed size of the bgzf blocks in which this txt file is compressed (in BGEN16)
    Buffer bgzf_starts;             // ZIP: offset in txt_file of each BGZF block
    Buffer bgzf_plausible_levels;   // ZIP: discovering library/level. .count = number of BGZF blocks tested so far
    struct FlagsBgzf bgzf_flags;    // corresponds to SectionHeader.flags in SEC_BGZF
    uint8_t bgzf_signature[3];      // PIZ: 3 LSB of size of source BGZF-compressed file, as passed in SectionHeaderTxtHeader.codec_info

    // TXT file: IGZIP stuff reading and writing compressed txt files 
    Buffer igzip_data;              // ZIP TXT GZ (with igzip): yet-uncompressed data read from disk
    Buffer igzip_state;             // ZIP TXT GZ (with igzip)
    uint64_t gzip_start_Ltxt;       // ZIP TXT GZ: Ltxt at the beginning for this gzip_section

    // TXT FILE: accounting for truncation when --truncate-partial-last-line is used
    bool bgzf_truncated_last_block;    // ZIP: detected a truncated last block
    uint32_t last_truncated_line_len;  // ZIP: bytes truncated due to incomplete final line. note that if file is BGZF, then this truncated data is contained in the final intact BGZF blocks, after already discarding the final incomplete BGZF block

    // TXT file: data used in --sex, --coverage and --idxstats
    Buffer coverage;
    Buffer read_count;
    Buffer unmapped_read_count;
    
    // Z_FILE: SAM/BAM SA stuff
    Buffer sag_grps;                   // Z_FILE: an SA group is a group of alignments, including the primary aligngment
    Buffer sag_grps_index;             // Z_FILE: index z_file->sag_grps by adler32(qname)
    Buffer sag_alns;                   // Z_FILE: array of {RNAME, STRAND, POS, CIGAR, NM, MAPQ} of the alignment
    Buffer sag_qnames;                 // Z_FILE
    Buffer sag_depn_index;             // Z_FILE: SAG_BY_FLAG: uniq-sorted hash(QNAME) of all depn alignments in the file
    
    union {
    Buffer sag_cigars;                 // Z_FILE: SAG_BY_SA: compressed CIGARs
    Buffer solo_data;                  // Z_FILE: SAG_BY_SOLO: solo data
    };
    Buffer sag_seq;                    // Z_FILE: bitmap of seqs in ACGT 2bit format
    Buffer sag_qual;                   // Z_FILE: compressed QUAL
    
    // Z_FILE: Deep
    Buffer vb_start_deep_line;         // Z_FILE: ZIP/PIZ: for each SAM VB, the first deepable_line_i of that VB (0-based, uint64_t)
    Buffer deep_ents;                  // Z_FILE: ZIP: entries of type ZipZDeep
                                       // Z_FILE: PIZ: an array of Buffers, one for each SAM VB, containing QNAME,SEQ,QUAL of all reconstructed lines
    
    Buffer deep_index;                 // Z_FILE: PIZ: an array of Buffers, one for each SAM VB, each containing an array of uint32_t - one for each primary line - index into deep_ents[vb_i] of PizZDeep of that line

    #define BY_SEQ   0
    #define BY_QNAME 1
    Buffer deep_index_by[2];           // Z_FILE: ZIP: hash table (BY_SEQ,BY_QNAME) - indices into deep ents, indexed by a subset of the hash(SEQ) bits 

    Mutex qname_huf_mutex;             // Z_FILE: PIZ: Deep: mutex to protect qname_huf
    HuffmanP qname_huf;                // Z_FILE: PIZ: Deep: Used for compressing QNAME in deep_ents
    char master_qname[SAM_MAX_QNAME_LEN+1]; // Z_FILE: PIZ: Deep: 
    int qnames_sampled;                // Z_FILE: PIZ: Deep: Number of QNAMEs sampled for producing the huffman compressor

    // Reconstruction plan, for reconstructing in sorted order if --sort: [0] is primary coords, [1] is luft coords
    Mutex recon_plan_mutex;            // TXT_FILE ZIP: VCF: protect vb_info and line_info during merging of VB data
    Buffer vb_info[2];                 // TXT_FILE ZIP: VCF: array of ZipVbInfo per VB, indexed by (vb_i-1), 0:PRIMARY, 1:LUFT
                                       // Z_FILE   ZIP: SAM: array of SamGcVbInfo for: 0:PRIM 1:DEPN
                                       // Z_FILE   PIZ: [0]: used by writer [1]: used to load SAM SA Groups - array of PlsgVbInfo - entry per PRIM vb 
    Buffer line_info[2];               // TXT_FILE ZIP: VCF: array of LineInfo per line or gapless range in txt_file
                                       //               SAM: array of uint32 - lengths of lines in PRIM/DEPN
    Buffer recon_plan;                 // TXT_FILE ZIP/PIZ: array of ReconPlanItem - order of reconstruction of ranges of lines, to achieve a sorted file. VCF: [0]=PRIM rendition [1]=LUFT rendition
                                       // Z_FILE   PIZ: plan for entire z_file, txt_file.recon_plan is assigned a portion of this plan
    Buffer recon_plan_index;           // TXT_FILE ZIP / Z_FILE PIZ: An array of BufWord, one for each VB: start and length of VB in recon_plan
    Buffer txt_header_info;            // Z_FILE   PIZ: used by writer
    uint64_t lines_written_so_far;     // TXT_FILE PIZ: number of textual lines (excluding the header) that passed all filters except downsampling, and is to be written to txt_file, or downsampled-out
    uint32_t num_vbs_verified;         // Z_FILE   PIZ: number of VBs verified so far
    
    // Z_FILE 
    Buffer bound_txt_names;            // ZIP: Stats data: a concatenation of all bound txt_names that contributed to this genozip file
    Buffer aliases;                    // ZIP/PIZ
    struct timespec start_time;        // Z_FILE: For stats: time z_file object was created in memory 
    Mutex ctx_mutex[MAX_DICTS];        // Z_FILE ZIP: Context z_file (only) is protected by a mutex 
    Mutex custom_merge_mutex;          // Z_FILE: ZIP: used to merge deep, but in the future could be used for other custom merges

    // Information content stats 
    CompIType num_components;          // ZIP/PIZ z_file: number of components in this file (inc. generated components)
    uint32_t num_vbs;                  // ZIP: z_file/txt_file PIZ: txt_file: number of VBs processed z_file: total VBs in file
    uint32_t num_vbs_dispatched;       // ZIP: txt_file
    uint32_t num_preproc_vbs_joined;   // Z_FILE: PIZ
    uint32_t max_conc_writing_vbs;     // Z_FILE: PIZ: the maximal value conc_writing_vbs across all SEC_RECON_PLAN sections in this z_file
    uint64_t deep_stats[NUM_DEEP_STATS]; // Z_FILE: ZIP/PIZ: SAM: stats collection on Deep performance
    uint64_t saggy_near_count, mate_line_count, prim_far_count; // Z_FILE ZIP: SAM: for stats
    union {
    uint64_t num_sequences;            // Z_FILE: ZIP: FASTA: num "DESC" lines in this file. for stats
    uint64_t num_perfect_matches;      // Z_FILE: ZIP: SAM/BAM/FASTQ: number of perfect matches found by aligner. for stats 
    Ploidy max_ploidy;                 // Z_FILE: ZIP: VCF 
    Ploidy max_ploidy_for_mux;         // Z_FILE: PIZ: VCF: copied from SectionHeaderGenozipHeader.max_ploidy_for_mux
    };
    uint64_t num_aligned;              // Z_FILE: ZIP: SAM/BAM/FASTQ: number of alignments successfully found by aligner. for stats 
    uint64_t domq_lines[MAX_NUM_COMPS];// Z_FILE: ZIP: number of lines per components of the various QUAL codecs (show in codec_qual_show_stats)
    uint64_t divr_lines[MAX_NUM_COMPS];
    uint64_t homp_lines[MAX_NUM_COMPS];
    uint64_t pacb_lines[MAX_NUM_COMPS];
    uint64_t longr_lines[MAX_NUM_COMPS];
    uint64_t normq_lines[MAX_NUM_COMPS];

    // per-component data (ZIP)
    uint64_t comp_num_lines[MAX_NUM_COMPS];             // Z_FILE: PIZ/ZIP: number of lines in each component 
    int64_t txt_file_disk_sizes_sum;                    // Z_FILE ZIP: sum of txt_file_disk_sizes[*]
    int64_t txt_file_disk_sizes[MAX_NUM_COMPS];         // Z_FILE ZIP: size of original txt file (whether or not compressed) of each component. 
    int64_t disk_so_far_comp[MAX_NUM_COMPS];            // Z_FILE ZIP: per-component size if z_file VB sections (note: global area sections, including SEC_DICT, are not accounted for)
    int64_t txt_data_so_far_bind_comp[MAX_NUM_COMPS];   // Z_FILE ZIP: per-component txt_size after modifications (due to --optimzie etc)
    int64_t txt_data_so_far_bind_0_comp[MAX_NUM_COMPS]; // Z_FILE ZIP: per-component txt_size before modifications 
    Codec comp_codec[MAX_NUM_COMPS];                    // Z_FILE ZIP: codec used for every txt file component (i.e. excluding generated components)
    Codec comp_source_codec[MAX_NUM_COMPS];             // Z_FILE ZIP: source codec used for every txt file component (i.e. excluding generated components)
    FlagsBgzf comp_bgzf[MAX_NUM_COMPS];                 // Z_FILE ZIP BGZF: library and level of BGZF of each comp
    bool gzip_section_size_single_block[MAX_NUM_COMPS]; // Z_FILE ZIP GZ: if true, disregard gzip_section_size as the entire file is just one GZ block
    uint64_t gzip_section_size[MAX_NUM_COMPS];          // Z_FILE ZIP GZ: size of one gzip section (in uncompressed terms) in case of concatenated gzip. -1 if they are not equal size.
} File;

#define z_has_gencomp (z_file && z_file->z_flags.has_gencomp) // ZIP/PIZ
#define z_sam_gencomp (z_file && (Z_DT(SAM) || Z_DT(BAM)) && z_has_gencomp) // note: is BAM file in piz are Z_DT(SAM) and in zip are Z_DT(BAM)

// methods
extern FileP file_open_z_read (rom filename);
extern FileP file_open_z_write (rom filename, FileMode mode, DataType data_type, Codec source_codec);
extern StrText file_get_z_run_time (FileP file);
extern FileP file_open_txt_read (rom filename);
extern FileP file_open_txt_write (rom filename, DataType data_type, BgzfLevel bgzf_level);
extern void file_close (FileP *file_p);
extern void file_write_txt (const void *data, unsigned len);
extern bool file_seek (FileP file, int64_t offset, int whence, rom mode, FailType soft_fail); // SEEK_SET, SEEK_CUR or SEEK_END
extern int64_t file_tell_do (FileP file, FailType soft_fail, FUNCLINE);
#define file_tell(file,soft_fail) file_tell_do ((file), (soft_fail), __FUNCLINE) 
extern FileType file_get_type (rom filename);
extern DataType file_get_data_type_of_input_file (FileType ft);
extern DataType file_piz_get_dt_of_out_filename (void);
extern Codec file_get_codec_by_txt_ft (DataType dt, FileType txt_ft, bool source);
extern uint32_t file_get_raw_name_and_type (rom filename, rom *raw_name, FileType *ft);
extern void file_assert_ext_decompressor (void);
extern void file_kill_external_compressors (void);
extern FileType file_get_z_ft_by_txt_in_ft (DataType dt, FileType txt_ft);
extern DataType file_get_dt_by_z_filename (rom z_filename);
extern FileType file_get_default_z_ft_of_data_type (DataType dt);
extern rom file_plain_ext_by_dt (DataType dt);
extern rom ft_name (FileType ft);

// wrapper operations for operating system files
typedef enum { VERIFY_NONE, VERIFY_ASCII, VERIFY_UTF8 } FileContentVerificationType; 
extern void file_get_file (VBlockP vb, rom filename, BufferP buf, rom buf_name, uint64_t max_size, FileContentVerificationType ver_type, bool add_string_terminator);

extern bool file_put_data (rom filename, const void *data, uint64_t len, mode_t mode);
extern void file_put_data_abort (void);
extern void file_put_data_reset_after_fork (void);

typedef struct { char s[64]; } PutLineFn;
extern PutLineFn file_put_line (VBlockP vb, STRp(line), rom msg);
extern bool file_exists (rom filename);
extern bool file_is_fifo (rom filename);
extern bool file_is_dir (rom filename);
extern void file_mkdir (rom dirname);
extern void file_remove (rom filename, bool fail_quietly);
extern bool file_rename (rom old_name, rom new_name, bool fail_quietly);
extern void file_gzip (char *filename);
extern void file_mkfifo (rom filename);
extern uint64_t file_get_size (rom filename);
extern StrTextLong file_compressible_extensions (bool plain_only);
extern bool file_buf_locate (FileP file, ConstBufferP buf);

#define FILENAME_STDIN  "(stdin)"
#define FILENAME_STDOUT "(stdout)"

#define file_printname(file) (!file             ? "(none)"       \
                           : (file)->name       ? (file)->name   \
                           : (file)->mode==READ ? FILENAME_STDIN \
                           :                      FILENAME_STDOUT)
#define txt_name file_printname(txt_file)
#define z_name   file_printname(z_file)

#define CLOSE(fd,name,quiet) ({ ASSERTW (!close (fd) || (quiet),  "Warning in %s:%u: Failed to close %s: %s",  __FUNCLINE, (name), strerror(errno));})
#define FCLOSE(fp,name) ({ if (fp) { ASSERTW (!fclose (fp), "Warning in %s:%u: Failed to fclose %s: %s", __FUNCLINE, (name), strerror(errno)); fp = NULL; } })
 
 // ---------------------------
// tests for compression types
// ---------------------------

#define SRC_CODEC(x) (txt_file->source_codec == CODEC_##x)

#define SC(x) (file->source_codec == CODEC_##x)
static inline bool is_read_via_ext_decompressor(ConstFileP file) { return SC(XZ)|| SC(ZIP) || SC(BCF)|| SC(CRAM) || SC(ORA); }
#undef SC

#define FC(x) (codec == CODEC_##x)
static inline bool is_read_via_int_decompressor(Codec codec)  { return FC(GZ)|| FC(BGZF)|| FC(BZ2); }
static inline bool is_written_via_ext_compressor(Codec codec) { return FC(BCF) || FC(CRAM); }
#undef FC

// read the contents and a newline-separated text file and split into lines - creating char **lines, unsigned *line_lens and unsigned n_lines
#define file_split_lines(fn, name, ver_type)                                                \
    static Buffer data = {};                                                                \
    ASSINP0 (!data.len, "only one instance of a " name " option can be used");              \
    file_get_file (evb, fn, &data, "file_split_lines__" name, 0, ver_type, false);          \
    ASSINP (data.len, "File %s is empty", (fn));                                            \
    ASSINP (*BLSTc(data) == '\n', "File %s expected to end with a newline", (fn));          \
                                                                                            \
    str_split_enforce (data.data, data.len, 0, '\n', line, true, (name));                   \
    ASSINP (!line_lens[n_lines-1], "Expecting %s to end with a newline", (fn));             \
                                                                                            \
    n_lines--; /* remove final empty line */                                                \
    str_remove_CR (line); /* note: fopen with non-binary mode (no "b") doesn't help, because it won't remove Windows-created \r when running on Unix */ \
    str_nul_separate (line);                                                                \
    /* note: its up to the caller to free data */

// platform compatibility stuff
#ifdef _WIN32
#define stat64  _stat64
#define fstat64 _fstat64
#else // this needs more work - there are more cases, depending if gcc is 32 or 64
#define stat64  stat
#define fstat64 fstat
#endif

#ifdef __APPLE__
#define fseeko64 fseeko
#define ftello64 ftello
#endif
