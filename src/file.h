// ------------------------------------------------------------------
//   file.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "file_types.h"
#include "context.h"
#include "aes.h"
#include "data_types.h"
#include "mutex.h"
#include "buffer.h"

#define MAX_GEN_COMP 2

typedef rom FileMode;
extern FileMode READ, WRITE, WRITEREAD;// this are pointers to static strings - so they can be compared eg "if (mode==READ)"

typedef struct File {
    void *file;
    char *name;                        // allocated by file_open(), freed by file_close()
    rom basename;                      // basename of name
    FileMode mode;
    FileSupertype supertype;            
    FileType type;
    bool is_remote;                    // true if file is downloaded from a url
    bool redirected;                   // txt_file: true if this file is redirected from stdin/stdout or a pipe
    bool is_eof;                       // we've read the entire file
    bool header_only;                  // ZIP txt_file: file has only the data-type header and no data
    DataType data_type;
    Codec source_codec;                // ZIP - txt_file: codec of txt file before redirection (eg CRAM, XZ, ZIP...)
    Codec codec;                       // ZIP - txt_file: generic codec used with this file (in PIZ we use flag.bgzf instead). If redirected - as read by txtfile (eg for cram files this is BGZF)

    // these relate to actual bytes on the disk
    int64_t disk_size;                 // 0 if not known (eg stdin or http stream). 
    int64_t disk_so_far;               // data read/write to/from "disk" (using fread/fwrite) (possibley gz/bz2 compressed)
    int64_t est_seggable_size;         // ZIP txt_file, access via txtfile_get_seggable_size(). Estimated size of txt_data in file, i.e. excluding the header. It is exact for plain files, or based on test_vb if the file has source compression

    // this relate to the textual data represented. In case of READ - only data that was picked up from the read buffer.
    int64_t txt_data_so_far_single;    // txt_file: data read (ZIP) or written (PIZ) to/from txt file so far
                                       // z_file: txt data represented in the GENOZIP data written (ZIP) or read (PIZ) to/from the genozip file so far for the current txt file
    int64_t header_size;               // txt_file ZIP: size of txt header  z_file ZIP: size of MAIN txt header
    int64_t header_size_bgzf;          // txt_file ZIP: compressed size of header - or 0 if not whole BGZF blocks
    
    int64_t txt_data_so_far_bind;      // z_file only: uncompressed txt data represented in the GENOZIP data written so far for all bound files
                                       // note: txt_data_so_far_single/bind accounts for txt modifications due to --optimize or --chain or compressing a Luft file
    int64_t txt_data_so_far_single_0;  // txt_file PIZ: accounting for data as it was in the original source file, as reading TxtHeader and VbHeader sections from the genozip file
                                       // z_file & ZIP only: same as txt_data_so_far_single/bind, but original sizes without modifications due to --chain/--optimize/Luft
    int64_t txt_data_so_far_bind_0;    // z_file & ZIP only: similar to txt_data_so_far_single_0, but for bound
    int64_t txt_disk_so_far_bind;      // z_file & ZIP only: compressed (with txt_file.codec - eg bgzf) txt data represented in the GENOZIP data written so far for all bound files
    int64_t num_lines;                 // z_file: number of lines in all txt files bound into this z_file
                                       // txt_file: number of lines, in source file terms, (read so far) in single txt file
    // per-component data (ZIP)
    int64_t disk_so_far_comp[MAX_NUM_COMPS];
    int64_t txt_data_so_far_bind_comp[MAX_NUM_COMPS];
    int64_t txt_data_so_far_bind_0_comp[MAX_NUM_COMPS];

    // Digest stuff - stored in z_file (ZIP & PIZ)
    Serializer digest_serializer;      // ZIP/PIZ: used for serializing VBs so they are MD5ed in order (not used for Adler32)
    DigestContext digest_ctx;          // ZIP/PIZ: z_file: digest context of txt file being compressed / reconstructed (used for MD5 and, in v9-13, for Adler32)
    Digest digest;                     // ZIP: z_file: digest of txt data read from input file  PIZ: z_file: as read from TxtHeader section (used for MD5 and, in v9-13, for Adler32)
    bool vb_digest_failed;             // PIZ: txt_file: At least one VB has an unexpected digest when decompressing
    
    // Used for READING & WRITING txt files - but stored in the z_file structure for zip to support bindenation (and in the txt_file structure for piz)
    uint32_t max_lines_per_vb;         // ZIP & PIZ - in ZIP, discovered while segmenting, in PIZ - given by SectionHeaderTxtHeader
    bool piz_header_init_has_run;      // PIZ: true if we called piz_header_init to initialize (only once per outputted txt_file, even if concatenated)

    // Used for READING GENOZIP files
    uint8_t genozip_version;           // GENOZIP_FILE_FORMAT_VERSION of the genozip file being read
        
    CompIType num_txt_files;          // PIZ z_file: set from genozip header 

    struct FlagsGenozipHeader z_flags; // z_file: genozip file flags as read from SectionHeaderGenozipHeader.flags
    struct FlagsTxtHeader txt_flags;   // txt_file PIZ: genozip file flags as read from SectionHeaderTxtHeader.flags

    Buffer vb_sections_index;          // PIZ z_file: an index into VB sections
    Buffer comp_sections_index;        // PIZ z_file: an index into Txt sections

    // Used for WRITING GENOZIP files
    uint64_t disk_at_beginning_of_this_txt_file;     // z_file: the value of disk_size when starting to read this txt file
    
    uint8_t txt_header_enc_padding[AES_BLOCKLEN-1];  // just so we can overwrite txt_header with encryption padding

    SectionHeaderTxtHeader txt_header_single;        // ZIP: store the txt header of single component in bound mode
                                                     // PIZ: z_file: header of COMP_MAIN
    uint8_t txt_header_enc_padding2[AES_BLOCKLEN-1]; // same

    bool z_closes_after_me;            // Z_FILE ZIP: this z_file will close after completing the current txt_file
    
    // dictionary information used for writing GENOZIP files - can be accessed only when holding mutex
    Mutex dicts_mutex;                 // this mutex protects contexts and num_contexts from concurrent adding of a new dictionary
    Did num_contexts;                  // length of populated subfield_ids and mtx_ctx;
    
    DictIdtoDidMap d2d_map; // map for quick look up of did_i from dict_id : 64K for key_map, 64K for alt_map 
    ContextArray contexts;             // a merge of dictionaries of all VBs
    Buffer ra_buf;                     // ZIP/PIZ:  RAEntry records: ZIP: of DC_PRIMARY ; PIZ - PRIMARY or LUFT depending on flag.luft
    Buffer ra_buf_luft;                // ZIP only: RAEntry records of DC_LUFT
    
    Buffer chrom2ref_map;              // ZIP/PIZ: mapping from user file chrom to alternate chrom in reference file
    Buffer ref2chrom_map;              // ZIP: an inverted chrom2ref_map, created by ref_compress_ref

    // section list - used for READING and WRITING genozip files
    Buffer section_list_buf;           // section list to be written as the payload of the genotype header section
    Buffer section_list_save;          // a copy of the section_list in case it is modified due to recon plan.
    CompIType num_txts_so_far;         // ZIP z_file: number of txt files compressed into this z_file - each becomes a component
                                       // PIZ z_file: number of txt files written from this z_file - each generated from one or more components 

    // TXT file: reading
    Buffer unconsumed_txt;             // ZIP: excess data read from the txt file - moved to the next VB

    // TXT file: stuff reading and writing txt files compressed with BGZF
    Buffer unconsumed_bgzf_blocks;     // ZIP: unconsumed or partially consumed bgzf blocks - moved to the next VB
    Buffer bgzf_isizes;                // ZIP/PIZ: uncompressed size of the bgzf blocks in which this txt file is compressed (in BGEN16)
    Buffer bgzf_starts;                // ZIP: offset in txt_file of each BGZF block
    struct FlagsBgzf bgzf_flags;       // correspond to SectionHeader.flags in SEC_BGZF
    uint8_t bgzf_signature[3];         // PIZ: 3 LSB of size of source BGZF-compressed file, as passed in SectionHeaderTxtHeader.codec_info

    // TXT file: data used in --sex, --coverage and --idxstats
    Buffer coverage;
    Buffer read_count;
    Buffer unmapped_read_count;
    
    // Z_FILE: SAM/BAM SA stuff
    Buffer sag_grps;                   // Z_FILE: an SA group is a group of alignments, including the primary aligngment
    Buffer sag_gps_index;              // Z_FILE: index z_file->sag_grps by adler32(qname)
    Buffer sag_alns;                   // Z_FILE: array of {RNAME, STRAND, POS, CIGAR, NM, MAPQ} of the alignment
    Buffer sag_qnames;                 // Z_FILE
    Buffer sag_depn_index;             // Z_FILE: SAG_BY_FLAG: uniq-sorted hash(QNAME) of all depn alignments in the file
    
    union {
    Buffer sag_cigars;                 // Z_FILE: SAG_BY_SA: compressed CIGARs
    Buffer solo_data;                  // Z_FILE: SAG_BY_SOLO: solo data
    };
    Buffer sag_seq;                    // Z_FILE: bitmap of seqs in ACGT 2bit format
    Buffer sag_qual;                   // Z_FILE: compressed QUAL
    
    uint64_t saggy_near_count, mate_line_count, prim_far_count; // Z_FILE ZIP: SAM: for stats

    // Z_FILE: Deep
    Buffer vb_start_deep_line;         // Z_FILE: ZIP/PIZ: for each SAM VB, the first prim_line_i of that VB (0-based, uint64_t)
    Buffer deep_ents;                  // Z_FILE: ZIP: entries of type ZipZDeep
                                       // Z_FILE: PIZ: an array of Buffers, one for each SAM VB, containing QNAME,SEQ,QUAL of all reconstructed lines
    union {
    Buffer deep_hash;                  // Z_FILE: ZIP: hash table - indices into deep ents 
    Buffer deep_index;                 // Z_FILE: PIZ: an array of Buffers, one for each SAM VB, each containing an array of uint32_t - one for each primary line - index into deep_ents[vb_i] of PizZDeep of that line
    };
    Mutex custom_merge_mutex;          // Z_FILE: ZIP: used to merge deep, but in the future could be used for other custom merges

    // Z_FILE: DVCF stuff
    Buffer rejects_report;             // Z_FILE ZIP --chain: human readable report about rejects
    Buffer apriori_tags;               // Z_FILE ZIP DVCF: used for INFO/FORMAT tag renaming. Data from command line options if --chain, or VCF header if DVCF
    
    // TXT_FILE: DVCF stuff
    uint8_t coords;                    // TXT_FILE ZIP: DVCF: Set from ##dual_coordinates and immutable thereafter
    uint64_t reject_bytes;             // ZIP DVCF: number of bytes in lines originating from ##primary_only/##luft_only, not yet assigned to a VB

    // Z_FILE: FASTA
    uint64_t num_sequences;            // ZIP: for stats

    // Reconstruction plan, for reconstructing in sorted order if --sort: [0] is primary coords, [1] is luft coords
    Mutex recon_plan_mutex[2];         // TXT_FILE ZIP: VCF: protect vb_info and line_info during merging of VB data
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
    
    // Information content stats - how many bytes and how many sections does this file have in each section type
    uint32_t num_vbs;                  // ZIP: z_file/txt_file PIZ: txt_file: number of VBs processed z_file: total VBs in file
    uint32_t num_vbs_dispatched;       // ZIP: txt_file
    uint32_t num_preproc_vbs_joined;   // PIZ: z_file
    uint32_t max_conc_writing_vbs;     // PIZ z_file: the maximal value conc_writing_vbs across all SEC_RECON_PLAN sections in this z_file
} File;

#define z_has_gencomp (z_file && z_file->z_flags.has_gencomp)
#define z_is_dvcf (z_file && (Z_DT(VCF) || Z_DT(BCF)) && z_has_gencomp)
#define z_sam_gencomp (z_file && (Z_DT(SAM) || Z_DT(BAM)) && z_has_gencomp) // note: is BAM file in piz are Z_DT(SAM) and in zip are Z_DT(BAM)

extern DataType last_z_dt; // data_type of last z_file opened

// methods
extern FileP file_open (rom filename, FileMode mode, FileSupertype supertype, DataType data_type /* only needed for WRITE */);
extern void file_close (FileP *file_p, bool index_txt, bool cleanup_memory /* optional */);
extern void file_write (FileP file, const void *data, unsigned len);
extern bool file_seek (FileP file, int64_t offset, int whence, int soft_fail); // SEEK_SET, SEEK_CUR or SEEK_END
extern int64_t file_tell_do (FileP file, FailType soft_fail, FUNCLINE);
#define file_tell(file,soft_fail) file_tell_do ((file), (soft_fail), __FUNCLINE) 
extern void file_set_input_type (rom type_str);
extern void file_set_input_size (rom size_str);
extern FileType file_get_type (rom filename);
extern DataType file_get_data_type (FileType ft, bool is_input);
extern Codec file_get_codec_by_txt_ft (DataType dt, FileType txt_ft, bool source);
extern rom file_plain_text_ext_of_dt (DataType dt);
extern void file_get_raw_name_and_type (char *filename, char **raw_name, FileType *ft);
extern void file_assert_ext_decompressor (void);
extern void file_kill_external_compressors (void);
extern FileType file_get_z_ft_by_txt_in_ft (DataType dt, FileType txt_ft);
extern DataType file_get_dt_by_z_ft (FileType z_ft);
extern FileType file_get_z_ft_by_dt (DataType dt);
extern rom file_plain_ext_by_dt (DataType dt);
extern void file_remove_codec_ext (char *filename, FileType ft);
extern rom ft_name (FileType ft);

// wrapper operations for operating system files
extern void file_get_file (VBlockP vb, rom filename, BufferP buf, rom buf_name, uint64_t max_size, bool verify_textual, bool add_string_terminator);
extern bool file_put_data (rom filename, const void *data, uint64_t len, mode_t mode);
extern void file_put_data_abort (void);

typedef struct { char s[64]; } PutLineFn;
extern PutLineFn file_put_line (VBlockP vb, STRp(line), rom msg);
extern bool file_exists (rom filename);
extern bool file_is_fifo (rom filename);
extern bool file_is_dir (rom filename);
extern void file_mkdir (rom dirname);
extern void file_remove (rom filename, bool fail_quietly);
extern bool file_rename (rom old_name, rom new_name, bool fail_quietly);
extern void file_mkfifo (rom filename);
extern uint64_t file_get_size (rom filename);
extern char *file_compressible_extensions (bool plain_only);
extern char *file_make_unix_filename (char *filename);

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

static inline bool file_is_read_via_ext_decompressor(ConstFileP file) { 
    return file->supertype == TXT_FILE && (file->source_codec == CODEC_XZ || file->source_codec == CODEC_ZIP || file->source_codec == CODEC_BCF || file->source_codec == CODEC_CRAM);
}

static inline bool file_is_read_via_int_decompressor(ConstFileP file) {
    return file->supertype == TXT_FILE && (file->codec == CODEC_GZ || file->codec == CODEC_BGZF || file->codec == CODEC_BZ2);
}

static inline bool file_is_written_via_ext_compressor(ConstFileP file) {
    return file->supertype == TXT_FILE && (file->codec == CODEC_BCF || file->codec == CODEC_GZ);
}

// read the contents and a newline-separated text file and split into lines - creating char **lines, unsigned *line_lens and unsigned n_lines
#define file_split_lines(fn, name, verify_textual) \
    static Buffer data = EMPTY_BUFFER; \
    ASSINP0 (!data.len, "only one instance of a " name " option can be used"); \
    file_get_file (evb, fn, &data, "file_split_lines__" name, 0, verify_textual, false); \
    \
    str_split_enforce (data.data, data.len, 0, '\n', line, true, (name)); \
    ASSINP (!line_lens[n_lines-1], "Expecting %s to end with a newline", (fn)); \
    \
    n_lines--; /* remove final empty line */\
    str_remove_CR (line); /* note: fopen with non-binary mode (no "b") doesn't help, because it won't remove Windows-created \r when running on Unix */ \
    str_nul_separate (line); \
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
