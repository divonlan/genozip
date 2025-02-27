// ------------------------------------------------------------------
//   file.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
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

typedef rom FileMode;
extern FileMode READ, WRITE, WRITEREAD;// this are pointers to static strings - so they can be compared eg "if (mode==READ)"

//                             number of sam alignments that are non-deepable for each of these reasons                             | number of fastq reads that are deepable or non-deepable for each of these reasons
typedef enum {                 RSN_DEEPABLE, RSN_SECONDARY, RSN_SUPPLEMENTARY, RSN_NO_SEQ, RSN_CONSENSUS, RSN_OVERFLOW, NDP_FQ_READS=16,           NDP_DEEPABLE, NDP_DEEPABLE_TRIM, NDP_NO_ENTS, NDP_BAD_N_QUAL, NDP_SAM_DUP, NDP_SAM_DUP_TRIM, NDP_NO_MATCH, NUM_DEEP_STATS_ZIP } DeepStatsZip; 
#define DEEP_STATS_NAMES_ZIP { "Deepable",   "Secondary",   "Supplementary",   "No_SEQ",   "Consensus",   "Overflow",   [NDP_FQ_READS]="FQ_reads", "Deepable",   "Dpable_trim",     "No_ents",   "Bad_N_qual",   "SAM_dup",   "SAM_dup_trm",    "No_match"  }

// bytes consumed by QNAME/SEQ/QUAL in z_file->deep_ents (SEQ_* in same order as PizZDeepSeqEncoding)                                                  | reasons why SEQ is stored explicitly in z_file->deep_ents and not vs the reference
typedef enum {                 QNAME_BYTES, SEQ_PACKED_BYTES, SEQ_BY_REF_BYTES, SEQ_ARITH_BYTES, QUAL_MONOCHAR_BYTES, QUAL_BYTES,  EXPL_DEEPABLE,    EXPL_SEQ_AS_REF, EXPL_TOO_MANY_MISMATCHES, EXPL_SEQ_NONREF_NON_ACGT, EXPL_SEQ_END_OF_CONTIG, EXPL_SEQ_COPY_VERBATIM, EXPL_BISULFITE, NUM_DEEP_STATS_PIZ } DeepStatsPiz;
#define DEEP_STATS_NAMES_PIZ { "QNAME",     "SEQ_packed",     "SEQ_by_ref",     "SEQ_arith",     "QUAL_monochar",     "QUAL_huff", "Total_deepable", "Stored_as_ref", "Too_many_mismatches",    "nonref_has_non-ACGT",    "End_of_contig",        "Copy_verbatim",        "Bisulfite"}

#define NUM_DEEP_STATS ((int)NUM_DEEP_STATS_ZIP > (int)NUM_DEEP_STATS_PIZ ? (int)NUM_DEEP_STATS_ZIP : (int)NUM_DEEP_STATS_PIZ)

//                             stats regarding sam alignments                                                                               
typedef enum {                 BA_TOTAL,     BA_USABLE, BA_UNMAPPED, BA_NO_SEQ, BA_NOT_PRIMARY, BA_NOT_IN_REF, BA_OVERFLOW, BA_AFTER } BamassStats; 
#define BAMASS_STATS_NAMES   { "Alignments", "Usable",  "Unmapped",  "No_SEQ",  "Not_primary",  "Not_in_ref",  "Overflow",  ""       }     

#define LIBDEFLATE_MAX_LEVEL 12
#define ZLIB_MAX_LEVEL 9
#define NUM_PLAUSIBLE_LEVELS ((LIBDEFLATE_MAX_LEVEL+1)+LIBDEFLATE_MAX_LEVEL+ZLIB_MAX_LEVEL)

typedef struct File {
    void *file;
    char *name;                        // allocated by file_open_z(), freed by file_close()
    rom basename;                      // basename of name
    FileMode mode;
    FileSupertype supertype;            
    FileType type;
    bool is_remote;                    // true if file is downloaded from a url
    bool redirected;                   // TXT_FILE: true if this file is redirected from stdin/stdout or a pipe
    bool no_more_blocks;               // ZIP TXT_FILE: txtfile_read_block has completed returning all the file's data (note: it is possible that we read all the data on the disk file so feof(fp)=true, but no_more_blocks=false, because some of it is waiting in buffers: gz_data or state->avail_in)
    bool is_in_tar;                    // z_file: file is embedded in tar file
    DataType data_type;
    Codec src_codec;                   // TXT_FILE ZIP/PIZ: internal or external codec of txt file (eg CRAM, BAM, XZ, ZIP, BGZF, MGZF, NONE...). Passed in SectionHeaderTxtHeader.src_codec.
                                       // Z_FILE PIZ: set to CODEC_BCF or CODEC_CRAM iff GenozipHeader.data_type is DT_BCF/DT_CRAM
    Codec effective_codec;             // TXT_FILE ZIP: method with which we actually uncompress txt_file: can be different than .codec, e.g. using GZ instead of a specialized gzip codec, or using BGZF/NONE when streaming
                                       // TXT_FILE PIZ: reconstruction codec (possibile values: NONE, BGZF, BCF, CRAM)
    uint32_t num_EOF_blocks;           // TXT_FILE ZIP MGZIP: number of EOF blocks encountered

    // these relate to actual bytes on the disk
    int64_t disk_size;                 // ZIP: size of actual file on disk. 0 if not known (eg stdin or http stream). 
    int64_t disk_so_far;               // ZIP: Z/TXT_FILE: data actually read/write to/from "disk" (using fread/fwrite), (TXT_FILE: possibley bgzf/gz/bz2 compressed ; 0 if external compressor is used for reading).
    int64_t disk_gz_uncomp_or_trunc;   // ZIP: TXT_FILE: gz-compressed data actually either decompressed or discarded due to truncate
    int64_t gz_blocks_so_far;          // ZIP: TXT_FILE: number of gz blocks read from disk
    int64_t est_seggable_size;         // TXT_FILE ZIP, Estimated size of txt_data in file, i.e. excluding the header. It is exact for plain files, or based on test_vb if the file has source compression
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

    Buffer vb_sections_index;          // ZIP/PIZ Z_FILE: an index into VB sections
    Buffer comp_sections_index;        // ZIP/PIZ Z_FILE: an index into Txt sections
    Buffer first_last_sec_index;       // ZIP/PIZ Z_FILE: sec_i of first section in the file of each section type

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
    Buffer section_list;               // Z_FILE ZIP/PIZ section list (payload of the GenozipHeader section)
    Buffer piz_reading_list;           // Z_FILE PIZ: a subset of z_file->section_list - only the TXT_HEADER and VB_HEADER sections of txt file which need to be read by PIZ. Accessed only by the main thread.
    
    // TXT file: reading
    Buffer unconsumed_txt;             // ZIP: excess uncompressed data read from the txt file - moved to the next VB: the final part of vb->txt_data that was not consumed
    Buffer unconsumed_mgzip_blocks;    // ZIP TXT MGZIP codecs: unconsumed or partially consumed MGZIP blocks - moved to the next VB
    Buffer gz_data;                    // ZIP TXT GZ: yet-unconsumed gz data read from disk. .comp_len/.uncomp_len refer to the first block in the buffer (in IL1M, but not BGZF, there might be additional data after the first block) 
    Buffer igzip_state;                // ZIP TXT GZ (with igzip).
    uint64_t start_gz_block;           // ZIP TXT GZ (with igzip): offset in file of start of a gz block
    uint32_t mgsp_vb_isize;            // ZIP TXT GZ (with MGSP): isize of each gz block in the VB (except for the last block that might be slightly more)
    uint32_t num_mgsp_blocks_in_vb;    // ZIP TXT GZ (with MGSP): number of gz blocks in this VB
    uint32_t max_mgsp_blocks_in_vb;    // ZIP TXT GZ (with MGSP): max gz blocks in a VB so far
    uint32_t max_mgzip_isize;          // ZIP TXT GZ: largest MGZIP gz block 
    bool discover_during_segconf;      // ZIP TXT GZ: gz discovery during segconf instead of file_open: for FASTQ files
    Codec vb_header_codec;             // ZIP TXT: codec used for compressing VB_HEADER payload

    // TXT_FILE: MGZIP stuff reading and writing compressed txt files 
    Mutex bgzf_discovery_mutex;        // TXT_FILE: ZIP: used to discover BGZF level
    Buffer mgzip_isizes;               // ZIP/PIZ: MGZIP: uncompressed size of the MGZIP blocks in which this txt file is compressed
    Buffer mgzip_starts;               // ZIP: offset in txt_file of each BGZF block
    struct FlagsMgzip bgzf_plausible_levels[NUM_PLAUSIBLE_LEVELS];      // ZIP: discovering library/level. .count = number of BGZF blocks tested so far
    uint8_t num_plausible_levels;       // ZIP: number of bgzf_plausible_levels that have not been marked as implausible
    uint8_t num_bgzf_blocks_tested_for_level; // ZIP: number of bgzf blocks tested for discovering level
    struct FlagsMgzip mgzip_flags;     // corresponds to SectionHeader.flags in SEC_MGZIP
    uint8_t bgzf_signature[3];         // PIZ: 3 LSB of size of source BGZF-compressed file, as passed in SectionHeaderTxtHeader.codec_info
    bool non_EOF_zero_block_found;     // ZIP: file contains an isize=0 block that is not identical to the EOF block, therefore we won't be able to reconstruct exactly
     
    // TXT_FILE: accounting for truncation when --truncate-partial-last-line is used
    uint32_t last_truncated_line_len;  // ZIP: bytes truncated due to incomplete final line. note that if file is BGZF, then this truncated data is contained in the final intact BGZF blocks, after already discarding the final incomplete BGZF block

    // TXT_FILE: data used in --sex, --coverage and --idxstats
    Buffer coverage;
    Buffer read_count;
    Buffer unmapped_read_count;
    
    // Z_FILE: SAM/BAM SAG stuff
    Buffer sag_grps;                   // Z_FILE ZIP/PIZ: an SA group is a group of alignments, including the primary aligngment
    Buffer sag_grps_index;             // Z_FILE ZIP: index z_file->sag_grps by qname hash
    Buffer sag_alns;                   // Z_FILE ZIP/PIZ: array of {RNAME, STRAND, POS, CIGAR, NM, MAPQ} of the alignment
    Buffer sag_qnames;                 // Z_FILE ZIP/PIZ: compressed qnames
    Buffer sag_depn_index;             // Z_FILE ZIP: SAG_BY_FLAG: uniq-sorted hash(QNAME) of all depn alignments in the file
    uint32_t SA_CIGAR_chewing_vb_i;    // Z_FILE ZIP: vb_i of the PRIM VB that will be chewing cigars (the first prim vb)

    union {
    Buffer sag_cigars;                 // Z_FILE ZIP: SAG_BY_SA: compressed CIGARs
    Buffer sag_solo_data;              // Z_FILE ZIP/PIZ: SAG_BY_SOLO: compressed solo data
    };
    Buffer sag_seq;                    // Z_FILE ZIP/PIZ: bitmap of seqs in ACGT 2bit format
    Buffer sag_qual;                   // Z_FILE ZIP/PIZ: compressed QUAL
    
    // Z_FILE: Deep
    Buffer vb_start_deep_line;         // Z_FILE: ZIP/PIZ: for each SAM VB, the first deepable_line_i of that VB (0-based, uint64_t)
    Buffer vb_num_deep_lines;          // Z_FILE: ZIP: for each SAM VB, the number of deepable lines
    union {
    Buffer deep_ents;                  // Z_FILE: ZIP: entries of type ZipZDeep
                                       // Z_FILE: PIZ: an array of Buffers, one for each SAM VB, containing QNAME,SEQ,QUAL of all reconstructed lines
    };
    union {
    Buffer deep_index;                 // Z_FILE: ZIP: hash table  - indices into deep_ents, indexed by a subset of the hash.qname bits 
                                       // Z_FILE: PIZ: an array of Buffers, one for each SAM VB, each containing an array of uint32_t - one for each primary line - index into deep_ents[vb_i] of PizZDeep of that line
    };
    // Reconstruction plan, for reconstructing in sorted order if --sort: [0] is primary coords, [1] is luft coords
    Buffer vb_info[3];                 // Z_FILE   ZIP: SAM: array of SamMainVbInfo for MAIN and SamGcVbInfo for PRIM, DEPN
                                       // Z_FILE   PIZ: [0]: used by writer [1]: used to load SAM SA Groups - array of PlsgVbInfo - entry per PRIM vb 
    Buffer recon_plan;                 // Z_FILE   PIZ: array of ReconPlanItem - reconstruction plan for the current txt file 
    Buffer txt_header_info;            // Z_FILE   PIZ: used by writer
    uint64_t lines_written_so_far;     // TXT_FILE PIZ: number of textual lines (excluding the header) that passed all filters except downsampling, and is to be written to txt_file, or downsampled-out
    uint32_t num_vbs_verified;         // Z_FILE   PIZ: number of VBs verified so far
    
    // Z_FILE 
    Buffer bound_txt_names;            // ZIP: Stats data: a concatenation of all bound txt_names that contributed to this genozip file
    Buffer aliases;                    // ZIP/PIZ
    struct timespec start_time;        // Z_FILE: For stats: time z_file object was created in memory 
    Mutex ctx_mutex[MAX_DICTS];        // Z_FILE ZIP: Context z_file (only) is protected by a mutex 
    Mutex custom_merge_mutex;          // Z_FILE: ZIP: used to merge deep, but in the future could be used for other custom merges
    Mutex test_abbrev_mutex;           // Z_FILE: ZIP: used to test if CIGAR_SA is abbreviated
    Buffer R1_txt_data_lens;           // Z_FILE: ZIP: FASTQ GZ: info regarding R1 VBs: txt_data.len32 of each VB 
    Buffer R1_last_qname_index;        // Z_FILE: ZIP: FASTQ GZ: info regarding R1 VBs: last qname of each VB, canonical form, nul-separated: index into R1_last_qname. Note: only accessible in main thread (bc may realloc)
    Buffer R1_last_qname;              // Z_FILE: ZIP: FASTQ GZ: data of R1_last_qname
    VBIType R1_first_vb_i;             // Z_FILE: ZIP: first vb_i of R1. Always 1 for --pair, more for --deep.

    // Information content stats 
    CompIType num_components;          // ZIP/PIZ z_file: number of components in this file (inc. generated components)
    uint32_t num_vbs;                  // ZIP: z_file/txt_file PIZ: txt_file: number of VBs processed z_file: total VBs in file
    uint32_t num_preproc_vbs_joined;   // Z_FILE: PIZ
    uint32_t max_conc_writing_vbs;     // Z_FILE: PIZ: SAM: the maximal number of hand-over VBs writer might have in memory. Passed through SectionHeaderGenozipHeader since 15.0.64, and through SectionHeaderReconPlan before that.
    uint64_t deep_stats[NUM_DEEP_STATS];  // Z_FILE: ZIP/PIZ: SAM: stats collection on Deep performance
    uint64_t secondary_count, supplementary_count, saggy_near_count, mate_line_count, depn_far_count; // Z_FILE ZIP: SAM: for stats
    union {
    uint64_t num_sequences;            // Z_FILE: ZIP: FASTA: num "DESC" lines in this file. for stats
    Ploidy max_ploidy;                 // Z_FILE: ZIP: VCF 
    Ploidy max_ploidy_for_mux;         // Z_FILE: PIZ: VCF: copied from SectionHeaderGenozipHeader.max_ploidy_for_mux
    };
    uint64_t sam_num_seq_by_aln;       // Z_FILE: ZIP: SAM/BAM: number of alignments segged vs reference by rname/pos/cigar (i.e. not aligner, not copy from prim/saggy, not verbatim)
    uint64_t sam_num_perfect_matches;  // Z_FILE: ZIP: SAM/BAM/FASTQ: number of perfect matches found by aligner. for stats 
    uint64_t sam_num_aligned;          // Z_FILE: ZIP: SAM/BAM: number of alignments successfully found by aligner. for stats 
    uint64_t sam_num_verbatim;         // Z_FILE: ZIP: SAM/BAM: number of alignments segged verbatim
    uint64_t sam_num_by_prim;          // Z_FILE: ZIP: SAM/BAM: number of alignments segged against prim / saggy
    uint64_t fq_num_aligned;           // Z_FILE: ZIP: FASTQ: number of alignments successfully found by aligner. for stats 
    uint64_t fq_num_perfect_matches;   // Z_FILE: ZIP: SAM/BAM/FASTQ: number of perfect matches found by aligner. for stats 
    uint64_t fq_num_verbatim;          // Z_FILE: ZIP: FASTQ: number of lines segged verbatim
    uint64_t fq_num_monochar;          // Z_FILE: ZIP: FASTQ: number of lines segged monochar
    uint64_t vcf_num_samples_copied;   // Z_FILE: ZIP: VCF: number of samples copied
    uint64_t domq_lines[MAX_NUM_COMPS];// Z_FILE: ZIP: number of lines per components of the various QUAL codecs (show in codec_qual_show_stats)
    uint64_t divr_lines[MAX_NUM_COMPS];
    uint64_t homp_lines[MAX_NUM_COMPS];
    uint64_t pacb_lines[MAX_NUM_COMPS];
    uint64_t longr_lines[MAX_NUM_COMPS];
    uint64_t normq_lines[MAX_NUM_COMPS];

    // per-component data 
    uint64_t comp_num_lines[MAX_NUM_COMPS];             // Z_FILE: PIZ/ZIP: number of lines in each component 
    QnameFlavorProp flav_prop[MAX_NUM_COMPS][NUM_QTYPES]; // Z_FILE: PIZ: qname properties as read from SectionHeaderTxtHeader.flav_prop. Used in SAM for Deep and Gencomp, and used in genocat SAM/FASTQ for qname filterning.
    int64_t txt_file_disk_sizes_sum;                    // Z_FILE ZIP: sum of txt_file_disk_sizes[*]
    int64_t txt_file_disk_sizes[MAX_NUM_COMPS];         // Z_FILE ZIP: size of original txt file (whether or not compressed) of each component. 
    int64_t disk_so_far_comp[MAX_NUM_COMPS];            // Z_FILE ZIP: per-component size if z_file VB sections (note: global area sections, including SEC_DICT, are not accounted for)
    int64_t txt_data_so_far_bind_comp[MAX_NUM_COMPS];   // Z_FILE ZIP: per-component txt_size after modifications (due to --optimzie etc)
    int64_t txt_data_so_far_bind_0_comp[MAX_NUM_COMPS]; // Z_FILE ZIP: per-component txt_size before modifications 
    Codec comp_src_codec[MAX_NUM_COMPS];                // Z_FILE ZIP: source codec used for every txt file component (i.e. excluding generated components)
    Codec comp_eff_codec[MAX_NUM_COMPS];                // Z_FILE ZIP: effective_codec used for every txt file component with a GZIP codec
    FlagsMgzip comp_bgzf[MAX_NUM_COMPS];                // Z_FILE ZIP BGZF: library and level of BGZF of each comp
    uint64_t gz_isize[MAX_NUM_COMPS][2];                // Z_FILE ZIP GZ: isize(=uncomp_size) of the first two MGZIP blocks (excluding known MGZIP codecs). 
    uint32_t comp_num_EOF_blocks[MAX_NUM_COMPS];        // Z_FILE ZIP MGZIP: number of EOF blocks encountered in the component
    
    #define GZ_HEADER_LEN 100
    union {
    uint8_t comp_gz_header[MAX_NUM_COMPS][GZ_HEADER_LEN]; // Z_FILE ZIP gzip codecs: first (usually all) bytes of the gz header 
    uint8_t gz_header[GZ_HEADER_LEN];                   // TXT_FILE ZIP: first gz_header of file (copied later to comp_gz_header)
    };
    uint32_t gz_header_len;                             // TXT_FILE ZIP: gz_header length of first gz block
} File;

#define z_has_gencomp (z_file && z_file->z_flags.has_gencomp) // ZIP/PIZ
#define z_sam_gencomp (z_file && (Z_DT(SAM) || Z_DT(BAM)) && z_has_gencomp) // note: is BAM file in piz are Z_DT(SAM) and in zip are Z_DT(BAM)

// methods
extern FileP file_open_z_read (rom filename);
extern FileP file_open_z_write (rom filename, FileMode mode, DataType data_type, Codec src_codec);
extern StrText file_get_z_run_time (FileP file);
extern FileP file_open_txt_read (rom filename);
extern FileP file_open_txt_write (rom filename, DataType data_type, MgzipLevel bgzf_level);
extern void file_close (FileP *file_p);
extern bool file_seek (FileP file, int64_t offset, int whence, rom mode, FailType soft_fail); // SEEK_SET, SEEK_CUR or SEEK_END
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

extern bool file_put_data (rom filename, const void *data, uint64_t len, mode_t chmod_to);
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

#define SRC_CODEC(x) (txt_file->src_codec == CODEC_##x)

#define SC(x) (file->src_codec == CODEC_##x)
static inline bool is_read_via_ext_decompressor(ConstFileP file) { return SC(XZ)|| SC(ZIP) || SC(BCF)|| SC(CRAM) || SC(ORA); }
#undef SC

#define FC(x) (codec == CODEC_##x)
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
