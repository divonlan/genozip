// ------------------------------------------------------------------
//   sections.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "digest.h"

// note: the numbering of the sections cannot be modified, for backward compatibility
typedef enum __attribute__ ((__packed__)) { // 1 byte
    SEC_NONE            = -1, // doesn't appear in the file 

    SEC_RANDOM_ACCESS   = 0,
    SEC_REFERENCE       = 1,
    SEC_REF_IS_SET      = 2,
    SEC_REF_HASH        = 3,
    SEC_REF_RAND_ACC    = 4,
    SEC_REF_CONTIGS     = 5,
    SEC_GENOZIP_HEADER  = 6, // SEC_GENOZIP_HEADER has been 6 since v2, so we can always read old versions' genozip header
    SEC_DICT_ID_ALIASES = 7,
    SEC_TXT_HEADER      = 8, 
    SEC_VB_HEADER       = 9,
    SEC_DICT            = 10, 
    SEC_B250            = 11, 
    SEC_LOCAL           = 12, 
    SEC_CHROM2REF_MAP   = 13,
    SEC_STATS           = 14,
    SEC_BGZF            = 15, // optionally appears per component (txt header) and contains the uncompressed sizes of the source file bgzf block
    SEC_RECON_PLAN      = 16, // introduced v12
    SEC_COUNTS          = 17, // introduced v12
    SEC_REF_IUPACS      = 18, // introduced v12

    NUM_SEC_TYPES // fake section for counting
} SectionType;

// this data must be perfectly aligned with SectionType.
#define SECTIONTYPE_ABOUT {    \
    {"SEC_RANDOM_ACCESS",   sizeof (SectionHeader)              }, \
    {"SEC_REFERENCE",       sizeof (SectionHeaderReference)     }, \
    {"SEC_REF_IS_SET",      sizeof (SectionHeaderReference)     }, \
    {"SEC_REF_HASH",        sizeof (SectionHeaderRefHash)       }, \
    {"SEC_REF_RAND_ACC",    sizeof (SectionHeader)              }, \
    {"SEC_REF_CONTIGS",     sizeof (SectionHeader)              }, \
    {"SEC_GENOZIP_HEADER",  sizeof (SectionHeaderGenozipHeader) }, \
    {"SEC_DICT_ID_ALIASES", sizeof (SectionHeader)              }, \
    {"SEC_TXT_HEADER",      sizeof (SectionHeaderTxtHeader)     }, \
    {"SEC_VB_HEADER",       sizeof (SectionHeaderVbHeader)      }, \
    {"SEC_DICT",            sizeof (SectionHeaderDictionary)    }, \
    {"SEC_B250",            sizeof (SectionHeaderCtx)           }, \
    {"SEC_LOCAL",           sizeof (SectionHeaderCtx)           }, \
    {"SEC_CHROM2REF_MAP",   sizeof (SectionHeader)              }, \
    {"SEC_STATS",           sizeof (SectionHeader)              }, \
    {"SEC_BGZF",            sizeof (SectionHeader)              }, \
    {"SEC_RECON_PLAN",      sizeof (SectionHeaderReconPlan)     }, \
    {"SEC_COUNTS",          sizeof (SectionHeaderCounts)        }, \
    {"SEC_REF_IUPACS",      sizeof (SectionHeader)              }, \
}

// Section headers - big endian

#define GENOZIP_MAGIC 0x27052012

#pragma pack(1) // structures that are part of the genozip format are packed.

// section headers are encoded in Big Endian (see https://en.wikipedia.org/wiki/Endianness)
// the reason for selecting big endian is that I am developing on little endian CPU (Intel) so
// endianity bugs will be discovered more readily this way

// note: #pragma pack doesn't affect enums
typedef enum __attribute__ ((__packed__)) { BGZF_LIBDEFLATE, BGZF_ZLIB, NUM_BGZF_LIBRARIES } BgzfLibraryType; // constants for BGZF FlagsBgzf.library
typedef enum __attribute__ ((__packed__)) { STORE_NONE, STORE_INT, STORE_FLOAT, STORE_INDEX } StoreType; // values for SectionFlags.ctx.store

// goes into SectionHeader.flags and also SectionEnt.flags
typedef union SectionFlags {  
    uint8_t flags;

    struct FlagsGenozipHeader {
        // note: if updating dts_* flags, update in zfile_compress_genozip_header, sections_show_header too
        #define dts_ref_internal dt_specific // SAM: REF_INTERNAL was used for compressing (i.e. SAM file without reference) (introduced v6)
        #define dts_paired       dt_specific // FASTQ: This z_file contains one or more pairs of FASTQs compressed with --pair (introduced v9.0.13)
        #define dts_mismatch     dt_specific // CHAIN: This chain file's contigs mismatch its references, so it cannot be used with --chain (introduced v12.0.35)
        uint8_t dt_specific      : 1; // this flag has a different meaning depending on the data_type, may be one of the above ^ 
        uint8_t aligner          : 1; // SAM, FASTQ: our aligner may have used to align sequences to the reference (always with FASTQ, sometimes with SAM)
        uint8_t txt_is_bin       : 1; // BAM: Source file is binary (BAM)
        uint8_t bgzf             : 1; // Reconstruct as BGZF (user may override) (determined by the last component)
        uint8_t adler            : 1; // true if Adler32 is used, false if MD5 is used (>= v9) or (either MD5 or nothing) (v8)
        uint8_t has_gencomp      : 1; // VCF: file supports dual coordinates - last two components are the "liftover rejects" data (v12)
                                      // SAM/BAM: PRIM and/or DEPN components exist 
        uint8_t has_taxid        : 1; // each line in the file has Taxonomic ID information (v12)
        uint8_t unused           : 1;
    } genozip_header;

    struct FlagsTxtHeader {
        uint8_t gencomp_num      : 2; // DVCF: DC_PRIMARY/DC_LUFT contains "##primary_only"/"##luft_only" variants or DC_NONE if not a rejects component (added v12)
                                      // SAM/BAM: CT_SA_PRIM/CT_SA_DEPN - added 13.0.9
        uint8_t is_txt_luft      : 1; // true if original source file was a dual-coordinates file in Luft rendition (v12)
        uint8_t unused           : 5;
    } txt_header;

    union FlagsVbHeader {
        struct FlagsVbHeaderVcf {
            uint8_t coords             : 2; // DC_PRIMARY if it contains TOPLEVEL container, DC_LUFT if LUFT toplevel container, or DC_BOTH if both (DC_NONE prior to v12)
            uint8_t use_null_DP_method : 1; // since 13.0.5: "null_DP" mehod is used
            uint8_t unused             : 5;
        } vcf;
        struct FlagsVbHeaderSam {
            uint8_t unused             : 2; // for now, reserved for coords, should we have a dual-coord SAM in the future
            uint8_t is_sorted          : 1; // Likely sorted (copied from segconf.sam_is_sorted)     - introduced 13.0.3
            uint8_t is_collated        : 1; // Likely collated (copied from segconf.sam_is_collated) - introduced 13.0.3
            uint8_t unused2            : 4;
        } sam;
    } vb_header;

    struct FlagsBgzf {
        uint8_t has_eof_block    : 1;
        uint8_t level            : 4; // 0-12 for libdeflate or 0-9 for zlib level: 15 means unknown
        BgzfLibraryType library  : 3; // ignored if level=15 (introduced 9.0.16)
    } bgzf;

    struct FlagsCtx {
        StoreType store          : 2; // after reconstruction of a snip, store it in ctx.last_value
        uint8_t paired           : 1; // reconstruction of this context requires access to the same section from the same vb of the previous (paired) file
        #define v8_container     store_delta  // in v8 files - if the context contains 1 or more containers
        uint8_t store_delta      : 1; // introduced v12.0.41: after reconstruction of a snip, store last_delta. notes: 1. last_delta also stored in case of a delta snip. 2. if using this, store=STORE_INT must be set.
        uint8_t copy_local_param : 1; // copy ctx.b250/local.param from SectionHeaderCtx.param
        uint8_t all_the_same     : 1; // SEC_B250: the b250 data contains only one element, and should be used to reconstruct any number of snips from this context
        #define delta_peek ctx_specific_flag // v13.0.5: Valid for contexts that use SNIP_OTHER_DELTA: whether reconstruct_from_delta should peek a value or use last_value
        uint8_t ctx_specific_flag: 1; // flag specific a context (introduced 10.0.3)
        uint8_t store_per_line   : 1; // introduced v12.0.41: store last_txt for each line - in context->txt_per_prev        
    } ctx;

    struct FlagsRandomAccess {
        uint8_t luft             : 1; // this section is for reconstructing the file in the Luft coordinates (introduced v12)
        uint8_t unused           : 7;
    } random_access;
    
    struct FlagsReconPlan {
        uint8_t luft             : 1; // this section is for reconstructing the file in the Luft coordinates (introduced v12)
        uint8_t unused           : 7;
    } recon_plan;

} SectionFlags;

#define SECTION_FLAGS_NONE ((SectionFlags){ .flags = 0 })

typedef struct SectionHeader {
    uint32_t     magic; 
    #define      uncomp_adler32 magic // used in --verify-codec
    uint32_t     compressed_offset;   // number of bytes from the start of the header that is the start of compressed data (sizeof header + header encryption padding)
    uint32_t     data_encrypted_len;  // = data_compressed_len + padding if encrypted, 0 if not
    uint32_t     data_compressed_len;
    uint32_t     data_uncompressed_len;
    uint32_t     vblock_i;            // VB with in file starting from 1 ; 0 for non-VB sections
    SectionType  section_type;        // 1 byte
    Codec        codec;               // 1 byte - primary codec in which this section is compressed
    Codec        sub_codec;           // 1 byte - sub codec, in case primary codec invokes another codec
    SectionFlags flags;                
} SectionHeader; 

typedef struct {
    SectionHeader h;
    uint8_t  genozip_version;
    EncryptionType encryption_type;   // one of ENC_TYPE_*
    uint16_t data_type;               // one of DATA_TYPE_*
    uint64_t recon_size_prim;         // data size of reconstructed file, if uncompressing as a single file in primary coordinates
    uint64_t num_lines_bound;         // number of lines in a bound file. "line" is data_type-dependent. For FASTQ, it is a read.
    uint32_t num_sections;            // number sections in this file (including this one)
    uint32_t num_components;          // number of txt bound components in this file (1 if no binding)
    Digest   digest_bound;
    uint8_t  password_test[16];       // short encrypted block - used to test the validy of a password
#define FILE_METADATA_LEN 72
    char     created[FILE_METADATA_LEN];  // nul-terminated metadata
    Digest   license_hash;            // MD5(license_num)
#define REF_FILENAME_LEN 256
    char     ref_filename[REF_FILENAME_LEN]; // external reference filename, nul-terimated. ref_filename[0]=0 if there is no external reference.
    Digest   ref_file_md5;            // SectionHeaderGenozipHeader.digest_bound.md5 of the reference FASTA genozip file
    union {
        struct {
            char prim_filename[REF_FILENAME_LEN]; // external primary coordinates reference file, nul-terimated. added v12.
            Digest prim_file_md5;     // SectionHeaderGenozipHeader.digest_bound.md5 of the primary reference file. added v12.
        } chain;
    } dt_specific;    
} SectionHeaderGenozipHeader;

// this footer appears AFTER the genozip header data, facilitating reading the genozip header in reverse from the end of the file
typedef struct {
    uint64_t genozip_header_offset;
    uint32_t magic;
} SectionFooterGenozipHeader;

// The text file header section appears once in the file (or multiple times in case of bound file), and includes the txt file header 
typedef struct {
    SectionHeader h;
    uint64_t txt_data_size;           // number of bytes in the original txt file
    uint64_t txt_num_lines;           // number of data (non-header) lines in the original txt file. Concat mode: entire file for first SectionHeaderTxtHeader, and only for that txt if not first
    uint32_t max_lines_per_vb;        // upper bound on how many data lines a VB can have in this file
    Codec    codec;                   // codec of original txt file (none, bgzf, gz, bz2...)
    uint8_t  codec_info[3];           // codec specific info: for CODEC_BGZF, these are the LSB, 2nd-LSB, 3rd-LSB of the source BGZF-compressed file size
    Digest   digest_single;           // digest of original single txt file. non-0 only if this genozip file is a result of binding. MD5 if --md5 or Adler32 otherwise. 0 if compressed in v8 without --md5. 
    Digest   digest_header;           // MD5 or Adler32 of header
#define TXT_FILENAME_LEN 256
    char     txt_filename[TXT_FILENAME_LEN]; // filename of this single component. without path, 0-terminated. always in base form like .vcf or .sam, even if the original is compressed .vcf.gz or .bam
    uint64_t txt_header_size;         // size of header in original txt file (likely different than reconstructed size if dual-coordinates)  (v12)
} SectionHeaderTxtHeader; 

typedef struct {
    SectionHeader h;            
    uint32_t unused;                  // "unused" since v12 (up to v11 it was "uint32_t first_line; // if 0, this is the terminating section of the components")
    uint32_t top_level_repeats;       // repeats of TOPLEVEL container in this VB (was called num_lines before v12)
    uint32_t recon_size_prim;         // size of vblock as it appears in the default PRIMARY reconstruction
    uint32_t z_data_bytes;            // total bytes of this vblock in the genozip file including all sections and their headers 
    uint32_t longest_line_len;        // length of the longest line in this vblock 
    Digest   digest_so_far;           // partial calculation of MD5 or Adler32 up to and including this VB 
    uint32_t num_lines_prim;          // number of lines in default reconstruction in PRIMARY coords (v12)
    uint32_t num_lines_luft;          // number of lines in default reconstruction in LUFT coords (v12)
    uint32_t recon_size_luft;         // size of vblock as it appears in the default LUFT reconstruction (v12)
} SectionHeaderVbHeader; 

typedef struct {
    SectionHeader h;           
    uint32_t num_snips;               // number of items in dictionary
    DictId   dict_id;           
} SectionHeaderDictionary;    

typedef struct {
    SectionHeader h;           
    int64_t  nodes_param;             // an extra piece of data transferred to/from Context.counts_extra
    DictId   dict_id;           
} SectionHeaderCounts;     

// LT_* values are consistent with BAM optional 'B' types (and extend them)
typedef enum __attribute__ ((__packed__)) { // 1 byte
    LT_TEXT      = 0,   // 0-seperated snips
    LT_INT8      = 1,    
    LT_UINT8     = 2,
    LT_INT16     = 3,
    LT_UINT16    = 4,
    LT_INT32     = 5,
    LT_UINT32    = 6,
    LT_INT64     = 7,   // ffu
    LT_UINT64    = 8,   // ffu
    LT_FLOAT32   = 9,   
    LT_FLOAT64   = 10,  // ffu
    LT_SEQUENCE  = 11,  // length of data extracted is determined by vb->seq_len
    LT_BITMAP    = 12,  // a bitmap
    LT_CODEC     = 13,  // codec specific type with its codec specific reconstructor
    LT_UINT8_TR  = 14,  // transposed array - number of columns in original array is in param (up to 255 columns)
    LT_UINT16_TR = 15,  // "
    LT_UINT32_TR = 16,  // "
    LT_UINT64_TR = 17,  // "
    NUM_LOCAL_TYPES
} LocalType;

typedef struct LocalTypeDesc {
    const char *name;
    const char sam_type;
    unsigned width;
    bool is_signed;
    int64_t min_int, max_int; // relevant for integer fields only
    BgEnBuf file_to_native;
} LocalTypeDesc;

extern const LocalTypeDesc lt_desc[NUM_LOCAL_TYPES];
#define LOCALTYPE_DESC { \
/*   name   sam  wid signed min_int                max_int                file_to_native */ \
   { "TXT", 0,   1,  0,     0,                     0,                     0                        }, \
   { "I8 ", 'c', 1,  1,     -0x80LL,               0x7fLL,                BGEN_deinterlace_d8_buf  }, \
   { "U8 ", 'C', 1,  0,     0,                     0xffLL,                BGEN_u8_buf              }, \
   { "I16", 's', 2,  1,     -0x8000LL,             0x7fffLL,              BGEN_deinterlace_d16_buf }, \
   { "U16", 'S', 2,  0,     0,                     0xffffLL,              BGEN_u16_buf             }, \
   { "I32", 'i', 4,  1,     -0x80000000LL,         0x7fffffffLL,          BGEN_deinterlace_d32_buf }, \
   { "U32", 'I', 4,  0,     0,                     0xffffffffLL,          BGEN_u32_buf             }, \
   { "I64", 0,   8,  1,     -0x8000000000000000LL, 0x7fffffffffffffffLL,  BGEN_deinterlace_d64_buf }, \
   { "U64", 0,   8,  0,     0,                     0x7fffffffffffffffLL,  BGEN_u64_buf             }, /* note: our internal representation is int64_t so max is limited by that */ \
   { "F32", 'f', 4,  0,     0,                     0,                     BGEN_u32_buf             }, \
   { "F64", 0,   8,  0,     0,                     0,                     BGEN_u64_buf             }, \
   { "SEQ", 0,   1,  0,     0,                     0,                     0                        }, \
   { "BMP", 0,   8,  0,     0,                     0,                     0                        }, \
   { "COD", 0,   1,  0,     0,                     0,                     0                        }, \
   { "T8 ", 0,   1,  0,     0,                     0xffLL,                BGEN_transpose_u8_buf    }, \
   { "T16", 0,   2,  0,     0,                     0xffffLL,              BGEN_transpose_u16_buf   }, \
   { "T32", 0,   4,  0,     0,                     0xffffffffLL,          BGEN_transpose_u32_buf   }, \
   { "T64", 0,   8,  0,     0,                     0x7fffffffffffffffLL,  BGEN_transpose_u64_buf   }, \
}

// used for SEC_LOCAL and SEC_B250
typedef struct {
    SectionHeader h;
    LocalType ltype; // populated in both SEC_B250 and SEC_LOCAL: goes into ctx.ltype - type of data for the ctx.local buffer
    uint8_t param;   // Three options: 1. goes into ctx.b250/local.param if flags.copy_local_param. (4/4/2021: actually NOT b250 bc not implemented in zfile_compress_b250_data) 
                     //                2. given to comp_uncompress as a codec parameter
                     //                3. starting 9.0.11 for ltype=LT_BITMAP: number of unused bits in top bitarray word
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
    PosType pos;               // first pos within chrom (1-based) of this range         
    PosType gpos;              // first pos within genome (0-based) of this range
    uint32_t num_bases;        // number of bases (nucleotides) in this range
    uint32_t chrom_word_index; // index in contexts[CHROM].word_list of the chrom of this reference range    
} SectionHeaderReference;

typedef struct {
    SectionHeader h;
    uint8_t num_layers;        // total number of layers
    uint8_t layer_i;           // layer number of this section (0=base layer, with the most bits)
    uint8_t layer_bits;        // number of bits in layer
    uint8_t ffu;
    uint32_t start_in_layer;   // start index within layer
} SectionHeaderRefHash;

// SEC_RECON_PLAN, contains ar array of ReconPlanItem
typedef struct SectionHeaderReconPlan {
    SectionHeader h;
    uint32_t conc_writing_vbs; // max number of concurrent VBs in possesion of the writer thread needed to execute this plan    
    uint32_t vblock_mb;        // size of vblock in MB
} SectionHeaderReconPlan;

// special values of ReconPlanItem.num_lines
#define PLAN_END_OF_VB   0xffffffff
#define PLAN_FULL_VB     0xfffffffe
#define PLAN_INTERLEAVE  0xfffffffd
#define PLAN_TXTHEADER   0xfffffffc
#define PLAN_DOWNSAMPLE  0xfffffffb
#define MIN_PLAN_TYPE    0xfffffffb
typedef union {
    struct {
        uint32_t vb_i;               
        uint32_t start_line; // 0-based line within vb_i
        uint32_t num_lines;
    } range;

    struct {
        uint32_t vb_i;               
        uint32_t unused;    
        uint32_t plan_type;  // must be PLAN_END_OF_VB
    } end_of_vb;

    struct {
        uint32_t vb_i;               
        uint32_t unused;    
        uint32_t plan_type;  // must be PLAN_FULL_VB
    } full_vb;

    struct {
        uint32_t vb_i;               
        uint32_t vb2_i;    
        uint32_t plan_type;  // must be PLAN_INTERLEAVE
    } interleave;

    struct {
        uint32_t vb_i;               
        uint32_t rp_comp_i;    
        uint32_t plan_type;  // must be PLAN_TXTHEADER
    } txt_header;

    struct {
        uint32_t vb_i;               
        uint32_t num_lines;  // copied from SectionHeaderVbHeader.num_lines_[prim|luft]   
        uint32_t plan_type;  // must be PLAN_DOWNSAMPLE
    } downsample;

    struct {
        uint32_t vb_i;               
        uint32_t unused;    
        uint32_t plan_type;  
    } x; // generic

} ReconPlanItem;

// the data of SEC_SECTION_LIST is an array of the following type, as is the z_file->section_list_buf
typedef const struct SectionEnt {
    uint64_t offset;            // offset of this section in the file
    DictId dict_id;             // used if this section is a DICT, LOCAL or a B250 section
    uint32_t vblock_i;
    SectionType st;             // 1 byte
    SectionFlags flags;         // same flags as in section header, since v12 (before was "unused")
    uint8_t unused[2];         
} SectionEnt;

// the data of SEC_RANDOM_ACCESS is an array of the following type, as is the z_file->ra_buf and vb->ra_buf
// we maintain one RA entry per vb per every chrom in the the VB
typedef struct RAEntry {
    uint32_t vblock_i;           // the vb_i in which this range appears
    WordIndex chrom_index;       // before merge: node index into chrom context nodes, after merge - word index in CHROM dictionary
    PosType min_pos, max_pos;    // POS field value of smallest and largest POS value of this chrom in this VB (regardless of whether the VB is sorted)
} RAEntry; 

// the data of SEC_REF_IUPACS (added v12)
typedef struct Iupac {
    PosType gpos;
    char iupac;
} Iupac;

typedef union {
    SectionHeader common;
    SectionHeaderGenozipHeader genozip_header;
    SectionHeaderTxtHeader txt_header;
    SectionHeaderVbHeader vb_header;
    SectionHeaderDictionary dict;
    SectionHeaderCounts counts;
    SectionHeaderCtx ctx;
    SectionHeaderReference reference;
    SectionHeaderRefHash ref_hash;
    SectionHeaderReconPlan recon_plan;
    char padding[sizeof(SectionHeaderGenozipHeader) + 15]; // SectionHeaderGenozipHeader is the largest, 15=crypt_max_padding_len()
} SectionHeaderUnion;

#pragma pack()

// ---------
// ZIP stuff
// ---------

extern void sections_add_to_list (VBlockP vb, const SectionHeader *header);
extern void sections_remove_from_list (VBlockP vb, uint64_t offset, uint64_t len);
extern void sections_list_concat (VBlockP vb);

// ---------
// PIZ stuff
// ---------

extern Section section_next (Section sec);

extern Section sections_first_sec (SectionType st, bool soft_fail);
extern Section sections_last_sec (SectionType st, bool soft_fail);
extern bool sections_next_sec2 (Section *sl_ent, SectionType st1, SectionType st2);
#define sections_next_sec(sl_ent,st) sections_next_sec2((sl_ent),(st),SEC_NONE)
extern bool sections_prev_sec2 (Section *sl_ent, SectionType st1, SectionType st2);
#define sections_prev_sec(sl_ent,st) sections_prev_sec2((sl_ent),(st),SEC_NONE)

extern Section sections_last_sec4 (Section sl, SectionType st1, SectionType st2, SectionType st3, SectionType st4);
#define sections_component_last(any_sl_in_component) sections_last_sec4 ((any_sl_in_component), SEC_B250, SEC_LOCAL, SEC_VB_HEADER, SEC_RECON_PLAN)

extern uint32_t sections_count_sections (SectionType st);
extern Section sections_vb_first (uint32_t vb_i, bool soft_fail);
#define sections_vb_last(any_sl_in_vb) sections_last_sec4 ((any_sl_in_vb), SEC_B250, SEC_LOCAL, SEC_NONE, SEC_NONE)

extern void sections_count_component_vbs (Section sl, uint32_t *num_vbs, uint32_t *first_vb);
extern Section sections_pull_vb_up (uint32_t vb_i, Section sl);
extern Section sections_pull_component_up (Section txtfile_sl_after_me, Section txtfile_sl_move_me);

extern void BGEN_sections_list(void);
extern const char *st_name (SectionType sec_type);
#define sections_has_dict_id(st) ((st) == SEC_B250 || (st) == SEC_LOCAL || (st) == SEC_DICT || (st) == SEC_COUNTS)
extern SectionType sections_st_by_name (char *name);
extern uint32_t st_header_size (SectionType sec_type);

extern void sections_get_refhash_details (uint32_t *num_layers, uint32_t *base_layer_bits);

// z_file sizes
extern int64_t sections_get_section_size (Section sl);
extern int64_t sections_get_vb_size (Section sl);
extern int64_t sections_get_vb_skipped_sections_size (Section vb_header_sl);
extern int64_t sections_get_ref_size (void);

// display functions
extern void sections_show_header (const SectionHeader *header, VBlockP vb /* optional if output to buffer */, uint64_t offset, char rw);
extern void sections_show_gheader (const SectionHeaderGenozipHeader *header);
extern void sections_show_section_list (VBlockP vb/*if VB*/, const SectionHeaderGenozipHeader *header/*if z_file*/);
extern const char *lt_name (LocalType lt);
