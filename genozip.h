// ------------------------------------------------------------------
//   genozip.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#pragma once

#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <unistd.h> 
#include <string.h> // must be after inttypes

#include "website.h"

#pragma GCC diagnostic ignored "-Wunknown-pragmas"    // needed for our #pragma GENDICT
#pragma GCC diagnostic ignored "-Wmicrosoft-anon-tag" // a clang warning

// we defined these ourselves (normally defined in stdbool.h), as not always available on all platforms (namely issues with Docker Hub)
#ifndef __bool_true_false_are_defined
typedef _Bool bool;
#define true 1
#define false 0
#endif

typedef enum __attribute__ ((__packed__)) { no=0, yes=1, unknown=2 } thool; // three-way "bool"

// -----------------
// system parameters
// -----------------
#define GENOZIP_EXT ".genozip"

#define MAX_POS ((PosType)UINT32_MAX) // maximum allowed value for POS (constraint: fit into uint32 ctx.local). Note: in SAM the limit is 2^31-1

#define MAX_DICTS 2048

#define MAX_FIELDS 2048  // Maximum number of fields in a line (eg VCF variant, SAM line etc), including VCF/FORMAT fields, VCF/INFO fields GVF/ATTR fields, SAM/AUX fields etc. 

#define DEFAULT_MAX_THREADS 8 // used if num_cores is not discoverable and the user didn't specifiy --threads

#define MEMORY_WARNING_THREASHOLD 0x100000000  // (4 GB) warning in some cases that we predict that user choices would cause us to consume more than this

#define NO_CALLBACK NULL

typedef enum { POOL_MAIN,  // used for all VBs, except non-pool VBs and BGZF VBs
               POOL_BGZF,  // PIZ only: used for BGZF-compression dispatcher
               NUM_POOL_TYPES } VBlockPoolType;
#define POOL_NAMES { "MAIN", "BGZF" };

// ------------------------------------------------------------------------------------------------------------------------
// pointers used in header files - so we don't need to include the whole .h (and avoid cyclicity and save compilation time)
// ------------------------------------------------------------------------------------------------------------------------
typedef struct VBlock *VBlockP;
typedef const struct VBlock *ConstVBlockP;
#define NULL_VB ((VBlockP)0)
typedef struct VBlockPool *VBlockPoolP;
typedef struct File *FileP;
typedef const struct File *ConstFileP;
typedef struct Buffer *BufferP;
typedef const struct Buffer *ConstBufferP;
typedef struct Container *ContainerP;
typedef const struct Container *ConstContainerP;
typedef struct SmallContainer *SmallContainerP;
typedef const struct SmallContainer *ConstSmallContainerP;
typedef struct MediumContainer *MediumContainerP;
typedef const struct MediumContainer *ConstMediumContainerP;
typedef struct ContainerItem *ContainerItemP;
typedef const struct ContainerItem *ConstContainerItemP;
typedef struct Context *ContextP;
#define CTX_NONE ((ContextP)0)
typedef const struct Context *ConstContextP;
typedef struct CtxNode *CtxNodeP;
typedef const struct CtxNode *ConstMtfNodeP;
typedef struct SectionHeader *SectionHeaderP;
typedef const struct SectionEnt *Section;
typedef struct Range *RangeP;
typedef const struct Range *ConstRangeP;
typedef struct Bits *BitsP;
typedef const struct Bits *ConstBitsP;
typedef struct RAEntry *RAEntryP;
typedef const struct RAEntry *ConstRAEntryP;
typedef struct Mutex *MutexP;
typedef struct Serializer *SerializerP;
typedef struct RefStruct *Reference;
typedef struct Contig *ContigP;
typedef const struct Contig *ConstContigP;
typedef struct ContigPkg *ContigPkgP;
typedef const struct ContigPkg *ConstContigPkgP;
typedef const struct QnameFlavorStruct *QnameFlavor; 
typedef struct DispatcherData *Dispatcher;

#define VB ((VBlockP)(vb))

typedef enum { EXE_GENOZIP, EXE_GENOUNZIP, EXE_GENOLS, EXE_GENOCAT, NUM_EXE_TYPES } ExeType;
#define is_genocat   (exe_type == EXE_GENOCAT)
#define is_genozip   (exe_type == EXE_GENOZIP)
#define is_genounzip (exe_type == EXE_GENOUNZIP)
#define is_genols    (exe_type == EXE_GENOLS)

// IMPORTANT: DATATYPES GO INTO THE FILE FORMAT - THEY CANNOT BE CHANGED
typedef enum { DT_NONE=-1, // used in the code logic, never written to the file
               DT_REF=0, DT_VCF=1, DT_SAM=2, DT_FASTQ=3, DT_FASTA=4, DT_GFF3=5, DT_ME23=6, // these values go into SectionHeaderGenozipHeader.data_type
               DT_BAM=7, DT_BCF=8, DT_GENERIC=9, DT_PHYLIP=10, DT_CHAIN=11, DT_KRAKEN=12, 
               DT_LOCS=13, DT_BCL=14, NUM_DATATYPES } DataType; 
#define Z_DT(dt)   (z_file->data_type   == (dt))
#define TXT_DT(dt) (txt_file->data_type == (dt))
#define VB_DT(dt)  (vb->data_type       == (dt))

typedef enum { DTYPE_FIELD, DTYPE_1, DTYPE_2 } DictIdType;

#pragma pack(1) // structures that are part of the genozip format are packed.
#define DICT_ID_LEN    ((int)sizeof(uint64_t))    // VCF/GFF3 specs don't limit the field name (tag) length, we limit it to 8 chars. zero-padded. (note: if two fields have the same 8-char prefix - they will just share the same dictionary)
typedef union DictId {
    uint64_t num;             // num is just for easy comparisons - it doesn't have a numeric value and endianity should not be changed
    uint8_t id[DICT_ID_LEN];  // \0-padded IDs 
    uint16_t map_key[4];      // we use the first two bytes ([0]) as they key into vb/z_file->dict_id_mapper

    struct {
        #define ALT_KEY(d) (0x10000 | ((d).alt_key.b0_4 << 11) | ((d).alt_key.b5_9 << 6) | ((d).alt_key.b10_14 << 1) | (d).alt_key.b15)
        uint64_t unused1 : 3;
        uint64_t b0_4    : 5; // 5 LSb from 1st character
        uint64_t unused2 : 3;
        uint64_t b5_9    : 5; // 5 LSb from 2nd character
        uint64_t unused3 : 3;
        uint64_t b10_14  : 5; // 5 LSb from 3rd character
        uint64_t unused4 : 7;
        uint64_t b15     : 1; // 1 LSb from 4th character
        uint64_t unused5 : 32;
    } alt_key;
} DictId __attribute__((__transparent_union__)); // DictId function parameters can receive arguments that are of the type of members of the union
#pragma pack()

typedef uint16_t Did;    // index of a context in vb->contexts or z_file->contexts / a counter of contexts
#define DID_NONE ((Did)0xFFFF)
#define DID_EOL  ((Did)0xFFFe)

typedef uint8_t CompIType;    // comp_i 
#define COMP_MAIN ((CompIType)0)
#define COMP_NONE ((CompIType)255)
#define MAX_NUM_COMPS 3       // can be increased if needed

typedef uint32_t VBIType;     // vblock_i
typedef uint64_t CharIndex;   // index within dictionary
typedef int32_t WordIndex;    // used for word and node indices
typedef int64_t PosType;      // used for position coordinate within a genome
typedef int32_t LineIType;
typedef const char *rom;      // "read-only memory"
typedef const uint8_t *bytes; // read-only array of bytes

// a reference into txt_data
typedef struct __attribute__ ((__packed__)) { uint32_t index, len; } TxtWord; // 32b as VBs are limited to 2GB (usually used as reference into txt_data)
#define TXTWORD(snip) ((TxtWord){ .index = BNUMtxt (snip),    .len = snip##_len }) // get coordinates in txt_data
#define TXTWORDi(x,i) ((TxtWord){ .index = BNUMtxt (x##s[i]), .len = x##_lens[i] }) 

typedef struct __attribute__ ((__packed__)) { uint64_t index; uint32_t len; } BufWord; // see also ZWord

// a reference into data in z_file 
#define ZWORD_MAX_INDEX ((1ULL << 40) - 1ULL) // 1 TB
#define ZWORD_MAX_LEN   ((1ULL << 24) - 1ULL) // 16 MB
typedef struct __attribute__ ((__packed__)) { uint64_t index : 40;  uint64_t len : 24; } ZWord; 
#define ZWORDtxt(snip) ((ZWord){ .index = BNUMtxt (snip), .len = snip##_len }) // get coordinates in txt_data

typedef union { // 64 bit
    int64_t i;
    double f;
    TxtWord; // index into in txt_data (note: gcc/clang flag -fms-extensions is needed for this type of anonymous struct use)
} ValueType __attribute__((__transparent_union__));

// global parameters - set before any thread is created, and never change
extern uint32_t global_max_threads;
#define MAX_GLOBAL_MAX_THREADS 10000

extern rom global_cmd;            // set once in main()
extern ExeType exe_type;

// global files (declared in file.c)
extern FileP z_file, txt_file; 

// IMPORTANT: This is part of the genozip file format. Also update codec.h/codec_args
// If making any changes, update arrays in 1. codec.h 2. (for codecs that have a public file format, eg .zip) txtfile_set_seggable_size
typedef enum __attribute__ ((__packed__)) { // 1 byte
    CODEC_UNKNOWN=0, 
    CODEC_NONE=1, CODEC_GZ=2, CODEC_BZ2=3, CODEC_LZMA=4, CODEC_BSC=5, 
    CODEC_RANS8=6, CODEC_RANS32=7, CODEC_RANS8_pack=8, CODEC_RANS32_pack=9, 
    
    CODEC_ACGT    = 10, CODEC_XCGT = 11, // compress sequence data - slightly better compression LZMA, 20X faster (these compress NONREF and NONREF_X respectively)
    CODEC_HAPM    = 12, // compress a VCF haplotype matrix - transpose, then sort lines, then bz2. 
    CODEC_DOMQ    = 13, // compress SAM/FASTQ quality scores, if dominated by a single character
    CODEC_GTSHARK = 14, // compress VCF haplotype matrix with gtshark (discontinued in v12)
    CODEC_PBWT    = 15, // compress VCF haplotype matrix with pbwt

    CODEC_ARITH8=16, CODEC_ARITH32=17, CODEC_ARITH8_pack=18, CODEC_ARITH32_pack=19, 

    // external compressors (used by executing an external application)
    CODEC_BGZF=20, CODEC_XZ=21, CODEC_BCF=22, 
    CODEC_BAM=23,       // in v8 BAM was a codec which was compressed using samtools as external compressor. Since v14 we use the codec name for displaying "BAM" in stats total line.
    CODEC_CRAM=24, CODEC_ZIP=25,

    CODEC_LONGR=26, CODEC_NORMQ=27,

    NUM_CODECS,
} Codec; 

// note: the numbering of the sections cannot be modified, for backward compatibility
typedef enum __attribute__ ((__packed__)) { // 1 byte
    SEC_NONE            = -1, // doesn't appear in the file 

    SEC_RANDOM_ACCESS   = 0,  // Global section
    SEC_REFERENCE       = 1,  // Global section
    SEC_REF_IS_SET      = 2,  // Global section
    SEC_REF_HASH        = 3,  // Global section
    SEC_REF_RAND_ACC    = 4,  // Global section
    SEC_REF_CONTIGS     = 5,  // Global section
    SEC_GENOZIP_HEADER  = 6,  // Global section: SEC_GENOZIP_HEADER has been 6 since v2, so we can always read old versions' genozip header
    SEC_DICT_ID_ALIASES = 7,  // Global section
    SEC_TXT_HEADER      = 8,  // Per-component section
    SEC_VB_HEADER       = 9,  // Per-VB section
    SEC_DICT            = 10, // Global section
    SEC_B250            = 11, // Multiple per-VB sections
    SEC_LOCAL           = 12, // Multiple per-VB sections
    SEC_CHROM2REF_MAP   = 13, // Global section
    SEC_STATS           = 14, // Global section
    SEC_BGZF            = 15, // Per-component section (optional): contains the uncompressed sizes of the source file bgzf block
    SEC_RECON_PLAN      = 16, // Per-component section (optional): introduced v12
    SEC_COUNTS          = 17, // Global section: introduced v12
    SEC_REF_IUPACS      = 18, // Global section: introduced v12

    NUM_SEC_TYPES // fake section for counting
} SectionType;

typedef enum { DATA_EXHAUSTED, READY_TO_COMPUTE, MORE_DATA_MIGHT_COME } DispatchStatus;

typedef void BgEnBufFunc (BufferP buf, uint8_t *lt); // we use uint8_t instead of LocalType (which 1 byte) to avoid #including sections.h
typedef BgEnBufFunc (*BgEnBuf);

typedef enum { SKIP_PURPOSE_RECON,    // true if this section should be skipped when reading VB data for reconstruction
               SKIP_PURPOSE_PREPROC   // true if this section should be skipped when reading data for preprocessing (SAM: SA Loading FASTA: grep)
} SkipPurpose;
#define IS_SKIP(func) bool func (SectionType st, CompIType comp_i, DictId dict_id, SkipPurpose purpose)
typedef IS_SKIP ((*Skip));

// PIZ / ZIP inspired by "We don't sell Duff. We sell Fudd"
typedef enum { NO_COMMAND=-1, ZIP='z', PIZ='d' /* this is unzip */, LIST='l', LICENSE='L', VERSION='V', HELP=140, SHOW_HEADERS=10, TEST_AFTER_ZIP } CommandType;
extern CommandType command, primary_command;
#define IS_ZIP (command == ZIP)
#define IS_PIZ (command == PIZ)

// external vb - used when an operation is needed outside of the context of a specific VB;
extern VBlockP evb;

// threads
typedef int ThreadId;
#define THREAD_ID_NONE ((ThreadId)-1)

typedef unsigned __int128 uint128_t;
typedef          __int128 int128_t;

#define VER(n) (z_file->genozip_version >= (n))

// macros with arguments that evaluate only once 
#define MIN_(a, b) ({ __typeof__(a) _a_=(a); __typeof__(b) _b_=(b); (_a_ < _b_) ? _a_ : _b_; }) // GCC / clang "statement expressions" extesion: https://gcc.gnu.org/onlinedocs/gcc/Statement-Exprs.html#Statement-Exprs
#define MAX_(a, b) ({ __typeof__(a) _a_=(a); __typeof__(b) _b_=(b); (_a_ > _b_) ? _a_ : _b_; })
#ifndef ABS
#define ABS(a) ({ __typeof__(a) _a_=(a); (_a_ >= 0) ? _a_ : -_a_; })
#endif

// round up or down to the nearest
#define ROUNDUP2(x)    (((x) + 1)       & ~(typeof(x))0x1) 
#define ROUNDUP4(x)    (((x) + 3)       & ~(typeof(x))0x3)
#define ROUNDUP8(x)    (((x) + 7)       & ~(typeof(x))0x7)
#define ROUNDUP16(x)   (((x) + 0xf)     & ~(typeof(x))0xf)
#define ROUNDUP32(x)   (((x) + 0x1f)    & ~(typeof(x))0x1f)    
#define ROUNDDOWN32(x) ((x)             & ~(typeof(x))0x1f)

#define ROUNDUP64(x)   (((x) + 0x3f)    & ~(typeof(x))0x3f)
#define ROUNDDOWN64(x) ((x)             & ~(typeof(x))0x3f)
#define ROUNDUP512(x)  (((x) + 0x1ff)   & ~(typeof(x))0x1ff)    
#define ROUNDUP1M(x)   (((x) + 0xfffff) & ~(typeof(x))0xfffff) 

#define ARRAY_LEN(array) ((unsigned)(sizeof(array) / sizeof(array[0])))

#define IS_FLAG(flag, mask) (((flag) & (mask)) == (mask))

#define SWAP(a,b)    ({ typeof(a) tmp = a; a = b; b = tmp; })
#define SWAPbit(a,b) ({ uint8_t   tmp = a; a = b; b = tmp; })  // meant for bit fields 

// used for qsort sort function - receives two integers of any type and returns -1/0/1 as required to sort in ascending order
#define SORTER(func) int func (const void *a, const void *b)
typedef SORTER ((*Sorter));
#define ASCENDING_RAW(a,b) (((a) > (b)) ? 1 : (a) < (b) ? -1 : 0)
#define DESCENDING_RAW(a,b) (-ASCENDING_RAW((a), (b)))
#define ASCENDING(struct_type,struct_field) ASCENDING_RAW (((struct_type *)a)->struct_field, ((struct_type *)b)->struct_field)
#define DESCENDING(struct_type,struct_field) (-ASCENDING(struct_type, struct_field))

#define ASCENDING_SORTER(func_name,struct_type,struct_field) \
    SORTER (func_name) { return ASCENDING (struct_type, struct_field); }

#define DESCENDING_SORTER(func_name,struct_type,struct_field) \
    SORTER (func_name) { return DESCENDING (struct_type, struct_field); }

#define DO_ONCE static uint64_t do_once=0; if (!__atomic_test_and_set (&do_once, __ATOMIC_RELAXED))  

// Strings - declarations
#define SNIP(len) uint32_t snip_len=(len); char snip[len]

#define SNIPi1(op,n)            \
    char snip[24];              \
    snip[0] = (op);             \
    unsigned snip_len = 1 + str_int ((n), &snip[1]);

#define SNIPi2(op0,op1,n)       \
    char snip[24];              \
    snip[0] = (op0);            \
    snip[1] = (op1);            \
    unsigned snip_len = 2 + str_int ((n), &snip[2]);

#define SNIPi3(op0,op1,op2,n)   \
    char snip[24];              \
    snip[0] = (op0);            \
    snip[1] = (op1);            \
    snip[2] = (op2);            \
    unsigned snip_len = 3 + str_int ((n), &snip[3]);

#define STR(x)   rom x;        uint32_t x##_len
#define eSTR(x)  extern rom x; extern uint32_t x##_len
#define STR0(x)  rom x=NULL;   uint32_t x##_len=0
#define STRw(x)  char *x;      uint32_t x##_len    // writeable
#define STRw0(x) char *x=NULL; uint32_t x##_len=0  // writeable, initializedx
#define STRl(name,len) char name[len]; uint32_t name##_len
#define eSTRl(x) extern char x[]; extern uint32_t x##_len

#define STRlast(name,ctx) rom name = last_txtx((VBlockP)(vb), (ctx)); unsigned name##_len = (ctx)->last_txt.len
#define CTXlast(name,ctx)  ({ name = last_txtx((VBlockP)(vb), (ctx));          name##_len = (ctx)->last_txt.len; })

// Strings - function parameters
#define STRp(x)  rom x,   uint32_t x##_len    
#define STRc(x)  char *x, uint32_t x##_len  // string is fixed-length, but editable 
#define pSTRp(x) rom *x,  uint32_t *x##_len  
#define qSTRp(x) char *x, uint32_t *x##_len // function populates a string and updates its length
#define STRe(x)  char *x, uint32_t *x##_len // string and its length are editable 
#define STRps(x) uint32_t n_##x##s, rom *x##s, const uint32_t *x##_lens  

// Strings - function arguments
#define STRa(x)     x, x##_len                       
#define STRas(x)    n_##x##s, x##s, x##_lens                       
#define STRasi(x,i) (n_##x##s-(i)), &x##s[i], &x##_lens[i] // subset strating from item i
#define STRd(x)     x##_str, x##_len                   
#define STRb(buf)   (buf).data, (buf).len                  
#define STRi(x,i)   x##s[i], x##_lens[i]             
#define pSTRa(x)    &x, &x##_len                      
#define qSTRa(x)    x, &x##_len                      
#define cSTR(x)     x, sizeof x-1                 // a use with a string literal

// for printf %.*s argument list
#define STRf(x)    ((int)x##_len), x          
#define STRfi(x,i) x##_lens[i], x##s[i]
#define STRfb(buf) (int)(buf).len, (buf).data 
#define STRfw(txtword) (txtword).len, Btxt ((txtword).index) // used with TxtWord
#define STRfBw(buf,txtword) (txtword).len, Bc ((buf), (txtword).index) // used with TxtWord or ZWord

#define STRcpy(dst,src)    ({ if (src##_len) { memcpy(dst,src,src##_len) ; dst##_len = src##_len; } })
#define STRcpyi(dst,i,src) ({ if (src##_len) { memcpy(dst##s[i],src,src##_len) ; dst##_lens[i] = src##_len; } })
#define STRset(dst,src)    ({ dst=src; dst##_len=src##_len; })
#define STRLEN(string_literal) (sizeof string_literal - 1)
#define STRtxt(x) Btxt (x), x##_len
#define STRtxtw(txtword) Btxt ((txtword).index), (txtword).len // used with TxtWord
#define STRBw(buf,txtword) Bc ((buf), (txtword).index), (txtword).len // used with TxtWord
#define FUNCLINE rom func, uint32_t code_line
#define __FUNCLINE __FUNCTION__, __LINE__
#define ARRAYp(name) uint32_t n_##name##s, rom name##s[], uint32_t name##_lens[] // function parameters
#define ARRAYa(name) n_##name##s, name##s, name##_lens // function arguments

#define SAVE_VALUE(var) typeof(var) save_##var = var 
#define TEMP_VALUE(var,temp) typeof(var) save_##var = var ; var = (temp)
#define RESET_VALUE(var) SAVE_VALUE(var) ; var=(typeof(var))(uint64_t)0
#define RESTORE_VALUE(var) var = save_##var

// returns true if new_value has been set
typedef enum { NO_NEW_VALUE, HAS_NEW_VALUE } HasNewValue;
#define SPECIAL_RECONSTRUCTOR(func) HasNewValue func (VBlockP vb, ContextP ctx, rom snip, uint32_t snip_len, ValueType *new_value, bool reconstruct)
#define SPECIAL_RECONSTRUCTOR_DT(func) HasNewValue func (VBlockP vb_/*vb_ instead of vb*/, ContextP ctx, rom snip, uint32_t snip_len, ValueType *new_value, bool reconstruct)
typedef SPECIAL_RECONSTRUCTOR ((*PizSpecialReconstructor));

#define SPECIAL(dt,num,name,func) \
    extern SPECIAL_RECONSTRUCTOR(func); \
    enum { dt##_SPECIAL_##name = (num + 32) }; // define constant - +32 to make it printable ASCII that can go into a snip 

// translations of Container items - when genounzipping one format translated to another
typedef uint8_t TranslatorId;
#define TRANS_ID_NONE    ((TranslatorId)0)
#define TRANS_ID_UNKNOWN ((TranslatorId)255)

#define TRANSLATOR_FUNC(func) int32_t func(VBlockP vb, ContextP ctx, char *recon, int32_t recon_len, uint32_t item_prefix_len, bool validate_only)
#define TRANSLATOR(src_dt,dst_dt,num,name,func)\
    extern TRANSLATOR_FUNC(func); \
    enum { src_dt##2##dst_dt##_##name = num }; // define constant

// filter is called before reconstruction of a repeat or an item, and returns false if item should 
// not be processed. if not processed, contexts are not consumed. if we need the contexts consumed,
// the filter can either set *reconstruct=false and return true, or use a callback instead which is called after reconstruction,
// and erase the reconstructed txt_data.
// NOTE: for a callback to be called on items of a container, the Container.callback flag needs to be set
#define CONTAINER_FILTER_FUNC(func) bool func(VBlockP vb, DictId dict_id, ConstContainerP con, unsigned rep, int item, bool *reconstruct)

// called after reconstruction of each repeat, IF Container.callback or Container.is_top_level is set
#define CONTAINER_CALLBACK(func) void func(VBlockP vb, DictId dict_id, bool is_top_level, unsigned rep, ConstContainerP con, \
                                           char *recon, int32_t recon_len, rom prefixes, uint32_t prefixes_len)

// called after reconstruction of an item, if CI1_ITEM_CB is specified
#define CONTAINER_ITEM_CALLBACK(func) void func(VBlockP vb, ConstContainerItemP con_item, rom recon, uint32_t recon_len)

#define TXTHEADER_TRANSLATOR(func) void func (VBlockP comp_vb, BufferP txtheader_buf)

// IMPORTANT: This is part of the genozip file format. 
typedef enum __attribute__ ((__packed__)) {  ENC_NONE = 0, ENC_AES256 = 1, NUM_ENCRYPTION_TYPES } EncryptionType;
#define ENC_NAMES { "NO_ENC", "AES256" }

// note: #pragma pack doesn't affect enums
typedef enum __attribute__ ((__packed__)) { BGZF_LIBDEFLATE, BGZF_ZLIB, NUM_BGZF_LIBRARIES         } BgzfLibraryType; // constants for BGZF FlagsBgzf.library
#define BGZF_LIB_NAMES                    { "libdeflate",    "zlib"                                }

#define COMPRESSOR_CALLBACK(func)                                   \
void func (VBlockP vb,                                              \
           ContextP ctx, /* NULL if not compressing a context */    \
           LineIType vb_line_i,                                     \
           char **line_data, uint32_t *line_data_len,               \
           uint32_t maximum_size, /* might be less than the size available if we're sampling in codec_assign_best_codec() */ \
           bool *is_rev)          // QUAL and SEQ callbacks - must be returned in SAM/BAM and set to 0 elsewhere
#define CALLBACK_NO_SIZE_LIMIT 0xffffffff // for maximum_size

typedef COMPRESSOR_CALLBACK (LocalGetLineCB);

#define SAFE_ASSIGNx(addr,char_val,x) /* we are careful to evaluate addr, char_val only once, lest they contain eg ++ */ \
    char *__addr##x = (char*)(addr); \
    char __save##x  = *__addr##x; \
    *__addr##x= (char_val)
#define SAFE_RESTOREx(x) *__addr##x = __save##x

#define SAFE_ASSIGN(addr,char_val) SAFE_ASSIGNx ((addr), (char_val), _)

#define SAFE_NUL(addr) SAFE_ASSIGN((addr), 0)
#define SAFE_NULT(str) SAFE_ASSIGN((&str[str##_len]), 0)

#define SAFE_RESTORE SAFE_RESTOREx(_)

// sanity checks
extern void main_exit (bool show_stack, bool is_error);
static inline void exit_on_error(bool show_stack) { main_exit (show_stack, true); }
static inline void exit_ok(void) { main_exit (false, false); }

extern FILE *info_stream;
extern bool is_info_stream_terminal; // is info_stream going to a terminal

static inline void iputc(char c) { fputc ((c), info_stream); } // no flushing

#define iprintf(format, ...)     ( { fprintf (info_stream, (format), __VA_ARGS__); fflush (info_stream); } )
static inline void iprint0 (rom str) { fprintf (info_stream, "%s", str); fflush (info_stream); } 

// bring the cursor down to a newline, if needed
extern bool progress_newline_since_update;
static inline void progress_newline(void) {
    if (!progress_newline_since_update) { 
        fputc ('\n', stderr);
        progress_newline_since_update = true;
    }
}

// check for a user error
#define ASSINP(condition, format, ...)       ( { if (!(condition)) { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); if (flags_command_line()) fprintf (stderr, "\n\ncommand: %s\n", flags_command_line()); else fprintf (stderr, "\n"); fflush (stderr); exit_on_error(false); }} )
#define ASSINP0(condition, string)           ASSINP (condition, string "%s", "")
#define ABORTINP(format, ...)                ( { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); exit_on_error(false);} )
#define ABORTINP0(string)                    ABORTINP (string "%s", "")

// check for a bug - prints stack
typedef struct { char s[80]; } StrText;
extern StrText str_time (void);

#define SUPPORT "\nIf this is unexpected, please contact "EMAIL_SUPPORT".\n"
#define ASSERT(condition, format, ...)       ( { if (!(condition)) { progress_newline(); fprintf (stderr, "%s Error in %s:%u: ", str_time().s, __FUNCLINE); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, SUPPORT); fflush (stderr); exit_on_error(true); }} )
#define ASSERT0(condition, string)           ASSERT (condition, string "%s", "")
#define ASSERTISNULL(p)                      ASSERT0 (!p, "expecting "#p" to be NULL")
#define ASSERTNOTNULL(p)                     ASSERT0 (p, #p" is NULL")
#define ASSERTNOTZERO(n,name)                ASSERT ((n), "%s: %s=0", (name), #n)
#define ASSERTISZERO(n)                      ASSERT (!(n), "%s!=0", #n)
#define ASSERTW(condition, format, ...)      ( { if (!(condition) && !flag.quiet) { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); }} )
#define ASSERTW0(condition, string)          ASSERTW (condition, string "%s", "")
#define ASSRET(condition, ret, format, ...)  ( { if (!(condition)) { progress_newline(); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); return ret; }} )
#define ASSRET0(condition, ret, string)      ASSRET (condition, ret, string "%s", "")
#define ASSERTRUNONCE(string)                ( { static bool once = false; /* this code path should run only once */ \
                                                 ASSINP0 (!__atomic_test_and_set (&once, __ATOMIC_RELAXED), string); } )
#define RETURNW(condition, ret, format, ...) ( { if (!(condition)) { if (!flag.quiet) { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); } return ret; }} )
#define RETURNW0(condition, ret, string)     ( { if (!(condition)) { if (!flag.quiet) { progress_newline(); fprintf (stderr, "%s: %s\n", global_cmd, string); fflush (stderr); } return ret; } } )
#define ABORT(format, ...)                   ( { progress_newline(); fprintf (stderr, "%s Error in %s:%u: ", str_time().s, __FUNCLINE); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, SUPPORT); fflush (stderr); exit_on_error(true);} )
#define ABORT_R(format, ...) /*w/ return 0*/ ( { progress_newline(); fprintf (stderr, "%s Error in %s:%u: ", str_time().s, __FUNCLINE); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, SUPPORT); fflush (stderr); exit_on_error(true); return 0;} )
#define ABORT0(string)                       ABORT (string "%s", "")
#define ABORT0_R(string)                     ABORT_R (string "%s", "")
#define WARN(format, ...)                    ( { if (!flag.quiet) { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); } } )
#define WARN0(string)                        WARN (string "%s", "")

#define WARN_ONCE(format, ...)               ( { static bool warning_shown = false; \
                                                 if (!flag.quiet && !__atomic_test_and_set (&warning_shown, __ATOMIC_RELAXED)) \
                                                    { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); } \
                                               } ) 

#define WARN_ONCE0(string)                   WARN_ONCE (string "%s", "")

#define ASSERTGOTO(condition, format, ...)   ( { if (!(condition)) { progress_newline(); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); goto error; }} )

// exit codes
#define EXIT_OK                   0
#define EXIT_GENERAL_ERROR        1
#define EXIT_INVALID_GENOZIP_FILE 2
#define EXIT_DOWNSTREAM_LOST      3
#define EXIT_STREAM               4
#define EXIT_SIGHUP               5
#define EXIT_SIGSEGV              6
#define EXIT_ABNORMAL             7
#define NUM_EXIT_CODES            8