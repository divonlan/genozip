// ------------------------------------------------------------------
//   genozip.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <unistd.h> 
#include <string.h> // must be after inttypes
#include <stdnoreturn.h>

#include "website.h"

#pragma GCC diagnostic ignored "-Wunknown-pragmas"    // needed for our #pragma GENDICT
#ifdef __clang__ 
#pragma GCC diagnostic ignored "-Wmicrosoft-anon-tag" // a clang warning
#endif

// we defined these ourselves (normally defined in stdbool.h), as not always available on all platforms (namely issues with Docker Hub)
#ifndef __bool_true_false_are_defined
typedef _Bool bool;
#define true 1
#define false 0
#endif

#define packed_enum enum __attribute__ ((packed))

typedef packed_enum { no=0, yes=1, unknown=2 } thool; // three-way "bool"

typedef packed_enum { 
    HARD_FAIL,   // display error and abort 
    SOFT_FAIL,   // silently return error
    WARNING_FAIL // display warning and return error
} FailType;

typedef packed_enum { RECON_OFF, RECON_ON, RECON_PREFIX_ONLY } ReconType;

typedef uint16_t __attribute__((aligned(1))) unaligned_uint16_t;
typedef uint32_t __attribute__((aligned(1))) unaligned_uint32_t;
typedef uint64_t __attribute__((aligned(1))) unaligned_uint64_t;

// uint40_t stuff
#define MAX_UINT40 0xffffffffffULL
typedef struct __attribute__((packed)) { uint32_t lo; uint8_t hi; } uint40_t;
// macros are good if lo is 32 bit and hi is 1 to 32 bits
#define U40to64(x40) (((uint64_t)(x40).hi << 32) | (uint64_t)(x40).lo)         
#define U64to40(x64) ((uint40_t){ .lo = (uint32_t)(x64), .hi = (x64) >> 32 })

typedef unsigned __int128 uint128_t;
typedef          __int128 int128_t;

typedef uint128_t Timestamp;

typedef int32_t TaxonomyId;

// -----------------
// system parameters
// -----------------
#define GENOZIP_EXT ".genozip"

#define MAX_POS ((PosType64)UINT32_MAX) // maximum allowed value for POS (constraint: fit into uint32 ctx.local). Note: in SAM the limit is 2^31-1

#define MAX_DICTS 2048

#define MAX_FIELDS 2048  // Maximum number of fields in a line (eg VCF variant, SAM line etc), including VCF/FORMAT fields, VCF/INFO fields GVF/ATTR fields, SAM/AUX fields etc. 

#define MAX_TAG_LEN 64   // including terminating nul (must be divisible by 8 for Tag struct)

#define DEFAULT_MAX_THREADS 8 // used if num_cores is not discoverable and the user didn't specifiy --threads

#define MEMORY_WARNING_THREASHOLD 0x100000000  // (4 GB) warning in some cases that we predict that user choices would cause us to consume more than this

#define MAX_QNAME_ITEMS 16 // mate + normal (matching Q?NAME defined in sam.h, fastq.h) 

#define NO_CALLBACK NULL

#define KB *((uint64_t)1<<10)
#define MB *((uint64_t)1<<20)
#define GB *((uint64_t)1<<30)
#define TB *((uint64_t)1<<40)
#define PB *((uint64_t)1<<50)

typedef packed_enum { 
    NO_POOL=-1, 
    POOL_MAIN,  // used for all VBs, except non-pool VBs and BGZF VBs
    POOL_BGZF,  // PIZ only: used for BGZF-compression dispatcher
    NUM_POOL_TYPES } VBlockPoolType;

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
typedef struct MiniContainer *MiniContainerP;
typedef const struct MiniContainer *ConstMiniContainerP;
typedef struct SmallContainer *SmallContainerP;
typedef const struct SmallContainer *ConstSmallContainerP;
typedef struct MediumContainer *MediumContainerP;
typedef const struct MediumContainer *ConstMediumContainerP;
typedef struct ContainerItem *ContainerItemP;
typedef const struct ContainerItem *ConstContainerItemP;
typedef struct Context *ContextP;
#define CTX_NONE ((ContextP)0)
typedef const struct Context *ConstContextP;
typedef struct SectionHeader *SectionHeaderP;
typedef const struct SectionHeader *ConstSectionHeaderP;
typedef const struct SectionEnt *Section;
typedef struct Range *RangeP;
typedef const struct Range *ConstRangeP;
typedef struct Bits *BitsP;
typedef const struct Bits *ConstBitsP;
typedef struct RAEntry *RAEntryP;
typedef const struct RAEntry *ConstRAEntryP;
typedef struct Mutex *MutexP;
typedef struct Serializer *SerializerP;
typedef struct Contig *ContigP;
typedef const struct Contig *ConstContigP;
typedef struct ContigPkg *ContigPkgP;
typedef const struct ContigPkg *ConstContigPkgP;
typedef const struct QnameFlavorStruct *QnameFlavor; 
typedef struct DispatcherData *Dispatcher;
typedef struct HuffmanChewer *HuffmanChewerP;
typedef struct HuffmanCodes *HuffmanCodesP;
typedef union BamCigarOp *BamCigarOpP;
typedef const union BamCigarOp *ConstBamCigarOpP;

typedef struct { char s[80];    } StrText;
typedef struct { char s[1 KB];  } StrTextLong;
typedef struct { char s[4 KB];  } StrTextSuperLong;
typedef struct { char s[16 KB]; } StrTextMegaLong; 
typedef struct { char s[128 KB]; } StrTextUltraLong; 

#include "version.h" // must be after StrText definition

#define VB ((VBlockP)(vb))

// IMPORTANT: DATATYPES GO INTO THE FILE FORMAT - THEY CANNOT BE CHANGED (note: if updating, also update single-user-stats-helper.c)
typedef packed_enum { 
    DT_NONE=-1,   // used in the code logic, never written to the file
    DT_REF=0, DT_VCF=1, DT_SAM=2, DT_FASTQ=3, DT_FASTA=4, DT_GFF=5, DT_ME23=6, // these values go into SectionHeaderGenozipHeader.data_type
    DT_BAM=7, DT_GNRIC=9, 
    DT_PHYLIP=10, // OBSOLETE: 9.0.11 - 15.0.41 
    DT_CHAIN=11,  // OBSOLETE: 11.0.8 - 15.0.41
    DT_KRAKEN=12, // OBSOLETE: 12.0.2 - 15.0.41
    DT_LOCS=13, DT_BED=14, 
    DT_BCF=8, DT_CRAM=15, // used only for OUT_DT in piz and set in GenozipHeader.data_type. It is NOT used in Z_DT/TXT_DT/VB_DT.
    NUM_DATATYPES 
} DataType; 

#define Z_DT(dt)   (z_file->data_type   == (DT_##dt))
#define TXT_DT(dt) (txt_file->data_type == (DT_##dt))
#define VB_DT(dt)  (vb->data_type       == (DT_##dt))
#define OUT_DT(dt) (flag.out_dt         == (DT_##dt))

typedef enum { DTYPE_FIELD, DTYPE_1, DTYPE_2 } DictIdType;

#define DICT_ID_LEN    ((int)sizeof(uint64_t))    // VCF/GFF3 specs don't limit the field name (tag) length, we limit it to 8 chars. zero-padded. (note: if two fields have the same 8-char prefix - they will just share the same dictionary)
typedef union DictId {
    uint64_t num;             // num is just for easy comparisons - it doesn't have a numeric value and endianity should not be changed
    uint8_t id[DICT_ID_LEN];  // \0-padded ID
    uint16_t map_key[4];      // we use the map_key[0] as the key into vb/z_file->d2d_map

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

typedef uint16_t Did;    // index of a context in vb->contexts or z_file->contexts / a counter of contexts
#define DID_NONE ((Did)0xFFFF)
#define DID_EOL  ((Did)0xFFFe)

typedef uint8_t CompIType;    // comp_i 
#define COMP_MAIN ((CompIType)0)
#define COMP_ALL  ((CompIType)254)
#define COMP_NONE ((CompIType)255)
#define MAX_NUM_COMPS 254 
#define MAX_NUM_TXT_FILES_IN_ZFILE (MAX_NUM_COMPS-2) // A high number to support many FASTQs in Deep. This is the maximum bc CompIType is uint8_t. 2 components reserved for generated components (in SAM).

typedef uint32_t VBIType;       // vblock_i
typedef uint64_t CharIndex;     // index within dictionary
typedef int32_t WordIndex;      // used for word and node indices
typedef int64_t PosType64;      // used for position coordinate within a genome
typedef int32_t PosType32;    
#define MAX_POS32 ((PosType64)0x7fffffff)  // max value of a POS field in VCF and SAM per their specifications
#define MAX_POS40 ((PosType64)0xffffffffff)// max gpos for deep and bamass

typedef int32_t LineIType;
#define NO_LINE (-1)

typedef const char *rom;        // "read-only memory"
typedef const uint8_t *bytes;   // read-only array of bytes

typedef uint8_t Ploidy;         // ploidy of a genotype

// a reference into txt_data
typedef struct { uint32_t index, len; } TxtWord; // 32b as VBs are limited to 2GB (usually used as reference into txt_data)
#define TXTWORD(snip) ((TxtWord){ .index = BNUMtxt (snip),    .len = snip##_len }) // get coordinates in txt_data
#define TXTWORDi(x,i) ((TxtWord){ .index = BNUMtxt (x##s[i]), .len = x##_lens[i] }) 

// a reference into a specific line
typedef struct { unaligned_uint16_t index/*in line*/; uint8_t len; } LineWord;  // short fields in short-read data
#define LONG_READ_LINE_INDEX_BITS 24
#define MAX_LONG_READ_LINE_INDEX ((1<<LONG_READ_LINE_INDEX_BITS) - 1)
typedef struct { uint32_t index : LONG_READ_LINE_INDEX_BITS; uint32_t len : 8; } LineWordL; // short fields in long-read data (line length is possibly > 64K) 
typedef struct { uint16_t index/*in line*/, len; } LineWordS;                   // long fields in short-read data
#define set_LineWord(dl, fld, str/*in vb*/, str_len) ({ if (IN_RANGE((int32_t)BNUMtxt(str) - (int32_t)(dl)->line_start, 0, 64 KB) && (str_len) < 256) \
                                                        (dl)->fld = (LineWord){ BNUMtxt(str) - (dl)->line_start, (str_len) }; })
#define set_LineWord_str(dl, fld, str/*in vb*/) set_LineWord((dl), fld, str, str##_len) 
#define LINEWORD(dl,snip)  ((LineWord) { .index = BNUMtxt (snip) - dl->line_start, .len = snip##_len }) // get coordinates in txt_data
#define LINEWORDL(dl,snip) ((LineWordL){ .index = BNUMtxt (snip) - dl->line_start, .len = snip##_len }) // get coordinates in txt_data

typedef struct __attribute__ ((packed,aligned(4))) { uint64_t index; uint32_t len; } BufWord; // see also ZWord

// a reference into data in z_file 
typedef struct { uint64_t index : 40; uint64_t len : 24; } ZWord; 
#define ZWORDtxt(snip) ((ZWord){ .index = BNUMtxt (snip), .len = snip##_len }) // get coordinates in txt_data

typedef union { // 64 bit
    int64_t i;
    double f;
    struct { float f, unused; } f32;  // used by bam_get_one_aux (note: struct due to clang)

    TxtWord;    // index into in txt_data (note: gcc/clang flag -fms-extensions is needed for this type of anonymous struct use)
    void *p; 
} ValueType __attribute__((__transparent_union__));
#define NO_VALUE ((ValueType){})

static inline bool is_p_in_range (const void *p, const void *start, uint64_t bytes)
{
    return p && start && (rom)p >= (rom)start && (rom)p < (rom)start + bytes;
}

// global parameters - set before any thread is created, and never change
extern uint32_t global_max_threads;
#define MAX_GLOBAL_MAX_THREADS 10000

extern char global_cmd[];            // set once in main()
extern unsigned file_i, n_files;

typedef enum { EXE_GENOZIP, EXE_GENOUNZIP, EXE_GENOCAT, EXE_GENOLS, NUM_EXE_TYPES } ExeType;
extern ExeType exe_type;

#define is_genozip   (exe_type == EXE_GENOZIP)
#define is_genounzip (exe_type == EXE_GENOUNZIP)
#define is_genocat   (exe_type == EXE_GENOCAT)
#define is_genols    (exe_type == EXE_GENOLS)


// global files (declared in file.c)
extern FileP z_file, txt_file; 

// IMPORTANT: This is part of the genozip file format. Also update codec.h/codec_args
// If making any changes, update arrays in 
// 1. CODEC_ARGS in codec.h 
// 2. for codecs that have a public file format, eg .zip: txtfile_set_seggable_size
// 3. for codecs used to compress sections: codec_show_time
typedef packed_enum { // 1 byte
    CODEC_UNKNOWN=0, CODEC_NONE=1, 
    
    // gzip codecs
    CODEC_GZ=2, CODEC_BGZF=20,    // General GZIP codecs 
    CODEC_IL1M=34,                // Illumina GZIP codecs 
    CODEC_MGZF=35, CODEC_MGSP=36, // MGI GZIP codecs
    CODEC_EMFL=37, CODEC_EMVL=38, // Element GZIP codecs
    CODEC_BAM=23,                 // treated as a gzip codec, with effective_codec=BGZF
    
    // other internal source codecs
    CODEC_BZ2=3, 

    // external source codecs (used by executing an external application)
    CODEC_XZ=21, CODEC_BCF=22, CODEC_CRAM=24, CODEC_ZIP=25, CODEC_ORA=32,

    // simple codecs used in genozip files
    CODEC_LZMA=4, CODEC_BSC=5, /*CODEC_BZ2=3,*/
    CODEC_RANB=6 /*8 bit*/, CODEC_RANW=7 /*32 bit*/, CODEC_RANb=8 /*8 bit packed*/, CODEC_RANw=9 /*32 bit packed*/, 
    CODEC_ARTB=16/*8 bit*/, CODEC_ARTW=17/*32 bit*/, CODEC_ARTb=18/*8 bit packed*/, CODEC_ARTw=19/*32 bit packed*/, 

    // complex codecs used in genozip files
    CODEC_ACGT    = 10, CODEC_XCGT = 11, // compress sequence data - slightly better compression LZMA, 20X faster (these compress NONREF and NONREF_X respectively)
    CODEC_HAPM    = 12, // compress a VCF haplotype matrix - transpose, then sort lines, then bz2. 
    CODEC_GTSHARK = 14, // compress VCF haplotype matrix with gtshark (discontinued in v12)
    CODEC_PBWT    = 15, // compress VCF haplotype matrix with pbwt
    CODEC_DOMQ    = 13, // compress QUAL data dominated by a single character (mostly binned Illumina)
    CODEC_LONGR   = 26, // compress Nanopore QUAL data 
    CODEC_PACB    = 30, // compress PacBio QUAL data 
    CODEC_HOMP    = 28, // compress Ultima QUAL data  
    CODEC_SMUX    = 31, // compress some types of MGI QUAL data
    CODEC_NORMQ   = 27, // compress QUAL data
    CODEC_T0      = 29, // compress the T0:Z field (Ultima) 
    CODEC_OQ      = 33, // compress the OQ:Z field (mostly generated by GATK BQSR) 
     
    NUM_CODECS    = 39,
} Codec; 

// note: the numbering of the sections cannot be modified, for backward compatibility
typedef packed_enum { // 1 byte
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
    SEC_B250            = 11, // Multiple sections per-VB
    SEC_LOCAL           = 12, // Multiple sections per-VB
    SEC_CHROM2REF_MAP   = 13, // Global section
    SEC_STATS           = 14, // Global section
    SEC_MGZIP           = 15, // Per-component section (optional): contains the uncompressed sizes of the source file mgzip block (called SEC_BGZF until 15.0.62)
    SEC_RECON_PLAN      = 16, // Per-component section (optional): used for DVCF v12-15.0.42, and for SAM gencomp v14-15.0.63
    SEC_COUNTS          = 17, // Global section: introduced v12
    SEC_REF_IUPACS      = 18, // Global section: introduced v12
    SEC_SUBDICTS        = 19, // Global section: introduced 15.0.25
    SEC_USER_MESSAGE    = 20, // Global section: introduced 15.0.30
    SEC_GENCOMP         = 21, // Section belonging to the SAM component (optional): used for SAM gencomp since v15.0.64
    SEC_HUFFMAN         = 22, // Section belonging to the SAM component (optional): huffman compression codes of QNAME (could be used for other data in the future)
    NUM_SEC_TYPES 
} SectionType;

typedef enum { DATA_EXHAUSTED, READY_TO_COMPUTE, MORE_DATA_MIGHT_COME } DispatchStatus;

typedef enum { SKIP_PURPOSE_RECON,    // true if this section should be skipped when reading VB data for reconstruction
               SKIP_PURPOSE_PREPROC   // true if this section should be skipped when reading data for preprocessing (SAM: SA Loading FASTA: grep)
} SkipPurpose;
#define IS_SKIP(func) bool func (SectionType st, CompIType comp_i, DictId dict_id, uint8_t f8/*not SectionFlags due to #include dependencies*/, SkipPurpose purpose)
typedef IS_SKIP ((*Skip));

// PIZ / ZIP inspired by "We don't sell Duff. We sell Fudd"
typedef enum { NO_COMMAND=-1, ZIP='z', PIZ='d' /* this is unzip */, LIST='l', LICENSE='L', VERSION='V', HELP=140, SHOW_HEADERS=10, TEST_AFTER_ZIP } CommandType;
extern CommandType command, primary_command;
#define IS_ZIP  (command == ZIP)
#define IS_PIZ  (command == PIZ)
#define IS_LIST (command == LIST)
#define IS_SHOW_HEADERS (command == SHOW_HEADERS)

typedef enum { VB_ID_EVB=-1, VB_ID_WRITER=-2, VB_ID_SEGCONF=-3, VB_ID_SCAN_VB=-4, VB_ID_COMPRESS_DEPN=-5, VB_ID_NONE=-999 } VBID;
#define NUM_NONPOOL_VBs 5

extern VBlockP evb; // External VB

// threads
typedef int ThreadId;
#define THREAD_ID_NONE ((ThreadId)-1)

#define VER(n) (z_file->genozip_version >= (n))
#define VER2(major,minor) ((z_file->genozip_version == (major) && (z_file->genozip_minor_ver >= (minor))) || \
                           (z_file->genozip_version > (major)))
                           
#define EXACT_VER(n) (z_file->genozip_version == (n))

#define load_relaxed(var)        __atomic_load_n    (&var, __ATOMIC_RELAXED)
#define load_acquire(var)        __atomic_load_n    (&var, __ATOMIC_ACQUIRE)
#define store_relaxed(var,value) __atomic_store_n   (&var, (value), __ATOMIC_RELAXED); 
#define store_release(var,value) __atomic_store_n   (&var, (value), __ATOMIC_RELEASE); 
#define increment_relaxed(var)   __atomic_fetch_add (&var, 1, __ATOMIC_RELAXED) // returns value before incrementing
#define decrement_relaxed(var)   __atomic_sub_fetch (&var, 1, __ATOMIC_RELAXED) // returns value after decrementing

// macros with arguments that evaluate only once 
#define MIN_(a, b) ({ __typeof__(a) _a_=(a); __typeof__(b) _b_=(b); (_a_ < _b_) ? _a_ : _b_; }) // GCC / clang "statement expressions" extesion: https://gcc.gnu.org/onlinedocs/gcc/Statement-Exprs.html#Statement-Exprs
#define MAX_(a, b) ({ __typeof__(a) _a_=(a); __typeof__(b) _b_=(b); (_a_ > _b_) ? _a_ : _b_; })
#ifndef ABS
#define ABS(a) ({ __typeof__(a) _a_=(a); (_a_ >= 0) ? _a_ : -_a_; })
#endif
#ifndef SQR 
#define SQR(x) ((x)*(x)) 
#endif

#define MINIMIZE(var,new_val) ({ if ((new_val) < (var)) var = (new_val); })
#define MAXIMIZE(var,new_val) ({ if ((new_val) > (var)) var = (new_val); })

#define IN_RANGE(x,min,after) ((x) >= (min) && (x) < (after)) // half_open [min,after)
#define IN_RANGX(x,min,max)   ((x) >= (min) && (x) <= (max))  // close [min,max]

#define MAXB64(x) ((uint64_t)((1ULL<<(x))-1))
#define MAXB(x) ((uint32_t)MAXB64(x))  // eg: MAXB(3) == 0b111 == 7

// round up or down to the nearest
#define ROUNDUP2(x)    (((x) + 1)       & ~(typeof(x))0x1) 
#define ROUNDDOWN2(x)  ((x)             & ~(typeof(x))0x1)
#define ROUNDUP4(x)    (((x) + 3)       & ~(typeof(x))0x3)
#define ROUNDDOWN4(x)  ((x)             & ~(typeof(x))0x3)
#define ROUNDUP8(x)    (((x) + 7)       & ~(typeof(x))0x7)
#define ROUNDDOWN8(x)  ((x)             & ~(typeof(x))0x7)
#define ROUNDUP16(x)   (((x) + 0xf)     & ~(typeof(x))0xf)
#define ROUNDDOWN16(x) ((x)             & ~(typeof(x))0xf)
#define ROUNDUP32(x)   (((x) + 0x1f)    & ~(typeof(x))0x1f)    
#define ROUNDDOWN32(x) ((x)             & ~(typeof(x))0x1f)
#define ROUNDUP64(x)   (((x) + 0x3f)    & ~(typeof(x))0x3f)
#define ROUNDDOWN64(x) ((x)             & ~(typeof(x))0x3f)
#define ROUNDUP512(x)  (((x) + 0x1ff)   & ~(typeof(x))0x1ff)    
#define ROUNDUP1M(x)   (((x) + 0xfffff) & ~(typeof(x))0xfffff) 

#define ARRAY_LEN(array) ((unsigned)(sizeof(array) / sizeof(array[0])))

#define IS_FLAG(flag, mask) (((flag) & (mask)) == (mask))

#define SWAP(a,b)     ({ typeof(a) tmp = a; a = b; b = tmp; })
#define SWAPbits(a,b) ({ uint64_t  tmp = a; a = b; b = tmp; })  // meant for bit fields (of any type uint8_t -> uint64_t)

// safe snprintf for multiple-step string building 
#define SNPRINTF(out/*StrText* */, format, ...) \
    ({ out##_len += snprintf (&out.s[out##_len], sizeof(out.s)-out##_len, (format), __VA_ARGS__); out##_len = MIN_(out##_len, sizeof(out.s)); })

#define SNPRINTF0(out/*StrText* */, str) \
    ({ out##_len += snprintf (&out.s[out##_len], sizeof(out.s)-out##_len, (str)); out##_len = MIN_(out##_len, sizeof(out.s)); })

// note: second caller will skip block, even if first caller is still executing it
#define DO_ONCE static bool do_once=0; if (!__atomic_test_and_set (&do_once, __ATOMIC_RELAXED))  

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

#define SNIPi2_2(op0,op1,n1,n2) \
    char snip[48];              \
    snip[0] = (op0);            \
    snip[1] = (op1);            \
    unsigned snip_len = 2 + str_int ((n1), &snip[2]); \
    snip[snip_len++] = ',';     \
    snip_len += str_int ((n2), &snip[snip_len])

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
#define STRw0(x) char *x=NULL; uint32_t x##_len=0  // writeable, initialized
#define sSTRl(name,len) static char name[len]; static uint32_t name##_len = (len)
#define STRl(name,len) char name[len]; uint32_t name##_len
#define mSTRl(name,multi,len) char name##s[multi][len]; uint32_t name##_len##s[multi]
#define STRli(name,len) uint32_t name##_len = (len) ; char name[name##_len] // avoid evaluating len twice
#define STRlic(name,len) uint32_t name##_len = len ; char name[len]         // integer constant len
#define eSTRl(x) extern char x[]; extern uint32_t x##_len
#define txtSTR(x,txtword) rom x = Btxt ((txtword).index); uint32_t x##_len = (txtword).len
#define ASSERT_LAST_TXT_VALID(ctx) ASSERT (is_last_txt_valid(ctx), "%s.last_txt is INVALID", (ctx)->tag_name)
#define STRlast(name,did_i)    ASSERT_LAST_TXT_VALID(CTX(did_i)); rom name = last_txt((VBlockP)(vb), did_i); uint32_t name##_len = CTX(did_i)->last_txt.len
#define SETlast(name,did_i) ({ ASSERT_LAST_TXT_VALID(CTX(did_i));     name = last_txt((VBlockP)(vb), did_i);          name##_len = CTX(did_i)->last_txt.len; })

#define STR_ARRAY_(x,n) rom x##s[n]; uint32_t x##_lens[n]
#define STR_ARRAY(x,n) STR_ARRAY_(x,n); uint32_t n_##x##s __attribute__((unused))
#define STRl_ARRAY(x,n,l) char x##s[n][l]; uint32_t x##_lens[n]
#define eSTRl_ARRAY(x,n,l) extern char x##s[n][l]; extern uint32_t x##_lens[n]
#define sSTRl_ARRAY(x,n,l) static char x##s[n][l]; static uint32_t x##_lens[n]

// Strings - function parameters
#define STRp(x)  rom x,   uint32_t x##_len    
#define STR8p(x) bytes x, uint32_t x##_len    
#define STRc(x)  char *x, uint32_t x##_len  // string is fixed-length, but editable 
#define STR8c(x) uint8_t *x, uint32_t x##_len  // string is fixed-length, but editable 
#define pSTRp(x) rom *x,  uint32_t *x##_len  
#define qSTRp(x) char *x, uint32_t *x##_len // function populates a string and updates its length
#define qSTR8p(x) uint8_t *x, uint32_t *x##_len 
#define STRps(x) uint32_t n_##x##s, rom *x##s, const uint32_t *x##_lens  

// Strings - function arguments
#define STRa(x)     x, x##_len                       
#define STRas(x)    n_##x##s, x##s, x##_lens                       
#define STRasi(x,i) (n_##x##s-(i)), &x##s[i], &x##_lens[i] // subset strating from item i
#define STRd(x)     x##_str, x##_len                   
#define STRb(buf)   (buf).data, (buf).len                  
#define STRi(x,i)   x##s[i], x##_lens[i]             
#define qSTRi(x,i)  x##s[i], &x##_lens[i]             
#define pSTRa(x)    &x, &x##_len                      
#define qSTRa(x)    x, &x##_len                      
#define cSTR(x)     x, sizeof x-1                 // a use with a string literal
#define STRlst(did_i) last_txt (VB, (did_i)), vb->last_txt_len (did_i)
#define STRlstf(did_i) vb->last_txt_len (did_i), last_txt (VB, (did_i))
#define STRlst_(ctx) last_txtx (VB, (ctx)), (ctx)->last_txt.len

// for printf %.*s argument list
#define STRf(x)    ((int)x##_len), x
#define STRfNUL(x) ((int)x##_len), ((x) ? (x) : "(null)")
#define STRfi(x,i) x##_lens[i], x##s[i]
#define STRfb(buf) (int)(buf).len, (buf).data 
#define STRfw(txtword) (txtword).len, Btxt ((txtword).index) // used with TxtWord
#define STRfBw(buf,txtword) (txtword).len, Bc ((buf), (txtword).index) // used with TxtWord or ZWord
#define STRfN(x) (x), ((x)!=1 ? "s" : "") // to match format eg "%u file%s" - singular or plural 

#define STRtxt(txtword) Btxt ((txtword).index), (txtword).len // used with TxtWord
#define STRline(dl,fld) Btxt ((dl)->line_start + (dl)->fld.index), (dl)->fld.len // used with LineWord

#define STRcpy(dst,src)    ({ if (src##_len) { memcpy(dst,src,src##_len) ; dst##_len = src##_len; } })
#define STRcpyi(dst,i,src) ({ if (src##_len) { memcpy(dst##s[i],src,src##_len) ; dst##_lens[i] = src##_len; } })
#define STRset(dst,src)    ({ dst=src; dst##_len=src##_len; })
#define STRinc(x,n)        ({ typeof(n) my_n = (n)/*eval once*/; x += my_n; x##_len -= my_n; })
#define STRdec(x,n)        ({ typeof(n) my_n = (n)/*eval once*/; x -= my_n; x##_len += my_n; })
#define STRLEN(string_literal) ((unsigned)(sizeof string_literal - 1))
#define _S(x) x, STRLEN(x)
#define _8(x) (bytes)x, STRLEN(x)
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
#define SPECIAL_RECONSTRUCTOR(func) HasNewValue func (VBlockP vb, ContextP ctx, rom snip, uint32_t snip_len, ValueType *new_value, ReconType reconstruct)
#define SPECIAL_RECONSTRUCTOR_DT(func) HasNewValue func (VBlockP vb_/*vb_ instead of vb*/, ContextP ctx, rom snip, uint32_t snip_len, ValueType *new_value, ReconType reconstruct)
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
    
typedef struct { uint64_t qname; uint32_t seq, qual; } DeepHash;

typedef enum { QNONE   = -6,
               QSAM2   = -5, // SAM QNAME2, while segging deep (for example, consensus reads flavor)
               QSAM    = -4, // SAM qname, while segging deep
               Q2orSAM = -3, // may appear as QNAME2 (eg in line3 if QNAME1 is NCBI) or in SAM
               Q1or3   = -2, // used in "only_q"
               QANY    = -1, 
               QNAME1  = 0,  // QNAME is SAM/BAM, line1 of FASTQ up to first space
               QNAME2  = 1,  // FASTQ: either: the remainder of line1, but excluding AUX (name=value) data, OR the second (original) QNAME on an NCBI line3, but only if different than line1  
                             // SAM/BAM: A second flavor (eg of consensus reads)
               QLINE3  = 2,  // FASTQ: The NCBI QNAME on line3, but only if different than line1 
               NUM_QTYPES/*NUM_QTYPES is part of the file format*/ } QType; 
#define QTYPE_NAME { "QNAME", "QNAME2", "LINE3" }
#define QTYPE_NEG_NAMES { "", "QANY", "Q1or3", "Q2orSAM", "QSAM", "QSAM2", "QNONE" }

// filter is called before reconstruction of a repeat or an item, and returns false if item should 
// not be processed. if not processed, contexts are not consumed. if we need the contexts consumed,
// the filter can either set *reconstruct=false and return true, or use a callback instead which is called after reconstruction,
// and erase the reconstructed txt_data.
// NOTE: for a callback to be called on items of a container, the Container.callback flag needs to be set
#define CONTAINER_FILTER_FUNC(func) bool func(VBlockP vb, DictId dict_id, ConstContainerP con, unsigned rep, int item, ReconType *reconstruct)

// called after reconstruction of each repeat, IF Container.callback or Container.is_top_level is set
#define CONTAINER_CALLBACK(func) void func(VBlockP vb_, DictId dict_id, bool is_top_level, unsigned rep, ConstContainerP con, \
                                           char *recon, int32_t recon_len, rom prefixes, uint32_t prefixes_len)

// called after reconstruction of an item, if CI1_ITEM_CB is specified
#define CONTAINER_ITEM_CALLBACK(func) void func(VBlockP vb_, ConstContainerItemP con_item, rom recon, uint32_t recon_len)

#define TXTHEADER_TRANSLATOR(func) void func (VBlockP comp_vb, BufferP txtheader_buf)

// IMPORTANT: This is part of the genozip file format. 
typedef packed_enum {  ENC_NONE = 0, ENC_AES256 = 1, NUM_ENCRYPTION_TYPES } EncryptionType;
#define ENC_NAMES { "NO_ENC", "AES256" }

// note: #pragma pack doesn't affect enums
typedef packed_enum { BGZF_LIBDEFLATE7=0, BGZF_ZLIB=1, BGZF_LIBDEFLATE19=2, BGZF_IGZIP=3, NUM_BGZF_LIBRARIES,
                      // the following are not part of the file format: used only in PIZ
                      BGZF_EXTERNAL_LIB, // level is sent to external compressor (used for BCF)
                      BGZF_NO_LIBRARY,
                      NUM_ALL_BGZF_LIBRARIES
                    } MgzipLibraryType; // constants for BGZF FlagsMgzip.library
#define BGZF_LIB_NAMES_LONG { "libdeflate_1.7", "zlib", "libdeflate_1.19", "igzip", "invalid", "external", "no_bgzf" }
#define BGZF_LIB_NAMES_SHRT { "libdef7",        "zlib", "libdef19",        "igzip", "invalid", "external", "no_bgzf" }

typedef packed_enum { BGZF_NOT_INITIALIZED    = -100,  
                      BGZF_BY_ZFILE           = -1,   // PIZ: use BGZF compression recorded in .genozip file 
                      // 0->14 are valid compression levels in a library-dependent level-space
                      BGZF_MAX_LEVEL          = 14,
                      BGZF_COMP_LEVEL_UNKNOWN = 15 
                      #define BGZF_NO_BGZF      15    // meaning if bgzf_library=BGZF_NO_LIBRARY
} MgzipLevel;

#define COMPRESSOR_CALLBACK(func)                                   \
void func (VBlockP vb,                                              \
           ContextP ctx, /* NULL if not compressing a context */    \
           LineIType vb_line_i,                                     \
           char **line_data, uint32_t *line_data_len,               \
           uint32_t maximum_size, /* might be less than the size available if we're sampling in codec_assign_best_codec() */ \
           bool *is_rev)          // QUAL and SEQ callbacks - must be returned in SAM/BAM and set to 0 elsewhere

#define COMPRESSOR_CALLBACK_DT(func)                                \
void func (VBlockP vb_,                                             \
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
    *__addr##x = (char_val)
#define SAFE_RESTOREx(x) *__addr##x = __save##x

#define SAFE_ASSIGN(addr,char_val) SAFE_ASSIGNx ((addr), (char_val), _)

#define SAFE_NUL(addr) SAFE_ASSIGN((addr), 0)
#define SAFE_NULT(str) SAFE_ASSIGN((&str[str##_len]), 0)
#define SAFE_NULB(buf) SAFE_ASSIGN((BAFTc((buf))), 0)

#define SAFE_RESTORE SAFE_RESTOREx(_)

#define RESAFE_NUL(addr) ({ SAFE_RESTORE ; __addr_ = (char*)(addr); __save_ = *__addr_ ; *__addr_ = 0; })

// exit codes of the genozip process
typedef enum {             EXIT_OK=0, EXIT_GENERAL_ERROR=1, EXIT_INVALID_GENOZIP_FILE=2, EXIT_DOWNSTREAM_LOST=3, EXIT_SIGHUP=5, EXIT_SIGSEGV=6, EXIT_ABNORMAL=7, EXIT_STREAM=127/*value that doesn't overlap curl or wget exit codes*/, NUM_EXIT_CODES } ExitCode;
#define EXIT_CODE_NAMES { "EXIT_OK", "EXIT_GENERAL_ERROR", "EXIT_INVALID_GENOZIP_FILE", "EXIT_DOWNSTREAM_LOST", "EXIT_SIGHUP", "EXIT_SIGSEGV", "EXIT_ABNORMAL", [EXIT_STREAM]="EXIT_STREAM" }

extern void noreturn main_exit (bool show_stack, bool is_error);
#define exit_on_error(show_stack) main_exit (show_stack, true)
#define exit_ok                   main_exit (false, false)

extern FILE *info_stream;
extern bool is_info_stream_terminal; // is info_stream going to a terminal

static inline void iputc(char c) { fputc ((c), info_stream); } // no flushing

#define iprintf(format, ...)     ( { fprintf (info_stream, (format), __VA_ARGS__); fflush (info_stream); } )
static inline void iprint0 (rom str) { fprintf (info_stream, "%s", str); fflush (info_stream); } 

extern noreturn void stall (void);

#if !defined(__GNUC__) // && !__has_builtin(__builtin_expect)
#define __builtin_expect(exp,c) (exp)
#endif

extern void progress_newline(void);

// check for a user error
#define ASSINP(condition, format, ...)       ( { if (__builtin_expect(!(condition), 0)) { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); if (flags_command_line()) fprintf (stderr, "\n\ncommand: %s\n", flags_command_line()); else fprintf (stderr, "\n"); fflush (stderr); exit_on_error(false); }} )
#define ASSINP0(condition, string)           ASSINP (condition, string "%s", "")
#define ABORTINP(format, ...)                ( { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); exit_on_error(false);} )
#define ABORTINP0(string)                    ABORTINP (string "%s", "")

extern StrTextLong str_time (void);

extern StrText license_get_number (void);

extern rom report_support (void);
extern rom report_support_if_unexpected (void);

#define ASSERT(condition, format, ...)       ( { if (__builtin_expect (!(condition), 0)) { progress_newline(); fprintf (stderr, "%s Error in %s:%u %s%s: ", str_time().s, __FUNCLINE, version_str().s, license_get_number().s); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "%s", report_support_if_unexpected()); fflush (stderr); exit_on_error(true); }} )
#define ASSERT0(condition, string)           ASSERT (condition, string "%s", "")
#define ASSERTISNULL(p)                      ASSERT0 (!p, "expecting "#p" to be NULL")
#define ASSERTNOTNULL(p)                     ASSERT0 (p, #p" is NULL")
#define ASSERTNOTZEROn(n,name)               ASSERT ((n), "%s: %s=0", (name), #n)
#define ASSERTNOTZERO(n)                     ASSERT ((n), "%s=0", #n)
#define ASSERTISZERO(n)                      ASSERT (!(n), "%s!=0", #n)
#define ASSERTINRANGE(n, min, after)         ASSERT (IN_RANGE((n), (min), (after)), "%s=%"PRId64" ∉ [%"PRId64",%"PRId64")", #n, (int64_t)(n), (int64_t)(min), ((int64_t)(after)))
#define ASSERTINRANGX(n, min, max)           ASSERT (IN_RANGX((n), (min), (max)),   "%s=%"PRId64" ∉ [%"PRId64",%"PRId64"]", #n, (int64_t)(n), (int64_t)(min), ((int64_t)(max)))

#define ASSERTW(condition, format, ...)      ( { if (__builtin_expect (!(condition), 0) && !flag.quiet) { progress_newline(); fprintf (stderr, "\n%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n\n"); fflush (stderr); }} )
#define WARN_IF(condition, format, ...)      ( { if (__builtin_expect ((condition), 0) && !flag.explicit_quiet) { progress_newline(); fprintf (stderr, "%s: WARNING: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n\n"); fflush (stderr); }} )
#define ASSERTW0(condition, string)          ASSERTW ((condition), string "%s", "")
#define WARN_IF0(condition, string)          WARN_IF ((condition), string "%s", "")
#define ASSERTWD(condition, format, ...)     ( { if (__builtin_expect (!(condition), 0) && flag.debug && !flag.quiet) { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); }} )
#define ASSERTWD0(condition, string)         ASSERTWD (condition, string "%s", "")
#define ASSRET(condition, ret, format, ...)  ( { if (__builtin_expect (!(condition), 0)) { progress_newline(); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); return ret; }} )
#define ASSRET0(condition, ret, string)      ASSRET (condition, ret, string "%s", "")
#define ASSERTRUNONCE(string)                ( { static bool once = false; /* this code path should run only once */ \
                                                 ASSINP0 (!__atomic_test_and_set (&once, __ATOMIC_RELAXED), string); } )
#define RETURNW(condition, ret, format, ...) ( { if (__builtin_expect (!(condition), 0)) { if (!flag.quiet) { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); } return ret; }} )
#define RETURNW0(condition, ret, string)     ( { if (__builtin_expect (!(condition), 0)) { if (!flag.quiet) { progress_newline(); fprintf (stderr, "%s: %s\n", global_cmd, string); fflush (stderr); } return ret; } } )
#define ABORT(format, ...)                   ( { progress_newline(); fprintf (stderr, "%s Error in %s:%u: ", str_time().s, __FUNCLINE); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "%s", report_support_if_unexpected()); fflush (stderr); exit_on_error(true);} )
#define ABORT0(string)                       ABORT (string "%s", "")
#define WARN(format, ...)                    ( { if (!flag.quiet) { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); } } )
#define WARN0(string)                        WARN (string "%s", "")
#define NOISYWARN(format, ...)               ( { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); fflush (stderr); } )
#define NOISYWARN0(string)                   NOISYWARN (string "%s", "")

#define WARN_ONCE(format, ...)               ( { static bool warning_shown = false; \
                                                 if (!flag.quiet && !__atomic_test_and_set (&warning_shown, __ATOMIC_RELAXED)) \
                                                    { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); } \
                                               } ) 

#define WARN_ONCE0(string)                   WARN_ONCE (string "%s", "")

#define ASSGOTO(condition, format, ...)      ( { if (__builtin_expect (!(condition), 0)) { progress_newline(); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); goto error; }} )
#define ASSGOTO0(condition, string)          ASSGOTO ((condition), string "%s", "")

extern rom main_input_filename (int file_i);
