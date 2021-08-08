// ------------------------------------------------------------------
//   genozip.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef GENOZIP_INCLUDED
#define GENOZIP_INCLUDED

#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <unistd.h> 
#include <string.h> // must be after inttypes

// -----------------
// system parameters
// -----------------
#define GENOZIP_EXT ".genozip"

#define MAX_POS ((PosType)UINT32_MAX) // maximum allowed value for POS (constraint: fit into uint32 ctx.local). Note: in SAM the limit is 2^31-1

#define MAX_FIELDS 2048  // Maximum number of fields in a line (eg VCF variant, SAM line etc), including VCF/FORMAT fields, VCF/INFO fields GVF/ATTR fields, SAM/OPTIONAL fields etc. 

#define DEFAULT_MAX_THREADS 8 // used if num_cores is not discoverable and the user didn't specifiy --threads

#define MEMORY_WARNING_THREASHOLD  0x100000000  // (4 GB) warning in some cases that we predict that user choices would cause us to consume more than this

// ------------------------------------------------------------------------------------------------------------------------
// pointers used in header files - so we don't need to include the whole .h (and avoid cyclicity and save compilation time)
// ------------------------------------------------------------------------------------------------------------------------
typedef struct VBlock *VBlockP;
typedef const struct VBlock *ConstVBlockP;
typedef struct ZipDataLine *ZipDataLineP;
typedef const struct ZipDataLine *ConstZipDataLineP;
typedef struct File *FileP;
typedef const struct File *ConstFileP;
typedef struct Buffer *BufferP;
typedef const struct Buffer *ConstBufferP;
typedef struct Container *ContainerP;
typedef const struct Container *ConstContainerP;
typedef struct Context *ContextP;
typedef const struct Context *ConstContextP;
typedef struct CtxNode *CtxNodeP;
typedef const struct CtxNode *ConstMtfNodeP;
typedef struct SectionHeader *SectionHeaderP;
typedef const struct SectionEnt *Section;
typedef struct Range *RangeP;
typedef const struct Range *ConstRangeP;
typedef struct BitArray *BitArrayP;
typedef const struct BitArray *ConstBitArrayP;
typedef struct RAEntry *RAEntryP;
typedef const struct RAEntry *ConstRAEntryP;
typedef struct Mutex *MutexP;
typedef struct RefStruct *Reference;

typedef void BgEnBufFunc (BufferP buf, uint8_t *lt); // we use uint8_t instead of LocalType (which 1 byte) to avoid #including sections.h
typedef BgEnBufFunc (*BgEnBuf);

typedef enum { EXE_GENOZIP, EXE_GENOUNZIP, EXE_GENOLS, EXE_GENOCAT, NUM_EXE_TYPES } ExeType;

// IMPORTANT: DATATYPES GO INTO THE FILE FORMAT - THEY CANNOT BE CHANGED
typedef enum { DT_NONE=-1, // used in the code logic, never written to the file
               DT_REF=0, DT_VCF=1, DT_SAM=2, DT_FASTQ=3, DT_FASTA=4, DT_GFF3=5, DT_ME23=6, // these values go into SectionHeaderGenozipHeader.data_type
               DT_BAM=7, DT_BCF=8, DT_GENERIC=9, DT_PHYLIP=10, DT_CHAIN=11, DT_KRAKEN=12, 
               NUM_DATATYPES 
             } DataType; 

typedef enum { DTYPE_FIELD, DTYPE_1, DTYPE_2 } DictIdType;

#pragma pack(1) // structures that are part of the genozip format are packed.
#define DICT_ID_LEN    ((int)sizeof(uint64_t))    // VCF/GFF3 specs don't limit the field name (tag) length, we limit it to 8 chars. zero-padded. (note: if two fields have the same 8-char prefix - they will just share the same dictionary)
typedef union DictId {
    uint64_t num;             // num is just for easy comparisons - it doesn't have a numeric value and endianity should not be changed
    uint8_t id[DICT_ID_LEN];  // \0-padded IDs 
    uint16_t map_key;         // we use the first two bytes as they key into vb/z_file->dict_id_mapper
} DictId;
#pragma pack()

typedef uint16_t DidIType;    // index of a context in vb->contexts or z_file->contexts / a counter of contexts
#define DID_I_NONE ((DidIType)0xFFFF)

typedef uint64_t CharIndex;   // index within dictionary
typedef int32_t WordIndex;    // used for word and node indices
typedef int64_t PosType;      // used for position coordinate within a genome

typedef union { // 64 bit
    int64_t i;
    double f;
} LastValueType __attribute__((__transparent_union__));

// global parameters - set before any thread is created, and never change
extern uint32_t global_max_threads;
extern const char *global_cmd;            // set once in main()
extern ExeType exe_type;

// global files (declared in file.c)
extern FileP z_file, txt_file; 

// IMPORTANT: This is part of the genozip file format. Also update codec.h/codec_args
// If making any changes, update arrays in 1. codec.h 2. txtfile_estimate_txt_data_size
typedef enum __attribute__ ((__packed__)) { // 1 byte
    CODEC_UNKNOWN=0, 
    CODEC_NONE=1, CODEC_GZ=2, CODEC_BZ2=3, CODEC_LZMA=4, CODEC_BSC=5, // internal compressors
    
    CODEC_ACGT    = 10, CODEC_XCGT = 11, // compress sequence data - slightly better compression LZMA, 20X faster (these compress NONREF and NONREF_X respectively)
    CODEC_HAPM    = 12, // compress a VCF haplotype matrix - transpose, then sort lines, then bz2. 
    CODEC_DOMQ    = 13, // compress SAM/FASTQ quality scores, if dominated by a single character
    CODEC_GTSHARK = 14, // compress VCF haplotype matrix with gtshark (discontinued in v12)
    CODEC_PBWT    = 15, // compress VCF haplotype matrix with pbwt
    
    // external compressors (used by executing an external application)
    CODEC_BGZF=20, CODEC_XZ=21, CODEC_BCF=22, 
    V8_CODEC_BAM=23,    // in v8 BAM was a codec which was compressed using samtools as external compressor. since v9 it is a full data type, and no longer a codec.
    CODEC_CRAM=24, CODEC_ZIP=25,  

    NUM_CODECS
} Codec; 

// PIZ / ZIP inspired by "We don't sell Duff. We sell Fudd"
typedef enum { NO_COMMAND=-1, ZIP='z', PIZ='d' /* this is unzip */, LIST='l', LICENSE='L', VERSION='V', HELP='h', TEST_AFTER_ZIP } CommandType;
extern CommandType command, primary_command;

// external vb - used when an operation is needed outside of the context of a specific variant block;
extern VBlockP evb;

// threads
typedef int ThreadId;
#define THREAD_ID_NONE ((ThreadId)-1)

// macros
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b) )
#define MAX(a, b) (((a) > (b)) ? (a) : (b) )
#endif

#define IS_FLAG(flag, mask) (((flag) & (mask)) == (mask))

#define SWAP(a,b) do { typeof(a) tmp = a; a = b; b = tmp; } while(0)

// we defined these ourselves (normally defined in stdbool.h), as not always available on all platforms (namely issues with Docker Hub)
#ifndef __bool_true_false_are_defined
typedef _Bool bool;
#define true 1
#define false 0
#endif

#define STR(x)   const char *x; unsigned x##_len
#define STR0(x)  const char *x=NULL; unsigned x##_len=0
#define STRp(x)  const char *x, unsigned x##_len // for function definitions 
#define pSTRp(x) const char **x, unsigned *x##_len // for function definitions 
#define STRa(x) x, x##_len // for function call arguments
#define STRi(x,i) x##s[i], x##_lens[i] // for function call arguments
#define pSTRa(x) &x, &x##_len // for function call arguments
#define STRf(x) x##_len, x  // for printf %.*s argument list
#define STRfi(x,i) x##_lens[i], x##s[i]  // for printf %.*s argument list
#define cSTR(x) x, sizeof x-1 // a constant string and its length
#define STRcpy(dst,src) do { if (src##_len) { memcpy(dst,src,src##_len) ; dst##_len = src##_len; } } while(0)
#define STRcpyi(dst,i,src) do { if (src##_len) { memcpy(dst##s[i],src,src##_len) ; dst##_lens[i] = src##_len; } } while(0)

#define ARRAYp(name) unsigned n_##name##s, const char *name##s[], unsigned name##_lens[] // function parameters
#define ARRAYa(name) n_##name##s, name##s, name##_lens // function arguments

#define SAVE_VALUE(flag) typeof(flag) save_##flag = flag 
#define TEMP_VALUE(flag,temp) typeof(flag) save_##flag = flag ; flag = (temp)
#define RESET_VALUE(flag) SAVE_VALUE(flag) ; flag=(typeof(flag))(uint64_t)0
#define RESTORE_VALUE(flag) flag = save_##flag

// returns true if new_value has been set
#define SPECIAL_RECONSTRUCTOR(func) bool func (VBlockP vb, ContextP ctx, const char *snip, unsigned snip_len, LastValueType *new_value, bool reconstruct)
typedef SPECIAL_RECONSTRUCTOR ((*PizSpecialReconstructor));

#define SPECIAL(dt,num,name,func) \
    extern SPECIAL_RECONSTRUCTOR(func); \
    enum { dt##_SPECIAL_##name = (num + 32) }; // define constant - +32 to make it printable ASCII that can go into a snip 

// translations of Container items - when genounzipping one format translated to another
typedef uint8_t TranslatorId;
#define TRANS_ID_NONE    ((TranslatorId)0)
#define TRANS_ID_UNKNOWN ((TranslatorId)255)

#define TRANSLATOR_FUNC(func) int32_t func(VBlockP vb, ContextP ctx, char *recon, int32_t recon_len, bool validate_only)
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
#define CONTAINER_CALLBACK(func) void func(VBlockP vb, DictId dict_id, bool is_top_level, unsigned rep, unsigned num_reps, \
                                           char *recon, int32_t recon_len, const char *prefixes, uint32_t prefixes_len)

#define TXTHEADER_TRANSLATOR(func) void func (VBlockP comp_vb, BufferP txtheader_buf)

// IMPORTANT: This is part of the genozip file format. 
typedef enum __attribute__ ((__packed__)) { // 1 byte
    ENC_NONE   = 0,
    ENC_AES256 = 1,
    NUM_ENCRYPTION_TYPES
} EncryptionType;

#define ENC_NAMES { "NO_ENC", "AES256" }

#define COMPRESSOR_CALLBACK(func) \
void func (VBlockP vb, uint64_t vb_line_i, \
           char **line_data, uint32_t *line_data_len,\
           uint32_t maximum_size); // might be less than the size available if we're sampling in zip_assign_best_codec()
#define CALLBACK_NO_SIZE_LIMIT 0xffffffff // for maximum_size

typedef COMPRESSOR_CALLBACK (LocalGetLineCB);

#define SAFE_ASSIGNx(addr,char_val,x) /* we are careful to evaluate addr, char_val only once, lest they contain eg ++ */ \
    char *__addr##x = (char*)(addr); \
    char __save##x  = *__addr##x; \
    *__addr##x= (char_val)

#define SAFE_ASSIGN(addr,char_val) SAFE_ASSIGNx ((addr), (char_val), _)

#define SAFE_NUL(addr) SAFE_ASSIGN((addr), 0)

#define SAFE_RESTOREx(x) *__addr##x = __save##x
#define SAFE_RESTORE SAFE_RESTOREx(_)

// sanity checks
extern void main_exit (bool show_stack, bool is_error);
static inline void exit_on_error(bool show_stack) { main_exit (show_stack, true); }
static inline void exit_ok(void) { main_exit (false, false); }

extern FILE *info_stream;
extern bool is_info_stream_terminal; // is info_stream going to a terminal

static inline void iputc(char c) { fputc ((c), info_stream); } // no flushing

#define iprintf(format, ...)     do { fprintf (info_stream, (format), __VA_ARGS__); fflush (info_stream); } while(0)
static inline void iprint0 (const char *str) { fprintf (info_stream, "%s", str); fflush (info_stream); } 

// bring the cursor down to a newline, if needed
extern bool progress_newline_since_update;
static inline void progress_newline(void) {
    if (!progress_newline_since_update) { 
        fputc ('\n', stderr);
        progress_newline_since_update = true;
    }
}

// check for a user error
#define ASSINP(condition, format, ...)       do { if (!(condition)) { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); if (flags_command_line()) fprintf (stderr, "\n\ncommand: %s\n", flags_command_line()); else fprintf (stderr, "\n"); exit_on_error(false); }} while(0)
#define ASSINP0(condition, string)           do { if (!(condition)) { progress_newline(); fprintf (stderr, "%s: %s\n", global_cmd, string); if (flags_command_line()) fprintf (stderr, "\ncommand: %s\n", flags_command_line()); exit_on_error(false); }} while(0)
#define ABORTINP(format, ...)                do { progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(false);} while(0)
#define ABORTINP0(string)                    do { progress_newline(); fprintf (stderr, "%s: %s\n", global_cmd, string); exit_on_error(false);} while(0)

// check for a bug - prints stack
#define ASSERT(condition, format, ...)       do { if (!(condition)) { progress_newline(); fprintf (stderr, "Error in %s:%u: ", __FUNCTION__, __LINE__); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); }} while(0)
#define ASSERT0(condition, string)           do { if (!(condition)) { progress_newline(); fprintf (stderr, "Error in %s:%u: %s\n", __FUNCTION__, __LINE__, string); exit_on_error(true); }} while(0)
#define ASSERTNOTNULL(p)                     ASSERT0 (p, #p" is NULL")
#define ASSERTW(condition, format, ...)      do { if (!(condition) && !flag.quiet) { progress_newline(); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); }} while(0)
#define ASSERTW0(condition, string)          do { if (!(condition) && !flag.quiet) { progress_newline(); fprintf (stderr, "%s\n", string); } } while(0)
#define ASSRET(condition, ret, format, ...)  do { if (!(condition)) { progress_newline(); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); return ret; }} while(0)
#define ASSRET0(condition, ret, string)      do { if (!(condition)) { progress_newline(); fprintf (stderr, "%s\n", string); return ret; } } while(0)
#define RETURNW(condition, ret, format, ...) do { if (!(condition)) { if (!flag.quiet) { progress_newline(); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); } return ret; }} while(0)
#define RETURNW0(condition, ret, string)     do { if (!(condition)) { if (!flag.quiet) { progress_newline(); fprintf (stderr, "%s\n", string); } return ret; } } while(0)
#define ABORT(format, ...)                   do { progress_newline(); fprintf (stderr, "Error in %s:%u: ", __FUNCTION__, __LINE__); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true);} while(0)
#define ABORT_R(format, ...) /*w/ return 0*/ do { progress_newline(); fprintf (stderr, "Error in %s:%u: ", __FUNCTION__, __LINE__); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); return 0;} while(0)
#define ABORT0(string)                       do { progress_newline(); fprintf (stderr, "Error in %s:%u: ", __FUNCTION__, __LINE__); fprintf (stderr, "%s\n", string); exit_on_error(true);} while(0)
#define ABORT0_R(string)                     do { progress_newline(); fprintf (stderr, "Error in %s:%u: ", __FUNCTION__, __LINE__); fprintf (stderr, "%s\n", string); exit_on_error(true); return 0; } while(0)
#define WARN(format, ...)                    do { if (!flag.quiet) { progress_newline(); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); } } while(0)
#define WARN0(string)                        do { if (!flag.quiet) { progress_newline(); fprintf (stderr, "%s\n", string); } } while(0)

#define WARN_ONCE(format, ...)               do { static bool warning_shown = false; \
                                                  if (!flag.quiet && !warning_shown) { \
                                                      progress_newline(); fprintf (stderr, "%s: ", global_cmd); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); \
                                                      warning_shown = true; \
                                                  } \
                                             } while(0) 

#define WARN_ONCE0(string)                   do { static bool warning_shown = false; \
                                                  if (!flag.quiet && !warning_shown) { \
                                                      progress_newline(); fprintf (stderr, "%s\n", string); \
                                                      warning_shown = true; \
                                                  } \
                                             } while(0) 

#define ABORTINP0_ONCE(string)                  do { static bool once = false; \
                                                     ASSINP0 (!once, string); \
                                                     once = true; } while (0)

#define ASSERTGOTO(condition, format, ...)   do { if (!(condition)) { progress_newline(); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); goto error; }} while(0)

// exit codes
#define EXIT_OK                   0
#define EXIT_GENERAL_ERROR        1
#define EXIT_INVALID_GENOZIP_FILE 2
#define EXIT_DOWNSTREAM_LOST      3
#define EXIT_STREAM               4
#define EXIT_SIGHUP               5
#define EXIT_SIGSEGV              6
#define EXIT_ABNORMAL             7
#endif