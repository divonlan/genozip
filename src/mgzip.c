// ------------------------------------------------------------------
//   mgzip.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <errno.h>
#include <math.h>

#include "igzip/igzip_lib.h"
#include "libdeflate_1.19/libdeflate.h"
#include "libdeflate_1.7/libdeflate.h"
#include "zlib/zlib.h"
#include "xxhash/xxhash.h"
#include "mgzip.h"
#include "arch.h"
#include "file.h"
#include "txtfile.h"
#include "threads.h"
#include "dispatcher.h"
#include "writer.h"
#include "gencomp.h"
#include "strings.h"
#include "codec.h"

typedef union {
    struct {
        uint8_t text    : 1;
        uint8_t hcrc    : 1;
        uint8_t extra   : 1;
        uint8_t name    : 1;
        uint8_t comment : 1;
        uint8_t unused  : 3;
    };
    uint8_t value;
} GzipFlags;

typedef struct __attribute__ ((packed, aligned(1))) GzipHeader { // 10 bytes
    uint8_t id1;    // Gzip id            - must be 31  (0x1f)
    uint8_t id2;    // Gzip id            - must be 139 (0x8b)
    uint8_t cm;     // Compression Method - must be 8
    GzipFlags flg;  // Flags              
    uint32_t mtime; // Modification Time  
    uint8_t xfl;    // eXtra Flags. Standard values: 0, 2=best compression 4=fastest compression
    uint8_t os;     // Operating System   - (e.g. Unix=3)
} GzipHeader;

#define GZ_HEADER_LEN 10

// all data in Little Endian. Defined in https://datatracker.ietf.org/doc/html/rfc1952 and https://samtools.github.io/hts-specs/SAMv1.pdf
typedef struct __attribute__ ((packed, aligned(1))) BgzfHeader { // 18 bytes
    GzipHeader;     // note: flg=4 (extra)
    uint16_t xlen;  // Size of extra fields - 6 if contain only BGZF (may be more)
    uint8_t si1;    // bsize field id - must be 66  (0x42)
    uint8_t si2;    // bsize field id - must be 67  (0x43)
    uint16_t slen;  // bsize field length - must be 2
    uint16_t bsize; // bsize field field - (compressed block size -1)
} BgzfHeader;

// see: https://docs.google.com/document/d/11yeGa1HzXi96D3VTeMwReW2BzG9hmspyKFxx-yRpLPo/edit
typedef struct __attribute__ ((packed, aligned(1))) MgzfHeader { // 29 bytes
    GzipHeader;     // note: flg=20 (0x14) (extra + comment); mtime=0; xfl=0; os=0xff
    uint16_t xlen;  // Size of extra fields - must be 8 
    uint8_t si1;    // bsize field id - must be (0x49)
    uint8_t si2;    // bsize field id - must be (0x47)
    uint16_t slen;  // bsize field extra field length - must be 4
    uint32_t bsize; // bsize field - compressed block size - header + body
    char comment[9];// nul-terminated string in the format "C001R015"
} MgzfHeader;

typedef struct __attribute__ ((packed, aligned(1))) GzipFooter {
    uint32_t crc32; // CRC32 of uncompressed data
    uint32_t isize; // Input (i.e. uncompressed) Size
} GzipFooter;

#define GZIP_FOOTER_LEN ((int)sizeof(GzipFooter))

static const uint8_t mgzip_header_len[NUM_CODECS] = { 
    [CODEC_BGZF] = BGZF_HEADER_LEN,
    [CODEC_MGZF] = MGZF_HEADER_LEN,
    [CODEC_EMFL] = GZ_HEADER_LEN,
    [CODEC_IL1M] = GZ_HEADER_LEN,
    [CODEC_IL4M] = GZ_HEADER_LEN,
    [CODEC_MGSP] = GZ_HEADER_LEN,
    [CODEC_EMVL] = GZ_HEADER_LEN,
    [CODEC_GZBL] = GZ_HEADER_LEN,
};

static const int igzip_level_buf_lens[ISAL_DEF_MAX_LEVEL+1] = 
    { 1/*avoid 0*/+ISAL_DEF_LVL0_DEFAULT, ISAL_DEF_LVL1_DEFAULT, ISAL_DEF_LVL2_DEFAULT, ISAL_DEF_LVL3_DEFAULT };

static inline uint8_t num_bgzf_blocks_tested_for_level (void)
{
    return load_relaxed (txt_file->num_bgzf_blocks_tested_for_level);
}

#define LIB(x) (mgzip_flags.library == BGZF_##x)

//-----------------------------------------------------
// Shared ZIP/PIZ functions
//-----------------------------------------------------

rom gzstatus_name (GzStatus st)
{
    return IN_RANGE(st, 0, NUM_GZ_STATUSES) ? (rom[])GZSTATUS_NAMES[st] : "InvalidGzStatus";
}

// possible return values, see libdeflate_result in libdeflate.h
static rom libdeflate_error (int err)
{
    switch (err) {
        case LIBDEFLATE_SUCCESS            : return "Success";
        case LIBDEFLATE_BAD_DATA           : return "BadData";
        case LIBDEFLATE_SHORT_OUTPUT       : return "ShortOutput";
        case LIBDEFLATE_INSUFFICIENT_SPACE : return "InsufficientSpace";
        case LIBDEFLATE_INSUFFICIENT_DATA  : return "InsufficientData";
        case LIBDEFLATE_NOT_GZIP           : return "NotGzip";
        case LIBDEFLATE_BAD_CRC            : return "BadCRC";
	    case LIBDEFLATE_WRONG_ISIZE        : return "WrongIsize";

        default                            : return "Undefined libdeflate error";
    }
}

rom isal_error (int ret)
{
    switch (ret) {
        case ISAL_DECOMP_OK          : return "Ok";
        case ISAL_END_INPUT          : return "EndInput";
        case ISAL_OUT_OVERFLOW       : return "OutOverflow";
        case ISAL_NAME_OVERFLOW      : return "NameOverflow";
        case ISAL_COMMENT_OVERFLOW   : return "CommentOverflow";
        case ISAL_EXTRA_OVERFLOW     : return "ExtraOverflow";
        case ISAL_NEED_DICT          : return "NeedDict";
        case ISAL_INVALID_BLOCK  	 : return "InvalidBlock";
        case ISAL_INVALID_SYMBOL     : return "InvalidSymbol";
        case ISAL_INVALID_LOOKBACK   : return "InvalidLookback";
        case ISAL_INVALID_WRAPPER    : return "InvalidWrapper";
        case ISAL_UNSUPPORTED_METHOD : return "UnsupportedMethod";
        case ISAL_INCORRECT_CHECKSUM : return "IncorrectChecksum";
        default                      : return "InvalidReturnCode";
    }
}

rom zlib_error (int ret)
{
    switch (ret) {
        case Z_OK            : return "Ok";
        case Z_STREAM_END    : return "StreamEnd";
        case Z_NEED_DICT     : return "NeedDict";
        case Z_ERRNO         : return "Errno";
        case Z_STREAM_ERROR  : return "StreamError";
        case Z_DATA_ERROR    : return "DataError";
        case Z_MEM_ERROR     : return "MemError";
        case Z_BUF_ERROR     : return "BufError";
        case Z_VERSION_ERROR : return "VersionError";
        default              : return "InvalidReturnCode";
    }
}

typedef struct { char s[100]; } BgzfBlockStr;
static BgzfBlockStr display_bb (GzBlockZip *bb)
{
    BgzfBlockStr s;
    snprintf (s.s, sizeof (s.s), "{txt_index=%u txt_size=%u gz_index=%u gz_size=%u is_uncompressed=%u}",
             bb->txt_index, bb->txt_size, bb->gz_index, bb->gz_size, bb->is_uncompressed);
    return s;
}

static void *gz_alloc (void *vb_, unsigned items, unsigned size, FUNCLINE)
{
    return codec_alloc_do ((VBlockP )vb_, (uint64_t)items * (uint64_t)size, 1, NULL, func, code_line); // all bzlib buffers are constant in size between subsequent compressions
}

static void gz_alloc_deflator (VBlockP vb, FlagsMgzip mgzip_flags)
{
    ASSERTISNULL (vb->gz_deflate_mem);

    vb->gz_deflate_mem = LIB(LIBDEFLATE19) ? libdeflate_alloc_compressor (vb, mgzip_flags.level, __FUNCLINE)
                       : LIB(LIBDEFLATE7)  ? libdeflate_alloc_compressor_1_7 (mgzip_flags.level, vb)
                       : LIB(ZLIB)         ? gz_alloc (vb, 1, sizeof (z_stream), __FUNCLINE)
                       : LIB(IGZIP)        ? gz_alloc (vb, 1, igzip_level_buf_lens[mgzip_flags.level], __FUNCLINE)
                       :                     ({ ABORT ("unrecognized library=%u", mgzip_flags.library); NULL; }); // invalid library

    if (LIB(ZLIB))
        *(z_stream *)vb->gz_deflate_mem = (z_stream){ .zalloc=gz_alloc, .zfree=codec_free_do, .opaque=vb };
}

static void gz_free_deflator (VBlockP vb, FlagsMgzip mgzip_flags)
{
    ASSERTNOTNULL (vb->gz_deflate_mem);

    if (LIB(LIBDEFLATE19))   libdeflate_free_compressor (vb->gz_deflate_mem, __FUNCLINE);
    else if LIB(LIBDEFLATE7) libdeflate_free_compressor_1_7 (vb->gz_deflate_mem);
    else                     codec_free (vb, vb->gz_deflate_mem);

    vb->gz_deflate_mem = NULL;
}

static uint32_t bgzf_igzip_deflate (VBlockP vb, FlagsMgzip mgzip_flags, STRp(uncomp), STRp(comp))
{
    ASSERTNOTNULL (vb->gz_deflate_mem);
    
    struct isal_zstream strm = {};
    isal_deflate_stateless_init (&strm);

    strm.gzip_flag      = ISAL_DEFLATE; 
    strm.level          = mgzip_flags.level;
    strm.level_buf_size = igzip_level_buf_lens[mgzip_flags.level];
    strm.level_buf      = vb->gz_deflate_mem; // caller responsible for allocating and freeing
    strm.next_in        = (uint8_t *)uncomp;
    strm.avail_in       = uncomp_len;
    strm.next_out       = (uint8_t *)comp;
    strm.avail_out      = comp_len;

    int ret = isal_deflate_stateless (&strm);
    ASSERT (ret == COMP_OK || ret == STATELESS_OVERFLOW, 
            "%s: isal_deflate_stateless: %s. isize=%u", VB_NAME, isal_error (ret), uncomp_len);

    return (ret == COMP_OK) ? (comp_len - strm.avail_out) : 0/*output buffer too small*/;
}

static uint32_t bgzf_zlib_deflate (VBlockP vb, FlagsMgzip mgzip_flags, STRp(uncomp), STRp(comp))
{
    z_stream strm = { .zalloc=gz_alloc, .zfree=codec_free_do, .opaque=vb };
    // deflateInit2 with the default zlib parameters, which is also the same as htslib does
    int ret = deflateInit2 (&strm, mgzip_flags.level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY);
    ASSERT (ret == Z_OK, "deflateInit2 failed: %s", zlib_error(ret));

    strm.next_in   = (uint8_t *)uncomp;
    strm.avail_in  = uncomp_len;
    strm.next_out  = (uint8_t *)comp;
    strm.avail_out = comp_len;
    
    ret = deflate (&strm, Z_FINISH);
    ASSERT (ret != Z_STREAM_ERROR, "zlib deflate failed: ret=%s msg=%s", zlib_error(ret), strm.msg);

    comp_len = (ret == Z_STREAM_END) ? (comp_len - strm.avail_out) : 0/*output buffer too small*/;
    
    ret = deflateEnd (&strm);
    ASSERT (ret == Z_OK, "deflateEnd failed: %s", zlib_error(ret));
    return comp_len;
}

// returns comp_len, or 0 if output buffer is too small
// caller is responsible for calling gz_alloc_deflator / gz_free_deflator
static uint32_t gz_deflate (VBlockP vb, FlagsMgzip mgzip_flags, STRp(uncomp), STRc(comp), bool is_bgzf)
{
    // for an empty BGZF block, always use BGZF EOF 03.00, not library native (igzip produces 01.00.00.FF.FF)
    // Notes: assuming that any 0-block in a BGZF file is BGZF_EOF (possibly mid-file if concatenated files)
    // TO DO: if in the future we gz-recompress to non-BGZF codecs, we need to use the empty-block of that codec
    if (!uncomp_len && is_bgzf) {
        ASSERT (comp_len >= STRLEN(BGZF_EOF_CDATA), "comp_len=%u too small", comp_len);
        memcpy (comp, BGZF_EOF_CDATA, STRLEN(BGZF_EOF_CDATA));
        return STRLEN(BGZF_EOF_CDATA);
    }

    else {
        bool allocate_here = !vb->gz_deflate_mem;
        if (allocate_here) gz_alloc_deflator (vb, mgzip_flags);

        comp_len = LIB(LIBDEFLATE19) ? (uint32_t)libdeflate_deflate_compress     (vb->gz_deflate_mem, STRa(uncomp), STRa(comp))
                 : LIB(LIBDEFLATE7)  ? (uint32_t)libdeflate_deflate_compress_1_7 (vb->gz_deflate_mem, STRa(uncomp), STRa(comp))
                 : LIB(IGZIP)        ? bgzf_igzip_deflate (vb, mgzip_flags, STRa(uncomp), STRa(comp))
                 : LIB(ZLIB)         ? bgzf_zlib_deflate  (vb, mgzip_flags, STRa(uncomp), STRa(comp))
                 :                     ({ ABORT ("unrecognized library=%u", mgzip_flags.library); 0; }); // invalid library

        if (allocate_here) gz_free_deflator (vb, mgzip_flags);
    }

    return comp_len;
}

//--------------------------------------------------------------------
// ZIP SIDE - library⁀level discovery
//--------------------------------------------------------------------

#define BGZF_DISCOVERY_MAX_TESTS 10       // maximum number of BGZF blocks to be tested

void bgzf_initialize_discovery (FileP file)
{
    // note: tested example files of MGZF, MGSP and IL1M and they don't match any of these libraries.
    if (IS_BGZF(file->effective_codec)) {
        for (int l=0; l <= LIBDEFLATE_MAX_LEVEL; l++) // level=0 only here, bc it would be the same in all libraries
            file->bgzf_plausible_levels[file->num_plausible_levels++] = (FlagsMgzip){ .library = BGZF_LIBDEFLATE19, .level = l};

        for (int l=1; l <= LIBDEFLATE_MAX_LEVEL; l++)
            file->bgzf_plausible_levels[file->num_plausible_levels++] = (FlagsMgzip){ .library = BGZF_LIBDEFLATE7,  .level = l};

        for (int l=1; l <= ZLIB_MAX_LEVEL; l++)
            file->bgzf_plausible_levels[file->num_plausible_levels++] = (FlagsMgzip){ .library = BGZF_ZLIB,         .level = l};

        for (int l=0; l <= IGZIP_MAX_LEVEL; l++)
            file->bgzf_plausible_levels[file->num_plausible_levels++] = (FlagsMgzip){ .library = BGZF_IGZIP,        .level = l};

        ASSERT (file->num_plausible_levels == NUM_PLAUSIBLE_LEVELS, 
                "Expecting file->num_plausible_levels=%u == NUM_PLAUSIBLE_LEVELS=%u", file->num_plausible_levels, NUM_PLAUSIBLE_LEVELS);

        mutex_initialize (file->bgzf_discovery_mutex);
    }

    else if (IS_MGZIP(file->effective_codec)) {
        txt_file->mgzip_flags.level = BGZF_COMP_LEVEL_UNKNOWN; // don't write MGZIP section (15.0.80)

        // bug 1101: we don't yet know the plausible levels for other MGZIP codecs
        mgzip_set_is_exactable (txt_file, false, "Not BGZF"); 
    }

    else 
        ABORT ("Unsupported effective_codec=%s for discovery", codec_name (file->effective_codec));
}

// ZIP main thread
static void bgzf_discover_finalize_testing (MgzipLibraryType lib, MgzipLevel level)
{
    txt_file->mgzip_flags.library = lib;  // assign field-by-field to not modify other fields
    txt_file->mgzip_flags.level   = level;

    if (flag.zip_comp_i < MAX_NUM_COMPS) // for stats
        z_file->comp_bgzf[flag.zip_comp_i] = txt_file->mgzip_flags;

    mgzip_set_is_exactable (txt_file, level != BGZF_COMP_LEVEL_UNKNOWN, "BGZF library unrecognized"); // usually already set in mgzip_uncompress_one_block(), but setting again here (at end of zip) as need for short files, with less than BGZF_DISCOVERY_MAX_TESTS BGZF blocks to test
}

// ZIP main thread
void bgzf_finalize_discovery (void)
{    
    if (txt_file->is_exactable == no) return;

    for (FlagsMgzip *ll = txt_file->bgzf_plausible_levels; ll < &txt_file->bgzf_plausible_levels[NUM_PLAUSIBLE_LEVELS]; ll++)
        // case: at least one library⁀level was verified with all test bgzf blocks (10 blocks, unless file is shorter) 
        if (!ll->implausible) { // plausible level
            bgzf_discover_finalize_testing (ll->library, ll->level);

            if (flag_show_bgzf) 
                iprintf ("Discover: %s is a %s file, %s %s level %u\n", txt_name, codec_name (txt_file->effective_codec),
                         (txt_file->num_plausible_levels == 1) ? "identified as generated with" : "with multiple plausible levels, arbitrarily selecting",
                         bgzf_library_name (txt_file->mgzip_flags.library, true), txt_file->mgzip_flags.level);
            return;
        }    

    // case: there is no library⁀level for which we can decompress with bgzf=exact
    bgzf_discover_finalize_testing (0, BGZF_COMP_LEVEL_UNKNOWN); // has BGZF, but cannot identify level

    if (flag_show_bgzf) 
        iprintf ("Discover:%s: is a %s file, generated by an unidentified library\n", txt_name, codec_name (txt_file->effective_codec));
}

// ZIP: test a BGZF block against all the remaining plausible levels, and eliminate those that don't match. 
static void bgzf_discover_library_and_level (VBlockP vb, STRp(comp), STRp(uncomp))
{   
    mutex_lock (txt_file->bgzf_discovery_mutex);

    if (!txt_file->num_plausible_levels) // already determined that no levels were discovered 
        goto done;

    if (num_bgzf_blocks_tested_for_level() == BGZF_DISCOVERY_MAX_TESTS)
        goto done; // tested all we need

    uint32_t header_len = mgzip_header_len[txt_file->effective_codec];
    ASSERT (header_len, "%s_HEADER_LEN missing in mgzip_header_len", codec_name (txt_file->effective_codec));

    if (comp_len <= header_len + GZIP_FOOTER_LEN) {
        for (int l=0; l < NUM_PLAUSIBLE_LEVELS; l++) 
            txt_file->bgzf_plausible_levels[l].implausible = true;

        txt_file->num_plausible_levels = 0;

        if (flag_show_bgzf) 
            iprintf ("%s: Block too small - could not identify compression library and level\n", txt_name);

        goto done;
    }

    // ignore the header and footer of the block
    comp     += header_len;
    comp_len -= header_len + GZIP_FOOTER_LEN;

    // compress with each of the remaining plausible levels - testing if the compression is identical to the actual
    uint32_t recomp_size = uncomp_len * 1.1 + 64 KB; // guessing the max compressed size in the worst case scenario of very bad compression
    char *recomp = MALLOC (recomp_size);

    for (FlagsMgzip *ll = txt_file->bgzf_plausible_levels; ll < &txt_file->bgzf_plausible_levels[NUM_PLAUSIBLE_LEVELS]; ll++) {
        if (ll->implausible) continue; // already determined in a previous gz block

        // for large test blocks, skip high compression levels which are not common anyway, as testing is too slow
        if (comp_len > 100 KB && ll->level >= 8)
            continue;

        uint32_t recomp_len = gz_deflate (vb, *ll, STRa(uncomp), recomp, recomp_size, true); // 0 if recomp_size is too small

        if (!str_issame (comp, recomp)) {
            ll->implausible = true;
            txt_file->num_plausible_levels--;
        }
        
        if (flag_show_bgzf) 
            iprintf ("Discover[%d]: library %s level %u: size_in_file=%u size_in_test=%u plausible=%s\n", 
                     num_bgzf_blocks_tested_for_level(), bgzf_library_name (ll->library, true), ll->level, comp_len, recomp_len, YN(!ll->implausible));
    }

    increment_relaxed (txt_file->num_bgzf_blocks_tested_for_level);

    FREE (recomp);

done:
    mutex_unlock (txt_file->bgzf_discovery_mutex);
}


//--------------------------------------------------------------------
// ZIP SIDE - decompress MGZIP-compressed file and prepare BGZF section
//--------------------------------------------------------------------

uint32_t mgzip_get_max_block_size (void)
{
    switch (txt_file->effective_codec) {
        case CODEC_BGZF : return 64 KB;
        case CODEC_IL1M : return 1  MB;
        case CODEC_IL4M : return 4  MB;
        case CODEC_EMFL : return txt_file->max_mgzip_isize; // determined during discovery
        default         : return 1;
    }
}

void inc_disk_gz_uncomp_or_trunc_(FileP file, uint64_t inc, FUNCLINE)
{
    add_relaxed (file->disk_gz_uncomp_or_trunc, inc);

    if (flag.show_gz_uncomp)
        iprintf ("%s:%u: disk_gz_uncomp_or_trunc + %"PRIu64"\t= %"PRIu64"\n", func, code_line, inc, file->disk_gz_uncomp_or_trunc); 
}

static bool is_header_prefix (bytes header, STR8p(prefix))
{
    for (uint32_t i=0; i < prefix_len; i++)
        if (prefix[i] != header[i] && prefix[i] != '\b')
            return false;

    return true;
}

static GzStatus mgzip_block_verify_header (FileP file, // option 1
                                           FILE *fp, STR8p(h), // option 2
                                           bool discovering, GzipFlags expected_flg, STR8p(prefix))
{
    if (file) {
        fp = (FILE *)file->file;
        h = B1ST8 (file->gz_data);
        h_len = file->gz_data.len32;
    }

    GzipHeader *gz_header = (GzipHeader *)h;

    // no data at all
    if (h_len == 0 && feof (fp)) {
        ASSERT0 (!discovering, "unexpected end of file when discovering");
        return GZ_EOF_WITHOUT_EOF_BLOCK; // end of file
    }
    
    // truncated mid-way through header
    if (h_len < GZ_HEADER_LEN) {
        if (discovering)
            return GZ_NOT_GZIP; // file smaller than a gzip header - its not GZIP
        
        else if (flag.truncate && feof (fp)) 
            return GZ_TRUNCATED;   // truncated file      
        
        else
            ABORT ("%s file %s appears truncated - it ends with a partial gzip block header. offset=%"PRIu64". If you expect this file to be truncated, use --truncate", 
                   codec_name (file->effective_codec), file->basename, (uint64_t)ftello64 (fp) - h_len); // less than the minimal gz block header size
    }

    // case: this is not a GZIP block at all (see: https://tools.ietf.org/html/rfc1952)
    else if (gz_header->id1 != 31 || gz_header->id2 != 139 || gz_header->cm != 8) {
        if (discovering)
            return GZ_NOT_GZIP;
        else
            ABORT ("expecting %s file %s to be compressed with gzip format, but it is not. offset=%"PRIu64, 
                   codec_name (file->effective_codec), file->basename, (int64_t)ftello64 (fp) - h_len);
    }

    // case: this GZIP block is NOT of the structure expected by the specific gz codec
    else if (expected_flg.value != gz_header->flg.value || // this works in discovery too, important for codecs with a prefix that is only known after discovery (EMVL, GZBL)
             !is_header_prefix (h, STRa(prefix))) {
        if (discovering)
            return GZ_IS_OTHER_FORMAT;
        else 
            RESTART ("--no-bgzf", "Encountered a GZIP block that unexpectedly is not %s in %s offset=%"PRIu64" file_size=%"PRIu64" found=%s expected=%s", 
                     codec_name (file->effective_codec), file->basename, (int64_t)ftello64 (fp) - h_len, file->disk_size, display_gz_header (h, h_len, false).s, display_gz_header (STRa(prefix), false).s);
    }

    return GZ_SUCCESS;
}

// read and uncompress the last two bgzf blocks of a file. Used for testing a BAM has unmapped aligments.
bool bgzf_read_and_uncomp_final_block (rom filename, qSTRp(uncomp))
{
    FILE *fp = fopen (filename, "rb");
    if (!fp) return false;

    uint8_t comp[BGZF_MAX_BLOCK_SIZE + BGZF_EOF_LEN];

    if (fseeko64 (fp, -(int)sizeof(comp), SEEK_END)) 
        goto fail; // note: fails if file is too short, but that's ok 

    int32_t comp_len = fread (comp, 1, sizeof(comp), fp);

    // search for BGZF block
    for (int32_t i=comp_len - BGZF_EOF_LEN - 1; i >= 0; i--) // block must contain at least 1 byte bigger than an empty EOF block
        if (comp[i] == 31 && comp[i+1] == 139 && comp[i+2] == 8 && 
            mgzip_block_verify_header (NULL, fp, comp+i, comp_len-i, true, (GzipFlags){ .extra=1 }, _8(BGZF_PREFIX)) == GZ_SUCCESS) {
            fclose (fp);

            struct libdeflate_decompressor *gz_inflate_mem = libdeflate_alloc_decompressor (evb, __FUNCLINE);

            size_t uncomp_len_size_t; 

            enum libdeflate_result ret = libdeflate_deflate_decompress (
                gz_inflate_mem, 
                comp+i + BGZF_HEADER_LEN, comp_len-i - BGZF_HEADER_LEN - GZIP_FOOTER_LEN, // compressed
                uncomp, *uncomp_len, &uncomp_len_size_t);  // uncompressed

            *uncomp_len = uncomp_len_size_t;

            libdeflate_free_decompressor (&gz_inflate_mem, __FUNCLINE);

            return (ret == LIBDEFLATE_SUCCESS);
        }

fail:
    fclose (fp);
    return false;
}

static void mgzip_show_truncated (FileP file, uint32_t comp_len_truncated)
{
    iprintf ("TRUNCATED  %s thread=%s comp_len_truncated=%u truncated incomplete final %s block\n", 
             codec_name (file->effective_codec), threads_am_i_main_thread() ? "MAIN" : "COMPUTE", comp_len_truncated, codec_name (file->effective_codec));
} 

#define ILxM_is_valid_isize(mb) \
static bool IL##mb##M_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool discovering, bool *is_end_of_vb) \
{   \
    /* normally Il1M blocks are exactly (mb) MB of txt data, but can be shorter: */ \
    bool allow_smaller_blocks = is_eof        /* the final block is always allowed to be shorter (inc. in discovery) */ \
                             || !discovering; /* all other blocks are also allowed to be shorter (happens when multiple ILxM files are concatenated) - except the first block during discovery (unless its the eof block) */ \
    \
    return proposed_isize == mb MB || (allow_smaller_blocks && proposed_isize < mb MB); \
}
ILxM_is_valid_isize(1);
ILxM_is_valid_isize(4);

static bool EMFL_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool discovering, bool *is_end_of_vb)
{
    bool allow_smaller_blocks = !discovering || is_eof; // see comment in IL1M

    return !file->max_mgzip_isize || // not set yet because this is the first block  
           proposed_isize == file->max_mgzip_isize || 
           (allow_smaller_blocks && proposed_isize < file->max_mgzip_isize);
}

static bool EMVL_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool discovering, bool *is_end_of_vb)
{
    return (proposed_isize < 512 MB); // sanity check
}

static bool GZBL_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool discovering, bool *is_end_of_vb)
{
    return (proposed_isize < 512 MB); // sanity check
}

static bool MGSP_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool discovering, bool *is_end_of_vb)
{
    if (proposed_isize > 64 MB) return false; // sanity

    // case: first gz block in this VB
    if (!file->num_mgsp_blocks_in_vb) 
        file->mgsp_vb_isize = proposed_isize;

    if (proposed_isize == file->mgsp_vb_isize || // gz block with same isize as first
        proposed_isize == 0                   || // EOF block is always accepted
        proposed_isize - file->mgsp_vb_isize <= file->num_mgsp_blocks_in_vb) {  // last gz block in a group allowed to be a bit longer, by the remainder of dividing the total isizes by the number of blocks: i.e. an integer from zero to the number of block minus one.

        file->num_mgsp_blocks_in_vb++;
        file->max_mgsp_blocks_in_vb = MAX_(file->max_mgsp_blocks_in_vb, file->num_mgsp_blocks_in_vb);
    }

    // proposed block is starting a new group - it will go into the next VB, not this one
    else {
        file->num_mgsp_blocks_in_vb = file->mgsp_vb_isize = 0; // initialize for next VB
        *is_end_of_vb = true;
    }

    return true;
}

// read a gz block of a codec that does not contain bsize in the gz header
// returns: discovering: GZ_SUCCESS, GZ_IS_OTHER_FORMAT
//          otherwise:   GZ_SUCCESS
GzStatus mgzip_read_block_no_bsize (FileP file, bool discovering, Codec codec)
{
    START_TIMER;

    struct { 
        bool valid_3_blocks_isize;  // the first 3 gz blocks are expected to have the same isize (used for discovery) 
        uint16_t gz_hdr_len;
        uint32_t max_bsize;         // an upper limit we set on compressed size (bsize) of a block based on observation (after discovery)
        bytes gz_hdr;               // fixed gz header
        IsValidSize is_valid_isize; // isize validation function
    } params;

    switch (codec) {
        #define SET_IF(c,valid_3,max_bsize,gz_hdr,gz_prefix_len) case CODEC_##c: params = (typeof(params)){ valid_3, ((gz_prefix_len) ? (gz_prefix_len) : mgzip_header_len[CODEC_##c]), max_bsize, (bytes)gz_hdr, c##_is_valid_isize }; break
        SET_IF (MGSP, true,  4  MB, MGSP_HEADER,     0);
        SET_IF (IL1M, true,  1  MB, ILxM_PREFIX,     ILxM_PREFIX_LEN);
        SET_IF (IL4M, true,  4  MB, ILxM_PREFIX,     ILxM_PREFIX_LEN);
        SET_IF (EMVL, false, 32 MB, EMVL_HEADER,     0);
        SET_IF (EMFL, true,  4  MB, file->gz_header, 0);
        SET_IF (GZBL, false, 32 MB, file->gz_header, 0);
        default: ABORT0 ("codec missing");
    }

    if (discovering && params.valid_3_blocks_isize)
        params.max_bsize *= 3; // read 3 blocks - txtfile_discover_specific_gz will use this data to verify they all have the same isize

    FILE *fp = (FILE *)file->file;
    file->gz_data.comp_len = file->gz_data.uncomp_len = 0; // init

    // top up gz_data to max_comp_size (or less if EOF)
    txtfile_fread (file, fp, NULL, (int32_t)params.max_bsize - (int32_t)file->gz_data.len32, &file->disk_so_far);
    
    GzStatus status = mgzip_block_verify_header (file, 0,0,0, discovering, 
                                                 (GzipFlags){}, // note: all "no bsize" flavors have flgs=0 (so far)
                                                 STRa(params.gz_hdr));
    if (status == GZ_EOF_WITHOUT_EOF_BLOCK) {
        if (!discovering) file->no_more_blocks = true;
        return GZ_SUCCESS;
    }

    if (status != GZ_SUCCESS) return status;
        
    // search for block size by beginning of next block, and an isize that "makes sense" relative to bsize
    // note: we do this even if EOF, because gz_data might contain several gz blocks. note: also NULL if data is too short.
    uint8_t *next_blk = B8(file->gz_data, params.gz_hdr_len) - 1;
    uint32_t bsize=0, isize=0, n_blks=0;
    do {
        next_blk = memmem (next_blk + 1, BAFT8(file->gz_data) - next_blk, params.gz_hdr, params.gz_hdr_len);
        n_blks++;
    }
    // if codec is variable length (is_valid_isize doesn't work well), we take the extra precaution of verifying its bsize makes sense  
    // mgzip "no_bsize" is only used for FASTQ for which the file is quite uniform, so we can bind the ratio relatively narrowly
#   define RATIO_MULT 1.7
    while (next_blk &&  
           (({  isize = GET_UINT32 (next_blk - 4);
                bsize = BNUM (file->gz_data, next_blk);
                double ratio = (double)isize / (double)bsize;
                if (discovering && flag.show_bgzf) 
                    iprintf ("Discovery %s: isize=%u bsize=%u ratio=%1.1f ∈ [%1.1f, %1.1f]\n", 
                             codec_name (codec), isize, bsize, ratio, segconf.gz_comp_ratio / RATIO_MULT, segconf.gz_comp_ratio * RATIO_MULT);   
                
                // case: while attempting to skip bogus gz headers, we inadvertedly skipped a real one (because didn't match the expected format or ratio)
                // symptom: ratio just gets worse and worse with every block we read (bsize keeps growing but not isize), 
                // and if we don't stop, we will read he entire file
                if (ratio < 0.95 && !discovering)
                    RESTART ("--no-bgzf", "Failed to identify next gz-block's gz header, marking the close of this block. file=%s codec=%s offset=%"PRId64" file_size=%"PRIu64, 
                             file->basename, codec_name (codec), (int64_t)ftello64 ((FILE *)file->file) - file->gz_data.len, file->disk_size);
                            
                // final while condition: loop back to extend the block, if the ratio doesn't make sense, but not if bsize is too small to judge ratio
                bsize > 2 KB && 
                (ratio < segconf.gz_comp_ratio / RATIO_MULT || ratio > segconf.gz_comp_ratio * RATIO_MULT); })));
    
    bool is_end_of_vb = false;

    // case: a block was found, and it is not the last block
    if (next_blk && params.is_valid_isize (file, isize, false/* there IS a next block so not EOF*/, discovering, &is_end_of_vb)) {
        
        file->gz_data.uncomp_len = isize;    // isize is the last field of the block ending before next_blk
        file->gz_data.comp_len   = bsize;
    }

    // case: remaining data could be a final gz block, or could be a truncated block, 
    // we will know for sure when trying to uncompress it
    else if (!next_blk && feof (fp) && 
             !(discovering && n_blks == 1) && // discovery: likely this is a single gz-block file - not MGZIP at all (uncompress in main, to allow correct divvying up to VBs)
             file->gz_data.len32 >= params.gz_hdr_len + GZIP_FOOTER_LEN &&
             params.is_valid_isize (file, GET_UINT32 (BAFT8 (file->gz_data) - 4), true, discovering, &is_end_of_vb)) {
        
        file->gz_data.uncomp_len = GET_UINT32 (BAFT8(file->gz_data) - 4);
        file->gz_data.comp_len = file->gz_data.len32;
    }
    
    // case: data in gz_data is does not contain a gz block - either not the right codec file or is truncated
    else {
        if (discovering)
            return GZ_IS_OTHER_FORMAT;

        // data is not e.g. IL1M somewhere in the middle of the file...
        if (!feof (fp)) 
            RESTART ("--no-bgzf", "Encountered a GZIP block that unexpectedly is not %s in %s offset=%"PRIu64" file_size=%"PRIu64, 
                     codec_name (codec), file->basename, (uint64_t)ftello64 ((FILE *)file->file) - file->gz_data.len, file->disk_size);

        // case: final data in file is not a full gz block and truncation allowed: 
        // account and then ignore the data that will not be gz-decompressed 
        if (flag.truncate) {
            WARN (_FYI "%s is truncated - its final %s block in incomplete. Dropping final %u bytes of the GZ data.", 
                  txt_name, codec_name (codec), file->gz_data.len32);

            if (flag_show_bgzf) mgzip_show_truncated (file, file->gz_data.len32);

            mgzip_set_is_exactable (file, false, "File is truncated");

            inc_disk_gz_uncomp_or_trunc (file, file->gz_data.len);
            file->gz_data.len32 = file->gz_data.uncomp_len = 0;
            segconf.zip_txt_modified = true;
            file->no_more_blocks = true;
        }

        else 
            ABORTINP ("%s is truncated mid-way through %s block. " _TIP "If this is expected, use --truncate to discard the final partial %s block", 
                      txt_name, codec_name (codec), codec_name (codec));
    }

    COPY_TIMER_EVB (mgzip_read_block_no_bsize);

    return (is_end_of_vb && !discovering) ? GZ_SUCCESS_END_OF_VB : GZ_SUCCESS;
}

static GzStatus bgzf_mgzf_set_block_lens (FileP file, uint32_t bsize, bool discovering)
{ 
    FILE *fp = (FILE *)file->file;

    if (file->gz_data.len32 >= bsize) {
        file->gz_data.comp_len   = bsize;
        file->gz_data.uncomp_len = GET_UINT32 (B8(file->gz_data, bsize-4)); 
        return GZ_SUCCESS;
    }

    if (discovering)
        return GZ_NOT_GZIP; 
    
    if (flag.truncate && feof (fp)) 
        return GZ_TRUNCATED;    // truncated file      

    int save_errno = errno; // we want to report errno of fread, not ftell.

    ABORT ("%s %s (ftell=%"PRId64" err=\"%s\" gz_data.len=%u but expecting=%u filesystem=%s). %s\n", 
            feof (fp) ? "Unexpected end of file while reading" : "Failed to read file", 
            file->basename, ftello64 (fp), 
            (file->is_remote && save_errno == ESPIPE) ? "Disconnected from remote host" : strerror (save_errno),
            file->gz_data.len32, bsize, arch_get_txt_filesystem().s,
            feof (fp) ? _TIP "If the file is expected to be truncated, you use --truncate to disregard the final partial BGZF block." : "");
}

static void bgzf_mgzf_verify_eof_block (FileP file, STRp(eof_block))
{
    // case: valid EOF block
    if (str_issame_(B1STc(file->gz_data), file->gz_data.comp_len, eof_block, eof_block_len))
        file->num_EOF_blocks++;
 
    // case: an isize=0 block that is not the EOF block (there are multiple representations of an empty block in deflate)
    // example: special.non-EOF-zero-bgzf-block.vcf.gz
    else {
        mgzip_set_is_exactable (file, false, "Non-EOF zero block");
    
        if (flag_show_bgzf)
            iprintf ("DETECTED non-EOF zero block - --bgzf=exact will not be available: %s\n", 
                     str_to_hex (B1ST8(file->gz_data), file->gz_data.comp_len).s);
    }
}

static bool mgzf_get_bsize (FileP file, uint32_t *bsize) 
{
    if (*bsize) 
        return true; // already set

    if (file->gz_data.len32 < MGZF_PREFIX_LEN + sizeof (uint32_t)) 
        return false; 

    if (memcmp (B1ST8 (file->gz_data), MGZF_PREFIX, MGZF_PREFIX_LEN))
        return false;

    MgzfHeader *h = B1ST(MgzfHeader, file->gz_data);
    *bsize = LTEN32 (h->bsize); 

    ASSERT (*bsize <= 128 MB, "bsize=%u seems too big", *bsize); // sanity

    return true;
}

// ZIP: reads and validates a MGZF block
// returns: discoverying: GZ_SUCCESS, GZ_NOT_GZIP, GZ_IS_OTHER_FORMAT
//          otherwise:    GZ_SUCCESS, GZ_EOF_WITHOUT_EOF_BLOCK, GZ_TRUNCATED
static GzStatus mgzf_read_block_do (FileP file, // txt_file is not yet assigned when called from txtfile_discover_specific_gz
                                    bool discovering)
{
    #define MGZF_CHUCK_SIZE ((uint32_t)(16 MB)) // max amount we read from disk at a time  
    
    FILE *fp = (FILE *)file->file;
    file->gz_data.comp_len = file->gz_data.uncomp_len = 0; // init
    uint32_t bsize = 0;

    // top-up if needed - in rare cases - twice (this happens very large block where gz_size>MGZF_CHUCK_SIZE and not known yet - first read of MGZF_CHUCK_SIZE includes the header)
    for (int i=0; i < 2; i++) 
        if ((!mgzf_get_bsize (file, &bsize) || bsize > file->gz_data.len32) && !feof (fp)) {
            int32_t chunk_size = MAX_(MGZF_CHUCK_SIZE, bsize);

            txtfile_fread (file, fp, NULL, chunk_size - (int32_t)file->gz_data.len32, &file->disk_so_far);
        }

    GzStatus status = mgzip_block_verify_header (file, 0,0,0, discovering, (GzipFlags){ .extra=1, .comment=1 }, _8(MGZF_PREFIX));
    if (status != GZ_SUCCESS) return status;
    
    status = bgzf_mgzf_set_block_lens (file, bsize, discovering); 

    if (status == GZ_SUCCESS) {
        // verify comment, length=8, format "C001R015"
        MgzfHeader *h = B1ST(MgzfHeader, file->gz_data);
        #define c h->comment
        bool is_valid_comment = strnlen (c, 9) == 8 &&
                                c[0]=='C' && IS_DIGIT(c[1]) && IS_DIGIT(c[2]) && IS_DIGIT(c[3]) &&
                                c[4]=='R' && IS_DIGIT(c[5]) && IS_DIGIT(c[6]) && IS_DIGIT(c[7]);
        #undef c
        if (discovering && !is_valid_comment) return GZ_IS_OTHER_FORMAT;

        if (h->bsize != MGZF_EOF_LEN && !is_valid_comment)
            RESTART ("--no-bgzf", "Invalid MGZF comment: gz_header={ %s }", 
                     display_gz_header ((bytes)STRb(file->gz_data), false).s);

        // if this is an isize=0 block, verify that it is an EOF block, else we won't be able reconstruct exact
        if (!file->gz_data.uncomp_len)
            bgzf_mgzf_verify_eof_block (file, MGZF_EOF, MGZF_EOF_LEN);
    }

    return status;
}

// ZIP: reads and validates a BGZF block
// returns: discoverying: GZ_SUCCESS, GZ_NOT_GZIP, GZ_IS_OTHER_FORMAT
//          otherwise:    GZ_SUCCESS, GZ_EOF_WITHOUT_EOF_BLOCK, GZ_TRUNCATED
static GzStatus bgzf_read_block_do (FileP file, // txt_file is not yet assigned when called from txtfile_discover_specific_gz
                                    bool discovering)
{
    FILE *fp = (FILE *)file->file;
    file->gz_data.comp_len = file->gz_data.uncomp_len = 0; // init

    // top-up if needed 
    if (file->gz_data.len32 < BGZF_MAX_BLOCK_SIZE && !feof (fp)) {
        int32_t chunk_size = txtfile_uncompress_mgzip_at_read() ? 150 KB // a bit more than the default block-device read-ahead buffer (128KB) for best parallelization between disk read-ahead and CPU decompression
                                                                : BGZF_MAX_CHUCK_SIZE; // bigger block is faster if we are prepared to yield the CPU when waiting for the disk 
        txtfile_fread (file, fp, NULL, chunk_size - (int32_t)file->gz_data.len32, &file->disk_so_far);
    }

    uint32_t bsize = (uint32_t)LTEN16 (B1ST(BgzfHeader, file->gz_data)->bsize) + 1;

    GzStatus status = mgzip_block_verify_header (file, 0,0,0, discovering, (GzipFlags){ .extra=1 }, _8(BGZF_PREFIX));
    if (status != GZ_SUCCESS) return status;

    status = bgzf_mgzf_set_block_lens (file, bsize, discovering); 

    if (status == GZ_SUCCESS) {
        ASSERT (file->gz_data.uncomp_len <= 65536, "isize=%u ∉ [0,65536] in %s offset=%"PRIu64, file->gz_data.uncomp_len, file->basename, (uint64_t)ftello64 (fp) - file->gz_data.len32);

        // if this is an isize=0 block, verify that it is an EOF block, else we won't be able reconstruct exact
        if (!file->gz_data.uncomp_len)
            bgzf_mgzf_verify_eof_block (file, _S(BGZF_EOF));
    }

    return status;
}

// ZIP main thread: reads a BGZF block into gz_data 
GzStatus mgzip_read_block_with_bsize (FileP file, bool discovering, Codec codec)
{
    START_TIMER;

    // with BGZF, gz_data is either empty or contains exactly 1 bgzf block
    if (file->gz_data.comp_len) goto success; // we already have 1 block

    GzStatus ret = (codec == CODEC_BGZF) ? bgzf_read_block_do (file, discovering)
                                         : mgzf_read_block_do (file, discovering);
    switch (ret) {
        case GZ_SUCCESS: success: // successful read of a BGZF block, or using a block already read during discovery             
            ret = GZ_SUCCESS;
            break;

        case GZ_NOT_GZIP: 
        case GZ_IS_OTHER_FORMAT:
            ASSERT (discovering, "ret=%d expected only if discovering", ret);
            break; // file->gz_data contains data that is not BGZF data

        case GZ_EOF_WITHOUT_EOF_BLOCK: // file ended without EOF block: that's fine
            ret = GZ_SUCCESS;
            break; // note: if file was not entirely read, we will detect that at the end of zip_one_file
            
        case GZ_TRUNCATED: // file ended mid-way through a BGZF block
            ASSERT0 (!discovering, "GZ_TRUNCATED unexpected when discovering");

            // case: truncation allowed: account and then discard the data that will not be gz-decompressed 
            if (flag.truncate) { 
                WARN (_FYI "%s is truncated - its final BGZF block in incomplete. Dropping final %u bytes of the GZ data.", txt_name, file->gz_data.len32);

                mgzip_set_is_exactable (file, false, "File is truncated");
                if (flag_show_bgzf) mgzip_show_truncated (file, file->gz_data.len32);

                inc_disk_gz_uncomp_or_trunc (file, file->gz_data.len);
                file->gz_data.len32 = file->gz_data.comp_len = file->gz_data.uncomp_len = 0; // discard partial BGZF block
                segconf.zip_txt_modified = true;

                ret = GZ_SUCCESS;
                break;
            }   
            
            else
                ABORTINP ("%s is truncated mid-way through BGZF block. " _TIP "If this is expected, use --truncate to discard the final partial BGZF block", txt_name);

        default:
            ABORT ("Unexpected ret=%s", gzstatus_name (ret));
    }

    if (!discovering)
        file->no_more_blocks = (file->gz_data.comp_len == file->gz_data.len32 && feof ((FILE*)file->file));

    COPY_TIMER_EVB (mgzip_read_block_with_bsize);
    return ret;
}

// uncompresses a BGZF block in vb->comp_txt_data referred to by bb, into its place in vb->txt_data as prescribed by bb
// might be called from main thread or compute threads
// caller is responsible for allocating and freeing vb->gz_deflate_mem
void mgzip_uncompress_one_block (VBlockP vb, GzBlockZip *bb, Codec codec)
{
    if (bb->is_uncompressed) return; // already uncompressed - nothing to do 
    ASSERTNOTNULL (vb->gz_inflate_mem);

    int header_len = mgzip_header_len[codec]; // note: EOF block's header length maybe different, eg in MGZF

    ASSERT (header_len, "%s_HEADER_LEN missing in mgzip_header_len[CODEC_%s]", codec_name (codec), codec_name (codec));

    uint8_t *h = B8(vb->comp_txt_data, bb->gz_index);

    // verify that entire block is within vb->comp_txt_data
    ASSERT (bb->gz_index + header_len < vb->comp_txt_data.len && // we have at least the header - we can access bsize
            bb->gz_index + bb->gz_size <= vb->comp_txt_data.len, 
            "%s: %s block size goes past the end of in vb->comp_txt_data: bb=%s gz_index=%u vb->comp_txt_data.len=%"PRIu64, 
            VB_NAME, codec_name (codec), display_bb (bb).s, bb->gz_index, vb->comp_txt_data.len);

    ASSERT (h[0]==31 && h[1]==139 && h[2]==8, "%s: invalid %s block in vb->comp_txt_data: gz_index=%u", VB_NAME, codec_name (codec), bb->gz_index);

    // possibly grow txt_data: can happen data is MGZIP and its length exceeds vb_size due to last block going over
    buf_alloc (vb, &vb->txt_data, 0/*don't use "more" bc Ltxt already incremented*/, bb->txt_index + bb->txt_size + TXTFILE_READ_VB_PADDING, char, 0, "txt_data");

    uint32_t this_block_digest=0;
    if (txt_file->is_exactable || flag_show_bgzf) {
        this_block_digest = (uint32_t)XXH3_64bits (h, bb->gz_size); // 32 LSb
        bb->gz_digest     = this_block_digest;
    }

    enum libdeflate_result ret = bb->txt_size  
        ? libdeflate_deflate_decompress (vb->gz_inflate_mem, 
                                         h + header_len, bb->gz_size - header_len - GZIP_FOOTER_LEN, // compressed
                                         Btxt (bb->txt_index), bb->txt_size, NULL)  // uncompressed
        : LIBDEFLATE_SUCCESS; // don't uncompress if empty block: 1. not needed 2. header size of EOF block might differ (as in MGZF) so arithmetic will be wrong

    // account for the case of decompression, and also the case bb is discarded due to a certain truncate situation (see below).
    inc_disk_gz_uncomp_or_trunc (txt_file, bb->gz_size);

    // case: wrong isize, likely becuase the gz block is actually multiple gz blocks - we missed the header 
    // of the pervious block(s) (e.g. because gz compression was outside of ratio, or header differed from first header)
    // note: normally we catch this in mgzip_read_block_no_bsize, but there could be edge cases where bsize/isize ratio appears ok and we arrive here
    while (ret == LIBDEFLATE_INSUFFICIENT_SPACE && !TXT_GZ_HEADER_HAS_BSIZE)
        RESTART ("--no-bgzf", "Main thread mis-divided gz blocks. codec=%s bb=%s gz_index=%u vb->comp_txt_data.len=%"PRIu64, 
                 codec_name (txt_file->effective_codec), display_bb (bb).s, bb->gz_index, vb->comp_txt_data.len);

    // case: final ILxM block, which is truncated, but we have --truncate, and the garbage last word
    // unluckily < xMB so it went undetected as a legimiate block in ILxM_is_valid_isize. we drop this block now.
    if (ret != LIBDEFLATE_SUCCESS && (TXT_IS(IL1M) || TXT_IS(IL4M)) && bb->is_eof) {
        if (flag.truncate) {
            // receive updates made by main thread to mgzip_isizes, mgzip_starts: no more are going to happen as we reached eof
            __atomic_thread_fence (__ATOMIC_ACQ_REL); 

            txt_file->mgzip_isizes.len--; // remove truncated block from isizes
            txt_file->mgzip_starts.len--;
            mgzip_show_truncated (txt_file, bb->gz_size);
            return; // with bb->is_uncompressed=false
        }

        else {
            ABORT ("Failed to uncompress the final %s block of the file: %s. " _TIP "If it is expected that the file is truncated, use --truncate to ignore the defective final block.", 
                   codec_name (vb->txt_codec), libdeflate_error(ret));
        }
    }

    ASSERT (ret == LIBDEFLATE_SUCCESS, "libdeflate_deflate_decompress failed: %s", libdeflate_error(ret));

    bb->is_uncompressed = true;

    if (flag_show_bgzf)
        iprintf ("UNCOMPRESS %s thread=%-7s%-11s digest=%08x bb_i=%-3u comp_index=%-9u comp_len=%-8u txt_index=%-9u txt_len=%-8u eof=%-5s%s%s\n",
                 codec_name (codec), threads_am_i_main_thread() ? "MAIN" : "COMPUTE", 
                 cond_str (vb->vblock_i, " vb=", VB_NAME), this_block_digest, 
                 BNUM (vb->gz_blocks, bb), bb->gz_index, bb->gz_size, bb->txt_index, bb->txt_size, TF(bb->is_eof),
                 cond_str (bb->txt_size, " uncomp[5]=", str_to_printable_(Btxt(bb->txt_index), MIN_(5, Ltxt - bb->txt_index)).s),
                 bb->gz_size == BGZF_EOF_LEN ? " EOF" : "");

    // discover which gzip library and compression level were used (testing the first few BGZF blocks)
    if (num_bgzf_blocks_tested_for_level() < BGZF_DISCOVERY_MAX_TESTS && txt_file->num_plausible_levels) { // fail fast: no thread issue - if we are finished testing in another thread but num_bgzf_blocks_tested_for_level or num_plausible are not visible yet - we harmlessly call bgzf_discover_library_and_level again   
        bgzf_discover_library_and_level (vb, (rom)h, bb->gz_size, Btxt (bb->txt_index), bb->txt_size);   
        
        if ((!txt_file->num_plausible_levels || num_bgzf_blocks_tested_for_level() >= BGZF_DISCOVERY_MAX_TESTS)
        &&  !test_and_set_relaxed (txt_file->once_set_is_exactable)) // first threads wins
            mgzip_set_is_exactable (txt_file, txt_file->num_plausible_levels > 0, "BGZF library unrecognized");
    }
}

// ZIP compute thread: from zip_compress_one_vb
// SAM_SCAN / BAMASS compute thread:  pre-processing reading of a file (on a separate txt_file object)
void mgzip_uncompress_vb (VBlockP vb, Codec codec)
{
    START_TIMER;
    ASSERTNOTEMPTY (vb->gz_blocks);

    ASSERTISNULL (vb->gz_inflate_mem);
    vb->gz_inflate_mem = libdeflate_alloc_decompressor(vb, __FUNCLINE);

    uint32_t total_vb_isizes = 0;
    for_buf (GzBlockZip, bb, vb->gz_blocks) {
        mgzip_uncompress_one_block (vb, bb, codec);
        total_vb_isizes += bb->txt_size;
    }

    // sanity - total_vb_isizes is at least the size of the vb (can be more, bc the first/last bb can be shared with previous/next vb)
    ASSERT (total_vb_isizes >= Ltxt, "%s: Expecting total_vb_isizes=%u >= Ltxt=%u. codec=%s", 
            VB_NAME, total_vb_isizes, Ltxt, codec_name (txt_file->effective_codec));

    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gz_inflate_mem, __FUNCLINE); // also sets gz_inflate_mem to NULL

    buf_destroy (vb->comp_txt_data); // now that we are finished decompressing we can release the memory (we won't need this buffer for a while so better destroy)

    if (flag.show_time) {
        if (threads_am_i_main_thread ()) COPY_TIMER (bgzf_io_thread);
        else                             COPY_TIMER (bgzf_compute_thread);
    }

    COPY_TIMER (mgzip_uncompress_vb);
}

// ZIP: decompresses a prescribed BGZF block when re-reading DEPN lines
static inline void bgzf_uncompress_one_prescribed_block (VBlockP vb, STRp(bgzf_block), STRc(uncomp_block), uint64_t bb_i)
{
    START_TIMER;

    BgzfHeader *h = (BgzfHeader *)bgzf_block;

    if (flag_show_bgzf)
        iprintf ("REREAD  %s reread bb_i=%"PRIu64" gz_size=%u uncomp_size=%u ", 
                 VB_NAME, bb_i, bgzf_block_len, uncomp_block_len);

    enum libdeflate_result ret = 
        libdeflate_deflate_decompress (vb->gz_inflate_mem, 
                                       h+1, bgzf_block_len - BGZF_HEADER_LEN - GZIP_FOOTER_LEN, // compressed
                                       STRa(uncomp_block), NULL);  // uncompressed

    ASSERT (ret == LIBDEFLATE_SUCCESS, "%s: libdeflate_deflate_decompress failed: %s. bgzf_block_len=%u uncomp_block_len=%u bb_i=%"PRIu64, 
            VB_NAME, libdeflate_error(ret), bgzf_block_len, uncomp_block_len, bb_i);

    if (flag_show_bgzf)
        #define C(i) (i < uncomp_block_len ? char_to_printable (uncomp_block[i]).s : "") 
        iprintf ("txt_data[5]=%1s%1s%1s%1s%1s\n", C(0), C(1), C(2), C(3), C(4));
        #undef C

    COPY_TIMER (bgzf_uncompress_one_prescribed_block);
}

// ZIP: re-reads and validates one BGZF block
static void bgzf_reread_one_prescribed_block (FILE *fp, uint64_t offset, qSTRp (bgzf_block))
{
    ASSERT (!fseeko64 (fp, offset, SEEK_SET),
            "fseeko64(%s, %"PRIu64") failed while rereading BGZF depn lines: %s", txt_name, offset, strerror(errno));

    // read the header
    uint32_t header_bytes = txtfile_fread (txt_file, fp, bgzf_block, BGZF_HEADER_LEN, NULL);

    // failed to read as prescribed
    ASSERT (header_bytes == BGZF_HEADER_LEN && !memcmp (bgzf_block, BGZF_PREFIX, BGZF_PREFIX_LEN),
            "failed to re-read a BGZF block header as perscribed BGZF: offset=%"PRIu64" bytes_read=%u header=%s", offset, header_bytes, str_to_hex ((bytes)bgzf_block, header_bytes).s);
    
    uint32_t body_size = (LTEN16 (((BgzfHeader*)bgzf_block)->bsize) + 1) - BGZF_HEADER_LEN;
    uint32_t body_bytes = txtfile_fread (txt_file, fp, bgzf_block + BGZF_HEADER_LEN, body_size, NULL);

    ASSERT (body_bytes == body_size, "failed to re-read a BGZF block body as perscribed BGZF: offset=%"PRIu64" bytes_read=%u expected=%u", 
            offset, body_bytes, body_size);

    *bgzf_block_len = BGZF_HEADER_LEN + body_size; 
}

// ZIP: SAM/BAM: compute thread of a DEPN VB: actually re-reading data into txt_data according to vb->reread_prescription
void bgzf_reread_uncompress_vb_as_prescribed (VBlockP vb, FILE *fp)
{
    uint64_t last_offset = -1LL;
    char uncomp_block[BGZF_MAX_BLOCK_SIZE];

    ASSERTISNULL (vb->gz_inflate_mem);
    vb->gz_inflate_mem = libdeflate_alloc_decompressor(vb, __FUNCLINE);

    for_buf (RereadLine, line, vb->reread_prescription) {
        
        // a line might span 1 or more BGZF blocks
        while (line->line_len) { 
            ASSERT (line->offset.bb_i < txt_file->mgzip_starts.len32, "Expecting bb_i=%"PRIu64" < mgzip_starts.len=%"PRIu64, 
                    (uint64_t)line->offset.bb_i, txt_file->mgzip_starts.len);

            uint64_t offset = *B64(txt_file->mgzip_starts, line->offset.bb_i);
            uint32_t isize  = *B32(txt_file->mgzip_isizes, line->offset.bb_i);

            if (offset != last_offset) {
                STRl (bgzf_block, BGZF_MAX_BLOCK_SIZE);

                bgzf_reread_one_prescribed_block (fp, offset, qSTRa(bgzf_block)); 
                bgzf_uncompress_one_prescribed_block (vb, STRa(bgzf_block), uncomp_block, isize, line->offset.bb_i);
            
                last_offset = offset;
            }

            uint32_t subline_len = MIN_(line->line_len, isize - line->offset.uoffset);
            memcpy (BAFTtxt, &uncomp_block[line->offset.uoffset], subline_len);
            Ltxt += subline_len;
            
            // if this line continues to next BGZF block - it starts from the beginning of that block, its remainder is subline_len shorter
            line->line_len -= subline_len;
            line->offset.bb_i++;
            line->offset.uoffset = 0;      
        }
    }

    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gz_inflate_mem, __FUNCLINE); // also sets gz_inflate_mem to NULL
}

void bgzf_libdeflate_1_7_initialize (void)
{
    libdeflate_set_memory_allocator_1_7 (gz_alloc, codec_free_do);
}

// ZIP: called by Seg to set the bgzf index of the next line
void mgzip_zip_advance_index (VBlockP vb, uint32_t line_len)
{
    if (!vb->gz_blocks.len) return; // no MGZIP blocks in this VB - all data came from "unconsumed_txt"

    vb->line_bgzf_uoffset += line_len;

    // udpate current_bb_i and bgzf_offset (note: line_len might span multiple bgzf blocks)
    GzBlockZip *bb;
    for (bb = B(GzBlockZip, vb->gz_blocks, vb->gz_blocks.current_bb_i); 
         vb->line_bgzf_uoffset && vb->line_bgzf_uoffset >= bb->txt_size; // note: careful to also terminate on the edge case that line_bgzf_uoffset==0 and in the final VB block bb->txt_size==0 
         bb++) 
      
        vb->line_bgzf_uoffset -= bb->txt_size; // index into the next BGZF block

    vb->gz_blocks.current_bb_i = BNUM(vb->gz_blocks, bb);
}

// ZIP: after reading data for a txt_header or VB, copy unconsumed gz_blocks to txt_file->unconsumed_mgzip_blocks
// The first block might be partially consumed.
int64_t mgzip_copy_unconsumed_blocks (VBlockP vb)
{
    START_TIMER;
    ASSERTISZERO (txt_file->unconsumed_mgzip_blocks.len32);

    if (!vb->gz_blocks.len) return 0; // not a BGZF-compressed file

    int32_t consumed = // amount of data in vb->gz_blocks that does NOT need to be copied to next VB bc it was consumed by this VB or the previous one 
        Ltxt +         // amount of data consumed by this VB
        vb->gz_blocks.consumed_by_prev_vb; // amount of data in first BGZF block was consumed by the previous VB

    ARRAY (GzBlockZip, bb, vb->gz_blocks);

    bool done = false;
    bool consumed_full_bgzf_blocks=false;
    int64_t compressed_size = 0;

    for (uint32_t i=0; i < bb_len; i++) {
        // if some of the BGZF blocks are not consumed (the first of them might be partially consumed) - move the blocks
        // to unconsumed_mgzip_blocks - to be moved to the next VB
        if (consumed - bb[i].txt_size < 0 && !done/*enter only once*/) {
            
            consumed_full_bgzf_blocks = (consumed == 0); // no partially-consumed block

            // block i might be partially consumed or not consumed at all, subsequent blocks are not consumed at all
            buf_append (evb, txt_file->unconsumed_mgzip_blocks, GzBlockZip, 
                        B(GzBlockZip, vb->gz_blocks, i), vb->gz_blocks.len32 - i, "txt_file->unconsumed_mgzip_blocks");

            txt_file->unconsumed_mgzip_blocks.consumed_by_prev_vb = consumed; // part of first BGZF block already consumed
            done = true;
        }
        else if (!done)
            compressed_size += bb[i].gz_size;

        consumed -= bb[i].txt_size;
    }

    // sanity check
    ASSERT (-consumed - (int32_t)txt_file->last_truncated_line_len == txt_file->unconsumed_txt.len32, "Expecting (-consumed)=%d - last_truncated_line_len=%u == unconsumed_txt.len=%u", 
            -consumed, txt_file->last_truncated_line_len, txt_file->unconsumed_txt.len32);

    // update bb.txt_index for next VB
    // note: first bb.txt_data of the next VB is possibly negative if some of its data was consumed by the current VB
    int32_t txt_index = -txt_file->unconsumed_mgzip_blocks.consumed_by_prev_vb;
    for_buf (GzBlockZip, bb, txt_file->unconsumed_mgzip_blocks) {
        bb->txt_index = txt_index; 
        txt_index += bb->txt_size;
    }

    COPY_TIMER (mgzip_copy_unconsumed_blocks);
    return consumed_full_bgzf_blocks ? compressed_size : 0;
} 

// return blocks used by the segconf VB to the unconsumed blocks
void mgzip_return_segconf_blocks (VBlockP vb)
{
    buf_copy (evb, &txt_file->unconsumed_mgzip_blocks, &vb->gz_blocks, GzBlockZip, 0, 0, 0);
    txt_file->unconsumed_mgzip_blocks.consumed_by_prev_vb = vb->gz_blocks.consumed_by_prev_vb;
}

// ZIP: before reading data for a VB, populate gz_blocks with some or all of the unconsumed blocks passed
// from the previous VB or txt_header
void mgzip_zip_init_vb (VBlockP vb)
{
    vb->vb_mgzip_i = txt_file->mgzip_isizes.len; // index of first mgzip block in vb->gz_blocks

    if (!txt_file->unconsumed_mgzip_blocks.len) return; // happens when either unconsumed_bytes=0 or not a BGZF-compressed file

    // data in the first MGZIP block already consumed by previous VB or txt_header 
    vb->gz_blocks.consumed_by_prev_vb = vb->line_bgzf_uoffset = txt_file->unconsumed_mgzip_blocks.consumed_by_prev_vb;

    // copy all unconsumed (or partially consumed) MGZIP blocks - we might not need all of them - the unconsumed ones will moved back in mgzip_copy_unconsumed_blocks
    buf_copy (vb, &vb->gz_blocks, &txt_file->unconsumed_mgzip_blocks, GzBlockZip, 0, 0, "gz_blocks"); 

    vb->vb_mgzip_i -= txt_file->unconsumed_mgzip_blocks.len32;

    txt_file->unconsumed_mgzip_blocks.len32 = txt_file->unconsumed_mgzip_blocks.consumed_by_prev_vb = 0;

    // sanity check
    int32_t available = -vb->gz_blocks.consumed_by_prev_vb; // possibly start negative
    for_buf (GzBlockZip, bb, vb->gz_blocks) 
        available += bb->txt_size;

    ASSERT (available >= Ltxt, "%s blocks in txt_file->unconsumed_mgzip_blocks cover only %d bytes, less than the needed unconsumed_bytes=%d", 
            codec_name (txt_file->effective_codec), available, Ltxt);
}

//-----------------------------------------------------
// PIZ SIDE - compressing txt_file with BGZF
//-----------------------------------------------------

static void bgzf_show_compress (VBlockP vb, int32_t gz_blk_i/*# of block in file*/, uint32_t comp_index, uint32_t comp_len, rom txt, uint32_t txt_index, uint32_t isize, uint32_t gz_digest/*32 LSBb*/)
{
    iprintf ("COMPRESS thread=%-7s %-7s gz_blk_i=%-3d comp_txt_data.i=%-8u bsize=%-5u txt_data.i=%-8d isize=%-5u gz_digest=%08x%s%s\n",
             threads_am_i_main_thread() ? "MAIN" : threads_am_i_writer_thread() ? "WRITER" : "COMPUTE", VB_NAME, gz_blk_i,
             comp_index, comp_len, txt_index, isize, gz_digest,
             cond_str (isize, " uncomp[5]=", str_to_printable_(txt, MIN_(isize, 5)).s),
             comp_len == BGZF_EOF_LEN ? " EOF" : "");
}

static void bgzf_failed_exact_verification (void)
{
    if (flag.bgzf == BGZF_EXACT_STRICT) 
        ABORTINP ("%s (technical details: library=%s) %s", 
                  NON_EXACT_ERROR, bgzf_lib_name_level (txt_file->mgzip_flags).s, WEBSITE_GZ);

    flag.bgzf = BGZF_EXACT_FAILED; // no point continuing to test + avoid full-file verification and error message
    
    WARN ("%s A different gz-recompression method was used.", NON_EXACT_ERROR);
}

static void bgzf_failed_exact_verification_one_block (VBlockP vb, const BgzfBlockPiz *restrict block, uint32_t this_block_digest, STRp(gz_data))
{
    if (!test_and_set_relaxed (txt_file->exact_failed_verification)) { // once per file
        uint64_t gz_blk_i = vb->vb_mgzip_i + BNUM64(vb->gz_blocks, block);
        if (flag_show_bgzf) 
            iprintf ("VERIFY exact gz-recompression FAILED for gz_blk_i=%"PRIu64" isize=%u bsize=%u zip_digest=%08x piz_digest=%08x library=%s\n", 
                     gz_blk_i, block->txt_size, gz_data_len, block->gz_digest, (uint32_t)this_block_digest, bgzf_lib_name_level (txt_file->mgzip_flags).s);
                
        if (flag.dump_gz_block == BGZF_EXACT_STRICT) {
            char dump_fn[128];
            snprintf (dump_fn, sizeof (dump_fn), "gz-block.piz.%"PRIu64".%s.gz", 
                    gz_blk_i, bgzf_lib_name_level (txt_file->mgzip_flags).s);
            
            if (file_put_data (dump_fn, STRa(gz_data), 0))
                fprintf (stderr, "\nDumped file %s\n", dump_fn);
        }

        bgzf_failed_exact_verification();
    }
}

// BGZF compute thread (in PIZ)
static void bgzf_compress_one_block (VBlockP vb, const BgzfBlockPiz *restrict block)
{
    START_TIMER;

    ASSERTNOTNULL (vb->gz_deflate_mem); // caller should allocate

    #define BGZF_MAX_CDATA_SIZE (BGZF_MAX_BLOCK_SIZE - BGZF_HEADER_LEN - GZIP_FOOTER_LEN)

    buf_alloc (vb, &vb->comp_txt_data, BGZF_MAX_BLOCK_SIZE, 0, char, 1.2, "comp_txt_data");
    uint32_t comp_index = vb->comp_txt_data.len32, recomp_len;

    BgzfHeader *header = (BgzfHeader *)BAFTc (vb->comp_txt_data);
    buf_add (&vb->comp_txt_data, BGZF_EOF, BGZF_HEADER_LEN); // template of header - only bsize needs updating

    recomp_len = gz_deflate (vb, txt_file->mgzip_flags, Btxt (block->txt_index), block->txt_size, BAFTc (vb->comp_txt_data), BGZF_MAX_CDATA_SIZE, true);

    // in case the compressed data doesn't fit in one BGZF block, move to compressing at the maximum level. this can
    // happen theoretically (maybe) if the original data was compressed with a higher level, and an uncompressible 64K block was
    // compressed to just under 64K while in our compression level it is just over 64K.
    if (!recomp_len) { 
        void *high_compressor = libdeflate_alloc_compressor (vb, LIBDEFLATE_MAX_LEVEL, __FUNCLINE); // libdefate's highest level
        recomp_len = libdeflate_deflate_compress (high_compressor, Btxt (block->txt_index), block->txt_size, BAFTc (vb->comp_txt_data), BGZF_MAX_CDATA_SIZE);
        libdeflate_free_compressor (high_compressor, __FUNCLINE);
    }

    ASSERT (recomp_len, "cannot compress block with %u bytes into a BGZF block with %u bytes", block->txt_size, BGZF_MAX_BLOCK_SIZE);
    vb->comp_txt_data.len32 += recomp_len;

    header->bsize = LTEN16 ((uint16_t)(BGZF_HEADER_LEN + recomp_len + GZIP_FOOTER_LEN - 1));

    GzipFooter footer = { .crc32 = LTEN32 (crc32 (0, Btxt(block->txt_index), block->txt_size)),
                          .isize = LTEN32 (block->txt_size) };
    buf_add (&vb->comp_txt_data, (rom)&footer, GZIP_FOOTER_LEN);

    recomp_len += BGZF_HEADER_LEN + GZIP_FOOTER_LEN;
    
    uint32_t this_block_digest=0;
    #define print_show_bgzf if (flag_show_bgzf) bgzf_show_compress (vb, BNUM(vb->gz_blocks, block), comp_index, recomp_len, Btxt(block->txt_index), block->txt_index, block->txt_size, this_block_digest)

    if (IS_PIZ_EXACT && txt_file->mgzip_digests.len > 0) {
        this_block_digest = (uint32_t)XXH3_64bits (Bc(vb->comp_txt_data, comp_index), recomp_len); // keep 32 LSb
        
        print_show_bgzf; // before checking for digest mismatch

        if (this_block_digest != block->gz_digest) 
            bgzf_failed_exact_verification_one_block (vb, block, this_block_digest, Bc(vb->comp_txt_data, comp_index), recomp_len);            
    }
    else
        print_show_bgzf;

    COPY_TIMER (bgzf_compress_one_block);
} 

// PIZ writer thread: appends file data to wvb->comp_txt_data
void bgzf_write_finalize (void)
{
    // exact verification for files up to 15.0.79. note: since 15.0.80, gz blocks are verified by block-level digests    
    if (IS_PIZ_EXACT && !VER2(15,80)) {
        bool success = (txt_file->disk_so_far % (16 MB) == txt_file->orig_gz_file_size // verification up to 15.0.79 was simply comparing the 3 LSB of the gz_file_size
                     || (Z_DT(FASTQ) && VER2(15,64))); // defect 2026-03-08 affecting FASTQ files 15.0.64-79 - we can't verify these files

        // case: gz file produced differs in size vs original file 
        if (!success) { // note: all blocks that were written had their gz_digest compared successfully, so a failure here can only happen due a genozip bug, e.g. blocks were omitted or write order switched 
            if (flag_show_bgzf) 
                iprintf ("VERIFY exact gz-recompression FAILED: gz_file_size=%"PRIu64" (24LSB==%"PRIu64") but original was gz_file_size.24LSB=%"PRIu64"\n", 
                         txt_file->disk_so_far, txt_file->disk_so_far % (16 MB), txt_file->orig_gz_file_size);
            
            bgzf_failed_exact_verification();
        }

        else if (flag_show_bgzf) // verified
            iprintf ("VERIFY exact gz-recompression SUCCEEDED: gz_file_size=%"PRIu64" 24LSB=%"PRIu64" same as original\n", 
                     txt_file->disk_so_far, txt_file->disk_so_far % (16 MB));
    }

    // if we're not using the gz blocks from mgzip_isizes, add EOF block now. 
    // note: if recomrpession-verification of a block failed, IS_PIZ_EXACT=false but we still continue to use mgzip_isizes and hence outputted an EOF already
    if (!txt_file->mgzip_isizes.len) {
        buf_add_more (wvb, &wvb->comp_txt_data, BGZF_EOF, BGZF_EOF_LEN, "comp_txt_data");

        if (flag_show_bgzf) 
            bgzf_show_compress (wvb, 0, wvb->comp_txt_data.len32, BGZF_EOF_LEN, NULL, 0, 0, 0);
    }
}

// Entry point of BGZF compression compute thread. This also calculates the BAI index for BAM files.
// bgzf-compress vb->txt_data into vb->comp_txt_data - using BGZF blocks as prescribed in vb->gz_blocks. 
// Note: we hope to reconstruct the exact same byte-level BGZF blocks, as the original files, but that 
// will only happen if the GZIP library (eg libdeflate), version and parameters are the same 
static void bgzf_compress_vb (VBlockP vb)
{
    START_TIMER;

    if (flag.make_bai && !vb->is_txt_header)
        bai_calculate_one_vb (vb);

    if (flag_show_bgzf)
        iprintf ("COMPRESS thread=%s %s initialized <%sexact> re-compression with %s(%s[%u])\n",
                 threads_am_i_main_thread() ? "MAIN" : threads_am_i_writer_thread() ? "WRITER" : "COMPUTE", VB_NAME,
                 IS_PIZ_EXACT ? "" : "non-",
                 codec_name (txt_file->effective_codec), bgzf_library_name (txt_file->mgzip_flags.library, true), txt_file->mgzip_flags.level);

    ASSERTNOTEMPTY (vb->gz_blocks);

    buf_alloc (vb, &vb->comp_txt_data, 0, vb->gz_blocks.len32 * BGZF_MAX_BLOCK_SIZE/2, uint8_t, 1, "comp_txt_data"); // alloc based on estimated size
    gz_alloc_deflator (vb, txt_file->mgzip_flags);

    for_buf2 (BgzfBlockPiz, block, bb_i, vb->gz_blocks) {
        ASSERT (block->txt_index + block->txt_size <= Ltxt, 
                "bb_i=%u out of range: expecting txt_index=%u txt_size=%u <= txt_data.len=%u",
                bb_i, block->txt_index, block->txt_size, Ltxt);

        block->gz_index = vb->comp_txt_data.len32;

        bgzf_compress_one_block (vb, block);
    }

    // add an "after" pseudo-block (e.g. for ease of looping in bai_update_vfo_to_bgzf_offsets)
    buf_append_one (vb->gz_blocks, (BgzfBlockPiz){ .gz_index = vb->comp_txt_data.len32 } );

    gz_free_deflator (vb, txt_file->mgzip_flags);

    vb_set_is_processed (vb); /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ 

    if (vb != evb) COPY_TIMER (bgzf_compute_thread);
}

#define BGZF_CREATED_BLOCK_SIZE 65280 // same size as observed in htslib-created files

// PIZ: calculate the BGZF blocks within this VB
static uint32_t bgzf_calculate_blocks_one_vb (VBlockP vb, bool is_last)
{
    // reconstruct --exact-ly based on mgzip_isizes
    // note: we continue to reconstruct with the prescribed gz_blocks even if a block fails verification
    if (txt_file->mgzip_isizes.len) {
        ARRAY (uint32_t, isize_arr,  txt_file->mgzip_isizes);
        ARRAY (uint32_t, digest_arr, txt_file->mgzip_digests);

        // move i to the index in mgzip_isizes after the gz_blocks needed for this vb
        uint32_t i, txt_index=0;
        for (i=txt_file->mgzip_isizes.next; i < txt_file->mgzip_isizes.len; i++)
            if (txt_index + isize_arr[i] <= Ltxt/*<= to include EOF in preceding VB */)
                txt_index += isize_arr[i];
            else
                break;

        buf_alloc_exact (vb, vb->gz_blocks, i - txt_file->mgzip_isizes.next, BgzfBlockPiz, "gz_blocks");

        // create vb->gz_blocks
        txt_index = 0; // index into txt_data of this vb
        for_buf (BgzfBlockPiz, blk, vb->gz_blocks) {
            *blk = (BgzfBlockPiz){ .txt_index = txt_index, 
                                   .txt_size  = isize_arr [txt_file->mgzip_isizes.next],
                                   .gz_digest = digest_arr ? digest_arr[txt_file->mgzip_isizes.next] : 0 };
            txt_file->mgzip_isizes.next++;
            txt_index += blk->txt_size;
        }

        int32_t remaining = Ltxt - txt_index;
        ASSERT (IN_RANGE(remaining, 0, BGZF_MAX_BLOCK_SIZE), "mgzip_isizes exhausted prematurely: remaining=%d", remaining); // if we have 65536 or more remaining, there should have been more isizes

        return remaining;
    }

    // create our own equal-isize blocks
    else {
        buf_alloc_exact (vb, vb->gz_blocks, ceill ((double)Ltxt / BGZF_CREATED_BLOCK_SIZE), BgzfBlockPiz, "gz_blocks");
        
        for_buf2 (BgzfBlockPiz, blk, i, vb->gz_blocks)
            *blk = (BgzfBlockPiz){ .txt_index = i * BGZF_CREATED_BLOCK_SIZE, .txt_size = BGZF_CREATED_BLOCK_SIZE };
        
        BLST(BgzfBlockPiz, vb->gz_blocks)->txt_size -= (BGZF_CREATED_BLOCK_SIZE * vb->gz_blocks.len - Ltxt); // remove excessive length from last block
    
        return 0; // our gz_blocks perfectly cover the VB - no data remaining
    }
}

// PIZ: Writer thread
void bgzf_dispatch_compress (Dispatcher dispatcher, STRp (uncomp), CompIType comp_i, bool is_last, bool is_txt_header)
{
    // uncompressed data to be dealt with by next call to this function (buffer belongs to writer thread)
    static Buffer intercall_txt = {}; // belongs to wvb
    buf_alloc (wvb, &intercall_txt, 0, BGZF_MAX_BLOCK_SIZE, char, 0, "intercall_txt");

    uint32_t next_isize = txt_file->mgzip_isizes.len ? *B32(txt_file->mgzip_isizes, txt_file->mgzip_isizes.next) 
                                                     : BGZF_CREATED_BLOCK_SIZE;

    // case: uncomp is not enough to fill a block, just store it to next call (but not if TXT_HEADER because txt_header and BAM data in the same VB will confuse BAI)
    if (!is_last && (uncomp_len + intercall_txt.len32 < next_isize) && !is_txt_header) {
        memcpy (BAFTc(intercall_txt), uncomp, uncomp_len);
        intercall_txt.len32 += uncomp_len;
        return;
    }

    if (uncomp_len ||          // might be 0 if is_last, in some cases
        intercall_txt.len) {   // non-zero only possible if --bgzf=exact

        VBlockP vb = dispatcher_generate_next_vb (dispatcher, wvb->vblock_i, COMP_NONE);
        vb->comp_i = comp_i;

        // build uncompressed data for this VB - some data left over from previous VB + data from wvb
        buf_alloc_exact (vb, vb->txt_data, intercall_txt.len + uncomp_len, char, "txt_data");
        if (intercall_txt.len32) memcpy (B1STtxt, intercall_txt.data, intercall_txt.len32);
        memcpy (Btxt (intercall_txt.len32), uncomp, uncomp_len);

        // calculate BGZF blocks - and trim data that doesn't fill a block - to be moved to next VB
        if ((intercall_txt.len32 = bgzf_calculate_blocks_one_vb (vb, is_last))) { // non-zero is possible only if bgzf=exact
            Ltxt -= intercall_txt.len32;
            memcpy (B1STc(intercall_txt), BAFTtxt, intercall_txt.len32);
        }

        vb->is_txt_header = is_txt_header; // not BAI indexing on txt_header VB

        // BGZF-compress vb->txt_data in a separate thread
        dispatcher_compute (dispatcher, bgzf_compress_vb);
    }

    if (is_last) {
        dispatcher_set_no_data_available (dispatcher, false, DATA_EXHAUSTED);
        buf_destroy (intercall_txt);
    }
}

// used by test/Makefile: stdin to stdout
void generate_il1m (void)
{
    FlagsMgzip mgzip_flags = { .library=BGZF_IGZIP, .level=3 };
    gz_alloc_deflator (evb, mgzip_flags);
    
    char *in = MALLOC (2 MB), *out = MALLOC (4 MB);
    uint32_t in_len;

    // note: if combined with --no-bgzf - 3rd block is artifically non-compliant
    for (int i=0; (in_len = fread (in, 1, 1 MB, stdin)); i++) {
        GzipFooter footer = { .crc32 = LTEN32 (crc32 (0, in, in_len)),
                              .isize = LTEN32 (in_len) };

        uint32_t out_len = gz_deflate (evb, mgzip_flags, in, in_len, out, 4 MB, false);

        
        GzipHeader h = { .id1=31, .id2=139, .cm=8, .xfl=flag.best?2 : flag.fast?4 : 0, .os=3 };
        if (i==2 && flag.no_bgzf) 
            h.mtime = 1971; // user requested that 3rd block is non-compliant to ILxM - we modify mtime

        ASSERT0 (1 == fwrite (&h, sizeof (h), 1, stdout), "fwrite failed #1");
        ASSERT  (1 == fwrite (STRa(out), 1, stdout), "fwrite failed: #3 out_len=%u", out_len);
        ASSERT0 (1 == fwrite (&footer, sizeof (footer), 1, stdout), "fwrite failed #4");
    }

    fflush (stdout);
}

// reads block and removes it from z_data
#define show_gz_next      B8(evb->z_data, evb->z_data.next)
#define show_gz_remaining (evb->z_data.len - evb->z_data.next)
static bool show_gz_uncompress_gz_block (rom filename, size_t *block_txt_len, size_t *block_gz_len)
{
    enum libdeflate_result ret;
    
    while ((ret = libdeflate_gzip_decompress_ex (evb->gz_inflate_mem, show_gz_next, show_gz_remaining, STRb(evb->txt_data), block_gz_len, block_txt_len))
            == LIBDEFLATE_INSUFFICIENT_SPACE)
            // double the size of txt_data if not enough txt space
            buf_alloc_exact (evb, evb->txt_data, evb->txt_data.len * 2, char, NULL); 
    
    if (ret == LIBDEFLATE_INSUFFICIENT_DATA 
     || ret == LIBDEFLATE_BAD_DATA) // can happen if we read the entire data and its stop abruptly because it is truncated
        return false; // out of data - we're done
        
    if (ret != LIBDEFLATE_SUCCESS) {
        printf ("%s: libdeflate_gzip_decompress_ex failed: %s (comp_len=%"PRIu64")\n", 
                filename, libdeflate_error(ret), evb->z_data.len - evb->z_data.next);
        return false;
    }

    return true;    
}

// called from --show-gz
void show_gz (rom filename) 
{
    // read a bunch of gz data from file
    bool is_analyzed = flag.explicit_quiet, is_emvl_empty_block = false;

    ASSINP (file_exists (filename), "Error: file not found: %s", filename);
    
    uint64_t file_size = file_get_size (filename), file_remaining=file_size, offset=0;
    ASSINP (file_size, "Error: file is empty: %s", filename);

    FILE *fp = fopen (filename, "rb");
    ASSERT (fp, "Failed to open %s: %s", filename, arch_str_error());

    buf_alloc (evb, &evb->z_data, 0, MIN_(100 MB, file_size), char, 0, "z_data");

    // file_get_file (evb, filename, &evb->z_data, "z_data", 100 MB, VERIFY_NONE, false);

    // uncompress first gz block of gz_data
    evb->gz_inflate_mem = libdeflate_alloc_decompressor (evb, __FUNCLINE);
    buf_alloc_exact (evb, evb->txt_data, 5 * evb->z_data.size, char, "txt_data");
    int size_width=0;

    int gz_blk_i=0; // sequential number of gz block in file
    while (file_remaining) {
        
        // move previously read final partial block to beginning of z_data
        if (evb->z_data.next) {
            memmove (B1ST8(evb->z_data), show_gz_next, show_gz_remaining);
            evb->z_data.len  = show_gz_remaining;
            evb->z_data.next = 0;
        }
        
        // top up z_data as much as possible
        uint64_t bytes_to_read = MIN_(evb->z_data.size - evb->z_data.len, file_remaining);
        ASSERT (fread (BAFT8(evb->z_data), bytes_to_read, 1, fp) == 1, "Failed to read %"PRIu64" bytes from %s\n", bytes_to_read, filename);
        evb->z_data.len += bytes_to_read;
        file_remaining -= bytes_to_read;

        // show all gz blocks currently in z_data 
        while (true) { 
            size_t block_txt_len, block_gz_len;
            
            // case: only partial block is available at end of z_data
            if (!show_gz_uncompress_gz_block (filename, &block_txt_len, &block_gz_len)) {
                bool is_full = (evb->z_data.len == evb->z_data.size && evb->z_data.next == 0);

                // case: we haven't read the entire file yet - top up for data and try again (unless entire data is full)
                if (file_remaining && !is_full)
                    break;

                if (is_full && !offset)
                    printf ("SINGLE GZ BLOCK (judging by first %s): gz_header=%s\n", 
                            str_size (evb->z_data.size).s, display_gz_header (show_gz_next, show_gz_remaining, false).s);

                else if (show_gz_remaining > 40 && !flag.explicit_quiet)
                    printf ("PARTIAL GZ BLOCK #%u offset=%"PRIu64": gz_header=%s\n", 
                            gz_blk_i, offset, display_gz_header (show_gz_next, show_gz_remaining, false).s);

                if (!gz_blk_i && !flag.explicit_quiet) {
                    if (show_gz_remaining >= STRLEN(MGSB_HEADER) && !memcmp (show_gz_next, MGSB_HEADER, STRLEN(MGSB_HEADER)))
                        printf ("Codec: Single gz-block MGI (based on header)\n"); // TODO: verify qname flavor
                    else
                        printf ("Cannot find a complete GZ block in the first %s - likely single-block GZIP\n", str_size (evb->z_data.len).s);
                }

                goto done;
            }

            if (!size_width && block_txt_len) 
                size_width = -str_get_uint_textual_len (block_txt_len);

            uint32_t h_len;
            StrText1K header_str = display_gz_header_ex (show_gz_next, block_gz_len, false, &h_len);

            printf ("i=%-2u: offset=%-10"PRIu64" isize=%*"PRIu64" (%-7s) bsize=%*"PRIu64" digest=%08x gz_header=%s.L=%u %s%s%s\n", 
                    gz_blk_i, offset, size_width, (uint64_t)block_txt_len, str_size (block_txt_len).s, size_width, (uint64_t)block_gz_len, 
                    (uint32_t)XXH3_64bits (show_gz_next, block_gz_len), // 32b digest
                    header_str.s, h_len,
                    (!block_txt_len && str_issame_((rom)show_gz_next, block_gz_len, BGZF_EOF, BGZF_EOF_LEN)) ? "⇐ BGZF_EOF" : "",
                    (!block_txt_len && str_issame_((rom)show_gz_next, block_gz_len, MGZF_EOF, MGZF_EOF_LEN)) ? "⇐ MGZF_EOF" : "",
                    (!block_txt_len && str_issame_((rom)show_gz_next, block_gz_len, MGSP_EOF, MGSP_EOF_LEN)) ? "⇐ MGSP_EOF" : "");
                    

            if (!block_txt_len) 
                is_emvl_empty_block = str_issame_(EMVL_FIRST_BLOCK, STRLEN(EMVL_FIRST_BLOCK), (rom)show_gz_next, block_gz_len);

            // analyze first block of at least 1 KB
            if (!is_analyzed && block_gz_len > 1 KB) {

                // analyze header (possibly multiple options)
                GzipHeader *h = (GzipHeader *)show_gz_next;
                bool is_bgzf = false;

                if (block_gz_len >= BGZF_PREFIX_LEN && !memcmp (h, BGZF_PREFIX, BGZF_PREFIX_LEN)) {
                    printf ("Codec: BGZF (based on header)\n");
                    is_bgzf = true;
                }

                if (block_gz_len >= ILxM_PREFIX_LEN && !memcmp (h, ILxM_PREFIX, ILxM_PREFIX_LEN) && 
                    h->os == 3 && (h->xfl==0 || h->xfl==2 || h->xfl==4)) {
                    if      (block_txt_len == 1 MB) printf ("Codec: IL1M (based on header and isize)\n");
                    else if (block_txt_len == 4 MB) printf ("Codec: IL4M (based on header and isize)\n");
                }

                if (block_gz_len > MGZF_PREFIX_LEN && !memcmp (h, MGZF_PREFIX, MGZF_PREFIX_LEN))
                    printf ("Codec: MGZF (based on header)\n");

                if (is_emvl_empty_block)
                    printf ("Codec: EMVL (based on first block being the EMVL empty block)\n");

                // TODO: MGSP, EMFL a bit more tricky to identify - would be helpful to get the qname flavor

                // detect library⁀level : recompress and compare
                FlagsMgzip ll[200] = {};
                int n_lls = 0;
                
                for (int l=0; l <= IGZIP_MAX_LEVEL; l++) // first because fastest
                    ll[n_lls++] = (FlagsMgzip){ .library = BGZF_IGZIP,        .level = l};

                for (int l=0; l <= LIBDEFLATE_MAX_LEVEL; l++) // level=0 only here, bc it would be the same in all libraries
                    ll[n_lls++] = (FlagsMgzip){ .library = BGZF_LIBDEFLATE19, .level = l};

                for (int l=0; l <= LIBDEFLATE_MAX_LEVEL; l++)
                    ll[n_lls++] = (FlagsMgzip){ .library = BGZF_LIBDEFLATE7,  .level = l};

                for (int l=0; l <= ZLIB_MAX_LEVEL; l++)
                    ll[n_lls++] = (FlagsMgzip){ .library = BGZF_ZLIB,         .level = l};

                ASSERTNOTINUSE (evb->scratch);
                buf_alloc_exact (evb, evb->scratch, evb->txt_data.len * 1.2, uint8_t, stringfy(evb->scratch)); // 1.2 in case deflate actually slightly grows rather than shrinks (unlikely edge cases)
                
                bytes payload = B8(evb->z_data, h_len);
                uint32_t payload_len = block_gz_len - h_len - sizeof(GzipFooter);

                #define eraser "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                #define spaces "                                                                                                                                                                "

                bool found=false;
                for (int i=0; i < n_lls; i++) {
                    // for large blocks, skip the super slow compression levels that are unlikely to be used by the software the geenerated these files
                    if (block_txt_len > 8 MB) {
                        if (ll[i].library==BGZF_LIBDEFLATE7  && ll[i].level >= 9)  continue; 
                        if (ll[i].library==BGZF_LIBDEFLATE19 && ll[i].level >= 10) continue; 
                        if (ll[i].library==BGZF_ZLIB         && ll[i].level >= 7)  continue; 
                    }

                    int len = printf ("Testing %-9s level=%-2u  ", bgzf_library_name (ll[i].library, false), ll[i].level);
                    fflush (stdout);

                    uint32_t test_bsize = gz_deflate (evb, ll[i], B1STc(evb->txt_data), block_txt_len, STRb(evb->scratch), is_bgzf);

                    if (str_issame_((rom)STRa(payload), B1STc(evb->scratch), test_bsize)) {
                        printf ("Plausible\n"); 
                        found = true;
                    }
                    else 
                        printf ("%.*s%.*s%.*s", len, eraser, len, spaces, len, eraser);
                    fflush (stdout);
                }

                if (!found) printf ("This file not compressed with any of the library⁀levels tested\n");

                is_analyzed = true;
            }

            gz_blk_i++;
            evb->z_data.next += block_gz_len;
            offset += block_gz_len;
            buf_free (evb->scratch);  
        }
    }

done:
    buf_free (evb->z_data);
    buf_free (evb->txt_data);
    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&evb->gz_inflate_mem, __FUNCLINE);
}


// called from --dump_gz_block=𝑛 
void dump_gz_block (rom filename) 
{
    ASSINP (file_exists (filename), "Error: file not found: %s", filename);
    
    uint64_t file_size = file_get_size (filename);
    ASSINP (file_size, "Error: file is empty: %s", filename);

    FILE *fp = fopen (filename, "rb");
    ASSERT (fp, "WARNING: Failed to open %s. fopen: %s", filename, strerror (errno));

    char data[64 KB];
    BgzfHeader *h = (BgzfHeader *)data;

    uint64_t offset=0;
    for (int i=0; i <= flag.dump_gz_block; i++) {
        fseeko64 (fp, offset, SEEK_SET);
        ASSERT (fread (h, BGZF_HEADER_LEN, 1, fp) == 1, "fread of block %u failed: %s", i, strerror (errno));
        
        ASSERT (!memcmp (h, BGZF_PREFIX, BGZF_PREFIX_LEN), "block %u is not a BGZF gzip block", i);
        offset += h->bsize + 1;
    }

    ASSERT (fread (&data[BGZF_HEADER_LEN], (h->bsize+1) - BGZF_HEADER_LEN, 1, fp) == 1, 
                   "fread of block %u failed (2): %s", flag.dump_gz_block, arch_str_error());
 
    char dump_fn[128];
    snprintf (dump_fn, sizeof (dump_fn), "gz-block.zip.%u.gz", flag.dump_gz_block);

    ASSERT (file_put_data (dump_fn, data, h->bsize+1, 0), "Failed to write file %s", dump_fn);

    fprintf (stderr, "\nDumped file %s\n", dump_fn);
    exit (0);
}


// Main thread: compress TBI: evb->txt_data to evb->comp_txt_data.
// NOTE: destroys txt_data fields, so must be called only after VCF is fully written
void bgzf_compress_tbi (void)
{
    START_TIMER;
    SAVE_FLAGS;

    txt_file->mgzip_isizes.len = 0; // not "exact" bgzf-blocks
    txt_file->mgzip_flags = bgzf_recompression_levels[3];

    bgzf_calculate_blocks_one_vb (evb, true);

    flag.make_bai  = false;
    flag.show_bgzf = false;
    bgzf_compress_vb (evb);

    buf_free (evb->gz_blocks);

    RESTORE_FLAGS;
    COPY_TIMER_EVB (bgzf_compress_tbi);
}
