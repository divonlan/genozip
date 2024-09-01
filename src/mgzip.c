// ------------------------------------------------------------------
//   mgzip.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
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
#include "mgzip.h"
#include "arch.h"
#include "file.h"
#include "zfile.h"
#include "zip.h"
#include "arch.h"
#include "txtfile.h"
#include "codec.h"
#include "threads.h"
#include "dispatcher.h"
#include "writer.h"
#include "gencomp.h"
#include "filename.h"
#include "strings.h"

#define LIBDEFLATE_MAX_LEVEL 12
#define ZLIB_MAX_LEVEL 9

// all data in Little Endian. Defined in https://datatracker.ietf.org/doc/html/rfc1952 and https://samtools.github.io/hts-specs/SAMv1.pdf
typedef struct __attribute__ ((packed, aligned(2))) BgzfHeader { // 18 bytes
    uint8_t id1;    // Gzip id - must be 31  (0x1f)
    uint8_t id2;    // Gzip id - must be 139 (0x8b)
    uint8_t cm;     // Compression Method - must be 8
    uint8_t flg;    // Flags - must be 4 (FEXTRA)
    uint32_t mtime; // Modification Time
    uint8_t xfl;    // eXtra Flags
    uint8_t os;     // Operating System
    uint16_t xlen;  // Size of extra fields - 6 if contain only BGZF (may be more)
    uint8_t si1;    // bsize field id - must be 66  (0x42)
    uint8_t si2;    // bsize field id - must be 67  (0x43)
    uint16_t slen;  // bsize field length - must be 2
    uint16_t bsize; // bsize field field - (compressed block size -1)
} BgzfHeader;

// see: https://docs.google.com/document/d/11yeGa1HzXi96D3VTeMwReW2BzG9hmspyKFxx-yRpLPo/edit
typedef struct __attribute__ ((packed)) MgzfHeader { // 29 bytes
    uint8_t id1;    // Gzip id - must be 31  (0x1f)
    uint8_t id2;    // Gzip id - must be 139 (0x8b)
    uint8_t cm;     // Compression Method - must be 8
    uint8_t flg;    // Flags - must be 20 (0x14) (FEXTRA | FCOMMENT)
    uint32_t mtime; // Modification Time - must be 0
    uint8_t xfl;    // eXtra Flags - must be 0 
    uint8_t os;     // Operating System - must be 0xff
    uint16_t xlen;  // Size of extra fields - must be 8 
    uint8_t si1;    // bsize field id - must be (0x49)
    uint8_t si2;    // bsize field id - must be (0x47)
    uint16_t slen;  // bsize field extra field length - must be 4
    uint32_t bsize; // bsize field - compressed block size - header + body
    char comment[9];// nul-terminated string in the format "C001R015"
} MgzfHeader;

typedef struct GzipFooter {
    uint32_t crc32; // CRC32 of uncompressed data
    uint32_t isize; // Input (i.e. uncompressed) Size
} GzipFooter;

#define GZIP_FOOTER_LEN ((int)sizeof(GzipFooter))

typedef struct __attribute__ ((packed, aligned(2))) GzipHeader { // 10 bytes
    uint8_t id1;    // Gzip id            - must be 31  (0x1f)
    uint8_t id2;    // Gzip id            - must be 139 (0x8b)
    uint8_t cm;     // Compression Method - must be 8
    uint8_t flg;    // Flags              - must be 0
    uint32_t mtime; // Modification Time  - must be 0
    uint8_t xfl;    // eXtra Flags        - must be 0
    uint8_t os;     // Operating System   - must be 3
} GzipHeader;

static FlagsMgzip bgzf_recompression_levels[1+MAX_FLAG_BGZF] = {
    { .library = BGZF_LIBDEFLATE19, .level = 0 },  // --bgzf=0 : BGZF blocks with no compression     
    { .library = BGZF_IGZIP,        .level = 1 },  // --bgzf=1 : note: this is IGZIP LVL0 
    { .library = BGZF_IGZIP,        .level = 2 },  // --bgzf=2 : note: this is IGZIP LVL1
    { .library = BGZF_LIBDEFLATE19, .level = 1 },  // --bgzf=3 
    { .library = BGZF_LIBDEFLATE19, .level = 7 },  // --bgzf=4 
    { .library = BGZF_LIBDEFLATE19, .level = 9 },  // --bgzf=5
};

static const uint8_t mgzip_header_len[NUM_CODECS] = MGZIP_HEADER_LEN_BY_CODEC;

#define bgzf_no_recompression (FlagsMgzip){ .library = BGZF_NO_LIBRARY, .level = BGZF_NO_BGZF }

rom gzstatus_name (GzStatus st)
{
    return IN_RANGE(st, 0, NUM_GZ_STATUSES) ? (rom[])GZSTATUS_NAMES[st] : "InvalidGzStatus";
}

// possible return values, see libdeflate_result in libdeflate.h
static rom libdeflate_error (int err)
{
    switch (err) {
        case LIBDEFLATE_SUCCESS            : return "SUCCESS";
        case LIBDEFLATE_BAD_DATA           : return "BAD DATA";
        case LIBDEFLATE_SHORT_OUTPUT       : return "SHORT OUTPUT";
        case LIBDEFLATE_INSUFFICIENT_SPACE : return "INSUFFICIENT SPACE";
        default                            : return "Undefined libdeflate error";
    }
}

typedef struct { char s[100]; } BgzfBlockStr;
static BgzfBlockStr display_bb (GzBlockZip *bb)
{
    BgzfBlockStr s;
    snprintf (s.s, sizeof (s.s), "{txt_index=%u txt_size=%u compressed_index=%u comp_size=%u is_uncompressed=%u}",
             bb->txt_index, bb->txt_size, bb->compressed_index, bb->comp_size, bb->is_uncompressed);
    return s;
}

static void *bgzf_alloc (void *vb_, unsigned items, unsigned size, FUNCLINE)
{
    return codec_alloc_do ((VBlockP )vb_, (uint64_t)items * (uint64_t)size, 1, func, code_line); // all bzlib buffers are constant in size between subsequent compressions
}

//--------------------------------------------------------------------
// ZIP SIDE - library/level discovery
//--------------------------------------------------------------------

#define BGZF_DISCOVERY_MAX_TESTS 10       // maximum number of BGZF blocks to be tested

void bgzf_initialize_discovery (FileP file)
{
    ASSERTNOTINUSE (file->bgzf_plausible_levels);
    
    // note: tested example files of MGZF, MGSP and IL1M and they don't match any of these libraries.
    if (IS_BGZF(file->effective_codec)) {
        ARRAY_alloc (FlagsMgzip, ll, (LIBDEFLATE_MAX_LEVEL+1)+LIBDEFLATE_MAX_LEVEL+ZLIB_MAX_LEVEL, 
                    false, file->bgzf_plausible_levels, evb, "txt_file->bgzf_plausible_levels");

        int next=0;
        for (int l=0; l <= LIBDEFLATE_MAX_LEVEL; l++) // level=0 only here, bc it would be the same in all libraries
            ll[next++] = (FlagsMgzip){ .library = BGZF_LIBDEFLATE19, .level = l};

        for (int l=1; l <= LIBDEFLATE_MAX_LEVEL; l++)
            ll[next++] = (FlagsMgzip){ .library = BGZF_LIBDEFLATE7,  .level = l};

        for (int l=1; l <= ZLIB_MAX_LEVEL; l++)
            ll[next++] = (FlagsMgzip){ .library = BGZF_ZLIB,         .level = l};
    }

    else if (IS_MGZIP(file->effective_codec)) {
        // bug 1101: we don't yet know the plausible levels for other MGZIP codecs
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
}

// ZIP main thread
void bgzf_finalize_discovery (void)
{    
    int n_levels = txt_file->bgzf_plausible_levels.len32;

    // case: there is no library/level combination for which we can decompress with bgzf=exact
    if (n_levels == 0) {
        bgzf_discover_finalize_testing (0, BGZF_COMP_LEVEL_UNKNOWN); // has BGZF, but cannot identify level

        if (flag.show_bgzf) 
            iprintf ("Discover:%s: is a %s file, generated by an unidentified library\n", txt_name, codec_name (txt_file->effective_codec));
    }

    // case: one or more library/level combinations was verified with all test bgzf blocks (10 blocks, unless file is shorter) 
    else {
        bgzf_discover_finalize_testing (B1ST(FlagsMgzip, txt_file->bgzf_plausible_levels)->library, B1ST(FlagsMgzip, txt_file->bgzf_plausible_levels)->level);

        if (flag.show_bgzf) 
            iprintf ("Discover: %s is a %s file, %s %s level %u\n", txt_name, codec_name (txt_file->effective_codec),
                     (n_levels == 1) ? "identified as generated with" : "with multiple plausible levels, arbitrarily selecting",
                     bgzf_library_name (txt_file->mgzip_flags.library, true), txt_file->mgzip_flags.level);
    }    
}

// ZIP: test a BGZF block against all the remaining plausible levels, and eliminate those that don't match. 
static void bgzf_discover_library_and_level (VBlockP vb, int test_block_i, STRp(comp), STRp(uncomp))
{
    uint32_t header_len = mgzip_header_len[txt_file->effective_codec];
    ASSERT (header_len, "%s_HEADER_LEN missing in mgzip_header_len", codec_name (txt_file->effective_codec));

    if (comp_len <= header_len + GZIP_FOOTER_LEN) {
        txt_file->bgzf_plausible_levels.len = 0;

        if (flag.show_bgzf || flag.show_gz) 
            iprintf ("%s: Block too small - could not identify compression library and level\n", txt_name);

        if (flag.show_gz) exit_ok;

        return;
    }

    // ignore the header and footer of the block
    comp     += header_len;
    comp_len -= header_len + GZIP_FOOTER_LEN;

    // compress with each of the remaining plausible levels - testing if the compression is identical to the actual
    uint32_t recomp_size = uncomp_len * 1.1 + 64 KB; // guessing the max compressed size in the worst case scenario of very bad compression
    char *recomp = MALLOC (recomp_size);
    uint32_t recomp_len;

    for_buf (FlagsMgzip, ll, txt_file->bgzf_plausible_levels) {

        // for large test blocks, skip high compression levels which are not common anyway, as testing is too slow
        if (comp_len > 100 KB && ll->level >= 8)
            continue;

        switch (ll->library) {
            case BGZF_LIBDEFLATE19 : {
                void *compressor = libdeflate_alloc_compressor (vb, ll->level, __FUNCLINE);
                recomp_len = (uint32_t)libdeflate_deflate_compress (compressor, STRa(uncomp), recomp, recomp_size);

                libdeflate_free_compressor (compressor, __FUNCLINE);
                break;
            }

            case BGZF_LIBDEFLATE7 : {
                void *compressor = libdeflate_alloc_compressor_1_7 (ll->level, vb);
                recomp_len = (uint32_t)libdeflate_deflate_compress_1_7 (compressor, STRa(uncomp), recomp, recomp_size);

                libdeflate_free_compressor_1_7 (compressor);
                break;
            }

            case BGZF_ZLIB : {
                z_stream strm = { .zalloc = bgzf_alloc, .zfree  = codec_free_do, .opaque = vb };
                // deflateInit2 with the default zlib parameters, which is also the same as htslib does
                ASSERT0 (deflateInit2 (&strm, ll->level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) == Z_OK, "deflateInit2 failed");

                strm.next_in   = (uint8_t *)uncomp;
                strm.avail_in  = uncomp_len;
                strm.next_out  = (uint8_t *)recomp;
                strm.avail_out = recomp_size;
                ASSERT (deflate (&strm, Z_FINISH) == Z_STREAM_END, "zlib deflate failed: msg=%s", strm.msg);

                recomp_len = recomp_size - strm.avail_out;
                
                ASSERT0 (deflateEnd (&strm) == Z_OK, "deflateEnd failed");
                break;
            }

            default: ABORT ("Invalid library=%d", ll->library);
        }

        bool plausible = str_issame (comp, recomp);

        if (flag.show_bgzf) 
            iprintf ("Discover[%d]: library %s level %u: size_in_file=%u size_in_test=%u plausible=%s\n", 
                     test_block_i, bgzf_library_name (ll->library, true), ll->level, comp_len, recomp_len, YN(plausible));

        if (!plausible) {
            buf_remove (txt_file->bgzf_plausible_levels, FlagsMgzip, BNUM(txt_file->bgzf_plausible_levels, ll), 1);
            ll--; fb_after--; // hack for_buf loop
        }
    }

    FREE (recomp);
}

//--------------------------------------------------------------------
// ZIP SIDE - decompress MGZIP-compressed file and prepare BGZF section
//--------------------------------------------------------------------

uint32_t mgzip_get_max_block_size (void)
{
    uint32_t max_block_size = (uint32_t[NUM_CODECS])MAX_ISIZE_BY_CODEC[txt_file->effective_codec];
    if (!max_block_size) max_block_size = 1;

    return max_block_size;
}

void inc_disk_gz_uncomp_or_trunc_(FileP file, uint64_t inc, FUNCLINE)
{
    __atomic_add_fetch (&file->disk_gz_uncomp_or_trunc, inc, __ATOMIC_RELAXED);

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

static GzStatus mgzip_block_verify_header (FileP file, bool discovering, int header_len, STR8p(prefix))
{
    FILE *fp = (FILE *)file->file;

    // no data at all
    if (file->gz_data.len32 == 0 && feof (fp)) {
        ASSERT0 (!discovering, "unexpected end of file when discovering");
        return GZ_EOF_WITHOUT_EOF_BLOCK; // end of file
    }
    
    // truncated mid-way through header
    ARRAY(uint8_t, h, file->gz_data);
    if (h_len < header_len) {
        if (discovering)
            return GZ_NOT_GZIP; // file smaller than a gzip header - its not GZIP
        
        else if (flag.truncate && feof (fp)) 
            return GZ_TRUNCATED;   // truncated file      
        
        else
            ABORT ("%s file %s appears truncated - it ends with a partial gzip block header. offset=%"PRIu64". If you expect this file to be truncated, use --truncate", 
                   codec_name (file->effective_codec), file->basename, (uint64_t)ftello64 (fp) - h_len); // less than the minimal gz block header size
    }

    // case: this is not a GZIP block at all (see: https://tools.ietf.org/html/rfc1952)
    else if (h[0] != 31 || h[1] != 139 || h[2] != 8) {
        if (discovering)
            return GZ_NOT_GZIP;
        else
            ABORT ("expecting %s file %s to be compressed with gzip format, but it is not. offset=%"PRIu64, 
                   codec_name (file->effective_codec), file->basename, (uint64_t)ftello64 (fp) - h_len);
    }

    // case: this is GZIP block (by the magic) but it is NOT a valid BGZF block (see: https://samtools.github.io/hts-specs/SAMv1.pdf)
    else if (!is_header_prefix (h, prefix, prefix_len)) {
        if (discovering)
            return GZ_IS_OTHER_FORMAT;
        else 
            ABORT ("Encountered a GZIP block that unexpectedly is not %s in %s offset=%"PRIu64" found=%s expected=%s. Solution: use --no-bgzf", 
                   codec_name (file->effective_codec), file->basename, (uint64_t)ftello64 (fp) - h_len, display_gz_header (h, h_len, false).s, display_gz_header (STRa(prefix), false).s);
    }

    return GZ_SUCCESS;
}

static void mgzip_update_file_isizes (FileP file)
{
    // changes since 15.0.63
    // - previously, we didn't add EOF blocks to mgzip_isizes and instead set txt_file->mgzip_flags.has_eof_block, incorrectly assuming that an EOF block can only occur at the end of the file
    // - previously, mgzip_isizes was an array of uint16_t whose elements were (file->gz_data.uncomp_len-1)

    // add isize to buffer that will be written to SEC_MGZIP
    if (!file->mgzip_isizes.len32) { 
        uint64_t est_n_blocks = (file->gz_data.comp_len > 10000) ? ((double)file->disk_size / (double)file->gz_data.comp_len * 1.1) : 0; // >10000 to avoid over-allocating due to a randomly small block
        buf_alloc (evb, &file->mgzip_isizes, 0, MAX_(10000, est_n_blocks), uint32_t, 0, "txt_file->mgzip_isizes"); 
        buf_alloc (evb, &file->mgzip_starts, 0, MAX_(10000, est_n_blocks), uint64_t, 0, "txt_file->mgzip_starts");
    }

    buf_append_one (file->mgzip_isizes, file->gz_data.uncomp_len);     
    buf_append_one (file->mgzip_starts, file->disk_so_far - file->gz_data.len); 
}

static void mgzip_show_truncated (FileP file, uint32_t comp_len_truncated)
{
    iprintf ("TRUNCATED  %s thread=%s comp_len_truncated=%u truncated incomplete final %s block\n", 
             codec_name (file->effective_codec), threads_am_i_main_thread() ? "MAIN" : "COMPUTE", comp_len_truncated, codec_name (file->effective_codec));
} 

bool il1m_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool *is_end_of_vb)
{
    return proposed_isize == 1 MB || (is_eof && proposed_isize < 1 MB);
}

bool emfl_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool *is_end_of_vb)
{
    return !file->max_mgzip_isize || // not set yet because this is the first block  
           proposed_isize == file->max_mgzip_isize || (is_eof && proposed_isize < file->max_mgzip_isize);
}

bool emvl_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool *is_end_of_vb)
{
    return proposed_isize < 512 MB; // sanity check
}

bool mgsp_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool *is_end_of_vb)
{
    if (proposed_isize > 64 MB) return false; // sanity

    // case: first gz block in this VB
    if (!file->num_mgsp_blocks_in_vb) 
        file->mgsp_vb_isize = proposed_isize;

    if (proposed_isize == file->mgsp_vb_isize || // gz block with same isize as first
        proposed_isize == 0                   || // EOF block is always accepted
        proposed_isize - file->mgsp_vb_isize <= file->num_mgsp_blocks_in_vb) {  // last block allowed to be the remainder of dividing the total isizes by the number of blocks: i.e. an integer from zero to the number of block minus one.

        file->num_mgsp_blocks_in_vb++;
        file->max_mgsp_blocks_in_vb = MAX_(file->max_mgsp_blocks_in_vb, file->num_mgsp_blocks_in_vb);
        return true;
    }

    // proposed block doesn't belong to this VB
    else {
        file->num_mgsp_blocks_in_vb = file->mgsp_vb_isize = 0; // initialize for next VB
        *is_end_of_vb = true;
        return false;
    }
}

// read a gz block of a codec that does not contain bsize in the gz header
// returns: discovering: GZ_SUCCESS, GZ_IS_OTHER_FORMAT
//          otherwise:   GZ_SUCCESS
GzStatus mgzip_read_block_no_bsize (FileP file, bool discovering, Codec codec)
{
    START_TIMER;
    static NoBsizeCodecParams no_bsize_codec_params[NUM_CODECS] = (NoBsizeCodecParams[])NO_BSIZE_CODECS_PARAMS;

    NoBsizeCodecParams params = no_bsize_codec_params[codec];
    if (!params.gz_hdr) params.gz_hdr = file->gz_header;

    ASSERT (params.max_bsize, "NO_BSIZE_CODECS_PARAMS missing data for %s", codec_name (codec));

    if (discovering && params.valid_3_blocks_isize)
        params.max_bsize *= 3; // read 3 blocks - txtfile_discover_specific_gz will use this data to verify they all have the same isize

    FILE *fp = (FILE *)file->file;
    file->gz_data.comp_len = file->gz_data.uncomp_len = 0; // init

    // top up gz_data to max_comp_size (or less if EOF)
    txtfile_fread (file, fp, NULL, (int32_t)params.max_bsize - (int32_t)file->gz_data.len32, &file->disk_so_far);
    
    GzStatus status = mgzip_block_verify_header (file, discovering, params.gz_hdr_len, STRa(params.gz_hdr));
    if (status == GZ_EOF_WITHOUT_EOF_BLOCK) {
        if (!discovering) file->no_more_blocks = true;
        return GZ_SUCCESS;
    }

    if (status != GZ_SUCCESS) return status;
        
    // search for block size by beginning of next block. note: we do this even if EOF,
    // because gz_data might contain several gz blocks. note: also NULL if data is too short.
    uint8_t *next_blk = memmem (B8(file->gz_data, params.gz_hdr_len), file->gz_data.len32 - params.gz_hdr_len, params.gz_hdr, params.gz_hdr_len);
    bool end_of_vb = false;

    // case: a block was found, and it is not the last block
    if (next_blk && params.is_valid_isize (file, GET_UINT32 (next_blk - 4), false/* there IS a next block so not EOF*/, &end_of_vb)) {
        
        file->gz_data.uncomp_len = GET_UINT32 (next_blk - 4);
        file->gz_data.comp_len = BNUM (file->gz_data, next_blk);
    }

    // case: remaining data could be a final gz block, or could be a truncated block, 
    // we will know for sure when trying to uncompress it
    else if (feof (fp) && file->gz_data.len32 >= params.gz_hdr_len + GZIP_FOOTER_LEN &&
             params.is_valid_isize (file, GET_UINT32 (BAFT8 (file->gz_data) - 4), true, &end_of_vb)) {
        
        file->gz_data.uncomp_len = GET_UINT32 (BAFT8(file->gz_data) - 4);
        file->gz_data.comp_len = file->gz_data.len32;
        if (!discovering) file->no_more_blocks = true;
    }
    
    // case: end of group of gz blocks that are designated for the current VB
    else if (end_of_vb) {
        // gz_data: keep data comp_len, uncomp_len remain 0, keep len remains >0 : this means end-of-vb
    }

    // case: data in gz_data is does not contain a gz block - either not the right codec file or is truncated
    else {
        if (discovering)
            return GZ_IS_OTHER_FORMAT;

        // data is not IL1M somewhere in the middle of the file...
        ASSERT (feof (fp), "Encountered a GZIP block that unexpectedly is not %s in %s offset=%"PRIu64"\nSolution: use --no-bgzf", 
                codec_name (codec), file->basename, (uint64_t)ftello64 ((FILE *)file->file) - file->gz_data.len);

        // case: final data in file is not a full gz block and truncation allowed: 
        // account and then ignore the data that will not be gz-decompressed 
        if (flag.truncate) {
            WARN ("FYI: %s is truncated - its final %s block in incomplete. Dropping final %u bytes of the GZ data.", 
                  txt_name, codec_name (codec), file->gz_data.len32);

            if (flag.show_bgzf) mgzip_show_truncated (file, file->gz_data.len32);

            inc_disk_gz_uncomp_or_trunc (file, file->gz_data.len);
            file->gz_data.len32 = file->gz_data.uncomp_len = 0;
            segconf.zip_txt_modified = true;
            file->no_more_blocks = true;
        }

        else 
            ABORTINP ("%s is truncated mid-way through %s block. Tip: If this is expected, use --truncate to discard the final partial %s block", 
                      txt_name, codec_name (codec), codec_name (codec));
    }

    if (file->gz_data.comp_len && !discovering) 
        mgzip_update_file_isizes (file);

    COPY_TIMER_EVB (mgzip_read_block_no_bsize);
    return GZ_SUCCESS;
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
            feof (fp) ? "Tip: If the file is expected to be truncated, you use --truncate to disregard the final partial BGZF block." : "");
}

static void bgzf_mgzf_verify_eof_block (FileP file, STRp(eof_block))
{
    // case: valid EOF block
    if (str_issame_(B1STc(file->gz_data), file->gz_data.comp_len, eof_block, eof_block_len))
        file->num_EOF_blocks++;
 
    // case: an isize=0 block that is not our EOF block (gz header differs or minimal payload differs - its possible in gzip)
    else {
        file->non_EOF_zero_block_found = true; // we won't be able to reconstruct exactly, as PIZ reconstructs the EOF block if isize=0
    
        if (flag.show_bgzf)
            iprintf ("DETECTED non-EOF zero block - we can't reconstruct this file --exact-ly: %s\n", 
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

// ZIP: reads and validates a BGZF block
// returns: discoverying: GZ_SUCCESS, GZ_NOT_GZIP, GZ_IS_OTHER_FORMAT
//          otherwise:    GZ_SUCCESS, GZ_EOF_WITHOUT_EOF_BLOCK, GZ_TRUNCATED
static GzStatus mgzf_read_block_do (FileP file, // txt_file is not yet assigned when called from txtfile_discover_specific_gz
                                    bool discovering)
{
    #define MGZF_CHUCK_SIZE ((uint32_t)(16 MB)) // max amount we read from disk at a time  
    
    FILE *fp = (FILE *)file->file;
    file->gz_data.comp_len = file->gz_data.uncomp_len = 0; // init
    uint32_t bsize = 0;

    // top-up if needed - in rare cases - twice (this happens very large block where comp_size>MGZF_CHUCK_SIZE and not known yet - first read of MGZF_CHUCK_SIZE includes the header)
    for (int i=0; i < 2; i++) 
        if ((!mgzf_get_bsize (file, &bsize) || bsize > file->gz_data.len32) && !feof (fp)) {
            int32_t chunk_size = MAX_(MGZF_CHUCK_SIZE, bsize);

            txtfile_fread (file, fp, NULL, chunk_size - (int32_t)file->gz_data.len32, &file->disk_so_far);
        }

    GzStatus status = mgzip_block_verify_header (file, discovering, MGZF_HEADER_LEN, _8(MGZF_PREFIX));
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

        ASSERT (h->bsize == MGZF_EOF_LEN || is_valid_comment, "Invalid MGZF comment: gz_header={ %s }. Solution: run with --no-bgzf", display_gz_header ((bytes)STRb(file->gz_data), false).s);

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
        int32_t chunk_size = flag.zip_uncompress_source_during_read ? 150 KB // a bit more than the default block-device read-ahead buffer (128KB) for best parallelization between disk read-ahead and CPU decompression
                                                                    : BGZF_MAX_CHUCK_SIZE; // bigger block is faster if we are prepared to yield the CPU when waiting for the disk 
        txtfile_fread (file, fp, NULL, chunk_size - (int32_t)file->gz_data.len32, &file->disk_so_far);
    }

    uint32_t bsize = (uint32_t)LTEN16 (B1ST(BgzfHeader, file->gz_data)->bsize) + 1;

    GzStatus status = mgzip_block_verify_header (file, discovering, BGZF_HEADER_LEN, _8(BGZF_PREFIX));
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
    if (file->gz_data.comp_len) return GZ_SUCCESS; // we already have 1 block

    GzStatus ret = (codec == CODEC_BGZF) ? bgzf_read_block_do (file, discovering)
                                         : mgzf_read_block_do (file, discovering);
    switch (ret) {
        case GZ_SUCCESS: // successful read of a BGZF block            
            mgzip_update_file_isizes (file);
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
                WARN ("FYI: %s is truncated - its final BGZF block in incomplete. Dropping final %u bytes of the GZ data.", txt_name, file->gz_data.len32);

                if (flag.show_bgzf) mgzip_show_truncated (file, file->gz_data.len32);

                inc_disk_gz_uncomp_or_trunc (file, file->gz_data.len);
                file->gz_data.len32 = file->gz_data.comp_len = file->gz_data.uncomp_len = 0; // discard partial BGZF block
                segconf.zip_txt_modified = true;
                
                ret = GZ_SUCCESS;
                break;
            }   
            
            else
                ABORTINP ("%s is truncated mid-way through BGZF block. Tip: If this is expected, use --truncate to discard the final partial BGZF block", txt_name);

        default:
            ABORT ("Unexpected ret=%s", gzstatus_name (ret));
    }

    if (!discovering)
        file->no_more_blocks = (file->gz_data.comp_len == file->gz_data.len32 && feof ((FILE*)file->file));

    COPY_TIMER_EVB (mgzip_read_block_with_bsize);
    return ret;
}

// ZIP: BGZF section per txt_file component
void mgzip_compress_mgzip_section (void)
{
    // cases where we don't write the BGZF blocks section
    if (!txt_file->mgzip_isizes.len ||        // this txt file is not compressed with MGZIP - we don't need a MGZIP section
        txt_file->mgzip_flags.level == BGZF_COMP_LEVEL_UNKNOWN ||  // we don't know the level - so PIZ will reconstruct at default level
        txt_file->non_EOF_zero_block_found || // we don't know how to reconstructed-exactly a zero block other than an EOF block, so no point storing isizes
        segconf.zip_txt_modified)             // the file has changed and we can't reconstruct to the same blocks
        return;
    
    // sanity check
    uint64_t total_isize=0;
    for_buf (uint32_t, isize_p, txt_file->mgzip_isizes)
        total_isize += *isize_p;
    
    ASSERT (total_isize == txt_file->txt_data_so_far_single + txt_file->last_truncated_line_len, 
            "Expecting total_isize=%"PRId64" == txt_data_so_far_single=%"PRId64,
            total_isize, txt_file->txt_data_so_far_single);

    BGEN_u32_buf (&txt_file->mgzip_isizes, NULL);
    txt_file->mgzip_isizes.len *= sizeof (uint32_t);

    Codec codec = codec_assign_best_codec (evb, NULL, &txt_file->mgzip_isizes, SEC_MGZIP);

    evb->comp_i = flag.zip_comp_i; // this goes into SectionEntFileFormat.comp_i via sections_add_to_list
    zfile_compress_section_data_ex (evb, NULL, SEC_MGZIP, &txt_file->mgzip_isizes, NULL, 0, codec, (SectionFlags)txt_file->mgzip_flags, NULL);
    txt_file->mgzip_isizes.len /= sizeof (uint32_t); // restore

    z_file->comp_num_EOF_blocks[flag.zip_comp_i] = txt_file->num_EOF_blocks;
}

// uncompresses a BGZF block in vb->comp_txt_data referred to by bb, into its place in vb->txt_data as prescribed by bb
// might be called from main thread or compute threads
void mgzip_uncompress_one_block (VBlockP vb, GzBlockZip *bb, Codec codec)
{
    if (bb->is_uncompressed) return; // already decompressed (or an empty (e.g. EOF) block) - nothing to do

    ASSERT0 (vb->gzip_compressor, "vb->gzip_compressor=NULL");

    int header_len = mgzip_header_len[codec];
    ASSERT (header_len, "%s_HEADER_LEN missing in mgzip_header_len", codec_name (codec));

    uint8_t *h = B8(vb->comp_txt_data, bb->compressed_index);

    // verify that entire block is within vb->comp_txt_data
    ASSERT (bb->compressed_index + header_len < vb->comp_txt_data.len && // we have at least the header - we can access bsize
            bb->compressed_index + bb->comp_size <= vb->comp_txt_data.len, 
            "%s: %s block size goes past the end of in vb->comp_txt_data: bb=%s compressed_index=%u vb->comp_txt_data.len=%"PRIu64, 
            VB_NAME, codec_name (codec), display_bb (bb).s, bb->compressed_index, vb->comp_txt_data.len);

    ASSERT (h[0]==31 && h[1]==139 && h[2]==8, "%s: invalid %s block in vb->comp_txt_data: compressed_index=%u", VB_NAME, codec_name (codec), bb->compressed_index);

    // possibly grow txt_data: can happen data is MGZIP and its length exceeds vb_size due to last block going over
    buf_alloc (vb, &vb->txt_data, 0/*don't use "more" bc Ltxt already incremented*/, bb->txt_index + bb->txt_size + TXTFILE_READ_VB_PADDING, char, 0, "txt_data");

    enum libdeflate_result ret = 
        libdeflate_deflate_decompress (vb->gzip_compressor, 
                                       h + header_len, bb->comp_size - header_len - GZIP_FOOTER_LEN, // compressed
                                       Btxt (bb->txt_index), bb->txt_size, NULL);  // uncompressed

    // account for the case of decompression, and also the case bb is discarded due to a certain truncate situation (see below).
    inc_disk_gz_uncomp_or_trunc (txt_file, bb->comp_size);

    // case: final IL1M block, which is truncated, but we have --truncate, and the garbage last word
    // unluckily < 1MB so it went undetected as a legimiate block in il1m_is_valid_isize. we drop this block now.
    if (ret != LIBDEFLATE_SUCCESS && TXT_IS_IL1M && bb->is_eof) {
        if (flag.truncate) {
            // receive updates made by main thread to mgzip_isizes, mgzip_starts: no more are going to happen as we reached eof
            __atomic_thread_fence (__ATOMIC_ACQ_REL); 

            txt_file->mgzip_isizes.len--; // remove truncated block from isizes
            txt_file->mgzip_starts.len--;
            mgzip_show_truncated (txt_file, bb->comp_size);
            return; // with bb->is_uncompressed=false
        }

        else {
            ABORT ("Failed to decompress the final %s block of the file: %s. Tip: If it is expected that the file is truncated, use --truncate to ignore the defective final block.", 
                   codec_name (vb->txt_codec), libdeflate_error(ret));
        }
    }

    ASSERT (ret == LIBDEFLATE_SUCCESS, "libdeflate_deflate_decompress failed: %s", libdeflate_error(ret));

    bb->is_uncompressed = true;

    if (flag.show_bgzf)
        iprintf ("UNCOMPRESS %s thread=%s%s i=%u comp_index=%u comp_len=%u txt_index=%u txt_len=%u eof=%s%s%s\n",
                 codec_name (codec), threads_am_i_main_thread() ? "MAIN" : "COMPUTE", 
                 cond_str (vb->vblock_i, " vb=", VB_NAME),
                 BNUM (vb->gz_blocks, bb), bb->compressed_index, bb->comp_size, bb->txt_index, bb->txt_size, TF(bb->is_eof),
                 cond_str (bb->txt_size, " uncomp[5]=", str_to_printable_(Btxt(bb->txt_index), MIN_(5, Ltxt - bb->txt_index)).s),
                 bb->comp_size == BGZF_EOF_LEN ? " EOF" : "");

    // discover which gzip library and compression level were used (testing the first few BGZF blocks)
    int test_block_i=0;
    if (txt_file->bgzf_plausible_levels.len32 && // only >0 in BGZF
        txt_file->bgzf_plausible_levels.count < BGZF_DISCOVERY_MAX_TESTS && // fail fast without atomic
        (test_block_i = __atomic_fetch_add (&txt_file->bgzf_plausible_levels.count, 1, __ATOMIC_RELAXED)) < BGZF_DISCOVERY_MAX_TESTS) { // note: if multiple threads test in parallel, count might be incremented beyond BGZF_DISCOVERY_MAX_TESTS - that's ok
        
        bgzf_discover_library_and_level (vb, test_block_i, (rom)h, bb->comp_size, Btxt (bb->txt_index), bb->txt_size);        
    }

    // in case of --show_gz: report and exit (otherwise, we finalize in the main thread to avoid thread issues with updating txt_file->mgzip_flags)
    if (flag.show_gz && (test_block_i == BGZF_DISCOVERY_MAX_TESTS-1 || !txt_file->bgzf_plausible_levels.len32)) 
        bgzf_finalize_discovery();
}

// ZIP: called from the compute thread: zip_compress_one_vb and main thread: txtfile_read_block_mgzip
void mgzip_uncompress_vb (VBlockP vb, Codec codec)
{
    START_TIMER;
    ASSERTNOTEMPTY (vb->gz_blocks);

    vb->gzip_compressor = libdeflate_alloc_decompressor(vb, __FUNCLINE);

    uint32_t total_vb_isizes = 0;
    for_buf (GzBlockZip, bb, vb->gz_blocks) {
        mgzip_uncompress_one_block (vb, bb, codec);
        total_vb_isizes += bb->txt_size;
    }

    // sanity - total_vb_isizes is at least the size of the vb (can be more, bc the first/last bb can be shared with previous/next vb)
    ASSERT (total_vb_isizes >= Ltxt, "%s: Expecting total_vb_isizes=%u >= Ltxt=%u. codec=%s", 
            VB_NAME, total_vb_isizes, Ltxt, codec_name (txt_file->effective_codec));

    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor, __FUNCLINE);

    buf_free (vb->comp_txt_data); // now that we are finished decompressing we can free it

    if (flag.show_time) {
        if (threads_am_i_main_thread ()) COPY_TIMER (bgzf_io_thread)
        else                             COPY_TIMER (bgzf_compute_thread);
    }

    COPY_TIMER (mgzip_uncompress_vb);
}

// ZIP: decompresses a prescribed BGZF block when re-reading DEPN lines
static inline void bgzf_uncompress_one_prescribed_block (VBlockP vb, STRp(bgzf_block), STRc(uncomp_block), uint64_t bb_i)
{
    START_TIMER;

    BgzfHeader *h = (BgzfHeader *)bgzf_block;

    if (flag.show_bgzf)
        iprintf ("REREAD  %s reread bb_i=%"PRIu64" comp_size=%u uncomp_size=%u ", 
                 VB_NAME, bb_i, bgzf_block_len, uncomp_block_len);

    enum libdeflate_result ret = 
        libdeflate_deflate_decompress (vb->gzip_compressor, 
                                       h+1, bgzf_block_len - BGZF_HEADER_LEN - GZIP_FOOTER_LEN, // compressed
                                       STRa(uncomp_block), NULL);  // uncompressed

    ASSERT (ret == LIBDEFLATE_SUCCESS, "%s: libdeflate_deflate_decompress failed: %s. bgzf_block_len=%u uncomp_block_len=%u bb_i=%"PRIu64, 
            VB_NAME, libdeflate_error(ret), bgzf_block_len, uncomp_block_len, bb_i);

    if (flag.show_bgzf)
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
    ASSERT (header_bytes == BGZF_HEADER_LEN && !memcmp (bgzf_block, BGZF_PREFIX, STRLEN(BGZF_PREFIX)),
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

    vb->gzip_compressor = libdeflate_alloc_decompressor(vb, __FUNCLINE);

    for_buf (RereadLine, line, vb->reread_prescription) {
        
        // a line might span 1 or more BGZF blocks
        while (line->line_len) { 
            ASSERT (line->offset.bb_i < txt_file->mgzip_starts.len32, "Expecting bb_i=%"PRIu64" < mgzip_starts.len=%"PRIu64, 
                    (uint64_t)line->offset.bb_i, txt_file->mgzip_starts.len);

            uint64_t offset = *B64 (txt_file->mgzip_starts, line->offset.bb_i);
            uint32_t isize  = *B32 (txt_file->mgzip_isizes, line->offset.bb_i);

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

    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor, __FUNCLINE);
}

void bgzf_libdeflate_1_7_initialize (void)
{
    libdeflate_set_memory_allocator_1_7 (bgzf_alloc, codec_free_do);
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
            compressed_size += bb[i].comp_size;

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
    vb->vb_mgzip_i = txt_file->mgzip_isizes.len; // index of first bgzf block to be used by the VB

    if (!txt_file->unconsumed_mgzip_blocks.len) return; // happens when either unconsumed_bytes=0 or not a BGZF-compressed file

    // data in the first BGZF block already consumed by previous VB or txt_header 
    vb->gz_blocks.consumed_by_prev_vb = vb->line_bgzf_uoffset = txt_file->unconsumed_mgzip_blocks.consumed_by_prev_vb;

    // copy all unconsumed BGZF blocks - we might not need all of them - the unconsumed ones will moved back in mgzip_copy_unconsumed_blocks
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
// PIZ SIDE - setting up BGZF for a particular txt file
//-----------------------------------------------------

static Buffer isizes = {}; // Will be grabbed into txt_file->mgzip_isizes.

static inline FlagsMgzip recompression_template (int bgzf_level)
{
    return (FlagsMgzip){ .level   = bgzf_recompression_levels[bgzf_level].level, // a 4-bit bitfield 
                         .library = bgzf_recompression_levels[bgzf_level].library };
}

// PIZ, after calling bgzf_load_isizes
static inline bool is_exact (void)
{
    return txt_file->mgzip_isizes.len > 0;
}

static FlagsMgzip bgzf_load_isizes (CompIType comp_i, bool show_only) 
{
    Section sec = sections_get_comp_bgzf_sec (comp_i);
    if (!sec) ignore: {
        ASSERTW0 (show_only, "FYI: --bgzf=exact ignored, because when compressing, genozip could not identify parameters of the .gz file");
        goto fallback; // this component doesn't contain a BGZF section
    }

    int32_t offset = zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", SEC_MGZIP, sec);

    SectionHeaderP header = (SectionHeaderP)Bc(evb->z_data, offset);

    // if we don't know the compression level (in older Genozip versions we wrote the SEC_MGZIP even 
    // if level discovery failed)
    if (header->flags.mgzip.level == BGZF_COMP_LEVEL_UNKNOWN) 
        goto ignore;

    zfile_uncompress_section (evb, header, &isizes, "txt_file->mgzip_isizes", 0, SEC_MGZIP);

    if (show_only) {
        buf_destroy (isizes);
        goto fallback;
    }

    if (VER2(15,63)) {
        isizes.len /= sizeof (uint32_t);
        BGEN_u32_buf (&isizes, NULL);
    }

    // up to 15.0.62 buffer was 16 bit, values were (isize-1), and EOF was indicated by header.has_eof_block
    else {
        isizes.len /= sizeof (uint16_t);
        buf_alloc (evb, &isizes, 0, isizes.len + 1, uint32_t, 0, NULL);

        for_buf_tandem_back (uint16_t, isize16, isizes, uint32_t, isize32, isizes)
            *isize32 = (uint32_t)BGEN16 (*isize16) + 1;

        if (header->flags.mgzip.OLD_has_eof_block)
            BNXT32(isizes) = 0; // append EOF block
    }

    return header->flags.mgzip; // mgzip_isizes successfully loaded

fallback:
    return recompression_template (BGZF_DEFAULT_LEVEL);
}          

// PIZ: called from main thread after reading txt_header's header
FlagsMgzip mgzip_piz_calculate_mgzip_flags (CompIType comp_i, Codec src_codec)
{
    #define C(cdc) (src_codec == CODEC_##cdc)
    FlagsMgzip mgzip_flags;

    #define HAS_EXT(x) filename_has_ext (flag.out_filename, #x)
    bool bgzf_implied_by_out_filename = flag.out_filename && (HAS_EXT(.gz) || HAS_EXT(.bgz) || HAS_EXT(.bam));
    bool no_bgzf_implied_by_out_filename = file_piz_get_dt_of_out_filename() == flag.out_dt && !(HAS_EXT(.gz) || HAS_EXT(.bgz) || HAS_EXT(.bam));
    bool isizes_loaded = false;

    // cases where there is no BGZF re-compression 
    if (flag.test    || 
        OUT_DT(CRAM) ||
        (flag.bgzf == 0 && !OUT_DT(BAM) && !OUT_DT(BCF)) || // note: in BCF and BAM --bgzf=0 means BGZF blocks with no compression (as opposed to no BGZF at all)
        (flag.bgzf == BGZF_BY_ZFILE && C(NONE) && flag.reconstruct_as_src))  // case: --bgzf=exact and source codec was CODEC_NONE
        
        mgzip_flags = bgzf_no_recompression; 
    
    // case: reconstructing BCF: piz sends VCF to bcftools in CODEC_NONE, and bcftools compressed by the level given in mgzip_flags
    else if (OUT_DT(BCF))
        mgzip_flags = (FlagsMgzip){ .library = BGZF_EXTERNAL_LIB, 
                                  .level = (flag.bgzf < 0) ? 4 : (int[]){0, 2, 4, 6, 8, 9 }[flag.bgzf] }; // convert Genozip level 0-5 to bcftools level 0-9
    
    // case: --bgzf=exact and source codec is other than CODEC_NONE 
    else if (flag.bgzf == BGZF_BY_ZFILE && !C(NONE)) {
        mgzip_flags = bgzf_load_isizes (comp_i, false); 
        isizes_loaded = true;
    }
    
    // case: --bgzf=0 to 5
    else if (flag.bgzf >= 0) {
        mgzip_flags = recompression_template (flag.bgzf); // set to --bgzf command line value

        // if user specified --bgzf and --output - make sure output filename is .gz, .bam or .bcf
        ASSINP (flag.force || !flag.out_filename || bgzf_implied_by_out_filename || HAS_EXT(.bcf) || mgzip_flags.level==0, 
                "using %s in combination with %s for outputting a %s file, requires the output filename to end with %s (override with --force)", 
                OT("output", "o"), OT("bgzf", "z"), dt_name(flag.out_dt), OUT_DT(BAM)?".bam" : OUT_DT(BCF)?".bcf" : ".gz");
    
        ASSINP0 (!OUT_DT(BCF) || flag.bgzf != BGZF_BY_ZFILE, "cannot use --bgzf=exact when outputing a BCF file"); // because we have no control over bcftools' BGZF block generation
    }

    // case: genocat to stdout without --bgzf: - no re-compression. 
    else if (is_genocat && !flag.out_filename)
        mgzip_flags = OUT_DT(BAM) ? bgzf_recompression_levels[0] : bgzf_no_recompression; // file_open_txt_write interprets level=0 as CODEC_BGZF without compression for BAM, and CODEC_NONE for other types 

    // case: genocat or genounzip out_filename and no --bgzf: - determine by file name (except BAM - bgzf regardless of filename)
    else if (flag.out_filename)
        mgzip_flags = (bgzf_implied_by_out_filename || OUT_DT(BAM) || (!no_bgzf_implied_by_out_filename && !C(NONE))) ? bgzf_recompression_levels[BGZF_DEFAULT_LEVEL] : bgzf_no_recompression;
    
    // case: genounzip without explicit filename, and no --bgzf: default compression or no compression
    else
        // note: for bz2, xz, and zip - we reconstruct as gz too. better choice than plain.
        mgzip_flags = (IS_GZIP(src_codec) || C(BAM) || C(BZ2) || C(XZ) || C(ZIP)) ? bgzf_recompression_levels[BGZF_DEFAULT_LEVEL] : bgzf_no_recompression; // note: similar logic to txtheader_piz_get_filename
    
    // case: user wants to see this section header, despite not needing BGZF data
    if (!isizes_loaded && (flag.only_headers == SEC_MGZIP+1 || flag.only_headers == SHOW_ALL_HEADERS))
        bgzf_load_isizes (comp_i, true); 
                
    if (flag.show_bgzf)
        iprintf ("comp_i=%u with src_codec=%s out_dt=%s: calculated mgzip_flags={%s, %d}\n",
                 comp_i, codec_name (src_codec), dt_name (flag.out_dt), bgzf_library_name (mgzip_flags.library, true), mgzip_flags.level);

    return mgzip_flags;
    #undef C
}

// PIZ main thread: update txt_file with BGZF info calculated earlier
void bgzf_piz_set_txt_file_bgzf_info (FlagsMgzip mgzip_flags, bytes codec_info)
{
    memcpy (txt_file->bgzf_signature, codec_info, 3);
    
    if (isizes.len) 
        buf_grab (evb, txt_file->mgzip_isizes, "txt_file->mgzip_isizes", isizes);
        
    txt_file->mgzip_flags = mgzip_flags;

    // sanity        
    ASSERT (txt_file->mgzip_flags.level >= 0 && txt_file->mgzip_flags.level <= BGZF_MAX_LEVEL, "txt_file->mgzip_flags.level=%u ∉ [0,%u]", 
            txt_file->mgzip_flags.level, BGZF_MAX_LEVEL);

    ASSERT (txt_file->mgzip_flags.library >= 0 && txt_file->mgzip_flags.library < NUM_BGZF_LIBRARIES, "txt_file->mgzip_flags.library=%u ∉ [0,%u]", 
            txt_file->mgzip_flags.level, NUM_BGZF_LIBRARIES-1);
}

//-----------------------------------------------------
// PIZ SIDE - compressing txt_file with BGZF
//-----------------------------------------------------

static void bgzf_alloc_compressor (VBlockP vb, FlagsMgzip mgzip_flags)
{
    ASSERT0 (!vb->gzip_compressor, "expecting vb->gzip_compressor=NULL");

    switch (mgzip_flags.library) {
        case BGZF_LIBDEFLATE19:
            vb->gzip_compressor = libdeflate_alloc_compressor (vb, mgzip_flags.level, __FUNCLINE);
            break;

        case BGZF_LIBDEFLATE7:
            vb->gzip_compressor = libdeflate_alloc_compressor_1_7 (mgzip_flags.level, vb);
            break;

        case BGZF_ZLIB:
            vb->gzip_compressor = bgzf_alloc (vb, 1, sizeof (z_stream), __FUNCLINE);
            *(z_stream *)vb->gzip_compressor = (z_stream){ .zalloc = bgzf_alloc, .zfree = codec_free_do, .opaque = vb };
            break;

        case BGZF_IGZIP:
            ASSERT (mgzip_flags.level==1 || mgzip_flags.level==2, "igzip: expecting mgzip_flags.level=%u ∈[1,2]", mgzip_flags.level);

            vb->gzip_compressor = bgzf_alloc (vb, 1, (int[]){ 1+ISAL_DEF_LVL0_DEFAULT, ISAL_DEF_LVL1_DEFAULT }[mgzip_flags.level-1], __FUNCLINE); // 1+ to avoid 0
            break;

        default:             
            ABORT ("Invalid mgzip_flags.library=%d", mgzip_flags.library);
    }
}

static void bgzf_free_compressor (VBlockP vb, FlagsMgzip mgzip_flags)
{
    switch (mgzip_flags.library) {  
        case BGZF_LIBDEFLATE7  :
            libdeflate_free_compressor_1_7 (vb->gzip_compressor);
            break;

        case BGZF_LIBDEFLATE19 :
            libdeflate_free_compressor (vb->gzip_compressor, __FUNCLINE);
            break;

        case BGZF_IGZIP :
        case BGZF_ZLIB  :
            codec_free (vb, vb->gzip_compressor);
            break;

        default:
            ABORT ("Invalid mgzip_flags.library=%d", mgzip_flags.library);
    }

    vb->gzip_compressor = NULL;
}

static void bgzf_show_compress (VBlockP vb, int32_t block_i, uint32_t comp_index, uint32_t comp_len, rom txt, uint32_t txt_index, uint32_t isize)
{
    iprintf ("COMPRESS thread=%s %s i=%d comp_txt_data.i=%u bsize=%u txt_data.i=%d isize=%u%s%s\n",
             threads_am_i_main_thread() ? "MAIN" : threads_am_i_writer_thread() ? "WRITER" : "COMPUTE", VB_NAME, block_i,
             comp_index, comp_len + BGZF_HEADER_LEN + GZIP_FOOTER_LEN, txt_index, isize,
             cond_str (isize, " uncomp[5]=", str_to_printable_(txt, MIN_(isize, 5)).s),
             comp_len == BGZF_EOF_LEN ? " EOF" : "");
}

static void bgzf_compress_one_block (VBlockP vb, rom in, uint32_t isize, 
                                     int32_t block_i, int32_t txt_index) // for show_bgzf (both may be negative - indicating previous VB)
{
    START_TIMER;

    ASSERT0 (vb->gzip_compressor, "vb->gzip_compressor=NULL");

    #define BGZF_MAX_CDATA_SIZE (BGZF_MAX_BLOCK_SIZE - BGZF_HEADER_LEN - GZIP_FOOTER_LEN)

    buf_alloc (vb, &vb->comp_txt_data, BGZF_MAX_BLOCK_SIZE, 0, char, 1.2, "comp_txt_data");
    uint32_t comp_index = vb->comp_txt_data.len32, comp_len;

    // if this is a isize=0 block, we reconstruct it as an EOF block, not through the compressor that might generate a different isize=0 block
    if (!isize) {
        buf_add (&vb->comp_txt_data, _S(BGZF_EOF));
        comp_len = BGZF_EOF_LEN;
    }

    else { // not EOF block
        BgzfHeader *header = (BgzfHeader *)BAFTc (vb->comp_txt_data);
        buf_add (&vb->comp_txt_data, BGZF_EOF, BGZF_HEADER_LEN); // template of header - only bsize needs updating

        if (txt_file->mgzip_flags.library == BGZF_IGZIP) {
            struct isal_zstream strm;
            isal_deflate_stateless_init (&strm);
            strm.gzip_flag      = ISAL_DEFLATE; 
            strm.flush          = NO_FLUSH;
            strm.level          = txt_file->mgzip_flags.level - 1; // note: level 1,2 in mgzip_flags corrsponds to IGZIP level 0,1
            strm.level_buf_size = (int[]){ ISAL_DEF_LVL0_DEFAULT, ISAL_DEF_LVL1_DEFAULT }[strm.level];
            strm.level_buf      = vb->gzip_compressor;
            strm.next_in        = (uint8_t *)in;
            strm.avail_in       = isize;
            strm.next_out       = BAFT8 (vb->comp_txt_data);
            strm.avail_out      = BGZF_MAX_CDATA_SIZE + GZIP_FOOTER_LEN;
            
            int ret = isal_deflate_stateless (&strm);
            ASSERT (ret == ISAL_DECOMP_OK, "%s: isal_deflate_stateless: %s. isize=%u", VB_NAME, isal_error (ret), isize);

            comp_len = BGZF_MAX_CDATA_SIZE + GZIP_FOOTER_LEN - strm.avail_out; 
        }

        else if (txt_file->mgzip_flags.library == BGZF_LIBDEFLATE19) { // libdeflate 1.19

            comp_len = (int)libdeflate_deflate_compress (vb->gzip_compressor, in, isize, BAFTc (vb->comp_txt_data), BGZF_MAX_CDATA_SIZE);

            // in case the compressed data doesn't fit in one BGZF block, move to compressing at the maximum level. this can
            // happen theoretically (maybe) if the original data was compressed with a higher level, and an uncompressible 64K block was
            // compressed to just under 64K while in our compression level it is just over 64K.
            if (!comp_len) {
                void *high_compressor = libdeflate_alloc_compressor (vb, LIBDEFLATE_MAX_LEVEL, __FUNCLINE); // libdefate's highest level
                comp_len = libdeflate_deflate_compress (high_compressor, in, isize, BAFTc (vb->comp_txt_data), BGZF_MAX_CDATA_SIZE);
                libdeflate_free_compressor (high_compressor, __FUNCLINE);
            }
        }

        else if (txt_file->mgzip_flags.library == BGZF_LIBDEFLATE7) { // libdeflate 1.7

            comp_len = (int)libdeflate_deflate_compress_1_7 (vb->gzip_compressor, in, isize, BAFTc (vb->comp_txt_data), BGZF_MAX_CDATA_SIZE);

            // see comment in BGZF_LIBDEFLATE19 above
            if (!comp_len) {
                void *high_compressor = libdeflate_alloc_compressor_1_7 (LIBDEFLATE_MAX_LEVEL, vb); // libdefate's highest level
                comp_len = libdeflate_deflate_compress_1_7 (high_compressor, in, isize, BAFTc (vb->comp_txt_data), BGZF_MAX_CDATA_SIZE);
                libdeflate_free_compressor_1_7 (high_compressor);
            }
        }

        else { // zlib
            #define strm ((z_stream *)vb->gzip_compressor)

            ASSERT0 (deflateInit2 (vb->gzip_compressor, txt_file->mgzip_flags.level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) == Z_OK, 
                    "deflateInit2 failed");

            strm->next_in   = (uint8_t *)in;
            strm->avail_in  = isize;
            strm->next_out  = BAFT8 (vb->comp_txt_data);
            strm->avail_out = BGZF_MAX_CDATA_SIZE;
            ASSERT (deflate (vb->gzip_compressor, Z_FINISH) == Z_STREAM_END, "deflate failed: msg=%s", strm->msg);

            comp_len = BGZF_MAX_CDATA_SIZE - strm->avail_out;
            
            ASSERT0 (deflateEnd (vb->gzip_compressor) == Z_OK, "deflateEnd failed");
            #undef strm
        }

        ASSERT (comp_len, "cannot compress block with %u bytes into a BGZF block with %u bytes", isize, BGZF_MAX_BLOCK_SIZE);
        vb->comp_txt_data.len32 += comp_len;

        header->bsize = LTEN16 ((uint16_t)(BGZF_HEADER_LEN + comp_len + GZIP_FOOTER_LEN - 1));

        GzipFooter footer = { .crc32 = LTEN32 (crc32 (0, in, isize)),
                              .isize = LTEN32 (isize) };
        buf_add (&vb->comp_txt_data, (rom)&footer, GZIP_FOOTER_LEN);
    }

    if (flag.show_bgzf)
        bgzf_show_compress (vb, block_i, comp_index, comp_len, in, txt_index, isize);

    COPY_TIMER (bgzf_compress_one_block);
} 

// appends file data to wvb->comp_txt_data
void bgzf_write_finalize (void)
{
    // if we attempted to reconstruct the BGZF block to the original file's mgzip_isizes - warn if we were unlucky and failed
    // note: EOF block(s) were already added according to mgzip_isizes
    if (is_exact()) {
        uint8_t signature[3];
        bgzf_sign (txt_file->disk_so_far, signature);
        
        bool verified = !memcmp (signature, txt_file->bgzf_signature, 3);
        
        // verify that we were successful in recompressing --exact-ly
        ASSERTW (verified, "FYI: %s is recompressed with %s (.gz). However, it seems that the original file was compressed with a different compression library than genozip uses, resulting in a slightly different level of compression. Rest assured that the actual data is identical.", 
                 txt_name, codec_name (txt_file->effective_codec));

        if (flag.show_bgzf) {
            #define INT_SIGN(s) ((uint32_t)s[0] | ((uint32_t)s[1] << 8) | ((uint32_t)s[2] << 16))
            if (verified) iprintf ("VERIFY <exact> recompression SUCCEEDED: (file_size=%"PRIu64" %% 16777216) = %u, same as source file\n",   txt_file->disk_so_far, INT_SIGN (signature));
            else          iprintf ("VERIFY <exact> recompression SUCCEEDED: (file_size=%"PRIu64" %% 16777216) = %u but source was file %u\n", txt_file->disk_so_far, INT_SIGN (signature), INT_SIGN (txt_file->bgzf_signature));
        }
    }

    // add EOF block when not reconstructing --exact-ly
    else {
        if (flag.show_bgzf) 
            bgzf_show_compress (wvb, 0, wvb->comp_txt_data.len32, BGZF_EOF_LEN, NULL, 0, 0);

        buf_add_more (wvb, &wvb->comp_txt_data, BGZF_EOF, BGZF_EOF_LEN, "comp_txt_data");
    }
}

void bgzf_sign (uint64_t disk_size, uint8_t *signature)
{
    signature[0] = (disk_size      ) & 0xff; // LSB of size
    signature[1] = (disk_size >> 8 ) & 0xff;
    signature[2] = (disk_size >> 16) & 0xff;
}

// Entry point of BGZF compression compute thread.
// bgzf-compress vb->txt_data into vb->comp_txt_data - using BGZF blocks as prescribed in vb->gz_blocks. 
// Note: we hope to reconstruct the exact same byte-level BGZF blocks, as the original files, but that 
// will only happen if the GZIP library (eg libdeflate), version and parameters are the same 
static void bgzf_compress_vb (VBlockP vb)
{
    START_TIMER;

    if (flag.show_bgzf)
        iprintf ("COMPRESS thread=%s %s initialized <%sexact> re-compression with %s(%s[%u])\n",
                 threads_am_i_main_thread() ? "MAIN" : threads_am_i_writer_thread() ? "WRITER" : "COMPUTE", VB_NAME,
                 is_exact() ? "" : "non-",
                 codec_name (txt_file->effective_codec), bgzf_library_name (txt_file->mgzip_flags.library, true), txt_file->mgzip_flags.level);

    ASSERTNOTEMPTY (vb->gz_blocks);

    buf_alloc (vb, &vb->comp_txt_data, 0, vb->gz_blocks.len32 * BGZF_MAX_BLOCK_SIZE/2, uint8_t, 1, "comp_txt_data"); // alloc based on estimated size
    bgzf_alloc_compressor (vb, txt_file->mgzip_flags);

    for_buf2 (BgzfBlockPiz, block, i, vb->gz_blocks) {
        ASSERT (block->txt_index + block->txt_size <= Ltxt, 
                "block=%u out of range: expecting txt_index=%u txt_size=%u <= txt_data.len=%u",
                i, block->txt_index, block->txt_size, Ltxt);

        bgzf_compress_one_block (vb, Btxt (block->txt_index), block->txt_size, i, block->txt_index);
    }

    bgzf_free_compressor (vb, txt_file->mgzip_flags);

    vb_set_is_processed (vb); /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ 
    COPY_TIMER (bgzf_compute_thread);
}

#define BGZF_CREATED_BLOCK_SIZE 65280 // same size as observed in htslib-created files

// PIZ: calculate the BGZF blocks within this VB
static uint32_t bgzf_calculate_blocks_one_vb (VBlockP vb, bool is_last)
{
    // create our own equal-isize blocks
    if (!is_exact()) {
        buf_alloc_exact (vb, vb->gz_blocks, ceill ((double)Ltxt / BGZF_CREATED_BLOCK_SIZE), BgzfBlockPiz, "gz_blocks");
        
        for_buf2 (BgzfBlockPiz, blk, i, vb->gz_blocks)
            *blk = (BgzfBlockPiz){ .txt_index = i * BGZF_CREATED_BLOCK_SIZE, .txt_size = BGZF_CREATED_BLOCK_SIZE };
        
        BLST(BgzfBlockPiz, vb->gz_blocks)->txt_size -= (BGZF_CREATED_BLOCK_SIZE * vb->gz_blocks.len - Ltxt); // remove excessive length from last block
    
        return 0; // our gz_blocks perfectly cover the VB - no data remaining
    }

    // reconstruct --exact-ly based on mgzip_isizes
    else {
        ARRAY (uint32_t, isizes, txt_file->mgzip_isizes);
        uint32_t i, index=0;
        for (i=txt_file->mgzip_isizes.next; i < txt_file->mgzip_isizes.len; i++)
            if (index + isizes[i] <= Ltxt/*<= to include EOF in preceding VB */)
                index += isizes[i];
            else
                break;

        buf_alloc_exact (vb, vb->gz_blocks, i - txt_file->mgzip_isizes.next, BgzfBlockPiz, "gz_blocks");

        index = 0;
        for_buf (BgzfBlockPiz, blk, vb->gz_blocks) {
            *blk = (BgzfBlockPiz){ .txt_index = index, .txt_size = isizes[txt_file->mgzip_isizes.next++] };
            index += blk->txt_size;
        }

        int32_t remaining = Ltxt - index;
        ASSERT (IN_RANGE(remaining, 0, BGZF_MAX_BLOCK_SIZE), "mgzip_isizes exhausted prematurely: remaining=%d", remaining); // if we have 65536 or more remaining, there should have been more isizes

        return remaining;    
    }
}

// PIZ
void bgzf_dispatch_compress (Dispatcher dispatcher, STRp (uncomp), CompIType comp_i, bool is_last)
{
    // uncompressed data to be dealt with by next call to this function (buffer belongs to writer thread)
    static Buffer intercall_txt = {}; // belongs to wvb
    buf_alloc (wvb, &intercall_txt, 0, BGZF_MAX_BLOCK_SIZE, char, 1.5, "intercall_txt");

    uint32_t next_isize = txt_file->mgzip_isizes.len ? *B32(txt_file->mgzip_isizes, txt_file->mgzip_isizes.next) 
                                                     : BGZF_CREATED_BLOCK_SIZE;

    // case: uncomp is not enough to fill a block, just store it to next call
    if (!is_last && (uncomp_len + intercall_txt.len32 < next_isize)) {
        memcpy (BAFTc(intercall_txt), uncomp, uncomp_len);
        intercall_txt.len32 += uncomp_len;
        return;
    }

    if (uncomp_len || intercall_txt.len) { // might be 0 if is_last, in some cases

        VBlockP vb = dispatcher_generate_next_vb (dispatcher, wvb->vblock_i, COMP_NONE);
        vb->comp_i = comp_i;

        // build uncompressed data for this VB - some data left over from previous VB + data from wvb
        buf_alloc_exact (vb, vb->txt_data, intercall_txt.len + uncomp_len, char, "txt_data");
        if (intercall_txt.len32) memcpy (B1STtxt, intercall_txt.data, intercall_txt.len32);
        memcpy (Btxt (intercall_txt.len32), uncomp, uncomp_len);

        // calculate BGZF blocks - and trim data that doesn't fill a block - to be moved to next VB
        if ((intercall_txt.len32 = bgzf_calculate_blocks_one_vb (vb, is_last))) {
            Ltxt -= intercall_txt.len32;
            memcpy (B1STc(intercall_txt), BAFTtxt, intercall_txt.len32);
        }

        // BGZF-compress vb->txt_data in a separate thread
        dispatcher_compute (dispatcher, bgzf_compress_vb);
    }

    if (is_last) {
        dispatcher_set_no_data_available (dispatcher, false, DATA_EXHAUSTED);
        buf_destroy (intercall_txt);
    }
}

rom bgzf_library_name (MgzipLibraryType library, bool long_name)
{
    return (library < 0 || library >= NUM_ALL_BGZF_LIBRARIES) ? "INVALID_BGZF_LIBRARY"
         : long_name ? (rom[])BGZF_LIB_NAMES_LONG[library]
         :             (rom[])BGZF_LIB_NAMES_SHRT[library];
}

// used by test/Makefile
void il1m_compress (void)
{
    void *compressor = libdeflate_alloc_compressor (evb, flag.fast?1 : flag.best?9 : 5, __FUNCLINE); 

    uint8_t *in = MALLOC (1 MB), *out = MALLOC (2 MB);
    uint32_t in_len;
    for (int i=0; (in_len = fread (in, 1, 1 MB, stdin)); i++) {
        GzipFooter footer = { .crc32 = LTEN32 (crc32 (0, in, in_len)),
                              .isize = LTEN32 (in_len) };

        uint32_t out_len = libdeflate_deflate_compress (compressor, in, in_len, out, 2 MB);
        ASSERT (out_len, "deflate failed: in_len=%u block_i=%u", in_len, i);

        ASSERT0 (1 == fwrite (_S(IL1M_HEADER), 1, stdout), "fwrite failed #1");
        ASSERT0 (1 == fwrite (flag.best?"\x02\x03" : flag.fast?"\x04\x03" : "\x00\x03", 2, 1, stdout), "fwrite failed #2");        
        ASSERT  (1 == fwrite (STRa(out), 1, stdout), "fwrite failed: #3 out_len=%u", out_len);
        ASSERT0 (1 == fwrite (&footer, sizeof (footer), 1, stdout), "fwrite failed #4");
    }

    fflush (stdout);
    exit (0);
}
