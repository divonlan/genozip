// ------------------------------------------------------------------
//   bgzf.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <errno.h>

#include "igzip/igzip_lib.h"
#include "libdeflate_1.19/libdeflate.h"
#include "libdeflate_1.7/libdeflate.h"
#include "zlib/zlib.h"
#include "bgzf.h"
#include "endianness.h"
#include "buffer.h"
#include "vblock.h"
#include "arch.h"
#include "strings.h"
#include "file.h"
#include "zfile.h"
#include "zip.h"
#include "arch.h"
#include "txtfile.h"
#include "codec.h"
#include "threads.h"
#include "segconf.h"
#include "dispatcher.h"
#include "writer.h"
#include "gencomp.h"
#include "filename.h"

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
    uint8_t si1;    // BGZF id - must be 66  (0x42)
    uint8_t si2;    // BGZF id - must be 67  (0x43)
    uint16_t slen;  // BGZF extra field length - must be 2
    uint16_t bsize; // BGZF extra field - (compressed block size -1)
} BgzfHeader;

typedef struct BgzfFooter {
    uint32_t crc32; // CRC32 of uncompressed data
    uint32_t isize; // Input (i.e. uncompressed) Size
} BgzfFooter;

static FlagsBgzf bgzf_recompression_levels[1+MAX_FLAG_BGZF] = {
    {                               .level = 0  },  // --bgzf=0 : CODEC_NONE, not BGZF     
    { .library = BGZF_IGZIP,        .level = 1  },  // --bgzf=1 : note: this is IGZIP LVL0 
    { .library = BGZF_IGZIP,        .level = 2  },  // --bgzf=2 : note: this is IGZIP LVL1
    { .library = BGZF_LIBDEFLATE19, .level = 1  },  // --bgzf=3 
    { .library = BGZF_LIBDEFLATE19, .level = 7  },  // --bgzf=4 
    { .library = BGZF_LIBDEFLATE19, .level = 9  },  // --bgzf=5
};
#define BGZF_DEFAULT_LEVEL 2 // used if --bgzf is not specified (it is actually faster than 1 if also writing to disk)

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
static BgzfBlockStr display_bb (BgzfBlockZip *bb)
{
    BgzfBlockStr s;
    sprintf (s.s, "{txt_index=%u txt_size=%u compressed_index=%u comp_size=%u is_decompressed=%u}",
             bb->txt_index, bb->txt_size, bb->compressed_index, bb->comp_size, bb->is_decompressed);
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

    if (file->codec == CODEC_GZ) {
        if (flag.show_gz) {
            iprintf ("%s: is GZIP but not BGZF\n", file->name); fflush (info_stream);
            exit_ok;
        }
        else return;
    }
 
    else if (file->codec != CODEC_BGZF) {
        if (flag.show_gz) {
            iprintf ("%s: is not GZIP, it is %s\n", file->name, codec_name (file->source_codec)); 
            exit_ok;
        }
        else return;
    }

    ARRAY_alloc (FlagsBgzf, ll, 13+12+9, false, file->bgzf_plausible_levels, evb, "txt_file->bgzf_plausible_levels");

    int next=0;
    for (int l=0; l <= 12; l++) // level=0 only here, bc it would be the same in all libraries
        ll[next++] = (FlagsBgzf){ .library = BGZF_LIBDEFLATE19, .level = l};

    for (int l=1; l <= 12; l++)
        ll[next++] = (FlagsBgzf){ .library = BGZF_LIBDEFLATE7,  .level = l};

    for (int l=1; l <= 9; l++)
        ll[next++] = (FlagsBgzf){ .library = BGZF_ZLIB,         .level = l};
}

// ZIP main thread
static void bgzf_discover_finalize_testing (BgzfLibraryType lib, BgzfLevel level)
{
    txt_file->bgzf_flags.library = lib;  // assign field-by-field: careful not to modify bgzf_flags.has_eof_block 
    txt_file->bgzf_flags.level   = level;

    if (flag.zip_comp_i < MAX_NUM_COMPS) // for stats
        z_file->comp_bgzf[flag.zip_comp_i] = txt_file->bgzf_flags;
}

// ZIP main thread
void bgzf_finalize_discovery (void)
{    
    int n_levels = txt_file->bgzf_plausible_levels.len32;

    // case: there is no library/level combination for which we can decompress with bgzf=exact
    if (n_levels == 0) {
        bgzf_discover_finalize_testing (0, BGZF_COMP_LEVEL_UNKNOWN);

        if (flag.show_bgzf || flag.show_gz) 
            iprintf ("%s: is a BGZF file, generated by an unidentified library\n", txt_name);
    }

    // case: one or more library/level combinations was verified with all test bgzf blocks (10 blocks, unless file is shorter) 
    else {
        bgzf_discover_finalize_testing (B1ST(FlagsBgzf, txt_file->bgzf_plausible_levels)->library, B1ST(FlagsBgzf, txt_file->bgzf_plausible_levels)->level);

        if (flag.show_bgzf || flag.show_gz) 
            iprintf ("%s: %s %s level %u\n", txt_name, 
                     (n_levels == 1) ? "Identified as generated with" : "Multiple plausible levels, arbitrarily selecting",
                     bgzf_library_name (txt_file->bgzf_flags.library, true), txt_file->bgzf_flags.level);
    }    

    if (flag.show_gz) exit_ok;
}

// ZIP: test a BGZF block against all the remaining plausible levels, and eliminate those that don't match. 
static void bgzf_discover_library_and_level (VBlockP vb, int test_block_i, STRp(comp), STRp(uncomp))
{
    if (comp_len <= sizeof (BgzfHeader) + sizeof (BgzfFooter)) {
        txt_file->bgzf_plausible_levels.len = 0;

        if (flag.show_bgzf || flag.show_gz) 
            iprintf ("%s: Block too small - could not identify compression library and level\n", txt_name);

        if (flag.show_gz) exit_ok;

        return;
    }

    // ignore the header and footer of the block
    comp     += sizeof (BgzfHeader);
    comp_len -= sizeof (BgzfHeader) + sizeof (BgzfFooter);

    // compress with each of the remaining plausible levels - testing if the compression is identical to the actual
    STRl (recomp, BGZF_MAX_BLOCK_SIZE); 

    for_buf (FlagsBgzf, ll, txt_file->bgzf_plausible_levels) {

        switch (ll->library) {
            case BGZF_LIBDEFLATE19 : {
                void *compressor = libdeflate_alloc_compressor (vb, ll->level, __FUNCLINE);
                recomp_len = (uint32_t)libdeflate_deflate_compress (compressor, STRa(uncomp), recomp, BGZF_MAX_BLOCK_SIZE);

                libdeflate_free_compressor (compressor, __FUNCLINE);
                break;
            }

            case BGZF_LIBDEFLATE7 : {
                void *compressor = libdeflate_alloc_compressor_1_7 (ll->level, vb);
                recomp_len = (uint32_t)libdeflate_deflate_compress_1_7 (compressor, STRa(uncomp), recomp, BGZF_MAX_BLOCK_SIZE);

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
                strm.avail_out = sizeof (recomp);
                ASSERT (deflate (&strm, Z_FINISH) == Z_STREAM_END, "deflate failed: msg=%s", strm.msg);

                recomp_len = sizeof (recomp) - strm.avail_out;
                
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
            buf_remove (txt_file->bgzf_plausible_levels, FlagsBgzf, BNUM(txt_file->bgzf_plausible_levels, ll), 1);
            ll--; fb_after--; // hack for_buf loop
        }
    }
}

//--------------------------------------------------------------------
// ZIP SIDE - decompress BGZF-compressed file and prepare BGZF section
//--------------------------------------------------------------------

// ZIP: reads and validates a BGZF block, and returns the uncompressed size or (only if soft_fail) an error
static int32_t bgzf_read_block_raw (FILE *file, // txt_file is not yet assigned when called from file_open_txt_read
                                    uint8_t *block /* must be BGZF_MAX_BLOCK_SIZE in size */, uint32_t *block_size /* out */,
                                    rom basename, bool is_remote, FailType soft_fail,
                                    int64_t *disk_so_far) 
{
    BgzfHeader *h = (BgzfHeader *)block;

    // read the header
    *block_size = fread (h, 1, sizeof (struct BgzfHeader), file);

    if (disk_so_far) *disk_so_far += *block_size;   

    if (*block_size != sizeof (struct BgzfHeader)) {
        ASSERT (!ferror (file), //|| (disk_so_far && *disk_so_far == *block_size), 
                "Error while reading %s: %s", basename, strerror (errno));
        return (disk_so_far && *disk_so_far == *block_size) ? BGZF_BLOCK_IS_NOT_GZIP : BGZF_ABRUBT_EOF ; // EOF without an EOF block (possibly a very short non-GZ file)
    }
    
    if (*block_size < 12) {
        ASSERT (soft_fail, "file %s appears truncated - it ends with a partial gzip block header", basename); // less than the minimal gz block header size
        return BGZF_BLOCK_IS_NOT_GZIP;
    } 

    // case: this is not a GZ / BGZF block at all (see: https://tools.ietf.org/html/rfc1952)
    if (h->id1 != 31 || h->id2 != 139) {
        ASSERT (soft_fail, "expecting %s to be compressed with gzip format, but it is not", basename);
        return BGZF_BLOCK_IS_NOT_GZIP;
    }

    // case: this is GZIP block that is NOT a valid BGZF block (see: https://samtools.github.io/hts-specs/SAMv1.pdf)
    if (flag.no_bgzf || // user instructed us to treat BGZF data as normal GZIP data
        (!(*block_size == sizeof (struct BgzfHeader) && !memcmp (h, BGZF_PREFIX, BGZF_PREFIX_LEN)))) {
        ASSINP (soft_fail, "Encountered a GZIP block that unexpectedly is not BGZF in %s offset=%"PRIu64"\nSolution: use --no-bgzf", 
                basename, (uint64_t)ftello64 (file) - *block_size);

        return BGZF_BLOCK_GZIP_NOT_BGZIP;
    }

    *block_size = LTEN16 (h->bsize) + 1;
    
    uint32_t body_size = *block_size - sizeof (struct BgzfHeader);
    uint32_t bytes = fread (h+1, 1, body_size, file);

    if (disk_so_far) *disk_so_far += bytes;   

    int save_errno = errno; // we want to report errno of fread, not ftell.

    // if failed, always error, even if soft_fail
    ASSERT (bytes == body_size || flag.truncate, 
            "%s %s (ftell=%"PRId64" err=\"%s\" bytes_read=%u but expecting=%u filesystem=%s). %s\n", 
            feof (file) ? "Unexpected end of file while reading" : "Failed to read body", 
            basename, ftello64 (file), 
            (is_remote && save_errno == ESPIPE) ? "Disconnected from remote host" : strerror (save_errno),
            bytes, body_size, arch_get_filesystem_type().s,
            feof (file) ? "If file is expected to be truncated, you may use --truncate-partial-last-line to disregard the final partial BGZF block." : "");
    
    return (bytes == body_size) ? BGZF_BLOCK_SUCCESS : BGZF_BLOCK_TRUNCATED;
}

// ZIP main thread: reads and validates a BGZF block, and returns the uncompressed size or (only if soft_fail) an error
int32_t bgzf_read_block (FileP file, // txt_file is not yet assigned when called from file_open_txt_read
                         uint8_t *block /* must be BGZF_MAX_BLOCK_SIZE in size */, uint32_t *block_size /* out */,
                         FailType soft_fail)
{
    START_TIMER;

    int ret = bgzf_read_block_raw ((FILE *)file->file, block, block_size, file->basename, file->is_remote, soft_fail, &file->disk_so_far);
    if (ret == BGZF_BLOCK_IS_NOT_GZIP || ret == BGZF_BLOCK_GZIP_NOT_BGZIP) 
        return ret; // happens only if soft_fail
    
    if (ret == BGZF_ABRUBT_EOF) {
        ASSERT (!file->disk_size                ||    // redirected or remote 
                flag.truncate ||    // possibly compressing while downloading
                file->disk_so_far == file->disk_size, // entire file was read
                "Abrupt EOF in BGZF file %s: disk_so_far=%s disk_size=%s filesystem=%s", 
                file->name, str_int_commas (file->disk_so_far).s, str_int_commas (file->disk_size).s, arch_get_filesystem_type().s);

        return 0; // no EOF block, that's fine
    }

    else if (ret == BGZF_BLOCK_TRUNCATED) {
        if (txt_file->bgzf_truncated_last_block) // we arrive here twice - show warning only on the second time
            WARN ("FYI: %s is truncated - its final BGZF block in incomplete. Dropping this defective BGZF block.", txt_name);
        txt_file->bgzf_truncated_last_block = true;
        return 0;
    }

    uint32_t isize_lt32 = GET_UINT32 (&block[*block_size - 4]);
    uint32_t isize = LTEN32 (isize_lt32); // 0...65536 per spec: "input size" = uncompressed data length
    ASSERT (isize <= 65536, "isize=%u out of range [0,65536]", isize);
    
    // add isize to buffer that will be written to SEC_BGZF
    if (isize) { // don't store EOF block (bc isize=0 cannot be represented as (isize-1) )
        #define BGZF_INITIAL_ALLOC 16 // just of the sake of a bit of effeciency: 16 chosen carefully so 16*63000 < 1MB min vb_size but over segconf size
        if (file->bgzf_isizes.len32 <= BGZF_INITIAL_ALLOC) { // entered thrice: when called from file_open_txt_read, segconf, and in first VB 
            buf_alloc (evb, &file->bgzf_isizes, 0, MAX_(BGZF_INITIAL_ALLOC, segconf.vb_size / 63000), uint16_t, 0, "txt_file->bgzf_isizes");
            buf_alloc (evb, &file->bgzf_starts, 0, MAX_(BGZF_INITIAL_ALLOC, segconf.vb_size / 63000), uint64_t, 0, "txt_file->bgzf_starts");
        }

        buf_append_one (file->bgzf_isizes, BGEN16 ((uint16_t)(isize - 1))); // -1 to make the range 0..65535
        buf_append_one (file->bgzf_starts, txt_file ? txt_file->disk_so_far - *block_size : 0); // not BGEN bc not written to z_file. note: first block is read from file_open_txt_read before txt_file is assigned
    }
    else {
        // if isize is 0, we're expecting an EOF block
        ASSERT (str_issame_((rom)block, *block_size, BGZF_EOF, BGZF_EOF_LEN),
                "Corrupt BGZF block in %s offset=%"PRIu64" bgzf_block_size=%u: isize=0 but this is not an EOF block",
                file->name, file->disk_so_far - *block_size, *block_size);

        ASSERT (file->disk_so_far == file->disk_size || // expected
                !file->disk_size || file_is_read_via_ext_decompressor (file) || file->redirected || file->is_remote, // cases in which we can't reliably test this condition
                "Corrupt BGZF file %s (size=%"PRIu64"): BGZF EOF block encountered at offset=%"PRIu64" length=%u, but this is not the end of the file",
                file->name, file->disk_size, file->disk_so_far - *block_size, *block_size);

        if (txt_file) 
            txt_file->bgzf_flags.has_eof_block = true;
    }

    COPY_TIMER_EVB (bgzf_read_block);
    return isize;
}

// ZIP: BGZF section per component
void bgzf_compress_bgzf_section (void)
{
    // cases where we don't write the BGZF blocks section
    if (!txt_file->bgzf_isizes.len ||  // this txt file is not compressed with BGZF - we don't need a BGZF section
        txt_file->bgzf_flags.level == BGZF_COMP_LEVEL_UNKNOWN ||  // we don't know the level - so PIZ will reconstruct at default level
        flag.zip_txt_modified)         // the file has changed and we can't reconstruct to the same blocks
        return;

    // sanity check
    int64_t total_isize = 0;
    for_buf (uint16_t, isize_p, txt_file->bgzf_isizes)
        total_isize += BGEN16 (*isize_p) + 1; // values 0-65535 correspond to isize 1-65536
    
    ASSERT (total_isize == txt_file->txt_data_so_far_single + txt_file->last_truncated_line_len, "Expecting total_isize=%"PRId64" == txt_file->txt_data_so_far_single=%"PRId64,
            total_isize, txt_file->txt_data_so_far_single);

    // get the best codec for the SEC_BGZF section
    txt_file->bgzf_isizes.len *= sizeof (uint16_t);
    Codec codec = codec_assign_best_codec (evb, NULL, &txt_file->bgzf_isizes, SEC_BGZF);

    evb->comp_i = flag.zip_comp_i; // this goes into SectionEntFileFormat.comp_i via sections_add_to_list
    zfile_compress_section_data_ex (evb, NULL, SEC_BGZF, &txt_file->bgzf_isizes, NULL, 0, codec, (SectionFlags)txt_file->bgzf_flags, NULL);
    txt_file->bgzf_isizes.len /= sizeof (uint16_t); // restore
}

// uncompresses a BGZF block in vb->scratch referred to by bb, into its place in vb->txt_data as prescribed by bb
// might be called from main thread or compute threads
void bgzf_uncompress_one_block (VBlockP vb, BgzfBlockZip *bb)
{
    if (bb->is_decompressed) return; // already decompressed - nothing to do

    ASSERT0 (vb->gzip_compressor, "vb->gzip_compressor=NULL");

    BgzfHeader *h = (BgzfHeader *)Bc(vb->scratch, bb->compressed_index);

    // verify that entire block is within vb->scratch
    ASSERT (bb->compressed_index + sizeof (BgzfHeader) < vb->scratch.len && // we have at least the header - we can access bsize
            bb->compressed_index + (uint32_t)LTEN16 (h->bsize) + 1 <= vb->scratch.len, 
            "%s: BGZF block size goes past the end of in vb->scratch: bb=%s compressed_index=%u vb->scratch.len=%"PRIu64, 
            VB_NAME, display_bb (bb).s, bb->compressed_index, vb->scratch.len);

    ASSERT (h->id1==31 && h->id2==139, "%s: invalid BGZF block in vb->scratch: compressed_index=%u", VB_NAME, bb->compressed_index);

    if (flag.show_bgzf)
        iprintf ("%-7s %s i=%u compressed_index=%u size=%u txt_index=%u size=%u ",
                 threads_am_i_main_thread() ? "MAIN" : "COMPUTE", VB_NAME, 
                 BNUM (vb->bgzf_blocks, bb), bb->compressed_index, bb->comp_size, bb->txt_index, bb->txt_size);

    enum libdeflate_result ret = 
        libdeflate_deflate_decompress (vb->gzip_compressor, 
                                       h+1, bb->comp_size - sizeof(BgzfHeader) - sizeof (BgzfFooter), // compressed
                                       Btxt (bb->txt_index), bb->txt_size, NULL);  // uncompressed

    ASSERT (ret == LIBDEFLATE_SUCCESS, "libdeflate_deflate_decompress failed: %s", libdeflate_error(ret));

    bb->is_decompressed = true;

    if (flag.show_bgzf)
        #define C(i) ((bb->txt_index + i < Ltxt) ? char_to_printable (*Btxt (bb->txt_index + (i))).s : "") 
        iprintf ("txt_data[5]=%1s%1s%1s%1s%1s %s\n", C(0), C(1), C(2), C(3), C(4), bb->comp_size == BGZF_EOF_LEN ? "EOF" : "");
        #undef C

    // discover which gzip library and compression level were used (testing the first few BGZF blocks)
    int test_block_i=0;
    if (txt_file->bgzf_plausible_levels.len32 &&
        txt_file->bgzf_plausible_levels.count < BGZF_DISCOVERY_MAX_TESTS && // fail fast without atomic
        (test_block_i = __atomic_fetch_add (&txt_file->bgzf_plausible_levels.count, 1, __ATOMIC_RELAXED)) < BGZF_DISCOVERY_MAX_TESTS) { // note: if multiple threads test in parallel, count might be incremented beyond BGZF_DISCOVERY_MAX_TESTS - that's ok
        
        bgzf_discover_library_and_level (vb, test_block_i, (rom)h, bb->comp_size, Btxt (bb->txt_index), bb->txt_size);        
    }

    // in case of --show_gz: report and exit (otherwise, we finalize in the main thread to avoid thread issues with updating txt_file->bgzf_flags)
    if (flag.show_gz && (test_block_i == BGZF_DISCOVERY_MAX_TESTS-1 || !txt_file->bgzf_plausible_levels.len32)) 
        bgzf_finalize_discovery();
}

// ZIP: called from the compute thread: zip_compress_one_vb and main thread: txtfile_read_block_bgzf
void bgzf_uncompress_vb (VBlockP vb)
{
    START_TIMER;

    vb->gzip_compressor = libdeflate_alloc_decompressor(vb, __FUNCLINE);

    for_buf (BgzfBlockZip, bb, vb->bgzf_blocks)
        bgzf_uncompress_one_block (vb, bb);

    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor, __FUNCLINE);

    buf_free (vb->scratch); // now that we are finished decompressing we can free it

    if (flag.show_time) {
        if (threads_am_i_main_thread ()) COPY_TIMER (bgzf_io_thread)
        else                             COPY_TIMER (bgzf_compute_thread);
    }

    COPY_TIMER (bgzf_uncompress_vb);
}

// ZIP: decompresses a prescribed BGZF block when re-reading DEPN lines
static inline void bgzf_uncompress_one_prescribed_block (VBlockP vb, STRp(bgzf_block), STRc (uncomp_block), uint64_t bb_i)
{
    START_TIMER;

    BgzfHeader *h = (BgzfHeader *)bgzf_block;

    // verify that entire block is within vb->scratch
    ASSERT ((uint32_t)LTEN16 (h->bsize) + 1 == bgzf_block_len, 
            "%s: BGZF reread: expecting block_size=%u but found %u", VB_NAME, bgzf_block_len, LTEN16 (h->bsize) + 1);

    ASSERT (h->id1==31 && h->id2==139, "%s: invalid BGZF block", VB_NAME);

    if (flag.show_bgzf)
        iprintf ("REREAD  vb=%s reread bb_i=%"PRIu64" comp_size=%u uncomp_size=%u ", 
                 VB_NAME, bb_i, bgzf_block_len, uncomp_block_len);

    enum libdeflate_result ret = 
        libdeflate_deflate_decompress (vb->gzip_compressor, 
                                       h+1, bgzf_block_len - sizeof(BgzfHeader) - sizeof (BgzfFooter), // compressed
                                       STRa(uncomp_block), NULL);  // uncompressed

    ASSERT (ret == LIBDEFLATE_SUCCESS, "%s: libdeflate_deflate_decompress failed: %s. bgzf_block_len=%u uncomp_block_len=%u bb_i=%"PRIu64, 
            VB_NAME, libdeflate_error(ret), bgzf_block_len, uncomp_block_len, bb_i);

    if (flag.show_bgzf)
        #define C(i) (i < uncomp_block_len ? char_to_printable (uncomp_block[i]).s : "") 
        iprintf ("txt_data[5]=%1s%1s%1s%1s%1s\n", C(0), C(1), C(2), C(3), C(4));
        #undef C

    COPY_TIMER (bgzf_uncompress_one_prescribed_block);
}

// ZIP: compute thread of a DEPN VB: actually re-reading data into txt_data according to vb->reread_prescription
void bgzf_reread_uncompress_vb_as_prescribed (VBlockP vb, FILE *file)
{
    uint64_t last_offset = -1LL;
    char uncomp_block[BGZF_MAX_BLOCK_SIZE];

    vb->gzip_compressor = libdeflate_alloc_decompressor(vb, __FUNCLINE);

    for_buf (RereadLine, line, vb->reread_prescription) {
        
        // a line might span 1 or more BGZF blocks
        while (line->line_len) { 

            uint64_t offset = *B64 (txt_file->bgzf_starts, line->offset.bb_i);
            uint32_t isize  = BGEN16 (*B16 (txt_file->bgzf_isizes, line->offset.bb_i)) + 1; // maximum isize is 65536 (not 65535)

            if (offset != last_offset) {
                ASSERT (!fseeko64 (file, offset, SEEK_SET),
                        "%s: fseeko64 on %s failed while rereading BGZF depn lines: %s", VB_NAME, txt_name, strerror(errno));

                STRl (bgzf_block, BGZF_MAX_BLOCK_SIZE);
                int32_t ret = bgzf_read_block_raw (file, (uint8_t*)qSTRa(bgzf_block), txt_file->basename, false, HARD_FAIL, NULL); 
                ASSERT (ret != BGZF_ABRUBT_EOF, "Unexpected BGZF_ABRUBT_EOF while re-reading BGZF block in %s: filesystem=%s offset=%"PRIu64" uncomp_block_size=%u", 
                        txt_name, arch_get_filesystem_type().s, offset, isize);

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
void bgzf_zip_advance_index (VBlockP vb, uint32_t line_len)
{
    if (!vb->bgzf_blocks.len) return; // no BGZF blocks in this VB - all data came from "unconsumed_txt"

    vb->line_bgzf_uoffset += line_len;

    // udpate current_bb_i and bgzf_offset (note: line_len might span multiple bgzf blocks)
    BgzfBlockZip *bb;
    for (bb = B(BgzfBlockZip, vb->bgzf_blocks, vb->bgzf_blocks.current_bb_i); 
         vb->line_bgzf_uoffset && vb->line_bgzf_uoffset >= bb->txt_size; // note: careful to also terminate on the edge case that line_bgzf_uoffset==0 and in the final VB block bb->txt_size==0 
         bb++) 
      
        vb->line_bgzf_uoffset -= bb->txt_size; // index into the next BGZF block

    vb->bgzf_blocks.current_bb_i = BNUM(vb->bgzf_blocks, bb);
}

// ZIP: after reading data for a txt_header or VB, copy unconsumed bgzf_blocks to txt_file->unconsumed_bgzf_blocks
// The first block might be partially consumed.
int64_t bgzf_copy_unconsumed_blocks (VBlockP vb)
{
    START_TIMER;

    if (!vb->bgzf_blocks.len) return 0; // not a BGZF-compressed file

    int32_t consumed = Ltxt +   // amount of data consumed by this VB
                       vb->bgzf_blocks.consumed_by_prev_vb; 

    ARRAY (BgzfBlockZip, bb, vb->bgzf_blocks);

    bool done = false;
    bool consumed_full_bgzf_blocks=false;
    int64_t compressed_size = 0;

    for (uint32_t i=0; i < bb_len; i++) {
        // if some of the BGZF blocks are not consumed (the first of them might be partially consumed) - move the blocks
        // to unconsumed_bgzf_blocks - to be moved to the next VB
        if (consumed - bb[i].txt_size < 0 && !done/*enter only once*/) {
            
            consumed_full_bgzf_blocks = (consumed == 0); // no partially-consumed block

            // block i might be partially consumed or not consumed at all, subsequent blocks are not consumed at all
            buf_append (evb, txt_file->unconsumed_bgzf_blocks, BgzfBlockZip, 
                        B(BgzfBlockZip, vb->bgzf_blocks, i), vb->bgzf_blocks.len32 - i, "txt_file->unconsumed_bgzf_blocks");

            txt_file->unconsumed_bgzf_blocks.consumed_by_prev_vb = consumed; // part of first BGZF block already consumed
            done = true;
        }
        else if (!done)
            compressed_size += bb[i].comp_size;

        consumed -= bb[i].txt_size;
    }

    // sanity check
    ASSERT (-consumed == txt_file->unconsumed_txt.len32, "Expecting (-consumed)=%d == unconsumed_txt.len=%u", -consumed, txt_file->unconsumed_txt.len32);

    COPY_TIMER (bgzf_copy_unconsumed_blocks);
    return consumed_full_bgzf_blocks ? compressed_size : 0;
} 

// return blocks used by the segconf VB to the unconsumed blocks
void bgzf_return_segconf_blocks (VBlockP vb)
{
    buf_copy (evb, &txt_file->unconsumed_bgzf_blocks, &vb->bgzf_blocks, BgzfBlockZip, 0, 0, 0);
    txt_file->unconsumed_bgzf_blocks.consumed_by_prev_vb = vb->bgzf_blocks.consumed_by_prev_vb;
}

// ZIP: before reading data for a VB, populate bgzf_blocks with some or all of the unconsumed blocks passed
// from the previous VB or txt_header
void bgzf_zip_init_vb (VBlockP vb)
{
    vb->vb_bgzf_i = txt_file->bgzf_isizes.len + txt_file->bgzf_flags.has_eof_block; // index of first bgzf block to be used by the VB

    if (!txt_file->unconsumed_bgzf_blocks.len) return; // happens when either unconsumed_bytes=0 or not a BGZF-compressed file

    // data in the first BGZF block already consumed by previous VB or txt_header 
    vb->bgzf_blocks.consumed_by_prev_vb = vb->line_bgzf_uoffset = txt_file->unconsumed_bgzf_blocks.consumed_by_prev_vb;

    // copy all unconsumed BGZF blocks - we might not need all of them - the unconsumed ones will moved back in bgzf_copy_unconsumed_blocks
    buf_copy (vb, &vb->bgzf_blocks, &txt_file->unconsumed_bgzf_blocks, BgzfBlockZip, 0, 0, "bgzf_blocks"); 

    vb->vb_bgzf_i -= txt_file->unconsumed_bgzf_blocks.len32;

    txt_file->unconsumed_bgzf_blocks.len32 = txt_file->unconsumed_bgzf_blocks.consumed_by_prev_vb = 0;

    // sanity check
    int32_t available = -vb->bgzf_blocks.consumed_by_prev_vb; // possibly start negative
    for_buf (BgzfBlockZip, bb, vb->bgzf_blocks) 
        available += bb->txt_size;

    ASSERT (available >= Ltxt, "BGZF blocks in txt_file->unconsumed_bgzf_blocks cover only %d bytes, less than the needed unconsumed_bytes=%d", 
            available, Ltxt);
}

//-----------------------------------------------------
// PIZ SIDE - setting up BGZF for a particular txt file
//-----------------------------------------------------

static Buffer isizes = {}; // will be grabbed into txt_file->bgzf_isizes;

static FlagsBgzf recompression_template (int bgzf_level)
{
    if (bgzf_level < 0 || bgzf_level > MAX_FLAG_BGZF) 
        bgzf_level = BGZF_DEFAULT_LEVEL;

    return (FlagsBgzf){ .has_eof_block = true, 
                        .level         = bgzf_recompression_levels[bgzf_level].level, // a 4-bit bitfield 
                        .library       = bgzf_recompression_levels[bgzf_level].library };
}

static FlagsBgzf bgzf_load_isizes (CompIType comp_i, bool show_only) 
{
    Section sec = sections_get_comp_bgzf_sec (comp_i);
    if (!sec) goto fallback; // this component doesn't contain a BGZF section

    int32_t offset = zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", SEC_BGZF, sec);

    SectionHeaderP header = (SectionHeaderP)Bc(evb->z_data, offset);

    // if we don't know the compression level (in older Genozip versions we wrote the SECTION_BGZF even 
    // if level discover failed), or if original file had compression level 0 (no compression), go with the default
    if (header->flags.bgzf.level == BGZF_COMP_LEVEL_UNKNOWN || !header->flags.bgzf.level) 
        goto fallback;

    zfile_uncompress_section (evb, header, &isizes, "txt_file->bgzf_isizes", 0, SEC_BGZF);
    isizes.len /= 2;

    if (show_only) buf_destroy (isizes);

    // convert to native endianity from big endian
    for_buf (uint16_t, isize_p, isizes)
        *isize_p = BGEN16 (*isize_p); // now it contains isize-1 - a value 0->65535 representing an isize 1->65536

    return header->flags.bgzf; // bgzf_isizes successfully loaded

fallback:
    return recompression_template (BGZF_DEFAULT_LEVEL);
}          

#define HAS_EXT(ext) filename_has_ext (flag.out_filename, (ext))
static bool flags_out_filename_implies_bgzf (void)
{
    return flag.out_filename && (HAS_EXT(".gz") || HAS_EXT(".bgz") || (Z_DT(SAM) && HAS_EXT(".bam")) || (Z_DT(VCF) && HAS_EXT(".bcf")));
}

// PIZ: called from main thread after reading txt_header's header
FlagsBgzf bgzf_piz_calculate_bgzf_flags (CompIType comp_i, Codec src_codec)
{
    FlagsBgzf bgzf_flags = recompression_template (flag.bgzf); // set to --bgzf command line value (set to BGZF_COMP_LEVEL_UNKNOWN if BGZF_NOT_INITIALIZED or BGZF_BY_ZFILE)

    if (flag.test)
        {} // bgzf_flags.level remains BGZF_COMP_LEVEL_UNKNOWN (conflict detected if -z is used with -t)

    // case: already set by command line (to a level or "exact")
    else if (flag.bgzf != BGZF_NOT_INITIALIZED) {
        // if user specified --bgzf and --output - make sure output filename is .gz or .bam
        ASSINP (!flag.out_filename || flags_out_filename_implies_bgzf() || bgzf_flags.level==0, 
                "using %s in combination with %s for outputting a %s file, requires the output filename to end with %s", 
                OT("output", "o"), OT("bgzf", "z"), dt_name(flag.out_dt), OUT_DT(BAM)?".bam" : OUT_DT(BCF)?".bcf" : ".gz");
    
        ASSINP0 (!OUT_DT(BCF) || flag.bgzf != BGZF_BY_ZFILE, "cannot use --bgzf=exact when outputing a BCF file"); // because we have no control over bcftools' BGZF block generation
    }

    // case: genocat to stdout - no BGZF, except BAM (bc some downstream tools (eg GATK) expect BAM to always be BGZF-compressed)
    else if (is_genocat && !flag.out_filename)
        bgzf_flags = bgzf_recompression_levels[OUT_DT(BAM) ? 1 : 0];

    // case: out_filename - determine by file name
    else if (flag.out_filename)
        bgzf_flags = bgzf_recompression_levels[flags_out_filename_implies_bgzf() ? BGZF_DEFAULT_LEVEL : 0];
    
    // case: genounzip without out_filename - determined by source file (0 or BGZF_DEFAULT_LEVEL) - 
    // set in txtheader_piz_read_and_reconstruct
    else 
        {}

    // case: attempt to reconstruct the BGZF following the instructions from the z_file
    if (flag.bgzf == BGZF_BY_ZFILE) // note: flags_update_piz_one_z_file enforce that !flag.zip_txt_modified in this case
        bgzf_flags = bgzf_load_isizes (comp_i, false); 
    
    // case: user wants to see this section header, despite not needing BGZF data
    else if (flag.only_headers == SEC_BGZF+1 || flag.only_headers == SHOW_ALL_HEADERS)
        bgzf_load_isizes (comp_i, true); 
        
    // case: bgzf_flags.level is to left to be determined by whether the source file was compressed (see bgzf_piz_calculate_bgzf_flags)
    if (bgzf_flags.level == BGZF_COMP_LEVEL_UNKNOWN)
        #define C(cdc) (src_codec == CODEC_##cdc)
        // note: for bz2, xz, and zip - we reconstruct as gz too. better choice than plain.
        bgzf_flags = bgzf_recompression_levels[(C(BGZF) || C(BAM) || C(GZ) || C(BZ2) || C(XZ) || C(ZIP)) ? BGZF_DEFAULT_LEVEL : 0]; // note: similar logic to txtheader_get_txt_filename_from_section
        #undef C
        
    return bgzf_flags;
}

// PIZ main thread: update txt_file with BGZF info calculated earlier
void bgzf_piz_set_txt_file_bgzf_info (FlagsBgzf bgzf_flags, bytes codec_info)
{
    memcpy (txt_file->bgzf_signature, codec_info, 3);
    
    if (isizes.len) 
        buf_grab (evb, txt_file->bgzf_isizes, "txt_file->bgzf_isizes", isizes);
        
    txt_file->bgzf_flags = bgzf_flags;

    // sanity        
    ASSERT (txt_file->bgzf_flags.level >= 0 && txt_file->bgzf_flags.level <= 12, "txt_file->bgzf_flags.level=%u out of range [0,12]", 
            txt_file->bgzf_flags.level);

    ASSERT (txt_file->bgzf_flags.library >= 0 && txt_file->bgzf_flags.library < NUM_BGZF_LIBRARIES, "txt_file->bgzf_flags.library=%u out of range [0,%u]", 
            txt_file->bgzf_flags.level, NUM_BGZF_LIBRARIES-1);
}

//-----------------------------------------------------
// PIZ SIDE - compressing txt_file with BGZF
//-----------------------------------------------------

static void bgzf_alloc_compressor (VBlockP vb, FlagsBgzf bgzf_flags)
{
    ASSERT0 (!vb->gzip_compressor, "expecting vb->gzip_compressor=NULL");

    switch (bgzf_flags.library) {
        case BGZF_LIBDEFLATE19:
            vb->gzip_compressor = libdeflate_alloc_compressor (vb, bgzf_flags.level, __FUNCLINE);
            break;

        case BGZF_LIBDEFLATE7:
            vb->gzip_compressor = libdeflate_alloc_compressor_1_7 (bgzf_flags.level, vb);
            break;

        case BGZF_ZLIB:
            vb->gzip_compressor = bgzf_alloc (vb, 1, sizeof (z_stream), __FUNCLINE);
            *(z_stream *)vb->gzip_compressor = (z_stream){ .zalloc = bgzf_alloc, .zfree = codec_free_do, .opaque = vb };
            break;

        case BGZF_IGZIP:
            ASSERT (bgzf_flags.level==1 || bgzf_flags.level==2, "igzip: expecting bgzf_flags.level=%u ∈[1,2]", bgzf_flags.level);

            vb->gzip_compressor = bgzf_alloc (vb, 1, (int[]){ 1+ISAL_DEF_LVL0_DEFAULT, ISAL_DEF_LVL1_DEFAULT }[bgzf_flags.level-1], __FUNCLINE); // 1+ to avoid 0
            break;

        default:             
            ABORT ("Invalid bgzf_flags.library=%d", bgzf_flags.library);
    }

    if (flag.show_bgzf)
        iprintf ("BGZF: initialized compressor %s level %u%s\n", 
                 bgzf_library_name (bgzf_flags.library, true), bgzf_flags.level,
                 flag.bgzf == BGZF_BY_ZFILE ? " EXACT" : ""); 
}

static void bgzf_free_compressor (VBlockP vb, FlagsBgzf bgzf_flags)
{
    switch (bgzf_flags.library) {  
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
            ABORT ("Invalid bgzf_flags.library=%d", bgzf_flags.library);
    }

    vb->gzip_compressor = NULL;
}

static void bgzf_compress_one_block (VBlockP vb, rom in, uint32_t isize, 
                                     int32_t block_i, int32_t txt_index, // for show_bgzf (both may be negative - indicating previous VB)
                                     BufferP compressed)
{
    START_TIMER;

    ASSERT0 (vb->gzip_compressor, "vb->gzip_compressor=NULL");

    #define BGZF_MAX_CDATA_SIZE (BGZF_MAX_BLOCK_SIZE - sizeof (BgzfHeader) - sizeof (BgzfFooter))

    buf_alloc (vb, compressed, BGZF_MAX_BLOCK_SIZE, 0, char, 1.2, "scratch");

    BgzfHeader *header = (BgzfHeader *)BAFTc (*compressed);
    buf_add (compressed, BGZF_EOF, sizeof (BgzfHeader)); // template of header - only bsize needs updating

    uint32_t comp_index = compressed->len;
    int out_size;

    if (txt_file->bgzf_flags.library == BGZF_IGZIP) {
        struct isal_zstream strm;
        isal_deflate_stateless_init (&strm);
        strm.gzip_flag      = ISAL_DEFLATE; 
        strm.flush          = NO_FLUSH;
        strm.level          = txt_file->bgzf_flags.level - 1; // note: level 1,2 in bgzf_flags corrsponds to IGZIP level 0,1
        strm.level_buf_size = (int[]){ ISAL_DEF_LVL0_DEFAULT, ISAL_DEF_LVL1_DEFAULT }[strm.level];
        strm.level_buf      = vb->gzip_compressor;
        strm.next_in        = (uint8_t *)in;
        strm.avail_in       = isize;
        strm.next_out       = BAFT8 (*compressed);
        strm.avail_out      = BGZF_MAX_CDATA_SIZE + sizeof (BgzfFooter);
        
        int ret = isal_deflate_stateless (&strm);
        ASSERT (ret == ISAL_DECOMP_OK, "%s: isal_deflate_stateless: %s. isize=%u", VB_NAME, isal_error (ret), isize);

        out_size = BGZF_MAX_CDATA_SIZE + sizeof (BgzfFooter) - strm.avail_out; 
    }

    else if (txt_file->bgzf_flags.library == BGZF_LIBDEFLATE19) { // libdeflate 1.19

        out_size = (int)libdeflate_deflate_compress (vb->gzip_compressor, in, isize, BAFTc (*compressed), BGZF_MAX_CDATA_SIZE);

        // in case the compressed data doesn't fit in one BGZF block, move to compressing at the maximum level. this can
        // happen theoretically (maybe) if the original data was compressed with a higher level, and an uncompressible 64K block was
        // "scratch" to just under 64K while in our compression level it is just over 64K.
        if (!out_size) {
            void *high_compressor = libdeflate_alloc_compressor (vb, 12, __FUNCLINE); // libdefate's highest level
            out_size = libdeflate_deflate_compress (vb->gzip_compressor, in, isize, BAFTc (*compressed), BGZF_MAX_CDATA_SIZE);
            libdeflate_free_compressor (high_compressor, __FUNCLINE);
        }
    }

    else if (txt_file->bgzf_flags.library == BGZF_LIBDEFLATE7) { // libdeflate 1.7

        out_size = (int)libdeflate_deflate_compress_1_7 (vb->gzip_compressor, in, isize, BAFTc (*compressed), BGZF_MAX_CDATA_SIZE);

        // in case the compressed data doesn't fit in one BGZF block, move to compressing at the maximum level. this can
        // happen theoretically (maybe) if the original data was compressed with a higher level, and an uncompressible 64K block was
        // "scratch" to just under 64K while in our compression level it is just over 64K.
        if (!out_size) {
            void *high_compressor = libdeflate_alloc_compressor_1_7 (12, vb); // libdefate's highest level
            out_size = libdeflate_deflate_compress_1_7 (vb->gzip_compressor, in, isize, BAFTc (*compressed), BGZF_MAX_CDATA_SIZE);
            libdeflate_free_compressor_1_7 (high_compressor);
        }
    }

    else { // zlib
        #define strm ((z_stream *)vb->gzip_compressor)

        ASSERT0 (deflateInit2 (vb->gzip_compressor, txt_file->bgzf_flags.level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) == Z_OK, 
                 "deflateInit2 failed");

        strm->next_in   = (uint8_t *)in;
        strm->avail_in  = isize;
        strm->next_out  = BAFT8 (*compressed);
        strm->avail_out = BGZF_MAX_CDATA_SIZE;
        ASSERT (deflate (vb->gzip_compressor, Z_FINISH) == Z_STREAM_END, "deflate failed: msg=%s", strm->msg);

        out_size = BGZF_MAX_CDATA_SIZE - strm->avail_out;
        
        ASSERT0 (deflateEnd (vb->gzip_compressor) == Z_OK, "deflateEnd failed");
        #undef strm
    }

    ASSERT (out_size, "cannot compress block with %u bytes into a BGZF block with %u bytes", isize, BGZF_MAX_BLOCK_SIZE);
    compressed->len += out_size;

    header->bsize = LTEN16 ((uint16_t)(sizeof (BgzfHeader) + out_size + sizeof (BgzfFooter) - 1));

    BgzfFooter footer = { .crc32 = LTEN32 (crc32 (0, in, isize)),
                          .isize = LTEN32 (isize) };
    buf_add (compressed, &footer, sizeof (BgzfFooter));

    if (flag.show_bgzf)
        #define C(i) (i < isize ? char_to_printable (in[i]).s : "")
        iprintf ("%-7s %s i=%d compressed_index=%u size=%u txt_index=%d size=%u txt_data[5]=%1s%1s%1s%1s%1s %s\n",
                threads_am_i_main_thread() ? "MAIN" : threads_am_i_writer_thread() ? "WRITER" : "COMPUTE", VB_NAME, block_i,
                comp_index - (int)sizeof (BgzfHeader), (unsigned)out_size + (int)sizeof (BgzfHeader) + (int)sizeof (BgzfFooter), txt_index, isize, C(0), C(1), C(2), C(3), C(4),
                out_size == BGZF_EOF_LEN ? "EOF" : "");
        #undef C

    COPY_TIMER (bgzf_compress_one_block);
} 

// appends file data to wvb->z_data
void bgzf_write_finalize (void)
{
    // write EOF block if needed
    if (txt_file->bgzf_flags.has_eof_block) {
        buf_add_more (wvb, &wvb->z_data, BGZF_EOF, BGZF_EOF_LEN, "z_data");
    
        if (flag.show_bgzf) iprintf ("%-7s vb=%u   EOF\n", "IO", 0);
    }

    // if we attempted to reconstruct the BGZF block to the original file's bgzf_isizes - warn if we were unlucky and failed
    if (txt_file->bgzf_isizes.len) {
        uint8_t signature[3];
        bgzf_sign (txt_file->disk_so_far + (txt_file->bgzf_flags.has_eof_block ? BGZF_EOF_LEN : 0), signature);

        // commented out because of an inconsistency: here we give an FYI in case we thought we could regenerate the 
        // blocks and failed, however we don't give a warning in the more common case that ZIP couldn't find the BGZF
        // level in the first place
        // ASSERTW (!memcmp (signature, txt_file->bgzf_signature, 3), 
        //          "FYI: %s is recompressed with BGZF (.gz). However, it seems that the original file was compressed with a different compression library than genozip uses, resulting in a slightly different level of compression. Rest assured that the actual data is identical.", 
        //          txt_name);
    }
}

void bgzf_sign (uint64_t disk_size, uint8_t *signature)
{
    signature[0] = (disk_size      ) & 0xff; // LSB of size
    signature[1] = (disk_size >> 8 ) & 0xff;
    signature[2] = (disk_size >> 16) & 0xff;
}

// Entry point of BGZF compression compute thread.
// bgzf-compress vb->txt_data into vb->z_data - using BGZF blocks as prescribed in vb->bgzf_blocks. 
// Note: we hope to reconstruct the exact same byte-level BGZF blocks, as the original files, but that 
// will only happen if the GZIP library (eg libdeflate), version and parameters are the same 
static void bgzf_compress_vb (VBlockP vb)
{
    START_TIMER;

    ASSERTNOTEMPTY (vb->bgzf_blocks);

    buf_alloc (vb, &vb->z_data, 0, vb->bgzf_blocks.len32 * BGZF_MAX_BLOCK_SIZE/2, uint8_t, 1, "z_data"); // alloc based on estimated size
    bgzf_alloc_compressor (vb, txt_file->bgzf_flags);

    for_buf2 (BgzfBlockPiz, block, i, vb->bgzf_blocks) {

        ASSERT (block->txt_index + block->txt_size <= Ltxt, 
                "block=%u out of range: expecting txt_index=%u txt_size=%u <= txt_data.len=%u",
                i, block->txt_index, block->txt_size, Ltxt);

        bgzf_compress_one_block (vb, Btxt (block->txt_index), block->txt_size, i, block->txt_index, &vb->z_data);
    }

    bgzf_free_compressor (vb, txt_file->bgzf_flags);

    vb_set_is_processed (vb); /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ 
    COPY_TIMER (bgzf_compute_thread);
}

static inline uint32_t bgzf_next_isize (void)
{
    #define BGZF_CREATED_BLOCK_SIZE 65280 // same size as observed in htslib-created files

    return txt_file->bgzf_isizes.len ? *B16(txt_file->bgzf_isizes, txt_file->bgzf_isizes.next) + 1 // +1 bc the array values are (isize-1)
                                     : BGZF_CREATED_BLOCK_SIZE; // case: no prescriped isizes 
}

// PIZ: calculate the BGZF blocks within this VB
static uint32_t bgzf_calculate_blocks_one_vb (VBlockP vb, bool is_last)
{
    uint32_t index = 0;

    while ((!txt_file->bgzf_isizes.len || txt_file->bgzf_isizes.next < txt_file->bgzf_isizes.len) && index < Ltxt) { 
        
        uint32_t isize = bgzf_next_isize();
        
        if (index + isize > Ltxt) {
            if (is_last) isize = Ltxt - index; // last BGZF block might be shorter
        else
            break; // the data at the end of this VB doesn't fill a whole BGZF block - pass it down to next vb 
        }

        buf_alloc (vb, &vb->bgzf_blocks, 1, Ltxt / 63000, BgzfBlockPiz, 1.5, "bgzf_blocks");

        BNXT (BgzfBlockPiz, vb->bgzf_blocks) = (BgzfBlockPiz){ .txt_index = index, .txt_size = isize };

        index += isize; 
        txt_file->bgzf_isizes.next++;
    }

    uint32_t remaining = Ltxt - index;
    ASSERT0 (remaining < BGZF_MAX_BLOCK_SIZE, "bgzf_isizes exhausted prematurely"); // if we have 65536 or more remaining, there should have been more isizes

    return remaining;    
}

void bgzf_dispatch_compress (Dispatcher dispatcher, STRp (uncomp), CompIType comp_i, bool is_last)
{
    // uncompressed data to be dealt with by next call to this function (buffer belongs to writer thread)
    static Buffer intercall_txt = {}; // belongs to wvb
    buf_alloc (wvb, &intercall_txt, 0, BGZF_MAX_BLOCK_SIZE, char, 0, "intercall_txt");

    // case: uncomp is not enough to fill a block, just store it to next call
    if (!is_last && (uncomp_len + intercall_txt.len32 < bgzf_next_isize())) {
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

rom bgzf_library_name (BgzfLibraryType library, bool long_name)
{
    return (library < 0 || library >= NUM_BGZF_LIBRARIES) ? "INVALID_BGZF_LIBRARY"
         : long_name ? (rom[])BGZF_LIB_NAMES_LONG[library]
         :             (rom[])BGZF_LIB_NAMES_SHRT[library];
}
