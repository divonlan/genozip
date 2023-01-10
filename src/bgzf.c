// ------------------------------------------------------------------
//   bgzf.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <errno.h>

#include "libdeflate/libdeflate.h"
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

// all data in Little Endian. Defined in https://datatracker.ietf.org/doc/html/rfc1952 and https://samtools.github.io/hts-specs/SAMv1.pdf
typedef struct __attribute__ ((__packed__)) BgzfHeader {
    uint8_t id1;    // Gzip id - must be 31
    uint8_t id2;    // Gzip id - must be 139
    uint8_t cm;     // Compression Method - must be 8
    uint8_t flg;    // Flags - must be 4 (FEXTRA)
    uint32_t mtime; // Modification Time
    uint8_t xfl;    // eXtra Flags
    uint8_t os;     // Operating System
    uint16_t xlen;  // Size of extra fields - 6 if contain only BGZF (may be more)
    uint8_t si1;    // BGZF id - must be 66
    uint8_t si2;    // BGZF id - must be 67
    uint16_t slen;  // BGZF extra field length - must be 2
    uint16_t bsize; // BGZF extra field - (compressed block size -1)
} BgzfHeader;

typedef struct __attribute__ ((__packed__)) BgzfFooter {
    uint32_t crc32; // CRC32 of uncompressed data
    uint32_t isize; // Input (i.e. uncompressed) Size
} BgzfFooter;

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

//--------------------------------------------------------------------
// ZIP SIDE - decompress BGZF-compressed file and prepare BGZF section
//--------------------------------------------------------------------

// ZIP: reads and validates a BGZF block, and returns the uncompressed size or (only if soft_fail) an error
static int32_t bgzf_read_block_raw (FILE *file, // txt_file is not yet assigned when called from file_open_txt_read
                                    uint8_t *block /* must be BGZF_MAX_BLOCK_SIZE in size */, uint32_t *block_size /* out */,
                                    rom basename, bool is_remote, bool soft_fail) 
{
    BgzfHeader *h = (BgzfHeader *)block;

    // read the header
    *block_size = fread (h, 1, sizeof (struct BgzfHeader), file);
    if (! *block_size) return BGZF_ABRUBT_EOF; // EOF without an EOF block

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
    if (!(*block_size == 18 && !memcmp (h, BGZF_PREFIX, BGZF_PREFIX_LEN))) {
        ASSERT (soft_fail, "invalid BGZF block while reading %s", basename);
        return BGZF_BLOCK_GZIP_NOT_BGZIP;
    }

    *block_size = LTEN16 (h->bsize) + 1;
    
    uint32_t body_size = *block_size - sizeof (struct BgzfHeader);
    uint32_t bytes = fread (h+1, 1, body_size, file);

    int save_errno = errno; // we want to report errno of fread, not ftell.

    // if failed, always error, even if soft_fail
    ASSERT (bytes == body_size, "%s Failed to read BGZF block of %s (ftell=%"PRId64" err=\"%s\")", 
            feof (file) ? "Unexpected end of file while reading" : "Failed to read body of", 
            basename, ftello64 (file), 
            (is_remote && save_errno == ESPIPE) ? "Disconnected from remote host" : strerror (save_errno));
    
    return 0; // success
}

// ZIP: reads and validates a BGZF block, and returns the uncompressed size or (only if soft_fail) an error
int32_t bgzf_read_block (File *file, // txt_file is not yet assigned when called from file_open_txt_read
                         uint8_t *block /* must be BGZF_MAX_BLOCK_SIZE in size */, uint32_t *block_size /* out */,
                         bool soft_fail)
{
    int ret = bgzf_read_block_raw ((FILE *)file->file, block, block_size, file->basename, file->is_remote, soft_fail);
    if (ret == BGZF_BLOCK_IS_NOT_GZIP || ret == BGZF_BLOCK_GZIP_NOT_BGZIP) return ret; // happens only if soft_fail
    if (ret == BGZF_ABRUBT_EOF) return 0;

    uint32_t isize_lt32 = *(uint32_t *)&block[*block_size - 4];
    uint32_t isize = LTEN32 (isize_lt32); // 0...65536 per spec
    ASSERT (isize <= 65536, "isize=%u out of range [0,65536]", isize);
    
    // add isize to buffer that will be written to SEC_BGZF
    if (isize) { // don't store EOF block (bc isize=0 cannot be represented as (isize-1) )
        #define BGZF_INITIAL_ALLOC 16 // just of the sake of a bit of effeciency: 16 chosen carefully so 16*63000 < 1MB min vb_size but over segconf size
        if (file->bgzf_isizes.len32 <= BGZF_INITIAL_ALLOC) { // entered thrice: when called from file_open_txt_read, segconf, and in first VB 
            buf_alloc (evb, &file->bgzf_isizes, 0, MAX_(BGZF_INITIAL_ALLOC, segconf.vb_size / 63000), uint16_t, 0, "txt_file->bgzf_isizes");
            buf_alloc (evb, &file->bgzf_starts, 0, MAX_(BGZF_INITIAL_ALLOC, segconf.vb_size / 63000), uint64_t, 0, "txt_file->bgzf_starts");
        }

        buf_append_one (file->bgzf_isizes, BGEN16 ((uint16_t)(isize - 1))); // -1 to make the range 0..65535
        buf_append_one (file->bgzf_starts, txt_file ? txt_file->disk_so_far : 0); // not BGEN bc not written to z_file. note: first block is read from file_open_txt_read before txt_file is assigned
    }
    else 
        if (txt_file) txt_file->bgzf_flags.has_eof_block = true;
    
    return isize;
}

// ZIP: BGZF section per component
void bgzf_compress_bgzf_section (void)
{
    // cases where we don't write the BGZF blocks section
    if (!txt_file->bgzf_isizes.len ||  // this txt file is not compressed with BGZF - we don't need a BGZF section
        txt_file->bgzf_flags.level == BGZF_COMP_LEVEL_UNKNOWN ||  // we don't know the level - so PIZ will reconstruct at default level
        flag.data_modified) return;     // we have data_modified-  the file has changed and we can't reconstruct to the same blocks

    // sanity check
    int64_t total_isize = 0;
    for_buf (uint16_t, isize_p, txt_file->bgzf_isizes)
        total_isize += BGEN16 (*isize_p) + 1; // values 0-65535 correspond to isize 1-65536
    
    ASSERT (total_isize == txt_file->txt_data_so_far_single, "Expecting total_isize=%"PRId64" == txt_file->txt_data_so_far_single=%"PRId64,
            total_isize, txt_file->txt_data_so_far_single);

    // get the best codec for the SEC_BGZF section
    txt_file->bgzf_isizes.len *= sizeof (uint16_t);
    Codec codec = codec_assign_best_codec (evb, NULL, &txt_file->bgzf_isizes, SEC_BGZF);

    evb->comp_i = flag.zip_comp_i; // this goes into SectionEntFileFormat.comp_i via sections_add_to_list
    zfile_compress_section_data_ex (evb, NULL, SEC_BGZF, &txt_file->bgzf_isizes, NULL, 0, codec, (SectionFlags)txt_file->bgzf_flags, NULL);
    txt_file->bgzf_isizes.len /= sizeof (uint16_t); // restore
}

// decompresses a BGZF block in vb->scratch referred to by bb, into its place in vb->txt_data as prescribed by bb
void bgzf_uncompress_one_block (VBlockP vb, BgzfBlockZip *bb)
{
    if (bb->is_decompressed) return; // already decompressed - nothing to do

    ASSERT0 (vb->gzip_compressor, "vb->gzip_compressor=NULL");

    BgzfHeader *h = (BgzfHeader *)Bc (vb->scratch, bb->compressed_index);

    // verify that entire block is within vb->scratch
    ASSERT (bb->compressed_index + sizeof (BgzfHeader) < vb->scratch.len && // we have at least the header - we can access bsize
            bb->compressed_index + (uint32_t)LTEN16 (h->bsize) + 1 <= vb->scratch.len, 
            "%s: BGZF block size goes past the end of in vb->scratch: bb=%s compressed_index=%u vb->scratch.len=%"PRIu64, 
            VB_NAME, display_bb (bb).s, bb->compressed_index, vb->scratch.len);

    ASSERT (h->id1==31 && h->id2==139, "%s: invalid BGZF block in vb->scratch: compressed_index=%u", VB_NAME, bb->compressed_index);

    if (flag.show_bgzf)
        iprintf ("%-7s vb=%s i=%u compressed_index=%u size=%u txt_index=%u size=%u ",
                 threads_am_i_main_thread() ? "MAIN" : "COMPUTE", VB_NAME, 
                 BNUM (vb->bgzf_blocks, bb), bb->compressed_index, bb->comp_size, bb->txt_index, bb->txt_size);

    enum libdeflate_result ret = 
        libdeflate_deflate_decompress (vb->gzip_compressor, 
                                       h+1, bb->comp_size - sizeof(BgzfHeader) - sizeof (BgzfFooter), // compressed
                                       Bc (vb->txt_data, bb->txt_index), bb->txt_size, NULL);  // uncompressed

    ASSERT (ret == LIBDEFLATE_SUCCESS, "libdeflate_deflate_decompress failed: %s", libdeflate_error(ret));

    bb->is_decompressed = true;

    if (flag.show_bgzf)
        #define C(i) ((bb->txt_index + i < vb->txt_data.len) ? char_to_printable (*Bc (vb->txt_data, bb->txt_index + (i))).s : "") 
        iprintf ("txt_data[5]=%1s%1s%1s%1s%1s %s\n", C(0), C(1), C(2), C(3), C(4), bb->comp_size == BGZF_EOF_LEN ? "EOF" : "");
        #undef C
}

// ZIP: called from the compute thread: zip_compress_one_vb and main thread: txtfile_read_block_bgzf
void bgzf_uncompress_vb (VBlockP vb)
{
    START_TIMER;

    vb->gzip_compressor = libdeflate_alloc_decompressor(vb);

    for_buf (BgzfBlockZip, bb, vb->bgzf_blocks)
        bgzf_uncompress_one_block (vb, bb);

    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor);

    buf_free (vb->scratch); // now that we are finished decompressing we can free it

    if (flag.show_time) {
        if (threads_am_i_main_thread ()) COPY_TIMER (bgzf_io_thread)
        else                             COPY_TIMER (bgzf_compute_thread);
    }
}

// ZIP: decompresses a prescribed BGZF block when re-reading DEPN lines
static inline void bgzf_uncompress_one_prescribed_block (VBlockP vb, STRp(bgzf_block), STRc (uncomp_block), uint64_t bb_i)
{
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

    ASSERT (ret == LIBDEFLATE_SUCCESS, "libdeflate_deflate_decompress failed: %s", libdeflate_error(ret));

    if (flag.show_bgzf)
        #define C(i) (i < uncomp_block_len ? char_to_printable (uncomp_block[i]).s : "") 
        iprintf ("txt_data[5]=%1s%1s%1s%1s%1s\n", C(0), C(1), C(2), C(3), C(4));
        #undef C
}

// ZIP: compute thread of a DEPN VB: actually re-reading data into txt_data according to vb->reread_prescription
void bgzf_reread_uncompress_vb_as_prescribed (VBlockP vb, FILE *file)
{
    uint64_t last_offset = -1LL;
    char uncomp_block[BGZF_MAX_BLOCK_SIZE];

    vb->gzip_compressor = libdeflate_alloc_decompressor(vb);

    for_buf (RereadLine, line, vb->reread_prescription) {
        
        // a line might span 1 or more BGZF blocks
        while (line->line_len) { 

            uint64_t offset = *B64 (txt_file->bgzf_starts, line->offset.bb_i);
            uint16_t isize  = BGEN16 (*B16 (txt_file->bgzf_isizes, line->offset.bb_i)) + 1;

            if (offset != last_offset) {
                ASSERT (!fseeko64 (file, offset, SEEK_SET),
                        "%s: fseeko64 on %s failed while rereading BGZF depn lines: %s", VB_NAME, txt_file->name, strerror(errno));

                STRl (bgzf_block, BGZF_MAX_BLOCK_SIZE);
                bgzf_read_block_raw (file, (uint8_t*)qSTRa(bgzf_block), txt_file->basename, false, false);
            
                bgzf_uncompress_one_prescribed_block (vb, STRa(bgzf_block), uncomp_block, isize, line->offset.bb_i);
            
                last_offset = offset;
            }

            uint32_t subline_len = MIN_(line->line_len, isize - line->offset.uoffset);
            memcpy (BAFTtxt, &uncomp_block[line->offset.uoffset], subline_len);
            vb->txt_data.len32 += subline_len;
            
            // if this line continues to next BGZF block - it starts from the beginning of that block, its remainder is subline_len shorter
            line->line_len -= subline_len;
            line->offset.bb_i++;
            line->offset.uoffset = 0;      
        }
    }

    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor);
}

static void *bgzf_alloc (void *vb_, unsigned items, unsigned size, FUNCLINE)
{
    return codec_alloc_do ((VBlockP )vb_, (uint64_t)items * (uint64_t)size, 1, func, code_line); // all bzlib buffers are constant in size between subsequent compressions
}

void bgzf_libdeflate_initialize (void)
{
    libdeflate_set_memory_allocator (bgzf_alloc, codec_free_do);
}

// ZIP: tests a BGZF block against libdeflate's level 0-12
// returns the level 0-12 if detected, or BGZF_COMP_LEVEL_UNKNOWN if not
// NOTE: even if we detected the level based on the first block of the file - it's not certain that the file was compressed with this 
// level. This is because multiple levels might result in the same compression for this block, but maybe not for other blocks in the file.
struct FlagsBgzf bgzf_get_compression_level (rom filename, bytes comp_block, uint32_t comp_block_size, uint32_t uncomp_block_size)
{
    static const struct { BgzfLibraryType library; int level; } levels[] =  // test in the order of likelihood of observing them
        { {BGZF_LIBDEFLATE,6}, {BGZF_ZLIB,6},  {BGZF_ZLIB,4},  {BGZF_LIBDEFLATE,9},  {BGZF_LIBDEFLATE,8}, {BGZF_LIBDEFLATE,7}, 
          {BGZF_LIBDEFLATE,5}, {BGZF_LIBDEFLATE,4}, {BGZF_LIBDEFLATE,3}, {BGZF_LIBDEFLATE,2}, {BGZF_LIBDEFLATE,1}, 
          {BGZF_LIBDEFLATE,0}, {BGZF_LIBDEFLATE,12}, {BGZF_LIBDEFLATE,11}, {BGZF_LIBDEFLATE,10}, {BGZF_ZLIB,9}, {BGZF_ZLIB,7}, 
          {BGZF_ZLIB,8}, {BGZF_ZLIB,5}, {BGZF_ZLIB,3}, {BGZF_ZLIB,2}, {BGZF_ZLIB,1} };

    // ignore the header and footer of the block
    comp_block      += sizeof (BgzfHeader);
    comp_block_size -= sizeof (BgzfHeader) + sizeof (BgzfFooter);

    // decompress block
    if (!evb->gzip_compressor) evb->gzip_compressor = libdeflate_alloc_decompressor(evb);

    uint8_t uncomp_block[uncomp_block_size]; // at most 64K - fits on stack
    enum libdeflate_result ret = 
        libdeflate_deflate_decompress (evb->gzip_compressor, comp_block, comp_block_size, uncomp_block, uncomp_block_size, NULL);

    ASSERT (ret == LIBDEFLATE_SUCCESS, "unable to read file %s. It appears to be compressed with BGZF, however decompression failed: %s", filename, libdeflate_error(ret));
    
    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&evb->gzip_compressor);

    // now, re-compress with each level until the compressed size matches
    uint32_t recomp_size;
    uint8_t recomp_block[BGZF_MAX_BLOCK_SIZE]; 

    for (int level_i=0; level_i < ARRAY_LEN (levels); level_i++) { 
        struct FlagsBgzf l = { .library = levels[level_i].library, .level = levels[level_i].level };

        if (l.library == BGZF_LIBDEFLATE) {
            void *compressor = libdeflate_alloc_compressor (l.level, evb);
            recomp_size = (uint32_t)libdeflate_deflate_compress (compressor, uncomp_block, uncomp_block_size, recomp_block, BGZF_MAX_BLOCK_SIZE);

            libdeflate_free_compressor (compressor);
        }
        else { // BGZF_ZLIB
            z_stream strm = { .zalloc = bgzf_alloc, .zfree  = codec_free_do, .opaque = evb };
            // deflateInit2 with the default zlib parameters, with is also the same as htslib does
            ASSERT0 (deflateInit2 (&strm, l.level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) == Z_OK, "deflateInit2 failed");

            strm.next_in   = uncomp_block;
            strm.avail_in  = uncomp_block_size;
            strm.next_out  = recomp_block;
            strm.avail_out = sizeof (recomp_block);
            ASSERT (deflate (&strm, Z_FINISH) == Z_STREAM_END, "deflate failed: msg=%s", strm.msg);

            recomp_size = sizeof (recomp_block) - strm.avail_out;
            
            ASSERT0 (deflateEnd (&strm) == Z_OK, "deflateEnd failed");
        }

        bool identical = recomp_size == comp_block_size && !memcmp (comp_block, recomp_block, comp_block_size);

        if (flag.show_bgzf) 
            iprintf ("Testing library %s level %u: size_in_file=%u size_in_test=%u identical=%s\n", 
                     l.library == BGZF_LIBDEFLATE ? "libdeflate" : "zlib", l.level, comp_block_size, recomp_size, identical ? "Yes" : "No");

        if (identical) {
            if (flag.show_bgzf) 
                iprintf ("File %s: Identified as compressed with %s level %u\n", 
                         filename, l.library == BGZF_LIBDEFLATE ? "libdeflate" : "zlib", l.level);
            return l; // this still might be wrong, see comment in function header ^
        }       
    }

    if (flag.show_bgzf) iprintf ("File %s: Could not identify compression library and level\n", filename);
    return (struct FlagsBgzf){ .level = BGZF_COMP_LEVEL_UNKNOWN };
}

// ZIP: called by Seg to set the bgzf index of the next line
void bgzf_zip_advance_index (VBlockP vb, uint32_t line_len)
{
    vb->line_bgzf_uoffset += line_len;

    // udpate current_bb_i and bgzf_offset (note: line_len might span multiple bgzf blocks)
    BgzfBlockZip *bb;
    for (bb = B(BgzfBlockZip, vb->bgzf_blocks, vb->bgzf_blocks.current_bb_i); vb->line_bgzf_uoffset >= bb->txt_size; bb++) 
        vb->line_bgzf_uoffset -= bb->txt_size; // index into the next BGZF block

    vb->bgzf_blocks.current_bb_i = BNUM(vb->bgzf_blocks, bb);
}

// ZIP: after reading data for a txt_header or VB, copy unconsumed bgzf_blocks to txt_file->unconsumed_bgzf_blocks
// The first block might be partially consumed.
int64_t bgzf_copy_unconsumed_blocks (VBlockP vb)
{
    if (!vb->bgzf_blocks.len) return 0; // not a BGZF-compressed file

    int32_t consumed = vb->txt_data.len32 +   // amount of data consumed by this VB
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
                        B(BgzfBlockZip, vb->bgzf_blocks, i), vb->bgzf_blocks.len32 - i, NULL);

            txt_file->unconsumed_bgzf_blocks.consumed_by_prev_vb = consumed; // part of first BGZF block already consumed
            done = true;
        }
        else if (!done)
            compressed_size += bb[i].comp_size;

        consumed -= bb[i].txt_size;
    }

    // sanity check
    ASSERT (-consumed == txt_file->unconsumed_txt.len32, "Expecting (-consumed)=%d == unconsumed_txt.len=%u", -consumed, txt_file->unconsumed_txt.len32);

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

    ASSERT (available >= vb->txt_data.len32, "BGZF blocks in txt_file->unconsumed_bgzf_blocks cover only %d bytes, less than the needed unconsumed_bytes=%d", 
            available, vb->txt_data.len32);
}

//---------
// PIZ SIDE
//---------

bool bgzf_load_isizes (Section sl_ent) 
{
    // skip passed all VBs of this component, and read SEC_BGZF - the last section of the component - if it exists
    // but stop if we encounter the next component's SEC_TXT_HEADER before seeing the BGZF
    if (!sections_next_sec2 (&sl_ent, SEC_BGZF, SEC_TXT_HEADER) // updates sl_ent, but its a local var, doesn't affect our caller
        || sl_ent->st == SEC_TXT_HEADER)
        return false; // this component doesn't contain a BGZF section

    int32_t offset = zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", SEC_BGZF, sl_ent);

    SectionHeaderP header = (SectionHeaderP)Bc (evb->z_data, offset);
    txt_file->bgzf_flags = header->flags.bgzf;

    // if we don't know the compression level, or if original file had compression level 0 (no compression), go with the default
    if (txt_file->bgzf_flags.level == BGZF_COMP_LEVEL_UNKNOWN ||
        !txt_file->bgzf_flags.level) { // the user can override this behaviour with --bgzf
        txt_file->bgzf_flags.level   = BGZF_COMP_LEVEL_DEFAULT;
        txt_file->bgzf_flags.library = BGZF_LIBDEFLATE;
    }

    zfile_uncompress_section (evb, header, &txt_file->bgzf_isizes, "txt_file->bgzf_isizes", 0, SEC_BGZF);
    txt_file->bgzf_isizes.len /= 2;

    // convert to native endianity from big endian
    for_buf (uint16_t, isize_p, txt_file->bgzf_isizes)
        *isize_p = BGEN16 (*isize_p); // now it contains isize-1 - a value 0->65535 representing an isize 1->65536

    return true; // bgzf_isizes successfully loaded
}                

static void bgzf_alloc_compressor (VBlockP vb, struct FlagsBgzf bgzf_flags)
{
    ASSERT0 (!vb->gzip_compressor, "expecting vb->gzip_compressor=NULL");

    if (bgzf_flags.library == BGZF_LIBDEFLATE)  // libdeflate
        vb->gzip_compressor = libdeflate_alloc_compressor (bgzf_flags.level, vb);

    else { // zlib
        vb->gzip_compressor = bgzf_alloc (vb, 1, sizeof (z_stream), __FUNCLINE);
        *(z_stream *)vb->gzip_compressor = (z_stream){ .zalloc = bgzf_alloc, .zfree  = codec_free_do, .opaque = vb };
    }
}

static void bgzf_free_compressor (VBlockP vb, struct FlagsBgzf bgzf_flags)
{
    if (bgzf_flags.library == BGZF_LIBDEFLATE)  // libdeflate
        libdeflate_free_compressor (vb->gzip_compressor);
    else
        codec_free (vb, vb->gzip_compressor);

    vb->gzip_compressor = NULL;
}

static uint32_t bgzf_compress_one_block (VBlockP vb, rom in, uint32_t isize, 
                                         int32_t block_i, int32_t txt_index, // for show_bgzf (both may be negative - indicating previous VB)
                                         BufferP compressed)
{
    ASSERT0 (vb->gzip_compressor, "vb->gzip_compressor=NULL");

    #define BGZF_MAX_CDATA_SIZE (BGZF_MAX_BLOCK_SIZE - sizeof (BgzfHeader) - sizeof (BgzfFooter))

    buf_alloc (vb, compressed, BGZF_MAX_BLOCK_SIZE, 0, char, 1.2, "scratch");

    BgzfHeader *header = (BgzfHeader *)BAFTc (*compressed);
    buf_add (compressed, BGZF_EOF, sizeof (BgzfHeader)); // template of header - only bsize needs updating

    uint32_t comp_index = compressed->len;
    int out_size;

    if (txt_file->bgzf_flags.library == BGZF_LIBDEFLATE) { // libdeflate

        out_size = (int)libdeflate_deflate_compress (vb->gzip_compressor, in, isize, BAFTc (*compressed), BGZF_MAX_CDATA_SIZE);

        // in case the compressed data doesn't fit in one BGZF block, move to compressing at the maximum level. this can
        // happen theoretically (maybe) if the original data was compressed with a higher level, and an uncompressible 64K block was
        // "scratch" to just under 64K while in our compression level it is just over 64K.
        if (!out_size) {
            void *high_compressor = libdeflate_alloc_compressor (12, vb); // libdefate's highest level
            out_size = libdeflate_deflate_compress (vb->gzip_compressor, in, isize, BAFTc (*compressed), BGZF_MAX_CDATA_SIZE);
            libdeflate_free_compressor (high_compressor);
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

    if (flag.show_bgzf)
        #define C(i) (i < isize ? char_to_printable (in[i]).s : "")
        iprintf ("%-7s vb=%s i=%d compressed_index=%u size=%u txt_index=%d size=%u txt_data[5]=%1s%1s%1s%1s%1s %s\n",
                threads_am_i_main_thread() ? "MAIN" : threads_am_i_writer_thread() ? "WRITER" : "COMPUTE", VB_NAME, block_i,
                comp_index, (unsigned)out_size, txt_index, isize, C(0), C(1), C(2), C(3), C(4),
                out_size == BGZF_EOF_LEN ? "EOF" : "");
        #undef C

    ASSERT (out_size, "cannot compress block with %u bytes into a BGZF block with %u bytes", isize, BGZF_MAX_BLOCK_SIZE);
    compressed->len += out_size;

    header->bsize = LTEN16 ((uint16_t)(sizeof (BgzfHeader) + out_size + sizeof (BgzfFooter) - 1));

    BgzfFooter footer = { .crc32 = LTEN32 (crc32 (0, in, isize)),
                          .isize = LTEN32 (isize) };
    buf_add (compressed, &footer, sizeof (BgzfFooter));

    return (uint32_t)out_size;
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
        //          txt_file->name);
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

        ASSERT (block->txt_index + block->txt_size <= vb->txt_data.len32, 
                "block=%u out of range: expecting txt_index=%u txt_size=%u <= txt_data.len=%u",
                i, block->txt_index, block->txt_size, vb->txt_data.len32);

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

    while ((!txt_file->bgzf_isizes.len || txt_file->bgzf_isizes.next < txt_file->bgzf_isizes.len) && index < vb->txt_data.len32) { 
        
        uint32_t isize = bgzf_next_isize();
        
        if (index + isize > vb->txt_data.len32) {
            if (is_last) isize = vb->txt_data.len32 - index; // last BGZF block might be shorter
        else
            break; // the data at the end of this VB doesn't fill a whole BGZF block - pass it down to next vb 
        }

        buf_alloc (vb, &vb->bgzf_blocks, 1, vb->txt_data.len32 / 63000, BgzfBlockPiz, 1.5, "bgzf_blocks");

        BNXT (BgzfBlockPiz, vb->bgzf_blocks) = (BgzfBlockPiz){ .txt_index = index, .txt_size = isize };

        index += isize; 
        txt_file->bgzf_isizes.next++;
    }

    uint32_t remaining = vb->txt_data.len32 - index;
    ASSERT0 (remaining < BGZF_MAX_BLOCK_SIZE, "bgzf_isizes exhausted prematurely"); // if we have 65536 or more remaining, there should have been more isizes

    return remaining;    
}

void bgzf_dispatch_compress (Dispatcher dispatcher, STRp (uncomp), bool is_last)
{
    // uncompressed data to be dealt with by next call to this function (buffer belongs to writer thread)
    static Buffer intercall_txt = EMPTY_BUFFER; // belongs to wvb
    buf_alloc (wvb, &intercall_txt, 0, BGZF_MAX_BLOCK_SIZE, char, 0, "intercall_txt");

    // case: uncomp is not enough to fill a block, just store it to next call
    if (!is_last && (uncomp_len + intercall_txt.len32 < bgzf_next_isize())) {
        memcpy (BAFTc(intercall_txt), uncomp, uncomp_len);
        intercall_txt.len32 += uncomp_len;
        return;
    }

    if (uncomp_len || intercall_txt.len) { // might be 0 if is_last, in some cases

        VBlockP vb = dispatcher_generate_next_vb (dispatcher, wvb->vblock_i, COMP_NONE);

        // build uncompressed data for this VB - some data left over from previous VB + data from wvb
        buf_alloc_exact (vb, vb->txt_data, intercall_txt.len + uncomp_len, char, "txt_data");
        if (intercall_txt.len32) memcpy (B1STtxt, intercall_txt.data, intercall_txt.len32);
        memcpy (Btxt(intercall_txt.len32), uncomp, uncomp_len);

        // calculate BGZF blocks - and trim data that doesn't fill a block - to be moved to next VB
        if ((intercall_txt.len32 = bgzf_calculate_blocks_one_vb (vb, is_last))) {
            vb->txt_data.len32 -= intercall_txt.len32;
            memcpy (B1STc(intercall_txt), BAFTtxt, intercall_txt.len32);
        }

        // BGZF-compress vb->txt_data in a separate thread
        dispatcher_compute (dispatcher, bgzf_compress_vb);
    }

    if (is_last) 
        dispatcher_set_no_data_available (dispatcher, false, DATA_EXHAUSTED);
}

rom bgzf_library_name (BgzfLibraryType library)
{
    rom names[] = BGZF_LIB_NAMES;

    return (library >= 0 && library < NUM_BGZF_LIBRARIES) ? names[library] : "INVALID_BGZF_LIBRARY";
}
