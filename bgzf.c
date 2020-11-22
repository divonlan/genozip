#include <errno.h>

#include <libdeflate.h>
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

// all data in Little Endian
typedef struct __attribute__ ((__packed__)) BgzfHeader {
    uint8_t id1, id2, cm, flg;
    uint32_t mtime;
    uint8_t xfl, os;
    uint16_t xlen;
    uint8_t si1, si2;
    uint16_t slen, bsize;
} BgzfHeader;

typedef struct __attribute__ ((__packed__)) BgzfFooter {
    uint32_t crc32;
    uint32_t isize;
} BgzfFooter;

// possible return values, see libdeflate_result in libdeflate.h
static const char *libdeflate_error (int err)
{
    switch (err) {
        case LIBDEFLATE_SUCCESS            : return "SUCCESS";
        case LIBDEFLATE_BAD_DATA           : return "BAD DATA";
        case LIBDEFLATE_SHORT_OUTPUT       : return "SHORT OUTPUT";
        case LIBDEFLATE_INSUFFICIENT_SPACE : return "INSUFFICIENT SPACE";
        default                            : return "Undefined libdeflate error";
    }
}

//--------------------------------------------------------------------
// ZIP SIDE - decompress BGZF-compressed file and prepare BGZF section
//--------------------------------------------------------------------

// unfortunately, we can't use fread(), since file_open_txt_read, who calls us, might later try to open the file
// with gzdopen which then gzread uses read() - causing skipping of the data that already is read into fread() buffers
// here, we read the amount of data requested, unless EOF is encountered first
static uint32_t bgzf_fread (File *file, void *dst_buf, uint32_t count)
{
    uint32_t bytes=0;
    while (bytes < count) {
        uint32_t bytes_once = read (fileno((FILE *)file->file), &((char*)dst_buf)[bytes], count - bytes);
        if (!bytes_once) return bytes;  // EOF

        bytes += bytes_once;
    }

    return bytes; // == count
}

// reads and validates a BGZF block, and returns the uncompressed size, or, in case of soft_fail:
int32_t bgzf_read_block (File *file, // txt_file is not yet assigned when called from file_open_txt_read
                         uint8_t *block /* must be BGZF_MAX_BLOCK_SIZE in size */, uint32_t *block_size /* out */,
                         bool soft_fail)
{
    BgzfHeader *h = (BgzfHeader *)block;

    *block_size = bgzf_fread (file, h, sizeof (struct BgzfHeader)); // read the header
    if (! *block_size) return 0; // EOF without an EOF block

    if (*block_size < 12) {
        ASSERT0 (soft_fail, "Error in bgzf_read_block: block_size read is < 12"); // less than the minimal gz block header size
        return BGZF_BLOCK_IS_NOT_GZIP;
    } 

    // case: this is not a GZ / BGZF block at all (see: https://tools.ietf.org/html/rfc1952)
    if (h->id1 != 31 || h->id2 != 139) {
        ASSERT (soft_fail, "Error: expecting %s to be compressed with gzip format, but it is not", file->name);
        return BGZF_BLOCK_IS_NOT_GZIP;
    }

    // case: this is GZIP block that is NOT a valid BGZF block (see: https://samtools.github.io/hts-specs/SAMv1.pdf)
    if (!(*block_size == 18 && h->cm==8 && h->flg==4 && h->si1==66 && h->si2==67)) {
        ASSERT (soft_fail, "Error in bgzf_read_block: invalid BGZF block while reading %s", file->name);
        return BGZF_BLOCK_GZIP_NOT_BGZIP;
    }

    *block_size = LTEN16 (h->bsize) + 1;
    
    uint32_t body_size = *block_size - sizeof (struct BgzfHeader);
    ASSERT (bgzf_fread (file, h+1, body_size) == body_size, 
            "Error in bgzf_read_block: failed to read body of BGZF block in %s: %s", file->name, strerror (errno));

    uint32_t isize_lt32 = *(uint32_t *)&block[*block_size - 4];
    uint32_t isize = LTEN32 (isize_lt32); // 0...65536 per spec

    // add isize to buffer that will be written to SEC_BGZF
    if (isize) { // don't store EOF block (with isize=0)
        buf_alloc_more (evb, &file->bgzf_isizes, 1, global_max_memory_per_vb / 63000, uint16_t, 2, "bgzf_isizes");
        NEXTENT (uint16_t, file->bgzf_isizes) = BGEN16 ((uint16_t)(isize - 1)); // -1 to make the range 0..65535
    }
    
    return isize;
}

void bgzf_compress_bgzf_section (void)
{
    if (!txt_file->bgzf_isizes.len) return; // this txt file is not compressed with BGZF - we don't need a BGZF section

    txt_file->bgzf_isizes.len *= sizeof (uint16_t);

    // usually BZ2 is the best, but if the sizes are highly random (happens in long read SAM/BAM), then BZ2 can be bloated
    // and even produce an error. that's why we test.
    Codec codec = codec_assign_best_codec (evb, NULL, &txt_file->bgzf_isizes, SEC_BGZF);

    zfile_compress_section_data_codec (evb, SEC_BGZF, &txt_file->bgzf_isizes, NULL, 0, codec);

    txt_file->bgzf_isizes.len /= sizeof (uint16_t); // restore
}

void bgzf_uncompress_one_block (VBlock *vb, BgzfBlockZip *bb)
{
    if (bb->is_decompressed) return; // already decompressed - nothing to do

    BgzfHeader *h = (BgzfHeader *)ENT (char, vb->compressed, bb->compressed_index);

    // verify that entire block is within vb->compressed
    ASSERT (bb->compressed_index + sizeof (BgzfHeader) < vb->compressed.len && // we have at least the header - we can access bsize
            bb->compressed_index + (uint32_t)LTEN16 (h->bsize) + 1 <= vb->compressed.len, 
            "Error in bgzf_uncompress_vb: bgzf block size goes past the end of in vb->compressed: vb=%u compressed_index=%u vb->compressed.len=%"PRIu64, 
            vb->vblock_i, bb->compressed_index, vb->compressed.len);

    ASSERT (h->id1==31 && h->id2==139, "Error in bgzf_uncompress_vb: not a valid bgzf block in vb->compressed: vb=%u compressed_index=%u", vb->vblock_i, bb->compressed_index);

    if (!vb->libdeflate) vb->libdeflate = libdeflate_alloc_decompressor();

    if (flag.show_bgzf)
        fprintf (stderr, "%-7s vb=%u i=%u compressed_index=%u size=%u txt_index=%u size=%u ",
                 arch_am_i_io_thread() ? "IO" : "COMPUTE", vb->vblock_i, 
                 ENTNUM (vb->bgzf_blocks, bb), bb->compressed_index, bb->comp_size, bb->txt_index, bb->txt_size);

    enum libdeflate_result ret = 
        libdeflate_deflate_decompress (vb->libdeflate, 
                                       h+1, bb->comp_size - sizeof(*h),   // compressed
                                       ENT (char, vb->txt_data, bb->txt_index), bb->txt_size, NULL); // uncompressed

    ASSERT (ret == LIBDEFLATE_SUCCESS, "Error in bgzf_uncompress_vb: libdeflate_deflate_decompress failed: %s", libdeflate_error(ret));

    bb->is_decompressed = true;

    if (flag.show_bgzf)
        #define C(i) ((bb->txt_index + i < vb->txt_data.len) ? char_to_printable (*ENT (char, vb->txt_data, bb->txt_index + (i))) : 0) 
        fprintf (stderr, "txt_data[5]=%c%c%c%c%c\n", C(0), C(1), C(2), C(3), C(4));
}

// ZIP: called from the compute thread: zip_compress_one_vb and I/O thread: txtfile_read_block_bgzf
void bgzf_uncompress_vb (VBlock *vb)
{
    START_TIMER;

    for (uint32_t block_i=0; block_i < vb->bgzf_blocks.len ; block_i++) {
        BgzfBlockZip *bb = ENT (BgzfBlockZip, vb->bgzf_blocks, block_i);
        bgzf_uncompress_one_block (vb, bb);
    }

    buf_free (&vb->compressed); // now that we are finished decompressing we can free it

    if (flag.show_time) {
        if (arch_am_i_io_thread ()) COPY_TIMER (bgzf_io_thread)
        else                        COPY_TIMER (bgzf_compute_thread);
    }
}

//---------
// PIZ SIDE
//---------

void bgzf_read_and_uncompress_isizes (const SectionListEntry *sl_ent) 
{
    // skip passed all VBs of this component, and read SEC_BGZF - the last section of the component - if it exists
    ASSERT0 (sections_get_next_section_of_type (&sl_ent, SEC_BGZF, false, false), // updates sl_ent, but its a local var, doesn't affect our caller
             "Error in bgzf_read_and_uncompress_isizes: file does not contain a SEC_BGZF section");

    int32_t offset = zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", SEC_BGZF, sl_ent);

    zfile_uncompress_section (evb, ENT (char, evb->z_data, offset), &txt_file->bgzf_isizes, "txt_file->bgzf_isizes", 0, SEC_BGZF);
    txt_file->bgzf_isizes.len /= 2;

    // convert to native endianity from big endian
    ARRAY (uint16_t, isizes, txt_file->bgzf_isizes);
    for (uint32_t i=0; i < txt_file->bgzf_isizes.len; i++) 
        isizes[i] = BGEN16 (isizes[i]); // now it contains isize-1 - a value 0->65535 representing an isize 1->65536
}                

// I/O thread ahead of dispatching - calculate the BGZF blocks within this VB that need to be compressed by
// the compute thread - i.e. excluding the flanking regions of txt_data that might share a BGZF block with the adjacent VB
// and will be compressed by bgzf_compress_and_write_split_blocks.
void bgzf_calculate_blocks_one_vb (VBlock *vb, uint32_t vb_txt_data_len)
{
    ARRAY (uint16_t, isizes, txt_file->bgzf_isizes);
    #define next_isize txt_file->bgzf_isizes.param // we use bgzf_isizes.param as "next" - iterator on bgzf_isizes
    
    // skip the initial region of txt_data that has a split block with the previous VB - 
    // the previous VBs final data is now in txt_file->unconsumed_txt
    uint32_t index = 0;
    if (next_isize && txt_file->unconsumed_txt.len)
        index = (uint32_t)(isizes[next_isize-1] + 1) - (uint32_t)txt_file->unconsumed_txt.len;

    while (next_isize < txt_file->bgzf_isizes.len) { 
        
        uint32_t isize = (uint32_t)isizes[next_isize] + 1;// +1 bc the array values are (isize-1)

        if (index + isize > vb_txt_data_len) break; // this block doesn't fit in the VB

        buf_alloc_more (vb, &vb->bgzf_blocks, 1, global_max_memory_per_vb / 63000, BgzfBlockPiz, 1.5, "bgzf_blocks");

        NEXTENT (BgzfBlockPiz, vb->bgzf_blocks) = (BgzfBlockPiz){ .txt_index = index, .txt_size = isize };

        index += isize; 
        next_isize++;
    }

    #undef next_isize
}

static uint32_t bgzf_compress_one_block (VBlock *vb, const char *in, uint32_t isize)
{
    START_TIMER;

    // use level 6 - same as samtools and bgzip default level. if we're lucky (same gzip library and user used default level),
    // we will reconstruct precisely even at the .gz level
    #define BGZF_DEFAULT_COMPRESSION_LEVEL 6 
    
    #define BGZF_MAX_CDATA_SIZE (BGZF_MAX_BLOCK_SIZE - sizeof (BgzfHeader) - sizeof (BgzfFooter))

    if (!vb->libdeflate) vb->libdeflate = libdeflate_alloc_compressor (BGZF_DEFAULT_COMPRESSION_LEVEL);

    if (flag.show_bgzf)
        fprintf (stderr, "%-7s vb=%u", arch_am_i_io_thread() ? "IO" : "COMPUTE", vb->vblock_i);

    buf_alloc_more (vb, &vb->compressed, BGZF_MAX_BLOCK_SIZE, 0, char, 1.2, "compressed");

    BgzfHeader *header = (BgzfHeader *)AFTERENT (char, vb->compressed);
    buf_add (&vb->compressed, BGZF_EOF, sizeof (BgzfHeader)); // template of header - only bsize needs updating

    size_t out_size = libdeflate_deflate_compress (vb->libdeflate, in, isize, AFTERENT (char, vb->compressed), BGZF_MAX_CDATA_SIZE);

    // in case the compressed data doesn't fit in one BGZF block, move to compressing at the maximum level. this can
    // happen theoretically (maybe) if the original data was compressed with a higher level, and an uncompressible 64K block was
    // "compressed" to just under 64K while in our compression level it is just over 64K.
    if (!out_size) {
        void *high_compressor = libdeflate_alloc_compressor (12); // libdefate's highest level
        out_size = libdeflate_deflate_compress (vb->libdeflate, in, isize, AFTERENT (char, vb->compressed), BGZF_MAX_CDATA_SIZE);
        libdeflate_free_compressor (high_compressor);
    }

    ASSERT (out_size, "Error in bgzf_uncompress_vb: cannot compress block with %u bytes into a BGZF block with %u bytes", isize, BGZF_MAX_BLOCK_SIZE);
    vb->compressed.len += out_size;

    header->bsize = LTEN16 ((uint16_t)(sizeof (BgzfHeader) + out_size + sizeof (BgzfFooter) - 1));

    BgzfFooter footer = { .crc32 = LTEN32 (libdeflate_crc32 (0, in, isize)),
                          .isize = LTEN32 (isize) };
    buf_add (&vb->compressed, &footer, sizeof (BgzfFooter));

    if (flag.show_time) {
        if (arch_am_i_io_thread ()) COPY_TIMER (bgzf_io_thread)
        else                        COPY_TIMER (bgzf_compute_thread);
    }

    return (uint32_t)out_size;
} 

// Called in the Compute Thread for VBs and in I/O thread for the Txt Header.
// bgzf-compress vb->txt_data into vb->compressed - reconstructing the same-isize BGZF blocks as the original txt file
// note: data at the beginning and end of txt_data that doesn't fit into a whole BGZF block (i.e. the block is shared
// with the adjacent VB) - we don't compress here, but rather in bgzf_compress_flanking().
// Note: we hope to reconstruct the exact same BGZF blocks, but that will only happen if the GZIP library, version and
// parameters are the same 
void bgzf_compress_vb (VBlock *vb)
{
    ASSERT0 (!vb->compressed.len, "Error in bgzf_compress_vb: expecting vb->compressed to be free, but its not");

    buf_alloc (vb, &vb->compressed, vb->bgzf_blocks.len * BGZF_MAX_BLOCK_SIZE/2, 1, "compressed"); // alloc based on estimated size

    ARRAY (BgzfBlockPiz, blocks, vb->bgzf_blocks);
    for (uint32_t i=0; i < vb->bgzf_blocks.len; i++) {

        ASSERT (blocks[i].txt_index + blocks[i].txt_size <= vb->txt_data.len, 
                "Error in bgzf_compress_vb: block=%u out of range: expecting txt_index=%u txt_size=%u <= txt_data.len=%u",
                i, blocks[i].txt_index, blocks[i].txt_size, (uint32_t)vb->txt_data.len);

        bgzf_compress_one_block (vb, ENT (char, vb->txt_data, blocks[i].txt_index), blocks[i].txt_size);
    }
}

// Called by I/O thread to complete the work Compute Thread cannot do - see 1,2,3 below
void bgzf_write_to_disk (VBlock *vb)
{
    if (!txt_file->unconsumed_txt.len && !vb->compressed.len) return; // nothing to do

    ASSERT (vb->bgzf_blocks.len, "Error in bgzf_write_to_disk: vb->bgzf_blocks.len=0 in vb=%u", vb->vblock_i);
    
    // get the regions at the beginning and end of txt_data which the compute thread did NOT compress
    uint32_t uncompressed_at_vb_start = FIRSTENT (BgzfBlockPiz, vb->bgzf_blocks)->txt_index;
    BgzfBlockPiz *last_block = LASTENT (BgzfBlockPiz, vb->bgzf_blocks);
    uint32_t uncompressed_index_at_vb_end = last_block->txt_index + last_block->txt_size;

    // Note: a BGZF block can be split by maximum between two VBs. In this, we are counting on VBs (apart for the last one)
    // being long enough to not allow a 64K bgzf block to span 3 VBs.

    // 1. bgzf-compress the BGZF block that is split between end the previous VB (data currently in txt_file->unconsumed_txt)
    //    and the beginning of this VB, and write the compressed data to disk
    if (txt_file->unconsumed_txt.len) {
    
        ASSERT (txt_file->unconsumed_txt.len + uncompressed_at_vb_start <= BGZF_MAX_BLOCK_SIZE, 
                "Error in bgzf_write_to_disk: too much uncompressed data at start: txt_file->unconsumed_txt.len=%u + uncompressed_at_vb_start=%u > %u",
                (uint32_t)txt_file->unconsumed_txt.len, uncompressed_at_vb_start, BGZF_MAX_BLOCK_SIZE);

        char block[BGZF_MAX_BLOCK_SIZE];
        uint32_t block_len = txt_file->unconsumed_txt.len + uncompressed_at_vb_start;
        memcpy (block, txt_file->unconsumed_txt.data, txt_file->unconsumed_txt.len);
        memcpy (&block[txt_file->unconsumed_txt.len], FIRSTENT (char, vb->txt_data), uncompressed_at_vb_start);

        bgzf_compress_one_block (evb, block, block_len);

        if (!flag.test) file_write (txt_file, evb->compressed.data, evb->compressed.len);

        txt_file->txt_data_so_far_single += block_len;
        txt_file->disk_so_far            += evb->compressed.len;

        buf_free (&txt_file->unconsumed_txt);
        buf_free (&evb->compressed);
    }

    // 2. Write all the "middle" blocks of this VB, compressed by bgzf_compress_vb, currently in vb->compressed, to disk
    if (!flag.test) file_write (txt_file, vb->compressed.data, vb->compressed.len);

    txt_file->txt_data_so_far_single += uncompressed_index_at_vb_end - uncompressed_at_vb_start;
    txt_file->disk_so_far            += evb->compressed.len;
    buf_free (&vb->compressed);

    // 3. move the final part (if any) of vb->txt_data, not compressed because its split with the subsequent VB, 
    //    into txt_file->unconsumed_txt, to be compressed in #1 of the next call to this function
    buf_copy (evb, &txt_file->unconsumed_txt, &vb->txt_data, 1, uncompressed_index_at_vb_end, 0, "txt_file->unconsumed_txt");
}
