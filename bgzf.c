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

typedef struct __attribute__ ((__packed__)) BgzfHeader {
    uint8_t id1, id2, cm, flg;
    uint32_t mtime;
    uint8_t xfl, os;
    uint16_t xlen;
    uint8_t si1, si2;
    uint16_t slen, bsize;
} BgzfHeader;

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
    if (! *block_size) return BGZF_BLOCK_IS_NOT_GZIP;

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
    zfile_compress_section_data (evb, SEC_BGZF, &txt_file->bgzf_isizes);
    txt_file->bgzf_isizes.len /= sizeof (uint16_t); // restore

    // actually write to disk
    zip_output_processed_vb (evb, NULL, false, PD_VBLOCK_DATA);  
}

void bgzf_uncompress_one_block (VBlock *vb, BgzfBlock *bb)
{
    if (bb->is_decompressed) return; // already decompressed - nothing to do

    BgzfHeader *h = (BgzfHeader *)ENT (char, vb->compressed, bb->compressed_index);

    ASSERT (h->id1==31 && h->id2==139, "Error in bgzf_uncompress_vb: corrupt bgzf data in vb->compressed: vb=%u compressed_index=%u", vb->vblock_i, bb->compressed_index);

    if (!vb->bgzf_decompressor) vb->bgzf_decompressor = libdeflate_alloc_decompressor();

    if (flag.show_bgzf)
        fprintf (stderr, "%-7s vb=%u i=%u compressed_index=%u size=%u txt_data_index=%u size=%u ",
                 arch_am_i_io_thread() ? "IO" : "COMPUTE", vb->vblock_i, 
                 ENTNUM (vb->bgzf_blocks, bb), bb->compressed_index, bb->comp_size, bb->txt_data_index, bb->uncomp_size);

//fprintf (stderr, "\nbb->txt_data_index=%u bb->uncomp_size=%u  xxx %s\n", bb->txt_data_index, bb->uncomp_size, buf_desc (&vb->compressed));
    enum libdeflate_result ret = 
        libdeflate_deflate_decompress (vb->bgzf_decompressor, 
                                       h+1, bb->comp_size - sizeof(*h),   // compressed
                                       ENT (char, vb->txt_data, bb->txt_data_index), bb->uncomp_size, NULL); // uncompressed

    ASSERT (ret == LIBDEFLATE_SUCCESS, "Error in bgzf_uncompress_vb: libdeflate_deflate_decompress failed: %s", libdeflate_error(ret));

    bb->is_decompressed = true;

    if (flag.show_bgzf)
        #define C(i) ((bb->txt_data_index + i < vb->txt_data.len) ? char_to_printable (*ENT (char, vb->txt_data, bb->txt_data_index + i)) : 0) 
        fprintf (stderr, "txt_data[5]=%c%c%c%c%c\n", C(0), C(1), C(2), C(3), C(4));
}

// ZIP: called from the compute thread: zip_compress_one_vb and I/O thread: txtfile_read_block_bgzf
void bgzf_uncompress_vb (VBlock *vb)
{
    START_TIMER;

    for (uint32_t block_i=0; block_i < vb->bgzf_blocks.len ; block_i++) 
        bgzf_uncompress_one_block (vb, ENT (BgzfBlock, vb->bgzf_blocks, block_i));

    buf_free (&vb->compressed); // now that we are finished decompressing we can free it

    COPY_TIMER (bgzf_uncompress_vb);
}

//---------
// PIZ SIDE
//---------

