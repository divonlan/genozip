// ------------------------------------------------------------------
//   txtfile.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
 
#include "profiler.h"

#ifdef __APPLE__
#define off64_t __int64_t // needed for for conda mac - otherwise zlib.h throws compilation errors
#endif
#define Z_LARGE64
#include <errno.h>
#include <bzlib.h>
#include "genozip.h"
#include "txtfile.h"
#include "vblock.h"
#include "vcf.h"
#include "zfile.h"
#include "file.h"
#include "strings.h"
#include "endianness.h"
#include "crypt.h"
#include "progress.h"
#include "codec.h"
#include "bgzf.h"
#include "mutex.h"
#include "zlib/zlib.h"

static bool is_first_txt = true; 
static uint32_t total_bound_txt_headers_len = 0;

MUTEX (vb_md5_mutex);   // ZIP: used for serializing MD5ing of VBs
static uint32_t vb_md5_last=0; // last vb to be MD5ed 

static const char *txtfile_dump_filename (VBlockP vb, const char *base_name, const char *ext) 
{
    char *dump_filename = MALLOC (strlen (base_name) + 100); // we're going to leak this allocation
    sprintf (dump_filename, "%s.vblock-%u.start-%"PRIu64".len-%u.%s", 
             base_name, vb->vblock_i, vb->vb_position_txt_file, (uint32_t)vb->txt_data.len, ext);
    return dump_filename;
}

// dump bad vb to disk
const char *txtfile_dump_vb (VBlockP vb, const char *base_name)
{
    const char *dump_filename = txtfile_dump_filename (vb, base_name, "bad");
    file_put_buffer (dump_filename, &vb->txt_data, 1);

    return dump_filename;
}

uint32_t txtfile_get_bound_headers_len(void) { return total_bound_txt_headers_len; }

static inline uint32_t txtfile_read_block_plain (VBlock *vb, uint32_t max_bytes)
{
    char *data = AFTERENT (char, vb->txt_data);
    int32_t bytes_read;

    // case: we have data passed to us from file_open_txt_read - handle it first
    if (!vb->txt_data.len && evb->compressed.len) {
        memcpy (data, evb->compressed.data, (bytes_read = evb->compressed.len));
        buf_free (&evb->compressed);
    }

    // case: normal read
    else {
        bytes_read = read (fileno((FILE *)txt_file->file), data, max_bytes); // -1 if error in libc
        ASSERT (bytes_read >= 0, "Error: read failed from %s: %s", txt_name, strerror(errno));

        // bytes_read=0 and we're using an external decompressor - it is either EOF or
        // there is an error. In any event, the decompressor is done and we can suck in its stderr to inspect it
        if (!bytes_read && file_is_read_via_ext_decompressor (txt_file)) {
            file_assert_ext_decompressor();
            txt_file->is_eof = true;
            return 0; // all is good - just a normal end-of-file
        }
    }

    txt_file->disk_so_far += (int64_t)bytes_read;

#ifdef _WIN32
    // in Windows using Powershell, the first 3 characters on an stdin pipe are BOM: 0xEF,0xBB,0xBF https://en.wikipedia.org/wiki/Byte_order_mark
    // these charactes are not in 7-bit ASCII, so highly unlikely to be present natrually in a textual txt file
    if (txt_file->redirected && 
        txt_file->disk_so_far == (int64_t)bytes_read &&  // start of file
        bytes_read >= 3  && 
        (uint8_t)data[0] == 0xEF && 
        (uint8_t)data[1] == 0xBB && 
        (uint8_t)data[2] == 0xBF) {

        // Bomb the BOM
        bytes_read -= 3;
        memcpy (data, data + 3, bytes_read);
        txt_file->disk_so_far -= 3;
    }
#endif
    vb->txt_data.len += bytes_read;

    return (uint32_t)bytes_read;
}

static inline uint32_t txtfile_read_block_gz (VBlock *vb, uint32_t max_bytes)
{
    uint32_t bytes_read = gzfread (AFTERENT (char, vb->txt_data), 1, max_bytes, (gzFile)txt_file->file);
    vb->txt_data.len += bytes_read;

    if (bytes_read)
        txt_file->disk_so_far = gzconsumed64 ((gzFile)txt_file->file); // this is actually all the data uncompressed so far, some of it not yet read by us and still waiting in zlib's output buffer
    else
        txt_file->is_eof = true;

    return bytes_read;
}

static inline uint32_t txtfile_read_block_bz2 (VBlock *vb, uint32_t max_bytes)
{
    uint32_t bytes_read = BZ2_bzread ((BZFILE *)txt_file->file, AFTERENT (char, vb->txt_data), max_bytes);
    vb->txt_data.len += bytes_read;

    if (bytes_read)
        txt_file->disk_so_far = BZ2_consumed ((BZFILE *)txt_file->file); 
    else
        txt_file->is_eof = true;

    return bytes_read;
}

// BGZF: we read *compressed* data into vb->compressed - that will be decompressed later. we read
// data with a *decompressed* size up to max_uncomp. vb->compressed always contains only full BGZF blocks
static inline uint32_t txtfile_read_block_bgzf (VBlock *vb, int32_t max_uncomp /* must be signed */, bool uncompress)
{
    #define uncomp_len param // we use vb->compress.param to hold the uncompressed length of the bgzf data in vb->compress

    if (uncompress) buf_alloc (vb, &vb->txt_data, max_uncomp, 0, "txt_data");

    uint32_t block_comp_len, block_uncomp_len, this_uncomp_len=0;
    while (vb->compressed.uncomp_len < max_uncomp - BGZF_MAX_BLOCK_SIZE) {

        buf_alloc (vb, &vb->compressed, MAX (max_uncomp/4, vb->txt_data.len + BGZF_MAX_BLOCK_SIZE), 1.5, "compressed")

        // case: we have data passed to us from file_open_txt_read - handle it first
        if (!vb->txt_data.len && evb->compressed.len) {
            block_uncomp_len = evb->compressed.uncomp_len;
            block_comp_len   = evb->compressed.len;

            // if we're reading a VB (not the txt header) - copy the compressed data from evb to vb
            if (evb != vb) {
                buf_copy (vb, &vb->compressed, &evb->compressed, 0,0,0,0);
                buf_free (&evb->compressed);
            }

            // add block to list
            buf_alloc_more (vb, &vb->bgzf_blocks, 1, 1.2 * max_uncomp / BGZF_MAX_BLOCK_SIZE, BgzfBlockZip, 2, "bgzf_blocks");
            NEXTENT (BgzfBlockZip, vb->bgzf_blocks) = (BgzfBlockZip)
                { .txt_index        = 0,
                  .compressed_index = 0,
                  .txt_size         = block_uncomp_len,
                  .comp_size        = block_comp_len,
                  .is_decompressed  = false };           
        }
        else {
            block_uncomp_len = (uint32_t)bgzf_read_block (txt_file, AFTERENT (uint8_t, vb->compressed), &block_comp_len, false);

            // case EOF - verify the BGZF end of file marker
            if (!block_uncomp_len) {
                ASSERT (block_comp_len == BGZF_EOF_LEN && !memcmp (AFTERENT (uint8_t, vb->compressed), BGZF_EOF, BGZF_EOF_LEN),
                        "Error: unexpected premature end of bgzf-compressed file %s", txt_name);
                txt_file->is_eof = true;
                break;
            }

            // add block to list
            buf_alloc_more (vb, &vb->bgzf_blocks, 1, 1.2 * max_uncomp / BGZF_MAX_BLOCK_SIZE, BgzfBlockZip, 2, "bgzf_blocks");
            NEXTENT (BgzfBlockZip, vb->bgzf_blocks) = (BgzfBlockZip)
                { .txt_index        = vb->txt_data.len, // after passed-down data and all previous blocks
                  .compressed_index = vb->compressed.len,
                  .txt_size         = block_uncomp_len,
                  .comp_size        = block_comp_len,
                  .is_decompressed  = false };           

            vb->compressed.len   += block_comp_len;   // compressed size
        }

        this_uncomp_len           += block_uncomp_len; // total uncompressed length of data read by this function call
        vb->compressed.uncomp_len += block_uncomp_len; // total uncompressed length of data in vb->compress
        vb->txt_data.len          += block_uncomp_len; // total length of txt_data after adding decompressed vb->compressed (may also include pass-down data)
        txt_file->disk_so_far     += block_comp_len;   

        // we decompress one block a time in the loop so that the decompression is parallel with the disk reading into cache
        if (uncompress) bgzf_uncompress_one_block (vb, LASTENT (BgzfBlockZip, vb->bgzf_blocks));  
    }

    return this_uncomp_len;
#undef param
}

// performs a single I/O read operation - returns number of bytes read
// data is placed in vb->txt_data, except if its BGZF and uncompress=false - compressed data is placed in vb->compressed
static uint32_t txtfile_read_block (VBlock *vb, uint32_t max_bytes,
                                    bool uncompress) // in BGZF, uncompress the data. ignored if not BGZF
{
    START_TIMER;

    if (txt_file->is_eof) return 0; // nothing more to read

    uint32_t bytes_read=0;

    if (file_is_plain_or_ext_decompressor (txt_file)) 
        bytes_read = txtfile_read_block_plain (vb, max_bytes);
    
    // BGZF: we read *compressed* data into vb->compressed - that will be decompressed later. we read
    // data with a *decompressed* size up to max_bytes. vb->compressed always contains only full BGZF blocks
    else if (txt_file->codec == CODEC_BGZF) 
        bytes_read = txtfile_read_block_bgzf (vb, max_bytes, uncompress); // bytes_read is in uncompressed terms

    else if (txt_file->codec == CODEC_GZ) 
        bytes_read = txtfile_read_block_gz (vb, max_bytes);

    else if (txt_file->codec == CODEC_BZ2) 
        bytes_read = txtfile_read_block_bz2 (vb, max_bytes);
    
    else 
        ABORT ("txtfile_read_block: Invalid file type %s (codec=%s)", ft_name (txt_file->type), codec_name (txt_file->codec));
    
    COPY_TIMER_VB (evb, read);
    return bytes_read;
}

// ZIP: called by compute thread to calculate MD5 for one VB - need to serialize VBs using a mutex
void txtfile_md5_one_vb (VBlock *vb)
{
    // wait for our turn
    while (1) {
        mutex_lock (vb_md5_mutex);
        if (vb_md5_last == vb->vblock_i - 1) break; // its our turn now

        // not our turn, wait 10ms and try again
        mutex_unlock (vb_md5_mutex);
        usleep (10000);
    }

    if (command == ZIP) {
        // MD5 of all data up to and including this VB is just the total MD5 of the file so far (as there is no unconsumed data)
        if (flag.bind) md5_update (&z_file->md5_ctx_bound, vb->txt_data.data, vb->txt_data.len);
        md5_update (&z_file->md5_ctx_single, vb->txt_data.data, vb->txt_data.len);

        // take a snapshot of MD5 as per the end of this VB - this will be used to test for errors in piz after each VB  
        vb->md5_hash_so_far = md5_snapshot (flag.bind ? &z_file->md5_ctx_bound : &z_file->md5_ctx_single);
    }
    
    else {
        static bool failed = false; // note: when testing multiple files, if a file fails the test, we don't test subsequent files, so no need to reset this variable

        md5_update (&txt_file->md5_ctx_bound, vb->txt_data.data, vb->txt_data.len);

        // if testing, compare MD5 file up to this VB to that calculated on the original file and transferred through SectionHeaderVbHeader
        // note: we cannot test this unbind mode, because the MD5s are commulative since the beginning of the bound file
        if (!failed && !flag.unbind && !md5_is_zero (vb->md5_hash_so_far)) {
            Md5Hash piz_hash_so_far = md5_snapshot (&txt_file->md5_ctx_bound);

            // warn if VB is bad, but don't exit, so file reconstruction is complete and we can debug it
            if (!md5_is_equal (vb->md5_hash_so_far, piz_hash_so_far)) {

                // dump bad vb to disk
                WARN ("MD5 of reconstructed vblock=%u (%s) differs from original file (%s).\n"
                    "Bad reconstructed vblock has been dumped to: %s\n"
                    "To see the same data in the original file:\n"
                    "   %s %s | head -c %"PRIu64" | tail -c %u > %s",
                    vb->vblock_i, md5_display (piz_hash_so_far), md5_display (vb->md5_hash_so_far), txtfile_dump_vb (vb, z_name),
                    codec_args[txt_file->codec].viewer, file_guess_original_filename (txt_file),
                    vb->vb_position_txt_file + vb->txt_data.len, (uint32_t)vb->txt_data.len, txtfile_dump_filename (vb, z_name, "good"));

                failed = true; // no point in test the rest of the vblocks as they will all fail - MD5 is commulative
            }
        }
    }

    vb_md5_last++; // next please
    mutex_unlock (vb_md5_mutex);
}

// default callback from DataTypeProperties.is_header_done: 
// returns header length if header read is complete + sets lines.len, 0 if complete but not existant, -1 not complete yet 
int32_t def_is_header_done (void)
{
    ARRAY (char, header, evb->txt_data);
    evb->lines.len = 0; // reset in case we've called this function a number of times (in case of a very large header)
    char prev_char = '\n';

    // check stop condition - a line not beginning with a 'first_char'
    for (int i=0; i < evb->txt_data.len; i++) { // start from 1 back just in case it is a newline, and end 1 char before bc our test is 2 chars
        if (header[i] == '\n') 
            evb->lines.len++;   

        if (prev_char == '\n' && header[i] != DTPT (txt_header_1st_char)) {
            // if we have no header, its an error if we require one
            ASSERT (i || DTPT (txt_header_required) != HDR_MUST, "Error: %s is missing a %s header", txt_name, dt_name (txt_file->data_type));
            return i; // return header length
        }
        prev_char = header[i];
    }

    return -1; // not end of header yet
}

// ZIP I/O thread: returns the hash of the header
static Md5Hash txtfile_read_header (bool is_first_txt)
{
    START_TIMER;

    ASSERT_DT_FUNC (txt_file, is_header_done);

    int32_t header_len;
    uint32_t bytes_read=1 /* non-zero */;

    // read data from the file until either 1. EOF is reached 2. end of txt header is reached
    #define HEADER_BLOCK (256*1024) // we have no idea how big the header will be... read this much at time
    while ((header_len = (DT_FUNC (txt_file, is_header_done)())) < 0) { // we might have data here from txtfile_test_data

        ASSERT (bytes_read, "Error in txtfile_read_header: %s: %s file too short - unexpected end-of-file", txt_name, dt_name(txt_file->data_type));

        buf_alloc_more (evb, &evb->txt_data, HEADER_BLOCK, 0, char, 1.15, "txt_data");    
        bytes_read = txtfile_read_block (evb, evb->txt_data.len + HEADER_BLOCK, true);
    }

    buf_free (&evb->compressed); // in case it was bgzf-compressed

    // the excess data is for the next vb to read 
    buf_copy (evb, &txt_file->unconsumed_txt, &evb->txt_data, 1, header_len, 0, "txt_file->unconsumed_txt");

    txt_file->txt_data_so_far_single = evb->txt_data.len = header_len; // trim

    // md5 header - always md5_ctx_single, md5_ctx_bound only if first component 
    if (flag.md5) {
        if (flag.bind && is_first_txt) md5_update (&z_file->md5_ctx_bound, evb->txt_data.data, evb->txt_data.len);
        md5_update (&z_file->md5_ctx_single, evb->txt_data.data, evb->txt_data.len);
    }

    Md5Hash header_md5 = md5_snapshot (&z_file->md5_ctx_single);

    COPY_TIMER_VB (evb, txtfile_read_header); // same profiler entry as txtfile_read_header

    return header_md5;
}

// default "unconsumed" function file formats where we need to read whole \n-ending lines. returns the unconsumed data length
int32_t def_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    ASSERT (*i >= 0 && *i < vb->txt_data.len, "Error in def_unconsumed: *i=%d is out of range [0,%"PRIu64"]", *i, vb->txt_data.len);

    for (; *i >= (int32_t)first_i; (*i)--) 
        if (vb->txt_data.data[*i] == '\n') 
            return vb->txt_data.len-1 - *i;

    return -1; // cannot find \n in the data starting first_i
}

static uint32_t txtfile_get_unconsumed_to_pass_up (VBlock *vb)
{
    int32_t passed_up_len;
    int32_t i=vb->txt_data.len-1; // next index to test (going backwards)

    // case: the data is BGZF-compressed in vb->compressed, except for passed down data from prev VB        
    // uncompress one block at a time to see if its sufficient. usually, one block is enough
    if (vb->compressed.len)
        for (int block_i=vb->bgzf_blocks.len-1; block_i >= 0; block_i--) {
            BgzfBlockZip *bb = ENT (BgzfBlockZip, vb->bgzf_blocks, block_i);
            bgzf_uncompress_one_block (vb, bb);

            passed_up_len = DT_FUNC(txt_file, unconsumed)(vb, bb->txt_index, &i);
            if (passed_up_len >= 0) goto done; // we have the answer (callback returns -1 if no it needs more data)
        }
        // if not found - fall through to test the passed-down data too now

    // test remaining txt_data including passed-down data from previous VB
    passed_up_len = DT_FUNC(txt_file, unconsumed)(vb, 0, &i);

    ASSERT (passed_up_len >= 0, "Error in txtfile_get_unconsumed_to_pass_up: failed to find a single complete line in the entire vb in vb=%u. VB dumped: %s", 
            vb->vblock_i, txtfile_dump_vb (vb, txt_name));

done:
    return (uint32_t)passed_up_len;
}

// ZIP I/O threads
void txtfile_read_vblock (VBlock *vb)
{
    START_TIMER;

    ASSERT_DT_FUNC (txt_file, unconsumed);

    uint64_t pos_before = 0;
    if (file_is_read_via_int_decompressor (txt_file))
        pos_before = file_tell (txt_file);

    buf_alloc (vb, &vb->txt_data, global_max_memory_per_vb, 1, "txt_data");    

    // start with using the data passed down from the previous VB (note: copy & free and not move! so we can reuse txt_data next vb)
    if (buf_is_allocated (&txt_file->unconsumed_txt)) {
        buf_copy (vb, &vb->txt_data, &txt_file->unconsumed_txt, 0 ,0 ,0, "txt_data");
        buf_free (&txt_file->unconsumed_txt);
    }

    // read data from the file until either 1. EOF is reached 2. end of block is reached
    uint64_t max_memory_per_vb = global_max_memory_per_vb;
    uint32_t passed_up_len=0;

    bool always_uncompress = flag.pair == PAIR_READ_2 || // if we're reading the 2nd paired file, fastq_txtfile_have_enough_lines needs the whole data
                             flag.make_reference;        // unconsumed callback for make-reference needs to inspect the whole data

    int32_t block_i=0 ; for (; vb->txt_data.len < max_memory_per_vb; block_i++) {  

        uint32_t len = txtfile_read_block (vb, max_memory_per_vb - vb->txt_data.len, always_uncompress);
        if (!len)  break; // EOF
            
        // case: this is the 2nd file of a fastq pair - make sure it has at least as many fastq "lines" as the first file
        if (flag.pair == PAIR_READ_2 &&  // we are reading the second file of a fastq file pair (with --pair)
            vb->txt_data.len >= max_memory_per_vb && // we are about to exit the loop
            !fastq_txtfile_have_enough_lines (vb, &passed_up_len)) { // we don't yet have all the data we need

            // if we need more lines - increase memory and keep on reading
            max_memory_per_vb *= 1.1; 
            buf_alloc (vb, &vb->txt_data, max_memory_per_vb, 1, "txt_data");    
        }
    }

    if (always_uncompress) buf_free (&vb->compressed); // tested by txtfile_get_unconsumed_to_pass_up

    // callback to decide what part of txt_data to pass up to the next VB (usually partial lines, but sometimes more)
    if (!passed_up_len && vb->txt_data.len) 
        passed_up_len = txtfile_get_unconsumed_to_pass_up (vb);

    // make sure file isn't truncated - if we reached EOF there should be no data to be passed up   
    ASSERT (!passed_up_len || !txt_file->is_eof,
            "Error: input file %s ends abruptly after reading %" PRIu64 " bytes in vb=%u", txt_name, vb->txt_data.len, vb->vblock_i);

    // if we have some unconsumed data, pass it up to the next vb
    if (passed_up_len) {
        buf_copy (evb, &txt_file->unconsumed_txt, &vb->txt_data, 1, // evb, because dst buffer belongs to File
                  vb->txt_data.len - passed_up_len, passed_up_len, "txt_file->unconsumed_txt");

        // now, if our data is bgzf-compressed, txt_data.len becomes shorter than indicated by vb->bgzf_blocks. that's ok - all that data
        // is decompressed and passed-down to the next VB. because it has been decompressed, the compute thread won't try to decompress it again
        vb->txt_data.len -= passed_up_len; 
    }

    vb->vb_position_txt_file = txt_file->txt_data_so_far_single;

    txt_file->txt_data_so_far_single += vb->txt_data.len;
    vb->vb_data_size = vb->txt_data.len; // initial value. it may change if --optimize is used.
    
    if (file_is_read_via_int_decompressor (txt_file))
        vb->vb_data_read_size = file_tell (txt_file) - pos_before; // gz/bz2 compressed bytes read

    if (DTPT(zip_read_one_vb)) DTPT(zip_read_one_vb)(vb);

    // if this is vb=1, we lock the mutex here in the I/O thread before any compute threads start running.
    // this will cause vb>=2 to block on merge, until vb=1 has completed its merge and unlocked it in zip_compress_one_vb()
    if (vb->vblock_i == 1) ctx_vb_1_lock(vb); 

    COPY_TIMER (txtfile_read_vblock);
}

// read num_lines of the txtfile (after the header), and call test_func for each line. true iff the proportion of lines
// that past the test is at least success_threashold
bool txtfile_test_data (char first_char,            // first character in every header line
                        unsigned num_lines_to_test, // number of lines to test
                        double success_threashold,  // proportion of lines that need to pass the test, for this function to return true
                        TxtFileTestFunc test_func)
{
    uint32_t line_start_i = 0;
    unsigned num_lines_so_far = 0; // number of data (non-header) lines
    unsigned successes = 0;

    #define TEST_BLOCK_SIZE (256 * 1024)

    while (1) {      // read data from the file until either 1. EOF is reached 2. we pass the header + num_lines_to_test lines
        buf_alloc_more (evb, &evb->txt_data, TEST_BLOCK_SIZE + 1 /* for \0 */, 0, char, 1.2, "txt_data");    

        uint64_t start_read = evb->txt_data.len;
        txtfile_read_block (evb, TEST_BLOCK_SIZE, true);
        if (start_read == evb->txt_data.len) break; // EOF

        ARRAY (char, str, evb->txt_data); // declare here, in case of a realloc ^ 
        for (uint64_t i=start_read; i < evb->txt_data.len; i++) {

            if (str[i] == '\n') { 
                if (str[line_start_i] != first_char) {  // data line
                    successes += test_func (&str[line_start_i], i - line_start_i);
                    num_lines_so_far++;

                    if (num_lines_so_far == num_lines_to_test) goto done;
                }
                line_start_i = i+1; 
            }
        }
    }
    // note: read data is left in evb->txt_data for the use of txtfile_read_header

done:
    return (double)successes / (double)num_lines_so_far >= success_threashold;
}

// PIZ
void txtfile_write_to_disk (Buffer *buf)
{
    if (!buf->len) return;
    
    if (!flag.test) file_write (txt_file, buf->data, buf->len);

    txt_file->txt_data_so_far_single += buf->len;
    txt_file->disk_so_far            += buf->len;
}

void txtfile_write_one_vblock (VBlockP vb)
{
    START_TIMER;

    if (txt_file->codec == CODEC_BGZF) 
        bgzf_write_to_disk (vb); 
    else
        txtfile_write_to_disk (&vb->txt_data);

    char s1[20], s2[20];

    ASSERTW (vb->txt_data.len == vb->vb_data_size || // files are the same size, expected
             exe_type == EXE_GENOCAT ||              // many genocat flags modify the output file, so don't compare
             !dt_get_translation().is_src_dt,        // we are translating between data types - the source and target txt files have different sizes
             "Warning: vblock_i=%u (num_lines=%u vb_start_line_in_file=%u) had %s bytes in the original %s file but %s bytes in the reconstructed file (diff=%d)", 
             vb->vblock_i, (uint32_t)vb->lines.len, vb->first_line,
             str_uint_commas (vb->vb_data_size, s1), dt_name (txt_file->data_type), str_uint_commas (vb->txt_data.len, s2), 
             (int32_t)vb->txt_data.len - (int32_t)vb->vb_data_size);

    COPY_TIMER (write);
}

// ZIP only - estimate the size of the txt data in this file. affects the hash table size and the progress indicator.
void txtfile_estimate_txt_data_size (VBlock *vb)
{
    uint64_t disk_size = txt_file->disk_size; 

    // case: we don't know the disk file size (because its stdin or a URL where the server doesn't provide the size)
    if (!disk_size) { 
        if (flag.stdin_size) disk_size = flag.stdin_size; // use the user-provided size, if there is one
        else return; // we're unable to estimate if the disk size is not known
    } 
    
    double ratio=1;

    bool is_no_ht_vcf = (txt_file->data_type == DT_VCF && vcf_vb_has_haplotype_data(vb));

    switch (txt_file->codec) {
        // if we decomprssed gz/bz2 data directly - we extrapolate from the observed compression ratio
        case CODEC_GZ:
        case CODEC_BGZF:
        case CODEC_BZ2:  ratio = (double)vb->vb_data_size / (double)vb->vb_data_read_size; break;

        // for compressed files for which we don't have their size (eg streaming from an http server) - we use
        // estimates based on a benchmark compression ratio of files with and without genotype data

        // note: .bcf files might be compressed or uncompressed - we have no way of knowing as 
        // "bcftools view" always serves them to us in plain VCF format. These ratios are assuming
        // the bcf is compressed as it normally is.
        case CODEC_BCF:  ratio = is_no_ht_vcf ? 55 : 8.5; break;

        case CODEC_XZ:   ratio = is_no_ht_vcf ? 171 : 12.7; break;

        case CODEC_BAM:  ratio = 7; break;

        case CODEC_CRAM: ratio = 9; break;

        case CODEC_ZIP:  ratio = 3; break;

        case CODEC_NONE: ratio = 1; break;

        default: ABORT ("Error in txtfile_estimate_txt_data_size: unspecified txt_file->codec=%s (%u)", codec_name (txt_file->codec), txt_file->codec);
    }

    txt_file->txt_data_size_single = disk_size * ratio;
}


// PIZ: called before reading each genozip file
void txtfile_header_initialize(void)
{
    is_first_txt = true;
    vcf_header_initialize(); // we don't yet know the data type, but we initialize the VCF stuff just in case, no harm.
}

// ZIP: reads txt header and writes its compressed form to the GENOZIP file
bool txtfile_header_to_genozip (uint32_t *txt_line_i)
{    
    Md5Hash header_md5 = MD5HASH_NONE;

    // intialize MD5 stuff
    mutex_initialize (vb_md5_mutex);
    if (!z_file->num_txt_components_so_far) vb_md5_last = 0; // reset if we're starting a new z_file

    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    if (DTPT(txt_header_required) == HDR_MUST || DTPT (txt_header_required) == HDR_OK)
        header_md5 = txtfile_read_header (is_first_txt); // reads into evb->txt_data and evb->lines.len
    
    *txt_line_i += (uint32_t)evb->lines.len;

    // for VCF, we need to check if the samples are the same before approving binding (other data types can bind without restriction)
    // for SAM, we check that the contigs specified in the header are consistent with the reference given in --reference/--REFERENCE
    if (!(DT_FUNC_OPTIONAL (txt_file, zip_inspect_txt_header, true)(&evb->txt_data))) { 
        // this is the second+ file in a bind list, but its samples are incompatible
        buf_free (&evb->txt_data);
        return false;
    }

    if (z_file && !flag.test_seg)       
        // we always write the txt_header section, even if we don't actually have a header, because the section
        // header contains the data about the file
        zfile_write_txt_header (&evb->txt_data, header_md5, is_first_txt); // we write all headers in bound mode too, to support --unbind

    // for stats: combined length of txt headers in this bound file, or only one file if not bound
    if (!flag.bind) total_bound_txt_headers_len=0;
    total_bound_txt_headers_len += evb->txt_data.len; 

    z_file->num_txt_components_so_far++; // when compressing

    buf_free (&evb->txt_data);
    
    is_first_txt = false;

    return true; // everything's good
}

// PIZ: reads the txt header from the genozip file and outputs it to the reconstructed txt file
void txtfile_genozip_to_txt_header (const SectionListEntry *sl, uint32_t unbind_component_i, Md5Hash *digest) // NULL if we're just skipped this header (2nd+ header in bound file)
{
    bool show_headers_only = (flag.show_headers && exe_type == EXE_GENOCAT);

    // initialize md5 stuff
    if (unbind_component_i == 0) {
        mutex_initialize (vb_md5_mutex);
        vb_md5_last = 0;
    }

    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    zfile_read_section (z_file, evb, 0, &evb->z_data, "header_section", SEC_TXT_HEADER, sl);

    // handle the GENOZIP header of the txt header section
    SectionHeaderTxtHeader *header = (SectionHeaderTxtHeader *)evb->z_data.data;

    ASSERT (!digest || BGEN32 (header->h.compressed_offset) == crypt_padded_len (sizeof(SectionHeaderTxtHeader)), 
            "Error: invalid txt header's header size: header->h.compressed_offset=%u, expecting=%u", BGEN32 (header->h.compressed_offset), (unsigned)sizeof(SectionHeaderTxtHeader));

    // 1. in unbind mode - we open the output txt file of the component
    // 2. when reading a reference file - we create txt_file here (but don't actually open the physical file)
    if (flag.unbind || flag.reading_reference) {
        ASSERT0 (!txt_file, "Error: not expecting txt_file to be open already in unbind mode or when reading reference");
        
        const char *filename = txtfile_piz_get_filename (header->txt_filename, flag.unbind, false);
        txt_file = file_open (filename, WRITE, TXT_FILE, z_file->data_type);
        FREE (filename); // file_open copies the names
    }

    txt_file->txt_data_size_single = BGEN64 (header->txt_data_size); 
    txt_file->max_lines_per_vb     = BGEN32 (header->max_lines_per_vb);
    txt_file->codec                = header->codec;
    
    if (is_first_txt || flag.unbind) 
        z_file->num_lines = BGEN64 (header->num_lines);

    if (flag.unbind) *digest = header->md5_hash_single; // override md5 from genozip header

    // now get the text of the txt header itself
    if (!show_headers_only)
        zfile_uncompress_section (evb, header, &evb->txt_data, "txt_data", 0, SEC_TXT_HEADER);

    if (z_file->data_type == DT_VCF) {

        if (!show_headers_only) vcf_header_set_globals(z_file->name, &evb->txt_data, false);

        if (flag.drop_genotypes) vcf_header_trim_header_line (&evb->txt_data); // drop FORMAT and sample names

        if (flag.header_one) vcf_header_keep_only_last_line (&evb->txt_data);  // drop lines except last (with field and samples name)
    }
    
    // if we're translating from one data type to another (SAM->BAM, BAM->FASTQ, ME23->VCF etc) translate the txt header 
    DtTranslation trans = dt_get_translation();
    if (trans.txtheader_translator && !show_headers_only) trans.txtheader_translator (&evb->txt_data); 

    // get SEC_BGZF data if needed
    if (flag.bgzf) {
        bgzf_read_and_uncompress_isizes (sl);     
        header = (SectionHeaderTxtHeader *)evb->z_data.data; // re-assign after possible realloc of z_data in bgzf_read_and_uncompress_isizes
    }
    
    // write txt header if not in bound mode, or, in bound mode, we write the txt header, only for the first genozip file
    if ((is_first_txt || flag.unbind) && !flag.no_header && !flag.reading_reference && !flag.genocat_info_only) {

        if (flag.md5) md5_update (&txt_file->md5_ctx_bound, evb->txt_data.data, evb->txt_data.len);

        if (txt_file->codec == CODEC_BGZF) {
            bgzf_calculate_blocks_one_vb (evb, evb->txt_data.len);
            bgzf_compress_vb (evb);
            bgzf_write_to_disk (evb);
        } 
        else
            txtfile_write_to_disk (&evb->txt_data);

        if (!md5_is_zero (header->md5_header)) {
            Md5Hash reconstructed_header_len = md5_do (evb->txt_data.data, evb->txt_data.len);

            if (!md5_is_equal (reconstructed_header_len, header->md5_header)) {
                WARN ("MD5 of reconstructed %s header (%s) differs from original file (%s)\n"
                      "Bad reconstructed vblock has been dumped to: %s\n",
                      dt_name (z_file->data_type), md5_display (reconstructed_header_len), md5_display (header->md5_header),
                      txtfile_dump_vb (evb, z_name));
            }
        }
    }
    
    buf_free (&evb->z_data);
    buf_free (&evb->txt_data);

    z_file->num_txt_components_so_far++;
    is_first_txt = false;
}

DataType txtfile_get_file_dt (const char *filename)
{
    FileType ft = file_get_stdin_type(); // check for --input option

    if (ft == UNKNOWN_FILE_TYPE) // no --input - get file type from filename
        ft = file_get_type (filename, false);

    return file_get_data_type (ft, true);
}

// get filename of output txt file in genounzip if user didn't specific it with --output
// case 1: outputing a single file - generate txt_filename based on the z_file's name
// case 2: unbinding a genozip into multiple txt files - generate txt_filename of a component file from the
//         component name in SEC_TXT_HEADER 
const char *txtfile_piz_get_filename (const char *orig_name,const char *prefix, bool is_orig_name_genozip)
{
    unsigned fn_len = strlen (orig_name);
    unsigned genozip_ext_len = is_orig_name_genozip ? strlen (GENOZIP_EXT) : 0;
    char *txt_filename = (char *)MALLOC(fn_len + 10);

    #define EXT2_MATCHES_TRANSLATE(from,to,ext)  \
        ((z_file->data_type==(from) && flag.out_dt==(to) && strcmp (&txt_filename[fn_len-genozip_ext_len-strlen(ext)], (ext))) ? (int)strlen(ext) : 0) 

    // length of extension to remove if translating, eg remove ".sam" if .sam.genozip->.bam */
    int old_ext_removed_len = EXT2_MATCHES_TRANSLATE (DT_SAM,  DT_BAM,   ".sam") +
                              EXT2_MATCHES_TRANSLATE (DT_SAM,  DT_SAM,   ".bam") +
                              EXT2_MATCHES_TRANSLATE (DT_SAM,  DT_FASTQ, ".sam") +
                              EXT2_MATCHES_TRANSLATE (DT_SAM,  DT_FASTQ, ".bam") +
                              EXT2_MATCHES_TRANSLATE (DT_VCF,  DT_BCF,   ".vcf") +
                              EXT2_MATCHES_TRANSLATE (DT_ME23, DT_VCF,   ".txt");

    sprintf ((char *)txt_filename, "%s%.*s%s%s", prefix,
                fn_len - genozip_ext_len - old_ext_removed_len, orig_name,
                old_ext_removed_len ? file_plain_ext_by_dt (flag.out_dt) : "", // add translated extension if needed
                (flag.bgzf && flag.out_dt != BAM) ? ".gz" : ""); // add .gz if --bgzip (except in BAM where it is implicit)

    return txt_filename;
}