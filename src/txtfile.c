// ------------------------------------------------------------------
//   txtfile.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.
 
#ifdef __APPLE__
#define off64_t __int64_t // needed for for conda mac - otherwise zlib.h throws compilation errors
#endif
#define Z_LARGE64
#include <errno.h>
#include "genozip.h"
#include "txtfile.h"
#include "vblock.h"
#include "filename.h"
#include "file.h"
#include "strings.h"
#include "progress.h"
#include "codec.h"
#include "bgzf.h"
#include "profiler.h"
#include "segconf.h"
#include "biopsy.h"
#include "zip.h"
#include "dispatcher.h"
#include "zlib/zlib.h"
#include "libdeflate/libdeflate.h"
#include "bzlib/bzlib.h"

#define MAX_TXT_HEADER_LEN ((uint64_t)0xffffffff) // maximum length of txt header - one issue with enlarging it is that we digest it in one go, and the digest module is 32 bit

// PIZ: dump bad vb to disk
rom txtfile_dump_vb (VBlockP vb, rom base_name)
{
    char *dump_filename = MALLOC (strlen (base_name) + 100); // we're going to leak this allocation
    sprintf (dump_filename, "%s.vblock-%u.start-%"PRIu64".len-%u.bad", 
            base_name, vb->vblock_i, vb->vb_position_txt_file, Ltxt);

    if (flag.is_windows) str_replace_letter (dump_filename, strlen(dump_filename), '/', '\\');

    buf_dump_to_file (dump_filename, &vb->txt_data, 1, false, false, false, true);

    return dump_filename;
}

static inline uint32_t txtfile_read_block_plain (VBlockP vb, uint32_t max_bytes)
{
    char *data = BAFTtxt;
    int32_t bytes_read;

    // case: we have data passed to us from file_open_txt_read - handle it first
    if (!Ltxt && evb->scratch.len) {
        memcpy (data, evb->scratch.data, (bytes_read = evb->scratch.len));
        buf_free (evb->scratch);
    }

    // case: normal read
    else {
        bytes_read = read (fileno((FILE *)txt_file->file), data, max_bytes); // -1 if error in libc
        ASSERT (bytes_read >= 0, "read failed from %s: %s", txt_name, strerror(errno));

        if (!bytes_read) { 
            // case external decompressor: inspect its stderr to make sure this is just an EOF and not an error 
            if (file_is_read_via_ext_decompressor (txt_file)) 
                file_assert_ext_decompressor();
            
            txt_file->is_eof = true;
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
        memmove (data, data + 3, bytes_read);
        txt_file->disk_so_far -= 3;
    }
#endif
    Ltxt += bytes_read;

    return (uint32_t)bytes_read;
}

static inline uint32_t txtfile_read_block_gz (VBlockP vb, uint32_t max_bytes)
{
    START_TIMER;

    uint32_t bytes_read = gzfread (BAFTtxt, 1, max_bytes, (gzFile)txt_file->file);
    Ltxt += bytes_read;

    if (bytes_read)
        txt_file->disk_so_far = gzconsumed64 ((gzFile)txt_file->file); // this is actually all the data uncompressed so far, some of it not yet read by us and still waiting in zlib's output buffer
    else
        txt_file->is_eof = true;

    COPY_TIMER (txtfile_read_block_gz);

    return bytes_read;
}

static inline uint32_t txtfile_read_block_bz2 (VBlockP vb, uint32_t max_bytes)
{
    START_TIMER;

    uint32_t bytes_read = BZ2_bzread ((BZFILE *)txt_file->file, BAFTtxt, max_bytes);
    Ltxt += bytes_read;

    if (bytes_read)
        txt_file->disk_so_far = BZ2_consumed ((BZFILE *)txt_file->file); 
    else
        txt_file->is_eof = true;

    COPY_TIMER (txtfile_read_block_bz2);

    return bytes_read;
}

// BGZF: we read *compressed* data into vb->scratch - that will be decompressed now or later, depending on "uncompress". 
// We read data with a *decompressed* size up to max_uncomp. vb->scratch always contains only full BGZF blocks
static inline uint32_t txtfile_read_block_bgzf (VBlockP vb, int32_t max_uncomp /* must be signed */, bool uncompress)
{
    START_TIMER;

    #define uncomp_len prm32[0] // we use vb->compress.param to hold the uncompressed length of the bgzf data in vb->compress

    uint32_t block_comp_len, block_uncomp_len, this_uncomp_len=0;

    if (uncompress)
        vb->gzip_compressor = libdeflate_alloc_decompressor(vb);
        
    int64_t start_uncomp_len = vb->scratch.uncomp_len;

    while (vb->scratch.uncomp_len - start_uncomp_len < max_uncomp - BGZF_MAX_BLOCK_SIZE) {

        buf_alloc (vb, &vb->scratch, BGZF_MAX_BLOCK_SIZE, max_uncomp/4, char, 1.5, "scratch");

        // case: we have data passed to us from file_open_txt_read - handle it first
        if (!Ltxt && evb->scratch.len) {
            block_uncomp_len = evb->scratch.uncomp_len;
            block_comp_len   = evb->scratch.len32;

            // if we're reading a VB (not the txt header) - copy the compressed data from evb to vb
            if (evb != vb) {
                buf_copy (vb, &vb->scratch, &evb->scratch, char,0,0,0);
                buf_free (evb->scratch);
            }

            // add block to list
            buf_alloc (vb, &vb->bgzf_blocks, 1, 1.2 * max_uncomp / BGZF_MAX_BLOCK_SIZE, BgzfBlockZip, 2, "bgzf_blocks");
            BNXT (BgzfBlockZip, vb->bgzf_blocks) = (BgzfBlockZip)
                { .txt_index        = 0,
                  .compressed_index = 0,
                  .txt_size         = block_uncomp_len,
                  .comp_size        = block_comp_len,
                  .is_decompressed  = false };           

            // note: this is the first BGZF block of the file, so vb->bgzf_i remains 0
        }
        else {
            if (!vb->bgzf_blocks.len)
                vb->vb_bgzf_i = txt_file->bgzf_isizes.len; // first bgzf block number for this VB

            block_uncomp_len = (uint32_t)bgzf_read_block (txt_file, BAFT8 (vb->scratch), &block_comp_len, HARD_FAIL);
            
            // check for corrupt data - at this point we've already confirm the file is BGZF so not expecting a different block
            if (block_uncomp_len == BGZF_BLOCK_GZIP_NOT_BGZIP || block_uncomp_len == BGZF_BLOCK_IS_NOT_GZIP) {
                // dump to file
                char dump_fn[strlen(txt_name)+100];
                sprintf (dump_fn, "%s.vb-%u.bad-bgzf.bad-offset-0x%X", txt_name, vb->vblock_i, vb->scratch.len32);
                Buffer dump_buffer = vb->scratch;    // a copy
                dump_buffer.len32 += block_comp_len; // compressed size
                buf_dump_to_file (dump_fn, &dump_buffer, 1, false, false, true, false);

                ABORT ("%s: Invalid BGZF block: block_comp_len=%u. Entire BGZF data of this vblock dumped to %s, bad block stats at offset 0x%X",
                       VB_NAME, block_comp_len, dump_fn, vb->scratch.len32);
            }

            // add block to list - including the EOF block (block_comp_len=BGZF_EOF_LEN block_uncomp_len=0)
            if (block_comp_len && block_uncomp_len/* 0 if truncated */) {
                buf_alloc (vb, &vb->bgzf_blocks, 1, 1.2 * max_uncomp / BGZF_MAX_BLOCK_SIZE, BgzfBlockZip, 2, "bgzf_blocks");
                BNXT (BgzfBlockZip, vb->bgzf_blocks) = (BgzfBlockZip)
                    { .txt_index        = Ltxt,  // after passed-down data and all previous blocks
                      .compressed_index = vb->scratch.len32,
                      .txt_size         = block_uncomp_len,
                      .comp_size        = block_comp_len,
                      .is_decompressed  = !block_uncomp_len }; // EOF block is always considered decompressed           

                vb->scratch.len32 += block_comp_len;   // compressed size
            }

            // case EOF - happens in 2 cases: 1. EOF block (block_comp_len=BGZF_EOF_LEN) or 2. no EOF block (block_comp_len=0)
            if (!block_uncomp_len) {
                txt_file->is_eof = true;
                if (flag.show_bgzf && txt_file->bgzf_flags.has_eof_block) 
                    iprint0 ("IO      vb=0 EOF\n");
                break;
            }
        }

        this_uncomp_len        += block_uncomp_len; // total uncompressed length of data read by this function call
        vb->scratch.uncomp_len += block_uncomp_len; // total uncompressed length of data in vb->compress
        Ltxt                   += block_uncomp_len; // total length of txt_data after adding decompressed vb->scratch (may also include pass-down data)
        txt_file->disk_so_far  += block_comp_len;   

        // we decompress one block a time in the loop so that the decompression is parallel with the disk reading into cache
        if (uncompress) {
            START_TIMER;
            bgzf_uncompress_one_block (vb, BLST (BgzfBlockZip, vb->bgzf_blocks));  
            COPY_TIMER(txtfile_read_block_bgzf_uncompress);
        }
    }

    if (uncompress) {
        buf_free (evb->scratch); 
        libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor);
    }

    COPY_TIMER (txtfile_read_block_bgzf);
    return this_uncomp_len;
#undef param
}

// performs a single I/O read operation - returns number of bytes read
// data is placed in vb->txt_data, except if its BGZF and uncompress=false - compressed data is placed in vb->scratch
static uint32_t txtfile_read_block (VBlockP vb, uint32_t max_bytes,
                                    bool uncompress) // in BGZF, whether to uncompress the data. ignored if not BGZF
{
    START_TIMER;

    if (txt_file->is_eof) return 0; // nothing more to read

    uint32_t bytes_read=0;

    // BGZF note: we read *compressed* data into vb->scratch - that will be decompressed later. we read
    // data with a *decompressed* size up to max_bytes. vb->scratch always contains only full BGZF blocks

    switch (txt_file->codec) {
        case CODEC_NONE : bytes_read = txtfile_read_block_plain (vb, max_bytes); break;
        case CODEC_BGZF : bytes_read = txtfile_read_block_bgzf  (vb, max_bytes, uncompress); break; // bytes_read is in uncompressed terms
        case CODEC_GZ   : bytes_read = txtfile_read_block_gz    (vb, max_bytes); break;
        case CODEC_BZ2  : bytes_read = txtfile_read_block_bz2   (vb, max_bytes); break;

        default: ABORT ("txtfile_read_block: Invalid file type %s (codec=%s)", ft_name (txt_file->type), codec_name (txt_file->codec));
    }

    COPY_TIMER_EVB (read);
    return bytes_read;
}

// iterator on a buffer containing newline-terminated lines
// false means continue iterating, true means stop
char *txtfile_foreach_line (BufferP txt_header,
                            bool reverse, // iterate backwards
                            TxtIteratorCallback callback, 
                            void *cb_param1, void *cb_param2, unsigned cb_param3, // passed as-is to callback
                            int64_t *line_len) // out
{
    if (line_len) *line_len = 0;

    if (!txt_header->len) return NULL;

    char *firstbuf = txt_header->data;
    char *afterbuf = BAFTc (*txt_header);

    char *first = !reverse ? firstbuf : 0;
    char *after = !reverse ? 0 : afterbuf;

    while (1) {
            
        // get one line - searching forward or backwards
        if (!reverse) {
            for (after=first ; after < afterbuf && *after != '\n' ; after++);
            after++; // skip newline
        }
        else {
            for (first=after-2 /* skip final \n */; first >= firstbuf && *first != '\n'; first--);
            first++; // after detected \n or at start of line
        }

        if (!reverse && after > afterbuf) return NULL; // we don't call callback if after>afterbuf - beyond end of line
            
        if (callback (first, after - first, cb_param1, cb_param2, cb_param3)) {
            if (line_len) *line_len = after - first;
            return first;
        }

        if (reverse && first == firstbuf) return NULL; // beginning of line - we called the cb

        if (!reverse) first=after;
        else          after=first;
    }

    return 0; // never reaches here
}   

// default callback from DataTypeProperties.is_header_done: 
// returns header length if header read is complete + sets lines.len, 0 if complete but not existant, -1 not complete yet 
int32_t def_is_header_done (bool is_eof)
{
    ARRAY (char, header, evb->txt_data);
    evb->lines.len = 0; // reset in case we've called this function a number of times (in case of a very large header)
    char prev_char = '\n';

    // check stop condition - a line not beginning with a 'first_char'
    for (int i=0; i < header_len; i++) { // start from 1 back just in case it is a newline, and end 1 char before bc our test is 2 chars
        if (header[i] == '\n') 
            evb->lines.len++;   

        if (prev_char == '\n' && header[i] != DTPT (txt_header_1st_char)) {
            // if we have no header, its an error if we require one
            TxtHeaderRequirement req = DTPT (txt_header_required); 
            ASSINP (i || (req != HDR_MUST && !(req == HDR_MUST_0 && evb->comp_i==0)), 
                    "Error: %s is missing a %s header", txt_name, dt_name (txt_file->data_type));
            return i; // return header length
        }
        prev_char = header[i];
    }

    // case: the entire file is just a header
    if (is_eof && prev_char == '\n') 
        return header_len;

    return -1; // not end of header yet
}

// ZIP main thread: reads txt header into evb->txt_data
void txtfile_read_header (bool is_first_txt)
{
    START_TIMER;

    ASSERT_DT_FUNC (txt_file, is_header_done);

    int32_t header_len;
    uint32_t bytes_read=1 /* non-zero */;

    // read data from the file until either 1. EOF is reached 2. end of txt header is reached
    #define HEADER_BLOCK (256 KB) // we have no idea how big the header will be... read this much at a time
    while ((header_len = (DT_FUNC (txt_file, is_header_done)(bytes_read==0))) < 0) { // we might have data here from txtfile_test_data
        
        if (!bytes_read) {
            if (flags_pipe_in_process_died()) // only works for Linux
                ABORTINP ("Pipe-in process %s (pid=%u) died before the %s header was fully read; only %"PRIu64" bytes were read",
                          flags_pipe_in_process_name(), flags_pipe_in_pid(), dt_name(txt_file->data_type), evb->txt_data.len);
            else
                ABORT ("Unexpected end-of-file while reading the %s header of %s (header_len=%"PRIu64")", 
                       dt_name(txt_file->data_type), txt_name, evb->txt_data.len);
        }

        ASSERT (evb->txt_data.len + HEADER_BLOCK <= MAX_TXT_HEADER_LEN, "%s header is larger than the maximum Genozip supports: %"PRIu64, 
                dt_name (txt_file->data_type), MAX_TXT_HEADER_LEN);

        buf_alloc (evb, &evb->txt_data, HEADER_BLOCK, 0, char, 2, "txt_data");    
        
        if (header_len != HEADER_DATA_TYPE_CHANGED) // note: if HEADER_DATA_TYPE_CHANGED - no need to read more data - we just process the same data again, with a different data type
           bytes_read = txtfile_read_block (evb, HEADER_BLOCK, true);
    }

    // the excess data is for the next vb to read 
    if (evb->txt_data.len > header_len) { 
        buf_copy (evb, &txt_file->unconsumed_txt, &evb->txt_data, char, header_len, 0, "txt_file->unconsumed_txt");
        evb->txt_data.len = header_len; // trim to uncompressed length of txt header

        txt_file->header_size_bgzf = bgzf_copy_unconsumed_blocks (evb); // copy unconsumed or partially consumed bgzf_blocks to txt_file->unconsumed_bgzf_blocks
    }

    txt_file->txt_data_so_far_single = txt_file->header_size = header_len; 

    biopsy_take (evb);

    COPY_TIMER_EVB (txtfile_read_header); // same profiler entry as txtfile_read_header
}

// default "unconsumed" function file formats where we need to read whole \n-ending lines. returns the unconsumed data length
int32_t def_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    ASSERT (*i >= 0 && *i < Ltxt, "*i=%d is out of range [0,%u]", *i, Ltxt);

    int32_t j; for (j=*i; j >= (int32_t)first_i; j--) // use j - automatic var - for speed 
        if (*Btxt (j) == '\n') {
            *i = j;
            return Ltxt -1 - j;
        }

    *i = j;
    return -1; // cannot find \n in the data starting first_i
}

static uint32_t txtfile_get_unconsumed_to_pass_to_next_vb (VBlockP vb)
{
    START_TIMER;

    int32_t pass_to_next_vb_len;
    int32_t last_i = Ltxt-1; // next index to test (going backwards)

    // case: the data is BGZF-compressed in vb->scratch, except for passed down data from prev VB        
    // uncompress one block at a time to see if its sufficient. usually, one block is enough
    if (txt_file->codec == CODEC_BGZF && vb->scratch.len) {

        vb->gzip_compressor = libdeflate_alloc_decompressor(vb);

        for (int block_i=vb->bgzf_blocks.len32 - 1; block_i >= 0; block_i--) {
            BgzfBlockZip *bb = B(BgzfBlockZip, vb->bgzf_blocks, block_i);
            bgzf_uncompress_one_block (vb, bb);

            pass_to_next_vb_len = (DT_FUNC(txt_file, unconsumed)(vb, bb->txt_index, &last_i));
            if (pass_to_next_vb_len >= 0) goto done; // we have the answer (callback returns -1 if it needs more data)
        }
    }

    // test remaining txt_data including passed-down data from previous VB
    pass_to_next_vb_len = (DT_FUNC(vb, unconsumed)(vb, 0, &last_i));

    if (flag.truncate_partial_last_line && pass_to_next_vb_len < 0 && !segconf.running) {
        WARN ("FYI: %s is truncated - its final %s in incomplete. Dropping this final %s.", txt_name, DTPT(line_name), DTPT(line_name));
        txt_file->last_truncated_line_len = Ltxt; 
        Ltxt = pass_to_next_vb_len = 0; // truncate last partial line
    }

    ASSERT (pass_to_next_vb_len >= 0 ||
            segconf.running, // case: we're testing memory and this VB is too small for a single line - return and caller will try again with a larger VB 
            "Reason: failed to find a full line %sin vb=%s data_type=%s txt_data.len=%u txt_file->codec=%s.\n"
            "Known possible causes:\n"
            "- The file is %s %s.\n"
            "- The file is not a %s file.\n"
            "VB dumped: %s\n",  
            DTPT(is_binary) ? "" : "(i.e. newline-terminated) ",
            VB_NAME, dt_name (txt_file->data_type), Ltxt, codec_name (txt_file->codec),
            DTPT(is_binary) ? "truncated but not on the boundary of the" : "missing a newline on the last", DTPT(line_name),
            TXT_DT(REF) ? "FASTA" : dt_name (txt_file->data_type),
            txtfile_dump_vb (vb, txt_name));

done:
    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor);
    
    COPY_TIMER (txtfile_get_unconsumed_to_pass_to_next_vb);
    return (uint32_t)pass_to_next_vb_len;
}

static bool seggable_size_is_modifiable (void)
{
    Codec c = txt_file->source_codec;
    return c==CODEC_GZ || c==CODEC_BGZF || c==CODEC_CRAM || c==CODEC_BZ2 || c==CODEC_BAM;
}
// estimate the size of the txt_data of the file - i.e. the uncompressed data excluding the header - 
// based on the observed or assumed compression ratio of the source compression so far
static void txtfile_set_seggable_size (void)
{
    uint64_t disk_size = txt_file->disk_size ? txt_file->disk_size 
                       : flag.stdin_size     ? flag.stdin_size // user-provided size
                       :                       0; // our estimate will be 0 
    double source_comp_ratio=1;
    switch (txt_file->source_codec) {
        case CODEC_GZ:   // for internal compressors, we use the observed source-compression ratio
        case CODEC_BGZF: 
        case CODEC_BAM:
        case CODEC_BZ2: {
            if (txt_file->is_remote || txt_file->redirected)
                source_comp_ratio = 4;
            else {    
                double plain_len  = txt_file->txt_data_so_far_single + txt_file->unconsumed_txt.len;
                double gz_bz2_len = file_tell (txt_file, HARD_FAIL); // should always work for bz2 or gz. For BZ2 this includes up to 64K read from disk but still in its internal buffers
                
                // case: header is whole BGZF blocks - remove header from calculation to get a better estimate of the seggable compression ratio
                if (txt_file->header_size_bgzf) { 
                    plain_len  -= txt_file->header_size;
                    gz_bz2_len -= txt_file->header_size_bgzf;
                }

                source_comp_ratio = plain_len / gz_bz2_len;
            }
            break;
        }
        
        // external decompressors
        case CODEC_BCF:  source_comp_ratio = 10; break; // note: .bcf files might be compressed or uncompressed - we have no way of knowing as "bcftools view" always serves them to us in plain VCF format. These ratios are assuming the bcf is compressed as it normally is.
        case CODEC_XZ:   source_comp_ratio = 15; break;
        case CODEC_CRAM: source_comp_ratio = 25; break;
        case CODEC_ZIP:  source_comp_ratio = 3;  break;

        case CODEC_NONE: source_comp_ratio = 1;  break;

        default: ABORT ("unspecified txt_file->codec=%s (%u)", codec_name (txt_file->codec), txt_file->codec);
    }
        
    txt_file->est_seggable_size = MAX_(0.0, (double)disk_size * source_comp_ratio - (double)txt_file->header_size);

    if (segconf.running)
        txt_file->txt_data_so_far_single = txt_file->header_size; // roll back as we will re-account for this data in VB=1
}

int64_t txtfile_get_seggable_size (void)
{
    // note: this changes each time a VB is read by the main thread, but since its a single 64B word,
    // a compute thread reading this value will always get a value that makes sense (not using atomic to avoid unnecessary overhead)
    return txt_file->est_seggable_size; 
}

uint64_t txtfile_max_memory_per_vb (void)
{
    return segconf.vb_size - TXTFILE_READ_VB_PADDING;
}

// ZIP main thread
void txtfile_read_vblock (VBlockP vb)
{
    START_TIMER;

    ASSERT_DT_FUNC (txt_file, unconsumed);
    ASSERT ((segconf.vb_size >= ABSOLUTE_MIN_VBLOCK_MEMORY && segconf.vb_size <= ABSOLUTE_MAX_VBLOCK_MEMORY) || segconf.running,
            "Invalid vb_size=%"PRIu64" comp_i(0-based)=%u", segconf.vb_size, z_file->num_txts_so_far-1);

    buf_alloc (vb, &vb->txt_data, 0, segconf.vb_size, char, 1, "txt_data");    

    // read data from the file until either 1. EOF is reached 2. end of block is reached
    uint64_t max_memory_per_vb = txtfile_max_memory_per_vb();
    uint32_t pass_to_next_vb_len = 0;

    // start with using the data passed down from the previous VB (note: copy & free and not move! so we can reuse txt_data next vb)
    if (txt_file->unconsumed_txt.len) {
        uint64_t bytes_moved = MIN_(txt_file->unconsumed_txt.len, max_memory_per_vb);
        buf_copy (vb, &vb->txt_data, &txt_file->unconsumed_txt, char, 0, bytes_moved, "txt_data");
        buf_remove (txt_file->unconsumed_txt, char, 0, bytes_moved);
    }

    bgzf_zip_init_vb (vb); 
    vb->comp_i = flag.zip_comp_i;  // needed for VB_NAME

    bool always_uncompress = flag.pair == PAIR_R2 || // if we're reading the 2nd paired file, fastq_txtfile_have_enough_lines needs the whole data
                             flag.make_reference      || // unconsumed callback for make-reference needs to inspect the whole data
                             segconf.running          ||
                             flag.optimize_DESC       || // fastq_zip_init_vb needs to count lines
                             flag.add_line_numbers    || // vcf_zip_init_vb   needs to count lines
                             flag.biopsy;

    txt_file->header_only = false; // optimistic initialization

    for (bool first=true; ; first=false) {

        bool is_bgzf = (txt_file->codec == CODEC_BGZF);

        uint32_t bytes_requested = MIN_(max_memory_per_vb - Ltxt, 1 GB /* read() can't handle more */);
        bool no_read_expected = is_bgzf && (bytes_requested <= BGZF_MAX_BLOCK_SIZE); // in this case, txtfile_read_block is expected to return 0 

        uint32_t len = (max_memory_per_vb > Ltxt) ? txtfile_read_block (vb, bytes_requested, always_uncompress) : 0;
        
        if (!len && first && !Ltxt) {
            if (vb->vblock_i <= 1) txt_file->header_only = true; // header-only file (the header was already read, and there is no additional data)
                                                                 // note: this doesn't capture header-only 2nd+ bound file
            return;
        }

        // when reading BGZF, we might be filled up even without completely filling max_memory_per_vb 
        // if there is room left for only a partial BGZF block (we can't read partial blocks)
        uint32_t filled_up = max_memory_per_vb - (is_bgzf ? (BGZF_MAX_BLOCK_SIZE - 1): 0);

        if (len && Ltxt < filled_up) continue;  // continue filling up txt_data...

        // case: this is the 2nd file of a fastq pair - make sure it has at least as many fastq "lines" as the first file
        uint32_t my_lines, pair_num_lines, pair_txt_data_len; 
        VBIType pair_vb_i;
        if (flag.pair == PAIR_R2 &&  // we are reading the second file of a fastq file pair (with --pair)
            !fastq_txtfile_have_enough_lines (vb, &pass_to_next_vb_len, &my_lines, &pair_vb_i, &pair_num_lines, &pair_txt_data_len)) { // we don't yet have all the data we need

            // note: the opposite case where R2 has more reads than R1 is caught in fastq_txtfile_have_enough_lines or zip_prepare_one_vb_for_dispatching
            ASSINP ((len || no_read_expected) && Ltxt, "Error: File %s has less FASTQ reads than its R1 mate (vb=%s has %u lines while its pair_vb_i=%d has pair_txt_data_len=%u pair_num_lines=%u; vb=%s Ltxt=%u bytes_requested=%u bytes_read=%u eof=%s max_memory_per_vb=%"PRIu64")", 
                    txt_name, VB_NAME, my_lines, pair_vb_i, pair_txt_data_len/*only set if flag.debug*/, pair_num_lines, VB_NAME, Ltxt, bytes_requested, len, TF(txt_file->is_eof), max_memory_per_vb);

            // if we need more lines - increase memory and keep on reading
            max_memory_per_vb += MAX_((is_bgzf ? 2 * BGZF_MAX_BLOCK_SIZE : 0), max_memory_per_vb / 16); 

            buf_alloc (vb, &vb->txt_data, 0, max_memory_per_vb, char, 1, "txt_data");    
        }
        else
            break;
    }

    if (always_uncompress) buf_free (vb->scratch); // tested by txtfile_get_unconsumed_to_pass_to_next_vb

    // callback to decide what part of txt_data to pass up to the next VB (usually partial lines, but sometimes more)
    // note: even if we haven't read any new data (everything was passed down), we still might data to pass up - eg
    // in FASTA with make-reference if we have a lots of small contigs, each VB will take one contig and pass up the remaining
    if (!pass_to_next_vb_len && Ltxt) {
        pass_to_next_vb_len = txtfile_get_unconsumed_to_pass_to_next_vb (vb);

        // case: return if we're testing memory, and there is not even one line of text  
        if (segconf.running && pass_to_next_vb_len == (uint32_t)-1) {
            buf_copy (evb, &txt_file->unconsumed_txt, &vb->txt_data, char, 0, 0, "txt_file->unconsumed_txt"); 
            buf_free (vb->txt_data);
            return;
        }
    }

    if (pass_to_next_vb_len) {

        // note: we might some unconsumed data, pass it up to the next vb. possibly we still have unconsumed data (can happen if DVCF reject
        // data was passed down from the txt header, greater than max_memory_per_vb)
        buf_insert (evb, txt_file->unconsumed_txt, char, 0, Btxt (Ltxt - pass_to_next_vb_len), pass_to_next_vb_len, "txt_file->unconsumed_txt");
        Ltxt -= pass_to_next_vb_len; 

        // copy unconsumed or partially consumed bgzf_blocks to txt_file->unconsumed_bgzf_blocks
        bgzf_copy_unconsumed_blocks (vb); 

        // if is possible we reached eof but still have pass_up_data - this happens eg in make-reference when a
        // VB takes only one contig from txt_data and pass up the rest - reset eof so that we come back here to process the rest
        txt_file->is_eof = false;
    }

    vb->vb_position_txt_file = txt_file->txt_data_so_far_single;
    vb->is_eof = txt_file->is_eof;
    txt_file->txt_data_so_far_single += Ltxt;

    zip_init_vb (vb);

    if (!txt_file->est_seggable_size || seggable_size_is_modifiable())
        txtfile_set_seggable_size();

    if (!segconf.running) {
        biopsy_take (vb);
        dispatcher_increment_progress ("read", txt_file->est_num_lines ? (Ltxt / MAX_(segconf.line_len,1)) : Ltxt);
    }

    COPY_TIMER (txtfile_read_vblock);
}

DataType txtfile_get_file_dt (rom filename)
{
    FileType ft = flag.stdin_type; // check for --input option

    if (ft == UNKNOWN_FILE_TYPE) // no --input - get file type from filename
        ft = file_get_type (filename);

    return file_get_data_type (ft, true);
}

