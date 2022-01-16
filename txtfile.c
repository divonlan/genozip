// ------------------------------------------------------------------
//   txtfile.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt
 
#ifdef __APPLE__
#define off64_t __int64_t // needed for for conda mac - otherwise zlib.h throws compilation errors
#endif
#define Z_LARGE64
#include <errno.h>
#include "genozip.h"
#include "txtfile.h"
#include "vblock.h"
#include "vcf.h"
#include "file.h"
#include "strings.h"
#include "endianness.h"
#include "crypt.h"
#include "progress.h"
#include "codec.h"
#include "bgzf.h"
#include "mutex.h"
#include "digest.h"
#include "profiler.h"
#include "segconf.h"
#include "biopsy.h"
#include "zlib/zlib.h"
#include "libdeflate/libdeflate.h"
#include "bzlib/bzlib.h"

const char *txtfile_dump_filename (VBlockP vb, const char *base_name, const char *ext) 
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
    buf_dump_to_file (dump_filename, &vb->txt_data, 1, false, false, false, false);

    return dump_filename;
}

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
        ASSERT (bytes_read >= 0, "read failed from %s: %s", txt_name, strerror(errno));

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
        memmove (data, data + 3, bytes_read);
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

// BGZF: we read *compressed* data into vb->compressed - that will be decompressed now or later, depending on uncompress. 
// We read data with a *decompressed* size up to max_uncomp. vb->compressed always contains only full BGZF blocks
static inline uint32_t txtfile_read_block_bgzf (VBlock *vb, int32_t max_uncomp /* must be signed */, bool uncompress)
{
    #define uncomp_len param // we use vb->compress.param to hold the uncompressed length of the bgzf data in vb->compress

    uint32_t block_comp_len, block_uncomp_len, this_uncomp_len=0;

    if (uncompress)
        vb->gzip_compressor = libdeflate_alloc_decompressor(vb);
        
    int64_t start_uncomp_len = vb->compressed.uncomp_len;

    while (vb->compressed.uncomp_len - start_uncomp_len < max_uncomp - BGZF_MAX_BLOCK_SIZE) {

        buf_alloc (vb, &vb->compressed, BGZF_MAX_BLOCK_SIZE, max_uncomp/4, char, 1.5, "compressed");

        // case: we have data passed to us from file_open_txt_read - handle it first
        if (!vb->txt_data.len && evb->compressed.len) {
            block_uncomp_len = evb->compressed.uncomp_len;
            block_comp_len   = evb->compressed.len;

            // if we're reading a VB (not the txt header) - copy the compressed data from evb to vb
            if (evb != vb) {
                buf_copy (vb, &vb->compressed, &evb->compressed, char,0,0,0);
                buf_free (&evb->compressed);
            }

            // add block to list
            buf_alloc (vb, &vb->bgzf_blocks, 1, 1.2 * max_uncomp / BGZF_MAX_BLOCK_SIZE, BgzfBlockZip, 2, "bgzf_blocks");
            NEXTENT (BgzfBlockZip, vb->bgzf_blocks) = (BgzfBlockZip)
                { .txt_index        = 0,
                  .compressed_index = 0,
                  .txt_size         = block_uncomp_len,
                  .comp_size        = block_comp_len,
                  .is_decompressed  = false };           
        }
        else {
            block_uncomp_len = (uint32_t)bgzf_read_block (txt_file, AFTERENT (uint8_t, vb->compressed), &block_comp_len, false);
            
            // check for corrupt data - at this point we've already confirm the file is BGZF so not expecting a different block
            if (block_uncomp_len == BGZF_BLOCK_GZIP_NOT_BGZIP || block_uncomp_len == BGZF_BLOCK_IS_NOT_GZIP) {
                // dump to file
                char dump_fn[strlen(txt_name)+100];
                sprintf (dump_fn, "%s.vb-%u.bad-bgzf.bad-offset-0x%X", txt_name, vb->vblock_i, (uint32_t)vb->compressed.len);
                Buffer dump_buffer = vb->compressed; // a copy
                dump_buffer.len   += block_comp_len; // compressed size
                buf_dump_to_file (dump_fn, &dump_buffer, 1, false, false, true, false);

                ABORT ("Error in txtfile_read_block_bgzf: Invalid BGZF block in vb=%u block_comp_len=%u. Entire BGZF data of this vblock dumped to %s, bad block stats at offset 0x%X",
                       vb->vblock_i, block_comp_len, dump_fn, (uint32_t)vb->compressed.len);
            }

            // add block to list - including the EOF block (block_comp_len=BGZF_EOF_LEN block_uncomp_len=0)
            if (block_comp_len) {
                buf_alloc (vb, &vb->bgzf_blocks, 1, 1.2 * max_uncomp / BGZF_MAX_BLOCK_SIZE, BgzfBlockZip, 2, "bgzf_blocks");
                NEXTENT (BgzfBlockZip, vb->bgzf_blocks) = (BgzfBlockZip)
                    { .txt_index        = vb->txt_data.len, // after passed-down data and all previous blocks
                      .compressed_index = vb->compressed.len,
                      .txt_size         = block_uncomp_len,
                      .comp_size        = block_comp_len,
                      .is_decompressed  = !block_uncomp_len }; // EOF block is always considered decompressed           

                vb->compressed.len += block_comp_len;   // compressed size
            }

            // case EOF - happens in 2 cases: 1. EOF block (block_comp_len=BGZF_EOF_LEN) or 2. no EOF block (block_comp_len=0)
            if (!block_uncomp_len) {
                txt_file->is_eof = true;
                if (flag.show_bgzf && txt_file->bgzf_flags.has_eof_block) 
                    iprint0 ("IO      vb=0 EOF\n");
                break;
            }
        }

        this_uncomp_len           += block_uncomp_len; // total uncompressed length of data read by this function call
        vb->compressed.uncomp_len += block_uncomp_len; // total uncompressed length of data in vb->compress
        vb->txt_data.len          += block_uncomp_len; // total length of txt_data after adding decompressed vb->compressed (may also include pass-down data)
        txt_file->disk_so_far     += block_comp_len;   

        // we decompress one block a time in the loop so that the decompression is parallel with the disk reading into cache
        if (uncompress) bgzf_uncompress_one_block (vb, LASTENT (BgzfBlockZip, vb->bgzf_blocks));  
    }

    if (uncompress) {
        buf_free (&evb->compressed); 
        libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor);
    }

    return this_uncomp_len;
#undef param
}

// performs a single I/O read operation - returns number of bytes read
// data is placed in vb->txt_data, except if its BGZF and uncompress=false - compressed data is placed in vb->compressed
static uint32_t txtfile_read_block (VBlock *vb, uint32_t max_bytes,
                                    bool uncompress) // in BGZF, whether to uncompress the data. ignored if not BGZF
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

// iterator on a buffer containing newline-terminated lines
// false means continue iterating, true means stop
char *txtfile_foreach_line (Buffer *txt_header,
                            bool reverse, // iterate backwards
                            TxtIteratorCallback callback, 
                            void *cb_param1, void *cb_param2, unsigned cb_param3, // passed as-is to callback
                            int64_t *line_len) // out
{
    if (line_len) *line_len = 0;

    if (!txt_header->len) return NULL;

    char *firstbuf = txt_header->data;
    char *afterbuf = AFTERENT (char, *txt_header);

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
    for (int i=0; i < evb->txt_data.len; i++) { // start from 1 back just in case it is a newline, and end 1 char before bc our test is 2 chars
        if (header[i] == '\n') 
            evb->lines.len++;   

        if (prev_char == '\n' && header[i] != DTPT (txt_header_1st_char)) {
            // if we have no header, its an error if we require one
            ASSINP (i || DTPT (txt_header_required) != HDR_MUST, "Error: %s is missing a %s header", txt_name, dt_name (txt_file->data_type));
            return i; // return header length
        }
        prev_char = header[i];
    }

    // case: the entire file is just a header
    if (is_eof && prev_char == '\n') 
        return evb->txt_data.len;

    return -1; // not end of header yet
}

// ZIP main thread: returns the hash of the header
Digest txtfile_read_header (bool is_first_txt)
{
    START_TIMER;

    ASSERT_DT_FUNC (txt_file, is_header_done);

    int32_t header_len;
    uint32_t bytes_read=1 /* non-zero */;

    // read data from the file until either 1. EOF is reached 2. end of txt header is reached
    #define HEADER_BLOCK (256*1024) // we have no idea how big the header will be... read this much at a time
    while ((header_len = (DT_FUNC (txt_file, is_header_done)(bytes_read==0))) < 0) { // we might have data here from txtfile_test_data
        
        if (!bytes_read) {
            if (flags_pipe_in_process_died()) // only works for Linux
                ABORTINP ("Pipe-in process %s (pid=%u) died before the %s header was fully read; only %"PRIu64" bytes were read",
                          flags_pipe_in_process_name(), flags_pipe_in_pid(), dt_name(txt_file->data_type), evb->txt_data.len);
            else
                ABORT ("Error in txtfile_read_header: unexpected end-of-file while reading the %s header of %s (header_len=%"PRIu64")", 
                       dt_name(txt_file->data_type), txt_name, evb->txt_data.len);
        }

        buf_alloc (evb, &evb->txt_data, HEADER_BLOCK, 0, char, 1.15, "txt_data");    
        bytes_read = txtfile_read_block (evb, HEADER_BLOCK, true);
    }

    // the excess data is for the next vb to read 
    if (evb->txt_data.len) 
        buf_copy (evb, &txt_file->unconsumed_txt, &evb->txt_data, char, header_len, 0, "txt_file->unconsumed_txt");
    
    txt_file->txt_data_so_far_single = txt_file->header_size = evb->txt_data.len = header_len; // trim to uncompressed length of txt header

    // md5 header - always digest_ctx_single, digest_ctx_bound only if first component 
    Digest header_digest = DIGEST_NONE;

    if (!flag.data_modified && !flag.gencomp_num) { // note: we don't add generated component's header to the digest
        if (flag.bind && is_first_txt) digest_update (&z_file->digest_ctx_bound, &evb->txt_data, "txt_header:digest_ctx_bound");
        digest_update (&z_file->digest_ctx_single, &evb->txt_data, "txt_header:digest_ctx_single");

        header_digest = digest_snapshot (&z_file->digest_ctx_single);
    }

    biopsy_take (evb);

    COPY_TIMER_VB (evb, txtfile_read_header); // same profiler entry as txtfile_read_header

    return header_digest;
}

// default "unconsumed" function file formats where we need to read whole \n-ending lines. returns the unconsumed data length
int32_t def_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    ASSERT (*i >= 0 && *i < vb->txt_data.len, "*i=%d is out of range [0,%"PRIu64"]", *i, vb->txt_data.len);

    for (; *i >= (int32_t)first_i; (*i)--) 
        if (vb->txt_data.data[*i] == '\n') 
            return vb->txt_data.len-1 - *i;

    return -1; // cannot find \n in the data starting first_i
}

static uint32_t txtfile_get_unconsumed_to_pass_to_next_vb (VBlock *vb)
{
    int32_t pass_to_next_vb_len;
    int32_t last_i = vb->txt_data.len-1; // next index to test (going backwards)

    // case: the data is BGZF-compressed in vb->compressed, except for passed down data from prev VB        
    // uncompress one block at a time to see if its sufficient. usually, one block is enough
    if (txt_file->codec == CODEC_BGZF && vb->compressed.len) {

        vb->gzip_compressor = libdeflate_alloc_decompressor(vb);

        for (int block_i=vb->bgzf_blocks.len-1; block_i >= 0; block_i--) {
            BgzfBlockZip *bb = ENT (BgzfBlockZip, vb->bgzf_blocks, block_i);
            bgzf_uncompress_one_block (vb, bb);

            pass_to_next_vb_len = (DT_FUNC(txt_file, unconsumed)(vb, bb->txt_index, &last_i));
            if (pass_to_next_vb_len >= 0) goto done; // we have the answer (callback returns -1 if no it needs more data)
        }

        libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor);

        // if not found - fall through to test the passed-down data too now
    }

    // test remaining txt_data including passed-down data from previous VB
    pass_to_next_vb_len = (DT_FUNC(txt_file, unconsumed)(vb, 0, &last_i));

    // case: we're testing memory and this VB is too small for a single line - return and caller will try again with a larger VB
    if (segconf.running && pass_to_next_vb_len < 0) return (uint32_t)-1;

    ASSERT (pass_to_next_vb_len >= 0, "Reason: failed to find a full line %sin vb=%u data_type=%s codec=%s.\n"
            "Known possible causes:\n"
            "- The file is %s %s.\n"
            "- The file is not a %s file.\n"
            "VB dumped: %s\n",  
            DTPT(is_binary) ? "" : "(i.e. newline-terminated) ",
            vb->vblock_i, dt_name (txt_file->data_type), codec_name (txt_file->codec),
            DTPT(is_binary) ? "truncated but not on the boundary of the" : "missing a newline on the last", DTPT(line_name),
            TXT_DT(DT_REF) ? "FASTA" : dt_name (txt_file->data_type),
            txtfile_dump_vb (vb, txt_name));

done:
    libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor);
    return (uint32_t)pass_to_next_vb_len;
}

// estimate the size of the txt_data of the file - i.e. the uncompressed data excluding the header - 
// based on the observed or assumed compression ratio of the source compression so far
static void txtfile_set_seggable_size (void)
{
    uint64_t disk_size = txt_file->disk_size ? txt_file->disk_size 
                       : flag.stdin_size     ? flag.stdin_size // user-provided size
                       :                       0; // our estimate will be 0 

    float source_comp_ratio=1;
    switch (txt_file->codec) {
        case CODEC_GZ:   // for internal compressors, we use the observed source-compression ratio
        case CODEC_BGZF:
        case CODEC_BZ2: {
            if (txt_file->is_remote || txt_file->redirected)
                source_comp_ratio = 4;
            else {    
                float plain_len  = txt_file->txt_data_so_far_single + txt_file->unconsumed_txt.len;
                float gz_bz2_len = file_tell (txt_file, false); // should always work for bz2 or gz. For BZ2 this includes up to 64K read from disk but still in its internal buffers
                source_comp_ratio = plain_len / gz_bz2_len;
            }
            break;
        }
        
        case CODEC_BCF:  source_comp_ratio = 10; break; // note: .bcf files might be compressed or uncompressed - we have no way of knowing as "bcftools view" always serves them to us in plain VCF format. These ratios are assuming the bcf is compressed as it normally is.
        case CODEC_XZ:   source_comp_ratio = 15; break;
        case CODEC_CRAM: source_comp_ratio = 25; break;
        case CODEC_ZIP:  source_comp_ratio = 3;  break;
        case CODEC_NONE: source_comp_ratio = 1;  break;

        default: ABORT ("Error in txtfile_set_seggable_size: unspecified txt_file->codec=%s (%u)", codec_name (txt_file->codec), txt_file->codec);
    }
        
    int64_t est_seggable_size = MAX_(0.0, (float)disk_size * source_comp_ratio - (float)txt_file->header_size);
    __atomic_store_n (&txt_file->est_seggable_size, est_seggable_size, __ATOMIC_RELAXED); // atomic loading for thread safety

    if (segconf.running)
        txt_file->txt_data_so_far_single = txt_file->header_size; // roll back as we will re-account for this data in VB=1
}

int64_t txtfile_get_seggable_size (void)
{
    return __atomic_load_n (&txt_file->est_seggable_size, __ATOMIC_RELAXED);
}

uint64_t txtfile_max_memory_per_vb (void)
{
    return segconf.vb_size - TXTFILE_READ_VB_PADDING;
}

// ZIP main threads
void txtfile_read_vblock (VBlock *vb)
{
    START_TIMER;

    ASSERT_DT_FUNC (txt_file, unconsumed);
    ASSERT ((segconf.vb_size >= (MIN_VBLOCK_MEMORY << 20) && segconf.vb_size <= ((uint64_t)MAX_VBLOCK_MEMORY << 20)) || segconf.running,
            "Invalid vb_size=%"PRIu64" component_i(1-based)=%u", segconf.vb_size, z_file->num_txt_components_so_far);

    buf_alloc (vb, &vb->txt_data, 0, segconf.vb_size, char, 1, "txt_data");    

    // read data from the file until either 1. EOF is reached 2. end of block is reached
    uint64_t max_memory_per_vb = txtfile_max_memory_per_vb();
    uint32_t pass_to_next_vb_len=0;

    // start with using the data passed down from the previous VB (note: copy & free and not move! so we can reuse txt_data next vb)
    if (txt_file->unconsumed_txt.len) {
        uint64_t bytes_moved = MIN_(txt_file->unconsumed_txt.len, max_memory_per_vb);
        buf_copy (vb, &vb->txt_data, &txt_file->unconsumed_txt, char, 0, bytes_moved, "txt_data");
        buf_copy (evb, &txt_file->unconsumed_txt, &txt_file->unconsumed_txt, char, bytes_moved, txt_file->unconsumed_txt.len - bytes_moved, NULL);
    }

    bool always_uncompress = flag.pair == PAIR_READ_2 || // if we're reading the 2nd paired file, fastq_txtfile_have_enough_lines needs the whole data
                             flag.make_reference      || // unconsumed callback for make-reference needs to inspect the whole data
                             segconf.running          ||
                             flag.optimize_DESC       || // fastq_zip_read_one_vb needs to count lines
                             flag.add_line_numbers    || // vcf_zip_read_one_vb   needs to count lines
                             flag.biopsy;

    while (1) {

        uint32_t len = max_memory_per_vb > vb->txt_data.len ? txtfile_read_block (vb, MIN_(max_memory_per_vb - vb->txt_data.len, 1<<30 /* read() can't handle more */), always_uncompress) 
                                                            : 0;

        // when reading BGZF, we might be filled up even without completely filling max_memory_per_vb 
        // if there is room left for only a partial BGZF block (we can't read partial blocks)
        uint32_t filled_up = max_memory_per_vb - (txt_file->codec == CODEC_BGZF) * BGZF_MAX_BLOCK_SIZE;

        if (!len || vb->txt_data.len >= filled_up) {  // EOF or we have filled up the allocted memory

            // case: this is the 2nd file of a fastq pair - make sure it has at least as many fastq "lines" as the first file
            uint32_t my_lines, her_lines;
            if (flag.pair == PAIR_READ_2 &&  // we are reading the second file of a fastq file pair (with --pair)
                !fastq_txtfile_have_enough_lines (vb, &pass_to_next_vb_len, &my_lines, &her_lines)) { // we don't yet have all the data we need

                ASSINP (len, "File %s has less FASTQ reads than its R1 counterpart (vb=%u has %u lines while counterpart has %u lines)", 
                        txt_name, vb->vblock_i, my_lines, her_lines);

                ASSERT (vb->txt_data.len, "txt_data.len=0 when reading pair_2 vb=%u", vb->vblock_i);

                // if we need more lines - increase memory and keep on reading
                max_memory_per_vb *= 1.1; 
                buf_alloc (vb, &vb->txt_data, 0, max_memory_per_vb, char, 1, "txt_data");    
            }
            else
                break;
        }
    }

    if (always_uncompress) buf_free (&vb->compressed); // tested by txtfile_get_unconsumed_to_pass_to_next_vb

    // callback to decide what part of txt_data to pass up to the next VB (usually partial lines, but sometimes more)
    // note: even if we haven't read any new data (everything was passed down), we still might data to pass up - eg
    // in FASTA with make-reference if we have a lots of small contigs, each VB will take one contig and pass up the remaining
    if (!pass_to_next_vb_len && vb->txt_data.len) {
        pass_to_next_vb_len = txtfile_get_unconsumed_to_pass_to_next_vb (vb);

        // case: return if we're testing memory, and there is not even one line of text  
        if (segconf.running && pass_to_next_vb_len == (uint32_t)-1) {
            buf_copy (evb, &txt_file->unconsumed_txt, &vb->txt_data, char, 0, 0, "txt_file->unconsumed_txt"); 
            buf_free (&vb->txt_data);
            return;
        }
    }

    if (pass_to_next_vb_len) {

        // note: we might some unconsumed data, pass it up to the next vb. possibly we still have unconsumed data (can happen if DVCF reject
        // data was passed down from the txt header, greater than max_memory_per_vb)
        buf_alloc (evb, &txt_file->unconsumed_txt, pass_to_next_vb_len, 0, char, 1, "txt_file->unconsumed_txt");
        memmove (ENT (char, txt_file->unconsumed_txt, pass_to_next_vb_len), FIRSTENT (char, txt_file->unconsumed_txt), txt_file->unconsumed_txt.len);
        memcpy (FIRSTENT (char, txt_file->unconsumed_txt), ENT (char, vb->txt_data, vb->txt_data.len - pass_to_next_vb_len), pass_to_next_vb_len);
        txt_file->unconsumed_txt.len += pass_to_next_vb_len;

        // now, if our data is bgzf-compressed, txt_data.len becomes shorter than indicated by vb->bgzf_blocks. that's ok - all that data
        // is decompressed and passed-down to the next VB. because it has been decompressed, the compute thread won't try to decompress it again
        vb->txt_data.len -= pass_to_next_vb_len; 

        // if is possible we reached eof but still have pass_up_data - this happens eg in make-reference when a
        // VB takes only one contig from txt_data and pass up the rest - reset eof so that we come back here to process the rest
        txt_file->is_eof = false;
    }

    vb->vb_position_txt_file = txt_file->txt_data_so_far_single;

    vb->recon_size = vb->recon_size_luft = vb->txt_data.len; // initial value. it may change if --optimize / --chain are used, or if dual coordintes - for the other coordinate
    vb->txt_size = vb->txt_data.len; // this copy doesn't change with --optimize / --chain.

    // ZIP of a dual-coordinates file: calculate how much of the VB is rejected lines originating from ##primary_only/##luft_only
    vb->reject_bytes = MIN_(vb->recon_size, txt_file->reject_bytes);
    txt_file->reject_bytes -= vb->reject_bytes;

    txt_file->txt_data_so_far_single += vb->txt_data.len;

    if (DTPT(zip_read_one_vb)) DTPT(zip_read_one_vb)(vb);

    txtfile_set_seggable_size();

    if (!segconf.running)
        biopsy_take (vb);

    COPY_TIMER (txtfile_read_vblock);
}

// read num_lines of the txtfile (after the header), and call test_func for each line. true iff the proportion of lines
// that past the test is at least success_threashold
bool txtfile_test_data (char first_char,            // first character in every header line
                        unsigned num_lines_to_test, // number of lines to test
                        float success_threashold,  // proportion of lines that need to pass the test, for this function to return true
                        TxtFileTestFunc test_func)
{
    uint32_t line_start_i = 0;
    unsigned num_lines_so_far = 0; // number of data (non-header) lines
    unsigned successes = 0;

    #define TEST_BLOCK_SIZE (256 * 1024)

    while (1) {      // read data from the file until either 1. EOF is reached 2. we pass the header + num_lines_to_test lines
        buf_alloc (evb, &evb->txt_data, TEST_BLOCK_SIZE + 1 /* for \0 */, 0, char, 1.2, "txt_data");    

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
    return (float)successes / (float)num_lines_so_far >= success_threashold;
}

DataType txtfile_get_file_dt (const char *filename)
{
    FileType ft = file_get_stdin_type(); // check for --input option

    if (ft == UNKNOWN_FILE_TYPE) // no --input - get file type from filename
        ft = file_get_type (filename);

    return file_get_data_type (ft, true);
}

// get filename of output txt file in genounzip if user didn't specific it with --output
// case 1: outputing a single file - generate txt_filename based on the z_file's name
// case 2: unbinding a genozip into multiple txt files - generate txt_filename of a component file from the
//         component name in SEC_TXT_HEADER 
const char *txtfile_piz_get_filename (const char *orig_name,const char *prefix, bool is_orig_name_genozip)
{
    unsigned fn_len = strlen (orig_name);
    unsigned genozip_ext_len = is_orig_name_genozip ? (sizeof GENOZIP_EXT - 1) : 0;
    char *txt_filename = (char *)MALLOC(fn_len + 10);

    #define EXT2_MATCHES_TRANSLATE(from,to,ext)  \
        ((z_file->data_type==(from) && flag.out_dt==(to) && \
         fn_len >= genozip_ext_len+strlen(ext) && \
         !strcmp (&txt_filename[fn_len-genozip_ext_len- (sizeof ext - 1)], (ext))) ? (int)(sizeof ext - 1) : 0) 

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
                (z_file->z_flags.bgzf && flag.out_dt != DT_BAM) ? ".gz" : ""); // add .gz if --bgzf (except in BAM where it is implicit)

    return txt_filename;
}