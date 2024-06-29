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
#include "txtfile.h"
#include "file.h"
#include "codec.h"
#include "bgzf.h"
#include "biopsy.h"
#include "zip.h"
#include "arch.h"
#include "dispatcher.h"
#include "zlib/zlib.h"
#include "libdeflate_1.19/libdeflate.h"
#include "bzlib/bzlib.h"
#include "igzip/igzip_lib.h"

#define MAX_TXT_HEADER_LEN ((uint64_t)0xffffffff) // maximum length of txt header - one issue with enlarging it is that we digest it in one go, and the digest module is 32 bit

#define TXTFILE_READ_VB_PADDING 16 // we need this quantity of unused bytes at the end of vb.txt_data

// PIZ: dump bad vb to disk
StrTextLong txtfile_dump_vb (VBlockP vb, rom base_name)
{
    StrTextLong dump_filename;

    snprintf (dump_filename.s, sizeof (dump_filename.s), "%.*s.vblock-%u.start-%"PRIu64".len-%u.bad", 
              (int)sizeof (dump_filename.s)-100, base_name, vb->vblock_i, vb->vb_position_txt_file, Ltxt);

    if (flag.is_windows) str_replace_letter (dump_filename.s, strlen(dump_filename.s), '/', '\\');

    buf_dump_to_file (dump_filename.s, &vb->txt_data, 1, false, false, false, true);

    return dump_filename;
}

// returns the requested number of bytes, except if eof in which case it could be less.
uint32_t txtfile_fread (FileP file, 
                        FILE *fp, // note: non-NULL if different from file->file (when re-reading) 
                        void *addr, uint32_t size, int64_t *disk_so_far)
{
    ASSERTNOTNULL (addr);
    if (!fp) fp = (FILE *)file->file;

    uint32_t bytes = fread (addr, 1, size, fp);
    if (disk_so_far) *disk_so_far += bytes;   

    ASSERT (bytes == size || !ferror (fp), "Error while reading %s codec=%s on filesystem=%s - requested %u bytes but read only %u: (%u)%s", 
            file->basename, codec_name (file->codec), arch_get_filesystem_type (file).s, size, bytes, errno, strerror (errno));

    // note: since we now took care of errors, we know that it feof iff bytes < size    
    return bytes;
}

void txtfile_fwrite (const void *data, uint32_t size)
{
    if (!size) return; // nothing to do

    ASSERTNOTNULL (txt_file);
    ASSERTNOTNULL (txt_file->file);
    ASSERTNOTNULL (data);
    
    uint32_t bytes = fwrite (data, 1, size, (FILE *)txt_file->file); 

    // if we're streaming our txt output to another process and that process has ended prematurely or 
    // otherwise closed the pipe, then exit quietly (note: sometimes the shell will kill us before we reach here) 
    if (bytes < size && errno == EPIPE) exit (EXIT_DOWNSTREAM_LOST);

    // error if failed to write to file
    ASSERT (bytes == size, "Error writing to %s on filesystem=%s - requested %u bytes but wrote only %u: (%u)%s", 
            txt_file->basename, arch_get_filesystem_type (txt_file).s, size, bytes, errno, strerror (errno));
}

static inline uint32_t txtfile_read_block_plain (VBlockP vb, uint32_t max_bytes)
{
    char *data = BAFTtxt;
    int32_t bytes_read;

    // case: we have data passed to us from txtfile_discover_gz_codec - handle it first (possibly txt_data already contains data passed down from previous VBs)
    if (txt_file->gz_data.len) {
        bytes_read = MIN_(txt_file->gz_data.len32, max_bytes);
        memcpy (BAFTtxt, B1STc (txt_file->gz_data), bytes_read);
        buf_remove (txt_file->gz_data, char, 0, bytes_read);
    }

    // case: normal read
    else {
        bytes_read = txtfile_fread (txt_file, NULL, data, max_bytes, &txt_file->disk_so_far);

        if (!bytes_read) { 
            // case external decompressor: inspect its stderr to make sure this is just an EOF and not an error 
            if (is_read_via_ext_decompressor (txt_file)) 
                file_assert_ext_decompressor();
            
            txt_file->no_more_blocks = true;
        }
    }

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
    }
#endif
    Ltxt += bytes_read;

    return (uint32_t)bytes_read;
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

// chuck size chosen to be equal to the default disk read-ahead buffer in Linux (eg /sys/block/sda/queue/read_ahead_kb) 
// so that decompression is parallelized with disk read-ahead buffer filling (a bigger buffer would cause the disk to be idle while we are still decompressing)
#define IGZIP_CHUNK (128 KB) 

static void txtfile_initialize_igzip (FileP file)
{
    ASSERTNOTINUSE (file->igzip_state);

    buf_alloc_exact_zero (evb, file->igzip_state, 1, struct inflate_state, "txt_file->igzip_state");
    struct inflate_state *state = B1ST (struct inflate_state, file->igzip_state);

    isal_inflate_init (state);
    state->crc_flag = ISAL_GZIP;
}

void txtfile_discover_gz_codec (FileP file)
{
    buf_alloc (evb, &file->gz_data, 0, MAX_(IGZIP_CHUNK, GZIL_MAX_BLOCK_SIZE), char, 0, "gz_data");
    
    // read the first potential BGZF block to test if this is GZ or BGZF 
    // note: we read if --no-bgzf to capture the data for z_file->gz_header
    GzStatus status = bgzf_read_block (file, true);

    // case: this is a BGZF block
    // note: we keep the still-compressed data in vb->scratch for later consumption
    if (!flag.no_bgzf && status == GZ_SUCCESS && file->gz_data.uncomp_len > 0) {
        if (file->source_codec != CODEC_CRAM && file->source_codec != CODEC_BAM && file->source_codec != CODEC_BCF) 
            file->source_codec = CODEC_BGZF;

        file->codec = CODEC_BGZF;
        bgzf_initialize_discovery (file);
    }

    // for regulars files, we already skipped 0 size files. This can happen in STDIN
    else if (status == GZ_SUCCESS && file->gz_data.uncomp_len == 0) {
        ASSINP (!flags_pipe_in_process_died(), // only works for Linux
                "Pipe-in process %s (pid=%u) died without sending any data",
                flags_pipe_in_process_name(), flags_pipe_in_pid());

        ABORTINP ("No data exists in input file %s", file->name ? file->name : FILENAME_STDIN);
    }
    
    // case: this is non-BGZF GZIP format 
    else if (flag.no_bgzf || status == GZ_IS_GZIP_NOT_BGZF) {      
        // case: this is FASTQ (judged by the filename) that is GZIL
        bool is_eof = false;
        if (!flag.no_bgzf && file->data_type == DT_FASTQ && 
            gzil_read_block (file, true, &is_eof) != GZ_IS_NOT_GZIL) 
            
            file->codec = file->source_codec = CODEC_GZIL;

        // case: neither BGZF or GZIL - treat as normal GZ
        else 
            file->codec = file->source_codec = CODEC_GZ;        
    } 

    // case: this is not GZIP format at all. treat as a plain file, and put the data read in vb->scratch 
    // for later consumption is txtfile_read_block_plain
    else if (status == GZ_IS_NOT_GZIP) {

        #define BZ2_MAGIC "BZh"
        #define XZ_MAGIC  (char[]){ 0xFD, '7', 'z', 'X', 'Z', 0 }
        #define ZIP_MAGIC (char[]){ 0x50, 0x4b, 0x03, 0x04 }
        #define ORA_MAGIC (char[]){ 0x49, 0x7c } // https://support-docs.illumina.com/SW/ORA_Format_Specification/Content/SW/ORA/ORAFormatSpecification.htm

        // we already open the file, so not easy to re-open with BZ2_bzopen as it would require injecting the read data into the BZ2 buffers
        if (str_isprefix_(STRb(file->gz_data), BZ2_MAGIC, 3)) 
            ABORTINP0 ("The data seems to be in bz2 format. Please use --input to specify the type (eg: \"genozip --input sam.bz2\")");

        else if (str_isprefix_(STRb(file->gz_data), XZ_MAGIC, 6)) {
            if (file->redirected) ABORTINP0 ("Compressing piped-in data in xz format is not currently supported");
            if (file->is_remote) ABORTINP0 ("The data seems to be in xz format. Please use --input to specify the type (eg: \"genozip --input sam.xz\")");
            ABORTINP0 ("The data seems to be in xz format. Please use --input to specify the type (eg: \"genozip --input sam.xz\")");
        }

        else if (str_isprefix_(STRb(file->gz_data), ZIP_MAGIC, 4)) {
            if (file->redirected) ABORTINP0 ("Compressing piped-in data in zip format is not currently supported");
            if (file->is_remote) ABORTINP0 ("The data seems to be in zip format. Please use --input to specify the type (eg: \"genozip --input generic.zip\")");
            ABORTINP0 ("The data seems to be in zip format. Please use --input to specify the type (eg: \"genozip --input generic.zip\")");
        }

        else if (str_isprefix_(STRb(file->gz_data), ORA_MAGIC, 2)) {
            if (file->redirected) ABORTINP0 ("Compressing piped-in data in ora format is not currently supported");
            if (file->is_remote) ABORTINP0 ("The data seems to be in ora format. Please use --input to specify the type (eg: \"genozip --input fastq.ora\")");
            ABORTINP0 ("The data seems to be in ora format. Please use --input to specify the type (eg: \"genozip --input fastq.ora\")");
        }

        file->codec = CODEC_NONE;
    }

    else
        ABORT ("Invalid status=%u", status);

    // if this is R2 we are going to uncompress in the main thread. IGZIP is a faster method for doing so 
    // than BGZF or GZIP, bc it is better at parallelizing disk read-aheads and decompression. The only reason
    // to keep BGZF is if we want to store BGZF isizes for exact reconstruction, which is only possible if we discovered
    // the library. BGZF library discovery has not yet occurred for R2, so we take the R1 results as proxy 
    // (if the proxying is wrong - either we will compress unneccsary slowly with BGZF instead of IGZIP, or we will 
    // incorrectly compress with IGZIP and drop the BGZF isizes preventing exact reconstruction - that's ok)
    bool is_pair2 = flag.pair && ((flag.zip_comp_i == FQ_COMP_R2 && Z_DT(FASTQ)) ||  // note: flag.pair is not incremented yet; z_file only exists if this 2nd+ component so test that first
                                  (flag.zip_comp_i == SAM_COMP_FQ01 && (Z_DT(BAM) || Z_DT(SAM)))); 
    if ((is_pair2 && (file->codec == CODEC_GZIL || (z_file->comp_codec[flag.zip_comp_i-1] == CODEC_BGZF && z_file->comp_bgzf[flag.zip_comp_i-1].level == BGZF_COMP_LEVEL_UNKNOWN))) ||    
        // likewise for --make-reference: we uncompress by main thread, and we don't care about retaining BGZF isizes
        (flag.make_reference && (file->codec == CODEC_GZIL || file->codec == CODEC_BGZF))) {

        file->gunzip_method = CODEC_GZ;
    }

    else 
        file->gunzip_method = file->codec;

    if (file->gunzip_method == CODEC_GZ)
        txtfile_initialize_igzip (file);
}

// ZIP main thread: called after txt and z are open, and txt codecs have been discovered.
void txtfile_zip_finalize_codecs (void)
{
    ASSERTNOTNULL (z_file);
    ASSERTNOTNULL (txt_file);
    
    if (flag.zip_comp_i < MAX_NUM_COMPS) { // for stats
        z_file->comp_codec[flag.zip_comp_i]         = txt_file->codec;
        z_file->comp_source_codec[flag.zip_comp_i]  = txt_file->source_codec;
        z_file->comp_gunzip_method[flag.zip_comp_i] = txt_file->gunzip_method;
    
        // copy GZ header (but not if BGZF, GZIL): data should be in gz_data txtfile_discover_gz_codec
        if (TXT_IS_GZ && txt_file->gz_data.len >= 12)
            memcpy (z_file->gz_header, B1ST8 (txt_file->gz_data), 12);
    }

    // note: for BGZF, we report in bgzf_finalize_discovery as we don't yet know the library/level here
    if ((flag.show_gz || flag.show_bgzf) && txt_file->gunzip_method != CODEC_BGZF) {
        iprintf ("%s: txt_codec=%s", txt_file->basename, txtfile_codec_name (z_file, flag.zip_comp_i).s); 
        if (flag.show_gz) { iprint0 ("\n"); exit_ok; };

        iprintf (" gunzip_method=%s\n", codec_name (txt_file->gunzip_method));
    }
}

// runs in main thread, reads and uncompressed GZ, and populates txt_data for vb
static uint32_t txtfile_read_block_igzip (VBlockP vb, uint32_t max_bytes)
{
    START_TIMER;
    ASSERTISALLOCED (txt_file->gz_data);

    struct inflate_state *state = B1ST (struct inflate_state, txt_file->igzip_state);

    // top up gz_data
    int32_t bytes_read = (txt_file->gz_data.len32 < IGZIP_CHUNK) 
        ? txtfile_fread (txt_file, NULL, BAFTc(txt_file->gz_data), IGZIP_CHUNK - txt_file->gz_data.len32, &txt_file->disk_so_far)
        : 0;
    
    ASSERT (!ferror((FILE *)txt_file->file) && bytes_read >= 0, "Error reading GZ file %s on filesystem=%s: %s", 
            txt_name, arch_get_txt_filesystem().s, strerror (errno));

    txt_file->gz_data.len32 += bytes_read; // yet-uncompressed data read from disk

    { START_TIMER

    state->next_in   = B1ST8 (txt_file->gz_data);
    state->avail_in  = txt_file->gz_data.len32;
    state->next_out  = BAFT8 (vb->txt_data);
    state->avail_out = max_bytes;

    // case: happens in blocked-GZ: in the previous call to this function we read the entire GZ-block, 
    // but unluckily just short of reading the checksum data. Now that we read more data from disk, 
    // we can verify the checksum, and if successful, the state will change to ISAL_BLOCK_FINISH
    if (state->block_state == ISAL_CHECKSUM_CHECK) {
        int ret = isal_inflate (state); // new gzip header in a file that has concatented gzip compressions
        ASSERT (ret == ISAL_DECOMP_OK, "isal_inflate failed checksum: %s avail_in=%u avail_out=%u",
                isal_error (ret), txt_file->gz_data.len32, max_bytes);
    }

    // case: happens in blocked-GZ: we decompressed an entire GZ-block and verified the 
    // checksum (either in isal_inflate below in the previous call to this function, 
    // or in isal_inflate above in this call). now we need to move on to the next GZ block.
    if (state->block_state == ISAL_BLOCK_FINISH) 
        isal_inflate_reset (state);

    int ret = isal_inflate (state); // new gzip header in a file that has concatented gzip compressions
    ASSERT (ret == ISAL_DECOMP_OK || ret == ISAL_END_INPUT, "isal_inflate error: %s avail_in=%u avail_out=%u",
            isal_error (ret), txt_file->gz_data.len32, max_bytes);

    COPY_TIMER(igzip_uncompress_during_read); }

    uint32_t gz_data_consumed = BNUM (txt_file->gz_data, state->next_in);

    // for stats: read and save the isize from the gzip footer (of the first 2 gzip blocks) 
    if (state->block_state == ISAL_BLOCK_FINISH && z_file && gz_data_consumed >= 4) 
        for (int i=0; i <= 1; i++)
            if (!z_file->gz_isize[flag.zip_comp_i][i]) {
                z_file->gz_isize[flag.zip_comp_i][i] = LTEN32 (GET_UINT32 (B8(txt_file->gz_data, gz_data_consumed - 4)));
                break;
            }

    buf_remove (txt_file->gz_data, char, 0, gz_data_consumed);

    inc_disk_gz_uncomp_or_trunc (txt_file, gz_data_consumed);

    txt_file->no_more_blocks = (!state->avail_in && feof ((FILE *)txt_file->file));

    Ltxt = BNUMtxt (state->next_out);

    COPY_TIMER (txtfile_read_block_igzip);

    return max_bytes - state->avail_out; // uncompressed data length
}

static inline uint32_t txtfile_read_block_bz2 (VBlockP vb, uint32_t max_bytes)
{
    START_TIMER;

    uint32_t bytes_read = BZ2_bzread ((BZFILE *)txt_file->file, BAFTtxt, max_bytes);
    Ltxt += bytes_read;

    if (bytes_read) 
        txt_file->disk_so_far = BZ2_consumed ((BZFILE *)txt_file->file); 
    else
        txt_file->no_more_blocks = true;

    COPY_TIMER (txtfile_read_block_bz2);

    return bytes_read;
}

// BGZF: we read *compressed* data into vb->scratch - that will be decompressed now or later, depending on "uncompress". 
// We read data with a *decompressed* size up to max_uncomp. vb->scratch always contains only full BGZF blocks
static inline uint32_t txtfile_read_block_bgz (VBlockP vb, int32_t max_uncomp /* must be signed */, bool uncompress)
{
    START_TIMER;

    uint32_t this_uncomp_len=0;

    if (uncompress)
        vb->gzip_compressor = libdeflate_alloc_decompressor(vb, __FUNCLINE);
        
    int64_t start_uncomp_len = vb->scratch.uncomp_len;
    int32_t max_block_size = TXT_IS_BGZF ? BGZF_MAX_BLOCK_SIZE : GZIL_MAX_BLOCK_SIZE;

    // scratch contains gz-compressed data; we use .uncomp_len to track its uncompress length
    buf_alloc (vb, &vb->scratch, 0, max_uncomp/2, char, 0, "scratch");

    while (vb->scratch.uncomp_len - start_uncomp_len <= max_uncomp - max_block_size && 
           !txt_file->no_more_blocks) {
        
        bool is_eof = false; // only used for GZIL
        GzStatus status = TXT_IS_BGZF ? bgzf_read_block (txt_file, false)
                                      : gzil_read_block (txt_file, false, &is_eof);

        uint32_t this_block_start = vb->scratch.len32;
        buf_add_more (vb, &vb->scratch, txt_file->gz_data.data, txt_file->gz_data.comp_len, "scratch");

        // check for corrupt data - at this point we've already confirm the file is BGZF so not expecting a different block
        if (status != GZ_SUCCESS) {
            // dump to file
            char dump_fn[strlen(txt_name)+100];
            snprintf (dump_fn, sizeof (dump_fn), "%s.vb-%u.bad-%s.bad-offset-0x%x", 
                      txt_name, vb->vblock_i, codec_name (txt_file->codec), this_block_start);
            
            buf_dump_to_file (dump_fn, &vb->scratch, 1, false, false, true, false);

            ABORT ("%s: Invalid %s block: block_comp_len=%u. Entire data of this vblock dumped to %s, bad block stats at offset 0x%x",
                    VB_NAME, codec_name (txt_file->codec), txt_file->gz_data.comp_len, dump_fn, this_block_start);
        }

        // add block to list - including the EOF block (block_comp_len=BGZF_EOF_LEN block_uncomp_len=0)
        if (txt_file->gz_data.comp_len/* note: if is 0 if truncated or EOF with no EOF block */) {
            buf_alloc (vb, &vb->gz_blocks, 1, MAX_(1000, 1.2 * max_uncomp / max_block_size), GzBlockZip, 2, "gz_blocks");
            BNXT (GzBlockZip, vb->gz_blocks) = (GzBlockZip)
                { .txt_index        = Ltxt,  // after passed-down data and all previous blocks
                  .compressed_index = this_block_start,
                  .txt_size         = txt_file->gz_data.uncomp_len,
                  .comp_size        = txt_file->gz_data.comp_len,
                  .is_decompressed  = !txt_file->gz_data.uncomp_len, // EOF block is always considered decompressed
                  .is_eof           = is_eof }; 

            // case EOF block: are not going to decompress the block, so account for it here
            if (!txt_file->gz_data.uncomp_len) 
                inc_disk_gz_uncomp_or_trunc (txt_file, txt_file->gz_data.comp_len);
        }
        
        // case EOF - happens in 2 cases: 1. EOF block (block_comp_len=BGZF_EOF_LEN) or 2. no EOF block (block_comp_len=0)
        if (!txt_file->gz_data.uncomp_len) {
            txt_file->no_more_blocks = true;

            if (flag.show_bgzf && txt_file->bgzf_flags.has_eof_block) 
                iprint0 ("IO      vb=0 EOF\n");
        }

        else {
            this_uncomp_len        += txt_file->gz_data.uncomp_len; // total uncompressed length of data read by this function call
            vb->scratch.uncomp_len += txt_file->gz_data.uncomp_len; // total uncompressed length of data in vb->compress
            Ltxt                   += txt_file->gz_data.uncomp_len; // total length of txt_data after adding decompressed vb->scratch (may also include pass-down data)

            // we decompress one block a time in the loop so that the decompression is parallel with the disk reading into cache
            if (uncompress) {
                START_TIMER;
                bgz_uncompress_one_block (vb, BLST (GzBlockZip, vb->gz_blocks), txt_file->codec);  
                COPY_TIMER(bgz_uncompress_during_read); 
            }
        }

        buf_remove (txt_file->gz_data, uint8_t, 0, txt_file->gz_data.comp_len);
        txt_file->gz_data.comp_len = txt_file->gz_data.uncomp_len = 0;
    }

    if (uncompress) {
        buf_free (vb->scratch); 
        libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor, __FUNCLINE);
    }

    COPY_TIMER (txtfile_read_block_bgz);
    return this_uncomp_len;
}

// performs a single I/O read operation - returns number of bytes read
// data is placed in vb->txt_data, except if its BGZF and uncompress=false - compressed data is placed in vb->scratch
static uint32_t txtfile_read_block (VBlockP vb, uint32_t max_bytes,
                                    bool uncompress) // in BGZF/GZIL, whether to uncompress the data. ignored if not BGZF/GZIL
{
    START_TIMER;

    if (txt_file->no_more_blocks) return 0; // nothing more to read

    uint32_t uncomp_len=0;

    // BGZF note: we read *compressed* data into vb->scratch - that will be decompressed later. we read
    // data with a *decompressed* size up to max_bytes. vb->scratch always contains only full BGZF blocks

    switch (txt_file->codec) {
        case CODEC_NONE : uncomp_len = txtfile_read_block_plain (vb, max_bytes); break;
        case CODEC_GZIL : 
        case CODEC_BGZF : uncomp_len = txtfile_read_block_bgz   (vb, max_bytes, uncompress); break; 
        case CODEC_GZ   : uncomp_len = txtfile_read_block_igzip    (vb, max_bytes); break;
        case CODEC_BZ2  : uncomp_len = txtfile_read_block_bz2   (vb, max_bytes); break;

        default: ABORT ("txtfile_read_block: Invalid file type %s (codec=%s)", ft_name (txt_file->type), codec_name (txt_file->codec));
    }

    COPY_TIMER_EVB (read);
    return uncomp_len;
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

        txt_file->header_size_bgzf = bgz_copy_unconsumed_blocks (evb); // copy unconsumed or partially consumed gz_blocks to txt_file->unconsumed_bgz_blocks
    }

    txt_file->txt_data_so_far_single = txt_file->header_size = header_len; 

    biopsy_take (evb);
    
    COPY_TIMER_EVB (txtfile_read_header); // same profiler entry as txtfile_read_header
}

// default "unconsumed" function file formats where we need to read whole \n-ending lines. returns the unconsumed data length
int32_t def_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    ASSERT (*i >= 0 && *i < Ltxt, "*i=%d ∉ [0,%u]", *i, Ltxt);

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

    // case: the data is BGZF-compressed in vb->scratch, except for passed down data from prev VB        
    // uncompress one block at a time to see if its sufficient. usually, one block is enough
    if ((TXT_IS_BGZF || TXT_IS_GZIL) && vb->scratch.len) {

        vb->gzip_compressor = libdeflate_alloc_decompressor (vb, __FUNCLINE);

        for (int block_i=vb->gz_blocks.len32 - 1; block_i >= 0; block_i--) {
            GzBlockZip *bb = B(GzBlockZip, vb->gz_blocks, block_i);
           
            START_TIMER;
            bgz_uncompress_one_block (vb, bb, txt_file->codec);
            COPY_TIMER(bgz_uncompress_during_read);

            // case: we dropped the bb: happens only for the final block in GZIL is truncated, and it was not detected earlier in gzil_read_block.
            if (!bb->is_decompressed) {
                vb->gz_blocks.len32--;
                Ltxt -= bb->txt_size;
                segconf.zip_txt_modified = true;

                WARN ("FYI: %s is truncated - its final GZIL block in incomplete. Dropping final %u bytes of the GZ data.", txt_name, bb->comp_size);
            }

            else  {
                START_TIMER;
                int32_t last_i = Ltxt-1; // test from end of data
                pass_to_next_vb_len = (DT_FUNC(txt_file, unconsumed)(vb, MAX_(bb->txt_index, 0), &last_i)); // note: bb->txt_index might be negative if part of this bb was consumed by the previous VB
                COPY_TIMER (txtfile_get_unconsumed_callback);

                if (pass_to_next_vb_len >= 0) goto done; // we have the answer (callback returns -1 if it needs more data)
            }
        }
    }

    // test remaining txt_data including passed-down data from previous VB
    {
    START_TIMER;
    int32_t last_i = Ltxt-1; // test from end of data
    pass_to_next_vb_len = (DT_FUNC(vb, unconsumed)(vb, 0, &last_i));
    COPY_TIMER (txtfile_get_unconsumed_callback);
    }

    // case: callback doesn't have enough data for even one line, but file has no more data
    if (flag.truncate && pass_to_next_vb_len < 0 && !segconf.running) {
        WARN ("FYI: %s is truncated - its final %s in incomplete. Dropping this partial final %s of %u bytes.", 
              txt_name, DTPT(line_name), DTPT(line_name), Ltxt);
        txt_file->last_truncated_line_len = Ltxt; 
        Ltxt = pass_to_next_vb_len = 0; // truncate last partial line
        segconf.zip_txt_modified = true;
    }

    ASSERT (pass_to_next_vb_len >= 0 ||
            segconf.running, // case: we're testing memory and this VB is too small for a single line - return and caller will try again with a larger VB 
            "Reason: failed to find a full line %sin vb=%s data_type=%s txt_data.len=%u txt_file->codec=%s is_last_vb_in_txt_file=%s interleaved=%s.\n"
            "Known possible causes:\n"
            "- The file is %s %s. Tip: try running with --truncate\n"
            "- The file is not a %s file.\n"
            "VB dumped: %s\n",  
            DTPT(is_binary) ? "" : "(i.e. newline-terminated) ",
            VB_NAME, dt_name (txt_file->data_type), Ltxt, codec_name (txt_file->codec), TF(vb->is_last_vb_in_txt_file), TF(segconf.is_interleaved),
            DTPT(is_binary) ? "truncated but not on the boundary of the" : "missing a newline on the last", DTPT(line_name),
            TXT_DT(REF) ? "FASTA" : dt_name (txt_file->data_type),
            txtfile_dump_vb (vb, txt_name).s);

done:
    if (vb->gzip_compressor)
        libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gzip_compressor, __FUNCLINE);
    
    return (uint32_t)pass_to_next_vb_len;
}

static bool seggable_size_is_modifiable (void)
{
    return !is_read_via_ext_decompressor (txt_file) && !TXT_IS_PLAIN;
}
// estimate the size of the txt_data of the file - i.e. the uncompressed data excluding the header - 
// based on the observed or assumed compression ratio of the source compression so far
static void txtfile_set_seggable_size (void)
{
    uint64_t disk_size = txt_file->disk_size ? txt_file->disk_size 
                       : flag.stdin_size     ? flag.stdin_size // user-provided size
                       :                       0; // our estimate will be 0 
    double source_comp_ratio=1;

    if (!is_read_via_ext_decompressor (txt_file)) {
        if (txt_file->is_remote || txt_file->redirected)
            source_comp_ratio = 4;

        else if (TXT_IS_PLAIN) 
            source_comp_ratio = 1; 

        else {    
            double plain_len = txt_file->txt_data_so_far_single + txt_file->unconsumed_txt.len; //  all data that has been decompressed
            double comp_len  = TXT_IS_BZ2 ? file_tell (txt_file, HARD_FAIL)
                                          : txt_file->disk_so_far - txt_file->gz_data.len; // data read from disk, excluding data still awaiting decompression
            
            // case: header is whole BGZF blocks - remove header from calculation to get a better estimate of the seggable compression ratio
            if (txt_file->header_size_bgzf) { 
                plain_len -= txt_file->header_size;
                comp_len  -= txt_file->header_size_bgzf;
            }

            source_comp_ratio = plain_len / MAX_(comp_len, 1);
        }
    }
        
    // external decompressors
    else switch (txt_file->source_codec) {
        case CODEC_BCF:  source_comp_ratio = 10; break; // note: .bcf files might be compressed or uncompressed - we have no way of knowing as "bcftools view" always serves them to us in plain VCF format. These ratios are assuming the bcf is compressed as it normally is.
        case CODEC_XZ:   source_comp_ratio = 15; break;
        case CODEC_CRAM: source_comp_ratio = 25; break;
        case CODEC_ORA:  source_comp_ratio = 25; break;
        case CODEC_ZIP:  source_comp_ratio = 3;  break;
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

static uint32_t txt_data_alloc_size (uint32_t vb_size)
{
    return vb_size + 
           TXTFILE_READ_VB_PADDING; // we need this quantity of unused bytes at the end of vb.txt_data
}

// ZIP main thread
void txtfile_read_vblock (VBlockP vb)
{
    START_TIMER;

    ASSERTNOTNULL (txt_file);
    ASSERT_DT_FUNC (txt_file, unconsumed);
    ASSERT (IN_RANGE (segconf.vb_size, ABSOLUTE_MIN_VBLOCK_MEMORY, ABSOLUTE_MAX_VBLOCK_MEMORY) || segconf.running,
            "Invalid vb_size=%"PRIu64" comp_i(0-based)=%u", segconf.vb_size, z_file->num_txts_so_far-1);

    if (txt_file->no_more_blocks && !txt_file->unconsumed_txt.len) return; // we're done
    
    uint32_t my_vb_size = segconf.vb_size; // might grow to match a FASTQ R2 vb to its R1 pair

    buf_alloc (vb, &vb->txt_data, 0, txt_data_alloc_size (my_vb_size), char, 1, "txt_data");    

    // read data from the file until either 1. EOF is reached 2. end of vb is reached
    uint32_t pass_to_next_vb_len = 0;

    // start with using the data passed down from the previous VB (note: copy & free and not move! so we can reuse txt_data next vb)
    if (txt_file->unconsumed_txt.len) {
        uint64_t bytes_moved = MIN_(txt_file->unconsumed_txt.len, segconf.vb_size);
        buf_copy (vb, &vb->txt_data, &txt_file->unconsumed_txt, char, 0, bytes_moved, "txt_data");
        buf_remove (txt_file->unconsumed_txt, char, 0, bytes_moved);
    }

    bool is_bgz = TXT_IS_BGZF || TXT_IS_GZIL;

    if (is_bgz) bgz_zip_init_vb (vb); 
    
    vb->comp_i = flag.zip_comp_i;  // needed for VB_NAME

    bool always_uncompress = flag.zip_uncompress_source_during_read || segconf.running;

    // case: compute thread should decompress
    if (!always_uncompress && (TXT_IS_BGZF || TXT_IS_GZIL))
        vb->txt_codec = txt_file->codec;

    uint32_t max_block_size = TXT_IS_BGZF ? BGZF_MAX_BLOCK_SIZE : GZIL_MAX_BLOCK_SIZE;

    for (bool first=true; ; first=false) {        
        uint32_t bytes_requested = MIN_(my_vb_size - Ltxt, 1 GB /* read() can't handle more */);
        bool no_read_expected = is_bgz && (bytes_requested <= max_block_size); // in this case, txtfile_read_block is expected to return 0 

        uint32_t len = (my_vb_size > Ltxt) ? txtfile_read_block (vb, bytes_requested, always_uncompress) : 0;
        
        if (!len && first && !Ltxt) goto done; // case: no data read nor pass up from prev vb (and hence also no data to pass down to next vb)

        // when reading BGZF, we might be filled up even without completely filling my_vb_size 
        // if there is room left for only a partial BGZF block (we can't read partial blocks)
        uint32_t filled_up = my_vb_size - (is_bgz ? (max_block_size - 1) : 0);

        if (len && Ltxt < filled_up) continue;  // continue filling up txt_data...

        // case: this is the 2nd file of a fastq pair - make sure it has at least as many fastq "lines" as the first file
        uint32_t my_lines, pair_num_lines, pair_txt_data_len; 
        VBIType pair_vb_i;
        if (flag.pair == PAIR_R2 &&  // we are reading the second file of a fastq file pair (with --pair)
            !fastq_txtfile_have_enough_lines (vb, &pass_to_next_vb_len, &my_lines, &pair_vb_i, &pair_num_lines, &pair_txt_data_len)) { // we don't yet have all the data we need

            // note: the opposite case where R2 has more reads than R1 is caught in fastq_txtfile_have_enough_lines or zip_prepare_one_vb_for_dispatching
            ASSINP ((len || no_read_expected) && Ltxt, "Error: File %s has less FASTQ reads than its R1 mate (vb=%s has %u lines while its pair_vb_i=%d num_R1_VBs=%u has pair_txt_data_len=%u pair_num_lines=%u; vb=%s Ltxt=%u bytes_requested=%u bytes_read=%u no_more_blocks=%s my_vb_size=%u vb_size=%s src_codec=%s disk_so_far=%"PRIu64").%s", 
                    txt_name, VB_NAME, my_lines, pair_vb_i, sections_get_num_vbs (FQ_COMP_R1), pair_txt_data_len/*only set if flag.debug*/, pair_num_lines, VB_NAME, Ltxt, bytes_requested, len, TF(txt_file->no_more_blocks), my_vb_size, str_size (segconf.vb_size).s, txtfile_codec_name (z_file, vb->comp_i).s, txt_file->disk_so_far,
                    (flag.truncate && (TXT_IS_BGZF || TXT_IS_GZIL || z_file->comp_codec[0] == CODEC_BGZF || z_file->comp_codec[0] == CODEC_GZIL)) ? " Tip: this might due to --truncate. Try adding --no-bgzf" : "");


            // if we need more lines - increase memory and keep on reading
            my_vb_size *= 1.1 * ((double)pair_num_lines / (double)my_lines);

            buf_alloc (vb, &vb->txt_data, 0, txt_data_alloc_size (my_vb_size), char, 1, "txt_data");    
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
            goto done;
        }
    }

    if (pass_to_next_vb_len) {
        // note: we might some unconsumed data, pass it up to the next vb. possibly we still have unconsumed data (can happen if DVCF reject
        // data was passed down from the txt header, greater than my_vb_size)
        buf_insert (evb, txt_file->unconsumed_txt, char, 0, Btxt (Ltxt - pass_to_next_vb_len), pass_to_next_vb_len, "txt_file->unconsumed_txt");
        Ltxt -= pass_to_next_vb_len; 

        // copy unconsumed or partially consumed gz_blocks to txt_file->unconsumed_bgz_blocks
        if (is_bgz) 
            bgz_copy_unconsumed_blocks (vb);
    }

    vb->vb_position_txt_file = txt_file->txt_data_so_far_single;
    vb->is_last_vb_in_txt_file = txt_file->no_more_blocks && !txt_file->unconsumed_txt.len;
    txt_file->txt_data_so_far_single += Ltxt;

    zip_init_vb (vb);

    if (!txt_file->est_seggable_size || seggable_size_is_modifiable())
        txtfile_set_seggable_size();

    if (!segconf.running) {
        biopsy_take (vb);
        dispatcher_increment_progress ("read", txt_file->est_num_lines ? (Ltxt / MAX_(segconf.line_len,1)) : Ltxt);
    }

done:
    if (flag_is_show_vblocks (ZIP_TASK_NAME)) 
        iprintf ("VB_READ(id=%d) vb=%s Ltxt=%u vb_position_txt_file=%"PRIu64" unconsumed_txt.len=%u is_last_vb_in_txt_file=%s\n", 
                 vb->id, VB_NAME, Ltxt, vb->vb_position_txt_file, txt_file->unconsumed_txt.len32, TF(vb->is_last_vb_in_txt_file));

    COPY_TIMER (txtfile_read_vblock);
}

DataType txtfile_zip_get_file_dt (rom filename)
{
    FileType ft = flag.stdin_type; // check for --input option

    if (ft == UNKNOWN_FILE_TYPE) // no --input - get file type from filename
        ft = file_get_type (filename);

    return file_get_data_type_of_input_file (ft);
}

// outputs details on txt_file->codec of a component, as stored in z_file
StrText txtfile_codec_name (FileP z_file/*obscures global*/, CompIType comp_i) 
{
    StrText s;

    if (!IN_RANGE (comp_i, 0, MAX_NUM_COMPS-1))
        snprintf (s.s, sizeof (s.s), "comp_i=%u out_of_range", comp_i);

    else if (z_file->comp_codec[comp_i] == CODEC_BGZF) {
            if (z_file->comp_bgzf[comp_i].level < BGZF_COMP_LEVEL_UNKNOWN)
                snprintf (s.s, sizeof (s.s), "BGZF(%s[%d])", bgzf_library_name (z_file->comp_bgzf[comp_i].library, false), z_file->comp_bgzf[comp_i].level);
            else
                strcpy (s.s, "BGZF(unknown_lib)");
        }
        
    else if (z_file->comp_codec[comp_i]==CODEC_GZ) {
        bool fextra = ((z_file->gz_header[comp_i][3] & 4) == 4); // FEXTRA is bit 2 of FLG
        
        snprintf (s.s, sizeof (s), "GZ(%.24s%.20s%.20s)", 
                  str_to_hex (z_file->gz_header[comp_i], fextra ? 12 : 10).s,
                  cond_str (z_file->gz_isize[comp_i][0], "-", str_size (z_file->gz_isize[comp_i][0]).s),
                  cond_str (z_file->gz_isize[comp_i][1], "-", str_size (z_file->gz_isize[comp_i][1]).s));
    }

    else 
        strcpy (s.s, codec_name (z_file->comp_codec[comp_i]));

    return s;
}