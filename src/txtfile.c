// ------------------------------------------------------------------
//   txtfile.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
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
#include "mgzip.h"
#include "biopsy.h"
#include "zip.h"
#include "arch.h"
#include "dispatcher.h"
#include "threads.h"
#include "zlib/zlib.h"
#include "libdeflate_1.19/libdeflate.h"
#include "bzlib/bzlib.h"
#include "igzip/igzip_lib.h"

#define MAX_TXT_HEADER_LEN ((uint64_t)0xffffffff) // maximum length of txt header - one issue with enlarging it is that we digest it in one go, and the digest module is 32 bit

// PIZ: dump bad vb to disk
StrText txtfile_dump_vb (VBlockP vb, rom base_name, BufferP txt_data)
{
    if (!txt_data) 
        return (StrText){"(txt_data=NULL)"};

    StrText dump_filename;
    
    snprintf (dump_filename.s, sizeof (dump_filename.s), "%u.bad-recon%s", vb->vblock_i, (txt_file ? file_plain_ext_by_dt (txt_file->data_type) : ""));

    if (flag.is_windows) str_replace_letter (dump_filename.s, strlen(dump_filename.s), '/', '\\');

    if (!txt_data) txt_data = &vb->txt_data;
    buf_dump_to_file (dump_filename.s, txt_data, 1, false, false, false, txt_data->len > 8 MB);

    return dump_filename;
}

// returns the requested number of bytes, except if eof in which case it could be less.
uint32_t txtfile_fread (FileP file, 
                        FILE *fp,   // note: non-NULL if different from file->file (when re-reading) 
                        void *addr, // NULL means append to gz_data (reallocing if needed) 
                        int32_t size, int64_t *disk_so_far)
{
    if (size <= 0) return 0;

    if (!fp) fp = (FILE *)file->file;
    bool is_gz_data = (addr == NULL);

    if (is_gz_data) {
        buf_alloc (evb, &file->gz_data, size, 0, uint8_t, 1.1, "txt_file->gz_data");
        addr = BAFT8 (file->gz_data); 
    }

    uint32_t bytes = fread (addr, 1, size, fp);
    if (disk_so_far) *disk_so_far += bytes;   

    int save_errno = errno;
    ASSERT (bytes == size || !ferror (fp) || 
            "Error while reading %s codec=%s on filesystem=%s - requested %u bytes but read only %u: (%u)%s", 
            file->basename, codec_name (file->effective_codec), arch_get_filesystem_type (file).s, size, bytes, save_errno, strerror (save_errno));

    if (is_gz_data) file->gz_data.len32 += bytes;

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
            txt_file->basename, arch_get_filesystem_type (txt_file).s, size, bytes, errno, arch_str_error());

    txt_file->disk_so_far += bytes;
}

static inline uint32_t txtfile_read_block_plain (VBlockP vb, uint32_t max_bytes)
{
    char *data = BAFTtxt;
    int32_t bytes_read;

    // case: we have data passed to us from txtfile_discover_specific_gz - handle it first (possibly txt_data already contains data passed down from previous VBs)
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

// chuck size chosen to be equal to the default disk read-ahead buffer in Linux (eg /sys/block/sda/queue/read_ahead_kb) 
// so that decompression is parallelized with disk read-ahead buffer filling (a bigger buffer would cause the disk to be idle while we are still decompressing)
#define IGZIP_CHUNK (128 KB) 

static StrText display_gz_flags (uint32_t flg)
{
    StrText s={};
    uint32_t s_len=0;

    rom names[] = { "TEXT", "HCRC", "XTRA", "NAME", "CMNT", "RES1", "RES2", "RES3" };

    for (int i=0; i < 8; i++) 
        if ((flg >> i) & 1) 
            SNPRINTF (s, "%s|", names[i]);
    if (s_len) 
        s.s[s_len-1] = 0; // remove final '|'
    else
        s.s[0] = '0';

    return s;
}

static StrText display_gz_xfl (uint8_t xfl)
{
    StrText s={};

    if      (xfl==2) strcpy (s.s, "BEST");
    else if (xfl==4) strcpy (s.s, "FAST");
    else snprintf (s.s, sizeof(s.s), "%u", xfl);

    return s;
}

StrText1K display_gz_header_ex (STR8p(h), bool obscure_fname, uint32_t *out_h_len/*out*/)
{
    StrText1K s = {};
    uint32_t s_len = 0;
    bytes save_h = h;
    #define ADVANCE_h(n) STRinc(h, n)

    if (h_len < 3) goto fail;
    if (memcmp (h, BGZF_PREFIX, 3))
        return (StrText1K){ "<not GZIP>" };

    if (h_len < 10) goto fail;
    uint8_t flg = h[3];
    time_t mtime = GET_UINT32(&h[4]);
    uint8_t xfl = h[8];
    uint8_t os  = h[9];

    struct tm *tm_info = localtime (&mtime);
    char time_str[32];
    if (IN_RANGX (mtime, 1, time(NULL)))
        strftime (time_str, sizeof(time_str), "%Y-%m-%d_%H:%M:%S", tm_info);
    
    // in Ultima, mtime holds the total_txt_size (R1+R2)
    else if (TECH(ULTIMA) && mtime && !IN_RANGX(tm_info->tm_year, 20, 50)) // note: tech not set in initial call, but set later (e.g. when reporting to stats)
        snprintf (time_str, sizeof(time_str), "ULTIMA=%u", (uint32_t)mtime);

    else
        snprintf (time_str, sizeof(time_str), "MTIME=%u", (uint32_t)mtime);

    #define OS(x,s) os==x?s :
    // note: for brevity, we don't display ID and CM: we verify above that they are BGZF_PREFIX.
    SNPRINTF (s, "{ FLG=%s %s XFL=%s OS=%s", display_gz_flags(flg).s, time_str, display_gz_xfl(xfl).s, 
              OS(3,"Unix") OS(7,"Mac") OS(11,"Windows") OS(255,"Unknown") str_int_s(h[9]).s);
    ADVANCE_h(10);

    if (IS_FLAG (flg, 4)) { // FEXTRA
        if (h_len < 2) goto fail;
        uint16_t xlen = GET_UINT16(h);
        ADVANCE_h(2);
        SNPRINTF (s, " XLEN=%u X=[", xlen);

        if (h_len < xlen) goto fail;
        
        while (xlen) {
            if (h_len < 4) goto fail;
            uint32_t f_len = GET_UINT16(&h[2]);
            if (h_len < 4 + f_len) goto fail;

            uint16_t id = GET_UINT16 (h);
            if (xlen == 6 && id == 0x4342 && f_len == 2)
                SNPRINTF (s, " { BGZF BSIZE-1=%-5u }", GET_UINT16(&h[4]));
            
            else if (xlen == 8 && id == 0x4749 && f_len == 4)
                SNPRINTF (s, " { MGZF BSIZE=%u }", GET_UINT32(&h[4]));

            else 
                SNPRINTF (s, " { ID=%02X%02X LEN=%u F=%s }", h[0], h[1], f_len, str_to_hex (&h[4], f_len).s);
            
            ADVANCE_h(f_len + 4);
            xlen -= f_len + 4;
        }

        SNPRINTF0 (s, " ]");
    }

    if (IS_FLAG(flg, 8)) { // FNAME 
        int name_len = strnlen ((rom)h, 1 KB);
        if (h_len < name_len+1 || name_len == 1 KB) goto fail; // somewhat safety 

        if (obscure_fname)
            SNPRINTF0 (s, " NAME=<hidden>");
        else
            SNPRINTF (s, " NAME=\"%s\"", h);
        ADVANCE_h (name_len+1);
    }

    if (IS_FLAG(flg, 16)) { // FCOMMENT 
        int comment_len = strnlen ((rom)h, 1 KB);
        if (h_len < comment_len+1 || comment_len == 1 KB) goto fail; // somewhat safety 

        SNPRINTF (s, " CMNT=\"%s\"", h);
        ADVANCE_h (comment_len+1);
    }

    if (IS_FLAG(flg, 1)) { // FHCRC
        if (h_len < 2) goto fail;
        SNPRINTF (s, " HCRC=%u", GET_UINT16(h));
        ADVANCE_h(2);
    } 

    SNPRINTF0 (s, " }");

    if (flg == 8/*name*/ && mtime && xfl == 0 && os == 3)  SNPRINTF0 (s, " (gzip/pigz)");

    if (out_h_len) *out_h_len = h - save_h;

    return s;

fail:
    return (StrText1K){ "<not enough data>" };
}

StrText1K display_gz_header (STR8p(h), bool obscure_fname)
{
    return display_gz_header_ex (STRa(h), obscure_fname, NULL);
}

uint32_t gzip_header_length (FileP file)
{
    ARRAY (char, h, file->gz_data);
    if (h_len < 10) return 0;
    
    if (h[0] != '\x1f' || h[1] != '\x8b' || h[2] != '\x08') return 0; // not gzip

    uint8_t flg = h[3];
    uint32_t xlen=0, name_len=0, comment_len=0, hcrc_len = 0;
    
    if (IS_FLAG (flg, 4)) { // FEXTRA
        if (h_len < 12) return 0;
        xlen = 2 + GET_UINT16(h+10);
        if (h_len < 10 + xlen) return 0;
    }

    if (IS_FLAG (flg, 8)) { // FNAME
        int offset = 10 + xlen;
        name_len = strnlen (h + offset, h_len - offset - 1);
        if (h[offset + name_len]) return 0; // expecting \0
        name_len++; // inc. the \0
    }

    if (IS_FLAG (flg, 16)) { // FCOMMENT
        int offset = 10 + xlen + name_len;
        comment_len = strnlen (h + offset, h_len - offset - 1);
        if (h[offset + comment_len]) return 0; // expecting \0
        comment_len++; // inc. the \0
    }

    if (IS_FLAG(flg, 1)) { // FHCRC
        int offset = 10 + xlen + name_len + comment_len;
        if (h_len < offset + 2) return 0;
        hcrc_len = 2;
    }

    return 10 + xlen + name_len + comment_len + hcrc_len;
}

void txtfile_initialize_igzip (FileP file)
{
    ASSERTNOTINUSE (file->igzip_state);

    buf_alloc_exact_zero (evb, file->igzip_state, 1, struct inflate_state, "txt_file->igzip_state");
    struct inflate_state *state = B1ST (struct inflate_state, file->igzip_state);

    isal_inflate_init (state);
    state->crc_flag = ISAL_GZIP;
}

bool txtfile_is_gzip (FileP file)
{
    txtfile_fread (file, NULL, NULL, 3, &file->disk_so_far);
    
    ARRAY (char, h, file->gz_data);
    return h[0] == '\x1f' && h[1] == '\x8b' && h[2] == '\x08';
}

// check if the first 3 gz blocks have the same isize. 
static bool txtfile_segconf_discover_constant_isize (STR8p(header))
{
    #define MAX_CONSTANT_SIZE_HEADER_LEN 255 // must be large enough for all the constant-size codecs we recognize
    if (MAX_CONSTANT_SIZE_HEADER_LEN > 255) return false; 

    // we search for the isize of the first gz block (which is its last field) followed by the constant header of a new block
    uint8_t signature[4 + header_len];
    PUT_UINT32 (signature, txt_file->gz_data.uncomp_len); //  this is the fixed-isize validated in mgzip_read_block_no_bsize 
    memcpy (signature + 4, header, header_len);

    // find 3 more instances of isize⁀header in the data with same isize
    bytes next_gz_block = B8(txt_file->gz_data, txt_file->gz_data.comp_len); // points to start of 2nd gz block in file
    bytes next  = next_gz_block - 4; // point to isize of first gz block
    bytes after = BAFT8(txt_file->gz_data);

    for (int i=0 ; i < 3; i++) {
        bytes next_isize;
        // can't find another isize⁀header
        if (!(next_isize = memmem (next, after - next, signature, 4 + header_len)))

            // cannot find another isize⁀header: this is ok iff we read all the data of the file, and the reason memmem 
            // can't find another isize⁀header is that there are no more headers i.e. not that the next header is preceded by an invalid isize
            return feof ((FILE *)txt_file->file) &&      // we read the entire file
                   after - next_gz_block <= txt_file->gz_data.uncomp_len &&  // the remaining data is less than uncomp_len (i.e. the "standard" isize) - used an upper-bound for comp_len
                   !memmem (next, after - next, STRa(header)); // there are no more headers (note: a rare case of false positive identification of header will cause this file not to be exactable - that's ok)

        next_gz_block = next_isize + 4;
        next          = next_isize + sizeof (signature); // continue search from after the signature
    }

    return true; // there are at least 4 gz blocks, and the first 3 have the same isize
}

// run when opening a txt file (before reading the txt_header), except for FASTQ for which it is 
// run at the end of segconf after segconf.tech is known.
void txtfile_discover_specific_gz (FileP file)
{   
    START_TIMER;
    bool keep_src_codec = (file->src_codec == CODEC_CRAM || file->src_codec == CODEC_BAM || file->src_codec == CODEC_BCF);
        
    // read the first potential BGZF block to test if this is GZ or BGZF 
    // note: we read if even if --no-bgzf to capture the data for z_file->gz_header
    GzStatus status = mgzip_read_block_with_bsize (file, true, CODEC_BGZF);

    // copy GZ header: data should be in gz_data 
    file->gz_header_len = gzip_header_length (file); 
    memcpy (file->gz_header, B1ST8 (file->gz_data), MIN_(MAX_GZ_HEADER_LEN, file->gz_header_len));

    // case: this is a BGZF block
    // note: we keep the still-compressed data in vb->comp_txt_data for later consumption
    if (!flag.no_bgzf && status == GZ_SUCCESS && file->gz_data.uncomp_len > 0) {
        if (!keep_src_codec) file->src_codec = CODEC_BGZF;
        file->effective_codec = CODEC_BGZF;
    }

    // for regulars files, we already skipped 0 size files. This can happen in STDIN
    else if (status == GZ_SUCCESS && file->gz_data.uncomp_len == 0) {
        ASSINP (!flags_pipe_in_process_died(), // only works for Linux
                "Pipe-in process %s (pid=%u) died without sending any data",
                flags_pipe_in_process_name(), flags_pipe_in_pid());

        ABORTINP ("No data exists in input file %s", file->name ? file->name : FILENAME_STDIN);
    }
    
    // case: this is non-BGZF GZIP format: test for one of the FASTQ codec. Run after segconf.
    else if (status == GZ_IS_OTHER_FORMAT && file->data_type == DT_FASTQ) {      
        #define SET_CODEC(codec) ({ file->src_codec = CODEC_##codec; \
                                    if (!flag.no_bgzf) file->effective_codec = CODEC_##codec; })

        if (TECH(ILLUMINA) && mgzip_read_block_no_bsize (file, true, CODEC_IL1M) == GZ_SUCCESS
        &&  txtfile_segconf_discover_constant_isize (_8(ILxM_PREFIX)))  // lightweight
            SET_CODEC(IL1M);
        
        else if (TECH(ILLUMINA) && mgzip_read_block_no_bsize (file, true, CODEC_IL4M) == GZ_SUCCESS
        &&  txtfile_segconf_discover_constant_isize (_8(ILxM_PREFIX)))           
            SET_CODEC(IL4M);
        
        else if (TECH(MGI) && mgzip_read_block_no_bsize (file, true, CODEC_MGSP) == GZ_SUCCESS &&
            txtfile_segconf_discover_constant_isize (_8(MGSP_HEADER))) {
            SET_CODEC(MGSP);
            file->num_mgsp_blocks_in_vb = 0; // reset
        }

        else if (TECH(MGI) && mgzip_read_block_with_bsize (file, true, CODEC_MGZF) == GZ_SUCCESS) 
            SET_CODEC(MGZF);

        else if (TECH(ELEMENT) && mgzip_read_block_no_bsize (file, true, CODEC_EMFL) == GZ_SUCCESS &&
            txtfile_segconf_discover_constant_isize (STRa(file->gz_header))) 
            SET_CODEC(EMFL);

        else if (TECH(ELEMENT) && mgzip_read_block_no_bsize (file, true, CODEC_EMVL) == GZ_SUCCESS &&
            str_issame_(B1STc(file->gz_data), file->gz_data.comp_len, _S(EMVL_FIRST_BLOCK)))  // first block is an empty block
            SET_CODEC(EMVL);

        // blocked-GZ of unknown type
        else if (mgzip_read_block_no_bsize (file, true, CODEC_GZBL) == GZ_SUCCESS)
            SET_CODEC (GZBL);

        else 
            goto generic_gz;
        #undef SET_CODEC
    } 

    // unblocked GZ, or blocked GZ with a header longer than 10 bytes
    else if (status != GZ_NOT_GZIP) generic_gz: {
        if (!keep_src_codec)file->src_codec = CODEC_GZ;
        file->effective_codec = CODEC_GZ;    
        
        mgzip_set_is_exactable (file, false, "Not an MGZIP codec");
    }

    // case: this is not GZIP format at all. treat as a plain file, and put the data read in vb->comp_txt_data 
    // for later consumption is txtfile_read_block_plain
    else {
        #define BZ2_MAGIC "BZh"
        #define XZ_MAGIC  (char[]){ 0xFD, '7', 'z', 'X', 'Z', 0 }
        #define ZIP_MAGIC (char[]){ 0x50, 0x4b, 0x03, 0x04 }
        #define ORA_MAGIC (char[]){ 0x49, 0x7c } // https://support-docs.illumina.com/SW/ORA_Format_Specification/Content/SW/ORA/ORAFormatSpecification.htm
        #define CRAM_MAGIC "CRAM"  // first 4 bytes of a CRAM file

        // we already open the file, so not easy to re-open with BZ2_bzopen as it would require injecting the read data into the BZ2 buffers
        if (str_isprefix_(STRb(file->gz_data), BZ2_MAGIC, 3)) 
            ABORTINP0 ("The data seems to be in bz2 format. "_TIP"Use --input to specify the type (eg: \"genozip --input sam.bz2\")");

        else if (str_isprefix_(STRb(file->gz_data), XZ_MAGIC, 6)) {
            if (file->redirected) ABORTINP0 ("Compressing piped-in data in xz format is not currently supported");
            if (file->is_remote) ABORTINP0 ("The data seems to be in xz format. "_TIP"Use --input to specify the type (eg: \"genozip --input sam.xz\")");
            ABORTINP0 ("The data seems to be in xz format. Please use --input to specify the type (eg: \"genozip --input sam.xz\")");
        }

        else if (str_isprefix_(STRb(file->gz_data), ZIP_MAGIC, 4)) {
            if (file->redirected) ABORTINP0 ("Compressing piped-in data in zip format is not currently supported");
            if (file->is_remote) ABORTINP0 ("The data seems to be in zip format. "_TIP"Use --input to specify the type (eg: \"genozip --input generic.zip\")");
            ABORTINP0 ("The data seems to be in zip format. Please use --input to specify the type (eg: \"genozip --input generic.zip\")");
        }

        else if (str_isprefix_(STRb(file->gz_data), ORA_MAGIC, 2)) {
            if (file->redirected) ABORTINP0 ("Compressing piped-in data in ora format is not currently supported");
            if (file->is_remote) ABORTINP0 ("The data seems to be in ora format. "_TIP"Use '--input fastq.ora'");
            ABORTINP0 ("The data seems to be in ora format. "_TIP" use '--input fastq.ora'");
        }

        // case: a CRAM file without a .cram file name extension
        else if (str_isprefix_(STRb(file->gz_data), CRAM_MAGIC, 4)) {
            if (file->is_remote || file->redirected) ABORTINP0 ("The data seems to be in CRAM format. "_TIP"Use '--input cram'");
            file->src_codec   = file->effective_codec = CODEC_CRAM;
            file->type        = CRAM;
            file->data_type   = DT_BAM;
            file->disk_so_far = 0;
            buf_free (file->gz_data); // we will re-read through external decompressor
            return;
        }

        if (!keep_src_codec) file->src_codec = CODEC_NONE;
        file->effective_codec = CODEC_NONE;
    }

    if (flag.make_reference && IS_MGZIP(file->src_codec))
        file->effective_codec = CODEC_GZ;
    
    // case: R2 codec differ than R1's - just use GZ (decompressing in main thread). This is a little
    // blunt as most codecs will interact well in most cases (e.g. if vb_size is large enough), but there
    // are many edge cases that need to be handled, and this is not a common use case worth dealing with.
    if (IS_R2 && file->effective_codec != z_file->comp_eff_codec[flag.zip_comp_i-1])
        file->effective_codec = CODEC_GZ; 

    if (file->effective_codec == CODEC_GZ)
        txtfile_initialize_igzip (file);

    else if (IS_MGZIP(file->effective_codec)) { 
        file->max_mgzip_isize = file->gz_data.uncomp_len; // note: will be 0 for EMVL, bc first block is 0

        bgzf_initialize_discovery (file); // to discover library and level
    }

    COPY_TIMER_EVB (txtfile_discover_specific_gz);
}

// ZIP main thread: called after txt and z are open, and txt codecs have been discovered.
void txtfile_zip_finalize_codecs (void)
{
    ASSERTNOTNULL (z_file);
    ASSERTNOTNULL (txt_file);
    
    if (flag.zip_comp_i < MAX_NUM_COMPS) { // for stats
        z_file->comp_src_codec[flag.zip_comp_i] = txt_file->src_codec;
        z_file->comp_eff_codec[flag.zip_comp_i] = txt_file->effective_codec;
    }

    memcpy (z_file->comp_gz_header[flag.zip_comp_i], txt_file->gz_header, MAX_GZ_HEADER_LEN);

    if (flag_show_bgzf) 
        iprintf ("%s: src_codec=%s effective_codec=%s gz_header=%s\n", txt_file->basename, // same format as in txtfile_zip_finalize_codecs
                 codec_name (txt_file->src_codec), codec_name (txt_file->effective_codec),
                 display_gz_header (z_file->comp_gz_header[flag.zip_comp_i], MAX_GZ_HEADER_LEN, false).s);
}

// runs in main thread, reads and uncompressed GZ, and populates txt_data for vb
static uint32_t txtfile_read_block_igzip (VBlockP vb, uint32_t max_bytes, bool *is_data_read) 
{
    START_TIMER;
    ASSERTISALLOCED (txt_file->gz_data);

    struct inflate_state *state = B1ST (struct inflate_state, txt_file->igzip_state);
    sSTRl (last_gz_header, MAX_GZ_HEADER_LEN)*0;

    // save in case of realloc in txtfile_fread
    uint32_t next_in_before = state->next_in ? BNUM(txt_file->gz_data, state->next_in) : 0; 
    uint64_t txt_index = vb->txt_data.len;

    // top up gz_data
    if (txt_file->gz_data.len32 < IGZIP_CHUNK) 
        txtfile_fread (txt_file, NULL, NULL, (int32_t)IGZIP_CHUNK - (int32_t)txt_file->gz_data.len32, &txt_file->disk_so_far);

    { START_TIMER
    state->next_in   = B8(txt_file->gz_data, next_in_before);
    state->avail_in  = BAFT8(txt_file->gz_data) - state->next_in;
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

    if (state->block_state == ISAL_BLOCK_NEW_HDR && flag_show_bgzf) {
        last_gz_header_len = MIN_(MAX_GZ_HEADER_LEN, state->avail_in);
        memcpy (last_gz_header, state->next_in, last_gz_header_len);
    }

    int ret = isal_inflate (state); // new gzip header in a file that has concatented gzip compressions
    ASSERT (ret == ISAL_DECOMP_OK || ret == ISAL_END_INPUT, "isal_inflate error: %s avail_in=%u avail_out=%u",
            isal_error (ret), txt_file->gz_data.len32, max_bytes);

    COPY_TIMER(igzip_uncompress_during_read); }

    uint32_t gz_data_consumed = BNUM (txt_file->gz_data, state->next_in);
    
    if (state->block_state == ISAL_BLOCK_FINISH && z_file && gz_data_consumed >= 4) {
        
        // note: this is actually (isize % 2^32) bc gzip isize field is 32bit
        uint32_t isize = GET_UINT32 (B8(txt_file->gz_data, gz_data_consumed - 4)); 

        bool is_single_block_gz = (!txt_file->max_mgzip_isize && feof((FILE *)txt_file->file));

        // for stats: save the isize from the gzip footer (of the first 2 gzip blocks) 
        if (!is_single_block_gz)
            for (int i=0; i <= 1; i++)
                if (!z_file->gz_isize[flag.zip_comp_i][i]) {
                    z_file->gz_isize[flag.zip_comp_i][i] = isize;
                    break;
                }

        uint64_t decompressed_so_far = txt_file->disk_so_far - txt_file->gz_data.len + BNUM (txt_file->gz_data, state->next_in);

        if (flag_show_bgzf) {
            uint64_t bsize = decompressed_so_far - txt_file->start_gz_block;
            iprintf ("UNCOMPRESS GZ   thread=MAIN    vb=%-7s block_i=%-16"PRIu64" comp_index=%-9"PRIu64" comp_len=%-8"PRIu64" txt_index=%-9"PRIu64" txt_len=%-8u%s gz_header=%s\n", 
                     VB_NAME, txt_file->gz_blocks_so_far, txt_file->start_gz_block, bsize, 
                     txt_index, isize, (bsize > 500 MB && (double) isize / (double) bsize < 0.98) ? "(overflow)" : "", // note: this doesn't catch all overflows, for example it won't catch: bsize=1GB and isize=3GB where the true isize is 7GB
                     display_gz_header ((bytes)STRa(last_gz_header), false).s);
        }

        txt_file->gz_blocks_so_far++;
        txt_file->start_gz_block = decompressed_so_far;

        if (!is_single_block_gz)
            MAXIMIZE (txt_file->max_mgzip_isize, isize);
    }

    *is_data_read = BNUM(txt_file->gz_data, state->next_in) > next_in_before;

    // case FASTQ: don't remove segconf data as we will re-uncompress it after effective_codec discovery 
    if (!txt_file->discover_during_segconf) {
        buf_remove (txt_file->gz_data, char, 0, gz_data_consumed);
        state->next_in = NULL;
        
        inc_disk_gz_uncomp_or_trunc (txt_file, gz_data_consumed);

        txt_file->no_more_blocks = (!state->avail_in && feof ((FILE *)txt_file->file));
    }

    else
        segconf.gz_comp_size = gz_data_consumed; // segconf, for use of txtfile_set_seggable_size 

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

static noreturn void txtfile_dump_comp_txt_data (VBlockP vb, uint32_t this_block_start, FUNCLINE)
{
    char dump_fn[strlen(txt_name)+100];
    snprintf (dump_fn, sizeof (dump_fn), "%s.vb-%u.bad-%s.bad-offset-0x%x", 
              txt_name, vb->vblock_i, codec_name (txt_file->effective_codec), this_block_start);
    
    buf_dump_to_file (dump_fn, &vb->comp_txt_data, 1, false, false, true, false);

    ABORT ("called from %s:%u: %s: Invalid %s block in: block_comp_len=%u. Entire data of this vblock dumped to %s, bad block stats at offset 0x%x",
           func, code_line,  VB_NAME, codec_name (txt_file->effective_codec), txt_file->gz_data.comp_len, dump_fn, this_block_start);
}

// ZIP main thread: add one gz_block to mgzip_isizes and mgzip_starts
static void txtfile_update_isizes (void)
{
    // since 15.0.63: EOF blocks are now in mgzip_isizes and instead set, and entries are "uint32_t isize". 
    // before:        entries were "uint16_t isize", txt_file->mgzip_flags.has_eof_block used to indicate an a file has an EOF block, incorrectly assuming that an EOF block can only occur at the end of the file
    // add isize to buffer that will be written to SEC_GZ_ISIZES
    if (!txt_file->mgzip_isizes.len32) { 
        uint64_t est_n_blocks = (txt_file->gz_data.comp_len > 10000) ? ((double)txt_file->disk_size / (double)txt_file->gz_data.comp_len * 1.1) : 0; // >10000 to avoid over-allocating due to a randomly small block
        buf_alloc (evb, &txt_file->mgzip_isizes,  0, MAX_(10000, est_n_blocks), uint32_t, 0, "txt_file->mgzip_isizes"); 
        buf_alloc (evb, &txt_file->mgzip_digests, 0, MAX_(10000, est_n_blocks), uint32_t, 0, "txt_file->mgzip_digests"); 
        buf_alloc (evb, &txt_file->mgzip_starts,  0, MAX_(10000, est_n_blocks), uint64_t, 0, "txt_file->mgzip_starts");
    }

    buf_append_one (txt_file->mgzip_isizes, txt_file->gz_data.uncomp_len); // add one isize
    buf_append_one (txt_file->mgzip_starts, txt_file->disk_so_far - txt_file->gz_data.len);
    // note: digests are added below ⇩ 
}

// ZIP main thread: called after VB is done uncompressing all gz blocks (during which the digest is calculated)
// might be called out-of-order in respect to VB order
void txtfile_update_gz_digests (VBlock𐤐 vb)
{
    if (!txt_file->is_exactable || 
        !vb->gz_blocks.len/*e.g. generated component VBs*/) return;

    // make sure the pre-allocated array is long enough
    MAXIMIZE (txt_file->mgzip_digests.len, vb->vb_mgzip_i + vb->gz_blocks.len);
    buf_alloc (evb, &txt_file->mgzip_digests, 0, txt_file->mgzip_digests.len, uint32_t, 0, "txt_file->mgzip_digests"); 
    
    // copy gz block digests from vb to txt_file
    uint32_t *restrict txt_digests = B32(txt_file->mgzip_digests, vb->vb_mgzip_i);
    for_buf2𐤐 (GzBlockZip, bb, bb_i, vb->gz_blocks)
        txt_digests[bb_i] = bb->gz_digest;
}

// ZIP: Multi-block GZIP (MGZIP): we read gz data into vb->comp_txt_data - that will be decompressed now or later, depending on "uncompress". 
// We read data with a *decompressed* size up to max_uncomp. vb->comp_txt_data always contains only full MGZIP blocks
static inline uint32_t txtfile_read_block_mgzip (VBlockP vb, 
                                                 int32_t requested_bytes, // 0 means read exactly one mgzip block
                                                 bool uncompress, 
                                                 bool *is_data_read)
{
    START_TIMER;

    uint32_t this_uncomp_len=0;
    *is_data_read = false; // initialize

    if (uncompress) {
        ASSERTISNULL (vb->gz_inflate_mem);
        vb->gz_inflate_mem = libdeflate_alloc_decompressor(vb, __FUNCLINE);
    }

    int64_t start_uncomp_len = vb->comp_txt_data.uncomp_len;
    int32_t max_block_size = mgzip_get_max_block_size();
   
    // comp_txt_data contains gz-compressed data; we use .uncomp_len to track its uncompress length
    buf_alloc (vb, &vb->comp_txt_data, 0, requested_bytes / 4/*typical gz ratio is 4.5-6*/, char, 0, "comp_txt_data");

    bool end_of_vb = false; // mostly MGSP uses this mechanism at the momemt
    while (!txt_file->no_more_blocks) {

        bool have_enough_data =
            ( end_of_vb               ? true // the last block read belongs to the next VB (decided by *_is_valid_isize)
            : (requested_bytes == 0)  ? true // no data requested, so we certainly have enough :)
            : TXT_IS(MGSP)            ? (vb->comp_txt_data.uncomp_len && !txt_file->num_mgsp_blocks_in_vb) // MGSP: we have reached the end of the VB and hence MGSP_is_valid_isize has set num_mgsp_blocks_in_vb to 0
            : TXT_IS_VB_SIZE_BY_BLOCK ? (vb->comp_txt_data.uncomp_len > 0) // one gz-block per VB (+ empty blocks) is all we need if TXT_IS_VB_SIZE_BY_BLOCK
            : /*requested_bytes > 0*/   (vb->comp_txt_data.uncomp_len - start_uncomp_len > requested_bytes - max_block_size));

        // we break if we have enough data, but only after we have verified that the next block is not empty (if it is, we don't break so we can include it in this VB)
        if (have_enough_data           && 
            txt_file->gz_data.comp_len && // there is a next block waiting in gz_data...
            txt_file->gz_data.uncomp_len) // ...and it is not an "empty block" (eg EOF)
            break;

        // read more gz data from disk into txt_file->gz_data
        GzStatus status = txt_file->gz_data.comp_len ? GZ_SUCCESS // we have a block waiting for already
                        : TXT_GZ_HEADER_HAS_BSIZE    ? mgzip_read_block_with_bsize (txt_file, false, txt_file->effective_codec)
                        :                              mgzip_read_block_no_bsize   (txt_file, false, txt_file->effective_codec);

        uint32_t this_block_start = vb->comp_txt_data.len32;
        
        // case: this block will be the first block of the next vb. In the mean time, we leave it in txt_file->gz_data
        if (status == GZ_SUCCESS_END_OF_VB)
            end_of_vb = true; 

        // check for corrupt data - at this point we've already confirm the file's codec so not expecting it to change
        else if (status != GZ_SUCCESS) {
            buf_add_more (vb, &vb->comp_txt_data, txt_file->gz_data.data, txt_file->gz_data.comp_len, "comp_txt_data");
            txtfile_dump_comp_txt_data (vb, this_block_start, __FUNCLINE);
        }

        // add block to list: if we need more data OR if the block is empty (e.g. EOF block)
        else if (txt_file->gz_data.comp_len && // this is 0 at the end of the file or if truncated
            (!have_enough_data || txt_file->gz_data.uncomp_len == 0)) {

            buf_add_more (vb, &vb->comp_txt_data, txt_file->gz_data.data, txt_file->gz_data.comp_len, "comp_txt_data");

            buf_alloc (vb, &vb->gz_blocks, 1, MAX_(1000, txt_file->max_mgzip_isize ? 1.2 * requested_bytes / txt_file->max_mgzip_isize : 0), GzBlockZip, 2, "gz_blocks");
            
            BNXT (GzBlockZip, vb->gz_blocks) = (GzBlockZip)
                { .txt_index  = Ltxt,  // after passed-down data and all previous blocks
                  .gz_index   = this_block_start,
                  .txt_size   = txt_file->gz_data.uncomp_len,
                  .gz_size    = txt_file->gz_data.comp_len,
                  .is_eof     = txt_file->no_more_blocks }; 

            // we add the isize at the same as adding the block to vb->gz_blocks. this is so we can rely on
            // txt_file->mgzip_isizes.len to count the number blocks added to VBs (for calculating vb_mgzip_i)
            txtfile_update_isizes();

            *is_data_read = true;

            if (flag_show_bgzf) {
                GzBlockZip *bb = BLST (GzBlockZip, vb->gz_blocks);
                iprintf ("READ       %-4s thread=MAIN   %-11s block_i=%-7"PRIu64" bb_i=%-3u comp_index=%-9u comp_len=%-8u txt_index=%-9u txt_len=%-8u eof=%-5s%s disk_so_far=%"PRIu64"\n",
                         codec_name (txt_file->effective_codec), 
                         cond_str (vb->vblock_i, " vb=", VB_NAME),
                         txt_file->gz_blocks_so_far, 
                         BNUM (vb->gz_blocks, bb), bb->gz_index, bb->gz_size, bb->txt_index, bb->txt_size, TF(bb->is_eof),
                         str_issame_(Bc(vb->comp_txt_data, this_block_start), txt_file->gz_data.comp_len, _S(BGZF_EOF)) ? " BGZF_EOF  " : "",
                         txt_file->disk_so_far);
            }
            txt_file->gz_blocks_so_far++; // counts blocks in entire file (note: bb_i counts within VB)

            MAXIMIZE (txt_file->max_mgzip_isize, txt_file->gz_data.uncomp_len);

            this_uncomp_len              += txt_file->gz_data.uncomp_len; // total uncompressed length of data read by this function call
            vb->comp_txt_data.uncomp_len += txt_file->gz_data.uncomp_len; // total uncompressed length of data in vb->compress
            Ltxt                         += txt_file->gz_data.uncomp_len; // total length of txt_data after adding decompressed vb->comp_txt_data (may also include pass-down data)

            // we uncompress one block a time in the loop so that the decompression is parallel with the disk read-ahead into cache
            if (uncompress) {
                START_TIMER;
                mgzip_uncompress_one_block (vb, BLST (GzBlockZip, vb->gz_blocks), txt_file->effective_codec);  
                COPY_TIMER(mgzip_uncompress_during_read); 
            }

            // remove the first MGZIP block from the gz_data now that its added it to vb->comp_txt_data and vb->gz_blocks
            buf_remove (txt_file->gz_data, uint8_t, 0, txt_file->gz_data.comp_len);
            txt_file->gz_data.comp_len = txt_file->gz_data.uncomp_len = 0; // note: these refer to the first block, gz_data.len might still be >0        
        }
        
        // case: no more data in the file
        else if (!txt_file->gz_data.len) {
            // previous block is eof. note: we miss it if previous block was in the previous VB
            if (vb->gz_blocks.len) 
                BLST (GzBlockZip, vb->gz_blocks)->is_eof = true;
            
            txt_file->no_more_blocks = true; // EOF without EOF block (usually set it in mgzip_read_block_[with|no]_bsize )
        }
    }

    if (uncompress) {
        buf_free (vb->comp_txt_data); 
        libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gz_inflate_mem, __FUNCLINE); // also sets gz_inflate_mem to NULL
    }

    COPY_TIMER (txtfile_read_block_mgzip);
    return this_uncomp_len;
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
                    "Error: %s is missing a %s header. "_TIP"Use --input=generic to compress as a generic file.", 
                    txt_name, dt_name (txt_file->data_type));

            return i; // return header length
        }
        prev_char = header[i];
    }

    // case: the entire file is just a header
    if (is_eof && prev_char == '\n') 
        return header_len;

    return -1; // not end of header yet
}

static uint32_t txtfile_read_block (VBlockP vb, uint32_t bytes_requested, bool uncompress, bool *is_data_read);

// ZIP main thread: reads txt header into evb->txt_data
void txtfile_read_header (bool is_first_txt)
{
    START_TIMER;

    ASSERT_DT_FUNC (txt_file, is_header_done);

    int32_t header_len;
    bool is_data_read = true;

    evb->vb_mgzip_i = 0; // header always starts with the first gz_block

    // read data from the file until either 1. EOF is reached 2. end of txt header is reached
    #define HEADER_BLOCK (256 KB) // we have no idea how big the header will be... read this much at a time
    while ((header_len = (DT_FUNC (txt_file, is_header_done)(txt_file->no_more_blocks))) < 0) { // we might have data here from txtfile_test_data
        
        if (!is_data_read) {
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
           txtfile_read_block (evb, HEADER_BLOCK, true, &is_data_read);
    }

    // the excess data is for the next vb to read 
    if (evb->txt_data.len > header_len) { 
        buf_copy (evb, &txt_file->unconsumed_txt, &evb->txt_data, char, header_len, 0, "txt_file->unconsumed_txt");
        evb->txt_data.len = header_len; // trim to uncompressed length of txt header

        txt_file->header_size_bgzf = mgzip_copy_unconsumed_blocks (evb); // copy unconsumed or partially consumed gz_blocks to txt_file->unconsumed_mgzip_blocks
    }

    txt_file->txt_data_so_far_single = txt_file->header_size = header_len; 

    txtfile_update_gz_digests (evb);

    biopsy_take (evb);
    
    COPY_TIMER_EVB (txtfile_read_header); // same profiler entry as txtfile_read_header
}

// default "unconsumed" function file formats where we need to read whole \n-ending lines. returns the unconsumed data length
int32_t def_unconsumed (VBlockP vb, uint32_t first_i)
{
    ASSERTNOTZERO (Ltxt);

    int32_t j; for (j=Ltxt-1; j >= (int32_t)first_i; j--) // use j - automatic var - for speed 
        if (*Btxt (j) == '\n') 
            return Ltxt -1 - j;

    return UNCONSUMED_NEED_MORE_DATA; // cannot find \n in the data starting first_i
}

static void txt_file_truncate_final_bytes (VBlockP vb, int32_t *n_bytes)
{
    bool is_R2_missing_R1 = (*n_bytes == -2);

    // case: VB consists of a single truncated line, or it is a R2 VB without an R1 counterpart
    if (*n_bytes < 0) *n_bytes = Ltxt; // entire VB 

    if (!is_R2_missing_R1)
        WARN (_FYI "File is truncated - its final %s in incomplete. Dropping this partial final %s of %u bytes.", 
              DTPT(line_name), DTPT(line_name), *n_bytes);
    else
        WARN (_FYI "This is an R2 file that is longer than its R1 counterpart: dropping %u bytes. (vb=%s)", *n_bytes, VB_NAME);

    txt_file->last_truncated_line_len = *n_bytes; 
    Ltxt -= *n_bytes; 
    *n_bytes = 0; // truncate last partial line
    segconf.zip_txt_modified = true;
    mgzip_set_is_exactable (txt_file, false, "File is truncated");
}

static bool txtfile_get_unconsumed_to_pass_to_next_vb (VBlockP vb, bool *R2_vb_truncated_away)
{
    START_TIMER;

    int32_t final_unconsumed_len = 0;
    
    // case: the data is multiblock-gzip-compressed in vb->comp_txt_data, except for passed down data from prev VB        
    // uncompress one block at a time to see if its sufficient. usually, one block is enough
    if (TXT_IS_MGZIP && vb->comp_txt_data.len) {
        ASSERTISNULL (vb->gz_inflate_mem);
        vb->gz_inflate_mem = libdeflate_alloc_decompressor (vb, __FUNCLINE);

        for_buf_back (GzBlockZip, bb, vb->gz_blocks) {
            START_TIMER;
            mgzip_uncompress_one_block (vb, bb, txt_file->effective_codec);
            COPY_TIMER(mgzip_uncompress_during_read);

            // case: we dropped the bb: happens only for the final block in ILxM is truncated, and it was not detected earlier in ILxM_is_valid_isize.
            if (!bb->is_uncompressed) {
                vb->gz_blocks.len32--;
                Ltxt -= bb->txt_size;
                segconf.zip_txt_modified = true;
                mgzip_set_is_exactable (txt_file, false, "ILxM File is truncated");
                
                ASSERT0 ((TXT_IS(IL1M) || TXT_IS(IL4M)) && flag.truncate, "not expecting to be here");

                WARN (_FYI "%s is truncated - its final %s block in incomplete. Dropping final %u bytes of the GZ data.", 
                      txt_name, codec_name (txt_file->effective_codec), bb->gz_size);
            }

            else  {
                START_TIMER;
                final_unconsumed_len = (DT_FUNC(txt_file, unconsumed)(vb, MAX_(bb->txt_index, 0))); // note: bb->txt_index might be negative if part of this bb was consumed by the previous VB
                COPY_TIMER (txtfile_get_unconsumed_callback);

                if (final_unconsumed_len >= 0) {
                    if (final_unconsumed_len && flag.truncate && txt_file->no_more_blocks) 
                        txt_file_truncate_final_bytes (vb, &final_unconsumed_len);

                    goto done; // we have the answer (callback returns -1 if it needs more data)
                }
            }
        }
    }

    // case: full line not detected in all gz-compressed data: test remaining txt_data including passed-down data from previous VB
    // case: codec is not BGZIP (i.e. it is NONE, GZ or BZ2) and hence already fully uncompressed
    {
    START_TIMER;
    final_unconsumed_len = (DT_FUNC(vb, unconsumed)(vb, 0));
    COPY_TIMER (txtfile_get_unconsumed_callback);
    }

    // case: truncate entire VB (requested by fastq_unconsumed in case this R2 VB doesn't have an R1 counterpart and we are allowed to truncate)
    if (final_unconsumed_len == UNCONSUMED_TRUNCATE_VB) {
        txt_file_truncate_final_bytes (vb, &final_unconsumed_len);       
        *R2_vb_truncated_away = true;
    }

    // case: there is not even one full line of text  
    else if (final_unconsumed_len == UNCONSUMED_NEED_MORE_DATA) {
        
        // case: segconf - segconf will try again, increasing the vb_size
        if (segconf_running && !txt_file->no_more_blocks)
            {}

        // case: R2 doesn't have enough data to sync with R1 (confirmed by counting lines) get more data
        else if (IS_R2) {
            uint32_t r1_num_lines = fastq_get_R1_num_lines (vb);
            uint32_t r2_num_lines = str_count_char (B1STtxt, Ltxt, '\n') / 4;

            // error if we did in fact read the same amount of lines as R1, uncompress everything, and still didn't find the QNAME
            ASSINP (r2_num_lines < r1_num_lines,
                    NO_PAIR_FMT_PREFIX "read name \"%s\" is missing in %s. codec=%s R2_vb_i=%u R2_num_lines=%u R1_vb_i=%u R1_num_lines=%u)%s",
                    txt_name, fastq_get_R1_last_qname(vb), txt_name, codec_name (txt_file->effective_codec),
                    vb->vblock_i, (uint32_t)r2_num_lines, fastq_get_R1_vb_i (vb), r1_num_lines, NO_PAIR_FMT_SUFFIX);
            
            // error if there isn't enough data in the file
            ASSINP (!txt_file->no_more_blocks || txt_file->unconsumed_txt.len,
                    NO_PAIR_FMT_PREFIX "read name \"%s\" is missing in %s, because the file appears shorter than its R1. codec=%s vb_i=%u R2_num_lines=%u R1_vb_i=%u R1_num_lines=%u)%s",
                    txt_name, fastq_get_R1_last_qname(vb), txt_name, codec_name (txt_file->effective_codec),
                    vb->vblock_i, (uint32_t)r2_num_lines, fastq_get_R1_vb_i (vb), r1_num_lines, NO_PAIR_FMT_SUFFIX);
        } 
        
        // case: file is truncated
        else if (flag.truncate) 
            txt_file_truncate_final_bytes (vb, &final_unconsumed_len);       

        // case: file is truncated, but user didn't specify --truncate 
        else 
            ABORT ("Reason: failed to find a full line %sin vb=%s data_type=%s txt_data.len=%u txt_file->effective_codec=%s is_last_vb_in_txt_file=%s interleaved=%s.\n"
                   "Known possible causes:\n"
                   "- The file is %s %s. " _TIP "try running with --truncate\n"
                   "- The file is not a %s file.\n"
                   "VB dumped: %s\n",  
                   DTPT(is_binary) ? "" : "(i.e. newline-terminated) ",
                   VB_NAME, dt_name (txt_file->data_type), Ltxt, codec_name (txt_file->effective_codec), TF(vb->is_last_vb_in_txt_file), TF(segconf.is_interleaved),
                   DTPT(is_binary) ? "truncated but not on the boundary of the" : "missing a newline on the last", DTPT(line_name),
                   TXT_DT(REF) ? "FASTA" : dt_name (txt_file->data_type),
                   txtfile_dump_vb (vb, txt_name, NULL).s);
    }

done:
    if (vb->gz_inflate_mem)
        libdeflate_free_decompressor ((struct libdeflate_decompressor **)&vb->gz_inflate_mem, __FUNCLINE); // also sets gz_inflate_mem to NULL
    
    // pass any unconsumed data at the end of txt_data to the next vb
    if (final_unconsumed_len > 0) {
        if (final_unconsumed_len < 0) final_unconsumed_len = Ltxt; // entire VB is unconsumed

        ASSERT (final_unconsumed_len <= Ltxt, "expecting final_unconsumed_len=%d <= Ltxt=%u", final_unconsumed_len, Ltxt);
        
        buf_insert (evb, txt_file->unconsumed_txt, char, 0, Btxt (Ltxt - final_unconsumed_len), final_unconsumed_len, "txt_file->unconsumed_txt");
        Ltxt -= final_unconsumed_len; 
    }

    return final_unconsumed_len >= 0; // false means more data is needed
}

// reads some bytes from beginning of a file into evb->txt_data
void txtfile_query_first_bytes_in_file (rom filename, uint32_t len)
{
    ASSERTNOTINUSE (evb->txt_data);
    ASSERTISNULL (txt_file);
    ASSERTMAINTHREAD;

    // note: used for FASTQ/A files, in which case gz discovery is skipped in file_open_txt_read_gz.
    txt_file = file_open_txt_read (filename);
    ASSINP (txt_file, "failed to open file %s", filename);
    
    uint32_t save_vb_size = segconf.vb_size;
    segconf.vb_size = SAM_MAX_QNAME_LEN+1;
    
    txtfile_read_vblock (evb);
        
    file_close (&txt_file);

    segconf.vb_size = save_vb_size;
}

static bool seggable_size_is_modifiable (void)
{
    return !is_read_via_ext_decompressor (txt_file) && !TXT_IS(NONE);
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
        if (TXT_IS(NONE)) 
            source_comp_ratio = 1; 

        else {    
            double plain_len = txt_file->txt_data_so_far_single + txt_file->unconsumed_txt.len; //  all data that has been decompressed
            double comp_len  = TXT_IS(BZ2)                       ? BZ2_consumed ((BZFILE *)txt_file->file)
                             : txt_file->discover_during_segconf ? segconf.gz_comp_size 
                             :                                     txt_file->disk_so_far - txt_file->gz_data.len; // data read from disk, excluding data still awaiting decompression
            
            // if we are uncompressing with igzip, we might have some unconsumed uncompressed txt still in the igzip buffers
            if (txt_file->igzip_state.len) {
                struct inflate_state *state = B1ST (struct inflate_state, txt_file->igzip_state);
                plain_len += (state->tmp_out_valid - state->tmp_out_processed);
            }
            
            // case: header is whole BGZF blocks - remove header from calculation to get a better estimate of the seggable compression ratio
            if (txt_file->header_size_bgzf) { 
                plain_len -= txt_file->header_size;
                comp_len  -= txt_file->header_size_bgzf;
            }

            source_comp_ratio = plain_len / MAX_(comp_len, 1.0);
        }
    }
        
    // external decompressors
    else switch (txt_file->src_codec) {
        case CODEC_BCF:  source_comp_ratio = 10; break; // note: .bcf files might be compressed or uncompressed - we have no way of knowing as "bcftools view" always serves them to us in plain VCF format. These ratios are assuming the bcf is compressed as it normally is.
        case CODEC_XZ:   source_comp_ratio = 15; break;
        case CODEC_CRAM: source_comp_ratio = 25; break;
        case CODEC_ORA:  source_comp_ratio = 25; break;
        case CODEC_ZIP:  source_comp_ratio = 3;  break;
        default: ABORT ("unspecified txt_file->src_codec=%s (%u)", codec_name (txt_file->src_codec), txt_file->src_codec);
    }
        
    txt_file->est_seggable_size = MAX_(0.0, (double)disk_size * source_comp_ratio - (double)txt_file->header_size);
        
    if (segconf_running) {
        txt_file->txt_data_so_far_single = txt_file->header_size; // roll back as we will re-account for this data in VB=1

        if (IS_GZIP(txt_file->src_codec)) 
            segconf.gz_comp_ratio = source_comp_ratio;
    }
}

uint32_t txt_data_alloc_size (uint32_t vb_size) 
{
    return TXT_IS(MGSP) ? MAX_(24, txt_file->max_mgsp_blocks_in_vb) * (txt_file->max_mgzip_isize + 1/*1 + for last gz block extra length*/) + TXTFILE_READ_VB_PADDING  
         : vb_size     ? vb_size + (TXT_IS_MGZIP ? txt_file->max_mgzip_isize : 0) + TXTFILE_READ_VB_PADDING 
         :               0;
}

// performs a single I/O read operation - returns number of bytes read
// data is placed in vb->txt_data, except if its BGZF and uncompress=false - compressed data is placed in vb->comp_txt_data
static uint32_t txtfile_read_block (VBlockP vb, uint32_t bytes_requested,
                                    bool uncompress,    // MGZIP codecs: whether to uncompress the data. ignored if not MGZIP
                                    bool *is_data_read) // out: true if read any data, including an isize=0 gz block
{   
    START_TIMER;

    if (txt_file->no_more_blocks) return 0; // nothing more to read

    uint32_t uncomp_len = 0;

    if (TXT_IS_MGZIP)     
        uncomp_len = txtfile_read_block_mgzip (vb, bytes_requested, uncompress, is_data_read);  // note: will possibly read more bytes than requested if last mgzip block goes over
    
    else if (TXT_IS(GZ))   
        uncomp_len = txtfile_read_block_igzip (vb, bytes_requested, is_data_read);

    else if (TXT_IS(NONE)) {
        uncomp_len = txtfile_read_block_plain (vb, bytes_requested); 
        *is_data_read = !!uncomp_len; 
    }

    else if (TXT_IS(BZ2)) {
        uncomp_len = txtfile_read_block_bz2 (vb, bytes_requested); 
        *is_data_read = !!uncomp_len; 
    }

    else
        ABORT ("unsupported codec %s", codec_name (txt_file->effective_codec));

    COPY_TIMER_EVB (read);
    return uncomp_len;
}

void txtfile_enforce_vb_large_enough_for_mgzip (uint32_t vb_size, EnforceVbSizeAction action)
{
    // max_block_size exists for fixed-block-size codecs: VB data read will be >= this size
    uint32_t max_block_size = mgzip_get_max_block_size();

    // VB needs be large enough to accomodate a full MGZIP block (since we can only read full blocks) plus a partial line passed 
    // down from the previous VB. Only relevant to FASTQ files with large MGZIP blocks (which might approach small VB sizes).
    uint32_t large_enough = max_block_size + 32 KB/*enough for a partial FASTQ read passed down*/;

    if (vb_size >= large_enough) return;

    switch (action) {
        case RESTART_IF_NOT:
            RESTART ("--no-bgzf", "%s: vblock_size=%s too small for codec=%s, needs to be at least %u bytes",
                     txt_name, str_size(vb_size).s, codec_name(txt_file->effective_codec), large_enough);
        
        case RESIZE_IF_NOT: segconf.vb_size = ROUNDUP1M (large_enough); break;
        
        default: ABORT0 ("unknown action");
    }
}

bool txtfile_uncompress_mgzip_at_read (void)
{
    return segconf_running     // segconf doesn't have a compute thread
        || IS_SHOW_BAI         // uncompressing a TBI file
        || IS_BIOPSY_BYTES
        || flag.biopsy         
        || flag.make_reference // unconsumed callback for make-reference needs to inspect the whole data
        || flag.zip_lines_counted_at_init_vb; // *_zip_init_vb needs to count lines (set at segconf)
}

// ZIP main thread
void txtfile_read_vblock (VBlockP vb)
{
    START_TIMER;

    ASSERTNOTNULL (txt_file);
    ASSERT_DT_FUNC (txt_file, unconsumed);
    ASSERT (vb == evb || IN_RANGX (segconf.vb_size, ABSOLUTE_MIN_VBLOCK_MEMORY, ABSOLUTE_MAX_VBLOCK_MEMORY) || segconf_running,
            "Invalid vb_size=%"PRIu64" comp_i(0-based)=%u", segconf.vb_size, z_file->num_txts_so_far-1);

    if (txt_file->no_more_blocks    && // nothing left on disk
        !txt_file->unconsumed_txt.len) // no data passed down to us from previous vb
        return; // we're done

    bool is_mgzip = TXT_IS_MGZIP;

    bool always_uncompress = !is_mgzip ||  // GZ, BZ2 and NONE always return uncompressed data anyway 
                             txtfile_uncompress_mgzip_at_read();

    vb->comp_i = flag.zip_comp_i;  // needed for VB_NAME

    // Note: VB might grow 1. if 0 (for large variable length MGZIP blocks) and 2. to match a FASTQ R2 vb to its R1 pair
    uint32_t my_vb_size = IS_R2 ? MAX_(fastq_get_R1_txt_data_len (vb), segconf.vb_size) // note: if no correspoding VB we go ahead and try to read data anyway, to make sure there is none
                                : segconf.vb_size; 
    ASSERTNOTZERO (my_vb_size);

    buf_alloc (vb, &vb->txt_data, 0, txt_data_alloc_size (my_vb_size), char, 1.05, "txt_data");    

    if (is_mgzip) {
        mgzip_zip_init_vb (vb); 
        txtfile_enforce_vb_large_enough_for_mgzip (my_vb_size, RESTART_IF_NOT); 
    }
    
    uint32_t max_block_size = mgzip_get_max_block_size(); // 1 if not a fixed-block-size mgzip codec

    while (1) {     
        uint32_t bytes_requested = IS_R2 ? ((double)my_vb_size * 1.03 + (max_block_size - 1)) // add 3% vs R1 (VB might be slightly bigger if reads on average are a bit longer) and round up to the next full block
                                         : MIN_(my_vb_size, 1 GB /* read() can't handle more */);
        bytes_requested -= MIN_(Ltxt, bytes_requested); // reduce data read, if we already have some from unconsumed_txt or previous iterations

        buf_alloc (vb, &vb->txt_data, bytes_requested, 0, char, 1.05, NULL); 

        bool is_data_passed_down=false, is_data_read_from_disk=false; 

        // start with using the data passed down from the previous VB
        if (bytes_requested && txt_file->unconsumed_txt.len) {
            uint64_t bytes_moved = MIN_(txt_file->unconsumed_txt.len, bytes_requested);
            buf_append (vb, vb->txt_data, char, B1STc(txt_file->unconsumed_txt), bytes_moved, "txt_data");
            buf_remove (txt_file->unconsumed_txt, char, 0, bytes_moved); // possibly some data remains
        
            is_data_passed_down = true;
            bytes_requested -= bytes_moved;
        }
        
        if (bytes_requested)
            txtfile_read_block (vb, bytes_requested, always_uncompress, &is_data_read_from_disk);

        bool is_data_read = is_data_passed_down || is_data_read_from_disk;

        // with an MGZIP codec, we might be filled up even without completely filling my_vb_size 
        // if there is room left for only a partial MGZIP block (we can't read partial blocks)
        uint32_t filled_up = my_vb_size - (is_mgzip ? (max_block_size - 1) : 0);

        // case: each VB contains one gz block (or group of blocks)
        if (TXT_IS_VB_SIZE_BY_MGZIP)
            break;

        // check if we're done: if can't read any more data, or if the VB appears full
        // note: is_data_read can be true even if len=0: when reading an isize=0 gz block
        else if (is_data_read && Ltxt < filled_up) 
            continue; // continue reading

        // case no more data read for R2 (for a subsequent read after previous txtfile_get_unconsumed_to_pass_to_next_vb couldn't find a matching QNAME)
        else if (!is_data_read && !Ltxt && IS_R2 && fastq_get_R1_vb_i(vb))
            ABORTINP ("--pair: file %s has less reads than its R1 counterpart. For example, read \"%s\" is missing for vb=%s R1_vb_1=%u",
                      txt_name, fastq_get_R1_last_qname (vb), VB_NAME, fastq_get_R1_vb_i (vb));

        // case: no data read nor pass up from prev vb (and hence also no data to pass down to next vb)            
        else if (!is_data_read && !Ltxt)
            goto done; // cancel this VB - we have no more data
    
        // callback to decide what part of txt_data to pass up to the next VB (usually partial lines, but sometimes more)
        // note: even if we haven't read any new data (everything was passed down), we still might data to pass up - eg
        // in FASTA with make-reference if we have a lots of small contigs, each VB will take one contig and pass up the remaining
        if (!vb->gz_blocks.len        ||  
            !txt_file->no_more_blocks || // usually, there's no need to check at the end of the file as there is no next VB
            flag.truncate             || // user informs us that there might be a partial line to be truncates, so proceed even if EOF
            IS_R2) {                     // In R2 we need to search for matching R1 QNAME

            bool R2_vb_truncated_away=false;

            // case where we don't need to check for unconusmed and/or sync with R1
            if (TXT_IS_IN_SYNC && !IS_R2 && !flag.truncate && txt_file->effective_codec == z_file->comp_eff_codec[flag.zip_comp_i-1])
                break; 
            
            else if (flag.biopsy_bytes.len)
                break; // --biopsy-bytes: we return everything read

            // case: we found the unconsumed text and (if needed) synced our R2 vb to its R1 counterpart. all good.
            else if (txtfile_get_unconsumed_to_pass_to_next_vb (vb, &R2_vb_truncated_away)) {
                // case: this is R2 with --truncate, this VB consisted of data that entirely goes beyond R1 - entire VB was eliminated
                if (R2_vb_truncated_away && !txt_file->no_more_blocks)
                    continue; // continue to read more data and elimitate it
                else
                    break;
            } 

            // case 1: very long single line in segconf 
            // case 2: cannot find matching QNAME in R2
            else { 
                my_vb_size *= 1.25;
                ASSERT (my_vb_size <= ABSOLUTE_MAX_VBLOCK_MEMORY, "%s: VBlock too big, > %s, when trying to grow vb", VB_NAME, str_size(ABSOLUTE_MAX_VBLOCK_MEMORY).s);
                continue;
            }
        }

        // case: we read the entire file, and this VB has some data
        else if (!is_data_read && Ltxt)
            break;
    };

    // case: compute thread should uncompress
    if (is_mgzip && !always_uncompress)
        vb->txt_codec = txt_file->effective_codec;
    
    // copy unconsumed or partially consumed gz_blocks to txt_file->unconsumed_mgzip_blocks.
    // note: this might happen even if final_unconsumed_len=0: in case we inherited gz_blocks from the 
    // previous VB (passed to us in mgzip_zip_init_vb), and we didn't use all of them because we took only
    // part of the txt_file->unconsumed_txt. This can happen, for example, if segconf 
    // generated lots of unconsumed_txt and unconsumed_mgzip_blocks, but segconf.vb_size is smaller for vb>1.
    if (is_mgzip) 
        mgzip_copy_unconsumed_blocks (vb);

    vb->vb_position_txt_file = txt_file->txt_data_so_far_single;
    vb->is_last_vb_in_txt_file = txt_file->no_more_blocks && !txt_file->unconsumed_txt.len;
    txt_file->txt_data_so_far_single += Ltxt;

    zip_init_vb (vb);

    if (!txt_file->est_seggable_size || seggable_size_is_modifiable())
        txtfile_set_seggable_size();

    if (!segconf_running) {
        // case: R1, and codec might or might not have mgzip blocks in sync with R2: store the info R2 will need for deciding this
        if (IS_R1 && Ltxt) 
            buf_append_one (z_file->R1_txt_data_lens, vb->txt_data.len32);

        biopsy_take (vb);
        dispatcher_increment_progress ("read", txt_file->est_num_lines ? (Ltxt / MAX_(segconf.line_len,1)) : Ltxt);
    }

done:
    if (flag_is_show_vblocks (TASK_ZIP)) 
        iprintf ("VB_READ(id=%d) vb=%s Ltxt=%u vb_position_txt_file=%"PRIu64" unconsumed_txt.len=%u is_last_vb_in_txt_file=%s\n", 
                 vb->id, VB_NAME, Ltxt, vb->vb_position_txt_file, txt_file->unconsumed_txt.len32, TF(vb->is_last_vb_in_txt_file));

    if (always_uncompress) 
        buf_free (vb->comp_txt_data); 

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
StrText1K txtfile_codec_name (FileP z_file/*obscures global*/, CompIType comp_i, bool obscure_fname) 
{
    StrText1K s;

    Codec src_codec = (z_file->comp_src_codec[comp_i] == CODEC_BAM) ? z_file->comp_eff_codec[comp_i] : z_file->comp_src_codec[comp_i]; // avoid BAM➤BGZF

    if (!IN_RANGE (comp_i, 0, MAX_NUM_COMPS))
        snprintf (s.s, sizeof (s.s), "comp_i=%u out_of_range", comp_i);

    else if (src_codec == CODEC_BGZF) {
            if (z_file->comp_bgzf[comp_i].level < BGZF_COMP_LEVEL_UNKNOWN)
                snprintf (s.s, sizeof (s.s), "BGZF(%s[%d])", bgzf_library_name (z_file->comp_bgzf[comp_i].library, false), z_file->comp_bgzf[comp_i].level);
            else
                strcpy (s.s, "BGZF(unknown_lib)");
        }
        
    else if (src_codec == CODEC_GZ) 
        snprintf (s.s, sizeof (s), "GZ(%.800s%.12s%.12s)", 
                  display_gz_header (z_file->comp_gz_header[comp_i], MAX_GZ_HEADER_LEN, obscure_fname).s,
                  cond_str (z_file->gz_isize[comp_i][0], "-", str_size (z_file->gz_isize[comp_i][0]).s),
                  cond_str (z_file->gz_isize[comp_i][1], "-", str_size (z_file->gz_isize[comp_i][1]).s));

    else 
        strcpy (s.s, codec_name (src_codec));

    return s;
}