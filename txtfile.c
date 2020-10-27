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
#include "zlib/zlib.h"

static bool is_first_txt = true; 
static uint32_t total_bound_txt_headers_len = 0;

uint32_t txtfile_get_bound_headers_len(void) { return total_bound_txt_headers_len; }

void txtfile_update_md5 (const char *data, uint32_t len, bool is_2ndplus_txt_header)
{
    if (flag_md5) {
        if (flag_bind && !is_2ndplus_txt_header)
            md5_update (&z_file->md5_ctx_bound, data, len);
        
        md5_update (&z_file->md5_ctx_single, data, len);
    }
}

// peformms a single I/O read operation - returns number of bytes read 
uint32_t txtfile_read_block (char *data, uint32_t max_bytes)
{
    START_TIMER;

    int32_t bytes_read=0;

    if (file_is_plain_or_ext_decompressor (txt_file)) {
        
        bytes_read = read (fileno((FILE *)txt_file->file), data, max_bytes); // -1 if error in libc
        ASSERT (bytes_read >= 0, "Error: read failed from %s: %s", txt_name, strerror(errno));

        // bytes_read=0 and we're using an external decompressor - it is either EOF or
        // there is an error. In any event, the decompressor is done and we can suck in its stderr to inspect it
        if (!bytes_read && file_is_read_via_ext_decompressor (txt_file)) {
            file_assert_ext_decompressor();
            goto finish; // all is good - just a normal end-of-file
        }
        
        txt_file->disk_so_far += (int64_t)bytes_read;

#ifdef _WIN32
        // in Windows using Powershell, the first 3 characters on an stdin pipe are BOM: 0xEF,0xBB,0xBF https://en.wikipedia.org/wiki/Byte_order_mark
        // these charactes are not in 7-bit ASCII, so highly unlikely to be present natrually in a VCF file
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
    }
    else if (txt_file->codec == CODEC_GZ || txt_file->codec == CODEC_BGZ) {
        bytes_read = gzfread (data, 1, max_bytes, (gzFile)txt_file->file);
        
        if (bytes_read)
            txt_file->disk_so_far = gzconsumed64 ((gzFile)txt_file->file); 
    }
    else if (txt_file->codec == CODEC_BZ2) { 
        bytes_read = BZ2_bzread ((BZFILE *)txt_file->file, data, max_bytes);

        if (bytes_read)
            txt_file->disk_so_far = BZ2_consumed ((BZFILE *)txt_file->file); 
    } 
    else {
        ABORT ("txtfile_read_block: Invalid file type %s (codec=%s)", ft_name (txt_file->type), codec_name (txt_file->codec));
    }
    
finish:
    COPY_TIMER_VB (evb, read);

    return bytes_read;
}

// ZIP: returns the hash of the header
Md5Hash txtfile_read_header (bool is_first_txt, bool header_required,
                             char first_char)  // first character in every header line
{
    START_TIMER;

    int32_t bytes_read;
    char prev_char='\n';

    // read data from the file until either 1. EOF is reached 2. end of txt header is reached
    for (bool first_block=true;; first_block=false) { 

        // enlarge if needed        
        if (!evb->txt_data.data || evb->txt_data.size - evb->txt_data.len < READ_BUFFER_SIZE) 
            buf_alloc (evb, &evb->txt_data, evb->txt_data.size + READ_BUFFER_SIZE + 1 /* for \0 */, 1.2, "txt_data", 0);    

        // case: first batch of data was already read in txtfile_read_header_and_100_lines(), use it
        if (first_block && evb->txt_data.len) {
            bytes_read = (int32_t)evb->txt_data.len;
            evb->txt_data.len = 0;
        }
        else
            bytes_read = txtfile_read_block (AFTERENT (char, evb->txt_data), READ_BUFFER_SIZE);

        if (!bytes_read) { // EOF
            ASSERT (!evb->txt_data.len || evb->txt_data.data[evb->txt_data.len-1] == '\n', 
                    "Error: invalid %s header in %s - expecting it to end with a newline", dt_name (txt_file->data_type), txt_name);
            goto finish;
        }

        const char *this_read = &evb->txt_data.data[evb->txt_data.len];

        ASSERT (!header_required || evb->txt_data.len || this_read[0] == first_char,
                "Error: %s is missing a %s header - expecting first character in file to be %c", 
                txt_name, dt_name (txt_file->data_type), first_char);

        // check stop condition - a line not beginning with a 'first_char'
        for (int i=0; i < bytes_read; i++) { // start from 1 back just in case it is a newline, and end 1 char before bc our test is 2 chars
            if (this_read[i] == '\n') 
                evb->lines.len++;   

            if (prev_char == '\n' && this_read[i] != first_char) {  

                uint32_t txt_header_len = evb->txt_data.len + i;
                evb->txt_data.len += bytes_read; // increase all the way - just for buf_copy

                // the excess data is for the next vb to read 
                buf_copy (evb, &txt_file->unconsumed_txt, &evb->txt_data, 1, txt_header_len,
                          bytes_read - i, "txt_file->unconsumed_txt", 0);

                txt_file->txt_data_so_far_single += i; 
                evb->txt_data.len = txt_header_len;

                goto finish;
            }
            prev_char = this_read[i];
        }

        evb->txt_data.len += bytes_read;
        txt_file->txt_data_so_far_single += bytes_read;
    }

finish:        
    // md5 header - with logic related to is_first
    txtfile_update_md5 (evb->txt_data.data, evb->txt_data.len, !is_first_txt);
    Md5Hash header_md5 = md5_snapshot (&z_file->md5_ctx_single);

    // md5 unconsumed_txt - always
    txtfile_update_md5 (txt_file->unconsumed_txt.data, txt_file->unconsumed_txt.len, false);

    *AFTERENT (char, evb->txt_data) = 0; // null-terminate

    COPY_TIMER_VB (evb, txtfile_read_header);

    return header_md5;
}

// returns true if txt_data[txt_i] is the end of a FASTQ line (= block of 4 lines in the file)
static inline bool txtfile_fastq_is_end_of_line (VBlock *vb, 
                                                 int32_t txt_i) // index of a \n in txt_data
{
#   define IS_NL_BEFORE_QUAL_LINE(i) \
        ((i > 3) && ((txt[i-2] == '\n' && txt[i-1] == '+') || /* \n line ending case */ \
                     (txt[i-3] == '\n' && txt[i-2] == '+' && txt[i-1] == '\r'))) /* \r\n line ending case */
    
    ARRAY (char, txt, vb->txt_data);

    // if we're not at the end of the data - we can just look at the next character
    // if it is a @ then that @ is a new record EXCEPT if the previous character is + and then
    // @ is actually a quality value... (we check two previous characters, as the previous one might be \r)
    if (txt_i < vb->txt_data.len-1)
        return txt[txt_i+1] == '@' && !IS_NL_BEFORE_QUAL_LINE(txt_i);

    // if we're at the end of the line, we scan back to the previous \n and check if it is at the +
    for (int32_t i=txt_i-1; i >= 0; i--) 
        if (txt[i] == '\n') 
            return IS_NL_BEFORE_QUAL_LINE(i);

    // we can't find a complete FASTQ block in the entire vb data
    ABORT ("Error when reading %s: last FASTQ record appears truncated, or the record is bigger than vblock", txt_name);
    return 0; // squash compiler warning ; never reaches here
}

// returns true if txt_data[txt_i] is the end of a FASTA read (= next char is > or end-of-file)
static inline bool txtfile_fasta_is_end_of_line (VBlock *vb, 
                                                 int32_t txt_i) // index of a \n in txt_data
{
    ARRAY (char, txt, vb->txt_data);

    // if we're not at the end of the data - we can just look at the next character
    // if it is a @ then that @ is a new record EXCEPT if the previous character is + and then
    // @ is actually a quality value... (we check two previous characters, as the previous one might be \r)
    if (txt_i < vb->txt_data.len-1)
        return txt[txt_i+1] == '>';

    // if we're at the end of the line, we scan back to the previous \n and check if it NOT the >
    for (int32_t i=txt_i-1; i >= 0; i--) 
        if (txt[i] == '\n') 
            return txt[txt_i+1] != '>'; // this row is a sequence row, not a description row

    // we can't find a complete FASTQ block in the entire vb data
    ABORT ("Error when reading %s: last FASTQ record appears truncated, or the record is bigger than vblock", txt_name);
    return 0; // squash compiler warning ; never reaches here
}

// ZIP I/O threads
void txtfile_read_vblock (VBlock *vb)
{
    START_TIMER;

    static Buffer block_md5_buf = EMPTY_BUFFER, block_start_buf = EMPTY_BUFFER, block_len_buf = EMPTY_BUFFER;

    uint64_t pos_before = 0;
    if (file_is_read_via_int_decompressor (txt_file))
        pos_before = file_tell (txt_file);

    buf_alloc (vb, &vb->txt_data, global_max_memory_per_vb, 1, "txt_data", vb->vblock_i);    

    // start with using the unconsumed data from the previous VB (note: copy & free and not move! so we can reuse txt_data next vb)
    if (buf_is_allocated (&txt_file->unconsumed_txt)) {
        buf_copy (vb, &vb->txt_data, &txt_file->unconsumed_txt, 0 ,0 ,0, "txt_data", vb->vblock_i);
        buf_free (&txt_file->unconsumed_txt);
    }

    // read data from the file until either 1. EOF is reached 2. end of block is reached
    uint64_t max_memory_per_vb = global_max_memory_per_vb;
    uint32_t unconsumed_len=0;
    int32_t block_i=0 ; for (; vb->txt_data.len < max_memory_per_vb; block_i++) {  

        buf_alloc (evb, &block_md5_buf,   sizeof (Md5Context) * MAX (100, block_md5_buf.len+1), 2, "block_md5_buf",   0); // static buffer added to evb list - we are in I/O thread now
        buf_alloc (evb, &block_start_buf, sizeof (char *)     * MAX (100, block_md5_buf.len+1), 2, "block_start_buf", 0); 
        buf_alloc (evb, &block_len_buf,   sizeof (uint32_t)   * MAX (100, block_md5_buf.len+1), 2, "block_len_buf",   0); 

        char *start = AFTERENT (char, vb->txt_data);
        uint32_t len = txtfile_read_block (start, MIN (READ_BUFFER_SIZE, max_memory_per_vb - vb->txt_data.len));

        if (!len) { // EOF - we're expecting to have consumed all lines when reaching EOF (this will happen if the last line ends with newline as expected)
            ASSERT (!vb->txt_data.len || vb->txt_data.data[vb->txt_data.len-1] == '\n', "Error: invalid input file %s - expecting it to end with a newline", txt_name);
            break;
        }
            
        // note: we md_udpate after every block, rather on the complete data (vb or txt header) when its done
        // because this way the OS read buffers / disk cache get pre-filled in parallel to our md5
        // Note: we md5 everything we read - even unconsumed data
        txtfile_update_md5 (start, len, false);

        NEXTENT (char *, block_start_buf)   = start;
        NEXTENT (uint32_t, block_len_buf)   = len;
        NEXTENT (Md5Context, block_md5_buf) = flag_bind ? z_file->md5_ctx_bound : z_file->md5_ctx_single; // MD5 of entire file up to and including this block

        vb->txt_data.len += len;

        // case: this is the 2nd file of a fastq pair - make sure it has at least as many fastq "lines" as the first file
        if (flag_pair == PAIR_READ_2 &&  // we are reading the second file of a fastq file pair (with --pair)
            vb->txt_data.len >= max_memory_per_vb && // we are about to exit the loop
            !fastq_txtfile_have_enough_lines (vb, &unconsumed_len)) { // we don't yet have all the data we need

            // if we need more lines - increase memory and keep on reading
            max_memory_per_vb *= 1.1; 
            buf_alloc (vb, &vb->txt_data, max_memory_per_vb, 1, "txt_data", vb->vblock_i);    
        }
    }
    
    // case: FASTA reference file - we allow only one contig (or part of it) per VB - move second contig onwards to next vb
    if (flag_make_reference) {
        bool data_found = false;
        ARRAY (char, txt, vb->txt_data);
        for (unsigned i=0; i < vb->txt_data.len; i++) {
            // just don't allow now-obsolete ';' rather than trying to disentangle comments from descriptions
            ASSERT (txt[i] != ';', "Error: %s contains a ';' character - this is not supported for reference files. Contig descriptions must begin with a >", txt_name);
        
            // if we've encountered a new DESC line after already seeing sequence data, move this DESC line and
            // everything following to the next VB
            if (data_found && txt[i]=='>' && txt[i-1]=='\n') {
                unconsumed_len = vb->txt_data.len - i;
                break;
            }                

            if (!data_found && (txt[i] != '\n' && txt[i] != '\r')) data_found = true; // anything, except for empty lines, is considered data
        }
    }

    // we move the final partial line to the next vb (unless we are already moving more, due to a reference file)
    if (!unconsumed_len) {
        for (int32_t i=vb->txt_data.len-1; i >= 0; i--) {

            if (vb->txt_data.data[i] == '\n') {

                // in FASTQ - an "end of line" is one that the next character is @, or it is the end of the file
                if (txt_file->data_type == DT_FASTQ && !txtfile_fastq_is_end_of_line (vb, i)) continue;

                // when compressing FASTA with a reference - an "end of line" is one that the next character is >, or it is the end of the file
                // note: when compressing FASTA with a reference (eg long reads stored in a FASTA instead of a FASTQ), line cannot be too long - they 
                // must fit in a VB
                if (txt_file->data_type == DT_FASTA && (flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE) && !txtfile_fasta_is_end_of_line (vb, i)) continue;

                unconsumed_len = vb->txt_data.len-1 - i;
                break;
            }
        }
    }

    // if we have some unconsumed data, pass it to the next vb
    if (unconsumed_len) {
        buf_copy (evb, &txt_file->unconsumed_txt, &vb->txt_data, 1, // evb, because dst buffer belongs to File
                  vb->txt_data.len - unconsumed_len, unconsumed_len, "txt_file->unconsumed_txt", vb->vblock_i);

        vb->txt_data.len -= unconsumed_len;

        // complete the MD5 of all data up to and including this VB based on the last full block + the consumed part of the 
        // last part-consumed-part-not-consumed block
        if (flag_md5 && !flag_make_reference) {
            // find the block that is part-consumed-part-not-consumed block (note: blocks can have any size, not necessarily READ_BUFFER_SIZE)
            int32_t partial_block = block_i - 1;
            uint32_t unincluded_len = unconsumed_len;
            ARRAY (uint32_t, block_len, block_len_buf);

            while (unincluded_len > block_len[partial_block] && partial_block >= 0) unincluded_len -= block_len[partial_block--];
            
            // in the extremely unlikely case where unconsumed block goes back and includes part of the unconsumed data passed on
            // from the previous VB, we will just not have a vb->md5_hash_so_far for this VB. no harm.
            if (partial_block >= 0) {
                // use the Md5Context of the last fully-consumed block as the basis, and update it with the consumed part of the partial block
                Md5Context md5_ctx = partial_block ? *ENT (Md5Context, block_md5_buf, partial_block-1) : MD5CONTEXT_NONE; // copy context
                md5_update (&md5_ctx, *ENT (char *, block_start_buf, partial_block), block_len[partial_block] - unincluded_len);

                vb->md5_hash_so_far = md5_snapshot (&md5_ctx);  
            }
        }
    }
    else if (flag_md5 && !flag_make_reference)
        // MD5 of all data up to and including this VB is just the total MD5 of the file so far (as there is no unconsumed data)
        vb->md5_hash_so_far = md5_snapshot (flag_bind ? &z_file->md5_ctx_bound : &z_file->md5_ctx_single);

    vb->vb_position_txt_file = txt_file->txt_data_so_far_single;

    txt_file->txt_data_so_far_single += vb->txt_data.len;
    vb->vb_data_size = vb->txt_data.len; // initial value. it may change if --optimize is used.
    
    if (file_is_read_via_int_decompressor (txt_file))
        vb->vb_data_read_size = file_tell (txt_file) - pos_before; // gz/bz2 compressed bytes read

    if (DTPT(zip_read_one_vb)) DTPT(zip_read_one_vb)(vb);
        
    buf_free (&block_md5_buf);
    buf_free (&block_start_buf);
    buf_free (&block_len_buf);

    COPY_TIMER (txtfile_read_vblock);
}

// read num_lines of the txtfile (after the header), and call test_func for each line. true iff the proportion of lines
// that past the test is at least success_threashold
bool txtfile_test_data (char first_char,    // first character in every header line
                        unsigned num_lines_to_test, // number of lines to test
                        double success_threashold, // proportion of lines that need to pass the test, for this function to return true
                        TxtFileTestFunc test_func)
{
    uint32_t line_start_i = 0;
    unsigned num_lines_so_far = 0; // number of data (non-header) lines
    unsigned successes = 0;

    while (1) {      // read data from the file until either 1. EOF is reached 2. we pass the header + 100 lines
        // enlarge if needed        
        if (!evb->txt_data.data || evb->txt_data.size - evb->txt_data.len < READ_BUFFER_SIZE) 
            buf_alloc (evb, &evb->txt_data, evb->txt_data.size + READ_BUFFER_SIZE + 1 /* for \0 */, 1.2, "txt_data", 0);    

        uint64_t start_read = evb->txt_data.len;
        evb->txt_data.len += txtfile_read_block (AFTERENT (char, evb->txt_data), READ_BUFFER_SIZE);
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

done:
    return (double)successes / (double)num_lines_so_far >= success_threashold;
}

// PIZ
static void txtfile_write_to_disk (Buffer *buf)
{
    if (flag_md5) md5_update (&txt_file->md5_ctx_bound, buf->data, buf->len);

    // if this is SAM/BAM make sure the output header is SAM or BAM according to flag_reconstruct_binary, and regardless of the source file
    if (z_file->data_type == DT_SAM) bam_prepare_txt_header (buf); 

    if (!flag_test) file_write (txt_file, buf->data, buf->len);

    txt_file->txt_data_so_far_single += buf->len;
    txt_file->disk_so_far            += buf->len;
}

void txtfile_write_one_vblock (VBlockP vb)
{
    START_TIMER;

    txtfile_write_to_disk (&vb->txt_data);

    // warn if VB is bad, but don't exit, so that we can continue and reconstruct the rest of the file for better debugging

    char s1[20], s2[20];
    ASSERTW (vb->txt_data.len == vb->vb_data_size || exe_type == EXE_GENOCAT, 
            "Warning: vblock_i=%u (num_lines=%u vb_start_line_in_file=%u) had %s bytes in the original %s file but %s bytes in the reconstructed file (diff=%d)", 
            vb->vblock_i, (uint32_t)vb->lines.len, vb->first_line,
            str_uint_commas (vb->vb_data_size, s1), dt_name (txt_file->data_type), str_uint_commas (vb->txt_data.len, s2), 
            (int32_t)vb->txt_data.len - (int32_t)vb->vb_data_size);

    // if testing, compare MD5 file up to this VB to that calculated on the original file and transferred through SectionHeaderVbHeader
    // note: we cannot test this unbind mode, because the MD5s are commulative since the beginning of the bound file
    if (flag_md5 && !flag_unbind && !md5_is_zero (vb->md5_hash_so_far)) {
        Md5Hash piz_hash_so_far = md5_snapshot (&txt_file->md5_ctx_bound);
        ASSERTW (md5_is_equal (vb->md5_hash_so_far, piz_hash_so_far), 
                "MD5 of reconstructed vblock=%u (%s) differs from original file (%s). To see bad vblock:\n"
                "%s %s | head -n%u | tail -n%u > bug%s", vb->vblock_i, md5_display (piz_hash_so_far), md5_display (vb->md5_hash_so_far), 
                file_viewer (txt_file), file_guess_original_filename (txt_file), 
                DTP (line_height) * (vb->first_line + (uint32_t)vb->lines.len - 1), 
                DTP (line_height) * (uint32_t)vb->lines.len,
                file_plain_text_ext_of_dt (vb->data_type));
    }

    COPY_TIMER (write);
}

// ZIP only - estimate the size of the txt data in this file. affects the hash table size and the progress indicator.
void txtfile_estimate_txt_data_size (VBlock *vb)
{
    uint64_t disk_size = txt_file->disk_size; 

    // case: we don't know the disk file size (because its stdin or a URL where the server doesn't provide the size)
    if (!disk_size) { 
        if (flag_stdin_size) disk_size = flag_stdin_size; // use the user-provided size, if there is one
        else return; // we're unable to estimate if the disk size is not known
    } 
    
    double ratio=1;

    bool is_no_ht_vcf = (txt_file->data_type == DT_VCF && vcf_vb_has_haplotype_data(vb));

    switch (txt_file->codec) {
        // if we decomprssed gz/bz2 data directly - we extrapolate from the observed compression ratio
        case CODEC_GZ:
        case CODEC_BGZ:
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

    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;

    if (DTPT(txt_header_required) == HDR_MUST || DTPT(txt_header_required) == HDR_OK)
        header_md5 = (DTPT(read_txt_header) ? DTPT(read_txt_header) : txtfile_read_header) (is_first_txt, DTPT(txt_header_required) == HDR_MUST, DTPT(txt_header_1st_char)); // reads into evb->txt_data and evb->lines.len
    
    *txt_line_i += (uint32_t)evb->lines.len;

    // for vcf, we need to check if the samples are the same before approving binding (other data types can bind without restriction)
    // for sam, we check that the contigs specified in the header are consistent with the reference given in --reference/--REFERENCE
    bool txt_file_ok = (DTPT (zip_inspect_txt_header) ? DTPT (zip_inspect_txt_header) (&evb->txt_data) : true);
    if (!txt_file_ok) { 
        // this is the second+ file in a bind list, but its samples are incompatible
        buf_free (&evb->txt_data);
        return false;
    }

    // we always write the txt_header section, even if we don't actually have a header, because the section
    // header contains the data about the file
    if (z_file && !flag_test_seg) zfile_write_txt_header (&evb->txt_data, header_md5, is_first_txt); // we write all headers in bound mode too, to support --unbind

    // for stats: combined length of txt headers in this bound file, or only one file if not bound
    if (!flag_bind) total_bound_txt_headers_len=0;
    total_bound_txt_headers_len += evb->txt_data.len; 

    z_file->num_txt_components_so_far++; // when compressing

    buf_free (&evb->txt_data);
    
    is_first_txt = false;

    return true; // everything's good
}

// PIZ: returns offset of header within data, EOF if end of file
bool txtfile_genozip_to_txt_header (const SectionListEntry *sl, Md5Hash *digest) // NULL if we're just skipped this header (2nd+ header in bound file)
{
    z_file->disk_at_beginning_of_this_txt_file = z_file->disk_so_far;
    static Buffer header_section = EMPTY_BUFFER;

    int header_offset = zfile_read_section (z_file, evb, 0, &header_section, "header_section", 
                                            sizeof(SectionHeaderTxtHeader), SEC_TXT_HEADER, sl);
    if (header_offset == EOF) {
        buf_free (&header_section);
        return false; // empty file (or in case of unbind mode - no more components) - not an error
    }

    // handle the GENOZIP header of the txt header section
    SectionHeaderTxtHeader *header = (SectionHeaderTxtHeader *)header_section.data;

    ASSERT (!digest || BGEN32 (header->h.compressed_offset) == crypt_padded_len (sizeof(SectionHeaderTxtHeader)), 
            "Error: invalid txt header's header size: header->h.compressed_offset=%u, expecting=%u", BGEN32 (header->h.compressed_offset), (unsigned)sizeof(SectionHeaderTxtHeader));

    // 1. in unbind mode - we open the output txt file of the component
    // 2. when reading a reference file - we create txt_file here (but don't actually open the physical file)
    if (flag_unbind || flag_reading_reference) {
        ASSERT0 (!txt_file, "Error: not expecting txt_file to be open already in unbind mode or when reading reference");
        
        char *filename = MALLOC (strlen (header->txt_filename) + strlen (flag_unbind) + 1);
        sprintf (filename, "%s%s", flag_unbind, header->txt_filename);

        txt_file = file_open (filename, WRITE, TXT_FILE, z_file->data_type);
        FREE (filename); // file_open copies the names
    }

    txt_file->txt_data_size_single = BGEN64 (header->txt_data_size); 
    txt_file->max_lines_per_vb     = BGEN32 (header->max_lines_per_vb);
    txt_file->codec                = header->compression_type;
    
    if (is_first_txt || flag_unbind) 
        z_file->num_lines = BGEN64 (header->num_lines);

    if (flag_unbind) *digest = header->md5_hash_single; // override md5 from genozip header

    // now get the text of the txt header itself
    static Buffer header_buf = EMPTY_BUFFER;
    
    if (!(flag_show_headers && exe_type == EXE_GENOCAT))
        zfile_uncompress_section (evb, header, &header_buf, "header_buf", 0, SEC_TXT_HEADER);

    bool is_vcf = (z_file->data_type == DT_VCF);

    bool can_bind = (is_vcf && !(flag_show_headers && exe_type == EXE_GENOCAT)) ? vcf_header_set_globals(z_file->name, &header_buf) : true;
    if (!can_bind) {
        buf_free (&header_section);
        buf_free (&header_buf);
        return false;
    }

    if (is_vcf && flag_drop_genotypes) vcf_header_trim_header_line (&header_buf); // drop FORMAT and sample names

    if (is_vcf && flag_header_one) vcf_header_keep_only_last_line (&header_buf);  // drop lines except last (with field and samples name)

    // write txt header if not in bound mode, or, in bound mode, we write the txt header, only for the first genozip file
    if ((is_first_txt || flag_unbind) && !flag_no_header && !flag_reading_reference && !flag_show_headers) {
        txtfile_write_to_disk (&header_buf);

        if (flag_md5 && !md5_is_zero (header->md5_header)) {
            Md5Hash reconstructed_header_len = md5_do (header_buf.data, header_buf.len);

            ASSERTW (md5_is_equal (reconstructed_header_len, header->md5_header), 
                     "MD5 of reconstructed %s header (%s) differs from original file (%s)",
                     dt_name (z_file->data_type), md5_display (reconstructed_header_len), md5_display (header->md5_header));
        }
    }
    
    buf_free (&header_section);
    buf_free (&header_buf);

    z_file->num_txt_components_so_far++;
    is_first_txt = false;

    return true;
}

