// ------------------------------------------------------------------
//   vcf_header.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <sys/types.h>
#include <sys/stat.h>

#include "genozip.h"
#include "zfile.h"
#include "vcffile.h"
#include "vcf_header.h"
#include "vb.h"
#include "crypt.h"
#include "version.h"
#include "endianness.h"
#include "file.h"

unsigned global_num_samples = 0; // a global param - assigned once before any thread is crated, and never changes - so no thread safety issue
static Buffer global_vcf_header_line = EMPTY_BUFFER; // header line of first VCF file read - use to compare to subsequent files to make sure they have the same header during concat

static bool vcf_header_set_globals(VariantBlock *vb, const char *filename, Buffer *vcf_header)
{
    static const char *vcf_header_line_filename = NULL; // file from which the header line was taken

    // count tabs in last line which should be the field header line
    unsigned tab_count = 0;
    for (int i=vcf_header->len-1; i >= 0; i--) {
        
        if (vcf_header->data[i] == '\t')
            tab_count++;
        
        // if this is the beginning of field header line 
        else if (vcf_header->data[i] == '#' && (i==0 || vcf_header->data[i-1] == '\n' || vcf_header->data[i-1] == '\r')) {
        
            // if first vcf file - copy the header to the global
            if (!buf_is_allocated (&global_vcf_header_line)) {
                buf_copy (vb, &global_vcf_header_line, vcf_header, 1, i, vcf_header->len - i, "global_vcf_header_line", 0);
                vcf_header_line_filename = filename;
            }

            // subsequent files - if we're in concat mode just compare to make sure the header is the same
            else if (flag_concat_mode && 
                     (vcf_header->len-i != global_vcf_header_line.len || memcmp (global_vcf_header_line.data, &vcf_header->data[i], global_vcf_header_line.len))) {

                fprintf (stderr, "%s: skipping %s: it has a different VCF header line than %s, see below:\n"
                                 "========= %s =========\n"
                                 "%.*s"
                                 "========= %s ==========\n"
                                 "%.*s"
                                 "=======================================\n", 
                         global_cmd, filename, vcf_header_line_filename,
                         vcf_header_line_filename, global_vcf_header_line.len, global_vcf_header_line.data,
                         filename, vcf_header->len-i, &vcf_header->data[i]);
                return false;
            }

            //count samples
            if (tab_count >= 9) {
                global_num_samples = tab_count-8; // number of samples;
                return true;
            }

            ASSERT0 (tab_count != 8, "Error: invalid VCF file - field header line contains a FORMAT field but no samples");

            ASSERT (tab_count == 7, "Error: invalid VCF file - field header line contains only %d fields, expecting at least 8", tab_count+1);

            global_num_samples = 0; // a VCF file without samples - that's perfectly fine
            return true; 
        }
    }

    ABORT ("Error: invalid VCF file - it does not contain a field header line; tab_count=%u", tab_count+1);
    return false; // avoid complication warnings
}

// reads VCF header and writes its compressed form to the GENOZIP file. returns num_samples.
bool vcf_header_vcf_to_genozip (VariantBlock *vb, unsigned *line_i, Buffer **first_data_line)
{    
    static Buffer vcf_header_line = EMPTY_BUFFER; // serves to read the header, then its the first line in the data, and again the header when starting the next vcf file
    static Buffer vcf_header_text = EMPTY_BUFFER;

    // in concat mode, we write the header to the genozip file, only for the first vcf file
    //bool use_vcf_header = !flag_concat_mode || !buf_is_allocated (&global_vcf_header_line) /* first vcf */;
    
    // line might be used as the first line, so it cannot be free at the end of this function
    // however, by the time we come here again, at the next VCF file, it is no longer needed
    if (buf_is_allocated (&vcf_header_line)) buf_free (&vcf_header_line); 

    *first_data_line = NULL;

    const unsigned INITIAL_BUF_SIZE = 65536;

    buf_alloc (vb, &vcf_header_text, INITIAL_BUF_SIZE, 0, "vcf_header_text", 0);

    bool is_first_vcf = !buf_is_allocated (&global_vcf_header_line); 

    bool skip_md5_vcf_header = flag_concat_mode && !is_first_vcf /* not first vcf */;
    
    while (1) {
        bool success = vcffile_get_line (vb, *line_i + 1, skip_md5_vcf_header, &vcf_header_line, "vcf_header_line");
        if (!success) break; // end of header - no data lines in this VCF file

        (*line_i)++;

        // case : we reached the end of the header - we exit here

        if (vcf_header_line.data[0] != '#') { // end of header - we have read the first data line
            vcf_header_line.data[vcf_header_line.len] = '\0'; // remove newline - data line reader is expecting lines without the newline

            *first_data_line = &vcf_header_line;

            break;
        }

        buf_alloc (vb, &vcf_header_text, vcf_header_line.len + vcf_header_text.len + 1, 2, "vcf_header_text", 1); // +1 for terminating \0

        memcpy (&vcf_header_text.data[vcf_header_text.len], vcf_header_line.data, vcf_header_line.len);
        
        vcf_header_text.len += vcf_header_line.len;
        vcf_header_text.data[vcf_header_text.len] = '\0';
    }

    // case - vcf header was found 
    if (vcf_header_text.len) {

        bool can_concatenate = vcf_header_set_globals(vb, vb->vcf_file->name, &vcf_header_text);
        if (!can_concatenate) { 
            // this is the second+ file in a concatenation list, but its samples are incompatible
            buf_free (&vcf_header_text);
            return false;
        }

        if (vb->z_file) {
            //if (use_vcf_header)
            zfile_write_vcf_header (vb, &vcf_header_text, is_first_vcf); // we write all headers in concat mode too, to support --split
            //else
            //    vb->z_file->vcf_data_so_far  += vcf_header_text.len; // length of the original VCF header
        }

        vb->vcf_file->section_bytes[SEC_VCF_HEADER] = vcf_header_text.len;
        vb->z_file  ->section_bytes[SEC_VCF_HEADER] = vb->z_section_bytes[SEC_VCF_HEADER]; // comes from zfile_compress
    }

    // case : header not found, but data line found
    else 
        ASSERT0 (! *first_data_line, "Error: file has no VCF header");

    // case : empty file - not an error (caller will see that *first_data_line is NULL)
    buf_free (&vcf_header_text);

    return true; // everything's good
}

// returns true if there's a file or false if its an empty file
bool vcf_header_genozip_to_vcf (VariantBlock *vb, Md5Hash *digest)
{
    // read vcf header if its not already read. it maybe already read, if we're in flag_split mode, for second component onwards
    // (it would have been read by previous component and stored for us)
    if (!buf_is_allocated (&vb->z_file->next_vcf_header)) {
        int ret = zfile_read_one_section (vb, &vb->z_file->next_vcf_header, "z_file->next_vcf_header", 
                                        sizeof(SectionHeaderVCFHeader), SEC_VCF_HEADER, true);
        if (ret == EOF) {
            buf_free (&vb->z_file->next_vcf_header);
            return false; // empty file (or in case of split mode - no more components) - not an error
        }
    }

    // handle the GENOZIP header of the VCF header section
    SectionHeaderVCFHeader *header = (SectionHeaderVCFHeader *)vb->z_file->next_vcf_header.data;

    ASSERT (header->genozip_version == GENOZIP_FILE_FORMAT_VERSION, "Error: file version %u is newer than the latest version supported %u. Please upgrade.",
            header->genozip_version, GENOZIP_FILE_FORMAT_VERSION);

    ASSERT (BGEN32 (header->h.compressed_offset) == crypt_padded_len (sizeof(SectionHeaderVCFHeader)), "Error: invalid VCF header's header size: header->h.compressed_offset=%u, expecting=%u", BGEN32 (header->h.compressed_offset), (unsigned)sizeof(SectionHeaderVCFHeader));

    // in split mode - we open the output VCF file of the component
    if (flag_split) {
        ASSERT0 (!vb->vcf_file, "Error: not expecting vb->vcf_file to be open already in split mode");
        vb->vcf_file = file_open (header->vcf_filename, WRITE, VCF);
    }

    bool first_vcf = !buf_is_allocated (&global_vcf_header_line);

    if (first_vcf || !flag_concat_mode) {
        vb->z_file->num_lines_concat     = vb->vcf_file->num_lines_concat     = BGEN64 (header->num_lines);
        vb->z_file->vcf_data_size_concat = vb->vcf_file->vcf_data_size_concat = BGEN64 (header->vcf_data_size);
    }

    *digest = flag_split ? header->md5_hash_single : header->md5_hash_concat;
    vb->vcf_file->has_md5 = digest->ulls[0] || digest->ulls[1]; // has_md5 iff not all 0. note: a chance of 1 in about 10^38 that we will get all-0 by chance in which case will won't perform the md5 comparison
        
    // now get the text of the VCF header itself
    static Buffer vcf_header_buf = EMPTY_BUFFER;
    zfile_uncompress_section (vb, header, &vcf_header_buf, SEC_VCF_HEADER);

    bool can_concatenate = vcf_header_set_globals (vb, vb->z_file->name, &vcf_header_buf);
    if (!can_concatenate) {
        buf_free (&vb->z_file->next_vcf_header);
        buf_free (&vcf_header_buf);
        return false;
    }

    // write vcf header if not in concat mode, or, in concat mode, we write the vcf header, only for the first genozip file
    if (first_vcf || !flag_concat_mode)
        vcffile_write_to_disk (vb->vcf_file, &vcf_header_buf);
    
    // if we didn't write the header (bc 2nd+ file in concat mode) - just account for it in MD5 if needed (this is normally done by vcffile_write_to_disk())
    else if (vb->vcf_file->has_md5)
        md5_update (&vb->vcf_file->md5_ctx_concat, vcf_header_buf.data, vcf_header_buf.len, true);

    buf_free (&vb->z_file->next_vcf_header);
    buf_free (&vcf_header_buf);

    return true;
}

// returns the the VCF header section's header from a GENOZIP file - used by main_list
bool vcf_header_get_vcf_header (File *z_file, SectionHeaderVCFHeader *vcf_header_header, bool *encrypted)
{
    int bytes = fread ((char*)vcf_header_header, 1, crypt_padded_len (sizeof(SectionHeaderVCFHeader)), (FILE *)z_file->file);
    
    if (bytes < sizeof(SectionHeaderVCFHeader)) {
        *encrypted = false; // we can't read for some reason - but its not an encryption issue
        return false;
    }

    if (BGEN32 (vcf_header_header->h.magic) == GENOZIP_MAGIC)
        return true; // not encrypted

    if (crypt_have_password()) {
        VariantBlock fake_vb;
        memset (&fake_vb, 0, sizeof(fake_vb));
        crypt_do (&fake_vb, (uint8_t *)vcf_header_header, crypt_padded_len (sizeof (SectionHeaderVCFHeader)), 0, -1);

        if (BGEN32 (vcf_header_header->h.magic) == GENOZIP_MAGIC)
            return true; // we successfully decrypted it with the user provided password
    }

    *encrypted = true;
    return false; // we either don't have a password, or the password we have didn't work
}
