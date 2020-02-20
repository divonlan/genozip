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

unsigned global_num_samples      = 0; // a global param - assigned once before any thread is created, and never changes - so no thread safety issue
unsigned global_max_lines_per_vb = 0; // for ZIP this is determined by as a function of the number of samples, and for PIZ by the SectionHeaderVCFHeader.max_lines_per_vb

Buffer global_vcf_header_line = EMPTY_BUFFER; // header line of first VCF file read - use to compare to subsequent files to make sure they have the same header during concat

// global - this names go into the dictionary names on disk. to preserve backward compatability, they should not be changed.
const char *vcf_field_names[] = { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT" };

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
            global_num_samples = (tab_count >= 9) ? tab_count-8 : 0; // note: a VCF file without samples would have tab_count==7 (8 fields) and is perfectly legal

            // set global_max_lines_per_vb when compressing. we decompressing, it will be determined by the input file.
            if (!global_max_lines_per_vb) {
                // if we have fewer samples, we can afford (memory and speed wise) larger variant blocks
                if      (global_num_samples >= 1024) global_max_lines_per_vb = 4   * 1024;
                else if (global_num_samples >= 512)  global_max_lines_per_vb = 8   * 1024;
                else if (global_num_samples >= 256)  global_max_lines_per_vb = 16  * 1024;
                else if (global_num_samples >= 64)   global_max_lines_per_vb = 32  * 1024;
                else if (global_num_samples >= 32)   global_max_lines_per_vb = 64  * 1024;
                else                                 global_max_lines_per_vb = 128 * 1024;
            }

            ASSERT0 (tab_count != 8, "Error: invalid VCF file - field header line contains a FORMAT field but no samples");

            ASSERT (tab_count >= 7, "Error: invalid VCF file - field header line contains only %d fields, expecting at least 8", tab_count+1);

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
        vb->z_file  ->num_sections [SEC_VCF_HEADER]++;
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
    vb->z_file->disk_at_beginning_of_this_vcf_file = vb->z_file->disk_so_far;

    static Buffer vcf_header_section = EMPTY_BUFFER;

    // note: for v1, we will use this function only for the very first VCF header (which will tell us this is v1)
    int ret = zfile_read_one_section (vb, &vcf_header_section, "vcf_header_section", 
                                      sizeof(SectionHeaderVCFHeader), SEC_VCF_HEADER);
    if (ret == EOF) {
        buf_free (&vcf_header_section);
        return false; // empty file (or in case of split mode - no more components) - not an error
    }

    // handle the GENOZIP header of the VCF header section
    SectionHeaderVCFHeader *header = (SectionHeaderVCFHeader *)vcf_header_section.data;

    ASSERT (BGEN32 (header->h.compressed_offset) == crypt_padded_len (sizeof(SectionHeaderVCFHeader)), "Error: invalid VCF header's header size: header->h.compressed_offset=%u, expecting=%u", BGEN32 (header->h.compressed_offset), (unsigned)sizeof(SectionHeaderVCFHeader));

    // in split mode - we open the output VCF file of the component
    if (flag_split) {
        ASSERT0 (!vb->vcf_file, "Error: not expecting vb->vcf_file to be open already in split mode");
        vb->vcf_file = file_open (header->vcf_filename, WRITE, VCF);
    }

    bool first_vcf = !buf_is_allocated (&global_vcf_header_line);

    uint32_t max_lines_per_vb = BGEN32 (header->max_lines_per_vb);

    if (first_vcf || !flag_concat_mode) {
        vb->z_file->num_lines_concat     = vb->vcf_file->num_lines_concat     = BGEN64 (header->num_lines);
        vb->z_file->vcf_data_size_concat = vb->vcf_file->vcf_data_size_concat = BGEN64 (header->vcf_data_size);

        global_max_lines_per_vb = max_lines_per_vb;
    }
    else {
        // 2nd+ concatenated file is expected to have the same number of max_lines_per_vb
        ASSERT (max_lines_per_vb == global_max_lines_per_vb, 
                "ERROR: invalidly, a 2nd+ concatenated file has max_lines_per_vb=%u different than global_max_lines_per_vb=%u",
                max_lines_per_vb, global_max_lines_per_vb)
    }

    if (flag_split) {
        vb->vcf_file->has_md5 = !md5_is_zero (header->md5_hash_single); // has_md5 iff not 0. note: a chance of 1 in about 10^38 that we will get all-0 by chance in which case will won't perform the md5 comparison
        *digest = header->md5_hash_single; // override md5 from genozip header
    }
        
    // now get the text of the VCF header itself
    static Buffer vcf_header_buf = EMPTY_BUFFER;
    zfile_uncompress_section (vb, header, &vcf_header_buf, "vcf_header_buf", SEC_VCF_HEADER);

    bool can_concatenate = vcf_header_set_globals (vb, vb->z_file->name, &vcf_header_buf);
    if (!can_concatenate) {
        buf_free (&vcf_header_section);
        buf_free (&vcf_header_buf);
        return false;
    }

    // write vcf header if not in concat mode, or, in concat mode, we write the vcf header, only for the first genozip file
    if (first_vcf || !flag_concat_mode)
        vcffile_write_to_disk (vb->vcf_file, &vcf_header_buf);
    
    // if we didn't write the header (bc 2nd+ file in concat mode) - just account for it in MD5 if needed (this is normally done by vcffile_write_to_disk())
    else if (vb->vcf_file->has_md5)
        md5_update (&vb->vcf_file->md5_ctx_concat, vcf_header_buf.data, vcf_header_buf.len, true);

    buf_free (&vcf_header_section);
    buf_free (&vcf_header_buf);

    return true;
}

#define V1_VCF_HEADER // select the vcf_header functions of v1.c
#include "v1.c"
