// ------------------------------------------------------------------
//   vcf_header.c
//   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "vczip.h"

unsigned global_num_samples = 0; // a global param - assigned once before any thread is crated, and never changes - so no thread safety issue
static Buffer global_vcf_header_line = EMPTY_BUFFER; // header line of first VCF file read - use to compare to subsequent files to make sure they have the same header during concat

static bool vcf_header_set_globals(VariantBlock *vb, const char *filename, Buffer *vcf_header, bool concat_mode)
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
                buf_copy (vb, &global_vcf_header_line, vcf_header, i, vcf_header->len - i);
                global_vcf_header_line.len = vcf_header->len - i;
                vcf_header_line_filename = filename;
            }

            // subsequent files - if we're in concat mode just compare to make sure the header is the same
            else if (concat_mode && 
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

            ASSERT (tab_count != 8, "Error: invalid VCF file - field header line contains a FORMAT field but no samples%s", "");

            ASSERT (tab_count == 7, "Error: invalid VCF file - field header line contains only %d fields, expecting at least 8", tab_count+1);

            global_num_samples = 0; // a VCF file without samples - that's perfectly fine
            return true; 
        }
    }

    ABORT ("Error: invalid VCF file - it does not contain a field header line; tab_count=%u", tab_count+1);
    return false; // avoid complication warnings
}

// reads VCF header and writes its compressed form to the VCZ file. returns num_samples.
bool vcf_header_vcf_to_vcz (VariantBlock *vb, 
                            bool concat_mode, unsigned *line_i, Buffer **first_data_line)
{    
    static Buffer line            = EMPTY_BUFFER; // serves to read the header, then its the first line in the data, and again the header when starting the next vcf file
    static Buffer vcf_header_text = EMPTY_BUFFER;

    // line might be used as the first line, so it cannot be free at the end of this function
    // however, by the time we come here again, at the next VCF file, it is no longer needed
    if (buf_is_allocated (&line)) buf_free (&line); 

    *first_data_line = NULL;

    const unsigned INITIAL_BUF_SIZE = 65536;

    buf_alloc (vb, &vcf_header_text, INITIAL_BUF_SIZE, 0, "vcf_header_do", 0);

    while (1) 
    {
        bool success = vcffile_get_line(vb, *line_i + 1, &line);
        if (!success) break; // end of header - no data lines in this VCF file

        (*line_i)++;

        // case : we reached the end of the header - we exit here

        if (line.data[0] != '#') { // end of header - we have read the first data line
            line.data[line.len] = '\0'; // remove newline - data line reader is expecting lines without the newline

            *first_data_line = &line;

            break;
        }

        buf_alloc (vb, &vcf_header_text, line.len + vcf_header_text.len + 1, 2, "vcf_header_do", 1); // +1 for terminating \0

        memcpy (&vcf_header_text.data[vcf_header_text.len], line.data, line.len);
        
        vcf_header_text.len += line.len;
        vcf_header_text.data[vcf_header_text.len] = '\0';
    }

    // case - vcf header was found 
    if (vcf_header_text.len) {

        bool first_vcf = !buf_is_allocated (&global_vcf_header_line);

        bool can_concatenate = vcf_header_set_globals(vb, vb->vcf_file->name, &vcf_header_text, concat_mode);
        if (!can_concatenate) { 
            // this is the second+ file in a concatenation list, but its samples are incompatible
            buf_free (&vcf_header_text);
            return false;
        }

        // in concat mode, we write the header to the dv file, only for the first vcf file
        if (vb->z_file && (!concat_mode || first_vcf)) 
            zfile_write_vcf_header (vb, &vcf_header_text); 

        vb->vcf_file->section_bytes[SEC_VCF_HEADER] = vcf_header_text.len;
        vb->z_file  ->section_bytes[SEC_VCF_HEADER] = vb->z_section_bytes[SEC_VCF_HEADER]; // comes from zfile_compress
    }

    // case : header not found, but data line found
    else 
        ASSERT(*first_data_line, "Error: file has no VCF header%s", "");

    // case : empty file - not an error (caller will see that *first_data_line is NULL)
    buf_free (&vcf_header_text);

    return true; // everything's good
}

bool vcf_header_vcz_to_vcf (VariantBlock *vb, bool concat_mode)
{
    static Buffer compressed_vcf_section = EMPTY_BUFFER;

    int ret = zfile_read_one_section(vb, &compressed_vcf_section, sizeof(SectionHeaderVCFHeader), SEC_VCF_HEADER, true);
    if (ret == EOF) {
        buf_free (&compressed_vcf_section);
        return false; // empty file - not an error
    }

    // handle the VCZ header of the VCF header section
    SectionHeaderVCFHeader *header = (SectionHeaderVCFHeader *)compressed_vcf_section.data;

    ASSERT (header->vczip_version == VCZIP_VERSION, "Error: file version %u is newer than the latest version supported %u. Please upgrade.",
            header->vczip_version, VCZIP_VERSION);

    ASSERT (header->compression_alg == COMPRESSION_ALG_BZLIB, "Unrecognized compression_alg=%u", header->compression_alg);

    ASSERT (ENDN32 (header->h.compressed_offset) == sizeof(SectionHeaderVCFHeader), "Error: invalid VCF header's header size: header->h.compressed_offset=%u, expecting=%u", ENDN32 (header->h.compressed_offset), (unsigned)sizeof(SectionHeaderVCFHeader));

    unsigned data_compressed_len = compressed_vcf_section.len - sizeof(SectionHeaderVCFHeader);
    ASSERT (data_compressed_len == ENDN32 (header->h.data_compressed_len), "Error: failed to read VCF header's data. Read %u bytes, expecting %u bytes", data_compressed_len, header->h.data_compressed_len);

    vb->z_file->num_lines     = vb->vcf_file->num_lines     = ENDN64 (header->num_lines);
    vb->z_file->vcf_data_size = vb->vcf_file->vcf_data_size = ENDN64 (header->vcf_data_size);

    // now get the text of the VCF header itself
    static Buffer vcf_header_buf = EMPTY_BUFFER;
    zfile_uncompress_section(vb, header, &vcf_header_buf, SEC_VCF_HEADER);

    bool first_vcf = !buf_is_allocated (&global_vcf_header_line);

    bool can_concatenate = vcf_header_set_globals (vb, vb->z_file->name, &vcf_header_buf, concat_mode);
    if (!can_concatenate) {
        buf_free(&compressed_vcf_section);
        buf_free(&vcf_header_buf);
        return false;
    }

    // in concat mode, we write the vcf header, only for the first dv file
    if (first_vcf || !concat_mode)
        vcffile_write_to_disk (vb->vcf_file, &vcf_header_buf);

    buf_free(&compressed_vcf_section);
    buf_free(&vcf_header_buf);

    return true;
}

// returns the the VCF header section's header from a VCZ file
bool vcf_header_get_vcf_header (File *z_file, SectionHeaderVCFHeader *vcf_header_header)
{
    int bytes = fread ((char*)vcf_header_header, 1, sizeof(SectionHeaderVCFHeader), z_file->file);
    
    return (bytes == sizeof(SectionHeaderVCFHeader) && ENDN32 (vcf_header_header->magic) == VCZIP_MAGIC);
}
