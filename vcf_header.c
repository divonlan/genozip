// ------------------------------------------------------------------
//   vcf_header.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "zfile.h"
#include "vcffile.h"
#include "vcf_header.h"
#include "vb.h"
#include "crypt.h"
#include "version.h"
#include "endianness.h"
#include "file.h"
#include "samples.h"

unsigned global_num_samples              = 0; // number of samples in the file
unsigned global_number_displayed_samples = 0; // PIZ only: number of samples to be displayed - might be less that global_num_samples if --samples is used
unsigned global_max_memory_per_vb        = DEFAULT_MAX_MEMORY_PER_VB; // ZIP only: used for reading VCF data

Buffer global_vcf_header_line = EMPTY_BUFFER; // header line of first VCF file read - use to compare to subsequent files to make sure they have the same header during concat

// global - this names go into the dictionary names on disk. to preserve backward compatability, they should not be changed.
const char *vcf_field_names[] = { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT" };

void vcf_header_initialize()
{
    global_num_samples              = 0;
    global_number_displayed_samples = 0;
    global_max_memory_per_vb        = DEFAULT_MAX_MEMORY_PER_VB;
    buf_free (&global_vcf_header_line);
}

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
            else if (flag_concat && 
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
            global_number_displayed_samples = global_num_samples;

            ASSERT0 (tab_count != 8, "Error: invalid VCF file - field header line contains a FORMAT field but no samples");

            ASSERT (tab_count >= 7, "Error: invalid VCF file - field header line contains only %d fields, expecting at least 8", tab_count+1);

            // if --samples is used, update vcf_header and global_number_displayed_samples
            if (flag_samples) samples_digest_vcf_header (vcf_header);

            return true; 
        }
    }

    ABORT ("Error: invalid VCF file - it does not contain a field header line; tab_count=%u", tab_count+1);
    return false; // avoid complication warnings
}

// reads VCF header and writes its compressed form to the GENOZIP file. returns num_samples.
bool vcf_header_vcf_to_genozip (uint32_t *vcf_line_i)
{    
    evb->z_file->disk_at_beginning_of_this_vcf_file = evb->z_file->disk_so_far;

    bool is_first_vcf = !buf_is_allocated (&global_vcf_header_line); 

    vcffile_read_vcf_header (is_first_vcf); // reads into evb->vcf_data and evb->num_lines
    
    *vcf_line_i += evb->num_lines;

    // case - vcf header was found 
    if (evb->vcf_data.len) {

        bool can_concatenate = vcf_header_set_globals(evb, evb->vcf_file->name, &evb->vcf_data);
        if (!can_concatenate) { 
            // this is the second+ file in a concatenation list, but its samples are incompatible
            buf_free (&evb->vcf_data);
            return false;
        }

        if (evb->z_file) zfile_write_vcf_header (&evb->vcf_data, is_first_vcf); // we write all headers in concat mode too, to support --split

        evb->vcf_file->section_bytes[SEC_VCF_HEADER] = evb->vcf_data.len;
        evb->z_file  ->section_bytes[SEC_VCF_HEADER] = evb->z_section_bytes[SEC_VCF_HEADER]; // comes from zfile_compress
        evb->z_file  ->num_sections [SEC_VCF_HEADER]++;
        evb->z_file  ->num_vcf_components_so_far++; // when compressing
    }

    if (flag_md5) {
        if (flag_concat && is_first_vcf) 
            md5_update (&evb->z_file->md5_ctx_concat, evb->vcf_data.data, evb->vcf_data.len);

        md5_update (&evb->z_file->md5_ctx_single, evb->vcf_data.data, evb->vcf_data.len);
    }

    buf_free (&evb->vcf_data);
    
    return true; // everything's good
}

// genocat: remove trip the vcf header line, in case of --drop-genotypes
void vcf_header_trim_header_line (Buffer *vcf_header_buf)
{
    static const char *standard = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    for (int i=vcf_header_buf->len-2; i >= 0; i--) // -2 - skip last newline
        if (vcf_header_buf->data[i] == '\n') { 
            if (!memcmp (&vcf_header_buf->data[i+1], standard, MIN (strlen (standard), vcf_header_buf->len-(i+1)))) {
                vcf_header_buf->len = i + strlen(standard) + 2; // fixed length of standard VCF header up to INFO
                vcf_header_buf->data[vcf_header_buf->len-1] = '\n';   
            }                  
            return;
        }
    // no newline found - nothing to do
}

// returns offset of header within data, EOF if end of file
bool vcf_header_genozip_to_vcf (Md5Hash *digest) // NULL if we're just skipped this header (2nd+ header in concatenated file)
{
    evb->z_file->disk_at_beginning_of_this_vcf_file = evb->z_file->disk_so_far;
    static Buffer vcf_header_section = EMPTY_BUFFER;

    int header_offset = zfile_read_one_section (evb, 0, &vcf_header_section, "vcf_header_section", 
                                                sizeof(SectionHeaderVCFHeader), SEC_VCF_HEADER);
    if (header_offset == EOF) {
        buf_free (&vcf_header_section);
        return false; // empty file (or in case of split mode - no more components) - not an error
    }

    // handle the GENOZIP header of the VCF header section
    SectionHeaderVCFHeader *header = (SectionHeaderVCFHeader *)vcf_header_section.data;

    ASSERT (!digest || BGEN32 (header->h.compressed_offset) == crypt_padded_len (sizeof(SectionHeaderVCFHeader)), 
            "Error: invalid VCF header's header size: header->h.compressed_offset=%u, expecting=%u", BGEN32 (header->h.compressed_offset), (unsigned)sizeof(SectionHeaderVCFHeader));

    // in split mode - we open the output VCF file of the component
    if (flag_split) {
        ASSERT0 (!evb->vcf_file, "Error: not expecting evb->vcf_file to be open already in split mode");
        evb->vcf_file = file_open (header->vcf_filename, WRITE, VCF);
        evb->z_file->vcf_data_size_single = BGEN64 (header->vcf_data_size);
    }

    bool first_vcf = !buf_is_allocated (&global_vcf_header_line);

    evb->vcf_file->max_lines_per_vb = BGEN32 (header->max_lines_per_vb);

    if (first_vcf || flag_split) {
        evb->z_file->num_lines_concat     = evb->vcf_file->num_lines_concat     = BGEN64 (header->num_lines);
        evb->z_file->vcf_data_size_concat = evb->vcf_file->vcf_data_size_concat = BGEN64 (header->vcf_data_size);
    }

    if (flag_split) *digest = header->md5_hash_single; // override md5 from genozip header
        
    // now get the text of the VCF header itself
    static Buffer vcf_header_buf = EMPTY_BUFFER;
    zfile_uncompress_section (evb, header, &vcf_header_buf, "vcf_header_buf", SEC_VCF_HEADER);

    bool can_concatenate = vcf_header_set_globals (evb, evb->z_file->name, &vcf_header_buf);
    if (!can_concatenate) {
        buf_free (&vcf_header_section);
        buf_free (&vcf_header_buf);
        return false;
    }

    if (flag_drop_genotypes) vcf_header_trim_header_line (&vcf_header_buf);

    // write vcf header if not in concat mode, or, in concat mode, we write the vcf header, only for the first genozip file
    if ((first_vcf || flag_split) && !flag_no_header)
        vcffile_write_to_disk (evb->vcf_file, &vcf_header_buf);
    
    buf_free (&vcf_header_section);
    buf_free (&vcf_header_buf);

    evb->z_file->num_vcf_components_so_far++;

    return true;
}

#define V1_VCF_HEADER // select the vcf_header functions of v1.c
#include "v1.c"
