// ------------------------------------------------------------------
//   header_vcf.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "vcf_private.h"
#include "zfile.h"
#include "txtfile.h"
#include "vblock.h"
#include "crypt.h"
#include "version.h"
#include "endianness.h"
#include "file.h"
#include "dispatcher.h"
#include "txtfile.h"
#include "strings.h"

// Globals
uint32_t vcf_num_samples           = 0; // number of samples in the file
uint32_t vcf_num_displayed_samples = 0; // PIZ only: number of samples to be displayed - might be less that vcf_num_samples if --samples is used
static Buffer vcf_header_line = EMPTY_BUFFER;  // header line of first VCF file read - use to compare to subsequent files to make sure they have the same header during bound

void vcf_header_initialize (void)
{
    vcf_num_samples           = 0;
    vcf_num_displayed_samples = 0;
    buf_free (&vcf_header_line);
}

bool vcf_inspect_txt_header (Buffer *txt_header)
{
    return vcf_header_set_globals (txt_file->name, txt_header);
}

bool vcf_header_set_globals(const char *filename, Buffer *vcf_header)
{
    static const char *vcf_header_line_filename = NULL; // file from which the header line was taken

    // count tabs in last line which should be the field header line
    unsigned tab_count = 0;
    for (int i=vcf_header->len-1; i >= 0; i--) {
        
        if (vcf_header->data[i] == '\t')
            tab_count++;
        
        // some times files have spaces instead of \t 
        else if (vcf_header->data[i] == ' ') {
            tab_count++;
            while (i >= 1 && vcf_header->data[i-1]==' ') i--; // skip run of spaces
        }

        // if this is the beginning of field header line 
        else if (vcf_header->data[i] == '#' && (i==0 || vcf_header->data[i-1] == '\n' || vcf_header->data[i-1] == '\r')) {
        
            // if first vcf file - copy the header to the global
            if (!buf_is_allocated (&vcf_header_line)) {
                buf_copy (evb, &vcf_header_line, vcf_header, 1, i, vcf_header->len - i, "vcf_header_line", 0);
                vcf_header_line_filename = filename;
            }

            // ZIP only: subsequent files - if we're in bound mode just compare to make sure the header is the same
            else if (flag_bind && 
                     (vcf_header->len-i != vcf_header_line.len || memcmp (vcf_header_line.data, &vcf_header->data[i], vcf_header_line.len))) {

                fprintf (stderr, "%s: skipping %s: it has a different VCF header line than %s, see below:\n"
                                 "========= %s =========\n"
                                 "%.*s"
                                 "========= %s ==========\n"
                                 "%.*s"
                                 "=======================================\n", 
                         global_cmd, filename, vcf_header_line_filename,
                         vcf_header_line_filename, (uint32_t)vcf_header_line.len, vcf_header_line.data,
                         filename, (uint32_t)vcf_header->len-i, &vcf_header->data[i]);
                return false;
            }

            //count samples
            vcf_num_samples = (tab_count >= 9) ? tab_count-8 : 0; 
            // note: a VCF file without samples may or may not have a "FORMAT" in the header, i.e. tab_count==7 or 8 (8 or 9 fields).
            // however, even if it has a FORMAT in the header, it won't have a FORMAT column in the data
            
            vcf_num_displayed_samples = vcf_num_samples;

            ASSERT (tab_count >= 7, "Error: invalid VCF file - field header line contains only %d fields, expecting at least 8", tab_count+1);

            // if --samples is used, update vcf_header and vcf_num_displayed_samples
            if (flag_samples) samples_digest_vcf_header (vcf_header);

            return true; 
        }
    }

    ABORT ("Error: invalid VCF file - it does not contain a field header line; tab_count=%u", tab_count+1);
    return false; // avoid complication warnings
}

// genocat: remove FORMAT and sample names from the vcf header line, in case of --drop-genotypes
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

// genocat: remove all lines but the last from the vcf header, in case of --header-one
void vcf_header_keep_only_last_line (Buffer *vcf_header_buf)
{
    for (int i=vcf_header_buf->len-2; i >= 0; i--) // -2 - skip last newline
        if (vcf_header_buf->data[i] == '\n') {
            vcf_header_buf->len = vcf_header_buf->len - i - 1;
            memcpy (vcf_header_buf->data, &vcf_header_buf->data[i+1], vcf_header_buf->len);
            break;
        }
}

uint32_t vcf_header_get_num_samples (void)
{
    return vcf_num_samples;
}
