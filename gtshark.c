// ------------------------------------------------------------------
//   gtshark.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
#include "gtshark.h"

#include <errno.h>
#include <sys/types.h>
#include <fcntl.h> 
#include <unistd.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif
#include "vb.h"
#include "buffer.h"
#include "file.h"
#include "endianness.h"
#include "stream.h"

static void gtshark_create_vcf_file (VariantBlock *vb, const Buffer *section_data, unsigned sb_i,
                                     const char *gtshark_vcf_name)
{
    FILE *file = fopen (gtshark_vcf_name, "wb");
    ASSERT (file, "Error: failed to create temporary file %s", gtshark_vcf_name);

    fprintf (file, "##fileformat=VCFv4.2\n");
    fprintf (file, "##contig=<ID=Z>\n");
    fprintf (file, "##FORMAT=<ID=GT>\n");
    fprintf (file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    unsigned num_haplotypes = vb->ploidy * vb_num_samples_in_sb (vb, sb_i); 
    for (unsigned i=0; i < num_haplotypes; i++)
        fprintf (file, "\t%u", i+1);
    fprintf (file, "\n");

    ASSERTGOTO (section_data->len == num_haplotypes * vb->num_lines, "Error: unexpected section_data->len=%u", section_data->len);

    // initialize allocation for exceptions
    buf_alloc (vb, &vb->gtshark_exceptions_line_i, MAX (vb->num_lines/100, 100) * sizeof(uint32_t), 1, "gtshark_exceptions_line_i", sb_i);
    buf_alloc (vb, &vb->gtshark_exceptions_ht_i, MAX (vb->num_lines/100, 100) * (global_num_samples / 3) * sizeof(uint16_t), 1, "gtshark_exceptions_ht_i", sb_i);
    buf_alloc (vb, &vb->gtshark_exceptions_allele, MAX (vb->num_lines/100, 100) * (global_num_samples / 3), 1, "gtshark_exceptions_allele", sb_i);
    
    for (unsigned vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        #define GTSHARK_CHROM_ID "Z"
        #define GTSHARK_VCF_LINE_VARDATA GTSHARK_CHROM_ID "\t.\t.\t.\t.\t.\t.\t.\tGT"

        fprintf (file, GTSHARK_VCF_LINE_VARDATA);
        unsigned num_exceptions_in_line = 0; 
        uint16_t last_exception_ht_i = 0;
        for (unsigned ht_i=0; ht_i < num_haplotypes; ht_i++) {
            char c = section_data->data[ht_i * vb->num_lines + vb_line_i];
            
            // case: gtshark can't handle alleles>2 (natively, it splits them to several lines).
            // we put this allele in the exception list and change it to '0' for gtshark.
            if (c < '0' || c > '2') {
                if (!num_exceptions_in_line) { // first exception for this line 
                    buf_alloc_more (vb, &vb->gtshark_exceptions_line_i, 1, uint32_t, 2);
                    *NEXTENT (uint32_t, &vb->gtshark_exceptions_line_i) = BGEN32 (vb_line_i);
                }

                buf_alloc_more (vb, &vb->gtshark_exceptions_ht_i, 2, uint16_t, 2); // room for terminator too
                *NEXTENT (uint16_t, &vb->gtshark_exceptions_ht_i) = BGEN16 (ht_i - last_exception_ht_i); // delta encoding
                last_exception_ht_i = ht_i;
                
                buf_alloc_more (vb, &vb->gtshark_exceptions_allele, 2, char, 2);   // room for terminator too
                *NEXTENT (char, &vb->gtshark_exceptions_allele) = c;

                num_exceptions_in_line++; 
                fprintf (file, "\t0");    
            }
            else 
                fprintf (file, "\t%c", c);
        }
        fprintf (file, "\n");

        if (num_exceptions_in_line) { // we have exceptions for this line - terminate the ht_i and allele arrays
            *NEXTENT (uint16_t, &vb->gtshark_exceptions_ht_i) = 0;
            *NEXTENT (char, &vb->gtshark_exceptions_allele) = 0;
        }
    }

    fclose (file);

    return;
error:
    if (file) {
        fclose (file);
        file_remove (gtshark_vcf_name, true);
    }
    my_exit();
}

#define PIPE_MAX_BYTES 32768

// check stdout / stderr output of the gtshark process for errors
static void gtshark_check_pipe_for_errors (char *data, FILE *fp, uint32_t vb_i, uint32_t sb_i, bool is_stderr) 
{
    unsigned bytes = fread (data, 1, PIPE_MAX_BYTES-1, fp); // -1 to leave room for \0
    data[bytes] = '\0';

    ASSERT (!strstr (data, "error")   && 
            !strstr (data, "E::")     && 
            !strstr (data, "Invalid") &&
            !strstr (data, "throw"),
            "Error: gtshark failed for vb_i=%u sb_i=%u. Here is its %s:\n%s\n", 
            vb_i, sb_i, is_stderr ? "STDERR" : "STDOUT", data);
}

// ZIP & PIZ
static bool gtshark_run (uint32_t vb_i, unsigned sb_i,
                         const char *command, const char *filename_1, const char *filename_2) 
{
    Stream gtshark = stream_create (DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0,
                                    "gtshark", command, filename_1, filename_2, NULL);

    // read pipe (up to 10000 characters)
    char stdout_data[PIPE_MAX_BYTES], stderr_data[PIPE_MAX_BYTES];
    gtshark_check_pipe_for_errors (stdout_data, gtshark.from_stream_stdout, vb_i, sb_i, false);
    gtshark_check_pipe_for_errors (stderr_data, gtshark.from_stream_stderr, vb_i, sb_i, true);

    // wait for gtshark to complete
    int exit_status = stream_wait_for_exit (gtshark);

    stream_close (&gtshark);

#ifndef _WIN32
    ASSERT (!WEXITSTATUS (exit_status), 
            "Error: gtshark exited with status=%u for vb_i=%u sb_i=%u. Here is its STDOUT:\n%s\nHere is the STDERR:\n%s\n", 
            WEXITSTATUS (exit_status), vb_i, sb_i, stdout_data, stderr_data);

    ASSERT (!WIFSIGNALED (exit_status),
            "Error: gtshark process was killed by a signal, it was running for vb_i=%u sb_i=%u. Here is its STDOUT:\n%s\nHere is the STDERR:\n%s\n", 
            vb_i, sb_i, stdout_data, stderr_data);
#endif

    ASSERT (!exit_status, 
            "Error: gtshark failed to exit normally for vb_i=%u sb_i=%u. Here is its STDOUT:\n%s\nHere is the STDERR:\n%s\n", 
            vb_i, sb_i, stdout_data, stderr_data);

    return true; // successfully executed gtshark
}

// ZIP
static void gtshark_run_compress (VariantBlock *vb, unsigned sb_i,
                                  const char *gtshark_vcf_name, const char *gtshark_db_name,
                                  const char *gtshark_db_db_name, const char *gtshark_db_gt_name)
{
    // remove in case of leftovers from previous run
    file_remove (gtshark_db_db_name, true);
    file_remove (gtshark_db_gt_name, true);

    gtshark_run (vb->variant_block_i, sb_i, "compress-db", gtshark_vcf_name, gtshark_db_name);

    // read both gtshark output files
    file_get_file (vb, gtshark_db_db_name, &vb->gtshark_db_db_data, "gtshark_db_db_data", vb->variant_block_i, false);
    file_get_file (vb, gtshark_db_gt_name, &vb->gtshark_db_gt_data, "gtshark_db_gt_data", vb->variant_block_i, false);
    
    file_remove (gtshark_vcf_name, false);
    file_remove (gtshark_db_db_name, false);
    file_remove (gtshark_db_gt_name, false);
}

// ZIP
void gtshark_compress_haplotype_data (VariantBlock *vb, const Buffer *section_data, unsigned sb_i)
{
    char *gtshark_vcf_name = malloc (strlen (z_file->name) + 20);
    sprintf (gtshark_vcf_name, "%s.%u.%u.db.vcf", file_printname (z_file), vb->variant_block_i, sb_i);

    char *gtshark_db_name = malloc (strlen (z_file->name) + 20);
    sprintf (gtshark_db_name, "%s.%u.%u.db", file_printname (z_file), vb->variant_block_i, sb_i);

    char *gtshark_db_db_name = malloc (strlen (z_file->name) + 20);
    sprintf (gtshark_db_db_name, "%s.%u.%u.db_db", file_printname (z_file), vb->variant_block_i, sb_i);

    char *gtshark_db_gt_name = malloc (strlen (z_file->name) + 20);
    sprintf (gtshark_db_gt_name, "%s.%u.%u.db_gt", file_printname (z_file), vb->variant_block_i, sb_i);

    gtshark_create_vcf_file (vb, section_data, sb_i, gtshark_vcf_name);

    gtshark_run_compress (vb, sb_i, gtshark_vcf_name, gtshark_db_name, gtshark_db_db_name, gtshark_db_gt_name);

    free (gtshark_vcf_name);
    free (gtshark_db_name);
    free (gtshark_db_db_name);
    free (gtshark_db_gt_name);
}

// PIZ
static char *gtshark_write_db_file (uint32_t vb_i, uint16_t section_i, const char *file_ext, 
                                    const Buffer *buf)
{
    char *filename = malloc (strlen (file_printname (z_file)) + 50);
    sprintf (filename, "%s.%u.%u.%s", file_printname (z_file), vb_i, section_i, file_ext);

    FILE *file = fopen (filename, "wb");

    size_t bytes_written = fwrite (buf->data, 1, buf->len, file);
    ASSERT (bytes_written == buf->len, 
            "Error in uncompressing vb_i=%u: failed to write file %s - only %u out of %u bytes written: %s",
            vb_i, filename, (unsigned)bytes_written, buf->len, strerror (errno));

    fclose (file);
    return filename;
}

// PIZ
static void gtshark_run_decompress (VariantBlock *vb, unsigned sb_i)
{
    char *gtshark_db_name = malloc (strlen (file_printname (z_file)) + 50);
    sprintf (gtshark_db_name, "%s.%u.%u.db", file_printname (z_file), vb->variant_block_i, sb_i);

    char *gtshark_vcf_name = malloc (strlen (file_printname (z_file)) + 50);
    sprintf (gtshark_vcf_name, "%s.%u.%u.vcf", file_printname (z_file), vb->variant_block_i, sb_i);

    // remove in case of leftovers from previous run
    file_remove (gtshark_vcf_name, true);

    gtshark_run (vb->variant_block_i, sb_i, "decompress-db", gtshark_db_name, gtshark_vcf_name);

    file_get_file (vb, gtshark_vcf_name, &vb->gtshark_vcf_data, "gtshark_vcf_data", vb->variant_block_i, true);

    file_remove (gtshark_vcf_name, false);

    free (gtshark_db_name);
    free (gtshark_vcf_name);
}    

// PIZ: convert the vcf generated by gtshark when decompressing the db, into our haplotype_data
static void gtshark_generate_haplotype_data (VariantBlock *vb, unsigned sb_i)
{
    uint32_t num_lines = vb->num_lines;
    unsigned num_hts = vb->ploidy * vb_num_samples_in_sb (vb, sb_i);

    buf_alloc (vb, &vb->haplotype_sections_data[sb_i], num_lines * num_hts, 1, "haplotype_sections_data", sb_i);
    vb->haplotype_sections_data[sb_i].len = num_lines * num_hts;
    
    uint8_t *haplotype_data = (uint8_t *)vb->haplotype_sections_data[sb_i].data;

    // look for the start of our data - a newline followed by the chrom
    const char *substr = strstr (vb->gtshark_vcf_data.data, "\n" GTSHARK_CHROM_ID);
    ASSERT (substr, "Error: cannot locate start of data within gtshark-produced vcf data for vb_i=%u", vb->variant_block_i);

    // build the transposed matrix from the vcf data - this will include alleles 0, 1 or 2 with higher
    // alleles showing as 0
    
    const char *next = substr + 1; // point to the CHROM = '1'
    unsigned prefix_len = strlen (GTSHARK_VCF_LINE_VARDATA) + 1; // +1 for the tab after GT
    
    for (uint32_t vb_line_i=0; vb_line_i < num_lines; vb_line_i++) {
        ASSERT (*next == 'Z', "Error: expecting vb_line_i=%u to start with 'Z'", vb_line_i);
        next += prefix_len; // skipping the line fields 1-9, arriving at the first haplotype
    
        for (uint32_t ht_i=0; ht_i < num_hts; ht_i++) {
            haplotype_data[ht_i * num_lines + vb_line_i] = *next; // haplotype matrix is transposed
            next += 2; // skip past this ht and also the following \t or \n
        }
    }

    // now enter the higher alleles from the exception list

    const uint32_t *exceptions_line_i_data = FIRSTENT (uint32_t, &vb->gtshark_exceptions_line_i);
    const uint16_t *next_ht_i_delta        = FIRSTENT (uint16_t, &vb->gtshark_exceptions_ht_i);
    const uint16_t *after_ht_i_delta       = AFTERENT (uint16_t, &vb->gtshark_exceptions_ht_i);
    const uint8_t  *next_allele            = FIRSTENT (uint8_t,  &vb->gtshark_exceptions_allele);
    const uint8_t  *after_allele           = AFTERENT (uint8_t,  &vb->gtshark_exceptions_allele);
    
    for (unsigned ex_line_i=0; ex_line_i < vb->gtshark_exceptions_line_i.len; ex_line_i++) {

        uint32_t vb_line_i = BGEN32 (exceptions_line_i_data[ex_line_i]);
        ASSERT (vb_line_i < num_lines, "Error processing exceptions for vb_i=%u: vb_line_i=%u is out of range (num_lines=%u)", 
                vb->variant_block_i, vb_line_i, num_lines);

        uint16_t last_ht_i = 0;
        
        for (unsigned ex_ht_i=0 ; !ex_ht_i || *next_ht_i_delta; ex_ht_i++) { // the list of ht_i's for this line is terminated with a 0

            ASSERT (next_ht_i_delta < after_ht_i_delta, "Error processing exceptions for vb_i=%u: next_ht_i_delta is out of range", vb->variant_block_i);
            ASSERT (next_allele < after_allele, "Error processing exceptions for vb_i=%u: next_allele is out of range", vb->variant_block_i);

            uint16_t delta = BGEN16 (*(next_ht_i_delta++));
            uint16_t ht_i = last_ht_i + delta; // decode delta encoding
            ASSERT (ht_i < num_hts, "Error processing exceptions for vb_i=%u: ht_i=%u is out of range (num_hts=%u)", 
                    vb->variant_block_i, ht_i, num_hts);

            ASSERT (*next_allele <= '0' + MAX_ALLELE_VALUE, "Error processing exceptions for vb_i=%u: allele is out of range (ascii(%u))", vb->variant_block_i, (unsigned)*next_allele);
            haplotype_data[ht_i * num_lines + vb_line_i] = *(next_allele++);
            last_ht_i = ht_i;
        }

        next_ht_i_delta++; // end-of-line separator in both ht_i and allele arrays
        next_allele++;
    }

    buf_free (&vb->gtshark_vcf_data);
}

void gtshark_uncompress_haplotype_data (VariantBlock *vb, unsigned sb_i)
{
//printf ("line_i_len=%u ht_i_len=%u allele_len=%u db_db_len=%u db_gt_len=%u %%-exceptions: %u\n", 
//exceptions_line_i_len, exceptions_ht_i_len, exceptions_allele_len, db_db_data_len, db_gt_data_len,
//(100*(exceptions_line_i_len+exceptions_ht_i_len+exceptions_allele_len) / compressed_len));
    char *filename_db_db = gtshark_write_db_file (vb->variant_block_i, sb_i, "db_db", &vb->gtshark_db_db_data);
    
    char *filename_db_gt = gtshark_write_db_file (vb->variant_block_i, sb_i, "db_gt", &vb->gtshark_db_gt_data);

    gtshark_run_decompress (vb, sb_i);                            

    file_remove (filename_db_db, false);
    file_remove (filename_db_gt, false);
    free (filename_db_db);
    free (filename_db_gt);

    gtshark_generate_haplotype_data (vb, sb_i); 

    // free buffers - we will need them for the next section
    buf_free (&vb->gtshark_exceptions_line_i);
    buf_free (&vb->gtshark_exceptions_ht_i);
    buf_free (&vb->gtshark_exceptions_allele);
    buf_free (&vb->gtshark_db_db_data);
    buf_free (&vb->gtshark_db_gt_data);
}