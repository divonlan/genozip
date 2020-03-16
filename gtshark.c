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

static void gtshark_create_vcf_file (VariantBlock *vb, const Buffer *section_data, unsigned sb_i,
                                     const char *gtshark_vcf_name)
{
    FILE *file = fopen (gtshark_vcf_name, "wb");
    ASSERT (file, "Error: failed to create temporary file %s", gtshark_vcf_name);

    fprintf (file, "##fileformat=VCFv4.2\n");
    fprintf (file, "##contig=<ID=1,length=1>\n");
    fprintf (file, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
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
        fprintf (file, "1\t1\t.\t.\t.\t.\t.\t.\tGT");
        unsigned num_exceptions_in_line = 0; 
        uint16_t last_exception_ht_i = 0;
        for (unsigned ht_i=0; ht_i < num_haplotypes; ht_i++) {
            char c = section_data->data[ht_i * vb->num_lines + vb_line_i];
            
            // case: gtshark can't handle alleles>2 (natively, it splits them to several lines).
            // we put this allele in the exception list and change it to '0' for gtshark.
            if (c > '2') {
                if (!num_exceptions_in_line) { // first exception for this line 

                    if (vb->gtshark_exceptions_ht_i.len) { // we have exception for a previous line - terminate it
                        *NEXTENT (uint16_t, &vb->gtshark_exceptions_ht_i) = 0;
                        *NEXTENT (char, &vb->gtshark_exceptions_allele) = 0;
                    }
                    buf_alloc_more (vb, &vb->gtshark_exceptions_line_i, 1, uint32_t, 2);
                    *NEXTENT (uint32_t, &vb->gtshark_exceptions_line_i) = vb_line_i;
                }

                buf_alloc_more (vb, &vb->gtshark_exceptions_ht_i, 2, uint16_t, 2); // room for terminator too
                *NEXTENT (uint16_t, &vb->gtshark_exceptions_ht_i) = ht_i - last_exception_ht_i; // delta encoding
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
//        if (num_exceptions_in_line) // DEBUG
//            printf("%u exceptions: vb_i=%u vb_line_i=%u sb_i=%u \n", // DEBUG
//                   num_exceptions_in_line, vb->variant_block_i, vb_line_i, sb_i);
    }

    fclose (file);

    return;
error:
    if (file) {
        fclose (file);
        file_remove (gtshark_vcf_name);
    }
    my_exit();
}

static void gtshark_run_gtshark (VariantBlock *vb, unsigned sb_i,
                                 const char *gtshark_vcf_name, const char *gtshark_db_name,
                                 const char *gtshark_db_db_name, const char *gtshark_db_gt_name)
{
#ifndef _WIN32
    int stdout_pipe[2];
    int ret = pipe (stdout_pipe);
    ASSERT (!ret, "Error: failed to created stdout_pipe pipe: %s", strerror (errno));
    
    pid_t child_pid = fork();
    if (!child_pid) {
        // redirect stdout a pipe
        dup2 (stdout_pipe[1], STDOUT_FILENO);
        close (stdout_pipe[0]); // close reading side of the pipes    

        // i am the child
        const char *argv[30];
        int argc = 0;
        argv[argc++] = "gtshark";
        argv[argc++] = "compress-db";
        argv[argc++] = gtshark_vcf_name;
        argv[argc++] = gtshark_db_name;
        argv[argc] = NULL;
        
        execvp ("gtshark", (char * const *)argv); // if successful, doesn't return
        ABORT ("Error executing gtshark: %s", strerror (errno)); 
    }
    // close writing side of the pipes
    close (stdout_pipe[1]);

    // read pipe (up to 10000 characters)
    #define PIPE_MAX_BYTES 10000
    char stdout_data[PIPE_MAX_BYTES];
    char *stdout_next = stdout_data;
//printf ("BEFORE READ FROM PIPE\n"); // DEBUG
    // get stdout output of the gtshark process
    while (read (stdout_pipe[0], stdout_next, 1) > 0 && (stdout_next - stdout_data < PIPE_MAX_BYTES-1)) stdout_next++;
    *stdout_next = 0; // string terminator
    close (stdout_pipe[0]);
//printf ("STDOUT:\n%s\n", stdout_data); // DEBUG

    ASSERT (!strstr (stdout_data, "error") && !strstr (stdout_data, "E::") & !strstr (stdout_data, "Invalid"),
            "Error: gtshark failed for vb_i=%u sb_i=%u. Here is its stdout:\n%s", vb->variant_block_i, sb_i, stdout_data);

    // wait for gtshark to complete
    int wstatus;
    waitpid (child_pid, &wstatus, 0); // I am the parent - wait for child, so that the terminal doesn't print the prompt until the child is done
    if (wstatus) fprintf (stderr, "gtshark exited\n"); // DEBUG

    ASSERT (!WEXITSTATUS (wstatus), 
            "Error: gtshark exited with status=%u for vb_i=%u sb_i=%u. Here is its stdout:\n%s", 
            WEXITSTATUS (wstatus), vb->variant_block_i, sb_i, stdout_data);

    ASSERT (WIFEXITED (wstatus), 
            "Error: gtshark failed to exist normally for vb_i=%u sb_i=%u. Here is its stdout:\n%s", 
            WEXITSTATUS (wstatus), vb->variant_block_i, sb_i, stdout_data);
#endif
//printf ("AFTER WAIT\n");
    // read both gtshark output files
    file_get_file (vb, gtshark_db_db_name, &vb->gtshark_db_db_data, "gtshark_db_db_data", vb->variant_block_i);
    file_get_file (vb, gtshark_db_gt_name, &vb->gtshark_db_gt_data, "gtshark_db_gt_data", vb->variant_block_i);
    
    file_remove (gtshark_vcf_name);
    file_remove (gtshark_db_db_name);
    file_remove (gtshark_db_gt_name);
}

void gtshark_compress (VariantBlock *vb, const Buffer *section_data, unsigned sb_i)
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

    gtshark_run_gtshark (vb, sb_i, gtshark_vcf_name, gtshark_db_name, gtshark_db_db_name, gtshark_db_gt_name);

    free (gtshark_vcf_name);
    free (gtshark_db_name);
    free (gtshark_db_db_name);
    free (gtshark_db_gt_name);
}

void gtshark_copy_to_z_data (VariantBlock *vb, Buffer *z_data, uint32_t compressed_offset)
{
    uint32_t save_len = z_data->len; // save

    z_data->len += compressed_offset;
    buf_add (z_data, vb->gtshark_exceptions_line_i.data, vb->gtshark_exceptions_line_i.len * sizeof (uint32_t));
    buf_add (z_data, vb->gtshark_exceptions_ht_i.data, vb->gtshark_exceptions_ht_i.len * sizeof (uint16_t));
    buf_add (z_data, vb->gtshark_exceptions_allele.data, vb->gtshark_exceptions_allele.len);
    buf_add (z_data, vb->gtshark_db_db_data.data, vb->gtshark_db_db_data.len);
    buf_add (z_data, vb->gtshark_db_gt_data.data, vb->gtshark_db_gt_data.len);
    
    z_data->len = save_len; // restore
}