// ------------------------------------------------------------------
//   codec_gtshark.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
#include <errno.h>
#include <sys/types.h>
#include <fcntl.h> 
#include <unistd.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif
#include "genozip.h"
#include "buffer.h"
#include "file.h"
#include "endianness.h"
#include "stream.h"
#include "vblock.h"
#include "codec.h"
#include "dict_id.h"
#include "strings.h"
#include "piz.h"

// -------------
// ZIP & PIZ
// -------------

#define PIPE_MAX_BYTES 32768

// check stdout / stderr output of the gtshark process for errors
static void gtshark_check_pipe_for_errors (char *data, FILE *fp, uint32_t vb_i, bool is_stderr) 
{
    unsigned bytes = fread (data, 1, PIPE_MAX_BYTES-1, fp); // -1 to leave room for \0
    data[bytes] = '\0';

    ASSERT (!strstr (data, "error")   && 
            !strstr (data, "E::")     && 
            !strstr (data, "Invalid") &&
            !strstr (data, "throw"),
            "Error: gtshark failed for vb_i=%u. Here is its %s:\n%s\n", 
            vb_i, is_stderr ? "STDERR" : "STDOUT", data);
}

static bool codec_gtshark_run (uint32_t vb_i, const char *command, 
                               const char *filename_1, const char *filename_2) 
{
    StreamP gtshark = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 0,
                                     "To use the --gtshark option",
                                     "gtshark", command, filename_1, filename_2, NULL);

    // read pipe (up to 10000 characters)
    char stdout_data[PIPE_MAX_BYTES], stderr_data[PIPE_MAX_BYTES];
    gtshark_check_pipe_for_errors (stdout_data, stream_from_stream_stdout (gtshark), vb_i, false);
    gtshark_check_pipe_for_errors (stderr_data, stream_from_stream_stderr (gtshark), vb_i, true);

    int exit_status = stream_close (&gtshark, STREAM_WAIT_FOR_PROCESS);  

#ifndef _WIN32
    ASSERT (!WEXITSTATUS (exit_status), 
            "Error: gtshark exited with status=%u for vb_i=%u. Here is its STDOUT:\n%s\nHere is the STDERR:\n%s\n", 
            WEXITSTATUS (exit_status), vb_i, stdout_data, stderr_data);

    ASSERT (!WIFSIGNALED (exit_status),
            "Error: gtshark process was killed by a signal, it was running for vb_i=%u. Here is its STDOUT:\n%s\nHere is the STDERR:\n%s\n", 
            vb_i, stdout_data, stderr_data);
#endif

    ASSERT (!exit_status, 
            "Error: gtshark failed to exit normally for vb_i=%u. Here is its STDOUT:\n%s\nHere is the STDERR:\n%s\n", 
            vb_i, stdout_data, stderr_data);

    return true; // successfully executed gtshark
}

// -------------
// ZIP side
// -------------

void codec_gtshark_comp_init (VBlock *vb)
{
    vb->ht_matrix_ctx    = mtf_get_ctx (vb, dict_id_FORMAT_GT_HT);
    vb->ht_matrix_ctx->lcodec = CODEC_GTSHARK; // this will trigger codec_gtshark_compress even though the section is not written to the file
    
    vb->gtshark_db       = mtf_get_ctx (vb, dict_id_FORMAT_GT_SHARK_DB);
    vb->gtshark_gt       = mtf_get_ctx (vb, dict_id_FORMAT_GT_SHARK_GT);
    vb->gtshark_x_line   = mtf_get_ctx (vb, dict_id_FORMAT_GT_SHARK_X_LINE);
    vb->gtshark_x_ht     = mtf_get_ctx (vb, dict_id_FORMAT_GT_SHARK_X_HT);
    vb->gtshark_x_allele = mtf_get_ctx (vb, dict_id_FORMAT_GT_SHARK_X_ALLELE);
}

static void codec_gtshark_create_vcf_file (VBlock *vb, const char *gtshark_vcf_name)
{
    FILE *file = fopen (gtshark_vcf_name, "wb");
    ASSERT (file, "Error: failed to create temporary file %s", gtshark_vcf_name);

    fprintf (file, "##fileformat=VCFv4.2\n");
    fprintf (file, "##contig=<ID=Z>\n");
    fprintf (file, "##FORMAT=<ID=GT>\n");

    #define GTSHARK_NUM_HT_PER_LINE "##num_haplotypes_per_line="
    fprintf (file, GTSHARK_NUM_HT_PER_LINE "%u\n", vb->num_haplotypes_per_line);
    fprintf (file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for (unsigned i=0; i < vb->num_haplotypes_per_line; i++)
        fprintf (file, "\t%u", i+1);
    fprintf (file, "\n");

    // initialize allocation for exceptions
    Buffer *x_line_buf   = &vb->gtshark_x_line->local;
    Buffer *x_ht_buf     = &vb->gtshark_x_ht->local;
    Buffer *x_allele_buf = &vb->gtshark_x_allele->local;

    buf_alloc (vb, x_line_buf,   MAX (vb->lines.len/100, 100) * sizeof(uint32_t), 1, "context->local", 0);
    buf_alloc (vb, x_ht_buf,     MAX (vb->lines.len/100, 100) * (vb->num_haplotypes_per_line / 3) * sizeof(uint16_t), 1, "context->local", 0);
    buf_alloc (vb, x_allele_buf, MAX (vb->lines.len/100, 100) * (vb->num_haplotypes_per_line / 3), 1, "context->local", 0);
    
    ARRAY (char, ht_matrix, vb->ht_matrix_ctx->local);

    for (unsigned vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        #define GTSHARK_CHROM_ID "Z"
        #define GTSHARK_VCF_LINE_VARDATA GTSHARK_CHROM_ID "\t.\t.\t.\t.\t.\t.\t.\tGT"

        fprintf (file, GTSHARK_VCF_LINE_VARDATA);
        unsigned num_exceptions_in_line = 0; 
        uint16_t last_exception_ht_i = 0;
        for (unsigned ht_i=0; ht_i < vb->num_haplotypes_per_line; ht_i++) {
            char c = ht_matrix[ht_i * vb->lines.len + vb_line_i];
            
            // case: gtshark can't handle alleles>2 (natively, it splits them to several lines).
            // we put this allele in the exception list and change it to '0' for gtshark.
            if (c < '0' || c > '2') {
                if (!num_exceptions_in_line) { // first exception for this line 
                    buf_alloc_more (vb, x_line_buf, 1, 1, uint32_t, 2);
                    NEXTENT (uint32_t, *x_line_buf) = BGEN32 (vb_line_i);
                }

                buf_alloc_more (vb, x_ht_buf, 2, 2, uint16_t, 2); // room for terminator too
                NEXTENT (uint16_t, *x_ht_buf) = BGEN16 (ht_i - last_exception_ht_i); // delta encoding
                last_exception_ht_i = ht_i;
                
                buf_alloc_more (vb, x_allele_buf, 2, 2, char, 2);   // room for terminator too
                NEXTENT (char, *x_allele_buf) = c;

                num_exceptions_in_line++; 
                fprintf (file, "\t0");    
            }
            else 
                fprintf (file, "\t%c", c);
        }
        fprintf (file, "\n");

        if (num_exceptions_in_line) { // we have exceptions for this line - terminate the ht_i and allele arrays
            NEXTENT (uint16_t, *x_ht_buf) = 0;
            NEXTENT (char, *x_allele_buf) = 0;
        }
    }

    FCLOSE (file, gtshark_vcf_name);

    return;
}

// ZIP
static void codec_gtshark_run_compress (VBlock *vb, 
                                        const char *gtshark_vcf_name, const char *gtshark_db_name,
                                        const char *gtshark_db_db_name, const char *gtshark_db_gt_name)
{
    // remove in case of leftovers from previous run
    file_remove (gtshark_db_db_name, true);
    file_remove (gtshark_db_gt_name, true);

    codec_gtshark_run (vb->vblock_i, "compress-db", gtshark_vcf_name, gtshark_db_name);

    // read both gtshark output files
    file_get_file ((VBlockP)vb, gtshark_db_db_name, &vb->gtshark_db->local, "context->local", vb->vblock_i, false);
    file_get_file ((VBlockP)vb, gtshark_db_gt_name, &vb->gtshark_gt->local, "context->local", vb->vblock_i, false);
    
    file_remove (gtshark_vcf_name,   false);
    file_remove (gtshark_db_db_name, false);
    file_remove (gtshark_db_gt_name, false);
}

bool codec_gtshark_compress (VBlock *vb, 
                             SectionHeader *header,    
                             const char *uncompressed, // option 1 - compress contiguous data
                             uint32_t *uncompressed_len, 
                             LocalGetLineCB callback,  // option 2 - not supported
                             char *compressed, uint32_t *compressed_len /* in/out */, 
                             bool soft_fail)           // soft fail not supported
{
    char gtshark_vcf_name[strlen (z_file->name) + 20];
    sprintf (gtshark_vcf_name, "%s.%u.db.vcf", z_name, vb->vblock_i);

    char gtshark_db_name[strlen (z_file->name) + 20];
    sprintf (gtshark_db_name, "%s.%u.db", z_name, vb->vblock_i);

    char gtshark_db_db_name[strlen (z_file->name) + 20];
    sprintf (gtshark_db_db_name, "%s.%u.db_db", z_name, vb->vblock_i);

    char gtshark_db_gt_name[strlen (z_file->name) + 20];
    sprintf (gtshark_db_gt_name, "%s.%u.db_gt", z_name, vb->vblock_i);

    codec_gtshark_create_vcf_file (vb, gtshark_vcf_name);

    codec_gtshark_run_compress (vb, gtshark_vcf_name, gtshark_db_name, gtshark_db_db_name, gtshark_db_gt_name);

    vb->gtshark_db      ->ltype  = LT_UINT8;
    vb->gtshark_gt      ->ltype  = LT_UINT8;
    vb->gtshark_x_line  ->ltype  = LT_UINT32;
    vb->gtshark_x_ht    ->ltype  = LT_UINT16;
    vb->gtshark_x_allele->ltype  = LT_UINT8;

    vb->gtshark_db->lcodec = CODEC_NONE; // these are already compressed by gtshark and not further compressible
    vb->gtshark_gt->lcodec = CODEC_NONE;
    // note: codec of gtshark_x_* remains CODEC_UNKNOWN for automatic assignment
    
    // put a gtshark codec for uncompression on the last section to be written which would trigger codec_gtshark_uncompress
    // after all sections have been uncompressed with their main (simple) codec
    if (vb->gtshark_x_allele->local.len) // if we have exceptions - this is the last sections
        vb->gtshark_x_allele->lsubcodec_piz = CODEC_GTSHARK; 

    else // if we have no exceptions - this is the last section 
        vb->gtshark_gt->lsubcodec_piz = CODEC_GTSHARK;

    // note: we create the data in the 5 contexts gtshark_*, which will be compressed subsequently
    // with simple codecs, but for this ht_matrix context - no section should be created for it in the file
    buf_free (&vb->ht_matrix_ctx->local); 
    *compressed_len = 0;

    return true;
}

// PIZ
static char *codec_gtshark_write_db_file (uint32_t vb_i, const char *file_ext, const Buffer *buf)
{
    char *filename = malloc (strlen (z_name) + 50);
    sprintf (filename, "%s.%u.%s", z_name, vb_i, file_ext);

    FILE *file = fopen (filename, "wb");

    size_t bytes_written = fwrite (buf->data, 1, buf->len, file);
    ASSERT (bytes_written == buf->len, 
            "Error in uncompressing vb_i=%u: failed to write file %s - only %u out of %u bytes written: %s",
            vb_i, filename, (unsigned)bytes_written, (uint32_t)buf->len, strerror (errno));

    FCLOSE (file, filename);
    return filename;
}

// ----------
// PIZ side
// ----------
static void codec_gtshark_run_decompress (VBlock *vb)
{
    char gtshark_db_name[strlen (z_name) + 50];
    sprintf (gtshark_db_name, "%s.%u.db", z_name, vb->vblock_i);

    char gtshark_vcf_name[strlen (z_name) + 50];
    sprintf (gtshark_vcf_name, "%s.%u.vcf", z_name, vb->vblock_i);

    // remove in case of leftovers from previous run
    file_remove (gtshark_vcf_name, true);

    codec_gtshark_run (vb->vblock_i, "decompress-db", gtshark_db_name, gtshark_vcf_name);

    file_get_file ((VBlockP)vb, gtshark_vcf_name, &vb->compressed, "compressed", vb->vblock_i, true);

    file_remove (gtshark_vcf_name, false);
}    

static void codec_gtshark_piz_apply_exceptions_to_ht_matrix (VBlock *vb, uint8_t *ht_matrix, uint32_t num_lines, uint32_t num_hts)
{
    const uint32_t *exceptions_line_i_data = FIRSTENT (uint32_t, vb->gtshark_x_line->local);
    const uint16_t *next_ht_i_delta        = FIRSTENT (uint16_t, vb->gtshark_x_ht->local);
    const uint16_t *after_ht_i_delta       = AFTERENT (uint16_t, vb->gtshark_x_ht->local);
    const uint8_t  *next_allele            = FIRSTENT (uint8_t,  vb->gtshark_x_allele->local);
    const uint8_t  *after_allele           = AFTERENT (uint8_t,  vb->gtshark_x_allele->local);
    
    uint32_t num_x_lines = vb->gtshark_x_line->local.len;
    for (unsigned x_line_i=0; x_line_i < num_x_lines; x_line_i++) {

        uint32_t vb_line_i = exceptions_line_i_data[x_line_i];
        ASSERT (vb_line_i < num_lines, "Error in codec_gtshark_piz_apply_exceptions_to_ht_matrix: vb_i=%u: vb_line_i=%u is out of range (num_lines=%u)", 
                vb->vblock_i, vb_line_i, num_lines);

        uint16_t last_ht_i = 0;
        
        for (unsigned ex_ht_i=0 ; !ex_ht_i || *next_ht_i_delta; ex_ht_i++) { // the list of ht_i's for this line is terminated with a 0

            ASSERT (next_ht_i_delta < after_ht_i_delta, "Error in codec_gtshark_piz_apply_exceptions_to_ht_matrix: vb_i=%u: next_ht_i_delta is out of range", vb->vblock_i);
            ASSERT (next_allele < after_allele, "Error in codec_gtshark_piz_apply_exceptions_to_ht_matrix: vb_i=%u: next_allele is out of range", vb->vblock_i);

            uint16_t delta = *(next_ht_i_delta++);
            uint16_t ht_i = last_ht_i + delta; // decode delta encoding
            ASSERT (ht_i < num_hts, "Error in codec_gtshark_piz_apply_exceptions_to_ht_matrix: vb_i=%u: ht_i=%u is out of range (num_hts=%u)", 
                    vb->vblock_i, ht_i, num_hts);

            ASSERT (*next_allele <= '0' + VCF_MAX_ALLELE_VALUE, "Error in codec_gtshark_piz_apply_exceptions_to_ht_matrix: vb_i=%u: allele is out of range (ascii(%u))", vb->vblock_i, (unsigned)*next_allele);
            ht_matrix[ht_i * num_lines + vb_line_i] = *next_allele++; // transposed matrix
            last_ht_i = ht_i;
        }

        next_ht_i_delta++; // end-of-line separator in both ht_i and allele arrays
        next_allele++;
    }
}

// PIZ: convert the vcf generated by gtshark when decompressing the db, into our haplotype matrix
static void codec_gtshark_piz_reconstruct_ht_matrix (VBlock *vb)
{
    // get vb->num_haplotypes_per_line
    const char *substr = strstr (vb->compressed.data, GTSHARK_NUM_HT_PER_LINE);
    ASSERT (substr, "Error: cannot locate \"" GTSHARK_NUM_HT_PER_LINE "\" within gtshark-produced vcf data for vb_i=%u", vb->vblock_i);
    uint32_t num_hts   = vb->num_haplotypes_per_line = atoi (substr + strlen (GTSHARK_NUM_HT_PER_LINE));

    uint32_t num_lines = vb->lines.len;

    ASSERT (num_lines && num_hts, "Error in codec_gtshark_piz_reconstruct_ht_matrix: Expecting num_lines=%u and num_hts=%u to be >0", num_lines, num_hts);
    
    buf_alloc (vb, &vb->ht_matrix_ctx->local, num_lines * num_hts, 1, "context->local", 0);
    vb->ht_matrix_ctx->local.len = num_lines * num_hts;
    
    ARRAY (uint8_t, ht_matrix, vb->ht_matrix_ctx->local);
    
    // look for the start of our data - a newline followed by the chrom
    substr = strstr (vb->compressed.data, "\n" GTSHARK_CHROM_ID);
    ASSERT (substr, "Error in codec_gtshark_piz_reconstruct_ht_matrix: cannot locate start of data within gtshark-produced vcf data for vb_i=%u", vb->vblock_i);

    // build the transposed matrix from the vcf data - this will include alleles 0, 1 or 2 with higher
    // alleles showing as 0
    
    const char *next = substr + 1; // point to the CHROM = '1'
    const unsigned prefix_len = strlen (GTSHARK_VCF_LINE_VARDATA) + 1; // +1 for the tab after GT
    
    for (uint32_t vb_line_i=0; vb_line_i < num_lines; vb_line_i++) {
        ASSERT (*next == 'Z', "Error in codec_gtshark_piz_reconstruct_ht_matrix: expecting vb_line_i=%u to start with 'Z'", vb_line_i);
        next += prefix_len; // skipping the line fields 1-9, arriving at the first haplotype
    
        for (uint32_t ht_i=0; ht_i < num_hts; ht_i++) {
            ht_matrix[ht_i * num_lines + vb_line_i] = *next; // haplotype matrix is transposed
            next += 2; // skip past this ht and also the following \t or \n
        }
    }

    // now enter the higher alleles from the exception list
    if (vb->gtshark_x_line) // has exceptions
        codec_gtshark_piz_apply_exceptions_to_ht_matrix (vb, ht_matrix, num_lines, num_hts);

    buf_free (&vb->compressed);
}

// After all 5 contexts are uncompressed, this function is called, as CODEC_GTSHARK is set as the subcodec
// of the last section (gtshark_x_allele)
void codec_gtshark_uncompress (VBlock *vb, Codec codec,
                               const char *compressed, uint32_t compressed_len,
                               Buffer *uncompressed_buf, uint64_t uncompressed_len,
                               Codec sub_codec)
{
    vb->ht_matrix_ctx    = mtf_get_ctx (vb, dict_id_FORMAT_GT_HT); // create context, as it doesn't exist in the file
    vb->ht_matrix_ctx->lcodec = CODEC_GTSHARK;
    vb->ht_matrix_ctx->ltype  = LT_CODEC; // reconstruction will go to codec_gtshark_reconstruct as defined in codec_args for CODEC_GTSHARK

    vb->gtshark_db       = mtf_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_DB);
    vb->gtshark_gt       = mtf_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_GT);
    vb->gtshark_x_line   = mtf_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_X_LINE); // note: gtshark_x_* data exists only if there are exceptions
    vb->gtshark_x_ht     = mtf_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_X_HT);
    vb->gtshark_x_allele = mtf_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_X_ALLELE);

    char *filename_db_db = codec_gtshark_write_db_file (vb->vblock_i, "db_db", &vb->gtshark_db->local);
    char *filename_db_gt = codec_gtshark_write_db_file (vb->vblock_i, "db_gt", &vb->gtshark_gt->local);

    codec_gtshark_run_decompress (vb);                            

    file_remove (filename_db_db, false);
    file_remove (filename_db_gt, false);
    free (filename_db_db);
    free (filename_db_gt);

    codec_gtshark_piz_reconstruct_ht_matrix (vb); 
}

void codec_gtshark_reconstruct (VBlock *vb, Codec codec, Context *ctx)
{
    if (vb->dont_show_curr_line) return;

    // find next allele - skipping unused spots ('*')
    uint8_t ht = '*';
    do { 
        ht = NEXTLOCAL (uint8_t, vb->ht_matrix_ctx);
    } while (ht == '*' && vb->ht_matrix_ctx->next_local < vb->ht_matrix_ctx->local.len);

    if (ht == '.' || IS_DIGIT(ht)) 
        RECONSTRUCT1 (ht);
    
    else if (ht == '*') 
        ABORT ("Error in codec_gtshark_reconstruct: reconstructing txt_line=%u vb_i=%u: unexpected end of ctx->local data in %s (len=%u)", 
               vb->line_i, vb->vblock_i, ctx->name, (uint32_t)ctx->local.len)
    
    else { // allele 10 to 99 (ascii 58 to 147)
        RECONSTRUCT_INT (ht - '0');
    }
}
