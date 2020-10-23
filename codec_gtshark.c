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
    
    vb->gtshark_db_ctx = mtf_get_ctx (vb, dict_id_FORMAT_GT_SHARK_DB);
    vb->gtshark_gt_ctx = mtf_get_ctx (vb, dict_id_FORMAT_GT_SHARK_GT);
    vb->gtshark_ex_ctx = mtf_get_ctx (vb, dict_id_FORMAT_GT_SHARK_EX);
}

static void codec_gtshark_create_vcf_file (VBlock *vb, const char *gtshark_vcf_name)
{
    const uint32_t num_hts = vb->num_haplotypes_per_line;

    FILE *file = fopen (gtshark_vcf_name, "wb");
    ASSERT (file, "Error: failed to create temporary file %s", gtshark_vcf_name);

    fputs ("##fileformat=VCFv4.2\n", file);
    fputs ("##contig=<ID=Z>\n", file);
    fputs ("##FORMAT=<ID=GT>\n", file);

    #define GTSHARK_NUM_HT_PER_LINE "##num_haplotypes_per_line="
    fprintf (file, GTSHARK_NUM_HT_PER_LINE "%u\n", num_hts);
    fputs ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", file);

    for (unsigned i=0; i < num_hts; i++)
        fprintf (file, "\t%u", i+1);
    fputc ('\n', file);

    // exceptions - one byte per matrix byte, ASCII 0 matrix is '0' or '1', or the matrix value if not
    buf_alloc (vb, &vb->gtshark_ex_ctx->local, vb->lines.len * num_hts, 1.1, "context->local", 0);
    buf_zero (&vb->gtshark_ex_ctx->local);
    vb->gtshark_ex_ctx->local.len = vb->lines.len * num_hts;
    bool has_ex = false;

    ARRAY (char, ht_matrix, vb->ht_matrix_ctx->local);
    ARRAY (char, gtshark_ex, vb->gtshark_ex_ctx->local);

    for (unsigned vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        #define GTSHARK_CHROM_ID "Z"
        #define GTSHARK_VCF_LINE_VARDATA GTSHARK_CHROM_ID "\t.\t.\t.\t.\t.\t.\t.\tGT"

        fprintf (file, GTSHARK_VCF_LINE_VARDATA);

        for (unsigned ht_i=0; ht_i < num_hts; ht_i++) {
            char s[2] = { '\t', ht_matrix[vb_line_i * num_hts + ht_i] }; 
            
            // case: gtshark can't handle alleles>2 (natively, it splits them to several lines).
            // we put this allele in the exception matrix and change it to '0' for gtshark.
            if (s[1] < '0' || s[1] > '2') {
                gtshark_ex[vb_line_i * num_hts + ht_i] = s[1];
                s[1] = '0';
                has_ex = true;
            }

            fwrite (s, 1, 2, file);
        }
        fputc ('\n', file);
    }

    FCLOSE (file, gtshark_vcf_name);

    if (!has_ex) buf_free (&vb->gtshark_ex_ctx->local); // no exceptions - no need to write an exceptions section
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
    file_get_file ((VBlockP)vb, gtshark_db_db_name, &vb->gtshark_db_ctx->local, "context->local", vb->vblock_i, false);
    file_get_file ((VBlockP)vb, gtshark_db_gt_name, &vb->gtshark_gt_ctx->local, "context->local", vb->vblock_i, false);
    
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

    vb->gtshark_db_ctx->ltype  = LT_UINT8;
    vb->gtshark_gt_ctx->ltype  = LT_UINT8;
    vb->gtshark_ex_ctx->ltype  = LT_UINT8;

    vb->gtshark_db_ctx->lcodec = CODEC_NONE; // these are already compressed by gtshark and not further compressible
    vb->gtshark_gt_ctx->lcodec = CODEC_NONE;
    vb->gtshark_ex_ctx->lcodec = CODEC_BZ2;  // a sparse matrix 
    
    // put a gtshark codec for uncompression on the last section to be written which would trigger codec_gtshark_uncompress
    // after all sections have been uncompressed with their main (simple) codec
    if (vb->gtshark_ex_ctx->local.len) // if we have exceptions - this is the last sections
        vb->gtshark_ex_ctx->lsubcodec_piz = CODEC_GTSHARK; 
    else 
        vb->gtshark_gt_ctx->lsubcodec_piz = CODEC_GTSHARK;

    // note: we created the data in the 3 contexts gtshark_*, which will be compressed subsequently.
    // For this ht_matrix context - no section should be created for it in the file
    buf_free (&vb->ht_matrix_ctx->local); 
    *compressed_len = 0;

    return true;
}

// ----------
// PIZ side
// ----------

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

// PIZ: convert the vcf generated by gtshark when decompressing the db, into our haplotype matrix
static void codec_gtshark_piz_reconstruct_ht_matrix (VBlock *vb)
{
    // get vb->num_haplotypes_per_line
    const char *substr = strstr (vb->compressed.data, GTSHARK_NUM_HT_PER_LINE);
    ASSERT (substr, "Error: cannot locate \"" GTSHARK_NUM_HT_PER_LINE "\" within gtshark-produced vcf data for vb_i=%u", vb->vblock_i);
    uint32_t num_hts = vb->num_haplotypes_per_line = atoi (substr + strlen (GTSHARK_NUM_HT_PER_LINE));

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
            ht_matrix[vb_line_i * num_hts + ht_i] = *next; // haplotypes it gtshark vcf file are transposed vs ht_matrix
            next += 2; // skip past this ht and also the following \t or \n
        }
    }

    // case: we have exceptions - update the ht_matrix with the "exceptional" alleles
    if (vb->gtshark_ex_ctx) { // has exceptions

        ARRAY (uint8_t, gtshark_ex, vb->gtshark_ex_ctx->local);

        for (uint64_t i=0; i < num_lines * num_hts; i++)
            if (gtshark_ex[i]) ht_matrix[i] = gtshark_ex[i];
    }

    buf_free (&vb->compressed); // free the gtshark-generated reconstructed vcf file data
}

// After all 5 contexts are uncompressed, this function is called, as CODEC_GTSHARK is set as the subcodec
// of the last section (gtshark_ex_allele)
void codec_gtshark_uncompress (VBlock *vb, Codec codec,
                               const char *compressed, uint32_t compressed_len,
                               Buffer *uncompressed_buf, uint64_t uncompressed_len,
                               Codec sub_codec)
{
    vb->ht_matrix_ctx    = mtf_get_ctx (vb, dict_id_FORMAT_GT_HT); // create context, as it doesn't exist in the file
    vb->ht_matrix_ctx->lcodec = CODEC_GTSHARK;
    vb->ht_matrix_ctx->ltype  = LT_CODEC; // reconstruction will go to codec_gtshark_reconstruct as defined in codec_args for CODEC_GTSHARK

    vb->gtshark_db_ctx = mtf_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_DB);
    vb->gtshark_gt_ctx = mtf_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_GT);
    vb->gtshark_ex_ctx = mtf_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_EX); // note: gtshark_ex_ctx will be set only if there was a gtshark_ex section

    char *filename_db_db = codec_gtshark_write_db_file (vb->vblock_i, "db_db", &vb->gtshark_db_ctx->local);
    char *filename_db_gt = codec_gtshark_write_db_file (vb->vblock_i, "db_gt", &vb->gtshark_gt_ctx->local);

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
