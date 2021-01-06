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
#include <pthread.h>
#include "genozip.h"
#include "codec.h"
#include "buffer.h"
#include "file.h"
#include "endianness.h"
#include "stream.h"
#include "vblock.h"
#include "dict_id.h"
#include "strings.h"
#include "reconstruct.h"
#include "vcf_private.h"

// -------------
// ZIP & PIZ
// -------------

#define PIPE_MAX_BYTES 32768

// gtshark vcf file stuff
#define GTSHARK_CHROM_ID "Z"
#define GTSHARK_NUM_HT_PER_LINE "##num_hts="
#define GTSHARK_VCF_LINE_VARDATA GTSHARK_CHROM_ID "\t.\t.\t.\t.\t.\t.\t.\tGT\t"
static const unsigned vardata_len = 19;

// check stdout / stderr output of the gtshark process for errors
static void gtshark_check_pipe_for_errors (char *data, FILE *fp, uint32_t vb_i, bool is_stderr) 
{
    unsigned bytes = fread (data, 1, PIPE_MAX_BYTES-1, fp); // -1 to leave room for \0
    data[bytes] = '\0';

    ASSERTE (!strstr (data, "error")   && 
             !strstr (data, "E::")     && 
             !strstr (data, "Invalid") &&
             !strstr (data, "throw"),
             "gtshark failed for vb_i=%u. Here is its %s:\n%s\n", 
             vb_i, is_stderr ? "STDERR" : "STDOUT", data);
}

static bool codec_gtshark_run (uint32_t vb_i, const char *command, 
                               const char *filename_1, const char *filename_2) 
{
    StreamP gtshark = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 0, 0,
                                     "To use the --gtshark option",
                                     "gtshark", command, filename_1, filename_2, NULL);

    // read pipe (up to 10000 characters)
    char stdout_data[PIPE_MAX_BYTES], stderr_data[PIPE_MAX_BYTES];
    gtshark_check_pipe_for_errors (stdout_data, stream_from_stream_stdout (gtshark), vb_i, false);
    gtshark_check_pipe_for_errors (stderr_data, stream_from_stream_stderr (gtshark), vb_i, true);

    int exit_status = stream_close (&gtshark, STREAM_WAIT_FOR_PROCESS);  

#ifndef _WIN32

    ASSERTE (!WEXITSTATUS (exit_status), 
             "gtshark exited with status=%u for vb_i=%u. Here is its STDOUT:\n%s\nHere is the STDERR:\n%s\n", 
             WEXITSTATUS (exit_status), vb_i, stdout_data, stderr_data);

    ASSERTE (!WIFSIGNALED (exit_status),
             "gtshark process was killed by a signal, it was running for vb_i=%u. Here is its STDOUT:\n%s\nHere is the STDERR:\n%s\n", 
             vb_i, stdout_data, stderr_data);
#endif

    ASSERTE (!exit_status, "gtshark failed to exit normally for vb_i=%u. Here is its STDOUT:\n%s\nHere is the STDERR:\n%s\n", 
             vb_i, stdout_data, stderr_data);

    return true; // successfully executed gtshark
}

typedef struct { VBlockVCFP vb; const char *filename; Buffer *buf; } RWThreadArg;

static void *codec_gtshark_read_gtshark_output_file (void *arg)
{
    const char *filename = ((RWThreadArg *)arg)->filename;
    Buffer *buf          = ((RWThreadArg *)arg)->buf;
    VBlockVCFP vb        = ((RWThreadArg *)arg)->vb;

    FILE *file = fopen (filename, "rb");
    ASSERTE (file, "cannot open %s: %s (errno=%u)", filename, strerror (errno), errno);

    #define CHUNK 100000

    uint32_t bytes_read;
    do {
        buf_alloc (vb, buf, buf->len + CHUNK, 2, "contexts->local"); 
        bytes_read = fread (AFTERENT (char, *buf), 1, CHUNK, file);
        buf->len += bytes_read;
    } while (bytes_read == CHUNK); // its EOF if its smaller

    ASSERTE (!fclose (file), "fclose failed: %s", strerror (errno));

    file_remove (filename, false);

    return NULL;
}

// -------------
// ZIP side
// -------------

void codec_gtshark_comp_init (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    vb->ht_matrix_ctx  = ctx_get_ctx (vb, dict_id_FORMAT_GT_HT);
    vb->ht_matrix_ctx->lcodec = CODEC_GTSHARK; // this will trigger codec_gtshark_compress even though the section is not written to the file
    
    vb->gtshark_db_ctx = ctx_get_ctx (vb, dict_id_FORMAT_GT_SHARK_DB);
    vb->gtshark_gt_ctx = ctx_get_ctx (vb, dict_id_FORMAT_GT_SHARK_GT);
    vb->gtshark_ex_ctx = ctx_get_ctx (vb, dict_id_FORMAT_GT_SHARK_EX);

    // in --stats, consolidate stats into GT
    vb->gtshark_db_ctx->st_did_i = vb->gtshark_gt_ctx->st_did_i = vb->gtshark_ex_ctx->st_did_i =
    vb->ht_matrix_ctx->st_did_i = ctx_get_ctx (vb, dict_id_FORMAT_GT)->did_i;
}

typedef struct { VBlockVCFP vb; const char *fifo; } VcfThreadArg;

static void *codec_gtshark_zip_create_vcf_file (void *arg)
{
    VBlockVCFP vb = ((VcfThreadArg *)arg)->vb;
    const char *gtshark_vcf_name = ((VcfThreadArg *)arg)->fifo;

    const uint32_t num_hts = vb->num_haplotypes_per_line;

    FILE *file = fopen (gtshark_vcf_name, "wb");
    ASSERTE (file, "failed to create temporary file %s", gtshark_vcf_name);

    // write gtshark vcf header
    fprintf (file, "##fileformat=VCFv4.2\n" 
                   "##contig=<ID=" GTSHARK_CHROM_ID ">\n" 
                   "##FORMAT=<ID=GT>\n"
                   GTSHARK_NUM_HT_PER_LINE "%u\n"
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", num_hts);

    for (unsigned i=0; i < num_hts; i++) fprintf (file, "\t%x", i+1);
    fputc ('\n', file);

    // exceptions - one byte per matrix byte, ASCII 0 matrix is '0' or '1', or the matrix value if not
    buf_alloc (vb, &vb->gtshark_ex_ctx->local, vb->lines.len * num_hts, 1.1, "contexts->local");
    buf_zero (&vb->gtshark_ex_ctx->local);
    vb->gtshark_ex_ctx->local.len = vb->lines.len * num_hts;
    bool has_ex = false;

    ARRAY (char, ht_matrix, vb->ht_matrix_ctx->local);
    ARRAY (char, gtshark_ex, vb->gtshark_ex_ctx->local);

    // prepare line template
    buf_alloc (vb, &vb->compressed, vardata_len + num_hts*2, 1.2, "compressed");
    buf_add (&vb->compressed, GTSHARK_VCF_LINE_VARDATA, vardata_len);
    memset (AFTERENT (char, vb->compressed), '\t', num_hts*2);  
    vb->compressed.len += num_hts*2;
    *LASTENT (char, vb->compressed) = '\n';

    // add the haplotypes and write lines - one line at a time
    for (unsigned vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        for (unsigned ht_i=0; ht_i < num_hts; ht_i++) {
            char c = ht_matrix[vb_line_i * num_hts + ht_i]; 
            
            // case: gtshark can't handle alleles>2 (natively, it splits them to several lines).
            // we put this allele in the exception matrix and change it to '0' for gtshark.
            if (c < '0' || c > '2') {
                gtshark_ex[vb_line_i * num_hts + ht_i] = c;
                c = '0';
                has_ex = true;
            }

            *ENT (char, vb->compressed, vardata_len + ht_i*2) = c;
        }
        fwrite (vb->compressed.data, 1, vb->compressed.len, file);
    }

    buf_free (&vb->compressed);
    FCLOSE (file, gtshark_vcf_name);

    if (!has_ex) buf_free (&vb->gtshark_ex_ctx->local); // no exceptions - no need to write an exceptions section

    return NULL;
}

#define GET_FILENAMES_FIFOS(vb_i) \
    const char *tmpdir = getenv ("TMPDIR") ? getenv ("TMPDIR") : "/tmp"; \
    const char *base = file_basename (z_name, 0,0,0,0); \
    unsigned fn_len = strlen (z_file->name) + 20; \
    char gtshark_vcf_name[fn_len], gtshark_base_name[fn_len], gtshark_db_name[fn_len], gtshark_gt_name[fn_len]; \
    sprintf (gtshark_vcf_name,  "%s/%s.%u.vcf", tmpdir, base, vb->vblock_i); \
    sprintf (gtshark_base_name, "%s/%s.%u", tmpdir, base, vb->vblock_i); \
    sprintf (gtshark_db_name,   "%s_db", gtshark_base_name); \
    sprintf (gtshark_gt_name,   "%s_gt", gtshark_base_name); \
    FREE (base); \
    file_mkfifo (gtshark_vcf_name); \
    file_mkfifo (gtshark_db_name); \
    file_mkfifo (gtshark_gt_name); 

bool codec_gtshark_compress (VBlock *vb_, 
                             SectionHeader *header,    
                             const char *uncompressed, // option 1 - compress contiguous data
                             uint32_t *uncompressed_len, 
                             LocalGetLineCB callback,  // option 2 - not supported
                             char *compressed, uint32_t *compressed_len /* in/out */, 
                             bool soft_fail)           // soft fail not supported
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    pthread_t vcf_thread, db_thread, gt_thread;
    GET_FILENAMES_FIFOS (vb->vblock_i);

    // create the VCF file to be consumed by gtshark - using a separate thread that will write to a FIFO (named pipe)
    // while gtshark reads from it
    VcfThreadArg vcf_thread_arg = { vb, gtshark_vcf_name };
    int err = pthread_create (&vcf_thread, NULL, codec_gtshark_zip_create_vcf_file, &vcf_thread_arg);
    ASSERTE (!err, "failed to create thread vcf_thread: %s", strerror (err));

    // Create reader threads to consume both gtshark output files, and remove them when done
    RWThreadArg db_read_thread_arg = { vb, gtshark_db_name, &vb->gtshark_db_ctx->local };
    err = pthread_create (&db_thread, NULL, codec_gtshark_read_gtshark_output_file, &db_read_thread_arg);
    ASSERTE (!err, "failed to create thread db_thread: %s", strerror (err));

    RWThreadArg gt_read_thread_arg = { vb, gtshark_gt_name, &vb->gtshark_gt_ctx->local };
    err = pthread_create (&gt_thread, NULL, codec_gtshark_read_gtshark_output_file, &gt_read_thread_arg);
    ASSERTE (!err, "failed to create thread gt_thread: %s", strerror (err));

    // Run gtshark - reading and writing data to fifos - i.e. in memory rather than disk
    codec_gtshark_run (vb->vblock_i, "compress-db", gtshark_vcf_name, gtshark_base_name);
    
    file_remove (gtshark_vcf_name, false);

    // wait for threads to complete (writer first, then readers)
    pthread_join (vcf_thread, NULL);
    pthread_join (db_thread,  NULL);
    pthread_join (gt_thread,  NULL);

    // set parameters for output sections
    vb->gtshark_db_ctx->ltype  = LT_UINT8;
    vb->gtshark_gt_ctx->ltype  = LT_UINT8;
    vb->gtshark_ex_ctx->ltype  = LT_UINT8;

    vb->gtshark_db_ctx->lcodec = CODEC_NONE; // these are already compressed by gtshark and not further compressible
    vb->gtshark_gt_ctx->lcodec = CODEC_NONE;

    // put a gtshark codec for uncompression on the last section to be written which would trigger codec_gtshark_uncompress
    // after all sections have been uncompressed with their main (simple) codec
    if (vb->gtshark_ex_ctx->local.len) {
        // since codecs were already assigned to contexts before compression of all contexts begun, but
        // we just created this context now, we assign a codec manually
        codec_assign_best_codec ((VBlockP)vb, vb->gtshark_ex_ctx, NULL, SEC_LOCAL);

        vb->gtshark_ex_ctx->lsubcodec_piz = CODEC_GTSHARK; // we have exceptions - EX is the last section
    }
    else 
        vb->gtshark_gt_ctx->lsubcodec_piz = CODEC_GTSHARK; // no exceptions - GT is the last section

    // note: we created the data in the 3 contexts gtshark_*, which will be compressed subsequently.
    // For this ht_matrix context - no section should be created for it in the file
    buf_free (&vb->ht_matrix_ctx->local); 
    *compressed_len = 0;

    return true;
}

// ----------
// PIZ side
// ----------

static void *codec_gtshark_write_gtshark_input_file (void *arg)
{
    const char *filename = ((RWThreadArg *)arg)->filename;
    Buffer *buf          = ((RWThreadArg *)arg)->buf;
    VBlockVCFP vb        = ((RWThreadArg *)arg)->vb;

    FILE *file = fopen (filename, "wb");

    size_t bytes_written = fwrite (buf->data, 1, buf->len, file);
    ASSERTE (bytes_written == buf->len, 
             "while uncompressing vb_i=%u: failed to write file %s - only %u out of %u bytes written: %s",
             vb->vblock_i, filename, (unsigned)bytes_written, (uint32_t)buf->len, strerror (errno));

    FCLOSE (file, filename);
    return NULL;
}

// PIZ: convert the vcf generated by gtshark when decompressing the db, into our haplotype matrix
static void codec_gtshark_reconstruct_ht_matrix (VBlockVCF *vb)
{
    // get vb->num_haplotypes_per_line
    const char *substr = strstr (vb->compressed.data, GTSHARK_NUM_HT_PER_LINE);
    ASSERTE (substr, "cannot locate \"" GTSHARK_NUM_HT_PER_LINE "\" within gtshark-produced vcf data for vb_i=%u", vb->vblock_i);
    uint32_t num_hts = vb->num_haplotypes_per_line = atoi (substr + strlen (GTSHARK_NUM_HT_PER_LINE));

    uint32_t num_lines = vb->lines.len;

    ASSERTE (num_lines && num_hts, "Expecting num_lines=%u and num_hts=%u to be >0", num_lines, num_hts);
    
    buf_alloc (vb, &vb->ht_matrix_ctx->local, num_lines * num_hts, 1, "contexts->local");
    vb->ht_matrix_ctx->local.len = num_lines * num_hts;
    
    ARRAY (uint8_t, ht_matrix, vb->ht_matrix_ctx->local);
    
    // look for the start of our data - a newline followed by the chrom
    substr = strstr (vb->compressed.data, "\n" GTSHARK_CHROM_ID);
    ASSERTE (substr, "cannot locate start of data within gtshark-produced vcf data for vb_i=%u", vb->vblock_i);

    // build the transposed matrix from the vcf data - this will include alleles 0, 1 or 2 with higher
    // alleles showing as 0
    const char *next = substr + 1; // point to the CHROM = '1'
    
    for (uint32_t vb_line_i=0; vb_line_i < num_lines; vb_line_i++) {
        ASSERTE (*next == GTSHARK_CHROM_ID[0], "expecting vb_line_i=%u to start with '" GTSHARK_CHROM_ID "' but seeing '%c' (ASCII %u)", 
                 vb_line_i, *next, *next);

        next += vardata_len; // skipping the line fields 1-9, arriving at the first haplotype
    
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
void codec_gtshark_uncompress (VBlock *vb_, Codec codec,
                               const char *compressed, uint32_t compressed_len,
                               Buffer *uncompressed_buf, uint64_t uncompressed_len,
                               Codec sub_codec)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;
    pthread_t vcf_thread, db_thread, gt_thread;
    GET_FILENAMES_FIFOS (vb->vblock_i);
    
    vb->ht_matrix_ctx    = ctx_get_ctx (vb, dict_id_FORMAT_GT_HT); // create context, as it doesn't exist in the file
    vb->ht_matrix_ctx->lcodec = CODEC_GTSHARK;
    vb->ht_matrix_ctx->ltype  = LT_CODEC; // reconstruction will go to codec_gtshark_reconstruct as defined in codec_args for CODEC_GTSHARK

    vb->gtshark_db_ctx = ctx_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_DB);
    vb->gtshark_gt_ctx = ctx_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_GT);
    vb->gtshark_ex_ctx = ctx_get_existing_ctx (vb, dict_id_FORMAT_GT_SHARK_EX); // note: gtshark_ex_ctx will be set only if there was a gtshark_ex section

    ASSERTE0 (vb->gtshark_db_ctx, "gtshark_db_ctx is NULL, perhaps the section is missing in the file?");
    ASSERTE0 (vb->gtshark_gt_ctx, "gtshark_gt_ctx is NULL, perhaps the section is missing in the file?");

    // create threads for writing gtshark files for its decompressor to consume 
    RWThreadArg db_write_thread_arg = { vb, gtshark_db_name, &vb->gtshark_db_ctx->local };
    int err = pthread_create (&db_thread, NULL, codec_gtshark_write_gtshark_input_file, &db_write_thread_arg);
    ASSERTE (!err, "failed to create thread db_thread: %s", strerror (err));

    RWThreadArg gt_write_thread_arg = { vb, gtshark_gt_name, &vb->gtshark_gt_ctx->local };
    err = pthread_create (&gt_thread, NULL, codec_gtshark_write_gtshark_input_file, &gt_write_thread_arg);
    ASSERTE (!err, "failed to create thread gt_thread: %s", strerror (err));

    // thread from reading VCF file outputed by gtshark decompressor
    RWThreadArg vcf_read_thread_arg = { vb, gtshark_vcf_name, &vb->compressed };
    err = pthread_create (&vcf_thread, NULL, codec_gtshark_read_gtshark_output_file, &vcf_read_thread_arg);
    ASSERTE (!err, "failed to create thread vcf_thread: %s", strerror (err));

    codec_gtshark_run (vb->vblock_i, "decompress-db", gtshark_base_name, gtshark_vcf_name);

    // wait for threads to complete (writers first, then reader)
    pthread_join (db_thread,  NULL);
    pthread_join (gt_thread,  NULL);
    pthread_join (vcf_thread, NULL);
    
    file_remove (gtshark_db_name, false);
    file_remove (gtshark_gt_name, false);

    codec_gtshark_reconstruct_ht_matrix (vb); 
}

void codec_gtshark_reconstruct (VBlock *vb_, Codec codec, Context *ctx)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    if (vb->dont_show_curr_line) return;

    // find next allele - skipping unused spots ('*')
    uint8_t ht = '*';
    do { 
        ht = NEXTLOCAL (uint8_t, vb->ht_matrix_ctx);
    } while (ht == '*' && vb->ht_matrix_ctx->next_local < vb->ht_matrix_ctx->local.len);

    if (ht == '.' || IS_DIGIT(ht)) 
        RECONSTRUCT1 (ht);
    
    else if (ht == '\b') 
        vb->txt_data.len--; // ploidy padding (starting 10.0.2) - appears as the 2nd+ HT - counted in GT.repeats: don't reconstruct anything, just remove previous phase character
    
    else if (ht == '*') 
        ABORT ("Error in codec_gtshark_reconstruct: reconstructing txt_line=%u vb_i=%u: unexpected end of ctx->local data in %s (len=%u)", 
               vb->line_i, vb->vblock_i, ctx->name, (uint32_t)ctx->local.len)
    
    else { // allele 10 to 99 (ascii 58 to 147)
        RECONSTRUCT_INT (ht - '0');
    }
}
