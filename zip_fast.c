// ------------------------------------------------------------------
//   zip_fast.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "vblock.h"
#include "buffer.h"
#include "zfile.h"
#include "header.h"
#include "zip.h"
#include "seg.h"
#include "txtfile.h"
#include "optimize.h"

#define DATA_LINE(vb,i) (&((ZipDataLineFAST *)((vb)->data_lines))[(i)])

// get lengths of sequence data (SEQ+E2) and quality data (QUAL+U2)
static void zip_fast_get_seq_len (VBlockFAST *vb, uint32_t *seq_len)
{
    *seq_len = 0;

    for (uint32_t vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) 
        *seq_len  += DATA_LINE (vb, vb_line_i)->seq_len;
}

// callback function for compress to get data of one line (called by comp_lzma_data_in_callback)
static void zip_fast_get_start_len_line_i_seq (VBlock *vb, uint32_t vb_line_i, 
                                               char **line_seq_data, uint32_t *line_seq_len,  // out 
                                               char **unused_data,  uint32_t *unused_len)
{
    ZipDataLineFAST *dl = DATA_LINE (vb, vb_line_i);
    *line_seq_data = ENT (char, vb->txt_data, dl->seq_data_start);
    *line_seq_len  = dl->seq_len;
    *unused_data   = NULL;
    *unused_len    = 0;
}   

// callback function for compress to get data of one line (called by comp_compress_bzlib)
static void zip_fast_get_start_len_line_i_qual (VBlock *vb, uint32_t vb_line_i, 
                                                char **line_qual_data, uint32_t *line_qual_len, // out
                                                char **unused_data,   uint32_t *unused_len) 
{
    ZipDataLineFAST *dl = DATA_LINE (vb, vb_line_i);
     
    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);
    *line_qual_len  = dl->seq_len;
    *unused_data    = NULL;
    *unused_len     = 0;

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    if (flag_optimize) optimize_phred_quality_string (*line_qual_data, *line_qual_len);
}

// this function receives all lines of a variant block and processes them
// in memory to the compressed format. This thread then terminates, and the I/O thread writes the output.
void zip_fastq_compress_one_vb (VBlockP vb_)
{ 
    START_TIMER;

    VBlockFAST *vb = (VBlockFAST *)vb_;
    bool is_fastq = (vb->data_type == DT_FASTQ);

    // if we're vb_i=1 lock, and unlock only when we're done merging. all other vbs need
    // to wait for our merge. that is because our dictionaries are sorted
    if (vb->vblock_i == 1) mtf_vb_1_lock(vb_);

    // allocate memory for the final compressed data of this vb. allocate 20% of the 
    // vb size on the original file - this is normally enough. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, vb->vb_data_size / 3, 1.2, "z_data", 0);

    // clone global dictionaries while granted exclusive access to the global dictionaries
    mtf_clone_ctx (vb_);

    // split each line in this variant block to its components
    seg_all_data_lines (vb_, is_fastq ? seg_fastq_data_line : seg_fasta_data_line, false, 
                        is_fastq ? sizeof (ZipDataLineFAST) : 0);

    // for the first vb only - sort dictionaries so that the most frequent entries get single digit
    // base-250 indices. This can be done only before any dictionary is written to disk, but likely
    // beneficial to all vbs as they are likely to more-or-less have the same frequent entries
    if (vb->vblock_i == 1) {
        mtf_sort_dictionaries_vb_1(vb_);

        // for a FASTQ/FASTA file that's compressed (and hence we don't know the size of content a-priori) AND we know the
        // compressed file size (eg a local gz/bz2 file or a compressed file on a ftp server) - we now estimate the 
        // txt_data_size_single that will be used for the global_hash and the progress indicator
        txtfile_estimate_txt_data_size (vb_);
    }

    zfile_compress_generic_vb_header (vb_); // vblock header

    // merge new words added in this vb into the z_file.mtf_ctx, ahead of zip_generate_b250_section() and
    // zip_vcf_generate_genotype_one_section(). writing indices based on the merged dictionaries. dictionaries are compressed. 
    // all this is done while holding exclusive access to the z_file dictionaries.
    mtf_merge_in_vb_ctx(vb_);

    // All fields 
    zip_generate_and_compress_fields (vb_);

    // DESC subfields
    zip_generate_and_compress_subfields (vb_, &vb->desc_mapper);

    // SEQ
    uint32_t seq_len;
    zip_fast_get_seq_len (vb, &seq_len);
    zfile_compress_section_data_alg (vb_, SEC_SEQ_DATA,  NULL, zip_fast_get_start_len_line_i_seq,  seq_len, COMPRESS_LZMA);

    // FASTQ only: QUAL
    if (is_fastq) 
        zfile_compress_section_data_alg (vb_, SEC_QUAL_DATA, NULL, zip_fast_get_start_len_line_i_qual, seq_len, COMPRESS_BZLIB);

    // FASTA obly: COMMENT
    else 
        zfile_compress_section_data (vb_, SEC_FASTA_COMMENT_DATA, &vb->comment_data);

    // tell dispatcher this thread is done and can be joined. 
    // thread safety: this isn't protected by a mutex as it will just be false until it at some point turns to true
    // this this operation needn't be atomic, but it likely is anyway
    vb->is_processed = true; 

    COPY_TIMER (vb->profile.compute);
  }
