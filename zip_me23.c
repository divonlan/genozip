// ------------------------------------------------------------------
//   zip_me2e.c
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

#define DATA_LINE(vb,i) (&((ZipDataLineME23 *)((vb)->data_lines))[(i)])

// this function receives all lines of a variant block and processes them
// in memory to the compressed format. This thread then terminates, and the I/O thread writes the output.
void zip_me23_compress_one_vb (VBlockP vb_)
{ 
    START_TIMER;

    VBlockME23 *vb = (VBlockME23 *)vb_;

    // if we're vb_i=1 lock, and unlock only when we're done merging. all other vbs need
    // to wait for our merge. that is because our dictionaries are sorted
    if (vb->vblock_i == 1) mtf_vb_1_lock(vb_);

    // allocate memory for the final compressed data of this vb. allocate 20% of the 
    // vb size on the original file - this is normally enough. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, vb->vb_data_size / 5, 1.2, "z_data", 0);

    // clone global dictionaries while granted exclusive access to the global dictionaries
    mtf_clone_ctx (vb_);

    // split each line in this variant block to its components
    seg_all_data_lines (vb_, seg_me23_data_line, seg_me23_initialize, 0);

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

    // generate & write b250 data for all fields 
    zip_generate_and_compress_fields (vb_);
    
    // generate & compress the ID and Genotype data
    zfile_compress_section_data_alg (vb_, SEC_ID_DATA,  &vb->rsid_data, NULL, 0, COMPRESS_LZMA);
    zfile_compress_section_data_alg (vb_, SEC_GT_DATA,  &vb->genotype_data, NULL, 0, COMPRESS_BZLIB);
    
    // tell dispatcher this thread is done and can be joined. 
    // thread safety: this isn't protected by a mutex as it will just be false until it at some point turns to true
    // this this operation needn't be atomic, but it likely is anyway
    vb->is_processed = true; 

    COPY_TIMER (vb->profile.compute);
  }
