// ------------------------------------------------------------------
//   zip_gff3.c
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
#include "random_access.h"

#define DATA_LINE(i) ENT (ZipDataLineGFF3, vb->lines, i)

// this function receives all lines of a vblock and processes them
// in memory to the compressed format. This thread then terminates, and the I/O thread writes the output.
void zip_gff3_compress_one_vb (VBlockP vb_)
{ 
    START_TIMER;

    VBlockGFF3 *vb = (VBlockGFF3 *)vb_;

    // if we're vb_i=1 lock, and unlock only when we're done merging. all other vbs need
    // to wait for our merge. that is because our dictionaries are sorted
    if (vb->vblock_i == 1) mtf_vb_1_lock(vb_);

    // allocate memory for the final compressed data of this vb. allocate 20% of the 
    // vb size on the original file - this is normally enough. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, vb->vb_data_size / 5, 1.2, "z_data", 0);

    // clone global dictionaries while granted exclusive access to the global dictionaries
    mtf_clone_ctx (vb_);

    // split each line in this variant block to its components
    seg_all_data_lines (vb_, seg_gff3_data_line, seg_gff3_initialize, sizeof (ZipDataLineGFF3));

    // for the first vb only - sort dictionaries so that the most frequent entries get single digit
    // base-250 indices. This can be done only before any dictionary is written to disk, but likely
    // beneficial to all vbs as they are likely to more-or-less have the same frequent entries
    if (vb->vblock_i == 1) {
        mtf_sort_dictionaries_vb_1(vb_);

        // for a file that's compressed (and hence we don't know the size of content a-priori) AND we know the
        // compressed file size (eg a local gz/bz2 file or a compressed file on a ftp server) - we now estimate the 
        // txt_data_size_single that will be used for the global_hash and the progress indicator
        txtfile_estimate_txt_data_size (vb_);
    }

    zfile_compress_generic_vb_header (vb_); // vblock header

    // merge new words added in this vb into the z_file.mtf_ctx, ahead of zip_generate_b250_section() and
    // zip_vcf_generate_genotype_one_section(). writing indices based on the merged dictionaries. dictionaries are compressed. 
    // all this is done while holding exclusive access to the z_file dictionaries.
    mtf_merge_in_vb_ctx(vb_);

    // now, we merge vb->ra_buf into z_file->ra_buf
    random_access_merge_in_vb (vb_);

    // generate & write b250 data for all fields 
    zip_generate_and_compress_fields (vb_);

    // generate & write b250 data for all ATTRIBUTES subfields
    uint8_t num_info_subfields=0;
    for (unsigned did_i=0; did_i < MAX_DICTS; did_i++) {
                
        MtfContext *ctx = &vb->mtf_ctx[did_i];
        
        if (ctx->dict_section_type == SEC_GFF3_ATTRS_SF_DICT && ctx->mtf_i.len) {
            zip_generate_b250_section (vb_, ctx);
            zfile_compress_b250_data (vb_, ctx);
            num_info_subfields++;
        }
    }
    ASSERT (num_info_subfields <= MAX_SUBFIELDS, "Error: vb_i=%u has %u ATTRIBUTES subfields, which exceeds the maximum of %u",
            vb->vblock_i, num_info_subfields, MAX_SUBFIELDS);

    // generate & compress the ID and Genotype data
    zfile_compress_section_data_alg (vb_, SEC_NUMERIC_ID_DATA,  &vb->dbxref_numeric_data, NULL, 0, COMP_LZMA);
    
    // tell dispatcher this thread is done and can be joined. 
    // thread safety: this isn't protected by a mutex as it will just be false until it at some point turns to true
    // this this operation needn't be atomic, but it likely is anyway
    vb->is_processed = true; 

    COPY_TIMER (vb->profile.compute);
}
