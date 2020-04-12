// ------------------------------------------------------------------
//   zip_sam.c
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

#define DATA_LINE(vb,i) (&((ZipDataLineSAM *)((vb)->data_lines))[(i)])

// get total number of bases in this VB
static uint32_t zip_sam_get_num_bases (VBlockSAM *vb)
{
    // calculate length
    uint32_t num_bases=0;
    for (uint32_t vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {
        ZipDataLineSAM *dl = DATA_LINE (vb, vb_line_i);
        num_bases += dl->seq_len; // length of SEQ and QUAL - they must be the same per SAM file specification
    }

    return num_bases;
}

static void zip_sam_get_start_len_line_i_seq (VBlock *vb, uint32_t vb_line_i, 
                                              char **line_data, uint32_t *line_data_len) // out
{
    ZipDataLineSAM *dl = DATA_LINE (vb, vb_line_i);
    uint32_t start_index = dl->seq_data_start;
    *line_data = ENT (char, vb->txt_data, start_index);
    *line_data_len   = dl->seq_len;
}   

static void zip_sam_get_start_len_line_i_qual (VBlock *vb, uint32_t vb_line_i, 
                                               char **line_data, uint32_t *line_data_len) // out
{
    ZipDataLineSAM *dl = DATA_LINE (vb, vb_line_i);
    uint32_t start_index = dl->qual_data_start;
    *line_data = ENT (char, vb->txt_data, start_index);
    *line_data_len   = dl->seq_len;
}   

// this function receives all lines of a variant block and processes them
// in memory to the compressed format. This thread then terminates the I/O thread writes the output.
void zip_sam_compress_one_vb (VBlockP vb_)
{ 
    START_TIMER;

    VBlockSAM *vb = (VBlockSAM *)vb_;

    // if we're vb_i=1 lock, and unlock only when we're done merging. all other vbs need
    // to wait for our merge. that is because our dictionaries are sorted
    if (vb->vblock_i == 1) mtf_vb_1_lock(vb_);

    // allocate memory for the final compressed data of this vb. allocate 20% of the 
    // vb size on the original file - this is normally enough. if not, we will realloc downstream
    buf_alloc (vb, &vb->z_data, vb->vb_data_size / 5, 1.2, "z_data", 0);

    // clone global dictionaries while granted exclusive access to the global dictionaries
    mtf_clone_ctx (vb_);

    // split each line in this variant block to its components
    seg_all_data_lines (vb_, seg_sam_data_line, sizeof (ZipDataLineSAM),
                        SAM_QNAME, SAM_OPTIONAL, sam_field_names, SEC_SAM_QNAME_DICT);

    // for the first vb only - sort dictionaries so that the most frequent entries get single digit
    // base-250 indices. This can be done only before any dictionary is written to disk, but likely
    // beneficial to all vbs as they are likely to more-or-less have the same frequent entries
    if (vb->vblock_i == 1) {
        mtf_sort_dictionaries_vb_1(vb_);

        // for a SAM file that's compressed (and hence we don't know the size of SAM content a-priori) AND we know the
        // compressed file size (eg a local gz/bz2 file or a compressed file on a ftp server) - we now estimate the 
        // txt_data_size_single that will be used for the global_hash and the progress indicator
        txtfile_estimate_txt_data_size (vb_);
    }

    //unsigned variant_data_header_pos = vb->z_data.len;
    zfile_sam_compress_vb_header (vb); // vblock header

    // merge new words added in this vb into the z_file.mtf_ctx, ahead of zip_generate_b250_section() and
    // zip_vcf_generate_genotype_one_section(). writing indices based on the merged dictionaries. dictionaries are compressed. 
    // all this is done while holding exclusive access to the z_file dictionaries.
    mtf_merge_in_vb_ctx(vb_);

    // generate & write b250 data for all fields (FLAG to OPTIONAL)
    for (SamFields f=SAM_QNAME ; f <= SAM_OPTIONAL ; f++) {
        MtfContext *ctx = &vb->mtf_ctx[f];
        if (ctx->mtf_i.len) {
            zip_generate_b250_section (vb_, ctx);
            zfile_compress_b250_data (vb_, ctx);
        }
    }

    // generate & write b250 data for all QNAME and OPTIONAL subfields
    for (unsigned did_i=0; did_i < MAX_DICTS; did_i++) {
                
        MtfContext *ctx = &vb->mtf_ctx[did_i];
        
        if (ctx->mtf_i.len && 
            (ctx->dict_section_type == SEC_SAM_QNAME_SF_DICT || ctx->dict_section_type == SEC_SAM_OPTNL_SF_DICT)) {
            zip_generate_b250_section (vb_, ctx);
            zfile_compress_b250_data (vb_, ctx);
        }
    }

    // generate & compress the SEQ & QUAL data
    uint32_t num_bases = zip_sam_get_num_bases (vb); 
    //zfile_compress_section_data_alg (vb_, SEC_SAM_SEQ_DATA,  NULL, zip_sam_get_start_len_line_i_seq,  num_bases, COMPRESS_LZMA);
    zfile_compress_section_data_alg (vb_, SEC_SAM_SEQ_DATA,  NULL, zip_sam_get_start_len_line_i_seq,  num_bases, COMPRESS_BZLIB);
    zfile_compress_section_data_alg (vb_, SEC_SAM_QUAL_DATA, NULL, zip_sam_get_start_len_line_i_qual, num_bases, COMPRESS_BZLIB);

    // tell dispatcher this thread is done and can be joined. 
    // thread safety: this isn't protected by a mutex as it will just be false until it at some point turns to true
    // this this operation needn't be atomic, but it likely is anyway
    vb->is_processed = true; 

    COPY_TIMER (vb->profile.compute);
  }
