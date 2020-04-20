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
#include "txtfile.h"
#include "optimize.h"
#include "random_access.h"

#define DATA_LINE(vb,i) (&((ZipDataLineSAM *)((vb)->data_lines))[(i)])

// get lengths of sequence data (SEQ+E2) and quality data (QUAL+U2)
static void zip_sam_get_seq_qual_len (VBlockSAM *vb, uint32_t *seq_len, uint32_t *qual_len)
{
    for (uint32_t vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {
        ZipDataLineSAM *dl = DATA_LINE (vb, vb_line_i);
        *seq_len  += dl->seq_len * (1 + !!dl->e2_data_start); // length SEQ and E2 - they must be the same per SAM file specification
        *qual_len += dl->seq_len * (1 + !!dl->u2_data_start); // length QUAL and U2 - they must be the same per SAM file specification        
    }
}

// callback function for compress to get data of one line (called by comp_lzma_data_in_callback)
static void zip_sam_get_start_len_line_i_seq (VBlock *vb, uint32_t vb_line_i, 
                                              char **line_seq_data, uint32_t *line_seq_len,  // out 
                                              char **line_e2_data,  uint32_t *line_e2_len)
{
    ZipDataLineSAM *dl = DATA_LINE (vb, vb_line_i);
    *line_seq_data = ENT (char, vb->txt_data, dl->seq_data_start);
    *line_seq_len  = dl->seq_len;
    *line_e2_data  = dl->e2_data_start ? ENT (char, vb->txt_data, dl->e2_data_start) : NULL;
    *line_e2_len   = dl->e2_data_start ? dl->seq_len : 0;
}   

// callback function for compress to get data of one line (called by comp_compress_bzlib)
static void zip_sam_get_start_len_line_i_qual (VBlock *vb, uint32_t vb_line_i, 
                                               char **line_qual_data, uint32_t *line_qual_len, // out
                                               char **line_u2_data,   uint32_t *line_u2_len) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb, vb_line_i);
     
    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);
    *line_qual_len  = dl->seq_len;
    *line_u2_data   = dl->u2_data_start ? ENT (char, vb->txt_data, dl->u2_data_start) : NULL;
    *line_u2_len    = dl->u2_data_start ? dl->seq_len : 0;

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    if (flag_optimize) {
        optimize_phred_quality_string (*line_qual_data, *line_qual_len);
        if (*line_u2_data) optimize_phred_quality_string (*line_u2_data, *line_u2_len);
    }
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
    seg_sam_initialize (vb);  

    seg_all_data_lines (vb_, seg_sam_data_line, sizeof (ZipDataLineSAM));

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

    // now, we merge vb->ra_buf into z_file->ra_buf
    random_access_merge_in_vb (vb_);

    // generate & write b250 data for all fields (FLAG to OPTIONAL)
    for (SamFields f=SAM_QNAME ; f <= SAM_OPTIONAL ; f++) {
        MtfContext *ctx = &vb->mtf_ctx[f];
        zip_generate_b250_section (vb_, ctx);
        zfile_compress_b250_data (vb_, ctx);
    }

    // generate & write b250 data for all QNAME subfields
    for (unsigned qname_sf_i=0; qname_sf_i < vb->qname_mapper.num_subfields; qname_sf_i++) {
        MtfContext *ctx = &vb->mtf_ctx[vb->qname_mapper.did_i[qname_sf_i]];
        zip_generate_b250_section (vb_, ctx);
        zfile_compress_b250_data  (vb_, ctx);
    }

    // generate & write b250 data for all OPTIONAL subfields
    for (unsigned did_i=datatype_last_field[DATA_TYPE_SAM]+1; did_i < vb->num_dict_ids; did_i++) {
                
        MtfContext *ctx = &vb->mtf_ctx[did_i];
        
        if (ctx->mtf_i.len && ctx->dict_section_type == SEC_SAM_OPTNL_SF_DICT) {
            zip_generate_b250_section (vb_, ctx);
            zfile_compress_b250_data (vb_, ctx);
        }
    }

    // generate & compress the MD and Random POS data
    zfile_compress_section_data_alg (vb_, SEC_SAM_RAND_POS_DATA,  &vb->random_pos_data, NULL, 0, COMPRESS_LZMA);
    zfile_compress_section_data_alg (vb_, SEC_SAM_MD_DATA, &vb->md_data, NULL, 0, COMPRESS_BZLIB);

    // generate & compress the SEQ & QUAL data
    uint32_t seq_len=0, qual_len=0;
    zip_sam_get_seq_qual_len (vb, &seq_len, &qual_len);
    zfile_compress_section_data_alg (vb_, SEC_SAM_SEQ_DATA,  NULL, zip_sam_get_start_len_line_i_seq,  seq_len,  COMPRESS_LZMA);
    zfile_compress_section_data_alg (vb_, SEC_SAM_QUAL_DATA, NULL, zip_sam_get_start_len_line_i_qual, qual_len, COMPRESS_BZLIB);

    // tell dispatcher this thread is done and can be joined. 
    // thread safety: this isn't protected by a mutex as it will just be false until it at some point turns to true
    // this this operation needn't be atomic, but it likely is anyway
    vb->is_processed = true; 

    COPY_TIMER (vb->profile.compute);
  }
