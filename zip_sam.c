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

#define DATA_LINE(i) ENT (ZipDataLineSAM, vb->lines, i)

// get lengths of sequence data (SEQ+E2) and quality data (QUAL+U2)
static void zip_sam_get_data_lengths (VBlockSAM *vb, 
                                      uint32_t *seq_len, uint32_t *qual_len, 
                                      uint32_t *bd_len, uint32_t *bi_len)
{
    *seq_len = *qual_len = *bi_len = *bd_len = 0;

    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {
        ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
        *seq_len  += dl->seq_len * (1 + !!dl->e2_data_start); // length SEQ and E2 - they must be the same per SAM file specification
        *qual_len += dl->seq_len * (1 + !!dl->u2_data_start); // length QUAL and U2 - they must be the same per SAM file specification        
        *bd_len   += dl->bd_data_len;
        *bi_len   += dl->bi_data_len;
    }
}

// callback function for compress to get data of one line (called by comp_lzma_data_in_callback)
static void zip_sam_get_start_len_line_i_seq (VBlock *vb, uint32_t vb_line_i, 
                                              char **line_seq_data, uint32_t *line_seq_len,  // out 
                                              char **line_e2_data,  uint32_t *line_e2_len)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    *line_seq_data = ENT (char, vb->txt_data, dl->seq_data_start);
    *line_seq_len  = dl->seq_data_len;
    *line_e2_data  = dl->e2_data_start ? ENT (char, vb->txt_data, dl->e2_data_start) : NULL;
    *line_e2_len   = dl->e2_data_len;

    // if SEQ is just "*" (i.e. unavailable) replace it by " " for consistency with QUAL (= simplify PIZ code)
    if (dl->seq_data_len == 1 && (*line_seq_data)[0] == '*') 
        *line_seq_data = " "; // pointer to static string
}   

// callback function for compress to get data of one line (called by comp_compress_bzlib)
static void zip_sam_get_start_len_line_i_qual (VBlock *vb, uint32_t vb_line_i, 
                                               char **line_qual_data, uint32_t *line_qual_len, // out
                                               char **line_u2_data,   uint32_t *line_u2_len) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
     
    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);
    *line_qual_len  = dl->qual_data_len;
    *line_u2_data   = dl->u2_data_start ? ENT (char, vb->txt_data, dl->u2_data_start) : NULL;
    *line_u2_len    = dl->u2_data_len;

    // if QUAL is just "*" (i.e. unavailable) replace it by " " because '*' is a legal PHRED quality value that will confuse PIZ
    if (dl->qual_data_len == 1 && (*line_qual_data)[0] == '*') 
        *line_qual_data = " "; // pointer to static string

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    else if (flag_optimize) {
        optimize_phred_quality_string (*line_qual_data, *line_qual_len);
        if (*line_u2_data) optimize_phred_quality_string (*line_u2_data, *line_u2_len);
    }
}

// callback function for compress to get data of one line
static void zip_sam_get_start_len_line_i_bd (VBlock *vb, uint32_t vb_line_i, 
                                             char **line_bd_data, uint32_t *line_bd_len,  // out 
                                             char **unused1,  uint32_t *unused2)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    *line_bd_data = dl->bd_data_len ? ENT (char, vb->txt_data, dl->bd_data_start) : NULL;
    *line_bd_len  = dl->bd_data_len;
    *unused1 = NULL;
    *unused2 = 0;
}   

// callback function for compress to get BI data of one line
// if BD data is available for this line too - we output the character-wise delta as they are correlated
// we expect the length of BI and BD to be the same, the length of the sequence. this is enforced by seg_sam_optional_field.
static void zip_sam_get_start_len_line_i_bi (VBlock *vb, uint32_t vb_line_i, 
                                             char **line_bi_data, uint32_t *line_bi_len,  // out 
                                             char **unused1,  uint32_t *unused2)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    if (dl->bi_data_len && dl->bd_data_len) {

        ASSERT (dl->bi_data_len == dl->bd_data_len, "Error: expecting dl->bi_data_len=%u to be equal to dl->bd_data_len=%u",
                dl->bi_data_len, dl->bd_data_len);

        char *bi       = ENT (char, vb->txt_data, dl->bi_data_start);
        const char *bd = ENT (char, vb->txt_data, dl->bd_data_start);

        // calculate character-wise delta
        for (unsigned i=0; i < dl->bi_data_len; i++) *(bi++) -= *(bd++);
    }

    *line_bi_data = dl->bi_data_len ? ENT (char, vb->txt_data, dl->bi_data_start) : NULL;
    *line_bi_len  = dl->bi_data_len;
    *unused1 = NULL;
    *unused2 = 0;
}   

// this function receives all lines of a vblock and processes them
// in memory to the compressed format. This thread then terminates, and the I/O thread writes the output.
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
    seg_all_data_lines (vb_, seg_sam_data_line, seg_sam_initialize, sizeof (ZipDataLineSAM));

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

    zfile_compress_generic_vb_header (vb_); // vblock header

    // merge new words added in this vb into the z_file.mtf_ctx, ahead of zip_generate_b250_section() and
    // zip_vcf_generate_genotype_one_section(). writing indices based on the merged dictionaries. dictionaries are compressed. 
    // all this is done while holding exclusive access to the z_file dictionaries.
    mtf_merge_in_vb_ctx(vb_);

    // now, we merge vb->ra_buf into z_file->ra_buf
    random_access_merge_in_vb (vb_);

    // generate & write b250 data for all fields (FLAG to OPTIONAL)
    zip_generate_and_compress_fields (vb_);;

    // generate & write b250 data for all QNAME subfields
    zip_generate_and_compress_subfields (vb_, &vb->qname_mapper);

    // generate & write b250 data for all OPTIONAL subfields
    for (unsigned did_i=datatype_last_field[DT_SAM]+1; did_i < vb->num_dict_ids; did_i++) {
                
        MtfContext *ctx = &vb->mtf_ctx[did_i];
        
        if (ctx->mtf_i.len && ctx->dict_section_type == SEC_SAM_OPTNL_SF_DICT) {
            zip_generate_b250_section (vb_, ctx);
            zfile_compress_b250_data (vb_, ctx);
        }
    }

    // generate & compress Random POS data
    vb->random_pos_data.len *= sizeof (uint32_t);
    zfile_compress_section_data_alg (vb_, SEC_SAM_RAND_POS_DATA,  &vb->random_pos_data, NULL, 0, COMP_LZMA);
        
    // generate & compress the MD, BD, BI data
    if (vb->md_data.len)
        zfile_compress_section_data_alg (vb_, SEC_SAM_MD_DATA, &vb->md_data, NULL, 0, COMP_BZ2);

    // generate & compress the BD, BI, SEQ & QUAL data
    uint32_t seq_len, qual_len, bd_len, bi_len;
    zip_sam_get_data_lengths (vb, &seq_len, &qual_len, &bd_len, &bi_len);

    if (bd_len)
        zfile_compress_section_data_alg (vb_, SEC_SAM_BD_DATA, NULL, zip_sam_get_start_len_line_i_bd, bd_len, COMP_LZMA);
    
    if (bi_len)
        zfile_compress_section_data_alg (vb_, SEC_SAM_BI_DATA, NULL, zip_sam_get_start_len_line_i_bi, bi_len, COMP_LZMA);
    
    zfile_compress_section_data_alg (vb_, SEC_SEQ_DATA,  NULL, zip_sam_get_start_len_line_i_seq,  seq_len,  COMP_LZMA);
    zfile_compress_section_data_alg (vb_, SEC_QUAL_DATA, NULL, zip_sam_get_start_len_line_i_qual, qual_len, COMP_BZ2);

    // tell dispatcher this thread is done and can be joined. 
    // thread safety: this isn't protected by a mutex as it will just be false until it at some point turns to true
    // this this operation needn't be atomic, but it likely is anyway
    vb->is_processed = true; 

    COPY_TIMER (vb->profile.compute);
  }
