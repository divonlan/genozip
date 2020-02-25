// ------------------------------------------------------------------
//   random_access.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "buffer.h"
#include "random_access.h"
#include "vb.h"
#include "file.h"
#include "endianness.h"

// called by ZIP compute thread, while holding the z_file mutex: merge in the VB's ra_buf in the global z_file one
void random_access_merge_in_vb (VariantBlock *vb)
{
    buf_alloc (vb, &vb->z_file->ra_buf, (vb->z_file->ra_buf.len + vb->ra_buf.len) * sizeof(RAEntry), 2, "z_file->ra_buf", 0);

    RAEntry *dst_ra = &((RAEntry *)vb->z_file->ra_buf.data)[vb->z_file->ra_buf.len];
    RAEntry *src_ra = ((RAEntry *)vb->ra_buf.data);

    MtfNode *chrom_mtf = (MtfNode *)vb->mtf_ctx[CHROM].mtf.data;

    for (unsigned i=0; i < vb->ra_buf.len; i++) {
        dst_ra->is_sorted       = src_ra->is_sorted;
        dst_ra->variant_block_i = BGEN32 (vb->variant_block_i);
        dst_ra->chrom           = BGEN32 (chrom_mtf[src_ra->chrom].word_index.n);
        dst_ra->first_pos       = BGEN32 (src_ra->first_pos);
        dst_ra->last_pos        = BGEN32 (src_ra->last_pos - src_ra->first_pos); // delta compresses better ?? check
        dst_ra->start_vb_line   = BGEN32 (src_ra->start_vb_line);
        dst_ra->num_vb_lines    = BGEN32 (src_ra->num_vb_lines);
    }
}