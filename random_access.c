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

    MtfNode *chrom_mtf = (MtfNode *)vb->mtf_ctx[CHROM].ol_mtf.data;

    for (unsigned i=0; i < vb->ra_buf.len; i++) {
        dst_ra[i].is_sorted       = src_ra[i].is_sorted;
        dst_ra[i].variant_block_i = BGEN32 (vb->variant_block_i);
        dst_ra[i].chrom           = BGEN32 (chrom_mtf[src_ra[i].chrom].word_index.n);
        dst_ra[i].first_pos       = BGEN32 (src_ra[i].first_pos);
        dst_ra[i].last_pos        = BGEN32 (src_ra[i].last_pos - src_ra[i].first_pos); // delta compresses better ?? check
        dst_ra[i].start_vb_line   = BGEN32 (src_ra[i].start_vb_line);
        dst_ra[i].num_vb_lines    = BGEN32 (src_ra[i].num_vb_lines);
    }

    vb->z_file->ra_buf.len += vb->ra_buf.len;
}

void BGEN_random_access (Buffer *ra_buf)
{
    RAEntry *ra = (RAEntry *)ra_buf->data;

    for (unsigned i=0; i < ra_buf->len; i++) {
        ra[i].variant_block_i = BGEN32 (ra[i].variant_block_i);
        ra[i].chrom           = BGEN32 (ra[i].chrom);
        ra[i].first_pos       = BGEN32 (ra[i].first_pos);
        ra[i].last_pos        = BGEN32 (ra[i].last_pos) + ra[i].first_pos; // restore from delta
        ra[i].start_vb_line   = BGEN32 (ra[i].start_vb_line);
        ra[i].num_vb_lines    = BGEN32 (ra[i].num_vb_lines);
    }
}