// ------------------------------------------------------------------
//   random_access.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "buffer.h"
#include "random_access.h"
#include "vb.h"
#include "file.h"
#include "endianness.h"
#include "regions.h"

// we pack this, because it gets written to disk in SEC_RANDOM_ACCESS. On disk, it is stored in Big Endian.

#pragma pack(push, 1) // structures that are part of the genozip format are packed.


#pragma pack(pop)

// ZIP only: called from seg_chrom_field when the CHROM changed - this might be a new chrom, or
// might be an exiting chrom if the VB is not sorted. we maitain one ra field per chrom per vb
void random_access_update_chrom (VariantBlock *vb, uint32_t vb_line_i, int32_t chrom_node_index)
{
    if (vb->curr_ra_ent && vb->curr_ra_ent->chrom == chrom_node_index) return; // all good - current chrom continues

    RAEntry *ra = ((RAEntry *)vb->ra_buf.data);

    // search for existing chrom
    unsigned i=0; for (; i < vb->ra_buf.len; i++)
        if (ra[i].chrom == chrom_node_index) { // found existing chrom
            vb->curr_ra_ent = &ra[i];
            return;
        }

    // this is a new chrom - we need a new entry
    buf_alloc (vb, &vb->ra_buf, sizeof (RAEntry) * (vb->ra_buf.len + 1), 2, "ra_buf", vb->variant_block_i);

    vb->curr_ra_ent = &((RAEntry *)vb->ra_buf.data)[vb->ra_buf.len++];
    memset (vb->curr_ra_ent, 0, sizeof(RAEntry));

    vb->curr_ra_ent->chrom           = chrom_node_index;
    vb->curr_ra_ent->start_vb_line   = vb_line_i;
    vb->curr_ra_ent->variant_block_i = vb->variant_block_i;
    vb->curr_ra_ent_is_initialized   = false; // we will finish the initialization or POS on the first call to random_access_update_pos
}

// ZIP only: called from seg_pos_field - update the pos in the existing chrom entry
void random_access_update_pos (VariantBlock *vb, int32_t this_pos)
{
    RAEntry *ra_ent = &((RAEntry *)vb->ra_buf.data)[vb->ra_buf.len-1];

    if (!vb->curr_ra_ent_is_initialized) {
        ra_ent->min_pos = ra_ent->max_pos = this_pos;
        vb->curr_ra_ent_is_initialized = true;
    }    

    else if (this_pos < ra_ent->min_pos) ra_ent->min_pos = this_pos; 
    
    else if (this_pos > ra_ent->max_pos) ra_ent->max_pos = this_pos;

    ra_ent->num_vb_lines++;
}

// called by ZIP compute thread, while holding the z_file mutex: merge in the VB's ra_buf in the global z_file one
void random_access_merge_in_vb (VariantBlock *vb)
{
    buf_alloc (vb, &vb->z_file->ra_buf, (vb->z_file->ra_buf.len + vb->ra_buf.len) * sizeof(RAEntry), 2, "z_file->ra_buf", 0);

    RAEntry *dst_ra = &((RAEntry *)vb->z_file->ra_buf.data)[vb->z_file->ra_buf.len];
    RAEntry *src_ra = ((RAEntry *)vb->ra_buf.data);

    MtfContext *chrom_ctx = &vb->mtf_ctx[CHROM];
    ASSERT0 (chrom_ctx, "Error: cannot find chrom_ctx");

    for (unsigned i=0; i < vb->ra_buf.len; i++) {
        MtfNode *chrom_node = mtf_node (chrom_ctx, src_ra[i].chrom, NULL, NULL);

        dst_ra[i].variant_block_i = vb->variant_block_i;
        dst_ra[i].chrom           = chrom_node->word_index.n; // note: in the VB we store the node index, while in zfile we store tha word index
        dst_ra[i].min_pos         = src_ra[i].min_pos;
        dst_ra[i].max_pos         = src_ra[i].max_pos;
        dst_ra[i].start_vb_line   = src_ra[i].start_vb_line;
        dst_ra[i].num_vb_lines    = src_ra[i].num_vb_lines;
    }

    vb->z_file->ra_buf.len += vb->ra_buf.len;
}

// PIZ I/O thread: check if for the given VB,
// the ranges in random access (from the file) overlap with the ranges in regions (from the command line -r or -R)
bool random_access_is_vb_included (uint32_t vb_i)
{
    if (!flag_regions) return true; // if no -r/-R was specified, all VBs are included

    static const RAEntry *next_ra = NULL;

    ASSERT0 ((vb_i==1) == !next_ra, "Error: expecting next_ra==NULL iff vb_i==1");

    const RAEntry *first_ra = (const RAEntry *)evb->z_file->ra_buf.data;
    if (!next_ra) next_ra = first_ra;

    bool included=false;
    while ((next_ra - first_ra < evb->z_file->ra_buf.len) && (next_ra->variant_block_i == vb_i)) {

        // we call regions_is_region_included_in_requested_regions() until a first included is found. then we just loop until the end of this vb
        if (!included && regions_is_region_included_in_requested_regions (next_ra->chrom, next_ra->min_pos, next_ra->max_pos))
            included = true; 

        next_ra++;
    }   

    return included; 
}

// Called by PIZ I/O thread (piz_read_global_area) and ZIP I/O thread (zip_write_global_area)
void BGEN_random_access()
{
    RAEntry *ra = (RAEntry *)evb->z_file->ra_buf.data;

    for (unsigned i=0; i < evb->z_file->ra_buf.len; i++) {
        ra[i].variant_block_i = BGEN32 (ra[i].variant_block_i);
        ra[i].chrom           = BGEN32 (ra[i].chrom);
        ra[i].min_pos         = BGEN32 (ra[i].min_pos);
        ra[i].max_pos         = BGEN32 (ra[i].max_pos);
        ra[i].start_vb_line   = BGEN32 (ra[i].start_vb_line);
        ra[i].num_vb_lines    = BGEN32 (ra[i].num_vb_lines);
    }
}

unsigned random_access_sizeof_entry()
{
    return sizeof (RAEntry);
}

void random_access_show_index ()
{
    printf ("Random-access index contents (result of --show-index):\n");
    
    const Buffer *ra_buf = &evb->z_file->ra_buf;

    for (unsigned i=0; i < ra_buf->len; i++) {
        
        RAEntry ra_ent = ((RAEntry *)ra_buf->data)[i]; // make a copy in case we need to BGEN
        
        printf ("vb_i=%u chrom_node_index=%u min_pos=%u max_pos=%u start_vb_line=%u num_vb_lines=%u\n",
                ra_ent.variant_block_i, ra_ent.chrom, ra_ent.min_pos, ra_ent.max_pos, ra_ent.start_vb_line, 
                ra_ent.num_vb_lines);
    }
}

