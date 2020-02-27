// ------------------------------------------------------------------
//   random_access.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "buffer.h"
#include "random_access.h"
#include "vb.h"
#include "file.h"
#include "endianness.h"

// we pack this, because it gets written to disk in SEC_RANDOM_ACCESS. On disk, it is stored in Big Endian.

#pragma pack(push, 1) // structures that are part of the genozip format are packed.

typedef struct {
    uint32_t chrom;                       // before merge: node index into chrom context mtf, after merge - word index in CHROM dictionary
    uint32_t min_pos;                     // POS field value of first position
    int32_t  max_pos;                     // in VB, this is in an absolute value. On disk, it is the detla vs first_pos
    uint32_t variant_block_i;             // the vb_i in which this range appears
    uint32_t start_vb_line, num_vb_lines; // corresponds to the line with this VB
    uint32_t is_sorted      : 1;          // is this range sorted in a non-decreasing order?
    uint32_t for_future_use : 31;         // put this 1 bit in 32 bits to avoid aliasing issues
} RAEntry; 

#pragma pack(pop)

int32_t random_access_get_last_chrom_node_index (VariantBlock *vb)
{
    ASSERT (vb->ra_buf.len > 0, "Error: vb->ra_buf.len is 0, can't get last, vb_i=%u", vb->variant_block_i);

    return ((RAEntry *)vb->ra_buf.data)[vb->ra_buf.len-1].chrom;
}

// ZIP only: called from seg_chrom_field
void random_access_new_entry (VariantBlock *vb, uint32_t vb_line_i, int32_t chrom_node_index)
{
    buf_alloc (vb, &vb->ra_buf, sizeof (RAEntry) * (vb->ra_buf.len + 1), 2, "ra_buf", vb->variant_block_i);

    RAEntry *ra_ent = &((RAEntry *)vb->ra_buf.data)[vb->ra_buf.len++];
    memset (ra_ent, 0, sizeof(RAEntry));

    ra_ent->chrom           = (chrom_node_index >= 0) ? chrom_node_index : (ra_ent-1)->chrom; // copy previous chrom if -1 (=unchanged)
    ra_ent->start_vb_line   = vb_line_i;
    ra_ent->is_sorted       = true; // sorted until proven otherwise
    ra_ent->variant_block_i = vb->variant_block_i;
}

// ZIP only: called from seg_pos_field
void random_access_update_last_entry (VariantBlock *vb, int32_t this_pos)
{
    RAEntry *ra_ent = &((RAEntry *)vb->ra_buf.data)[vb->ra_buf.len-1];

    if (this_pos < ra_ent->max_pos) ra_ent->is_sorted = false;
    
    if (ra_ent->min_pos == 0 || this_pos < ra_ent->min_pos) 
        ra_ent->min_pos = this_pos; // new vb or new new chrom 
    
    if (this_pos > ra_ent->max_pos)
        ra_ent->max_pos = this_pos;

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

        dst_ra[i].is_sorted       = src_ra[i].is_sorted;
        dst_ra[i].variant_block_i = BGEN32 (vb->variant_block_i);
        dst_ra[i].chrom           = BGEN32 (chrom_node->word_index.n); // note: in the VB we store the node index, while in zfile we store tha word index
        dst_ra[i].min_pos         = BGEN32 (src_ra[i].min_pos);
        dst_ra[i].max_pos         = BGEN32 (src_ra[i].max_pos - src_ra[i].min_pos); // delta compresses better ?? check
        dst_ra[i].start_vb_line   = BGEN32 (src_ra[i].start_vb_line);
        dst_ra[i].num_vb_lines    = BGEN32 (src_ra[i].num_vb_lines);
    }

    vb->z_file->ra_buf.len += vb->ra_buf.len;
}

// converts back to native from big endian, and also recovers max_pos
static inline void BGEN_random_access_entry (RAEntry *ra_ent)
{
    ra_ent->variant_block_i = BGEN32 (ra_ent->variant_block_i);
    ra_ent->chrom           = BGEN32 (ra_ent->chrom);
    ra_ent->min_pos         = BGEN32 (ra_ent->min_pos);
    ra_ent->max_pos         = BGEN32 (ra_ent->max_pos) + ra_ent->min_pos; // restore from delta
    ra_ent->start_vb_line   = BGEN32 (ra_ent->start_vb_line);
    ra_ent->num_vb_lines    = BGEN32 (ra_ent->num_vb_lines);
}

// Called by PIZ I/O thread
void BGEN_random_access (Buffer *ra_buf)
{
    RAEntry *ra = (RAEntry *)ra_buf->data;

    for (unsigned i=0; i < ra_buf->len; i++) 
        BGEN_random_access_entry (&ra[i]);
}

unsigned random_access_sizeof_entry()
{
    return sizeof (RAEntry);
}

void random_access_show (const Buffer *ra_buf, bool is_big_endian)
{
    fprintf (stderr, "Random-access index contents (result of --show-index):\n");

    for (unsigned i=0; i < ra_buf->len; i++) {
        
        RAEntry ra_ent = ((RAEntry *)ra_buf->data)[i]; // make a copy in case we need to BGEN
        
        if (is_big_endian) BGEN_random_access_entry (&ra_ent);

        const char *chrom_str = is_big_endian ? "chrom_WORD_index" : "chrom_NODE_index";

        fprintf (stderr, "vb_i=%u %s=%u min_pos=%u max_pos=%u start_vb_line=%u num_vb_lines=%u is_sorted=%u\n",
                 ra_ent.variant_block_i, chrom_str, ra_ent.chrom, ra_ent.min_pos, ra_ent.max_pos, ra_ent.start_vb_line, 
                 ra_ent.num_vb_lines, ra_ent.is_sorted);
    }
}