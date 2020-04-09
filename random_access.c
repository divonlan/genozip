// ------------------------------------------------------------------
//   random_access.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <pthread.h>
#include "buffer.h"
#include "random_access.h"
#include "vcf_vb.h"
#include "file.h"
#include "endianness.h"
#include "regions.h"

static pthread_mutex_t ra_mutex;
static bool ra_mutex_initialized = false;

void random_access_initialize(void)
{
    if (ra_mutex_initialized) return;

    unsigned ret = pthread_mutex_init (&ra_mutex, NULL);
    ASSERT0 (!ret, "Error: pthread_mutex_init failed for ra_mutex");

    ra_mutex_initialized = true;
}

// ZIP only: called from seg_chrom_field when the CHROM changed - this might be a new chrom, or
// might be an exiting chrom if the VB is not sorted. we maitain one ra field per chrom per vb
void random_access_update_chrom (VariantBlock *vb, uint32_t vb_line_i, int32_t chrom_node_index)
{
    ASSERT (chrom_node_index >= 0, "Error in random_access_update_chrom: chrom_node_index=%d in vb_i=%u", 
            chrom_node_index, vb->variant_block_i);

    if (vb->curr_ra_ent && vb->curr_ra_ent->chrom_index == chrom_node_index) return; // all good - current chrom continues

    RAEntry *ra = ((RAEntry *)vb->ra_buf.data);

    // search for existing chrom
    unsigned i=0; for (; i < vb->ra_buf.len; i++)
        if (ra[i].chrom_index == chrom_node_index) { // found existing chrom
            vb->curr_ra_ent = &ra[i];
            return;
        }

    // this is a new chrom - we need a new entry
    buf_alloc (vb, &vb->ra_buf, sizeof (RAEntry) * (vb->ra_buf.len + 1), 2, "ra_buf", vb->variant_block_i);

    vb->curr_ra_ent = &((RAEntry *)vb->ra_buf.data)[vb->ra_buf.len++];
    memset (vb->curr_ra_ent, 0, sizeof(RAEntry));

    vb->curr_ra_ent->chrom_index     = chrom_node_index;
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
}

// called by ZIP compute thread, while holding the z_file mutex: merge in the VB's ra_buf in the global z_file one
// note: the order of the merge is not necessarily the sequential order of VBs
void random_access_merge_in_vb (VariantBlock *vb)
{
    pthread_mutex_lock (&ra_mutex);

    buf_alloc (evb, &z_file->ra_buf, (z_file->ra_buf.len + vb->ra_buf.len) * sizeof(RAEntry), 2, "z_file->ra_buf", 0);

    RAEntry *dst_ra = &((RAEntry *)z_file->ra_buf.data)[z_file->ra_buf.len];
    RAEntry *src_ra = ((RAEntry *)vb->ra_buf.data);

    MtfContext *chrom_ctx = &vb->mtf_ctx[CHROM];
    ASSERT0 (chrom_ctx, "Error: cannot find chrom_ctx");

    for (unsigned i=0; i < vb->ra_buf.len; i++) {
        MtfNode *chrom_node = mtf_node (chrom_ctx, src_ra[i].chrom_index, NULL, NULL);

        dst_ra[i].variant_block_i = vb->variant_block_i;
        dst_ra[i].chrom_index     = chrom_node->word_index.n; // note: in the VB we store the node index, while in zfile we store tha word index
        dst_ra[i].min_pos         = src_ra[i].min_pos;
        dst_ra[i].max_pos         = src_ra[i].max_pos;
    }

    z_file->ra_buf.len += vb->ra_buf.len;

    pthread_mutex_unlock (&ra_mutex);
}

// PIZ I/O thread: check if for the given VB,
// the ranges in random access (from the file) overlap with the ranges in regions (from the command line -r or -R)
bool random_access_is_vb_included (uint32_t vb_i,
                                   Buffer *region_ra_intersection_matrix) // out - a bytemap - rows are ra's of this VB, columns are regions, a cell is 1 if there's an intersection
{
    if (!flag_regions) return true; // if no -r/-R was specified, all VBs are included

    ASSERT0 (region_ra_intersection_matrix, "Error: region_ra_intersection_matrix is NULL");

    // allocate bytemap. note that it allocate in evb, and will be buf_moved to the vb after it is generated 
    ASSERT (!buf_is_allocated (region_ra_intersection_matrix), "Error: expecting region_ra_intersection_matrix to be unallcoated vb_i=%u", vb_i);

    unsigned num_regions = regions_max_num_chregs();
    buf_alloc (evb, region_ra_intersection_matrix, z_file->ra_buf.len * num_regions, 1, "region_ra_intersection_matrix", vb_i);
    buf_zero (region_ra_intersection_matrix);

    static const RAEntry *next_ra = NULL;
    ASSERT0 ((vb_i==1) == !next_ra, "Error: expecting next_ra==NULL iff vb_i==1");
    if (!next_ra) next_ra = (const RAEntry *)z_file->ra_buf.data; // initialize static on first call


    bool vb_is_included=false;
    for (unsigned ra_i=0; 
         ra_i < z_file->ra_buf.len && next_ra->variant_block_i == vb_i;
         ra_i++, next_ra++) {

        if (regions_get_ra_intersection (next_ra->chrom_index, next_ra->min_pos, next_ra->max_pos,
                                         &region_ra_intersection_matrix->data[ra_i * num_regions]))  // the matrix row for this ra
            vb_is_included = true; 
    }   

    if (!vb_is_included) buf_free (region_ra_intersection_matrix);

    return vb_is_included; 
}

// PIZ I/O threads: get last vb_i that is included in the regions requested with --regions, or -1 if no vb includes regions.
// used for reading dictionaries from a genozip file - there's no need to read dictionaries beyond this vb_i
int32_t random_access_get_last_included_vb_i (void)
{
    int32_t last_vb_i = -1;
    for (unsigned ra_i=0; ra_i < z_file->ra_buf.len; ra_i++) { // note that all entries of the same vb_i are together, but vb_i's are not necessarily in seqential order
        
        const RAEntry *ra = ENT (RAEntry, &z_file->ra_buf, ra_i);
        if ((int32_t)ra->variant_block_i <= last_vb_i) continue; // we already decided to include this vb_i - no need to check further

        if (regions_get_ra_intersection (ra->chrom_index, ra->min_pos, ra->max_pos, NULL))
            last_vb_i = (int32_t)ra->variant_block_i; 
    }   
    return last_vb_i;
}

// Called by PIZ I/O thread (piz_read_global_area) and ZIP I/O thread (zip_write_global_area)
void BGEN_random_access()
{
    RAEntry *ra = (RAEntry *)z_file->ra_buf.data;

    for (unsigned i=0; i < z_file->ra_buf.len; i++) {
        ra[i].variant_block_i = BGEN32 (ra[i].variant_block_i);
        ra[i].chrom_index     = BGEN32 (ra[i].chrom_index);
        ra[i].min_pos         = BGEN32 (ra[i].min_pos);
        ra[i].max_pos         = BGEN32 (ra[i].max_pos);
    }
}

unsigned random_access_sizeof_entry()
{
    return sizeof (RAEntry);
}

void random_access_show_index ()
{
    fprintf (stderr, "Random-access index contents (result of --show-index):\n");
    
    const Buffer *ra_buf = &z_file->ra_buf;

    for (unsigned i=0; i < ra_buf->len; i++) {
        
        RAEntry ra_ent = ((RAEntry *)ra_buf->data)[i]; // make a copy in case we need to BGEN
        
        fprintf (stderr, "vb_i=%u chrom_node_index=%u min_pos=%u max_pos=%u\n",
                 ra_ent.variant_block_i, ra_ent.chrom_index, ra_ent.min_pos, ra_ent.max_pos);
    }
}

