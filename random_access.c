// ------------------------------------------------------------------
//   random_access.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <pthread.h>
#include "buffer.h"
#include "random_access.h"
#include "vblock.h"
#include "file.h"
#include "endianness.h"
#include "regions.h"
#include "sections.h"

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
void random_access_update_chrom (VBlock *vb, uint32_t vb_line_i, int32_t chrom_node_index)
{
    ASSERT (chrom_node_index >= 0, "Error in random_access_update_chrom: chrom_node_index=%d in vb_i=%u", 
            chrom_node_index, vb->vblock_i);

    // make sure ra_buf is big enough
    uint64_t old_len = vb->ra_buf.len;
    uint64_t new_len = chrom_node_index+1;
    if (new_len > old_len) {
        buf_alloc (vb, &vb->ra_buf, sizeof (RAEntry) * MAX (new_len, 500), 2, "ra_buf", vb->vblock_i);
        memset (ENT (RAEntry, vb->ra_buf, old_len), 0, (new_len - old_len) * sizeof (RAEntry));
        vb->ra_buf.len = new_len;
    }

    RAEntry *ra_ent = ENT (RAEntry, vb->ra_buf, chrom_node_index);
    if (!ra_ent->vblock_i) { // first occurance of this chrom
        ra_ent->chrom_index = chrom_node_index;
        ra_ent->vblock_i    = vb->vblock_i;
    }

    vb->chrom_node_index = chrom_node_index;
}

// ZIP only: called from seg_pos_field - update the pos in the existing chrom entry
void random_access_update_pos (VBlock *vb, int32_t this_pos)
{
    if (!this_pos) return; // ignore pos=0 (in SAM, it means unmapped POS)

    RAEntry *ra_ent = ENT (RAEntry, vb->ra_buf, vb->chrom_node_index);

    if (!ra_ent->min_pos) // first line this chrom is encountered in this vb - initialize pos
        ra_ent->min_pos = ra_ent->max_pos = this_pos;

    else if (this_pos < ra_ent->min_pos) ra_ent->min_pos = this_pos; 
    
    else if (this_pos > ra_ent->max_pos) ra_ent->max_pos = this_pos;
}

// called by ZIP compute thread, while holding the z_file mutex: merge in the VB's ra_buf in the global z_file one
// note: the order of the merge is not necessarily the sequential order of VBs
void random_access_merge_in_vb (VBlock *vb)
{
    pthread_mutex_lock (&ra_mutex);

    buf_alloc (evb, &z_file->ra_buf, (z_file->ra_buf.len + vb->ra_buf.len) * sizeof(RAEntry), 2, "z_file->ra_buf", 0);

    ARRAY (RAEntry, src_ra, vb->ra_buf);

    MtfContext *chrom_ctx = &vb->mtf_ctx[chrom_did_i_by_dt[vb->data_type]];
    ASSERT0 (chrom_ctx, "Error in random_access_merge_in_vb: cannot find chrom_ctx");

    for (unsigned i=0; i < vb->ra_buf.len; i++) {
        
        if (!src_ra[i].min_pos) continue; // chrom node_index=i has no range in this vb

        RAEntry *dst_ra = &NEXTENT (RAEntry, z_file->ra_buf);

        MtfNode *chrom_node = mtf_node (chrom_ctx, src_ra[i].chrom_index, NULL, NULL);

        dst_ra->vblock_i    = vb->vblock_i;
        dst_ra->chrom_index = chrom_node->word_index.n; // note: in the VB we store the node index, while in zfile we store tha word index
        dst_ra->min_pos     = src_ra[i].min_pos;
        dst_ra->max_pos     = src_ra[i].max_pos;
    }

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
         ra_i < z_file->ra_buf.len && next_ra->vblock_i == vb_i;
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
        
        const RAEntry *ra = ENT (RAEntry, z_file->ra_buf, ra_i);
        if ((int32_t)ra->vblock_i <= last_vb_i) continue; // we already decided to include this vb_i - no need to check further

        if (regions_get_ra_intersection (ra->chrom_index, ra->min_pos, ra->max_pos, NULL))
            last_vb_i = (int32_t)ra->vblock_i; 
    }   
    return last_vb_i;
}

// Called by PIZ I/O thread (piz_read_global_area) and ZIP I/O thread (zip_write_global_area)
void BGEN_random_access()
{
    ARRAY (RAEntry, ra, z_file->ra_buf);

    for (unsigned i=0; i < z_file->ra_buf.len; i++) {
        ra[i].vblock_i    = BGEN32 (ra[i].vblock_i);
        ra[i].chrom_index = BGEN32 (ra[i].chrom_index);
        ra[i].min_pos     = BGEN32 (ra[i].min_pos);
        ra[i].max_pos     = BGEN32 (ra[i].max_pos);
    }
}

unsigned random_access_sizeof_entry()
{
    return sizeof (RAEntry);
}

void random_access_show_index (bool from_zip)
{
    fprintf (stderr, "Random-access index contents (result of --show-index):\n");
    
    ARRAY (RAEntry, ra, z_file->ra_buf);

    MtfContext *ctx = &z_file->mtf_ctx[chrom_did_i_by_dt[z_file->data_type]];

    for (unsigned i=0; i < z_file->ra_buf.len; i++) {
        
        const char *chrom_snip; unsigned chrom_snip_len;
        if (from_zip) {
            MtfNode *chrom_node = mtf_get_node_by_word_index (ctx, ra[i].chrom_index);
            chrom_snip    = ENT (char, ctx->dict, chrom_node->char_index);
            chrom_snip_len = chrom_node->snip_len;
        }
        else {
            MtfWord *chrom_word = ENT (MtfWord, ctx->word_list, ra[i].chrom_index);
            chrom_snip = ENT (char, ctx->dict, chrom_word->char_index);
            chrom_snip_len = chrom_word->snip_len;
        }
        fprintf (stderr, "vb_i=%u chrom='%.*s' (chrom_word_index=%u) min_pos=%u max_pos=%u\n",
                    ra[i].vblock_i, chrom_snip_len, chrom_snip, ra[i].chrom_index, ra[i].min_pos, ra[i].max_pos);
    }
}

