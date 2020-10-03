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
#include "mutex.h"
#include "zfile.h"

MUTEX (ra_mutex);

void random_access_initialize(void)
{    
    mutex_initialize (ra_mutex);
}

// --------------------
// ZIP stuff
// --------------------

// ZIP only
void random_access_alloc_ra_buf (VBlock *vb, int32_t chrom_node_index)
{
    uint64_t old_len = vb->ra_buf.len;
    uint64_t new_len = chrom_node_index + 2; // +2 because we store for chrom_node_index [-1, chrom_node_index]
    if (new_len > old_len) {
        buf_alloc (vb, &vb->ra_buf, sizeof (RAEntry) * MAX (new_len, 500), 2, "ra_buf", vb->vblock_i);
        memset (ENT (RAEntry, vb->ra_buf, old_len), 0, (new_len - old_len) * sizeof (RAEntry));
        vb->ra_buf.len = new_len;
    }
}

// ZIP only: called from vcf_seg_chrom_field when the CHROM changed - this might be a new chrom, or
// might be an exiting chrom (for example, in an unsorted SAM or VCF). we maitain one ra field per chrom per vb
void random_access_update_chrom (VBlock *vb, WordIndex chrom_node_index, const char *chrom_name, unsigned chrom_name_len)
{
    // note: when FASTA calls this for a sequence that started in the previous vb, and hence chrom is unknown, chrom_node_index==-1.
    ASSERT (chrom_node_index >= -1, "Error in random_access_update_chrom: chrom_node_index=%d in vb_i=%u", 
            chrom_node_index, vb->vblock_i);

    random_access_alloc_ra_buf (vb, chrom_node_index); // make sure ra_buf is big enough

    RAEntry *ra_ent = ENT (RAEntry, vb->ra_buf, chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc
    if (!ra_ent->vblock_i) { // first occurance of this chrom
        ra_ent->chrom_index = chrom_node_index;
        ra_ent->vblock_i    = vb->vblock_i;
    }

    vb->chrom_node_index = chrom_node_index;

    if (chrom_node_index != WORD_INDEX_NONE) {
        vb->chrom_name       = chrom_name;
        vb->chrom_name_len   = chrom_name_len;
    }
}

// ZIP only: update the pos in the existing chrom entry
void random_access_update_pos (VBlock *vb, DidIType did_i_pos)
{
    PosType this_pos = vb->contexts[did_i_pos].last_value.i;

    if (this_pos <= 0) return; // ignore pos<=0 (in SAM, 0 means unmapped POS)

    RAEntry *ra_ent = ENT (RAEntry, vb->ra_buf, vb->chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc

    if (!ra_ent->min_pos) // first line this chrom is encountered in this vb - initialize pos
        ra_ent->min_pos = ra_ent->max_pos = this_pos;

    else if (this_pos < ra_ent->min_pos) ra_ent->min_pos = this_pos; 
    
    else if (this_pos > ra_ent->max_pos) ra_ent->max_pos = this_pos;
}

// ZIP: update increment reference pos
void random_access_increment_last_pos (VBlockP vb, PosType increment)
{
    RAEntry *ra_ent = ENT (RAEntry, vb->ra_buf, vb->chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc

    if (!ra_ent->min_pos) ra_ent->min_pos = 1; 
    ra_ent->max_pos += increment;
}

// ZIP: update last reference pos
void random_access_update_last_pos (VBlock *vb, PosType last_pos)
{
    RAEntry *ra_ent = ENT (RAEntry, vb->ra_buf, vb->chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc
    if (last_pos > ra_ent->max_pos) ra_ent->max_pos = last_pos;
}

void random_access_update_to_entire_chrom (VBlockP vb, PosType first_pos_of_chrom, PosType last_pos_of_chrom)
{
    RAEntry *ra_ent = ENT (RAEntry, vb->ra_buf, vb->chrom_node_index + 1); // chrom_node_index=-1 goes into entry 0 etc
    ra_ent->min_pos = first_pos_of_chrom;
    ra_ent->max_pos = last_pos_of_chrom;
}

// called by ZIP compute thread, while holding the z_file mutex: merge in the VB's ra_buf in the global z_file one
// note: the order of the merge is not necessarily the sequential order of VBs
void random_access_merge_in_vb (VBlock *vb)
{
    mutex_lock (ra_mutex);

    buf_alloc (evb, &z_file->ra_buf, (z_file->ra_buf.len + vb->ra_buf.len) * sizeof(RAEntry), 2, "z_file->ra_buf", 0); 

    ARRAY (RAEntry, src_ra, vb->ra_buf);

    Context *chrom_ctx = &vb->contexts[CHROM];
    ASSERT0 (chrom_ctx, "Error in random_access_merge_in_vb: cannot find chrom_ctx");

    for (unsigned i=0; i < vb->ra_buf.len; i++) {
        
        if (!src_ra[i].min_pos) continue; // chrom node_index=i has no range in this vb

        RAEntry *dst_ra = &NEXTENT (RAEntry, z_file->ra_buf);

        dst_ra->vblock_i = vb->vblock_i;
        dst_ra->min_pos  = src_ra[i].min_pos;
        dst_ra->max_pos  = src_ra[i].max_pos;

        if (src_ra[i].chrom_index != WORD_INDEX_NONE) {
            MtfNode *chrom_node = mtf_node_vb (chrom_ctx, (WordIndex)src_ra[i].chrom_index, NULL, NULL);
            dst_ra->chrom_index = chrom_node->word_index.n; // note: in the VB we store the node index, while in zfile we store tha word index
        }
        else 
            dst_ra->chrom_index = WORD_INDEX_NONE; // to be updated in random_access_finalize_entries()
    }

    mutex_unlock (ra_mutex);
}

// ZIP
static Buffer *ra_buf_being_sorted;
int random_access_sort_by_vb_i (const void *a_, const void *b_)
{
    ARRAY (RAEntry, ra, *ra_buf_being_sorted);

    int32_t a = *(int32_t *)a_;
    int32_t b = *(int32_t *)b_;
    
    int vb_i_diff = (int)ra[a].vblock_i - (int)ra[b].vblock_i;

    if (vb_i_diff)
        return vb_i_diff;   // from lower to higher vb_i
    else
        return a - b; // for same vb_i, keep order as it currently is in array (order will be from lower to higher address)
}

// ZIP (I/O thread) sort RA, update overflowing chroms, create and merge in evb ra
void random_access_finalize_entries (Buffer *ra_buf)
{
    // build an index into ra_buf that we will sort. we need that, because for same-vb entries we need to 
    // maintain their current order - sorted by the index
    int32_t *sorter = malloc (ra_buf->len * sizeof (int32_t));
    ASSERT (sorter, "Error in random_access_finalize_entries: failed to malloc %u bytes", (uint32_t)(ra_buf->len * sizeof (uint32_t)));

    for (int32_t i=0; i < ra_buf->len; i++) sorter[i] = i;

    // at this point, VBs in RA might be out of order, but all entries of a specific VB are grouped together in order
    // we sort so that the VBs are in order, and the order within a VB is unchanged
    ra_buf_being_sorted = ra_buf;
    qsort (sorter, ra_buf->len, sizeof (uint32_t), random_access_sort_by_vb_i);

    // use sorter to consturct a sorted RA
    static Buffer sorted_ra_buf = EMPTY_BUFFER; // must be static because its added to buf_list
    buf_alloc (evb, &sorted_ra_buf, sizeof (RAEntry) * ra_buf->len, 1, ra_buf->name, 0);
    sorted_ra_buf.len = ra_buf->len;

    for (uint32_t i=0; i < ra_buf->len; i++) 
        *ENT (RAEntry, sorted_ra_buf, i) = *ENT (RAEntry, *ra_buf, sorter[i]);

    // replace ra_buf with sorted one
    buf_destroy (ra_buf);
    buf_move (evb, ra_buf, evb, &sorted_ra_buf);

    FREE (sorter);

    // now that the VBs are in order, we can updated the "CONTINUED FROM PREVIOUS VB" ra's to their final values
    for (uint32_t i=0; i < ra_buf->len; i++) {
        
        RAEntry *ra = ENT (RAEntry, *ra_buf, i);

        // case: chrom is unknown - if the sequence started in the previous VB (this happens in FASTA) - we update now
        if (ra->chrom_index == WORD_INDEX_NONE) {
            // we expect this ra to be the first in its VB, and the previous ra to be of the previous VB
            ASSERT (i && (ra->vblock_i == (ra-1)->vblock_i+1), "Error in random_access_finalize_entries: currupt ra[%u]: chrom_index=WORD_INDEX_NONE but vb_i=%u and (ra-1)->vb_i=%u (expecting it to be %u)",
                    i, ra->vblock_i, (ra-1)->vblock_i, ra->vblock_i-1);

            ra->chrom_index = (ra-1)->chrom_index;
            ra->min_pos    += (ra-1)->max_pos;
            ra->max_pos    += (ra-1)->max_pos;
        }
    }
}

// --------------------
// PIZ stuff
// --------------------

// PIZ: binary search for first entry of vb_i, using the fact that in random_access_finalize_entries (in ZIP) we sorted the entries by vb_i
static const RAEntry *random_access_get_first_ra_of_vb (uint32_t vb_i, const RAEntry *first_ra, const RAEntry *last_ra)
{
    if (first_ra > last_ra) return NULL; // vb_i not found

    const RAEntry *mid_ra = first_ra + ((last_ra - first_ra) / 2);
    bool first_ra_in_vb = (mid_ra == FIRSTENT (RAEntry, z_file->ra_buf)) || (mid_ra->vblock_i != (mid_ra-1)->vblock_i);

    if (mid_ra->vblock_i < vb_i) 
        return random_access_get_first_ra_of_vb (vb_i, mid_ra+1, last_ra);

    if (mid_ra->vblock_i > vb_i || !first_ra_in_vb)
        return random_access_get_first_ra_of_vb (vb_i, first_ra, mid_ra-1);

    return mid_ra;
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

    const RAEntry *ra = random_access_get_first_ra_of_vb (vb_i, FIRSTENT (RAEntry, z_file->ra_buf), LASTENT (RAEntry, z_file->ra_buf));

    bool vb_is_included = false;
    for (unsigned ra_i=0; ra->vblock_i == vb_i; ra_i++, ra++) {
        if (regions_get_ra_intersection (ra->chrom_index, ra->min_pos, ra->max_pos,
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
    for (unsigned ra_i=0; ra_i < z_file->ra_buf.len; ra_i++) { 
        
        const RAEntry *ra = ENT (RAEntry, z_file->ra_buf, ra_i);
        if ((int32_t)ra->vblock_i <= last_vb_i) continue; // we already decided to include this vb_i - no need to check further

        if (regions_get_ra_intersection (ra->chrom_index, ra->min_pos, ra->max_pos, NULL))
            last_vb_i = (int32_t)ra->vblock_i; 
    }   
    return last_vb_i;
}

// PIZ I/O thread: gets min/max pos value for a particular chrom, across the entire file, by looking at the RA entries
void random_access_pos_of_chrom (WordIndex chrom_word_index, PosType *min_pos, PosType *max_pos)
{
    typedef struct { PosType min_pos, max_pos, gpos; } MinMax;

    // first time here for this z_file - we initialize ra_min_max_by_chrom
    if (!buf_is_allocated (&z_file->ra_min_max_by_chrom)) {
        uint64_t num_chroms = z_file->contexts[CHROM].word_list.len;

        buf_alloc (evb, &z_file->ra_min_max_by_chrom, num_chroms * sizeof (MinMax), 1, "z_file->ra_min_max_by_chrom", 0);
        buf_zero (&z_file->ra_min_max_by_chrom); // safety
        z_file->ra_min_max_by_chrom.len = num_chroms;

        // initialize
        for (unsigned i=0; i < num_chroms; i++) { 
            MinMax *mm = ENT (MinMax, z_file->ra_min_max_by_chrom, i);
            mm->min_pos = RA_MISSING_RA_MIN;
            mm->max_pos = RA_MISSING_RA_MAX;
        }

        // calculate min and max for each chrom by traversing the whole index (one entry per chrome per vb)
        for (unsigned i=0; i < z_file->ra_buf.len; i++) {
            RAEntry *ra = ENT (RAEntry, z_file->ra_buf, i);
            MinMax *mm = ENT (MinMax, z_file->ra_min_max_by_chrom, ra->chrom_index);
            mm->max_pos = MAX (mm->max_pos, ra->max_pos);
            mm->min_pos = MIN (mm->min_pos, ra->min_pos);
        }
    }

    MinMax *mm = ENT (MinMax, z_file->ra_min_max_by_chrom, chrom_word_index);
    if (min_pos) *min_pos = mm->min_pos;
    if (max_pos) *max_pos = mm->max_pos;
}

// FASTA PIZ compute thread when consuming a reference FASTA
// sets vb->{chrom_name,chrom_name_len,chrom_node_index} and returns start_pos
void random_access_get_first_chrom_of_vb (VBlockP vb, PosType *first_pos, PosType *last_pos)
{
    Context *ctx = &z_file->contexts[CHROM];
    ASSERT (ctx->word_list.len, "Error in random_access_get_first_chrom_of_vb: word_list of %s is empty", ctx->name);

    const RAEntry *ra = random_access_get_first_ra_of_vb (vb->vblock_i, FIRSTENT (RAEntry, z_file->ra_buf), LASTENT (RAEntry, z_file->ra_buf));
    ASSERT (ra, "Error in random_access_get_first_chrom_of_vb: vb_i=%u not found in random access index", vb->vblock_i);

    MtfWord *chrom_word  = ENT (MtfWord, ctx->word_list, ra->chrom_index);            
    vb->chrom_name       = ENT (const char, ctx->dict, chrom_word->char_index);
    vb->chrom_name_len   = chrom_word->snip_len;
    vb->chrom_node_index = ra->chrom_index;
    *first_pos           = ra->min_pos;
    *last_pos            = ra->max_pos;
}

// FASTA PIZ
bool random_access_does_last_chrom_continue_in_next_vb (uint32_t vb_i)
{
    const RAEntry *ra = random_access_get_first_ra_of_vb (vb_i+1, FIRSTENT (RAEntry, z_file->ra_buf), LASTENT (RAEntry, z_file->ra_buf));
    if (!ra) return false; // vb_i is the last vb (there is no vb_i+1) - so it doesn't continue in the next vb...

    return ra->chrom_index == (ra-1)->chrom_index;
}

// Called by PIZ I/O thread (piz_read_global_area) and ZIP I/O thread (zip_write_global_area)
void BGEN_random_access (Buffer *ra_buf)
{
    ARRAY (RAEntry, ra, *ra_buf);

    for (unsigned i=0; i < ra_buf->len; i++) {
        ra[i].vblock_i    = BGEN32 (ra[i].vblock_i);
        ra[i].chrom_index = BGEN32 (ra[i].chrom_index);
        ra[i].min_pos     = BGEN64 (ra[i].min_pos);
        ra[i].max_pos     = BGEN64 (ra[i].max_pos);
    }
}

void random_access_show_index (const Buffer *ra_buf, bool from_zip, const char *msg)
{
    fprintf (stderr, "\n%s:\n", msg);
    
    ARRAY (const RAEntry, ra, *ra_buf);

    Context *ctx = &z_file->contexts[CHROM];

    for (unsigned i=0; i < ra_buf->len; i++) {
        
        const char *chrom_snip; unsigned chrom_snip_len;
        if (from_zip) {
            if (ra[i].chrom_index != WORD_INDEX_NONE) {
                MtfNode *chrom_node = mtf_get_node_by_word_index (ctx, ra[i].chrom_index);
                chrom_snip     = ENT (char, ctx->dict, chrom_node->char_index);
                chrom_snip_len = chrom_node->snip_len;
            }
            else {
                chrom_snip     = "CONT_FROM_PREV_VB";
                chrom_snip_len = strlen (chrom_snip);
            }
        }
        else {
            MtfWord *chrom_word = ENT (MtfWord, ctx->word_list, ra[i].chrom_index);
            chrom_snip = ENT (char, ctx->dict, chrom_word->char_index);
            chrom_snip_len = chrom_word->snip_len;
        }
        fprintf (stderr, "vb_i=%u chrom='%.*s' (chrom_word_index=%d) min_pos=%"PRId64" max_pos=%"PRId64"\n",
                 ra[i].vblock_i, chrom_snip_len, chrom_snip, ra[i].chrom_index, ra[i].min_pos, ra[i].max_pos);
    }
}

void random_access_get_ra_info (uint32_t vblock_i, WordIndex *chrom_index, PosType *min_pos, PosType *max_pos)
{
    const RAEntry *ra = random_access_get_first_ra_of_vb (vblock_i, FIRSTENT (RAEntry, z_file->ra_buf), LASTENT (RAEntry, z_file->ra_buf));

    *chrom_index = ra->chrom_index;
    *min_pos     = ra->min_pos;
    *max_pos     = ra->max_pos;
}

void random_access_load_ra_section (SectionType sec_type, Buffer *ra_buf, const char *buf_name, const char *show_index_msg)
{
    SectionListEntry *ra_sl = sections_get_offset_first_section_of_type (sec_type, false);
    zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", sizeof (SectionHeader), sec_type, ra_sl);

    zfile_uncompress_section (evb, evb->z_data.data, ra_buf, buf_name, 0, sec_type);

    ra_buf->len /= sizeof (RAEntry);
    BGEN_random_access (ra_buf);

    if (show_index_msg) {
        random_access_show_index (ra_buf, false, show_index_msg);
        if (exe_type == EXE_GENOCAT) exit_ok; // in genocat --show-index, we only show the index, not the data
    }

    buf_free (&evb->z_data);
}

void random_access_compress (Buffer *ra_buf, SectionType sec_type, const char *msg)
{
    if (msg) random_access_show_index (ra_buf, true, msg);
    
    BGEN_random_access (ra_buf); // make ra_buf into big endian

    ra_buf->len *= sizeof (RAEntry); // change len to count bytes
    zfile_compress_section_data_codec (evb, sec_type, ra_buf, 0,0, CODEC_LZMA); // ra data compresses better with LZMA than BZLIB
    ra_buf->len /= sizeof (RAEntry); // restore
}