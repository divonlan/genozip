// ------------------------------------------------------------------
//   hash.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <math.h>
#include "genozip.h"
#include "buffer.h"
#include "move_to_front.h"
#include "file.h"
#include "vblock.h"
#include "hash.h"
#include "dict_id.h"

#define NO_NEXT 0xffffffff
typedef struct {        
    int32_t mtf_i;            // index into MtfContext.ol_mtf (if < ol_mtf.len) or MtfContext.mtf or NIL
    uint32_t next;            // linked list - index into MtfContext.global/local_hash or NIL
                              //               local_hash indeces started at LOCAL_HASH_OFFSET
} LocalHashEnt;

#pragma pack(push, hash, 4)
typedef struct {        
    int32_t mtf_i;            // index into MtfContext.ol_mtf (if < ol_mtf.len) or MtfContext.mtf or NIL
    uint32_t next;            // linked list - index into MtfContext.global/local_hash or NIL
    int32_t merge_num;        // the merge_num in which the "mtf_i" field was set. when this global hash is overlayed 
                              // to a vb_ctx, that vb_ctx is permitted use the mtf_i value if this merge_num is <= vb_ctx->merge_num,
                              // otherwise, it should treat it as NIL.
} GlobalHashEnt;
#pragma pack (pop)

// get the size of the hash table - a primary number between roughly 0.5K and 8M that is close to size, or a bit bigger
static uint32_t hash_next_size_up (uint64_t size)
{
    // primary numbers just beneath the powers of 2
    static uint32_t hash_sizes[] = { 509, 1021, 2039, 4039, 8191, 16381, 32749, 65521, 131071, 
                                     262139, 524287, 1048573, 2097143, 4194301, 8388593, 16777213, 33554393 };
    #define NUM_HASH_SIZES (sizeof(hash_sizes) / sizeof(hash_sizes[0]))

    for (int i=0; i < NUM_HASH_SIZES; i++)
        if (size < (uint64_t)hash_sizes[i]) return hash_sizes[i];

    return hash_sizes[NUM_HASH_SIZES-1]; // the maximal size
}

// this is called when the first VB is merging - after its mtf and dict have been moved to zf_ctx.
// Now we need to populate hash. 
static void hash_populate_from_mtf (MtfContext *zf_ctx)
{
    for (uint32_t i=0; i < zf_ctx->mtf.len; i++) {
        MtfNode *node = ENT (MtfNode, zf_ctx->mtf, i);
        const char *snip = &zf_ctx->dict.data[node->char_index];
        hash_get_entry_for_merge (zf_ctx, snip, node->snip_len, i, NULL /* the snip is not in the hash for sure */); 
    }
}

// This is called when the VB encounters a first entry that's not in the global dictionary (possibly there is no global dict)
// allocation algorithm:
// 1. If we got info on the size of this dict with the previous merged vb - use that size
// 2. If not - use either num_lines for the size, or the smallest size for dicts that are typically small
void hash_alloc_local (VBlock *segging_vb, MtfContext *vb_ctx)
{
    vb_ctx->local_hash_prime = 0; // initialize

    // if known from previously merged vb - use those values
    if (vb_ctx->num_new_entries_prev_merged_vb)
        // 3X the expected number of entries to reduce spill-over
        vb_ctx->local_hash_prime = hash_next_size_up (vb_ctx->num_new_entries_prev_merged_vb * 3);

    else if (z_file->data_type == DATA_TYPE_VCF) {

        // if typically small - use minimal hash table
        if (vb_ctx->dict_id.num == dict_id_vcf_fields[VCF_CHROM]  ||
            vb_ctx->dict_id.num == dict_id_vcf_fields[VCF_FORMAT] ||
            vb_ctx->dict_id.num == dict_id_vcf_fields[VCF_INFO]   ||
            vb_ctx->dict_id.num == dict_id_vcf_fields[VCF_REFALT] ||
            vb_ctx->dict_id.num == dict_id_vcf_fields[VCF_FILTER] ||
            vb_ctx->dict_id.num == dict_id_INFO_AC ||
            vb_ctx->dict_id.num == dict_id_INFO_AF ||
            vb_ctx->dict_id.num == dict_id_INFO_AN ||
            vb_ctx->dict_id.num == dict_id_INFO_DP)
            
            vb_ctx->local_hash_prime = hash_next_size_up(1);

        // typically big - use large hash table
        else 
        if (vb_ctx->dict_id.num == dict_id_INFO_VQSLOD ||
            vb_ctx->dict_id.num == dict_id_FORMAT_GL   ||
            vb_ctx->dict_id.num == dict_id_FORMAT_PL)

            vb_ctx->local_hash_prime = hash_next_size_up(segging_vb->num_lines);
    }

    else if (z_file->data_type == DATA_TYPE_SAM) {
        // typically small - use minimal hash table ~ 500
        if (vb_ctx->dict_id.num == dict_id_sam_fields[SAM_FLAG]      ||
            vb_ctx->dict_id.num == dict_id_sam_fields[SAM_MAPQ]      ||
            vb_ctx->dict_id.num == dict_id_sam_fields[SAM_QNAME]     ||
            vb_ctx->dict_id.num == dict_id_sam_fields[SAM_OPTIONAL]  ||

            // standard tags, see here: https://samtools.github.io/hts-specs/SAMtags.pdf
            vb_ctx->dict_id.num == dict_id_OPTION_AM ||
            vb_ctx->dict_id.num == dict_id_OPTION_AS ||
            vb_ctx->dict_id.num == dict_id_OPTION_CM ||
            vb_ctx->dict_id.num == dict_id_OPTION_LB ||
            vb_ctx->dict_id.num == dict_id_OPTION_FI ||
            vb_ctx->dict_id.num == dict_id_OPTION_H0 ||
            vb_ctx->dict_id.num == dict_id_OPTION_H1 ||
            vb_ctx->dict_id.num == dict_id_OPTION_H2 ||
            vb_ctx->dict_id.num == dict_id_OPTION_MQ ||
            vb_ctx->dict_id.num == dict_id_OPTION_NH ||
            vb_ctx->dict_id.num == dict_id_OPTION_NM ||
            vb_ctx->dict_id.num == dict_id_OPTION_OC ||
            vb_ctx->dict_id.num == dict_id_OPTION_PG ||
            vb_ctx->dict_id.num == dict_id_OPTION_PQ ||
            vb_ctx->dict_id.num == dict_id_OPTION_PU ||
            vb_ctx->dict_id.num == dict_id_OPTION_RG ||
            vb_ctx->dict_id.num == dict_id_OPTION_SA ||
            vb_ctx->dict_id.num == dict_id_OPTION_SM ||
            vb_ctx->dict_id.num == dict_id_OPTION_TC ||
            vb_ctx->dict_id.num == dict_id_OPTION_UQ ||
            
            // bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
            vb_ctx->dict_id.num == dict_id_OPTION_X0 ||
            vb_ctx->dict_id.num == dict_id_OPTION_X1 ||
            vb_ctx->dict_id.num == dict_id_OPTION_XA ||
            vb_ctx->dict_id.num == dict_id_OPTION_XN ||
            vb_ctx->dict_id.num == dict_id_OPTION_XM ||
            vb_ctx->dict_id.num == dict_id_OPTION_XO ||
            vb_ctx->dict_id.num == dict_id_OPTION_XG ||
            vb_ctx->dict_id.num == dict_id_OPTION_XS ||
            vb_ctx->dict_id.num == dict_id_OPTION_XE ||
            
            vb_ctx->dict_id.num == dict_id_OPTION_STRAND)
            
            vb_ctx->local_hash_prime = hash_next_size_up(500);

        // typically smallish - use hash table ~ 2000
        else 
        if (vb_ctx->dict_id.num == dict_id_sam_fields[SAM_RNAME]  ||
            vb_ctx->dict_id.num == dict_id_OPTION_CC)

            vb_ctx->local_hash_prime = hash_next_size_up(2000);

        // typically medium - use hash table ~ 50000
        else 
        if (vb_ctx->dict_id.num == dict_id_sam_fields[SAM_CIGAR]  ||
            vb_ctx->dict_id.num == dict_id_OPTION_CG ||
            vb_ctx->dict_id.num == dict_id_OPTION_MC)

            vb_ctx->local_hash_prime = hash_next_size_up(50000);
    }

    // default: it could be big - start with num_lines / 10 (this is an estimated num_lines that is likely inflated)
    if (!vb_ctx->local_hash_prime) 
        vb_ctx->local_hash_prime = hash_next_size_up(segging_vb->num_lines / 10);

    // note: we can't be too generous with the initial allocation because this memory is usually physically allocated
    // to ALL VB structures before any of them merges. Better start smaller for vb_i=1 and let it extend if needed
    buf_alloc (segging_vb, &vb_ctx->local_hash, (vb_ctx->local_hash_prime * 1.2) * sizeof (LocalHashEnt) /* room for expansion */, 1, 
               "mtf_ctx->local_hash", vb_ctx->did_i);
    vb_ctx->local_hash.len = vb_ctx->local_hash_prime;
    memset (vb_ctx->local_hash.data, 0xff, vb_ctx->local_hash_prime * sizeof (LocalHashEnt)); // initialize core table
//printf ("Seg vb_i=%u: local hash: dict=%.8s size=%u\n", segging_vb->vblock_i, dict_id_printable (vb_ctx->dict_id).id, vb_ctx->local_hash_prime); 
}

// ZIP merge: allocating the global cache for a dictionary, when merging the first VB that encountered it
// an attempt is made to set the cache size according to the expected needs. for this, we use data
// about how many entries were added in the first half of the vb vs the second half, as well as the file size,
// to extrapolate the expected growth
// it is very important to get this as accurate as possible: merge is our bottleneck for core-count scalability
// as the merge is protected by a per-dictionary mutex. if the global hash table size is too small, search time
// goes up (because of the need to traverse linked lists) during the bottleneck time. Coversely, if the hash
// table size is too big, it both consumes a lot memory, as well as slows down the search time as the dictionary
// is less likely to fit into the CPU memory caches
void hash_alloc_global (VBlock *merging_vb, MtfContext *zf_ctx, const MtfContext *first_merging_vb_ctx)
{
    // note on txt_data_size_single: if its a physical plain VCF file - this is the file size. 
    // if not - its an estimate done after the first VB by txtfile_estimate_txt_data_size
    double estimated_num_vbs = MAX (1, (double)txt_file->txt_data_size_single / (double)merging_vb->txt_data.len);
    double estimated_num_lines = estimated_num_vbs * (double)merging_vb->num_lines;

    double n1 = first_merging_vb_ctx->mtf_len_at_half;
    double n2 = (int)first_merging_vb_ctx->ol_mtf.len - (int)first_merging_vb_ctx->mtf_len_at_half;

    static struct { uint64_t vbs, factor; } growth_plan[] =
    { 
        { 0, 1 },
        { 1, 2 },
        { 5, 4 },
        { 9, 8 },
        { 17, 16 },
        { 44, 23 },
        { 114, 32 },
        { 295, 45 },
        { 765, 64 },
        { 1981, 91 },
        { 5132, 128 },
        { 13291, 181 },
        { 34423, 256 },
        { 89155, 362 },
        { 230912, 512 },
        { 598063, 724 },
        { 1548983, 1024 },
    };

    double estimated_entries=0;
    double n_ratio = n2 ? n1/n2 : 0;    unsigned max_growth_plan=0;

    if (n2 == 0) 
        estimated_entries = zf_ctx->mtf.len;

    // To do - take into account optional fields partial appearance (e.g. if a field appears 1% of the time...)

    else if (n_ratio > 0.8 && n_ratio < 1.2)  // looks like almost linear growth
        estimated_entries = estimated_num_lines * ((n1+n2 )/ merging_vb->num_lines) * 0.75;
    
    else {
        if      (n_ratio > 2.5 || !n1) max_growth_plan = 2;
        else if (n_ratio > 2.1) max_growth_plan = 3;
        else if (n_ratio > 1.8) max_growth_plan = 4;
        else if (n_ratio > 1.5) max_growth_plan = 5;
        else if (n_ratio > 1.2) max_growth_plan = 7;
        else if (n_ratio > 1.1) max_growth_plan = 9;
        else                    max_growth_plan = sizeof(growth_plan) / sizeof(growth_plan[0]);
        
        for (int i=max_growth_plan; i >= 0 ; i--)
            if (estimated_num_vbs > growth_plan[i].vbs) {
                estimated_entries = zf_ctx->mtf.len * growth_plan[i].factor;
                break;
            }
    }

    if (!estimated_entries) estimated_entries = 100000000000; // very very big

    zf_ctx->global_hash_prime = hash_next_size_up (estimated_entries * 5);
    //printf ("dict=%.8s n1=%2.2lf n2=%2.2lf n1/n2=%2.2lf max_growth_plan=%u vbs=%u num_lines=%u zf_ctx->mtf.len=%u entries=%2.2lf hashsize=%u\n", 
    //        dict_id_printable(zf_ctx->dict_id).id, n1, n2, n_ratio, max_growth_plan, (unsigned)estimated_num_vbs, (unsigned)estimated_num_lines, (uint32_t)zf_ctx->mtf.len, estimated_entries, zf_ctx->global_hash_prime); 

    buf_alloc (evb, &zf_ctx->global_hash, sizeof(GlobalHashEnt) * zf_ctx->global_hash_prime * 1.5, 1,  // 1.5 - leave some room for extensions
               "z_file->mtf_ctx->global_hash", zf_ctx->did_i);
    buf_set_overlayable (&zf_ctx->global_hash);

    zf_ctx->global_hash.len = zf_ctx->global_hash_prime; // global_hash.len can get longer over time as extension links are added

    // we set all entries to {NO_NEXT, NIL, NIL} == {0xffffffff x 3} (note: GlobalHashEnt is packed)
    memset (zf_ctx->global_hash.data, 0xff, sizeof(GlobalHashEnt) * zf_ctx->global_hash.len);

    hash_populate_from_mtf (zf_ctx);
}

// tested hash table sizes up to 5M. turns out smaller tables (up to a point) are faster, despite having longer
// average linked lists. probably bc the CPU can store the entire hash and mtf arrays in L1 or L2
// memory cache during segmentation
static inline uint32_t hash_do (uint32_t hash_len, const char *snip, unsigned snip_len)
{
    if (!hash_len) return NO_NEXT; // hash table does not exist
    
    // spread the snip throughout the 64bit word before taking a mod - to ensure about-even distribution 
    // across the hash table
    uint64_t result=0;
    for (unsigned i=0; i < snip_len; i++) 
        result = ((result << 23) | (result >> 41)) ^ (uint64_t)((uint8_t)snip[i]);

    return (uint32_t)(result % hash_len);
}

// gets the mtf_i if the snip is already in the hash table, or puts a new one in the hash table in not
int32_t hash_get_entry_for_merge (MtfContext *zf_ctx, const char *snip, unsigned snip_len, 
                                  int32_t new_mtf_i_if_no_old_one,
                                  MtfNode **node)        // out - node if node is found, NULL if not
{
    GlobalHashEnt g_head, *g_hashent = &g_head;
    g_hashent->next = hash_do (zf_ctx->global_hash_prime, snip, snip_len); // entry in hash table determined by hash function on snip
    int32_t hashent_i = NO_NEXT; // squash compiler warning (note: global hash table is always allocated in the first merge)

    while (g_hashent->next != NO_NEXT) {

        ASSERT (g_hashent->next < zf_ctx->global_hash.len, "Error in hash_get_entry_for_merge: g_hashent->next=%d out of range, hash.len=%"PRIu64,  
                g_hashent->next, zf_ctx->global_hash.len);

        hashent_i = g_hashent->next;
        g_hashent = ENT(GlobalHashEnt, zf_ctx->global_hash, hashent_i);

        // case: snip is not in core hash table and also no other snip occupies the slot (mtf_i==NIL happens only in the core table)
        if (g_hashent->mtf_i == NIL) { // unoccupied space in core hash table
            // thread safety: VB threads with merge_num < ours, might be segmenting right now, and have this global hash overlayed 
            // and accessing it. we make sure that the setting of g_hashent->merge_num is atomic and the other threads will at all times either
            // see NIL or merge_num - both of which indicate the this entry is effectively NIL as they have an older merge_num

            g_hashent->next = NO_NEXT;
            g_hashent->mtf_i = new_mtf_i_if_no_old_one;
            __atomic_store_n (&g_hashent->merge_num, zf_ctx->merge_num, __ATOMIC_RELAXED); // stamp our merge_num as the ones that set the mtf_i
            if (node) *node = NULL;
            return NIL;
        }

        if (node) {  // if node=NULL, caller is telling us it is not in MTF for sure
            const char *snip_in_dict;
            uint32_t snip_len_in_dict;

            *node = mtf_node (zf_ctx, g_hashent->mtf_i, &snip_in_dict, &snip_len_in_dict);
        
            // case: snip is in the hash table 
            if (snip_len == snip_len_in_dict && !memcmp (snip, snip_in_dict, snip_len)) 
                return g_hashent->mtf_i;
        }
    }

    // case: not found in hash table, and we are required to provide a new hash entry on the linked list
    buf_alloc (evb, &zf_ctx->global_hash, sizeof (GlobalHashEnt) * (1 + zf_ctx->global_hash.len), 2, 
               "z_file->mtf_ctx->global_hash", zf_ctx->did_i);

    g_hashent = ENT (GlobalHashEnt, zf_ctx->global_hash, hashent_i); // might have changed after realloc

    // thread safetey:  VB threads with merge_num < ours, might be segmenting right now, and have this global hash overlayed 
    // and accessing it. We make sure to first prepare the new entry including the merge_num which will prohibit old
    // VBs from using it, before we atomically set the "next"
    uint32_t next = zf_ctx->global_hash.len++;
    GlobalHashEnt *new_hashent = ENT (GlobalHashEnt, zf_ctx->global_hash, next);
    new_hashent->merge_num = zf_ctx->merge_num; // stamp our merge_num as the ones that set the mtf_i
    new_hashent->mtf_i     = new_mtf_i_if_no_old_one;
    new_hashent->next      = NO_NEXT;
    
    // now, with the new g_hashent set, we can atomically update the "next"
    __atomic_store_n (&g_hashent->next, next, __ATOMIC_RELAXED);

    if (node) *node = NULL; // we don't have an old node
    return NIL;
}

// gets the mtf_i if the snip is already in the hash table, or puts a new one in the hash table in not
// 1. if its in the global hash table, with merge_num lower or equal to ours - i.e. added by an earler thread - we take it 
// 2. if its in the local hash table - i.e. added by us (this vb) earlier - we take it
// 3. if not found - we add it to the local hash table
int32_t hash_get_entry_for_seg (VBlock *segging_vb, MtfContext *vb_ctx,
                                const char *snip, unsigned snip_len, 
                                int32_t new_mtf_i_if_no_old_one,
                                MtfNode **node)        // out - node if node is found, NULL if not
{
    // first, search for the snip in the global table
    GlobalHashEnt g_head, *g_hashent = &g_head;
    g_hashent->next = hash_do (vb_ctx->global_hash_prime, snip, snip_len); // entry in hash table determined by hash function on snip

    for (unsigned depth=0; ; depth++) {
        
        // four cases for "next":
        // 1. NO_NEXT - the linked list is terminated and we have not found our snip
        // 2. next points to a valid index with our overlaid global_cache - we follow it
        // two cases where another thread sets next during merge with a higher merge number. two case:
        // 3. "next" is within the length cloned by us - we will follow it, but arrive at a node 
        //    which its merge_num prohibits us from consuming
        // 4. "next" falls outside the length cloned by us - entries added to hash, but without reallocing
        //    after we cloned. in this case, we will detected that "next" is out of range.
        uint32_t next = __atomic_load_n (&g_hashent->next, __ATOMIC_RELAXED);
        
        if (next == NO_NEXT || /* case 1 */ next >= vb_ctx->global_hash.len) // case 4
            break;

        g_hashent = ENT(GlobalHashEnt, vb_ctx->global_hash, next);

        // case: snip is not in core hash table (at least it wasn't there when we cloned and set our maximum merge_num we accept)
        uint32_t merge_num = __atomic_load_n (&g_hashent->merge_num, __ATOMIC_RELAXED);
        if (g_hashent->mtf_i == NIL || merge_num > vb_ctx->merge_num) break; // case 3

        const char *snip_in_dict;
        uint32_t snip_len_in_dict;
        *node = mtf_node (vb_ctx, g_hashent->mtf_i, &snip_in_dict, &snip_len_in_dict);

        // case: snip is in the global hash table - we're done
        if (snip_len == snip_len_in_dict && !memcmp (snip, snip_in_dict, snip_len)) 
            return g_hashent->mtf_i; // case 2
    }

    // snip was not found in the global hash table (as it was at the time we cloned), we now search
    // in our local hash table - and if not found there - we will add it

    // allocate hash table, if not already allocated, based on experience of previous VBs, or pre-set default if there isn't any
    if (!buf_is_allocated (&vb_ctx->local_hash)) 
        hash_alloc_local (segging_vb, vb_ctx);

    LocalHashEnt l_head, *l_hashent = &l_head;
    l_hashent->next = hash_do (vb_ctx->local_hash_prime, snip, snip_len); // entry in hash table determined by hash function on snip
    int32_t l_hashent_i = NO_NEXT; // initialize to squash compiler warning

    while (l_hashent->next != NO_NEXT) {

        ASSERT (l_hashent->next < vb_ctx->local_hash.len, 
                "Error: l_hashent->next=%d out of range, local_hash.len=%u", l_hashent->next, (uint32_t)vb_ctx->local_hash.len);
        l_hashent_i = l_hashent->next;
        l_hashent = ENT (LocalHashEnt, vb_ctx->local_hash, l_hashent_i);

        // case: snip is not in hash table and also no other snip occupies the slot (mtf_i==NIL happens only in the core table)
        if (l_hashent->mtf_i == NIL) { // unoccupied space in core hash table
            l_hashent->next = NO_NEXT;
            l_hashent->mtf_i = new_mtf_i_if_no_old_one;
            if (node) *node = NULL;
            return NIL;
        }

        if (node) { // if the caller doesn't provide "node", he is telling us that with certainly the snip is not in the hash table
            const char *snip_in_dict;
            uint32_t snip_len_in_dict;
            *node = mtf_node (vb_ctx, l_hashent->mtf_i, &snip_in_dict, &snip_len_in_dict);

            // case: snip is in the hash table - we're done
            if (snip_len == snip_len_in_dict && !memcmp (snip, snip_in_dict, snip_len)) 
                return l_hashent->mtf_i;
        }
    }

    // case: not found in hash table, and we are required to provide a new hash entry on the linked list
    buf_alloc (segging_vb, &vb_ctx->local_hash, sizeof (LocalHashEnt) * (1 + vb_ctx->local_hash.len), // realloc if needed
                1.5, "mtf_ctx->local_hash", vb_ctx->did_i);

    l_hashent = ENT (LocalHashEnt, vb_ctx->local_hash, l_hashent_i);  // might have changed after realloc
    l_hashent->next = vb_ctx->local_hash.len++;

    LocalHashEnt *new_l_hashent = ENT (LocalHashEnt, vb_ctx->local_hash, l_hashent->next);
    new_l_hashent->next = NO_NEXT;
    new_l_hashent->mtf_i = new_mtf_i_if_no_old_one;

    if (node) *node = NULL;
    return NIL;
}
