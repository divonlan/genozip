// ------------------------------------------------------------------
//   hash.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <math.h>
#include "genozip.h"
#include "buffer.h"
#include "move_to_front.h"
#include "file.h"
#include "vb.h"
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
static uint32_t hash_next_size_up (uint32_t size)
{
    // primary numbers just beneath the powers of 2
    #define NUM_HASH_SIZES 19
    static uint32_t hash_sizes[NUM_HASH_SIZES] = { 509, 1021, 2039, 4039, 8191, 16381, 32749, 65521, 131071, 
                                                   262139, 524287, 1048573, 2097143, 4194301, 8388593, 16777213, 
                                                   33554393, 67108859, 134217689 };

    for (int i=0; i < NUM_HASH_SIZES; i++)
        if (size < hash_sizes[i]) return hash_sizes[i];

    return hash_sizes[NUM_HASH_SIZES-1]; // the maximal size
}

// this is called when the first VB is merging - after its mtf and dict have been moved to zf_ctx.
// Now we need to populate hash. 
static void hash_populate_from_mtf (MtfContext *zf_ctx)
{
    for (uint32_t i=0; i < zf_ctx->mtf.len; i++) {
        MtfNode *node = ENT (MtfNode, &zf_ctx->mtf, i);
        const char *snip = &zf_ctx->dict.data[node->char_index];
        hash_get_entry_for_merge (zf_ctx, snip, node->snip_len, i, NULL /* the snip is not in the hash for sure */); 
    }
}

// This is called when the VB encounters a first entry that's not in the global dictionary (possibly there is no global dict)
// allocation algorithm:
// 1. If we got info on the size of this dict with the previous merged vb - use that size
// 2. If not - use either num_lines for the size, or the smallest size for dicts that are typically small
void hash_alloc_local (VariantBlock *segging_vb, MtfContext *vb_ctx)
{
    // if known from previously merged vb - use those values
    if (vb_ctx->num_new_entries_prev_merged_vb)
        // 3X the expected number of entries to reduce spill-over
        vb_ctx->local_hash_prime = hash_next_size_up (vb_ctx->num_new_entries_prev_merged_vb * 3);

    // if typically small - use minimal hash table
    else if (vb_ctx->dict_id.num == dict_id_vardata_fields[CHROM]  ||
             vb_ctx->dict_id.num == dict_id_vardata_fields[FORMAT] ||
             vb_ctx->dict_id.num == dict_id_vardata_fields[INFO]   ||
             vb_ctx->dict_id.num == dict_id_vardata_fields[REFALT] ||
             vb_ctx->dict_id.num == dict_id_vardata_fields[FILTER] ||
             vb_ctx->dict_id.num == dict_id_INFO_AC ||
             vb_ctx->dict_id.num == dict_id_INFO_AF ||
             vb_ctx->dict_id.num == dict_id_INFO_AN ||
             vb_ctx->dict_id.num == dict_id_INFO_DP)
             
        vb_ctx->local_hash_prime = hash_next_size_up(1);

    // typically big - use large hash table
    else if (vb_ctx->dict_id.num == dict_id_INFO_VQSLOD ||
             vb_ctx->dict_id.num == dict_id_FORMAT_GL ||
             vb_ctx->dict_id.num == dict_id_FORMAT_PL)

        vb_ctx->local_hash_prime = hash_next_size_up(segging_vb->num_lines);

    // if could be big - start with num_lines / 10 (this is an estimated num_lines that is likely inflated)
    else 
        vb_ctx->local_hash_prime = hash_next_size_up(segging_vb->num_lines / 10);

             
    // note: we can't be too generous with the initial allocation because this memory is usually physically allocated
    // to ALL VB structures before any of them merges. Better start smaller for vb_i=1 and let it extend if needed
    buf_alloc (segging_vb, &vb_ctx->local_hash, (vb_ctx->local_hash_prime * 1.2) * sizeof (LocalHashEnt) /* room for expansion */, 1, 
               "mtf_ctx->local_hash", vb_ctx->did_i);
    vb_ctx->local_hash.len = vb_ctx->local_hash_prime;
    memset (vb_ctx->local_hash.data, 0xff, vb_ctx->local_hash_prime * sizeof (LocalHashEnt)); // initialize core table
//printf ("Seg vb_i=%u: local hash: dict=%.8s size=%u\n", segging_vb->variant_block_i, dict_id_printable (vb_ctx->dict_id).id, vb_ctx->local_hash_prime); 
}

// ZIP merge: allocating the global cache for a dictionary, when merging the first VB that encountered it
// an attempt is made to set the cache size according to the expected needs. for this, we use data
// about how many entries were added in the first half of the vb vs the second half, as well as the file size,
// to extrapolate the expected growth
void hash_alloc_global (VariantBlock *merging_vb, MtfContext *zf_ctx, const MtfContext *first_merging_vb_ctx)
{
    // an approximation of the number of VBs we expect with similar data. It is inaccurate in these cases:
    // 1. when the file is not roughly uniform - i.e. if other VBs look differently
    // 2. if there are more concatenated VCF files coming - we don't take them into account
    // 3. if we don't have the file size because it is coming from redirected stdin
    // 4. if we don't have the file size because it is a compressed (gz/bz2) file - we attempt to correct for this
    uint64_t estimated_vcf_file_size = vcf_file->disk_size;
    if (vcf_file->type == VCF_BZ2) estimated_vcf_file_size *= 15; // compression ratio of bzip2 of a "typical" vcf file
    if (vcf_file->type == VCF_GZ)  estimated_vcf_file_size *= 8;  // compression ratio of gzip of a "typical" vcf file

    uint64_t expected_num_vbs = MAX (1, estimated_vcf_file_size / (uint64_t)merging_vb->vcf_data.len);

    double n1 = first_merging_vb_ctx->mtf_len_at_half;
    double n2 = (int)first_merging_vb_ctx->ol_mtf.len - (int)first_merging_vb_ctx->mtf_len_at_half;

    // method 1 - we got it with trial and error - more accurate than method 2 with large dictionaries,
    // but tends to under-estimate with small ones
    double reduction_factor_per_vb = MIN (1, ((n1 > 0) ? (n2 / n1) * (n2 / n1) : 0.000001)); // a number (0,1] )

    double method_1_total = n1+n2; 
    double growth_by = n1+n2; 

    for (unsigned vb_i=1; vb_i < expected_num_vbs && growth_by > 500; vb_i++) {
        growth_by *= reduction_factor_per_vb;
        reduction_factor_per_vb = pow (reduction_factor_per_vb, 0.818); // magic number
        method_1_total += growth_by; 
    }

    // method 2 - somewhat more scientific - accurate results for smaller dictionaries but underestimates
    // for larger ones 

    // some math: definitions:
    // simplifying assumptions: 
    // 1. we assume each value is random across a uniform space (i.e. values are spread across the file without structure)
    // 2. we assume that when the second half of the lines were checking if entries exist, they only checked those
    //    entries from the first half. this is of course false, but it is an ok approximation in most cases
    // 3. we assume that VBs and the number of lines in each VB is similar, so we will take the first VB data
    //    to be representative of all VBs to come
    //
    // n1, n2 - the number of entries, as counted in the vb, in each one of approximately half of the lines of the vb
    // L1, L2 - the number of lines acorss which n1 and n2 were counted
    // N - the unknown size of the uniform space
    // calculation regarding the H2 (second half of the vb), when n1 out of the possible (unkown) N values are already in ctx:
    // approximation of EXPECTED number of rows to be NOT in ctx in H2: (L2 - (((n1+2)/N) * L2). 
    // We equate this is ACTUAL = n2 and extract N:
    // N = (n1*L2) / (L2-n2)
    // For the next VB with L3 lines, and np previous entires (starting with total_so_far = n1+n2), the expected 
    //     % of hits = (1 - (total_so_far / N)) * L3

    // note: we do the calculations in "double" for more accurate divisions, and to support large intermediate numbers
    double line_factor = dict_id_is_format_subfield (zf_ctx->dict_id) ? global_num_samples : 1; // for FORMAT subfields - the number of "slots" for items is lines x samples
    double L2 = ((int)merging_vb->num_lines - (int)first_merging_vb_ctx->num_lines_at_half) * line_factor;    
    double N = (L2 != n2) ? (((n1+n2)*L2) / (L2-n2)) : 0;
    double method_2_total = n1+n2; 

    if (N > 0) {
        // now, we iterate to calculate how many entries is each VB expected to add, and what do we
        // expect the total number of entries to be
        for (unsigned vb_i=1; vb_i < expected_num_vbs; vb_i++) 
            method_2_total += (MIN (N, 1 - MIN (1, method_2_total / N)) * merging_vb->num_lines * line_factor);  // merging_vb->num_lines is an approximation of next VBs num_lines
    }

    // our final estimate is the maximum of method 1 and method 2

    // we allocate 5X total_so_far to leave the hash table a bit less crowded - less spill over
    // (this is the global table - only 1 copy of it shared by all threads (unless there are reallocs) -
    // so we can be a bit more generous with space). reallocs are also less likely with larger spaces
    // (less entries spill over), so it might even not add up that much more memory.
    zf_ctx->global_hash_prime = hash_next_size_up ((uint32_t)MAX(method_1_total, method_2_total) * 5);

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

        ASSERT (g_hashent->next < zf_ctx->global_hash.len, "Error in hash_get_entry_for_merge: g_hashent->next=%d out of range, hash.len=%u",  
                g_hashent->next, zf_ctx->global_hash.len);

        hashent_i = g_hashent->next;
        g_hashent = ENT(GlobalHashEnt, &zf_ctx->global_hash, hashent_i);

        // case: snip is not in core hash table and also no other snip occupies the slot (mtf_i==NIL happens only in the core table)
        if (g_hashent->mtf_i == NIL) { // unoccupied space in core hash table
            // thread safety: VB threads with merge_num < ours, might be segmenting right now, and have this global hash overlayed 
            // and accessing it. we make sure that the setting of g_hashent->merge_num is atomic and the other threads will at all times either
            // see NIL or merge_num - both of which indicate the this entry is effectively NIL as they have an older merge_num
// clang: check __c11_atomic_store: https://clang.llvm.org/docs/LanguageExtensions.html#langext-c11-atomic
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

    g_hashent = ENT (GlobalHashEnt, &zf_ctx->global_hash, hashent_i); // might have changed after realloc

    // thread safetey:  VB threads with merge_num < ours, might be segmenting right now, and have this global hash overlayed 
    // and accessing it. We make sure to first prepare the new entry including the merge_num which will prohibit old
    // VBs from using it, before we atomically set the "next"
    uint32_t next = zf_ctx->global_hash.len++;
    GlobalHashEnt *new_hashent = ENT (GlobalHashEnt, &zf_ctx->global_hash, next);
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
int32_t hash_get_entry_for_seg (VariantBlock *segging_vb, MtfContext *vb_ctx,
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

        g_hashent = ENT(GlobalHashEnt, &vb_ctx->global_hash, g_hashent->next);

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
                "Error: l_hashent->next=%d out of range, local_hash.len=%u", l_hashent->next, vb_ctx->local_hash.len);
        l_hashent_i = l_hashent->next;
        l_hashent = ENT (LocalHashEnt, &vb_ctx->local_hash, l_hashent_i);

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

    l_hashent = ENT (LocalHashEnt, &vb_ctx->local_hash, l_hashent_i);  // might have changed after realloc
    l_hashent->next = vb_ctx->local_hash.len++;

    LocalHashEnt *new_l_hashent = ENT (LocalHashEnt, &vb_ctx->local_hash, l_hashent->next);
    new_l_hashent->next = NO_NEXT;
    new_l_hashent->mtf_i = new_mtf_i_if_no_old_one;

    if (node) *node = NULL;
    return NIL;
}
