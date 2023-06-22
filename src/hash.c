// ------------------------------------------------------------------
//   hash.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "genozip.h"
#include "buffer.h"
#include "context.h"
#include "file.h"
#include "vblock.h"
#include "hash.h"
#include "strings.h"
#include "dict_id.h"
#include "segconf.h"

typedef struct {        
    WordIndex node_index;     // index into Context.ston_nodes (if < ston_nodes.len) or Context.nodes or NODE_INDEX_NONE
    uint32_t next;            // linked list - index into Context.global/local_ents or NO_NEXT
                              //               local_ents indices start at LOCAL_HASH_OFFSET
} LocalHashEnt;

#pragma pack(4)
typedef struct {        
    WordIndex node_index;     // index into Context.nodes or NODE_INDEX_NONE
    uint32_t next;            // linked list - index into Context.global/local_ents or NO_NEXT
    int32_t merge_num;        // the merge_num in which the "node_index" field was set. when this global hash is overlayed 
                              // to a vctx, that vctx is permitted use the node_index value if this merge_num is <= vctx->merge_num,
                              // otherwise, it should treat it as NODE_INDEX_NONE.
} GlobalHashEnt;
#pragma pack()

// get the size of the hash table - a primary number between roughly 0.5K and 8M that is close to size, or a bit bigger
uint32_t hash_next_size_up (uint64_t size, bool allow_huge)
{
    // if user set a low --vblock, use this as the limit for the hash table size too, but don't restrict tighter than 16MB
    size = MIN_(size, MAX_(16000000, segconf.vb_size));

    // primary numbers just beneath the powers of 2^0.5 (and 2^0.25 for the larger numbers)
    // minimum ~64K to prevent horrible miscalculations in edge cases that result in dramatic slow down
    static uint32_t hash_sizes[] = { 65521, 92681, 131071, 
                                     185363, 262139, 370723, 524287, 741431, 1048573, 1482907, 2097143, 2965819, 4194301, 5931641, 8388593, 
                                     11863279, 16777213, 19951579, 23726561, 28215799, 33554393, 39903161, 47453111, 56431601, 67108859,
                                     94906265, 134217757, 189812533, 268435459, 379625083, 536870923 };
    #define NUM_HASH_REGULAR 25
    #define NUM_HASH_HUGE ARRAY_LEN(hash_sizes)
    #define NUM_HASH_SIZES (allow_huge ? NUM_HASH_HUGE : NUM_HASH_REGULAR)

    for (int i=0; i < NUM_HASH_SIZES; i++)
        if (size < (uint64_t)hash_sizes[i]) return hash_sizes[i];

    return hash_sizes[NUM_HASH_SIZES-1]; // the maximal size
}

// this is called when allocing memory for global hash - it copies pre-exiting nodes data, for example, when copying in a
// reference contig dictionary  
static void hash_populate_from_nodes (ContextP zctx)
{
    uint32_t len = zctx->nodes.len32;
    zctx->nodes.len32 = 0; // hash_global_get_entry will increment it back to its original value

    for (uint32_t i=0; i < len; i++) {
        CtxNode *node = B(CtxNode, zctx->nodes, i);
        rom snip = Bc (zctx->dict, node->char_index);
        hash_global_get_entry (zctx, snip, node->snip_len, HASH_NEW_OK_NOT_SINGLETON, NULL /* the snip is not in the hash for sure */); 
    }
}

// This is called when the VB encounters a first snip that's not in the ol_dict 
// allocation algorithm:
// 1. If we got info on the size of this dict with the previous merged vb - use that size
// 2. If not - use either num_lines for the size, or the smallest size for dicts that are typically small
static void hash_alloc_local (VBlockP segging_vb, ContextP vctx)
{
    vctx->local_hash_prime = 0; // initialize

    // if known from previously merged vb - use those values
    if (vctx->num_new_entries_prev_merged_vb)
        // 3X the expected number of entries to reduce hash contention
        vctx->local_hash_prime = hash_next_size_up (vctx->num_new_entries_prev_merged_vb * 3, false);

    // if known to small, use hash table of ~ 64K
    else if (DT_(segging_vb, seg_is_small)(segging_vb, vctx->dict_id))
        vctx->local_hash_prime = hash_next_size_up(1, false);
    
    // default: it could be big - start with num_lines / 10 (this is an estimated num_lines that is likely inflated)
    else
        vctx->local_hash_prime = hash_next_size_up (segging_vb->lines.len32 / 10, false);

    // note: we can't be too generous with the initial allocation because this memory is usually physically allocated
    // to ALL VB structures before any of them merges. Better start smaller for vb_i=1 and let it extend if needed
    buf_alloc_exact_255 (segging_vb, vctx->local_hash, vctx->local_hash_prime, uint32_t, "contexts->local_hash");

    buf_alloc_exact_255 (segging_vb, vctx->local_ents, vctx->local_hash_prime / 3, LocalHashEnt, "contexts->local_ents");
    vctx->local_ents.len = 0;
}

// ZIP merge: allocating the global hash for a context, when merging the first VB that encountered it
// an attempt is made to set the hash size according to the expected needs. for this, we use data
// about how many entries were added in the first half of the vb vs the second half, as well as the file size,
// to extrapolate the expected growth
// it is very important to get this as accurate as possible: merge is our bottleneck for core-count scalability
// as the merge is protected by a per-dictionary mutex. if the global hash table size is too small, search time
// goes up (because of the need to traverse linked lists) during the bottleneck time. Coversely, if the hash
// table size is too big, it both consumes a lot memory, as well as slows down the search time as the dictionary
// is less likely to fit into the CPU memory caches
uint32_t hash_get_estimated_entries (VBlockP vb, ContextP zctx, ConstContextP vctx)
{
    double effective_num_vbs   = 0; 
    double estimated_num_vbs   = MAX_(1, (double)txtfile_get_seggable_size() / (double)vb->txt_data.len);
    double estimated_num_lines = estimated_num_vbs * (double)vb->lines.len;

    if (flag.show_hash && vctx->did_i==0) 
        iprintf ("\n\nOutput of --show-hash:\n"
                 "est_vbs=%u vb_1_num_lines=%s est_total_lines=%s\n", 
                 (unsigned)ceil(estimated_num_vbs), str_int_commas (vb->lines.len).s, str_int_commas ((uint64_t)estimated_num_lines).s);

    // if known to small, use hash table of ~ 64K
    if (DT_(vb, seg_is_small) (vb, zctx->dict_id)) {

        if (flag.show_hash)
            iprintf ("dict=%s : known to be small. hashsize=%s\n", 
                     vctx->tag_name, str_int_commas (hash_next_size_up (1, false)).s); 

        return 1; // will yield smallest hash table - around 64K
    }

    // for growth purposes, we discard the first 1/3 of VB1, and compare the growth of the 2nd vs the 3rd 1/3. 
    // this is because many fields display the charactistics of having a relatively small number of high frequency entries
    // which mostly show up in n1, and then a long tail of low frequency entries. The comparison of n2 to n3,
    // assuming that the new snips first introduced in them are mostly low frequency ones, will give us a more accurate predication
    // of the gradient appearance of new snips
    double n1 = vctx->nodes_len_at_1_3;
    double n2 = vctx->nodes_len_at_2_3 ? ((double)vctx->nodes_len_at_2_3 - (double)vctx->nodes_len_at_1_3) : 0;
    double n3 = (double)vctx->nodes.len - n1 - n2;

    double n1_lines      = (double)vb->num_lines_at_1_3;
    double n2_lines      = (double)vb->num_lines_at_2_3 - n1_lines;
    double n2_n3_lines   = (double)vb->lines.len - n1_lines;
    double n3_lines      = (double)vb->lines.len - n1_lines - n2_lines;

    double n1_density    = n1_lines    ? (n1 / n1_lines)         : 0; // might be more than 1 if multiple per line, eg VCF sample subfields
    double n2_density    = n2_lines    ? (n2 / n2_lines)         : 0; // might be more than 1 if multiple per line, eg VCF sample subfields
    double n2_n3_density = n2_n3_lines ? ((n2+n3) / n2_n3_lines) : 0;  
    double n3_density    = n3_lines    ? (n3 / n3_lines)         : 0;

    double n2n3_density_ratio = n3_density ? n2_density/n3_density : 0;    

    // my secret formula for arriving at these numbers: https://docs.google.com/spreadsheets/d/1UijOuPgquZ71kEdB7kUf1OYUMs-Xhu2gAOOHgy6QCIQ/edit#gid=0
    static struct { 
    uint64_t   vbs,     factor; } growth_plan[] = { 
    /* 0  */ { 0,       1       },
    /* 1  */ { 1,       2       },
    /* 2  */ { 5,       4       },
    /* 3  */ { 9,       8       },
    /* 4  */ { 17,      16      },
    /* 5  */ { 44,      29      },
    /* 6  */ { 114,     49      },
    /* 7  */ { 295,     78      },
    /* 8  */ { 765,     117     },
    /* 9  */ { 1981,    164     },
    /* 10 */ { 5132,    230     },
    /* 11 */ { 13291,   322     },
    /* 12 */ { 34423,   451     },
    /* 13 */ { 89155,   631     },
    /* 14 */ { 230912,  884     },
    /* 15 */ { 598063,  1238    },
    /* 16 */ { 1548983, 1744    },
    };

    double estimated_entries=0;
    unsigned gp=0;

    if (n3 == 0) 
        estimated_entries = vctx->nodes.len; // we don't expect new entries coming from the next VBs

    // case: growth from n2 to n3 looks like is almost linear growth. we distiguish between 2 cases
    // 1. this is truly a uniform introduction of new entries (for example, an ID per line) - in which case the
    //    density in n1 should be about the same
    // 2. the density of n2 and n3 being so similar is a fluke - we can see a big decrease in densify from n1 - 
    //    we artificially update them to be a bit less similar
#   define UNIFORM_THREASHOLD 1.05    
    else if (n2n3_density_ratio <= UNIFORM_THREASHOLD) {
        if (n2_n3_density * 1.2 < n1_density) // case 2 
            n2n3_density_ratio = 1.06; // density in n2/n3 fell more than 20% from n1 - this is not a uniform case. arbitrarily set to 1.06
        
        else // case 1
            estimated_entries = estimated_num_lines * n2_n3_density; // growth is going to be uniform - very large hash table
    }
    
    if (n3 && n2n3_density_ratio > UNIFORM_THREASHOLD) { // if statement, not "else", because we might have updated n2n3_density_ratio

        // depending on the ratio, the number of new entries might drop fast. we cap the number of estimated vbs
        // at the VB where the expected number of new entries drops to 1 (dropping at a rate of n2n3_density_ratio every 1/3 VB)
        // n3 * [ (1/n2n3_density_ratio) ^ (3*effective_num_vbs)] = 1
        // which we develop to: effective_num_vbs = log (1/n3) / (3 * log (1/n2n3_density_ratio) )

        effective_num_vbs = ceil (log (1/n3) / (3 * log (1/n2n3_density_ratio) )); // note that n2n3_density_ratio > 1.05 and n3 >= 1 due to the if statements above

        for (gp = ARRAY_LEN(growth_plan) - 1; gp >= 0 ; gp--)
            if (MIN_(effective_num_vbs + 1 /* +1 for first vb */, estimated_num_vbs) > growth_plan[gp].vbs) {
                estimated_entries = vctx->nodes.len * n2_n3_density * growth_plan[gp].factor;
                break;
            }
    
        if (!estimated_entries) estimated_entries = 100000000000.0; // very very big - force largest available hash table
    }

    // add "stochastic noise" - to cover cases where we have a very large file, and new words continue to be introduced
    // at a low rate throughout. We add words at 10% of what we viewed in n3 - for the entire file
    if (n3_lines) estimated_entries += n3_density * estimated_num_lines * 0.10;

    #define BIG_SIZE (1 MB)
    if (estimated_entries < BIG_SIZE && (DT_FUNC(vb, seg_is_big)(vb, vctx->dict_id, (vctx->st_did_i != DID_NONE ? CTX(vctx->st_did_i)->dict_id : DICT_ID_NONE))))
        estimated_entries = BIG_SIZE; 

    if (flag.show_hash) {
 
        if (vctx->did_i==0) 
            iprintf ("\n\nOutput of --show-hash:\n"
                     "est_vbs=%u vb_1_num_lines=%s est_total_lines=%s\n", 
                     (unsigned)ceil(estimated_num_vbs), str_int_commas (vb->lines.len).s, str_int_commas ((uint64_t)estimated_num_lines).s);
        
        iprintf ("dict=%s n1=%d n2=%d n3=%d n2/n3=%2.2lf growth_plan=%u effc_vbs=%u "
                 "n2_n3_lines=%s vctx->nodes.len=%u est_entries=%d hashsize=%s\n", 
                 vctx->tag_name, (int)n1, (int)n2, (int)n3, n2n3_density_ratio, gp, (unsigned)effective_num_vbs, 
                 str_int_commas ((uint64_t)n2_n3_lines).s, vctx->nodes.len32, (int)estimated_entries, 
                 str_int_commas (hash_next_size_up (estimated_entries * 5, false)).s); 
    }

    return (uint32_t)estimated_entries;
}

void hash_alloc_global (ContextP zctx, uint32_t estimated_entries)
{
    if (!estimated_entries) estimated_entries = 1000; // this happens when populating contigs outside of seg
    
    zctx->global_hash_prime = hash_next_size_up (estimated_entries * 5, false);

    // note: we set all entries to NO_NEXT == 0xffffffff
    buf_alloc_exact_255 (evb, zctx->global_hash, zctx->global_hash_prime, uint32_t, "zctx->global_hash");
    buf_set_shared (&zctx->global_hash);

    // note: we set all entries to {NODE_INDEX_NONE, NO_NEXT, NODE_INDEX_NONE} == {0xffffffff x 3} (note: GlobalHashEnt is packed)
    buf_alloc_exact_255 (evb, zctx->global_ents, estimated_entries, GlobalHashEnt, "zctx->global_ents");
    buf_set_shared (&zctx->global_ents);
    zctx->global_ents.len = 0;

    hash_populate_from_nodes (zctx);
}


// creates a node in the hash_ents, unless the snip is already there. 
// the old node is in node, and NULL if its a new node.
// returns the node_index (positive if in nodes and negative-2 if in ston_nodes - the singleton buffer)
WordIndex hash_global_get_entry (ContextP zctx, STRp(snip), HashGlobalGetEntryMode mode,
                                 CtxNode **old_node)        // out - node if node is found, NULL if not
{    
    uint32_t *hash_p = B32(zctx->global_hash, hash_do (zctx->global_hash_prime, STRa(snip)));
    GlobalHashEnt g_head = { .next = *hash_p }, *g_hashent = &g_head;
    g_hashent->next = *hash_p; 
    int32_t hashent_i = NO_NEXT; // squash compiler warning (note: global hash table is always allocated in the first merge)
    bool singleton_encountered = false;

    while (g_hashent->next != NO_NEXT) {

        ASSERT (g_hashent->next < zctx->global_ents.len, "g_hashent->next=%d out of range in context=%s, global_ents.len=%"PRIu64,  
                g_hashent->next, zctx->tag_name, zctx->global_ents.len);

        hashent_i = g_hashent->next;

        g_hashent = B(GlobalHashEnt, zctx->global_ents, hashent_i);

        // case: snip is not in core hash table and also no other snip occupies the slot (node_index==NODE_INDEX_NONE happens only in the core table)
        if (g_hashent->node_index == NODE_INDEX_NONE) { // unoccupied space in core hash table
            // thread safety: VB threads with merge_num < ours, might be segmenting right now, and have this global hash overlayed 
            // and accessing it. we make sure that the setting of g_hashent->merge_num is atomic and the other threads will at all times either
            // see NIL or merge_num - both of which indicate the this entry is effectively NODE_INDEX_NONE as they have an older merge_num

            if (mode != HASH_READ_ONLY) {
                g_hashent->next = NO_NEXT;

                if (mode == HASH_NEW_OK_SINGLETON_IN_VB) {
                    g_hashent->node_index = -((int32_t)zctx->ston_nodes.len++) - 2; // -2 because: 0 is mapped to -2, 1 to -3 etc (as 0 is ambiguius and -1 is NODE_INDEX_NONE)
                    ASSERT (zctx->ston_nodes.len <= MAX_WORDS_IN_CTX, "zctx->ston_nodes is full: too many singletons in zctx=%s len=%u. snip=\"%.*s\"", 
                            zctx->tag_name, zctx->ston_nodes.len32, STRf(snip));
                }
                else {
                    g_hashent->node_index = zctx->nodes.len32++; 
                    ASSERT (zctx->nodes.len <= MAX_WORDS_IN_CTX, "zctx->nodes is full: too many nodes in zctx=%s len=%u. snip=\"%.*s\"", 
                            zctx->tag_name, zctx->nodes.len32, STRf(snip));
                }                

                __atomic_store_n (&g_hashent->merge_num, zctx->merge_num, __ATOMIC_RELEASE); // stamp our merge_num as the ones that set the node_index
            }

            if (old_node) *old_node = NULL; // no old node
            return g_hashent->node_index;
        }

        if (old_node) {  // note: if node=NULL, caller is telling us it is not in nodes for sure
            STR (snip_in_dict);
            *old_node = ctx_node_zf (zctx, g_hashent->node_index, &snip_in_dict, &snip_in_dict_len);
        
            // case: snip is in the hash table 
            if (str_issame (snip, snip_in_dict)) {

                // case: it is already in the hash table - as a singleton. we note that we encountered a singleton
                // and continue searching
                if (g_hashent->node_index < 0) {
                    singleton_encountered = true;
                    continue;
                }

                return g_hashent->node_index;
            }
        }
    }

    // case: not found - return if this was a read only request
    if (mode == HASH_READ_ONLY) {
        if (old_node) *old_node = NULL; // we don't have an old node
        return NODE_INDEX_NONE;
    }

    // case: not found, and we are required to provide a new entry in global_ents
    buf_alloc (evb, &zctx->global_ents, 1, 16384, GlobalHashEnt, 2, "zctx->global_ents");

    if (g_hashent != &g_head)
        g_hashent = B(GlobalHashEnt, zctx->global_ents, hashent_i); // might have changed after realloc

    // thread safetey:  VB threads with merge_num < ours, might be segmenting right now, and have this global hash overlayed 
    // and accessing it. We make sure to first prepare the new entry including the merge_num which will prohibit old
    // VBs from using it, before we atomically set the "next"
    ASSERT (zctx->global_ents.len32 < 0xffffffff, "no more room in global_ents of context %s", zctx->tag_name);
    uint32_t next = zctx->global_ents.len32++;

    GlobalHashEnt *new_hashent = B(GlobalHashEnt, zctx->global_ents, next);
    __atomic_store_n (&new_hashent->merge_num, zctx->merge_num, __ATOMIC_RELEASE); // stamp our merge_num as the ones that set the node_index
    
    // we enter the node as a singleton (=in ston_nodes) if this was a singleton in this VB and not found in any previous VB 
    // (the second occurance in the file isn't a singleton anymore)
    bool is_singleton_global   = ((mode == HASH_NEW_OK_SINGLETON_IN_VB) && !singleton_encountered);
    new_hashent->node_index    = is_singleton_global ? (-((int32_t)zctx->ston_nodes.len32++) - 2) : zctx->nodes.len32++; // -2 because: 0 is mapped to -2, 1 to -3 etc (as 0 is ambiguius and -1 is NODE_INDEX_NONE)
    new_hashent->next          = NO_NEXT;

    if (is_singleton_global) 
        ASSERT (zctx->ston_nodes.len32 <= MAX_NODE_INDEX, "number of singleton in context %s exceeded the maximum of %u. snip=\"%.*s\"", 
                zctx->tag_name, MAX_NODE_INDEX, STRf(snip));
    else
        ASSERT (zctx->nodes.len32 <= MAX_NODE_INDEX, "number of nodes in context %s exceeded the maximum of %u. snip=\"%.*s\"", 
                zctx->tag_name, MAX_NODE_INDEX, STRf(snip));

#ifdef DEBUG
    #define HASH_OCC_WARNING 2ULL
#else
    #define HASH_OCC_WARNING 10ULL
#endif
    if (zctx->nodes.len > HASH_OCC_WARNING * (uint64_t)zctx->global_hash_prime) {
        if (txt_file->redirected)
            WARN_ONCE ("Unusually slow compression due to Genozip under-allocating resources because the input file is streaming through a pipe preventing it from knowing the file size. To overcome this, please use --input-size (value in bytes, can be approximate) to inform Genozip of the file size. ctx=%s hash_prime=%u snip=\"%.*s\"", 
                       zctx->tag_name, zctx->global_hash_prime, STRf(snip));
        else
            WARN_ONCE ("Unexpected structure of file is causing unusually slow compression. ctx=%s hash_prime=%u snip=\"%.*s\"", 
                       zctx->tag_name, zctx->global_hash_prime, STRf(snip));
    }

    // now, with the new_hashent set, we can atomically update the "next"
    if (g_hashent == &g_head)
        __atomic_store_n (hash_p, next, __ATOMIC_RELEASE);
    else
        __atomic_store_n (&g_hashent->next, next, __ATOMIC_RELEASE);

    if (singleton_encountered)  // a snip that was previously counted as a singleton is encountered for the second time. it is therefore a failed singleton
        zctx->num_failed_singletons++;

    if (old_node) *old_node = NULL; // we don't have an old node
    return new_hashent->node_index;
}

// gets the node_index if the snip is already in the hash table, or puts a new one in the hash table in not
// 1. if its in the global hash table, with merge_num lower or equal to ours - i.e. added by an earler thread - we take it 
// 2. if its in the local hash table - i.e. added by us (this vb) earlier - we take it
// 3. if not found - we add it to the local hash table
WordIndex hash_get_entry_for_seg (VBlockP segging_vb, ContextP vctx, STRp(snip), 
                                  WordIndex node_index_if_new, // NODE_INDEX_NONE means "read only in global hash (ol_nodes)"
                                  CtxNode **node)        // out - node if node is found, NULL if not
{
    // first, search for the snip in the global table
    GlobalHashEnt g_head, *g_hashent = &g_head;
    
    if (vctx->global_hash.len32) {
        uint32_t hash = hash_do (vctx->global_hash_prime, STRa(snip)); // entry in hash table determined by hash function on snip
        g_hashent->next = __atomic_load_n (B32(vctx->global_hash, hash), __ATOMIC_ACQUIRE); // ACQUIRE means: compiler may not re-order load/store appearing in the code after this line, to before it
    }
    else
        g_hashent->next = NO_NEXT;

    for (unsigned depth=0; ; depth++) {
        
        // four cases for "next":
        // 1. NO_NEXT - the linked list is terminated and we have not found our snip
        // 2. next points to a valid index with our overlaid global_cache - we follow it
        // two cases where another thread sets next during merge with a higher merge number. two case:
        // 3. "next" is within the length cloned by us - we will follow it, but arrive at a node 
        //    which its merge_num prohibits us from consuming
        // 4. "next" falls outside the length cloned by us - entries added to hash, but without reallocing
        //    after we cloned. in this case, we will detected that "next" is out of range.
        uint32_t next = __atomic_load_n (&g_hashent->next, __ATOMIC_ACQUIRE); 
        
        if (next == NO_NEXT || /* case 1 */ next >= vctx->global_ents.len32) // case 4
            break;

        g_hashent = B(GlobalHashEnt, vctx->global_ents, next);
                
        // case: snip is not in core hash table (at least it wasn't there when we cloned and set our maximum merge_num we accept)
        uint32_t merge_num = __atomic_load_n (&g_hashent->merge_num, __ATOMIC_ACQUIRE);
        if (g_hashent->node_index == NODE_INDEX_NONE || merge_num > vctx->merge_num) break; // case 3

        // we skip singletons and continue searching
        if (g_hashent->node_index < 0) continue;

        STR(snip_in_dict);
        *node = ctx_node_vb (vctx, g_hashent->node_index, &snip_in_dict, &snip_in_dict_len);

        // case: snip is in the global hash table - we're done
        if ((*node)->node_index != WORD_INDEX_NONE && str_issame (snip, snip_in_dict))  // note: WORD_INDEX_NONE if node was canceled in ctx_rollback
            return g_hashent->node_index; // case 2
    }

    // caller requested read-only
    if (node_index_if_new == NODE_INDEX_NONE)
        goto not_found;

    // snip was not found in the global hash table (as it was at the time we cloned), we now search
    // in our local hash table - and if not found there - we will add it

    // allocate hash table, if not already allocated, based on experience of previous VBs, or pre-set default if there isn't any
    if (!buf_is_alloc (&vctx->local_hash)) 
        hash_alloc_local (segging_vb, vctx);

    uint32_t hash = hash_do (vctx->local_hash_prime, STRa(snip));
    LocalHashEnt l_head = { .next = *B32(vctx->local_hash, hash) }, *l_hashent = &l_head;

    int32_t l_hashent_i = NO_NEXT; // initialize to squash compiler warning

    while (l_hashent->next != NO_NEXT) {

        ASSERT (l_hashent->next < vctx->local_ents.len32, 
                "l_hashent->next=%d out of range, local_ents.len=%u", l_hashent->next, vctx->local_ents.len32);
        l_hashent_i = l_hashent->next;
        l_hashent = B(LocalHashEnt, vctx->local_ents, l_hashent_i);

        // case: snip is not in hash table and also no other snip occupies the slot (node_index==NODE_INDEX_NONE happens only in the core table)
        if (l_hashent->node_index == NODE_INDEX_NONE) { // unoccupied space in core hash table
            l_hashent->next = NO_NEXT;
            l_hashent->node_index = node_index_if_new;
            if (node) *node = NULL;
            return NODE_INDEX_NONE;
        }

        if (node) { // if the caller doesn't provide "node", he is telling us that with certainly the snip is not in the hash table
            STR(snip_in_dict);
            *node = ctx_node_vb (vctx, l_hashent->node_index, &snip_in_dict, &snip_in_dict_len);

            // case: snip is in the hash table - we're done
            if ((*node)->word_index != WORD_INDEX_NONE && str_issame (snip, snip_in_dict)) // note: WORD_INDEX_NONE if this node was canceled in ctx_rollback
                return l_hashent->node_index;
        }
    }

    // case: not found in hash table, and we are required to provide a new hash entry on the linked list
    buf_alloc (segging_vb, &vctx->local_ents, 1, 0, LocalHashEnt, 1.5, "contexts->local_ents");

    // case: first entry in local_ent for this hash value - link from local_hash
    if (l_hashent == &l_head) 
        *B32(vctx->local_hash, hash) = vctx->local_ents.len32++;

    // case: not first entry in local_ent for this hash value - link from previous entry with same hash
    else
        B(LocalHashEnt, vctx->local_ents, l_hashent_i)->next = vctx->local_ents.len32++; // note: don't use l_hashent as it might have changed after realloc

    LocalHashEnt *new_l_hashent = BLST(LocalHashEnt, vctx->local_ents);
    new_l_hashent->next = NO_NEXT;
    new_l_hashent->node_index = node_index_if_new;

not_found:
    if (node) *node = NULL;
    return NODE_INDEX_NONE;
}
