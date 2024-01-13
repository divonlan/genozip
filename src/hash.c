// ------------------------------------------------------------------
//   hash.c
//   Copyright (C) 2020-2024 Genozip Limited
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

#pragma pack(4)
typedef struct {              // 8 bytes
    WordIndex node_index;     // index into vctx->nodes or WORD_INDEX_NONE
    uint32_t next;            // linked list - index into local_ents or NO_NEXT
                              //               local_ents indices start at LOCAL_HASH_OFFSET
} LocalHashEnt, StonHashEnt;

#pragma pack()

#define ENT_NONE 0xffffffffUL // no global_hash/stons_hash entry
#define ENT_HEAD 0xfffffffeUL // the head rather than ents buffer

#ifdef DEBUG
    #define HASH_OCC_WARNING 2ULL
#else
    #define HASH_OCC_WARNING 10ULL
#endif

size_t sizeof_GlobalHashEnt (void) { return sizeof (GlobalHashEnt); }
size_t sizeof_StonHashEnt   (void) { return sizeof (StonHashEnt);   }

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

    for (CtxWord *node = B1ST(CtxWord, zctx->nodes); node < B(CtxWord, zctx->nodes, len); node++) {
        rom snip = Bc(zctx->dict, node->char_index);

        CtxWord *please_update_index;
        hash_global_get_entry (zctx, snip, node->snip_len, false, true, &please_update_index); 
        please_update_index->char_index = node->char_index;
    }
}

// This is called when the VB encounters a first snip that's not in the ol_dict 
// allocation algorithm:
// 1. If we got info on the size of this dict with the previous merged vb - use that size
// 2. If not - use either num_lines for the size, or the smallest size for dicts that are typically small
static void hash_alloc_local (VBlockP vb, ContextP vctx)
{
    vctx->local_hash_prime = 0; // initialize

    // if known from previously merged vb - use those values
    if (vctx->num_new_entries_prev_merged_vb)
        // 3X the expected number of entries to reduce hash contention
        vctx->local_hash_prime = hash_next_size_up (vctx->num_new_entries_prev_merged_vb * 3, false);

    // if known to small, use hash table of ~ 64K
    else if (DT_(vb, seg_is_small)(vb, vctx->dict_id))
        vctx->local_hash_prime = hash_next_size_up(1, false);
    
    // default: it could be big - start with num_lines / 10 (this is an estimated num_lines that is likely inflated)
    else
        vctx->local_hash_prime = hash_next_size_up (vb->lines.len32 / 10, false);

    // note: we can't be too generous with the initial allocation because this memory is usually physically allocated
    // to ALL VB structures before any of them merges. Better start smaller for vb_i=1 and let it extend if needed
    buf_alloc_exact_255 (vb, vctx->local_hash, vctx->local_hash_prime, uint32_t, "contexts->local_hash");

    buf_alloc_exact_255 (vb, vctx->local_ents, vctx->local_hash_prime / 3, LocalHashEnt, "contexts->local_ents");
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

    // note: we set all entries to {WORD_INDEX_NONE, NO_NEXT, WORD_INDEX_NONE} == {0xffffffff x 3} (note: GlobalHashEnt is packed)
    buf_alloc_exact_255 (evb, zctx->global_ents, estimated_entries, GlobalHashEnt, "zctx->global_ents");
    buf_set_shared (&zctx->global_ents);
    zctx->global_ents.len = 0;

    hash_populate_from_nodes (zctx);
}

// search for snip in singletons - returns true and result if found
// ston_hash  - fixed size hash table, type uint32_t points to head of linked-list - index into ston_ents
// ston_ents  - linked lists associated with each hash value, type StonHashEnt
// ston_nodes - Nodes of singletons - of type CtxWord (contains snip index and length)
// stons      - buffer containing the snips. each snip is nul-terminated.
// Note: these ^ buffers are protected by the zctx mutex and are not overlayed to the compute threads, so no thread safety issues.

// returns: false - singleton not found ; true - singleton was found and removed
static inline bool hash_stons_remove_singleton (ContextP zctx, uint32_t hash, STRp(snip))
{
    if (!zctx->ston_ents.len32) return false; // this zctx has no singletons

    uint32_t *head = B32(zctx->ston_hash, hash);
    uint32_t *prevs_next = head;

    StonHashEnt dummy = { .next = *head }, *g_stonent = &dummy;

    while (g_stonent->next != NO_NEXT) {
        ASSERT (g_stonent->next < zctx->ston_ents.len32, "g_stonent->next=%d out of range in context=%s, ston_ents.len=%"PRIu64,  
                g_stonent->next, zctx->tag_name, zctx->ston_ents.len);

        g_stonent = B(StonHashEnt, zctx->ston_ents, g_stonent->next);

        // get node and snip from ston_nodes/stons and check if it is our snip. 
        CtxWord ston_node = *B(CtxWord, zctx->ston_nodes, g_stonent->node_index);
        rom snip_in_stons = Bc(zctx->stons, ston_node.char_index);
            
        // case: this snip encountered in the current VB is already a singleton entered by a prior VB,
        // therefore it is no longer a singleton - and we remove it
        if (str_issame_(STRa(snip), snip_in_stons, ston_node.snip_len)) {
            // remove singleton from its hash linked list
            *prevs_next = g_stonent->next;
            
            // add the singleton to the linked list of all decommissioned singletons waiting to recycle -
            // headed by the last entry in ston_hash
            head = BLST32 (zctx->ston_hash); 

            // LIFO - insert newly decommissioned singleton at the head (LIFO is just easier to implement that FIFO)
            g_stonent->next = *head;
            *head = BNUM (zctx->ston_ents, g_stonent);

            // We call this a failed singleton, because the initial VB thoughts its a singleton and segged
            // it in local, and now we're moving it to dict - so it appears in z twice.
            zctx->num_failed_singletons++; 
            
            return true; // singleton found and removed
        }

        prevs_next = &g_stonent->next;
    }

    return false; // we've reached the end of the linked list and the singleton was not found
}

static void hash_verify_nodes_space (ContextP zctx, BufferP nodes, rom item_name, STRp(snip))
{
    ASSERT (nodes->len < 0xffffffff, "number of %s in context %s exceeded the maximum of %u. snip=%s", 
            item_name, zctx->tag_name, 0xfffffffe, str_snip (STRa(snip), true).s);

    if (nodes->len > HASH_OCC_WARNING * (uint64_t)zctx->global_hash_prime) {
        if (txt_file->redirected)
            WARN_ONCE ("(%s) Unusually slow compression due to Genozip under-allocating resources because the input file is streaming through a pipe preventing it from knowing the file size. To overcome this, please use --input-size (value in bytes, can be approximate) to inform Genozip of the file size. ctx=%s hash_prime=%u snip=\"%s\"", 
                       item_name, zctx->tag_name, zctx->global_hash_prime, str_snip (STRa(snip), true).s);
        else
            WARN_ONCE ("(%s) Unexpected structure of file is causing unusually slow compression. ctx=%s hash_prime=%u snip=\"%s\"%s", 
                       item_name, zctx->tag_name, zctx->global_hash_prime, str_snip (STRa(snip), true).s, SUPPORT);
    }
}

// add a singleton to the singleton buffers
static inline void hash_stons_add_singleton (ContextP zctx, uint32_t hash, STRp(snip), CtxWord **please_update_index)
{
    // case: first singleton in this zctx - allocate - set all entries to NO_NEXT == 0xffffffff
    if (!zctx->ston_hash.len32) 
        // global_hash_prime linked-lists - one for each hash value, +1 for decommissions singletons
        buf_alloc_exact_255 (evb, zctx->ston_hash, zctx->global_hash_prime + 1, uint32_t, "zctx->ston_hash"); 

    // first, try to find a decommissioned singleton of the same length
    uint32_t *decommissioned_head = BLST32(zctx->ston_hash);
    uint32_t *decommissioned_prev_next = decommissioned_head;
    StonHashEnt dummy = { .next = *decommissioned_head }, *g_stonent = &dummy;

    while (g_stonent->next != NO_NEXT) {
        ASSERT (g_stonent->next < zctx->ston_ents.len32, "g_stonent->next=%d out of range in context=%s, ston_ents.len=%"PRIu64,  
                g_stonent->next, zctx->tag_name, zctx->ston_ents.len);

        g_stonent = B(StonHashEnt, zctx->ston_ents, g_stonent->next);
        CtxWord *ston_node = B(CtxWord, zctx->ston_nodes, g_stonent->node_index);

        // case: found a decommissioned singleton of the same length - update the snip string and move between linked-lists
        if (ston_node->snip_len == snip_len) {
            // update snip string
            *please_update_index = NULL; // ctx_commit_node doesn't need to update since we're updating here
            memcpy (Bc(zctx->stons, ston_node->char_index), snip, snip_len);

            // remove entry from decommissioned singletons linked list: connect previous entry to next one
            *decommissioned_prev_next = g_stonent->next; 

            // connect recommissioned entry to the start of the linked-list of hash (just because its easier than to add to the end)
            // before: head -> old_ent1 -> old_ent2...
            // after : head -> recom_ent -> old_ent1 -> old_ent2...
            g_stonent->next = *B32(zctx->ston_hash, hash);
            *B32(zctx->ston_hash, hash) = BNUM (zctx->ston_ents, g_stonent);

            return; // we re-commissioned the decommissioned singleton with the new snip - nothing more to do
        }

        decommissioned_prev_next = &g_stonent->next;
    }

    // case: no decommissioned singleton of same length. creating new singleton entry
    buf_alloc (evb, &zctx->ston_ents, 1, 16 KB, StonHashEnt, 2, "zctx->ston_ents");
    ASSERT (zctx->ston_ents.len32 < 0xffffffff, "no more room in ston_ents of context %s", zctx->tag_name);

    // connect new entry at the start of the linked-list of hash
    // before: head -> old_ent1 -> old_ent2...
    // after : head -> new_ent -> old_ent1 -> old_ent2...
    BNXT(StonHashEnt, zctx->ston_ents) = (StonHashEnt){ .node_index = zctx->ston_nodes.len32,
                                                        .next = *B32(zctx->ston_hash, hash) };
         
    *B32(zctx->ston_hash, hash) = zctx->ston_ents.len32 - 1; 

    hash_verify_nodes_space (zctx, &zctx->ston_nodes, "singletons", STRa(snip));
    buf_alloc (evb, &zctx->ston_nodes, 1, INITIAL_NUM_NODES, CtxWord, CTX_GROWTH, "zctx->ston_nodes");
    
    CtxWord *ston_node = &BNXT(CtxWord, zctx->ston_nodes);
    ston_node->snip_len = snip_len;
    *please_update_index = ston_node; // ctx_commit_node needs to update char_index
}

// search for snip in (non-singleton) nodes
static inline bool hash_global_find_in_nodes (ContextP zctx, uint32_t hash, STRp(snip), bool compare_to_snip, WordIndex *existing_node_index, uint32_t *last_ent_on_linked_list)
{
    uint32_t *head = B32(zctx->global_hash, hash);
    GlobalHashEnt dummy = { .next = *head }, *g_hashent = &dummy;
    g_hashent->next = *head; 

    while (g_hashent->next != NO_NEXT) {
        ASSERT (g_hashent->next < zctx->global_ents.len32, "g_hashent->next=%d out of range in context=%s, global_ents.len=%"PRIu64,  
                g_hashent->next, zctx->tag_name, zctx->global_ents.len);

        g_hashent = B(GlobalHashEnt, zctx->global_ents, g_hashent->next);

        if (compare_to_snip) { // we're actually trying to find the snip, not just traverse to the end of the linked list
            CtxWord node = *B(CtxWord, zctx->nodes, g_hashent->word_index);
            rom snip_in_dict = Bc(zctx->dict, node.char_index);
        
            // case: found snip's node 
            if (str_issame_(STRa(snip), snip_in_dict, node.snip_len)) {
                if (existing_node_index) *existing_node_index = g_hashent->word_index;
                return true;
            }
        }
    }

    *last_ent_on_linked_list = (g_hashent == &dummy) ? ENT_HEAD : BNUM (zctx->global_ents, g_hashent);

    return false; // we've reached the end of the linked list and the node was not found
}

// add a new (non-singleton) node
static inline WordIndex hash_global_add_node (ContextP zctx, uint32_t hash, uint32_t last_ent_on_linked_list, STRp(snip), 
                                              CtxWord **please_update_index)
{
    // realloc
    hash_verify_nodes_space (zctx, &zctx->nodes, "nodes", STRa(snip));

    buf_alloc (evb, &zctx->global_ents, 1, 16 KB, GlobalHashEnt, 2, "zctx->global_ents");
    buf_alloc (evb, &zctx->nodes, 1, INITIAL_NUM_NODES, CtxWord, CTX_GROWTH, "zctx->nodes");

    CtxWord *new_node = &BNXT(CtxWord, zctx->nodes);
    *new_node = (CtxWord){ .snip_len = snip_len }; 
    *please_update_index = new_node;

    // thread safetey:  VB threads with merge_num < ours, might be segmenting right now, and have this global hash overlayed 
    // and accessing it. We make sure to first prepare the new entry including the merge_num which will prohibit old
    // VBs from using it, before we atomically set the "next"
    ASSERT (zctx->global_ents.len32 < 0xffffffff, "no more room in global_ents of context %s", zctx->tag_name);
    uint32_t new_ent_index = zctx->global_ents.len32++;

    // note: when traversing the linked list, the merge_num is monotonically increasing, allowing the
    // compute thread to stop at the first merge_num that is beyond their's.
    GlobalHashEnt *new_hashent = B(GlobalHashEnt, zctx->global_ents, new_ent_index);
    __atomic_store_n (&new_hashent->merge_num, zctx->merge_num, __ATOMIC_RELEASE); // stamp our merge_num as the ones that set the word_index
    
    new_hashent->word_index = zctx->nodes.len32 - 1; 
    new_hashent->next       = NO_NEXT;

    // now, with the new_hashent set, we can atomically update link the previous entry to the new one.
    // note: we must add end the end (rather than start) of the linked-list, to maintain the merge_nums in order.
    uint32_t *prevs_next = (last_ent_on_linked_list == ENT_HEAD) 
        ? B32(zctx->global_hash, hash) : &B(GlobalHashEnt, zctx->global_ents, last_ent_on_linked_list)->next;

    __atomic_store_n (prevs_next, new_ent_index, __ATOMIC_RELEASE);

    return zctx->nodes.len32 - 1; // word_index
}

// finds or creates a node in the hash_ents or ston_ents
// returns the word_index if node and WORD_INDEX_NONE if singleton.
// returns char_index_p - if not NULL - caller should allocate space in dict/stons and copy snip 
// (if NULL this is a re-commissioned singleton and snip already copied)
WordIndex hash_global_get_entry (ContextP zctx, STRp(snip), 
                                 bool allow_singleton,
                                 bool snip_is_definitely_new,
                                 CtxWord **please_update_index)  // out: in not NULL, caller should update index 
{    
    uint32_t hash = hash_do (zctx->global_hash_prime, STRa(snip));
    uint32_t last_ent_on_linked_list=ENT_NONE;

    bool is_previously_a_ston = !snip_is_definitely_new && hash_stons_remove_singleton (zctx, hash, STRa(snip));
    
    WordIndex existing_node_index = WORD_INDEX_NONE;
    bool is_in_nodes = !snip_is_definitely_new && 
                       !is_previously_a_ston && 
                       hash_global_find_in_nodes (zctx, hash, STRa(snip), true, &existing_node_index, &last_ent_on_linked_list);
    
    // case: new singleton 
    if (!is_in_nodes && !is_previously_a_ston && allow_singleton) {
        hash_stons_add_singleton (zctx, hash, STRa(snip), please_update_index);
        return WORD_INDEX_NONE; // means "singleton"
    }

    // case: new node (not a singleton in this VB, or seen exactly once in a previous VB and stored as a singleton)
    else if (!is_in_nodes) {
        // get the end of the linked list, if we don't already have it
        if (last_ent_on_linked_list == ENT_NONE)
            hash_global_find_in_nodes (zctx, hash, 0, 0, false, NULL, &last_ent_on_linked_list);

        return hash_global_add_node (zctx, hash, last_ent_on_linked_list, STRa(snip), please_update_index);
    }

    // case: existing node
    else {
        *please_update_index = NULL; // no need to update - this node already exists
        return existing_node_index;
    }
}

// gets the node_index if the snip is already in the hash table, or puts a new one in the hash table in not
// 1. if its in the global hash table, with merge_num lower or equal to ours - i.e. added by an earler thread - we take it 
// 2. if its in the local hash table - i.e. added by us (this vb) earlier - we take it
// 3. if not found - we add it to the local hash table, and return WORD_INDEX_NONE
// note: during seg, compute threads don't have access to global singletons. 
WordIndex hash_get_entry_for_seg (VBlockP vb, ContextP vctx, STRp(snip), 
                                  WordIndex node_index_if_new, // WORD_INDEX_NONE means "read only in global hash (ol_nodes)"
                                  rom *snip_in_dict_out)       // optional out - snip - pointer into vctx->dict or vctx->ol_dict (only if existing, NULL if not)
{
    // first, search for the snip in the global table
    GlobalHashEnt dummy, *g_hashent = &dummy;
    
    if (vctx->global_hash.len32) {
        uint32_t hash = hash_do (vctx->global_hash_prime, STRa(snip)); // entry in hash table determined by hash function on snip
        g_hashent->next = __atomic_load_n (B32(vctx->global_hash, hash), __ATOMIC_ACQUIRE); // ACQUIRE means: compiler may not re-order load/store appearing in the code after this line, to before it
    }
    else
        g_hashent->next = NO_NEXT;

    for (unsigned __attribute__((unused)) depth=0; ; depth++) {
        
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
        if (g_hashent->word_index == WORD_INDEX_NONE || merge_num > vctx->merge_num) break; // case 3

        STR(snip_in_dict);
        ctx_node_vb (vctx, g_hashent->word_index, &snip_in_dict, &snip_in_dict_len);

        // case: snip is in the global hash table - we're done
        if (str_issame (snip, snip_in_dict)) { 
            if (snip_in_dict_out) *snip_in_dict_out = snip_in_dict; // pointer into vctx->ol_dict
            return g_hashent->word_index; // case 2
        }
    }

    // caller requested read-only and only from ol_nodes
    if (node_index_if_new == WORD_INDEX_NONE)
        goto not_found;

    // snip was not found in the global hash table (as it was at the time we cloned), we now search
    // in our local hash table - and if not found there - we will add it

    // allocate hash table, if not already allocated, based on experience of previous VBs, or pre-set default if there isn't any
    if (!buf_is_alloc (&vctx->local_hash)) 
        hash_alloc_local (vb, vctx);

    uint32_t hash = hash_do (vctx->local_hash_prime, STRa(snip));
    LocalHashEnt l_head = { .next = *B32(vctx->local_hash, hash) }, *l_hashent = &l_head;

    int32_t l_hashent_i = NO_NEXT; 

    while (l_hashent->next != NO_NEXT) {

        ASSERT (l_hashent->next < vctx->local_ents.len32, 
                "l_hashent->next=%d out of range, local_ents.len=%u", l_hashent->next, vctx->local_ents.len32);
        l_hashent_i = l_hashent->next;
        l_hashent = B(LocalHashEnt, vctx->local_ents, l_hashent_i);

        // case: snip is not in hash table and also no other snip occupies the slot (node_index==WORD_INDEX_NONE happens only in the core table)
        if (l_hashent->node_index == WORD_INDEX_NONE) { // unoccupied space in core hash table
            l_hashent->next = NO_NEXT;
            l_hashent->node_index = node_index_if_new;
            if (snip_in_dict_out) *snip_in_dict_out = NULL;
            return WORD_INDEX_NONE; // new node
        }

        STR(snip_in_dict);
        CtxVbNode node = ctx_node_vb (vctx, l_hashent->node_index, &snip_in_dict, &snip_in_dict_len);

        // case: snip is in the local hash table - we're done
        if (node.canceled != VB_NODE_CANCELED && str_issame (snip, snip_in_dict)) {
            if (snip_in_dict_out) *snip_in_dict_out = snip_in_dict; // pointer into vctx->dict
            return l_hashent->node_index;
        }
    }

    // case: not found in either hash table, and we are required to provide a new hash entry on the linked list
    buf_alloc (vb, &vctx->local_ents, 1, 0, LocalHashEnt, 1.5, "contexts->local_ents");

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
    if (snip_in_dict_out) *snip_in_dict_out = NULL;
    return WORD_INDEX_NONE;
}
