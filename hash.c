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
#include "strings.h"

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
    // primary numbers just beneath the powers of 2^0.5 (and 2^0.25 for the larger numbers)
    // minimum ~64K to prevent horrible miscalculations in edge cases that result in dramatic slow down
    static uint32_t hash_sizes[] = { /* 509, 719, 1021, 1447, 2039, 2887, 4039, 5791, 8191, 11579, 16381, 23167, 32749, 46337,*/ 65521, 92681, 131071, 
                                     185363, 262139, 370723, 524287, 741431, 1048573, 1482907, 2097143, 2965819, 4194301, 5931641, 8388593, 
                                     11863279, 16777213, 19951579, 23726561, 28215799, 33554393, 39903161, 47453111, 56431601, 67108859 };
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

    else switch (segging_vb->data_type) {
    
    case DT_VCF:
        // if typically small - use minimal hash table
        if (vb_ctx->dict_id.num == dict_id_fields[VCF_CHROM]  ||
            vb_ctx->dict_id.num == dict_id_fields[VCF_FORMAT] ||
            vb_ctx->dict_id.num == dict_id_fields[VCF_INFO]   ||
            vb_ctx->dict_id.num == dict_id_fields[VCF_REFALT] ||
            vb_ctx->dict_id.num == dict_id_fields[VCF_FILTER] ||
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

            vb_ctx->local_hash_prime = hash_next_size_up((uint32_t)segging_vb->lines.len);
        break;

    case DT_SAM:
        // typically small - use minimal hash table ~ 500
        if (vb_ctx->dict_id.num == dict_id_fields[SAM_FLAG]      ||
            vb_ctx->dict_id.num == dict_id_fields[SAM_MAPQ]      ||
            vb_ctx->dict_id.num == dict_id_fields[SAM_QNAME]     ||
            vb_ctx->dict_id.num == dict_id_fields[SAM_OPTIONAL]  ||

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
            
            // biobambam tags
            
            vb_ctx->dict_id.num == dict_id_OPTION_STRAND)
            
            vb_ctx->local_hash_prime = hash_next_size_up(500);

        // typically smallish - use hash table ~ 2000
        else 
        if (vb_ctx->dict_id.num == dict_id_fields[SAM_RNAME]  ||
            vb_ctx->dict_id.num == dict_id_OPTION_CC)

            vb_ctx->local_hash_prime = hash_next_size_up(2000);

        // typically medium - use hash table ~ 50000
        else 
        if (vb_ctx->dict_id.num == dict_id_fields[SAM_CIGAR]  ||
            vb_ctx->dict_id.num == dict_id_OPTION_MC)

            vb_ctx->local_hash_prime = hash_next_size_up(50000);
        break;

    case DT_FASTQ:
    case DT_FASTA:
        if (vb_ctx->dict_id.num == dict_id_fields[FAST_DESC] ||
            vb_ctx->dict_id.num == dict_id_fields[FAST_LINEMETA])
            
            vb_ctx->local_hash_prime = hash_next_size_up(500);
        break;

    case DT_GFF3:
        if (vb_ctx->dict_id.num == dict_id_fields[GFF3_SEQID] ||
            vb_ctx->dict_id.num == dict_id_fields[GFF3_SOURCE] ||
            vb_ctx->dict_id.num == dict_id_fields[GFF3_TYPE] ||
            vb_ctx->dict_id.num == dict_id_fields[GFF3_END] ||
            vb_ctx->dict_id.num == dict_id_fields[GFF3_SCORE] ||
            vb_ctx->dict_id.num == dict_id_fields[GFF3_STRAND] ||
            vb_ctx->dict_id.num == dict_id_fields[GFF3_PHASE] ||
            vb_ctx->dict_id.num == dict_id_fields[GFF3_ATTRS])
            
            vb_ctx->local_hash_prime = hash_next_size_up(500);
        break;

    case DT_ME23:
        if (vb_ctx->dict_id.num == dict_id_fields[ME23_CHROM])
            
            vb_ctx->local_hash_prime = hash_next_size_up(500);
        break;


    default: ABORT ("hash_alloc_local: unknown data_type=%s", dt_name (segging_vb->data_type));
    }

    // default: it could be big - start with num_lines / 10 (this is an estimated num_lines that is likely inflated)
    if (!vb_ctx->local_hash_prime) 
        vb_ctx->local_hash_prime = hash_next_size_up ((uint32_t)segging_vb->lines.len / 10);

    // note: we can't be too generous with the initial allocation because this memory is usually physically allocated
    // to ALL VB structures before any of them merges. Better start smaller for vb_i=1 and let it extend if needed
    buf_alloc (segging_vb, &vb_ctx->local_hash, (vb_ctx->local_hash_prime * 1.2) * sizeof (LocalHashEnt) /* room for expansion */, 1, 
               "mtf_ctx->local_hash", vb_ctx->did_i);
    vb_ctx->local_hash.len = vb_ctx->local_hash_prime;
    memset (vb_ctx->local_hash.data, 0xff, vb_ctx->local_hash_prime * sizeof (LocalHashEnt)); // initialize core table
//printf ("Seg vb_i=%u: local hash: dict=%.8s size=%u\n", segging_vb->vblock_i, err_dict_id (vb_ctx->dict_id), vb_ctx->local_hash_prime); 
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
    // note on txt_data_size_single: if its a physical plain txt file - this is the file size. 
    // if not - its an estimate done after the first VB by txtfile_estimate_txt_data_size
    double effective_num_vbs=0, estimated_num_vbs = MAX (1, (double)txt_file->txt_data_size_single / (double)merging_vb->txt_data.len);
    double estimated_num_lines = estimated_num_vbs * (double)merging_vb->lines.len;

    // for growth purposes, we discard the first 1/3 of VB1, and compare the growth of the 2nd vs the 3rd 1/3. 
    // this is because many fields display the charactistics of having a relatively small number of high frequency entries
    // which mostly show up in n1, and then a long tail of low frequency entries. The comparison of n2 to n3,
    // assuming that the new snips first introduced in them are mostly low frequency ones, will give us a more accurate predication
    // of the gradient appearance of new snips
    double n1 = first_merging_vb_ctx->mtf_len_at_1_3;
    double n2 = first_merging_vb_ctx->mtf_len_at_2_3 ? ((double)first_merging_vb_ctx->mtf_len_at_2_3 - (double)first_merging_vb_ctx->mtf_len_at_1_3) : 0;
    double n3 = (double)first_merging_vb_ctx->ol_mtf.len - n1 - n2;

    double n1_lines      = (double)merging_vb->num_lines_at_1_3;
    double n2_lines      = (double)merging_vb->num_lines_at_2_3 - n1_lines;
    double n2_n3_lines   = (double)merging_vb->lines.len - n1_lines;
    double n3_lines      = (double)merging_vb->lines.len - n1_lines - n2_lines;

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
        estimated_entries = zf_ctx->mtf.len; // we don't expect new entries coming from the next VBs

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

        for (gp = sizeof(growth_plan) / sizeof(growth_plan[0]) - 1; gp >= 0 ; gp--)
            if (MIN (effective_num_vbs + 1 /* +1 for first vb */, estimated_num_vbs) > growth_plan[gp].vbs) {
                estimated_entries = zf_ctx->mtf.len * n2_n3_density * growth_plan[gp].factor;
                break;
            }
    
        if (!estimated_entries) estimated_entries = 100000000000; // very very big - force largest available hash table
    }

    // add "stochastic noise" - to cover cases where we have a very large file, and new words continue to be introduced
    // at a low rate throughout. We add words at 10% of what we viewed in n3 - for the entire file
    if (n3_lines) estimated_entries += n3_density * estimated_num_lines * 0.10;

    zf_ctx->global_hash_prime = hash_next_size_up (estimated_entries * 5);

    if (flag_show_hash) {
        char s1[30], s2[30];
 
        if (zf_ctx->did_i==0) {
            fprintf (stderr, "\n\nOutput of --show-hash:\n");
            fprintf (stderr, "est_vbs=%u vb_1_num_lines=%s est_total_lines=%s\n", 
                     (unsigned)ceil(estimated_num_vbs), str_uint_commas (merging_vb->lines.len, s1), str_uint_commas ((uint64_t)estimated_num_lines, s2));
        }
        
        fprintf (stderr, "dict=%.8s n1=%d n2=%d n3=%d n2/n3=%2.2lf growth_plan=%u effc_vbs=%u "
                 "n2_n3_lines=%s zf_ctx->mtf.len=%u est_entries=%d hashsize=%s\n", 
                 dict_id_printable(zf_ctx->dict_id).id, (int)n1, (int)n2, (int)n3, n2n3_density_ratio, gp, (unsigned)effective_num_vbs, 
                 str_uint_commas ((uint64_t)n2_n3_lines, s2), (uint32_t)zf_ctx->mtf.len, (int)estimated_entries, 
                 str_uint_commas (zf_ctx->global_hash_prime, s1)); 
    }

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
