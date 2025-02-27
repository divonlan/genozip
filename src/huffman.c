// ------------------------------------------------------------------
//   huffman.c
//   Copyright (C) 2023-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "bits.h"
#include "huffman.h"
#include "buffer.h"
#include "flags.h"
#include "vblock.h"
#include "compressor.h"
#include "file.h"
#include "zfile.h"
#include "sam_friend.h"

#define ALPHABET 256

typedef struct { // 12 bytes
    uint32_t freq;
    uint8_t c;
    uint8_t unused;
    uint16_t parent, left, right;
} HuffmanNode;

#define HUFFMAN_MAX_MASTER_LEN 62 // note: with 62, HuffmanCodes fits in 4B words (feels safter...)
#define HUFFMAN_MAX_MAX_PREFIX_LEN 31 // MAX_MAX_ is intentional 
#define HUFFMAN_PREFIX_LEN_BITS 5 // max_prefix_len fits in this many bits.

// ZIP/PIZ: used for compressing data (this is the payload of the SEC_HUFFMAN section)
typedef struct HuffmanCodes { // 9548 bytes
    uint64_t codes[ALPHABET];      // huffman variable-length representation of each value (up to 64 bits)
    uint8_t code_n_bits[ALPHABET]; // number of bits in huffman representation of each value
    HuffmanNode nodes[2*ALPHABET]; // 256 leafs, 255 internal nodes, 1 nil node
    int roots[ALPHABET];           // indices into nodes
    Did did_i; 
    uint8_t longest_code_n_bits;
    uint8_t ununsed;
    uint32_t num_chewed;           // ZIP only
    uint32_t unused1;            
    uint8_t max_prefix_len;        // subset of master for which we reasonably expect identical prefixes
    uint8_t master_len;
    char master[HUFFMAN_MAX_MASTER_LEN]; // used for finding an identical prefix, and XORing the rest
} HuffmanCodes;

#define MAX_CHEWED 10000 // maximum number of samples after which we stop chewing

// ZIP: used for calculating codes (=huffman representation of values) 
typedef struct HuffmanChewer { // 1348 bytes
    uint32_t num_chewed;
    #define CHEW_UNIT 256 /*larger than adding the frequencies of all non-encountered values (which have freq=1)*/
    uint32_t freq_table[ALPHABET]; // number of occurances of each chewed value
    uint8_t max_prefix_len;
    uint8_t master_len;
    char master[HUFFMAN_MAX_MASTER_LEN];
    #define HUFFMAN_MAX_PREV_LEN 255
    uint8_t prev_sample_len;
    char prev_sample[HUFFMAN_MAX_PREV_LEN];
} HuffmanChewer;

typedef packed_enum { 
    NC_N_CIGAR_OP, NC_N_MEX, NC_N_SH, NC_N_ID, NC_N_NP, // channels containing a uint8_t
    NC_OP_FLANK, NC_OP_AFTER_SHNP, NC_OP_AFTER_ME, NC_OP_AFTER_X, NC_OP_AFTER_D, NC_OP_AFTER_I, // channels containing a BamCigarOpType
    NUM_NICO_CHANNELS } NicoChannels;

#define NICO_CHAN_NAMES { "N_CIGAR_OP", "N_MEX", "N_SH", "N_ID", "N_NP", "OP_FLANK", "OP_AFTER_SHNP", "OP_AFTER_ME", "OP_AFTER_X", "OP_AFTER_D", "OP_AFTER_I" }

#define ASSHUFFEXISTS ASSERT (huffman_exists(did_i), "%s.huffman doesn't exist", zctx->tag_name)
// -------------------------------------------------------------------------------------
// private - for use by huffman_compress and nico functions
// -------------------------------------------------------------------------------------

static inline int huffman_compress_init_bits (STR8c(comp), BitsP bits, uint64_t *next_bit)
{
    // shift back so Bits starts from full word
    int shift = (uint64_t)comp & 0b111; // bytes shift 

    *next_bit = shift * 8;

    *bits = (Bits){ .nbits  = (comp_len + shift) * 8,    // note: possibly only partially covering last word
                    .nwords = ROUNDUP8(comp_len + shift) / 8,
                    .words  = (uint64_t *)(comp - shift), // shift to bring address to word boundary
                    .type   = BITS_STANDALONE };

    return shift;
}

static inline void huffman_compress_finalize_bits (ContextP zctx, BitsP bits, uint64_t next_bit, int shift, qSTR8p(comp))
{
    ASSERT (next_bit <= bits->nbits, "bits overflow: expecting next_bit=%"PRIu64" <= nbits=%"PRIu64" for %s",
            next_bit, bits->nbits, zctx->tag_name);

    bits_truncate (bits, next_bit);

    uint32_t updated_comp_len = roundup_bits2bytes(next_bit) - shift; // number of bytes
    ASSERT (updated_comp_len <= *comp_len, "compression overflow: available %u bytes but generated %u bytes", *comp_len, updated_comp_len);

    *comp_len = updated_comp_len;

    // clear unused bits in last byte
    int bits_in_last_byte = next_bit % 8;
    if (bits_in_last_byte) 
        comp[*comp_len-1] &= bitmask8(bits_in_last_byte);
}

#define huffman_compress_add_int(word,n_bits,avail_comp_len) \
    ({ ASSERT (next_bit + (n_bits) <= bits.nbits, "%s: comp_len=%d too short when compressing %s", LN_NAME, (avail_comp_len), zctx->tag_name); \
       bits_set_wordn (&bits, next_bit, (word), (n_bits)); \
       next_bit += (n_bits); })

#define huffman_compress_add_one(x_,avail_comp_len) \
    ({ uint8_t x = (x_); /* evaluate x_ just once */\
       ASSERT (h->code_n_bits[x], "%s: value %u in not in the huffman alphabet of %s channel_i=%u", LN_NAME, x, zctx->tag_name,(int)(h - B1ST(HuffmanCodes, zctx->huffman))); \
       huffman_compress_add_int (h->codes[x], h->code_n_bits[x], (avail_comp_len)); })

#define huffman_uncompress_init_bits                                        \
    Bits bits = { .nbits  = 64000000000ULL, /* some very large number */    \
                  .nwords = 1000000000ULL,                                  \
                  .type   = BITS_STANDALONE,                                \
                  .words  = (uint64_t *)ROUNDDOWN8 ((uint64_t)comp) }; /* shift back to align with a 64b word boundary (counting on pointers being up to 64 bit) */ \
    uint64_t bit_i = (comp - (bytes)bits.words) * 8; /* note: if we shifted .words back, the initial value of bit_i will point to were the data actually starts */

#define huffman_uncompress_comp_len (roundup_bits2bytes(bit_i) - (comp - (bytes)bits.words))

#define huffman_uncompress_get_c ({                                         \
    int node_i = h->roots[0];                                               \
    while (h->nodes[node_i].left)  /* has children */                       \
        node_i = bits_get (&bits, bit_i++) ? h->nodes[node_i].right : h->nodes[node_i].left; \
    h->nodes[node_i].c; })

//-----------------------------------------
// plain vanilla Huffman
//-----------------------------------------

// main thread
void huffman_start_chewing (Did did_i, STRp(master), uint8_t max_prefix_len, int n_channels)
{
    decl_zctx(did_i);

    ASSERTINRANGX (max_prefix_len, 0, HUFFMAN_MAX_MAX_PREFIX_LEN);
    if (master) ASSERTNOTZERO (max_prefix_len);

    buf_alloc_exact_zero (evb, zctx->huffman, n_channels, HuffmanCodes, "contexts->huffman"); // used for both HuffmanChewer and HuffmanCodes - the latter is much bigger
    
    for_buf2 (HuffmanChewer, chewer, i, zctx->huffman) {
        if (i==0 && master_len) {
            chewer->master_len = MIN_(master_len, HUFFMAN_MAX_MASTER_LEN); 
            memcpy (chewer->master, master, chewer->master_len);

            chewer->max_prefix_len = MIN_(max_prefix_len, master_len);
        }

        for (int i=0; i < ALPHABET; i++)
            chewer->freq_table[i] = 1; // frequency cannot be 0, as the algorithm won't work - adding two roots 0 + 0 is not greater than 0
    }
}

void huffman_chew_one_sample (Did did_i, STRp(sample), bool skip_if_same_as_prev)
{
    HuffmanChewerP chewer = B1ST (HuffmanChewer, ZCTX(did_i)->huffman);

    if (chewer->num_chewed >= MAX_CHEWED) return;

    if (skip_if_same_as_prev && str_issame(sample, chewer->prev_sample))
        return;

    uint32_t i=0;

    // ignore any master of sample that is the same chewer->master (up to max_prefix_len)
    while (i < MIN_(sample_len, chewer->max_prefix_len) && sample[i] == chewer->master[i]) 
        i++;
       
    for (; i < sample_len; i++) 
        chewer->freq_table[(uint8_t)sample[i]] += CHEW_UNIT; // much larger than the "1" default value

    if (skip_if_same_as_prev && sample_len <= HUFFMAN_MAX_PREV_LEN) {
        chewer->prev_sample_len = sample_len;
        memcpy (chewer->prev_sample, sample, sample_len);
    }

    chewer->num_chewed++;
}

uint32_t huffman_get_num_chewed (Did did_i)
{
    decl_zctx (did_i);

    return zctx->huffman.param ? B1ST (HuffmanCodes,  zctx->huffman)->num_chewed
                               : B1ST (HuffmanChewer, zctx->huffman)->num_chewed;
}

#define channel_i BNUM(zctx->huffman, h)
#define channel_str cond_str (zctx->huffman.len > 1, " channel=", ((rom[])NICO_CHAN_NAMES)[channel_i])

static void huffman_print_tree (HuffmanCodesP h, int node_i, int level)
{
    decl_zctx (h->did_i);
    rom space = "                                                                                                                                ";
    HuffmanNode *n = &h->nodes[node_i];
    
    if (!level) 
        iprintf ("%s%s: Huffman binary tree:\n", zctx->tag_name, channel_str);

    if (!n->left) // leaf
        iprintf ("%s: %.*s level=%u node_i=%u freq=%u P=%u c=%u \n", 
                 zctx->tag_name, level*2, space, level, node_i, n->freq, n->parent, n->c);

    else {
        iprintf ("%s: %.*s level=%u node_i=%u freq=%u P=%u L=%u R=%u\n", 
                 zctx->tag_name, level*2, space, level, node_i, n->freq, n->parent, n->left, n->right);

        huffman_print_tree (h, n->left,  level+1);
        huffman_print_tree (h, n->right, level+1);
    }
}

// convert nodes to a binary tree by successively linking two tree roots, to replace the two trees with one new tree.
static uint32_t huffman_generate_binary_tree (ContextP zctx, HuffmanChewerP chewer, HuffmanCodesP h, const bool mask[ALPHABET])
{
    // initialize
    uint32_t num_nodes=1; // we skip node 0, bc "0" means "NIL" 
    uint32_t num_roots=0, num_leaves=0;

    // initialize nodes from frequency table - all initial nodes are the leafs of the eventual binary tree
    for (int c=0; c < ALPHABET; c++) {
        if (mask[c]) {
            num_leaves++;
            h->roots[num_roots++] = num_nodes;
            h->nodes[num_nodes++] = (HuffmanNode){ .c = c, .freq = chewer->freq_table[c] };
        }
    }
    
    while (num_roots > 1) {

        // find the two lowest-frequency roots
        uint32_t min1_freq=0xffffffff, min2_freq=0xffffffff;
        int min1_i=0, min2_i=0; 
        
        for (int i=0; i < num_roots; i++) 
            // check if roots[i] is the lowest freq so far 
            if (h->nodes[h->roots[i]].freq < min1_freq) { // lowest so far
                min2_i = min1_i; // demote min1 to now be the second loweest
                min2_freq = min1_freq;

                min1_i = i;      // record new lowest
                min1_freq = h->nodes[h->roots[i]].freq;
            }

            // not the lowest, maybe the 2nd lowest?
            else if (h->nodes[h->roots[i]].freq < min2_freq) { // lowest so far
                min2_i = i;      // record new 2nd lowest
                min2_freq = h->nodes[h->roots[i]].freq;
            }

        // new root for the 2nd lowest frequency roots, which now cease to be roots
        h->nodes[num_nodes++] = (HuffmanNode){ .freq  = min1_freq + min2_freq,
                                               .left  = h->roots[min1_i],
                                               .right = h->roots[min2_i] };
        
        // add parent to both previous roots
        h->nodes[h->roots[min1_i]].parent = h->nodes[h->roots[min2_i]].parent = num_nodes - 1;

        // remove old roots from root array
        if (min2_i < min1_i) SWAP (min1_i, min2_i); // min1_i is now the smaller

        memmove (&h->roots[min1_i],   &h->roots[min1_i+1], (min2_i - min1_i - 1) * sizeof (h->roots[0]));    // move roots between min1 and min2 one back, overwriting min1
        memmove (&h->roots[min2_i-1], &h->roots[min2_i+1], (num_roots - min2_i - 1) * sizeof (h->roots[0])); // move roots after min2 two back - overwriting the empty space caused by the first move + min2
        num_roots -= 2;

        // add new root
        h->roots[num_roots++] = num_nodes - 1;
    }

    if (flag.debug_huffman_dict_id.num && flag.debug_huffman_dict_id.num == dict_id_typeless (zctx->dict_id).num)
        huffman_print_tree (h, h->roots[0], 0);

    return num_leaves;
}

static void huffman_show (HuffmanCodesP h)
{
    decl_zctx(h->did_i);

    iprintf ("\n%s %s Huffman codes:%s%s%s\n", IS_ZIP ? "ZIP" : "PIZ", zctx->tag_name, channel_str,
             cond_int (IS_ZIP, " num_chewed=", h->num_chewed), 
             cond_int (VER2(15,58), " longest_code_n_bits=", h->longest_code_n_bits));

    bool is_nico = zctx->huffman.len > 1; 

    for (int c=0; c < ALPHABET; c++) 
        if (h->code_n_bits[c]) {
            if (is_nico && channel_i >= NC_OP_FLANK) 
                iprintf ("%s: code of '%c'", zctx->tag_name, cigar_op_to_char[c]); // nico
            
            else if (!is_nico && IN_RANGX(c,' ', '~')) 
                iprintf ("%s: code of '%c'", zctx->tag_name, c);
            
            else                           
                iprintf ("%s: code of %-3d", zctx->tag_name, c);

            Bits bits = { .nwords = 1, .nbits = h->code_n_bits[c], .words = &h->codes[c], .type = BUF_REGULAR };
            bits_print (&bits);    
        }
}

static void huffman_set_longest_code_n_bits (HuffmanCodesP h)
{
    h->longest_code_n_bits = 0;
    for (int c=0; c < ALPHABET; c++)
        MAXIMIZE (h->longest_code_n_bits, h->code_n_bits[c]);
}

static void huffman_generate_codes (ContextP zctx, HuffmanCodesP h, uint32_t num_leaves)
{
    // initialize
    memset (h->codes, 0, sizeof (h->codes));
    memset (h->code_n_bits, 0, sizeof (h->code_n_bits));

    // generate codes - variable number of bits for c (i.e. each leaf of the binary tree)
    for (int i=1; i <= num_leaves; i++) { // stop at first non-leaf

        int c = h->nodes[i].c; // also == i-1

        Bits bits = { .nwords = 1, .nbits = 64, .words = &h->codes[c], .type = BUF_REGULAR };
        
        int node, parent;
        for (node = i; node != h->roots[0]; node = parent) {
            parent = h->nodes[node].parent;

            if (h->nodes[parent].right == node)      // node is its parent's left child
                bits_set (&bits, h->code_n_bits[c]); // note: stays 0 if left child

            h->code_n_bits[c]++;

            // edge case - code is more than 64 bits. TO DO: except for the highest 32 frequncies, set all others to equal frequency "1", and recalculate.
            if (h->code_n_bits[c] == 64 && parent) {
                huffman_print_tree (h, h->roots[0], 0); 
                ABORT ("While generating huffman codes of %s: code for c=%u exceeded 64 bits", zctx->tag_name, c);
            }
        }

        bits.nbits = h->code_n_bits[c];

        bits_reverse (&bits);
    }

    huffman_set_longest_code_n_bits (h);

    if (flag.show_huffman_dict_id.num && flag.show_huffman_dict_id.num == dict_id_typeless (zctx->dict_id).num) 
        huffman_show (h);
}

void huffman_produce_compressor (Did did_i, const HuffmanMask *mask) // caller guarantees that all chewed and future data is in this range - huffman does not test this
{
    decl_zctx(did_i);

    BufferP huff = &zctx->huffman;
    if (!buf_is_alloc (huff)) return; // context not encountered in segconf (possibly header-only file)

    // we're going to convert huff from HuffmanChewer to HuffmanCodes - save chewer data
    HuffmanChewer chewer[huff->len];
    memcpy (chewer, B1ST (HuffmanChewer, *huff), huff->len * sizeof (HuffmanChewer));

    buf_zero (huff);

    // we store the "longest_code_n_bits" across all channels in channel 0 
    uint8_t *longest_code_n_bits = &B1ST(HuffmanCodes, *huff)->longest_code_n_bits;

    for_buf2 (HuffmanCodes, h, i, *huff) {
        h->did_i = did_i;
        
        if (i == 0) { // only first chewer can have master (since it is only used for QNAME)
            h->max_prefix_len = chewer[i].max_prefix_len;
            h->master_len     = chewer[i].master_len;
            memcpy (h->master, chewer[i].master, chewer[i].master_len);
        }

        uint32_t num_leaves = huffman_generate_binary_tree (zctx, &chewer[i], h, mask[i]);

        h->num_chewed = chewer->num_chewed; // for show-huffman

        huffman_generate_codes (zctx, h, num_leaves);

        // zero entries not needed any more, to reduce the size of the SEC_HUFFMAN section
        memset (&h->roots[1], 0, (ALPHABET-1) * sizeof (h->roots[0]));
        for (int i=0; i < ARRAY_LEN(h->nodes); i++) {
            h->nodes[i].freq   = 0;
            h->nodes[i].parent = 0;
        }

        uint8_t this_longest = h->longest_code_n_bits;
        h->longest_code_n_bits = 0; // 0 to all channels but first, to improve section compression 
        h->num_chewed = 0;

        MAXIMIZE (*longest_code_n_bits, this_longest);
    }

    __atomic_thread_fence (__ATOMIC_ACQ_REL); // make sure codes are visible to all threads before setting param

    store_release (huff->param, IS_ZIP ? HUFF_PRODUCED_BY_ZIP : HUFF_PRODUCED_BY_PIZ); // indicate that huffman exists and is produced
}

void huffman_get_master (Did did_i, pSTRp(master))
{
    decl_zctx (did_i);
    ASSHUFFEXISTS;

    HuffmanCodesP h = B1ST(HuffmanCodes, zctx->huffman);
    STRset (*master, h->master);
} 

// maximum comp_len for first (or only) channel
uint32_t huffman_get_theoretical_max_comp_len (Did did_i, uint32_t uncomp_len) 
{
    decl_zctx (did_i);
    ASSHUFFEXISTS;
    
    return roundup_bits2bytes64 (uncomp_len * B1ST(HuffmanCodes, zctx->huffman)->longest_code_n_bits);
}

// note: sag_load expects compressed length in piz to be identical to zip, so changing anything
// in this compression will require supporting the old one as well in sag_load for back comp
uint32_t huffman_compress (VBlockP vb, Did did_i, STRp(uncomp), 
                           qSTR8p(comp)) // must be allocated by caller to (at least) maximum theoretical length 
{
    decl_zctx (did_i);
    ASSHUFFEXISTS;

    ARRAY (HuffmanCodes, h, zctx->huffman);

    Bits bits;
    uint64_t next_bit;
    int shift = huffman_compress_init_bits (comp, *comp_len, &bits, &next_bit);

    // ignore any master of sample that is the same chewer->master - just store its length (up to max_prefix_len, which is not more than 31)
    if (h->max_prefix_len) { // this implies single channel
        uint32_t actual_max_prefix_len = MIN_(uncomp_len, h->max_prefix_len);
        uint32_t i=0;
        while (i < actual_max_prefix_len && uncomp[i] == h->master[i]) 
            i++;

        bits_set_wordn (&bits, next_bit, i, HUFFMAN_PREFIX_LEN_BITS); // i is prefix_len
        next_bit += HUFFMAN_PREFIX_LEN_BITS;

        if (VER2(15,68))
            for (;i < uncomp_len; i++) 
                huffman_compress_add_one ((uint8_t)uncomp[i], *comp_len);

        else // backcomp for QNAME 15.0.65-15.0.67
            for (;i < uncomp_len; i++) {
                uint8_t c = (i < h->master_len) ? ((uint8_t)uncomp[i] ^ (uint8_t)h->master[i])
                          : h->master_len       ? ((uint8_t)uncomp[i] ^ (uint8_t)h->master[h->master_len-1]) // xor with last master char as to hopefully bring the range of the uncomp characters beyond the master to the similar range as the characters xored with the master
                          :                       (uint8_t)uncomp[i];

                huffman_compress_add_one (c, *comp_len);
            }
    }

    // shortcut if no master
    else 
        for (uint32_t i=0 ;i < uncomp_len; i++) 
            huffman_compress_add_one ((uint8_t)uncomp[i], *comp_len);

    huffman_compress_finalize_bits (zctx, &bits, next_bit, shift, STRa(comp));

    if (flag.debug_huffman_dict_id.num && flag.debug_huffman_dict_id.num == dict_id_typeless (zctx->dict_id).num)
        iprintf ("%s: uncomp_len=%u comp_len=%u ratio=%.1f\n", zctx->tag_name, uncomp_len, *comp_len, (double)uncomp_len / (double)*comp_len);

    return *comp_len;
}

// get length (in bytes) of compressed data, without generating the compressed data itself
static uint32_t huffman_compressed_len_(HuffmanCodesP h, STRp(uncomp))
{
    uint32_t i=0;
    uint64_t count_bits = 0;

    // ignore any master of sample that is the same chewer->master - just store its length (up to max_prefix_len, which is not more than 31)
    if (h->max_prefix_len) {
        uint32_t actual_max_prefix_len = MIN_(uncomp_len, h->max_prefix_len);
        while (i < actual_max_prefix_len && uncomp[i] == h->master[i]) 
            i++;

        count_bits += HUFFMAN_PREFIX_LEN_BITS;

        bool new_codec = VER2(15,68);

        for (;i < uncomp_len; i++) {
            int c = new_codec           ? (uint8_t)uncomp[i]
                  // backcomp for QNAME 15.0.65-15.0.67
                  : (i < h->master_len) ? ((uint8_t)uncomp[i] ^ (uint8_t)h->master[i])
                  : h->master_len       ? ((uint8_t)uncomp[i] ^ (uint8_t)h->master[h->master_len-1]) // xor with last master char as to hopefully bring the range of the uncomp characters beyond the master to the similar range as the characters xored with the master
                  :                       (uint8_t)uncomp[i];

            count_bits += h->code_n_bits[c];
        }
    }

    // shortcut if single-channel no master
    else 
        for (;i < uncomp_len; i++) 
            count_bits += h->code_n_bits[(uint8_t)uncomp[i]];

    return roundup_bits2bytes (count_bits); // number of bytes
}

uint32_t huffman_compressed_len (Did did_i, STRp(uncomp))
{
    decl_zctx (did_i);
    ASSHUFFEXISTS;

    if (!uncomp_len) return 0; // quick short cut
    
    ASSERT0 (zctx->huffman.len == 1, "multi channel not supported yet");

    return huffman_compressed_len_(B1ST(HuffmanCodes, zctx->huffman), STRa(uncomp));
}

// returns comp_len
int huffman_uncompress (Did did_i, bytes comp, STRc(uncomp))
{
    decl_zctx (did_i);
    ASSHUFFEXISTS;

    HuffmanCodesP h = B1ST(HuffmanCodes, zctx->huffman);

    huffman_uncompress_init_bits;

    // if we use master, then first 5 bits indicate the subset of master used
    uint32_t prefix_len = 0;
    if (h->max_prefix_len) {
        prefix_len = bits_get5 (&bits, bit_i); // assumes HUFFMAN_PREFIX_LEN_BITS=5
        ASSERT (prefix_len <= h->master_len, "expected prefix_len=%u <= h->master_len=%u", prefix_len, h->master_len);

        memcpy (uncomp, h->master, prefix_len);
        bit_i += 5;
    }

    // new codec starting 15.0.68. Files up to 15.0.64 didn't have SEC_HUFFMAN and a trival new codec huffman was produced in huffman_piz_read_all
    bool new_codec = VER2(15,68) || !VER2(15,65); 

    for (int i=prefix_len; i < uncomp_len; i++) {
        uint8_t c = huffman_uncompress_get_c;
                    
        uncomp[i] = new_codec           ? c
                  // backcomp for QNAME 15.0.65-15.0.67
                  : (!h->master_len)    ? c
                  : (i < h->master_len) ? (c ^ (uint8_t)h->master[i])
                  :                       (c ^ (uint8_t)h->master[h->master_len-1]); // xor with last master char as to hopefully bring the range of the uncomp characters beyond the master to the similar range as the characters xored with the master                        
    }

    return huffman_uncompress_comp_len;
}

int RECONSTRUCT_huffman (VBlockP vb, Did did_i, uint32_t uncomp_len, bytes comp) 
{
    int comp_len = huffman_uncompress (did_i, comp, BAFTtxt, uncomp_len);
    Ltxt += uncomp_len;

    return comp_len;
}

// get length (in bytes) of compressed data, without generating the uncompressed data itself
int huffman_uncompress_len (Did did_i, bytes comp, uint32_t uncomp_len)
{
    decl_zctx(did_i);
    ASSHUFFEXISTS;

    HuffmanCodesP h = B1ST(HuffmanCodes, zctx->huffman);

    Bits bits = { .nbits  = 64000000000ULL, // some very large number
                  .nwords = 1000000000ULL,
                  .type   = BITS_STANDALONE,
                  .words  = (uint64_t *)ROUNDDOWN8 ((uint64_t)comp) }; // shift back to align with a 64b word boundary (counting on pointers being up to 64 bit)

    // note: if we shifted .words back, the initial value of bit_i will point to were the data actually starts
    uint64_t bit_i = (comp - (bytes)bits.words) * 8;

    // if we use master, then first 5 bits indicate the subset of master used
    uint32_t prefix_len = 0;
    if (h->max_prefix_len) {
        prefix_len = bits_get5 (&bits, bit_i); // assumes HUFFMAN_PREFIX_LEN_BITS=5
        ASSERT (prefix_len <= h->master_len, "expected prefix_len=%u <= h->master_len=%u", prefix_len, h->master_len);

        bit_i += 5;
    }

    for (int i=prefix_len; i < uncomp_len; i++) {
        int node_i = h->roots[0];
        while (h->nodes[node_i].left)  // has children
            node_i = bits_get (&bits, bit_i++) ? h->nodes[node_i].right : h->nodes[node_i].left;
    }

    return roundup_bits2bytes(bit_i) - (comp - (bytes)bits.words);
}

bool huffman_issame (Did did_i, bytes comp, uint32_t uncomp_len, STRp(uncomp_compare_to))
{
    if (uncomp_len != uncomp_compare_to_len) return false;

    char uncomp[uncomp_len];
    huffman_uncompress (did_i, comp, STRa(uncomp));

    return !memcmp (uncomp, uncomp_compare_to, uncomp_len);
}

static void bgen_huffman_codes (ContextP zctx)
{
    for_buf (HuffmanCodes, h, zctx->huffman) {
        for (int i=0; i < ARRAY_LEN(h->codes); i++)
            h->codes[i] = BGEN64 (h->codes[i]);

        h->did_i = IS_ZIP ? 0 : zctx->did_i; // recalculated in piz, as might be different than in ZIP
                          
        // LTEN due to back comp: in 15.0.65-67 left, right and roots[0] were not modified (not BGENed).
        for (int i=0; i < ARRAY_LEN(h->nodes); i++) { 
            h->nodes[i].left  = LTEN16 (h->nodes[i].left);
            h->nodes[i].right = LTEN16 (h->nodes[i].right);
        }

        h->roots[0] = LTEN32 (h->roots[0]);
    }
}

// prepare and output the SEC_HUFFMAN section
void huffman_compress_section (Did did_i)
{
    if (!huffman_exists (did_i)) return;

    decl_zctx(did_i);

    bgen_huffman_codes (zctx);

    zctx->huffman.len *= sizeof (HuffmanCodes);
    Codec codec = codec_assign_best_codec (evb, NULL, &zctx->huffman, SEC_HUFFMAN);

    SectionHeaderHuffman header = (SectionHeaderHuffman){ 
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .section_type          = SEC_HUFFMAN,
        .data_uncompressed_len = BGEN32 (zctx->huffman.len32),
        .codec                 = codec,
        .dict_id               = zctx->dict_id,   
        .vblock_i              = 0,
    };

    comp_compress (evb, zctx, &evb->z_data, &header, zctx->huffman.data, NO_CALLBACK, "SEC_HUFFMAN");

    zctx->huffman.len /= sizeof (HuffmanCodes); // restore
}

static Mutex backcomp_huffman_qual_mutex = {};

// in case SEC_HUFFMAN of QUAL is missing, generate it based on first qual string reconstructed
void huffman_piz_backcomp_produce_qual (STRp(qual))
{
    mutex_lock (backcomp_huffman_qual_mutex);
    if (huffman_exists (SAM_QUAL)) goto done; // another thread beat up to producing - all good

    huffman_start_chewing (SAM_QUAL, 0, 0, 0, 1);
    huffman_chew_one_sample (SAM_QUAL, STRa(qual), false);
    huffman_produce_compressor (SAM_QUAL, (HuffmanMask[1]){ {[33 ... 126] = true} });   

    done:
    mutex_unlock (backcomp_huffman_qual_mutex);
}

// PIZ: called by the main thread from piz_read_global_area
void huffman_piz_read_all (void)
{
    if (flag.reading_reference) return;

    if (VER2(15,65)) {
        Section sec = NULL;
        while (sections_next_sec (&sec, SEC_HUFFMAN)) {
            ContextP zctx = ctx_get_existing_zctx (sec->dict_id);

            zfile_get_global_section (SectionHeaderHuffman, sec, &zctx->huffman, "huffman");
            if (flag.only_headers || !zctx->huffman.len) continue; // only show headers, or section skipped

            zctx->huffman.len /= sizeof (HuffmanCodes);
            zctx->huffman.param = HUFF_PRODUCED_BY_ZIP; 

            bgen_huffman_codes (zctx);

            HuffmanCodes *h = B1ST(HuffmanCodes, zctx->huffman);

            if (!VER2(15,68))  // in files compressed since 15.0.68, ZIP calculates this 
                huffman_set_longest_code_n_bits (h); 

            if (flag.show_huffman_dict_id.num == dict_id_typeless (sec->dict_id).num) {
                huffman_show (h);
                if (is_genocat) exit_ok;
            }
        }
    }

    // generate trivial huffmans for solo fields for which they are missing 
    if (IS_SAG_SOLO) // note: Solo huffmans start to be available in 15.0.68
        sam_piz_produce_trivial_solo_huffmans();

    if ((z_sam_gencomp || flag.deep) && !huffman_exists (SAM_QNAME)) { // QNAME huffman available since 15.0.65
        huffman_start_chewing (SAM_QNAME, 0, 0, 0, 1);
        huffman_produce_compressor (SAM_QNAME, (HuffmanMask[1]){ {[32 ... 126] = true} });     
    }

    if (flag.deep && !VER2(15,69)) {
        buf_free (ZCTX(SAM_CIGAR)->huffman); // note: in 15.0.68 CIGAR huffman was a now-obsolete single-channel huffman
        nico_start_chewing (SAM_CIGAR);
        nico_produce_compressor (SAM_CIGAR, true);
    }

    if (IS_SAG_SA && !huffman_exists(OPTION_SA_CIGAR)) { 
        nico_start_chewing (OPTION_SA_CIGAR);
        nico_produce_compressor (OPTION_SA_CIGAR, false);        
    }

    if ((z_sam_gencomp || flag.deep) && !huffman_exists (SAM_QUAL)) {
        // we will generate the qual in huffman_piz_backcomp_produce_qual based on the first QUAL string reconstructed
        mutex_initialize (backcomp_huffman_qual_mutex);

        buf_set_promiscuous (&ZCTX(SAM_QUAL)->huffman, "contexts->huffman");
    }
}

//--------------------------------------------------------------------
// Nico - multi-channel huffman for in-memory CIGAR compression
//--------------------------------------------------------------------

static uint8_t op_channel_by_prev_op[16] = { 
    [BC_NONE]=NC_OP_FLANK, 
    [BC_S]=NC_OP_AFTER_SHNP, [BC_H]=NC_OP_AFTER_SHNP, [BC_N]=NC_OP_AFTER_SHNP, [BC_P]=NC_OP_AFTER_SHNP,
    [BC_M]=NC_OP_AFTER_ME,   [BC_E]=NC_OP_AFTER_ME,
    [BC_X]=NC_OP_AFTER_X,    [BC_D]=NC_OP_AFTER_D,    [BC_I]=NC_OP_AFTER_I };

static uint8_t n_channel_by_op[9] = { 
    [BC_M]=NC_N_MEX, [BC_E]=NC_N_MEX, [BC_X]=NC_N_MEX, 
    [BC_S]=NC_N_SH,  [BC_H]=NC_N_SH,
    [BC_D]=NC_N_ID,  [BC_I]=NC_N_ID,
    [BC_N]=NC_N_NP,  [BC_P]=NC_N_NP };

// main thread
void nico_start_chewing (Did did_i)
{
    huffman_start_chewing (did_i, 0, 0, 0, NUM_NICO_CHANNELS);
}

// any thread - but only one designated thread as there is no thread-saftey 
void nico_chew_one_cigar (Did did_i, BamCigarOpP cigar, uint32_t n_cigar_op)
{
    HuffmanChewerP chewer = B1ST (HuffmanChewer, ZCTX(did_i)->huffman);

    if (chewer->num_chewed >= MAX_CHEWED) return;

    chewer[NC_N_CIGAR_OP].freq_table[MIN_(n_cigar_op, 255)] += CHEW_UNIT; 
    BamCigarOpType prev_op = BC_NONE;

    for (uint32_t i=0; i < n_cigar_op; i++) {
        BamCigarOpType op = cigar[i].op;

        int op_channel = (i==n_cigar_op-1) ? NC_OP_FLANK : op_channel_by_prev_op[prev_op]; 
        chewer[op_channel].freq_table[op] += 10; // much larger than the "1" default value        
        chewer[n_channel_by_op[op]].freq_table[MIN_((uint32_t)cigar[i].n, 255)] += CHEW_UNIT;
        
        prev_op = op;
    }
    
    chewer[0].num_chewed++; // num_chewed is counted on the [0] 
}

void nico_chew_one_textual_cigar (VBlockP vb, Did did_i, STRp(textual_cigar))
{
    if (huffman_get_num_chewed (did_i) >= MAX_CHEWED) return; // we don't need more samples

    ASSERTNOTINUSE (vb->scratch);
    sam_cigar_textual_to_binary (vb, STRa(textual_cigar), &vb->scratch, "scratch");

    nico_chew_one_cigar (did_i, CIG(vb->scratch));

    buf_free (vb->scratch);
}

// main thread
void nico_produce_compressor (Did did_i, bool is_deep_cigar)
{
    if (is_deep_cigar)
        huffman_produce_compressor (did_i, (HuffmanMask[NUM_NICO_CHANNELS]){ [NC_OP_FLANK ... NC_OP_AFTER_I]      = { [BC_M ... BC_D]=1, [BC_S]=1 }, // we only need S for bamass with insert contexts 
                                                                             [NC_N_CIGAR_OP]                      = { [1    ... 255 ]=1 }, 
                                                                             [NC_N_MEX ... NC_N_ID]               = { [0    ... 255 ]=1 },
                                                                             [NC_N_NP]                            = { [0    ... 1   ]=1 } }); // N,P not use in deep cigars - trivial channel
    else
        huffman_produce_compressor (did_i, (HuffmanMask[NUM_NICO_CHANNELS]){ [NC_OP_FLANK]                        = { [BC_M ... BC_X]=1 }, // the first or last op can be any op 
                                                                             [NC_OP_AFTER_SHNP ... NC_OP_AFTER_I] = { [BC_M ... BC_S]=1, [BC_P ... BC_X]=1 }, // the non-first/last cannot be H (but can be S!) 
                                                                             [NC_N_CIGAR_OP]                      = { [1    ... 255 ]=1 }, 
                                                                             [NC_N_MEX ... NC_N_NP]               = { [0    ... 255 ]=1 } });
}

#define INLINE_N_CIGAR_OP_BITS 20 // our limit
#define INLINE_OP_N_BITS       28 // 28 bits is the size of n per BAM spec 

// any thread: compress one cigar with big numbers inline
uint32_t nico_compress_cigar (VBlockP vb, Did did_i, BamCigarOpP cigar, uint32_t n_cigar_op, STR8c(comp)) 
{
    decl_zctx (did_i);
    ASSHUFFEXISTS;

    ARRAY (HuffmanCodes, channels, zctx->huffman);

    ASSERT (n_cigar_op < 1 MB, "cannot inline-compress a CIGAR with n_cigar_op=%u >= %u", n_cigar_op, (uint32_t)(1 MB));

    Bits bits;
    uint64_t next_bit;
    int shift = huffman_compress_init_bits (STRa(comp), &bits, &next_bit);

    HuffmanCodesP h = &channels[NC_N_CIGAR_OP];
    huffman_compress_add_one(MIN_(n_cigar_op, 255), comp_len);
    if (n_cigar_op >= 255)
        huffman_compress_add_int (n_cigar_op, INLINE_N_CIGAR_OP_BITS, comp_len); 

    BamCigarOpType prev_op = BC_NONE;

    for (uint32_t i=0; i < n_cigar_op; i++) {
        BamCigarOpType op = cigar[i].op;
        uint32_t n = cigar[i].n;

        h = &channels[(i==n_cigar_op-1) ? NC_OP_FLANK : op_channel_by_prev_op[prev_op]];
        huffman_compress_add_one (op, comp_len);
        prev_op = op;

        h = &channels[n_channel_by_op[op]]; 
        huffman_compress_add_one (MIN_(n,255), comp_len);        

        if (n >= 255)
            huffman_compress_add_int (n, INLINE_OP_N_BITS, comp_len); 
    }

    huffman_compress_finalize_bits (zctx, &bits, next_bit, shift, qSTRa(comp));

    if (flag.debug_huffman_dict_id.num && flag.debug_huffman_dict_id.num == dict_id_typeless (zctx->dict_id).num)
        iprintf ("%s: n_cigar_op=%u comp_len=%u ratio=%.1f\n", 
                 zctx->tag_name, n_cigar_op, comp_len, ((double)(n_cigar_op*4)+2) / (double)comp_len);

    return comp_len;
}

uint32_t nico_compress_textual_cigar (VBlockP vb, Did did_i, STRp(textual_cigar), STR8c(comp)) 
{
    BufferP binary_cigar = &vb->codec_bufs[6]; 
    ASSERTNOTINUSE (*binary_cigar);

    sam_cigar_textual_to_binary (vb, STRa(textual_cigar), binary_cigar, "scratch");

    comp_len = nico_compress_cigar (vb, did_i, B1ST(BamCigarOp, *binary_cigar), binary_cigar->len32, comp, comp_len);

    buf_free (*binary_cigar);

    return comp_len;
}

static uint32_t nico_compressed_len (Did did_i, BamCigarOpP cigar, uint32_t n_cigar_op)
{
    decl_zctx(did_i);
    ASSHUFFEXISTS;

    ARRAY (HuffmanCodes, channels, zctx->huffman);

    ASSERT (n_cigar_op < 1 MB, "cannot inline-compress a CIGAR with n_cigar_op=%u >= %u", n_cigar_op, (uint32_t)(1 MB));

    // bits due to n_cigar_op
    HuffmanCodesP h = &channels[NC_N_CIGAR_OP];
    uint32_t n_bits = h->code_n_bits[MIN_(n_cigar_op, 255)];
    if (n_cigar_op >= 255) n_bits += INLINE_N_CIGAR_OP_BITS;

    BamCigarOpType prev_op = BC_NONE;

    for (uint32_t i=0; i < n_cigar_op; i++) {
        BamCigarOpType op = cigar[i].op;
        uint32_t n = cigar[i].n;

        // bits due to OP
        h = &channels[(i==n_cigar_op-1) ? NC_OP_FLANK : op_channel_by_prev_op[prev_op]];
        n_bits += h->code_n_bits[op];
        prev_op = op;

        // bits due to N
        h = &channels[n_channel_by_op[op]]; 
        n_bits += h->code_n_bits[MIN_(n, 255)];
        if (n >= 255) n_bits += INLINE_N_CIGAR_OP_BITS;
    }

    return roundup_bits2bytes (n_bits);
}

uint32_t nico_compressed_len_textual_cigar (VBlockP vb, Did did_i, STRp(textual_cigar)) 
{
    ASSERTNOTINUSE (vb->scratch);
    sam_cigar_textual_to_binary (vb, STRa(textual_cigar), &vb->scratch, "scratch");

    uint32_t comp_len = nico_compressed_len (did_i, CIG(vb->scratch));

    buf_free (vb->scratch);

    return comp_len;
}

#define nico_uncompress_cigar_get_int(x, big_num_bits)      \
    x = huffman_uncompress_get_c;                           \
    if (x == 255) {                                         \
        x = bits_get_wordn (&bits, bit_i, (big_num_bits));  \
        bit_i += (big_num_bits);                            \
    }

// loose upper bound on compression size of a cigar with this many ops
uint32_t nico_max_comp_len (Did did_i, uint32_t n_cigar_op)
{
    return huffman_get_theoretical_max_comp_len (did_i, 1 + 2 * n_cigar_op);
} 

uint32_t nico_uncompress_cigar (VBlockP vb, Did did_i, bytes comp, BufferP cigar, rom buf_name)
{
    decl_zctx(did_i);
    ASSHUFFEXISTS;

    ARRAY (HuffmanCodes, channels, zctx->huffman);

    huffman_uncompress_init_bits;

    HuffmanCodesP h = &channels[NC_N_CIGAR_OP]; 
    nico_uncompress_cigar_get_int (cigar->len32, INLINE_N_CIGAR_OP_BITS);

    buf_alloc_exact (vb, *cigar, cigar->len, BamCigarOp, buf_name);
    BamCigarOpType prev_op = BC_NONE;

    for_buf2 (BamCigarOp, op, i, *cigar) {        
        h = &channels[(i == cigar->len32-1) ? NC_OP_FLANK : op_channel_by_prev_op[prev_op]];
        op->op = huffman_uncompress_get_c;
        prev_op = op->op;

        h = &channels[n_channel_by_op[op->op]]; 
        nico_uncompress_cigar_get_int (op->n, INLINE_OP_N_BITS);
    }        

    return huffman_uncompress_comp_len;
}

uint32_t nico_uncompress_textual_cigar (Did did_i, bytes comp, BufferP textual_cigar, 
                                        uint32_t textual_len, bool do_htos)
{
    decl_zctx(did_i);
    ASSHUFFEXISTS;

    ASSERT (huffman_exists(did_i), "%s.huffman doesn't exist", zctx->tag_name);

    ARRAY (HuffmanCodes, channels, zctx->huffman);

    huffman_uncompress_init_bits;

    HuffmanCodesP h = &channels[NC_N_CIGAR_OP]; 
    uint32_t n_cigar_op;
    nico_uncompress_cigar_get_int (n_cigar_op, INLINE_N_CIGAR_OP_BITS);

    char *next = BAFTc(*textual_cigar); // appends textual cigar to this buffer 

    BamCigarOpType prev_op = BC_NONE;

    for (int i=0; i < n_cigar_op; i++) {
        h = &channels[(i == n_cigar_op-1) ? NC_OP_FLANK : op_channel_by_prev_op[prev_op]];
        BamCigarOpType op = huffman_uncompress_get_c;
        
        if (do_htos && op == BC_S) op = BC_H;

        prev_op = op;

        h = &channels[n_channel_by_op[op]]; 
        uint32_t n;
        nico_uncompress_cigar_get_int (n, INLINE_OP_N_BITS);

        next += str_int (n, next);
        *next++ = cigar_op_to_char[op];
    }        

    // sanity check. note: this is our only use for cigar.piz.len in the code 
    uint32_t uncomp_len = next - BAFTc(*textual_cigar);
    ASSERT (uncomp_len == textual_len, "unexpected uncomp_len=%u, expected textual_len=%u", uncomp_len, textual_len);

    textual_cigar->len32 += uncomp_len;

    return huffman_uncompress_comp_len;
}

//--------------------------------------------------
// Unit test
//--------------------------------------------------

#ifdef DEBUG
void huffman_unit_test (void)
{
    flag.debug_huffman_dict_id.num = _SAM_QNAME;

    rom training_set[] = {
        "A00925:74:H25J5DSXY:1:1645:20383:5603",
        "A00910:90:H2V3LDSXY:4:1577:7961:20134",
        "A00910:90:H2V3LDSXY:4:1577:7961:20134",
        "A00925:74:H25J5DSXY:3:2424:22688:36996",
        "A00966:52:H2FYVDSXY:1:2221:24099:25473",
        "A00910:90:H2V3LDSXY:4:2144:20654:1407",
        "A00966:41:H25GJDSXY:3:1578:32615:24940",
        "A00910:85:HYGWJDSXX:4:1257:2899:21543",
        "A00910:90:H2V3LDSXY:3:1137:15094:26162",       
    };

    rom validation_set[] = {
        "A00966:52:H2FYVDSXY:1:2520:24288:19633",
        "A00966:41:H25GJDSXY:1:1103:7048:1485",
        "A00966:41:H25GJDSXY:1:1103:7048:1485",
        "A00966:52:H2FYVDSXY:1:2221:24099:25473",
        "A00910:90:H2V3LDSXY:4:2144:20654:1407",
        "A00910:85:HYGWJDSXX:2:1201:28203:14606"       
    };

    if (!z_file) z_file = CALLOC (sizeof(File)); // huffman  functions assign to z_file->contexts

    huffman_start_chewing (SAM_QNAME, training_set[0], strlen(training_set[0]), 31, 1);

    for (int i=0; i < ARRAY_LEN(training_set); i++)
        huffman_chew_one_sample (SAM_QNAME, training_set[i], strlen (training_set[i]), true);

    huffman_produce_compressor (SAM_QNAME, (HuffmanMask[1]){ {[32 ... 126] = true} });

    for (int i=0; i < ARRAY_LEN(validation_set); i++) {
        int uncomp_len = strlen(validation_set[i]);
        STRli(comp, huffman_get_theoretical_max_comp_len (SAM_QNAME, uncomp_len)); 

        huffman_compress (evb, SAM_QNAME, validation_set[i], uncomp_len, (uint8_t *)comp, &comp_len);
    
        char uncomp[uncomp_len];
        huffman_uncompress (SAM_QNAME, (bytes)comp, STRa(uncomp));

        printf ("i=%-2d before=\"%.*s\"\tafter=\"%.*s\"\tsame=%u\n", i, uncomp_len, validation_set[i], uncomp_len, uncomp, 
                !memcmp (validation_set[i], uncomp, uncomp_len));
    }
}
#endif
