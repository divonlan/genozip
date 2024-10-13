// ------------------------------------------------------------------
//   huffman.c
//   Copyright (C) 2023-2024 Genozip Limited. Patent Pending.
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
    uint32_t unused1[2];            
    uint8_t max_prefix_len;        // subset of master for which we reasonably expect identical prefixes
    uint8_t master_len;
    char master[HUFFMAN_MAX_MASTER_LEN]; // used for finding an identical prefix, and XORing the rest
} HuffmanCodes;

// ZIP: used for calculating codes (=huffman representation of values) 
typedef struct HuffmanChewer { 
    uint32_t freq_table[ALPHABET]; // number of occurances of each chewed value
    uint8_t max_prefix_len;
    uint8_t master_len;
    char master[HUFFMAN_MAX_MASTER_LEN];
    #define HUFFMAN_MAX_PREV_LEN 255
    uint8_t prev_sample_len;
    char prev_sample[HUFFMAN_MAX_PREV_LEN];
} HuffmanChewer;

void huffman_start_chewing (Did did_i, STRp(master), uint8_t max_prefix_len)
{
    ASSERTINRANGX (max_prefix_len, 0, HUFFMAN_MAX_MAX_PREFIX_LEN);
    if (master) ASSERTNOTZERO (max_prefix_len);

    buf_alloc_exact_zero (evb, ZCTX(did_i)->huffman, MAX_(sizeof (HuffmanChewer), sizeof(HuffmanCodes)), char, "huffman"); // alloc enough for HuffmanCodes so we don't need to realloc later
    
    HuffmanChewer *chewer = B1ST (HuffmanChewer, ZCTX(did_i)->huffman);
    
    if (master_len) {
        chewer->master_len = MIN_(master_len, HUFFMAN_MAX_MASTER_LEN); 
        memcpy (chewer->master, master, chewer->master_len);

        chewer->max_prefix_len = MIN_(max_prefix_len, master_len);
    }

    for (int i=0; i < ALPHABET; i++)
        chewer->freq_table[i] = 1; // frequency cannot be 0, as the algorithm won't work - adding two roots 0 + 0 is not greater than 0
}

void huffman_chew_one_sample (Did did_i, STRp(sample), bool skip_if_same_as_prev)
{
    HuffmanChewer *chewer = B1ST (HuffmanChewer, ZCTX(did_i)->huffman);

    if (skip_if_same_as_prev && str_issame(sample, chewer->prev_sample))
        return;

    uint32_t i=0;

    // ignore any master of sample that is the same chewer->master (up to max_prefix_len)
    while (i < MIN_(sample_len, chewer->max_prefix_len) && sample[i] == chewer->master[i]) 
        i++;
       
    for (; i < sample_len; i++) 
        chewer->freq_table[(uint8_t)sample[i]] += 10; // much larger than the "1" default value

    if (skip_if_same_as_prev && sample_len <= HUFFMAN_MAX_PREV_LEN) {
        chewer->prev_sample_len = sample_len;
        memcpy (chewer->prev_sample, sample, sample_len);
    }
}

static void huffman_print_tree (HuffmanCodesP h, int node_i, int level)
{
    rom space = "                                                                                                                                ";
    HuffmanNode *n = &h->nodes[node_i];
    
    if (!level) iprintf ("%s: Huffman binary tree:\n", ZCTX(h->did_i)->tag_name);

    if (!n->left) // leaf
        iprintf ("%s: %.*s level=%u node_i=%u freq=%u P=%u c=%u \n", 
                 ZCTX(h->did_i)->tag_name, level*2, space, level, node_i, n->freq, n->parent, n->c);

    else {
        iprintf ("%s: %.*s level=%u node_i=%u freq=%u P=%u L=%u R=%u\n", 
                 ZCTX(h->did_i)->tag_name, level*2, space, level, node_i, n->freq, n->parent, n->left, n->right);

        huffman_print_tree (h, n->left,  level+1);
        huffman_print_tree (h, n->right, level+1);
    }
}

// convert nodes to a binary tree by successively linking two tree roots, to replace the two trees with one new tree.
static uint32_t huffman_generate_binary_tree (HuffmanChewerP chewer, HuffmanCodesP h, const bool mask[ALPHABET])
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

        memmove (&h->roots[min1_i],   &h->roots[min1_i+1], (min2_i - min1_i - 1) * sizeof (h->roots[0]));      // move roots between min1 and min2 one back, overwriting min1
        memmove (&h->roots[min2_i-1], &h->roots[min2_i+1], (num_roots - min2_i - 1) * sizeof (h->roots[0])); // move roots after min2 two back - overwriting the empty space caused by the first move + min2
        num_roots -= 2;

        // add new root
        h->roots[num_roots++] = num_nodes - 1;
    }

    if (flag.debug_huffman)
        huffman_print_tree (h, h->roots[0], 0);

    return num_leaves;
}

static void huffman_show_one (HuffmanCodesP h, int c)
{
    if (IN_RANGX(c,' ', '~')) iprintf ("%s: code of '%c'", ZCTX(h->did_i)->tag_name, c);
    else                      iprintf ("%s: code of %-3d", ZCTX(h->did_i)->tag_name, c);

    Bits bits = { .nwords = 1, .nbits = h->code_n_bits[c], .words = &h->codes[c], .type = BUF_REGULAR };
    bits_print (&bits);    
}

static void huffman_show (HuffmanCodesP h)
{
    iprintf ("\n%s Huffman codes: %s\n", ZCTX(h->did_i)->tag_name,
             cond_int (VER2(15,58), " longest_code_n_bits=", h->longest_code_n_bits));

    for (int c=0; c < ALPHABET; c++) 
        if (h->code_n_bits[c]) 
            huffman_show_one (h, c);
}

static void huffman_generate_codes (HuffmanCodesP h, uint32_t num_leaves)
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
                ABORT ("While generating huffman codes of %s: code for c=%u exceeded 64 bits", ZCTX(h->did_i)->tag_name, c);
            }
        }

        bits.nbits = h->code_n_bits[c];

        bits_reverse (&bits);
    }

    if (flag.debug_huffman || flag.show_huffman) 
        huffman_show (h);
}

static void huffman_set_longest_code_n_bits (HuffmanCodesP h)
{
    h->longest_code_n_bits = 0;
    for (int c=0; c < ALPHABET; c++)
        MAXIMIZE (h->longest_code_n_bits, h->code_n_bits[c]);
}

void huffman_produce_compressor (Did did_i, 
                                 const bool mask[ALPHABET]) // caller guarantees that all chewed and future data is in this range - huffman does not test this
{
    decl_zctx(did_i);

    BufferP huff = &zctx->huffman;
    if (!buf_is_alloc (huff)) return; // context not encountered in segconf (possibly header-only file)

    // we're going to convert huff from HuffmanChewer to HuffmanCodes - save chewer data
    HuffmanChewer chewer = *B1ST (HuffmanChewer, *huff); // make a copy

    buf_zero (huff);
    HuffmanCodesP h = B1ST(HuffmanCodes, *huff);
    h->did_i          = did_i;
    h->max_prefix_len = chewer.max_prefix_len;
    h->master_len     = chewer.master_len;
    memcpy (h->master, chewer.master, chewer.master_len);

    uint32_t num_leaves = huffman_generate_binary_tree (&chewer, h, mask);
    
    huffman_generate_codes (h, num_leaves);

    huff->len = 1; // indicate that huffman exists and is produced

    // zero entries not needed any more, to reduce the size of the SEC_HUFFMAN section
    memset (&h->roots[1], 0, (ALPHABET-1) * sizeof (h->roots[0]));
    for (int i=0; i < ARRAY_LEN(h->nodes); i++) {
        h->nodes[i].freq   = 0;
        h->nodes[i].parent = 0;
    }

    huffman_set_longest_code_n_bits (h);
}

void huffman_get_master (Did did_i, pSTRp(master))
{
    ASSERT (huffman_exists(did_i), "%s.huffman doesn't exist", ZCTX(did_i)->tag_name);

    HuffmanCodesP h = B1ST(HuffmanCodes, ZCTX(did_i)->huffman);
    STRset (*master, h->master);
} 

uint32_t huffman_get_theoretical_max_comp_len (Did did_i, uint32_t uncomp_len) 
{
    return roundup_bits2bytes64 (uncomp_len * B1ST(HuffmanCodes, ZCTX(did_i)->huffman)->longest_code_n_bits);
}

// note: sag_load expects compressed length in piz to be identical to zip, so changing anything
// in this compression will require supporting the old one as well in sag_load for back comp
void huffman_compress (Did did_i, STRp(uncomp), 
                       qSTR8p(comp)) // must be allocated by caller to (at least) maximum theoretical length 
{
    HuffmanCodesP h = B1ST(HuffmanCodes, ZCTX(did_i)->huffman);

    // shift back so Bits starts from full word
    int shift = (uint64_t)comp & 0b111; // bytes shift 

    Bits bits = { .nbits  = (*comp_len + shift) * 8,    // note: possibly only partially covering last word
                  .nwords = ROUNDUP8(*comp_len + shift) / 8,
                  .words  = (uint64_t *)(comp - shift), // shift to bring address to word boundary
                  .type   = BITS_STANDALONE };

    uint64_t next_bit = shift * 8;

    // ignore any master of sample that is the same chewer->master - just store its length (up to max_prefix_len, which is not more than 31)
    if (h->max_prefix_len) {
        uint32_t actual_max_prefix_len = MIN_(uncomp_len, h->max_prefix_len);
        uint32_t i=0;
        while (i < actual_max_prefix_len && uncomp[i] == h->master[i]) 
            i++;

        bits_set_wordn (&bits, next_bit, i, HUFFMAN_PREFIX_LEN_BITS); // i is prefix_len
        next_bit += HUFFMAN_PREFIX_LEN_BITS;

        bool new_codec = VER2(15,68);
        
        for (;i < uncomp_len; i++) {
            int c = new_codec           ? (uint8_t)uncomp[i]
                // backcomp for QNAME 15.0.65-15.0.67
                : (i < h->master_len)   ? ((uint8_t)uncomp[i] ^ (uint8_t)h->master[i])
                : h->master_len         ? ((uint8_t)uncomp[i] ^ (uint8_t)h->master[h->master_len-1]) // xor with last master char as to hopefully bring the range of the uncomp characters beyond the master to the similar range as the characters xored with the master
                :                         (uint8_t)uncomp[i];

            ASSERT (next_bit + h->code_n_bits[c] <= bits.nbits, "%s: comp_len=%u too short for compressing \"%.*s\"", 
                    ZCTX(did_i)->tag_name, *comp_len, STRf(uncomp));

            bits_set_wordn (&bits, next_bit, h->codes[c], h->code_n_bits[c]);
            next_bit += h->code_n_bits[c];
        }
    }

    // shortcut if no master
    else 
        for (uint32_t i=0 ;i < uncomp_len; i++) {
            int c = (uint8_t)uncomp[i];

            ASSERT (next_bit + h->code_n_bits[c] <= bits.nbits, "%s: comp_len=%d too short for compressing \"%.*s\"", 
                    ZCTX(did_i)->tag_name, *comp_len, STRf(uncomp));

            bits_set_wordn (&bits, next_bit, h->codes[c], h->code_n_bits[c]);

            next_bit += h->code_n_bits[c];
        }

    ASSERT (next_bit <= bits.nbits, "bits overflow: expecting next_bit=%"PRIu64" <= nbits=%"PRIu64" for %s",
            next_bit, bits.nbits, ZCTX(did_i)->tag_name);

    bits_truncate (&bits, next_bit);

    *comp_len = roundup_bits2bytes(next_bit) - shift; // number of bytes

    // clear unused bits in last byte
    int bits_in_last_byte = next_bit % 8;
    if (bits_in_last_byte) 
        comp[*comp_len-1] &= bitmask8(bits_in_last_byte);

    if (flag.debug_huffman)
        iprintf ("%s: uncomp_len=%u comp_len=%u ratio=%.1f\n", ZCTX(did_i)->tag_name, uncomp_len, *comp_len, (double)uncomp_len / (double)*comp_len);
}

uint32_t huffman_compress_or_copy (Did did_i, STRp(uncomp), 
                                   STR8c(comp)) // must be allocated by caller to (at least) maximum theoretical length 
{
    if (huffman_exists (did_i)) { 
        huffman_compress (SAM_CIGAR, STRa(uncomp), comp, &comp_len);
        return comp_len;
    }
    else {
        memcpy (comp, uncomp, uncomp_len);
        return uncomp_len;
    }

    return comp_len;
}

// get length (in bytes) of compressed data, without generating the compressed data itself
uint32_t huffman_compress_len (Did did_i, STRp(uncomp))
{
    if (!uncomp_len) return 0; // quick short cut
    
    HuffmanCodesP h = B1ST(HuffmanCodes, ZCTX(did_i)->huffman);

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

    // shortcut if no master
    else 
        for (;i < uncomp_len; i++) 
            count_bits += h->code_n_bits[(uint8_t)uncomp[i]];

    return roundup_bits2bytes (count_bits); // number of bytes
}

// returns comp_len
int huffman_uncompress (Did did_i, bytes comp, STRc(uncomp))
{
    ASSERT (huffman_exists(did_i), "%s.huffman doesn't exist", ZCTX(did_i)->tag_name);

    HuffmanCodesP h = B1ST(HuffmanCodes, ZCTX(did_i)->huffman);

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

        memcpy (uncomp, h->master, prefix_len);
        bit_i += 5;
    }

    bool new_codec = VER2(15,68);

    for (int i=prefix_len; i < uncomp_len; i++) {
        int node_i = h->roots[0];
        while (h->nodes[node_i].left)  // has children
            node_i = bits_get (&bits, bit_i++) ? h->nodes[node_i].right : h->nodes[node_i].left;
                    
        uncomp[i] = new_codec           ? h->nodes[node_i].c
                  // backcomp for QNAME 15.0.65-15.0.67
                  : (!h->master_len)    ? h->nodes[node_i].c
                  : (i < h->master_len) ? (h->nodes[node_i].c ^ (uint8_t)h->master[i])
                  :                       (h->nodes[node_i].c ^ (uint8_t)h->master[h->master_len-1]); // xor with last master char as to hopefully bring the range of the uncomp characters beyond the master to the similar range as the characters xored with the master                        
    }

    return roundup_bits2bytes(bit_i) - (comp - (bytes)bits.words);
}

// uncompress from huffman compression, also supporting pre-15.0.65 files where data is not compressed
int huffman_uncompress_or_copy (Did did_i, bytes comp, STRc(uncomp)) 
{
    if (huffman_exists (did_i)) 
        return huffman_uncompress (did_i, comp, STRa(uncomp));

    else { 
        memcpy (uncomp, comp, uncomp_len);
        return uncomp_len;
    }
}

int RECONSTRUCT_huffman_or_copy (VBlockP vb, Did did_i, uint32_t uncomp_len, bytes comp) 
{
    int comp_len = huffman_uncompress_or_copy (did_i, comp, BAFTtxt, uncomp_len);
    Ltxt += uncomp_len;

    return comp_len;
}

// get length (in bytes) of compressed data, without generating the uncompressed data itself
int huffman_uncompress_len (Did did_i, bytes comp, uint32_t uncomp_len)
{
    ASSERT (huffman_exists(did_i), "%s.huffman doesn't exist", ZCTX(did_i)->tag_name);

    HuffmanCodesP h = B1ST(HuffmanCodes, ZCTX(did_i)->huffman);

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
    HuffmanCodes *h = B1ST (HuffmanCodes, zctx->huffman);

    for (int i=0; i < ARRAY_LEN(h->codes); i++)
        h->codes[i] = BGEN64 (h->codes[i]);

    h->did_i = 0; // recalculated in piz, as might be different than in ZIP

    // LTEN due to back comp: in 15.0.65-67 left, right and roots[0] were not modified (not BGENed).
    for (int i=0; i < ARRAY_LEN(h->nodes); i++) { 
        h->nodes[i].left  = LTEN16 (h->nodes[i].left);
        h->nodes[i].right = LTEN16 (h->nodes[i].right);
    }

    h->roots[0] = LTEN32 (h->roots[0]);
}

// prepare and output the SEC_HUFFMAN section
void huffman_compress_section (Did did_i)
{
    if (!huffman_exists (did_i)) return;

    ContextP zctx = ZCTX(did_i);

    zctx->huffman.len = sizeof (HuffmanCodes);
    Codec codec = codec_assign_best_codec (evb, NULL, &zctx->huffman, SEC_HUFFMAN);

    bgen_huffman_codes (zctx);

    SectionHeaderHuffman header = (SectionHeaderHuffman){ 
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .section_type          = SEC_HUFFMAN,
        .data_uncompressed_len = BGEN32 (sizeof (HuffmanCodes)),
        .codec                 = codec,
        .dict_id               = zctx->dict_id,   
        .vblock_i              = 0,
    };

    comp_compress (evb, zctx, &evb->z_data, &header, zctx->huffman.data, NO_CALLBACK, "SEC_HUFFMAN");

    zctx->huffman.len = 1; // restore
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

            bgen_huffman_codes (zctx);
            zctx->huffman.len = 1;

            HuffmanCodes *h = B1ST(HuffmanCodes, zctx->huffman);

            h->did_i = zctx->did_i; // might be different did_i than during zip

            if (!VER2(15,68))  // in files compressed since 15.0.68, ZIP calculates this 
                huffman_set_longest_code_n_bits (h); 

            if (flag.show_huffman) 
                huffman_show (h);
        }
    }

    if (is_genocat && flag.show_huffman) exit_ok;
}

#ifdef DEBUG
void huffman_unit_test (void)
{
    flag.debug_huffman = true;

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

    huffman_start_chewing (SAM_QNAME, training_set[0], strlen(training_set[0]), 31);

    for (int i=0; i < ARRAY_LEN(training_set); i++)
        huffman_chew_one_sample (SAM_QNAME, training_set[i], strlen (training_set[i]), true);

    huffman_produce_compressor (SAM_QNAME, (bool[256]){[32 ... 126] = true});

    for (int i=0; i < ARRAY_LEN(validation_set); i++) {
        int uncomp_len = strlen(validation_set[i]);
        STRli(comp, huffman_get_theoretical_max_comp_len (SAM_QNAME, uncomp_len)); 

        huffman_compress (SAM_QNAME, validation_set[i], uncomp_len, (uint8_t *)comp, &comp_len);
    
        char uncomp[uncomp_len];
        huffman_uncompress (SAM_QNAME, (bytes)comp, STRa(uncomp));

        printf ("i=%-2d before=\"%.*s\"\tafter=\"%.*s\"\tsame=%u\n", i, uncomp_len, validation_set[i], uncomp_len, uncomp, 
                !memcmp (validation_set[i], uncomp, uncomp_len));
    }
}
#endif
