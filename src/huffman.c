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
    uint32_t num_roots, num_nodes, num_leaves;
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
       
    for (; i < sample_len; i++) {
        int c = (i < chewer->master_len) ? ((uint8_t)sample[i] ^ (uint8_t)chewer->master[i])
              : chewer->master_len       ? ((uint8_t)sample[i] ^ (uint8_t)chewer->master[chewer->master_len-1]) // xor with last master char as to hopefully bring the range of the uncomp characters beyond the master to the similar range as the characters xored with the master
              :                            (uint8_t)sample[i];
        
        chewer->freq_table[c] += 10; // much larger than the "1" default value
    }

    if (skip_if_same_as_prev && sample_len <= HUFFMAN_MAX_PREV_LEN) {
        chewer->prev_sample_len = sample_len;
        memcpy (chewer->prev_sample, sample, sample_len);
    }
}

static void huffman_print_tree (HuffmanCodesP h, int node_i, int level)
{
    rom space = "                                                                                                                                ";
    HuffmanNode *n = &h->nodes[node_i];
    
    if (!level) printf ("Huffman binary tree:\n");

    if (!n->left) // leaf
        printf ("%.*s level=%u node_i=%u freq=%u P=%u c=%u \n", level*2, space, level, node_i, n->freq, n->parent, n->c);

    else {
        printf ("%.*s level=%u node_i=%u freq=%u P=%u L=%u R=%u\n", 
                level*2, space, level, node_i, n->freq, n->parent, n->left, n->right);

        huffman_print_tree (h, n->left,  level+1);
        huffman_print_tree (h, n->right, level+1);
    }
}

// convert nodes to a binary tree by successively linking two tree roots, to replace the two trees with one new tree.
static void huffman_generate_binary_tree (HuffmanChewerP chewer, HuffmanCodesP h)
{
    // initialize
    h->num_nodes = 1; // we skip node 0, bc "0" means "NIL" 
    h->num_roots = 0;
    h->num_leaves = ALPHABET;

    // initialize nodes from frequency table - all initial nodes are the leafs of the eventual binary tree
    for (int c=0; c < ALPHABET; c++) {
        h->roots[h->num_roots++] = h->num_nodes;
        h->nodes[h->num_nodes++] = (HuffmanNode){ .c = c, .freq = chewer->freq_table[c] };
    }
    
    while (h->num_roots > 1) {

        // find the two lowest-frequency roots
        uint32_t min1_freq=0xffffffff, min2_freq=0xffffffff;
        int min1_i=0, min2_i=0; 
        
        for (int i=0; i < h->num_roots; i++) 
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
        h->nodes[h->num_nodes++] = (HuffmanNode){ .freq  = min1_freq + min2_freq,
                                                  .left  = h->roots[min1_i],
                                                  .right = h->roots[min2_i] };
        
        // add parent to both previous roots
        h->nodes[h->roots[min1_i]].parent = h->nodes[h->roots[min2_i]].parent = h->num_nodes - 1;

        // remove old roots from root array
        if (min2_i < min1_i) SWAP (min1_i, min2_i); // min1_i is now the smaller

        memmove (&h->roots[min1_i],   &h->roots[min1_i+1], (min2_i - min1_i - 1) * sizeof (h->roots[0]));      // move roots between min1 and min2 one back, overwriting min1
        memmove (&h->roots[min2_i-1], &h->roots[min2_i+1], (h->num_roots - min2_i - 1) * sizeof (h->roots[0])); // move roots after min2 two back - overwriting the empty space caused by the first move + min2
        h->num_roots -= 2;

        // add new root
        h->roots[h->num_roots++] = h->num_nodes - 1;
    }

    if (flag.debug_huffman)
        huffman_print_tree (h, h->roots[0], 0);
}

static bool huffman_generate_codes (HuffmanCodesP h)
{
    // initialize
    memset (h->codes, 0, sizeof (h->codes));
    memset (h->code_n_bits, 0, sizeof (h->code_n_bits));

    if (flag.debug_huffman) 
        iprint0 ("Huffman codes:\n");

    // generate codes - variable number of bits for c (i.e. each leaf of the binary tree)
    for (int i=1; i <= h->num_leaves; i++) { // stop at first non-leaf

        int c = h->nodes[i].c; // also == i-1

        Bits bits = { .nwords = 1, .nbits = 64, .words = &h->codes[c], .type = BUF_REGULAR };
        
        int node, parent;
        for (node = i; node != h->roots[0]; node = parent) {
            parent = h->nodes[node].parent;

            if (h->nodes[parent].right == node)      // node is its parent's left chile
                bits_set (&bits, h->code_n_bits[c]); // note: stays 0 if left child

            h->code_n_bits[c]++;

            // edge case - code is more than 64 bits.
            if (h->code_n_bits[c] == 64 && parent) return false;
        }

        bits.nbits = h->code_n_bits[c];

        bits_reverse (&bits);

        if (flag.debug_huffman) {
            iprintf ("code of %-3d: ", c);
            bits_print (&bits);        
        }
    }

    return true;
}

static void huffman_flatten_frequencies (HuffmanCodesP h)
{
    // TO DO
    ABORT0 ("Huffman produces codes with more than 64 bits");
}

void huffman_produce_compressor (Did did_i)
{
    // we're going to convert huff from HuffmanChewer to HuffmanCodes - save chewer data
    BufferP huff = &ZCTX(did_i)->huffman;
    if (!buf_is_alloc (huff)) return; // context not encountered in segconf (possibly header-only file)

    HuffmanChewer chewer = *B1ST (HuffmanChewer, *huff); // make a copy

    buf_zero (huff);
    HuffmanCodesP h = B1ST(HuffmanCodes, *huff);
    h->max_prefix_len = chewer.max_prefix_len;
    h->master_len     = chewer.master_len;
    memcpy (h->master, chewer.master, chewer.master_len);

    huffman_generate_binary_tree (&chewer, h);
    
    if (!huffman_generate_codes (h)) {
        // edge case - code is more than 64 bits. Solution: except for the highest 32 frequncies, set all others to equal frequency "1", and recalculate.
        huffman_flatten_frequencies (h);

        huffman_generate_binary_tree (&chewer, h);

        // now that all but 32 frequenies are equal, the tree should not be able to reach a height of 64
        ASSERT0 (huffman_generate_codes (h), "Failed to generate Huffman codes");
    }

    huff->len = 1; // exists
}

void huffman_get_master (Did did_i, pSTRp(master))
{
    ASSERT (huffman_exists(did_i), "%s.huffman doesn't exist", ZCTX(did_i)->tag_name);

    HuffmanCodesP h = B1ST(HuffmanCodes, ZCTX(did_i)->huffman);
    STRset (*master, h->master);
} 

void huffman_compress (Did did_i, STRp(uncomp), 
                       qSTR8p(comp)) // must be allocated by caller to huffman_comp_len_required_allocation() 
{
    HuffmanCodesP h = B1ST(HuffmanCodes, ZCTX(did_i)->huffman);

    uint32_t i=0;

    Bits bits = { .nbits  = *comp_len * 8,
                  .nwords = *comp_len / 8,
                  .words  = (uint64_t *)comp,
                  .type   = BUF_REGULAR      };

    uint64_t next_bit = 0;

    // ignore any master of sample that is the same chewer->master - just store its length (up to max_prefix_len, which is not more than 31)
    if (h->max_prefix_len) {
        uint32_t actual_max_prefix_len = MIN_(uncomp_len, h->max_prefix_len);
        while (i < actual_max_prefix_len && uncomp[i] == h->master[i]) 
            i++;

        bits_set_wordn (&bits, 0, i, HUFFMAN_PREFIX_LEN_BITS); // i is prefix_len
        next_bit += HUFFMAN_PREFIX_LEN_BITS;
    }

    for (;i < uncomp_len; i++) {
        int c = (i < h->master_len) ? ((uint8_t)uncomp[i] ^ (uint8_t)h->master[i])
              : h->master_len       ? ((uint8_t)uncomp[i] ^ (uint8_t)h->master[h->master_len-1]) // xor with last master char as to hopefully bring the range of the uncomp characters beyond the master to the similar range as the characters xored with the master
              :                        (uint8_t)uncomp[i];

        bits_set_wordn (&bits, next_bit, h->codes[c], h->code_n_bits[c]);
        next_bit += h->code_n_bits[c];
    }
    
    *comp_len = roundup_bits2bytes(next_bit); // number of bytes

    // hygine: clear unused bits in last byte
    int bits_in_last_byte = next_bit % 8;
    if (bits_in_last_byte) 
        comp[*comp_len-1] &= bitmask8(bits_in_last_byte);

    if (flag.debug_huffman)
        iprintf ("uncomp_len=%u comp_len=%u ratio=%.1f\n", uncomp_len, *comp_len, (double)uncomp_len / (double)*comp_len);
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

    for (int i=prefix_len; i < uncomp_len; i++) {
        int node_i = h->roots[0];
        while (h->nodes[node_i].left)  // has children
            node_i = bits_get (&bits, bit_i++) ? h->nodes[node_i].right : h->nodes[node_i].left;

        uncomp[i] = (i < h->master_len) ? (h->nodes[node_i].c ^ (uint8_t)h->master[i])
                  : h->master_len       ? (h->nodes[node_i].c ^ (uint8_t)h->master[h->master_len-1]) // xor with last master char as to hopefully bring the range of the uncomp characters beyond the master to the similar range as the characters xored with the master
                  :                        h->nodes[node_i].c;
    }

    return roundup_bits2bytes(bit_i) - (comp - (bytes)bits.words);
}

// reconstructs from huffman compression, also supporting pre-15.0.65 files where data is not compressed
int RECONSTRUCT_huffman (VBlockP vb, Did did_i, uint32_t uncomp_len, bytes comp) 
{
    int ret; 

    if (huffman_exists (did_i)) // files since 15.0.65
        ret = huffman_uncompress ((did_i), (comp), BAFTtxt, (uncomp_len));

    else { // files up to 15.0.64
        memcpy (BAFTtxt, (comp), (uncomp_len));
        ret = (uncomp_len);
    }
    Ltxt += uncomp_len;
    return ret;
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
    uint64_t *codes = B1ST (HuffmanCodes, zctx->huffman)->codes;
    for (int i=0; i < ALPHABET; i++)
        codes[i] = BGEN64 (codes[i]);
}

// prepare and output the SEC_HUFFMAN section
void huffman_compress_section (Did did_i)
{
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
}

// PIZ: called by the main thread from piz_read_global_area
void huffman_piz_read_all (void)
{
    if (VER2(15,65)) {
        Section sec = NULL;
        while (sections_next_sec (&sec, SEC_HUFFMAN))  {
            ContextP zctx = ctx_get_existing_zctx (sec->dict_id);

            zfile_get_global_section (SectionHeaderHuffman, sec, &zctx->huffman, "huffman");
            if (flag.only_headers || !zctx->huffman.len) continue; // only show headers, or section skipped

            bgen_huffman_codes (zctx);
            zctx->huffman.len = 1;
        }
    }
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

    huffman_produce_compressor (SAM_QNAME);

    for (int i=0; i < ARRAY_LEN(validation_set); i++) {
        int uncomp_len = strlen(validation_set[i]);
        STRli(comp, huffman_comp_len_required_allocation(uncomp_len));

        huffman_compress (SAM_QNAME, validation_set[i], uncomp_len, (uint8_t *)comp, &comp_len);
    
        char uncomp[uncomp_len];
        huffman_uncompress (SAM_QNAME, (bytes)comp, STRa(uncomp));

        printf ("i=%-2d before=\"%.*s\"\tafter=\"%.*s\"\tsame=%u\n", i, uncomp_len, validation_set[i], uncomp_len, uncomp, 
                !memcmp (validation_set[i], uncomp, uncomp_len));
    }
}
#endif
