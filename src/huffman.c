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

#define ALPHABET 256

typedef struct {
    uint32_t freq;
    uint8_t c;
    uint16_t parent, left, right;
} HuffmanNode;

typedef struct Huffman {
    uint32_t freq_table[ALPHABET];
    uint64_t codes[ALPHABET];
    uint8_t code_n_bits[ALPHABET];
    HuffmanNode nodes[2*ALPHABET]; // 256 leafs, 255 internal nodes, 1 nil node
    int roots[ALPHABET]; // indices into nodes
    uint32_t num_roots, num_nodes, num_leaves;
    int samples_chewed;
    bool have_codes;
} Huffman;

HuffmanP huffman_initialize (void)
{
    HuffmanP h = CALLOC (sizeof (Huffman));

    for (int i=0; i < ALPHABET; i++)
        h->freq_table[i] = 1; // frequency cannot be 0, as the algorithm won't work - adding two roots 0 + 0 is not greater than 0

    return h;
}

void huffman_destroy (HuffmanP *h)
{
    FREE (*h);
}

int huffman_chew_one_sample (HuffmanP h, bytes sample, uint32_t sample_len)
{
    ASSERTISZERO (h->have_codes);

    for (uint32_t i=0; i < sample_len; i++)
        h->freq_table[(int)sample[i]] += 10; // much larger than the "1" default value

    h->samples_chewed++;

    return h->samples_chewed;
}

static void huffman_print_tree (HuffmanP h, int node_i, int level)
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
static void huffman_generate_binary_tree (HuffmanP h)
{
    // initialize
    h->num_nodes = 1; // we skip node 0, bc "0" means "NIL" 
    h->num_roots = 0;
    h->num_leaves = ALPHABET;

    // initialize nodes from frequency table - all initial nodes are the leafs of the eventual binary tree
    for (int c=0; c < ALPHABET; c++) {
        h->roots[h->num_roots++] = h->num_nodes;
        h->nodes[h->num_nodes++] = (HuffmanNode){ .c = c, .freq = h->freq_table[c] };
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

static bool huffman_generate_codes (HuffmanP h)
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

static void huffman_flatten_frequencies (HuffmanP h)
{
    // TO DO
    ABORT0 ("Huffman produces codes with more than 64 bits");
}

void huffman_produce_compressor (HuffmanP h)
{
    ASSERTISZERO (h->have_codes);

    huffman_generate_binary_tree (h);
    
    if (!huffman_generate_codes (h)) {
        // edge case - code is more than 64 bits. Solution: except for the highest 32 frequncies, set all others to equal frequency "1", and recalculate.
        huffman_flatten_frequencies (h);

        huffman_generate_binary_tree (h);

        // now that all but 32 frequenies are equal, the tree should not be able to reach a height of 64
        ASSERT0 (huffman_generate_codes (h), "Failed to generate Huffman codes");
    }

    // We need to set have_codes only after we have them, so other threads can rely on it. 
    __atomic_store_n (&h->have_codes, (bool)true, __ATOMIC_RELEASE); 
}

bool huffman_is_produced (HuffmanP h)
{
    return h && h->have_codes;
}

void huffman_compress (HuffmanP h, 
                       bytes uncomp,  uint32_t uncomp_len, 
                       uint8_t *comp, uint32_t *comp_len) // must be allocated by caller to 8*uncomp_len (hypothetical worse case scenario is 64bit per character)
{
    Bits bits = { .nbits  = *comp_len * 8,
                  .nwords = *comp_len / 8,
                  .words  = (uint64_t *)comp,
                  .type   = BUF_REGULAR      };

    uint64_t next_bit = 0;
    for (uint32_t i = 0 ; i < uncomp_len; i++) {
        int c = uncomp[i];
        bits_set_wordn (&bits, next_bit, h->codes[c], h->code_n_bits[c]);
        next_bit += h->code_n_bits[c];
    }

    *comp_len = roundup_bits2bytes(next_bit); // number of bytes

    if (flag.debug_huffman)
        iprintf ("uncomp_len=%u comp_len=%u ratio=%.1f\n", uncomp_len, *comp_len, (double)uncomp_len / (double)*comp_len);
}

// returns comp_len
int huffman_decompress (HuffmanP h, bytes comp, uint8_t *uncomp, uint32_t uncomp_len)
{
    Bits bits = { .nbits  = 64000000000ULL, // some very large number
                  .nwords = 1000000000ULL,
                  .type   = BITS_STANDALONE,
                  .words  = (uint64_t *)ROUNDDOWN8 ((uint64_t)comp) }; // shift back to align with a 64b word boundary (counting on pointers being up to 64 bit)

    // note: if we shifted .words back, the initial value of bit_i will point to were the data actually starts
    uint64_t bit_i = (comp - (bytes)bits.words) * 8;

    for (int i=0; i < uncomp_len; i++) {
        int node_i = h->roots[0];
        while (h->nodes[node_i].left)  // has children
            node_i = bits_get (&bits, bit_i++) ? h->nodes[node_i].right : h->nodes[node_i].left;

        uncomp[i] = h->nodes[node_i].c;
    }

    return roundup_bits2bytes(bit_i) - (comp - (bytes)bits.words);
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

    HuffmanP h = huffman_initialize();

    for (int i=0; i < ARRAY_LEN(training_set); i++)
        huffman_chew_one_sample (h, (bytes)training_set[i], strlen (training_set[i]));

    huffman_produce_compressor (h);

    for (int i=0; i < ARRAY_LEN(validation_set); i++) {
        int uncomp_len = strlen(validation_set[i]);
        STRli(comp, huffman_comp_len_required_allocation(uncomp_len));

        huffman_compress (h, (bytes)validation_set[i], uncomp_len, (uint8_t *)comp, &comp_len);
    
        char decomp[uncomp_len];
        huffman_decompress (h, (bytes)comp, (uint8_t *)decomp, uncomp_len);

        printf ("i=%-2d before=\"%.*s\"\tafter=\"%.*s\"\tsame=%u\n", i, uncomp_len, validation_set[i], uncomp_len, decomp, 
                !memcmp (validation_set[i], decomp, uncomp_len));
    }
}
#endif
