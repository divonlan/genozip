// ------------------------------------------------------------------
//   move_to_front.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef MOVE_TO_FRONT_INCLUDED
#define MOVE_TO_FRONT_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "base250.h"
#include "dict_id.h"

// fake mtf index values that go into genotype_data after segregation if subfields are missing
#define SEG_MAX_INDEX  0xfffffffdUL // the number just smaller than all the special values below
#define SEG_EMPTY_SF   0xfffffffeUL // subfield is missing, terminating : present
#define SEG_MISSING_SF 0xffffffffUL // subfield is missing at end of cell, no :

#define NIL -1
typedef struct {
    uint32_t char_index;       // character index into dictionary array
    uint32_t snip_len;         // not including \t terminator present in dictionary array
    Base250 word_index;        // word index into dictionary to be written to gt section
} MtfNode;

typedef struct {        
    int32_t mtf_i;             // index into MtfContext.ol_mtf (if < ol_mtf.len) or MtfContext.mtf or NIL
    int32_t next;              // linked list - index into MtfContext.hash or NIL
} HashEnt;

// used by variant_block_i=1 to collect frequency statistics
typedef struct {
    int32_t mtf_i;             // index into MtfContext.mtf
    uint32_t count;            // number of times this snip has been encoutered so far
} SorterEnt;

typedef struct {
    DictIdType dict_id;        // which dict_id is this MTF dealing with

    Buffer ol_dict;            // tab-delimited list of all unique snips - overlayed all previous VB dictionaries
    Buffer ol_mtf;             // MTF nodes - overlayed all previous VB dictionaries. char/word indeces are into ol_dict.
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    Buffer mtf;                // array of MtfNode - in this VB that don't exist in ol_mtf. char/word indeces are into dict.
    Buffer hash;               // hash table of entries HashEnt - initialized as a copy of all previous VBs. For entries are
                               // obtained by hash function hash(snip) and the rest of linked to them by linked list
    Buffer sorter;             // used by the first VB only of ZIP to sort the dictionary - entries of SorterEnt
    Buffer word_list;          // word list - used for PIZ only
} MtfContext;

extern int32_t mtf_evaluate_snip (VariantBlockP vb, MtfContext *ctx, const char *snip, uint32_t snip_len, bool overlayable, MtfNode **node /* out */);
extern void mtf_get_snip_by_word_index (VariantBlockP vb, MtfContext *ctx, const uint8_t *word_index_base250, char **snip, uint32_t *snip_len);
extern void mtf_clone_ctx (VariantBlockP vb);
extern MtfNode *mtf_node (const MtfContext *ctx, uint32_t mtf_i, const char **snip_in_dict /* optional out */);
extern unsigned mtf_merge_in_vb_ctx (VariantBlockP vb);
extern unsigned mtf_get_did_i_by_dict_id (MtfContext *mtf_ctx, unsigned *num_dict_ids, unsigned *num_subfields, DictIdType dict_id);
extern void mtf_integrate_dictionary_fragment (VariantBlockP vb, char *data);
extern void mtf_overlay_dictionaries_to_vb (VariantBlockP vb);
extern void mtf_sort_dictionaries_vb_1(VariantBlockP vb);

extern void mtf_initialize_mutex (FileP z_file, unsigned next_variant_i_to_merge);

extern void mtf_free_context (MtfContext *ctx);
#ifdef DEBUG
extern void mtf_tree_test (const MtfContext *ctx);
#endif

#endif