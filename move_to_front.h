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
#include "base250.h"
#include "sections.h"

// fake mtf index values that go into genotype_data after segregation if subfields are missing
#define WORD_INDEX_MAX_INDEX  0xfffffffcUL // the number just smaller than all the special values below
#define WORD_INDEX_ONE_UP     0xfffffffdUL // the value is the one more than the previous value
#define WORD_INDEX_EMPTY_SF   0xfffffffeUL // subfield is missing, terminating : present
#define WORD_INDEX_MISSING_SF 0xffffffffUL // subfield is missing at end of cell, no :

#define NIL -1
typedef struct mtfnode_ {
    uint32_t char_index;      // character index into dictionary array
    uint32_t snip_len;        // not including \t terminator present in dictionary array
    Base250 word_index;       // word index into dictionary 
} MtfNode;

typedef struct {
    uint32_t char_index;
    uint32_t snip_len;
} MtfWord;

typedef struct {        
    int32_t mtf_i;            // index into MtfContext.ol_mtf (if < ol_mtf.len) or MtfContext.mtf or NIL
    int32_t next;             // linked list - index into MtfContext.hash or NIL
} HashEnt;

// used by variant_block_i=1 to collect frequency statistics
typedef struct {
    int32_t mtf_i;             // index into MtfContext.mtf
    int32_t count;             // number of times this snip has been encoutered so far
} SorterEnt;

typedef struct {
    const uint8_t *next_b250;  // Pointer into b250 of the next b250 to be read (must be initialized to NULL)
    int32_t prev_word_index;   // When decoding, if word_index==BASE250_ONE_UP, then make it prev_word_index+1 (must be initalized to -1)
} SnipIterator;

typedef struct mtfcontext_ {
    unsigned did_i;            // the index of this ctx within the array vb->mtf_ctx
    DictIdType dict_id;        // ZIP & PIZ. which dict_id is this MTF dealing with
    SectionType b250_section_type; // section type where the the corresponding b250 data goes
    SectionType dict_section_type; // section type for dictionary statistics (this is a "fake" section type, not one that we write to disk)
    Buffer ol_dict;            // ZIP only. tab-delimited list of all unique snips - overlayed all previous VB dictionaries
    Buffer ol_mtf;             // ZIP only. MTF nodes - overlayed all previous VB dictionaries. char/word indeces are into ol_dict.
    Buffer dict;               // ZIP & PIZ. tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    Buffer mtf;                // ZIP only. array of MtfNode - in this VB that don't exist in ol_mtf. char/word indeces are into dict.
    Buffer hash;               // ZIP only. hash table of entries HashEnt - initialized as a copy of all previous VBs. For entries are
                               // obtained by hash function hash(snip) and the rest of linked to them by linked list
    Buffer sorter;             // ZIP only. used by the first VB only of ZIP to sort the dictionary - entries of SorterEnt
    Buffer word_list;          // PIZ only. word list. an array of MtfWord - referring to the snips in dictionary
    
    Buffer mtf_i;              // ZIP only: contains 32bit indeces into the ctx->mtf - this is an intermediate step before generating fields_sections_data 
    
    Buffer b250;               // ZIP & PIZ: The buffer of b250 data containing indeces (in b250) to word_list
    
    // PIZ only: these two fields are used to iterate on the context, reading one b250 word_index at a time
    SnipIterator iterator;
} MtfContext;

static inline void mtf_init_iterator (MtfContext *ctx) { ctx->iterator.next_b250 = NULL ; ctx->iterator.prev_word_index = -1; }

extern int32_t mtf_evaluate_snip (VariantBlockP vb, MtfContext *ctx, const char *snip, uint32_t snip_len, MtfNode **node, bool *is_new);
extern uint32_t mtf_get_next_snip (VariantBlockP vb, MtfContext *ctx, SnipIterator *override_iterator, const char **snip, uint32_t *snip_len, uint32_t vcf_line);
extern void mtf_clone_ctx (VariantBlockP vb);
extern MtfNode *mtf_node (const MtfContext *ctx, uint32_t mtf_i, const char **snip_in_dict /* optional */, uint32_t *snip_len /* optional */);
extern void mtf_merge_in_vb_ctx (VariantBlockP vb);
extern MtfContext *mtf_get_ctx_by_dict_id (MtfContext *mtf_ctx, unsigned *num_dict_ids, unsigned *num_subfields, DictIdType dict_id, SectionType dict_section_type);
extern int mtf_get_existing_did_i_by_dict_id (VariantBlockP vb, DictIdType dict_id);
extern void mtf_integrate_dictionary_fragment (VariantBlockP vb, char *data);
extern void mtf_overlay_dictionaries_to_vb (VariantBlockP vb);
extern void mtf_sort_dictionaries_vb_1(VariantBlockP vb);
extern void mtf_zero_all_sorters (VariantBlockP vb);

extern void mtf_initialize_mutex (FileP z_file, unsigned next_variant_i_to_merge);

extern void mtf_free_context (MtfContext *ctx);

#endif