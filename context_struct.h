// ------------------------------------------------------------------
//   context_struct.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef CONTEXT_STRUCT_INCLUDED
#define CONTEXT_STRUCT_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "mutex.h"
#include "sections.h"

typedef struct { // initialize with ctx_init_iterator()
    const uint8_t *next_b250;  // Pointer into b250 of the next b250 to be read (must be initialized to NULL)
    WordIndex prev_word_index; // When decoding, if word_index==BASE250_ONE_UP, then make it prev_word_index+1 (must be initalized to -1)
} SnipIterator;

typedef struct Context {
    // ----------------------------
    // common fields for ZIP & PIZ
    // ----------------------------
    #define MAX_TAG_LEN 64     // including terminating nul (must be divisible by 8 for Tag struct)
    char tag_name[MAX_TAG_LEN];// nul-terminated tag name 
    DidIType did_i;            // the index of this ctx within the array vb->contexts
    DidIType st_did_i;         // in --stats, consolidate this context into st_did_i
    LocalType ltype;           // LT_* - type of local data - included in the section header
    struct FlagsCtx flags;     // flags to be included in section header
    struct FlagsCtx pair_flags;// Used if this file is a PAIR_2 - contains ctx->flags of the PAIR_1
    DictId dict_id;            // which dict_id is this MTF dealing with
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict

    #define num_ctx_words param // b250.param, when it contains b250 data, holds the number of words (len is the number of bytes)
    Buffer b250;               // ZIP: During Seg, .data contains 32b indices into context->nodes. In zip_generate_b250_section, 
                               //      the "node indices" are converted into "word indices" - indices into the future 
                               //      context->word_list, in base-250. the number of words is moved from .len to .param. 
                               // PIZ: .data contains the word indices (i.e. indices into word_list) in base-250
    Buffer local;              // VB: Data private to this VB that is not in the dictionary
    Buffer pair;               // Used if this file is a PAIR_2 - contains a copy of either b250 or local of the PAIR_1 (if inst.pair_b250 or inst.pair_local is set)
    
    int64_t compressor_time;   // Used when --show-time - time for compressing / decompressing this context

    // rollback point
    uint64_t rback_b250_len, rback_local_len, rback_txt_len; // ZIP: data to roll back the last seg
    uint32_t rback_num_singletons, rback_last_txt_index, rback_last_txt_len;
    LastValueType rback_last_value;
    int64_t rback_last_delta;
    
    // ----------------------------
    // ZIP only fields
    // ----------------------------
    Buffer ol_dict;            // ZIP VB: tab-delimited list of all unique snips - overlayed all previous VB dictionaries
                               // ZIP zfile: singletons are stored here
                               // PIZ: counts are read to here (from SEC_COUNTS) - aligned to the words in word_list/dict
    Buffer ol_nodes;           // ZIP array of CtxNode - overlayed all previous VB dictionaries. char/word indices are into ol_dict.
#define ston_nodes ol_nodes    // ZIP zfile: nodes of singletons
    Buffer nodes;              // ZIP: array of CtxNode - in this VB that don't exist in ol_nodes. char/word indices are into dict.
                               // PIZ: in kraken's KRAKEN_QNAME context - contains qname_nodes 
                               // ZIP->PIZ zctx.nodes.param is transferred via SectionHeaderCounts.nodes_param if counts_section=true
    Buffer counts;             // ZIP/PIZ: counts of snips (array of int64_t)
    
    // Seg: snip (in dictionary) and node_index the last non-empty ("" or NULL) snip evaluated
    const char *last_snip;     
    unsigned last_snip_len;
    WordIndex last_snip_ni;   

    int tag_i;                 // ZIP dual-coordinates, VB only: index into vb->tags for tag renaming 

    // settings
    Codec lcodec, bcodec;      // codec used to compress local and b250
    Codec lsubcodec_piz;       // piz to decompress with this codec, AFTER decompressing with lcodec

    // ZIP-only instructions NOT written to the genozip file
    bool no_stons;             // don't attempt to move singletons to local (singletons are never moved anyway if ltype!=LT_TEXT)
    bool pair_local;           // ZIP: this is the 2nd file of a pair - compare vs the first file, and set flags.paired in the header of SEC_LOCAL
                               // PIZ: pair local data is loaded to context.pair
    bool pair_b250;            // ZIP: this is the 2nd file of a pair - compare vs the first file, and set flags.paired in the header of SEC_B250
                               // PIZ: pair b250 data is loaded to context.pair
    bool stop_pairing;         // this is the 2nd file of a pair - don't use SNIP_PAIR_LOOKUP/DELTA anymore until the end of this VB
    bool no_callback;          // don't use LOCAL_GET_LINE_CALLBACK for compressing, despite it being defined
    bool local_param;          // copy local.param to SectionHeaderCtx
    bool no_vb1_sort;          // don't sort the dictionary in ctx_sort_dictionaries_vb_1
    bool no_all_the_same;      // the b250 section cannot be optimized away in no_all_the_same (eg if we need section header to carry a param)
    bool local_always;         // always create a local section in zfile, even if it is empty 
    bool dynamic_size_local;   // resize LT_UINT32 according to data during generate (also do BGEN)
    bool numeric_only;         // if both numeric_only and dynamic_size_local are set, 
    bool is_stats_parent;      // other contexts have this context in st_did_i
    bool counts_section;       // output a SEC_COUNTS section for this context
    bool line_is_luft_trans;   // Seg: true if current line, when reconstructed with --luft, should be translated with luft_trans (false if no
                               //      trans_luft exists for this context, or it doesn't trigger for this line, or line is already in LUFT coordinates)
    TranslatorId luft_trans;   // ZIP: Luft translator for the context, set at context init and immutable thereafter
    
    // ZIP only: hash stuff 
    Buffer local_hash;         // hash table for entries added by this VB that are not yet in the global (until merge_number)
                               // obtained by hash function hash(snip) and the rest of linked to them by linked list
    uint32_t local_hash_prime; // prime number - size of the core (without extensions) has table 
    int32_t num_new_entries_prev_merged_vb; // zctx: updated in every merge - how many new words happened in this VB
                               // vctx: copied from zctx during clone, and used to initialize the size of local_hash
                               //         0 means no VB merged yet with this. if a previous vb had 0 new words, it will still be 1.
    Buffer global_hash;        // global hash table that is populated during merge in zctx and is overlayed to vctx during clone.
    uint32_t global_hash_prime; // prime number - size of the core (without extensions) has table 

    uint32_t merge_num;        // in vctx: the merge_num when global_hash was cloned. only entries with merge_num <= this number 
                               // are valid. other entries may be added by later merges and should be ignored.
                               // in zctx: incremented with every merge into this ctx.
    // the next 2 are used in merge to set the size of the global hash table, when the first vb to create a ctx does so
    uint32_t nodes_len_at_1_3, nodes_len_at_2_3;  // value of nodes->len after an estimated 1/3 + 2/3 of the lines have been segmented
    
    // ZIP: stats
    uint64_t txt_len;          // How many characters in reconstructed text are accounted for by snips in this ctx (for stats), when file reconstructed in PRIMARY coordinates (i.e. PRIMARY reconstruction for regular VBs, LUFT reconstruction for ##luft_only VBs, and no reconstruction for ##primary_only VBs)
    uint32_t num_singletons;   // True singletons that appeared exactly once in the entire file

    // PIZ-only
    Buffer word_list;          // PIZ only: word list. an array of CtxWord - listing the snips in dictionary
    bool semaphore;            // valid within the context of reconstructing a single line. MUST be reset ahead of completing the line.

    // ----------------------------
    // ZIP in z_file only
    // ----------------------------
    uint32_t num_failed_singletons;// (for stats) Words that we wrote into local in one VB only to discover later that they're not a singleton, and wrote into the global dict too
    Mutex mutex;               // Context in z_file (only) is protected by a mutex 
    
    // ------------------------------------------------------------------------------------------------
    // START: RECONSTRUCT STATE : copied in reconstruct_peek 
    #define reconstruct_state_start(ctx) ((char*)&(ctx)->last_value)
    #define reconstruct_state_size(ctx)  ((char*)(&(ctx)->pair_b250_iter + 1) - (char*)(&(ctx)->last_value))

    LastValueType last_value;  // ZIP/PIZ: last value of this context (it can be a basis for a delta, used for BAM translation, and other uses)
    int64_t last_delta;        // last delta value calculated
    
    #define INVALID_LAST_TXT_INDEX ((uint32_t)-1)
    uint32_t last_txt_index;   // ZIP/PIZ: index into vb->txt_data of last seg/reconstruction (always in PIZ, sometimes in Seg) (introduced 10.0.5)
    uint32_t last_txt_len;     // ZIP/PIZ: length (in vb->txt_data) of last seg/reconstruction (always in PIZ, sometimes in Seg)

    #define LAST_LINE_I_INIT -0x7fffffff
    int32_t last_line_i;       // ZIP/PIZ: =vb->line_i this line, so far, generated a valid last_value that can be used by downstream fields 
                               //          =(-vb->line_i-1) means ctx encountered in this line (so far) but last_value was not set 
    int32_t ctx_specific;      // ZIP/PIZ: for context-specific usage 
    uint32_t next_local;       // PIZ: iterator on Context.local
    SnipIterator iterator;     // PIZ: used to iterate on the context, reading one b250 word_index at a time
    SnipIterator pair_b250_iter; // PIZ: Iterator on pair, if it contains b250 data  <--- LAST in RECONSTRUCT START
    // END: RECONSTRUCT STATE 
    // ----------------------------------------------------------------------------------------

    // Container cache 
    Buffer con_cache;          // PIZ: Handled by container_reconstruct - an array of Container which includes the did_i. Each struct is truncated to used items, followed by prefixes. 
                               // ZIP: Each context is free to use it on its own
    Buffer con_index;          // Array of uint32_t - PIZ: index into con_cache - Each item corresponds to word_index. ZIP: context-specific.
    Buffer con_len;            // Array of uint16_t - length of item in cache
    
} Context;

#endif
