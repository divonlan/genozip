// ------------------------------------------------------------------
//   context_struct.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "buffer.h"
#include "mutex.h"
#include "sections.h"

typedef struct { // initialize with ctx_init_iterator()
    bytes next_b250;           // Pointer into b250 of the next b250 to be read (must be initialized to NULL)
    WordIndex prev_word_index; // When decoding, if word_index==BASE250_ONE_UP, then make it prev_word_index+1 (must be initialized to -1)
} SnipIterator;

typedef enum { DYN_DEC, DYN_hex, DYN_HEX } DynType;

typedef struct Context {
    // ----------------------------
    // common fields for ZIP & PIZ
    // ----------------------------
    #define MAX_TAG_LEN 64     // including terminating nul (must be divisible by 8 for Tag struct)
    char tag_name[MAX_TAG_LEN];// nul-terminated tag name 
    Did did_i;            // the index of this ctx within the array vb->contexts
    union {
    Did st_did_i;         // ZIP: in --stats, consolidate this context into st_did_i
    Did other_did_i;      // PIZ: cache the other context needed for reconstructing this one
    };
    LocalType ltype;           // LT_* - type of local data - included in the section header
    struct FlagsCtx flags;     // flags to be included in section header
    struct FlagsCtx pair_flags;// Used if this file is a PAIR_2 - contains ctx->flags of the PAIR_1
    B250Size b250_size;        // Max size of element in b250 data (PIZ and ZIP after generation) v14
    B250Size pair_b250_size;
    DictId dict_id;            // which dict_id is this MTF dealing with
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    Buffer b250;               // ZIP: During Seg, .data contains 32b indices into context->nodes. In zip_generate_b250_section, 
                               //      the "node indices" are converted into "word indices" - indices into the future 
                               //      context->word_list, in base-250. the number of words is moved from .len to .count. 
                               // PIZ: .data contains the word indices (i.e. indices into word_list) in base-250
    Buffer local;              // VB: Data private to this VB that is not in the dictionary
                               // ZIP Z_FILE .len  # fields of this type segged in the file (for stats)
    Buffer pair;               // Used if this file is a PAIR_2 - contains a copy of either b250 or local of the PAIR_1 (if pair_b250 or pair_local is set)
    
    int64_t compressor_time;   // Used when --show-time - time for compressing / decompressing this context

    // rollback point - used for rolling back during Seg
    int64_t rback_id;          // ZIP: rollback data valid only if ctx->rback_id == vb->rback_id
    TxtWord rback_last_txt;
    uint32_t rback_b250_len, rback_local_len, rback_nodes_len, rback_txt_len; // ZIP: data to roll back the last seg
    uint32_t rback_b250_count, rback_local_num_words;
    ValueType rback_last_value;// also used in PIZ for rolling back VCF_POS.last_value after INFO/END
    int64_t rback_last_delta, rback_ctx_spec_param;
    
    union {
    Buffer ol_dict;            // ZIP VB: tab-delimited list of all unique snips - overlayed zctx->dict (i.e. all previous VB dictionaries)
    Buffer stons;              // ZIP zfile: singletons are stored here
    };

    union {
    Buffer ol_nodes;           // ZIP array of CtxNode - overlayed all previous VB dictionaries. char/word indices are into ol_dict.
    Buffer ston_nodes;         // ZIP z_file: nodes of singletons

    // PIZ: context-specific buffer
    Buffer qname_nodes;        // PIZ: used in KRAKEN_QNAME
    Buffer cigar_anal_history; // PIZ: used in SAM_CIGAR - items of type CigarAnalItem
    Buffer line_sqbitmap;      // PIZ: used in SAM_SQBITMAP
    Buffer domq_denorm;        // PIZ SAM/BAM/FASTQ: DomQual codec denormalization table for contexts with QUAL data 
    Buffer piz_lookback_buf;   // PIZ: SAM: used by contexts with lookback 
    };

    union {
    Buffer nodes;              // ZIP: array of CtxNode - in this VB that don't exist in ol_nodes. char/word indices are into dict.
                               // ZIP->PIZ zctx.nodes.param is transferred via SectionHeaderCounts.nodes_param if counts_section=true
    Buffer piz_word_list_hash; // PIZ: used by ctx_get_word_index_by_snip
    };

    Buffer counts;             // ZIP/PIZ: counts of snips (VB:uint32_t, z_file:uint64_t)
                               // ZIP: counts.param is a context-specific global counter that gets accumulated in zctx during merge (e.g. OPTION_SA_CIGAR)
    
    // Seg: snip (in dictionary) and node_index the last non-empty ("" or NULL) snip evaluated
    rom last_snip;             
    unsigned last_snip_len;
    WordIndex last_snip_ni;   

    int tag_i;                 // ZIP dual-coordinates VB only: index into vb->tags for tag renaming 

    // codecs
    uint8_t lcodec_count, bcodec_count; // ZIP z_file, --best: approximate number of VBs in a row that selected this codec
    Codec lcodec;              // codec used to compress local and dict
    Codec bcodec;              // codec used to compress b250
    Codec dcodec;              // codec used to compress dict
    Codec lsubcodec_zip;       // zip to compress with this codec AFTER compressing with lcodec
    Codec lsubcodec_piz;       // piz to decompress with this codec, AFTER decompressing with lcodec
    Codec lcodec_non_inherited;// ZIP z_file: non-inherited lcodec - used only for submitting stats
    
    // ZIP-only instructions NOT written to the genozip file
    bool no_stons;             // ZIP: don't attempt to move singletons to local (singletons are never moved anyway if ltype!=LT_TEXT)
    bool pair_local;           // ZIP: this is the 2nd file of a pair - compare vs the first file, and set flags.paired in the header of SEC_LOCAL
                               // PIZ: pair local data is loaded to context.pair
    bool pair_b250;            // ZIP: this is the 2nd file of a pair - compare vs the first file, and set flags.paired in the header of SEC_B250
                               // PIZ: pair b250 data is loaded to context.pair
    bool no_callback;          // don't use callback for compressing, despite it being defined
    bool local_is_lten;        // if true local data is LTEN, otherwise it is the machine (native) endianity
    bool local_param;          // copy local.param to SectionHeaderCtx
    bool no_vb1_sort;          // don't sort the dictionary in ctx_sort_dictionaries_vb_1
    bool no_drop_b250;         // the b250 section cannot be optimized away in zip_generate_b250_section (eg if we need section header to carry a param)
    bool local_always;         // always create a local section in zfile, even if it is empty 
    bool is_stats_parent;      // other contexts have this context in st_did_i
    bool counts_section;       // output a SEC_COUNTS section for this context
    bool line_is_luft_trans;   // Seg: true if current line, when reconstructed with --luft, should be translated with luft_trans (false if no
                               //      trans_luft exists for this context, or it doesn't trigger for this line, or line is already in LUFT coordinates)
    bool lcodec_hard_coded;    // ZIP: lcodec is hard-coded and should not be reassigned

    enum __attribute__ ((__packed__)) { DEP_L0, DEP_L1, DEP_L2, NUM_LOCAL_DEPENDENCY_LEVELS } local_dep; // ZIP: this local is created when another local is compressed (each NONREF_X is created with NONREF is compressed) (value=0,1,2)
    bool is_loaded;            // PIZ: either dict or local or b250 are loaded (not skipped) so context can be reconstructed
    bool is_initialized;       // ZIP / PIZ: context-specific initialization has been done
    union {
    bool local_compressed;     // ZIP: VB: local has been compressed
    bool local_uncompressed;   // PIZ: VB: local has been uncompressed
    };
    union {
    bool b250_compressed;      // ZIP: VB: b250 has been compressed
    bool b250_uncompressed;    // PIZ: VB: b250 has been uncompressed
    };
    bool dict_merged;          // ZIP: VB: dict has been merged into zctx
    bool please_remove_dict;   // ZFILE: one or more of the VBs request NOT compressing this dict (will be dropped unless another VB insists on keeping it)
    bool dict_len_excessive;   // ZFILE: dict is very big, indicating an ineffecient segging of this context
    
    TranslatorId luft_trans;   // ZIP: VCF: Luft translator for the context, set at context init and immutable thereafter
    
    int16_t sf_i;              // ZIP VCF FORMAT fields: 0-based index of this context within the FORMAT of this line (only for fields defined in vcf.h); -1 if not context present in this line

    uint32_t local_in_z;       // ZIP: index and len into z_data where local compressed data is
    uint32_t local_in_z_len;   
    uint32_t b250_in_z;        // ZIP: index and len into z_data where b250 compressed data is
    uint32_t b250_in_z_len;    

    union {
    Buffer local_hash;         // ZIP: hash table for entries added by this VB that are not yet in the global (until merge_number)
                               // obtained by hash function hash(snip) and the rest of linked to them by linked list
    Buffer history;            // PIZ: used if FlagsCtx.store_per_line and also for lookback (for files compressed starting with v12.0.41) - contains an array of either int64_t (if STORE_INT) or HistoryWord
    };

    uint32_t local_hash_prime; // prime number - size of the core (without extensions) has table 
    int32_t num_new_entries_prev_merged_vb; // zctx: updated in every merge - how many new words happened in this VB
                               // vctx: copied from zctx during clone, and used to initialize the size of local_hash
                               //         0 means no VB merged yet with this. if a previous vb had 0 new words, it will still be 1.
    union {
    Buffer global_hash;        // ZIP: global hash table that is populated during merge in zctx and is overlayed to vctx during clone.
    Buffer per_line;           // PIZ: data copied from txt_data for fields with textual store_per_line, used in if the line was dropped
    };

    uint32_t global_hash_prime;// prime number - size of the core (without extensions) hash table 

    uint32_t merge_num;        // in vctx: the merge_num when global_hash was cloned. only entries with merge_num <= this number 
                               // are valid. other entries may be added by later merges and should be ignored.
                               // in zctx: incremented with every merge into this ctx.
    // the next 2 are used in merge to set the size of the global hash table, when the first vb to create a ctx does so
    uint32_t nodes_len_at_1_3, nodes_len_at_2_3;  // value of nodes->len after an estimated 1/3 + 2/3 of the lines have been segmented
    
    union {
    ConstContainerP parent_container; // PIZ: last container that invoked reconstruction of this context

    // ZIP: stats
    uint64_t txt_len;          // How many characters in reconstructed text are accounted for by snips in this ctx (for stats), when file reconstructed in PRIMARY coordinates (i.e. PRIMARY reconstruction for regular VBs, LUFT reconstruction for ##luft_only VBs, and no reconstruction for ##primary_only VBs)
    };
    uint32_t local_num_words;  // ZIP: number of words (segs) that went into local. If a field is segged into multiple contexts - this field is incremented in each of them. If the context also uses b250, this field is ignored by stats which uses count instead.

    union {
    Buffer word_list;          // PIZ z_file: word list. an array of CtxWord - listing the snips in dictionary
    Buffer zip_lookback_buf;   // ZIP VB: lookback_buf for contexts that use lookback
    };

    bool semaphore;            // valid within the context of reconstructing a single line. MUST be reset ahead of completing the line.
    bool value_is_missing;     // PIZ: set by a SPECIAL function, as if there was a WORD_INDEX_MISSING b250
    
    // ----------------------------
    // ZIP in z_file only
    // ----------------------------
    uint32_t num_failed_singletons;// (for stats) Words that we wrote into local in one VB only to discover later that they're not a singleton, and wrote into the global dict too
    Mutex mutex;               // Context in z_file (only) is protected by a mutex 
    
    // ------------------------------------------------------------------------------------------------
    // START: RECONSTRUCT STATE : copied in reconstruct_peek 
    #define reconstruct_state_start(ctx) ((char*)&(ctx)->last_value)
    #define reconstruct_state_size_formula  ((char*)(&evb->contexts[0].pair_b250_iter + 1) - (char*)(&evb->contexts[0].last_value))

    ValueType last_value;      // ZIP/PIZ: last value of this context (it can be a basis for a delta, used for BAM translation, and other uses)
    int64_t last_delta;        // last delta value calculated
    
    #define INVALID_LAST_TXT_INDEX ((uint32_t)-1)
    TxtWord last_txt;          // ZIP/PIZ: index/len into vb->txt_data of last seg/reconstruction (always in PIZ, sometimes in Seg) (introduced 10.0.5)

    #define LAST_LINE_I_INIT -0x7fffffff
    LineIType last_line_i;     // ZIP/PIZ: =vb->line_i this line, so far, generated a valid last_value that can be used by downstream fields 
                               //          =(-vb->line_i-1) means ctx encountered in this line (so far) but last_value was not set 
    int32_t last_sample_i;     // ZIP/PIZ: Current sample in VCF/FORMAT ; must be set to 0 if not VCF/FORMAT
   
    union { // 32 bit
        int32_t ctx_specific;
        bool last_is_alt;         // CHROM (all DTs): ZIP: last CHROM was an alt
        bool last_is_new;         // SAM_QNAME:       ZIP: used in segconf.running
        int32_t sum_dp_this_line; // INFO_DP:         ZIP/PIZ: sum of FORMAT/DP of samples in this line ('.' counts as 0)
        int32_t last_end_line_i;  // INFO_END:        PIZ: last line on which INFO/END was encountered 
             
        enum   __attribute__ ((__packed__)) { PAIR1_ALIGNED_UNKNOWN=-1, PAIR1_NOT_ALIGNED=0, PAIR1_ALIGNED=1 } pair1_is_aligned;  // FASTQ_SQBITMAP:  PIZ: used when reconstructing pair-2
        struct __attribute__ ((__packed__)) { uint16_t gt_prev_ploidy; char gt_prev_phase; }; // FORMAT_GT: ZIP/PIZ
        struct __attribute__ ((__packed__)) { enum __attribute__ ((__packed__)) { PS_NONE, PS_POS, PS_POS_REF_ALT, PS_UNKNOWN } ps_type; }; // FORMAT_PS and FORMAT_PID
    };

    bool last_encounter_was_reconstructed; // PIZ: only valid if ctx_encountered() is true. Means last encountered was also reconstructed.
    uint32_t next_local;       // PIZ: iterator on Context.local 
    SnipIterator iterator;     // PIZ: used to iterate on the context, reading one b250 word_index at a time
    SnipIterator pair_b250_iter; // PIZ: Iterator on pair, if it contains b250 data  <--- LAST in RECONSTRUCT START 
                               // ZIP FASTQ paired: iterating on pair-1 data while compressing pair-2
    // END: RECONSTRUCT STATE 
    // ----------------------------------------------------------------------------------------

    // Container cache 
    Buffer con_cache;          // PIZ: Handled by container_reconstruct - an array of Container which includes the did_i. Each struct is truncated to used items, followed by prefixes. 
                               //      also used to cached Multiplexers in vcf_piz_special_MUX_BY_DOSAGE , piz_special_MINUS
                               // ZIP: Each context is free to use it on its own
    Buffer con_index;          // Array of uint32_t - PIZ: index into con_cache - Each item corresponds to word_index. 
                               // ZIP: used by: 1. seg_array_of_struct
    Buffer con_len;            // Array of uint16_t - length of item in cache
    
} Context;

typedef struct {
    DictId dict_id;
    Did did_i;
} ContextIndex;
