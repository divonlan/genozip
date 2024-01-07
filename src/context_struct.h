// ------------------------------------------------------------------
//   context_struct.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
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

typedef enum __attribute__ ((__packed__)) { IDT_UNKNOWN, IDT_ALPHA_INT, IDT_ALPHA_NUM, IDT_ALPHA_INT_DOT_INT, IDT_ALPHA_NUM_DOT_INT, IDT_OTHER } IdType;

typedef enum __attribute__ ((__packed__)) { DEP_L0, DEP_L1, DEP_L2, NUM_LOCAL_DEPENDENCY_LEVELS } LocalDepType;

typedef struct Context {
    // ----------------------------
    // common fields for ZIP & PIZ
    // ----------------------------
    char tag_name[MAX_TAG_LEN];// nul-terminated tag name 
    Did did_i;                 // the index of this ctx within the array vb->contexts. PIZ: if this context is an ALIAS_CTX, did_i contains the destination context did_i
    union {
    Did st_did_i;              // ZIP: in --stats, consolidate this context into st_did_i
    Did other_did_i;           // PIZ: cache the other context needed for reconstructing this one
    };
    Did dict_did_i;            // ZIP/PIZ: zctx only: normally ==did_i, but if context is a ALIAS_DICT, did_i of its destination (shared dictionary between otherwise independent contexts)
    LocalType ltype;           // LT_* - type of local data - included in the section header
    LocalType pair_ltype;      // LT_* - Used if this file is a PAIR_2 - type of local data of PAIR_1
    struct FlagsCtx flags;     // flags to be included in section header
    struct FlagsCtx pair_flags;// Used if this file is a PAIR_2 - contains ctx->flags of the PAIR_1
    struct FlagsDict dict_flags;  // ZIP: zctx only ; PIZ: flags included in Dictionary section header (v15)
    SectionType pair_assist_type; // PIZ FASTQ R2: SEC_LOCAL, SEC_B250 is pair-assist, SEC_NONE if not.
    B250Size b250_size;        // Max size of element in b250 data (PIZ and ZIP after generation) v14
    B250Size pair_b250_size;
    DictId dict_id;            // the dict_id of this context
    #define FIRST_BUFFER_IN_Context dict
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    #define CTX_TAG_B250  "contexts->b250"
    Buffer b250;               // ZIP: During Seg, .data contains 32b indices into context->nodes. In zip_generate_b250_section, 
                               //      the "node indices" are converted into "word indices" - indices into the future 
                               //      context->word_list, in base-250. the number of words is moved from .len to .count. 
                               // PIZ: .data contains the word indices (i.e. indices into word_list) in base-250
    #define CTX_TAG_LOCAL "contexts->local"
    Buffer local;              // VB: Data private to this VB that is not in the dictionary
                               // ZIP Z_FILE .len  # fields of this type segged in the file (for stats)
    Buffer b250R1;             // ZIP/PIZ: used by PAIR_2 FASTQ VBs (inc. in Deep SAM), for paired contexts: PAIR_1 b250 data from corresponding VB (in PIZ: only if CTX_PAIR_LOAD)
    
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
    Buffer piz_ctx_specific_buf;
    Buffer qname_nodes;        // PIZ: used in KRAKEN_QNAME
    Buffer cigar_anal_history; // PIZ: used in SAM_CIGAR - items of type CigarAnalItem
    Buffer line_sqbitmap;      // PIZ: used in SAM_SQBITMAP
    Buffer domq_denorm;        // PIZ SAM/BAM/FASTQ: DomQual codec denormalization table for contexts with QUAL data 
    Buffer piz_lookback_buf;   // PIZ: SAM: used by contexts with lookback 
    Buffer channel_data;       // PIZ: SAM: QUAL/OPTION_iq_Z/OPTION_dq_Z/OPTION_sq_Z : used by PACB codec
    Buffer homopolymer;        // PIZ: SAM: OPTION_tp_B_c
    Buffer columns_data;       // PIZ: VCF: FORMAT_GT_HT with CODEC_HAPMAT
    Buffer column_of_zeros;    // PIZ: VCF: FORMAT_GT_HT_INDEX with CODEC_HAPMAT
    Buffer one_array;          // PIZ: VCF: FORMAT_PBWT_RUNS with CODEC_HAPMAT
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
    Codec lcodec;              // codec used to compress local
    Codec bcodec;              // codec used to compress b250
    Codec dcodec;              // codec used to compress dict
    Codec lsubcodec_zip;       // zip to compress with this codec AFTER compressing with lcodec
    Codec lsubcodec_piz;       // piz to decompress with this codec, AFTER decompressing with lcodec
    Codec lcodec_non_inherited;// ZIP z_file: non-inherited lcodec - used only for submitting stats
    
    // ZIP-only instructions NOT written to the genozip file
    union {
    bool no_stons;             // ZIP: don't attempt to move singletons to local even if ltype=LT_SINGLETON or if local is not used for anything else
    bool is_ctx_alias;         // PIZ: context is an alias            
    };

    bool no_callback;          // don't use callback for compressing, despite it being defined
    bool local_is_lten;        // if true local data is LTEN, otherwise it is the machine (native) endianity
    bool local_param;          // copy local.param to SectionHeaderCtx
    bool no_vb1_sort;          // don't sort the dictionary in ctx_sort_dictionaries_vb_1
    bool no_drop_b250;         // ZIP: the b250 section cannot be optimized away in zip_generate_b250_section (eg if we need section header to carry a param)
    bool z_data_exists;        // ZIP/PIZ: z_file has SEC_DICT, SEC_B250 and/or SEC_LOCAL sections of this context (not necessarily loaded)
    bool local_always;         // always create a local section in zfile, even if it is empty 
    bool is_stats_parent;      // other contexts have this context in st_did_i
    bool counts_section;       // ZIP: output ctx->counts to SEC_COUNTS section for this context
    bool subdicts_section;     // ZIP: output ctx->subdicts to SEC_SUBDICTS section for this context
    bool line_is_luft_trans;   // Seg: true if current line, when reconstructed with --luft, should be translated with luft_trans (false if no
                               //      trans_luft exists for this context, or it doesn't trigger for this line, or line is already in LUFT coordinates)
    union {
    bool lcodec_hard_coded;    // ZIP: lcodec is hard-coded and should not be reassigned
    bool empty_lookup_ok;      // PIZ: 
    };
    LocalDepType local_dep;    // ZIP: this local is created when another local is compressed (each NONREF_X is created with NONREF is compressed) (value=0,1,2)
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
                               // PIZ: local section of this VB found in file (used to determined if pair-identical R1 section should be loaded)
    uint32_t local_in_z_len;   
    uint32_t b250_in_z;        // ZIP: index and len into z_data where b250 compressed data is
                               // PIZ: b250 section of this VB found in file (used to determined if pair-identical R1 section should be loaded)
    uint32_t b250_in_z_len;    

    union {
    Buffer local_hash;         // ZIP: hash table for entries added by this VB that are not yet in the global (until merge_number)
                               // obtained by hash function hash(snip) and contains indices into local_ents
    Buffer history;            // PIZ: used if FlagsCtx.store_per_line and also for lookback (for files compressed starting with v12.0.41) - contains an array of either int64_t (if STORE_INT) or HistoryWord
    };
    Buffer local_ents;         // ZIP: linked entries - pointed to from local_hash

    uint32_t local_hash_prime; // prime number - size of the core (without extensions) has table 
    int32_t num_new_entries_prev_merged_vb; // zctx: updated in every merge - how many new words happened in this VB
                               // vctx: copied from zctx during clone, and used to initialize the size of local_hash
                               //         0 means no VB merged yet with this. if a previous vb had 0 new words, it will still be 1.
    union {
    Buffer global_hash;        // ZIP: global hash table that is populated during merge in zctx and is overlayed to vctx during clone. contains indices into global_ents.
    Buffer per_line;           // PIZ: data copied from txt_data for fields with textual store_per_line, used in if the line was dropped
    };
    Buffer global_ents;        // ZIP: pointed to by global cache

    uint32_t global_hash_prime;// prime number - size of the core (without extensions) hash table 

    uint32_t merge_num;        // in vctx: the merge_num when global_hash was cloned. only entries with merge_num <= this number 
                               // are valid. other entries may be added by later merges and should be ignored.
                               // in zctx: incremented with every merge into this ctx.
    // the next 2 are used in merge to set the size of the global hash table, when the first vb to create a ctx does so
    uint32_t nodes_len_at_1_3, nodes_len_at_2_3;  // value of nodes->len after an estimated 1/3 + 2/3 of the lines have been segmented
    
    union {
    ConstContainerP parent_container; // PIZ: last container that invoked reconstruction of this context

    // ZIP: stats
    uint64_t txt_len;          // ZIP: number of characters in reconstructed text are accounted for by snips in this ctx (for stats), when file reconstructed in PRIMARY coordinates (i.e. PRIMARY reconstruction for regular VBs, LUFT reconstruction for ##luft_only VBs, and no reconstruction for ##primary_only VBs) (note: seg_seg_long_CIGAR assumes this is uint64_t)
    };
    uint64_t local_num_words;  // ZIP: number of words (segs) that went into local. If a field is segged into multiple contexts - this field is incremented in each of them. If the context also uses b250, this field is ignored by stats which uses count instead.

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
    
    // ------------------------------------------------------------------------------------------------
    // START: RECONSTRUCT STATE : copied in reconstruct_peek 
    #define reconstruct_state_start(ctx) ((char*)&(ctx)->last_value)
    #define reconstruct_state_size_formula  ((char*)(&evb->contexts[0].pair_b250_iter + 1) - (char*)(&evb->contexts[0].last_value))

    ValueType last_value;      // ZIP/PIZ: last value of this context (it can be a basis for a delta, used for BAM translation, and other uses)
    union {
    int64_t last_delta;        // ZIP/PIZ: last delta value calculated (always in PIZ, sometimes in ZIP)
    WordIndex last_con_wi;     // PIZ: word index of last container retrieved from this ctx
    };

    #define INVALID_LAST_TXT_INDEX ((uint32_t)-1)
    TxtWord last_txt;          // ZIP/PIZ: index/len into vb->txt_data of last seg/reconstruction (always in PIZ, sometimes in Seg) (introduced 10.0.5)

    #define LAST_LINE_I_INIT -0x7fffffff
    LineIType last_line_i;     // ZIP/PIZ: =vb->line_i this line, so far, generated a valid last_value that can be used by downstream fields 
                               //          =(-vb->line_i-1) means ctx encountered in this line (so far) but last_value was not set 
    int32_t last_sample_i;     // ZIP/PIZ: Current sample in VCF/FORMAT ; must be set to 0 if not VCF/FORMAT
   
    union { // 32 bit
        int32_t ctx_specific;
        uint32_t segconf_max;       // maximum value during segconf
        bool last_is_alt;           // CHROM (all DTs): ZIP: last CHROM was an alt
        bool last_is_new;           // SAM_QNAME:       ZIP: used in segconf.running
        bool has_len;               // ZIP: INFO_ANN subfields of cDNA, CDS, AA
        thool use_special_sf;       // ZIP: INFO_SF
        thool line_has_RGQ;         // ZIP/PIZ: FORMAT_RGQ : GVCF

        struct {                    // SAM_QUAL, SAM_CQUAL, OPTION_OQ_Z, FASTQ_QUAL: 
            bool longr_bins_calculated; // ZIP zctx: codec_longr: were LONGR bins calculated in segconf
            Codec qual_codec;       // ZIP zctx: for displaying in stats
        };
        struct ctx_tp {             // PIZ: OPTION_tp_B_ARR (15.0.10-15.0.27: OPTION_tp_B_c)
            #define TP_LEN_BITS 11
            #define HP_LEN_BITS 9
            uint32_t repeat    : TP_LEN_BITS; // repeat within tp:B
            uint32_t hp_start  : TP_LEN_BITS; // start of current condensed homopolymer (-1 if this is not a condensed homopolymer)
            uint32_t hp_len    : HP_LEN_BITS; // length of current homopolymer
            uint32_t condensed : 1;           // current homopolymer is condensed
        } tp;
        struct {                    // INFO_DP:
            int32_t by_format_dp        : 1;   // ZIP/PIZ: segged vs sum of FORMAT/DP
            int32_t reconstruct         : 1;   // PIZ: valid if by_format_dp=1: does this INFO/DP need to be reconstructed
            int32_t sum_format_dp       : 30;  // ZIP/PIZ: sum of FORMAT/DP of samples in this line ('.' counts as 0).
        } dp;
        struct {                    // INFO_QD:         ZIP/PIZ: 
            uint32_t sum_dp_with_dosage : 28;  // sum of FORMAT/DP of samples in this line and dosage >= 1
            uint32_t pred_type          : 4;   // predictor type;
        } qd;
        int32_t last_end_line_i;    // INFO_END:        PIZ: last line on which INFO/END was encountered 

        IdType id_type;             // ZIP: type of ID in fields segged with seg_id_field        
        enum   __attribute__ ((__packed__)) { PAIR1_ALIGNED_UNKNOWN=-1, PAIR1_NOT_ALIGNED=0, PAIR1_ALIGNED=1 } pair1_is_aligned;  // FASTQ_SQBITMAP:  PIZ: used when reconstructing pair-2
        struct __attribute__ ((__packed__)) { Ploidy gt_prev_ploidy, gt_actual_last_ploidy; char gt_prev_phase; }; // FORMAT_GT: ZIP/PIZ
        struct __attribute__ ((__packed__)) { enum __attribute__ ((__packed__)) { PS_NONE, PS_POS, PS_POS_REF_ALT, PS_UNKNOWN } ps_type; }; // FORMAT_PS, FORMAT_PID, FORMAT_IPSphased
    };

    bool last_encounter_was_reconstructed; // PIZ: only valid if ctx_encountered() is true. Means last encountered was also reconstructed.
    uint32_t next_local;       // PIZ: iterator on Context.local 
    SnipIterator iterator;     // PIZ: used to iterate on the context, reading one b250 word_index at a time
    SnipIterator pair_b250_iter; // PIZ: Iterator on pair, if it contains b250 data  <--- LAST in RECONSTRUCT START 
                               // ZIP FASTQ paired: iterating on pair-1 data while compressing pair-2
    // END: RECONSTRUCT STATE 
    // ----------------------------------------------------------------------------------------
    
    // ZIP/PIZ: context specific #1
    union {
    #define CTX_TAG_CON_CACHE "contexts->con_cache"        
    Buffer con_cache;          // PIZ: vctx: use by contexts that might have containers: Handled by container_reconstruct - an array of Container which includes the did_i. 
                               //      Each struct is truncated to used items, followed by prefixes. 
                               // ZIP: vctx: seg_array, sam_seg_array_field_get_con cache a container.
    Buffer ctx_cache;          // PIZ: vctx: used to cached Contexts of Multiplexers and other dict_id look ups
    Buffer packed;             // PIZ: vctx: used by contexts that compressed CODEC_ACTG               
    Buffer zip_ctx_specific_buf;
    Buffer subdicts;           // ZIP/PIZ: zctx: Used by contexts that set ctx->subdicts_section: QUAL with PACB codec, iq:Z
    Buffer value_to_bin;       // ZIP: Used by LONGR codec on *_DOMQRUNS contexts
    Buffer longr_state;        // ZIP: Used by LONGR codec on QUAL contexts
    Buffer chrom2ref_map;      // ZIP (vctx & zctx), PIZ(zctx): Used by CHROM and contexts with a dict alias to it. Mapping from user file chrom to alternate chrom in reference file (for ZIP-VB: new chroms in this VB) - incides match ctx->nodes
    Buffer qual_line;          // ZIP: used by DOMQ codec on *_DOMQRUNS contexts
    Buffer normalize_buf;      // ZIP: used by DOMQ codec on QUAL contexts
    Buffer interlaced;         // ZIP: SAM: used to interlace BD/BI and iq/dq/sq line data
    Buffer mi_history;         // ZIP: SAM: used by OPTION_MI_Z in Ultima
    Buffer info_items;         // ZIP: VCF: VCF_INFO
    Buffer sf_snip;            // ZIP/PIZ: VCF: INFO_SF
    };

    // ZIP/PIZ: context specific #2
    union {
    #define CTX_TAG_CON_LEN "contexts->con_len"
    Buffer con_len;            // PIZ vctx: use by contexts that might have containers: Array of uint16_t - length of item in cache
    Buffer localR1;            // ZIP/PIZ vctx: PAIR_2 FASTQ VBs (inc. in Deep SAM): for paired contexts: PAIR_1 local data from corresponding VB (in PIZ: only if fastq_use_pair_assisted). Note: contexts with containers are always no_stons, so they have no local - therefore no issue with union conflict.
    Buffer ol_chrom2ref_map;   // ZIP vctx: SAM/BAM/VCF/CHAIN: CHROM: mapping from user file chrom to alternate chrom in reference file (cloned) - indices match vb->contexts[CHROM].ol_nodes. New nodes are stored in ctx->chrom2ref_map.
    Buffer ref2chrom_map;      // ZIP zctx: SAM/BAM/VCF/CHAIN: CHROM: reverse mapping from ref_index to chrom, created by ref_compress_ref
    };

    #define CTX_TAG_CON_INDEX "contexts->con_index"
    Buffer con_index;          // PIZ: use by contexts that might have containers: Array of uint32_t - index into con_cache - Each item corresponds to word_index. 
                               // ZIP: used by: 1. seg_array_of_struct
} Context;

typedef struct {
    DictId dict_id;
    Did did_i;
} ContextIndex;

typedef Did DictIdtoDidMap[65536 * 2];
static inline void init_dict_id_to_did_map (DictIdtoDidMap d2d_map)
{
    memset (d2d_map, 0xff, sizeof (DictIdtoDidMap)); // DID_NONE
}

typedef Context ContextArray[MAX_DICTS];

