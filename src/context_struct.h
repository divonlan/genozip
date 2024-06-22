// ------------------------------------------------------------------
//   context_struct.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "buffer.h"
#include "mutex.h"
#include "sections.h"

typedef struct { // initialize with ctx_init_iterator()
    bytes next_b250;           // Pointer into b250 of the next b250 to be read (must be initialized to NULL)
    WordIndex prev_word_index; // When decoding, if word_index==BASE250_ONE_UP, then make it prev_word_index+1 (must be initialized to -1)
} SnipIterator;

typedef enum { DYN_DEC, DYN_hex, DYN_HEX } DynType;

typedef packed_enum { IDT_UNKNOWN, IDT_ALPHA_INT, IDT_ALPHA_NUM, IDT_ALPHA_INT_DOT_INT, IDT_ALPHA_NUM_DOT_INT, IDT_OTHER } IdType;

typedef packed_enum { DEP_L0, DEP_L1, DEP_L2, NUM_LOCAL_DEPENDENCY_LEVELS } LocalDepType;

typedef packed_enum {  // PIZ: set by a SPECIAL function, as if there was a WORD_INDEX_MISSING b250
    SPEC_RES_OK,        
    SPEC_RES_IS_MISSING,  // value missing, as if there was a WORD_INDEX_MISSING b250
    SPEC_RES_DEFERRED     // not reconstructed - will be reconstructed from top level container callback
} SpecialResult;              

#define CTX_TAG_B250       "contexts->b250"
#define CTX_TAG_LOCAL      "contexts->local"
#define CTX_TAG_DICT       "contexts->dict"
#define CTX_TAG_COUNTS     "contexts->counts"
#define CTX_TAG_NODES      "contexts->nodes"
#define CTX_TAG_LOCAL_HASH "contexts->local_hash"

typedef struct Context {
    // ------ common fields for ZIP & PIZ ------ 
    char tag_name[MAX_TAG_LEN];// nul-terminated tag name 
    DictId dict_id;            // the dict_id of this context

    Did did_i;                 // the index of this ctx within the array vb->contexts. PIZ: if this context is an ALIAS_CTX, did_i contains the destination context did_i
    Did dict_did_i;            // ZIP/PIZ: zctx only: normally ==did_i, but if context is a ALIAS_DICT, did_i of its destination (shared dictionary between otherwise independent contexts)

    uint32_t local_in_z;       // ZIP: index and len into z_data where local compressed data is
                               // PIZ: local section of this VB found in file (used to determined if pair-identical R1 section should be loaded)
    uint32_t local_in_z_len;   
    uint32_t b250_in_z;        // ZIP: index and len into z_data where b250 compressed data is
                               // PIZ: b250 section of this VB found in file (used to determined if pair-identical R1 section should be loaded)
    uint32_t b250_in_z_len;    

    union {
    uint8_t dict_helper;       // ZIP zctx / PIZ zctx+vctx: context-specific value passed through SectionHeaderDictionary.dict_helper (since 15.0.42)
    uint8_t con_rep_special;   // ZIP/PIZ: zctx: SPECIAL for getting container repeats in case of CON_REPEATS_IS_SPECIAL
    };

    LocalType ltype;           // LT_* - type of local data - included in the section header
    LocalType pair_ltype;      // LT_* - Used if this file is a PAIR_2 - type of local data of PAIR_1
    struct FlagsCtx flags;     // flags to be included in section header
    struct FlagsCtx pair_flags;// Used if this file is a PAIR_2 - contains ctx->flags of the PAIR_1
    struct FlagsDict dict_flags;  // ZIP zctx ; PIZ: zctx+vctx . Tramsmiited via SectionFlags.dictinonary (v15)
    B250Size b250_size;        // Size type of element in b250 data (PIZ and ZIP after generation) v14
    B250Size pair_b250_size;
    Codec lcodec;              // ZIP/PIZ: vctx/zctx: codec used to compress local
    union {
    Codec lsubcodec_piz;       // ZIP/PIZ: vctx: piz to decompress with this codec, AFTER decompressing with lcodec
    Codec qual_codec;          // ZIP zctx: QUAL codec selected in codec_assign_best_qual_codec
    };
    bool z_data_exists;        // ZIP/PIZ: z_file has SEC_DICT, SEC_B250 and/or SEC_LOCAL sections of this context (not necessarily loaded)
    bool is_initialized;       // ZIP / PIZ: context-specific initialization has been done
    
    uint8_t nothing_char;      // ZIP/PIZ: if non-zero, if local integer == max_int (for its ltype), nothing_char will be reconstructed instead. In PIZ, 0xff means fallback to pre-15.0.39

    #define FIRST_BUFFER_IN_Context dict
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    Buffer b250;               // ZIP: During Seg, .data contains variable-length indices into context->nodes. .count contains the number 
                               //      of b250s (still >1 if multiple b250s collapse with all-the-same) and .len contains the length in bytes.
                               //      In b250_zip_generate_section the "node indices" are converted into "word indices" - indices into the future context->word_list.
                               // PIZ: .data contains the word indices (i.e. indices into word_list) in base-250
    Buffer local;              // ZIP/PIZ vctx: Data private to this VB that is not in the dictionary
                               // ZIP zctx - only .len - number of fields of this type segged in the file (for stats)
    Buffer b250R1;             // ZIP/PIZ: used by PAIR_2 FASTQ VBs (inc. in Deep SAM), for paired contexts: PAIR_1 b250 data from corresponding VB (in PIZ: only if CTX_PAIR_LOAD)    

    Buffer counts;             // ZIP/PIZ: counts of snips (VB:uint32_t, z_file:uint64_t)
                               // ZIP: counts.param is a context-specific global counter that gets accumulated in zctx during merge (e.g. OPTION_SA_CIGAR)

    #define CTX_TAG_CON_INDEX "contexts->con_index"
    Buffer con_index;          // PIZ: use by contexts that might have containers: Array of uint32_t - index into con_cache - Each item corresponds to word_index. 
                               // ZIP: used by: 1. seg_array_of_struct

    // ZIP/PIZ: context specific buffer #1
    union {
        // GENERAL
        #define CTX_TAG_CON_CACHE "contexts->con_cache"        
        Buffer con_cache;          // PIZ: vctx: use by contexts that might have containers: Handled by container_reconstruct - an array of Container which includes the did_i. 
                                   //      Each struct is truncated to used items, followed by prefixes. 
                                   // ZIP: vctx: seg_array, sam_seg_array_field_get_con cache a container.
        Buffer ctx_cache;          // PIZ: vctx: used to cached Contexts of Multiplexers and other dict_id look ups
        Buffer chrom2ref_map;      // ZIP (vctx & zctx), PIZ(zctx): Used by CHROM and contexts with a dict alias to it. Mapping from user file chrom to alternate chrom in reference file (for ZIP-VB: new chroms in this VB) - incides match ctx->nodes
        Buffer snip_cache;         // ZIP: vctx: used by contexts that call seg_delta_vs_other_local*

        // CODECs
        Buffer packed;             // PIZ: vctx: used by contexts that compressed CODEC_ACTG               
        Buffer subdicts;           // ZIP/PIZ: zctx: Used by contexts that set ctx->subdicts_section: QUAL with PACB codec, iq:Z
        Buffer value_to_bin;       // ZIP: Used by LONGR codec on *_DOMQRUNS contexts
        Buffer longr_state;        // ZIP: Used by LONGR codec on QUAL contexts
        Buffer qual_line;          // ZIP: used by DOMQ codec on *_DOMQRUNS contexts
        Buffer normalize_buf;      // ZIP: used by DOMQ codec on QUAL contexts
        // SAM/BAM
        Buffer qname_hash;         // ZIP: vctx: SAM_QNAME: each entry i contains a line number for which the hash(qname)=i (or -1). prm8[0] is log2(len) (i.e., the number of bits)
        Buffer interlaced;         // ZIP: SAM: used to interlace BD/BI and iq/dq/sq line data
        Buffer mi_history;         // ZIP: SAM: used by OPTION_MI_Z in Ultima
        Buffer XG;                 // ZIP/PIZ: OPTION_XG_Z in bsseeker2: XG:Z field with the underscores removed. ZIP: revcomped if FLAG.revcomp. PIZ: as reconstructed when peeking during XM:Z special recon
        // VCF
        Buffer format_mapper_buf;  // ZIP: vctx: VCF_SAMPLES: an array of type Container - one entry per entry in CTX(VCF_FORMAT)->nodes   
        Buffer last_format;        // ZIP: vctx: VCF_FORMAT: cache previous line's FORMAT string
        Buffer id_hash;            // ZIP: vctx: VCF_ID (BND mates): each entry i contains a line number for which the hash(ID₀)=i (or -1). prm8[0] is log2(len) (i.e., the number of bits)
        Buffer info_items;         // ZIP: vctx: VCF_INFO
        Buffer deferred_snip;      // ZIP/PIZ: VCF: snip of a field whose seg/recon is postponed to after samples: INFO_SF, INFO_DPB
    };

    // ZIP/PIZ: context specific #2
    union {
        // GENERAL
        Buffer ol_chrom2ref_map;   // ZIP: vctx: SAM/BAM/VCF: CHROM: mapping from user file chrom to alternate chrom in reference file (cloned) - indices match vb->contexts[CHROM].ol_nodes. New nodes are stored in ctx->chrom2ref_map.
        Buffer ref2chrom_map;      // ZIP: zctx: SAM/BAM/VCF: CHROM: reverse mapping from ref_index to chrom, created by ref_compress_ref
        Buffer con_len;            // PIZ: vctx: use by contexts that might have containers: Array of uint16_t - length of item in cache
        // FASTQ
        Buffer localR1;            // ZIP/PIZ vctx: PAIR_2 FASTQ VBs (inc. in Deep SAM): for paired contexts: PAIR_1 local data from corresponding VB (in PIZ: only if fastq_use_pair_assisted). Note: contexts with containers are always no_stons, so they have no local - therefore no issue with union conflict.
        // VCF
        Buffer format_contexts;    // ZIP: vctx: VCF_SAMPLES: an array of format_mapper_buf.len of ContextPBlock
        Buffer insertion;          // PIZ: vctx: INFO_SF: inserted INFO fields reconstructed after samples
    };
                
    // ------------------------------------------------------------------------------------------------
    // START: RECONSTRUCT STATE : copied in reconstruct_peek 
    #define reconstruct_state_start(ctx) ((char*)&(ctx)->last_value)
    #define reconstruct_state_size_formula  ((char*)(&evb->contexts[0].last_encounter_was_reconstructed + 1) - (char*)(&evb->contexts[0].last_value))

    ValueType last_value;          // ZIP/PIZ: last value of this context (it can be a basis for a delta, used for BAM translation, and other uses)
    union {
        int64_t last_delta;        // ZIP/PIZ: last delta value calculated (always in PIZ, sometimes in ZIP)
        WordIndex last_con_wi;     // PIZ: word index of last container retrieved from this ctx
    };

    #define INVALID_LAST_TXT_INDEX ((uint32_t)-1)
    TxtWord last_txt;              // ZIP/PIZ: index/len into vb->txt_data of last seg/reconstruction (always in PIZ, sometimes in Seg) (introduced 10.0.5)

    #define LAST_LINE_I_INIT -0x7fffffff
    LineIType last_line_i;         // ZIP/PIZ: =vb->line_i this line, so far, generated a valid last_value that can be used by downstream fields 
                                   //          =(-vb->line_i-1) means ctx encountered in this line (so far) but last_value was not set 
    int32_t last_sample_i;         // ZIP/PIZ: Current sample in VCF/FORMAT ; must be set to 0 if not VCF/FORMAT
    
    union { // 64 bit
        int64_t ctx_specific;
        uint32_t segconf_max;      // maximum value during segconf
        bool last_is_alt;          // CHROM (all DTs): ZIP: last CHROM was has an alternative name
        IdType id_type;            // ZIP: type of ID in fields segged with seg_id_field        

        // SAM / BAM
        bool last_is_new;          // SAM_QNAME:       ZIP: used in segconf.running
        thool XG_inc_S;            // ZIP: bsseeker2 OPTION_XG_Z: whether to include soft_clip[0]
        struct {                   // SAM_QUAL, SAM_CQUAL, OPTION_OQ_Z, FASTQ_QUAL: 
            bool longr_bins_calculated; // ZIP zctx: codec_longr: were LONGR bins calculated in segconf
        };
        struct ctx_tp {            // PIZ: OPTION_tp_B_ARR (15.0.10-15.0.27: OPTION_tp_B_c)
            #define TP_LEN_BITS 23
            #define HP_LEN_BITS 17
            uint64_t repeat    : TP_LEN_BITS; // repeat within tp:B
            uint64_t hp_start  : TP_LEN_BITS; // start of current condensed homopolymer (-1 if this is not a condensed homopolymer)
            uint64_t hp_len    : HP_LEN_BITS; // length of current homopolymer
            uint64_t condensed : 1;           // current homopolymer is condensed
        } tp;

        // VCF
        struct {                    // ZIP/PIZ: FORMAT_GT 
            Ploidy prev_ploidy, actual_last_ploidy; 
            char prev_phase; 
        } gt; 
        struct {                    // ZIP/PIZ: FORMAT_GT_HT
            uint32_t use_HT_matrix : 1;  // ZIP: GT data from this line is placed in the HT matrix
            uint32_t ht_per_line   : 31; // number of haplotypes (columns) in the matrix
            uint32_t HT_n_lines;         // number of lines included in the matrix (not included: lines with no GT, lines copied in vcf_seg_sv_SAMPLES)
        };
        PosType32 pos_last_value;   // PIZ: VCF_POS: value for rolling back last_value after INFO/END
        bool has_len;               // ZIP: INFO_ANN subfields of cDNA, CDS, AA
        // char deferred;              // PIZ: !=0 if reconstruction deferred. possibly a seg-passed parameter

        struct {                    // ZIP: INFO_SF
            uint32_t next;
            uint32_t SF_by_GT; 
        } sf;    
        thool line_has_RGQ;         // ZIP/PIZ: FORMAT_RGQ : GVCF
        struct {                    // 
            int32_t sum_format_dp;  //   ZIP/PIZ: sum of FORMAT/DP of samples in this line ('.' counts as 0).
        } dp;
        struct {
            uint32_t count_ht;      // ZIP/PIZ: INFO/AN: sum of non-. haplotypes in FORMAT/GT, used to calculate INFO/AN
        } an;
        struct {                    // ZIP/PIZ: INFO_QD:  
            uint32_t sum_dp_with_dosage; // sum of FORMAT/DP of samples in this line and dosage >= 1
            uint32_t pred_type;     // predictor type
        } qd;
        struct {                    // PIZ: VCF_QUAL
            bool by_GP;             // QUAL_BY_GP used for this line 
            uint8_t decimals;
            bool is_rounddown;
            bool truncate_trailing_zeros;     
        } QUAL;
        uint16_t sum_sb[4];         // ZIP/PIZ: FORMAT_SB: sum_sb[i] is the sum of SBᵢ across all samples, where SBᵢ is the i'th component of a bi-allelic FORMAT/SB. 
        int32_t last_end_line_i;    // PIZ: INFO_END: last line on which INFO/END was encountered 
        TxtWord predicted_RU;       // PIZ/ZIP: INFO_RU: pointer into REFALT field in this line in txt_data
        bool saggy_seq_needs_fq_reversal; // PIZ: SAM_SQBITMAP: true if saggy copied is a reverse         
        struct { packed_enum { PS_NONE, PS_POS, PS_POS_REF_ALT, PS_UNKNOWN } ps_type; }; // FORMAT_PS, FORMAT_PID, FORMAT_IPSphased
        ContextP other_ctx;         // ZIP: used by FORMAT/RO, FORMAT/AO

        // FASTQ
        packed_enum { PAIR1_ALIGNED_UNKNOWN=-1, PAIR1_NOT_ALIGNED=0, PAIR1_ALIGNED=1 } r1_is_aligned;  // FASTQ_SQBITMAP:  PIZ: used when reconstructing pair-2
        TxtWord last_line1;         // ZIP segconf QNAME: entire line1
    };

    SnipIterator iterator;     // PIZ: used to iterate on the ctx->b250, reading one b250 word_index at a time
    SnipIterator pair_b250_iter; // PIZ: Iterator on pair, if it contains b250 data 
                               // ZIP FASTQ paired: iterating on pair-1 data while compressing pair-2
    uint32_t next_local;       // PIZ: iterator on Context.local 
    bool last_encounter_was_reconstructed; // PIZ: only valid if ctx_encountered() is true. Means last encountered was also reconstructed.
    // END: RECONSTRUCT STATE 
    // ----------------------------------------------------------------------------------------
    
    union {
    
    // ------ ZIP-only fields - common zctx and vctx ------ 
    struct {
    Buffer nodes;              // ZIP: array of CtxNode - in this VB that don't exist in ol_nodes. char/word indices are into dict.
                               // ZIP->PIZ zctx.nodes.param is transferred via SectionHeaderCounts.nodes_param if counts_section=true
    Buffer global_hash;        // ZIP: zctx/vctx: global hash table that is populated during merge in zctx and is overlayed to vctx during clone. contains indices into global_ents.

    uint64_t txt_len;          // ZIP zctx/vctx: number of characters in reconstructed (possibly modified) text are accounted for by snips in this ctx (for stats)  (note: seg_seg_long_CIGAR assumes this is uint64_t)
    int64_t txt_shrinkage;     // ZIP zctx/vctx: number of characters removed from txt due to modifications (can be negative)
    uint64_t local_num_words;  // ZIP zctx/vctx: number of words (segs) that went into local. If a field is segged into multiple contexts - this field is incremented in each of them. If the context also uses b250, this field is ignored by stats which uses count instead.

    int32_t num_new_entries_prev_merged_vb; // zctx: updated in every merge - how many new words happened in this VB
                               // vctx: copied from zctx during clone, and used to initialize the size of local_hash
                               //         0 means no VB merged yet with this. if a previous vb had 0 new words, it will still be 1.

    Did st_did_i;              // ZIP: in --stats, consolidate this context into st_did_i

    Codec bcodec;              // ZIP zctx/vctx: codec used to compress b250
    Codec dcodec;              // ZIP zctx/vctx: codec used to compress dict

    bool no_stons;             // ZIP: zctx/vctx: don't attempt to move singletons to local even if ltype=LT_SINGLETON or if local is not used for anything else
    bool local_param;          // copy local.param to SectionHeaderCtx
    bool counts_section;       // ZIP: zctx/vctx: output ctx->counts to SEC_COUNTS section for this context
    bool subdicts_section;     // ZIP: zctx/vctx: output ctx->subdicts to SEC_SUBDICTS section for this context
    bool lcodec_hard_coded;    // ZIP: zctx/vctx: lcodec is hard-coded and should not be reassigned
    bool is_stats_parent;      // other contexts have this context in st_did_i
    StoreType seg_to_local;    // ZIP: zctx/vctx: seg_array: this Int/Float field should be segged to local 
    
    union {
        struct {
            packed_enum { NOT_IN_HEADER=0, NUMBER_R=-1, NUMBER_A=-2, NUMBER_G=-3, NUMBER_VAR=-4, NUMBER_LARGE=-5/*Number∉[1,127]*/ } Number; // contains a value 1->127 or one of enumerated values
            packed_enum { VCF_Unknown_Type, VCF_Float, VCF_Integer, VCF_Character, VCF_Flag, VCF_String } Type;
        } vcf;
    } header_info;

    union {

    // ------ ZIP-only fields - vctx only ------ 
    struct { 
    Buffer ol_dict;            // ZIP vctx: tab-delimited list of all unique snips - overlayed zctx->dict (i.e. all previous VB dictionaries)
    Buffer ol_nodes;           // ZIP vctx: array of CtxNode - overlayed all previous VB dictionaries. char/word indices are into ol_dict.
    Buffer local_hash;         // ZIP: vctx: hash table for entries added by this VB that are not yet in the global (until merge_number)
                               // obtained by hash function hash(snip) and contains indices into vctx->nodes
    Buffer zip_lookback_buf;   // ZIP vctx: lookback_buf for contexts that use lookback

    // rollback point - used for rolling back during Seg (64b fields first and 32b fields after)
    int64_t rback_id;          // ZIP: rollback data valid only if ctx->rback_id == vb->rback_id
    TxtWord rback_last_txt;
    ValueType rback_last_value;
    int64_t rback_last_delta, rback_ctx_spec_param;
    
    int64_t dyn_int_min, dyn_int_max; // ZIP vctx: if ltype=LT_DYN* - min and max values encountered in this VB so far

    STR (last_snip);           // Seg: snip (in dictionary) and node_index the last non-empty ("" or NULL) snip evaluated             

    uint32_t rback_b250_count, rback_local_num_words, rback_local_len, rback_nodes_len, rback_txt_len; // ZIP: data to roll back the last seg

    uint32_t nodes_len_at_1_3, nodes_len_at_2_3;  // used in merge to set the size of the global hash table, when the first vb to create a ctx does so: value of nodes->len after an estimated 1/3 + 2/3 of the lines have been segmented

    bool local_always;         // ZIP vctx: always create a local section in zfile, even if it is empty 
    bool local_is_lten;        // ZIP vctx: if true local data is LTEN, otherwise it is the machine (native) endianity
    uint8_t dyn_lt_order;      // ZIP vctx: if ltype=LT_DYN*, the current ltype order of the data in local
    bool dyn_transposed;       // ZIP vctx: matrix should be transposed, if possible
    bool no_drop_b250;         // ZIP: the b250 section cannot be optimized away in b250_zip_generate_section (eg if we need section header to carry a param)
    bool no_callback;          // ZIP vctx: don't use callback for compressing, despite it being defined
    LocalDepType local_dep;    // ZIP: this local is created when another local is compressed (each NONREF_X is created with NONREF is compressed) (value=0,1,2)
    bool local_compressed;     // ZIP: VB: local has been compressed
    bool b250_compressed;      // ZIP: VB: b250 has been compressed
    bool nodes_converted;      // ZIP vctx: nodes have been converted in ctx_merge_in_one_vctx from index/len to word_index
    bool dict_merged;          // ZIP vctx: dict has been merged into zctx
    int16_t sf_i;              // ZIP VCF sample fields: 0-based index of this context within the FORMAT of this line (only for fields defined in vcf.h); -1 if not context present in this line
    int tag_i;                 // ZIP VCF: dual-coordinates VB only: index into vb->tags for tag renaming 
    WordIndex last_snip_ni;    // ZIP
    }; // ------ End of ZIP-only fields - vctx only ------

    // ------ ZIP-only fields - zctx only ------ 
    struct { 
    Buffer ston_hash;          // ZIP zctx: hash table for global singletons - each entry is a head of linked-list - index into ston_ents
    Buffer ston_ents;          // ZIP zctx: ents of hash of singletons - of type LocalHashEnt. contains link lists for each hash entry - headed from ston_hash
    uint32_t num_failed_singletons;// zctx: (for stats) Words that we wrote into local in one VB only to discover later that they're not a singleton, and wrote into the global dict too
    Codec lcodec_non_inherited;// ZIP zctx: non-inherited lcodec - used only for submitting stats
    uint8_t lcodec_count, bcodec_count; // ZIP zctx --best: approximate number of VBs in a row that selected this codec
    bool dict_len_excessive;   // ZIP zctx: dict is very big, indicating an ineffecient segging of this context
    bool please_remove_dict;   // zctx: one or more of the VBs request NOT compressing this dict (will be dropped unless another VB insists on keeping it)
    bool rm_dict_all_the_same; // zctx: we can remove the dict bc it is all-the-same
    bool override_rm_dict_ats; // zctx: don't remove dict, even if rm_dict_all_the_same is set
    bool all_the_same_wi_is_set; // zctx: zctx->dict_flags.all_the_same_wi is set 
    int8_t vb_1_pending_merges;// ZIP zctx: count of vb=1 merges still pending for this context (>1 if it has aliases). Other VBs can merge only if this is 0.
    }; // ------ End of ZIP-only fields - zctx only -------
    };
    };

    // ------ PIZ-only fields ------ 
    struct {
    Buffer word_list;          // PIZ zctx: word list. an array of CtxWord - listing the snips in dictionary
    Buffer per_line;           // PIZ: data copied from txt_data for fields with textual store_per_line, used in if the line was dropped
    Buffer piz_word_list_hash; // PIZ: used by ctx_get_word_index_by_snip
    Buffer history;            // PIZ: used if FlagsCtx.store_per_line and also for lookback (for files compressed starting with v12.0.41) - contains an array of either int64_t (if STORE_INT) or HistoryWord

    // PIZ: context-specific buffer
    union {
        Buffer piz_ctx_specific_buf;
        Buffer cigar_anal_history; // PIZ: used in SAM_CIGAR - items of type CigarAnalItem
        Buffer line_sqbitmap;      // PIZ: used in SAM_SQBITMAP
        Buffer domq_denorm;        // PIZ SAM/BAM/FASTQ: DomQual codec denormalization table for contexts with QUAL data 
        Buffer piz_lookback_buf;   // PIZ: SAM: used by contexts with lookback 
        Buffer channel_data;       // PIZ: SAM: QUAL/OPTION_iq_Z/OPTION_dq_Z/OPTION_sq_Z : used by PACB codec
        Buffer homopolymer;        // PIZ: SAM: OPTION_tp_B_c
    };

    ConstContainerP curr_container;// PIZ: current container in this context currently in the stack. NULL if none.

    WordIndex last_wi;         // PIZ: last word_index retrieved from b250 
    LineIType recon_insertion; // PIZ VCF: for deferred fields (mostly INFO fields inserted after samples) - whether to reconstruct. if to reconstruct - set to vb->line_i+1. any other value means "don't reconstruct"
    Did other_did_i;           // PIZ: cache the other context needed for reconstructing this one
    SectionType pair_assist_type; // PIZ FASTQ R2: SEC_LOCAL, SEC_B250 is pair-assist, SEC_NONE if not.
    bool is_ctx_alias;         // PIZ: context is an alias            
    bool local_uncompressed;   // PIZ: VB: local has been uncompressed
    bool b250_uncompressed;    // PIZ: VB: b250 has been uncompressed
    bool empty_lookup_ok;      // PIZ: 
    bool is_loaded;            // PIZ: either dict or local or b250 are loaded (not skipped) so context can be reconstructed
    bool semaphore;            // valid within the context of reconstructing a single line. MUST be reset ahead of completing the line.
    SpecialResult special_res; // PIZ: set by a SPECIAL function in case of result for which the reconstructor needs to take further action
    }; // ------ End of PIZ-only fields ------- 
    };

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

