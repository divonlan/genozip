// ------------------------------------------------------------------
//   fastq_private.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "file.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "deep.h"
#include "aligner.h"
#include "sam_friend.h"

#define DTYPE_QNAME     DTYPE_1
#define DTYPE_FASTQ_AUX DTYPE_2

#define MAX_DESC_FIELDS (MAX_FIELDS-100)

typedef struct {
    TxtWord seq;
    TxtWord qual;               // start within vb->txt_data (qual.len==seq.len except if condensed by codec_homp_compress)
    uint32_t sam_seq_len;       // Deep / BamAss: seq_len of matching sam alignment
    uint32_t sam_seq_offset;    // Deep / BamAss: start of matching SAM SEQ within FASTQ SEQ
    bool dont_compress_QUAL;    // true in case of Deep and fully copied from SAM
    bool deep_qual;             // QUAL is to be fully or partially copied from Deep
    bool monochar;              // sequence is entirely of the same character (eg: NNNNN)
    bool qual_is_trims_only;    // Deep: fastq_zip_qual has modified the QUAL data to be trims only 
} ZipDataLineFASTQ, *ZipDataLineFASTQP;

// IMPORTANT: if changing fields in VBlockFASTQ, also update vb_fast_release_vb 
typedef struct VBlockFASTQ {
    VBLOCK_COMMON_FIELDS

    // --------- current line - reset before every line by fastq_reset_line() ----------
    #define first_fastq_vb_zero_per_line deep_seq_len
    uint32_t deep_seq_len;          // PIZ Deep
    uint32_t sam_seq_offset;        // PIZ Deep: offset of start of SEQ / QUAL copied from SAM with the FASTQ SEQ / QUAL
   
    PosType64 gpos;                 // ZIP bamass: gpos of lowest-gpos base consumed in the reference. note: if rev_comp, this correspond to the *last* ref-consumed base of the FASTQ SEQ 
    uint32_t ref_consumed;          // ZIP bamass
    uint32_t ref_and_seq_consumed;  // ZIP bamass
    uint32_t insertions;            // ZIP bamass

    PizDeepSeqFlags piz_deep_flags; // PIZ Deep
    bool is_forward;                // ZIP bamass: true if the FASTQ sequence is forward relative to the reference

    #define after_fastq_vb_zero_per_line R1_vb_i
    // --------- END OF current line ----------

    // data used for segging R2
    uint32_t R1_vb_i;               // ZIP/PIZ R2: the equivalent vb_i in the R1 (vb_i >= 1)
    STR (R1_last_qname);            // ZIP R2: pointer into z_file->R1_last_qname: last reversed, nul-terminated canonical qname of the correspondining R1 VB
    uint32_t R1_num_lines;          // ZIP R2: number of reads (FASTQ lines) in the corresponding R1 vb
    int32_t R2_lowest_read;         // ZIP: in fastq_unconsumed: lowest/highest index in txt_data of reads that we tested against R1_last_qname
    int32_t R2_highest_read; 

    uint64_t first_line;            // ZIP R2: used for optimizing QNAME  

    bool has_extra;                 // ZIP: a VB-private copy of segconf.has_extra

    #define FASTQ_NUM_TOP_LEVEL_FIELDS 16
    bool item_filter[FASTQ_NUM_TOP_LEVEL_FIELDS];  // PIZ: item filter - true if to reconstrut, false if not. index is item in toplevel container. length must be >= nitems_lo of all genozip verions

    // stats
    uint32_t deep_stats[NUM_DEEP_STATS_ZIP];  // ZIP: stats collection regarding deep / bamass
    uint32_t seq_comp_len;          // total comrpessed len of sections contributing to SEQ compression. Use for calculating bamass_trims.
    uint32_t num_consumed;          // number of consumed deep/bamass ents by this VB
    uint32_t num_monochar;          // number of lines segged monochar
    Multiplexer2 mux_ultima_c;
} VBlockFASTQ, *VBlockFASTQP;

#define VB_FASTQ ((VBlockFASTQP)vb)

#define DATA_LINE(i) B(ZipDataLineFASTQ, vb->lines, i)
 
typedef struct __attribute__((packed, aligned(4))) { // 20 bytes
    uint32_t seq_hash;
    uint32_t aln_lo;
    union { 
        uint64_t vb_qname_hash;            // before linking: full 64b hash (note: 54 MSb overlap with z_qname_hash_hi). 
        struct { 
#ifdef __BIG_ENDIAN__                      // see: https://gcc.gnu.org/legacy-ml/gcc/2004-09/msg00581.html
            uint64_t z_qname_hash_hi : 54; // after linking: high 54 bits of vb_qname_hash (first field, because its big endian)
#endif 
            uint64_t z_aln_hi        : 8;  // after linking: high byte of 40-bit aln
            uint64_t z_is_forward    : 1;  // after linking: true if the FASTQ seq is forward relative to the reference 
            uint64_t z_is_long_gpos  : 1;  // after linking: gpos is 5 bytes instead of 4
#ifdef __LITTLE_ENDIAN__
            uint64_t z_qname_hash_hi : 54; // after linking: high 54 bits of vb_qname_hash (little endian) (low bits are implied by linked list)
#endif 
        };
    };
    union {
        uint32_t next;                     // after linking: linked list of entries with the with the same (hash.qname & mask)
        struct {
            bool vb_is_forward;            // before linking
            bool vb_is_long_gpos;          // before linking
        };     
    };
} BamAssEnt, *BamAssEntP; 
#define NO_NEXT             0xffffffff
#define MAX_BAMASS_ENTS     0xfffffffe // our current implementation is limited to 4G reads because the linked list is 32b
#define MAX_BAMASS_ALNS_LEN 0xfffffffffe

extern uint64_t global_num_consumed;

// DESC
extern void fastq_seg_QNAME (VBlockFASTQP vb, STRp(qname), uint32_t line1_len, bool deep, uint32_t uncanonical_suffix_len);
extern void fastq_seg_DESC (VBlockFASTQP vb, STRp(desc), bool deep_qname2, uint32_t uncanonical_suffix_len);
extern void fastq_seg_LINE3 (VBlockFASTQP vb, STRp(qline3), STRp(qname), STRp(desc));
extern void fastq_segconf_analyze_DESC (VBlockFASTQP vb, STRp(desc));
extern bool fastq_is_line3_copy_of_line1 (STRp(qname), STRp(line3), uint32_t desc_len);

// SAUX
extern bool fastq_segconf_analyze_saux (VBlockFASTQP vb, STRp(saux));
extern void fastq_seg_saux (VBlockFASTQP vb, STRp(saux));

// AUX
extern int fastq_seg_aux (VBlockFASTQP vb, STRps(item));

// Agilent stuff
extern void agilent_seg_initialize (VBlockP vb);
extern void agilent_seg_RX (VBlockP vb, ContextP ctx, STRp(rx), unsigned add_bytes); // RX and QX are also in sam_private.h.
extern void agilent_seg_QX (VBlockP vb, ContextP ctx, STRp(qx), unsigned add_bytes);

// SEQ
extern void fastq_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), bool deep);

// QUAL
extern void fastq_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qual));

// Pairing stuff
extern bool fastq_piz_R1_test_aligned (VBlockFASTQP vb);
extern int64_t reconstruct_from_pair_int (VBlockFASTQP vb, ContextP ctx);

// Deep stuff
extern void fastq_deep_zip_initialize (void);
extern bool fastq_seg_deep (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qname), STRp(qname2), STRp(seq), STRp(qual), uint32_t *uncanonical_suffix_len);
extern void fastq_deep_seg_finalize_segconf (uint32_t n_lines);
extern void fastq_deep_seg_initialize (VBlockFASTQP vb);
extern void fastq_deep_seg_QNAME (VBlockFASTQP vb, Did did_i, STRp(qname), uint32_t uncanonical_suffix_len, unsigned add_bytes);
extern void fastq_deep_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), ContextP bitmap_ctx, ContextP nonref_ctx);
extern void fastq_deep_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, ContextP qual_ctx, uint32_t qual_len);
extern void fastq_deep_zip_after_compute (VBlockFASTQP vb);
extern void fastq_deep_zip_finalize (void);
extern void fastq_deep_piz_wait_for_deep_data (void);
extern bool fastq_deep_seg_find_subseq (VBlockFASTQP vb, STRp (fastq_seq), uint32_t sam_seq_len, uint32_t sam_seq_hash, bool allow_offset, uint32_t *sam_seq_offset);

// BAM-Assist stuff
extern void fastq_bamass_populate (void);
extern void fastq_bamass_seg_initialize (VBlockFASTQP vb);
extern void fastq_bamass_seg_finalize (VBlockFASTQP vb);
extern void fastq_bamass_zip_comp_cb (VBlockFASTQP vb, ContextP ctx, SectionType st, uint32_t comp_len);
extern void fastq_bamass_zip_after_compress (VBlockFASTQP vb);
extern void fastq_bamass_zip_finalize (bool is_last_fastq);
extern void fastq_bamass_retrieve_ent (VBlockP vb, const BamAssEnt *e, bool get_cigar, bool *is_fwd, PosType64 *gpos, uint32_t *seq_len, uint32_t *ref_consumed, uint32_t *ref_and_seq_consumed, uint32_t *insertions);
extern MappingType fastq_bamass_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), bool is_pair_2, PosType64 pair_gpos, bool pair_is_forward);
extern void fastq_bamass_seg_CIGAR (VBlockFASTQP vb);
extern DeepStatsZip fastq_seg_find_bamass (VBlockFASTQP vb, ZipDataLineFASTQ *dl, DeepHash *deep_hash, STRp(seq), BamAssEnt **matching_ent);
extern StrTextLong bamass_dis_ent (VBlockP vb, const BamAssEnt *e, uint64_t qname_hash);
extern void fastq_bamass_consider_stopping_aligner (VBlockFASTQP vb);

extern Buffer bamass_ents, bamass_heads;
eSTRl(copy_qname_snip);
