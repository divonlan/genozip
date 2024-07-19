// ------------------------------------------------------------------
//   fastq.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "sections.h"
#include "multiplexer.h"

// SAM and FASTQ share the same Dids and DictIds
#define FASTQ_CONTIG        SAM_RNAME
#define FASTQ_QNAME         SAM_QNAME
#define FASTQ_Q0NAME        SAM_Q0NAME
#define FASTQ_QmNAME        SAM_QmNAME
#define FASTQ_QNAME2        SAM_QNAME2
#define FASTQ_QmNAME2       SAM_QmNAME2
#define FASTQ_AUX           SAM_AUX
#define FASTQ_SQBITMAP      SAM_SQBITMAP
#define FASTQ_NONREF        SAM_NONREF
#define FASTQ_NONREF_X      SAM_NONREF_X
#define FASTQ_GPOS          SAM_GPOS
#define FASTQ_GPOS_DELTA    SAM_GPOS_DELTA
#define FASTQ_GPOS_R2       SAM_GPOS_R2    // 15.0.58: used for interleaved files
#define FASTQ_STRAND        SAM_STRAND
#define FASTQ_STRAND_R2     SAM_STRAND_R2  // 15.0.58: used for interleaved files
#define FASTQ_SEQMIS_A      SAM_SEQMIS_A
#define FASTQ_SEQMIS_C      SAM_SEQMIS_C
#define FASTQ_SEQMIS_G      SAM_SEQMIS_G
#define FASTQ_SEQMIS_T      SAM_SEQMIS_T
#define FASTQ_QUAL          SAM_QUAL
#define FASTQ_DOMQRUNS      SAM_DOMQRUNS
#define FASTQ_QUALMPLX      SAM_QUALMPLX
#define FASTQ_DIVRQUAL      SAM_DIVRQUAL
#define FASTQ_TOPLEVEL      SAM_TOPLEVEL
#define FASTQ_BUDDY         SAM_BUDDY
#define FASTQ_TAXID         SAM_TAXID
#define FASTQ_DEBUG_LINES   SAM_DEBUG_LINES

#define _FASTQ_CONTIG       _SAM_RNAME
#define _FASTQ_QNAME        _SAM_QNAME
#define _FASTQ_Q0NAME       _SAM_Q0NAME
#define _FASTQ_Q1NAME       _SAM_Q1NAME
#define _FASTQ_QmNAME       _SAM_QmNAME
#define _FASTQ_QNAME2       _SAM_QNAME2
#define _FASTQ_Q0NAME2      _SAM_Q0NAME2
#define _FASTQ_Q1NAME2      _SAM_Q1NAME2
#define _FASTQ_QmNAME2      _SAM_QmNAME2
#define _FASTQ_AUX          _SAM_AUX
#define _FASTQ_SQBITMAP     _SAM_SQBITMAP
#define _FASTQ_NONREF       _SAM_NONREF
#define _FASTQ_NONREF_X     _SAM_NONREF_X
#define _FASTQ_GPOS         _SAM_GPOS
#define _FASTQ_GPOS_DELTA   _SAM_GPOS_DELTA
#define _FASTQ_GPOS_R2      _SAM_GPOS_R2    
#define _FASTQ_STRAND       _SAM_STRAND
#define _FASTQ_STRAND_R2    _SAM_STRAND_R2  
#define _FASTQ_SEQMIS_A     _SAM_SEQMIS_A
#define _FASTQ_SEQMIS_C     _SAM_SEQMIS_C
#define _FASTQ_SEQMIS_G     _SAM_SEQMIS_G
#define _FASTQ_SEQMIS_T     _SAM_SEQMIS_T
#define _FASTQ_QUAL         _SAM_QUAL
#define _FASTQ_DOMQRUNS     _SAM_DOMQRUNS
#define _FASTQ_QUALMPLX     _SAM_QUALMPLX
#define _FASTQ_DIVRQUAL     _SAM_DIVRQUAL
#define _FASTQ_TOPLEVEL     _SAM_TOPLEVEL
#define _FASTQ_BUDDY        _SAM_BUDDY
#define _FASTQ_TAXID        _SAM_TAXID
#define _FASTQ_DEBUG_LINES  _SAM_DEBUG_LINES

#define FASTQ_PREDEFINED    SAM_PREDEFINED

#define NUM_FASTQ_FIELDS    NUM_SAM_FIELDS

#define NO_PAIR_FMT_PREFIX "--pair cannot be used because %s is not perfectly paired with its counterpart (read names differ or are not aligned) (technical: "
#define NO_PAIR_FMT_SUFFIX (flag.deep ? " Solution: add --not-paired" : "")

// Txtfile stuff
extern int32_t fastq_unconsumed (VBlockP vb, uint32_t first_i);
extern bool fastq_txtfile_sync_to_R1_by_num_lines (VBlockP vb, uint32_t bytes_requested, uint32_t len, bool no_read_expected, uint32_t *my_vb_size);
extern bool is_fastq (STRp(header), bool *need_more);
extern bool is_fastq_pair_2 (VBlockP vb);
extern VBIType fastq_get_R1_vb_i (VBlockP vb);
extern uint32_t fastq_get_R1_num_lines (VBlockP vb);
extern uint32_t fastq_get_R1_txt_data_len (VBlockP vb);
extern rom fastq_get_R1_last_qname (VBlockP vb);
extern void fastq_zip_set_txt_header_flags (struct FlagsTxtHeader *f);

// ZIP Stuff
extern void fastq_zip_initialize (void);
extern rom fastq_zip_modify (VBlockP vb, rom line_start, uint32_t remaining);
extern void fastq_segconf_set_r1_or_r2 (void);
extern void fastq_zip_after_segconf (void);
extern void fastq_zip_finalize (bool is_last_user_txt_file);
extern void fastq_zip_init_vb (VBlockP vb);
extern void fastq_zip_after_compute (VBlockP vb);
extern bool fastq_zip_use_pair_assisted (DictId dict_id, SectionType st);
extern bool fastq_zip_use_pair_identical (DictId dict_id);
extern uint32_t fastq_zip_get_seq_len (VBlockP vb, uint32_t line_i) ;

COMPRESSOR_CALLBACK (fastq_zip_seq);
COMPRESSOR_CALLBACK(fastq_zip_qual); // used by codec_longr_compress

// SEG Stuff
extern void fastq_seg_initialize (VBlockP vb);
extern void fastq_segconf_finalize (VBlockP vb);
extern void fastq_seg_finalize (VBlockP vb);
extern bool fastq_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern rom fastq_seg_txt_line (VBlockP vb, rom line_start, uint32_t remaining, bool *has_13);
extern rom fastq_assseg_line (VBlockP vb);
extern void fastq_seg_r2_gpos (VBlockP vb, PosType64 r1_pos, PosType64 r2_gpos);
extern void fastq_seg_interleaved_gpos (VBlockP vb, PosType64 pair_gpos/*only if we are R2*/, PosType64 gpos);
extern void fastq_update_qual_len (VBlockP vb, uint32_t line_i, uint32_t new_len);
extern Multiplexer2P fastq_get_ultima_c_mux (VBlockP vb);

// PIZ Stuff
extern bool fastq_piz_initialize (CompIType comp_i);
extern void fastq_piz_header_init (void);
extern void fastq_piz_process_recon (VBlockP vb);
extern void fastq_piz_before_read (VBlockP vb);
extern bool fastq_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header);
CONTAINER_FILTER_FUNC (fastq_piz_filter);
extern IS_SKIP (fastq_piz_is_skip_section);
extern void fastq_recon_aligned_SEQ (VBlockP vb, STRp(seq_len_str), ReconType reconstruct);
extern CONTAINER_CALLBACK (fastq_piz_container_cb);

// VBlock stuff
extern unsigned fastq_vb_size (DataType dt);
extern unsigned fastq_vb_zip_dl_size (void);
extern void fastq_reset_line (VBlockP vb);

// file pairing (--pair) stuff
extern void fastq_read_R1_data (VBlockP vb, VBIType R1_vb_i);

// FASTQ-specific fields in genozip header
extern void fastq_zip_genozip_header (SectionHeaderGenozipHeaderP header);
extern void fastq_piz_genozip_header (ConstSectionHeaderGenozipHeaderP header);

#define FASTQ_LOCAL_GET_LINE_CALLBACKS           \
    { DT_FASTQ, _FASTQ_QUAL,   fastq_zip_qual }, 
 // { DT_FASTQ, _FASTQ_NONREF, fastq_zip_seq  },  this callback is called directly from codec_longr, it is not added to the list, since it should NOT be called in other codecs

#define FASTQ_DICT_ID_ALIASES                               \
    /*           type        alias        maps to      */   \
    { DT_FASTQ,  ALIAS_DICT, _FASTQ_E1L, _FASTQ_E2L     },  /* each of E1L and E2L *may* have different b250s if lines differ in their EOL, but normally they are both all-the-same, so only the dict of E2L survives, and no b250 sections */

typedef enum { FQ_COMP_R1, FQ_COMP_R2 } FastqComponentType;
#define FASTQ_COMP_NAMES { "FQR1", "FQR2" }
 
// SPECIALs
SPECIAL (FASTQ, 0,  unaligned_SEQ,     fastq_special_unaligned_SEQ);          // v14
SPECIAL (FASTQ, 1,  PAIR2_GPOS,        fastq_special_PAIR2_GPOS);             // v14
SPECIAL (FASTQ, 2,  mate_lookup,       fastq_special_mate_lookup);            // v14
SPECIAL (FASTQ, 3,  set_deep,          fastq_special_set_deep);               // v15
SPECIAL (FASTQ, 4,  deep_copy_QNAME,   fastq_special_deep_copy_QNAME);        // v15
SPECIAL (FASTQ, 5,  deep_copy_SEQ,     fastq_special_deep_copy_SEQ);          // v15
SPECIAL (FASTQ, 6,  deep_copy_QUAL,    fastq_special_deep_copy_QUAL);         // v15
SPECIAL (FASTQ, 7,  backspace,         fastq_special_backspace);              // v15
SPECIAL (FASTQ, 8,  copy_line1,        fastq_special_copy_line1);             // v15
SPECIAL (FASTQ, 9,  monochar_QUAL,     fastq_special_monochar_QUAL);          // v15
SPECIAL (FASTQ, 10, ULTIMA_C,          ultima_c_piz_special_DEMUX_BY_Q4NAME); // introduced 15.0.15
SPECIAL (FASTQ, 11, AGENT_RX,          agilent_special_AGENT_RX);             // introduced 15.0.23
SPECIAL (FASTQ, 12, AGENT_QX,          agilent_special_AGENT_QX);             // introduced 15.0.23
SPECIAL (FASTQ, 13, qname_rng2seq_len, special_qname_rng2seq_len);            // introduced 15.0.26
SPECIAL (FASTQ, 14, DEMUX_BY_R,        fastq_special_DEMUX_BY_R);             // introduced 15.0.58
