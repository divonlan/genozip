// ------------------------------------------------------------------
//   vcf.h
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VCF_INCLUDED
#define VCF_INCLUDED

#include "genozip.h"
#include "digest.h"
#include "sections.h"

// default max number of samples in each sample block within a variant block. user configurable with --sblock
#define VCF_SAMPLES_PER_VBLOCK "4096" 

#define VCF_MAX_PLOIDY 100  // set to a reasonable 100 to avoid memory allocation explosion in case of an error in the VCF file
#if VCF_MAX_PLOIDY > 65535
#error "VCF_MAX_PLOIDY cannot go beyond 65535 VBlockVCF.ploidy are uint16_t"
#endif

typedef uint8_t Allele; // elements of ht_matrix: values 48->147 for allele 0 to 99, '*' for unused, '%', '-'

// ZIP stuff
extern void vcf_zip_initialize (void);

// SEG stuff
extern const char *vcf_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void vcf_seg_initialize (VBlockP vb_);
extern void vcf_zip_after_compute (VBlockP vb);
extern void vcf_seg_finalize (VBlockP vb_);
extern bool vcf_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern TranslatorId vcf_lo_luft_trans_id (DictId dict_id, char number);

// PIZ stuff
extern bool vcf_piz_read_one_vb (VBlockP vb, Section sl);
extern bool vcf_vb_is_luft (VBlockP vb);
extern bool vcf_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
CONTAINER_FILTER_FUNC (vcf_piz_filter);
CONTAINER_CALLBACK (vcf_piz_container_cb);

// VCF Header stuff
extern void vcf_header_piz_init (void);
extern bool vcf_inspect_txt_header (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);
extern uint32_t vcf_header_get_num_samples (void);

// VBlock stuff
extern void vcf_vb_release_vb();
extern void vcf_vb_destroy_vb();
extern void vcf_vb_cleanup_memory();
extern unsigned vcf_vb_size (void);
extern unsigned vcf_vb_zip_dl_size (void);
extern bool vcf_vb_has_haplotype_data (VBlockP vb);

// Liftover - INFO fields
#define INFO_LUFT  "LUFT"
#define INFO_PRIM  "PRIM"
#define INFO_LREJ  "Lrej"
#define INFO_PREJ  "Prej" // lower case so Prej doesn't have the same first 2 chars as PRIM (to not conflict in dict_id_to_did_i_map)
#define INFO_DVCF_LEN (sizeof INFO_LUFT - 1) // these 4 must be the same length

#define VCF_CONTIG_FMT "##contig=<ID=%.*s,length=%"PRId64">"

// Samples stuff
extern void vcf_samples_add  (const char *samples_str);

#define VCF_SPECIAL { vcf_piz_special_main_REFALT, vcf_piz_special_FORMAT, vcf_piz_special_INFO_AC, vcf_piz_special_INFO_SVLEN, \
                      vcf_piz_special_FORMAT_DS, vcf_piz_special_INFO_BaseCounts, vcf_piz_special_INFO_SF, vcf_piz_special_MINUS,  \
                      vcf_piz_special_LIFT_REF, vcf_piz_special_COPYSTAT, vcf_piz_special_other_REFALT, vcf_piz_special_COPYPOS }
SPECIAL (VCF, 0,  main_REFALT,  vcf_piz_special_main_REFALT);
SPECIAL (VCF, 1,  FORMAT,       vcf_piz_special_FORMAT)
SPECIAL (VCF, 2,  AC,           vcf_piz_special_INFO_AC);
SPECIAL (VCF, 3,  SVLEN,        vcf_piz_special_INFO_SVLEN);
SPECIAL (VCF, 4,  DS,           vcf_piz_special_FORMAT_DS);
SPECIAL (VCF, 5,  BaseCounts,   vcf_piz_special_INFO_BaseCounts);
SPECIAL (VCF, 6,  SF,           vcf_piz_special_INFO_SF);
SPECIAL (VCF, 7,  MINUS,        vcf_piz_special_MINUS);        // added v12
SPECIAL (VCF, 8,  LIFT_REF,     vcf_piz_special_LIFT_REF);     // added v12
SPECIAL (VCF, 9,  COPYSTAT,     vcf_piz_special_COPYSTAT);     // added v12
SPECIAL (VCF, 10, other_REFALT, vcf_piz_special_other_REFALT); // added v12
SPECIAL (VCF, 11, COPYPOS,      vcf_piz_special_COPYPOS);      // added v12
#define NUM_VCF_SPECIAL 12

// Translators for Luft (=secondary coordinates)
TRANSLATOR (VCF, VCF,   1, G,      vcf_piz_luft_G)       // same order as LiftOverStatus starting LO_CANT_G
TRANSLATOR (VCF, VCF,   2, R,      vcf_piz_luft_R)
TRANSLATOR (VCF, VCF,   3, R2,     vcf_piz_luft_R2)
TRANSLATOR (VCF, VCF,   4, A_AN,   vcf_piz_luft_A_AN)
TRANSLATOR (VCF, VCF,   5, A_1,    vcf_piz_luft_A_1)
TRANSLATOR (VCF, VCF,   6, GT,     vcf_piz_luft_GT)      
TRANSLATOR (VCF, VCF,   7, END,    vcf_piz_luft_END)      
TRANSLATOR (VCF, VCF,   8, XREV,   vcf_piz_luft_XREV)      

#define NUM_VCF_TRANS   9 // including "none"
#define VCF_TRANSLATORS { NULL /* none */, vcf_piz_luft_G, vcf_piz_luft_R, vcf_piz_luft_R2, vcf_piz_luft_A_AN, \
                          vcf_piz_luft_A_1, vcf_piz_luft_GT, vcf_piz_luft_END, vcf_piz_luft_XREV }

typedef struct {
    const char *alg_name;
    enum { TW_NEVER, TW_ALWAYS, TW_REF_ALT_SWITCH, TW_XSTRAND } upon;
} LuftTransLateProp;

// names of INFO / FORMAT algorithms, goes into VCF header's ##INFO / ##FORMAT "RenderAlg" attribute
                           /* Algorithm   When-triggered */
#define DVCF_TRANS_PROPS { { "NONE",      TW_NEVER          },   /* never translate */\
                           { "G",         TW_REF_ALT_SWITCH },   /* reshuffle a  'G' vector (one element per genotype) if REF<>ALT changed */\
                           { "R",         TW_REF_ALT_SWITCH },   /* reshuffle an 'R' vector (one element per ref/alt allele) if REF<>ALT changed */\
                           { "R2",        TW_REF_ALT_SWITCH },   /* reshuffle a vector with 2 elements per ref/alt allele, if REF<>ALT changed */\
                           { "A_AN",      TW_REF_ALT_SWITCH },   /* recalculate an 'A' vector (one element per ALT allele) if REF<>ALT changed, who's elements, including a missing element for REF, add up to AN (example: AC). */ \
                           { "A_1",       TW_REF_ALT_SWITCH },   /* recalculate an 'A' vector (one element per ALT allele) if REF<>ALT changed, who's elements, including a missing element for REF, add up to 1 (example: AF). */ \
                           { "GT",        TW_REF_ALT_SWITCH },   /* recalculate the allele numbers FORMAT/GT if REF<>ALT changed */ \
                           { "END",       TW_ALWAYS         },   /* recalculate INFO/END */\
                           { "XREV",      TW_XSTRAND        } }  /* reverse the elements of a vector if XSTRAND. Example: INFO/BaseCounts */ 
                           
extern const LuftTransLateProp ltrans_props[NUM_VCF_TRANS];

#define needs_translation(ctx)  (z_dual_coords && (ctx)->luft_trans && \
    ((ltrans_props[(ctx)->luft_trans].upon == TW_REF_ALT_SWITCH && last_ostatus == LO_OK_REF_ALT_SWTCH) || \
     (ltrans_props[(ctx)->luft_trans].upon == TW_ALWAYS         && LO_IS_OK (last_ostatus))             || \
     (ltrans_props[(ctx)->luft_trans].upon == TW_XSTRAND        && LO_IS_OK (last_ostatus) && *vb->contexts[VCF_oXSTRAND].last_snip == 'X')))

#define VCF_DICT_ID_ALIASES \
    /*         alias                           maps to this ctx          */  \
    { DT_VCF,  &dict_id_INFO_END,              &dict_id_fields[VCF_POS]    }, \

#define VCF_LOCAL_GET_LINE_CALLBACKS

#define dict_id_is_vcf_info_sf   dict_id_is_type_1
#define dict_id_is_vcf_format_sf dict_id_is_type_2

#define DTYPE_VCF_INFO   DTYPE_1
#define DTYPE_VCF_FORMAT DTYPE_2

enum { oCHROM, oPOS, oREF, oXSTRAND, oSTATUS }; // order of fields as defined in data_types.h
#define ODID(offset) (DTFZ(ochrom)+(DidIType)(offset))

#endif
