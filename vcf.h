// ------------------------------------------------------------------
//   vcf.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
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

// ZIP stuff
extern void vcf_zip_initialize (void);

// SEG stuff
extern const char *vcf_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void vcf_seg_initialize (VBlockP vb_);
extern void vcf_zip_after_compute (VBlockP vb);
extern void vcf_seg_finalize (VBlockP vb_);
extern bool vcf_seg_is_small (ConstVBlockP vb, DictId dict_id);

// PIZ stuff
extern bool vcf_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
CONTAINER_FILTER_FUNC (vcf_piz_filter);
CONTAINER_CALLBACK (vcf_piz_container_cb);

// VCF Header stuff
extern void vcf_header_piz_init (void);
extern bool vcf_inspect_txt_header (BufferP txt_header);
extern uint32_t vcf_header_get_num_samples (void);

// VBlock stuff
extern void vcf_vb_release_vb();
extern void vcf_vb_destroy_vb();
extern void vcf_vb_cleanup_memory();
extern unsigned vcf_vb_size (void);
extern unsigned vcf_vb_zip_dl_size (void);
extern bool vcf_vb_has_haplotype_data (VBlockP vb);

// Liftover - INFO fields
#define INFO_LIFTOVER  "LIFTOVER"
#define INFO_LIFTOVER_LEN 8
#define INFO_LIFTBACK  "LIFTBACK"
#define INFO_LIFTBACK_LEN 8
#define INFO_LIFTREJD  "LIFTREJD"
#define INFO_LIFTREJD_LEN 8

// Liftover - header keys
#define HK_GENOZIP_CMD "##genozip_command="
#define HK_LO_CONTIG   "##liftover_contig="     
#define HK_LB_CONTIG   "##liftback_contig="
#define HK_LO_REF      "##liftover_reference="
#define HK_LB_REF      "##liftback_reference="
#define HK_LB_REJECT   "##liftback_reject="
#define HK_CHAIN       "##chain="
#define HK_DC          "##dual_coordinates"
#define HK_DC_PRIMARY  HK_DC"=PRIMARY"
#define HK_DC_LUFT     HK_DC"=LUFT"
#define KH_INFO_LO     "##INFO=<ID=" INFO_LIFTOVER ",Number=5,Type=String,Description=\"dual-coordinates VCF: Information for lifting over the variant to luft coordinates\",Source=\"genozip\",Version=\"%s\">"
#define KH_INFO_LB     "##INFO=<ID=" INFO_LIFTBACK ",Number=5,Type=String,Description=\"dual-coordinates VCF: Information for retrieving the variant in the primary coordinates\",Source=\"genozip\",Version=\"%s\">"
#define KH_INFO_LR     "##INFO=<ID=" INFO_LIFTREJD ",Number=1,Type=String,Description=\"dual-coordinates VCF: Reason variant was rejected for lift over\",Source=\"genozip\",Version=\"%s\">"

extern void vcf_seg_ref_alt (VBlockP vb, const char *ref, unsigned ref_len, const char *alt, unsigned alt_len);

// Samples stuff
extern void vcf_samples_add  (const char *samples_str);

#define VCF_SPECIAL { vcf_piz_special_REFALT, vcf_piz_special_FORMAT, vcf_piz_special_AC, vcf_piz_special_SVLEN, \
                      vcf_piz_special_DS, vcf_piz_special_BaseCounts, vcf_piz_special_SF, \
                      vcf_piz_special_OREF, vcf_piz_special_LIFTREJD, vcf_piz_special_LIFTBACK }
SPECIAL (VCF, 0, REFALT,     vcf_piz_special_REFALT);
SPECIAL (VCF, 1, FORMAT,     vcf_piz_special_FORMAT)
SPECIAL (VCF, 2, AC,         vcf_piz_special_AC);
SPECIAL (VCF, 3, SVLEN,      vcf_piz_special_SVLEN);
SPECIAL (VCF, 4, DS,         vcf_piz_special_DS);
SPECIAL (VCF, 5, BaseCounts, vcf_piz_special_BaseCounts);
SPECIAL (VCF, 6, SF,         vcf_piz_special_SF);
SPECIAL (VCF, 7, OREF,       vcf_piz_special_OREF);     // added v12
SPECIAL (VCF, 8, LIFTREJD,   vcf_piz_special_LIFTREJD); // added v12 - maybe be used by other data types too in the future
SPECIAL (VCF, 9, LIFTBACK,   vcf_piz_special_LIFTBACK); // added v12 - maybe be used by other data types too in the future

#define NUM_VCF_SPECIAL 10

// Translators for Luft (=secondary coordinates)
TRANSLATOR (VCF, VCF,   1,  CHROM,  vcf_piz_luft_CHROM)
TRANSLATOR (VCF, VCF,   2,  POS,    vcf_piz_luft_POS)
TRANSLATOR (VCF, VCF,   3,  REFALT, vcf_piz_luft_REFALT)
TRANSLATOR (VCF, VCF,   4,  AC,     vcf_piz_luft_AC)
TRANSLATOR (VCF, VCF,   5,  AF,     vcf_piz_luft_AF)
TRANSLATOR (VCF, VCF,   6,  AD,     vcf_piz_luft_AD)
TRANSLATOR (VCF, VCF,   7,  END,    vcf_piz_luft_END)
TRANSLATOR (VCF, VCF,   8,  GT,     vcf_piz_luft_GT)
TRANSLATOR (VCF, VCF,   9,  GL,     vcf_piz_luft_GL)

#define NUM_VCF_TRANS   10 // including "none"
#define VCF_TRANSLATORS { NULL /* none */, vcf_piz_luft_CHROM, vcf_piz_luft_POS, vcf_piz_luft_REFALT, \
                          vcf_piz_luft_AC, vcf_piz_luft_AF, \
                          vcf_piz_luft_AD, vcf_piz_luft_END, vcf_piz_luft_GT, vcf_piz_luft_GL }

#define VCF_DICT_ID_ALIASES \
    /*         alias                           maps to this ctx          */  \
    { DT_VCF,  &dict_id_INFO_END,              &dict_id_fields[VCF_POS]    }, \

#define VCF_LOCAL_GET_LINE_CALLBACKS

#define dict_id_is_vcf_info_sf   dict_id_is_type_1
#define dict_id_is_vcf_format_sf dict_id_is_type_2

#define DTYPE_VCF_INFO   DTYPE_1
#define DTYPE_VCF_FORMAT DTYPE_2

#endif
