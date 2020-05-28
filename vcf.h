// ------------------------------------------------------------------
//   vcf.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VCF_INCLUDED
#define VCF_INCLUDED

#include "genozip.h"
#include "md5.h"

// SEG stuff
extern const char *vcf_seg_txt_line (VBlockP vb_, const char *field_start_line, bool *has_special_eol);
extern void vcf_seg_initialize (VBlockP vb_);

// ZIP stuff
extern void vcf_zip_set_global_samples_per_block (const char *num_samples_str);
extern void vcf_zip_initialize (void);
extern void vcf_zip_compress_one_vb  (VBlockP vb_);
extern void vcf_zip_generate_ht_gt_compress_vb_header (VBlockP vb_);

// PIZ stuff
extern bool vcf_piz_read_one_vb (VBlockP vb, SectionListEntryP sl);
extern bool vcf_v1_piz_read_one_vb();// v1 compatibility
extern void vcf_piz_uncompress_vb(); // no parameter - implicit casting of VBlockP to VBlockVCFP
extern bool vcf_piz_is_skip_section (VBlockP vb, SectionType st, DictIdType dict_id);

// ZFILE stuff
extern void vcf_zfile_compress_vb_header (VBlockP vb);
extern void vcf_zfile_update_compressed_vb_header (VBlockP vb, uint32_t vcf_first_line_i);

// VCF Header stuff
extern void vcf_header_initialize (void);
extern bool vcf_header_set_globals (const char *filename, BufferP vcf_header);
extern void vcf_header_trim_header_line (BufferP vcf_header_buf);
extern void vcf_header_keep_only_last_line (BufferP vcf_header_buf);

// VBlock stuff
extern void vcf_vb_release_vb();
extern void vcf_vb_destroy_vb();
extern void vcf_vb_cleanup_memory();
extern unsigned vcf_vb_size (void);
extern unsigned vcf_vb_zip_dl_size (void);
extern bool vcf_vb_has_haplotype_data (VBlockP vb);

// Samples stuff
extern void vcf_samples_add  (const char *samples_str);
extern bool vcf_is_sb_included (void *vb_, uint32_t sb_i);

// GL-optimize stuff
extern const char *gl_optimize_dictionary (VBlockP vb, BufferP dict, MtfNodeP nodes, uint64_t dict_start_char, unsigned num_words);
extern void gl_deoptimize_dictionary (char *data, int len);
extern bool gl_optimize (const char *snip, unsigned len, char *updated_snip, unsigned *updated_len);

// v1 compatibility
extern bool vcf_v1_header_genozip_to_vcf (Md5Hash *digest);
extern bool vcf_v1_header_get_vcf_header (uint64_t *uncompressed_data_size, uint32_t *num_samples,uint64_t *num_items_concat,
                                          Md5Hash  *md5_hash_concat, char *created, unsigned created_len);

// global parameters - set before any thread is created
extern uint32_t global_vcf_num_samples, global_vcf_num_displayed_samples;

#define VCF_DICT_ID_ALIASES \
    /*         alias                           maps to this ctx          */  \
    { DT_VCF,  &dict_id_INFO_END,              &dict_id_fields[VCF_POS]    }, \

#define VCF_LOCAL_COMPRESSOR_CALLBACKS

#define dict_id_is_vcf_info_sf   dict_id_is_type_1
#define dict_id_is_vcf_format_sf dict_id_is_type_2

#define dict_id_vcf_info_sf      dict_id_type_1
#define dict_id_vcf_format_sf    dict_id_type_2

#endif
