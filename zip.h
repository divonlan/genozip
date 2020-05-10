// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ZIP_INCLUDED
#define ZIP_INCLUDED

#include "genozip.h"

extern void zip_dispatcher (const char *vcf_basename, unsigned max_threads, bool is_last_file);

extern void zip_vcf_set_global_samples_per_block (const char *num_samples_str);

extern void zip_generate_b250_section (VBlockP vb, MtfContextP ctx);
extern void zip_generate_and_compress_fields (VBlockP vb);
extern void zip_generate_and_compress_subfields (VBlockP vb, ConstSubfieldMapperP mapper);
extern void zip_generate_and_compress_subfields2 (VBlockP vb, SectionType sec_dict);

extern void zip_vcf_initialize (void);

extern void zip_vcf_compress_one_vb  (VBlockP vb_);
extern void zip_sam_compress_one_vb  (VBlockP vb_);
extern void zip_fast_compress_one_vb (VBlockP vb_);
extern void zip_gff3_compress_one_vb (VBlockP vb_);
extern void zip_me23_compress_one_vb (VBlockP vb_);

extern void zip_vcf_generate_ht_gt_compress_vb_header (VBlockP vb_);
extern uint32_t global_vcf_samples_per_block;

// --------------------------------------------------
// utilities for use by zip_*_compress_one_vb
// --------------------------------------------------

#define COMPRESS_START(vblock_type) \
    START_TIMER; \
    vblock_type *vb = (vblock_type *)vb_; \
    /* if we're vb_i=1 lock, and unlock only when we're done merging. all other vbs need  \
       to wait for our merge. that is because our dictionaries are sorted */ \
    if (vb->vblock_i == 1) mtf_vb_1_lock(vb_); \
    /* allocate memory for the final compressed data of this vb. allocate 20% of the  \
       vb size on the original file - this is normally enough. if not, we will realloc downstream */ \
    buf_alloc (vb, &vb->z_data, vb->vb_data_size / 5, 1.2, "z_data", 0); \
    /* clone global dictionaries while granted exclusive access to the global dictionaries */ \
    mtf_clone_ctx (vb_);

#define COMPRESS_DATA_SECTION(sec,vb_buf_name,data_type,comp_alg,is_optional) { \
    if (vb->vb_buf_name.len) { \
        vb->vb_buf_name.len *= sizeof (data_type); \
        zfile_compress_section_data_alg ((VBlockP)vb, (sec),  &vb->vb_buf_name, NULL, 0, comp_alg); \
        vb->vb_buf_name.len /= sizeof (data_type); /* restore */ \
    } \
}

#define COMPRESS_DATA_SECTION_CALLBACK(sec,callback,total_len,comp_alg,is_optional) { \
    if (!is_optional || total_len) \
        zfile_compress_section_data_alg (vb_, (sec), NULL, callback, total_len, comp_alg);\
}

#define COMPRESS_DONE \
    /* tell dispatcher this thread is done and can be joined.  \
       thread safety: this isn't protected by a mutex as it will just be false until it at some point turns to true \
       this this operation needn't be atomic, but it likely is anyway */ \
    vb->is_processed = true; \
    COPY_TIMER (vb->profile.compute);

#endif
