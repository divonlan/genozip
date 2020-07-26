// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ZIP_INCLUDED
#define ZIP_INCLUDED

#include "genozip.h"

extern void zip_dispatcher (const char *vcf_basename, bool is_last_file);
extern void zip_generate_b250_section (VBlockP vb, ContextP ctx);

typedef enum { PD_VBLOCK_DATA, PD_DICT_DATA, PD_REFERENCE_DATA } ProcessedDataType;
extern void zip_output_processed_vb (VBlockP vb, BufferP section_list_buf, bool update_txt_file, ProcessedDataType pd_type);

extern void zip_initialize_binary_dump (const char *field, DictId *dict_id, const char *filename_ext);

// --------------------------------------------------
// utilities for use by zip_*_compress_one_vb
// --------------------------------------------------

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

#endif
