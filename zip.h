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

// --------------------------------------------------
// utilities for use by zip_*_compress_one_vb
// --------------------------------------------------

#define COMPRESS_DATA_SECTION(sec,vb_buf_name,data_type,codec,is_optional) { \
    if (vb->vb_buf_name.len) { \
        vb->vb_buf_name.len *= sizeof (data_type); \
        zfile_compress_section_data_codec ((VBlockP)vb, (sec),  &vb->vb_buf_name, NULL, 0, codec); \
        vb->vb_buf_name.len /= sizeof (data_type); /* restore */ \
    } \
}

#define COMPRESS_DATA_SECTION_CALLBACK(sec,callback,total_len,codec,is_optional) { \
    if (!is_optional || total_len) \
        zfile_compress_section_data_codec (vb_, (sec), NULL, callback, total_len, codec);\
}

#endif
