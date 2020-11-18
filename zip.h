// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ZIP_INCLUDED
#define ZIP_INCLUDED

#include "genozip.h"

extern void zip_one_file (const char *vcf_basename, bool is_last_file);

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
