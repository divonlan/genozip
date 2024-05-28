// ------------------------------------------------------------------
//   zfile.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "sections.h"
#include "digest.h"

// --------
// ZIP side
// --------

extern void zfile_compress_genozip_header (void);

extern void zfile_compress_section_data_ex (VBlockP vb, ContextP ctx, SectionType section_type, 
                                            BufferP section_data, LocalGetLineCB callback, uint32_t total_len, 
                                            Codec codec, SectionFlags flags,
                                            rom context_name); // NULL if not context data

#define zfile_compress_section_data(vb, section_type, section_data) \
    zfile_compress_section_data_ex ((vb), NULL, (section_type), (section_data), NULL, 0, CODEC_BZ2, SECTION_FLAGS_NONE, NULL)

extern uint32_t zfile_compress_b250_data  (VBlockP vb, ContextP ctx);
extern uint32_t zfile_compress_local_data (VBlockP vb, ContextP ctx, uint32_t sample_size);

extern void zfile_compress_vb_header (VBlockP vb);
extern void zfile_update_compressed_vb_header (VBlockP vb);

extern void zfile_update_txt_header_section_header (uint64_t pos_of_current_vcf_header);

extern void zfile_remove_ctx_group_from_z_data (VBlockP vb, Did did_i);

extern void zfile_output_processed_vb (VBlockP vb);
extern void zfile_output_processed_vb_ext (VBlockP vb, bool background);

// --------
// PIZ side
// --------

extern uint64_t zfile_read_genozip_header_get_offset (bool as_is);
extern bool zfile_read_genozip_header (SectionHeaderGenozipHeaderP header, FailType fail_type);

extern SectionHeaderUnion zfile_read_section_header_do (VBlockP vb, Section sec, SectionType expected_sec_type, FUNCLINE);
#define zfile_read_section_header(vb, sec, expected_sec_type) \
    zfile_read_section_header_do ((VBlockP)(vb), (sec), (expected_sec_type), __FUNCLINE)

#define SECTION_SKIPPED ((int32_t)-1)
extern int32_t zfile_read_section_do (FileP file, VBlockP vb, uint32_t original_vb_i, 
                                      BufferP data /* buffer to append */, rom buf_name,
                                      SectionType expected_sec_type, 
                                      Section sec, FUNCLINE); 
#define zfile_read_section(file,vb,original_vb_i,data,buf_name,expected_sec_type,sec) \
    zfile_read_section_do ((file),(VBlockP)(vb),(original_vb_i),(data),(buf_name),(expected_sec_type),(sec), __FUNCLINE)

extern void zfile_uncompress_section (VBlockP vb, SectionHeaderUnionP section_header, BufferP uncompressed_data, rom uncompressed_data_buf_name, uint32_t expected_vb_i, SectionType expected_section_type);
extern void zfile_uncompress_section_into_buf (VBlockP vb, SectionHeaderUnionP section_header_p, uint32_t expected_vb_i, SectionType expected_section_type, BufferP dst_buf, char *dst);

#define zfile_get_global_section(HeaderType, sec, out_buf, out_buf_name) \
    bool skipped = (SECTION_SKIPPED == zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", (sec)->st, (sec))); \
    if (!skipped && (!flag.only_headers || (sec)->st == SEC_RANDOM_ACCESS)) \
        zfile_uncompress_section (evb, B1ST(SectionHeader, evb->z_data), (out_buf), (out_buf_name), 0, (sec)->st); \
    HeaderType header __attribute__((unused)) = !skipped ? *(HeaderType *)evb->z_data.data : (HeaderType){}; /* make a copy of the header */ \
    if (!skipped) buf_free (evb->z_data); 

extern bool zfile_advance_to_next_header (uint64_t *offset, uint64_t *gap);

extern DataType zfile_piz_get_file_dt (rom filename);

#ifdef __APPLE__
#define off64_t __int64_t // needed for conda mac - otherwise zlib.h throws compilation errors
#endif
