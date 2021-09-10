// ------------------------------------------------------------------
//   zfile.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef ZFILE_INCLUDED
#define ZFILE_INCLUDED

#include "genozip.h"
#include "sections.h"
#include "digest.h"

// --------
// ZIP side
// --------

extern void zfile_compress_genozip_header (Digest single_component_digest);

extern void zfile_compress_section_data_ex (VBlockP vb, SectionType section_type, 
                                            BufferP section_data, LocalGetLineCB callback, uint32_t total_len, 
                                            Codec codec, SectionFlags flags);

#define zfile_compress_section_data(vb, section_type, section_data) \
    zfile_compress_section_data_ex ((vb), (section_type), (section_data), NULL, 0, CODEC_BZ2, SECTION_FLAGS_NONE)

extern uint32_t zfile_compress_b250_data  (VBlockP vb, ContextP ctx);
extern uint32_t zfile_compress_local_data (VBlockP vb, ContextP ctx, uint32_t sample_size);

extern LocalGetLineCB *zfile_get_local_data_callback (DataType dt, ContextP ctx);

extern void zfile_compress_vb_header (VBlockP vb);
extern void zfile_update_compressed_vb_header (VBlockP vb);

extern void zfile_write_txt_header (BufferP vcf_header_text, uint64_t unmodified_txt_header_len, Digest header_md5, bool is_first_vcf);
extern bool zfile_update_txt_header_section_header (uint64_t pos_of_current_vcf_header, uint32_t max_lines_per_vb, Digest *md5);

extern void zfile_output_processed_vb (VBlockP vb);

// --------
// PIZ side
// --------

extern bool zfile_read_genozip_header (SectionHeaderGenozipHeader *header);

extern SectionHeaderUnion zfile_read_section_header (VBlockP vb, uint64_t offset, uint32_t original_vb_i, SectionType expected_sec_type);

#define SECTION_SKIPPED ((int32_t)-1)
extern int32_t zfile_read_section_do (FileP file, VBlockP vb, uint32_t original_vb_i, 
                                      BufferP data /* buffer to append */, const char *buf_name,
                                      SectionType expected_sec_type, 
                                      Section sl, uint32_t header_size, const char *func, uint32_t code_line); 
#define zfile_read_section(file,vb,original_vb_i,data,buf_name,expected_sec_type,sl) \
    zfile_read_section_do ((file),(VBlockP)(vb),(original_vb_i),(data),(buf_name),(expected_sec_type),(sl), st_header_size (expected_sec_type), __FUNCTION__, __LINE__)

extern void zfile_uncompress_section (VBlockP vb, void *section_header, 
                                      BufferP uncompressed_data, 
                                      const char *uncompressed_data_buf_name,
                                      uint32_t expected_vb_i, SectionType expected_section_type);

#define zfile_get_global_section(HeaderType, sec,sl,out_buf,out_buf_name) \
    zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", (sec), (sl)); \
    if (!flag.show_headers || (sec) == SEC_RANDOM_ACCESS) zfile_uncompress_section (evb, evb->z_data.data, (out_buf), (out_buf_name), 0, (sec)); \
    HeaderType header __attribute__((unused)) = *(HeaderType *)evb->z_data.data; /* make a copy of the header */ \
    buf_free (&evb->z_data); 

extern DataType zfile_get_file_dt (const char *filename);

// --------
// General
// --------

extern void zfile_show_header (const SectionHeader *header, VBlockP vb /* optional if output to buffer */, uint64_t offset, char rw);

#ifdef __APPLE__
#define off64_t __int64_t // needed for conda mac - otherwise zlib.h throws compilation errors
#endif

#endif
