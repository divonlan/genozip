// ------------------------------------------------------------------
//   zfile.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ZFILE_INCLUDED
#define ZFILE_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "sections.h"
#include "dict_id.h"

extern void zfile_write_vcf_header (VariantBlockP vb, Buffer *vcf_header_text, bool is_first_vcf);
extern void zfile_compress_variant_data (VariantBlockP vb);
extern void zfile_update_compressed_variant_data_header (VariantBlockP vb, unsigned pos, unsigned num_dictionary_sections);
extern void zfile_compress_section_data (VariantBlockP vb, SectionType section_type, Buffer *section_data);
extern void zfile_compress_dictionary_data (VariantBlockP vb, DictIdType dict_id, 
                                            uint32_t num_words, const char *data, uint32_t num_chars);

extern bool zfile_read_one_vb (VariantBlockP vb);

// returns offset of header within data, EOF if end of file (or end of VCF component in the case of flag_split)
extern int zfile_read_one_section (VariantBlockP vb, 
                                   Buffer *data /* buffer to append */, const char *buf_name,
                                   unsigned header_size, SectionType expected_sec_type,
                                   bool allow_eof);

extern void zfile_uncompress_section (VariantBlockP vb, void *section_header, Buffer *uncompressed_data, SectionType expected_section_type);

#ifdef __APPLE__
#define off64_t __int64_t // needed for conda mac - otherwise zlib.h throws compilation errors
#endif
extern void zfile_update_vcf_header_section_header (VariantBlockP vb, off64_t vcf_header_header_pos_single, bool final_for_concat);

#endif
