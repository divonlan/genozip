// ------------------------------------------------------------------
//   zfile.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ZFILE_INCLUDED
#define ZFILE_INCLUDED

#include "genozip.h"
#include "sections.h"
#include "dict_id.h"
#include "md5.h"
#include "compressor.h"

extern int16_t zfile_read_genozip_header (Md5Hash *digest);
extern void zfile_compress_genozip_header (const Md5Hash *single_component_md5);
extern bool zfile_get_genozip_header (uint64_t *uncompressed_data_size, uint32_t *num_samples,
                                      uint64_t *num_items_concat, Md5Hash *md5_hash_concat, char *created, unsigned created_len);

extern void zfile_compress_section_data_alg (VBlockP vb, SectionType section_type, 
                                             BufferP section_data, CompGetLineCallback callback, uint32_t total_len, 
                                             CompressionAlg comp_alg);
#define zfile_compress_section_data(vb, section_type, section_data) \
    zfile_compress_section_data_alg ((vb), (section_type), (section_data), NULL, 0, COMP_BZ2)

typedef enum {DICTREAD_ALL, DICTREAD_CHROM_ONLY, DICTREAD_EXCEPT_CHROM} ReadChromeType;
extern void zfile_read_all_dictionaries (uint32_t last_vb_i /* 0 means all VBs */, ReadChromeType read_chrom);

extern void zfile_compress_dictionary_data (VBlockP vb, MtfContextP ctx, 
                                            uint32_t num_words, const char *data, uint32_t num_chars);
extern void zfile_compress_b250_data (VBlockP vb, MtfContextP ctx);

// returns offset of header within data, EOF if end of file (or end of VCF component in the case of flag_split)
#define SEEK_NONE ((uint64_t)-1)
#define NO_SB_I ((uint32_t)-1)
extern int zfile_read_section (VBlockP vb, uint32_t original_vb_i, uint32_t sb_i, /* NO_SB_I if not a sample related section */
                               BufferP data /* buffer to append */, const char *buf_name,
                               unsigned header_size, SectionType expected_sec_type,
                               ConstSectionListEntryP sl); 

extern void zfile_uncompress_section (VBlockP vb, void *section_header, 
                                      BufferP uncompressed_data, const char *uncompressed_data_buf_name,
                                      SectionType expected_section_type);

extern void zfile_show_header (const SectionHeader *header, VBlockP vb /* optional if output to buffer */);

extern bool zfile_is_skip_section (void *vb, SectionType st, DictIdType dict_id);

extern void zfile_write_txt_header (BufferP vcf_header_text, bool is_first_vcf);
extern bool zfile_update_txt_header_section_header (uint64_t pos_of_current_vcf_header, uint32_t max_lines_per_vb, Md5Hash *md5);

#ifdef __APPLE__
#define off64_t __int64_t // needed for conda mac - otherwise zlib.h throws compilation errors
#endif

// -----------------------------
// VCF stuff
// -----------------------------
extern void zfile_vcf_compress_vb_header (VBlockP vb);
extern void zfile_vcf_update_compressed_vb_header (VBlockP vb, uint32_t vcf_first_line_i);
extern void zfile_vcf_compress_haplotype_data_gtshark (VBlockVCFP vb, ConstBufferP haplotype_sections_data, unsigned sb_i);
extern void zfile_vcf_read_one_vb (VBlockP vb);

// -----------------------------
// v1 compatibility (VCF only)
// -----------------------------

extern bool v1_zfile_vcf_read_one_vb (VBlockVCFP vb);
extern int v1_zfile_read_section (VBlockP vb, BufferP data, const char *buf_name, unsigned header_size, SectionType expected_sec_type, bool allow_eof);

// -----------------------------
// SAM stuff
// -----------------------------

extern void zfile_sam_read_one_vb (VBlockP vb);
extern void zfile_compress_generic_vb_header (VBlockP vb);
extern void zfile_update_compressed_vb_header (VBlockP vb, uint32_t vcf_first_line_i);

// -----------------------------
// FASTQ/FASTA stuff
// -----------------------------

extern void zfile_fast_read_one_vb (VBlockP vb);
extern void zfile_fast_compress_vb_header (VBlockP vb);
extern void zfile_fast_update_compressed_vb_header (VBlockP vb, uint32_t vcf_first_line_i);

// -----------------------------
// ME23 stuff
// -----------------------------

extern void zfile_me23_read_one_vb (VBlockP vb);
extern void zfile_me23_compress_vb_header (VBlockP vb);
extern void zfile_me23_update_compressed_vb_header (VBlockP vb, uint32_t vcf_first_line_i);

#endif
