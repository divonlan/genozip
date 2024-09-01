// ------------------------------------------------------------------
//   txtfile.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "digest.h"

#define TXTFILE_READ_VB_PADDING 16 // we need this quantity of unused bytes at the end of vb.txt_data

extern uint32_t txtfile_fread (FileP file, FILE *fp, void *addr, int32_t size, int64_t *disk_so_far);
extern void txtfile_fwrite (const void *data, uint32_t size);
extern void txtfile_query_first_bytes_in_file (rom filename, uint32_t len);
extern void txtfile_initialize_igzip (FileP file);
extern StrText txtfile_dump_vb (VBlockP vb, rom base_name, BufferP txt_data);
extern StrTextLong txtfile_codec_name (FileP z_file, CompIType comp_i, bool obscure_fname);
extern void txtfile_zip_finalize_codecs (void);
extern void txtfile_read_header (bool is_first_txt);
extern uint32_t txt_data_alloc_size (uint32_t vb_size) ;
extern void txtfile_read_vblock (VBlockP vb);
extern void txtfile_set_seggable_size (void);
extern int64_t txtfile_get_seggable_size (void);
extern bool txtfile_is_gzip (FileP file);
extern void txtfile_discover_specific_gz (FileP file);
extern rom isal_error (int ret);

// callbacks
extern int32_t def_unconsumed (VBlockP vb, uint32_t first_i);
extern int32_t def_is_header_done (bool is_eof);

extern DataType txtfile_zip_get_file_dt (rom filename);

// misc
extern StrTextLong display_gz_header (STR8p(h), bool obscure_fname);