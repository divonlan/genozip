// ------------------------------------------------------------------
//   tar.h
//   Copyright (C) 2021-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"

extern void tar_set_tar_name (rom tar_filename);
extern rom tar_get_tar_name (void);
extern void tar_initialize (BufferP input_files_buf);
extern FILE *tar_open_file (rom fn_on_disk, rom fn_in_tar);
extern void tar_copy_file (rom fn_on_disk, rom fn_in_tar);
extern bool tar_zip_is_tar (void);
extern int64_t tar_file_offset (void);
extern void tar_close_file (void **file);
extern void tar_finalize (void);
