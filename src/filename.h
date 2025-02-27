// ------------------------------------------------------------------
//   filename.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "file_types.h"
#include "flags.h"

extern rom filename_z_by_flag (void);
extern rom filename_z_normal (rom txt_filename, DataType dt, FileType txt_ft);
extern char *filename_z_pair (rom fn1, rom fn2, bool test_only);
extern char *filename_z_deep (rom sam_name);

extern rom filename_guess_original (ConstFileP file);

extern bool filename_has_ext (rom filename, rom extension);
extern void filename_remove_codec_ext (char *filename, FileType ft);
extern rom filename_base (rom filename, bool remove_exe, rom default_basename, char *basename, unsigned basename_size);

extern char *filename_make_unix (char *filename);
