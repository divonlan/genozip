// ------------------------------------------------------------------
//   tar.h
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

extern void tar_set_tar_name (rom tar_filename);
extern void tar_initialize (void);
extern FILE *tar_open_file (rom z_filename);
extern bool tar_is_tar (void);
extern int64_t tar_file_offset (void);
extern void tar_close_file (void **file);
extern void tar_copy_file (rom z_filename);
extern void tar_finalize (void);
