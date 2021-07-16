// ------------------------------------------------------------------
//   tar.h
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef TAR_INCLUDED
#define TAR_INCLUDED

#include "genozip.h"

extern void tar_set_tar_name (const char *tar_filename);
extern void tar_initialize (void);
extern FILE *tar_open_file (const char *z_filename);
extern bool tar_is_tar (void);
extern int64_t tar_file_offset (void);
extern void tar_close_file (void **file);
extern void tar_copy_file (const char *z_filename);
extern void tar_finalize (void);

#endif
