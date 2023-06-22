// ------------------------------------------------------------------
//   license.h
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern void license_register (void);
extern bool license_is_registered (void);
extern void license_set_filename (rom filename);
extern void license_load (void);
extern rom license_get_one_line (void);
extern void license_display (void);
extern bool license_allow_tip (void);
extern void license_one_file_compressed (DataType dt);
extern void license_show_deep_notice (void);
extern bool license_is_eval (void);

typedef struct { char s[16384]; } StrNotice; 
extern StrNotice license_print_default_notice (void);

