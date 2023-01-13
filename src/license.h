// ------------------------------------------------------------------
//   license.h
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

typedef enum __attribute__ ((__packed__))/*1 byte*/ { LIC_TYPE_NONE, LIC_TYPE_ACADEMIC, LIC_TYPE_EVAL, LIC_TYPE_STANDARD, LIC_TYPE_DEEP, NUM_LIC_TYPES } LicenseType; 

extern void license_register (void);
extern bool license_is_registered (void);
extern void license_set_filename (rom filename);
extern void license_load (void);
extern uint32_t license_get_number (void);
extern rom license_get_one_line (void);
extern void license_display (void);
extern bool license_allow_stats (void);
extern void license_one_file_compressed (DataType dt);
extern void license_show_deep_notice (void);

typedef struct { char s[16384]; } StrNotice; 
extern StrNotice license_print_default_notice (void);

