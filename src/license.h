// ------------------------------------------------------------------
//   license.h
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

typedef enum __attribute__ ((__packed__))/*1 byte*/ { LIC_TYPE_NONE, LIC_TYPE_ACADEMIC, LIC_TYPE_EVAL, LIC_TYPE_STANDARD, NUM_LIC_TYPES } LicenseType; 

extern void license_register (void);
extern bool license_is_registered (void);
extern void license_set_filename (rom filename);
extern uint32_t license_get_number (void);
extern rom license_get_one_line (void);
extern void license_display (void);
extern LicenseType license_get_type (void);
extern bool license_allow_stats (void);
extern rom license_print_default_notice (void);
extern void license_one_file_compressed (void);