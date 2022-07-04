// ------------------------------------------------------------------
//   license.h
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

typedef enum __attribute__ ((__packed__))/*1 byte*/ { LIC_TYPE_NONE, LIC_TYPE_ACADEMIC, LIC_TYPE_EVAL, LIC_TYPE_PAID, NUM_LIC_TYPES } LicenseType; 

extern void license_register (void);
extern void license_set_filename (rom filename);
extern uint32_t license_get_number (void);
extern rom license_get_one_line (void);
extern void license_display (void);
extern LicenseType license_get_type (void);
extern void license_print_tip (void);
