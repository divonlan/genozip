// ------------------------------------------------------------------
//   license.h
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

extern void license_register (void);
extern void license_set_filename (const char *filename);
extern uint32_t license_get_number (void);
extern const char *license_get_one_line (void);
extern bool license_has_details (void);
extern void license_display (void);

