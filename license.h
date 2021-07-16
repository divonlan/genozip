// ------------------------------------------------------------------
//   license.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef LICENSE_INCLUDED
#define LICENSE_INCLUDED

extern void license_register (void);
extern void license_set_filename (const char *filename);
extern uint32_t license_get_number (void);
extern const char *license_get_one_line (void);
extern bool license_has_details (void);
extern void license_display (void);
#endif
