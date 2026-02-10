// ------------------------------------------------------------------
//   license.h
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "sections.h"

extern void license_activate (void);
extern bool license_is_activated (void);
extern void license_set_filename (rom filename);
extern void license_load (void);
extern StrText license_get_number (void);
extern rom lic_type_name (int8_t lic_type);
extern rom license_get_one_line (void);
extern void license_display (bool show_latest);
extern bool license_allow_tip (void);
extern bool license_allow_distribution (void);
extern void license_eval_notice (void);
extern bool license_is_eval (void);
extern bool license_is_standard (void);
extern bool license_is_enterprise (void);
extern bool license_is_academic (void);
extern void license_prepare (rom arg);
extern bool license_piz_prepare_genozip_header (SectionHeaderGenozipHeaderP header, FailType fail_type);
extern StrTextLong license_academic_tip (void);
extern bool am_i_submit (void);