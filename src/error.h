// ------------------------------------------------------------------
//   error.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

extern void error_initialize (int argc, char *argv[]);
extern void error_print_call_stack (void);
