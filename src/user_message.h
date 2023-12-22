// ------------------------------------------------------------------
//   user_message.h
//   Copyright (C) 2023-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern void user_message_init (rom filename);
extern void user_message_compress (void);
extern void user_message_display (void);
extern bool has_user_msg (void);
