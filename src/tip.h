// ------------------------------------------------------------------
//   tip.h
//   Copyright (C) 2022-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern void tip_dt_encountered (DataType dt);
extern void tip_print (void);

// protect from spamming the user with more than one tip
// Note: use this for learning-curve tips, not for failure cases which should always display
#define TIP(format, ...)  ( { if (!flag.no_tip) iprintf ("Tip: "format "\n", __VA_ARGS__); flag.no_tip = true; } ) 
