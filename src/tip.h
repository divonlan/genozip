// ------------------------------------------------------------------
//   tip.h
//   Copyright (C) 2022-2024 Genozip Limited. Patent Pending.
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
#define TIP(format, ...)  ( { if (!flag.no_tip) iprintf ("\nTip: "format "\n\n", __VA_ARGS__); flag.no_tip = true; } ) 
#define TIP0(str)  ( { if (!flag.no_tip) iprint0 ("\nTip: " str "\n\n"); flag.no_tip = true; } ) 
