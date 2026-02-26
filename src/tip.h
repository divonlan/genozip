// ------------------------------------------------------------------
//   tip.h
//   Copyright (C) 2022-2026 Genozip Limited. Patent Pending.
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
#define TIP(format, ...)  ( { if (!flag.no_tip) iprintf ("\n\n" _TIP format "\n\n", __VA_ARGS__); flag.no_tip = true; } ) 
#define TIP0(str)  ( { if (!flag.no_tip) iprint0 ("\n\n" _TIP str "\n\n"); flag.no_tip = true; } ) 

// strings for tips to ensure consistency (italics from https://en.wikipedia.org/wiki/Mathematical_Alphanumeric_Symbols)
#define _TIP "💡 " 
#define _REFFILE "𝑟𝑒𝑓-𝑓𝑖𝑙𝑒" // better than sans-serif: "𝘳𝘦𝘧-𝘧𝘪𝘭𝘦"
#define _BAMFILE "𝑏𝑎𝑚-𝑓𝑖𝑙𝑒"
