#!/bin/sed -Ef

# ------------------------------------------------------------------
#   html-to-md.sed
#   Copyright (C) 2020-2022 Genozip Limited
#   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

s/<br>$/  /g
s/<b>/\*\*/g
s/<\/b>/\*\*/g
s/<i>/\*/g
s/<\/i>/\*/g
s/<h1>//g
s/<\/h1>/\n=======/g
s/&nbsp/ /g
