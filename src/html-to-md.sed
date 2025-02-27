#!/bin/sed -Ef

# ------------------------------------------------------------------
#   html-to-md.sed
#   Copyright (C) 2020-2025 Genozip Limited. Patent Pending.
#   Please see terms and conditions in the file LICENSE.txt
#
#   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
#   and subject to penalties specified in the license.

s/<br>$/  /g
s/<b>/\*\*/g
s/<\/b>/\*\*/g
s/<i>/\*/g
s/<\/i>/\*/g
s/<h1>//g
s/<\/h1>/\n=======/g
s/&nbsp/ /g
