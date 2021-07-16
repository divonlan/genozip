#!/bin/sed -Ef

# ------------------------------------------------------------------
#   html-to-md.sed
#   Copyright (C) 2020-2021 Black Paw Ventures Limited
#   Please see terms and conditions in the file LICENSE.txt

s/<br>$/  /g
s/<b>/\*\*/g
s/<\/b>/\*\*/g
s/<i>/\*/g
s/<\/i>/\*/g
s/<h1>//g
s/<\/h1>/\n=======/g
s/&nbsp/ /g
