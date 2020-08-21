#!/bin/sed -Ef

# ------------------------------------------------------------------
#   html-to-md.sed
#   Copyright (C) 2020 Divon Lan <divon@genozip.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

s/<br>$/  /g
s/<b>/\*\*/g
s/<\/b>/\*\*/g
s/<i>/\*/g
s/<\/i>/\*/g
s/<h1>//g
s/<\/h1>/\n=======/g
s/&nbsp/ /g
