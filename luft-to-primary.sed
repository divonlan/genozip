#!/bin/sed -Erf

# ------------------------------------------------------------------
#   luft-to-primary.sed
#   Copyright (C) 2021-2021 Divon Lan <divon@genozip.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
#
#   Convert a LUFT VCF to a PRIMARY VCF (just renaming stuff, no liftover). Used for testing during development.
#
#   See also: primary-to-luft.sed

s/##dual_coordinates=LUFT/##dual_coordinates=PRIMARY/g
s/##primary_/##luft_/g
s/LIFTBACK=/LIFTOVER=/g
s/REJTBACK=/REJTOVER=/g
s/^(##luft_only=.*)REJTOVER(.*)$/\1REJTBACK\2/g  # note: ##primary_only already replaced to ##luft_only
