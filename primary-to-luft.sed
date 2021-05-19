#!/bin/sed -Erf

# ------------------------------------------------------------------
#   primary-to-luft.sed
#   Copyright (C) 2021-2021 Divon Lan <divon@genozip.com>
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
#
#   Convert a PRIMARY VCF to a LUFT VCF (just renaming strings, no liftover). Used for testing during development.
#
#   See also: luft-to-primary.sed

s/##dual_coordinates=PRIMARY/##dual_coordinates=LUFT/g
s/##luft_/##primary_/g
s/LIFTOVER=/LIFTBACK=/g
s/REJTOVER=/REJTBACK=/g
s/^(##primary_only=.*)REJTBACK(.*)$/\1REJTOVER\2/g  # note: ##luft_only already replaced to ##primary_only

