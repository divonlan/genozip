#!/bin/sed -Erf

# ------------------------------------------------------------------
#   primary-to-luft.sed
#   Copyright (C) 2021-2022 Black Paw Ventures Limited
#   Please see terms and conditions in the file LICENSE.txt
#
#   Convert a PRIMARY VCF to a LUFT VCF (just renaming strings, no liftover). Used for testing during development.
#
#   See also: luft-to-primary.sed

s/##dual_coordinates=PRIMARY/##dual_coordinates=LUFT/g
s/##luft_/##primary_/g
s/LIFTOVER=/LIFTBACK=/g
s/REJTOVER=/REJTBACK=/g
s/^(##primary_only=.*)REJTBACK(.*)$/\1REJTOVER\2/g  # note: ##luft_only already replaced to ##primary_only

