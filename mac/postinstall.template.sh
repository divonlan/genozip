#!/bin/bash --norc

# ------------------------------------------------------------------
#   postinstall
#   Copyright (C) 2020-2022 Genozip Limited 
#   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

FILES=(__FILES__)

[ -d /usr/local/bin ] || mkdir /usr/local/bin

chmod -R 755 /Library/genozip

for file in "${FILES[@]}"
do
    rm -f /usr/local/bin/${file}
    ln -s /Library/genozip/${file} /usr/local/bin/${file}
done

