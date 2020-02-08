#!/bin/bash --norc

# ------------------------------------------------------------------
#   postinstall
#   Copyright (C) 2020 Divon Lan <divon@genozip.com> 
#   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

FILES=(__FILES__)

[ -d /usr/local/bin ] || mkdir /usr/local/bin

chmod -R 755 /Library/genozip

for file in "${FILES[@]}"
do
    rm -f /usr/local/bin/${file}
    ln -s /Library/genozip/${file} /usr/local/bin/${file}
done

