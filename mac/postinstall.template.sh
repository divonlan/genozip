#!/bin/bash --norc

# ------------------------------------------------------------------
#   postinstall
#   Copyright (C) 2020-2021 Black Paw Ventures Limited 
#   Please see terms and conditions in the file LICENSE.txt

FILES=(__FILES__)

[ -d /usr/local/bin ] || mkdir /usr/local/bin

chmod -R 755 /Library/genozip

for file in "${FILES[@]}"
do
    rm -f /usr/local/bin/${file}
    ln -s /Library/genozip/${file} /usr/local/bin/${file}
done

