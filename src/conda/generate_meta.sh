#/bin/bash --norc

# ------------------------------------------------------------------
#   generate-meta.sh
#   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
#   Please see terms and conditions in the file LICENSE.txt
#
#   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
#   and subject to penalties specified in the license.

# This script is invoked by Makefile

version=$(head -n1 version.h |cut -d\" -f2)

(  sed /__README_MD__/Q conda/meta.template.yaml ; \
    cat ../README.md | sed "s/<.\{1,3\}>//g"| sed "s/\&nbsp/ /g" | grep -v "<!" | sed 's/^/    /' ; \
    sed '0,/__README_MD__/d' conda/meta.template.yaml \
    ) | \
    sed s/__SHA256__/$(openssl sha256 .archive.tar.gz | cut -d= -f2 | cut -c2-)/ | \
    sed s/__VERSION__/${version}/g | \
    grep -v "^#" 
