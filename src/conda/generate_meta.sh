#/bin/bash --norc

# This script is invoked by Makefile

version=$(head -n1 version.h |cut -d\" -f2)

(  sed /__README_MD__/Q conda/meta.template.yaml ; \
    cat README.md | sed "s/<.\{1,3\}>//g"| sed "s/\&nbsp/ /g" | grep -v "<!" | sed 's/^/    /' ; \
    sed '0,/__README_MD__/d' conda/meta.template.yaml \
    ) | \
    sed s/__SHA256__/$(openssl sha256 .archive.tar.gz | cut -d= -f2 | cut -c2-)/ | \
    sed s/__VERSION__/${version}/g | \
    grep -v "^#" 
