#/bin/bash --norc

# This script is invoked by Makefile

(  sed /__README_MD__/Q conda/README.template.md ; \
   ./html-to-md.sed ../README.md | grep -v "<!" ; \
   sed '0,/__README_MD__/d' conda/README.template.md \
   ) | \
   grep -v "^#" 
   