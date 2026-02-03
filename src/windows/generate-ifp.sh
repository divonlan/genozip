#/bin/bash --norc

# ------------------------------------------------------------------
#   generate-ifp.sh
#   Copyright (C) 2024-2026 Genozip Limited. Patent Pending.
#   Please see terms and conditions in the file LICENSE.txt
#
#   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
#   and subject to penalties specified in the license.

# This script is invoked by Makefile

if [[ `pwd` != *windows* ]]; then
    echo "$0 must be run from the src/windows directory"
    exit 1
fi

template=genozip-installer.template.ifp

if [ ! -f $template ]; then
    echo "$template: file not found"
    exit 1
fi

version=$(head -n1 ../version.h |cut -d\" -f2)

lines=`wc -l $template | cut -d" " -f1`
head=$(( `grep -n __LICENSE__ $template | cut -d: -f1` - 1 )) # __LICENSE__ must be on a stand-alone line
tail=$(( $lines - $head - 1 ))
windows_dir=`echo $PWD | sed "s/^\/mnt\/c/C:/g" | sed "s/\//\\\\\\\\\\\\\\\\/g"`

cd ../..
root_dir=`echo $PWD | sed "s/^\/mnt\/c/C:/g" | sed "s/\//\\\\\\\\\\\\\\\\/g"`
cd - > /dev/null

cd ../../installers
installers_dir=`echo $PWD | sed "s/^\/mnt\/c/C:/g" | sed "s/\//\\\\\\\\\\\\\\\\/g"`
cd - > /dev/null

exe_mb=5.4 # not sure if this needs to be accurate
license_kb=11.0

rm -f genozip-installer.ifp

# NOTE: if making any change, verify that the Registry section shows up as expected in 
# InstallForge. It is very sensitive to any change in the License section.

(head -$head $template ; sed "s/'/\\\\f1\\\\rquote\\\\f0 /g" LICENSE.for-installer.txt | sed "s/$/\\\\par\r/g"; tail -$tail $template) \
| sed s/__VERSION__/${version}/g \
| sed "s/__ROOT_DIR__/${root_dir}/g" \
| sed "s/__WINDOWS_DIR__/${windows_dir}/g" \
| sed "s/__INSTALLERS__/${installers_dir}/g" \
| sed "s/__EXE_MB__/${exe_mb}/g" \
| sed "s/__LICENSE_KB__/${license_kb}/g" \
> genozip-installer.ifp
