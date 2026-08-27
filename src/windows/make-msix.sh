#/bin/bash --norc

# ------------------------------------------------------------------
#   generate-manifest.sh
#   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
#   Please see terms and conditions in the file LICENSE.txt
#
#   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
#   and subject to penalties specified in the license.

# This script is invoked by src/Makefile

if [[ -z "$GENOZIP_HOME" ]]; then # definition in /home/divon/.bashrc overrides definition in Windows Settings->Environment Variables
    echo "GENOZIP_HOME is not set"
    exit 1
fi

if [[ `uname -o` != Msys ]]; then
    echo "This script must be run from Windows (Msys terminal)"
    exit 1
fi

cd $GENOZIP_HOME/src
make -j windows/genozip.exe

cd $GENOZIP_HOME/src/windows

template=AppxManifest.template.xml

if [ ! -f $template ]; then
    echo "$template: file not found"
    exit 1
fi

version=$(head -n1 ../version.h |cut -d\" -f2)

windows_build=$(powershell.exe -NoProfile -Command "[System.Environment]::OSVersion.Version.Build" 2>/dev/null)

rm -f msix/*
mkdir -p msix

cat $template \
| sed s/__VERSION__/${version}/g \
| sed s/__WINDOWS_BUILD__/${windows_build}/g \
> msix/AppxManifest.xml

# copy from windows dir (made with -DDISTRIBUTION=MicrosoftStore)
cp genozip.exe msix/genozip.exe
cp genozip.exe msix/genounzip.exe
cp genozip.exe msix/genocat.exe
cp genozip.exe msix/genols.exe
cp ../../LICENSE.txt msix/
cp genozip.png StoreLogo.png Square150x150Logo.png Square44x44Logo.png msix/

winapp pack ./msix
