#!/bin/bash --norc

# ------------------------------------------------------------------
#   uninstall.sh
#   Copyright (C) 2020-2024 Genozip Limited where applies
#   Please see terms and conditions in the file LICENSE.txt
#
#   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
#   and subject to penalties specified in the license.

# loosely based on https://github.com/KosalaHerath/macos-installer-builder which is licensed under Apache 2.0 license

if (( $EUID != 0 )); then
    echo "Please run as root."
    exit
fi

while true; do
    read -p "Do you wish to uninstall genozip [Y/n]?" answer
    [[ $answer == "y" || $answer == "Y" || $answer == "" ]] && break
    [[ $answer == "n" || $answer == "N" ]] && exit 0
    echo "Please answer with 'y' or 'n'"
done


(cd /usr/local/bin; rm -f __FILES__)
 
pkgutil --forget "org.genozip.__VERSION__" > /dev/null 2>&1

[ -e "/Library/genozip" ] && rm -rf "/Library/genozip"

echo Done
exit 0
