#!/bin/bash --norc

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
