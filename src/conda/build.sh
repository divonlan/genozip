#!/bin/bash 

make -C src || exit $?

echo Installing in $PREFIX/bin
mkdir -p $PREFIX/bin 
cp -f genozip $PREFIX/bin/genozip || exit $?
chmod a+x $PREFIX/bin/genozip

for exe in genounzip genocat genols; do
    ln $PREFIX/bin/genozip $PREFIX/bin/$exe || exit $?
done
