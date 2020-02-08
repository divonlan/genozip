#!/bin/bash --norc

# loosely based on https://github.com/KosalaHerath/macos-installer-builder which is licensed under Apache 2.0 license

# reference: http://thegreyblog.blogspot.com/2014/06/os-x-creating-packages-from-command_2.html

MAC_DIR=mac-pkg
TARGET_DIR=${MAC_DIR}/target
PRODUCT=genozip
VERSION=`head -n1 version.h |cut -d\" -f2`
FILES=(genozip genounzip genols genocat) # array
FILES_STR=${FILES[@]} # string

# create installation directory
rm -rf $TARGET_DIR
mkdir $TARGET_DIR
if [[ $? != 0 ]]; then
    echo "Failed to create $TARGET_DIR directory" $?
    exit 1
fi
mkdir -p ${TARGET_DIR}/darwinpkg/Library/genozip ${TARGET_DIR}/Resources ${TARGET_DIR}/scripts 

# copy and adjust files
cp ${FILES[@]} ${TARGET_DIR}/darwinpkg/Library/genozip
cp ${MAC_DIR}/Distribution ${TARGET_DIR}
cp ${MAC_DIR}/welcome.html ${MAC_DIR}/conclusion.html ${MAC_DIR}/banner.png ${MAC_DIR}/uninstall.sh ${TARGET_DIR}/Resources 
sed -e "s/__FILES__/${FILES_STR}/g" ${MAC_DIR}/postinstall > ${TARGET_DIR}/scripts/postinstall

cp ${MAC_DIR}/uninstall.sh ${TARGET_DIR}/darwinpkg/Library/genozip
sed -i '' -e "s/__VERSION__/${VERSION}/g" "${TARGET_DIR}/darwinpkg/Library/genozip/uninstall.sh"
sed -i '' -e "s/__FILES__/${FILES_STR}/g" "${TARGET_DIR}/darwinpkg/Library/genozip/uninstall.sh"

chmod -R 755 ${TARGET_DIR}

# build package
echo Building package
pkgbuild --identifier org.genozip.${VERSION} --version ${VERSION} --scripts ${TARGET_DIR}/scripts --root ${TARGET_DIR}/darwinpkg ${MAC_DIR}/genozip.pkg > /dev/null 2>&1

# build product
echo Building product
PRODUCT=${MAC_DIR}/genozip-macos-installer-x64-${VERSION}.pkg
productbuild --distribution ${TARGET_DIR}/Distribution --resources ${TARGET_DIR}/Resources --package-path ${MAC_DIR} ${PRODUCT} > /dev/null 2>&1

# sign product - IF we have a certificate from Apple
if [ -f apple_developer_certificate_id ]; then
    echo Signing the product
    APPLE_DEVELOPER_CERTIFICATE_ID=`cat apple_developer_certificate_id`
    productsign --sign "Developer ID Installer: ${APPLE_DEVELOPER_CERTIFICATE_ID}" ${PRODUCT} ${PRODUCT}.signed
    pkgutil --check-signature ${PRODUCT}.signed
    mv -f ${PRODUCT}.signed ${PRODUCT}
fi

exit 0
