#/bin/bash --norc

old_version=`cut -d\" -f2 version.h|head -n1`
old_small=`echo $old_version | cut -d. -f3`
new_small=`expr $old_small + 1`
new_version=`echo $old_version | cut -d. -f1-2`.$new_small
new_big=`echo $new_version | cut -d. -f1`

echo \#define GENOZIP_CODE_VERSION \"$new_version\" > version.h   # override previous
echo \#define GENOZIP_FILE_FORMAT_VERSION $new_big >> version.h
echo New version is $new_version

