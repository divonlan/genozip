#/bin/bash --norc

github_version=`curl -Ls -o /dev/null -w %{url_effective} https://github.com/divonlan/genozip/releases/latest | rev | cut -d- -f1 | rev`
curr_version=`cut -d\" -f2 version.h|head -n1`

if [ $curr_version != $github_version ] ; then
    echo Version is already incremented: in version.h \"$curr_version\" vs latest github release: \"$github_version\"
    exit 0
fi

curr_small=`echo $curr_version | cut -d. -f3`
new_small=`expr $curr_small + 1`
new_version=`echo $curr_version | cut -d. -f1-2`.$new_small
new_big=`echo $new_version | cut -d. -f1`

echo \#define GENOZIP_CODE_VERSION \"$new_version\" > version.h   # override previous
echo \#define GENOZIP_FILE_FORMAT_VERSION $new_big >> version.h
echo New version is $new_version

