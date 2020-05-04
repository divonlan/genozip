#!/bin/bash

output=test-output

is_windows=`uname|grep -i mingw`

# -----------------
# platform settings
# -----------------
if [ -n "$is_windows" ]; then
    path=`pwd| cut -c3-|tr / '\\\\'`\\
else
    path=$PWD/
fi

if `command -v md5 >& /dev/null`; then
    md5=md5 # mac
else
    md5=md5sum 
fi

cmp_2_files() {
    if [[ `$md5 $1 | cut -d" " -f1` != `$md5 ${2%.*} | cut -d" " -f1` ]]; then
        echo FAILED!
        exit 1
    fi
}

test_header() {
    printf "\n=================================================================\n"
    printf   "TESTING $1 \n"
    printf   "=================================================================\n"
}

files=(test-file.vcf test-file.sam test-file.fq test-file.fa genome_23andme_Full_test-file.txt)
for file in ${files[@]}; do
    test_header "$file - basic test - Unix-style end-of-line"
    cat $file | tr -d "\r" > unix-nl.$file
    ./genozip unix-nl.$file -ft -o ${output}.genozip || exit 1

    test_header "$file - Window-style end-of-line"
    sed 's/$/\r/g' unix-nl.$file > windows-nl.$file
    ./genozip windows-nl.$file -ft -o ${output}.genozip || exit 1
    rm unix-nl.$file windows-nl.$file

    test_header "$file - as URL"
    ./genozip file://${path}$file -ft -o ${output}.genozip || exit 1

    test_header "$file - encrypted"
    ./genozip $file --password abc -ft -o ${output}.genozip || exit 1

    test_header "$file - redirected from stdin"
    cat $file | ./genozip --test --force --output ${output}.genozip --input-type ${file#*.} - || exit 1

    if [ -z "$is_windows" ]; then # windows can't redirect binary data
        test_header "$file - redirecting stdout"
        ./genozip ${file} --stdout > ${output}.genozip || exit 1
        ./genounzip ${output}.genozip -f || exit 1
        cmp_2_files $file $output
    fi

    test_header "$file - concat & split"
    file1=copy1.$file
    file2=copy2.$file
    cp $file $file1
    cp $file $file2
    ./genozip $file1 $file2 -ft -o ${output}.genozip || exit 1
    ./genounzip ${output}.genozip -O -t || exit 1
    ls  cop*
    rm $file1 $file2

    test_header "$file --optimize - NOT checking correctness, just that it doesn't crash"
    ./genozip $file -f --optimize -o ${output}.genozip || exit 1

    test_header "$file --strip - NOT checking correctness, just that it doesn't crash"
    ./genocat ${output}.genozip --strip > /dev/null || exit 1
done

echo "Backward compatability tests"
files=`ls backward-compatibility-test/*.vcf` 
for file in $files; do
    ./genounzip ${file}.genozip -fo $output || exit 1
    cmp_2_files $file $output
    
done
 
echo " "
echo "ALL GOOD!"

rm -f $output ${output}.genozip
