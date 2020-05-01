#!/bin/bash

output=test-output

files=(test-file.vcf test-file.sam test-file.fq test-file.fa genome_23andme_Full_test-file.txt)
for file in ${files[@]}; do
    printf "\nTESTING $file - basic test - Unix-style end-of-line\n"
    cat $file | tr -d "\r" > unix-nl.$file
    ./genozip unix-nl.$file -ft -o ${output}.genozip || exit 1

    printf "\nTESTING $file - Window-style end-of-line\n"
    sed.exe 's/$/\r/g' unix-nl.$file > windows-nl.$file
    ./genozip windows-nl.$file -ft -o ${output}.genozip || exit 1
    rm unix-nl.$file windows-nl.$file

    printf "\nTESTING $file - encrypted\n"
    ./genozip $file --password abc -ft -o ${output}.genozip || exit 1

    printf "\nTESTING $file - redirected from stdin\n"
    cat $file | ./genozip --test --force --output ${output}.genozip --input-type ${file#*.} - || exit 1

    printf "\nTESTING $file - concat & split\n"
    file1=copy1.$file
    file2=copy2.$file
    cp $file $file1
    cp $file $file2
    ./genozip $file1 $file2 -ft -o ${output}.genozip || exit 1
    ./genounzip ${output}.genozip -O -t
    ls  cop*
    rm $file1 $file2

    printf "\nTESTING $file --optimize - NOT checking correctness, just that it doesn't crash\n"
    ./genozip $file -f --optimize -o ${output}.genozip || exit 1

    printf "\nTESTING $file --strip - NOT checking correctness, just that it doesn't crash\n"
    ./genocat ${output}.genozip --strip > /dev/null || exit 1
done

if `command -v md5 >& /dev/null`; then
    md5=md5 # mac
else
    md5=md5sum 
fi

echo "Backward compatability tests"
files=`ls backward-compatibility-test/*.vcf` 
for file in $files; do
    ./genounzip ${file}.genozip -fo $output
    if [[ `$md5 $file | cut -d" " -f1` != `$md5 ${output%.*} | cut -d" " -f1` ]]; then
        echo FAILED!
        exit 1
    fi
done
 
echo " "
echo "ALL GOOD!"

rm -f $output ${output}.genozip
