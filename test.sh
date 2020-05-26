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
        echo "MD5 comparison FAILED!"
        exit 1
    fi
}

test_header() {
    sep="=======================================================================================================\n"
    printf "\n${sep}TESTING $1 \n${sep}"
}

files=(test-file.vcf test-file.sam test-file.fq test-file.fa test-file.gvf genome_23andme_Full_test-file.txt)
for file in ${files[@]}; do
    test_header "$file - basic test - Unix-style end-of-line"
    cat $file | tr -d "\r" > unix-nl.$file
    ./genozip unix-nl.$file -ft -o ${output}.genozip || exit 1

    test_header "$file - Window-style end-of-line"
    sed 's/$/\13/g' unix-nl.$file > windows-nl.$file # note: sed on mac doesn't recognize \r
    ./genozip windows-nl.$file -ft -o ${output}.genozip || exit 1
    rm unix-nl.$file windows-nl.$file

    test_header "$file - as URL"
    ./genozip file://${path}$file -ft -o ${output}.genozip || exit 1

    test_header "$file - encrypted"
    ./genozip $file --password abc -ft -o ${output}.genozip || exit 1

    test_header "$file - redirected from stdin"
    cat $file | ./genozip --test --force --output ${output}.genozip --input-type ${file#*.} - || exit 1

    if [ $file != test-file.sam ] && [ $file != genome_23andme_Full_test-file.txt ]; then
        allow_compressed=1;
    else 
        allow_compressed=0;
    fi

    if `command -v gzip >& /dev/null` && [ $allow_compressed == 1 ]; then
        test_header "${file} - with gzip"
        cp $file copy.$file
        gzip copy.$file
        ./genozip copy.${file}.gz -ft -o ${output}.genozip || exit 1
        rm copy.${file}.gz
    fi
    
    if `command -v bzip2 >& /dev/null` && [ $allow_compressed == 1 ]; then
        test_header "${file} - with bzip2"
        cp $file copy.$file
        bzip2 copy.$file
        ./genozip copy.${file}.bz2 -ft -o ${output}.genozip || exit 1
        rm copy.${file}.bz2
    fi
    
    if `command -v xz >& /dev/null` && [ $allow_compressed == 1 ]; then
        test_header "${file} - with xz"
        cp $file copy.$file
        xz copy.$file
        ./genozip copy.${file}.xz -ft -o ${output}.genozip || exit 1
        rm copy.${file}.xz
    fi
        
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

done

files=`ls backward-compatibility-test/*.vcf` 
for file in $files; do
    test_header "$file - backward compatability test"
    ./genounzip ${file}.genozip -fo $output || exit 1
    cmp_2_files $file $output
done

if `command -v gtshark >& /dev/null`; then
    test_header "test-file.vcf --gtshark"
    ./genozip test-file.vcf --gtshark -ft -o ${output}.genozip || exit 1
fi

test_header "test-file.vcf without FORMAT or samples"
cut -f1-8 test-file.vcf > test-input.vcf
./genozip test-input.vcf -ft -o ${output}.genozip || exit 1
rm test-input.vcf

if `command -v samtools >& /dev/null`; then
    test_header "test_file.sam - input and output as BAM"
    samtools view test-file.sam -OBAM -h > bam-test.input.bam    
    ./genozip bam-test.input.bam -fto ${output}.genozip || exit 1
    ./genounzip ${output}.genozip --force --output bam-test.output.bam
    cmp_2_files bam-test.input.bam bam-test.output.bam.fake-extension
    rm bam-test.input.bam bam-test.output.bam
fi

printf "\nALL GOOD!\n"

rm -f $output ${output}.genozip
