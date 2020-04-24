#!bash

output=regression-output

files=(test-file.vcf test-file.sam test-file.fq test-file.fa test-file.23andme)
for file in ${files[@]}; do
    echo " "
    echo TESTING $file
    ./genozip $file -ft -o ${output}.genozip || exit 1
done

echo " "
echo "Backward compatability tests"
files=`ls backward-compatibility-test/*.vcf` 
for file in $files; do
    ./genounzip ${file}.genozip -fo $output
    if [[ `md5sum $file | cut -d" " -f1` != `md5sum regression-output | cut -d" " -f1` ]]; then
        echo FAILED!
        exit 1
    fi
done
 
echo " "
echo "ALL GOOD!"

rm $output ${output}.genozip
