#!bash

set -e # exit when any command fails

echo " "
echo "VCF test file"
./genozip test-file.vcf -ft

echo " "
echo "SAM test file"
./genozip test-file.sam -ft

echo " "
echo "FASTQ test file"
./genozip test-file.fq -ft

echo " "
echo "FASTA test file"
./genozip test-file.fa -ft

echo " "
echo "23andMe test file"
./genozip test-file.23andme -ft

echo "ALL GOOD!"