#!/bin/bash 

TESTDIR=test
OUTDIR=$TESTDIR/tmp

cleanup() { 
    rm -f $OUTDIR/* 
}

cmp_2_files() {
#    if (( `$md5 $1 ${2%.*} | cut -d" " -f1 | uniq | wc -l` != 1 )) ; then

    if (( `$md5 $1 $2 | cut -d" " -f1 | uniq | wc -l` != 1 )) ; then
        echo "MD5 comparison FAILED:"
#        $md5 $1 ${2%.*}
        $md5 $1 $2
        exit 1
    fi
}

test_header() {
    sep="=======================================================================================================\n"
    printf "\n${sep}TESTING ${FUNCNAME[2]}: "
    echo $1 | tr "\\\\" "/" # \ -> \\ coming from $path string on Windows
    printf "$sep"
}

test_count_genocat_lines() {
    local cmd="$genocat $output $2 $arg1"
    test_header "$cmd"
    $genozip $arg1 $1 -fo $output || exit 1
    local wc=`$cmd | wc -l`
    if (( $wc != $3 )); then
        echo "FAILED - expected $3 lines, but getting $wc"
        exit 1
    fi  
}

test_standard()  # $1 genozip args $2 genounzip args $3... filenames 
{
    local zip_args=( $1 )
    local unzip_args=( $2 )
    local args=( "$@" )

    local copies=()
    if [[ ${zip_args[0]} == "COPY" ]]; then
        zip_args=( ${zip_args[@]:1} ) # remove COPY
        for file in ${args[@]:2}; do 
            cp $TESTDIR/$file $OUTDIR/copy.$file
            copies+=( $OUTDIR/copy.$file )
            args+=( tmp/copy.$file ) # adding the test/ prefix in a sec
        done
    fi

    local files=( ${args[@]:2} )
    if [[ ${zip_args[0]} == "NOPREFIX" ]]; then
        zip_args=( ${zip_args[@]:1} ) # remove NOPREFIX
    else
        files=( "${files[@]/#/${TESTDIR}/}" )
    fi

    local single_output=0
    if [[ ${zip_args[0]} == "CONCAT" ]]; then
        zip_args=( ${zip_args[@]:1} ) # remove CONCAT
        local single_output=1
    fi

    test_header "genozip $arg1 ${zip_args[*]} ${files[*]}" # after COPY, NOPREFIX and CONCAT modifications occurred
    
    if (( ${#files[@]} == 1 || $single_output )); then # test with Adler32, unless caller specifies --md5
        $genozip $arg1 ${zip_args[@]} ${files[@]} -o $output -f || exit 1
        $genounzip $arg1 ${unzip_args[@]} $output -t || exit 1
    else 
        $genozip $arg1 ${zip_args[@]} ${unzip_args[@]} ${files[@]} -ft || exit 1
    fi
    
    if [ ! -n "$single_output" ]; then
        local count=`ls -1 $TESTDIR/*.genozip | wc -l`   # unfortunately, these go to TESTDIR not OUTDIR
        local num_files=$(( $# - 1 ))
        if (( $count != $num_files )); then
            echo "Error: compressed $num_files files, but only $count genozip files found in td. Files compressed: "
            echo ${files[@]}
            exit 1
        fi

        rm -f $TESTDIR/*.genozip
    fi

    local list=`ls $OUTDIR`
    rm -f ${list/test\/tmp\/basic-ref.ref.genozip}
}

test_redirected() { # $1=filename  $2...$N=optional extra genozip arg
    test_header "$1 - redirected from stdin"
    local file=$TESTDIR/$1
    local args=( "$@" )
    cat $file | $genozip $arg1 ${args[@]:1} --test --force --output $output --input ${file#*.} - || exit 1
    cleanup
}

test_unix_style() {  # $1=filename  ; optional $2=rm
    test_header "$1 - basic test - Unix-style end-of-line"
    local file=$TESTDIR/$1

    if [ ! -f $file ] ; then echo "$1: File $file not found"; exit 1; fi

    cat $file | tr -d "\r" > $OUTDIR/unix-nl.$1
    $genozip $OUTDIR/unix-nl.$1 -ft -o $output || exit 1    

    # no cleanup as we need the output for test_windows_style
}

test_windows_style() {  # $1=filename 
    if [ -n "$is_mac" ]; then return ; fi  # note: can't run this test on mac - sed on mac doesn't recognize \r

    test_header "$1 - Window-style end-of-line"
    local file=$TESTDIR/$1

    if [ ! -f $file ] ; then echo "$1: File $file not found"; exit 1; fi

    sed 's/$/\r/g' $OUTDIR/unix-nl.$1 > $OUTDIR/windows-nl.$1 || exit 1 # note: sed on mac doesn't recognize \r
    $genozip $arg1 $OUTDIR/windows-nl.$1 -ft -o $output || exit 1
    cleanup
}

test_stdout()
{
    test_header "$1 - redirecting stdout"
    local file=$TESTDIR/$1
    
    $genozip $arg1 ${file} -fo $output || exit 1
    $genounzip $arg1 $output --stdout | tr -d "\r" > $OUTDIR/unix-nl.$1 || exit 1

    cmp_2_files $file $OUTDIR/unix-nl.$1
    cleanup
}

test_multi_bound() # $1=filename $2=REPLACE (optional)
{
    test_header "$1 - bind & unbind (2 files with 2 components each)"
    local file=$TESTDIR/$1
    local file1=$OUTDIR/copy1.$1
    local file2=$OUTDIR/copy2.$1

    cp -f $file $file1
    if [[ $2 == "REPLACE" ]]; then
        cat $file | sed 's/PRFX/FILE2/g' > $file2
    else
        cp -f $file $file2
    fi

    $genozip $arg1 $file1 $file2 -ft -o $output || exit 1 # test as bound
    local output2=$OUTDIR/output2.genozip
    cp -f $output $output2
    $genounzip $arg1 $output $output2 -u -t || exit 1 # test unbind 2x2
    cleanup
}

test_optimize()
{
    test_header "$1 --optimize - NOT checking correctness, just that it doesn't crash"
    local file=$TESTDIR/$1
    $genozip $arg1 $file -f --optimize -o $output || exit 1
    cleanup
}

test_translate_bam_to_sam() # $1 bam file 
{
    test_header "$1 - translate BAM to SAM"

    local bam=$TESTDIR/$1
    local sam=${bam%.bam}.sam
    local copy=$OUTDIR/copy.sam
    if [ ! -f $bam ] ; then echo "$bam: File not found"; exit 1; fi
    if [ ! -f $sam ] ; then echo "$sam: File not found"; exit 1; fi

    $genozip $arg1 -f $bam -o $output  || exit 1
    $genocat $arg1 $output --no-PG -fo $copy  || exit 1
    cmp_2_files $sam $copy
    cleanup
}

test_translate_sam_to_bam() # $1 bam file 
{
    test_header "$1 - translate SAM to BAM"

    local bam=$TESTDIR/$1
    local sam=${bam%.bam}.sam
    local new_bam=$OUTDIR/copy_new.bam
    local new_sam=$OUTDIR/copy_new.sam
    if [ ! -f $bam ] ; then echo "$bam: File not found"; exit 1; fi
    if [ ! -f $sam ] ; then echo "$sam: File not found"; exit 1; fi

    $genozip $arg1 -f $sam -o $output  || exit 1
    $genounzip $arg1 $output --bam --no-PG -fo $new_bam || exit 1
    
    # we compare the BAMs on a textual basis as the original BAM might be produced with a different version
    samtools view --no-PG -h $new_bam > $new_sam || exit 1
    cmp_2_files $sam $new_sam
    cleanup
}

# note: only runs it to see that it doesn't crash, doesn't validate results
test_translate_sambam_to_fastq() # $1 sam or bam file 
{
    test_header "$1 - translate SAM/BAM to FASTQ"

    local sambam=$TESTDIR/$1
    local fastq=$OUTDIR/copy.fastq.gz
    if [ ! -f $sambam ] ; then echo "$sambam: File not found"; exit 1; fi

    $genozip $arg1 -f $sambam -o $output  || exit 1
    $genounzip $arg1 $output -fo $fastq   || exit 1

    cleanup
}

# note: only runs it to see that it doesn't crash, doesn't validate results
test_translate_23andMe_to_vcf() # $1 23andMe file
{
    test_header "$1 - translate 23andMe to VCF"

    local me23=$TESTDIR/$1
    local vcf=$OUTDIR/copy.vcf.gz
    if [ ! -f $sambam ] ; then echo "$sambam: File not found"; exit 1; fi

    $genozip $arg1 -f $me23 -o $output  || exit 1
    $genounzip $arg1 $output -fo $vcf -e $hg19  || exit 1

    cleanup
}

test_backward_compatability()
{
    test_header "$1 - backward compatability test"
    $genounzip $arg1 -t $1 || exit 1
}

# minimal files - expose edge cases where fields have only 1 instance
batch_minimal()
{
    local files=(minimal.vcf minimal.sam minimal.fq minimal.fa minimal.gvf minimal.genome_Full.me23.txt)
    local file
    for file in ${files[@]}; do
        test_standard " " " " $file
    done
}

# basic files
batch_basic()
{
    local files=(basic.vcf basic.sam basic.fq basic.fa basic.gvf basic.genome_Full.me23.txt)
    local file
    for file in ${files[@]}; do
        
        test_unix_style $file
        test_windows_style $file
        test_standard "NOPREFIX CONCAT" " " file://${path}${TESTDIR}/$file
        test_standard "-p123" "--password 123" $file
        test_redirected $file
        test_stdout $file
        test_standard "COPY" " " $file
        test_multi_bound $file REPLACE # REPLAEC to adjust the contig name for .fa as we can't have two contigs with the same name
        test_optimize $file
    done
}

# pre-compressed files (except BGZF) and non-precompressed BAM
batch_precompressed()
{
    local files=(basic-gzip.sam.gz basic-bz2.sam.bz2 basic-xz.sam.xz basic-nobgzip.bam) 
    local file
    for file in ${files[@]}; do

        if [ -x "$(command -v xz)" -o "${file##*.}" != xz ] ; then # skip .xz files if xz is not installed
            test_standard " " " " "$file"
            test_standard "NOPREFIX CONCAT" " " file://${path}${TESTDIR}/$file
            test_standard "-p123" "--password 123" $file
        fi
    done
}

# bgzf files
batch_bgzf()
{
    local files=(basic.bam basic-bgzf-6.sam.gz basic-bgzf-9.sam.gz basic-bgzf-6-no-eof.sam.gz)
    local file
    for file in ${files[@]}; do
        test_standard " " " " $file
        test_standard "NOPREFIX CONCAT" " " file://${path}${TESTDIR}/$file
        test_standard "-p123" "--password 123" $file
        if [ -z "$is_windows" ]; then # windows can't redirect binary data
            test_redirected $file
        fi
        test_standard "COPY" " " $file
        test_multi_bound $file
    done
}
        
# files represent cases that cannot be written into the test files because they would conflict
batch_special_algs()
{
    local files=(basic-domqual.fq basic-domqual.sam basic-unaligned.sam basic-no-samples.vcf)
    local file
    for file in ${files[@]}; do
        test_unix_style $file                # standard
        test_standard "-p123" "-p 123" $file # encrypted
        test_standard "COPY" " " $file       # multiple files unbound
        test_multi_bound $file               # multiple files bound
        test_optimize $file                  # optimize - only compress to see that it doesn't error
    done
}

# Test translations
batch_translations()
{
    # note: we have these files in both sam and bam versions generated with samtools
    local files=(test.NA12878.chr22.1x.bam 
                 test.m64136_200621_234916.ccs.10k.bam  # unaligned SAM/BAM with no SQ records
                 test.human2.bam)
    local file
    for file in ${files[@]}; do
        test_translate_bam_to_sam $file

        if command -v samtools &> /dev/null; then # test_translate_sam_to_bam requires samtools
            test_translate_sam_to_bam $file
        fi

        test_translate_sambam_to_fastq $file
        test_translate_sambam_to_fastq ${file%.bam}.sam
    done

    test_translate_23andMe_to_vcf test.genome_Full.txt
}

batch_genocat_tests()
{
    # FASTA genocat tests
    local file=$TESTDIR/basic.fa
    test_count_genocat_lines $file "--sequential" 9
    test_count_genocat_lines $file "--header-only" 3
    test_count_genocat_lines $file "--header-one" 18
    test_count_genocat_lines $file "--no-header" 15
    test_count_genocat_lines $file "--no-header --sequential" 6
    test_count_genocat_lines $file "--grep cytochrome" 6
    test_count_genocat_lines $file "--grep cytochrome --sequential " 2
    test_count_genocat_lines $file "--grep cytochrome --sequential --no-header " 1

    # Multi-FASTA to Phylip
    local file=$TESTDIR/basic-multifasta.fa
    test_count_genocat_lines $file "--phylip" 4

    # FASTQ genocat tests
    file=$TESTDIR/basic.fq
    test_count_genocat_lines $file "--header-only" `grep @ $file | wc -l` 
    test_count_genocat_lines $file "--grep line5" 4
    test_count_genocat_lines $file "--grep line5 --header-only" 1
}

batch_backward_compatability()
{
    local files=`ls $TESTDIR/back-compat/*.genozip` 
    local file
    for file in $files; do
        test_backward_compatability $file
    done
}

batch_real_world_subsets()
{
    rm -f $TESTDIR/*.genozip # unfortunately, these go to TESTDIR not OUTDIR

    if [ -x "$(command -v xz)" ] ; then # xz available
        local files=( `cd test; ls -1 test.*vcf* test.*sam* test.*bam* test.*fq* test.*fastq* test.*fa* test.*fasta* test.*vcf* test.*gvf* test.*txt*` )
    else
        local files=( `cd test; ls -1 test.*vcf* test.*sam* test.*bam* test.*fq* test.*fastq* test.*fa* test.*fasta* test.*vcf* test.*gvf* test.*txt*|grep -v xz` )
    fi
    
    echo "subsets (~3 VBs) or real world files"
    test_standard "-m" " " ${files[@]}
}

batch_misc_cases()
{
    # Test binding SAM files with lots of contigs (no reference)
    echo "binding SAM files with lots of contigs (no reference)"
    test_multi_bound test.human-unsorted.sam
}

batch_external_tools()
{
    # VCF gtshark test
    if `command -v gtshark >& /dev/null`; then
        test_standard --gtshark " " basic.vcf
    fi

    # CRAM hg19
    if `command -v samtools >& /dev/null`; then
        echo "CRAM" 
        test_standard "-E$hg19" " " test.human2.cram   
    fi

    # BCF
    if `command -v bcftools >& /dev/null`; then
        echo "BCF"
        test_standard " " " " test.human2.filtered.snp.bcf    
    fi

    # unzip
    if `command -v unzip >& /dev/null`; then
        echo "unzip"
        test_standard " " " " test.genome_Full.zip    
    fi
}

batch_make_reference()
{
    cleanup

    # Making a reference
    echo "Making a reference"
    local fa_file=$TESTDIR/basic-ref.fa 
    local ref_file=$OUTDIR/basic-ref.ref.genozip
    $genozip $arg1 --make-reference $fa_file --force -o $ref_file || exit 1

    local ref="--reference $ref_file"
    local REF="--REFERENCE $ref_file"

    echo "unaligned SAM with --reference"
    test_standard "$ref" "$ref" basic-unaligned.sam 

    echo "unaligned SAM with --REFERENCE"
    test_standard "$REF" " " basic-unaligned.sam

    echo "unaligned BAM with --reference"
    test_standard "$ref" "$ref" basic-unaligned.bam

    echo "FASTQ with --REFERENCE"
    test_standard "$REF" " " basic.fq 

    echo "unaligned SAM with --REFERENCE - from stdin"
    test_redirected basic-unaligned.sam "$REF"

    cleanup
}

batch_reference()
{
    echo "paired FASTQ with --reference, --password"
    test_standard "CONCAT -e$GRCh38 -p 123 --pair" "-p123" test.human2-R1.100K.fq.bz2 test.human2-R2.100K.fq.bz2

    echo "4 paired FASTQ with --REFERENCE (not bound)"
    test_standard "COPY -E$GRCh38 --pair" " " test.human2-R1.100K.fq.bz2 test.human2-R2.100K.fq.bz2

    echo "4 paired FASTQ with --REFERENCE (bound)"
    test_standard "COPY CONCAT -E$GRCh38 -2" "-u" test.human2-R1.100K.fq.bz2 test.human2-R2.100K.fq.bz2

    echo "command line with mixed SAM and FASTQ files with --reference"
    echo "Note: '$GRCh38' needs to be up to date with the latest genozip format"
    test_standard "-me$GRCh38" "-e$GRCh38" test.human-unsorted.sam test.human.fq test.human-sorted.sam

    echo "multiple bound SAM with --REFERENCE" 
    test_standard "-mE$GRCh38" " " test.human-unsorted.sam test.human-sorted.sam
    
    echo "SAM with --reference and --password" 
    test_standard "-me$GRCh38 --password 123" "-p123 -e$GRCh38" test.human-unsorted.sam
    
    echo "SAM with --REFERENCE and --password" 
    test_standard "-E$GRCh38 --password 123" "-p123" test.human-unsorted.sam
    
    echo "multiple bound VCF with --reference, --md5 using hg19, and unbind"
    test_standard "COPY CONCAT -me$hg19" "-u" test.human2-R1.100K.fq.bz2 test.human2-R2.100K.fq.bz2

    echo "multiple VCF with --REFERENCE using hg19" 
    test_standard "-mE$hg19" " " test.ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf test.human2.filtered.snp.vcf
}

output=${OUTDIR}/output.genozip

is_windows=`uname|grep -i mingw`
is_mac=`uname|grep -i Darwin`

hg19=data/hs37d5.ref.genozip
GRCh38=data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip

if (( $# < 1 )); then
    echo "Usage: tesh.sh <start-test> [optional-genozip-arg]"
    exit 0
fi

arg1=$2 # optional arg for genozip/genounzip

# debug
is_debug=`echo $1|grep debug`
if [ -n "$is_debug" ]; then 
    debug=-debug; 
    shift
fi

# -----------------
# platform settings
# -----------------
if [ -n "$is_windows" ]; then
    genozip=./genozip${debug}.exe
    genounzip=./genounzip${debug}.exe
    genocat=./genocat${debug}.exe
    genols=./genols${debug}.exe
    path=`pwd| cut -c3-|tr / '\\\\'`\\
else
    genozip=./genozip${debug}
    genounzip=./genounzip${debug}
    genocat=./genocat${debug}
    genols=./genols${debug}
    path=$PWD/
fi

exes=($genozip $genounzip $genocat $genols)
for exe in ${exes[@]}; do
    if [ ! -x $exe ]; then
        echo "Error: $exe does not exist"
        exit 1
    fi
done

if `command -v md5 >& /dev/null`; then
    md5=md5 # mac
else
    md5=md5sum 
fi

mkdir $OUTDIR >& /dev/null
cleanup

start=$1
if (( $start <=  1 )); then batch_minimal                  ; fi
if (( $start <=  2 )); then batch_basic                    ; fi
if (( $start <=  3 )); then batch_precompressed            ; fi
if (( $start <=  4 )); then batch_bgzf                     ; fi
if (( $start <=  5 )); then batch_special_algs             ; fi
if (( $start <=  6 )); then batch_translations             ; fi
if (( $start <=  7 )); then batch_genocat_tests            ; fi
if (( $start <=  8 )); then batch_backward_compatability   ; fi
if (( $start <=  9 )); then batch_real_world_subsets       ; fi
if (( $start <= 10 )); then batch_misc_cases               ; fi
if (( $start <= 11 )); then batch_external_tools           ; fi
if (( $start <= 12 )); then batch_reference                ; fi
if (( $start <= 13 )); then batch_make_reference           ; fi

printf "\nALL GOOD!\n"

