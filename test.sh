#!/bin/bash 

cmp_2_files() {
    if (( `$md5 $1 ${2%.*} | cut -d" " -f1 | uniq | wc -l` != 1 )) ; then
        echo "MD5 comparison FAILED:"
        $md5 $1 ${2%.*}
        exit 1
    fi
}

test_header() {
    sep="=======================================================================================================\n"
    printf "\n${sep}TESTING ${FUNCNAME[1]} $1 \n${sep}"
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
    rm -f test/*.genozip
    test_header "genozip $*"

    local zip_args=( $1 )
    local unzip_args=( $2 )
    local args=( "$@" )
    local files=() # ${args[@]:2} )

    local copies=()
    if [[ ${zip_args[0]} == "COPY" ]]; then
        zip_args=( ${zip_args[@]:1} ) # remove COPY
        for file in ${files[@]}; do 
            cp test/$file test/copy.$file
            copies+=( test/copy.$file )
            args+=( copy.$file ) # adding the test/ prefix in a sec
        done
    fi

    if [[ ${zip_args[0]} == "NOPREFIX" ]]; then
        zip_args=( ${zip_args[@]:1} ) # remove NOPREFIX
        files=( ${args[@]:2} )
    else
        mapfile -t -d $'\0' files < <(printf "test/%s\0" "${args[@]:2}")
    fi

    local single_output=0
    if [[ ${zip_args[0]} == "CONCAT" ]]; then
        zip_args=( ${zip_args[@]:1} ) # remove CONCAT
        local single_output=1
    fi

    if (( ${#files[@]} == 1 || $single_output )); then # test with Adler32, unless caller specifies --md5
        $genozip $arg1 ${zip_args[@]} ${files[@]} -o $output -f || exit 1
        $genounzip $arg1 ${unzip_args[@]} $output -t || exit 1
    else
        $genozip $arg1 ${zip_args[@]} ${unzip_args[@]} ${files[@]} -ft || exit 1
    fi
    
    if [ ! -n "$single_output" ]; then
        local count=`ls -1 test/*.genozip | wc -l`
        local num_files=$(( $# - 1 ))
        if (( $count != $num_files )); then
            echo "Error: compressed $num_files files, but only $count genozip files found in td. Files compressed: "
            echo ${files[@]}
            exit 1
        fi
    fi

    rm -f test/*.genozip ${copies[@]}
}

test_redirected() { # $1=filename  $2...$N=optional extra genozip arg
    local file=test/$1
    test_header "$1 - redirected from stdin"
    local args=( "$@" )
    cat $file | $genozip $arg1 ${args[@]:1} --test --force --output $output --input-type ${file#*.} - || exit 1
    rm -f $output
}

test_eval() { # $1=filename $2=eval expression
    local file=test/$1
    eval `echo $2 | sed 's/FILE/$file/g'` > $1.eval
    test_header "$1 $2"
    $genozip $arg1 $1.eval -ft -o $output || exit 1
    rm -f $output $1.eval
}

test_unix_style() {  # $1=filename  ; optional $2=rm
    local file=test/$1
    test_header "$1 - basic test - Unix-style end-of-line"

    if [ ! -f $file ] ; then echo "$1: File $file not found"; exit 1; fi

    cat $file | tr -d "\r" > unix-nl.$1
    $genozip unix-nl.$1 -ft -o $output || exit 1    

    if [ -n "$3" ]; then rm -f unix-nl.$1 $output; fi
}

test_windows_style() {  # $1=filename 
    local file=test/$1
    test_header "$1 - Window-style end-of-line"

    if [ ! -f $file ] ; then echo "$1: File $file not found"; exit 1; fi

    if [ ! -n "$is_mac" ]; then  # note: sed on mac doesn't recognize \r
        sed 's/$/\r/g' unix-nl.$1 > windows-nl.$1 || exit 1 # note: sed on mac doesn't recognize \r
        $genozip $arg1 windows-nl.$1 -ft -o $output || exit 1
        rm -f unix-nl.$1 windows-nl.$1
    fi
}

test_stdout()
{
    if [ -z "$is_windows" ]; then # windows can't redirect binary data
        local file=test/$1
        test_header "$$1 - redirecting stdout"
        $genozip $arg1 ${file} --stdout > $output || exit 1
        $genounzip $arg1 $output -f || exit 1
        cmp_2_files $file $output
        rm -f $output
    fi
}

test_multi_bound()
{
    local file=test/$1
    test_header "$1 - bind & unbind (2 files with 2 components each)"
    local file1=test/copy1.$1
    local file2=test/copy2.$1
    cp -f $file $file1
    cat $file | sed 's/PRFX/FILE2/g' > $file2
    $genozip $arg1 $file1 $file2 -ft -o $output || exit 1 # test as bound
    local output2=test/output2.genozip
    cp -f $output $output2
    $genounzip $arg1 $output $output2 -u -t || exit 1 # test unbind 2x2
    rm -f $file1 $file2 $output $output2
}

test_optimize()
{
    local file=test/$1
    test_header "$1 --optimize - NOT checking correctness, just that it doesn't crash"
    $genozip $arg1 $file -f --optimize -o $output || exit 1
}

test_translate_bam_to_sam() # $1 bam file 
{
    local bam=test/$1
    local sam=${bam%.bam}.sam
    if [ ! -f $bam ] ; then echo "$bam: File not found"; exit 1; fi
    if [ ! -f $sam ] ; then echo "$sam: File not found"; exit 1; fi

    test_header "$1 - translate BAM to SAM"

    $genozip $arg1 -f $bam -o $output  || exit 1
    $genocat $output > copy.sam  || exit 1
    cmp_2_files $sam copy.sam
    rm -f copy.sam $output
}

test_translate_sam_to_bam() # $1 bam file 
{
    local bam=test/$1
    local sam=${bam%.bam}.sam
    if [ ! -f $bam ] ; then echo "$bam: File not found"; exit 1; fi
    if [ ! -f $sam ] ; then echo "$sam: File not found"; exit 1; fi

    test_header "$1 - translate SAM to BAM"

    $genozip $arg1 -f $sam -o $output  || exit 1
    $genounzip $output --bam -o copy.bam   || exit 1
    cmp_2_files $bam copy.bam
    rm -f copy.bam $output
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
        test_unix_style $file minimal rm
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
        test_standard "NOPREFIX" " " file://${path}test/$file
        test_standard "-p123" "--password 123" $file
        test_redirected $file
        test_stdout $file
        test_standard "COPY" " " $file
        test_multi_bound $file
        test_optimize $file
    done
}

# pre-compressed files (except BGZF) and non-precompressed BAM
batch_precompressed()
{
    local files=(basic-gzip.sam.gz basic-bz2.sam.bz2 basic-xz.sam.xz basic-nobgzip.bam) 
    local file
    for file in ${files[@]}; do
        test_standard " " " " $file
        test_standard "NOPREFIX" " " file://${path}test/$file
        test_standard "-p123" "--password 123" $file
    done
}

# bgzf files
batch_bgzf()
{
    local files=(basic.bam basic-bgzf-6.sam.gz basic-bgzf-9.sam.gz basic-bgzf-6-no-eof.sam.gz)
    local file
    for file in ${files[@]}; do
        test_standard " " " " $file
        test_standard "NOPREFIX" " " file://${path}$file
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
    local files=(basic-domqual.fq basic-domqual.sam basic-unaligned.sam)
    local file
    for file in ${files[@]}; do
        test_unix_style $file
        test_standard "-p123" "-p 123" $file
        test_standard "COPY" " " $file
        test_multi_bound $file
        test_optimize $file
    done
}

# Test SAM<->BAM translation
batch_translate_bam_sam()
{
    # note: we have these files in both sam and bam versions generated with samtools
    local files=(test/test.NA12878.chr22.1x.bam test/test.m64136_200621_234916.ccs.10k.bam test/test.GFX0241869.bam)
    local file
    for file in ${files[@]}; do
        test_translate_bam_to_sam $file
        test_translate_sam_to_bam $file
    done
}

batch_genocat_tests()
{
    # FASTA genocat tests
    test_count_genocat_lines basic.fa "--sequential" 9
    test_count_genocat_lines basic.fa "--header-only" 3
    test_count_genocat_lines basic.fa "--header-one" 3
    test_count_genocat_lines basic.fa "--no-header" 15
    test_count_genocat_lines basic.fa "--no-header --sequential" 6
    test_count_genocat_lines basic.fa "--grep cytochrome" 6
    test_count_genocat_lines basic.fa "--grep cytochrome --sequential " 2
    test_count_genocat_lines basic.fa "--grep cytochrome --sequential --no-header " 1

    # FASTQ genocat tests
    test_count_genocat_lines basic.fq "--header-only" `grep @ basic.fq | wc -l` 
    test_count_genocat_lines basic.fq "--header-one" `grep @ basic.fq | wc -l`
    test_count_genocat_lines basic.fq "--grep line5" 4
    test_count_genocat_lines basic.fq "--grep line5 --header-only" 1
}

batch_backward_compatability()
{
    local files=`ls backward-compatibility-test/*.genozip` 
    local file
    for file in $files; do
        test_backward_compatability $file
    done
}

batch_real_world_subsets()
{
    echo "subsets (~3 VBs) or real world files"
    test_standard "-m" " " test/test.*
}

batch_misc_cases()
{
    # Test binding SAM files with lots of contigs (no reference)
    echo "binding SAM files with lots of contigs (no reference)"
    test_multi_bound test.transfly-unsorted.sam

    # VCF gtshark test
    test_standard --gtshark " " basic.vcf

    # VCF without samples
    echo "basic.vcf without FORMAT or samples"
    test_eval basic.vcf "cut -f1-8 FILE"
}

batch_small_reference()
{
    # Making a reference
    echo "Making a reference"
    file=basic-ref.fa 
    $genozip $arg1 --make-reference $file --force -o copy.${file}.ref.genozip || exit 1

    ref="--reference copy.${file}.ref.genozip"
    REF="--REFERENCE copy.${file}.ref.genozip"

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

    rm -f copy.${file}.ref.genozip $output
}

batch_real_reference()
{
    echo "command line with mixed SAM and FASTQ files with --reference"
    echo "Note: '$GRCh38' needs to be up to date with the latest genozip format"
    test_standard "-me$GRCh38" "-e$GRCh38" test.transfly-unsorted.sam test.transfly.fq test.transfly-sorted.sam

    echo "multiple bound SAM with --REFERENCE" 
    test_standard "-mE$GRCh38" " " test.transfly-unsorted.sam test.transfly-sorted.sam
    
    echo "SAM with --reference and --password" 
    test_standard "-me$GRCh38 --password 123" "-p123 -e$GRCh38" test.transfly-unsorted.sam
    
    echo "SAM with --reference and --password" 
    test_standard "-E$GRCh38 --password 123" "-p123" test.transfly-unsorted.sam
    
    echo "paired FASTQ with --reference, --password"
    test_standard "CONCAT -e$GRCh38 -p 123 --pair" "-p123" test.divon-R1.100K.fq.bz2 test.divon-R2.100K.fq.bz2

    echo "4 paired FASTQ with --REFERENCE"
    test_standard "COPY CONCAT -E$GRCh38 -2" "-u" test.divon-R1.100K.fq.bz2 test.divon-R2.100K.fq.bz2

    echo "multiple bound VCF with --reference, --md5 using hg19, and unbind"
    test_standard "COPY CONCAT -me$hg19" "-u" test.divon-R1.100K.fq.bz2 test.divon-R2.100K.fq.bz2

    echo "multiple VCF with --REFERENCE using hg19" 
    test_standard "-mE$hg19" " " test.ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf test.GFX0241869.filtered.snp.vcf
}

output=test/output.genozip

is_windows=`uname|grep -i mingw`
is_mac=`uname|grep -i Darwin`

hg19=data/hg19.p13.plusMT.full_analysis_set.ref.genozip
GRCh38=data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip

arg1=$1

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

batch_minimal
batch_basic
batch_precompressed
batch_bgzf
batch_special_algs
batch_translate_bam_sam
batch_genocat_tests
batch_backward_compatability
batch_real_world_subsets
batch_misc_cases
batch_small_reference
batch_real_reference

printf "\nALL GOOD!\n"

