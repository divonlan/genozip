#!/bin/bash 

# ------------------------------------------------------------------
#   test.sh
#   Copyright (C) 2019-2021 Black Paw Ventures Limited
#   Please see terms and conditions in the file LICENSE.txt

shopt -s extglob  # Enables extglob - see https://mywiki.wooledge.org/glob
export GENOZIP_TEST="Yes" # Causes output of debugger arguments
unset GENOZIP_REFERENCE   # initialize

ulimit -c unlimited # enable core dumps

TESTDIR=private/test
OUTDIR=$TESTDIR/tmp
REFDIR=../genozip/data

cleanup() { 
    rm -fR $OUTDIR/* $TESTDIR/*.bad $TESTDIR/*.rejects.* 
}

cmp_2_files() {
    if [ ! -f $1 ] ; then echo "File $1 not found while in cmp_2_files()"; exit 1; fi
    if [ ! -f $2 ] ; then echo "File $2 not found while in cmp_2_files()"; exit 1; fi

    if (( `$md5 $1 $2 | cut -d" " -f1 | uniq | wc -l` != 1 )) ; then
        echo "MD5 comparison FAILED: $1 $2"
        $md5 $1 $2
        exit 1
    fi
}

test_header() {
    sep="=======================================================================================================\n"
    printf "\n${sep}TESTING ${FUNCNAME[2]} (batch_id=${batch_id}): "
    echo $1 | tr "\\\\" "/" # \ -> \\ coming from $path string on Windows
    printf "$sep"
}

test_count_genocat_lines() { # $1 - genozip arguments $2 - genocat arguments $3 - expected number of output lines
    test_header "genozip $1 ; genocat $2"
    $genozip $1 -fo $output || exit 1
    $genocat $output $2 -fo $recon || exit 1
    local wc=`cat $recon |wc -l`

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
            args+=( tmp/copy.$file ) # adding the ${TESTDIR}/ prefix in a sec
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

    test_header "$genozip ${zip_args[*]} ${files[*]}" # after COPY, NOPREFIX and CONCAT modifications occurred
    
    if (( ${#files[@]} == 1 || $single_output )); then # test with Adler32, unless caller specifies --md5
        $genozip ${zip_args[@]} ${files[@]} -o $output -f || exit 1
        $genounzip ${unzip_args[@]} $output -t || exit 1
    else 
        $genozip ${zip_args[@]} ${unzip_args[@]} ${files[@]} -ft || exit 1
    fi
    
    if [ ! -n "$single_output" ]; then
        local count=`ls -1 $TESTDIR/*.genozip | wc -l`   # unfortunately, these go to TESTDIR not OUTDIR
        local num_files=$(( $# - 1 ))
        if (( $count != $num_files )); then
            echo "Error: compressed $num_files files, but only $count genozip files found. Files compressed: "
            echo ${files[@]}
            exit 1
        fi
        cleanup
    fi
}

test_redirected() { # $1=filename  $2...$N=optional extra genozip arg
    test_header "$1 - redirecting $file to genozip via stdin"
    local file=$TESTDIR/$1
    local args=( "$@" )

    # instead of passing the input, we pass the filename. that works, because the code compares extensions
    local input=${file#*.}
    if [[ $input == genome_Full.me23.txt ]]; then input=23andme; fi

    cat $file | $genozip ${args[@]:1} --test --force --output $output --input $input - || exit 1
    cleanup
}

test_unix_style() {  # $1=filename  ; optional $2=rm
    test_header "$1 - basic test - Unix-style end-of-line"
    local file=$TESTDIR/$1

    if [ ! -f $file ] ; then echo "$1: File $file not found"; exit 1; fi

    cat $file | tr -d "\r" > $OUTDIR/unix-nl.$1 || exit 1
    $genozip $OUTDIR/unix-nl.$1 -ft -o $output || exit 1    

    # no cleanup as we need the output for test_windows_style
}

test_windows_style() {  # $1=filename 
    if [ -n "$is_mac" ]; then return ; fi  # note: can't run this test on mac - sed on mac doesn't recognize \r

    test_header "$1 - Window-style end-of-line"
    local file=$TESTDIR/$1

    if [ ! -f $file ] ; then echo "$1: File $file not found"; exit 1; fi

    sed 's/$/\r/g' $OUTDIR/unix-nl.$1 > $OUTDIR/windows-nl.$1 || exit 1 # note: sed on mac doesn't recognize \r
    $genozip $OUTDIR/windows-nl.$1 -ft -o $output || exit 1
    cleanup
}

test_stdout()
{
    test_header "$1 - redirecting stdout"
    local file=$TESTDIR/$1
    local ext=${file#*.}
    local arg;

    if [ "$ext" = "bam" ]; then arg='--bam -z0'; fi

    $genozip ${file} -fo $output || exit 1
    ($genocat --no-pg $output $arg || exit 1) | tr -d "\r" > $OUTDIR/unix-nl.$1 

    cmp_2_files $file $OUTDIR/unix-nl.$1
    cleanup
}

test_optimize()
{
    test_header "$1 --optimize - NOT checking correctness, just that it doesn't crash"
    local file=$TESTDIR/$1
    $genozip $file -f --optimize -o $output || exit 1
    cleanup
}

test_md5()
{
    test_header "$1 --md5 - see that it is the correct MD5"
    local file=$TESTDIR/$1

    $genozip $file -f --md5 -o $output || exit 1
    genozip_md5=`$genols $output | grep $output | cut -c 51-82`
    real_md5=`$md5 $file | cut -d" " -f1`

    if [[ "$genozip_md5" != "$real_md5" ]]; then echo "FAILED - expected $file to have MD5=\"$real_md5\" but genozip calculated MD5=\"$genozip_md5\""; exit 1; fi

    cleanup
}

test_translate_sam_to_bam_to_sam() # $1 bam file $2 genozip options $3 genocat options
{
    test_header "$1 - translate SAM to BAM to SAM \"$2\" \"$3\""

    local bam=$TESTDIR/$1
    local sam=${bam%.bam}.sam
    local new_bam=$OUTDIR/copy_new.bam
    local new_sam=$OUTDIR/copy_new.sam
    if [ ! -f $bam ] ; then echo "$bam: File not found"; exit 1; fi
    if [ ! -f $sam ] ; then echo "$sam: File not found"; exit 1; fi

    # SAM -> BAM
    echo "STEP 1: sam -> sam.genozip"
    $genozip -f $sam -o $output $2 || exit 1

    echo "STEP 2: sam.genozip -> bam"
    $genocat $output --bam --no-PG -fo $new_bam $3 || exit 1

    # BAM -> SAM
    echo "STEP 3: bam -> bam.genozip"
    $genozip -f $new_bam -o $output $2 || exit 1

    echo "STEP 4: bam.genozip -> sam"
    $genocat $output --sam --no-PG -fo $new_sam $3 || exit 1

    # compare original SAM and SAM created via sam->sam.genozip->bam->bam.genozip->sam
    echo "STEP 5: compare original and output SAMs"
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

    $genozip -f $sambam -o $output || exit 1
    $genocat $output -fo $fastq    || exit 1

    cleanup
}

view_file()
{
    if [[ $1 =~ \.gz$ ]]; then  
        gunzip -c $1 || exit 1
    else
        cat $1
    fi 
}

batch_print_header()
{
    batch_id=$((batch_id + 1))
    echo "***************************************************************************"
    echo "******* ${FUNCNAME[1]} (batch_id=${batch_id}) " $1
    echo "***************************************************************************"
}

# minimal files - expose edge cases where fields have only 1 instance
batch_minimal()
{
    batch_print_header
    local files=(minimal.vcf minimal.sam minimal.fq minimal.fa minimal.gvf minimal.genome_Full.me23.txt minimal.kraken)
    local file
    for file in ${files[@]}; do
        test_standard " " " " $file
    done
}

# basic files
batch_basic()
{
    batch_print_header

    local file replace
    file=$1

    # if [ $file == basic.chain ]; then
    #     export GENOZIP_REFERENCE=${hs37d5}:${GRCh38}
    # else
    #     unset GENOZIP_REFERENCE
    # fi

    test_standard "" "" $file

    test_md5 $file # note: basic.bam needs to be non-BGZF for this to pass

    if [ $file != basic.bam ] && [ $file != basic.generic ]; then # binary files have no \n 
        test_unix_style $file
        test_windows_style $file
        replace=REPLACE
    else
        replace=
    fi

    test_standard "NOPREFIX CONCAT $ref" " " file://${path}${TESTDIR}/$file
    test_standard "-p123 $ref" "--password 123" $file
    if [ -z "$is_windows" ] || [ $file != basic.bam ]; then # can't redirect binary files in Windows
    # if [ -z "$is_windows" ]; then # in windows, we don't support redirecting stdin (bug 339)
        test_redirected $file
    # fi
#    if [ -z "$is_windows" ] || [ $file != basic.bam ]; then # can't redirect binary files in Windows
        test_stdout $file
    fi
    test_standard "COPY $ref" " " $file

    test_optimize $file
    unset GENOZIP_REFERENCE
}

# pre-compressed files (except BGZF) and non-precompressed BAM
batch_precompressed()
{
    batch_print_header
    local files=(basic-gzip.sam.gz basic-bz2.sam.bz2 basic-xz.sam.xz) 
    local file
    for file in ${files[@]}; do

        if [ -x "$(command -v xz)" -o "${file##*.}" != xz ] ; then # skip .xz files if xz is not installed
            test_standard " " " " "$file"
            test_standard "NOPREFIX CONCAT" " " file://${path}${TESTDIR}/$file
            test_standard "-p123" "--password 123" $file
        fi
    done
}

verify_bgzf() # $1 file that we wish to inspect $2 expected result (1 bgzf 0 not-bgzf)
{
    if [ "$(head -c4 $1 | od -x | head -1 | awk '{$2=$2};1')" == "0000000 8b1f 0408" ]; then 
        if [ $2 -eq 0 ]; then
            echo $1 is unexpectedly BGZF-compressed
            exit 1
        fi
    else
        if [ $2 -eq 1 ]; then
            echo $1 is not BGZF-compressed
            exit 1
        fi
    fi
}

# bgzf files
batch_bgzf()
{
    batch_print_header
    local files=(basic-bgzf.bam basic-bgzf-6.sam.gz basic-bgzf-9.sam.gz basic-bgzf-6-no-eof.sam.gz basic-1bgzp_block.bam)
    #local files=()
    local file
    for file in ${files[@]}; do
        test_standard " " " " $file
        test_standard "NOPREFIX CONCAT" " " file://${path}${TESTDIR}/$file
        test_standard "-p123" "--password 123" $file
#        if [ -z "$is_windows" ]; then # in windows, we don't support redirecting stdin
            test_redirected $file
#        fi
        test_standard "COPY" " " $file
    done

    test_header "sam -> sam.genozip -> genocat to sam.gz - see that it is BGZF"
    local sam_gz=${OUTDIR}/bgzf_test.sam.gz
    $genozip ${TESTDIR}/basic.sam -fo $output || exit 1
    $genocat --no-pg $output -fo $sam_gz || exit 1
#    if [ "$(head -c4 $sam_gz | od -x | head -1 | awk '{$2=$2};1')" != "0000000 8b1f 0408" ] ; then  # the awk converts the Mac output to be the same as Linux (removing redundant spaces)
    verify_bgzf $sam_gz 1

    test_header "sam-> sam.genozip -> genocat to bam - see that it is BGZF"
    local bam=${OUTDIR}/bgzf_test.bam
    $genocat --no-pg $output -fo $bam || exit 1
    verify_bgzf $bam 1

    test_header "sam.gz -> sam.genozip -> genocat to sam - see that it is not BGZF"
    local sam=${OUTDIR}/bgzf_test.sam
    $genozip $sam_gz -fo $output || exit 1
    $genocat --no-pg $output -fo $sam || exit 1
    verify_bgzf $sam 0

    test_header "sam.gz -> sam.genozip -> genounzip to sam.gz - see that it is BGZF"
    local sam_gz2=${OUTDIR}/bgzf_test2.sam.gz
    $genozip $sam_gz -fo $output || exit 1
    $genounzip $output -fo $sam_gz2 || exit 1
    verify_bgzf $sam_gz2 1

    test_header "bam -> bam.genozip -> genounzip to bam - see that it is BGZF"
    local bam2=${OUTDIR}/bgzf_test2.bam
    $genozip $bam -fo $output || exit 1
    $genounzip $output -fo $bam2 || exit 1
    verify_bgzf $bam2 1

    test_header "bam -> bam.genozip -> genounzip -z0 to bam - see that it is not BGZF"
    $genounzip $output -z0 -fo $bam2 || exit 1
    verify_bgzf $bam2 0

    # test with gencomp
    file=special.sag-by-sa.bam.gz
    $genozip ${TESTDIR}/$file -fo $output --force-gencomp
    $genounzip $output -fo ${OUTDIR}/$file
    verify_bgzf ${OUTDIR}/$file 1
}

batch_subdirs()
{
    batch_print_header
    local files=(${TESTDIR}/minimal.sam ${TESTDIR}/basic-subdirs)

    cleanup
    $genozip -Dft ${files[*]} || exit 1

    # verify that the file residing in the subdir was compressed
    $genounzip -f ${TESTDIR}/basic-subdirs/basic.subdirs.txt.genozip || exit
    cleanup
}   
        
# files represent cases that cannot be written into the test files because they would conflict
batch_special_algs()
{
    batch_print_header
    local files=(basic-domqual.fq basic-domqual.sam basic-unaligned.sam basic-no-samples.vcf)
    local file
    for file in ${files[@]}; do
        test_unix_style $file                # standard
        test_standard "-p123" "-p 123" $file # encrypted
        test_standard "COPY" " " $file       # multiple files unbound
        test_optimize $file                  # optimize - only compress to see that it doesn't error
    done
}

batch_dvcf()
{
    batch_print_header

    local files=(minimal.vcf basic-dvcf-source.vcf basic-dvcf-luft.vcf test.NA12878.sorted.vcf test.clinvar37.vcf.gz test.1KG-37.indels.vcf test.chr17.SS6004478.vcf test.ExAC.vcf.gz)
    local file

    # prepare chain file
    test_header "${files[0]} - DVCF test - preparing chain file"
    $genozip -e $hs37d5 -e $GRCh38 ${chain37_38%%.genozip} -fqo $chain --match-chrom || exit 1

    # test explicit reference
    test_header "${files[0]} - DVCF test - explicit reference"
    $genozip ${TESTDIR}/${files[0]} -fo $output -C $chain -e $hs37d5 -e $GRCh38 || exit 1

    for file in ${files[@]}; do
        test_header "$file - DVCF test"

        local src=$OUTDIR/${file}
        local primary=$OUTDIR/primary.vcf
        local luft=$OUTDIR/luft.vcf
        local primary2=$OUTDIR/primary2.vcf
        local dvcf=${primary%%.vcf}.d.vcf.genozip
        
        # compare src to primary (ignoring header and INFO fields)
        echo -n "Step 1: make $dvcf from $file : " 
        if [ "$file" == basic-dvcf-luft.vcf ]; then # this file is already a Luft rendition
            $genozip ${TESTDIR}/$file -fo $dvcf || exit 1
        else
            $genozip -C $chain ${TESTDIR}/$file -fo $dvcf --dvcf-rename="FORMAT/QDF:STRAND>QDR,QDR:STRAND>QDF" --dvcf-drop="INFO/CLN:REFALT" || exit 1
        fi
    
        echo -n "Step 2: make ${primary} from $dvcf : " 
        $genocat $dvcf --no-pg -fo ${primary} || exit 1

        # convert primary -> luft -> primary
        echo -n "Step 3: make ${luft} from $dvcf : " 
        $genocat --luft --no-pg $dvcf -fo ${luft} || exit 1
        echo -n "Step 4: make ${luft}.genozip from ${luft} : " 
        $genozip $luft -fo ${luft}.genozip || exit 1
        echo -n "Step 5: make ${primary2} from ${luft}.genozip : " 
        $genocat ${luft}.genozip --no-pg -fo ${primary2} || exit 1
        echo "Step 6: compare $primary to $primary2" 
        cmp_2_files $primary $primary2 

        rm -f ${src}.noinfo $dvcf ${primary}.noinfo ${luft} ${luft}.genozip ${primary2}
    done
    cleanup
}

batch_match_chrom()
{
    batch_print_header
    
    # hg19 tests
    local files=(basic-dvcf-source.vcf basic.genome_Full.me23.txt special.match.sam special.match.bam test.homo_sapiens_incl_consequences-chrY.gvf)
    local file f
    for f in ${files[@]}; do

        test_header "$f --match-chrom"

        file=$TESTDIR/$f # has contigs from both styles
        one=$OUTDIR/one.$f
        two=$OUTDIR/two.$f
        three=$OUTDIR/three.$f

        # convert to CHROM_STYLE_22
        $genozip --match-chrom $file -fo ${one}.genozip -e $hg19 || exit 1
        $genounzip -z0 ${one}.genozip || exit 1

        #convert CHROM_STYLE_chr22 and then to CHROM_STYLE_22
        $genozip --match-chrom $file -fo ${two}.genozip -e $hs37d5 || exit 1
        $genounzip -z0 -f ${two}.genozip || exit 1

        $genozip --match-chrom $two -fo ${three}.genozip -e $hg19 || exit 1
        $genounzip -z0 -f ${three}.genozip || exit 1

        cmp_2_files $three $one 
    done

    cleanup
}

test_kraken() { # $1 file ; $2 1st genocat arguments ; $3 2nd genocat arguments
    test_header "$file - kraken test (genozip $1 ; genocat $2 ; genocat $3)"

    local zip_args=$1
    local cat1_args=$2
    local cat2_args=$3

    $genozip $1 -fo $output || exit 1

    local lines_plus=`$genocat_no_echo $output -Hq ${cat1_args[*]} --count`
    if [ "$lines_plus" == "" ] || [ "$lines_plus" -eq 0 ]; then echo "genocat error - \$lines_plus=\"$lines_plus\""; exit 1; fi

    local lines_minus=`$genocat_no_echo $output -Hq $3 --count`
    if [ "$lines_minus" == "" ] || [ "$lines_minus" -eq 0 ]; then echo "genocat error - \$lines_minus=\"$lines_minus\""; exit 1; fi
    
    local lines=`$genocat_no_echo -Hq $output --count`
    if [ "$lines" == "" ] || [ "$lines" -eq 0 ]; then echo "genocat error - \$lines=\"$lines\""; exit 1; fi
    
    echo "$file : lines_plus=$lines_plus lines_minus=$lines_minus lines=$lines"
    if [ $(($lines_plus + $lines_minus)) -ne $lines ]; then
        echo "$file: adding up kraken positive- and negative- filtered data, isn't the size of the original file"
        exit 1
    fi        
}

batch_kraken() # $1 genozip arguments #2 genocat (one of them must include --kraken)
{
    batch_print_header

    local files=(${TESTDIR}/basic.bam ${TESTDIR}/basic-test-kraken.fq ${TESTDIR}/basic-test-kraken.fa ${TESTDIR}/basic.kraken \
                 ${TESTDIR}/basic.sam)
    local file

    $genozip ${TESTDIR}/basic.kraken -fo $kraken

    # testing filtering FASTA, FASTQ, SAM and KRAKEN itself with --taxid 
    for file in ${files[@]}; do
        test_kraken "$file $1" "$2 -k570" "$2 -k^570"
    done

    # testing filtering SAM translated to FASTQ with --taxid 
    test_kraken "${TESTDIR}/basic.bam $1" \
                "--fastq -k570 $2" \
                "--fastq -k^570 $2"
    
    # testing filtering SAM translated to FASTQ with --taxid and +0 
    test_kraken "${TESTDIR}/basic.bam $1" \
                "--fastq -k570+0 $2" \
                "--fastq -k^570+0 $2"
    
    # testing filtering SAM with a concatenated KRAKEN file (representing slightly different classifications
    # originating from separate kraken2 of R1 and R2 FASTQ files)
    cat ${TESTDIR}/basic.kraken ${TESTDIR}/basic-2nd-file.kraken | $genozip -i kraken -fo $kraken

    test_kraken "${TESTDIR}/basic.bam $1"\
                "-k570 $2" \
                "-k^570 $2"

    test_kraken "${TESTDIR}/basic.bam $1"\
                "-k570+0 $2" \
                "-k^570+0 $2"
    cleanup
}

# unit test for ref_copy_compressed_sections_from_reference_file. 
batch_copy_ref_section()
{
    batch_print_header

    #created with -r1:9660000-10650000, and contains 99% of vb=11 of hs37d5.ref.genozip which is 3867649-4834572
    local file=${TESTDIR}/unit-test.-E.copy-ref-section.sam.gz

    $genozip -E $hs37d5 -p 123 -ft $file
}    

# test -@1 - different code paths
batch_single_thread()
{
    batch_print_header

    # note -@1 will override previous -@
    test_standard "-@1" "-@1" basic.vcf 
    
    # with reference
    test_standard "-@1 -e$hs37d5 -e$GRCh38" "-@1" basic.chain 
    
    # with chain and reference (note: cannot --test dual-coord files)
    $genozip -@1 -C $chain37_38 ${TESTDIR}/basic-dvcf-source.vcf -f || exit 1

    # with kraken aux file
    $genozip -@1 ${TESTDIR}/basic.kraken -fo $kraken || exit 1
    test_kraken "-@1 -K$kraken ${TESTDIR}/basic.bam" "-@1 -k570" "-@1 -k^570" 
}

batch_iupac() 
{
    batch_print_header

    # SAM - genocat and verifying with wc
    non_iupac_lines=$(( `grep -v "^@" ${TESTDIR}/basic.sam | wc -l` - 1 )) # we have 1 IUPAC line (E100020409L1C001R0030000234)

    test_header "genocat --bases ACGTN (SAM)"
    test_count_genocat_lines ${TESTDIR}/basic.sam "-H --bases=ACGTN" $non_iupac_lines

    test_header "genocat --bases ^ACGTN (SAM)"
    test_count_genocat_lines ${TESTDIR}/basic.sam "-H --bases=^ACGTN" 1

    # SAM - using --count
    test_header "genocat --bases ACGTN --count (SAM)"
    local count=`$genocat_no_echo $output -H --bases ACGTN --count -q`
    if [ "$count" == "" ]; then echo genocat error; exit 1; fi

    if [ "$count" -ne $non_iupac_lines ]; then echo "bad count = $count, expecting $non_iupac_lines"; exit 1; fi

    test_header "genocat --bases ^ACGTN --count (SAM)"
    local count=`$genocat_no_echo $output -H --bases ^ACGTN --count -q`
    if [ "$count" -ne 1 ]; then echo "bad count = $count"; exit 1; fi

    # BAM - using --count
    test_header "genocat --bases ACGTN --count --bam"
    local count=`$genocat_no_echo $output --bam --bases ACGTN --count -q`
    if [ "$count" == "" ]; then echo genocat error; exit 1; fi

    if [ "$count" -ne $non_iupac_lines ]; then echo "bad count = $count, expecting $non_iupac_lines"; exit 1; fi

    test_header "genocat --bases ^ACGTN --count --bam"
    local count=`$genocat_no_echo $output --bam --bases ^ACGTN --count -q`
    if [ "$count" -ne 1 ]; then echo "bad count = $count"; exit 1; fi

    # FASTQ (verifying with wc)
    test_count_genocat_lines ${TESTDIR}/basic.fq "-H --IUPAC=ACGTN" 20
    test_count_genocat_lines ${TESTDIR}/basic.fq "-H --IUPAC=^ACGTN" 4
    non_iupac_lines=5
    
    # FASTQ - using --count
    test_header "genocat --bases ACGTN --count (FASTQ)"
    local count=`$genocat_no_echo $output -H --bases ACGTN --count -q`
    if [ "$count" == "" ]; then echo genocat error; exit 1; fi

    if [ "$count" -ne $non_iupac_lines ]; then echo "bad count = $count, expecting $non_iupac_lines"; exit 1; fi

    test_header "genocat --bases ^ACGTN --count (FASTQ)"
    local count=`$genocat_no_echo $output -H --bases ^ACGTN --count -q`
    if [ "$count" -ne 1 ]; then echo "bad count = $count"; exit 1; fi
}

# Test SAM/BAM translations
batch_sam_bam_translations()
{
    batch_print_header

    # test different buddy code path for subsetted file
    test_translate_sam_to_bam_to_sam special.buddy.bam " " -r22

    # test with gencomp
    test_translate_sam_to_bam_to_sam special.depn.bam --force-gencomp

    # note: we have these files in both sam and bam versions generated with samtools
    local files=(special.buddy.bam 
                 special.depn.bam            # depn/prim with/without QUAL
                 special.NA12878.bam 
                 special.pacbio.ccs.bam      # unaligned SAM/BAM with no SQ records
                 special.human2.bam          
                 special.bsseeker2-rrbs.bam) # sam_piz_special_BSSEEKER2_XM sensitive to SAM/BAM
    local file
    for file in ${files[@]}; do
        test_translate_sam_to_bam_to_sam $file
    done
}

# Test SAM/BAM->FQ translations
batch_sam_fq_translations()
{
    batch_print_header

    # note: we have these files in both sam and bam versions generated with samtools
    local files=(special.buddy.bam 
                 special.depn.bam        # depn/prim with/without QUAL
                 special.NA12878.bam 
                 special.pacbio.ccs.bam  # unaligned SAM/BAM with no SQ records
                 special.human2.bam)
    local file
    for file in ${files[@]}; do
        test_translate_sambam_to_fastq $file
        test_translate_sambam_to_fastq ${file%.bam}.sam
    done
}

# Test 23andMe translations
# note: only runs it to see that it doesn't crash, doesn't validate results
batch_23andMe_translations()
{
    batch_print_header

    local file=test.genome_Full.txt
 
    test_header "$file - translate 23andMe to VCF"

    local me23=$TESTDIR/$file
    local vcf=$OUTDIR/copy.vcf.gz

    $genozip -f $me23 -o $output       || exit 1
    $genocat $output -fo $vcf -e $hs37d5 || exit 1

    cleanup
}

batch_phylip_translations()
{
    batch_print_header

    local multifasta=$TESTDIR/basic-multifasta.fa
    local phylip=$OUTDIR/phylip.phy
    local seq=$OUTDIR/seq.fa
    local multifasta2=$OUTDIR/multifasta2.fa

    test_header "$file - translate multifasta to phylip and back"
 
    $genozip $multifasta -fo $output           || exit 1 # compress multifasta
    $genocat $output --sequential --header-one \
        | tr -d "\r" > $seq                    || exit 1 # reconstruct as phylip
    $genocat $output --phylip -fo $phylip      || exit 1 # reconstruct as phylip
    $genozip $phylip -fo $output2              || exit 1 # compress the phylip
    $genocat $output2 --fasta -fo $multifasta2 || exit 1 # reconstruct as multifasta
    cmp_2_files $multifasta2 $seq                        # compare

    cleanup
}

batch_genocat_tests()
{
    batch_print_header

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

    # Phylip to Multi-fasta
    local file=$TESTDIR/basic.phy
    test_count_genocat_lines $file "--fasta" 9

    # FASTQ genocat tests
    file=$TESTDIR/basic.fq
    local num_lines=`grep + $file | wc -l`
    test_count_genocat_lines $file "--header-only" $num_lines 
    test_count_genocat_lines $file "--downsample 2" $(( 4 * $num_lines / 2 )) 
    test_count_genocat_lines "--pair -E $GRCh38 $file $file" "--interleave" $(( 4 * $num_lines * 2 )) 
    test_count_genocat_lines "--pair -E $GRCh38 $file $file" "--interleave --downsample=5,4" $(( 4 * $num_lines / 5 * 2 )) 
    test_count_genocat_lines "--pair -E $GRCh38 $file $file" "--grep PRFX --header-only" 2
    test_count_genocat_lines "--pair -E $GRCh38 $file $file" "--R1" $(( 4 * $num_lines )) 
    test_count_genocat_lines "--pair -E $GRCh38 $file $file" "--R2" $(( 4 * $num_lines ))

    # test --interleave and with --grep
    sed "s/PRFX/prfx/g" $file > $OUTDIR/prfx.fq
    test_count_genocat_lines "--pair -E $GRCh38 $file $OUTDIR/prfx.fq" "--interleave=either --grep PRFX" 8
    test_count_genocat_lines "--pair -E $GRCh38 $file $OUTDIR/prfx.fq" "--interleave=both --grep PRFX" 0

    # grep without pairing
    test_count_genocat_lines $file "--grep line5 --header-only" 1
}

# test --grep, --count, --lines
batch_grep_count_lines()
{
    batch_print_header

    local file
    for file in ${basics[@]}; do
        if [ $file == basic.fa ] || [ $file == basic.bam ] || [ $file == basic.locs ] || [ $file == basic.generic ]; then continue; fi

        if [ $file == basic.chain ]; then
            export GENOZIP_REFERENCE=${hs37d5}:${GRCh38}
        else
            unset GENOZIP_REFERENCE
        fi

        # number of expected lines
        local lines=1
        if [ $file == basic.fq ];    then lines=4; fi
        if [ $file == basic.chain ]; then lines=3; fi

        # grep        
        test_count_genocat_lines $TESTDIR/$file "--grep PRFX --no-header" $lines
        test_count_genocat_lines $TESTDIR/$file "--grep NONEXISTANT --no-header" 0

        # count
        $genozip $TESTDIR/$file -fo $output
        local count=`$genocat_no_echo --quiet --count $output || exit`
        if [ "$count" == "" ]; then echo genocat error; exit 1; fi

        # lines
        test_count_genocat_lines $TESTDIR/$file "--no-header --lines=${count}-${count}" $lines
        test_count_genocat_lines $TESTDIR/$file "--no-header --lines=100000-" 0

        unset GENOZIP_REFERENCE
    done

    # regions-file
    test_count_genocat_lines "$TESTDIR/basic.vcf" "-R $TESTDIR/basic.vcf.regions -H" 7

    # regions
    test_count_genocat_lines "$TESTDIR/basic.vcf" "--regions 13:207237509-207237510,1:207237250 -H" 7
}

batch_backward_compatability()
{
    batch_print_header
    local files=( `ls -r $TESTDIR/back-compat/[0-9]*/*.genozip` )
    local file
    for file in ${files[@]}; do
        test_header "$file - backward compatability test"
        $genounzip -t $file || exit 1
    done
}

batch_prod_compatability()
{
    if [ "$i_am_prod" == "1" ]; then return; fi 

    if [ ! -d ../genozip-prod ]; then return; fi
    
    save_genozip=$genozip
    genozip=$genozip_prod -m

    if (( $1 <= $2 + 0  )) ; then batch_basic basic.vcf     ; fi
    if (( $1 <= $2 + 1  )) ; then batch_dvcf                ; fi
    if (( $1 <= $2 + 2  )) ; then batch_reference_fastq     ; fi
    if (( $1 <= $2 + 3  )) ; then batch_reference_sam       ; fi
    if (( $1 <= $2 + 4  )) ; then batch_basic basic.bam     ; fi
    if (( $1 <= $2 + 5  )) ; then batch_basic basic.fq      ; fi
    if (( $1 <= $2 + 6  )) ; then batch_basic basic.fa      ; fi
    if (( $1 <= $2 + 7  )) ; then batch_basic basic.chain   ; fi
    if (( $1 <= $2 + 8  )) ; then batch_basic basic.gvf     ; fi
    if (( $1 <= $2 + 9  )) ; then batch_basic basic.genome_Full.me23.txt ; fi
    if (( $1 <= $2 + 10 )) ; then batch_kraken " " "-K$kraken"           ; fi
    if (( $1 <= $2 + 11 )) ; then batch_basic basic.phy     ; fi
    if (( $1 <= $2 + 12 )) ; then batch_basic basic.generic ; fi

    genozip=$save_genozip
}
num_batch_prod_compatability_tests=13
    
batch_real_world_1_adler32() # $1 extra genozip argument
{
    batch_print_header

    cleanup # note: cleanup doesn't affect TESTDIR, but we shall use -f to overwrite any existing genozip files

    local filter_out=nothing
    if [ ! -x "$(command -v xz)" ] ; then # xz unavailable
        local filter_out=.xz
    fi

    # without reference
    local files=( `cd $TESTDIR; ls -1 test.*vcf* test.*sam* test.*bam \
                   test.*fq* test.*fa* \
                   basic.phy* test.*gvf* test.*gff* test.*locs* \
                   test.*txt* test.*kraken* | \
                   grep -v "$filter_out" | grep -v headerless | grep -v .genozip | grep -v .md5 | grep -v .bad ` ) 

    for f in ${files[@]}; do rm -f ${f}.genozip; done

    # test genozip and genounzip --test
    echo "subsets of real world files (without reference)"
    test_standard "-f $1 --show-filename" " " ${files[*]}

    # don't remove .genozip files as we will use them in batch_real_world_genounzip
    #for f in ${files[@]}; do rm -f ${f}.genozip; done
}

batch_real_world_genounzip() # $1 extra genozip argument
{
    batch_print_header

    cleanup # note: cleanup doesn't affect TESTDIR, but we shall use -f to overwrite any existing genozip files

    local filter_out=nothing
    if [ ! -x "$(command -v xz)" ] ; then # xz unavailable
        local filter_out=.xz
    fi

    # without reference
    local files=( `cd $TESTDIR; ls -1 test.*vcf* test.*sam* test.*bam \
                   test.*fq* test.*fa* \
                   basic.phy* test.*gvf* test.*gff* test.*locs* \
                   test.*txt* test.*kraken* | \
                   grep -v "$filter_out" | grep -v headerless | grep -v .genozip | grep -v .md5 | grep -v .bad | grep -v .xz | grep -v .bz2` )
    
    # test full genounzip (not --test), including generation of BZGF
    for f in ${files[@]}; do 

        local genozip_file=${TESTDIR}/${f%.gz}.genozip

        # note: normally, the test runs on files compressed in batch_real_world_1_adler32 - we compress them here if not
        if [ ! -f $genozip_file ]; then
            $genozip ${TESTDIR}/$f -fX || exit 1
        fi

        local recon=${OUTDIR}/$f
        $genounzip $genozip_file -o $recon || exit 1

        # same is in private/test/Makefile
        if [[ `head -c2 $recon | od -x | head -1 | cut -d" " -f2 || exit 1` == 8b1f ]]; then
            local actual_md5=`gzip -dc < $recon | md5sum | cut -d" " -f1`
        else
            local actual_md5=`md5sum $recon | cut -d" " -f1`
        fi

        local expected_md5=`cat ${TESTDIR}/${f}.md5` # calculated in test/Makefile
        if [[ "$actual_md5" != "$expected_md5" ]] ; then
            echo "${TESTDIR}/$f has MD5=$expected_md5 but reconstructed file ${OUTDIR}/$f has MD5=$actual_md5"
            exit 1
        fi

        rm -f ${OUTDIR}/$f
    done

    cleanup 
}

batch_real_world_1_backcomp()
{
    batch_print_header

    cleanup # note: cleanup doesn't affect TESTDIR, but we shall use -f to overwrite any existing genozip files

    local filter_out=nothing
    if [ ! -x "$(command -v xz)" ] ; then # xz unavailable
        local filter_out=.xz
    fi

    # without reference
    local files=(  `cd $TESTDIR; ls -1 test.*vcf test.*vcf.gz test.*sam test.*sam.gz test.*bam \
                   test.*fq test.*fq.gz test.*fa test.*fa.gz test.*fasta test.*fasta.gz \
                   basic.phy test.*gvf test.*gvf.gz test.*gff test.*gff.gz test.*locs \
                   test.*txt test.*txt.gz test.*kraken test.*kraken.gz \
                   | grep -v headerless.sam` ) 

    local i=0
    for f in ${files[@]}; do 
        i=$(( i + 1 ))

        # exceptions - supported in 14.0.0 but not earlier
        # if [[ $f == test.3rd-line-sra.fq ]]; then continue; fi 
        # if [[ $f == test.bsseeker2-pysam-qual.bam ]]; then continue; fi
        # if [[ $f == test.covid.gisaid.nuke.fasta.gz ]]; then continue; fi
        # if [[ $f == test.transcriptome.bam ]]; then continue; fi

        # exceptions - supported in 14.0.5 but not earlier
        # if [[ $f == test.novoalign.bam ]]; then continue; fi

        # exceptions - supported in 14.0.7 but not earlier
        # if [[ $f == test.ubam.bam ]]; then continue; fi

        # exceptions - supported in 14.0.9 but not earlier
        if [[ $f == test.MC_Z_in_prim.sam ]]; then continue; fi
        if [[ $f == test.10M-contigs.sam.bz2 ]]; then continue; fi
        if [[ $f == test.unmapped-is-saggy.sam ]]; then continue; fi

        test_header "$f - backward compatability with prod ($i/${#files[@]})"
        $genozip_prod private/test/$f --md5 -fo $output || exit 1
        $genounzip -t $output || exit 1
    done
}

batch_real_world_with_ref_md5() # $1 extra genozip argument
{
    batch_print_header $1

    cleanup # note: cleanup doesn't affect TESTDIR, but we shall use -f to overwrite any existing genozip files

    # with two references
    test_standard "-mf $1 -e $GRCh38 -e $hs37d5" " " test.GRCh38_to_GRCh37.chain 

    local files37=( test.IonXpress.sam.gz \
                    test.human.fq.gz test.human2.bam test.pacbio.clr.bam \
                    test.human2-R1.100K.fq.bz2 test.pacbio.ccs.10k.bam test.unmapped.sam.gz \
                    test.NA12878.chr22.1x.bam test.NA12878-R1.100k.fq test.pacbio-blasr.bam \
                    test.human2.filtered.snp.vcf test.solexa-headerless.sam )

    local files38=( test.1KG-38.vcf.gz test.human-collated-headerless.sam test.human-sorted-headerless.sam )

    local filesT2T1_1=( test.nanopore.t2t_v1_1.bam )

    test_standard "-mf $1 -e $hs37d5 --show-filename" " " ${files37[*]}
    test_standard "-mf $1 -e $GRCh38 --show-filename" " " ${files38[*]}
    test_standard "-mf $1 -e $T2T1_1 --show-filename" " " ${filesT2T1_1[*]}

    for f in ${files37[@]} ${files38[@]} ${filesT2T1_1[@]} test.GRCh38_to_GRCh37.chain; do rm -f ${TESTDIR}/${f}.genozip ; done
}


batch_real_world_with_ref_backcomp()
{
    if [ "$i_am_prod" == "1" ]; then return; fi 

    batch_print_header

    cleanup # note: cleanup doesn't affect TESTDIR, but we shall use -f to overwrite any existing genozip files

    # with a reference
    local files37=( test.IonXpress.sam.gz \
                    test.human.fq.gz test.human2.bam test.pacbio.clr.bam \
                    test.human2-R1.100K.fq.bz2 test.pacbio.ccs.10k.bam \
                    test.NA12878.chr22.1x.bam test.NA12878-R1.100k.fq  \
                    test.human2.filtered.snp.vcf test.solexa-headerless.sam )

    local files38=( test.1KG-38.vcf.gz test.human-collated-headerless.sam )

    local filesT2T1_1=( test.nanopore.t2t_v1_1.bam )

    local total=$(( ${#files37[@]} + ${#files38[@]} + ${#filesT2T1_1[@]} ))

    local i=0
    for f in ${files37[@]}; do     
        i=$(( i + 1 ))
        test_header "$f - backward compatability with prod (with reference) - 37 ($i/$total)"
        $genozip_prod private/test/$f -mf -e $hs37d5 -o $output || exit 1
        $genounzip -t $output || exit 1
    done

    for f in ${files38[@]}; do 
        i=$(( i + 1 ))
        test_header "$f - backward compatability with prod (with reference) - 38 ($i/$total)"
        $genozip_prod private/test/$f -mf -e $GRCh38 -o $output || exit 1
        $genounzip -t $output || exit 1
    done

    for f in ${filesT2T1_1[@]}; do 
        i=$(( i + 1 ))
        test_header "$f - backward compatability with prod (with reference) - T2T ($i/$total)"
        $genozip_prod private/test/$f -mf -e $T2T1_1 -o $output || exit 1
        $genounzip -t $output || exit 1
    done
}

batch_real_world_small_vbs()
{
    batch_print_header

    cleanup # note: cleanup doesn't affect TESTDIR, but we shall use -f to overwrite any existing genozip files

    local filter_out=nothing
    if [ ! -x "$(command -v xz)" ] ; then # xz unavailable
        local filter_out=.xz
    fi

    # lots of small VBs
    local files=( test.IonXpress.sam.gz                                 \
                  test.human.fq.gz test.human2.bam                      \
                  test.human2-R1.100K.fq.bz2 test.pacbio.ccs.10k.bam    \
                  test.pacbio.clr.bam `# multiple PRIM and DEPN vbs`    \
                  test.NA12878.chr22.1x.bam test.NA12878-R1.100k.fq     \
                  test.sequential.fa.gz )

    if [ -x "$(command -v xz)" ] ; then # skip .xz files if xz is not installed
        files+=( test.pacbio.10k.fasta.xz )
    fi

    echo "subsets of real world files (lots of small VBs -B1)"
    test_standard "-mf $1 -B1 --show-filename" " " ${files[*]}

    for f in ${files[@]}; do rm -f ${f}.genozip; done
}

batch_multiseq()
{
    batch_print_header
    test_standard "--multiseq" " " test.coronavirus.fasta

    # regions
    test_count_genocat_lines "$TESTDIR/test.coronavirus.fasta" "--regions MW362225.1" 22
    test_count_genocat_lines "$TESTDIR/test.coronavirus.fasta" "--regions ^MW362225.1" 99978

    test_standard "--multiseq" " " test.nanopore-virus.fq
}

# CRAM hs37d5
batch_external_cram()
{
    batch_print_header
    if `command -v samtools >& /dev/null`; then
        test_standard "-E$hs37d5" " " test.human2.cram   
    fi
}

# BCF
batch_external_bcf()
{
    batch_print_header
    if `command -v bcftools >& /dev/null`; then
        test_standard " " " " test.human2.filtered.snp.bcf    
    fi
}

# unzip
batch_external_unzip()
{
    batch_print_header
    if `command -v unzip >& /dev/null`; then
        test_standard " " " " test.genome_Full.zip    
    fi
}

batch_reference_fastq()
{
    batch_print_header

    echo "paired FASTQ with --reference, --password (BZ2), --md5, -B1"
    test_standard "CONCAT -e$GRCh38 -p 123 --pair -mB1" "-p123" test.human2-R1.100K.fq.bz2 test.human2-R2.100K.fq.bz2

    echo "4 paired FASTQ with --REFERENCE (BGZF, decompress concatenated, password)"
    test_standard "COPY -E$GRCh38 -2 -p 123" " " test.human2-R1.100K.fq.gz test.human2-R2.100K.fq.gz

    # solexa read style
    test_standard "-e$GRCh38 --pair" "" special.solexa-R1.fq special.solexa-R2.fq
}

batch_reference_sam()
{
    batch_print_header

    echo "command line with mixed SAM and FASTQ files with --reference"
    echo "Note: '$GRCh38' needs to be up to date with the latest genozip format"
    test_standard "-me$GRCh38" " " test.human2.bam test.human.fq.gz test.human3-collated.bam

    echo "multiple SAM with --REFERENCE" 
    test_standard "-mE$GRCh38" " " test.human-collated-headerless.sam test.human3-collated.bam
    
    echo "SAM with --REFERENCE and --password" 
    test_standard "-E$GRCh38 --password 123" "-p123" test.human-collated-headerless.sam

    echo "BAM with --reference and --password, alternate chrom names" 
    test_standard "-me$hg19 --password 123" "-p123 -e$hg19" test.human2.bam  

    echo "SAM with large (>4GB) plant genome"  
    test_standard "-me$chinese_spring --password 123" "-p123 -e$chinese_spring" test.bsseeker2-wgbs.sam.gz  
}

batch_reference_vcf()
{
    batch_print_header

    echo "multiple VCF with --reference, --md5 using hs37d5 ; alternate chroms names"
    test_standard "COPY -me$hg19" " " test.human2.filtered.snp.vcf

    echo "GVCF with --reference, --md5 using GRCh38"
    test_standard "-me$GRCh38" " " test.g.vcf.gz

    echo "multiple VCF with --REFERENCE using hs37d5, password" 
    test_standard "-mE$hs37d5 -p123" "--password 123" test.1KG-37.vcf test.human2.filtered.snp.vcf
}

batch_many_small_files()
{
    batch_print_header

    cleanup

    local num_files=300
    echo "Generating $num_files small files"
    mkdir ${OUTDIR}/smalls
    for i in `seq $num_files`; do cp ${TESTDIR}/minimal.sam ${OUTDIR}/smalls/minimal.${i}.sam; done

    $genozip -ft -D ${OUTDIR}/smalls || exit 1

    cleanup
}

batch_make_reference()
{
    batch_print_header

    cleanup

    # Making a reference
    echo "Making a reference"
    local fa_file=$REFDIR/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz 
    local ref_file=$OUTDIR/output.ref.genozip
    test_header "$genozip --make-reference $fa_file"
    $genozip --make-reference $fa_file --force -o $ref_file || exit 1

    local ref="--reference $ref_file"
    local REF="--REFERENCE $ref_file"

    echo "unaligned SAM with --reference"
    test_standard "$ref" "$ref" basic-unaligned.sam 

    echo "unaligned SAM with --REFERENCE"
    test_standard "$REF" " " basic-unaligned.sam

    echo "unalignable SAM with --REFERENCE"
    test_standard "$REF" " " basic-unalignable.sam

    echo "unaligned BAM with --reference"
    test_standard "$ref" "$ref" basic-unaligned.bam

    echo "FASTQ with --REFERENCE"
    test_standard "$REF" " " basic.fq 

#    if [ -z "$is_windows" ]; then # in windows, we don't support redirecting stdin
        echo "unaligned SAM with --REFERENCE - from stdin"
        test_redirected basic-unaligned.sam "$REF"
#    fi

    #cleanup - no cleanup, we need the reference for batch_reference_backcomp
}

# compress with prod and prod-created ref file, uncompress with new version and new version-created ref file
# ref files are expected have the same MD5
batch_reference_backcomp()
{
    batch_print_header    
    
    local fa_file=$REFDIR/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz 
    local ref_file=$OUTDIR/output.ref.genozip
    local prod_ref_file=$OUTDIR/output.prod.ref.genozip

    # new ref file normally created by batch_make_reference, but we create if for some reason its not
    if [ ! -f $ref_file ]; then
        echo "Making reference"
        $genozip --make-reference $fa_file --force -o $ref_file || exit 1
    fi

    if [ ! -f $prod_ref_file ]; then
        echo "Making prod reference"
        $genozip_prod --make-reference $fa_file --force -o $prod_ref_file || exit 1
    fi

    local files38=( test.human.fq.gz test.human-collated-headerless.sam test.1KG-38.vcf.gz )

    for f in ${files38[@]}; do 
        test_header "$f - reference file backward compatability with prod"
        $genozip_prod private/test/$f -mf -e $prod_ref_file -o $output || exit 1
        $genounzip -t $output -e $ref_file || exit 1
    done

    cleanup
}

# compress headerless SAM with wrong ref. Expected to work (with bad compression ratio)
batch_headerless_wrong_ref()
{
    batch_print_header

    $genozip -ft ${TESTDIR}/test.human-collated-headerless.sam -e $hs37d5 || exit 1 
    
    cleanup
}


batch_genols()
{
    batch_print_header

    $genozip ${TESTDIR}/basic.fq ${TESTDIR}/basic.fq -2 -e $GRCh38 -fo $output -p abcd || exit 1
    $genols $output -p abcd || exit 1
    rm -f $output
}

batch_tar_files_from()
{
    batch_print_header
    cleanup

    local tar=${OUTDIR}/output.tar

    $genozip -D -T ${TESTDIR}/basic-files-from -f --tar $tar  || exit 1
    tar xvf $tar || exit 1

    cat ${TESTDIR}/basic-files-from-genozip | $genounzip --files-from - -t || exit 1
    $genols --files-from ${TESTDIR}/basic-files-from-genozip || exit 1
    $genocat --files-from ${TESTDIR}/basic-files-from-genozip > $output || exit 1
    
    cleanup
}

batch_gencomp_depn_methods() # note: use --debug-gencomp for detailed tracking
{
    batch_print_header

    # invoke REREAD with plain file (test.pacbio.clr.bam is NOT compressed with BGZF)
    # -B1 forces multiple depn VBs
    # -@3 sets DEPN queue length to 3 forcing other VBs to be REREAD
    $genozip -fB1 -@3 -t $TESTDIR/test.pacbio.clr.bam || exit 1 

    # invoke REREAD with BGZF file (test.pacbio.clr.bam.gz is generated by test/Makefile)
    $genozip -fB1 -@3 -t $TESTDIR/test.pacbio.clr.bam.gz || exit 1

    # invoke OFFLOAD method
    $genozip -fB1 -@3 -t file://${path}${TESTDIR}/test.pacbio.clr.bam || exit 1
    cat $TESTDIR/test.pacbio.clr.bam.gz | $genozip -i bam -fB1 -@3 -t - -o $output || exit 1

    cleanup
}

output=${OUTDIR}/output.genozip
output2=${OUTDIR}/output2.genozip
recon=${OUTDIR}/recon.txt
kraken=${OUTDIR}/kraken.genozip
chain=${OUTDIR}/chain.genozip

is_windows="`uname|grep -i mingw``uname|grep -i MSYS`"
is_mac=`uname|grep -i Darwin`

# reference and chain files
hg19=$REFDIR/hg19.p13.plusMT.full_analysis_set.ref.genozip
hs37d5=$REFDIR/hs37d5.ref.genozip
GRCh38=$REFDIR/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
T2T1_1=$REFDIR/chm13.draft_v1.1.ref.genozip
chinese_spring=$REFDIR/161010_Chinese_Spring_v1.0_pseudomolecules_parts.ref.genozip
chain37_38=$REFDIR/GRCh37_to_GRCh38.chain.genozip

if (( $# < 1 )); then
    echo "Usage: test.sh [debug|opt|prod] <batch_id-test> [optional-genozip-arg]"
    exit 0
fi

# debug, opt, prod
is_debug=`echo $1|grep debug`
if [ -n "$is_debug" ]; then 
    debug_zip=-debug
    debug_nonzip=-debug
    shift
fi

is_opt=`echo $1|grep opt`
if [ -n "$is_opt" ]; then 
    debug_zip=-opt
    debug_nonzip=-opt
    shift
fi

# test prod (for a maintainence release)
is_prod=`echo $1|grep prod`
if [ -n "$is_prod" ]; then 
    dir=../genozip-prod
    shift
fi

if [ ! -n "$dir" ]; then 
    dir=.
fi

if (( `pwd | grep genozip-prod | wc -l` == 1 )); then
    i_am_prod=1;
fi

# -----------------
# platform settings
# -----------------
if [ -n "$is_windows" ]; then
    genozip_exe=$dir/genozip${debug_zip}.exe
    genozip_prod_exe=$dir/genozip-prod.exe
    genounzip_exe=$dir/genounzip${debug_nonzip}.exe
    genocat_exe=$dir/genocat${debug_nonzip}.exe
    genols_exe=$dir/genols${debug_nonzip}.exe 
#    zip_threads="-@3"
#    piz_threads="-@5"
    path=`pwd| cut -c3-|tr / '\\\\'`\\
else
    genozip_exe=$dir/genozip${debug_zip}
    genozip_prod_exe=$dir/genozip-prod
    genounzip_exe=$dir/genounzip${debug_nonzip}
    genocat_exe=$dir/genocat${debug_nonzip}
    genols_exe=$dir/genols${debug_nonzip} 
    path=$PWD/
fi

genozip="$genozip_exe --echo $2 $zip_threads"
genozip_prod="$genozip_prod_exe --echo $2 $zip_threads"
genounzip="$genounzip_exe --echo $2 $piz_threads"
genocat_no_echo="$genocat_exe $2 $piz_threads"
genocat="$genocat_exe --echo $2 $piz_threads"
genols=$genols_exe 

basics=(basic.vcf basic.chain basic.sam basic.vcf basic.bam basic.fq basic.fa basic.gvf basic.genome_Full.me23.txt \
        basic.kraken basic.phy basic.locs basic.generic)

exes=($genozip_exe $genounzip_exe $genocat_exe $genols_exe)
for exe in ${exes[@]}; do
    if [ ! -x $exe ]; then
        echo "Error: $exe does not exist"
        exit 1
    fi
done

if `command -v md5 >& /dev/null`; then
    md5="md5 -q" # mac
else
    md5=md5sum 
fi

mkdir $OUTDIR >& /dev/null
cleanup

make -j -C $TESTDIR --quiet || exit 1

# only if doing a full test (starting from 0) - delete genome and hash caches
sparkling_clean()
{
    rm -f ${hg19}.*cache* ${hs37d5}.*cache* ${GRCh38}.*cache* ${TESTDIR}/*.genozip ${TESTDIR}/basic-subdirs/*.genozip ${TESTDIR}/*rejects* ${TESTDIR}/*.DEPN
}

# unfortunately Mac's bash doesn't support "case" with fall-through ( ;& )
batch_id=$1
batch_id=$((batch_id - 1))

if (( $1 <= 0  )) ; then  sparkling_clean              ; fi
if (( $1 <= 1  )) ; then  batch_minimal                ; fi
if (( $1 <= 2  )) ; then  batch_basic basic.vcf        ; fi
if (( $1 <= 3  )) ; then  batch_basic basic.bam        ; fi
if (( $1 <= 4  )) ; then  batch_basic basic.sam        ; fi
if (( $1 <= 5  )) ; then  batch_basic basic.fq         ; fi
if (( $1 <= 6  )) ; then  batch_basic basic.fa         ; fi
if (( $1 <= 7  )) ; then  batch_basic basic.chain      ; fi
if (( $1 <= 8  )) ; then  batch_basic basic.gvf        ; fi
if (( $1 <= 9  )) ; then  batch_basic basic.genome_Full.me23.txt ; fi
if (( $1 <= 10 )) ; then  batch_basic basic.kraken     ; fi
if (( $1 <= 11 )) ; then  batch_basic basic.phy        ; fi
if (( $1 <= 12 )) ; then  batch_basic basic.generic    ; fi
if (( $1 <= 13 )) ; then  batch_precompressed          ; fi
if (( $1 <= 14 )) ; then  batch_bgzf                   ; fi
if (( $1 <= 15 )) ; then  batch_subdirs                ; fi
if (( $1 <= 16 )) ; then  batch_special_algs           ; fi
if (( $1 <= 17 )) ; then  batch_dvcf                   ; fi
if (( $1 <= 18 )) ; then  batch_sam_bam_translations   ; fi
if (( $1 <= 19 )) ; then  batch_sam_fq_translations    ; fi
if (( $1 <= 20 )) ; then  batch_23andMe_translations   ; fi
if (( $1 <= 21 )) ; then  batch_phylip_translations    ; fi
if (( $1 <= 22 )) ; then  batch_genocat_tests          ; fi
if (( $1 <= 23 )) ; then  batch_grep_count_lines       ; fi
if (( $1 <= 24 )) ; then  batch_backward_compatability ; fi
if (( $1 <= 25 )) ; then  batch_match_chrom            ; fi
if (( $1 <= 26 )) ; then  batch_kraken " " "-K$kraken" ; fi   # genocat loads kraken data
if (( $1 <= 27 )) ; then  batch_kraken "-K$kraken" " " ; fi   # genozip loads kraken data
if (( $1 <= 28 )) ; then  batch_single_thread          ; fi 
if (( $1 <= 29 )) ; then  batch_copy_ref_section       ; fi 
if (( $1 <= 30 )) ; then  batch_iupac                  ; fi 
if (( $1 <= 31 )) ; then  batch_genols                 ; fi
if (( $1 <= 32 )) ; then  batch_tar_files_from         ; fi
if (( $1 <= 33 )) ; then  batch_gencomp_depn_methods   ; fi 
if (( $1 <= 34 )) ; then  batch_real_world_small_vbs   ; fi 
if (( $1 <= 35 )) ; then  batch_real_world_1_adler32   ; fi 
if (( $1 <= 36 )) ; then  batch_real_world_genounzip   ; fi 
if (( $1 <= 37 )) ; then  batch_real_world_1_adler32 "--best -f" ; fi 
if (( $1 <= 38 )) ; then  batch_real_world_1_adler32 --fast    ; fi 
if (( $1 <= 39 )) ; then  batch_real_world_with_ref_md5; fi 
if (( $1 <= 40 )) ; then  batch_real_world_with_ref_md5 --best ; fi 
if (( $1 <= 41 )) ; then  batch_multiseq               ; fi
if (( $1 <= 42 )) ; then  batch_external_cram          ; fi
if (( $1 <= 43 )) ; then  batch_external_bcf           ; fi
if (( $1 <= 44 )) ; then  batch_external_unzip         ; fi
if (( $1 <= 45 )) ; then  batch_reference_fastq        ; fi
if (( $1 <= 46 )) ; then  batch_reference_sam          ; fi
if (( $1 <= 47 )) ; then  batch_reference_vcf          ; fi
if (( $1 <= 48 )) ; then  batch_many_small_files       ; fi
if (( $1 <= 49 )) ; then  batch_make_reference         ; fi
if (( $1 <= 50 )) ; then  batch_reference_backcomp     ; fi
if (( $1 <= 51 )) ; then  batch_headerless_wrong_ref   ; fi

if [[ `basename $PWD` != genozip-prod ]]; then
    if (( $1 <= 52 )) ; then  batch_real_world_1_backcomp  ; fi 
    if (( $1 <= 53 )) ; then  batch_real_world_with_ref_backcomp ; fi 
    next=54
    if (( $1 <= $next + $num_batch_prod_compatability_tests )) ; then batch_prod_compatability $1 $next ; fi
fi

printf "\nALL GOOD!\n"
