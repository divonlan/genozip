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

cleanup() { 
    rm -f $OUTDIR/* $TESTDIR/*.bad $TESTDIR/*.rejects.* 
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
            echo "Error: compressed $num_files files, but only $count genozip files found in td. Files compressed: "
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
    
    $genozip ${file} -fo $output || exit 1
    ($genocat --no-pg $output || exit 1) | tr -d "\r" > $OUTDIR/unix-nl.$1 

    cmp_2_files $file $OUTDIR/unix-nl.$1
    cleanup
}

test_multi_bound() # $1=filename $2=REPLACE (optional)
{
    test_header "$1 - test_multi_bound - bind & unbind (2 files with 2 components each)"
    local file=$TESTDIR/$1
    local file1=$OUTDIR/copy1.$1
    local file2=$OUTDIR/copy2.$1

    cp -f $file $file1 || exit 1
    if [[ $2 == "REPLACE" ]]; then
        cat $file | sed 's/PRFX/FIL2/g' > $file2 || exit 1 # note - FIL2 needs to be the same length of PRFX or basic.phy will break
    else
        cp -f $file $file2 || exit 1
    fi

    $genozip $file1 $file2 -ft -o $output || exit 1 # test as bound
    cp -f $output $output2 || exit 1
    $genounzip $output $output2 -t || exit 1 # test unbind 2x2

    # test --component
    if [[ $2 == "REPLACE" ]]; then

        $genocat $output --component 1 -fo $recon || exit 1
        local wc=`cat $recon | grep PRFX | wc -l`
        if (( "$wc" == 0 )); then echo "FAILED --component 1 - expected 1 lines, but getting $wc" ; exit 1; fi

        $genocat $output --component 2 -fo $recon || exit 1
        local wc=`cat $recon | grep FIL2 | wc -l`
        if (( "$wc" == 0 )); then echo "FAILED --component 2 - expected 1 lines, but getting $wc" ; exit 1; fi
    fi  

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

test_translate_bam_to_sam() # $1 bam file 
{
    test_header "$1 - translate BAM to SAM"

    local bam=$TESTDIR/$1
    local sam=${bam%.bam}.sam
    local copy=$OUTDIR/copy.sam
    if [ ! -f $bam ] ; then echo "$bam: File not found"; exit 1; fi
    if [ ! -f $sam ] ; then echo "$sam: File not found"; exit 1; fi

    $genozip -f $bam -o $output  || exit 1
    $genocat $output --no-PG -fo $copy  || exit 1
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

    $genozip -f $sam -o $output  || exit 1
    $genocat $output --bam --no-PG -fo $new_bam || exit 1
    
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

test_backward_compatability()
{
    test_header "$1 - backward compatability test"
    $genounzip -t $1 || exit 1
}

batch_print_header()
{
    batch_id=$((batch_id + 1))
    echo "***************************************************************************"
    echo "******* ${FUNCNAME[1]} (batch_id=${batch_id}) "
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
    if [ $file != basic.bam ] && [ $file != basic.generic ] && [ $file != basic.phy ]; then # issue with redirection on Windows of Phylip files (bug 339)
        test_redirected $file
        test_stdout $file
    fi
    test_standard "COPY $ref" " " $file

    test_multi_bound $file $replace # REPLACE to adjust the contig name for .fa as we can't have two contigs with the same name
    test_optimize $file
    unset GENOZIP_REFERENCE
}

# pre-compressed files (except BGZF) and non-precompressed BAM
batch_precompressed()
{
    batch_print_header
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
    batch_print_header
    local files=(basic-bgzf.bam basic-bgzf-6.sam.gz basic-bgzf-9.sam.gz basic-bgzf-6-no-eof.sam.gz basic-1bgzp_block.bam)
    #local files=()
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

    test_header "sam -> sam.genozip -> sam.gz - see that it is BGZF"
    local sam_gz=${OUTDIR}/bgzf_test.sam.gz
    $genozip ${TESTDIR}/basic.sam -fo $output || exit 1
    $genocat --no-pg $output -fo $sam_gz || exit 1
    if [ "$(head -c4 $sam_gz | od -x | head -1 | awk '{$2=$2};1')" != "0000000 8b1f 0408" ] ; then  # the awk converts the Mac output to be the same as Linux (removing redundant spaces)
        echo $sam_gz is not BGZF-compressed
        exit 1
    fi

    test_header "sam-> sam.genozip -> bam - see that it is BGZF"
    local bam=${OUTDIR}/bgzf_test.bam
    $genocat --no-pg $output -fo $bam || exit 1
    if [ "$(head -c4 $bam | od -x | head -1 | awk '{$2=$2};1')" != "0000000 8b1f 0408" ] ; then 
        echo $bam is not BGZF-compressed
        exit 1
    fi

    test_header "sam.gz -> sam.genozip -> sam - see that it is not BGZF"
    local sam=${OUTDIR}/bgzf_test.sam
    $genozip $sam_gz -fo $output || exit 1
    $genocat --no-pg $output -fo $sam || exit 1
    if [ "$(head -c4 $sam | od -x | head -1 | awk '{$2=$2};1')" == "0000000 8b1f 0408" ] ; then 
        echo $sam is unexpectedly BGZF-compressed
        exit 1
    fi
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
        test_multi_bound $file               # multiple files bound
        test_optimize $file                  # optimize - only compress to see that it doesn't error
    done
}

batch_dvcf()
{
    batch_print_header

    local files=(basic-dvcf-source.vcf basic-dvcf-luft.vcf test.NA12878.sorted.vcf test.clinvar37.vcf.gz test.chr22.indels.vcf test.chr17.SS6004478.vcf test.ExAC.vcf.gz)
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
        $genocat $dvcf --no-pg -fo ${primary}

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
        $genounzip ${one}.genozip || exit 1

        #convert CHROM_STYLE_chr22 and then to CHROM_STYLE_22
        $genozip --match-chrom $file -fo ${two}.genozip -e $hs37d5 || exit 1
        $genounzip -f ${two}.genozip || exit 1

        $genozip --match-chrom $two -fo ${three}.genozip -e $hg19 || exit 1
        $genounzip -f ${three}.genozip || exit 1

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

    local lines_plus=`$genocat $output -Hq ${cat1_args[*]} --count`
    if [ "$lines_plus" == "" ] || [ "$lines_plus" -eq 0 ]; then echo "genocat error - \$lines_plus=\"$lines_plus\""; exit 1; fi

    local lines_minus=`$genocat $output -Hq $3 --count`
    if [ "$lines_minus" == "" ] || [ "$lines_minus" -eq 0 ]; then echo "genocat error - \$lines_minus=\"$lines_minus\""; exit 1; fi
    
    local lines=`$genocat -Hq $output --count`
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
    $genozip ${TESTDIR}/basic.kraken ${TESTDIR}/basic-2nd-file.kraken -fo $kraken

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

    # SAM
    non_iupac_lines=$(( `grep -v "^@" ${TESTDIR}/basic.sam | wc -l` - 1 ))
    test_count_genocat_lines ${TESTDIR}/basic.sam "-H --bases=AGCTN" $non_iupac_lines
    test_count_genocat_lines ${TESTDIR}/basic.sam "-H --bases=^AGCTN" 1

    # BAM
    test_header "genocat --bases AGCTN --count --bam"
    local count=`$genocat $output --echo -H --bam --bases AGCTN --count -q`
    if [ "$count" == "" ]; then echo genocat error; exit 1; fi

    if [ "$count" -ne $non_iupac_lines ]; then echo "bad count = $count"; exit 1; fi

    test_header "genocat --bases ^AGCTN --count --bam"
    local count=`$genocat $output --echo -H --bam --bases ^AGCTN --count -q`
    if [ "$count" -ne 1 ]; then echo "bad count = $count"; exit 1; fi

    # FASTQ
    test_count_genocat_lines ${TESTDIR}/basic.fq "-H --IUPAC=AGCTN" 20
    test_count_genocat_lines ${TESTDIR}/basic.fq "-H --IUPAC=^AGCTN" 4
}

# Test SAM/BAM translations
batch_sam_translations()
{
    batch_print_header

    # note: we have these files in both sam and bam versions generated with samtools
    local files=(test.NA12878.chr22.1x.bam 
                 test.pacbio.ccs.10k.bam  # unaligned SAM/BAM with no SQ records
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
    test_count_genocat_lines $file "--header-only" `grep + $file | wc -l` 
    test_count_genocat_lines $file "--downsample 2" $(( 4 * `grep + $file | wc -l` / 2 )) 
    test_count_genocat_lines "--pair -E $GRCh38 $file $file" "--interleave" $(( 4 * `grep + $file | wc -l` * 2 )) 
    test_count_genocat_lines "--pair -E $GRCh38 $file $file" "--interleave --downsample=5,4" $(( 4 * `grep + $file | wc -l` / 5 * 2 )) 
    test_count_genocat_lines "--pair -E $GRCh38 $file $file" "--grep PRFX --header-only" 2

    # test --interleave and with --grep
    sed "s/PRFX/prfx/g" $file > $OUTDIR/prfx.fq
    test_count_genocat_lines "--pair -E $GRCh38 $file $OUTDIR/prfx.fq" "--interleave=either --grep PRFX" 8
    test_count_genocat_lines "--pair -E $GRCh38 $file $OUTDIR/prfx.fq" "--interleave=both --grep PRFX" 0

    test_count_genocat_lines $file "--grep line5 --header-only" 1
}

# test --grep, --count, --lines
batch_grep_count_lines()
{
    batch_print_header

    local file
    for file in ${basics[@]}; do
        if [ $file == basic.fa ] || [ $file == basic.bam ] || [ $file == basic.generic ]; then continue; fi

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
        local count=`$genocat --quiet --count $output || exit`
        if [ "$count" == "" ]; then echo genocat error; exit 1; fi

        # lines
        test_count_genocat_lines $TESTDIR/$file "--no-header --lines=${count}-${count}" $lines
        test_count_genocat_lines $TESTDIR/$file "--no-header --lines=100000-" 0

        unset GENOZIP_REFERENCE
    done

    # regions-file
    test_count_genocat_lines "$TESTDIR/basic.vcf" "-R $TESTDIR/basic.vcf.regions -H" 6

    # regions
    test_count_genocat_lines "$TESTDIR/basic.vcf" "--regions 13:207237509-207237510,1:207237250 -H" 6
}

batch_backward_compatability()
{
    batch_print_header
    local files=`ls $TESTDIR/back-compat/[0-9]*/*.genozip` 
    local file
    for file in $files; do
        test_backward_compatability $file
    done
}

batch_prod_compatability()
{
    if [ ! -d genozip-prod ]; then return; fi

    batch_print_header

    (cd ../genozip-prod; make ${debug:1})
    
    save_genozip=$genozip
    genozip=../genozip-prod/$genozip

    batch_basic basic.vcf        
    #batch_dvcf
    batch_reference_fastq
    batch_reference_sam
    batch_basic basic.bam        
    batch_basic basic.fq         
    batch_basic basic.fa         
    batch_basic basic.chain      
    batch_basic basic.gvf        
    batch_basic basic.genome_Full.me23.txt 
    batch_kraken " " "-K$kraken"     
    batch_basic basic.phy        
    batch_basic basic.generic    

    genozip=$save_genozip
}
    
batch_real_world_1()
{
    batch_print_header

    cleanup # note: cleanup doesn't affect TESTDIR, but we shall use -f to overwrite any existing genozip files

    local filter_out=nothing
    if [ ! -x "$(command -v xz)" ] ; then # xz unavailable
        local filter_out=.xz
    fi

    # without reference
    local files=( `cd $TESTDIR; ls -1 test.*vcf* test.*sam* test.*bam* \
                   test.*fq* test.*fa* \
                   basic.phy* test.*gvf* test.*gff* \
                   test.*txt* test.*kraken* | \
                   grep -v "$filter_out" | grep -v .genozip` )

    for f in $files; do rm -f ${f}.genozip; done

    echo "subsets of real world files (without reference)"
    test_standard "-mf $1" " " ${files[*]}

    for f in $files; do rm -f ${f}.genozip; done
}

batch_real_world_with_ref()
{
    batch_print_header

    cleanup # note: cleanup doesn't affect TESTDIR, but we shall use -f to overwrite any existing genozip files

    # with two references
    test_standard "-mf $1 -e $GRCh38 -e $hs37d5" " " test.GRCh38_to_GRCh37.chain 

    # with a reference
    local files=( test.IonXpress.sam \
                  test.human.fq.gz test.human2.bam test.human2.sam \
                  test.human2-R1.100K.fq.bz2 test.pacbio.ccs.10k.bam test.pacbio.ccs.10k.sam \
                  test.NA12878.chr22.1x.bam test.NA12878-R1.100k.fq test.pacbio.10k.hg19.sam.gz )

    echo "subsets of real world files (with reference)"
    test_standard "-mf $1 -e $hs37d5" " " ${files[*]}

    for f in $files test.GRCh38_to_GRCh37.chain; do rm -f ${TESTDIR}/${f}.genozip ; done
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
    local files=( test.IonXpress.sam \
                  test.human.fq.gz test.human2.bam test.human2.sam \
                  test.human2-R1.100K.fq.bz2 test.pacbio.ccs.10k.bam \
                  test.NA12878.chr22.1x.bam test.NA12878-R1.100k.fq \
                  test.sequential.fa.gz )

    if [ -x "$(command -v xz)" ] ; then # skip .xz files if xz is not installed
        files+=( test.pacbio.10k.fasta.xz )
    fi

    echo "subsets of real world files (lots of small VBs -B1)"
    test_standard "-mf $1 -B1" " " ${files[*]}

    for f in $files; do rm -f ${f}.genozip; done
}

batch_multifasta()
{
    batch_print_header
    test_standard "--multifasta" " " test.coronavirus.fasta

    # regions
    test_count_genocat_lines "$TESTDIR/test.coronavirus.fasta" "--regions MW362225.1" 22
    test_count_genocat_lines "$TESTDIR/test.coronavirus.fasta" "--regions ^MW362225.1" 99978
}

batch_misc_cases()
{
    batch_print_header

    # Test binding SAM files with lots of contigs (no reference)
    echo "binding SAM files with lots of contigs (no reference)"
    test_multi_bound test.human-collated.sam
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

    echo "paired FASTQ with --reference, --password"
    test_standard "CONCAT -e$GRCh38 -p 123 --pair" "-p123" test.human2-R1.100K.fq.bz2 test.human2-R2.100K.fq.bz2

    echo "4 paired FASTQ with --REFERENCE (BGZF, decompress concatenated)"
    test_standard "COPY -E$GRCh38 --pair" " " test.human2-R1.100K.fq.gz test.human2-R2.100K.fq.gz

    echo "4 paired FASTQ with --REFERENCE (BZ2, decompress unbound) and password"
    test_standard "COPY CONCAT -E$GRCh38 -2 -p 123" "-p123" test.human2-R1.100K.fq.bz2 test.human2-R2.100K.fq.bz2
}

batch_reference_sam()
{
    batch_print_header

    echo "command line with mixed SAM and FASTQ files with --reference"
    echo "Note: '$GRCh38' needs to be up to date with the latest genozip format"
    test_standard "-me$GRCh38" " " test.human-collated.sam test.human.fq.gz test.human-sorted.sam

    echo "multiple bound SAM with --REFERENCE" 
    test_standard "-mE$GRCh38" " " test.human-collated.sam test.human-sorted.sam
    
    echo "SAM with --REFERENCE and --password" 
    test_standard "-E$GRCh38 --password 123" "-p123" test.human-collated.sam

    echo "SAM with --reference and --password, alternate chrom names" 
    test_standard "-me$hg19 --password 123" "-p123 -e$hg19" test.human2.sam    
}

batch_reference_vcf()
{
    batch_print_header

    echo "multiple bound VCF with --reference, --md5 using hs37d5, and unbind ; alternate chroms names"
    test_standard "COPY CONCAT -me$hg19" " " test.human2.filtered.snp.vcf

    echo "multiple VCF with --REFERENCE using hs37d5" 
    test_standard "-mE$hs37d5" " " test.ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf test.human2.filtered.snp.vcf
}

batch_make_reference()
{
    batch_print_header

    cleanup

    # Making a reference
    echo "Making a reference"
    local fa_file=data/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz 
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

    echo "unaligned SAM with --REFERENCE - from stdin"
    test_redirected basic-unaligned.sam "$REF"

    cleanup
}

batch_genols()
{
    batch_print_header

    $genozip ${TESTDIR}/basic.sam ${TESTDIR}/minimal.sam -fo $output -p abcd || exit 1
    $genols $output -p abcd || exit 1
    rm -f $output
}

batch_tar_files_from()
{
    batch_print_header

    local tar=${OUTDIR}/output.tar

    $genozip -T ${TESTDIR}/basic-files-from -f --tar $tar  || exit 1
    tar xvf $tar || exit 1
    
    cat ${TESTDIR}/basic-files-from-genozip | $genounzip --files-from - -t || exit 1
    $genols --files-from ${TESTDIR}/basic-files-from-genozip || exit 1
    $genocat --files-from ${TESTDIR}/basic-files-from-genozip -fo $output || exit 1
    rm -fR $tar
}

output=${OUTDIR}/output.genozip
output2=${OUTDIR}/output2.genozip
recon=${OUTDIR}/recon.txt
kraken=${OUTDIR}/kraken.genozip
chain=${OUTDIR}/chain.genozip

is_windows=`uname|grep -i mingw`
is_mac=`uname|grep -i Darwin`

if [ -n "$is_windows" ]; then
    make --quiet testfiles || exit 1
fi

# standard file - test.sh should change these
hg19=data/hg19.p13.plusMT.full_analysis_set.ref.genozip
hs37d5=data/hs37d5.ref.genozip
GRCh38=data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
chain37_38=data/GRCh37_to_GRCh38.chain.genozip

if (( $# < 1 )); then
    echo "Usage: test.sh [debug|opt|prod] <batch_id-test> [optional-genozip-arg]"
    exit 0
fi

# debug and opt
is_debug=`echo $1|grep debug`
if [ -n "$is_debug" ]; then 
    debug=-debug; 
    shift
fi

is_opt=`echo $1|grep opt`
if [ -n "$is_opt" ]; then 
    debug=-opt; 
    shift
fi

is_prod=`echo $1|grep prod`
if [ -n "$is_prod" ]; then 
    dir=../genozip-prod
    shift
fi

if [ ! -n "$dir" ]; then 
    dir=.
fi

# -----------------
# platform settings
# -----------------
if [ -n "$is_windows" ]; then
    genozip_exe=$dir/genozip${debug}.exe
    genounzip_exe=$dir/genounzip${debug}.exe
    genocat_exe=$dir/genocat${debug}.exe
    genols_exe=$dir/genols${debug}.exe 
#    zip_threads="-@3"
#    piz_threads="-@5"
    path=`pwd| cut -c3-|tr / '\\\\'`\\
else
    genozip_exe=$dir/genozip${debug}
    genounzip_exe=$dir/genounzip${debug}
    genocat_exe=$dir/genocat${debug}
    genols_exe=$dir/genols${debug} 
    path=$PWD/
fi

genozip="$genozip_exe --echo $2 $zip_threads"
genounzip="$genounzip_exe --echo $2 $piz_threads"
genocat="$genocat_exe --echo $2 $piz_threads"
genols=$genols_exe 

basics=(basic.vcf basic.chain basic.sam basic.vcf basic.bam basic.fq basic.fa basic.gvf basic.genome_Full.me23.txt \
        basic.kraken basic.phy basic.generic)

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

# only if doing a full test (starting from 0) - delete genome and hash caches
sparkling_clean()
{
    rm -f ${hg19}.*cache* ${hs37d5}.*cache* ${GRCh38}.*cache* ${TESTDIR}/*.genozip ${TESTDIR}/*rejects*
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
if (( $1 <= 15 )) ; then  batch_special_algs           ; fi
if (( $1 <= 16 )) ; then  batch_dvcf                   ; fi
if (( $1 <= 17 )) ; then  batch_sam_translations       ; fi
if (( $1 <= 18 )) ; then  batch_23andMe_translations   ; fi
if (( $1 <= 19 )) ; then  batch_phylip_translations    ; fi
if (( $1 <= 20 )) ; then  batch_genocat_tests          ; fi
if (( $1 <= 21 )) ; then  batch_grep_count_lines       ; fi
if (( $1 <= 22 )) ; then  batch_backward_compatability ; fi
if (( $1 <= 23 )) ; then  batch_match_chrom            ; fi
if (( $1 <= 24 )) ; then  batch_kraken " " "-K$kraken" ; fi   # genocat loads kraken data
if (( $1 <= 25 )) ; then  batch_kraken "-K$kraken" " " ; fi   # genozip loads kraken data
if (( $1 <= 26 )) ; then  batch_single_thread          ; fi 
if (( $1 <= 27 )) ; then  batch_copy_ref_section       ; fi 
if (( $1 <= 28 )) ; then  batch_iupac                  ; fi 
if (( $1 <= 29 )) ; then  batch_real_world_small_vbs   ; fi 
if (( $1 <= 30 )) ; then  batch_real_world_1           ; fi 
if (( $1 <= 31 )) ; then  batch_real_world_with_ref    ; fi 
if (( $1 <= 32 )) ; then  batch_multifasta             ; fi
if (( $1 <= 33 )) ; then  batch_misc_cases             ; fi
if (( $1 <= 34 )) ; then  batch_external_cram          ; fi
if (( $1 <= 35 )) ; then  batch_external_bcf           ; fi
if (( $1 <= 36 )) ; then  batch_external_unzip         ; fi
if (( $1 <= 37 )) ; then  batch_reference_fastq        ; fi
if (( $1 <= 38 )) ; then  batch_reference_sam          ; fi
if (( $1 <= 39 )) ; then  batch_reference_vcf          ; fi
if (( $1 <= 40 )) ; then  batch_genols                 ; fi
if (( $1 <= 41 )) ; then  batch_tar_files_from         ; fi
if (( $1 <= 42 )) ; then  batch_make_reference         ; fi
if (( $1 <= 43 )) ; then  batch_prod_compatability     ; fi

printf "\nALL GOOD!\n"
