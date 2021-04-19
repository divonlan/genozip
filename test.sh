#!/bin/bash 

shopt -s extglob  # Enables extglob - see https://mywiki.wooledge.org/glob

TESTDIR=test
OUTDIR=$TESTDIR/tmp

cleanup() { 
    rm -f $OUTDIR/* $TESTDIR/*.bad $TESTDIR/*.rejects.* 
}

cmp_2_files() {
#    if (( `$md5 $1 ${2%.*} | cut -d" " -f1 | uniq | wc -l` != 1 )) ; then

    if (( `$md5 $1 $2 | cut -d" " -f1 | uniq | wc -l` != 1 )) ; then
        echo "MD5 comparison FAILED: $1 $2"
#        $md5 $1 ${2%.*}
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
    test_header "$1 - redirected from stdin"
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
    test_header "$1 - bind & unbind (2 files with 2 components each)"
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

    # test --comoponent
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
    echo "******* ${FUNCNAME[1]} (batch_id=${batch_id}) *******"
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
    local files=(basic.vcf basic.sam basic.fq basic.fa basic.gvf basic.genome_Full.me23.txt basic.chain basic.kraken basic.phy)
    local file
    for file in ${files[@]}; do
        
        test_unix_style $file
        test_windows_style $file
        test_standard "NOPREFIX CONCAT" " " file://${path}${TESTDIR}/$file
        test_standard "-p123" "--password 123" $file
        if [ $file != basic.phy ]; then # issue with redirection on Windows of Phylip files (bug 339)
            test_redirected $file
            test_stdout $file
        fi
        test_standard "COPY" " " $file
        test_multi_bound $file REPLACE # REPLACE to adjust the contig name for .fa as we can't have two contigs with the same name
        test_optimize $file
    done
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
    local files=(basic.bam basic-bgzf-6.sam.gz basic-bgzf-9.sam.gz basic-bgzf-6-no-eof.sam.gz basic-1bgzp_block.bam)
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


batch_dual_coordinates()
{
    batch_print_header

    local files=(test.human2.filtered.snp.vcf test.NA12878.sorted.vcf test.ExAC.vcf.gz)
    local chain=data/GRCh37_to_GRCh38.chain.genozip
    local file

    for file in ${files[@]}; do
        test_header "$file - dual-coordintes test"

        local src=$OUTDIR/${file}
        local primary=$OUTDIR/primary.vcf
        local luft=$OUTDIR/luft.vcf
        local primary2=$OUTDIR/primary2.vcf

        # compare src to primary (ignoring header and INFO fields)
        echo -n "Step 1: make ${primary}.genozip from $file : " 
        $genozip -C $chain test/$file -fo ${primary}.genozip || exit 1
        echo -n "Step 2: make ${primary} from ${primary}.genozip : " 
        $genocat ${primary}.genozip --no-pg -fo ${primary}

        # convert primary -> luft -> primary
        echo -n "Step 3: make ${luft} from ${primary}.genozip : " 
        $genocat --luft --no-pg ${primary}.genozip -fo ${luft} || exit 1
        echo -n "Step 4: make ${luft}.genozip from ${luft} : " 
        $genozip $luft -fo ${luft}.genozip || exit 1
        echo -n "Step 5: make ${primary2} from ${luft}.genozip : " 
        $genocat ${luft}.genozip --no-pg -fo ${primary2} || exit 1
        echo "Step 6: compare $primary to $primary2" 
        cmp $primary $primary2 || exit 1

        cleanup
        rm -f ${src}.noinfo ${primary}.genozip  ${primary}.noinfo ${luft} ${luft}.genozip ${primary2}
    done
}

batch_kraken() # $1 genozip arguments #2 genocat (one of them must include --kraken)
{
    batch_print_header

    local files=(test/basic-test-kraken.fa test/basic-test-kraken.fq test/basic.kraken test/basic.sam) # SAM must be last
    local file

    $genozip test/basic.kraken -i kraken -fo $kraken

    # testing filtering FASTA, FASTQ, SAM and KRAKEN itself with --taxid 
    for file in ${files[@]}; do
        test_header "$file - kraken test (genozip $1 ; genocat $2)"

        echo -n genozip $file : 
        $genozip $1 $file -fo $output || exit 1

        local chars_plus=`($genocat $2 -Hq $output -k570 || exit 1) | wc -c`
        local chars_minus=`($genocat $2 -Hq $output -k^570 || exit 1) | wc -c`
        local chars=`($genocat -Hq $output || exit 1) | wc -c`
        echo "$file : chars_plus=$chars_plus chars_minus=$chars_minus chars=$chars"
        if [ $(($chars_plus + $chars_minus)) -ne $chars ]; then
            echo "$file: adding up kraken positive- and negative- filtered data, isn't the size of the original file"
            exit 1
        fi        
    done

    # testing filtering SAM translated to FASTQ with --taxid 
    file=test/basic.sam
    test_header "$file --fastq - kraken test"

    local chars_plus=`($genocat $2 --fastq -Hq $output -k570 || exit 1) | wc -c`
    local chars_minus=`($genocat $2 --fastq -Hq $output -k^570 || exit 1) | wc -c`
    local chars=`($genocat --fastq -Hq $output || exit 1) | wc -c`
    echo "$file --fastq : chars_plus=$chars_plus chars_minus=$chars_minus chars=$chars"
    if [ $(($chars_plus + $chars_minus)) -ne $chars ]; then
        echo "$file: adding up kraken positive- and negative- filtered data, isn't the size of the original file"
        exit 1
    fi        
    
    # testing filtering SAM with a concatenated KRAKEN file (representing slightly different classifications
    # originating from separate kraken2 of R1 and R2 FASTQ files)
    file=test/basic.sam
    test_header "$file concatenated kraken filter test"

    $genozip test/basic.kraken test/basic-2nd-file.kraken -fo $kraken

    local chars_plus=`($genocat $2 -Hq $output -k570 || exit 1) | wc -c`
    local chars_minus=`($genocat $2 -Hq $output -k^570 || exit 1) | wc -c`
    local chars=`($genocat -Hq $output || exit 1) | wc -c`
    echo "$file : chars_plus=$chars_plus chars_minus=$chars_minus chars=$chars"
    if [ $(($chars_plus + $chars_minus)) -ne $chars ]; then
        echo "$file: adding up kraken positive- and negative- filtered data, isn't the size of the original file"
        exit 1
    fi        
    
    cleanup
}

# Test SAM/BAM translations
batch_sam_translations()
{
    batch_print_header

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
    if [ ! -f $sambam ] ; then echo "$sambam: File not found"; exit 1; fi

    $genozip -f $me23 -o $output       || exit 1
    $genocat $output -fo $vcf -e $hg19 || exit 1

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
    cmp $multifasta2 $seq                      || exit 1 # compare

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
    test_count_genocat_lines $file "--header-only" `grep @ $file | wc -l` 
    test_count_genocat_lines $file "--grep line5" 4
    test_count_genocat_lines $file "--grep line5 --header-only" 1
    test_count_genocat_lines $file "--downsample 2" $(( 4 * `grep @ $file | wc -l` / 2 )) 
    test_count_genocat_lines "--pair -E $GRCh38 $file $file" "--interleave" $(( 4 * `grep @ $file | wc -l` * 2 )) 
    test_count_genocat_lines "--pair -E $GRCh38 $file $file" "--interleave --downsample=5,4" $(( 4 * `grep @ $file | wc -l` / 5 * 2 )) 
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

batch_real_world_subsets()
{
    batch_print_header

    cleanup # note: cleanup doesn't affect TESTDIR, but we shall use -f to overwrite any existing genozip files

    filter_out=nothing
    if [ ! -x "$(command -v xz)" ] ; then # xz unavailable
        filter_out=.xz
    fi

    local files=( `cd test; ls -1 test.*vcf* test.*sam* test.*bam* \
                   test.*fq* test.*fastq* test.*fa* test.*fasta* \
                   basic.phy* test.*chain* test.*gvf* \
                   test.*txt* test.*kraken* | \
                   grep -v "$filter_out" | grep -v .genozip` )

    echo "subsets (~3 VBs) or real world files"
    test_standard "-mf $1" " " ${files[*]}
}

batch_multifasta()
{
    batch_print_header
    test_standard "--multifasta" " " test.coronavirus.fasta
}

batch_misc_cases()
{
    batch_print_header

    # Test binding SAM files with lots of contigs (no reference)
    echo "binding SAM files with lots of contigs (no reference)"
    test_multi_bound test.human-unsorted.sam
}

# CRAM hg19
batch_external_cram()
{
    batch_print_header
    if `command -v samtools >& /dev/null`; then
        batch_print_header
        test_standard "-E$hg19" " " test.human2.cram   
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

batch_reference()
{
    batch_print_header

    echo "paired FASTQ with --reference, --password"
    test_standard "CONCAT -e$GRCh38 -p 123 --pair" "-p123" test.human2-R1.100K.fq.bz2 test.human2-R2.100K.fq.bz2

    echo "4 paired FASTQ with --REFERENCE (BGZF, decompress concatenated)"
    test_standard "COPY -E$GRCh38 --pair" " " test.human2-R1.100K.fq.gz test.human2-R2.100K.fq.gz

    echo "4 paired FASTQ with --REFERENCE (BZ2, decompress unbound)"
    test_standard "COPY CONCAT -E$GRCh38 -2" " " test.human2-R1.100K.fq.bz2 test.human2-R2.100K.fq.bz2

    echo "command line with mixed SAM and FASTQ files with --reference"
    echo "Note: '$GRCh38' needs to be up to date with the latest genozip format"
    test_standard "-me$GRCh38" "-e$GRCh38" test.human-unsorted.sam test.human.fq.gz test.human-sorted.sam

    echo "multiple bound SAM with --REFERENCE" 
    test_standard "-mE$GRCh38" " " test.human-unsorted.sam test.human-sorted.sam
    
    echo "SAM with --reference and --password" 
    test_standard "-me$GRCh38 --password 123" "-p123 -e$GRCh38" test.human-unsorted.sam
    
    echo "SAM with --REFERENCE and --password" 
    test_standard "-E$GRCh38 --password 123" "-p123" test.human-unsorted.sam
    
    echo "multiple bound VCF with --reference, --md5 using hg19, and unbind"
    test_standard "COPY CONCAT -me$hg19" " " test.human2-R1.100K.fq.bz2 test.human2-R2.100K.fq.bz2

    echo "multiple VCF with --REFERENCE using hg19" 
    test_standard "-mE$hg19" " " test.ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf test.human2.filtered.snp.vcf
}

batch_make_reference()
{
    batch_print_header

    cleanup

    # Making a reference
    echo "Making a reference"
    local fa_file=data/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz 
    local ref_file=$OUTDIR/output.ref.genozip
    $genozip --make-reference $fa_file --force -o $ref_file || exit 1

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

batch_genols()
{
    batch_print_header

    $genozip test/basic.sam test/minimal.sam -fo $output -p abcd || exit 1
    $genols $output -p abcd || exit 1
    rm -f $output
}

make --quiet testfiles

output=${OUTDIR}/output.genozip
output2=${OUTDIR}/output2.genozip
recon=${OUTDIR}/recon.txt
kraken=${OUTDIR}/kraken.genozip

is_windows=`uname|grep -i mingw`
is_mac=`uname|grep -i Darwin`

hg19=data/hs37d5.ref.genozip
GRCh38=data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip

if (( $# < 1 )); then
    echo "Usage: test.sh [debug] <batch_id-test> [optional-genozip-arg]"
    exit 0
fi

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
    genozip_exe=./genozip${debug}.exe
    genounzip_exe=./genounzip${debug}.exe
    genocat_exe=./genocat${debug}.exe
    genols_exe=./genols${debug}.exe 
    path=`pwd| cut -c3-|tr / '\\\\'`\\
else
    genozip_exe=./genozip${debug}
    genounzip_exe=./genounzip${debug}
    genocat_exe=./genocat${debug}
    genols_exe=./genols${debug} 
    path=$PWD/
fi

genozip="$genozip_exe -@3 --echo $2"
genounzip="$genounzip_exe --echo -@5 $2"
genocat="$genocat_exe --echo -@5 $2"
genols=$genols_exe 

exes=($genozip_exe $genounzip_exe $genocat_exe $genols_exe)
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

# only if doing a full test - delete genome and hash caches
if (( $1 == 1 )) ; then
    rm -f ${hg19}.*cache* ${GRCh38}.*cache* 
fi

# unfortunately Mac's bash doesn't support "case" with fall-through ( ;& )
batch_id=$1
batch_id=$((batch_id - 1))

if (( $1 <= 1  )) ; then  batch_minimal                ; fi
if (( $1 <= 2  )) ; then  batch_basic                  ; fi
if (( $1 <= 3  )) ; then  batch_precompressed          ; fi
if (( $1 <= 4  )) ; then  batch_bgzf                   ; fi
if (( $1 <= 5  )) ; then  batch_special_algs           ; fi
if (( $1 <= 6  )) ; then  batch_dual_coordinates       ; fi
if (( $1 <= 7  )) ; then  batch_sam_translations       ; fi
if (( $1 <= 8  )) ; then  batch_23andMe_translations   ; fi
if (( $1 <= 9  )) ; then  batch_phylip_translations    ; fi
if (( $1 <= 10 )) ; then  batch_genocat_tests          ; fi
if (( $1 <= 11 )) ; then  batch_backward_compatability ; fi
if (( $1 <= 12 )) ; then  batch_kraken " " "-K$kraken" ;      # genocat loads kraken data
                          batch_kraken "-K$kraken" " " ; fi   # genozip loads kraken data
if (( $1 <= 13 )) ; then  batch_real_world_subsets     ; fi ; # natural VB size
if (( $1 <= 14 )) ; then  batch_real_world_subsets -B1 ; fi ; # many VBs
if (( $1 <= 15 )) ; then  batch_multifasta             ; fi
if (( $1 <= 16 )) ; then  batch_misc_cases             ; fi
if (( $1 <= 17 )) ; then  batch_external_cram          ; fi
if (( $1 <= 18 )) ; then  batch_external_bcf           ; fi
if (( $1 <= 19 )) ; then  batch_external_unzip         ; fi
if (( $1 <= 20 )) ; then  batch_reference              ; fi
if (( $1 <= 21 )) ; then  batch_make_reference         ; fi
if (( $1 <= 22 )) ; then  batch_genols                 ; fi

printf "\nALL GOOD!\n"


