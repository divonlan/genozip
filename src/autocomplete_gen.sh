#!/usr/bin/env bash
# ^ finds bash according to $PATH
 
# ------------------------------------------------------------------
#   autocomplete_gen.sh
#   Copyright (C) 2024-2024 Genozip Limited. Patent Pending.
#   Please see terms and conditions in the file LICENSE.txt

# Tutorial: https://www.gnu.org/software/gnuastro/manual/html_node/Bash-TAB-completion-tutorial.html
# $1 = The name of the command
# $2 = The current word being completed (empty unless we are in the middle of typing a word).
# $3 = The word before the word being completed.

if [[ `basename $PWD` != src ]]; then
    echo "autocomplete_gen.sh must be run from the src directory"
    exit 1
fi

if [[ "$1" != debug ]]; then
    # remove diagnostic options, and also alternative names for options (eg list-chroms==contigs, license==licence etc)
    diagnostics_genozip="|^debug|^show-|^dump-|^biopsy|^bgzf|head|check-latest|force-[PLy|deep|domq|homp|longr|normq|pacb|smux]|hold-cache|licence|^optimise|license-prepare|no-[header|domq|eval|homp|longr|pacb|smux|zriter|faf]|STATS|analyze-insertions|skip-segconf|submit-stats|^t_|verify-codec|xthreads|log-digest|chroms|no-FAF|no-interleaved"
    diagnostics_genounzip="|^debug|^show-|^dump-|STATS|chroms|licence|log-digest|make-reference|no-eval|^recover$|^t_|verify-codec|xthreads|decompress|hold-cache|tar"
    diagnostics_genocat="|^debug|^show-|^dump-|STATS|chroms|licence|log-digest|one-vb|^FASTQ|^FQ|^fq|^BAM|^SAM|^VCF|^BCF|decompress|gpos|hold-cache|make-reference|no-[eval|pg]|^qname$|qname-file|^recover$|verify-codec|xthreads"
    diagnostics_genols="|^debug|^show-|licence|no-eval"
    out=autocomplete.sh
else
    out=autocomplete-debug.sh
    func_suffix=_debug
fi

is_windows="`uname|grep -i mingw``uname|grep -i MSYS`"
if [ -n "$is_windows" ]; then exe=.exe ; fi

if type ./genozip$exe >& /dev/null ; then
    unset GENOZIP_TEST
    fasta_exts=`./genozip$exe --help=input | grep FASTA | cut -d: -f2 | xargs | tr " " "\n" | sed "s/^/./g" | sed "s/$/$/g" | tr "\n" "|" ; echo -n "/$"`
    inputs=`./genozip$exe --help=input | egrep "\.gz|\.zip|generic" | cut -d: -f2 | xargs`
else
    fasta_exts=".fasta$|.fasta.gz$|.fasta.bz2$|.fasta.xz$|.fasta.zip$|.faa$|.faa.gz$|.faa.bz2$|.faa.xz$|.faa.zip$|.ffn$|.ffn.gz$|.ffn.bz2$|.ffn.xz$|.ffn.zip$|.fnn$|.fnn.gz$|.fnn.bz2$|.fnn.xz$|.fnn.zip$|.fna$|.fna.gz$|.fna.bz2$|.fna.xz$|.fna.zip$|.fna$|.frn.gz$|.frn.bz2$|.frn.xz$|.frn.zip$|.fas$|.fas.gz$|.fas.bz2$|.fas.xz$|.fas.zip$|.fsa$|.fsa.gz$|.fsa.bz2$|.fsa.xz$|.fsa.zip$|.fa$|.fa.gz$|.fa.bz2$|.fa.xz$|.fa.zip$|/$"
    inputs="vcf vcf.gz vcf.bgz vcf.bz2 vcf.xz sam sam.gz sam.bgz sam.bz2 sam.xz fastq fastq.gz fastq.bz2 fastq.xz fastq.ora fq fq.gz fq.bz2 fq.xz fq.ora fasta fasta.gz fasta.bz2 fasta.xz fasta.zip faa faa.gz faa.bz2 faa.xz faa.zip ffn ffn.gz ffn.bz2 ffn.xz ffn.zip fnn fnn.gz fnn.bz2 fnn.xz fnn.zip fna fna.gz fna.bz2 fna.xz fna.zip fna frn.gz frn.bz2 frn.xz frn.zip fas fas.gz fas.bz2 fas.xz fas.zip fsa fsa.gz fsa.bz2 fsa.xz fsa.zip fa fa.gz fa.bz2 fa.xz fa.zip gff3 gff3.gz gff3.bz2 gff3.xz gff gff.gz gff.bz2 gff.xz gvf gvf.gz gvf.bz2 gvf.xz gtf gtf.gz gtf.bz2 gtf.xz bam bam.gz bam.bgz cram bcf bcf.gz bcf.bgz locs locs.gz locs.bz2 locs.xz bed bed.gz bed.bz2 bed.xz track track.gz track.bz2 track.xz 23andme 23andme.zip generic"
fi

rm -f $out

echo "# bash autocompletion script for Genozip commands" >> $out
echo "# " >> $out
echo "# Copyright (C) 2024-2024 Genozip Limited." >> $out
echo >> $out
echo "# For this to work, bash autocomplete must be already installed on your computer (it usually already is): " >> $out
echo "#    Ubuntu: https://www.cyberciti.biz/faq/add-bash-auto-completion-in-ubuntu-linux/" >> $out
echo "#    Mac: https://www.simplified.guide/macos/bash-completion" >> $out
echo "#" >> $out
echo "# To enable Genozip autocompletion, add the following line to your ~/.bash_completion file (this is done automatically when you register to genozip):"  >> $out
echo '# source <...your-path-to-genozip...>/autocomplete.sh ; complete -F _genozip genozip; complete -F _genounzip genounzip ; complete -F _genocat genocat ; complete -F _genols genols'  >> $out
echo '# this will run when a new shell is opened, but you can also run it explicitly: source ~/.bash_completion' >> $out
echo >> $out

for cmd in genozip genounzip genocat genols; do    
    diagnostics=$( eval "echo \$diagnostics_$cmd" )
    filter=$( eval "echo \$filter_$cmd" )

    # function header
    echo "_${cmd}${func_suffix}() {" >> $out

    # get last short option (which might have a argument)
    unset x >> $out
    echo 'x="$2"' >> $out
    echo 'if [[ "${x:0:1}" == '"'-'"' && "${x:1:1}" != '"'-'"' ]]; then' >> $out
    echo '     short_option=${x:0-1}' >> $out
    echo 'fi' >> $out
    
    # show command line options relevant for the cmd after --
    echo -n '    options="' >> $out
    egrep_arg=$( grep ${cmd}_lo flags.c | tr " " "\n" | egrep "^_"  | sed "s/,$//g" | tr "\n" "|" ; echo "@@@@@" ) # @@@@@ is the 2nd argument to final | 

    grep "#define" flags.c | egrep -w "$egrep_arg" | cut -d\" -f2 | egrep -v "_00$diagnostics" | sed "s/^/--/g" | tr "\n" " " >> $out    
    echo '"' >> $out

    echo '    if [[ "$2" =~ ^--* ]]; then' >> $out 
    echo '        COMPREPLY=( `compgen -W "$options" -- "$2" `) ' >> $out 

    # show reference files after e.g. --reference (except genols): files from current dir and also $GENOZIP_REFERENCE dir
    if [ $cmd != genols ] ; then 
    echo '    elif [[ "$3" =~ "-".*[eE]$ ]] || [ "$3" == "--reference" ] || [ "$3" == "--REFERENCE" ]; then' >> $out 
    echo '        COMPREPLY=( `ls -d1p ${2}* 2>/dev/null | egrep ".ref.genozip$|/$" ; if [ -v GENOZIP_REFERENCE ] && [ -d $GENOZIP_REFERENCE ]; then cd $GENOZIP_REFERENCE; ls -d1p ${2}* 2>/dev/null | egrep ".ref.genozip$" ; cd - > /dev/null; fi` )' >> $out
    fi

    if [ $cmd == genozip ] ; then 
    # show input options after --input or -i
    echo '    elif [[ "$3" =~ "-".*i || "$3" == "--input" ]]; then' >> $out 
    echo '        COMPREPLY=( `compgen -W "'$inputs'" -- "$2"` )' >> $out 
    
    # show tar files after --tar
    echo '    elif [[ "$3" == "--tar" ]]; then' >> $out 
    echo '        COMPREPLY=( `ls -dp ${2}* 2>/dev/null | xargs | tr " " "\n" | egrep "\.tar$|/$"` )' >> $out

    # show FASTA files after --make-reference
    echo '    elif [[ "$3" == "--make-reference" ]]; then' >> $out 
    echo '        COMPREPLY=( `ls -dp ${2}* 2>/dev/null | xargs | tr " " "\n" | egrep "'$fasta_exts'"` )' >> $out

    # show directories and all files except .genozip files (good for --user-message argument too)
    echo '    else' >> $out 
    echo '        COMPREPLY=( `ls -d1p ${2}* 2>/dev/null | grep -v ".genozip$"` )' >> $out

    # genounzip: show directories and .genozip files, excluding .ref.genozip
    elif [ $cmd == genounzip ] ; then 
    echo '    else' >> $out 
    echo '        COMPREPLY=( `ls -d1p ${2}* 2>/dev/null | egrep "\.genozip$|/$" | grep -Fv  ".ref.genozip"` )' >> $out

    # genocat: show directories and .genozip files, including .ref.genozip
    elif [ $cmd == genocat ] ; then 
    echo '    else' >> $out 
    echo '        COMPREPLY=( `ls -d1p ${2}* 2>/dev/null | egrep "\.genozip$|/$"` )' >> $out
    fi

    echo '    fi' >> $out 

    # if completing a directory, remove trailing space
    echo '    [[ $COMPREPLY == */ ]] && compopt -o nospace ' >> $out

    echo '}' >> $out
    echo >> $out
done
