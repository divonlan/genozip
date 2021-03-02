Example of a FASTQ to BAM pipeline using Genozip
================================================

This pipeline starts with a ``.fq.genozip`` file, cleans it up with ``fastp``, maps it with ``bwa mem``, sorts it with ``bamsort``, removes duplicates with ``bamstreamingmarkduplicates`` and finally outputs a ``.bam.genozip`` file.

This pipeline has Genozip files on both ends, and uses no intermediate files, so it is efficient in storage usage. 

Genozip adds very little overhead, as its CPU consumption is insignificant in comparison with the tools used in this pipeline, in particular ``bwa mem``.

**Notes**
  | • The input ``fastq.genozip`` files were originally created from pairs of fq.gz files:
  |   ``genozip sample01-R1.fq.gz sample01-R2.fq.gz --pair``
  |
  | • ``genocat --interleave`` is used to generate FASTQ data consisting of reads from R1 interleaved with reads from R2, which feed into ``fastp --interleaved_in``. This works because we used ``--pair`` to generate the ``.fq.genozip`` file.
  |
  | • ``genozip`` uses ``GRCh38_full_analysis_set_plus_decoy_hla.fa`` which was prepared with:
  |   ``genozip --make-reference GRCh38_full_analysis_set_plus_decoy_hla.fa``
  |
  | • This script is designed to have multiple instances of it running in parallel. Each instance will work on different files. 
  |   *Technical note*: ``${out}.doing_now`` is used as a mutex and ``touch``/``rm`` as mutex_lock/unlock. Since ``touch`` is not atomic, exclusivity is not 100% guarateed, but in practice it seems to work. 

::

    #!/bin/bash

    ref=GRCh38_full_analysis_set_plus_decoy_hla.fa
    study=mystudy
    fastq=myfastqdir
    mapped=mymappeddir

    files=($fastq/*.genozip)

    processed=1

    while (( processed == 1 )); do
        processsed=0
        for file in ${files[@]}
        do
            sample=`echo $file |grep -o -E sample'[[:digit:]]{2}'`  # convert the file name to a sample name
            out=$mapped/${sample}.bam.genozip

            if [ -f $out ]; then continue; fi # already processed
            if [ -f ${out}.doing_now ]; then continue; fi # another instance of this script is working on it

            processed=1
            touch ${out}.doing_now

            echo =========================================
            echo Sample $sample
            echo =========================================

            ( genocat --interleave $file -e ${ref%.fa}.ref.genozip                                                  || >&2 echo "genocat exit=$?" )|\
            ( fastp --stdin --interleaved_in --stdout --html ${fastq}/${sample}.html --json ${fastq}/${sample}.json || >&2 echo "fastp exit=$?"   )|\
            ( bwa mem $ref - -p -t 54 -T 0 -R "@RG\tID:$sample\tSM:$study\tPL:Illumina"                             || >&2 echo "bwa exit=$?"     )|\
            ( samtools view -h -OSAM                                                                                || >&2 echo "samtools exit=$?")|\
            ( bamsort fixmates=1 adddupmarksupport=1 inputformat=sam outputformat=sam inputthreads=5 outputthreads=5 sortthreads=30 level=1  || >&2 echo "bamsort exit=$?" )|\
            ( bamstreamingmarkduplicates inputformat=sam inputthreads=3 outputthreads=3 level=1                     || >&2 echo "bamstreamingmarkduplicates exit=$?" )|\
            ( genozip -e ${ref%.fa}.ref.genozip -i bam -o $out -t                                                   || >&2 echo "genozip exit=$?" )

            rm ${out}.doing_now
        done
    done

