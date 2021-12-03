.. _archiving:

Archiving
=========

Genozip has two separate capabilities of aggregating files, each with its pros and cons:

**1. tar with --tar**

Using the ``genozip --tar`` option, genozip compresses files directly into a standard `tar file <https://en.wikipedia.org/wiki/Tar_(computing)>`_. 

Each file is compressed independently and written directly into a standard tar file as it is being formed. This is faster and consumes less disk space than first genozipping files and then packaging them into a tar file, since no separate .genozip files are created - just the tar file. 

Example 1:

::

    > # Compressing
    > genozip --tar mydata.tar sample1.bam sample2.bam variants.vcf

    > # Listing the contents of the tar file
    > tar tvf mydata.tar
    -rw-rw-rw- USER/USER   3424847 2021-06-01 11:34 sample1.bam.genozip
    -rw-rw-rw- USER/USER   6765323 2021-03-04 22:04 sample2.bam.genozip
    -rw-rw-rw- USER/USER    765323 2021-03-04 22:08 variants.vcf.genozip
    
    > # Unarchiving and decompressing all files
    > tar xvf mydata.tar |& genounzip --files-from - --replace

Example 2: compress all files in a directory and its sub-directories, using ``--subdirs``:

::

    > genozip --tar mydata.tar --subdirs my-data-dir

Example 3: compress and archive all BAM files in the current directory and its sub-directories, preserving the directory struture:

::

    > find . -name "*.bam" | genozip --tar mydata.tar -T-

Implementation note: Genozip implements the IEEE 1003.1-1988 ("ustar") standard of tar files, with the size field in binary format for files 8GB or larger. This is compatible with most modern tar implementations, including GNU tar.

**2. Binding with --output**

| Binding files only works for files that are similar: 
| - The files must be of the same type (eg VCF, SAM/BAM, FASTQ etc) 
| - In case of VCF files: they must have the same samples 
| - In case of SAM/BAM files: they must have identical @SQ header lines (contigs). 

Unlike ``--tar``, binding does not record the directory structure.

Unlike tar that merely aggregates files, binding leverages similarities and data redundancies between files to improve compression - typically, the size of a bound file is smaller than the combined sizes of files compressed separately.

Binding also allows analyzing the combined file as a concatenated data set.

Example:

::

    > # binding
    > genozip chr1.vcf chr2.vcf --output variants.vcf.genozip

    > # unbinding back to the original files
    > genounzip variants.vcf.genozip

    > # viewing the concatenated data
    > genocat variants.vcf.genozip

    > # extracting the first component (chr1.vcf in this case) from the bound file
    > genocat --component 1 variants.vcf.genozip

    > # listing the components of a bound file
    > genols variants.vcf.genozip


In the following example, we bind all the BAM files in the current directory:

::

    > # binding
    > genozip *.bam --output mydata.bam.genozip 

    > # unbinding back to the original files
    > genounzip mydata.bam.genozip

    > # unbinding into a different directory
    > genounzip mydata.bam.genozip --prefix=mydir/

     