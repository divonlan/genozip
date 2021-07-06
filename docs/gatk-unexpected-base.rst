.. _gatk-unexpected-base:

Handling GATK's "Unexpected base in allele bases" error
=======================================================

When running GATK's HaplotypeCaller (and perhaps other commands?) on a BAM files that contains bases other than A,C,T,G (and N?), GATK throws an exception:

``java.lang.IllegalArgumentException: Unexpected base in allele bases``

This was observed both in GATK 3.5 and 4.1.

Genozip can directly filter the offending lines out of a BAM file:

Step 1: Compress the file with Genozip:

::

    genozip myfile.bam

Step 2: Filter it:

- Keep only lines in which SEQ consists of only A,C,T,G,N:

::

    genocat myfile.bam.genozip --bases ACGTN

- See the offending lines:

::

    genocat myfile.bam.genozip --bases ^ACGTN

- Count the number of offending lines:

::

    genocat myfile.bam.genozip --bases ^ACGTN --count

| This also works for SAM and FASTQ files.
|
| The list of IUPAC chacacters can be found here: `IUPAC codes <https://www.bioinformatics.org/sms/bases.html>`_




