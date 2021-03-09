.. _idxstats:

idxstats
========

Data Types: SAM, BAM and FASTQ

**Usage**

``genocat --idxstats my-file.bam.genozip`` 

``genocat --idxstats my-file.fq.genozip`` 

**Description**

Generates the list of contigs, along with number of mapped and unmapped reads for each contig. Reads with an undefined contig are grouped under "*".

The output format and contents are identical to `samtools idxstats <http://www.htslib.org/doc/samtools-idxstats.html>`_.

This works both on SAM/BAM and on FASTQ. 

For FASTQ, the mapping to contigs is as reported by the Genozip Aligner. The Genozip Aligner maps reads for compression purposes and does not attempt to map them according to the biological truth. However, usually the large majority of reads are in fact mapped to their correct position, so this can give a reasonable approximation of idxstats of the data directly from FASTQ without needing to map it to BAM.

**Output**

A tab-seperated list:

The columns are:
   #. Contig name
   #. Contig length
   #. Number of mapped reads
   #. Number of unmapped reads

::

    genocat --idxstats my-file.bam.genozip

    chr1    248956422       22502721        242067
    chr2    242193529       22856412        237624
    chr3    198295559       18923298        182644
    chr4    190214555       18505304        178501
    chr5    181538259       16884321        169270

**Example**

Alternative sex assignment with `sexassign.py <https://github.com/grahamgower/sexassign>`_ (on Linux):

``genocat mysample.bam.genozip --idxstats | sexassign.py /dev/fd/0``
