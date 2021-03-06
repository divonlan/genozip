.. _coverage:

Coverage and Depth
==================

Data Types: SAM, BAM and FASTQ

**Usage**

``genocat --show-coverage my-file.bam.genozip``

``genocat --show-coverage=all my-file.bam.genozip`` 

``genocat --show-coverage my-file.fq.genozip`` 

**Description**

Displays the per-chromosome Number-of-reads, Number-of-bases, % of total bases that are on this chromosome and Depth. With the paramater =all, shows all contigs, not just primary chromosomes.

Fine details: 

  *Number-of-bases* is the number of bases mapped to a contig, excluding bases with a ``S`` (Soft Clip) CIGAR (see 1.4.6 `here <https://samtools.github.io/hts-specs/SAMv1.pdf>`_).

  *Number-of-reads* is the number of reads mapped to a contig.
  
  *Depth* is the Coverage divided by the Length (LN) of the contig.
  
  *Length* is defined in the header of the SAM/BAM, or in the reference file used to compress it, or in the case of a header-less SAM file that was compressed without a reference - the highest coordinate of the contig that appears in the file.

  *Contigs that are chromosomes* are contigs with a name of up to 5 characters. For example ``chr22`` is, but ``chr22_KI270731v1_random`` is not.

  *Number-of-bases* and *Number-of-reads* excludes unmapped reads, and in SAM/BAM also excludes reads with a FLAG indicating Unmapped, Duplicate, Secondary and Failed Filters. These are reported separately.

  *Unmapped* reads (e.g. with FLAG 0x4) are counted as "Unmapped" and not as reads of a contig, even if their RNAME contains a value.
  
**Output**
    
Output is space-seperated if sent to a terminal, and tab-seperated if redirected to a file or pipe.

::

    $ genocat.exe --show-coverage my-sample.bam.genozip
    Contig         LN        Reads        Bases       Bases  Depth
    chr1           249.0 Mb  24,907,479   3312.8 Mb   7.0 %  13.31
    chr2           242.2 Mb  23,976,680   3240.4 Mb   6.9 %  13.38
    chr3           198.3 Mb  19,735,323   2677.6 Mb   5.7 %  13.50
    chr4           190.2 Mb  19,445,892   2612.9 Mb   5.5 %  13.74
    chr5           181.5 Mb  17,740,719   2401.5 Mb   5.1 %  13.23
    chr6           170.8 Mb  16,544,998   2240.0 Mb   4.7 %  13.11
    chr7           159.3 Mb  15,806,289   2144.2 Mb   4.5 %  13.46
    chr8           145.1 Mb  14,174,411   1923.1 Mb   4.1 %  13.25
    chr9           138.4 Mb  12,030,410   1639.8 Mb   3.5 %  11.85
    chr10          133.8 Mb  13,489,393   1822.2 Mb   3.9 %  13.62
    chr11          135.1 Mb  13,204,566   1794.4 Mb   3.8 %  13.28
    chr12          133.3 Mb  13,022,071   1767.1 Mb   3.7 %  13.26
    chr13          114.4 Mb  10,215,665   1381.9 Mb   2.9 %  12.08
    chr14          107.0 Mb  8,628,259    1167.9 Mb   2.5 %  10.91
    chr15          102.0 Mb  8,018,740    1092.0 Mb   2.3 %  10.71
    chr16          90.3 Mb   8,643,329    1182.0 Mb   2.5 %  13.08
    chr17          83.3 Mb   7,788,629    1055.1 Mb   2.2 %  12.67
    chr18          80.4 Mb   7,898,796    1071.2 Mb   2.3 %  13.33
    chr19          58.6 Mb   5,306,793    725.5 Mb    1.5 %  12.38
    chr20          64.4 Mb   6,554,184    886.3 Mb    1.9 %  13.75
    chr21          46.7 Mb   4,505,754    610.7 Mb    1.3 %  13.08
    chr22          50.8 Mb   3,634,924    498.7 Mb    1.1 %  9.81
    chrX           156.0 Mb  15,120,488   2061.5 Mb   4.4 %  13.21
    chrY           57.2 Mb   677,995      38.0 Mb     0.1 %  0.66
    chrM           16.6 Kb   119,173      17.8 Mb     0.0 %  1075.42
    Soft clip                             1509.6 Mb   3.2 %
    Unmapped                 2,063,908    309.6 Mb    0.7 %
    Duplicate                31,525,370   4741.9 Mb   10.1%
    Other contigs            9,593,918    1252.8 Mb   2.7 %
    TOTAL                    334,374,156  47178.5 Mb  100.0%
    
**Using with FASTQ**

  This works both on SAM/BAM and on FASTQ. For FASTQ, the mapping to contigs is as reported by the Genozip Aligner. The Genozip Aligner maps reads for compression purposes and does not attempt to map them according to the biological truth. However, usually the large majority of reads are in fact mapped to their correct position, so this can give a reasonable approximation of coverage of the data directly from FASTQ without needing to map it to BAM. 
