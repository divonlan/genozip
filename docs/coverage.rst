.. _coverage:

Coverage and Depth
==================

Data Types: SAM, BAM and FASTQ

``genocat --show-coverage=all my-file.bam.genozip`` displays the per-contig Number-of-reads, Number-of-bases and Depth.

``genocat --show-coverage my-file.bam.genozip`` does the same, but just for contigs that are chromosomes.

Fine details: 

  *Number-of-bases* is the number of bases mapped to a contig, excluding bases with a ``S`` (Soft Clip) CIGAR (see 1.4.6 `here <https://samtools.github.io/hts-specs/SAMv1.pdf>`_).

  *Number-of-reads* is the number of reads mapped to a contig.
  
  *Depth* is the Coverage divided by the Length (LN) of the contig.
  
  *Length* is defined in the header of the SAM/BAM, or in the reference file used to compress it, or in the case of a header-less SAM file that was compressed without a reference - the highest coordinate of the contig that appears in the file, rounded up to the next 1MB.

  *Contigs that are chromosomes* are contigs with a name of up to 5 characters. For example ``chr22`` is, but ``chr22_KI270731v1_random`` is not.

  *Number-of-bases* and *Number-of-reads* excludes unmapped reads, and in SAM/BAM also excludes reads with a FLAG indicating Duplicate, Secondary and Failed Filters. These are reported separately.
  
**Output**
    
::

    $ genocat.exe --show-coverage my-sample.bam.genozip
    Contig         LN        Reads        Bases       Depth
    1              249.3 Mb  61,261,021   9000.6 Mb   36.11 
    2              243.2 Mb  63,534,042   9314.4 Mb   38.30
    3              198.0 Mb  50,851,120   7479.4 Mb   37.77
    4              191.2 Mb  48,167,951   7082.8 Mb   37.05
    5              180.9 Mb  46,524,643   6842.5 Mb   37.82
    6              171.1 Mb  43,689,854   6423.2 Mb   37.54
    7              159.1 Mb  41,213,705   6050.3 Mb   38.02
    8              146.4 Mb  38,109,245   5602.2 Mb   38.28 
    9              141.2 Mb  31,779,723   4668.8 Mb   33.06
    10             135.5 Mb  34,988,073   5137.7 Mb   37.91
    11             135.0 Mb  35,262,312   5183.1 Mb   38.39
    12             133.9 Mb  34,551,733   5075.2 Mb   37.92
    13             115.2 Mb  24,543,420   3608.4 Mb   31.33
    14             107.3 Mb  23,546,204   3459.0 Mb   32.22
    15             102.5 Mb  22,632,125   3323.9 Mb   32.42
    16             90.4 Mb   23,703,488   3471.9 Mb   38.43
    17             81.2 Mb   22,169,323   3249.4 Mb   40.02
    18             78.1 Mb   19,555,374   2874.6 Mb   36.82
    19             59.1 Mb   16,392,666   2395.8 Mb   40.52
    20             63.0 Mb   16,544,349   2428.3 Mb   38.53
    21             48.1 Mb   9,823,711    1440.9 Mb   29.94
    22             51.3 Mb   10,247,743   1501.8 Mb   29.27
    X              155.3 Mb  20,446,847   2996.6 Mb   19.30
    Y              59.4 Mb   3,984,478    582.7 Mb    9.81
    MT             16.6 Kb   106,209      15.5 Mb     938.41
    Soft clip                             1185.0 Mb
    Unmapped                 17,071,958   2433.6 Mb
    Duplicate                71,068,615   10556.8 Mb
    Other contigs            29,259,649   4138.6 Mb

**Using with FASTQ**

    The data reported is NOT the true coverages and depths when running on FASTQ. Rather, it reports data derived from the assignment of reads to contigs by the Genozip Aligner. The Genozip Aligner is designed to be very fast at the expense of inaccuracies, and does not make the claim of mapping reads consistent with the biological truth. 
