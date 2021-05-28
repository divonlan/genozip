.. _coverage:

Coverage and Depth
==================

**EXPERIMENTAL FEATURE**

Data Types: SAM, BAM and FASTQ

**Usage**

``genocat --show-coverage my-file.bam.genozip``

``genocat --show-coverage=all my-file.bam.genozip`` 

``genocat --show-coverage=one *.bam.genozip`` 

``genocat --show-coverage my-file.fq.genozip`` 

**Description**

Displays the per-chromosome reads, bases, % of total bases that are on this chromosome and Depth. With the argument =all, shows all contigs, not just primary chromosomes. With the argument =one, shows just a one number summary of coverage (the same number that appears as the *Depth* on the *All contigs* line).

**Output**
    
Output is space-seperated if sent to a terminal, and tab-seperated if redirected to a file or pipe.

::

    $ genocat.exe --coverage my-sample.bam.genozip
    
    --coverage for: my-sample.bam.genozip
    Contig         LN        Reads        Bases       Bases   Depth
    1              249.3 Mb  61,215,549   9.0 Gb      7.1 %   36.10
    2              243.2 Mb  63,342,488   9.3 Gb      7.3 %   38.26
    3              198.0 Mb  50,819,357   7.5 Gb      5.9 %   37.76
    4              191.2 Mb  48,131,726   7.1 Gb      5.6 %   37.04
    5              180.9 Mb  46,495,225   6.8 Gb      5.4 %   37.81
    6              171.1 Mb  43,659,445   6.4 Gb      5.0 %   37.53
    7              159.1 Mb  41,171,496   6.0 Gb      4.7 %   38.01
    8              146.4 Mb  38,081,676   5.6 Gb      4.4 %   38.27
    9              141.2 Mb  31,755,198   4.7 Gb      3.7 %   33.05
    10             135.5 Mb  34,954,121   5.1 Gb      4.0 %   37.90
    11             135.0 Mb  35,234,129   5.2 Gb      4.1 %   38.38
    12             133.9 Mb  34,523,219   5.1 Gb      4.0 %   37.91
    13             115.2 Mb  24,526,664   3.6 Gb      2.8 %   31.32
    14             107.3 Mb  23,525,557   3.5 Gb      2.7 %   32.21
    15             102.5 Mb  22,613,925   3.3 Gb      2.6 %   32.41
    16             90.4 Mb   23,675,707   3.5 Gb      2.7 %   38.41
    17             81.2 Mb   22,148,439   3.2 Gb      2.5 %   40.01
    18             78.1 Mb   19,542,517   2.9 Gb      2.3 %   36.81
    19             59.1 Mb   16,374,766   2.4 Gb      1.9 %   40.50
    20             63.0 Mb   16,530,922   2.4 Gb      1.9 %   38.52
    21             48.1 Mb   9,811,790    1.4 Gb      1.1 %   29.93
    22             51.3 Mb   10,238,901   1.5 Gb      1.2 %   29.26
    X              155.3 Mb  20,421,365   3.0 Gb      2.3 %   19.29
    Y              59.4 Mb   3,971,007    582.0 Mb    0.5 %    9.80
    MT             16.6 Kb   105,906      15.5 Mb     0.0 %  937.41
    Other contigs            28,989,258   4.1 Gb      3.2 %
    -----
    All contigs              771,860,353  113.3 Gb    88.8%   36.60
    Soft clip                             1.2 Gb      0.9 %
    Unmapped                 17,071,958   2.4 Gb      1.9 %
    Supplementary            1,028,655    52.2 Mb     0.0 %
    Duplicate                71,068,615   10.6 Gb     8.3 %
    TOTAL                    861,029,581  127.5 Gb    100.0%

**Details** 

  *Length* is defined in the header of the SAM/BAM, or in the reference file used to compress it, or in the case of a header-less SAM file that was compressed without a reference - the highest coordinate of the contig that appears in the file.

  *Contig* includes all contigs present in the file if run with ``--coverage=all``, or if not, just contigs with a name of up to 5 characters, and all other contigs are summed up in the *Other contings* line. The *All contigs* line sums up the values for all contigs present in the file.

  *Unmapped*, *Failed filters*, *Duplicate*, *Secondary*, *Supplementary* (together: *Excluded reads*) and reads with the appropriate SAM FLAG set (i.e. 0x4, 0x200, 0x400, 0x100, 0x800 respectively). Reads missing an RNAME are also counted as unmapped irrespective of their FLAG. Reads can count in only one category even if multiple flags are set, with this order of precedence.

  *Soft clip* are bases with an ``S`` (Soft Clip) CIGAR (see 1.4.6 `here <https://samtools.github.io/hts-specs/SAMv1.pdf>`_), but not counting those on *Excluded reads*.

  The *Reads* and *Bases* reported for each *Contig* exclude *Excluded reads* and the *Bases* also exclude *Soft clip* bases.
  
  *Depth* is *Bases* divided by the Length (LN) of the contig.
  
  
**Using with FASTQ**

  This works both on SAM/BAM and on FASTQ. For FASTQ, the mapping to contigs is as reported by the Genozip Aligner. The Genozip Aligner maps reads for compression purposes and does not attempt to map them according to the biological truth. However, usually the large majority of reads are in fact mapped to their correct position, so this can give a reasonable approximation of coverage of the data directly from FASTQ without needing to map it to BAM. 
