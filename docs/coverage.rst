Coverage and Depth
==================

Data Types: SAM, BAM and FASTQ

``genocat --show-coverage my-file.bam.genozip`` displays the per-contig Coverage and Depth.

``genocat --show-coverage-chrom my-file.bam.genozip`` does the same, but just for contigs that are chromosomes.

Fine details: 

  *Coverage* is the number of bases mapped to a contig. 
  
  *Depth* is the Coverage divided by the Length (LN) of the contig.
  
  *Length* is defined in the header of the SAM/BAM, or in the reference file used to compress it, or in the case of a header-less SAM file that was compressed without a reference - the highest coordinate of the contig that appears in the file, rounded up to the next 1MB.

  *Contigs that are chromosomes* are contigs with a name of up to 5 characters. For example ``chr22`` is, but ``chr22_KI270731v1_random`` is not.

  In SAM and BAM (but not FASTQ) we exclude:
  
  #. We exclude bases that have an ``S`` CIGAR ("soft clipping") (see 1.4.6 `here <https://samtools.github.io/hts-specs/SAMv1.pdf>`_).
  
  #. We exclude reads that are not mapped (FLAG=0x4), failed filters (FLAG=0x200) are secondary (FLAG=0x100) or are duplicate (FLAG=0x400)

**Output**
  | A tab-seperated table, example:
    
::

    $ genocat myfile.bam.genozip --show-coverage-chrom

    Contig  LN      Coverage        Depth
    chr1    248956422       2667178287       10.71
    chr2    242193529       2641166618       10.91
    chr3    198295559       2191671750       11.05
    chr4    190214555       2105061023       11.07
    chr5    181538259       1946975026       10.72
    chr6    170805979       1815243563       10.63
    chr7    159345973       1733006366       10.88
    chr8    145138636       1568148827       10.80
    chr9    138394717       1318887120        9.53
    chr10   133797422       1481588178       11.07
    chr11   135086622       1465298907       10.85
    chr12   133275309       1443045061       10.83
    chr13   114364328       1130334663        9.88
    chr14   107043718       953895279         8.91
    chr15   101991189       881547642         8.64
    chr16   90338345        975530338        10.80
    chr17   83257441        870265768        10.45
    chr18   80373285        884126928        11.00
    chr19   58617616        601881274        10.27
    chr20   64444167        737073875        11.44
    chr21   46709983        483872201        10.36
    chr22   50818468        414579819         8.16
    chrX    156040895       882316415         5.65
    chrY    57227415        215579644         3.77
    chrM    16569   16851259        1017.04

**Using with FASTQ**

    The data reported is NOT the true coverages and depths when running on FASTQ. Rather, it reports data derived from the assignment of reads to contigs by the Genozip Aligner. The Genozip Aligner is designed to be very fast at the expense of inaccuracies, and does not make the claim of mapping reads consistent with the biological truth. 
