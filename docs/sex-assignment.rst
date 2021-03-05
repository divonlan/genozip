.. _sex:

Sex assignment
==============

.. toctree::
   :hidden:

   sam <sex-assignment-alg-sam>
   fastq <sex-assignment-alg-fastq>

Data Types: SAM, BAM and FASTQ

**Usage**

``genocat --show-sex my-file.bam.genozip`` 

``genocat --show-sex my-file.fq.genozip`` 

**Description**

Determines the genetic *Sex* of the individual from which the data originated.

**Output**

Output is space-seperated if sent to a terminal, and tab-seperated if redirected to a file or pipe.

::

  $ genocat *.genozip --show-sex

  Sex     File                        DP_1    DP_X    DP_Y    1/X   X/Y
  Male    sample01.10x.bam.genozip    13.188  6.921   5.887   1.9   1.2
  Male    sample02.10x.bam.genozip    14.644  7.662   7.432   1.9   1.0
  Male    sample03.10x.bam.genozip    10.993  5.614   6.135   2.0   0.9
  Female  sample04.10x.bam.genozip    4.274   4.238   0.215   1.0   19.7
  Male    sample05.10x.bam.genozip    12.497  6.565   4.320   1.9   1.5
  Female  sample06.10x.bam.genozip    10.728  10.866  0.534   1.0   20.4


**Classifier algorithm**
  | See :doc:`sex-assignment-alg-sam` and :doc:`sex-assignment-alg-fastq`

**Advantages**
  | • This method works fast, because just the relevant fields are read from disk - a small subset of the file.
  |
  | • Since this algorithm is based on counting bases rather than counting reads, it would work just as well with data that has highly variable read lengths, as common, for example, in long-read technologies.
  |
  | • Since the algorithm uses both X/Y and Autosome/X, it can detect a XXY.

**Limitations**
  | • This feature has been tested on human data. It may or may not work for other species.
  |
  | • Rare sexes such as XXXY or XYY are not detected.
  |
  | • In case of a bound genozip file (i.e. one created with eg ``genozip file1.sam file2.sam file3.sam --output my-bound.sam.genozip``) calculation will be done on the entire file leading non-sensical results.
  |
  | • SAM/BAM lines with an undefined RNAME or CIGAR are ignored, therefore this will not work for unaligned SAM/BAM files.
  |
  | • This feature should never be used for clinical diagnostics.

**Additional limitations when used with FASTQ**
  | • Sex assignment of a FASTQ is a an experimental feature, and the results should be taken with a grain of salt.
  |
  | • Assignment of reads to contigs is based on the Genozip Aligner which is designed to be very fast at the expense of accuracy. This feature works because the Male vs. Female signal is usually stronger than the Genozip Aligner inaccuracies, however when running on FASTQ, the algorithm will be more conservative than when running on SAM/BAM and will report "Unassigned" in some cases where it would confidently call a Sex when running on SAM/BAM data.
  |
  | • This has been tested on fastq.genozip files which consist of human data compressed with the GRCh38 reference genome. It yet unknown if the algorithm would work with other species data or other reference genomes.  
