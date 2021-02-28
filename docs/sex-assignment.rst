Sex assignment
==============

Data Types: SAM and BAM

``genocat --show-sex my-file.bam.genozip`` determines *Sex* to the data.

**Output**
  | A tab-seperated table, example:

::

  $ genocat *.genozip --show-sex

  File                  1_Depth   X_Depth  Y_Depth  1/X   X/Y  Sex
  sample01.bam.genozip  10.71343  5.65439  3.76707  1.9   1.5  Male
  sample02.bam.genozip  14.99745  7.63580  5.30143  2.0   1.4  Male
  sample03.bam.genozip  14.55365  7.29276  5.71087  2.0   1.3  Male
  sample04.bam.genozip  17.32741  9.14707  8.88780  1.9   1.0  Male
  sample05.bam.genozip  14.56035  7.61892  6.32542  1.9   1.2  Male
  sample06.bam.genozip  16.48654  8.53986  8.08142  1.9   1.1  Male
  sample07.bam.genozip  12.17690  6.20695  6.63175  2.0   0.9  Male
  sample08.bam.genozip   4.68230  4.64500  0.22333  1.0  20.8  Female
  sample09.bam.genozip  14.10800  7.31946  4.66873  1.9   1.6  Male
  sample10.bam.genozip  12.11650 12.26979  0.56171  1.0  21.8  Female
  sample11.bam.genozip  16.48835 16.00353  0.86917  1.0  18.4  Female
  sample12.bam.genozip  17.15595  8.95983  7.67593  1.9   1.2  Male
  sample13.bam.genozip  17.39297  8.76055  9.27279  2.0   0.9  Male


**Classifier algorithm**
  | 1) Calculate the per-contig Depth as described in :doc:`coverage`.
  |
  | 2) Calculate the ratio 1_Depth / X_Depth and X_Depth / Y_Depth.
  |
  | 3) Decision matrix: 

    =========================== ================================= =============
    **Assigned Sex**            **1 / X**                         **X / Y**
    *Male*                      > 1.75 (single-X)                 < 9 (Has Y)
    *Female*                    < 1.1 (double-X)                  > 5 (Not enough Y for Male)
    *Female*                    < 1.3 (not quite double-X)        > 9 (No Y)
    *Male-XXY*                  < 1.1 (double-X)                  ∈ (1.8,5) (XXY)
    *Male-XXY or XY/XXY mosaic* ∈ (1.1, 1.3) (not quite double-X) ∈ (1.8,5) (XXY)
    *Male-XXY or XXYY*          < 1.1 (double-X)                  < 1.8 (Has Y, possibly similar ratio Y to X) 
    *Unassigned*                All other combinations
    =========================== ================================= =============

  | *Definitions*:
  | • Chromosome X is the contig named "X", "chrX" or "ChrX". Likewise for Y and 1.
  |
  | • Depth is the total number of bases mapped to the chromosome which are marked by the CIGAR Op ``M``, ``X``, ``=`` or ``I`` (see 1.4.6 `here <https://samtools.github.io/hts-specs/SAMv1.pdf>`_), divided by the length (LN) of the chromosome as defined in the SAM/BAM header.
  
**Advantages**
  | • This method works extremely fast, because just the RNAME and CIGAR data is read from disk - a small subset of the file.
  |
  | • Since this algorithm is based on counting bases rather than counting reads, it would work just as well with data that has highly variable read lengths, as common, for example, in long-read technologies.
  |
  | • Since the algorithm uses both X/Y and 1/X, it can detect a XXY.

**Limitations**
  | • This feature has been tested on human data. It may or may not work for other species.
  |
  | • Rare sexes such as XXXY or XYY are not detected.
  |
  | • In case of a bound genozip file (i.e. one created with eg ``genozip file1.sam file2.sam file3.sam --output my-bound.sam.genozip``) calculation will be done on the entire file leading non-sensical results.
  |
  | • Lines with an undefined RNAME or CIGAR are ignored, therefore this will not work for unaligned SAM/BAM files.
  |
  | • This feature should never be used for clinical diagnostics.
