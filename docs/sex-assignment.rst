.. _sex:

Sex assignment
==============

.. toctree::
   :hidden:

   sam <sex-assignment-alg-sam>
   fastq <sex-assignment-alg-fastq>

Data Types: SAM, BAM and FASTQ

``genocat --show-sex my-file.bam.genozip`` determines the genetic *Sex* of the individual from which the data originated.

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
