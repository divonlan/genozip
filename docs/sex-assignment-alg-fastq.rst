Sex Classifier algorithm (FASTQ)
================================

  | 1) Calculate the per-contig Depth as described in :doc:`coverage`.
  |
  | 2) Calculate X_Depth / Y_Depth and Autosome_Depth / X_Depth. 
  |
  | 3) Multiply Autosome_Depth / X_Depth by a *correction factor* of 1.333 - this is to correct for observed Genozip Aligner biases in favor of X.
  |
  | 4) Decision matrix (for SAM / BAM): 

    =========================== ================================= =============
    **Assigned Sex**            **AS / X**                         **X / Y**
    *Male*                      > 1.75 (single-X)                 < 9 (Has Y)
    *Female*                    < 1.3 (double-X or almsot)   > 5 (Not enough Y for Male)
    *Male-XXY or XXYY*          < 1.1 (double-X)                  < 1.8 (Has Y, possibly similar ratio Y to X) 
    *Unassigned*                All other combinations
    =========================== ================================= =============

  | *Definitions*:
  | • *Depth* is defined here :doc:`coverage`.
  |
  | • *Autosome_Depth* means combined coverage of all *autosome contigs* divided by combined length of all *autosome contigs*.
  |
  | • *Autosome contigs* are all contigs excluding X, Y, MT and excluding non-primary contigs like ``chr22_KI270731v1_random``.
  |
  | • *Chromosome X* is the contig named "X", "chrX" or "ChrX", and similarly for Y. For MT, contig names based on "M" and "MT" are accepted.
  |
  | • The *correction factor* of 1.333 is set based on empirical tests on human data with the GRCh38 reference genome. It is not yet known if this factor will remain the same for other species or other reference genomes.
  

