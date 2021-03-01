Sex Classifier algorithm (SAM/BAM)
==================================

  | 1) Calculate the per-contig Depth as described in :doc:`coverage`.
  |
  | 2) Calculate X_Depth / Y_Depth and chr1_Depth / X_Depth. 
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
  | • *Depth* is defined here :doc:`coverage`.
  |
  | • *Chromosome X* is the contig named "X", "chrX" or "ChrX", and similarly for Y and chr1. 
  

