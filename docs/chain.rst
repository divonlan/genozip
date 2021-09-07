.. _chain:

Chain files
===========

.. include:: dvcf-see-also.rst

**Obtaining a chain file**

Here is a non-exhaustive list of chain files (also called LiftOver files). The creators of these chain files are not associated with Genozip in any way, and Genozip does not endorse, recommend or warrant the correctness of any particular file.

`From Ensembl <https://ftp.ensembl.org/pub/assembly_mapping/>`_

`From UCSC <https://hgdownload.soe.ucsc.edu/downloads.html>`_

**Chain file format**

The file format of chain files can be found here: `UCSC chain file format <https://genome.ucsc.edu/goldenPath/help/chain.html>`_. 

**Chain files and Genozip**

Chain files are compressable with Genozip, and indeed, need to be compressed before used with ``genozip --chain`` to :ref:`lift a VCF file<dvcf>`. To compress a chain file, two reference files are required - the reference genome in Primary coordinates, and the reference genome in Luft coordintes, for example:

::

    genozip GRCh37_to_GRCh38.chain --reference hs37d5.ref.genozip --reference GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip

**Viewing chain file alignments**

To view the alignment in the chain, use ``--show-chain``:

.. _show_chain:
    
::

    > genocat --show-chain GRCh37_to_GRCh38.chain.genozip

    ##fileformat=GENOZIP-CHAIN
    #ALN_I  PRIM_CONTIG     PRIM_START      PRIM_END        LUFT_CONTIG     LUFT_START      LUFT_ENDS       XSTRAND ALN_OVERLAP
    1       1       10001   177417  chr1    10001   177417  -
    2       1       227418  267719  chr1    257667  297968  -
    3       1       317720  471368  chr1    347969  501617  X
    4       1       521369  1566075 chr1    585989  1630695 -
    5       1       1566076 1569784 chr1    1630697 1634405 -
    6       1       1569785 1570918 chr1    1634409 1635542 -
    7       1       1570919 1570922 chr1    1635547 1635550 -
    8       1       1570923 1574299 chr1    1635561 1638937 -
    9       1       1574300 1583669 chr1    1638939 1648308 -

Column 1 is the sequential number of this alignment.

Columns 2,3,4 are coordinates in the Primary reference genome and columns 5,6,7 are coordinates in the Luft reference genome. Note that these are 1-based coordinates (as in VCF), whereas the Chain file format has 0-based coordinates. 

Column 8 contains an ``X`` if the range in the Luft reference is reverse complimented.

Column 9 contains the ALN_I of an alignment with an overlapping Luft range: This is for cases in which multiple alignments contain overlapping Luft ranges. If this alignment's Luft range overlaps with the Luft range of more than one other alignment, the column will contain one of them.

**Viewing chain file contig information**

::

    genocat GRCh37_to_GRCh38.chain.genozip --show-chain-contigs

    PRIMARY chain file contigs that also exist in the reference file:
    PRIMARY 1 length=249250621
    PRIMARY 2 length=243199373
    PRIMARY 3 length=198022430
    PRIMARY 4 length=191154276
    PRIMARY 5 length=180915260

    LUFT chain file contigs that also exist in the reference file:
    LUFT chr1 length=248956422
    LUFT chr2 length=242193529
    LUFT chr3 length=198295559
    LUFT chr4 length=190214555
    LUFT chr5 length=181538259

Note: All contigs of the respective reference files are shown, even if not referred to in the chain data.
