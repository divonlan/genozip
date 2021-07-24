.. _multifasta2phylip.rst:

Converting MultiFASTA to Phylip and back
========================================

Data Types: FASTA, Phylip

Converting a MultiFASTA file to a Phylip file - all sequences must be the same length:

::

    genozip mydata.fa.gz
    genocat mydata.fa.genozip --phylip --output mydata.phy

Converting a Phylip file to a MultiFASTA file:

::

    genozip mydata.phy
    genocat mydata.phy.genozip --fasta --output mydata.fa.gz
    
Note: the input files can be plain files, or compressed with .gz .bz2 or .xz. The output files may be plain files or .gz.

Note: if the Phylip input file doesn't have a .phy (or .phy.gz / .phy.bz2 / .phy.xz) file name extenion, you can tell genozip that this is a Phylip file with ``--input phy``.