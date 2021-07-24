.. _23andMe2vcf:

Converting a 23andMe *Raw Genetic File* to VCF
==============================================

Data Types: 23andMe

`23andMe <https://www.23andme.com>`_ customers can download their raw genetic data, following these `instructions <https://customercare.23andme.com/hc/en-us/articles/212196868-Accessing-Your-Raw-Genetic-Data>`_.

However, this data comes in propietary 23andMe format. 

The 23andMe file is called something like `genome_John_Doe_v3_Full_20190101201010.zip` (the exact file name format may vary). 

Here, we explain how to convert the file to the standard VCF format.

Step 1: Download a reference file - any version of hg19 or GRCh37 will do, for example this one: `hs37d5.fa.gz <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz>`_. This file is quite large: appoximately 900MB.

Step 2: Create a Genozip reference file: This takes about 10 minutes to run:

::

    genozip --make-reference hs37d5.fa.gz

Step 3: Compress your 23andMe file with Genozip: 

::

    genozip genome_John_Doe_v3_Full_20190101201010.zip

Step 4: Convert the file to VCF: 

::

    genocat --reference hs37d5.ref.genozip --vcf genome_John_Doe_v3_Full_20190101201010.genozip --output mydata.vcf

Note: INDEL genotypes ('DD' 'DI' 'II') as well as uncalled sites ('--') are discarded
