.. _dvcf:

Dual-coordinates VCF files
==========================

See also:

    | :ref:`Rendering a DVCF <dvcf-rendering>`
    |
    | :ref:`Dual-coordinates VCF Specification <dvcf-spec>`
    |
    | :ref:`Chain files <dvcf-chain-files>`

.. toctree::
    :hidden:
 
    Rendering a DVCF <dvcf-rendering>
    Dual-coordinates VCF Specification <dvcf-spec>
    Chain files <dvcf-chain-files>

Dual coordinates VCF files (or DVCFs), are VCF files that contain coordinates in two coordinate systems concurrently, for example, GRCh37 and GRCh38. A DVCF can be *rendered* in either *Primary* coordinates or *Luft* coordinates (*Luft* being our made-up past particile of *Lift*). 

Both *renditions* are standard-compliant VCF files, and both contain precisely the same information, with the only difference being how the information is presented.

Crucially, since both renditions contain all the information needed for rendering the file in both coordinate systems, the file in either rendition can be processed through any bioinformatics tool or pipeline, and still maintain its dual coordinates. In other words, an analytics pipeline may now have some steps that operate on the file in one coordinate system, and other steps operate in other coordinate system, and the same data can flow through these steps seemlessly.

Let's walk through an example:

Consider this small VCF file, which is in the coordinates defined by GRCh37, a popular Homo Sapiens reference:

::

    ##fileformat=VCFv4.2
    ##reference=file:///references/grch37/reference.bin
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##contig=<ID=1,length=249250621>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1	Person2
    1	10285	1	T	C	4.4	PASS	AC=3;AN=4	GT:AD:AF:PL	0/1:31,18:0.367:37,0,46	1/1
    1	329162	2	A	T	4.6	PASS	AC=3;AN=4	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1
    1	366043	3	CA	A	100	PASS	.	GT	1|0	0|0

We shall now convert this file to a Dual-coordinates VCF, containing both GRCh37 coordinates, as well as coordinates from another popular reference, GRCh38. This process is called *lifting*. We shall call GRCh37 the *Primary* coordinates, and GRCh38 the *Luft* coordinates.

Steps:

1. Prepare a reference file in Luft coordinates in genozip format. The following results in the file named ``GRCh38.ref.genozip`` :

::

    genozip --make-reference GRCh38.fa 


2. Prepare a chain file in genozip format. Genozip uses chain files in the popular `UCSC chain file format <https://genome.ucsc.edu/goldenPath/help/chain.html>`_. Here are some :ref:`links to available chain files<dvcf-chain-files>`. The following command line generates the file ``GRCh37_to_GRCh38.chain.genozip`` :

::

    genozip GRCh37_to_GRCh38.chain.gz --reference GRCh38.ref.genozip

Note that the first two steps are preparation steps that need to be executed only once. Now that we have the chain file, we can convert any number of VCF files to DVCF:

:: 

    genozip mydata.vcf.gz --chain GRCh37_to_GRCh38.chain.genozip

Note that as usual with ``genozip``, .vcf .vcf.gz .vcf.bz2 and .vcf.xz are all accepted. Genozip does not require the VCF file to be indexed.

Now, we have a dual-coordinate VCF file: ``mydata.vcf.genozip``. We can render it in either coordinate:

Rendering in Primary coordinates (GRCh37), followed by Luft coordinates (GRCh38):

::

    > genocat mydata.vcf.genozip | tee mydata.37.vcf

    ##fileformat=VCFv4.2
    ##dual_coordinates=PRIMARY
    ##chain=file:///data/GRCh37_to_GRCh38.chain.genozip
    ##reference=file:///references/grch37/reference.bin
    ##luft_reference=file:///data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles",RendAlg=R>
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed",RendAlg=A_1>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",RendAlg=GT>
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes",RendAlg=G>
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele",RendAlg=A_AN>
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes",RendAlg=NONE>
    ##INFO=<ID=LUFT,Number=4,Type=String,Description="Info for rendering variant in LUFT coords",RendAlg=NONE>
    ##INFO=<ID=PRIM,Number=4,Type=String,Description="Info for rendering variant in PRIMARY coords",RendAlg=NONE>
    ##INFO=<ID=Lrej,Number=1,Type=String,Description="Reason variant was rejected for LUFT coords",RendAlg=NONE>
    ##INFO=<ID=Prej,Number=1,Type=String,Description="Reason variant was rejected for PRIMARY coords",RendAlg=NONE>
    ##contig=<ID=1,length=249250621>
    ##luft_contig=<ID=chr1,length=248956422>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1	Person2
    1	10285	1	T	C	4.4	PASS	AC=3;AN=4;LUFT=chr1,10285,T,-	GT:AD:AF:PL	0/1:31,18:0.367:37,0,46	1/1
    1	329162	3	A	T	4.6	PASS	AC=3;AN=4;LUFT=chr1,248466248,T,-	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1
    1	366043	4	CA	A	100	PASS	Lrej=RefTooLong	GT	1|0	0|0

    > genocat mydata.vcf.genozip --luft | tee mydata.38.vcf

    ##fileformat=VCFv4.2
    ##primary_only=1	366043	4	CA	A	100	PASS	Lrej=RefTooLong	GT	1|0	0|0
    ##dual_coordinates=LUFT
    ##chain=file:///data/GRCh37_to_GRCh38.chain.genozip
    ##reference=file:///data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##primary_reference=file:///references/grch37/reference.bin
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles",RendAlg=R>
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed",RendAlg=A_1>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",RendAlg=GT>
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes",RendAlg=G>
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele",RendAlg=A_AN>
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes",RendAlg=NONE>
    ##INFO=<ID=LUFT,Number=4,Type=String,Description="Info for rendering variant in LUFT coords",RendAlg=NONE>
    ##INFO=<ID=PRIM,Number=4,Type=String,Description="Info for rendering variant in PRIMARY coords",RendAlg=NONE>
    ##INFO=<ID=Lrej,Number=1,Type=String,Description="Reason variant was rejected for LUFT coords",RendAlg=NONE>
    ##INFO=<ID=Prej,Number=1,Type=String,Description="Reason variant was rejected for PRIMARY coords",RendAlg=NONE>
    ##primary_contig=<ID=1,length=249250621>
    ##contig=<ID=chr1,length=248956422>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1	Person2
    chr1	10285	1	T	C	4.4	PASS	AC=3;AN=4;PRIM=1,10285,T,-	GT:AD:AF:PL	0/1:31,18:0.367:37,0,46	1/1
    chr1	248466248	3	T	A	4.6	PASS	AC=1;AN=4;PRIM=1,329162,A,-	GT:AD:AF:PL	1/0:9,28:0.7:0,0,36	0/0

So what do we see here?

The first thing to appreciate, is that both ``mydata.37.vcf`` and ``mydata.38.vcf`` contain precisely the same information. In other words, you can ``genozip`` either file, and then render it again in both Primary and Luft coordinates.

Take a look at the first file - the Primary rendition - at the variant with ID=3. The INFO/LUFT field contains the information of the Luft coordinate, which, in this case, is mapped to the Luft reference at CHROM=chr1 POS=248466248, and the REF (the third value in INFO/LUFT) is T - which would be change from the current (GRCh37) REF that is A, and identical to the current ALT which is a T. This is common situation which we call a REF⇆ALT switch.

Now take a look at the second file - the Luft rendition - at ID=3. Indeed, the REF and ALT have switched - and also the INFO/AC, FORMAT/GT, FORMAT/AD, FORMAT/AF and FORMAT/PL fields have changed too, to reflect this REF⇆ALT switch. In addition, the INFO/PRIM field contains the Primary coordinates of this variant.

Now back to the Primary rendition - take a look the variant with ID=4. It has an INFO/Lrej field indicating that Genozip cannot render it in Luft coordinates and hence it is a Primary-only variant. The reason given is ``RefTooLong`` - in other words, Genozip can't cross-render this variant because its REF is ``CA`` and Genozip can only handle *lifting* of single-base REF fields. To understand more how Genozip cross-renders each field, see :ref:`dvcf-rendering`.

In the Luft rendition, this Primary-only variant, ID=4, is indeed missing in the VCF data lines, which now contain only two variants. However, notice that it still survives as a meta-information line with the key ``##primary_only``.

Now, let's take our Luft rendition VCF ``mydata.38.vcf`` and merge it (using ``bcftools merge`` in this example) with another VCF file, ``person3.vcf``, which is in GRCh38 coordinates and has a single sample, Person3. Note that ``person3.vcf`` is *not* a DVCF. ``person3.vcf`` contains 2 variants on chr1: ID=3 (also present in ``mydata.38.vcf``), and ID=2 (not present in ``mydata.38.vcf``):

::

    > cat person3.vcf
    
    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person3
    chr1	20197381	2	G	C	100	PASS	AC=1;AN=2	GT:AD	1/0:20,15
    chr1	248466248	3	T	A	4.6	PASS	AC=0;AN=2	GT:AD:GP	0/0:30,3:37.9,14.9,0.14
    
Notice that Person3 has the field FORMAT/GP in variant ID=3 that doesn't exist in ``mydata.38.vcf``. 

The merged file will be a DVCF file looking like this:

::

    > bgzip person3.vcf 
    > bgzip mydata.38.vcf
    > bcftools index mydata.38.vcf.gz 
    > bcftools index person3.vcf.gz
    > bcftools merge mydata.38.vcf.gz person3.vcf.gz | tee merged.38.vcf

    ##fileformat=VCFv4.2
    ##primary_only=1        366043  4       CA      A       100     PASS    Lrej=RefTooLong GT      1|0     0|0
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##dual_coordinates=LUFT
    ##chain=file:///data/GRCh37_to_GRCh38.chain.genozip
    ##reference=file:///data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##primary_reference=file:///references/grch37/reference.bin
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles",RendAlg=R>
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed",RendAlg=A_1>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",RendAlg=GT>
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes",RendAlg=G>
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele",RendAlg=A_AN>
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes",RendAlg=NONE>
    ##INFO=<ID=LUFT,Number=4,Type=String,Description="Info for rendering variant in LUFT coords",RendAlg=NONE>
    ##INFO=<ID=PRIM,Number=4,Type=String,Description="Info for rendering variant in PRIMARY coords",RendAlg=NONE>
    ##INFO=<ID=Lrej,Number=1,Type=String,Description="Reason variant was rejected for LUFT coords",RendAlg=NONE>
    ##INFO=<ID=Prej,Number=1,Type=String,Description="Reason variant was rejected for PRIMARY coords",RendAlg=NONE>
    ##primary_contig=<ID=1,length=249250621>
    ##contig=<ID=chr1,length=248956422>
    ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification">
    ##bcftools_mergeVersion=1.11+htslib-1.11
    ##bcftools_mergeCommand=merge mydata.38.vcf.gz person3.vcf.gz; Date=Wed May 26 19:32:48 2021
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Person1 Person2 Person3
    chr1    10285   1       T       C       4.4     PASS    PRIM=1,10285,T,-;AN=4;AC=3      GT:AD:AF:PL     0/1:31,18:0.367:37,0,46 1/1:.:.:.       ./.:.:.:.
    chr1    20197381        2       G       C       100     PASS    AN=2;AC=1       GT:AD   ./.:.   ./.:.   1/0:20,15
    chr1    248466248       3       T       A       4.6     PASS    PRIM=1,329162,A,-;AN=6;AC=1     GT:AD:AF:PL:GP  1/0:9,28:0.7:0,0,36:.   0/0:.:.:.:.     0/0:30,3:.:.:37.9,14.9,0.14

The merged file, as expected, is in Luft (i.e. GRCh38) coordinates, and the new variant, ID=2, has no DVCF INFO tag (i.e. no INFO/PRIM).

Now let us genozip the merged file and render it first in Luft coordinates:

::

    > genozip merged.38.vcf
    genozip merged.38.vcf : 0%
    FYI: the number of samples in variant CHROM=1 POS=366043 is 2, different than the VCF column header line which has 3 samples
    Done (0 seconds)

    > genocat merged.38.vcf.genozip --luft
    
    ##fileformat=VCFv4.2
    ##primary_only=1        366043  4       CA      A       100     PASS    Lrej=RefTooLong GT      1|0     0|0
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##dual_coordinates=LUFT
    ##chain=file:///data/GRCh37_to_GRCh38.chain.genozip
    ##reference=file:///data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##primary_reference=file:///references/grch37/reference.bin
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles",RendAlg=R>
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed",RendAlg=A_1>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",RendAlg=GT>
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes",RendAlg=G>
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele",RendAlg=A_AN>
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes",RendAlg=NONE>
    ##INFO=<ID=LUFT,Number=4,Type=String,Description="Info for rendering variant in LUFT coords",RendAlg=NONE>
    ##INFO=<ID=PRIM,Number=4,Type=String,Description="Info for rendering variant in PRIMARY coords",RendAlg=NONE>
    ##INFO=<ID=Lrej,Number=1,Type=String,Description="Reason variant was rejected for LUFT coords",RendAlg=NONE>
    ##INFO=<ID=Prej,Number=1,Type=String,Description="Reason variant was rejected for PRIMARY coords",RendAlg=NONE>
    ##primary_contig=<ID=1,length=249250621>
    ##contig=<ID=chr1,length=248956422>
    ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification",RendAlg=G>
    ##bcftools_mergeVersion=1.11+htslib-1.11
    ##bcftools_mergeCommand=merge mydata.38.vcf.gz person3.vcf.gz; Date=Wed May 26 19:32:48 2021
    ##genozip_command="genocat --luft newdata.38.vcf.genozip" 2021-05-26 19:36:01 ACDT
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Person1 Person2 Person3
    chr1    10285   1       T       C       4.4     PASS    AN=4;AC=3;PRIM=1,10285,T,-      GT:AD:AF:PL     0/1:31,18:0.367:37,0,46 1/1:.:.:.       ./.:.:.:.
    chr1    20197381        2       G       C       100     PASS    AN=2;AC=1;Prej=AddedVariant     GT:AD   ./.:.   ./.:.   1/0:20,15
    chr1    248466248       3       T       A       4.6     PASS    AN=6;AC=1;PRIM=1,329162,A,-     GT:AD:AF:PL:GP  1/0:9,28:0.7:0,0,36:.   0/0:.:.:.:.     0/0:30,3:.:.:37.9,14.9,0.14

Looking at the Luft rendition of the merged file, we notice that the new variant, ID=2, now contains a INFO/Prej field. Therefore, this is a Luft-only variant that cannot be expressed in Primary coordinates. This is always the case for variants that are added to an existing DVCF file, unless the new variants also originate from a DVCF file.

Also, note that when running ``genozip``, we got a warning about the number of samples of variant ID=4. This is our Primary-only variant, that didn't participate in the merge in Luft coordinates, and remains with two samples only. This is ok, as long as the DVCF file is the first file on the ``bctools merge`` command line, so that the missing sample aligns with the new sample Person3. It is perfectly acceptable to have missing samples in VCF.

Let us now take a look at the Primary rendition of the merged file:

::

    > genocat merged.38.vcf.genozip
    
    ##fileformat=VCFv4.2
    ##luft_only=chr1        20197381        2       G       C       100     PASS    AN=2;AC=1;Prej=AddedVariant     GT:AD   ./.:.   ./.:.   1/0:20,15
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##dual_coordinates=PRIMARY
    ##chain=file:///data/GRCh37_to_GRCh38.chain.genozip
    ##luft_reference=file:///data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##reference=file:///references/grch37/reference.bin
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles",RendAlg=R>
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed",RendAlg=A_1>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",RendAlg=GT>
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes",RendAlg=G>
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele",RendAlg=A_AN>
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes",RendAlg=NONE>
    ##INFO=<ID=LUFT,Number=4,Type=String,Description="Info for rendering variant in LUFT coords",RendAlg=NONE>
    ##INFO=<ID=PRIM,Number=4,Type=String,Description="Info for rendering variant in PRIMARY coords",RendAlg=NONE>
    ##INFO=<ID=Lrej,Number=1,Type=String,Description="Reason variant was rejected for LUFT coords",RendAlg=NONE>
    ##INFO=<ID=Prej,Number=1,Type=String,Description="Reason variant was rejected for PRIMARY coords",RendAlg=NONE>
    ##contig=<ID=1,length=249250621>
    ##luft_contig=<ID=chr1,length=248956422>
    ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification",RendAlg=G>
    ##bcftools_mergeVersion=1.11+htslib-1.11
    ##bcftools_mergeCommand=merge mydata.38.vcf.gz person3.vcf.gz; Date=Wed May 26 19:32:48 2021
    ##genozip_command="genocat newdata.38.vcf.genozip" 2021-05-26 19:37:29 ACDT
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Person1 Person2 Person3
    1       10285   1       T       C       4.4     PASS    AN=4;AC=3;LUFT=chr1,10285,T,-   GT:AD:AF:PL     0/1:31,18:0.367:37,0,46 1/1:.:.:.       ./.:.:.:.
    1       329162  3       A       T       4.6     PASS    AN=6;AC=5;LUFT=chr1,248466248,T,-       GT:AD:AF:PL:GP  0/1:28,9:0.3:36,0,0:.   1/1:.:.:.:.     1/1:3,30:.:.:0.14,14.9,37.9
    1       366043  4       CA      A       100     PASS    Lrej=RefTooLong GT      1|0     0|0
    
Notice that our new variant, ID=2, which is a Luft-only variant, is indeed absent from the variant data lines when rendering in Primary coordinates, and is instead present in the ``##luft_only`` meta-information line at the beginning of the file.

Also, look at the variant ID=3 - this is the variant with the REF⇆ALT switch. Notice how the GT, AD and GP fields were cross-rendered for the new sample Person3 to reflect the REF⇆ALT switch even though this sample was not present when we *lifted* the file with ``genozip --chain``. For the GP field, this happened even though the GP field was entirely absent from the file at the time of lifting:

  - GT: The values were flipped - 0/0 became 1/1
  - AD: The order of the values was reversed: 30,3 became 3,30
  - GP: The order of the values was reversed: 37.9,14.9,0.14 became 0.14,14.9,37.9 

This demonstrates how a VCF file can continue to evolve, adding and dropping variants, adding, removing or modifying INFO and FORMAT fields, all while continuing to maintain both coordinates, renderable in either Primary or Luft coordinates as required by any particular step in an analysis process.

Finally, let's look at ``--show-dvcf`` - a useful tool for getting visibility into how Genozip handled each variant. It may also be used in combination with ``--luft``:

::

    > genocat newdata.38.vcf.genozip --show-dvcf --header-one

    #COORD  oSTATUS  CHROM   POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Person1 Person2 Person3
    BOTH    OKRefSame       1       10285   1       T       C       4.4     PASS    AN=4;AC=3;LUFT=chr1,10285,T,-   GT:AD:AF:PL     0/1:31,18:0.367:37,0,46 1/1:.:.:.       ./.:.:.:.
    BOTH    OkRefAltSwitch  1       329162  3       A       T       4.6     PASS    AN=6;AC=5;LUFT=chr1,248466248,T,-       GT:AD:AF:PL:GP  0/1:28,9:0.3:36,0,0:.   1/1:.:.:.:.       1/1:3,30:.:.:0.14,14.9,37.9
    PRIM    RefTooLong      1       366043  4       CA      A       100     PASS    Lrej=RefTooLong GT      1|0     0|0

    > genocat newdata.38.vcf.genozip --show-dvcf --header-one --luft

    #COORD  oSTATUS  CHROM   POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Person1 Person2 Person3
    BOTH    OKRefSame       chr1    10285   1       T       C       4.4     PASS    AN=4;AC=3;PRIM=1,10285,T,-      GT:AD:AF:PL     0/1:31,18:0.367:37,0,46 1/1:.:.:.       ./.:.:.:.
    LUFT    AddedVariant    chr1    20197381        2       G       C       100     PASS    AN=2;AC=1;Prej=AddedVariant     GT:AD   ./.:.   ./.:.   1/0:20,15
    BOTH    OkRefAltSwitch  chr1    248466248       3       T       A       4.6     PASS    AN=6;AC=1;PRIM=1,329162,A,-     GT:AD:AF:PL:GP  1/0:9,28:0.7:0,0,36:.   0/0:.:.:.:.       0/0:30,3:.:.:37.9,14.9,0.14

Two columns were added at the start of each variant line - the first states the coordinates of the variant - Primary-only, Luft-only or Both coordinates. The second is the oSTATUS, states how the variant was handled. Dual-coordinate variants ("BOTH") have an oStatus of either ``OKRefSame`` indicating that the REF is the same (or reverse-complemented in XSTRAND=X) in both coordinates, or ``OkRefAltSwitch`` indicating the REF and ALT are switched between the two coordinates. At present, a variant with a REF change other than REF⇆ALT switch, cannot be lifted by Genozip and will be a Primary-only or Luft-only variant.

More inforation about how Genozip *lifts* and *cross-renders* each field can found here: :ref:`Rendering a DVCF <dvcf-rendering>`
