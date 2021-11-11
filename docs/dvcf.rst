.. _dvcf:

Dual-coordinate VCF files
==========================

.. include:: dvcf-see-also.rst

.. toctree::
    :hidden:
 
    Rendering a DVCF <dvcf-rendering>
    Renaming and dropping annotations in a DVCF <dvcf-renaming>
    Chain files <chain>
    Genozip DVCF lifting limitations <dvcf-limitations>
    
**At a glance**

Dual coordinates VCF files (or DVCFs), are VCF files that contain coordinates in two coordinate systems concurrently, for example, GRCh37 and GRCh38. A DVCF can be *rendered* in either *Primary* coordinates or *Luft* coordinates (*Luft* being our made-up past particile of *Lift*). 

Both *renditions* are standard-compliant VCF files, and both contain precisely the same information, with the only difference being how the information is presented.

Crucially, since both renditions contain all the information needed for rendering the file in both coordinate systems, the file in either rendition can be processed through any bioinformatics tool or pipeline, and still maintain its dual coordinates. In other words, an analytics pipeline may now have some steps that operate on the file in one coordinate system, and other steps operate in the other coordinate system, and the same data can flow through these steps seamlessly.

|

**An example - (1) using --chain to create a DVCF**

Let's walk through an example:

Consider this small VCF file, called ``mydata.vcf.gz``, which is in the coordinates defined by GRCh37, a popular Homo Sapiens reference:

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
    ##INFO=<ID=OBSCURE,Number=1,Type=Float,Description="Obscure annotation",RendAlg="A_1"
    ##contig=<ID=1,length=249250621>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1	Person2
    1	10285	.	T	C	4.4	PASS	AC=3;AN=4	GT:AD:AF:PL	0/1:31,18:0.367:37,0,46	1/1
    1	329162	.	A	T	4.6	PASS	AC=3;AN=4	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1
    1	366042	.	TA	A	100	PASS	.	GT	1|0	0|0
    1	20159588	.	C	A	100	PASS	.	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1
                    
We shall now convert this file to a Dual-coordinate VCF, containing both GRCh37 coordinates, as well as coordinates from another popular reference, GRCh38. This process is called *lifting*. We shall call GRCh37 the *Primary* coordinates, and GRCh38 the *Luft* coordinates.

Steps:

1. Prepare two reference files in genozip format: for Primary and Luft coordinates. The following results in the files named ``GRCh37.ref.genozip`` and ``GRCh38.ref.genozip`` :

::

    genozip --make-reference GRCh37.fa 
    genozip --make-reference GRCh38.fa 


2. Prepare a chain file in genozip format. Genozip uses chain files in the `UCSC chain file format <https://genome.ucsc.edu/goldenPath/help/chain.html>`_. Here are some :ref:`links to available chain files<chain>`. The following command line generates the file ``GRCh37_to_GRCh38.chain.genozip`` :

::

    genozip GRCh37_to_GRCh38.chain.gz --reference GRCh37.ref.genozip --reference GRCh38.ref.genozip

Note that the first two steps are preparation steps that need to be executed only once. 

1. Now that we have the chain file, we can convert any number of VCF files to DVCF:

:: 

    genozip mydata.vcf.gz --chain GRCh37_to_GRCh38.chain.genozip --add-line-numbers

Note that as usual with ``genozip``, .vcf .vcf.gz .vcf.bz2 and .vcf.xz are all accepted. Genozip does not require the VCF file to be indexed.

|

**An example - (2) rendering the DVCF in Primary and Luft coordinates**

Now, we have a dual-coordinate VCF file: ``mydata.vcf.genozip``. We can render it in either coordinate:

Rendering in Primary coordinates (GRCh37), followed by Luft coordinates (GRCh38):

::

    > genocat mydata.d.vcf.genozip | tee mydata.37.vcf

    ##fileformat=VCFv4.2
    ##original_reference=file:///references/grch37/reference.bin
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles",RendAlg="R"
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed",RendAlg="A_1"
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",RendAlg="GT"
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes",RendAlg="G"
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele",RendAlg="A_AN"
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes",RendAlg="NONE"
    ##INFO=<ID=OBSCURE,Number=1,Type=Float,Description="Obscure annotation",RendAlg="A_1"
    ##contig=<ID=1,length=249250621>
    ##dual_coordinates=PRIMARY
    ##chain=file:///C/Users/USER/projects/genozip/data/GRCh37_to_GRCh38.chain.genozip
    ##luft_reference=file://data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##reference=file://data/hs37d5.ref.genozip
    ##INFO=<ID=LUFT,Number=4,Type=String,Description="Info for rendering variant in LUFT coords. See https://genozip.com/dvcf.html",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=PRIM,Number=4,Type=String,Description="Info for rendering variant in PRIMARY coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=Lrej,Number=1,Type=String,Description="Reason variant was rejected for LUFT coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=Prej,Number=1,Type=String,Description="Reason variant was rejected for PRIMARY coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##luft_contig=<ID=chr1,length=248956422>
    ##genozip_command="C:\Users\USER\projects\genozip\genocat-debug.exe -o docs/dvcf-example-files/mydata.37.vcf -f docs/dvcf-example-files/mydata.d.vcf.genozip" 2021-06-15 11:36:11 Cen. Australia Daylight Time
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1	Person2
    1	10285	LN=1	T	C	4.4	PASS	AC=3;AN=4;LUFT=chr1,10285,T,-	GT:AD:AF:PL	0/1:31,18:0.367:37,0,46	1/1
    1	329162	LN=2	A	T	4.6	PASS	AC=3;AN=4;LUFT=chr1,490175,T,X	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1
    1	366042	LN=3	TA	A	100	PASS	Lrej=RefNewAlleleIndel	GT	1|0	0|0
    1	20159588	LN=4	C	A	100	PASS	AC=3;AN=4;AF=0.75;OBSCURE=0.002;LUFT=chr1,19833095,A,-	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1
                
    > genocat mydata.d.vcf.genozip --luft | tee mydata.38.vcf

    ##fileformat=VCFv4.2
    ##primary_only=1	366042	LN=1	TA	A	100	PASS	Lrej=RefNewAlleleIndel	GT	1|0	0|0
    ##original_reference=file:///references/grch37/reference.bin
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles",RendAlg="R"
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed",RendAlg="A_1"
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",RendAlg="GT"
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes",RendAlg="G"
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele",RendAlg="A_AN"
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes",RendAlg="NONE"
    ##INFO=<ID=OBSCURE,Number=1,Type=Float,Description="Obscure annotation",RendAlg="A_1"
    ##primary_contig=<ID=1,length=249250621>
    ##dual_coordinates=LUFT
    ##chain=file:///C/Users/USER/projects/genozip/data/GRCh37_to_GRCh38.chain.genozip
    ##reference=file://data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##primary_reference=file://data/hs37d5.ref.genozip
    ##INFO=<ID=LUFT,Number=4,Type=String,Description="Info for rendering variant in LUFT coords. See https://genozip.com/dvcf.html",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=PRIM,Number=4,Type=String,Description="Info for rendering variant in PRIMARY coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=Lrej,Number=1,Type=String,Description="Reason variant was rejected for LUFT coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=Prej,Number=1,Type=String,Description="Reason variant was rejected for PRIMARY coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##contig=<ID=chr1,length=248956422>
    ##genozip_command="C:\Users\USER\projects\genozip\genocat-debug.exe -vo docs/dvcf-example-files/mydata.38.vcf -f docs/dvcf-example-files/mydata.d.vcf.genozip" 2021-06-15 11:36:09 Cen. Australia Daylight Time
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1	Person2
    chr1	10285	LN=1	T	C	4.4	PASS	AC=3;AN=4;PRIM=1,10285,T,-	GT:AD:AF:PL	0/1:31,18:0.367:37,0,46	1/1
    chr1	490175	LN=2	T	A	4.6	PASS	AC=3;AN=4;PRIM=1,329162,A,X	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1
    chr1	19833095	LN=4	A	C	100	PASS	AC=1;AN=4;AF=0.25;OBSCURE=0.998;PRIM=1,20159588,C,-	GT:AD:AF:PL	1/0:9,28:0.7:0,0,36	0/0
                    
So what do we see here?

The first thing to appreciate, is that both ``mydata.37.vcf`` and ``mydata.38.vcf`` contain precisely the same information. In other words, you can ``genozip`` either file, and then render it again in both Primary and Luft coordinates.

Notice the ID field contains a sequential line number (starting with ``LN=``). This is a result of the ``genozip --add-line-numbers`` option and it will make it easier for us to correspond variants in the Primary rendition to the Luft rendition.

Take a look at the first file - the Primary rendition - at the variant with LN=4. The INFO/LUFT field contains the information of the Luft coordinate, which, in this case, is mapped to the Luft reference at CHROM=chr1 POS=19833095, and the REF (the third value in INFO/LUFT) is A - which would be change from the Primary (GRCh37) REF that is C, and identical to the Primary ALT which is an A. This is common situation which we call a REF⇆ALT switch.

Now take a look at the second file - the Luft rendition - at LN=4. Indeed, the REF and ALT have switched - and also the INFO/AC, INFO/AF,INFO/OBSCURE, FORMAT/GT, FORMAT/AD, FORMAT/AF and FORMAT/PL fields have changed too, to reflect this REF⇆ALT switch. In addition, the INFO/PRIM field contains the Primary coordinates of this variant. 

The algorithms determining how to *cross-rendered* each INFO and FORMAT field, are defined in the ``RendAlg`` attribute of the ``##INFO`` and ``##FORMAT`` meta-information lines, and you may modify them as you wish. Notice, for example, that INFO/OBSCURE, an annotation made up for the purpose of this example, was given a RendAlg=A_1 in the original ``mydata.vcf.gz`` file, an algorithm that cross-renders the value to (1-value) in case of a REF⇆ALT switch. To learn more about RendAlg options, see :ref:`dvcf-rendering`. Notice that we don't bother setting RendAlg for the other INFO and FORMAT fields in mydata.vcf.gz as they are all well known annotations for which Genozip's RendAlg choices behave as expected. These default choices can be seen in the ##INFO and ##FORMAT meta-information lines of both renditions. 

In contrast, take a look at variant LN=2. In the Primary rendition has REF=A ALT=T and it in the Luft rendition it has REF=T ALT=A. This appears to be a REF⇆ALT switch, however, the INFO and FORMAT fields were not cross-rendered. This is because this variant is actually *not* a REF⇆ALT switch - rather, the chain file mapped the alignment on which this variant resides to the reverse strand in the Luft reference, as indicated by the ``X`` in the fourth value of the INFO/PRIM field. The nominal base change is due to the required reverse-complementing of REF and ALT (as in the VCF file format variants are always shown in forward strand terms), and this variant has actually not changed between the references.

Let us now take a look the variant with LN=3. It has an INFO/Lrej field indicating that Genozip cannot render it in Luft coordinates and hence it is a Primary-only variant (or in other words, it was *rejected* from lifting). The reason given in this case is ``RefNewAlleleIndel`` - which means that neither the REF nor the ALT of the Primary indel variant match the Luft reference at this position. There are many reasons lifting can fail - they are listed in the table below.

In the Luft rendition, this Primary-only variant, LN=3, is indeed missing in the VCF data lines, which now contain only three variants. However, notice that it still survives as a meta-information line with the key ``##primary_only``.

|

**An example - (3) applying bioinformatics tools to a DVCF rendition**

Now, let us take our Luft rendition VCF ``mydata.38.vcf`` and merge it (using ``bcftools merge`` in this example) with another VCF file, ``person3.vcf``, which is in GRCh38 coordinates and has a single sample, Person3. Note that ``person3.vcf`` is *not* a DVCF. ``person3.vcf`` contains 2 variants on chr1: LN=1 (also present in ``mydata.38.vcf``), and LN=5 (not present in ``mydata.38.vcf``):

::

    > cat person3.vcf
    
    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person3
    chr1	10285	LN=1	T	C	100	PASS	AC=1;AN=2	GT:AD	1/0:20,15
    chr1	22000000	LN=5	A	C	4	PASS	AC=0;AN=2	GT:AD:GP	0/0:30,3:37.9,14.9,0.14
            
Notice that Person3 has the field FORMAT/GP in variant LN=5 that doesn't exist in ``mydata.38.vcf``. 

The merged file will be a DVCF file looking like this:

::

    > bgzip person3.vcf 
    > bgzip mydata.38.vcf
    > bcftools index mydata.38.vcf.gz 
    > bcftools index person3.vcf.gz
    > bcftools merge mydata.38.vcf.gz person3.vcf.gz | tee merged.38.vcf

    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##primary_only=1	366042	LN=1	TA	A	100	PASS	Lrej=RefNewAlleleIndel	GT	1|0	0|0
    ##original_reference=file:///references/grch37/reference.bin
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles",RendAlg="R"
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed",RendAlg="A_1"
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",RendAlg="GT"
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes",RendAlg="G"
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele",RendAlg="A_AN"
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes",RendAlg="NONE"
    ##INFO=<ID=OBSCURE,Number=1,Type=Float,Description="Obscure annotation",RendAlg="A_1"
    ##primary_contig=<ID=1,length=249250621>
    ##dual_coordinates=LUFT
    ##chain=file:///C/Users/USER/projects/genozip/data/GRCh37_to_GRCh38.chain.genozip
    ##reference=file://data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##primary_reference=file://data/hs37d5.ref.genozip
    ##INFO=<ID=LUFT,Number=4,Type=String,Description="Info for rendering variant in LUFT coords. See https://genozip.com/dvcf.html",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=PRIM,Number=4,Type=String,Description="Info for rendering variant in PRIMARY coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=Lrej,Number=1,Type=String,Description="Reason variant was rejected for LUFT coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=Prej,Number=1,Type=String,Description="Reason variant was rejected for PRIMARY coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##contig=<ID=chr1,length=248956422>
    ##genozip_command="C:\Users\USER\projects\genozip\genocat-debug.exe -vo docs/dvcf-example-files/mydata.38.vcf -f docs/dvcf-example-files/mydata.d.vcf.genozip" 2021-06-15 11:36:09 Cen. Australia Daylight Time
    ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification">
    ##bcftools_mergeVersion=1.11+htslib-1.11
    ##bcftools_mergeCommand=merge mydata.38.vcf.gz person3.vcf.gz; Date=Tue Jun 15 12:13:51 2021
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1	Person2	Person3
    chr1	10285	LN=1	T	C	100	PASS	PRIM=1,10285,T,-;AN=6;AC=4	GT:AD:AF:PL	0/1:31,18:0.367:37,0,46	1/1:.:.:.	1/0:20,15:.:.
    chr1	490175	LN=2	T	A	4.6	PASS	PRIM=1,329162,A,X;AN=4;AC=3	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1:.:.:.	./.:.:.:.
    chr1	19833095	LN=4	A	C	100	PASS	AF=0.25;OBSCURE=0.998;PRIM=1,20159588,C,-;AN=4;AC=1	GT:AD:AF:PL	1/0:9,28:0.7:0,0,36	0/0:.:.:.	./.:.:.:.
    chr1	22000000	LN=5	A	C	4	PASS	AN=2;AC=0	GT:AD:GP	./.:.:.	./.:.:.	0/0:30,3:37.9,14.9,0.14
    
The merged file, as expected, is in Luft (i.e. GRCh38) coordinates, and the new variant, LN=5, has no DVCF INFO tag (i.e. no INFO/PRIM).

|

**An example - (4) the modified file is still a DVCF**

Now let us genozip the merged file and render it first in Luft coordinates:

::

    > genozip merged.38.vcf
    genozip merged.38.vcf : 0%
    FYI: the number of samples in variant CHROM=1 POS=366043 is 2, different than the VCF column header line which has 3 samples
    Done (0 seconds)

    > genocat merged.38.vcf.genozip --luft
    
    ##fileformat=VCFv4.2
    ##primary_only=1	366042	LN=1	TA	A	100	PASS	Lrej=RefNewAlleleIndel	GT	1|0	0|0
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##original_reference=file:///references/grch37/reference.bin
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles",RendAlg="R"
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed",RendAlg="A_1"
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",RendAlg="GT"
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes",RendAlg="G"
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele",RendAlg="A_AN"
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes",RendAlg="NONE"
    ##INFO=<ID=OBSCURE,Number=1,Type=Float,Description="Obscure annotation",RendAlg="A_1"
    ##primary_contig=<ID=1,length=249250621>
    ##dual_coordinates=LUFT
    ##chain=file:///C/Users/USER/projects/genozip/data/GRCh37_to_GRCh38.chain.genozip
    ##reference=file://data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##primary_reference=file://data/hs37d5.ref.genozip
    ##INFO=<ID=LUFT,Number=4,Type=String,Description="Info for rendering variant in LUFT coords. See https://genozip.com/dvcf.html",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=PRIM,Number=4,Type=String,Description="Info for rendering variant in PRIMARY coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=Lrej,Number=1,Type=String,Description="Reason variant was rejected for LUFT coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=Prej,Number=1,Type=String,Description="Reason variant was rejected for PRIMARY coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##contig=<ID=chr1,length=248956422>
    ##genozip_command="C:\Users\USER\projects\genozip\genocat-debug.exe -vo docs/dvcf-example-files/mydata.38.vcf -f docs/dvcf-example-files/mydata.d.vcf.genozip" 2021-06-15 11:36:09 Cen. Australia Daylight Time
    ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification",RendAlg="G"
    ##bcftools_mergeVersion=1.11+htslib-1.11
    ##bcftools_mergeCommand=merge mydata.38.vcf.gz person3.vcf.gz; Date=Tue Jun 15 12:13:51 2021
    ##genozip_command="C:\Users\USER\projects\genozip\genocat-debug.exe -vo ../../junk.vcf docs/dvcf-example-files/merged.38.vcf.genozip" 2021-06-15 12:18:35 Cen. Australia Daylight Time
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1	Person2	Person3
    chr1	10285	LN=1	T	C	100	PASS	AN=6;AC=4;PRIM=1,10285,T,-	GT:AD:AF:PL	0/1:31,18:0.367:37,0,46	1/1:.:.:.	1/0:20,15:.:.
    chr1	490175	LN=2	T	A	4.6	PASS	AN=4;AC=3;PRIM=1,329162,A,X	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1:.:.:.	./.:.:.:.
    chr1	19833095	LN=4	A	C	100	PASS	AF=0.25;OBSCURE=0.998;AN=4;AC=1;PRIM=1,20159588,C,-	GT:AD:AF:PL	1/0:9,28:0.7:0,0,36	0/0:.:.:.	./.:.:.:.
    chr1	22000000	LN=5	A	C	4	PASS	AN=2;AC=0;Prej=AddedVariant	GT:AD:GP	./.:.:.	./.:.:.	1/1:3,30:37.9,14.9,0.14
    
Looking at the Luft rendition of the merged file, we notice that the new variant, LN=5, now contains a INFO/Prej field. Therefore, this is a Luft-only variant that cannot be expressed in Primary coordinates. This is always the case for variants that are added to an existing DVCF file, unless the new variants also originate from a DVCF file.

Also, note that when running ``genozip``, we got a warning about the number of samples of variant LN=1. This is our Primary-only variant, that didn't participate in the merge in Luft coordinates, and remains with two samples only. This is ok, as long as the DVCF file is the first file on the ``bctools merge`` command line, so that the missing sample aligns with the new sample Person3. It is perfectly acceptable to have missing samples in VCF.

Let us now take a look at the Primary rendition of the merged file:

::

    > genocat merged.38.vcf.genozip
    
    ##fileformat=VCFv4.2
    ##luft_only=chr1	22000000	LN=5	A	C	4	PASS	AN=2;AC=0;Prej=AddedVariant	GT:AD:GP	./.:.:.	./.:.:.	0/0:30,3:37.9,14.9,0.14
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##original_reference=file:///references/grch37/reference.bin
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles",RendAlg="R"
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed",RendAlg="A_1"
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",RendAlg="GT"
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes",RendAlg="G"
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele",RendAlg="A_AN"
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes",RendAlg="NONE"
    ##INFO=<ID=OBSCURE,Number=1,Type=Float,Description="Obscure annotation",RendAlg="A_1"
    ##contig=<ID=1,length=249250621>
    ##dual_coordinates=PRIMARY
    ##chain=file:///C/Users/USER/projects/genozip/data/GRCh37_to_GRCh38.chain.genozip
    ##luft_reference=file://data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##reference=file://data/hs37d5.ref.genozip
    ##INFO=<ID=LUFT,Number=4,Type=String,Description="Info for rendering variant in LUFT coords. See https://genozip.com/dvcf.html",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=PRIM,Number=4,Type=String,Description="Info for rendering variant in PRIMARY coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=Lrej,Number=1,Type=String,Description="Reason variant was rejected for LUFT coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##INFO=<ID=Prej,Number=1,Type=String,Description="Reason variant was rejected for PRIMARY coords",Source="genozip",Version="12.0.25",RendAlg="NONE"
    ##luft_contig=<ID=chr1,length=248956422>
    ##genozip_command="C:\Users\USER\projects\genozip\genocat-debug.exe -vo docs/dvcf-example-files/mydata.38.vcf -f docs/dvcf-example-files/mydata.d.vcf.genozip" 2021-06-15 11:36:09 Cen. Australia Daylight Time
    ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification",RendAlg="G"
    ##bcftools_mergeVersion=1.11+htslib-1.11
    ##bcftools_mergeCommand=merge mydata.38.vcf.gz person3.vcf.gz; Date=Tue Jun 15 12:13:51 2021
    ##genozip_command="C:\Users\USER\projects\genozip\genocat-debug.exe -o ../../junk.vcf docs/dvcf-example-files/merged.38.vcf.genozip" 2021-06-15 12:20:53 Cen. Australia Daylight Time
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1	Person2	Person3
    1	10285	LN=1	T	C	100	PASS	AN=6;AC=4;LUFT=chr1,10285,T,-	GT:AD:AF:PL	0/1:31,18:0.367:37,0,46	1/1:.:.:.	1/0:20,15:.:.
    1	329162	LN=2	A	T	4.6	PASS	AN=4;AC=3;LUFT=chr1,490175,T,X	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1:.:.:.	./.:.:.:.
    1	366042	LN=1	TA	A	100	PASS	Lrej=RefNewAlleleIndel	GT	1|0	0|0
    1	20159588	LN=4	C	A	100	PASS	AF=0.75;OBSCURE=0.002;AN=4;AC=3;LUFT=chr1,19833095,A,-	GT:AD:AF:PL	0/1:28,9:0.3:36,0,0	1/1:.:.:.	./.:.:.:.
        
Notice that our new variant, LN=5, which is a Luft-only variant, is indeed absent from the variant data lines when rendering in Primary coordinates, and is instead present in the ``##luft_only`` meta-information line at the beginning of the file.

Also, look at the variant LN=4 - this is the variant with the REF⇆ALT switch. Notice how the INFO and FORMAT fields were cross-rendered for the new sample Person3 to reflect the REF⇆ALT switch even though this sample was not present when we *lifted* the file with ``genozip --chain``. For the GP field, this happened even though the GP field was entirely absent from the file at the time of lifting:

  - GT: The values were flipped - 0/0 became 1/1
  - AD: The order of the values was reversed: 30,3 became 3,30
  - GP: The order of the values was reversed: 37.9,14.9,0.14 became 0.14,14.9,37.9 
  - INFO/AC: The new value is (AN-value)
  - INFO/AF, INFO/AC and INFO/OBSCURE: The value is (1-value) 

This demonstrates how a VCF file can continue to evolve, adding and dropping variants, adding, removing or modifying INFO and FORMAT fields, all while continuing to maintain both coordinates, renderable in either Primary or Luft coordinates as required by any particular step in an analysis process.

|

**Annotation dropping and renaming**

In some cases, an annotation tag (rather than value) changes between the Primary and Luft rendition, and in other cases it is desireable to drop annotations in the Luft rendition.

Genozip renames and drops certain annotations by default (which may be overridden if needed) and additional renaming or dropping may be specified using the command line options ``--dvcf-rename`` and ``--dvcf-drop``.

See: :ref:`Renaming and dropping annotations in a DVCF <dvcf-renaming>`

|

**--single-coord: Converting back to a single coordinate (i.e. normal) VCF file**

| ``genocat --single-coord myfile.d.vcf.genozip``
| ``genocat --single-coord --luft myfile.d.vcf.genozip``

The option ``--single-coord`` removes all DVCF-specific lines from the VCF header, and removes the DVCF INFO annotations, leaving the file as a normal VCF file in single coordinates - either the Primary coordinates, or Luft coordinates (when combined with ``--luft``).

**--show-dvcf and --show-ostatus - variant by variant info**

Now let's look at a interesting analytics tool: ``--show-dvcf`` is a useful tool for getting visibility into how Genozip handled each variant. It may also be used in combination with ``--luft``:

::

    > genocat merged.38.vcf.genozip --show-dvcf --header-one

    #COORD  oSTATUS CHROM   POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Person1 Person2 Person3
    BOTH    OkRefSameSNP    1       10285   LN=1    T       C       100     PASS    AN=6;AC=4;LUFT=chr1,10285,T,-   GT:AD:AF:PL     0/1:31,18:0.367:37,0,46 1/1:.:.:.       1/0:20,15:.:.
    BOTH    OkRefSameSNP    1       329162  LN=2    A       T       4.6     PASS    AN=4;AC=3;LUFT=chr1,490175,T,X  GT:AD:AF:PL     0/1:28,9:0.3:36,0,0     1/1:.:.:.       ./.:.:.:.
    PRIM    RefNewAlleleIndel       1       366042  LN=1    TA      A       100     PASS    Lrej=RefNewAlleleIndel  GT      1|0     0|0
    BOTH    OkRefAltSwitchSNP       1       20159588        LN=4    C       A       100     PASS    AF=0.75;OBSCURE=0.002;AN=4;AC=3;LUFT=chr1,19833095,A,-  GT:AD:AF:PL     0/1:28,9:0.3:36,0,0   
            1/1:.:.:.       ./.:.:.:.
            
    > genocat merged.38.vcf.genozip --show-dvcf --header-one --luft

    #COORD  oSTATUS CHROM   POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Person1 Person2 Person3
    BOTH    OkRefSameSNP    chr1    10285   LN=1    T       C       100     PASS    AN=6;AC=4;PRIM=1,10285,T,-      GT:AD:AF:PL     0/1:31,18:0.367:37,0,46 1/1:.:.:.       1/0:20,15:.:.
    BOTH    OkRefSameSNP    chr1    490175  LN=2    T       A       4.6     PASS    AN=4;AC=3;PRIM=1,329162,A,X     GT:AD:AF:PL     0/1:28,9:0.3:36,0,0     1/1:.:.:.       ./.:.:.:.
    BOTH    OkRefAltSwitchSNP       chr1    19833095        LN=4    A       C       100     PASS    AF=0.25;OBSCURE=0.998;AN=4;AC=1;PRIM=1,20159588,C,-     GT:AD:AF:PL     1/0:9,28:0.7:0,0,36   
            0/0:.:.:.       ./.:.:.:.
    LUFT    AddedVariant    chr1    22000000        LN=5    A       C       4       PASS    AN=2;AC=0;Prej=AddedVariant     GT:AD:GP        ./.:.:. ./.:.:. 1/1:3,30:37.9,14.9,0.14

Two columns were added at the start of each variant line - the first states the coordinates of the variant - Primary-only, Luft-only or Both coordinates. The second is the oSTATUS, states how the variant was handled. Dual-coordinate variants ("BOTH") have an oStatus starting with Ok, while Primary-only and Luft-only variants have a status indicating the reason for being rejected from the opposite rendition.

Alternatively, ``--show-ostatus`` may be used to add the oSTATUS to the INFO field:

::

    > genocat merged.38.vcf.genozip --show-ostatus --header-one

    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Person1 Person2 Person3
    1       10285   LN=1    T       C       100     PASS    AN=6;AC=4;LUFT=chr1,10285,T,-;oSTATUS=OkRefSameSNP      GT:AD:AF:PL     0/1:31,18:0.367:37,0,46 1/1:.:.:.       1/0:20,15:.:.       
    1       329162  LN=2    A       T       4.6     PASS    AN=4;AC=3;LUFT=chr1,490175,T,X;oSTATUS=OkRefSameSNP     GT:AD:AF:PL     0/1:28,9:0.3:36,0,0     1/1:.:.:.       ./.:.:.:.
    1       366042  LN=1    TA      A       100     PASS    Lrej=RefNewAlleleIndel;oSTATUS=RefNewAlleleIndel        GT      1|0     0|0
    1       20159588        LN=4    C       A       100     PASS    AF=0.75;OBSCURE=0.002;AN=4;AC=3;LUFT=chr1,19833095,A,-;oSTATUS=OkRefAltSwitchSNP        GT:AD:AF:PL     0/1:28,9:0.3:36,0,0 
            1/1:.:.:.       ./.:.:.:.
    
More inforation about how Genozip *lifts* and *cross-renders* each field can found here: :ref:`Rendering a DVCF <dvcf-rendering>`

|

**The rejects file and --show-lifts**

When running ``genozip --chain``, a human-readable *rejects file*, ``mydata.d.vcf.genozip.rejects`` in our case, is generated containing additional information regarding the rejection reasons related to each rejected variant:

::

    > cat mydata.d.vcf.genozip.rejects

    ##fileformat=GENOZIP-REJECTSv12.0.33
    #oSTATUS          CHROM   POS      REF     ALT     REASON
    RefNewAlleleIndel 1       366042   TA      A       Genozip limitation: new indel allele: PRIM=TAA LUFT(453295)=TAAAA

Here we can see why TA would be a new allele in the Luft reference: it is because the Primary REF allele relates to T followed by two As, and hence the ALT deletion allele refers to a T followed by a single A. Since the Luft reference contains four repeating As, and therefore the corresponding variant in Luft coordinates would therefore be "REF=TAAA ALT=TA,T", however at present Genozip is limited in that it cannot change a bi-allelic variant to a tri-allelic one, and hence the variant is rejected.

The option ``--show-lifts`` in combination with ``genozip --chain`` causes all lifts (not just rejected ones) to be reported in the rejects file:

::

    > genozip mydata.vcf.gz --chain GRCh37_to_GRCh38.chain.genozip --show-lifts
    > cat mydata.d.vcf.genozip.rejects 

    ##fileformat=GENOZIP-REJECTSv12.0.33
    #oSTATUS          CHROM   POS      REF     ALT     REASON
    OkRefSameSNP      1       10285    T       .       SNP: REF unchanged
    OkRefSameSNP      1       329162   A       .       SNP: REF unchanged
    RefNewAlleleIndel 1       366042   TA      A       Genozip limitation: new indel allele: PRIM=TAA LUFT(453295)=TAAAA
    OkRefAltSwitchSNP 1       20159588 C       A       SNP: REF<>ALT switch
    
|

**Overlapping variants and the .overlaps file** 

There are cases where two or more variants which have distinct Primary coordinates, are mapped to the same Luft coordinates, due to overlapping alignments in the chain file. When this occurs, Genozip generates a file with the ``.overlaps`` extension, containing the CHROM and POS of overlapping variants in Luft coordintes. This file can be used to view these variants:

``genocat myfile.d.vcf.genozip --luft --regions-file myfile.d.vcf.genozip.overlaps``

Filtering out these variants can be done with:

``genocat myfile.d.vcf.genozip --luft --regions-file ^myfile.d.vcf.genozip.overlaps``

The format of the ``.overlaps`` file, after removing the comment lines, is also compatible with ``bcftools view --regions-file``.

**Viewing summary statistics**

To see summary statistics of how variants were handled, we can use --show-counts (this works both with genozip and genocat):

::

    > genocat --show-counts=o\$TATUS merged.38.vcf.genozip
    
    Showing counts of o$TATUS (did_i=16). Total items=5 Number of categories=4
    OkRefSameSNP    2       40.00%
    RefNewAlleleIndel       1       20.00%
    AddedVariant    1       20.00%
    OkRefAltSwitchSNP       1       20.00%
    
    > genocat --show-counts=COORDS merged.38.vcf.genozip                                                                                                                             

    Showing counts of COORDS (did_i=15). Total items=5 Number of categories=3
    BOTH    3       60.00%
    LUFT    1       20.00%
    PRIM    1       20.00%

|

**Other useful options**

- ``genocat --show-chain <my-chain.chain.genozip>`` - displays all the chain file alignments, including overlap information.

    |

- ``genocat --snps-only <myfile.vcf.genozip>`` - drops non-SNP variants.

    | 

- ``genocat --indels-only <myfile.vcf.genozip>`` - drops non-indel variants.

    | 

- ``genocat --reference <ref-file> --regions chr22:17000000-17000100`` (or equivalently: ``chr22:17000000+101``) - displays a region of a reference.

    |

- ``genocat --reference <ref-file> --regions chr22:17000000-16999900`` (or equivalently: ``chr22:17000000-101``) - displays a region in reverse-complement.

    |

- ``genozip --chain <my-chain.chain.genozip> --add-line-numbers <myfile.vcf>`` - replaces the ID field with a sequential variant number, making it easier to locate corresponding variants between Primary and Luft renditions.

|

**oSTATUS list**

The oSTATUSes produced are as follows:

========================= ==================================
oSTATUS                   Occurs when...
========================= ==================================
OkRefSameSNP              Lift succeeded - REF unchanged for a SNP variant
OkRefSameIndel            Lift succeeded - REF unchanged for an indel variant
OkRefSameNDNIRev          Lift succeeded - Lift succeeded - REF unchanged for an left-anchored non-Ins non-Del indel variant with strand reversal
OkRefSameDelRev           Lift succeeded - Lift succeeded - REF unchanged for an left-anchored Deletion variant with strand reversal
OkRefSameInsRev           Lift succeeded - Lift succeeded - REF unchanged for an left-anchored Insertion variant with strand reversal
OkRefSameNotLeftAnc       Lift succeeded - REF unchanged for a variant that is neither a SNP nor left-anchored
OkRefSameStructVariant    Lift succeeded - REF unchanged for a variant with a symbolic ALT allele
OkRefAltSwitchSNP         Lift succeeded - REF⇆ALT switch for a SNP variant
OkRefAltSwitchIndel       Lift succeeded - REF⇆ALT switch for an left-anchored indel variant
OkRefAltSwitchDelToIns    Lift succeeded - Indel REF⇆ALT switch: The Deletion variant in Primary is incorporated in the Luft reference
OkRefAltSwitchIndelRpts   Lift succeeded - Indel REF⇆ALT switch: Switched number of payload repeats in reference
OkRefAltSwitchIndelFlank  Lift succeeded - Indel REF⇆ALT switch: REF bases the same, but switch called based on flanking regions
OkRefAltSwitchNDNI        Lift succeeded - Indel REF⇆ALT switch: Left-anchored non-Ins non-Del indel
OkRefAltSwitchWithGap     Lift succeeded - Indel REF⇆ALT switch: Deletion with payload in chain file gap
OkRefAltSwitchNotLeftAnc  Lift succeeded - REF⇆ALT switch for a variant that is neither a SNP nor left-anchored
OkNewRefSNP               Lift succeeded - REF changed (not to ALT) for a SNP variant with AF=1 or AC=AN
ChromNotInPrimReference   Lift failed - Primary coordinates reference file doesn't contain CHROM 
ChromNotInChainFile       Lift failed - chain file doesn't have a mapping for CHROM 
RefNotMappedInChain       Lift failed - Chain file doesn't contain a mapping covering REF
NewAnchorNotInChrom       Lift failed - New left-anchor base (after reverse-complementing) is before beginning of the chromosome
RefSplitInChain           Lift failed - REF is not fully within a single alignment in the chain file
RefMismatchesReference    Lift failed - REF in VCF file mismatches the Primary coordinates reference file. There is an error in the VCF file, or the wrong reference file is provided
RefMultiAltSwitchSNP      Lift failed - A REF⇆ALT switch for a multi-allelic SNP
RefMultiAltSwitchIndel    Lift failed - A REF⇆ALT switch for a multi-allelic indel
RefNewAlleleSNP           Lift failed - The Luft reference would introduce a new allele which is neither REF or ALT
RefNewAlleleDelRefChgHas* Lift failed - REF changed in a Deletion variant that has a "*" ALT
RefNewAlleleDelRefChanged Lift failed - REF changed in Deletion variant, but not REF<>ALT switch (i.e. Deletion not integrated into new reference)
RefNewAlleleDelSameRef    Lift failed - REF bases match, but this is a new Deletion allele based on context
RefNewAlleleInsSameRef    Lift failed - REF bases match, but this is a new Insertion allele based on context
RefNewAlleleIndelNoSwitch Lift failed - REF switched with one of the ALTs, but flanking regions differ
RefNewAlleleNDNI          Lift failed - REF is a new allele a left-anchored non-Ins non-Del indel variant
RefNewAlleleNotLeftAnc    Lift failed - The Luft reference would introduce a new allele for a variant that is neither a SNP nor left-anchored
RefNewAlleleSV            Lift failed - The Luft reference is different than REF for a variant with a symbolic ALT allele
XstrandNotLeftAnc         Lift failed - A variant that is neither a SNP nor left-anchored mapped to the reverse strand
XstrandSV                 Lift failed - A variant with a symbolic ALT allele mapped to the reverse strand
ComplexRearrangements     Lift failed - Variant is a Complex Rearrangement
INFO/END                  Lift failed - POS and INFO/END are not on the same chain file alignment
AddedVariant              Failed to cross-render a variant - Non-dual-coordinate variant added to a DVCF 
UnsupportedRefAlt         Failed to cross-render a variant - Combination of REF/ALT in Primary and Luft coordinates not supported by Genozip
INFO/AC                   Failed to cross-render AC: INFO/AN is missing for this variant 
INFO/*tag*                Failed to cross-render this INFO tag
FORMAT/*tag*              Failed to cross-render this FORMAT tag
========================= ==================================
