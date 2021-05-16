.. _dc-vcf-spec:

The Variant Call Format (VCF) - Dual Coordinates Extension
==========================================================

Version 1.0, March 1, 2021

Written by: Divon Lan, divon@shablife.com

**Contents**

1. Dual Coordinates VCF concept

2. An example
   
3. Meta information lines

4. Liftable variants

5. Non-liftable variants

6. Liftover of CHROM, POS, REF, ALT
   
7. Subfield liftover algorithms
   
8. Sorting
   
9.  File names


**1. Dual Coordinates VCF concept**

This specificiation is fully compatible with the `VCFv4.3 specification <https://samtools.github.io/hts-specs/VCFv4.3.pdf>`_, and extends it. It is also fully compatible with VCF v4.1 and v4.2.

The specification defines two formats of VCF, both fully compliant with the VCF specification, which are called the Primary Coordinates VCF file and the Luft Coordintes VCF file (or *PVCF* and *LVCF* in short). A pair of PVCF and LVCF files contain information about the same genetic variants in two different coordinate systems. The key feature of these files is that they contain precisely the same information, merely *rendered* in two different ways.

Once a VCF file is *lifted* to a Dual Coordinate VCF file (i.e. PVCF and/or LVCF) - it can be processed through an analytical pipeline, and since the data can be rendered in either coordinate systems, each stage of the pipeline can arbitrarily operate on either coordinate system. Importantly, the rendering continues to work as fields and annotations are added, removed or modified as the data works its way down the pipeline. Rendering is a fast operation that does not require a reference or chain file.

Definitions:

- A *Source VCF* is a VCF file that is not a PVCF or LVCF. 

- The *Primary coordinate system* is the coordinate system of the Source VCF.

- The *Luft coordinate system* is the other coordinate system in which the variant data will be expressed ("Luft" being a made-up past participle of "Lift").

- A *Primary VCF* (or *PVCF*) and a *Luft VCF* (or *LVCF*) are VCF files expressed in the Primary or Luft coordinates repectively, which are equivalent to each other and contain all the information of a Source VCF along with all the information needed to render the PVCF as a LVCF and vice versa.

- A *Lifter* is a software functionality that converts, or *Lifts*, a Source VCF to a PVCF and/or LVCF. It typically uses additional information such as a reference file in the Luft coordinates and a chain file.

- A *Renderer* is a software functionality that renders a PVCF as an LVCF (referred to as *lifting over*) or renders a LVCF as a PVCF ("*lifting back*"). It does not require any external information beyond the input PVCF or LVCF file itself. Specifically, it does not require a reference file or a chain file.

- An *implementation* means a particular software package including the functionalities of a Lifter and Renderer.
   
This specification defines the LVCF and PVCF formats. It does not define the algorithms of the Lifter or Renderer. Indeed, different implementations of Lifters and Translators might operate differently to address different needs. Adhering to this specification will ensure that the resulting files are inter-operable across different systems.

**2. An example** 

A PVCF:

::

    ##fileformat=VCFv4.3
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ###contig=<ID=1,length=249250621>
    ##contig=<ID=2,length=243199373>
    ##dual_coordinates=PRIMARY
    ##chain=file:///data/GRCh37_to_GRCh38.chain.genozip
    ##liftover_reference=file:///data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##INFO=<ID=LIFTOVER,Number=5,Type=String,Description="dual-coordinates VCF: Information for lifting over the variant to luft coordinates">
    ##INFO=<ID=LIFTBACK,Number=5,Type=String,Description="dual-coordinates VCF: Information for retrieving the variant in the primary coordinates">
    ##INFO=<ID=LIFTREJT,Number=1,Type=String,Description="dual-coordinates VCF: Reason variant was rejected for lift over">
    ##liftover_contig=<ID=chr1,length=248956422>
    ##liftover_contig=<ID=chr2,length=242193529>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1
    1	00329162	.	G	A	100	.	AN=2;AC=1;AF=0.500;LIFTREJT=REF_CHANGE_NOT_ALT	GT:AD:PL	1|0:11,41:81,17,0
    1	00348307	.	c	G	100	.	AN=2;AC=1;AF=0.50;LIFTOVER=chr1,248485393,X,1,g	GT:AD:PL	1|0:20,10:60,20,1
    1	20159588	.	C	a	100	.	AN=4;AC=1;AF=0.25;LIFTOVER=chr1,19833095,-,1,a	GT:AD:PL	1|0:11,41:109,60,0

The LVCF file corresponding to the PVCF above:

::

    ##fileformat=VCFv4.3
    ##liftback_reject=1	329162	.	G	A	100	.	AN=2;AC=1;AF=0.500;LIFTREJT=REF_CHANGE_NOT_ALT	GT:AD:PL	1|0:11,41:81,17,0
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ###contig=<ID=1,length=249250621>
    ##liftback_contig=<ID=2,length=243199373>
    ##dual_coordinates=LUFT
    ##chain=file:///data/GRCh37_to_GRCh38.chain.genozip
    ##reference=file:///data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    ##INFO=<ID=LIFTOVER,Number=5,Type=String,Description="dual-coordinates VCF: Information for lifting over the variant to luft coordinates">
    ##INFO=<ID=LIFTBACK,Number=5,Type=String,Description="dual-coordinates VCF: Information for retrieving the variant in the primary coordinates">
    ##INFO=<ID=LIFTREJT,Number=1,Type=String,Description="dual-coordinates VCF: Reason variant was rejected for lift over">
    ##contig=<ID=chr1,length=248956422>
    ##contig=<ID=chr2,length=242193529>
    ##genozip_command="C:\Users\USER\projects\genozip\genocat-debug.exe -v -fo example.l.vcf example.p.vcf.genozip" 2021-05-07 22:07:44 Cen. Australia Daylight Time
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1
    chr1	19833095	.	a	C	100	.	AN=4;AC=3;AF=0.75;LIFTBACK=1,20159588,-,1,C	GT:AD:PL	0|1:41,11:0,60,109
    chr1	248485393	.	g	C	100	.	AN=2;AC=1;AF=0.50;LIFTBACK=1,348307,X,1,c	GT:AD:PL	1|0:20,10:60,20,1

**3. Meta-information lines**

The following meta information lines are added or modified. Their order is in the file has no significance.

3.1 Coordinates

``##dual_coordinates=PRIMARY``

This field is required.

Permitted values: ``PRIMARY``, ``LUFT``. Defines the coordinates of the current rendering.

3.2 Chain file URL

This field is recommended.

``##chain=file:///data/GRCh37_to_GRCh38.chain.genozip``

The URL of the chain file used for generating this Dual Coordintes VCF. The file format is implementation-specific and out of scope of this specification.

3.3 Reference files' URLs

These field are recommended.

``##reference=file:///data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip``
``##primary_reference=file:///data/hg19.p13.plusMT.full_analysis_set.ref.genozip``
``##luft_reference=file:///data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip``

In a PVCF, ##reference contains the URL of the reference file of coordintes currently rendered (the Primary coordinates), and ##luft_reference contains the URL of the reference file of the Luft coordinates.

Similarly, in LVCF, we have ##reference and ##primary_reference.

3.4 Contigs

The ##contig key is as defined in the VCF specification. It refers to contigs of the current coordinates.

In a PVCF file, meta information lines with a ##luft_contig key may exist, and have an identical format ##contig. They describe the contigs that appear in the LVCF. It is recommended that a PVCF includes a ##luft-contig line for each contigs that appears in the LVCF file. 

Similarly, a LVCF file may contain ##primary_contig keys, describing the PVCF's contigs.

3.5 RenderAlg attribute of ##INFO and ##FORMAT

|``##FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype Likelihoods",RenderAlg=G>``
|``##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed",RenderAlg=A_AN>``

The RenderAlg attribute *must* be present in all INFO and FORMAT fields listed in the header. While it is not required that all INFO and FORMAT fields present in the file are listed in the header, only those listed in the header may be re-rendered by the Renderer. The Renderer *must* add RenderAlg to INFO or FORMAT lines in the header that are missing them (this is possible, for example, if the file acquired additional INFO or FORMAT fields in an analysis step). A user may edit the RenderAlg attribute and the Renderer *must* use the algorithm specified by the user, or error.

3.6 LIFTOVER, LIFTBACK and LIFTREJT


**7. Setting of oCHROM, oPOS, oREF and XSTRAND**

An implementation, when generating an PVCF and/or LVCF file **must** set XSTRAND to one of two values:
    - '-' if this variant has the same orientation in the reference data of both coordinate systems.
    - 'X' if this variant has the opposite orientations in the reference data of the two coordinate systems.
  
An implementation, when generating an PVCF and/or LVCF file from a source VCF, for each variant, must set the values of oCHROM, oPOS and oREF to be the genomically correct values in the Luft coordinates, generate (for LVCF) or be able to generate (for PVCF) the oALT based oCHROM, oPOS, oREF and XSTRAND. If it cannot do so, it **must** reject the variant.

**6. Subfield liftover algorithms**

If the REF changes, the implementation **must** either correct the REF and ALT fields, or reject the variant. The implementation **may** fix some of the INFO and FORMAT subfields, so that they are correctly reflect the changed REF and ALT fields.

An INFO for FORMAT subfield that an implementation fixes in case of a REF change, must contain an additional ``Liftover`` key its meta information line, indicating the algorithm used to lift it over.

::

    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",Liftover=GT>

Any algorithm **must** conform with the following three rules:

   1. To calculate the lifting, the algorithm may only use information contained within the same VCF line (variant). There is no limitation to the modifications that the algorithm may conduct on the VCF line and they may include other VCF fields or INFO or FORMAT subfields.

   2. The algorithm must be invertible, so that the inverse algorithm precisely undoes all the modifications to the VCF line.

   3. The inverse algorithm must adhere to these same three rules. 


The following algorithms are defined in this specification, and additional algorithms may be added by implementations. 

====  =======================  =========================================================================================================
ID    Relevant for (at least)  Algorithm
====  =======================  =========================================================================================================
GT    FORMAT/GT                Only if REF/ALT change: Change the allele numbers in-place (i.e. without changing their order) to reflect the REF and ALT changes
R     Subfields with Number=R  Only if REF/ALT change: Change the order of the elements in the array to reflect the REF / ALT changes
G     Subfields with Number=G  Only if REF/ALT change: Change the order of the elements in the array to reflect the REF / ALT changes
A_1   FORMAT/AF, INFO/AF       Only if REF/ALT change: For a Number=A field in which the sum of all values, including the implicit value relating to the REF allele, equals 1 - change to reflect the REF / ALT changes
A_AN  INFO/AC                  Only if REF/ALT change: For a Number=A field in which the sum of all values, including the implicit value relating to the REF allele, equals INFO/AN - change to reflect the REF / ALT changes
END   INFO/END                 Keep the same difference vs POS
====  =======================  =========================================================================================================

With the exception of END that applies to all lines with an INFO/END, the other algorithms listed do nothing to variants other than those with a REF change which is not marely a change in orientation.

An implementation, when generating the PVCF and LVCF files from a source VCF, **must** include an ##INFO or ##FORMAT meta information line for each subfield that is handled, even if such meta information lines are missing in the source VCF. Then, for each occurance of the subfield, it **must** either handle it according to the specified algorithm *or* reject the variant with the reason code LO_CANT_xxx where *xxx* is the ID of the algorithm, for example LO_CANT_GT. 

An implementation, when lifting over a PVCF file or lifting back a LVCF file, for each subfield that has a ``Liftover`` key in its meta information line: the implementation **must** lift-over or lift-back all non-rejected lines according to the specified algorithm. If an implementation is not capable of doing so (perhaps because the file was generated by another implementation with different capabilities), then it **must** error.

**7. Sorting**

The PVCF and LVCF files are both sorted, per the VCF specification, each in its respective coordinates. The ##liftback_reject meta-information lines need not be sorted.
