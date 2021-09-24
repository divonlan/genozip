.. _match-chrom:

Matching contig names of the file to those in the reference file
================================================================

The unfortunate reality of the bioinformatics world is that contigs may appear with different names in different reference files, causing a problems in analysis pipelines.

Examples: 

| - chr22 ⇆ 22
| - M ⇆ chrM ⇆ MT ⇆ chrMT
| - chr21_gl000210_random ⇆ GL000210.1

Genozip offers a command line option, ``--match-chrom-to-reference``, that updates the contigs of a file to match those of the reference file.

``--match-chrom-to-reference`` works for all file types that have contigs: SAM/BAM, VCF, GFF3/GVF, 23andMe, and chain files.

**Example**

Notice that in the example below, the contig name **1** was updated to **chr1** both in the SAM header and the 3rd column of the data line, which is the RNAME field.

::
    
    > cat example.sam

    @HD	VN:1.4	SO:coordinate
    @SQ	SN:1	LN:249250621
    A00910:85:HYGWJDSXX:2:2502:31647:7701	99	1	9997	34	28M1I6M1I39M4D68M7S	=	10159	324	CCCTTAACCCTAACCCTAACCCTAACCCTTAACCCTTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAAACCCTAACCCTAACCCTAACCCTAACCCAAACCAAACCCTAACCCTAACCCTAACCCTAACCCTAACACCCAAA	FFFFFFFFFFF:FFFFFFF:FFFFFFFFF:F:FFFF:FFFFFFFFF:FFFFFFFFF:FF,:FFFFFFFFFFF,FFFFFFFF:FFF:::FFFF,F::FF:FFFFF::,FF,::FFF,:,FFF,,,,FF,::FFF:F,FF,,:FF:FFF,:,	AS:i:99	XS:i:96	MD:Z:0N0N0N0N69^CCCT29T4C33	NM:i:12	RG:Z:1

    > genozip example.sam --reference hg19.p13.plusMT.full_analysis_set.ref.genozip --match-chrom-to-reference
    genozip example.sam : Done (1 second, SAM compression ratio: 14.4)

    > genocat example.sam.genozip

    @HD	VN:1.4	SO:coordinate
    @SQ	SN:chr1	LN:249250621
    A00910:85:HYGWJDSXX:2:2502:31647:7701	99	chr1	9997	34	28M1I6M1I39M4D68M7S	=	10159	324	CCCTTAACCCTAACCCTAACCCTAACCCTTAACCCTTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAAACCCTAACCCTAACCCTAACCCTAACCCAAACCAAACCCTAACCCTAACCCTAACCCTAACCCTAACACCCAAA	FFFFFFFFFFF:FFFFFFF:FFFFFFFFF:F:FFFF:FFFFFFFFF:FFFFFFFFF:FF,:FFFFFFFFFFF,FFFFFFFF:FFF:::FFFF,F::FF:FFFFF::,FF,::FFF,:,FFF,,,,FF,::FFF:F,FF,,:FF:FFF,:,	AS:i:99	XS:i:96	MD:Z:0N0N0N0N69^CCCT29T4C33	NM:i:12	RG:Z:1

**Which fields are updated?**

========= ==========================================================================================
File type Fields updated
========= ==========================================================================================
SAM, BAM  @SQ lines in the file header; RNAME (column 3) and RNEXT (column 7); the contig name in the optional SA, OA and XA fields
VCF       ##contig lines in the file header; CHROM (column 1)
Chain     qName and tName fields
GFF3, GVF SequenceId (column 1)
23andMe   Chromosome (column 2)
========= ==========================================================================================


**How are contig names converted?**

Each contig name that appears in the file, is searched for in the reference file, in its unmodified form as well as variations of its name. Once a match is found in the reference file, the length of the contig, if known, is compared to make sure it is indeed the same contig. If a match is not found in the reference file, the contig name is not changed.

The variations of the contig names considered are:

- For chromosomes with a numeric suffix (eg chr22) as well as chrX, chrY, chrW and chrZ - with or without the *chr* prefix - the name with or without a lower case *chr* prefix is considered.

- For the Mitochondria chromosome, four options are considered - M, MT, chrM, chrMT

- For contigs that contain an embedded `Accession Number <https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/>`_, the following variations are considered:

========================== =======================================
Example contig name        Interpreted as...
========================== =======================================
chr4_gl383528_alt          Accession Number GL383528 version 1
chrUn_JTFH01001867v2_decoy Accession Number JTFH01001867 version 2
GL000192.1                 Accession Number GL000192 version 1
========================== =======================================

