Dual-coordinates VCF files
==========================

Working with dual coordinates consists of three processes:

1. *Alignment* of a VCF file to a laft reference. This results in a *dual-coordinates VCF*, which is the original VCF with one additional INFO subfield per variant, which can be either LIFTOVER or LIFTREJD. The INFO/LIFTOVER subfield contains information a Liftover tool may use to liftover the VCF to the laft coordinate system, while the INFO/LIFTREJD subfield indicates that reason as to why this variant cannot be lifted over. The method of aligning a VCF to the laft reference is determined by the Aligner tool and is out of scope of this specification.

2. *Liftover* of a dual-coordinates VCF file. This results in a *Laft VCF* as follows:
   
    2.1 For variants that have a INFO/LIFTOVER subfield in the dual-coordinates VCF: These variants are lifted-over and some fields, INFO subfields and CHROM subfields are modified. The INFO/LIFTOVER subfield is replaced with an INFO/LIFTBACK subfield, described below, containing all the information a liftover software needs to liftback this variant to the primary coordinates. The selection of which fields or subfields to modify in the liftover process, and the algorithm to do so, are determined by the Liftover tool used and is out of scope of this specification.

    
    2.2 Variants with a INFO/LIFTREJT subfield in the dual-coordinates VCF: The entire variant line is moved, as is, to the VCF header, and prefixed with "##LIFTOVER_REJECT=". 
 
 
3. *Liftback* of a Laft VCF file. This results in a *dual-coordinates VCF*. 

    3.1 All variants are lifted-backed to the primary coordinate system using the information in INFO/LIFTBACK. The INFO/LIFTBACK subfield is replaced by an INFO/LIFTOVER subfield which is expected to be identical to the INFO/LIFTOVER subfield of the original dual-coordinates VCF prior to liftover. 
    
    3.2 Variants contained in ##LIFTOVER_REJECT header lines are extracted and added back as variants, and the header lines removed. Liftover tools are encouraged to offer the user an option to drop these variants altogether.


**Extensions to the primary VCF file** 

1. INFO/LIFTOVER subfield:

LIFTOVER=*CHROM*,*POS*,*REF*,*STRAND* 

Example: ``LIFTOVER=chr1,24143532,A,+``

============ ==================================================================================================
*OCHROM*     The value of the CHROM field in the *Laft VCF*
*OPOS*       The value of the POS field in the *Laft VCF*
*OREF*       The value of the REF field in the *Laft VCF*
*OSTRAND*    "-" if the ALT alleles should be reverse-complemented during lift-over (after any changes in REF and ALT), and "+" if not
*OALTRULE*   One of: Y or N: determines how to calculate OALT
             Y: if *OREF* equals to REF, then OALT=ALT
             if *OREF* is equal to one of the alleles in the ALT, then OALT=ALT with that allele replaced by REF
             if *OREF* is not equal to any REF or ALT allele, then OALT=ALT concatenated with REF
             N: OALT=ALT
============ ==================================================================================================

The corresponding header line is: 
``##INFO=<ID=LIFTOVER,Number=4,Type=String,Description="dual-coordinates VCF: Information for lifting over the variant to laft coordinates">``

1. VCF header ##DUAL_COORDINATES key 

This key must exist for a dual coordinates file, and must contain one of the values: PRIMARY or LAFT, indicating the coordinates of the file.

2. VCF header ##LIFTOVER_REJECT header key
   
Exists in LAFT files only, one per variant rejected by lift over. The value is the variant line as-is.

3. INFO/LIFTREJD subfield: 
   
LIFTREJD=*REASON*

Example: ``LIFTREJD=NO_CHROM``. 

It represents the reasons for liftover rejection. The following reasons are defined, and a Liftover tool may add additional ones:

=========== ==================================================================================================
*String*    *Rejection reason*
OK          Liftover successful (reserved, but usually doesn't appear in INFO/LIFTREJD)
NO_CHROM    CHROM is not available in lift-over reference
NO_MAPPING  POS does not map to the lift-over reference
REF_SPLIT   The first and last base in a multi-base REF map to different alignments in the lift-over reference
UNSUPPORTED The variant is mappable in principal, but this case is not supported by the liftover tool
=========== ==================================================================================================

The corresponding header line is: 
``##INFO=<ID=LIFTREJD,Number=1,Type=String,Description="dual-coordinates VCF: Reason variant was rejected for lift over">``

4. INFO/LIFTBACK subfield:

LIFTBACK=*CHROM*,*POS*,*STRAND*,*REF*,*ALTRULE* 

Example: ``LIFTBACK=chr1,24143532,A,+,Y``

=========== ==================================================================================================
*CHROM*     The value of the CHROM field in the *Primary VCF*
*POS*       The value of the POS field in the *Primary VCF*
*REF*       The value of the REF field in the *Primary VCF*
*STRAND*    "-" if Laft's ALT alleles should be reverse-complemented before any changes to REF and ALT, "+" if not
*ALTRULE*   One of: Y or N: determines how to calculate ALT
            Y: if *REF* equals to OREF, then ALT=OALT
            if *REF* is equal to one of the alleles in the OALT, then ALT=OALT with that allele replaced by OREF
            if *REF* is not equal to any OREF or OALT allele, then ALT=OALT concatenated with OREF
            N: ALT=OALT
=========== ==================================================================================================

The corresponding header line is: 
``##INFO=<ID=LIFTBACK,Number=5,Type=String,Description="dual-coordinates VCF: Information for retrieving the variant in the primary coordinates">``

