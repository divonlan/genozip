.. _dual_coordinates:

Dual-coordinates VCF files
==========================

At the core of the Dual Coordinate concept, is the fact that both these views are equivalent and contain the exact same information - they may be losslessly converted back and forth, without needing any external information (such as reference file or a chain file). Importantly, they may be converted back and forth even as they are modified - annotations added, changed or removed, variant order changed, variants filtered out etc - so long as the dual-coordinate information described below needed for dual-coordinate conversion is retained.

Working with dual coordinates consists of three processes:

1. *Alignment* of a VCF file to a luft reference. This results in a *dual-coordinates VCF*, which is the original VCF with one additional INFO subfield per variant, which can be either LIFTOVER or LIFTREJT. The INFO/LIFTOVER subfield contains information a Liftover tool may use to liftover the VCF to the luft coordinate system, while the INFO/LIFTREJT subfield indicates that reason as to why this variant cannot be lifted over. The method of aligning a VCF to the luft reference is determined by the Aligner tool and is out of scope of this specification.

2. *Liftover* of a dual-coordinates VCF file. This results in a *Luft VCF* as follows:
   
    2.1 For variants that have a INFO/LIFTOVER subfield in the dual-coordinates VCF: These variants are lifted-over and some fields, INFO subfields and CHROM subfields are modified. The INFO/LIFTOVER subfield is replaced with an INFO/LIFTBACK subfield, described below, containing all the information a liftover software needs to liftback this variant to the primary coordinates. The selection of which fields or subfields to modify in the liftover process, and the algorithm to do so, are determined by the Liftover tool used and is out of scope of this specification.

    
    2.2 Variants with a INFO/LIFTREJT subfield in the dual-coordinates VCF: The entire variant line is moved, as is, to the VCF header, and prefixed with "##LIFTOVER_REJECT=". 
 
 
3. *Liftback* of a Luft VCF file. This results in a *dual-coordinates VCF*. 

    3.1 All variants are lifted-backed to the primary coordinate system using the information in INFO/LIFTBACK. The INFO/LIFTBACK subfield is replaced by an INFO/LIFTOVER subfield which is expected to be identical to the INFO/LIFTOVER subfield of the original dual-coordinates VCF prior to liftover. 
    
    3.2 Variants contained in ##LIFTOVER_REJECT header lines are extracted and added back as variants, and the header lines removed. Liftover tools are encouraged to offer the user an option to drop these variants altogether.


**Extensions to the primary VCF file** 

1. INFO/LIFTOVER subfield:

LIFTOVER=*CHROM*,*POS*,*REF*,*STRAND* 

Example: ``LIFTOVER=chr1,24143532,A,+``

============ ==================================================================================================
*OCHROM*     The value of the CHROM field in the *Luft VCF*
*OPOS*       The value of the POS field in the *Luft VCF*
*XSTRAND*    "X" if the ALT alleles should be reverse-complemented during lift-over (after any changes in REF and ALT), and "-" if not
*OREF*       The value of the REF field in the *Luft VCF*
============ ==================================================================================================

The corresponding header line is: 
``##INFO=<ID=LIFTOVER,Number=4,Type=String,Description="dual-coordinates VCF: Information for lifting over the variant to luft coordinates">``

1. VCF header ##DUAL_COORDINATES key 

This key must exist for a dual coordinates file, and must contain one of the values: PRIMARY or LUFT, indicating the coordinates of the file.

2. VCF header ##LIFTOVER_REJECT header key
   
Exists in LUFT files only, one per variant rejected by lift over. The value is the variant line as-is.

3. INFO/LIFTREJT subfield: 
   
LIFTREJT=*REASON*

Example: ``LIFTREJT=NO_CHROM``. 

It represents the reasons for liftover rejection. The following reasons are defined, and a Liftover tool may add additional ones:

=========== ==================================================================================================
*String*    *Rejection reason*
OK          Liftover successful (reserved, but usually doesn't appear in INFO/LIFTREJT)
NO_CHROM    CHROM is not available in lift-over reference
NO_MAPPING  POS does not map to the lift-over reference
REF_SPLIT   The first and last base in a multi-base REF map to different alignments in the lift-over reference
UNSUPPORTED The variant is mappable in principal, but this case is not supported by the liftover tool
=========== ==================================================================================================

The corresponding header line is: 
``##INFO=<ID=LIFTREJT,Number=1,Type=String,Description="dual-coordinates VCF: Reason variant was rejected for lift over">``

4. INFO/LIFTBACK subfield:

LIFTBACK=*CHROM*,*POS*,*STRAND*,*REF* 

Example: ``LIFTBACK=chr1,24143532,=,Y,A``

=========== ==================================================================================================
*CHROM*     The value of the CHROM field in the *Primary VCF*
*POS*       The value of the POS field in the *Primary VCF*
*XSTRAND*   "X" if Luft's ALT alleles should be reverse-complemented before any changes to REF and ALT, "-" if not
*REF*       The value of the REF field in the *Primary VCF*
=========== ==================================================================================================

The corresponding header line is: 
``##INFO=<ID=LIFTBACK,Number=5,Type=String,Description="dual-coordinates VCF: Information for retrieving the variant in the primary coordinates">``

