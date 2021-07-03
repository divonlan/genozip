.. _dvcf-rendering:

Rendering a DVCF
======================

See also:

    | :ref:`Dual-coordinates VCF files <dvcf>`
    |
    | :ref:`Dual-coordinates VCF Specification <dvcf-spec>`
    |
    | :ref:`Chain files <dvcf-chain-files>`

This page describes how Genozip *lifts* a VCF file to be a :ref:`dual-coordinates VCF file<dvcf>` (DVCF), and how it *cross-renders* each field between Primary and Luft VCF renditions of a DVCF.

*Notation*: when ambiguous, we will use a subscript to indicate whether the field is in the Primary or Luft rendition, e.g. REF\ :subscript:`prim` and REF\ :subscript:`luft`

**Lifting vs cross-rendering**

When *lifting* a VCF file to a DVCF file using ``genozip --chain``, Genozip loads the chain file, the Primary-coordinates reference file, and the Luft-coordinates reference file. It then inspects each variant, and *lifts* its values of CHROM, POS and REF. By *lift* we mean that it calculates these values relative to the Luft reference file. In addition, it calculates a value called XSTRAND (pronounced: cross-strand), which indicates whether the chain file maps this coordinate to the forward or reverse strand of the Luft reference (see: `UCSC chain file format <https://genome.ucsc.edu/goldenPath/help/chain.html>`_).

When *rendering* a DVCF file in Luft coordinates with ``genocat --luft``, Genozip uses the lifted values of CHROM, POS and REF, and for all other fields, it *cross-renders* their values from Primary to Luft coordinates.

When *compressing* a Luft-coordinate DVCF rendition (initially created with ``genocat --luft`` and possibly modified after), Genozip *cross-renders* the Luft values in the Luft rendition file back to their Primary values, as values are always stored in the ``.vcf.genozip`` file in Primary-coordinate terms. The exception to this is if a Luft variant cannot be cross-rendered back to Primary. This can happen if it is a new variant added to the Luft file, or an existing Luft variant that had an INFO or FORMAT field added or modified in a way that cannot be cross-rendered. In this case, it becomes a Luft-only variant, and is stored in the .vcf.genozip file in Luft-coordinate terms. 

Important to note that only the initial lifting operation requires loading of the chain and reference files. Once a DVCF file is prepared, any rendering of it in Primary or Luft renditions, and any subsequent compression of these renditions with ``genozip`` (possibly after modifications with bioinformatics tools), does *not* require loading of the chain and reference files, as the DVCF file and its renditions already contain values in both coordinates.

Both *lifting* and *cross-rendering* may fail for any field, if the conditions for that particular operation are not met. A failure in any field causes the variant to become single-coordinate in the current coordinates - i.e. a Primary-only variant or a Luft-only variant. In this case, an INFO/Lrej will be created instead of an INFO/LUFT (for a Primary-only variant) or INFO/Prej will be created instead of INFO/PRIM (for a Luft-only variant). The INFO/Lrej or INFO/Prej field will contain the reason for rejection.

**CHROM and POS**

The CHROM\ :subscript:`luft` and POS\ :subscript:`luft` are initially obtained from the chain file when *lifting* a VCF file to a DVCF file using ``genozip --chain``, and placed in the INFO/LUFT field, as the first and second values respectively. Thereafter, when the file is rendered with ``--luft``, these Luft coordinates are presented in the VCF's CHROM and POS fields, while CHROM\ :subscript:`prim` and POS\ :subscript:`prim` are placed in the INFO/PRIM field.

**XSTRAND**

When *lifting* a variant during ``genozip --chain``, XSTRAND is a value placed as the fourth value in INFO/LUFT, and has the following value:
  - If qStrand in the chain file is '+', XSTRAND is set to ``-`` (hyphen)
  - If qStrand in the chain file is '-', XSTRAND is set to ``X`` (capital X)
  - Genozip requires tStrand to always be '+', which is the normal case for chain files.
  
When a DVCF file is rendered with ``--luft``, XSTRAND is placed, unmodified, as the fourth value in INFO/PRIM.

**REF**

Genozip checks whether the Luft reference matches REF\ :subscript:`prim` or one of the ALT\ :subscript:`prim` alleles. For Indels it does so by counting the number of repeats to the right of variant. It reverse-complements if XSTRAND=X.

If an Indel is reverse-complemented, it is also left-anchored, for example:

  | ``POS=10 REF=ACG ALT=A`` Original Deletion 

  | ``POS=8 REF=CGT ALT=T`` Step 1: Reverse complementing: The original ``A`` left-hand anchor base (POS=10) now appears as a right-hand anchor base ``T``. POS is shifted to point to the left-hand base 'C' (POS=8) rather than the 'T'
  
  | ``POS=7 REF=GCG ALT=G`` Step 2: Left-anchoring: The base to the left of REF (as it appears in the reference) - ``G`` in this example - is added as a new left-anchor, and the right-side anchor ``T`` is discarded. We also decrement POS to reflect the added base.

*Note*: Left-anchoring is *not* the same as left-aligning. Had the variant been left aligned, the POS would have been rolled all the way back to the beginning of the repeats of the *payload* of the Indel (the payload in the example above being ``CG``). Genozip does not left-align Indel variants when lifting because it would have been incorrect - the bases that appear before the a variant, for the specific samples in a particular VCF file, are not knowable from the VCF data alone.

If the mapped sequence in the Luft reference matches neither REF\ :subscript:`prim` nor ALT\ :subscript:`prim`, or if it matches one of the ALTs but the variant is not bi-allelic, the variant is rejected from lifting with RefMultiAltSwitchSNP or RefMultiAltSwitchIndel. This is because Genozip only supports REF⇆ALT switches of bi-allelic variants.

Complex Indels (those where both REF and at least one of the ALTs are multi-base) - are lifted only if the REF matches exactly in both references (with or without XSTRAND). 

Variants with symbolic allele ALTs (those with ALT values such as <DEL> or <INS>) are lifted only if the REF matches exactly in both references. They are rejected with XstrandSV if the mapping is with XSTRAND, and they are also rejected with INFO/END, if INFO/END does not map to the same chain file alignment as POS.

Complex Rearrangement variants (whose ALT contains the characters ``[`` or ``]``) are always rejected, with ComplexRearrangements.

Variants with a multi-base REF are lifted only if the entire REF fits into a single chain file alignment, with the exception of Deletion variants where a Deletion variant is allowed to have its anchor base on an alignment and its payload in the gap following the alignment (this would be a REF⇆ALT switch).

Variants are rejected with REFMismatchesReference if their REF does not match the Primary reference (this would be an indication of an error in the VCF file, or usage of the wrong reference file). In case the reference contains `IUPAC "bases" <http://www.bioinformatics.org/sms/iupac.html>`_ (other than A,C,T,G,N), a base is considered matching if it matches one of the IUPAC's "base" constituent bases.

When a DVCF file is rendered with ``--luft``, REF\ :subscript:`luft` is placed in the REF field, and REF\ :subscript:`prim` is placed as the third value of INFO/PRIM.

**ALT**

The ALT field is not *lifted*, rather, it is *cross-rendered* - its value is calculated by ``genocat`` based on REF in Primary and Luft coordinates, ALT in Primary coordinates and XSTRAND.

**INFO and FORMAT fields**

Like ALT, the INFO and FORMAT subfields are not *lifted*, they are only *cross-rendered*.

Genozip implements the 10 Rendering Algorithms (or *RendAlgs*) listed in the table below. They are similar to those defined in the :ref:`DVCF specification <dvcf-spec>`, with the following refinements:
  - The trigger defined as *REF Change* in the DVCF specification is implemented in Genozip only in the case of a REF⇆ALT switch in a bi-allelic variant. 
  - The A_tag RendAlg defined in the DVCF specification is implemented only for AN and is hence defined as A_AN in the table below.
  - All Genozip's RendAlgs are not only losslessly invertible (as required by the DVCF specification), but they are involutions - i.e. they are the inverse of themselves - applying them twice on the intended INFO or FORMAT value, results in getting the same value back. Hence, the same algorithm is applied to convert a field from Primary to Luft, and from Luft to Primary.

Genozip selects the algorithm to apply each field, by the RendAlg attribute of the corresponding ##INFO or ##FORMAT meta-information line in the VCF header. The default RendAlg algorithm selected for each INFO and FORMAT tag is listed in the *Applied to* column of the table below, and is based on the tag's ``ID`` and ``Number`` attributes. 

The RendAlg attribute, if it is missing from any particular ##INFO or ##FORMAT line, is added by Genozip when *lifting* a source VCF to a DVCF (using ``genozip --chain``) or when compressing a Primary-coordinates or Luft-coordinates DVCF file, according to the table below. You may set the RendAlg attribute yourself or modify the default setting of Genozip, to one of the supported RendAlgs, or turn it off, by setting it to NONE.

For FORMAT or INFO subfields that are lacking a corresponding ##INFO or ##FORMAT line, Genozip selects the default RendAlg algorithm defined in the table below.

When cross-rendering (either in ``genocat --luft`` or when executing ``genozip`` on a Luft-rendition VCF file), an INFO or FORMAT field is unmodified, unless the trigger defined for its RendAlg activated (i.e. the conditions for its activation defined in the table below apply). If the trigger is activated, the data transformation described in the Action column of the table below is applied to the field's value, to obtain the cross-rendered value.

+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| RendAlg   | Triggered upon    | Action                                     | Rejected if     | Applied to     |    
+===========+===================+============================================+=================+================+
| G         | REF⇆ALT switch    | Re-order the per-genotype values of the    | Not bi-allelic  | Number=G       |
|           |                   | subfield                                   |                 | FORMAT/GL      |
|           |                   |                                            | Ploidy > 2      | FORMAT/PL      |
|           |                   |                                            |                 | FORMAT/PRI     |
|           |                   |                                            |                 | FORMAT/GP      |
|           |                   |                                            |                 | FORMAT/PP      |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| R         | REF⇆ALT switch    | Re-order the per-allele values of the      | Not bi-allelic  | Number=R       |
|           |                   | subfield                                   |                 | FORMAT/AD      |
|           |                   |                                            |                 | FORMAT/ADF     |
|           |                   |                                            |                 | FORMAT/ADR     |
|           |                   |                                            |                 | FORMAT/ADALL   |
|           |                   |                                            |                 | FORMAT/F1R2    |
|           |                   |                                            |                 | FORMAT/F2R1    |
|           |                   |                                            |                 | FORMAT/DP_HIST |
|           |                   |                                            |                 | FORMAT/GQ_HIST |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| R2        | REF⇆ALT switch    | Switch the first 2 values with the last 2  | Not bi-allelic  | FORMAT/SB      |
|           |                   | for these 4-value fields                   |                 | FORMAT/MB      |
|           |                   |                                            |                 | FORMAT/SAC     |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| A_1       | REF⇆ALT switch    | Applied to AF-like fields:                 | Value ∉ [0,1]   | INFO/AF        |
|           |                   | new value is (1.0-value)                   |                 | INFO/AF\_\*    |
|           |                   |                                            |                 | INFO/\*_AF     |
|           |                   |                                            |                 | (not MAX_AF)   |
|           |                   |                                            |                 | INFO/MLEAF     |
|           |                   |                                            |                 | INFO/LDAF      |
|           |                   |                                            |                 | FORMAT/AF      |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| A_AN      | REF⇆ALT switch    | Applied to AC-like fields:                 | Not bi-allelic  | INFO/AC        |
|           |                   | new value is (AN-value)                    |                 | INFO/MLEAC     |
|           |                   |                                            | No INFO/AN      |                |
|           |                   |                                            |                 |                |
|           |                   |                                            | AC ∉ [0,AN]     |                |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| PLOIDY    | REF⇆ALT switch    | New value is (ploidy-value).               | Not bi-allelic  | FORMAT/DS      |
|           |                   | ploidy=number of haplotypes in GT of this  |                 |                |
|           |                   | sample.                                    | No FORMAT/GT    |                |
|           |                   |                                            |                 |                |
|           |                   |                                            | val ∉ [0,ploidy]|                |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| GT        | REF⇆ALT switch    | In a GT subfield: switch 0⇆1               | Not bi-allelic  | FORMAT/GT      |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| XREV      | XSTRAND=X         | Reverse the order of values                |                 | INFO/BaseCounts|
|           |                   |                                            |                 |                |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| END       | Always            | Recalculate the value of END so that       | POS and END not | INFO/END       |
|           |                   | END-POS remains unchanged                  | on same chain   |                |
|           |                   |                                            | file alignment  |                |
|           |                   |                                            |                 |                |
|           |                   |                                            | Strand change   |                |
|           |                   |                                            | (XSTRAND=X)     |                |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| ALLELE    | Always            | Value is identical to one of the REF or    | Value is REF,   | INFO/AA        | 
|           |                   | ALT alleles. It shall remain identical to  | and reference   |                |
|           |                   | that allele even if it changes order or    | changed but not |                |
|           |                   | is reverse-complemented and shifted        | to an ALT       |                |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+
| NONE      | Never             | Do nothing                                 |                 | All other      |
|           |                   |                                            |                 | subfields      |
+-----------+-------------------+--------------------------------------------+-----------------+----------------+

|

* Note: fields with RendAlg=A_1 or PLOIDY which contain a value in scientific notation, e.g. 2.3e-04, who's variant is a REF⇆ALT switch, are converted to standard notation during *lifting* (i.e. ``genozip --chain``), so ``2.3e-04`` is changed to ``0.00023``. This is to ensure losslessness when cross-rendering variants with a REF⇆ALT switch.
  