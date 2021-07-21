.. _dvcf-limitations:

Genozip DVCF lifting limitations
================================

See also:

    | :ref:`Dual-coordinate VCF files <dvcf>`
    |
    | :ref:`Rendering a DVCF <dvcf-rendering>`
    |
    | `Dual-coordinate VCF Specification <https://www.researchgate.net/publication/351904893_The_Variant_Call_Format_Dual_Coordinates_Extension_DVCF_Specification>`_
    |
    | :ref:`Chain files <dvcf-chain-files>`

Here is a list of the current limitations of Genozip DVCF lifting:

1. Cases in which variants won't be lifted and will remain primary-only, despite it being theoretically possible (yet complicated) to lift the variant:

    | a. Bi-allelic variants in which the sequence in the Luft reference genome represents an allele that is neither the REF allele nor the ALT one. However, if this variant is a SNP in which AF=1 or AC=AN the variant is nevertheless lifted.
    |
    | b. Variants with more than 2 alleles (i.e. more than one ALT allele), in which the sequence in the Luft reference and Primary reference differ.
    |
    | c. Complex re-arrangement variants (as defined in section 5.4 of the `VCF specification <https://samtools.github.io/hts-specs/VCFv4.3.pdf>`_).
    | 
    | d. Variants with a symbolic ALT allele (eg <DEL>) where the Luft reference and Primary reference differ and/or there is a strand reversal.
    | 
    | e. An indel variant that is not left-anchored (example: in "A ACCT", 'A' is the *anchor base* and it appears on the left, hence the variant is left-anchored, while "G CTTG" is not left-anchored).
    | 
    | Note: When running ``genozip --chain``, a human-readable *rejects file*, eg ``mydata.d.vcf.genozip.rejects.txt`` will describe the cause of each rejected variant.
    
2. Genozip does not yet implement the ``MAX_tag`` *RendAlg* defined in the `DVCF Specification <https://www.researchgate.net/publication/351904893_The_Variant_Call_Format_Dual_Coordinates_Extension_DVCF_Specification>`_, so it can't cross-render INFO/MAX_AF

3. Genozip does not yet implement *Tag Renaming* defined in the `DVCF Specification <https://www.researchgate.net/publication/351904893_The_Variant_Call_Format_Dual_Coordinates_Extension_DVCF_Specification>`_, so can't get handle tag switching (for example ``ADF`` ⇆ ``ADR``) in case of strand reversal.


    