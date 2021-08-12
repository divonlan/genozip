.. _dvcf-limitations:

Genozip DVCF lifting limitations
================================

.. include:: dvcf-see-also.rst

Cases in which variants won't be lifted and will remain primary-only, despite it being theoretically possible (yet complicated) to lift the variant:

| 1. Bi-allelic variants in which the sequence in the Luft reference genome represents an allele that is neither the REF allele nor the ALT one. However, if this variant is a SNP in which AF=1 or AC=AN the variant is nevertheless lifted.
|
| 2. Variants with more than 2 alleles (i.e. more than one ALT allele), in which the sequence in the Luft reference and Primary reference differ.
|
| 3. Complex re-arrangement variants (as defined in section 5.4 of the `VCF specification <https://samtools.github.io/hts-specs/VCFv4.3.pdf>`_).
| 
| 4. Variants with a symbolic ALT allele (eg <DEL>) where the Luft reference and Primary reference differ and/or there is a strand reversal.
| 
| 5. An indel variant that is not left-anchored and has a strand reversal ("left-anchored" means for example "A ACCT" - 'A' is the *anchor base* and it appears on the left, while "G CTTG" is not left-anchored).
| 
| Note: When running ``genozip --chain``, a human-readable *rejects file*, eg ``mydata.d.vcf.genozip.rejects.txt`` will describe the cause of each rejected variant.
    


    