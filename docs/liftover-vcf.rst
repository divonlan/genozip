.. _liftover_vcf:

How Genozip lifts-over a VCF file
=================================

As described in :ref:`dual_coordinates`, a dual coordinates VCF file may be viewed in two ways: by default, ``genocat`` shows it the Primary coordinates, and when using ``genocat --luft``, it is shown in the "Luft" coordinates. In case you're wondering, Luft is our made-up alternative past paticiple of the verb Lift, similar to "swim"->"swum", "begin"->"begun".

The Primary and Luft views of the file contain the same infomation, and maybe losslessly converted back and forth. Here we describe. We call "Lifting over" the processing of converting a Primary coordinates file to the Luft coordinates, and we call "Lifting back" the inverse process.

**Lift over of the VCF header**


**Cases where Genozip rejects lifting over a variant**

**Handling INFO and FORMAT fields in case of a REF <> ALT switch**

*Format GT*



