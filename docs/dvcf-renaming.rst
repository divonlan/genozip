.. _dvcf-renaming:

Renaming and dropping annotations in a DVCF
===========================================

.. include:: dvcf-see-also.rst


**At a glance**

In some cases, an annotation's name (rather than value) changes between the Primary and Luft renditions. This happens in case of a REF⇆ALT switch for annotations with a name that contains a reference to the REF or ALT allele (for example: ALT_F1R2), and in case of a strand reversal where the annotation name contains a reference to the strand, for example ADF. In some other cases, the annotation makes no sense in the Luft coordinates, and should be dropped entirely. Genozip implements annotation dropping by adding a "DROP\_" prefix to their name.

**Genozip default annotation renaming and dropping**

These are the annotations that are renamed by default:

============= ====== ============ =================
Annotation    Type   Renamed to   Upon
============= ====== ============ =================
MAX_AF        INFO   DROP_MAX_AF  REF⇆ALT switch 
CLNHGVS       INFO   DROP_CLNHGVS Always
ADF           FORMAT ADR          Strand reversal
ADR           FORMAT ADF          Strand reversal
F1R1          FORMAT F2R1         Strand reversal
F2R1          FORMAT F1R2         Strand reversal
REF_F1R2      FORMAT REF_F2R1     Strand reversal
                     ALT_F1R2     REF⇆ALT switch
                     ALT_F2R1     REF⇆ALT + Strand
ALT_F1R2      FORMAT ALT_F2R1     Strand reversal
                     REF_F1R2     REF⇆ALT switch
                     REF_F2R1     REF⇆ALT + Strand
REF_F2R1      FORMAT REF_F1R2     Strand reversal
                     ALT_F2R1     REF⇆ALT switch
                     ALT_F1R2     REF⇆ALT + Strand
ALT_F2R1      FORMAT ALT_F1R2     Strand reversal
                     REF_F2R1     REF⇆ALT switch
                     REF_F1R2     REF⇆ALT + Strand
============= ====== ============ =================

**The --dvcf-rename option**

Annotations may be renamed with the ``--dvcf-rename`` command line option, for example:

``--dvcf-rename=FORMAT/ADF:STRAND>ADR|REFALT>DROP_ADF``

The argument is a comma-separated list of all annotations that need to be renamed (this example contains only one annotation - ``FORMAT/ADF``):

- The annotation name (``FORMAT/ADF`` in this case) can be the name only (eg ``ADF``), or prefixed with ``INFO/`` or ``FORMAT/`` to resolve ambiguity. 

- The rules for renaming the particular annotation are specified as a \| (pipe)-seperated list to the right of the colon. In the example above, we have two rules: ``STRAND>ADR`` and ``REFALT>DROP_ADF``.

- Each rule consists of an event and a destination annotation name, seperated by a > (greater-than) character. The event can be one of the four: 

========== ============================
Rule       Rule activated upon
========== ============================
``STRAND`` Strand reversal
``REFALT`` REF⇆ALT switch
``TLAFER`` Concurrent strand reversal and REF⇆ALT switch
``ALWAYS`` Always
========== ============================

**The --dvcf-drop option**

Annotations maybe dropped with the ``--dvcf-drop`` command line option, for example:

``--dvcf-drop=INFO/MAX_AF:REFALT``

This is equivalent of:

``--dvcf-rename=INFO/MAX_AF:REFALT>DROP_MAX_AF``

To override Genozip's default renaming, just rename the tag to itself, for example:

``--dvcf-rename=INFO/MAX_AF:Always>MAX_AF``

**The --show-rename-tags option**

Shows the list of tags that are to be renamed.
