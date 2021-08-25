.. _gff3:

Genozip GFF3 files
==================

**Compressing a GFF3 file**

::

    $ genozip myfile.gff3 
    genozip myfile.gff3 : Done (4 seconds, GFF3 compression ratio: 6.6)
    
    $ ls -lh myfile.gff3*
    -rwxrwxrwx 1 divon divon  26M Aug 22 22:48 myfile.gff3
    -rwxrwxrwx 1 divon divon 3.9M Aug 22 22:49 myfile.gff3.genozip

This creates a compressed file, without modifying the original file. This also works with ``.gff3.gz``, ``.gff3.bz2``, ``.gff3.xz`` as well as ``.gff`` and ``.gvf`` in all these combinations.

Some useful command line options (for a full list, see :ref:`genozip manual<genozip>`):

``genozip --test myfile.gff3``: after completing the compression, the file is uncompressed in memory, and its MD5 is compared to that of the original file.

``genozip --replace myfile.gff3``: the original file is removed after successful compression

**Compressing multiple files into a tar archive**

``genozip *.gff3 --tar mydata.tar``. See details: :ref:`archiving`.

**Optimizing compression**

These are options that modify the file in ways that improve compression. ``--optimize`` is an umbrella option that activates all optimization options.

``genozip --optimize-sort myfile.gff3.gz``

Sorts ATTR subfields alphabetically.

``genozip --optimize-Vf myfile.gff3.gz``

The value of Variant_freq is rounded to 2 significant digits

The option ``--stats`` can be used in ``genozip``, ``genounzip`` or ``genocat`` to get a better understanding of the information content of the file. For example:
   
::

    $ genocat --stats myfile.gff3.genozip

    GFF3 file: myfile.gff3
    Sequences: 40,849   Dictionaries: 55   Vblocks: 1 x 25 MB  Sections: 75
    Genozip version: 12.0.30 github
    Date compressed: 2021-08-22 22:49:04 ACDT
    License v12.0.11 granted to: ***** accepted by: ***** on 2021-07-23 14:33:51 ACDT from IP=*****
    
    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    Gap                     2.1 MB  55.4%    9.2 MB  36.2%    4.3X
    ID                    734.5 KB  18.6%    1.4 MB   5.5%    2.0X
    Target                192.6 KB   4.9%    1.2 MB   4.6%    6.2X
    START                 119.0 KB   3.0%  350.4 KB   1.3%    2.9X
    pct_coverage          114.6 KB   2.9%  360.3 KB   1.4%    3.1X
    pct_coverage_hiqual   114.6 KB   2.9%  360.3 KB   1.4%    3.1X
    END                    90.0 KB   2.3%  350.5 KB   1.3%    3.9X
    num_ident              82.0 KB   2.1%  157.1 KB   0.6%    1.9X
    pct_identity_gap       78.1 KB   2.0%  252.0 KB   1.0%    3.2X
    pct_identity_gapopen   75.5 KB   1.9%  251.9 KB   1.0%    3.3X
    pct_identity_ungap     73.5 KB   1.9%  250.3 KB   1.0%    3.4X
    num_mismatch           40.5 KB   1.0%   80.2 KB   0.3%    2.0X
    gap_count              27.6 KB   0.7%   60.6 KB   0.2%    2.2X
    ATTRS                   8.6 KB   0.2%    9.8 MB  38.6% 1169.3X
    pct_ident_quantized     5.7 KB   0.1%   79.3 KB   0.3%   14.0X
    reciprocity             2.4 KB   0.1%   39.9 KB   0.2%   16.5X
    same_unit_reciprocit    1.7 KB   0.0%   31.4 KB   0.1%   18.9X
    STRAND                  1.3 KB   0.0%   79.8 KB   0.3%   59.2X
    Other                    587 B   0.0%         -   0.0%    0.0X
    TXT_HEADER               417 B   0.0%      69 B   0.0%    0.2X
    SEQID                    320 B   0.0%  438.8 KB   1.7% 1404.2X
    SOURCE                    47 B   0.0%  279.2 KB   1.1% 6083.9X
    TYPE                      46 B   0.0%  239.3 KB   0.9% 5328.1X
    SCORE                     42 B   0.0%   79.8 KB   0.3% 1945.2X
    PHASE                     42 B   0.0%   79.8 KB   0.3% 1945.2X
    chromosome_to_chromo      42 B   0.0%   21.3 KB   0.1%  518.5X
    mitochondrion_to_mit      42 B   0.0%       3 B   0.0%    0.1X
    best_on_query_same_u      42 B   0.0%    9.2 KB   0.0%  225.1X
    best_on_subject           42 B   0.0%    5.8 KB   0.0%  140.9X
    genomic_to_genomic        42 B   0.0%   39.9 KB   0.2%  972.5X
    best_on_query             42 B   0.0%    6.3 KB   0.0%  154.7X
    salvage                   42 B   0.0%    3.4 KB   0.0%   82.7X
    chromosome_to_alt         42 B   0.0%    5.8 KB   0.0%  140.4X
    fix_patch_align           42 B   0.0%    1.2 KB   0.0%   28.5X
    COMMENT                   41 B   0.0%         -   0.0%    0.0X
    TOTAL                   3.9 MB 100.0%   25.5 MB 100.0%    6.6X

**Uncompressing**

``genounzip myfile.gff3.genozip``

Uncompresses a file

``genocat myfile.gff3.genozip``

Uncompresses a file into stdout (i.e. the terminal).

``genocat --bgzf 6 myfile.gff3.genozip`` 
``genounzip --bgzf 6 myfile.gff3.genozip`` 

Sets the level BGZF compression (for .gff3.gz output format) - from 0 (no compression) to 12 (best yet slowest compression). Absent this option, ``genounzip`` attemps to recover the BGZF compression level of the original file, while ``genocat`` uncompresses without BGZF compression. 
    
**Using in a pipeline**

| Compressing piped input: 
| ``my-pipeline | genozip - --input gff3 --output myfile.gff3.genozip`` 

| Uncompressing to a pipe: 
| ``genocat myfile.gff3.genozip | my-pipeline``

**Grepping**

``genocat --grep-w rs1357314184 myfile.gff3.genozip`` 

Displays the lines containing "rs1357314184" (strings that match exactly).

``genocat --grep Dbxref=dbSNP_152:rs myfile.gff3.genozip`` 

Displays the lines containing "Dbxref=dbSNP_152:rs" (possibly a substring of a longer string).

**Filtering specific regions of the genome**

Examples of using ``--regions`` (or its shortcut ``-r``):

=============================================== =============================================
``genocat myfile.gff3.genozip -r 22:1000-2000`` Positions 1000 to 2000 on contig 22
``genocat myfile.gff3.genozip -r 22:1000+151``  151 bases, starting pos 1000, on contig 22
``genocat myfile.gff3.genozip -r -2000,2500-``  Two ranges on all contigs
``genocat myfile.gff3.genozip -r chr21,chr22``  Contigs chr21 and chr22 in their entirety
``genocat myfile.gff3.genozip -r ^MT,Y``        All contigs, excluding MT and Y
``genocat myfile.gff3.genozip -r ^-1000``       All contigs, excluding positions up to 1000
``genocat myfile.gff3.genozip -r chrM``         Contig chrM
=============================================== =============================================

``genocat --regions-file <filename> myfile.gff3.genozip`` 

Get regions from a tab-separated file. An example of a valid file:

::

   chr22	17000000	17000099
   chr22	17000000	+100
   chr22	17000000

**Multi-threading**

By default, Genozip attempts to utilize as many cores as available. For that, it sets the number of threads to be a bit more than the number of cores (a practice known as "over-subscription"), as at any given moment some threads might be idle, waiting for a resource to become available. The ``--threads <number>`` option allows explicit specification of the number of "compute threads" to be used (in addition a small number of I/O threads is used too, usually 1 or 2).

**Memory (RAM) consumption**

In ``genozip``, each compute thread is assigned a segment of the input file, known as a VBlock. By default, the size of the VBlock is set automatically to balance memory consumption and compression ratio for the particular input file, however it may be set explicitly with ``genozip --vblock <megabytes>`` (<megabytes> is an integer between 1 and 2048). A larger VBlock usually results in better compression while a smaller VBlock causes ``genozip`` to consume less RAM. The VBlock size can be observed at the top of the ``--stats`` report. ``genozip``'s memory consumption is linear with (VBlock-size X number-of-threads). 

``genocat`` and ``genounzip`` also consume memory linearly with (VBlock-size X number-of-threads), where VBlock-size is the value used by ``genozip`` of the particular file (it cannot be modified ``genocat`` or ``genounzip``). Usually, ``genocat`` and ``genounzip`` consume significantly less memory compared to ``genozip``.

