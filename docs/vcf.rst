.. _vcf:

Genozip VCF files
=================

**Compressing a VCF or BCF file**

::

    $ genozip myfile.vcf.gz
    genozip myfile.vcf.gz : Done (2 seconds, VCF compression ratio: 21.9 - better than .vcf.gz by a factor of 2.2)

    $ ls -lh myfile.vcf*
    -rwxrwxrwx 1 divon divon 1.9M Aug 22 00:15 myfile.vcf.genozip
    -rwxrwxrwx 1 divon divon 4.0M Aug 22 00:14 myfile.vcf.gz

This creates a compressed file, without modifying the original file. This also works with ``.vcf.gz``, ``.vcf.bz2``, ``.vcf.xz`` and ``.bcf``.

Some useful command line options (for a full list, see :ref:`genozip manual<genozip>`):

``genozip --test myfile.vcf``: after completing the compression, the file is uncompressed in memory, and its MD5 is compared to that of the original file.

``genozip --replace myfile.vcf``: the original file is removed after successful compression

**Compressing multiple files into a tar archive**

``genozip *.vcf --tar mydata.tar``. See details: :ref:`archiving`.

**Optimizing compression**

These are options that modify the file in ways that improve compression. ``--optimize`` is an umbrella option that activates all optimization options.

``genozip --optimize-sort myfile.vcf.gz``

Sorts INFO subfields alphabetically.

``genozip --optimize-phred myfile.vcf.gz``

Applied to FORMAT/PL FORMAT/PRI FORMAT/PP and (VCF v4.2 or earlier) FORMAT/GL - Phred scores are rounded to the nearest integer and capped at 60.

``genozip --GL-to-PL myfile.vcf.gz``

The FORMAT/GL field is converted to PL and Phred values are capped at 60.

``genozip --GP-to-PP myfile.vcf.gz``

Applicable to VCF v4.3 and later: The FORMAT/GP field is converted to PP and Phred values are capped at 60.

``genozip --optimize-VQSLOD myfile.vcf.gz``

The VQSLOD value is rounded to 2 significant digits. 

``genozip --reference reference-file.ref.genozip myfile.vcf.gz``

Compresses against a reference. This option is *not* included in ``--optimize``. It improves compression in files in which the "REF+ALT" field consumes significant part of the GENOZIP content (see ``--stats`` below).

The option ``--stats`` can be used in ``genozip``, ``genounzip`` or ``genocat`` to get a better understanding of the information content of the file. For example:
   
::

    $ genocat --stats myfile.vcf.genozip

    VCF file: myfile.vcf.gz
    Samples: 1211   Variants: 1,500   Dictionaries: 93   Vblocks: 3 x 16 MB  Sections: 210
    Genozip version: 12.0.30 github
    Date compressed: 2021-08-22 00:15:08 ACDT
    License v12.0.11 granted to: ***** accepted by: ***** on 2021-07-23 14:33:51 ACDT from IP=*****

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    PL                      1.1 MB  59.9%   10.7 MB  26.8%    9.8X
    GQ                    288.4 KB  15.4%    2.3 MB   5.6%    8.0X
    AD                    281.7 KB  15.1%    5.5 MB  13.6%   19.9X
    GT                     70.7 KB   3.8%    5.2 MB  13.0%   75.3X
    TXT_HEADER             21.0 KB   1.1%  145.6 KB   0.4%    6.9X
    PID                    17.4 KB   0.9%    1.4 MB   3.4%   80.6X
    PGT                    10.3 KB   0.6%    1.3 MB   3.3%  130.4X
    QUAL                    4.4 KB   0.2%    9.6 KB   0.0%    2.2X
    InbreedingCoeff         4.2 KB   0.2%    9.7 KB   0.0%    2.3X
    AF                      4.2 KB   0.2%   11.0 KB   0.0%    2.6X
    MQ                      4.0 KB   0.2%    7.0 KB   0.0%    1.8X
    MLEAF                   3.9 KB   0.2%   10.0 KB   0.0%    2.6X
    QD                      3.8 KB   0.2%    6.4 KB   0.0%    1.7X
    ExcessHet               3.7 KB   0.2%    7.8 KB   0.0%    2.1X
    DP                      3.6 KB   0.2%    5.9 KB   0.0%    1.6X
    AN                      3.5 KB   0.2%    5.6 KB   0.0%    1.6X
    SOR                     3.4 KB   0.2%    7.2 KB   0.0%    2.1X
    FS                      2.9 KB   0.2%    4.8 KB   0.0%    1.6X
    MQRankSum               2.8 KB   0.2%    5.8 KB   0.0%    2.1X
    BaseQRankSum            2.8 KB   0.1%    5.5 KB   0.0%    2.0X
    ReadPosRankSum          2.7 KB   0.1%    5.6 KB   0.0%    2.0X
    MLEAC                   1.9 KB   0.1%    2.1 KB   0.0%    1.1X
    POS                     1.6 KB   0.1%    8.7 KB   0.0%    5.4X
    DP                      1.5 KB   0.1%    2.0 MB   5.0% 1338.7X
    Other                   1.2 KB   0.1%   11.3 MB  28.1% 9658.3X
    REF+ALT                 1.1 KB   0.1%    5.9 KB   0.0%    5.3X
    CHROM                    915 B   0.0%    2.9 KB   0.0%    3.3X
    INFO                     729 B   0.0%  179.4 KB   0.4%  252.0X
    AC                       555 B   0.0%    1.9 KB   0.0%    3.4X
    FORMAT                   526 B   0.0%   30.7 KB   0.1%   59.9X
    COORDS                   391 B   0.0%         -   0.0%    0.0X
    BGZF                      56 B   0.0%         -   0.0%    0.0X
    oXSTRAND                  44 B   0.0%         -   0.0%    0.0X
    ID                        42 B   0.0%    2.9 KB   0.0%   71.4X
    FILTER                    42 B   0.0%    2.9 KB   0.0%   71.4X
    ClippingRankSum           42 B   0.0%    1.2 KB   0.0%   29.4X
    GENOZIP vs BGZF         1.8 MB 100.0%    4.0 MB 100.0%    2.2X
    GENOZIP vs TXT          1.8 MB 100.0%   40.1 MB 100.0%   21.9X

In this paritcular example, we observe that the PL field consumes a whopping 59.9% of the total compressed file size. Therefore, we can expect that ``--optimize-phred`` will significantly reduce the compressed file size. In contrast, the REF+ALT field, in this case, consumes only 0.1% of the compressed file size. Therefore, we can expect that using ``--reference`` will *not* significantly reduce the compressed file size.

**Uncompressing**

``genounzip myfile.vcf.genozip``

Uncompresses a file

``genocat myfile.vcf.genozip``

Uncompresses a file into stdout (i.e. the terminal).

``genounzip --index myfile.vcf.genozip``

Uncompresses a file and also generates a CSI index file, using `bcftools index <http://samtools.github.io/bcftools/bcftools.html#index>`_. bcftools needs to be installed for this option to work. 

| ``genocat --bgzf 6 myfile.vcf.genozip``
| ``genounzip --bgzf 6 myfile.vcf.genozip`` 

Sets the level BGZF compression (for .vcf.gz output format) - from 0 (no compression) to 12 (best yet slowest compression). Absent this option, ``genounzip`` attemps to recover the BGZF compression level of the original file, while ``genocat`` uncompresses without BGZF compression. 
    
**Using in a pipeline**

| Compressing piped input: 
| ``my-pipeline | genozip - --input vcf --output myfile.vcf.genozip`` 

| Uncompressing to a pipe: 
| ``genocat myfile.vcf.genozip | my-pipeline``

**Downsampling** 

``genocat --downsample 10,0 myfile.vcf.genozip`` 

Displays only the first (#0) variant in every 10 variants.

**Grepping**

``genocat --grep-w AC=2 myfile.vcf.genozip`` 

Displays the variants containing "AC=2" (strings that match exactly).

``genocat --grep ACCTTAAT myfile.vcf.genozip`` 

Displays the variants containing "ACCTTAAT" (possibly a substring of a longer string).

**Selecting samples**

``genocat myfile.vcf.genozip --samples HG00255,HG00256``  

Shows two samples.

``genocat myfile.vcf.genozip --samples ^HG00255,HG00256`` 

Shows all samples except these two.

``genocat myfile.vcf.genozip --samples 5``                

Shows the first 5 samples.

``genocat myfile.vcf.genozip --drop-genotypes``                

Drops all samples and the FORMAT columns. ``--drop-genotypes`` is the same as ``--samples 0``, ``-s 0`` and ``-G``.

``genocat myfile.vcf.genozip --GT-only``

Within samples, outputs only genotype (GT) data - dropping the other subfields.

**SNPs or indels only**

``genocat myfile.vcf.genozip --snps-only``

Drops variants that are not a Single Nucleotide Polymorphism (SNP).

``genocat myfile.vcf.genozip --indels-only``

Drops variants that are not Insertions or Deletions (indel).

**The VCF header**

``genocat --header-only myfile.vcf.genozip``

Displays only the VCF header.

``genocat --no-header myfile.vcf.genozip`` 

Displays the file without the VCF header.

``genocat --header-one myfile.vcf.genozip`` 

Displays the file without the VCF header, except for the #CHROM line.

``genocat --no-PG myfile.vcf.genozip`` 

When modifying the data in a file using genocat, Genozip normally adds a "##genozip_command" line to the VCF header. With this option it doesn't.

**Filtering specific regions of the genome**

Examples of using ``--regions`` (or its shortcut ``-r``):

============================================== =============================================
``genocat myfile.vcf.genozip -r 22:1000-2000`` Positions 1000 to 2000 on contig 22
``genocat myfile.vcf.genozip -r 22:1000+151``  151 bases, starting pos 1000, on contig 22
``genocat myfile.vcf.genozip -r -2000,2500-``  Two ranges on all contigs
``genocat myfile.vcf.genozip -r chr21,chr22``  Contigs chr21 and chr22 in their entirety
``genocat myfile.vcf.genozip -r ^MT,Y``        All contigs, excluding MT and Y
``genocat myfile.vcf.genozip -r ^-1000``       All contigs, excluding positions up to 1000
``genocat myfile.vcf.genozip -r chrM``         Contig chrM
============================================== =============================================

``genocat --regions-file <filename> myfile.vcf.genozip`` 

Get regions from a tab-separated file. An example of a valid file:

::

   chr22	17000000	17000099
   chr22	17000000	+100
   chr22	17000000

**Sorting**

``genozip --sort myfile.vcf``

Variants are sorted by CHROM and POS. This works for "mildly unsorted" files. This is the default with ``--chain`` is used, unless ``--unsorted`` is specified.

``genocat --unsorted myfile.vcf.genozip``

Shows the variants in their original order.

**Adding line numbers**

``genozip --add-line-numbers myfile.vcf``

Replaces the ID field in each variant with a sequential line number starting from 1.

**Flat coordinates (GPOS)**

``genocat --gpos --reference reference-file.ref.genozip myfile.vcf.genozip`` 

Replaces (CHROM,POS) with a coordinate in GPOS (Global POSition) terms. GPOS is a single genome-wide coordinate defined by a reference file, in which contigs appear in the order of the original FASTA data used to generate the reference file. 

``genocat --show-ref-contigs reference-file.ref.genozip``

Shows the mapping of CHROM to GPOS.

**BCF files**

Genozip does not support BCF natively - it uses `bcftools <http://samtools.github.io/bcftools/bcftools.html>`_ to convert BCF files to/from the VCF format, and as such it requires bcftools to be installed for the BCF features to work.

``genozip myfile.bcf``

Compresses a BCF file.

``genocat --bcf myfile.vcf.genozip`` 

Outputs the file in BCF format.

**Dual-coordinate VCF files**

Genozip has the unique ability to represent a VCF file with coordinates in two different reference genomes concurrently. See :ref:`dvcf`.

``genozip --chain mychainfile.chain.genozip myfile.vcf``

Lifts a VCF file to a dual-coordinate VCF (DVCF) - this generates ``myfile.d.vcf.genozip``.

| When using ``--chain``, additional options may be combined:
| - ``--dvcf-rename``, ``--dvcf-drop`` - specify annotations that should be renamed or dropped when cross rendering Primary➝Luft or Luft➝Primary. See :ref:`dvcf-renaming`.
| - ``--show-rename-tags`` -  shows tags that are to be renamed. Used when compressing a DVCF or in combination with --chain.
| - ``--show-lifts`` - output successful lifts to the rejects file too, not only rejected lifts.
| - ``--show-counts=o\$TATUS``, ``--show-counts=COORDS`` - see below
| - ``--show-chain`` - displays all chain file alignments.


``genocat myfile.d.vcf.genozip``

Displays the file in the *Primary* coordinates.

``genocat --luft myfile.d.vcf.genozip``

Displays the file in the *Luft* coordinates.

``genocat --show-ostatus myfile.d.vcf.genozip``

Adds oSTATUS to the INFO field - the status of the variant relative to the lift process. 

``genocat --show-counts=o\$TATUS myfile.d.vcf.genozip``

Shows summary statistics of variant lift outcome (also works with ``genozip --chain``).

``genocat --show-counts=COORDS myfile.d.vcf.genozip``

Shows summary statistics of variant coordinates (also works with ``genozip --chain``).

``genocat --show-dvcf myfile.d.vcf.genozip``

For each variant, shows its coordinate system (Primary or Luft or Both) and its oStatus. May be used with or without --luft.

**Multi-threading**

By default, Genozip attempts to utilize as many cores as available. For that, it sets the number of threads to be a bit more than the number of cores (a practice known as "over-subscription"), as at any given moment some threads might be idle, waiting for a resource to become available. The ``--threads <number>`` option allows explicit specification of the number of "compute threads" to be used (in addition a small number of I/O threads is used too, usually 1 or 2).

**Memory (RAM) consumption**

In ``genozip``, each compute thread is assigned a segment of the input file, known as a VBlock. By default, the size of the VBlock is set automatically to balance memory consumption and compression ratio for the particular input file, however it may be set explicitly with ``genozip --vblock <megabytes>`` (<megabytes> is an integer between 1 and 2048). A larger VBlock usually results in better compression while a smaller VBlock causes ``genozip`` to consume less RAM. The VBlock size can be observed at the top of the ``--stats`` report. ``genozip``'s memory consumption is linear with (VBlock-size X number-of-threads). 

``genocat`` and ``genounzip`` also consume memory linearly with (VBlock-size X number-of-threads), where VBlock-size is the value used by ``genozip`` of the particular file (it cannot be modified ``genocat`` or ``genounzip``). Usually, ``genocat`` and ``genounzip`` consume significantly less memory compared to ``genozip``.

