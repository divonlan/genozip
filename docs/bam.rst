.. _bam:

Genozip BAM, SAM or CRAM files
==============================

**Compressing a SAM or BAM file**

::

    $ genozip myfile.bam
    genozip myfile.bam : Done (2 seconds, BAM compression ratio: 2.2)    
    
    $ ls -lh myfile.bam*
    -rwxrwxrwx 1 divon divon 5.7M Aug 20 13:59 myfile.bam
    -rwxrwxrwx 1 divon divon 2.6M Aug 20 13:59 myfile.bam.genozip

This creates a compressed file, without modifying the original file. This also works with ``.sam``, ``.sam.gz``, ``.sam.bz2`` and ``.sam.xz`` files.

Some useful command line options (for a full list, see :ref:`genozip manual<genozip>`):

``genozip --test myfile.bam``: after completing the compression, the file is uncompressed in memory, and its MD5 is compared to that of the original file.

``genozip --replace myfile.bam``: the original file is removed after successful compression

``genozip --optimize-QUAL myfile.bam``: the QUAL field (base quality scores) is optimized to improve compression. See :ref:`genozip manual<genozip>` for details.

**Compressing multiple files into a tar archive**

``genozip *.bam --tar mydata.tar``. See details: :ref:`archiving`.

**Uncompressing**

``genounzip myfile.bam.genozip``

Uncompresses a file.

``genocat myfile.bam.genozip``

Uncompresses a file into stdout (i.e. the terminal).

``genounzip --index myfile.bam.genozip``

Uncompresses a file and also generates a BAI index file, using `samtools index <http://www.htslib.org/doc/samtools-index.html>`_. samtools needs to be installed for this option to work. 

``genounzip --output newname.sam.gz myfile.sam.genozip``
``genocat --output newname.bam myfile.bam.genozip``

Uncompresses to a particular name. The file format produced depends on whether the output file name extension is ``.bam``, ``.sam`` or ``.sam.gz``. Converting between SAM and BAM is possible with ``genocat`` but not ``genounzip``.

``genocat --sam myfile.bam.genozip`` 
``genocat --bam myfile.sam.genozip`` 

Specifies a SAM or BAM file output (implicit if ``--output`` is used).

``genocat --bgzf 6 myfile.bam.genozip`` 
``genounzip --bgzf 6 myfile.bam.genozip`` 

Sets the level BGZF compression (for .bam and .sam.gz output format) - from 0 (no compression) to 12 (best yet slowest compression). Absent this option, ``genounzip`` attemps to recover the BGZF compression level of the original file, while ``genocat`` uncompresses without BGZF compression for SAM and with BGZF compression level 1 for BAM (this is because some popular bioinforamtics tools error on BAM files with BGZF compression level 0). 
    
**Using in a pipeline**

| Compressing piped input: (use ``--input sam`` for a SAM file)
| ``my-pipeline | genozip - --input bam --output myfile.bam.genozip`` 

| Uncompressing to a pipe: 
| ``genocat myfile.bam.genozip | my-pipeline``

**Using a reference file**

Compressing against a reference file improves compression:

::

    $ genozip --make-reference hs37d5.fa.gz  # one-time preparation of the reference. should be similar to the reference used to create the BAM file
    $ genozip myfile.bam --reference hs37d5.ref.genozip
    
Note: the reference file is also needed for uncompressing. Alternatively, use ``--REFERENCE`` to copy the relevant parts of the reference file to myfile.bam.genozip, so that the reference file is not needed for decompression. Compression with a reference file works regardless of whether the SAM/BAM file is aligned.

**Compressing CRAM files**

Genozip does not support CRAM natively - it uses `samtools <https://en.wikipedia.org/wiki/SAMtools>`_ to convert CRAM files to/from the SAM format, and as such it requires samtools to be installed for CRAM compression to work.

``genozip --reference hs37d5.ref.genozip myfile.cram``

Compresses a CRAM file into ``myfile.sam.genozip`` in this case. 

Note: use of a reference is mandatory when compressing a CRAM file.

Note: When decompressing the file, it is decompressed into SAM or BAM format. Genozip does not support decompressing to CRAM format.

**Compressing Ion Torrent BAM files**

Genozip offers the ``--optimize-ZM`` option to optimize the ZM (flow signal) data for improved compression. See :ref:`genozip manual<genozip>` for details.

Example:

::

    $ wget ftp://ftp-trace.ncbi.nih.gov/HG002_NA24385_SRR1767406_IonXpress_020_rawlib_24028.30k.b37.bam

    $ genozip IonXpress_020_rawlib.b37.bam

    $ genozip IonXpress_020_rawlib.b37.bam --optimize-ZM -o IonXpress_020_rawlib.b37.optimized.bam.genozip

    $ ls -Ggh IonXpress_020_rawlib.b37*
    -rw-rw-r--+ 1 26G Aug 13 23:53 IonXpress_020_rawlib.b37.bam
    -rw-rw-r--+ 1 17G Aug 14 00:10 IonXpress_020_rawlib.b37.bam.genozip
    -rw-rw-r--+ 1 12G Aug 14 00:17 IonXpress_020_rawlib.b37.optimized.bam.genozip


**Converting to a FASTQ** 

``genocat --fastq myfile.bam.genozip`` or ``genozip --fastq=all myfile.bam.genozip`` may be used to output the data in FASTQ format. See :ref:`sam2fq` for details.

**Downsampling** 

``genocat --downsample 10,0 myfile.bam.genozip`` 

Displays only the first (#0) read in every 10 reads.

**Grepping**

``genocat --grep-w MC:Z:151M myfile.bam.genozip`` 

Displays the lines containing "MC:Z:151M" (strings that match exactly).

``genocat --grep ACCTTAAT myfile.bam.genozip`` 

Displays the lines containing "ACCTTAAT" (possibly a substring of a longer string).

**The SAM header**

``genocat --header-only myfile.bam.genozip``

Displays only the SAM header.

``genocat --no-header myfile.bam.genozip`` 

Displays the file without the SAM header.

``genocat --no-PG myfile.bam.genozip`` 

When modifying the data in a file using genocat Genozip normally adds a @PG line to the header with information about the modification. With this option it doesn't.

**Filtering specific regions of the genome**

Examples of using ``--regions`` (or its shortcut ``-r``):

============================================== =============================================
``genocat myfile.bam.genozip -r 22:1000-2000`` Positions 1000 to 2000 on contig 22
``genocat myfile.sam.genozip -r 22:1000+151``  151 bases, starting pos 1000, on contig 22
``genocat myfile.bam.genozip -r -2000,2500-``  Two ranges on all contigs
``genocat myfile.sam.genozip -r chr21,chr22``  Contigs chr21 and chr22 in their entirety
``genocat myfile.bam.genozip -r ^MT,Y``        All contigs, excluding MT and Y
``genocat myfile.bam.genozip -r ^-1000``       All contigs, excluding positions up to 1000
``genocat myfile.bam.genozip -r chrM``         Contig chrM
============================================== =============================================

``genocat --regions-file <filename>`` 

Get regions from a tab-separated file. An example of a valid file:

::

   chr22	17000000	17000099
   chr22	17000000	+100
   chr22	17000000

**Filtering reads based on FLAG**

``genocat --FLAG *{+-^}value* myfile.bam.genozip``.  Filter lines based on the FLAG value: <value> is a decimal or hexadecimal value and should be prefixed by + - or ^: 

    ==  =======================================================================
    \+  INCLUDES lines in which ALL flags in *value* are set in the line's FLAG
    \-  INCLUDES lines in which NO flags in *value* are set in the line's FLAG
    ^   EXCLUDES lines in which ALL flags in *value* are set in the line's FLAG
    ==  =======================================================================

*Example*: ``genocat --FLAG -192`` includes only lines in which neither FLAG 64 nor 128 are set. This can also be expressed as ``--FLAG -0xC0``

The FLAGs are defined in the `SAM specification <https://samtools.github.io/hts-specs/SAMv1.pdf>`_ as follows:

    ======= ===== =================================================================== 
    Decimal Hex   Meaning
    ======= ===== =================================================================== 
    1       0x1   template having multiple segments in sequencing
    2       0x2   each segment properly aligned according to the aligner
    4       0x4   segment unmapped
    8       0x8   next segment in the template unmapped
    16      0x10  SEQ being reverse complemented
    32      0x20  SEQ of the next segment in the template being reverse complemented
    64      0x40  the first segment in the template
    128     0x80  the last segment in the template
    256     0x100 secondary alignment
    512     0x200 not passing filters, such as platform/vendor quality controls
    1024    0x400 PCR or optical duplicate
    2048    0x800 supplementary alignment
    ======= ===== =================================================================== 
 
    |

**Filtering reads based on MAPQ**

``genocat --MAPQ [^]value myfile.bam.genozip`` 

Filters lines based on the MAPQ value: INCLUDE (or EXCLUDE if *value* is prefixed with ^) lines with a MAPQ greater or equal to *value*. 

**Filtering non-ACTGN "bases"**

``genocat --bases ACGTN myfile.bam.genozip``  

Displays only lines in which all characters of the SEQ are one of A,C,G,T,N

``genocat --bases ^ACGTN myfile.bam.genozip`` 

Displays only lines in which NOT all characters of the SEQ are one of A,C,G,T,N

Note: In all lines missing a sequence (i.e. SEQ=*) are included in positive --bases filters (the first example above) and excluded in negative ones.

Note: The list of IUPAC chacacters can be found here: `IUPAC codes <https://www.bioinformatics.org/sms/iupac.html>`_

**Filtering reads by species**

Genozip has the unique ability to filter SAM/BAM files by species (taxonomy id). This is useful, for example, for filtering out bacterial contamination by directly removing reads that map to bacterial genomes rather than just removing reads with low mapping quality, assuming they represent contamination. See :ref:`kraken`.

**Inspecting field-level compression statistics**

If optimizing the compressed file size is important, the option ``--stats`` can be used in ``genozip``, ``genounzip`` or ``genocat`` to get a better understanding of the information content of the individual fields. For example:
   
::

    $ genocat --stats myfile.bam.genozip
    
    BAM file: myfile.bam
    Alignment lines: 99,909   Dictionaries: 50   Vblocks: 2 x 16 MB  Sections: 143
    Genozip version: 12.0.11 conda
    Date compressed: 2021-08-20 17:28:44 ACDT
    License v12.0.11 granted to: ***** accepted by: ***** on 2021-07-23 14:33:51 ACDT from IP=*****
    
    Sections (sorted by % of genozip file):
    NAME              GENOZIP      %      TXT       %   RATIO
    QUAL             978.8 KB  38.1%   14.1 MB  45.1%   14.8X
    QNAME            606.2 KB  23.6%    3.8 MB  12.0%    6.4X
    SEQ              357.4 KB  13.9%    7.5 MB  23.9%   21.4X
    MD:Z             131.9 KB   5.1%  598.7 KB   1.9%    4.5X
    TLEN             122.1 KB   4.8%  390.3 KB   1.2%    3.2X
    PNEXT            118.8 KB   4.6%  390.3 KB   1.2%    3.3X
    XS:i              46.5 KB   1.8%   97.4 KB   0.3%    2.1X
    CIGAR             43.8 KB   1.7%  646.3 KB   2.0%   14.7X
    POS               40.1 KB   1.6%  390.3 KB   1.2%    9.7X
    FLAG              29.6 KB   1.2%  195.1 KB   0.6%    6.6X
    AS:i              29.3 KB   1.1%   97.4 KB   0.3%    3.3X
    NM:i              23.0 KB   0.9%   97.4 KB   0.3%    4.2X
    MAPQ              21.0 KB   0.8%   97.6 KB   0.3%    4.6X
    TXT_HEADER         8.2 KB   0.3%   24.7 KB   0.1%    3.0X
    SA:Z               5.4 KB   0.2%   13.3 KB   0.0%    2.5X
    Other              2.8 KB   0.1%    1.8 MB   5.8%  666.7X
    RNEXT              1.4 KB   0.1%  390.3 KB   1.2%  284.0X
    XQ:i                991 B   0.0%     522 B   0.0%    0.5X
    BGZF                792 B   0.0%         -   0.0%    0.0X
    RNAME               297 B   0.0%  390.3 KB   1.2% 1345.6X
    BAM_BIN              43 B   0.0%  195.1 KB   0.6% 4646.9X
    RG:Z                 42 B   0.0%  195.1 KB   0.6% 4757.6X
    GENOZIP vs BGZF    2.5 MB 100.0%    5.7 MB 100.0%    2.3X
    GENOZIP vs TXT     2.5 MB 100.0%   31.3 MB 100.0%   12.5X
    
In this paritcular example, we observe that the QUAL field consumes 38.1% of the total compressed file size. Therefore, we can expect that ``--optimize-QUAL`` will significantly reduce the compressed file size. In contrast, NM:i, for example, consumes only 0.9% of the compressed file size. Therefore, we can expect that getting rid of NM:i will *not* significantly reduce the compressed file size.

**idxstats**

``genocat --idxstats myfile.bam.genozip``

Calculates idxstats, similar to `samtools idxstats <http://www.htslib.org/doc/samtools-idxstats.html>`_. See :ref:`idxstats`.

**Per-contig coverage and depth**

``genocat --show-coverage myfile.bam.genozip``

An experimental feature for calculating coverage and depth, see :ref:`coverage`.

**Sex assignment**

``genocat --show-sex myfile.bam.genozip``

An experimental feature for determining the sex of a sample, see :ref:`sex`.

**Multi-threading**

By default, Genozip attempts to utilize as many cores as available. For that, it sets the number of threads to be a bit more than the number of cores (a practice known as "over-subscription"), as at any given moment some threads might be idle, waiting for a resource to become available. The ``--threads <number>`` option allows explicit specification of the number of "compute threads" to be used (in addition a small number of I/O threads is used too, usually 1 or 2).

**Memory (RAM) consumption**

In ``genozip``, each compute thread is assigned a segment of the input file, known as a VBlock. By default, the size of the VBlock is set automatically to balance memory consumption and compression ratio for the particular input file, however it may be set explicitly with ``genozip --vblock <megabytes>`` (<megabytes> is an integer between 1 and 2048). A larger VBlock usually results in better compression while a smaller VBlock causes ``genozip`` to consume less RAM. The VBlock size can be observed at the top of the ``--stats`` report. ``genozip``'s memory consumption is linear with (VBlock-size X number-of-threads). 

``genocat`` and ``genounzip`` also consume memory linearly with (VBlock-size X number-of-threads), where VBlock-size is the value used by ``genozip`` of the particular file (it cannot be modified ``genocat`` or ``genounzip``). Usually, ``genocat`` and ``genounzip`` consume significantly less memory compared to ``genozip``.

When using a reference file, usually only the required parts of it are actually loaded to RAM - for example, if the BAM file contains only data of one specific chromosome, then only the reference file data related to that chromosome will be in RAM. If multiple ``genozip``/ ``genocat`` / ``genounzip`` processes are running in parallel, only one copy of the reference file is loaded to memory and shared between all processes, and depending on how busy the computer is, that reference file data might persist in RAM even *between* consecutive runs, saving Genozip the need to load it again from disk. All this all happens behind the scenes.

