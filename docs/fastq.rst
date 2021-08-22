.. _fastq:

Genozip FASTQ files
===================

**Compressing a FASTQ file**

While Genozip is technically capable of compressing FASTQ files without using a reference, in practice, to achieve good compression ratios, compression should always be done against a reference genome.

Preparing the reference file is a one-time step (per-species), where the input file is a FASTA file representing the reference genome.

Notes:

- If the species has multiple versions of its reference genome FASTA, any one of them will work. 
- If there is no reference genome for the target species, a reference genome of a closely related species may be used.
- For meta-genomic applications, it is possible to use a reference FASTA that contains sequences from multiple species, up to a total of 4 Gbp.

::

    $ genozip --make-reference hs37d5.fa.gz  # this generates the file hs37d5.ref.genozip
    
Once a reference file is prepared, we can compress the FASTQ files:

Compressing a single FASTQ file:

::

    $ genozip --reference hs37d5.ref.genozip myfile.fq.gz
    genozip myfile.fq.gz : Done (31 seconds, FASTQ compression ratio: 22.5 - better than .fq.gz by a factor of 4.5)

    $ ls -l myfile*
    -rwxrwxrwx 1 divon divon 1640641 Aug 21 12:28 myfile.fq.genozip
    -rwxrwxrwx 1 divon divon 7338338 Aug 21 12:02 myfile.fq.gz

Note: supported input file extensions include ``.fq`` ``.fq.gz`` ``.fq.bz2`` ``.fq.xz`` and also ``.fastq`` ``.fastq.gz`` ``.fastq.bz2`` ``.fastq.xz``. For FASTQ files with a different extension, use ``--input fastq`` to inform Genozip that this is FASTQ data. 

Tip: Use --REFERENCE instead of --reference to store the reference data as part of the compressed file, obliviating the need for a separate reference file when uncompressing. This is in particular beneficial when binding multiple files together with --output, see :ref:`archiving`. 

**Compressing a paired-end FASTQ files**

For paired-end FASTQ files, it is advisable to compress the two files together, as this improves the compression ratio:

::

    $ genozip --reference hs37d5.ref.genozip --pair myfile-R1.fq.gz myfile-R2.fq.gz
    genozip myfile-R1.fq.gz : Done (2 seconds)
    genozip myfile-R2.fq.gz : Done (8 seconds, FASTQ compression ratio: 20.3 - better than .fq.gz by a factor of 4.3)

    $ ls -l myfile*
    -rwxrwxrwx 1 divon divon 3624227 Aug 21 12:50 myfile-R1+2.fastq.genozip
    -rwxrwxrwx 1 divon divon 7338338 Aug 21 12:02 myfile-R1.fq.gz
    -rwxrwxrwx 1 divon divon 8232187 Aug 21 12:02 myfile-R2.fq.gz
    
| Accessing R1 data: ``genocat myfile-R1+2.fq.genozip --component 1``
| Accessing R2 data: ``genocat myfile-R1+2.fq.genozip --component 2``

Many downstream bioinforamtics tools can accept paired-end FASTQ data in *interleaved* format, for example `bwa mem -p <http://bio-bwa.sourceforge.net/bwa.shtml>`_. To access the data in interleaved format: ``genocat myfile-R1+2.fq.genozip --interleaved``

Some useful command line options (for a full list, see :ref:`genozip manual<genozip>`):

``genozip --test myfile.fq.gz``: after completing the compression, the file is uncompressed in memory, and its MD5 is compared to that of the original file.

``genozip --replace myfile.fq.gz``: the original file is removed after successful compression

**Compressing multiple files into a tar archive**

| ``genozip *.fq.gz --reference hs37d5.ref.genozip --tar mydata.tar``. 
| See details: :ref:`archiving`.

**Optimizing compression**

These are options that modify the file in ways that improve compression. ``--optimize`` is an umbrella option that activates all optimization options.

``genozip myfile.fq.gz --reference hs37d5.ref.genozip --optimize-DESC`` 

Replaces the description line with @filename:read_number. Also - if the 3rd line (the '+' line) contains a copy of the description it is shortened to just '+'.

``genozip myfile.fq.gz --reference hs37d5.ref.genozip --optimize-QUAL`` 

The quality data is optimized as follows:

    ============ ======
    *Old values* *New value*                 
    2-9          6
    10-10        15
    20-24        22
    25-29        27
    \.\.\.
    85-89        87
    90-92        91
    93           Unchanged
    ============ ======

The option ``--stats`` can be used in ``genozip``, ``genounzip`` or ``genocat`` to get a better understanding of the information content of the file. For example:
   
::

    $ genocat --stats myfile-R1+2.fastq.genozip 
    
    FASTQ files (paired): myfile-R1.fq.gz myfile-R2.fq.gz
    Reference: hs37d5.ref.genozip
    Sequences: 200,000   Dictionaries: 25   Vblocks: 6 x 16 MB  Sections: 132
    Genozip version: 12.0.30 github
    Date compressed: 2021-08-21 18:34:02 Cen. Australia Daylight Time
    License v12.0.29 granted to: ***** accepted by:***** on 2021-08-18 20:30:56 Cen. Australia Daylight Time from IP=*****
    
    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    QUAL                    2.0 MB  58.3%   28.3 MB  40.3%   14.1X
    SEQ                     1.3 MB  37.5%   28.3 MB  40.3%   21.9X
    DESC                  144.9 KB   4.1%   12.5 MB  17.8%   88.2X
    Other                   1.1 KB   0.0%    1.1 MB   1.6% 1097.9X
    TXT_HEADER               696 B   0.0%         -   0.0%    0.0X
    LINE3                    246 B   0.0%         -   0.0%    0.0X
    BGZF                     112 B   0.0%         -   0.0%    0.0X
    GENOZIP vs BGZF         3.5 MB 100.0%   14.8 MB 100.0%    4.3X
    GENOZIP vs TXT          3.5 MB 100.0%   70.3 MB 100.0%   20.4X

In this paritcular example, we observe that the quality line consumes 58.3% of the total compressed file size. Therefore, we can expect that ``--optimize-QUAL`` will significantly reduce the compressed file size. In contrast, the description line, in this case, consumes only 4.1% of the compressed file size. Therefore, we can expect that ``--optimize-DESC`` will *not* significantly reduce the compressed file size.

**Uncompressing**

``genounzip myfile.fq.genozip``

Uncompresses a file.

``genocat myfile.fq.genozip``

Uncompresses a file into stdout (i.e. the terminal).

``genounzip --index myfile.fq.genozip``

Uncompresses a file and also generates a FAI index file, using `samtools faidx <http://www.htslib.org/doc/samtools-faidx.html>`_. samtools needs to be installed for this option to work. 

``genounzip --output newname.fq.gz myfile.fq.genozip``

Uncompressing to a particular name. Whether or not the name has a ``.gz`` extension detemines whether the output file is BGZF-compressed.

``genocat --bgzf 6 myfile.fq.genozip`` 
``genounzip --bgzf 6 myfile.fq.genozip`` 

Sets the level BGZF compression (for .fq.gz output format) - from 0 (no compression) to 12 (best yet slowest compression). Absent this option, ``genounzip`` attemps to recover the BGZF compression level of the original file, while ``genocat`` uncompresses without BGZF compression. 
    
**Using in a pipeline**

| Compressing piped input: 
| ``my-pipeline | genozip - --input fastq --output myfile.fq.genozip`` 

| Uncompressing to a pipe: 
| ``genocat myfile.fq.genozip | my-pipeline                     # not paired-end`` 
| ``genocat myfile.R1+2.fq.genozip --interleaved | my-pipeline  # paired-end`` 

**Showing only the descripion lines**

``genocat --header-only myfile.fq.genozip``

**Downsampling**

``genocat --downsample 10,0 myfile.fq.genozip`` 

Displays only the first (#0) read in every 10 reads.

**Grepping**

``genocat --grep ACCTTAAT myfile.fq.genozip`` 

Displays reads with the string "ACCTTAAT" anywhere in the read (description, seqeuence or quality lines) - possibly a substring of a longer string.

``genocat --grep-w ACCTTAAT myfile.fq.genozip`` 

Displays reads with the string "ACCTTAAT" exactly matching a component of the description, or the entire sequence line or the entire quality line.

**Filtering non-ACTGN "bases"**

``genocat --bases ACGTN myfile.fq.genozip``  

Displays only reads in which all characters of the sequence are one of A,C,G,T,N

``genocat --bases ^ACGTN myfile.fq.genozip`` 

Displays only reads in which NOT all characters of the sequence are one of A,C,G,T,N

Note: The list of IUPAC chacacters can be found here: `IUPAC codes <https://www.bioinformatics.org/sms/iupac.html>`_

**Filtering reads by species**

Genozip has the ability to filter FASTQ files by species (taxonomy id). See :ref:`kraken`.

**idxstats**

``genocat --idxstats myfile.fq.genozip``

Calculates approximate idxstats, directly from the FASTQ data, expected to be outputted by `samtools idxstats <http://www.htslib.org/doc/samtools-idxstats.html>`_ after this FASTQ file is mapped against the reference genome. See :ref:`idxstats`.

**Per-contig coverage and depth**

``genocat --show-coverage myfile.fq.genozip``

An experimental feature for calculating coverage and depth directly from a FASTQ file, see :ref:`coverage`.

**Sex assignment**

``genocat --show-sex myfile.fq.genozip``

An experimental feature for determining the sex of a sample from a FASTQ file, see :ref:`sex`.

**Multi-threading**

By default, Genozip attempts to utilize as many cores as available. For that, it sets the number of threads to be a bit more than the number of cores (a practice known as "over-subscription"), as at any given moment some threads might be idle, waiting for a resource to become available. The ``--threads <number>`` option allows explicit specification of the number of "compute threads" to be used (in addition a small number of I/O threads is used too, usually 1 or 2).

**Memory (RAM) consumption**

In ``genozip``, each compute thread is assigned a segment of the input file, known as a VBlock. By default, the size of the VBlock for most FASTQ files is 16MB, however it may be set explicitly with ``genozip --vblock <megabytes>`` (<megabytes> is an integer between 1 and 2048). A larger VBlock usually results in better compression while a smaller VBlock causes ``genozip`` to consume less RAM. The VBlock size can be observed at the top of the ``--stats`` report. ``genozip``'s memory consumption is linear with (VBlock-size X number-of-threads). 

``genocat`` and ``genounzip`` also consume memory linearly with (VBlock-size X number-of-threads), where VBlock-size is the value used by ``genozip`` of the particular file (it cannot be modified ``genocat`` or ``genounzip``). Usually, ``genocat`` and ``genounzip`` consume significantly less memory compared to ``genozip``.

When using a reference file, it is loaded to memory too. If multiple ``genozip``/ ``genocat`` / ``genounzip`` processes are running in parallel, only one copy of the reference file is loaded to memory and shared between all processes, and depending on how busy the computer is, that reference file data might persist in RAM even *between* consecutive runs, saving Genozip the need to load it again from disk. All this all happens behind the scenes.

