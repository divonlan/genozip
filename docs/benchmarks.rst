.. _benchmarks:

Benchmarks
==========

Below are some benchmarks measuring Genozip's performance on a variety of file types. In the case of BAM, CRAM, FASTQ and VCF, Genozip is shown here compressing already-compressed files.

To see more, peer-reviewed, benchmarks, see :ref:`Publications <publications>`.

As can be appreciated from the :ref:`detailed compression reports<benchmark-details>` following the table, the compression gains strongly depend on the specific fields contained in a file, and can therefore vary significantly. 

================================ ================ ========= ========= ========= ============ ============ ================================= ==============================
Details                          Type             Size      .gz Size  .genozip  Genozip gain Gain vs .gz  Data                              Link to data
================================ ================ ========= ========= ========= ============ ============ ================================= ==============================
:ref:`details<benchmark-FASTQ1>` FASTQ (Illumina) 294.7 GB  61.2 GB   13.5 GB   21.8X        4.5X         30x Illumina NovaSeq (R1+R2)      Unpublished
:ref:`details<benchmark-FASTQ2>` FASTQ (Nanopore) 538 MB    268 MB    169.5 MB  3.2X         1.6X         Nanopore (virus) :sup:`1`         `NCBI <https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos/sra-pub-run-1/ERR2708427/ERR2708427.1>`_
:ref:`details<benchmark-BAM1>`   BAM (Illumina)   41 GB               10 GB     4.1X                      30x Illumina NovaSeq + Dragen     Unpublished
:ref:`details<benchmark-BAM2>`   BAM (PacBio CLR) 53.8 GB             28 GB     1.9X                      PacBio CLR (mapped)               `1000 Genozip Project <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20131209 na12878 pacbio/si/NA12878.pacbio.bwa-sw.20140202.bam>`_
:ref:`details<benchmark-BAM3>`   BAM (PacBio CCS) 73.8 GB             38.7 GB   1.9X                      PacBio CCS + Winnowmap (mapped)   `Telomere-to-telomere consortium <https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.1.hifi_20k.wm_2.01.pri.bam>`_
:ref:`details<benchmark-BAM4>`   BAM (Nanopore)   384.1 GB            209.3 GB  1.8X                      Nanopore + Winnowmap (mapped)     `Telomere-to-telomere consortium <https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.1.ont_guppy_3.6.0.wm_2.01.pri.bam>`_
:ref:`details<benchmark-CRAM1>`  CRAM (Illumina)  15 GB               9.8 GB    1.5X                      Illumina NovaSeq + bwa mem        `The European Bioinformatics Institute <ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram>`_
:ref:`details<benchmark-VCF1>`   VCF              175.3 GB  27 GB     7.8 GB    22.4X        3.3X         3202 human samples                `1000 Genome Project <ftp://ftp=trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/20201028_CCDG_14151_B01_GRM_WGS_2020=08=05_chr22.recalibrated_variants.vcf.gz>`_
:ref:`details<benchmark-VCF2>`   VCF              1.2 TB    132.4 GB  31.6 GB   37.9X        4.2X         1135 plant samples                `1001 Genomes - Arabidopsis Thaliana <https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz>`_
:ref:`details<benchmark-VCF3>`   VCF              26 GB     1.9 GB    315.5 MB  84.3X        6.2X         Rice                              `3K Rice Genome <https://3kricegenome.s3.amazonaws.com/9311/IRIS_313-10000.snp.vcf.gz>`_
:ref:`details<benchmark-FASTA1>` FASTA            1.2 GB    254.9 MB  1.5 MB    838.4X       170X         Covid-19 multi-FASTA :sup:`2`     `coronavirus.innar.com <https://coronavirus.innar.com/coronavirus.unwrapped.fasta.zip>`_
:ref:`details<benchmark-GFF1>`   GFF3             3.7 GB    91.8 MB   32.3 MB   117.9X       2.8X         Gene annotation                   `Telomere-to-telomere consortium <https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v1.0.gene_annotation.v4.gff3.gz>`_
:ref:`details<benchmark-ME1>`    23andMe          23.6 MB             4.2 MB    5.7X                      Consumer DNA test "raw data"      Unpublished
:ref:`details<benchmark-LOCS1>`  LOCS (Illumina)  32 MB     12MB      5.9 KB    5430X        2006X        Illumina s.locs file              `BaseSpace Demo Data <https://basespace.illumina.com/run/199406228/NVS1p5-10X-Multi-ATAC-001/files/tree>`_
:ref:`details<benchmark-LOCS2>`  LOCS (Illumina)  6.2 MB    4.4MB     3.3 MB    1.9X         1.3X         Illumina s_X_XXXX.locs file       `BaseSpace Demo Data <https://basespace.illumina.com/run/207258100/MiSeq-v2-0-Sars-CoV-2-samples/files/tree>`_
================================ ================ ========= ========= ========= ============ ============ ================================= ==============================

Notes:

    - The tests were conducted with the ``--best`` option. For BAM, CRAM, Illumina FASTQ files, the ``--reference`` option was used to specify the appropriate reference file. For Illumina FASTQ, the ``--pair`` option was used.
    
    - :sup:`1` the data was converted to fastq with ``fastq-dump`` and compressed to .gz with ``bgzip``. genozip was used with the ``--multiseq`` option.
  
    - :sup:`2` the data was unzipped with ``unzip`` and then compressed to .gz with ``bgzip``. genozip was used with the ``--multiseq`` option.
  
.. _benchmark-details:

**Detailed compression reports**

*The following reports can be produced during compression with* ``genozip --stats`` *or after compression with* ``genocat --stats <myfile>.genozip``.


.. _benchmark-FASTQ1:

**FASTQ - 30x Illumina NovaSeq**

::

    FASTQ files (paired): <redacted>_R1_001.fastq.gz <redacted>_R2_001.fastq.gz
    Reference: hs37d5.ref.genozip
    Sequences: 860,000,926   Dictionaries: 36   Vblocks: 590 x 512 MB  Sections: 10823
    Read name style: Illumina-fastq
    Genozip version: 13.0.6 github
    Date compressed: 2021-12-08 00:49:27 ACDT
    Command line: genozip --best --pair -ft -e hs37d5.ref.genozip --stats <redacted>_R1_001.fastq.gz <redacted>_R2_001.fastq.gz

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %       TXT      %   RATIO
    QUAL                    7.7 GB  56.9%  118.7 GB  40.3%   15.4X
    SEQ                     5.2 GB  38.1%  118.7 GB  40.3%   23.0X
    DESC                  682.4 MB   4.9%   52.4 GB  17.8%   78.7X
    Other                  71.1 KB   0.0%    4.8 GB   1.6% 70888.0X
    LINE3                  23.6 KB   0.0%         -   0.0%    0.0X
    GENOZIP vs BGZF        13.5 GB 100.0%   61.2 GB 100.0%    4.5X
    GENOZIP vs TXT         13.5 GB 100.0%  294.7 GB 100.0%   21.8X


.. _benchmark-FASTQ2:

**FASTQ - Nanopore virus data**

::

    FASTQ file: ERR2708427.1.fastq.gz
    Sequences: 218,903   Dictionaries: 26   Vblocks: 2 x 512 MB  Sections: 37
    Genozip version: 13.0.5 github
    Date compressed: 2021-12-03 10:31:26 ACDT
    Command line: genozip --stats --best --multiseq ERR2708427.1.fastq.gz

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %       TXT      %   RATIO
    QUAL                  123.2 MB  72.7%  254.1 MB  47.3%    2.1X
    SEQ                    41.8 MB  24.6%  254.1 MB  47.3%    6.1X
    DESC                    4.6 MB   2.7%   14.0 MB   2.6%    3.1X
    Other                    470 B   0.0%    1.3 MB   0.2% 2794.5X
    TXT_HEADER               348 B   0.0%         -   0.0%    0.0X
    LINE3                     54 B   0.0%   14.0 MB   2.6% 271727.3X
    GENOZIP vs BGZF       169.5 MB 100.0%  267.7 MB 100.0%    1.6X
    GENOZIP vs TXT        169.5 MB 100.0%  537.5 MB 100.0%    3.2X


.. _benchmark-BAM1:

**BAM - 30x Illumina NovaSeq + Dragen:**

::

    BAM file: <redacted>.bam
    Reference: hs37d5.ref.genozip
    Alignments: 636,660,181   Dictionaries: 161   Vblocks: 389 x 512 MB  Sections: 17340
    Sorting: Sorted by POS
    Read name style: Illumina
    Genozip version: 13.0.2 github
    Date compressed: 2021-11-09 17:14:19 Cen. Australia Daylight Time

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    QUAL                    5.7 GB  57.2%   87.9 GB  45.3%   15.3X
    QNAME                   1.9 GB  19.3%   23.4 GB  12.1%   12.1X
    SEQ                   890.2 MB   8.7%   46.5 GB  24.0%   53.5X
    PNEXT                 441.8 MB   4.3%    2.4 GB   1.2%    5.5X
    POS                   211.0 MB   2.1%    2.4 GB   1.2%   11.5X
    CIGAR                 183.2 MB   1.8%    3.8 GB   1.9%   21.1X
    FLAG                  151.9 MB   1.5%    1.2 GB   0.6%    8.0X
    XS:i                  139.6 MB   1.4%  973.5 MB   0.5%    7.0X
    AS:i                  124.6 MB   1.2%    2.4 GB   1.2%   19.5X
    TLEN                  119.6 MB   1.2%    2.4 GB   1.2%   20.3X
    XQ:i                   70.4 MB   0.7%  580.0 MB   0.3%    8.2X
    Other                  41.7 MB   0.4%   10.1 GB   5.2%  248.6X
    MAPQ                   20.3 MB   0.2%  607.2 MB   0.3%   29.9X
    SA:Z                    7.8 MB   0.1%   31.2 MB   0.0%    4.0X
    RNEXT                   6.6 MB   0.1%    2.4 GB   1.2%  367.5X
    RNAME                  24.4 KB   0.0%    2.4 GB   1.2% 101926.8X
    NM:i                   12.6 KB   0.0%    2.4 GB   1.2% 196672.2X
    TXT_HEADER              5.6 KB   0.0%   16.8 KB   0.0%    3.0X
    BAM_BIN                   43 B   0.0%    1.2 GB   0.6% 29612100.0X
    RG:Z                      42 B   0.0%    1.2 GB   0.6% 30317150.0X
    GENOZIP vs BGZF        10.0 GB 100.0%   41.0 GB 100.0%    4.1X
    GENOZIP vs TXT         10.0 GB 100.0%  194.0 GB 100.0%   19.4X


.. _benchmark-BAM2:

**BAM - PacBio CLR (mapped)**

::

    BAM file: NA12878.pacbio.bwa-sw.20140202.bam
    Reference: hs37d5.ref.genozip
    Alignments: 25,968,256   Dictionaries: 163   Vblocks: 215 x 512 MB  Sections: 9825
    Sorting: Sorted by POS
    Read name style: PacBio-Range
    Genozip version: 13.0.5 github
    Date compressed: 2021-11-30 21:05:14 ACDT
    Command line: genozip --reference hs37d5.ref.genozip --best --stats NA12878.pacbio.bwa-sw.20140202.bam

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    QUAL                   19.3 GB  69.2%   41.1 GB  38.3%    2.1X
    SEQ                     3.7 GB  13.3%   20.6 GB  19.2%    5.6X
    CIGAR                   2.4 GB   8.6%   22.6 GB  21.1%    9.4X
    SA:Z                    2.3 GB   8.2%   19.1 GB  17.8%    8.4X
    QNAME                 103.6 MB   0.4%    1.9 GB   1.8%   18.7X
    AS:i                   32.0 MB   0.1%   38.0 MB   0.0%    1.2X
    XS:i                   31.1 MB   0.1%   29.2 MB   0.0%    0.9X
    POS                    24.5 MB   0.1%   99.1 MB   0.1%    4.0X
    Other                  11.6 MB   0.0%  520.0 MB   0.5%   44.9X
    MAPQ                   10.6 MB   0.0%   24.8 MB   0.0%    2.3X
    FLAG                    7.1 MB   0.0%   49.5 MB   0.0%    6.9X
    PNEXT                  51.7 KB   0.0%   99.1 MB   0.1% 1962.5X
    NM:i                   49.2 KB   0.0%   29.9 MB   0.0%  622.7X
    RNAME                  10.8 KB   0.0%   99.1 MB   0.1% 9369.7X
    RNEXT                   9.5 KB   0.0%   99.1 MB   0.1% 10686.5X
    TXT_HEADER              4.0 KB   0.0%   19.7 KB   0.0%    5.0X
    TLEN                      83 B   0.0%   99.1 MB   0.1% 1251482.2X
    RG:Z                      56 B   0.0%  396.2 MB   0.4% 7419501.5X
    PG:Z                      55 B   0.0%  371.5 MB   0.3% 7082251.5X
    BAM_BIN                   43 B   0.0%   49.5 MB   0.0% 1207825.9X
    GENOZIP vs BGZF        28.0 GB 100.0%   53.8 GB 100.0%    1.9X
    GENOZIP vs TXT         28.0 GB 100.0%  107.3 GB 100.0%    3.8X


.. _benchmark-BAM3:

**BAM - PacBio CCS (mapped)**

::

    BAM file: chm13.draft_v1.1.hifi_20k.wm_2.01.pri.bam
    Reference: chm13.draft_v1.1.ref.genozip
    Alignments: 5,575,318   Dictionaries: 164   Vblocks: 286 x 512 MB  Sections: 11865
    Sorting: Sorted by POS
    Read name style: PacBio-Label
    Genozip version: 13.0.5 conda
    Date compressed: 2021-12-08 13:35:00 ACDT
    Command line: genozip -ftw -e chm13.draft_v1.1.ref.genozip --best chm13.draft_v1.1.hifi_20k.wm_2.01.pri.bam
    
    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %       TXT      %   RATIO
    QUAL                   38.3 GB  98.9%   93.5 GB  65.5%    2.4X
    CIGAR                 273.3 MB   0.7%    1.4 GB   1.0%    5.4X
    SEQ                    75.4 MB   0.2%   46.8 GB  32.7%  635.7X
    QNAME                  16.7 MB   0.0%  188.1 MB   0.1%   11.2X
    AS:i                    9.4 MB   0.0%   10.6 MB   0.0%    1.1X
    ms:i                    9.0 MB   0.0%   10.6 MB   0.0%    1.2X
    s1:i                    8.2 MB   0.0%   10.6 MB   0.0%    1.3X
    POS                     7.4 MB   0.0%   21.3 MB   0.0%    2.9X
    Other                   5.9 MB   0.0%  196.8 MB   0.1%   33.1X
    cm:i                    5.9 MB   0.0%   10.6 MB   0.0%    1.8X
    rl:i                    5.4 MB   0.0%    6.9 MB   0.0%    1.3X
    de:f                    4.0 MB   0.0%   21.3 MB   0.0%    5.3X
    s2:i                  842.4 KB   0.0%    5.7 MB   0.0%    6.9X
    FLAG                  721.6 KB   0.0%   10.6 MB   0.0%   15.1X
    SA:Z                  330.1 KB   0.0%  596.4 KB   0.0%    1.8X
    MD:Z                   35.2 KB   0.0%  533.8 MB   0.4% 15524.7X
    MAPQ                   33.7 KB   0.0%    5.3 MB   0.0%  161.6X
    PNEXT                  32.5 KB   0.0%   21.3 MB   0.0%  670.5X
    nn:i                   13.2 KB   0.0%    5.3 MB   0.0%  412.3X
    RNEXT                  12.4 KB   0.0%   21.3 MB   0.0% 1753.8X
    RNAME                  12.1 KB   0.0%   21.3 MB   0.0% 1801.1X
    TXT_HEADER              2.3 KB   0.0%   24.6 KB   0.0%   10.8X
    zd:i                    1.2 KB   0.0%      69 B   0.0%    0.1X
    tp:A                     233 B   0.0%    5.3 MB   0.0% 23928.4X
    NM:i                      85 B   0.0%    5.3 MB   0.0% 65785.6X
    TLEN                      83 B   0.0%   21.3 MB   0.0% 268690.0X
    BAM_BIN                   43 B   0.0%   10.6 MB   0.0% 259317.1X
    GENOZIP vs BGZF        38.7 GB 100.0%   73.8 GB 100.0%    1.9X
    GENOZIP vs TXT         38.7 GB 100.0%  142.8 GB 100.0%    3.7X


.. _benchmark-BAM4:

**BAM - Nanopore (mapped)**

::

    BAM file: chm13.draft_v1.1.ont_guppy_3.6.0.wm_2.01.pri.bam
    Reference: chm13.draft_v1.1.ref.genozip
    Alignments: 13,364,876   Dictionaries: 168   Vblocks: 1321 x 512 MB  Sections: 65852
    Sorting: Sorted by POS
    Read name style: Nanopore
    Genozip version: 13.0.6 github
    Date compressed: 2021-12-08 02:36:21 ACDT
    Command line: genozip --best -wtf -e chm13.draft_v1.1.ref.genozip chm13.draft_v1.1.ont_guppy_3.6.0.wm_2.01.pri.bam

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %       TXT      %   RATIO
    QUAL                  186.1 GB  88.9%  351.6 GB  53.3%    1.9X
    CIGAR                  11.8 GB   5.6%   84.2 GB  12.8%    7.2X
    SEQ                    11.0 GB   5.3%  175.9 GB  26.6%   15.9X
    QNAME                 214.8 MB   0.1%  484.3 MB   0.1%    2.3X
    MD:Z                   39.6 MB   0.0%   46.7 GB   7.1% 1208.4X
    AS:i                   27.8 MB   0.0%   30.6 MB   0.0%    1.1X
    SA:Z                   26.9 MB   0.0%  154.8 MB   0.0%    5.8X
    ms:i                   25.0 MB   0.0%   30.6 MB   0.0%    1.2X
    s1:i                   21.9 MB   0.0%   23.1 MB   0.0%    1.1X
    Other                  20.9 MB   0.0%  475.1 MB   0.1%   22.8X
    de:f                   17.4 MB   0.0%   51.0 MB   0.0%    2.9X
    cm:i                   16.3 MB   0.0%   18.6 MB   0.0%    1.1X
    POS                    15.3 MB   0.0%   51.0 MB   0.0%    3.3X
    CG:B:I                 14.3 MB   0.0%  118.9 MB   0.0%    8.3X
    rl:i                    6.7 MB   0.0%   13.6 MB   0.0%    2.0X
    s2:i                    3.7 MB   0.0%   13.4 MB   0.0%    3.6X
    FLAG                    2.1 MB   0.0%   25.5 MB   0.0%   12.3X
    MAPQ                    1.6 MB   0.0%   12.7 MB   0.0%    8.0X
    PNEXT                 167.9 KB   0.0%   51.0 MB   0.0%  311.0X
    zd:i                  126.6 KB   0.0%  109.7 KB   0.0%    0.9X
    nn:i                   59.4 KB   0.0%   12.7 MB   0.0%  219.6X
    RNEXT                  56.9 KB   0.0%   51.0 MB   0.0%  917.7X
    RNAME                  52.6 KB   0.0%   51.0 MB   0.0%  993.3X
    NM:i                   39.2 KB   0.0%   21.8 MB   0.0%  569.9X
    tp:A                   28.1 KB   0.0%   12.7 MB   0.0%  464.4X
    TXT_HEADER              1.8 KB   0.0%   12.0 KB   0.0%    6.8X
    TLEN                      83 B   0.0%   51.0 MB   0.0% 644090.4X
    BAM_BIN                   43 B   0.0%   25.5 MB   0.0% 621622.1X
    GENOZIP vs BGZF       209.3 GB 100.0%  384.1 GB 100.0%    1.8X
    GENOZIP vs TXT        209.3 GB 100.0%  660.1 GB 100.0%    3.2X


.. _benchmark-CRAM1:

**CRAM - Illumina NovaSeq + bwa mem**

::

    SAM file: NA12878.final.cram
    Reference: GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    Alignments: 768,580,569   Dictionaries: 160   Vblocks: 749 x 512 MB  Sections: 44749
    Sorting: Sorted by POS
    Read name style: Illumina
    Genozip version: 13.0.2 conda
    Date compressed: 2021-11-09 17:32:34 ACDT

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    QUAL                    2.7 GB  27.7%  108.1 GB  28.9%   39.9X
    QNAME                   2.4 GB  24.7%   27.5 GB   7.4%   11.4X
    XA:Z                    1.4 GB  14.8%   13.7 GB   3.7%    9.5X
    SEQ                   925.9 MB   9.3%  108.1 GB  28.9%  119.5X
    PNEXT                 630.8 MB   6.3%    6.6 GB   1.8%   10.7X
    XS:i                  305.0 MB   3.1%    2.1 GB   0.6%    7.2X
    POS                   275.7 MB   2.8%    6.6 GB   1.8%   24.4X
    RG:Z                  275.0 MB   2.8%   29.3 GB   7.8%  109.3X
    FLAG                  213.7 MB   2.1%    2.6 GB   0.7%   12.6X
    AS:i                  129.6 MB   1.3%    2.8 GB   0.8%   22.4X
    SA:Z                  115.3 MB   1.2%  838.7 MB   0.2%    7.3X
    MAPQ                  112.9 MB   1.1%    2.1 GB   0.6%   19.1X
    MC:Z                  106.3 MB   1.1%    3.8 GB   1.0%   36.5X
    CIGAR                  91.8 MB   0.9%    3.8 GB   1.0%   42.2X
    MQ:i                   37.0 MB   0.4%    2.1 GB   0.6%   57.5X
    RNEXT                  20.6 MB   0.2%    1.5 GB   0.4%   76.5X
    Other                  17.6 MB   0.2%   29.0 GB   7.8% 1692.4X
    TLEN                   14.0 MB   0.1%    3.2 GB   0.9%  233.7X
    pa:f                    6.9 MB   0.1%   38.5 MB   0.0%    5.6X
    MD:Z                    1.2 MB   0.0%    4.0 GB   1.1% 3464.2X
    RNAME                 374.0 KB   0.0%    4.1 GB   1.1% 11515.0X
    TXT_HEADER             72.5 KB   0.0%  626.9 KB   0.0%    8.6X
    NM:i                   66.8 KB   0.0%    1.4 GB   0.4% 22624.0X
    Reference                112 B   0.0%         -   0.0%    0.0X
    PG:Z                      55 B   0.0%   10.7 GB   2.9% 209612880.0X
    BAM_BIN                   43 B   0.0%         -   0.0%    0.0X
    TOTAL                   9.8 GB 100.0%  374.2 GB 100.0%   38.3X


.. _benchmark-VCF1:

**VCF - 3202 samples from the 1000 Genome Project**

::

    VCF file: 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz
    Samples: 3202   Variants: 1,927,372   Dictionaries: 401   Vblocks: 351 x 512 MB  Sections: 158051
    Genozip version: 13.0.3 github
    Date compressed: 2021-11-14 08:56:28 ACDT

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    FORMAT/PL               5.0 GB  64.3%   59.4 GB  33.9%   11.8X
    FORMAT/AD               2.4 GB  30.9%   24.7 GB  14.1%   10.2X
    FORMAT/GT              67.2 MB   0.8%   17.2 GB   9.8%  262.8X
    FORMAT/GQ              65.3 MB   0.8%   11.2 GB   6.4%  176.0X
    FORMAT/PID             63.5 MB   0.8%    3.2 GB   1.8%   51.6X
    FORMAT/PGT             38.4 MB   0.5%    2.3 GB   1.3%   60.1X
    FORMAT/DP              19.6 MB   0.2%   11.4 GB   6.5%  593.2X
    FORMAT/AB              14.9 MB   0.2%    5.7 GB   3.2%  388.8X
    INFO/AC_Het_EUR_unre    5.9 MB   0.1%   26.9 MB   0.0%    4.6X
    QUAL                    5.1 MB   0.1%   13.7 MB   0.0%    2.7X
    INFO/DP                 4.4 MB   0.1%   10.5 MB   0.0%    2.4X
    INFO/AF_AMR_unrel       3.4 MB   0.0%   20.4 MB   0.0%    6.0X
    INFO/VQSLOD             3.3 MB   0.0%    9.6 MB   0.0%    2.9X
    INFO/FS                 3.0 MB   0.0%    7.3 MB   0.0%    2.5X
    INFO/AF                 2.9 MB   0.0%   21.9 MB   0.0%    7.6X
    INFO/MQRankSum          2.8 MB   0.0%    9.3 MB   0.0%    3.3X
    INFO/SOR                2.8 MB   0.0%    9.0 MB   0.0%    3.2X
    INFO/BaseQRankSum       2.8 MB   0.0%    9.2 MB   0.0%    3.4X
    INFO/QD                 2.7 MB   0.0%    8.5 MB   0.0%    3.1X
    INFO/ReadPosRankSum     2.7 MB   0.0%    9.0 MB   0.0%    3.3X
    INFO/AF_EUR_unrel       2.7 MB   0.0%   16.9 MB   0.0%    6.2X
    INFO/ClippingRankSum    2.7 MB   0.0%    9.3 MB   0.0%    3.4X
    INFO/AC_AMR_unrel       2.3 MB   0.0%    5.8 MB   0.0%    2.5X
    INFO/MLEAF              2.3 MB   0.0%   17.3 MB   0.0%    7.7X
    INFO/AF_AFR             2.2 MB   0.0%   11.9 MB   0.0%    5.4X
    INFO/ExcHet             2.2 MB   0.0%   10.8 MB   0.0%    4.9X
    INFO/AC_Het_AFR         2.1 MB   0.0%    5.7 MB   0.0%    2.8X
    INFO/InbreedingCoeff    2.0 MB   0.0%   10.8 MB   0.0%    5.5X
    INFO/ExcHet_AFR         1.9 MB   0.0%    7.8 MB   0.0%    4.0X
    REF+ALT                 1.9 MB   0.0%   12.4 MB   0.0%    6.4X
    INFO/MLEAC              1.9 MB   0.0%    3.5 MB   0.0%    1.8X
    INFO/HWE                1.9 MB   0.0%    6.2 MB   0.0%    3.2X
    INFO/AC_EUR_unrel       1.9 MB   0.0%    5.6 MB   0.0%    2.9X
    INFO/AC_Het             1.9 MB   0.0%    3.4 MB   0.0%    1.8X
    INFO/AF_SAS             1.6 MB   0.0%    9.2 MB   0.0%    5.6X
    INFO/AF_AMR             1.6 MB   0.0%    8.9 MB   0.0%    5.5X
    INFO/AC_AFR             1.6 MB   0.0%    3.1 MB   0.0%    1.9X
    INFO/AF_EUR             1.6 MB   0.0%    8.6 MB   0.0%    5.4X
    INFO/AF_SAS_unrel       1.6 MB   0.0%    8.8 MB   0.0%    5.6X
    INFO/AC_Het_SAS         1.5 MB   0.0%    5.4 MB   0.0%    3.5X
    INFO/AF_EAS             1.5 MB   0.0%    8.5 MB   0.0%    5.6X
    INFO/AC_Het_AMR         1.5 MB   0.0%    5.3 MB   0.0%    3.5X
    POS                     1.5 MB   0.0%   16.5 MB   0.0%   11.2X
    INFO/AC_Het_EUR         1.5 MB   0.0%    5.4 MB   0.0%    3.7X
    INFO/HWE_AFR            1.4 MB   0.0%    5.0 MB   0.0%    3.5X
    INFO/AC_Het_EAS         1.4 MB   0.0%    5.3 MB   0.0%    3.8X
    INFO/ExcHet_AMR         1.4 MB   0.0%    6.2 MB   0.0%    4.4X
    INFO/ExcHet_SAS         1.4 MB   0.0%    5.9 MB   0.0%    4.3X
    INFO/MQ                 1.4 MB   0.0%    5.8 MB   0.0%    4.3X
    INFO/ExcHet_EUR         1.3 MB   0.0%    5.7 MB   0.0%    4.2X
    INFO/ExcHet_EAS         1.3 MB   0.0%    5.4 MB   0.0%    4.3X
    INFO/AC_SAS             1.2 MB   0.0%    2.8 MB   0.0%    2.3X
    INFO/AC_AMR             1.2 MB   0.0%    2.8 MB   0.0%    2.3X
    INFO/AC_EUR             1.2 MB   0.0%    2.8 MB   0.0%    2.4X
    INFO/AC_SAS_unrel       1.2 MB   0.0%    2.8 MB   0.0%    2.4X
    INFO/AC_EAS             1.1 MB   0.0%    2.8 MB   0.0%    2.5X
    INFO/HWE_SAS            1.1 MB   0.0%    4.2 MB   0.0%    4.0X
    INFO/HWE_AMR            1.0 MB   0.0%    4.1 MB   0.0%    4.0X
    INFO/HWE_EUR            1.0 MB   0.0%    4.1 MB   0.0%    4.0X
    INFO/HWE_EAS          981.2 KB   0.0%    4.0 MB   0.0%    4.1X
    INFO/AC_Hom           936.7 KB   0.0%    2.8 MB   0.0%    3.1X
    INFO/ME               926.4 KB   0.0%    4.5 MB   0.0%    4.9X
    INFO/AC               778.4 KB   0.0%    3.5 MB   0.0%    4.6X
    INFO/AN_AMR_unrel     537.4 KB   0.0%   12.8 MB   0.0%   24.5X
    INFO                  479.4 KB   0.0%    1.6 GB   0.9% 3600.3X
    INFO/AN_EUR_unrel     472.6 KB   0.0%   14.4 MB   0.0%   31.2X
    INFO/AN               457.6 KB   0.0%    7.4 MB   0.0%   16.4X
    Other                 391.4 KB   0.0%   38.1 GB  21.7% 102014.1X
    FORMAT                347.0 KB   0.0%   37.9 MB   0.0%  111.8X
    INFO/culprit          340.1 KB   0.0%    5.5 MB   0.0%   16.7X
    INFO/AN_AFR           337.9 KB   0.0%    7.3 MB   0.0%   22.2X
    INFO/AN_EUR           285.7 KB   0.0%    7.3 MB   0.0%   26.2X
    INFO/AN_SAS           283.2 KB   0.0%    7.3 MB   0.0%   26.4X
    INFO/AN_EAS           273.6 KB   0.0%    7.3 MB   0.0%   27.3X
    INFO/AN_AMR           272.8 KB   0.0%    5.5 MB   0.0%   20.7X
    INFO/AN_SAS_unrel     266.6 KB   0.0%    5.5 MB   0.0%   21.2X
    FILTER                147.0 KB   0.0%   15.0 MB   0.0%  104.6X
    INFO/NEGATIVE_TRAIN_   14.1 KB   0.0%         -   0.0%    0.0X
    TXT_HEADER             13.2 KB   0.0%  201.7 KB   0.0%   15.2X
    INFO/POSITIVE_TRAIN_   12.2 KB   0.0%         -   0.0%    0.0X
    COORDS                   536 B   0.0%         -   0.0%    0.0X
    CHROM                    139 B   0.0%   11.0 MB   0.0% 83195.9X
    ID                        42 B   0.0%    3.7 MB   0.0% 91779.6X
    INFO/MQ0                  42 B   0.0%    1.8 MB   0.0% 45889.8X
    GENOZIP vs BGZF         7.8 GB 100.0%   26.0 GB 100.0%    3.3X
    GENOZIP vs TXT          7.8 GB 100.0%  175.3 GB 100.0%   22.4X


.. _benchmark-VCF2:

**VCF - 1135 Arabidopsis Thaliana samples**

::

    VCF file: 1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz
    Samples: 1135   Variants: 119,146,348   Dictionaries: 190   Vblocks: 2398 x 512 MB  Sections: 56957
    Genozip version: 13.0.5 github
    Date compressed: 2021-11-22 11:28:13 ACDT
    Command line: genozip --stats --best 1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    FORMAT/DP              22.6 GB  71.6%  206.5 GB  17.2%    9.1X
    FORMAT/GQ               8.7 GB  27.5%  231.7 GB  19.3%   26.6X
    FORMAT/GT             237.1 MB   0.7%  377.8 GB  31.5% 1631.6X
    REF+ALT                41.4 MB   0.1%  460.5 MB   0.0%   11.1X
    QUAL                  928.5 KB   0.0%  340.9 MB   0.0%  376.0X
    Other                 296.1 KB   0.0%  377.8 GB  31.5% 1337972.4X
    POS                   236.6 KB   0.0%  969.7 MB   0.1% 4197.5X
    CHROM                  74.6 KB   0.0%  227.3 MB   0.0% 3119.3X
    TXT_HEADER              1.5 KB   0.0%    6.0 KB   0.0%    3.9X
    COORDS                   547 B   0.0%         -   0.0%    0.0X
    INFO                      76 B   0.0%  454.5 MB   0.0% 6270860.5X
    FORMAT                    51 B   0.0% 1022.6 MB   0.1% 21025826.0X
    INFO/DP                   49 B   0.0%  545.3 MB   0.0% 11669053.0X
    FILTER                    45 B   0.0%  568.1 MB   0.0% 13238482.0X
    ID                        42 B   0.0%  227.3 MB   0.0% 5673636.0X
    GENOZIP vs BGZF        31.6 GB 100.0%  132.4 GB 100.0%    4.2X
    GENOZIP vs TXT         31.6 GB 100.0%    1.2 TB 100.0%   37.9X


.. _benchmark-VCF3:

**VCF - 3K Rice Genome**

::

    VCF file: IRIS_313-10000.snp.vcf.gz
    Samples: 1   Variants: 409,606,670   Dictionaries: 302   Vblocks: 52 x 512 MB  Sections: 3288
    Genozip version: 13.0.4 github
    Date compressed: 2021-11-20 20:51:00 ACDT

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    REF+ALT                93.3 MB  29.6%    1.5 GB   5.9%   16.7X
    QUAL                   70.0 MB  22.2%    2.0 GB   7.8%   29.7X
    INFO/MQ                68.2 MB  21.6%    1.7 GB   6.4%   25.0X
    INFO/DP                52.3 MB  16.6%  593.1 MB   2.2%   11.3X
    INFO/MQ0                6.6 MB   2.1%  360.3 MB   1.4%   54.5X
    INFO                    3.4 MB   1.1%    5.8 GB  22.4% 1757.9X
    FORMAT/DP               2.7 MB   0.8%  561.5 MB   2.1%  211.5X
    FORMAT/PL               2.6 MB   0.8%   11.2 MB   0.0%    4.3X
    FORMAT/GT               2.6 MB   0.8%    1.1 GB   4.4%  454.9X
    Other                   2.6 MB   0.8%  716.6 MB   2.7%  280.6X
    FORMAT                  2.5 MB   0.8%    2.1 GB   8.1%  848.9X
    INFO/QD                 2.1 MB   0.7%    6.8 MB   0.0%    3.3X
    FORMAT/AD               1.1 MB   0.3%    5.2 MB   0.0%    4.7X
    FILTER                902.5 KB   0.3%  785.2 MB   3.0%  890.9X
    FORMAT/GQ             811.5 KB   0.3%    2.7 MB   0.0%    3.4X
    POS                   715.0 KB   0.2%    3.2 GB  12.2% 4627.1X
    INFO/MQRankSum        647.9 KB   0.2%    2.3 MB   0.0%    3.6X
    INFO/ReadPosRankSum   643.2 KB   0.2%    2.2 MB   0.0%    3.5X
    INFO/BaseQRankSum     638.8 KB   0.2%    2.2 MB   0.0%    3.5X
    INFO/HaplotypeScore   443.5 KB   0.1%    7.7 MB   0.0%   17.8X
    INFO/FS               373.2 KB   0.1%    7.2 MB   0.0%   19.7X
    CHROM                 137.7 KB   0.0%    4.4 GB  17.1% 33717.5X
    INFO/MLEAF             82.6 KB   0.0%    6.1 MB   0.0%   75.2X
    INFO/MLEAC             82.6 KB   0.0%    1.4 MB   0.0%   17.7X
    INFO/AF                82.1 KB   0.0%    6.1 MB   0.0%   75.6X
    INFO/RPA               77.6 KB   0.0%  317.9 KB   0.0%    4.1X
    INFO/RU                64.3 KB   0.0%  138.0 KB   0.0%    2.1X
    TXT_HEADER             46.7 KB   0.0%  546.0 KB   0.0%   11.7X
    INFO/Dels              35.0 KB   0.0%    5.1 MB   0.0%  149.8X
    INFO/AC                 4.7 KB   0.0%    1.4 MB   0.0%  313.4X
    INFO/STR                2.1 KB   0.0%         -   0.0%    0.0X
    COORDS                   547 B   0.0%         -   0.0%    0.0X
    INFO/AN                   83 B   0.0%  321.7 MB   1.2% 4064608.5X
    ID                        42 B   0.0%  781.3 MB   2.9% 19505078.0X
    GENOZIP vs BGZF       315.5 MB 100.0%    1.9 GB 100.0%    6.2X
    GENOZIP vs TXT        315.5 MB 100.0%   26.0 GB 100.0%   84.3X


.. _benchmark-FASTA1:

**FASTA - Covid-19 multi-FASTA**

::

    FASTA file: coronavirus.unwrapped.fasta.gz
    Lines: 89,836   Dictionaries: 11   Vblocks: 3 x 512 MB  Sections: 32
    Sequence type: Nucleotide bases
    Genozip version: 13.0.5 github
    Date compressed: 2021-12-03 11:25:06 ACDT
    Command line: genozip --best --stats --multiseq coronavirus.unwrapped.fasta.gz

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %       TXT      %   RATIO
    NONREF                  1.4 MB  92.6%    1.2 GB 100.0%  904.8X
    Other                  95.3 KB   6.2%   87.7 KB   0.0%    0.9X
    DESC                   18.1 KB   1.2%  526.3 KB   0.0%   29.1X
    GENOZIP vs BGZF         1.5 MB 100.0%  254.9 MB 100.0%  170.0X
    GENOZIP vs TXT          1.5 MB 100.0%    1.2 GB 100.0%  838.4X


.. _benchmark-GFF1:

**GFF3 - Gene annotation**

::

    GFF3 file: chm13.draft_v1.0.gene_annotation.v4.gff3.gz
    Sequences: 3,873,663   Dictionaries: 79   Vblocks: 8 x 512 MB  Sections: 562
    Genozip version: 13.0.5 conda
    Date compressed: 2021-12-08 12:07:18 ACDT
    Command line: genozip -ft --best --stats chm13.draft_v1.0.gene_annotation.v4.gff3.gz

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %       TXT      %   RATIO
    ID                      4.5 MB  13.9%   78.6 MB   2.1%   17.5X
    END                     3.0 MB   9.1%   34.0 MB   0.9%   11.5X
    alignment_id            2.0 MB   6.1%   69.1 MB   1.8%   35.0X
    source_transcript       1.9 MB   5.9%   60.8 MB   1.6%   31.9X
    source_gene_common_n    1.9 MB   5.8%   56.2 MB   1.5%   29.9X
    START                   1.7 MB   5.2%   34.0 MB   0.9%   20.3X
    transcript_name         1.7 MB   5.1%   35.9 MB   0.9%   21.7X
    Parent                  1.6 MB   4.9%   50.9 MB   1.3%   32.1X
    havana_transcript       1.5 MB   4.7%   64.2 MB   1.7%   42.5X
    transcript_id           1.3 MB   4.1%   50.9 MB   1.3%   38.6X
    intron_annotation_su    1.2 MB   3.8%  101.0 MB   2.7%   83.3X
    exon_anotation_suppo    1.0 MB   3.2%  108.3 MB   2.8%  104.2X
    protein_id            997.4 KB   3.0%   45.0 MB   1.2%   46.2X
    adj_stop              825.9 KB   2.5%   22.3 MB   0.6%   27.6X
    ATTRS                 818.0 KB   2.5%    2.1 GB  56.8% 2711.0X
    adj_start             804.7 KB   2.4%   22.3 MB   0.6%   28.3X
    tag                   689.3 KB   2.1%   84.8 MB   2.2%  125.9X
    transcript_support_l  569.0 KB   1.7%    3.1 MB   0.1%    5.7X
    transcript_biotype    420.6 KB   1.3%   53.9 MB   1.4%  131.2X
    Name                  394.6 KB   1.2%   21.8 MB   0.6%   56.6X
    TYPE                  386.3 KB   1.2%   22.4 MB   0.6%   59.4X
    gene_name             326.5 KB   1.0%   51.3 MB   1.3%  160.9X
    PHASE                 321.0 KB   1.0%    7.4 MB   0.2%   23.6X
    SCORE                 315.6 KB   1.0%    8.1 MB   0.2%   26.2X
    source_gene           293.2 KB   0.9%   65.8 MB   1.7%  229.7X
    ccdsid                258.3 KB   0.8%   14.6 MB   0.4%   57.9X
    havana_gene           241.3 KB   0.7%   71.0 MB   1.9%  301.3X
    proper_orf            205.7 KB   0.6%   14.9 MB   0.4%   74.4X
    valid_stop            187.7 KB   0.6%   14.9 MB   0.4%   81.3X
    level                 176.8 KB   0.5%    3.5 MB   0.1%   20.4X
    valid_start           170.5 KB   0.5%   14.9 MB   0.4%   89.2X
    hgnc_id               164.9 KB   0.5%   32.1 MB   0.8%  199.6X
    gene_id               129.5 KB   0.4%   51.7 MB   1.4%  408.8X
    frameshift            127.2 KB   0.4%   11.1 MB   0.3%   89.4X
    unfiltered_paralogy   115.3 KB   0.3%   11.7 MB   0.3%  103.9X
    paralogy               76.8 KB   0.2%   11.1 MB   0.3%  147.7X
    transcript_class       76.4 KB   0.2%   30.6 MB   0.8%  409.8X
    gene_biotype           54.5 KB   0.2%   49.9 MB   1.3%  937.7X
    SEQID                  45.4 KB   0.1%   20.2 MB   0.5%  456.6X
    STRAND                 33.0 KB   0.1%    7.4 MB   0.2%  229.2X
    transcript_modes       32.5 KB   0.1%   29.3 MB   0.8%  921.5X
    rna_support            31.0 KB   0.1%   10.3 MB   0.3%  340.7X
    collapsed_gene_names   22.2 KB   0.1%   15.1 MB   0.4%  695.0X
    collapsed_gene_ids     20.8 KB   0.1%   18.7 MB   0.5%  917.8X
    gene_alternate_conti   18.2 KB   0.1%    1.6 MB   0.0%   90.2X
    alternative_source_t    5.7 KB   0.0%   11.1 MB   0.3% 1982.3X
    ont                     4.6 KB   0.0%  679.7 KB   0.0%  147.0X
    Other                   1.7 KB   0.0%         -   0.0%    0.0X
    possible_split_gene_    1.3 KB   0.0%   12.5 KB   0.0%    9.5X
    novel_5p_cap            1.0 KB   0.0%  394.9 KB   0.0%  383.7X
    TXT_HEADER               364 B   0.0%      16 B   0.0%    0.0X
    COMMENT                  328 B   0.0%         -   0.0%    0.0X
    SOURCE                   160 B   0.0%   14.8 MB   0.4% 97251.2X
    extra_paralog            159 B   0.0%   18.5 MB   0.5% 121710.3X
    novel_poly_a             130 B   0.0%  390.6 KB   0.0% 3076.8X
    reference_support         45 B   0.0%   13.6 MB   0.4% 317191.9X
    GENOZIP vs GZ          32.3 MB 100.0%   91.8 MB 100.0%    2.8X
    GENOZIP vs TXT         32.3 MB 100.0%    3.7 GB 100.0%  117.9X


.. _benchmark-ME1:

**23andMe - Consumer DNA test "raw data"**

::

    23ANDME file: genome_<redacted>.txt
    SNPs: 960,613   Dictionaries: 7   Vblocks: 2 x 16 MB  Sections: 27
    Genozip version: 13.0.2 github
    Date compressed: 2021-11-09 18:05:24 Cen. Australia Daylight Time

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    ID                      2.3 MB  55.2%    9.3 MB  39.5%    4.1X
    POS                     1.5 MB  37.1%    8.4 MB  35.8%    5.5X
    GENOTYPE              327.5 KB   7.7%    2.7 MB  11.5%    8.5X
    CHROM                   1.9 KB   0.0%    2.2 MB   9.4% 1200.4X
    TXT_HEADER               931 B   0.0%     940 B   0.0%    1.0X
    Other                    804 B   0.0%  938.1 KB   3.9% 1194.8X
    TOTAL                   4.2 MB 100.0%   23.6 MB 100.0%    5.7X


.. _benchmark-LOCS1:

**Illumina s.locs file**

::

    LOCS file: s.locs
    Clusters: 4,091,904   Dictionaries: 3   Vblocks: 1 x 512 MB  Sections: 9
    Genozip version: 13.0.5
    Date compressed: 2021-11-28 16:51:17 Cen. Australia Daylight Time
    Command line: genozip --best --stats s.locs

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    X                       3.0 KB  79.9%   15.6 MB  50.0% 5410.8X
    TXT_HEADER               360 B   9.5%      12 B   0.0%    0.0X
    Y                        239 B   6.3%   15.6 MB  50.0% 68483.8X
    Other                    162 B   4.3%         -   0.0%    0.0X
    TOTAL                   3.7 KB 100.0%   31.2 MB 100.0% 8646.4X

    Note: ZIP total file size excludes overhead of 2.2 KB


.. _benchmark-LOCS2:

**Illumina locs file (for example: s_1_1101.locs)**

::

    LOCS file: s_1_1101.locs
    Clusters: 812,166   Dictionaries: 3   Vblocks: 1 x 512 MB  Sections: 9
    Genozip version: 13.0.5
    Date compressed: 2021-11-28 16:57:58 Cen. Australia Daylight Time
    Command line: genozip --best --stats s_1_1101.locs

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    X                       2.6 MB  78.5%    3.1 MB  50.0%    1.2X
    Y                     720.1 KB  21.5%    3.1 MB  50.0%    4.4X
    TXT_HEADER               360 B   0.0%      12 B   0.0%    0.0X
    Other                    162 B   0.0%         -   0.0%    0.0X
    TOTAL                   3.3 MB 100.0%    6.2 MB 100.0%    1.9X

    