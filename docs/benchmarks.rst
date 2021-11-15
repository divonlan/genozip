.. _benchmarks:

Benchmarks
==========

Below are some benchmarks measuring Genozip's performance on a variety of file types. In the case of BAM, CRAM, FASTQ and VCF, Genozip is shown here compressing already-compressed files.

To see more, peer-reviewed, benchmarks, see :ref:`Publications <publications>`.

As can be appreciated from the *detailed compression reports* following the table, the compression gains strongly depend on the specific fields contained in a file, and can therefore vary significantly. 

================================ =========== ========= ========= ============ ================================= ==============================
Details                          Type        Size      .genozip  Genozip gain Data                              Link to data
================================ =========== ========= ========= ============ ================================= ==============================
:ref:`details<benchmark-BAM1>`   BAM         41 GB     10 GB     4.1X         30x Illumina NovaSeq + Dragen     Unpublished
:ref:`details<benchmark-BAM2>`   BAM         53.8 GB   30.2 GB   1.8X         PacBio CLR (mapped)               `1000 Genozip Project <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20131209 na12878 pacbio/si/NA12878.pacbio.bwa-sw.20140202.bam>`_
:ref:`details<benchmark-CRAM1>`  CRAM        15 GB     9.8 GB    1.5X         Illumina NovaSeq + bwa mem        `The European Bioinformatics Institute <ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram>`_
:ref:`details<benchmark-FASTQ1>` FASTQ.gz    61.2 GB   14.2 GB   4.3X         30x Illumina NovaSeq (R1+R2)      Unpublished
:ref:`details<benchmark-VCF1>`   VCF.gz      241.8 MB  85.8 MB   2.8X         30x Illumina NovaSeq + Dragen     Unpublished
:ref:`details<benchmark-VCF2>`   VCF.gz      27 GB     7.8 GB    3.3X         3202 samples                      `1000 Genome Project <ftp://ftp=trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/20201028_CCDG_14151_B01_GRM_WGS_2020=08=05_chr22.recalibrated_variants.vcf.gz>`_
:ref:`details<benchmark-FASTA1>` FASTA       1.2 GB    5.9 MB    212.2X       Covid-19 multi-FASTA              `coronavirus.innar.com <https://coronavirus.innar.com/coronavirus.unwrapped.fasta.zip>`_
:ref:`details<benchmark-GFF1>`   GFF3        714 KB    32 KB     22.2X        GRCh38 issues                     `NCBI <https://ftp.ncbi.nlm.nih.gov/pub/grc/human/GRC/Issue_Mapping/GRCh38.p9_issues.gff3>`_
:ref:`details<benchmark-ME1>`    23andMe     23.6 MB   4.2 MB    5.7X         Consumer DNA test "raw data"      Unpublished
================================ =========== ========= ========= ============ ================================= ==============================

Notes:

    - The tests were conducted with the ``--best`` option. For BAM, CRAM, FASTQ files, the ``--reference`` option was used to specify the appropriate reference file. For FASTQ, the ``--pair`` was used.
    
  
**Detailed compression reports**

*The following reports can be produced during compression with* ``genozip --stats`` *or after compression with* ``genocat --stats <myfile>.genozip``.

.. _benchmark-BAM1:

**BAM - 30x Illumina NovaSeq + Dragen:**

::

    BAM file: <redacted>.bam
    Reference: data/hs37d5.ref.genozip
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
    Reference: data/hs37d5.ref.genozip
    Alignments: 25,968,256   Dictionaries: 159   Vblocks: 215 x 512 MB  Sections: 9624
    Sorting: Sorted by POS
    Read name style: PacBio-Range
    Genozip version: 13.0.2 conda
    Date compressed: 2021-11-09 18:10:37 ACDT

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    QUAL                   21.6 GB  71.5%   41.1 GB  38.3%    1.9X
    SEQ                     3.7 GB  12.3%   20.6 GB  19.2%    5.6X
    CIGAR                   2.4 GB   8.0%   22.6 GB  21.1%    9.4X
    SA:Z                    2.3 GB   7.6%   19.1 GB  17.8%    8.4X
    QNAME                 103.6 MB   0.3%    1.9 GB   1.8%   18.7X
    AS:i                   32.0 MB   0.1%   38.0 MB   0.0%    1.2X
    XS:i                   31.1 MB   0.1%   29.2 MB   0.0%    0.9X
    POS                    24.5 MB   0.1%   99.1 MB   0.1%    4.0X
    Other                  11.6 MB   0.0%  520.0 MB   0.5%   44.9X
    MAPQ                   10.6 MB   0.0%   24.8 MB   0.0%    2.3X
    FLAG                    7.1 MB   0.0%   49.5 MB   0.0%    6.9X
    PNEXT                  50.9 KB   0.0%   99.1 MB   0.1% 1991.1X
    NM:i                   49.3 KB   0.0%   29.9 MB   0.0%  622.6X
    RNAME                  11.4 KB   0.0%   99.1 MB   0.1% 8877.3X
    RNEXT                   9.5 KB   0.0%   99.1 MB   0.1% 10686.5X
    TXT_HEADER              4.0 KB   0.0%   19.7 KB   0.0%    5.0X
    TLEN                      83 B   0.0%   99.1 MB   0.1% 1251482.2X
    RG:Z                      56 B   0.0%  396.2 MB   0.4% 7419501.5X
    PG:Z                      55 B   0.0%  371.5 MB   0.3% 7082251.5X
    BAM_BIN                   43 B   0.0%   49.5 MB   0.0% 1207825.9X
    GENOZIP vs BGZF        30.2 GB 100.0%   53.8 GB 100.0%    1.8X
    GENOZIP vs TXT         30.2 GB 100.0%  107.3 GB 100.0%    3.6X

.. _benchmark-CRAM1:

**CRAM - Illumina NovaSeq + bwa mem**

::

    SAM file: NA12878.final.cram
    Reference: data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
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

.. _benchmark-FASTQ1:

**FASTQ - 30x Illumina NovaSeq**

::

    FASTQ files (paired): <redacted>_R1_001.fastq.gz <redacted>_R2_001.fastq.gz
    Reference: GRCh38_full_analysis_set_plus_decoy_hla.ref.genozip
    Sequences: 860,000,926   Dictionaries: 23   Vblocks: 590 x 512 MB  Sections: 11689
    Read name style: Illumina-fastq
    Genozip version: 13.0.2 conda
    Date compressed: 2021-11-09 16:57:31 ACDT

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    QUAL                    7.7 GB  54.1%  118.7 GB  40.3%   15.4X
    SEQ                     5.2 GB  36.8%  118.7 GB  40.3%   22.7X
    DESC                    1.3 GB   9.1%   52.4 GB  17.8%   40.4X
    Other                  73.3 KB   0.0%    4.8 GB   1.6% 68760.6X
    LINE3                  23.6 KB   0.0%         -   0.0%    0.0X
    TXT_HEADER               696 B   0.0%         -   0.0%    0.0X
    GENOZIP vs BGZF        14.2 GB 100.0%   61.2 GB 100.0%    4.3X
    GENOZIP vs TXT         14.2 GB 100.0%  294.7 GB 100.0%   20.7X


.. _benchmark-VCF1:

**VCF - 30x Illumina NovaSeq + Dragen**

::

    VCF file: <redacted>.vcf.gz
    Samples: 1   Variants: 3,866,255   Dictionaries: 249   Vblocks: 3 x 512 MB  Sections: 400
    Genozip version: 13.0.2 conda
    Date compressed: 2021-11-09 16:32:06 ACDT

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    INFO/VQSLOD            11.1 MB  13.0%   26.1 MB   2.4%    2.3X
    FORMAT/GP               8.6 MB  10.0%   68.5 MB   6.3%    8.0X
    POS                     5.0 MB   5.9%   33.9 MB   3.1%    6.8X
    INFO/R2_5P_bias         4.6 MB   5.3%   12.8 MB   1.2%    2.8X
    QUAL                    4.6 MB   5.3%   20.1 MB   1.8%    4.4X
    INFO/SOR                4.4 MB   5.1%   18.1 MB   1.7%    4.1X
    INFO/ReadPosRankSum     3.7 MB   4.3%   11.5 MB   1.1%    3.1X
    FORMAT/AD               3.7 MB   4.3%   16.7 MB   1.5%    4.6X
    FORMAT/SB               3.6 MB   4.2%   32.7 MB   3.0%    9.1X
    FORMAT/F2R1             3.6 MB   4.2%   14.5 MB   1.3%    4.1X
    FORMAT/F1R2             3.6 MB   4.1%   14.5 MB   1.3%    4.1X
    INFO/QD                 3.6 MB   4.1%   14.3 MB   1.3%    4.0X
    INFO/MQRankSum          3.2 MB   3.7%   11.6 MB   1.1%    3.6X
    INFO/MQ                 3.1 MB   3.7%   15.6 MB   1.4%    5.0X
    FORMAT/MB               3.1 MB   3.6%   32.8 MB   3.0%   10.6X
    FORMAT/PL               2.7 MB   3.2%   28.1 MB   2.6%   10.2X
    INFO/FS                 2.7 MB   3.2%    9.9 MB   0.9%    3.6X
    FORMAT/AF               2.6 MB   3.0%   12.3 MB   1.1%    4.8X
    INFO/DP                 2.4 MB   2.8%    7.3 MB   0.7%    3.0X
    REF+ALT                 1.5 MB   1.8%   14.8 MB   1.4%    9.7X
    CHROM                   1.4 MB   1.6%    9.1 MB   0.8%    6.5X
    INFO/FractionInforma  953.3 KB   1.1%    6.7 MB   0.6%    7.2X
    FORMAT/GQ             679.8 KB   0.8%    8.2 MB   0.8%   12.3X
    FORMAT/GT             496.3 KB   0.6%   11.1 MB   1.0%   22.8X
    INFO                  289.0 KB   0.3%  336.8 MB  31.0% 1193.7X
    INFO/AF               282.4 KB   0.3%    8.2 MB   0.8%   29.6X
    Other                 183.1 KB   0.2%   44.7 MB   4.1%  249.9X
    FORMAT                181.8 KB   0.2%  152.5 MB  14.0%  858.8X
    FORMAT/PS              76.7 KB   0.1%    3.5 MB   0.3%   46.9X
    FILTER                 44.3 KB   0.1%   19.4 MB   1.8%  448.8X
    INFO/AC                 2.9 KB   0.0%    3.7 MB   0.3% 1292.4X
    TXT_HEADER              2.8 KB   0.0%    9.0 KB   0.0%    3.2X
    COORDS                   476 B   0.0%         -   0.0%    0.0X
    INFO/LOD                 413 B   0.0%     364 B   0.0%    0.9X
    FORMAT/PRI               274 B   0.0%   48.0 MB   4.4% 183538.6X
    FORMAT/DP                 96 B   0.0%    7.3 MB   0.7% 80128.0X
    INFO/AN                   83 B   0.0%    3.7 MB   0.3% 46580.6X
    ID                        42 B   0.0%    7.4 MB   0.7% 184107.4X
    GENOZIP vs BGZF        85.8 MB 100.0%  241.8 MB 100.0%    2.8X
    GENOZIP vs TXT         85.8 MB 100.0%    1.1 GB 100.0%   12.7X


.. _benchmark-VCF2:

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


.. _benchmark-FASTA1:

**FASTA - Covid-19 multi-FASTA**

::

    FASTA file: coronavirus.unwrapped.fasta
    Lines: 89,914   Dictionaries: 11   Vblocks: 79 x 16 MB  Sections: 586
    Sequence type: Nucleotide bases
    Genozip version: 13.0.2 conda
    Date compressed: 2021-11-09 19:26:11 ACDT

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    NONREF                  5.8 MB  97.1%    1.2 GB 100.0%  218.4X
    Other                 142.6 KB   2.4%   87.7 KB   0.0%    0.6X
    DESC                   33.9 KB   0.6%  526.3 KB   0.0%   15.5X
    TXT_HEADER               348 B   0.0%         -   0.0%    0.0X
    TOTAL                   5.9 MB 100.0%    1.2 GB 100.0%  212.2X


.. _benchmark-GFF1:

**GFF3 - GRCh38 issues**

::

    GFF3 file: https://ftp.ncbi.nlm.nih.gov/pub/grc/human/GRC/Issue_Mapping/GRCh38.p9_issues.gff3
    Sequences: 5,256   Dictionaries: 41   Vblocks: 1 x 16 MB  Sections: 36
    Genozip version: 13.0.2 conda
    Date compressed: 2021-11-09 19:21:06 ACDT

    Sections (sorted by % of genozip file):
    NAME                   GENOZIP      %      TXT       %   RATIO
    END                    10.1 KB  31.6%   39.9 KB   5.6%    3.9X
    START                   9.1 KB  28.4%   33.5 KB   4.7%    3.7X
    Name                    2.9 KB   9.0%   34.6 KB   4.8%   12.0X
    SEQID                   2.6 KB   8.0%   67.2 KB   9.4%   26.0X
    Other                   2.1 KB   6.4%         -   0.0%    0.0X
    fixVersion              1.5 KB   4.6%   42.5 KB   6.0%   29.0X
    type                    1.0 KB   3.2%   54.1 KB   7.6%   52.1X
    affectVersion            964 B   2.9%   34.6 KB   4.8%   36.7X
    status                   717 B   2.2%   44.2 KB   6.2%   63.2X
    chr                      392 B   1.2%    8.7 KB   1.2%   22.8X
    TXT_HEADER               389 B   1.2%      41 B   0.0%    0.1X
    ATTRS                    151 B   0.5%  266.9 KB  37.4% 1810.0X
    TYPE                      47 B   0.1%   35.9 KB   5.0%  782.8X
    SOURCE                    44 B   0.1%   20.5 KB   2.9%  477.8X
    SCORE                     42 B   0.1%   10.3 KB   1.4%  250.3X
    STRAND                    42 B   0.1%   10.3 KB   1.4%  250.3X
    PHASE                     42 B   0.1%   10.3 KB   1.4%  250.3X
    COMMENT                   41 B   0.1%         -   0.0%    0.0X
    TOTAL                  32.1 KB 100.0%  713.6 KB 100.0%   22.2X


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

