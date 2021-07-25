.. highlight:: none

genocat
=======
Display contents or metadata of a file compressed with ``genozip``.

Usage: ``genocat`` [options]... [files]...

One or more file names must be given.

**General options**

.. include:: opt-piz.rst

.. option:: -z, --bgzf level.  Compress the output to the BGZF format (.gz extension). Note that by default genocat does not compress with BGZF (except for BAM that is compressed with level 1). Use this option if downstream tools require it.

.. include:: opt-quiet.rst

.. include:: opt-threads.rst
.. include:: opt-stats.rst                        

.. option:: --validate[=valid]  Validates that the file(s) are valid genozip files. By default reports files that are invalid. With --validate=valid reports files that are valid, and if run on a single exit code indicates validity.
    
          |

.. include:: opt-help.rst

.. option:: -e, --reference filename.  Load a reference file prior to decompressing. Used for files compressed with --reference. 

   | Note: this is equivalent of setting the environment variable $GENOZIP_REFERENCE with the reference filename.
   |

.. option:: --show-reference  Show the name and MD5 of the reference file that needs to be provided to uncompress this file.

          |

**Subsetting options (options resulting in modified display of the data)**

.. option:: --downsample rate[,shard].  Show only one in every <rate> lines (or reads in the case of FASTQ). The optional <shard> parameter indicates which of the shards is shown. Other subsetting options (if any) will be applied to the surviving lines only.

          |

.. option:: --component component-number.  View a specific component of a genozip file. <component-number> is the number of the component as it appears in the genols list - the first component being number 1.

          |

.. option:: -r, --regions [^]chr|chr:pos|pos|chr:from-to|chr:from-|chr:-to|from-to|from-|-to|from+len[,...].  (VCF SAM/BAM GFF3/GVF FASTA 23andMe Chain Reference) Show one or more regions of the file. 

   |
   | *Examples*: 

   ============================================== ======================================
   ``genocat myfile.vcf.genozip -r 22:1000-2000`` Positions 1000 to 2000 on contig 22
   ``genocat myfile.ref.genozip -r 22:2000-1000`` Reverse complement of positions 1000 to 2000 on contig 22 (reference file only)
   ``genocat myfile.sam.genozip -r 22:1000+151``  151 bases, starting pos 1000, on contig 22
   ``genocat myfile.ref.genozip -r 22:1000-151``  Reverse complement of 151 bases, from 1000 to 850, on contig 22 (reference file only)
   ``genocat myfile.vcf.genozip -r -2000,2500-``  Two ranges on all contigs
   ``genocat myfile.sam.genozip -r chr21,chr22``  Contigs chr21 and chr22 in their entirety
   ``genocat myfile.vcf.genozip -r ^MT,Y``        All contigs, excluding MT and Y
   ``genocat myfile.vcf.genozip -r ^-1000``       All contigs, excluding positions up to 1000
   ``genocat myfile.fa.genozip  -r chrM``         Contig chrM
   ============================================== ======================================

   | *Note*: genozip files are indexed automatically during compression. There is no separate indexing step or separate index file.
   |
   | *Note*: Indels are considered part of a region if their start position is.
   |
   | *Note*: Multiple ``-r`` arguments may be specified - this is equivalent to chaining their regions with a comma separator in a single argument.
   |
   | *Note*: For FASTA and Chain files, only whole-contig regions are possible.
   |
   | *Note*: For Chain files this applies to the Primary contig (qName).
   |

.. option:: -R, --regions-file [^]filename.  (VCF SAM/BAM GFF3/GVF FASTA 23andMe Chain Reference) Show regions from a list in tab-separated file. To include all regions EXCEPT those in the file٫ prefix the filename with ^. If filename is - (or ^-) data is taken from stdin rather than a file.

|
| *Example of a valid file: The first two rows produce the same 100-base region, and the third row is a single base*:

::

   chr22	17000000	17000099
   chr22	17000000	+100
   chr22	17000000

.. option:: --grep string.  Show only lines (FASTA: sequences ; FASTQ: reads ; CHAIN: sets) in which <string> is a case-sensitive substring of the lines (FASTA: description). This does not affect showing the file header.

          |

.. option:: -g, --grep-w string.  Same as --grep, but restrict to whole words.

          |

.. option:: -n, --lines [first]-[last] or [first].  Show a certain range of lines. <first> and <last> are numbers of lines in the file (starting from 1). 

   |
   | *Examples*: 

================================ ========================================================
``genocat --lines 1000-2000``    displays the 1001 lines between 1000 and 2000
``genocat --lines=1000-``        displays all lines starting from 1000 (*optional* =)
``genocat -n -2000``             displays lines 1 to 2000 (``-n`` *instead of* ``--lines``)
``genocat -n 1000``              displays 10 lines starting from line 1000
================================ ========================================================


   | *Note on outputting as BAM*: The numbering excludes the BAM header. 
   | *Note on FASTQ*: The numbering is of reads rather than lines. 
   | *Note*: The entire file header is included if any part of it is.
   | *Note*: Line numbers are taken before any additional filters are applied.

.. option:: --head [num_lines].  Show <num_lines> lines from the start of the file.

          |

.. option:: --tail [num_lines].  Show <num_lines> lines from the end of the file.

          |

**VCF options**

.. option:: -s, --samples [^]sample[,...] or num_samples.  Show a subset of samples (individuals). No other fields (such as AF, AC) are updated.

   |
   | *Examples*:

   ================================================== ======================================
   ``genocat myfile.vcf.genozip -s HG00255,HG00256``  show two samples
   ``genocat myfile.vcf.genozip -s ^HG00255,HG00256`` show all samples except these two
   ``genocat myfile.vcf.genozip -s 5``                show the first 5 samples
   ================================================== ======================================
   
   | *Note*: This does not change the INFO data (including the AC and AN tags).
   |
   | *Note*: Sample names are case-sensitive.
   |
   | *Note*: Multiple ``-s`` arguments may be specified - this is equivalent to chaining their samples with a comma separator in a single argument.
   |

.. option:: -G, --drop-genotypes.  Output the data without the samples and FORMAT column. No other fields (such as AF, AC) are updated.
   
          |

.. option:: --GT-only.  Within samples output only genotype (GT) data - dropping the other subfields.

          |

.. option:: --snps-only.  Drops variants that are not a Single Nucleotide Polymorphism (SNP).

          |

.. option:: --indels-only.  Drops variants that are not Insertions or Deletions (INDEL).

          |

.. option:: --unsorted.  If a file contains a "reconstruction plan" (see genozip --sort) the file will be displayed sorted by default. --unsorted overrides this behavior and shows the file in its unsorted form. This is useful if the file was highly unsorted causing sorting during genocat to consume a lot of memory.

          |

.. option:: -1, --header-one.  Output only the last line on the header (the line with the field and sample names).
          
          |

.. option:: --bcf  Output as BCF. Note: bcftools needs to be installed for this option to work.

          |

.. option:: --luft.  Render a DVCF file in Luft coordinates (absent this option, a DVCF will be rendered in Primary coordinates).
      
   | See: :ref:`dvcf`
   |

.. option:: -y, --show-dvcf.  For each variant show its coordinate system (Primary or Luft or Both) and its oStatus. May be used with or without --luft.
   
   | See: :ref:`dvcf`
   |
          
.. option:: --show-ostatus.  Add oSTATUS to the INFO field. May be used with or without of --luft.
   
   | See: :ref:`dvcf`
   |
          
.. option:: --show-counts=o\$TATUS.  Show summary statistics of variant lift outcome.

   | See: :ref:`dvcf`
   |
      
.. option:: --show-counts=COORDS.  Show summary statistics of variant coordinates.

   | See: :ref:`dvcf`
   |

.. option:: --no-PG.  When modifying the data in a file using genocat Genozip normally adds a "##genozip_command" line to the VCF header. With this option it doesn't.

          |

.. option:: --gpos.  Replaces (CHROM,POS) with a coordinate in GPOS (Global POSition) terms. GPOS is a single genome-wide coordinate defined by a reference file, in which contigs appear in the order of the original FASTA data used to generate the reference file. Must be used in combination with --reference. The mapping of CHROM to GPOS can be viewed with "genocat --show-ref-contigs <reference-file.ref.genozip>". 

          |

**SAM and BAM options**

.. option:: --FLAG {+-^}value.  Filter lines based on the FLAG value: <value> is a decimal or hexadecimal value and should be prefixed by + - or ^: 

   ==  =======================================================================
   \+  INCLUDES lines in which ALL flags in *value* are set in the line's FLAG
   \-  INCLUDES lines in which NO flags in *value* are set in the line's FLAG
   ^   EXCLUDES lines in which ALL flags in *value* are set in the line's FLAG
   ==  =======================================================================

   | *Example*: --FLAG -192 includes only lines in which neither FLAG 64 nor 128 are set. This can also be expressed as --FLAG -0xC0 
   |
   | More information of FLAGs can be found in section 1.4 of the `SAM spec <https://samtools.github.io/hts-specs/SAMv1.pdf>`_.
   |

.. option:: --MAPQ [^]value.  Filter lines based on the MAPQ value: INCLUDE (or EXCLUDE if <value> is prefixed with ^) lines with a MAPQ greater or equal to <value> 
   
          |

.. option:: --bases [^]value.  Filter lines based on the IUPAC characters (bases) of the sequence data. 
   
   |
   | *Examples*: 

   ========================== ===============================================================================
   ``genocat --bases ACGTN``  displays only lines in which all characters of the SEQ are one of A,C,G,T,N
   ``genocat --bases ^ACGTN`` displays only lines in which NOT all characters of the SEQ are one of A,C,G,T,N
   ========================== ===============================================================================

   | Note: In SAM/BAM, all lines missing a sequence (i.e. SEQ=*) are included in positive --bases filters (the first example above) and excluded in negative ones.
   | Note: The list of IUPAC chacacters can be found here: `IUPAC codes <https://www.bioinformatics.org/sms/iupac.html>`_
   |

.. option:: --bam  Output as BAM. Note: this option is implicit if --output specifies a filename ending with .bam

          |

.. option:: --sam  Output as SAM. This option is the default in genocat on SAM and BAM data and is implicit if --output specifies a filename ending with .sam

          |

.. option:: --fastq  Output as FASTQ. Note: this option is implicit if --output specifies a filename ending with .fq or .fastq

   | see more details: :ref:`sam2fq`
   |

.. option:: --no-PG.  When modifying the data in a file using genocat Genozip normally adds information about the modification in the file header. With this option it doesn't.

          |

**FASTQ options**

.. option:: --interleaved[=both|either]  For FASTQ data compressed with --pair: Show every pair of paired-end FASTQ files with their reads interleaved: first one read of the first file ; then a read from the second file ; then the next read from the first file and so on. Optional argument 'both' (default) or 'either' determines whether both reads of a pair or only one is required for the pair to survive when combining with a subsetting option such as --grep.

          |

.. option:: --bases [^]value.  Filter lines based on the IUPAC characters (bases) of the sequence data (see SAM/BAM options).

          |

**FASTA options**

.. option:: -1, --header-one.  Output the sequence name up to the first space or tab.

          |

.. option:: -H, --no-header.  Don't output the header lines.

          |

.. option:: --header-only.  Output only the header lines.

          |

.. option:: --sequential.  Output in sequential format - each sequence in a single line.
   
          |

.. option:: --phylip  Output a Multi-FASTA in Phylip format. All sequences must be the same length.

          |

**Phylip options**

.. option:: --fasta  (Phylip only) Output as Multi-FASTA.

          |

**Reference file options**

.. option:: --reference <file> --regions <regions> [--header-only] View one or more regions of a reference file

   | Note: For reverse complement, use a reverse range, eg -r1000000-999995 or equivalently -r1000000-6
   | Note: --regions-file maybe used intead of --regions
   | Note: Combine with --no-header to suppress output of the chromosome name.
   | Note: Short forms of the options (eg -e instead of --reference) are fine too. 

          |

.. option:: --gpos. In combination with --reference and --regions or --regions-file - shows coordinates in GPOS (Global POSition) terms - a single genome-wide numeric coordinate - rather than (CHROM,POS).

          |

.. option:: --show-ref-contigs. Show the details of the reference file contigs. 

          |

.. option:: --show-ref-iupacs. Show non-ACTGN `IUPAC <http://www.bioinformatics.org/sms/iupac.html>`_ pseudo-bases in the reference file. 

          |

**Chain file options**

.. option:: --with-chr.  rewrites the qName - single-character or double-digit chromosome names are prefixed with "chr" and "MT" is rewritten as "chrM".

          |

.. option:: --show-chain  Show chain file alignments. 

          |

.. option:: --show-chain-contigs  Show the details of the chain file contigs. 

          |

**23andMe options**

.. option:: --vcf  Output as VCF. --vcf must be used in combination with --reference to specify the reference file as listed in the header of the 23andMe file (usually this is GRCh37). Note: Indel variants ('DD' 'DI' 'II') as well as uncalled sites ('--') are discarded.

          |

**Filtering using kraken data**

| Works for SAM, BAM, FASTQ and FASTA - See: :ref:`kraken`

.. option:: -K, --kraken filename. Load a .kraken.genozip file for use with --taxid. 

          |

.. option:: -k, --taxid [^]taxid[+0].  Show only lines than match the Taxonomy ID <taxid>. ^ for a negative search. +0 means <taxid> AND unclassified. Requires either using in combination with --kraken or for the file to have been compressed with genozip --kraken.

          |

.. option:: --show-kraken[=INCLUDED|EXCLUDED]. In combination with --taxid reports whether each line is included or excluded. =INCLUDED or =EXCLUDED reports only a subset of lines accordingly. Combine with --count for a fast report without display the file itself. 

          |

**Analysis options**

.. option:: --contigs.  (VCF SAM BAM FASTA GFF3/GVF 23andMe) List the names of the chromosomes (or contigs) included in the file. Altnernative option names: --list-chroms --chroms.
    
          |

.. option:: --sex.  (SAM BAM FASTQ) Determine whether a SAM/BAM is a Male or a Female. Limitations when using on FASTQ. See "Sex assignment" for  details.
    
          |

.. option:: --coverage[=all|=one].  (SAM BAM FASTQ) Shows the coverage and depth of each contig. Approximate values when using on FASTQ. "Coverage and Depth" for details. 
    
          |

.. option:: --idxstats.  (SAM BAM FASTQ) Shows the count of mapped and unmapped reads by contig. Approximate values when using on FASTQ. Same output format as samtools idxstats. 
    
          |

.. option:: --count.  Rather than displaying the file content just report the number of lines (FASTQ: reads ; CHAIN: sets) (excluding the header) that would have been displayed. Useful in combination with filtering options.
    
          |

