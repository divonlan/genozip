genocat
=======
Display contents or metadata of a file compressed with ``genozip``.

Usage: ``genocat`` [options]... [files]...

One or more file names must be given.

**Reference-file related options**

.. option:: -e, --reference filename.  Load a reference file prior to decompressing. Required only for files compressed with --reference. When no non-reference file is specified display the reference data itself (typically used in combination with --regions).

          |

.. option:: -E, --REFERENCE filename.  With no non-reference file specified. Display the reverse complement of the reference data itself. Typically used in combination with --regions.

          |

.. option:: --show-reference  Show the name and MD5 of the reference file that needs to be provided to uncompress this file.

          |

**Subsetting (aka filtering) options (options resulting in modified display of the data)**

.. option:: --downsample rate[,shard].  Show only one in every <rate> lines (or reads in the case of FASTQ), optional <shard> parameter indicates which of the shards is shown. Other subsetting options (if any) will be applied to the surviving lines only.

          |

.. option:: --component component-number.  View a specific component of a genozip file. <component-number> is the number of the component as it appears in the genols list - the first component being number 1.

          |
.. option:: --interleaved  For FASTQ data compressed with --pair: Show every pair of paired-end FASTQ files with their reads interleaved: first one read of the first file ; then a read from the second file ; then the next read from the first file and so on.

          |

.. option:: -r, --regions [^]chr|chr:pos|pos|chr:from-to|chr:from-|chr:-to|from-to|from-|-to|from+len[,...].  (FASTA SAM/BAM GVF 23andMe Chain) Show one or more regions of the file. Examples:

   ============================================== ======================================
   ``genocat myfile.vcf.genozip -r 22:1000-2000`` Positions 1000 to 2000 on contig 22
   ``genocat myfile.sam.genozip -r 22:1000+151``  151 bases, starting pos 1000, on contig 22
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
   | *Note*: For Chain files this applies to the source contig (qName).
   |

.. option:: -s, --samples [^]sample[,...].  (VCF) Show a subset of samples (individuals). Examples:

   ================================================== ======================================
   ``genocat myfile.vcf.genozip -s HG00255,HG00256``  show two samples
   ``genocat myfile.vcf.genozip -s ^HG00255,HG00256`` show all samples except these two
   ================================================== ======================================
   
   | *Note*: This does not change the INFO data (including the AC and AN tags).
   |
   | *Note*: Sample names are case-sensitive.
   |
   | *Note*: Multiple ``-s`` arguments may be specified - this is equivalent to chaining their samples with a comma separator in a single argument.
   |

.. option:: --FLAG {+-^}value.  (SAM BAM) Filter lines based on the FLAG value: <value> is a decimal or hexadecimal value and should be prefixed by + - or ^: 

==  =======================================================================
\+  INCLUDES lines in which ALL flags in *value* are set in the line's FLAG
\-  INCLUDES lines in which NO flags in *value* are set in the line's FLAG
^   EXCLUDES lines in which ALL flags in *value* are set in the line's FLAG
==  =======================================================================

| *Example*: --FLAG -192 includes only lines in which neither FLAG 64 nor 128 are set. This can also be expressed as --FLAG -0xC0 
|
| More information of FLAGs can be found in section 1.4 of the `SAM spec <https://samtools.github.io/hts-specs/SAMv1.pdf>`_.

   |

.. option:: --MAPQ [^]value.  (SAM BAM) Filter lines based on the MAPQ value: INCLUDE (or EXCLUDE if <value> is prefixed with ^) lines with a MAPQ greater or equal to <value> 
   
          |

.. option:: -g, --grep string.  Show only lines (FASTA: sequences ; FASTQ: reads ; CHAIN: sets) in which <string> is a case-sensitive substring of the lines (FASTA and FASTQ: description). This does not affect showing the file header.

          |

.. option:: -n, --lines [first]-[last] or [first].  Show a certain range of lines. <first> and <last> are numbers of lines in the file (starting from 1). 

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

.. option:: --head [num_lines].  Show a certain number of lines from the start of the file.

          |

.. option:: --tail [num_lines].  Show a certain number of lines from the end of the file.

          |

.. option:: --iupac [^]value.  (SAM BAM FASTQ) Filter lines based on the IUPAC characters of the sequence data. Examples:
   
| *Examples*: 

========================== ===============================================================================
``genocat --iupac ACGTN``  displays only lines in which all characters of the SEQ are one of A,C,G,T,N
``genocat --iupac ^ACGTN`` displays only lines in which NOT all characters of the SEQ are one of A,C,G,T,N
========================== ===============================================================================

| Note: In SAM/BAM, all lines missing a sequence (i.e. SEQ=*) are included in positive iupac filters (the first example above) and excluded in negative ones.
| Note: The full list of IUPAC chacacters is here: `IUPAC codes <https://www.bioinformatics.org/sms/iupac.html>`_


.. option:: -K, --kraken filename. Load a .kraken.genozip file for use with --taxid. 

| See: :ref:`kraken`

.. option:: -k, --taxid [^]taxid[+0].  Show only lines than match the Taxonomy ID <taxid>. ^ for a negative search. +0 means <taxid> AND unclassified. Requires either using in combination with --kraken or for the file to have been compressed with genozip --kraken.

| See: :ref:`kraken`

.. option:: -G, --drop-genotypes.  (VCF) Output the data without the samples and FORMAT column.
   
          |

.. option:: -H, --no-header.  Don't output the header lines.

          |

.. option:: -1, --header-one.  (VCF FASTA) VCF: Output only the last line on the header (the line with the field and sample names). FASTA: Output the sequence name up to the first space or tab.

          |

.. option:: --header-only.  Output only the header lines.

          |

.. option:: --GT-only.  (VCF) Within samples output only genotype (GT) data - dropping the other subfields.

          |

.. option:: --unsorted.  If a file contains a "reconstruction plan" (see genozip --sort) the file will be displayed sorted by default. --unsorted overrides this behavior and shows the file in its unsorted form. This is useful if the file was highly unsorted causing sorting during genocat to consume a lot of memory.

          |

.. option:: --sequential.  (FASTA) Output in sequential format - each sequence in a single line.
   
          |

**Analysis options**

.. option:: --show-chroms.  (VCF SAM BAM FASTA GVF 23andMe) List the names of the chromosomes (or contigs) included in the file. Altnernative option names: --show-contigs --list-chroms --list-contigs.
    
          |

.. option:: --sex.  (SAM BAM FASTQ) Determine whether a SAM/BAM is a Male or a Female. Limitations when using on FASTQ. See "Sex assignment" for  details.
    
          |

.. option:: --coverage[=all|=one].  (SAM BAM FASTQ) Shows the coverage and depth of each contig. Approximate values when using on FASTQ. "Coverage and Depth" for details. 
    
          |

.. option:: --idxstats.  (SAM BAM FASTQ) Shows the count of mapped and unmapped reads by contig. Approximate values when using on FASTQ. Same output format as samtools idxstats. 
    
          |

.. option:: --count.  Rather than displaying the file content just report the number of lines (FASTQ: reads ; CHAIN: sets) (excluding the header) that would have been displayed. Useful in combination with filtering options.
    
          |

.. include:: opt-translation.rst

**General options**

.. include:: opt-piz.rst

.. option:: -z, --bgzf level.  Compress the output to the BGZF format (.gz extension). Note that by default genocat does not compress with BGZF (except for BAM that is compressed with level 1). Use this option if downstream tools require it.

.. include:: opt-quiet.rst

.. option:: --no-PG.  (VCF SAM BAM) When modifying the data in a file using genocat Genozip normally adds information about the modification in the file header. With this option it doesn't.
    
|

.. include:: opt-threads.rst
.. include:: opt-stats.rst

.. option:: --show-kraken[=INCLUDED|EXCLUDED]. In combination with --taxid reports whether each line is included or excluded. =INCLUDED or =EXCLUDED reports only a subset of lines accordingly. Combine with --count for a fast report without display the file itself. 

| See: :ref:`kraken`
                        

.. option:: --validate[=valid]  Validates that the file(s) are valid genozip files. By default reports files that are invalid. With --validate=valid reports files that are valid, and if run on a single exit code indicates validity.
    
          |

.. include:: opt-help.rst
