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

**Subsetting options (options resulting in modified display of the data)**

.. option:: --downsample rate[,shard].  Show only one in every <rate> lines (or reads in the case of FASTQ), optional <shard> parameter indicates which of the shards is shown. Other subsetting options (if any) will be applied to the surviving lines only.

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

.. option:: -g, --grep string.  (FASTQ FASTA) Show only records in which <string> is a case-sensitive substring of the description.

          |

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

.. option:: --list-chroms.  (VCF SAM BAM FASTA GVF 23andMe) List the names of the chromosomes (or contigs) included in the file.
    
          |

.. option:: --show-sex.  (SAM BAM) Determine whether a SAM/BAM is a Male or a Female. See "Sex assignment" use case.
    
          |

.. option:: --show-coverage[=all].  (SAM BAM) Shows the coverage and depth of each contig. Without =all it shows only contigs that are chromosomes and groups the other contigs under "Other contigs". See "Coverage and Depth" use case.
    
          |

.. option:: --show-coverage-chrom.  (SAM BAM) Same as --show-coverage but shows only contigs that are chromosomes.
    
          |

.. include:: opt-translation.rst

**General options**

.. include:: opt-piz.rst
.. include:: opt-quiet.rst

.. option:: --no-PG.  (VCF SAM BAM) When modifying the data in a file using genocat Genozip normally adds information about the modification in the file header. With this option it doesn't.
    
|

.. include:: opt-threads.rst
.. include:: opt-stats.rst

.. option:: --validate  Validates that the file(s) are valid genozip files.
    
          |

.. include:: opt-help.rst
