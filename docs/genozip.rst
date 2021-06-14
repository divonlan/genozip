genozip
=======
Compress files. 

``genozip`` can compress any file, but is optimally designed to compress the following file types: VCF/BCF, SAM/BAM/CRAM, FASTQ, FASTA, GVF, PHYLIP, Chain, Kraken and 23andMe

Usage: ``genozip`` [options]... [files or urls]...

One or more file names or URLs may be given, or if omitted, standard input is used instead. ``-`` means standard input.

Supported input file types, as recognized by their listed filename extension(s): 

   ======== ==========================================================
   *Type*   *Filename extensions*
   ======== ==========================================================
   FASTA    fasta, fa, faa, ffn, fnn, fna (possibly .gz .bgz .bz2 .xz)
   FASTQ    fastq, fq (possibly .gz .bgz .bz2 .xz)
   SAM      sam (possibly .gz .bgz .bz2 .xz)
   BAM      bam
   CRAM     cram
   VCF      vcf (possibly .gz .bgz .bz2 .xz)
   BCF      bcf (possibly .gz .bgz)
   GVF      gvf (possibly .gz .bgz .bz2 .xz)
   PHYLIP   phy (possibly .gz .bgz .bz2 .xz)
   Chain    chain (possibly .gz .bgz .bz2 .xz)
   Kraken   kraken (possibly .gz .bgz .bz2 .xz)
   23andMe  genome\*Full\*.txt (possibly zip)
   Generic  any other file (possibly .gz .bgz .bz2 .xz)
   ======== ==========================================================

Note: compressing .bcf, .cram or .xz files requires bcftools, samtools or xz, respectively, to be installed.

Examples: 

    ``genozip sample.bam``

    ``genozip sample.R1.fq.gz sample.R2.fq.gz --pair --reference hg19.ref.genozip -o sample.genozip``      

    ``genozip --optimize -password 12345 ftp://ftp.ncbi.nlm.nih.gov/file2.vcf.gz``

**Flags**

.. option:: -i, --input data-type.  <data-type> is one of the supported file extensions listed in the table above (eg bam vcf.gz fq.xz. See "genozip --help=input" for full list of accepted file types. This flag should be used when redirecting input data with a < or | or if the input file type cannot be determined by its file name.

                     |
                     
.. option:: -f, --force  Force overwrite of the output file or force writing the compressed data to standard output.

                     |
                     
.. option:: -^, --replace  Replace the source file with the result file rather than leaving it unchanged.

                     |
                     
.. option:: -o, --output output-filename.  This option can also be used to bind multiple input files into a single genozip file. genounzip will unbind the file back to its components while genocat will concatenate them. To bind files they must be of the same type (e.g. VCF or SAM) and if they are VCF files they must contain the same samples. genozip takes advantage of similarities between the input files so that the bound file is usually smaller than the combined size of individually compressed files.

                     |
                     
.. option:: --best  Best compression but slower than --fast mode. This is the default mode of genozip - this flag has no additional effect.

                     |
                     
.. option:: -F, --fast  Fast compression but lower compression ratio than --best. Files compressed with this option also uncompress faster. Compressing with this option also consumes less memory.

                     |
                     
.. option:: -p, --password password.  Password-protected - encrypted with 256-bit AES.

                     |
                     
.. option:: -m, --md5  Use MD5 (rather than the default Adler32) to calculate the digest of the original textual file. The MD5 digest is also viewable with genols. Note: for compressed files (e.g. myfile.vcf.gz or myfile.bam) the MD5 calculated is that of the original uncompressed textual file - myfile.vcf or myfile.sam respectively.
             
                     |
                     
.. option:: -I, --input-size file-size-in-bytes.  genozip configures its internal data structures to optimize execution speed based on the file size. When redirecting the input file with < or | genozip cannot determine its size and this might result in slower execution. This problem can be overcome by using this flag to inform genozip of the file size.

                     |
                     
.. option:: -t, --test  After compressing normally decompresss in memory (i.e. without writing the decompressed file to disk) - comparing the MD5 of the resulting textual decompressed file to that of the original textual file. This option also activates --md5.

                     |
                     
.. include:: opt-quiet.rst
.. include:: opt-threads.rst
                     
.. option:: -B, --vblock number.  Set the maximum size of data (between 1 and 2048 in megabytes) of the textual input data that a thread processes at any given time. By default genozip sets this value dynamically based on the characateristics of the file and it is reported in --show-stats (but capped at 32MB on Windows and MacOS). Smaller values will result in faster subsetting with genocat --regions and --grep. Larger values will result in better compression. Note that memory consumption of both genozip and genounzip is linear with the vblock value used for compression.

                     |
                     
.. option:: -e, --reference filename.  Use a reference file (filename extension .ref.genozip) - this is a FASTA file genozipped with the --make-reference option. The same reference needs to be provided to genounzip or genocat. While genozip is capabale of compressing without a reference it can utilize a reference file to improve compression of FASTQ SAM/BAM and VCF files. The improvement for FASTQ files is substantial; for SAM/BAM it may be significant; for VCF if it is significant only if REFALT content is a significant percentage of the zip content (see "% of zip" in --show-stats)

   | Note: this is equivalent of setting the environment variable $GENOZIP_REFERENCE with the reference filename.
   |
                     
.. option:: -E, --REFERENCE filename.  Similar to --reference except genozip copies the reference (or part of it) to the output file so there is no need to specify --reference in genounzip and genocat. Note on using with --password: the copy of the reference file stored in the compressed file is never encrypted.

                     |

.. option:: -K, --kraken filename. Incorporate the Taxonomy ID of each line into the file. For use with genocat --taxid. For SAM/BAM it also adds a TX:i field. 

| See: :ref:`kraken`
                        
.. include:: opt-stats.rst

.. option:: --register  Register (or re-register) a non-commericial license to use genozip.

                     |

**VCF-specific options**

.. option:: --chain  chain-file.  Lifts a VCF to be a dual-coordinate VCF (DVCF).

   | See: :ref:`dvcf`
   |

.. option:: --allow-ambiguous.  Used in combination with --chain - allow Insertion Deletion and Structural variants that Genozip isn't sure whether they are a REF⇆ALT switch or not. Absent this flag these variants would not be lifted and would remain Primary-only variants.

   | See: :ref:`dvcf`
   |

.. option:: --reject-ref-alt-switches.  Used in combination with --chain - Bi-allelic SNP variants with a REF⇆ALT switch are normally lifted and their INFO and FORMAT subfields modified to reflect the REF⇆ALT switch. With this option these variants are rejected from lifting and remain Primary-only variants.


   | See: :ref:`dvcf`
   |

.. option:: --show-counts=o\$TATUS.  Show summary statistics of variant lift outcome. This is set by default when using --chain.

   | See: :ref:`dvcf`
   |
      
.. option:: --show-counts=COORDS.  Show summary statistics of variant coordinates.

   | See: :ref:`dvcf`
   |

.. option:: --show-chain.  Used in combination with --chain - displays all chain file alignments.

   |

.. option:: --sort.  Causes genozip to generate a "reconstruction plan" that will allow genocat to show the file sorted. This is designed for mildly-unsorted files. If the file is highly unsorted this might result in genocat loading a big portion of the uncompressed file to memory (genocat --unsorted can be used to prevent sorting). This option is always set for dual-coordinates files unless overridden with --unsorted.

   |

.. option:: --unsorted.  Don't generate a "reconstruction plan".

   |

.. option:: --add-line-numbers.  Replaces the ID field in each variant with a sequential line number starting from 1.

   |

**FASTQ-specific options (ignored for other file types)**

.. option:: -2, --pair  Compress pairs of paired-end FASTQ files resulting in compression ratios better than compressing the files individually. When using this option every two consecutive files on the file list should be paired-end FASTQ files with an identical number of reads and consistent file names and --reference or --REFERENCE must be specified. The resulting genozip file is a bound file. To display it interleaved use genocat --interleaved. To unbind the genozip file back to its original FASTQ files use genounzip.

                     |
                     
**FASTA-specific options (ignored for other file types):**

.. option:: --make-reference  Compresss a FASTA file to be used as a reference in --reference or --REFERENCE.

                     |
                     
.. option:: --multifasta  All contigs in the FASTA file are variations of a the same contig (i.e. they are somewhat similar to each other). genozip uses this information to improve the compression.

                     |
                     
**Optimizing**

.. option:: -9, --optimize, --optimise  Modify the file in ways that are likely insignificant for analytical purposes but significantly improve compression and somewhat improve the speed of genocat --regions. This option activates all optimizations. Alternatively they can be activated individually as listed below.

   Note: files compressed with this option are NOT identical to the original file after decompression. For this reason, it is not possible to use this option in combination with --test or --md5

*VCF optimizations* 

.. option:: --optimize-sort  INFO subfields are sorted alphabetically.               

                              | Example: ``AN=21;AC=3`` -> ``AC=3;AN=21``
                              |


.. option:: --optimize-PL  PL data: Phred values of over 60 are changed to 60.     

                              | Example: ``0,18,270`` -> ``0,18,60``
                              |

.. option:: --optimize-GL  GL data: Numbers are rounded to 2 significant digits.   

                              | Example: ``-2.61618,-0.447624,-0.193264`` -> ``-2.6,-0.45,-0.19``
                              |

.. option:: --optimize-GP  GP data: Numbers are rounded to 2 significant digits as with GL.

                              |

.. option:: --optimize-VQSLOD  VQSLOD data: Number is rounded to 2 significant digits. 

                              | Example: ``-4.19494`` -> ``-4.2``

*SAM and BAM optimizations*
   
.. option:: --optimize-QUAL  The QUAL quality field and the secondary U2 quality field (if it exists) are modified to group quality scores into a smaller number of bins:

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

                     | This assumes a standard Sanger format of `Phred quality scores <https://en.wikipedia.org/wiki/Phred_quality_score>`_ 0->93 encoded in ASCII 33->126
                     | Note: this follows `Illumina's quality bins <https://sapac.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf>`_ for values up to Phred 39, and extends with additional similar bins for values of 40 and above common in some non-Illumina technologies.
                     | Example: ``LSVIHINKHK`` -> ``IIIIFIIIFI``
                     |

.. option:: --optimize-ZM  ZM:B:s data: negative Ion Torrent flow signal values are changed to zero, and positives are rounded to the nearest 10.  

                     Example: ``-20,212,427`` -> ``0,210,430``

*FASTQ optimizations*
   
.. option:: --optimize-DESC  Replaces the description line with @filename:read_number.
                     
                     | Example: ``@A00488:61:HMLGNDSXX:4:1101:1561:1000 2:N:0:CTGAAGCT+ATAGAGGC`` -> ``@sample.fq.gz:100`` (100 is the read sequential number within this fastq file)
                     |

.. option:: --optimize-QUAL  The quality data is optimized as described for SAM above.

*GVF optimizations*
   
.. option:: --optimize-sort  Attributes are sorted alphabetically.
                     
                     | Example: ``Notes=hi;ID=rs12`` -> ``ID=rs12;Notes=hi``
                     |

.. option:: --optimize-Vf  Variant_freq data: Number is rounded to 2 significant digits. 
                     
                     | Example: ``0.006351`` -> ``0.0064``
                     |

**Other actions - uses other than compressing a file**

.. option:: -d, --decompress  Same as running genounzip.
    
                     |

.. option:: -l, --list  Same as running genols.

                     |
                     
.. include:: opt-help.rst
                     
