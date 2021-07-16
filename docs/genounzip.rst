.. highlight:: none

genounzip
=========
Uncompress files compressed with ``genozip``.

Usage: ``genounzip`` [options]... [files]...

One or more file names must be given

          |

**Examples** 

     | ``genounzip file1.vcf.genozip file2.sam.genozip``
     |
     | ``genounzip file.vcf.genozip --output file.vcf.gz``
     |
     | ``genounzip bound.vcf.genozip --prefix=new_directory/``

**Options**

.. option:: -e, --reference filename.  Load a reference file prior to decompressing. Required only for files compressed with --reference

     | Note: this is equivalent of setting the environment variable $GENOZIP_REFERENCE with the reference filename.
     |

.. include:: opt-piz.rst

.. option:: -z, --bgzf level.  Compress the output to the BGZF format (.gz extension) using libdeflate at the compression level specified by the argument. Argument specifies the compression level from 0 (no compression) to 12 (best yet slowest compression). If you are not sure what value to choose - 6 is a popular option. Note: by default (absent this option) genozip will attempt to re-create the same BGZF compression as in the original file. Whether genozip succeeds in re-creating the exact same BGZF compression ratio depends on the compression library used by the application that generated the original file.

          |

.. option:: -u --prefix <prefix> Specify a prefix that is added to each file component name. A prefix may include a directory.

          |

.. option:: -m, --md5  Show the digest of the decompressed file - MD5 if the file was compressed with --md5 or --test and Adler32 if not. Note: for compressed files (e.g. myfile.vcf.gz) the digest calculated is that of the original uncompressed file. 

          |

.. option:: -t, --test  Decompress in memory (i.e. without writing the decompressed file to disk) and use the digest (MD5 or Adler32) to verify that the resulting decompressed file is identical to the original file.

          |

.. option:: -^, --replace  Replace the source file with the result file rather than leaving it unchanged.

          |


.. include:: opt-quiet.rst
.. include:: opt-threads.rst
.. include:: opt-stats.rst

.. option:: --no-PG.  (SAM BAM) When translating the file Genozip normally adds information in a @PG header line. With this option it doesn't.

**General options**

.. include:: opt-help.rst

