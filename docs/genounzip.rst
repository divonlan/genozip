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

          |

.. include:: opt-piz.rst

.. option:: -u --prefix <prefix> Specify a prefix that is added to each file component name. A prefix may include a directory.

          |

.. option:: -m, --md5  Show the digest of the decompressed file - MD5 if the file was compressed with --md5 or --test and Adler32 if not. Note: for compressed files (e.g. myfile.vcf.gz) the digest calculated is that of the original uncompressed file. 

          |

.. option:: -t, --test  Decompress in memory (i.e. without writing the decompressed file to disk) and use the digest (MD5 or Adler32) to verify that the resulting decompressed file is identical to the original file.

          |

.. option:: -c, --stdout  (obsolete) use genocat instead          

          |

.. include:: opt-quiet.rst
.. include:: opt-threads.rst
.. include:: opt-stats.rst
.. include:: opt-translation.rst

.. option:: --no-PG.  (SAM BAM) When translating the file Genozip normally adds information in a @PG header line. With this option it doesn't.

**Other actions - uses other than uncompressing a file**

.. include:: opt-help.rst

