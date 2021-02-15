.. option:: -c, --stdout  Send output to standard output instead of a file.

          |

.. option:: -f, --force  Force overwrite of the output file.

          |

.. option:: -z, --bgzf level.  Compress the output to the BGZF format (.gz extension) using libdeflate at the compression level specified by the argument. Argument specifies the compression level from 0 (no compression) to 12 (best yet slowest compression). If you are not sure what value to choose - 6 is a popular option. Note: by default (absent this option) genozip will attempt to re-create the same BGZF compression as in the original file. Whether genozip succeeds in re-creating the exact same BGZF compression ratio depends on the compression library used by the application that generated the original file.

          |

.. option:: -^, --replace  Replace the source file with the result file rather than leaving it unchanged.

          |

.. option:: -o, --output output-filename.  Output to this filename.

          |

.. option:: -p, --password password. Provide password to access file(s) that were compressed with ``--password``.

          |

.. option:: -x, --index  Create an index file alongside the decompressed file. The index file is created as described:

          |

     ============ ============
     *Data type*  *Tool used*
     SAM/BAM      ``samtools index``
     FASTQ        ``samtools faidx``
     FASTA        ``samtools faidx``
     VCF          ``bcftools index``
     Other types  ``--index`` not supported
     ============ ============

|

