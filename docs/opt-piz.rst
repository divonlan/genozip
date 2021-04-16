.. option:: -c, --stdout  Send output to standard output instead of a file.

          |

.. option:: -f, --force  Force overwrite of the output file.

          |

.. option:: -^, --replace  Replace the source file with the result file rather than leaving it unchanged.

          |

.. option:: -o, --output output-filename.  Output to this filename.

          |

.. option:: -p, --password password.  Provide password to access file(s) that were compressed with --password.

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

