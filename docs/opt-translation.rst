**Translation options (convertion from one format to another)**

.. option:: --bam  (SAM and BAM only) Output as BAM. Note: this option is implicit if --output specifies a filename ending with .bam

          |

.. option:: --sam  (SAM and BAM only) Output as SAM. This option is the default in genocat on SAM and BAM data.

          |

.. option:: --fastq  (SAM and BAM only) Output as FASTQ. 

    | see more details: :ref:`sam2fq`
    |

.. option:: --bcf  (VCF only) Output as BCF. Note: bcftools needs to be installed for this option to work.

          |

.. option:: --phylip  (FASTA only) Output a Multi-FASTA in Phylip format. All sequences must be the same length.

          |

.. option:: --fasta  (Phylip only) Output as Multi-FASTA.

          |

.. option:: --vcf  (23andMe only) Output as VCF. --vcf must be used in combination with --reference to specify the reference file as listed in the header of the 23andMe file (usually this is GRCh37). Note: INDEL genotypes ('DD' 'DI' 'II') as well as uncalled sites ('--') are discarded.

          |

