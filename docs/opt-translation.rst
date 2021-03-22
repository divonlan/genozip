**Translation options (convertion from one format to another)**

.. option:: --bam  (SAM and BAM only) Output as BAM. Note: this option is implicit if --output specifies a filename ending with .bam

          |

.. option:: --sam  (SAM and BAM only) Output as SAM. This option is the default in genocat on SAM and BAM data.

          |

.. option:: --fastq  (SAM and BAM only) Output as FASTQ. The alignments are outputted as FASTQ reads in the order they appear in the SAM/BAM file. Alignments with FLAG 16 (reverse complimented) have their SEQ reverse complimented and their QUAL reversed. Alignments with FLAG 4 (unmapped) or 256 (secondary) are dropped. Alignments with FLAG 64 (or 128) (the first (or last) segment in the template) have a '1' (or '2') added after the read name. Usually (if the original order of the SAM/BAM file has not been tampered with) this would result in a valid interleaved FASTQ file. Note: this option is implicit if --output specifies a filename ending with .fq[.gz] or .fastq[.gz]

          |

.. option:: --bcf  (VCF only) Output as BCF. Note: bcftools needs to be installed for this option to work.

          |

.. option:: --phylip  (FASTA only) Output a Multi-FASTA in Phylip format. All sequences must be the same length.

          |

.. option:: --fasta  (Phylip only) Output as Multi-FASTA.

          |

.. option:: --vcf  (23andMe only) Output as VCF. --vcf must be used in combination with --reference to specify the reference file as listed in the header of the 23andMe file (usually this is GRCh37). Note: INDEL genotypes ('DD' 'DI' 'II') as well as uncalled sites ('--') are discarded.

          |

