.. _losslessness:

Losslessness
============

Genozip compression is `lossless <https://en.wikipedia.org/wiki/Lossless_compression>`_ relative to the underlying data, which means that the data reconstructed during uncompression is **exactly identical** to the original data compressed.

**Verification of Losslessness**

When uncompressing with ``genounzip`` or ``genocat``, these tools verify that the reconstructed data is exactly identical to the source data, using a digest. Some exceptions apply. See here: :ref:`digest`.

**Exceptions to Losslessness**

1. There are some cases in which you may request ``genozip`` to change the source data before compressing it. In these cases, the digest is not calculated. These cases are:

- Using ``--optimize`` or any of the ``--optimize-*`` options 
- Using ``--GL-to-PL`` or ``--GP-to-PP``
- Generating a dual-coordinates file with ``--chain``
- Compressing a Luft file (a lifted-over dual-coordinates file)
- Compressing a SAM or BAM file with ``--taxid``

        |

2. Compressing already-compressed files (BGZF): Many tools that generate genomic files use `BGZF <https://www.htslib.org/doc/bgzip.html#BGZF_FORMAT>`_ to compress them. These typically have a .gz extension (although not all files a .gz extension are compressed with BGZF), and BGZF compression is also used internally in .bam and .bcf files. When ``genozip`` compresses a BGZF-compressed file, it first decompresses it to recover the original underlying data, and then compresses the data with Genozip. Likewise, when ``genounzip`` uncompresses a .genozip file, it recompresses the data back to BGZF if the original file was compressed with BGZF. Genozip records the parameters of the BGZF compression in the .genozip file (such as estimated compression level), and ``genounzip`` attempts to recompress to BGZF format using the same parameters. However, since there are many BGZF compression libraries, each with dozens of parameter combinations, it is possible that ``genounzip`` BGZF compression will achieve a slightly different compression level than the BGZF compression of the original file, resulting in the final BGZF-compressed file differing than the original BGZF file. However, the underlying data is still exactly identical and verified with the *digest*. 

        |

3. Compressing already-compressed files (.bz2 .xz and non-BGZF .gz): ``genozip`` is capable of compressing files that are already compressed with these methods. It does so by first uncompressing the file, and then recompressing with Genozip. However, in these cases, ``genounzip`` does not recompress the files back to .bz2 .xz or .gz.
   