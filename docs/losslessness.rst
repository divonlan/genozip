.. _losslessness:

Losslessness
============

Genozip compression is `lossless <https://en.wikipedia.org/wiki/Lossless_compression>`_ relative to the underlying data, which means that the data reconstructed during uncompression is **exactly identical** to the original data compressed.

**Verification of Losslessness**

When reading the data for compression, ``genozip`` calculates its *digest* using either `Adler32 <https://en.wikipedia.org/wiki/Adler-32>`_ (default) or `MD5 <https://en.wikipedia.org/wiki/MD5>`_ (if ``--md5`` is specified), and stores it along with the compressed data in the .genozip file.

After uncompressing a .genozip file with ``genounzip`` or ``genocat`` or when testing with ``genozip --test``, the *digest* of the reconstructed data is calculated and compared to the *digest* stored in the .genozip file. If it differs (this should never happen), an error is displayed.

**Exceptions to Losslessness**

1. When compressing with ``--optimize`` or one of the ``--optimize-*`` options, Genozip modifies the file data in ways that are designed to be insignificant for downstream analysis but result in better compression. Since these options are designed to modify the data, losslessness is not preserved.

    |

2. Compressing already-compressed files (BGZF): Many tools that generate genomic files use `BGZF <http://www.htslib.org/doc/bgzip.html#BGZF_FORMAT>`_ to compress them. These typically have a .gz extension (although not all files a .gz extension are compressed with BGZF), and BGZF compression is also used internally in .bam and .bcf files. When ``genozip`` compresses a BGZF-compressed file, it first decompresses it to recover the original underlying data, and then compresses the data with Genozip. Likewise, when ``genounzip`` uncompresses a .genozip file, it recompresses the data back to BGZF if the original file was compressed with BGZF. Genozip records the parameters of the BGZF compression in the .genozip file (such as estimated compression level), and ``genounzip`` attempts to recompress to BGZF format using the same parameters. However, since there are many BGZF compression libraries, each with dozens of parameter combinations, it is possible that ``genounzip`` BGZF compression will achieve a slightly different compression level than the BGZF compression of the original file, resulting in the final BGZF-compressed file differing than the original BGZF file. However, the underlying data is still exactly identical and verified with the *digest*. 

    |

3. Compressing already-compressed files (.bz2 .xz and non-BGZF .gz): ``genozip`` is capable of compressing files that are already compressed with these methods. It does so by first uncompressing the file, and then recompressing with Genozip. However, in these cases, ``genounzip`` does not recompress the files back to .bz2 .xz or .gz.
   