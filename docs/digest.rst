.. _digest:

Verifying file integrity
========================

When compressing with ``genozip``, as the source file is being read from disk (or standard input), Genozip calculates a numeric value which is a function of the entire source file's content - a value known as the *digest* of the file. If the source file is compressed in .gz .bz2 or .xz - the digest is calculated on the source data after decompressing it.

The digest value is then stored in the compressed genozip file. 

When the genozip file is decompressed ``genounzip`` or ``genocat``, a digest is again calculated on the final output data. The output data digest is then compared to the digest stored in the genozip file. If these two digests are identical, then we know with a very high probability (see below), that the uncompressed file is exactly identical to the source file. If it is not, ``genounzip`` or ``genocat`` will report an error. Short of intentional tampering of the file, or a bug in Genozip, this should never happen.

When running ``genozip`` with the ``--test`` option (which is always a good idea), ``genounzip --test`` automatically is run on the resulting output genozip file, to ensure its integrity. ``genounzip --test``, which may also be run separately, decompresses the file in memory just in order to compare the digest - the resulting reconstructed data is then discarded and not written to disk.

**Probabilites**

By default, the digest is calculated using the `Alder32 algorithm <https://en.wikipedia.org/wiki/Adler-32>`_, resulting in a 32 bit digest. This would mean that the probability of a corrupt file incorrectly passing the integrity test being about 1 in 4 billion. 

Genozip also offers a stronger verification with the `MD5 algorithm <https://en.wikipedia.org/wiki/MD5>`_, which may be invoked by compressing with ``genozip --md5``. MD5 results in a 128 bit digest, bringing down the probability of a corrupt file incorrectly passing the integrity test to about 2^(-128) or approximately 1 in 34,000,000,000,000,000,000,000,000,000,000,000,000. 

MD5, unlike Adler32, also has cryptographic usage - an attacker interested in maliciously replacing a file with another file, will find it practically impossible to purposely generate the other file so that it has the same MD5 digest as the original file. Therefore, recording the MD5 digest of the original file, and comparing it to the MD5 of the final file, guarantees, for all practical purposes, that this is indeed the same file.  

**Cases in which genozip doesn't calculate a digest**

There are some cases in which you may request ``genozip`` to change the source data before compressing it. In these cases, the digest is not calculated. These cases are:

- Using ``--optimize`` or any of the ``--optimize-*`` options 
- Generating a dual-coordinate file with ``--chain``
- Compressing a Luft file (a lifted-over dual-coordinate file)
- Compressing a SAM or BAM file with ``--taxid``

**Cases in which genocat doesn't verify the digest**

When using any ``genocat`` option which results in intentional modification of the output data, the digest is not verified. These options include (partial list): ``--regions`` ``--grep`` ``--downsample`` ``--taxid`` ``--drop-genotypes`` ``--gt-only`` ``--no-header``, ``--header-only``, ``header-one``, ``--MAPQ``, ``--FLAG``, ``--bases``, ``--samples``, ``--sequential``, ``--luft``, ``--show-dvcf``, ``--head``, ``--tail``, ``--lines``, ``--one-vb``, ``--one-component``

**Concatenated files**

When concatenating multiple files into a single genozip file, two types of digests are calculated - 

1) A digest for each component separately

2) A digest for entire concatenated file as it will be viewed using ``genocat``. For files formats that include a header (eg VCF, SAM), this will include the header of the first component only.
   
**Command line examples**

Compressing with Adler32 digest:

::

    genozip myfile.sam

Compressing with MD5 digest:

::

    genozip --md5 myfile.sam

Compressing with MD5 and verifying:

::
    
    genozip --test myfile.sam

Verifying a compressed file by decompressing in memory without output:

::
    
    genounzip --test myfile.sam.genozip

Viewing the MD5 digest of a single file or the entire directory:

::

    genols myfile.sam.genozip
    genols

Creating a concatenated file, uncompressing it back to the individual componenets verifying against the individual components' digests, and outputting the concatenated data - verifying against the concatenated digest:

::

    genozip file1.sam file2.sam file3.sam --output myfiles.sam.genozip
    genounzip --prefix newdir/ myfiles.sam.genozip
    genocat myfiles.sam.genozip --output concatenated.sam

Viewing the flow of creating the digest (useful mostly for Genozip developers):

::
    
    genozip --show-digest myfile.sam
    genounzip --test --show-digest myfile.sam
    genocat --show-digest myfile.sam

Viewing the digest of the entire as it appears in the GENOZIP_HEADER section of the genozip file, and the digest of the individual components as it appears in the TXT_HEADER sections (useful mostly for Genozip developers):

::

    genocat --show-header=GENOZIP_HEADER myfile.sam.genozip
    genocat --show-header=TXT_HEADER myfile.sam.genozip


