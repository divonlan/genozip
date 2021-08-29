.. _compression:

Compression
===========

Genozip can run on any type of file, but it is optimized to compress genomic file formats.

**Simple compression and decompression**

    | ``genozip sample.bam``
    |
    | ``genounzip sample.bam.genozip``

    
**Compressing and decompressing FASTQ with paired-end reads with a reference** 

    | First, create a refrence file. This is a one-time step, that may take up to 10 minutes: 
    | ``genozip --make-reference myfasta.fa``
    |   
    | Second, compress or decompresses your file(s) using the reference:
    | ``genozip --reference myfasta.ref.genozip --pair mysample-R1.fastq.gz mysample-R2.fastq.gz``
    |
    | ``genounzip --reference myfasta.ref.genozip mysample-R1+2.fastq.genozip``
    |
    | :ref:`Details <fastq>`
    
**Compressing and decompressing multiple files into a tar file (preserving directory structure)**

    | ``genozip *.bam --tar mysamples.tar``  
    |
    | ``tar xvf mydata.tar | genounzip --files-from - --replace``
    |
    | :ref:`Details <archiving>`

**Using genozip in a pipline**

    | ``genocat mysample.sam.genozip | samtools - .....``
    |
    | ``my-sam-outputing-method | genozip - --input sam --output mysample.sam.genozip``

**Viewing compression stats**

    | ``genocat sample.bam.genozip --stats``
     
**Lookups, downsampling and other subsets**

    | ``genocat --regions chr1:10000-20000 mysamples.vcf.genozip``  
    | Displays a specific region.
    |
    | ``genocat --regions ^Y,MT mysample.bam.genozip``
    | Displays all alignments except Y and MT contigs.
    |
    | ``genocat --regions chrM GRCh38.fa.genozip``  
    | Dislays the sequence of chrM.
    |
    | ``genocat --samples SMPL1,SMPL2 mysamples.vcf.genozip``   
    | Displays 2 samples.
    |
    | ``genocat --grep 1101:2392 myreads.fq.genozip``   
    | Displays reads that have "1101:2392" anywhere in the description.
    |
    | ``genocat --downsample 10 mysample.fq.genozip``   
    | Displays 1 in 10 reads.

    | *Note*: These are just some examples - there are many more subsetting options see :doc:`genocat`.

**Compressing even better, with some minor modifications of the data**

    | ``genozip file.bam --optimize``
    | *Note*: compression with ``--optimize`` is *not* lossless - see :doc:`genozip` for details.

**Compressing faster, sacrificing a bit of compression ratio**

    | ``genozip file.bam --fast``

**Encrypting (256 bit AES)**

    | ``genozip file.vcf --password abc``
    | ``genounzip file.vcf.genozip --password abc``
    |
    | :ref:`Details <encryption>`

**Converting SAM/BAM to FASTQ**

    | ``genounzip file.bam.genozip --fastq``
    |
    | :ref:`Details <bam>`

**Generating a samtools/bcftools index file when uncompressing**
    | ``genounzip file.bam.genozip --index``

**Calculating the MD5 of the underlying textual file (also included in --test)**

    | ``genozip file.vcf --md5``
    | ``genounzip file.vcf.genozip --md5``
    | ``genols file.vcf.genozip``

**Compressing and then verifying that the compressed file decompresses correctly**

    | ``genozip file.vcf --test``
