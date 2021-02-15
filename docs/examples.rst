Examples
========

**Simple compression and uncompression**

    | ``genozip sample.bam``
    |
    | ``genounzip sample.bam.genozip``

**Creating a reference file**
    
    ``genozip --make-reference myfasta.fa``

**Compressing a FASTQ, SAM/BAM or VCF file(s) with a reference**

    | ``genozip --reference myfasta.ref.genozip mysample1.fq mysample2.fq mysample3.fq``
    |
    | ``genozip --reference myfasta.ref.genozip mysample.bam``
    |
    | ``genozip --reference myfasta.ref.genozip mysamples.vcf.gz``
    |
    | ``genozip --reference myfasta.ref.genozip *``
    | compresses all files in the current directory

    | *Notes*:
    |
    | 1. Genozip can compress with or without a reference - using a reference achieves much better compression when compressing FASTQ or unaligned SAM/BAM, and modestly better compression in other cases.
    |
    | 2. SAM/BAM - compression of aligned or unaligned SAM/BAM files is possible. Sorting makes no difference.
    |
    | 3. Long reads - compression of long reads (Pac Bio / Nanopore) achieves signficantly better results when compressing an aligned BAM vs an unaligned BAM or FASTQ.
    |
    | 4. Compression of CRAM (but not SAM or BAM) files requires samtools to be installed.
    |
    | 5. Use ``--REFERENCE`` instead of ``--reference`` to store the relevant parts of the reference file as part of the compressed file itself, which will then allow decompression with genounzip without need of the reference file.

**Compressing and uncompressing paired-end reads with --pair** 

    | ``genozip --reference myfasta.ref.genozip --pair mysample-R1.fastq.gz mysample-R2.fastq.gz``
    |
    | ``genounzip --reference myfasta.ref.genozip --unbind mysample-R1+2.fastq.genozip``

    | *Note*: with ``--pair`` genozip uses similarities between the files to enhance compression.

**Using genozip in a pipline**

    | ``genocat mysample.sam.genozip | samtools - .....``
    |
    | ``my-sam-outputing-method | genozip - --input sam --output mysample.sam.genozip``

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

    | *Notes*:
    |
    | 1. ``--regions`` works with VCF, SAM/BAM, FASTA, 23andMe, GVF and reference files ; ``--grep`` works with FASTQ, FASTA ; ``--samples`` works with VCF ; ``--downsample`` works with all types
    |
    | 2. There is no need for a separate indexing step or index file
    |
    | 3. Many more options (see :doc:`genocat` for full list): ``--no-header`` ; ``--header-only`` ; ``--header-one`` ; ``--sequential`` ; ``--list-chroms`` ; ``--drop-genotypes`` ; ``--GT-only``

**Binding mutiple files into a single genozip file and unbinding**

    | ``genozip *.fq.gz -o all-samples.fq.genozip``  
    | Binds all .fq.gz files in the current directory.
    |
    | ``genounzip my-project.fq.genozip --unbind``

**Compressing even better, with some minor modifications of the data**

    | ``genozip file.bam --optimize``
    | *Note*: compression with ``--optimize`` is *not* lossless - see :doc:`genozip` for details.

**Compressing faster, sacrificing a bit of compression ratio**

    | ``genozip file.bam --fast``

**Encrypting (256 bit AES)**

    | ``genozip file.vcf --password abc``
    | ``genounzip file.vcf.genozip --password abc``

**Converting SAM/BAM to FASTQ**

    | ``genounzip file.bam.genozip --fastq``

**Converting 23andMe to VCF**

    | ``genounzip genome_mydata-Full.txt.genozip --vcf -e GRCh37.ref.genozip``

**Generating a samtools/bcftools index file when uncompressing**
    | ``genounzip file.bam.genozip --index``

**Calculating the MD5 of the underlying textual file (also included in --test)**

    | ``genozip file.vcf --md5``
    | ``genounzip file.vcf.genozip --md5``
    | ``genols file.vcf.genozip``

**Compressing and then verifying that the compressed file decompresses correctly**

    | ``genozip file.vcf --test``
