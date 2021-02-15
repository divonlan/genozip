<!DOCTYPE html>
<!--                                                                                                      -->
<!-- README.md                                                                                            -->
<!-- Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>                                                -->
<!-- Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt   -->
<!--                                                                                                      -->
<!-- This file needs to be compliant to both Markdown and HTML. It is:                                    -->
<!-- 1. rendered as README.md by github                                                                   -->
<!-- 2. copied as HTML to the Mac installer                                                               -->
<!-- 3. copied into meta.yaml, after removing all the HTML stuff                                          -->
<!-- 4. rendered as README.md in Docker Hub                                                               -->
<!-- 5. converted to Markdown and embedded in conda/README.template.md to generate conda feedstock README -->
<!--                                                                                                      -->
<!-- To preview in Visual Studio Code: Ctrl+Shift+V with the "HTML Preview" extension                     -->
<h1>Genozip</h1><br>
<br>
(available on <b>Conda</b>, <b>Docker Hub</b> and https://github.com/divonlan/genozip ; Documentation: http://genozip.com)<br>
<br>
<b>Genozip</b> is a compressor for genomic files - while it can compress any file (i.e. not only genomic files), it is optimized to compress FASTQ, SAM/BAM/CRAM, VCF/BCF, FASTA, GVF, Phylip and 23andMe files. It can even compress them if they are already compressed with .gz .bz2 .xz (for full list of supported file types see 'genozip --help=input').<br>
<br>
The compression ratio depends on the data being compressed, and you can usually expect about a 1.5-3X ratio when compressing .bam, 2X-5X for .fastq.gz files (i.e. compressing already-compressed files), and up to 200X when compressing an uncompressed high-sample-count .vcf file with only GT data.<br>
<br>
<b>Sign up</b> to receive low-frequency updates related to Genozip: https://tinyurl.com/genozip<br>
<br>
The compression is lossless - the decompressed file is 100% identical to the original file.<br>
Notes: <br>
1. Losslessness when the original file is already compressed (with .gz, .bz2 or .xz as well as .bam) is relative to the underlying uncompressed file. If the original file was compressed with BGZF (as in BAM and most genomic files with a .gz extension) - genounzip will re-compresses the file with BGZF upon decompression (unless --plain is specified), and, further, it will attempt to recover the exact same BGZF compression as in the original file. However, sometimes exact-same BGZF compression is not possible due to different libraries used.<br>
2. The one exception to Genozip's strict losslessness is when using the --optimize option which is lossy (see --help for details)<br>
<br>
The command line options are similar to gzip and samtools/bcftools, so if you are familiar with these, it works pretty much the same. To get started, try: <b>genozip</b> --help<br>
<br>
<b><i>Commands</i></b>: <br>
<b>genozip</b>   - compress one or more files <br>
<b>genounzip</b> - decompress one or more files <br>
<b>genols</b>    - show metadata of one or more files or the entire directory <br>
<b>genocat</b>   - view one or more files <br>
<br>
<b><i>Some examples:</i></b><br>
<br>
<b><i>Creating a reference file:</i></b><br>
<b>genozip</b> --make-reference <i>myfasta.fa</i><br>
<br>
<b><i>Compressing a FASTQ, SAM/BAM or VCF file(s) with a reference:</i></b><br>
<b>genozip</b> --reference <i>myfasta.ref.genozip</i> mysample1.fq mysample2.fq mysample3.fq<br>
<b>genozip</b> --reference <i>myfasta.ref.genozip</i> mysample.bam<br>
<b>genozip</b> --reference <i>myfasta.ref.genozip</i> mysamples.vcf.gz<br>
<b>genozip</b> --reference <i>myfasta.ref.genozip</i> *  &nbsp&nbsp&nbsp←compresses all files in the current directory<br>
<br>
Notes:<br>
1. Genozip can compress with or without a reference - using a reference achieves much better compression when compressing FASTQ or unaligned SAM/BAM, and modestly better compression in other cases<br>
2. SAM/BAM - compression of aligned or unaligned SAM/BAM files is possible. Sorting makes no difference<br>
3. Long reads - compression of long reads (Pac Bio / Nanopore) achieves signficantly better results when compressing an aligned BAM vs an unaligned BAM or FASTQ<br>
4. Compression of CRAM (but not SAM or BAM) files requires samtools to be installed<br>
5. Use --REFERENCE instead of --reference to store the relevant parts of the reference file as part of the compressed file itself, which will then allow decompression with genounzip without need of the reference file.<br>
<br>
<b><i>Compressing and uncompressing paired-end reads with --pair - better than compressing FASTQs individually</i></b><br>
<b>genozip</b> --reference <i>myfasta.ref.genozip</i> --pair mysample-R1.fastq.gz mysample-R2.fastq.gz<br>
<b>genounzip</b> --reference <i>myfasta.ref.genozip</i> --unbind mysample-R1+2.fastq.genozip<br>
<br>
<b><i>Using</i> genozip <i>in a pipline:</i></b><br>
<b>genocat</b> <i>mysample.sam.genozip</i> | samtools - .....<br>
my-sam-outputing-method | <b>genozip</b> - --input sam --output <i>mysample.sam.genozip</i><br>
<br>
<b><i>Lookups, downsampling and other subsets:</i></b><br>
<b>genocat</b> --regions chr1:10000-20000 <i>mysamples.vcf.genozip</i>&nbsp&nbsp&nbsp←displays a specific region<br>
<b>genocat</b> --regions ^Y,MT <i>mysample.bam.genozip</i>&nbsp&nbsp&nbsp←displays all alignments except Y and MT contigs<br>
<b>genocat</b> --regions chrM <i>GRCh38.fa.genozip</i>&nbsp&nbsp&nbsp←displays the sequence of chrM<br>
<b>genocat</b> --samples SMPL1,SMPL2 <i>mysamples.vcf.genozip</i>&nbsp&nbsp&nbsp←displays 2 samples<br>
<b>genocat</b> --grep 1101:2392 <i>myreads.fq.genozip</i>&nbsp&nbsp&nbsp←displays reads that have "1101:2392" anywhere in the description<br>
<b>genocat</b> --downsample 10 <i>mysample.fq.genozip</i>&nbsp&nbsp&nbsp←displays 1 in 10 reads<br>
Notes:<br>
1. --regions works with VCF, SAM/BAM, FASTA, 23andMe, GVF and reference files ; --grep works with FASTQ, FASTA ; --samples works with VCF ; --downsample works with all types<br>
2. There is no need for a separate indexing step or index file<br>
3. Many more options (see --help for full list): --no-header ; --header-only ; --header-one ; --sequential ; --list-chroms ; --drop-genotypes ; --GT-only<br>
<br>
<b><i>Binding mutiple files into a single genozip file and unbinding:</i></b><br>
<b>genozip</b> <i>*.fq.gz</i> -o <i>all-samples.fq.genozip</i>&nbsp&nbsp&nbsp←binds all .fq.gz files in the current directory<br>
<b>genounzip</b> <i>my-project.fq.genozip</i> --unbind <br>
<br>
<b><i>Compressing even better, with some minor modifications of the data (therefore not lossless, see --help for details):</i></b><br>
<b>genozip</b> <i>file.bam</i> --optimize <br>
<br>
<b><i>Compressing faster, sacrificing a bit of compression ratio:</i></b><br>
<b>genozip</b> <i>file.bam</i> --fast <br>
<br>
<b><i>Encrypting (256 bit AES):</i></b><br>
<b>genozip</b> <i>file.vcf</i> --password <i>abc</i> <br>
<b>genounzip</b> <i>file.vcf.genozip</i> --password <i>abc</i> <br>
<br>
<b><i>Converting SAM/BAM to FASTQ:</i></b><br>
<b>genounzip</b> <i>file.bam.genozip</i> --fastq<br>
<br>
<b><i>Converting 23andMe to VCF:</i></b><br>
<b>genounzip</b> <i>genome_mydata-Full.txt.genozip</i> --vcf -e GRCh37.ref.genozip<br>
<br>
<b><i>Generating a samtools/bcftools index file when uncompressing:</i></b><br>
<b>genounzip</b> <i>file.bam.genozip</i> --index<br>
<br>
<b><i>Calculating the MD5 of the underlying textual file (also included in </i>--test<i>):</i></b><br>
<b>genozip</b> <i>file.vcf</i> --md5 <br>
<b>genounzip</b> <i>file.vcf.genozip</i> --md5 <br>
<b>genols</b> <i>file.vcf.genozip</i><br>
<br>
<b><i>Compressing and then verifying that the compressed file decompresses correctly:</i></b><br>
<b>genozip</b> <i>file.vcf</i> --test <br>
<br>
<b><i>Citing</i></b><br>
Do you find Genozip useful? Please support continued development by citing:<br>
Lan, D., et al. <i>Bioinformatics</i>, 36, 4091–4092, July 2020, https://doi.org/10.1093/bioinformatics/btaa290<br>
<br> 
Feature requests and bug reports: <b>bugs@genozip.com</b> <br>
<br>
<b>Genozip</b> is free for non-commercial use. For a commercial license, please contact <b>sales@genozip.com</b> <br>
<br>
Usage is subject to terms and conditions. The non-commercial license can be viewed with <b>genozip</b> --license<br>
<br>
THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.<br>
