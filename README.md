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
<h1>genozip</h1><br>
<br>
(also available on <b>Conda</b> and <b>Docker Hub</b>)<br>
<br>
<b>genozip</b> is a compressor for genomic files - it compresses FASTQ, SAM/BAM/CRAM, VCF/BCF, FASTA, GVF and 23andMe files. If can even compress them if they are already compressed with .gz .bz2 .xz (for full list of supported file types see 'genozip --input --help').<br>
<br>
It achieves x2 to x5 better compression ratios than gzip because it leverages some properties specific to genomic data to compress better. It is also a lot faster than gzip.<br>
<br>
The compression is lossless - the decompressed file is 100% identical to the original file.<br>
Notes: <br>
1. Losslessness is relative to the underlying textual file - for example, when compressing a .bam or a .fq.gz file - the losslessness is relative to the underlying SAM or FASTQ file respectively<br>
2. The one exception is when using the --optimize option which is lossy (see --help for details)<br>
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
<b>genozip</b> --reference <i>myfasta.ref.genozip</i> *    ← compresses all files in the current directory<br>
<br>
Notes:<br>
1. genozip can compress with or without a reference - using a reference often achieves much better compression<br>
2. SAM/BAM - compression of aligned or unaligned SAM/BAM files is possible. Sorting makes no difference<br>
3. Long reads - compression of long reads (Pac Bio / Nanopore) achieves signficantly better results when compressing an aligned BAM vs an unaligned BAM or FASTQ<br>
4. Compression of BAM and CRAM (but not SAM) files requires samtools to be installed<br>
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
<b><i>Lookups:</i></b><br>
<b>genocat</b> --regions ^Y,MT <i>mysamples.vcf.genozip</i>  ← displays all chromosomes except Y and MT<br>
<b>genocat</b> --regions -10000 <i>mysample.sam.genozip</i>  ← displays positions up to 10000<br>
<b>genocat</b> --samples SMPL1,SMPL2 <i>mysamples.vcf.genozip</i>  ← displays 2 samples<br>
<b>genocat</b> --grep August-10 <i>myfasta.fa.genozip</i>  ← displays contigs/reads that have "August-10" in the header<br>
Notes:<br>
1. --regions works with VCF, SAM/BAM, FASTA ; --grep works with FASTQ, FASTA ; --samples works with VCF<br>
2. There is no need for a separate indexing step or index file<br>
3. Many more options (see --help for full list): --no-header ; --header-only ; --sequential ; --list-chroms ; --drop-genotypes ; --GT-only<br>
<br>
<b><i>Binding mutiple files into a single genozip file & unbinding:</i></b><br>
<b>genozip</b> <i>*.fq.gz</i> -o <i>all-samples.fq.genozip</i> ← binds all .fq.gz files in the current directory<br>
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
<b><i>Calculating the MD5 of the underlying textual file (also included in </i>--test<i>):</i></b><br>
<b>genozip</b> <i>file.vcf</i> --md5 <br>
<b>genounzip</b> <i>file.vcf.genozip</i> --md5 <br>
<br>
<b><i>Compressing and then verifying that the compressed file decompresses correctly:</i></b><br>
<b>genozip</b> <i>file.vcf</i> --test <br>
<br>
<b><i>Citing</i></b><br>
Do you find genozip useful? Please support continued development by citing:<br>
Lan, D., et al. (2020) <i>Bioinformatics</i>, 36, 4091–4092, https://doi.org/10.1093/bioinformatics/btaa290<br>
<br>
Feature requests and bug reports: <b>bugs@genozip.com</b> <br>
<br>
<b>genozip</b> is free for non-commercial use. For a commercial license, please contact <b>sales@genozip.com</b> <br>
<br>
Usage is subject to terms and conditions. The non-commercial license can be viewed with <b>genozip</b> --license<br>
<br>
THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.<br>
