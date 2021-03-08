<!DOCTYPE html>
<!--                                                                                                      -->
<!-- README.md                                                                                            -->
<!-- Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>                                                -->
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
<b>Genozip</b> is a compressor for genomic files - while it can compress any file (i.e. not only genomic files), it is optimized to compress FASTQ, SAM/BAM/CRAM, VCF/BCF, FASTA, GVF, PHYLIP, Chain and 23andMe files.<br>
<br>
<b>Citing</b> Do you find Genozip useful? Please cite:<br>
Lan, D., et al. (2021) <b>Genozip: a universal extensible genomic data compressor</b>. <i>Bioinformatics</i>, https://doi.org/10.1093/bioinformatics/btab102<br>
Lan, D., et al. (2020) <b>genozip: a fast and efficient compression tool for VCF files</b> <i>Bioinformatics</i>, 36, 4091–4092, https://doi.org/10.1093/bioinformatics/btaa290<br>
</p>
<br> 
Typically, a <b>2X-5X improvement over the existing compression</b> is achieved when compressing already-compressed files like .fastq.gz .bam vcf.gz and much higher ratios in some other cases.<br> 
<br> 
<b>Yes</b>, Genozip can compress already-compressed files (.gz .bz2 .xz .bam .cram).<br> 
<br> 
The compression is <b>lossless</b> - the decompressed file is 100% identical to the original file (see documentation for exceptions).<br> 
<b>Sign up</b> to receive low-frequency updates related to Genozip: https://tinyurl.com/genozip<br>
<br>
The command line options are similar to gzip and samtools/bcftools, so if you are familiar with these, it works pretty much the same. To get started, see: http://genozip.com<br>
<br>
<b>Genozip</b> is free for non-commercial use. For a commercial license, please contact <b>sales@genozip.com</b> <br>
<br>
Usage is subject to terms and conditions. The non-commercial license can be viewed on http://genozip.com/license.html<br>
<br>
THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.<br>
