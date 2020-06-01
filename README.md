<!DOCTYPE html>
<!--                                                                                                    -->
<!-- README.md                                                                                          -->
<!-- Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>                                              -->
<!-- Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt -->
<!--                                                                                                    -->
<!-- This file needs to be compliant to both Markdown and HTML. It is:                                  -->
<!-- 1. rendered as README.md by github                                                                 -->
<!-- 2. copied as HTML to the Mac installer                                                             -->
<!-- 3. copied into meta.yaml, after removing all the HTML stuff                                        -->
<!--                                                                                                    -->
<h1>genozip</h1><br> 
<br>
(also available on <b>conda</b> and <b>Docker Hub</b>)
<br>
<b>genozip</b> is a compressor for genomic files - it compresses VCF/BCF, SAM/BAM, fastq, fasta, GVF and 23andMe files. If can even compress them if they are already compressed with .gz .bz2 .xz (for full list of supported file types see 'genozip --input --help').<br>
<br>
It achieves x2 to x5 better compression ratios than gzip because it leverages some properties specific to genomic data to compress better. It is also a lot faster than gzip.<br>
<br>
The compression is lossless - the decompressed file is 100% identical to the original file.<br>
<br>
The command line options are similar to gzip and bcftools, so if you're familiar with these, it works pretty much the same. To get started, try: <b>genozip</b> --help<br>
<br>
<b><i>Commands</i></b>: <br>
<b>genozip</b>   - compress one or more files <br>
<b>genounzip</b> - decompress one or more files <br>
<b>genols</b>    - show metadata of one or more files or the entire directory <br>
<b>genocat</b>   - view one or more files <br>
<br>
<b><i>Some advanced options:</i></b><br>
<br>
<b><i>Lookups:</i></b><br>
<b>genocat</b> -r ^Y,MT <i>file1.vcf</i>       -- displays all chromosomes except Y and MT<br>
<b>genocat</b> -r -10000 <i>file1.vcf</i>      -- displays positions up to 10000<br>
<b>genocat</b> -s SMPL1,SMPL2 <i>file1.vcf</i> -- displays 2 samples<br>
Note: there is no need for a separate indexing step or index file<br>
<br>
<b><i>Concatenating & splitting:</i></b><br>
<b>genozip</b> <i>file1.vcf file2.vcf</i> -o <i>concat.vcf.genozip</i> <br>
<b>genounzip</b> <i>concat.vcf.genozip</i> -O <br>
<br>
<b><i>Calculating the MD5:</i></b><br>
<b>genozip</b> <i>file.vcf</i> --md5 <br>
<br>
<b><i>Encryption:</i></b><br>
<b>genozip</b> <i>file.vcf</i> --password <i>abc</i> <br>
<br>
<b><i>Even better compression, with some minor modifications of the data:</i></b><br>
<b>genozip</b> <i>file.vcf</i> --optimize <br>
<br>
<b><i>Compress and then verify that the compressed file decompresses correctly:</i></b><br>
<b>genozip</b> <i>file.vcf</i> --test <br>
<br>
Do you find genozip to be helpful in your research? Please be so kind as to support continued development by citing
<b>Citing:</b> https://doi.org/10.1093/bioinformatics/btaa290
<br>
Feature requests and bug reports: <b>bugs@genozip.com</b> <br>
<br>
<b>genozip</b> is free for non-commercial use. For a commercial license, please contact <b>sales@genozip.com</b> <br>
<br>
Usage is subject to terms and conditions. The non-commercial license can be viewed with <b>genozip</b> --license<br>
<br>
THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.<br>
