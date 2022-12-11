<!DOCTYPE html>
<!--                                                                                                      -->
<!-- README.md                                                                                            -->
<!-- Copyright (C) 2019-2022 Genozip Limited. Patent Pending.                                                -->
<!-- Please see terms and conditions in the file LICENSE.txt
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
<b>Genozip</b> is a compressor for BAM / SAM / CRAM, FASTQ, VCF and other genomic files - see https://genozip.com<br>
<br>
Genozip is also available on <b>Conda</b> and binary downloads, see <a href="https://genozip.com/installing">installation options</a>.<br>
<br>
For Illumina data <i>.bam</i> and <i>.fastq.gz</i> files, the typical gain over gzip is around 4X. For PacBio and Oxford Nanopore data aligned <i>.bam</i> files, the gain is typically around 2X. For <i>.vcf.gz</i> files, the gain over gzip is typically 3-6X. Here are some <a href="https://genozip.com/benchmarks">benchmarks</a>.<br>
<br>
<b>Yes</b>, Genozip can compress already-compressed files (.gz .bz2 .xz .bam .cram).<br> 
<br> 
The compression is <a href="https://genozip.com/losslessness">lossless</a> - the decompressed file is 100% identical to the original file (see documentation for exceptions), verifiable by MD5.<br>
<br>
<b>Genozip</b> is free for certain academic research use. See <a href="https://genozip.com/get-genozip">eligibility and other licensing options</a> or contact <a href="mailto://sales@genozip.com">sales@genozip.com</a> <br>
<br>
Usage is subject to terms and conditions. See our <a href="https://genozip.com/license">license</a> for details.<br>
<br>
<b>THIS SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS, COPYRIGHT HOLDERS OR DISTRIBUTORS OF THIS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</b><br>
