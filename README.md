# genozip
<b>genozip</b> is a compressor for VCF genomic files (it compresses .vcf or .vcf.gz or .vcf.bz2 files). 

It achieves x2 to x5 better compression ratios than gzip because it leverages some properties of the genomic data, such as linkage disequilibrium, to compress better. It is also a lot faster than gzip. 

The compression is lossless - the decompressed VCF file is 100% identical to the original VCF file.

The command line options are very much similar to gzip, so if you're familiar with these, it works pretty much the same. To get started, try: <b>genozip</b> --help

<b><i>Commands</b></i>: \
<b>genozip</b>   - compress one or more files \
<b>genounzip</b> - decompress one or more files \
<b>genols</b>    - show metadata of files or the entire directory \
<b>genocat</b>   - view one or more files

<b><i>Some advanced options</b></i>:
 
<b><i>Concatenating & splitting</b></i>: \
<b>genozip</b> <i>file1.vcf file2.vcf</i> -o <i>concat.vcf.genozip</i> \
<b>genounzip</b> <i>concat.vcf.genozip</i> -O 

<b><i>Calculate MD5 of the VCF file during compression</b></i>: \
<b>genozip</b> <i>file.vcf</i> --md5 \
<b>genols</b> <i>file.vcf.genozip</i> --md5 \
Note: When compressing with --md5, the md5 is automatically verified during <b>genounzip</b>

<b><i>Encryption</b></i>: \
<b>genozip</b> <i>file.vcf</i> --password <i>abc</i>

Feature requests and bug reports: <b>bugs@genozip.com</b>

<b>genozip</b> is free for non-commercial use. For a commercial license, please contact <b>sales@genozip.com</b>

Usage is subject to terms and conditions, see the LICENSE.commercial.txt and LICENSE.non-commercial.txt file.

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
