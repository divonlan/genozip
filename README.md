# genozip
__genozip__ is a compressor for VCF genomic files (it compresses .vcf or .vcf.gz or .vcf.bz2 files). 

It achieves x2 to x5 better compression ratios than gzip because it leverages some properties of the genomic data, such as linkage disequilibrium, to compress better. It is also a lot faster than gzip. 

The compression is lossless - the decompressed VCF file is 100% identical to the original VCF file.

The command line options are very much similar to gzip, so if you're familiar with these, it works pretty much the same. To get started, try: __genozip__ --help

***Commands***: \
__genozip__   - compress one or more files \
__genounzip__ - decompress one or more files \
__genols__    - show metadata of files or the entire directory \
__genocat__   - view one or more files

***Some advanced options***:
 
***Concatenating & splitting***: \
__genozip__ _file1.vcf file2.vcf_ -o _concat.vcf.genozip_ \
__genounzip__ _concat.vcf.genozip_ -O 

***MD5 of the VCF file***: \
__genozip__ --md5 _file.vcf_ \
__genols__ --md5 _file.vcf.genozip_ \
Note: When compressing with --md5, the md5 is automatically verified during _genounzip_

***Encryption***: \
__genozip__ _file.vcf_ --password _abc_

Feature requests and bug reports: __bugs@genozip.com__

__genozip__ is free for non-commercial use. For a commercial license, please contact __sales@genozip.com__

Usage is subject to terms and conditions, see the LICENSE.commercial.txt and LICENSE.non-commercial.txt file.

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
