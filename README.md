# genozip
_genozip_ is a compressor for VCF genomic files (it compresses .vcf or .vcf.gz or .vcf.bz2 files). 

It achieves x2 to x5 better compression ratios than gzip because it leverages some properties of the genomic data, such as linkage disequilibrium, to compress better. It is also a lot faster than gzip. 

It is very easy to use - in fact, if you're familiar with gzip, it works pretty much the same.
    
The compression is lossless - the decompressed VCF file is 100% identical to the original VCF file.

The command line options are very much similar to gzip, so if you're familiar with these, it works pretty much the same. Try _genozip --help_ to get started.

__Commands__:
_genozip_   - compress one or more files
_genounzip_ - decompress one ore more files
_genols_    - show metadata of files or the entire directory
_genocat_   - view one or more files

Some advanced options:

__Concatenating & splitting__: 
_genozip_ file1.vcf file2.vcf _-o_ concat.vcf.genozip
_genounzip_ concat.vcf.genozip _-O_ 

__MD5 hash__:
_genozip --md5_ file.vcf 
_genols --md5_ file.vcf.genozip
Note: The md5 hash is automatically verified during _genounzip_

__Encryption__: 
_genozip --password_ abc file.vcf

Feedback and bug reports: _bugs@genozip.com_
Genozip is free for non-commercial use. For a commercial license, please contact _sales@genozip.com_

Usage is subject to terms and conditions, see the LICENSE.commercial.txt and LICENSE.non-commercial.txt file.

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
