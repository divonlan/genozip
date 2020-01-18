
Compress or uncompress VCF (Variant Call Format) files

Usage: genozip [options]... [files]...
       genounzip [options]... [files]...
       genocat [options]... [files]...
       genols files

Actions - use at most one of these actions:
   -z --compress     Compress a .vcf, .vcf.gz or .vcf.bz2 file (Yes! You can compress an already-compressed file).
                     The source file is left unchanged. This is the default action of genozip
   -d --decompress   Decompress a .vcf.genozip file. The .vcf .genozip file is left unchanged.
                     This is the default action of genounzip
   -l --list         List the compression ratios of the .vcf.genozip files
                     This is the default action of genols
   -t --test         Test genozip. Compress the .vcf file(s), uncompress, and then compare the
                     result to the original .vcf - all in memory without writing to any file
   -h --help         Show this help page
   -L --license      Show the license terms and conditions for this product
   -V --version      Display version number

Flags:
   -c --stdout       Send output to standard output instead of a file. -dc is the default action of genocat
   -f --force        Force overwrite of the output file, or force writing .vcf.genozip data to standard output
   -R --replace      Replace the source file with the result file, rather than leaving it unchanged
   -o --output       Output file name. This option can also be used to concatenate multiple input files
                     with the same individuals, into a single concatenated output file
   -p --password     Password-protected - encrypted with 256-bit AES
   -m --md5          When compressing - records the MD5 hash of the VCF file in the genozip file header
                     When listing (--list) - shows the MD5 of each file
                     Decompress always compares the MD5 to the uncompressed VCF, if compress was done with --md5
   -q --quiet        Don't show the progress indicator
   -@ --threads      Specify the number of threads to use. By default, genozip uses all cores available to it
   --show-content    Show the information content of VCF files and the compression ratios of each component
   --show-alleles    Output allele values to stdout. Each row corresponds to a row in the VCF file
                     Mixed-ploidy regions are padded, and 2-digit allele values are replaced by an ascii character
   --show-time       Show what functions are consuming the most time (useful mostly for developers of genozip)
   --show-memory     Show what buffers are consuming the most memory (useful mostly for developers of genozip)

One or more file names may be given, or if omitted, standard input/output is used instead

Genozip is available for free for non-commercial use. Commercial use requires a commercial license.

For bug reports: bugs@genozip.com and license inquiries: sales@genozip.com

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

