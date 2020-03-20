
Compress VCF (Variant Call Format) files

Usage: genozip [options]... [files]...

See also: genounzip genocat genols

Actions - use at most one of these actions:
   -z --compress     This is the default action. Compress a .vcf, .vcf.gz or .vcf.bz2 file (Yes! You can compress an
                     already-compressed file). The source file is left unchanged. This is the default action of genozip

   -d --decompress   Same as running genounzip. For more details, run: genounzip --help

   -l --list         Same as running genols. For more details, run: genols --help

   -h --help         Show this help page. Use with -f to see developer options.

   -L --license      Show the license terms and conditions for this product

   -V --version      Display version number

Flags:
   -c --stdout       Send output to standard output instead of a file

   -f --force        Force overwrite of the output file, or force writing .vcf.genozip data to standard output

   -^ --replace      Replace the source file with the result file, rather than leaving it unchanged

   -o --output       <output-filename>. This option can also be used to concatenate multiple input files with the same
                     individuals, into a single concatenated output file

   -p --password     <password>. Password-protected - encrypted with 256-bit AES

   -m --md5          Calculates the MD5 hash of the VCF file. When the resulting file is decompressed, this MD5 will be
                     compared to the MD5 of the decompressed VCF.
                     Note: for compressed files, e.g. myfile.vcf.gz, the MD5 calculated is that of the original,
                     uncompressed file. 
                     In addition, if the VCF file has Windows-style \r\n line endings, the md5 will be that of the
                     modified file with the \r removed

   -q --quiet        Don't show the progress indicator or warnings

   -Q --noisy        The --quiet is turned on by default when outputting to the terminal. --noisy stops the suppression
                     of warnings

   -t --test         After compressing normally, decompresss in memory (i.e. without writing the decompressed file to
                     disk) - comparing the MD5 of the resulting decompressed file to that of the original VCF. This
                     option also activates --md5.

   -@ --threads      <number>. Specify the maximum number of threads. By default, this is set to the number of cores
                     available. The number of threads actually used may be less, if sufficient to balance CPU and I/O

   --show-content    Show the information content of VCF files and the compression ratios of each component

Optimizing:
   -9 --optimize     Modify the VCF file in ways that are likely insignificant for analytical purposes, but make a
                     signficant difference for compression. At the moment, these optimizations include:
                     - PL data: Phred values of over 60 are changed to 60.     Example: '0,18,270' -> '0,18,60'
                     - GL data: Numbers are rounded to 2 significant digits.   Example: '-2.61618,-0.447624,-0.193264'
                     -> '-2.6,-0.45,-0.19'
                     - GP data: Numbers are rounded to 2 significant digits, as with GL.
                     - VQSLOD data: Number is rounded to 2 significant digits. Example: '-4.19494' -> '-4.2'
                     Note: due to these data modifications, files compressed with --optimized are NOT identical as the
                     original VCF after decompression. For this reason, it is not possible to use this option in
                     combination with --test or --md5

   -B --vblock       <number between 1 and 2048>. Sets the maximum size of memory (in megabytes) of VCF file data that
                     can go into one variant block. By default, this is set to 128 MB. The variant block is the basic
                     unit of data on which genozip and genounzip operate. This value affects a number of things: 1.
                     Memory consumption of both compression and decompression are linear with the variant block size.
                     2. Compression is sometimes better with larger block sizes, in particular if the number of samples
                     is small. 3. Smaller blocks will result in faster 'genocat --regions' lookups

   -S --sblock       <number>. Sets the number of samples per sample block. By default, it is set to 4096. When
                     compressing or decompressing a variant block, the samples within the block are divided to sample
                     blocks which are compressed separately. A higher value will result in a better compression ratio,
                     while a lower value will result in faster 'genocat --samples' lookups

   --gtshark         Use gtshark instead of the default bzlib as the final compression step for allele data (the GT
                     subfield in the sample data). 
                     Note: For this to work, gtshark needs to be installed and in the execution path - it is a separate
                     software package that is not affliated with genozip in any way. It can be found here:
                     https://github.com/refresh-bio/GTShark
                     Note: gtshark also needs to be installed for decompressing files that were compressed with this
                     option. 
                     Note: This option isn't supported on Windows

One or more file names may be given, or if omitted, standard input is used instead

genozip is available for free for non-commercial use. Commercial use requires a commercial license

Bug reports and feature requests: bugs@genozip.com
Commercial license inquiries: sales@genozip.com

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

------------------------------------------------------------------------------------------------------------------------
Uncompress VCF (Variant Call Format) files previously compressed with genozip

Usage: genounzip [options]... [files]...

See also: genozip genocat genols

Options:
   -c --stdout       Send output to standard output instead of a file

   -f --force        Force overwrite of the output file

   -^ --replace      Replace the source file with the result file, rather than leaving it unchanged

   -O --split        Split a concatenated file back to its original components

   -o --output       <output-filename>. Output to this filename instead of the default one

   -p --password     <password>. Provide password to access file(s) that were compressed with --password

   -m --md5          Shows the MD5 hash of the decompressed VCF file. If the file was originally compressed with --md5,
                     it also verifies that the MD5 of the original VCF file is identical to the MD5 of the decompressed
                     VCF.
                     Note: for compressed files, e.g. myfile.vcf.gz, the MD5 calculated is that of the original,
                     uncompressed file. 
                     Note: if the VCF file has Windows-style \r\n line endings, the md5 will be that of the modified
                     file with the \r removed

   -q --quiet        Don't show the progress indicator or warnings

   -Q --noisy        The --quiet is turned on by default when outputting to the terminal. --noisy stops the suppression
                     of warnings

   -t --test         Decompresss in memory (i.e. without writing the decompressed file to disk) - comparing the MD5 of
                     the resulting decompressed file to that of the original VCF. Works only if the file was compressed
                     with --md5

   -@ --threads      <number>. Specify the maximum number of threads. By default, this is set to the number of cores
                     available. The number of threads actually used may be less, if sufficient to balance CPU and I/O

   -h --help         Show this help page. Use with -f to see developer options.

   -L --license      Show the license terms and conditions for this product

   -V --version      Display version number

One or more file names must be given

Bug reports and feature requests: bugs@genozip.com
Commercial license inquiries: sales@genozip.com

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

------------------------------------------------------------------------------------------------------------------------
View metadata of VCF (Variant Call Format) files previously compressed with genozip

Usage: genols [options]... [files or directories]...

See also: genozip genounzip genocat

Options:
   -q --quiet        Don't show warnings

   -h --help         Show this help page

   -L --license      Show the license terms and conditions for this product

   -V --version      Display version number

One or more file or directory names may be given, or if omitted, genols runs on the current directory

Bug reports and feature requests: bugs@genozip.com
Commercial license inquiries: sales@genozip.com

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

------------------------------------------------------------------------------------------------------------------------
Print VCF (Variant Call Format) file(s) previously compressed with genozip

Usage: genocat [options]... [files]...

See also: genozip genounzip genols

Options:
   -r --regions      [^]chr|chr:pos|pos|chr:from-to|chr:from-|chr:-to|from-to|from-|-to[,...]
                     Shows one or more regions of the file. Examples:
                               genocat myfile.vcf.genozip -r22:1000000-2000000  (A range of chromosome 22)
                               genocat myfile.vcf.genozip -r-2000000,2500000-   (Two ranges of all chromosomes)
                               genocat myfile.vcf.genozip -r21,22               (All of chromosome 21 and 22)
                               genocat myfile.vcf.genozip -r^MT,Y               (All of chromosomes except for MT and Y)
                               genocat myfile.vcf.genozip -r^-10000             (All sites on all chromosomes, except
                     positions up to 10000)
                     Note: genozip files are indexed automatically during compression. There is no separate indexing
                     step or separate index file
                     Note: Indels are considered part of a region if their start position is
                     Note: Multiple -r arguments may be specified - this is equivalent to chaining their regions with a
                     comma separator in a single argument

   -t --targets      Identical to --regions, provided for pipeline compatability

   -s --samples      [^]sample[,...]
                     Shows a subset of samples (individuals). Examples:
                               genocat myfile.vcf.genozip -s HG00255,HG00256    (show two samples)
                               genocat myfile.vcf.genozip -s ^HG00255,HG00256   (show all samples except these two)
                     Note: This does not change the INFO data (including the AC and AN tags)
                     Note: sample names are case-sensitive

                     Note: Multiple -s arguments may be specified - this is equivalent to chaining their samples with a
                     comma separator in a single argument

   -G --drop-genotypes Output the data without the individual genotypes and FORMAT column

   -H --no-header    Don't output the VCF header

      --header-only  Output only the VCF header

      --strip        Don't output values for ID, QUAL, FILTER, INFO; FORMAT is only GT (at most); Samples include
                     allele values (i.e. GT subfield) only

   -o --output       <output-filename>. Output to this filename instead of stdout

   -p --password     Provide password to access file(s) that were compressed with --password

   -@ --threads      Specify the maximum number of threads. By default, this is set to the number of cores available.
                     The number of threads actually used may be less, if sufficient to balance CPU and I/O

   -q --quiet        Don't show warnings

   -Q --noisy        The --quiet is turned on by default when outputting to the terminal. --noisy stops the suppression
                     of warnings

   -h --help         Show this help page. Use with -f to see developer options. Use --header-only if that's what you're
                     looking for

   -L --license      Show the license terms and conditions for this product

   -V --version      Display version number

One or more file names must be given

Bug reports and feature requests: bugs@genozip.com
Commercial license inquiries: sales@genozip.com

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

