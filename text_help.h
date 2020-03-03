// ------------------------------------------------------------------
//   help-text.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"

static const char *help_genozip[] = {
    "",
    "Compress VCF (Variant Call Format) files",
    "",
    "Usage: genozip [options]... [files]...",
    "",
    "See also: genounzip genocat genols",
    "",
    "Actions - use at most one of these actions:",
    "   -z --compress     This is the default action. Compress a .vcf, .vcf.gz or .vcf.bz2 file (Yes! You can compress an already-compressed file). The source file is left unchanged. This is the default action of genozip",
    "   -d --decompress   Same as running genounzip. For more details, run: genounzip --help",
    "   -l --list         Same as running genols. For more details, run: genols --help",
    "   -h --help         Show this help page. Use with -f to see developer options.",
    "   -L --license      Show the license terms and conditions for this product",
    "   -V --version      Display version number",
    "",    
    "Flags:",    
    "   -c --stdout       Send output to standard output instead of a file",    
    "   -f --force        Force overwrite of the output file, or force writing .vcf" GENOZIP_EXT " data to standard output",    
    "   -^ --replace      Replace the source file with the result file, rather than leaving it unchanged",    
    "   -o --output       <output-filename>. This option can also be used to concatenate multiple input files with the same individuals, into a single concatenated output file",
    "   -p --password     <password>. Password-protected - encrypted with 256-bit AES",
    "   -m --md5          Shows the MD5 hash of the VCF file. Note: genozip always calculates and stores the md5, and when the resulting file is decompressed, this MD5 will be compared to the MD5 of the decompressed VCF. Note: for compressed files, e.g. myfile.vcf.gz, the MD5 calculated is that of the original, uncompressed file",
    "   -q --quiet        Don't show the progress indicator or warnings",    
    "   -t --test         After compressing normally, decompresss in memory (i.e. without writing the decompressed file to disk) - comparing the MD5 of the resulting decompressed file to that of the original VCF",
    "   -@ --threads      <number>. Specify the maximum number of threads. By default, this is set to the number of cores available. The number of threads actually used may be less, if sufficient to balance CPU and I/O",
    "   --show-content    Show the information content of VCF files and the compression ratios of each component",
    "   --show-alleles    Output allele values to stdout. Each row corresponds to a row in the VCF file. Mixed-ploidy regions are padded, and 2-digit allele values are replaced by an ascii character",
    "",
    "Optimizing:",    
    "   -B --vblock       <number>. Sets the number of variants (sites) in a variant block. Memory consumption of both compression and decompression are linear with the variant block size. Compression is sometimes better with larger block sizes, in particular if the number of samples is small. By default, the variant block size is set to a value between 4096 (for VCF files with 1024 samples or more) and 131,072 (for VCF files with 31 samples or less). The genozip file is indexed on the variant block level, so smaller blocks will result in faster lookup using -r or -R, but at the cost of a larger index",    
    "",
    "One or more file names may be given, or if omitted, standard input is used instead",
    "",
    "genozip is available for free for non-commercial use. Commercial use requires a commercial license",
};

static const char *help_genozip_developer[] = {
    "Options useful mostly for developers of genozip:",
    "   --show-time       Show what functions are consuming the most time",
    "   --show-memory     Show what buffers are consuming the most memory",
    "   --show-sections   Shows the section types of the output genozip file and the compression ratios of each component",
    "   --show-dict       Show dictionary fragments written for each variant block (works for genounzip too)",
    "   --show-one-dict   <field-name> Show the dictionary for this field in a tab-seperated list - <field-name> may be one of the fields 1-9 (CHROM to FORMAT) or a subfield of INFO or a subfield of FORMAT (works for genounzip too)",
    "   --show-gt-nodes   Show transposed GT matrix - each value is an index into its dictionary",
    "   --show-b250       Show fields 1-9 (CHROM to FORMAT) as well as the subfields of INFO - each value shows the line (counting from 1) and the index into its dictionary (note: REF and ALT are compressed together as they are correlated). This also works with genounzip, but without the line numbers.",
    "   --show-one-b250   <field-name> Show the values for this field - maybe one of the fields 1-9 (CHROM to FORMAT) or a subfield of INFO",
    "   --show-headers    Show the sections headers (works for genounzip too)",
    "   --show-index      Show the content of the random access index"
    "   --show-gheader    Show the content of the genozip header (which also includes the list of all sections in the file)",
};

static const char *help_genounzip[] = {
    "",
    "Uncompress VCF (Variant Call Format) files previously compressed with genozip",
    "",
    "Usage: genounzip [options]... [files]...",
    "",
    "See also: genozip genocat genols",
    "",
    "Options:",
    "   -c --stdout       Send output to standard output instead of a file",
    "   -f --force        Force overwrite of the output file",
    "   -^ --replace      Replace the source file with the result file, rather than leaving it unchanged",    
    "   -O --split        Split a concatenated file back to its original components",
    "   -o --output       <output-filename>. This option can also be used to concatenate multiple input files with the same individuals, into a single concatenated output file",
    "   -p --password     <password>. Provide password to access file(s) that were compressed with --password",
    "   -m --md5          Shows the MD5 hash of the VCF file. Note: genozip always compares the MD5 of the original VCF file to the MD5 of the decompressed VCF. Note: for compressed files, e.g. myfile.vcf.gz, the MD5 calculated is that of the original, uncompressed file",
    "   -q --quiet        Don't show the progress indicator or warnings",    
    "   -t --test         Decompresss in memory (i.e. without writing the decompressed file to disk) - comparing the MD5 of the resulting decompressed file to that of the original VCF",
    "   -@ --threads      <number>. Specify the maximum number of threads. By default, this is set to the number of cores available. The number of threads actually used may be less, if sufficient to balance CPU and I/O",
    "   -h --help         Show this help page. Use with -f to see developer options.",
    "   -L --license      Show the license terms and conditions for this product",
    "   -V --version      Display version number",
    "",
    "One or more file names must be given"
};

static const char *help_genols[] = {
    "",
    "View metadata of VCF (Variant Call Format) files previously compressed with genozip",
    "",
    "Usage: genols [options]... [files or directories]...",
    "",
    "See also: genozip genounzip genocat",
    "",
    "Options:",
    "   -q --quiet        Don't show warnings",    
    "   -h --help         Show this help page",
    "   -L --license      Show the license terms and conditions for this product",
    "   -V --version      Display version number",
    "",
    "One or more file or directory names may be given, or if omitted, genols runs on the current directory",
};

static const char *help_genocat[] = {
    "",
    "Print VCF (Variant Call Format) file(s) previously compressed with genozip",
    "",
    "Usage: genocat [options]... [files]...",
    "",
    "See also: genozip genounzip genols",
    "",
    "Options:",    
    "   -r --regions      chr|chr:pos|pos|chr:from-to|chr:from-|chr:-to|from-to|from-|-to[,...]",
    "                     Shows one or more regions of the file. Examples:",
    "                               genocat myfile.vcf.genozip -r22:1000000-2000000  (A range of chromosome 22)",
    "                               genocat myfile.vcf.genozip -r-2000000,2500000-   (Two ranges of all chromosomes)",
    "                               genocat myfile.vcf.genozip -r21,22               (All of chromosome 21 and 22)",            
    "   -p --password     Provide password to access file(s) that were compressed with --password",
    "   -@ --threads      Specify the maximum number of threads. By default, this is set to the number of cores available. The number of threads actually used may be less, if sufficient to balance CPU and I/O",
    "   -q --quiet        Don't show warnings",    
    "   -h --help         Show this help page. Use with -f to see developer options.",
    "   -L --license      Show the license terms and conditions for this product",
    "   -V --version      Display version number",
    "",
    "One or more file names must be given",
};

static const char *help_footer[] = {
    "",
    "Bug reports and feature requests: bugs@genozip.com",
    "Commercial license inquiries: sales@genozip.com",
    "",
    "THIS SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.",
    ""
};

