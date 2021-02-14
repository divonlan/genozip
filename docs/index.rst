Welcome to Genozip
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Welcome to Genozip
   genozip
   genounzip
   genocat
   genols
   developer options <developer>
   Publications
   License <license>
     
Genozip is a universal compressor for genomic files - it is optimized to compress FASTQ, SAM/BAM/CRAM, VCF/BCF, FASTA, GVF, Phylip and 23andMe files, but it can also compress any other file (including non-genomic files). 

It can even compress files that are already compressed with .gz .bz2 .xz (for full list of supported file types see ``genozip --help=input``).

The compression ratio one can expect to see depends on the data being compressed - it is typically about a 1.5-3X ratio when compressing .bam, 2X-5X for .fastq.gz files (i.e. compressing already-compressed files), and up to 200X when compressing an uncompressed high-sample-count .vcf file with only GT data or a Multi-FASTA file with similar contigs.

The compression is lossless - the decompressed file is 100% identical to the original file (except when using the ``--optimize`` option).

While Genozip was designed primarily to allow cost-effective storage and fast transmission of genomic files, it also includes extensive options that allow direct use of genozip-compressed files in analytical pipelines.

Genozip consists of four command line tools:
  
   - :doc:`genozip` compresses files
  
   - :doc:`genounzip` decompresses files

   - :doc:`genocat` is the workhorse for using genozip in analytical pipelines:
      - Display the contents of a compressed file - possibly piping it into a downstream tool
      - Subset a compressed file - show a specific part of its contents
      - Translate a compressed file to another format (eg BAM to FASTQ or Multi-FASTA to Phylip)
      - Show various statistics and metadata related to a compressed file
  
   - :doc:`genols` is similar to ls, just for compressed files

Installing
==========

Citing
======
Lan, D., et al. Bioinformatics, 36, 4091–4092, July 2020, https://doi.org/10.1093/bioinformatics/btaa290

Contact
=======
Bug reports and feature requests: bugs@genozip.com

Commercial license inquiries: sales@genozip.com

Requests for support for compression of additional public or proprietary genomic file formats: sales@genozip.com

License
=======
Genozip's :doc:`license` allows for free non-commerical use, subject to certain conditions. For a commercial license, please contact sales@genozip.com.
   
THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
