Genozip
=======

*A universal compressor for genomic files*
      
.. toctree::
   :maxdepth: 2

   Installing <installing> 
   Examples <examples> 
   Use cases <use-cases>
   Using genozip <genozip>
   Using genounzip <genounzip>
   Using genocat <genocat>
   Using genols <genols>
   Publications & Citing <publications>
   Source code <source>
   Developer options <developer>
   License <license>
   Contact <contact>
  
|
      
About Genozip
=============

`Sign up <https://tinyurl.com/genozip>`_ to receive low-frequency updates related to Genozip.

Genozip is a universal compressor for genomic files - it is optimized to compress FASTQ, SAM/BAM/CRAM, VCF/BCF, FASTA, GVF, Phylip and 23andMe files, but it can also compress any other file (including non-genomic files). 

Typically, a **2X-5X improvement over the existing compression** is achieved when compressing already-compressed files like .fastq.gz .bam vcf.gz, and much higher ratios in some other cases.

Yes, Genozip can compress already-compressed files (.gz .bz2 .xz .bam .cram).

The compression is lossless - the decompressed file is 100% identical to the original file (except when using the ``--optimize`` option).

Genozip consists of four command line tools:
  
   * :doc:`genozip` compresses files
  
   * :doc:`genounzip` decompresses files

   * :doc:`genols` shows metadata of compressed files and directories

   * :doc:`genocat` is the workhorse for using genozip in analytical pipelines:
      - Display the contents of a compressed file - possibly piping it into a downstream tool
      - Subset a compressed file - show a specific part of its contents
      - Translate a compressed file to another format (eg BAM to FASTQ or Multi-FASTA to Phylip)
      - Show various statistics and metadata related to a compressed file
  
|
 
.. include:: installing.rst
.. include:: publications.rst
.. include:: contact.rst

License
=======
Genozip's :doc:`license` allows for free non-commerical use, subject to certain conditions. For a commercial license, please contact sales@genozip.com.
   
THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
