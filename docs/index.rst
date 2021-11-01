Genozip
=======

*A universal compressor for genomic files*
      
.. toctree::
   :maxdepth: 3

   Installing <installing> 
   File Types <data-types>
   Manual <manual>
   Capabilities <capabilities>
   Applications <applications> 
   Algorithms <algorithms>
   Source code <source>
   Publications & Citing <publications>
   License <license>
   Release Notes <release-notes>
   Contact <contact>
  
|
      
About Genozip
=============

`Sign up <https://tiny.cc/genozip>`_ to receive low-frequency updates related to Genozip.

Genozip is a universal compressor for genomic files - it is optimized to compress FASTQ, SAM/BAM/CRAM, VCF/BCF, FASTA, GFF3/GVF, PHYLIP, Chain, Kraken and 23andMe files, but it can also compress any other file (including non-genomic files). 

Typically, a **2X-5X improvement over the existing compression** is achieved when compressing already-compressed files like .fastq.gz .bam vcf.gz, and much higher ratios in some other cases.

**Yes**, Genozip can compress already-compressed files (.gz .bz2 .xz .bam .cram).

The compression is **lossless** - the decompressed file is 100% identical to the original file (some :ref:`exceptions<losslessness>` apply).

Genozip consists of four command line tools:
  
   * :doc:`genozip` compresses files
  
   * :doc:`genounzip` decompresses files

   * :doc:`genols` shows metadata of compressed files and directories

   * :doc:`genocat` is the workhorse for using genozip in analytical pipelines:
      - Display the contents of a compressed file - possibly piping it into a downstream tool
      - Subset a compressed file - show a specific part of its contents
      - Translate a compressed file to another format (eg BAM to FASTQ or Multi-FASTA to Phylip)
      - Analyze a compressed file (eg showing the :ref:`sex<sex>`, :ref:`coverage<coverage>` or compression statistics)
  
|
 
.. include:: installing.rst
.. include:: publications.rst
.. include:: contact.rst

License
=======
Genozip's :doc:`license` allows for free non-commerical use, subject to certain conditions. For a commercial license, please contact sales@genozip.com.
   
THIS SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS, COPYRIGHT HOLDERS OR DISTRIBUTORS OF THIS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
