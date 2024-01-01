// ------------------------------------------------------------------
//   ref_make.c
//   Copyright (C) 2020-2024 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "data_types.h"
#include "codec.h" // must be included before reference.h
#include "reference.h"
#include "vblock.h"
#include "fasta_private.h"
#include "ref_private.h"
#include "mutex.h"
#include "refhash.h"
#include "random_access.h"
#include "context.h"
#include "filename.h"
#include "file.h"
#include "ref_iupacs.h"
#include "contigs.h"
#include "stream.h"
#include "arch.h"
#include "digest.h"

SPINLOCK (make_ref_spin);
#define MAKE_REF_NUM_RANGES 1000000 // should be more than enough (in GRCh38 we have 6389)

static Buffer contig_metadata = {}; // contig header of each contig, except for chrom name

Serializer make_ref_merge_serializer = {};

bool is_ref (STRp(data), bool *need_more/*optional*/)
{
    // this file can be segged into DT_REF only if --make-reference is specified and it is FASTA
    return flag.make_reference && is_fasta (STRa(data), need_more);
}

void ref_make_seg_initialize (VBlockP vb)
{
    START_TIMER;

    ASSINP (vb->vblock_i > 1 || *B1STtxt == '>' || *B1STtxt == ';',
            "Error: expecting FASTA file %s to start with a '>' or a ';'", txt_name);

    CTX(FASTA_CONTIG)->no_vb1_sort = true; // keep contigs in the order of the reference, i.e. in the order they would appear in BAM header created with this reference 
    CTX(FASTA_CONTIG)->no_stons = true; // needs b250 node_index for reference

    if (segconf.running) segconf.fasta_has_contigs = true; // initialize optimistically

    COPY_TIMER (seg_initialize);
}

// called from ref_make_create_range, returns the range for this fasta VB. note that we have exactly one range per VB
// as txtfile_read_vblock makes sure we have only one full or partial contig per VB (if flag.make_reference)
static Range *ref_make_ref_get_range (VBIType vblock_i)
{
    // access ranges.len under the protection of the mutex
    spin_lock (make_ref_spin);
    gref->ranges.len32 = MAX_(gref->ranges.len32, vblock_i); // note that this function might be called out order (called from ref_make_create_range - FASTA ZIP compute thread)
    ASSERT (gref->ranges.len <= MAKE_REF_NUM_RANGES, "reference file too big - number of ranges exceeds %u", MAKE_REF_NUM_RANGES);
    spin_unlock (make_ref_spin);

    return B(Range, gref->ranges, vblock_i-1);
}

// called during REF ZIP compute thread, from zip_compress_one_vb (as "zip_after_compress" defined in data_types.h)
// converts the vb sequence into a range
void ref_make_create_range (VBlockP vb)
{
    Range *r = ref_make_ref_get_range (vb->vblock_i);
    uint64_t seq_len = CTX(FASTA_NONREF)->local.len;

    // at this point, we don't yet know the first/last pos or the chrom - we just create the 2bit sequence array.
    // the missing details will be added during ref_make_prepare_one_range_for_compress
    r->ref = bits_alloc (seq_len * 2, false); // 2 bits per base
    r->range_id = vb->vblock_i-1;

    uint64_t bit_i=0;
    for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) {
        
        uint32_t seq_data_start, seq_len;
        fasta_get_data_line (vb, line_i, &seq_data_start, &seq_len);

        bytes line_seq = B8 (vb->txt_data, seq_data_start);
        for (uint64_t base_i=0; base_i < seq_len; base_i++, bit_i += 2) {
            char base = line_seq[base_i];
            bits_assign2 (&r->ref, bit_i, acgt_encode[(int)base]); 

            // store very rare IUPAC bases (GRCh38 has 94 of them)
            ref_iupacs_add (vb, bit_i/2, base);
        }
    }

    ASSERT (seq_len * 2 == bit_i, "Expecting SEQ.local.len (x2 = %"PRId64") == bit_i (%"PRId64")", seq_len * 2, bit_i);
}

// in make_ref, each VB gets it own range indexed by vb->vblock_i - so they can work on them in parallel without
// worrying about byte overlap. called from zip_one_file as zip_initialize
void ref_make_ref_init (void)
{
    ASSERT0 (flag.make_reference, "Expecting flag.make_reference=true");

    buf_alloc (evb, &gref->ranges, 0, MAKE_REF_NUM_RANGES, Range, 1, "ranges"); // must be allocated by main thread as its evb
    gref->ranges.rtype = RT_MAKE_REF;
    
    buf_zero (&gref->ranges);

    refhash_initialize_for_make();

    spin_initialize (make_ref_spin);

    serializer_initialize (make_ref_merge_serializer);
}

void ref_make_prepare_ranges_for_compress (void)
{
    for_buf2 (Range, r, i, gref->ranges) { 
        // we have exactly one contig for each VB (but possibly multiple VBs with the same contig), one one RAEntry for that contig
        // during seg we didn't know the chrom,first,last_pos, so we add them now, from the RA
        random_access_get_ra_info (i+1, &r->chrom, &r->first_pos, &r->last_pos);

        // set gpos "global pos" - a single 0-based coordinate spanning all ranges in order        
        // each chrom's gpos be 64-base aligned (so that is_set, which has 1 bit per base, can be word-aligned for each chrom) 
        r->gpos = !i                         ? 0                                         // first range
                : (r->chrom != (r-1)->chrom) ? ROUNDUP64 ((r-1)->gpos + ref_size (r-1))  // first range of a contig
                :                              (r-1)->gpos + ref_size (r-1);             // 2nd+ range of contig
    }

    PosType64 max_gpos = BLST(Range, gref->ranges)->gpos;

    // since our refhash entries are 32 bit, we cannot use the reference data beyond the first 4Gbp for creating the refhash
    // TO DO: make the hash entries 40bit (or 64 bit?) if genome size > 4Gbp (bug 150)
    ASSERTW (max_gpos <= MAX_ALIGNER_GPOS,
             "FYI: %s contains %s bases. When compressing a FASTQ or unaligned (i.e. missing RNAME, POS) "
             "SAM/BAM file using the reference being generated, only the first %s bases of the reference will be used, "
             "possibly affecting the compression ratio. This limitation doesn't apply to aligned SAM/BAM files and VCF files. "
             "If you need to use reference files larger than %s, let us know at " EMAIL_SUPPORT ".", 
             txt_name, str_bases (max_gpos).s, str_bases (MAX_ALIGNER_GPOS).s, str_bases (MAX_ALIGNER_GPOS).s);
}

// the "read" part of reference-compressing dispatcher, called from ref_compress_ref
void ref_make_prepare_one_range_for_compress (VBlockP vb)
{
    if (vb->vblock_i-1 == gref->ranges.len32) return; // we're done

    RangeP r = B(Range, gref->ranges, vb->vblock_i-1); // vb_i=1 goes to ranges[0] etc

    vb->range          = r; // range to compress
    vb->range->num_set = r->ref.nbits / 2;
    vb->dispatch       = READY_TO_COMPUTE;
}

// ZIP main thread: mimic the logic of loading a reference, and calculate the digest
void ref_make_calculate_digest (void)
{
    START_TIMER;

    // mirrors the logic in ref_initialize_ranges
    Range *last_r = BLST(Range, gref->ranges);
    gref->genome_nbases = ROUNDUP64 (last_r->gpos + ref_size (last_r)) + 64;

    gref->genome_buf.can_be_big = true; // supress warning in case of an extra large genome (eg plant genomes)
    gref->genome = buf_alloc_bits (evb, &gref->genome_buf, gref->genome_nbases * 2, CLEAR, 0, "ref->genome_buf");

    for_buf (Range, r, gref->ranges) 
        bits_copy (gref->genome, r->gpos * 2, &r->ref, 0, ref_size(r) * 2);

    z_file->digest = digest_do (STRb(gref->genome_buf), !flag.md5, "genome");
    
    buf_free (gref->genome_buf);

    COPY_TIMER_EVB (ref_make_calculate_digest);
}

// zip_after_compute callback: make-reference called by main thread after completing compute thread of VB.
// must be called in order of VBs so that contigs are in the same order as the FASTA, resulting in GPOS being allocated
// to contigs consistently across executions.
void ref_make_after_compute (VBlockP vb_)
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;

    // add contig metadata
    if (vb->has_contig_metadata) {
        buf_alloc (evb, &contig_metadata, 1, 1000, ContigMetadata, 2, "contig_metadata");
        BNXT (ContigMetadata, contig_metadata) = vb->contig_metadata;
    }

    // collect iupacs
    ref_iupacs_after_compute (vb_);
}

ConstBufferP ref_make_get_contig_metadata (void) 
{ 
    return &contig_metadata; 
}

// callback from zfile_compress_genozip_header
void ref_make_genozip_header (SectionHeaderGenozipHeaderP header)
{
    header->genome_digest = z_file->digest;
}

// zip_finalize callback
void ref_make_finalize (bool is_last_user_txt_file) 
{ 
    buf_destroy (contig_metadata); 
    serializer_destroy (make_ref_merge_serializer);
}

// Get reference file name from FASTA name, and if reference file does not exist, run a separate process to --make-reference
rom ref_fasta_to_ref (FileP file)
{
    rom ref_filename = filename_z_normal (file->name, DT_REF, file->type);

    // if file reference doesn't exist yet - --make-reference now, in a separate process
    if (!file_exists (ref_filename)) {

        WARN ("FYI: cannot find reference file %s: generating it now from %s", ref_filename, file->name);

        StreamP make_ref = stream_create (NULL, 0, 0, 0, 0, 0, 0, "Make reference",
                                          arch_get_executable().s, "--make-reference", file->name, 
                                          flag.stats_submit ? "--submit" : SKIP_ARG,
                                          "--no-tip", NULL);
        
        // wait for child process to finish
        ASSINP (!stream_wait_for_exit (make_ref), "Failed to make reference file %s. Try making it explicitly with \"genozip --make-reference %s\"", 
                ref_filename, file->name); 
    }

    FREE (file->name);
    file->name = (char *)ref_filename;

    Reference ref = flag.reading_reference ? flag.reading_reference : gref;

    REALLOC ((char **)&ref->filename, strlen (ref_filename) + 1, "ref->filename");
    strcpy ((char*)ref->filename, ref_filename);

    file->disk_size = file_get_size (ref_filename);

    return ref->filename;
}
