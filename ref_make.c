// ------------------------------------------------------------------
//   ref_make.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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
#include "file.h"
#include "ref_iupacs.h"

SPINLOCK (make_ref_spin);
#define MAKE_REF_NUM_RANGES 1000000 // should be more than enough (in GRCh38 we have 6389)

static Buffer contig_metadata = {}; // contig header of each contig, except for chrom name

// called from ref_make_create_range, returns the range for this fasta VB. note that we have exactly one range per VB
// as txtfile_read_vblock makes sure we have only one full or partial contig per VB (if flag.make_reference)
static Range *ref_make_ref_get_range (uint32_t vblock_i)
{
    // access ranges.len under the protection of the mutex
    spin_lock (make_ref_spin);
    gref->ranges.len = MAX (gref->ranges.len, (uint64_t)vblock_i); // note that this function might be called out order (called from ref_make_create_range - FASTA ZIP compute thread)
    ASSERT (gref->ranges.len <= MAKE_REF_NUM_RANGES, "reference file too big - number of ranges exceeds %u", MAKE_REF_NUM_RANGES);
    spin_unlock (make_ref_spin);

    return ENT (Range, gref->ranges, vblock_i-1);
}

// called during REF ZIP compute thread, from zip_compress_one_vb (as "compress" defined in data_types.h)
// converts the vb sequence into a range
void ref_make_create_range (VBlockP vb)
{
    Range *r = ref_make_ref_get_range (vb->vblock_i);
    uint64_t seq_len = CTX(FASTA_NONREF)->local.len;

    // at this point, we don't yet know the first/last pos or the chrom - we just create the 2bit sequence array.
    // the missing details will be added during ref_make_prepare_range_for_compress
    r->ref = bit_array_alloc (seq_len * 2, false); // 2 bits per base
    r->range_id = vb->vblock_i-1;

    uint64_t bit_i=0;
    for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) {
        
        uint32_t seq_data_start, seq_len;
        fasta_get_data_line (vb, line_i, &seq_data_start, &seq_len);

        const uint8_t *line_seq = ENT (uint8_t, vb->txt_data, seq_data_start);
        for (uint64_t base_i=0; base_i < seq_len; base_i++) {
            char base = line_seq[base_i];
            uint8_t encoding = acgt_encode[(int)base];
            bit_array_assign (&r->ref, bit_i, encoding & 1);
            bit_array_assign (&r->ref, bit_i + 1, (encoding >> 1) & 1);

            // store very rare IUPAC bases (GRCh38 has 94 of them)
            ref_iupacs_add (vb, bit_i/2, base);

            bit_i += 2;
        }
    }

    ASSERT (seq_len * 2 == bit_i, "Expecting SEQ.local.len (x2 = %"PRId64") == bit_i (%"PRId64")", seq_len * 2, bit_i);
}

// in make_ref, each VB gets it own range indexed by vb->vblock_i - so they can work on them in parallel without
// worrying about byte overlap. called from zip_one_file as zip_initialize
void ref_make_ref_init (void)
{
    ASSERT0 (flag.make_reference, "Expecting flag.make_reference=true");

    // remove old cache files
    ref_remove_cache (gref);
    refhash_remove_cache();

    buf_alloc (evb, &gref->ranges, 0, MAKE_REF_NUM_RANGES, Range, 1, "ranges"); // must be allocated by main thread as its evb
    ranges_type (gref) = RT_MAKE_REF;
    
    buf_zero (&gref->ranges);

    refhash_initialize (NULL);

    spin_initialize (make_ref_spin);
}


// the "read" part of reference-compressing dispatcher, called from ref_compress_ref
void ref_make_prepare_range_for_compress (VBlockP vb)
{
    if (vb->vblock_i-1 == gref->ranges.len) return; // we're done

    Range *r = ENT (Range, gref->ranges, vb->vblock_i-1); // vb_i=1 goes to ranges[0] etc

    // we have exactly one contig for each VB (but possibly multiple VBs with the same contig), one one RAEntry for that contig
    // during seg we didn't know the chrom,first,last_pos, so we add them now, from the RA
    random_access_get_ra_info (vb->vblock_i, &r->chrom, &r->first_pos, &r->last_pos);
    
    // set gpos "global pos" - a single 0-based coordinate spanning all ranges in order
    r->gpos = vb->vblock_i > 1 ? (r-1)->gpos + ref_size (r-1) : 0;

    // each chrom's gpos must start on a 64bit aligned word
    if (r->gpos % 64 && r->chrom != (r-1)->chrom) // each new chrom needs to have a GPOS aligned to 64, so that we can overload is_set bits between the whole genome and individual chroms
        r->gpos = ROUNDUP64 (r->gpos);

    vb->range             = r; // range to compress
    vb->range->num_set    = r->ref.nbits / 2;
    vb->ready_to_dispatch = true;
}

// make-refernece called by main thread after completing compute thread of VB 
void ref_make_after_compute (VBlockP vb_)
{
    VBlockFASTA *vb = (VBlockFASTA *)vb_;

    // add contig metadata
    if (vb->has_contig_metadata) {
        buf_alloc (evb, &contig_metadata, 1, 1000, ContigMetadata, 2, "contig_metadata");
        NEXTENT (ContigMetadata, contig_metadata) = vb->contig_metadata;
    }

    // collect iupacs
    ref_iupacs_after_compute (vb_);
}

ConstBufferP ref_make_get_contig_metadata (void) 
{ 
    return &contig_metadata; 
}

void ref_make_finalize (void) 
{ 
    buf_free (&contig_metadata); 
}

