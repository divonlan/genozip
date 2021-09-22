// ------------------------------------------------------------------
//   reference.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <errno.h>
#include "reference.h"
#include "buffer.h"
#include "strings.h"
#include "dict_id.h"
#include "dispatcher.h"
#include "zfile.h"
#include "endianness.h"
#include "random_access.h"
#include "seg.h"
#include "piz.h"
#include "vblock.h"
#include "context.h"
#include "mutex.h"
#include "sections.h"
#include "file.h"
#include "regions.h"
#include "refhash.h"
#include "ref_private.h"
#include "compressor.h"
#include "threads.h"
#include "website.h"
#include "ref_iupacs.h"
#include "chrom.h"
#include "crypt.h"

static RefStruct refs[2] = { { .ctgs.name = "gref" }, { .ctgs.name = "prim_ref"} };
Reference gref     = &refs[0]; // global reference 
Reference prim_ref = &refs[1]; // chain file primary coordinates reference

#define CHROM_GENOME 0
#define CHROM_NAME_GENOME "GENOME"

#define CHROM_GENOME_REV 1
#define CHROM_NAME_GENOME_REV "GENOME_REV"

#define ref_is_range_used(r) ((r)->ref.nbits && ((r)->is_set.nbits || flag.make_reference))

// get / set functions 
const char *ref_get_filename (const Reference ref) { return ref->filename; }
uint8_t ref_get_genozip_version (const Reference ref) { return ref->genozip_version; }
BufferP ref_get_stored_ra (Reference ref) { return &ref->stored_ra; }
Digest ref_get_file_md5 (const Reference ref) { return ref->file_md5; }
ConstContigPkgP ref_get_ctgs (const Reference ref) { return &ref->ctgs; }
uint32_t ref_num_contigs (Reference ref) { return (uint32_t)ref->ctgs.contigs.len; }

void ref_get_genome (Reference ref, const BitArray **genome, const BitArray **emoneg, PosType *genome_nbases)
{
    ASSERT0 (ref->genome, "Reference file not loaded");
    
    *genome = ref->genome;
    *emoneg = ref->emoneg;
    *genome_nbases = ref->genome_nbases;
}

void ref_set_genome_is_used (Reference ref, PosType gpos, uint32_t len)
{
    if (len == 1)
        bit_array_set (ref->genome_is_set, gpos); 
    else 
        bit_array_set_region (ref->genome_is_set, gpos, len);
}

static inline bool ref_has_is_set (void)
{
    return primary_command == PIZ || (primary_command == ZIP && flag.reference == REF_EXT_STORE);
}

static void ref_free_denovo_ranges (Reference ref)
{
    if (!buf_is_alloc (&ref->ranges)) return;

    ARRAY (Range, r, ref->ranges);
    for (unsigned i=0; i < ref->ranges.len ; i++) {
        FREE (r[i].ref.words);
        FREE (r[i].is_set.words);
        if (primary_command == ZIP) 
            FREE (r[i].chrom_name); // allocated only in ZIP/REF_INTERNAL - otherwise a pointer into another Buffer
    }
}

// free memory allocations between files, when compressing multiple non-bound files or decompressing multiple files
void ref_unload_reference (Reference ref)
{
    ref_create_cache_join (ref, true); // wait for cache to finish writing, if applicable

    if (ranges_type(ref) == RT_DENOVO) 
        ref_free_denovo_ranges (ref);
    
    // case: the reference has been modified and we can't use it for the next file
    if (flag.reference & REF_STORED) {
        buf_free (&ref->genome_buf);
        buf_free (&ref->emoneg_buf);
        buf_free (&ref->genome_cache);
        buf_free (&ref->ranges);
    }

    // case: free, unless these buffers are immutable so the next file can use them
    if (flag.reference == REF_INTERNAL || flag.reference == REF_STORED /* PIZ */) {
        buf_free (&ref->ref_external_ra);
        buf_free (&ref->ref_file_section_list);
        buf_free (&ref->genome_is_set_buf);
        contigs_free (&ref->ctgs);
        ref->genome_nbases = 0;

        // locks
        ref_lock_free (ref);
    }
    
     // in ZIP with REF_EXTERNAL, we just cleanup is_set
    if (flag.reference == REF_EXTERNAL && command == ZIP) 
        buf_zero (&ref->genome_is_set_buf);
    
    buf_free (&ref->region_to_set_list);
    buf_free (&ref->stored_ra);

    ref->external_ref_is_loaded = false;
}

void ref_destroy_reference (Reference ref, bool destroy_only_if_not_mmap)
{
    if (destroy_only_if_not_mmap && ranges_type(ref) == RT_CACHED) return;

    ref_create_cache_join (ref, true); // wait for cache to finish writing, if applicable

    if (ranges_type(ref) == RT_DENOVO) ref_free_denovo_ranges (ref);

    buf_destroy (&ref->ranges);
    buf_destroy (&ref->genome_buf);
    buf_destroy (&ref->emoneg_buf);
    buf_destroy (&ref->genome_cache);
    buf_destroy (&ref->genome_is_set_buf);
    buf_destroy (&ref->region_to_set_list);
    buf_destroy (&ref->ref_external_ra);
    buf_destroy (&ref->stored_ra);
    buf_destroy (&ref->ref_file_section_list);
    FREE (ref->filename);
    FREE (ref->cache_fn);

    // ref contig stuff
    contigs_destroy (&ref->ctgs);

    // iupacs stuff
    buf_destroy (&ref->iupacs_buf);

    ref_lock_free (ref);
    
    refhash_destroy();

    memset (ref, 0, sizeof (RefStruct));
}

// account for memory allocations NOT through Buffers 
MemStats ref_memory_consumption (Reference ref)
{
    MemStats stats = { "reference", 0, 0 };

    if (ranges_type(ref) == RT_DENOVO) {
        ARRAY (Range, r, ref->ranges);
        for (unsigned i=0; i < ref->ranges.len; i++) {
            if (r[i].ref.nwords) {
                stats.bytes += r[i].ref.nwords * sizeof (word_t);
                stats.buffers++;
            }
            if (r[i].is_set.nwords) {
                stats.bytes += r[i].is_set.nwords * sizeof (word_t);
                stats.buffers++;
            }
        }
    }

    return stats;
}

// PIZ: returns a range which is the entire contig
const Range *ref_piz_get_range (VBlockP vb, Reference ref, PosType first_pos_needed, uint32_t num_nucleotides_needed)
{
    ASSERTISALLOCED (ref->ranges);

    // caching
    if (vb->prev_range[0] && vb->prev_range_chrom_node_index[0] == vb->chrom_node_index)
        return vb->prev_range[0];

    // gets the index of the matching chrom in the reference - either its the chrom itself, or one with an alternative name
    // eg 'chr22' instead of '22'
    WordIndex ref_contig_index = chrom_2ref_piz_get (vb->chrom_node_index);

    ASSPIZ (ref_contig_index != WORD_INDEX_NONE, "Requested a reference range of a contig \"%.s\" with no reference.", ctx_get_words_snip (CHROM, vb->chrom_node_index));

    ASSPIZ (ref_contig_index < ref->ranges.len, "Expecting ref_contig_index=%d < ref->ranges.len=%"PRIu64 " in %s. FYI: z_file->chrom2ref_map.len=%"PRIu64, 
            ref_contig_index, ref->ranges.len, ref_get_filename(ref), z_file->chrom2ref_map.len);

    Range *r = ENT (Range, ref->ranges, ref_contig_index);
    if (!r->ref.nwords) return NULL; // this can ligitimately happen if entire chromosome is verbatim in SAM, eg. unaligned (pos=4) or SEQ or CIGAR are unavailable

    if (first_pos_needed + num_nucleotides_needed - 1 <= r->last_pos) {
        // this can happen if the reference originated from REF_INTERNAL, and the latter part of the range requested
        // originated from a missing range due to hash contention, and this missing part also happens to be
        // at the end of the chromosome, thereby causing r->last_pos to be less than the length of the chromosome. 
        // we return the range anyway, as the missing parts will have is_set=0 and seq_bitmap=0, retrieving from nonref

        // TODO: r->last_pos in REF_INTERNAL should include the missing ranges at the end of the chromosome
    }

    vb->prev_range[0] = r;
    vb->prev_range_chrom_node_index[0] = vb->chrom_node_index;

    return r;
}

// -------------------------------------------------------------------------------------------------------
// PIZ: read and uncompress stored ranges (originally produced with --REFERENCE or SAM internal reference)
// -------------------------------------------------------------------------------------------------------

// PIZ: uncompact a region within ref - called by compute thread of reading the reference
static void ref_uncompact_ref (Range *r, int64_t first_bit, int64_t last_bit, const BitArray *compacted)
{
    uint64_t start_1_offset=first_bit, start_0_offset, len_1; // coordinates into r->is_set (in nucleotides)
    uint64_t next_compacted=0; // coordinates into compacted (in nucleotides)

    while (1) {
        // find length of set region
        bool has_any_bit = bit_array_find_next_clear_bit (&r->is_set, start_1_offset, &start_0_offset);
        if (!has_any_bit || start_0_offset > last_bit) 
            start_0_offset = last_bit + 1; // this is the last region of 1s

        len_1 = start_0_offset - start_1_offset;
        ASSERT (len_1 > 0, "len_1 is not positive: start_0_offset=%"PRId64" start_1_offset=%"PRId64" first_bit=%"PRId64" last_bit=%"PRId64, 
                start_0_offset, start_1_offset, first_bit, last_bit);

        // do actual uncompacting
        bit_array_copy (&r->ref, start_1_offset * 2, compacted, next_compacted * 2, len_1 * 2);
        next_compacted += len_1;

        if (start_0_offset > last_bit) break; // we're done (we always end with a region of 1s because we removed the flanking 0s during compacting)

        // skip the clear region
        has_any_bit = bit_array_find_next_set_bit (&r->is_set, start_0_offset, &start_1_offset); 
        ASSERT0 (has_any_bit, "cannot find next set bit");
        ASSERT (start_1_offset <= last_bit, "expecting start_1_offset(%"PRId64") <= last_bit(%"PRId64")",
                start_1_offset, last_bit); // we removed the flanking regions, so there is always an 1 after a 0 within the region
    }

    ASSERT (next_compacted * 2 == compacted->nbits, "expecting next_compacted(%"PRId64") * 2 == compacted->nbits(%"PRId64")",
            next_compacted, compacted->nbits);
}

// Compute thread: called by ref_uncompress_one_range
Range *ref_get_range_by_chrom (Reference ref, WordIndex chrom, const char **chrom_name)
{
    Context *ctx = ZCTX(CHROM);
    ASSERT (chrom >= 0 && chrom < ctx->word_list.len, "chrom=%d out of range - ctx->word_list.len=%u",
            chrom, (uint32_t)ctx->word_list.len);

    if (chrom_name)
        *chrom_name = ctx_get_words_snip (ctx, chrom);

    ASSERT (chrom < ref->ranges.len, "expecting chrom=%d < ranges.len=%"PRIu64, chrom, ref->ranges.len);
    
    Range *r = ENT (Range, ref->ranges, chrom); // in PIZ, we have one range per chrom
    return r;
}

Range *ref_get_range_by_ref_index (VBlockP vb, Reference ref, WordIndex ref_contig_index)
{
    if (ref_contig_index == WORD_INDEX_NONE || ref_contig_index >= ref->ranges.len) return NULL;

    return ENT (Range, ref->ranges, ref_contig_index); // in PIZ, we have one range per chrom
}

// Print this array to a file stream.  Prints '0's and '1'.  Doesn't print newline.
static void ref_print_bases (FILE *file, const BitArray *bitarr, 
                             bit_index_t start_base, bit_index_t num_of_bases, bool is_forward)
{
    static const char fwd[2][2] = { { 'A', 'C' }, {'G', 'T'} };
    static const char rev[2][2] = { { 'T', 'G' }, {'C', 'A'} };

#define BASES_PER_LINE 100

    if (is_forward)
        for (bit_index_t i=start_base*2; i < (start_base + num_of_bases)*2; i+=2) {
            if (!flag.sequential && (i-start_base*2) % (BASES_PER_LINE*2) == 0)
                fprintf (file, "%8"PRIu64": ", i/2);
            fputc (fwd[bit_array_get(bitarr, i+1)][bit_array_get(bitarr, i)], file);
            if (!flag.sequential && ((i-start_base*2) % (BASES_PER_LINE*2) == 2*(BASES_PER_LINE-1))) fputc ('\n', file);
        }
    
    else 
        for (int64_t i=(start_base+num_of_bases-1)*2; i >= start_base*2; i-=2) { // signed type
            fputc (rev[bit_array_get(bitarr, i+1)][bit_array_get(bitarr, i)], file);
            if (!flag.sequential && (((start_base+num_of_bases-1)*2-i) % (BASES_PER_LINE*2) == (BASES_PER_LINE-1)*2)) fputc ('\n', file);
        }
    
    fputc ('\n', file);
}

static void ref_show_sequence (Reference ref)
{
    for (uint32_t range_i=0; range_i < ref->ranges.len; range_i++) {
        Range *r = ENT (Range, ref->ranges, range_i);

        // get first pos and last pos, potentially modified by --regions
        PosType first_pos, last_pos;
        bool revcomp; // to do: implement
        if (!r->ref.nbits ||
            !regions_get_range_intersection (r->chrom, r->first_pos, r->last_pos, 0, &first_pos, &last_pos, &revcomp)) continue;

        if (r->ref.nbits) {
            iprintf ("%.*s\n", r->chrom_name_len, r->chrom_name);
            ref_print_bases (info_stream, &r->ref, first_pos-1, last_pos-first_pos, true);
        }
    }

    if (exe_type == EXE_GENOCAT) exit_ok();  // in genocat this, not the data
}

// entry point of compute thread of reference decompression. this is called when pizzing a file with a stored reference,
// including reading the reference file itself.
// vb->z_data contains a SEC_REFERENCE section and sometimes also a SEC_REF_IS_SET section
static void ref_uncompress_one_range (VBlockP vb)
{
    START_TIMER;

    if (!buf_is_alloc (&vb->z_data) || !vb->z_data.len) goto finish; // we have no data in this VB because it was skipped due to --regions or genocat --show-headers

    SectionHeaderReference *header = (SectionHeaderReference *)vb->z_data.data;

    WordIndex chrom          = (WordIndex)BGEN32 (header->chrom_word_index);
    uint32_t uncomp_len      = BGEN32 (header->h.data_uncompressed_len);
    PosType ref_sec_pos      = (PosType)BGEN64 (header->pos);
    PosType ref_sec_gpos     = (PosType)BGEN64 (header->gpos); // this is equal to sec_start_gpos. However up to 12.0.3 we had a bug in case of compacted ranges in a SAM/BAM DENOVO reference with start flanking regions - GPOS in the section header of a didn't reflect the flanking removal, so header->gpos cannot be trusted as correct for older SAM/BAM DENOVO reference files
    PosType ref_sec_len      = (PosType)BGEN32 (header->num_bases);
    PosType ref_sec_last_pos = ref_sec_pos + ref_sec_len - 1;
    PosType compacted_ref_len=0, initial_flanking_len=0, final_flanking_len=0; 

    ASSERT0 (chrom != WORD_INDEX_NONE, "Unexpected reference section with chrom=WORD_INDEX_NONE");

    const char *chrom_name;
    Range *r = ref_get_range_by_chrom (vb->ref, chrom, &chrom_name);
    PosType sec_start_within_contig = ref_sec_pos - r->first_pos;
    PosType sec_start_gpos          = r->gpos + sec_start_within_contig;
    PosType sec_end_within_contig   = sec_start_within_contig + ref_sec_len - 1;
    
    bool is_compacted = (header->h.section_type == SEC_REF_IS_SET); // we have a SEC_REF_IS_SET if SEC_REFERENCE was compacted

    if (flag.show_reference && primary_command == PIZ && r)  // in ZIP, we show the compression of SEC_REFERENCE into z_file, not the uncompression of the reference file
        iprintf ("vb_i=%u Uncompressing %-14s chrom=%u ('%.*s') gpos=%"PRId64" pos=%"PRId64" num_bases=%u comp_bytes=%u\n", 
                 vb->vblock_i, st_name (header->h.section_type), BGEN32 (header->chrom_word_index), r->chrom_name_len, r->chrom_name, BGEN64 (header->gpos), 
                 BGEN64 (header->pos), BGEN32 (header->num_bases), BGEN32 (header->h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));

    // initialization of is_set:
    // case 1: ZIP (reading an external reference) - we CLEAR is_set, and let seg set the bits that are to be
    //    needed from the reference for pizzing (so we don't store the others)
    //    note: in case of ZIP with REF_INTERNAL, we CLEAR the bits in ref_seg_get_locked_range
    // case 2: PIZ, reading an uncompacted (i.e. complete) reference section - which is always the case when
    //    reading an external reference and sometimes when reading a stored one - we SET all the bits as they are valid for pizzing
    //    we do this in ref_load_stored_reference AFTER all the SEC_REFERENCE/SEC_REF_IS_SEC sections are uncompressed,
    //    so that in case this was an REF_EXT_STORE compression, we first copy the contig-wide IS_SET sections (case 3) 
    //    (which will have 0s in the place of copied FASTA sections), and only after do we set these regions to 1.
    //    note: in case of PIZ, entire contigs are initialized to clear in ref_initialize_ranges as there might be
    //    regions missing (not covered by SEC_REFERENCE sections)
    // case 3: PIZ, reading a compacted reference - we receive the correct is_set in the SEC_REF_IS_SET section and don't change it


    // case: if compacted, this SEC_REF_IS_SET sections contains r->is_set and its first/last_pos contain the coordinates
    // of the range, while the following SEC_REFERENCE section contains only the bases for which is_set is 1, 
    // first_pos=0 and last_pos=(num_1_bits_in_is_set-1)
    if (is_compacted) {

        // if compacted, the section must be within the boundaries of the contig (this is not true if the section was copied with ref_copy_one_compressed_section)
        ASSERT (sec_start_within_contig >= 0 && ref_sec_last_pos <= r->last_pos, 
                "section range out of bounds for chrom=%d \"%s\": in SEC_REFERENCE being uncompressed: first_pos=%"PRId64" last_pos=%"PRId64" but in reference contig as initialized: first_pos=%"PRId64" last_pos=%"PRId64,
                chrom, r->chrom_name, ref_sec_pos, ref_sec_last_pos, r->first_pos, r->last_pos);

        ASSERT (uncomp_len == roundup_bits2bytes64 (ref_sec_len), "when uncompressing SEC_REF_IS_SET: uncomp_len=%u inconsistent with ref_sec_len=%"PRId64" (roundup_bits2bytes64 (ref_sec_len)=%"PRId64")", 
                uncomp_len, ref_sec_len, roundup_bits2bytes64 (ref_sec_len)); 

        // uncompress into r->is_set, via vb->compressed
        ASSERTNOTINUSE (vb->compressed);
        zfile_uncompress_section (vb, (SectionHeaderP)header, &vb->compressed, "compressed", 0, SEC_REF_IS_SET);

        BitArray *is_set = buf_zfile_buf_to_bitarray (&vb->compressed, ref_sec_len);

        // note on locking: while different threads uncompress regions of the range that are non-overlapping, 
        // there might be a 64b word that is split between two ranges
        RefLock lock = ref_lock (vb->ref, sec_start_gpos, ref_sec_len + 63); // +63 to ensure lock covers entire last word

        bit_array_copy (&r->is_set, sec_start_within_contig, is_set, 0, ref_sec_len); // initialization of is_set - case 3
        ref_unlock (vb->ref, lock);

        buf_free (&vb->compressed);

        // display contents of is_set if user so requested
        if (flag.show_is_set && !strcmp (chrom_name, flag.show_is_set)) 
            ref_print_is_set (r, -1, info_stream);

        // prepare for uncompressing the next section - which is the SEC_REFERENCE
        header = (SectionHeaderReference *)&vb->z_data.data[*ENT (uint32_t, vb->z_section_headers, 1)];

        if (flag.show_reference && primary_command == PIZ && r) 
            iprintf ("vb_i=%u Uncompressing %-14s chrom=%u ('%.*s') gpos=%"PRId64" pos=%"PRId64" num_bases=%u comp_bytes=%u\n", 
                     vb->vblock_i, st_name (header->h.section_type), BGEN32 (header->chrom_word_index), r->chrom_name_len, r->chrom_name, BGEN64 (header->gpos), 
                     BGEN64 (header->pos), BGEN32 (header->num_bases), BGEN32 (header->h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));

        compacted_ref_len  = (PosType)BGEN32(header->num_bases);
        uncomp_len         = BGEN32 (header->h.data_uncompressed_len);

        ASSERT (uncomp_len == roundup_bits2bytes64 (compacted_ref_len*2), 
                "uncomp_len=%u inconsistent with compacted_ref_len=%"PRId64, uncomp_len, compacted_ref_len); 

        ASSERT0 (BGEN32 (header->chrom_word_index) == chrom && BGEN64 (header->pos) == ref_sec_pos && BGEN64 (header->gpos) == ref_sec_gpos, // chrom should be the same between the two sections
                  "header mismatch between SEC_REF_IS_SET and SEC_REFERENCE sections");
    }
    
    // case: not compacted means that entire range is set
    else {
        ASSERT (uncomp_len == roundup_bits2bytes64 (ref_sec_len*2), "uncomp_len=%u inconsistent with ref_len=%"PRId64, uncomp_len, ref_sec_len); 

        if (primary_command == ZIP && flag.reference == REF_EXT_STORE) { // initialization of is_set - case 1
            RefLock lock = ref_lock (vb->ref, sec_start_gpos, ref_sec_len + 63); // +63 to ensure lock covers entire last word
            bit_array_clear_region (&r->is_set, sec_start_within_contig, ref_sec_len); // entire range is cleared
            ref_unlock (vb->ref, lock);
        }

        else if (primary_command == PIZ) { // initialization of is_set - case 2

            // it is possible that the section goes beyond the boundaries of the contig, this can happen when we compressed with --REFERENCE
            // and the section was copied in its entirety from the .ref.genozip file (in ref_copy_one_compressed_section)
            // even though a small amount of flanking regions is not set. in this case, we copy from the section only the part needed
            initial_flanking_len = (sec_start_within_contig < 0)    ? -sec_start_within_contig       : 0; // nucleotides in the section that are before the start of our contig
            final_flanking_len   = (ref_sec_last_pos > r->last_pos) ? ref_sec_last_pos - r->last_pos : 0; // nucleotides in the section that are after the end of our contig

            bit_index_t start = MAX_(sec_start_within_contig, 0);
            bit_index_t len   = ref_sec_len - initial_flanking_len - final_flanking_len;
            ASSERT (len >= 0 && len <= ref_sec_len, "expecting ref_sec_len=%"PRIu64" >= initial_flanking_len=%"PRIu64" + final_flanking_len=%"PRIu64,
                    ref_sec_len, initial_flanking_len, final_flanking_len);

            RefLock lock = ref_lock (vb->ref, start + r->gpos, len + 63); 
            bit_array_set_region (&r->is_set, start, len);
            ref_unlock (vb->ref, lock);

            // save the region we need to set, we will do the actual setting in ref_load_stored_reference
            spin_lock (vb->ref->region_to_set_list_spin);
            RegionToSet *rts = &NEXTENT (RegionToSet, vb->ref->region_to_set_list);
            spin_unlock (vb->ref->region_to_set_list_spin);
            rts->is_set    = &r->is_set;
            rts->first_bit = MAX_(sec_start_within_contig, 0);
            rts->len       = ref_sec_len - initial_flanking_len - final_flanking_len;
        }
   
        if (!uncomp_len) return;  // empty header - if it appears, it is the final header (eg in case of an unaligned SAM file)
    }

    // uncompress into r->ref, via vb->compressed
    ASSERTNOTINUSE (vb->compressed);
    zfile_uncompress_section (vb, (SectionHeaderP)header, &vb->compressed, "compressed", 0, SEC_REFERENCE);

    // lock - while different threads uncompress regions of the range that are non-overlapping, they might overlap at the bit level
    RefLock lock = ref_lock (vb->ref, sec_start_gpos, ref_sec_len + 63); // +63 to ensure lock covers entire last word

    if (is_compacted) {
        const BitArray *compacted = buf_zfile_buf_to_bitarray (&vb->compressed, compacted_ref_len * 2);
        ref_uncompact_ref (r, sec_start_within_contig, sec_end_within_contig, compacted);
    }

    else {
        BitArray *ref = buf_zfile_buf_to_bitarray (&vb->compressed, ref_sec_len * 2);

        // copy the section, excluding the flanking regions
        bit_array_copy (&r->ref, MAX_(sec_start_within_contig, 0) * 2, // dst
                        ref, initial_flanking_len * 2, // src
                        (ref_sec_len - initial_flanking_len - final_flanking_len) * 2); // len
    }

    ref_unlock (vb->ref, lock);

    buf_free (&vb->compressed);

finish:
    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. 

    COPY_TIMER (ref_uncompress_one_range);
}

static Reference ref_load_stored_reference_ref; // ref_read_one_range is run from the main thread, so no thread safetey issues
static Section sec = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
static void ref_read_one_range (VBlockP vb)
{
    START_TIMER;

    vb->ref = ref_load_stored_reference_ref;

    if (!sections_next_sec2 (&sec, SEC_REFERENCE, SEC_REF_IS_SET)) return;  // no more reference sections
        
    unsigned header_len = (Z_DT(DT_SAM) && z_file->z_flags.dts_ref_internal) ? crypt_padded_len (sizeof (SectionHeaderReference)) // REF_INTERNAL sections are encrypted if --password
                                                                             : sizeof (SectionHeaderReference); // data originating from a reference file is never encrypted
    if (((sec+1)->offset - sec->offset) == header_len) return; // final, terminating header-only section sometimes exists (see ref_compress_ref)
    
    if (sec->vblock_i == 0) // section was created with ref_copy_one_compressed_section
        z_file->num_copied_ref_sections++;
    else        
        ASSERT (sec->vblock_i + z_file->num_copied_ref_sections == vb->vblock_i, 
                "mismatch: sec->vblock_i=%u but vb->vblock_i=%u, z_file->num_copied_ref_sections=%u",
                sec->vblock_i, vb->vblock_i, z_file->num_copied_ref_sections);

    // if the user specified --regions, check if this ref range is needed
    bool range_is_included = true;
    RAEntry *ra = NULL;
    if (flag.regions) { 
        if (vb->vblock_i > vb->ref->stored_ra.len) return; // we're done - no more ranges to read, per random access (this is the empty section)

        ra = ENT (RAEntry, vb->ref->stored_ra, vb->vblock_i-1);
        ASSERT (ra->vblock_i == vb->vblock_i, "expecting ra->vblock_i(%u) == vb->vblock_i(%u)", ra->vblock_i, vb->vblock_i);

        range_is_included = regions_is_ra_included (ra);
    }

    if (range_is_included) { 
        buf_alloc (vb, &vb->z_section_headers, 0, 2, int32_t, 0, "z_section_headers"); // room for 2 section headers  
        ASSERT0 (vb->z_section_headers.len < 2, "unexpected 3rd recursive entry");

        NEXTENT (int32_t, vb->z_section_headers) = zfile_read_section (z_file, vb, sec->vblock_i, &vb->z_data, "z_data", sec->st, sec);
    }

    // if this is SEC_REF_IS_SET, read the SEC_REFERENCE section now (even if its not included - we need to advance the cursor)
    if (sec->st == SEC_REF_IS_SET) 
        ref_read_one_range (vb);

    if (flag.only_headers) 
        vb->z_data.len = 0; // roll back if we're only showing headers

    vb->ready_to_dispatch = true; // to simplify the code, we will dispatch the thread even if we skip the data, but we will return immediately. 

    COPY_TIMER (ref_read_one_range);
}

// PIZ: loading a reference stored in the genozip file - this could have been originally stored as REF_INTERNAL or REF_EXT_STORE
// or this could be a .ref.genozip file (called from load_external->piz_one_txt_file)
void ref_load_stored_reference (Reference ref)
{
    START_TIMER;

    ASSERTNOTINUSE (ref->ranges);

    if (!flag.only_headers)
        ref_initialize_ranges (ref, RT_LOADED);
    
    sec = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 

    spin_initialize (ref->region_to_set_list_spin);
    buf_alloc (evb, &ref->region_to_set_list, 0, sections_count_sections (SEC_REFERENCE), RegionToSet, 1, "region_to_set_list");
    
    // decompress reference using Dispatcher
    bool external = (bool)(flag.reference != REF_STORED);

    ref_load_stored_reference_ref = ref;
    dispatcher_fan_out_task ("load_ref", external ? ref->filename : z_file->basename, 
                             external ? PROGRESS_PERCENT : PROGRESS_NONE, 
                             external ? "Reading reference file" : NULL, 
                             flag.test, true, true, false, 0, 100,
                             ref_read_one_range, 
                             ref_uncompress_one_range, 
                             NULL);

    if (flag.only_headers) return;

    if (flag.show_ref_seq) ref_show_sequence (ref);

    // now we can safely set the is_set regions originating from non-compacted ranges. we couldn't do it before, because
    // copied-from-FASTA ranges appear first in the genozip file, and after them could be compacted ranges that originate
    // from a full-contig range in EXT_STORE, whose regions copied-from-FASTA are 0s.
    ARRAY (RegionToSet, rts, ref->region_to_set_list);
    for (uint32_t i=0; i < ref->region_to_set_list.len; i++)
        bit_array_set_region (rts[i].is_set, rts[i].first_bit, rts[i].len);

    COPY_TIMER_VB (evb, ref_load_stored_reference);
}

// ---------------------
// Cache stuff
// ---------------------

static inline const char *ref_get_cache_fn (Reference ref)
{
    if (!ref->cache_fn) {
        ref->cache_fn = MALLOC (strlen (z_name) + 20);
        sprintf ((char *)ref->cache_fn, "%s.gcache", z_name); // constant thereafter
    }

    return ref->cache_fn;
}

void ref_remove_cache (Reference ref)
{
    file_remove (ref_get_cache_fn(ref), true);
}

// mmap the reference cached file, as copy-on-write - modifications are private to process and not written to the file
bool ref_mmap_cached_reference (Reference ref)
{  
    if (!buf_is_alloc (&ref->ranges)) {  // possibly already loaded from previous file
        if (!file_exists (ref_get_cache_fn(ref))) return false; // file doesn't exist

        ref_initialize_ranges (ref, RT_CACHED); // also does the actual buf_mmap
    }

    if (ref_has_is_set()) buf_zero (&ref->genome_is_set_buf);

    // PIZ: all ranges of contigs are "set", i.e. the genome is valid in this location (it is not set in the short space between contigs)
    if (primary_command == PIZ) 
        for (uint32_t chrom=0; chrom < ref->ranges.len; chrom++) {
            Range *r = ENT (Range, ref->ranges, chrom);
            bit_array_set_region (&r->is_set, 0, ref_size (r));
        }

    if (flag.show_ref_seq) ref_show_sequence (ref);

    return true;
}

static void ref_create_cache (VBlockP cache_create_vb) 
{
    buf_dump_to_file (ref_get_cache_fn(cache_create_vb->ref), &cache_create_vb->ref->genome_cache, 1, true, false, false);
}

// initiate creating the the genome cache in a background thread
void ref_create_cache_in_background (Reference ref)
{
    if (flag.regions) return; // we can't create cache as we haven't loaded the entire reference

    ref->cache_create_vb = vb_initialize_nonpool_vb (VB_ID_GCACHE_CREATE, DT_NONE, "ref_create_cache_in_background");
    ref->cache_create_vb->ref = ref;
    ref->cache_create_vb->compute_task = "create_reference_cache";

    ref_get_cache_fn (ref); // generate name before closing z_file
    ref->ref_cache_creation_thread_id = threads_create (ref_create_cache, ref->cache_create_vb);
}

void ref_create_cache_join (Reference ref, bool free_mem)
{
    if (!ref->cache_create_vb) return;

    threads_join (&ref->ref_cache_creation_thread_id, NULL);

    if (free_mem) vb_destroy_vb (&ref->cache_create_vb);
}


// ------------------------------------
// ZIP side
// ------------------------------------

// case REF_INTERNAL - we get a range_id by hashing the chrom name and the range_i - we can't use the chrom_node_index
// because we have multiple threads in parallel that might have the same node_index for different chroms
static inline uint32_t ref_range_id_by_hash (VBlockP vb, uint32_t range_i)
{
    ASSERT0 (vb->chrom_name_len > 0, "vb->chrom_name_len==0");

    uint32_t value, n=0;
    bool is_major_chrom=false;

    // step 1: get number embedded in the chrom name
    for (unsigned i=0; i < vb->chrom_name_len; i++)
        if (IS_DIGIT (vb->chrom_name[i])) 
            n = n*10 + (vb->chrom_name[i] - '0');

    // short name - possibly a major chromosome
    if (vb->chrom_name_len <= 5) {
        is_major_chrom=true; // possible major...

        // chromosome name contains a number of 1-124 and chromosome name is short - major chromosomes
        if (n >= 1 && n <= 124) {
            // we're good
        }
        // other major chromosomes - a number 125-127 (note: these are for human, we can add here others for other popular species)
        #define IS_CHROM(s) (vb->chrom_name_len == strlen(s) && !memcmp (vb->chrom_name, s, strlen(s))) // hopefully the compiler optimizes away strlen(const s) 
        else if (IS_CHROM ("X") || IS_CHROM ("chrX")) n = 125;
        else if (IS_CHROM ("Y") || IS_CHROM ("chrY")) n = 126;
        else if (IS_CHROM ("M") || IS_CHROM ("chrM") || IS_CHROM ("MT") || IS_CHROM ("chrMT")) n = 127; // note: even though MT is short, it is major, as we might have many reads for it and we don't want hash contention

        else is_major_chrom = false; // not major
    }

    // non-major chromosomse - if n is too small (perhaps indicating non-uniqueness of the number) 
    // get a new n (a number 0->28668) derived from the last 8 characters of the chromosome name
    if (!is_major_chrom && n < 10000) 
        n = (uint32_t)( ((                          ((uint64_t)vb->chrom_name[vb->chrom_name_len-1]) << 0)       |   
                        (vb->chrom_name_len >= 2 ? (((uint64_t)vb->chrom_name[vb->chrom_name_len-2]) << 8)  : 0) |
                        (vb->chrom_name_len >= 3 ? (((uint64_t)vb->chrom_name[vb->chrom_name_len-3]) << 16) : 0) | 
                        (vb->chrom_name_len >= 4 ? (((uint64_t)vb->chrom_name[vb->chrom_name_len-4]) << 24) : 0) |
                        (vb->chrom_name_len >= 5 ? (((uint64_t)vb->chrom_name[vb->chrom_name_len-5]) << 32) : 0) |
                        (vb->chrom_name_len >= 6 ? (((uint64_t)vb->chrom_name[vb->chrom_name_len-6]) << 40) : 0) |
                        (vb->chrom_name_len >= 7 ? (((uint64_t)vb->chrom_name[vb->chrom_name_len-7]) << 48) : 0) |
                        (vb->chrom_name_len >= 8 ? (((uint64_t)vb->chrom_name[vb->chrom_name_len-8]) << 56) : 0) ) % 28669); // 28669 is prime -> even distribution

    // a non-major chromosome with a number that appears unique - just mod it into 0-28668
    else if (!is_major_chrom)
        n %= 28669;

    // major chromosomes - we have up to 127 of them, each can have 1024 ranges of 1Mbp each
    if (is_major_chrom)
        value = (0b111 << 17)  |   // top 3 MSb 111
                (n     << 10)  |   // chromosome component - 7 bits (values 1-127)
                (range_i & 0x3ff); // range_i component - 10 bits (up to 1024 ranges of 1 Mbp each =REF_NUM_DENOVO_SITES_PER_RANGE)
    
    // non-major chromosomes - we have 28669 slots for for them - fit into 15 bit with top 3 MSB being 000->110 
    // (as 111 is major chromosomes) each can have 32 ranges of 1Mbp each
    else
        value = (n << 5) |         // chromosome component - 15 bits, top 3 MSB are NOT 111 (values 0 - 28668)
                (range_i & 0x1f);  // range_i component - 5 bits (up to 32 ranges of 1 Mbp each)

    return value; 
}

static Range *ref_seg_get_locked_range_denovo (VBlockP vb, Reference ref, WordIndex chrom, PosType pos, const char *field /* used for ASSSEG */, RefLock *lock)  
{
    uint32_t range_i = pos2range_i (pos); // range within contig 

    // case: we're asking for the same range as the previous one (for example, subsequent line in a sorted SAM)
    if (vb && vb->prev_range[0] && vb->prev_range_chrom_node_index[0] == chrom && vb->prev_range_range_i == range_i) {
        *lock = ref_lock_range (ref, ENTNUM (ref->ranges, vb->prev_range[0]));
        return vb->prev_range[0];
    }

    uint32_t range_id = ref_range_id_by_hash (vb, range_i);
    ASSSEG (range_id < ref->ranges.len, field, "range_id=%u expected to be smaller than ranges.len=%u", range_id, (uint32_t)ref->ranges.len);

    Range *range = ENT (Range, ref->ranges, range_id);
    *lock = ref_lock_range (ref, range_id);

    if (range->ref.nbits) {

        // check for hash conflict 
        if (range->range_i != range_i || !str_issame (vb->chrom_name, range->chrom_name)) {
            *lock = ref_unlock (ref, *lock);

            ASSERTW (!flag.seg_only, "Warning: ref range contention: chrom=%.*s pos=%u (this slightly affects compression ratio, but is harmless)", 
                     vb->chrom_name_len, vb->chrom_name, (uint32_t)pos); // only show this in --seg-only

            range = NULL;  // no soup for you
        }

        return range; // range already initialized or cannot be initialized ^ - we're done
    }

    *range = (Range){
        .range_id       = range_id,
        .range_i        = range_i,
        .first_pos      = range_i2pos (range_i),
        .last_pos       = range_i2pos (range_i) + REF_NUM_DENOVO_SITES_PER_RANGE - 1,
        .chrom_name_len = vb->chrom_name_len,
        .chrom          = WORD_INDEX_NONE,  // Note: in REF_INTERNAL the chrom index is private to the VB prior to merge, so we can't use it
        .chrom_name     = MALLOC (vb->chrom_name_len),
        .ref            = bit_array_alloc (REF_NUM_DENOVO_SITES_PER_RANGE * 2, false),
        .is_set         = bit_array_alloc (REF_NUM_DENOVO_SITES_PER_RANGE, true), // nothing is set yet - in ZIP, bits get set as they are encountered in the compressing txt file
    };

    // vb->chrom_name points into vb->txt_data that will disappear at the end of this VB, 
    memcpy ((char *)range->chrom_name, vb->chrom_name, vb->chrom_name_len);

    if (vb) {
        vb->prev_range[0] = range;
        vb->prev_range_chrom_node_index[0] = chrom;
        vb->prev_range_range_i = range_i;
    }

    return range; // returning locked range
}

// Seg of VCF and SAM - compressing against a reference when chrom is known
static Range *ref_seg_get_locked_range_loaded (VBlockP vb, Reference ref, WordIndex chrom, STRp(chrom_name), 
                                               PosType pos, uint32_t seq_len, 
                                               WordIndex ref_index, // if known (mandatory if not prim_chrom), WORD_INDEX_NONE if not
                                               const char *field /* used for ASSSEG */, 
                                               RefLock *lock) // optional - range locked if provided
{
    if (ref_index == WORD_INDEX_NONE)
        ref_index = chrom_2ref_seg_get (ref, vb, chrom); // only works for prim_chrom

    if (ref_index == WORD_INDEX_NONE) return NULL;

    Range *range = ENT (Range, ref->ranges, ref_index);

    if (lock) {
        PosType gpos = range->gpos + (pos - range->first_pos);
        *lock = ref_lock (ref, gpos, seq_len);
    }

    // when using an external refernce, pos has to be within the reference range
    // note: in SAM, if a read starts within the valid range, it is allowed to overflow beyond it - and we will circle
    // around to the beginning of the range assuming its a circular chromosome (see in sam_seg_seq_field)
    ASSSEG ((pos >= range->first_pos && pos <= range->last_pos), field,
            "POS=%"PRId64" for contig \"%.*s\" (mapped to reference contig \"%.*s\" ref_index=%d), but this contig's range is %"PRId64" - %"PRId64". Likely this is because %s was created using a reference file other than %s.",
            pos, vb->chrom_name_len, vb->chrom_name, range->chrom_name_len, range->chrom_name, ref_index, range->first_pos, range->last_pos, txt_name, ref->filename);

    vb->prev_range[ref==prim_ref] = range;
    vb->prev_range_chrom_node_index[ref==prim_ref] = chrom; // the chrom that started this search, leading to this range

    return range; 
}

// ZIP: returns a range that includes pos and is locked for (pos,len)
// case 1: ZIP: in SAM with REF_INTERNAL, when segging a SEQ field ahead of committing it to the reference
// case 2: ZIP: SAM and VCF with REF_EXTERNAL: when segging a SAM_SEQ or VCF_REFALT field
// if range is found, returns a locked range, and its the responsibility of the caller to unlock it. otherwise, returns NULL
Range *ref_seg_get_locked_range (VBlockP vb, Reference ref, WordIndex chrom, STRp(chrom_name), 
                                 PosType pos, uint32_t seq_len, 
                                 WordIndex ref_index, // if known (mandatory if not prim_chrom), WORD_INDEX_NONE if not                                
                                 const char *field /* used for ASSSEG */, 
                                 RefLock *lock) // optional if RT_LOADED/RT_CACHED
{
    // sanity checks
    ASSERT0 (vb->chrom_name, "vb->chrom_name=NULL");

    switch (ranges_type(ref)) {
        case RT_DENOVO : return ref_seg_get_locked_range_denovo (vb, ref, chrom, pos, field, lock);
        case RT_CACHED :
        case RT_LOADED : return ref_seg_get_locked_range_loaded (vb, ref, chrom, STRa(chrom_name), pos, seq_len, ref_index, field, lock);
        default        : ABORT ("Error in ref_seg_get_locked_range: invalid ranges_type=%"PRId64, ranges_type(ref)); return 0;
    }
}

// ----------------------------------------------
// Compressing ranges into SEC_REFERENCE sections
// ----------------------------------------------

// ZIP main thread
static void ref_copy_one_compressed_section (Reference ref, File *ref_file, const RAEntry *ra, Section *sl)
{
    // get section list entry from ref_file_section_list - which will be used by zfile_read_section to seek to the correct offset
    while (*sl < AFTERENT (SectionEnt, ref->ref_file_section_list) && 
           !((*sl)->vblock_i == ra->vblock_i && (*sl)->st == SEC_REFERENCE)) 
        (*sl)++;

    ASSERT (*sl < AFTERENT (SectionEnt, ref->ref_file_section_list), "cannot find FASTA_NONREF of vb_i=%u in section list of reference file", ra->vblock_i);

    static Buffer ref_seq_section = EMPTY_BUFFER;

    CLEAR_FLAG (show_headers);
    zfile_read_section (ref_file, evb, ra->vblock_i, &ref_seq_section, "ref_seq_section", SEC_REFERENCE, *sl);
    RESTORE_FLAG (show_headers);

    SectionHeaderReference *header = (SectionHeaderReference *)ref_seq_section.data;

    WordIndex ref_index = BGEN32 (header->chrom_word_index);
    ASSERT0 (ref_index == ra->chrom_index && BGEN64 (header->pos) == ra->min_pos,
              "RA and Section don't agree on chrom or pos");

    // update header->chrom_word_index, currently refering to ref->contigs, to refer to ZCTX(CHROM)
    STR(contig_name);
    contig_name = ref_contigs_get_name (ref, ref_index, &contig_name_len);

    WordIndex chrom_word_index = chrom_get_by_name (STRa(contig_name));
    ASSERT (chrom_word_index != WORD_INDEX_NONE, "Cannot find reference contig \"%.*s\" in CHROM context", STRf(contig_name));

    header->chrom_word_index = BGEN32 (chrom_word_index);

    // some minor changes to the header...
    header->h.vblock_i  = 0; // we don't belong to any VB and there is no encryption of external ref

    // "manually" add the reference section to the section list - normally it is added in comp_compress()
    sections_add_to_list (evb, &header->h);
    sections_list_concat (evb); // must be called before disk_so_far is updated

    // Write header and body of the reference to z_file
    // Note on encryption: reference sections originating from an external reference are never encrypted - not
    // by us here, and not in the source reference fasta (because with disallow --make-reference in combination with --password)
    START_TIMER;
    file_write (z_file, ref_seq_section.data, ref_seq_section.len);
    COPY_TIMER_VB (evb, write);

    z_file->disk_so_far += ref_seq_section.len;   // length of GENOZIP data writen to disk

    if (flag.show_reference) {
        Context *ctx = ZCTX(CHROM);
        CtxNode *node = ENT (CtxNode, ctx->nodes, BGEN32 (header->chrom_word_index));
        iprintf ("Copying SEC_REFERENCE from %s: chrom=%u (%s) gpos=%"PRId64" pos=%"PRId64" num_bases=%u section_size=%u\n", 
                 ref->filename, BGEN32 (header->chrom_word_index), 
                 ENT (char, ctx->dict, node->char_index), 
                 BGEN64 (header->gpos), BGEN64 (header->pos), 
                 BGEN32 (header->num_bases), 
                 BGEN32 (header->h.data_compressed_len) + BGEN32 (header->h.compressed_offset));
    }

    buf_free (&ref_seq_section);
}

// ZIP copying parts of external reference to fine - called by main thread from zip_write_global_area->ref_compress_ref
static void ref_copy_compressed_sections_from_reference_file (Reference ref)
{
    ASSERT (primary_command == ZIP && flag.reference == REF_EXT_STORE, 
            "not expecting to be here: primary_command=%u flag.reference=%u", primary_command, flag.reference);

    File *ref_file = file_open (ref->filename, READ, Z_FILE, DT_FASTA);

    // note: in a FASTA file compressed with --make-reference, there is exactly one RA per VB (a contig or part of a contig)
    // and, since this is ZIP with EXT_STORE, also exactly one range per contig. We loop one RA at a time and:
    // 1. If 95% of the ref file RA is set in the zfile contig range - we copy the compressed reference section directly from the ref FASTA
    // 2. If we copied from the FASTA, we mark those region covered by the RA as "is_set=0", so that we don't compress it later
    Section sl = FIRSTENT (SectionEnt, ref->ref_file_section_list);
    ARRAY (RAEntry, fasta_sec, ref->ref_external_ra);

    chrom_index_by_name (CHROM);

    for (uint32_t i=0; i < ref->ref_external_ra.len; i++) {

        Range *contig_r = ENT (Range, ref->ranges, fasta_sec[i].chrom_index);
        PosType fasta_sec_start_in_contig_r = fasta_sec[i].min_pos - contig_r->first_pos; // the start of the FASTA section (a bit less than 1MB) within the full-contig range

        PosType fasta_sec_len = fasta_sec[i].max_pos - fasta_sec[i].min_pos + 1;
        PosType bits_is_set   = bit_array_num_bits_set_region (&contig_r->is_set, fasta_sec_start_in_contig_r, fasta_sec_len);

        // if this at least 95% of the RA is covered, just copy the corresponding FASTA section to our file, and
        // mark all the ranges as is_set=false indicating that they don't need to be compressed individually
        if ((double)bits_is_set / (double)fasta_sec_len >= 0.95) {
            ref_copy_one_compressed_section (ref, ref_file, &fasta_sec[i], &sl);
            bit_array_clear_region (&contig_r->is_set, fasta_sec_start_in_contig_r, fasta_sec_len);
        }
    }

    file_close (&ref_file, false, false);
}

// remove the unused parts of a range and the beginning and end of the range, and update first/last_pos.
static bool ref_remove_flanking_regions (Reference ref, Range *r, uint64_t *start_flanking_region_len /* out */)
{
// threshold - if the number of clear bits (excluding flanking regions) is beneath this, we will not copmact, as the cost in
// z_file size of a SEC_REF_IS_SET section needed if compacting will be more than what we save in compression of r->ref
#define THRESHOLD_FOR_COMPACTING 470 

    uint64_t end_flanking_region_len, last_1;
    
    char has_any_bit = bit_array_find_first_set_bit (&r->is_set, start_flanking_region_len);
    ASSERT (has_any_bit, "range %u (%s) has no bits set in r->is_set", ENTNUM (ref->ranges, r), r->chrom_name); // ref_prepare_range_for_compress is responsible not to send us 0-bit ranges

    has_any_bit = bit_array_find_prev_set_bit (&r->is_set, r->is_set.nbits, &last_1);
    ASSERT (has_any_bit, "range %u (%s) has no bits set in r->is_set (#2)", ENTNUM (ref->ranges, r), r->chrom_name); // this should definitely never happen, since we already know the range has bits
    end_flanking_region_len = r->is_set.nbits - last_1 - 1;

    uint64_t num_clear_bits_excluding_flanking_regions = 
        r->is_set.nbits - r->num_set - *start_flanking_region_len - end_flanking_region_len;

    // remove flanking regions - will allow a smaller allocation for the reference in PIZ 
    r->gpos      += *start_flanking_region_len;
    r->first_pos += *start_flanking_region_len;
    r->last_pos  -= end_flanking_region_len;

    ASSERT (r->last_pos >= r->first_pos, "bad removal of flanking regions: r->first_pos=%"PRId64" r->last_pos=%"PRId64,
            r->first_pos, r->last_pos);

    bit_array_remove_flanking (&r->is_set, *start_flanking_region_len, end_flanking_region_len);

    // if all we're doing is removing the flanking regions, we update ref now. if not, ref_compact_ref will update it
    bool is_compact_needed = num_clear_bits_excluding_flanking_regions >= THRESHOLD_FOR_COMPACTING;
    if (!is_compact_needed) 
        bit_array_remove_flanking (&r->ref, *start_flanking_region_len * 2, end_flanking_region_len * 2);

    // return true if compacting is needed
    return is_compact_needed;
}

// we compact one range by squeezing together all the bases that have is_set=1. return true if compacted
static bool ref_compact_ref (Reference ref, Range *r)
{
    if (!r || !r->num_set) return false;

    ASSERT0 (r->is_set.nbits, "r->is_set.nbits=0");

    // remove flanking regions
    uint64_t start_flanking_region_len;
    bool is_compact_needed = ref_remove_flanking_regions (ref, r, &start_flanking_region_len);
    if (!is_compact_needed) return false;

    uint64_t start_1_offset=0, start_0_offset=0, compact_len=0;
    while (1) {
        
        // find length of set region
        bool has_any_bit = bit_array_find_next_clear_bit (&r->is_set, start_1_offset, &start_0_offset);
        uint64_t len_1 = (has_any_bit ? start_0_offset : r->is_set.nbits) - start_1_offset;

        // do actual compacting - move set region to be directly after the previous set region (or at the begining if its the first)
        bit_array_copy (&r->ref, compact_len * 2, &r->ref, (start_flanking_region_len + start_1_offset) * 2, len_1 * 2);
        compact_len += len_1;

        if (!has_any_bit) break; // case: we're done- this 1 region goes to the end the range - there are no more clear regions

        // find length of clear region
        has_any_bit = bit_array_find_next_set_bit (&r->is_set, start_0_offset, &start_1_offset);
        ASSERT0 (has_any_bit, "cannot find set bits"); // this should never happen as we removed the flanking regions
    }

    // set length of ref - this is the data that will be compressed
    r->ref.nbits  = compact_len * 2;
    r->ref.nwords = roundup_bits2words64 (r->ref.nbits); 

    return true;
}

static void ref_compress_one_range (VBlockP vb)
{
    Range *r = vb->range; // will be NULL if we're being asked to write a final, empty section

    // remove flanking regions, and if beneficial also compact it further by removing unused nucleotides 
    bool is_compacted = flag.make_reference ? false : ref_compact_ref (gref, r); // true if it is compacted beyong just the flanking regions

    // get the index into the ZCTX(CHROM) dictionary
    WordIndex chrom_word_index = !r                                   ? WORD_INDEX_NONE
                               : (flag.reference & REF_ZIP_CHROM2REF) ? *ENT (WordIndex, z_file->ref2chrom_map, r->chrom)
                               :                                        r->chrom;

    SectionHeaderReference header = { .h.vblock_i          = BGEN32 (vb->vblock_i),
                                      .h.magic             = BGEN32 (GENOZIP_MAGIC),
                                      .h.compressed_offset = BGEN32 (sizeof(header)),
                                      .chrom_word_index    = BGEN32 (chrom_word_index),
                                      .pos                 = r ? BGEN64 (r->first_pos) : 0,
                                      .gpos                = r ? BGEN64 (r->gpos)      : 0 };

    vb->z_data.param = vb->vblock_i;

    // First, SEC_REF_IS_SET section (but not needed if the entire range is used)
    if (r && is_compacted) {

        LTEN_bit_array (&r->is_set);

        header.h.section_type          = SEC_REF_IS_SET;  // most of the header is the same as ^
        header.h.codec                 = CODEC_BZ2;
        header.h.data_uncompressed_len = BGEN32 (r->is_set.nwords * sizeof (uint64_t));
        header.num_bases               = BGEN32 ((uint32_t)ref_size (r)); // full length, after flanking regions removed
        comp_compress (vb, &vb->z_data, (SectionHeader*)&header, (char *)r->is_set.words, NULL);

        if (flag.show_reference && r) 
            iprintf ("vb_i=%u Compressing SEC_REF_IS_SET chrom=%u (%.*s) gpos=%"PRIu64" pos=%"PRIu64" num_bases=%u section_size=%u bytes\n", 
                     vb->vblock_i, BGEN32 (header.chrom_word_index), r->chrom_name_len, r->chrom_name,
                     BGEN64 (header.gpos), BGEN64 (header.pos), BGEN32 (header.num_bases), 
                     BGEN32 (header.h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));
    }

    // Second, SEC_REFERENCE
    if (r) LTEN_bit_array (&r->ref);

    header.h.section_type          = SEC_REFERENCE;
    header.h.codec                 = CODEC_LZMA; // better than BSC: slightly better compression and compression speed, 2.5X faster decompression
    header.h.compressed_offset     = BGEN32 (sizeof(header)); // reset compressed offset - if we're encrypting - REF_IS_SET was encrypted and compressed_offset padded, by REFERENCE is never encrypted
    header.h.data_uncompressed_len = r ? BGEN32 (r->ref.nwords * sizeof (uint64_t)) : 0;
    header.num_bases               = r ? BGEN32 (r->ref.nbits / 2) : 0; // less than ref_size(r) if compacted
    comp_compress (vb, &vb->z_data, (SectionHeader*)&header, r ? (char *)r->ref.words : NULL, NULL);

    if (flag.show_reference && r) 
        iprintf ("vb_i=%u Compressing SEC_REFERENCE chrom=%u (%.*s) %s gpos=%"PRIu64" pos=%"PRIu64" num_bases=%u section_size=%u bytes\n", 
                 vb->vblock_i, BGEN32 (header.chrom_word_index), r->chrom_name_len, r->chrom_name, is_compacted ? "compacted " : "",
                 BGEN64 (header.gpos), BGEN64 (header.pos), BGEN32 (header.num_bases), 
                 BGEN32 (header.h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));

    // store the stored_ra data for this range
    if (r) {
        spin_lock (gref->stored_ra_spin);
        NEXTENT (RAEntry, gref->stored_ra) = (RAEntry){ .vblock_i    = vb->vblock_i, 
                                                        .chrom_index = r->chrom,
                                                        .min_pos     = r->first_pos,
                                                        .max_pos     = r->last_pos };
        spin_unlock (gref->stored_ra_spin);
    }

    // insert this range sequence into the ref_hash (included in the reference file, for use to compress of FASTQ, unaligned SAM and FASTA)
    if (flag.make_reference)
        refhash_calc_one_range (r, ISLASTENT (gref->ranges, r) ? NULL : r+1);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

#define MAX_BASES_FOR_RANGE_MERGE REF_NUM_DENOVO_SITES_PER_RANGE
#define MIN_BASES_TO_NOT_MERGE (MAX_BASES_FOR_RANGE_MERGE/2)
#define MAX_MERGED_RANGES 32 // if this number is too high, then this function becomes a bottleneck and we can't scale threads

// look forward and count how many additional ranges are going to get merged 
static inline unsigned ref_prepare_expected_more_merges (const Range *this_r, int64_t num_set_this_merge, uint32_t max_ranges)
{
    unsigned more=0;
    for (Range *r=(Range *)this_r+1; r < AFTERENT (Range, gref->ranges) && more <= max_ranges; r++, more++) {

        if (r->num_set == -1) // calcualte num_set if not already calculated
            r->num_set = bit_array_num_bits_set (&r->is_set);

        num_set_this_merge += r->num_set;

        // cases where more < max_ranges:
        if (r->chrom != this_r->chrom ||          // next range is on a different chrom
            r->first_pos != (r-1)->last_pos+1 ||  // next range is not consecutive pos
            num_set_this_merge > MAX_BASES_FOR_RANGE_MERGE) // we reached the maximum number of set bases
            break;
    }

    return MIN_(more, max_ranges);
}

// compress the reference - one section at the time, using Dispatcher to do them in parallel 
// note: this is not called in make_reference - instead, ref_make_prepare_range_for_compress is called
// note: this function is run in the main thread, so no thread safety issues with static variables
static void ref_prepare_range_for_compress (VBlockP vb)
{
    static uint32_t next_range_i=0;
    unsigned num_merged=0;

    // initialize if new file
    if (vb->vblock_i == 1) next_range_i = 0;

    // find next occupied range
    for (; next_range_i < gref->ranges.len ; next_range_i++) {
        Range *r = ENT (Range, gref->ranges, next_range_i);

        // case: we're done merging, as this range has a different chrom or not consecutive, therefore we're not consuming ranges[next_range_i]
        if (vb->range && (num_merged >= MAX_MERGED_RANGES || r->chrom != vb->range->chrom || r->first_pos != vb->range->last_pos+1)) break; 

        if (r->num_set == -1) // calcualte num_set if not already calculated
            r->num_set = bit_array_num_bits_set (&r->is_set);

        if (!r->num_set) {
            r->is_set.nbits = 0;
            continue; // nothing to with this range - perhaps copied and cleared in ref_copy_compressed_sections_from_reference_file
        }

        // case: we're done merging - we have all the bases we need already, therefore we're not consuming ranges[next_range_i]
        if (vb->range && (vb->range->num_set + r->num_set > MAX_BASES_FOR_RANGE_MERGE))
            break; 

        // case: first range in merge set
        else if (!vb->range) {
            vb->range              = r; // range to compress
            vb->ready_to_dispatch  = true;

            // case: we have enough bases, no need to attempt merging
            if (r->num_set > MIN_BASES_TO_NOT_MERGE) {
                next_range_i++; // next call to this function we move on to the next range
                break; 
            }
        }

        // case: merge r into vb->range 
        else {
            num_merged++;

            vb->range->num_set += r->num_set;
            vb->range->last_pos = r->last_pos;

            static unsigned expected_more_merges=0;
            expected_more_merges = (num_merged==1) ? ref_prepare_expected_more_merges (r, vb->range->num_set, MAX_MERGED_RANGES-1)
                                                   : expected_more_merges - 1;

            bit_array_concat (&vb->range->ref, &r->ref, expected_more_merges);            
            bit_array_concat (&vb->range->is_set, &r->is_set, expected_more_merges);
            bit_array_free (&r->ref);
            bit_array_free (&r->is_set);            
        }
    }
  }

static int ref_contigs_range_sorter (const void *a, const void *b)
{
    const Range *r_a = a;
    const Range *r_b = b;
    
    bool is_used_a = ref_is_range_used (r_a);
    bool is_used_b = ref_is_range_used (r_b);
    
    // if one of them is not initialized - it comes after
    if (!is_used_a && !is_used_b) return 0;
    if (!is_used_a) return 1;
    if (!is_used_b) return -1;

    if (r_a->chrom != r_b->chrom)
        return r_a->chrom - r_b->chrom; // ascending order
    
    return r_a->last_pos - r_b->last_pos; // ascending order
}

// sorts denovo ranges by chrom and pos, removes unused ranges from ranges, 
static void ref_finalize_denovo_ranges (void)
{
    // create contigs from CHROM dictionary
    ref_contigs_generate_data_if_denovo (gref);

    // calculate all chrom indices
    for (uint32_t range_i=0; range_i < gref->ranges.len; range_i++) {
        Range *r = ENT (Range, gref->ranges, range_i);
        
        if (ref_is_range_used (r))
            r->chrom = ref_contigs_get_by_name (gref, r->chrom_name, r->chrom_name_len, false, false);
    }

    // sort by chrom then pos, and place the unused ranges at the end
    qsort (gref->ranges.data, gref->ranges.len, sizeof (Range), ref_contigs_range_sorter);

    // shorten the array to only used ranges
    Range *r ; for (r=FIRSTENT (Range, gref->ranges); r < AFTERENT (Range, gref->ranges) && ref_is_range_used (r); r++) {};
    gref->ranges.len = ENTNUM (gref->ranges, r);

    // shorten the last range in each to contig according to usage
    for (uint32_t range_i=0; range_i < gref->ranges.len; range_i++) {
        Range *r = ENT (Range, gref->ranges, range_i);
        
        if (range_i < gref->ranges.len - 1 && r->chrom == (r+1)->chrom) continue; // not last range of the contig

        bit_index_t effective_len = bit_array_effective_length (&r->is_set);
        r->last_pos = r->first_pos + effective_len - 1;
        bit_array_truncate (&r->ref, 2 * effective_len);
        bit_array_truncate (&r->is_set, effective_len);
    }
}

// ZIP: compress and write reference sections. either compressed by us, or copied directly from the reference file.
void ref_compress_ref (void)
{
    if (!buf_is_alloc (&gref->ranges)) return;

    ref_create_cache_join (gref, true); // finish dumping reference to cache before we modify it via compacting

    // calculate z_file->ref2chrom_map, the inverse of z_file->chrom2ref_map
    if (flag.reference & REF_ZIP_CHROM2REF)
        chrom_calculate_ref2chrom (ref_num_contigs(gref));

    if ((ranges_type(gref) == RT_DENOVO) &&
        buf_is_alloc (&ZCTX(CHROM)->dict)) // did we have an aligned lines? (to do: this test is not enough)
        ref_finalize_denovo_ranges(); // assignes chroms; sorts ranges by chrom, pos; gets rid of unused ranges

    if (ranges_type(gref) == RT_DENOVO)
        ref_contigs_compress_internal (gref); // also assigns gpos to de-novo ranges 

    else if (ranges_type(gref) != RT_MAKE_REF)
        ref_contigs_compress_ext_store (gref);  

    // copy already-compressed SEC_REFERENCE sections from the genozip reference file, but only such sections that are almost entirely
    // covered by ranges with is_accessed=true. we mark these ranges affected as is_accessed=false.
    if (ranges_type(gref) == RT_LOADED || ranges_type(gref) == RT_CACHED)
        ref_copy_compressed_sections_from_reference_file (gref);

    buf_alloc (evb, &gref->stored_ra, 0, gref->ranges.len, RAEntry, 1, "stored_ra");
    gref->stored_ra.len = 0; // re-initialize, in case we read the external reference into here
    
    spin_initialize (gref->stored_ra_spin);

    // compression of reference doesn't output % progress
    SAVE_FLAGS;
    if (flag.show_reference) flag.quiet = true; // show references instead of progress

    // initialize Range.num_set (used by ref_prepare_range_for_compress)
    if (!flag.make_reference)
        for (unsigned range_i=0; range_i < gref->ranges.len; range_i++)
            if (! ENT (Range, gref->ranges, range_i)->num_set)  // if not already set by ref_contigs_populate_aligned_chroms
                ENT (Range, gref->ranges, range_i)->num_set = -1;

    // proceed to compress all ranges that have still have data in them after copying
    Dispatcher dispatcher = 
        dispatcher_fan_out_task ("compress_ref", NULL, PROGRESS_MESSAGE, "Writing reference...", false, true, true, false, 0, 5000,
                                 flag.make_reference ? ref_make_prepare_range_for_compress : ref_prepare_range_for_compress, 
                                 ref_compress_one_range, 
                                 zfile_output_processed_vb);

    uint32_t num_vbs_dispatched = dispatcher_get_next_vb_i (dispatcher);

    RESTORE_FLAGS;
    
    // SAM require at least one reference section, but if the SAM is unaligned, there will be none - create one empty section
    // (this will also happen if SAM has just only reference section, we will just needlessly write another tiny section - no harm)
    // incidentally, this empty section will also be written in case of a small (one vb) reference - no harm
    if (Z_DT(DT_SAM) && num_vbs_dispatched==1) {
        evb->range = NULL;
        ref_compress_one_range (evb); // written with vb_i=0, section header only (no body)
    }

    // compress reference random access (note: in case of a reference file, SEC_REF_RAND_ACC will be identical to SEC_RANDOM_ACCESS. That's ok, its tiny)
    random_access_finalize_entries (&gref->stored_ra); // sort in order of vb_i

    // range data needed for contigs is only set in ref_make_prepare_range_for_compress which happens as part of the dispatcher_fan_out_task
    if (ranges_type (gref) == RT_MAKE_REF) {
        // sort by chrom then pos - so later piz can binary-search by chrom index
        qsort (gref->ranges.data, gref->ranges.len, sizeof (Range), ref_contigs_range_sorter);

        ref_contigs_compress_internal (gref); 
    }
}

// -------------------------------
// Loading an external reference
// -------------------------------

void ref_set_reference (Reference ref, const char *filename, ReferenceType ref_type, bool is_explicit)
{
    if (!is_explicit && ref->filename) return; // already set explicitly

    unsigned filename_len;
    if (!filename) {
        const char *env = getenv ("GENOZIP_REFERENCE");
        if (!env) return; // nothing to set
        str_split (env, strlen (env), 2, ':', ref_fn, false)
        ASSERT (n_ref_fns, "Invalid value in $GENOZIP_REFERENCE=\"%s\"- expecting a reference file name or two file names separated by a ':'", env);

        ASSERT (ref == gref || n_ref_fns==2, "Invalid value in $GENOZIP_REFERENCE=\"%s\"- expecting two reference file names (Primary and Luft) separated by a ':'. See" WEBSITE_DVCF, env);
        
        filename = ref_fns[ref==gref]; // first name is primary, second is luft
        filename_len = ref_fn_lens[ref==gref];

        WARN ("Note: using the reference file %.*s set in $GENOZIP_REFERENCE. You can override this with --reference", filename_len, filename);
    }
    else
        filename_len = strlen (filename);

    static int num_explicit = 0; // user can have up to 2 --reference arguments
    ASSINP0 (!is_explicit || (++num_explicit <= 2), "More than two --reference arguments");
    
    // no need for a reference if we're just doing "genocat --show-chain myfile.chain.genozip" (--show-chain only works for genocat and chain files)
    if (!flag.show_chain) {

        // case two --reference arguments: we move the first to prim_ref, and the second will be gref (destpination ref)
        // note: we we read the first argument, we didn't yet know if there is another one
        if (num_explicit == 2) {
            SWAP (gref, prim_ref);
            ref = gref;
        }

        // case: pizzing subsequent files with implicit reference (reference from file header)
        if (!is_explicit && ref->filename) {
            if (!strcmp (filename, ref->filename)) return; // same file - we're done

            // in case a different reference is loaded - destroy it
            ref_destroy_reference (ref, false);
        }
    
        flag.reference    = ref_type; 
        flag.explicit_ref = is_explicit;
    }

    ref->filename = CALLOC (filename_len + 1);
    memcpy ((char*)ref->filename, filename, filename_len);
}

// called when loading an external reference
void ref_set_ref_file_info (Reference ref, Digest md5, const char *fasta_name, uint8_t genozip_version)
{
    ref->file_md5 = md5;
    ref->genozip_version = genozip_version;

    if (fasta_name[0]) {
        ref->ref_fasta_name = MALLOC (strlen (fasta_name) + 1); 
        strcpy (ref->ref_fasta_name, fasta_name);
    }
}

// display the reference:
// show a subset of the reference if --regions is specified - but only up to one region per chromosome
void ref_display_ref (Reference ref)
{
    Context chrom_ctx = {};
    ref_load_external_reference (ref, &chrom_ctx);

    if (flag.regions) regions_make_chregs (&chrom_ctx);

    for (const Range *r = FIRSTENT (Range, ref->ranges); r < AFTERENT (Range, ref->ranges); r++) {

        unsigned num_intersections = regions_get_num_range_intersections (r->chrom);
        if (!num_intersections) continue;

        if (!flag.no_header && !flag.gpos) printf ("%.*s\n", r->chrom_name_len, r->chrom_name); // note: this runs stand-alone, so we always output to stdout, not info_stream
        
        for (unsigned inter_i=0; inter_i < num_intersections; inter_i++) {
         
            PosType display_first_pos, display_last_pos;
            bool revcomp;
            
            ASSERT0 (regions_get_range_intersection (r->chrom, r->first_pos, r->last_pos, inter_i, &display_first_pos, &display_last_pos, &revcomp), 
                     "unexpectedly, no intersection");

            int64_t adjust = flag.gpos ? r->gpos - r->first_pos : 0;

            // case: normal sequence
            if (!revcomp) {
                if (display_first_pos == display_last_pos)
                    iprintf ("%"PRIu64"\t", display_first_pos + adjust);
                else
                    iprintf ("%"PRIu64"-%"PRIu64"\t", display_first_pos + adjust, display_last_pos + adjust);

                for (PosType pos=display_first_pos, next_iupac_pos=display_first_pos ; pos <= display_last_pos ; pos++) {
                    char iupac = (pos==next_iupac_pos) ? ref_iupacs_get (gref, r, pos, false, &next_iupac_pos) : 0;
                    char base = iupac ? iupac : ref_base_by_pos (r, pos);
                    fputc (base, stdout);
                }
            }

            // case: reverse complement
            else {
                if (display_first_pos == display_last_pos)
                    iprintf ("COMPLEM %"PRIu64"\t", display_first_pos + adjust);
                else 
                    iprintf ("COMPLEM %"PRIu64"-%"PRIu64"\t", display_last_pos + adjust, display_first_pos + adjust);

                for (PosType pos=display_last_pos, next_iupac_pos=display_last_pos; pos >= display_first_pos; pos--) {
                    char iupac = (pos==next_iupac_pos) ? ref_iupacs_get (gref, r, pos, true, &next_iupac_pos) : 0;
                    char base = iupac ? iupac : ref_base_by_pos (r, pos);
                    fputc (COMPLEM[(int)base], stdout);
                }        
            }

            fputc ('\n', stdout);
        }
    }

    ctx_free_context (&chrom_ctx);
}

#define REV_CODEC_GENOME_BASES_PER_THREAD (1 << 27) // 128Mbp

static Reference ref_reverse_compliment_genome_ref = 0; // ref_generate_reverse_complement_genome is called from the main thread so no thread safety issues
static void ref_reverse_compliment_genome_prepare (VBlock *vb)
{
    vb->ref = ref_reverse_compliment_genome_ref;
    vb->ready_to_dispatch = (vb->vblock_i-1) * REV_CODEC_GENOME_BASES_PER_THREAD < vb->ref->genome_nbases;
}

static void ref_reverse_compliment_genome_do (VBlock *vb)
{
    bit_array_reverse_complement_all (vb->ref->emoneg, vb->ref->genome, 
                                      (vb->vblock_i-1) * REV_CODEC_GENOME_BASES_PER_THREAD, 
                                      REV_CODEC_GENOME_BASES_PER_THREAD);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

void ref_generate_reverse_complement_genome (Reference ref)
{
    START_TIMER;
    ref_reverse_compliment_genome_ref = ref;
    dispatcher_fan_out_task ("generate_rev_comp_genome", NULL, PROGRESS_NONE, 0, false, true, true, false, 0, 10,
                             ref_reverse_compliment_genome_prepare, 
                             ref_reverse_compliment_genome_do, 
                             NULL);
    COPY_TIMER_VB (evb, generate_rev_complement_genome);
}

bool ref_is_external_loaded (Reference ref)
{
    return ref->external_ref_is_loaded;
}

// do we have a reference (EXTERNAL, STORED, ....)
bool ref_is_loaded (Reference ref)
{
    return ref->genome && ref->genome->nbits;
}

// ZIP & PIZ: import external reference
void ref_load_external_reference (Reference ref, ContextP chrom_ctx)
{
    ASSERTNOTNULL (ref->filename);
    SAVE_FLAGS_AUX;

    flag.reading_reference = ref; // tell file.c and fasta.c that this is a reference
    flag.no_writer = true;
    flag.luft = false;
    flag.list_chroms = false;
    flag.show_gheader = false;

    z_file = file_open (ref->filename, READ, Z_FILE, DT_FASTA);    
    z_file->basename = file_basename (ref->filename, false, "(reference)", NULL, 0);

    TEMP_VALUE (command, PIZ);

    // the reference file has no components or VBs - it consists of only a global area including reference sections
    piz_read_global_area (ref);
    
    // case: we are requested the chrom context - word_list and dict
    if (chrom_ctx) {
        buf_copy (evb, &chrom_ctx->word_list, &ZCTX(REF_CONTIG)->word_list, CtxWord, 0, 0, "contexts->word_list");
        buf_copy (evb, &chrom_ctx->dict, &ZCTX(REF_CONTIG)->dict, char, 0, 0, "contexts->word_list");
    }

    file_close (&z_file, false, false);
    file_close (&txt_file, false, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtheader_piz_read_and_reconstruct.

    ref->external_ref_is_loaded = true;

    // recover globals
    RESTORE_VALUE (command);
    RESTORE_FLAGS;
}

static void ref_initialize_loaded_ranges (Reference ref, RangesType type)
{
    if (flag.reading_reference) {
        buf_copy (evb, &ref->ref_external_ra, &z_file->ra_buf, RAEntry, 0, 0, "ref_external_ra");
        buf_copy (evb, &ref->ref_file_section_list, &z_file->section_list_buf, SectionEnt, 0, 0, "ref_file_section_list");
    }

    // notes on ranges.len:
    // 1. in case stored reference originates from REF_INTERNAL, we have a contig for every chrom for which there was a sequence,
    // but there might be additional chroms in the SAM file (eg in RNEXT) that don't have a sequence. Since we're indexing ranges by chrom,
    // we need a range for each chrom, even if it doesn't have data 
    // 2. in case stored reference originates from REF_EXT_STORE, we have contigs for all the chroms in the original reference file,
    // but our CHROM contexts also includes alternate chrom names that aren't in the original reference that appear after the reference
    // chroms. we need to make sure ranges.len does not include alternate chroms as that's how we know a chrom is alternate in ref_piz_get_range
    // 3. in case loading from a reference file, the number of contigs will match the number of chroms, so no issues.
    ref->ranges.len = IS_REF_INTERNAL (z_file) ? ZCTX(CHROM)->word_list.len : ref->ctgs.contigs.len;

    buf_alloc (evb, &ref->ranges, ref->ranges.len, 0, Range, 1, "ranges");     
    buf_zero (&ref->ranges);
    ranges_type(ref) = type;

    Context *chrom_ctx = ZCTX(CHROM);

    for (uint32_t range_id=0; range_id < ref->ranges.len; range_id++) {
        Range *r = ENT (Range, ref->ranges, range_id);
        r->range_id = r->chrom = range_id;

        if (flag.reference == REF_STORED) // PIZ
            ctx_get_snip_by_word_index (chrom_ctx, r->chrom, &r->chrom_name, &r->chrom_name_len);
        else
            ref_contigs_get_name_by_ref_index (ref, r->chrom, &r->chrom_name, &r->chrom_name_len);
    }

    // we don't need is_set if we're compressing with REF_EXTERNAL 

    // note: genome_nbases must be full words, so that bit_array_reverse_complement doesn't need to shift
    ref->genome_nbases = ROUNDUP64 (ref_contigs_get_genome_nbases (ref)) + 64; // round up to the nearest 64 bases, and add one word, needed by aligner_get_match_len for bit shifting overflow

    if (ref_has_is_set()) 
        ref->genome_is_set = buf_alloc_bitarr (evb, &ref->genome_is_set_buf, ref->genome_nbases, "genome_is_set_buf");

    // we protect genome->ref while uncompressing reference data, and genome->is_set while segging
    ref_lock_initialize_loaded_genome (ref);
}

static void overlay_ranges_on_loaded_genome (Reference ref)
{
    // overlay all chromosomes (range[i] goes to chrom_index=i) - note some chroms might not have a contig in 
    // which case their range is not initialized
    for (Range *r = FIRSTENT (Range, ref->ranges) ; r < AFTERENT (Range, ref->ranges); r++) {
        r->chrom = ENTNUM (ref->ranges, r);
        const Contig *rc = ref_contigs_get_contig_by_ref_index (ref, r->chrom, true);

        if (rc) { // this chromosome has reference data 
            r->gpos      = rc->gpos;
            r->first_pos = rc->min_pos;
            r->last_pos  = rc->max_pos;
            ref_contigs_get_name_by_ref_index (ref, r->chrom, &r->chrom_name, &r->chrom_name_len);
            
            PosType nbases = rc->max_pos - rc->min_pos + 1;

            ASSERT (r->gpos + nbases <= ref->genome->nbits / 2, "adding range \"%s\": r->gpos(%"PRId64") + nbases(%"PRId64") (=%"PRId64") is beyond genome_nbases=%"PRId64,
                    r->chrom_name, r->gpos, nbases, r->gpos+nbases, ref->genome_nbases);

            bit_array_overlay (&r->ref, ref->genome, r->gpos*2, nbases*2);

            if (ref_has_is_set()) 
                bit_array_overlay (&r->is_set, ref->genome_is_set, r->gpos, nbases);
        }
    }
}

// case 1: in case of ZIP with external reference, called by ref_load_stored_reference during piz_read_global_area of the reference file
// case 2: in case of PIZ: also called from ref_load_stored_reference with RT_LOADED
// case 3: in case of ZIP of SAM using internal reference - called from sam_zip_initialize
// note: ranges allocation must be called by the main thread as it adds a buffer to evb buf_list
void ref_initialize_ranges (Reference ref, RangesType type)
{
    if (type == RT_LOADED || type == RT_CACHED) {
        ref_initialize_loaded_ranges (ref, type);

        if (type == RT_LOADED) 
            buf_alloc (evb, &ref->genome_cache, 0, ref->genome_nbases / 4 * 2, uint8_t, 1, "genome_cache"); // contains both forward and rev. compliment
        
        else  // RT_CACHED 
            ASSERT0 (buf_mmap (evb, &ref->genome_cache, ref_get_cache_fn(ref), false, "genome_cache"),  // we map the entire file (forward and revese complement genomes) onto genome_cache
                     "failed to map cache. Please try again");

        // overlay genome and emoneg. we do it this was so we can use just a single file
        buf_set_overlayable (&ref->genome_cache);
        ref->genome = buf_overlay_bitarr (evb, &ref->genome_buf, &ref->genome_cache, 0, ref->genome_nbases * 2, "genome_buf");
        ref->emoneg = buf_overlay_bitarr (evb, &ref->emoneg_buf, &ref->genome_cache, ref->genome_nbases / 4, ref->genome_nbases * 2, "emoneg_buf");

        overlay_ranges_on_loaded_genome (ref);
    }

    else { // RT_DENOVO
        if (buf_is_alloc (&ref->ranges)) return; // case: 2nd+ bound file

        ref->ranges.len = REF_NUM_DENOVO_RANGES;
        ranges_type(ref) = RT_DENOVO;
        buf_alloc (evb, &ref->ranges, 0, ref->ranges.len, Range, 1, "ranges"); 
        buf_zero (&ref->ranges);
        ref_lock_initialize_denovo_genome(ref);
    }
}

typedef struct { PosType min_pos, max_pos; } MinMax;

//---------------------------------------
// Printing
//---------------------------------------
RangeStr ref_display_range (const Range *r)
{
    RangeStr s;

    sprintf (s.s, "range_id=%u ref.num_bits=%"PRIu64" is_set.num_bits=%"PRIu64" chrom_name=%.*s ref_index=%d range_i=%u first_pos=%"PRId64
             " last_pos=%"PRId64" gpos=%"PRId64,
             r->range_id, r->ref.nbits, r->is_set.nbits, r->chrom_name_len, r->chrom_name, r->chrom, r->range_i, 
             r->first_pos, r->last_pos, r->gpos);
    return s;
}

void ref_print_subrange (const char *msg, const Range *r, PosType start_pos, PosType end_pos, FILE *file) /* start_pos=end_pos=0 if entire ref */
{
    uint64_t start_idx = start_pos ? start_pos - r->first_pos : 0;
    uint64_t end_idx   = (end_pos ? MIN_(end_pos, r->last_pos) : r->last_pos) - r->first_pos;

    fprintf (file, "%s: %.*s %"PRId64" - %"PRId64" (len=%u): ", msg, r->chrom_name_len, r->chrom_name, start_pos, end_pos, (uint32_t)(end_pos - start_pos + 1));
    for (uint64_t idx = start_idx; idx <= end_idx; idx++) 
        fputc (ref_base_by_idx (r, idx) + (32 * !ref_is_nucleotide_set (r, idx)), file); // uppercase if set, lowercase if not

    fputc ('\n', file);
}

// outputs in seq, a nul-terminated string of up to (len-1) bases
char *ref_dis_subrange (Reference ref, const Range *r, PosType start_pos, PosType len, char *seq, bool revcomp) // in_reference allocated by caller to length len
{
    if (!revcomp) {
        PosType next_iupac_pos, pos, end_pos = MIN_(start_pos + len - 1 - 1, r->last_pos);  // -1 to leave room for \0
        for (pos=start_pos, next_iupac_pos=start_pos ; pos <= end_pos; pos++) {
            char iupac = (pos==next_iupac_pos) ? ref_iupacs_get (ref, r, pos, false, &next_iupac_pos) : 0;
            seq[pos - start_pos] = iupac ? iupac : ref_base_by_pos (r, pos);
        }
        seq[pos - start_pos] = 0;
    }

    // revcomp: display the sequence starting at start_pos and going backwards - complemented
    else {
        PosType next_iupac_pos, pos, end_pos = MAX_(start_pos - (len - 1) + 1, r->first_pos);  // -1 to leave room for \0
        for (pos=start_pos, next_iupac_pos=start_pos ; pos >= end_pos; pos--) {
            char iupac = (pos==next_iupac_pos) ? ref_iupacs_get (gref, r, pos, true, &next_iupac_pos) : 0;
            seq[start_pos - pos] = iupac ? iupac : COMPLEM[(int)ref_base_by_pos (r, pos)];
        }
        seq[start_pos - pos] = 0;
    }

    return seq;
}

void ref_print_is_set (const Range *r,
                       PosType around_pos,  // display around this neighborhoud ; -1 means entire range
                       FILE *file)
{
#   define neighborhood (PosType)10000

    fprintf (file, "\n\nRegions set for ref_index %d \"%.*s\" [%"PRId64"-%"PRId64"] according to range.is_set (format- \"first_pos-last_pos (len)\")\n", 
             r->chrom, r->chrom_name_len, r->chrom_name, r->first_pos, r->last_pos);
    fprintf (file, "In the neighborhood of about %u bp around pos=%"PRId64"\n", (unsigned)neighborhood, around_pos);

    if (!r->is_set.nbits)
        fprintf (file, "No data: r->is_set.nbits=0\n");

    if (around_pos < r->first_pos || around_pos > r->last_pos)
        fprintf (file, "No data: pos=%"PRId64" is outside of [first_pos=%"PRId64" - last_pos=%"PRId64"]\n", around_pos, r->first_pos, r->last_pos);

    uint64_t next;
    for (uint64_t i=0; i < r->is_set.nbits; ) {
        
        bool found = bit_array_find_next_clear_bit (&r->is_set, i, &next);
        if (!found) next = r->is_set.nbits;

        bool in_neighborhood = (around_pos - (PosType)(r->first_pos+i) > -neighborhood) && (around_pos - (PosType)(r->first_pos+i) < neighborhood);
        if (next > i && (around_pos == -1 || in_neighborhood)) {
            if (next - i > 1)
                fprintf (file, "%"PRId64"-%"PRIu64"(%u)\t", r->first_pos + i, r->first_pos + next-1, (uint32_t)(next - i));
            else
                fprintf (file, "%"PRId64"(1)\t", r->first_pos + i);
        }                   
        if (!found) break;

        i = next;

        found = bit_array_find_next_set_bit (&r->is_set, i, &next);
        if (!found) next = r->is_set.nbits;

        i = next;
    }
    fprintf (file, "\n");
}

// returns the reference file name for CRAM, derived from the genozip reference name
const char *ref_get_cram_ref (Reference ref)
{
    static char *samtools_T_option = NULL;
    if (samtools_T_option) goto done; // already calculated

    ASSINP0 (ref->filename, "when compressing a CRAM file, --reference or --REFERENCE must be specified");

    // if we're attempting to open a cram file, just to check whether it is aligned (in main_load_reference),
    // then we haven't loaded the reference file yet, and hence we don't know ref_fasta_name.
    // in that case, we will just load the reference file's header
    z_file = file_open (ref->filename, READ, Z_FILE, DT_FASTA);    
    flag.reading_reference = gref;
    zfile_read_genozip_header (0);
    flag.reading_reference = NULL;
    file_close (&z_file, false, true);


    ASSINP (ref->ref_fasta_name, "cannot compress a CRAM file because %s is lacking the name of the source fasta file - likely because it was created by piping a fasta from from stdin, or because the name of the fasta provided exceed %u characters",
            ref->filename, REF_FILENAME_LEN-1);

    samtools_T_option = MALLOC (MAX_(strlen (ref->ref_fasta_name), strlen (ref->filename)) + 10);

    // case: fasta file is in its original location
    if (file_exists (ref->ref_fasta_name)) 
        sprintf (samtools_T_option, "-T%s", ref->ref_fasta_name);

    // try: fasta file is in directory of reference file
    else {
        const char *slash = strrchr (ref->ref_fasta_name, '/');
        if (!slash) slash = strrchr (ref->ref_fasta_name, '\\'); 
        const char *basename = slash ? slash+1 : ref->ref_fasta_name;

        slash = strrchr (ref->filename, '/');
        if (!slash) slash = strrchr (ref->filename, '\\'); 
        unsigned dirname_len = slash ? ((slash+1) - ref->filename) : 0;

        sprintf (samtools_T_option, "-T%.*s%s", dirname_len, ref->filename, basename);

        ASSINP (file_exists (&samtools_T_option[2]), 
                "Searching of the fasta file used to create %s. It was not found in %s or %s. Note: it is needed for passing to samtools as a reference file (-T option) for reading the CRAM file", 
                ref->filename, ref->ref_fasta_name, &samtools_T_option[2]);        
    }

done:
    return samtools_T_option;
}
