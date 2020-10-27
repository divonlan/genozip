// ------------------------------------------------------------------
//   reference.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <pthread.h>
#include <errno.h>
#include "reference.h"
#include "buffer.h"
#include "strings.h"
#include "dict_id.h"
#include "dispatcher.h"
#include "zip.h"
#include "zfile.h"
#include "endianness.h"
#include "random_access.h"
#include "seg.h"
#include "piz.h"
#include "vblock.h"
#include "context.h"
#include "hash.h"
#include "mutex.h"
#include "sections.h"
#include "profiler.h"
#include "bit_array.h"
#include "file.h"
#include "regions.h"
#include "refhash.h"
#include "ref_private.h"
#include "bit_array.h"
#include "codec.h"
#include "compressor.h"

// ZIP and PIZ, internal or external reference ranges. If in ZIP-INTERNAL we have REF_NUM_DENOVO_RANGES Range's - each allocated on demand. In all other cases we have one range per contig.
Buffer ranges                       = EMPTY_BUFFER; 

Range genome = {}, genome_rev = {}; // used in ZIP/PIZ of fasta and fastq with a reference - we have only one range that is the entire genome (this points to the first entry in ranges)
Buffer genome_ref_copy              = EMPTY_BUFFER; // saved copy of genome.ref, used when compressing multiple files with REF_EXT_STORE - which compacts ranges - recovering between files
Buffer ranges_copy                  = EMPTY_BUFFER;     // saved copy of ranges in initialize state, for same purpose

PosType genome_size = 0;

static Buffer ref_external_ra       = EMPTY_BUFFER; // Random Access data of the external reference file
static Buffer ref_file_section_list = EMPTY_BUFFER; // Section List of the external reference file

// PIZ: an array of RegionToSet - list of non-compacted regions, for which we should set is_set bit to 1
static Buffer region_to_set_list    = EMPTY_BUFFER; 
SPINLOCK (region_to_set_list_spin);

// ZIP/PIZ: random_access data for the reference sections stored in a target genozip file
Buffer ref_stored_ra                = EMPTY_BUFFER;
SPINLOCK (ref_stored_ra_spin); // ZIP only

typedef struct {
    BitArray *is_set;
    PosType first_bit, len;
} RegionToSet;

static SectionListEntry *sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
static uint32_t ref_range_cursor = 0;

static char *ref_fasta_name = NULL;

// globals
const char *ref_filename = NULL; // filename of external reference file
Md5Hash ref_md5 = MD5HASH_NONE;

#define CHROM_GENOME 0
#define CHROM_NAME_GENOME "GENOME"

#define CHROM_GENOME_REV 1
#define CHROM_NAME_GENOME_REV "GENOME_REV"

#define ref_is_range_used(r) ((r)->ref.num_of_bits && ((r)->is_set.num_of_bits || flag_make_reference))

// free memory allocations between files, when compressing multiple non-bound files or decompressing multiple files
void ref_unload_reference (bool force_clean_all)
{
    // in case of REF_EXTERNAL/REF_EXT_STORE, we keep the reference for other files on the command line
    // note: in PIZ we never cleanup external references as they can never change (the user can only specify one --reference)
    if ((flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE) && !force_clean_all) {  
        
        // zip sets 1 to nucleotides used. we clear it for the next file.
        if (primary_command == ZIP)
            bit_array_clear_all (&genome.is_set); 

        // REF_EXT_STORE compacts ranges. We recover them from the copy
        if (buf_is_allocated (&genome_ref_copy)) { // won't be allocated if we are compressing only one file or if not REF_EXT_STORE
            memcpy (genome.ref.words, genome_ref_copy.data, genome.ref.num_of_words * sizeof (word_t));
            buf_copy (evb, &ranges, &ranges_copy, sizeof (Range), 0, 0, "ranges", 0);
        }
    }

    // case ZIP: REF_INTERNAL - we're done with ranges - free so next non-bound file can build its own
    // case ZIP: REF_EXT_STORE - unfortunately we will need to re-read the FASTA as we damaged it with removing flanking regions
    // case PIZ: REF_STORED - cleanup for next file with its own reference data
    // case force_clean_all: happens before start the test process (in ZIP with --test) for the last file

    else {
        if (buf_is_allocated (&ranges) && ranges.param == RT_DENOVO) {
            ARRAY (Range, rng, ranges);
            for (unsigned i=0; i < ranges.len ; i++) {
                FREE (rng[i].ref.words);
                FREE (rng[i].is_set.words);
                if (primary_command == ZIP) 
                    FREE (rng[i].chrom_name); // allocated only in ZIP/REF_INTERNAL - otherwise a pointer into another Buffer
            }
        }
        
        if (ranges.param == RT_LOADED) {
            FREE (genome.ref.words);
            FREE (genome.is_set.words);
            FREE (genome_rev.ref.words);
            genome = genome_rev = (Range){};
            genome_size = 0;
        }

        buf_free (&ranges);
        buf_free (&region_to_set_list);
        buf_free (&ref_external_ra);
        buf_free (&ref_stored_ra);
        buf_free (&ref_file_section_list);
        buf_free (&genome_ref_copy);
        buf_free (&ranges_copy);

        ref_contigs_free();
        ref_lock_free();
        refhash_free();
    }
}

MemStats ref_memory_consumption (void)
{
    MemStats stats = { "reference", 0, 0 };

    if (ranges.param == RT_LOADED) {
        stats.bytes += (genome.ref.num_of_words + genome.is_set.num_of_words + genome_rev.ref.num_of_words) * sizeof (word_t);
        stats.buffers = !!genome.ref.num_of_words + !!genome.is_set.num_of_words + !!genome_rev.ref.num_of_words;
    }

    else if (ranges.param == RT_DENOVO) {
        ARRAY (Range, r, ranges);
        for (unsigned i=0; i < ranges.len; i++) {
            if (r[i].ref.num_of_words) {
                stats.bytes += r[i].ref.num_of_words * sizeof (word_t);
                stats.buffers++;
            }
            if (r[i].is_set.num_of_words) {
                stats.bytes += r[i].is_set.num_of_words * sizeof (word_t);
                stats.buffers++;
            }
        }
    }

    return stats;
}

// PIZ: returns a range which is the entire contig
const Range *ref_piz_get_range (VBlockP vb, PosType first_pos_needed, uint32_t num_nucleotides_needed)
{
    // caching
    if (vb->prev_range && vb->prev_range_chrom_node_index == vb->chrom_node_index)
        return vb->prev_range;

    // gets the index of the matching chrom in the reference - either its the chrom itself, or one with an alternative name
    // eg 'chr22' instead of '22'
    uint32_t index = buf_is_allocated (&z_file->alt_chrom_map) ? *ENT (WordIndex, z_file->alt_chrom_map, vb->chrom_node_index)
                                                               : vb->chrom_node_index;
    Range *r = ENT (Range, ranges, index);
    if (!r->ref.num_of_words) return NULL; // this can ligitimately happen if entire chromosome is verbatim in SAM, eg. unaligned (pos=4) or SEQ or CIGAR are unavailable

    if (first_pos_needed + num_nucleotides_needed - 1 <= r->last_pos) {
        // this can happen if the reference originated from REF_INTERNAL, and the latter part of the range requested
        // originated from a missing range due to hash contention, and this missing part also happens to be
        // at the end of the chromosome, thereby causing r->last_pos to be less than the length of the chromosome. 
        // we return the range anyway, as the missing parts will have is_set=0 and seq_bitmap=0, retrieving from nonref

        // TODO: r->last_pos in REF_INTERNAL should include the missing ranges at the end of the chromosome
    }

    vb->prev_range = r;
    vb->prev_range_chrom_node_index = vb->chrom_node_index;

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
        ASSERT (len_1 > 0, "Error in ref_uncompact_ref: len_1 is not positive: start_0_offset=%"PRId64" start_1_offset=%"PRId64" first_bit=%"PRId64" last_bit=%"PRId64, 
                start_0_offset, start_1_offset, first_bit, last_bit);

        // do actual uncompacting
        bit_array_copy (&r->ref, start_1_offset * 2, compacted, next_compacted * 2, len_1 * 2);
        next_compacted += len_1;

        if (start_0_offset > last_bit) break; // we're done (we always end with a region of 1s because we removed the flanking 0s during compacting)

        // skip the clear region
        has_any_bit = bit_array_find_next_set_bit (&r->is_set, start_0_offset, &start_1_offset); 
        ASSERT0 (has_any_bit, "Error in ref_uncompact_ref: cannot find next set bit");
        ASSERT (start_1_offset <= last_bit, "Error in ref_uncompact_ref: expecting start_1_offset(%"PRId64") <= last_bit(%"PRId64")",
                start_1_offset, last_bit); // we removed the flanking regions, so there is always an 1 after a 0 within the region
    }

    ASSERT (next_compacted * 2 == compacted->num_of_bits, "Error in ref_uncompact_ref: expecting next_compacted(%"PRId64") * 2 == compacted->num_of_bits(%"PRId64")",
            next_compacted, compacted->num_of_bits);
}

// Compute thread: called by ref_uncompress_one_range
static Range *ref_get_range_by_chrom (WordIndex chrom, const char **chrom_name)
{
    Context *ctx = &z_file->contexts[CHROM];
    ASSERT (chrom >= 0 && chrom < ctx->word_list.len, "Error in ref_get_range_by_chrom: chrom=%d out of range - ctx->word_list.len=%u",
            chrom, (uint32_t)ctx->word_list.len);

    // if chrom is not in reference, get the reference chrom to which it is mapped (eg. 'chr22' -> '22')
    //chrom = buf_is_allocated (&z_file->alt_chrom_map) ? *ENT (WordIndex, z_file->alt_chrom_map, chrom) : chrom; <--- need reverse lookup - file chrom by alt chrom

    uint32_t chrom_name_len;
    mtf_get_snip_by_word_index (&ctx->word_list, &ctx->dict, chrom, chrom_name, &chrom_name_len);
    Range *r = ENT (Range, ranges, chrom); // in PIZ, we have one range per chrom

    return r;
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
            if (!flag_sequential && (i-start_base*2) % (BASES_PER_LINE*2) == 0)
                fprintf (stderr, "%8"PRIu64": ", i/2);
            fputc (fwd[bit_array_get(bitarr, i+1)][bit_array_get(bitarr, i)], file);
            if (!flag_sequential && ((i-start_base*2) % (BASES_PER_LINE*2) == 2*(BASES_PER_LINE-1))) fputc ('\n', file);
        }
    
    else 
        for (int64_t i=(start_base+num_of_bases-1)*2; i >= start_base*2; i-=2) { // signed type
            fputc (rev[bit_array_get(bitarr, i+1)][bit_array_get(bitarr, i)], file);
            if (!flag_sequential && (((start_base+num_of_bases-1)*2-i) % (BASES_PER_LINE*2) == (BASES_PER_LINE-1)*2)) fputc ('\n', file);
        }
    
    fputc ('\n', file);
}

static void ref_show_sequence (void)
{
    for (uint32_t range_i=0; range_i < ranges.len; range_i++) {
        Range *r = ENT (Range, ranges, range_i);

        // get first pos and last pos, potentially modified by --regions
        PosType first_pos, last_pos;
        if (!r->ref.num_of_bits ||
            !regions_get_range_intersection (r->chrom, r->first_pos, r->last_pos, &first_pos, &last_pos)) continue;

        if (r->ref.num_of_bits) {
            fprintf (stderr, "%.*s\n", r->chrom_name_len, r->chrom_name);
            ref_print_bases (stderr, &r->ref, first_pos, last_pos-first_pos+1, true);
        }
    }

    if (exe_type == EXE_GENOCAT) exit_ok;  // in genocat this, not the data
}

// entry point of compute thread of reference decompression. this is called when pizzing a file with a stored reference,
// including reading the reference file itself.
// vb->z_data contains a SEC_REFERENCE section and sometimes also a SEC_REF_IS_SET sections
static void ref_uncompress_one_range (VBlockP vb)
{
    if (!buf_is_allocated (&vb->z_data) || !vb->z_data.len) goto finish; // we have no data in this VB because it was skipped due to --regions or genocat --show-headers

    SectionHeaderReference *header = (SectionHeaderReference *)vb->z_data.data;

    WordIndex chrom          = (WordIndex)BGEN32 (header->chrom_word_index);
    uint32_t uncomp_len      = BGEN32 (header->h.data_uncompressed_len);
    PosType ref_sec_pos      = (PosType)BGEN64 (header->pos);
    PosType ref_sec_gpos     = (PosType)BGEN64 (header->gpos);
    PosType ref_sec_len      = (PosType)BGEN32 (header->num_bases);
    PosType ref_sec_last_pos = ref_sec_pos + ref_sec_len - 1;
    PosType compacted_ref_len=0, initial_flanking_len=0, final_flanking_len=0; 

    const char *chrom_name;
    Range *r = ref_get_range_by_chrom (chrom, &chrom_name);
    PosType sec_start_within_contig = ref_sec_pos - r->first_pos;
    PosType sec_start_gpos          = r->gpos + sec_start_within_contig;
    PosType sec_end_within_contig   = sec_start_within_contig + ref_sec_len - 1;
    
    bool is_compacted = (header->h.section_type == SEC_REF_IS_SET); // we have a SEC_REF_IS_SET if  SEC_REFERENCE was compacted

    if (flag_show_reference && primary_command == PIZ && r)  // in ZIP, we show the compression of SEC_REFERENCE into z_file, not the uncompression of the reference file
        fprintf (stderr, "vb_i=%u Uncompressing %-14s chrom=%u (%.*s) gpos=%"PRId64" pos=%"PRId64" num_bases=%u Bytes=%u\n", 
                vb->vblock_i, st_name (header->h.section_type), BGEN32 (header->chrom_word_index), r->chrom_name_len, r->chrom_name, BGEN64 (header->gpos), 
                BGEN64 (header->pos), BGEN32 (header->num_bases), BGEN32 (header->h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));

    // initialization of is_set:
    // case 1: ZIP (reading an external reference) - we CLEAR is_set, and let seg set the bits that are to be
    //    needed from the reference for pizzing (so we don't store the others)
    //    note: in case of ZIP with REF_INTERNAL, we CLEAR the bits in ref_seg_get_locked_range
    // case 2: PIZ, reading an uncompacted (i.e. complete) reference section - which is always the case when
    //    reading an external reference and sometimes when reading an stored one - we SET all the bits as they are valid for pizzing
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
                "Error in ref_uncompress_one_range: section range out of bounds for chrom=%d \"%s\": in SEC_REFERENCE being uncompressed: first_pos=%"PRId64" last_pos=%"PRId64" but in reference contig as initialized: first_pos=%"PRId64" last_pos=%"PRId64,
                chrom, r->chrom_name, ref_sec_pos, ref_sec_last_pos, r->first_pos, r->last_pos);

        ASSERT (uncomp_len == roundup_bits2bytes64 (ref_sec_len), "Error in ref_uncompress_one_range: when uncompressing SEC_REF_IS_SET: uncomp_len=%u inconsistent with len=%"PRId64, uncomp_len, ref_sec_len); 

        // uncompress into r->is_set, via vb->compressed
        ASSERT0 (!vb->compressed.len, "Error in ref_uncompress_one_range: expecting vb->compressed to be free, but its not");
        zfile_uncompress_section (vb, (SectionHeaderP)header, &vb->compressed, "compressed", 0, SEC_REF_IS_SET);

        BitArray *is_set = buf_zfile_buf_to_bitarray (&vb->compressed, ref_sec_len);

        // note on locking: while different threads uncompress regions of the range that are non-overlapping, 
        // there might be a 64b word that is split between two ranges
        RefLock lock = ref_lock (sec_start_gpos, ref_sec_len); 
        bit_array_copy (&r->is_set, sec_start_within_contig, is_set, 0, ref_sec_len); // initialization of is_set - case 3
        ref_unlock (lock);

        buf_free (&vb->compressed);

        // display contents of is_set if user so requested
        if (flag_show_is_set && !strcmp (chrom_name, flag_show_is_set)) 
            ref_print_is_set (r, -1);

        // prepare for uncompressing the next section - which is the SEC_REFERENCE
        header = (SectionHeaderReference *)&vb->z_data.data[*ENT (uint32_t, vb->z_section_headers, 1)];

        if (flag_show_reference && primary_command == PIZ && r) 
            fprintf (stderr, "vb_i=%u Uncompressing %-14s chrom=%u (%.*s) gpos=%"PRId64" pos=%"PRId64" num_bases=%u Bytes=%u\n", 
                    vb->vblock_i, st_name (header->h.section_type), BGEN32 (header->chrom_word_index), r->chrom_name_len, r->chrom_name, BGEN64 (header->gpos), 
                    BGEN64 (header->pos), BGEN32 (header->num_bases), BGEN32 (header->h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));

        compacted_ref_len  = (PosType)BGEN32(header->num_bases);
        uncomp_len         = BGEN32 (header->h.data_uncompressed_len);

        ASSERT (uncomp_len == roundup_bits2bytes64 (compacted_ref_len*2), 
                "Error: uncomp_len=%u inconsistent with compacted_ref_len=%"PRId64, uncomp_len, compacted_ref_len); 

        ASSERT0 (BGEN32 (header->chrom_word_index) == chrom && BGEN64 (header->pos) == ref_sec_pos && BGEN64 (header->gpos) == ref_sec_gpos, // chrom should be the same between the two sections
                "Error in ref_uncompress_one_range: header mismatch between SEC_REF_IS_SET and SEC_REFERENCE sections");
    }
    
    // case: not compacted means that entire range is set
    else {
        ASSERT (uncomp_len == roundup_bits2bytes64 (ref_sec_len*2), "Error: uncomp_len=%u inconsistent with ref_len=%"PRId64, uncomp_len, ref_sec_len); 

        if (primary_command == ZIP && flag_reference == REF_EXT_STORE) { // initialization of is_set - case 1
            RefLock lock = ref_lock (sec_start_gpos, ref_sec_len);
            bit_array_clear_region (&r->is_set, sec_start_within_contig, ref_sec_len); // entire range is cleared
            ref_unlock (lock);
        }

        else if (primary_command == PIZ) { // initialization of is_set - case 2

            // it is possible that the section goes beyond the boundaries of the contig, this can happen when we compressed with --REFERENCE
            // and the section was copied in its entirety from the .ref.genozip file (in ref_copy_one_compressed_section)
            // even though a small amount of flanking regions is not set. in this case, we copy from the section only the part needed
            initial_flanking_len = (sec_start_within_contig < 0)    ? -sec_start_within_contig       : 0; // nucleotides in the section that are before the start of our contig
            final_flanking_len   = (ref_sec_last_pos > r->last_pos) ? ref_sec_last_pos - r->last_pos : 0; // nucleotides in the section that are after the end of our contig

            bit_index_t start = MAX (sec_start_within_contig, 0);
            bit_index_t len   = ref_sec_len - initial_flanking_len - final_flanking_len;
            
            RefLock lock = ref_lock (start + r->gpos, len);
            bit_array_set_region (&r->is_set, start, len);
            ref_unlock (lock);

            // save the region we need to set, we will do the actual setting in ref_load_stored_reference
            spin_lock (region_to_set_list_spin);
            RegionToSet *rts = &NEXTENT (RegionToSet, region_to_set_list);
            spin_unlock (region_to_set_list_spin);
            rts->is_set    = &r->is_set;
            rts->first_bit = MAX (sec_start_within_contig, 0);
            rts->len       = ref_sec_len - initial_flanking_len - final_flanking_len;
        }
   
        if (!uncomp_len) return;  // empty header - if it appears, it is the final header (eg in case of an unaligned SAM file)
    }

    // uncompress into r->ref, via vb->compressed
    ASSERT0 (!vb->compressed.len, "Error in ref_uncompress_one_range: expecting vb->compressed to be free, but its not");
    zfile_uncompress_section (vb, (SectionHeaderP)header, &vb->compressed, "compressed", 0, SEC_REFERENCE);

    // lock - while different threads uncompress regions of the range that are non-overlapping, they might overlap at the bit level
    RefLock lock = ref_lock (ref_sec_gpos, ref_sec_len); 

    if (is_compacted) {
        const BitArray *compacted = buf_zfile_buf_to_bitarray (&vb->compressed, compacted_ref_len * 2);
        ref_uncompact_ref (r, sec_start_within_contig, sec_end_within_contig, compacted);
    }

    else {
        BitArray *ref = buf_zfile_buf_to_bitarray (&vb->compressed, ref_sec_len * 2);

        // copy the section, excluding the flanking regions
        bit_array_copy (&r->ref, MAX (sec_start_within_contig, 0) * 2, // dst
                        ref, initial_flanking_len * 2, // src
                        (ref_sec_len - initial_flanking_len - final_flanking_len) * 2); // len
    }

    ref_unlock (lock);
    buf_free (&vb->compressed);

finish:
    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. 
}

static void ref_read_one_range (VBlockP vb)
{
    if (!sections_get_next_section_of_type (&sl_ent, &ref_range_cursor, SEC_REFERENCE, SEC_REF_IS_SET) || // no more reference sections
        ((sl_ent+1)->offset - sl_ent->offset) == sizeof (SectionHeaderReference)) // final, header-only section sometimes exists (see ref_compress_ref)
        return; // we're done
    
    if (sl_ent->vblock_i == 0) // section was created with ref_copy_one_compressed_section
        z_file->num_copied_ref_sections++;
    else        
        ASSERT (sl_ent->vblock_i + z_file->num_copied_ref_sections == vb->vblock_i, 
                "Error in ref_read_one_range: mismatch: sl_ent->vblock_i=%u but vb->vblock_i=%u, z_file->num_copied_ref_sections=%u",
                sl_ent->vblock_i, vb->vblock_i, z_file->num_copied_ref_sections);

    // if the user specified --regions, check if this ref range is needed
    bool range_is_included = true;
    RAEntry *ra = NULL;
    if (flag_regions) {
        if (vb->vblock_i > ref_stored_ra.len) return; // we're done - no more ranges to read, per random access (this is the empty section)

        ra = ENT (RAEntry, ref_stored_ra, vb->vblock_i-1);
        ASSERT (ra->vblock_i == vb->vblock_i, "Error in ref_read_one_range: expecting ra->vblock_i(%u) == vb->vblock_i(%u)", 
                ra->vblock_i, vb->vblock_i);

        range_is_included = regions_is_ra_included (ra); 
    }

    if (range_is_included) { 

        buf_alloc (vb, &vb->z_section_headers, 2 * sizeof(int32_t), 0, "z_section_headers", 2); // room for 2 section headers  

        ASSERT0 (vb->z_section_headers.len < 2, "Error in ref_read_one_range: unexpected 3rd recursive entry");

        int32_t section_offset = 
            zfile_read_section (z_file, vb, sl_ent->vblock_i, &vb->z_data, "z_data", 
                                sizeof(SectionHeaderReference), sl_ent->section_type, sl_ent);    

        ASSERT (section_offset != EOF, "Error in ref_read_one_range: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);

        NEXTENT (int32_t, vb->z_section_headers) = section_offset;

        // allocate memory for entire chrom reference if this is the first range of this chrom
        SectionHeaderReference *header = (SectionHeaderReference *)&vb->z_data.data[section_offset];
        
        WordIndex chrom = BGEN32 (header->chrom_word_index);
        if (chrom == NODE_INDEX_NONE) return; // we're done - terminating empty section that sometimes appears (eg in unaligned SAM that don't have any reference and yet are REF_INTERNAL)
    }

    // if this is SEC_REF_IS_SET, read the SEC_REFERENCE section now (even if its not included - we need to advance the cursor)
    if (sl_ent->section_type == SEC_REF_IS_SET) 
        ref_read_one_range (vb);

    if (flag_show_headers && exe_type == EXE_GENOCAT) 
        vb->z_data.len = 0; // roll back if we're only showing headers

    vb->ready_to_dispatch = true; // to simplify the code, we will dispatch the thread even if we skip the data, but we will return immediately. 
}

// PIZ: loading a reference stored in the genozip file - this could have been originally stored as REF_INTERNAL or REF_EXT_STORE
// or this could be a .ref.genozip file (called from load_external->piz_one_file)
void ref_load_stored_reference (void)
{
    ASSERT0 (!buf_is_allocated (&ranges), "Error in ref_load_stored_reference: expecting ranges to be unallocated");
    
    if (!(flag_show_headers && exe_type == EXE_GENOCAT)) {
        random_access_pos_of_chrom (0, 0, 0); // initialize if not already initialized
        
        if (flag_reading_reference) {
            buf_copy (evb, &ref_external_ra, &z_file->ra_buf, sizeof (RAEntry), 0, 0, "ref_external_ra", 0);
            buf_copy (evb, &ref_file_section_list, &z_file->section_list_buf, sizeof (SectionListEntry), 0, 0, "ref_file_section_list", 0);
        }

        ref_initialize_ranges (RT_LOADED);
        
        sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
        ref_range_cursor = 0;

        spin_initialize (region_to_set_list_spin);
        buf_alloc (evb, &region_to_set_list, sections_count_sections (SEC_REFERENCE) * sizeof (RegionToSet), 1, "region_to_set_list", 0);
    }
    
    // decompress reference using Dispatcher
    bool external = flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE;
    dispatcher_fan_out_task (external ? ref_filename     : z_file->basename, 
                             external ? PROGRESS_MESSAGE : PROGRESS_NONE, 
                             external ? "Reading reference file..." : NULL, 
                             flag_test, 
                             ref_read_one_range, 
                             ref_uncompress_one_range, 
                             NULL);

    if (flag_show_ref_seq) 
        ref_show_sequence();
        
    // now we can safely set the is_set regions originating from non-compacted ranges. we couldn't do it before, because
    // copied-from-FASTA ranges appear first in the genozip file, and after them could be compacted ranges that originate
    // from a full-contig range in EXT_STORE, whose regions copied-from-FASTA are 0s.
    ARRAY (RegionToSet, rts, region_to_set_list);
    for (uint32_t i=0; i < region_to_set_list.len; i++)
        bit_array_set_region (rts[i].is_set, rts[i].first_bit, rts[i].len);

    buf_test_overflows_all_vbs ("ref_load_stored_reference");
}

// ------------------------------------
// ZIP side
// ------------------------------------

// case REF_INTERNAL - we get a range_id by hashing the chrom name and the range_i - we can't use the chrom_node_index
// because we have multiple threads in parallel that might have the same node_index for different chroms
static inline uint32_t ref_range_id_by_hash (VBlockP vb, uint32_t range_i)
{
    ASSERT0 (vb->chrom_name_len > 0, "Error in ref_range_id_by_hash: vb->chrom_name_len==0");

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

static Range *ref_seg_get_locked_range_denovo (VBlockP vb, PosType pos, const char *field /* used for ASSSEG */, RefLock *lock)  
{
    uint32_t range_i = pos2range_i (pos); // range within chromosome 

    // case: we're asking for the same range as the previous one (for example, subsequent line in a sorted SAM)
    if (vb && vb->prev_range && vb->prev_range_chrom_node_index == vb->chrom_node_index && vb->prev_range_range_i == range_i) {
        *lock = ref_lock_range (vb->prev_range - FIRSTENT (Range, ranges));
        return vb->prev_range;
    }

    uint32_t range_id = ref_range_id_by_hash (vb, range_i);
    ASSSEG (range_id < ranges.len, field, "Error in ref_seg_get_locked_range: range_id=%u expected to be smaller than ranges.len=%u", range_id, (uint32_t)ranges.len);

    Range *range = ENT (Range, ranges, range_id);
    *lock = ref_lock_range (range_id);

    if (range->ref.num_of_bits) {

        // check for hash conflict 
        if ((range->range_i != range_i || vb->chrom_name_len != range->chrom_name_len || memcmp (vb->chrom_name, range->chrom_name, vb->chrom_name_len))) {
            *lock = ref_unlock (*lock);

            ASSERTW (!flag_test_seg, "Warning: ref range contention: chrom=%.*s pos=%u (this slightly affects compression ratio, but is harmless)", 
                     vb->chrom_name_len, vb->chrom_name, (uint32_t)pos); // only show this in --test-seg

            range = NULL;  // no soup for you
        }

        return range; // range already initialized or cannot be initialized ^ - we're done
    }

    range->range_id       = range_id;
    range->range_i        = range_i;
    range->first_pos      = range_i2pos (range_i);
    range->last_pos       = range->first_pos + REF_NUM_DENOVO_SITES_PER_RANGE - 1;
    range->chrom_name_len = vb->chrom_name_len;
    range->chrom          = WORD_INDEX_NONE;  // Note: in REF_INTERNAL the chrom index is private to the VB prior to merge, so we can't use it
    
    // chrom_name points into vb->txt_data that will disappear at the end of this VB, 
    range->chrom_name = MALLOC (vb->chrom_name_len);
    memcpy ((char *)range->chrom_name, vb->chrom_name, vb->chrom_name_len);

    // nothing is set yet - in ZIP, bits get set as they are encountered in the compressing txt file
    bit_array_alloc (&range->ref, REF_NUM_DENOVO_SITES_PER_RANGE * 2);
    bit_array_alloc (&range->is_set, REF_NUM_DENOVO_SITES_PER_RANGE);
    bit_array_clear_region (&range->is_set, 0, REF_NUM_DENOVO_SITES_PER_RANGE);

    if (vb) {
        vb->prev_range = range;
        vb->prev_range_chrom_node_index = vb->chrom_node_index;
        vb->prev_range_range_i = range_i;
    }

    return range; // returning locked range
}

static Range *ref_seg_get_locked_range_loaded (VBlockP vb, PosType pos, uint32_t seq_len, const char *field /* used for ASSSEG */, RefLock *lock)  
{
    // case: we're asking for the same range as the previous one (for example, subsequent line in a sorted SAM)
    if (vb && vb->prev_range && vb->prev_range_chrom_node_index == vb->chrom_node_index) {
        PosType gpos = vb->prev_range->gpos + (pos - vb->prev_range->first_pos);
        *lock = ref_lock (gpos, seq_len);
        return vb->prev_range;
    }

    if (!flag_reading_reference) { // segging VCF or SAM with external reference

        // if the chrom is not in the reference and it is numeric only, attempt to change it eg "22"->"chr22"
        uint32_t num_contigs = ref_contigs_num_contigs();
        if (vb->chrom_node_index >= num_contigs) 
            vb->chrom_node_index = ref_alt_chroms_zip_get_alt_index (vb->chrom_name, vb->chrom_name_len, WI_REF_CONTIG, vb->chrom_node_index); // change temporarily just for ref_range_id_by_word_index()

        // case: the contig is not in the reference - we will just consider it an unaligned line 
        // (we already gave a warning for this in ref_contigs_verify_identical_chrom, so no need for another one)
        if (vb->chrom_node_index >= num_contigs) return NULL;
    }

    ASSSEG (vb->chrom_node_index < ranges.len, field, "Error in ref_seg_get_locked_range: vb->chrom_node_index=%u expected to be smaller than ranges.len=%u", 
            vb->chrom_node_index, (uint32_t)ranges.len);

    Range *range = ENT (Range, ranges, vb->chrom_node_index);
    PosType gpos = range->gpos + (pos - range->first_pos);

    *lock = ref_lock (gpos, seq_len);

    // when using an external refernce, pos has to be within the reference range
    // note: in SAM, if a read starts within the valid range, it is allowed to overflow beyond it - and we will circle
    // around to the beginning of the range assuming its a circular chromosome (see in sam_seg_seq_field)
    ASSSEG ((pos >= range->first_pos && pos <= range->last_pos), field,
            "Error: file has POS=%"PRId64" for contig \"%.*s\", but this contig's range is %"PRId64" - %"PRId64". Likely this is because %s was created using a reference file other than %s.",
            pos, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos, txt_name, ref_filename);

    return range; // returning locked range
}

// ZIP: returns a range that includes pos and is locked for (pos,len)
// case 1: ZIP: in SAM with REF_INTERNAL, when segging a SEQ field ahead of committing it to the reference
// case 2: ZIP: SAM and VCF with REF_EXTERNAL: when segging a SAM_SEQ or VCF_REFALT field
// if range is found, returns a locked range, and its the responsibility of the caller to unlock it. otherwise, returns NULL
Range *ref_seg_get_locked_range (VBlockP vb, PosType pos, uint32_t seq_len, const char *field /* used for ASSSEG */, RefLock *lock)  
{
    // sanity checks
    ASSSEG0 (vb->chrom_name, field, "Error in ref_seg_get_locked_range: vb->chrom_name=NULL");

    uint32_t save_chrom_index = vb->chrom_node_index; // might temporarily change with alt chrom names
    
    Range *range = NULL;
    switch (ranges.param) {
        case RT_DENOVO : range = ref_seg_get_locked_range_denovo (vb, pos, field, lock); break;
        case RT_LOADED : range = ref_seg_get_locked_range_loaded (vb, pos, seq_len, field, lock); break;
        default        : ABORT ("Error in ref_seg_get_locked_range: invalid ranges.param=%"PRId64, ranges.param);
    }

    vb->chrom_node_index = save_chrom_index; // restore
    return range; // returning locked range
}

// ----------------------------------------------
// Compressing ranges into SEC_REFERENCE sections
// ----------------------------------------------

static void ref_copy_one_compressed_section (File *ref_file, const RAEntry *ra, SectionListEntry **sl)
{
    // get section list entry from ref_file_section_list - which will be used by zfile_read_section to seek to the correct offset
    while (*sl < AFTERENT (SectionListEntry, ref_file_section_list) && 
           !((*sl)->vblock_i == ra->vblock_i && (*sl)->section_type == SEC_REFERENCE)) 
        (*sl)++;

    ASSERT (*sl < AFTERENT (SectionListEntry, ref_file_section_list), "Error in ref_copy_one_compressed_section: cannot find FASTA_NONREF of vb_i=%u in section list of reference file", ra->vblock_i);

    static Buffer ref_seq_section = EMPTY_BUFFER;

    RESET_FLAG (flag_show_headers);
    zfile_read_section (ref_file, evb, ra->vblock_i, &ref_seq_section, "ref_seq_section", 
                        sizeof (SectionHeaderReference), SEC_REFERENCE, *sl);
    RESTORE_FLAG (flag_show_headers);

    SectionHeaderReference *header = (SectionHeaderReference *)ref_seq_section.data;

    ASSERT0 (BGEN32 (header->chrom_word_index) == ra->chrom_index && BGEN64 (header->pos) == ra->min_pos,
            "Error in ref_copy_one_compressed_section: RA and Section don't agree on chrom or pos");

    // some minor changes to the header...
    header->h.vblock_i  = 0; // we don't belong to any VB and there is no encryption of external ref

    // Write header and body of the reference to z_file
    // Note on encryption: reference sections originating from an external reference are never encrypted - not
    // by us here, and not in the source reference fasta (because with disallow --make-reference in combination with --password)
    START_TIMER;
    file_write (z_file, ref_seq_section.data, ref_seq_section.len);
    COPY_TIMER_VB (evb, write);

    // "manually" add the reference section to the section list - normally it is added in comp_compress()
    sections_add_to_list (evb, &header->h);

    z_file->disk_so_far += ref_seq_section.len;   // length of GENOZIP data writen to disk

    if (flag_show_reference) {
        Context *ctx = &z_file->contexts[CHROM];
        MtfNode *node = ENT (MtfNode, ctx->mtf, BGEN32 (header->chrom_word_index));
        fprintf (stderr, "Copying SEC_REFERENCE from %s: chrom=%u (%s) gpos=%"PRId64" pos=%"PRId64" num_bases=%u section_size=%u\n", 
                 ref_filename, BGEN32 (header->chrom_word_index), 
                 ENT (char, ctx->dict, node->char_index), 
                 BGEN64(header->gpos), BGEN64(header->pos), 
                 BGEN32 (header->num_bases), 
                 BGEN32 (header->h.data_compressed_len) + BGEN32 (header->h.compressed_offset));
    }

    buf_free (&ref_seq_section);
}

// ZIP copying parts of external reference to fine - called by I/O thread from zip_write_global_area->ref_compress_ref
static void ref_copy_compressed_sections_from_reference_file (void)
{
    ASSERT (primary_command == ZIP && flag_reference == REF_EXT_STORE, 
            "Error in ref_copy_compressed_sections_from_reference_file: not expecting to be here: primary_command=%u flag_reference=%u", primary_command, flag_reference);

    File *ref_file = file_open (ref_filename, READ, Z_FILE, DT_FASTA);

    // note: in a FASTA file compressed with --make-reference, there is exactly one RA per VB (a contig or part of a contig)
    // and, since this is ZIP with EXT_STORE, also exactly one range per contig. We loop one RA at a time and:
    // 1. If 95% of the ref file RA is set in the zfile contig range - we copy the compressed reference section directly from the ref FASTA
    // 2. If we copied from the FASTA, we mark those region covered by the RA as "is_set=0", so that we don't compress it later
    SectionListEntry *sl = FIRSTENT (SectionListEntry, ref_file_section_list);
    ARRAY (RAEntry, fasta_sec, ref_external_ra);

    for (uint32_t i=0; i < ref_external_ra.len; i++) {

        Range *contig_r = ENT (Range, ranges, fasta_sec[i].chrom_index);
        PosType fasta_sec_start_in_contig_r = fasta_sec[i].min_pos - contig_r->first_pos; // the start of the FASTA section (a bit less than 1MB) within the full-contig range

        PosType fasta_sec_len = fasta_sec[i].max_pos - fasta_sec[i].min_pos + 1;
        PosType bits_is_set   = bit_array_num_bits_set_region (&contig_r->is_set, fasta_sec_start_in_contig_r, fasta_sec_len);

        // if this at least 95% of the RA is covered, just copy the corresponding FASTA section to our file, and
        // mark all the ranges as is_set=false indicating that they don't need to be compressed individually
        if ((double)bits_is_set / (double)fasta_sec_len >= 0.95) {
            ref_copy_one_compressed_section (ref_file, &fasta_sec[i], &sl);
            bit_array_clear_region (&contig_r->is_set, fasta_sec_start_in_contig_r, fasta_sec_len);
        }
    }

    file_close (&ref_file, false);
}

// remove the unused parts of a range and the beginning and end of the range, and update first/last_pos.
static bool ref_remove_flanking_regions (Range *r, uint64_t r_num_set_bits, uint64_t *start_flanking_region_len /* out */)
{
// threshold - if the number of clear bits (excluding flanking regions) is beneath this, we will not copmact, as the cost in
// z_file size of a SEC_REF_IS_SET section needed if compacting will be more than what we save in compression of r->ref
#define THRESHOLD_FOR_COMPACTING 470 

    uint64_t end_flanking_region_len, last_1;
    
    char has_any_bit = bit_array_find_first_set_bit (&r->is_set, start_flanking_region_len);
    ASSERT (has_any_bit, "Error in ref_remove_flanking_regions: range %u (%s) has no bits set in r->is_set", (uint32_t)(r-(Range*)ranges.data), r->chrom_name); // ref_prepare_range_for_compress is responsible not to send us 0-bit ranges

    has_any_bit = bit_array_find_prev_set_bit (&r->is_set, r->is_set.num_of_bits, &last_1);
    ASSERT (has_any_bit, "Error in ref_remove_flanking_regions: range %u (%s) has no bits set in r->is_set (#2)", (uint32_t)(r-(Range*)ranges.data), r->chrom_name); // this should definitely never happen, since we already know the range has bits
    end_flanking_region_len = r->is_set.num_of_bits - last_1 - 1;

    uint64_t num_clear_bits_excluding_flanking_regions = 
        r->is_set.num_of_bits - r_num_set_bits - *start_flanking_region_len - end_flanking_region_len;

    // remove flanking regions - will allow a smaller allocation for the reference in PIZ 
    r->first_pos += *start_flanking_region_len;
    r->last_pos  -= end_flanking_region_len;

    if (ranges.param == RT_LOADED) // note: for RT_DENOVO, gpos is still 0 at this point
        r->gpos  += *start_flanking_region_len;

    ASSERT (r->last_pos >= r->first_pos, "Error in ref_compact_ref: bad removal of flanking regions: r->first_pos=%"PRId64" r->last_pos=%"PRId64,
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
static bool ref_compact_ref (Range *r, uint64_t r_num_set_bits)
{
    if (!r || !r_num_set_bits) return false;

    ASSERT0 (r->is_set.num_of_bits, "Error in ref_compact_ref: r->is_set.num_of_bits=0");

    // remove flanking regions
    uint64_t start_flanking_region_len;
    bool is_compact_needed = ref_remove_flanking_regions (r, r_num_set_bits, &start_flanking_region_len);
    if (!is_compact_needed) return false;

    uint64_t start_1_offset=0, start_0_offset=0, compact_len=0;
    while (1) {
        
        // find length of set region
        bool has_any_bit = bit_array_find_next_clear_bit (&r->is_set, start_1_offset, &start_0_offset);
        uint64_t len_1 = (has_any_bit ? start_0_offset : r->is_set.num_of_bits) - start_1_offset;

        // do actual compacting - move set region to be directly after the previous set region (or at the begining if its the first)
        bit_array_copy (&r->ref, compact_len * 2, &r->ref, (start_flanking_region_len + start_1_offset) * 2, len_1 * 2);
        compact_len += len_1;

        if (!has_any_bit) break; // case: we're done- this 1 region goes to the end the range - there are no more clear regions

        // find length of clear region
        has_any_bit = bit_array_find_next_set_bit (&r->is_set, start_0_offset, &start_1_offset);
        ASSERT0 (has_any_bit, "Error in ref_compact_ref: cannot find set bits"); // this should never happen as we removed the flanking regions
    }

    // set length of ref - this is the data that will be compressed
    r->ref.num_of_bits  = compact_len * 2;
    r->ref.num_of_words = roundup_bits2words64 (r->ref.num_of_bits); 

    return true;
}

static void ref_compress_one_range (VBlockP vb)
{
    Range *r = vb->range; // will be NULL if we're being asked to write a final, empty section

    // remove flanking regions, and if beneficial also compact it further by removing unused nucleotides 
    bool is_compacted = flag_make_reference ? false : ref_compact_ref (r, vb->range_num_set_bits); // true if it is compacted beyong just the flanking regions

    SectionHeaderReference header = { .h.vblock_i          = BGEN32 (vb->vblock_i),
                                      .h.magic             = BGEN32 (GENOZIP_MAGIC),
                                      .h.compressed_offset = BGEN32 (sizeof(header)),
                                      .chrom_word_index    = r ? BGEN32 (r->chrom) : WORD_INDEX_NONE,
                                      .pos                 = r ? BGEN64 (r->first_pos) : 0,
                                      .gpos                = r ? BGEN64 (r->gpos) : 0 };

    vb->z_data.name  = "z_data"; // comp_compress requires that these are pre-set
    vb->z_data.param = vb->vblock_i;

    // First, SEC_REF_IS_SET section (but not needed if the entire range is used)
    if (r && is_compacted) {

        LTEN_bit_array (&r->is_set);

        header.h.section_type          = SEC_REF_IS_SET;  // most of the header is the same as ^
        header.h.codec                 = CODEC_BZ2;
        header.h.data_uncompressed_len = BGEN32 (r->is_set.num_of_words * sizeof (uint64_t));
        header.num_bases               = BGEN32 ((uint32_t)ref_size (r)); // full length, after flanking regions removed
        comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, (char *)r->is_set.words, NULL);

        if (flag_show_reference && r) 
            fprintf (stderr, "vb_i=%u Compressing SEC_REF_IS_SET chrom=%u (%.*s) gpos=%"PRIu64" pos=%"PRIu64" num_bases=%u section_size=%u bytes\n", 
                    vb->vblock_i, BGEN32 (header.chrom_word_index), r->chrom_name_len, r->chrom_name,
                    BGEN64 (header.gpos), BGEN64 (header.pos), BGEN32 (header.num_bases), 
                    BGEN32 (header.h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));
    }

    // Second. SEC_REFERENCE
    if (r) LTEN_bit_array (&r->ref);

    header.h.section_type          = SEC_REFERENCE;
    header.h.codec                 = CODEC_LZMA; // better than BSC: slightly better compression and compression speed, 2.5X faster decompression
    header.h.compressed_offset     = BGEN32 (sizeof(header)); // reset compressed offset - if we're encrypting - REF_IS_SET was encrypted and compressed_offset padded, by REFERENCE is never encrypted
    header.h.data_uncompressed_len = r ? BGEN32 (r->ref.num_of_words * sizeof (uint64_t)) : 0;
    header.num_bases               = r ? BGEN32 (r->ref.num_of_bits / 2) : 0; // less than ref_size(r) if compacted
    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, r ? (char *)r->ref.words : NULL, NULL);

    if (flag_show_reference && r) 
        fprintf (stderr, "vb_i=%u Compressing SEC_REFERENCE chrom=%u (%.*s) %s gpos=%"PRIu64" pos=%"PRIu64" num_bases=%u section_size=%u bytes\n", 
                 vb->vblock_i, BGEN32 (header.chrom_word_index), r->chrom_name_len, r->chrom_name, is_compacted ? "compacted " : "",
                 BGEN64 (header.gpos), BGEN64 (header.pos), BGEN32 (header.num_bases), 
                 BGEN32 (header.h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));

    // store the ref_stored_ra data for this range
    if (r) {
        spin_lock (ref_stored_ra_spin);
        NEXTENT (RAEntry, ref_stored_ra) = (RAEntry){ .vblock_i    = vb->vblock_i, 
                                                      .chrom_index = r->chrom,
                                                      .min_pos     = r->first_pos,
                                                      .max_pos     = r->last_pos };
        spin_unlock (ref_stored_ra_spin);
    }

    // insert this range sequence into the ref_hash (included in the reference file, for use to compress of FASTQ, unaligned SAM and FASTA)
    if (flag_make_reference)
        refhash_calc_one_range (r, ISLASTENT (ranges, r) ? NULL : r+1);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

void ref_output_vb (VBlockP vb) 
{ 
    zip_output_processed_vb (vb, &vb->section_list_buf, false, PD_REFERENCE_DATA); 
}

// compress the reference - one section at the time, using Dispatcher to do them in parallel 
// note: this is not called in make_reference - instead, ref_make_prepare_range_for_compress is called
static void ref_prepare_range_for_compress (VBlockP vb)
{
    static uint32_t next_range_i=0;
    if (vb->vblock_i == 1) next_range_i=0; // initialize once per non-bound file

    // find next occupied range
    for (; !vb->ready_to_dispatch && next_range_i < ranges.len ; next_range_i++) {

        Range *r = ENT (Range, ranges, next_range_i);

        uint64_t num_set = bit_array_num_bits_set (&r->is_set);
        if (!num_set) {
            r->is_set.num_of_bits = 0;
            continue; // nothing to with this range - perhaps copied and cleared in ref_copy_compressed_sections_from_reference_file
        }

        vb->range              = r; // range to compress
        vb->range_num_set_bits = num_set;
        vb->ready_to_dispatch  = true;
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
    ref_contigs_generate_data_if_denovo();

    // calculate all chrom indices
    for (uint32_t range_i=0; range_i < ranges.len; range_i++) {
        Range *r = ENT (Range, ranges, range_i);
        
        if (ref_is_range_used (r))
            r->chrom = ref_contigs_get_word_index (r->chrom_name, r->chrom_name_len, WI_REF_CONTIG, false);
    }

    // sort by chrom then pos, and place the unused ranges at the end
    qsort (ranges.data, ranges.len, sizeof (Range), ref_contigs_range_sorter);

    // shorten the array to only used ranges
    Range *r ; for (r=FIRSTENT (Range, ranges); r < AFTERENT (Range, ranges) && ref_is_range_used (r); r++) {};
    ranges.len = r - FIRSTENT (Range, ranges);
}

// ZIP: compress and write reference sections. either compressed by us, or copied directly from the reference file.
void ref_compress_ref (void)
{
    if (!buf_is_allocated (&ranges)) return;

    if ((ranges.param == RT_DENOVO) &&
        buf_is_allocated (&z_file->contexts[CHROM].dict)) // did we have an aligned lines? (to do: this test is not enough)
        ref_finalize_denovo_ranges(); // assignes chroms; sorts ranges by chrom, pos; gets rid of unused ranges

    if (ranges.param != RT_MAKE_REF) {
        ref_contigs_compress(); // also assigns gpos to de-novo ranges 
        zip_output_processed_vb (evb, &evb->section_list_buf, false, PD_REFERENCE_DATA); 
    }

    // copy already-compressed SEQ sections from the genozip reference file, but only such sections that are almost entirely
    // covered by ranges with is_accessed=true. we mark these ranges affected as is_accessed=false.
    if (ranges.param == RT_LOADED)
        ref_copy_compressed_sections_from_reference_file ();

    buf_alloc (evb, &ref_stored_ra, sizeof (RAEntry) * ranges.len, 1, "ref_stored_ra", 0);
    ref_stored_ra.len = 0; // re-initialize, in case we read the external reference into here
    
    spin_initialize (ref_stored_ra_spin);

    // compression of reference doesn't output % progress
    SAVE_FLAG (flag_quiet);
    if (flag_show_reference) flag_quiet = true; // show references instead of progress

    // proceed to compress all ranges that have still have data in them after copying
    uint32_t num_vbs_dispatched = 
        dispatcher_fan_out_task (NULL, PROGRESS_MESSAGE, "Writing reference...", false, 
                                 flag_make_reference ? ref_make_prepare_range_for_compress : ref_prepare_range_for_compress, 
                                 ref_compress_one_range, 
                                 ref_output_vb);

    RESTORE_FLAG (flag_quiet);
    
    // SAM require at least one reference section, but if the SAM is unaligned, there will be none - create one empty section
    // (this will also happen if SAM has just only reference section, we will just needlessly write another tiny section - no harm)
    // incidentally, this empty section will also be written in case of a small (one vb) reference - no harm
    if (z_file->data_type == DT_SAM && num_vbs_dispatched==1) {
        evb->range = NULL;
        ref_compress_one_range (evb); // written with vb_i=0, section header only (no body)
    }

    // compress reference random access (note: in case of a reference file, SEC_REF_RAND_ACC will be identical to SEC_RANDOM_ACCESS. That's ok, its tiny)
    random_access_finalize_entries (&ref_stored_ra); // sort in order of vb_i

    // range data needed for contigs is only set in ref_make_prepare_range_for_compress which happens as part of the dispatcher_fan_out_task
    if (ranges.param == RT_MAKE_REF) {
        // sort by chrom then pos - so later piz can binary-search by chrom index
        qsort (ranges.data, ranges.len, sizeof (Range), ref_contigs_range_sorter);

        ref_contigs_compress(); 
        zip_output_processed_vb (evb, &evb->section_list_buf, false, PD_REFERENCE_DATA); 
    }
}

// -------------------------------
// Loading an external reference
// -------------------------------

void ref_set_reference (const char *filename)
{
    ASSERT0 (filename, "Error in ref_set_reference: filename is NULL");

    ref_filename = MALLOC (strlen (filename) + 1);

    strcpy ((char*)ref_filename, filename);
}

// called when loading an external reference
void ref_set_ref_file_info (Md5Hash md5, const char *fasta_name)
{
    ref_md5 = md5;

    if (fasta_name[0]) {
        ref_fasta_name = MALLOC (strlen (fasta_name) + 1); 
        strcpy (ref_fasta_name, fasta_name);
    }
}

// display the reference:
// if --reference is used, normal reference is shown, if --REFERENCE - reverse complement is shown
// show a subset of the reference if --regions is specified - but only up to one region per chromosome
static void ref_display_ref (void)
{
    for (const Range *r = FIRSTENT (Range, ranges); r < AFTERENT (Range, ranges); r++) {

        PosType display_first_pos, display_last_pos;
        if (!regions_get_range_intersection (r->chrom, r->first_pos, r->last_pos, &display_first_pos, &display_last_pos))
            continue;

        printf ("%.*s\n", r->chrom_name_len, r->chrom_name);

        // case: normal sequence
        if (flag_reference == REF_EXTERNAL)
            for (PosType pos=display_first_pos; pos <= display_last_pos; pos++)
                fputc (ref_get_nucleotide (r, pos - r->first_pos), stdout);

        // case: reverse complement
        else
            for (PosType pos=display_last_pos; pos >= display_first_pos; pos--) {
                char base = ref_get_nucleotide (r, pos - r->first_pos);
                switch (base) {
                    case 'G': fputc ('C', stdout); break;
                    case 'C': fputc ('G', stdout); break;
                    case 'A': fputc ('T', stdout); break;
                    case 'T': fputc ('A', stdout); break;
                    case 'g': fputc ('c', stdout); break;
                    case 'c': fputc ('g', stdout); break;
                    case 'a': fputc ('t', stdout); break;
                    case 't': fputc ('a', stdout); break;
                    default : fputc (base, stdout); 
                }
            }        
        
        fputc ('\n', stdout);
    }
}

#define REV_CODEC_GENOME_BASES_PER_THREAD (1 << 27) // 128Mbp

static void ref_reverse_compliment_genome_prepare (VBlock *vb)
{
    vb->ready_to_dispatch = (vb->vblock_i-1) * REV_CODEC_GENOME_BASES_PER_THREAD < genome_size;
}

static void ref_reverse_compliment_genome_do (VBlock *vb)
{
    bit_array_reverse_complement_all (&genome_rev.ref, &genome.ref, 
                                      (vb->vblock_i-1) * REV_CODEC_GENOME_BASES_PER_THREAD, 
                                      REV_CODEC_GENOME_BASES_PER_THREAD);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

void ref_generate_reverse_complement_genome (void)
{
    START_TIMER;
    dispatcher_fan_out_task (NULL, PROGRESS_NONE, 0, false, 
                                ref_reverse_compliment_genome_prepare, 
                                ref_reverse_compliment_genome_do, 
                                NULL);
    COPY_TIMER_VB (evb, generate_rev_complement_genome);
}

static void ref_save_genome_copy_if_needed (bool is_last_file)
{
    if (flag_reference != REF_EXT_STORE) return; // relevant only for REF_EXT_STORE

    ASSERT0 (genome.ref.words, "Error in ref_save_genome_copy_if_needed: genome.ref is not allocated");

    // make a copy of genome and ranges
    if (!is_last_file && !buf_is_allocated (&genome_ref_copy)) {
        buf_alloc (evb, &genome_ref_copy, genome.ref.num_of_words * sizeof (word_t), 1, "genome_ref_copy", 0);
        BitArray *bar = buf_get_bitarray (&genome_ref_copy);
        bar->num_of_bits  = genome.ref.num_of_bits;
        bar->num_of_words = genome.ref.num_of_words;

        // copy initial state of ranges, before they are modifed by compacting
        buf_copy (evb, &ranges_copy, &ranges, sizeof (Range), 0, 0, "ranges_copy", 0);
    }
}

bool ref_is_reference_loaded (void)
{
    return buf_is_allocated (&ranges);
}

// ZIP & PIZ: import external reference
void ref_load_external_reference (bool display, bool is_last_file)
{
    ASSERT0 (ref_filename, "Error: ref_filename is NULL");

    flag_reading_reference = true; // tell file.c and fasta.c that this is a reference

    z_file = file_open (ref_filename, READ, Z_FILE, DT_FASTA);    
    z_file->basename = file_basename (ref_filename, false, "(reference)", NULL, 0);
        
    // save and reset flags that are intended to operate on the compressed file rather than the reference file
    RESET_FLAG (flag_test);
    RESET_FLAG (flag_unbind);
    RESET_FLAG (flag_md5);
    RESET_FLAG (flag_show_time);
    RESET_FLAG (flag_show_memory);
    RESET_FLAG (flag_show_stats);
    RESET_FLAG (flag_no_header);
    RESET_FLAG (flag_header_one);
    RESET_FLAG (flag_header_only);
    RESET_FLAG (flag_grep);
    RESET_FLAG (flag_regions);
    RESET_FLAG (flag_show_index);
    RESET_FLAG (flag_show_dict);
    RESET_FLAG (flag_show_b250);
    RESET_FLAG (flag_show_ref_contigs);
    RESET_FLAG (dict_id_show_one_b250);
    RESET_FLAG (dict_id_show_one_dict);
    RESET_FLAG (dump_one_b250_dict_id);
    RESET_FLAG (dump_one_local_dict_id);
    RESET_FLAG (flag_list_chroms);
    TEMP_FLAG (command, PIZ);

    bool piz_successful = piz_one_file (true, false);
    ASSERT (piz_successful, "Error: failed to uncompress reference file %s", ref_filename);

    // recover globals
    flag_reading_reference = false;
    RESTORE_FLAG (command);
    RESTORE_FLAG (flag_test);
    RESTORE_FLAG (flag_unbind);
    RESTORE_FLAG (flag_md5);
    RESTORE_FLAG (flag_show_time);
    RESTORE_FLAG (flag_show_memory);
    RESTORE_FLAG (flag_show_stats);
    RESTORE_FLAG (flag_no_header);
    RESTORE_FLAG (flag_header_one);
    RESTORE_FLAG (flag_header_only);
    RESTORE_FLAG (flag_grep);
    RESTORE_FLAG (flag_regions);
    RESTORE_FLAG (flag_show_index);
    RESTORE_FLAG (flag_show_dict);
    RESTORE_FLAG (flag_show_b250);
    RESTORE_FLAG (flag_show_ref_contigs);
    RESTORE_FLAG (dict_id_show_one_b250);
    RESTORE_FLAG (dict_id_show_one_dict);
    RESTORE_FLAG (dump_one_b250_dict_id);
    RESTORE_FLAG (dump_one_local_dict_id);
    RESTORE_FLAG (flag_list_chroms);

    file_close (&z_file, false);
    file_close (&txt_file, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtfile_genozip_to_txt_header.

    // If we're zipping a FASTQ or unaligned SAM, create a reverse-complement genome too
    if (flag_ref_use_aligner && primary_command == ZIP) 
        ref_generate_reverse_complement_genome();

    ref_save_genome_copy_if_needed (is_last_file);

    if (display) ref_display_ref();
}

static void ref_allocate_loaded_genome (void)
{
    // notes on ranges.len:
    // 1. in case stored reference originates from REF_INTERNAL, we have a contig for every chrom for which there was a sequence,
    // but there might be additional chroms in the SAM file (eg in RNEXT) that don't have a sequence. Since we're indexing ranges by chrom,
    // we need a range for each chrom, even if it doesn't have data 
    // 2. in case stored reference originates from REF_EXT_STORE, we have contigs for all the chroms in the original reference file,
    // but our CHROM contexts also includes alternate chrom names that aren't in the original reference that appear after the reference
    // chroms. we need to make sure ranges.len does not include alternate chroms as that's how we know a chrom is alternate in ref_piz_get_range
    // 3. in case loading from a reference file, the number of contigs will match the number of chroms, so no issues.
    ranges.len = (z_file->flags & GENOZIP_FL_REF_INTERNAL) ? z_file->contexts[CHROM].word_list.len :
                                                             ref_contigs_num_contigs();

    buf_alloc (evb, &ranges, ranges.len * sizeof (Range), 1, "ranges", RT_LOADED); 
    buf_zero (&ranges);
    
    for (uint32_t range_id=0; range_id < ranges.len; range_id++) 
        ENT (Range, ranges, range_id)->range_id = range_id;

    // note: genome_size must be full words, so that bit_array_reverse_complement doesn't need to shift
    genome_size = ROUNDUP64 (ref_contigs_get_genome_size()) + 64; // round up to the nearest 64, and add one word, needed by aligner_get_match_len for bit shifting overflow
    genome = (Range){ .chrom_name = CHROM_NAME_GENOME, 
                      .chrom_name_len = strlen (CHROM_NAME_GENOME),
                      .chrom = CHROM_GENOME, 
                      .last_pos = genome_size - 1 };

    bit_array_alloc (&genome.ref, genome_size * 2); // 2 bits per base
    bit_array_clear_all (&genome.ref); // make sure the gaps between the contigs are zero as "matches" can covered gaps as well

    // we don't need is_set if we're compressing with REF_EXTERNAL 
    bool has_is_set = primary_command == PIZ || (primary_command == ZIP && flag_reference == REF_EXT_STORE);
    if (has_is_set) {
        bit_array_alloc (&genome.is_set, genome_size); 
        bit_array_clear_all (&genome.is_set); 
    }
    
    genome_rev = (Range){ .chrom_name = CHROM_NAME_GENOME_REV, 
                          .chrom_name_len = strlen (CHROM_NAME_GENOME_REV),
                          .chrom = CHROM_GENOME_REV, 
                          .last_pos = genome_size - 1 };
    bit_array_alloc (&genome_rev.ref, genome_size * 2); 
    bit_array_clear_all (&genome_rev.ref); // make sure the gaps between the contigs are zero

    // we protect genome->ref while uncompressing reference data, and genome->is_set while segging
    ref_lock_initialize_loaded_genome();

    // overlay all chromosomes (range[i] goes to chrom_index=i) - note some chroms might not have a contig in 
    // which case their range is not initialized
    for (Range *r = FIRSTENT (Range, ranges) ; r < AFTERENT (Range, ranges); r++) {
        r->chrom = r - FIRSTENT (Range, ranges);
        const RefContig *rc = ref_contigs_get_contig (r->chrom, true);

        if (rc) { // this chromosome has reference data 
            r->gpos      = rc->gpos;
            r->first_pos = rc->min_pos;
            r->last_pos  = rc->max_pos;
            ref_contigs_get_chrom_snip (r->chrom, &r->chrom_name, &r->chrom_name_len);
            
            PosType len = rc->max_pos - rc->min_pos + 1;

            ASSERT (r->gpos + len <= genome.ref.num_of_bits / 2, "Error in ref_allocate_loaded_genome: adding range \"%s\": r->gpos(%"PRId64") + len(%"PRId64") (=%"PRId64") is beyond genome_size=%"PRId64,
                    r->chrom_name, r->gpos, len, r->gpos+len, genome_size);

            bit_array_overlay (&r->ref, &genome.ref, r->gpos*2, len*2);

            if (has_is_set) 
                bit_array_overlay (&r->is_set, &genome.is_set, r->gpos, len);
        }
    }
}

 
// Initialize the overlay ranges, one per chromosome for a loaded genome. The remaining fields will be set in ref_read_one_range 
static void ref_create_contig_ranges_for_loaded_genome (void)
{
    for (unsigned chrom_index=0; chrom_index < ranges.len; chrom_index++) {
        Range *r = ENT (Range, ranges, chrom_index); // ranges is indexed by chrom_index

        r->chrom = chrom_index;      
        
        if (flag_reference == REF_STORED) {
            Context *ctx = &z_file->contexts[CHROM];
            mtf_get_snip_by_word_index (&ctx->word_list, &ctx->dict, r->chrom, &r->chrom_name, &r->chrom_name_len);
        }
        else
            ref_contigs_get_chrom_snip (r->chrom, &r->chrom_name, &r->chrom_name_len);

        // case: this chrom has no RA - for example, it is a "=" or "*", it doesn't have a reference region
        if (r->last_pos < 0) continue;
    }
}

// case 1: in case of ZIP with external reference, called by ref_load_stored_reference during piz_read_global_area of the reference file
// case 2: in case of PIZ: also called from ref_load_stored_reference with RT_LOADED
// case 3: in case of ZIP of SAM using internal reference - called from sam_zip_initialize
// note: ranges allocation must be called by the I/O thread as it adds a buffer to evb buf_list
void ref_initialize_ranges (RangesType ranges_type)
{
    if (ranges_type == RT_LOADED) {
        ref_allocate_loaded_genome();
        ref_create_contig_ranges_for_loaded_genome();
    }
    else { // RT_DENOVO
        if (buf_is_allocated (&ranges)) return; // case: 2nd+ bound file

        ranges.len = REF_NUM_DENOVO_RANGES;
    
        buf_alloc (evb, &ranges, ranges.len * sizeof (Range), 1, "ranges", RT_DENOVO); 
        buf_zero (&ranges);
        ref_lock_initialize_denovo_genome();
    }
}

typedef struct { PosType min_pos, max_pos; } MinMax;

//---------------------------------------
// Printing
//---------------------------------------

void ref_print_subrange (const char *msg, const Range *r, PosType start_pos, PosType end_pos) /* start_pos=end_pos=0 if entire ref */
{
    uint64_t start_idx = start_pos ? start_pos - r->first_pos : 0;
    uint64_t end_idx   = (end_pos ? MIN (end_pos, r->last_pos) : r->last_pos) - r->first_pos;

    fprintf (stderr, "%s: %.*s %"PRId64" - %"PRId64" (len=%u): ", msg, r->chrom_name_len, r->chrom_name, start_pos, end_pos, (uint32_t)(end_pos - start_pos + 1));
    for (uint64_t idx = start_idx; idx <= end_idx; idx++) 
        fputc (ref_get_nucleotide (r, idx) + (32 * !ref_is_nucleotide_set (r, idx)), stderr); // uppercase if set, lowercase if not

    fputc ('\n', stderr);
}

void ref_print_is_set (const Range *r,
                       PosType around_pos) // display around this neighborhoud ; -1 means entire range
{
#   define neighborhood (PosType)10000

    fprintf (stderr, "\n\nRegions set for chrom %u \"%.*s\" [%"PRId64"-%"PRId64"] according to range.is_set (format- \"first_pos-last_pos (len)\")\n", 
             r->chrom, r->chrom_name_len, r->chrom_name, r->first_pos, r->last_pos);
    fprintf (stderr, "In the neighborhood of about %u bp around pos=%"PRId64"\n", (unsigned)neighborhood, around_pos);

    if (!r->is_set.num_of_bits)
        fprintf (stderr, "No data: r->is_set.num_of_bits=0\n");

    if (around_pos < r->first_pos || around_pos > r->last_pos)
        fprintf (stderr, "No data: pos=%"PRId64" is outside of [first_pos=%"PRId64" - last_pos=%"PRId64"]\n", around_pos, r->first_pos, r->last_pos);

    uint64_t next;
    for (uint64_t i=0; i < r->is_set.num_of_bits; ) {
        
        bool found = bit_array_find_next_clear_bit (&r->is_set, i, &next);
        if (!found) next = r->is_set.num_of_bits;

        bool in_neighborhood = (around_pos - (PosType)(r->first_pos+i) > -neighborhood) && (around_pos - (PosType)(r->first_pos+i) < neighborhood);
        if (next > i && (around_pos == -1 || in_neighborhood)) {
            if (next - i > 1)
                fprintf (stderr, "%"PRId64"-%"PRIu64"(%u)\t", r->first_pos + i, r->first_pos + next-1, (uint32_t)(next - i));
            else
                fprintf (stderr, "%"PRId64"(1)\t", r->first_pos + i);
        }                   
        if (!found) break;

        i = next;

        found = bit_array_find_next_set_bit (&r->is_set, i, &next);
        if (!found) next = r->is_set.num_of_bits;

        i = next;
    }
    fprintf (stderr, "\n");
}

// returns the reference file name for CRAM, derived from the genozip reference name
const char *ref_get_cram_ref (void)
{
    static char *samtools_T_option = NULL;
    if (samtools_T_option) goto done; // already calculated

    ASSINP (ref_filename, "%s: when compressing a CRAM file, --reference or --REFERENCE must be specified", global_cmd);

    ASSINP (ref_fasta_name, "%s: cannot compress a CRAM file because %s is lacking the name of the source fasta file - likely because it was created by piping a fasta from from stdin, or because the name of the fasta provided exceed %u characters",
            global_cmd, ref_filename, REF_FILENAME_LEN-1);

    bool file_exists = (access (ref_fasta_name, F_OK) == 0);
    ASSINP (file_exists, "%s: cannot find the fasta file %s. Note: this is the file that was used to create %s, and it needs to exist in this name, in order to be passed to samtools as a reference file (-T option) for reading the CRAM file", 
            global_cmd, ref_fasta_name, ref_filename);

    samtools_T_option = MALLOC (strlen (ref_fasta_name) + 10);
    sprintf (samtools_T_option, "-T%s", ref_fasta_name);

done:
    return samtools_T_option;
}
