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
#include "arch.h"
#include "sections.h"
#include "profiler.h"
#include "bit_array.h"
#include "file.h"
#include "regions.h"

// ZIP and PIZ, internal or external reference ranges. If in ZIP-INTERNAL we have REF_NUM_RANGES Range's - each allocated on demand. In all other cases we have one range per contig.
static Buffer ranges                    = EMPTY_BUFFER; 

// data derived from the external reference FASTA's CONTIG context
static Buffer contig_dict               = EMPTY_BUFFER;
static Buffer contig_words              = EMPTY_BUFFER;
static Buffer contig_maxes              = EMPTY_BUFFER; // an array of int64 - the max_pos of each contig
static Buffer contig_words_sorted_index = EMPTY_BUFFER; // an array of uint32 of indexes into contig_words - sorted by alphabetical order of the snip in contig_dict
static Buffer ref_external_ra           = EMPTY_BUFFER; // Random Access data of the external reference file
static Buffer ref_file_section_list     = EMPTY_BUFFER; // Section List of the external reference file

// PIZ: an array of RegionToSet - list of non-compacted regions, for which we should set is_set bit to 1
static Buffer region_to_set_list        = EMPTY_BUFFER; 
static pthread_mutex_t region_to_set_list_mutex;
static bool region_to_set_list_mutex_initialized=false;

// ZIP/PIZ: random_access data for the reference sections stored in a target genozip file
Buffer ref_stored_ra                    = EMPTY_BUFFER;
static pthread_mutex_t ref_stored_ra_mutex; // ZIP only
static bool ref_stored_ra_mutex_initialized = false;

typedef struct {
    BitArray *is_set;
    int64_t first_bit, len;
} RegionToSet;

static SectionListEntry *sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
static uint32_t ref_range_cursor = 0;

// globals
const char *ref_filename = NULL; // filename of external reference file
Md5Hash ref_md5 = MD5HASH_NONE;

#define SAVE_FLAG(flag) int save_##flag = flag ; flag=0
#define RESTORE_FLAG(flag) flag = save_##flag

// forward declarations
static int32_t ref_get_index_of_chrom_with_alt_name (VBlockP vb); 
static void ref_copy_chrom_data_from_z_file (void);

static inline int64_t ref_size (const Range *r) { return r ? (r->last_pos - r->first_pos + 1) : 0; }

// called by mtf_copy_reference_contig_to_chrom_ctx when initializing ZIP for a new file using a pre-loaded external reference
void ref_get_contigs (ConstBufferP *out_contig_dict, ConstBufferP *out_contig_words)
{
    *out_contig_dict  = &contig_dict;
    *out_contig_words = &contig_words;
}

// free memory allocations between files, when compressing multiple non-bound files or decompressing multiple files
void ref_cleanup_memory(void)
{
    if (flag_reference == REF_EXTERNAL) {
        if (primary_command == ZIP) {
            for (unsigned i=0; i < ranges.len ; i++) 
                bit_array_clear_all (&ENT (Range, ranges, i)->is_set); // zip sets 1 to nucleotides used. we clear it for the next file.
        }
        else {
            // PIZ we never cleanup external references as they can never change (the user can only specify one --reference)
        }
    }
    // case ZIP: REF_INTERNAL - we're done with ranges - free so next non-bound file can build its own
    // case ZIP: REF_EXT_STORE - unfortunately we will need to re-read the FASTA as we damaged it with removing flanking regions
    // case PIZ: REF_STORED - cleanup for next file with its own reference data

    else {
        if (buf_is_allocated (&ranges)) {
            ARRAY (Range, rng, ranges);
            for (unsigned i=0; i < ranges.len ; i++) {
                if (rng[i].ref.words)    FREE (rng[i].ref.words);
                if (rng[i].is_set.words) FREE (rng[i].is_set.words);
                if (primary_command == ZIP && flag_reference == REF_INTERNAL) 
                    FREE ((char *)rng[i].chrom_name); // allocated only in ZIP/REF_INTERNAL - otherwise a pointer into another Buffer
            }
            buf_free (&ranges);
        }
        buf_free (&region_to_set_list);

        buf_free (&contig_dict);          
        buf_free (&contig_words);
        buf_free (&contig_maxes);
        buf_free (&contig_words_sorted_index);
        buf_free (&ref_external_ra);
        buf_free (&ref_stored_ra);
        buf_free (&ref_file_section_list);
    }
}

MemStats ref_memory_consumption (void)
{
    MemStats stats = { "reference", 0, 0 };

    ARRAY (Range, r, ranges);
    if (buf_is_allocated (&ranges)) {
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
const Range *ref_piz_get_range (VBlockP vb, int64_t first_pos_needed, uint32_t num_nucleotides_needed)
{
    // caching
    if (vb->prev_range && vb->prev_range_chrom_node_index == vb->chrom_node_index)
        return vb->prev_range;

    // if the chrom is not in the reference and it is numeric only, attempt to change it eg "22"->"chr22"
    uint32_t index = vb->chrom_node_index < ranges.len ? vb->chrom_node_index : ref_get_index_of_chrom_with_alt_name (vb);

    ASSERT (index < ranges.len, "Error in ref_piz_get_range: vb->chrom_node_index=%u is out of range: ranges.len=%u",
            index, (uint32_t)ranges.len);

    Range *r = ENT (Range, ranges, index);
    if (!r->ref.num_of_words) return NULL; // this can ligitimately happen if entire chromosome is verbatim in SAM, eg. unaligned (pos=4) or SEQ or CIGAR are unavailable

    ASSERT (first_pos_needed + num_nucleotides_needed - 1 <= r->last_pos, "Error in ref_piz_get_range: out of range reconstructing txt_line_i=%u: pos=%"PRId64" num_nucleotides_needed=%u but range->last_pos=%"PRId64,
            vb->line_i, first_pos_needed, num_nucleotides_needed, r->last_pos);

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

        // do actual compacting - move set region to be directly after the previous set region (or at the begining if its the first)
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

// entry point of compute thread of reference decompression. this is called when pizzing a file with a stored reference.
// vb->z_data contains a SEC_REFERENCE section and sometimes also a SEC_REF_IS_SET sections
static void ref_uncompress_one_range (VBlockP vb)
{
    if (!buf_is_allocated (&vb->z_data)) goto finish; // we have no data in this VB because it was skipped due to --regions

    SectionHeaderReference *header = (SectionHeaderReference *)vb->z_data.data;

    int32_t chrom             = (int32_t)BGEN32 (header->chrom_word_index);
    uint32_t uncomp_len       = BGEN32 (header->h.data_uncompressed_len);
    int64_t ref_sec_first_pos = (int64_t)BGEN64 (header->first_pos);
    int64_t ref_sec_last_pos  = (int64_t)BGEN64 (header->last_pos);
    int64_t ref_sec_len       = ref_sec_last_pos - ref_sec_first_pos + 1;
    int64_t compacted_last_pos=0, compacted_ref_len=0, initial_flanking_len=0, final_flanking_len=0; 

    Context *ctx = &z_file->contexts[DTFZ(chrom)];
    ASSERT (chrom >= 0 && chrom < ctx->word_list.len, "Error in ref_uncompress_one_range: chrom=%d out of range - ctx->word_list.len=%u",
            chrom, (uint32_t)ctx->word_list.len);

    const char *chrom_name = mtf_get_snip_by_word_index (&ctx->word_list, &ctx->dict, chrom, 0, 0);
    Range *r = ENT (Range, ranges, chrom); // in PIZ, we have one range per chrom

    // if we can't find the chrom in the file data, try a reduced name "chr22" -> "22"
    int32_t alt_chrom = ctx->word_list.len;
    if (r->last_pos == RA_MISSING_RA_MAX && !strncmp (chrom_name, "chr", 3)) {

        // ineffecient search. consider adding a sorter to the chrom word_list and binary searching
        for (alt_chrom=0 ; alt_chrom < ctx->word_list.len; alt_chrom++) {
            char *alt_chrom_name = ENT (char, ctx->dict, ENT (MtfWord, ctx->word_list, alt_chrom)->char_index);
            if (!strcmp (&chrom_name[3], alt_chrom_name)) 
                r = ENT (Range, ranges, alt_chrom);
        }
    }

    ASSERT (r->last_pos != RA_MISSING_RA_MAX, "Error in ref_uncompress_one_range #2: chrom=%d \"%s\" appears as a reference section, but it doesn't appear in the file data",
            chrom, chrom_name);

    int64_t sec_start_within_contig = ref_sec_first_pos - r->first_pos;
    int64_t sec_end_within_contig   = ref_sec_last_pos  - r->first_pos;

    bool is_compacted = (header->h.section_type == SEC_REF_IS_SET); // we have a SEC_REF_IS_SET if  SEC_REFERENCE was compacted

    if (flag_show_reference && primary_command == PIZ && r)  // in ZIP, we show the compression of SEC_REFERENCE into z_file, not the uncompression of the reference file
        fprintf (stderr, "Uncompressing %-14s chrom=%u (%.*s) %"PRId64" - %"PRId64" (size=%"PRId64"). Bytes=%u\n", 
                st_name (header->h.section_type), BGEN32 (header->chrom_word_index), r->chrom_name_len, r->chrom_name, BGEN64 (header->first_pos), BGEN64 (header->last_pos), BGEN64 (header->last_pos) - BGEN64 (header->first_pos) + 1, 
                BGEN32 (header->h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));

    // initialization of is_set:
    // case 1: ZIP (reading an external reference) - we CLEAR is_set, and let seg set the bits that are to be
    //    needed from the reference for pizzing (so we don't store the others)
    //    note: in case of ZIP with REF_INTERNAL, we CLEAR the bits in ref_zip_get_locked_range
    // case 2: PIZ, reading an uncompacted (i.e. complete) reference section - which is always the case when
    //    reading an external reference and sometimes when reading an stored one - we SET all the bits as they are valid for pizzing
    //    we do this in ref_load_stored_reference AFTER all the SEC_REFERENCE/SEC_REF_IS_SEC sections are uncompressed,
    //    so that in case this was an REF_EXT_STORE compression, we first copy the contig-wide IS_SET sections (case 3) 
    //    (which will have 0s in the place of copied FASTA sections), and only after do we set these regions to 1.
    //    note: in case of PIZ, entire contigs are initialized to clear in ref_initialize_ranges as there might be
    //    regions missing (not covered by SEC_REFERENCE sections)
    // case 3: PIZ, reading a compacted reference (this only happens with stored references that original from 
    //    internal references) - we receive the correct is_set in the SEC_REF_IS_SET section and don't change it


    // case: if compacted, this SEC_REF_IS_SET sections contains r->is_set and its first/last_pos contain the coordinates
    // of the range, while the following SEC_REFERENCE section contains only the bases for which is_set is 1, 
    // first_pos=0 and last_pos=(num_1_bits_in_is_set-1)
    if (is_compacted) {

        // if compacted, the section must be within the boundaries of the contig (this is not true if the section was copied with ref_copy_one_compressed_section)
        ASSERT (sec_start_within_contig >= 0 && ref_sec_last_pos <= r->last_pos, 
                "Error in ref_uncompress_one_range: section range out of bounds for chrom=%d \"%s\": in SEC_REFERENCE being uncompressed: first_pos=%"PRId64" last_pos=%"PRId64" but in reference contig as initialized: first_pos=%"PRId64" last_pos=%"PRId64,
                chrom, chrom_name, ref_sec_first_pos, ref_sec_last_pos, r->first_pos, r->last_pos);

        ASSERT (uncomp_len == roundup_bits2words64 (ref_sec_len) * sizeof (uint64_t), "Error in ref_uncompress_one_range: when uncompressing SEC_REF_IS_SET: uncomp_len=%u inconsistent with len=%"PRId64, uncomp_len, ref_sec_len); 

        // uncompress into r->is_set, via vb->compress
        zfile_uncompress_section (vb, (SectionHeaderP)header, &vb->compressed, "compressed", SEC_REF_IS_SET);

        BitArray *is_set = buf_zfile_buf_to_bitarray (&vb->compressed, ref_sec_len);

        mutex_lock (r->mutex); // while different threads uncompress regions of the range that are non-overlapping, they might overlap at the bit level
        bit_array_copy (&r->is_set, sec_start_within_contig, is_set, 0, ref_sec_len); // initialization of is_set - case 3
        mutex_unlock (r->mutex);

        if (flag_show_is_set && !strcmp (chrom_name, flag_show_is_set)) ref_print_is_set (r);

        // prepare for uncompressing the next section - which is the SEC_REFERENCE
        header = (SectionHeaderReference *)&vb->z_data.data[*ENT (unsigned, vb->z_section_headers, 1)];

        if (flag_show_reference && primary_command == PIZ && r) 
            fprintf (stderr, "Uncompressing %-14s chrom=%u (%.*s) %"PRId64" - %"PRId64" (size=%"PRId64"). Bytes=%u\n", 
                    st_name (header->h.section_type), BGEN32 (header->chrom_word_index), r->chrom_name_len, r->chrom_name, BGEN64 (header->first_pos), BGEN64 (header->last_pos), BGEN64 (header->last_pos) - BGEN64 (header->first_pos) + 1, 
                    BGEN32 (header->h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));

        compacted_last_pos = BGEN64 (header->last_pos);
        compacted_ref_len  = compacted_last_pos - ref_sec_first_pos + 1;
        uncomp_len         = BGEN32 (header->h.data_uncompressed_len);

        ASSERT (uncomp_len == roundup_bits2words64 (compacted_ref_len*2) * sizeof (uint64_t), 
                "Error: uncomp_len=%u inconsistent with compacted_ref_len=%"PRId64, uncomp_len, compacted_ref_len); 

        ASSERT0 (BGEN32 (header->chrom_word_index) == chrom && BGEN64 (header->first_pos) == ref_sec_first_pos, // chrom should be the same between the two sections
                "Error in ref_uncompress_one_range: header mismatch between SEC_REF_IS_SET and SEC_REFERENCE sections");
    }
    
    // case: not compacted means that entire range is set
    else {
        ASSERT (uncomp_len == roundup_bits2words64 (ref_sec_len*2) * sizeof (uint64_t), "Error: uncomp_len=%u inconsistent with ref_len=%"PRId64, uncomp_len, ref_sec_len); 

        if (primary_command == ZIP) { // initialization of is_set - case 1
            mutex_lock (r->mutex); // while different threads uncompress regions of the range that are non-overlapping, they might overlap at the bit level
            bit_array_clear_region (&r->is_set, sec_start_within_contig, ref_sec_len); // entire range is cleared
            mutex_unlock (r->mutex);
        }

        else if (primary_command == PIZ) { // initialization of is_set - case 2

            // it is possible that the section goes beyond the boundaries of the contig, this can happen when we compressed with --REFERENCE
            // and the section was copied in its entirety from the reference FASTA file (in ref_copy_one_compressed_section)
            // in this case, we copy from the section only the part needed
            // IS THIS THIS REALLY POSSIBLE? NOT SURE. anyway, no harm in this code. divon 5/7/2020
            initial_flanking_len = (sec_start_within_contig < 0)    ? -sec_start_within_contig       : 0; // nucleotides in the section that are before the start of our contig
            final_flanking_len   = (ref_sec_last_pos > r->last_pos) ? ref_sec_last_pos - r->last_pos : 0; // nucleotides in the section that are after the end of our contig

            bit_array_set_region (&r->is_set, MAX (sec_start_within_contig, 0), ref_sec_len - initial_flanking_len - final_flanking_len);
            // save the region we need to set, we will do the actual setting in ref_load_stored_reference
            mutex_lock (region_to_set_list_mutex);
            RegionToSet *rts = &NEXTENT (RegionToSet, region_to_set_list);
            mutex_unlock (region_to_set_list_mutex);
            rts->is_set    = &r->is_set;
            rts->first_bit = MAX (sec_start_within_contig, 0);
            rts->len       = ref_sec_len - initial_flanking_len - final_flanking_len;
        }
   
        if (!uncomp_len) return;  // empty header - if it appears, it is the final header (eg in case of an unaligned SAM file)
    }

    // uncompress into r->ref, via vb->compress
    zfile_uncompress_section (vb, (SectionHeaderP)header, &vb->compressed, "compressed", SEC_REFERENCE);

    mutex_lock (r->mutex); // while different threads uncompress regions of the range that are non-overlapping, they might overlap at the bit level

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

    mutex_unlock (r->mutex);

finish:
    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. 
}

static void ref_read_one_range (VBlockP vb)
{
    if (!sections_get_next_section_of_type (&sl_ent, &ref_range_cursor, SEC_REFERENCE, SEC_REF_IS_SET))
        return; // no more reference sections

    // if the user specified --regions, check if this ref range is needed
    bool range_is_included = true;
    RAEntry *ra = NULL;
    if (flag_regions) {
        ASSERT (vb->vblock_i <= ref_stored_ra.len, "Error in ref_read_one_range: expecting vb->vblock_i(%u) <= ref_stored_ra.len(%u)", vb->vblock_i, (uint32_t)ref_stored_ra.len);

        ra = ENT (RAEntry, ref_stored_ra, vb->vblock_i-1);
        range_is_included = regions_is_ra_included (ra);
    }

    if (range_is_included) { 

        buf_alloc (vb, &vb->z_section_headers, 2 * sizeof(int32_t), 0, "z_section_headers", 2); // room for 2 section headers  

        ASSERT0 (vb->z_section_headers.len < 2, "Error in ref_read_one_range: unexpected 3rd recursive entry");

        int32_t section_offset = 
            zfile_read_section (z_file, vb, sl_ent->vblock_i, NO_SB_I, &vb->z_data, "z_data", 
                                sizeof(SectionHeaderReference), sl_ent->section_type, sl_ent);    

        ASSERT (*LASTENT(int32_t, vb->z_section_headers) != EOF, "Error in ref_read_one_range: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);

        NEXTENT (int32_t, vb->z_section_headers) = section_offset;

        // allocate memory for entire chrom reference if this is the first range of this chrom
        SectionHeaderReference *header = (SectionHeaderReference *)&vb->z_data.data[section_offset];
        
        int32_t chrom = BGEN32 (header->chrom_word_index);
        if (chrom == NIL) return; // we're done - terminating empty section that sometimes appears (eg in unaligned SAM that don't have any reference and yet are REF_INTERNAL)

        ASSERT (!flag_regions || (ra->chrom_index == chrom && ra->min_pos == BGEN64 (header->first_pos)), // note we don't compare max_pos because in the header it will be different if ref range is compacted
                "Error in ref_read_one_range: inconsistency between RA and ref header: RA=(chrom=%d min=%"PRId64") header==(chrom=%d min=%"PRId64")",
                ra->chrom_index, ra->min_pos, chrom, BGEN64 (header->first_pos));

        // if this is SEC_REF_IS_SET, read the SEC_REFERENCE section now
        if (header->h.section_type == SEC_REF_IS_SET) 
            ref_read_one_range (vb);
    }

    vb->ready_to_dispatch = true; // to simplify the code, we will dispatch the thread even if we skip the data, but we will return immediately. 
}

// PIZ: loading a reference stored in the genozip file - this could have been originally stored as REF_INTERNAL or REF_EXT_STORE
// or this could be a reference fasta file
void ref_load_stored_reference (void)
{
    ASSERT0 (!buf_is_allocated (&ranges), "Error in ref_load_stored_reference: expecting ranges to be unallocated");
    
    ref_initialize_ranges (true);
    
    sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
    ref_range_cursor = 0;

    mutex_initialize (&region_to_set_list_mutex, &region_to_set_list_mutex_initialized);
    buf_alloc (evb, &region_to_set_list, sections_count_sections (SEC_REFERENCE) * sizeof (RegionToSet), 1, "region_to_set_list", 0);

    // decompress reference using Dispatcher
    dispatcher_fan_out_task (flag_reference == REF_EXTERNAL ? ref_filename : z_file->basename, 
                             flag_reference == REF_EXTERNAL ? "Reading reference file..." : "Reading stored reference", flag_test, 
                             ref_read_one_range, 
                             ref_uncompress_one_range, 
                             NULL);

    // now we can safely set the is_set regions originating from non-compacted ranges. we couldn't do it before, because
    // copied-from-FASTA ranges appear first in the genozip file, and after them could be compacted ranges that originate
    // from a full-contig range in EXT_STORE, whose regions copied-from-FASTA are 0s.
    ARRAY (RegionToSet, rts, region_to_set_list);
    for (uint32_t i=0; i < region_to_set_list.len; i++)
        bit_array_set_region (rts[i].is_set, rts[i].first_bit, rts[i].len);

    buf_test_overflows_all_vbs ("ref_load_stored_reference");

    // exit now if all we wanted was just to see the reference
    if ((flag_show_reference || flag_show_is_set) && exe_type == EXE_GENOCAT) exit(0);
}

// ------------------------------------
// ZIP side
// ------------------------------------

// case REF_INTERNAL - we get a range_id by hashing the chrom name and the range_i - we can't use the chrom_node_index
// because we have multiple threads in parallel that might have the same node_index for different chroms
static inline uint32_t ref_range_id_by_hash (VBlockP vb, uint32_t range_i)
{
    ASSERT0 (vb->chrom_name_len > 0, "Error in ref_range_id_by_hash: vb->chrom_name_len==0");

    // step 1: get number embedded in the chrom name
    uint32_t n=0;
    for (unsigned i=0; i < vb->chrom_name_len; i++)
        if (IS_DIGIT (vb->chrom_name[i])) 
            n = n*10 + (vb->chrom_name[i] - '0');

    // if there are no digits, take the last 4 characters
    if (!n)
        n = (                           (((uint32_t)vb->chrom_name[vb->chrom_name_len-1]) ^ 0x1f))            ^ 
            (vb->chrom_name_len >= 2 ? ((((uint32_t)vb->chrom_name[vb->chrom_name_len-2]) ^ 0x1f) << 3)  : 0) ^
            (vb->chrom_name_len >= 3 ? ((((uint32_t)vb->chrom_name[vb->chrom_name_len-3]) ^ 0x1f) << 4)  : 0) ^ 
            (vb->chrom_name_len >= 4 ? ((((uint32_t)vb->chrom_name[vb->chrom_name_len-4]) ^ 0x1f) << 5)  : 0) ;

    // step: calculate the hash - 10 bit for n, 10 bit for range_i
    uint32_t chr_component = (vb->chrom_name_len <= 6) ? (n & 0x1f) : ((n % 896) + 128); // first 128 entries are reserved for the major chromosomes, heuristically identified as name length 6 or less

    uint32_t value = chr_component | (((range_i & 0x3ff) ^ ((n*3) & 0x3ff)) << 10);
    return value; 
}

// binary search for this chrom in contig_sorted_words. we count on gcc tail recursion optimization to keep this fast.
static int32_t ref_get_contig_word_index_do (const char *chrom_name, unsigned chrom_name_len, 
                                              int32_t first_sorted_index, int32_t last_sorted_index)
{
    if (first_sorted_index > last_sorted_index) return NIL; // not found

    uint32_t mid_sorted_index = (first_sorted_index + last_sorted_index) / 2;
    uint32_t word_index = *ENT (uint32_t, contig_words_sorted_index, mid_sorted_index);
    MtfWord *mid_word = ENT (MtfWord, contig_words, word_index);
    const char *snip = ENT (const char, contig_dict, mid_word->char_index);

    int cmp = strncmp (snip, chrom_name, chrom_name_len);
    if (!cmp && mid_word->snip_len != chrom_name_len) // identical prefix but different length
        cmp = mid_word->snip_len - chrom_name_len;

    if (cmp < 0) return ref_get_contig_word_index_do (chrom_name, chrom_name_len, mid_sorted_index+1, last_sorted_index);
    if (cmp > 0) return ref_get_contig_word_index_do (chrom_name, chrom_name_len, first_sorted_index, mid_sorted_index-1);

    return (int32_t)word_index;
}             

// binary search for this chrom in contig_sorted_words. we count on gcc tail recursion optimization to keep this fast.
static int32_t ref_get_contig_word_index (const char *chrom_name, unsigned chrom_name_len)
{
    int32_t word_index = ref_get_contig_word_index_do (chrom_name, chrom_name_len, 0, contig_words_sorted_index.len-1);

    ASSERT (word_index != NIL, "Error: contig \"%.*s\" is observed in %s but is not found in the reference %s",
            chrom_name_len, chrom_name, txt_name, ref_filename);

    return word_index;
}

// returns a pointer to the snip in contig_dict of the chrom_name given
void ref_contig_word_index_from_name (const char *chrom_name, unsigned chrom_name_len, 
                                      const char **snip, int32_t *word_index) // out
{
    *word_index = ref_get_contig_word_index (chrom_name, chrom_name_len);
    MtfWord *word = ENT (MtfWord, contig_words, *word_index);
    *snip = ENT (const char, contig_dict, word->char_index);
}


// in case of REF_EXTERNAL - we use the FASTA_CONTIG word_index from the reference FASTA file to build the range_id
// we allow the first 64 contigs (assumingly the chromosomes) to have up to MAX_POS bp and the rest to have
// up to 16Mbp (assumingly alt contigs, decoys etc)
// Our 2^20 ranges are divided: 
//  - 128 primary contigs * MAX_POS (2^32) / REF_NUM_SITES_PER_RANGE (2^20) = 512K ranges
//  - 32K secondary contigs * 16Mbp (2^24) / REF_NUM_SITES_PER_RANGE (2^20) = 512K ranges
#define REF_NUM_PRIMARY_CONTIGS       (1<<7 ) // 128
#define REF_NUM_SECONDARY_CONTIGS     (1<<15) // 32K
#define REF_MAX_POS_SECONDARY_CONTIG  (1<<24) // 16 Mbp
/*static inline uint32_t ref_range_id_by_word_index (VBlock *vb, uint32_t range_i)
{
    // note: before starting ZIP, we copied the pre-loaded external reference FASTA CONTIG dictionary into our chrom
    // dictionary, ensuring the that chrom_node_index of the chroms in the file, if they exist in reference, are the same
    ASSERT (vb->chrom_node_index < contig_words.len, "Error in ref_range_id_by_word_index: chrom \"%.*s\" in %s does not appear in the reference %s",
            vb->chrom_name_len, vb->chrom_name, txt_name, ref_filename);

    if (vb->chrom_node_index < REF_NUM_PRIMARY_CONTIGS) {
        ASSERT (range_i < pos2range_i (MAX_POS+1), "Error in ref_range_id_by_word_index: range_i=%u is out of range", range_i);

        // 20 bits id: MSb=1 ; 7 bits are the contig, 12 MSb are the range_i 
        return (1 << 19) | (vb->chrom_node_index << 12) | range_i;
    }

    else {
        ASSERT (range_i < REF_MAX_POS_SECONDARY_CONTIG / REF_NUM_SITES_PER_RANGE, "Error: this FASTA cannot be used as a reference: contig %.*s, which is the #%u contig in the FASTA file, has %u bp or more, but contigs beyond the first %u contigs are only allowed up to %u bps.", 
                vb->chrom_name_len, vb->chrom_name,
                vb->chrom_node_index+1, (uint32_t)range_i2pos (range_i), REF_NUM_PRIMARY_CONTIGS, (uint32_t)range_i2pos (range_i)-1);
    
        // 20 bits: MSb=0 ; 15 bit contig ; 4 bit range_i
        return (0 << 19) | ((vb->chrom_node_index - REF_NUM_PRIMARY_CONTIGS) << 4) | range_i;
    }
}
*/
static int32_t ref_get_index_of_chrom_with_alt_name (VBlockP vb)
{
    if (vb->chrom_name_len > 2 || !str_is_int (vb->chrom_name, vb->chrom_name_len))  // we only handle 1 or 2 digit chrom names
        goto fail;

    char chr_chrom[5] = "chr";
    chr_chrom[3] = vb->chrom_name[0];
    chr_chrom[4] = (vb->chrom_name_len == 2 ? vb->chrom_name[0] : 0);

    int32_t alternative_chrom_word_index = ref_get_contig_word_index (chr_chrom, vb->chrom_name_len+3); 
    if (alternative_chrom_word_index == NIL) goto fail;

    return (uint32_t)alternative_chrom_word_index;

fail:
    return vb->chrom_node_index;
}

// ZIP: Allocated and initializes the ref and mutex buffers for the given chrom/range
// case 1: ZIP, reading an external reference: called ahead of committing a fasta vb contig to the reference (REF_EXTERNAL/REF_EXT_STORE)
// case 2: ZIP: in SAM with REF_INTERNAL, when segging a SEQ field ahead of committing it to the reference
// case 3: ZIP: SAM and VCF with REF_EXTERNAL: when segging a SAM_SEQ or VCF_REFALT field, called with RGR_MUST_EXIST
// if range is found, returns a locked range, and its the responsibility of the caller to unlock it. otherwise, returns NULL
Range *ref_zip_get_locked_range (VBlockP vb, int64_t pos)  // out - index of pos within range
{
    uint32_t range_i = flag_reference == REF_INTERNAL ? pos2range_i (pos) : 0; // range within chromosome (only non-0 for REF_INTERNAL)

    // case: we're asking for the same range as the previous one (for example, subsequent line in a sorted SAM)
    if (vb && vb->prev_range && vb->prev_range_chrom_node_index == vb->chrom_node_index && vb->prev_range_range_i == range_i) {
        mutex_lock (vb->prev_range->mutex);
        return vb->prev_range;
    }

    // sanity checks
    ASSERT0 (vb->chrom_name, "Error in ref_zip_get_locked_range: vb->chrom_name=NULL");

    uint32_t save_chrom_index = vb->chrom_node_index;

    if (!flag_reading_reference && (flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE)) { // segging VCF or SAM with external reference

        // if the chrom is not in the reference and it is numeric only, attempt to change it eg "22"->"chr22"
        if (vb->chrom_node_index >= contig_words.len) 
            vb->chrom_node_index = ref_get_index_of_chrom_with_alt_name (vb); // change temporarily just for ref_range_id_by_word_index()
        
        ASSERT (vb->chrom_node_index < contig_words.len, "Error in ref_zip_get_locked_range: chrom \"%.*s\" is not found in the reference file %s",
                MIN (vb->chrom_name_len, 100), vb->chrom_name, ref_filename);
    }

    // one_chrom_ranges_buf is an array of indeces into ranges - one index per range of size REF_NUM_SITES_PER_RANGE
    // note: in case of REF_INTERNAL, we can get hash conflicts, but in case of REF_EXTERNAL, it is guaranteed that there are no conflicts
    uint32_t range_id = (flag_reference == REF_INTERNAL) ? ref_range_id_by_hash (vb, range_i) : vb->chrom_node_index;

    ASSERT (range_id < ranges.len, "Error in ref_zip_get_locked_range: range_id=%u expected to be smaller than ranges.len=%u", range_id, (uint32_t)ranges.len);
    
    Range *range = ENT (Range, ranges, range_id);
    mutex_lock (range->mutex);

    if (range->ref.num_of_bits) {

        // when using an external refernce, pos has to be within the reference range
        ASSERT (flag_reference == REF_INTERNAL || (pos >= range->first_pos && pos <= range->last_pos), 
                "Error: file has POS=%"PRId64" for contig \"%.*s\", but this contig's range is %"PRId64" - %"PRId64,
                pos, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);

        // check for hash conflict (can only happen in REF_INTERNAL) 
        if (flag_reference == REF_INTERNAL && 
            (range->range_i != range_i || vb->chrom_name_len != range->chrom_name_len || memcmp (vb->chrom_name, range->chrom_name, vb->chrom_name_len))) {
            mutex_unlock (range->mutex);
            range = NULL;  // no soup for you
        }

        goto done; // range already initialized or cannot be initialized ^ - we're done
    }

    // check that the range exists, we we're expecting it to (when zipping against an external reference)
    ASSERT (flag_reference == REF_INTERNAL, "Error in ref_zip_get_locked_range: range_id=%u in reference (pos=%"PRId64" chrom=%u \"%.*s\") does not contain any data",
            range_id, pos, vb->chrom_node_index, vb->chrom_name_len, vb->chrom_name);

    // the following applies only to REF_INTERNAL
    uint64_t ref_size = REF_NUM_SITES_PER_RANGE;

    range->range_i          = range_i;
    range->first_pos        = range_i2pos (range_i);
    range->last_pos         = range->first_pos + ref_size - 1;
    range->chrom_name_len   = vb->chrom_name_len;
    
    // the index included in the reference - either the original vb->chrom_node_index or the alternative. 
    // Note: in REF_INTERNAL the chrom index is private to the VB prior to merge, so we can't use it
    range->chrom            = (flag_reference != REF_INTERNAL) ? vb->chrom_node_index : NIL; 
    
    // in case 1, hrom_name points into z_file->contexts of fasta that we be closed when we're done reading.
    // in case 2, chrom_name points into vb->txt_data that will disappear at the end of this VB, 
    // case 3 doesn't reach here. therefore, we make a copy of chrom_name.
    range->chrom_name = malloc (vb->chrom_name_len);
    memcpy ((char *)range->chrom_name, vb->chrom_name, vb->chrom_name_len);

    // nothing is set yet - in ZIP, bits get set as they are encountered in the compressing txt file
    bit_array_alloc (&range->ref, ref_size * 2);
    bit_array_alloc (&range->is_set, ref_size);
    bit_array_clear_region (&range->is_set, 0, ref_size);

    if (vb) {
        vb->prev_range = range;
        vb->prev_range_chrom_node_index = vb->chrom_node_index;
        vb->prev_range_range_i = range_i;
    }

done:
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

    ASSERT (*sl < AFTERENT (SectionListEntry, ref_file_section_list), "Error in ref_copy_one_compressed_section: cannot find FASTA_SEQ of vb_i=%u in section list of reference file", ra->vblock_i);

    static Buffer ref_seq_section = EMPTY_BUFFER;

    SAVE_FLAG (flag_show_headers);
    zfile_read_section (ref_file, evb, ra->vblock_i, NO_SB_I, &ref_seq_section, "ref_seq_section", 
                        sizeof (SectionHeaderReference), SEC_REFERENCE, *sl);
    RESTORE_FLAG (flag_show_headers);

    SectionHeaderReference *header = (SectionHeaderReference *)ref_seq_section.data;

    // some minor changes to the header...
    header->h.vblock_i  = 0; // we don't belong to any VB and there is no encryption of external ref
    header->h.section_i = 0; // we don't belong to any VB

    // Write header and body of the reference to z_file
    // Note on encryption: reference sections originating from an external reference are never encrypted - not
    // by us here, and not in the source reference fasta (because with disallow --make-reference in combination with --password)
    START_TIMER;
    file_write (z_file, ref_seq_section.data, ref_seq_section.len);
    COPY_TIMER (evb->profile.write);

    // "manually" add the reference section to the section list - normally it is added in comp_compress()
    sections_add_to_list (evb, &header->h);

    z_file->disk_so_far += ref_seq_section.len;   // length of GENOZIP data writen to disk

    if (flag_show_reference) {
        Context *ctx = &z_file->contexts[DTFZ(chrom)];
        MtfNode *node = ENT (MtfNode, ctx->mtf, ra->chrom_index);
        fprintf (stderr, "Copying reference    %.20s %"PRId64" - %"PRId64"\n", ENT (char, ctx->dict, node->char_index), ra->min_pos, ra->max_pos);
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
        int64_t fasta_sec_start_in_contig_r = fasta_sec[i].min_pos - contig_r->first_pos; // the start of the FASTA section (a bit less than 1MB) within the full-contig range
        int64_t fasta_sec_len = fasta_sec[i].max_pos - fasta_sec[i].min_pos + 1;
        int64_t bits_is_set   = bit_array_num_bits_set_region (&contig_r->is_set, fasta_sec_start_in_contig_r, fasta_sec_len);

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
static bool ref_remove_flanking_regions (Range *r, uint64_t r_num_set_bits, uint64_t *start_flanking_region_len)
{
// threshold - if the number of clear bits (excluding flanking regions) is beneath this, we will not copmact, as the cost in
// z_file size of a SEC_REF_IS_SET section needed if compacting will be more than what we save in compression of r->ref
#define THRESHOLD_FOR_COMPACTING 470 

    uint64_t end_flanking_region_len, last_1;
    
    char has_any_bit = bit_array_find_first_set_bit (&r->is_set, start_flanking_region_len);
    ASSERT0 (has_any_bit, "Error in ref_remove_flanking_regions: range has no bits set"); // ref_prepare_range_for_compress is responsible not to send us 0-bit ranges

    has_any_bit = bit_array_find_prev_set_bit (&r->is_set, r->is_set.num_of_bits, &last_1);
    ASSERT0 (has_any_bit, "Error in ref_remove_flanking_regions: range has no bits set #2"); // this should definitely never happen, since we already know the range has bits
    end_flanking_region_len = r->is_set.num_of_bits - last_1 - 1;

    uint64_t num_clear_bits_excluding_flanking_regions = 
        r->is_set.num_of_bits - r_num_set_bits - *start_flanking_region_len - end_flanking_region_len;

    // remove flanking regions - will allow a smaller allocation for the reference in PIZ 
    r->first_pos += *start_flanking_region_len;
    r->last_pos  -= end_flanking_region_len;

    ASSERT0 (r->last_pos >= r->first_pos, "Error in ref_compact_ref: bad removal of flanking regions");

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

    // in case of REF_INTERNAL, chrom_word_index might stil be NIL. We get it now from the z_file data
    if (r && r->chrom == NIL)
        r->chrom = ref_get_contig_word_index (r->chrom_name, r->chrom_name_len); // this changes first/last_pos, removing flanking regions

    SectionHeaderReference header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.vblock_i          = BGEN32 (vb->vblock_i); 
    header.h.section_i         = 0; 
    header.h.magic             = BGEN32 (GENOZIP_MAGIC);
    header.h.compressed_offset = BGEN32 (sizeof(header));
    header.chrom_word_index    = r ? BGEN32 (r->chrom) : NIL;
    header.first_pos           = r ? BGEN64 ((uint64_t)r->first_pos) : 0;

    vb->z_data.name  = "z_data"; // comp_compress requires that these are pre-set
    vb->z_data.param = vb->vblock_i;

    // First, SEC_REF_IS_SET section (but not needed if the entire range is used)
    int64_t compacted_ref_size = r ? (r->ref.num_of_bits / 2) : 0; 

    if (r && is_compacted) {

        LTEN_bit_array (&r->is_set, true);

        header.h.section_type          = SEC_REF_IS_SET;  // most of the header is the same as ^
        header.h.sec_compression_alg   = COMP_BZ2;
        header.h.data_uncompressed_len = BGEN32 (r->is_set.num_of_words * sizeof (uint64_t));
        header.last_pos                = BGEN64 ((uint64_t)r->last_pos); // full length, after flanking regions removed
        comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, (char *)r->is_set.words, NULL);

        if (flag_show_reference && r) 
            fprintf (stderr, "vb_i=%u Compressing SEC_REF_IS_SET chrom=%u (%.*s) pos=[%"PRId64",%"PRId64"] (%"PRId64" bases) section_size=%u bytes\n", 
                    vb->vblock_i, BGEN32 (header.chrom_word_index), r->chrom_name_len, r->chrom_name, BGEN64 (header.first_pos), BGEN64 (header.last_pos), BGEN64 (header.last_pos) - BGEN64 (header.first_pos) + 1, 
                    BGEN32 (header.h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));
    }

    // Second. SEC_REFERENCE
    if (r) LTEN_bit_array (&r->ref, true);

    header.h.section_type          = SEC_REFERENCE;
    header.h.sec_compression_alg   = COMP_LZMA;
    header.h.data_uncompressed_len = r ? BGEN32 (r->ref.num_of_words * sizeof (uint64_t)) : 0;
    // we shorten the length if compacted - keeping first_pos intact, and last_pos reflecting the number of bits
    header.last_pos                = r ? BGEN64 ((uint64_t)(r->first_pos + compacted_ref_size - 1))  : 0;  // this is actually the same as r->last_pos if not compacted
    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, r ? (char *)r->ref.words : NULL, NULL);

    if (flag_show_reference && r) 
        fprintf (stderr, "vb_i=%u Compressing SEC_REFERENCE chrom=%u (%.*s) %s=[%"PRId64",%"PRId64"] (%"PRId64" bases) section_size=%u bytes\n", 
                 vb->vblock_i, BGEN32 (header.chrom_word_index), r->chrom_name_len, r->chrom_name, is_compacted ? "compacted" : "pos",
                 BGEN64 (header.first_pos), BGEN64 (header.last_pos), BGEN64 (header.last_pos) - BGEN64 (header.first_pos) + 1, 
                 BGEN32 (header.h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderReference));

    // store the ref_stored_ra data for this range
    if (r) {
        mutex_lock (ref_stored_ra_mutex);
        NEXTENT (RAEntry, ref_stored_ra) = (RAEntry){ .vblock_i    = vb->vblock_i, 
                                                      .chrom_index = r->chrom,
                                                      .min_pos     = r->first_pos,
                                                      .max_pos     = r->last_pos   };
        mutex_unlock (ref_stored_ra_mutex);
    }

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

static void ref_output_one_range (VBlockP vb) 
{ 
    zip_output_processed_vb (vb, &vb->section_list_buf, false, PD_REFERENCE_DATA); 
}

static inline bool ref_is_range_used (const Range *r)
{
    // note: in make_ref mode, is_set is always 0, and ref is 0 for unused ranges
    return r->ref.num_of_bits && (r->is_set.num_of_bits || flag_make_reference);
}

// compress the reference - one section at the time, using Dispatcher to do them in parallel 
static void ref_prepare_range_for_compress (VBlockP vb)
{
    static uint32_t next_range_i=0;
    if (vb->vblock_i == 1) next_range_i=0; // initialize for once for non-bound file

    // find next occupied range
    for (; !vb->ready_to_dispatch && next_range_i < ranges.len ; next_range_i++) {

        Range *r = ENT (Range, ranges, next_range_i);

        if (!ref_is_range_used (r)) continue;
//        if ((!flag_make_reference && !r->is_set.num_of_bits) || !r->ref.num_of_bits) 
//            continue; // unused range (note: in make_ref mode, is_set is always 0, and ref is 0 for unused ranges)

        uint64_t num_set = flag_make_reference ? (r->ref.num_of_bits / 2) : bit_array_num_bits_set (&r->is_set);
        if (!num_set) {
            r->is_set.num_of_bits = 0;
            continue; // nothing to with this range - perhaps copied and cleared in ref_copy_compressed_sections_from_reference_file
        }

        // in make_reference, we have exactly one contig for each VB, one one RAEntry for that contig
        // during seg we didn't know the chrom,first,last_pos, so we add them now, from the RA
        if (flag_make_reference) 
            random_access_get_ra_info (vb->vblock_i, &r->chrom, &r->first_pos, &r->last_pos);

        vb->range = r; // range to compress
        vb->range_num_set_bits = num_set;
        vb->ready_to_dispatch = true;
    }
}

// ZIP: compress and write reference sections. either compressed by us, or copied directly from the reference file.
void ref_compress_ref (void)
{
    if (!buf_is_allocated (&ranges)) return;

    // copy already-compressed SEQ sections from the FASTA genozip reference, but only such sections that are entirely
    // covered by ranges with is_accessed=true. we mark these ranges affected as is_accessed=false.
    if (flag_reference == REF_EXT_STORE)
        ref_copy_compressed_sections_from_reference_file ();

    // make a copy of the file's own chrom dictionary into the reference chrom dictionary, for code consistency with
    // the case where it is read from the external reference
    else if (flag_reference == REF_INTERNAL)
        ref_copy_chrom_data_from_z_file();

    // initialize ref_stored_ra for the creation of SEC_REF_RANDOM_ACC section
    uint32_t count_used_ranges=0;
    for (int32_t i=0; i < ranges.len; i++)
        if (ref_is_range_used (ENT (Range, ranges, i)))
            count_used_ranges++;

    buf_alloc (evb, &ref_stored_ra, sizeof (RAEntry) * count_used_ranges, 1, "ref_stored_ra", 0);
    mutex_initialize (&ref_stored_ra_mutex, &ref_stored_ra_mutex_initialized);

    // compression of reference doesn't output % progress
    SAVE_FLAG (flag_quiet);
    if (flag_show_reference) flag_quiet = true; // show references instead of progress

    // proceed to compress all ranges that have still have data in them after copying
    uint32_t num_vbs_dispatched = 
        dispatcher_fan_out_task (NULL, "Writing reference...", false, 
                                 ref_prepare_range_for_compress, 
                                 ref_compress_one_range, 
                                 ref_output_one_range);

    RESTORE_FLAG (flag_quiet);
    
    // SAM require at least one reference section, but if the SAM is unaligned, there will be none - create one empty section
    // (this will also happen if SAM has just only reference section, we will just needless write another tiny section - no harm)
    if (z_file->data_type == DT_SAM && num_vbs_dispatched==1) {
        evb->range = NULL;
        ref_compress_one_range (evb); // incidentally, will also be written in case of a small (one vb) reference - no harm
    }

    if (flag_show_ref_index) 
        random_access_show_index (&ref_stored_ra, true, "Reference random-access index contents (result of --show-ref-index)");

    // compress reference random access
    random_access_finalize_entries (&ref_stored_ra); // sort in order of vb_i

    BGEN_random_access (&ref_stored_ra); // make ra_buf into big endian
    ref_stored_ra.len *= sizeof (RAEntry); // change len to count bytes
    zfile_compress_section_data_alg (evb, SEC_REF_RANDOM_ACC, &ref_stored_ra, 0,0, COMP_LZMA); // ra data compresses better with LZMA than BZLIB
}

// -------------------------------
// Importing an external reference
// -------------------------------

void ref_set_reference (const char *filename)
{
    ref_filename = filename;
}

void ref_set_md5 (Md5Hash md5)
{
    ref_md5 = md5;
}

// ZIP & PIZ: import external reference
void ref_load_external_reference (void)
{
    ASSERT0 (ref_filename, "Error: ref_filename is NULL");

    z_file = file_open (ref_filename, READ, Z_FILE, DT_FASTA);    
    z_file->basename = file_basename (ref_filename, false, "(reference)", NULL, 0);
        
    // save and reset some globals that after pizzing FASTA
    int save_flag_test         = flag_test        ; flag_test        = 0;
    int save_flag_md5          = flag_md5         ; flag_md5         = 0;
    int save_flag_show_time    = flag_show_time   ; flag_show_time   = 0;
    int save_flag_show_memory  = flag_show_memory ; flag_show_memory = 0;
    int save_flag_no_header    = flag_no_header   ; flag_no_header   = 0;
    int save_flag_header_one   = flag_header_one  ; flag_header_one  = 0;
    int save_flag_header_only  = flag_header_only ; flag_header_only = 0;
    char *save_flag_grep       = flag_grep        ; flag_grep        = 0;
    CommandType save_command   = command;         ; command          = PIZ;
    flag_reading_reference = true; // tell fasta.c that this is a reference
    
    bool piz_successful = piz_dispatcher (true, false);

    // recover globals
    flag_reading_reference = false;
    command          = save_command;
    flag_test        = save_flag_test;
    flag_md5         = save_flag_md5;
    flag_show_time   = save_flag_show_time;
    flag_show_memory = save_flag_show_memory;
    flag_no_header   = save_flag_no_header;
    flag_header_only = save_flag_header_only;
    flag_header_one  = save_flag_header_one;
    flag_grep        = save_flag_grep;

    ASSERT (piz_successful, "Error: failed to uncompress reference file %s", ref_filename);
    
    file_close (&z_file, false);
    file_close (&txt_file, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtfile_genozip_to_txt_header.
}

// case 1: in case of ZIP with external reference, called by ref_load_stored_reference during piz_read_global_area of the reference file
// case 2: in case of PIZ: also called from ref_load_stored_reference one_range_per_contig=true
// case 3: in case of ZIP of SAM using internal reference - called from sam_zip_initialize
// note: ranges allocation must be called by the I/O thread as it adds a buffer to evb buf_list
void ref_initialize_ranges (bool one_range_per_contig)
{
    ranges.len = one_range_per_contig ? (flag_reference == REF_STORED ? z_file->contexts[DTFZ(chrom)].word_list.len : contig_words.len)
                                      : REF_NUM_RANGES;
    buf_alloc (evb, &ranges, ranges.len * sizeof (Range), 1, "ranges", 0); 
    buf_zero (&ranges);

    for (unsigned i=0; i < ranges.len; i++) {
        Range *r = ENT (Range, ranges, i);

        pthread_mutex_init (&r->mutex, NULL);

        // in PIZ or when reading an external range, we allocate and initialize one range per chromosome (in ZIP from REF_INTERNAL, we will allocate ranges as they are encoutered in the data)
        if (one_range_per_contig) { // PIZ always and ZIP with external reference

            r->chrom = i;      

            random_access_pos_of_chrom (i, &r->first_pos, &r->last_pos); // in reference fasta.genozip if REF_EXTERNAL, or in file itself is REF_STORED
            
            if (flag_reference == REF_STORED) {
                Context *ctx = &z_file->contexts[DTFZ(chrom)];
                mtf_get_snip_by_word_index (&ctx->word_list, &ctx->dict, r->chrom, &r->chrom_name, &r->chrom_name_len);
            }
            else
                mtf_get_snip_by_word_index (&contig_words, &contig_dict, r->chrom, &r->chrom_name, &r->chrom_name_len);

            // case: this chrom has no RA - for example, it is a "=" or "*". No need to allocate memory for it (we do get the chrom name ^ to ease debugging)
            if (r->last_pos < 0) continue;

            bit_array_alloc (&r->ref,    ref_size(r) * 2); 
            bit_array_alloc (&r->is_set, ref_size(r));

            // in PIZ, we start by clearing the entire contig, and then in ref_uncompress_one_range we set those regions
            // that are read from SEC_REFERENCE sections in a stored or external reference
            if (primary_command == PIZ) 
                bit_array_clear_all (&r->is_set); // entire contig is cleared
        }
    }
}

static int ref_sort_words_alphabetically (const void *a, const void *b)
{
    uint32_t index_a = *(uint32_t *)a;
    uint32_t index_b = *(uint32_t *)b;

    MtfWord *word_a = ENT (MtfWord, contig_words, index_a);
    MtfWord *word_b = ENT (MtfWord, contig_words, index_b);
    
    return strcmp (ENT (char, contig_dict, word_a->char_index),
                   ENT (char, contig_dict, word_b->char_index));
}

static void ref_copy_chrom_data_from_z_file (void)
{
    // copy data from the reference FASTA's CONTIG context, so it survives after we finish reading the reference and close z_file
    Context *chrom_ctx = &z_file->contexts[DTFZ(chrom)];

    ASSERT (flag_reference == REF_INTERNAL || (buf_is_allocated (&chrom_ctx->dict) && buf_is_allocated (&chrom_ctx->word_list)),
            "Error: cannot use %s as a reference as it is missing a CONTIG dictionary", z_name);

    buf_copy (evb, &contig_dict, &chrom_ctx->dict, 1, 0, 0, "contig_dict", 0);
    
    // case: in REF_INTERNAL, we copy from the z_file context after we completed segging a file
    // note: we can't rely on chrom_ctx->mtf for the correct order as vb_i=1 resorted. rather, by mimicking the
    // word_list generation as done in PIZ, we guarantee that we will get the same chrom_index
    // in case of multiple bound files, we re-do this in every file in case of additional chroms (not super effecient, but good enough because the context is small)
    if (flag_reference == REF_INTERNAL) {
        contig_words.len = chrom_ctx->mtf.len;
        buf_alloc (evb, &contig_words, contig_words.len * sizeof (MtfWord), 1, "contig_words", 0);

        // similar logic to mtf_integrate_dictionary_fragment
        char *start = contig_dict.data;
        for (uint32_t snip_i=0; snip_i < contig_words.len; snip_i++) {

            MtfWord *word = ENT (MtfWord, contig_words, snip_i);

            char *c=start; while (*c != SNIP_SEP) c++;
            word->snip_len   = c - start;
            word->char_index = start - contig_dict.data;

            start = c+1; // skip over the SNIP_SEP
        }
    }
    else
        buf_copy (evb, &contig_words, &chrom_ctx->word_list, sizeof (MtfWord), 0, 0, "contig_words", 0);

    // contig_words_sorted_index - an array of uint32 of indexes into contig_words - sorted by alphabetical order of the snip in contig_dict
    buf_alloc (evb, &contig_words_sorted_index, sizeof(uint32_t) * contig_words.len, 1, "contig_words_sorted_index", 0);
    for (uint32_t i=0; i < contig_words.len; i++)
        NEXTENT (uint32_t, contig_words_sorted_index) = i;

    qsort (contig_words_sorted_index.data, contig_words_sorted_index.len, sizeof(uint32_t), ref_sort_words_alphabetically);
}

// ZIP & PIZ: called by piz_read_global_area during PIZ of the external reference fa.genozip file
void ref_consume_ref_fasta_global_area (void)
{
    ref_copy_chrom_data_from_z_file();

    random_access_pos_of_chrom (0, 0, 0); // initialize if not already initialized
    buf_copy (evb, &contig_maxes, &z_file->ra_min_max_by_chrom, sizeof (int64_t), 0, 0, "contig_maxes", 0);
    buf_copy (evb, &ref_external_ra, &z_file->ra_buf, sizeof (RAEntry), 0, 0, "ref_external_ra", 0);
    buf_copy (evb, &ref_file_section_list, &z_file->section_list_buf, sizeof (SectionListEntry), 0, 0, "ref_file_section_list", 0);
}

// ZIP & PIZ: used when zipping or pizzing the target file, after the external reference is read
int64_t ref_min_max_of_chrom (int32_t chrom, bool get_max)
{
    ASSERT (chrom >= 0 && chrom < contig_maxes.len, "Error in ref_min_max_of_chrom: chrom=%d out of range, contig_maxes.len=%u",
            chrom, (uint32_t)contig_maxes.len);

    ASSERT0 (flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE, "Error in ref_min_max_of_chrom: function can only be called with external reference");

    typedef struct { int64_t min_pos, max_pos; } MinMax;
    MinMax *mm = ENT (MinMax, contig_maxes, chrom);
    return get_max ? mm->max_pos : mm->min_pos;
}

//---------------------------------------
// Reference creation stuff
//---------------------------------------

static pthread_mutex_t make_ref_mutex;

#define MAKE_REF_NUM_RANGES 1000000 // should be more than enough (in GRCh38 we have 6389)

// in make_ref, each VB gets it own range indexed by vb->vblock_i - so they can work on them in parallel without
// worrying about byte overlap.
void ref_make_ref_init (void)
{
    ASSERT0 (flag_make_reference, "Expecting flag_make_reference=true");

    buf_alloc (evb, &ranges, MAKE_REF_NUM_RANGES * sizeof (Range), 1, "ranges", 0); // must be allocated by I/O thread as its evb
    buf_zero (&ranges);

    pthread_mutex_init (&make_ref_mutex, NULL);
}

Range *ref_make_ref_get_range (uint32_t vblock_i)
{
    mutex_lock (make_ref_mutex);
    ranges.len = MAX (ranges.len, (uint64_t)vblock_i + 1); // note that this function might be called out order (called from FASTA ZIP compute thread)
    mutex_unlock (make_ref_mutex);

    ASSERT (ranges.len <= MAKE_REF_NUM_RANGES, "Error in ref_make_ref_get_range: reference file too big - number of ranges exceeds %u", MAKE_REF_NUM_RANGES);

    return ENT (Range, ranges, vblock_i);
}

//---------------------------------------
// Printing
//---------------------------------------

void ref_print_subrange (const char *msg, const Range *r, int64_t start_pos, int64_t end_pos) /* start_pos=end_pos=0 if entire ref */
{
    uint64_t start_idx = start_pos ? start_pos - r->first_pos : 0;
    uint64_t end_idx   = (end_pos ? MIN (end_pos, r->last_pos) : r->last_pos) - r->first_pos;

    fprintf (stderr, "%s: %.*s %"PRId64" - %"PRId64" (len=%u): ", msg, r->chrom_name_len, r->chrom_name, start_pos, end_pos, (uint32_t)(end_pos - start_pos + 1));
    for (uint64_t idx = start_idx; idx <= end_idx; idx++) 
        fputc (ref_get_nucleotide (r, idx) + (32 * !ref_is_nucleotide_set (r, idx)), stderr); // uppercase if set, lowercase if not

    fputc ('\n', stderr);
}

void ref_print_is_set (const Range *r)
{
    fprintf (stderr, "\n\nRegions set for chrom %u \"%.*s\" [%"PRId64"-%"PRId64"] according to range.is_set (format- \"first_pos-last_pos (len)\")\n", 
             r->chrom, r->chrom_name_len, r->chrom_name, r->first_pos, r->last_pos);

    if (!r->is_set.num_of_bits)
        fprintf (stderr, "No data: r->is_set.num_of_bits=0\n");

    uint64_t next;
    for (uint64_t i=0; i < r->is_set.num_of_bits; ) {
        
        bool found = bit_array_find_next_clear_bit (&r->is_set, i, &next);
        if (!found) next = r->is_set.num_of_bits;
        if (next > i) {
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