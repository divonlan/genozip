
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

static Buffer ranges = EMPTY_BUFFER;      // ZIP: a buffer containing a array of REF_NUM_RANGES Range's - each allocated on demand. PIZ - containing one element per chrom.
static CommandType ref_primary_command;   // the command that caused us to read the reference - either ZIP or UNZIP

// data derived from the external reference FASTA's CONTIG context
static Buffer contig_dict               = EMPTY_BUFFER;
static Buffer contig_words              = EMPTY_BUFFER;
static Buffer contig_maxes              = EMPTY_BUFFER; // an array of int64 - the max_pos of each contig
static Buffer contig_words_sorted_index = EMPTY_BUFFER; // an array of uint32 of indexes into contig_words - sorted by alphabetical order of the snip in contig_dict
static Buffer ref_file_ra               = EMPTY_BUFFER;
static Buffer ref_file_section_list     = EMPTY_BUFFER;

// globals
bool ref_flag_reading_reference = false;
const char *ref_filename = NULL; // filename of external reference file
Md5Hash ref_md5 = MD5HASH_NONE;

#define SAVE_FLAG(flag) int save_##flag = flag ; flag=0
#define RESTORE_FLAG(flag) flag = save_##flag

// called by mtf_copy_reference_contig_to_chrom_ctx when initializing ZIP for a new file using a pre-loaded external reference
void ref_get_contigs (ConstBufferP *out_contig_dict, ConstBufferP *out_contig_words)
{
    *out_contig_dict  = &contig_dict;
    *out_contig_words = &contig_words;
}

// free memory allocations between files, when compressing multiple non-concatenated files or decompressing multiple files
void ref_cleanup_memory(void)
{
    if (flag_reference != REF_INTERNAL && flag_reference != REF_STORED) return; // we never cleanup external references as they can never change (the user can only specify one --reference)

    if (buf_is_allocated (&ranges)) {
        ARRAY (Range, rng, ranges);
        for (unsigned i=0; i < ranges.len ; i++) {
            if (rng[i].ref) FREE (rng[i].ref);
            if (flag_reference == REF_INTERNAL && rng[i].chrom_name) FREE ((char *)rng[i].chrom_name); // allocated only in REF_INTERNAL
        }

        buf_free (&ranges);
    }
}

MemStats ref_memory_consumption (void)
{
    MemStats stats = { "reference", 0, 0 };

    ARRAY (Range, r, ranges);
    if (buf_is_allocated (&ranges)) {
        for (unsigned i=0; i < ranges.len; i++)
            if (r[i].ref_size) {
                stats.bytes += r[i].ref_size;
                stats.buffers++;
            }
    }

    return stats;
}

static uint32_t ref_get_alternative_chrom_name (VBlockP vb); // forward

const char *ref_get_ref (VBlockP vb, uint64_t pos, uint32_t ref_consumed)
{
    // if the chrom is not in the reference and it is numeric only, attempt to change it eg "22"->"chr22"
    uint32_t index = vb->chrom_node_index < ranges.len ? vb->chrom_node_index : ref_get_alternative_chrom_name (vb);

    ASSERT (index < ranges.len, "Error in ref_get_ref: vb->chrom_node_index=%u is out of range: ranges.len=%u",
            index, (uint32_t)ranges.len);

    Range *r = ENT (Range, ranges, index);
    if (!r->ref) return NULL; // this can if entire chromosome is verbatim, eg. unaligned (pos=4) or SEQ or CIGAR are unavailable

    ASSERT (pos + ref_consumed - 1 <= r->last_pos, "Error in ref_get_ref: out of range reconstructing txt_line_i=%u: pos=%"PRId64" ref_consumed=%u but range->last_pos=%"PRId64,
            vb->line_i, pos, ref_consumed, r->last_pos);

    return &r->ref[pos];
}

// ----------------------------------------------------------
// PIZ: read and uncompress stored ranges (if no --reference)
// ----------------------------------------------------------

static void ref_uncompress_one_stored_range (VBlockP vb)
{
    SectionHeaderReference *header = (SectionHeaderReference *)vb->z_data.data;

    uint32_t ref_i      = BGEN32 (header->chrom_word_index);
    int64_t first_pos   = BGEN64 (header->first_pos);
    int64_t last_pos    = BGEN64 (header->last_pos);
    uint32_t uncomp_len = BGEN32 (header->h.data_uncompressed_len);

    Context *ctx = &z_file->contexts[DTFZ(chrom)];
    const MtfWord *word = ENT (const MtfWord, ctx->word_list, ref_i);

    Range *r = ENT (Range, ranges, ref_i);
    ASSERT (last_pos <= r->last_pos, "Error in ref_uncompress_one_stored_range: ref range out of bounds for ref_i=%u chrom=%.*s: first_pos=%"PRId64" last_pos=%"PRId64" but range->last_pos=%"PRId64,
            ref_i, word->snip_len, ENT (char, ctx->dict, word->char_index), first_pos, last_pos, r->last_pos);

    char *uncompressed_data = &r->ref[first_pos];

    zfile_uncompress_section (vb, (SectionHeaderP)header, uncompressed_data, NULL, SEC_REFERENCE);

    if (!header->h.data_uncompressed_len) return;  // reference contains a single empty header, in case of an unaligned SAM file

    // case: ref has been compacted - we uncompact it
    // note: the unused gap areas of reference contain undefined data and are not used. we are careful not to zero
    // them as this way the OS lazy-allocation algorithm doesn't actually consume memory blocks that are entirely unused
    if (uncomp_len < last_pos - first_pos + 1) {

        // traverse backwards and "stretch" the reference 
        char *next_ref = &r->ref[last_pos]; // the last char of the uncompacted ref
        for (char *next_com = &uncompressed_data[uncomp_len-1]; next_com >= uncompressed_data; next_com--) { 

            if (*next_com == '\t') {
                uint32_t gap_len = (uint32_t)(uint8_t)next_com[-3] | ((uint32_t)(uint8_t)next_com[-2] << 8) | ((uint32_t)(uint8_t)next_com[-1] << 16);
                next_ref -= gap_len;
                next_com -= 3; // skip gap length
            }
            else {
                ASSERT (next_ref >= r->ref, "Error in ref_uncompress_one_stored_range: next_ref out of range when uncompacting chrom=%.*s, first_pos=%"PRId64" last_pos=%"PRId64,
                        word->snip_len, ENT (char, ctx->dict, word->char_index), first_pos, last_pos);

                *(next_ref--) = *next_com;
            }
        } 
    }

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. 
}

static SectionListEntry *sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
static uint32_t ref_range_cursor = 0;

static uint64_t ref_calc_chrom_ref_size (uint64_t max_pos)
{
    return (((max_pos+1) + (REF_NUM_SITES_PER_RANGE-1)) / REF_NUM_SITES_PER_RANGE) * REF_NUM_SITES_PER_RANGE;
}

static void ref_read_one_stored_range (VBlockP vb)
{
    if (sections_get_next_section_of_type (&sl_ent, &ref_range_cursor, SEC_REFERENCE)) {
        
        zfile_read_section (z_file, evb, sl_ent->vblock_i, NO_SB_I, &vb->z_data, "z_data", 
                            sizeof(SectionHeaderReference), sl_ent->section_type, sl_ent);    
        
        // allocate memory for entire chrom reference if this is the first range of this chrom
        uint32_t chrom_word_index = BGEN32 (((SectionHeaderReference *)vb->z_data.data)->chrom_word_index);
        Range *r = ENT (Range, ranges, chrom_word_index);

        if (!r->ref) {
            r->last_pos = ref_calc_chrom_ref_size (random_access_max_pos_of_chrom (chrom_word_index)) - 1;
            // note: the OS lazy-allocation will only actually allocate memory blocks that are actually used 
            r->ref_size = r->last_pos + 1;
            r->ref = malloc (r->ref_size);
            ASSERT (r->ref, "Error in ref_read_one_stored_range: failed to allocated %"PRId64" bytes", r->ref_size); 
        }

        vb->ready_to_dispatch = true;
    }
}

void ref_uncompress_all_stored_ranges (void)
{
    ASSERT0 (!buf_is_allocated (&ranges), "Error in ref_uncompress_all_stored_ranges: expecting ranges to be unallocated");
    
    ranges.len = z_file->contexts[DTFZ(chrom)].word_list.len;
    buf_alloc (evb, &ranges, ranges.len * sizeof(Range), 1, "ranges", 0);
    buf_zero (&ranges);
    
    sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
    ref_range_cursor = 0;

    // decompress reference using Dispatcher
    dispatcher_fan_out_task ("Reading reference from genozip file", flag_test, 
                             ref_read_one_stored_range, 
                             ref_uncompress_one_stored_range, 
                             NULL);

    buf_test_overflows_all_vbs ("ref_uncompress_all_stored_ranges");
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
void ref_get_contig (const char *chrom_name, unsigned chrom_name_len, 
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
static inline uint32_t ref_range_id_by_word_index (VBlock *vb, uint32_t range_i)
{
    // note: before starting ZIP, we copied the pre-loaded external reference FASTA CONTIG dictionary into our chrom
    // dictionary, ensuring the that chrom_node_index of the chroms in the file, if they exist in reference, are the same
    ASSERT (vb->chrom_node_index < contig_words.len, "Error in ref_range_id_by_word_index: chrom \"%.*s\" in %s does not appear in the reference %s",
            vb->chrom_name_len, vb->chrom_name, txt_name, ref_filename);

    if (vb->chrom_node_index < REF_NUM_PRIMARY_CONTIGS) {
        ASSERT (range_i < (MAX_POS+1) / REF_NUM_SITES_PER_RANGE, "Error in ref_range_id_by_word_index: range_i=%u is out of range", range_i);

        // 20 bits id: MSb=1 ; 7 bits are the contig, 12 MSb are the range_i 
        return (1 << 19) | (vb->chrom_node_index << 12) | range_i;
    }

    else {
        ASSERT (range_i < REF_MAX_POS_SECONDARY_CONTIG / REF_NUM_SITES_PER_RANGE, "Error: this FASTA cannot be used as a reference: contig %.*s, which is the #%u contig in the FASTA file, has %u bp or more, but contigs beyond the first %u contigs are only allowed up to %u bps.", 
                vb->chrom_name_len, vb->chrom_name,
                vb->chrom_node_index+1, range_i * REF_NUM_SITES_PER_RANGE, REF_NUM_PRIMARY_CONTIGS, range_i * REF_NUM_SITES_PER_RANGE-1);
    
        // 20 bits: MSb=0 ; 15 bit contig ; 4 bit range_i
        return (0 << 19) | ((vb->chrom_node_index - REF_NUM_PRIMARY_CONTIGS) << 4) | range_i;
    }
}

// PIZ: set a subset of the reference of a single chromosome, from external data
static void ref_set_ref_from_external_data_piz (VBlock *vb, uint64_t start_pos, ConstBufferP data)
{
    Range *r = ENT (Range, ranges, vb->chrom_node_index);

    ASSERT (start_pos + data->len - 1 <= r->ref_size, "Error in ref_set_ref_from_external_data_piz: range exceeds allocated ref: chrom=%.*s, start_pos=%"PRId64", data->len=%"PRId64" range->ref_size=%"PRId64,
            vb->chrom_name_len, vb->chrom_name, start_pos, data->len, r->ref_size);

    memcpy (&r->ref[start_pos], data->data, data->len);   
}

// ZIP: set a subset of the reference of a single chromosome, from external data
static void ref_set_ref_from_external_data_zip (VBlock *vb, int64_t start_pos,
                                                ConstBufferP data, uint64_t data_start, bool might_span_multiple_vbs)
{
    ASSERT (data_start <= data->len, "Error in ref_set_ref_from_external_data: expected data_start=%"PRId64" <= data->len=%"PRId64, 
            data_start, data->len);

    if (!data->len) return; // nothing to do

    uint64_t range_i = start_pos / REF_NUM_SITES_PER_RANGE;
    uint64_t this_range_start = start_pos % REF_NUM_SITES_PER_RANGE;
    uint64_t this_range_data_len = MIN ((range_i+1) * REF_NUM_SITES_PER_RANGE - start_pos, data->len - data_start);

    // get range - note possibly two consecutive VBs threads might be writing to the same range at the same time (but different pos range)
    Range *range = ref_get_range (vb, range_i, RGR_CAN_BE_NEW_OR_EXISTING,
                                  might_span_multiple_vbs ? 0 : this_range_data_len); // full size for a range that might be multi-vb

    // copy data (no need for mutex protection as threads write to non-overlapping regions)
    memcpy (&range->ref[this_range_start], &data->data[data_start], this_range_data_len);

    // set range info (note: in case of two consecutive VB threads writing to the same range, we don't know who writes first)
    mutex_lock (range->mutex);
    range->first_pos = MIN (range->first_pos, start_pos);
    range->last_pos  = MAX (range->last_pos,  start_pos + this_range_data_len - 1);
    range->num_set  += this_range_data_len;
    mutex_unlock (range->mutex);

    if (data_start + this_range_data_len < data->len) 
        // continue recursively - tail recursion
        ref_set_ref_from_external_data (vb, start_pos + this_range_data_len, 
                                        data, data_start + this_range_data_len, might_span_multiple_vbs);    
}

void ref_set_ref_from_external_data (VBlock *vb, uint64_t start_pos,
                                     ConstBufferP data, uint64_t data_start, bool might_span_multiple_vbs)
{
    if (ref_primary_command == ZIP)
        ref_set_ref_from_external_data_zip (vb, start_pos, data, data_start, might_span_multiple_vbs);
    else
        ref_set_ref_from_external_data_piz (vb, start_pos, data);
}

static uint32_t ref_get_alternative_chrom_name (VBlockP vb)
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
// case 1: ZIP: called when reading an external sequence ahead of committing it to the reference (REF_EXTERNAL/REF_EXT_STORE)
// case 2: ZIP: in SAM, when reading an SEQ field ahead of committing it to the refernece (REF_INTERNAL)
Range *ref_get_range (VBlockP vb, uint32_t range_i, RgrMode mode, uint32_t reduced_size /* 0=default */)
{
    // case: we're asking for the same range as the previous one (for example, subsequent line in a sorted SAM)
    if (vb && vb->prev_range && vb->prev_range_chrom_node_index == vb->chrom_node_index && vb->prev_range_range_i == range_i)
        return vb->prev_range;

    // sanity checks
    ASSERT0 (vb->chrom_name, "Error in ref_get_range: vb->chrom_name=NULL");

    uint32_t save_chrom_index = vb->chrom_node_index;

    if (!ref_flag_reading_reference && (flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE)) {

        // if the chrom is not in the reference and it is numeric only, attempt to change it eg "22"->"chr22"
        if (vb->chrom_node_index >= contig_words.len)
            vb->chrom_node_index = ref_get_alternative_chrom_name (vb); // change temporarily just for ref_range_id_by_word_index()
        
        ASSERT (vb->chrom_node_index < contig_words.len, "Error in ref_get_range: chrom \"%.*s\" is not found in the reference file %s",
                MIN (vb->chrom_name_len, 100), vb->chrom_name, ref_filename);
    }

    // one_chrom_ranges_buf is an array of indeces into ranges - one index per range of size REF_NUM_SITES_PER_RANGE
    // note: in case of REF_INTERNAL, we can get hash conflicts, but in case of REF_EXTERNAL, it is guaranteed that there are no conflicts
    uint32_t range_id = (flag_reference == REF_INTERNAL) ? ref_range_id_by_hash (vb, range_i)
                                                         : ref_range_id_by_word_index (vb, range_i);

    vb->chrom_node_index = save_chrom_index; // restore

    ASSERT (range_id < REF_NUM_SITES_PER_RANGE, "Error in ref_get_range: range_id=%u expected to be smaller than %u", range_id, REF_NUM_SITES_PER_RANGE);
    Range *range = ENT (Range, ranges, range_id);

    // if ref is set, it is guaranteed that the range is good to go. This way, we can complete this function most of the
    // time without needing to lock a mutex
    char *ref = __atomic_load_n (&range->ref, __ATOMIC_RELAXED); // we load atomically, as we don't have the mutex locked
    if (ref || mode == RGR_MUST_EXIST) goto found;

    // if not - we lock the mutex and test again (it might have changed between the first and second test)
    mutex_lock (range->mutex);

    // case: it has been initialized since we first tested above 
    if (range->ref) {
        mutex_unlock (range->mutex);
        goto found;
    }

    range->range_i          = range_i;
    range->num_set          = 0;
    range->first_pos        = MAX_POS;
    range->last_pos         = 0;
    range->chrom_name_len   = vb->chrom_name_len;
    range->ref_size         = reduced_size ? reduced_size : REF_NUM_SITES_PER_RANGE;

    // case REF_INTERNAL: chrom_name points into vb->txt_data, so we make a copy of it so it survives the VB
    if (flag_reference == REF_INTERNAL) {
        range->chrom_name   = malloc (vb->chrom_name_len);
        memcpy ((char *)range->chrom_name, vb->chrom_name, vb->chrom_name_len);
    }
    // case external reference: chrom_name points into contig_dict or z_file->contexts[chrom].dict, so we just copy the pointer
    else
        range->chrom_name   = vb->chrom_name;

    ref = (flag_reference == REF_INTERNAL) ? calloc (range->ref_size + 1, 1) // +1 bc reference starts at pos=1
                                           : malloc (range->ref_size + 1);   // (no need to zero for external references)

    __atomic_store_n (&range->ref, ref, __ATOMIC_RELAXED); // the very last thing - store atomically as readers might not have mutex locked

    mutex_unlock (range->mutex);
    goto finalize;

found:
    ASSERT (mode != RGR_MUST_BE_NEW, "Error in ref_get_range: range for chrom=%.*s range_i=%u is unexpectedly already allocated",
            vb->chrom_name_len, vb->chrom_name, range_i);

    // check for hash conflict (can only happen in REF_INTERNAL) - range properties are immutable, we can safety test them
    if (flag_reference == REF_INTERNAL &&
        (range->range_i != range_i || vb->chrom_name_len != range->chrom_name_len || memcmp (vb->chrom_name, range->chrom_name, vb->chrom_name_len)))
        range = NULL;  // no soup for you

finalize:
    if (vb) {
        vb->prev_range = range;
        vb->prev_range_chrom_node_index = vb->chrom_node_index;
        vb->prev_range_range_i = range_i;
    }
    
    return range;
}

// ----------------------------------------------
// Compressing ranges into SEC_REFERENCE sections
// ----------------------------------------------

static void ref_copy_one_compressed_section (File *ref_file, const RAEntry *ra, SectionListEntry **sl)
{
    // get section list entry from ref_file_section_list - which will be used by zfile_read_section to seek to the correct offset
    while (*sl < AFTERENT (SectionListEntry, ref_file_section_list) && 
           !((*sl)->vblock_i == ra->vblock_i && (*sl)->section_type == SEC_LOCAL && (*sl)->dict_id.num == dict_id_FASTA_SEQ)) 
        (*sl)++;

    ASSERT (*sl < AFTERENT (SectionListEntry, ref_file_section_list), "Error in ref_copy_one_compressed_section: cannot find FASTA_SEQ of vb_i=%u in section list of reference file", ra->vblock_i);

    static Buffer ref_seq_section = EMPTY_BUFFER;

    SAVE_FLAG (flag_show_headers);
    zfile_read_section (ref_file, evb, ra->vblock_i, NO_SB_I, &ref_seq_section, "ref_seq_section", sizeof (SectionHeaderCtx), SEC_LOCAL, *sl);
    RESTORE_FLAG (flag_show_headers);

    SectionHeaderCtx *ref_seq_section_header = (SectionHeaderCtx *)ref_seq_section.data;

    SectionHeaderReference z_ref_section_header;
    memset (&z_ref_section_header, 0, sizeof(z_ref_section_header)); // safety

    z_ref_section_header.h.vblock_i              = 0; // we don't belong to any VB and there is no encryption of external ref
    z_ref_section_header.h.section_i             = 0; // we don't belong to any VB
    z_ref_section_header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    z_ref_section_header.h.section_type          = SEC_REFERENCE;
    z_ref_section_header.h.data_uncompressed_len = ref_seq_section_header->h.data_uncompressed_len; // already big endian
    z_ref_section_header.h.data_compressed_len   = ref_seq_section_header->h.data_compressed_len;
    z_ref_section_header.h.compressed_offset     = BGEN32 (sizeof(z_ref_section_header));
    z_ref_section_header.h.sec_compression_alg   = ref_seq_section_header->h.sec_compression_alg;
    z_ref_section_header.chrom_word_index        = BGEN32 (ra->chrom_index);
    z_ref_section_header.first_pos               = BGEN64 (ra->min_pos);
    z_ref_section_header.last_pos                = BGEN64 (ra->max_pos); 

    const char *body_start = ref_seq_section.data + BGEN32 (ref_seq_section_header->h.compressed_offset);
    uint32_t body_len = ref_seq_section.len       - BGEN32 (ref_seq_section_header->h.compressed_offset);

    // Write header and body of the reference to z_file
    // Note on encryption: reference sections originating from an external reference are never encrypted - not
    // by us here, and not in the source reference fasta (because with disallow --make-reference in combination with --password)
    START_TIMER;
    file_write (z_file, &z_ref_section_header, sizeof (z_ref_section_header));
    file_write (z_file, body_start, body_len); // write body (already compressed from the reference fasta genozip)
    COPY_TIMER (evb->profile.write);

    // "manually" add the reference section to the section list - normally it is added in comp_compress()
    sections_add_to_list (evb, &z_ref_section_header.h);

    z_file->disk_so_far += sizeof (z_ref_section_header) + body_len;   // length of GENOZIP data writen to disk

    buf_free (&ref_seq_section);
}

static void ref_copy_compressed_sections_from_reference_file (void)
{
    SectionListEntry *sl = FIRSTENT (SectionListEntry, ref_file_section_list);
    ARRAY (RAEntry, ra, ref_file_ra);

    File *ref_file = file_open (ref_filename, READ, Z_FILE, DT_FASTA);

    // note: in a FASTA file compressed with --make-reference, there is exactly one RA per VB (a contig or part of a contig)
    // we go one RA at a time and:
    // 1. If the entire RA is covered by accessed ranges - we copy the compressed reference section directly from the ref FASTA
    // 2. If we copied from the FASTA, we mark those region covered by the RA as "unaccessed", so that we don't compress it later
    for (uint32_t i=0; i < ref_file_ra.len; i++) {
        uint32_t first_range_i = (uint32_t)(ra[i].min_pos / REF_NUM_SITES_PER_RANGE);
        uint32_t last_range_i = (uint32_t)(ra[i].max_pos / REF_NUM_SITES_PER_RANGE);

        bool all_ranges_accessed = true;
        for (uint32_t range_i=first_range_i ; range_i <= last_range_i ; range_i++) {
            
            // check this range was used by this file
            evb->chrom_node_index = ra[i].chrom_index;
            uint32_t range_id = ref_range_id_by_word_index (evb, range_i);
            Range *r = ENT (Range, ranges, range_id);

            if (!r->is_accessed) {
                all_ranges_accessed = false;
                break;
            }
        }

        // if this entire RA is covered, just copy the corresponding FASTA section to our file, and
        // mark all the ranges as is_accessed=false indicated they don't need to be compressed individually
        if (all_ranges_accessed) {
            ref_copy_one_compressed_section (ref_file, &ra[i], &sl);

            for (uint32_t range_i=first_range_i ; range_i <= last_range_i ; range_i++) {
                evb->chrom_node_index = ra[i].chrom_index;
                uint32_t range_id = ref_range_id_by_word_index (evb, range_i);
                Range *r = ENT (Range, ranges, range_id);

                // if the remaining set area of range is entirely included in this ra - we don't need to compress it at all
                if (r->first_pos >= ra[i].min_pos && r->last_pos <= ra[i].max_pos)
                    r->is_accessed = false;

                // if the remaining set area of range is only patially included in this ra - shrink the set area to only the remaining.
                // note: that ref contigs always start at pos 1 and the RAs are in order
                else if (r->last_pos > ra[i].max_pos)
                    r->first_pos = ra[i].max_pos + 1;

                else if (r->first_pos < ra[i].min_pos)
                    r->last_pos = ra[i].min_pos - 1;
            }
        }
    }

    file_close (&ref_file, false);
}

static uint32_t ref_compact_ref (Range *r)
{
    uint32_t first_offset = (uint32_t)(r->first_pos - r->range_i * REF_NUM_SITES_PER_RANGE);
    uint32_t last_offset  = (uint32_t)(r->last_pos  - r->range_i * REF_NUM_SITES_PER_RANGE);

    uint32_t gap_start=first_offset, stretch_start=first_offset, i=first_offset;
    char *next = r->ref;
    do {
        while (i <= last_offset && !r->ref[i]) i++; // first first non-zero (or end of ref)
        uint32_t gap_len = i - gap_start;
        
        // if the gap between this start and the previous stretch is big enough to justify it, we enter the length
        if (gap_len >= 32) {
            *(next++) = gap_len         & 0xff; // LSB (3 bytes are enough as our ranges are 1MB long)
            *(next++) = (gap_len >> 8)  & 0xff;
            *(next++) = (gap_len >> 16) & 0xff; // MSB
            *(next++) = '\t';                   // gap indicator (appears last bc in piz we scan backwards)

            stretch_start = i;
        }
        else {} // no change in stretch_start 

        while (i <= last_offset && r->ref[i]) i++; // find the first 0 (or end of ref)
        uint32_t stretch_len = i - stretch_start; // this includes the preceding short gap (< 32 chars) if there is one

        memmove (next, &r->ref[stretch_start], stretch_len);
        next += stretch_len;
        gap_start = stretch_start = i;
    
    } while (i <= last_offset);

    return (uint32_t)(next - r->ref);
}

static void ref_compress_one_range (VBlockP vb, const char *data,
                                    uint32_t uncompressed_len, uint32_t section_i, uint32_t chrom_word_index, 
                                    int64_t first_pos, int64_t last_pos)
{
    SectionHeaderReference header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.vblock_i              = 0; // reference doesn't belong to any VB (consequence: in case of encryption of internal reference, all reference sections are encrypted with the same AES key, that's fine.)
    header.h.section_i             = 0; 
    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = SEC_REFERENCE;
    header.h.data_uncompressed_len = BGEN32 (uncompressed_len);
    header.h.compressed_offset     = BGEN32 (sizeof(header));
    header.h.sec_compression_alg   = COMP_LZMA;
    header.chrom_word_index        = BGEN32 (chrom_word_index);
    header.first_pos               = BGEN64 (first_pos);
    header.last_pos                = BGEN64 (last_pos); 

    vb->z_data.name  = "z_data"; // comp_compress requires that these are pre-set
    vb->z_data.param = vb->vblock_i;
    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, data, NULL);
}

// thread entry point - compresses as many ranges as fit in global_max_memory_per_vb
static void ref_compress_some_ranges (VBlockP vb)
{
    uint32_t total_uncompressed = 0;
    const Range *after_ranges = AFTERENT (const Range, ranges);

    // compress ranges up to global_max_memory_per_vb of uncompressed_len of data
    for (; vb->range < after_ranges && (total_uncompressed + vb->range->uncompressed_len <= global_max_memory_per_vb);
         vb->range++) {

        if (!vb->range->is_accessed) continue; // skip unaccessed range

        int64_t first_index = vb->range->is_compacted ? 0 : vb->range->first_pos - vb->range->range_i * REF_NUM_SITES_PER_RANGE;
    
        // at this point, after the compute threads are all done, z_file->contexts is immutable, so we can safely access it
        MtfNode *node;
        hash_global_get_entry (&z_file->contexts[DTFZ(chrom)], vb->range->chrom_name, vb->range->chrom_name_len, HASH_READ_ONLY, &node);

        ref_compress_one_range (vb, &vb->range->ref[first_index], vb->range->uncompressed_len, vb->z_next_header_i++, 
                                node->word_index.n, vb->range->first_pos, vb->range->last_pos);

        total_uncompressed += vb->range->uncompressed_len;

        //fprintf (stderr, "Compressing range_i=%u with chrom=%.*s pos=[%u,%u]\n", (uint32_t)(vb->range - (Range *)ranges.data),
        //         vb->range->chrom_name_len, vb->range->chrom_name, (uint32_t)vb->range->first_pos, (uint32_t)vb->range->last_pos);
    }

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

static void ref_output_one_range (VBlockP vb) 
{ 
    zip_output_processed_vb (vb, &vb->section_list_buf, false, PD_REFERENCE_DATA); 
}

// compress the reference - one section at the time, using Dispatcher to do them in parallel 
static uint32_t next_range_i=0;
static void ref_prepare_range_for_compress (VBlockP vb)
{
    // find next occupied range
    ARRAY (Range, rng, ranges);
    vb->range = &rng[next_range_i]; // first range to compress
    uint32_t total_uncompressed = 0;

    if (next_range_i == REF_NUM_RANGES) {
        vb->ready_to_dispatch = false; // we're done
        return;
    }

    // ref_compress_some_ranges will consume data up to global_max_memory_per_vb, so skip ranges up to this amount
    for (; next_range_i < REF_NUM_RANGES ; next_range_i++) {

        Range *r = &rng[next_range_i];
        if (!r->is_accessed) continue; // skip ranges not used by this file

        if (!r->uncompressed_len) { // possibly already handled in previous VB, if it was the last range and didn't make it due to going oversize
        
            // if this ref is sparse (eg due to a low coverage sample), we compact it (only happens for REF_INTERNAL, bc for external r->num_set is always the entire range)
            if ((double)(r->last_pos - r->first_pos + 1) / (double)r->num_set >= 1.5) { 
                r->uncompressed_len = ref_compact_ref (r);
                r->is_compacted = true;
            }
            else {
                r->uncompressed_len = r->last_pos - r->first_pos + 1;                
                r->is_compacted = false;
            }
        }

        total_uncompressed += r->uncompressed_len;

        if (total_uncompressed > global_max_memory_per_vb) break; // break before incrementing next_range_i
    }
    vb->ready_to_dispatch = true;
}

void ref_compress_ref (void)
{
    if (!buf_is_allocated (&ranges)) return;

    // copy already-compressed SEQ sections from the FASTA genozip reference, but only such sections that are entirely
    // covered by ranges with is_accessed=true. we mark these ranges affected as is_accessed=false.
    if (flag_reference == REF_EXT_STORE)
         ref_copy_compressed_sections_from_reference_file ();

    // proceed to compress all ranges that remain with is_accessed=true
    next_range_i=0; // can be initialized multiple times if we're compressing multiple non-concatenated SAMs - but not concurrently

    uint32_t num_ref_sections = 
        dispatcher_fan_out_task ("Writing reference to genozip file", false, 
                                 ref_prepare_range_for_compress, 
                                 ref_compress_some_ranges, 
                                 ref_output_one_range);

    // SAM require at least one reference section, but if the SAM is unaligned, there will be none - create one empty section
    if (z_file->data_type == DT_SAM && !num_ref_sections)
        ref_compress_one_range (evb, NULL, 0, 0, 0, 0, 0);
}

void ref_set_reference (const char *filename)
{
    ref_filename = filename;
}

void ref_set_md5 (Md5Hash md5)
{
    ref_md5 = md5;
}

// ZIP & PIZ: import external reference
void ref_read_external_reference (void)
{
    ASSERT0 (ref_filename, "Error: ref_filename is NULL");

    z_file = file_open (ref_filename, READ, Z_FILE, DT_FASTA);    
    const char *basename = file_basename (ref_filename, false, "(reference)", NULL, 0);
        
    // save and reset some globals that after pizzing FASTA
    int save_flag_test         = flag_test        ; flag_test        = 0;
    int save_flag_md5          = flag_md5         ; flag_md5         = 0;
    int save_flag_show_time    = flag_show_time   ; flag_show_time   = 0;
    int save_flag_show_memory  = flag_show_memory ; flag_show_memory = 0;
    int save_flag_no_header    = flag_no_header   ; flag_no_header   = 0;
    int save_flag_header_one   = flag_header_one  ; flag_header_one  = 0;
    int save_flag_header_only  = flag_header_only ; flag_header_only = 0;
    char *save_flag_grep       = flag_grep        ; flag_grep        = 0;
    ref_primary_command        = command          ; command          = UNZIP;
    ref_flag_reading_reference = true; // tell fasta.c that this is a reference
    
    bool piz_successful = piz_dispatcher (basename, true, false);

    // recover globals
    ref_flag_reading_reference = false;
    command          = ref_primary_command;
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

// case 1: in case of ZIP with --reference or --REFERENCE, called by ref_consume_ref_fasta_global_area during piz_read_global_area of FASTA ref file
// case 2: in case of ZIP of SAM using internal reference - called from sam_zip_initialize
// note: ranges allocation must be called by the I/O thread as it adds a buffer to evb buf_list
void ref_zip_initialize_ranges (void)
{
    buf_alloc (evb, &ranges, REF_NUM_RANGES * sizeof (Range), 1, "ranges", 0); 
    ranges.len = REF_NUM_RANGES;
    buf_zero (&ranges);

    // create all 1M muteces, hopefully the OS can handle it
    for (unsigned i=0; i < REF_NUM_RANGES; i++)
        pthread_mutex_init (&ENT (Range, ranges, i)->mutex, NULL);
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

// ZIP & PIZ: called by piz_read_global_area during PIZ of the external reference FASTA file
void ref_consume_ref_fasta_global_area (void)
{
    // copy data from the reference FASTA's CONTIG context, so it survives after we finish reading the reference and close z_file
    Context *fasta_contig_ctx = &z_file->contexts[DTFZ(chrom)];

    ASSERT (buf_is_allocated (&fasta_contig_ctx->dict) && buf_is_allocated (&fasta_contig_ctx->word_list),
            "Error: cannot use %s as a reference as it is missing a CONTIG dictionary", z_name);

    buf_copy (evb, &contig_dict,  &fasta_contig_ctx->dict,      1               , 0, 0, "contig_dict",  0);
    buf_copy (evb, &contig_words, &fasta_contig_ctx->word_list, sizeof (MtfWord), 0, 0, "contig_words", 0);
    
    random_access_max_pos_of_chrom (0); // initialize if not already initialized
    buf_copy (evb, &contig_maxes, &z_file->ra_max_pos_by_chrom, sizeof (int64_t), 0, 0, "contig_maxes", 0);
    buf_copy (evb, &ref_file_ra, &z_file->ra_buf, sizeof (RAEntry), 0, 0, "ref_file_ra", 0);

    buf_copy (evb, &ref_file_section_list, &z_file->section_list_buf, sizeof (SectionListEntry), 0, 0, "ref_file_section_list", 0);

    // contig_words_sorted_index - an array of uint32 of indexes into contig_words - sorted by alphabetical order of the snip in contig_dict
    buf_alloc (evb, &contig_words_sorted_index, sizeof(uint32_t) * contig_words.len, 1, "contig_words_sorted_index", 0);
    for (uint32_t i=0; i < contig_words.len; i++)
        NEXTENT (uint32_t, contig_words_sorted_index) = i;

    qsort (contig_words_sorted_index.data, contig_words_sorted_index.len, sizeof(uint32_t), ref_sort_words_alphabetically);

    if (ref_primary_command == ZIP) 
        ref_zip_initialize_ranges(); // in ZIP, we allocate 1M range, 1MB each (both internal and external ref) 
    
    else { // in PIZ using an external reference, we allocate 1 range per contig (chrom)
        buf_alloc (evb, &ranges, contig_words.len * sizeof (Range), 1, "ranges", 0); 
        ranges.len = contig_words.len;
        buf_zero (&ranges);

        ARRAY (MtfWord, words, contig_words);

        for (uint32_t i=0; i < ranges.len; i++) {
            Range *r = ENT (Range, ranges, i);
            r->chrom_name     = ENT (char, contig_dict, words[i].char_index);
            r->chrom_name_len = words[i].snip_len;
            r->first_pos      = 1;
            r->last_pos       = random_access_max_pos_of_chrom (i);
            r->num_set        = r->last_pos;
            r->ref_size       = r->last_pos + 1;
            r->ref            = malloc (r->ref_size);

            ASSERT (r->ref, "Error in ref_consume_ref_fasta_global_area: failed to allocate %"PRId64" bytes for reference chrom %.*s",
                    r->ref_size, r->chrom_name_len, r->chrom_name);
        }
    }
}

// ZIP & PIZ: used when zipping or pizing the target file, after the external reference is read
int64_t ref_max_pos_of_chrom (uint32_t chrom_word_index)
{
    if (command == ZIP) {
        ASSERT (chrom_word_index < contig_maxes.len, "Error in ref_max_pos_of_chrom: chrom_word_index=%u out of range, contig_maxes.len=%u",
                chrom_word_index, (uint32_t)contig_maxes.len);

        ASSERT0 (flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE, "Error in ref_max_pos_of_chrom: function can only be called with external reference");

        return *ENT (int64_t, contig_maxes, chrom_word_index);
    }
    else { // UNZIP
        ASSERT (chrom_word_index < ranges.len, "Error in ref_max_pos_of_chrom: chrom_word_index=%u out of range, ranges.len=%u",
                chrom_word_index, (uint32_t)ranges.len);

        return ENT (RAEntry, ranges, chrom_word_index)->max_pos;
    }
}

/*                // case: ref exists but is different than SEQ - replace base with 1, 2 or 3 - 1 being the most frequent SNP etc. This reduces the alphabet from
                // ~4 letters to ~3, and differentiates them more probability-wise (this saves about 1% of SEQ)
                TO DO: handle N   
                else if (ref_value == 'A') {
                    if      (seq[i] == 'G') seq[i] = '1'; // 633K A->G in Divon's VCF
                    else if (seq[i] == 'C') seq[i] = '2'; // 162K A->C in Divon's VCF
                    else if (seq[i] == 'T') seq[i] = '3'; // 141K A->T in Divon's VCF
                }
                
                else if (ref_value == 'C') {
                    if      (seq[i] == 'T') seq[i] = '1'; // 660K C->T in Divon's VCF
                    else if (seq[i] == 'G') seq[i] = '2'; // 168K C->G in Divon's VCF
                    else if (seq[i] == 'A') seq[i] = '3'; // 166K C->A in Divon's VCF
                }
                
                else if (ref_value == 'G') {
                    if      (seq[i] == 'A') seq[i] = '1'; // 660K A->G in Divon's VCF
                    else if (seq[i] == 'C') seq[i] = '2'; // 168K A->C in Divon's VCF
                    else if (seq[i] == 'T') seq[i] = '3'; // 166K A->T in Divon's VCF (yes! exactly the same histogram as C)
                }
                
                else if (ref_value == 'T') {
                    if      (seq[i] == 'C') seq[i] = '1'; // 633K A->G in Divon's VCF
                    else if (seq[i] == 'G') seq[i] = '2'; // 161K A->C in Divon's VCF
                    else if (seq[i] == 'A') seq[i] = '3'; // 140K A->T in Divon's VCF
                }
                // case: ref exists, but is not A,C,G or T (eg N) - do nothing - SEQ is unchanged for this site
     */           
                // case: ref exists, but is different that seq - store the seq base in nonref_ctx and indicate '.'
