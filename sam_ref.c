
// ------------------------------------------------------------------
//   sam_ref.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <pthread.h>
#include "sam_private.h"
#include "buffer.h"
#include "strings.h"
#include "dict_id.h"
#include "dispatcher.h"
#include "zip.h"
#include "zfile.h"
#include "endianness.h"
#include "hash.h"
#include "random_access.h"
#include "seg.h"
#include "piz.h"

typedef struct Range {
    char *ref; // once malloced, never changes
    pthread_mutex_t mutex;
    char *chrom_name;
    unsigned chrom_name_len;
    uint32_t chrom_node_index;
    uint32_t range_i;
    uint32_t num_set; // number of bases set in this range
    uint32_t first_pos, last_pos;
    uint32_t uncompressed_len; // final length of range ready for compression
    bool is_compacted;
} Range;

// reference sequences - one per range of 1MB. ranges (chrom, pos) are mapped here with a hash function. In the rare case two unrelated ranges
// are mapped to the same entry - only the first range will have a reference, and the other ones will not. this will hurt the compression ratio,
// but not the correctness.
// Thread safety: "ref" is set atomically as the last set in initialization. If its set, then it is correct and will never change. If it appears to be
// not set yet, it is necesarry to lock the mutex and test again, and initialize if still not set.

#define NUM_RANGES (1 << 20)
#define NUM_SITES_PER_RANGE (1 << 20) // 1 MB
static Buffer ranges = EMPTY_BUFFER; // ZIP: a buffer containing a array of NUM_RANGES Range's - each allocated on demand. 

// PIZ: an array of references for each chrom (each reference in a Buffer) - indexed by word_index as in the Chrom context. 
static Buffer *refs = NULL; 
static unsigned refs_len = 0;

// called between files when compressing or decompressing multiple genozip files
void sam_ref_cleanup_memory (VBlock *vb)
{
    if (vb != evb) return;

    if (refs) {
        for (unsigned i=0; i < refs_len; i++) buf_destroy (&refs[i]);
        free (refs);
        refs = NULL;
        refs_len = 0;
    }

    buf_free (&ranges);
}


// import external reference
void sam_ref_import (const char *filename)
{

}

// ------------------------------------
// PIZ side
// ------------------------------------

static const char *sam_ref_get_ref (VBlockSAMP vb)
{
    int64_t pos = vb->contexts[SAM_POS].last_value;
    int32_t ref_i = vb->chrom_node_index;

    if (!buf_is_allocated (&refs[ref_i])) return NULL; // this can if entire chromosome is verbatim, eg. unaligned (pos=4) or SEQ or CIGAR are unavailable

    if (*vb->last_cigar == '*') return NULL; // no cigar, means no reference - eg unaligned sequenece

    ASSERT (pos + vb->ref_consumed <= refs[ref_i].len, "Error in sam_ref_get_ref: out of range reconstructing txt_line_i=%u: pos=%u ref_consumed=%u but refs[ref_i].len=%u",
            vb->line_i, (uint32_t)pos, vb->ref_consumed, (uint32_t)refs[ref_i].len);

    return ENT (const char, refs[ref_i], pos);
}

// PIZ: SEQ contains : 
// '-' - data should be taken from the reference (-) 
// '.' - data should be taken from SQnonref.local (.)
// other - verbatim (this happens if there is no reference, eg unaligned BAM)
void sam_ref_reconstruct (VBlock *vb_, Context *seq_ctx)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    ASSERT0 (seq_ctx, "Error in sam_ref_reconstruct: seq_ctx is NULL");

    if (piz_is_skip_section (vb, SEC_LOCAL, seq_ctx->dict_id)) return; // if case we need to skip the SEQ field (for the entire file)

    Context *nonref_ctx = mtf_get_ctx (vb, (DictId)dict_id_SAM_SQnonref);
    const char *nonref = &nonref_ctx->local.data[nonref_ctx->next_local]; // possibly, this VB has no nonref (i.e. everything is ref), in which cse nonref would be an invalid pointer. That's ok, as it will not be accessed.
    const char *nonref_start = nonref;
    const char *seq = &seq_ctx->local.data[seq_ctx->next_local];
    unsigned subcigar_len = 0;
    char cigar_op;
    
    // case where seq is '*' (rewritten as ' ' by the zip callback)
    if (seq[0] == ' ') {
        RECONSTRUCT1 ('*');
        seq_ctx->last_value = seq_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence
        seq_ctx->next_local++; // only 1
        return;
    }

    const char *ref = sam_ref_get_ref (vb); // pointer to ref @ last chrom and pos. ref[0] is the base at POS.

    unsigned seq_consumed=0, ref_consumed=0;
    while (seq_consumed < vb->seq_len) {
        
        if (!subcigar_len) {
            subcigar_len = strtod (vb->last_cigar, (char **)&vb->last_cigar); // get number and advance next_cigar
            cigar_op = *(vb->last_cigar++);

            // case: Deletion or Skipping - skip some of the reference
            if (cigar_op == 'D' || cigar_op == 'N') {
                ref_consumed += subcigar_len;
                subcigar_len = 0;
                continue;
            } 

            // case: hard clipping - just ignore this subcigar
            if (cigar_op == 'H' || cigar_op == 'P') {
                subcigar_len = 0;
                continue;
            }
        }

        char c = seq[seq_consumed++];

        if (c == '-') {
            ASSERT0 (ref, "SEQ shows -, but ref is unavailable (eg bc entire chromosome is POS=0 or CIGAR=* or SEQ=*");
            RECONSTRUCT1 (ref[ref_consumed++]); 
        }

        else if (c == '.') {
            RECONSTRUCT1 (*nonref);
            nonref++; // consume non-ref even if not reconstructing this seq

            // advance ref if this is a SNP (but not in case of 'I' or 'S')
            if (cigar_op == 'M' || cigar_op == '=' || cigar_op == 'X') ref_consumed++; 
        }
        
        // in case of SAM having CIGAR='*' or POS=0 (eg unaligned BAM), we just copy verbatim
        else RECONSTRUCT1 (c); 


        subcigar_len--;
    }

    seq_ctx->last_value = seq_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence
    seq_ctx->next_local += vb->seq_len;
    
    nonref_ctx->next_local += (uint32_t)(nonref - nonref_start);
}


static void sam_ref_uncompress_one_range (VBlockP vb)
{
    SectionHeaderSAMReference *header = (SectionHeaderSAMReference *)vb->z_data.data;

    uint32_t ref_i      = BGEN32 (header->chrom_word_index);
    uint32_t first_pos  = BGEN32 (header->first_pos);
    uint32_t last_pos   = BGEN32 (header->last_pos);
    uint32_t uncomp_len = BGEN32 (header->h.data_uncompressed_len);

    Context *ctx = &z_file->contexts[SAM_RNAME];
    const MtfWord *word = ENT (const MtfWord, ctx->word_list, ref_i);

    ASSERT (last_pos < refs[ref_i].len, "Error in sam_ref_uncompress_one_range: ref range out of bounds for chrom=%.*s: first_pos=%u last_pos=%u but ref_len=%u",
            word->snip_len, ENT (char, ctx->dict, word->char_index), first_pos, last_pos, (uint32_t)refs[ref_i].len);

    char *uncompressed_data = ENT (char, refs[ref_i], first_pos);

    zfile_uncompress_section (vb, (SectionHeaderP)header, uncompressed_data, NULL, SEC_SAM_REFERENCE);

    // case: ref has been compacted - we uncompact it
    // note: the unused gap areas of reference contain undefined data and are not used. we are careful not to zero
    // them as this way the OS lazy-allocation algorithm doesn't actually consume memory blocks that are entirely unused
    if (uncomp_len < last_pos - first_pos + 1) {

        // traverse backwards and "stretch" the reference 
        char *next_ref = ENT (char, refs[ref_i], last_pos); // the last char of the uncompacted ref
        for (char *next_com = &uncompressed_data[uncomp_len-1]; next_com >= uncompressed_data; next_com--) { 

            if (*next_com == '\t') {
                uint32_t gap_len = (uint32_t)(uint8_t)next_com[-3] | ((uint32_t)(uint8_t)next_com[-2] << 8) | ((uint32_t)(uint8_t)next_com[-1] << 16);
                next_ref -= gap_len;
                next_com -= 3; // skip gap length
            }
            else {
                ASSERT (next_ref >= refs[ref_i].data, "Error in sam_ref_uncompress_one_range: next_ref out of range when uncompacting chrom=%.*s, first_pos=%u last_pos=%u",
                        word->snip_len, ENT (char, ctx->dict, word->char_index), first_pos, last_pos);

                *(next_ref--) = *next_com;
            }
        } 
    }

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. 
}

static bool sam_ref_is_ref_range_section (SectionType st) { return st == SEC_SAM_REFERENCE; }

static SectionListEntry *sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
static uint32_t ref_range_cursor = 0;

static void sam_ref_read_one_range (VBlockP vb)
{
    if (sections_get_next_section_of_type (&sl_ent, &ref_range_cursor, sam_ref_is_ref_range_section)) {
        
        zfile_read_section (evb, sl_ent->vblock_i, NO_SB_I, &vb->z_data, "z_data", sizeof(SectionHeaderSAMReference), sl_ent->section_type, sl_ent);    
        
        // allocate memory for entire chrom reference if this is the first range of this chrom
        uint32_t chrom_word_index = BGEN32 (((SectionHeaderSAMReference *)vb->z_data.data)->chrom_word_index);
        if (!buf_is_allocated (&refs[chrom_word_index])) {
            refs[chrom_word_index].len = 1 + (uint64_t)random_access_get_vb0_max_pos (chrom_word_index); // cast to 64b to allow len=2^32 if case of MAX_POS
            // note: the OS lazy-allocation will only actually allocate memory blocks that are actually used 
            buf_alloc (evb, &refs[chrom_word_index], refs[chrom_word_index].len, 1, "refs", chrom_word_index);
        }

        vb->ready_to_dispatch = true;
    }
}

void sam_ref_read_all_ranges (void)
{
    ASSERT0 (!refs, "Error in sam_ref_read_all_ranges: expecting refs to be NULL");
    
    unsigned num_chroms = z_file->contexts[SAM_RNAME].word_list.len;
    refs = calloc (num_chroms, sizeof (Buffer));
    
    sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
    ref_range_cursor = 0;

    // decompress reference using Dispatcher
    dispatcher_fan_out_task ("Internal reference sequence", sam_ref_read_one_range, sam_ref_uncompress_one_range, NULL);
}

// ------------------------------------
// ZIP side
// ------------------------------------

// this must be called by the I/O thread as it adds buffers to evb buf_list
void sam_ref_zip_initialize (void)
{
    buf_alloc (evb, &ranges, NUM_RANGES * sizeof (Range), 1, "ranges", 0); // Never freed or realloced.
    buf_zero (&ranges);

    // create all 1M muteces, hopefully the OS can handle it
    for (unsigned i=0; i < NUM_RANGES; i++)
        pthread_mutex_init (&ENT (Range, ranges, i)->mutex, NULL);
}

static inline uint32_t sam_ref_hash_range (const char *chrom_name, unsigned chrom_name_len, uint32_t range_i)
{
    ASSERT0 (chrom_name_len > 0, "Error in sam_ref_hash_range: chrom_name_len==0");

    // step 1: get number embedded in the chrom name
    uint32_t n=0;
    for (unsigned i=0; i < chrom_name_len; i++)
        if (IS_DIGIT (chrom_name[i])) 
            n = n*10 + (chrom_name[i] - '0');

    // if there are no digits, take the last 4 characters
    if (!n)
        n = (                       (((uint32_t)chrom_name[chrom_name_len-1]) ^ 0x1f))            ^ 
            (chrom_name_len >= 2 ? ((((uint32_t)chrom_name[chrom_name_len-2]) ^ 0x1f) << 3)  : 0) ^
            (chrom_name_len >= 3 ? ((((uint32_t)chrom_name[chrom_name_len-3]) ^ 0x1f) << 4)  : 0) ^ 
            (chrom_name_len >= 4 ? ((((uint32_t)chrom_name[chrom_name_len-4]) ^ 0x1f) << 5)  : 0) ;

    // step: calculate the hash - 10 bit for n, 10 bit for range_i
    uint32_t chr_component = (chrom_name_len <= 6) ? (n & 0x1f) : ((n % 896) + 128); // first 128 entries are reserved for the major chromosomes, heuristically identified as name length 6 or less

    uint32_t value = chr_component | (((range_i & 0x3ff) ^ ((n*3) & 0x3ff)) << 10);
    return value; 
}

// Allocated and initializes the ref and mutex buffers for the given chrom/range
// Note: this is the only function that accesses the ranges buffer directly during the compute threads stage
static inline Range *sam_ref_get_range (const char *chrom_name, unsigned chrom_name_len, uint32_t range_i)
{
    // one_chrom_ranges_buf is an array of indeces into ranges - one index per range of size NUM_SITES_PER_RANGE
    uint32_t range_hash = sam_ref_hash_range (chrom_name, chrom_name_len, range_i);
    ASSERT (range_hash < NUM_SITES_PER_RANGE, "Error in sam_ref_get_range: range_hash=%u expected to be smaller than %u", range_hash, NUM_SITES_PER_RANGE);

    Range *range = ENT (Range, ranges, range_hash);

    // if ref is set, it is guaranteed that the range is good to go. This way, we can complete this function most of the
    // time without needing to lock a mutex
    char *ref = __atomic_load_n (&range->ref, __ATOMIC_RELAXED); // we load atomically, as we don't have the mutex locked
    if (ref) goto found;

    // if not - we lock the mutex and test again (it might have changed between the first and second test)
    pthread_mutex_lock (&range->mutex);

    // case: it has been initialized since we first tested above 
    if (range->ref) {
        pthread_mutex_unlock (&range->mutex);
        goto found;
    }

    range->range_i          = range_i;
    range->num_set          = 0;
    range->first_pos        = MAX_POS;
    range->last_pos         = 0;
    range->chrom_node_index = NIL;
    range->chrom_name_len   = chrom_name_len;
    range->chrom_name       = malloc (chrom_name_len);
    memcpy (range->chrom_name, chrom_name, chrom_name_len);

    ref = calloc (NUM_SITES_PER_RANGE, 1); // first prepare...
    __atomic_store_n (&range->ref, ref, __ATOMIC_RELAXED); // the very last thing - store atomically as readers might not have mutex locked

    pthread_mutex_unlock (&range->mutex);
    return range;

found:
    // range properties are immutable, we can safety test them
    if (range->range_i == range_i && range->chrom_name_len == chrom_name_len && !memcmp (range->chrom_name, chrom_name, chrom_name_len)) 
        return range; // we found a range
    else 
        return NULL;  // sorry, no soup for you, hash conflict
}

void sam_ref_refy_one_seq (VBlock *vb, char *seq, uint32_t seq_len, uint32_t pos, const char *cigar, 
                           const char *chrom_name, unsigned chrom_name_len)
{
    uint32_t range_i = (uint32_t)(pos / NUM_SITES_PER_RANGE);
    bool range_mutex_is_locked = false;

    Context *nonref_ctx = mtf_get_ctx (vb, (DictId)dict_id_SAM_SQnonref);
    buf_alloc (vb, &nonref_ctx->local, MAX (nonref_ctx->local.len + seq_len, vb->lines.len * seq_len / 5), CTX_GROWTH, "mtf->local", nonref_ctx->did_i);
    Range *range = sam_ref_get_range (chrom_name, chrom_name_len, range_i);
    
    // this range cannot be diffed against a reference, as the hash entry for this range is unfortunately already occupied by another range
    if (!range) {
        buf_add (&nonref_ctx->local, seq, seq_len);
        memset (seq, '.', seq_len);
        return; 
    }

    uint32_t next_ref  = pos - range_i * NUM_SITES_PER_RANGE;

    const char *next_cigar = cigar;
    uint32_t i=0;
    int subcigar_len=0;
    char cigar_op;

    while (i < seq_len && next_ref < NUM_SITES_PER_RANGE) {

        subcigar_len = strtod (next_cigar, (char **)&next_cigar); // get number and advance next_cigar
        cigar_op = *(next_cigar++);

        if (cigar_op == 'M' || cigar_op == '=' || cigar_op == 'X') { // alignment match or sequence match or mismatch

            ASSERT (subcigar_len > 0 && subcigar_len <= (seq_len - i), 
                    "Error in sam_ref_refy_one_seq: CIGAR %s implies seq_len longer than actual seq_len=%u", cigar, seq_len);

            while (subcigar_len && next_ref < NUM_SITES_PER_RANGE) {

                // note that our "ref" might be different than the alignment reference, therefore X and = are not an indication
                // whether it matches our ref or not. However, if it is an X we don't enter it into our ref, and we wait
                // for a read with a = or M for that site
            
                // tread saftey: if range->ref[next_ref] contains a non-zero value, then that value is final and immutable. if it does not yet contain
                // a value, we need to grab the mutex to verify and potentially set it
                char ref_value = range->ref[next_ref];

                if (!ref_value && !range_mutex_is_locked) {
                    pthread_mutex_lock (&range->mutex); // lock on our first change to this ref
                    range_mutex_is_locked = true;

                    ref_value = range->ref[next_ref]; // check again - this time we are sure that if it is 0, it will not change as we locked the mutex
                }

                if (ref_value == seq[i]) 
                    seq[i] = '-'; // ref
                
                else if (!ref_value) { // no ref yet for this site - we will be the ref
                    range->ref[next_ref] = seq[i];
                    range->num_set++;
                    seq[i] = '-'; 
                }

                else {
                    NEXTENT (char, nonref_ctx->local) = seq[i];
                    seq[i] = '.';
                } 

                subcigar_len--;
                next_ref++;
                i++;
            }
        } // end if 'M', '=', 'X'

        // for Insertion or Soft clipping - this SEQ segment doesn't align with the ref - we leave it as is 
        else if (cigar_op == 'I' || cigar_op == 'S') {

            ASSERT (subcigar_len > 0 && subcigar_len <= (seq_len - i), 
                    "Error in sam_ref_refy_one_seq: CIGAR %s implies seq_len longer than actual seq_len=%u", cigar, seq_len);

            buf_add (&nonref_ctx->local, &seq[i], subcigar_len);
            memset (&seq[i], '.', subcigar_len);
            i += subcigar_len;
        }

        // for Deletion or Skipping - we move the next_ref ahead
        else if (cigar_op == 'D' || cigar_op == 'N') {
            unsigned ref_consumed = MIN (subcigar_len, NUM_SITES_PER_RANGE - next_ref);

            next_ref += ref_consumed;
            subcigar_len -= ref_consumed;
        }

        // in other cases - Hard clippping (H) or padding (P) we do nothing
        else {} 
    }

    if (range_mutex_is_locked)
        pthread_mutex_unlock (&range->mutex);       

    // update first, last pos of this range
    if (pos < range->first_pos) range->first_pos = pos;

    uint32_t last_pos = (range_i * NUM_SITES_PER_RANGE) + next_ref - 1;
    if (last_pos > range->last_pos) range->last_pos = last_pos;

    // case: we have reached the end of the current ref range, but we still have sequence left - 
    // call recursively with remaining sequence and next ref range 
    if (i < seq_len) {

        ASSSEG (last_pos < MAX_POS, cigar, "%s: Error: reference pos, considering POS=%u and the consumed reference impied by CIGAR=%s, exceeds MAX_POS=%u",
                global_cmd, pos, cigar, MAX_POS);

        char updated_cigar[100];
        if (subcigar_len) sprintf (updated_cigar, "%u%c%s", subcigar_len, cigar_op, next_cigar);

        sam_ref_refy_one_seq (vb, seq + i, seq_len - i, (range_i+1) * NUM_SITES_PER_RANGE,
                              subcigar_len ? updated_cigar : next_cigar, 
                              chrom_name, chrom_name_len);
    }
    else // update RA of the VB with last pos of this line as implied by the CIGAR string
        random_access_update_last_pos (vb, (uint32_t)last_pos);
}

static uint32_t sam_ref_compact_ref (Range *r)
{
    uint32_t first_offset = r->first_pos - r->range_i * NUM_SITES_PER_RANGE;
    uint32_t last_offset  = r->last_pos  - r->range_i * NUM_SITES_PER_RANGE;

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

// thread entry point
static void sam_ref_compress_one_range (VBlockP vb_)
{
    VBlockSAM *vb = (VBlockSAM *)vb_;

    uint32_t first_index = vb->range->is_compacted ? 0 : vb->range->first_pos - vb->range->range_i * NUM_SITES_PER_RANGE;
 
    // at this point, after the compute threads are all done, z_file->contexts is immutable, so we can safely access it
    MtfNode *node;
    hash_global_get_entry (&z_file->contexts[SAM_RNAME], vb->range->chrom_name, vb->range->chrom_name_len, HASH_MUST_EXIST, &node);

    SectionHeaderSAMReference header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = SEC_SAM_REFERENCE;
    header.h.data_uncompressed_len = BGEN32 (vb->range->uncompressed_len);
    header.h.compressed_offset     = BGEN32 (sizeof(header));
    header.h.sec_compression_alg   = COMP_LZMA;
    header.h.section_i             = BGEN16 (vb->z_next_header_i++);
    header.chrom_word_index        = BGEN32 (node->word_index.n);
    header.first_pos               = BGEN32 (vb->range->first_pos);
    header.last_pos                = BGEN32 (vb->range->last_pos); 

    vb->z_data.name  = "z_data"; // comp_compress requires that these are pre-set
    vb->z_data.param = vb->vblock_i;
    comp_compress (vb_, &vb->z_data, false, (SectionHeader*)&header, &vb->range->ref[first_index], NULL);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

// compress the reference - one section at the time, using Dispatcher to do them in parallel 
static uint32_t next_range_i=0;
static void sam_ref_output_one_range (VBlockP vb) { zip_output_processed_vb (vb, &vb->section_list_buf, false, PD_SAM_REF_DATA); }
static void sam_ref_prepare_range_for_compress (VBlockP vb)
{
    // find next occupied range
    ARRAY (Range, rng, ranges);
    while (!rng[next_range_i].num_set && next_range_i < NUM_RANGES) next_range_i++;

    if (next_range_i < NUM_RANGES) {  // we found some data 
        Range *r = &rng[next_range_i];

        // if this ref is sparse (eg due to a low coverage sample), we compact it
        if (((double)r->last_pos - (double)r->first_pos + 1) / (double)r->num_set >= 1.5) { // case last/first to double to prevent last-first+1 overflowing uint32 in case of [0,MAX_POS]
            r->uncompressed_len = sam_ref_compact_ref (r);
            r->is_compacted = true;
        }
        else {
            r->uncompressed_len = r->last_pos - r->first_pos + 1;                
            r->is_compacted = false;
        }

        ((VBlockSAMP)vb)->range = r;
        vb->ready_to_dispatch = true;
        next_range_i++;
    }
}

void sam_ref_compress_ref (void)
{
    if (!buf_is_allocated (&ranges)) return;

    next_range_i=0; // can be initialized multiple times if we're compressing multiple SAMs - but not concurrently

    dispatcher_fan_out_task ("Internal reference sequence", sam_ref_prepare_range_for_compress, sam_ref_compress_one_range, sam_ref_output_one_range);
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
