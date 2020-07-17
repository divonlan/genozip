// ------------------------------------------------------------------
//   refhash.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "arch.h"
#include "reference.h"
#include "sections.h"
#include "buffer.h"
#include "vblock.h"
#include "file.h"
#include "endianness.h"
#include "dispatcher.h"
#include "zfile.h"
#include "refhash.h"
#include "compressor.h"
#include "bit_array.h"

// ref_hash logic:
// we use the 28 bits (14 nucleotides) following a "G" hook, as the hash value, to index into a hash table with multiple
// layers - each layer being of size of half of the size of the previous layer.
// the full hash value is used to index into the first hash table layer, the 27 LSb of the hash value into the 2nd layer etc.
//
// The size of the first hash layer is 256M entries (28 bits) * sizeof(uint32_t) = 1GB, and all layers < 2 GB. 
//
// With the simplified assumption that G appears 1/4 of all bases, then a human genome will have ~ 3B / 4 =~ 750M hash values 
// competing to enter these 256M (one layer) - 512M slots (28 layers).
//
// The expected number of times any particular 29 bit ("G" + 28) string will appear in a 3B genome is 3B / 750M = 4.
// On a typical 150 bases short read, we expected 150/4 =~ 37.5 'G' hooks - hopefully one of them will point to the correct
// reference gpos
//
#define HOOK 'G'

// values used in make-reference. when loading a reference file, we get it from the file
// note: the code supports modifying MAKE_REF_BASE_LAYER_BITS and MAKE_REF_NUM_LAYERS without further code changes.
// this only affects make-reference. These values were chosen to balance RAM and almost optimimum "near-perfect short read matches"
// (47-49% of reads in our fastq benchmarks are 95%+ matched to the reference), with memory use of refhash - 1.875 GB
#define MAKE_REF_BASE_LAYER_BITS 28 // number of bits of layer_i=0. each subsequent layer has 1 bit less 
#define MAKE_REF_NUM_LAYERS 4
static uint32_t make_ref_vb_size = 0; // max bytes in a refhash vb

// each rehash_buf entry is a Buffer containing a hash layer containing 256M gpos values (4B each) = 1GB
static unsigned num_layers=0;
static uint32_t bits_per_hash=0;  // = layer_bits[0]
static uint32_t nukes_per_hash=0; // = layer_bits[0] / 2
static uint32_t layer_bits[64]; // number of bits in each layer - layer_bits[0] is the base (widest) layer
static uint32_t layer_size[64]; // size (in bytes) of each layer - layer_size[0] is the base (biggest) layer
static uint32_t layer_bitmask[64]; // 1s in the layer_bits[] LSbs
static Buffer *refhash_bufs = NULL; // array of buffers, one for each layer
static uint32_t **refhashs = NULL; // array of pointers to the beginning of each layer

static SectionListEntry *sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
static uint32_t ref_hash_cursor = 0;


// used for parallelizing read / write of the refhash
static uint32_t next_task_layer = 0;
static uint32_t next_task_start_within_layer = 0;

static inline uint32_t refhash_get_word (const Range *r, int64_t idx)
{
    // num_bits_this_range will be:
    // * 0 if HOOK is in this range, but all bits in next range
    // * larger than 0, smaller than BITS_PER_HASH if HOOK is in this range, 
    //   and some bits in this range and some in next
    // * BITS_PER_HASH if HOOK and all bits are in this range
    int64_t num_bits_this_range = MIN (bits_per_hash, (int64_t)r->ref.num_of_bits - idx*2);

    uint32_t refhash_word = 0;

    // if the are any of the 29 bits in this range, use them
    if (num_bits_this_range > 0)
        refhash_word = bit_array_get_wordn (&r->ref, idx * 2, num_bits_this_range);

    // if there are any bits in the next range, they are the more significant bits in the word we are forming
    if (num_bits_this_range < bits_per_hash)
        // note: refhash_calc_one_range guarantees us that if the word overflows to the next range, then there is a valid next range with sufficient bits
        refhash_word |= bit_array_get_wordn (&(r+1)->ref, 0, bits_per_hash - num_bits_this_range) << num_bits_this_range;

    return refhash_word;
}

// get base - might be in this range or next range
#define GET_BASE(idx) (idx           < num_bases        ? ref_get_nucleotide (r, idx)                : \
                       idx-num_bases < ref_size(next_r) ? ref_get_nucleotide (next_r, idx-num_bases) : \
                       /* else */                         'X' /* no more bases */)

// make-reference: generate the refhash data for one range of the reference. called by ref_compress_one_range (compute thread)
// towards the end of the range, we might have hash values that start in this range and end in the next range.
void refhash_calc_one_range (const Range *r, const Range *next_r /* NULL if r is the last range */)
{
    int64_t this_range_size = ref_size (r);
    int64_t next_range_size = ref_size (next_r);
    
    // number of bases - considering the availability of bases in the next range, as we will overflow to it at the
    // end of this one (note: we only look at one next range - even if it is very short, we will not overflow to the next one after)
    int64_t num_bases = this_range_size - (nukes_per_hash - MIN (next_range_size, nukes_per_hash) ); // take up to NUKES_PER_HASH bases from the next range, if available

    for (int64_t i=0; i < num_bases; i++)

        // take only the final hook in a polymer string of hooks (i.e. the last G in a e.g. GGGGG)
        if (GET_BASE(i) == HOOK && GET_BASE(i+1) != HOOK) {
            
            uint32_t refhash_word = refhash_get_word (r, i+1); // +1 so not to include the hook

            // since our refhash entries are 32 bit, we cannot use the reference data beyond the first 4Gbp for creating the refhash
            // TO DO: make the hash entries 40bit (or 64 bit?) if genome size > 4Gbp (bug 150)
            if (r->gpos + i > ((uint64_t)1 << 32)) {
                static bool warning_given = false;

                ASSERTW (warning_given, "Warning: %s contains more than 2^32 (~4 billion) nucleaotides. When compressing a FASTQ or FASTA file using the reference being generated, only the first 2^32 nucleotides of the reference will be used (no such limitation when compressing other file types)", txt_name);
                warning_given = true; // display this warning only once
                return;
            }

            // we enter our own gpos to the hash table at refhash_word - we place it in the first empty layer.
            // if all the layers are full - then with probability 50%, we replace one of the layer in random.
            // (unlikely edge case: an entry with gpos=0 will always be replaced)
            bool set=false;
            for (unsigned layer_i=0; layer_i < num_layers; layer_i++) {
                
                uint32_t idx = refhash_word & layer_bitmask[layer_i];

                if (refhashs[layer_i][idx] == 0) {
                    refhashs[layer_i][idx] = BGEN32 (r->gpos + i);
                    set=true;
                    break;
                }
            }

            if (!set && ((rand() & 3) == 0)) { // if all layers are already set, and we are on the lucky side of a 25% chance
                uint32_t layer_i = rand() % num_layers;
                uint32_t idx = refhash_word & layer_bitmask[layer_i];
                refhashs[layer_i][idx] = BGEN32 (r->gpos + i); // replace one layer in random
            }
        }
}

// compress the reference - one section at the time, using Dispatcher to do them in parallel (make_ref version)
static void refhash_prepare_for_compress (VBlockP vb)
{
    if (next_task_layer == num_layers) return; // we're done

    // tell this vb what to do
    vb->refhash_layer = next_task_layer;
    vb->refhash_start_in_layer = next_task_start_within_layer;
    vb->ready_to_dispatch = true;

    // incremenet parameters for the next vb
    next_task_start_within_layer += make_ref_vb_size;

    if (next_task_start_within_layer >= layer_size[next_task_layer]) {
        next_task_layer++;
        next_task_start_within_layer = 0;
    }
}

static void refhash_compress_one_vb (VBlockP vb)
{
    uint32_t uncompressed_size = MIN (make_ref_vb_size, layer_size[vb->refhash_layer] - vb->refhash_start_in_layer);

    // tradeoff between LZMA and BZ2 in this case:
    // BZ2 pros : BZ2 is about 8x faster than LZMA in generating the reference (56 sec vs 7:42 min on my PC),
    // LZMA pros: LZMA it saves 4-6 seconds (on my PC) in loading the generated reference for compressing fasta/fastq due to faster decompression. The reference file is ~1.5% smaller.
    SectionHeaderRefHash header = { .h.section_type          = SEC_REF_HASH, 
                                    .h.sec_compression_alg   = COMP_LZMA, 
                                    .h.data_uncompressed_len = BGEN32 (uncompressed_size),
                                    .h.vblock_i              = BGEN32 (vb->vblock_i),
                                    .h.magic                 = BGEN32 (GENOZIP_MAGIC),
                                    .h.compressed_offset     = BGEN32 (sizeof(header)),
                                    .num_layers              = (uint8_t)num_layers,
                                    .layer_i                 = (uint8_t)vb->refhash_layer,
                                    .layer_bits              = (uint8_t)layer_bits[vb->refhash_layer],
                                    .start_in_layer          = BGEN32 (vb->refhash_start_in_layer)     };

    vb->z_data.name  = "z_data"; // comp_compress requires that these are pre-set
    vb->z_data.param = vb->vblock_i;
    
    comp_compress (vb, &vb->z_data, false, (SectionHeaderP)&header, ENT (char, refhash_bufs[vb->refhash_layer], vb->refhash_start_in_layer), NULL);

    if (flag_show_ref_hash) 
        fprintf (stderr, "vb_i=%u Compressing SEC_REF_HASH num_layers=%u layer_i=%u layer_bits=%u start=%u size=%u bytes size_of_disk=%u bytes\n", 
                 vb->vblock_i, header.num_layers, header.layer_i, header.layer_bits, vb->refhash_start_in_layer, uncompressed_size, BGEN32 (header.h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderRefHash));

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

// ZIP-FASTA-make-reference: called by I/O thread in zip_write_global_area
void refhash_compress_refhash (void)
{
    next_task_layer = 0;
    next_task_start_within_layer = 0;

    dispatcher_fan_out_task (NULL, "Writing hash table used for compressing fastq and fasta files...", false, 
                             refhash_prepare_for_compress, 
                             refhash_compress_one_vb, 
                             ref_output_vb);
}

// ----------------------------------------------------------------------------------------------
// stuff related to ZIPping fasta and fastq using an external reference that includes the refhash
// ----------------------------------------------------------------------------------------------

// entry point of compute thread of refhash decompression
static void refhash_uncompress_one_vb (VBlockP vb)
{
    SectionHeaderRefHash *header = (SectionHeaderRefHash *)vb->z_data.data;
    uint32_t start = BGEN32 (header->start_in_layer);
    uint32_t size  = BGEN32 (header->h.data_uncompressed_len);
    uint32_t layer_i = header->layer_i;

    if (flag_show_ref_hash) // before the sanity checks
        fprintf (stderr, "vb_i=%u Uncompressing SEC_REF_HASH num_layers=%u layer_i=%u layer_bits=%u start=%u size=%u bytes size_of_disk=%u bytes\n", 
                 vb->vblock_i, header->num_layers, layer_i, header->layer_bits, start, size, BGEN32 (header->h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderRefHash));

    // sanity checks
    ASSERT (layer_i < num_layers, "Error in refhash_uncompress_one_vb: expecting header->layer_i=%u < num_layers=%u", layer_i, num_layers);

    ASSERT (header->layer_bits == layer_bits[layer_i], "Error in refhash_uncompress_one_vb: expecting header->layer_bits=%u to be %u",
            header->layer_bits, layer_bits[layer_i]);

    ASSERT (start + size <= layer_size[layer_i], "Error in refhash_uncompress_one_vb: expecting start=%u + size=%u <= layer_size=%u",
            start, size, layer_size[layer_i]);

    zfile_uncompress_section (vb, header, &refhashs[layer_i][start / sizeof(uint32_t)], NULL, SEC_REF_HASH);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

static void refhash_read_one_vb (VBlockP vb)
{
    buf_alloc (vb, &vb->z_section_headers, 1 * sizeof(int32_t), 0, "z_section_headers", 0); // room for 1 section header

    if (!sections_get_next_section_of_type (&sl_ent, &ref_hash_cursor, SEC_REF_HASH, SEC_NONE))
        return; // no more refhash sections

    int32_t section_offset = zfile_read_section (z_file, vb, sl_ent->vblock_i, NO_SB_I, &vb->z_data, 
                                                 "z_data", sizeof(SectionHeaderRefHash), sl_ent->section_type, sl_ent);

    if (((SectionHeaderRefHash *)vb->z_data.data)->layer_i >= num_layers)
        return; // don't read the high layers if beyond the requested num_layers

    ASSERT (section_offset != EOF, "Error in ref_read_one_range: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);

    NEXTENT (int32_t, vb->z_section_headers) = section_offset;

    vb->ready_to_dispatch = true;
}

// called by the I/O thread - piz_read_global_area when reading the reference file, ahead of compressing a fasta or fastq file. 
void refhash_load(void)
{
    refhash_initialize();

    sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
    ref_hash_cursor = 0;

    dispatcher_fan_out_task (ref_filename,
                             "Reading reference hash table...", flag_test, 
                             refhash_read_one_vb, 
                             refhash_uncompress_one_vb, 
                             NULL);

    buf_test_overflows_all_vbs ("refhash_load");
}

// strict encoding of A,C,G,T - everything else in non-encodable (a 4 here)
static const uint32_t nuke_encode[256] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 4
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 16
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 32
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 48
                                           4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,   // 64  A(65)->0 C(67)->1 G(71)->2
                                           4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 84  T(84)->3
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 96  
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 112 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   // 128
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

// calculates a refhash word from 14 nucleotides following a 'G' (only last G in a sequenece of GGGG...)
static inline bool refhash_get_word_from_seq (const char *seq, uint32_t *refhash_word)
{   
    *refhash_word = 0;

    for (unsigned i=0; i < nukes_per_hash; i++) {   
        uint32_t base = nuke_encode[(uint8_t)seq[i]];
        if (base == 4) return false; // not a A,C,G,T

        *refhash_word = *refhash_word + (base << (i*2));
    }

    *refhash_word &= layer_bitmask[0]; // remove MSb in case of an odd number of base bits

    return true;
}

int64_t refhash_get_final_match_len (const char *seq, int64_t gpos, int64_t fwd_len)
{
    int64_t count=0;

    for (int64_t i=0; i < fwd_len; i++) 
        if (seq[i] == ACTG_DECODE (&genome->ref, gpos + i)) count++;

    return count;
}

static inline int64_t refhash_get_match_len (const char *seq, int64_t gpos, int64_t fwd_len, int64_t bwd_len,
                                             int64_t *fwd_match, int64_t *bwd_match, int64_t *snps)
{
    *snps=0;
    int64_t i;

    // search forwards
    for (i=0; i < fwd_len; i++) 
        if (seq[i] != ACTG_DECODE (&genome->ref, gpos + i)) {
            if (i < fwd_len-1 && seq[i+1] == ACTG_DECODE (&genome->ref, gpos + i+1))
                (*snps)++; // tolerate a SNP - this position is different, but next position in the same again
            else
                break;            
        }
    *fwd_match = i;

    // search backwards
    for (i=-1; i >= -bwd_len; i--)
        if (seq[i] != ACTG_DECODE (&genome->ref, gpos + i)) {
            if (i >= -bwd_len+1 && seq[i-1] == ACTG_DECODE (&genome->ref, gpos + i-1))
                (*snps)++; // tolerate a SNP - this position is different, but next position in the same again
            else
                break;            
        }
    *bwd_match = -1 - i;

    return *fwd_match + *bwd_match; 
}

// returns gpos aligned with seq, containing the longest match to the reference, or REFHASH_NOMATCH if no match was found
int64_t refhash_best_match (const char *seq, int64_t seq_len)
{
    int64_t longest_len = 0, num_matches=0, total_len=0, hooks=0, fwd_match, bwd_match, snps, final_match_len,
            best_fwd_match=0, best_bwd_match=0, best_snps=0, best_gpos=0, best_i=0;

    bool almost_perfect = false;
    for (int64_t i=0; i < seq_len - nukes_per_hash && !almost_perfect; i++) {
        
        if (seq[i] != HOOK || seq[i+1] == HOOK) continue; // take the G - if there is a polymer GGGG... take the last one
        hooks++;

        uint32_t refhash_word;
        if (!refhash_get_word_from_seq (&seq[i+1], &refhash_word)) continue; 

        for (unsigned layer_i=0; layer_i < num_layers; layer_i++) {

            int64_t gpos = (int64_t)BGEN32 (refhashs[layer_i][refhash_word & layer_bitmask[layer_i]]); // position of the start of the G... sequence in the genome
            
            // case: hash table is empty for this word - and the remaining layers do not have an entry that originated from refhash_word.
            // note that they might have an entry from another word that shares LSb with refhash_word, but that by definition doesnt much our seq)
            if (!gpos) break; 

            int64_t match_len = refhash_get_match_len (&seq[i], gpos, MIN (seq_len-i, genome_size-gpos), MIN(i, gpos), &fwd_match, &bwd_match, &snps);
            if (match_len < 20) continue;

            if (match_len > longest_len) {
                longest_len    = match_len;
                best_gpos      = gpos;
                best_i         = i;
                best_fwd_match = fwd_match;
                best_bwd_match = bwd_match;
                best_snps      = snps;
            }

            total_len += match_len;
            num_matches++;

            // case: we have over 95% match, no need to continue searching
            if (match_len >= (seq_len * 95) / 100) {
                almost_perfect=true;
                break; 
            }
        }
    }

    return longest_len ? (best_gpos - best_i) : REFHASH_NOMATCH; // can be negative if seq is aligned to before the start of the reference
/*
    if (longest_len > 0 && longest_len < seq_len) {
        int64_t start_gpos = MAX (best_gpos - best_i, 0);
        int64_t room_fwd = MIN (seq_len, genome_size-start_gpos);
        final_match_len = refhash_get_final_match_len (seq, start_gpos, room_fwd); 
    }
    else 
        final_match_len = longest_len;

    printf ("hooks=%u num_matches=%u longest_match_len=%u(%u%%) fwd=%u bwd=%u snps=%u final_match_len=%u(%u%%)\n", 
            (uint32_t)hooks, (uint32_t)num_matches, (uint32_t)longest_len, (uint32_t)((longest_len*100)/seq_len), 
            (uint32_t)best_fwd_match, (uint32_t)best_bwd_match, (uint32_t)best_snps,
            (uint32_t)final_match_len, (uint32_t)((final_match_len*100)/seq_len));
*/
    return almost_perfect;
}

// ----------------------------
// general stuff
// ----------------------------

void refhash_initialize (void)
{
    uint32_t base_layer_bits;

    // case 1: called from ref_make_ref_init - initialize for making a reference file
    if (flag_make_reference) {
        num_layers      = MAKE_REF_NUM_LAYERS;
        base_layer_bits = MAKE_REF_BASE_LAYER_BITS;
        // we use the default vb size (16MB) not the reduced make-ref size (1MB), unless user overrides with --vblock
        make_ref_vb_size = flag_vblock ? global_max_memory_per_vb : (atoi (TXT_DATA_PER_VB_DEFAULT) << 20);
    }

    // case 2: piz_read_global_area from piz_read_global_area -> refhash_load - initialize when reading an external reference for ZIP of fasta or fastq
    else 
        sections_get_refhash_details (num_layers ? NULL : &num_layers, &base_layer_bits);

    for (unsigned layer_i=0; layer_i < num_layers; layer_i++) {
        layer_bits[layer_i]    = base_layer_bits - layer_i;
        layer_bitmask[layer_i] = bitmask32 (layer_bits[layer_i]);
        layer_size[layer_i]    = ((1 << layer_bits[layer_i]) * sizeof (uint32_t));
    }

    bits_per_hash  = layer_bits[0];
    nukes_per_hash = (1 + bits_per_hash) / 2; // round up

    // allocate memory
    refhash_bufs = calloc (num_layers, sizeof (Buffer)); // array of Buffer
    ASSERT0 (refhash_bufs, "Error in refhash_initialize: failed to calloc refhash_bufs");

    refhashs = calloc (num_layers, sizeof (uint32_t *)); // array of pointers
    ASSERT0 (refhashs, "Error in refhash_initialize: failed to calloc refhashs");

    // base layer size is 1GB, and every layer is half the size of its predecessor, so total less than 2GB
    for (unsigned layer_i=0; layer_i < num_layers; layer_i++) {
        buf_alloc (evb, &refhash_bufs[layer_i], layer_size[layer_i], 1, "refhash_bufs", layer_i);
        if (flag_make_reference) buf_zero (&refhash_bufs[layer_i]); // no need to zero in ZIP, as we will be reading the data from the refernce file

        refhashs[layer_i] = FIRSTENT (uint32_t, refhash_bufs[layer_i]);
    }
}

void refhash_free (void)
{
    if (refhash_bufs) {
        for (unsigned layer_i=0; layer_i < num_layers; layer_i++) 
            buf_destroy (&refhash_bufs[layer_i]);

        FREE (refhash_bufs); refhash_bufs = NULL;
        FREE (refhashs);     refhashs     = NULL;
    }
}
