// ------------------------------------------------------------------
//   refhash.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "mutex.h"
#include "codec.h" // must be included before reference.h
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
#include "bits.h"
#include "profiler.h"
#include "threads.h"
#include "segconf.h"
#include "strings.h"

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

// values used in make-reference. when loading a reference file, we get it from the file
// note: the code supports modifying MAKE_REF_BASE_LAYER_BITS and MAKE_REF_NUM_LAYERS without further code changes.
// this only affects make-reference. These values were chosen to balance RAM and almost optimimum "near-perfect short read matches"
// (47-49% of reads in our fastq benchmarks are 95%+ matched to the reference), with memory use of refhash - 1.875 GB
#define MAKE_REF_BASE_LAYER_BITS 28   // number of bits of layer_i=0. each subsequent layer has 1 bit less 
#define MAKE_REF_NUM_LAYERS 4
static uint32_t make_ref_vb_size = 0; // max bytes in a refhash vb

// each rehash_buf entry is a Buffer containing a hash layer containing 256M gpos values (4B each) = 1GB
unsigned num_layers=0;
bool bits_per_hash_is_odd=0; // true bits_per_hash is odd
uint32_t bits_per_hash=0;    // = layer_bits[0]
uint32_t nukes_per_hash=0;   // = layer_bits[0] / 2
uint32_t layer_bits[64];     // number of bits in each layer - layer_bits[0] is the base (widest) layer
uint32_t layer_size[64];     // size (in bytes) of each layer - layer_size[0] is the base (biggest) layer
uint32_t layer_bitmask[64];  // 1s in the layer_bits[] LSbs
static Buffer refhash_buf = {}; // One buffer that includes all layers
uint32_t **refhashs = NULL;  // array of pointers to into refhash_buf.data - beginning of each layer 

static Section sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 

// used for parallelizing read / write of the refhash
static uint32_t next_task_layer = 0;
static uint32_t next_task_start_within_layer = 0;

// lookup table for base complement
const char complement[256] =  { ['A']='T', ['C']='G', ['G']='C', ['T']='A',  // complement A,C,G,T, others are 4
                                [0 ...'@']=4, ['B']=4, ['D'...'F']=4, ['U'...255]=0 };

// cache stuff
static ThreadId refhash_cache_creation_thread_id;
static VBlockP cache_create_vb = NULL;

// ------------------------------------------------------
// stuff related to creating the refhash
// ------------------------------------------------------

static inline uint32_t refhash_get_word (const Range *r, int64_t base_i)
{
    // num_bits_this_range will be:
    // * 0 if HOOK is in this range, but all bits in next range
    // * larger than 0, smaller than BITS_PER_HASH if HOOK is in this range, 
    //   and some bits in this range and some in next
    // * BITS_PER_HASH if HOOK and all bits are in this range
    PosType num_bits_this_range = MIN_(bits_per_hash, (PosType)r->ref.nbits - base_i*2);
    uint32_t refhash_word = 0;

    // if the are any of the 29 bits in this range, use them
    if (num_bits_this_range > 0)
        refhash_word = bits_get_wordn (&r->ref, base_i * 2, num_bits_this_range);

    // if there are any bits in the next range, they are the more significant bits in the word we are forming
    if (num_bits_this_range < bits_per_hash)
        // note: refhash_calc_one_range guarantees us that if the word overflows to the next range, then there is a valid next range with sufficient bits
        refhash_word |= bits_get_wordn (&(r+1)->ref, 0, bits_per_hash - num_bits_this_range) << num_bits_this_range;

    return refhash_word;
}

// get base - might be in this range or next range
#define GET_BASE(idx) (idx           < num_bases        ? ref_base_by_idx (r, idx)                : \
                       idx-num_bases < ref_size(next_r) ? ref_base_by_idx (next_r, idx-num_bases) : \
                       /* else */                         'X' /* no more bases */)

// make-reference: generate the refhash data for one range of the reference. called by ref_compress_one_range (compute thread)
// towards the end of the range, we might have hash values that start in this range and end in the next range.
void refhash_calc_one_range (const Range *r, const Range *next_r /* NULL if r is the last range */)
{
    PosType this_range_size = ref_size (r);
    PosType next_range_size = ref_size (next_r);
    
    ASSERT (this_range_size * 2 == r->ref.nbits, 
            "mismatch between this_range_size=%"PRId64" (x2 = %"PRId64") and r->ref.nbits=%"PRIu64". Expecting the latter to be exactly double the former. chrom=%s r->first_pos=%"PRId64" r->last_pos=%"PRId64" r->range_id=%u", 
            this_range_size, this_range_size*2, r->ref.nbits, Bc (ZCTX(0)->dict, B(CtxNode, ZCTX(0)->nodes, r->chrom)->char_index), 
            r->first_pos, r->last_pos, r->range_id);
            
    // number of bases - considering the availability of bases in the next range, as we will overflow to it at the
    // end of this one (note: we only look at one next range - even if it is very short, we will not overflow to the next one after)
    PosType num_bases = this_range_size - (nukes_per_hash - MIN_(next_range_size, nukes_per_hash) ); // take up to NUKES_PER_HASH bases from the next range, if available

    for (PosType base_i=0; base_i < num_bases; base_i++)

        // take only the final hook in a polymer string of hooks (i.e. the last G in a e.g. GGGGG)
        if (GET_BASE(base_i) == HOOK && GET_BASE(base_i+1) != HOOK) {
            
            uint32_t refhash_word = refhash_get_word (r, base_i+1); // +1 so not to include the hook

            // since our refhash entries are 32 bit, we cannot use the reference data beyond the first 4Gbp for creating the refhash
            // TO DO: make the hash entries 40bit (or 64 bit?) if genome size > 4Gbp (bug 150)
            if (r->gpos + base_i > MAX_ALIGNER_GPOS) {
                WARN_ONCE ("FYI: %s contains more than %s bases. When compressing a FASTQ or unaligned (i.e. missing RNAME, POS) SAM/BAM file using the reference being generated, only the first %s bases of the reference will be used (no such limitation when compressing other file types). This might affect the compression ratio.", 
                            txt_name, str_int_commas (MAX_ALIGNER_GPOS).s, str_int_commas (MAX_ALIGNER_GPOS).s);
                return;
            }

            // we enter our own gpos to the hash table at refhash_word - we place it in the first empty layer.
            // if all the layers are full - then with probability 50%, we replace one of the layer in random.
            // (unlikely edge case: an entry with gpos=0 will always be replaced)
            bool set=false;
            for (unsigned layer_i=0; layer_i < num_layers; layer_i++) {
                
                uint32_t idx = refhash_word & layer_bitmask[layer_i];

                if (refhashs[layer_i][idx] == NO_GPOS) {
                    refhashs[layer_i][idx] = BGEN32 (r->gpos + base_i);
                    set=true;
                    break;
                }
            }

            // if all layers are already set, and we are on the lucky side of a 25% chance we 
            // overwrite one of the layers
            if (!set && ((rand() & 3) == 0)) { 
                uint32_t layer_i = rand() % num_layers;
                uint32_t idx = refhash_word & layer_bitmask[layer_i];
                refhashs[layer_i][idx] = BGEN32 (r->gpos + base_i); // replace one layer in random
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
    vb->dispatch = READY_TO_COMPUTE;

    // incremenet parameters for the next vb
    next_task_start_within_layer += make_ref_vb_size;

    if (next_task_start_within_layer >= layer_size[next_task_layer]) {
        next_task_layer++;
        next_task_start_within_layer = 0;
    }
}

// part of --make-reference - compute thread for compressing part of the hash
static void refhash_compress_one_vb (VBlockP vb)
{
    uint32_t uncompressed_size = MIN_(make_ref_vb_size, layer_size[vb->refhash_layer] - vb->refhash_start_in_layer);
    const uint32_t *hash_data = &refhashs[vb->refhash_layer][vb->refhash_start_in_layer / sizeof (uint32_t)];

    // calculate density to decide on compression codec
    uint32_t num_zeros=0;
    for (uint32_t i=0 ; i < uncompressed_size / sizeof (uint32_t); i++)
        if (!hash_data[i]) num_zeros++;

    SectionHeaderRefHash header = { .h.section_type          = SEC_REF_HASH, 
                                    .h.codec                 = CODEC_RANS32, // Much!! faster than LZMA (compress and uncompress), 8% worse compression of human refs, MUCH better on small refs
                                    .h.data_uncompressed_len = BGEN32 (uncompressed_size),
                                    .h.vblock_i              = BGEN32 (vb->vblock_i),
                                    .h.magic                 = BGEN32 (GENOZIP_MAGIC),
                                    .h.compressed_offset     = BGEN32 (sizeof(header)),
                                    .num_layers              = (uint8_t)num_layers,
                                    .layer_i                 = (uint8_t)vb->refhash_layer,
                                    .layer_bits              = (uint8_t)layer_bits[vb->refhash_layer],
                                    .start_in_layer          = BGEN32 (vb->refhash_start_in_layer)     };

    comp_compress (vb, NULL, &vb->z_data, (SectionHeaderP)&header, (char*)hash_data, NO_CALLBACK, "SEC_REF_HASH");

    if (flag.show_ref_hash) 
        iprintf ("vb_i=%u Compressing SEC_REF_HASH num_layers=%u layer_i=%u layer_bits=%u start=%u size=%u bytes size_of_disk=%u bytes\n", 
                 vb->vblock_i, header.num_layers, header.layer_i, header.layer_bits, vb->refhash_start_in_layer, uncompressed_size, BGEN32 (header.h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderRefHash));

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

// ZIP-FASTA-make-reference: called by main thread in zip_write_global_area
void refhash_compress_refhash (void)
{
    next_task_layer = 0;
    next_task_start_within_layer = 0;

    dispatcher_fan_out_task ("compress_refhash", NULL, 0, "Writing hash table (this can take several minutes)...", 
                             false, false, false, 0, 100,
                             refhash_prepare_for_compress, 
                             refhash_compress_one_vb, 
                             zfile_output_processed_vb);
}

// -----------------------------------
// stuff related to refhash cache file
// -----------------------------------

static inline rom refhash_get_cache_fn (void)
{
    static char *cache_fn = NULL;

    if (!cache_fn) {
        cache_fn = MALLOC (strlen (z_name) + 20);
        sprintf (cache_fn, "%s.hcache", z_name);
    }

    return cache_fn;
}

void refhash_remove_cache (void)
{
    file_remove (refhash_get_cache_fn(), true);
}

// thread entry for creating refhash cache
static void refhash_create_cache (VBlockP unused)
{
    buf_dump_to_file (refhash_get_cache_fn(), &refhash_buf, 1, true, false, false, false);
}

static void refhash_create_cache_in_background (void)
{
    if (flag.regions) return; // can't create cache as reference isn't fully loaded

    cache_create_vb = vb_initialize_nonpool_vb (VB_ID_HCACHE_CREATE, DT_NONE, "refhash_create_cache_in_background");
    cache_create_vb->compute_task = "create_refhash_cache";

    refhash_get_cache_fn(); // generate name before we close z_file
    refhash_cache_creation_thread_id = threads_create (refhash_create_cache, cache_create_vb);
}

void refhash_create_cache_join (bool free_mem)
{
    if (!cache_create_vb) return;

    threads_join (&refhash_cache_creation_thread_id);   

    if (free_mem) vb_destroy_vb (&cache_create_vb);
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

    if (flag.show_ref_hash) // before the sanity checks
        iprintf ("vb_i=%u Uncompressing SEC_REF_HASH num_layers=%u layer_i=%u layer_bits=%u start=%u size=%u bytes size_of_disk=%u bytes\n", 
                 vb->vblock_i, header->num_layers, layer_i, header->layer_bits, start, size, BGEN32 (header->h.data_compressed_len) + (uint32_t)sizeof (SectionHeaderRefHash));

    // sanity checks
    ASSERT (layer_i < num_layers, "expecting header->layer_i=%u < num_layers=%u", layer_i, num_layers);

    ASSERT (header->layer_bits == layer_bits[layer_i], "expecting header->layer_bits=%u to be %u", header->layer_bits, layer_bits[layer_i]);

    ASSERT (start + size <= layer_size[layer_i], "expecting start=%u + size=%u <= layer_size=%u", start, size, layer_size[layer_i]);

    // a hack for uncompressing to a location withing the buffer - while multiple threads are uncompressing into 
    // non-overlappying regions in the same buffer in parallel
    Buffer copy = refhash_buf;
    copy.data = (char *)refhashs[layer_i] + start; // refhashs[layer_i] points to the start of the layer within refhash_buf.data
    zfile_uncompress_section (vb, header, &copy, NULL, 0, SEC_REF_HASH);

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}

static void refhash_read_one_vb (VBlockP vb)
{
    buf_alloc (vb, &vb->z_section_headers, 0, 1, int32_t, 0, "z_section_headers"); // room for 1 section header

    if (!sections_next_sec (&sl_ent, SEC_REF_HASH))
        return; // no more refhash sections

    int32_t section_offset = zfile_read_section (z_file, vb, sl_ent->vblock_i, &vb->z_data, "z_data", sl_ent->st, sl_ent);

    if (((SectionHeaderRefHash *)vb->z_data.data)->layer_i >= num_layers)
        return; // don't read the high layers if beyond the requested num_layers

    BNXT (int32_t, vb->z_section_headers) = section_offset;

    vb->dispatch = READY_TO_COMPUTE;
}

void refhash_load_standalone (void)
{
    flag.reading_reference = gref; // tell file.c and fasta.c that this is a reference

    TEMP_VALUE (command, PIZ);
    CLEAR_FLAG (test);

    z_file = file_open (ref_get_filename (gref), READ, Z_FILE, DT_FASTA);    
    z_file->basename = file_basename (ref_get_filename (gref), false, "(reference)", NULL, 0);

    zfile_read_genozip_header (0);

    refhash_initialize (NULL);

    file_close (&z_file, false, false);
    file_close (&txt_file, false, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtheader_piz_read_and_reconstruct.
    
    RESTORE_FLAG (test);
    RESTORE_VALUE (command);

    flag.reading_reference = NULL;
}

// ----------------------------
// general stuff
// ----------------------------

static void refhash_initialize_refhashs_array (void)
{
    refhashs = CALLOC (num_layers * sizeof (uint32_t *)); // array of pointers

    // set layer pointers
    uint64_t offset=0;
    for (unsigned layer_i=0; layer_i < num_layers; layer_i++) {
        refhashs[layer_i] = (uint32_t*)B8 (refhash_buf, offset);
        offset += layer_size[layer_i];
    }
}

// called by the main thread - piz_read_global_area when reading the reference file, ahead of compressing a fasta or fastq file. 
// returns true if mapped cache
void refhash_initialize (bool *dispatcher_invoked)
{
    uint32_t base_layer_bits;
    
    if (dispatcher_invoked) *dispatcher_invoked = false; // initialize

    if (buf_is_alloc (&refhash_buf)) return; // already loaded from a previous file

    // case 1: called from ref_make_ref_init - initialize for making a reference file
    if (flag.make_reference) {
        #define VBLOCK_MEMORY_REFHASH  (16 << 20) // VB memory with --make-reference - refhash data (overridable with --vblock)
        num_layers       = MAKE_REF_NUM_LAYERS;
        base_layer_bits  = MAKE_REF_BASE_LAYER_BITS;
        // we use the default vb size (16MB) not the reduced make-ref size (1MB), unless user overrides with --vblock
        make_ref_vb_size = flag.vblock ? segconf.vb_size : VBLOCK_MEMORY_REFHASH;
    }

    // case 2: piz_read_global_area from piz_read_global_area -> refhash_load - initialize when reading an external reference for ZIP of fasta or fastq
    else 
        sections_get_refhash_details (num_layers ? NULL : &num_layers, &base_layer_bits);

    uint64_t refhash_size = 0;
    for (unsigned layer_i=0; layer_i < num_layers; layer_i++) {
        layer_bits[layer_i]    = base_layer_bits - layer_i;
        layer_bitmask[layer_i] = bitmask32 (layer_bits[layer_i]);
        layer_size[layer_i]    = ((1 << layer_bits[layer_i]) * sizeof (uint32_t));
        refhash_size          += layer_size[layer_i];
    }

    bits_per_hash  = layer_bits[0];
    bits_per_hash_is_odd = bits_per_hash % 2;
    nukes_per_hash = (1 + bits_per_hash) / 2; // round up

    // if not making reference - we try to load - first from cache, then from reference file
    if (!flag.make_reference) {

        // attempt to mmap the cache, but if it doesn't exist read from the reference and create the cache
        bool mapped_cache = buf_mmap (evb, &refhash_buf, refhash_get_cache_fn(), true, "refhash_buf");
        if (!mapped_cache) { 
            // allocate memory - base layer size is 1GB, and every layer is half the size of its predecessor, so total less than 2GB
            buf_alloc (evb, &refhash_buf, 0, refhash_size, uint8_t, 1, "refhash_buf"); 
            refhash_initialize_refhashs_array();

            sl_ent = NULL; // NULL -> first call to this sections_get_next_ref_range() will reset cursor 
            dispatcher_fan_out_task ("load_refhash", ref_get_filename (gref),
                                     0, "Reading and caching reference hash table...",
                                     true, flag.test, false, 0, 100,
                                     refhash_read_one_vb, 
                                     refhash_uncompress_one_vb, 
                                     NO_CALLBACK);
 
            refhash_create_cache_in_background();
            
            if (dispatcher_invoked) *dispatcher_invoked = true;
        }
    }

    else { // make_reference
        // set all entries to NO_GPOS. note: no need to set in ZIP, as we will be reading the data from the refernce file
        // NOTE: setting NO_GPOS to 0xff rather than 0x00 causes make-ref to take ~8 min on my PC instead of < 1 min
        // due to a different LZMA internal mode when compressing the hash. However, the resulting file is MUCH smaller,
        // and loading of refhash during zip is MUCH faster
        buf_alloc (evb, &refhash_buf, 0, refhash_size, uint8_t, 1, "refhash_buf"); 
        buf_set (&refhash_buf, 0xff); 
    }

    // if we haven't prepared refhashs yet, do it now.
    if (!refhashs) refhash_initialize_refhashs_array();
}

void refhash_destroy (bool destroy_only_if_not_mmap)
{
    if (refhash_buf.type == BUF_UNALLOCATED || 
        (refhash_buf.type == BUF_MMAP_RO && destroy_only_if_not_mmap)) return;

    refhash_create_cache_join (true); // wait for cache writing, if we're writing

    buf_destroy (refhash_buf);
    FREE (refhashs);

    flag.aligner_available = false;
}

void refhash_verify_before_exit (void)
{
    // verify read-only mmap'ed reference integrity before existing process
    buf_destroy (refhash_buf); // also verifies loaded mmap
}

bool refhash_has_refhash (void)
{
    return buf_is_alloc (&refhash_buf);
}
