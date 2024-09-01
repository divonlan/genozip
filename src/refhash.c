// ------------------------------------------------------------------
//   refhash.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "compressor.h" // codec.h included in compressor.h must be included before reference.h
#include "file.h"
#include "dispatcher.h"
#include "zfile.h"
#include "refhash.h"
#include "threads.h"

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
uint32_t num_layers=0;
bool bits_per_hash_is_odd=0;     // true bits_per_hash is odd
static uint32_t bits_per_hash=0; // = layer_bits[0]
uint32_t nukes_per_hash=0;       // = layer_bits[0] / 2
static uint32_t layer_bits[64];  // number of bits in each layer - layer_bits[0] is the base (widest) layer
static uint32_t layer_size[64];  // size (in bytes) of each layer - layer_size[0] is the base (biggest) layer
uint32_t layer_bitmask[64];      // 1s in the layer_bits[] LSbs
Buffer refhash_buf = {};         // One buffer that includes all layers
uint32_t **refhashs = NULL;      // array of pointers to into refhash_buf.data - beginning of each layer 
static Digest refhash_digest = DIGEST_NONE;

// used for parallelizing read / write of the refhash
static uint32_t next_task_layer = 0;
static uint32_t next_task_start_within_layer = 0;

static Buffer emoneg_buf = {}; // reverse compliment of genome (2-bit). only used by aligner.

// -------------------------------------------
// shared logic (make refhash + use refhash)
// -------------------------------------------

uint64_t refhash_get_refhash_size (void)
{
    Section sec = sections_last_sec (SEC_REF_HASH, HARD_FAIL);

    SectionHeaderRefHash header = zfile_read_section_header (evb, sec, SEC_REF_HASH).ref_hash;
    uint32_t base_layer_bits = header.layer_bits + header.layer_i; // layer_i=0 is the base layer, layer_i=1 has 1 bit less etc

    uint64_t refhash_size = 0;
    for (uint32_t layer_i=0; layer_i < header.num_layers; layer_i++) 
        refhash_size += ((1 << (base_layer_bits - layer_i)) * sizeof (uint32_t));

    return refhash_size;
}

static uint64_t refhash_initialize_layers (uint32_t base_layer_bits)
{
    uint64_t refhash_size = 0;
    for (uint32_t layer_i=0; layer_i < num_layers; layer_i++) {
        layer_bits[layer_i]    = base_layer_bits - layer_i;
        layer_bitmask[layer_i] = bitmask32 (layer_bits[layer_i]);
        layer_size[layer_i]    = ((1 << layer_bits[layer_i]) * sizeof (uint32_t));
        refhash_size          += layer_size[layer_i];
    }

    bits_per_hash  = layer_bits[0];
    bits_per_hash_is_odd = bits_per_hash % 2;
    nukes_per_hash = (1 + bits_per_hash) / 2; // round up

    // note: base layer size is 1GB, and every layer is half the size of its predecessor, so total less than 2GB
    if (refhash_buf.type != BUF_SHM) {
        buf_alloc_exact (evb, refhash_buf, refhash_size, uint8_t, "refhash_buf"); 
        
        if (flag.make_reference) 
            memset (refhash_buf.data, 255/*NO_GPOS*/, refhash_size);
    }

    refhashs = CALLOC (num_layers * sizeof (uint32_t *)); // array of pointers

    // set layer pointers
    uint64_t offset=0;
    for (unsigned layer_i=0; layer_i < num_layers; layer_i++) {
        refhashs[layer_i] = (uint32_t*)B8 (refhash_buf, offset);
        offset += layer_size[layer_i];
    }

    __atomic_thread_fence (__ATOMIC_RELEASE); 

    return refhash_size;
}

// update addresses of refhashs if address of refhash.data has changed
void refhash_update_layers (int64_t delta_bytes)
{
    ASSERTNOTNULL (refhashs);

    ASSERT (delta_bytes % 4 == 0, "expecting delta_bytes=%"PRIu64" to be a multiple of 4", delta_bytes);

    for (unsigned layer_i=0; layer_i < num_layers; layer_i++) 
        refhashs[layer_i] += delta_bytes / sizeof (uint32_t);

    if (flag.show_cache) iprint0 ("show-cache: updating refhash layers after address changed\n");
}

// ------------------------------------------------------
// stuff related to creating the refhash
// ------------------------------------------------------

// called from ref_make_ref_init - initialize for making a reference file
void refhash_initialize_for_make (void)
{
    num_layers = MAKE_REF_NUM_LAYERS;

    #define VBLOCK_MEMORY_REFHASH (16 MB) // VB memory with --make-reference - refhash data (overridable with --vblock)
    make_ref_vb_size = flag.vblock ? segconf.vb_size : VBLOCK_MEMORY_REFHASH;

    refhash_initialize_layers (MAKE_REF_BASE_LAYER_BITS);
}

static inline uint32_t refhash_get_word (ConstRangeP r, int64_t base_i)
{
    // num_bits_this_range will be:
    // * 0 if HOOK is in this range, but all bits in next range
    // * larger than 0, smaller than BITS_PER_HASH if HOOK is in this range, 
    //   and some bits in this range and some in next
    // * BITS_PER_HASH if HOOK and all bits are in this range
    PosType64 num_bits_this_range = MIN_(bits_per_hash, (PosType64)r->ref.nbits - base_i*2);
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
#define GET_BASE(idx) ((idx)           < num_bases        ? ref_base_by_idx (r, (idx))              : \
                       (idx)-num_bases < ref_size(next_r) ? ref_base_by_idx (next_r, idx-num_bases) : \
                       /* else */                         'X' /* no more bases */)

// "random" number 0..3, treating the reference itself as a random number generator (but mod 13, to break the uneven frequencies of nucleotides in the reference)
#define RANDOM4(idx) ((bits_get4 (&r->ref, (idx)) % 13) & 3) // a *lot* faster than rand(). note: idx, not 2*idx, so never goes beyond num_bases*2

// make-reference: generate the refhash data for one range of the reference. called by ref_compress_one_range (compute thread)
// towards the end of the range, we might have hash values that start in this range and end in the next range.
void refhash_calc_one_range (VBlockP vb, // VB of reference compression dispatcher
                             ConstRangeP r, ConstRangeP next_r /* NULL if r is the last range */)
{
    START_TIMER;

    PosType64 this_range_size = ref_size (r);
    PosType64 next_range_size = ref_size (next_r);
    
    ASSERT (this_range_size * 2 == r->ref.nbits, 
            "mismatch between this_range_size=%"PRId64" (x2 = %"PRId64") and r->ref.nbits=%"PRIu64". Expecting the latter to be exactly double the former. chrom=%s r->first_pos=%"PRId64" r->last_pos=%"PRId64" r->range_id=%u", 
            this_range_size, this_range_size*2, r->ref.nbits, Bc (ZCTX(0)->dict, B(CtxNode, ZCTX(0)->nodes, r->chrom)->char_index), 
            r->first_pos, r->last_pos, r->range_id);
            
    // number of bases - considering the availability of bases in the next range, as we will overflow to it at the
    // end of this one (note: we only look at one next range - even if it is very short, we will not overflow to the next one after)
    PosType64 num_bases = this_range_size - (nukes_per_hash - MIN_(next_range_size, nukes_per_hash) ); // take up to NUKES_PER_HASH bases from the next range, if available

    // we cannot use the reference data beyond the MAX_ALIGNER_GPOS for creating the refhash
    if (r->gpos + num_bases - 1 > MAX_ALIGNER_GPOS)
        return; // note: warning already displayed in ref_compress_ref

    decl_acgt_decode;
    for (PosType64 base_i=0; base_i < num_bases; base_i++)

        // take only the final hook in a polymer string of hooks (i.e. the last G in a e.g. GGGGG)
        if (GET_BASE(base_i) == HOOK && GET_BASE(base_i+1) != HOOK) {
            uint32_t refhash_word = refhash_get_word (r, base_i+1); // +1 so not to include the hook

            // we enter our own gpos to the hash table at refhash_word - we place it in the first empty layer.
            // if all the layers are full - then with probability 50%, we replace one of the layer in random.
            // (unlikely edge case: an entry with gpos=0 will always be replaced)
            bool set=false;
            for (unsigned layer_i=0; layer_i < num_layers; layer_i++) {    
                uint32_t idx = refhash_word & layer_bitmask[layer_i];

                if (refhashs[layer_i][idx] == (uint32_t)NO_GPOS) {
                    refhashs[layer_i][idx] = BGEN32 (r->gpos + base_i);
                    set=true;
                    break;
                }
            }

            // if all layers are already set, and we are on the lucky side of a 25% chance, we overwrite one of the layers
            if (!set && (RANDOM4(base_i+2) == 0)) { 
                uint32_t layer_i = RANDOM4(base_i+3); // choose layer to overwrite. RANDOM4 bc MAKE_REF_NUM_LAYERS=4
                uint32_t idx = refhash_word & layer_bitmask[layer_i];
                refhashs[layer_i][idx] = BGEN32 (r->gpos + base_i); // replace one layer selected randomly
            }
        }

    COPY_TIMER(refhash_calc_one_range);
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
    START_TIMER;

    uint32_t uncompressed_size = MIN_(make_ref_vb_size, layer_size[vb->refhash_layer] - vb->refhash_start_in_layer);
    uint32_t *hash_data = &refhashs[vb->refhash_layer][vb->refhash_start_in_layer / sizeof (uint32_t)];

    // // convert NO_GPOS to 0 - compresses a lot better and faster - we convert back in refhash_uncompress_one_vb
    // for (uint32_t i=0 ; i < uncompressed_size / sizeof (uint32_t); i++)
    //     if (hash_data[i] == 0xffffffff) 
    //         hash_data[i] = 0; 

    SectionHeaderRefHash header = { .section_type          = SEC_REF_HASH, 
                                    .codec                 = CODEC_RANW, // Much!! faster than LZMA (compress and uncompress), 8% worse compression of human refs, MUCH better on small refs
                                    .data_uncompressed_len = BGEN32 (uncompressed_size),
                                    .vblock_i              = BGEN32 (vb->vblock_i),
                                    .magic                 = BGEN32 (GENOZIP_MAGIC),
                                    .num_layers            = (uint8_t)num_layers,
                                    .layer_i               = (uint8_t)vb->refhash_layer,
                                    .layer_bits            = (uint8_t)layer_bits[vb->refhash_layer],
                                    .start_in_layer        = BGEN32 (vb->refhash_start_in_layer)     };

    comp_compress (vb, NULL, &vb->z_data, &header, (char*)hash_data, NO_CALLBACK, "SEC_REF_HASH");

    if (flag.show_ref_hash) 
        iprintf ("vb_i=%u Compressing SEC_REF_HASH num_layers=%u layer_i=%u layer_bits=%u start=%u size=%u bytes size_of_disk=%u bytes\n", 
                 vb->vblock_i, header.num_layers, header.layer_i, header.layer_bits, vb->refhash_start_in_layer, uncompressed_size, BGEN32 (header.data_compressed_len) + (uint32_t)sizeof (SectionHeaderRefHash));

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.

    COPY_TIMER(refhash_compress_one_vb);
}

// ZIP-FASTA-make-reference: called by main thread in zip_write_global_area
void refhash_compress_refhash (void)
{
    START_TIMER;

    next_task_layer = 0;
    next_task_start_within_layer = 0;

    dispatcher_fan_out_task ("compress_refhash", NULL, 0, "Writing hash table (this can take several minutes)...", 
                             false, false, false, 0, 100, true,
                             refhash_prepare_for_compress, 
                             refhash_compress_one_vb, 
                             zfile_output_processed_vb);

    {   
        START_TIMER;
        refhash_digest = digest_do (STRb(refhash_buf), !flag.md5, "refhash"); 
        COPY_TIMER_EVB (refhash_compress_digest); 
    }

    COPY_TIMER_EVB (refhash_compress_refhash);
}

Digest refhash_get_digest (void)
{
    return refhash_digest;
}

// ----------------------------------------------------------------------------------------------
// stuff related to ZIPping fasta and fastq using an external reference that includes the refhash
// ----------------------------------------------------------------------------------------------

ConstBitsP refhash_get_emoneg (void)
{
    ASSERTNOTEMPTY (emoneg_buf);
    return (BitsP)&emoneg_buf;
}

// private variables used for generating emoneg
static Reference ref_reverse_compliment_genome_ref = 0; // ref_generate_reverse_complement_genome is called from the main thread so no thread safety issues
static ConstBitsP genome;
static PosType64 genome_nbases;
#define REV_CODEC_GENOME_BASES_PER_THREAD (1ULL << 27) // 128Mbp

static void refhash_revcomp_genome_prepare (VBlockP vb)
{
    vb->ref = ref_reverse_compliment_genome_ref;
    if ((uint64_t)(vb->vblock_i-1) * REV_CODEC_GENOME_BASES_PER_THREAD < genome_nbases)
        vb->dispatch = READY_TO_COMPUTE;
    else
        vb->dispatch = DATA_EXHAUSTED;
}

static void refhash_revcomp_genome_do (VBlockP vb)
{
    bits_reverse_complement_aligned ((BitsP)&emoneg_buf, genome, 
                                     (uint64_t)(vb->vblock_i-1) * REV_CODEC_GENOME_BASES_PER_THREAD, 
                                     REV_CODEC_GENOME_BASES_PER_THREAD);

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.
}


// Generate an in-memory revcomp copy of genome. Only needed for aligner in ZIP.
static void refhash_generate_emoneg (Reference ref)
{
    START_TIMER;
    if (buf_is_alloc (&emoneg_buf)) return; // already generated

    SAVE_FLAGS_AUX((rom)0); // silence some flags if --xthreads

    ref_get_genome (ref, &genome, NULL, &genome_nbases);
    emoneg_buf.can_be_big = true; 
    buf_alloc_bits (evb, &emoneg_buf, genome_nbases * 2, NOINIT, 0, "emoneg_buf");

    ref_reverse_compliment_genome_ref = ref;
    dispatcher_fan_out_task ("generate_rev_comp_genome", NULL, 0, 0, true, false, false, 0, 10, true,
                             refhash_revcomp_genome_prepare, 
                             refhash_revcomp_genome_do, 
                             NO_CALLBACK);

    
    // cleanup 
    ref_reverse_compliment_genome_ref = 0;
    genome = 0;
    genome_nbases = 0;
    
    RESTORE_FLAGS;

    COPY_TIMER_EVB (refhash_generate_emoneg);
}

void refhash_set_digest (Digest digest)
{
    refhash_digest = digest;
}

// entry point of compute thread of refhash decompression
static void refhash_uncompress_one_vb (VBlockP vb)
{
    START_TIMER;

    SectionHeaderRefHashP header = (SectionHeaderRefHashP )vb->z_data.data;
    uint32_t start = BGEN32 (header->start_in_layer);
    uint32_t size  = BGEN32 (header->data_uncompressed_len);
    uint32_t layer_i = header->layer_i;

    if (flag.show_ref_hash) // before the sanity checks
        iprintf ("vb_i=%u Uncompressing SEC_REF_HASH num_layers=%u layer_i=%u layer_bits=%u start=%u size=%u bytes size_of_disk=%u bytes\n", 
                 vb->vblock_i, header->num_layers, layer_i, header->layer_bits, start, size, BGEN32 (header->data_compressed_len) + (uint32_t)sizeof (SectionHeaderRefHash));

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

    COPY_TIMER (refhash_uncompress_one_vb);
}

static void refhash_read_one_vb (VBlockP vb)
{
    START_TIMER;

    buf_alloc (vb, &vb->z_section_headers, 0, 1, int32_t, 0, "z_section_headers"); // room for 1 section header

    static Section sec;    
    if (vb->vblock_i == 1) sec = NULL; // reset for each file

    if (!sections_next_sec (&sec, SEC_REF_HASH))
        return; // no more refhash sections

    int32_t section_offset = zfile_read_section (z_file, vb, sec->vblock_i, &vb->z_data, "z_data", sec->st, sec);

    // TO DO: if refhash_load sets num_layers to a value other than header.num_layers, we don't need
    // to even read this section (currently never happens) (bug 348)
    if (((SectionHeaderRefHashP)vb->z_data.data)->layer_i >= num_layers)
        return; 

    BNXT (int32_t, vb->z_section_headers) = section_offset;

    vb->dispatch = READY_TO_COMPUTE;

    if (flag.debug_or_test) buflist_test_overflows(vb, __FUNCTION__); 

    COPY_TIMER (refhash_read_one_vb);
}

// ZIP: needed before using aligner
void refhash_load (Reference ref)
{
    START_TIMER;
    if (refhash_buf.count) return; // already loaded

    // emoneg is needed iff aligner is used (in ZIP)
    if (flag.aligner_available)
        refhash_generate_emoneg (ref);

    Section sec = sections_last_sec (SEC_REF_HASH, HARD_FAIL);
    SectionHeaderRefHash header = zfile_read_section_header (evb, sec, SEC_REF_HASH).ref_hash;
    num_layers = header.num_layers; // set global

    uint32_t base_layer_bits = header.layer_bits + header.layer_i; // layer_i=0 is the base layer, layer_i=1 has 1 bit less etc
    refhash_initialize_layers (base_layer_bits);

    if (!ref_cache_is_cached (ref)) { // no cache, or loading to cache, but not already cached
        dispatcher_fan_out_task ("load_refhash", ref_get_filename (ref),
                                0, (ref_cache_is_populating(ref) ? "Caching reference file" : "Reading reference file"), // same message as in ref_load_stored_reference
                                true, flag.test, false, 0, 100, true,
                                refhash_read_one_vb, 
                                refhash_uncompress_one_vb, 
                                NO_CALLBACK);
    
        if (flag.show_cache) iprint0 ("show-cache: done reading refhash from disk\n");
    }

    // calculate in-memory digest of loaded genome, and compare it to the value calculated by
    // --make-reference, stored in GENOZIP_HEADER.refhash_digest
    // This is slow: for a human genome, digest takes about 0.7 seconds (adler) or 7 seconds (md5), so we test the digest sparingly. 
    // This is low risk: if refhash is corrupted, it will affect compression, but not data integrity (see bug 825)
    if (ref_get_genozip_version(ref) >= 15 &&              // GENOZIP_HEADER.refhash_digest is available since v15
        (!ref_cache_is_cached (ref) || flag.show_cache) && // test only on loading from disk or on show-cache
        !flag.fast) {                                      // skip if --fast. 

        START_TIMER;

        if (flag.show_cache) iprintf ("show-cache: calculating refhash digest %s\n", ref_cache_is_cached (ref) ? "(done only due to --show-cache)" : "");

        bool is_adler = ref_is_digest_adler (ref);
        Digest digest = digest_do (STRb(refhash_buf), is_adler, "refhash"); 

        ASSERT (digest_is_equal (digest, refhash_digest), "Bad reference file: In-memory digest of refhash is %s, different than calculated by make-reference: %s",
                digest_display_(digest, is_adler).s, digest_display_(refhash_digest, is_adler).s);

        if (flag.show_cache) iprint0 ("show-cache: verified refhash digest\n");

        COPY_TIMER_EVB (refhash_load_digest);
    }
    else
        if (flag.show_cache) iprintf ("show-cache: loaded refhash without verifying digest (version %u)\n", ref_get_genozip_version(ref));

    refhash_buf.count = true; // loaded

    COPY_TIMER_EVB (refhash_load);
}

void refhash_load_standalone (void)
{
    if (ref_cache_is_cached (gref))  // already loaded 
        refhash_generate_emoneg (gref);

    flag.reading_reference = gref; // tell file.c and fasta.c that this is a reference

    TEMP_VALUE (command, PIZ);
    TEMP_VALUE (z_file, NULL);   // save z_file and txt_file in case we are called from sam_seg_finalize_segconf
    TEMP_VALUE (txt_file, NULL);
    CLEAR_FLAG (test);

    z_file = file_open_z_read (ref_get_filename (gref));    

    zfile_read_genozip_header (0, HARD_FAIL);

    refhash_load (gref);

    file_close (&z_file);
    file_close (&txt_file); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open_z called from txtheader_piz_read_and_reconstruct.
    
    RESTORE_FLAG (test);
    RESTORE_VALUE (command);
    RESTORE_VALUE (z_file);
    RESTORE_VALUE (txt_file);
    
    flag.reading_reference = NULL;
}

void refhash_destroy (void)
{
    if (refhash_buf.type == BUF_UNALLOCATED) return;

    if (flag.show_cache && refhash_buf.type == BUF_SHM) 
        iprint0 ("show-cache: destroy rehash_buf attached to shm\n");

    buf_destroy (refhash_buf);
    buf_destroy (emoneg_buf);

    FREE (refhashs);

    flag.aligner_available = false;
}

