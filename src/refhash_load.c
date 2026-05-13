// ------------------------------------------------------------------
//   refhash_load.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "compressor.h" // codec.h included in compressor.h must be included before reference.h
#include "file.h"
#include "dispatcher.h"
#include "zfile.h"
#include "refhash_friend.h"
#include "threads.h"
#include "reference.h"

bool refhash_is_flat = true;     // true iff version is at least 15.0.81
int gpos_bytes=0;                // 4 or 5
int bits_per_hash_out=0;         // FLAT: determined by genome size. // hash table has 2^ entries of this, each entry is gpos_bytes. this can be changed, and hash tables over 4 GB are ok. LAYERED: =layer_bits[0]
int bases_per_hash=0;            // length of input to the hash function
int bits_per_hash_in=0;          // bases_per_hash * 2
int hash_shift=0;                // (64 - bits_per_hash_out)
MakeRefSize make_ref_size;       // of reference currently loaded
Buffer refhash_buf = {};         // One buffer that includes all layers
Digest refhash_digest = DIGEST_NONE;

// ---------------------------------------------------------------
// loading refhash stuff stuff
// ---------------------------------------------------------------

Digest refhash_get_digest (void)        { return refhash_digest;   }

// main thread
static void refhash_read_one_vb (VBlockP vb)
{
    START_TIMER;

    buf_alloc (vb, &vb->z_section_headers, 0, 1, int32_t, 0, "z_section_headers"); // room for 1 section header

    static Section sec;    
    if (vb->vblock_i == 1) sec = NULL; // reset for each file

    if (!sections_next_sec (&sec, SEC_REF_HASH))
        return; // no more refhash sections

    int32_t section_offset = zfile_read_section (z_file, vb, sec->vblock_i, &vb->z_data, "z_data", sec->st, sec);

    BNXT (int32_t, vb->z_section_headers) = section_offset;

    vb->dispatch = READY_TO_COMPUTE;

    if (flag.debug_or_test) buflist_test_overflows(vb, __FUNCTION__); 

    COPY_TIMER (refhash_read_one_vb);
}

// entry point of compute thread of refhash decompression
static void refhash_uncompress_one_vb (VBlockP vb)
{
    START_TIMER;

    SectionHeaderRefHashP header = (SectionHeaderRefHashP )vb->z_data.data;

    // a hack for uncompressing to a location withing the buffer - while multiple threads are uncompressing into 
    // non-overlappying regions in the same buffer in parallel
    Buffer copy = refhash_buf;

    if (refhash_is_flat) {
        ASSERTNOTZERO (gpos_bytes);
        uint64_t start = header->first_ent * gpos_bytes;
        copy.data = Bc(refhash_buf, start);
    }

    else  // up to 15.0.80
        copy.data = (char *)B32(refhash_buf, LAYER_START[header->OLD.layer_i]) + BGEN32 (header->OLD.start_in_layer);
    
    zfile_uncompress_section (vb, header, &copy, NULL, 0, SEC_REF_HASH);

    vb_set_is_processed (vb); // tell dispatcher this thread is done and can be joined.

    COPY_TIMER (refhash_uncompress_one_vb);
}

// when reading a reference file: set info from reference file's GENOZIP_HEADER section
void refhash_set_ref_file_info (Digest digest, uint8_t ref_bases_per_hash, uint8_t ref_bits_per_hash_out, uint8_t ref_gpos_bytes, MakeRefSize my_make_ref_size)
{
    refhash_digest  = digest;
    make_ref_size   = my_make_ref_size;
    refhash_is_flat = (make_ref_size != MAKE_REF_OLD); // true since 15.0.81
    
    if (refhash_is_flat) { 
        gpos_bytes        = ref_gpos_bytes;
        bases_per_hash    = ref_bases_per_hash;
        bits_per_hash_out = ref_bits_per_hash_out;
        hash_shift        = (64 - bits_per_hash_out);
    }
    else { // up to 15.0.80
        gpos_bytes        = 4;
        bases_per_hash    = 14;
        bits_per_hash_out = 28;
    }

    bits_per_hash_in = bases_per_hash * 2;
}

// check if this reference file has refhash
bool refhash_exists (void)
{
    ASSERTNOTZERO (make_ref_size); // this function only works after reference is loaded

    return make_ref_size != MAKE_REF_MINIMAL;
}

// ZIP: needed before using aligner
void refhash_load (void)
{
    if (refhash_buf.count) return; // already loaded

    Version ref_ver = ref_get_genozip_ver();
    
    refhash_buf.can_be_big = true; // avoid warning
    
    if (refhash_buf.type != BUF_SHM) { // if cached refhash_buf is initialized in ref_cache_initialize_genome()
        if (!refhash_is_flat)     buf_alloc_exact (evb, refhash_buf, 256 MB + 128 MB + 64 MB + 32 MB, uint32_t, "refhash_buf"); // layered refhash - up to 15.0.80
        else if (gpos_bytes == 4) buf_alloc_exact (evb, refhash_buf, ref_hash_len, uint32_t, "refhash_buf"); 
        else                      buf_alloc_exact (evb, refhash_buf, ref_hash_len, uint40_t, "refhash_buf"); 
    
        refhash_buf.len *= gpos_bytes; // len counts bytes
    }
    
    if (!ref_cache_is_cached ()) { // no cache, or loading to cache, but not already cached
        dispatcher_fan_out_task (TASK_LOAD_REFHASH, ref_get_filename(),
                                 0, JOIN_OUT_OF_ORDER, false, 0, 100,
                                 refhash_read_one_vb, 
                                 refhash_uncompress_one_vb, 
                                 NO_CALLBACK);
    
        if (flag.show_cache) iprint0 ("show-cache: done reading refhash from disk\n");
    }

    // calculate in-memory digest of loaded genome, and compare it to the value calculated by
    // --make-reference, stored in GENOZIP_HEADER.refhash_digest
    // This is slow: for a human genome, digest takes about 0.7 seconds (adler) or 7 seconds (md5), so we test the digest sparingly. 
    // This is low risk: if refhash is corrupted, it will affect compression, but not data integrity (see bug 825)
    if (ref_ver.major >= 15 &&                         // GENOZIP_HEADER.refhash_digest is available since v15
        (!ref_cache_is_cached() || flag.show_cache) && // test only on loading from disk or on show-cache
        !flag.fast) {                                  // skip if --fast. 

        START_TIMER;

        if (flag.show_cache) iprintf ("show-cache: calculating refhash digest %s\n", ref_cache_is_cached() ? "(done only due to --show-cache)" : "");

        DigestAlg alg = ref_get_genome_digest_alg();

        Digest digest = digest_do (STRb(refhash_buf), alg, "refhash"); 
        
        // a different refhash_digest is an indication that the reference file is not the file that is cached (perhaps is was 
        // re-made: every make results in a different refhash_digest), so we remove the shm to avoid a debuggability nightmare
        if (!digest_is_equal (digest, refhash_digest)) {
            ref_cache_remove_do (true, false); 
            ABORT ("Bad reference file: In-memory %s digest of refhash is %s, different than calculated by make-reference: %s",
                   digest_alg_name(alg), digest_display_(digest, alg).s, digest_display_(refhash_digest, alg).s);
        }

        if (flag.show_cache) iprint0 ("show-cache: verified refhash digest\n");

        COPY_TIMER_EVB (refhash_load_digest);
    }
    else
        if (flag.show_cache) iprintf ("show-cache: loaded refhash without verifying digest (version %s)\n", STRver(ref_ver).s);

    refhash_buf.count = true; // loaded
}

void refhash_load_standalone (void)
{
    flag.reading_reference = true; // tell file.c and fasta.c that this is a reference

    TEMP_VALUE (command, PIZ);
    TEMP_VALUE (z_file, NULL);   // save z_file and txt_file in case we are called from sam_seg_finalize_segconf
    TEMP_VALUE (txt_file, NULL);
    CLEAR_FLAG (test);

    z_file = file_open_z_read (ref_get_filename());    

    zfile_read_genozip_header (0, HARD_FAIL);

    refhash_load();

    file_close (&z_file);
    file_close (&txt_file); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open_z called from txtheader_piz_read_and_reconstruct.
    
    RESTORE_FLAG (test);
    RESTORE_VALUE (command);
    RESTORE_VALUE (z_file);
    RESTORE_VALUE (txt_file);
    
    flag.reading_reference = false;
}

void refhash_destroy (void)
{
    if (refhash_buf.type == BUF_UNALLOCATED) return;

    if (flag.show_cache && refhash_buf.type == BUF_SHM) 
        iprint0 ("show-cache: destroy rehash_buf attached to shm\n");

    buf_destroy (refhash_buf);

    flag.aligner_available = false;
}

// returns size in bytes of refhash
uint64_t refhash_get_refhash_size (void)
{
    if (refhash_is_flat) // data from GenozipHeader
        return ref_hash_len * gpos_bytes;

    else  // up 15.0.80 - in original code controled by variables, but variables never changed
        return 1 GB + 512 MB + 256 MB + 128 MB;
}
