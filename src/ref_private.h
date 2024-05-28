// ------------------------------------------------------------------
//   ref_private.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "reference.h"
#include "sections.h"
#include "contigs.h"

typedef struct {
    Bits *is_set;
    PosType64 first_bit, len;
} RegionToSet;

typedef enum { CHROM_STYLE_UNKNOWN, CHROM_STYLE_chr22, CHROM_STYLE_22 } RefChromeStyle;

typedef enum { CACHE_INITITAL, CACHE_READY/*shm read-only*/, CACHE_POPULATING/*shm read-write*/, CACHE_NONE } RefCacheState;

// a reference cache is a shared memory segment consisting of a RefCache struct, followed by the genome, follewed by the refhash
typedef struct ref_cache { 
    uint32_t magic;               // set to GENOZIP_MAGIC - used to detect whether this shm segment is a Genozip reference cache
    uint32_t creator_pid;
    uint64_t creation_ts;         // if set - cache is ready OR currently being populated
    uint64_t shm_size;                
    uint8_t genozip_version;      // Genozip version that created this cache (NOT genozip version of reference file)
    bool is_populated;            // if set - cache is ready
    bool terminate_holder;        // Windows: a message to the holder process that it can terminate now
    uint8_t unused[5];            // start genome_data 64b-word-aligned
    char ref_basename[256];       // nul-terminated basename of reference filename, or <unused> if too long
    char genome_data[0];
} RefCache;

typedef struct RefStruct {
    // file 
    rom filename;                 // filename of external reference file
    Digest genome_digest;         // v15: digest of genome as it is loaded to memory. Up to v14: MD5 of original FASTA file (buggy)
    bool is_adler;                // true if genome_digest is Adler32, false if it MD5
    uint8_t genozip_version;
    bool is_primary;
    bool is_filename_allocated;
    
    // features
    RefChromeStyle chrom_style;

    // ZIP and PIZ, internal or external reference ranges. If in ZIP-INTERNAL we have REF_NUM_DENOVO_RANGES Range's - each allocated on demand. In all other cases we have one range per contig.
    Buffer ranges; 
    #define rtype param
    
    Buffer genome_buf, genome_is_set_buf;
    BitsP genome,                 // the genome in 2-bit representation. attached to shared memory or allocated privately 
          genome_is_set;          // 1 bit per reference base, indicates if base is needed for reconstructing current file. 

    PosType64 genome_nbases;

    Buffer ref_external_ra;       // Random Access data of the external reference file
    Buffer ref_file_section_list; // Section List of the external reference file

    // PIZ: an array of RegionToSet - list of non-compacted regions, for which we should set is_set bit to 1
    Buffer region_to_set_list; 
    SPINLOCK (region_to_set_list_spin);

    // ZIP/PIZ: random_access data for the reference sections stored in a target genozip file
    Buffer stored_ra;
    SPINLOCK (stored_ra_spin);    // ZIP only

    char *ref_fasta_name;

    bool external_ref_is_loaded;

    // contigs loaded from a reference file
    ContigPkg ctgs;

    // lock stuff
    Buffer genome_muteces;        // one mutex per GENOME_BASES_PER_MUTEX (64K) bases - protects genome->is_set
    Buffer genome_mutex_names;

    // iupac stuff
    Buffer iupacs_buf; 

    // reference cache
    RefCacheState cache_state;
#ifndef _WIN32
    #define CACHE_SHM_NONE -1
    int cache_shm;
#else
    #define CACHE_SHM_NONE NULL
    void *cache_shm; // Windows HANDLE is defined as void *
#endif
    RefCache *cache; // cache consists of a block memory containing a RefCache struct, followed by the genome, follewed by the refhash

} RefStruct;

extern void ref_make_prepare_ranges_for_compress (void);
extern void ref_make_prepare_one_range_for_compress (VBlockP vb);

// lock stuff
extern void ref_lock_initialize (Reference ref);
extern void ref_lock_free (Reference ref);

// contigs stuff
extern rom ref_contigs_get_name_by_ref_index (Reference ref, WordIndex chrom_index, pSTRp(snip), PosType64 *gpos);
extern const Contig *ref_contigs_get_contig_by_ref_index (Reference ref, WordIndex chrom_index, FailType soft_fail);
extern WordIndex ref_seg_get_alt_chrom (VBlockP vb);
extern void ref_contigs_compress_ref_make (Reference ref);
extern void ref_contigs_compress_stored (Reference ref);
extern void ref_make_calculate_digest (void);

// cache stuff
extern bool ref_cache_initialize_genome (Reference ref);
extern void ref_cache_done_populating (Reference ref);
extern void ref_cache_remove_do (Reference ref, bool cache_exists, bool verbose);
