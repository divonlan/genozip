// ------------------------------------------------------------------
//   ref_private.h
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "reference.h"
#include "sections.h"
#include "mutex.h"
#include "contigs.h"

typedef struct {
    Bits *is_set;
    PosType64 first_bit, len;
} RegionToSet;

typedef enum { CHROM_STYLE_UNKNOWN, CHROM_STYLE_chr22, CHROM_STYLE_22 } RefChromeStyle;

typedef struct RefStruct {

    // file 
    rom filename; // filename of external reference file
    Digest genome_digest; // v15: digest of genome as it is loaded to memory. Up to v14: MD5 of original FASTA file (buggy)
    bool is_adler;        // true if genome_digest is Adler32, false if it MD5
    uint8_t genozip_version;
    bool is_primary;
    bool is_filename_allocated;
    
    // features
    RefChromeStyle chrom_style;

    // ZIP and PIZ, internal or external reference ranges. If in ZIP-INTERNAL we have REF_NUM_DENOVO_RANGES Range's - each allocated on demand. In all other cases we have one range per contig.
    Buffer ranges; 
    #define rtype param
    
    Buffer genome_buf, emoneg_buf, genome_is_set_buf;
    BitsP genome, emoneg/*reverse compliment*/, genome_is_set;

    PosType64 genome_nbases;

    Buffer ref_external_ra;       // Random Access data of the external reference file
    Buffer ref_file_section_list; // Section List of the external reference file

    // PIZ: an array of RegionToSet - list of non-compacted regions, for which we should set is_set bit to 1
    Buffer region_to_set_list; 
    SPINLOCK (region_to_set_list_spin);

    // ZIP/PIZ: random_access data for the reference sections stored in a target genozip file
    Buffer stored_ra;
    SPINLOCK (stored_ra_spin); // ZIP only

    char *ref_fasta_name;

    bool external_ref_is_loaded;

    // contigs loaded from a reference file
    ContigPkg ctgs;

    // lock stuff
    Buffer genome_muteces; // one mutex per 64K bases - protects genome->is_set
    Buffer genome_mutex_names;

    // iupac stuff
    Buffer iupacs_buf; 

    // cache
    enum { CACHE_INITITAL, CACHE_OK, CACHE_LOADING, CACHE_NONE } cache_state;
    int cache_shm;

    // cache consists of a block memory containing a ref_cache struct, followed by the genome, follewed by the refhash
    struct ref_cache {
        uint64_t creation_ts;
        uint64_t genome_size, refhash_size;
        uint8_t ref_genozip_ver;
        bool is_populated;
        uint8_t base_layer_bits;
        uint8_t unused[61];    // start genome_data word-aligned
        char genome_data[0];
    } *cache;
    void *cache_refhash_data;  // pointer into cache data

} RefStruct;

extern void ref_make_prepare_range_for_compress (VBlockP vb);

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
