// ------------------------------------------------------------------
//   ref_private.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef REF_PRIVATE_INCLUDED
#define REF_PRIVATE_INCLUDED

#include "genozip.h"
#include "reference.h"
#include "sections.h"
#include "mutex.h"

#pragma GENDICT_PREFIX REF
#pragma GENDICT REF_CONTIG=DTYPE_FIELD=CONTIG 

typedef struct {
    BitArray *is_set;
    PosType first_bit, len;
} RegionToSet;

typedef struct RefStruct {

    // file
    const char *filename; // filename of external reference file
    const char *cache_fn;
    Digest file_md5;
    uint8_t genozip_version;
    
    // ZIP and PIZ, internal or external reference ranges. If in ZIP-INTERNAL we have REF_NUM_DENOVO_RANGES Range's - each allocated on demand. In all other cases we have one range per contig.
    Buffer ranges; 

    Buffer genome_buf, emoneg_buf, genome_is_set_buf, genome_cache;
    BitArray *genome, *emoneg /* reverse compliment */, *genome_is_set;

    PosType genome_nbases;

    Buffer ref_external_ra;       // Random Access data of the external reference file
    Buffer ref_file_section_list; // Section List of the external reference file

    // PIZ: an array of RegionToSet - list of non-compacted regions, for which we should set is_set bit to 1
    Buffer region_to_set_list; 
    SPINLOCK (region_to_set_list_spin);

    // ZIP/PIZ: random_access data for the reference sections stored in a target genozip file
    Buffer stored_ra;
    SPINLOCK (stored_ra_spin); // ZIP only

    char *ref_fasta_name;

    ThreadId ref_cache_creation_thread_id;
    VBlockP cache_create_vb;
    bool external_ref_is_loaded;

    // contigs loaded from a reference file
    Buffer loaded_contigs; // array of RefContig
    Buffer loaded_contigs_dict;
    Buffer loaded_contigs_by_name; // an array of uint32 of indexes into loaded_contigs - sorted by alphabetical order of the snip in contig_dict
    Buffer loaded_contigs_by_LN;

    // lock stuff
    Mutex *genome_muteces; // one spinlock per 16K bases - protects genome->is_set
    uint32_t genome_num_muteces;
    char *genome_mutex_names;

    // iupac stuff
    Buffer iupacs_buf; 

} RefStruct;

extern void ref_make_prepare_range_for_compress (VBlockP vb);

// lock stuff
extern void ref_lock_initialize_loaded_genome (Reference ref);
extern void ref_lock_initialize_denovo_genome (Reference ref);
extern void ref_lock_free (Reference ref);

// contigs stuff
extern const char *ref_contigs_get_chrom_snip (Reference ref, WordIndex chrom_index, const char **snip, uint32_t *snip_len);
extern const RefContig *ref_contigs_get_contig (Reference ref, WordIndex chrom_index, bool soft_fail);
extern PosType ref_contigs_get_genome_nbases (Reference ref);
extern void ref_contigs_generate_data_if_denovo (Reference ref);
extern WordIndex ref_seg_get_alt_chrom (VBlockP vb);
extern void ref_contigs_compress (Reference ref);
extern WordIndex ref_contigs_get_by_accession_number (Reference ref, const char *ac, unsigned ac_len);
extern WordIndex ref_contigs_get_by_uniq_len (Reference ref, PosType LN);

#define ranges_type(ref) (ref)->ranges.param

#define IS_REF_INTERNAL(f) (((f)->data_type == DT_SAM || (f)->data_type == DT_BAM) && (f)->z_flags.dts_ref_internal)

#define ROUNDUP32(x) (((x) + 32) & ~(typeof(x))0x1f) // round up to the nearest 32
#define ROUNDDOWN32(x) ((x)      & ~(typeof(x))0x1f) // round down to the nearest 32

#define ROUNDUP64(x) (((x) + 63) & ~(typeof(x))0x3f) // round up to the nearest 64
#define ROUNDDOWN64(x) ((x)      & ~(typeof(x))0x3f) // round down to the nearest 64

#endif
