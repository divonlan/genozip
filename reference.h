// ------------------------------------------------------------------
//   reference.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REFERENCE_INCLUDED
#define REFERENCE_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "md5.h"

// reference sequences - one per range of 1MB. ranges (chrom, pos) are mapped here with a hash function. In the rare case two unrelated ranges
// are mapped to the same entry - only the first range will have a reference, and the other ones will not. this will hurt the compression ratio,
// but not the correctness.
// Thread safety: "ref" is set atomically as the last set in initialization. If its set, then it is correct and will never change. If it appears to be
// not set yet, it is necesarry to lock the mutex and test again, and initialize if still not set.

#define REF_NUM_RANGES (1 << 20)
#define REF_NUM_SITES_PER_RANGE (1 << 20) // 1 Mbp

// values of flag_reference
#define REF_NONE      0
#define REF_INTERNAL  1 // ZIP: reference is calculated from file(s) data
#define REF_EXTERNAL  2 // ZIP & PIZ: user specified --reference
#define REF_EXT_STORE 3 // ZIP: user specified --REFERENCE
#define REF_STORED    4 // PIZ: file contains REFERENCE sections (user cannot specify --reference)

typedef struct Range {
    char *ref;                    // actual reference data
    const char *chrom_name;
    unsigned chrom_name_len;
    uint32_t range_i;
    uint32_t num_set;             // number of bases set in this range
    uint32_t uncompressed_len;    // final length of range ready for compression
    int64_t first_pos, last_pos;  // [first_pos,last_pos] include all set locii in ZIP-INTERNAL (might also include locii that were not set) or entire data read from the reference in ZIP-EXTERNAL and PIZ
    uint64_t ref_size;            // number of bytes allocated to ref (at most REF_NUM_SITES_PER_RANGE)
    bool is_compacted;
    pthread_mutex_t mutex;

    bool is_accessed;             // has this range been encountered in the file (if not, we won't store it)
} Range;

extern void ref_compress_ref (void);
extern void ref_uncompress_all_ranges (void);
extern void ref_set_reference (const char *filename);
extern void ref_set_md5 (Md5Hash md5);
extern void ref_read_external_reference (void);
extern void ref_cleanup_memory (void);
extern MemStats ref_memory_consumption (void);
extern const char *ref_get_ref (VBlockP vb, uint64_t pos, uint32_t ref_consumed);
extern void ref_consume_ref_fasta_global_area (void);
extern void ref_get_contig (const char *chrom_name, unsigned chrom_name_len, const char **snip, int32_t *word_index);
extern void ref_set_ref_from_external_data (VBlockP vb, uint64_t start_pos, ConstBufferP data, uint64_t data_start, bool might_span_multiple_vbs);
extern void ref_get_contigs (ConstBufferP *out_contig_dict, ConstBufferP *out_contig_words);
extern int64_t ref_max_pos_of_chrom (uint32_t chrom_word_index);

typedef enum { RGR_MUST_BE_NEW, RGR_MUST_EXIST, RGR_CAN_BE_NEW_OR_EXISTING } RgrMode;
extern Range *ref_get_range (VBlockP vb, uint32_t range_i, RgrMode mode, uint32_t size);

// globals
extern bool ref_flag_reading_reference;
extern const char *ref_filename;
extern Md5Hash ref_md5;

#endif
