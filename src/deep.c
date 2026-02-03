// ------------------------------------------------------------------
//   deep.c
//   Copyright (C) 2023-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "deep.h"
#include "vblock.h"
#include "libdeflate_1.19/libdeflate.h"

uint8_t num_hash_bits; // ZIP: number of bits of seq_hash used for deep_hash_by_* (i.e. hash table is of size 2^num_hash_bits) (note: cannot be in z_file because bamass may span several z_files)

#define MAX_AUTO_READ_LEN 20000 // not too long, so that the "longer reads" code path also gets some mileage

uint64_t deep_qname_hash (VBlockP vb, QType q, STRp(qname), thool is_last, uint32_t *uncanonical_suffix_len) 
{
    if (is_last != unknown && segconf.is_interleaved && vb->line_i % 2) is_last = !is_last;

    return qname_calc_hash (q, COMP_NONE, STRa(qname), is_last, true, CRC64, uncanonical_suffix_len); 
}

// hash of a SEQ field in the forward direction 
// note: I tested crc32 after converting seq to 2-bit. No advantage - Almost identical linked-list-length histogram.
uint32_t deep_seq_hash (VBlockP vb, STRp(seq), bool is_revcomp)
{
    char short_read_data[MAX_AUTO_READ_LEN]; 

    if (is_revcomp) {
        char *rc; 
        
        // short enough reads - use automatic allocation
        if (seq_len <= MAX_AUTO_READ_LEN) 
            rc = short_read_data;

        // longer reads - use heap allocation
        else {
            ASSERTNOTINUSE (vb->scratch);
            buf_alloc (vb, &vb->scratch, 0, seq_len, char, 0, "scratch");
            rc = B1STc (vb->scratch);
        }
        str_revcomp (rc, STRa(seq));
        seq = rc;
    }

    uint32_t hash = crc32 (0, seq, seq_len);

    if (is_revcomp && seq_len > MAX_AUTO_READ_LEN)
        buf_free (vb->scratch);

    return hash;
}

// hash of a QUAL field in the forward direction 
uint32_t deep_qual_hash (VBlockP vb, STRp(qual), bool is_revcomp)
{
    char short_read_data[MAX_AUTO_READ_LEN] = {};

    if (is_revcomp) {
        char *rev; 
        
        // short enough reads - use automatic allocation
        if (qual_len <= MAX_AUTO_READ_LEN) 
            rev = short_read_data;

        // longer reads - use heap allocation
        else {
            ASSERTNOTINUSE (vb->scratch);
            buf_alloc (vb, &vb->scratch, 0, qual_len, char, 0, "scratch");
            rev = B1STc (vb->scratch);
        }
        str_reverse (rev, STRa(qual));
        qual = rev;
    }

    uint32_t hash = crc32 (0, STRa(qual));
    
    // note: qual_hash=0 means "QUAL not hashed", so both crc32=0 and crc32=1 get mapped to hash=1
    if (hash == 0) hash = 1;

    if (is_revcomp && qual_len > MAX_AUTO_READ_LEN)
        buf_free (vb->scratch);

    return hash;
}
