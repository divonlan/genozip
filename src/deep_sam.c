// ------------------------------------------------------------------
//   deep_sam.c
//   Copyright (C) 2023-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "libdeflate/libdeflate.h"

// -------------------------------------------------
// SAM seg: calculate hash values into vb->deep_data
// -------------------------------------------------

// ZIP: Per-line Deep data (excluding SECONDARY and SUPPLEMENTARY lines)
typedef struct {
    uint32_t qname_hash;
    uint32_t seq_hash;
    uint32_t qual_hash;
} ZipVbDeep;

typedef struct {
    uint32_t qname_hash;
    uint32_t seq_hash;
    uint32_t qual_hash;
    VBIType  sam_vb_i;   // SAM VB in which of alignment. TO DO: store global line instead of VB - and z_file->num_lines_per_sam_vb
    uint32_t sam_line_i; // line_i within SAM VB
    uint32_t next;       // linked list
} ZipZDeep;

#define NO_NEXT     0xffffffff
#define MAX_ENTRIES 0xfffffffe // our current implementation is limited to 4G reads

void deep_sam_seg_initialize (VBlockSAMP vb)
{
    buf_alloc (vb, &vb->deep_data, 0, vb->lines.len, ZipVbDeep, 0, "deep_data");
}

void deep_sam_set_QNAME_hash (VBlockSAMP vb, ZipDataLineSAM *dl, uint32_t qname_hash)
{
    buf_alloc (vb, &vb->deep_data, 1, 0, ZipVbDeep, CTX_GROWTH, "deep_data");

    // note: qname_hash=0 means "line is not primary", so both crc32=0 and crc32=1 get mapped to hash=1
    BNXT(ZipVbDeep, vb->deep_data).qname_hash = FLAG_IS_PRIMARY(dl->FLAG) ? MAX_(1, qname_hash) : 0;
}

// set hash of forward SEQ (i.e. as it would appear in the FASTQ file) 
void deep_sam_set_SEQ_hash (VBlockSAMP vb,ZipDataLineSAM *dl, STRp(textual_seq))
{
    if (!FLAG_IS_PRIMARY(dl->FLAG)) return;

    if (dl->FLAG.rev_comp) {
        ASSERTNOTINUSE (vb->scratch);
        buf_alloc (vb, &vb->scratch, 0, textual_seq_len, char, 0, "scratch");

        str_revcomp_in_out (B1STc (vb->scratch), STRa(textual_seq));
        textual_seq = B1STc (vb->scratch);
    }

    BLST(ZipVbDeep, vb->deep_data)->seq_hash = crc32 (0, STRa(textual_seq));

    if (dl->FLAG.rev_comp)
        buf_free (vb->scratch);
}

// hash of forward QUAL (i.e. as it would appear in the FASTQ file) for Deep 
void deep_sam_set_QUAL_hash (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qual))
{
    if (!FLAG_IS_PRIMARY(dl->FLAG)) return;

    if (vb->has_qual && !segconf.has_bqsr) {
        if (dl->FLAG.rev_comp) {
            ASSERTNOTINUSE (vb->scratch);
            buf_alloc (vb, &vb->scratch, 0, qual_len, char, 0, "scratch");

            str_reverse (B1STc (vb->scratch), STRa(qual));
            qual = B1STc (vb->scratch);
        }

        // note: qual_hash=0 means "QUAL not hashed", so both crc32=0 and crc32=1 get mapped to hash=1
        BLST(ZipVbDeep, vb->deep_data)->qual_hash = MAX_(1, crc32 (0, STRa(qual)));

        if (dl->FLAG.rev_comp)
            buf_free (vb->scratch);
    }

    else
        BLST(ZipVbDeep, vb->deep_data)->qual_hash = 0;
}

void deep_sam_merge_deep_data (VBlockSAMP vb)
{
    #define num_hash_bits z_file->deep_hash.prm8[0]

    ASSERT (vb->lines.len == vb->deep_data.len, "Expecting vb->lines.len=%u == vb->deep_data.len=%u", vb->lines.len32, vb->deep_data.len32);

    mutex_lock (z_file->deep_populate_ents_mutex);

    if (!z_file->deep_hash.len) {
        double est_num_vbs = MAX_(1, (double)txtfile_get_seggable_size() / (double)vb->txt_data.len);
        uint32_t est_num_lines = MIN_(MAX_ENTRIES, (uint64_t)(est_num_vbs * (double)vb->lines.len)); 

        // number of bits of hash table size (maximum 31)
        num_hash_bits = MIN_(16, 32 - (int)__builtin_clzl (MIN_(0x7fffffffULL, est_num_lines * 2))); // num hash bits

        buf_alloc_exact_255 (evb, z_file->deep_hash, (1 << num_hash_bits), uint32_t, "z_file->deep_hash"); // pointer into deep_ents
        
        buf_alloc_exact (evb, z_file->deep_ents, MAX_(est_num_lines, vb->deep_data.len), ZipZDeep, "z_file->deep_ents");
        z_file->deep_ents.len = 0;
    }

    buf_alloc (evb, &z_file->deep_ents, vb->deep_data.len32, MIN_(MAX_ENTRIES, CTX_GROWTH * z_file->deep_ents.len), 
               ZipZDeep, 0, "z_file->deep_ents"); // limit growth to MAX_ENTRIES

    uint32_t ent_i = z_file->deep_ents.len32;

    ARRAY (uint32_t, hash_table, z_file->deep_hash);
    ARRAY (ZipZDeep, deep_ents, z_file->deep_ents);

    uint32_t mask = bitmask32 (num_hash_bits);
    for_buf2 (ZipVbDeep, line_deep_data, line_i, vb->deep_data) {
        if (!line_deep_data->qname_hash) continue; // case: secondary or supplementary line - qname_hash=0

        // we can only handle 4G entries, because the hash table as well as ZipZDeep.next is 32 bit.
        if (ent_i == 0xffffffff) {
            WARN_ONCE ("The number of alignments in %s exceeds the maximum supported for optimal Deep compression. The compression will not be as good as it could have. Please report this to " EMAIL_SUPPORT, txt_name);
            break;
        }

        uint32_t hash = line_deep_data->seq_hash & mask; 
        
        // case: first entry for this hash value
        if (hash_table[hash] == NO_NEXT) 
            hash_table[hash] = ent_i;

        // case: add entry add end of linked list
        uint32_t linked_ent_i = hash_table[hash];
        while (deep_ents[linked_ent_i].next != NO_NEXT) 
            linked_ent_i = deep_ents[linked_ent_i].next;

        deep_ents[linked_ent_i].next = ent_i; // extend linked list

        // create new entry
        deep_ents[ent_i++] = (ZipZDeep){ .next       = NO_NEXT, 
                                         .seq_hash   = line_deep_data->seq_hash,
                                         .qual_hash  = line_deep_data->qual_hash,
                                         .sam_vb_i   = vb->vblock_i,
                                         .sam_line_i = line_i };
    }

    z_file->deep_ents.len32 = ent_i;

    mutex_unlock (z_file->deep_populate_ents_mutex);
} 
