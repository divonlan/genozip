// ------------------------------------------------------------------
//   cram.c
//   Copyright (C) 2023-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "buffer.h"
#include "sam_private.h"
#include "bits.h"

// see: https://samtools.github.io/hts-specs/CRAMv3.pdf
typedef struct {
    int32_t length;    // including header
    int32_t ref_seq_id;
    int32_t start_pos;
    int32_t aln_span;
    int32_t n_records;
    int64_t record_counter;
    int64_t bases;
    int32_t n_blocks;
    int32_t n_landmarks;
    int32_t crc32; 
} CramContainerHeader; // note: excluding the landmark array

static bool get_int32 (int32_t *value)
{
    if (evb->scratch.next + 4 > evb->scratch.len) return false;

    memcpy (value, Bc(evb->scratch, evb->scratch.next), 4);
    
    evb->scratch.next += 4;
    return true;
}

static bool get_itf8 (int32_t *value)
{
    if (evb->scratch.next + 5 > evb->scratch.len) return false;

    uint8_t *b = B8(evb->scratch, evb->scratch.next);
    
    int n_bytes = MIN_(4, __builtin_clz ((uint32_t)~b[0] << 24)); // 0 to 4: number of leading 1s followed by a 0 (except if 4 1s - no following 0)
    
    uint32_t v = b[0] & bitmask8(7 - MIN_(n_bytes,3)); // note: a bitmask is a value with 1s in the N lower bits and 0 in the higher bits

    for (int i=1; i <= n_bytes; i++)
        v = (v << 8) | b[i];

    if (value) *value = (int32_t)v;

    evb->scratch.next += 1 + n_bytes;
    return true;
}

static bool get_ltf8 (int64_t *value)
{
    if (evb->scratch.next + 9 > evb->scratch.len) return false;

    uint8_t *b = B8(evb->scratch, evb->scratch.next);
    
    int n_bytes = __builtin_clz ((uint32_t)~b[0] << 24); // 0 to 8: number of leading 1s followed by a 0 (except if 8 1s - no following 0)
    
    uint64_t v = (n_bytes < 7) ? (b[0] & bitmask8(7 - n_bytes)) : 0;

    for (int i=1; i <= n_bytes; i++)
        v = (v << 8) | b[i];

    if (value) *value = (int64_t)v;

    evb->scratch.next += 1 + n_bytes;
    return true;
}

// populates a CramContainerHeader and skips to end of container
static bool cram_get_container_header (CramContainerHeader *h)
{
    uint64_t before = evb->scratch.next;

    if (!get_int32 (&h->length))        return false;
    if (!get_itf8 (&h->ref_seq_id))     return false;
    if (!get_itf8 (&h->start_pos))      return false;
    if (!get_itf8 (&h->aln_span))       return false;
    if (!get_itf8 (&h->n_records))      return false;
    if (!get_ltf8 (&h->record_counter)) return false;
    if (!get_ltf8 (&h->bases))          return false;
    if (!get_itf8 (&h->n_blocks))       return false;
    if (!get_itf8 (&h->n_landmarks))    return false;

    // skip landmarks
    for (int i=0; i < h->n_landmarks; i++)
        if (!get_itf8 (NULL))           return false;

    if (!get_int32 (&h->crc32))         return false;

    // skip blocks    
    evb->scratch.next += h->length; // this might cause next to be more than len

    h->length = evb->scratch.next - before; // update length to include header size

    return true; // true even if the data blocks are not included in scratch, as we don't need them
}

static bool cram_verify_file_definition (void)
{
    if (evb->scratch.len < 26) return false;

    rom c = B1STc (evb->scratch);
    evb->scratch.next = 26;

    return c[0]=='C' && c[1]=='R' && c[2]=='A' && c[3]=='M' && c[4]==3; // CRAM version 3
}

uint64_t cram_estimated_num_alignments (rom filename)
{
    FILE *fp = fopen (filename, "rb");
    if (!fp) return 0;

    ASSERTNOTINUSE (evb->scratch);
    buf_alloc_exact (evb, evb->scratch, 8 MB, char, "scratch"); 

    int64_t est_num_records = 0; // return 0 if failed

    evb->scratch.len32 = fread (B1STc(evb->scratch), 1, evb->scratch.len32, fp);

    if (!cram_verify_file_definition()) goto done;

    // parse the SAM Header container
    CramContainerHeader h;
    if (!cram_get_container_header (&h)) goto done; // cram header
    
    int64_t cram_header_size = evb->scratch.next;   // including file definition header
    int64_t n_records=0, disk_consumed=0;

    while (cram_get_container_header (&h)) {
        n_records     += h.n_records; // records are alignments
        disk_consumed += h.length;
    }

    if (!n_records) goto done;

    int64_t file_size = file_get_size (filename);

    double one_record_size = (double)disk_consumed / (double)n_records; // this also allocates part of the container header to each record
    
    est_num_records = (double)(file_size - cram_header_size) / one_record_size;

done:
    fclose (fp);
    buf_free (evb->scratch);

    return est_num_records;
}
