// ------------------------------------------------------------------
//   cram.c
//   Copyright (C) 2023-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "igzip/igzip_lib.h"

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

typedef packed_enum { // 1 byte 
    CRAM_CODEC_NONE=0, CRAM_CODEC_GZIP=1, CRAM_CODEC_BZ2=2, CRAM_CODEC_LZMA=3,
    CRAM_CODEC_RANS4x8=4, CRAM_CODEC_RANS4x16=5, CRAM_CODEC_ARITH=6, 
    CRAM_CODEC_FQZCOMP=7, CRAM_CODEC_TOKENIZER=8 
} CramCodec;

typedef packed_enum { // 1 byte 
    CRAM_FILE_HEADER=0, CRAM_COMPRESSION_HEADER=1, CRAM_SLICE_HEADER=2, EXTERNAL_DATA=4, CORE_DATA=5
} CramBlockContentType;

typedef struct {
    CramCodec codec; 
    CramBlockContentType content_type;
    int32_t content_id;
    int32_t compressed_size;
    int32_t uncompressed_size;
} CramBlockHeader;

static bool get_uint8 (uint8_t *value)
{
    if (evb->scratch.next + 1 > evb->scratch.len) return false;

    *value = *B8(evb->scratch, evb->scratch.next++);
    return true;
}

static bool get_int32 (int32_t *value)
{
    if (evb->scratch.next + 4 > evb->scratch.len) return false;

    memcpy (value, Bc(evb->scratch, evb->scratch.next), 4);
    *value = LTEN32 (*value);

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

// populate evb->txt_data, evb->lines.len, file->header_size, file->txt_data_so_far_single
static bool cram_read_sam_header (FileP file, bool read_sam_header)
{
    uint64_t before_sam_header = evb->scratch.next;
    int32_t sam_header_len32 = 0;

    // note: the Header container may contain multiple blocks, but per the spec the entire
    // SAM header is in the first block, and the remaining blocks are empty
    CramBlockHeader b;
    if (!get_uint8 (&b.codec) || (b.codec != CRAM_CODEC_NONE && b.codec != CRAM_CODEC_GZIP)) return false;
    if (!get_uint8 (&b.content_type) || b.content_type != CRAM_FILE_HEADER) return false;
    if (!get_itf8 (&b.content_id) || b.content_id != 0) return false;
    if (!get_itf8 (&b.compressed_size))   return false;
    if (!get_itf8 (&b.uncompressed_size) || b.uncompressed_size < 4) return false;

    uint32_t remaining = evb->scratch.len - evb->scratch.next;

    if (b.codec == CRAM_CODEC_NONE) {
        if (!get_int32 (&sam_header_len32) ||
             sam_header_len32 > remaining) // we haven't read the entire SAM header from disk
            return false;

        buf_alloc_exact (evb, evb->txt_data, sam_header_len32 + 1, char, "txt_data");
        memcpy (B1STc(evb->txt_data), Bc(evb->scratch, evb->scratch.next), sam_header_len32);
    }

    else { // CRAM_CODEC_GZIP
        if (b.compressed_size > remaining) // we haven't read the entire SAM header from disk
            return false;

        struct inflate_state state = {};
        isal_inflate_init (&state);

        state.crc_flag  = ISAL_GZIP;
        state.next_in   = B8(evb->scratch, evb->scratch.next);
        state.avail_in  = b.compressed_size;
        state.next_out  = (uint8_t *)&sam_header_len32;
        state.avail_out = sizeof (sam_header_len32);

        // read length of header
        int ret = isal_inflate (&state); 
        ASSERT (ret == ISAL_DECOMP_OK, "%s: Failed to read SAM header #%d (ret=%d)", file->name, 1, ret);

        buf_alloc_exact (evb, evb->txt_data, sam_header_len32 + 1, char, "txt_data");

        state.next_out  = B1ST8(evb->txt_data);
        state.avail_out = sam_header_len32;

        ret = isal_inflate (&state); 
        ASSERT (ret == ISAL_DECOMP_OK, "%s: Failed to read SAM header #%d (ret=%d)", file->name, 2, ret);
    }

    file->header_size = file->txt_data_so_far_single = sam_header_len32;

    ASSERT (*Bc(evb->txt_data, evb->txt_data.len32-2) == '\n', "%s: expecting SAM header to end with a \\n", file->name);

    // count lines    
    evb->lines.len32 = str_count_char (Bc(evb->txt_data, 8), evb->txt_data.len32 - 8, '\n');
    
    // verify that every header contig appears in the reference too, with same name and length (alt contigs not supported for CRAM)
    *BLSTc(evb->txt_data) = 0;
    sam_header_zip_inspect_SQ_lines_in_cram (file->basename);
    buf_free (evb->txt_data);
    
    evb->scratch.next = before_sam_header; // rollback
    return true;
}

// populates a CramContainerHeader and skips to end of container
static bool cram_get_container_header (FileP file, CramContainerHeader *h, bool read_sam_header)
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

    // case: get textual length of SAM header 
    if (read_sam_header && 
        !cram_read_sam_header (file, true)) return false;

    // skip blocks    
    evb->scratch.next += h->length; // this might cause next to be more than len

    h->length = evb->scratch.next - before; // update length to include header size

    return true; // true even if the data blocks are not included in scratch, as we don't need them
}

static bool cram_inspect_file_definition_data (FileP file)
{
    ARRAY(uint8_t, c, evb->scratch);

    // case: this is actually a GZ file (possibly BAM)
    if (c_len >= 2 && c[0] == 0x1f && c[1] == 0x8b) {
        file->codec = file->source_codec = CODEC_GZ;
        file->data_type = DT_GNRIC; // generic_is_header_done will figure out the true data type
        file->type = GNRIC_GZ;
        return false;
    }

    // case: this is not CRAM, but not a GZ file (perhaps SAM) 
    else if (c_len < 26 || memcmp (c, CRAM_MAGIC, STRLEN(CRAM_MAGIC))) {
        file->codec = file->source_codec = CODEC_NONE;
        file->data_type = DT_GNRIC; 
        file->type = GNRIC;
        return false;
    }

    // case: CRAM of a version other than 3 - we don't know how to parse it here, but samtools might work and hence we will be able to compress it
    if (c[4] != 3) return false;

    evb->scratch.next = 26; // past file definition
    return true;
}

// read and inspect the first (CRAM_INSPECTION_SIZE bytes) of the CRAM file: 
// 1. File Definition data (possibly changing file->codec/source_code/data_type/type to Generic)
// 2. The Header Container which contains the SAM header: populate evb->txt_data, evb->lines.len, file->header_size, file->txt_data_so_far_single
// 3. Some Data Containers which contain CRAM data: populate file->est_num_lines
void cram_inspect_file (FileP file)
{
    START_TIMER;

    #define CRAM_INSPECTION_SIZE (8 MB)

    FILE *fp = fopen (file->name, "rb");
    if (!fp) return;

    ASSERTNOTINUSE (evb->scratch);
    buf_alloc_exact (evb, evb->scratch, CRAM_INSPECTION_SIZE, char, "scratch"); 

    evb->scratch.len32 = fread (B1STc(evb->scratch), 1, evb->scratch.len32, fp);
    ASSERT (!ferror (fp), "Failed to read file %s", file->name);

    fclose (fp);

    // file definition
    if (!cram_inspect_file_definition_data (file)) goto done;

    // the Header container (containing the SAM header)
    CramContainerHeader h;
    if (!cram_get_container_header (file, &h, true)) goto done; // first block is the SAM header
    
    int64_t header_container_size = evb->scratch.next;   // including file definition header

    // Data containers
    int64_t n_records=0, disk_consumed=0;
    while (cram_get_container_header (file, &h, false)) {
        n_records     += h.n_records; // records are alignments
        disk_consumed += h.length;
    }

    if (!n_records) goto done; // not even one full data container - we can't set est_num_lines

    int64_t file_size = file_get_size (file->name);

    double one_record_size = (double)disk_consumed / (double)n_records; // this also allocates part of the container header to each record
    
    file->est_num_lines = (double)(file_size - header_container_size) / one_record_size;

done:
    buf_free (evb->scratch);
    COPY_TIMER_EVB (cram_inspect_file);
}

// returns the -T (reference) option for CRAM, derived from the genozip reference name
StrTextSuperLong cram_get_samtools_option_T (Reference ref)
{
    if (!ref_is_loaded (ref)) return (StrTextSuperLong){};

    StrTextSuperLong samtools_T_option;
    uint32_t samtools_T_option_len = 0;
    rom ref_filename = ref_get_filename (ref);
    rom ref_fasta_name = ref_get_fasta_name (ref);

    ASSINP0 (ref_filename, "when compressing a CRAM file, --reference or --REFERENCE must be specified");

    // cases where the FASTA name is not in SEC_GENOZIP_HEADER
    ASSINP (ref_fasta_name, "Genozip limitation: Reference file %s cannot be used to compress CRAM files because it was created by piping a fasta file from from a url or stdin, or because the name of the fasta file exceeds %u characters",
            ref_filename, REF_FILENAME_LEN-1);

    int samtools_T_option_size = MAX_(strlen (ref_fasta_name), strlen (ref_filename)) + 10;
    ASSERT (samtools_T_option_size < sizeof (StrTextSuperLong), "samtools_T_option_size=%u too large: ref_fasta_name=\"%s\" ref_filename=\"%s\"", 
            samtools_T_option_size, ref_fasta_name, ref_filename);

    // case: fasta file is in its original location
    if (file_exists (ref_fasta_name)) 
        SNPRINTF (samtools_T_option, "-T%s", ref_fasta_name);

    // try: fasta file is in directory of reference file
    else {
        rom slash = strrchr (ref_fasta_name, '/');
        if (!slash) slash = strrchr (ref_fasta_name, '\\'); 
        rom basename = slash ? slash+1 : ref_fasta_name;

        slash = strrchr (ref_filename, '/');
        if (!slash) slash = strrchr (ref_filename, '\\'); 
        unsigned dirname_len = slash ? ((slash+1) - ref_filename) : 0;

        SNPRINTF (samtools_T_option, "-T%.*s%s", dirname_len, ref_filename, basename);

        ASSINP (file_exists (&samtools_T_option.s[2]), 
                "Searching of the fasta file used to create %s. It was not found in %s or %s. Note: it is needed for passing to samtools as a reference file (-T option) for reading the CRAM file", 
                ref_filename, ref_fasta_name, &samtools_T_option.s[2]);        
    }

    return samtools_T_option;
}
