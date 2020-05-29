// ------------------------------------------------------------------
//   zfile.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <errno.h>
#include <time.h>
#include <math.h>
#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "crypt.h"
#include "vblock.h"
#include "move_to_front.h"
#include "file.h"
#include "endianness.h"
#include "version.h"
#include "sections.h"
#include "compressor.h"
#include "piz.h"
#include "zip.h"
#include "license.h"
#include "arch.h"
#include "strings.h"
#include "dict_id.h"

bool is_v2_or_above=true, is_v3_or_above=true, is_v4_or_above=true, is_v5_or_above=true; // default to 'true' for ZIP, modifed when PIZ reads the genozip header

static const char *password_test_string = "WhenIThinkBackOnAllTheCrapIlearntInHighschool";

void zfile_show_header (const SectionHeader *header, VBlock *vb /* optional if output to buffer */)
{
    DictIdType dict_id = {0};
    char flags[10], ltype[10];
    bool has_flags = false;

    if (section_type_is_dictionary (header->section_type)) dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    
    else if (header->section_type == SEC_LOCAL || section_type_is_b250  (header->section_type)) {
        SectionHeaderCtx *header_ctx = (SectionHeaderCtx *)header;
        dict_id = header_ctx->dict_id;
        str_int (header_ctx->flags, flags);
        str_int (header_ctx->ltype, ltype);
        has_flags = true;
    }

    char str[1000];

    sprintf (str, "%-19s %*.*s %6s%-3s %11s%-3s alg=%u vb_i=%-3u sec_i=%-2u comp_offset=%-6u uncomp_len=%-6u comp_len=%-6u enc_len=%-6u magic=%8.8x\n",
             st_name(header->section_type), -DICT_ID_LEN, DICT_ID_LEN, dict_id.num ? dict_id_printable (dict_id).id : dict_id.id,
             has_flags ? "flags=" : "", has_flags ? flags : "", has_flags ? "ltype=" : "", has_flags ? ltype : "", 
             header->sec_compression_alg,
             BGEN32 (header->vblock_i), BGEN16 (header->section_i), 
             BGEN32 (header->compressed_offset), BGEN32 (header->data_uncompressed_len),
             BGEN32 (header->data_compressed_len), BGEN32 (header->data_encrypted_len), BGEN32 (header->magic));

    if (vb) {
        unsigned len = strlen (str);
        buf_alloc (vb, &vb->show_headers_buf, vb->show_headers_buf.len + len + 1, 2, "show_headers_buf", 0);
        strcpy (&vb->show_headers_buf.data[vb->show_headers_buf.len], str);
        vb->show_headers_buf.len += len;
    }
    else 
        printf ("%s", str);
}

static void zfile_show_b250_section (void *section_header_p, Buffer *b250_data)
{
    SectionHeaderCtx *header = (SectionHeaderCtx *)section_header_p;

    if (!flag_show_b250 && dict_id_printable (header->dict_id).num != dict_id_show_one_b250.num) return;

    static bool first = true;
    if (flag_show_b250 && first) { 
        fprintf (stderr, "Base-250 data for VB %u (result of '--show-b250'):\n", BGEN32 (header->h.vblock_i));
        first = false;
    }

    fprintf (stderr, "  %*.*s: ", -DICT_ID_LEN-1, DICT_ID_LEN, dict_id_printable (header->dict_id).id);

    const uint8_t *data  = FIRSTENT (const uint8_t, *b250_data);
    const uint8_t *after = AFTERENT (const uint8_t, *b250_data);

    while (data < after) {
        uint32_t word_index = base250_decode (&data);
        switch (word_index) {
            case WORD_INDEX_ONE_UP     : fprintf (stderr, "ONE_UP "); break;
            case WORD_INDEX_EMPTY_SF   : fprintf (stderr, "EMPTY "); break;
            case WORD_INDEX_MISSING_SF : fprintf (stderr, "MISSING "); break;
            default: fprintf (stderr, "%u ", word_index);
        }
    }
    fprintf (stderr, "\n");
}

// uncompressed a block and adds a \0 at its end. Returns the length of the uncompressed block, without the \0.
// when we get here, the header is already unencrypted zfile_`one_section
void zfile_uncompress_section (VBlock *vb,
                               void *section_header_p,
                               Buffer *uncompressed_data,
                               const char *uncompressed_data_buf_name,
                               SectionType expected_section_type) 
{
    START_TIMER;

    DictIdType dict_id = {0};
    if (section_type_is_dictionary (expected_section_type))
        dict_id = ((SectionHeaderDictionary *)section_header_p)->dict_id;
    else if (section_type_is_b250 (expected_section_type))
        dict_id = ((SectionHeaderCtx *)section_header_p)->dict_id;
    else if (expected_section_type == SEC_LOCAL)
        dict_id = ((SectionHeaderCtx *)section_header_p)->dict_id;

    if (piz_is_skip_section (vb, expected_section_type, dict_id)) return; // we skip some sections based on flags

    SectionHeader *section_header = (SectionHeader *)section_header_p;
    
    uint32_t compressed_offset     = BGEN32 (section_header->compressed_offset);
    uint32_t data_encrypted_len    = BGEN32 (section_header->data_encrypted_len);
    uint32_t data_compressed_len   = BGEN32 (section_header->data_compressed_len);
    uint32_t data_uncompressed_len = BGEN32 (section_header->data_uncompressed_len);
    uint32_t vblock_i              = BGEN32 (section_header->vblock_i);
    uint16_t section_i             = BGEN16 (section_header->section_i);
    CompressionAlg comp_alg        = (CompressionAlg)section_header->sec_compression_alg;

    // prior to version 5, the algorithm was hard coded, and this field was called "unused"
    if (!is_v5_or_above) {
        if (expected_section_type == SEC_HT_GTSHARK_DB_DB || expected_section_type == SEC_HT_GTSHARK_DB_GT)
            comp_alg = COMP_PLN;
        else
            comp_alg = COMP_BZ2;
    }

    // sanity checks
    ASSERT (section_header->section_type == expected_section_type, "Error in zfile_uncompress_section: expecting section type %s but seeing %s", st_name(expected_section_type), st_name(section_header->section_type));
    
    bool expecting_vb_i = !section_type_is_dictionary (expected_section_type) && expected_section_type != SEC_TXT_HEADER;
    ASSERT (vblock_i == vb->vblock_i || !expecting_vb_i, // dictionaries are uncompressed by the I/O thread with pseduo_vb (vb_i=0) 
             "Error: bad vblock_i: in file=%u in vb=%u", vblock_i, vb->vblock_i);

    // decrypt data (in-place) if needed
    if (data_encrypted_len) {
        if (is_v2_or_above)
            crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, vblock_i, section_header->section_type, false);
        else
            vcf_v1_crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, vblock_i, section_i);
    }

    if (data_uncompressed_len > 0) { // FORMAT, for example, can be missing in a sample-less file

        buf_alloc (vb, uncompressed_data, data_uncompressed_len + 1, 1.1, uncompressed_data_buf_name, 0); // +1 for \0
        uncompressed_data->len = data_uncompressed_len;

        comp_uncompress (vb, comp_alg, (char*)section_header + compressed_offset, data_compressed_len, uncompressed_data);
    }

    if (flag_show_b250 && section_type_is_b250 (expected_section_type)) 
        zfile_show_b250_section (section_header_p, uncompressed_data);
    
    if (vb) COPY_TIMER(vb->profile.zfile_uncompress_section);
}

// generate the metadata string eg "user@domain.com 2019-11-16 16:48:19"
static void zfile_get_metadata(char *metadata)
{
    time_t timer;
    char time_buf[20];
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", tm_info);

    const char *user = getenv ("USER");
    const char *host = getenv ("HOSTNAME");
    sprintf(metadata, "%.19s %.20s%s%.30s", time_buf,  // 71 chars + the string termintor \0 = 72
            user ? user : "",
            user && host ? "@" : "",
            host ? host : "");

    // sanity
    ASSERT0 (strlen (metadata) < FILE_METADATA_LEN, "Error: metadata too long");
}

// ZIP: called by compute threads, while holding the compress_dictionary_data_mutex
void zfile_compress_dictionary_data (VBlock *vb, MtfContext *ctx, 
                                     uint32_t num_words, const char *data, uint32_t num_chars)
{
    START_TIMER;

    SectionHeaderDictionary header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = SEC_DICT;
    header.h.data_uncompressed_len = BGEN32 (num_chars);
    header.h.compressed_offset     = BGEN32 (sizeof(SectionHeaderDictionary));
    header.h.sec_compression_alg  = COMP_BZ2;
    header.h.vblock_i              = BGEN32 (vb->vblock_i);
    header.h.section_i             = BGEN16 (vb->z_next_header_i++);
    header.num_snips               = BGEN32 (num_words);
    header.dict_id                 = ctx->dict_id;

    if (flag_show_dict) {
        fprintf (stderr, "%s (vb_i=%u, did=%u, num_snips=%u):\t", 
                 ctx->name, vb->vblock_i, ctx->did_i, num_words);
        str_print_null_seperated_data (data, num_chars, true);
    }
    
    if (dict_id_printable (ctx->dict_id).num == dict_id_show_one_dict.num)
        str_print_null_seperated_data (data, num_chars, false);

    z_file->dict_data.name  = "z_file->dict_data"; // comp_compress requires that it is set in advance
    comp_compress (vb, &z_file->dict_data, true, (SectionHeader*)&header, data, NULL);

    COPY_TIMER (vb->profile.zfile_compress_dictionary_data)    
//printf ("End compress dict vb_i=%u did_i=%u\n", vb->vblock_i, ctx->did_i);
}

void zfile_compress_b250_data (VBlock *vb, MtfContext *ctx, CompressionAlg comp_alg)
{
    SectionHeaderCtx header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = SEC_B250;
    header.h.data_uncompressed_len = BGEN32 (ctx->b250.len);
    header.h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCtx));
    header.h.sec_compression_alg   = comp_alg;
    header.h.vblock_i              = BGEN32 (vb->vblock_i);
    header.h.section_i             = BGEN16 (vb->z_next_header_i++);
    header.dict_id                 = ctx->dict_id;
    header.flags                   = ctx->flags;
    header.ltype              = ctx->ltype;
            
    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, ctx->b250.data, NULL);
}

static CompGetLineCallback *zfile_get_local_data_callback (DataType dt, DictIdType dict_id)
{
    static struct { DataType dt; uint64_t *dict_id_num; CompGetLineCallback *func; } callbacks[] = LOCAL_COMPRESSOR_CALLBACKS;

    for (unsigned i=0; i < sizeof(callbacks)/sizeof(callbacks[0]); i++)
        if (callbacks[i].dt == dt && *callbacks[i].dict_id_num == dict_id.num) 
            return callbacks[i].func;

    return NULL;
}

void zfile_compress_local_data (VBlock *vb, MtfContext *ctx)
{   
    SectionHeaderCtx header;
    memset (&header, 0, sizeof(header)); // safety

    unsigned local_len_multiplier  = ctx_lt_sizeof_one[ctx->ltype];

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = SEC_LOCAL;
    header.h.data_uncompressed_len = BGEN32 (ctx->local.len * local_len_multiplier); 
    header.h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCtx));
    header.h.sec_compression_alg   = (ctx->flags & CTX_FL_LOCAL_LZMA) ? COMP_LZMA : COMP_BZ2;
    header.h.vblock_i              = BGEN32 (vb->vblock_i);
    header.h.section_i             = BGEN16 (vb->z_next_header_i++);
    header.dict_id                 = ctx->dict_id;
    header.flags                   = ctx->flags;
    header.ltype                   = ctx->ltype;

    CompGetLineCallback *callback = zfile_get_local_data_callback (vb->data_type, ctx->dict_id);

    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, 
                   callback ? NULL : ctx->local.data, 
                   callback);
}

// compress section - two options for input data - 
// 1. contiguous data in section_data 
// 2. line by line data - by providing a callback + total_len
void zfile_compress_section_data_alg (VBlock *vb, SectionType section_type, 
                                      Buffer *section_data,          // option 1 - compress contiguous data
                                      CompGetLineCallback callback, uint32_t total_len, // option 2 - compress data one line at a time
                                      CompressionAlg comp_alg)
{
    SectionHeader header;
    memset (&header, 0, sizeof(header)); // safety

    header.magic                   = BGEN32 (GENOZIP_MAGIC);
    header.section_type            = section_type;
    header.data_uncompressed_len   = BGEN32 (section_data ? section_data->len : total_len);
    header.compressed_offset       = BGEN32 (sizeof(header));
    header.sec_compression_alg = comp_alg;
    header.vblock_i                = BGEN32 (vb->vblock_i);
    header.section_i               = BGEN16 (vb->z_next_header_i++);
    
    vb->z_data.name  = "z_data"; // comp_compress requires that these are pre-set
    vb->z_data.param = vb->vblock_i;
    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, 
                   section_data ? section_data->data : NULL, 
                   callback);
}

// reads exactly the length required, error otherwise. manages read buffers to optimize I/O performance.
// this doesn't make a big difference for SSD, but makes a huge difference for HD
// return true if data was read as requested, false if file has reached EOF and error otherwise
void *zfile_read_from_disk (VBlock *vb, Buffer *buf, unsigned len, bool fail_quietly_if_not_enough_data)
{
    START_TIMER;

    ASSERT0 (len, "Error: in zfile_read_from_disk, len is 0");

    void *start = (void *)&buf->data[buf->len];

    unsigned memcpyied = 0;
    unsigned len_save  = len;

    while (len) {

        if (z_file->z_next_read == z_file->z_last_read) { // nothing left in read_buffer - replenish it from disk

            if (z_file->z_last_read != READ_BUFFER_SIZE && !memcpyied) return NULL; // EOF reached last time, nothing more to read

            if (z_file->z_last_read != READ_BUFFER_SIZE && memcpyied && fail_quietly_if_not_enough_data) return NULL; // partial read = quiet error as requested
            
            ASSERT (z_file->z_last_read == READ_BUFFER_SIZE, 
                    "Error: end-of-file while reading %s, read %u bytes, but wanted to read %u bytes", 
                    z_name, memcpyied, len_save);

            z_file->z_last_read = fread (z_file->read_buffer, 1, READ_BUFFER_SIZE, (FILE *)z_file->file);
            z_file->z_next_read = 0;
            z_file->disk_so_far += z_file->z_last_read;

            ASSERT (z_file->z_last_read || !memcpyied, "Error: data requested could not be read bytes_so_far=%"PRIu64"", z_file->disk_so_far);
            if (!z_file->z_last_read) {
                COPY_TIMER (vb->profile.read);
                return NULL; // file is exhausted - nothing read
            }
        }

        unsigned memcpy_len = MIN (len, z_file->z_last_read - z_file->z_next_read);

        buf_add (buf, z_file->read_buffer + z_file->z_next_read, memcpy_len);
        len                 -= memcpy_len;
        z_file->z_next_read += memcpy_len;

        memcpyied += memcpy_len;
    }

    COPY_TIMER (vb->profile.read);

    return start;
}

// read section header - called from the I/O thread, but for a specific VB
// returns offset of header within data, EOF if end of file
int zfile_read_section (VBlock *vb, 
                        uint32_t original_vb_i, // the vblock_i used for compressing. this is part of the encryption key. dictionaries are compressed by the compute thread/vb, but uncompressed by the I/O thread (vb=0)
                        uint32_t sb_i,          // sample block number, NO_SB_I if this section type is not related to vcf samples
                        Buffer *data, const char *buf_name, // buffer to append 
                        unsigned header_size, SectionType expected_sec_type,
                        const SectionListEntry *sl)   // NULL for no seeking
{
    ASSERT (!sl || expected_sec_type == sl->section_type, "Error in zfile_read_section: expected_sec_type=%s but encountered sl->section_type=%s. vb_i=%u, sb_i=%d",
            st_name (expected_sec_type), st_name(sl->section_type), vb->vblock_i, sb_i);

    if (sl && piz_is_skip_section (vb, expected_sec_type, sl->dict_id)) return 0; // skip if this section is not needed according to flags

    if (sb_i != NO_SB_I && !vcf_is_sb_included(vb, sb_i)) return 0; // skip section if this sample block is excluded by --samples

    // note: the first section is always read by zfile_read_section(). if it is a v1, or
    // if it is not readable (perhaps v1 encrypted) - we assume its a v1 VCF header
    // in v2+ the first section is always an unencrypted SectionHeaderGenozipHeader
    ASSERT0 (is_v2_or_above || expected_sec_type==SEC_GENOZIP_HEADER, "Error: zfile_read_section cannot read v1 data");

    bool is_encrypted = (expected_sec_type != SEC_GENOZIP_HEADER) &&
                         crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    unsigned header_offset = data->len;
    buf_alloc (vb, data, header_offset + header_size, 2, buf_name, 1);
    
    // move the cursor to the section. file_seek is smart not to cause any overhead if no moving is needed
    if (sl) file_seek (z_file, sl->offset, SEEK_SET, false);

    SectionHeader *header = zfile_read_from_disk (vb, data, header_size, false); // note: header in file can be shorter than header_size if its an earlier version

    // case: this is a v5+ genozip header - read the extra field
    if (header->section_type == SEC_GENOZIP_HEADER &&
        ((v2v3v4_SectionHeaderGenozipHeader *)header)->genozip_version >= 5) {
        zfile_read_from_disk (vb, data, sizeof (SectionHeaderGenozipHeader) - sizeof (v2v3v4_SectionHeaderGenozipHeader), false);
        header_size += sizeof (Md5Hash);
    }

    // case: we're done! no more concatenated files
    if (!header && expected_sec_type == SEC_TXT_HEADER) return EOF; 

    ASSERT (header, "Error in zfile_read_section: Failed to read data from file %s while expecting section type %s: %s", 
            z_name, st_name(expected_sec_type), strerror (errno));
    
    if (flag_show_headers) zfile_show_header (header, NULL);

    bool is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;

    // decrypt header (note: except for SEC_GENOZIP_HEADER - this header is never encrypted)
    if (is_encrypted) {
        ASSERT (BGEN32 (header->magic) != GENOZIP_MAGIC, 
                "Error: password provided, but file %s is not encrypted", z_name);

        crypt_do (vb, (uint8_t*)header, header_size, original_vb_i, expected_sec_type, true); // negative section_i for a header
    
        is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC; // update after decryption
    }

    ASSERT (is_magical, "Error: corrupt data (magic is wrong) when attempting to read section %s of vblock_i=%u in file %s", 
            st_name (expected_sec_type), vb->vblock_i, z_name);

    unsigned compressed_offset   = BGEN32 (header->compressed_offset);
    ASSERT (compressed_offset, "Error: header.compressed_offset is 0 when reading section_type=%s", st_name(expected_sec_type));

    unsigned data_compressed_len = BGEN32 (header->data_compressed_len);
    unsigned data_encrypted_len  = BGEN32 (header->data_encrypted_len);

    unsigned data_len = MAX (data_compressed_len, data_encrypted_len);
    
    // check that we received the section type we expect, 
    ASSERT (expected_sec_type == header->section_type,
            "Error: Unexpected section type: expecting %s, found %s", st_name(expected_sec_type), st_name(header->section_type));

    ASSERT (compressed_offset == header_size,
            "Error: invalid header - expecting compressed_offset to be %u but found %u. section_type=%s", 
            header_size, compressed_offset, st_name(header->section_type));

    // allocate more memory for the rest of the header + data (note: after this realloc, header pointer is no longer valid)
    buf_alloc (vb, data, header_offset + compressed_offset + data_len, 2, "zfile_read_section", 2);
    header = (SectionHeader *)&data->data[header_offset]; // update after realloc

    // read section data
    if (data_len) {
        ASSERT (zfile_read_from_disk (vb, data, data_len, false), 
                "Error: failed to read section data, section_type=%s: %s", st_name(header->section_type), strerror (errno));
    }

    return header_offset;
}

void zfile_read_all_dictionaries (uint32_t last_vb_i /* 0 means all VBs */, ReadChromeType read_chrom)
{
    SectionListEntry *sl_ent = NULL; // NULL -> first call to this sections_get_next_dictionary() will reset cursor 

    mtf_initialize_primary_field_ctxs (z_file->contexts, z_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_dict_ids);

    while (sections_get_next_dictionary (&sl_ent)) {

        if (last_vb_i && sl_ent->vblock_i > last_vb_i) break;

        // cases where we can skip reading these dictionaries because we don't be using them
        bool is_chrom = (sl_ent->dict_id.num == dict_id_fields[DTFZ(chrom)]);
        if (read_chrom == DICTREAD_CHROM_ONLY  && !is_chrom) continue;
        if (read_chrom == DICTREAD_EXCEPT_CHROM && is_chrom) continue;

        if (piz_is_skip_sectionz (sl_ent->section_type, sl_ent->dict_id)) continue;
        
        zfile_read_section (evb, sl_ent->vblock_i, NO_SB_I, &evb->z_data, "z_data", sizeof(SectionHeaderDictionary), sl_ent->section_type, sl_ent);    

        // update dictionaries in z_file->contexts with dictionary data 
        mtf_integrate_dictionary_fragment (evb, evb->z_data.data);

        buf_free (&evb->z_data);
    }

    // output the dictionaries if we're asked to
    if (flag_show_dict || dict_id_show_one_dict.num) {
        for (uint32_t did_i=0; did_i < z_file->num_dict_ids; did_i++) {
            MtfContext *ctx = &z_file->contexts[did_i];

#define MAX_PRINTABLE_DICT_LEN 100000

            if (dict_id_printable (ctx->dict_id).num == dict_id_show_one_dict.num) 
                fprintf (stderr, "%.*s\t", (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), ctx->dict.data);
            
            if (flag_show_dict)
                fprintf (stderr, "%s (did_i=%u, num_snips=%u):\t%.*s\n", 
                         ctx->name, did_i, (uint32_t)ctx->word_list.len, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), ctx->dict.data);
        }
        fprintf (stderr, "\n");

        if (exe_type == EXE_GENOCAT) exit(0); // if this is genocat - we're done
    }
}

// returns the read data IF this is an invalid v2+ file, and hence might be v1
int16_t zfile_read_genozip_header (Md5Hash *digest) // out
{
    // read the footer from the end of the file
    file_seek (z_file, -sizeof(SectionFooterGenozipHeader), SEEK_END, false);

    SectionFooterGenozipHeader footer;
    int ret = fread (&footer, sizeof (footer), 1, (FILE *)z_file->file);
    ASSERTW (ret == 1, "Skipping empty file %s", z_name);
    if (!ret) return DT_NONE;
    
    // case: this is not a valid genozip v2+ file... maybe its v1?
    if (BGEN32 (footer.magic) != GENOZIP_MAGIC) {
        is_v2_or_above = is_v3_or_above = is_v4_or_above = is_v5_or_above = false;
        file_seek (z_file, 0, SEEK_SET, false);
        return DT_VCF_V1;
    }

    // read genozip header

    // note: for v1, we will use this function only for the very first VCF header (which will tell us this is v1)
    SectionListEntry dummy_sl = { .section_type = SEC_GENOZIP_HEADER,
                                  .offset       = BGEN64 (footer.genozip_header_offset) };

    // first, read the older section len, if its version 5+, read the extra field
    ret = zfile_read_section (evb, 0, NO_SB_I, &evb->z_data, "genozip_header", sizeof(v2v3v4_SectionHeaderGenozipHeader), 
                              SEC_GENOZIP_HEADER, &dummy_sl);

    ASSERT0 (ret != EOF, "Error: unexpected EOF when reading genozip header");
    
    SectionHeaderGenozipHeader *header = (SectionHeaderGenozipHeader *)evb->z_data.data;

    DataType data_type = (DataType)(BGEN16 (header->data_type)); 
    ASSERT ((unsigned)data_type < NUM_DATATYPES, "Error in zfile_read_genozip_header: unrecognized data_type=%d", data_type);
    z_file->data_type = data_type; // update in case type was not know from file extension
    if (txt_file) txt_file->data_type = data_type; // txt_file is still NULL in case of --split
    
    ASSERT (header->genozip_version <= GENOZIP_FILE_FORMAT_VERSION, 
            "Error: %s cannot be openned because it was compressed with a newer version of genozip (version %u.x.x) while the version you're running is older (version %s).\n"
            "You might want to consider upgrading genozip to the newest version.\n"
            "Note: genozip is backward-compatible: newer versions of genozip can always decompress files compressed with older versions",
            z_name, header->genozip_version, GENOZIP_CODE_VERSION);

    ASSERT (header->encryption_type != ENCRYPTION_TYPE_NONE || !crypt_have_password(), 
            "Error: password provided, but file %s is not encrypted", z_name);

    // get & test password, if file is encrypted
    if (header->encryption_type != ENCRYPTION_TYPE_NONE) {

        if (!crypt_have_password()) crypt_prompt_for_password();

        crypt_do (evb, header->password_test, sizeof(header->password_test), 0, SEC_NONE, true); // decrypt password test

        ASSERT (!memcmp (header->password_test, password_test_string, sizeof(header->password_test)),
                "Error: password is wrong for file %s", z_name);
    }

    global_vcf_num_samples    = BGEN32 (header->num_samples); // possibly 0, if genozip header was not rewritten. in this case, piz will get it from the first VCF header, but genols will show 0
    z_file->genozip_version   = header->genozip_version;
    z_file->num_components    = BGEN32 (header->num_components);
    *digest                   = header->md5_hash_concat; 

    // global bools to help testing
    is_v2_or_above = (z_file->genozip_version >= 2);
    is_v3_or_above = (z_file->genozip_version >= 3);
    is_v4_or_above = (z_file->genozip_version >= 4);
    is_v5_or_above = (z_file->genozip_version >= 5);
     
    zfile_uncompress_section (evb, header, &z_file->section_list_buf, "z_file->section_list_buf", SEC_GENOZIP_HEADER);
    z_file->section_list_buf.len /= sizeof (SectionListEntry); // fix len
    BGEN_sections_list();

    if (flag_show_gheader) {
        sections_show_gheader (header);
        if (exe_type == EXE_GENOCAT) exit(0); // in genocat, exit after showing the requested data
    }

    buf_free (&evb->z_data);

    return data_type;
}

void zfile_compress_genozip_header (const Md5Hash *single_component_md5)
{
    SectionHeaderGenozipHeader header;

    // start with just the fields needed by sections_add_to_list
    memset (&header, 0, sizeof(SectionHeaderGenozipHeader)); // safety
    header.h.section_type = SEC_GENOZIP_HEADER;

    // "manually" add the genozip section to the section list - normally it is added in comp_compress()
    // but in this case the genozip section containing the list will already be ready...
    sections_add_to_list (evb, &header.h);

    bool is_encrypted = crypt_have_password();

    uint32_t num_sections = z_file->section_list_buf.len;

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.compressed_offset     = BGEN32 (sizeof (SectionHeaderGenozipHeader));
    header.h.data_uncompressed_len = BGEN32 (z_file->section_list_buf.len * sizeof (SectionListEntry));
    header.h.sec_compression_alg   = COMP_BZ2;
    header.genozip_version         = GENOZIP_FILE_FORMAT_VERSION;
    header.data_type               = BGEN16 ((uint16_t)z_file->data_type);
    header.encryption_type         = is_encrypted ? ENCRYPTION_TYPE_AES256 : ENCRYPTION_TYPE_NONE;
    header.uncompressed_data_size  = BGEN64 (z_file->txt_data_so_far_concat);
    header.num_samples             = BGEN32 (global_vcf_num_samples);
    header.num_items_concat        = BGEN64 (z_file->num_lines);
    header.num_sections            = BGEN32 (num_sections); 
    header.num_components          = BGEN32 (z_file->num_txt_components_so_far);

    uint32_t license_num_bgen = BGEN32 (license_get());
    md5_do (&license_num_bgen, sizeof (int32_t), &header.license_hash);

    if (flag_md5) {
        if (flag_concat) {
            md5_finalize (&z_file->md5_ctx_concat, &header.md5_hash_concat);
            if (flag_md5 && z_file->num_txt_components_so_far > 1 && !flag_quiet) 
                fprintf (stderr, "Concatenated VCF MD5 = %s\n", md5_display (&header.md5_hash_concat, false));
        } 
        else 
            header.md5_hash_concat = *single_component_md5; // if not in concat mode - just copy the md5 of the single file
    }

    zfile_get_metadata (header.created);

    if (is_encrypted) {
        memcpy (header.password_test, password_test_string, sizeof(header.password_test));
        crypt_do (evb, header.password_test, sizeof(header.password_test), 0, SEC_NONE, true);
    }

    Buffer *z_data = &evb->z_data;

    uint64_t genozip_header_offset = z_file->disk_so_far + z_data->len; // capture before comp_compress that increases len

    // prepare section list for disk - Big Endian and length in bytes
    BGEN_sections_list();
    z_file->section_list_buf.len *= sizeof (SectionListEntry); // change to counting bytes

    // compress section into z_data - to be eventually written to disk by the I/O thread
    evb->z_data.name  = "z_data"; // comp_compress requires that these are pre-set
    comp_compress (evb, z_data, true, &header.h, z_file->section_list_buf.data, NULL);

    // restore
    z_file->section_list_buf.len /= sizeof (SectionListEntry); 
    BGEN_sections_list();

    if (flag_show_gheader) sections_show_gheader (&header); 

    // add a footer to this section - this footer appears AFTER the genozip header data, 
    // facilitating reading the genozip header in reverse from the end of the file
    SectionFooterGenozipHeader footer;
    footer.magic                 = BGEN32 (GENOZIP_MAGIC);
    footer.genozip_header_offset = BGEN64 (genozip_header_offset);

    buf_alloc (evb, z_data, z_data->len + sizeof(SectionFooterGenozipHeader), 1.5, "z_data", 0);
    memcpy (&z_data->data[z_data->len], &footer, sizeof(SectionFooterGenozipHeader));
    z_data->len += sizeof(SectionFooterGenozipHeader);
}

// reads the the genozip header section's header from a GENOZIP file - used by main_list, returns true if successful
bool zfile_get_genozip_header (uint64_t *uncompressed_data_size,
                               uint32_t *num_samples,
                               uint64_t *num_items_concat,
                               Md5Hash  *md5_hash_concat,
                               char *created, unsigned created_len /* caller allocates space */,
                               Md5Hash  *license_hash)
{
    // read the footer from the end of the file
    if (!file_seek (z_file, -sizeof(SectionFooterGenozipHeader), SEEK_END, true))
        return false;

    SectionFooterGenozipHeader footer;
    int ret = fread (&footer, sizeof (footer), 1, (FILE *)z_file->file);
    ASSERTW (ret == 1, "Skipping empty file %s", z_name);    
    if (!ret) return false; // empty file / cannot read
    
    // case: this is not a valid genozip v2+ file... maybe its v1?
    if (BGEN32 (footer.magic) != GENOZIP_MAGIC) {
        file_seek (z_file, 0, SEEK_SET, false);
        return vcf_v1_header_get_vcf_header (uncompressed_data_size, num_samples, num_items_concat, 
                                             md5_hash_concat, created, created_len);
    }

    // read genozip header
    uint64_t genozip_header_offset = BGEN64 (footer.genozip_header_offset);
    if (!file_seek (z_file, genozip_header_offset, SEEK_SET, true))
        return false;

    SectionHeaderGenozipHeader header;
    int bytes = fread ((char*)&header, 1, sizeof(SectionHeaderGenozipHeader), (FILE *)z_file->file);
    if (bytes < sizeof(SectionHeaderGenozipHeader)) return false;

    ASSERTW (BGEN32 (header.h.magic) == GENOZIP_MAGIC, "Error reading %s: corrupt data", z_name);
    if (BGEN32 (header.h.magic) != GENOZIP_MAGIC) return false;

    *uncompressed_data_size = BGEN64 (header.uncompressed_data_size);
    *num_samples            = BGEN32 (header.num_samples);
    *num_items_concat       = BGEN64 (header.num_items_concat);
    *md5_hash_concat        = header.md5_hash_concat;

    if (header.genozip_version >= 5) // is_v5_or_above is for genols
        *license_hash = header.license_hash;
    else
        license_hash->ulls[0] = license_hash->ulls[1] = 0;

    memcpy (created, header.created, MIN (FILE_METADATA_LEN, created_len));

    return true;
}

// ZIP
void zfile_write_txt_header (Buffer *txt_header_text, bool is_first_txt)
{
    SectionHeaderTxtHeader header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = SEC_TXT_HEADER;
    header.h.data_uncompressed_len = BGEN32 (txt_header_text->len);
    header.h.compressed_offset     = BGEN32 (sizeof (SectionHeaderTxtHeader));
    header.h.sec_compression_alg  = COMP_BZ2;
    header.num_samples             = BGEN32 (global_vcf_num_samples);
    header.num_lines               = NUM_LINES_UNKNOWN; 
    header.compression_type        = (uint8_t)txt_file->comp_alg; 

    file_basename (txt_file->name, false, "(stdin)", header.txt_filename, TXT_FILENAME_LEN);
    header.txt_filename[strlen(header.txt_filename)- (strlen(file_exts[txt_file->type])-4)] = '\0'; // remove the .gz/.bgz/.bz2
    
    static Buffer txt_header_buf = EMPTY_BUFFER;

    buf_alloc (evb, &txt_header_buf, sizeof (SectionHeaderTxtHeader) + txt_header_text->len / 3, // generous guess of compressed size
            1, "txt_header_buf", 0); 

    comp_compress (evb, &txt_header_buf, true, (SectionHeader*)&header, 
                   txt_header_text->len ? txt_header_text->data : NULL, // actual header may be missing (eg in SAM it is permitted to not have a header)
                   NULL);

    file_write (z_file, txt_header_buf.data, txt_header_buf.len);

    z_file->disk_so_far            += txt_header_buf.len;   // length of GENOZIP data writen to disk
    z_file->txt_data_so_far_single += txt_header_text->len; // length of the original VCF header
    z_file->txt_data_so_far_concat += txt_header_text->len;

    buf_free (&txt_header_buf); 

    // copy it to z_file - we might need to update it at the very end in zfile_update_txt_header_section_header()
    if (is_first_txt)
        memcpy (&z_file->txt_header_first, &header, sizeof (header));

    memcpy (&z_file->txt_header_single, &header, sizeof (header));
}

// updating the VCF bytes of a GENOZIP file. If we're compressing a simple VCF file, we will know
// the bytes upfront, but if we're concatenating or compressing a VCF.GZ, we will need to update it
// when we're done. num_lines can only be known after we're done with this VCF component.
// if we cannot update the header - that's fine, these fields are only used for the progress indicator on --list
bool zfile_update_txt_header_section_header (uint64_t pos_of_current_vcf_header, uint32_t max_lines_per_vb,
                                             Md5Hash *md5 /* out */)
{
    // rewind to beginning of current (latest) vcf header - nothing to do if we can't
    if (!file_seek (z_file, pos_of_current_vcf_header, SEEK_SET, true)) return false;

    unsigned len = crypt_padded_len (sizeof (SectionHeaderTxtHeader));

    // update the header of the single (current) vcf. 
    SectionHeaderTxtHeader *curr_header = &z_file->txt_header_single;
    curr_header->txt_data_size    = BGEN64 (txt_file->txt_data_size_single);
    curr_header->num_lines        = BGEN64 (txt_file->num_lines);
    curr_header->max_lines_per_vb = BGEN32 (max_lines_per_vb);

    md5_finalize (&z_file->md5_ctx_single, &curr_header->md5_hash_single);
    *md5 = curr_header->md5_hash_single;
    if (flag_md5 && !flag_quiet) 
        fprintf (stderr, "MD5 = %s\n", md5_display (&curr_header->md5_hash_single, false));

    if (pos_of_current_vcf_header == 0) 
        z_file->txt_header_first.md5_hash_single = curr_header->md5_hash_single; // first vcf - update the stored header 

    // encrypt if needed
    if (crypt_have_password()) 
        crypt_do (evb, (uint8_t *)curr_header, len, 0, curr_header->h.section_type, true);

    file_write (z_file, curr_header, len);
    fflush ((FILE*)z_file->file); // its not clear why, but without this fflush the bytes immediately after the first header get corrupted (at least on Windows with gcc)
    
    file_seek (z_file, 0, SEEK_END, false); // return to the end of the file

    return true; // success
}

// ZIP compute thread - called from sam_zip_compress_one_vb()
void zfile_compress_generic_vb_header (VBlock *vb)
{
    uint32_t sizeof_header = sizeof (SectionHeaderVbHeader);

    SectionHeaderVbHeader vb_header;
    memset (&vb_header, 0, sizeof(vb_header)); // safety
    
    vb_header.h.magic               = BGEN32 (GENOZIP_MAGIC);
    vb_header.h.section_type        = SEC_VB_HEADER;
    vb_header.h.compressed_offset   = BGEN32 (sizeof_header);
    vb_header.h.vblock_i            = BGEN32 (vb->vblock_i);
    vb_header.h.section_i           = BGEN16 (vb->z_next_header_i++); // always 0
    vb_header.h.sec_compression_alg = COMP_PLN;
    vb_header.num_lines             = BGEN32 ((uint32_t)vb->lines.len);
    vb_header.vb_data_size          = BGEN32 (vb->vb_data_size);
    vb_header.longest_line_len      = BGEN32 (vb->longest_line_len);

    // copy section header into z_data - to be eventually written to disk by the I/O thread. this section doesn't have data.
    vb->z_data.name  = "z_data"; // this is the first allocation of z_data - comp_compress requires that it is pre-named
    vb->z_data.param = vb->vblock_i;
    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&vb_header, NULL, NULL);
}

// ZIP only: called by the I/O thread in the sequential order of VBs: updating of the already compressed
// variant data section (compressed by the compute thread in zfile_compress_generic_vb_header) just before writing it to disk
// note: this updates the z_data in memory (not on disk)
void zfile_update_compressed_vb_header (VBlock *vb, uint32_t txt_first_line_i)
{
    SectionHeaderVbHeader *vb_header = (SectionHeaderVbHeader *)vb->z_data.data;
    vb_header->z_data_bytes = BGEN32 ((uint32_t)vb->z_data.len);
    vb_header->first_line   = BGEN32 (txt_first_line_i);

    if (flag_show_vblocks) {
        fprintf (stderr, "vb_i=%u first_line=%u num_lines=%u uncomprssed=%u compressed=%u longest_line_len=%u\n",
                 vb->vblock_i, txt_first_line_i, BGEN32 (vb_header->num_lines), 
                 BGEN32 (vb_header->vb_data_size), BGEN32 (vb_header->z_data_bytes), 
                 BGEN32 (vb_header->longest_line_len));
    }
    // now we can finally encrypt the header - if needed
    if (crypt_have_password())  
        crypt_do (vb, (uint8_t*)vb_header, BGEN32 (vb_header->h.compressed_offset),
                  BGEN32 (vb_header->h.vblock_i), vb_header->h.section_type, true);
}
