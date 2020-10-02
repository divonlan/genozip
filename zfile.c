// ------------------------------------------------------------------
//   zfile.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <errno.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "crypt.h"
#include "vblock.h"
#include "context.h"
#include "file.h"
#include "endianness.h"
#include "version.h"
#include "sections.h"
#include "compressor.h"
#include "piz.h"
#include "zip.h"
#include "license.h"
#include "mutex.h"
#include "strings.h"
#include "dict_id.h"
#include "reference.h"

bool is_v6_or_above=true;

static const char *password_test_string = "WhenIThinkBackOnAllTheCrapIlearntInHighschool";

void zfile_show_header (const SectionHeader *header, VBlock *vb /* optional if output to buffer */)
{
    static const char *comp_names[NUM_COMPRESSION_ALGS] = CODEC_NAMES;

    if (flag_reading_reference) return; // don't show headers of reference file
    
    DictId dict_id = {0};
    char flags[10] = "", param[10] = "";
    const char *ltype="";
    bool has_ltype = false;

    if (header->section_type == SEC_DICT) 
        dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    
    else if (header->section_type == SEC_LOCAL || header->section_type == SEC_B250) {
        SectionHeaderCtx *header_ctx = (SectionHeaderCtx *)header;
        dict_id = header_ctx->dict_id;
        str_int (header_ctx->param, param);
        ltype = lt_desc[header_ctx->ltype].name;
        has_ltype = true;
    }

    if (header->flags) 
        str_int (header->flags, flags);

    char str[1000];

    sprintf (str, "%-19s %*.*s %6s%-3s %6s%-3s %7s%-3s codec=%-4.4s vb_i=%-3u comp_offset=%-6u uncomp_len=%-7u comp_len=%-6u enc_len=%-6u magic=%8.8x\n",
             st_name(header->section_type), -DICT_ID_LEN, DICT_ID_LEN, dict_id.num ? dict_id_printable (dict_id).id : dict_id.id,
             header->flags ? "flags=" : "", flags, has_ltype ? "ltype=" : "", ltype, has_ltype ? "param=" : "", param, 
             comp_names[header->codec],
             BGEN32 (header->vblock_i), BGEN32 (header->compressed_offset), BGEN32 (header->data_uncompressed_len),
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

static void zfile_show_b250_section (void *section_header_p, const Buffer *b250_data)
{
    MUTEX (show_b250_mutex); // protect so compute thread's outputs don't get mix

    SectionHeaderCtx *header = (SectionHeaderCtx *)section_header_p;

    if (!flag_show_b250 && dict_id_printable (header->dict_id).num != dict_id_show_one_b250.num) return;

    mutex_initialize (show_b250_mutex); // possible unlikely race condition on initializing - good enough for debugging purposes
    mutex_lock (show_b250_mutex);

    fprintf (stderr, "vb_i=%u %*.*s: ", BGEN32 (header->h.vblock_i), -DICT_ID_LEN-1, DICT_ID_LEN, dict_id_printable (header->dict_id).id);

    const uint8_t *data  = FIRSTENT (const uint8_t, *b250_data);
    const uint8_t *after = AFTERENT (const uint8_t, *b250_data);

    while (data < after) {
        WordIndex word_index = base250_decode (&data);
        switch (word_index) {
            case WORD_INDEX_ONE_UP     : fprintf (stderr, "ONE_UP "); break;
            case WORD_INDEX_EMPTY_SF   : fprintf (stderr, "EMPTY "); break;
            case WORD_INDEX_MISSING_SF : fprintf (stderr, "MISSING "); break;
            default: fprintf (stderr, "%u ", word_index);
        }
    }
    fprintf (stderr, "\n");

    mutex_unlock (show_b250_mutex);
}

// uncompressed a block and adds a \0 at its end. Returns the length of the uncompressed block, without the \0.
// when we get here, the header is already unencrypted zfile_one_section
void zfile_uncompress_section (VBlock *vb,
                               void *section_header_p,
                               void *uncompressed_data, // Buffer * or char *
                               const char *uncompressed_data_buf_name, // a name if Buffer, NULL if char *
                               SectionType expected_section_type) 
{
    START_TIMER;

    DictId dict_id = {0};
    if (expected_section_type      == SEC_DICT)
        dict_id = ((SectionHeaderDictionary *)section_header_p)->dict_id;
    else if (expected_section_type == SEC_B250)
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
    Codec codec        = (Codec)section_header->codec;

    // sanity checks
    ASSERT (section_header->section_type == expected_section_type, "Error in zfile_uncompress_section: expecting section type %s but seeing %s", st_name(expected_section_type), st_name(section_header->section_type));
    
    bool expecting_vb_i = expected_section_type != SEC_DICT && 
                          expected_section_type != SEC_TXT_HEADER && expected_section_type != SEC_REFERENCE && expected_section_type != SEC_REF_IS_SET;
                          
    ASSERT (vblock_i == vb->vblock_i || !expecting_vb_i, // dictionaries are uncompressed by the I/O thread with pseduo_vb (vb_i=0) 
             "Error in zfile_uncompress_section: bad vblock_i: vblock_i in file=%u but expecting it to be %u (section_type=%s)", 
             vblock_i, vb->vblock_i, st_name (expected_section_type));

    // decrypt data (in-place) if needed
    if (data_encrypted_len) 
        crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, vblock_i, section_header->section_type, false);

    BufferP uncompressed_buf = NULL;
    if (data_uncompressed_len > 0) { // FORMAT, for example, can be missing in a sample-less file

        if (uncompressed_data_buf_name) { // we're decompressing into a buffer
            uncompressed_buf = (BufferP)uncompressed_data;
            buf_alloc (vb, uncompressed_buf, data_uncompressed_len + sizeof (uint64_t), 1.1, uncompressed_data_buf_name, 0); // add a 64b word for safety in case this buffer will be converted to a bitarray later
            uncompressed_buf->len = data_uncompressed_len;
            uncompressed_data = uncompressed_buf->data;
        }

        comp_uncompress (vb, codec, (char*)section_header + compressed_offset, data_compressed_len, 
                         uncompressed_data, data_uncompressed_len);
    }
 
    if (uncompressed_buf && flag_show_b250 && expected_section_type == SEC_B250) 
        zfile_show_b250_section (section_header_p, uncompressed_buf);
    
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
void zfile_compress_dictionary_data (VBlock *vb, Context *ctx, 
                                     uint32_t num_words, const char *data, uint32_t num_chars)
{
    START_TIMER;

    SectionHeaderDictionary header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = SEC_DICT;
    header.h.data_uncompressed_len = BGEN32 (num_chars);
    header.h.compressed_offset     = BGEN32 (sizeof(SectionHeaderDictionary));
    header.h.codec   = CODEC_BZ2;
    header.h.vblock_i              = BGEN32 (vb->vblock_i);
    header.h.flags                 = (uint8_t)(ctx->flags & 0xff); // 8 bits
    header.num_snips               = BGEN32 (num_words);
    header.dict_id                 = ctx->dict_id;

    if (flag_show_dict) {
        fprintf (stderr, "%s (vb_i=%u, did=%u, num_snips=%u):\t", 
                 ctx->name, vb->vblock_i, ctx->did_i, num_words);
        str_print_null_seperated_data (data, num_chars, true, false);
    }
    
    if (dict_id_printable (ctx->dict_id).num == dict_id_show_one_dict.num)
        str_print_null_seperated_data (data, num_chars, false, false);

    if (flag_list_chroms && ctx->did_i == CHROM)
        str_print_null_seperated_data (data, num_chars, false, vb->data_type == DT_SAM);


    z_file->dict_data.name  = "z_file->dict_data"; // comp_compress requires that it is set in advance
    comp_compress (vb, &z_file->dict_data, true, (SectionHeader*)&header, data, NULL);

    COPY_TIMER (vb->profile.zfile_compress_dictionary_data)    
//printf ("End compress dict vb_i=%u did_i=%u\n", vb->vblock_i, ctx->did_i);
}

void zfile_compress_b250_data (VBlock *vb, Context *ctx, Codec codec)
{
    SectionHeaderCtx header = (SectionHeaderCtx) { 
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_B250,
        .h.data_uncompressed_len = BGEN32 (ctx->b250.len),
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCtx)),
        .h.codec   = codec,
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .h.flags                 = (uint8_t)(ctx->flags & 0xff) | ((ctx->inst & CTX_INST_PAIR_B250) ? CTX_FL_PAIRED : 0), // 8 bit
        .dict_id                 = ctx->dict_id,
        .ltype                   = ctx->ltype
    };
    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, ctx->b250.data, NULL);
}

static LocalGetLineCallback *zfile_get_local_data_callback (DataType dt, Context *ctx)
{
    static struct { DataType dt; const uint64_t *dict_id_num; LocalGetLineCallback *func; } callbacks[] = LOCAL_GET_LINE_CALLBACKS;

    for (unsigned i=0; i < sizeof(callbacks)/sizeof(callbacks[0]); i++)
        if (callbacks[i].dt == dt && *callbacks[i].dict_id_num == ctx->dict_id.num && !(ctx->inst & CTX_INST_NO_CALLBACK)) 
            return callbacks[i].func;

    return NULL;
}

void zfile_compress_local_data (VBlock *vb, Context *ctx)
{   
    vb->has_non_agct = false;
    
    uint8_t flags = ctx->flags | 
                    ((ctx->inst & CTX_INST_PAIR_LOCAL)  ? CTX_FL_PAIRED     : 0) |
                    ((ctx->inst & CTX_INST_LOCAL_PARAM) ? CTX_FL_COPY_PARAM : 0);
                    
    SectionHeaderCtx header = (SectionHeaderCtx) {
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_LOCAL,
        .h.data_uncompressed_len = BGEN32 (ctx->local.len * lt_desc[ctx->ltype].width),
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCtx)),
        .h.codec   = ctx->lcodec,
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .h.flags                 = flags,
        .dict_id                 = ctx->dict_id,
        .ltype                   = ctx->ltype,
        .param                   = (flags & CTX_FL_COPY_PARAM) ? (uint8_t)ctx->local.param : 0
    };

    LocalGetLineCallback *callback = zfile_get_local_data_callback (vb->data_type, ctx);

    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, 
                   callback ? NULL : ctx->local.data, 
                   callback);

    // for CODEC_ACGT, if we discovered any non-A,G,C,T, then we need an additional section with CODEC_NON_ACGT 
    // that will have the correction needed
    if ((ctx->lcodec == CODEC_ACGT) && vb->has_non_agct) {
        header.h.codec = CODEC_NON_ACGT;

        comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, 
                       callback ? NULL : ctx->local.data, 
                       callback);
    }
}

// compress section - two options for input data - 
// 1. contiguous data in section_data 
// 2. line by line data - by providing a callback + total_len
void zfile_compress_section_data_codec (VBlock *vb, SectionType section_type, 
                                        Buffer *section_data,          // option 1 - compress contiguous data
                                        LocalGetLineCallback callback, uint32_t total_len, // option 2 - compress data one line at a time
                                        Codec codec)
{
    SectionHeader header;
    memset (&header, 0, sizeof(header)); // safety

    header.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.section_type          = section_type;
    header.data_uncompressed_len = BGEN32 (section_data ? section_data->len : total_len);
    header.compressed_offset     = BGEN32 (sizeof(header));
    header.codec   = codec;
    header.vblock_i              = BGEN32 (vb->vblock_i);
    
    vb->z_data.name  = "z_data"; // comp_compress requires that these are pre-set
    vb->z_data.param = vb->vblock_i;
    comp_compress (vb, &vb->z_data, false, &header, 
                   section_data ? section_data->data : NULL, 
                   callback);
}

// reads exactly the length required, error otherwise. manages read buffers to optimize I/O performance.
// this doesn't make a big difference for SSD, but makes a huge difference for HD
// return true if data was read as requested, false if file has reached EOF and error otherwise
static void *zfile_read_from_disk (File *file, VBlock *vb, Buffer *buf, unsigned len, bool fail_quietly_if_not_enough_data)
{
    START_TIMER;

    ASSERT0 (len, "Error: in zfile_read_from_disk, len is 0");

    void *start = (void *)&buf->data[buf->len];

    unsigned memcpyied = 0;
    unsigned len_save  = len;

    while (len) {

        if (file->z_next_read == file->z_last_read) { // nothing left in read_buffer - replenish it from disk

            if (file->z_last_read != READ_BUFFER_SIZE && !memcpyied) return NULL; // EOF reached last time, nothing more to read

            if (file->z_last_read != READ_BUFFER_SIZE && memcpyied && fail_quietly_if_not_enough_data) return NULL; // partial read = quiet error as requested
            
            ASSERT (file->z_last_read == READ_BUFFER_SIZE, 
                    "Error: end-of-file while reading %s, read %u bytes, but wanted to read %u bytes", 
                    z_name, memcpyied, len_save);

            file->z_last_read = fread (file->read_buffer, 1, READ_BUFFER_SIZE, (FILE *)file->file);

            file->z_next_read = 0;
            file->disk_so_far += file->z_last_read;

            ASSERT (file->z_last_read || !memcpyied, "Error: data requested could not be read bytes_so_far=%"PRIu64"", file->disk_so_far);
            if (!file->z_last_read) {
                COPY_TIMER (vb->profile.read);
                return NULL; // file is exhausted - nothing read
            }
        }

        unsigned memcpy_len = MIN (len, file->z_last_read - file->z_next_read);

        buf_add (buf, file->read_buffer + file->z_next_read, memcpy_len);
        len                 -= memcpy_len;
        file->z_next_read += memcpy_len;

        memcpyied += memcpy_len;
    }

    COPY_TIMER (vb->profile.read);

    return start;
}

// read section header - called from the I/O thread, but for a specific VB
// returns offset of header within data, EOF if end of file
int32_t zfile_read_section (File *file,
                            VBlock *vb, 
                            uint32_t original_vb_i, // the vblock_i used for compressing. this is part of the encryption key. dictionaries are compressed by the compute thread/vb, but uncompressed by the I/O thread (vb=0)
                            Buffer *data, const char *buf_name, // buffer to append 
                            unsigned header_size, SectionType expected_sec_type,
                            const SectionListEntry *sl)   // NULL for no seeking
{
    ASSERT (!sl || expected_sec_type == sl->section_type, "Error in zfile_read_section: expected_sec_type=%s but encountered sl->section_type=%s. vb_i=%u",
            st_name (expected_sec_type), st_name(sl->section_type), vb->vblock_i);

    if (sl && file == z_file && piz_is_skip_section (vb, expected_sec_type, sl->dict_id)) return 0; // skip if this section is not needed according to flags

    unsigned unencrypted_header_size = header_size;

    bool is_encrypted =  (z_file->data_type != DT_REF) && 
                         (expected_sec_type != SEC_GENOZIP_HEADER) &&
                         crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    unsigned header_offset = data->len;
    buf_alloc (vb, data, header_offset + header_size, 2, buf_name, 1);
    
    // move the cursor to the section. file_seek is smart not to cause any overhead if no moving is needed
    if (sl) file_seek (file, sl->offset, SEEK_SET, false);

    SectionHeader *header = zfile_read_from_disk (file, vb, data, header_size, false); // note: header in file can be shorter than header_size if its an earlier version

    // case: we're done! no more bound files
    if (!header && expected_sec_type == SEC_TXT_HEADER) return EOF; 

    ASSERT (header, "Error in zfile_read_section: Failed to read data from file %s while expecting section type %s: %s", 
            z_name, st_name(expected_sec_type), strerror (errno));
    
    if (flag_show_headers) zfile_show_header (header, NULL);

    bool is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;

    // SEC_REFERENCE is never encrypted when originating from a reference file, it is encrypted (if the file is encrypted) if it originates from REF_INTERNAL 
    if (is_encrypted && header->section_type == SEC_REFERENCE && !header->data_encrypted_len) {
        is_encrypted = false;
        header_size = unencrypted_header_size;
    }

    // decrypt header (note: except for SEC_GENOZIP_HEADER - this header is never encrypted)
    if (is_encrypted) {
        ASSERT (BGEN32 (header->magic) != GENOZIP_MAGIC, 
                "Error: password provided, but file %s is not encrypted (sec_type=%s)", z_name, st_name (header->section_type));

        crypt_do (vb, (uint8_t*)header, header_size, original_vb_i, expected_sec_type, true); 
    
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
    ASSERT (expected_sec_type == header->section_type || 
            (expected_sec_type == SEC_GENOZIP_HEADER && header->flags == SEC_GENOZIP_HEADER), // in v2-5, the section_type field was located where flags is now
            "Error: Unexpected section type when reading %s: expecting %s, found %s", 
            z_name, st_name(expected_sec_type), st_name(header->section_type));

    ASSERT (compressed_offset == header_size || expected_sec_type == SEC_GENOZIP_HEADER || // we allow SEC_GENOZIP_HEADER of other sizes, for older versions
            "Error: invalid header when reading %s - expecting compressed_offset to be %u but found %u. section_type=%s", 
            z_name, header_size, compressed_offset, st_name(header->section_type));

    // allocate more memory for the rest of the header + data (note: after this realloc, header pointer is no longer valid)
    buf_alloc (vb, data, header_offset + compressed_offset + data_len, 2, "zfile_read_section", 2);
    header = (SectionHeader *)&data->data[header_offset]; // update after realloc

    // read section data - but only if header size is as expected
    if (data_len && compressed_offset == header_size) {
        ASSERT (zfile_read_from_disk (file, vb, data, data_len, false), 
                "Error: failed to read section data, section_type=%s: %s", st_name(header->section_type), strerror (errno));
    }

    return header_offset;
}

// Read one section header - NOT thread-safe, should only be called by the I/O thread
void *zfile_read_section_header (uint64_t offset, uint32_t size)
{
    // get the uncompressed size from one of the headers - they are all the same size, and the reference file is never encrypted
    file_seek (z_file, offset, SEEK_SET, false);

    static Buffer one_header_buf = EMPTY_BUFFER;
    buf_alloc (evb, &one_header_buf, size, 4, "one_header_buf", 0);
    one_header_buf.len = 0;

    return zfile_read_from_disk (z_file, evb, &one_header_buf, size, false); 
}

// PIZ
void zfile_read_all_dictionaries (uint32_t last_vb_i /* 0 means all VBs */, ReadChromeType read_chrom)
{
    SectionListEntry *sl_ent = NULL; // NULL -> first call to this sections_get_next_section_of_type() will reset cursor 

    mtf_initialize_primary_field_ctxs (z_file->contexts, z_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_contexts);

    while (sections_get_next_section_of_type (&sl_ent, &z_file->sl_dir_cursor, SEC_DICT, SEC_NONE)) {

        if (last_vb_i && sl_ent->vblock_i > last_vb_i) break;

        // cases where we can skip reading these dictionaries because we don't be using them
        bool is_chrom = (sl_ent->dict_id.num == dict_id_fields[CHROM]);
        if (read_chrom == DICTREAD_CHROM_ONLY  && !is_chrom) continue;
        if (read_chrom == DICTREAD_EXCEPT_CHROM && is_chrom) continue;

        if (piz_is_skip_sectionz (sl_ent->section_type, sl_ent->dict_id)) continue;
        
        zfile_read_section (z_file, evb, sl_ent->vblock_i, &evb->z_data, "z_data", sizeof(SectionHeaderDictionary), sl_ent->section_type, sl_ent);    

        // update dictionaries in z_file->contexts with dictionary data 
        mtf_integrate_dictionary_fragment (evb, evb->z_data.data);

        buf_free (&evb->z_data);
    }

    // output the dictionaries if we're asked to
    if (flag_show_dict || dict_id_show_one_dict.num || flag_list_chroms) {
        for (uint32_t did_i=0; did_i < z_file->num_contexts; did_i++) {
            Context *ctx = &z_file->contexts[did_i];

#define MAX_PRINTABLE_DICT_LEN 100000

            if (dict_id_printable (ctx->dict_id).num == dict_id_show_one_dict.num) 
                str_print_null_seperated_data (ctx->dict.data, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), false, false);
            
            if (flag_list_chroms && ctx->did_i == CHROM)
                str_print_null_seperated_data (ctx->dict.data, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), false, z_file->data_type == DT_SAM);
            
            if (flag_show_dict) {
                fprintf (stderr, "%s (did_i=%u, num_snips=%u):\t", ctx->name, did_i, (uint32_t)ctx->word_list.len);
                str_print_null_seperated_data (ctx->dict.data, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), true, false);
            }
        }
        fprintf (stderr, "\n");

        if (exe_type == EXE_GENOCAT) exit_ok; // if this is genocat - we're done
    }
}

// returns the file's data_type
int16_t zfile_read_genozip_header (Md5Hash *digest) // out
{
    DataType data_type = DT_NONE;

    // read the footer from the end of the file
    file_seek (z_file, -sizeof(SectionFooterGenozipHeader), SEEK_END, false);

    SectionFooterGenozipHeader footer;
    int ret = fread (&footer, sizeof (footer), 1, (FILE *)z_file->file);
    ASSERTW (ret == 1, "Skipping empty file %s", z_name);
    if (!ret) goto final;
    
    // case: there is no genozip header. this can happen if the file was truncated (eg because compression did not complete)
    // note: this can also happen if the file is genozip v1, but I don't think there are any real v1 files in the wild
    // so I will keep the error message simple and not mention it
    ASSERT (BGEN32 (footer.magic) == GENOZIP_MAGIC, "Error: failed to read file %s - the file appears to be incomplete.", z_name);

    // read genozip header
    uint64_t footer_offset = BGEN64 (footer.genozip_header_offset);
    
    SectionListEntry dummy_sl = { .section_type = SEC_GENOZIP_HEADER,
                                  .offset       = footer_offset };

    // header might be smaller for older versions - we limit our reading of it to the entire section size so we don't
    // fail due to end-of-file. This is just so we can observe the section number, and give a proper error message
    // for unsupported version<=5 files.
    unsigned sizeof_genozip_header = MIN (sizeof (SectionHeaderGenozipHeader),
                                          (unsigned)(z_file->disk_size - footer_offset - sizeof(SectionFooterGenozipHeader)));
    
    ret = zfile_read_section (z_file, evb, 0, &evb->z_data, "genozip_header", sizeof_genozip_header, SEC_GENOZIP_HEADER, &dummy_sl);

    ASSERT0 (ret != EOF, "Error: unexpected EOF when reading genozip header");
    
    SectionHeaderGenozipHeader *header = (SectionHeaderGenozipHeader *)evb->z_data.data;

    ASSERT (header->genozip_version <= GENOZIP_FILE_FORMAT_VERSION, 
            "Error: %s cannot be openned because it was compressed with a newer version of genozip (version %u.x.x) while the version you're running is older (version %s).\n"
            "You might want to consider upgrading genozip to the newest version.\n",
            z_name, header->genozip_version, GENOZIP_CODE_VERSION);

    // in version 6, we canceled backward compatability with v1-v5
    ASSERT (header->genozip_version >= 6, "Error: %s was compressed with an older version of genozip - version %u.\nIt may be uncompressed with genozip versions %u to 5",
            z_name, header->genozip_version, header->genozip_version);

    // in version 7, we canceled backward compatability with v6
    ASSERT (header->genozip_version >= 7, "Error: %s was compressed with an older version of genozip - version %u.\nIt may be uncompressed with genozip versions 6",
            z_name, header->genozip_version);

    data_type = (DataType)(BGEN16 (header->data_type)); 
    ASSERT ((unsigned)data_type < NUM_DATATYPES, "Error in zfile_read_genozip_header: unrecognized data_type=%d", data_type);

    if (z_file->data_type == DT_NONE) {
        z_file->data_type = data_type;
        z_file->type      = file_get_z_ft_by_dt (z_file->data_type);  
    }
    else
        ASSERT (z_file->data_type == data_type, "Error: %s - file extension indicates this is a %s file, but according to its contents it is a %s", 
                z_name, dt_name (z_file->data_type), dt_name (data_type));

    if (txt_file) txt_file->data_type = data_type; // txt_file is still NULL in case of --1d
    
    ASSERT (header->encryption_type != ENCRYPTION_TYPE_NONE || !crypt_have_password() || z_file->data_type == DT_REF, 
            "Error: password provided, but file %s is not encrypted", z_name);

    ASSERT (BGEN32 (header->h.compressed_offset) == sizeof (SectionHeaderGenozipHeader),
            "Error: invalid genozip header - expecting compressed_offset to be %u but found %u", 
            (unsigned)sizeof (SectionHeaderGenozipHeader), BGEN32 (header->h.compressed_offset));

    // get & test password, if file is encrypted
    if (header->encryption_type != ENCRYPTION_TYPE_NONE) {

        if (!crypt_have_password()) crypt_prompt_for_password();

        crypt_do (evb, header->password_test, sizeof(header->password_test), 0, SEC_NONE, true); // decrypt password test

        ASSERT (!memcmp (header->password_test, password_test_string, sizeof(header->password_test)),
                "Error: password is wrong for file %s", z_name);
    }

    global_vcf_num_samples    = BGEN32 (header->num_samples); // possibly 0, if genozip header was not rewritten. in this case, piz will get it from the first VCF header, but genols will show 0
    z_file->num_components    = BGEN32 (header->num_components);
    z_file->genozip_version   = header->genozip_version;
    z_file->flags             = header->h.flags;
    if (digest) *digest       = header->md5_hash_bound; 
    
    // global bools to help testing
    is_v6_or_above = (z_file->genozip_version >= 6);
         
    zfile_uncompress_section (evb, header, &z_file->section_list_buf, "z_file->section_list_buf", SEC_GENOZIP_HEADER);
    z_file->section_list_buf.len /= sizeof (SectionListEntry); // fix len
    BGEN_sections_list();

    if (flag_show_gheader) {
        sections_show_gheader (header);
        if (exe_type == EXE_GENOCAT) exit_ok; // in genocat, exit after showing the requested data
    }

    // case: we are reading a file expected to be the reference file itself
    if (flag_reading_reference) {
        ASSERT (header->data_type == DT_REF, "Error: %s is not a reference file. To create a reference file, use 'genozip --make-reference <fasta-file.fa>'",
                ref_filename);

        ref_set_ref_file_info (header->md5_hash_bound, header->ref_filename); // in the reference file itself, header->ref_filename is the original fasta used to create this reference
    }

    // case: we are reading a file that is not expected to be a reference file
    else {
        // case: we are attempting to decompress a reference file - this is not supported
        if (data_type == DT_REF && !(flag_genocat_info_only && exe_type == EXE_GENOCAT)) { // we will stop a bit later in this case
            WARN ("%s is a reference file - it cannot be decompressed. Skipping it.", z_name);
            data_type = DT_NONE;
            goto final;
        }

        if (flag_show_reference && !md5_is_zero (header->ref_file_md5)) {
            fprintf (stderr, "%s was compressed using the reference file:\nName: %s\nMD5: %s\n",
                     z_name, header->ref_filename, md5_display (header->ref_file_md5));
            if (exe_type == EXE_GENOCAT) exit_ok; // in genocat --show-reference, we only show the reference, not the data
        }

        ASSERT (md5_is_zero (header->ref_file_md5) || flag_reference == REF_EXTERNAL || flag_genocat_info_only, 
                "Error reading %s: please use --reference to specify the reference filename.\nNote: the reference used to compress this file was %s",
                z_name, header->ref_filename);

        if (!(flag_reference == REF_NONE || flag_reference == REF_INTERNAL || md5_is_equal (header->ref_file_md5, ref_md5))) {
    
            // fail if user loaded an external reference, but file does not require one
            ASSINP (!md5_is_zero (header->ref_file_md5), "Error: %s does not need --reference to decompress", txt_name);

            // just warn, don't fail - there are use cases where the user might do this on purpose
            ASSERTW (false, 
                    "WARNING: The reference file has a different MD5 than the reference file used to compress %s\n"
                    "If these two files contain a different sequence in the genomic regions contained in the compressed file, then \n"
                    "THE UNCOMPRESSED FILE WILL BE DIFFERENT THAN ORIGINAL FILE\n"
                    "Reference you are using now: %s MD5=%s\n"
                    "Reference used to compress the file: %s MD5=%s\n", 
                    z_name, ref_filename, md5_display (ref_md5), 
                    header->ref_filename, md5_display (header->ref_file_md5));
        }
    }
     
final:
    buf_free (&evb->z_data);
    return data_type;
}

void zfile_compress_genozip_header (Md5Hash single_component_md5)
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
    header.h.codec   = CODEC_BZ2;
    header.h.flags                 = ((flag_reference == REF_INTERNAL ? GENOZIP_FL_REF_INTERNAL : 0) |
                                      (flag_ref_use_aligner           ? GENOZIP_FL_ALIGNER      : 0) );
    header.genozip_version         = GENOZIP_FILE_FORMAT_VERSION;
    header.data_type               = BGEN16 ((uint16_t)z_file->data_type);
    header.encryption_type         = is_encrypted ? ENCRYPTION_TYPE_AES256 : ENCRYPTION_TYPE_NONE;
    header.uncompressed_data_size  = BGEN64 (z_file->txt_data_so_far_bind);
    header.num_samples             = BGEN32 (global_vcf_num_samples);
    header.num_items_bound         = BGEN64 (z_file->num_lines);
    header.num_sections            = BGEN32 (num_sections); 
    header.num_components          = BGEN32 (z_file->num_txt_components_so_far);
    
    // when decompressing will require an external reference, we set header.ref_filename to the name of the genozip reference file
    if (flag_reference == REF_EXTERNAL) {   
        strncpy (header.ref_filename, ref_filename, REF_FILENAME_LEN-1);
        header.ref_file_md5 = ref_md5;
    }

    // in --make-ref, we set header.ref_filename to the original fasta file, to be used later in ref_get_cram_ref
    // (unless the fasta is piped from stdin, or its name is too long)
    else if (flag_make_reference && strcmp (txt_name, FILENAME_STDIN) && strlen (txt_name) <= REF_FILENAME_LEN-1) {
#ifndef WIN32
        realpath (txt_name, header.ref_filename);
#else
        _fullpath (header.ref_filename, txt_name, REF_FILENAME_LEN);
#endif
    }

    uint32_t license_num_bgen = BGEN32 (license_get());
    header.license_hash = md5_do (&license_num_bgen, sizeof (int32_t));

    if (flag_md5) {
        if (flag_bind) {
            // get the hash from a copy of the md5 context, because we still need it for displaying the compression ratio later
            Md5Context copy_md5_ctx_bound = z_file->md5_ctx_bound;
            header.md5_hash_bound = md5_finalize (&copy_md5_ctx_bound);
        } 
        else 
            header.md5_hash_bound = single_component_md5; // if not in bound mode - just copy the md5 of the single file
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

// reads the the genozip header section's header from a GENOZIP file - used by main_list and main_is_fastq, returns true if successful
bool zfile_get_genozip_header (File *file,
                               DataType *dt,
                               uint64_t *uncomp_data_size,
                               uint32_t *num_samples,
                               uint64_t *num_items_bound,
                               Md5Hash  *md5_hash_bound,
                               char *created, unsigned created_len, // caller allocates space 
                               Md5Hash  *license_hash,
                               char *ref_file_name, unsigned ref_file_name_len,  // caller allocates space 
                               Md5Hash *ref_file_md5)
{
    // read the footer from the end of the file
    if (!file_seek (file, -sizeof(SectionFooterGenozipHeader), SEEK_END, true))
        return false;

    SectionFooterGenozipHeader footer;
    int ret = fread (&footer, sizeof (footer), 1, (FILE *)file->file);
    ASSERTW (ret == 1, "Skipping empty file %s", z_name);    
    if (!ret) return false; // empty file / cannot read
    
    // case: this is not a valid genozip v2+ file
    if (BGEN32 (footer.magic) != GENOZIP_MAGIC) return false;

    // read genozip header
    uint64_t genozip_header_offset = BGEN64 (footer.genozip_header_offset);
    if (!file_seek (file, genozip_header_offset, SEEK_SET, true))
        return false;

    SectionHeaderGenozipHeader header;
    int bytes = fread ((char*)&header, 1, sizeof(SectionHeaderGenozipHeader), (FILE *)file->file);
    if (bytes < sizeof(SectionHeaderGenozipHeader)) return false;

    ASSERTW (BGEN32 (header.h.magic) == GENOZIP_MAGIC, "Error reading %s: corrupt data", z_name);
    if (BGEN32 (header.h.magic) != GENOZIP_MAGIC) return false;

    if (dt)               *dt               = (DataType)BGEN16 (header.data_type);
    if (uncomp_data_size) *uncomp_data_size = BGEN64 (header.uncompressed_data_size);
    if (num_samples)      *num_samples      = BGEN32 (header.num_samples);
    if (num_items_bound)  *num_items_bound  = BGEN64 (header.num_items_bound);
    if (md5_hash_bound)   *md5_hash_bound   = header.md5_hash_bound;
    if (license_hash)     *license_hash     = header.license_hash;
    if (ref_file_md5)     *ref_file_md5     = header.ref_file_md5;

    memset (created, 0, created_len);
    memcpy (created, header.created, MIN (FILE_METADATA_LEN, created_len));

    memset (ref_file_name, 0, ref_file_name_len);
    memcpy (ref_file_name, header.ref_filename, MIN (REF_FILENAME_LEN, ref_file_name_len));

    return true;
}

// ZIP
void zfile_write_txt_header (Buffer *txt_header_text, Md5Hash header_md5, bool is_first_txt)
{
    SectionHeaderTxtHeader header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = SEC_TXT_HEADER;
    header.h.data_uncompressed_len = BGEN32 (txt_header_text->len);
    header.h.compressed_offset     = BGEN32 (sizeof (SectionHeaderTxtHeader));
    header.h.codec   = CODEC_BZ2;
    header.num_samples             = BGEN32 (global_vcf_num_samples);
    header.num_lines               = NUM_LINES_UNKNOWN; 
    header.compression_type        = (uint8_t)txt_file->codec; 
    header.md5_header              = header_md5;
    
    file_basename (txt_file->name, false, FILENAME_STDIN, header.txt_filename, TXT_FILENAME_LEN);
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
    z_file->txt_data_so_far_bind += txt_header_text->len;

    buf_free (&txt_header_buf); 

    // copy it to z_file - we might need to update it at the very end in zfile_update_txt_header_section_header()
    if (is_first_txt)
        memcpy (&z_file->txt_header_first, &header, sizeof (header));

    memcpy (&z_file->txt_header_single, &header, sizeof (header));
}

// updating the VCF bytes of a GENOZIP file. If we're compressing a simple VCF file, we will know
// the bytes upfront, but if we're binding or compressing a VCF.GZ, we will need to update it
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
    curr_header->md5_hash_single  = flag_md5 ? md5_finalize (&z_file->md5_ctx_single) : MD5HASH_NONE;

    *md5 = curr_header->md5_hash_single;

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
void zfile_compress_vb_header (VBlock *vb)
{
    uint32_t sizeof_header = sizeof (SectionHeaderVbHeader);

    SectionHeaderVbHeader vb_header = {
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_VB_HEADER,
        .h.compressed_offset     = BGEN32 (sizeof_header),
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .h.codec   = CODEC_NONE,
        .num_lines               = BGEN32 ((uint32_t)vb->lines.len),
        .vb_data_size            = BGEN32 (vb->vb_data_size),
        .longest_line_len        = BGEN32 (vb->longest_line_len),
        .md5_hash_so_far         = vb->md5_hash_so_far,
        .num_samples             = BGEN32 (global_vcf_num_samples), // 0 if not VCF
        .num_haplotypes_per_line = BGEN32 (vb->num_haplotypes_per_line), // 0 if not VCF
    };

    // copy section header into z_data - to be eventually written to disk by the I/O thread. this section doesn't have data.
    vb->z_data.name  = "z_data"; // this is the first allocation of z_data - comp_compress requires that it is pre-named
    vb->z_data.param = vb->vblock_i;
    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&vb_header, NULL, NULL);
}

// ZIP only: called by the I/O thread in the sequential order of VBs: updating of the already compressed
// variant data section (compressed by the compute thread in zfile_compress_vb_header) just before writing it to disk
// note: this updates the z_data in memory (not on disk)
void zfile_update_compressed_vb_header (VBlock *vb, uint32_t txt_first_line_i)
{
    SectionHeaderVbHeader *vb_header = (SectionHeaderVbHeader *)vb->z_data.data;
    vb_header->z_data_bytes = BGEN32 ((uint32_t)vb->z_data.len);
    vb_header->first_line   = BGEN32 (txt_first_line_i);

    if (flag_show_vblocks) 
        fprintf (stderr, "vb_i=%u first_line=%u num_lines=%u txt_file=%u genozip_size=%u longest_line_len=%u\n",
                 vb->vblock_i, txt_first_line_i, BGEN32 (vb_header->num_lines), 
                 BGEN32 (vb_header->vb_data_size), BGEN32 (vb_header->z_data_bytes), 
                 BGEN32 (vb_header->longest_line_len));

    // now we can finally encrypt the header - if needed
    if (crypt_have_password())  
        crypt_do (vb, (uint8_t*)vb_header, BGEN32 (vb_header->h.compressed_offset),
                  BGEN32 (vb_header->h.vblock_i), vb_header->h.section_type, true);
}

