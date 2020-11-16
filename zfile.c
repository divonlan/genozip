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
#include "codec.h"
#include "piz.h"
#include "zip.h"
#include "license.h"
#include "mutex.h"
#include "strings.h"
#include "dict_id.h"
#include "reference.h"

static const char *password_test_string = "WhenIThinkBackOnAllTheCrapIlearntInHighschool";

void zfile_show_header (const SectionHeader *header, VBlock *vb /* optional if output to buffer */, uint64_t offset, char rw)
{
    if (flag.reading_reference) return; // don't show headers of reference file
    
    DictId dict_id = {0};
    char flags[10] = "", param[10] = "";
    const char *ltype="";
    bool has_ltype = (header->section_type == SEC_LOCAL);
    bool has_param = (header->section_type == SEC_LOCAL || header->section_type == SEC_B250);
    bool is_dict_offset = (header->section_type == SEC_DICT && rw == 'W'); // at the point calling this function in zip, SEC_DICT offsets are not finalized yet and are relative to the beginning of the dictionary area in the genozip file
    if (header->section_type == SEC_DICT) 
        dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    
    else if (header->section_type == SEC_LOCAL || header->section_type == SEC_B250) {
        SectionHeaderCtx *header_ctx = (SectionHeaderCtx *)header;
        dict_id = header_ctx->dict_id;
        str_int (header_ctx->param, param);
        if (has_ltype) ltype = lt_desc[header_ctx->ltype].name;
    }

    if (header->flags) 
        str_int (header->flags, flags);

    char str[1000];
    #define PRINT if (vb) buf_add_string (vb, &vb->show_headers_buf, str); else printf ("%s", str)

    sprintf (str, "%c %s%-*"PRIu64" %-19s %*.*s %4s%-3s %5s%-3s %4s%-3s %-4.4s %-4.4s vb=%-3u z_off=%-6u txt_len=%-7u z_len=%-7u enc_len=%-7u mgc=%8.8x\n",
             rw, 
             is_dict_offset ? "~" : "", 9-is_dict_offset, offset, 
             st_name(header->section_type), 
             -DICT_ID_LEN, DICT_ID_LEN, dict_id.num ? dict_id_printable (dict_id).id : dict_id.id,
             header->flags ? "flg=" : "", flags, 
             has_ltype ? "type=" : "", ltype, 
             has_param ? "prm=" : "", param, 
             codec_name (header->codec), codec_name (header->sub_codec),
             BGEN32 (header->vblock_i), 
             BGEN32 (header->compressed_offset), 
             BGEN32 (header->data_uncompressed_len), 
             BGEN32 (header->data_compressed_len), 
             BGEN32 (header->data_encrypted_len), 
             BGEN32 (header->magic));
    PRINT;

#define SEC_TAB "            ++  "
    if (header->section_type == SEC_GENOZIP_HEADER) {
        SectionHeaderGenozipHeader *h = (SectionHeaderGenozipHeader *)header;
        sprintf (str, SEC_TAB "ver=%u enc=%s dt=%s usize=%"PRIu64" lines=%"PRIu64" secs=%u txts=%u md5bound=%s md5ref=%s \n" 
                      SEC_TAB "created=\"%.*s\" ref=\"%.*s\"\n",
                 h->genozip_version, encryption_name (h->encryption_type), dt_name (BGEN16 (h->data_type)), 
                 BGEN64 (h->uncompressed_data_size), BGEN64 (h->num_items_bound), BGEN32 (h->num_sections), BGEN32 (h->num_components),
                 md5_display (h->md5_hash_bound), md5_display (h->ref_file_md5), FILE_METADATA_LEN, h->created, REF_FILENAME_LEN, h->ref_filename);
        PRINT;
    }

    else if (header->section_type == SEC_TXT_HEADER) {
        SectionHeaderTxtHeader *h = (SectionHeaderTxtHeader *)header;
        sprintf (str, SEC_TAB "txt_size=%"PRIu64" lines=%"PRIu64" max_lines_per_vb=%u md5_single=%s md5_header=%s\n" 
                      SEC_TAB "txt_codec=%s (args=0x%02X.%02X.%02X) txt_filename=\"%.*s\"\n",
                 BGEN64 (h->txt_data_size), BGEN64 (h->num_lines), BGEN32 (h->max_lines_per_vb), 
                 md5_display (h->md5_hash_single), md5_display (h->md5_header), 
                 codec_name (h->codec), h->codec_args[0], h->codec_args[1], h->codec_args[2], TXT_FILENAME_LEN, h->txt_filename);
        PRINT;
    }

    else if (header->section_type == SEC_VB_HEADER) {
        SectionHeaderVbHeader *h = (SectionHeaderVbHeader *)header;
        sprintf (str, SEC_TAB "first_line=%u lines=%u longest_line=%u vb_data_size=%u z_data_bytes=%u md5_hash_so_far=%s\n",
                 BGEN32 (h->first_line), BGEN32 (h->num_lines), BGEN32 (h->longest_line_len), BGEN32 (h->vb_data_size), 
                 BGEN32 (h->z_data_bytes), md5_display (h->md5_hash_so_far));
        PRINT;
    }

    else if (header->section_type == SEC_REFERENCE || header->section_type == SEC_REF_IS_SET) {
        SectionHeaderReference *h = (SectionHeaderReference *)header;
        sprintf (str, SEC_TAB "pos=%"PRIu64" gpos=%"PRIu64" num_bases=%u chrom_word_index=%u\n",
                 BGEN64 (h->pos), BGEN64 (h->gpos), BGEN32 (h->num_bases), BGEN32 (h->chrom_word_index)); 
        PRINT;
    }
    
    else if (header->section_type == SEC_REF_HASH) {
        SectionHeaderRefHash *h = (SectionHeaderRefHash *)header;
        sprintf (str, SEC_TAB "num_layers=%u layer_i=%u layer_bits=%u start_in_layer=%u\n",
                 h->num_layers, h->layer_i, h->layer_bits, BGEN32 (h->start_in_layer)); 
        PRINT;
    }
    

}

static void zfile_show_b250_section (void *section_header_p, const Buffer *b250_data)
{
    MUTEX (show_b250_mutex); // protect so compute thread's outputs don't get mix

    SectionHeaderCtx *header = (SectionHeaderCtx *)section_header_p;

    if (!flag.show_b250 && dict_id_printable (header->dict_id).num != flag.dict_id_show_one_b250.num) return;

    mutex_initialize (show_b250_mutex); // possible unlikely race condition on initializing - good enough for debugging purposes
    mutex_lock (show_b250_mutex);

    fprintf (stderr, "vb_i=%u %*.*s: ", BGEN32 (header->h.vblock_i), -DICT_ID_LEN-1, DICT_ID_LEN, dict_id_printable (header->dict_id).id);

    const uint8_t *data  = FIRSTENT (const uint8_t, *b250_data);
    const uint8_t *after = AFTERENT (const uint8_t, *b250_data);

    while (data < after) {
        WordIndex word_index = base250_decode (&data, true);
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

// Write uncompressed, unencrypted section to <section-type>.<vb>.<dict_id>.[header|body]. 
// Note: header includes encryption padding if it was encrypted
static void zfile_dump_section (Buffer *uncompressed_data, SectionHeader *section_header, unsigned section_header_len, DictId dict_id)
{
    char filename[100];
    uint32_t vb_i = BGEN32 (section_header->vblock_i);

    // header
    sprintf (filename, "%s.%u.%s.header", st_name (section_header->section_type), vb_i, dict_id_print (dict_id));
    file_put_data (filename, section_header, section_header_len);

    // body
    sprintf (filename, "%s.%u.%s.body", st_name (section_header->section_type), vb_i, dict_id_print (dict_id));
    file_put_buffer (filename, uncompressed_data, 1);
}

// uncompressed a block and adds a \0 at its end. Returns the length of the uncompressed block, without the \0.
// when we get here, the header is already unencrypted zfile_one_section
void zfile_uncompress_section (VBlock *vb,
                               void *section_header_p,
                               Buffer *uncompressed_data, 
                               const char *uncompressed_data_buf_name, // a name if Buffer, NULL ok if buffer need not be realloced
                               uint32_t expected_vb_i,
                               SectionType expected_section_type) 
{
    START_TIMER;

    DictId dict_id = DICT_ID_NONE;
    if (expected_section_type == SEC_DICT)
        dict_id = ((SectionHeaderDictionary *)section_header_p)->dict_id;
    else if (expected_section_type == SEC_B250 || expected_section_type == SEC_LOCAL)
        dict_id = ((SectionHeaderCtx *)section_header_p)->dict_id;

    if (piz_is_skip_section (vb, expected_section_type, dict_id)) return; // we skip some sections based on flags

    SectionHeader *section_header  = (SectionHeader *)section_header_p;
    uint32_t compressed_offset     = BGEN32 (section_header->compressed_offset);
    uint32_t data_encrypted_len    = BGEN32 (section_header->data_encrypted_len);
    uint32_t data_compressed_len   = BGEN32 (section_header->data_compressed_len);
    uint32_t data_uncompressed_len = BGEN32 (section_header->data_uncompressed_len);
    uint32_t vblock_i              = BGEN32 (section_header->vblock_i);

    // sanity checks
    ASSERT (section_header->section_type == expected_section_type, "Error in zfile_uncompress_section: expecting section type %s but seeing %s", st_name(expected_section_type), st_name(section_header->section_type));
    
    ASSERT (vblock_i == expected_vb_i || !expected_vb_i, // dictionaries are uncompressed by the I/O thread with pseduo_vb (vb_i=0) 
            "Error in zfile_uncompress_section: bad vblock_i: vblock_i in file=%u but expecting it to be %u (section_type=%s)", 
            vblock_i, expected_vb_i, st_name (expected_section_type));

    // decrypt data (in-place) if needed
    if (data_encrypted_len) 
        crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, vblock_i, section_header->section_type, false);

    if (data_uncompressed_len > 0) { // FORMAT, for example, can be missing in a sample-less file

        if (uncompressed_data_buf_name) { 
            buf_alloc (vb, uncompressed_data, data_uncompressed_len + sizeof (uint64_t), 1.1, uncompressed_data_buf_name); // add a 64b word for safety in case this buffer will be converted to a bitarray later
            uncompressed_data->len = data_uncompressed_len;
        }

        comp_uncompress (vb, section_header->codec, section_header->sub_codec, (char*)section_header + compressed_offset, data_compressed_len, 
                         uncompressed_data, data_uncompressed_len);
    }
 
    if (flag.show_b250 && expected_section_type == SEC_B250) 
        zfile_show_b250_section (section_header_p, uncompressed_data);
    
    if (flag.dump_section && !strcmp (st_name (expected_section_type), flag.dump_section))
        zfile_dump_section (uncompressed_data, section_header, compressed_offset, dict_id);

    if (vb) COPY_TIMER (zfile_uncompress_section);
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
    Codec codec = CODEC_BSC;

    SectionHeaderDictionary header = (SectionHeaderDictionary){ 
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_DICT,
        .h.data_uncompressed_len = BGEN32 (num_chars),
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderDictionary)),
        .h.codec                 = codec,
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .num_snips               = BGEN32 (num_words),
        .dict_id                 = ctx->dict_id
    };

    if (flag.show_dict) {
        fprintf (stderr, "%s (vb_i=%u, did=%u, num_snips=%u):\t", 
                 ctx->name, vb->vblock_i, ctx->did_i, num_words);
        str_print_null_seperated_data (data, num_chars, true, false);
    }
    
    if (dict_id_printable (ctx->dict_id).num == flag.dict_id_show_one_dict.num)
        str_print_null_seperated_data (data, num_chars, false, false);

    if (flag.list_chroms && ctx->did_i == CHROM)
        str_print_null_seperated_data (data, num_chars, false, vb->data_type == DT_SAM);

    if (flag.show_time) codec_show_time (vb, "DICT", ctx->name, codec);

    z_file->dict_data.name  = "z_file->dict_data"; // comp_compress requires that it is set in advance
    comp_compress (vb, &z_file->dict_data, true, (SectionHeader*)&header, data, NULL);

    COPY_TIMER (zfile_compress_dictionary_data)    
//printf ("End compress dict vb_i=%u did_i=%u\n", vb->vblock_i, ctx->did_i);
}

uint32_t zfile_compress_b250_data (VBlock *vb, Context *ctx)
{
    SectionHeaderCtx header = (SectionHeaderCtx) { 
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_B250,
        .h.data_uncompressed_len = BGEN32 (ctx->b250.len),
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCtx)),
        .h.codec                 = ctx->bcodec == CODEC_UNKNOWN ? CODEC_BZ2 : ctx->bcodec,
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .h.flags                 = ctx->flags | ((ctx->inst & CTX_INST_PAIR_B250) ? CTX_FL_PAIRED : 0), 
        .dict_id                 = ctx->dict_id,
        .ltype                   = ctx->ltype
    };
    
    return comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, ctx->b250.data, NULL);
}

static LocalGetLineCB *zfile_get_local_data_callback (DataType dt, Context *ctx)
{
    static struct { DataType dt; const uint64_t *dict_id_num; LocalGetLineCB *func; } callbacks[] = LOCAL_GET_LINE_CALLBACKS;

    for (unsigned i=0; i < sizeof(callbacks)/sizeof(callbacks[0]); i++)
        if (callbacks[i].dt == dt && *callbacks[i].dict_id_num == ctx->dict_id.num && !(ctx->inst & CTX_INST_NO_CALLBACK)) 
            return callbacks[i].func;

    return NULL;
}

// returns compressed size
uint32_t zfile_compress_local_data (VBlock *vb, Context *ctx, uint32_t sample_size /* 0 means entire local buffer */)
{   
    vb->has_non_agct = false;
    
    SectionFlags flags = ctx->flags | 
                         ((ctx->inst & CTX_INST_PAIR_LOCAL)  ? CTX_FL_PAIRED     : 0) |
                         ((ctx->inst & CTX_INST_LOCAL_PARAM) ? CTX_FL_COPY_PARAM : 0);
                    
    uint32_t uncompressed_len = ctx->local.len * lt_desc[ctx->ltype].width;
    
    // case: we're just testing a small sample
    if (sample_size && uncompressed_len > sample_size) 
        uncompressed_len = sample_size;

    SectionHeaderCtx header = (SectionHeaderCtx) {
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_LOCAL,
        .h.data_uncompressed_len = BGEN32 (uncompressed_len),
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCtx)),
        .h.codec                 = ctx->lcodec == CODEC_UNKNOWN ? CODEC_BZ2 : ctx->lcodec, // if codec has not been decided yet, fall back on BZ2
        .h.sub_codec             = ctx->lsubcodec_piz ? ctx->lsubcodec_piz : codec_args[ctx->lcodec].sub_codec,
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .h.flags                 = flags,
        .dict_id                 = ctx->dict_id,
        .ltype                   = ctx->ltype,
        .param                   = (flags & CTX_FL_COPY_PARAM) ? (uint8_t)ctx->local.param : 0
    };

    LocalGetLineCB *callback = zfile_get_local_data_callback (vb->data_type, ctx);

    return comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, 
                          callback ? NULL : ctx->local.data, 
                          callback);
}

// compress section - two options for input data - 
// 1. contiguous data in section_data 
// 2. line by line data - by providing a callback + total_len
void zfile_compress_section_data_codec (VBlock *vb, SectionType section_type, 
                                        Buffer *section_data,          // option 1 - compress contiguous data
                                        LocalGetLineCB callback, uint32_t total_len, // option 2 - compress data one line at a time
                                        Codec codec)
{
    SectionHeader header = (SectionHeader){
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .section_type          = section_type,
        .data_uncompressed_len = BGEN32 (section_data ? section_data->len : total_len),
        .compressed_offset     = BGEN32 (sizeof(header)),
        .codec                 = codec,
        .vblock_i              = BGEN32 (vb->vblock_i)
    };

    if (flag.show_time) codec_show_time (vb, st_name (section_type), NULL, codec);

    vb->z_data.name  = "z_data"; // comp_compress requires that these are pre-set
    comp_compress (vb, &vb->z_data, false, &header, 
                   section_data ? section_data->data : NULL, 
                   callback);
}

// reads exactly the length required, error otherwise. manages read buffers to optimize I/O performance.
// this doesn't make a big difference for SSD, but makes a huge difference for HD
// return a pointer to the data read
static void *zfile_read_from_disk (File *file, VBlock *vb, Buffer *buf, uint32_t len, SectionType st)
{
    START_TIMER;

    ASSERT (len, "Error in zfile_read_from_disk reading %s: len is 0", st_name (st));
    ASSERT (buf_has_space (buf, len), "Error in zfile_read_from_disk reading %s: buf is out of space: len=%u but remaining space in buffer=%u (tip: run with --show-headers to see where it fails)",
            st_name (st), len, (uint32_t)(buf->size - buf->len));

    char *start = AFTERENT (char, *buf);
    uint32_t bytes = fread (start, 1, len, (FILE *)file->file);
    ASSERT (bytes == len, "Error in zfile_read_from_disk reading %s: read only %u bytes out of len=%u", st_name (st), bytes, len);

    buf->len += bytes;

    COPY_TIMER (read);

    return start;
}


// read section header - called from the I/O thread. returns offset of header within data
int32_t zfile_read_section_do (File *file,
                               VBlock *vb, 
                               uint32_t original_vb_i, // the vblock_i used for compressing. this is part of the encryption key. dictionaries are compressed by the compute thread/vb, but uncompressed by the I/O thread (vb=0)
                               Buffer *data, const char *buf_name, // buffer to append 
                               SectionType expected_sec_type,
                               const SectionListEntry *sl,
                               uint32_t header_size)   // NULL for no seeking
{
    ASSERT (!sl || expected_sec_type == sl->section_type, "Error in zfile_read_section: expected_sec_type=%s but encountered sl->section_type=%s. vb_i=%u",
            st_name (expected_sec_type), st_name(sl->section_type), vb->vblock_i);

    if (sl && file == z_file && piz_is_skip_section (vb, expected_sec_type, sl->dict_id)) return 0; // skip if this section is not needed according to flags
    uint32_t unencrypted_header_size = header_size;

    // note: for an encrypted file, while reading the reference, we don't yet know until getting the header whether it
    // will be an SEC_REF_IS_SET (encrypted) or SEC_REFERENCE (not encrypted if originating from external, encryptd if de-novo)
    bool is_encrypted =  (z_file->data_type != DT_REF) && 
                         (expected_sec_type != SEC_GENOZIP_HEADER) &&
                         crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    uint32_t header_offset = data->len;
    buf_alloc (vb, data, header_offset + header_size, 2, buf_name);
    data->param = 1;
    
    // move the cursor to the section. file_seek is smart not to cause any overhead if no moving is needed
    if (sl) file_seek (file, sl->offset, SEEK_SET, false);

    SectionHeader *header = zfile_read_from_disk (file, vb, data, header_size, expected_sec_type); // note: header in file can be shorter than header_size if its an earlier version
    uint32_t bytes_read = header_size;

    ASSERT (header, "Error in zfile_read_section: Failed to read data from file %s while expecting section type %s: %s", 
            z_name, st_name(expected_sec_type), strerror (errno));
    
    bool is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;

    // SEC_REFERENCE is never encrypted when originating from a reference file, it is encrypted (if the file is encrypted) if it originates from REF_INTERNAL 
    if (is_encrypted && header->section_type == SEC_REFERENCE && !header->data_encrypted_len) {
        is_encrypted = false;
        header_size  = unencrypted_header_size;
    }

    // decrypt header (note: except for SEC_GENOZIP_HEADER - this header is never encrypted)
    if (is_encrypted) {
        ASSERT (BGEN32 (header->magic) != GENOZIP_MAGIC, 
                "Error: password provided, but file %s is not encrypted (sec_type=%s)", z_name, st_name (header->section_type));

        crypt_do (vb, (uint8_t*)header, header_size, original_vb_i, expected_sec_type, true); 
    
        is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC; // update after decryption
    }

    if (flag.show_headers) {
        zfile_show_header (header, NULL, sl ? sl->offset : 0, 'R');
        if (exe_type == EXE_GENOCAT && (expected_sec_type == SEC_B250 || expected_sec_type == SEC_LOCAL || expected_sec_type == SEC_DICT))
             return header_offset; // in genocat --show-header - we only show headers, nothing else
    }

    ASSERT (is_magical, "Error in zfile_read_section: corrupt data (magic is wrong) when attempting to read section=%s dict_id=%s of vblock_i=%u component=%u in file %s", 
            st_name (expected_sec_type), sl ? err_dict_id (sl->dict_id) : "(no sl)", vb->vblock_i, z_file->num_txt_components_so_far, z_name);

    uint32_t compressed_offset   = BGEN32 (header->compressed_offset);
    ASSERT (compressed_offset, "Error: header.compressed_offset is 0 when reading section_type=%s", st_name(expected_sec_type));

    uint32_t data_compressed_len = BGEN32 (header->data_compressed_len);
    uint32_t data_encrypted_len  = BGEN32 (header->data_encrypted_len);

    uint32_t data_len = MAX (data_compressed_len, data_encrypted_len);

    // in case where we already read part of the body (eg if is_encrypted was initially set and then unset) (remaining_data_len might be negative)
    int32_t remaining_data_len = (int32_t)data_len - (int32_t)(bytes_read - header_size); 
    
    // check that we received the section type we expect, 
    char s[30];
    ASSERT (expected_sec_type == header->section_type || 
            (expected_sec_type == SEC_GENOZIP_HEADER && header->flags == SEC_GENOZIP_HEADER), // in v2-5, the section_type field was located where flags is now
            "Error: Unexpected section type when reading %s: expecting %s, found %s sl(expecting)=(offset=%s, dict_id=%s)",
            z_name, st_name(expected_sec_type), st_name(header->section_type), 
            sl ? str_uint_commas (sl->offset, s) : "N/A", sl ? err_dict_id (sl->dict_id) : "N/A");

    ASSERT (compressed_offset == header_size || expected_sec_type == SEC_GENOZIP_HEADER || // we allow SEC_GENOZIP_HEADER of other sizes, for older versions
            "Error: invalid header when reading %s - expecting compressed_offset to be %u but found %u. section_type=%s", 
            z_name, header_size, compressed_offset, st_name(header->section_type));

    // allocate more memory for the rest of the header + data (note: after this realloc, header pointer is no longer valid)
    buf_alloc (vb, data, header_offset + compressed_offset + data_len, 2, "zfile_read_section");
    data->param = 2;
    header = (SectionHeader *)&data->data[header_offset]; // update after realloc

    // read section data - but only if header size is as expected
    if (remaining_data_len > 0 && compressed_offset == header_size)
        zfile_read_from_disk (file, vb, data, remaining_data_len, expected_sec_type);

    return header_offset;
}

// Read one section header - returns the header in vb->compressed - caller needs to free vb->compressed
void *zfile_read_section_header (VBlockP vb, uint64_t offset, 
                                 uint32_t original_vb_i, // the vblock_i used for compressing. this is part of the encryption key. dictionaries are compressed by the compute thread/vb, but uncompressed by the I/O thread (vb=0)
                                 SectionType expected_sec_type)
{
    uint32_t header_size = st_header_size (expected_sec_type);

    // get the uncompressed size from one of the headers - they are all the same size, and the reference file is never encrypted
    file_seek (z_file, offset, SEEK_SET, false);

    bool is_encrypted =  (z_file->data_type != DT_REF) && 
                         (expected_sec_type != SEC_GENOZIP_HEADER) &&
                         crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
        
    ASSERT (!vb->compressed.len, "Error in zfile_read_section_header vb_i=%u expected_sec_type=%s: expecting vb->compressed to be free, but its not",
            vb->vblock_i, st_name (expected_sec_type));
    
    buf_alloc (evb, &vb->compressed, header_size, 4, "compressed"); 

    SectionHeader *header = zfile_read_from_disk (z_file, evb, &vb->compressed, header_size, expected_sec_type); 

    ASSERT (header, "Error in zfile_read_section_header: Failed to read header of section type %s from file %s: %s", 
            st_name(expected_sec_type), z_name, strerror (errno));

    bool is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;

    // decrypt header 
    if (is_encrypted) {
        ASSERT (BGEN32 (header->magic) != GENOZIP_MAGIC, 
                "Error in zfile_read_section_header: password provided, but file %s is not encrypted (sec_type=%s)", z_name, st_name (header->section_type));

        crypt_do (vb, (uint8_t*)header, header_size, original_vb_i, expected_sec_type, true); 
    
        is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC; // update after decryption
    }

    ASSERT (is_magical, "Error in zfile_read_section_header: corrupt data (magic is wrong) when attempting to read header of section %s in file %s", 
            st_name (expected_sec_type), z_name);

    return header;
}

// PIZ
void zfile_read_all_dictionaries (uint32_t last_vb_i /* 0 means all VBs */, ReadChromeType read_chrom)
{
    ctx_initialize_primary_field_ctxs (z_file->contexts, z_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_contexts);

    const SectionListEntry *sl_ent = NULL; 
    while (sections_get_next_section_of_type (&sl_ent, SEC_DICT, true, false)) {

        if (last_vb_i && sl_ent->vblock_i > last_vb_i) break;

        // cases where we can skip reading these dictionaries because we don't be using them
        bool is_chrom = (sl_ent->dict_id.num == dict_id_fields[CHROM]);
        if (read_chrom == DICTREAD_CHROM_ONLY  && !is_chrom) continue;
        if (read_chrom == DICTREAD_EXCEPT_CHROM && is_chrom) continue;

        if (piz_is_skip_sectionz (sl_ent->section_type, sl_ent->dict_id)) continue;
        
        zfile_read_section (z_file, evb, sl_ent->vblock_i, &evb->z_data, "z_data", sl_ent->section_type, sl_ent);    

        // update dictionaries in z_file->contexts with dictionary data 
        if (!(flag.show_headers && exe_type == EXE_GENOCAT))
            ctx_integrate_dictionary_fragment (evb, evb->z_data.data);

        buf_free (&evb->z_data);
    }

    // output the dictionaries if we're asked to
    if (flag.show_dict || flag.dict_id_show_one_dict.num || flag.list_chroms) {
        for (uint32_t did_i=0; did_i < z_file->num_contexts; did_i++) {
            Context *ctx = &z_file->contexts[did_i];

#define MAX_PRINTABLE_DICT_LEN 100000

            if (dict_id_printable (ctx->dict_id).num == flag.dict_id_show_one_dict.num) 
                str_print_null_seperated_data (ctx->dict.data, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), false, false);
            
            if (flag.list_chroms && ctx->did_i == CHROM)
                str_print_null_seperated_data (ctx->dict.data, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), false, z_file->data_type == DT_SAM);
            
            if (flag.show_dict) {
                fprintf (stderr, "%s (did_i=%u, num_snips=%u):\t", ctx->name, did_i, (uint32_t)ctx->word_list.len);
                str_print_null_seperated_data (ctx->dict.data, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), true, false);
            }
        }
        fprintf (stderr, "\n");

        if (exe_type == EXE_GENOCAT) exit_ok; // if this is genocat - we're done
    }
}

// returns false if file should be skipped
bool zfile_read_genozip_header (Md5Hash *digest, uint64_t *txt_data_size, uint64_t *num_items_bound, char *created) // optional outs
{
    bool success=false;

    // read the footer from the end of the file
    if (!file_seek (z_file, -sizeof(SectionFooterGenozipHeader), SEEK_END, 2))
        goto final; // likely an empty file 

    SectionFooterGenozipHeader footer;
    int ret = fread (&footer, sizeof (footer), 1, (FILE *)z_file->file);
    ASSERTW (ret == 1, "Skipping empty file %s", z_name);
    if (!ret) goto final;
    
    // case: there is no genozip header. this can happen if the file was truncated (eg because compression did not complete)
    // note: this can also happen if the file is genozip v1, but I don't think there are any real v1 files in the wild
    // so I will keep the error message simple and not mention it
    ASSERTGOTO (BGEN32 (footer.magic) == GENOZIP_MAGIC, "Error: failed to read file %s - the file appears to be incomplete.", z_name);

    // read genozip header
    uint64_t footer_offset = BGEN64 (footer.genozip_header_offset);
    
    SectionListEntry dummy_sl = { .section_type = SEC_GENOZIP_HEADER,
                                  .offset       = footer_offset };

    // header might be smaller for older versions - we limit our reading of it to the entire section size so we don't
    // fail due to end-of-file. This is just so we can observe the section number, and give a proper error message
    // for unsupported version<=5 files.
    uint32_t sizeof_genozip_header = MIN (sizeof (SectionHeaderGenozipHeader),
                                          (uint32_t)(z_file->disk_size - footer_offset - sizeof(SectionFooterGenozipHeader)));
    
    zfile_read_section_do (z_file, evb, 0, &evb->z_data, "genozip_header", SEC_GENOZIP_HEADER, &dummy_sl, sizeof_genozip_header);

    SectionHeaderGenozipHeader *header = (SectionHeaderGenozipHeader *)evb->z_data.data;

    ASSERT (header->genozip_version <= GENOZIP_FILE_FORMAT_VERSION, 
            "Error: %s cannot be openned because it was compressed with a newer version of genozip (version %u.x.x) while the version you're running is older (version %s).\n"
            "You might want to consider upgrading genozip to the newest version.\n",
            z_name, header->genozip_version, GENOZIP_CODE_VERSION);

    // in version 6, we canceled backward compatability with v1-v5
    ASSERT (header->genozip_version >= 6, "Error: %s was compressed with an older version of genozip - version %u.\nIt may be uncompressed with genozip versions %u to 5",
            z_name, header->genozip_version, header->genozip_version);

    // in version 7, we canceled backward compatability with v6
    ASSERT (header->genozip_version >= 7, "Error: %s was compressed with an version 6 of genozip - version %u.\nIt may be uncompressed with genozip versions 6",
            z_name, header->genozip_version);

    // in version 8, we canceled backward compatability with v7
    ASSERT (header->genozip_version >= 8, "Error: %s was compressed with an version 7 of genozip - version %u.\nIt may be uncompressed with genozip versions 7",
            z_name, header->genozip_version);

    DataType data_type = (DataType)(BGEN16 (header->data_type)); 
    ASSERT ((unsigned)data_type < NUM_DATATYPES, "Error in zfile_read_genozip_header: unrecognized data_type=%d", data_type);

    if (z_file->data_type == DT_NONE) {
        z_file->data_type = data_type;
        z_file->type      = file_get_z_ft_by_dt (z_file->data_type);  
    }
    else
        ASSERT (z_file->data_type == data_type, "Error: %s - file extension indicates this is a %s file, but according to its contents it is a %s", 
                z_name, dt_name (z_file->data_type), dt_name (data_type));

    if (txt_file) txt_file->data_type = data_type; // txt_file is still NULL in case of --1d
        
    ASSERT (header->encryption_type != ENC_NONE || !crypt_have_password() || z_file->data_type == DT_REF, 
            "Error: password provided, but file %s is not encrypted", z_name);

    ASSERT (BGEN32 (header->h.compressed_offset) == sizeof (SectionHeaderGenozipHeader),
            "Error: invalid genozip header - expecting compressed_offset to be %u but found %u", 
            (unsigned)sizeof (SectionHeaderGenozipHeader), BGEN32 (header->h.compressed_offset));

    // get & test password, if file is encrypted
    if (header->encryption_type != ENC_NONE) {

        if (!crypt_have_password()) crypt_prompt_for_password();

        crypt_do (evb, header->password_test, sizeof(header->password_test), 0, SEC_NONE, true); // decrypt password test

        ASSERT (!memcmp (header->password_test, password_test_string, sizeof(header->password_test)),
                "Error: password is wrong for file %s", z_name);
    }

    z_file->num_components    = BGEN32 (header->num_components);
    if (z_file->num_components < 2) flag.unbind = 0; // override user's --unbind if file has only 1 component

    z_file->genozip_version   = header->genozip_version;
    z_file->flags             = header->h.flags;
    if (digest) *digest       = header->md5_hash_bound; 
    if (txt_data_size) *txt_data_size = BGEN64 (header->uncompressed_data_size);
    if (num_items_bound) *num_items_bound = BGEN64 (header->num_items_bound); 
    if (created) memcpy (created, header->created, FILE_METADATA_LEN);

    zfile_uncompress_section (evb, header, &z_file->section_list_buf, "z_file->section_list_buf", 0, SEC_GENOZIP_HEADER);
    z_file->section_list_buf.len /= sizeof (SectionListEntry); // fix len
    BGEN_sections_list();

    if (flag.show_gheader) {
        sections_show_gheader (header);
        if (exe_type == EXE_GENOCAT) exit_ok; // in genocat, exit after showing the requested data
    }

    // case: we are reading a file expected to be the reference file itself
    if (flag.reading_reference) {
        ASSERT (data_type == DT_REF, "Error: %s is not a reference file. To create a reference file, use 'genozip --make-reference <fasta-file.fa>'",
                ref_filename);

        ref_set_ref_file_info (header->md5_hash_bound, header->ref_filename); // in the reference file itself, header->ref_filename is the original fasta used to create this reference
    }

    // case: we are reading a file that is not expected to be a reference file
    else {
        // case: we are attempting to decompress a reference file - this is not supported
        if (data_type == DT_REF && !(flag.genocat_info_only && exe_type == EXE_GENOCAT) && exe_type != EXE_GENOLS) { // we will stop a bit later in this case
            WARN ("%s is a reference file - it cannot be decompressed. Skipping it.", z_name);
            goto final;
        }

        if (flag.show_reference && !md5_is_zero (header->ref_file_md5)) {
            fprintf (stderr, "%s was compressed using the reference file:\nName: %s\nMD5: %s\n",
                     z_name, header->ref_filename, md5_display (header->ref_file_md5));
            if (exe_type == EXE_GENOCAT) exit_ok; // in genocat --show-reference, we only show the reference, not the data
        }

        if (!(flag.reference == REF_NONE || flag.reference == REF_INTERNAL || md5_is_equal (header->ref_file_md5, ref_md5))) {
    
            // just warn, don't fail - there are use cases where the user might do this on purpose
            ASSERTW (md5_is_zero (header->ref_file_md5), // its ok if file doesn't need a reference
                    "WARNING: The reference file has a different MD5 than the reference file used to compress %s\n"
                    "If these two files contain a different sequence in the genomic regions contained in the compressed file, then \n"
                    "THE UNCOMPRESSED FILE WILL BE DIFFERENT THAN ORIGINAL FILE\n"
                    "Reference you are using now: %s MD5=%s\n"
                    "Reference used to compress the file: %s MD5=%s\n", 
                    z_name, ref_filename, md5_display (ref_md5), 
                    header->ref_filename, md5_display (header->ref_file_md5));
        }

        // case: this file requires an external reference, but command line doesn't include --reference - attempt to use the
        // reference specified in the header. 
        // Note: this code will be executed when zfile_read_genozip_header is called from main_genounzip.
        if (!md5_is_zero (header->ref_file_md5) && !ref_filename && exe_type != EXE_GENOLS) {
            ASSERTW (flag.genocat_info_only, "Note: using the reference file %s. You can override this with --reference", header->ref_filename);
            ref_set_reference (header->ref_filename);
            flag.reference = REF_EXTERNAL;
        }
    }
    success = true;
     
error: // for ASSERTGOTO
final:
    buf_free (&evb->z_data);
    return success;
}

void zfile_compress_genozip_header (Md5Hash single_component_md5)
{
    SectionHeaderGenozipHeader header = {};

    // start with just the fields needed by sections_add_to_list
    header.h.section_type = SEC_GENOZIP_HEADER;

    // "manually" add the genozip section to the section list - normally it is added in comp_compress()
    // but in this case the genozip section containing the list will already be ready...
    sections_add_to_list (evb, &header.h);

    bool is_encrypted = crypt_have_password();

    uint32_t num_sections = z_file->section_list_buf.len;

    // BAM txt files result in SAM genozip files 
    DataType z_data_type =  (z_file->data_type == DT_BAM ? DT_SAM : z_file->data_type);

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.compressed_offset     = BGEN32 (sizeof (SectionHeaderGenozipHeader));
    header.h.data_uncompressed_len = BGEN32 (z_file->section_list_buf.len * sizeof (SectionListEntry));
    header.h.codec                 = CODEC_BZ2;
    header.h.flags                 = z_file->flags |
                                     (flag.reference == REF_INTERNAL ? GENOZIP_FL_REF_INTERNAL : 0) |
                                     (flag.ref_use_aligner           ? GENOZIP_FL_ALIGNER      : 0);
    header.genozip_version         = GENOZIP_FILE_FORMAT_VERSION;
    header.data_type               = BGEN16 ((uint16_t)z_data_type);
    header.encryption_type         = is_encrypted ? ENC_AES256 : ENC_NONE;
    header.uncompressed_data_size  = BGEN64 (z_file->txt_data_so_far_bind);
    header.num_items_bound         = BGEN64 (z_file->num_lines);
    header.num_sections            = BGEN32 (num_sections); 
    header.num_components          = BGEN32 (z_file->num_txt_components_so_far);
    
    // when decompressing will require an external reference, we set header.ref_filename to the name of the genozip reference file
    if (flag.reference == REF_EXTERNAL) {   
        strncpy (header.ref_filename, ref_filename, REF_FILENAME_LEN-1);
        header.ref_file_md5 = ref_md5;
    }

    // in --make-ref, we set header.ref_filename to the original fasta file, to be used later in ref_get_cram_ref
    // (unless the fasta is piped from stdin, or its name is too long)
    else if (flag.make_reference && strcmp (txt_name, FILENAME_STDIN) && strlen (txt_name) <= REF_FILENAME_LEN-1) {
#ifndef WIN32
        char *ref_filename = realpath (txt_name, NULL); // allocates memory
        ASSERT (ref_filename, "Error in zfile_compress_genozip_header: realpath() failed: %s", strerror (errno));

        strncpy (header.ref_filename, ref_filename, REF_FILENAME_LEN-1);
        FREE (ref_filename);
#else
        ASSERT0 (_fullpath (header.ref_filename, txt_name, REF_FILENAME_LEN), "Error in zfile_compress_genozip_header: _fullpath() failed");
#endif
    }

    uint32_t license_num_bgen = BGEN32 (license_get());
    header.license_hash = md5_do (&license_num_bgen, sizeof (int32_t));

    if (flag.md5) {
        if (flag.bind) {
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

    if (flag.show_gheader) sections_show_gheader (&header); 

    // add a footer to this section - this footer appears AFTER the genozip header data, 
    // facilitating reading the genozip header in reverse from the end of the file
    SectionFooterGenozipHeader footer;
    footer.magic                 = BGEN32 (GENOZIP_MAGIC);
    footer.genozip_header_offset = BGEN64 (genozip_header_offset);

    buf_alloc (evb, z_data, z_data->len + sizeof(SectionFooterGenozipHeader), 1.5, "z_data");
    memcpy (&z_data->data[z_data->len], &footer, sizeof(SectionFooterGenozipHeader));
    z_data->len += sizeof(SectionFooterGenozipHeader);
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
    header.h.codec                 = CODEC_BZ2;
    header.num_lines               = NUM_LINES_UNKNOWN; 
    header.codec                   = txt_file->codec; 
    header.md5_header              = header_md5;
    
    // In BGZF, we store the 3 least significant bytes of the file size, so check if the reconstructed BGZF file is likely the same
    if (txt_file->codec == CODEC_BGZF) {
        header.codec_args[0] = (txt_file->disk_size      ) & 0xff; // LSB of size
        header.codec_args[1] = (txt_file->disk_size >> 8 ) & 0xff;
        header.codec_args[2] = (txt_file->disk_size >> 16) & 0xff;
    }
        
    file_basename (txt_file->name, false, FILENAME_STDIN, header.txt_filename, TXT_FILENAME_LEN);
    file_remove_codec_ext (header.txt_filename, txt_file->type); // eg "xx.fastq.gz -> xx.fastq"
    
    static Buffer txt_header_buf = EMPTY_BUFFER;

    buf_alloc (evb, &txt_header_buf, sizeof (SectionHeaderTxtHeader) + txt_header_text->len / 3, // generous guess of compressed size
               1, "txt_header_buf"); 

    comp_compress (evb, &txt_header_buf, true, (SectionHeader*)&header, 
                   txt_header_text->len ? txt_header_text->data : NULL, // actual header may be missing (eg in SAM it is permitted to not have a header)
                   NULL);

    file_write (z_file, txt_header_buf.data, txt_header_buf.len);

    z_file->disk_so_far            += txt_header_buf.len;   // length of GENOZIP data writen to disk
    z_file->txt_data_so_far_single += txt_header_text->len; // length of the original VCF header
    z_file->txt_data_so_far_bind   += txt_header_text->len;

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

    uint32_t len = crypt_padded_len (sizeof (SectionHeaderTxtHeader));

    // update the header of the single (current) vcf. 
    SectionHeaderTxtHeader *curr_header = &z_file->txt_header_single;
    curr_header->txt_data_size    = BGEN64 (txt_file->txt_data_size_single);
    curr_header->num_lines        = BGEN64 (txt_file->num_lines);
    curr_header->max_lines_per_vb = BGEN32 (max_lines_per_vb);
    curr_header->md5_hash_single  = flag.md5 ? md5_finalize (&z_file->md5_ctx_single) : MD5HASH_NONE;

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

// ZIP compute thread - called from zip_compress_one_vb()
void zfile_compress_vb_header (VBlock *vb)
{
    uint32_t sizeof_header = sizeof (SectionHeaderVbHeader);

    SectionHeaderVbHeader vb_header = {
        .h.magic             = BGEN32 (GENOZIP_MAGIC),
        .h.section_type      = SEC_VB_HEADER,
        .h.compressed_offset = BGEN32 (sizeof_header),
        .h.vblock_i          = BGEN32 (vb->vblock_i),
        .h.codec             = CODEC_NONE,
        .num_lines           = BGEN32 ((uint32_t)vb->lines.len),
        .vb_data_size        = BGEN32 (vb->vb_data_size),
        .longest_line_len    = BGEN32 (vb->longest_line_len),
        .md5_hash_so_far     = vb->md5_hash_so_far
    };

    // copy section header into z_data - to be eventually written to disk by the I/O thread. this section doesn't have data.
    vb->z_data.name  = "z_data"; // this is the first allocation of z_data - comp_compress requires that it is pre-named
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

    if (flag.show_vblocks) 
        fprintf (stderr, "vb_i=%u first_line=%u num_lines=%u txt_file=%u genozip_size=%u longest_line_len=%u\n",
                 vb->vblock_i, txt_first_line_i, BGEN32 (vb_header->num_lines), 
                 BGEN32 (vb_header->vb_data_size), BGEN32 (vb_header->z_data_bytes), 
                 BGEN32 (vb_header->longest_line_len));

    // now we can finally encrypt the header - if needed
    if (crypt_have_password())  
        crypt_do (vb, (uint8_t*)vb_header, BGEN32 (vb_header->h.compressed_offset),
                  BGEN32 (vb_header->h.vblock_i), vb_header->h.section_type, true);
}

DataType zfile_get_file_dt (const char *filename)
{
    FileType ft = file_get_type (filename, false);
    DataType dt = file_get_dt_by_z_ft (ft);

    // case: we don't know yet what file type this is - we need to read the genozip header to determine
    if (dt == DT_NONE && filename) {
        File *file = file_open (filename, READ, Z_FILE, DT_NONE);

        // read the footer from the end of the file
        if (!file_seek (file, -sizeof(SectionFooterGenozipHeader), SEEK_END, true))
            goto done;

        SectionFooterGenozipHeader footer;
        int ret = fread (&footer, sizeof (footer), 1, (FILE *)file->file);
        ASSERTW (ret == 1, "Skipping empty file %s", z_name);    
        if (!ret) goto done; // empty file / cannot read
        
        // case: this is not a valid genozip v2+ file
        if (BGEN32 (footer.magic) != GENOZIP_MAGIC) goto done;

        // read genozip header
        uint64_t genozip_header_offset = BGEN64 (footer.genozip_header_offset);
        if (!file_seek (file, genozip_header_offset, SEEK_SET, true))
            goto done;

        SectionHeaderGenozipHeader header;
        int bytes = fread ((char*)&header, 1, sizeof(SectionHeaderGenozipHeader), (FILE *)file->file);
        if (bytes < sizeof(SectionHeaderGenozipHeader)) goto done;

        ASSERTW (BGEN32 (header.h.magic) == GENOZIP_MAGIC, "Error reading %s: corrupt data", z_name);
        if (BGEN32 (header.h.magic) != GENOZIP_MAGIC) goto done;

        dt = (DataType)BGEN16 (header.data_type);

        file_close (&file, false);
    }

done:
    return dt;
}

