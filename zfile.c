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
#include "dispatcher.h"
#include "bgzf.h"
#include "digest.h"
#include "md5.h"

static const char *password_test_string = "WhenIThinkBackOnAllTheCrapIlearntInHighschool";

void zfile_show_header (const SectionHeader *header, VBlock *vb /* optional if output to buffer */, uint64_t offset, char rw)
{
    if (flag.reading_reference || flag.reading_chain) return; // don't show headers of reference file
    
    if (  flag.show_headers   != -1 &&                 // we don't need to show all sections
          flag.show_headers-1 != header->section_type) // we don't need to show this section
        return;

    bool is_dict_offset = (header->section_type == SEC_DICT && rw == 'W'); // at the point calling this function in zip, SEC_DICT offsets are not finalized yet and are relative to the beginning of the dictionary area in the genozip file

    char str[1000];
    #define PRINT { if (vb) buf_add_string (vb, &vb->show_headers_buf, str); else printf ("%s", str); } 

    sprintf (str, "%c %s%-*"PRIu64" %-19s %-4.4s %-4.4s vb=%-3u z_off=%-6u txt_len=%-7u z_len=%-7u enc_len=%-7u mgc=%8.8x\n",
             rw, 
             is_dict_offset ? "~" : "", 9-is_dict_offset, offset, 
             st_name(header->section_type), 
             codec_name (header->codec), codec_name (header->sub_codec),
             BGEN32 (header->vblock_i), 
             BGEN32 (header->compressed_offset), 
             BGEN32 (header->data_uncompressed_len), 
             BGEN32 (header->data_compressed_len), 
             BGEN32 (header->data_encrypted_len), 
             BGEN32 (header->magic));
    PRINT;

#define SEC_TAB "            ++  "
    
    switch (header->section_type) {
    case SEC_GENOZIP_HEADER: {
        SectionHeaderGenozipHeader *h = (SectionHeaderGenozipHeader *)header;
        struct FlagsGenozipHeader f = h->h.flags.genozip_header;
        DataType dt = BGEN16 (h->data_type);
        sprintf (str, SEC_TAB "ver=%u enc=%s dt=%s usize=%"PRIu64" lines=%"PRIu64" secs=%u txts=%u digest_bound=%s\n" 
                      SEC_TAB "%s=%u aligner=%u txt_is_bin=%u bgzf=%u adler=%u dual_coords=%u ref=\"%.*s\" md5ref=%s\n"
                      SEC_TAB "created=\"%.*s\"\n",
                 h->genozip_version, encryption_name (h->encryption_type), dt_name (dt), 
                 BGEN64 (h->uncompressed_data_size), BGEN64 (h->num_items_bound), BGEN32 (h->num_sections), BGEN32 (h->num_components),
                 digest_display (h->digest_bound).s, 
                 (dt==DT_FASTQ ? "dts_paired" : "dts_ref_internal"), f.dt_specific, 
                 f.aligner, f.txt_is_bin, f.bgzf, f.adler, f.dual_coords,
                 REF_FILENAME_LEN, h->ref_filename, digest_display_ex (h->ref_file_md5, DD_MD5).s,
                 FILE_METADATA_LEN, h->created);
        break;
    }

    case SEC_TXT_HEADER: {
        SectionHeaderTxtHeader *h = (SectionHeaderTxtHeader *)header;
        sprintf (str, SEC_TAB "txt_data_size=%"PRIu64" lines=%"PRIu64" max_lines_per_vb=%u md5_single=%s digest_header=%s\n" 
                      SEC_TAB "txt_codec=%s (args=0x%02X.%02X.%02X) liftover_rejects=%u txt_filename=\"%.*s\"\n",
                 BGEN64 (h->txt_data_size), BGEN64 (h->num_lines), BGEN32 (h->max_lines_per_vb), 
                 digest_display (h->digest_single).s, digest_display (h->digest_header).s, 
                 codec_name (h->codec), h->codec_info[0], h->codec_info[1], h->codec_info[2], 
                 h->h.flags.txt_header.liftover_rejects, TXT_FILENAME_LEN, h->txt_filename);
        break;
    }

    case SEC_VB_HEADER: {
        SectionHeaderVbHeader *h = (SectionHeaderVbHeader *)header;
        sprintf (str, SEC_TAB "first_line=%u lines=%u longest_line=%u vb_data_size=%u z_data_bytes=%u digest_so_far=%s\n",
                 BGEN32 (h->first_line), BGEN32 (h->num_lines), BGEN32 (h->longest_line_len), BGEN32 (h->vb_data_size), 
                 BGEN32 (h->z_data_bytes), digest_display (h->digest_so_far).s);
        break;
    }

    case SEC_REFERENCE:
    case SEC_REF_IS_SET: {
        SectionHeaderReference *h = (SectionHeaderReference *)header;
        sprintf (str, SEC_TAB "pos=%"PRIu64" gpos=%"PRIu64" num_bases=%u chrom_word_index=%u\n",
                 BGEN64 (h->pos), BGEN64 (h->gpos), BGEN32 (h->num_bases), BGEN32 (h->chrom_word_index)); 
        break;
    }
    
    case SEC_REF_HASH: {
        SectionHeaderRefHash *h = (SectionHeaderRefHash *)header;
        sprintf (str, SEC_TAB "num_layers=%u layer_i=%u layer_bits=%u start_in_layer=%u\n",
                 h->num_layers, h->layer_i, h->layer_bits, BGEN32 (h->start_in_layer)); 
        break;
    }
    
    case SEC_RECON_PLAN: {
        SectionHeaderReconPlan *h = (SectionHeaderReconPlan *)header;
        sprintf (str, SEC_TAB "num_txt_data_bufs=%u laft=%u\n", h->num_txt_data_bufs, h->h.flags.recon_plan.laft); 
        break;
    }
    
    case SEC_BGZF: {
        SectionHeaderRefHash *h = (SectionHeaderRefHash *)header;
        sprintf (str, SEC_TAB "level=%u has_eof=%u\n", h->h.flags.bgzf.level, h->h.flags.bgzf.has_eof_block); 
        break;
    }
    
    case SEC_B250: {
        SectionHeaderCtx *h = (SectionHeaderCtx *)header;
        static const char *store[4] = { [STORE_NONE]="NONE", [STORE_INT]="INT", [STORE_FLOAT]="FLOAT", [STORE_INDEX]="INDEX"};

        sprintf (str, SEC_TAB "%s param=%u store=%s paired=%u copy_param=%u all_the_same=%u ctx_specific=%u\n",
                 dis_dict_id_name (h->dict_id).s, h->param, store[h->h.flags.ctx.store], 
                 h->h.flags.ctx.paired, h->h.flags.ctx.copy_param, h->h.flags.ctx.all_the_same, h->h.flags.ctx.ctx_specific); 
        break;
    }

    case SEC_LOCAL: {
        SectionHeaderCtx *h = (SectionHeaderCtx *)header;
        sprintf (str, SEC_TAB "%s ltype=%s param=%u paired=%u copy_param=%u ctx_specific=%u\n",
                 dis_dict_id_name (h->dict_id).s, lt_name (h->ltype), h->param, 
                 h->h.flags.ctx.paired, h->h.flags.ctx.copy_param, h->h.flags.ctx.ctx_specific); 
        break;
    }

    case SEC_DICT: {
        SectionHeaderDictionary *h = (SectionHeaderDictionary *)header;
        sprintf (str, SEC_TAB "%s num_snips=%u\n", dis_dict_id (h->dict_id).s, BGEN32 (h->num_snips)); 
        break;
    }

    default: 
        str[0] = 0; 
    }

    if (str[0]) PRINT;
}

static void zfile_show_b250_section (void *section_header_p, const Buffer *b250_data)
{
    static Mutex show_b250_mutex = {}; // protect so compute thread's outputs don't get mix

    SectionHeaderCtx *header = (SectionHeaderCtx *)section_header_p;

    if (!flag.show_b250 && dict_id_typeless (header->dict_id).num != flag.dict_id_show_one_b250.num) return;

    mutex_initialize (show_b250_mutex); // possible unlikely race condition on initializing - good enough for debugging purposes
    mutex_lock (show_b250_mutex);

    fprintf (info_stream, "vb_i=%u %*.*s: ", BGEN32 (header->h.vblock_i), -DICT_ID_LEN-1, DICT_ID_LEN, dict_id_typeless (header->dict_id).id);

    const uint8_t *data  = FIRSTENT (const uint8_t, *b250_data);
    const uint8_t *after = AFTERENT (const uint8_t, *b250_data);

    while (data < after) {
        WordIndex word_index = base250_decode (&data, true, "zfile_show_b250_section");
        switch (word_index) {
            case WORD_INDEX_ONE_UP     : fprintf (info_stream, "ONE_UP "); break;
            case WORD_INDEX_EMPTY_SF   : fprintf (info_stream, "EMPTY "); break;
            case WORD_INDEX_MISSING_SF : fprintf (info_stream, "MISSING "); break;
            default: fprintf (info_stream, "%u ", word_index);
        }
    }
    fprintf (info_stream, "\n");

    fflush (info_stream);

    mutex_unlock (show_b250_mutex);
}

// Write uncompressed, unencrypted section to <section-type>.<vb>.<dict_id>.[header|body]. 
// Note: header includes encryption padding if it was encrypted
static void zfile_dump_section (Buffer *uncompressed_data, SectionHeader *section_header, unsigned section_header_len, DictId dict_id)
{
    char filename[100];
    uint32_t vb_i = BGEN32 (section_header->vblock_i);

    // header
    sprintf (filename, "%s.%u.%s.header", st_name (section_header->section_type), vb_i, dis_dict_id (dict_id).s);
    file_put_data (filename, section_header, section_header_len);

    // body
    if (uncompressed_data->len) {
        sprintf (filename, "%s.%u.%s.body", st_name (section_header->section_type), vb_i, dis_dict_id (dict_id).s);
        buf_dump_to_file (filename, uncompressed_data, 1, false, false, true);
    }
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
    uint8_t param = 0;

    if (expected_section_type == SEC_DICT)
        dict_id = ((SectionHeaderDictionary *)section_header_p)->dict_id;
    else if (expected_section_type == SEC_B250 || expected_section_type == SEC_LOCAL) {
        dict_id = ((SectionHeaderCtx *)section_header_p)->dict_id;
        param   = ((SectionHeaderCtx *)section_header_p)->param;
    }

    if (piz_is_skip_section (vb, expected_section_type, dict_id)) return; // we skip some sections based on flags

    SectionHeader *section_header  = (SectionHeader *)section_header_p;
    uint32_t compressed_offset     = BGEN32 (section_header->compressed_offset);
    uint32_t data_encrypted_len    = BGEN32 (section_header->data_encrypted_len);
    uint32_t data_compressed_len   = BGEN32 (section_header->data_compressed_len);
    uint32_t data_uncompressed_len = BGEN32 (section_header->data_uncompressed_len);
    uint32_t vblock_i              = BGEN32 (section_header->vblock_i);

    // sanity checks
    ASSERTE (section_header->section_type == expected_section_type, "expecting section type %s but seeing %s", st_name(expected_section_type), st_name(section_header->section_type));
    
    ASSERTE (vblock_i == expected_vb_i || !expected_vb_i, // dictionaries are uncompressed by the I/O thread with pseduo_vb (vb_i=0) 
             "bad vblock_i: vblock_i in file=%u but expecting it to be %u (section_type=%s)", 
             vblock_i, expected_vb_i, st_name (expected_section_type));

    // decrypt data (in-place) if needed
    if (data_encrypted_len) 
        crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, vblock_i, section_header->section_type, false);

    if (data_uncompressed_len > 0) { // FORMAT, for example, can be missing in a sample-less file

        if (uncompressed_data_buf_name) {
            buf_alloc (vb, uncompressed_data, data_uncompressed_len + sizeof (uint64_t), 1.1, uncompressed_data_buf_name); // add a 64b word for safety in case this buffer will be converted to a bitarray later
            uncompressed_data->len = data_uncompressed_len;
        }

        comp_uncompress (vb, section_header->codec, section_header->sub_codec, param,
                         (char*)section_header + compressed_offset, data_compressed_len, 
                         uncompressed_data, data_uncompressed_len);
    }
 
    if (flag.show_b250 && expected_section_type == SEC_B250) 
        zfile_show_b250_section (section_header_p, uncompressed_data);
    
    if (flag.dump_section && !strcmp (st_name (expected_section_type), flag.dump_section)) {
        uint64_t save_len = uncompressed_data->len;
        uncompressed_data->len = data_uncompressed_len; // might be different, eg in the case of ref_hash
        zfile_dump_section (uncompressed_data, section_header, compressed_offset, dict_id);
        uncompressed_data->len = save_len; // restore
    }

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
    ASSERTE0 (strlen (metadata) < FILE_METADATA_LEN, "metadata too long");
}

uint32_t zfile_compress_b250_data (VBlock *vb, Context *ctx)
{
    struct FlagsCtx flags = ctx->flags; // make a copy
    flags.paired = ctx->pair_b250;

    SectionHeaderCtx header = (SectionHeaderCtx) { 
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_B250,
        .h.data_uncompressed_len = BGEN32 (ctx->b250.len),
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCtx)),
        .h.codec                 = ctx->bcodec == CODEC_UNKNOWN ? CODEC_BZ2 : ctx->bcodec,
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .h.flags.ctx             = flags,
        .dict_id                 = ctx->dict_id,
        .ltype                   = ctx->ltype
    };

    return comp_compress (vb, &vb->z_data, (SectionHeader*)&header, ctx->b250.data, NULL);
}

LocalGetLineCB *zfile_get_local_data_callback (DataType dt, Context *ctx)
{
    static struct { DataType dt; const uint64_t *dict_id_num; LocalGetLineCB *func; } callbacks[] = LOCAL_GET_LINE_CALLBACKS;

    for (unsigned i=0; i < sizeof(callbacks)/sizeof(callbacks[0]); i++)
        if (callbacks[i].dt == dt && *callbacks[i].dict_id_num == ctx->dict_id.num && !ctx->no_callback) 
            return callbacks[i].func;

    return NULL;
}

// returns compressed size
uint32_t zfile_compress_local_data (VBlock *vb, Context *ctx, uint32_t sample_size /* 0 means entire local buffer */)
{   
    vb->has_non_agct = false;
    
    struct FlagsCtx flags = ctx->flags; // make a copy
    flags.paired     = ctx->pair_local;
    flags.copy_param = ctx->local_param;
                    
    uint32_t uncompressed_len = ctx->local.len * lt_desc[ctx->ltype].width;
    
    // case: we're just testing a small sample
    if (sample_size && uncompressed_len > sample_size) 
        uncompressed_len = sample_size;

    uint8_t unused_bits = 0;
    if (ctx->ltype == LT_BITMAP) {
        BitArray *bm = buf_get_bitarray (&ctx->local);
        unused_bits = ((uint8_t)64 - (uint8_t)(bm->nbits % 64)) % (uint8_t)64;
    }

    SectionHeaderCtx header = (SectionHeaderCtx) {
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_LOCAL,
        .h.data_uncompressed_len = BGEN32 (uncompressed_len),
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCtx)),
        .h.codec                 = ctx->lcodec == CODEC_UNKNOWN ? CODEC_BZ2 : ctx->lcodec, // if codec has not been decided yet, fall back on BZ2
        .h.sub_codec             = ctx->lsubcodec_piz ? ctx->lsubcodec_piz : codec_args[ctx->lcodec].sub_codec,
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .h.flags.ctx             = flags,
        .dict_id                 = ctx->dict_id,
        .ltype                   = ctx->ltype,
        .param                   = flags.copy_param ? (uint8_t)ctx->local.param : unused_bits,
    };

    LocalGetLineCB *callback = zfile_get_local_data_callback (vb->data_type, ctx);

    return comp_compress (vb, &vb->z_data, (SectionHeader*)&header, 
                          callback ? NULL : ctx->local.data, 
                          callback);
}

// compress section - two options for input data - 
// 1. contiguous data in section_data 
// 2. line by line data - by providing a callback + total_len
void zfile_compress_section_data_ex (VBlock *vb, SectionType section_type, 
                                     Buffer *section_data,          // option 1 - compress contiguous data
                                     LocalGetLineCB callback, uint32_t total_len, // option 2 - compress data one line at a time
                                     Codec codec, SectionFlags flags)
{
    SectionHeader header = (SectionHeader){
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .section_type          = section_type,
        .data_uncompressed_len = BGEN32 (section_data ? section_data->len : total_len),
        .compressed_offset     = BGEN32 (sizeof(header)),
        .codec                 = codec,
        .vblock_i              = BGEN32 (vb->vblock_i),
        .flags                 = flags
    };

    if (flag.show_time) codec_show_time (vb, st_name (section_type), NULL, codec);

    comp_compress (vb, &vb->z_data, &header, 
                   section_data ? section_data->data : NULL, 
                   callback);
}

// reads exactly the length required, error otherwise. manages read buffers to optimize I/O performance.
// this doesn't make a big difference for SSD, but makes a huge difference for HD
// return a pointer to the data read
static void *zfile_read_from_disk (File *file, VBlock *vb, Buffer *buf, uint32_t len, SectionType st)
{
    START_TIMER;

    ASSERTE (len, "reading %s: len is 0", st_name (st));
    ASSERTE (buf_has_space (buf, len), "reading %s: buf is out of space: len=%u but remaining space in buffer=%u (tip: run with --show-headers to see where it fails)",
            st_name (st), len, (uint32_t)(buf->size - buf->len));

    char *start = AFTERENT (char, *buf);
    uint32_t bytes = fread (start, 1, len, (FILE *)file->file);
    ASSERTE (bytes == len, "reading %s: read only %u bytes out of len=%u", st_name (st), bytes, len);

    buf->len += bytes;
    file->disk_so_far += bytes;

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
    ASSERTE (!sl || expected_sec_type == sl->section_type, "expected_sec_type=%s but encountered sl->section_type=%s. vb_i=%u",
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

    ASSERTE (header, "Failed to read data from file %s while expecting section type %s: %s", 
             z_name, st_name(expected_sec_type), strerror (errno));
    
    bool is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;

    // SEC_REFERENCE is never encrypted when originating from a reference file, it is encrypted (if the file is encrypted) if it originates from REF_INTERNAL 
    if (is_encrypted && header->section_type == SEC_REFERENCE && !header->data_encrypted_len) {
        is_encrypted = false;
        header_size  = unencrypted_header_size;
    }

    // decrypt header (note: except for SEC_GENOZIP_HEADER - this header is never encrypted)
    if (is_encrypted) {
        ASSINP (BGEN32 (header->magic) != GENOZIP_MAGIC, 
                "password provided, but file %s is not encrypted (sec_type=%s)", z_name, st_name (header->section_type));

        crypt_do (vb, (uint8_t*)header, header_size, original_vb_i, expected_sec_type, true); 
    
        is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC; // update after decryption
    }

    if (flag.show_headers) {
        zfile_show_header (header, NULL, sl ? sl->offset : 0, 'R');
        if (exe_type == EXE_GENOCAT && (expected_sec_type == SEC_B250 || expected_sec_type == SEC_LOCAL || expected_sec_type == SEC_DICT))
             return header_offset; // in genocat --show-header - we only show headers, nothing else
    }

    ASSERTE (is_magical, "corrupt data (magic is wrong) when attempting to read section=%s dict_id=%s of vblock_i=%u component=%u in file %s", 
            st_name (expected_sec_type), sl ? dis_dict_id (sl->dict_id).s : "(no sl)", vb->vblock_i, z_file->num_txt_components_so_far, z_name);

    uint32_t compressed_offset   = BGEN32 (header->compressed_offset);
    ASSERTE (compressed_offset, "header.compressed_offset is 0 when reading section_type=%s", st_name(expected_sec_type));

    uint32_t data_compressed_len = BGEN32 (header->data_compressed_len);
    uint32_t data_encrypted_len  = BGEN32 (header->data_encrypted_len);

    uint32_t data_len = MAX (data_compressed_len, data_encrypted_len);

    // in case where we already read part of the body (eg if is_encrypted was initially set and then unset) (remaining_data_len might be negative)
    int32_t remaining_data_len = (int32_t)data_len - (int32_t)(bytes_read - header_size); 
    
    // check that we received the section type we expect, 
    ASSERTE (expected_sec_type == header->section_type || 
             (expected_sec_type == SEC_GENOZIP_HEADER && (SectionType)header->sub_codec == SEC_GENOZIP_HEADER), // in v2-5, the section_type field was located where sub_codec is now
             "Unexpected section type when reading %s: expecting %s, found %s sl(expecting)=(offset=%s, dict_id=%s)",
             z_name, st_name(expected_sec_type), st_name(header->section_type), 
             sl ? str_uint_commas (sl->offset).s : "N/A", sl ? dis_dict_id (sl->dict_id).s : "N/A");

    ASSERTE (compressed_offset == header_size || expected_sec_type == SEC_GENOZIP_HEADER, // we allow SEC_GENOZIP_HEADER of other sizes, for older versions
             "invalid header when reading %s - expecting compressed_offset to be %u but found %u. section_type=%s", 
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
SectionHeader *zfile_read_section_header (VBlockP vb, uint64_t offset, 
                                          uint32_t original_vb_i, // the vblock_i used for compressing. this is part of the encryption key. dictionaries are compressed by the compute thread/vb, but uncompressed by the I/O thread (vb=0)
                                          SectionType expected_sec_type)
{
    uint32_t header_size = st_header_size (expected_sec_type);

    // get the uncompressed size from one of the headers - they are all the same size, and the reference file is never encrypted
    file_seek (z_file, offset, SEEK_SET, false);

    bool is_encrypted =  (z_file->data_type != DT_REF) && 
                         (expected_sec_type != SEC_GENOZIP_HEADER) &&
                         crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
        
    ASSERTE (!vb->compressed.len, "vb_i=%u expected_sec_type=%s: expecting vb->compressed to be free, but it's not: %s",
             vb->vblock_i, st_name (expected_sec_type), buf_desc (&vb->compressed).s);
    
    buf_alloc (vb, &vb->compressed, header_size, 4, "compressed"); 

    SectionHeader *header = zfile_read_from_disk (z_file, evb, &vb->compressed, header_size, expected_sec_type); 

    ASSERTE (header, "Failed to read header of section type %s from file %s: %s", 
             st_name(expected_sec_type), z_name, strerror (errno));

    bool is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;

    // decrypt header 
    if (is_encrypted) {
        ASSERTE (BGEN32 (header->magic) != GENOZIP_MAGIC, 
                 "password provided, but file %s is not encrypted (sec_type=%s)", z_name, st_name (header->section_type));

        crypt_do (vb, (uint8_t*)header, header_size, original_vb_i, expected_sec_type, true); 
    
        is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC; // update after decryption
    }

    ASSERTE (is_magical, "corrupt data (magic is wrong) when attempting to read header of section %s in file %s", 
             st_name (expected_sec_type), z_name);

    return header;
}

static void zfile_read_genozip_header_handle_ref_info (const SectionHeaderGenozipHeader *header)
{
    if (flag.show_reference && !md5_is_zero (header->ref_file_md5)) {
        iprintf ("%s was compressed using the reference file:\nName: %s\nMD5: %s\n",
                    z_name, header->ref_filename, digest_display (header->ref_file_md5).s);
        if (exe_type == EXE_GENOCAT) exit_ok; // in genocat --show-reference, we only show the reference, not the data
    }

    if (flag.reading_reference)
        ref_file_md5 = header->ref_file_md5;

    else if (!md5_is_zero (header->ref_file_md5) && // reference file specified in the header
             exe_type != EXE_GENOLS) {              // we don't need the reference for genols

        // case: this file requires an external reference, but command line doesn't include --reference - attempt to use the
        // reference specified in the header. 
        // Note: this code will be executed when zfile_read_genozip_header is called from main_genounzip.
        if (!flag.explicit_ref && // reference NOT was specified on command line
            !(ref_filename && !strcmp (ref_filename, header->ref_filename))) { // ref_filename already set from a previous file with the same reference
            
            if (file_exists (header->ref_filename)) {
                ASSERTW (flag.genocat_no_reconstruct, "Note: using the reference file %s. You can override this with --reference", header->ref_filename);
                ref_set_reference (header->ref_filename, REF_EXTERNAL, false);
            }
            else 
                ASSINP (flag.genocat_no_ref_file, "%s: please use --reference specify the current path to reference file with which %s was compressed (original path was %s)",
                        global_cmd, z_name, header->ref_filename);
        }

        // test for matching MD5 between specified external reference and reference in the header
        ASSERTE (!flag.explicit_ref || (digest_is_zero (ref_file_md5) || digest_is_equal (header->ref_file_md5, ref_file_md5)),
                 "%s: Bad reference file:\n%s (MD5=%s) was used for compressing\n%s (MD5=%s) has a different MD5",
                 z_name, header->ref_filename, digest_display (header->ref_file_md5).s, ref_filename, digest_display (ref_file_md5).s);
    }
}

// returns false if file should be skipped
bool zfile_read_genozip_header (Digest *digest, uint64_t *txt_data_size, uint64_t *num_items_bound, char *created) // optional outs
{
    // read the footer from the end of the file
    if (file_get_size (z_file->name) < sizeof(SectionFooterGenozipHeader) ||
        !file_seek (z_file, -sizeof(SectionFooterGenozipHeader), SEEK_END, 2))
        goto error;

    SectionFooterGenozipHeader footer;
    int ret = fread (&footer, sizeof (footer), 1, (FILE *)z_file->file);
    ASSERTW (ret == 1, "Skipping empty file %s", z_name);
    if (!ret) goto error;
    
    // case: there is no genozip header. this can happen if the file was truncated (eg because compression did not complete)
    // note: this can also happen if the file is genozip v1, but I don't think there are any real v1 files in the wild
    // so I will keep the error message simple and not mention it
    ASSERTGOTO (BGEN32 (footer.magic) == GENOZIP_MAGIC, "Error: failed to read file %s - the file appears to be incomplete (it is missing the genozip footer).", z_name);

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

    ASSERTGOTO (header->genozip_version <= GENOZIP_FILE_FORMAT_VERSION, 
                "Error: %s cannot be opened because it was compressed with a newer version of genozip (version %u.x.x) while the version you're running is older (version %s).\n"
                "You might want to consider upgrading genozip to the newest version.\n",
                z_name, header->genozip_version, GENOZIP_CODE_VERSION);

    // in version 6, we canceled backward compatability with v1-v5
    ASSERTGOTO (header->genozip_version >= 6, "Skipping %s: it was compressed with an older version of genozip - version %u.\nIt may be uncompressed with genozip versions %u to 5",
                z_name, header->genozip_version, header->genozip_version);

    // in version 7, we canceled backward compatability with v6
    ASSERTGOTO (header->genozip_version >= 7, "Skipping %s: it was compressed with version 6 of genozip. It may be uncompressed with genozip version 6",
                z_name);

    // in version 8, we canceled backward compatability with v7
    ASSERTGOTO (header->genozip_version >= 8, "Skipping %s: it was compressed with version 7 of genozip. It may be uncompressed with genozip version 7",
                z_name);

    DataType data_type = (DataType)(BGEN16 (header->data_type)); 
    ASSERTE ((unsigned)data_type < NUM_DATATYPES, "unrecognized data_type=%d: please upgrade genozip to the latest version", data_type);

    if (z_file->data_type == DT_NONE || z_file->data_type == DT_GENERIC) {
        z_file->data_type = data_type;
        z_file->type      = file_get_z_ft_by_dt (z_file->data_type);  
    }
    else
        ASSINP (z_file->data_type == data_type, "%s - file extension indicates this is a %s file, but according to its contents it is a %s", 
                z_name, dt_name (z_file->data_type), dt_name (data_type));

    if (txt_file && header->h.flags.genozip_header.txt_is_bin)// txt_file is still NULL when called from main_genozip
        txt_file->data_type = DTPZ (bin_type);

    ASSINP (header->encryption_type != ENC_NONE || !crypt_have_password() || z_file->data_type == DT_REF, 
            "password provided, but file %s is not encrypted", z_name);

    ASSERTE (BGEN32 (header->h.compressed_offset) == sizeof (SectionHeaderGenozipHeader),
             "invalid genozip header - expecting compressed_offset to be %u but found %u", 
             (unsigned)sizeof (SectionHeaderGenozipHeader), BGEN32 (header->h.compressed_offset));

    // get & test password, if file is encrypted
    if (header->encryption_type != ENC_NONE) {

        if (!crypt_have_password()) crypt_prompt_for_password();

        crypt_do (evb, header->password_test, sizeof(header->password_test), 0, SEC_NONE, true); // decrypt password test

        ASSINP (!memcmp (header->password_test, password_test_string, sizeof(header->password_test)),
                "password is wrong for file %s", z_name);
    }

    z_file->num_components    = BGEN32 (header->num_components);
    if (z_file->num_components < 2) flag.unbind = 0; // override user's --unbind if file has only 1 component

    z_file->genozip_version   = header->genozip_version;
    
    int dts = z_file->z_flags.dt_specific; // save in case its set already (eg dts_paired is set in fastq_piz_is_paired)
    z_file->z_flags = header->h.flags.genozip_header;
    z_file->z_flags.dt_specific |= dts;

    if (digest) *digest       = header->digest_bound; 
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
        ASSINP (data_type == DT_REF, "Error: %s is not a reference file. To create a reference file, use 'genozip --make-reference <fasta-file.fa>'",
                ref_filename);

        ref_set_ref_file_info (header->digest_bound, header->ref_filename); // in the reference file itself, header->ref_filename is the original fasta used to create this reference
    }

    // case: we are reading a file that is not expected to be a reference file
    else {
        // case: we are attempting to decompress a reference file - this is not supported
        ASSERTGOTO (data_type != DT_REF || (flag.genocat_no_reconstruct && exe_type == EXE_GENOCAT) || exe_type == EXE_GENOLS,
                    "%s is a reference file - it cannot be decompressed. Skipping it.", z_name);

        // handle reference file info
        zfile_read_genozip_header_handle_ref_info (header);
    }
     
    buf_free (&evb->z_data);
    return true;

error: 
    ASSERTE0 (!flag.reading_reference, "failed to read reference file");
    buf_free (&evb->z_data);
    return false;
}

void zfile_compress_genozip_header (Digest single_component_digest)
{
    SectionHeaderGenozipHeader header = {};

    // start with just the fields needed by sections_add_to_list
    header.h.section_type = SEC_GENOZIP_HEADER;

    // "manually" add the genozip section to the section list - normally it is added in comp_compress()
    // but in this case the genozip section containing the list will already be ready...
    sections_add_to_list (evb, &header.h);
    sections_list_concat (evb); // usually done in output_process_vb, but the section list will be already compressed within the genozip header...

    bool is_encrypted = crypt_have_password();

    uint32_t num_sections = z_file->section_list_buf.len;

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.compressed_offset     = BGEN32 (sizeof (SectionHeaderGenozipHeader));
    header.h.data_uncompressed_len = BGEN32 (z_file->section_list_buf.len * sizeof (SectionListEntry));
    header.h.codec                 = CODEC_BZ2;
    header.h.flags.genozip_header  = (struct FlagsGenozipHeader) {
        .txt_is_bin   = DTPT (is_binary),
        .dt_specific  = DT_FUNC (z_file, zip_dts_flag)(),
        .aligner      = (flag.ref_use_aligner > 0),
        .bgzf         = (txt_file->codec == CODEC_BGZF),
        .adler        = !flag.md5,
        .dual_coords  = z_file->z_flags.dual_coords
    };
    header.genozip_version         = GENOZIP_FILE_FORMAT_VERSION;
    header.data_type               = BGEN16 ((uint16_t)dt_get_txt_dt (z_file->data_type));
    header.encryption_type         = is_encrypted ? ENC_AES256 : ENC_NONE;
    header.uncompressed_data_size  = BGEN64 (z_file->txt_data_so_far_bind);
    header.num_items_bound         = BGEN64 (z_file->num_lines);
    header.num_sections            = BGEN32 (num_sections); 
    header.num_components          = BGEN32 (z_file->num_txt_components_so_far);
    
    // when decompressing will require an external reference, we set header.ref_filename to the name of the genozip reference file
    if (flag.reference == REF_EXTERNAL) {   
        strncpy (header.ref_filename, ref_filename, REF_FILENAME_LEN-1);
        header.ref_file_md5 = ref_file_md5;
    }

    // in --make-ref, we set header.ref_filename to the original fasta file, to be used later in ref_get_cram_ref
    // (unless the fasta is piped from stdin, or its name is too long)
    else if (flag.make_reference && strcmp (txt_name, FILENAME_STDIN) && strlen (txt_name) <= REF_FILENAME_LEN-1) {
#ifndef WIN32
        char *ref_filename = realpath (txt_name, NULL); // allocates memory
        ASSERTE (ref_filename, "realpath() failed: %s", strerror (errno));

        strncpy (header.ref_filename, ref_filename, REF_FILENAME_LEN-1);
        FREE (ref_filename);
#else
        ASSERTE0 (_fullpath (header.ref_filename, txt_name, REF_FILENAME_LEN), "_fullpath() failed");
#endif
    }

    uint32_t license_num_bgen = BGEN32 (license_get());
    header.license_hash = md5_do (&license_num_bgen, sizeof (int32_t));
    
    header.digest_bound = flag.data_modified ? DIGEST_NONE
                        : flag.bind          ? digest_snapshot (&z_file->digest_ctx_bound)
                        :                      single_component_digest; // if not in bound mode - just copy the md5 of the single file

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
    comp_compress (evb, z_data, &header.h, z_file->section_list_buf.data, NULL);

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

    zfile_output_processed_vb (evb); // write footer
}

// ZIP
void zfile_write_txt_header (Buffer *txt_header, 
                             uint64_t unmodified_txt_header_len, // length of header before modifications, eg due to --chain or compressing a Laft file
                             Digest header_md5, bool is_first_txt)
{
    SectionHeaderTxtHeader header = {
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_TXT_HEADER,
        .h.data_uncompressed_len = BGEN32 (txt_header->len),
        .h.compressed_offset     = BGEN32 (sizeof (SectionHeaderTxtHeader)),
        .h.codec                 = CODEC_BZ2,
        .h.flags.txt_header.liftover_rejects = flag.processing_rejects,
        .codec                   = txt_file->codec, 
        .digest_header           = flag.data_modified ? DIGEST_NONE : header_md5,
    };

    // In BGZF, we store the 3 least significant bytes of the file size, so check if the reconstructed BGZF file is likely the same
    if (txt_file->codec == CODEC_BGZF) 
        bgzf_sign (txt_file->disk_size, header.codec_info);
        
    file_basename (txt_file->name, false, FILENAME_STDIN, header.txt_filename, TXT_FILENAME_LEN);
    file_remove_codec_ext (header.txt_filename, txt_file->type); // eg "xx.fastq.gz -> xx.fastq"
    
    static Buffer txt_header_buf = EMPTY_BUFFER;

    buf_alloc (evb, &txt_header_buf, sizeof (SectionHeaderTxtHeader) + txt_header->len / 3, // generous guess of compressed size
               1, "txt_header_buf"); 

    comp_compress (evb, &txt_header_buf, (SectionHeader*)&header, 
                   txt_header->len ? txt_header->data : NULL, // actual header may be missing (eg in SAM it is permitted to not have a header)
                   NULL);

    file_write (z_file, txt_header_buf.data, txt_header_buf.len);

    z_file->disk_so_far += txt_header_buf.len;   // length of GENOZIP data writen to disk

    // note: the liftover reject txt data is of course not counted as part of the file txt data for stats...
    if (!flag.processing_rejects) {        
        z_file->txt_data_so_far_single   += txt_header->len; // length of txt header as it would be reconstructed (possibly afer modifications)
        z_file->txt_data_so_far_bind     += txt_header->len;
        z_file->txt_data_so_far_single_0 += unmodified_txt_header_len; // length of the original txt header as read from the file
        z_file->txt_data_so_far_bind_0   += unmodified_txt_header_len;
        z_file->codec                     = txt_file->codec; // for stats (stats can't use txt_file.codec as it might be the rejects file)
    }

    buf_free (&txt_header_buf); 

    // copy it to z_file - we might need to update it at the very end in zfile_update_txt_header_section_header()
    if (is_first_txt)
        memcpy (&z_file->txt_header_first, &header, sizeof (header));

    memcpy (&z_file->txt_header_single, &header, sizeof (header));
}

// Update SEC_TXT_HEADER. If we're compressing a plain file, we will know
// the bytes upfront, but if we're binding or compressing a eg .GZ, we will need to update it
// when we're done. num_lines can only be known after we're done with this txt component.
// if we cannot update the header - that's fine, these fields are only used for the progress indicator on --list
bool zfile_update_txt_header_section_header (uint64_t offset_in_z_file, uint32_t max_lines_per_vb,
                                             Digest *md5 /* out */)
{
    // rewind to beginning of current (latest) vcf header - nothing to do if we can't
    if (!file_seek (z_file, offset_in_z_file, SEEK_SET, true)) return false;

    uint32_t len = crypt_padded_len (sizeof (SectionHeaderTxtHeader));

    // sanity check - we skip empty files, so data is expected
    ASSERTE (txt_file->txt_data_size_single > 0, "Expecting txt_file->txt_data_size_single=%"PRId64" > 0", txt_file->txt_data_size_single);
    
    // update the header of the single (current) vcf. 
    SectionHeaderTxtHeader *curr_header = &z_file->txt_header_single;
    curr_header->txt_data_size    = BGEN64 (txt_file->txt_data_size_single);
    curr_header->num_lines        = BGEN64 (txt_file->num_lines);
    curr_header->max_lines_per_vb = BGEN32 (max_lines_per_vb);
    curr_header->digest_single    = flag.data_modified ? DIGEST_NONE 
                                                       : digest_finalize (&z_file->digest_ctx_single, "component:digest_ctx_single");

    *md5 = curr_header->digest_single;

    if (offset_in_z_file == 0) 
        z_file->txt_header_first.digest_single = curr_header->digest_single; // first vcf - update the stored header 

    // encrypt if needed
    if (crypt_have_password()) 
        crypt_do (evb, (uint8_t *)curr_header, len, 0, curr_header->h.section_type, true);

    file_write (z_file, curr_header, len);
    fflush ((FILE*)z_file->file); // its not clear why, but without this fflush the bytes immediately after the first header get corrupted (at least on Windows with gcc)
    
    file_seek (z_file, 0, SEEK_END, false); // return to the end of the file

    if (flag.show_headers)
        zfile_show_header ((SectionHeader *)curr_header, NULL, offset_in_z_file, 'W'); 

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
        .h.flags.vb_header   = vb->vb_header_flags,
        .num_lines           = BGEN32 ((uint32_t)vb->lines.len),
        .vb_data_size        = BGEN32 (vb->vb_data_size),
        .longest_line_len    = BGEN32 (vb->longest_line_len),
        .digest_so_far       = flag.data_modified ? DIGEST_NONE : vb->digest_so_far
    };

    // copy section header into z_data - to be eventually written to disk by the I/O thread. this section doesn't have data.
    comp_compress (vb, &vb->z_data, (SectionHeader*)&vb_header, NULL, NULL);
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
        iprintf ("vb_i=%u component=%u first_line=%u num_lines=%u txt_file=%u genozip_size=%u longest_line_len=%u\n",
                 vb->vblock_i, z_file->num_txt_components_so_far, txt_first_line_i, BGEN32 (vb_header->num_lines), 
                 BGEN32 (vb_header->vb_data_size), BGEN32 (vb_header->z_data_bytes), 
                 BGEN32 (vb_header->longest_line_len));

    // now we can finally encrypt the header - if needed
    if (crypt_have_password())  
        crypt_do (vb, (uint8_t*)vb_header, BGEN32 (vb_header->h.compressed_offset),
                  BGEN32 (vb_header->h.vblock_i), vb_header->h.section_type, true);
}

void zfile_output_processed_vb (VBlock *vb)
{
    START_TIMER;

    sections_list_concat (vb);
    
    file_write (z_file, vb->z_data.data, vb->z_data.len);
    COPY_TIMER (write);

    z_file->disk_so_far += (int64_t)vb->z_data.len;
    vb->z_data.len = 0;

    // this function holds the mutex and hence has a non-trival performance penalty. we call
    // it only if the user specifically requested --show-stats
    if (flag.show_stats) ctx_update_stats (vb);

    if (flag.show_headers && buf_is_allocated (&vb->show_headers_buf))
        buf_print (&vb->show_headers_buf, false);

    if (flag.show_threads) dispatcher_show_time ("Write genozip data done", -1, vb->vblock_i);
}

DataType zfile_get_file_dt (const char *filename)
{
    FileType ft = file_get_type (filename);
    DataType dt = file_get_dt_by_z_ft (ft);
    File *file = NULL;

    // case: we don't know yet what file type this is - we need to read the genozip header to determine
    if (dt == DT_NONE && filename) {
        if (!(file = file_open (filename, READ, Z_FILE, DT_NONE)) || !file->file)
            goto done; // not a genozip file

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
    }

done:
    file_close (&file, false, false);
    return dt;
}

