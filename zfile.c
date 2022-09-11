// ------------------------------------------------------------------
//   zfile.c
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <errno.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include "vblock.h"
#include "zfile.h"
#include "crypt.h"
#include "context.h"
#include "file.h"
#include "endianness.h"
#include "version.h"
#include "sections.h"
#include "compressor.h"
#include "codec.h"
#include "piz.h"
#include "license.h"
#include "mutex.h"
#include "strings.h"
#include "dict_id.h"
#include "reference.h"
#include "bgzf.h"
#include "digest.h"
#include "md5.h"
#include "flags.h"
#include "website.h"
#include "gencomp.h"
#include "libdeflate/libdeflate.h"

static rom password_test_string = "WhenIThinkBackOnAllTheCrapIlearntInHighschool";

static void zfile_show_b250_section (void *section_header_p, ConstBufferP b250_data)
{
    static Mutex show_b250_mutex = {}; // protect so compute thread's outputs don't get mix

    SectionHeaderCtx *header = (SectionHeaderCtx *)section_header_p;

    if (!flag.show_b250 && dict_id_typeless (header->dict_id).num != flag.dict_id_show_one_b250.num) return;

    mutex_initialize (show_b250_mutex); // possible unlikely race condition on initializing - good enough for debugging purposes
    mutex_lock (show_b250_mutex);

    iprintf ("vb_i=%u %*.*s: ", BGEN32 (header->h.vblock_i), -DICT_ID_LEN-1, DICT_ID_LEN, dict_id_typeless (header->dict_id).id);

    bytes data  = B1ST (const uint8_t, *b250_data);
    bytes after = BAFT (const uint8_t, *b250_data);

    while (data < after) {
        WordIndex word_index = ctx_decode_b250 (&data, true, header->b250_size, "zfile_show_b250_section");
        switch (word_index) {
            case WORD_INDEX_ONE_UP  : iprint0 ("ONE_UP " ) ; break ;
            case WORD_INDEX_EMPTY   : iprint0 ("EMPTY "  ) ; break ;
            case WORD_INDEX_MISSING : iprint0 ("MISSING ") ; break ;
            default                 : iprintf ("%u ", word_index);
        }
    }
    iprint0 ("\n");

    mutex_unlock (show_b250_mutex);
}

// Write uncompressed, unencrypted section to <section-type>.<vb>.<dict_id>.[header|body]. 
// Note: header includes encryption padding if it was encrypted
static void zfile_dump_section (BufferP uncompressed_data, SectionHeader *section_header, unsigned section_header_len, DictId dict_id)
{
    char filename[100];
    VBIType vb_i = BGEN32 (section_header->vblock_i);

    // header
    sprintf (filename, "%s.%u.%s.header", st_name (section_header->section_type), vb_i, dis_dict_id (dict_id).s);
    file_put_data (filename, section_header, section_header_len, 0);

    // body
    if (uncompressed_data->len) {
        sprintf (filename, "%s.%u.%s.body", st_name (section_header->section_type), vb_i, dis_dict_id (dict_id).s);
        buf_dump_to_file (filename, uncompressed_data, 1, false, false, true, false);
    }
}

// uncompressed a block and adds a \0 at its end. Returns the length of the uncompressed block, without the \0.
// when we get here, the header is already unencrypted zfile_one_section
void zfile_uncompress_section (VBlockP vb,
                               void *section_header_p,
                               BufferP uncompressed_data, 
                               rom uncompressed_data_buf_name, // a name if Buffer, NULL ok if buffer need not be realloced
                               uint32_t expected_vb_i,
                               SectionType expected_section_type) 
{
    START_TIMER;
    ASSERTNOTNULL (section_header_p);
    
    DictId dict_id = DICT_ID_NONE;
    uint8_t codec_param = 0;

    if (expected_section_type == SEC_DICT)
        dict_id = ((SectionHeaderDictionary *)section_header_p)->dict_id;
    else if (expected_section_type == SEC_B250 || expected_section_type == SEC_LOCAL) {
        dict_id = ((SectionHeaderCtx *)section_header_p)->dict_id;
        codec_param = ((SectionHeaderCtx *)section_header_p)->param; 
    }
    else if (expected_section_type == SEC_COUNTS) 
        dict_id = ((SectionHeaderCounts *)section_header_p)->dict_id;

    ContextP ctx = NULL;
    if (IS_DICTED_SEC (expected_section_type)) { 
        ctx = ECTX(dict_id); 
        if (ctx && (!ctx->is_loaded && flag.pair!=2)) // pair=2 happens in ZIP of paired FASTQ when reading pair 1 z_file data 
            return;  // section was skipped 
    }
    else 
        if (piz_is_skip_undicted_section (expected_section_type)) return; // section was skipped 

    SectionHeader *section_header  = (SectionHeader *)section_header_p;
    uint32_t compressed_offset     = BGEN32 (section_header->compressed_offset);
    uint32_t data_encrypted_len    = BGEN32 (section_header->data_encrypted_len);
    uint32_t data_compressed_len   = BGEN32 (section_header->data_compressed_len);
    uint32_t data_uncompressed_len = BGEN32 (section_header->data_uncompressed_len);
    VBIType vblock_i               = BGEN32 (section_header->vblock_i);

    // sanity checks
    ASSERT (section_header->section_type == expected_section_type, "expecting section type %s but seeing %s", st_name(expected_section_type), st_name(section_header->section_type));
    
    ASSERT (vblock_i == expected_vb_i || !expected_vb_i, // dictionaries are uncompressed by the main thread with pseduo_vb (vb_i=0) 
            "bad vblock_i: section_header->vblock_i=%u but expecting it to be %u (section_type=%s dict_id=%s)", 
            vblock_i, expected_vb_i, st_name (expected_section_type), dis_dict_id(dict_id).s);

    if (flag.show_uncompress)
        iprintf ("Uncompress: %s %-9s %-8s comp_len=%-7u uncomp_len=%u\n", VB_NAME, 
                 st_name (expected_section_type), dict_id.num ? dis_dict_id (dict_id).s : "", data_compressed_len, data_uncompressed_len);
        
    // decrypt data (in-place) if needed
    if (data_encrypted_len) 
        crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, vblock_i, section_header->section_type, false);
 
    bool bad_compression = false;
    if (data_uncompressed_len > 0) { // FORMAT, for example, can be missing in a sample-less file

        if (uncompressed_data_buf_name) {
            buf_alloc (vb, uncompressed_data, 0, data_uncompressed_len + sizeof (uint64_t), char, 1.1, uncompressed_data_buf_name); // add a 64b word for safety in case this buffer will be converted to a bitarray later
            uncompressed_data->len = data_uncompressed_len;
        }

        comp_uncompress (vb, section_header->codec, section_header->sub_codec, codec_param,
                         (char*)section_header + compressed_offset, data_compressed_len, 
                         uncompressed_data, data_uncompressed_len,
                         dict_id.num ? dis_dict_id(dict_id).s : st_name(expected_section_type));

        //--verify-codec: verify that adler32 of the decompressed data is equal that of the original uncompressed data
        if (flag.verify_codec && uncompressed_data && data_uncompressed_len && 
            BGEN32 (section_header->magic) != GENOZIP_MAGIC &&
            section_header->uncomp_adler32 != adler32 (1, uncompressed_data->data, data_uncompressed_len)) {
        
            iprintf ("--verify-codec: BAD ADLER32 section decompressed incorrectly: codec=%s\n", codec_name(section_header->codec));
            sections_show_header (section_header, NULL, 0, 'R');
            bad_compression = true;
        }
    }
 
    if (flag.show_b250 && expected_section_type == SEC_B250) 
        zfile_show_b250_section (section_header_p, uncompressed_data);
    
    if ((flag.dump_section && !strcmp (st_name (expected_section_type), flag.dump_section)) || bad_compression) {
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

    rom user = getenv ("USER");
    rom host = getenv ("HOSTNAME");
    sprintf(metadata, "%.19s %.20s%s%.30s", time_buf,  // 71 chars + the string termintor \0 = 72
            user ? user : "",
            user && host ? "@" : "",
            host ? host : "");

    // sanity
    ASSERT0 (strlen (metadata) < FILE_METADATA_LEN, "metadata too long");
}

uint32_t zfile_compress_b250_data (VBlockP vb, ContextP ctx)
{
    struct FlagsCtx flags = ctx->flags; // make a copy
    
    if (VB_DT(DT_FASTQ))
        flags.paired = ctx->pair_b250;

    SectionHeaderCtx header = (SectionHeaderCtx) { 
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_B250,
        .h.data_uncompressed_len = BGEN32 (ctx->b250.len32),
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCtx)),
        .h.codec                 = ctx->bcodec == CODEC_UNKNOWN ? CODEC_RANS8 : ctx->bcodec,
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .h.flags.ctx             = flags,
        .dict_id                 = ctx->dict_id,
        .ltype                   = ctx->ltype,
        .b250_size               = ctx->b250_size
    };

    ctx->b250_in_z = vb->z_data.len32;

    uint32_t compressed_size = comp_compress (vb, ctx, &vb->z_data, (SectionHeader*)&header, ctx->b250.data, NO_CALLBACK, ctx->tag_name);

    ctx->b250_in_z_len = vb->z_data.len32 - ctx->b250_in_z;

    return compressed_size;
}

LocalGetLineCB *zfile_get_local_data_callback (DataType dt, ContextP ctx)
{
    static struct { DataType dt; const uint64_t dict_id_num; LocalGetLineCB *func; } callbacks[] = LOCAL_GET_LINE_CALLBACKS;

    if (ctx && !ctx->no_callback)
        for (unsigned i=0; i < ARRAY_LEN(callbacks); i++)
            if (callbacks[i].dt == dt && callbacks[i].dict_id_num == ctx->dict_id.num) 
                return callbacks[i].func;

    return NULL;
}

// returns compressed size
uint32_t zfile_compress_local_data (VBlockP vb, ContextP ctx, uint32_t sample_size /* 0 means entire local buffer */)
{   
    struct FlagsCtx flags = ctx->flags; // make a copy

    if (VB_DT(DT_FASTQ))
        flags.paired = ctx->pair_local;

    uint32_t uncompressed_len = ctx->local.len32 * lt_desc[ctx->ltype].width;
    
    // case: we're just testing a small sample
    if (sample_size && uncompressed_len > sample_size) 
        uncompressed_len = sample_size;

    SectionHeaderCtx header = (SectionHeaderCtx) {
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_LOCAL,
        .h.data_uncompressed_len = BGEN32 (uncompressed_len),
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCtx)),
        .h.codec                 = ctx->lcodec == CODEC_UNKNOWN ? CODEC_RANS8 : ctx->lcodec, // if codec has not been decided yet, fall back on RANS8
        .h.sub_codec             = ctx->lsubcodec_piz ? ctx->lsubcodec_piz : CODEC_UNKNOWN,
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .h.flags.ctx             = flags,
        .dict_id                 = ctx->dict_id,
        .ltype                   = ctx->ltype,
        .param                   = ctx->local_param ? ctx->local.prm8[0] : 0,
    };

    LocalGetLineCB *callback = zfile_get_local_data_callback (vb->data_type, ctx);

    ctx->local_in_z = vb->z_data.len;

    uint32_t compressed_size = comp_compress (vb, ctx, &vb->z_data, (SectionHeader*)&header, 
                                              callback ? NULL : ctx->local.data, callback, ctx->tag_name);

    ctx->local_in_z_len = vb->z_data.len32 - ctx->local_in_z;

    return compressed_size;
}

// compress section - two options for input data - 
// 1. contiguous data in section_data 
// 2. line by line data - by providing a callback + total_len
void zfile_compress_section_data_ex (VBlockP vb, 
                                     ContextP ctx, // NULL if not context data
                                     SectionType section_type, 
                                     BufferP section_data,          // option 1 - compress contiguous data
                                     LocalGetLineCB callback, uint32_t total_len, // option 2 - compress data one line at a time
                                     Codec codec, SectionFlags flags, 
                                     rom name) 
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

    if (flag.show_time) codec_show_time (vb, name ? name : st_name (section_type), NULL, codec);

    comp_compress (vb, ctx,
                   // note: when called from codec_assign_best_codec we use z_data_test. this is because codec_assign_best_codec can be
                   // called from within complex codecs for their subcodecs, and if we had used z_data, comp_compress could realloc it as it
                   // is being populated by complex codec
                   in_assign_codec ? &vb->z_data_test : &vb->z_data, 
                   &header, 
                   section_data ? section_data->data : NULL, 
                   callback, st_name (section_type));
}

typedef struct { uint64_t start, len; } RemovedSection;

static DESCENDING_SORTER (sort_removed_sections, RemovedSection, start)

// remove ctx and all other ctxs consolidated to it from z_data. akin of unscrambling an egg.
void zfile_remove_ctx_group_from_z_data (VBlockP vb, Did remove_did_i)
{
    unsigned num_rms=0;
    RemovedSection rm[vb->num_contexts * 2];

    // remove all contexts in the group
    CTX(remove_did_i)->st_did_i = remove_did_i; // so the loop catches it too
    for_ctx
        if (ctx->st_did_i == remove_did_i) {
            if (ctx->b250_in_z_len) 
                rm[num_rms++] = (RemovedSection){.start = ctx->b250_in_z, .len = ctx->b250_in_z_len };

            if (ctx->local_in_z_len) 
                rm[num_rms++] = (RemovedSection){.start = ctx->local_in_z, .len = ctx->local_in_z_len};

            vb->recon_size -= ctx->txt_len; // it won't be reconstructed after all

            ctx_substract_txt_len (vb, ctx); // substract txt_len added to zctx during merge

            ctx_free_context (ctx, ctx->did_i);
        }

    // update VB Header (always first in z_data) with reduced recon_size (re-encrypting it if encrypting)
    uint64_t save = vb->z_data.len;
    vb->z_data.len = 0;
    zfile_compress_vb_header (vb);
    vb->z_data.len = save;

    // sort indices to the to-be-removed sections in reverse order
    qsort (rm, num_rms, sizeof(RemovedSection), sort_removed_sections);

    bool is_encrypted = crypt_have_password(); // we can't (easily) test magic if header is encrypted

    for (unsigned i=0; i < num_rms; i++) {
        ASSERT (is_encrypted || ((SectionHeader*)B8 (vb->z_data, rm[i].start))->magic == BGEN32(GENOZIP_MAGIC),
                "Data to be cut out start=%"PRIu64" len=%"PRIu64" is not on section boundary", rm[i].start, rm[i].len);
        
        buf_remove (vb->z_data, char, rm[i].start, rm[i].len);
        sections_remove_from_list (vb, rm[i].start, rm[i].len);

        ASSERT (is_encrypted || rm[i].start == vb->z_data.len || ((SectionHeader*)B8 (vb->z_data, rm[i].start))->magic == BGEN32(GENOZIP_MAGIC),
                "Data cut out is not exactly one section start=%"PRIu64" len=%"PRIu64, rm[i].start, rm[i].len);
    }
}

// reads exactly the length required, error otherwise. 
// return a pointer to the data read
static void *zfile_read_from_disk (File *file, VBlockP vb, BufferP buf, uint32_t len, SectionType st)
{
    START_TIMER;

    ASSERT (len, "reading %s: len is 0", st_name (st));
    ASSERT (buf_has_space (buf, len), "reading %s: buf is out of space: len=%u but remaining space in buffer=%u (tip: run with --show-headers to see where it fails)",
            st_name (st), len, (uint32_t)(buf->size - buf->len));

    char *start = BAFTc (*buf);
    uint32_t bytes = fread (start, 1, len, (FILE *)file->file);
    ASSERT (bytes == len, "reading %s: read only %u bytes out of len=%u", st_name (st), bytes, len);

    buf->len += bytes;
    file->disk_so_far += bytes; // consumed by dispatcher_increment_progress

    COPY_TIMER (read);

    return start;
}


// read section header - called from the main thread. 
// returns offset of header within data, or SECTION_SKIPPED if section is skipped
int32_t zfile_read_section_do (File *file,
                               VBlockP vb, 
                               uint32_t original_vb_i, // the vblock_i used for compressing. this is part of the encryption key. dictionaries are compressed by the compute thread/vb, but uncompressed by the main thread (vb=0)
                               BufferP data, rom buf_name, // buffer to append 
                               SectionType expected_sec_type,
                               Section sec, // NULL for no seeking
                               uint32_t header_size,   
                               FUNCLINE)
{
    ASSERT (!sec || expected_sec_type == sec->st, "called from %s:%u: expected_sec_type=%s but encountered sec->st=%s. vb_i=%u",
            func, code_line, st_name (expected_sec_type), st_name(sec->st), vb->vblock_i);

    // skip if this section is not needed according to flags
    if (sec && file == z_file && 
        piz_is_skip_section (sec->st, (vb ? vb->comp_i : COMP_NONE), (IS_DICTED_SEC (sec->st) ? sec->dict_id : DICT_ID_NONE), 
                             (vb && vb->preprocessing) ? SKIP_PURPOSE_PREPROC : SKIP_PURPOSE_RECON)) 
        return SECTION_SKIPPED; 

    uint32_t unencrypted_header_size = header_size;

    // note: for an encrypted file, while reading the reference, we don't yet know until getting the header whether it
    // will be an SEC_REF_IS_SET (encrypted) or SEC_REFERENCE (not encrypted if originating from external, encryptd if de-novo)
    bool is_encrypted =  (z_file->data_type != DT_REF) && 
                         (expected_sec_type != SEC_GENOZIP_HEADER) &&
                         crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    uint32_t header_offset = data->len;
    buf_alloc (vb, data, 0, header_offset + header_size, uint8_t, 2, buf_name);
    data->param = 1;
    
    // move the cursor to the section. file_seek is smart not to cause any overhead if no moving is needed
    if (sec) file_seek (file, sec->offset, SEEK_SET, false);

    SectionHeader *header = zfile_read_from_disk (file, vb, data, header_size, expected_sec_type); // note: header in file can be shorter than header_size if its an earlier version
    uint32_t bytes_read = header_size;

    ASSERT (header, "called from %s:%u: Failed to read data from file %s while expecting section type %s: %s", 
            func, code_line, z_name, st_name(expected_sec_type), strerror (errno));
    
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
        sections_show_header (header, NULL, sec ? sec->offset : 0, 'R');
        if (is_genocat && (expected_sec_type == SEC_B250      || expected_sec_type == SEC_LOCAL     || 
                                        expected_sec_type == SEC_DICT      || expected_sec_type == SEC_COUNTS    || 
                                        expected_sec_type == SEC_REFERENCE || expected_sec_type == SEC_REF_IS_SET))
             return header_offset; // in genocat --show-header - we only show headers, nothing else
    }

    ASSERT (is_magical || flag.verify_codec, "called from %s:%u: corrupt data (magic is wrong) when attempting to read section=%s dict_id=%s of vblock_i=%u comp=%s in file %s", 
            func, code_line, st_name (expected_sec_type), sec ? dis_dict_id (sec->dict_id).s : "(no sec)", vb->vblock_i, comp_name(vb->comp_i), z_name);

    uint32_t compressed_offset   = BGEN32 (header->compressed_offset);
    ASSERT (compressed_offset, "called from %s:%u: header.compressed_offset is 0 when reading section_type=%s", func, code_line, st_name(expected_sec_type));

    uint32_t data_compressed_len = BGEN32 (header->data_compressed_len);
    uint32_t data_encrypted_len  = BGEN32 (header->data_encrypted_len);

    uint32_t data_len = MAX_(data_compressed_len, data_encrypted_len);

    // in case where we already read part of the body (eg if is_encrypted was initially set and then unset) (remaining_data_len might be negative)
    int32_t remaining_data_len = (int32_t)data_len - (int32_t)(bytes_read - header_size); 
    
    // check that we received the section type we expect, 
    ASSERT (expected_sec_type == header->section_type || 
            (expected_sec_type == SEC_GENOZIP_HEADER && (SectionType)header->sub_codec == SEC_GENOZIP_HEADER), // in v2-5, the section_type field was located where sub_codec is now
            "called from %s:%u: Unexpected section type when reading %s: expecting %s, found %s sec(expecting)=(offset=%s, dict_id=%s)",
            func, code_line, z_name, st_name(expected_sec_type), st_name(header->section_type), 
            sec ? str_int_commas (sec->offset).s : "N/A", sec ? dis_dict_id (sec->dict_id).s : "N/A");

    ASSERT (BGEN32 (header->vblock_i) == original_vb_i, 
            "Requested to read %s with vb_i=%u, but actual section has vb_i=%u",
            st_name(expected_sec_type), original_vb_i, BGEN32 (header->vblock_i));

    ASSERT (compressed_offset == header_size || expected_sec_type == SEC_GENOZIP_HEADER, // we allow SEC_GENOZIP_HEADER of other sizes, for older versions
            "called from %s:%u: invalid header when reading %s - expecting compressed_offset to be %u but found %u. section_type=%s", 
            func, code_line, z_name, header_size, compressed_offset, st_name(header->section_type));

    // allocate more memory for the rest of the header + data (note: after this realloc, header pointer is no longer valid)
    buf_alloc (vb, data, 0, header_offset + compressed_offset + data_len, uint8_t, 2, "zfile_read_section");
    data->param = 2;
    header = (SectionHeader *)Bc(*data, header_offset); // update after realloc

    // read section data - but only if header size is as expected
    if (remaining_data_len > 0 && compressed_offset == header_size)
        zfile_read_from_disk (file, vb, data, remaining_data_len, expected_sec_type);

    return header_offset;
}

// Read one section header - returns the header in vb->scratch - caller needs to free vb->scratch
SectionHeaderUnion zfile_read_section_header (VBlockP vb, uint64_t offset, 
                                              uint32_t original_vb_i, // the vblock_i used for compressing. this is part of the encryption key. dictionaries are compressed by the compute thread/vb, but uncompressed by the main thread (vb=0)
                                              SectionType expected_sec_type)
{
    uint32_t header_size = st_header_size (expected_sec_type);
    uint32_t unencrypted_header_size = header_size;

    // get the uncompressed size from one of the headers - they are all the same size, and the reference file is never encrypted
    file_seek (z_file, offset, SEEK_SET, false);

    bool is_encrypted = (z_file->data_type != DT_REF)             && 
                        (expected_sec_type != SEC_GENOZIP_HEADER) &&
                        crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    SectionHeaderUnion header;
    uint32_t bytes = fread (&header, 1, header_size, (FILE *)z_file->file);
    
    ASSERT (bytes == header_size, "Failed to read header of section type %s from file %s: %s", 
            st_name(expected_sec_type), z_name, strerror (errno));

    bool is_magical = BGEN32 (header.common.magic) == GENOZIP_MAGIC;

    // SEC_REFERENCE is never encrypted when originating from a reference file, it is encrypted (if the file is encrypted) if it originates from REF_INTERNAL 
    if (is_encrypted && header.common.section_type == SEC_REFERENCE && !header.common.data_encrypted_len) {
        is_encrypted = false;
        header_size  = unencrypted_header_size;
    }

    // decrypt header 
    if (is_encrypted) {
        ASSERT (BGEN32 (header.common.magic) != GENOZIP_MAGIC, 
                "password provided, but file %s is not encrypted (sec_type=%s)", z_name, st_name (header.common.section_type));

        crypt_do (vb, (uint8_t*)&header, header_size, original_vb_i, expected_sec_type, true); 
    
        is_magical = BGEN32 (header.common.magic) == GENOZIP_MAGIC; // update after decryption
    }

    ASSERT (is_magical, "corrupt data (magic is wrong) when attempting to read header of section %s in file %s", 
            st_name (expected_sec_type), z_name);

    ASSERT ((BGEN32 (header.common.vblock_i) == original_vb_i && header.common.section_type == expected_sec_type) ||
            (!VER(14) && expected_sec_type == SEC_REF_HASH), // in V<=13, REF_HASH didn't have a vb_i in the section list
            "Requested to read %s with vb_i=%u, but actual section is %s with vb_i=%u",
            st_name(expected_sec_type), original_vb_i, st_name(header.common.section_type), BGEN32 (header.common.vblock_i));

    return header;
}

// check if reference filename exists in the absolute or relative path from the chain header, and if not, 
// check the relative path from the chain file
static rom zfile_read_genozip_header_get_ref_filename (rom header_fn)
{
    // if header_filename exists, use it
    if (file_exists (header_fn)) return header_fn;

    // case absolute path and it doesn't exist 
    if (header_fn[0] == '/' || header_fn[0] == '\\') return NULL;

    rom slash = strrchr (z_name, '/');
    if (!slash && flag.is_windows) slash = strrchr (z_name, '\\');
    if (!slash) return NULL; // chain file is in the current dir

    unsigned dirname_len = slash - z_name + 1; // including slash
    char *fn = MALLOC (strlen (header_fn) + dirname_len + 1); // we're going to leak this memory (+1 for \0)
    sprintf (fn, "%.*s%s", dirname_len, z_name, header_fn);

    return file_exists (fn) ? fn :NULL;
}

// reference data when NOT reading a reference file
static void zfile_read_genozip_header_handle_ref_info (const SectionHeaderGenozipHeader *header)
{
    ASSERT0 (!flag.reading_reference, "we should not be here");

    if (md5_is_zero (header->ref_file_md5)) return; // no reference info in header - we're done

    if (flag.show_reference) {
        iprintf ("%s was compressed using the reference file:\nName: %s\nMD5: %s\n",
                    z_name, header->ref_filename, digest_display (header->ref_file_md5).s);
        if (is_genocat) exit_ok(); // in genocat --show-reference, we only show the reference, not the data
    }

    if (exe_type != EXE_GENOLS) { // note: we don't need the reference for genols

        rom gref_fn = ref_get_filename (gref);

        // case --chain which requires the prim_ref, but an external reference isn't provided (these fields didn't exist before v12)
        if (flag.reading_chain && !flag.explicit_ref && VER(12)) {

            rom prim_ref_filename = zfile_read_genozip_header_get_ref_filename (header->chain.prim_filename);
            rom luft_ref_filename = zfile_read_genozip_header_get_ref_filename (header->ref_filename);

            if (file_exists (prim_ref_filename) && file_exists (luft_ref_filename)) {

                if (!flag.show_chain) 
                    WARN_ONCE ("Note: using the reference files PRIMARY=%s LUFT=%s. You can override this with --reference, see: " WEBSITE_DVCF,
                            prim_ref_filename, luft_ref_filename);
                    
                ref_set_reference (gref,     luft_ref_filename, REF_LIFTOVER, false);
                ref_set_reference (prim_ref, prim_ref_filename, REF_LIFTOVER, false);
            }
            else 
                ASSINP (flag.genocat_no_ref_file, "Please use two --reference arguments to specify the paths to the PRIMAY and LUFT coordinates reference file. Original paths were: PRIMARY=%.*s LUFT=%.*s",
                        REF_FILENAME_LEN, header->chain.prim_filename, REF_FILENAME_LEN, header->ref_filename);
        }

        // case: this file requires an external reference, but command line doesn't include --reference - attempt to use the
        // reference specified in the header. 
        // Note: this code will be executed when zfile_read_genozip_header is called from main_genounzip.
        else if (!flag.explicit_ref && // reference NOT was specified on command line
            !(gref_fn && !strcmp (gref_fn, header->ref_filename))) { // ref_filename already set from a previous file with the same reference

            rom ref_filename = zfile_read_genozip_header_get_ref_filename (header->ref_filename);

            if (!flag.genocat_no_ref_file && ref_filename) {
                WARN ("Note: using the reference file %s. You can override this with --reference", ref_filename);
                ref_set_reference (gref, ref_filename, REF_EXTERNAL, false);
            }
            else 
                ASSINP (flag.genocat_no_ref_file, "Please use --reference to specify the path to the %sreference file. Original path was: %.*s",
                        (Z_DT(DT_CHAIN) ? "LUFT (target) coordinates " : ""), REF_FILENAME_LEN, header->ref_filename);
        }

        // test for matching MD5 between specified external reference and reference in the header
        // note: for --chain, we haven't read the reference yet, so we will check this separately in chain_load
        if (flag.explicit_ref && !flag.reading_chain) 
            digest_verify_ref_is_equal (gref, header->ref_filename, header->ref_file_md5);
    }
}

// gets offset to the beginning of the GENOZIP_HEADER section, and sets z_file->genozip_version
static uint64_t zfile_read_genozip_header_get_offset (void)
{
    // read the footer from the end of the file
    if (file_get_size (z_file->name) < sizeof(SectionFooterGenozipHeader) ||
        !z_file->file ||
        !file_seek (z_file, -sizeof(SectionFooterGenozipHeader), SEEK_END, 2))
        return 0; // failed

    SectionFooterGenozipHeader footer;
    int ret = fread (&footer, sizeof (footer), 1, (FILE *)z_file->file);
    ASSERTW (ret == 1, "Skipping empty file %s", z_name);
    if (!ret) return 0; // failed
    
    // case: there is no genozip header. this can happen if the file was truncated (eg because compression did not complete)
    // note: this can also happen if the file is genozip v1, but I don't think there are any real v1 files in the wild
    // so I will keep the error message simple and not mention it
    RETURNW (BGEN32 (footer.magic) == GENOZIP_MAGIC, 0, "Error: failed to read file %s - the file appears to be incomplete (it is missing the genozip footer).", z_name);
    uint64_t offset = BGEN64 (footer.genozip_header_offset);
    
    // read genozip_version directly, needed to determine the section header size
    RETURNW (file_seek (z_file, offset, SEEK_SET, true), 0, "Error in %s: corrupt offset=%"PRIu64" in genozip footer", z_name, offset);

    char top[sizeof (SectionHeader) + 1];
    RETURNW (fread (top, sizeof top, 1, z_file->file) == 1, 0, "Error in %s: failed to read genozip header", z_name);

    z_file->genozip_version = top[sizeof (SectionHeader)]; 

    ASSRET (z_file->genozip_version <= GENOZIP_FILE_FORMAT_VERSION, 0,
             "Error: %s cannot be opened because it was compressed with a newer version of genozip (version %u) while the version you're running is older (version %s).\n"
             "You must upgrade Genozip to open this file. See: %s\n",
             z_name, z_file->genozip_version, GENOZIP_CODE_VERSION, WEBSITE_INSTALLING);

    // in version 6, we canceled backward compatability with v1-v5
    ASSRET (VER(6), 0, "Skipping %s: it was compressed with an older version of genozip - version %u.\nIt may be uncompressed with genozip versions %u to 5",
                z_name, z_file->genozip_version, z_file->genozip_version);

    // in version 7, we canceled backward compatability with v6
    ASSRET (VER(7), 0, "Skipping %s: it was compressed with version 6 of genozip. It may be uncompressed with genozip version 6",
                z_name);

    // in version 8, we canceled backward compatability with v7
    ASSRET (VER(8), 0, "Skipping %s: it was compressed with version 7 of genozip. It may be uncompressed with genozip version 7",
                z_name);

    return offset;
}

// returns false if file should be skipped
bool zfile_read_genozip_header (SectionHeaderGenozipHeader *out_header) // optional outs
{
    SectionEnt sec = { .st = SEC_GENOZIP_HEADER, 
                       .offset = zfile_read_genozip_header_get_offset() };
    
    if (!sec.offset) goto error;

    zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", SEC_GENOZIP_HEADER, &sec);

    SectionHeaderGenozipHeader *header = (SectionHeaderGenozipHeader *)evb->z_data.data;
    if (out_header) *out_header = *header;

    DataType data_type = (DataType)(BGEN16 (header->data_type)); 
    ASSERT ((unsigned)data_type < NUM_DATATYPES, "unrecognized data_type=%d: please upgrade genozip to the latest version", data_type);

    if (Z_DT(DT_NONE) || Z_DT(DT_GENERIC)) {
        z_file->data_type = data_type;
        z_file->type      = file_get_z_ft_by_dt (z_file->data_type);  

        if (data_type != DT_REF) last_z_dt = data_type;
    }
    else
        ASSINP (z_file->data_type == data_type, "%s - file extension indicates this is a %s file, but according to its contents it is a %s", 
                z_name, dt_name (z_file->data_type), dt_name (data_type));

    flag.genocat_no_ref_file |= (Z_DT(DT_CHAIN) && !flag.reading_chain); // initialized in flags_update: we only need the reference when using the chain file with --chain

    //xxx if (txt_file && !flag.explicit_out_dt /*not explicitly specified by genocat user*/ && header->h.flags.genozip_header.txt_is_bin && 
    //     DTPZ (txt_type))
    //     txt_file->data_type = (is_genocat ? DTPZ(txt_type) : DTPZ (bin_type));

    ASSINP (header->encryption_type != ENC_NONE || !crypt_have_password() || Z_DT(DT_REF), 
            "password provided, but file %s is not encrypted", z_name);

    ASSERT (BGEN32 (header->h.compressed_offset) == st_header_size (SEC_GENOZIP_HEADER),
            "invalid genozip header of %s - expecting compressed_offset to be %u in genozip_version=%u but found %u", 
            z_name, st_header_size (SEC_GENOZIP_HEADER), header->genozip_version, BGEN32 (header->h.compressed_offset));

    // get & test password, if file is encrypted
    if (header->encryption_type != ENC_NONE) {

        if (!crypt_have_password()) crypt_prompt_for_password();

        crypt_do (evb, header->password_test, sizeof(header->password_test), 0, SEC_NONE, true); // decrypt password test

        ASSINP (!memcmp (header->password_test, password_test_string, sizeof(header->password_test)),
                "password is wrong for file %s", z_name);
    }

    z_file->num_components = header->num_components;
    if (z_file->num_components < 2) flag.unbind = 0; // override user's prefix if file has only 1 component (bug 326)

    int dts = z_file->z_flags.dt_specific; // save in case its set already (eg dts_paired is set in fastq_piz_is_paired)
    z_file->z_flags = header->h.flags.genozip_header;
    z_file->z_flags.dt_specific |= dts; 
    z_file->num_lines = BGEN64 (header->num_lines_bound);
    z_file->txt_data_so_far_bind = BGEN64 (header->recon_size_prim);
    segconf.vb_size = (uint64_t)BGEN16 (header->vb_size) << 20;

    DT_FUNC (z_file, piz_genozip_header)(header); // data-type specific processing of the Genozip Header

    if (!z_file->section_list_buf.param) { // not already initialized in a previous call to this function
        zfile_uncompress_section (evb, header, &z_file->section_list_buf, "z_file->section_list_buf", 0, SEC_GENOZIP_HEADER);

        sections_list_file_to_memory_format (header);

        if (flag.show_gheader==1) {
            DO_ONCE sections_show_gheader (header);
            if (is_genocat) exit_ok(); // in genocat, exit after showing the requested data
        }

        z_file->section_list_buf.param = 1;
    }

    // case: we are reading a file expected to be the reference file itself
    if (flag.reading_reference) {
        ASSINP (data_type == DT_REF, "Error: %s is not a reference file. To create a reference file, use 'genozip --make-reference <fasta-file.fa>'",
                ref_get_filename(gref));

        // note: in the reference file itself, header->ref_filename is the original fasta used to create this reference
        ref_set_ref_file_info (flag.reading_reference, header->REF_fasta_md5, header->ref_filename, header->genozip_version); 
    }

    // case: we are reading a file that is not expected to be a reference file
    else {
        // case: we are attempting to decompress a reference file - this is not supported
        ASSERTGOTO (data_type != DT_REF || (flag.genocat_no_reconstruct && is_genocat) || is_genols,
                    "%s is a reference file - it cannot be decompressed. Skipping it.", z_name);

        // handle reference file info
        if (!flag.genocat_no_ref_file)
            zfile_read_genozip_header_handle_ref_info (header);
    }
     
    buf_free (evb->z_data);
    return true;

error:
    buf_free (evb->z_data);
    ASSERT0 (!flag.reading_reference, "failed to read reference file");
    return false;
}

void zfile_compress_genozip_header (void)
{
    START_TIMER;

    SectionHeaderGenozipHeader header = {};

    // start with just the fields needed by sections_add_to_list
    header.h.section_type = SEC_GENOZIP_HEADER;

    header.h.flags.genozip_header  = (struct FlagsGenozipHeader) {
        .txt_is_bin   = DTPT (is_binary),
        .dt_specific  = (DT_FUNC (z_file, zip_dts_flag)()),
        .aligner      = (flag.aligner_available > 0),
        .bgzf         = (txt_file->codec == CODEC_BGZF || txt_file->codec == CODEC_GZ), // note: if txt file is compressed with GZ, we will reconstruct it with BGZF
        .adler        = !flag.md5,
        .has_gencomp  = z_file->z_flags.has_gencomp,
        .has_taxid    = kraken_is_loaded,
    };

    // "manually" add the genozip section to the section list - normally it is added in comp_compress()
    // but in this case the genozip section containing the list will already be ready...
    sections_add_to_list (evb, &header.h);
    sections_list_concat (evb); // usually done in output_process_vb, but the section list will be already compressed within the genozip header...

    bool is_encrypted = crypt_have_password();

    uint32_t num_sections = z_file->section_list_buf.len32;

    sections_list_memory_to_file_format (false); // into evb->scratch, just up to CODEC_ASSIGN_SAMPLE_SIZE
    Codec codec = codec_assign_best_codec (evb, NULL, &evb->scratch, SEC_GENOZIP_HEADER);
    buf_free (evb->scratch);

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.compressed_offset     = BGEN32 (sizeof (SectionHeaderGenozipHeader));
    header.h.data_uncompressed_len = BGEN32 (z_file->section_list_buf.len * sizeof (SectionEntFileFormat));
    header.h.codec                 = codec == CODEC_UNKNOWN ? CODEC_NONE : codec;
    header.genozip_version         = GENOZIP_FILE_FORMAT_VERSION;
    header.data_type               = BGEN16 ((uint16_t)dt_get_txt_dt (z_file->data_type));
    header.encryption_type         = is_encrypted ? ENC_AES256 : ENC_NONE;
    header.recon_size_prim         = BGEN64 (z_file->txt_data_so_far_bind);
    header.num_lines_bound         = BGEN64 (z_file->num_lines);
    header.num_sections            = BGEN32 (num_sections); 
    header.num_components          = z_file->num_txts_so_far;
    header.vb_size                 = BGEN16 (segconf.vb_size >> 20);

    // when decompressing will require an external reference, we set header.ref_filename to the name of the genozip reference file
    if (IS_REF_EXTERNAL || IS_REF_MAKE_CHAIN) {   
        strncpy (header.ref_filename, ref_get_filename (gref), REF_FILENAME_LEN-1);
        header.ref_file_md5 = ref_get_file_md5 (gref);
    }

    // in --make-ref, we set header.ref_filename to the original fasta file, to be used later in ref_get_cram_ref
    // (unless the fasta is piped from stdin, or its name is too long)
    else if (flag.make_reference && strcmp (txt_name, FILENAME_STDIN) && strlen (txt_name) <= REF_FILENAME_LEN-1) {
#ifndef WIN32
        char *ref_filename = realpath (txt_name, NULL); // allocates memory
        ASSERT (ref_filename, "realpath() failed: %s", strerror (errno));

        strncpy (header.ref_filename, ref_filename, REF_FILENAME_LEN-1);
        FREE (ref_filename);
#else
        ASSERT0 (_fullpath (header.ref_filename, txt_name, REF_FILENAME_LEN), "_fullpath() failed");
#endif
    }

    uint32_t license_num_bgen = BGEN32 (license_get_number());
    header.license_hash = md5_do (&license_num_bgen, sizeof (int32_t));
    
    zfile_get_metadata (header.created);

    if (is_encrypted) {
        memcpy (header.password_test, password_test_string, sizeof(header.password_test));
        crypt_do (evb, header.password_test, sizeof(header.password_test), 0, SEC_NONE, true);
    }

    DT_FUNC (z_file, zip_genozip_header)(&header);
    
    BufferP z_data = &evb->z_data;

    uint64_t genozip_header_offset = z_file->disk_so_far + z_data->len; // capture before comp_compress that increases len

    if (flag.show_gheader) sections_show_gheader (&header); 

    // prepare section list for disk - Big Endian and length in bytes
    sections_list_memory_to_file_format (true);

    // compress section into z_data - to be eventually written to disk by the main thread
    comp_compress (evb, NULL, z_data, &header.h, z_file->section_list_buf.data, NO_CALLBACK, "SEC_GENOZIP_HEADER");

    // add a footer to this section - this footer appears AFTER the genozip header data, 
    // facilitating reading the genozip header in reverse from the end of the file
    SectionFooterGenozipHeader footer = { .magic                 = BGEN32 (GENOZIP_MAGIC),
                                          .genozip_header_offset = BGEN64 (genozip_header_offset) };
    buf_add_more (evb, z_data, (char*)&footer, sizeof(SectionFooterGenozipHeader), "z_data");

    COPY_TIMER_VB (evb, zfile_compress_genozip_header);

    zfile_output_processed_vb (evb); // write footer
}

// ZIP
void zfile_write_txt_header (BufferP txt_header, 
                             uint64_t unmodified_txt_header_len, // length of header before modifications, eg due to --chain or compressing a Luft file
                             Digest header_md5, bool is_first_txt, CompIType comp_i)
{
    Codec codec = codec_assign_best_codec (evb, NULL, txt_header, SEC_TXT_HEADER);

    SectionHeaderTxtHeader header = {
        .h.magic                   = BGEN32 (GENOZIP_MAGIC),
        .h.section_type            = SEC_TXT_HEADER,
        .h.data_uncompressed_len   = BGEN32 (txt_header->len),
        .h.compressed_offset       = BGEN32 (sizeof (SectionHeaderTxtHeader)),
        .h.codec                   = (codec == CODEC_UNKNOWN) ? CODEC_NONE : codec,
        .codec                     = txt_file->codec, 
        .digest_header             = flag.data_modified ? DIGEST_NONE : header_md5,
        .txt_header_size           = BGEN64 (unmodified_txt_header_len),
    };

    // data type specific fields
    if (DTPZ(zip_set_txt_header_specific)) DTPZ(zip_set_txt_header_specific)(&header);

    // In BGZF, we store the 3 least significant bytes of the file size, so check if the reconstructed BGZF file is likely the same
    if (txt_file->codec == CODEC_BGZF) 
        bgzf_sign (txt_file->disk_size, header.codec_info);
        
    file_basename (txt_file->name, false, FILENAME_STDIN, header.txt_filename, TXT_FILENAME_LEN);
    file_remove_codec_ext (header.txt_filename, txt_file->type); // eg "xx.fastq.gz -> xx.fastq"
    
    static Buffer txt_header_buf = EMPTY_BUFFER;

    buf_alloc (evb, &txt_header_buf, 0, sizeof (SectionHeaderTxtHeader) + txt_header->len / 3, // generous guess of compressed size
               char, 1, "txt_header_buf"); 

    CompIType save_evb_comp_i = evb->comp_i;
    evb->comp_i = comp_i; // goes into SectionEntFileFormat.comp_i
    comp_compress (evb, NULL, &txt_header_buf, (SectionHeader*)&header, 
                   txt_header->len ? txt_header->data : NULL, // actual header may be missing (eg in SAM it is permitted to not have a header)
                   NO_CALLBACK, "SEC_TXT_HEADER");
    evb->comp_i = save_evb_comp_i; // restore

    file_write (z_file, STRb(txt_header_buf));

    z_file->disk_so_far += txt_header_buf.len;   // length of GENOZIP data writen to disk
    z_file->disk_so_far_comp[comp_i] += txt_header_buf.len;

    // VCF note: we don't account for DVCF rejects files - the added header lines are duplicates of the main header
    if (!z_is_dvcf || comp_i == VCF_COMP_MAIN) {        
        z_file->txt_data_so_far_single   += txt_header->len; // length of txt header as it would be reconstructed (possibly afer modifications)
        z_file->txt_data_so_far_bind     += txt_header->len;
        z_file->txt_data_so_far_single_0 += unmodified_txt_header_len; // length of the original txt header as read from the file
        z_file->txt_data_so_far_bind_0   += unmodified_txt_header_len;
    }

    z_file->txt_data_so_far_bind_comp[comp_i] += txt_header->len;
    z_file->txt_data_so_far_bind_0_comp[comp_i] += unmodified_txt_header_len;

    buf_free (txt_header_buf); 

    memcpy (&z_file->txt_header_single, &header, sizeof (header));
}

// Update SEC_TXT_HEADER. If we're compressing a plain file, we will know
// the bytes upfront, but if we're binding or compressing a eg .GZ, we will need to update it
// when we're done. num_lines can only be known after we're done with this txt component.
// if we cannot update the header - that's fine, these fields are only used for the progress indicator on --list
bool zfile_update_txt_header_section_header (uint64_t offset_in_z_file, uint32_t max_lines_per_vb)
{
    // rewind to beginning of current (latest) vcf header - nothing to do if we can't
    if (!file_seek (z_file, offset_in_z_file, SEEK_SET, true)) return false;

    uint32_t len = crypt_padded_len (sizeof (SectionHeaderTxtHeader));

    // sanity check - we skip empty files, so data is expected
    ASSERT (txt_file->txt_data_so_far_single > 0, "Expecting txt_file->txt_data_so_far_single=%"PRId64" > 0", txt_file->txt_data_so_far_single);
    
    // update the header of the single (current) vcf. 
    SectionHeaderTxtHeader *curr_header = &z_file->txt_header_single;
    curr_header->txt_data_size    = BGEN64 (txt_file->txt_data_so_far_single);
    curr_header->txt_num_lines    = BGEN64 (txt_file->num_lines);
    curr_header->max_lines_per_vb = BGEN32 (max_lines_per_vb);

    if (!flag.data_modified && gencomp_comp_eligible_for_digest(NULL))
        curr_header->digest = digest_snapshot (&z_file->digest_ctx, "file");

    if (flag.show_headers)
        sections_show_header ((SectionHeader *)curr_header, NULL, offset_in_z_file, 'W'); 

    // encrypt if needed
    if (crypt_have_password()) 
        crypt_do (evb, (uint8_t *)curr_header, len, 0, curr_header->h.section_type, true);

    file_write (z_file, curr_header, len);
    fflush ((FILE*)z_file->file); // its not clear why, but without this fflush the bytes immediately after the first header get corrupted (at least on Windows with gcc)
    
    file_seek (z_file, 0, SEEK_END, false); // return to the end of the file

    return true; // success
}

// ZIP compute thread - called from zip_compress_one_vb()
void zfile_compress_vb_header (VBlockP vb)
{
    uint32_t sizeof_header = sizeof (SectionHeaderVbHeader);

    SectionHeaderVbHeader vb_header = {
        .h.magic             = BGEN32 (GENOZIP_MAGIC),
        .h.section_type      = SEC_VB_HEADER,
        .h.compressed_offset = BGEN32 (sizeof_header),
        .h.vblock_i          = BGEN32 (vb->vblock_i),
        .h.codec             = CODEC_NONE,
        .h.flags.vb_header   = vb->flags,
        .recon_size_prim     = BGEN32 (vb->recon_size),
        .longest_line_len    = BGEN32 (vb->longest_line_len),
        .digest              = flag.data_modified ? DIGEST_NONE : vb->digest,
    };

    DT_FUNC (vb, zip_set_vb_header_specific)(vb, &vb_header);

    // copy section header into z_data - to be eventually written to disk by the main thread. this section doesn't have data.
    comp_compress (vb, NULL, &vb->z_data, (SectionHeader*)&vb_header, NULL, NO_CALLBACK, "SEC_VB_HEADER");
}

// ZIP only: called by the main thread in the sequential order of VBs: updating of the already compressed
// variant data section (compressed by the compute thread in zfile_compress_vb_header) just before writing it to disk
// note: this updates the z_data in memory (not on disk)
void zfile_update_compressed_vb_header (VBlockP vb)
{
    if (flag.biopsy) return; // we have no z_data in biopsy mode

    SectionHeaderVbHeader *vb_header = (SectionHeaderVbHeader *)vb->z_data.data;
    vb_header->z_data_bytes = BGEN32 (vb->z_data.len32);

    if (flag.show_vblocks) 
        iprintf ("UPDATE_VB_HEADER(id=%d) vb_i=%u comp_i=%u recon_size_prim=%u recon_size_luft=%u genozip_size=%u longest_line_len=%u\n",
                 vb->id, vb->vblock_i, vb->comp_i, 
                 BGEN32 (vb_header->recon_size_prim), BGEN32 (vb_header->dvcf_recon_size_luft), 
                 BGEN32 (vb_header->z_data_bytes), BGEN32 (vb_header->longest_line_len));

    // now we can finally encrypt the header - if needed
    if (crypt_have_password())  
        crypt_do (vb, (uint8_t*)vb_header, BGEN32 (vb_header->h.compressed_offset),
                  BGEN32 (vb_header->h.vblock_i), vb_header->h.section_type, true);
}

// ZIP
void zfile_output_processed_vb (VBlockP vb)
{
    START_TIMER;

    sections_list_concat (vb);
    
    file_write (z_file, STRb(vb->z_data));
    COPY_TIMER (write);

    z_file->disk_so_far += vb->z_data.len;
    z_file->disk_so_far_comp[vb->comp_i] += vb->z_data.len;
    vb->z_data.len = 0;

    ctx_update_stats (vb);

    if (flag.show_headers && buf_is_alloc (&vb->show_headers_buf))
        buf_print (&vb->show_headers_buf, false);
}

DataType zfile_get_file_dt (rom filename)
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

