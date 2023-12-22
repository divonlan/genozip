// ------------------------------------------------------------------
//   zfile.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

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
#include "zip.h"
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
#include "threads.h"
#include "refhash.h"
#include "seg.h"
#include "dispatcher.h"
#include "qname.h"
#include "libdeflate/libdeflate.h"

static void zfile_show_b250_section (SectionHeaderUnionP header_p, ConstBufferP b250_data)
{
    static Mutex show_b250_mutex = {}; // protect so compute thread's outputs don't get mix

    SectionHeaderCtxP header = header_p.ctx;

    if (!flag.show_b250 && dict_id_typeless (header->dict_id).num != flag.dict_id_show_one_b250.num) return;

    mutex_initialize (show_b250_mutex); // possible unlikely race condition on initializing - good enough for debugging purposes
    mutex_lock (show_b250_mutex);

    iprintf ("vb_i=%u %*.*s: ", BGEN32 (header->vblock_i), -DICT_ID_LEN-1, DICT_ID_LEN, dict_id_typeless (header->dict_id).id);

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
static void zfile_dump_section (BufferP uncompressed_data, SectionHeaderP header, unsigned section_len, DictId dict_id)
{
    char filename[100];
    VBIType vb_i = BGEN32 (header->vblock_i);

    // header
    sprintf (filename, "%s.%u.%s.header", st_name (header->section_type), vb_i, dis_dict_id (dict_id).s);
    file_put_data (filename, header, section_len, 0);

    // body
    if (uncompressed_data->len) {
        sprintf (filename, "%s.%u.%s.body", st_name (header->section_type), vb_i, dis_dict_id (dict_id).s);
        buf_dump_to_file (filename, uncompressed_data, 1, false, false, true, false);
    }
}

// uncompressed a block and adds a \0 at its end. Returns the length of the uncompressed block, without the \0.
// when we get here, the header is already unencrypted zfile_one_section
void zfile_uncompress_section (VBlockP vb,
                               SectionHeaderUnionP header_p,
                               BufferP uncompressed_data, 
                               rom uncompressed_data_buf_name, // a name if Buffer, NULL ok if buffer need not be realloced
                               uint32_t expected_vb_i,
                               SectionType expected_section_type) 
{
    START_TIMER;
    ASSERTNOTNULL (header_p.common);
    
    DictId dict_id = DICT_ID_NONE;
    uint8_t codec_param = 0;

    if (expected_section_type == SEC_DICT)
        dict_id = header_p.dict->dict_id;
    else if (expected_section_type == SEC_B250 || expected_section_type == SEC_LOCAL) {
        dict_id = header_p.ctx->dict_id;
        codec_param = header_p.ctx->param; 
    }
    else if (expected_section_type == SEC_COUNTS) 
        dict_id = header_p.counts->dict_id;
    else if (expected_section_type == SEC_SUBDICTS) 
        dict_id = header_p.subdicts->dict_id;

    ContextP ctx = NULL;
    if (IS_DICTED_SEC (expected_section_type)) { 
        ctx = ECTX(dict_id); 
        if (ctx && !ctx->is_loaded && IS_PIZ)  // note: never skip in ZIP (when an R2 VB uncompressed R1 sections)
            return;  // section was skipped 
    }
    else 
        if (piz_is_skip_undicted_section (expected_section_type)) return; // undicted section was skipped 

    SectionHeaderP header  = header_p.common;
    uint32_t data_encrypted_len    = BGEN32 (header->data_encrypted_len);
    uint32_t data_compressed_len   = BGEN32 (header->data_compressed_len);
    uint32_t data_uncompressed_len = BGEN32 (header->data_uncompressed_len);
    uint32_t expected_z_digest     = BGEN32 (header->z_digest);
    VBIType vblock_i               = BGEN32 (header->vblock_i);

    // sanity checks
    ASSERT (header->section_type == expected_section_type, "expecting section type %s but seeing %s", st_name(expected_section_type), st_name(header->section_type));
    
    ASSERT (vblock_i == expected_vb_i || !expected_vb_i, // dictionaries are uncompressed by the main thread with pseduo_vb (vb_i=0) 
            "bad vblock_i: header->vblock_i=%u but expecting it to be %u (section_type=%s dict_id=%s)", 
            vblock_i, expected_vb_i, st_name (expected_section_type), dis_dict_id(dict_id).s);

    if (flag.show_uncompress)
        iprintf ("Uncompress: %s %-9s %-8s comp_len=%-7u uncomp_len=%u\n", VB_NAME, 
                 st_name (expected_section_type), dict_id.num ? dis_dict_id (dict_id).s : "", data_compressed_len, data_uncompressed_len);
        
    uint32_t compressed_offset = st_header_size (header->section_type);
    if (data_encrypted_len) compressed_offset = ROUNDUP16 (compressed_offset);

    uint32_t actual_z_digest = adler32 (1, (uint8_t*)header + compressed_offset, MAX_(data_compressed_len, data_encrypted_len));

    if (VER(15) && expected_z_digest != actual_z_digest) {
        sections_show_header (header_p.common, vb, 0, 'E');
        ABORT ("%s:%s: Section %s data failed digest verification: expected_z_digest=%u != actual_z_digest=%u", 
               z_name, VB_NAME, st_name(header->section_type), expected_z_digest, actual_z_digest);
    }

    // decrypt data (in-place) if needed
    if (data_encrypted_len)        
        crypt_do (vb, (uint8_t*)header + compressed_offset, data_encrypted_len, vblock_i, header->section_type, false);
    
    bool bad_compression = false;

    if (data_uncompressed_len > 0) { // FORMAT, for example, can be missing in a sample-less file

        if (uncompressed_data_buf_name) {
            buf_alloc (vb, uncompressed_data, 0, data_uncompressed_len + sizeof (uint64_t), char, 1.1, uncompressed_data_buf_name); // add a 64b word for safety in case this buffer will be converted to a bits later
            uncompressed_data->len = data_uncompressed_len;
        }

        comp_uncompress (vb, ctx, header->codec, header->sub_codec, codec_param,
                         (char*)header + compressed_offset, data_compressed_len, 
                         uncompressed_data, data_uncompressed_len,
                         dict_id.num ? dis_dict_id(dict_id).s : st_name(expected_section_type));

        //--verify-codec: verify that adler32 of the decompressed data is equal that of the original uncompressed data
        if (flag.verify_codec && uncompressed_data && data_uncompressed_len && 
            BGEN32 (header->magic) != GENOZIP_MAGIC &&
            header->uncomp_adler32 != adler32 (1, uncompressed_data->data, data_uncompressed_len)) {
        
            iprintf ("--verify-codec: BAD ADLER32 section decompressed incorrectly: codec=%s\n", codec_name(header->codec));
            sections_show_header (header, NULL, 0, 'R');
            bad_compression = true;
        }
    }
 
    if (flag.show_b250 && expected_section_type == SEC_B250) 
        zfile_show_b250_section (header_p, uncompressed_data);
    
    if ((flag.dump_section && !strcmp (st_name (expected_section_type), flag.dump_section)) || bad_compression) {
        uint64_t save_len = uncompressed_data->len;
        uncompressed_data->len = data_uncompressed_len; // might be different, eg in the case of ref_hash
        zfile_dump_section (uncompressed_data, header, compressed_offset, dict_id);
        uncompressed_data->len = save_len; // restore
    }

    if (vb) COPY_TIMER (zfile_uncompress_section);
}

// uncompress into a specific offset in a pre-allocated buffer
void zfile_uncompress_section_into_buf (VBlockP vb, SectionHeaderUnionP header_p, uint32_t expected_vb_i, SectionType expected_section_type,
                                        BufferP dst_buf,
                                        char *dst) // pointer into dst_buf.data
{
    ASSERT (dst >= B1STc(*dst_buf) && dst <= BLSTc(*dst_buf), "expecting dst=%p to be within dst_buf=%s", dst, buf_desc(dst_buf).s);

    Buffer copy = *dst_buf;
    copy.data = dst; // somewhat of a hack
    zfile_uncompress_section (vb, header_p, &copy, NULL, expected_vb_i, expected_section_type); // NULL name prevents buf_alloc
}

uint32_t zfile_compress_b250_data (VBlockP vb, ContextP ctx)
{
    struct FlagsCtx flags = ctx->flags; // make a copy
    
    if (VB_DT(FASTQ))
        flags.paired = (flag.pair == PAIR_R1 && fastq_zip_use_pair_identical (ctx->dict_id)) ||        // "paired" flag in R1 means: "In R2, reconstruct R1 data IFF R2 data is absent" (v15)
                       (flag.pair == PAIR_R2 && fastq_zip_use_pair_assisted (ctx->dict_id, SEC_B250)); // "paired" flag in R2 means: "Reconstruction of R2 requires R2 data as well as R1 data"

    SectionHeaderCtx header = (SectionHeaderCtx) { 
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .section_type          = SEC_B250,
        .data_uncompressed_len = BGEN32 (ctx->b250.len32),
        .codec                 = ctx->bcodec == CODEC_UNKNOWN ? CODEC_RANS8 : ctx->bcodec,
        .vblock_i              = BGEN32 (vb->vblock_i),
        .flags.ctx             = flags,
        .dict_id               = ctx->dict_id,
        .b250_size             = ctx->b250_size
    };

    ctx->b250_in_z = vb->z_data.len32;

    uint32_t compressed_size = comp_compress (vb, ctx, &vb->z_data, &header, ctx->b250.data, NO_CALLBACK, ctx->tag_name);

    ctx->b250_in_z_len = vb->z_data.len32 - ctx->b250_in_z;

    ctx_zip_z_data_exist (ctx);

    return compressed_size;
}

// returns compressed size
uint32_t zfile_compress_local_data (VBlockP vb, ContextP ctx, uint32_t sample_size /* 0 means entire local buffer */)
{   
    struct FlagsCtx flags = ctx->flags; // make a copy

    if (VB_DT(FASTQ))
        flags.paired = (flag.pair == PAIR_R1 && fastq_zip_use_pair_identical (ctx->dict_id)) ||         // "paired" flag in R1 means: "Load R1 data in R2, if R2 data is absent" (v15)
                       (flag.pair == PAIR_R2 && fastq_zip_use_pair_assisted (ctx->dict_id, SEC_LOCAL)); // "paired" flag in R2 means: "Reconstruction of R2 requires R2 data as well as R1 data"

    uint32_t uncompressed_len = ctx->local.len32 * lt_width(ctx);
    
    // case: we're just testing a small sample
    if (sample_size && uncompressed_len > sample_size) 
        uncompressed_len = sample_size;

    SectionHeaderCtx header = (SectionHeaderCtx) {
        .magic                 = BGEN32 (GENOZIP_MAGIC),
        .section_type          = SEC_LOCAL,
        .data_uncompressed_len = BGEN32 (uncompressed_len),
        .codec                 = ctx->lcodec == CODEC_UNKNOWN ? CODEC_RANS8 : ctx->lcodec, // if codec has not been decided yet, fall back on RANS8
        .sub_codec             = ctx->lsubcodec_piz ? ctx->lsubcodec_piz : CODEC_UNKNOWN,
        .vblock_i              = BGEN32 (vb->vblock_i),
        .flags.ctx             = flags,
        .dict_id               = ctx->dict_id,
        .ltype                 = ctx->ltype,
        .param                 = ctx->local_param ? ctx->local.prm8[0] : 0,
    };

    LocalGetLineCB *callback = zip_get_local_data_callback (vb->data_type, ctx);

    ctx->local_in_z = vb->z_data.len32;

    uint32_t compressed_size = comp_compress (vb, ctx, &vb->z_data, &header, 
                                              callback ? NULL : ctx->local.data, callback, ctx->tag_name);

    ctx->local_in_z_len = vb->z_data.len32 - ctx->local_in_z;

    ctx_zip_z_data_exist (ctx);

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
    for_ctx_that (ctx->st_did_i == remove_did_i) {
        if (ctx->b250_in_z_len) 
            rm[num_rms++] = (RemovedSection){.start = ctx->b250_in_z, .len = ctx->b250_in_z_len };

        if (ctx->local_in_z_len) 
            rm[num_rms++] = (RemovedSection){.start = ctx->local_in_z, .len = ctx->local_in_z_len};

        vb->recon_size -= ctx->txt_len; // it won't be reconstructed after all

        ctx_substract_txt_len (vb, ctx); // substract txt_len added to zctx during merge

        buflist_free_ctx (vb, ctx);
    }

    // update VB Header (always first in z_data) with reduced recon_size (re-encrypting it if encrypting)
    uint64_t save = vb->z_data.len;
    vb->z_data.len = 0;
    zfile_compress_vb_header (vb);
    vb->z_data.len = save;

    // sort indices to the to-be-removed sections in reverse order
    qsort (rm, num_rms, sizeof(RemovedSection), sort_removed_sections);

    bool is_encrypted = has_password(); // we can't (easily) test magic if header is encrypted

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
static void *zfile_read_from_disk (FileP file, VBlockP vb, BufferP buf, uint32_t len, SectionType st, DictId dict_id)
{
    START_TIMER;

    ASSERT (len, "reading %s%s: len is 0", st_name (st), cond_str(dict_id.num, " dict_id=", dis_dict_id(dict_id).s));
    ASSERT (buf_has_space (buf, len), "reading %s: buf is out of space: len=%u but remaining space in buffer=%u (tip: run with --show-headers to see where it fails)",
            st_name (st), len, (uint32_t)(buf->size - buf->len));

    char *start = BAFTc (*buf);
    uint32_t bytes = fread (start, 1, len, (FILE *)file->file);
    ASSERT (bytes == len, "reading %s%s read only %u bytes out of len=%u: %s", 
            st_name (st), cond_str(dict_id.num, " dict_id=", dis_dict_id(dict_id).s), bytes, len, strerror(errno));

    buf->len += bytes;
    file->disk_so_far += bytes; // consumed by dispatcher_increment_progress

    COPY_TIMER (read);

    return start;
}


// read section header - called from the main thread. 
// returns offset of header within data, or SECTION_SKIPPED if section is skipped
int32_t zfile_read_section_do (FileP file,
                               VBlockP vb, 
                               uint32_t original_vb_i, // the vblock_i used for compressing. this is part of the encryption key. dictionaries are compressed by the compute thread/vb, but uncompressed by the main thread (vb=0)
                               BufferP data, rom buf_name, // buffer to append 
                               SectionType expected_sec_type,
                               Section sec, // NULL for no seeking
                               FUNCLINE)
{
    ASSERTMAINTHREAD;

    ASSERT (!sec || expected_sec_type == sec->st, "called from %s:%u: expected_sec_type=%s but encountered sec->st=%s. vb_i=%u",
            func, code_line, st_name (expected_sec_type), st_name(sec->st), vb->vblock_i);

    // skip if this section is not needed according to flags
    if (sec && file == z_file && 
        piz_is_skip_section ((vb ? vb->data_type : z_file->data_type), sec->st, (vb ? vb->comp_i : COMP_NONE), (IS_DICTED_SEC (sec->st) ? sec->dict_id : DICT_ID_NONE), 
                             sec->flags.flags,
                             (vb && vb->preprocessing) ? SKIP_PURPOSE_PREPROC : SKIP_PURPOSE_RECON)) 
        return SECTION_SKIPPED; 

    uint32_t header_size = st_header_size (expected_sec_type);
    uint32_t unencrypted_header_size = header_size;

    // note: for an encrypted file, while reading the reference, we don't yet know until getting the header whether it
    // will be an SEC_REF_IS_SET (encrypted) or SEC_REFERENCE (not encrypted if originating from external, encryptd if de-novo)
    bool is_encrypted =  !Z_DT(REF) && 
                         expected_sec_type != SEC_GENOZIP_HEADER &&
                         crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    uint32_t header_offset = data->len;
    buf_alloc (vb, data, 0, header_offset + header_size, uint8_t, 2, buf_name);
    data->param = 1;
    
    // move the cursor to the section. file_seek is smart not to cause any overhead if no moving is needed
    if (sec) file_seek (file, sec->offset, SEEK_SET, HARD_FAIL);

    SectionHeaderP header = zfile_read_from_disk (file, vb, data, header_size, expected_sec_type, IS_DICTED_SEC(sec->st) ? sec->dict_id : DICT_ID_NONE); 
    uint32_t bytes_read = header_size;

    ASSERT (header, "called from %s:%u: Failed to read data from file %s while expecting section type %s: %s", 
            func, code_line, z_name, st_name(expected_sec_type), strerror (errno));
    
    bool is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;

    // SEC_REFERENCE is never encrypted when originating from a reference file, it is encrypted (if the file is encrypted) if it originates from REF_INTERNAL 
    if (is_encrypted && HEADER_IS(REFERENCE) && !header->data_encrypted_len) {
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
        if (is_genocat && (IS_DICTED_SEC (expected_sec_type) || expected_sec_type == SEC_REFERENCE || expected_sec_type == SEC_REF_IS_SET))
             return header_offset; // in genocat --show-header - we only show headers, nothing else
    }

    ASSERT (is_magical || flag.verify_codec, "called from %s:%u: corrupt data (magic is wrong) when attempting to read section=%s dict_id=%s of vblock_i=%u comp=%s in file %s", 
            func, code_line, st_name (expected_sec_type), sec ? dis_dict_id (sec->dict_id).s : "(no sec)", vb->vblock_i, comp_name(vb->comp_i), z_name);

    uint32_t data_compressed_len = BGEN32 (header->data_compressed_len);
    uint32_t data_encrypted_len  = BGEN32 (header->data_encrypted_len);

    uint32_t data_len = MAX_(data_compressed_len, data_encrypted_len);

    // in case where we already read part of the body (eg if is_encrypted was initially set and then unset) (remaining_data_len might be negative)
    int32_t remaining_data_len = (int32_t)data_len - (int32_t)(bytes_read - header_size); 
    
    // check that we received the section type we expect, 
    ASSERT (expected_sec_type == header->section_type,
            "called from %s:%u: Unexpected section type when reading %s: expecting %s, found %s sec(expecting)=(offset=%s, dict_id=%s)",
            func, code_line, z_name, st_name(expected_sec_type), st_name(header->section_type), 
            sec ? str_int_commas (sec->offset).s : "N/A", sec ? dis_dict_id (sec->dict_id).s : "N/A");

    ASSERT (BGEN32 (header->vblock_i) == original_vb_i, 
            "Requested to read %s with vb_i=%u, but actual section has vb_i=%u",
            st_name(expected_sec_type), original_vb_i, BGEN32 (header->vblock_i));

    // up to v14, we had compressed_offset instead of z_digest. Since we have it, we might as well use it
    // as an extra verification of the SectionHeader integrity 
    ASSERT (VER(15) || BGEN32 (header->v14_compressed_offset) == header_size,
            "called from %s:%u: invalid header when reading %s - expecting compressed_offset to be %u but found %u. genozip_version=%u section_type=%s", 
            func, code_line, z_name, header_size, BGEN32 (header->v14_compressed_offset), z_file->genozip_version/*set from footer*/, st_name(header->section_type));

    // allocate more memory for the rest of the header + data 
    buf_alloc (vb, data, 0, header_offset + header_size + data_len, uint8_t, 2, "zfile_read_section");
    header = (SectionHeaderP)Bc(*data, header_offset); // update after realloc
    
    data->param = 2;

    // read section data 
    if (remaining_data_len > 0)
        zfile_read_from_disk (file, vb, data, remaining_data_len, expected_sec_type, sections_get_dict_id (header));

    return header_offset;
}

// Read one section header - returns the header in vb->scratch - caller needs to free vb->scratch
SectionHeaderUnion zfile_read_section_header_do (VBlockP vb, Section sec, 
                                                 SectionType expected_sec_type, // optional: if not SEC_NONE, also verifies section is of expected type
                                                 FUNCLINE)
{
    ASSERT (expected_sec_type == SEC_NONE || sec->st == expected_sec_type, 
            "called from %s:%u: expecting sec.st=%s to be %s", func, code_line, st_name (sec->st), st_name (expected_sec_type));

    uint32_t header_size = st_header_size (sec->st);
    uint32_t unencrypted_header_size = header_size;

    file_seek (z_file, sec->offset, SEEK_SET, HARD_FAIL);

    bool is_encrypted = (z_file->data_type != DT_REF)   && 
                        (sec->st != SEC_GENOZIP_HEADER) &&
                        crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    SectionHeaderUnion header;
    uint32_t bytes = fread (&header, 1, header_size, (FILE *)z_file->file);
    
    ASSERT (bytes == header_size, "called from %s:%u: Failed to read header of section type %s from file %s: %s (bytes=%u header_size=%u)", 
            func, code_line, st_name(sec->st), z_name, strerror (errno), bytes, header_size);

    bool is_magical = BGEN32 (header.common.magic) == GENOZIP_MAGIC;

    // SEC_REFERENCE is never encrypted in references files, or if REF_EXT_STORE is used. 
    // It is encrypted (if the file is encrypted) if REF_INTERNAL is used. 
    if (is_encrypted && header.common.section_type == SEC_REFERENCE && !header.common.data_encrypted_len) {
        is_encrypted = false;
        header_size  = unencrypted_header_size;
    }

    // decrypt header 
    if (is_encrypted) {
        ASSERT (BGEN32 (header.common.magic) != GENOZIP_MAGIC, 
                "called from %s:%u: password provided, but file %s is not encrypted (sec_type=%s)", func, code_line, z_name, st_name (header.common.section_type));

        crypt_do (vb, (uint8_t*)&header, header_size, sec->vblock_i, sec->st, true); 
    
        is_magical = BGEN32 (header.common.magic) == GENOZIP_MAGIC; // update after decryption
    }

    ASSERT (is_magical, "called from %s:%u: corrupt data (magic is wrong) when attempting to read header of section %s in file %s", 
            func, code_line, st_name (sec->st), z_name);

    ASSERT (expected_sec_type == SEC_NONE ||
            (BGEN32 (header.common.vblock_i) == sec->vblock_i && header.common.section_type == sec->st) ||
            (!VER(14) && sec->st == SEC_REF_HASH), // in V<=13, REF_HASH didn't have a vb_i in the section list
            "called from %s:%u: Requested to read %s with vb_i=%u, but actual section is %s with vb_i=%u",
            func, code_line, st_name(sec->st), sec->vblock_i, st_name(header.common.section_type), BGEN32 (header.common.vblock_i));

    return header;
}

// up to v14, we had no explicit "has_digest" flag - we calculate it here by searching for proof of digest.
// since a digest might be 0 by chance, a 0 is not a proof of non-digest, however several 0s are strong enough evidence.
static bool zfile_get_has_digest_up_to_v14 (SectionHeaderGenozipHeaderP header)
{
    // proof: a file was compressed with --md5 (zip verifies --md5 conflicts)
    if (!header->flags.genozip_header.adler) return true;

    // proof: a file is up to v13 with digest_bound 
    if (!VER(14) && !digest_is_zero (header->FASTQ_v13_digest_bound)) return true;

    // search for a non-0 digest in the first 3 TXT/VB headers
    Section sec = NULL;
    for (int i=0 ; i < 3 && sections_next_sec2 (&sec, SEC_TXT_HEADER, SEC_VB_HEADER); i++) {
        SectionHeaderUnion header = zfile_read_section_header (evb, sec, SEC_NONE);
        
        // proof: a TXT_HEADER has a digest of the txt_header (0 if file has no header) or 
        // digest of the entire file. 
        if (sec->st == SEC_TXT_HEADER &&              
               (!digest_is_zero (header.txt_header.digest) || !digest_is_zero (header.txt_header.digest_header)))
            return true; 
    
        // proof: a VB has a digest
        if (sec->st == SEC_VB_HEADER && !digest_is_zero (header.vb_header.digest)) return true;
    }

    return false; // no proof of digest
}

bool zfile_advance_to_next_header (uint64_t *offset, uint64_t *gap)
{
    uint64_t start_offset = *offset;
    file_seek (z_file, start_offset, SEEK_SET, HARD_FAIL);

    char data[128 KB + 4];
    while (1) {
        memset (data, 0, sizeof(data));
        
        uint32_t bytes;
        if (!(bytes = fread (data+4, 1, 128 KB, (FILE *)z_file->file)))
            return false; // possibly 4 bytes of the Footer magic remaining

        // note: we accept a magic in the final 4 bytes of data - this could be a Footer. We 
        // move those last 4 bytes to the next iteration
        for (int i=0; i < bytes; i++) 
            if (BGEN32(GET_UINT32 (&data[i])) == GENOZIP_MAGIC) {
                *offset += i - 4;
                *gap = *offset - start_offset; 
                return true;
            }

        *offset += 128 KB;
        memcpy (data, &data[128 KB], 4);
    }
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
static void zfile_read_genozip_header_handle_ref_info (ConstSectionHeaderGenozipHeaderP header)
{
    ASSERT0 (!flag.reading_reference, "we should not be here");

    if (digest_is_zero (header->ref_genome_digest)) return; // no reference info in header - we're done

    if (flag.show_reference) {
        if (flag.force)
            iprintf ("%s", header->ref_filename);
        else
            iprintf ("%s was compressed using the reference file:\nName: %s\nMD5: %s\n",
                        z_name, header->ref_filename, digest_display (header->ref_genome_digest).s);
        if (is_genocat) exit_ok; // in genocat --show-reference, we only show the reference, not the data
    }

    if (!is_genols) { // note: we don't need the reference for genols

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
        else if (!flag.explicit_ref && !getenv ("GENOZIP_REFERENCE") && // reference NOT was specified on command line
                 !Z_DT(REF) && // for reference files, this field is actual fasta_filename
                 !(gref_fn && !strcmp (gref_fn, header->ref_filename))) { // ref_filename already set from a previous file with the same reference

            rom ref_filename = zfile_read_genozip_header_get_ref_filename (header->ref_filename);

            if (!flag.genocat_no_ref_file && ref_filename) {
                WARN ("Note: using the reference file %s. You can override this with --reference or $GENOZIP_REFERENCE", ref_filename);
                ref_set_reference (gref, ref_filename, REF_EXTERNAL, false);
            }
            else 
                ASSINP (flag.genocat_no_ref_file, "Please use --reference to specify the path to the %sreference file. Original path was: %.*s",
                        (Z_DT(CHAIN) ? "LUFT (target) coordinates " : ""), REF_FILENAME_LEN, header->ref_filename);
        }

        // test for matching digest between specified external reference and reference in the header
        // note: for --chain, we haven't read the reference yet, so we will check this separately in chain_load
        if (flag.explicit_ref && !flag.reading_chain) 
            digest_verify_ref_is_equal (gref, header->ref_filename, header->ref_genome_digest);
    }
}

static uint64_t zfile_read_genozip_header_get_actual_offset (void)
{
    uint32_t size = MIN_(z_file->disk_size, 16 MB);
    file_seek (z_file, z_file->disk_size - size, SEEK_SET, HARD_FAIL);

    ASSERTNOTINUSE (evb->scratch);
    buf_alloc_exact_zero (evb, evb->scratch, size + 100, char, "scratch");
    evb->scratch.len -= 100; // extra allocated memory to ease the scan loop

    int ret = fread (evb->scratch.data, size, 1, (FILE *)z_file->file);
    ASSERT (ret == 1, "Failed to read %u bytes from the end of %s", size, z_name);

    for_buf_back (uint8_t, p, evb->scratch)
        if (BGEN32(GET_UINT32(p)) == GENOZIP_MAGIC && ((SectionHeaderP)p)->section_type == SEC_GENOZIP_HEADER) 
            return BNUM (evb->scratch, p) + (z_file->disk_size - size);

    ABORT ("Cannot locate the SEC_GENOZIP_HEADER in the final %u bytes of %s", size, z_name);
}

// gets offset to the beginning of the GENOZIP_HEADER section, and sets z_file->genozip_version
uint64_t zfile_read_genozip_header_get_offset (bool as_is)
{
    // read the footer from the end of the file
    if (z_file->disk_size < sizeof(SectionFooterGenozipHeader) ||
        !z_file->file ||
        !file_seek (z_file, -sizeof(SectionFooterGenozipHeader), SEEK_END, SOFT_FAIL))
        return 0; // failed

    TEMP_FLAG(quiet, false);

    SectionFooterGenozipHeader footer;
    int ret = fread (&footer, sizeof (footer), 1, (FILE *)z_file->file);
    ASSERTW (ret == 1, "Skipping empty file %s", z_name);
    if (!ret) return 0; // failed
    
    // case: there is no genozip header. this can happen if the file was truncated (eg because compression did not complete)
    // note: this can also happen if the file is genozip v1, but I don't think there are any real v1 files in the wild
    // so I will keep the error message simple and not mention it
    RETURNW (BGEN32 (footer.magic) == GENOZIP_MAGIC, 0, "Error in %s: the file appears to be incomplete (it is missing the Footer).", z_name);
    
    uint64_t offset = flag.recover ? zfile_read_genozip_header_get_actual_offset() // get correct offset in case of corruption
                                   : BGEN64 (footer.genozip_header_offset);
        
    if (as_is) return offset;

    // read genozip_version directly, needed to determine the section header size
    RETURNW (file_seek (z_file, offset, SEEK_SET, WARNING_FAIL), 0, "Error in %s: corrupt offset=%"PRIu64" in Footer  (file_size=%"PRIu64")", 
             z_name, offset, z_file->disk_size);

    SectionHeaderGenozipHeader top;
    RETURNW (fread (&top, sizeof (SectionHeaderGenozipHeader), 1, z_file->file) == 1, 0, "Error in %s: failed to read genozip header", z_name);

    RETURNW (BGEN32 (top.magic) == GENOZIP_MAGIC, 0, "Error in %s: offset=%"PRIu64" of the GENOZIP_HEADER section as it appears in the Footer appears to be wrong, or the GENOZIP_HEADER section has bad magic (file_size=%"PRIu64").%s", 
             z_name, offset, z_file->disk_size, flag.debug_or_test ? " Try again with --recover." : "");

    RESTORE_FLAG(quiet);

    z_file->genozip_version   = top.genozip_version;
    z_file->genozip_minor_ver = top.genozip_minor_ver; // 0 before 15.0.28

    z_file->data_type = BGEN16 (top.data_type);

    // check that file version is at most this executable version, except for reference file for which only major version is tested
    ASSINP (z_file->genozip_version < code_version_major() || 
            (z_file->genozip_version == code_version_major() && (z_file->genozip_minor_ver <= code_version_minor() || Z_DT(REF))),
            "Error: %s cannot be opened because it was compressed with genozip version %u.0.%u which is newer than the version running - %s.\n%s",
            z_name, z_file->genozip_version, z_file->genozip_minor_ver, GENOZIP_CODE_VERSION, genozip_update_msg());

    bool metadata_only = is_genocat && (flag.show_stats || flag.show_gheader || flag.show_headers || flag.show_aliases || flag.show_dict);

    #define MSG "Error: %s was compressed with version %u of genozip. It may be uncompressed with genozip versions %u to %u"

    // in version 6, we canceled backward compatability with v1-v5
    ASSINP (VER(6), MSG, z_name, z_file->genozip_version, z_file->genozip_version, 5);

    // in version 7, we canceled backward compatability with v6
    ASSINP (VER(7), MSG, z_name, z_file->genozip_version, 6, 6);

    // in version 8, we canceled backward compatability with v7
    ASSINP (VER(8), MSG, z_name, z_file->genozip_version, 7, 7);

    // in version 15, we canceled backward compatability with v8,9,10 (except reference files which continue to be supported back to v8, as they might be needed to decompress files of later versions)
    ASSINP (metadata_only || VER(11) || Z_DT(REF), MSG, z_name, z_file->genozip_version, z_file->genozip_version, 14);

    #undef MSG
    return offset;
}

// returns false if file should be skipped
bool zfile_read_genozip_header (SectionHeaderGenozipHeaderP out_header, FailType fail_type) // optional outs
{
    if (z_file->section_list_buf.len) return true; // header already read

    SectionEnt sec = { .st     = SEC_GENOZIP_HEADER, 
                       .offset = zfile_read_genozip_header_get_offset (false) };
    
    if (!sec.offset) goto error;

    zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", SEC_GENOZIP_HEADER, &sec);

    SectionHeaderGenozipHeaderP header = (SectionHeaderGenozipHeaderP)evb->z_data.data;
    if (out_header) *out_header = *header;

    DataType data_type = (DataType)(BGEN16 (header->data_type)); 
    ASSERT ((unsigned)data_type < NUM_DATATYPES, "unrecognized data_type=%d. %s", data_type, genozip_update_msg());

    if (Z_DT(NONE) || Z_DT(GENERIC)) {
        z_file->data_type = data_type;
        z_file->type      = file_get_z_ft_by_dt (z_file->data_type);  
    }
    else
        ASSINP (z_file->data_type == data_type, "%s - file extension indicates this is a %s file, but according to its contents it is a %s", 
                z_name, z_dt_name(), dt_name (data_type));

    flag.genocat_no_ref_file |= (Z_DT(CHAIN) && !flag.reading_chain); // initialized in flags_update: we only need the reference when using the chain file with --chain

    ASSINP (header->encryption_type != ENC_NONE || !has_password() || Z_DT(REF), 
            "password provided, but file %s is not encrypted", z_name);

    ASSERT (VER(15) || BGEN32 (header->v14_compressed_offset) == st_header_size (SEC_GENOZIP_HEADER),
            "invalid genozip header of %s - expecting compressed_offset to be %u in genozip_version=%u but found %u", 
            z_name, st_header_size (SEC_GENOZIP_HEADER), header->genozip_version, BGEN32 (header->v14_compressed_offset));

    // get & test password, if file is encrypted
    if (header->encryption_type != ENC_NONE) {

        if (!has_password()) crypt_prompt_for_password();

        crypt_do (evb, header->password_test, sizeof(header->password_test), 0, SEC_NONE, true); // decrypt password test

        ASSINP (!memcmp (header->password_test, PASSWORD_TEST, sizeof(header->password_test)),
                "password is wrong for file %s", z_name);
    }

    z_file->num_txt_files = header->num_txt_files;
    if (z_file->num_txt_files < 2) flag.unbind = 0; // override user's prefix if file has only 1 component (bug 326)

    int dts = z_file->z_flags.dt_specific; // save in case its set already (eg dts_paired is set in sections_is_paired)
    z_file->z_flags = header->flags.genozip_header;

    z_file->z_flags.dt_specific |= dts; 
    z_file->num_lines = BGEN64 (header->num_lines_bound);
    z_file->txt_data_so_far_bind = BGEN64 (header->recon_size_prim);
    segconf.vb_size = (uint64_t)BGEN16 (header->vb_size) MB;

    if (flag.show_data_type) {
        iprintf ("%s\n", z_dt_name());
        exit_ok;
    }

    DT_FUNC (z_file, piz_genozip_header)(header); // data-type specific processing of the Genozip Header

    bool has_section_list = true; 
    if (!z_file->section_list_buf.param) { // not already initialized in a previous call to this function
        
        has_section_list = license_piz_prepare_genozip_header (header, IS_LIST || (IS_SHOW_HEADERS && flag.force));

        if (has_section_list) {
            zfile_uncompress_section (evb, header, &z_file->section_list_buf, "z_file->section_list_buf", 0, SEC_GENOZIP_HEADER);

            sections_list_file_to_memory_format (header);
        }

        if (flag.show_gheader==1) {
            DO_ONCE sections_show_gheader (header);
            if (is_genocat) exit_ok; // in genocat, exit after showing the requested data
        }

        z_file->section_list_buf.param = 1;
    }

    if (!VER(15)) 
        z_file->z_flags.has_digest = zfile_get_has_digest_up_to_v14 (header); // overwrites v14_bgzf that is no longer used for PIZ

    // case: we are reading a file expected to be the reference file itself
    if (flag.reading_reference) {
        ASSINP (data_type == DT_REF, "Error: %s is not a reference file. To create a reference file, use 'genozip --make-reference <fasta-file.fa>'",
                ref_get_filename(gref));

        // note: in the reference file itself, header->ref_filename is the original fasta used to create this reference
        ref_set_ref_file_info (flag.reading_reference, header->genome_digest, header->flags.genozip_header.adler, 
                               header->fasta_filename, header->genozip_version); 

        refhash_set_digest (header->refhash_digest);
        
        buf_free (evb->z_data);
    }

    // case: we are reading a file that is not expected to be a reference file
    else {
        // case: we are attempting to decompress a reference file - this is not supported
        ASSGOTO (data_type != DT_REF || (flag.genocat_no_reconstruct && is_genocat) || is_genols,
                    "%s is a reference file - it cannot be decompressed - skipping it. Did you intend to use --reference?.", z_name);

        // handle reference file info
        flags_update_piz_no_ref_file();
        
        if (!flag.genocat_no_ref_file && data_type != DT_REF)
            zfile_read_genozip_header_handle_ref_info (header);

        buf_free (evb->z_data); // free before ctx_piz_initialize_zctxs that might read aliases - header not valid after freeing

        // create all contexts for B250/LOCAL/DICT data in the z_file (or predefined) - 
        // flags_update_piz_one_z_file and IS_SKIP functions may rely on Context.z_data_exists
        if (has_section_list) 
            ctx_piz_initialize_zctxs();
    }
     
    return true;

error:
    buf_free (evb->z_data);
    ASSERT (fail_type == SOFT_FAIL, "failed to read %s", z_name);
    return false;
}

// Update SEC_TXT_HEADER. If we're compressing a plain file, we will know
// the bytes upfront, but if we're binding or compressing a eg .GZ, we will need to update it
// when we're done. num_lines can only be known after we're done with this txt component.
// if we cannot update the header - that's fine, these fields are only used for the progress indicator on --list
bool zfile_update_txt_header_section_header (uint64_t offset_in_z_file)
{
    // rewind to beginning of current (latest) vcf header - nothing to do if we can't
    if (!file_seek (z_file, offset_in_z_file, SEEK_SET, WARNING_FAIL)) return false;

    uint32_t len = crypt_padded_len (sizeof (SectionHeaderTxtHeader));

    // sanity check - we skip empty files, so data is expected
    ASSERT (txt_file->txt_data_so_far_single > 0, "Expecting txt_file->txt_data_so_far_single=%"PRId64" > 0", txt_file->txt_data_so_far_single);
    
    // update the first txt_header fragment of the single (current) txt file. 
    SectionHeaderTxtHeaderP header = &z_file->txt_header_single;
    header->txt_data_size    = BGEN64 (txt_file->txt_data_so_far_single);
    header->txt_num_lines    = BGEN64 (txt_file->num_lines);
    header->max_lines_per_vb = BGEN32 (txt_file->max_lines_per_vb);

    // qname stuff
    for (QType q=0; q < NUM_QTYPES; q++)
        header->flav_prop[q] = segconf.flav_prop[q];

    if (flag.md5 && !flag.data_modified && gencomp_comp_eligible_for_digest(NULL))
        header->digest = digest_snapshot (&z_file->digest_ctx, "file");

    if (flag.show_headers)
        sections_show_header ((SectionHeaderP)header, NULL, offset_in_z_file, 'W'); 

    // encrypt if needed
    if (has_password()) 
        crypt_do (evb, (uint8_t *)header, len, 1 /*was 0 up to 14.0.8*/, header->section_type, true);

    file_write (z_file, header, len);
    fflush ((FILE*)z_file->file); // its not clear why, but without this fflush the bytes immediately after the first header get corrupted (at least on Windows with gcc)
    
    file_seek (z_file, 0, SEEK_END, HARD_FAIL); // return to the end of the file

    return true; // success
}

// ZIP compute thread - called from zip_compress_one_vb()
void zfile_compress_vb_header (VBlockP vb)
{

    SectionHeaderVbHeader vb_header = {
        .magic             = BGEN32 (GENOZIP_MAGIC),
        .section_type      = SEC_VB_HEADER,
        .vblock_i          = BGEN32 (vb->vblock_i),
        .codec             = CODEC_NONE,
        .flags.vb_header   = vb->flags,
        .recon_size_prim   = BGEN32 (vb->recon_size),
        .longest_line_len  = BGEN32 (vb->longest_line_len),
        .longest_seq_len   = BGEN32 (vb->longest_seq_len), // since v15 (non-0 for SAM, BAM, FASTQ)
        .digest            = flag.data_modified ? DIGEST_NONE : vb->digest,
    };

    DT_FUNC (vb, zip_set_vb_header_specific)(vb, &vb_header);

    // copy section header into z_data - to be eventually written to disk by the main thread. this section doesn't have data.
    comp_compress (vb, NULL, &vb->z_data, &vb_header, NULL, NO_CALLBACK, "SEC_VB_HEADER");
}

// ZIP only: called by the main thread in the sequential order of VBs: updating of the already compressed
// variant data section (compressed by the compute thread in zfile_compress_vb_header) just before writing it to disk
// note: this updates the z_data in memory (not on disk)
void zfile_update_compressed_vb_header (VBlockP vb)
{
    if (flag.biopsy) return; // we have no z_data in biopsy mode

    SectionHeaderVbHeaderP vb_header = (SectionHeaderVbHeaderP)vb->z_data.data;
    vb_header->z_data_bytes = BGEN32 (vb->z_data.len32);

    if (flag_is_show_vblocks (ZIP_TASK_NAME)) 
        iprintf ("UPDATE_VB_HEADER(id=%d) vb_i=%u comp_i=%u recon_size_prim=%u recon_size_luft=%u genozip_size=%u longest_line_len=%u\n",
                 vb->id, vb->vblock_i, vb->comp_i, 
                 BGEN32 (vb_header->recon_size_prim), BGEN32 (vb_header->dvcf_recon_size_luft), 
                 BGEN32 (vb_header->z_data_bytes), BGEN32 (vb_header->longest_line_len));

    // now we can finally encrypt the header - if needed
    if (has_password())  
        crypt_do (vb, (uint8_t*)vb_header, ROUNDUP16(sizeof(SectionHeaderVbHeader)),
                  BGEN32 (vb_header->vblock_i), vb_header->section_type, true);
}

// ZIP - main thread
void zfile_output_processed_vb (VBlockP vb)
{
    START_TIMER;
    ASSERTMAINTHREAD;

    sections_list_concat (vb);
    
    if (!flag.zip_no_z_file) 
        file_write (z_file, STRb(vb->z_data));

    COPY_TIMER (write);

    z_file->disk_so_far += vb->z_data.len;
    if (vb->comp_i != COMP_NONE) z_file->disk_so_far_comp[vb->comp_i] += vb->z_data.len;
    vb->z_data.len = 0;

    // sanity
    if (!flag.zip_no_z_file) {
        uint64_t actual_disk_so_far = file_tell (z_file, false);
        ASSERT (actual_disk_so_far == z_file->disk_so_far, "%s: Expecting actual_disk_so_far=%"PRIu64" == z_file->disk_so_far=%"PRIu64,
                VB_NAME, actual_disk_so_far, z_file->disk_so_far);
    }
    
    ctx_update_stats (vb);

    if (flag.show_headers && buf_is_alloc (&vb->show_headers_buf))
        buf_print (&vb->show_headers_buf, false);
}

DataType zfile_get_file_dt (rom filename)
{
    FileType ft = file_get_type (filename);
    DataType dt = file_get_dt_by_z_ft (ft);
    FileP file = NULL;

    // case: we don't know yet what file type this is - we need to read the genozip header to determine
    if (dt == DT_NONE && filename) {
        if (!(file = file_open_z_read (filename)) || !file->file)
            goto done; // not a genozip file

        // read the footer from the end of the file
        if (!file_seek (file, -sizeof(SectionFooterGenozipHeader), SEEK_END, WARNING_FAIL))
            goto done;

        SectionFooterGenozipHeader footer;
        int ret = fread (&footer, sizeof (footer), 1, (FILE *)file->file);
        ASSERTW (ret == 1, "Skipping empty file %s", z_name);    
        if (!ret) goto done; // empty file / cannot read
        
        // case: this is not a valid genozip v2+ file
        if (BGEN32 (footer.magic) != GENOZIP_MAGIC) goto done;

        // read genozip header
        uint64_t genozip_header_offset = BGEN64 (footer.genozip_header_offset);
        if (!file_seek (file, genozip_header_offset, SEEK_SET, WARNING_FAIL))
            goto done;

        SectionHeaderGenozipHeader header;
        int bytes = fread ((char*)&header, 1, sizeof(SectionHeaderGenozipHeader), (FILE *)file->file);
        if (bytes < sizeof(SectionHeaderGenozipHeader)) goto done;

        ASSERTW (BGEN32 (header.magic) == GENOZIP_MAGIC, "Error reading %s: corrupt data", z_name);
        if (BGEN32 (header.magic) != GENOZIP_MAGIC) goto done;

        dt = (DataType)BGEN16 (header.data_type);
    }

done:
    file_close (&file);
    return dt;
}

