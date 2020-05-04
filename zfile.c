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
#include "samples.h"
#include "squeeze_vcf.h"
#include "header.h"
#include "gtshark_vcf.h"
#include "compressor.h"
#include "piz.h"
#include "license.h"

static const char *password_test_string = "WhenIThinkBackOnAllTheCrapIlearntInHighschool";

void zfile_show_header (const SectionHeader *header, VBlock *vb /* optional if output to buffer */)
{
    DictIdType dict_id = {0};

    if (section_type_is_dictionary (header->section_type)) dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    if (section_type_is_b250       (header->section_type)) dict_id = ((SectionHeaderBase250    *)header)->dict_id;

    char str[1000];

    sprintf (str, "%-23s %*.*s vb_i=%-3u sec_i=%-2u comp_offset=%-6u uncomp_len=%-6u comp_len=%-6u enc_len=%-6u magic=%8.8x\n",
             st_name(header->section_type), -DICT_ID_LEN, DICT_ID_LEN, dict_id.num ? dict_id_printable (dict_id).id : dict_id.id,
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
        fprintf (stderr, "%s", str);
}

static void zfile_show_b250_section (void *section_header_p, Buffer *b250_data)
{
    SectionHeaderBase250 *header = (SectionHeaderBase250 *)section_header_p;

    if (!flag_show_b250 && dict_id_printable (header->dict_id).num != dict_id_show_one_b250.num) return;

    if (flag_show_b250 && header->h.section_type == SEC_CHROM_B250)
        fprintf (stderr, "Base-250 data for VB %u (result of '--show-b250'):\n", BGEN32 (header->h.vblock_i));
        
    fprintf (stderr, "  %*.*s: ", -DICT_ID_LEN-1, DICT_ID_LEN, dict_id_printable (header->dict_id).id);

    const uint8_t *data = (const uint8_t *)b250_data->data;
    for (unsigned i=0; i < BGEN32 (header->num_b250_items); i++) {
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

bool zfile_is_skip_section (void *vb_, SectionType st, DictIdType dict_id)
{
    VBlock *vb = (VBlock *)vb_;

    static const struct {char *name; bool strip;} abouts[NUM_SEC_TYPES] = SECTIONTYPE_ABOUT;

    if (flag_strip && abouts[st].strip) return true;

    // other situations not covered by abouts.strip
    switch (vb ? vb->data_type : z_file->data_type) {
        case DT_VCF:
            if ((flag_drop_genotypes || flag_gt_only) && 
                (st == SEC_VCF_FORMAT_B250 || st == SEC_VCF_FORMAT_DICT || st == SEC_VCF_GT_DATA || st == SEC_VCF_FRMT_SF_DICT))
                return true;

            break;

        case DT_SAM:
            // note: we need OPTIONAL to check for E2 that consumes SEQ, and we need OPTNL_SF for OA, SA, and XA only
            if (flag_strip && 
                (st == SEC_SAM_OPTNL_SF_DICT || st == SEC_SAM_OPTNL_SF_B250) && 
                (dict_id.num != dict_id_OPTION_SA && dict_id.num != dict_id_OPTION_OA && dict_id.num != dict_id_OPTION_XA))
                return true;

            break;

        case DT_FASTQ:
        case DT_FASTA:
            if (flag_header_one &&
                (st == SEC_SEQ_DATA || st == SEC_QUAL_DATA || st == SEC_FASTA_COMMENT_DATA))
                return true;

        default: break; // other sections don't have special logic not already covered in abouts
    }

    return false; // don't skip this section
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
        dict_id = ((SectionHeaderBase250 *)section_header_p)->dict_id;

    if (zfile_is_skip_section (vb, expected_section_type, dict_id)) return; // we skip some sections based on flags

    SectionHeader *section_header = (SectionHeader *)section_header_p;
    
    uint32_t compressed_offset     = BGEN32 (section_header->compressed_offset);
    uint32_t data_encrypted_len    = BGEN32 (section_header->data_encrypted_len);
    uint32_t data_compressed_len   = BGEN32 (section_header->data_compressed_len);
    uint32_t data_uncompressed_len = BGEN32 (section_header->data_uncompressed_len);
    uint32_t vblock_i              = BGEN32 (section_header->vblock_i);
    uint16_t section_i             = BGEN16 (section_header->section_i);
    CompressionAlg comp_alg        = (CompressionAlg)section_header->sec_compression_alg;

    // prior to version 5, the algorithm was hard coded, and this field was called "unused"
    if (z_file->genozip_version < 5) {
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
        if (z_file->genozip_version > 1)
            crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, vblock_i, section_header->section_type, false);
        else
            v1_crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, vblock_i, section_i);
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

    ASSERT (section_type_is_dictionary(ctx->dict_section_type),
            "Error: dict_section_type=%s is not a dictionary section", st_name(ctx->dict_section_type));

    SectionHeaderDictionary header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = ctx->dict_section_type;
    header.h.data_uncompressed_len = BGEN32 (num_chars);
    header.h.compressed_offset     = BGEN32 (sizeof(SectionHeaderDictionary));
    header.h.sec_compression_alg  = COMP_BZ2;
    header.h.vblock_i              = BGEN32 (vb->vblock_i);
    header.h.section_i             = BGEN16 (vb->z_next_header_i++);
    header.num_snips               = BGEN32 (num_words);
    header.dict_id                 = ctx->dict_id;

    if (flag_show_dict)
        fprintf (stderr, "%.*s (vb_i=%u, %s, did=%u, num_snips=%u):\t%.*s\n", 
                DICT_ID_LEN, dict_id_printable (ctx->dict_id).id, vb->vblock_i, st_name(ctx->dict_section_type), 
                ctx->did_i, num_words, num_chars, data);

    if (dict_id_printable (ctx->dict_id).num == dict_id_show_one_dict.num)
        fprintf (stderr, "%.*s\t", num_chars, data);

    z_file->dict_data.name  = "z_file->dict_data"; // comp_compress requires that it is set in advance
    comp_compress (vb, &z_file->dict_data, true, (SectionHeader*)&header, data, NULL);

    COPY_TIMER (vb->profile.zfile_compress_dictionary_data)    
//printf ("End compress dict vb_i=%u did_i=%u\n", vb->vblock_i, ctx->did_i);
}

void zfile_compress_b250_data (VBlock *vb, MtfContext *ctx)
{
    SectionHeaderBase250 header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = ctx->b250_section_type;
    header.h.data_uncompressed_len = BGEN32 (ctx->b250.len);
    header.h.compressed_offset     = BGEN32 (sizeof(SectionHeaderBase250));
    header.h.sec_compression_alg  = COMP_BZ2;
    header.h.vblock_i              = BGEN32 (vb->vblock_i);
    header.h.section_i             = BGEN16 (vb->z_next_header_i++);
    header.dict_id                 = ctx->dict_id;
    header.num_b250_items          = BGEN32 (ctx->mtf_i.len);
            
    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, ctx->b250.data, NULL);
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
static void *zfile_read_from_disk (VBlock *vb, Buffer *buf, unsigned len, bool fail_quietly_if_not_enough_data)
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

    if (sl && zfile_is_skip_section (vb, expected_sec_type, sl->dict_id)) return 0; // skip if this section is not needed according to flags

    if (sb_i != NO_SB_I && ! *ENT(bool, ((VBlockVCFP)vb)->is_sb_included, sb_i)) return 0; // skip section if this sample block is excluded by --samples

    // note: the first section is always read by zfile_read_section(). if it is a v1, or
    // if it is not readable (perhaps v1 encrypted) - we assume its a v1 VCF header
    // in v2+ the first section is always an unencrypted SectionHeaderGenozipHeader
    ASSERT0 (z_file->genozip_version != 1, "Error: zfile_read_section cannot read v1 data");

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
        zfile_read_from_disk (vb, data, sizeof (Md5Hash), false);
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

    mtf_initialize_primary_field_ctxs (z_file->mtf_ctx, z_file->data_type, 0, z_file->dict_id_to_did_i_map, &z_file->num_dict_ids);

    while (sections_get_next_dictionary (&sl_ent)) {

        static SectionType chrom_dict_sec[NUM_DATATYPES] = CHROM_FIELD_DICT_SECTION;

        if (last_vb_i && sl_ent->vblock_i > last_vb_i) break;

        // cases where we can skip reading these dictionaries because we don't be using them
        SectionType st = sl_ent->section_type; 
        if (read_chrom == DICTREAD_CHROM_ONLY   && st != chrom_dict_sec[z_file->data_type]) continue;
        if (read_chrom == DICTREAD_EXCEPT_CHROM && st == chrom_dict_sec[z_file->data_type]) continue;
        if (zfile_is_skip_section (NULL, st, sl_ent->dict_id)) continue;
        
        zfile_read_section (evb, sl_ent->vblock_i, NO_SB_I, &evb->z_data, "z_data", sizeof(SectionHeaderDictionary), st, sl_ent);    

        // update dictionaries in z_file->mtf_ctx with dictionary data 
        mtf_integrate_dictionary_fragment (evb, evb->z_data.data);

        buf_free (&evb->z_data);
    }

    if (flag_show_dict || dict_id_show_one_dict.num) {
        for (uint32_t did_i=0; did_i < z_file->num_dict_ids; did_i++) {
            MtfContext *ctx = &z_file->mtf_ctx[did_i];

#define MAX_PRINTABLE_DICT_LEN 100000

            if (dict_id_printable (ctx->dict_id).num == dict_id_show_one_dict.num) 
                fprintf (stderr, "%.*s\t", (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), ctx->dict.data);
            
            if (flag_show_dict)
                fprintf (stderr, "%.*s (%s, did=%u, num_snips=%u):\t%.*s\n", 
                         DICT_ID_LEN, dict_id_printable (ctx->dict_id).id, st_name(ctx->dict_section_type), 
                         did_i, (uint32_t)ctx->word_list.len, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), ctx->dict.data);
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

        crypt_do (evb, header->password_test, sizeof(header->password_test), 0, SEC_EOF, true); // decrypt password test

        ASSERT (!memcmp (header->password_test, password_test_string, sizeof(header->password_test)),
                "Error: password is wrong for file %s", z_name);
    }

    global_vcf_num_samples    = BGEN32 (header->num_samples); // possibly 0, if genozip header was not rewritten. in this case, piz will get it from the first VCF header, but genols will show 0
    z_file->genozip_version   = header->genozip_version;
    z_file->num_components    = BGEN32 (header->num_components);
    *digest                   = header->md5_hash_concat; 

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
        crypt_do (evb, header.password_test, sizeof(header.password_test), 0, SEC_EOF, true);
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
    evb->z_section_bytes[SEC_GENOZIP_HEADER] += sizeof(SectionFooterGenozipHeader);
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
        return v1_vcf_header_get_vcf_header (uncompressed_data_size, num_samples, num_items_concat, 
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

    if (header.genozip_version >= 5)
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

// ZIP compute thread - called from zip_sam_compress_one_vb()
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

//------------------------------------------------------------------------------
// VCF stuff
//------------------------------------------------------------------------------

// called by ZIP compute thread - called from zip_vcf_compress_one_vb()
void zfile_vcf_compress_vb_header (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    uint32_t my_squeeze_len = squeeze_len (vb->num_haplotypes_per_line);
    uint32_t sizeof_header = sizeof (SectionHeaderVbHeaderVCF);

    SectionHeaderVbHeaderVCF vb_header;
    memset (&vb_header, 0, sizeof(SectionHeaderVbHeaderVCF)); // safety
    
    vb_header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    vb_header.h.section_type          = SEC_VB_HEADER;
    vb_header.h.data_uncompressed_len = BGEN32 (my_squeeze_len);
    vb_header.h.compressed_offset     = BGEN32 (sizeof_header);
    vb_header.h.vblock_i              = BGEN32 (vb->vblock_i);
    vb_header.h.section_i             = BGEN16 (vb->z_next_header_i++); // always 0
    vb_header.h.sec_compression_alg  = COMP_PLN; // in versions 2-4 it was (needlessly) compressed with bzlib
    vb_header.num_lines               = BGEN32 ((uint32_t)vb->lines.len);
    vb_header.phase_type              = (char)vb->phase_type; 
    vb_header.has_genotype_data       = vb->has_genotype_data;
    vb_header.num_samples             = BGEN32 (global_vcf_num_samples);
    vb_header.num_haplotypes_per_line = BGEN32 (vb->num_haplotypes_per_line);
    vb_header.num_samples_per_block   = BGEN32 (vb->num_samples_per_block);
    vb_header.num_sample_blocks       = BGEN32 (vb->num_sample_blocks);
    vb_header.ploidy                  = BGEN16 (vb->ploidy);
    vb_header.vb_data_size            = BGEN32 (vb->vb_data_size);
    vb_header.max_gt_line_len         = BGEN32 (vb->max_gt_line_len);
    vb_header.is_gtshark              = flag_gtshark;

    // create squeezed index - IF we have haplotype data AND more than one haplotype per line (i.e. my_squeeze_len > 0)
    if (my_squeeze_len) {
        buf_alloc (vb, &vb->haplotype_permutation_index_squeezed, my_squeeze_len, 1, "haplotype_permutation_index_squeezed", 0);

        uint16_t haplotype_index_checksum;
        squeeze (vb, (uint8_t *)vb->haplotype_permutation_index_squeezed.data, &haplotype_index_checksum,
                (unsigned *)vb->haplotype_permutation_index.data, 
                vb->num_haplotypes_per_line);

        vb_header.haplotype_index_checksum = BGEN16(haplotype_index_checksum);
    }
    
    // compress section into z_data - to be eventually written to disk by the I/O thread
    vb->z_data.name = "z_data"; // this is the first allocation of z_data - comp_compress requires that it is pre-named
    vb->z_data.param = vb->vblock_i;
    comp_compress ((VBlockP)vb, &vb->z_data, false, (SectionHeader*)&vb_header, 
                    my_squeeze_len ? vb->haplotype_permutation_index_squeezed.data : NULL, NULL);

    // account for haplotype_index as part of the haplotype_data for compression statistics
    // even though we store it in the variant_data header
    vb->z_section_bytes[SEC_HT_DATA]   += my_squeeze_len;
    vb->z_section_bytes[SEC_VB_HEADER] -= my_squeeze_len;
}

// ZIP only: called by the I/O thread in the sequential order of VBs: updating of the already compressed
// variant data section (compressed by the compute thread in zfile_vcf_compress_vb_header) just before writing it to disk
// note: this updates the z_data in memory (not on disk)
void zfile_vcf_update_compressed_vb_header (VBlock *vb, uint32_t txt_first_line_i)
{
    SectionHeaderVbHeaderVCF *vb_header = (SectionHeaderVbHeaderVCF *)vb->z_data.data;
    vb_header->z_data_bytes = BGEN32 ((uint32_t)vb->z_data.len);
    vb_header->first_line   = BGEN32 (txt_first_line_i);

    if (flag_show_vblocks) {
        fprintf (stderr, "vb_i=%u first_line=%u num_lines=%u uncomprssed=%u compressed=%u num_sample_blocks=%u ht_per_line=%u max_gt_line_len=%u\n",
                 vb->vblock_i, txt_first_line_i, BGEN32 (vb_header->num_lines), 
                 BGEN32 (vb_header->vb_data_size), BGEN32 (vb_header->z_data_bytes), 
                 BGEN32 (vb_header->num_sample_blocks), BGEN32 (vb_header->num_haplotypes_per_line), 
                 BGEN32 (vb_header->max_gt_line_len));
    }
    // now we can finally encrypt the header - if needed
    if (crypt_have_password())  
        crypt_do (vb, (uint8_t*)vb_header, BGEN32 (vb_header->h.compressed_offset),
                  BGEN32 (vb_header->h.vblock_i), vb_header->h.section_type, true);
}

void zfile_vcf_compress_haplotype_data_gtshark (VBlockVCF *vb, const Buffer *haplotype_sections_data, unsigned sb_i)
{
    gtshark_compress_haplotype_data (vb, haplotype_sections_data, sb_i); // populates vb->gtshark_*
    
    vb->gtshark_exceptions_line_i.len *= sizeof (uint32_t);
    vb->gtshark_exceptions_ht_i.len *= sizeof (uint16_t);

    zfile_compress_section_data     ((VBlockP)vb, SEC_HT_GTSHARK_X_LINE,   &vb->gtshark_exceptions_line_i);
    zfile_compress_section_data     ((VBlockP)vb, SEC_HT_GTSHARK_X_HTI,    &vb->gtshark_exceptions_ht_i);
    zfile_compress_section_data     ((VBlockP)vb, SEC_HT_GTSHARK_X_ALLELE, &vb->gtshark_exceptions_allele);
    zfile_compress_section_data_alg ((VBlockP)vb, SEC_HT_GTSHARK_DB_DB,    &vb->gtshark_db_db_data, 0, 0, COMP_PLN);
    zfile_compress_section_data_alg ((VBlockP)vb, SEC_HT_GTSHARK_DB_GT,    &vb->gtshark_db_gt_data, 0, 0, COMP_PLN);

    // free buffers - they will be needed by the next section
    buf_free (&vb->gtshark_exceptions_line_i);
    buf_free (&vb->gtshark_exceptions_ht_i);
    buf_free (&vb->gtshark_exceptions_allele);
    buf_free (&vb->gtshark_db_db_data);
    buf_free (&vb->gtshark_db_gt_data);
}

#define READ_SB_SECTION(sec,header_type,sb_i) \
    { *ENT (unsigned, vb->z_section_headers, section_i++) = vb->z_data.len; \
      zfile_read_section ((VBlockP)vb, vb->vblock_i, sb_i, &vb->z_data, "z_data", sizeof(header_type), sec, sl++); }

#define READ_SECTION(sec,header_type) READ_SB_SECTION(sec, header_type, NO_SB_I)

void zfile_vcf_read_one_vb (VBlock *vb_)
{ 
    START_TIMER;

    VBlockVCF *vb = (VBlockVCF *)vb_;
    // The VB is read from disk here, in the I/O thread, and is decompressed in piz_uncompress_all_sections() in the 
    // Compute thread, with the exception of dictionaries that are handled here - this VBs dictionary fragments are
    // integrated into the global dictionaries.
    // Order of sections in a V2 VB:
    // 1. SEC_VB_HEADER - its data is the haplotype index
    // 2. SEC_VCF_INFO_SF_B250 - Fields 1-9 b250 data
    // 3. SEC_VCF_INFO_SF_B250 - All INFO subfield data
    // 4. All sample data: up 3 sections per sample block:
    //    4a. SEC_VCF_GT_DATA - genotype data
    //    4b. SEC_VCF_PHASE_DATA - phase data
    //    4c. SEC_HT_DATA  or SEC_HAPLOTYPE_GTSHARK - haplotype data

    SectionListEntry *sl = sections_vb_first (vb->vblock_i);

    int vb_header_offset = zfile_read_section ((VBlockP)vb, vb->vblock_i, NO_SB_I, &vb->z_data, "z_data",
                                               sizeof(SectionHeaderVbHeaderVCF), SEC_VB_HEADER, sl++);

    // note - use a macro and not a variable bc vb_header changes when z_data gets realloced as we read more data
    #define vb_header ((SectionHeaderVbHeaderVCF *)&vb->z_data.data[vb_header_offset])

    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);

    // overlay all dictionaries (not just those that have fragments in this variant block) to the vb
    mtf_overlay_dictionaries_to_vb ((VBlockP)vb);
    
    // iname mapper: first VB, we map all inname subfields to this global mapper, and each VB overlays it on its own map
    if (!flag_strip) {
    
        static Buffer global_iname_mapper_buf = EMPTY_BUFFER;
    
        if (vb->vblock_i == 1) {
            buf_free (&global_iname_mapper_buf); // in case it was used to piz a previous file
            (z_file->genozip_version >= 4 ? piz_vcf_map_iname_subfields : v2v3_piz_vcf_map_iname_subfields)(&global_iname_mapper_buf);
            buf_set_overlayable (&global_iname_mapper_buf);
        }
        
        buf_overlay (vb_, &vb->iname_mapper_buf, &global_iname_mapper_buf, "iname_mapper_buf", vb->vblock_i);    

        vb->end_did_i = mtf_get_existing_did_i_by_dict_id ((DictIdType)dict_id_INFO_END); 
    }

    // read the the data sections (fields, info sub fields, genotype, phase, haplotype)

    buf_alloc (vb, &vb->z_section_headers, 1000 /* arbitrary initial value */ * sizeof(char*), 0, "z_section_headers", 1);
    
    ((unsigned *)vb->z_section_headers.data)[0] = vb_header_offset; // variant data header is at index 0

    unsigned section_i=1;

    // read the 8 fields (CHROM to FORMAT)    
    for (VcfFields f=VCF_CHROM; f <= VCF_FORMAT; f++)
        READ_SECTION (FIELD_TO_B250_SECTION(DT_VCF, f), SectionHeaderBase250);

    // read the info subfield sections into memory (if any)
    vb->num_info_subfields = sections_count_sec_type (vb->vblock_i, SEC_VCF_INFO_SF_B250); // also used later in piz_uncompress_all_sections()
    for (uint8_t sf_i=0; sf_i < vb->num_info_subfields; sf_i++) 
        READ_SECTION (SEC_VCF_INFO_SF_B250, SectionHeaderBase250);

    // read the numberic data of the ID field (the non-numeric part is in SEC_ID_B250)
    if (z_file->genozip_version >= 5)
        READ_SECTION (SEC_NUMERIC_ID_DATA, SectionHeader);
        
    // read the sample data
    uint32_t num_sample_blocks     = BGEN32 (vb_header->num_sample_blocks);
    uint32_t num_samples_per_block = BGEN32 (vb_header->num_samples_per_block);

    buf_alloc (vb, &vb->is_sb_included, num_sample_blocks * sizeof(bool), 1, "is_sb_included", vb->vblock_i);

    for (unsigned sb_i=0; sb_i < num_sample_blocks; sb_i++) {

        // calculate whether this block is included. zfile_read_section will skip reading and piz_uncompress_all_sections
        // will skip uncompressing based on this value
        NEXTENT (bool, vb->is_sb_included) = samples_is_sb_included (num_samples_per_block, sb_i);
 
        // make sure we have enough space for the section pointers
        buf_alloc (vb, &vb->z_section_headers, sizeof (uint32_t) * (section_i + 3), 2, "z_section_headers", 2);

        if (vb_header->has_genotype_data)
            READ_SB_SECTION (SEC_VCF_GT_DATA,         SectionHeader, sb_i);

        if (vb_header->phase_type == PHASE_MIXED_PHASED) 
            READ_SB_SECTION (SEC_VCF_PHASE_DATA,      SectionHeader, sb_i);

        if (vb_header->num_haplotypes_per_line != 0 && !vb_header->is_gtshark) 
            READ_SB_SECTION (SEC_HT_DATA,             SectionHeader, sb_i);

        if (vb_header->num_haplotypes_per_line != 0 && vb_header->is_gtshark) {
            READ_SB_SECTION (SEC_HT_GTSHARK_X_LINE,   SectionHeader, sb_i);
            READ_SB_SECTION (SEC_HT_GTSHARK_X_HTI,    SectionHeader, sb_i);
            READ_SB_SECTION (SEC_HT_GTSHARK_X_ALLELE, SectionHeader, sb_i);
            READ_SB_SECTION (SEC_HT_GTSHARK_DB_DB,    SectionHeader, sb_i);
            READ_SB_SECTION (SEC_HT_GTSHARK_DB_GT,    SectionHeader, sb_i);
        }
    }
    
    vb->ready_to_dispatch = true; // all good

    COPY_TIMER (vb->profile.zfile_read_one_vb);

    #undef vb_header
}

#define V1_ZFILE // select the zfile functions of v1.c
#include "v1_vcf.c"

//------------------------------------------------------------------------------
// SAM stuff
//------------------------------------------------------------------------------

void zfile_sam_read_one_vb (VBlock *vb_)
{ 
    START_TIMER;

    VBlockSAM *vb = (VBlockSAM *)vb_;

    // The VB is read from disk here, in the I/O thread, and is decompressed in piz_uncompress_all_sections() in the 
    // Compute thread, with the exception of dictionaries that are handled here - this VBs dictionary fragments are
    // integrated into the global dictionaries.
    // Order of sections in a VB:
    // 1. SEC_VB_HEADER - its data is the haplotype index
    // 2. SEC_SAM_QNAME_B250 to SEC_SAM_OPTIONAL_B250- Main b250 data
    // 3. SEC_SAM_QNAME_SF_B250 - QNAME subfields
    // 4. SEC_SAM_OPTNL_SF_B250 - OPTIONAL subfields
    // 5. SEC_SAM_RAND_POS_DATA - Random POS data - for non-delta pos data from POS, PNEXT, SA, OA, XA
    // 6. SEC_SAM_MD_DATA       - Data for the optional MD field
    // 7. SEC_SEQ_DATA      - Sequences data from SEQ, E2
    // 8. SEC_QUAL_DATA     - Quality data from QUAL, U2    

    SectionListEntry *sl = sections_vb_first (vb->vblock_i);
    
    int vb_header_offset = zfile_read_section ((VBlockP)vb, vb->vblock_i, NO_SB_I, &vb->z_data, "z_data",
                                               sizeof(SectionHeaderVbHeader), SEC_VB_HEADER, sl++);

    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);

    // overlay all dictionaries (not just those that have fragments in this vblock) to the vb
    mtf_overlay_dictionaries_to_vb ((VBlockP)vb);

    // read the the data sections (fields, QNAME and OPTIONAL subfields, SEQ and QUAL data)

    // room for section headers (we have at most MAX_DICTS + 5 as all sections are b250 except for VB header, SEQ,QUAL,MD,RAND_POS)
    buf_alloc (vb, &vb->z_section_headers, (MAX_DICTS + 5)  * sizeof(char*), 0, "z_section_headers", 1);
    
    *FIRSTENT (unsigned, vb->z_section_headers) = vb_header_offset; // vb header is at index 0

    unsigned section_i=1;

    // read the fields
    for (SamFields f=SAM_QNAME; f <= SAM_OPTIONAL; f++) 
        READ_SECTION (FIELD_TO_B250_SECTION(z_file->data_type, f), SectionHeaderBase250);

    // read the QNAME subfields sections into memory (if any)
    vb->qname_mapper.num_subfields = sections_count_sec_type (vb->vblock_i, SEC_SAM_QNAME_SF_B250); 
    for (uint8_t sf_i=0; sf_i < vb->qname_mapper.num_subfields; sf_i++) 
        READ_SECTION (SEC_SAM_QNAME_SF_B250, SectionHeaderBase250);

    // read the OPTIONAL subfields sections into memory (if any)
    vb->num_optional_subfield_b250s = sections_count_sec_type (vb->vblock_i, SEC_SAM_OPTNL_SF_B250);
    for (uint8_t sf_i=0; sf_i < vb->num_optional_subfield_b250s; sf_i++) 
        READ_SECTION (SEC_SAM_OPTNL_SF_B250, SectionHeaderBase250);

    // read the RAND_POS, MD, SEQ and QUAL data
    READ_SECTION (SEC_SAM_RAND_POS_DATA, SectionHeader);
    READ_SECTION (SEC_SAM_MD_DATA,       SectionHeader);
    READ_SECTION (SEC_SEQ_DATA,          SectionHeader);
    READ_SECTION (SEC_QUAL_DATA,         SectionHeader);

    vb->ready_to_dispatch = true; // all good

    COPY_TIMER (vb->profile.zfile_read_one_vb);
}

//------------------------------------------------------------------------------
// FASTQ & FASTA stuff
//------------------------------------------------------------------------------

void zfile_fast_read_one_vb (VBlock *vb_)
{ 
    START_TIMER;

    VBlockFAST *vb = (VBlockFAST *)vb_;

    // The VB is read from disk here, in the I/O thread, and is decompressed in piz_uncompress_all_sections() in the 
    // Compute thread, with the exception of dictionaries that are handled here - this VBs dictionary fragments are
    // integrated into the global dictionaries.
    // Order of sections in a VB:
    // 1. SEC_VB_HEADER - its data is the haplotype index
    // 2. SEC_FAST_DESC_B250 and SEC_FAST_LINEMETA_B250
    // 3. SEC_FAST_DESC_SF_B250 - Description subfields
    // 7. SEC_SEQ_DATA      - Sequences data 
    // 8. SEC_QUAL_DATA     - Quality data (FASTQ only)    

    SectionListEntry *sl = sections_vb_first (vb->vblock_i);
    
    int vb_header_offset = zfile_read_section ((VBlockP)vb, vb->vblock_i, NO_SB_I, &vb->z_data, "z_data",
                                               sizeof(SectionHeaderVbHeader), SEC_VB_HEADER, sl++);

    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);

    // overlay all dictionaries (not just those that have fragments in this vblock) to the vb
    mtf_overlay_dictionaries_to_vb ((VBlockP)vb);

    // read the the data sections (fields, QNAME and OPTIONAL subfields, SEQ and QUAL data)

    // room for section headers (we have at most MAX_DICTS + 3 as all sections are b250 except for VB header, SEQ and QUAL/COMMENT)
    buf_alloc (vb, &vb->z_section_headers, (MAX_DICTS + 3)  * sizeof(char*), 0, "z_section_headers", 1);
    
    *FIRSTENT (unsigned, vb->z_section_headers) = vb_header_offset; // vb header is at index 0

    unsigned section_i=1;

    // read the field sections
    READ_SECTION (SEC_FAST_DESC_B250, SectionHeaderBase250);
    READ_SECTION (SEC_FAST_LINEMETA_B250, SectionHeaderBase250);
    
    // read the DESC subfields sections into memory (if any)
    vb->desc_mapper.num_subfields = sections_count_sec_type (vb->vblock_i, SEC_FAST_DESC_SF_B250); 
    for (uint8_t sf_i=0; sf_i < vb->desc_mapper.num_subfields; sf_i++) 
        READ_SECTION (SEC_FAST_DESC_SF_B250, SectionHeaderBase250);

    if (flag_grep && !piz_fast_test_grep (vb)) goto finish; // ususually, we uncompress and reconstruct the DESC from the I/O thread in case of --grep

    // read SEQ and QUAL data (FASTQ) or COMMENT data (FASTA)
    READ_SECTION (SEC_SEQ_DATA, SectionHeader);
    READ_SECTION (z_file->data_type == DT_FASTQ ? SEC_QUAL_DATA : SEC_FASTA_COMMENT_DATA, SectionHeader);

    vb->ready_to_dispatch = true; // all good

finish:
    COPY_TIMER (vb->profile.zfile_read_one_vb);
}

//------------------------------------------------------------------------------
// ME23 stuff
//------------------------------------------------------------------------------

void zfile_me23_read_one_vb (VBlock *vb_)
{ 
    START_TIMER;

    VBlockME23 *vb = (VBlockME23 *)vb_;

    SectionListEntry *sl = sections_vb_first (vb->vblock_i);
    
    int vb_header_offset = zfile_read_section ((VBlockP)vb, vb->vblock_i, NO_SB_I, &vb->z_data, "z_data",
                                               sizeof(SectionHeaderVbHeader), SEC_VB_HEADER, sl++);

    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);

    // overlay all dictionaries (not just those that have fragments in this vblock) to the vb
    mtf_overlay_dictionaries_to_vb ((VBlockP)vb);

    // read the the data sections (fields, QNAME and OPTIONAL subfields, SEQ and QUAL data)

    // room for section headers - we have 5 sections (inc. the VB header)
    buf_alloc (vb, &vb->z_section_headers, 5 * sizeof(char*), 0, "z_section_headers", 1);
    
    *FIRSTENT (unsigned, vb->z_section_headers) = vb_header_offset; // vb header is at index 0

    unsigned section_i=1;

    // read the sections
    READ_SECTION (SEC_CHROM_B250, SectionHeaderBase250);
    READ_SECTION (SEC_POS_B250, SectionHeaderBase250);
    READ_SECTION (SEC_ID_B250, SectionHeaderBase250);
    READ_SECTION (SEC_NUMERIC_ID_DATA, SectionHeader);
    READ_SECTION (SEC_HT_DATA, SectionHeader);

    COPY_TIMER (vb->profile.zfile_read_one_vb);

    vb->ready_to_dispatch = true; // all good
}