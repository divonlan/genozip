// ------------------------------------------------------------------
//   sections.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "sections.h"
#include "buffer.h"
#include "file.h"
#include "vb.h"
#include "endianness.h"

// ZIP only: create section list that goes into the genozip header, as we are creating the sections
void sections_add_to_list (VariantBlock *vb, const SectionHeader *header)
{
    DictIdType dict_id = { 0 };
    bool is_dict = (section_type_is_dictionary (header->section_type));

    if (is_dict)                                          dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    else if (section_type_is_b250 (header->section_type)) dict_id = ((SectionHeaderBase250    *)header)->dict_id;

    // 1. if this is a vcf_header, random_access or genozip_header - it goes directly into the z_file by the I/O thread
    //    before or after all the compute threads are operational
    // 2. if this is a dictionary - it goes directly into z_file by the Compute thread while merge holds the mutex:
    //    mtf_merge_in_vb_ctx_one_dict_id -> zfile_compress_dictionary_data
    // 3. if we this section is part of a VB (other than a dictionary), we store the entry within the VB and merge it to
    //    the zfile in the correct order of VBs by the I/O thread after the compute thread is finished.
    //
    // offsets in case 2 and 3 are relative to their buffer at this point, and will be updated later

    Buffer *section_list_buf;
    uint64_t offset;
    if (!vb->variant_block_i) {  // case 1
        section_list_buf = &vb->z_file->section_list_buf;
        offset = vb->z_file->disk_so_far + vb->z_data.len;
    }
    else if (is_dict) {          // case 2
        section_list_buf = &vb->z_file->section_list_dict_buf;
        offset = vb->z_file->dict_data.len;
    }
    else {                       // case 3
        section_list_buf = &vb->section_list_buf;
        offset = vb->z_data.len;
    }

    buf_alloc (vb, section_list_buf, MAX (section_list_buf->len + 1, 50) * sizeof(SectionListEntry), 2, 
               "section_list_buf", vb->variant_block_i);

    SectionListEntry *ent = &((SectionListEntry *)section_list_buf->data)[section_list_buf->len++];
    ent->section_type    = header->section_type;
    ent->variant_block_i = BGEN32 (header->variant_block_i); // big endian in header - convert back to native
    ent->dict_id         = dict_id;
    ent->offset          = offset;  // this is a partial offset (within d) - we will correct it later
    ent->include         = 0;
    ent->for_future_use  = 0;
}

// Called by ZIP I/O thread. concatenates a vb or dictionary section list to the z_file sectinon list - just before 
// writing those sections to the disk. we use the current disk position to update the offset
void sections_list_concat (VariantBlock *vb, BufferP section_list_buf)
{
    File *zfile = vb->z_file;

    buf_alloc (vb, &zfile->section_list_buf, 
              (zfile->section_list_buf.len + section_list_buf->len) * sizeof(SectionListEntry), 2, 
              "section_list_buf", 0);
  
    SectionListEntry *dst = &((SectionListEntry *)zfile->section_list_buf.data)[zfile->section_list_buf.len];
    SectionListEntry *src = ((SectionListEntry *)section_list_buf->data);

    // update the offset
    for (unsigned i=0; i < section_list_buf->len; i++)
        src[i].offset += zfile->disk_so_far;

    // copy all entries
    memcpy (dst, src, section_list_buf->len * sizeof(SectionListEntry));
    zfile->section_list_buf.len += section_list_buf->len;

    buf_free (section_list_buf);
}

// called by PIZ I/O thread: zfile_read_on_vb
uint32_t sections_count_info_b250s (File *z_file)
{
    SectionListEntry *sl = ((SectionListEntry *)z_file->section_list_buf.data);

    // skip to the first SEC_INFO_SUBFIELD_B250
    while (sl[z_file->sl_cursor].section_type != SEC_INFO_SUBFIELD_B250) z_file->sl_cursor++;

    // count the SEC_INFO_SUBFIELD_B250 sections
    uint32_t start = z_file->sl_cursor;
    while (sl[z_file->sl_cursor].section_type == SEC_INFO_SUBFIELD_B250) z_file->sl_cursor++;

    return z_file->sl_cursor - start;
}

// called by PIZ I/O to know if next up is a VB Header or VCF Header or EOF
SectionType sections_get_next_header_type (File *z_file)
{
    SectionListEntry *sl = ((SectionListEntry *)z_file->section_list_buf.data);

    // find the next VB or VCF header section
    while (z_file->sl_cursor < z_file->section_list_buf.len) {
        SectionType sec_type = sl[z_file->sl_cursor++].section_type;
        if (sec_type == SEC_VB_HEADER || sec_type == SEC_VCF_HEADER) return sec_type;
    }

    return SEC_EOF; // no more headers
}

// called by PIZ I/O when splitting a concatenated file - to know if there are any more VCF components remaining
bool sections_has_more_vcf_components (File *z_file)
{
    SectionListEntry *sl = ((SectionListEntry *)z_file->section_list_buf.data);

    return z_file->sl_cursor==0 || sl[z_file->sl_cursor-1].section_type == SEC_VCF_HEADER;
}

void sections_show_genozip_header (VariantBlock *pseudo_vb, SectionHeaderGenozipHeader *header)
{
    unsigned num_sections = BGEN32 (header->num_sections);
    char size_str[50];

    printf ("Below are the contents of the genozip header. This is the output of --show-gheader:\n");
    printf ("  genozip_version: %u\n",               header->genozip_version);
    printf ("  data_type: %s\n",                     dt_name (BGEN16 (header->data_type)));
    printf ("  encryption_type: %s\n",               encryption_name (header->encryption_type)); 
    printf ("  num_samples: %u\n",                   BGEN32 (header->num_samples));
    printf ("  uncompressed_data_size: %s\n",        buf_human_readable_uint (BGEN64 (header->uncompressed_data_size), size_str));
    printf ("  num_items_concat: %"PRIu64"\n",       BGEN64 (header->num_items_concat));
    printf ("  num_sections: %u\n",                  num_sections);
    printf ("  num_vcf_components: %u\n",            BGEN32 (header->num_vcf_components));
    printf ("  md5_hash_concat: %s\n",               md5_display (&header->md5_hash_concat, false));
    printf ("  created: %*s\n",                      -FILE_METADATA_LEN, header->created);

    printf ("  sections:\n");

    SectionListEntry *ents = (SectionListEntry *)pseudo_vb->z_file->section_list_buf.data;

    for (unsigned i=0; i < num_sections; i++) {
     
        uint64_t this_offset = BGEN64 (ents[i].offset);
        uint64_t next_offset = (i < num_sections-1) ? BGEN64 (ents[i+1].offset) : pseudo_vb->z_file->disk_so_far;

        printf ("    %3u. %-24.24s %*.*s vb_i=%u offset=%"PRIu64" size=%"PRId64"\n", 
                 i, st_name(ents[i].section_type), 
                 -DICT_ID_LEN, DICT_ID_LEN, ents[i].dict_id.num ? dict_id_printable (ents[i].dict_id).id : ents[i].dict_id.id, 
                 BGEN32 (ents[i].variant_block_i), this_offset, next_offset - this_offset);
    }
}

void BGEN_sections_list (Buffer *sections_list_buf)
{
    SectionListEntry *ent = (SectionListEntry *)sections_list_buf->data;

    // update only entries that belong to the same vb AND have the same class - dictionary or not dictionary
    for (unsigned i=0; i < sections_list_buf->len; i++) {
        ent[i].variant_block_i = BGEN32 (ent[i].variant_block_i);
        ent[i].offset          = BGEN64 (ent[i].offset);
    }
}

static const char *sections_type_name (unsigned item, const char **names, unsigned num_names)
{
    if (item > num_names) {
        static char str[50];
        sprintf (str, "%u (out of range)", item);
        return str;
    }
    
    return names[item];    
}

const char *st_name(unsigned sec_type)
{
    static const char *names[] = SECTIONTYPE_NAMES;
    return sections_type_name (sec_type, names, sizeof(names)/sizeof(names[0]));
}

const char *dt_name (unsigned data_type)
{
    static const char *names[] = DATATYPE_NAMES;
    return sections_type_name (data_type, names, sizeof(names)/sizeof(names[0]));
}

const char *encryption_name (unsigned encryption_type)
{
    static const char *names[] = ENCRYPTION_TYPE_NAMES;
    return sections_type_name (encryption_type, names, sizeof(names)/sizeof(names[0]));
}