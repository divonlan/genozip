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
    ent->include = ent->for_future_use = 0;
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

void sections_show_genozip_header (VariantBlock *pseudo_vb, SectionHeaderGenozipHeader *header)
{
    unsigned num_sections = BGEN32 (header->num_sections);

    fprintf (stderr, "The genozip header follows, the result of --show-gheader:\n");
    fprintf (stderr, "  genozip_version: %u\n",               header->genozip_version);
    fprintf (stderr, "  encryption_type: %u\n",               header->encryption_type); 
    fprintf (stderr, "  data_type: %u\n",                     BGEN16 (header->data_type));
    fprintf (stderr, "  num_samples: %u\n",                   BGEN32 (header->num_samples));
    fprintf (stderr, "  uncompressed_data_size: %"PRIu64"\n", BGEN64 (header->uncompressed_data_size));
    fprintf (stderr, "  num_items_concat: %"PRIu64"\n",       BGEN64 (header->num_items_concat));
    fprintf (stderr, "  num_sections: %u\n",                  num_sections);
    fprintf (stderr, "  num_vcf_components: %u\n",            BGEN32 (header->num_vcf_components));
    fprintf (stderr, "  num_info_dictionary_sections: %u\n",  BGEN32 (header->num_info_dictionary_sections));
    fprintf (stderr, "  num_gt_dictionary_sections: %u\n",    BGEN32 (header->num_gt_dictionary_sections));
    fprintf (stderr, "  md5_hash_concat: %s\n",               md5_display (&header->md5_hash_concat, false));
    fprintf (stderr, "  created: %*s\n",                      -FILE_METADATA_LEN, header->created);

    fprintf (stderr, "  sections:\n");

    SectionListEntry *ents = (SectionListEntry *)pseudo_vb->z_file->section_list_buf.data;

    for (unsigned i=0; i < num_sections; i++)
        fprintf (stderr, "    %3u. %-22.22s %*.*s vb_i=%u offset_on_disk=%"PRIu64"\n", 
                 i, st_name(ents[i].section_type), -DICT_ID_LEN, DICT_ID_LEN, 
                 ents[i].dict_id.num ? dict_id_printable (ents[i].dict_id).id : ents[i].dict_id.id, 
                 BGEN32 (ents[i].variant_block_i), BGEN64 (ents[i].offset));
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

const char *st_name(unsigned sec_type)
{
    static const char *st_names[] = SECTIONTYPE_NAMES;

    if (sec_type > sizeof(st_names)/sizeof(st_names[0])) {
        static char str[50];
        sprintf (str, "%u (out of range)", sec_type);
        return str;
    }
    
    return st_names[sec_type];
}

