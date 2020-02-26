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

void sections_add_to_list (VariantBlock *pseudo_vb, const SectionHeader *header, uint64_t offset)
{
    DictIdType dict_id = { 0 };
    if (section_type_is_dictionary (header->section_type)) dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    if (section_type_is_b250       (header->section_type)) dict_id = ((SectionHeaderBase250    *)header)->dict_id;

    buf_alloc (pseudo_vb, &pseudo_vb->z_file->section_list_buf, MAX (pseudo_vb->z_file->section_list_buf.len + 1, 100) * sizeof(SectionListEntry), 
               2, "z_file->section_list_buf", 0);

    SectionListEntry *ent = &((SectionListEntry *)pseudo_vb->z_file->section_list_buf.data)[pseudo_vb->z_file->section_list_buf.len++];
    ent->section_type    = header->section_type;
    ent->variant_block_i = BGEN32 (header->variant_block_i); // big endian in header - convert back to native
    ent->dict_id         = dict_id;
    ent->offset          = offset;
    ent->include = ent->for_future_use = 0;
}

void sections_update_list (VariantBlock *vb, bool is_z_data)
{
    File *zfile = vb->z_file;

    SectionListEntry *ent = (SectionListEntry *)zfile->section_list_buf.data;

    // update only entries that belong to the same vb AND have the same class - dictionary or not dictionary
    for (unsigned i=0; i < zfile->section_list_buf.len; i++) {
        bool is_dict = section_type_is_dictionary (ent[i].section_type);
        if ((ent[i].variant_block_i == vb->variant_block_i && is_z_data && !is_dict) || 
            (!is_z_data && is_dict)) {
        
            ent[i].offset += zfile->disk_so_far;
            printf ("New offset of %s vb_i=%u offset=%u\n", st_name(ent[i].section_type), ent[i].variant_block_i, (unsigned)ent[i].offset);
        }
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

