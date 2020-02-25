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

    buf_alloc (pseudo_vb, &pseudo_vb->z_file->section_list_buf, MAX (pseudo_vb->z_file->section_list_buf.len + 1, 100) * sizeof(SectionListEntry), 2, "", 0);

    SectionListEntry *ent = &((SectionListEntry *)pseudo_vb->z_file->section_list_buf.data)[pseudo_vb->z_file->section_list_buf.len++];
    ent->section_type    = header->section_type;
    ent->variant_block_i = header->variant_block_i; // already big endian
    ent->dict_id         = dict_id;
    ent->offset          = BGEN64 (offset);
    ent->include = ent->for_future_use = 0;
}
