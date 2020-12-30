// ------------------------------------------------------------------
//   sections.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "sections.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "endianness.h"
#include "random_access.h"
#include "strings.h"
#include "crypt.h"
#include "dict_id.h"
#include "zfile.h"

static const struct {const char *name; uint32_t header_size; } abouts[NUM_SEC_TYPES] = SECTIONTYPE_ABOUT;

const LocalTypeDesc lt_desc[NUM_LOCAL_TYPES] = LOCALTYPE_DESC;

// ZIP only: create section list that goes into the genozip header, as we are creating the sections. returns offset
void sections_add_to_list (VBlock *vb, const SectionHeader *header)
{
    DictId dict_id = DICT_ID_NONE;
    bool has_dict_id = true;

    if      (header->section_type == SEC_DICT ) dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    else if (header->section_type == SEC_B250 ) dict_id = ((SectionHeaderCtx        *)header)->dict_id;
    else if (header->section_type == SEC_LOCAL) dict_id = ((SectionHeaderCtx        *)header)->dict_id;
    else has_dict_id = false;

    ASSERTE0 (!has_dict_id || dict_id.num, "dict_id=0");
    
    buf_alloc_more (vb, &vb->section_list_buf, 1, 50, SectionListEntry, 2, "section_list_buf");
    
    NEXTENT (SectionListEntry, vb->section_list_buf) = (SectionListEntry) {
        .section_type = header->section_type,
        .vblock_i     = BGEN32 (header->vblock_i), // big endian in header - convert back to native
        .dict_id      = dict_id,
        .offset       = vb->z_data.len  // this is a partial offset (within d) - we will correct it later
    };
}

// Called by ZIP I/O thread. concatenates a vb or dictionary section list to the z_file section list - just before 
// writing those sections to the disk. we use the current disk position to update the offset
void sections_list_concat (VBlock *vb)
{
    if (!vb->section_list_buf.len) return;

    buf_alloc_more (evb, &z_file->section_list_buf, vb->section_list_buf.len, 0, SectionListEntry, 2, "z_file->section_list_buf");

    SectionListEntry *dst = AFTERENT (SectionListEntry, z_file->section_list_buf);
    SectionListEntry *src = FIRSTENT (SectionListEntry, vb->section_list_buf);

    // update the offset
    for (unsigned i=0; i < vb->section_list_buf.len; i++)
        src[i].offset += z_file->disk_so_far;

    // copy all entries
    memcpy (dst, src, vb->section_list_buf.len * sizeof(SectionListEntry));
    z_file->section_list_buf.len += vb->section_list_buf.len;

    buf_free (&vb->section_list_buf);
}

// section iterator. returns true if a section of this type was found.
bool sections_get_next_section_of_type2 (const SectionListEntry **sl_ent, // optional in/out. if NULL - search entire list
                                         SectionType st1, SectionType st2, 
                                         bool must_be_next_section,       // check only next section, not entire remaining list
                                         bool seek)                       // if true, seek if found
{
    const SectionListEntry *sl = sl_ent ? *sl_ent : NULL; 
    bool found = false;

    // case: first time or we are allowed skip sections
    if (!sl || !must_be_next_section) {

        while (sl < AFTERENT (const SectionListEntry, z_file->section_list_buf) - 1) {

            sl = sl ? (sl + 1) : FIRSTENT (const SectionListEntry, z_file->section_list_buf); 

            if (sl->section_type == st1 || sl->section_type == st2) {
                found = true;
                break;
            }
        }
    }
    // case: we aren't allowed to skip sections - check if the next section is what we were after
    else {
        sl++;
        found = sl->section_type == st1 || sl->section_type == st2;
    }

    if (found && seek)
        file_seek (z_file, sl->offset, SEEK_SET, false);

    if (sl_ent) *sl_ent = sl;
    return found;
}

// find the first and last vb_i of a the immediately previous bound file, start from an sl in this file
void sections_get_prev_component_vb_i (const SectionListEntry *sl, // any sl of the current file
                                  uint32_t *prev_file_first_vb_i, uint32_t *prev_file_last_vb_i) //out
{
#   define SAFTEY ASSERT0 ((char*)sl > z_file->section_list_buf.data, "Error in sections_get_prev_component_vb_i: cannot find previous file VBs")

    // search back to current file's txt header
    do { SAFTEY; sl--; } while (sl->section_type != SEC_TXT_HEADER);

    // search back to previous file's last VB header
    do { SAFTEY; sl--; } while (sl->section_type != SEC_VB_HEADER);
    *prev_file_last_vb_i = sl->vblock_i;

    // search back to the previous file's txt header
    do { SAFTEY; sl--; } while (sl->section_type != SEC_TXT_HEADER);

    sl++;
    *prev_file_first_vb_i = sl->vblock_i;

#   undef SAFETY
}

// count how many sections we have of a certain type
uint32_t sections_count_sections (SectionType st)
{
    uint32_t count=0;
    for (uint32_t i=0; i < z_file->section_list_buf.len; i++) 
        if (ENT (SectionListEntry, z_file->section_list_buf, i)->section_type == st)
            count++;

    return count;
}

// called by PIZ I/O : vcf_zfile_read_one_vb. Sets *sl_ent to the first section of this vb_i, and returns its offset
const SectionListEntry *sections_vb_first (uint32_t vb_i, bool soft_fail)
{
    const SectionListEntry *sl=NULL;
    unsigned i=0; for (; i < z_file->section_list_buf.len; i++) {
        sl = ENT (const SectionListEntry, z_file->section_list_buf, i);
        if (sl->vblock_i == vb_i) break; // found!
    }

    if (i >= z_file->section_list_buf.len) {
        if (soft_fail) return NULL;
        ABORT ("Error in sections_get_next_vb_section: cannot find any section for vb_i=%u", vb_i);
    }

    return sl;
}

// I/O thread: called by refhash_initialize - get details of the refhash ahead of loading it from the reference file 
void sections_get_refhash_details (uint32_t *num_layers, uint32_t *base_layer_bits) // optional outs
{
    ASSERT0 (flag.reading_reference || flag.show_ref_hash, 
             "Error in sections_get_refhash_details: can only be called while reading reference");

    for (int i=z_file->section_list_buf.len-1; i >= 0; i--) { // search backwards as the refhash sections are near the end
        const SectionListEntry *sl = ENT (const SectionListEntry, z_file->section_list_buf, i);
        if (sl->section_type == SEC_REF_HASH) {

            SectionHeaderRefHash *header = (SectionHeaderRefHash *)zfile_read_section_header (evb, sl->offset, 0, SEC_REF_HASH);
            if (num_layers) *num_layers = header->num_layers;
            if (base_layer_bits) *base_layer_bits = header->layer_bits + header->layer_i; // layer_i=0 is the base layer, layer_i=1 has 1 bit less etc

            buf_free (&evb->compressed); // allocated by zfile_read_section_header
            return;
        }
        else if (sl->section_type == SEC_REFERENCE)
            break; // we arrived at a SEC_REFERENCE - there won't be any more SEC_REF_HASH sections
    }

    ABORT ("Error in sections_get_refhash_details: can't find SEC_REF_HASH sections in %s", z_name);
}

void BGEN_sections_list()
{
    ARRAY (SectionListEntry, ent, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++) {
        ent[i].vblock_i = BGEN32 (ent[i].vblock_i);
        ent[i].offset   = BGEN64 (ent[i].offset);
    }
}

void sections_show_gheader (const SectionHeaderGenozipHeader *header)
{
    if (flag.reading_reference) return; // don't show gheaders of reference file
    
    unsigned num_sections = BGEN32 (header->num_sections);

    fprintf (info_stream, "Contents of the genozip header (output of --show-gheader) of %s:\n", z_name);
    fprintf (info_stream, "  genozip_version: %u\n",         header->genozip_version);
    fprintf (info_stream, "  data_type: %s\n",               dt_name (BGEN16 (header->data_type)));
    fprintf (info_stream, "  encryption_type: %s\n",         encryption_name (header->encryption_type)); 
    fprintf (info_stream, "  uncompressed_data_size: %s\n",  str_uint_commas (BGEN64 (header->uncompressed_data_size)).s);
    fprintf (info_stream, "  num_items_bound: %"PRIu64"\n", BGEN64 (header->num_items_bound));
    fprintf (info_stream, "  num_sections: %u\n",            num_sections);
    fprintf (info_stream, "  num_components: %u\n",          BGEN32 (header->num_components));
    fprintf (info_stream, "  digest_bound.md5: %s\n",        digest_display (header->digest_bound).s);
    fprintf (info_stream, "  created: %*s\n",                -FILE_METADATA_LEN, header->created);
    fprintf (info_stream, "  license_hash: %s\n",            digest_display (header->license_hash).s);
    fprintf (info_stream, "  reference filename: %*s\n",     -REF_FILENAME_LEN, header->ref_filename);
    fprintf (info_stream, "  reference file hash: %s\n",     digest_display (header->ref_file_md5).s);

    fprintf (info_stream, "  sections:\n");

    ARRAY (const SectionListEntry, ents, z_file->section_list_buf);

    for (unsigned i=0; i < num_sections; i++) {
     
        uint64_t this_offset = ents[i].offset;
        uint64_t next_offset;
        
        if (i < num_sections-1) 
            next_offset = ents[i+1].offset;

        else // we're at the last section genozip header+footer
            next_offset = this_offset + BGEN32 (header->h.data_compressed_len) + BGEN32 (header->h.compressed_offset) + sizeof (SectionFooterGenozipHeader);

        fprintf (info_stream, "    %3u. %-24.24s %*.*s vb_i=%u offset=%"PRIu64" size=%"PRId64"\n", 
                 i, st_name(ents[i].section_type), 
                 -DICT_ID_LEN, DICT_ID_LEN, ents[i].dict_id.num ? dict_id_printable (ents[i].dict_id).id : ents[i].dict_id.id, 
                 ents[i].vblock_i, this_offset, next_offset - this_offset);
    }
}

const char *st_name (SectionType sec_type)
{
    ASSERT (sec_type >= SEC_NONE && sec_type < NUM_SEC_TYPES, "Error in st_name: sec_type=%u out of range [-1,%u]", sec_type, NUM_SEC_TYPES-1);

    return (sec_type == SEC_NONE) ? "SEC_NONE" : type_name (sec_type, &abouts[sec_type].name , sizeof(abouts)/sizeof(abouts[0]));
}

// called to parse the optional argument to --show-headers. we accept eg "REFERENCE" or "SEC_REFERENCE" or even "ref"
SectionType sections_st_by_name (char *name)
{
    if (!name) return -2; // all sections

    str_toupper (name, name); // name is case insesitive (compare uppercase)

    for (SectionType st=0; st < NUM_SEC_TYPES; st++)
        if (strstr (abouts[st].name, name)) // found if name is a substring of the section name
            return st;

    ABORT_R ("%s: bad argument for --show-headers - \"%s\", is not a recognized section type", 
           global_cmd, name);
}

uint32_t st_header_size (SectionType sec_type)
{
    ASSERT (sec_type >= SEC_NONE && sec_type < NUM_SEC_TYPES, "Error in st_header_size: sec_type=%u out of range [-1,%u]", sec_type, NUM_SEC_TYPES-1);

    return (sec_type == SEC_NONE) ? 0 : abouts[sec_type].header_size;
}

// called by PIZ I/O
const SectionListEntry *sections_get_first_section_of_type (SectionType st, bool soft_fail)
{
    ARRAY (const SectionListEntry, sl, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++)
        if (sl[i].section_type == st) return &sl[i];

    ASSERT (soft_fail, "Error in sections_get_first_section_of_type: Cannot find section_type=%s in z_file", st_name (st));

    return NULL;
}

const char *lt_name (LocalType lt)
{
    if (lt >= 0 && lt < NUM_LOCAL_TYPES) 
        return lt_desc[lt].name;
    else
        return "INVALID_LT";
}
