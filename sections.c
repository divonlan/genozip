// ------------------------------------------------------------------
//   sections.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "sections.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "strings.h"
#include "crypt.h"
#include "dict_id.h"
#include "zfile.h"
#include "piz.h"
#include "endianness.h"
#include "codec.h"
#include "bgzf.h"

typedef struct SectionEnt SectionEntModifiable;

static const struct {rom name; uint32_t header_size; } abouts[NUM_SEC_TYPES] = SECTIONTYPE_ABOUT;

const LocalTypeDesc lt_desc[NUM_LOCAL_TYPES] = LOCALTYPE_DESC;

typedef struct {
    Section first_sec, last_sec;
    Section next_vb_same_comp_sec; // linked list of VBs of the same comp. final VB is NULL.
} SectionsVbIndexEnt;

typedef struct {
    Section txt_header_sec, bgzf_sec, recon_plan_sec[2];
    uint32_t first_vb_i, num_vbs;
    SectionsVbIndexEnt *first_vb_index_ent, *last_vb_index_ent; // first and last VB in this component
} SectionsCompIndexEnt;

// ZIP only: create section list that goes into the genozip header, as we are creating the sections. returns offset
void sections_add_to_list (VBlockP vb, const SectionHeader *header)
{
    // case: we're re-creaating a section already on this list - nothing to do
    if (vb->section_list_buf.len && vb->z_data.len <= BLST (SectionEnt, vb->section_list_buf)->offset)
        return;

    SectionType st = header->section_type;
    DictId dict_id = DICT_ID_NONE;

    if (IS_DICTED_SEC(st)) {
        switch (st) {
            case SEC_DICT   : dict_id = ((SectionHeaderDictionary *)header)->dict_id; break;
            case SEC_B250   : dict_id = ((SectionHeaderCtx        *)header)->dict_id; break;
            case SEC_LOCAL  : dict_id = ((SectionHeaderCtx        *)header)->dict_id; break;
            case SEC_COUNTS : dict_id = ((SectionHeaderCounts     *)header)->dict_id; break;
            default         : ABORT0 ("new section type with dict_id?");
        }
        ASSERT0 (dict_id.num, "dict_id=0");
    }
    
    buf_alloc (vb, &vb->section_list_buf, 1, 50, SectionEnt, 2, "section_list_buf");
    
    BNXT (SectionEntModifiable, vb->section_list_buf) = (SectionEntModifiable) {
        .st        = header->section_type,
        .vblock_i  = (IS_VB_SEC(st) || IS_FRAG_SEC(st)) ? BGEN32 (header->vblock_i) : 0, // big endian in header - convert back to native
        .comp_i    = IS_COMP_SEC(st) ? vb->comp_i : COMP_NONE,
        .offset    = vb->z_data.len,  // this is a partial offset (within d) - we will correct it later
        .flags     = header->flags,
        .num_lines = st==SEC_VB_HEADER ? vb->lines.len32 : 0, 
        .size      = BGEN32 (header->compressed_offset) + BGEN32 (header->data_compressed_len)
    };

    if (IS_DICTED_SEC(st)) 
        BLST (SectionEntModifiable, vb->section_list_buf)->dict_id = dict_id; // note: this is in union with num_lines
}

// ZIP: remove section from list when it is deleted from z_data by zfile_remove_ctx_group_from_z_data
void sections_remove_from_list (VBlockP vb, uint64_t offset, uint64_t len)
{
    SectionEntModifiable *sec;
    for (sec=BLST (SectionEntModifiable, vb->section_list_buf); sec->offset > offset; sec--) 
        sec->offset -= len;

    ASSERT (sec->offset == offset, "cannot find section with offset=%"PRIu64" in vb=%s", offset, VB_NAME);

    buf_remove (vb->section_list_buf, SectionEnt, BNUM (vb->section_list_buf, sec), 1);
}

// Called by ZIP main thread. concatenates a vb or dictionary section list to the z_file section list - just before 
// writing those sections to the disk. we use the current disk position to update the offset
void sections_list_concat (VBlockP vb)
{
    ARRAY (SectionEntModifiable, vb_sec, vb->section_list_buf);

    // update the offset
    for (uint32_t i=0; i < vb_sec_len; i++)
        vb_sec[i].offset += z_file->disk_so_far;

    // copy all entries
    buf_add_buf (evb, &z_file->section_list_buf, &vb->section_list_buf, SectionEntModifiable, "z_file->section_list_buf");

    buf_free (vb->section_list_buf);
}

// section iterator. returns true if a section of this type was found.
bool sections_next_sec2 (Section *sl_ent,   // optional in/out. if NULL - search entire list
                         SectionType st1, SectionType st2) // check only next section, not entire remaining list
{
    Section sec = sl_ent ? *sl_ent : NULL; 
    bool found = false;

    ASSERT (!sec || BISVALID (z_file->section_list_buf, sec), "Invalid sec: st1=%s st2=%s", st_name (st1), st_name (st2));

    while (sec < BAFT (SectionEnt, z_file->section_list_buf) - 1) {

        sec = sec ? (sec + 1) : B1ST (SectionEnt, z_file->section_list_buf); 

        if (sec->st == st1 || sec->st == st2) {
            found = true;
            break;
        }
    }

    if (sl_ent) *sl_ent = sec;
    return found;
}

// section iterator. returns true if a section of this type was found.
bool sections_prev_sec2 (Section *sl_ent,   // optional in/out. if NULL - search entire list
                         SectionType st1, SectionType st2)
{
    Section sec = sl_ent ? *sl_ent : NULL; 

    ASSERT (!sec || (sec >= B1ST(SectionEnt, z_file->section_list_buf) && sec <= BLST(SectionEnt, z_file->section_list_buf)),
           "Invalid sec: st1=%s st2=%s", st_name (st1), st_name (st2));

    while (!sec || sec >= B1ST (SectionEnt, z_file->section_list_buf)) {

        sec = sec ? (sec - 1) : BLST (SectionEnt, z_file->section_list_buf); 

        if (sec->st == st1 || sec->st == st2) {
            if (sl_ent) *sl_ent = sec;
            return true;
        }
    }

    return false;
}

// section iterator. skips to the last consecutive section from after sec, that is type st1-4
Section sections_last_sec4 (Section sec, SectionType st1, SectionType st2, SectionType st3, SectionType st4)
{
    ASSERT (!sec || (sec >= B1ST(SectionEnt, z_file->section_list_buf) && sec <= BLST(SectionEnt, z_file->section_list_buf)),
           "Invalid sec: st1=%s st2=%s st3=%s st4=%s", st_name (st1), st_name (st2), st_name (st3), st_name (st4));

    while (sec < BAFT (SectionEnt, z_file->section_list_buf)) {

        sec = sec ? (sec + 1) : B1ST (SectionEnt, z_file->section_list_buf); 

        if (!((st1 != SEC_NONE && sec->st == st1) ||
              (st2 != SEC_NONE && sec->st == st2) ||
              (st3 != SEC_NONE && sec->st == st3) ||
              (st4 != SEC_NONE && sec->st == st4))) break;
    } 
    
    return sec - 1; // sec is either on the first section that is NOT st1 or st2, or we are after the section list - so we decrement it
}

// count how many sections we have of a certain type
uint32_t sections_count_sections_until (SectionType st, Section first_sec, SectionType until_encountering)
{
    ASSERT0 (!first_sec || (first_sec >= B1ST(SectionEnt, z_file->section_list_buf) && first_sec <= BLST(SectionEnt, z_file->section_list_buf)),
             "Invalid sec");

    ARRAY (SectionEnt, secs, z_file->section_list_buf);
    uint32_t count=0;
    uint32_t first_i = first_sec ? BNUM (z_file->section_list_buf, first_sec) : 0;

    for (uint32_t i=first_i; i < secs_len; i++) 
        if (secs[i].st == st)
            count++;
        else if (secs[i].st == until_encountering)
            break;

    return count;
}

static CompIType sections_calculate_num_comps (void)
{
    ARRAY (SectionEnt, secs, z_file->section_list_buf);
    CompIType max_comp_i = 0;

    for (uint32_t i=0; i < secs_len; i++) 
        if ((secs[i].st == SEC_VB_HEADER || secs[i].st == SEC_TXT_HEADER) && secs[i].comp_i > max_comp_i)
            max_comp_i = secs[i].comp_i;

    return max_comp_i + 1;
}

// -------------------
// VBs
// -------------------

// PIZ
// note: not all VBs and components in the z_file are in the section list or its index - some might have been removed
// from the section list by the recon_plan (eg in DVCF we remove the unneeded component for the rendition) 
static void sections_create_index (bool force)
{
    if (!force && z_file->comp_sections_index.len) return; // already created

    if (!z_file->num_vbs)
        z_file->num_vbs = sections_count_sections (SEC_VB_HEADER);

    CompIType num_comps = sections_calculate_num_comps(); 

    ARRAY_alloc (SectionsVbIndexEnt, vb_index, z_file->num_vbs + 1, true, z_file->vb_sections_index, evb, "z_file->vb_sections_index");
    ARRAY_alloc (SectionsCompIndexEnt, comp_index, num_comps, true, z_file->comp_sections_index, evb, "z_file->comp_sections_index");
    ARRAY (SectionEnt, sec, z_file->section_list_buf);
    
    // initialize (or reset if initialized before)
    z_file->first_dict_section = NULL;
    
    for (uint32_t sec_i=0; sec_i < sec_len; sec_i++) {
    
        CompIType comp_i = sec[sec_i].comp_i;
        
        switch (sec[sec_i].st) {
            
            case SEC_TXT_HEADER : 
                comp_index[comp_i].txt_header_sec = &sec[sec_i]; 
                break;
            
            case SEC_BGZF : 
                comp_index[comp_i].bgzf_sec = &sec[sec_i]; 
                break;

            case SEC_RECON_PLAN : {
                bool is_luft_plan = sec[sec_i].flags.recon_plan.luft;
                if (!comp_index[comp_i].recon_plan_sec[is_luft_plan]) // save first fragment
                    comp_index[comp_i].recon_plan_sec[is_luft_plan] = &sec[sec_i]; 
                break;
            }

            case SEC_VB_HEADER : {
                VBIType vb_i = sec[sec_i].vblock_i;
                
                ASSERT (vb_i>=1 && vb_i <= z_file->num_vbs, "Bad vb_i=%u, max is %u", vb_i, z_file->num_vbs);
                
                vb_index[vb_i] = (SectionsVbIndexEnt){
                    .first_sec = &sec[sec_i],
                    .last_sec  = sections_last_sec4 (&sec[sec_i], SEC_B250, SEC_LOCAL, SEC_NONE, SEC_NONE)
                };

                if (!comp_index[comp_i].first_vb_i) {
                    comp_index[comp_i].first_vb_index_ent = &vb_index[vb_i];
                    comp_index[comp_i].first_vb_i = vb_i;
                }

                else // link to previous VB in this component (not necessarily consecutive in the file)
                    comp_index[comp_i].last_vb_index_ent->next_vb_same_comp_sec = &sec[sec_i];
                    
                comp_index[comp_i].last_vb_index_ent = &vb_index[vb_i]; // last so far
                comp_index[comp_i].num_vbs++;
                sec_i = BNUM (z_file->section_list_buf, vb_index[vb_i].last_sec);
                break;
            }

            case SEC_DICT :
                if (!z_file->first_dict_section) z_file->first_dict_section = &sec[sec_i];
                break;

            default: break;
        }
    }
}
static inline void sections_create_index_if_needed(void) 
{ 
    if (!z_file->comp_sections_index.len) sections_create_index (false);
}

// ZIP / PIZ: returns the SEC_VB_HEADER section, or NULL if this vb_i doesn't exist
Section sections_vb_header (VBIType vb_i, bool soft_fail)
{
    sections_create_index (false);
    
    if (IS_PIZ) {
        ASSERTNOTEMPTY (z_file->vb_sections_index);
        
        if (vb_i <= z_file->num_vbs) {
            Section sec = B(SectionsVbIndexEnt, z_file->vb_sections_index, vb_i)->first_sec;
            if (!sec) goto fail;

            ASSBISVALID (z_file->section_list_buf, sec);
            return sec;
        }
        else 
            return NULL;
    }

    else {
        Section sec=NULL;
        uint32_t i=0; for (; i < z_file->section_list_buf.len32; i++) {
            sec = B(SectionEnt, z_file->section_list_buf, i);
            if (sec->vblock_i == vb_i) break; // found!
        }

        if (i >= z_file->section_list_buf.len) goto fail;
        return sec;
    }

fail:
    if (soft_fail) return NULL;
    ABORT_R ("cannot find SEC_VB_HEADER for vb_i=%u", vb_i);
}

Section sections_one_before (Section sec) 
{ 
    return (sec == B1ST (SectionEnt, z_file->section_list_buf)) ? NULL : sec-1; 
}

// PIZ
Section sections_vb_last (VBIType vb_i)
{
    ASSERT (vb_i <= z_file->num_vbs, "Bad vb_i=%u, max is %u", vb_i, z_file->num_vbs);
    return B(SectionsVbIndexEnt, z_file->vb_sections_index, vb_i)->last_sec;
}

// main thread: called by refhash_initialize - get details of the refhash ahead of loading it from the reference file 
void sections_get_refhash_details (uint32_t *num_layers, uint32_t *base_layer_bits) // optional outs
{
    ASSERT0 (flag.reading_reference || flag.show_ref_hash, "can only be called while reading reference");

    for (int i=z_file->section_list_buf.len-1; i >= 0; i--) { // search backwards as the refhash sections are near the end
        Section sec = B(SectionEnt, z_file->section_list_buf, i);
        if (sec->st == SEC_REF_HASH) {

            SectionHeaderRefHash header = zfile_read_section_header (evb, sec->offset, sec->vblock_i, SEC_REF_HASH).ref_hash;
            if (num_layers) *num_layers = header.num_layers;
            if (base_layer_bits) *base_layer_bits = header.layer_bits + header.layer_i; // layer_i=0 is the base layer, layer_i=1 has 1 bit less etc
            return;
        }
        else if (sec->st == SEC_REFERENCE)
            break; // we arrived at a SEC_REFERENCE - there won't be any more SEC_REF_HASH sections
    }

    ABORT ("can't find SEC_REF_HASH sections in %s", z_name);
}

//---------------------------------------------
// PIZ: constructing a updated section list
//---------------------------------------------

void sections_new_list_add_vb (BufferP new_list, VBIType vb_i)
{
    SectionsVbIndexEnt *vbent = B(SectionsVbIndexEnt, z_file->vb_sections_index, vb_i);
    uint32_t num_sections = vbent->last_sec - vbent->first_sec + 1;

    memcpy (BAFT (SectionEntModifiable, *new_list), vbent->first_sec, num_sections * sizeof (SectionEnt));
    new_list->len += num_sections;
}

void sections_new_list_add_txt_header (BufferP new_list, CompIType comp_i)
{
    SectionsCompIndexEnt *compent = B(SectionsCompIndexEnt, z_file->comp_sections_index, comp_i);
    
    BNXT (SectionEntModifiable, *new_list) = *compent->txt_header_sec;
}

// PIZ: If any of the components has a SEC_BGZF add it
void sections_new_list_add_bgzf (BufferP new_list)
{
    for (CompIType comp_i=0; comp_i < z_file->comp_sections_index.len; comp_i++) {
        SectionsCompIndexEnt *comp = B(SectionsCompIndexEnt, z_file->comp_sections_index, comp_i);

        if (comp->bgzf_sec) 
            BNXT (SectionEntModifiable, *new_list) = *comp->bgzf_sec;
    } 
}

void sections_new_list_add_global_sections (BufferP new_list)
{
    uint32_t num_sections = BAFT(SectionEnt, z_file->section_list_buf) - z_file->first_dict_section;

    memcpy (BAFT (SectionEntModifiable, *new_list), z_file->first_dict_section, num_sections * sizeof (SectionEnt));
    new_list->len += num_sections;
}

// replace current section list with the new list (if one exists)
void sections_commit_new_list (BufferP new_list)
{
    if (!new_list->len) return;

    ASSERT (new_list->len <= z_file->section_list_buf.len, "Expecting new_list->len=%"PRIu64" <= z_file->section_list_buf.len=%"PRIu64, 
            new_list->len, z_file->section_list_buf.len);

    buf_destroy (z_file->section_list_buf);
    buf_move (evb, &z_file->section_list_buf, evb, new_list); // note: in writer, tf_ent[]->last_sec points into new_list (which won't work if we ever change buf_move to buf_copy)

    sections_create_index (true);
}

//---------------
// misc functions
//---------------

// ZIP
void sections_list_memory_to_file_format (bool in_place) // in place, or to evb->scratch for testing codec
{
    if (!in_place) { // just generating for codec testing
        ASSERTNOTINUSE(evb->scratch);
        buf_alloc_exact (evb, evb->scratch, (MIN_(z_file->section_list_buf.len * sizeof (SectionEntFileFormat), CODEC_ASSIGN_SAMPLE_SIZE) / sizeof(SectionEntFileFormat)),
                         SectionEntFileFormat, "scratch");
    }

    BufferP out = in_place ? &z_file->section_list_buf : &evb->scratch;
    ARRAY (const SectionEntModifiable, mem_sec, z_file->section_list_buf); // memory format (entries are larger)
    ARRAY (SectionEntFileFormat, file_sec, *out); // file format
        
    for (uint32_t i=0; i < file_sec_len; i++) {
        file_sec[i] = (SectionEntFileFormat){
            .vblock_i    = BGEN32 (mem_sec[i].vblock_i),
            .comp_i      = mem_sec[i].comp_i != COMP_NONE ? mem_sec[i].comp_i : 0, // since v14
            .offset      = BGEN64 (mem_sec[i].offset),
            .st_specific = mem_sec[i].st_specific, 
            .st          = mem_sec[i].st,
            .flags       = mem_sec[i].flags    // since v12
        };
    
        if (mem_sec[i].st == SEC_VB_HEADER) 
            file_sec[i].num_lines = BGEN32 (file_sec[i].num_lines);
    }

    out->len *= sizeof (SectionEntFileFormat); // change to counting bytes
}

// PIZ: create in-memory format of the section list - copy from z_file section, BGEN, and add sizes
void sections_list_file_to_memory_format (SectionHeaderGenozipHeader *genozip_header)
{
    struct FlagsGenozipHeader f = genozip_header->h.flags.genozip_header;
    int V = genozip_header->genozip_version;
    DataType dt = BGEN16 (genozip_header->data_type);

    // For files v13 and earlier, we can only read them if they are a single component, or two paired FASTQs, or DVCF
    // Other bound files need to be decompressed with Genozip v13
    uint32_t v13_num_components = BGEN32 (genozip_header->v13_num_components);
    ASSINP (V >= 14 || 
            v13_num_components == 1 || // single file
            (dt == DT_VCF && f.has_gencomp) || // DVCF
            (dt == DT_FASTQ && (f.dts_paired || V <= 9) && v13_num_components == 2) || // Paired FASTQs (the dts_paired flag was introduced in V9.0.13)
            is_genols ||
            flags_is_genocat_global_area_only(), // only show meta-data (can't use flag.genocat_global_area_only bc flags are not set yet)
            "%s is comprised of %u files bound together. The bound file feature was discontinued in Genozip v14. To decompress this file, use Genozip v13",
            z_name, v13_num_components);
    
    // treat V8/V9 genozip files containing two FASTQ components as paired (the dts_paired flag was introduced in V9.0.13)    
    if (dt == DT_FASTQ && V <= 9 && v13_num_components == 2) 
        z_file->z_flags.dts_paired = true; 

    z_file->section_list_buf.len /= sizeof (SectionEntFileFormat); // fix len

    buf_alloc (evb, &z_file->section_list_buf, 0, z_file->section_list_buf.len, SectionEnt, 0, NULL); // extend

    ARRAY (SectionEntModifiable, mem_sec,  z_file->section_list_buf); // memory format (entries are larger)
    ARRAY (const SectionEntFileFormat, file_sec, z_file->section_list_buf); // file format

    // note: we work backwards as mem_sec items are larger than file_sec
    for (int i=file_sec_len-1; i >= 0; i--) {
        SectionEntFileFormat sec = file_sec[i]; // make a copy as assignment is overlapping
        
        mem_sec[i] = (SectionEntModifiable){
            .offset      = BGEN64 (sec.offset),
            .st_specific = sec.st_specific,  
            .vblock_i    = BGEN32 (sec.vblock_i),
            .st          = sec.st,
            .comp_i      = (V >= 14 && IS_COMP_SEC(sec.st)) ? sec.comp_i : COMP_NONE, // note: in file format, COMP_NONE sections have 0, as comp_i is just 2 bit
            .flags       = (V >= 12) ? sec.flags  : (SectionFlags){} // flags were introduced in v12
        };

        if (mem_sec[i].st == SEC_VB_HEADER) 
            mem_sec[i].num_lines = BGEN32 (mem_sec[i].num_lines);        

        if (i < file_sec_len-1)
            mem_sec[i].size = mem_sec[i+1].offset - mem_sec[i].offset;
    }

    // comp_i was introduced V14, for V<=13, get comp_i by the component's relative position in the file. 
    if (V <= 13) {
        CompIType comp_i_by_consecutive = COMP_NONE; 

        for (int i=0; i < mem_sec_len; i++) {
            SectionEntModifiable sec = mem_sec[i]; // make a copy as assignment is overlapping
            
            if (sec.st == SEC_TXT_HEADER) comp_i_by_consecutive++; // note: we only read bound V<=13 files if the are a single FQ pair, or a DVCF

            if (IS_COMP_SEC(sec.st))
                mem_sec[i].comp_i = comp_i_by_consecutive; 
        }
    }

    mem_sec[file_sec_len-1].size = BGEN32 (genozip_header->h.data_compressed_len) + 
                                   BGEN32 (genozip_header->h.compressed_offset) + 
                                   sizeof (SectionFooterGenozipHeader);

    sections_create_index (true);

    // Up to V11, top_level_repeats existed (called num_lines) but was not populated
    // in V12-13 num_lines was transmitted through the VbHeader.top_level_repeats
    // Since V14, num_lines is transmitted through SectionEntFileFormat
    if (V >= 12 && V <= 13) 
        for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {
            SectionEntModifiable *sec = (SectionEntModifiable *)sections_vb_header (vb_i, false);
            SectionHeaderUnion header = zfile_read_section_header (evb, sec->offset, vb_i, SEC_VB_HEADER);
            sec->num_lines = BGEN32 (header.vb_header.v13_top_level_repeats);
        }

    // copy vb_i=1 sections to section_list_vb1, because they are needed for setting the default context
    // flags, and vb_i=1 may be eliminated by the recon_plan (for example, FASTQ with --R2)
    if (z_file->num_vbs >= 1) { // note: some files have no VBs, eg DT_REF, or a header-only file
        Section first_sec_vb1 = sections_vb_header (1, false);

        uint32_t first_sec_i = BNUM (z_file->section_list_buf, first_sec_vb1);
        uint32_t num_secs = BNUM (z_file->section_list_buf, sections_vb_last (1)) - first_sec_i + 1;

        buf_copy (evb, &z_file->section_list_vb1, &z_file->section_list_buf, SectionEnt, first_sec_i, num_secs, "z_file->section_list_vb1");
    }
}

rom st_name (SectionType sec_type)
{
    ASSERT (sec_type >= SEC_NONE && sec_type < NUM_SEC_TYPES, "sec_type=%u out of range [-1,%u]", sec_type, NUM_SEC_TYPES-1);

    return (sec_type == SEC_NONE) ? "SEC_NONE" : type_name (sec_type, &abouts[sec_type].name , ARRAY_LEN(abouts));
}

rom comp_name (CompIType comp_i)
{
    static int max_comps_by_dt[NUM_DATATYPES] = { [DT_VCF]=3, [DT_BCF]=3, [DT_SAM]=3, [DT_BAM]=3, [DT_FASTQ]=2 };
    static rom comp_names[NUM_DATATYPES][3] = { [DT_VCF] = VCF_COMP_NAMES, [DT_BCF] = VCF_COMP_NAMES,
                                                [DT_SAM] = SAM_COMP_NAMES, [DT_BAM] = SAM_COMP_NAMES,
                                                [DT_FASTQ] = FASTQ_COMP_NAMES };
    if (!z_file) return "EMEM"; // can happen if another thread is busy exiting the process
    
    DataType dt = z_file->data_type;

    if (!max_comps_by_dt[dt] && comp_i==0)
        return "MAIN";

    else if (comp_i==COMP_NONE)
        return "NONE";

    else if (comp_i >= max_comps_by_dt[dt]) {
        static char invalid_comp_number[24];
        sprintf (invalid_comp_number, "invalid-comp_i-%d", comp_i); // not thread safe, but should never happen
        return invalid_comp_number;
    }

    else
        return comp_names[dt][comp_i];
}

rom comp_name_ex (CompIType comp_i, SectionType st)
{
    if (IS_COMP_SEC (st))     return comp_name (comp_i);
    if (st == SEC_REFERENCE)  return "REFR";
    if (st == SEC_REF_IS_SET) return "REFR";
    if (st == SEC_REF_HASH)   return "REFH";
    if (st == SEC_DICT)       return "DICT";
    else                      return "----";
}

VbNameStr vb_name (VBlockP vb)
{
    VbNameStr s;
    if (vb && vb->vblock_i)
        sprintf (s.s, "%s/%u", comp_name (vb->comp_i), vb->vblock_i);
    else if (vb && segconf.running)
        strcpy (s.s, "SEGCONF/0");
    else
        strcpy (s.s, "NONCOMPUTE");

    return s;
}

LineNameStr line_name (VBlockP vb)
{
    LineNameStr s;
    sprintf (s.s, "%s/%u%s", vb_name(vb).s, vb->line_i, vb->preprocessing ? "(preproc)" : "");
    return s;
}

// called to parse the optional argument to --show-headers. we accept eg "REFERENCE" or "SEC_REFERENCE" or even "ref"
SectionType sections_st_by_name (char *name)
{
    if (!name) return -2; // all sections

    str_toupper (name, name); // name is case insesitive (compare uppercase)

    SectionType candidate=-1;
    unsigned candidate_len=100000;

    // choose the shortest section for which name is a case-insensitive substring (eg "dict" will result in SEC_DICT not SEC_DICT_ID_ALIASES)
    for (SectionType st=0; st < NUM_SEC_TYPES; st++)
        if (strstr (abouts[st].name, name)) { // found if name is a substring of the section name
            unsigned len = strlen (abouts[st].name);
            if (len < candidate_len) { // best candidate so far
                candidate     = st;
                candidate_len = len;
            }
        }

    ASSINP (candidate >= 0, "bad argument - \"%s\", is not a recognized section type (or a case-insensitive sub-string of a section type)", name);

    return candidate;
}

uint32_t st_header_size (SectionType sec_type)
{
    ASSERT (sec_type >= SEC_NONE && sec_type < NUM_SEC_TYPES, "sec_type=%u out of range [-1,%u]", sec_type, NUM_SEC_TYPES-1);

    return sec_type == SEC_NONE ? 0 
         : sec_type == SEC_VB_HEADER      && command != ZIP && !VER(12) ? sizeof (SectionHeaderVbHeader) - 4*sizeof(uint32_t) // in v8-11, SectionHeaderVbHeader was shorter by 4 32b words
         : sec_type == SEC_VB_HEADER      && command != ZIP && !VER(14) ? sizeof (SectionHeaderVbHeader) - 1*sizeof(uint32_t) // in v12-13, SectionHeaderVbHeader was shorter by 1 32b word
         : sec_type == SEC_TXT_HEADER     && command != ZIP && !VER(12) ? sizeof (SectionHeaderTxtHeader) -  sizeof(uint64_t) // in v8-11, SectionHeaderTxtHeader was shorter
         : sec_type == SEC_GENOZIP_HEADER && command != ZIP && !VER(12) ? sizeof (SectionHeaderGenozipHeader) -  REF_FILENAME_LEN - sizeof(Digest) // in v8-11, SectionHeaderTxtHeader was shorter
         : abouts[sec_type].header_size;
}

Section section_next (Section sec)
{
    if (!z_file->section_list_buf.len) return NULL;
    if (sec == NULL) return B1ST (SectionEnt, z_file->section_list_buf);
    if (sec < BLST (SectionEnt, z_file->section_list_buf)) return sec+1;
    return NULL;
}

// called by PIZ main thread
Section sections_first_sec (SectionType st, bool soft_fail)
{
    ARRAY (SectionEnt, sec, z_file->section_list_buf);

    for (unsigned i=0; i < sec_len; i++)
        if (sec[i].st == st) return &sec[i];

    ASSERT (soft_fail, "Cannot find section_type=%s in z_file", st_name (st));

    return NULL;
}

// called by ZIP, PIZ main thread
Section sections_last_sec (SectionType st, bool soft_fail)
{
    ARRAY (SectionEnt, sec, z_file->section_list_buf);

    for (int i=sec_len-1; i >= 0; i--)
        if (sec[i].st == st) return &sec[i];

    ASSERT (soft_fail, "Cannot find section_type=%s in z_file", st_name (st));

    return NULL;
}

CompIType sections_get_num_comps (void)
{
    sections_create_index_if_needed();

    return z_file->comp_sections_index.len;
}

static const SectionsCompIndexEnt *sections_get_comp_index_ent (CompIType comp_i)
{
    sections_create_index_if_needed();

    ASSERTNOTEMPTY(z_file->comp_sections_index);
    ASSERT (comp_i < z_file->comp_sections_index.len, "comp_i=%u out of range [0,%u]", comp_i, (int)z_file->comp_sections_index.len-1);

    const SectionsCompIndexEnt *comp_index_ent = B(SectionsCompIndexEnt, z_file->comp_sections_index, comp_i);
    return comp_index_ent;
}

VBIType sections_get_num_vbs (CompIType comp_i) 
{ 
    sections_create_index_if_needed();

    if (comp_i >= z_file->comp_sections_index.len) return 0; // this component doesn't exist (can happen eg when checking for gencomp that doesn't exist in this file)

    return sections_get_comp_index_ent (comp_i)->num_vbs;
}

VBIType sections_get_first_vb_i (CompIType comp_i) 
{ 
    return sections_get_comp_index_ent (comp_i)->first_vb_i;
}

Section sections_get_comp_txt_header_sec (CompIType comp_i)
{
    return sections_get_comp_index_ent (comp_i)->txt_header_sec;
}

Section sections_get_comp_recon_plan_sec (CompIType comp_i, bool is_luft_plan)
{
    return sections_get_comp_index_ent (comp_i)->recon_plan_sec[is_luft_plan];
}

Section sections_get_next_vb_of_comp_sec (CompIType comp_i, Section *vb_sec)
{
    const SectionsCompIndexEnt *comp_index_ent = sections_get_comp_index_ent (comp_i);

    if (! *vb_sec)
        *vb_sec = comp_index_ent->first_vb_index_ent->first_sec;
    
    else 
        *vb_sec = B(SectionsVbIndexEnt, z_file->vb_sections_index, (*vb_sec)->vblock_i)->next_vb_same_comp_sec;

    return *vb_sec;
}

rom lt_name (LocalType lt)
{
    if (lt >= 0 && lt < NUM_LOCAL_TYPES) 
        return lt_desc[lt].name;
    else
        return "INVALID_LT";
}

rom store_type_name (StoreType store)
{
    switch (store) {
        case STORE_NONE  : return "NONE";
        case STORE_INT   : return "INT";
        case STORE_FLOAT : return "FLOAT";
        case STORE_INDEX : return "INDEX";
        default          : return "InvalidStoreType";
    }
}

typedef struct { char s[128]; } FlagStr;
static FlagStr sections_dis_flags (SectionFlags f, SectionType st, DataType dt)
{
    static rom dts[NUM_DATATYPES] = { [DT_FASTQ]="dts_paired", [DT_SAM]="dts_ref_internal", [DT_CHAIN]="dts_mismatch" };

    FlagStr str = {};

    switch (st) {
        case SEC_GENOZIP_HEADER:
            sprintf (str.s, "%s=%u aligner=%u txt_is_bin=%u bgzf=%u adler=%u has_gencomp=%u has_taxid=%u",
                     dts[dt] ? dts[dt] : "dt_specitic", f.genozip_header.dt_specific, f.genozip_header.aligner, 
                     f.genozip_header.txt_is_bin, f.genozip_header.bgzf, f.genozip_header.adler, f.genozip_header.has_gencomp,
                     f.genozip_header.has_taxid);
            break;

        case SEC_TXT_HEADER:
            sprintf (str.s, "is_txt_luft=%u", f.txt_header.is_txt_luft);
            break;

        case SEC_VB_HEADER:
            switch (dt) {
                case DT_VCF: 
                    sprintf (str.s, "coords=%s null_DP=%u", vcf_coords_name (f.vb_header.vcf.coords), f.vb_header.vcf.use_null_DP_method);
                    break;
                case DT_SAM: case DT_BAM:
                    if (VER(13) && !VER(14)) sprintf (str.s, "sorted=%u collated=%u", f.vb_header.sam.v13_is_sorted, f.vb_header.sam.v13_is_collated);
                    break;
                default:
                    str.s[0] = 0;
            }
            break;

        case SEC_BGZF:
            sprintf (str.s, "library=%s level=%u has_eof=%u", bgzf_library_name(f.bgzf.library), f.bgzf.level, f.bgzf.has_eof_block); 
            break;

        case SEC_LOCAL:
            sprintf (str.s, "store=%-5s per_ln=%u delta=%u paired=%u spl_custom=%u specific=%u",
                     store_type_name(f.ctx.store), f.ctx.store_per_line, f.ctx.store_delta, f.ctx.paired, f.ctx.spl_custom, f.ctx.ctx_specific_flag); // note: we don't print ctx_specific as its not currently used
            break;

        case SEC_B250:
            sprintf (str.s, "store=%-5s per_ln=%u delta=%u paired=%u spl_custom=%u specific=%u same=%u",
                     store_type_name(f.ctx.store), f.ctx.store_per_line, f.ctx.store_delta, f.ctx.paired, f.ctx.spl_custom, f.ctx.ctx_specific_flag, f.ctx.all_the_same); 
            break;
 
        case SEC_RANDOM_ACCESS:
        case SEC_RECON_PLAN:
            sprintf (str.s, "is_luft=%u", f.recon_plan.luft);
            break;

        default: 
            str.s[0] = 0;
    }

    return str;
}

void sections_show_header (const SectionHeader *header, VBlockP vb /* optional if output to buffer */, uint64_t offset, char rw)
{
    if (flag_loading_auxiliary && !flag.debug_read_ctxs) return; // don't show headers of an auxiliary file in --show-header, but show in --debug-read-ctx
    
    if (  flag.show_headers   != -1 &&                   // we don't need to show all sections
          flag.show_headers-1 != header->section_type && // we don't need to show this section
          !flag.debug_read_ctxs)
        return;

    bool is_dict_offset = (header->section_type == SEC_DICT && rw == 'W'); // at the point calling this function in zip, SEC_DICT offsets are not finalized yet and are relative to the beginning of the dictionary area in the genozip file
    bool v12 = (IS_ZIP || VER(12));
    bool v14 = (IS_ZIP || VER(14));

    char str[2048];
    #define PRINT { if (vb) buf_append_string (vb, &vb->show_headers_buf, str); else iprintf ("%s", str); } 
    #define SEC_TAB "            ++  "

    sprintf (str, "%c %s%-*"PRIu64" %-19s %-4.4s %-4.4s vb=%-3u z_off=%-6u txt_len=%-7u z_len=%-7u enc_len=%-7u%s ",
             rw, 
             is_dict_offset ? "~" : "", 9-is_dict_offset, offset, 
             st_name(header->section_type), 
             codec_name (header->codec), codec_name (header->sub_codec),
             BGEN32 (header->vblock_i), 
             BGEN32 (header->compressed_offset), 
             BGEN32 (header->data_uncompressed_len), 
             BGEN32 (header->data_compressed_len), 
             BGEN32 (header->data_encrypted_len), 
             BGEN32 (header->magic) != GENOZIP_MAGIC ? " BAD-MAGIC" : ""); // usually, we don't want the magic to take up line real estatee
    PRINT;

    SectionFlags f = header->flags;
    SectionType st = header->section_type;
    DataType dt = z_file->data_type;

    switch (st) {

    case SEC_GENOZIP_HEADER: {
        SectionHeaderGenozipHeader *h = (SectionHeaderGenozipHeader *)header;
        z_file->z_flags.adler = h->h.flags.genozip_header.adler; // needed by digest_display_ex
        char dt_specific[REF_FILENAME_LEN + 200] = "";

        if (dt == DT_CHAIN)
            sprintf (dt_specific, SEC_TAB "prim=\"%.*s\" md5=%s\n", REF_FILENAME_LEN, h->chain.prim_filename, 
                     digest_display_ex (h->chain.prim_file_md5, DD_MD5).s);

        else if ((dt == DT_SAM || dt == DT_BAM) && v14)
            sprintf (dt_specific, SEC_TAB "segconf=(sorted=%u,collated=%u,seq_len=%u,seq_len_to_cm=%u,ms_type=%u,has_MD_or_NM=%u,bisulfite=%u,MD_NM_by_unconverted=%u,predict_meth=%u,is_paired=%u,sag_type=%s,sag_has_AS=%u,pysam_qual=%u,cellranger=%u,SA_HtoS=%u,seq_len_dict_id=%s)\n", 
                     h->sam.segconf_is_sorted, h->sam.segconf_is_collated, BGEN32 (h->sam.segconf_seq_len), h->sam.segconf_seq_len_cm, h->sam.segconf_ms_type, h->sam.segconf_has_MD_or_NM, 
                     h->sam.segconf_bisulfite, h->sam.segconf_MD_NM_by_un, h->sam.segconf_predict_meth, 
                     h->sam.segconf_is_paired, sag_type_name(h->sam.segconf_sag_type), h->sam.segconf_sag_has_AS, 
                     h->sam.segconf_pysam_qual, h->sam.segconf_cellranger, h->sam.segconf_SA_HtoS, dis_dict_id(h->sam.segconf_seq_len_dict_id).s);

        else if (dt == DT_REF)
            sprintf (dt_specific, SEC_TAB "fast_md5=%s\n", digest_display (h->REF_fasta_md5).s);

        else if (dt == DT_FASTQ && h->genozip_version <= 13) 
            sprintf (dt_specific, SEC_TAB "bound_digest=%s segconf_seq_len_dict_id=%s\n", 
                     digest_display (h->FASTQ_v13_digest_bound).s, dis_dict_id(h->fastq.segconf_seq_len_dict_id).s);

        sprintf (str, "\n"SEC_TAB "ver=%u enc=%s dt=%s usize=%"PRIu64" lines=%"PRIu64" secs=%u comps=%u vb_size=%u\n" 
                      SEC_TAB "%s ref=\"%.*s\" md5ref=%s\n"
                      "%s" // dt_specific, if there is any
                      SEC_TAB "created=\"%.*s\"\n",
                 h->genozip_version, encryption_name (h->encryption_type), dt_name (dt), 
                 BGEN64 (h->recon_size_prim), BGEN64 (h->num_lines_bound), BGEN32 (h->num_sections), BGEN32 (h->num_components),
                 BGEN16(h->vb_size), sections_dis_flags (f, st, dt).s,
                 REF_FILENAME_LEN, h->ref_filename, digest_display_ex (h->ref_file_md5, DD_MD5).s,
                 dt_specific, 
                 FILE_METADATA_LEN, h->created);
        break;
    }

    case SEC_TXT_HEADER: {
        SectionHeaderTxtHeader *h = (SectionHeaderTxtHeader *)header;
        sprintf (str, "\n"SEC_TAB "txt_data_size=%"PRIu64" txt_header_size=%"PRIu64" lines=%"PRIu64" max_lines_per_vb=%u digest=%s digest_header=%s\n" 
                      SEC_TAB "txt_codec=%s (args=0x%02X.%02X.%02X) %s txt_filename=\"%.*s\"\n",
                 BGEN64 (h->txt_data_size), v12 ? BGEN64 (h->txt_header_size) : 0, BGEN64 (h->txt_num_lines), BGEN32 (h->max_lines_per_vb), 
                 digest_display (h->digest).s, digest_display (h->digest_header).s, 
                 codec_name (h->codec), h->codec_info[0], h->codec_info[1], h->codec_info[2], 
                 sections_dis_flags (f, st, dt).s, TXT_FILENAME_LEN, h->txt_filename);
        break;
    }

    case SEC_VB_HEADER: {
        SectionHeaderVbHeader *h = (SectionHeaderVbHeader *)header;
        if (Z_DT(DT_VCF)) 
            sprintf (str, 
                    "\n"SEC_TAB "recon_size=(PRIM:%u, LUFT:%u) longest_line=%u z_data_bytes=%u digest=%s %s\n",
                    BGEN32 (h->recon_size_prim), v12 ? BGEN32 (h->dvcf_recon_size_luft) : 0, 
                    BGEN32 (h->longest_line_len), 
                    BGEN32 (h->z_data_bytes), digest_display (h->digest).s, sections_dis_flags (f, st, dt).s);
        else if (Z_DT(DT_SAM))
            sprintf (str, 
                    "\n"SEC_TAB "recon_size=%u longest_line=%u z_data_bytes=%u digest=%s prim=(seq=%u comp_qual=%u qname=%u num_alns=%u first_grp_i=%u %s=%u) %s\n", 
                    BGEN32 (h->recon_size_prim),  BGEN32 (h->longest_line_len), 
                    BGEN32 (h->z_data_bytes), digest_display (h->digest).s, 
                    v14 ? BGEN32 (h->sam_prim_seq_len)          : 0,
                    v14 ? BGEN32 (h->sam_prim_comp_qual_len)    : 0,
                    v14 ? BGEN32 (h->sam_prim_qname_len)        : 0,
                    v14 ? BGEN32 (h->sam_prim_num_sag_alns)     : 0,
                    v14 ? BGEN32 (h->sam_prim_first_grp_i)      : 0,
                    IS_SAG_SA?"comp_cigars" : IS_SAG_SOLO?"solo_data" : "unused",
                    v14 ? BGEN32 (h->sam_prim_comp_cigars_len)  : 0,
                    sections_dis_flags (f, st, dt).s);
        else
            sprintf (str, 
                    "\n"SEC_TAB "recon_size=%u longest_line=%u z_data_bytes=%u digest=%s %s\n",
                    BGEN32 (h->recon_size_prim),BGEN32 (h->longest_line_len), 
                    BGEN32 (h->z_data_bytes), digest_display (h->digest).s, sections_dis_flags (f, st, dt).s);

        break;
    }

    case SEC_REFERENCE:
    case SEC_REF_IS_SET: {
        SectionHeaderReference *h = (SectionHeaderReference *)header;
        sprintf (str, "pos=%-9"PRIu64" gpos=%-9"PRIu64" num_bases=%-6u chrom_word_index=%-4d\n",
                 BGEN64 (h->pos), BGEN64 (h->gpos), BGEN32 (h->num_bases), BGEN32 (h->chrom_word_index)); 
        break;
    }
    
    case SEC_REF_HASH: {
        SectionHeaderRefHash *h = (SectionHeaderRefHash *)header;
        sprintf (str, "num_layers=%u layer_i=%u layer_bits=%u start_in_layer=%u\n",
                 h->num_layers, h->layer_i, h->layer_bits, BGEN32 (h->start_in_layer)); 
        break;
    }
    
    case SEC_REF_CONTIGS: {
        sprintf (str, "sequential_ref_index=%u\n", header->flags.ref_contigs.sequential_ref_index);
        break;
    }
    
    case SEC_RECON_PLAN: {
        SectionHeaderReconPlan *h = (SectionHeaderReconPlan *)header;
        sprintf (str, "conc_writing_vbs=%u %s\n", BGEN32 (h->conc_writing_vbs), sections_dis_flags (f, st, dt).s); 
        break;
    }
        
    case SEC_BGZF:
    case SEC_RANDOM_ACCESS: {
        sprintf (str, SEC_TAB "%s\n", sections_dis_flags (f, st, dt).s); 
        break;
    }
    
    case SEC_B250: {
        SectionHeaderCtx *h = (SectionHeaderCtx *)header;
        sprintf (str, "%s/%-8s\tltype=%s b250_size=%u param=%u %s\n",
                 dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, lt_name (h->ltype), 
                 h->b250_size==B250_BYTES_1?1 : h->b250_size==B250_BYTES_2?2 : h->b250_size==B250_BYTES_3?3 : 4,
                 h->param, sections_dis_flags (f, st, dt).s);
        break;
    }

    case SEC_LOCAL: {
        SectionHeaderCtx *h = (SectionHeaderCtx *)header;
        sprintf (str, "%s/%-8s\tltype=%s param=%u %s\n",
                 dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, lt_name (h->ltype), h->param, sections_dis_flags (f, st, dt).s);
        break;
    }

    case SEC_DICT: {
        SectionHeaderDictionary *h = (SectionHeaderDictionary *)header;
        sprintf (str, "%s/%-8s\tnum_snips=%u\n", dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, BGEN32 (h->num_snips)); 
        break;
    }

    case SEC_COUNTS: {
        SectionHeaderCounts *h = (SectionHeaderCounts *)header;
        sprintf (str, "  %s/%-8s\t\n", dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s); 
        break;
    }

    default: 
        str[0] = '\n'; str[1] = 0; 
    }

    PRINT;
}

// called from main()
void genocat_show_headers (rom z_filename)
{
    z_file = file_open (z_filename, READ, Z_FILE, DT_NONE);    

    SectionHeaderGenozipHeader header;
    TEMP_FLAG (show_headers, 0);
    zfile_read_genozip_header (&header);
    RESTORE_FLAG (show_headers);

    ARRAY (SectionEnt, sec, z_file->section_list_buf);
    for (uint32_t i=0; i < sec_len; i++) {
        header = zfile_read_section_header (evb, sec[i].offset, sec[i].vblock_i, sec[i].st).genozip_header; // we assign the largest of the SectionHeader* types
        sections_show_header ((SectionHeader *)&header, NULL, sec[i].offset, 'R');
    }        

    exit_ok();
}

void sections_show_section_list (DataType dt) // optional - take data from z_data
{    
    for (Section s=B1ST (SectionEnt, z_file->section_list_buf); s < BAFT (SectionEnt, z_file->section_list_buf); s++)
        if (s->st == SEC_B250 || s->st == SEC_LOCAL || s->st == SEC_DICT)
            iprintf ("%5u %-20.20s %s%s%-8.8s\tvb=%s/%-4u offset=%-8"PRIu64"  size=%-6u  %s\n", 
                     BNUM(z_file->section_list_buf, s), st_name(s->st), 
                     s->dict_id.num ? dtype_name_z(s->dict_id) :"     ", 
                     s->dict_id.num ? "/" : "", 
                     s->dict_id.num ? dis_dict_id (s->dict_id).s : "", 
                     comp_name_ex (s->comp_i, s->st), s->vblock_i, s->offset, s->size, 
                     sections_dis_flags (s->flags, s->st, dt).s);
        
        else if (s->st == SEC_VB_HEADER)
            iprintf ("%5u %-20.20s\t\t\tvb=%s/%-4u offset=%-8"PRIu64"  size=%-6u  num_lines=%-8u%s\n", 
                     BNUM(z_file->section_list_buf, s), st_name(s->st), comp_name (s->comp_i), 
                     s->vblock_i, s->offset, s->size, s->num_lines, sections_dis_flags (s->flags, s->st, dt).s);

        else if (IS_FRAG_SEC(s->st) || s->st == SEC_BGZF)
            iprintf ("%5u %-20.20s\t\t\tvb=%s/%-4u offset=%-8"PRIu64"  size=%-6u  %s\n", 
                     BNUM(z_file->section_list_buf, s), st_name(s->st), 
                     comp_name_ex (s->comp_i, s->st), s->vblock_i, s->offset, s->size, sections_dis_flags (s->flags, s->st, dt).s);

        else
            iprintf ("%5u %-20.20s\t\t\t             offset=%-8"PRIu64"  size=%-6u  %s\n", 
                     BNUM(z_file->section_list_buf, s), st_name(s->st),
                     s->offset, s->size, sections_dis_flags (s->flags, s->st, dt).s);
}

void sections_show_gheader (const SectionHeaderGenozipHeader *header)
{
    if (flag_loading_auxiliary) return; // don't show gheaders of an auxiliary file
    
    DataType dt = BGEN16 (header->data_type);
    ASSERT (dt < NUM_DATATYPES, "Invalid data_type=%u", dt);
    
    if (header) {
        iprintf ("Contents of the SEC_GENOZIP_HEADER section (output of --show-gheader) of %s:\n", z_name);
        iprintf ("  genozip_version: %u\n",         header->genozip_version);
        iprintf ("  data_type: %s\n",               dt_name (dt));
        iprintf ("  encryption_type: %s\n",         encryption_name (header->encryption_type)); 
        iprintf ("  recon_size_prim: %s\n",         str_int_commas (BGEN64 (header->recon_size_prim)).s);
        iprintf ("  num_lines_bound: %"PRIu64"\n",  BGEN64 (header->num_lines_bound));
        iprintf ("  num_sections: %u\n",            z_file->section_list_buf.len32);
        iprintf ("  num_components: %u\n",          header->num_components);
        if (dt == DT_REF)
            iprintf ("  REF_fasta_md5: %s\n",           digest_display (header->REF_fasta_md5).s);
        iprintf ("  created: %*s\n",                -FILE_METADATA_LEN, header->created);
        iprintf ("  license_hash: %s\n",            digest_display (header->license_hash).s);
        if (header->ref_filename[0]) {
            iprintf ("  reference filename: %s\n",  header->ref_filename);
            iprintf ("  reference file hash: %s\n", digest_display (header->ref_file_md5).s);
        }
        iprintf ("  flags: %s\n",                   sections_dis_flags (header->h.flags, SEC_GENOZIP_HEADER, dt).s);

        switch (dt) {
            case DT_CHAIN:
                iprintf ("  primary-coordinates reference filename: %s\n", header->chain.prim_filename);
                iprintf ("  primary-coordinates reference file hash: %s\n", digest_display (header->chain.prim_file_md5).s);
                break;

            default: break;
        }
    }

    iprint0 ("  sections:\n");

    sections_show_section_list (dt);
}
