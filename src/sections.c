// ------------------------------------------------------------------
//   sections.c
//   Copyright (C) 2020-2024 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
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
#include "threads.h"

typedef struct SectionEnt SectionEntModifiable;

static const struct {rom name; uint32_t header_size; } abouts[NUM_SEC_TYPES] = SECTIONTYPE_ABOUT;

const LocalTypeDesc lt_desc[NUM_LOCAL_TYPES] = LOCALTYPE_DESC;

typedef struct SectionsVbIndexEnt {
    int32_t vb_header_sec_i, last_sec_i; // -1 if none
    struct SectionsVbIndexEnt *next_vb_of_comp; // linked list of VBs of the same comp. list is terminated with -1. VBs are in the order they appear in section list, not necessarily consecutive vb_i's.
} SectionsVbIndexEnt;

typedef struct {
    int32_t txt_header_sec_i, bgzf_sec_i, recon_plan_sec_i[2];
    VBIType first_vb_i, num_vbs;
    SectionsVbIndexEnt *first_vb_of_comp, *last_vb_of_comp; // first and last VB in this component - in the order the appear in the section list (not necessarily consecutive vb_i)
} SectionsCompIndexEnt;

static Section Bsec (int32_t sec_i)
{
    if (sec_i == -1)
        return NULL;

    else if (sec_i < z_file->section_list_buf.len32)
        return B(SectionEnt, z_file->section_list_buf, sec_i);

    else
        ABORT ("sec_i=%d out of range [0,%u]", sec_i, z_file->vb_sections_index.len32-1);
}

static const SectionsVbIndexEnt *Bvbindex (VBIType vb_i)
{
    ASSERT (vb_i >= 1 && vb_i < z_file->vb_sections_index.len32, "vb_i=%u out of range [1,%d] (this happens if flag.pair is wrong or other reasons)", vb_i, (int)z_file->vb_sections_index.len32-1);

    const SectionsVbIndexEnt *vb_index_ent = B(SectionsVbIndexEnt, z_file->vb_sections_index, vb_i);
    return vb_index_ent;
}

static const SectionsCompIndexEnt *Bcompindex (CompIType comp_i)
{
    ASSERT (comp_i < z_file->comp_sections_index.len32, "comp_i=%u out of range [0,%u]", comp_i, z_file->comp_sections_index.len32-1);

    const SectionsCompIndexEnt *comp_index_ent = B(SectionsCompIndexEnt, z_file->comp_sections_index, comp_i);
    return comp_index_ent;
}

DictId sections_get_dict_id (ConstSectionHeaderP header)
{
    if (!header) return DICT_ID_NONE;
    
    switch (header->section_type) {
        case SEC_DICT     : return ((SectionHeaderDictionaryP)header)->dict_id; break;
        case SEC_B250     : return ((SectionHeaderCtxP       )header)->dict_id; break;
        case SEC_LOCAL    : return ((SectionHeaderCtxP       )header)->dict_id; break;
        case SEC_COUNTS   : return ((SectionHeaderCountsP    )header)->dict_id; break;
        case SEC_SUBDICTS : return ((SectionHeaderSubDictsP  )header)->dict_id; break;
        default           : return DICT_ID_NONE;
    }
}

// ZIP only: create section list that goes into the genozip header, as we are creating the sections. returns offset
void sections_add_to_list (VBlockP vb, ConstSectionHeaderP header)
{
    // case: we're re-creaating a section already on this list - nothing to do
    if (vb->section_list_buf.len && vb->z_data.len <= BLST (SectionEnt, vb->section_list_buf)->offset)
        return;

    SectionType st = header->section_type;
    ASSERT (st >= SEC_NONE && st < NUM_SEC_TYPES, "sec_type=%u out of range [-1,%u]", st, NUM_SEC_TYPES-1);

    DictId dict_id = sections_get_dict_id (header);
    ASSERT0 (!IS_DICTED_SEC(st) || dict_id.num, "dict_id=0");
    
    buf_alloc (vb, &vb->section_list_buf, 1, 50, SectionEnt, 2, "section_list_buf");
    
    int header_size = st_header_size (header->section_type);
    if (header->data_encrypted_len) header_size = ROUNDUP16 (header_size);

    BNXT (SectionEntModifiable, vb->section_list_buf) = (SectionEntModifiable) {
        .st        = header->section_type,
        .vblock_i  = (IS_VB_SEC(st) || IS_FRAG_SEC(st)) ? BGEN32 (header->vblock_i) : 0, // big endian in header - convert back to native
        .comp_i    = IS_COMP_SEC(st) ? vb->comp_i : COMP_NONE,
        .offset    = vb->z_data.len,  // this is a partial offset (within vb) - we will correct it later
        .flags     = header->flags,
        .num_lines = ST(VB_HEADER) ? vb->lines.len32 : 0, 
        .num_deep_lines = (ST(VB_HEADER) && (VB_DT(SAM) || VB_DT(BAM))) ? vb->lines.count : 0, // alignments that are not SUPPLEMENTARY or SECONDARY
        .size      = header_size  + MAX_(BGEN32 (header->data_compressed_len), BGEN32 (header->data_encrypted_len))
        //up to v14: header_size  + BGEN32 (header->data_compressed_len)
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
void sections_list_concat (BufferP section_list_buf)
{
    ARRAY (SectionEntModifiable, vb_sec, *section_list_buf);

    // update the offset
    for (uint32_t i=0; i < vb_sec_len; i++) 
        vb_sec[i].offset += z_file->disk_so_far;

    // copy all entries. note: z_file->section_list_buf is protected by zriter_mutex 
    buf_add_buf (evb, &z_file->section_list_buf, section_list_buf, SectionEntModifiable, "z_file->section_list_buf");

    section_list_buf->len = 0;
}

// section iterator. returns true if a section of this type was found.
bool sections_next_sec3 (Section *sl_ent,   // optional in/out. if NULL - search entire list
                         SectionType st1, SectionType st2, SectionType st3) // check only next section, not entire remaining list
{
    ASSERTNOTEMPTY (z_file->section_list_buf);
    
    Section sec = sl_ent ? *sl_ent : NULL; 
    bool found = false;

    ASSERT (!sec || BISVALID (z_file->section_list_buf, sec), "Invalid sec: st1=%s st2=%s st3=%s", 
            st_name (st1), st_name (st2), st_name (st3));

    while (sec < BAFT (SectionEnt, z_file->section_list_buf) - 1) {

        sec = sec ? (sec + 1) : B1ST (SectionEnt, z_file->section_list_buf); 

        if (sec->st == st1 || sec->st == st2 || sec->st == st3) {
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

// -------------------
// VBs
// -------------------

// PIZ: called twice: once when z_file is opened, and once from sections_commit_new_list (after some VBs have been removed by the recon_plan)
// ZIP: called after each txt_file is zipped into a z_file
void sections_create_index (void)
{
    START_TIMER;

    ASSERTMAINTHREAD;

    if (Z_DT(REF)) return; // reference files have no components or VBs (save time)

    // free existing index in case it is being re-created (after each txt_file is compressed)
    if (IS_ZIP) {
        buf_free (z_file->comp_sections_index); 
        buf_free (z_file->vb_sections_index);
    }

    ASSERTNOTZERO (z_file->num_components);
    ARRAY_alloc (SectionsVbIndexEnt, vb_index, z_file->num_vbs + 1, true, z_file->vb_sections_index, evb, "z_file->vb_sections_index"); // +1 bc vb_i is 1-based
    ARRAY_alloc (SectionsCompIndexEnt, comp_index, z_file->num_components, true, z_file->comp_sections_index, evb, "z_file->comp_sections_index");

    // initialize
    for (int i=0; i < vb_index_len; i++)
        vb_index[i].vb_header_sec_i = vb_index[i].last_sec_i = -1;

    for (int i=0; i < comp_index_len; i++)
        comp_index[i].bgzf_sec_i = comp_index[i].txt_header_sec_i = 
        comp_index[i].recon_plan_sec_i[0] = comp_index[i].recon_plan_sec_i[1] = -1;

    // populate index
    SectionsVbIndexEnt *vbinx = NULL; // if non-NULL, we're within a VB section block
    for_buf2 (SectionEnt, sec, sec_i, z_file->section_list_buf) {
        if (sec->st == SEC_B250 || sec->st == SEC_LOCAL) continue; // short circuit

        SectionsCompIndexEnt *comp = &comp_index[sec->comp_i];

        // case: we passed the last B250/LOCAL section of the VB 
        if (vbinx) {
            vbinx->last_sec_i = sec_i - 1;
            vbinx = NULL;
        }

        switch (sec->st) {
            case SEC_TXT_HEADER : 
                if (comp->txt_header_sec_i == -1) // store first fragment of txt header
                    comp->txt_header_sec_i = sec_i; 
                break;
            
            case SEC_BGZF : 
                comp->bgzf_sec_i = sec_i; 
                break;

            case SEC_RECON_PLAN : {
                bool is_luft_plan = sec->flags.recon_plan.luft;
                if (comp->recon_plan_sec_i[is_luft_plan] == -1) // save first fragment
                    comp->recon_plan_sec_i[is_luft_plan] = sec_i; 
                break;
            }

            // note: in ZIP, VBs are written out-of-order, so some VBs might still not be written and hence their index entry will remain empty
            case SEC_VB_HEADER : {
                VBIType vb_i = sec->vblock_i;
                
                ASSERT (vb_i>=1 && vb_i < z_file->vb_sections_index.len32, "Bad vb_i=%u, z_file->vb_sections_index.len=%u", vb_i, z_file->vb_sections_index.len32);
                
                vbinx = &vb_index[vb_i];
                *vbinx = (SectionsVbIndexEnt){ .vb_header_sec_i = sec_i };

                if (!comp->first_vb_i || vb_i < comp->first_vb_i) 
                    comp->first_vb_i = vb_i; // first vb_i in numerical order of vb_i

                if (!comp->first_vb_of_comp)
                    comp->first_vb_of_comp = vbinx; // first VB in the order it appears in the section list (not necessarily the lowest vb_i)

                else // link to previous VB in this component (consecutive in section list, but not necessarilty consecutive vb_i)
                    comp->last_vb_of_comp->next_vb_of_comp = vbinx;
                    
                comp->last_vb_of_comp = vbinx; // last so far
                comp->num_vbs++;
                break;
            }

            default: break;
        }
    }

    COPY_TIMER_EVB (sections_create_index);
}

// ZIP / PIZ: returns the SEC_VB_HEADER section, or NULL if this vb_i doesn't exist
Section sections_vb_header (VBIType vb_i)
{    
    uint32_t vb_header_sec_i = Bvbindex(vb_i)->vb_header_sec_i;
    Section sec = Bsec(vb_header_sec_i);
    
    // sanity
    ASSERT (sec, "vb_i=%u fits in z_file->vb_sections_index, but it is not indexed", vb_i);
    ASSERT (sec->st == SEC_VB_HEADER, "Expecting indexed VB section of vb_i=%u to have st=VB_HEADER but it has st=%s", vb_i, st_name(sec->st));

    return sec;
}

Section sections_one_before (Section sec) 
{ 
    return (!sec || sec == B1ST (SectionEnt, z_file->section_list_buf)) ? NULL : sec-1; 
}

//---------------------------------------------
// PIZ: constructing a updated section list
//---------------------------------------------

void sections_new_list_add_vb (BufferP new_list, VBIType vb_i)
{
    const SectionsVbIndexEnt *vbent = Bvbindex(vb_i);
    uint32_t num_sections = vbent->last_sec_i - vbent->vb_header_sec_i + 1;

    memcpy (BAFT (SectionEntModifiable, *new_list), Bsec(vbent->vb_header_sec_i), num_sections * sizeof (SectionEnt));
    new_list->len += num_sections;
}

void sections_new_list_add_txt_header (BufferP new_list, CompIType comp_i)
{
    // add all fragments of the txt_header of this component
    for (Section sec = sections_get_comp_txt_header_sec (comp_i); 
         sec->st == SEC_TXT_HEADER && sec->comp_i == comp_i; 
         sec++)
        BNXT (SectionEntModifiable, *new_list) = *sec;
}

// PIZ: If any of the components has a SEC_BGZF add it
void sections_new_list_add_bgzf (BufferP new_list)
{
    for_buf2 (SectionsCompIndexEnt, comp, comp_i, z_file->comp_sections_index)
        if (comp->bgzf_sec_i != -1) 
            BNXT (SectionEntModifiable, *new_list) = *Bsec(comp->bgzf_sec_i);
}

void sections_new_list_add_global_sections (BufferP new_list)
{
    // get first section that's not TXT_HEADER/BGZF/VB_HEADER/LOCAL/B250/COUNT/DICT
    Section sec = NULL;
    for_buf_back (SectionEnt, s, z_file->section_list_buf) 
        if (IS_DICTED_SEC(s->st) || s->st == SEC_VB_HEADER || s->st == SEC_TXT_HEADER || s->st == SEC_BGZF) {
            sec = s; 
            break;
        }
    sec++;

    buf_append (NULL, *new_list, SectionEntModifiable, sec, BAFT(SectionEnt, z_file->section_list_buf) - sec, NULL);
}

// PIZ: replace current section list with the new list (if one exists)
void sections_commit_new_list (BufferP new_list)
{
    if (!new_list->len) return;

    ASSERT (new_list->len <= z_file->section_list_buf.len, "Expecting new_list->len=%"PRIu64" <= z_file->section_list_buf.len=%"PRIu64, 
            new_list->len, z_file->section_list_buf.len);

    buf_destroy (z_file->section_list_buf);
    buf_grab (evb, z_file->section_list_buf, "z_file->section_list_buf", *new_list);

    sections_create_index();
}

//---------------
// misc functions
//---------------

// up to v14
static void v14_sections_list_file_to_memory_format (void)
{
    ARRAY (const SectionEntFileFormatV14, file_sec, z_file->section_list_buf); // file format

    for_buf2_back (SectionEntModifiable, sec, i, z_file->section_list_buf) {
        SectionEntFileFormatV14 fsec = file_sec[i]; // make a copy as assignment is overlapping
        
        *sec = (SectionEntModifiable){
            .offset      = BGEN64 (fsec.offset),
            .st_specific = fsec.st_specific,  
            .vblock_i    = BGEN32 (fsec.vblock_i),
            .st          = fsec.st,
            .comp_i      = (VER(14) && IS_COMP_SEC(fsec.st)) ? fsec.comp_i : COMP_NONE, // note: this field was introduced in v14: COMP_NONE sections have 0, as comp_i is just 2 bit. 
            .flags       = VER(12) ? fsec.flags  : (SectionFlags){} // prior to v12, flags were stored only in SectionHeader. Since v12, they are also copied to SectionEntFileFormat.
        };

        if (sec->st == SEC_VB_HEADER) {
            sec->num_lines  = BGEN32 (sec->num_lines);        
            z_file->num_vbs = MAX_(z_file->num_vbs, sec->vblock_i);

            if (sec->comp_i != COMP_NONE) {
                z_file->comp_num_lines[sec->comp_i] += sec->num_lines;
                z_file->num_components = MAX_(z_file->num_components, sec->comp_i);
            }
        }    

        else if (sec->st == SEC_TXT_HEADER && sec->comp_i != COMP_NONE)
            z_file->num_components = MAX_(z_file->num_components, sec->comp_i);

        if (i < file_sec_len-1)
            sec->size = (sec+1)->offset - sec->offset;
    }

    // comp_i was introduced V14, for V<=13, get comp_i by the component's relative position in the file. 
    if (VER(14)) 
        z_file->num_components++; // one more than the max comp_i

    else {
        CompIType comp_i_by_consecutive = COMP_NONE; 

        for_buf (SectionEntModifiable, sec, z_file->section_list_buf) {
            if (sec->st == SEC_TXT_HEADER) comp_i_by_consecutive++; // note: we only read bound V<=13 files if the are a single FQ pair, or a DVCF

            if (IS_COMP_SEC(sec->st))
                sec->comp_i = comp_i_by_consecutive; 

            // Up to V11, flags were stored only in SectionHeader
            if (!VER(12)) 
                sec->flags = zfile_read_section_header (evb, sec, SEC_NONE).common.flags;
        }

        z_file->num_components = comp_i_by_consecutive + 1; // one more than the max comp_i
    }
}

// ZIP (see also bug 819)
void sections_list_memory_to_file_format (bool in_place) // in place, or to evb->scratch for testing codec
{
    if (!in_place) { // just generating for codec testing
        ASSERTNOTINUSE(evb->scratch);
        buf_alloc_exact (evb, evb->scratch, (MIN_(z_file->section_list_buf.len * sizeof (SectionEntFileFormat), CODEC_ASSIGN_SAMPLE_SIZE) / sizeof(SectionEntFileFormat)),
                         SectionEntFileFormat, "scratch");
    }

    BufferP out = in_place ? &z_file->section_list_buf : &evb->scratch;

    // replace dict_id with the the sec_i of its first appearahce. 
    int32_t first_appearance[z_file->num_contexts];
    memset (first_appearance, 255, z_file->num_contexts * sizeof (int32_t));

    SectionEntModifiable prev_sec = {};
    uint32_t prev_num_lines = 0;
    for_buf2 (SectionEntFileFormat, fsec, i, *out) {
        SectionEnt sec = *B(SectionEnt, z_file->section_list_buf, i); // copy before it gets overwritten
        
        int64_t offset_delta = (int64_t)sec.offset - (int64_t)prev_sec.offset;
        ASSERT (offset_delta >=0LL && offset_delta <= 0xffffffffLL,  // note: offset_delta is size of previous section
                "section_i=%u size=%"PRId64" st=%s is too big", i-1, offset_delta, st_name ((fsec-1)->st));

        int32_t vb_delta = INTERLACE(int32_t, (int32_t)sec.vblock_i - (int32_t)prev_sec.vblock_i);

        *fsec = (SectionEntFileFormat){
            .vblock_i_delta = BGEN32 (vb_delta),
            .comp_i_plus_1  = (i && sec.comp_i == prev_sec.comp_i) ? 0 
                            : (sec.comp_i == COMP_NONE)            ? COMP_NONE
                            :                                        (1 + sec.comp_i),  
            .offset_delta   = BGEN32 ((uint32_t)offset_delta),
            .st             = sec.st,
            .flags          = sec.flags    // since v12
        }; 

        if (IS_DICTED_SEC(sec.st)) {
            Did did_i = ctx_get_existing_zctx (sec.dict_id)->did_i;
            if (first_appearance[did_i] == -1) {
                fsec->dict_id = sec.dict_id;
                first_appearance[did_i] = i;
            }
            else 
                fsec->dict_sec_i = BGEN32 (first_appearance[did_i]);
        }

        // To do: num_lines - delta from prev_num_lines ; num_deep_lines non-zero only if deep&SAM - -delta vs num_lines
        // tested, but these make BZ2 compression worse
        else if (sec.st == SEC_VB_HEADER) {
            int32_t num_lines_delta = INTERLACE(int32_t, (int32_t)sec.num_lines - (int32_t)prev_num_lines);
            prev_num_lines = sec.num_lines;

            fsec->num_lines = BGEN32 (num_lines_delta);
            if (flag.deep && (sec.comp_i == SAM_COMP_MAIN || sec.comp_i == SAM_COMP_PRIM))
                fsec->num_non_deep_lines = BGEN32(sec.num_lines - sec.num_deep_lines);
        }

        prev_sec = sec;
    }

    out->len *= sizeof (SectionEntFileFormat); // change to counting bytes
}

// PIZ: create in-memory format of the section list - copy from z_file section, BGEN, and add sizes
void sections_list_file_to_memory_format (SectionHeaderGenozipHeaderP genozip_header)
{
    struct FlagsGenozipHeader f = genozip_header->flags.genozip_header;
    DataType dt = BGEN16 (genozip_header->data_type);

    // For files v13 or earlier, we can only read them if they are a single component, or two paired FASTQs, or DVCF
    // Other bound files need to be decompressed with Genozip v13
    uint32_t v13_num_components = BGEN32 (genozip_header->v13_num_components);
    ASSINP (VER(14)                         || 
            v13_num_components == 1         ||   // single file
            (dt == DT_VCF && f.has_gencomp) ||   // DVCF
            (dt == DT_FASTQ && f.v14_dts_paired && v13_num_components == 2) || // Paired FASTQs (the v14_dts_paired flag was introduced in V9.0.13)
            is_genols                       ||
            flags_is_genocat_global_area_only(), // only show meta-data (can't use flag.genocat_global_area_only bc flags are not set yet)
            "%s is comprised of %u %s files bound together. The bound file feature was discontinued in Genozip v14. To decompress this file, use Genozip v13",
            z_name, v13_num_components, dt_name (BGEN16(genozip_header->data_type)));
    
    z_file->section_list_buf.len /= VER(15) ? sizeof (SectionEntFileFormat) : sizeof (SectionEntFileFormatV14); // fix len

    buf_alloc (evb, &z_file->section_list_buf, 0, z_file->section_list_buf.len, SectionEnt, 0, NULL); // extend

    if (VER(15)) {
        ARRAY (const SectionEntFileFormat, file_sec, z_file->section_list_buf);  // file format

        // note: we work backwards as mem_sec items are larger than file_sec
        for_buf2_back (SectionEntModifiable, sec, i, z_file->section_list_buf) { // memory format (entries are larger)
            SectionEntFileFormat fsec = file_sec[i]; // make a copy as assignment is overlapping
            
            *sec = (SectionEntModifiable){
                .offset      = BGEN32(fsec.offset_delta),
                .st_specific = fsec.st_specific,  
                .vblock_i    = DEINTERLACE (int32_t, (int32_t)BGEN32 (fsec.vblock_i_delta)),
                .st          = fsec.st,
                .comp_i      = fsec.comp_i_plus_1, 
                .flags       = fsec.flags
            };
        }

        // finalize values  
        int32_t prev_num_lines = 0;
        for_buf2 (SectionEntModifiable, sec, i, z_file->section_list_buf) {
            if (i) (sec-1)->size = sec->offset;
            if (i) sec->offset += (sec-1)->offset;

            sec->comp_i = (sec->comp_i == 0)         ? (sec-1)->comp_i
                        : (sec->comp_i == COMP_NONE) ? COMP_NONE
                        :                              (sec->comp_i - 1);

            sec->vblock_i = i ? (int32_t)(sec-1)->vblock_i + (int32_t)sec->vblock_i : (int32_t)sec->vblock_i;

            if (IS_DICTED_SEC(sec->st) && !sec->is_dict_id) 
                sec->dict_id = (sec - i + BGEN32(sec->dict_sec_i))->dict_id; // copy from section of first appearance of this dict_id

            else if (sec->st == SEC_VB_HEADER) {
                sec->num_lines      = prev_num_lines + DEINTERLACE(int32_t, (int32_t)BGEN32 (sec->num_lines/*this is still delta_num_lines*/));
                prev_num_lines      = sec->num_lines;
                sec->num_deep_lines = sec->num_lines - BGEN32 (sec->num_deep_lines/*this is still num_non_deep_lines*/);
                z_file->num_vbs     = MAX_(z_file->num_vbs, sec->vblock_i);

                z_file->comp_num_lines[sec->comp_i] += sec->num_lines;
            }    

            if (sec->st == SEC_VB_HEADER || sec->st == SEC_TXT_HEADER)
                z_file->num_components = MAX_(z_file->num_components, sec->comp_i);
        }
    
        z_file->num_components++; // one more than the max comp_i
    }

    else  // up to v14
        v14_sections_list_file_to_memory_format();

    BLST(SectionEntModifiable, z_file->section_list_buf)->size =
        st_header_size (SEC_GENOZIP_HEADER) + BGEN32 (genozip_header->data_compressed_len) + sizeof (SectionFooterGenozipHeader);

    sections_create_index(); // PIZ initial indexing, we will index again after recon_plan potentially removes some VBs
    
    // Up to V11, top_level_repeats existed (called num_lines) but was not populated
    // in V12-13 num_lines was transmitted through the VbHeader.top_level_repeats
    // Since V14, num_lines is transmitted through SectionEntFileFormat
    if (VER(12) && !VER(14)) 
        for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {
            SectionEntModifiable *sec = (SectionEntModifiable *)sections_vb_header (vb_i);
            SectionHeaderUnion header = zfile_read_section_header (evb, sec, SEC_VB_HEADER);
            sec->num_lines = BGEN32 (header.vb_header.v13_top_level_repeats);
        }

    // save sections to section_list_save, because they are needed for setting the default context
    // flags, and vb_i=1 may be eliminated by the recon_plan (for example, FASTQ with --R2). Also, we need to recover between txt_files.
    buf_copy (evb, &z_file->section_list_save, &z_file->section_list_buf, SectionEnt, 0, 0, "z_file->section_list_save");

    if (z_file->num_vbs >= 1) {
        z_file->section_list_save.prm32[0] = Bvbindex(1)->vb_header_sec_i; 
        z_file->section_list_save.prm32[1] = Bvbindex(1)->last_sec_i;
    }
    else
        z_file->section_list_save.param = 0;
}

// check if there is any DICT, LOCAL, B250 or COUNT section of a certain dict_id
bool is_there_any_section_with_dict_id (DictId dict_id)
{
    for_buf (SectionEnt, sec, z_file->section_list_buf)
        if (sec->dict_id.num == dict_id.num && IS_DICTED_SEC (sec->st))
            return true;

    return false;
}

rom st_name (SectionType sec_type)
{
    static char invalid[24]; // not thread safe, but not expected except in error situations
    if (sec_type < SEC_NONE || sec_type >= NUM_SEC_TYPES) {
        sprintf (invalid, "INVALID(%d)", sec_type);
        return invalid;
    }

    return (sec_type == SEC_NONE) ? "SEC_NONE" : type_name (sec_type, &abouts[sec_type].name , ARRAY_LEN(abouts));
}

StrText comp_name_(CompIType comp_i)
{
    static int max_comps_by_dt[NUM_DATATYPES] = { [DT_VCF]=3, [DT_BCF]=3, [DT_SAM]=3/*except if --deep*/, [DT_BAM]=3, [DT_FASTQ]=2 };
    
    static rom comp_names[NUM_DATATYPES][5]   = { [DT_VCF] = VCF_COMP_NAMES, [DT_BCF] = VCF_COMP_NAMES,
                                                  [DT_SAM] = SAM_COMP_NAMES, [DT_BAM] = SAM_COMP_NAMES,
                                                  [DT_FASTQ] = FASTQ_COMP_NAMES };
    
    DataType dt = z_file ? z_file->data_type : DT_NONE;
    StrText s;

    if (!z_file) // can happen if another thread is busy exiting the process
        strcpy (s.s, "EMEM");
    
    else if (!max_comps_by_dt[dt] && comp_i==0)
        strcpy (s.s, "MAIN");

    else if (comp_i==COMP_NONE)
        strcpy (s.s, "NONE");

    else if ((dt==DT_SAM || dt==DT_BAM) && comp_i >= SAM_COMP_FQ00) 
        sprintf (s.s, "FQ%02u", comp_i - SAM_COMP_FQ00); 

    else if (comp_i >= max_comps_by_dt[dt]) 
        sprintf (s.s, "invalid-comp_i-%d", comp_i);

    else 
        strcpy (s.s, comp_names[dt][comp_i]);

    return s;
}

static StrText comp_name_ex (CompIType comp_i, SectionType st)
{
    if (ST(RECON_PLAN) || ST(TXT_HEADER)) {
        StrText s;
        sprintf (s.s, "%s%02u",ST(RECON_PLAN) ? "RP" : "TX", comp_i);
        return s;
    }
        
    else {
        if (IS_COMP_SEC(st)) return comp_name_(comp_i);
        if (ST(REFERENCE))   return (StrText){ "REFR" };
        if (ST(REF_IS_SET))  return (StrText){ "REFR" };
        if (ST(REF_HASH))    return (StrText){ "REFH" };
        if (ST(DICT))        return (StrText){ "DICT" };
        else                 return (StrText){ "----" };
    }
}

StrText vb_name (VBlockP vb)
{
    StrText s;
    if (vb && vb->vblock_i)
        sprintf (s.s, "%.10s/%u", comp_name (vb->comp_i), vb->vblock_i);
    else if (vb && segconf.running)
        strcpy (s.s, "SEGCONF/0");
    else
        strcpy (s.s, "NONCOMPUTE");

    return s;
}

StrText line_name (VBlockP vb)
{
    StrText s;
    sprintf (s.s, "%.10s/%u%s", vb_name(vb).s, vb->line_i, vb->preprocessing ? "(preproc)" : "");
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
         : sec_type == SEC_VB_HEADER      && command != ZIP && !VER(12) ? sizeof (SectionHeaderVbHeader) - 5*sizeof(uint32_t) // in v8-11, SectionHeaderVbHeader was shorter by 5 32b words
         : sec_type == SEC_VB_HEADER      && command != ZIP && !VER(14) ? sizeof (SectionHeaderVbHeader) - 2*sizeof(uint32_t) // in v12-13, SectionHeaderVbHeader was shorter by 2 32b word
         : sec_type == SEC_VB_HEADER      && command != ZIP && !VER(15) ? sizeof (SectionHeaderVbHeader) - 1*sizeof(uint32_t) // in v14, SectionHeaderVbHeader was shorter by 1 32b word
         : sec_type == SEC_TXT_HEADER     && command != ZIP && !VER(12) ? sizeof (SectionHeaderTxtHeader) - 3*sizeof(uint16_t) - sizeof(uint64_t) // in v8-11, SectionHeaderTxtHeader was shorter
         : sec_type == SEC_TXT_HEADER     && command != ZIP && !VER(15) ? sizeof (SectionHeaderTxtHeader) - 3*sizeof(uint16_t) // in v12-14, SectionHeaderTxtHeader was shorter by 3XQnameFlavorProp
         : sec_type == SEC_GENOZIP_HEADER && command != ZIP && !VER(12) ? sizeof (SectionHeaderGenozipHeader) - REF_FILENAME_LEN - sizeof(Digest) // in v8-11, SectionHeaderTxtHeader was shorter
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
Section sections_first_sec (SectionType st, FailType soft_fail)
{
    ARRAY (SectionEnt, sec, z_file->section_list_buf);

    for (unsigned i=0; i < sec_len; i++)
        if (sec[i].st == st) return &sec[i];

    ASSERT (soft_fail, "Cannot find section_type=%s in z_file", st_name (st));

    return NULL;
}

// called by ZIP, PIZ main thread
Section sections_last_sec (SectionType st, FailType soft_fail)
{
    ARRAY (SectionEnt, sec, z_file->section_list_buf);

    for (int i=sec_len-1; i >= 0; i--)
        if (sec[i].st == st) return &sec[i];

    ASSERT (soft_fail, "Cannot find section_type=%s in z_file", st_name (st));

    return NULL;
}

VBIType sections_get_num_vbs (CompIType comp_i) 
{ 
    if (comp_i >= z_file->comp_sections_index.len32) return 0; // this component doesn't exist (can happen eg when checking for gencomp that doesn't exist in this file)

    return Bcompindex(comp_i)->num_vbs;
}

VBIType sections_get_num_vbs_(CompIType first_comp_i, CompIType last_comp_i)
{
    ASSERT (first_comp_i != COMP_NONE && last_comp_i != COMP_NONE, "COMP_NONE unexpected: first_comp_i=%s last_comp_i=%s",
            comp_name (first_comp_i), comp_name (last_comp_i));
            
    VBIType num_vbs = 0;
    for (CompIType comp_i = first_comp_i; comp_i <= last_comp_i; comp_i++)
        num_vbs += sections_get_num_vbs (comp_i);

    return num_vbs;
}

VBIType sections_get_first_vb_i (CompIType comp_i) 
{ 
    VBIType first_vb_i = Bcompindex(comp_i)->first_vb_i;
    ASSERT (first_vb_i, "comp_i=%s has no VBs", comp_name (comp_i));
    return first_vb_i;
}

Section sections_get_comp_txt_header_sec (CompIType comp_i)
{
    return Bsec(Bcompindex(comp_i)->txt_header_sec_i);
}

Section sections_get_comp_bgzf_sec (CompIType comp_i)
{
    return Bsec(Bcompindex(comp_i)->bgzf_sec_i);
}

Section sections_get_comp_recon_plan_sec (CompIType comp_i, bool is_luft_plan)
{
    return Bsec(Bcompindex(comp_i)->recon_plan_sec_i[is_luft_plan]);
}

// get next VB_HEADER in the order they appear in the section list (not necessarily consecutive vb_i)
Section sections_get_next_vb_header_sec (CompIType comp_i, Section *vb_header_sec)
{
    const SectionsCompIndexEnt *comp_index_ent = Bcompindex(comp_i);

    if (! *vb_header_sec)
        *vb_header_sec = comp_index_ent->first_vb_of_comp ? 
            Bsec(comp_index_ent->first_vb_of_comp->vb_header_sec_i) : NULL/*component has no VBs*/;
    
    else {
        SectionsVbIndexEnt *next_vb_ent = Bvbindex((*vb_header_sec)->vblock_i)->next_vb_of_comp;
        *vb_header_sec = next_vb_ent ? Bsec (next_vb_ent->vb_header_sec_i) : NULL;
    }

    ASSERT (!(*vb_header_sec) || (*vb_header_sec)->st == SEC_VB_HEADER, 
            "expecting section to have a SEC_VB_HEADER but it has %s", st_name ((*vb_header_sec)->st));

    return *vb_header_sec;
}

// inspects z_file flags and if needed reads additional data, and returns true if the z_file consists of FASTQs compressed with --pair
bool sections_is_paired (void)
{
    bool is_paired;

    if (VER(15)) 
        is_paired = (Z_DT(FASTQ) && z_file->comp_sections_index.len32 == 2) ||
                    (Z_DT(SAM)   && z_file->comp_sections_index.len32 == 5 && 
                    Bsec(Bcompindex(SAM_COMP_FQ01)->txt_header_sec_i)->flags.txt_header.pair);

    // up to V14, only FASTQ files could be paired. back comp note: until v13, it was possible to have bound files with multiple pairs (is this true?). not supported here.
    else if (!Z_DT(FASTQ) || z_file->num_txt_files != 2)  
        is_paired = false;

    // from V10 to V14 - determine by v14_dts_paired (introduced in 9.0.13) 
    else 
        is_paired = z_file->z_flags.v14_dts_paired;
    
    // sanity check: R1 and R2 are expected to have the same number of VBs (only if index is already created)
    if (is_paired && z_file->comp_sections_index.len32) {
        uint32_t num_vbs_R1 = Bcompindex(z_file->comp_sections_index.len32-2)->num_vbs;
        uint32_t num_vbs_R2 = Bcompindex(z_file->comp_sections_index.len32-1)->num_vbs;
        ASSERT (num_vbs_R1 == num_vbs_R2, "R1 and R2 have a different number of VBs (%u vs %u respecitvely)", num_vbs_R1, num_vbs_R2); 
    }

    return is_paired;
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
    rom dts[NUM_DATATYPES]  = { [DT_FASTQ]=(!VER(15) ? "v14_dts_paired" : 0), [DT_SAM]="dts_ref_internal", [DT_CHAIN]="dts_mismatch" };
    rom dts2[NUM_DATATYPES] = { [DT_SAM]="dts2_deep" };

    FlagStr str = {};

    switch (st) {
        case SEC_GENOZIP_HEADER:
            sprintf (str.s, "%s=%u %s=%u aligner=%u txt_is_bin=%u %s=%u adler=%u has_gencomp=%u has_taxid=%u",
                     (dt != DT_NONE && dts[dt])  ? dts[dt]  : "dt_specific",  f.genozip_header.dt_specific, 
                     (dt != DT_NONE && dts2[dt]) ? dts2[dt] : "dt_specific2", f.genozip_header.dt_specific2, 
                     f.genozip_header.aligner, f.genozip_header.txt_is_bin, 
                     VER(15) ? "has_digest" : "v14_bgzf", f.genozip_header.has_digest, 
                     f.genozip_header.adler, f.genozip_header.has_gencomp,f.genozip_header.has_taxid);
            break;

        case SEC_TXT_HEADER: {
            char extra[64] = {};
            if (dt==DT_VCF && !VER(14)) sprintf (extra, " dvcf_comp_i=%u", f.txt_header.v13_dvcf_comp_i);
            if ((dt==DT_SAM || dt==DT_BAM || dt==DT_FASTQ) && VER(15)) sprintf (extra, " pair=%s", pair_type_name (f.txt_header.pair));
            sprintf (str.s, "is_txt_luft=%u no_gz_ext=%u %s", f.txt_header.is_txt_luft, f.txt_header.no_gz_ext, extra);
            break;
        }

        case SEC_VB_HEADER:
            switch (dt) {
                case DT_VCF: 
                    sprintf (str.s, "coords=%s null_DP=%u", vcf_coords_name (f.vb_header.vcf.coords), f.vb_header.vcf.use_null_DP_method);
                    break;
                
                case DT_SAM: case DT_BAM:
                    if (VER(13) && !VER(14)) sprintf (str.s, "sorted=%u collated=%u", f.vb_header.sam.v13_is_sorted, f.vb_header.sam.v13_is_collated);
                    break;
                
                case DT_GFF:
                    if (VER(15)) sprintf (str.s, "embedded_fasta=%u", f.vb_header.gff.embedded_fasta);
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

        case SEC_DICT:
            sprintf (str.s, "deep_sam=%u deep_fq=%u", f.dictionary.deep_sam, f.dictionary.deep_fastq);
            break;

        default: 
            str.s[0] = 0;
    }

    return str;
}

void sections_show_header (ConstSectionHeaderP header, VBlockP vb /* optional if output to buffer */, uint64_t offset, char rw)
{
    #define DT(x) ((dt) == DT_##x)

    if (flag_loading_auxiliary && !flag.debug_read_ctxs) return; // don't show headers of an auxiliary file in --show-header, but show in --debug-read-ctx
    
    if (  flag.show_headers   != SHOW_ALL_HEADERS &&     // we don't need to show all sections
          flag.show_headers-1 != header->section_type && // we don't need to show this section
          !flag.debug_read_ctxs)
        return;

    bool is_dict_offset = (HEADER_IS(DICT) && rw == 'W'); // at the point calling this function in zip, SEC_DICT offsets are not finalized yet and are relative to the beginning of the dictionary area in the genozip file
    bool v12 = (IS_ZIP || VER(12));
    bool v14 = (IS_ZIP || VER(14));
    bool v15 = (IS_ZIP || VER(15));

    char str[2048];
    #define PRINT { if (vb) buf_append_string (vb, &vb->show_headers_buf, str); else iprintf ("%s", str); } 
    
    rom SEC_TAB = isatty(1) ? "            ++  " : " "; // single line if not terminal - eg for grepping

    sprintf (str, "%c %s%-*"PRIu64" %-19s %-4.4s %-4.4s vb=%-3u %s=%*u txt_len=%-7u z_len=%-7u enc_len=%-7u%s ",
             rw, 
             is_dict_offset ? "~" : "", 9-is_dict_offset, offset, 
             st_name(header->section_type), 
             codec_name (header->codec), codec_name (header->sub_codec),
             BGEN32 (header->vblock_i), 
             VER(15) ? "z_digest" : "comp_off",
             VER(15) ? -10 : -4,
             VER(15) ? BGEN32 (header->z_digest) : BGEN32 (header->v14_compressed_offset), 
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
        SectionHeaderGenozipHeaderP h = (SectionHeaderGenozipHeaderP)header;
        z_file->z_flags.adler = h->flags.genozip_header.adler; // needed by digest_display_ex
        char dt_specific[REF_FILENAME_LEN + 200] = "";
        dt = BGEN16 (h->data_type); // for GENOZIP_HEADER, go by what header declares
        if (dt >= NUM_DATATYPES) dt = DT_NONE;

        if (DT(CHAIN))
            sprintf (dt_specific, "%sprim=\"%.*s\" md5=%s\n", 
                     SEC_TAB, REF_FILENAME_LEN, h->chain.prim_filename, 
                     digest_display_ex (h->chain.prim_genome_digest, DD_MD5).s);

        else if ((DT(VCF) || DT(BCF)) && v14)
            sprintf (dt_specific, "%ssegconf=(has_RGQ=%s,GQ_method=%s,FMT_DP_method=%s) width=(AC=%u,AN=%u,MLEAC=%u,DP=%u,QD=%u,SF=%u) max_ploidy_for_mux=%u\n", 
                     SEC_TAB, TF(h->vcf.segconf_has_RGQ), GQ_method_name (h->vcf.segconf_GQ_method), FMT_DP_method_name (h->vcf.segconf_FMT_DP_method), 
                     h->vcf.width.AC, h->vcf.width.AN, h->vcf.width.MLEAC, h->vcf.width.DP, h->vcf.width.QD, h->vcf.width.SF, 
                     h->vcf.max_ploidy_for_mux);

        else if ((DT(SAM) || DT(BAM)) && v14)
            sprintf (dt_specific, "%ssegconf=(sorted=%u,collated=%u,seq_len=%u,seq_len_to_cm=%u,ms_type=%u,has_MD_or_NM=%u,bisulfite=%u,MD_NM_by_unconverted=%u,predict_meth=%u,is_paired=%u,sag_type=%s,sag_has_AS=%u,pysam_qual=%u,cellranger=%u,SA_HtoS=%u,seq_len_dict_id=%s,deep_qname1=%u,deep_qname2=%u,deep_no_qual=%u,use_ins_ctx=%u,%uXsam_factor=%u)\n", 
                     SEC_TAB, h->sam.segconf_is_sorted, h->sam.segconf_is_collated, BGEN32 (h->sam.segconf_seq_len), h->sam.segconf_seq_len_cm, h->sam.segconf_ms_type, h->sam.segconf_has_MD_or_NM, 
                     h->sam.segconf_bisulfite, h->sam.segconf_MD_NM_by_un, h->sam.segconf_predict_meth, 
                     h->sam.segconf_is_paired, sag_type_name(h->sam.segconf_sag_type), h->sam.segconf_sag_has_AS, 
                     h->sam.segconf_pysam_qual, h->sam.segconf_cellranger, h->sam.segconf_SA_HtoS, dis_dict_id(h->sam.segconf_seq_len_dict_id).s,
                     h->sam.segconf_deep_qname1, h->sam.segconf_deep_qname2, h->sam.segconf_deep_no_qual, 
                     h->sam.segconf_use_ins_ctxs, SAM_FACTOR_MULT, h->sam.segconf_sam_factor);

        else if (DT(REF)) {
            if (v15) sprintf (dt_specific, "%sgenome_digest=%s\n", SEC_TAB, digest_display (h->genome_digest).s);
            else     sprintf (dt_specific, "%sfasta_md5=%s\n", SEC_TAB, digest_display (h->v14_REF_fasta_md5).s);
        }

        else if (DT(FASTQ) && v14) 
            sprintf (dt_specific, "%sFASTQ_v13_digest_bound=%s segconf_seq_len_dict_id=%s\n", 
                     SEC_TAB, digest_display (h->FASTQ_v13_digest_bound).s, dis_dict_id(h->fastq.segconf_seq_len_dict_id).s);

        sprintf (str, "\n%sver=%u.0.%u private=%u enc=%s dt=%s usize=%"PRIu64" lines=%"PRIu64" secs=%u txts=%u vb_size=%u\n" 
                      "%s%s %s=\"%.*s\" %s=%s\n"
                      "%s" // dt_specific, if there is any
                      "%screated=\"%.*s\"\n",
                 SEC_TAB, h->genozip_version, h->genozip_minor_ver/*15.0.28*/, h->private_file, encryption_name (h->encryption_type), dt_name (dt), 
                 BGEN64 (h->recon_size_prim), BGEN64 (h->num_lines_bound), BGEN32 (h->num_sections), h->num_txt_files,
                 BGEN16(h->vb_size), 
                 SEC_TAB, sections_dis_flags (f, st, dt).s,
                 DT(REF) ? "fasta" : "ref", REF_FILENAME_LEN, h->ref_filename, 
                 DT(REF) ? "refhash_digest" : "ref_genome_digest", 
                 DT(REF) ? digest_display(h->refhash_digest).s : digest_display_ex (h->ref_genome_digest, DD_MD5).s,
                 dt_specific, 
                 SEC_TAB, FILE_METADATA_LEN, h->created);
        break;
    }

    case SEC_TXT_HEADER: {
        SectionHeaderTxtHeaderP h = (SectionHeaderTxtHeaderP)header;
        if (VER(15))
            sprintf (str, "\n%stxt_data_size=%"PRIu64" txt_header_size=%"PRIu64" lines=%"PRIu64" max_lines_per_vb=%u digest=%s digest_header=%s\n" 
                    "%ssrc_codec=%s (args=0x%02X.%02X.%02X) %s txt_filename=\"%.*s\"\n",
                    SEC_TAB, BGEN64 (h->txt_data_size), v12 ? BGEN64 (h->txt_header_size) : 0, BGEN64 (h->txt_num_lines), BGEN32 (h->max_lines_per_vb), 
                    digest_display (h->digest).s, digest_display (h->digest_header).s, 
                    SEC_TAB, codec_name (h->src_codec), h->codec_info[0], h->codec_info[1], h->codec_info[2], 
                    sections_dis_flags (f, st, dt).s, TXT_FILENAME_LEN, h->txt_filename);
        else
            sprintf (str, "\n%stxt_data_size=%"PRIu64" txt_header_size=%"PRIu64" lines=%"PRIu64" max_lines_per_vb=%u digest=%s digest_header=%s\n" 
                    "%ssrc_codec=%s (args=0x%02X.%02X.%02X) %s txt_filename=\"%.*s\" flav_prop=(id,has_seq_len,is_mated,cnn)=[[%u,%u,%u],[%u,%u,%u],[%u,%u,%u]]\n",
                    SEC_TAB, BGEN64 (h->txt_data_size), v12 ? BGEN64 (h->txt_header_size) : 0, BGEN64 (h->txt_num_lines), BGEN32 (h->max_lines_per_vb), 
                    digest_display (h->digest).s, digest_display (h->digest_header).s, 
                    SEC_TAB, codec_name (h->src_codec), h->codec_info[0], h->codec_info[1], h->codec_info[2], 
                    sections_dis_flags (f, st, dt).s, TXT_FILENAME_LEN, h->txt_filename,
                    h->flav_prop[0].has_seq_len, h->flav_prop[0].is_mated, h->flav_prop[0].cnn,
                    h->flav_prop[1].has_seq_len, h->flav_prop[1].is_mated, h->flav_prop[1].cnn,
                    h->flav_prop[2].has_seq_len, h->flav_prop[2].is_mated, h->flav_prop[2].cnn);

        break;
    }

    case SEC_VB_HEADER: {
        SectionHeaderVbHeaderP h = (SectionHeaderVbHeaderP)header;
        if (Z_DT(VCF) || Z_DT(BCF)) 
            sprintf (str, 
                    "\n%srecon_size=(PRIM:%u, LUFT:%u) longest_line=%u z_data_bytes=%u digest=%s %s\n",
                    SEC_TAB, BGEN32 (h->recon_size_prim), v12 ? BGEN32 (h->dvcf_recon_size_luft) : 0, 
                    BGEN32 (h->longest_line_len), 
                    BGEN32 (h->z_data_bytes), digest_display (h->digest).s, sections_dis_flags (f, st, dt).s);
        else if (Z_DT(SAM))
            sprintf (str, 
                    "\n%srecon_size=%u longest_line=%u longest_seq=%u z_data_bytes=%u digest=%s prim=(seq=%u comp_qual=%u qname=%u num_alns=%u first_grp_i=%u %s=%u) %s\n", 
                    SEC_TAB, BGEN32 (h->recon_size_prim),  
                    BGEN32 (h->longest_line_len), BGEN32(h->longest_seq_len),
                    BGEN32 (h->z_data_bytes), digest_display (h->digest).s, 
                    v14 ? BGEN32 (h->sam_prim_seq_len)          : 0,
                    v14 ? BGEN32 (h->sam_prim_comp_qual_len)    : 0,
                    v14 ? BGEN32 (h->sam_prim_qname_len)        : 0,
                    v14 ? BGEN32 (h->sam_prim_num_sag_alns)     : 0,
                    v14 ? BGEN32 (h->sam_prim_first_grp_i)      : 0,
                    IS_SAG_SA?"comp_cigars" : IS_SAG_SOLO?"solo_data" : "unused",
                    v14 ? BGEN32 (h->sam_prim_comp_cigars_len)  : 0,
                    sections_dis_flags (f, st, dt).s);
        else if (Z_DT(FASTQ))
            sprintf (str, 
                    "\n%srecon_size=%u longest_line=%u longest_seq=%u z_data_bytes=%u digest=%s %s\n",
                    SEC_TAB, BGEN32 (h->recon_size_prim),BGEN32 (h->longest_line_len), BGEN32(h->longest_seq_len),
                    BGEN32 (h->z_data_bytes), digest_display (h->digest).s, sections_dis_flags (f, st, dt).s);
        else
            sprintf (str, 
                    "\n%srecon_size=%u longest_line=%u z_data_bytes=%u digest=%s %s\n",
                    SEC_TAB, BGEN32 (h->recon_size_prim),BGEN32 (h->longest_line_len), 
                    BGEN32 (h->z_data_bytes), digest_display (h->digest).s, sections_dis_flags (f, st, dt).s);

        break;
    }

    case SEC_REFERENCE:
    case SEC_REF_IS_SET: {
        SectionHeaderReferenceP h = (SectionHeaderReferenceP)header;
        sprintf (str, "pos=%-9"PRIu64" gpos=%-9"PRIu64" num_bases=%-6u chrom_word_index=%-4d\n",
                 BGEN64 (h->pos), BGEN64 (h->gpos), BGEN32 (h->num_bases), BGEN32 (h->chrom_word_index)); 
        break;
    }
    
    case SEC_REF_HASH: {
        SectionHeaderRefHashP h = (SectionHeaderRefHashP)header;
        sprintf (str, "num_layers=%u layer_i=%u layer_bits=%u start_in_layer=%u\n",
                 h->num_layers, h->layer_i, h->layer_bits, BGEN32 (h->start_in_layer)); 
        break;
    }
    
    case SEC_REF_CONTIGS: {
        sprintf (str, "sequential_ref_index=%u\n", header->flags.ref_contigs.sequential_ref_index);
        break;
    }
    
    case SEC_RECON_PLAN: {
        SectionHeaderReconPlanP h = (SectionHeaderReconPlanP)header;
        sprintf (str, "conc_writing_vbs=%u %s\n", BGEN32 (h->conc_writing_vbs), sections_dis_flags (f, st, dt).s); 
        break;
    }
        
    case SEC_BGZF:
    case SEC_RANDOM_ACCESS: {
        sprintf (str, "%s%s\n", SEC_TAB, sections_dis_flags (f, st, dt).s); 
        break;
    }
    
    case SEC_B250: {
        SectionHeaderCtxP h = (SectionHeaderCtxP)header;
        sprintf (str, "%s/%-8s\tb250_size=%u param=%u %s\n",
                 dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s,  
                 h->b250_size==B250_BYTES_1?1 : h->b250_size==B250_BYTES_2?2 : h->b250_size==B250_BYTES_3?3 : 4,
                 h->param, sections_dis_flags (f, st, dt).s);
        break;
    }

    case SEC_LOCAL: {
        SectionHeaderCtxP h = (SectionHeaderCtxP)header;
        sprintf (str, "%s/%-8s\tltype=%s param=%u %s\n",
                 dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, lt_name (h->ltype), h->param, sections_dis_flags (f, st, dt).s);
        break;
    }

    case SEC_DICT: {
        SectionHeaderDictionaryP h = (SectionHeaderDictionaryP)header;
        sprintf (str, "%s/%-8s\tnum_snips=%u %s\n", dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, BGEN32 (h->num_snips), 
                 sections_dis_flags (f, st, dt).s); 
        break;
    }

    case SEC_COUNTS: {
        SectionHeaderCountsP h = (SectionHeaderCountsP)header;
        sprintf (str, "  %s/%-8s param=%"PRId64"\t\n", dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, h->nodes_param); 
        break;
    }

    case SEC_SUBDICTS: {
        SectionHeaderSubDictsP h = (SectionHeaderSubDictsP)header;
        sprintf (str, "  %s/%-8s param=%"PRId64"\t\n", dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, h->param); 
        break;
    }

    default: 
        str[0] = '\n'; str[1] = 0; 
    }

    // if not going directly to the terminal, replace non-final newlines with \t, to allow grep etc
    if (!isatty(1)) str_replace_letter (str, strlen(str)-1, '\n', '\t');

    PRINT;
}

// called from main()
void noreturn genocat_show_headers (rom z_filename)
{
    z_file = file_open_z_read (z_filename);    
    SectionHeaderGenozipHeader header;

    flag.genocat_no_ref_file = true;

    TEMP_FLAG (show_headers, 0);
    TEMP_FLAG (genocat_no_reconstruct, 1);
    zfile_read_genozip_header (&header, HARD_FAIL); // also sets z_file->genozip_version
    RESTORE_FLAG (show_headers);
    RESTORE_FLAG (genocat_no_reconstruct);

    // normal --show-headers - go by the section list
    if (!flag.force) {
        ARRAY (SectionEnt, sec, z_file->section_list_buf);
        for (uint32_t i=0; i < sec_len; i++) {
            header = zfile_read_section_header (evb, &sec[i], SEC_NONE).genozip_header; // we assign the largest of the SectionHeader* types
            sections_show_header ((SectionHeaderP)&header, NULL, sec[i].offset, 'R');
        }        
    }

    // --show-headers --force - search for actual headers in case of file corruption
    else {
        file_seek (z_file, 0, SET, READ, HARD_FAIL);
        uint64_t gap; // gap before section
        uint64_t accumulated_gap = 0;
        SectionEntModifiable sec = { .st = SEC_GENOZIP_HEADER/*largest header*/ };

        for (int sec_i=0; zfile_advance_to_next_header (&sec.offset, &gap); sec_i++) {
            if (gap || accumulated_gap) 
                iprintf ("ERROR: unexpected of %"PRIu64" bytes before next section\n", gap + accumulated_gap);
            
            header = zfile_read_section_header (evb, &sec, SEC_NONE).genozip_header; // we assign the largest of the SectionHeader* types
            if (header.section_type < 0 || header.section_type >= NUM_SEC_TYPES) { // not true section - magic matches by chance
                sec.offset += 4;
                accumulated_gap += 4;
                continue;
            }
            
            int header_size = st_header_size (header.section_type);
            if (header.data_encrypted_len) header_size = ROUNDUP16(header_size);

            // up to v14 we verify v14_compressed_offset
            if (!VER(15) && header_size != BGEN32 (header.v14_compressed_offset)) {
                sec.offset += 4;
                accumulated_gap += 4;
                continue;
            }

            if (flag.show_headers == SHOW_ALL_HEADERS || flag.show_headers-1 == header.section_type) 
                iprintf ("%5u ", sec_i);
            
            sections_show_header ((SectionHeaderP)&header, NULL, sec.offset, 'R');

            sec.offset += header_size + BGEN32 (header.data_compressed_len);
            accumulated_gap = 0;
        }

        if (flag.show_headers == SHOW_ALL_HEADERS) {
            if (gap) iprintf ("ERROR: unexpected gap of %"PRIu64" bytes before Footer\n", gap);

            if ((sec.offset = zfile_read_genozip_header_get_offset (true)))
                iprintf ("R %9"PRIu64" FOOTER              genozip_header_offset=%"PRIu64"\n", 
                        z_file->disk_size - sizeof (SectionFooterGenozipHeader), sec.offset);
            else
                iprint0 ("ERROR: no valid Footer\n");
        }
    }

    exit_ok;
}

void sections_show_section_list (DataType dt) // optional - take data from z_data
{    
    for_buf (SectionEnt, s, z_file->section_list_buf)
        if (s->st == SEC_B250 || s->st == SEC_LOCAL || s->st == SEC_DICT)
            iprintf ("%5u %-20.20s %s%s%-8.8s\tvb=%s/%-4u offset=%-8"PRIu64"  size=%-6u  %s\n", 
                     BNUM(z_file->section_list_buf, s), st_name(s->st), 
                     s->dict_id.num ? dtype_name_z(s->dict_id) :"     ", 
                     s->dict_id.num ? "/" : "", 
                     s->dict_id.num ? dis_dict_id (s->dict_id).s : "", 
                     comp_name_ex (s->comp_i, s->st).s, s->vblock_i, s->offset, s->size, 
                     sections_dis_flags (s->flags, s->st, dt).s);
        
        else if (s->st == SEC_VB_HEADER)
            iprintf ("%5u %-20.20s\t\t\tvb=%s/%-4u offset=%-8"PRIu64"  size=%-6u  num_lines=%-8u%s\n", 
                     BNUM(z_file->section_list_buf, s), st_name(s->st), comp_name (s->comp_i), 
                     s->vblock_i, s->offset, s->size, s->num_lines, sections_dis_flags (s->flags, s->st, dt).s);

        else if (IS_FRAG_SEC(s->st) || s->st == SEC_BGZF)
            iprintf ("%5u %-20.20s\t\t\tvb=%s/%-4u offset=%-8"PRIu64"  size=%-6u  %s\n", 
                     BNUM(z_file->section_list_buf, s), st_name(s->st), 
                     comp_name_ex (s->comp_i, s->st).s, s->vblock_i, s->offset, s->size, sections_dis_flags (s->flags, s->st, dt).s);

        else if (s->st == SEC_COUNTS || s->st == SEC_SUBDICTS)
            iprintf ("%5u %-20.20s %s%s%-8.8s\t             offset=%-8"PRIu64"  size=%-6u  %s\n", 
                     BNUM(z_file->section_list_buf, s), st_name(s->st),
                     s->dict_id.num ? dtype_name_z(s->dict_id) :"     ", 
                     s->dict_id.num ? "/" : "", 
                     s->dict_id.num ? dis_dict_id (s->dict_id).s : "", 
                     s->offset, s->size, sections_dis_flags (s->flags, s->st, dt).s);
        else
            iprintf ("%5u %-20.20s\t\t\t             offset=%-8"PRIu64"  size=%-6u  %s\n", 
                     BNUM(z_file->section_list_buf, s), st_name(s->st),
                     s->offset, s->size, sections_dis_flags (s->flags, s->st, dt).s);
}

void sections_show_gheader (ConstSectionHeaderGenozipHeaderP header)
{
    bool v15 = (IS_ZIP || VER(15));

    if (flag_loading_auxiliary) return; // don't show gheaders of an auxiliary file
    
    DataType dt = BGEN16 (header->data_type);
    ASSERT (dt < NUM_DATATYPES, "Invalid data_type=%u", dt);
    
    if (header) {
        iprintf ("Contents of the SEC_GENOZIP_HEADER section (output of --show-gheader) of %s:\n", z_name);
        iprintf ("  genozip_version: %u.0.%u\n",    header->genozip_version, header->genozip_minor_ver); // note: minor version always 0 before 15.0.28
        iprintf ("  data_type: %s\n",               dt_name (dt));
        iprintf ("  encryption_type: %s\n",         encryption_name (header->encryption_type)); 
        iprintf ("  recon_size_prim: %s\n",         str_int_commas (BGEN64 (header->recon_size_prim)).s);
        iprintf ("  num_lines_bound: %"PRIu64"\n",  BGEN64 (header->num_lines_bound));
        iprintf ("  num_sections: %u\n",            z_file->section_list_buf.len32);
        iprintf ("  num_txt_files: %u\n",           header->num_txt_files);
        if (dt == DT_REF)
            iprintf ("  %s: %s\n", (v15 ? "genome_digest" : "REF_fasta_md5"), digest_display (header->genome_digest).s);
        iprintf ("  created: %*s\n",                -FILE_METADATA_LEN, header->created);
        iprintf ("  license_hash: %s\n",            digest_display (header->license_hash).s);
        iprintf ("  private_file: %s\n",            TF(header->private_file));
        if (header->ref_filename[0]) {
            iprintf ("  reference filename: %s\n",  header->ref_filename);
            iprintf ("  reference file hash: %s\n", digest_display (header->ref_genome_digest).s);
        }
        iprintf ("  flags: %s\n",                   sections_dis_flags (header->flags, SEC_GENOZIP_HEADER, dt).s);

        switch (dt) {
            case DT_CHAIN:
                iprintf ("  primary-coordinates reference filename: %s\n", header->chain.prim_filename);
                iprintf ("  primary-coordinates reference file hash: %s\n", digest_display (header->chain.prim_genome_digest).s);
                break;

            default: break;
        }
    }

    iprint0 ("  sections:\n");

    sections_show_section_list (dt);
}
