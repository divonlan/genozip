// ------------------------------------------------------------------
//   sections.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "sections.h"
#include "file.h"
#include "zfile.h"
#include "threads.h"
#include "codec.h"

typedef struct SectionEnt SectionEntModifiable;

const struct {rom name; uint32_t header_size; } abouts[NUM_SEC_TYPES] = {
    [SEC_RANDOM_ACCESS]   = {"SEC_RANDOM_ACCESS",   sizeof (SectionHeader)              }, \
    [SEC_REFERENCE]       = {"SEC_REFERENCE",       sizeof (SectionHeaderReference)     }, \
    [SEC_REF_IS_SET]      = {"SEC_REF_IS_SET",      sizeof (SectionHeaderReference)     }, \
    [SEC_REF_HASH]        = {"SEC_REF_HASH",        sizeof (SectionHeaderRefHash)       }, \
    [SEC_REF_RAND_ACC]    = {"SEC_REF_RAND_ACC",    sizeof (SectionHeader)              }, \
    [SEC_REF_CONTIGS]     = {"SEC_REF_CONTIGS",     sizeof (SectionHeader)              }, \
    [SEC_GENOZIP_HEADER]  = {"SEC_GENOZIP_HEADER",  sizeof (SectionHeaderGenozipHeader) }, \
    [SEC_DICT_ID_ALIASES] = {"SEC_DICT_ID_ALIASES", sizeof (SectionHeader)              }, \
    [SEC_TXT_HEADER]      = {"SEC_TXT_HEADER",      sizeof (SectionHeaderTxtHeader)     }, \
    [SEC_VB_HEADER]       = {"SEC_VB_HEADER",       sizeof (SectionHeaderVbHeader)      }, \
    [SEC_DICT]            = {"SEC_DICT",            sizeof (SectionHeaderDictionary)    }, \
    [SEC_B250]            = {"SEC_B250",            sizeof (SectionHeaderCtx)           }, \
    [SEC_LOCAL]           = {"SEC_LOCAL",           sizeof (SectionHeaderCtx)           }, \
    [SEC_CHROM2REF_MAP]   = {"SEC_CHROM2REF_MAP",   sizeof (SectionHeader)              }, \
    [SEC_STATS]           = {"SEC_STATS",           sizeof (SectionHeader)              }, \
    [SEC_GZ_ISIZES]       = {"SEC_GZ_ISIZES",       sizeof (SectionHeader)              }, \
    [SEC_GZ_DIGESTS]      = {"SEC_GZ_DIGESTS",      sizeof (SectionHeaderGzDigests)     }, \
    [SEC_RECON_PLAN]      = {"SEC_RECON_PLAN",      sizeof (SectionHeaderReconPlan)     }, \
    [SEC_COUNTS]          = {"SEC_COUNTS",          sizeof (SectionHeaderCounts)        }, \
    [SEC_REF_IUPACS]      = {"SEC_REF_IUPACS",      sizeof (SectionHeader)              }, \
    [SEC_SUBDICTS]        = {"SEC_SUBDICTS",        sizeof (SectionHeaderSubDicts)      }, \
    [SEC_USER_MESSAGE]    = {"SEC_USER_MESSAGE",    sizeof (SectionHeader)              }, \
    [SEC_GENCOMP]         = {"SEC_GENCOMP",         sizeof (SectionHeader)              }, \
    [SEC_HUFFMAN]         = {"SEC_HUFFMAN",         sizeof (SectionHeaderHuffman)       }, \
};

const LocalTypeDesc lt_desc[NUM_LOCAL_TYPES] = LOCALTYPE_DESC;

typedef struct SectionsVbIndexEnt {
    int32_t vb_header_sec_i, last_sec_i; // -1 if none
    VBIType next_vb_i; // linked list of VBs of the same comp in the order of vb_i. list is terminated with -1. VBs are in the order they appear in section list, not necessarily consecutive vb_i's.
} SectionsVbIndexEnt;

typedef struct {
    int32_t txt_header_sec_i, gz_isizes_sec_i, recon_plan_sec_i;
    VBIType first_vb_i, last_vb_i, num_vbs; // vb_i's for each component form a linked list first_vb_i->next_vb_i->last_vb_i - vb_i's in each list are monotonously increasing, but not necessarily consecutive
} SectionsCompIndexEnt;

typedef struct {
    int32_t first_sec_i, last_sec_i;
} SectionsFirstLastIndexEnt;

static Section Bsec (int32_t sec_i)
{
    if (sec_i == -1)
        return NULL;

    else if (sec_i < z_file->section_list.len32)
        return B(SectionEnt, z_file->section_list, sec_i);

    else
        ABORT ("sec_i=%d ∉ [0,%u]", sec_i, z_file->section_list.len32-1);
}

static const SectionsVbIndexEnt *Bvbindex (VBIType vb_i)
{
    ASSERTINRANGE(vb_i, 1, z_file->vb_sections_index.len32);

    const SectionsVbIndexEnt *vb_index_ent = B(SectionsVbIndexEnt, z_file->vb_sections_index, vb_i);
    return vb_index_ent;
}

static const SectionsCompIndexEnt *Bcompindex (CompIType comp_i)
{
    ASSERT (comp_i < z_file->comp_sections_index.len32, "comp_i=%u ∉ [0,%u]", comp_i, z_file->comp_sections_index.len32-1);

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
        case SEC_HUFFMAN  : return ((SectionHeaderHuffmanP   )header)->dict_id; break;
        default           : return DICT_ID_NONE;
    }
}

// ZIP only: create section list that goes into the genozip header, as we are creating the sections. returns offset
void sections_add_to_list (VBlockP vb, ConstSectionHeaderP header)
{
    // case: we're re-creaating a section already on this list - nothing to do
    if (vb->section_list.len && vb->z_data.len <= BLST (SectionEnt, vb->section_list)->offset)
        return;

    SectionType st = header->section_type;
    ASSERT (st >= SEC_NONE && st < NUM_SEC_TYPES, "sec_type=%u ∉ [-1,%u]", st, NUM_SEC_TYPES-1);

    DictId dict_id = sections_get_dict_id (header);
    ASSERT0 (!IS_DICTED_SEC(st) || dict_id.num, "dict_id=0");
    
    buf_alloc (vb, &vb->section_list, 1, 50, SectionEnt, 2, "section_list");
    
    int header_size = st_header_size (header->section_type);
    if (header->data_encrypted_len) header_size = ROUNDUP16 (header_size);

    BNXT (SectionEntModifiable, vb->section_list) = (SectionEntModifiable) {
        .st             = header->section_type,
        .vblock_i       = (IS_VB_SEC(st) || IS_FRAG_SEC(st)) ? BGEN32 (header->vblock_i) : 0, // big endian in header - convert back to native
        .comp_i         = IS_COMP_SEC(st) ? vb->comp_i : COMP_NONE,
        .offset         = vb->z_data.len,  // this is a partial offset (within vb) - we will correct it later
        .flags          = header->flags,
        .num_lines      = ST(VB_HEADER) ? vb->lines.len32 : 0, 
        .size           = header_size  + MAX_(BGEN32 (header->data_compressed_len), BGEN32 (header->data_encrypted_len))
        //up to v14: header_size  + BGEN32 (header->data_compressed_len)
    };

    if (IS_DICTED_SEC(st)) 
        BLST (SectionEntModifiable, vb->section_list)->dict_id = dict_id; // note: this is in union with num_lines
}

// ZIP (possibly compute thread): remove section from list when it is deleted from z_data by zfile_remove_ctx_group_from_z_data
void sections_remove_from_list (VBlockP vb, uint64_t offset, uint64_t len)
{
    SectionEntModifiable *sec;
    for (sec=BLST (SectionEntModifiable, vb->section_list); sec->offset > offset; sec--) 
        sec->offset -= len;

    ASSERT (sec->offset == offset, "cannot find section with offset=%"PRIu64" in vb=%s", offset, VB_NAME);

    buf_remove (vb->section_list, SectionEnt, BNUM (vb->section_list, sec), 1);
}

// Called by ZIP main thread. concatenates a vb or dictionary section list to the z_file section list - just before 
// writing those sections to the disk. we use the current disk position to update the offset
void sections_list_concat (BufferP section_list)
{
    ASSERTMAINTHREAD;

    ARRAY (SectionEntModifiable, vb_sec, *section_list);

    // update the offset
    for (uint32_t i=0; i < vb_sec_len; i++) 
        vb_sec[i].offset += z_file->disk_so_far;

    // copy all entries
    buf_append_buf (evb, &z_file->section_list, section_list, SectionEntModifiable, "z_file->section_list");

    section_list->len = 0;
}

// section iterator. returns true if a section of this type was found.
bool sections_next_sec3 (Section *sl_ent,   // optional in/out. if NULL - search entire list
                         SectionType st1, SectionType st2, SectionType st3) // check only next section, not entire remaining list
{
    ASSERTNOTEMPTY (z_file->section_list);

    ARRAY (SectionsFirstLastIndexEnt, fl, z_file->first_last_sec_index);

    Section sec = sl_ent ? *sl_ent : NULL; 
    bool found = false;

    if (!sec && fl_len) { // indexed
        uint32_t firsts[3] = { (st1 != SEC_NONE ? fl[st1].first_sec_i : 0xffffffff),
                               (st2 != SEC_NONE ? fl[st2].first_sec_i : 0xffffffff),
                               (st3 != SEC_NONE ? fl[st3].first_sec_i : 0xffffffff) };
        uint32_t first_sec_i = MIN_(MIN_(firsts[0], firsts[1]), firsts[2]);

        if (first_sec_i != 0xffffffff) {
            sec = B(SectionEnt, z_file->section_list, first_sec_i);
            found = true;
        }
        
        goto done;
    }

    SectionEnt *last_sec = BLST (SectionEnt, z_file->section_list);

    if (fl_len) { // note: if we have fl_len, we also definitely have last_sec_i, otherwise we also don't have first_sec_i and we already returned above 
        int32_t lasts[3] = { (st1 != SEC_NONE ? fl[st1].last_sec_i : 0),
                             (st2 != SEC_NONE ? fl[st2].last_sec_i : 0),
                             (st3 != SEC_NONE ? fl[st3].last_sec_i : 0) };
        int32_t last_sec_i = MAX_(MAX_(lasts[0], lasts[1]), lasts[2]);
    
        last_sec = B(SectionEnt, z_file->section_list, last_sec_i);
    }

    ASSERT (!sec || BISVALID (z_file->section_list, sec), "Invalid sec: st1=%s st2=%s st3=%s", 
            st_name (st1), st_name (st2), st_name (st3));

    while (sec < last_sec) {
        sec = sec ? (sec + 1) : B1ST (SectionEnt, z_file->section_list); 

        if (sec->st == st1 || sec->st == st2 || sec->st == st3) {
            found = true;
            break;
        }
    }

done:
    if (sl_ent) *sl_ent = sec;
    return found;
}

// section iterator. returns true if a section of this type was found.
bool sections_prev_sec2 (Section *sl_ent,   // optional in/out. if NULL - search entire list
                         SectionType st1, SectionType st2)
{
    Section sec = sl_ent ? *sl_ent : NULL; 

    ASSERT (!sec || IN_RANGX (sec, B1ST(SectionEnt, z_file->section_list), BLST(SectionEnt, z_file->section_list)),
           "Invalid sec: st1=%s st2=%s", st_name (st1), st_name (st2));

    while (!sec || sec >= B1ST (SectionEnt, z_file->section_list)) {

        sec = sec ? (sec - 1) : BLST (SectionEnt, z_file->section_list); 

        if (sec->st == st1 || sec->st == st2) {
            if (sl_ent) *sl_ent = sec;
            return true;
        }
    }

    return false;
}

// count how many sections we have of a certain type
uint32_t sections_count_sections_until (BufferP sec_list, SectionType st, Section first_sec, SectionType until_encountering)
{
    ASSERT0 (!first_sec || IN_BUF (first_sec, *sec_list), "Invalid sec");

    ARRAY (SectionEnt, secs, *sec_list);
    uint32_t count=0;
    uint32_t first_i = first_sec ? BNUM (*sec_list, first_sec) : 0;

    for (uint32_t i=first_i; i < secs_len; i++) 
        if (secs[i].st == st)
            count++;
        else if (secs[i].st == until_encountering)
            break;

    return count;
}

uint32_t sections_txt_header_get_num_fragments (void)
{
    uint32_t num_frags = 0;

    for (CompIType comp_i=0; comp_i < z_file->num_components; comp_i++) 
        for (Section sec = sections_get_comp_txt_header_sec (comp_i); sec && IS_TXT_HEADER(sec) && sec->comp_i == comp_i; sec++)
            num_frags++;

    return num_frags;
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

    // free existing index in case it is being re-created (after each txt_file is compressed)
    if (IS_ZIP) {
        buf_free (z_file->comp_sections_index); 
        buf_free (z_file->vb_sections_index);
    }

    ASSERTNOTZERO (z_file->num_components);
    ARRAY_alloc (SectionsVbIndexEnt, vb_index, z_file->num_vbs + 1, true, z_file->vb_sections_index, evb, "z_file->vb_sections_index"); // +1 bc vb_i is 1-based
    ARRAY_alloc (SectionsCompIndexEnt, comp_index, z_file->num_components, true, z_file->comp_sections_index, evb, "z_file->comp_sections_index");
    ARRAY_alloc (SectionsFirstLastIndexEnt, fl_sec_index, NUM_SEC_TYPES, false, z_file->first_last_sec_index, evb, "z_file->first_last_sec_index");

    // initialize
    for (int i=0; i < vb_index_len; i++)
        vb_index[i].vb_header_sec_i = vb_index[i].last_sec_i = -1;

    for (int i=0; i < comp_index_len; i++)
        comp_index[i].gz_isizes_sec_i = comp_index[i].txt_header_sec_i = 
        comp_index[i].recon_plan_sec_i = -1;

    memset (fl_sec_index, 0xff, fl_sec_index_len * sizeof (SectionsFirstLastIndexEnt)); // set to -1

    // populate index
    SectionsVbIndexEnt *vbinx = NULL; // if non-NULL, we're within a VB section block
    for_buf2 (SectionEnt, sec, sec_i, z_file->section_list) {
        if (fl_sec_index[sec->st].first_sec_i == -1)
            fl_sec_index[sec->st].first_sec_i = sec_i;

        fl_sec_index[sec->st].last_sec_i = sec_i;

        if (IS_B250(sec) || IS_LOCAL(sec) || Z_DT(REF)) continue; // short circuit

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
            
            case SEC_GZ_ISIZES : 
                comp->gz_isizes_sec_i = sec_i; 
                break;

            case SEC_RECON_PLAN : {
                if (comp->recon_plan_sec_i == -1) // save first fragment
                    comp->recon_plan_sec_i = sec_i; 
                break;
            }

            // note: in ZIP, VBs are written out-of-order, so some VBs might still not be written and hence their index entry will remain empty
            case SEC_VB_HEADER : {
                VBIType vb_i = sec->vblock_i;
                
                ASSERTINRANGE (vb_i, 1, z_file->vb_sections_index.len32);
                vbinx = &vb_index[vb_i];
                *vbinx = (SectionsVbIndexEnt){ .vb_header_sec_i = sec_i };
                break;
            }

            default: break;
        }
    }

    // create linked lists of VBs of each component in order of vb_i
    // note: VBs may appear in arbitrary order in the file, and in the case of gencomp, VBs don't necessarily appear in the file in the order of components
    for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {
        CompIType comp_i = Bsec(vb_index[vb_i].vb_header_sec_i)->comp_i;
        SectionsCompIndexEnt *comp = &comp_index[comp_i];

        if (!comp->first_vb_i) 
            comp->first_vb_i = comp->last_vb_i = vb_i; // first vb_i in numerical order of vb_i
        
        else {
            vb_index[comp->last_vb_i].next_vb_i = vb_i;
            comp->last_vb_i = vb_i;
        }
        
        comp->num_vbs++;
    }
    
    COPY_TIMER_EVB (sections_create_index);
}

// ZIP / PIZ: returns the SEC_VB_HEADER section, or NULL if this vb_i doesn't exist
Section sections_vb_header (VBIType vb_i)
{    
    ASSERTINRANGX (vb_i, 1, z_file->num_vbs);

    uint32_t vb_header_sec_i = Bvbindex(vb_i)->vb_header_sec_i;
    Section sec = Bsec(vb_header_sec_i);
    
    // sanity
    ASSERT (sec, "vb_i=%u fits in z_file->vb_sections_index, but it is not indexed", vb_i);

    if (!IS_VB_HEADER(sec)) {
        sections_show_section_list (z_file->data_type, &z_file->section_list, SEC_VB_HEADER);
        ABORT ("Expecting indexed VB section of vb_i=%u sec_i=%u to be a VB_HEADER but it is a %s.", vb_i, vb_header_sec_i, st_name(sec->st));
    }

    return sec;
}

Section sections_vb_last_section (VBIType vb_i)
{    
    uint32_t last_sec_i = Bvbindex(vb_i)->last_sec_i;
    Section sec = Bsec(last_sec_i);
    
    // sanity
    ASSERT (sec, "vb_i=%u fits in z_file->vb_sections_index, but it is not indexed", vb_i);
 
    return sec;
}

//---------------------------------------------
// PIZ: constructing a updated section list
//---------------------------------------------

void sections_reading_list_add_vb_header (VBIType vb_i)
{
    BNXT(SectionEntModifiable, z_file->piz_reading_list) = *Bsec(Bvbindex(vb_i)->vb_header_sec_i);
}

void sections_reading_list_add_txt_header (CompIType comp_i)
{
    // add all fragments of the txt_header of this component
    for (Section sec = sections_get_comp_txt_header_sec (comp_i); 
         IS_TXT_HEADER(sec) && sec->comp_i == comp_i; 
         sec++)
        BNXT (SectionEntModifiable, z_file->piz_reading_list) = *sec;
}

//---------------
// misc functions
//---------------

// up to v14
static void v14_sections_list_file_to_memory_format (void)
{
    ARRAY (const SectionEntFileFormatV14, file_sec, z_file->section_list); // file format

    for_buf2_back (SectionEntModifiable, sec, i, z_file->section_list) {
        SectionEntFileFormatV14 fsec = file_sec[i]; // make a copy as assignment is overlapping
        
        *sec = (SectionEntModifiable){
            .offset      = BGEN64 (fsec.offset),
            .st_specific = fsec.st_specific,  
            .vblock_i    = BGEN32 (fsec.vblock_i),
            .st          = fsec.st,
            .comp_i      = (VER(14) && IS_COMP_SEC(fsec.st)) ? fsec.comp_i : COMP_NONE, // note: this field was introduced in v14: COMP_NONE sections have 0, as comp_i is just 2 bit. 
            .flags       = VER(12) ? fsec.flags  : (SectionFlags){} // prior to v12, flags were stored only in SectionHeader. Since v12, they are also copied to SectionEntFileFormat.
        };

        if (IS_VB_HEADER(sec)) {
            sec->num_lines  = BGEN32 (sec->num_lines);        
            z_file->num_vbs = MAX_(z_file->num_vbs, sec->vblock_i);

            if (sec->comp_i != COMP_NONE) {
                z_file->comp_num_lines[sec->comp_i] += sec->num_lines;
                z_file->num_components = MAX_(z_file->num_components, sec->comp_i);
            }
        }    

        else if (IS_TXT_HEADER(sec) && sec->comp_i != COMP_NONE)
            z_file->num_components = MAX_(z_file->num_components, sec->comp_i);

        if (i < file_sec_len-1)
            sec->size = (sec+1)->offset - sec->offset;
    }

    // comp_i was introduced V14, for V<=13, get comp_i by the component's relative position in the file. 
    if (VER(14)) 
        z_file->num_components++; // one more than the max comp_i

    else {
        CompIType comp_i_by_consecutive = COMP_NONE; 

        for_buf (SectionEntModifiable, sec, z_file->section_list) {
            if (IS_TXT_HEADER(sec)) comp_i_by_consecutive++; // note: we only read bound V<=13 files if the are a single FQ pair, or a DVCF

            if (IS_COMP_SEC(sec->st))
                sec->comp_i = comp_i_by_consecutive; 

            // Up to V11, flags were stored only in SectionHeader
            if (!VER(12)) 
                sec->flags = zfile_read_section_header (evb, sec, SEC_NONE).common.flags;
        }

        z_file->num_components = comp_i_by_consecutive + 1; // one more than the max comp_i
    }
}

// ZIP: creates file-format section list in evb->scratch (see also bug 819)
void sections_list_memory_to_file_format (void) 
{
    ASSERTNOTINUSE(evb->scratch);
    buf_alloc_exact (evb, evb->scratch, (MAX_(z_file->section_list.len * sizeof (SectionEntFileFormat), CODEC_ASSIGN_SAMPLE_SIZE) / sizeof(SectionEntFileFormat)),
                        SectionEntFileFormat, "scratch");
    evb->scratch.len = z_file->section_list.len;
    
    // replace dict_id with the the sec_i of its first appearahce. 
    int32_t first_appearance[z_file->ca.num_contexts];
    memset (first_appearance, 255, z_file->ca.num_contexts * sizeof (int32_t));

    SectionEntModifiable prev_sec = {};
    uint32_t prev_num_lines = 0;
    for_buf2 (SectionEntFileFormat, fsec, i, evb->scratch) {
        SectionEnt sec = *B(SectionEnt, z_file->section_list, i); 
        
        int64_t offset_delta = (int64_t)sec.offset - (int64_t)prev_sec.offset;
        ASSERT (IN_RANGX (offset_delta, 0LL, 0xffffffffLL),  // note: offset_delta is size of previous section
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

        else if (IS_VB_HEADER(&sec)) {
            int32_t num_lines_delta = INTERLACE(int32_t, (int32_t)sec.num_lines - (int32_t)prev_num_lines);
            prev_num_lines = sec.num_lines;

            fsec->num_lines = BGEN32 (num_lines_delta);
        }

        prev_sec = sec;
    }

    evb->scratch.len *= sizeof (SectionEntFileFormat); // change to counting bytes
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
            (dt == DT_FASTQ && f.v14_is_paired && v13_num_components == 2) || // Paired FASTQs (the v14_is_paired flag was introduced in V9.0.13)
            is_genols                       ||
            flags_is_genocat_global_area_only(), // only show meta-data (can't use flag.genocat_global_area_only bc flags are not set yet)
            "%s is comprised of %u %s files bound together. The bound file feature was discontinued in Genozip v14. To decompress this file, use Genozip v13",
            z_name, v13_num_components, dt_name (BGEN16(genozip_header->data_type)));
    
    z_file->section_list.len /= VER(15) ? sizeof (SectionEntFileFormat) : sizeof (SectionEntFileFormatV14); // fix len

    buf_alloc (evb, &z_file->section_list, 0, z_file->section_list.len, SectionEnt, 0, NULL); // extend

    if (VER(15)) {
        ARRAY (const SectionEntFileFormat, file_sec, z_file->section_list);  // file format

        // note: we work backwards as mem_sec items are larger than file_sec
        for_buf2_back (SectionEntModifiable, sec, i, z_file->section_list) { // memory format (entries are larger)
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
        for_buf2 (SectionEntModifiable, sec, i, z_file->section_list) {
            if (i) (sec-1)->size = sec->offset;
            if (i) sec->offset += (sec-1)->offset;

            sec->comp_i = (sec->comp_i == 0)         ? (sec-1)->comp_i
                        : (sec->comp_i == COMP_NONE) ? COMP_NONE
                        :                              (sec->comp_i - 1);

            sec->vblock_i = i ? (int32_t)(sec-1)->vblock_i + (int32_t)sec->vblock_i : (int32_t)sec->vblock_i;

            if (IS_DICTED_SEC(sec->st) && !sec->is_dict_id) 
                sec->dict_id = (sec - i + BGEN32(sec->dict_sec_i))->dict_id; // copy from section of first appearance of this dict_id

            else if (IS_VB_HEADER(sec)) {
                sec->num_lines  = prev_num_lines + DEINTERLACE(int32_t, (int32_t)BGEN32 (sec->num_lines/*this is still delta_num_lines*/));
                prev_num_lines  = sec->num_lines;
                z_file->num_vbs = MAX_(z_file->num_vbs, sec->vblock_i);

                z_file->comp_num_lines[sec->comp_i] += sec->num_lines;
            }    

            if (IS_VB_HEADER(sec) || IS_TXT_HEADER(sec))
                z_file->num_components = MAX_(z_file->num_components, sec->comp_i);
        }
    
        z_file->num_components++; // one more than the max comp_i
    }

    else  // up to v14
        v14_sections_list_file_to_memory_format();

    BLST(SectionEntModifiable, z_file->section_list)->size =
        st_header_size (SEC_GENOZIP_HEADER) + BGEN32 (genozip_header->data_compressed_len) + sizeof (SectionFooterGenozipHeader);

    sections_create_index(); // PIZ indexing
    
    // Up to V11, top_level_repeats existed (called num_lines) but was not populated
    // in V12-13 num_lines was transmitted through the VbHeader.top_level_repeats
    // Since V14, num_lines is transmitted through SectionEntFileFormat
    if (VER(12) && !VER(14)) 
        for (VBIType vb_i=1; vb_i <= z_file->num_vbs; vb_i++) {
            SectionEntModifiable *sec = (SectionEntModifiable *)sections_vb_header (vb_i);
            SectionHeaderUnion header = zfile_read_section_header (evb, sec, SEC_VB_HEADER);
            sec->num_lines = BGEN32 (header.vb_header.v13_top_level_repeats);
        }
}

// check if there is any DICT, LOCAL, B250 or COUNT section of a certain dict_id
bool is_there_any_section_with_dict_id (DictId dict_id)
{
    for_buf (SectionEnt, sec, z_file->section_list)
        if (sec->dict_id.num == dict_id.num && IS_DICTED_SEC (sec->st))
            return true;

    return false;
}

// get a section with a specific vb_i, section type and dict_id, or return NULL
SectionEnt *section_get_section (VBIType vb_i, SectionType st, DictId dict_id)
{
    const SectionsVbIndexEnt *vbent = Bvbindex(vb_i);
    ARRAY (SectionEnt, sec, z_file->section_list);

    for (uint32_t i=vbent->vb_header_sec_i+1; i <= vbent->last_sec_i; i++) 
        if (sec[i].dict_id.num == dict_id.num && sec[i].st == st)
            return &sec[i];

    return NULL; // not found
}

uint32_t st_header_size (SectionType sec_type)
{
    ASSERTINRANGE (sec_type, SEC_NONE, NUM_SEC_TYPES);

    if (!IS_ZIP) 
        switch (sec_type) {
            case SEC_NONE : 
                return 0;
            
            case SEC_VB_HEADER:
                if      (!VER(12)) return sizeof (SectionHeaderVbHeader) - 5*sizeof(uint32_t);  // in v8-11, SectionHeaderVbHeader was shorter by 5 32b words
                else if (!VER(14)) return sizeof (SectionHeaderVbHeader) - 2*sizeof(uint32_t);  // in v12-13, SectionHeaderVbHeader was shorter by 2 32b word
                else if (!VER(15)) return sizeof (SectionHeaderVbHeader) - 1*sizeof(uint32_t);  // in v14, SectionHeaderVbHeader was shorter by 1 32b word
                break;

            case SEC_TXT_HEADER:
                if      (!VER(12)) return sizeof (SectionHeaderTxtHeader) - 3*sizeof(uint16_t) - sizeof(uint64_t); // in v8-11, SectionHeaderTxtHeader was shorter
                else if (!VER(15)) return sizeof (SectionHeaderTxtHeader) - 3*sizeof(uint16_t); // in v12-14, SectionHeaderTxtHeader was shorter by 3XQnameFlavorProp
                break;

            case SEC_GENOZIP_HEADER:
                if (!VER(12)) return sizeof (SectionHeaderGenozipHeader) - REF_FILENAME_LEN - sizeof(Digest); // in v8-11, SectionHeaderTxtHeader was shorter
                break;

            default: break;
        }

    return abouts[sec_type].header_size;
}

Section section_next (Section sec)
{
    if (!z_file->section_list.len) return NULL;
    if (sec == NULL) return B1ST (SectionEnt, z_file->section_list);
    if (sec < BLST (SectionEnt, z_file->section_list)) return sec+1;
    return NULL;
}

// called by PIZ main thread
Section sections_first_sec (SectionType st, FailType soft_fail)
{
    ASSERTINRANGE (st, 0, NUM_SEC_TYPES);

    int32_t first_sec_i = B(SectionsFirstLastIndexEnt, z_file->first_last_sec_index, st)->first_sec_i; 

    if (first_sec_i != -1)
        return B(SectionEnt, z_file->section_list, first_sec_i);
    
    ASSERT (soft_fail, "Cannot find section_type=%s in z_file", st_name (st));

    return NULL;
}

// called by ZIP, PIZ main thread
Section sections_last_sec (SectionType st, FailType soft_fail)
{
    ASSERTINRANGE (st, SEC_NONE, NUM_SEC_TYPES);
    ARRAY (SectionEnt, sec, z_file->section_list);

    // case SEC_NONE - get the very last section regardless of type
    if (st == SEC_NONE)
        return &sec[sec_len-1];

    // case: indexed (always true in PIZ)
    else if (z_file->first_last_sec_index.len) { 
        int32_t last_sec_i = B(SectionsFirstLastIndexEnt, z_file->first_last_sec_index, st)->last_sec_i; 

        if (last_sec_i != -1) return &sec[last_sec_i];
    }

    // case: not indexed (can only happen in ZIP)
    else 
        for (int i=sec_len-1; i >= 0; i--)
            if (sec[i].st == st) return &sec[i];

    ASSERT (soft_fail, "Cannot find section_type=%s in z_file", st_name (st));

    return NULL;
}

VBIType sections_get_num_vbs (CompIType comp_i) 
{ 
    if (comp_i >= z_file->comp_sections_index.len32) return 0; // this component doesn't exist (can happen eg when checking for gencomp that doesn't exist in this file)

    uint32_t num_vbs = Bcompindex(comp_i)->num_vbs;
    return num_vbs;
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
    // ZIP: works for current comp_i in certain cases, most importantly for a FQ component in deep
    if (comp_i > 0 && comp_i == z_file->comp_sections_index.len32 && IS_ZIP && flag.zip_comp_i == comp_i)
        return Bcompindex(comp_i-1)->first_vb_i + Bcompindex(comp_i-1)->num_vbs; 

    VBIType first_vb_i = Bcompindex(comp_i)->first_vb_i;
    ASSERT (first_vb_i, "comp_i=%s has no VBs", comp_name (comp_i));
    return first_vb_i;
}

// get number of VBs in comp_i, for which vblock_i is at most max_vb_i
VBIType sections_get_num_vbs_up_to (CompIType comp_i, VBIType max_vb_i)
{
    VBIType num_vbs = 0;

    for (VBIType vb_i = Bcompindex(comp_i)->first_vb_i; 
         vb_i != 0 && vb_i <= max_vb_i;
         vb_i = Bvbindex(vb_i)->next_vb_i) // 0 means end-of-list
        num_vbs++;

    return num_vbs;
}

Section sections_get_comp_txt_header_sec (CompIType comp_i)
{
    return Bsec(Bcompindex(comp_i)->txt_header_sec_i);
}

Section sections_get_comp_GZ_ISIZES_sec (CompIType comp_i)
{
    return Bsec(Bcompindex(comp_i)->gz_isizes_sec_i);
}

// note: only files v12-15.0.63 could have a recon sec
uint32_t sections_get_recon_plan (Section *recon_sec)
{
        // add all fragments of the txt_header of this component
    Section sec = Bsec(Bcompindex(SAM_COMP_MAIN)->recon_plan_sec_i);
    if (recon_sec) *recon_sec = sec;

    if (!sec) return 0; // no SEC_RECON_PLAN section

    uint32_t num_fragments = 0;
    while (sec->st == SEC_RECON_PLAN) sec++, num_fragments++;

    return num_fragments;
}

// get next VB_HEADER in the order they appear in the section list (not necessarily consecutive vb_i)
Section sections_get_next_vb_header_sec (CompIType comp_i, Section *vb_header_sec)
{
    const SectionsCompIndexEnt *comp_index_ent = Bcompindex(comp_i);

    if (!comp_index_ent->first_vb_i)
        *vb_header_sec = NULL; // component has no VBs
        
    else if (! *vb_header_sec) {
        const SectionsVbIndexEnt *first_vb_ent = Bvbindex(comp_index_ent->first_vb_i);
        *vb_header_sec = Bsec(first_vb_ent->vb_header_sec_i);
    }

    else {
        VBIType prev_vb_i = (*vb_header_sec)->vblock_i;
        if (prev_vb_i == comp_index_ent->last_vb_i)
            *vb_header_sec = NULL; // no more VBs in this component

        else {
            VBIType next_vb_i = Bvbindex(prev_vb_i)->next_vb_i;
            const SectionsVbIndexEnt *next_vb_ent = Bvbindex(next_vb_i);
            *vb_header_sec = Bsec (next_vb_ent->vb_header_sec_i);
        }
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
                    (Z_DT(SAM)   && z_file->comp_sections_index.len32 >= 5 && // >= bc possibly multiple pairs of FASTQs
                     Bsec(Bcompindex(SAM_COMP_FQ01)->txt_header_sec_i)->flags.txt_header.pair);

    // up to V14, only FASTQ files could be paired. back comp note: until v13, it was possible to have bound files with multiple pairs (is this true?). not supported here.
    else if (!Z_DT(FASTQ) || z_file->num_txt_files != 2)  
        is_paired = false;

    // from V10 to V14 - determine by v14_is_paired (introduced in 9.0.13) 
    else 
        is_paired = z_file->z_flags.v14_is_paired;
    
    // sanity check: R1 and R2 are expected to have the same number of VBs (only if index is already created)
    if (is_paired && z_file->comp_sections_index.len32) {
        uint32_t num_vbs_R1 = Bcompindex(z_file->comp_sections_index.len32-2)->num_vbs;
        uint32_t num_vbs_R2 = Bcompindex(z_file->comp_sections_index.len32-1)->num_vbs;
        ASSERT (num_vbs_R1 == num_vbs_R2, "R1 and R2 have a different number of VBs (%u vs %u respecitvely)", num_vbs_R1, num_vbs_R2); 
    }

    return is_paired;
}

