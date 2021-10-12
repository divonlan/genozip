// ------------------------------------------------------------------
//   sections.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

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

typedef struct SectionEnt SectionEntModifiable;

static const struct {const char *name; uint32_t header_size; } abouts[NUM_SEC_TYPES] = SECTIONTYPE_ABOUT;

const LocalTypeDesc lt_desc[NUM_LOCAL_TYPES] = LOCALTYPE_DESC;

// ZIP only: create section list that goes into the genozip header, as we are creating the sections. returns offset
void sections_add_to_list (VBlock *vb, const SectionHeader *header)
{
    DictId dict_id = DICT_ID_NONE;
    bool has_dict_id = true;

    if      (header->section_type == SEC_DICT  ) dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    else if (header->section_type == SEC_B250  ) dict_id = ((SectionHeaderCtx        *)header)->dict_id;
    else if (header->section_type == SEC_LOCAL ) dict_id = ((SectionHeaderCtx        *)header)->dict_id;
    else if (header->section_type == SEC_COUNTS) dict_id = ((SectionHeaderCounts     *)header)->dict_id;
    else has_dict_id = false;

    ASSERT0 (!has_dict_id || dict_id.num, "dict_id=0");
    
    buf_alloc (vb, &vb->section_list_buf, 1, 50, SectionEnt, 2, "section_list_buf");
    
    NEXTENT (SectionEntModifiable, vb->section_list_buf) = (SectionEntModifiable) {
        .st       = header->section_type,
        .vblock_i = BGEN32 (header->vblock_i), // big endian in header - convert back to native
        .dict_id  = dict_id,
        .offset   = vb->z_data.len,  // this is a partial offset (within d) - we will correct it later
        .flags    = header->flags
    };
}

// Called by ZIP main thread. concatenates a vb or dictionary section list to the z_file section list - just before 
// writing those sections to the disk. we use the current disk position to update the offset
void sections_list_concat (VBlock *vb)
{
    if (!vb->section_list_buf.len) return;

    buf_alloc (evb, &z_file->section_list_buf, vb->section_list_buf.len, 0, SectionEnt, 2, "z_file->section_list_buf");

    SectionEntModifiable *dst = AFTERENT (SectionEntModifiable, z_file->section_list_buf);
    SectionEntModifiable *src = FIRSTENT (SectionEntModifiable, vb->section_list_buf);

    // update the offset
    for (unsigned i=0; i < vb->section_list_buf.len; i++)
        src[i].offset += z_file->disk_so_far;

    // copy all entries
    memcpy (dst, src, vb->section_list_buf.len * sizeof(SectionEnt));
    z_file->section_list_buf.len += vb->section_list_buf.len;

    buf_free (&vb->section_list_buf);
}

// section iterator. returns true if a section of this type was found.
bool sections_next_sec2 (Section *sl_ent,   // optional in/out. if NULL - search entire list
                         SectionType st1, SectionType st2) // check only next section, not entire remaining list
{
    Section sl = sl_ent ? *sl_ent : NULL; 
    bool found = false;

    while (sl < AFTERENT (SectionEnt, z_file->section_list_buf) - 1) {

        sl = sl ? (sl + 1) : FIRSTENT (SectionEnt, z_file->section_list_buf); 

        if (sl->st == st1 || sl->st == st2) {
            found = true;
            break;
        }
    }

    if (sl_ent) *sl_ent = sl;
    return found;
}

// section iterator. returns true if a section of this type was found.
bool sections_prev_sec2 (Section *sl_ent,   // optional in/out. if NULL - search entire list
                         SectionType st1, SectionType st2)
{
    Section sl = sl_ent ? *sl_ent : NULL; 

    while (!sl || sl >= FIRSTENT (SectionEnt, z_file->section_list_buf)) {

        sl = sl ? (sl - 1) : LASTENT (SectionEnt, z_file->section_list_buf); 

        if (sl->st == st1 || sl->st == st2) {
            if (sl_ent) *sl_ent = sl;
            return true;
        }
    }

    return false;
}

// section iterator. skips to the last consecutive section from after sl, that is type st1, st2 or st3.
Section sections_last_sec4 (Section sl, SectionType st1, SectionType st2, SectionType st3, SectionType st4)
{
    while (sl < AFTERENT (SectionEnt, z_file->section_list_buf) - 1) {

        sl = sl ? (sl + 1) : FIRSTENT (SectionEnt, z_file->section_list_buf); 

        if (!((st1 != SEC_NONE && sl->st == st1) ||
              (st2 != SEC_NONE && sl->st == st2) ||
              (st3 != SEC_NONE && sl->st == st3) ||
              (st4 != SEC_NONE && sl->st == st4))) {
            break;
        }
    } 
    
    return sl - 1; // sl is either on the first section that is NOT st1 or st2, or we are after the section list - so we decrement it
}

// count how many sections we have of a certain type
uint32_t sections_count_sections (SectionType st)
{
    uint32_t count=0;
    for (uint32_t i=0; i < z_file->section_list_buf.len; i++) 
        if (ENT (SectionEnt, z_file->section_list_buf, i)->st == st)
            count++;

    return count;
}

// called by PIZ main thread : vcf_zfile_read_one_vb. returns the first section of this vb_i
Section sections_vb_first (uint32_t vb_i, bool soft_fail)
{
    Section sl=NULL;
    unsigned i=0; for (; i < z_file->section_list_buf.len; i++) {
        sl = ENT (SectionEnt, z_file->section_list_buf, i);
        if (sl->vblock_i == vb_i) break; // found!
    }

    if (i >= z_file->section_list_buf.len) {
        if (soft_fail) return NULL;
        ABORT ("Error in sections_get_next_vb_section: cannot find any section for vb_i=%u", vb_i);
    }

    return sl;
}

// main thread: called by refhash_initialize - get details of the refhash ahead of loading it from the reference file 
void sections_get_refhash_details (uint32_t *num_layers, uint32_t *base_layer_bits) // optional outs
{
    ASSERT0 (flag.reading_reference || flag.show_ref_hash, "can only be called while reading reference");

    for (int i=z_file->section_list_buf.len-1; i >= 0; i--) { // search backwards as the refhash sections are near the end
        Section sl = ENT (SectionEnt, z_file->section_list_buf, i);
        if (sl->st == SEC_REF_HASH) {

            SectionHeaderRefHash header = zfile_read_section_header (evb, sl->offset, 0, SEC_REF_HASH).ref_hash;
            if (num_layers) *num_layers = header.num_layers;
            if (base_layer_bits) *base_layer_bits = header.layer_bits + header.layer_i; // layer_i=0 is the base layer, layer_i=1 has 1 bit less etc
            return;
        }
        else if (sl->st == SEC_REFERENCE)
            break; // we arrived at a SEC_REFERENCE - there won't be any more SEC_REF_HASH sections
    }

    ABORT ("Error in sections_get_refhash_details: can't find SEC_REF_HASH sections in %s", z_name);
}

void BGEN_sections_list (void)
{
    ARRAY (SectionEntModifiable, ent, z_file->section_list_buf);

    for (unsigned i=0; i < ent_len; i++) {
        ent[i].vblock_i = BGEN32 (ent[i].vblock_i);
        ent[i].offset   = BGEN64 (ent[i].offset);
        if (command==PIZ && z_file->genozip_version < 12) // flags were introduced in v12
            ent[i].flags.flags = 0;
    }
}

const char *st_name (SectionType sec_type)
{
    ASSERT (sec_type >= SEC_NONE && sec_type < NUM_SEC_TYPES, "sec_type=%u out of range [-1,%u]", sec_type, NUM_SEC_TYPES-1);

    return (sec_type == SEC_NONE) ? "SEC_NONE" : type_name (sec_type, &abouts[sec_type].name , sizeof(abouts)/sizeof(abouts[0]));
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
         : sec_type == SEC_VB_HEADER      && command == PIZ && z_file->genozip_version < 12 ? sizeof (SectionHeaderVbHeader) - 3*sizeof(uint32_t) // in v8-11, SectionHeaderVbHeader was shorter
         : sec_type == SEC_TXT_HEADER     && command == PIZ && z_file->genozip_version < 12 ? sizeof (SectionHeaderTxtHeader) -  sizeof(uint64_t) // in v8-11, SectionHeaderTxtHeader was shorter
         : sec_type == SEC_GENOZIP_HEADER && command == PIZ && z_file->genozip_version < 12 ? sizeof (SectionHeaderGenozipHeader) -  REF_FILENAME_LEN - sizeof(Digest) // in v8-11, SectionHeaderTxtHeader was shorter
         : abouts[sec_type].header_size;
}

Section section_next (Section sec)
{
    if (!z_file->section_list_buf.len) return NULL;
    if (sec == NULL) return FIRSTENT (SectionEnt, z_file->section_list_buf);
    if (sec < LASTENT (SectionEnt, z_file->section_list_buf)) return sec+1;
    return NULL;
}

// called by PIZ main thread
Section sections_first_sec (SectionType st, bool soft_fail)
{
    ARRAY (SectionEnt, sl, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++)
        if (sl[i].st == st) return &sl[i];

    ASSERT (soft_fail, "Cannot find section_type=%s in z_file", st_name (st));

    return NULL;
}

// called by PIZ main thread
Section sections_last_sec (SectionType st, bool soft_fail)
{
    ARRAY (SectionEnt, sl, z_file->section_list_buf);

    for (int i=z_file->section_list_buf.len-1; i >= 0; i--)
        if (sl[i].st == st) return &sl[i];

    ASSERT (soft_fail, "Cannot find section_type=%s in z_file", st_name (st));

    return NULL;
}

// count VBs of this component
void sections_count_component_vbs (Section sl, // must be set to SEC_TXT_HEADER of the requested component
                                   uint32_t *first_vb, uint32_t *num_vbs) // out
{
    *first_vb = *num_vbs = 0;

    while (sections_next_sec2 (&sl, SEC_VB_HEADER, SEC_TXT_HEADER) &&
           sl->st == SEC_VB_HEADER) {
        (*num_vbs)++;
        if (! *first_vb) *first_vb = sl->vblock_i;
    }
}

// move all sections of vb_i to be immediately after sl ; returns last section of vb_i after move
// note: vb_i is expected to start strictly after sl
Section sections_pull_vb_up (uint32_t vb_i, Section sl) 
{
    Section new_vb_sl    = sl+1;
    Section vb_header_sl = sections_vb_first (vb_i, false);
    Section last_vb_sl   = sections_vb_last (vb_header_sl);

    // case: VB is already in its desired place
    if (sl+1 == vb_header_sl) return last_vb_sl; // return last section (B250/Local) of this VB

    // save all entries of the VB in temp
    uint32_t vb_sections_len  = (last_vb_sl - vb_header_sl) + 1;
    uint32_t vb_sections_size = vb_sections_len * sizeof (SectionEnt);
    SectionEntModifiable temp[vb_sections_len]; // size is bound by MAX_DICTS x 2 x sizeof (SectionEnt) = 2048x2x24=96K
    memcpy (temp, vb_header_sl, vb_sections_size);

    // push down all sections after sl and before vb_header_sl
    memmove ((SectionEntModifiable *)new_vb_sl + vb_sections_len, new_vb_sl, (vb_header_sl - new_vb_sl) * sizeof (SectionEnt));

    // copy our vb to its requested place at sl+1
    memcpy ((SectionEntModifiable *)new_vb_sl, temp, vb_sections_size);
             
    return new_vb_sl + vb_sections_len - 1; // last section of the VB
}

// move all sections of component txtfile_sl_move_me to be immediately after txtfile_sl_after_me, or txtfile_sl_after_me=NULL to move to the top of the section list
void sections_pull_component_up (Section txtfile_sl_after_me, // destination - after this component (NULL means move to beginning)
                                 Section txtfile_sl_move_me)  // this is the component to move
{    
    Section dst_sl = txtfile_sl_after_me ? sections_component_last (txtfile_sl_after_me) + 1
                                         : FIRSTENT (SectionEnt, z_file->section_list_buf);

    // case: component is already in its desired place
    if (dst_sl == txtfile_sl_move_me) return;

    // save all entries of the VB in temp
    Section last_moving_sl = sections_component_last (txtfile_sl_move_me);
    uint32_t vb_sections_len  = (last_moving_sl - txtfile_sl_move_me) + 1;
    uint32_t vb_sections_size = vb_sections_len * sizeof (SectionEnt);
    
    SectionEntModifiable *temp = MALLOC (vb_sections_size);
    
    // swap
    memcpy (temp, txtfile_sl_move_me, vb_sections_size); // moving component moves to temp
    memmove ((SectionEntModifiable *)dst_sl + vb_sections_len, dst_sl, (txtfile_sl_move_me - dst_sl) * sizeof (SectionEnt)); // push down all sections after dst_sl, making room for moving component
    memcpy ((SectionEntModifiable *)dst_sl, temp, vb_sections_size); // moving component moves to dst

    FREE (temp);
}

int64_t sections_get_section_size (Section sl)
{
    // last section is always the genozip header
    if (sl->st == SEC_GENOZIP_HEADER)
        return z_file->disk_size - sl->offset;

    // case: the sections are consecutive
    Section next_vb_sl;
    if ((sl->st != SEC_B250 && sl->st != SEC_LOCAL) ||              // this is not a B250/LOCAL section
        (sl+1)->st == SEC_B250 || (sl+1)->st == SEC_LOCAL ||        // next section is still a B250/LOCAL section
        ((sl+1)->st == SEC_VB_HEADER && (sl+1)->vblock_i == sl->vblock_i+1) || // next section is the consecutive VB header
        !(next_vb_sl = sections_vb_first (sl->vblock_i + 1, true))) // this is the last VB in the file
        return (sl+1)->offset - sl->offset;

    // at this point, we know sl is the last B250/LOCAL section of a VB (but not the last VB of the file) 
    // and the next VB doesn't follow immediately
    
    // case: next VB is in the same component (because it NOT the first in a new component)
    if ((next_vb_sl-1)->st != SEC_TXT_HEADER) 
        return (next_vb_sl->offset - sl->offset);

    // case: next VB is in a new component directly consecutive - only a TXT_HEADER or a RECON_PLAN+TXT_HEADER 
    // are between the two VBs - next section is consecutive
    if (next_vb_sl - sl <= 2) 
        return (sl+1)->offset - sl->offset;

    // case: this is the last VB of the component - and this component has a recon plan - the recon plan is the next section
    Section sl_after_b250_local = sl;
    if (sections_next_sec2 (&sl_after_b250_local, SEC_TXT_HEADER, SEC_RECON_PLAN) &&
        sl_after_b250_local->st == SEC_RECON_PLAN)
        return sl_after_b250_local->offset - sl->offset;

    // case: this is the last VB of the component - the next section is the TXT_HEADER of the next vb's component
    return (next_vb_sl-1)->offset - sl->offset;
}

int64_t sections_get_vb_size (Section sl)
{
    ASSERT0 (sl->st == SEC_VB_HEADER, "expecing sl to be of a VB header");

    Section last_sl_in_vb = sections_vb_last (sl);

    return last_sl_in_vb->offset - sl->offset +       // size of all B250/LOCAL sections except the last
           sections_get_section_size (last_sl_in_vb); // size of the last B250/LOCAL
}

int64_t sections_get_vb_skipped_sections_size (Section vb_header_sl)
{
    ASSERT0 (vb_header_sl->st == SEC_VB_HEADER, "expecing sl to be of a VB header");

    Section last_sl_in_vb = sections_vb_last (vb_header_sl);
    int64_t size=0;

    for (Section sl = vb_header_sl+1; sl <= last_sl_in_vb; sl++)
        if (piz_is_skip_sectionz (sl->st, sl->dict_id))
            size += (sl==last_sl_in_vb) ? sections_get_section_size (last_sl_in_vb)
                                        : (sl+1)->offset - sl->offset;

    ASSERT (size >= 0, "expecting size=%"PRId64 "to be >=0", size);

    return size;
}

// returns the total size SEC_REFERENCE, SEC_REF_IS_SET, SEC_REF_RAND_ACCESS and SEC_REF_CONTIGS occupy in z_file
int64_t sections_get_ref_size (void)
{
    Section first_sl = sections_first_sec (SEC_REFERENCE, true);
    if (!first_sl) return 0;

    Section sl = sections_last_sec4 (first_sl, SEC_REFERENCE, SEC_REF_IS_SET, SEC_NONE, SEC_NONE);
    
    int64_t size = (sl+1)->offset - first_sl->offset;

    if ((sl = sections_last_sec (SEC_REF_RAND_ACC, true))) // there is only one section of this type
        size += (sl+1)->offset - sl->offset;

    if ((sl = sections_last_sec (SEC_REF_CONTIGS, true))) // there is only one section of this type
        size += (sl+1)->offset - sl->offset;

    ASSERT (size >= 0, "expecting size=%"PRId64 "to be >=0", size);

    return size;
}

const char *lt_name (LocalType lt)
{
    if (lt >= 0 && lt < NUM_LOCAL_TYPES) 
        return lt_desc[lt].name;
    else
        return "INVALID_LT";
}

typedef struct { char s[128]; } FlagStr;
static FlagStr sections_dis_flags (SectionFlags f, SectionType st, DataType dt)
{
    static const char *store[4] = { [STORE_NONE]="NONE", [STORE_INT]="INT", [STORE_FLOAT]="FLOAT", [STORE_INDEX]="INDEX"};
    static const char *dts[NUM_DATATYPES] = { [DT_FASTQ]="dts_paired", [DT_SAM]="dts_ref_internal", [DT_CHAIN]="dts_mismatch" };

    FlagStr str;

    switch (st) {
        case SEC_GENOZIP_HEADER:
            sprintf (str.s, "%s=%u aligner=%u txt_is_bin=%u bgzf=%u adler=%u dual_coords=%u has_taxid=%u",
                     dts[dt] ? dts[dt] : "dt_specitic", f.genozip_header.dt_specific, f.genozip_header.aligner, 
                     f.genozip_header.txt_is_bin, f.genozip_header.bgzf, f.genozip_header.adler, f.genozip_header.dual_coords,
                     f.genozip_header.has_taxid);
            break;

        case SEC_TXT_HEADER:
            sprintf (str.s, "rejects_coord=%s is_txt_luft=%u", coords_name (f.txt_header.rejects_coord), f.txt_header.is_txt_luft);
            break;

        case SEC_VB_HEADER:
            sprintf (str.s, "coords=%s", coords_name (f.vb_header.coords));
            break;

        case SEC_BGZF:
            sprintf (str.s, "level=%u has_eof=%u", f.bgzf.level, f.bgzf.has_eof_block); 
            break;

        case SEC_LOCAL:
            sprintf (str.s, "store=%s,per_line=%u,delta=%u paired=%u cy_param=%u",
                     store[f.ctx.store], f.ctx.store_per_line, f.ctx.store_delta, f.ctx.paired, f.ctx.copy_local_param); // note: we don't print ctx_specific as its not currently used
            break;

        case SEC_B250:
            sprintf (str.s, "store=%s,per_line=%u,delta=%u paired=%u cy_param=%u all_same=%u",
                     store[f.ctx.store], f.ctx.store_per_line, f.ctx.store_delta, f.ctx.paired, f.ctx.copy_local_param, f.ctx.all_the_same); 
            break;
 
        case SEC_RANDOM_ACCESS:
        case SEC_RECON_PLAN:
            sprintf (str.s, "rejects_coord=%s is_txt_luft=%u", coords_name (f.txt_header.rejects_coord), f.txt_header.is_txt_luft);
            break;

        default: 
            str.s[0] = 0;
    }

    return str;
}

void sections_show_header (const SectionHeader *header, VBlock *vb /* optional if output to buffer */, uint64_t offset, char rw)
{
    if (flag_loading_auxiliary) return; // don't show headers of an auxiliary file
    
    if (  flag.show_headers   != -1 &&                 // we don't need to show all sections
          flag.show_headers-1 != header->section_type) // we don't need to show this section
        return;

    bool is_dict_offset = (header->section_type == SEC_DICT && rw == 'W'); // at the point calling this function in zip, SEC_DICT offsets are not finalized yet and are relative to the beginning of the dictionary area in the genozip file
    bool v12 = (command == ZIP || z_file->genozip_version >= 12);

    char str[2048];
    #define PRINT { if (vb) buf_add_string (vb, &vb->show_headers_buf, str); else iprintf ("%s", str); } 
    #define SEC_TAB "            ++  "

    sprintf (str, "%c %s%-*"PRIu64" %-19s %-4.4s %-4.4s vb=%-3u z_off=%-6u txt_len=%-7u z_len=%-7u enc_len=%-7u %s ",
             rw, 
             is_dict_offset ? "~" : "", 9-is_dict_offset, offset, 
             st_name(header->section_type), 
             codec_name (header->codec), codec_name (header->sub_codec),
             BGEN32 (header->vblock_i), 
             BGEN32 (header->compressed_offset), 
             BGEN32 (header->data_uncompressed_len), 
             BGEN32 (header->data_compressed_len), 
             BGEN32 (header->data_encrypted_len), 
             BGEN32 (header->magic) != GENOZIP_MAGIC ? "BAD-MAGIC" : ""); // usually, we don't want the magic to take up line real estatee
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
            sprintf (dt_specific, SEC_TAB "prim=\"%.*s\" md5=%s\n", REF_FILENAME_LEN, h->dt_specific.chain.prim_filename, 
                     digest_display_ex (h->dt_specific.chain.prim_file_md5, DD_MD5).s);

        sprintf (str, "\n"SEC_TAB "ver=%u enc=%s dt=%s usize=%"PRIu64" lines=%"PRIu64" secs=%u txts=%u digest_bound=%s\n" 
                      SEC_TAB "%s ref=\"%.*s\" md5ref=%s\n"
                      "%s" // dt_specific, if there is any
                      SEC_TAB "created=\"%.*s\"\n",
                 h->genozip_version, encryption_name (h->encryption_type), dt_name (dt), 
                 BGEN64 (h->recon_size_prim), BGEN64 (h->num_lines_bound), BGEN32 (h->num_sections), BGEN32 (h->num_components),
                 digest_display (h->digest_bound).s, sections_dis_flags (f, st, dt).s,
                 REF_FILENAME_LEN, h->ref_filename, digest_display_ex (h->ref_file_md5, DD_MD5).s,
                 dt_specific, 
                 FILE_METADATA_LEN, h->created);
        break;
    }

    case SEC_TXT_HEADER: {
        SectionHeaderTxtHeader *h = (SectionHeaderTxtHeader *)header;
        sprintf (str, "\n"SEC_TAB "txt_data_size=%"PRIu64" txt_header_size=%"PRIu64" lines=%"PRIu64" max_lines_per_vb=%u md5_single=%s digest_header=%s\n" 
                      SEC_TAB "txt_codec=%s (args=0x%02X.%02X.%02X) %s txt_filename=\"%.*s\"\n",
                 BGEN64 (h->txt_data_size), v12 ? BGEN64 (h->txt_header_size) : 0, BGEN64 (h->txt_num_lines), BGEN32 (h->max_lines_per_vb), 
                 digest_display (h->digest_single).s, digest_display (h->digest_header).s, 
                 codec_name (h->codec), h->codec_info[0], h->codec_info[1], h->codec_info[2], 
                 sections_dis_flags (f, st, dt).s, TXT_FILENAME_LEN, h->txt_filename);
        break;
    }

    case SEC_VB_HEADER: {
        SectionHeaderVbHeader *h = (SectionHeaderVbHeader *)header;
        sprintf (str, "\n"SEC_TAB "lines=(PRIM:%u, LUFT:%u) recon_size=(PRIM:%u, LUFT:%u) longest_line=%u z_data_bytes=%u digest_so_far=%s %s\n",
                 BGEN32 (h->num_lines_prim),  v12 ? BGEN32 (h->num_lines_luft)  : 0, 
                 BGEN32 (h->recon_size_prim), v12 ? BGEN32 (h->recon_size_luft) : 0, 
                 BGEN32 (h->longest_line_len), 
                 BGEN32 (h->z_data_bytes), digest_display (h->digest_so_far).s, sections_dis_flags (f, st, dt).s);
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
        sprintf (str, "%s/%-8s\tparam=%u %s\n",
                 dtype_name_z (h->dict_id), dis_dict_id (h->dict_id).s, h->param, sections_dis_flags (f, st, dt).s); 
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

void sections_show_gheader (const SectionHeaderGenozipHeader *header /* optional */)
{
    if (flag_loading_auxiliary) return; // don't show gheaders of an auxiliary file
    
    ARRAY (SectionEnt, ents, z_file->section_list_buf);

    DataType dt = BGEN16 (header->data_type);
    ASSERT (dt < NUM_DATATYPES, "Invalid data_type=%u", dt);
    
    if (header) {
        iprintf ("Contents of the SEC_GENOZIP_HEADER section (output of --show-gheader) of %s:\n", z_name);
        iprintf ("  genozip_version: %u\n",         header->genozip_version);
        iprintf ("  data_type: %s\n",               dt_name (dt));
        iprintf ("  encryption_type: %s\n",         encryption_name (header->encryption_type)); 
        iprintf ("  recon_size_prim: %s\n",         str_uint_commas (BGEN64 (header->recon_size_prim)).s);
        iprintf ("  num_lines_bound: %"PRIu64"\n",  BGEN64 (header->num_lines_bound));
        iprintf ("  num_sections: %u\n",            (unsigned)ents_len);
        iprintf ("  num_components: %u\n",          BGEN32 (header->num_components));
        iprintf ("  digest_bound.md5: %s\n",        digest_display (header->digest_bound).s);
        iprintf ("  created: %*s\n",                -FILE_METADATA_LEN, header->created);
        iprintf ("  license_hash: %s\n",            digest_display (header->license_hash).s);
        iprintf ("  reference filename: %s\n",      header->ref_filename);
        iprintf ("  reference file hash: %s\n",     digest_display (header->ref_file_md5).s);
        iprintf ("  flags: %s\n",                   sections_dis_flags (header->h.flags, SEC_GENOZIP_HEADER, dt).s);

        switch (dt) {
            case DT_CHAIN:
                iprintf ("  primary-coordinates reference filename: %s\n", header->dt_specific.chain.prim_filename);
                iprintf ("  primary-coordinates reference file hash: %s\n", digest_display (header->dt_specific.chain.prim_file_md5).s);
                break;

            default: break;
        }
    }

    iprint0 ("  sections:\n");

    for (unsigned i=0; i < ents_len; i++) {
     
        uint64_t this_offset = ents[i].offset;
        
        uint64_t next_offset = (i < ents_len-1) ? ents[i+1].offset
                             : header               ? this_offset + BGEN32 (header->h.data_compressed_len) + BGEN32 (header->h.compressed_offset) + sizeof (SectionFooterGenozipHeader) // we're at the last section genozip header+footer
                             :                        this_offset; // without the GenozipHeader, we can't know its length


        iprintf ("%3u\t%-20.20s\t%s%s%-8.8s\tvb_i=%u\toffset=%-8"PRIu64"\tsize=%-6"PRId64"\t%s\n", 
                 i, st_name(ents[i].st), 
                 ents[i].dict_id.num ? dtype_name_z(ents[i].dict_id) :"     ", 
                 ents[i].dict_id.num ? "/" : "", 
                 ents[i].dict_id.num ? dis_dict_id (ents[i].dict_id).s : "", 
                 ents[i].vblock_i, this_offset, next_offset - this_offset, sections_dis_flags (ents[i].flags, ents[i].st, dt).s);
    }
}
