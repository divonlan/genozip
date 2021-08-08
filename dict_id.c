// ------------------------------------------------------------------
//   dict_id.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "dict_id.h"
#include "data_types.h"
#include "file.h"
#include "zfile.h"
#include "sections.h"
#include "vblock.h"
#include "libdeflate/libdeflate.h"
#include "vcf.h"

// globals externed in dict_id.h and initialized in dict_id_initialize
static Buffer dict_id_aliases_buf  = EMPTY_BUFFER;
const DictIdAlias *dict_id_aliases = NULL;
uint32_t dict_id_num_aliases = 0;

DictId dict_id_make (const char *str, unsigned str_len, DictIdType dict_id_type) 
{ 
    DictId dict_id = DICT_ID_NONE; 

    if (!str_len) str_len = strlen (str);

    if (str_len <= DICT_ID_LEN) 
        memcpy (dict_id.id, str, str_len);
    
    else { 
        #define half1_len (DICT_ID_LEN/2)
        #define half2_len (DICT_ID_LEN - DICT_ID_LEN/2)

        memcpy (dict_id.id, str, half1_len); // take 1/2 from the start and 1/2 from then end (note: vcf_lo_is_INFO_AF_type depends on this)
        memcpy (dict_id.id + half1_len, str+str_len-half2_len, half2_len);
    }

    switch (dict_id_type) {
        case DTYPE_FIELD  : dict_id.id[0] = dict_id.id[0] & 0x3f; break;
        case DTYPE_1      : dict_id.id[0] = dict_id.id[0] | 0xc0; break;
        case DTYPE_2      : break;
        default: ABORT ("Error in dict_id_make: invalid type %d", dict_id_type);
    }
    return dict_id;
}

#define s0(s,dict_id_type) ((dict_id_type==DTYPE_1)     ? (s[0] | 0xc0) \
                          : (dict_id_type==DTYPE_2)     ?  s[0] \
                          :           /* DTYPE_FIELD */   (s[0] & 0x3f))

#define DICT_ID_MAKE(s, dict_id_type) ((sizeof s <= 2) ? (DictId){.id = { s0(s,dict_id_type), } } \
                                    :  (sizeof s == 3) ? (DictId){.id = { s0(s,dict_id_type), s[1] } } \
                                    :  (sizeof s == 4) ? (DictId){.id = { s0(s,dict_id_type), s[1], s[2], } } \
                                    :  (sizeof s == 5) ? (DictId){.id = { s0(s,dict_id_type), s[1], s[2], s[3] } } \
                                    :  (sizeof s == 6) ? (DictId){.id = { s0(s,dict_id_type), s[1], s[2], s[3], s[4] } } \
                                    :  (sizeof s == 7) ? (DictId){.id = { s0(s,dict_id_type), s[1], s[2], s[3], s[4], s[5] } } \
                                    :  (sizeof s == 8) ? (DictId){.id = { s0(s,dict_id_type), s[1], s[2], s[3], s[4], s[5], s[6] } } \
                                    :                    (DictId){.id = { s0(s,dict_id_type), s[1], s[2], s[3], s[4], s[5], s[6], s[7] } } )

static void dict_id_show_aliases (void)
{
    iprint0 ("Contents of SEC_DICT_ID_ALIASES section:\n");
    for (unsigned i=0; i < dict_id_num_aliases; i++) 
        iprintf ("alias=%s dst=%s\n", 
                 dis_dict_id (dict_id_aliases[i].alias).s, dis_dict_id (dict_id_aliases[i].dst).s);

    if (exe_type == EXE_GENOCAT) exit_ok();
}

// called by ZIP main thread for writing to global section
Buffer *dict_id_create_aliases_buf (void)
{
    static struct { DataType dt; uint64_t dict_id_alias; uint64_t dict_id_dst; } aliases_def[] = DICT_ID_ALIASES;

    // count
    dict_id_num_aliases = 0;
    for (unsigned i=0; i < sizeof(aliases_def)/sizeof(aliases_def[0]); i++)
        if (aliases_def[i].dt == z_file->data_type)
            dict_id_num_aliases++;

    // build global alias reference, which will be immutable until the end of this z_file
    dict_id_aliases_buf.len = dict_id_num_aliases * sizeof (DictIdAlias);
    buf_alloc (evb, &dict_id_aliases_buf, 0, dict_id_aliases_buf.len, char, 1, "dict_id_aliases_buf");

    DictIdAlias *next = FIRSTENT (DictIdAlias, dict_id_aliases_buf);
    for (unsigned i=0; i < sizeof(aliases_def)/sizeof(aliases_def[0]); i++)
        if (aliases_def[i].dt == z_file->data_type) {
            next->alias = (DictId)aliases_def[i].dict_id_alias;
            next->dst   = (DictId)aliases_def[i].dict_id_dst;
            next++;
        }

    if (flag.show_aliases) dict_id_show_aliases();

    return &dict_id_aliases_buf;
}

// PIZ main thread: read all dict_id aliaeses, if there are any
void dict_id_read_aliases (void) 
{ 
    Section sl = sections_last_sec (SEC_DICT_ID_ALIASES, true);
    if (!sl) return; // no aliases section

    buf_free (&dict_id_aliases_buf); // needed in case this is the 2nd+ file being pizzed

    zfile_get_global_section (SectionHeader, SEC_DICT_ID_ALIASES, sl, &dict_id_aliases_buf, "dict_id_aliases_buf");

    dict_id_aliases = FIRSTENT (DictIdAlias, dict_id_aliases_buf);
    dict_id_num_aliases = dict_id_aliases_buf.len / sizeof (DictIdAlias);

    for (unsigned i=0; i < dict_id_num_aliases; i++) 
        ASSERT0 (dict_id_aliases[i].dst.id[0] && dict_id_aliases[i].alias.id[0], "corrupted aliases buffer");
    
    if (flag.show_aliases) dict_id_show_aliases();
}

// template can be 0 - anything OR a type - must 2 MSb of id[0] are used OR a specific dict_id
// candidate is a specific dict_id that we test for its matching of the template
bool dict_id_is_match (DictId template, DictId candidate)
{
    if (!template.num) return true;

    if (template.num == candidate.num) return true;

    // template if it is 0 except for the 2 MSb of byte 0 - i.e. we're searching for a type of dict_id rather than a specific one
    DictId copy = template;
    copy.id[0] &= 0x3f; // remove the two MSb */

    if (!copy.num && ((template.id[0] & 0xc0) == (candidate.id[0] & 0xc0)))
        return true;

    return false;
}

const char *dict_id_display_type (DataType dt, DictId dict_id)
{
    if (dict_id_is_field  (dict_id)) return dt_props[dt].dtype_names[0];
    if (dict_id_is_type_1 (dict_id)) return dt_props[dt].dtype_names[1]; 
    if (dict_id_is_type_2 (dict_id)) return dt_props[dt].dtype_names[2]; 

    ABORT0_R ("Error in dict_id_display_type");
}

// print the dict_id - NOT thread safe, for use in execution-termination messages
DisplayPrintId dis_dict_id (DictId dict_id)
{
    DisplayPrintId s = {};

    if (!dict_id.num) return (DisplayPrintId){ .s = "<none>" };

    s.s[0] = (dict_id.id[0] & 0x7f) | 0x40;  // set 2 Msb to 01

    for (unsigned i=1; i < DICT_ID_LEN; i++) 
        s.s[i] = dict_id.id[i] == 0  ? ' '
               : dict_id.id[i] < 32  ? '?' // non printable char is converted to '?'
               : dict_id.id[i] > 126 ? '?'
               :                       dict_id.id[i];

    // trim final ' '
    for (int i=DICT_ID_LEN-1; i >= 0; i--) {
        if (s.s[i] != ' ') break;
        s.s[i] = 0;
    }

    return s;
}

