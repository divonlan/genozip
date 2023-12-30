// ------------------------------------------------------------------
//   dict_id.c
//   Copyright (C) 2020-2024 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <stdarg.h>
#include "genozip.h"
#include "dict_id.h"
#include "file.h"
#include "zfile.h"

// note: this function is near-identical to the one generated by dict_id_gen.sh
DictId dict_id_make (STRp(str), DictIdType dict_id_type) 
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

    // bc it would causes dict_id.id[0] to be 0 and id[0]=0 means lookback in SectionEntFileFormat
    ASSERT (dict_id_type != DTYPE_FIELD || dict_id.id[0] != '@', "A DTYPE_FIELD dict_id=%s cannot begin with '@'.", dis_dict_id (dict_id).s);

    switch (dict_id_type) {
        case DTYPE_FIELD  : dict_id.id[0] = dict_id.id[0] & 0x3f; break;
        case DTYPE_1      : dict_id.id[0] = dict_id.id[0] | 0xc0; break;
        case DTYPE_2      : break;
        default: ABORT ("invalid type %d", dict_id_type);
    }
    return dict_id;
}

#define s0(s,dict_id_type) ((dict_id_type==DTYPE_1)     ? (s[0] | 0xc0) \
                          : (dict_id_type==DTYPE_2)     ?  s[0] \
                          :           /* DTYPE_FIELD */   (s[0] & 0x3f))

#define DICT_ID_MAKE(s, dict_id_type) ((sizeof s <= 2) ? {.id = { s0(s,dict_id_type), } } \
                                    :  (sizeof s == 3) ? {.id = { s0(s,dict_id_type), s[1] } } \
                                    :  (sizeof s == 4) ? {.id = { s0(s,dict_id_type), s[1], s[2], } } \
                                    :  (sizeof s == 5) ? {.id = { s0(s,dict_id_type), s[1], s[2], s[3] } } \
                                    :  (sizeof s == 6) ? {.id = { s0(s,dict_id_type), s[1], s[2], s[3], s[4] } } \
                                    :  (sizeof s == 7) ? {.id = { s0(s,dict_id_type), s[1], s[2], s[3], s[4], s[5] } } \
                                    :  (sizeof s == 8) ? {.id = { s0(s,dict_id_type), s[1], s[2], s[3], s[4], s[5], s[6] } } \
                                    :                    {.id = { s0(s,dict_id_type), s[1], s[2], s[3], s[4], s[5], s[6], s[7] } } )

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

rom dict_id_display_type (DataType dt, DictId dict_id)
{
    if (dict_id_is_field  (dict_id)) return dt_props[dt].dtype_names[0];
    if (dict_id_is_type_1 (dict_id)) return dt_props[dt].dtype_names[1]; 
    if (dict_id_is_type_2 (dict_id)) return dt_props[dt].dtype_names[2]; 

    ABORT0 ("");
}

// return true is dict_id is contained in the DICT_ID_NONE-terminated list
bool dict_id_is_in (DictId dict_id, ...)
{
    va_list args;                                       
    va_start (args, dict_id);

    uint64_t test_dnum; 
    do {
        test_dnum = (uint64_t)va_arg (args, uint64_t);
    } while (test_dnum && test_dnum != dict_id.num);

    va_end (args);                              

    return test_dnum != 0;
}

// print the dict_id - NOT thread safe, for use in execution-termination messages
DisplayPrintId dis_dict_id (DictId dict_id)
{
    if (!dict_id.num) return (DisplayPrintId){ .s = "(none)" };

    DisplayPrintId s = { .s[0] = (dict_id.id[0] & 0x7f) | 0x40 }; // set 2 Msb to 01

    for (int i=1; i < DICT_ID_LEN; i++) 
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

// show a dict_id if the requested string is a subset of it, excluding unprintable characters
bool dict_id_is_show (DictId dict_id)
{
    if (!flag.show_one_dict) return false;
    
    char dict_id_str[9] = "";
    unsigned s_len=0;
    dict_id = dict_id_typeless (dict_id);

    for (unsigned i=0; i < DICT_ID_LEN; i++)
        if (IS_NON_WS_PRINTABLE(dict_id.id[i]))
            dict_id_str[s_len++] = dict_id.id[i];

    unsigned len = strlen (flag.show_one_dict);
    if (len <= 8)
        return !dict_id_str[MIN_(len,8)] && !memcmp (dict_id_str, flag.show_one_dict, len);
    else
        return !dict_id_str[MIN_(len,8)] && !memcmp (dict_id_str, flag.show_one_dict, 4) && !memcmp(&dict_id_str[4], &flag.show_one_dict[len-4], 4);
}

