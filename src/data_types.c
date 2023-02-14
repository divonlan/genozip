// ------------------------------------------------------------------
//   data_types.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "data_types.h"
#include "vblock.h"
#include "context.h"
#include "dict_id.h"
#include "file.h"
#include "strings.h"

DataTypeProperties dt_props [NUM_DATATYPES] = DATA_TYPE_PROPERTIES;
DataTypeProperties dt_props_def             = DATA_TYPE_FUNCTIONS_DEFAULT;
DataTypeFields     dt_fields[NUM_DATATYPES] = DATA_TYPE_FIELDS;

// PIZ: gets the toplevel and factor, and returns true if translating
const DtTranslation dt_get_translation (VBlockP vb) // vb=NULL relates to the txt header
{
    static DtTranslation translations[] = TRANSLATIONS;
    bool i_am_binary = z_file->z_flags.txt_is_bin;
    DataType z_dt = vb ? vb->data_type : z_file->data_type; // note: vb is NULL when called with txtheader 

    // exception: no SAM->FASTQ translation if deep
    for (unsigned i=0; i < ARRAY_LEN (translations); i++) {
/*if (i==ARRAY_LEN (translations)-1) { printf ("translations[i].src_z_non_bin_dt=%u z_dt=%u " 
            "translations[i].src_z_is_binary=%u i_am_binary=%u "
            "translations[i].dst_txt_dt=%u  flag.out_dt=%u "
            "translations[i].is_translation=%u *translations[i].is_translation=%u\n", 
            translations[i].src_z_non_bin_dt , z_dt ,
            translations[i].src_z_is_binary  , i_am_binary ,
            translations[i].dst_txt_dt       , flag.out_dt ,
            translations[i].is_translation , *translations[i].is_translation); exit_ok();}
*/
        if (translations[i].src_z_non_bin_dt == z_dt &&
            translations[i].src_z_is_binary  == i_am_binary &&
            translations[i].dst_txt_dt       == flag.out_dt &&
            (!translations[i].is_translation || translations[i].is_translation (vb))) // if non-NULL, this flag determines if it is a translation
            return translations[i + flag.extended_translation]; // get next one if extended
    }

    // translation not found - return a non-translation "translation"
    return (DtTranslation){ .dst_txt_dt      = z_dt, 
                            .factor          = 1, 
                            .is_src_dt       = true, 
                            .src_z_is_binary = i_am_binary, 
                            .toplevel        = DTFZ(toplevel) };
}

// if file is_txt_binary - return the equivalent textual type, or just the type if not
DataType dt_get_txt_dt (DataType dt)
{
    if (!dt_props[dt].is_binary) return dt;

    for (DataType txt_dt=0; txt_dt < NUM_DATATYPES; txt_dt++)
        if (dt_props[txt_dt].bin_type == dt)
            return txt_dt;

    ABORT_R ("cannot find textual type for binary data type %s", dt_name (dt));
}

rom dt_name (DataType dt)
{
    if (dt == DT_NONE) return "NONE";
    return type_name (dt, &dt_props[dt].name, NUM_DATATYPES);
}

