// ------------------------------------------------------------------
//   data_types.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

// gets the toplevel and factor, and returns true if translating
const DtTranslation dt_get_translation (void)
{
    static DtTranslation translations[] = TRANSLATIONS;
    #define NUM_TRANSLATIONS (sizeof (translations) / sizeof (translations[0]))

    bool i_am_binary = z_file->z_flags.txt_is_bin;

    for (unsigned i=0; i < NUM_TRANSLATIONS; i++)
        if (translations[i].src_z_non_bin_dt == z_file->data_type &&
            translations[i].src_z_is_binary  == i_am_binary &&
            translations[i].dst_txt_dt       == flag.out_dt) 
            return translations[i];

    // translation not found - return a non-translation "translation"
    return (DtTranslation){ .dst_txt_dt      = z_file->data_type, 
                            .factor          = 1, 
                            .is_src_dt       = true, 
                            .src_z_is_binary = i_am_binary, 
                            .toplevel        = DTFZ(toplevel) };
}

const char *dt_name (DataType dt)
{
    return type_name (dt, &dt_props[dt].name, NUM_DATATYPES);
}

