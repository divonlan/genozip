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
DataTypeProperties dt_props_def = DATA_TYPE_FUNCTIONS_DEFAULT;

DataTypeFields     dt_fields[NUM_DATATYPES] = DATA_TYPE_FIELDS;

static struct { DataType src_z_non_bin_dt; 
                uint8_t src_z_is_binary;  // same type as z_file.flags
                DataType dst_txt_dt;
                DidIType toplevel;
                float factor; 
                TXTHEADER_TRANSLATOR ((*txtheader_translator));
              } translations[] = TRANSLATIONS;

#define NUM_TRANSLATIONS (sizeof (translations) / sizeof (translations[0]))

// gets the toplevel and factor, and returns true if translating
bool dt_get_translation (DidIType *toplevel, float *factor) // optional outs
{
    DidIType my_toplevel;
    float my_factor;

    for (unsigned i=0; i < NUM_TRANSLATIONS; i++)
        if (translations[i].src_z_non_bin_dt == z_file->data_type &&
            translations[i].src_z_is_binary  == !!(z_file->flags & GENOZIP_FL_TXT_IS_BIN) &&
            translations[i].dst_txt_dt       == flag.out_dt) {
                my_factor   = translations[i].factor;
                my_toplevel = translations[i].toplevel;
                goto finish;
            }
            
    // not a translation
    my_factor   = 1;
    my_toplevel = DTFZ(toplevel);

finish:
    if (factor)   *factor   = my_factor;
    if (toplevel) *toplevel = my_toplevel;

    return (my_toplevel != DTFZ(toplevel)); // true we're using an alternative toplevel snip
}

// calls the txtheader translator
void dt_translate_txtheader (BufferP txt)
{
    for (unsigned i=0; i < NUM_TRANSLATIONS; i++)
        if (translations[i].src_z_non_bin_dt == z_file->data_type &&
            translations[i].src_z_is_binary  == !!(z_file->flags & GENOZIP_FL_TXT_IS_BIN) &&
            translations[i].dst_txt_dt       == flag.out_dt) {
                // we found a translation - it may or may not have a txtheader translator
                if (translations[i].txtheader_translator) translations[i].txtheader_translator (txt);
                break;
            }
}

const char *dt_name (DataType dt)
{
    return type_name (dt, &dt_props[dt].name, NUM_DATATYPES);
}

