// ------------------------------------------------------------------
//   data_types.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "data_types.h"
#include "vblock.h"
#include "seg.h"
#include "zip.h"
#include "zfile.h"
#include "piz.h"
#include "strings.h"

DataTypeProperties dt_props [NUM_DATATYPES] = DATA_TYPE_PROPERTIES;
DataTypeProperties dt_props_def = DATA_TYPE_FUNCTIONS_DEFAULT;

DataTypeFields     dt_fields[NUM_DATATYPES] = DATA_TYPE_FIELDS;

const char *dt_name (DataType dt)
{
    return type_name (dt, &dt_props[dt].name, NUM_DATATYPES);
}
