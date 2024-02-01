// ------------------------------------------------------------------
//   data_types.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
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

void dt_initialize (void)
{
    // since V15, FASTQ_PREDEFINED is a copy of SAM PREDEFINED, except for FASTQ_CONTIG that remains "CONTIG" - see also ctx_initialize_predefined_ctxs
    dt_fields[DT_FASTQ].predefined[FASTQ_CONTIG].dict_id = (DictId)DICT_ID_MAKEF_6 ("CONTIG");
    dt_fields[DT_FASTQ].predefined[FASTQ_CONTIG].tag_name = "CONTIG";
    dt_fields[DT_FASTQ].predefined[FASTQ_CONTIG].tag_name_len = 6;
}

// PIZ: gets the toplevel and factor, and returns true if translating
const DtTranslation dt_get_translation (VBlockP vb) // vb=NULL relates to the txt header
{
    static DtTranslation translations[] = TRANSLATIONS;
    bool i_am_binary = z_file->z_flags.txt_is_bin;
    DataType z_dt = vb ? vb->data_type : z_file->data_type; // note: vb is NULL when called with txtheader 

    // exception: no SAM->FASTQ translation if deep
    for (unsigned i=0; i < ARRAY_LEN (translations); i++) {
        // if (vb->vblock_i == 1)
        //     printf ("translations[i].src_z_non_bin_dt=%u z_dt=%u " 
        //             "translations[i].src_z_is_binary=%u i_am_binary=%u "
        //             "translations[i].dst_txt_dt=%u  flag.out_dt=%u "
        //             "translations[i].is_translation=%u *translations[i].is_translation=%u\n", 
        //             translations[i].src_z_non_bin_dt , z_dt ,
        //             translations[i].src_z_is_binary  , i_am_binary ,
        //             translations[i].dst_txt_dt       , flag.out_dt ,
        //             translations[i].is_translation , *translations[i].is_translation); 

        if (translations[i].src_z_non_bin_dt == z_dt &&
            (translations[i].src_z_is_binary == i_am_binary || translations[i].src_z_is_binary == unknown) &&
            translations[i].dst_txt_dt       == flag.out_dt &&
            (!translations[i].is_translation || translations[i].is_translation (vb))) {// if non-NULL, this flag determines if it is a translation
            
            // special case BAM->SAM: since 15.0.28 segconf reports the actual factor which may vary widely between different files
            if (translations[i].toplevel.num == _SAM_TOPLEVEL && segconf.est_sam_factor) {
                DtTranslation copy = translations[i];
                copy.factor = segconf.est_sam_factor * 1.20; // 20% safety margin above what segconf reported
                return copy;
            }
            else
                return translations[i]; // get one of the following entries if extended
        }
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

    ABORT ("cannot find textual type for binary data type %s", dt_name (dt));
}

rom dt_name (DataType dt)
{
    if (dt == DT_NONE) return "NONE";
    return type_name (dt, &dt_props[dt].name, NUM_DATATYPES);
}

rom z_dt_name (void)
{
    return IS_SRC_BAM_PIZ ? "BAM" 
         : z_file         ? dt_name (z_file->data_type)
         :                  "ERR_NULL_Z_FILE";
}
