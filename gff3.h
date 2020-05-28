// ------------------------------------------------------------------
//   gff3.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef GFF3_INCLUDED
#define GFF3_INCLUDED

#include "genozip.h"

// SEG Stuff
extern const char *gff3_seg_txt_line (VBlockP vb_, const char *field_start_line, bool *has_special_eol);
extern void gff3_seg_initialize (VBlockP vb_);

// PIZ Stuff
extern bool gff3_piz_read_one_vb (VBlockP vb, SectionListEntryP sl);
extern void gff3_piz_reconstruct_vb(); 

// VBlock stuff
extern void gff3_vb_release_vb();
extern void gff3_vb_destroy_vb();
extern unsigned gff3_vb_size (void);
extern unsigned gff3_vb_zip_dl_size (void);

#define GFF3_DICT_ID_ALIASES \
    /*          alias                           maps to this ctx          */  \
    { DT_GFF3, &dict_id_ATTR_Variant_seq  ,    &dict_id_ATTR_Reference_seq }, \
    { DT_GFF3, &dict_id_ATTR_ancestral_allele, &dict_id_ATTR_Reference_seq }, 

#define GFF3_LOCAL_COMPRESSOR_CALLBACKS

#define dict_id_is_gff3_attr_sf dict_id_is_type_1
#define dict_id_gff3_attr_sf dict_id_type_1

#endif

