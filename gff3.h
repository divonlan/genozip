// ------------------------------------------------------------------
//   gff3.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef GFF3_INCLUDED
#define GFF3_INCLUDED

#include "genozip.h"

// SEG Stuff
extern const char *gff3_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void gff3_seg_initialize (VBlockP vb_);
extern void gff3_seg_finalize (VBlockP vb);

// VBlock stuff
extern void gff3_vb_release_vb();
extern void gff3_vb_destroy_vb();
extern unsigned gff3_vb_size (void);

#define GFF3_DICT_ID_ALIASES \
    /*          alias                           maps to this ctx          */  \
    { DT_GFF3, &dict_id_ATTR_Variant_seq  ,    &dict_id_ATTR_Reference_seq }, \
    { DT_GFF3, &dict_id_ATTR_ancestral_allele, &dict_id_ATTR_Reference_seq }, 

#define GFF3_LOCAL_GET_LINE_CALLBACKS

#define DTYPE_GFF3_ATTR    DTYPE_1
#define DTYPE_GFF3_ENST DTYPE_2
#endif

