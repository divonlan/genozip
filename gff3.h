// ------------------------------------------------------------------
//   gff3.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef GFF3_INCLUDED
#define GFF3_INCLUDED

#include "genozip.h"

// SEG Stuff
extern const char *gff3_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void gff3_seg_initialize (VBlockP vb_);
extern void gff3_seg_finalize (VBlockP vb);
extern bool gff3_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern bool gff3_seg_special_info_subfields (VBlockP vb, DictId dict_id, const char **this_value, unsigned *this_value_len);

// VBlock stuff
extern void gff3_vb_release_vb();
extern void gff3_vb_destroy_vb();
extern unsigned gff3_vb_size (DataType dt);

// PIZ stuff
CONTAINER_FILTER_FUNC (gff3_piz_filter);

#define GFF3_DICT_ID_ALIASES \
    /*          alias                           maps to this ctx          */  \
    { DT_GFF3, &dict_id_ATTR_Variant_seq  ,    &dict_id_ATTR_Reference_seq }, \
    { DT_GFF3, &dict_id_ATTR_ancestral_allele, &dict_id_ATTR_Reference_seq }, 

#define GFF3_LOCAL_GET_LINE_CALLBACKS

#define DTYPE_GFF3_ATTR    DTYPE_1
#define DTYPE_GFF3_ENST DTYPE_2
#endif

