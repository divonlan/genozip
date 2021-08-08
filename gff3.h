// ------------------------------------------------------------------
//   gff3.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef GFF3_INCLUDED
#define GFF3_INCLUDED

#include "genozip.h"

#define DTYPE_GFF3_ATTR DTYPE_1
#define DTYPE_GFF3_ENST DTYPE_2

#define _GFF3_SEQID            DICT_ID_MAKEF_5 ("SEQID")
#define _GFF3_SOURCE           DICT_ID_MAKEF_6 ("SOURCE")
#define _GFF3_TYPE             DICT_ID_MAKEF_4 ("TYPE")
#define _GFF3_START            DICT_ID_MAKEF_5 ("START")
#define _GFF3_END              DICT_ID_MAKEF_3 ("END")
#define _GFF3_SCORE            DICT_ID_MAKEF_5 ("SCORE")
#define _GFF3_STRAND           DICT_ID_MAKEF_6 ("STRAND")
#define _GFF3_PHASE            DICT_ID_MAKEF_5 ("PHASE")
#define _GFF3_ATTRS            DICT_ID_MAKEF_5 ("ATTRS")
#define _GFF3_EOL              DICT_ID_MAKEF_3 ("EOL")
#define _GFF3_TOPLEVEL         DICT_ID_MAKEF_L (TOPLEVEL)
#define _GFF3_COMMENT          DICT_ID_MAKEF_7 ("COMMENT")

// standard GVF fields (ID is also a standard GFF3 field)
#define _ATTR_ID               DICT_ID_MAKE1_2 ("ID")
#define _ATTR_Variant_seq      DICT_ID_MAKE1_L ("Variant_seq")
#define _ATTR_Reference_seq    DICT_ID_MAKE1_L ("Reference_seq")
#define _ATTR_Variant_freq     DICT_ID_MAKE1_L ("Variant_freq")

// fields added in the GVFs of GRCh37/38
#define _ATTR_Dbxref           DICT_ID_MAKE1_6 ("Dbxref")
#define _ATTR_ancestral_allele DICT_ID_MAKE1_L ("ancestral_allele")
#define _ATTR_Variant_effect   DICT_ID_MAKE1_L ("Variant_effect")
#define _ATTR_sift_prediction  DICT_ID_MAKE1_L ("sift_prediction")
#define _ATTR_polyphen_prediction DICT_ID_MAKE1_L ("polyphen_prediction")
#define _ATTR_variant_peptide  DICT_ID_MAKE1_L ("variant_peptide")

#define _ENSTid                DICT_ID_MAKE2_6 ("ENSTid") 

typedef enum { GFF3_SEQID, GFF3_SOURCE, GFF3_TYPE, GFF3_START, GFF3_END, GFF3_SCORE, GFF3_STRAND, GFF3_PHASE, GFF3_ATTRS, GFF3_EOL, GFF3_TOPLEVEL, GFF3_COMMENT, 
               ATTR_ID, ATTR_Variant_effect, ATTR_sift_prediction, ATTR_polyphen_prediction, ATTR_variant_peptide, ATTR_Reference_seq,
               ATTR_Dbxref, ENSTid,
               NUM_GFF3_FIELDS } Gff3Fields;

#define GFF3_MAPPING { V(GFF3_SEQID), V(GFF3_SOURCE), V(GFF3_TYPE), V(GFF3_START), V(GFF3_END), V(GFF3_SCORE), \
                       V(GFF3_STRAND), V(GFF3_PHASE), V(GFF3_ATTRS), V(GFF3_EOL), V(GFF3_TOPLEVEL), V(GFF3_COMMENT), \
                       V(ATTR_ID), V(ATTR_Variant_effect), V(ATTR_sift_prediction), V(ATTR_polyphen_prediction), V(ATTR_variant_peptide), V(ATTR_Reference_seq), \
                       V(ATTR_Dbxref), V(ENSTid), }

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
    { DT_GFF3, _ATTR_Variant_seq  ,    _ATTR_Reference_seq }, \
    { DT_GFF3, _ATTR_ancestral_allele, _ATTR_Reference_seq }, 

#define GFF3_LOCAL_GET_LINE_CALLBACKS

#endif

