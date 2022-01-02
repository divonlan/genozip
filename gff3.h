// ------------------------------------------------------------------
//   gff3.h
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

#define DTYPE_GFF3_ATTR DTYPE_1
#define DTYPE_GFF3_ENST DTYPE_2

#pragma GENDICT_PREFIX GFF3

#pragma GENDICT GFF3_SEQID=DTYPE_FIELD=SEQID
#pragma GENDICT GFF3_SOURCE=DTYPE_FIELD=SOURCE
#pragma GENDICT GFF3_TYPE=DTYPE_FIELD=TYPE
#pragma GENDICT GFF3_START=DTYPE_FIELD=START
#pragma GENDICT GFF3_END=DTYPE_FIELD=END
#pragma GENDICT GFF3_SCORE=DTYPE_FIELD=SCORE
#pragma GENDICT GFF3_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT GFF3_PHASE=DTYPE_FIELD=PHASE
#pragma GENDICT GFF3_ATTRS=DTYPE_FIELD=ATTRS
#pragma GENDICT GFF3_EOL=DTYPE_FIELD=EOL
#pragma GENDICT GFF3_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT GFF3_COMMENT=DTYPE_FIELD=COMMENT

// standard GFF3 attributes defined in https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#pragma GENDICT ATTR_ID=DTYPE_1=ID
#pragma GENDICT ATTR_Name=DTYPE_1=Name
#pragma GENDICT ATTR_Alias=DTYPE_1=Alias
#pragma GENDICT ATTR_Parent=DTYPE_1=Parent
#pragma GENDICT ATTR_Target=DTYPE_1=Target
#pragma GENDICT ATTR_Target_ID=DTYPE_1=T1gtID
#pragma GENDICT ATTR_Target_POS=DTYPE_1=T2gtPOS
#pragma GENDICT ATTR_Target_STRAND=DTYPE_1=T3gtSTRAND
#pragma GENDICT ATTR_Gap=DTYPE_1=Gap
#pragma GENDICT ATTR_Derives_from=DTYPE_1=Derives_from
#pragma GENDICT ATTR_Note=DTYPE_1=Note
#pragma GENDICT ATTR_Dbxref=DTYPE_1=Dbxref
#pragma GENDICT ATTR_Ontology_term=DTYPE_1=Ontology_term
#pragma GENDICT ATTR_Is_circular=DTYPE_1=Is_circular

// standard GVF fields 
#pragma GENDICT ATTR_Variant_seq=DTYPE_1=Variant_seq
#pragma GENDICT ATTR_Reference_seq=DTYPE_1=Reference_seq
#pragma GENDICT ATTR_Variant_freq=DTYPE_1=Variant_freq

// fields added in the GVFs of GRCh37/38
#pragma GENDICT ATTR_ancestral_allele=DTYPE_1=ancestral_allele
#pragma GENDICT ATTR_Variant_effect=DTYPE_1=Variant_effect
#pragma GENDICT ATTR_sift_prediction=DTYPE_1=sift_prediction
#pragma GENDICT ATTR_polyphen_prediction=DTYPE_1=polyphen_prediction
#pragma GENDICT ATTR_variant_peptide=DTYPE_1=variant_peptide

// other fields
#pragma GENDICT ATTR_chr=DTYPE_1=chr

#pragma GENDICT ENSTid=DTYPE_2=ENSTid 
#pragma GENDICT EnNSTid=DTYPE_2=EnNSTid 

// SEG Stuff
extern void gff3_zip_initialize (void);
extern const char *gff3_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void gff3_seg_initialize (VBlockP vb_);
extern void gff3_seg_finalize (VBlockP vb);
extern bool gff3_seg_is_small (ConstVBlockP vb, DictId dict_id);

// VBlock stuff
extern void gff3_vb_release_vb();
extern void gff3_vb_destroy_vb();
extern unsigned gff3_vb_size (DataType dt);

// PIZ stuff
CONTAINER_FILTER_FUNC (gff3_piz_filter);

#define GFF3_DICT_ID_ALIASES \
    /*          alias                  maps to this ctx          */  \
    { DT_GFF3, _ATTR_Variant_seq  ,    _ATTR_Reference_seq }, \
    { DT_GFF3, _ATTR_ancestral_allele, _ATTR_Reference_seq }, 

#define GFF3_LOCAL_GET_LINE_CALLBACKS
