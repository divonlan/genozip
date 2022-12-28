// ------------------------------------------------------------------
//   gff.h
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#define DTYPE_GFF_ATTR DTYPE_1
#define DTYPE_GFF_ENST DTYPE_2

#pragma GENDICT_PREFIX GFF

#pragma GENDICT GFF_SEQID=DTYPE_FIELD=SEQID
#pragma GENDICT GFF_SOURCE=DTYPE_FIELD=SOURCE
#pragma GENDICT GFF_TYPE=DTYPE_FIELD=TYPE
#pragma GENDICT GFF_START=DTYPE_FIELD=START
#pragma GENDICT GFF_END=DTYPE_FIELD=END
#pragma GENDICT GFF_SCORE=DTYPE_FIELD=SCORE
#pragma GENDICT GFF_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT GFF_PHASE=DTYPE_FIELD=PHASE
#pragma GENDICT GFF_ATTRS=DTYPE_FIELD=ATTRS
#pragma GENDICT GFF_EOL=DTYPE_FIELD=EOL
#pragma GENDICT GFF_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT GFF_COMMENT=DTYPE_FIELD=COMMENT
#pragma GENDICT GFF_DEBUG_LINES=DTYPE_FIELD=DBGLINES      // used by --debug-lines

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
#pragma GENDICT ATTR_db_xref=DTYPE_1=db_xref
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

// GTF fields
#pragma GENDICT ATTR_gene_id=DTYPE_1=gene_id
#pragma GENDICT ATTR_gene_name=DTYPE_1=gene_name
#pragma GENDICT ATTR_transcript_id=DTYPE_1=transcript_id
#pragma GENDICT ATTR_transcript_name=DTYPE_1=transcript_name
#pragma GENDICT ATTR_transcript_name_gene=DTYPE_1=t0ranscript_name
#pragma GENDICT ATTR_transcript_name_num=DTYPE_1=t1ranscript_name

#pragma GENDICT ATTR_protein_id=DTYPE_1=protein_id
#pragma GENDICT ATTR_ccds_id=DTYPE_1=ccds_id
#pragma GENDICT ATTR_exon_id=DTYPE_1=exon_id
#pragma GENDICT ATTR_exon_number=DTYPE_1=exon_number

// other fields
#pragma GENDICT ATTR_chr=DTYPE_1=chr

#pragma GENDICT ENSTid=DTYPE_2=ENSTid 

// ZIP Stuff
extern void gff_zip_initialize (void);
extern bool is_gff (STRp(header), bool *need_more);
extern bool gff_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);
extern rom gff_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void gff_seg_initialize (VBlockP vb_);
extern void gff_seg_finalize (VBlockP vb);
extern bool gff_seg_is_small (ConstVBlockP vb, DictId dict_id);

// VBlock stuff
extern void gff_vb_release_vb();
extern void gff_vb_destroy_vb();
extern unsigned gff_vb_size (DataType dt);

// PIZ stuff
extern CONTAINER_FILTER_FUNC (gff_piz_filter);
extern CONTAINER_CALLBACK (gff_piz_container_cb);

#define GFF_SPECIAL {  gff_piz_special_exon_number }
SPECIAL (GFF, 0,  exon_number,         gff_piz_special_exon_number);
#define NUM_GFF_SPECIAL 1

#define GFF_DICT_ID_ALIASES \
    /*         alias                  maps to this ctx   */  \
    { DT_GFF, _ATTR_Variant_seq,      _ATTR_Reference_seq }, \
    { DT_GFF, _ATTR_ancestral_allele, _ATTR_Reference_seq }, 