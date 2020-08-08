// ------------------------------------------------------------------
//   sam.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SAM_INCLUDED
#define SAM_INCLUDED

#include "genozip.h"
#include "sections.h"

// ZIP Stuff
COMPRESSOR_CALLBACK(sam_zip_get_start_len_line_i_seq)
COMPRESSOR_CALLBACK(sam_zip_get_start_len_line_i_qual)
COMPRESSOR_CALLBACK(sam_zip_get_start_len_line_i_bd)
COMPRESSOR_CALLBACK(sam_zip_get_start_len_line_i_bi)
extern void sam_zip_initialize (void);
extern bool sam_zip_is_unaligned_line (const char *line, int len);
extern bool sam_inspect_txt_header (BufferP txt_header);

// SEG Stuff
extern void sam_seg_initialize (VBlockP vb);
extern void sam_seg_finalize (VBlockP vb);
extern const char *sam_seg_txt_line (VBlockP vb_, const char *field_start_line, bool *has_special_eol);
extern uint32_t sam_seg_seq_len_from_cigar (const char *cigar, unsigned cigar_len);
extern uint32_t sam_seg_get_seq_len_by_MD_field (const char *md_str, unsigned md_str_len);

// PIZ Stuff
extern bool sam_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
extern void sam_piz_reconstruct_vb ();
extern void sam_piz_reconstruct_seq (VBlockP vb, ContextP ctx, const char *unused, unsigned unused2);

// VB stuff
extern void sam_vb_release_vb();
extern void sam_vb_destroy_vb();
extern unsigned sam_vb_size (void);
extern unsigned sam_vb_zip_dl_size (void);

#define SAM_SPECIAL { sam_piz_special_CIGAR, sam_piz_special_TLEN, sam_piz_special_BI, sam_piz_special_AS, sam_piz_special_MD }
SPECIAL (SAM, 0, CIGAR, sam_piz_special_CIGAR);
SPECIAL (SAM, 1, TLEN,  sam_piz_special_TLEN);
SPECIAL (SAM, 2, BI,    sam_piz_special_BI);
SPECIAL (SAM, 3, AS,    sam_piz_special_AS);
SPECIAL (SAM, 4, MD,    sam_piz_special_MD);
#define NUM_SAM_SPECIAL 5

// SAM field types 
#define sam_dict_id_is_qname_sf  dict_id_is_type_1
#define sam_dict_id_is_optnl_sf  dict_id_is_type_2

#define sam_dict_id_qname_sf     dict_id_type_1
#define sam_dict_id_optnl_sf     dict_id_type_2

#define SAM_DICT_ID_ALIASES \
    /*         alias                           maps to this ctx          */   \
    { DT_SAM,  &dict_id_fields[SAM_RNEXT],     &dict_id_fields[SAM_RNAME]  }, \
    { DT_SAM,  &dict_id_OPTION_MC,             &dict_id_fields[SAM_CIGAR]  }, \
    { DT_SAM,  &dict_id_OPTION_OC,             &dict_id_fields[SAM_CIGAR]  }, \
    { DT_SAM,  &dict_id_OPTION_E2,             &dict_id_fields[SAM_SEQ_BITMAP] }, \
    { DT_SAM,  &dict_id_OPTION_U2,             &dict_id_fields[SAM_QUAL]   },

#define SAM_LOCAL_GET_LINE_CALLBACKS  \
    { DT_SAM,   &dict_id_OPTION_BD,          sam_zip_get_start_len_line_i_bd    }, \
    { DT_SAM,   &dict_id_OPTION_BI,          sam_zip_get_start_len_line_i_bi    }, \
    { DT_SAM,   &dict_id_fields[SAM_QUAL],   sam_zip_get_start_len_line_i_qual  }, 

#endif
