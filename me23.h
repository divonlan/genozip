// ------------------------------------------------------------------
//   me23.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt


#ifndef ME23_INCLUDED
#define ME23_INCLUDED

#include "genozip.h"

#define _ME23_CHROM     DICT_ID_MAKEF_5 ("CHROM")
#define _ME23_POS       DICT_ID_MAKEF_3 ("POS")
#define _ME23_ID        DICT_ID_MAKEF_2 ("ID")
#define _ME23_GENOTYPE  DICT_ID_MAKEF_L ("GENOTYPE")
#define _ME23_EOL       DICT_ID_MAKEF_3 ("EOL")
#define _ME23_TOPLEVEL  DICT_ID_MAKEF_L (TOPLEVEL)
#define _ME23_TOP2VCF   DICT_ID_MAKEF_7 ("TOP2VCF")

typedef enum { ME23_CHROM, ME23_POS, ME23_ID, ME23_GENOTYPE, ME23_EOL, ME23_TOPLEVEL, ME23_TOP2VCF, NUM_ME23_FIELDS } Me23Fields;  

#define ME23_MAPPING { V(ME23_CHROM), V(ME23_POS), V(ME23_ID), V(ME23_GENOTYPE), V(ME23_EOL), V(ME23_TOPLEVEL), V(ME23_TOP2VCF) }

extern const char *me23_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void me23_seg_initialize (VBlockP vb);
extern void me23_seg_finalize (VBlockP vb);
extern bool me23_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern bool me23_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);

#define ME23_DICT_ID_ALIASES

// translator numbers must start from 1 - 0 is reserved for "none"
TRANSLATOR (ME23, VCF, 1, GENOTYPE, sam_piz_m232vcf_GENOTYPE)   // reconstruct VCF GENOTYPE field as VCF - REF,ALT,QUAL,FILTER,INFO,FORMAT,Sample
#define NUM_ME23_TRANS 2 // including "none"
#define ME23_TRANSLATORS { NULL /* none */, sam_piz_m232vcf_GENOTYPE }

TXTHEADER_TRANSLATOR (txtheader_me232vcf);

#endif
