// ------------------------------------------------------------------
//   me23.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt


#pragma once

#include "genozip.h"

#pragma GENDICT_PREFIX ME23
#pragma GENDICT ME23_CHROM=DTYPE_FIELD=CHROM
#pragma GENDICT ME23_POS=DTYPE_FIELD=POS
#pragma GENDICT ME23_ID=DTYPE_FIELD=ID
#pragma GENDICT ME23_GENOTYPE=DTYPE_FIELD=GENOTYPE
#pragma GENDICT ME23_EOL=DTYPE_FIELD=EOL
#pragma GENDICT ME23_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT ME23_TOP2VCF=DTYPE_FIELD=TOP2VCF
#pragma GENDICT ME23_DEBUG_LINES=DTYPE_FIELD=DBGLINES      // used by --debug-lines

extern rom me23_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void me23_seg_initialize (VBlockP vb);
extern void me23_seg_finalize (VBlockP vb);
extern bool me23_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern bool me23_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);

// translator numbers must start from 1 - 0 is reserved for "none"
TRANSLATOR (ME23, VCF, 1, GENOTYPE, sam_piz_m232vcf_GENOTYPE)   // reconstruct VCF GENOTYPE field as VCF - REF,ALT,QUAL,FILTER,INFO,FORMAT,Sample
#define NUM_ME23_TRANS 2 // including "none"
#define ME23_TRANSLATORS { NULL /* none */, sam_piz_m232vcf_GENOTYPE }

TXTHEADER_TRANSLATOR (txtheader_me232vcf);
