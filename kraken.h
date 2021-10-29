// ------------------------------------------------------------------
//   kraken.h
//   Copyright (C) 2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

#pragma GENDICT_PREFIX KRAKEN

#pragma GENDICT KRAKEN_CU=DTYPE_FIELD=CU
#pragma GENDICT KRAKEN_QNAME=DTYPE_FIELD=QNAME
#pragma GENDICT KRAKEN_Q0NAME=DTYPE_1=Q0NAME // MAX_QNAME_ITEMS - fixed qname items must have a did_i directly after container's (MUST be the same dict_id as in sam.h)
#pragma GENDICT KRAKEN_Q1NAME=DTYPE_1=Q1NAME 
#pragma GENDICT KRAKEN_Q2NAME=DTYPE_1=Q2NAME
#pragma GENDICT KRAKEN_Q3NAME=DTYPE_1=Q3NAME
#pragma GENDICT KRAKEN_Q4NAME=DTYPE_1=Q4NAME
#pragma GENDICT KRAKEN_Q5NAME=DTYPE_1=Q5NAME
#pragma GENDICT KRAKEN_Q6NAME=DTYPE_1=Q6NAME 
#pragma GENDICT KRAKEN_Q7NAME=DTYPE_1=Q7NAME 
#pragma GENDICT KRAKEN_QmatNAME=DTYPE_1=QmatNAME // QmatNAME reserved for mate number (always the last dict_id in the container)

#pragma GENDICT KRAKEN_TAXID=DTYPE_FIELD=TAXID

#pragma GENDICT KRAKEN_SEQLEN=DTYPE_FIELD=SEQLEN
#pragma GENDICT KRAKEN_SEQLEN_1=DTYPE_2=S1EQLEN
#pragma GENDICT KRAKEN_SEQLEN_2=DTYPE_2=S2EQLEN

#pragma GENDICT KRAKEN_KMERS=DTYPE_FIELD=KMERS
#pragma GENDICT KRAKEN_KMERTAX=DTYPE_FIELD=KMERTAX
#pragma GENDICT KRAKEN_KMERLEN=DTYPE_FIELD=KMERLEN
#pragma GENDICT KRAKEN_EOL=DTYPE_FIELD=EOL
#pragma GENDICT KRAKEN_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT KRAKEN_TOP2TAXID=DTYPE_FIELD=TOP2HASH

typedef int32_t TaxonomyId;
#define TAXID_NONE         ((TaxonomyId)-1)
#define TAXID_UNCLASSIFIED ((TaxonomyId)-0)

// seg of a kraken file
extern void kraken_zip_initialize (void);
extern void kraken_seg_initialize (VBlockP vb);
extern void kraken_seg_finalize (VBlockP vb);
extern const char *kraken_seg_txt_line (VBlockP vb, const char *line, uint32_t remaining_txt_len, bool *has_13);
extern bool kraken_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern void kraken_zip_after_compute (VBlockP vb);

// piz of a kraken file
extern bool kraken_is_translation (VBlockP vb);
extern bool kraken_piz_initialize (void);
extern bool kraken_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
extern CONTAINER_CALLBACK (kraken_piz_container_cb);
extern void kraken_piz_handover_data (VBlockP vb);

// using the kraken data in genocat --kraken
extern void kraken_set_taxid (const char *optarg);
extern void kraken_load (void);
extern bool kraken_is_included_loaded (VBlockP vb, const char *qname, unsigned qname_len);
extern bool kraken_is_included_stored (VBlockP vb, DidIType did_i_taxid, bool already_reconstructed);
extern void kraken_destroy (void);

// using the kraken data in genozip --kraken
extern unsigned kraken_seg_taxid_do (VBlockP vb, DidIType did_i_taxid, const char *qname, unsigned qname_len, char *snip, bool fail_if_missing);
extern unsigned kraken_seg_taxid (VBlockP vb, DidIType did_i_taxid, const char *qname, unsigned qname_len, bool fail_if_missing);

// misc
extern void kraken_set_show_kraken (const char *optarg);
extern char *kraken_filename; // global
#define kraken_is_loaded ((bool)kraken_filename)
