// ------------------------------------------------------------------
//   kraken.h
//   Copyright (C) 2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef KRAKEN_INCLUDED
#define KRAKEN_INCLUDED

#include "genozip.h"

#define _KRAKEN_CU        DICT_ID_MAKEF_2 ("CU")
#define _KRAKEN_QNAME     DICT_ID_MAKEF_5 ("QNAME")
#define _KRAKEN_TAXID     DICT_ID_MAKEF_5 ("TAXID")
#define _KRAKEN_SEQLEN    DICT_ID_MAKEF_6 ("SEQLEN")
#define _KRAKEN_KMERS     DICT_ID_MAKEF_5 ("KMERS")
#define _KRAKEN_KMERTAX   DICT_ID_MAKEF_L ("KMERTAX")
#define _KRAKEN_KMERLEN   DICT_ID_MAKEF_L ("KMERLEN")
#define _KRAKEN_EOL       DICT_ID_MAKEF_3 ("EOL")
#define _KRAKEN_TOPLEVEL  DICT_ID_MAKEF_L (TOPLEVEL)
#define _KRAKEN_TOP2TAXID DICT_ID_MAKEF_L ("TOP2HASH")

typedef enum { KRAKEN_CU, KRAKEN_QNAME, KRAKEN_TAXID, KRAKEN_SEQLEN, KRAKEN_KMERS, KRAKEN_KMERTAX, KRAKEN_KMERLEN, 
               KRAKEN_EOL, KRAKEN_TOPLEVEL, KRAKEN_TOP2TAXID, NUM_KRAKEN_FIELDS } KrakenFields;

#define KRAKEN_MAPPING { V(KRAKEN_CU), V(KRAKEN_QNAME), V(KRAKEN_TAXID), V(KRAKEN_SEQLEN), V(KRAKEN_KMERS), V(KRAKEN_KMERTAX), V(KRAKEN_KMERLEN), V(KRAKEN_EOL), V(KRAKEN_TOPLEVEL), V(KRAKEN_TOP2TAXID) }

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

#endif
