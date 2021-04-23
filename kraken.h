// ------------------------------------------------------------------
//   kraken.h
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef KRAKEN_INCLUDED
#define KRAKEN_INCLUDED

#include "genozip.h"

typedef int32_t TaxonomyId;
#define TAXID_NONE ((TaxonomyId)-1)

// seg of a kraken file
extern void kraken_zip_initialize (void);
extern void kraken_seg_initialize (VBlockP vb);
extern void kraken_seg_finalize (VBlockP vb);
extern const char *kraken_seg_txt_line (VBlockP vb, const char *line, uint32_t remaining_txt_len, bool *has_13);
extern bool kraken_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern void kraken_zip_after_compute (VBlockP vb);

// piz of a kraken file
extern bool kraken_piz_initialize (void);
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
