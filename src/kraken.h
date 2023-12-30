// ------------------------------------------------------------------
//   kraken.h
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

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
#pragma GENDICT KRAKEN_Q8NAME=DTYPE_1=Q8NAME 
#pragma GENDICT KRAKEN_Q9NAME=DTYPE_1=Q9NAME 
#pragma GENDICT KRAKEN_QANAME=DTYPE_1=QANAME 
#pragma GENDICT KRAKEN_QBNAME=DTYPE_1=QBNAME 
#pragma GENDICT KRAKEN_QCNAME=DTYPE_1=QCNAME 
#pragma GENDICT KRAKEN_QDNAME=DTYPE_1=QDNAME 
#pragma GENDICT KRAKEN_QENAME=DTYPE_1=QENAME // if adding more Q*NAMEs - add to fastq.h and sam.h and update MAX_QNAME_ITEMS
#pragma GENDICT KRAKEN_QmNAME=DTYPE_1=QmNAME // QmNAME reserved for mate number (always the last dict_id in the container)

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
#pragma GENDICT KRAKEN_DEBUG_LINES=DTYPE_FIELD=DBGLINES      // used by --debug-lines

#define TAXID_NONE         ((TaxonomyId)-1)
#define TAXID_UNCLASSIFIED ((TaxonomyId)-0)

// VB management
extern unsigned kraken_vb_size (DataType dt);

// seg of a kraken file
extern void kraken_zip_initialize (void);
extern void kraken_seg_initialize (VBlockP vb);
extern void kraken_seg_finalize (VBlockP vb);
extern rom kraken_seg_txt_line (VBlockP vb, rom line, uint32_t remaining_txt_len, bool *has_13);
extern bool kraken_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern void kraken_zip_after_compute (VBlockP vb);

// piz of a kraken file
extern bool kraken_is_translation (VBlockP vb);
extern bool kraken_piz_initialize (CompIType comp_i);
extern IS_SKIP (kraken_piz_is_skip_section);
extern CONTAINER_CALLBACK (kraken_piz_container_cb);
extern void kraken_piz_process_recon (VBlockP vb);

// using the kraken data in genocat --kraken
extern void kraken_set_taxid (rom optarg);
extern void kraken_load (void);
extern bool kraken_is_included_loaded (VBlockP vb, STRp(qname));
extern bool kraken_is_included_stored (VBlockP vb, Did did_i_taxid, bool already_reconstructed);
extern void kraken_destroy (void);

// using the kraken data in genozip --kraken
extern unsigned kraken_seg_taxid_do (VBlockP vb, Did did_i_taxid, STRp(qname), char *snip, bool fail_if_missing);
extern unsigned kraken_seg_taxid (VBlockP vb, Did did_i_taxid, STRp(qname), bool fail_if_missing);

// misc
extern void kraken_set_show_kraken (rom optarg);
extern char *kraken_filename; // global
#define kraken_is_loaded ((bool)(kraken_filename != NULL))

#define KRAKEN_SPECIAL {  }

#define NUM_KRAKEN_SPECIAL 0
