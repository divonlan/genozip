// ------------------------------------------------------------------
//   txtheader.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"

//--------------
// General stuff
//--------------

void txtheader_finalize (void);

//----------
// ZIP stuff
//----------

extern bool txtheader_zip_read_and_compress (uint64_t *txt_header_size);
extern ConstBufferP txtheader_get_contigs (void);
extern const char *txtheader_get_contig_name (uint32_t index, uint32_t *snip_len);
extern void txtheader_add_contig (const char *chrom_name, unsigned chrom_name_len, PosType last_pos, void *liftover_);
extern void txtheader_verify_contig_init (void);
extern void txtheader_verify_contig (const char *chrom_name, unsigned chrom_name_len, PosType last_pos, void *liftover_);
extern uint32_t txtheader_get_bound_headers_len(void); // for stats
extern void txtheader_alloc_contigs (uint32_t more_contigs, uint32_t more_dict_len, bool liftover);
extern void txtheader_zip_prepopulate_contig_ctxs (void);

//----------
// PIZ stuff
//----------

extern Coords txtheader_piz_read_and_reconstruct (uint32_t component_i, Section txt_header_sl);
