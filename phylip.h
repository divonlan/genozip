// ------------------------------------------------------------------
//   phylip.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PHY_INCLUDED
#define PHY_INCLUDED

#include "genozip.h"

// ZIP side
COMPRESSOR_CALLBACK(phy_zip_id)
COMPRESSOR_CALLBACK(phy_zip_seq)
#define PHY_LOCAL_GET_LINE_CALLBACKS  \
    { DT_PHYLIP, &dict_id_fields[PHY_ID],  phy_zip_id  }, \
    { DT_PHYLIP, &dict_id_fields[PHY_SEQ], phy_zip_seq }, 

extern unsigned phy_vb_zip_dl_size (void);
extern int32_t phy_is_header_done (bool is_eof);
extern bool phy_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);
extern void phy_seg_initialize (VBlockP vb);
extern void phy_seg_finalize (VBlockP vb);
extern bool phy_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern const char *phy_seg_txt_line (VBlockP vb, const char *line, uint32_t remaining_txt_len, bool *has_13);

//---------------------------------------
// PHYLIP -> Multifasta translation stuff
//---------------------------------------

// Important: Numbers (and order) of translators cannot be changed, as they are part of the file format
// (they are included in the TOPLEVEL container)
// translator numbers must start from 1 - 0 is reserved for "none"
TRANSLATOR (PHYLIP, FASTA, 1, ID, phy_piz_phy2fasta_ID)   // remove final spaces

#define NUM_PHYLIP_TRANS 2 // including "none"
#define PHYLIP_TRANSLATORS { NULL /* none */, phy_piz_phy2fasta_ID }

TXTHEADER_TRANSLATOR (txtheader_phy2fa);

#endif
