// ------------------------------------------------------------------
//   phylip.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef PHY_INCLUDED
#define PHY_INCLUDED

#include "genozip.h"

#define _PHY_ID         DICT_ID_MAKEF_2 ("ID")  
#define _PHY_SEQ        DICT_ID_MAKEF_3 ("SEQ") 
#define _PHY_EOL        DICT_ID_MAKEF_3 ("EOL")
#define _PHY_TOPLEVEL   DICT_ID_MAKEF_L (TOPLEVEL)
#define _PHY_TOP2FASTA  DICT_ID_MAKEF_6 ("TOP2FA")

typedef enum { PHY_ID, PHY_SEQ, PHY_EOL, PHY_TOPLEVEL, PHY_TOP2FASTA, NUM_PHY_FIELDS } PhyFields;
#define PHY_MAPPING { V(PHY_ID), V(PHY_SEQ), V(PHY_EOL), V(PHY_TOPLEVEL), V(PHY_TOP2FASTA) }

// ZIP side
COMPRESSOR_CALLBACK(phy_zip_id)
COMPRESSOR_CALLBACK(phy_zip_seq)
#define PHY_LOCAL_GET_LINE_CALLBACKS  \
    { DT_PHYLIP, _PHY_ID,  phy_zip_id  }, \
    { DT_PHYLIP, _PHY_SEQ, phy_zip_seq }, 

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
