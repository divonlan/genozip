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
extern bool phy_header_inspect (BufferP txt_header);
extern void phy_seg_initialize (VBlockP vb);
extern void phy_seg_finalize (VBlockP vb);
extern bool phy_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern const char *phy_seg_txt_line (VBlockP vb, const char *line, uint32_t remaining_txt_len, bool *has_13);

// PIZ side
TXTHEADER_TRANSLATOR (txtheader_phy2fa);

#endif
