// ------------------------------------------------------------------
//   phylip.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "file.h"
#include "random_access.h"
#include "endianness.h"
#include "strings.h"
#include "piz.h"
#include "dict_id.h"
#include "stats.h"
#include "reference.h"
#include "codec.h"
#include "version.h"
#include "phylip.h"

#define PHY_ID_LEN 10 

static uint32_t phy_seq_len = 0; // read from txt header

//-----------------
// Header functions
//-----------------

// ZIP: called from txtfile_read_header
// returns header length if header read is complete + sets lines.len, -1 not complete yet 
int32_t phy_is_header_done (bool is_eof)
{
    ARRAY (char, header, evb->txt_data);

    for (uint32_t i=0; i < evb->txt_data.len32; i++)
        if (header[i] == '\n') 
            return i+1;

    // case: the entire file is just a header
    if (is_eof && *BLSTc (evb->txt_data) == '\n') 
        return evb->txt_data.len;

    return HEADER_NEED_MORE;
}

bool phy_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags)
{
    ARRAY (char, header, *txt_header);

    uint32_t num_seqs;
    int ret = sscanf (header, "%u %u", &num_seqs, &phy_seq_len);
    ASSINP (ret==2, "Error: invalid PHYLIP header line: \"%.*s\"", (int)txt_header->len, txt_header->data);

    ASSERTW0 (num_seqs > 0, "FYI: unusual PHYLIP header: number_of_sequences==0");
    ASSERTW0 (phy_seq_len > 0, "FYI: unusual PHYLIP header: sequence_length==0");

    return true;
}

//----------------------
// Compression functions
//----------------------

typedef struct {
    uint32_t line_start;  // start within vb->txt_data
} ZipDataLinePHY;
#define DATA_LINE(i) B(ZipDataLinePHY, vb->lines, i)

unsigned phy_vb_zip_dl_size (void) { return sizeof (ZipDataLinePHY); }

// callback functions for codec_* to get data of one line
COMPRESSOR_CALLBACK (phy_zip_id)
{
    ZipDataLinePHY *dl = DATA_LINE (vb_line_i);
    *line_data_len = PHY_ID_LEN;
    *line_data     = Bc (vb->txt_data, dl->line_start);
    if (is_rev) *is_rev = 0;
}

COMPRESSOR_CALLBACK (phy_zip_seq)
{
    ZipDataLinePHY *dl = DATA_LINE (vb_line_i);
    *line_data_len = phy_seq_len;
    *line_data     = Bc (vb->txt_data, dl->line_start) + PHY_ID_LEN;
    if (is_rev) *is_rev = 0;
}

//-----------------------
// Segmentation functions
//-----------------------

void phy_seg_initialize (VBlockP vb)
{
    CTX(PHY_SEQ)->ltype = LT_SEQUENCE;
    CTX(PHY_ID)->ltype  = LT_SEQUENCE;
    CTX(PHY_TOPLEVEL)->no_stons  = true; // keep in b250 so it can be eliminated as all_the_same
    CTX(PHY_TOP2FASTA)->no_stons = true;
}

void phy_seg_finalize (VBlockP vb)
{
    // top level snip
    SmallContainer top_level = { 
        .repeats     = vb->lines.len,
        .is_toplevel = true,
        .nitems_lo   = 3,
        .items       = { { .dict_id = { _PHY_ID }  },
                         { .dict_id = { _PHY_SEQ } },
                         { .dict_id = { _PHY_EOL } } }
    };

    container_seg (vb, CTX(PHY_TOPLEVEL), (ContainerP)&top_level, 0, 0, 0);

    SmallContainer top_level_to_fasta = { 
        .repeats   = vb->lines.len,
        .is_toplevel = true,
        .nitems_lo = 2,
        .items     = { { .dict_id = { _PHY_ID },  .separator = "\n", .translator = PHYLIP2FASTA_ID },
                       { .dict_id = { _PHY_SEQ }, .separator = "\n\n" } }
    };

    static const char fasta_prefix[] = { CON_PX_SEP,        // has prefix 
                                         CON_PX_SEP,        // end of (empty) container-wide prefix
                                         '>', CON_PX_SEP }; // sequence ID prefix in fasta

    container_seg (vb, CTX(PHY_TOP2FASTA), (ContainerP)&top_level_to_fasta, fasta_prefix, sizeof (fasta_prefix), 0);
}

bool phy_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return true; // contexts are expected to have small dictionaries
}

rom phy_seg_txt_line (VBlockP vb, rom line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    Context *id_ctx  = CTX(PHY_ID);
    Context *seq_ctx = CTX(PHY_SEQ);

    ASSSEG0 (remaining_txt_len >= PHY_ID_LEN + phy_seq_len + 1 /* newline */, line, "Phylip data ends abruptly");

    *has_13 = (line[PHY_ID_LEN + phy_seq_len] == '\r');

    ASSSEG (line[PHY_ID_LEN + phy_seq_len + *has_13] == '\n', line, "Expecting the line to be exactly %u characters long", PHY_ID_LEN + phy_seq_len);

    DATA_LINE (vb->line_i)->line_start = line - vb->txt_data.data;

    // ID
    seg_by_ctx (VB, (char[]){ SNIP_LOOKUP, '1', '0' }, 3, id_ctx, PHY_ID_LEN); // lookup 10 characters from local
    id_ctx ->local.len += PHY_ID_LEN;

    // SEQ
    seg_lookup_with_length (vb, seq_ctx, phy_seq_len, phy_seq_len);
    seq_ctx->local.len += phy_seq_len;

    // EOL
    SEG_EOL (PHY_EOL, true);
    
    return line + PHY_ID_LEN + phy_seq_len + *has_13 + 1;
}

//---------------------------------------
// PHYLIP -> Multifasta translation stuff
//---------------------------------------

// remove Phylip header line
TXTHEADER_TRANSLATOR (txtheader_phy2fa)
{
    txtheader_buf->len = 0;
}

// Translating PHYLIP->FASTA: remove redundant terminating spaces from ID (FASTA's DESC)
TRANSLATOR_FUNC (phy_piz_phy2fasta_ID)
{
    while (vb->txt_data.len && *BLSTtxt==' ')
        vb->txt_data.len--;

    return 0;
}

