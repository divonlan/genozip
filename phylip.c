// ------------------------------------------------------------------
//   phylip.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

// returns header length if header read is complete + sets lines.len, -1 not complete yet 
int32_t phy_is_header_done (bool is_eof)
{
    ARRAY (char, header, evb->txt_data);

    for (uint32_t i=0; i < evb->txt_data.len; i++)
        if (header[i] == '\n') 
            return i+1;

    // case: the entire file is just a header
    if (is_eof && *LASTENT (char, evb->txt_data) == '\n') 
        return evb->txt_data.len;

    return -1;
}

bool phy_header_inspect (BufferP txt_header)
{
    ARRAY (char, header, evb->txt_data);

    uint32_t num_seqs;
    int ret = sscanf (header, "%u %u", &num_seqs, &phy_seq_len);
    ASSINP (ret==2, "Error: invalid Phylip header line: \"%.*s\"", (int)txt_header->len, txt_header->data);

    return true;
}

//----------------------
// Compression functions
//----------------------

typedef struct {
    uint32_t line_start;  // start within vb->txt_data
} ZipDataLinePHY;
#define DATA_LINE(i) ENT (ZipDataLinePHY, vb->lines, i)

unsigned phy_vb_zip_dl_size (void) { return sizeof (ZipDataLinePHY); }

// callback functions for codec_* to get data of one line
void phy_zip_id (VBlock *vb, uint32_t vb_line_i, char **seq_data,  uint32_t *seq_len, uint32_t maximum_len) 
{
    ZipDataLinePHY *dl = DATA_LINE (vb_line_i);
    *seq_len  = PHY_ID_LEN;
    *seq_data = ENT (char, vb->txt_data, dl->line_start);
}

void phy_zip_seq (VBlock *vb, uint32_t vb_line_i, char **seq_data,  uint32_t *seq_len, uint32_t maximum_len) 
{
    ZipDataLinePHY *dl = DATA_LINE (vb_line_i);
    *seq_len  = phy_seq_len;
    *seq_data = ENT (char, vb->txt_data, dl->line_start) + PHY_ID_LEN;
}

//-----------------------
// Segmentation functions
//-----------------------

void phy_seg_initialize (VBlock *vb)
{
    vb->contexts[PHY_SEQ].ltype = LT_SEQUENCE;
    vb->contexts[PHY_ID].ltype  = LT_SEQUENCE;
    vb->contexts[PHY_TOPLEVEL].no_stons  = true; // keep in b250 so it can be eliminated as all_the_same
    vb->contexts[PHY_TOP2FASTA].no_stons = true;
}

void phy_seg_finalize (VBlockP vb)
{
    // top level snip
    SmallContainer top_level = { 
        .repeats     = vb->lines.len,
        .is_toplevel = true,
        .nitems_lo   = 3,
        .items       = { { .dict_id = (DictId)dict_id_fields[PHY_ID]  },
                         { .dict_id = (DictId)dict_id_fields[PHY_SEQ] },
                         { .dict_id = (DictId)dict_id_fields[PHY_EOL] } }
    };

    container_seg_by_ctx (vb, &vb->contexts[PHY_TOPLEVEL], (ContainerP)&top_level, 0, 0, 0);

    SmallContainer top_level_to_fasta = { 
        .repeats   = vb->lines.len,
        .is_toplevel = true,
        .nitems_lo = 2,
        .items     = { { .dict_id = (DictId)dict_id_fields[PHY_ID],  .seperator = "\n", .translator = PHYLIP2FASTA_ID },
                       { .dict_id = (DictId)dict_id_fields[PHY_SEQ], .seperator = "\n\n" } }
    };

    static const char fasta_prefix[] = { CON_PREFIX_SEP,        // has prefix 
                                         CON_PREFIX_SEP,        // end of (empty) container-wide prefix
                                         '>', CON_PREFIX_SEP }; // sequence ID prefix in fasta

    container_seg_by_ctx (vb, &vb->contexts[PHY_TOP2FASTA], (ContainerP)&top_level_to_fasta, fasta_prefix, sizeof (fasta_prefix), 0);
}

bool phy_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return true; // contexts are expected to have small dictionaries
}

const char *phy_seg_txt_line (VBlock *vb, const char *line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    Context *id_ctx  = &vb->contexts[PHY_ID];
    Context *seq_ctx = &vb->contexts[PHY_SEQ];

    ASSSEG0 (remaining_txt_len >= PHY_ID_LEN + phy_seq_len + 1 /* newline */, line, "Phylip data ends abruptly");

    *has_13 = (line[PHY_ID_LEN + phy_seq_len] == '\r');

    ASSSEG (line[PHY_ID_LEN + phy_seq_len + *has_13] == '\n', line, "Expecting the line to be exactly %u characters long", PHY_ID_LEN + phy_seq_len);

    DATA_LINE (vb->line_i)->line_start = line - vb->txt_data.data;

    // ID
    static char id_lookup[3] = { SNIP_LOOKUP, '1', '0' }; // lookup 10 characters from local
    seg_by_ctx (vb, id_lookup, sizeof (id_lookup), id_ctx, PHY_ID_LEN);
    id_ctx ->local.len += PHY_ID_LEN;

    // SEQ
    static char seq_lookup[10] = { SNIP_LOOKUP }; 
    unsigned snip_len = 1 + str_int (phy_seq_len, &seq_lookup[1]);
    seg_by_ctx (vb, seq_lookup, snip_len, seq_ctx, phy_seq_len);
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
    while (vb->txt_data.len && *LASTENT (char, vb->txt_data)==' ')
        vb->txt_data.len--;

    return 0;
}

