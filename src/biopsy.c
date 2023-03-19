// ------------------------------------------------------------------
//   txtheader.c
//   Copyright (C) 2o21-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "buffer.h"
#include "strings.h"
#include "vblock.h"
#include "file.h"

static Buffer biopsy_vb_i = { .name = "biopsy_vb_i" };
static Buffer biopsy_data = { .name = "biopsy_data", .can_be_big = true };
static char *biopsy_fn = NULL;

void biopsy_init (rom optarg)
{
    str_split (optarg, strlen(optarg), 0, ',', item, false);

    ASSINP (n_items > 0, "Invalid biopsy argument: \"%s\", expecting a comma-seperated list of VB numbers, with 0 meaning the txt header", optarg);

    for (int i=0; i < n_items; i++) {
        str_split (items[i], item_lens[i], 2, '-', startend, false);
        if (n_startends == 2) { // a range, eg "20-25"
            uint32_t first_vb_i = atoi (startends[0]);
            uint32_t last_vb_i  = atoi (startends[1]);

            for (VBIType vb_i = first_vb_i; vb_i <= last_vb_i; vb_i++)
                buf_add_int (evb, biopsy_vb_i, vb_i); 
        }
        else // single vb eg "1" or "0" (i.e. txt header only)
            buf_add_int (evb, biopsy_vb_i, ((uint32_t)atoi (items[i]))); 
    }

    biopsy_fn = MALLOC (strlen(optarg)+50);
    sprintf (biopsy_fn, "%s.biopsy", optarg);

    flag.biopsy   = true;
    flag.seg_only = true; // no need to write the genozip file (+ we don't even seg, see zip_compress_one_vb)
}

void biopsy_take (VBlockP vb)
{
    if (!biopsy_vb_i.len) return;

    ARRAY (uint32_t, vb_i, biopsy_vb_i);

    for (int i=0; i < vb_i_len; i++)
        if (vb_i[i] == vb->vblock_i) {
            memmove (&vb_i[i], &vb_i[i+1], sizeof(uint32_t) * (vb_i_len-(i+1))); // remove from list
            biopsy_vb_i.len--;

            goto start_biopsy;
        }

    if (vb->vblock_i == 0) goto start_biopsy; // always output the txt header

    return; // we were not requested to take a biopsy from this vb

start_biopsy:
    
    buf_add_buf (evb, &biopsy_data, &vb->txt_data, char, "biopsy_data");
    
    DO_ONCE sprintf (&biopsy_fn[strlen(biopsy_fn)], "%s", file_plain_ext_by_dt (txt_file->data_type));

    // case: biopsy is ready - dump it to a file and exit
    if (!biopsy_vb_i.len) {
        buf_dump_to_file (biopsy_fn, &biopsy_data, 1, false, false, true, true);
        exit_ok();
    }
}
