// ------------------------------------------------------------------
//   txtheader.c
//   Copyright (C) 2o21-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "buffer.h"
#include "strings.h"
#include "vblock.h"
#include "file.h"

static Buffer biopsy_vb_i = EMPTY_BUFFER;
static Buffer biopsy_data = EMPTY_BUFFER;
static char *biopsy_fn = NULL;

void biopsy_init (const char *optarg)
{
    str_split (optarg, strlen(optarg), 0, ',', item, false);

    ASSINP (n_items > 0, "Invalid biopsy argument: \"%s\", expecting a comma-seperated list of VB numbers, with 0 meaning the txt header", optarg);

    for (int i=0; i < n_items; i++) {
        str_split (items[i], item_lens[i], 2, '-', startend, false);
        if (n_startends == 2) { // a range, eg "20-25"
            uint32_t first_vb_i = atoi (startends[0]);
            uint32_t last_vb_i  = atoi (startends[1]);

            buf_alloc (evb, &biopsy_vb_i, last_vb_i-first_vb_i+1, 0, uint32_t, 0, "biopsy");

            for (uint32_t vb_i = first_vb_i; vb_i <= last_vb_i; vb_i++)
                NEXTENT (uint32_t, biopsy_vb_i) = vb_i;
        }
        else { // single vb eg "1"    
            buf_alloc (evb, &biopsy_vb_i, 1, n_items, uint32_t, 0, "biopsy");
            NEXTENT (uint32_t, biopsy_vb_i) = atoi (items[i]);
        }
    }

    biopsy_fn = malloc (strlen(optarg)+50);
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

    return; // we were not requested to take a biopsy from vb

start_biopsy:
    
    buf_add_more (evb, &biopsy_data, STRb(vb->txt_data), "biopsy_data");
    
    DO_ONCE sprintf (&biopsy_fn[strlen(biopsy_fn)], "%s", file_plain_ext_by_dt (txt_file->data_type));

    // case: biopsy is ready - dump it to a file and exit
    if (!biopsy_vb_i.len) {
        buf_dump_to_file (biopsy_fn, &biopsy_data, 1, false, false, true, true);
        exit_ok();
    }
}
