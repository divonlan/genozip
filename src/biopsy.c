// ------------------------------------------------------------------
//   txtheader.c
//   Copyright (C) 2o21-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include "genozip.h"
#include "buffer.h"
#include "strings.h"
#include "vblock.h"
#include "file.h"
#include "biopsy.h"

static Buffer biopsy_vb_i = { .name = "biopsy_vb_i" };
static StrTextSuperLong biopsy_fn;
static int biopsy_fn_len = 0;

void biopsy_init (rom optarg)
{
    str_split (optarg, strlen(optarg), 0, ',', item, false);

    ASSINP (n_items > 0, "Invalid biopsy argument: \"%s\", expecting a comma-seperated list of VB numbers, with 0 meaning the txt header", optarg);

    if (!item_lens[n_items-1]) n_items--; // remove final empty item in case of terminal ',' 

    for (int i=0; i < n_items; i++) {
        str_split (items[i], item_lens[i], 2, '-', startend, false);
        if (n_startends == 2) { // a range, eg "20-25"
            uint32_t first_vb_i = atoi (startends[0]);
            uint32_t last_vb_i  = atoi (startends[1]);

            for (VBIType vb_i = first_vb_i; vb_i <= last_vb_i; vb_i++)
                buf_add_int (evb, biopsy_vb_i, vb_i); 

            biopsy_vb_i.count += (last_vb_i - first_vb_i + 1); // counts VBs (not MAIN/PRIM/DEPN)        
        }

        else if (str_issame_(STRi(item,i), "MAIN", 4))
            buf_add_int (evb, biopsy_vb_i, (int32_t)-SAM_COMP_MAIN-1); // -1

        else if (str_issame_(STRi(item,i), "PRIM", 4))
            buf_add_int (evb, biopsy_vb_i, (int32_t)-SAM_COMP_PRIM-1); // -2

        else if (str_issame_(STRi(item,i), "DEPN", 4))
            buf_add_int (evb, biopsy_vb_i, (int32_t)-SAM_COMP_DEPN-1); // -3
             
        else { // single vb eg "1" or "0" (i.e. txt header only)
            buf_add_int (evb, biopsy_vb_i, ((int32_t)atoi (items[i]))); 
            biopsy_vb_i.count++;
        }
    }

    SNPRINTF (biopsy_fn, "%s.biopsy", optarg);

    flag.biopsy   = true;
    flag.seg_only = true; // no need to write the genozip file (+ we don't even seg, see zip_compress_one_vb)
}

static void biopsy_compress (void)
{
    iprint0 ("Biopsy: compressing biopsy file\n");   
    
    file_gzip (biopsy_fn.s);
    
    iprintf ("Biopsy: Done. Biopsy file is %s\n", biopsy_fn.s);    
}

static bool data_exhausted = false;
void biopsy_data_is_exhausted (void)
{
    data_exhausted = true;
}

void biopsy_take (VBlockP vb)
{
    if (!flag.biopsy || !Ltxt) return;

    ARRAY (int32_t, vb_i, biopsy_vb_i);

    for (int i=0; i < vb_i_len; i++)
        if (vb_i[i] == vb->vblock_i) {
            memmove (&vb_i[i], &vb_i[i+1], sizeof(uint32_t) * (vb_i_len-(i+1))); // remove from list
            biopsy_vb_i.count--;

            goto start_biopsy;
        }

        else if (-vb_i[i]-1 == vb->comp_i) 
            goto start_biopsy;


    // always output the txt header except if --no-header
    if (vb->vblock_i == 0 && !flag.no_header) 
        goto start_biopsy; 

    return; // we were not requested to take a biopsy from this vb

start_biopsy: {
    DO_ONCE {
        SNPRINTF (biopsy_fn, "%s", file_plain_ext_by_dt (txt_file->data_type));
        file_remove (biopsy_fn.s, true); // remove old file
    }

    if (biopsy_fn_len == sizeof (biopsy_fn.s)) goto too_long;
    
    // append to file
    FILE *fp = fopen (biopsy_fn.s, "ab");

    if (!fp && errno == ENAMETOOLONG) too_long: { 
        biopsy_fn_len = 0;
        SNPRINTF (biopsy_fn, "biopsy%s", file_plain_ext_by_dt (txt_file->data_type));
        file_remove (biopsy_fn.s, true); // remove old file
        fp = fopen (biopsy_fn.s, "ab");
    }

    ASSERT (fp, "failed to open biopsy file %s: %s", biopsy_fn.s, strerror (errno));

    ASSERT (fwrite (B1STtxt, 1, Ltxt, fp) == Ltxt, "failed to write %u bytes to biopsy file %s: %s", 
            Ltxt, biopsy_fn.s, strerror (errno));

    ASSERT (!fclose (fp), "failed to close biopsy file %s: %s", biopsy_fn.s, strerror (errno)); 
    
    if (vb->vblock_i == 0) iprint0 ("Biopsy: wrote txt_header\n");
    else                   iprintf ("Biopsy: wrote vblock=%s\n", VB_NAME);

    // case: biopsy is ready (works only with VB list, not components) - dump it to a file and exit
    if (!biopsy_vb_i.len) { 
        biopsy_compress();
        exit_ok;
    }
}
}

bool biopsy_is_done (void)
{
    bool ret = !biopsy_vb_i.len ||                       // biopsy of only VBs (no components) - and all done
               (!biopsy_vb_i.count && data_exhausted); // no VBs remaining, only components - and we're at EOF - done

    return ret;
}

// called if not already finalized after completing VB list (i.e. if not all VBs were encountered, or
// components were listed)
void biopsy_finalize (void)
{
    biopsy_compress();

    ASSINP0 (biopsy_is_done(), "FYI: One or more of biopsy VBs was not encountered");
}
