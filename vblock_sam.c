// ------------------------------------------------------------------
//   vblock_sam.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "vblock.h"

void vb_sam_release_vb (VBlockSAM *vb)
{
    vb->last_pos = 0;
    memset (&vb->qname_mapper, 0, sizeof (vb->qname_mapper));
    
    // note: vb->data_line is not freed but rather used by subsequent vbs
    if (command == ZIP && vb->data_lines)
        memset (vb->data_lines, 0, sizeof(ZipDataLineSAM) * vb->num_lines_alloced);
}

void vb_sam_cleanup_memory (VBlockSAM *vb)
{

}
