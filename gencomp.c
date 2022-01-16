// ------------------------------------------------------------------
//   gencomp.c - "generated component"
//   Copyright (C) 2021-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <errno.h>
#include "genozip.h"
#include "file.h"
#include "gencomp.h"

// ZIP: called when inspecting the txtheader to add header data, and after each VB to add rejected line
void gencomp_append_file (VBlockP vb, 
                          GenCompNum file_gc, // 1 or 2
                          const char *name)
{
    ASSERTNOTNULL (z_file);
    int file_i = file_gc - 1;

    // create gencomp file (DVCF: rejects SAM: Primary/Dependent) if not already open
    if (!z_file->gencomp_file[file_i]) {
        z_file->gencomp_file_name[file_i] = MALLOC (strlen (z_file->name) + 20);
        sprintf (z_file->gencomp_file_name[file_i], "%s.%s_ONLY%s", z_file->name, name, file_plain_ext_by_dt (z_file->data_type));
        
        z_file->gencomp_file[file_i] = fopen (z_file->gencomp_file_name[file_i], "wb+");
        ASSERT (z_file->gencomp_file[file_i], "fopen() failed to open %s: %s", z_file->gencomp_file_name[file_i], strerror (errno));
    }

    ASSERT0 (z_file->gencomp_file[file_i], "liftover rejects file is not open");

    fwrite (vb->gencomp[file_i].data, 1, vb->gencomp[file_i].len, z_file->gencomp_file[file_i]);
    z_file->gencomp_disk_size[file_i] += vb->gencomp[file_i].len;

    buf_free (&vb->gencomp[file_i]);

    flag.bind = BIND_GENCOMP; // we need to bind the generated components 
}
