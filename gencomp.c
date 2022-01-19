// ------------------------------------------------------------------
//   gencomp.c - "generated component"
//   Copyright (C) 2021-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <errno.h>
#include "genozip.h"
#include "file.h"
#include "gencomp.h"

// create gencomp file (DVCF: rejects SAM: Primary/Dependent) if not already open
void gencomp_initialize_file (GenCompNum file_gc, const char *name)
{
    int file_i = file_gc - 1;

    if (z_file->gencomp_file[file_i]) return; // already initialized

    z_file->gencomp_file_name[file_i] = MALLOC (strlen (z_file->name) + 20);
    sprintf (z_file->gencomp_file_name[file_i], "%s.%s_ONLY%s", z_file->name, name, file_plain_ext_by_dt (z_file->data_type));
    
    z_file->gencomp_file[file_i] = fopen (z_file->gencomp_file_name[file_i], "wb+");
    ASSERT (z_file->gencomp_file[file_i], "fopen() failed to open %s: %s", z_file->gencomp_file_name[file_i], strerror (errno));

    // SAM/BAM: we add a SAM header. One benefit is it forces a 3-component file in case one the gencomp components are missing. The
    // header will marked as not "in_plan" in writer_init_comp_info, and only reconstructed in case of --component 2 or --component 3
    static const char *sam_gencomp_header[2] = {
        "@CO	PRIMARY - All alignments that have a supplementary or secondary alignment\n",
        "@CO	DEPENDENT - All supplementary and secondary alignments\n"
    };
    uint32_t l_text = strlen(sam_gencomp_header[file_i]);
    if (Z_DT(DT_SAM)) {
        fwrite (sam_gencomp_header[file_i], 1, l_text, z_file->gencomp_file[file_i]);
        z_file->gencomp_disk_size[file_i] += l_text;
    }
    // write a valid BAM header
    else if (Z_DT (DT_BAM)) {
        fwrite ("BAM\1", 1, 4, z_file->gencomp_file[file_i]);
        uint32_t l_text_lten = LTEN32 (l_text);
        fwrite (&l_text_lten, 1, 4, z_file->gencomp_file[file_i]);
        fwrite (sam_gencomp_header[file_i], 1, l_text, z_file->gencomp_file[file_i]);
        fwrite ("\0\0\0\0", 1, 4, z_file->gencomp_file[file_i]);

        z_file->gencomp_disk_size[file_i] += l_text + 12;
    }

    z_file->z_flags.has_gencomp = true;
}

// ZIP: called when inspecting the txtheader to add header data, and after each VB to add rejected line
void gencomp_append_file (VBlockP vb, 
                          GenCompNum file_gc) // 1 or 2
{
    ASSERTNOTNULL (z_file);
    
    int file_i = file_gc - 1;

    ASSERT (z_file->gencomp_file[file_i], "gencomp[%u] file is not open", file_i);

    if (vb->gencomp[file_i].len) {
        fwrite (vb->gencomp[file_i].data, 1, vb->gencomp[file_i].len, z_file->gencomp_file[file_i]);
        z_file->gencomp_disk_size[file_i] += vb->gencomp[file_i].len;
    }

    buf_free (&vb->gencomp[file_i]);

    flag.bind = BIND_GENCOMP; // we need to bind the generated components 
}
