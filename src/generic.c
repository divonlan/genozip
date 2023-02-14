// ------------------------------------------------------------------
//   generic.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "vblock.h"
#include "buffer.h"
#include "seg.h"
#include "generic.h"
#include "dict_id.h"
#include "file.h"
#include "tar.h"

#define MAGIC_SIZE 32
static char magic[MAGIC_SIZE] = {}; // first 8 bytes of the generic file
static char ext[32] = {}; // nul-terminated txt filename extension

// all data is always consumed
int32_t generic_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    return 0;
}

// we use the header callback to try to detect if this generic file can be recognized as another
// data type by a signature. if not, it remains generic, and we set the header length to 0
int32_t generic_is_header_done (bool is_eof)
{
    ARRAY (char, header, evb->txt_data);
    DataType new_dt = DT_NONE;
    
    if (!header_len)
        return is_eof ? 0 : HEADER_NEED_MORE;
    
    SAFE_NUL(&header[header_len]);
    bool need_more = false;

    // search for a data type who's signature is in this header
    if (!flag.explicitly_generic)
        for (DataType dt=0; dt < NUM_DATATYPES; dt++)
            if (dt_props[dt].is_data_type && dt_props[dt].is_data_type (STRa(header), &need_more)) {
                new_dt = dt;
                break;
            }

    SAFE_RESTORE;

    if (new_dt != DT_NONE) {
        txt_file->data_type = z_file->data_type = new_dt;

        // recreate predefined contexts
        ctx_initialize_predefined_ctxs (z_file->contexts, new_dt, z_file->dict_id_to_did_i_map, &z_file->num_contexts);
        
        return HEADER_DATA_TYPE_CHANGED;
    }

    else if (need_more && !is_eof) // header too short to determine - we need more data
        return HEADER_NEED_MORE;

    else {
        // if we can't recognize a piped-in file, then we require the user to tell us what it is with --input
        ASSINP (!txt_file->redirected || flag.stdin_type, 
                "to pipe data in, please use --input (or -i) to specify its type, which can be one of the following:\n%s", file_compressible_extensions (true));
        
        if (!flag.stdin_type && !tar_is_tar()) 
            WARN_ONCE ("FYI: genozip doesn't recognize %s file's type, so it will be compressed as GENERIC. In the future, you may specify the type with \"--input <type>\". To suppress this warning, use \"--input generic\".", txt_name);

        return 0;
    }
}

void generic_seg_initialize (VBlockP vb)
{
    // capture the first MAGIC_SIZE bytes and the extension to be reported in stats
    if (vb->vblock_i == 1) {
        memset (magic, 0, MAGIC_SIZE);
        memcpy (magic, B1STtxt, MIN_(MAGIC_SIZE, vb->txt_data.len32));

        // copy return last component if it is not the whole filename, and short (we want the extension that indicates the file type, not the part of the filename that indicates the data)
        memset (ext, 0, sizeof (ext));
        rom last_dot = txt_file->name ? strrchr (txt_file->name, '.') : NULL;
        if (last_dot) { 

            // if extension is gz, bz2 or xz - add the prior filename component
            rom last_dot2 = NULL;
            if (!strcmp (last_dot+1, "gz") || !strcmp (last_dot+1, "bz2") || !strcmp (last_dot+1, "xz")) {
                SAFE_NUL(last_dot);
                last_dot2 = strrchr (txt_file->name, '.');
                SAFE_RESTORE;
            }
         
            strncpy (ext, (last_dot2 ? last_dot2 : last_dot) + 1, sizeof(ext)-1);
        }
    }
}

void generic_seg_finalize (VBlockP vb)
{
    Context *data_ctx = CTX(GNRIC_DATA);
    data_ctx->ltype = LT_UINT8;
    buf_move (vb, &data_ctx->local, vb, &vb->txt_data);
    data_ctx->txt_len += data_ctx->local.len;

    Context *toplevel_ctx = CTX(GNRIC_TOPLEVEL);
    toplevel_ctx->no_stons = true; // keep in b250 so it can be eliminated as all_the_same
    
    static const char snip[2] = { SNIP_SPECIAL, GNRIC_SPECIAL_TOPLEVEL };
    seg_by_ctx (VB, snip, 2, toplevel_ctx, 0); 
}

bool generic_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return true; // contexts are expected to have small dictionaries
}

SPECIAL_RECONSTRUCTOR (generic_piz_TOPLEVEL)
{
    buf_destroy (vb->txt_data);
    buf_move (vb, &vb->txt_data, vb, &CTX(GNRIC_DATA)->local);
    return NO_NEW_VALUE;
}

rom generic_get_magic (void)
{
    static char s[MAGIC_SIZE * 4 + 10];
    s[0] = '"';
    str_to_printable (magic, sizeof(magic), &s[1]);

    int len = strlen(s);
    s[len] = '"';
    s[len+1] = ' ';

    str_to_hex ((bytes)magic, sizeof(magic), &s[len+2], true);

    return s;
}

rom generic_get_ext (void)
{
    return ext;
}
 