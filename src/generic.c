// ------------------------------------------------------------------
//   generic.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "seg.h"
#include "file.h"
#include "tar.h"
#include "threads.h"
#include "tip.h"

#define MAGIC_SIZE 32
static char magic[MAGIC_SIZE+1] = {}; // first bytes of the generic file
static char ext[32] = {}; // nul-terminated txt filename extension

// all data is always consumed
int32_t generic_unconsumed (VBlockP vb, uint32_t first_i)
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

    // if new data type requires an external compressor, ask user to re-start (unelegant solution to a rare scenario)
    ASSINP0 (new_dt != DT_CRAM, "This is CRAM file, but Genozip got confused because the file's name doesn't end with .cram. Solution: re-run and add the command line option \"--input cram\"");
    ASSINP0 (new_dt != DT_BCF,  "This is BCF file, but Genozip got confused because the file's name doesn't end with .bcf. Solution: re-run and add the command line option \"--input bcf\"");
    
    if (new_dt != DT_NONE) {
        txt_file->data_type = z_file->data_type = new_dt;

        z_file->z_flags.txt_is_bin = DTPT (is_binary); // update

        // recreate predefined contexts
        ctx_initialize_predefined_ctxs (z_file->contexts, new_dt, z_file->d2d_map, &z_file->num_contexts);

        // note on FASTQ: effective_codec is currently GZ or BGZF. We will not try to re-discover in segconf as that is too complicated.        
        return HEADER_DATA_TYPE_CHANGED;
    }

    else if (need_more && !is_eof) // note: need_more=true, if no data type was identified, and at least one data type requested more data
        return HEADER_NEED_MORE;   

    if (flag.explicitly_generic)
        {} // no message

    else if (tar_zip_is_tar() || !txt_file->redirected) 
        WARN_ONCE ("FYI: genozip doesn't recognize %s file's type, so it will be compressed as GENERIC. In the future, you may specify the type with \"--input <type>\". To suppress this warning, use \"--input generic\".", txt_name);

    else 
        ABORTINP ("to pipe data in, please use --input (or -i) to specify its type, which can be one of the following:\n%s", file_compressible_extensions (true).s);

    return 0;
}

void generic_seg_initialize (VBlockP vb)
{
    // capture the first MAGIC_SIZE bytes and the extension to be reported in stats
    if (vb->vblock_i == 1) {
        memset (magic, 0, MAGIC_SIZE+1);
        memcpy (magic, B1STtxt, MIN_(MAGIC_SIZE, Ltxt));

        // copy return last component if it is not the whole filename, and short (we want the extension that indicates the file type, not the part of the filename that indicates the data)
        memset (ext, 0, sizeof (ext));
        rom last_dot = txt_file->name ? strrchr (txt_file->name, '.') : NULL;
        if (last_dot) { 

            // if extension is gz, bz2 or xz - add the prior filename component
            rom last_dot2 = NULL;
            if (!strcmp (last_dot+1, "gz") || !strcmp (last_dot+1, "bz2") || !strcmp (last_dot+1, "xz") || !strcmp (last_dot+1, "zip")) {
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
    ContextP data_ctx = CTX(GNRIC_DATA);
    data_ctx->ltype = LT_UINT8;
    buf_move (vb, data_ctx->local, CTX_TAG_LOCAL, vb->txt_data);
    data_ctx->txt_len += data_ctx->local.len;

    ContextP toplevel_ctx = CTX(GNRIC_TOPLEVEL);
    toplevel_ctx->no_stons = true; // keep in b250 so it can be eliminated as all_the_same
    
    static const char snip[2] = { SNIP_SPECIAL, GNRIC_SPECIAL_TOPLEVEL };
    seg_by_ctx (VB, snip, 2, toplevel_ctx, 0); 
}

rom generic_seg_txt_line (VBlockP vb, rom next_line, uint32_t remaining_txt_len, bool *has_13)
{
    return next_line + remaining_txt_len; // consumed entire VB
}

rom generic_assseg_line (VBlockP vb)
{
    return "GENERIC cannot display offending line"; 
}

bool generic_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return true; // contexts are expected to have small dictionaries
}

SPECIAL_RECONSTRUCTOR (generic_piz_TOPLEVEL)
{
    buf_destroy (vb->txt_data);
    buf_move (vb, vb->txt_data, "txt_data", CTX(GNRIC_DATA)->local);
    return NO_NEW_VALUE;
}

StrTextLong generic_get_magic (void)
{
    StrTextLong s = {};
    s.s[0] = '"';
    int len = 1 + str_to_printable (magic, strlen(magic), &s.s[1], sizeof(s.s) - MAGIC_SIZE*3 - 10);
    s.s[len++] = '"';
    s.s[len++] = ' ';

    str_to_hex_((bytes)magic, strlen(magic), &s.s[len], true);

    return s;
}

rom generic_get_ext (void)
{
    return ext;
}
 
// to be called from segconf of other data types, if it is discovered that the file is not actually of that data type
rom fallback_to_generic (VBlockP vb) // vb is optional
{
    ASSERT0 (threads_am_i_main_thread(), "fallback_to_generic can only be called by the main thread");

    SAVE_FLAG(quiet);

    if (!flag.explicit_quiet) flag.quiet = false; // cancel segconf's quietness
    WARN ("%s is not a valid %s file, compressing as GENERIC", txt_name, dt_name(txt_file->data_type));
    
    z_file->data_type = txt_file->data_type = DT_GNRIC;

    // recreate predefined contexts
    ctx_initialize_predefined_ctxs (z_file->contexts, DT_GNRIC, z_file->d2d_map, &z_file->num_contexts);

    RESTORE_FLAG(quiet);
    return vb ? BAFTtxt : NULL;
}
