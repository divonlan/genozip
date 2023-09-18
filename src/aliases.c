// ------------------------------------------------------------------
//   aliases.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "aliases.h"
#include "file.h"
#include "zfile.h"
#include "buffer.h"
#include "qname.h"

static Buffer aliases = {};

static void show_aliases (void)
{
    if (!aliases.len)
        iprint0 ("No aliases in this file\n");

    else {
        static rom names[] = ALIAS_TYPE_NAMES;

        iprintf ("Contents of SEC_DICT_ID_ALIASES section (num_aliases=%u):\n", aliases.len32);
        
        for_buf (DictIdAlias, alias, aliases)
            iprintf ("type=%-4s\talias=%s/%-8s\tdst=%s/%-8s\n", names[alias->alias_type],
                    dtype_name_z (alias->alias), dis_dict_id (alias->alias).s, 
                    dtype_name_z (alias->dst),   dis_dict_id (alias->dst).s);
    }

    if (is_genocat) exit_ok;
}

// return predefined aliases of z_file->data_type
static void aliases_zip_get_predefined (bool remove_if_no_data) 
{
    static struct { DataType dt; 
                    AliasType alias_type;
                    uint64_t dict_id_alias; 
                    uint64_t dict_id_dst; } aliases_def[] = DICT_ID_ALIASES;

    buf_alloc (evb, &aliases, 0, ARRAY_LEN(aliases_def), DictIdAlias, 0, "aliases");

    for (int i=0; i < ARRAY_LEN(aliases_def); i++)
        if ((aliases_def[i].dt == z_file->data_type && 
             (!remove_if_no_data || ctx_get_existing_zctx (aliases_def[i].dict_id_dst)->z_data_exists))) // to reduce size, keep only aliases used in this z_file: it would be better to test existance of the alias rather than dst, but that's not easy to do. testing dst will be a little harmlessly wasteful bc in some cases we will write an alias section when unneeded

            BNXT (DictIdAlias, aliases) = (DictIdAlias){
                .alias_type = aliases_def[i].alias_type,
                .alias      = (DictId)aliases_def[i].dict_id_alias,
                .dst        = (DictId)aliases_def[i].dict_id_dst
            };
}

// ZIP main thread: write to global section
void aliases_compress (void)
{
    ASSERTNOTINUSE (aliases);
    aliases_zip_get_predefined (true);

    // add qname alias if there is one
    buf_alloc (evb, &aliases, NUM_QTYPES, 0, DictIdAlias, 0, NULL);
    for (QType q=QNAME1; q < NUM_QTYPES; q++) 
        if (qname_get_alias(q).alias_type)
            BNXT (DictIdAlias, aliases) = qname_get_alias(q);

    if (flag.show_aliases) show_aliases();

    if (aliases.len) {
        aliases.len *= sizeof (DictIdAlias);
        zfile_compress_section_data (evb, SEC_DICT_ID_ALIASES, &aliases);
    }

    buf_destroy (aliases);
}

static void aliases_upgrade_to_v15 (void)
{
    typedef struct { DictId alias, dst; } DictIdAliasV14; // up to v14

    aliases.len /= sizeof (DictIdAliasV14);
    buf_alloc (evb, &aliases, aliases.len, 0, char, 0, NULL); // add one byte per alias

    ARRAY (DictIdAliasV14, v14, aliases);
    ARRAY (DictIdAlias,    v15, aliases);

    for (int32_t i=v14_len-1; i >= 0; i--) {
        // note: explictly create a struct, then assign it. If done in a single operation, the compiler optimizer
        // messes it up, as v14 and v15 overlap
        DictIdAlias new = (DictIdAlias){ .alias_type = ALIAS_CTX,
                                         .alias      = v14[i].alias,
                                         .dst        = v14[i].dst    };
        v15[i] = new;
    }
}

// PIZ main thread: read all dict_id aliases, if there are any. Caller should free buffer.
void aliases_piz_read_aliases_section (void) 
{ 
    Section sec = sections_last_sec (SEC_DICT_ID_ALIASES, SOFT_FAIL);
    if (!sec) return; // no aliases section

    zfile_get_global_section (SectionHeader, sec, &aliases, "aliases");

    if (!VER(15)) 
        aliases_upgrade_to_v15();
    else
        aliases.len /= sizeof (DictIdAlias);

    for_buf (DictIdAlias, alias, aliases)
        ASSERT0 (alias->dst.id[0] && alias->alias.id[0], "corrupted aliases buffer"); // id[0] can never be 0

    if (flag.show_aliases) show_aliases();
}

// gets buffer of aliases - caller needs to destroy it. NULL if there are no aliases.
BufferP aliases_get (void) 
{ 
    ASSERTNOTINUSE (aliases);

    if (IS_ZIP) aliases_zip_get_predefined (false);
    else        aliases_piz_read_aliases_section();

    return &aliases;
}

