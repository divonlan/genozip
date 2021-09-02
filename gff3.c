// ------------------------------------------------------------------
//   gff3.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

// GFF3 specification and examples: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
// Another one: http://gmod.org/wiki/GFF3 and https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "random_access.h"
#include "file.h"
#include "strings.h"
#include "optimize.h"
#include "piz.h"
#include "dict_id.h"
#include "codec.h"
#include "vcf.h"
#include "dict_id_gen.h"

#define MAX_ENST_ITEMS 10 // maximum number of items in an enst structure. this can be changed without impacting backward compatability.

// called from seg_all_data_lines
void gff3_seg_initialize (VBlock *vb)
{
    CTX(GFF3_SEQID)->flags.store   = STORE_INDEX; // since v12
    CTX(GFF3_START)->flags.store   = STORE_INT;   // since v12
    CTX(GFF3_COMMENT)->flags.store = STORE_INDEX; // COMMENT introduced in 12.0.12
    CTX(GFF3_COMMENT)->no_stons    = true; // required by STORE_INDEX (otherwise singletons don't get their index stored)
    CTX(GFF3_SEQID)->no_stons      = true; // needs b250 node_index for random access
    CTX(GFF3_ATTRS)->no_stons      = true;
    CTX(GFF3_TOPLEVEL)->no_stons   = true; // keep in b250 so it can be eliminated as all_the_same

    CTX(ATTR_Target_ID)->st_did_i = CTX(ATTR_Target_POS)->st_did_i = CTX(ATTR_Target_STRAND)->st_did_i = ATTR_Target;

}

void gff3_seg_finalize (VBlockP vb)
{
    // top level snip
    SmallContainer top_level = { 
        .repeats   = vb->lines.len,
        .is_toplevel  = true,
        .filter_items = true,
        .nitems_lo = 11,
        .items     = { { .dict_id = { _GFF3_COMMENT },                  },
                       { .dict_id = { _GFF3_SEQID },  .seperator = "\t" },
                       { .dict_id = { _GFF3_SOURCE }, .seperator = "\t" },
                       { .dict_id = { _GFF3_TYPE },   .seperator = "\t" },
                       { .dict_id = { _GFF3_START },  .seperator = "\t" },
                       { .dict_id = { _GFF3_END },    .seperator = "\t" },
                       { .dict_id = { _GFF3_SCORE },  .seperator = "\t" },
                       { .dict_id = { _GFF3_STRAND }, .seperator = "\t" },
                       { .dict_id = { _GFF3_PHASE },  .seperator = "\t" },
                       { .dict_id = { _GFF3_ATTRS },                    },
                       { .dict_id = { _GFF3_EOL },                      } }
    };

    container_seg (vb, CTX(GFF3_TOPLEVEL), (Container *)&top_level, 0, 0, 0);
}

bool gff3_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == _GFF3_TOPLEVEL  ||
           dict_id.num == _GFF3_SEQID     ||
           dict_id.num == _GFF3_SOURCE    ||
           dict_id.num == _GFF3_TYPE      ||
           dict_id.num == _GFF3_END       ||
           dict_id.num == _GFF3_SCORE     ||
           dict_id.num == _GFF3_STRAND    ||
           dict_id.num == _GFF3_PHASE     ||
           dict_id.num == _GFF3_ATTRS     ||
           dict_id.num == _GFF3_EOL;
}

// returns trus if successful
static bool gff3_seg_target (VBlockP vb, const char *value, unsigned value_len)
{
    str_split (value, value_len, 4, ' ', item, false);
    if (n_items != 3 && n_items != 4) return false; // not standard Target format

    SmallContainer con = {
        .repeats   = 1,
        .nitems_lo = n_items,
        .drop_final_item_sep = true,
        .items     = { { .dict_id = { _ATTR_Target_ID     }, .seperator = " " },
                       { .dict_id = { _ATTR_Target_POS    }, .seperator = " " }, // START
                       { .dict_id = { _ATTR_Target_POS    }, .seperator = " " }, // END
                       { .dict_id = { _ATTR_Target_STRAND }, .seperator = " " } }
    };

    seg_by_did_i (vb, items[0], item_lens[0], ATTR_Target_ID, item_lens[0]);
    seg_pos_field (vb, ATTR_Target_POS, ATTR_Target_POS, 0, 0, items[1], item_lens[1], 0, item_lens[1]);
    seg_pos_field (vb, ATTR_Target_POS, ATTR_Target_POS, 0, 0, items[2], item_lens[2], 0, item_lens[2]);

    if (n_items == 4)
        seg_by_did_i (vb, items[3], item_lens[3], ATTR_Target_STRAND, item_lens[3]);

    // TO DO: use seg_duplicate_last if n_items hasn't changed
    container_seg (vb, CTX (ATTR_Target), (ContainerP)&con, NULL, 0, n_items-1); // account for space separators

    return true;
}

static inline DictId gff3_seg_attr_subfield (VBlockP vb, const char *tag_name, unsigned tag_name_len, const char *value, unsigned value_len)
{
    DictId dict_id = dict_id_make (tag_name, tag_name_len, DTYPE_GFF3_ATTR);

    ContextP ctx = ctx_get_ctx_tag (vb, dict_id, tag_name, tag_name_len);
    
    switch (dict_id.num) {

    // ID - sometimes this is a sequential number (GRCh37/38)
    // sometimes it is something like this: c5312581-5d6e-4234-89d7-4974581f2993
    case _ATTR_ID: 
        if (str_is_int (value, value_len))
            seg_pos_field (vb, ATTR_ID, ATTR_ID, SPF_BAD_SNIPS_TOO, 0, value, value_len, 0, value_len);
        else
            seg_by_did_i (vb, value, value_len, ATTR_ID, value_len);
        break;

    // Dbxref (example: "dbSNP_151:rs1307114892") - we divide to the non-numeric part which we store
    // in a dictionary and the numeric part which store in a local
    case _ATTR_Dbxref:
        seg_id_field (vb, ATTR_Dbxref, value, value_len, false); // discard the const as seg_id_field modifies
        break;

    case _ATTR_Target:
        if (!gff3_seg_target (vb, value, value_len))
            goto plain_seg;
        break;

    case _ATTR_Name:
        seg_id_field (vb, ATTR_Name, value, value_len, false);
        break;

    // example: Parent=mRNA00001,mRNA00002,mRNA00003
    case _ATTR_Parent:
        seg_array (vb, CTX(ATTR_Parent), ATTR_Parent, value, value_len, ',', 0, false, false, false);
        break;

    //case _ATTR_Gap: // I tried: 1. array (no improvement) ; 2. string of op-codes in b250 + integers in local (negligible improvement)

    // subfields that are arrays of structs, for example:
    // "non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
    case _ATTR_Variant_effect: {
        static const SmallContainer Variant_effect = {
            .nitems_lo   = 4, 
            .drop_final_repeat_sep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="V0arEff" }, .seperator = {' '} },
                             { .dict_id={.id="V1arEff" }, .seperator = {' '} },
                             { .dict_id={.id="V2arEff" }, .seperator = {' '} },
                             { .dict_id={.num=_ENSTid  },                    } }
        };
        seg_array_of_struct (vb, CTX(ATTR_Variant_effect), Variant_effect, value, value_len, true);
        break;
    }

    case _ATTR_sift_prediction: {
        static const SmallContainer sift_prediction = {
            .nitems_lo   = 4, 
            .drop_final_repeat_sep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="S0iftPr" }, .seperator = {' '} },
                             { .dict_id={.id="S1iftPr" }, .seperator = {' '} },
                             { .dict_id={.id="S2iftPr" }, .seperator = {' '} },
                             { .dict_id={.num=_ENSTid  },                    } }
        };
        seg_array_of_struct (vb, CTX(ATTR_sift_prediction), sift_prediction, value, value_len, true);
        break;
    }

    case _ATTR_polyphen_prediction: {
        static const SmallContainer polyphen_prediction = {
            .nitems_lo   = 4, 
            .drop_final_repeat_sep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="P0olyPhP" }, .seperator = {' '} },
                             { .dict_id={.id="P1olyPhP" }, .seperator = {' '} },
                             { .dict_id={.id="P2olyPhP" }, .seperator = {' '} },
                             { .dict_id={.num=_ENSTid   },                    } }
        };
        seg_array_of_struct (vb, CTX(ATTR_polyphen_prediction), polyphen_prediction, value, value_len, true);
        break;
    }

    case _ATTR_variant_peptide: {
        static const SmallContainer variant_peptide = {
            .nitems_lo   = 3, 
            .drop_final_repeat_sep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="v0arPep"  }, .seperator = {' '} }, // small v to differentiate from Variant_effect, so that dict_id to did_i mapper can map both
                             { .dict_id={.id="v1arPep"  }, .seperator = {' '} },
                             { .dict_id={.num=_ENSTid   },                    } }
        };
        seg_array_of_struct (vb, CTX(ATTR_variant_peptide), variant_peptide, value, value_len, true);
        break;
    }

    // we store these 3 in one dictionary, as they are correlated and will compress better together
    case _ATTR_Variant_seq:
    case _ATTR_Reference_seq:
    case _ATTR_ancestral_allele: 
        // note: all three are stored together in _ATTR_Reference_seq as they are correlated
        seg_add_to_local_text (vb, CTX(ATTR_Reference_seq), value, value_len, value_len); 
        break;

    // Optimize Variant_freq
    case _ATTR_Variant_freq:
        if (flag.optimize_Vf) {
            unsigned optimized_snip_len = value_len + 20;
            char optimized_snip[optimized_snip_len]; // used for 1. fields that are optimized 2. fields translated luft->primary

            if (optimize_float_2_sig_dig (value, value_len, 0, optimized_snip, &optimized_snip_len)) {        
                vb->recon_size -= (int)(value_len) - (int)optimized_snip_len;
                value = optimized_snip;
                value_len = optimized_snip_len;
            }            
        }
        goto plain_seg; // proceed with adding to dictionary/b250

    default:
    plain_seg:
        seg_by_ctx (vb, value, value_len, ctx, value_len);
    }

    return dict_id;
}

static void gff3_seg_attrs_field (VBlock *vb, const char *field, unsigned field_len)
{
    // case: "." field
    if (field_len == 1 && *field == '.') {
        seg_by_did_i (vb, ".", 1, GFF3_ATTRS, 2);
        return;
    }

    Container con = { .repeats             = 1, 
                      .drop_final_item_sep = true };
    char prefixes[field_len + 2];
    prefixes[0] = prefixes[1] = CON_PREFIX_SEP;
    unsigned prefixes_len = 2;

    str_split (field, field_len, MAX_DICTS, ';', attr, false);
    ASSSEG (n_attrs, field, "Invalid attributes field: %.*s", STRf(field));
    
    con_set_nitems (con, n_attrs);

    for (unsigned i=0; i < n_attrs; i++) {

        str_split (attrs[i], attr_lens[i], 2, '=', tag_val, true);
        ASSSEG (n_tag_vals == 2 && tag_val_lens[0] > 0 && tag_val_lens[1] > 0, attrs[i], "Invalid attribute: %.*s", STRfi (attr, i));

        memcpy (&prefixes[prefixes_len], tag_vals[0], tag_val_lens[0] + 1);
        prefixes_len += tag_val_lens[0] + 2; // including the '=' and CON_PREFIX_SEP
        prefixes[prefixes_len-1] = CON_PREFIX_SEP;

        DictId dict_id = gff3_seg_attr_subfield (vb, STRi(tag_val, 0), STRi (tag_val, 1));
        con.items[i] = (ContainerItem){ .dict_id = dict_id, .seperator = { ';' } }; 
    }

    container_seg (vb, CTX(GFF3_ATTRS), &con, prefixes, prefixes_len, (prefixes_len-2) /* names inc. = and (; or \n) seperator */);
}

const char *gff3_seg_txt_line (VBlock *vb, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;
    int32_t len = (int32_t)remaining_txt_len;

    // handle a comment (#), directive (##) or empty line 
    if (next_field[0] == '#' || next_field[0] == '\n' || next_field[0] == '\r') {
        SEG_NEXT_ITEM_NL (GFF3_COMMENT);

        // TO DO: support embedded FASTA, bug 392
        ASSINP (!str_issame_ (GFF3_COMMENT_str, GFF3_COMMENT_len, "##FASTA", 7), "%s contains a ##FASTA directive. To compress it use --input=generic", txt_name);
        
        goto eol; // if we have a comment, then during piz, the other fields will be filtered out
    }
    else 
        seg_by_did_i (vb, NULL, 0, GFF3_COMMENT, 0); // missing comment field

    GET_NEXT_ITEM (GFF3_SEQID);
    seg_chrom_field (vb, field_start, field_len);

    SEG_NEXT_ITEM (GFF3_SOURCE);
    SEG_NEXT_ITEM (GFF3_TYPE);

    GET_NEXT_ITEM (GFF3_START);
    seg_pos_field (vb, GFF3_START, GFF3_START, 0, 0, field_start, field_len, 0, field_len+1);
    random_access_update_pos (vb, DC_PRIMARY, GFF3_START);

    GET_NEXT_ITEM (GFF3_END);
    seg_pos_field (vb, GFF3_END, GFF3_START, 0, 0, field_start, field_len, 0, field_len+1);

    SEG_NEXT_ITEM (GFF3_SCORE);
    SEG_NEXT_ITEM (GFF3_STRAND);
    SEG_MAYBE_LAST_ITEM (GFF3_PHASE);

    if (separator != '\n') { // "attrbutes" is an optional field per http://gmod.org/wiki/GFF3
        GET_LAST_ITEM (GFF3_ATTRS); 
        gff3_seg_attrs_field (vb, field_start, field_len);
    }
    else
        seg_by_did_i (vb, NULL, 0, GFF3_ATTRS, 0); // NULL=MISSING so previous \t is removed

eol:             
    SEG_EOL (GFF3_EOL, false);

    return next_field;
}

//----------
// PIZ stuff
//----------

// filter is called before reconstruction of a repeat or an item, and returns true if item should be reconstructed. if not reconstructed, contexts are not consumed.
CONTAINER_FILTER_FUNC (gff3_piz_filter)
{
    if (item < 1 || item == 10) return true; // we filter all items except COMMENT and EOL

    return CTX(GFF3_COMMENT)->last_value.i == WORD_INDEX_MISSING; // reconstruct non-comment contexts only if this isn't a comment line
}
