// ------------------------------------------------------------------
//   seg.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

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

// returns length of next expected item, and 0 if unsuccessful
static unsigned gff3_seg_get_aofs_item_len (const char *str, unsigned len, bool is_last_item)
{
    unsigned i=0; for (; i < len; i++) {
        bool is_comma = str[i] == ',';
        bool is_space = str[i] == ' ';

        if ((is_last_item && is_comma) || (!is_last_item && is_space)) return i; // we reached the end of the item - return its length

        if ((is_last_item && is_space) || (!is_last_item && is_comma)) return 0; // invalid string - unexpected seperator
    }

    if (is_last_item) return i; // last item of last entry
    else return 0; // the string ended prematurely - this is not yet the last item
}

// a field that looks like: "non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
// we have an array (2 entires in this example) of items (4 in this examples) - the entries are separated by comma and the items by space
// observed in Ensembel generated GVF: Variant_effect, sift_prediction, polyphen_prediction, variant_peptide
// The last item is treated as an ENST_ID (format: ENST00000399012) while the other items are regular dictionaries
// the names of the dictionaries are the same as the ctx, with the 2nd character replaced by 1,2,3...
// the field itself will contain the number of entries
static void gff3_seg_array_of_struct (VBlock *vb, Context *subfield_ctx, 
                                      SmallContainer con, 
                                      const char *snip, unsigned snip_len)
{
    bool is_last_entry = false;
    uint32_t num_items = con_nitems (con);

    // get ctx's
    Context *ctxs[MAX_ENST_ITEMS] = {}; // an array of length num_items_in_struct (pointer to start of sub-array in vb->contexts)
    for (unsigned i=0; i < num_items; i++) {
        ctxs[i] = ctx_get_ctx (vb, con.items[i].dict_id); 
        ctxs[i]->st_did_i = subfield_ctx->did_i;
    }

    // set roll back point
    uint64_t saved_node_i_len[MAX_ENST_ITEMS], saved_local_len[MAX_ENST_ITEMS], saved_txt_len[MAX_ENST_ITEMS];
    for (unsigned item_i=0; item_i < num_items; item_i++) {
        saved_node_i_len[item_i] = ctxs[item_i]->b250.len;
        saved_local_len[item_i]  = ctxs[item_i]->local.len;
        saved_txt_len  [item_i]  = ctxs[item_i]->txt_len;
    }
    const char *saved_snip = snip;
    unsigned saved_snip_len = snip_len;

    con.repeats = 0;

    while (snip_len) {
        
        for (unsigned item_i=0; item_i < num_items; item_i++) {
            bool is_last_item = (item_i == num_items-1);
            unsigned item_len = gff3_seg_get_aofs_item_len (snip, snip_len, is_last_item);
            if (!item_len) goto badly_formatted;

            if (!is_last_item)
                seg_by_ctx (vb, snip, item_len, ctxs[item_i], item_len);
            else {
                is_last_entry = (snip_len - item_len == 0);
                seg_id_field (vb, ENSTid, snip, item_len, false);
            }
    
            snip     += item_len + 1 - is_last_entry; // 1 for either the , or the ' ' (except in the last item of the last entry)
            snip_len -= item_len + 1 - is_last_entry;
        }

        if (!is_last_entry && snip[-1]!=',') goto badly_formatted; // expecting a , after the end of all items in this entry
        
        con.repeats++;

        ASSSEG (con.repeats <= CONTAINER_MAX_REPEATS, snip, "exceeded maximum repeats allowed (%u) while parsing %s",
                CONTAINER_MAX_REPEATS, subfield_ctx->tag_name);
    }

    // finally, the Container snip itself
    container_seg (vb, subfield_ctx, (ContainerP)&con, NULL, 0, 
                           con.repeats * (num_items-1) /* space seperators */ + con.repeats-1 /* comma separators */);

    return;

badly_formatted:
    // roll back all the changed data
    for (unsigned item_i=0; item_i < num_items ; item_i++) {
        ctxs[item_i]->b250.len  = saved_node_i_len[item_i];
        ctxs[item_i]->local.len = saved_local_len[item_i];
        ctxs[item_i]->txt_len   = saved_txt_len[item_i];
    }

    // now save the entire snip in the dictionary
    seg_by_dict_id (vb, saved_snip, saved_snip_len, subfield_ctx->dict_id, saved_snip_len); 
}                           

bool gff3_seg_special_info_subfields (VBlockP vb, DictId dict_id, const char **this_value, unsigned *this_value_len)
{
    // ID - this is a sequential number (at least in GRCh37/38)
    if (dict_id.num == _ATTR_ID) {
        seg_pos_field ((VBlockP)vb, ATTR_ID, ATTR_ID, SPF_BAD_SNIPS_TOO, 0, *this_value, *this_value_len, 0, *this_value_len);
        return false; // do not add to dictionary/b250 - we already did it
    }

    // Dbxref (example: "dbSNP_151:rs1307114892") - we divide to the non-numeric part which we store
    // in a dictionary and the numeric part which store in a NUMERICAL_ID_DATA section
    if (dict_id.num == _ATTR_Dbxref) {
        seg_id_field (vb, ATTR_Dbxref, *this_value, *this_value_len, false); // discard the const as seg_id_field modifies
        return false; // do not add to dictionary/b250 - we already did it
    }

    // subfields that are arrays of structs, for example:
    // "non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
    if (dict_id.num == _ATTR_Variant_effect) {
        static const SmallContainer Variant_effect = {
            .nitems_lo   = 4, 
            .drop_final_repeat_sep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="V0arEff" }, .seperator = {' '} },
                             { .dict_id={.id="V1arEff" }, .seperator = {' '} },
                             { .dict_id={.id="V2arEff" }, .seperator = {' '} },
                             { .dict_id={.id="ENSTid"  },                    } }
        };
        gff3_seg_array_of_struct (vb, CTX(ATTR_Variant_effect), Variant_effect, *this_value, *this_value_len);
        return false;
    }

    if (dict_id.num == _ATTR_sift_prediction) {
        static const SmallContainer sift_prediction = {
            .nitems_lo   = 4, 
            .drop_final_repeat_sep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="S0iftPr" }, .seperator = {' '} },
                             { .dict_id={.id="S1iftPr" }, .seperator = {' '} },
                             { .dict_id={.id="S2iftPr" }, .seperator = {' '} },
                             { .dict_id={.id="ENSTid"  },                    } }
        };
        gff3_seg_array_of_struct (vb, CTX(ATTR_sift_prediction), sift_prediction, *this_value, *this_value_len);
        return false;
    }

    if (dict_id.num == _ATTR_polyphen_prediction) {
        static const SmallContainer polyphen_prediction = {
            .nitems_lo   = 4, 
            .drop_final_repeat_sep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="P0olyPhP" }, .seperator = {' '} },
                             { .dict_id={.id="P1olyPhP" }, .seperator = {' '} },
                             { .dict_id={.id="P2olyPhP" }, .seperator = {' '} },
                             { .dict_id={.id="ENSTid"   },                    } }
        };
        gff3_seg_array_of_struct (vb, CTX(ATTR_polyphen_prediction), polyphen_prediction, *this_value, *this_value_len);
        return false;
    }

    if (dict_id.num == _ATTR_variant_peptide) {
        static const SmallContainer variant_peptide = {
            .nitems_lo   = 3, 
            .drop_final_repeat_sep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="v0arPep" }, .seperator = {' '} }, // small v to differentiate from Variant_effect, so that dict_id to did_i mapper can map both
                             { .dict_id={.id="v1arPep" }, .seperator = {' '} },
                             { .dict_id={.id="ENSTid"  },                    } }
        };
        gff3_seg_array_of_struct (vb, CTX(ATTR_variant_peptide), variant_peptide, *this_value, *this_value_len);
        return false;
    }

    // we store these 3 in one dictionary, as they are correlated and will compress better together
    if (dict_id.num == _ATTR_Variant_seq   ||
        dict_id.num == _ATTR_Reference_seq ||
        dict_id.num == _ATTR_ancestral_allele) {

        // note: all three are stored together in _ATTR_Reference_seq as they are correlated
        seg_add_to_local_text (vb, CTX(ATTR_Reference_seq), *this_value, *this_value_len, *this_value_len);
        
        return false; // do not add to dictionary/b250 - we already did it
    }

    // Optimize Variant_freq
    unsigned optimized_snip_len = *this_value_len + 20;
    char optimized_snip[optimized_snip_len]; // used for 1. fields that are optimized 2. fields translated luft->primary

    if (flag.optimize_Vf && (dict_id.num == _ATTR_Variant_freq) &&
        optimize_float_2_sig_dig (*this_value, *this_value_len, 0, optimized_snip, &optimized_snip_len)) {
        
        vb->recon_size -= (int)(*this_value_len) - (int)optimized_snip_len;
        *this_value = optimized_snip;
        *this_value_len = optimized_snip_len;
        return true; // proceed with adding to dictionary/b250
    }

    return true; // all other cases -  procedue with adding to dictionary/b250
}

typedef struct { const char *start; 
                 unsigned len; 
                 DictId dict_id;   } AttrsItem;

static int sort_by_subfield_name (const void *a, const void *b)  
{ 
    AttrsItem *ina = (AttrsItem *)a;
    AttrsItem *inb = (AttrsItem *)b;
    
    return strncmp (ina->start, inb->start, MIN (ina->len, inb->len));
}

static void gff3_seg_attrs_field (VBlock *vb, const char *info_str, unsigned info_len)
{
    Container con = { .repeats             = 1, 
                      .drop_final_item_sep = true };

    const char *this_name = info_str, *this_value = NULL;
    int this_name_len = 0, this_value_len=0; // int and not unsigned as it can go negative

    AttrsItem info_items[MAX_FIELDS];

    // get name / value pairs - and insert values to the "name" dictionary
    bool reading_name = true;
    for (unsigned i=0; i < info_len + 1; i++) {
        char c = (i==info_len) ? ';' : info_str[i]; // add an artificial ; at the end of the INFO data

        if (reading_name) {

            if (c == '=' || c == ';') {  // end of valueful or valueless name

                bool valueful = (c == '=');

                ASSSEG0 (this_name_len > 0, info_str, "Error: ATTRS contains a = or ; without a preceding subfield name");

                if (this_name_len > 0) { 
                    ASSSEG ((this_name[0] >= 64 && this_name[0] <= 127) || this_name[0] == '.', info_str,
                            "ATTRS contains a name %.*s starting with an illegal character '%c' (ASCII %u)", 
                            this_name_len, this_name, this_name[0], this_name[0]);

                    DictId dict_id = valueful ? dict_id_make (this_name, this_name_len, DTYPE_1) : DICT_ID_NONE;

                    info_items[con_nitems(con)] = (AttrsItem) {
                        .start   = this_name,
                        .len     = this_name_len + valueful, // include the '=' if there is one 
                        .dict_id = dict_id
                    };

                    // create context with tag, if it doesn't already exist
                    if (valueful) ctx_get_ctx_tag (vb, dict_id, this_name, this_name_len);

                    this_value = &info_str[i+1]; 
                    this_value_len = -valueful; // if there is a '=' to be skipped, start from -1
                    reading_name = false; 
                }
            }
            else this_name_len++; // don't count the = or ; in the len
        }
        
        if (!reading_name) {

            if (c == ';') { // end of value
                // If its a valueful item, seg it (either special or regular)
                DictId dict_id = info_items[con_nitems(con)].dict_id;
                if (dict_id.num && gff3_seg_special_info_subfields (vb, dict_id, &this_value, (unsigned *)&this_value_len))
                    seg_by_dict_id (vb, this_value, this_value_len, dict_id, this_value_len);

                reading_name = true;  // end of value - move to the next item
                this_name = &info_str[i+1]; // move to next field in info string
                this_name_len = 0;
                con_inc_nitems (con);

                ASSSEG (con_nitems(con) <= MAX_FIELDS, info_str, 
                        "A line has too many subfields (tags) in ATTRS - the maximum supported is %u", MAX_FIELDS);
            }
            else this_value_len++;
        }
    }

    // if requested, we will re-sort the info fields in alphabetical order. This will result less words in the dictionary
    // thereby both improving compression and improving --regions speed. 
    if (flag.optimize_sort && con_nitems(con) > 1) 
        qsort (info_items, con_nitems(con), sizeof(AttrsItem), sort_by_subfield_name);

    char prefixes[CONTAINER_MAX_PREFIXES_LEN]; // these are the Container prefixes
    prefixes[0] = prefixes[1] = CON_PREFIX_SEP; // initial CON_PREFIX_SEP follow by seperator of empty Container-wide prefix
    unsigned prefixes_len = 2;

    // Populate the Container 
    uint32_t total_names_len=0;
    for (unsigned i=0; i < con_nitems(con); i++) {
        // Set the Container item and find (or create) a context for this name
        AttrsItem *ii = &info_items[i];
        con.items[i] = (ContainerItem){ .dict_id   = ii->dict_id,
                                        .seperator = { ';' } }; 
            
        // add to the prefixes
        ASSSEG (prefixes_len + ii->len + 1 <= CONTAINER_MAX_PREFIXES_LEN, info_str, 
                "ATTRS contains tag names that, combined (including the '='), exceed the maximum of %u characters", CONTAINER_MAX_PREFIXES_LEN);

        memcpy (&prefixes[prefixes_len], ii->start, ii->len);
        prefixes_len += ii->len;
        prefixes[prefixes_len++] = CON_PREFIX_SEP;

        total_names_len += ii->len + 1; // +1 for ; \t or \n separator
    }

    container_seg (vb, CTX(GFF3_ATTRS), &con, prefixes, prefixes_len, total_names_len /* names inc. = and separator */);
}

// GFF3 file format: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
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
    SEG_NEXT_ITEM (GFF3_PHASE);

    if (separator != '\n') {
        GET_LAST_ITEM (GFF3_ATTRS); 
        gff3_seg_attrs_field (vb, field_start, field_len); // GFF3 ATTRS has the same format and VCF INFO
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
