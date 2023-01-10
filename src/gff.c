// ------------------------------------------------------------------
//   gff.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// GFF3 specification and examples: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
//                                  http://gmod.org/wiki/GFF3 
//                                  https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt
// GVF: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gvf.md
// GTF / GFF2: https://www.ensembl.org/info/website/upload/gff.html
//             http://gmod.org/wiki/GFF2

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "random_access.h"
#include "file.h"
#include "strings.h"
#include "optimize.h"
#include "piz.h"
#include "gff.h"
#include "dict_id.h"
#include "codec.h"
#include "vcf.h"
#include "dict_id_gen.h"
#include "chrom.h"
#include "stats.h"

#define MAX_ENST_ITEMS 10 // maximum number of items in an enst structure. this can be changed without impacting backward compatability.

// IMPORTANT: if changing fields in VBlockGFF, also update gff_vb_release_vb 
typedef struct VBlockGFF {    
    VBLOCK_COMMON_FIELDS
    STR(prev_type);
} VBlockGFF;

typedef VBlockGFF *VBlockGFFP;
#define VB_GFF ((VBlockGFFP)vb)

// cleanup vb (except common) and get it ready for another usage (without freeing memory held in the Buffers)
void gff_vb_release_vb (VBlockGFFP vb) 
{
    vb->prev_type     = 0;
    vb->prev_type_len = 0;
}

unsigned gff_vb_size (DataType dt) { return sizeof (VBlockGFF); }

sSTRl(copy_gene_name_snip,32);
sSTRl(transcript_name_container_snip,100);
void gff_zip_initialize (void)
{
    // transcript_name staff
    SmallContainer con = { 
        .repeats      = 1,
        .nitems_lo    = 2,
        .items        = { { .dict_id = { _ATTR_transcript_name_gene }, .separator = {'-'} },  // 0
                          { .dict_id = { _ATTR_transcript_name_num  }                     } } // 1
    };
    container_prepare_snip ((ContainerP)&con, 0, 0, qSTRa (transcript_name_container_snip));
    seg_prepare_snip_other (SNIP_COPY, _ATTR_gene_name, false, 0, copy_gene_name_snip);
}

// detect if a generic file is actually a GFF3/GVF/GTF (but cannot detect non-GTF GFF2)
bool is_gff (STRp(header), bool *need_more)
{
    SAFE_NULT (header);

    bool is_gff = (header_len >= 13 && !memcmp (header, "##gff-version", 13)) || // must appear for GFF3/GTF
                   strstr (header, "#!genome-build")                          || // normally appears in GTF / GVF
                   strstr (header, "#!genome-version");                          // normally appears in GTF / GVF

    SAFE_RESTORE;
    return is_gff;
}

bool gff_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags)
{
    SAFE_NUL (BAFTc (*txt_header)); // for atoi

    // case: gff version announced explicitly. otherwise we will deduce it during segconf
    if (txt_header->len >= 15 && !memcmp (B1STc(*txt_header), "##gff-version", 13))
        segconf.gff_version = atoi (Bc(*txt_header, 14));

    else if (txt_header->len >= 14 && !memcmp (B1STc(*txt_header), "#gtf-version", 12))
        segconf.gff_version = 2;

    SAFE_RESTORE;

    return true;
}

// called from seg_all_data_lines
void gff_seg_initialize (VBlockP vb)
{
    ctx_set_store (VB, STORE_INT, GFF_START/*since v12*/, GFF_END/*14.0.18*/, ATTR_Target_POS, ATTR_ID, ATTR_exon_number, DID_EOL);

    ctx_set_store (VB, STORE_INDEX, GFF_SEQID/*v12*/, GFF_COMMENT/*12.0.12*/, DID_EOL);

    ctx_set_no_stons (vb, GFF_COMMENT/*required by STORE_INDEX (otherwise singletons don't get their index stored)*/, 
                      GFF_SEQID/*needs b250 node_index for random access*/, 
                      GFF_ATTRS, ATTR_Variant_seq, ATTR_Reference_seq, ATTR_ancestral_allele,
                      ATTR_Target_POS, ATTR_ID, GFF_START, GFF_END, DID_EOL); // as requied by seg_pos_field

    ctx_set_ltype (vb, LT_DYN_INT, ATTR_transcript_name_num, DID_EOL);

    ctx_consolidate_stats (vb, ATTR_transcript_name, ATTR_transcript_name_gene, ATTR_transcript_name_num, DID_EOL);
    ctx_consolidate_stats (vb, ATTR_Target, ATTR_Target_ID, ATTR_Target_POS, ATTR_Target_STRAND, DID_EOL);

    seg_id_field_init (CTX(ATTR_Dbxref));
    seg_id_field_init (CTX(ATTR_db_xref));
    seg_id_field_init (CTX(ATTR_Name));
    seg_id_field_init (CTX(ATTR_exon_id));
    seg_id_field_init (CTX(ATTR_gene_id));
    seg_id_field_init (CTX(ATTR_protein_id));
    seg_id_field_init (CTX(ATTR_ccds_id));
    seg_id_field_init (CTX(ATTR_transcript_id));
}

void gff_seg_finalize (VBlockP vb)
{
    CTX(ENSTid)->st_did_i = DID_NONE; // cancel consolidatation as it goes into multiple attributes

    // top level snip
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .filter_items = true,
        .nitems_lo    = 11,
        .items        = { { .dict_id = { _GFF_COMMENT },                   },
                          { .dict_id = { _GFF_SEQID   }, .separator = "\t" },
                          { .dict_id = { _GFF_SOURCE  }, .separator = "\t" },
                          { .dict_id = { _GFF_TYPE    }, .separator = "\t" },
                          { .dict_id = { _GFF_START   }, .separator = "\t" },
                          { .dict_id = { _GFF_END     }, .separator = "\t" },
                          { .dict_id = { _GFF_SCORE   }, .separator = "\t" },
                          { .dict_id = { _GFF_STRAND  }, .separator = "\t" },
                          { .dict_id = { _GFF_PHASE   }, .separator = "\t" },
                          { .dict_id = { _GFF_ATTRS   },                   },
                          { .dict_id = { _GFF_EOL     },                   } }
    };

    container_seg (vb, CTX(GFF_TOPLEVEL), (Container *)&top_level, 0, 0, 0);
}

bool gff_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == _GFF_TOPLEVEL  ||
           dict_id.num == _GFF_SEQID     ||
           dict_id.num == _GFF_SOURCE    ||
           dict_id.num == _GFF_TYPE      ||
           dict_id.num == _GFF_END       ||
           dict_id.num == _GFF_SCORE     ||
           dict_id.num == _GFF_STRAND    ||
           dict_id.num == _GFF_PHASE     ||
           dict_id.num == _GFF_ATTRS     ||
           dict_id.num == _GFF_EOL;
}

// returns trus if successful
static bool gff_seg_target (VBlockP vb, STRp(value))
{
    str_split (value, value_len, 4, ' ', item, false);
    if (n_items != 3 && n_items != 4) return false; // not standard Target format

    SmallContainer con = {
        .repeats   = 1,
        .nitems_lo = n_items,
        .drop_final_item_sep = true,
        .items     = { { .dict_id = { _ATTR_Target_ID     }, .separator = " " },
                       { .dict_id = { _ATTR_Target_POS    }, .separator = " " }, // START
                       { .dict_id = { _ATTR_Target_POS    }, .separator = " " }, // END
                       { .dict_id = { _ATTR_Target_STRAND }, .separator = " " } }
    };

    seg_by_did (VB, items[0], item_lens[0], ATTR_Target_ID, item_lens[0]);
    seg_pos_field (vb, ATTR_Target_POS, ATTR_Target_POS, 0, 0, items[1], item_lens[1], 0, item_lens[1]);
    seg_pos_field (vb, ATTR_Target_POS, ATTR_Target_POS, 0, 0, items[2], item_lens[2], 0, item_lens[2]);

    if (n_items == 4)
        seg_by_did (VB, items[3], item_lens[3], ATTR_Target_STRAND, item_lens[3]);

    // TO DO: use seg_duplicate_last if n_items hasn't changed
    container_seg (vb, CTX (ATTR_Target), (ContainerP)&con, NULL, 0, n_items-1); // account for space separators

    return true;
}

static int exon_number_prediction (VBlockGFFP vb, ContextP ctx)
{
    bool is_exon = str_issame_(STRtxtw(CTX(GFF_TYPE)->last_txt), "exon", 4);
    bool prev_is_transcript = str_issame_(STRa(vb->prev_type), "transcript", 10);

    int prediction = (is_exon && !prev_is_transcript) ? ctx->last_value.i + 1
                   : (is_exon && prev_is_transcript)  ? 1
                   : /* not exon */                     ctx->last_value.i;

    return prediction;
}

// transcript name eg UBE2J2-219 consists of the gene name (UBE2J2) and an integer
static void gff_seg_transcript_name (VBlockGFFP vb, STRp(transcript_name))
{
    ContextP ctx = CTX(ATTR_transcript_name);

    str_split (transcript_name, transcript_name_len, 2, '-', item, true);
    
    if (n_items == 2 && 
        ctx_encountered_in_line (VB, ATTR_gene_name) &&
        str_issame_(STRi(item,0), STRtxtw(CTX(ATTR_gene_name)->last_txt)) &&
        str_is_int (STRi(item,1))) {

        seg_by_ctx (VB, STRa(copy_gene_name_snip), CTX(ATTR_transcript_name_gene), item_lens[0]);
        seg_integer (VB, CTX(ATTR_transcript_name_num), atoi (items[1]), true, item_lens[1]);
        seg_by_ctx (VB, STRa(transcript_name_container_snip), ctx, 1/*separator*/);
    }

    else // fallback
        seg_by_ctx (VB, STRa(transcript_name), ctx, transcript_name_len); 
}

static void gff_seg_exon_number (VBlockGFFP vb, STRp(exon_number))
{
    ContextP ctx = CTX(ATTR_exon_number);

    int64_t en;
    ASSSEG (str_get_int (STRa(exon_number), &en), exon_number, "expecting exon_number=\"%.*s\" to be an integer", STRf(exon_number));

    if (en == exon_number_prediction (vb, ctx))
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, GFF_SPECIAL_exon_number }, 2, ctx, exon_number_len);

    else 
        seg_by_ctx (VB, STRa(exon_number), ctx, exon_number_len);
    
    ctx_set_last_value (VB, ctx, en);
}

SPECIAL_RECONSTRUCTOR (gff_piz_special_exon_number)
{
    new_value->i = exon_number_prediction (VB_GFF, ctx);
    
    if (reconstruct) 
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

static inline DictId gff_seg_attr_subfield (VBlockP vb, STRp(tag), STRp(value))
{
    DictId dict_id = dict_id_make (STRa(tag), DTYPE_GFF_ATTR);

    ContextP ctx = ctx_get_ctx_tag (vb, dict_id, tag, tag_len);
    
    switch (dict_id.num) {

    // ID - sometimes this is a sequential number (GRCh37/38)
    // sometimes it is something like this: c5312581-5d6e-4234-89d7-4974581f2993
    case _ATTR_ID: 
        if (str_is_int (value, value_len))
            seg_pos_field (vb, ATTR_ID, ATTR_ID, SPF_BAD_SNIPS_TOO, 0, STRa(value), 0, value_len);
        else
            seg_by_did (VB, STRa(value), ATTR_ID, value_len);
        break;

    // Dbxref (example: "dbSNP_151:rs1307114892") - we divide to the non-numeric part which we store
    // in a dictionary and the numeric part which store in a local
    case _ATTR_Dbxref:
    case _ATTR_db_xref: 
        seg_id_field (vb, ctx, value, value_len, false); 
        break;

    case _ATTR_Target:
        if (!gff_seg_target (vb, STRa(value)))
            goto plain_seg;
        break;

    // example: Name=gene021872.2
    case _ATTR_Name:
        seg_id_field (vb, CTX(ATTR_Name), value, value_len, false);
        break;

    // example: Parent=mRNA00001,mRNA00002,mRNA00003
    case _ATTR_Parent:
        seg_array (vb, CTX(ATTR_Parent), ATTR_Parent, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len);
        break;

    //case _ATTR_Gap: // I tried: 1. array (no improvement) ; 2. string of op-codes in b250 + integers in local (negligible improvement)

    // subfields that are arrays of structs, for example:
    // "non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
    case _ATTR_Variant_effect: {
        static const MediumContainer Variant_effect = {
            .nitems_lo   = 4, 
            .drop_final_repsep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="V0arEff" }, .separator = {' '} },
                             { .dict_id={.id="V1arEff" }, .separator = {' '} },
                             { .dict_id={.id="V2arEff" }, .separator = {' '} },
                             { .dict_id={.num=_ENSTid  },                    } }
        };
        seg_array_of_struct (vb, CTX(ATTR_Variant_effect), Variant_effect, STRa(value), NULL, value_len);
        break;
    }

    case _ATTR_sift_prediction: {
        static const MediumContainer sift_prediction = {
            .nitems_lo   = 4, 
            .drop_final_repsep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="S0iftPr" }, .separator = {' '} },
                             { .dict_id={.id="S1iftPr" }, .separator = {' '} },
                             { .dict_id={.id="S2iftPr" }, .separator = {' '} },
                             { .dict_id={.num=_ENSTid  },                    } }
        };
        seg_array_of_struct (vb, CTX(ATTR_sift_prediction), sift_prediction, STRa(value), NULL, value_len);
        break;
    }

    case _ATTR_polyphen_prediction: {
        static const MediumContainer polyphen_prediction = {
            .nitems_lo   = 4, 
            .drop_final_repsep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="P0olyPhP" }, .separator = {' '} },
                             { .dict_id={.id="P1olyPhP" }, .separator = {' '} },
                             { .dict_id={.id="P2olyPhP" }, .separator = {' '} },
                             { .dict_id={.num=_ENSTid   },                    } }
        };
        seg_array_of_struct (vb, CTX(ATTR_polyphen_prediction), polyphen_prediction, STRa(value), NULL, value_len);
        break;
    }

    case _ATTR_variant_peptide: {
        static const MediumContainer variant_peptide = {
            .nitems_lo   = 3, 
            .drop_final_repsep = true,
            .repsep      = {','},
            .items       = { { .dict_id={.id="v0arPep"  }, .separator = {' '} }, // small v to differentiate from Variant_effect, so that dict_id to did_i mapper can map both
                             { .dict_id={.id="v1arPep"  }, .separator = {' '} },
                             { .dict_id={.num=_ENSTid   },                    } }
        };
        seg_array_of_struct (vb, CTX(ATTR_variant_peptide), variant_peptide, STRa(value), NULL, value_len);
        break;
    }

    // we store these 3 in one dictionary, as they are correlated and will compress better together
    case _ATTR_Variant_seq:
    case _ATTR_Reference_seq:
    case _ATTR_ancestral_allele: 
        // note: all three are stored together in _ATTR_Reference_seq as they are correlated
        seg_add_to_local_text (vb, CTX(ATTR_Reference_seq), STRa(value), LOOKUP_NONE, value_len); 
        break;

    case _ATTR_chr:
        chrom_seg_by_did_i (vb, ATTR_chr, STRa(value), value_len);
        break;

    // Optimize Variant_freq
    case _ATTR_Variant_freq:
        if (flag.optimize_Vf) {
            unsigned optimized_snip_len = value_len + 20;
            char optimized_snip[optimized_snip_len]; // used for 1. fields that are optimized 2. fields translated luft->primary

            if (optimize_float_2_sig_dig (STRa(value), 0, optimized_snip, &optimized_snip_len)) {        
                vb->recon_size -= (int)(value_len) - (int)optimized_snip_len;
                value = optimized_snip;
                value_len = optimized_snip_len;
            }            
        }
        goto plain_seg; // proceed with adding to dictionary/b250

    // GTF fields
    case _ATTR_gene_name:
        set_last_txt (ctx->did_i, value); // consumed by gff_seg_transcript_name
        goto plain_seg;

    case _ATTR_transcript_name:
        gff_seg_transcript_name (VB_GFF, STRa(value));
        break;

    case _ATTR_gene_id: 
    case _ATTR_transcript_id: 
    case _ATTR_ccds_id: 
    case _ATTR_protein_id: 
    case _ATTR_exon_id: 
        seg_id_field2 (vb, ctx, STRa(value), value_len);
        break;

    case _ATTR_exon_number:
        gff_seg_exon_number (VB_GFF, STRa(value));
        break;

    default:
    plain_seg:
        seg_by_ctx (VB, STRa(value), ctx, value_len);
    }

    ctx_set_encountered (vb, ctx);

    return dict_id;
}

// determine GFF version by ATTR key/value seperator
static void gff_segconf_set_gff_version (VBlockP vb, STRp(attr))
{
    for (int i=0; i < attr_len-1; i++)
        switch (attr[i]) {
            case ' ' : segconf.gff_version = 2; return;
            case '=' : segconf.gff_version = 3; return;
            case '"' : ASSSEG0 (false, attr, "Invalid attribute: quotation mark before separator");
            default  : break;
        }

    ASSSEG0 (false, attr, "Invalid attribute: separator not found");
}

// see: http://gmod.org/wiki/GFF2 and https://www.ensembl.org/info/website/upload/gff.html
static void gff_seg_gff2_attrs_field (VBlockP vb, STRp(attribute))
{
    // case: "." attribute
    if (attribute_len == 1 && *attribute == '.') {
        seg_by_did (VB, ".", 1, GFF_ATTRS, 2);
        return;
    }

    char prefixes[attribute_len + 1];
    prefixes[0] = prefixes[1] = CON_PX_SEP;
    unsigned prefixes_len = 2;

    bool extra_space = (attribute[attribute_len-1] == ' ');
    if (extra_space) attribute_len--;

    // a semicolon contained within a quote is replaced to avoid splitting by it
    bool in_quotes=false;
    for (int i=0; i < attribute_len; i++)
        if (attribute[i] == '"') in_quotes = !in_quotes;
        else if (in_quotes && attribute[i] == ';') ((char*)attribute)[i] = 1;

    str_split (attribute, attribute_len, MAX_DICTS, ';', attr, false);
    ASSSEG (n_attrs, attribute, "Invalid attributes field: %.*s", STRf(attribute));

    // replace back
    str_replace_letter ((char*)STRa(attribute), 1, ';');

    Container con = { .repeats = 1 };

    if (!attr_lens[n_attrs-1]) // last item is ends with a ; - creating a fake final item
        n_attrs--;
    else
        con.drop_final_item_sep = true; // last item doesn't not end with a semicolon

    con_set_nitems (con, n_attrs);

    // The prefix of the container includes the following stuff:
    // 1. tag
    // 2. (sometimes) leading space - space may appear before tag (but not in the first item)
    // 3. (sometimes) final space - space may appear after value
    // 4. (sometimes) 2 quotation marks - value may be enclosed in quotes
    int total_values_len = 0; 
    for (unsigned i=0; i < n_attrs; i++) {

        ASSSEG0 (attr_lens[i] >= 3, attrs[i], "Invalid GFF2 attributes field: Expecting attribute to have at least 3 characters");

        // case: space after semicolon of previous field: “gene_id "ENSG00000223972"; gene_version "5"; ”
        bool leading_space = (i && attrs[i][0] == ' ');
        if (leading_space) { attrs[i]++; attr_lens[i]--; }

        // case: final space before semicolon
        bool final_space = (attrs[i][attr_lens[i] - 1]== ' ');
        if (final_space) attr_lens[i]--;

        // find value - after first space - this is the tag/value separator (there might be additional spaces in a quotes-enclosed value)
        STR0(value); 
        STR0(tag);
        for (rom c = attrs[i]; c < attrs[i] + attr_lens[i] - 1; c++)
            if (*c == ' ') {
                tag       = attrs[i];
                tag_len   = c - tag;
                value     = c + 1;
                value_len = attrs[i] + attr_lens[i] - value;
                break;
            }

        ASSSEG0 (value_len, attrs[i], "Invalid GFF2 attributes field: No space separator");

        // case: value is enclosed in quotes: “gene_id "ENSG00000223972"; ”
        bool has_quotes = (value[0] == '"'); 
        ASSSEG (has_quotes == (value[value_len-1] == '"'), attrs[i], "Invalid GFF2 attributes field: expecting no quotes or a pair of quotes: %.*s", STRfi (attr, i));
        
        if (has_quotes) {  value += 1; value_len -= 2; } // remove quotes

        total_values_len += value_len;

        // we currently support either a final space or a final quote - we've yet to see an example of space following
        // a quote, and it would be tricky to implement as container item separators have only two characters
        ASSSEG0 (!(has_quotes && final_space), attrs[i], "Invalid GFF2 attributes field: a field has a space after the closing quote");
        
        // seg value
        con.items[i].dict_id = gff_seg_attr_subfield (vb, STRa(tag), STRa(value));

        // add leading_space, name, ' ' separator, opening quote to prefix 
        int item_prefix_len = leading_space + tag_len + 1 + has_quotes;
        memcpy (&prefixes[prefixes_len], tag - leading_space, item_prefix_len); 
        prefixes_len += item_prefix_len + 1/*CON_PX_SEP*/; 
        prefixes[prefixes_len-1] = CON_PX_SEP;

        // add (closing quote or final_space) and ';' separator to item separator 
        con.items[i].separator[0] = final_space ? ' ' : has_quotes ? '"' : ';';  
        con.items[i].separator[1] = (final_space || has_quotes) ? ';' : 0;
    }

    // add extra_space to the prefixes
    if (extra_space) {
        prefixes[prefixes_len++] = ' ';
        prefixes[prefixes_len++] = CON_PX_SEP;
    }

    container_seg (vb, CTX(GFF_ATTRS), &con, prefixes, prefixes_len, attribute_len - total_values_len + 1/*\n*/ + extra_space); 
}

static void gff_seg_gff3_attrs_field (VBlockP vb, STRp(field))
{
    // case: "." field
    if (field_len == 1 && *field == '.') {
        seg_by_did (VB, ".", 1, GFF_ATTRS, 2);
        return;
    }

    char prefixes[field_len + 2];
    prefixes[0] = prefixes[1] = CON_PX_SEP;
    unsigned prefixes_len = 2;

    str_split (field, field_len, MAX_DICTS, ';', attr, false);
    ASSSEG (n_attrs, field, "Invalid attributes field: %.*s", STRf(field));

    Container con = { .repeats = 1 };

    if (!attr_lens[n_attrs-1]) // last item is ends with a ; - creating a fake final item
        n_attrs--;
    else
        con.drop_final_item_sep = true; // last item doesn't not end with a semicolon

    con_set_nitems (con, n_attrs);

    for (unsigned i=0; i < n_attrs; i++) {

        str_split (attrs[i], attr_lens[i], 2, '=', tag_val, true);
        ASSSEG (n_tag_vals == 2 && tag_val_lens[0] > 0 && tag_val_lens[1] > 0, attrs[i], "Invalid attribute: %.*s", STRfi (attr, i));

        memcpy (&prefixes[prefixes_len], tag_vals[0], tag_val_lens[0] + 1);
        prefixes_len += tag_val_lens[0] + 2; // including the '=' and CON_PX_SEP
        prefixes[prefixes_len-1] = CON_PX_SEP;

        con.items[i].dict_id = gff_seg_attr_subfield (vb, STRi(tag_val, 0), STRi (tag_val, 1));
        con.items[i].separator[0] = ';'; 
    }

    container_seg (vb, CTX(GFF_ATTRS), &con, prefixes, prefixes_len, 
                   prefixes_len - 1 - con.drop_final_item_sep); // tags inc. = and (; or \n) separator 
}

static void gff_seg_START (VBlockGFFP vb, STRp(start))
{
    seg_pos_field (VB, GFF_START, GFF_START, 0, 0, STRa(start), 0, start_len+1);
    random_access_update_pos (VB, 0, GFF_START);
}

static void gff_seg_TYPE (VBlockGFFP vb, STRp(type))
{
    seg_by_did (VB, STRa(type), GFF_TYPE, type_len + 1);

    set_last_txt (GFF_TYPE, type);
}

rom gff_seg_txt_line (VBlockP vb, rom field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    rom next_field=field_start_line, field_start;
    unsigned field_len=0;
    char separator;
    int32_t len = (int32_t)remaining_txt_len;

    // handle a comment (#), directive (##), experimental directive (#!) or empty line 
    if (next_field[0] == '#' || next_field[0] == '\n' || next_field[0] == '\r' ||
        (remaining_txt_len >= 6 && !memcmp (next_field, "track ", 6))) { // Ensembl GTF track lines, eg: “track name=coords description="Chromosome coordinates list" priority=2"”
        SEG_NEXT_ITEM_NL (GFF_COMMENT);

        // TO DO: support embedded FASTA, bug 392
        ASSINP (!str_issame_ (STRd(GFF_COMMENT), "##FASTA", 7), "%s contains a ##FASTA directive. To compress it use --input=generic", txt_name);
        
        goto eol; // if we have a comment, then during piz, the other fields will be filtered out
    }
    else 
        seg_by_did (VB, NULL, 0, GFF_COMMENT, 0); // missing comment field

    GET_NEXT_ITEM (GFF_SEQID);
    chrom_seg (vb, STRd(GFF_SEQID));

    SEG_NEXT_ITEM (GFF_SOURCE);
    
    GET_NEXT_ITEM (GFF_TYPE);
    gff_seg_TYPE (VB_GFF, STRd(GFF_TYPE));

    GET_NEXT_ITEM (GFF_START);
    gff_seg_START (VB_GFF, STRd(GFF_START));

    GET_NEXT_ITEM (GFF_END);
    seg_pos_field (vb, GFF_END, GFF_START, 0, 0, STRd(GFF_END), 0, GFF_END_len+1);

    SEG_NEXT_ITEM (GFF_SCORE);
    SEG_NEXT_ITEM (GFF_STRAND);
    SEG_MAYBE_LAST_ITEM (GFF_PHASE);

    if (separator != '\n') { // "attributes" is an optional field per http://gmod.org/wiki/GFF3
        GET_LAST_ITEM (GFF_ATTRS); 

        if (segconf.running && !segconf.gff_version)
            gff_segconf_set_gff_version (vb, STRd(GFF_ATTRS));

        if (segconf.gff_version >= 3) gff_seg_gff3_attrs_field (vb, STRd(GFF_ATTRS));
        else                          gff_seg_gff2_attrs_field (vb, STRd(GFF_ATTRS));
    }
    else
        seg_by_did (VB, NULL, 0, GFF_ATTRS, 0); // NULL=MISSING so previous \t is removed

    VB_GFF->prev_type = GFF_TYPE_str;
    VB_GFF->prev_type_len = GFF_TYPE_len;

eol:             
    SEG_EOL (GFF_EOL, false);

    return next_field;
}

//----------
// PIZ stuff
//----------

// filter is called before reconstruction of a repeat or an item, and returns true if item should be reconstructed. if not reconstructed, contexts are not consumed.
CONTAINER_FILTER_FUNC (gff_piz_filter)
{
    if (item < 1 || item == 10) return true; // we filter all items except COMMENT and EOL

    return CTX(GFF_COMMENT)->last_value.i == WORD_INDEX_MISSING; // reconstruct non-comment contexts only if this isn't a comment line
}

CONTAINER_CALLBACK (gff_piz_container_cb)
{
    if (is_top_level) {
        VB_GFF->prev_type = last_txt (vb, GFF_TYPE);
        VB_GFF->prev_type_len = vb->last_txt_len (GFF_TYPE);
    }
}