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
#include "contigs.h"

#define MAX_ENST_ITEMS 10 // maximum number of items in an enst structure. this can be changed without impacting backward compatability.

typedef struct VBlockGFF {    
    VBLOCK_COMMON_FIELDS
    STR(prev_type);

    // current line
    bool has_parent;  
} VBlockGFF;

typedef VBlockGFF *VBlockGFFP;
#define VB_GFF ((VBlockGFFP)vb)

unsigned gff_vb_size (DataType dt) { return sizeof (VBlockGFF);  }

sSTRl(copy_gene_name_snip,32);
sSTRl(dbx_container_snip,100);
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

    con = (SmallContainer){
        .repeats      = 1,
        .nitems_lo    = 2,
        .items        = { { .dict_id = { _ATTR_DBXdb }, .separator = {':'} },  // 0
                          { .dict_id = { _ATTR_DBXid }                     } } // 1
    };
    container_prepare_snip ((ContainerP)&con, 0, 0, qSTRa (dbx_container_snip));
}

// detect if a generic file is actually a GFF3/GVF/GTF (but cannot detect non-GTF GFF2)
bool is_gff (STRp(header), bool *need_more)
{
    SAFE_NULT (header);

    bool is_gff = (str_isprefix_(STRa(header), _S("##gff-version"))) || // must appear for GFF3/GTF
                  (str_isprefix_(STRa(header), _S("## mirGFF3")))    || // mirGFF, a specific type of GFF3
                   strstr (header, "#!genome-build")                 || // normally appears in GTF / GVF
                   strstr (header, "#!genome-version");                 // normally appears in GTF / GVF

    SAFE_RESTORE;
    return is_gff;
}

bool gff_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags)
{
    SAFE_NULB (*txt_header); // for atoi

    // case: gff version announced explicitly. otherwise we will deduce it during segconf
    if (str_isprefix_(STRb(*txt_header), _S("##gff-version")) && txt_header->len32 >= 15)
        segconf.gff_version = atoi (Bc(*txt_header, 14));

    else if (str_isprefix_(STRb(*txt_header), _S("## mirGFF3")))
        segconf.gff_version = 3;

    else if (str_isprefix_(STRb(*txt_header), _S("#gtf-version")))
        segconf.gff_version = 2;

    SAFE_RESTORE;

    return true;
}

// search for last newline, and also search for embedded FASTA
int32_t gff_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    ASSERT (*i >= 0 && *i < Ltxt, "*i=%d is out of range [0,%u]", *i, Ltxt);

    int32_t final_i = *i;
    int32_t last_newline = -1;

    for (int32_t j=first_i; j <= final_i; j++)
        if (*Btxt (j) == '\n') {
            last_newline = j;

            if (j < final_i && *Btxt (j+1) == '>') {
                if (!segconf.running) {
                    // note: we don't run segconf on an embedded FASTA - we set the values here instead
                    segconf.has_embedded_fasta = true;
                    segconf.fasta_has_contigs  = false; // GFF3-embedded FASTA doesn't have contigs, because did=0 is reserved for GFF's SEQID
                    segconf.seq_type           = SQT_NUKE;
                }

                break; // terminate VB (and component) at this newline - next line is FASTA
            }
        }

    if (last_newline != -1) {
        *i = last_newline;
        return Ltxt-1 - last_newline;
    }

    else { // no newline found
        *i = (int32_t)first_i - 1;
        return -1; // cannot find \n in the data starting first_i
    }
}

// called from seg_all_data_lines
void gff_seg_initialize (VBlockP vb)
{
    ctx_set_store (VB, STORE_INT, GFF_START/*since v12*/, GFF_END/*14.0.18*/, ATTR_Target_POS, ATTR_ID, ATTR_exon_number, DID_EOL);

    ctx_set_store (VB, STORE_INDEX, GFF_SEQID/*v12*/, GFF_COMMENT/*12.0.12*/, DID_EOL);

    ctx_set_no_stons (vb, GFF_COMMENT/*required by STORE_INDEX (otherwise singletons don't get their index stored)*/, 
                      GFF_SEQID/*needs b250 node_index for random access*/, 
                      GFF_ATTRS, 
                      ATTR_Target_POS, ATTR_ID, GFF_START, GFF_END, DID_EOL); // as requied by seg_pos_field

    ctx_set_ltype (vb, LT_STRING, ATTR_Variant_seq, ATTR_Reference_seq, ATTR_ancestral_allele, DID_EOL);

    ctx_set_ltype (vb, LT_DYN_INT, ATTR_transcript_name_num, DID_EOL);

    ctx_consolidate_stats (vb, ATTR_transcript_name, ATTR_transcript_name_gene, ATTR_transcript_name_num, DID_EOL);
    ctx_consolidate_stats (vb, ATTR_Target, ATTR_Target_ID, ATTR_Target_POS, ATTR_Target_STRAND, DID_EOL);

    if (segconf.has[ATTR_Dbxref])
        ctx_consolidate_stats (vb, ATTR_Dbxref, ATTR_DBXid, ATTR_DBXdb, DID_EOL);
    else if (segconf.has[ATTR_db_xref])
        ctx_consolidate_stats (vb, ATTR_db_xref, ATTR_DBXid, ATTR_DBXdb, DID_EOL);
}

static void gff_seg_finalize_segconf (VBlockP vb)
{
    if (segconf.has[ATTR_score] && segconf.has[ATTR_cscore] && segconf.has[ATTR_sscore])
        stats_add_one_program (_S("Prodigal"));
}

void gff_seg_finalize (VBlockP vb)
{
    if (segconf.running) 
        gff_seg_finalize_segconf (vb);

    CTX(ENSTid)->st_did_i = DID_NONE; // cancel consolidatation as it goes into multiple attributes

    // top level snip
    SmallContainer top_level = { 
        .repeats      = vb->lines.len32,
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

// initialization of the line
void gff_reset_line (VBlockP vb_)
{
    VBlockGFFP vb = (VBlockGFFP)vb_;

    ASSERT (VB_DT(GFF), "VB has wrong data type: %s", dt_name (vb->data_type));

    if (IS_ZIP) {
        vb->has_parent = false;
    }
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

bool gff_seg_is_big (ConstVBlockP vb, DictId dict_id, DictId st_dict_id)
{
    return dict_id.num == _ATTR_DBXid           ||
           dict_id.num == _ATTR_transcript_name ||
           dict_id.num == _ATTR_gene_name;
}

// eg: "DGVa_201810:essv9162", "GeneID:100287102,HGNC:HGNC:37102"
static bool gff_seg_dbxref_one (VBlockP vb, ContextP ctx, STRp(value), uint32_t repeat)
{
    // if this file has no "Parent", we're better off segging as a simple ID
    if (!segconf.has[ATTR_Parent] || segconf.running)
        seg_id_field (vb, ctx, STRa(value), false, value_len);  

    else {
        // split by the first ':'. We don't use str_split bc the second item can contain : which we don't split
        rom colon = memchr (value, ':', value_len);
        if (colon) {
            int db_len = colon - value;
            int id_len = value_len - db_len - 1;

            // container
            seg_by_ctx (vb, STRa(dbx_container_snip), ctx, 1); // account for ':' 
            
            // DB name component
            seg_by_did (vb, value, db_len, ATTR_DBXdb, db_len); 

            // ID component
            seg_id_field (vb, CTX(ATTR_DBXid), colon+1, id_len, false, id_len); 
        }

        else  // fallback
            seg_by_ctx (vb, STRa(value), ctx, value_len); 
    }

    return true;
}

static void gff_seg_dbxref (VBlockGFFP vb, ContextP ctx, STRp(value))
{
    STRlast(last, ctx->did_i);

    // e.g. exons usually have the same DBXref as their transcript 
    if (str_issame(value, last)) 
        seg_by_ctx (VB, (char[]){ SNIP_COPY }, 1, ctx, value_len);

    else
        seg_array_by_callback (VB, ctx, STRa(value), ',', gff_seg_dbxref_one, value_len);

    seg_set_last_txt (VB, ctx, STRa(value));
}

// returns true if successful
static bool gff_seg_target (VBlockGFFP vb, STRp(value))
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
    seg_pos_field (VB, ATTR_Target_POS, ATTR_Target_POS, 0, 0, items[1], item_lens[1], 0, item_lens[1]);
    seg_pos_field (VB, ATTR_Target_POS, ATTR_Target_POS, 0, 0, items[2], item_lens[2], 0, item_lens[2]);

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
    decl_ctx (ATTR_transcript_name);

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
    decl_ctx (ATTR_exon_number);

    int64_t en;
    ASSSEG (str_get_int (STRa(exon_number), &en), "expecting exon_number=\"%.*s\" to be an integer", STRf(exon_number));

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

static inline DictId gff_seg_attr_subfield (VBlockGFFP vb, STRp(tag), STRp(value))
{
    DictId dict_id = dict_id_make (STRa(tag), DTYPE_GFF_ATTR);

    ContextP ctx = ctx_get_ctx_tag (vb, dict_id, tag, tag_len);
    
    if (segconf.running)
        segconf.has[ctx_get_ctx (VB, dict_id)->did_i]++;

    switch (dict_id.num) {

    // ID - sometimes this is a sequential number (GRCh37/38)
    // sometimes it is something like this: c5312581-5d6e-4234-89d7-4974581f2993, sometimes exon-NR_024540.1-2
    case _ATTR_ID: 
        if (str_is_int (value, value_len))
            seg_pos_field (VB, ATTR_ID, ATTR_ID, SPF_BAD_SNIPS_TOO, 0, STRa(value), 0, value_len);
        else
            seg_by_did (VB, STRa(value), ATTR_ID, value_len);
        break;

    // example: "GeneID:100287102,GenBank:NR_046018.2,HGNC:HGNC:37102"
    case _ATTR_Dbxref:
    case _ATTR_db_xref: 
        gff_seg_dbxref (vb, ctx, STRa(value));
        break;

    case _ATTR_Target:
        if (!gff_seg_target (vb, STRa(value)))
            goto plain_seg;
        break;

    // example: { "gene021872.2" } { "HG-1885", "HG-789" }
    case _ATTR_Name:
        seg_id_field (VB, CTX(ATTR_Name), STRa(value), false, value_len);
        break;

    // example: Parent=mRNA00001,mRNA00002,mRNA00003
    case _ATTR_Parent:
        seg_array (VB, CTX(ATTR_Parent), ATTR_Parent, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len);
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
        seg_array_of_struct (VB, CTX(ATTR_Variant_effect), Variant_effect, STRa(value), NULL, value_len);
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
        seg_array_of_struct (VB, CTX(ATTR_sift_prediction), sift_prediction, STRa(value), NULL, value_len);
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
        seg_array_of_struct (VB, CTX(ATTR_polyphen_prediction), polyphen_prediction, STRa(value), NULL, value_len);
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
        seg_array_of_struct (VB, CTX(ATTR_variant_peptide), variant_peptide, STRa(value), NULL, value_len);
        break;
    }

    // we store these 3 in one dictionary, as they are correlated and will compress better together
    case _ATTR_Variant_seq:
    case _ATTR_Reference_seq:
    case _ATTR_ancestral_allele: 
        // note: all three are stored together in _ATTR_Reference_seq as they are correlated
        seg_add_to_local_string (VB, CTX(ATTR_Reference_seq), STRa(value), LOOKUP_NONE, value_len); 
        break;

    case _ATTR_chr:
        chrom_seg_by_did_i (VB, ATTR_chr, STRa(value), value_len);
        break;

    // Optimize Variant_freq
    case _ATTR_Variant_freq:
        if (flag.optimize_Vf) {
            STRli (optimized_snip, value_len + 20); // used for 1. fields that are optimized 2. fields translated luft->primary

            if (optimize_float_2_sig_dig (STRa(value), 0, qSTRa(optimized_snip))) {        
                vb->recon_size -= (int)(value_len) - (int)optimized_snip_len;
                STRset (value, optimized_snip);
            }            
        }
        goto plain_seg; // proceed with adding to dictionary/b250

    // GTF fields
    case _ATTR_gene_name:
        set_last_txt (ctx->did_i, value); // consumed by gff_seg_transcript_name
        goto plain_seg;

    case _ATTR_transcript_name:
        gff_seg_transcript_name (vb, STRa(value));
        break;

    case _ATTR_gene_id:
    case _ATTR_havana_gene: 
    case _ATTR_transcript_id: 
    case _ATTR_havana_transcript: 
    case _ATTR_ccds_id: 
    case _ATTR_ccdsid: 
    case _ATTR_protein_id: 
    case _ATTR_exon_id: 
        seg_id_field (VB, ctx, STRa(value), true, value_len);
        break;

    case _ATTR_exon_number:
        gff_seg_exon_number (vb, STRa(value));
        break;

    default:
    plain_seg:
        seg_by_ctx (VB, STRa(value), ctx, value_len);
    }

    ctx_set_encountered (VB, ctx);

    return dict_id;
}

// determine GFF version by ATTR key/value seperator
static void gff_segconf_set_gff_version (VBlockGFFP vb, STRp(attr))
{
    for (int i=0; i < attr_len-1; i++)
        switch (attr[i]) {
            case ' ' : segconf.gff_version = 2; return;
            case '=' : segconf.gff_version = 3; return;
            case '"' : ABOSEG0 ("Invalid attribute: quotation mark before separator");
            default  : break;
        }

    ABOSEG0 ("Invalid attribute: separator not found");
}

// see: http://gmod.org/wiki/GFF2 and https://www.ensembl.org/info/website/upload/gff.html
static void gff_seg_gff2_attrs_field (VBlockGFFP vb, STRp(attribute))
{
    // case: "." attribute
    if (IS_PERIOD(attribute)) {
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
    ASSSEG (n_attrs, "Invalid attributes field: %.*s", STRf(attribute));

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

        ASSSEG0 (attr_lens[i] >= 3, "Invalid GFF2 attributes field: Expecting attribute to have at least 3 characters");

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

        ASSSEG0 (value_len, "Invalid GFF2 attributes field: No space separator");

        // case: value is enclosed in quotes: “gene_id "ENSG00000223972"; ”
        bool has_quotes = (value[0] == '"'); 
        ASSSEG (has_quotes == (value[value_len-1] == '"'), "Invalid GFF2 attributes field: expecting no quotes or a pair of quotes: %.*s", STRfi (attr, i));
        
        if (has_quotes) {  value += 1; value_len -= 2; } // remove quotes

        total_values_len += value_len;

        // we currently support either a final space or a final quote - we've yet to see an example of space following
        // a quote, and it would be tricky to implement as container item separators have only two characters
        ASSSEG0 (!(has_quotes && final_space), "Invalid GFF2 attributes field: a field has a space after the closing quote");
        
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

static void gff_seg_gff3_attrs_field (VBlockGFFP vb, STRp(field))
{
    // case: "." field
    if (IS_PERIOD (field)) {
        seg_by_did (VB, ".", 1, GFF_ATTRS, 2);
        return;
    }

    char prefixes[field_len + 2];
    prefixes[0] = prefixes[1] = CON_PX_SEP;
    unsigned prefixes_len = 2;

    str_split (field, field_len, MAX_DICTS, ';', attr, false);
    ASSSEG (n_attrs, "Invalid attributes field: %.*s", STRf(field));

    Container con = { .repeats = 1 };

    if (!attr_lens[n_attrs-1]) // last item is ends with a ; - creating a fake final item
        n_attrs--;
    else
        con.drop_final_item_sep = true; // last item doesn't not end with a semicolon

    // check for parent - hacky for now
    for (int i=n_attrs-1; i >= 0; i--) // Parent is usually towards the end
        if (str_isprefix_(STRi(attr,i), _S("Parent="))) {
            vb->has_parent = true;
            break;
        }

    con_set_nitems (con, n_attrs);

    for (unsigned i=0; i < n_attrs; i++) {

        str_split (attrs[i], attr_lens[i], 2, '=', tag_val, true);
        ASSSEG (n_tag_vals == 2 && tag_val_lens[0] > 0 && tag_val_lens[1] > 0, "Invalid attribute: %.*s", STRfi (attr, i));

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

rom gff_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    VBlockGFFP vb = (VBlockGFFP)vb_;

    rom next_field=field_start_line, field_start;
    unsigned field_len=0;
    char separator;
    int32_t len = (int32_t)remaining_txt_len;

    // handle a comment (#), directive (##), experimental directive (#!) or empty line 
    if (next_field[0] == '#' || next_field[0] == '\n' || next_field[0] == '\r' ||
        (remaining_txt_len >= 6 && !memcmp (next_field, "track ", 6))) { // Ensembl GTF track lines, eg: “track name=coords description="Chromosome coordinates list" priority=2"”
        GET_NEXT_ITEM_NL (GFF_COMMENT); 
        seg_by_did (VB, field_start, field_len, GFF_COMMENT, field_len+1);
        
        goto eol; // if we have a comment, then during piz, the other fields will be filtered out
    }
    else 
        seg_by_did (VB, NULL, 0, GFF_COMMENT, 0); // missing comment field

    GET_NEXT_ITEM (GFF_SEQID);
    chrom_seg (VB, STRd(GFF_SEQID));

    SEG_NEXT_ITEM (GFF_SOURCE);
    
    GET_NEXT_ITEM (GFF_TYPE);
    gff_seg_TYPE (vb, STRd(GFF_TYPE));

    GET_NEXT_ITEM (GFF_START);
    gff_seg_START (vb, STRd(GFF_START));

    GET_NEXT_ITEM (GFF_END);
    seg_pos_field (VB, GFF_END, GFF_START, 0, 0, STRd(GFF_END), 0, GFF_END_len+1);

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

bool gff_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header, uint32_t *txt_data_so_far_single_0_increment)
{
    if (VER(15) && vb->flags.gff.embedded_fasta) 
        vb->data_type = DT_FASTA;

    return true;
}

// filter is called before reconstruction of a repeat or an item, and returns true if item should be reconstructed. if not reconstructed, contexts are not consumed.
CONTAINER_FILTER_FUNC (gff_piz_filter)
{
    if (item < 1 || item == 10) return true; // we filter all items except COMMENT and EOL

    return CTX(GFF_COMMENT)->last_value.i == WORD_INDEX_MISSING; // reconstruct non-comment contexts only if this isn't a comment line
}

CONTAINER_CALLBACK (gff_piz_container_cb)
{
    if (is_top_level) 
        SETlast (VB_GFF->prev_type, GFF_TYPE);
}