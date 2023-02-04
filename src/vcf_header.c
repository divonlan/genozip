// ------------------------------------------------------------------
//   header_vcf.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "zfile.h"
#include "txtfile.h"
#include "vblock.h"
#include "crypt.h"
#include "version.h"
#include "endianness.h"
#include "file.h"
#include "dispatcher.h"
#include "txtfile.h"
#include "txtheader.h"
#include "strings.h"
#include "dict_id.h"
#include "website.h"
#include "contigs.h"
#include "gencomp.h"
#include "stats.h"
#include "tip.h"

// Globals
static VcfVersion vcf_version;
uint32_t vcf_num_samples = 0; // number of samples in the file
static uint32_t vcf_num_displayed_samples = 0; // PIZ only: number of samples to be displayed - might be less that vcf_num_samples if --samples is used
static Buffer vcf_field_name_line = EMPTY_BUFFER;  // header line of first VCF file read - use to compare to subsequent files to make sure they have the same header during bound

static bool vcf_has_file_format_vcf = false; 
bool vcf_header_get_has_fileformat (void) { return vcf_has_file_format_vcf; }

rom vcf_header_rename_attrs[NUM_RENAME_ATTRS] = { HK_RENAME_REFALT_ATTR, HK_RENAME_STRAND_ATTR, HK_RENAME_TLAFER_ATTR, HK_RENAME_ALWAYS_ATTR };
const unsigned vcf_header_rename_attr_lens[NUM_RENAME_ATTRS] = { sizeof HK_RENAME_REFALT_ATTR-1, sizeof HK_RENAME_STRAND_ATTR-1, sizeof HK_RENAME_TLAFER_ATTR-1, sizeof HK_RENAME_ALWAYS_ATTR-1 };

static Buffer vcf_header_liftover_dst_contigs = {};

#define FI_LEN (is_info ? 7 : 9) // length of ##FORMAT= or ##INFO=

// referring to sample strings from the --samples command line option
static Buffer cmd_samples_buf = EMPTY_BUFFER; // an array of (char *)
static bool cmd_is_negative_samples = false;

// referring to samples in the vcf file
char *vcf_samples_is_included;         // a bytemap indicating for each sample if it is included
static char **vcf_sample_names;        // an array of char * to nul-terminated names of samples 
static char *vcf_sample_names_data;    // vcf_sample_names point into here

#define LINEIS(s) (line_len > (sizeof s - 1) && !memcmp (line, (s), sizeof s - 1)) 

#define SUBST_LABEL(old,new) ({ buf_add_more (NULL, new_txt_header, (new), (sizeof new - 1), NULL); \
                                buf_add_more (NULL, new_txt_header, line + (sizeof old - 1), line_len - (sizeof old - 1), NULL); })

typedef struct { STR (key); STR (value); } Attr;

void vcf_header_piz_init (void)
{
    vcf_num_samples = 0;
    buf_free (vcf_field_name_line);

    // note: we don't re-initialize vcf_num_displayed_samples - this is calculated only once
}

static void vcf_header_subset_samples (BufferP vcf_header);
static bool vcf_header_set_globals (rom filename, BufferP vcf_header, FailType soft_fail);
static void vcf_header_trim_field_name_line (BufferP vcf_header);

// returns the length of the first line, if it starts with ##fileformat, or 0
// note: per VCF spec, the ##fileformat line must be the first in the file
static inline unsigned vcf_header_parse_fileformat_line (ConstBufferP txt_header)
{
    ARRAY (char, line, *txt_header);
    if (!LINEIS ("##fileformat")) return 0;

    rom newline = memchr (line, '\n', line_len);
    unsigned len = newline ? newline - line + 1 : 0;

    if      (!memcmp (line, "##fileformat=VCFv4.1", MIN_(20, len))) vcf_version = VCF_v4_1;
    else if (!memcmp (line, "##fileformat=VCFv4.2", MIN_(20, len))) vcf_version = VCF_v4_2;
    else if (!memcmp (line, "##fileformat=VCFv4.3", MIN_(20, len))) vcf_version = VCF_v4_3;
    else if (!memcmp (line, "##fileformat=VCFv4.4", MIN_(20, len))) vcf_version = VCF_v4_4;
    else if (!memcmp (line, "##fileformat=VCFv4.5", MIN_(20, len))) vcf_version = VCF_v4_5;

    return len;
}

static inline bool is_field_name_line (STRp(line))
{
    // compare to "#CHROM" - not to entire VCF_FIELD_NAMES as it may be a line maybe separated by spaces instead of tabs
    return LINEIS ("#CHROM");
}

// returns the length of the last line, if it starts with #CHROM, or 0
static bool vcf_header_get_last_line_cb (STRp(line), void *unused1, void *unused2, unsigned unused3)
{ return true; }
static unsigned vcf_header_get_last_line (BufferP txt_header, char **line_p)
{
    int64_t line_len;
    char *line = txtfile_foreach_line (txt_header, true, vcf_header_get_last_line_cb, 0, 0, 0, &line_len);
    if (line_p) *line_p = line;

    return is_field_name_line (STRa(line)) ? (unsigned)line_len : 0;
}

static void vcf_header_add_genozip_command (VBlockP txt_header_vb, BufferP txt_header)
{
    // the command line length is unbound, careful not to put it in a bufprintf
    buf_append_string (txt_header_vb, txt_header, HK_GENOZIP_CMD"\"");
    buf_append_string (txt_header_vb, txt_header, flags_command_line());
    bufprintf (txt_header_vb, txt_header, "\" %s\n", str_time().s);
}

static void vcf_header_get_attribute (STRp(line), unsigned key_len, STRp(attr), bool remove_quotes, bool enforce,// in
                                      rom *snip, unsigned *snip_len) // out
{
    SAFE_NUL (&line[line_len]);
    rom start = strstr (line + key_len, attr);
    SAFE_RESTORE;

    if (!start) goto missing; // attr doesn't exist or malformed line

    start += attr_len;

    rom comma   = memchr (start, ',', &line[line_len] - start);
    rom bracket = memchr (start, '>', &line[line_len] - start);

    if (!comma && !bracket) goto missing; // attr isn't terminated with a comma or >
    rom after = !comma ? bracket : !bracket ? comma : MIN_(comma,bracket);
    
    *snip_len = after - start ;
    *snip = *snip_len ? start : NULL;

    if (remove_quotes && *snip) {
        ASSINP ((*snip)[0]=='"' && (*snip)[*snip_len-1]=='"', "Expecting value an attribute %.*s to be enclosed in quotes. Line=\"%.*s\"", STRf(attr), STRf(line));
        (*snip)++;
        (*snip_len) -= 2;
    }
    return;

missing:
    *snip = NULL;
    *snip_len = 0;
    ASSINP (!enforce, "missing attr %s in header line %.*s", attr, line_len, line);
}

// ZIP: VCF main component (rejects components txt headers don't have contigs)
static void vcf_header_consume_contig (STRp (contig_name), PosType *LN, bool is_luft_contig)
{
    WordIndex ref_index = WORD_INDEX_NONE;

    // case: we have a reference - verify length
    if (!is_luft_contig && (flag.reference & REF_ZIP_LOADED))
        ref_index = ref_contigs_ref_chrom_from_header_chrom (IS_REF_LIFTOVER ? prim_ref : gref, 
                                                             STRa(contig_name), LN);

    // if sorting - we have to pre-populate header & reference contigs (in vcf_zip_initialize) so that vb->is_unsorted is calculated correctly 
    // in vcf_seg_evidence_of_unsorted() (if VBs have their own chrom node_index's, they might be in reverse order vs the eventual z_file chroms)
    if (vcf_is_sorting (VCF_COMP_MAIN)) 
        ctx_populate_zf_ctx (is_luft_contig ? VCF_oCHROM : VCF_CHROM, STRa(contig_name), ref_index, false); 

    if (flag.show_txt_contigs) 
        iprintf ("%s \"%.*s\" LN=%"PRId64" ref_index=%d\n", 
                 is_luft_contig ? "is_luft_contig: " : "", STRf(contig_name), *LN, ref_index);
}

static bool vcf_header_get_dual_coords (STRp(line), void *unused1, void *unused2, unsigned unused3)
{
    if (LINEIS (HK_DC_PRIMARY)) {
        ASSERT (!chain_is_loaded, "--chain cannot be used with %s - it is already a dual-coordinates VCF file - it contains \""HK_DC_PRIMARY"\" in its header", txt_name);
        txt_file->coords   = DC_PRIMARY;
        flag.bind          = BIND_DVCF;
    }

    else if (LINEIS (HK_DC_LUFT)) {
        ASSERT (!chain_is_loaded, "--chain cannot be used with %s - it is already a dual-coordinates VCF file - it contains \""HK_DC_LUFT"\" in its header", txt_name);
        txt_file->coords   = DC_LUFT;
        flag.bind          = BIND_DVCF;
        flag.data_modified = true; // we will be zipping this as a dual-coordinates VCF with default reconstruction as PRIMARY - this turns off digest
        z_file->digest_ctx = DIGEST_CONTEXT_NONE;
    }

    else return false; // continue iterating

    // --test doesn't work with DVCF
    ASSINP0 (!flag.explicit_test, "--test is not supported for dual-coordinates files");
    flag.test = false;

    return true; // found - no more iterating
}

// PIZ with --luft - one line: remove fileformat, update *contig, *reference, dual_coordinates keys to Luft format
static bool vcf_header_piz_liftover_header (STRp(line), void *new_txt_header_, void *unused1, unsigned unused2)
{
    BufferP new_txt_header = (BufferP )new_txt_header_;

    if      (LINEIS (HK_LUFT_CONTIG)) SUBST_LABEL (HK_LUFT_CONTIG, "##contig=");
    else if (LINEIS ("##contig="))    SUBST_LABEL ("##contig=", HK_PRIM_CONTIG);
    else if (LINEIS (HK_LUFT_REF))    SUBST_LABEL (HK_LUFT_REF, "##reference=");
    else if (LINEIS ("##reference=")) SUBST_LABEL ("##reference=", HK_PRIM_REF);
    else if (LINEIS (HK_DC_PRIMARY))  SUBST_LABEL (HK_DC_PRIMARY, HK_DC_LUFT);
    else                              buf_add_more (NULL, new_txt_header, line, line_len, NULL);

    return false; // continue iterating
}

// used with --single-coord: PIZ remove DVCF related lines
static bool vcf_header_piz_single_coord (STRp(line), void *new_txt_header_, void *unused1, unsigned unused2)
{
    BufferP new_txt_header = (BufferP )new_txt_header_;

    if (!LINEIS (HK_LUFT_CONTIG) && !LINEIS (HK_PRIM_CONTIG) &&
        !LINEIS (HK_LUFT_REF) && !LINEIS (HK_PRIM_REF) &&
        !LINEIS (HK_DC) && !LINEIS (HK_CHAIN) && !LINEIS (HK_ORIGINAL_REF) &&
        !LINEIS (KH_INFO INFO_LUFT_NAME) && !LINEIS (KH_INFO INFO_PRIM_NAME) && !LINEIS (KH_INFO INFO_LREJ_NAME) && !LINEIS (KH_INFO INFO_PREJ_NAME))
    
        buf_add_more (NULL, new_txt_header, line, line_len, NULL);

    return false; // continue iterating
}

// PIZ with --luft: remove fileformat, update *contig, *reference, dual_coordinates keys to Luft format
static void vcf_header_rewrite_header (VBlockP txt_header_vb, BufferP txt_header, TxtIteratorCallback callback, void *ret)
{
    #define new_txt_header txt_header_vb->codec_bufs[0]

    buf_alloc (txt_header_vb, &new_txt_header, 0, txt_header->len, char, 0, "codec_bufs[0]"); // initial allocation (might be a bit bigger due to label changes)

    txtfile_foreach_line (txt_header, false, callback, &new_txt_header, ret, 0, 0);

    // replace txt_header with lifted back one
    buf_free (*txt_header);
    buf_copy (txt_header_vb, txt_header, &new_txt_header, char, 0, 0, "txt_data");
    buf_free (new_txt_header);
    #undef new_txt_header
}

// squeeze in the liftover rejects from the header, with their label removed, before the existing unconsumed text
// (data that was read beyond the header) 
static void vcf_header_move_rejects_to_unconsumed_text (BufferP rejects)
{
    if (rejects->len) {
        buf_insert (evb, txt_file->unconsumed_txt, char, 0, rejects->data, rejects->len, "txt_file->unconsumed_txt");
        txt_file->reject_bytes = rejects->len;
    }

    buf_free (*rejects);
}

// ZIP of Luft line: update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
static bool vcf_header_zip_liftback_header_one_line (STRp(line), 
                                                     void *new_txt_header_, void *rejects_, unsigned unused)
{
    BufferP new_txt_header = (BufferP )new_txt_header_;
    
    if      (LINEIS (HK_PRIM_ONLY))   buf_add_more (evb, (BufferP )rejects_, line + (sizeof HK_PRIM_ONLY -1), line_len - (sizeof HK_PRIM_ONLY -1), NULL); // add ##primary_only variants to "rejects" buffer
    else if (LINEIS (HK_PRIM_CONTIG)) SUBST_LABEL (HK_PRIM_CONTIG, "##contig=");
    else if (LINEIS ("##contig="))    SUBST_LABEL ("##contig=", HK_LUFT_CONTIG);
    else if (LINEIS (HK_PRIM_REF))    SUBST_LABEL (HK_PRIM_REF, "##reference=");
    else if (LINEIS ("##reference=")) SUBST_LABEL ("##reference=", HK_LUFT_REF);
    else if (LINEIS (HK_DC_LUFT))     SUBST_LABEL (HK_DC_LUFT, HK_DC_PRIMARY);
    else                              buf_add_more (evb, new_txt_header, line, line_len, NULL);

    return false; // continue iterating
}

// ZIP of LUFT-coords file: update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
static void vcf_header_zip_liftback_header (BufferP txt_header)
{
    #define new_txt_header evb->codec_bufs[0]
    #define rejects        evb->codec_bufs[1]

    // case: Luft file - scan the header and for LIFTBACK_REJECT lines, and move them to txt_file->unconsumed (after removing the label) 
    // to be processsed as normal variants
    buf_alloc (evb, &new_txt_header, 0, txt_header->len + 16384, char, 0, "evb->codec_bufs[0]"); // initial allocation (might be a bit bigger due to label changes)
    buf_alloc (evb, &rejects, 0, 65536, char, 0, "evb->codec_bufs[1]"); // initial allocation

    txtfile_foreach_line (txt_header, false, vcf_header_zip_liftback_header_one_line, &new_txt_header, &rejects, 0, 0);

    // replace txt_header with lifted back one
    buf_free (*txt_header);
    buf_copy (evb, txt_header, &new_txt_header, char, 0, 0, "txt_data");

    // add rejects to the beginning of unconsumed_txt
    vcf_header_move_rejects_to_unconsumed_text (&rejects);

    buf_free (new_txt_header); 

    #undef new_txt_header
    #undef rejects
}

static bool vcf_header_zip_get_luft_only_lines_do (STRp(line), void *rejects_, void *unused1, unsigned unused2)
{
    // add ##luft_only variants to "rejects" buffer
    if (LINEIS (HK_LUFT_ONLY)) buf_add_more (evb, (BufferP )rejects_, line + sizeof HK_LUFT_ONLY -1, line_len - sizeof HK_LUFT_ONLY -1, NULL); 
    return false; // continue iterating
}

// ZIP of PRIMARY-coords file: update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
static void vcf_header_zip_get_luft_only_lines (BufferP txt_header)
{
    #define rejects evb->codec_bufs[1]

    buf_alloc (evb, &rejects, 0, 65536, char, 0, "evb->codec_bufs[1]"); // initial allocation

    txtfile_foreach_line (txt_header, false, vcf_header_zip_get_luft_only_lines_do, &rejects, 0, 0, 0);

    // add rejects to the beginning of unconsumed_txt
    vcf_header_move_rejects_to_unconsumed_text (&rejects);

    #undef rejects
}

// add a key=value to the end of an existing VCF header line
static void vcf_header_add_line_add_attrs (BufferP new_header, STRp(line), Attr *attrs, unsigned num_attr)
{
    bool has_nl, has_13, has_bt;
    line_len -= (has_nl = (line_len >= 1 && line[line_len-1] == '\n'));
    line_len -= (has_13 = (line_len >= 1 && line[line_len-1] == '\r'));
    line_len -= (has_bt = (line_len >= 1 && line[line_len-1] == '>'));
    
    ASSERT (has_nl && has_bt, "While adding attributes: malformed VCF meta-information line: \"%.*s\"", line_len, line);

    buf_add_more (NULL, new_header, line, line_len, NULL);

    for (unsigned i=0; i < num_attr; i++) {
        if (!attrs[i].value_len) attrs[i].value_len = strlen (attrs[i].value);
        bufprintf (evb, new_header, ",%.*s\"%.*s\"", STRf(attrs[i].key), STRf(attrs[i].value));
    }

    BNXTc (*new_header) = '>';
    if (has_13) BNXTc (*new_header) = '\r';
    BNXTc (*new_header) = '\n';
}

// remove an attribute from a line
static void vcf_header_remove_attribute (char *line, unsigned *line_len, STRp(key), unsigned value_len)
{
    SAFE_NUL (&line[*line_len]);
    char *start = strstr (line, key) -1; // -1 to remove the preceding comma
    SAFE_RESTORE;

    unsigned attr_len = 1 /* comma */ + key_len + value_len;
    memmove (start, start + attr_len, &line[*line_len] - (start + attr_len)); 
    
    *line_len -= attr_len;
}

static void vcf_header_zip_convert_header_to_dc_add_lines (BufferP new_header)
{
    // sanity, these should never happen
    ASSERT0 (ref_get_filename (gref), "can't find Luft reference file name - looks like reference was not loaded properly");
    ASSERT0 (ref_get_filename (prim_ref), "can't find Primary reference file name - looks like reference was not loaded properly");
    ASSERT0 (chain_filename, "can't find chain file name - looks like chain file was not loaded properly");

    bufprintf (evb, new_header, "%s\n", HK_DC_PRIMARY);
    bufprintf (evb, new_header, HK_CHAIN"file://%s\n", chain_filename);
    bufprintf (evb, new_header, HK_LUFT_REF"file://%s\n", ref_get_filename (gref));
    bufprintf (evb, new_header, "##reference=file://%s\n", ref_get_filename (prim_ref));

    bufprintf (evb, new_header, KH_INFO_LUFT"\n", GENOZIP_CODE_VERSION);
    bufprintf (evb, new_header, KH_INFO_PRIM"\n", GENOZIP_CODE_VERSION);
    bufprintf (evb, new_header, KH_INFO_LREJ"\n", GENOZIP_CODE_VERSION);
    bufprintf (evb, new_header, KH_INFO_PREJ"\n", GENOZIP_CODE_VERSION);

    // add all Luft contigs (note: we get them from liftover directly, as they zip_initialize that adds them to zctx has not been called yet)
    ARRAY (WordIndex, dst_contig, vcf_header_liftover_dst_contigs); // dst_contig is sorted, but might contain duplicate elements
    for (unsigned i=0; i < dst_contig_len; i++) {
        if (i && dst_contig[i] == dst_contig[i-1]) continue; // skip duplicate

        PosType length;
        rom luft_contig = chain_get_luft_contig (dst_contig[i], &length);

        if (length) bufprintf (evb, new_header, HK_LUFT_CONTIG"<ID=%s,length=%"PRId64">\n", luft_contig, length);
        else        bufprintf (evb, new_header, HK_LUFT_CONTIG"<ID=%s>\n", luft_contig);
    }
    buf_free (vcf_header_liftover_dst_contigs);
}

// Update zctx and tags after reading or creating an INFO / FORMAT RendAlg=value
static TranslatorId vcf_header_set_luft_trans (STRp(id), bool is_info, TranslatorId luft_trans, STRp(number))
{
    DictId dict_id = dict_id_make (id, id_len, is_info ? DTYPE_VCF_INFO : DTYPE_VCF_FORMAT);

    if (luft_trans == TRANS_ID_UNKNOWN) 
        luft_trans = (number_len==1) ? vcf_lo_luft_trans_id (dict_id, *number) : TRANS_ID_NONE;

    ctx_add_new_zf_ctx_from_txtheader (id, id_len, dict_id, luft_trans);
    
    return luft_trans;
}

// If Liftover doesn't exist - adds it according to Genozip's capabilities and also sets ctx->trans_luft
// If it does exist - verifies that it is recognized by Genozip, and sets ctx->trans_luft accordingly
static void vcf_header_update_INFO_FORMAT (STRp(line), bool is_info, BufferP new_header)
{
    STR (number); STR(type); STR (id); STR (rendalg);
    Attr attrs[NUM_RENAME_ATTRS + 1]; unsigned num_attrs=0;
    TranslatorId luft_trans;

    vcf_header_get_attribute (STRa(line), FI_LEN, cSTR("Number="),       false, false, pSTRa(number));
    vcf_header_get_attribute (STRa(line), FI_LEN, cSTR("Type="),         false, false, pSTRa(type));
    vcf_header_get_attribute (STRa(line), FI_LEN, cSTR("ID="),           false, false, pSTRa(id));
    vcf_header_get_attribute (STRa(line), FI_LEN, cSTR(HK_RENDALG_ATTR), true,  false, pSTRa(rendalg));

    // case: liftover algorithm is specified in header, see that Genozip supports it and set context, or error if not 
    if (rendalg) {
        SAFE_NUL (&rendalg[rendalg_len]);
        for (luft_trans=0; luft_trans < NUM_VCF_TRANS; luft_trans++)
            if (!strcmp (rendalg, ltrans_props[luft_trans].alg_name)) 
                goto found;

        ABORT ("Unrecongized Liftover value \"%s\". Offending line:\n%.*s", rendalg, line_len, line);
    
        found: 
        SAFE_RESTORE;
        vcf_header_set_luft_trans (id, id_len, is_info, luft_trans, 0, 0);        
    }

    // case: liftover algorithm not specified - get Genozip's default for this ID
    else {
        luft_trans = vcf_header_set_luft_trans (id, id_len, is_info, TRANS_ID_UNKNOWN, number, number_len);

        attrs[num_attrs++] = (Attr){ .key = HK_RENDALG_ATTR, .key_len = sizeof HK_RENDALG_ATTR-1, 
                                     .value = ltrans_props[luft_trans].alg_name };
    }

    // update renaming attributes
    char new_line[line_len+1]; // +1 for \0
    for (RenameAttr attr_i=0; attr_i < NUM_RENAME_ATTRS; attr_i++) {
        STR (value);
        vcf_header_get_attribute (STRa(line), FI_LEN, STRi(vcf_header_rename_attr,attr_i), true, false, pSTRa(value));
        unsigned header_value_len = value_len;

        // over-writes attrs if there is a pre-existing value, or adds attribute to apriori tag if not
        bool changed = vcf_tags_add_attr_from_header (is_info ? DTYPE_VCF_INFO : DTYPE_VCF_FORMAT, STRa(id), attr_i, 
                                                      STRa(number), STRa(type), ltrans_props[luft_trans].alg_name, strlen (ltrans_props[luft_trans].alg_name),
                                                      pSTRa(value), false);

        if (changed) {
            attrs[num_attrs++] = (Attr){ .key       = vcf_header_rename_attrs[attr_i], 
                                         .key_len   = vcf_header_rename_attr_lens[attr_i], 
                                         .value     = value,
                                         .value_len = value_len };

            if (line != new_line) {
                memcpy (new_line, line, line_len);
                line = new_line;
            }

            if (header_value_len)
                vcf_header_remove_attribute (new_line, &line_len, STRa(attrs[num_attrs-1].key), header_value_len);
        }
    }

    if (num_attrs) 
        vcf_header_add_line_add_attrs (new_header, STRa(line), attrs, num_attrs);
    else
        buf_add_more (NULL, new_header, line, line_len, NULL); // unchanged line
}

// --chain: convert one line at a time to dual-coordinate format, and add all additional lines before #CHROM
static bool vcf_header_zip_update_to_dual_coords_one_line (STRp(line), void *new_header_p, void *unused2, unsigned unused3)
{
    BufferP new_txt_header = (BufferP )new_header_p;
    bool is_info;
    
    if ((is_info = LINEIS ("##INFO=")) || LINEIS ("##FORMAT=")) 
        vcf_header_update_INFO_FORMAT (STRa(line), is_info, new_txt_header);
    
    else if (LINEIS ("##fileformat")) 
        { } // remove this line - it was added to the rejects component instead
        
    else if (chain_is_loaded && LINEIS ("##reference="))
        SUBST_LABEL ("##reference=", HK_ORIGINAL_REF);

    else {
        if (chain_is_loaded && LINEIS ("#CHROM"))
            vcf_header_zip_convert_header_to_dc_add_lines (new_txt_header);

        buf_add_more (evb, new_txt_header, line, line_len, NULL); // just add the line, unmodified
    }

    return false; // continue iterating
}

// ZIP with --chain: add tags that are missing bc they are only in the command line, or symmetrical tags of header tags
static void vcf_header_zip_add_missing_tags (BufferP txt_header)
{
    Tag *tag = NULL;
    txt_header->len -= vcf_field_name_line.len; // remove field name line

    while ((tag = vcf_tags_get_next_missing_tag (tag))) {

        // update number, type, rendalg if missing
        if (!tag->number_len)  { tag->number[0] = '.';             tag->number_len  = 1; }
        if (!tag->type_len)    { memcpy (tag->type, "String", 6);  tag->type_len    = 6; }
        if (!tag->rendalg_len) { memcpy (tag->rendalg, "NONE", 4); tag->rendalg_len = 4; }
        
        bufprintf (evb, txt_header, "##%s=<ID=%.*s,Number=%.*s,Type=%.*s,Description=\"Generated by Genozip due to DVCF tag renaming\","TAG_SOURCE, 
                   DTPT(dtype_names)[tag->dtype], STRf(tag->tag_name), STRf(tag->number), STRf (tag->type));

        if (tag->rendalg_len)
            bufprintf (evb, txt_header, ","HK_RENDALG_ATTR"\"%.*s\"", STRf(tag->rendalg));

        for (RenameAttr attr_i=0; attr_i < NUM_RENAME_ATTRS; attr_i++)
            if (tag->dest_lens[attr_i])
                bufprintf (evb, txt_header, ",%s\"%.*s\"", vcf_header_rename_attrs[attr_i], tag->dest_lens[attr_i], tag->dests[attr_i]);

        bufprint0 (evb, txt_header, ">\n");
    }

    // add back the field name (#CHROM) line
    buf_add_more (evb, txt_header, STRb(vcf_field_name_line), "txt_data");
}

// ZIP with --chain: add liftover_contig, liftover_reference, chain, dual_coordinates keys
static void vcf_header_zip_update_to_dual_coords (BufferP txt_header)
{
    #define new_txt_header evb->codec_bufs[0]
    
    buf_alloc (evb, &new_txt_header, 0, txt_header->len + 262144, char, 0, "evb->codec_bufs[0]"); // initial allocation

    txtfile_foreach_line (txt_header, false, vcf_header_zip_update_to_dual_coords_one_line, &new_txt_header, 0, 0, 0);

    vcf_header_zip_add_missing_tags (&new_txt_header);

    buf_free (*txt_header);
    buf_copy (evb, txt_header, &new_txt_header, char, 0, 0, "txt_header");
    buf_free (new_txt_header);

    vcf_tags_finalize_tags_from_vcf_header();

    #undef new_txt_header
}
 
static SORTER (dst_contigs_sorter) { return ASCENDING_RAW (*(WordIndex *)a, *(WordIndex *)b); }

// scan header for contigs - used to pre-populate z_file contexts in ctx_populate_zf_ctx_from_contigs
// in --chain - also builds vcf_header_liftover_dst_contigs 
static bool vcf_header_handle_contigs (STRp(line), void *new_txt_header_, void *num_contig_lines, unsigned unused2)
{
    bool printed = false;
    BufferP new_txt_header = (BufferP )new_txt_header_;

    const bool is_contig    = LINEIS (HK_CONTIG);
    const bool is_lo_contig = LINEIS (HK_LUFT_CONTIG);
    const bool is_lb_contig = LINEIS (HK_PRIM_CONTIG);

    if (is_contig) 
        (*(uint32_t *)num_contig_lines)++;

    if (!is_contig && !is_lo_contig && !is_lb_contig) goto copy_line; // not a contig line

    ASSERT0 (!chain_is_loaded || (!is_lo_contig && !is_lb_contig), "Cannot use --chain on this file because it already contains \""HK_LUFT_CONTIG"\" or \"" HK_PRIM_CONTIG"\" lines in its header");

    PosType LN = 0;
    unsigned key_len = is_contig    ? STRLEN (HK_CONTIG)
                     : is_lb_contig ? STRLEN (HK_PRIM_CONTIG)
                     :                STRLEN (HK_LUFT_CONTIG);

    // parse eg "##contig=<ID=NKLS02001838.1,length=29167>" and "##luft_contig=<ID=22>"
    SAFE_NUL (&line[line_len]);
    STR (contig_name);
    vcf_header_get_attribute (STRa(line), key_len, cSTR("ID="),     false, true,  pSTRa(contig_name));
    STR(length_str);
    vcf_header_get_attribute (STRa(line), key_len, cSTR("length="), false, false, pSTRa(length_str));
    str_get_int_range64 (STRa(length_str), 1, 1000000000000ULL, &LN); // length stays 0 if length_str_len=0
    SAFE_RESTORE;

    // case: Luft contigs. These only in DVCF files, and when compressing them there is no Luft reference loaded
    if ((is_contig && txt_file->coords == DC_LUFT) || is_lo_contig)
        vcf_header_consume_contig (STRa(contig_name), &LN, true); 

    // case: Primary contigs
    else if (is_contig || is_lb_contig) {

        // case match_chrom_to_reference: update contig_name to the matching name in the reference
        // Note: we match reference contigs by name only, not searching LN, we just verify the LN after the fact.
        // This is so the process is the same when matching CHROM fields while segging, which will result in the header and lines
        // behaving the same
        if (flag.match_chrom_to_reference) {
            Reference ref = chain_is_loaded ? prim_ref : gref;
            WordIndex ref_index = ref_contigs_get_by_name (ref , STRa(contig_name), true, true);
            if (ref_index != WORD_INDEX_NONE) {
                contig_name = ref_contigs_get_name (ref, ref_index, &contig_name_len); // this contig as its called in the reference
                LN = ref_contigs_get_contig_length (ref, ref_index, 0, 0, true);   // update - we check later that it is consistent
            }

            bufprintf (evb, new_txt_header, "%.*s<ID=%.*s,length=%"PRIu64">\n", key_len, line, STRf(contig_name), LN);
            printed = true;
        }

        vcf_header_consume_contig (STRa (contig_name), &LN, false); // also verifies / updates length
    }

    // case: --chain: collect all vcf_header_liftover_dst_contigs. non sorted/unique buf for now, will fix in vcf_header_zip_update_to_dual_coords
    if (chain_is_loaded && !evb->comp_i && is_contig)
        chain_append_all_luft_ref_index (STRa (contig_name), LN, &vcf_header_liftover_dst_contigs);

copy_line:
    if (!printed) 
        buf_add_more (NULL, new_txt_header, line, line_len, NULL); // unchanged
    
    return false; // continue iterating
}

static void vcf_header_add_FORMAT_lines (BufferP txt_header)
{
    txt_header->len -= vcf_field_name_line.len; // remove field name line

    if (flag.GP_to_PP) buf_add_moreC (evb, txt_header, "##FORMAT=<ID=PP,Number=G,Type=Integer,Description=\"Phred-scaled genotype posterior probabilities rounded to the closest integer\">\n", "txt_data");
    if (flag.GL_to_PL) buf_add_moreC (evb, txt_header, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods rounded to the closest integer\">\n", "txt_data");

    // add back the field name (#CHROM) line
    buf_add_more (evb, txt_header, STRb(vcf_field_name_line), "txt_data");
}

static void vcf_add_all_ref_contigs_to_header (BufferP txt_header)
{
    txt_header->len -= vcf_field_name_line.len; // remove field name line

    ConstContigPkgP prim_contigs = ref_get_ctgs (prim_ref);
    ARRAY (const Contig, ctgs, prim_contigs->contigs);

    for (uint64_t i=0; i < ctgs_len; i++) {
        rom contig_name = Bc (prim_contigs->dict, ctgs[i].char_index);
        bufprintf (evb, txt_header, VCF_CONTIG_FMT"\n", ctgs[i].snip_len, contig_name, ctgs[i].max_pos);

        chain_append_all_luft_ref_index (contig_name, ctgs[i].snip_len, ctgs[i].max_pos, &vcf_header_liftover_dst_contigs);
    }

    // add back the field name (#CHROM) line
    buf_add_more (evb, txt_header, STRb(vcf_field_name_line), "txt_data");
}

static bool vcf_header_build_stats_programs (STRp(line), void *unused1, void *unused2, unsigned unused3)
{
    if (LINEIS ("##source=")) {
        buf_append (evb, stats_programs, char, &line[9], line_len-9, "stats_programs");
        *BLSTc (stats_programs) = ';'; // replace \n with ; - the separator required by stats_programs
    }

    return false; // continue iterating
}

static bool vcf_inspect_txt_header_zip (BufferP txt_header)
{
    if (!vcf_header_set_globals (txt_file->name, txt_header, true)) return false; // samples are different than a previous concatented file

    txtfile_foreach_line (txt_header, false, vcf_header_build_stats_programs, 0, 0, 0, 0);
    if (stats_programs.len) *BAFTc (stats_programs) = 0; // nul-terminate

    SAFE_NUL (BAFTc (*txt_header));
    #define IF_IN_SOURCE(signature, segcf) if (stats_programs.len && strstr (stats_programs.data, signature)) segconf.segcf = true
    #define IF_IN_HEADER(signature, segcf) if (strstr (txt_header->data   , signature)) segconf.segcf = true
    
    IF_IN_SOURCE ("infiniumFinalReportConverter", vcf_infinium);
    IF_IN_SOURCE ("VarScan", vcf_is_varscan);
    IF_IN_SOURCE ("dbSNP", vcf_dbSNP);
    IF_IN_HEADER ("GenotypeGVCFs", vcf_is_gvcf);
    IF_IN_HEADER ("CombineGVCFs", vcf_is_gvcf);
    IF_IN_HEADER ("beagle", vcf_is_beagle);
    IF_IN_HEADER ("Illumina GenCall", vcf_illum_gtyping);
    IF_IN_HEADER ("Log R Ratio", vcf_illum_gtyping);
    IF_IN_HEADER ("Number of cases used to estimate genetic effect", vcf_is_gwas); // v1.0
    IF_IN_HEADER ("##trait", vcf_is_gwas); /*v1.2*/
    
    bool has_PROBE = !!strstr (txt_header->data, "##INFO=<ID=PROBE_A");
    SAFE_RESTORE;

    if (!flag.reference && segconf.vcf_is_gvcf)
        TIP ("compressing a GVCF file using a reference file can reduce the compressed file's size by 10%%-30%%.\nUse: \"genozip --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n",
            txt_file->name);

    if (!flag.reference && segconf.vcf_illum_gtyping && has_PROBE)
        TIP ("compressing an Illumina Genotyping VCF file using a reference file can reduce the compressed file's size by 20%%.\nUse: \"genozip --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n",
             txt_file->name);

    if (chain_is_loaded && !evb->comp_i)
        vcf_tags_populate_tags_from_command_line();
        
    // scan header for ##dual_coordinates - this sets txt_file->coords
    txtfile_foreach_line (txt_header, false, vcf_header_get_dual_coords, 0, 0, 0, 0);
    
    // handle contig lines
    ASSERTNOTINUSE (vcf_header_liftover_dst_contigs);
    if (chain_is_loaded && !evb->comp_i)
        buf_alloc (evb, &vcf_header_liftover_dst_contigs, 0, 100, WordIndex, 0, "codec_bufs[1]");
    
    uint32_t num_contig_lines=0;
    vcf_header_rewrite_header (evb, txt_header, vcf_header_handle_contigs, &num_contig_lines);

    if (chain_is_loaded && !evb->comp_i)     // sort (but not uniq)
        qsort (STRb(vcf_header_liftover_dst_contigs), sizeof (WordIndex), dst_contigs_sorter);

    // in liftover, if header has no contigs, add reference contigs and derived ocontigs
    if (chain_is_loaded && !evb->comp_i && !num_contig_lines) 
        vcf_add_all_ref_contigs_to_header (txt_header);

    // add PP and/or PL INFO lines if needed
    if ((flag.GP_to_PP || flag.GL_to_PL) && !evb->comp_i)
        vcf_header_add_FORMAT_lines (txt_header);

    // set vcf_version
    unsigned fileformat_line_len = vcf_header_parse_fileformat_line (txt_header);

    // if --chain: scan header for ##FORMAT and ##INFO - for translatable subfields - create contexts, and populate Context.luft_trans

    // if we're compressing a txt file with liftover, we prepare the two rejects headers: it consists of two lines:
    // 1. the first "##fileformat" line if one exists in our file 
    // 2. the last "#CHROM" line - needing for zipping the rejects file, but will be removed in reconstruction
    // in reconstruction without --luft, we will reconstruct the normal header 
    // in reconstruction with --luft, we will reconstruct the rejected header first, then the normal header without the first line
    if ((chain_is_loaded || txt_file->coords) && !evb->comp_i) {
                
        STRw(chrom_line) = vcf_header_get_last_line (txt_header, &chrom_line); // including sample names

        // write both PRIM_ONLY and LUFT_ONLY txt headers (they are the same)
        ASSERTNOTINUSE (evb->scratch);
        buf_add_more (evb, &evb->scratch, B1STc(*txt_header), fileformat_line_len, "scratch"); // add the ##fileformat line, if there is one
        buf_add_more (evb, &evb->scratch, STRa(chrom_line), "scratch"); // add the #CHROM line
 
        txtheader_compress (&evb->scratch, evb->scratch.len, DIGEST_NONE, false, VCF_COMP_PRIM_ONLY); 
        txtheader_compress (&evb->scratch, evb->scratch.len, DIGEST_NONE, false, VCF_COMP_LUFT_ONLY);
        
        buf_free (evb->scratch);
    }

    // case: Luft file - update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
    if (txt_file->coords == DC_LUFT)
        vcf_header_zip_liftback_header (txt_header);
    else if (txt_file->coords == DC_PRIMARY) 
        vcf_header_zip_get_luft_only_lines (txt_header);

    // case: --chain: add liftover_contig, liftover_reference, chain, dual_coordinates keys
    if ((chain_is_loaded || txt_file->coords) && evb->comp_i == VCF_COMP_MAIN)
        vcf_header_zip_update_to_dual_coords (txt_header);

    //xxx if (z_file && !z_file->num_txt_files)  // first component 
    if (z_file)
        z_file->z_flags.has_gencomp |= (txt_file->coords != DC_NONE);  // note: in case of --chain, z_flags.dual_coords is set in file_open_z

    return true; // all good
}

static bool vcf_inspect_txt_header_piz (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags)
{
    vcf_header_set_globals (z_file->name, txt_header, false);

    if (flag.genocat_no_reconstruct) return true;

    // remove #CHROM line (it is saved in vcf_field_name_line by vcf_header_set_globals()) - or everything
    // if we ultimately want the #CHROM line. We will add back the #CHROM line later.
    if (flag.header_one) 
        txt_header->len = 0; 
    else
        txt_header->len -= vcf_field_name_line.len; 

    // if genocat --samples is used, update vcf_header and vcf_num_displayed_samples
    // note: we subset reject samples (in the header) as well - we analyze only in the rejects section
    if (flag.samples) vcf_header_subset_samples (&vcf_field_name_line);
    else              vcf_num_displayed_samples = vcf_num_samples;

    // for the rejects part of the header - we're done
    if (evb->comp_i) 
        return true;

    // if needed, liftover the header
    if (flag.luft && !flag.header_one && !evb->comp_i) 
        vcf_header_rewrite_header (txt_header_vb, txt_header, vcf_header_piz_liftover_header, NULL);

    if (flag.single_coord)
        vcf_header_rewrite_header (txt_header_vb, txt_header, vcf_header_piz_single_coord, NULL);

    // if needed, add ##INFO oSTATUS line
    if (flag.show_ostatus)
        bufprintf (txt_header_vb, txt_header, KH_INFO_oSTATUS"\n", GENOZIP_CODE_VERSION);

    // add genozip command line
    if (!flag.header_one && is_genocat && !flag.genocat_no_reconstruct && !evb->comp_i
        && !flag.no_pg && (flag.data_modified || z_is_dvcf)) 
        vcf_header_add_genozip_command (txt_header_vb, txt_header);

    if (flag.drop_genotypes) 
        vcf_header_trim_field_name_line (&vcf_field_name_line); // drop FORMAT and sample names

    // if --show_dvcf, add two extra fields at beginning for field name line, and remove the # from #CHROM
    if (flag.show_dvcf)
        buf_add_moreC (txt_header_vb, txt_header,  "#COORD\toSTATUS\t", "txt_data");

    // add the (perhaps modified) field name (#CHROM) line
    buf_add_more (txt_header_vb, txt_header, &vcf_field_name_line.data[!!flag.show_dvcf], vcf_field_name_line.len - !!flag.show_dvcf, "txt_data");

    return true; // all good
}

bool vcf_inspect_txt_header (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags)
{
    return (IS_ZIP) ? vcf_inspect_txt_header_zip (txt_header)
                            : vcf_inspect_txt_header_piz (txt_header_vb, txt_header, txt_header_flags);
}

static bool vcf_header_set_globals (rom filename, BufferP vcf_header, FailType soft_fail)
{
    static rom vcf_field_name_line_filename = NULL; // file from which the header line was taken

    // check for ##fileformat=VCF
    vcf_has_file_format_vcf = vcf_header->len >= 16 && !memcmp (vcf_header->data, "##fileformat=VCF", 16);

    // count tabs in last line which should be the field header line
    unsigned tab_count = 0;
    for (int i=vcf_header->len-1; i >= 0; i--) {
        
        if (vcf_header->data[i] == '\t')
            tab_count++;
        
        // some times files have spaces instead of \t 
        else if (vcf_header->data[i] == ' ') {
            tab_count++;
            while (i >= 1 && vcf_header->data[i-1]==' ') i--; // skip run of spaces
        }

        // if this is the beginning of field header line 
        else if (vcf_header->data[i] == '#' && (i==0 || vcf_header->data[i-1] == '\n' || vcf_header->data[i-1] == '\r')) {
        
            // ZIP: if first vcf file ; PIZ: everytime - copy the header to the global
            if (!buf_is_alloc (&vcf_field_name_line) || IS_PIZ) {
                buf_copy (evb, &vcf_field_name_line, vcf_header, char, i, vcf_header->len - i, "vcf_field_name_line");
                vcf_field_name_line_filename = filename;
            }

            // ZIP only: subsequent files - if we're in bound mode just compare to make sure the header is the same
            else if (flag.bind && 
                     (vcf_header->len-i != vcf_field_name_line.len || memcmp (vcf_field_name_line.data, &vcf_header->data[i], vcf_field_name_line.len))) {

                WARN ("%s: skipping %s: it has a different VCF header line than %s, see below:\n"
                      "========= %s =========\n"
                      "%.*s"
                      "========= %s ==========\n"
                      "%.*s"
                      "=======================================\n", 
                      global_cmd, filename, vcf_field_name_line_filename,
                      vcf_field_name_line_filename, STRfb(vcf_field_name_line),
                      filename, vcf_header->len32-i, &vcf_header->data[i]);
                
                if (soft_fail) return false;
                else           exit_on_error (false);
            }

            // count samples
            vcf_num_samples = (tab_count >= 9) ? tab_count-8 : 0; 
            
            // note: a VCF file without samples may or may not have a "FORMAT" in the header, i.e. tab_count==7 or 8 (8 or 9 fields).
            // however, even if it has a FORMAT in the header, it won't have a FORMAT column in the data

            ASSINP (tab_count >= 7, "Error: invalid VCF field/samples header line - it contains only %d fields, expecting at least 8", tab_count+1);

            return true; 
        }
    }

    ABORT_R ("Error: invalid VCF file - it does not contain a field header line; tab_count=%u", tab_count+1);
}

// genocat: remove FORMAT and sample names from the vcf header line, in case of --drop-genotypes
// TODO - remove FORMAT lines from header itself
static void vcf_header_trim_field_name_line (BufferP vcf_header)
{
    char *line;
    unsigned line_len = vcf_header_get_last_line (vcf_header, &line);

    if (is_field_name_line (STRa(line))) {
        vcf_header->len = strstr (line, "INFO") + 5 - B1STc (*vcf_header);
        *BLSTc (*vcf_header) = '\n';
    }                  
}

uint32_t vcf_header_get_num_samples (void)
{
    if (Z_DT(VCF) || Z_DT(BCF))
        return vcf_num_samples;
    else
        return 0;
}

// -------------
// samples stuff
// -------------

// processes the vcf header sample line according to the --samples option, removing samples that are not required
// and building a bytemap PIZ filter the samples
static void vcf_header_subset_samples (BufferP vcf_field_name_line)
{
    // accept a sample from the vcf file's samples as consistent with the --samples requested
    #define samples_accept(sample_str) { \
        buf_append_string (evb, vcf_field_name_line, sample_str); \
        buf_append_string (evb, vcf_field_name_line, "\t"); \
        vcf_num_displayed_samples++; \
    }

    ARRAY (char, line, *vcf_field_name_line);

    flag.samples = is_field_name_line (STRa(line));
    RETURNW (flag.samples,, "Warning: found non-standard VCF sample header line. Ingoring --samples : \n%.*s", (int)line_len, line);

    int32_t num_samples=-8;
    for (unsigned i=0; i < line_len; i++)
        if (line[i] == '\t') num_samples++;
        
    ASSINP (num_samples >= 1, "Cannot use --samples with %s, as this file has no samples", z_name);
    
    vcf_samples_is_included = MALLOC (num_samples);
    memset (vcf_samples_is_included, cmd_is_negative_samples, num_samples); // 0 if not included unless list says so (positive) and vice versa

    unsigned vcf_names_start_index = (line + sizeof VCF_FIELD_NAMES_LONG-1 + 1) - B1STc (*vcf_field_name_line);
    unsigned vcf_names_data_len = vcf_field_name_line->len - vcf_names_start_index;
    vcf_sample_names_data = MALLOC (vcf_names_data_len);
    memcpy (vcf_sample_names_data, &vcf_field_name_line->data[vcf_names_start_index], vcf_names_data_len);
    vcf_sample_names_data[vcf_names_data_len-1] = '\t'; // change last separator from \n to \t

    vcf_sample_names = MALLOC (num_samples * sizeof (char *));

    vcf_field_name_line->len = vcf_names_start_index;
    char *next_token = vcf_sample_names_data;
    
    // go through the vcf file's samples and add those that are consistent with the --samples requested
    vcf_num_displayed_samples = 0;

    buf_copy (evb, &evb->codec_bufs[0], &cmd_samples_buf, char*, 0, 0, 0);
    ARRAY (rom , snames, evb->codec_bufs[0]);

    int fixed_len = (snames_len == 1 && atoi (snames[0]) >= 1) ? atoi (snames[0]) : 0;

    for (unsigned i=0; i < num_samples; i++) {
        vcf_sample_names[i] = strtok_r (next_token, "\t", &next_token);
        bool handled = false;

        // case: --samples <number>
        if (fixed_len && i < fixed_len) {
            vcf_samples_is_included[i] = true;
            samples_accept (vcf_sample_names[i]);
            (*(uint64_t *)&snames_len) = 0;
            continue;
        }
        
        for (unsigned s=0; s < snames_len; s++) 
            if (!strcmp (vcf_sample_names[i], snames[s])) { // found
                vcf_samples_is_included[i] = !cmd_is_negative_samples;
                if (!cmd_is_negative_samples) 
                    samples_accept (vcf_sample_names[i]);

                // remove this sample from the --samples list as we've found it already
                memcpy (&snames[s], &snames[s+1], (snames_len-s-1) * sizeof (char *));
                (*(uint64_t *)&snames_len)--; // override const
                handled = true;
                break;
            }
            
        if (!handled && cmd_is_negative_samples) 
            samples_accept (vcf_sample_names[i]);
    }

    *BLSTc (*vcf_field_name_line) = '\n'; // change last separator from \t

    // warn about any --samples items that were not found in the vcf file (all the ones that still remain in the buffer)
    for (unsigned s=0; s < snames_len; s++) 
        ASSERTW (false, "Warning: requested sample '%s' is not found in the VCF file, ignoring it", snames[s]);

    // if the user filtered out all samples, its equivalent of drop_genotypes
    if (!vcf_num_displayed_samples) flag.drop_genotypes = true;

    buf_free (evb->codec_bufs[0]);
}

// called from genozip.c for processing the --samples flag
void vcf_samples_add  (rom samples_str)
{
    ASSERTNOTNULL (samples_str);

    // --samples 0 is interpreted as --drop-genotypes
    if (samples_str[0]=='0' && !samples_str[1] && !cmd_samples_buf.len) {
        flag.drop_genotypes = true;
        flag.samples = false;
    }

    bool is_negated = samples_str[0] == '^';

    bool is_conflicting_negation = (cmd_samples_buf.len && (cmd_is_negative_samples != is_negated));
    ASSINP0 (!is_conflicting_negation, "Error: inconsistent negation - all samples listed must either be negated or not");

    cmd_is_negative_samples = is_negated;

    // make a copy of the string and leave the original one for error message. 
    // we don't free this memory as chrom fields in regions will be pointing to it
    char *next_region_token = MALLOC (strlen (samples_str)+1); // heap memory, as cmd_samples_buf elements point into this
    strcpy (next_region_token, samples_str + is_negated); // drop the ^ if there is one

    while (1) {
        char *one_sample = strtok_r (next_region_token, ",", &next_region_token);
        if (!one_sample) break;

        bool is_duplicate = false;
        for (unsigned s=0; s < cmd_samples_buf.len32; s++)
            if (!strcmp (one_sample, *B(char *, cmd_samples_buf, s))) {
                is_duplicate = true;
                break;
            }
        if (is_duplicate) continue; // skip duplicates "genocat -s sample1,sample2,sample1"

        buf_alloc (evb, &cmd_samples_buf, 1, 100, char*, 2, "cmd_samples_buf");

        BNXT (char *, cmd_samples_buf) = one_sample;
    }
}

VcfVersion vcf_header_get_version (void) { return vcf_version; }