// ------------------------------------------------------------------
//   header_vcf.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

// Globals
uint32_t vcf_num_samples = 0; // number of samples in the file
static uint32_t vcf_num_displayed_samples = 0; // PIZ only: number of samples to be displayed - might be less that vcf_num_samples if --samples is used
static Buffer vcf_field_name_line = EMPTY_BUFFER;  // header line of first VCF file read - use to compare to subsequent files to make sure they have the same header during bound

#define vcf_header_liftover_dst_contigs evb->codec_bufs[1]

// referring to sample strings from the --samples command line option
static Buffer cmd_samples_buf = EMPTY_BUFFER; // an array of (char *)
static bool cmd_is_negative_samples = false;

// referring to samples in the vcf file
char *vcf_samples_is_included;         // a bytemap indicating for each sample if it is included
static char **vcf_sample_names;        // an array of char * to nul-terminated names of samples 
static char *vcf_sample_names_data;    // vcf_sample_names point into here

#define LINEIS(s) (line_len > (sizeof s - 1) && !memcmp (line, (s), sizeof s - 1)) 

#define SUBST_LABEL(old,new) do { buf_add_more (NULL, new_txt_header, (new), (sizeof new - 1), NULL); \
                                  buf_add_more (NULL, new_txt_header, line + (sizeof old - 1), line_len - (sizeof old - 1), NULL); } while (0)

// use in --chain for generating LIFTREJT strings with full INFO and FORMAT ID names (but just the 8 chars in ctx->name)
typedef struct { DictId dict_id; char ID[MAX_VCF_ID_LEN]; } DictIdToID;
static Buffer dict_id_to_ID = EMPTY_BUFFER; // array of DictIdToID

void vcf_header_piz_init (void)
{
    vcf_num_samples = 0;
    buf_free (&vcf_field_name_line);

    // note: we don't re-initialize vcf_num_displayed_samples - this is calculated only once
}


static void vcf_header_subset_samples (Buffer *vcf_header);
static bool vcf_header_set_globals (const char *filename, Buffer *vcf_header, bool soft_fail);
static void vcf_header_trim_field_name_line (Buffer *vcf_header);

// returns the length of the first line, if it starts with ##fileformat, or 0
static inline unsigned vcf_header_get_first_line_len (const Buffer *txt_header)
{
    ARRAY (char, line, *txt_header);
    if (!LINEIS ("##fileformat")) return 0;

    const char *newline = memchr (line, '\n', line_len);
    return newline ? newline - line + 1 : 0;
}

static inline bool is_field_name_line (const char *line, unsigned line_len)
{
    // compare to "#CHROM" - not to entire VCF_FIELD_NAMES as it may be a line maybe separated by spaces instead of tabs
    return LINEIS ("#CHROM");
}

// returns the length of the last line, if it starts with #CHROM, or 0
static bool vcf_header_get_last_line_cb (const char *line, unsigned line_len, void *unused1, void *unused2, unsigned unused3)
{ return true; }
static unsigned vcf_header_get_last_line (Buffer *txt_header, char **line_p)
{
    int64_t line_len;
    char *line = txtfile_foreach_line (txt_header, true, vcf_header_get_last_line_cb, 0, 0, 0, &line_len);
    if (line_p) *line_p = line;

    return is_field_name_line (line, line_len) ? (unsigned)line_len : 0;
}

static void vcf_header_add_genozip_command (VBlockP txt_header_vb, Buffer *txt_header)
{
    // the command line length is unbound, careful not to put it in a bufprintf
    buf_add_string (txt_header_vb, txt_header, HK_GENOZIP_CMD"\"");
    buf_add_string (txt_header_vb, txt_header, flags_command_line());
    bufprintf (txt_header_vb, txt_header, "\" %s\n", str_time().s);
}

static void vcf_header_get_attribute (const char *line, unsigned line_len, unsigned key_len,
                                      const char *attr, unsigned attr_len, bool enforce,// in
                                      const char **snip, unsigned *snip_len) // out
{
    SAFE_NUL (&line[line_len]);
    const char *start = strstr (line + key_len, attr);
    SAFE_RESTORE;

    if (!start) goto missing; // attr doesn't exist or malformed line

    start += attr_len;

    const char *comma   = memchr (start, ',', &line[line_len] - start);
    const char *bracket = memchr (start, '>', &line[line_len] - start);

    if (!comma && !bracket) goto missing; // attr isn't terminated with a comma or >
    const char *after = !comma ? bracket : !bracket ? comma : MIN (comma,bracket);
    
    *snip_len = after - start;
    *snip = *snip_len ? start : NULL;
    return;

missing:
    *snip = NULL;
    *snip_len = 0;
    ASSINP (!enforce, "missing attr %s in header line %.*s", attr, line_len, line);
}

// add contig to buffer if line is a valid contig, and fail silently if not
static void vcf_header_consume_contig (const char *contig_name, unsigned contig_name_len, PosType length, bool liftover)
{
    // if its the 2nd+ file while binding in ZIP - we just check that the contigs are the same
    // (or a subset) of previous file's contigs (to do: merge headers, bug 327) 
    if (command == ZIP && flag.bind && z_file->num_txt_components_so_far /* 2nd+ file */) 
        txtheader_verify_contig (contig_name, contig_name_len, length, (void *)liftover);

    else {
        txtheader_alloc_contigs (1, contig_name_len+1, liftover);
        txtheader_add_contig (contig_name, contig_name_len, length, (void *)liftover);
    }
}

static bool vcf_header_extract_contigs (const char *line, unsigned line_len, void *dst_contigs_buf, void *unused2, unsigned unused3)
{
    bool is_contig    = LINEIS ("##contig=");
    bool is_lo_contig = LINEIS (HK_LUFT_CONTIG);
    bool is_lb_contig = LINEIS (HK_PRIM_CONTIG);

    if (!is_contig && !is_lo_contig && !is_lb_contig) goto done;

    ASSERT0 (!chain_is_loaded || (!is_lo_contig && !is_lb_contig), "Cannot use --chain on this file because it already contains \""HK_LUFT_CONTIG"\" or \"" HK_PRIM_CONTIG"\" lines in its header");

    const char *contig_name;
    unsigned contig_name_len;
    PosType length = 0;
    unsigned key_len = is_contig ? 9 : (sizeof HK_LUFT_CONTIG - 1);
    DECLARE_SNIP;

    // parse eg "##contig=<ID=NKLS02001838.1,length=29167>" and "##luft_contig=<ID=22>"
    SAFE_NUL (&line[line_len]);
    vcf_header_get_attribute (line, line_len, key_len, "ID=", 3, true, &contig_name, &contig_name_len);
    vcf_header_get_attribute (line, line_len, key_len, "length=", 7, false, &snip, &snip_len);
    str_get_int_range64 (snip, snip_len, 1, 1000000000000ULL, &length); // length stays 0 if snip_len=0
    SAFE_RESTORE;

    if ((is_contig && txt_file->coords == DC_LUFT) || is_lo_contig)
        vcf_header_consume_contig (contig_name, contig_name_len, length, true);

    else if (is_contig || is_lb_contig) 
        vcf_header_consume_contig (contig_name, contig_name_len, length, false);

    // case: --chain: collect all vcf_header_liftover_dst_contigs. non sorted/unique buf for now, will fix in vcf_header_zip_update_to_dual_coords
    if (chain_is_loaded && !flag.rejects_coord && LINEIS ("##contig="))
        chain_append_all_dst_contig_index (contig_name, contig_name_len, (Buffer *)dst_contigs_buf);

done:
    return false; // continue iterating
}


static bool vcf_header_get_dual_coords (const char *line, unsigned line_len, void *unused1, void *unused2, unsigned unused3)
{
    if (LINEIS (HK_DC_PRIMARY)) {
        ASSERT (!chain_is_loaded, "--chain cannot be used with %s - it is already a dual-coordinates VCF file - it contains \""HK_DC_PRIMARY"\" in its header", txt_name);
        txt_file->coords = DC_PRIMARY;
        flag.bind = BIND_REJECTS;
        z_file->digest_ctx_bound = z_file->digest_ctx_single; // we already calculated the digest for the header
    }

    else if (LINEIS (HK_DC_LUFT)) {
        ASSERT (!chain_is_loaded, "--chain cannot be used with %s - it is already a dual-coordinates VCF file - it contains \""HK_DC_LUFT"\" in its header", txt_name);
        txt_file->coords = DC_LUFT;
        flag.bind = BIND_REJECTS;
        flag.data_modified = true; // we will be zipping this as a dual-coordinates VCF with default reconstruction as PRIMARY - this turns off digest
        z_file->digest_ctx_bound  = (DigestContext){};
        z_file->digest_ctx_single = (DigestContext){};
    }

    else return false; // continue iterating

    // genozip of a dual coordinates file implies --sort, unless overridden with --unsorted
    if (!flag.unsorted) flag.sort = true;

    ASSINP0 (!flag.test, "--test is not supported for dual-coordinates files");

    return true; // found - no more iterating
}

// PIZ with --luft - one line: remove fileformat, update *contig, *reference, dual_coordinates keys to Luft format
static bool vcf_header_piz_liftover_header_one_line (const char *line, unsigned line_len, 
                                                     void *new_txt_header_, void *unused1, unsigned unused2)
{
    Buffer *new_txt_header = (Buffer *)new_txt_header_;

    if      (LINEIS (HK_LUFT_CONTIG))   SUBST_LABEL (HK_LUFT_CONTIG, "##contig=");
    else if (LINEIS ("##contig="))    SUBST_LABEL ("##contig=", HK_PRIM_CONTIG);
    else if (LINEIS (HK_LUFT_REF))      SUBST_LABEL (HK_LUFT_REF, "##reference=");
    else if (LINEIS ("##reference=")) SUBST_LABEL ("##reference=", HK_PRIM_REF);
    else if (LINEIS (HK_DC_PRIMARY))  SUBST_LABEL (HK_DC_PRIMARY, HK_DC_LUFT);
    else                              buf_add_more (NULL, new_txt_header, line, line_len, NULL);

    return false; // continue iterating
}

// PIZ with --luft: remove fileformat, update *contig, *reference, dual_coordinates keys to Luft format
static void vcf_header_piz_liftover_header (VBlockP txt_header_vb, Buffer *txt_header)
{
    #define new_txt_header txt_header_vb->codec_bufs[0]

    buf_alloc (txt_header_vb, &new_txt_header, 0, txt_header->len, char, 0, "codec_bufs[0]"); // initial allocation (might be a bit bigger due to label changes)

    txtfile_foreach_line (txt_header, false, vcf_header_piz_liftover_header_one_line, &new_txt_header, 0, 0, 0);

    // replace txt_header with lifted back one
    buf_free (txt_header);
    buf_copy (txt_header_vb, txt_header, &new_txt_header, char, 0, 0, "txt_data");
    buf_free (&new_txt_header);
    #undef new_txt_header
}

// squeeze in the liftover rejects from the header, with their label removed, before the existing unconsumed text
// (data that was read beyond the header) 
static void vcf_header_move_rejects_to_unconsumed_text (Buffer *rejects)
{
    if (rejects->len) {
        buf_alloc (evb, &txt_file->unconsumed_txt, rejects->len, 0, char, 0, "txt_file->unconsumed_txt");
        ARRAY (char, unconsumed_txt, txt_file->unconsumed_txt);
        
        memmove (&unconsumed_txt[rejects->len], unconsumed_txt, unconsumed_txt_len);
        memcpy (unconsumed_txt, rejects->data, rejects->len);
        
        txt_file->unconsumed_txt.len += rejects->len;
        txt_file->reject_bytes = rejects->len;
    }

    buf_free (rejects);
}

// ZIP of Luft line: update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
static bool vcf_header_zip_liftback_header_one_line (const char *line, unsigned line_len, 
                                                     void *new_txt_header_, void *rejects_, unsigned unused)
{
    Buffer *new_txt_header = (Buffer *)new_txt_header_;
    
    if      (LINEIS (HK_PRIM_ONLY))   buf_add_more (evb, (Buffer *)rejects_, line + (sizeof HK_PRIM_ONLY -1), line_len - (sizeof HK_PRIM_ONLY -1), NULL); // add ##primary_only variants to "rejects" buffer
    else if (LINEIS (HK_PRIM_CONTIG)) SUBST_LABEL (HK_PRIM_CONTIG, "##contig=");
    else if (LINEIS ("##contig="))    SUBST_LABEL ("##contig=", HK_LUFT_CONTIG);
    else if (LINEIS (HK_PRIM_REF))    SUBST_LABEL (HK_PRIM_REF, "##reference=");
    else if (LINEIS ("##reference=")) SUBST_LABEL ("##reference=", HK_LUFT_REF);
    else if (LINEIS (HK_DC_LUFT))     SUBST_LABEL (HK_DC_LUFT, HK_DC_PRIMARY);
    else                              buf_add_more (evb, new_txt_header, line, line_len, NULL);

    return false; // continue iterating
}

// ZIP of LUFT-coords file: update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
static void vcf_header_zip_liftback_header (Buffer *txt_header)
{
    #define new_txt_header evb->codec_bufs[0]
    #define rejects        evb->codec_bufs[1]

    // case: Luft file - scan the header and for LIFTBACK_REJECT lines, and move them to txt_file->unconsumed (after removing the label) 
    // to be processsed as normal variants
    buf_alloc (evb, &new_txt_header, 0, txt_header->len + 16384, char, 0, "evb->codec_bufs[0]"); // initial allocation (might be a bit bigger due to label changes)
    buf_alloc (evb, &rejects, 0, 65536, char, 0, "evb->codec_bufs[1]"); // initial allocation

    txtfile_foreach_line (txt_header, false, vcf_header_zip_liftback_header_one_line, &new_txt_header, &rejects, 0, 0);

    // replace txt_header with lifted back one
    buf_free (txt_header);
    buf_copy (evb, txt_header, &new_txt_header, char, 0, 0, "txt_data");

    // add rejects to the beginning of unconsumed_txt
    vcf_header_move_rejects_to_unconsumed_text (&rejects);

    buf_free (&new_txt_header); 

    #undef new_txt_header
    #undef rejects
}

static bool vcf_header_zip_get_luft_only_lines_do (const char *line, unsigned line_len, void *rejects_, void *unused1, unsigned unused2)
{
    // add ##luft_only variants to "rejects" buffer
    if (LINEIS (HK_LUFT_ONLY)) buf_add_more (evb, (Buffer *)rejects_, line + sizeof HK_LUFT_ONLY -1, line_len - sizeof HK_LUFT_ONLY -1, NULL); 
    return false; // continue iterating
}

// ZIP of PRIMARY-coords file: update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
static void vcf_header_zip_get_luft_only_lines (Buffer *txt_header)
{
    #define rejects evb->codec_bufs[1]

    buf_alloc (evb, &rejects, 0, 65536, char, 0, "evb->codec_bufs[1]"); // initial allocation

    txtfile_foreach_line (txt_header, false, vcf_header_zip_get_luft_only_lines_do, &rejects, 0, 0, 0);

    // add rejects to the beginning of unconsumed_txt
    vcf_header_move_rejects_to_unconsumed_text (&rejects);

    #undef rejects
}

// add a key=value to the end of an existing VCF header line
static void vcf_header_add_line_add_key (Buffer *new_header, const char *line, unsigned line_len, const char *key, const char *value)
{
    bool has_nl, has_13, has_bt;
    line_len -= (has_nl = (line_len >= 1 && line[line_len-1] == '\n'));
    line_len -= (has_13 = (line_len >= 1 && line[line_len-1] == '\r'));
    line_len -= (has_bt = (line_len >= 1 && line[line_len-1] == '>'));
    
    ASSERT (has_nl && has_bt, "While adding key=%s: malformed VCF header line: \"%.*s\"", key, line_len, line);

    buf_add_more (NULL, new_header, line, line_len, NULL);
    NEXTENT (char, *new_header) = ',';
    buf_add_more (NULL, new_header, key, strlen (key), NULL);
    NEXTENT (char, *new_header) = '=';
    buf_add_more (NULL, new_header, value, strlen (value), NULL);
    NEXTENT (char, *new_header) = '>';
    if (has_13) NEXTENT (char, *new_header) = '\r';
    NEXTENT (char, *new_header) = '\n';
}

static void vcf_header_zip_convert_header_to_dc_add_lines (Buffer *new_header)
{
    bufprintf (evb, new_header, "%s\n", HK_DC_PRIMARY);
    bufprintf (evb, new_header, HK_CHAIN"file://%s\n", chain_filename ? chain_filename : "<null>");
    bufprintf (evb, new_header, HK_LUFT_REF"file://%s\n", ref_filename ? ref_filename : "<null>");

    bufprintf (evb, new_header, KH_INFO_LUFT"\n", GENOZIP_CODE_VERSION);
    bufprintf (evb, new_header, KH_INFO_PRIM"\n", GENOZIP_CODE_VERSION);
    bufprintf (evb, new_header, KH_INFO_LREJ"\n", GENOZIP_CODE_VERSION);
    bufprintf (evb, new_header, KH_INFO_PREJ"\n", GENOZIP_CODE_VERSION);

    // add all Luft contigs (note: we get them from liftover directly, as they zip_initialize that adds them to zf_ctx has not been called yet)
    ARRAY (WordIndex, dst_contig, vcf_header_liftover_dst_contigs); // dst_contig is sorted, but might contain duplicate elements
    for (unsigned i=0; i < dst_contig_len; i++) {
        if (i && dst_contig[i] == dst_contig[i-1]) continue; // skip duplicate

        PosType length;
        const char *luft_contig = chain_get_dst_contig (dst_contig[i], &length);

        if (length) bufprintf (evb, new_header, HK_LUFT_CONTIG"<ID=%s,length=%"PRId64">\n", luft_contig, length);
        else        bufprintf (evb, new_header, HK_LUFT_CONTIG"<ID=%s>\n", luft_contig);
    }
    buf_free (&vcf_header_liftover_dst_contigs);
}

// Update zf_ctx and dict_id_to_ID after reading or creating an INFO / FORMAT RenderAlg=vale
static TranslatorId vcf_header_set_luft_trans (const char *id, unsigned id_len, bool is_info, TranslatorId luft_trans, const char *number, unsigned number_len)
{
    DictId dict_id = dict_id_make (id, id_len, is_info ? DTYPE_VCF_INFO : DTYPE_VCF_FORMAT);

    // add to dict_id_to_ID
    buf_alloc (NULL, &dict_id_to_ID, 1, 0, DictIdToID, 2, NULL);
    DictIdToID *ent = &NEXTENT (DictIdToID, dict_id_to_ID);
    *ent = (DictIdToID){.dict_id = dict_id};
    memcpy (LASTENT (DictIdToID, dict_id_to_ID)->ID, id, MIN (id_len, sizeof (ent->ID)-1)); // already nul-terminated as array is initialized to zeros

    if (luft_trans == TRANS_ID_UNKNOWN) 
        luft_trans = (number_len==1) ? vcf_lo_luft_trans_id (dict_id, *number) : TRANS_ID_NONE;

    ctx_add_new_zf_ctx_from_txtheader (dict_id, id, id_len, luft_trans);

    return luft_trans;
}

// If Liftover doesn't exist - adds it according to Genozip's capabilities and also sets ctx->trans_luft
// If it does exist - verifies that it is recognized by Genozip, and sets ctx->trans_luft accordingly
static void vcf_header_update_INFO_FORMAT (const char *line, unsigned line_len, bool is_info, Buffer *new_header)
{
    const char *number, *id, *alg;
    unsigned number_len=0, id_len=0, alg_len=0;
    TranslatorId luft_trans;

    vcf_header_get_attribute (line, line_len, is_info ? 7 : 9, "Number=", 7, false, &number, &number_len);
    vcf_header_get_attribute (line, line_len, is_info ? 7 : 9, "ID=",     3, false, &id, &id_len);
    vcf_header_get_attribute (line, line_len, is_info ? 7 : 9, HK_RENDERALG_ATTR"=", sizeof HK_RENDERALG_ATTR-1 +1, false, &alg, &alg_len);

    // case: liftover algorithm is specified in header, see that Genozip supports it and set context, or error if not 
    if (alg) {
        SAFE_NUL (&alg[alg_len]);
        for (luft_trans=0; luft_trans < NUM_VCF_TRANS; luft_trans++)
            if (!strcmp (alg, ltrans_props[luft_trans].alg_name)) 
                goto found;

        ABORT ("Unrecongized Liftover value \"%s\". Offending line:\n%.*s", alg, line_len, line);
    
        found: 
        SAFE_RESTORE;
        vcf_header_set_luft_trans (id, id_len, is_info, luft_trans, 0, 0);
        
        buf_add_more (NULL, new_header, line, line_len, NULL); // unchanged line
    }

    // case: liftover algorithm not specified - get Genozip's default for this ID
    else {
        luft_trans = vcf_header_set_luft_trans (id, id_len, is_info, TRANS_ID_UNKNOWN, number, number_len);

        vcf_header_add_line_add_key (new_header, line, line_len, HK_RENDERALG_ATTR, ltrans_props[luft_trans].alg_name);
    }
}

// --chain: convert one line at a time to dual-coordinate format, and add all additional lines before #CHROM
static bool vcf_header_zip_update_to_dual_coords_one_line (const char *line, unsigned line_len, void *new_header_p, void *unused2, unsigned unused3)
{
    Buffer *new_header = (Buffer *)new_header_p;
    bool is_info;
    
    if ((is_info = LINEIS ("##INFO=")) || LINEIS ("##FORMAT=")) 
        vcf_header_update_INFO_FORMAT (line, line_len, is_info, new_header);
    
    else if (LINEIS ("##fileformat")) {
        // remove this line - it was added to the rejects file instead
    }
    else {
        if (LINEIS ("#CHROM") && chain_is_loaded)
            vcf_header_zip_convert_header_to_dc_add_lines (new_header);

        buf_add_more (evb, new_header, line, line_len, NULL); // just add the line, unmodified
    }

    return false; // continue iterating
}

static int dict_id_to_ID_sorter (const void *a, const void *b) { return ((DictIdToID *)a)->dict_id.num - ((DictIdToID *)b)->dict_id.num; }

// ZIP with --chain: add liftover_contig, liftover_reference, chain, dual_coordinates keys
static void vcf_header_zip_update_to_dual_coords (Buffer *txt_header)
{
    #define new_txt_header evb->codec_bufs[0]

    if (z_file->num_txt_components_so_far == 0)  // first component of this z_file - allocation by previous z_file
        buf_free (&dict_id_to_ID);
    
    // a mapping of full ID names of FORMAT / INFO field for inclusion in LIFTREJT fields (we can't use 
    // ctx->name as its only 8 characters) 
    buf_alloc (evb, &dict_id_to_ID, 0, 100, DictIdToID, 0, "dict_id_to_ID");
    buf_alloc (evb, &new_txt_header, 0, txt_header->len + 262144, char, 0, "evb->codec_bufs[0]"); // initial allocation

    txtfile_foreach_line (txt_header, false, vcf_header_zip_update_to_dual_coords_one_line, &new_txt_header, 0, 0, 0);

    buf_free (txt_header);
    buf_copy (evb, txt_header, &new_txt_header, char, 0, 0, "txt_header");
    buf_free (&new_txt_header);

    // sort by dict_id.num
    qsort (dict_id_to_ID.data, dict_id_to_ID.len, sizeof (DictIdToID), dict_id_to_ID_sorter);

    #undef new_txt_header
}

static const char *vcf_header_get_VCF_ID_by_dict_id_do (DictId dict_id, int first, int last, bool must_exist)
{
    if (first > last) { // not found
        ASSERT (!must_exist, "Cannot find ID=%s in VCF header", dis_dict_id_ex (dict_id, true).s);
        return NULL;
    }
    
    int mid = (first + last) / 2;
    DictIdToID *ent = ENT (DictIdToID, dict_id_to_ID, mid);

    if (dict_id.num == ent->dict_id.num) return ent->ID;

    if (dict_id.num > ent->dict_id.num) return vcf_header_get_VCF_ID_by_dict_id_do (dict_id, mid+1, last, must_exist);
    else                                return vcf_header_get_VCF_ID_by_dict_id_do (dict_id, first, mid-1, must_exist);
}

// returns the ID of an INFO or FORMAT field as it appears in the VCF header - nul-terminated
const char *vcf_header_get_VCF_ID_by_dict_id (DictId dict_id, bool must_exist)
{
    return vcf_header_get_VCF_ID_by_dict_id_do (dict_id, 0, (int)(dict_id_to_ID.len)-1, must_exist);
}
 
static int dst_contigs_sorter (const void *a, const void *b) { return *(WordIndex *)a - *(WordIndex *)b; }

// scan header for contigs - used to pre-populate z_file contexts in txtheader_zip_prepopulate_contig_ctxs
// in --chain - also builds vcf_header_liftover_dst_contigs 
static void vcf_header_handle_contigs (Buffer *txt_header)
{
    ASSERTNOTINUSE (vcf_header_liftover_dst_contigs);

    if (chain_is_loaded && !flag.rejects_coord)
        buf_alloc (evb, &vcf_header_liftover_dst_contigs, 0, 100, WordIndex, 0, "codec_bufs[1]");

    txtheader_verify_contig_init();
    txtfile_foreach_line (txt_header, false, vcf_header_extract_contigs, &vcf_header_liftover_dst_contigs, 0, 0, 0);

    // sort (but not uniq)
    if (chain_is_loaded && !flag.rejects_coord) 
        qsort (vcf_header_liftover_dst_contigs.data, vcf_header_liftover_dst_contigs.len, sizeof (WordIndex), dst_contigs_sorter);
}

static bool vcf_inspect_txt_header_zip (Buffer *txt_header)
{
    if (!vcf_header_set_globals (txt_file->name, txt_header, true)) return false; // samples are different than a previous concatented file

    // scan header for ##dual_coordinates - this sets txt_file->coords
    txtfile_foreach_line (txt_header, false, vcf_header_get_dual_coords, 0, 0, 0, 0);
    
    // handle contig lines
    vcf_header_handle_contigs (txt_header);

    // if --chain: scan header for ##FORMAT and ##INFO - for translatable subfields - create contexts, and populate Context.luft_trans
    if (chain_is_loaded) flag.bind = BIND_REJECTS;

    // if we're compressing a txt file with liftover, we prepare the two rejects headers: it consists of two lines:
    // 1. the first "##fileformat" line if one exists in our file 
    // 2. the last "#CHROM" line - needing for zipping the rejects file, but will be removed in reconstruction
    // in reconstruction without --luft, we will reconstruct the normal header 
    // in reconstruction with --luft, we will reconstruct the rejected header first, then the normal header without the first line
    if ((chain_is_loaded || txt_file->coords) && !flag.rejects_coord) {
        buf_add_more (evb, &evb->lo_rejects[0], txt_header->data, vcf_header_get_first_line_len (txt_header), "lo_rejects");
        buf_add_more (evb, &evb->lo_rejects[1], txt_header->data, vcf_header_get_first_line_len (txt_header), "lo_rejects");
        
        char *line;
        unsigned len = vcf_header_get_last_line (txt_header, &line);
        buf_add_more (0, &evb->lo_rejects[0], line, len, 0);
        buf_add_more (0, &evb->lo_rejects[1], line, len, 0);

        vcf_lo_append_rejects_file (evb, DC_PRIMARY);
        vcf_lo_append_rejects_file (evb, DC_LUFT);
    }

    // case: Luft file - update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
    if (txt_file->coords == DC_LUFT)
        vcf_header_zip_liftback_header (txt_header);
    else if (txt_file->coords == DC_PRIMARY) 
        vcf_header_zip_get_luft_only_lines (txt_header);

    // case: --chain: add liftover_contig, liftover_reference, chain, dual_coordinates keys
    if ((chain_is_loaded || txt_file->coords) && !flag.rejects_coord)
        vcf_header_zip_update_to_dual_coords (txt_header);

    return true; // all good
}

static bool vcf_inspect_txt_header_piz (VBlock *txt_header_vb, Buffer *txt_header, struct FlagsTxtHeader txt_header_flags)
{
    vcf_header_set_globals (z_file->name, txt_header, false);

    if (flag.genocat_no_reconstruct) return true;

    // if genocat --samples is used, update vcf_header and vcf_num_displayed_samples
    // note: we subset reject samples (in the header) as well - we analyze only in the rejects section
    if (flag.samples) vcf_header_subset_samples (&vcf_field_name_line);
    else              vcf_num_displayed_samples = vcf_num_samples;

    // for the rejects part of the header - we just remove its #CHROM line as it will be taken from the primary component
    if (txt_header_flags.rejects_coord) {
        txt_header->len -= vcf_field_name_line.len;
        return true;
    }

    // remove #CHROM line (it is saved in vcf_field_name_line by vcf_header_set_globals()) - or everything
    // if we ultimately one want the #CHROM line. We will add back the #CHROM line later.
    if (flag.header_one) 
        txt_header->len = 0; 
    else
        txt_header->len -= vcf_field_name_line.len; 

    // if needed, liftover the header
    if (flag.luft && !flag.header_one && !flag.rejects_coord) 
        vcf_header_piz_liftover_header (txt_header_vb, txt_header);

    // add genozip command line
    if (!flag.header_one && exe_type == EXE_GENOCAT && !flag.genocat_no_reconstruct && !flag.rejects_coord
        && !flag.no_pg && (flag.data_modified || z_dual_coords)) 
        vcf_header_add_genozip_command (txt_header_vb, txt_header);

    if (flag.drop_genotypes) 
        vcf_header_trim_field_name_line (&vcf_field_name_line); // drop FORMAT and sample names

    // if --show_dvcf, add two extra fields at beginning for field name line, and remove the # from #CHROM
    if (flag.show_dvcf)
        #define ADDITIONAL_FIELDS "#COORD\tSTATUS\t"
        buf_add_more (txt_header_vb, txt_header, ADDITIONAL_FIELDS, sizeof ADDITIONAL_FIELDS-1, "txt_data");

    // add the (perhaps modified) field name (#CHROM) line
    buf_add_more (txt_header_vb, txt_header, &vcf_field_name_line.data[!!flag.show_dvcf], vcf_field_name_line.len - !!flag.show_dvcf, "txt_data");

    return true; // all good
}

bool vcf_inspect_txt_header (VBlockP txt_header_vb, Buffer *txt_header, struct FlagsTxtHeader txt_header_flags)
{
    return (command == ZIP) ? vcf_inspect_txt_header_zip (txt_header)
                            : vcf_inspect_txt_header_piz (txt_header_vb, txt_header, txt_header_flags);
}

static bool vcf_header_set_globals (const char *filename, Buffer *vcf_header, bool soft_fail)
{
    static const char *vcf_field_name_line_filename = NULL; // file from which the header line was taken

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
            if (!buf_is_alloc (&vcf_field_name_line) || command == PIZ) {
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
                      vcf_field_name_line_filename, (uint32_t)vcf_field_name_line.len, vcf_field_name_line.data,
                      filename, (uint32_t)vcf_header->len-i, &vcf_header->data[i]);
                
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
static void vcf_header_trim_field_name_line (Buffer *vcf_header)
{
    char *line;
    unsigned line_len = vcf_header_get_last_line (vcf_header, &line);

    if (is_field_name_line (line, line_len)) {
        vcf_header->len = strstr (line, "INFO") + 5 - FIRSTENT (char, *vcf_header);
        *LASTENT (char, *vcf_header) = '\n';
    }                  
}

uint32_t vcf_header_get_num_samples (void)
{
    if (z_file->data_type == DT_VCF || z_file->data_type == DT_BCF)
        return vcf_num_samples;
    else
        return 0;
}

// -------------
// samples stuff
// -------------

// processes the vcf header sample line according to the --samples option, removing samples that are not required
// and building a bytemap PIZ filter the samples
static void vcf_header_subset_samples (Buffer *vcf_field_name_line)
{
    // accept a sample from the vcf file's samples as consistent with the --samples requested
    #define samples_accept(sample_str) { \
        buf_add_string (evb, vcf_field_name_line, sample_str); \
        buf_add_string (evb, vcf_field_name_line, "\t"); \
        vcf_num_displayed_samples++; \
    }

    ARRAY (char, line, *vcf_field_name_line);

    flag.samples = is_field_name_line (line, line_len);
    RETURNW (flag.samples,, "Warning: found non-standard VCF sample header line. Ingoring --samples : \n%.*s", (int)line_len, line);

    int32_t num_samples=-8;
    for (unsigned i=0; i < line_len; i++)
        if (line[i] == '\t') num_samples++;
        
    vcf_samples_is_included = MALLOC (num_samples);
    memset (vcf_samples_is_included, cmd_is_negative_samples, num_samples); // 0 if not included unless list says so (positive) and vice versa

    unsigned vcf_names_start_index = (line + sizeof VCF_FIELD_NAMES_LONG-1 + 1) - FIRSTENT (char, *vcf_field_name_line);
    unsigned vcf_names_data_len = vcf_field_name_line->len - vcf_names_start_index;
    vcf_sample_names_data = MALLOC (vcf_names_data_len);
    memcpy (vcf_sample_names_data, &vcf_field_name_line->data[vcf_names_start_index], vcf_names_data_len);
    vcf_sample_names_data[vcf_names_data_len-1] = '\t'; // change last separator from \n to \t

    vcf_sample_names = MALLOC (num_samples * sizeof (char *));

    vcf_field_name_line->len = vcf_names_start_index;
    char *next_token = vcf_sample_names_data;
    
    // go through the vcf file's samples and add those that are consistent with the --samples requested
    vcf_num_displayed_samples = 0;

    #define working_buf evb->codec_bufs[0]

    buf_copy (evb, &working_buf, &cmd_samples_buf, char*, 0, 0, 0);

    for (unsigned i=0; i < num_samples; i++) {
        vcf_sample_names[i] = strtok_r (next_token, "\t", &next_token);
        bool handled = false;

        for (unsigned s=0; s < working_buf.len; s++) 
            if (!strcmp (vcf_sample_names[i], *ENT(char *, working_buf, s))) { // found
                vcf_samples_is_included[i] = !cmd_is_negative_samples;
                if (!cmd_is_negative_samples) 
                    samples_accept (vcf_sample_names[i]);

                // remove this sample from the --samples list as we've found it already
                memcpy (ENT(char *, working_buf, s), ENT(char *, working_buf, s+1), (working_buf.len-s-1) * sizeof (char *));
                working_buf.len--;
                handled = true;
                break;
            }
            
        if (!handled && cmd_is_negative_samples) 
            samples_accept (vcf_sample_names[i]);
    }

    *LASTENT (char, *vcf_field_name_line) = '\n'; // change last separator from \t

    // warn about any --samples items that were not found in the vcf file (all the ones that still remain in the buffer)
    for (unsigned s=0; s < working_buf.len; s++) 
        ASSERTW (false, "Warning: requested sample '%s' is not found in the VCF file, ignoring it", *ENT(char *, working_buf, s));

    // if the user filtered out all samples, its equivalent of drop_genotypes
    if (!vcf_num_displayed_samples) flag.drop_genotypes = true;

    buf_free (&working_buf);
    #undef working_buf
}

// called from genozip.c for processing the --samples flag
void vcf_samples_add  (const char *samples_str)
{
    ASSERTNOTNULL (samples_str);

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
        for (unsigned s=0; s < cmd_samples_buf.len; s++)
            if (!strcmp (one_sample, *ENT (char *, cmd_samples_buf, s))) {
                is_duplicate = true;
                break;
            }
        if (is_duplicate) continue; // skip duplicates "genocat -s sample1,sample2,sample1"

        buf_alloc_old (evb, &cmd_samples_buf, MAX (cmd_samples_buf.len + 1, 100) * sizeof (char*), 2, "cmd_samples_buf");

        NEXTENT (char *, cmd_samples_buf) = one_sample;
    }
}


