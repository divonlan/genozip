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

// Globals
uint32_t vcf_num_samples           = 0; // number of samples in the file
static uint32_t vcf_num_displayed_samples = 0; // PIZ only: number of samples to be displayed - might be less that vcf_num_samples if --samples is used
static Buffer vcf_field_name_line = EMPTY_BUFFER;  // header line of first VCF file read - use to compare to subsequent files to make sure they have the same header during bound

// referring to sample strings from the --samples command line option
static Buffer cmd_samples_buf = EMPTY_BUFFER; // an array of (char *)
static bool cmd_is_negative_samples = false;

// referring to samples in the vcf file
char *vcf_samples_is_included;         // a bytemap indicating for each sample if it is included
static char **vcf_sample_names;        // an array of char * to nul-terminated names of samples 
static char *vcf_sample_names_data;    // vcf_sample_names point into here

#define LINEIS(s) (line_len > strlen(s) && !memcmp (line, (s), strlen(s))) // note: gcc optimizer resolves strlen("CONST") in compilation time

#define SUBST_LABEL(old,new) do { buf_add_more (NULL, new_txt_header, (new), strlen (new), NULL); \
                                  buf_add_more (NULL, new_txt_header, line + strlen (old), line_len - strlen (old), NULL); } while (0)

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
    buf_add_string (txt_header_vb, txt_header, flags_command_line()->data);
    bufprintf (txt_header_vb, txt_header, "\" %s\n", str_time().s);
}

static void vcf_header_get_subvalue (const char *line, unsigned line_len, unsigned key_len,
                                     const char *subkey, unsigned subkey_len, bool enforce,// in
                                     const char **snip, unsigned *snip_len) // out
{
    const char *start = strstr (line + key_len, subkey);
    if (!start) goto missing; // subkey doesn't exist
    start += subkey_len;

    const char *comma   = memchr (start, ',', &line[line_len] - start);
    const char *bracket = memchr (start, '>', &line[line_len] - start);

    if (!comma && !bracket) goto missing; // subkey isn't terminated with a comma or >
    const char *after = !comma ? bracket : !bracket ? comma : MIN (comma,bracket);
    
    *snip = start;
    *snip_len = after - start;
    return;

missing:
    ASSINP (!enforce, "missing subkey %s in header line %.*s", subkey, line_len, line);
}

// add contig to buffer if line is a valid contig, and fail silently if not
static void vcf_header_parse_contig_line (const char *line, unsigned line_len, bool liftover)
{
    const char *id=0, *len_str=0;
    unsigned id_len=0, len_str_len=0;
    PosType length = 0; 

    #define HK_CONTIG "contig=<"
    #define HK_CONTIG_LEN 8

    // parse eg "##contig=<ID=NKLS02001838.1,length=29167>" and "##liftover_contig=<ID=22>"
    SAFE_NUL (&line[line_len]);
    if (!strstr (line, HK_CONTIG)) goto done;

    vcf_header_get_subvalue (line, line_len, HK_CONTIG_LEN, "ID=", 3, true, &id, &id_len);
    vcf_header_get_subvalue (line, line_len, HK_CONTIG_LEN, "length=", 7, false, &len_str, &len_str_len);

    str_get_int_range64 (len_str, len_str_len, 1, 1000000000000ULL, &length); // length stays 0 if len_str_len=0

    // if its the 2nd+ file while binding in ZIP - we just check that the contigs are the same
    // (or a subset) of previous file's contigs (to do: merge headers, bug 327) 
    if (command == ZIP && flag.bind && z_file->num_txt_components_so_far /* 2nd+ file */) 
        txtheader_verify_contig (id, id_len, length, (void *)liftover);

    else {
        txtheader_alloc_contigs (1, id_len+1, liftover);
        txtheader_add_contig (id, id_len, length, (void *)liftover);
    }

done:
    SAFE_RESTORE;
}

static bool vcf_header_extract_contigs (const char *line, unsigned line_len, void *unused1, void *unused2, unsigned unused3)
{
    if ((LINEIS ("##contig=") && txt_file->dual_coords == DC_LUFT) || 
         LINEIS (HK_LO_CONTIG))
        vcf_header_parse_contig_line (line, line_len, true);

    else if (LINEIS ("##contig=") || LINEIS (HK_LB_CONTIG)) 
        vcf_header_parse_contig_line (line, line_len, false);

    return false; // continue iterating
}


static bool vcf_header_get_dual_coords (const char *line, unsigned line_len, void *unused1, void *unused2, unsigned unused3)
{
    if (LINEIS (HK_DC_PRIMARY)) 
        txt_file->dual_coords = DC_PRIMARY;

    else if (LINEIS (HK_DC_LUFT)) {
        txt_file->dual_coords = DC_LUFT;
        flag.data_modified = true; // we will be zipping this as a dual-coordinates VCF with default reconstruction as PRIMARY - this turns off digest
    }

    else return false; // continue iterating

    // genozip of a dual coordinates file implies --sort, unless overridden with --unsorted
    if (!flag.unsorted) flag.sort = true;

    return true; // found - no more iterating
}

// PIZ with --luft - one line: remove fileformat, update *contig, *reference, dual_coordinates keys to Luft format
static bool vcf_header_piz_liftover_header_one_line (const char *line, unsigned line_len, 
                                                     void *new_txt_header_, void *unused1, unsigned unused2)
{
    Buffer *new_txt_header = (Buffer *)new_txt_header_;

    if      (LINEIS (HK_LO_CONTIG))   SUBST_LABEL (HK_LO_CONTIG, "##contig=");
    else if (LINEIS ("##contig="))    SUBST_LABEL ("##contig=", HK_LB_CONTIG);
    else if (LINEIS (HK_LO_REF))      SUBST_LABEL (HK_LO_REF, "##reference=");
    else if (LINEIS ("##reference=")) SUBST_LABEL ("##reference=", HK_LB_REF);
    else if (LINEIS (HK_DC_PRIMARY))  SUBST_LABEL (HK_DC_PRIMARY, HK_DC_LUFT);
    else if (LINEIS ("##fileformat=")) {}  // remove ##fileformat as it will be reconstructed from the rejects component
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

// ZIP of Luft line: update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
static bool vcf_header_zip_liftback_header_one_line (const char *line, unsigned line_len, 
                                                     void *new_txt_header_, void *rejects_, unsigned unused)
{
    Buffer *new_txt_header = (Buffer *)new_txt_header_;
    
    if      (LINEIS (HK_LB_REJECT))   buf_add_more (evb, (Buffer *)rejects_, line + strlen (HK_LB_REJECT), line_len - strlen (HK_LB_REJECT), NULL); // // add rejects to "rejects" buffer
    else if (LINEIS (HK_LB_CONTIG))   SUBST_LABEL (HK_LB_CONTIG, "##contig=");
    else if (LINEIS ("##contig="))    SUBST_LABEL ("##contig=", HK_LO_CONTIG);
    else if (LINEIS (HK_LB_REF))      SUBST_LABEL (HK_LB_REF, "##reference=");
    else if (LINEIS ("##reference=")) SUBST_LABEL ("##reference=", HK_LO_REF);
    else if (LINEIS (HK_DC_LUFT))     SUBST_LABEL (HK_DC_LUFT, HK_DC_PRIMARY);
    else                              buf_add_more (evb, new_txt_header, line, line_len, NULL);

    return false; // continue iterating
}

// ZIP of Luft file: update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
static void vcf_header_zip_liftback_header (Buffer *txt_header)
{
    #define new_txt_header evb->codec_bufs[0]
    #define rejects        evb->codec_bufs[1]

    // case: Luft file - scan the header and for LIFTBACK_REJECT lines, and move them to txt_file->unconsumed (after removing the label) 
    // to be processsed as normal variants
    buf_alloc_old (evb, &new_txt_header, txt_header->len + 1 /*+1 to make sure alloc allocates */, 0, "evb->codec_bufs[0]"); // initial allocation (might be a bit bigger due to label changes)
    buf_alloc_old (evb, &rejects, 65536, 0, "evb->codec_bufs[1]"); // initial allocation

    txtfile_foreach_line (txt_header, false, vcf_header_zip_liftback_header_one_line, &new_txt_header, &rejects, 0, 0);

    // replace txt_header with lifted back one
    buf_free (txt_header);
    buf_copy (evb, txt_header, &new_txt_header, char, 0, 0, "txt_data");

    // squeeze in the liftover rejects from the header, with their label removed, before the existing unconsumed text
    // (data that was read beyond the header) 
    if (rejects.len) {
        buf_alloc (evb, &txt_file->unconsumed_txt, rejects.len, 0, char, 0, "txt_file->unconsumed_txt");
        ARRAY (char, unconsumed_txt, txt_file->unconsumed_txt);
        
        memmove (&unconsumed_txt[rejects.len], unconsumed_txt, unconsumed_txt_len);
        memcpy (unconsumed_txt, rejects.data, rejects.len);
        
        txt_file->unconsumed_txt.len += rejects.len;
        txt_file->luft_reject_bytes = rejects.len;
    }

    buf_free (&new_txt_header);
    buf_free (&rejects);

    #undef new_txt_header
    #undef rejects
}

// ZIP with --chain: add liftover_contig, liftover_reference, chain, dual_coordinates keys
static void vcf_header_zip_convert_header_to_dual_coords (Buffer *txt_header)
{
    txt_header->len -= vcf_field_name_line.len; // remove field name line

    bufprintf (evb, txt_header, "%s\n", HK_DC_PRIMARY);
    bufprintf (evb, txt_header, HK_CHAIN"file://%s\n", chain_filename ? chain_filename : "<null>");
    bufprintf (evb, txt_header, HK_LO_REF"file://%s\n", ref_filename ? ref_filename : "<null>");

    bufprintf (evb, txt_header, KH_INFO_LO"\n", GENOZIP_CODE_VERSION);
    bufprintf (evb, txt_header, KH_INFO_LB"\n", GENOZIP_CODE_VERSION);
    bufprintf (evb, txt_header, KH_INFO_LR"\n", GENOZIP_CODE_VERSION);

    // add all Luft contigs (note: we get them from liftover directly, as they zip_initialize that adds them to zf_ctx has not been called yet)
    const char *luft_contig;
    PosType length;
    for (unsigned i=0; (luft_contig = chain_get_dst_contig(i, &length)); i++) 
        if (length) bufprintf (evb, txt_header, HK_LO_CONTIG"<ID=%s,length=%"PRId64">\n", luft_contig, length);
        else        bufprintf (evb, txt_header, HK_LO_CONTIG"<ID=%s>\n", luft_contig);

    buf_add_buf (evb, txt_header, &vcf_field_name_line, char, "txt_header"); // careful not to use bufprintf as it is unbound
}

static bool vcf_inspect_txt_header_zip (Buffer *txt_header)
{
    if (!vcf_header_set_globals (txt_file->name, txt_header, true)) return false; // samples are different than a previous concatented file

    // scan header for ##dual_coordinates - this sets txt_file->dual_coords
    txtfile_foreach_line (txt_header, false, vcf_header_get_dual_coords, 0, 0, 0, 0);
    
    // scan header for contigs - used to pre-populate z_file contexts in txtheader_zip_prepopulate_contig_ctxs
    txtheader_verify_contig_init();
    txtfile_foreach_line (txt_header, false, vcf_header_extract_contigs, 0, 0, 0, 0);

    ASSERT (!chain_is_loaded || !txt_file->dual_coords, 
            "--chain cannot be used with %s - it is already a dual-coordinates VCF file - it contains \""HK_DC"\" in its header", txt_name);

    // if we're compressing a txt file with liftover, we prepare the rejects header: it consists of two lines:
    // 1. the first "##fileformat" line if one exists in our file 
    // 2. the last "#CHROM" line - needing for zipping the rejects file, but will be removed in reconstruction
    // in reconstruction without --luft, we will reconstruct the normal header 
    // in reconstruction with --luft, we will reconstruct the rejected header first, then the normal header without the first line
    if ((chain_is_loaded /* --chain */ || txt_file->dual_coords /* dual coordinates file */) && !flag.processing_rejects) {
        buf_add_more (evb, &evb->liftover_rejects, txt_header->data, vcf_header_get_first_line_len (txt_header), "liftover_rejects");
        
        char *line;
        unsigned len = vcf_header_get_last_line (txt_header, &line);
        buf_add_more (evb, &evb->liftover_rejects, line, len, "liftover_rejects");
        liftover_append_rejects_file (evb);
    }

    // case: --chain: add liftover_contig, liftover_reference, chain, dual_coordinates keys
    if (chain_is_loaded && !flag.processing_rejects)
        vcf_header_zip_convert_header_to_dual_coords (txt_header);

    // case: Luft file - update *contig, *reference, dual_coordinates keys to Primary format, and move rejects to unconsumed_txt
    else if (txt_file->dual_coords == DC_LUFT)
        vcf_header_zip_liftback_header (txt_header);

    return true; // all good
}

static bool vcf_inspect_txt_header_piz (VBlock *txt_header_vb, Buffer *txt_header)
{
    vcf_header_set_globals (z_file->name, txt_header, false);

    if (flag.genocat_no_reconstruct) return true;

    // if genocat --samples is used, update vcf_header and vcf_num_displayed_samples
    // note: we subset reject samples (in the header) as well - we analyze only in the rejects section
    if (flag.samples) vcf_header_subset_samples (&vcf_field_name_line);
    else              vcf_num_displayed_samples = vcf_num_samples;

    // for the rejects part of the header - we just remove its #CHROM line as it will be taken from the primary component
    if (flag.processing_rejects) {
        txt_header->len -= vcf_field_name_line.len;
        return true;
    }

    // remove #CHROM line (it is saved in vcf_field_name_line by vcf_header_set_globals()) - or everything
    // if we ultimately one want the #CHROM line. We will add back the #CHROM line later.
    if (flag.header_one) 
        txt_header->len = 0; 
    else
        txt_header->len -= vcf_field_name_line.len; 
    
    // in case of --luft, we remove the first line of the regular file header (##fileformat), as they are 
    // duplicate in both files been output as part of the rejects file
    if (flag.luft && !flag.header_one && !flag.genocat_no_reconstruct) {
        unsigned len = vcf_header_get_first_line_len (txt_header);
        memmove (txt_header->data, txt_header->data + len, txt_header->len - len);
        txt_header->len -= len;
    }

    // if needed, liftover the header
    if (flag.luft && !flag.header_one && !flag.processing_rejects) 
        vcf_header_piz_liftover_header (txt_header_vb, txt_header);

    // add genozip command line
    if (!flag.header_one && exe_type == EXE_GENOCAT && !flag.genocat_no_reconstruct && !flag.processing_rejects
        && !flag.no_pg && (flag.data_modified || z_file->z_flags.dual_coords)) 
        vcf_header_add_genozip_command (txt_header_vb, txt_header);

    if (flag.drop_genotypes) 
        vcf_header_trim_field_name_line (&vcf_field_name_line); // drop FORMAT and sample names

    // add the (perhaps modified) header
    buf_add_buf (txt_header_vb, txt_header, &vcf_field_name_line, char, "txt_data");

    return true; // all good
}

bool vcf_inspect_txt_header (VBlockP txt_header_vb, Buffer *txt_header)
{
    return (command == ZIP) ? vcf_inspect_txt_header_zip (txt_header)
                            : vcf_inspect_txt_header_piz (txt_header_vb, txt_header);
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
        
            // if first vcf file - copy the header to the global
            if (!buf_is_alloc (&vcf_field_name_line)) {
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

            ASSINP (tab_count >= 7, "Error: invalid field header line - it contains only %d fields, expecting at least 8", tab_count+1);

            return true; 
        }
    }

    ABORT_R ("Error: invalid VCF file - it does not contain a field header line; tab_count=%u", tab_count+1);
}

// genocat: remove FORMAT and sample names from the vcf header line, in case of --drop-genotypes
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

    unsigned vcf_names_start_index = (line + strlen(VCF_FIELD_NAMES_LONG) + 1) - FIRSTENT (char, *vcf_field_name_line);
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


