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

#define working_buf evb->codec_bufs[0] // so we don't override cmd_samples_buf that may be used for the next files

void vcf_header_piz_init (void)
{
    vcf_num_samples = 0;
    buf_free (&vcf_field_name_line);

    // note: we don't re-initialize vcf_num_displayed_samples - this is calculated only once
}


static void vcf_header_subset_samples (Buffer *vcf_header);
static bool vcf_header_set_globals (const char *filename, Buffer *vcf_header, bool soft_fail);
static void vcf_header_trim_field_name_line (Buffer *vcf_header);

// false means continue iterating, true means stop
typedef bool (*VcfHeaderIteratorCallback)(const char *line, unsigned line_len, const void *cb_param1, unsigned cb_param2);
static char *vcf_header_foreach_line (Buffer *txt_header,
                                      bool reverse, // iterate backwards
                                      VcfHeaderIteratorCallback callback, const void *cb_param1, unsigned cb_param2,
                                      int64_t *line_len) // out
{
    if (line_len) *line_len = 0;

    if (!txt_header->len) return NULL;

    char *firstbuf = txt_header->data;
    char *afterbuf = AFTERENT (char, *txt_header);

    char *first = !reverse ? firstbuf : 0;
    char *after = !reverse ? 0 : afterbuf;

    while (1) {
            
        // get one line - searching forward or backwards
        if (!reverse) {
            for (after=first ; after < afterbuf && *after != '\n' ; after++);
            after++; // skip newline
        }
        else {
            for (first=after-2 /* skip final \n */; first >= firstbuf && *first != '\n'; first--);
            first++; // after detected \n or at start of line
        }

        if (!reverse && after > afterbuf) return NULL; // we don't call callback if after>afterbuf - beyond end of line
            
        if (callback (first, after - first, cb_param1, cb_param2)) {
            if (line_len) *line_len = after - first;
            return first;
        }

        if (reverse && first == firstbuf) return NULL; // beginning of line - we called the cb

        if (!reverse) first=after;
        else          after=first;
    }

    return 0; // never reaches here
}   

// false means continue iterating, true means stop
static bool vcf_header_is_key_cb (const char *line, unsigned line_len, const void *key, unsigned key_len)
{
    return line_len >= key_len && !memcmp (key, line, key_len);
}

static bool vcf_header_one_line_cb (const char *line, unsigned line_len, const void *key, unsigned key_len)
{
    return true;
}

// replaces a line that starts with key with new_line, or if none exists, appends new_line to the end of the buffer
static void vcf_header_add_or_replace_line (Buffer *txt_header, const char *key, const char *new_line) 
{
    int64_t old_line_len;
    int64_t new_line_len = strlen (new_line);
    char *line = vcf_header_foreach_line (txt_header, true, vcf_header_is_key_cb, key, strlen (key), &old_line_len);

    if (!line) line = AFTERENT (char, *txt_header);

    uint64_t index = ENTNUM (*txt_header, line);

    if (new_line_len > old_line_len)
        buf_alloc_more (evb, txt_header, new_line_len - old_line_len, 0, char, 1.5, "txt_data");

    line = ENT (char, *txt_header, index); // reset after realloc

    // make room (or shrink it)
    memmove (line + new_line_len, line + old_line_len, AFTERENT (char, *txt_header) - line - old_line_len);

    // assign new data
    memcpy (line, new_line, new_line_len);

    txt_header->len += new_line_len - old_line_len; // possitive or negative
}

// returns the length of the first line, if it starts with ##fileformat, or 0
static inline unsigned vcf_header_get_first_line_len (const Buffer *txt_header)
{
    if (txt_header->len < 12 || memcmp (txt_header->data, "##fileformat", 12)) return 0;

    const char *newline = memchr (txt_header->data, '\n', txt_header->len);
    
    return newline ? newline - txt_header->data + 1 : 0;
}

static inline bool is_field_name_line (const char *line, unsigned line_len)
{
    unsigned field_name_line_len = strlen (VCF_FIELD_NAMES);
    return line_len >= field_name_line_len && !memcmp (line, VCF_FIELD_NAMES, field_name_line_len);
}

// returns the length of the last line, if it starts with #CHROM, or 0
static unsigned vcf_header_get_last_line (Buffer *txt_header, char **line_p)
{
    int64_t line_len;
    char *line = vcf_header_foreach_line (txt_header, true, vcf_header_one_line_cb, 0, 0, &line_len);
    if (line_p) *line_p = line;

    return is_field_name_line (line, line_len) ? (unsigned)line_len : 0;
}

static void vcf_header_add_dual_coord_keys (Buffer *txt_header)
{
    #define LIFTOVER "##INFO=<ID=" INFO_LIFTOVER ",Number=.,Type=String,Description=\"Dual-coordinate VCF: Information for lifting over the variant to laft coordinates\",Source=\"genozip\",Version=\"%s\">\n"
    #define LIFTBACK "##INFO=<ID=" INFO_LIFTBACK ",Number=.,Type=String,Description=\"Dual-coordinate VCF: Information for retrieving the variant in the primary coordinates\",Source=\"genozip\",Version=\"%s\">\n"
    #define LIFTREJD "##INFO=<ID=" INFO_LIFTREJD ",Number=.,Type=String,Description=\"Dual-coordinate VCF: Reason variant was rejected for lift over\",Source=\"genozip\",Version=\"%s\">\n"

    // add INFO lines only if they're not there
    char line[strlen(LIFTOVER)+50];

    sprintf (line, HEADER_KEY_DC"%s\n", flag.laft ? HEADER_KEY_DC_LAFT : HEADER_KEY_DC_PRIMARY);
    vcf_header_add_or_replace_line (txt_header, HEADER_KEY_DC, line);

    sprintf (line, LIFTOVER, GENOZIP_CODE_VERSION);
    vcf_header_add_or_replace_line (txt_header, "##INFO=<ID=" INFO_LIFTOVER, line);
    
    sprintf (line, LIFTBACK, GENOZIP_CODE_VERSION);
    vcf_header_add_or_replace_line (txt_header, "##INFO=<ID=" INFO_LIFTBACK, line);

    sprintf (line, LIFTREJD, GENOZIP_CODE_VERSION);
    vcf_header_add_or_replace_line (txt_header, "##INFO=<ID=" INFO_LIFTREJD, line);
}

static void vcf_header_add_genozip_command (Buffer *txt_header)
{
    // the command line length is unbound, careful not to put it in a bufprintf
    buf_add_string (evb, txt_header, "##genozip_Command=\"");
    buf_add_string (evb, txt_header, flags_command_line()->data);
    bufprintf (evb, txt_header, "\" %s\n", str_time().s);
}

static bool vcf_header_extract_one_liftover_reject (const char *line, unsigned line_len, const void *first_reject_line, unsigned unused)
{
    // if not a reject line - return
    unsigned label_len = strlen (HEADER_KEY_LIFTOVER_REJECT);
    if (line_len <= label_len || memcmp (line, HEADER_KEY_LIFTOVER_REJECT, label_len))
        // we reached the end the ##LIFTOVER_REJECTS area if we already observed a liftover line, otherwise he haven't reached it
        // yet - continue iterating
        return (bool)*(const char **)first_reject_line; 

    buf_add_more (evb, &working_buf, line + label_len, line_len - label_len, NULL);

    if (! *(const char **)first_reject_line) 
        *(const char **)first_reject_line = line; 

    return false; // continue iterating
}

static bool vcf_header_get_dual_coords (const char *line, unsigned line_len, const void *unused1, unsigned unused2)
{
    // if not a DUAL_COORDINATES line - return and continue iterating
    unsigned label_len = strlen (HEADER_KEY_DC);
    unsigned prim_len  = strlen (HEADER_KEY_DC_PRIMARY);
    unsigned laft_len  = strlen (HEADER_KEY_DC_LAFT);

    if (line_len <= label_len || memcmp (line, HEADER_KEY_DC, label_len)) return false; // continue iterating

    const char *value = &line[label_len];
    unsigned value_len = line_len - label_len - 1; 
    if (value[value_len-1] == '\r') value_len--;

    if (value_len == prim_len && !memcmp (value, HEADER_KEY_DC_PRIMARY, value_len))
        txt_file->dual_coords = DC_PRIMARY;
    
    else if (value_len == laft_len && !memcmp (value, HEADER_KEY_DC_LAFT, value_len)) {
        txt_file->dual_coords = DC_LAFT;
        flag.data_modified = true; // we will be zipping this as a dual-coordinate VCF with default reconstruction as PRIMARY - this turns off digest
    }
    
    else
        ABORT ("Expecting \""HEADER_KEY_DC"\" header line to be "HEADER_KEY_DC_PRIMARY" or "HEADER_KEY_DC_LAFT", but it is %.*s", 
               value_len, value);

    return true; // found - not more iterating
}

static bool vcf_inspect_txt_header_zip (Buffer *txt_header)
{
    if (!vcf_header_set_globals (txt_file->name, txt_header, true)) return false; // samples are different than a previous concatented file

    // scan header for ##DUAL_COORDINATES
    vcf_header_foreach_line (txt_header, false, vcf_header_get_dual_coords, 0, 0, 0);

    // if we're compressing a txt file with liftover, we prepare the rejects header: it consists of two lines:
    // 1. the first "##fileformat" line if one exists in our file 
    // 2. the last "#CHROM" line - needing for zipping the rejects file, but will be removed in reconstruction
    // in reconstruction without --laft, we will reconstruct the normal header 
    // in reconstruction with --laft, we will reconstruct the rejected header first, then the normal header without the first line
    if ((chain_is_loaded /* --chain */ || txt_file->dual_coords /* dual coordinates file */) && !flag.processing_rejects) {
        buf_add_more (evb, &evb->liftover_rejects, txt_header->data, vcf_header_get_first_line_len (txt_header), "liftover_rejects");
        
        char *line;
        unsigned len = vcf_header_get_last_line (txt_header, &line);
        buf_add_more (evb, &evb->liftover_rejects, line, len, "liftover_rejects");
        liftover_append_rejects_file (evb);
    }

    // scan the header for LIFTOVER_REJECT lines, and move them to txt_file->unconsumed (after removing the lable) to be processsed as 
    // normal variants
    buf_alloc (evb, &working_buf, 500000, 0, "evb->codec_bufs[0]"); // initial allocation
    char *first_reject=0, *after_rejects=0;
    after_rejects = vcf_header_foreach_line (txt_header, false, vcf_header_extract_one_liftover_reject, &first_reject, 0, 0);

    // if we moved some LIFTOVER_REJECT to working_buf, we shrink the txt_header to remove that area
    if (first_reject) { 
        uint64_t rejects_bytes_inc_label = after_rejects - first_reject;
        uint64_t after_rejects_bytes = AFTERENT (char, *txt_header) - after_rejects;

        // remove rejects from header
        memmove (first_reject, after_rejects, after_rejects_bytes);
        txt_header->len -= rejects_bytes_inc_label;

        // squeeze in the liftover rejects from the header, with their label removed, befoer the existing unconsumed text
        // (data that was read beyond the header) 
        buf_alloc_more (evb, &txt_file->unconsumed_txt, working_buf.len, 0, char, 0, "txt_file->unconsumed_txt");
        ARRAY (char, unconsumed_txt, txt_file->unconsumed_txt);
        
        memmove (&unconsumed_txt[working_buf.len], unconsumed_txt, unconsumed_txt_len);
        memcpy (unconsumed_txt, working_buf.data, working_buf.len);
        
        txt_file->unconsumed_txt.len += working_buf.len;
        txt_file->laft_reject_bytes = working_buf.len;
    }

    buf_free (&working_buf);
    return true; // all good
}

static bool vcf_inspect_txt_header_piz (Buffer *txt_header)
{
    vcf_header_set_globals (z_file->name, &evb->txt_data, false);

    if (flag.genocat_no_reconstruct) return true;

    // if genocat --samples is used, update vcf_header and vcf_num_displayed_samples
    // note: we subset reject samples (in the header) as well - we analyze only in the rejects section
    if (flag.samples) vcf_header_subset_samples (&vcf_field_name_line);
    else              vcf_num_displayed_samples = vcf_num_samples;

    // for the rejects part of the header - we just remove its #CHROM line as it will be taken from the primary component
    if (flag.processing_rejects) {
        txt_header->len -= vcf_header_get_last_line (txt_header, 0); 
        return true;
    }

    // remove #CHROM line (it is saved in vcf_field_name_line by vcf_header_set_globals()) - or everything
    // if we ultimately one want the #CHROM line. We will add back the #CHROM line later.
    if (flag.header_one) 
        txt_header->len = 0; 
    else
        txt_header->len -= vcf_header_get_last_line (txt_header, 0); 
    
    // in case of --laft, we remove the first line of the regular file header (##fileformat), as they are 
    // duplicate in both files been output as part of the rejects file
    if (!flag.header_one && flag.laft) {
        unsigned len = vcf_header_get_first_line_len (txt_header);
        memmove (txt_header->data, txt_header->data + len, txt_header->len - len);
        txt_header->len -= len;
    }

    // add dual-coordinate keys to the header data
    if (z_file->z_flags.dual_coords && exe_type == EXE_GENOCAT && !flag.header_one && !flag.genocat_no_reconstruct && 
        !flag.processing_rejects)
        vcf_header_add_dual_coord_keys (txt_header);

    // add genozip command line
    if (!flag.header_one && exe_type == EXE_GENOCAT && !flag.genocat_no_reconstruct && !flag.processing_rejects
        && !flag.no_pg && (flag.data_modified || z_file->z_flags.dual_coords)) 
        vcf_header_add_genozip_command (txt_header);

    if (flag.drop_genotypes) 
        vcf_header_trim_field_name_line (&vcf_field_name_line); // drop FORMAT and sample names

    // add the (perhaps modified) header
    buf_add_more (evb, txt_header, vcf_field_name_line.data, vcf_field_name_line.len, "txt_data");

    return true; // all good
}

bool vcf_inspect_txt_header (Buffer *txt_header)
{
    return (command == ZIP) ? vcf_inspect_txt_header_zip (txt_header)
                            : vcf_inspect_txt_header_piz (txt_header);
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
            if (!buf_is_allocated (&vcf_field_name_line)) {
                buf_copy (evb, &vcf_field_name_line, vcf_header, 1, i, vcf_header->len - i, "vcf_field_name_line");
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

            ASSINP (tab_count >= 7, "Error: invalid VCF file - field header line contains only %d fields, expecting at least 8", tab_count+1);

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

    buf_copy (evb, &working_buf, &cmd_samples_buf, sizeof (char*), 0, 0, 0);

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
}

// called from genozip.c for processing the --samples flag
void vcf_samples_add  (const char *samples_str)
{
    ASSERTE0 (samples_str, "samples_str is NULL");

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

        buf_alloc (evb, &cmd_samples_buf, MAX (cmd_samples_buf.len + 1, 100) * sizeof (char*), 2, "cmd_samples_buf");

        NEXTENT (char *, cmd_samples_buf) = one_sample;
    }
}


