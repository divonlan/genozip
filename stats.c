// ------------------------------------------------------------------
//   stats.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "buffer.h"
#include "dict_id.h"
#include "strings.h"
#include "stats.h"
#include "sections.h"
#include "file.h"
#include "vblock.h"
#include "txtfile.h"
#include "reference.h"
#include "zfile.h"
#include "version.h"

static int *count_per_section = NULL;

// a concatenation of all bound txt_names that contributed to this genozip file
static Buffer bound_txt_names = EMPTY_BUFFER;

#define OVERHEAD_SEC_VB_HDR      -1
#define OVERHEAD_SEC_ALIASES     -2
#define OVERHEAD_SEC_TXT_HDR     -3
#define OVERHEAD_SEC_RA_INDEX    -4
#define OVERHEAD_SEC_REFERENCE   -5 // used only for REF_EXT_STORE
#define OVERHEAD_SEC_REF_HASH    -6 // used only in a reference file
#define OVERHEAD_SEC_REF_CONTIGS -7 
#define OVERHEAD_SEC_ALT_CHROMS  -8 
#define LAST_OVERHEAD            -8

static void stats_get_sizes (DictId dict_id /* option 1 */, int overhead_sec /* option 2*/, 
                             int64_t *dict_compressed_size, int64_t *b250_compressed_size, int64_t *local_compressed_size)
{
    *dict_compressed_size = *b250_compressed_size = *local_compressed_size = 0;

    for (unsigned i=0; i < z_file->section_list_buf.len; i++) {

        SectionListEntry *section = ENT (SectionListEntry, z_file->section_list_buf, i);
        int64_t after_sec = (i == z_file->section_list_buf.len - 1) ? z_file->disk_so_far : (section+1)->offset;
        int64_t sec_size = after_sec - section->offset;

        count_per_section[i]++; // we're optimistically assuming section will be count_per_section - we will revert if not
        
        // we put all the SEQ derived data in a single stats line for easy comparison (reference in dict, the +- section SEQ in b250, and the unique non-ref sequence stretchs in local)
        if (z_file->data_type == DT_SAM && dict_id.num == dict_id_fields[SAM_SEQ_BITMAP]) { // 'SEQ' line - must be before SEC_LOCAL
            if (section->dict_id.num == dict_id_fields[SAM_NONREF])
                *local_compressed_size += sec_size;
            else if ((section->section_type == SEC_REFERENCE || section->section_type == SEC_REF_IS_SET) && flag_reference == REF_INTERNAL)
                *dict_compressed_size += sec_size;
            else if (section->dict_id.num == dict_id_fields[SAM_SEQ_BITMAP])
                *b250_compressed_size += sec_size;
        }

        else if (overhead_sec == OVERHEAD_SEC_REFERENCE && 
                 (flag_reference == REF_EXT_STORE || z_file->data_type != DT_SAM) && 
                 (section->section_type == SEC_REFERENCE || section->section_type == SEC_REF_IS_SET))
            *dict_compressed_size += sec_size;

        else if (z_file->data_type == DT_SAM && dict_id.num == dict_id_fields[SAM_NONREF]) {} // ^ already handled above

        else if (section->dict_id.num == dict_id.num && section->section_type == SEC_DICT)
            *dict_compressed_size += sec_size;

        else if (section->dict_id.num == dict_id.num && section->section_type == SEC_B250)
            *b250_compressed_size += sec_size;

        else if ((section->dict_id.num == dict_id.num && section->section_type == SEC_LOCAL) ||
                 (section->section_type == SEC_VB_HEADER       && overhead_sec == OVERHEAD_SEC_VB_HDR) ||
                 (section->section_type == SEC_DICT_ID_ALIASES && overhead_sec == OVERHEAD_SEC_ALIASES) ||
                 (section->section_type == SEC_REF_CONTIGS     && overhead_sec == OVERHEAD_SEC_REF_CONTIGS) ||
                 (section->section_type == SEC_ALT_CHROMS      && overhead_sec == OVERHEAD_SEC_ALT_CHROMS) ||
                 (section->section_type == SEC_REF_HASH        && overhead_sec == OVERHEAD_SEC_REF_HASH) ||
                 (section->section_type == SEC_TXT_HEADER      && overhead_sec == OVERHEAD_SEC_TXT_HDR) ||
                 (section->section_type == SEC_RANDOM_ACCESS   && overhead_sec == OVERHEAD_SEC_RA_INDEX) ||
                 (section->section_type == SEC_REF_RAND_ACC    && overhead_sec == OVERHEAD_SEC_RA_INDEX))
            *local_compressed_size += sec_size;
             
        else count_per_section[i]--; // acually... not count_per_section!
    }
}

static void stats_check_count (uint64_t all_comp_total)
{
    if (all_comp_total == z_file->disk_so_far) return; // all good

    ARRAY (SectionListEntry, sections, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++) 
        if (!count_per_section[i]) 
            fprintf (stderr, "Section not counted: %s section_i=%u\n", st_name (sections[i].section_type), i);
        else if (count_per_section[i] > 1) 
            fprintf (stderr, "Section overcounted: %s section_i=%u dict=%.8s counted %u times\n", st_name (sections[i].section_type), i, dict_id_printable (sections[i].dict_id).id, count_per_section[i]);

    char s1[20], s2[20];
    ASSERTW (false, "Hmm... incorrect calculation for GENOZIP sizes: total section sizes=%s but file size is %s (diff=%d)", 
             str_uint_commas (all_comp_total, s1), str_uint_commas (z_file->disk_so_far, s2), (int32_t)(z_file->disk_so_far - all_comp_total));
}

static void stats_show_file_metadata (Buffer *buf)
{
    bufprintf (evb, buf, "%s", "\n\n");
    if (txt_file->name) 
        bufprintf (evb, buf, "%s file%s: %.*s\n", dt_name (z_file->data_type), 
                   memchr (bound_txt_names.data, ' ', bound_txt_names.len) ? "s" : "", // a space separator indicates more than one file
                   (int)bound_txt_names.len, bound_txt_names.data);
    
    if (flag_reference == REF_INTERNAL && z_file->data_type == DT_SAM) {
        bufprintf (evb, buf, "Reference: Internal ('SEQ' line: 'comp dict'=stored reference, 'comp 250'=bitmap, 'comp local'=non-refereance bases)%s\n", "");
    }
    else if (flag_reference == REF_EXTERNAL) {
        if (z_file->data_type == DT_SAM)
            bufprintf (evb, buf, "Reference: %s ('SEQ' line: 'comp 250'=bitmap, 'comp local'=non-refereance bases)\n", ref_filename)
        else
            bufprintf (evb, buf, "Reference: %s (not stored in genozip file)\n", ref_filename);
    }
    else if (flag_reference == REF_INTERNAL || flag_reference == REF_EXT_STORE)
        bufprintf (evb, buf, "Reference: %s (size appears in 'Reference' line)\n", ref_filename);

    char ls[30];
    if (z_file->data_type == DT_VCF) 
        bufprintf (evb, buf, "Samples: %u   ", vcf_header_get_num_samples());

    bufprintf (evb, buf, "%s: %s   Dictionaries: %u   Vblocks: %u   Sections: %u\n", 
               DTPZ (show_stats_line_name), str_uint_commas (z_file->num_lines, ls), z_file->num_contexts, 
               z_file->num_vbs, (uint32_t)z_file->section_list_buf.len);

    char timestr[100];
    time_t now = time (NULL);
    strftime (timestr, 100, "%Y-%m-%d %H:%M:%S", gmtime (&now));
    bufprintf (evb, buf, "Genozip version: %s Date compressed: %s UTC\n", GENOZIP_CODE_VERSION, timestr);
}

typedef struct {
    int64_t total_comp_size;
    char did_i[20], name[50], type[20], words[20], hash[20], uncomp_dict[20], comp_dict[20],
         comp_b250[20], comp_data[20], comp_total[20], txt[20];
         double comp_ratio, pc_txt, pc_genozip, pc_dict, pc_singletons, pc_failed_singletons, pc_hash_occupancy;
} StatsByLine;


static int stats_sort_by_total_comp_size(const void *a, const void *b)  
{ 
    return (((StatsByLine*)b)->total_comp_size > ((StatsByLine*)a)->total_comp_size) ? 1 : -1; // use comparison (>) and not minus (-) as the retun value is only 32 bit
}

// generate the stats text - all sections except genozip header and the stats section itself 
void stats_compress (void)
{
    stats_show_file_metadata(&z_file->stats_buf_1);
    buf_copy (evb, &z_file->stats_buf_2, &z_file->stats_buf_1, 0,0,0, "z_file->stats_buf_2", 0);

    int64_t all_comp_dict=0, all_uncomp_dict=0, all_comp_b250=0, all_comp_data=0, all_comp_total=0, all_txt=0;

    //prepare data
    StatsByLine sbl[MAX_DICTS+50], *s = sbl;
    memset (sbl, 0, sizeof(sbl));

    count_per_section = calloc (z_file->section_list_buf.len, sizeof (int));

    for (int i=LAST_OVERHEAD; i < (int)z_file->num_contexts; i++) { 

        Context *ctx = (i>=0) ? &z_file->contexts[i] : NULL;
   
        int64_t dict_compressed_size, b250_compressed_size, local_compressed_size, txt_len;
        stats_get_sizes (ctx ? ctx->dict_id : DICT_ID_NONE, ctx ? 0 : i, 
                         &dict_compressed_size, &b250_compressed_size, &local_compressed_size);

        s->total_comp_size = dict_compressed_size + b250_compressed_size + local_compressed_size;

        ASSERTW (s->total_comp_size >= 0, "Hmm... s->total_comp_size=%"PRId64" is negative for overhead=%d", s->total_comp_size, i);

        if (ctx && !ctx->mtf_i.len && !ctx->txt_len && !ctx->b250.len && !s->total_comp_size) continue;

        txt_len = ctx ? ctx->txt_len : (i==OVERHEAD_SEC_TXT_HDR ? txtfile_get_last_header_len() : 0);

        bool is_dict_a_reference = (i == OVERHEAD_SEC_REFERENCE) || // reference appreas in "comp_dict" of "Reference" line
                                   (ctx && z_file->data_type == DT_SAM && (ctx->dict_id.num == dict_id_fields[SAM_SEQ_BITMAP] || ctx->dict_id.num == dict_id_fields[SAM_NONREF])); // reference appreas in "comp_dict" of "SEQ" line
        
        if (!is_dict_a_reference) all_comp_dict += dict_compressed_size;
        all_uncomp_dict += ctx ? ctx->dict.len : 0;
        all_comp_b250   += b250_compressed_size;
        all_comp_data   += local_compressed_size;
        all_comp_total  += s->total_comp_size;
        all_txt         += txt_len;

        /* Name         */ 
        if (ctx)                              strcpy (s->name, ctx->name); 
        else if (i==OVERHEAD_SEC_ALIASES)     strcpy (s->name, "aliases");
        else if (i==OVERHEAD_SEC_VB_HDR)      strcpy (s->name, "VBlock hdrs");
        else if (i==OVERHEAD_SEC_TXT_HDR)     strcpy (s->name, "Txt file header");
        else if (i==OVERHEAD_SEC_RA_INDEX)    strcpy (s->name, "--regions index");
        else if (i==OVERHEAD_SEC_REFERENCE)   strcpy (s->name, "Reference");
        else if (i==OVERHEAD_SEC_REF_HASH)    strcpy (s->name, "Refhash");
        else if (i==OVERHEAD_SEC_REF_CONTIGS) strcpy (s->name, "Ref Contigs");
        else if (i==OVERHEAD_SEC_ALT_CHROMS)  strcpy (s->name, "Alt Chroms");

        /* Type           */ strcpy (s->type, ctx ? dict_id_display_type (z_file->data_type, ctx->dict_id) : "OTHER");

        if (ctx) {
        /* did_i          */ str_uint_commas ((uint64_t)ctx->did_i, s->did_i); 
        /* #Words in file */ str_uint_commas (ctx->mtf_i.len, s->words);
        /* % dict         */ s->pc_dict              = !ctx->mtf_i.len         ? 0 : 100.0 * (double)ctx->mtf.len / (double)ctx->mtf_i.len;
        /* % singletons   */ s->pc_singletons        = !ctx->mtf_i.len         ? 0 : 100.0 * (double)ctx->num_singletons / (double)ctx->mtf_i.len;
        /* % failed singl.*/ s->pc_failed_singletons = !ctx->mtf_i.len         ? 0 : 100.0 * (double)ctx->num_failed_singletons / (double)ctx->mtf_i.len;
        /* % hash occupn. */ s->pc_hash_occupancy    = !ctx->global_hash_prime ? 0 : 100.0 * (double)(ctx->mtf.len + ctx->ol_mtf.len) / (double)ctx->global_hash_prime;
        /* Hash           */ str_uint_commas (ctx->global_hash_prime, s->hash);
        /* uncomp dict    */ str_size (ctx->dict.len, s->uncomp_dict);
        /* comp dict      */ str_size (dict_compressed_size, s->comp_dict);
        }
        else 
            s->did_i[0] = s->words[0] = s->hash[0] = s->uncomp_dict[0] = '-';
        
        if (ctx || is_dict_a_reference)
        /* comp dict      */ str_size (dict_compressed_size, s->comp_dict);
        else 
            s->comp_dict[0] = '-';

        /* txt            */ str_size (txt_len, s->txt);
        /* comp b250      */ str_size (b250_compressed_size, s->comp_b250);
        /* comp data      */ str_size (local_compressed_size, s->comp_data);
        /* comp total     */ str_size (s->total_comp_size, s->comp_total);
        /* comp ratio     */ s->comp_ratio = (double)txt_len / (double)s->total_comp_size;
        /* % of txt       */ s->pc_txt     = 100.0 * (double)txt_len / (double)z_file->txt_data_so_far_bind;
        /* % of genozip   */ s->pc_genozip = 100.0 * (double)s->total_comp_size / (double)z_file->disk_so_far;

        s++;
    }
    unsigned num_stats = s - sbl;

    // sort by total compressed size
    qsort (sbl, num_stats, sizeof (sbl[0]), stats_sort_by_total_comp_size);

    bufprintf (evb, &z_file->stats_buf_1, "\nSections (sorted by %% of genozip file):%s\n", "");
    bufprintf (evb, &z_file->stats_buf_1, "NAME                  TXT      %%       ZIP      %%   RATIO%s\n", "");

    bufprintf (evb, &z_file->stats_buf_2, "\nSections (sorted by %% of genozip file):%s\n", "");
    bufprintf (evb, &z_file->stats_buf_2, "did_i Name            Type            #Words  Snips-(%% of #Words)        Hash-table    uncomp      comp      comp      comp      comp       txt    comp   %% of   %% of  %s\n", "");
    bufprintf (evb, &z_file->stats_buf_2, "                                     in file   Dict  Local   Both         Size Occp      dict      dict      b250     local     TOTAL             ratio    txt    zip%s\n", "");

    for (uint32_t i=0; i < num_stats; i++) { // don't show primary fields as they are already showed above

        StatsByLine *s = &sbl[i];
        if (!s->total_comp_size) continue;

#define PC(pc) ((pc==0 || pc>=10) ? 0 : (pc<1 ? 2:1))
        bufprintf (evb, &z_file->stats_buf_1, "%-15.15s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                   s->name, s->txt, s->pc_txt, s->comp_total, s->pc_genozip, s->comp_ratio);

        bufprintf (evb, &z_file->stats_buf_2, "%-2.2s    %-15.15s %-6.6s %15s  %4.*f%%  %4.*f%%  %4.*f%% %12s %3.0f%% %9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
                    s->did_i, s->name, s->type, s->words, 
                    PC (s->pc_dict), s->pc_dict, PC(s->pc_singletons), s->pc_singletons, PC(s->pc_failed_singletons), s->pc_failed_singletons, 
                    s->hash, s->pc_hash_occupancy, // Up to here - these don't appear in the total
                    s->uncomp_dict, s->comp_dict, s->comp_b250, s->comp_data, s->comp_total, s->txt, s->comp_ratio, s->pc_txt, s->pc_genozip);
    }

    double all_comp_ratio = (double)all_txt / (double)all_comp_total;
    double all_pc_txt     = 100.0 * (double)all_txt        / (double)z_file->txt_data_so_far_bind;
    double all_pc_genozip = 100.0 * (double)all_comp_total / (double)z_file->disk_so_far;
    
    char s1[20], s2[20], s3[20], s4[20], s5[20], s6[20];

    bufprintf (evb, &z_file->stats_buf_1, "TOTAL           "
                "%9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                str_size (all_txt, s6), all_pc_txt, str_size (all_comp_total, s5), all_pc_genozip, all_comp_ratio);

    bufprintf (evb, &z_file->stats_buf_2, "TOTAL                                                                               "
                "%9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
                str_size (all_uncomp_dict, s1), str_size (all_comp_dict, s2),  str_size (all_comp_b250, s3), 
                str_size (all_comp_data, s4),   str_size (all_comp_total, s5), str_size (all_txt, s6), 
                all_comp_ratio, all_pc_txt, all_pc_genozip);

    stats_check_count (all_comp_total);

    // note: we use txt_data_so_far_single and not txt_data_size_single, because the latter has estimated size if disk_so_far is 
    // missing, while txt_data_so_far_single is what was actually processed
    ASSERTW (all_txt == z_file->txt_data_so_far_bind || flag_optimize || flag_make_reference, 
             "Hmm... incorrect calculation for %s sizes: total section sizes=%s but file size is %s (diff=%d)", 
             dt_name (z_file->data_type), str_uint_commas (all_txt, s1), str_uint_commas (z_file->txt_data_so_far_bind, s2), 
             (int32_t)(z_file->txt_data_so_far_bind - all_txt)); 

    zfile_compress_section_data (evb, SEC_STATS, &z_file->stats_buf_1);
    zfile_compress_section_data (evb, SEC_STATS, &z_file->stats_buf_2);

    FREE (count_per_section);
}

void stats_display (void)
{
    Buffer *buf = flag_show_stats == 1 ? &z_file->stats_buf_1 : &z_file->stats_buf_2;

    if (!buf_is_allocated (buf)) return; // no stats available

    buf_print (buf , false);

    SectionListEntry *sl = sections_get_first_section_of_type (SEC_STATS, false);

    if (z_file->disk_size < (1<<20)) { // no need to print this note if total size > 1MB, as the ~2K of overhead is rounded off anyway
        char s[20];
        // stats text doesn't include SEC_STATS and SEC_GENOZIP_HEADER - the last 2 sections in the file - since stats text is generated before these sections are compressed
        printf ("\nNote: ZIP total file size excludes overhead of %s\n", str_size (z_file->disk_size - sl->offset, s));
    }

    printf ("\n");
}

void stats_read_and_display (void)
{
    SectionListEntry *sl = sections_get_first_section_of_type (SEC_STATS, true);
    if (!sl) return; // genozip file does not contain stats sections (SEC_STATS was introduced in v 7.0.5)

    // read and uncompress the requested stats section
    zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", sizeof (SectionHeader), SEC_STATS, sl + (flag_show_stats==2));
    zfile_uncompress_section (evb, evb->z_data.data, flag_show_stats == 1 ? &z_file->stats_buf_1 : &z_file->stats_buf_2, "z_file->stats_buf", 0, SEC_STATS);
    buf_free (&evb->z_data);
    
    stats_display();

    if (exe_type == EXE_GENOCAT) exit_ok; // if this is genocat - we're done
}

// concatenate txt names of bound files so we can show them all
void stats_add_txt_name (const char *fn)
{
    unsigned fn_len = strlen (fn);

    buf_alloc (evb, &bound_txt_names, bound_txt_names.len + fn_len + 1, 2, "bound_txt_names", 0);
    
    if (bound_txt_names.len) buf_add (&bound_txt_names, " ", 1);
    buf_add (&bound_txt_names, fn, fn_len);
}


