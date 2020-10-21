// ------------------------------------------------------------------
//   stats.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <stdarg.h>
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

static void stats_get_sizes (DictId dict_id /* option 1 */, SectionType non_ctx_sec /* option 2*/, 
                             int64_t *dict_compressed_size, int64_t *b250_compressed_size, int64_t *local_compressed_size)
{
    *dict_compressed_size = *b250_compressed_size = *local_compressed_size = 0;

    for (uint64_t i=0; i < z_file->section_list_buf.len; i++) {

        SectionListEntry *section = ENT (SectionListEntry, z_file->section_list_buf, i);
        int64_t after_sec = (i == z_file->section_list_buf.len - 1) ? z_file->disk_so_far : (section+1)->offset;
        int64_t sec_size = after_sec - section->offset;

        count_per_section[i]++; // we're optimistically assuming section will be count_per_section - we will revert if not

        if (section->dict_id.num == dict_id.num && section->section_type == SEC_DICT)
            *dict_compressed_size += sec_size;

        else if (section->dict_id.num == dict_id.num && section->section_type == SEC_B250)
            *b250_compressed_size += sec_size;

        else if ((section->dict_id.num == dict_id.num && section->section_type == SEC_LOCAL) ||
                 (section->section_type == non_ctx_sec))
            *local_compressed_size += sec_size;
             
        else count_per_section[i]--; // acually... not count_per_section!
    }
}

static void stats_check_count (uint64_t all_z_size)
{
    if (all_z_size == z_file->disk_so_far) return; // all good

    ARRAY (SectionListEntry, sections, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++) 
        if (!count_per_section[i]) 
            fprintf (stderr, "Section not counted: %s section_i=%u\n", st_name (sections[i].section_type), i);
        else if (count_per_section[i] > 1) 
            fprintf (stderr, "Section overcounted: %s section_i=%u dict=%.8s counted %u times\n", st_name (sections[i].section_type), i, dict_id_printable (sections[i].dict_id).id, count_per_section[i]);

    char s1[20], s2[20];
    ASSERTW (false, "Hmm... incorrect calculation for GENOZIP sizes: total section sizes=%s but file size is %s (diff=%d)", 
             str_uint_commas (all_z_size, s1), str_uint_commas (z_file->disk_so_far, s2), (int32_t)(z_file->disk_so_far - all_z_size));
}

static void stats_show_file_metadata (Buffer *buf)
{
    bufprintf (evb, buf, "%s", "\n\n");
    if (txt_file->name) 
        bufprintf (evb, buf, "%s file%s: %.*s\n", dt_name (z_file->data_type), 
                   memchr (bound_txt_names.data, ' ', bound_txt_names.len) ? "s" : "", // a space separator indicates more than one file
                   (int)bound_txt_names.len, bound_txt_names.data);
    
    if (flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE) 
        bufprintf (evb, buf, "Reference: %s\n", ref_filename);

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
    int64_t txt_size, z_size;
    const char *name;
    char did_i[20], type[20], words[20], hash[20], uncomp_dict[20], comp_dict[20],
    comp_b250[20], comp_data[20];
    double pc_of_txt, pc_of_z, pc_dict, pc_singletons, pc_failed_singletons, pc_hash_occupancy;
} StatsByLine;

static int stats_sort_by_z_size(const void *a, const void *b)  
{ 
    return (((StatsByLine*)b)->z_size > ((StatsByLine*)a)->z_size) ? 1 : -1; // use comparison (>) and not minus (-) as the retun value is only 32 bit
}

static void stats_consolidate_compound (StatsByLine *sbl, unsigned num_stats, const char *consolidated_name, const char *field)
{
    unsigned field_len = strlen (field);

    StatsByLine *survivor = NULL;    
    for (uint32_t i=0; i < num_stats; i++)
        if (sbl[i].name && !strcmp (sbl[i].name, field)) {
            survivor = &sbl[i];
            survivor->name = consolidated_name;
        }

    if (!survivor) return; // this data type doesn't have this compound field

    for (unsigned i=0; i < num_stats; i++) {
        
        // check for eg DESC -> D1ESC
        if (sbl[i].name && (field_len+1 == strlen (sbl[i].name)) && (field[0] == sbl[i].name[0]) && !strcmp (&field[1], &sbl[i].name[2])) {
            survivor->txt_size  += sbl[i].txt_size;
            survivor->z_size    += sbl[i].z_size;  
            survivor->pc_of_txt += sbl[i].pc_of_txt;
            survivor->pc_of_z   += sbl[i].pc_of_z;  

            sbl[i] = (StatsByLine){};
        }
    }
}

static void stats_consolidate_related (StatsByLine *sbl, unsigned num_stats, const char *consolidated_name, unsigned num_deps, ...)
{
    va_list args;
    va_start (args, num_deps);

    const char *deps[num_deps];
    for (unsigned d=0; d < num_deps; d++) 
        deps[d] = va_arg (args, const char *);

    StatsByLine *survivor = NULL;
    for (unsigned i=0; i < num_stats; i++) {
        
        if (!sbl[i].name) continue; // unused entry

        for (unsigned d=0; d < num_deps; d++)
            if (!strcmp (deps[d], sbl[i].name)) {
                if (!survivor) {
                    survivor = &sbl[i]; // first found gets to be the survivor
                }
                else {
                    survivor->txt_size  += sbl[i].txt_size;
                    survivor->z_size    += sbl[i].z_size;  
                    survivor->pc_of_txt += sbl[i].pc_of_txt;
                    survivor->pc_of_z   += sbl[i].pc_of_z; 
                    survivor->name       = consolidated_name; // rename only if at least one was consolidated

                    sbl[i] = (StatsByLine){}; 
                }
                break;
            } 
    }

    va_end (args);
}

static void stats_output_stats (StatsByLine *s, unsigned num_stats, 
                                int64_t all_txt_size, int64_t all_z_size, double all_pc_of_txt, double all_pc_of_z, double all_comp_ratio)
{
    char s1[20], s2[20];

    bufprintf (evb, &z_file->stats_buf_1, "\nSections (sorted by %% of genozip file):%s\n", "");
    bufprintf (evb, &z_file->stats_buf_1, "NAME                  TXT      %%       ZIP      %%   RATIO%s\n", "");

    for (uint32_t i=0; i < num_stats; i++, s++)
        if (s->z_size)
            bufprintf (evb, &z_file->stats_buf_1, "%-15.15s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                       s->name, str_size (s->txt_size, s1), s->pc_of_txt, str_size (s->z_size, s2), s->pc_of_z, 
                       (double)s->txt_size / (double)s->z_size);
    
    bufprintf (evb, &z_file->stats_buf_1, "TOTAL           "
                "%9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                str_size (all_txt_size, s1), all_pc_of_txt, str_size (all_z_size, s2), all_pc_of_z, all_comp_ratio);
}

static void stats_output_STATS (StatsByLine *s, unsigned num_stats,
                                int64_t all_txt_size, int64_t all_uncomp_dict, int64_t all_comp_dict, int64_t all_comp_b250, int64_t all_comp_data, 
                                int64_t all_z_size, double all_pc_of_txt, double all_pc_of_z, double all_comp_ratio)
{
#define PC(pc) ((pc==0 || pc>=10) ? 0 : (pc<1 ? 2:1))
    char s1[20], s2[20], s3[20], s4[20], s5[20], s6[20];

    bufprintf (evb, &z_file->stats_buf_2, "\nSections (sorted by %% of genozip file):%s\n", "");
    bufprintf (evb, &z_file->stats_buf_2, "did_i Name            Type            #Words  Snips-(%% of #Words)        Hash-table    uncomp      comp      comp      comp      comp       txt    comp   %% of   %% of  %s\n", "");
    bufprintf (evb, &z_file->stats_buf_2, "                                     in file   Dict  Local   Both         Size Occp      dict      dict      b250     local     TOTAL             ratio    txt    zip%s\n", "");

    for (uint32_t i=0; i < num_stats; i++, s++)
        if (s->z_size)
            bufprintf (evb, &z_file->stats_buf_2, "%-2.2s    %-15.15s %-6.6s %15s  %4.*f%%  %4.*f%%  %4.*f%% %12s %3.0f%% %9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
                        s->did_i, s->name, s->type, s->words, 
                        PC (s->pc_dict), s->pc_dict, PC(s->pc_singletons), s->pc_singletons, PC(s->pc_failed_singletons), s->pc_failed_singletons, 
                        s->hash, s->pc_hash_occupancy, // Up to here - these don't appear in the total
                        s->uncomp_dict, s->comp_dict, s->comp_b250, s->comp_data, str_size (s->z_size, s2), 
                        str_size (s->txt_size, s1), (double)s->txt_size / (double)s->z_size, s->pc_of_txt, s->pc_of_z);

    bufprintf (evb, &z_file->stats_buf_2, "TOTAL                                                                               "
                "%9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
                str_size (all_uncomp_dict, s1), str_size (all_comp_dict, s2),  str_size (all_comp_b250, s3), 
                str_size (all_comp_data, s4),   str_size (all_z_size, s5), str_size (all_txt_size, s6), 
                all_comp_ratio, all_pc_of_txt, all_pc_of_z);
}

// generate the stats text - all sections except genozip header and the stats section itself 
void stats_compress (void)
{
    stats_show_file_metadata(&z_file->stats_buf_1);
    buf_copy (evb, &z_file->stats_buf_2, &z_file->stats_buf_1, 0,0,0, "z_file->stats_buf_2", 0);

    int64_t all_comp_dict=0, all_uncomp_dict=0, all_comp_b250=0, all_comp_data=0, all_z_size=0, all_txt_size=0;

    //prepare data
    StatsByLine sbl[MAX_DICTS + NUM_SEC_TYPES] = { }, *s = sbl;

    count_per_section = CALLOC (z_file->section_list_buf.len * sizeof (int));

    #define SEC(i) (i<0 ? -(i)-1 : SEC_NONE) // i to section type
    #define ST_NAME(st) (&st_name(st)[4]) // cut off "SEC_" 

    for (int i=-NUM_SEC_TYPES; i < (int)z_file->num_contexts; i++) { // sections go into -1 to -NUM_SEC_TYPES (see SEC())

        if (SEC(i) == SEC_DICT || SEC(i) == SEC_B250 || SEC(i) == SEC_LOCAL) continue; // these are covered by indiviual contexts

        Context *ctx = (i>=0) ? &z_file->contexts[i] : NULL;
   
        int64_t dict_compressed_size, b250_compressed_size, local_compressed_size;
        stats_get_sizes (ctx ? ctx->dict_id : DICT_ID_NONE, SEC(i), 
                         &dict_compressed_size, &b250_compressed_size, &local_compressed_size);

        s->z_size = dict_compressed_size + b250_compressed_size + local_compressed_size;

        ASSERTW (s->z_size >= 0, "Hmm... s->z_size=%"PRId64" is negative for %s", s->z_size, s->name);

        if (ctx && !ctx->mtf_i.len && !ctx->txt_len && !ctx->b250.len && !s->z_size) continue;

        s->txt_size = ctx ? ctx->txt_len : (i==-SEC_TXT_HEADER ? txtfile_get_bound_headers_len() : 0);
        
        all_comp_dict   += dict_compressed_size;
        all_uncomp_dict += ctx ? ctx->dict.len : 0;
        all_comp_b250   += b250_compressed_size;
        all_comp_data   += local_compressed_size;
        all_z_size      += s->z_size;
        all_txt_size    += s->txt_size;

        /* Name           */ s->name = ctx ? ctx->name : ST_NAME (SEC(i)); 
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
        
        if (ctx)
        /* comp dict      */ str_size (dict_compressed_size, s->comp_dict);
        else 
            s->comp_dict[0] = '-';

        /* comp b250      */ str_size (b250_compressed_size, s->comp_b250);
        /* comp data      */ str_size (local_compressed_size, s->comp_data);
        /* % of txt       */ s->pc_of_txt = 100.0 * (double)s->txt_size / (double)z_file->txt_data_so_far_bind;
        /* % of genozip   */ s->pc_of_z   = 100.0 * (double)s->z_size / (double)z_file->disk_so_far;

        s++;
    }
    unsigned num_stats = s - sbl;

    double all_comp_ratio = (double)all_txt_size / (double)all_z_size;
    double all_pc_of_txt  = 100.0 * (double)all_txt_size        / (double)z_file->txt_data_so_far_bind;
    double all_pc_of_z    = 100.0 * (double)all_z_size / (double)z_file->disk_so_far;

    // long form stats from --show-STATS    
    qsort (sbl, num_stats, sizeof (sbl[0]), stats_sort_by_z_size);  // sort by compressed size
    stats_output_STATS (sbl, num_stats, all_txt_size, all_uncomp_dict, all_comp_dict, all_comp_b250, all_comp_data, all_z_size, all_pc_of_txt, all_pc_of_z, all_comp_ratio);
    
    // consolidates stats of related lines into one - for SEQ, QUAL, DESC and QNAME
    stats_consolidate_compound (sbl, num_stats, "Description", "DESC");
    
    stats_consolidate_compound (sbl, num_stats, "QNAME", "QNAME");

    stats_consolidate_related (sbl, num_stats, "Sequence", 5, "SQBITMAP", "GPOS", "NONREF", "NONREF_X", "STRAND");
    
    stats_consolidate_related (sbl, num_stats, "Quality",  2, "QUAL", "DOMQRUNS");
    
    stats_consolidate_related (sbl, num_stats, "Reference", 5, ST_NAME (SEC_REFERENCE), ST_NAME (SEC_REF_IS_SET), 
                               ST_NAME (SEC_REF_CONTIGS), ST_NAME (SEC_REF_RAND_ACC), ST_NAME (SEC_REF_ALT_CHROMS));

    stats_consolidate_related (sbl, num_stats, "Other", 12, "E1L", "E2L", "EOL", "SAMPLES", "OPTIONAL", "TOPLEVEL", "LINEMETA", "CONTIG",
                               ST_NAME (SEC_RANDOM_ACCESS), ST_NAME (SEC_DICT_ID_ALIASES), 
                               ST_NAME (SEC_TXT_HEADER), ST_NAME (SEC_VB_HEADER));

    stats_consolidate_related (sbl, num_stats, "BD_BI",  3, "BD_BI", "BD:Z", "BI:Z");

    // short form stats from --show-stats    
    qsort (sbl, num_stats, sizeof (sbl[0]), stats_sort_by_z_size);  // re-sort after consolidation
    stats_output_stats (sbl, num_stats, all_txt_size, all_z_size, all_pc_of_txt, all_pc_of_z, all_comp_ratio);

    stats_check_count (all_z_size);

    // note: we use txt_data_so_far_single and not txt_data_size_single, because the latter has estimated size if disk_so_far is 
    // missing, while txt_data_so_far_single is what was actually processed
    char s1[30], s2[30];
    ASSERTW (all_txt_size == z_file->txt_data_so_far_bind || flag_optimize || flag_make_reference, 
             "Hmm... incorrect calculation for %s sizes: total section sizes=%s but file size is %s (diff=%d)", 
             dt_name (z_file->data_type), str_uint_commas (all_txt_size, s1), str_uint_commas (z_file->txt_data_so_far_bind, s2), 
             (int32_t)(z_file->txt_data_so_far_bind - all_txt_size)); 

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
