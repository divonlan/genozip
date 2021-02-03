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
#include "arch.h"
#include "codec.h"

static void stats_get_sizes (DictId dict_id /* option 1 */, SectionType non_ctx_sec /* option 2*/, 
                             int64_t *dict_compressed_size, int64_t *b250_compressed_size, int64_t *local_compressed_size,
                             int *count_per_section)
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

static void stats_check_count (uint64_t all_z_size, const int *count_per_section)
{
    if (all_z_size == z_file->disk_so_far) return; // all good

    ARRAY (SectionListEntry, sections, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++) 
        if (!count_per_section[i]) 
            WARN ("stats_check_count: Section not counted: %s section_i=%u\n", st_name (sections[i].section_type), i);
        else if (count_per_section[i] > 1) 
            WARN ("stats_check_count: Section overcounted: %s section_i=%u dict=%s counted %u times\n", 
                   st_name (sections[i].section_type), i, dis_dict_id (sections[i].dict_id).s, count_per_section[i]);

    WARN ("Hmm... incorrect calculation for GENOZIP sizes: total section sizes=%s but file size is %s (diff=%d)", 
          str_uint_commas (all_z_size).s, str_uint_commas (z_file->disk_so_far).s, (int32_t)(z_file->disk_so_far - all_z_size));
}

static void stats_show_file_metadata (Buffer *buf)
{
    bufprintf (evb, buf, "%s", "\n\n");
    if (txt_file->name) 
        bufprintf (evb, buf, "%s file%s%s: %.*s\n", dt_name (z_file->data_type), 
                   z_file->bound_txt_names.param > 1 ? "s" : "", // param holds the number of txt files
                   flag.pair ? " (paired)" : "",
                   (int)z_file->bound_txt_names.len, z_file->bound_txt_names.data);
    
    if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) 
        bufprintf (evb, buf, "Reference: %s\n", ref_filename);

    if (z_file->data_type == DT_VCF) 
        bufprintf (evb, buf, "Samples: %u   ", vcf_header_get_num_samples());

    bufprintf (evb, buf, "%s: %s   Dictionaries: %u   Vblocks: %u x %u MB  Sections: %u\n", 
               DTPZ (show_stats_line_name), str_uint_commas (z_file->num_lines).s, z_file->num_contexts, 
               z_file->num_vbs, (uint32_t)(flag.vblock_memory >> 20), (uint32_t)z_file->section_list_buf.len);

    char timestr[100];
    time_t now = time (NULL);
    strftime (timestr, 100, "%Y-%m-%d %H:%M:%S", localtime (&now));
    bufprintf (evb, buf, "Genozip version: %s %s\nDate compressed: %s %s\n", 
               GENOZIP_CODE_VERSION, arch_get_distribution(), timestr, tzname[daylight]);
}

typedef struct {
    DidIType my_did_i, st_did_i;
    int64_t txt_size, z_size;
    const char *name;
    char  type[20];
    StrText did_i, words, hash, uncomp_dict, comp_dict, comp_b250, comp_data;
    double pc_of_txt, pc_of_z, pc_dict, pc_singletons, pc_failed_singletons, pc_hash_occupancy;
} StatsByLine;

static int stats_sort_by_z_size(const void *a, const void *b)  
{ 
    return (((StatsByLine*)b)->z_size > ((StatsByLine*)a)->z_size) ? 1 : -1; // use comparison (>) and not minus (-) as the retun value is only 32 bit
}

static void stats_consolidate_ctxs (StatsByLine *sbl, unsigned num_stats)
{
    for (unsigned parent=0; parent < num_stats; parent++) 

        // case: we might consolidate to this context
        if (sbl[parent].my_did_i != DID_I_NONE && sbl[parent].st_did_i == DID_I_NONE)  {
            for (unsigned child=0; child < num_stats; child++) {

                if (sbl[parent].my_did_i  == sbl[child].st_did_i) {
                    sbl[parent].txt_size  += sbl[child].txt_size;
                    sbl[parent].z_size    += sbl[child].z_size;  
                    sbl[parent].pc_of_txt += sbl[child].pc_of_txt;
                    sbl[parent].pc_of_z   += sbl[child].pc_of_z;  

                    sbl[child] = (StatsByLine){ .my_did_i = DID_I_NONE, .st_did_i = DID_I_NONE };
                }
            }

            if (!strcmp (sbl[parent].name, "SQBITMAP")) sbl[parent].name = "SEQ"; // rename
        }
}

void stats_set_consolidation (VBlock *vb, DidIType parent, unsigned num_deps, ...)
{
    va_list args;
    va_start (args, num_deps);

    for (unsigned d=0; d < num_deps; d++) {
        DidIType dep = (DidIType)va_arg (args, int); 
        vb->contexts[dep].st_did_i = parent;
    }

    vb->contexts[parent].is_stats_parent = true;
    
    va_end (args);
}

static void stats_consolidate_non_ctx (StatsByLine *sbl, unsigned num_stats, const char *consolidated_name, unsigned num_deps, ...)
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

static void stats_output_stats (StatsByLine *s, unsigned num_stats, double txt_ratio, Codec codec,
                                int64_t all_txt_size, int64_t all_z_size, double all_pc_of_txt, double all_pc_of_z, double all_comp_ratio)
{
    bufprintf (evb, &z_file->stats_buf, "\nSections (sorted by %% of genozip file):%s\n", "");
    bufprintf (evb, &z_file->stats_buf, "NAME              GENOZIP      %%      TXT       %%   RATIO\n%s", "");

    for (uint32_t i=0; i < num_stats; i++, s++)
        if (s->z_size)
            bufprintf (evb, &z_file->stats_buf, "%-15.15s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                       s->name, 
                       str_size (s->z_size).s, s->pc_of_z, // z size and % of total z that is in this line
                       str_size ((double)s->txt_size).s, s->pc_of_txt, // txt size and % of total txt which is in this line
                       (double)s->txt_size / (double)s->z_size); // ratio z vs txt
                       
    if (txt_ratio > 1)
        bufprintf (evb, &z_file->stats_buf, 
                   "TOTAL vs %-6s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                   codec_name (codec),
                   str_size (all_z_size).s, all_pc_of_z, // total z size and sum of all % of z (should be 10-=0)
                   str_size (all_txt_size / txt_ratio).s, all_pc_of_txt, // total txt fize and ratio z vs txt
                   all_comp_ratio / txt_ratio);
    
    bufprintf (evb, &z_file->stats_buf, 
               "%-15s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
               txt_ratio > 1 ? "TOTAL vs TXT" : "TOTAL",
               str_size (all_z_size).s, all_pc_of_z, // total z size and sum of all % of z (should be 10-=0)
               str_size (all_txt_size).s, all_pc_of_txt, // total txt fize and ratio z vs txt
               all_comp_ratio);
}

static void stats_output_STATS (StatsByLine *s, unsigned num_stats,
                                int64_t all_txt_size, int64_t all_uncomp_dict, int64_t all_comp_dict, int64_t all_comp_b250, int64_t all_comp_data, 
                                int64_t all_z_size, double all_pc_of_txt, double all_pc_of_z, double all_comp_ratio)
{
#define PC(pc) ((pc==0 || pc>=10) ? 0 : (pc<1 ? 2:1))

    // add diagnostic info
    bufprintf (evb, &z_file->STATS_buf, "Command line: %s", "");
    buf_add_string (evb, &z_file->STATS_buf, flags_command_line()->data); // careful not to use bufprintf with command_line as it can exceed the maximum length in bufprintf
    bufprintf (evb, &z_file->STATS_buf, "\nSystem info: OS=%s cores=%u endianity=%s\n", 
               arch_get_os(), arch_get_num_cores(), arch_get_endianity());
    bufprintf (evb, &z_file->STATS_buf, "\nSections (sorted by %% of genozip file):%s\n", "");
    bufprintf (evb, &z_file->STATS_buf, "did_i Name            Type            #Words  Snips-(%% of #Words)        Hash-table    uncomp      comp      comp      comp      comp       txt    comp   %% of   %% of  %s\n", "");
    bufprintf (evb, &z_file->STATS_buf, "                                     in file   Dict  Local   Both         Size Occp      dict      dict      b250     local     TOTAL             ratio    txt    zip%s\n", "");

    for (uint32_t i=0; i < num_stats; i++, s++)
        if (s->z_size)
            bufprintf (evb, &z_file->STATS_buf, "%-2.2s    %-15.15s %-6.6s %15s  %4.*f%%  %4.*f%%  %4.*f%% %12s %3.0f%% %9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
                       s->did_i.s, s->name, s->type, s->words.s, 
                       PC (s->pc_dict), s->pc_dict, PC(s->pc_singletons), s->pc_singletons, PC(s->pc_failed_singletons), s->pc_failed_singletons, 
                       s->hash.s, s->pc_hash_occupancy, // Up to here - these don't appear in the total
                       s->uncomp_dict.s, s->comp_dict.s, s->comp_b250.s, s->comp_data.s, str_size (s->z_size).s, 
                       str_size ((double)s->txt_size).s, (double)s->txt_size / (double)s->z_size, s->pc_of_txt, s->pc_of_z);

    bufprintf (evb, &z_file->STATS_buf, "TOTAL                                                                               "
               "%9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
               str_size (all_uncomp_dict).s, str_size (all_comp_dict).s,  str_size (all_comp_b250).s, 
               str_size (all_comp_data).s,   str_size (all_z_size).s, str_size (all_txt_size).s, 
               all_comp_ratio, all_pc_of_txt, all_pc_of_z);
}

// generate the stats text - all sections except genozip header and the two stats sections 
void stats_compress (void)
{
    stats_show_file_metadata(&z_file->stats_buf);
    buf_copy (evb, &z_file->STATS_buf, &z_file->stats_buf, 0,0,0, "z_file->STATS_buf");

    int64_t all_comp_dict=0, all_uncomp_dict=0, all_comp_b250=0, all_comp_data=0, all_z_size=0, all_txt_size=0;

    // prepare data
    StatsByLine sbl[MAX_DICTS + NUM_SEC_TYPES] = { }, *s = sbl;

    static Buffer count_per_section_buf = EMPTY_BUFFER;
    buf_alloc (evb, &count_per_section_buf, z_file->section_list_buf.len * sizeof (int), 1, "count_per_section");
    buf_zero (&count_per_section_buf);
    ARRAY (int, count_per_section, count_per_section_buf);

    #define SEC(i) (i<0 ? -(i)-1 : SEC_NONE) // i to section type
    #define ST_NAME(st) (&st_name(st)[4]) // cut off "SEC_" 

    for (int i=-NUM_SEC_TYPES; i < (int)z_file->num_contexts; i++) { // sections go into -1 to -NUM_SEC_TYPES (see SEC())

        Context *ctx = (i>=0) ? &z_file->contexts[i] : NULL;

        if (SEC(i) == SEC_DICT || SEC(i) == SEC_B250 || SEC(i) == SEC_LOCAL) continue; // these are covered by indiviual contexts

        int64_t dict_compressed_size, b250_compressed_size, local_compressed_size;
        stats_get_sizes (ctx ? ctx->dict_id : DICT_ID_NONE, SEC(i), 
                         &dict_compressed_size, &b250_compressed_size, &local_compressed_size, count_per_section);

        s->z_size = dict_compressed_size + b250_compressed_size + local_compressed_size;

        ASSERTW (s->z_size >= 0, "Hmm... s->z_size=%"PRId64" is negative for %s", s->z_size, s->name);

        if (ctx && !ctx->b250.num_b250_words && !ctx->txt_len && !ctx->b250.len && !ctx->is_stats_parent && !s->z_size) 
            continue;

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
                             s->my_did_i = ctx->did_i;
                             s->st_did_i = ctx->st_did_i;
        /* did_i          */ s->did_i = str_uint_commas ((uint64_t)ctx->did_i); 
        /* #Words in file */ s->words = str_uint_commas (ctx->b250.num_b250_words);
        /* % dict         */ s->pc_dict              = !ctx->b250.num_b250_words         ? 0 : 100.0 * (double)ctx->nodes.len / (double)ctx->b250.num_b250_words;
        /* % singletons   */ s->pc_singletons        = !ctx->b250.num_b250_words         ? 0 : 100.0 * (double)ctx->num_singletons / (double)ctx->b250.num_b250_words;
        /* % failed singl.*/ s->pc_failed_singletons = !ctx->b250.num_b250_words         ? 0 : 100.0 * (double)ctx->num_failed_singletons / (double)ctx->b250.num_b250_words;
        /* % hash occupn. */ s->pc_hash_occupancy    = !ctx->global_hash_prime  ? 0 : 100.0 * (double)(ctx->nodes.len + ctx->ol_nodes.len) / (double)ctx->global_hash_prime;
        /* Hash           */ s->hash = str_uint_commas (ctx->global_hash_prime);
        /* uncomp dict    */ s->uncomp_dict = str_size (ctx->dict.len);
        /* comp dict      */ s->comp_dict   = str_size (dict_compressed_size);
        }
        else {
            s->my_did_i = s->st_did_i = DID_I_NONE;
            s->did_i.s[0] = s->words.s[0] = s->hash.s[0] = s->uncomp_dict.s[0] = '-';
        }

        if (ctx)
        /* comp dict      */ s->comp_dict = str_size (dict_compressed_size);
        else 
            s->comp_dict.s[0] = '-';

        /* comp b250      */ s->comp_b250 = str_size (b250_compressed_size);
        /* comp data      */ s->comp_data = str_size (local_compressed_size);
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
    stats_output_STATS (sbl, num_stats, 
                        all_txt_size, all_uncomp_dict, all_comp_dict, all_comp_b250, all_comp_data, all_z_size, all_pc_of_txt, all_pc_of_z, all_comp_ratio);
    
    // consolidates stats of child contexts into the parent one
    stats_consolidate_ctxs (sbl, num_stats);
    
    stats_consolidate_non_ctx (sbl, num_stats, "Reference", 5, ST_NAME (SEC_REFERENCE), ST_NAME (SEC_REF_IS_SET), 
                               ST_NAME (SEC_REF_CONTIGS), ST_NAME (SEC_REF_RAND_ACC), ST_NAME (SEC_REF_ALT_CHROMS));

    stats_consolidate_non_ctx (sbl, num_stats, "Other", 16, "E1L", "E2L", "EOL", "SAMPLES", "OPTIONAL", 
                               TOPLEVEL, "TOP2BAM", "TOP2FQ", "TOP2VCF", "LINEMETA", "CONTIG",
                               ST_NAME (SEC_RANDOM_ACCESS), ST_NAME (SEC_DICT_ID_ALIASES), 
                               ST_NAME (SEC_TXT_HEADER), ST_NAME (SEC_VB_HEADER), ST_NAME (SEC_BGZF));

    // consolidate SAM arrays

    // short form stats from --show-stats    
    qsort (sbl, num_stats, sizeof (sbl[0]), stats_sort_by_z_size);  // re-sort after consolidation

    double txt_ratio = (double)z_file->txt_data_so_far_bind / (double)z_file->txt_disk_so_far_bind; // source compression, eg BGZF

    stats_output_stats (sbl, num_stats, txt_ratio, txt_file->codec, all_txt_size, all_z_size, all_pc_of_txt, all_pc_of_z, all_comp_ratio);

    stats_check_count (all_z_size, count_per_section);

    // note: we use txt_data_so_far_single and not txt_data_size_single, because the latter has estimated size if disk_so_far is 
    // missing, while txt_data_so_far_single is what was actually processed
    ASSERTW (all_txt_size == z_file->txt_data_so_far_bind || flag.optimize || flag.make_reference, 
             "Hmm... incorrect calculation for %s sizes: total section sizes=%s but file size is %s (diff=%d)", 
             dt_name (z_file->data_type), str_uint_commas (all_txt_size).s, str_uint_commas (z_file->txt_data_so_far_bind).s, 
             (int32_t)(z_file->txt_data_so_far_bind - all_txt_size)); 

    zfile_compress_section_data (evb, SEC_STATS, &z_file->stats_buf);
    zfile_compress_section_data (evb, SEC_STATS, &z_file->STATS_buf);

    buf_free (&count_per_section_buf);
}

void stats_display (void)
{
    Buffer *buf = flag.show_stats == 1 ? &z_file->stats_buf : &z_file->STATS_buf;

    if (!buf_is_allocated (buf)) return; // no stats available

    buf_print (buf , false);

    const SectionListEntry *sl = sections_get_first_section_of_type (SEC_STATS, false);

    if (z_file->disk_size < (1<<20))  // no need to print this note if total size > 1MB, as the ~2K of overhead is rounded off anyway
        // stats text doesn't include SEC_STATS and SEC_GENOZIP_HEADER - the last 2 sections in the file - since stats text is generated before these sections are compressed
        iprintf ("\nNote: ZIP total file size excludes overhead of %s\n", str_size (z_file->disk_size - sl->offset).s);

    iprint0 ("\n");
}

void stats_read_and_display (void)
{
    const SectionListEntry *sl = sections_get_first_section_of_type (SEC_STATS, true);
    if (!sl) return; // genozip file does not contain stats sections (SEC_STATS was introduced in v 7.0.5)

    // read and uncompress the requested stats section
    zfile_read_section (z_file, evb, 0, &evb->z_data, "z_data", SEC_STATS, sl + (flag.show_stats==2));
    zfile_uncompress_section (evb, evb->z_data.data, flag.show_stats == 1 ? &z_file->stats_buf : &z_file->STATS_buf, "z_file->stats_buf", 0, SEC_STATS);
    buf_free (&evb->z_data);
    
    stats_display();

    if (exe_type == EXE_GENOCAT) exit_ok; // if this is genocat - we're done
}

// concatenate txt names of bound files so we can show them all
void stats_add_txt_name (const char *fn)
{
    bufprintf (evb, &z_file->bound_txt_names, "%s%s", z_file->bound_txt_names.len ? " ": "", fn);
    z_file->bound_txt_names.param++; // we store the number of files in param
}
