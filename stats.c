// ------------------------------------------------------------------
//   stats.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "dict_id.h"
#include "strings.h"
#include "stats.h"
#include "section_types.h"
#include "file.h"
#include "vblock.h"
#include "txtfile.h"

// VCF: used for allocating compressed SEC_VCF_GT_DATA bytes to the ctxs that contributed to it, pro-rated by the number
// of words each ctx contributed
static struct { uint64_t words_total, words_so_far, comp_bytes_so_far; } vcf_gt_data_tracker;
static int *count_per_section = NULL;

#define OVERHEAD_SEC_VB_HDR      -1
#define OVERHEAD_SEC_GENOZIP_HDR -2
#define OVERHEAD_SEC_TXT_HDR     -3
#define OVERHEAD_SEC_RA_INDEX    -4
#define LAST_OVERHEAD            -4
static void stats_get_sizes (DictIdType dict_id /* option 1 */, int overhead_sec /* option 2*/, 
                             uint64_t *dict_compressed_size, uint64_t *b250_compressed_size, uint64_t *local_compressed_size,
                             double format_proportion /* only needed for FORMAT subfields - number of words of this dict_id in gt_data divided by all words in gt_data */)
{
    if (!count_per_section) count_per_section = calloc (z_file->section_list_buf.len, sizeof (int));

    *dict_compressed_size = *b250_compressed_size = *local_compressed_size = 0;
    uint64_t total_gt_data=0;

    for (unsigned i=0; i < z_file->section_list_buf.len; i++) {

        SectionListEntry *section = ENT (SectionListEntry, z_file->section_list_buf, i);

        count_per_section[i]++; // we're optimistically assuming section will be count_per_section - we will revert if not
        
        if (section->dict_id.num == dict_id.num && section->section_type == SEC_DICT)
            *dict_compressed_size += (section+1)->offset - section->offset;

        else if (section->dict_id.num == dict_id.num && section->section_type == SEC_B250)
            *b250_compressed_size += (section+1)->offset - section->offset;

        else if (section->dict_id.num == dict_id.num && section->section_type == SEC_LOCAL)
            *local_compressed_size += (section+1)->offset - section->offset;

        else if (section->section_type == SEC_VB_HEADER && overhead_sec == OVERHEAD_SEC_VB_HDR)
            *local_compressed_size += (z_file->data_type == DT_VCF) ? sizeof (SectionHeaderVbHeaderVCF) // payload is haplotype index which goes to VCF_GT
                                                                    : (section+1)->offset - section->offset;

        else if (section->section_type == SEC_GENOZIP_HEADER && overhead_sec == OVERHEAD_SEC_GENOZIP_HDR)
            *local_compressed_size += z_file->disk_size - section->offset;

        else if (section->section_type == SEC_DICT_ID_ALIASES && overhead_sec == OVERHEAD_SEC_GENOZIP_HDR)
            *local_compressed_size += z_file->disk_size - section->offset;

        else if (section->section_type == SEC_TXT_HEADER && overhead_sec == OVERHEAD_SEC_TXT_HDR)
            *local_compressed_size += (section+1)->offset - section->offset;

        else if (section->section_type == SEC_RANDOM_ACCESS && overhead_sec == OVERHEAD_SEC_RA_INDEX)
            *local_compressed_size += (section+1)->offset - section->offset;
        
        else if ((section->section_type == SEC_VCF_HT_DATA && dict_id.num == dict_id_fields[VCF_GT]) ||
                 (section->section_type == SEC_VCF_PHASE_DATA && dict_id.num == dict_id_fields[VCF_GT]))
            *local_compressed_size += (section+1)->offset - section->offset;
            
        // we allocate the size of SEC_VCF_GT_DATA by the proportion of words due to this FORMAT dict_id
        else if (section->section_type == SEC_VCF_GT_DATA && dict_id_is_vcf_format_sf (dict_id)) {

            if (vcf_gt_data_tracker.words_so_far != vcf_gt_data_tracker.words_total) {
                uint64_t increment = (uint64_t)(format_proportion * (double)((section+1)->offset - section->offset));
                *b250_compressed_size += increment;
                vcf_gt_data_tracker.comp_bytes_so_far += increment;
            }
            else
                total_gt_data += ((section+1)->offset - section->offset);
        }

        else if (section->section_type == SEC_VB_HEADER && dict_id.num == dict_id_fields[VCF_GT]) {// squeezed index is part of GT data
            *local_compressed_size += (section+1)->offset - section->offset - sizeof (SectionHeaderVbHeaderVCF); // section header already accounted for
            count_per_section[i]--; // we already counted SEC_VB_HEADER above. this is just accounting for the squeeze payload
        }
        
        else count_per_section[i]--; // acually... not count_per_section!
    }

    // case: the last dict_id contributing to gt_data. take everything that's remaining to make sure all the components add up
    if (z_file->data_type == DT_VCF && dict_id_is_vcf_format_sf (dict_id) && vcf_gt_data_tracker.words_so_far == vcf_gt_data_tracker.words_total)
        *b250_compressed_size = total_gt_data - vcf_gt_data_tracker.comp_bytes_so_far;
}

static void stats_check_count (uint64_t all_comp_total)
{
    if (all_comp_total == z_file->disk_size) return; // all good

    ARRAY (SectionListEntry, sections, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++) 
        if (!count_per_section[i]) 
            fprintf (stderr, "Section not counted: %s section_i=%u\n", st_name (sections[i].section_type), i);
        else if (count_per_section[i] > 1 && ENT (SectionListEntry, z_file->section_list_buf, i)->section_type != SEC_VCF_GT_DATA) // note: SEC_VCF_GT_DATA is expected to have 1 count per FORMAT dict. to do: verify this.
            fprintf (stderr, "Section overcounted: %s section_i=%u counted %u times\n", st_name (sections[i].section_type), i, count_per_section[i]);

    char s1[20], s2[20];
    ASSERTW (false, "Hmm... incorrect calculation for GENOZIP sizes: total section sizes=%s but file size is %s (diff=%d)", 
             str_uint_commas (all_comp_total, s1), str_uint_commas (z_file->disk_size, s2), (int32_t)(z_file->disk_size - all_comp_total));
}

static void stats_show_file_metadata (void)
{
    fprintf (stderr, "\n\n");
    if (txt_file->name) fprintf (stderr, "%s file name: %s\n", dt_name (z_file->data_type), txt_file->name);
    
    char ls[30];
    if (z_file->data_type == DT_VCF) 
        fprintf (stderr, "Individuals: %u   ", global_vcf_num_samples);

    fprintf (stderr, "%s: %s   Dictionaries: %u   Vblocks: %u   Sections: %u\n", 
                DTPZ (show_sections_line_name), str_uint_commas (z_file->num_lines, ls), z_file->num_dict_ids, 
                z_file->num_vbs, (uint32_t)z_file->section_list_buf.len);
}

static void initialize_vcf_gt_data_tracker (void)
{
    memset (&vcf_gt_data_tracker, 0, sizeof (vcf_gt_data_tracker));

    for (uint32_t i=0; i < z_file->num_dict_ids; i++) { 

        const MtfContext *ctx = &z_file->mtf_ctx[i];
        if (dict_id_is_vcf_format_sf (ctx->dict_id))
            vcf_gt_data_tracker.words_total += ctx->mtf_i.len;
    }
}

typedef struct {
    uint64_t total_comp_size;
    char did_i[20], name[50], type[20], words[20], hash[20], uncomp_dict[20], comp_dict[20],
         comp_b250[20], comp_data[20], comp_total[20], txt[20];
         double comp_ratio, pc_txt, pc_genozip, pc_dict, pc_singletons, pc_failed_singletons, pc_hash_occupancy;
} StatsByLine;


static int stats_sort_by_total_comp_size(const void *a, const void *b)  
{ 
    return ((StatsByLine*)b)->total_comp_size - ((StatsByLine*)a)->total_comp_size;
}

void stats_show_sections (void)
{
    stats_show_file_metadata();

    uint64_t all_comp_dict=0, all_uncomp_dict=0, all_comp_b250=0, all_comp_data=0, all_comp_total=0, all_txt=0;

    //prepare data
    StatsByLine sbl[MAX_DICTS+10], *s = sbl;
    memset (sbl, 0, sizeof(sbl));

    if (z_file->data_type == DT_VCF) initialize_vcf_gt_data_tracker();

    for (int i=LAST_OVERHEAD; i < (int)z_file->num_dict_ids; i++) { 

        MtfContext *ctx = (i>=0) ? &z_file->mtf_ctx[i] : NULL;

        // for FORMAT subfields, we compress all of them together in GT_DATA, because they are often correlated
        // we allocate their compressed size as the proportion of GT_DATA words that are of this subfield
        double format_proportion = 0;
        if (z_file->data_type == DT_VCF && ctx && dict_id_is_vcf_format_sf (ctx->dict_id)) {
            format_proportion = (double)ctx->mtf_i.len / (double)vcf_gt_data_tracker.words_total;
            vcf_gt_data_tracker.words_so_far += ctx->mtf_i.len;
        }
   
        uint64_t dict_compressed_size, b250_compressed_size, local_compressed_size, txt_len;
        stats_get_sizes (ctx ? ctx->dict_id : DICT_ID_NONE, ctx ? 0 : i, 
                         &dict_compressed_size, &b250_compressed_size, &local_compressed_size, format_proportion);

        s->total_comp_size = dict_compressed_size + b250_compressed_size + local_compressed_size;
        if (ctx && !ctx->mtf_i.len && !ctx->txt_len && !ctx->b250.len && !s->total_comp_size) continue;

        txt_len = ctx ? ctx->txt_len : (i==OVERHEAD_SEC_TXT_HDR ? txtfile_get_last_header_len() : 0);

        all_comp_dict   += dict_compressed_size;
        all_uncomp_dict += ctx ? ctx->dict.len : 0;
        all_comp_b250   += b250_compressed_size;
        all_comp_data   += local_compressed_size;
        all_comp_total  += s->total_comp_size;
        all_txt         += txt_len;

        /* Name         */ 
        if (ctx)                              strcpy (s->name, ctx->name); 
        else if (i==OVERHEAD_SEC_GENOZIP_HDR) strcpy (s->name, "genozip hdr");
        else if (i==OVERHEAD_SEC_VB_HDR)      strcpy (s->name, "VBlock hdrs");
        else if (i==OVERHEAD_SEC_TXT_HDR)     strcpy (s->name, "Txt file header");
        else if (i==OVERHEAD_SEC_RA_INDEX)    strcpy (s->name, "--regions index");

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
            s->did_i[0] = s->words[0] = s->hash[0] = s->uncomp_dict[0] = s->comp_dict[0] = '-';
        
        /* txt            */ str_size (txt_len, s->txt);
        /* comp b250      */ str_size (b250_compressed_size, s->comp_b250);
        /* comp data      */ str_size (local_compressed_size, s->comp_data);
        /* comp total     */ str_size (s->total_comp_size, s->comp_total);
        /* comp ratio     */ s->comp_ratio = (double)txt_len / (double)s->total_comp_size;
        /* % of txt       */ s->pc_txt     = 100.0 * (double)txt_len / (double)txt_file->txt_data_so_far_single;
        /* % of genozip   */ s->pc_genozip = 100.0 * (double)s->total_comp_size / (double)z_file->disk_size;

        s++;
    }
    unsigned num_stats = s - sbl;

    // sort by total compressed size
    qsort (sbl, num_stats, sizeof (sbl[0]), stats_sort_by_total_comp_size);

    fprintf (stderr, "\nSections (sorted by %% of genozip file):\n");
    fprintf (stderr, "did_i Name            Type            #Words  Snips-(%% of #Words)        Hash-table    uncomp      comp      comp      comp      comp       txt    comp   %% of   %% of  \n");
    fprintf (stderr, "                                     in file   Dict  Local   Both         Size Occp      dict      dict      b250     local     TOTAL             ratio    txt    zip\n");

    for (uint32_t i=0; i < num_stats; i++) { // don't show CHROM-FORMAT as they are already showed above

        StatsByLine *s = &sbl[i];
        if (!s->total_comp_size) continue;

#define PC(pc) ((pc==0 || pc>=10) ? 0 : (pc<1 ? 2:1))
        fprintf (stderr, "%-2.2s    %-15.15s %-6.6s %15s  %4.*f%%  %4.*f%%  %4.*f%% %12s %3.0f%% %9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
                 s->did_i, s->name, s->type, s->words, 
                 PC (s->pc_dict), s->pc_dict, PC(s->pc_singletons), s->pc_singletons, PC(s->pc_failed_singletons), s->pc_failed_singletons, 
                 s->hash, s->pc_hash_occupancy, // Up to here - these don't appear in the total
                 s->uncomp_dict, s->comp_dict, s->comp_b250, s->comp_data, s->comp_total, s->txt, s->comp_ratio, s->pc_txt, s->pc_genozip);
    }

    double all_comp_ratio = (double)all_txt / (double)all_comp_total;
    double all_pc_txt     = 100.0 * (double)all_txt        / (double)txt_file->txt_data_so_far_single;
    double all_pc_genozip = 100.0 * (double)all_comp_total / (double)z_file->disk_size;
    
    char s1[20], s2[20], s3[20], s4[20], s5[20], s6[20];
    fprintf (stderr, "TOTAL                                                                               "
             "%9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
             str_size (all_uncomp_dict, s1), str_size (all_comp_dict, s2),  str_size (all_comp_b250, s3), 
             str_size (all_comp_data, s4),   str_size (all_comp_total, s5), str_size (all_txt, s6), 
             all_comp_ratio, all_pc_txt, all_pc_genozip);

    stats_check_count(all_comp_total);

    // note: we use txt_data_so_far_single and not txt_data_size_single, because the latter has estimated size if disk_size is 
    // missing, while txt_data_so_far_single is what was actually processed
    ASSERTW (all_txt == txt_file->txt_data_so_far_single || flag_optimize, "Hmm... incorrect calculation for %s sizes: total section sizes=%s but file size is %s (diff=%d)", 
             dt_name (z_file->data_type), str_uint_commas (all_txt, s1), str_uint_commas (txt_file->txt_data_size_single, s2), 
             (int32_t)(txt_file->txt_data_so_far_single - all_txt)); 
}
