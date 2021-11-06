// ------------------------------------------------------------------
//   stats.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <stdarg.h>
#include "genozip.h"
#include "buffer.h"
#include "strings.h"
#include "stats.h"
#include "sections.h"
#include "file.h"
#include "vblock.h"
#include "txtheader.h"
#include "reference.h"
#include "zfile.h"
#include "version.h"
#include "arch.h"
#include "codec.h"
#include "license.h"
#include "segconf.h"
#include "qname.h"

static void stats_get_sizes (DictId dict_id /* option 1 */, SectionType non_ctx_sec /* option 2*/, 
                             int64_t *dict_compressed_size, int64_t *b250_compressed_size, int64_t *local_compressed_size,
                             int *count_per_section)
{
    *dict_compressed_size = *b250_compressed_size = *local_compressed_size = 0;

    unsigned i=0;
    for (Section sec = section_next(0); sec; sec = section_next (sec), i++) {

        int64_t after_sec = section_next (sec) ? (sec+1)->offset : z_file->disk_so_far;
        int64_t sec_size = after_sec - sec->offset;

        count_per_section[i]++; // we're optimistically assuming sec will be count_per_section - we will revert if not

        if (sec->dict_id.num == dict_id.num && (sec->st == SEC_DICT || sec->st == SEC_COUNTS))
            *dict_compressed_size += sec_size;

        else if (sec->dict_id.num == dict_id.num && sec->st == SEC_B250)
            *b250_compressed_size += sec_size;

        else if ((sec->dict_id.num == dict_id.num && sec->st == SEC_LOCAL) ||
                 (sec->st == non_ctx_sec))
            *local_compressed_size += sec_size;
             
        else count_per_section[i]--; // acually... not count_per_section!
    }
}

static void stats_check_count (uint64_t all_z_size, const int *count_per_section)
{
    if (all_z_size == z_file->disk_so_far) return; // all good

    ARRAY (SectionEnt, sections, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++) 
        if (!count_per_section[i]) 
            WARN ("stats_check_count: Section not counted: %s section_i=%u dict=%s\n", 
                  st_name (sections[i].st), i, dis_dict_id (sections[i].dict_id).s);
        else if (count_per_section[i] > 1) 
            WARN ("stats_check_count: Section overcounted: %s section_i=%u dict=%s counted %u times\n", 
                   st_name (sections[i].st), i, dis_dict_id (sections[i].dict_id).s, count_per_section[i]);

    WARN ("Hmm... incorrect calculation for GENOZIP sizes: total section sizes=%s but file size is %s (diff=%d)", 
          str_uint_commas (all_z_size).s, str_uint_commas (z_file->disk_so_far).s, (int32_t)(z_file->disk_so_far - all_z_size));
}

static void stats_output_file_metadata (Buffer *buf)
{
    bufprintf (evb, buf, "%s", "\n\n");
    if (txt_file->name) 
        bufprintf (evb, buf, "%s file%s%s: %.*s\n", dt_name (z_file->data_type), 
                   z_file->bound_txt_names.param > 1 ? "s" : "", // param holds the number of txt files
                   flag.pair ? " (paired)" : "",
                   (int)z_file->bound_txt_names.len, z_file->bound_txt_names.data);
    
    if (flag.reference == REF_MAKE_CHAIN) {
        bufprintf (evb, buf, "PRIM reference: %s\n", ref_get_filename (prim_ref));
        bufprintf (evb, buf, "LUFT reference: %s\n", ref_get_filename (gref));
    }

    else if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE || flag.reference == REF_LIFTOVER) 
        bufprintf (evb, buf, "Reference: %s\n", ref_get_filename (gref));

    if (Z_DT(DT_VCF)) 
        bufprintf (evb, buf, "Samples: %u   ", vcf_header_get_num_samples());

    bufprintf (evb, buf, "%ss: %s   Dictionaries: %u   Vblocks: %u x %u MB  Sections: %u\n", 
               DTPZ (line_name), str_uint_commas (z_file->num_lines).s, z_file->num_contexts, 
               z_file->num_vbs, (uint32_t)(segconf.vb_size >> 20), (uint32_t)z_file->section_list_buf.len);

    if (Z_DT(DT_KRAKEN)) {
        int64_t dominant_taxid_count;
        const char *dominant_taxid = ctx_get_snip_with_largest_count (KRAKEN_TAXID, &dominant_taxid_count);

        if (dominant_taxid_count != -1)
            bufprintf (evb, buf, "Dominant TaxID: %s  %s: %s (%-5.2f%%)\n", dominant_taxid, DTPZ (line_name),
                       str_uint_commas (dominant_taxid_count).s, 100.0 * (float)dominant_taxid_count / (float)z_file->num_lines); 
        else
            bufprint0 (evb, buf, "Dominant TaxID: No dominant species\n"); 
    }  
    
    else if (kraken_is_loaded) 
        bufprintf (evb, buf, "Features: Per-line taxonomy ID data\n%s", "");

    else if (Z_DT(DT_CHAIN) && flag.reference == REF_MAKE_CHAIN && !segconf.chain_mismatches_ref)
        bufprintf (evb, buf, "Features: Chain file suitable for use with genozip --chain\n%s", "");

    if (chain_is_loaded || txt_file->coords) 
        bufprintf (evb, buf, "Features: Dual-coordinates\n%s", "");

    if ((Z_DT(DT_SAM) || Z_DT(DT_BAM)) && segconf.sam_is_sorted)
        bufprintf (evb, buf, "Sorting: Sorted by POS\n%s", "");

    if ((Z_DT(DT_SAM) || Z_DT(DT_BAM)) && segconf.sam_is_collated)
        bufprintf (evb, buf, "Sorting: Collated by QNAME\n%s", "");

    if (Z_DT(DT_FASTA))
        bufprintf (evb, buf, "Sequence type: %s\n", segconf.seq_type==SQT_AMINO ? "Amino acids" : "Nucleotide bases");

    if (segconf.qname_flavor) 
        bufprintf (evb, buf, "Read name style: %s\n", qf_name(segconf.qname_flavor));

    bufprintf (evb, buf, "Genozip version: %s %s\nDate compressed: %s\n", 
               GENOZIP_CODE_VERSION, DISTRIBUTION, str_time().s);

    if (license_has_details())
        bufprintf (evb, buf, "%s\n", license_get_one_line());
}

typedef struct {
    DidIType my_did_i, st_did_i;
    int64_t txt_len, z_size;
    char name[100];
    const char *type;
    StrText did_i, words, hash, uncomp_dict, comp_dict, comp_b250, comp_data;
    float pc_of_txt, pc_of_z, pc_dict, pc_singletons, pc_failed_singletons, pc_hash_occupancy;
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
                    sbl[parent].txt_len   += sbl[child].txt_len;
                    sbl[parent].z_size    += sbl[child].z_size;  
                    sbl[parent].pc_of_txt += sbl[child].pc_of_txt;
                    sbl[parent].pc_of_z   += sbl[child].pc_of_z;  

                    if (flag.debug_stats)
                        iprintf ("Consolidated %s (did=%u) (txt_len=%"PRIu64" z_size=%"PRIu64") into %s (did=%u) (AFTER: txt_len=%"PRIu64" z_size=%"PRIu64")\n",
                                 sbl[child].name, sbl[child].my_did_i, sbl[child].txt_len, sbl[child].z_size, 
                                 sbl[parent].name, sbl[parent].my_did_i, sbl[parent].txt_len, sbl[parent].z_size);

                    sbl[child] = (StatsByLine){ .my_did_i = DID_I_NONE, .st_did_i = DID_I_NONE };
                }
            }

            if (!strcmp (sbl[parent].name, "SQBITMAP")) strcpy (sbl[parent].name, "SEQ"); // rename
        }
}

void stats_set_consolidation (VBlock *vb, DidIType parent, unsigned num_deps, ...)
{
    va_list args;
    va_start (args, num_deps);

    for (unsigned d=0; d < num_deps; d++) {
        DidIType dep = (DidIType)va_arg (args, int); 
        CTX(dep)->st_did_i = parent;
    }

    CTX(parent)->is_stats_parent = true;
    
    va_end (args);
}

static void stats_consolidate_non_ctx (StatsByLine *sbl, unsigned num_stats, const char *consolidated_name, 
                                       unsigned num_deps, ...)
{
    va_list args;
    va_start (args, num_deps);

    const char *deps[num_deps];
    for (unsigned d=0; d < num_deps; d++) 
        deps[d] = va_arg (args, const char *);

    // use existing SBL if it matches the consolidated name
    StatsByLine *survivor = NULL;
    for (unsigned i=0; i < num_stats; i++) 
        if (sbl[i].name[0] && !strcmp (consolidated_name, sbl[i].name)) {
            survivor = &sbl[i];
            break;
        }

    for (unsigned i=0; i < num_stats; i++) {
        
        if (!sbl[i].name[0]) continue; // unused entry

        for (unsigned d=0; d < num_deps; d++)
            if (!strcmp (deps[d], sbl[i].name)) {
                if (!survivor) {
                    survivor = &sbl[i]; // first found gets to be the survivor
                }
                else {
                    survivor->txt_len   += sbl[i].txt_len;
                    survivor->z_size    += sbl[i].z_size;  
                    survivor->pc_of_txt += sbl[i].pc_of_txt;
                    survivor->pc_of_z   += sbl[i].pc_of_z; 
                    strcpy (survivor->name, consolidated_name); // rename only if at least one was consolidated

                    if (flag.debug_stats)
                        iprintf ("Consolidated %s (txt_len=%"PRIu64" z_size=%"PRIu64") into %s (AFTER: txt_len=%"PRIu64" z_size=%"PRIu64")\n",
                                 sbl[i].name, sbl[i].txt_len, sbl[i].z_size, survivor->name, survivor->txt_len, survivor->z_size);

                    sbl[i] = (StatsByLine){}; 
                }
                break;
            } 
    }

    va_end (args);
}

static void stats_output_stats (StatsByLine *s, unsigned num_stats, float src_comp_ratio, Codec codec,
                                int64_t all_txt_len, int64_t all_txt_len_0, int64_t all_z_size, float all_pc_of_txt, float all_pc_of_z, float all_comp_ratio)
{
    bufprintf (evb, &z_file->stats_buf, "\nSections (sorted by %% of genozip file):%s\n", "");
    bufprintf (evb, &z_file->stats_buf, "NAME                   GENOZIP      %%      TXT       %%   RATIO\n%s", "");

    for (uint32_t i=0; i < num_stats; i++, s++)
        if (s->z_size)
            bufprintf (evb, &z_file->stats_buf, "%-20.20s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                       s->name, 
                       str_size (s->z_size).s, s->pc_of_z, // z size and % of total z that is in this line
                       str_size ((float)s->txt_len).s, s->pc_of_txt, // txt size and % of total txt which is in this line
                       (float)s->txt_len / (float)s->z_size); // ratio z vs txt

    if (src_comp_ratio != 1)
        bufprintf (evb, &z_file->stats_buf, 
                   "GENOZIP vs %-9s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                   codec_name (codec),
                   str_size (all_z_size).s, all_pc_of_z, // total z size and sum of all % of z (should be 100)
                   str_size (all_txt_len_0 / src_comp_ratio).s, all_pc_of_txt, // total txt fize and ratio z vs txt
                   all_comp_ratio / src_comp_ratio);
    
    bufprintf (evb, &z_file->stats_buf, 
               "%-20s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
               src_comp_ratio != 1 ? "GENOZIP vs TXT" : "TOTAL",
               str_size (all_z_size).s, all_pc_of_z, // total z size and sum of all % of z (should be 100)
               str_size (all_txt_len_0).s, all_pc_of_txt, // total txt fize and ratio z vs txt
               all_comp_ratio);
}

static void stats_output_STATS (StatsByLine *s, unsigned num_stats,
                                int64_t all_txt_len, int64_t all_uncomp_dict, int64_t all_comp_dict, int64_t all_comp_b250, int64_t all_comp_data, 
                                int64_t all_z_size, float all_pc_of_txt, float all_pc_of_z, float all_comp_ratio)
{
#define PC(pc) ((pc==0 || pc>=10) ? 0 : (pc<1 ? 2:1))

    // add diagnostic info
    bufprintf (evb, &z_file->STATS_buf, "Command line: %s", "");
    buf_add_string (evb, &z_file->STATS_buf, flags_command_line()); // careful not to use bufprintf with command_line as it can exceed the maximum length in bufprintf
    bufprintf (evb, &z_file->STATS_buf, "\nSystem info: OS=%s cores=%u endianity=%s\n", 
               arch_get_os(), arch_get_num_cores(), arch_get_endianity());
    bufprintf (evb, &z_file->STATS_buf, "\nSections (sorted by %% of genozip file):%s\n", "");
    bufprintf (evb, &z_file->STATS_buf, "did_i Name              Parent            #Words  Snips-(%% of #Words)    Hash-table    uncomp      comp      comp      comp      comp       txt    comp   %% of   %% of  %s\n", "");
    bufprintf (evb, &z_file->STATS_buf, "                                         in file   Dict  Local   Both     Size Occp      dict      dict      b250     local     TOTAL             ratio    txt    zip%s\n", "");

    for (uint32_t i=0; i < num_stats; i++, s++)
        if (s->z_size)
            bufprintf (evb, &z_file->STATS_buf, "%-2.2s    %-17.17s %-17.17s %6s  %4.*f%%  %4.*f%%  %4.*f%% %8s %3.0f%% %9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
                       s->did_i.s, s->name, s->type, s->words.s, 
                       PC (s->pc_dict), s->pc_dict, PC(s->pc_singletons), s->pc_singletons, PC(s->pc_failed_singletons), s->pc_failed_singletons, 
                       s->hash.s, s->pc_hash_occupancy, // Up to here - these don't appear in the total
                       s->uncomp_dict.s, s->comp_dict.s, s->comp_b250.s, s->comp_data.s, str_size (s->z_size).s, 
                       str_size ((float)s->txt_len).s, (float)s->txt_len / (float)s->z_size, s->pc_of_txt, s->pc_of_z);

    bufprintf (evb, &z_file->STATS_buf, "TOTAL                                                                               "
               "%9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
               str_size (all_uncomp_dict).s, str_size (all_comp_dict).s,  str_size (all_comp_b250).s, 
               str_size (all_comp_data).s,   str_size (all_z_size).s, str_size (all_txt_len).s, 
               all_comp_ratio, all_pc_of_txt, all_pc_of_z);
}

// generate the stats text - all sections except genozip header and the two stats sections 
void stats_compress (void)
{    
    stats_output_file_metadata(&z_file->stats_buf);
    buf_copy (evb, &z_file->STATS_buf, &z_file->stats_buf, char,0,0, "z_file->STATS_buf");

    int64_t all_comp_dict=0, all_uncomp_dict=0, all_comp_b250=0, all_comp_data=0, all_z_size=0, all_txt_len=0;

    // prepare data
    #define NUM_SBL (NUM_SEC_TYPES + z_file->num_contexts + 2) // 2 for consolidated groups
    StatsByLine sbl[NUM_SBL], *s = sbl;
    memset (sbl, 0, NUM_SBL * sizeof (StatsByLine)); // initialize

    static Buffer count_per_section_buf = EMPTY_BUFFER;
    buf_alloc (evb, &count_per_section_buf, 0, z_file->section_list_buf.len, int, 1, "count_per_section");
    buf_zero (&count_per_section_buf);
    ARRAY (int, count_per_section, count_per_section_buf);

    #define SEC(i) (i<0 ? -(i)-1 : SEC_NONE) // i to section type
    #define ST_NAME(st) (&st_name(st)[4]) // cut off "SEC_" 

    for (int i=-NUM_SEC_TYPES; i < (int)z_file->num_contexts; i++) { // sections go into -1 to -NUM_SEC_TYPES (see SEC())

        Context *ctx = (i>=0) ? ZCTX(i) : NULL;

        if (SEC(i) == SEC_DICT || SEC(i) == SEC_B250 || SEC(i) == SEC_LOCAL || SEC(i) == SEC_COUNTS) continue; // these are covered by individual contexts

        int64_t dict_compressed_size, b250_compressed_size, local_compressed_size;
        stats_get_sizes (ctx ? ctx->dict_id : DICT_ID_NONE, SEC(i), 
                         &dict_compressed_size, &b250_compressed_size, &local_compressed_size, count_per_section);

        s->z_size = dict_compressed_size + b250_compressed_size + local_compressed_size;

        ASSERTW (s->z_size >= 0, "Hmm... s->z_size=%"PRId64" is negative for %s", s->z_size, s->name);

        if (ctx && !ctx->b250.num_ctx_words && !ctx->txt_len && !ctx->b250.len && !ctx->is_stats_parent && !s->z_size) 
            continue;

        s->txt_len = SEC(i) == SEC_TXT_HEADER  ? txtheader_get_bound_headers_len()
                   : !ctx                      ? 0 
                   :                             ctx->txt_len;
        
        all_comp_dict   += dict_compressed_size;
        all_uncomp_dict += ctx ? ctx->dict.len : 0;
        all_comp_b250   += b250_compressed_size;
        all_comp_data   += local_compressed_size;
        all_z_size      += s->z_size;
        all_txt_len     += s->txt_len;

        if (ctx) {
            if (Z_DT(DT_VCF) && dict_id_type(ctx->dict_id))
                sprintf (s->name, "%s/%s", dtype_name_z(ctx->dict_id), ctx->tag_name);
            else 
                strcpy (s->name, ctx->tag_name);

            s->type  = ctx->st_did_i != DID_I_NONE ? ZCTX(ctx->st_did_i)->tag_name : ctx->tag_name;
        }
        else {
            strcpy (s->name, ST_NAME (SEC(i))); 
            s->type  = "Global";
        }

        uint64_t n_words = !ctx                     ? 0
                         : ctx->b250.num_ctx_words  ? ctx->b250.num_ctx_words
                         : ctx->ltype != LT_TEXT    ? ctx->local.num_ctx_words
                         :                            0; // no data in context

        if (ctx) {
            s->my_did_i             = ctx->did_i;
            s->st_did_i             = ctx->st_did_i;
            s->did_i                = str_uint_commas ((uint64_t)ctx->did_i); 
            s->words                = str_uint_commas_limit (n_words, 99999);
            s->pc_dict              = !ctx->b250.num_ctx_words ? 0 : 100.0 * (float)ctx->nodes.len / (float)ctx->b250.num_ctx_words;
            s->pc_singletons        = !ctx->b250.num_ctx_words ? 0 : 100.0 * (float)ctx->num_singletons / (float)ctx->b250.num_ctx_words;
            s->pc_failed_singletons = !ctx->b250.num_ctx_words ? 0 : 100.0 * (float)ctx->num_failed_singletons / (float)ctx->b250.num_ctx_words;
            s->pc_hash_occupancy    = !ctx->global_hash_prime   ? 0 : 100.0 * (float)(ctx->nodes.len + ctx->ol_nodes.len) / (float)ctx->global_hash_prime;
            s->hash                 = str_size (ctx->global_hash_prime);
            s->uncomp_dict          = str_size (ctx->dict.len);
            s->comp_dict            = str_size (dict_compressed_size);
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
        /* % of txt       */ s->pc_of_txt = z_file->txt_data_so_far_bind ? 100.0 * (float)s->txt_len / (float)z_file->txt_data_so_far_bind : 0;
        /* % of genozip   */ s->pc_of_z   = z_file->disk_so_far          ? 100.0 * (float)s->z_size / (float)z_file->disk_so_far           : 0;

        s++;
    }
    unsigned num_stats = s - sbl;

    // note: for txt size and compression ratio in the TOTAL line (all_comp_ratio) we use txt_data_so_far_bind_0 
    // (the original txt data size) and not all_txt_size (size after ZIP modifications like --optimize). 
    // Therefore, in case of ZIP-modified txt, the sum of the (modified) fields in the TXT column will NOT equal the
    // TOTAL in the TXT column. That's ok.
    float all_comp_ratio = (float)z_file->txt_data_so_far_bind_0 /* without modifications */ / (float)all_z_size;
    float all_pc_of_txt  = z_file->txt_data_so_far_bind ? 100.0 * (float)all_txt_len / (float)z_file->txt_data_so_far_bind : 0 /* with modifications */;
    float all_pc_of_z    = z_file->disk_so_far          ? 100.0 * (float)all_z_size  / (float)z_file->disk_so_far          : 0;

    // long form stats from --STATS    
    qsort (sbl, num_stats, sizeof (sbl[0]), stats_sort_by_z_size);  // sort by compressed size

    stats_output_STATS (sbl, num_stats, 
                        all_txt_len, all_uncomp_dict, all_comp_dict, all_comp_b250, all_comp_data, all_z_size, all_pc_of_txt, all_pc_of_z, all_comp_ratio);

    // consolidates stats of child contexts into the parent one
    stats_consolidate_ctxs (sbl, num_stats);
    
    stats_consolidate_non_ctx (sbl, num_stats, 
                               flag.reference == REF_INTERNAL ? "SEQ" : "Reference", // when compressing SAM with REF_INTERNAL, count the internal reference data as part of SEQ
                               6, ST_NAME (SEC_REFERENCE), ST_NAME (SEC_REF_IS_SET), 
                               ST_NAME (SEC_REF_CONTIGS), ST_NAME (SEC_REF_RAND_ACC), ST_NAME (SEC_CHROM2REF_MAP),
                               ST_NAME (SEC_REF_IUPACS));

    stats_consolidate_non_ctx (sbl, num_stats, "Other", 18, "E1L", "E2L", "EOL", "SAMPLES", "OPTIONAL", 
                               TOPLEVEL, "ToPLUFT", "TOP2BAM", "TOP2FQ", "TOP2FQEX", "TOP2VCF", "TOP2HASH", "LINEMETA", "CONTIG", 
                               ST_NAME (SEC_RANDOM_ACCESS), ST_NAME (SEC_DICT_ID_ALIASES), 
                               ST_NAME (SEC_VB_HEADER), ST_NAME (SEC_BGZF));
    
    ASSERTW (all_txt_len == z_file->txt_data_so_far_bind || flag.make_reference, // all_txt_len=0 in make-ref as there are no contexts
             "Expecting all_txt_len=%"PRId64" == z_file->txt_data_so_far_bind=%"PRId64, all_txt_len, z_file->txt_data_so_far_bind);

    // short form stats from --stats    
    qsort (sbl, num_stats, sizeof (sbl[0]), stats_sort_by_z_size);  // re-sort after consolidation

    // source compression, eg BGZF, against txt before any modifications
    float src_comp_ratio = (float)z_file->txt_data_so_far_bind_0 / (float)z_file->txt_disk_so_far_bind; 

    stats_output_stats (sbl, num_stats, src_comp_ratio, z_file->codec, all_txt_len, z_file->txt_data_so_far_bind_0, all_z_size, all_pc_of_txt, all_pc_of_z, all_comp_ratio);

    stats_check_count (all_z_size, count_per_section);

    // note: we use txt_data_so_far_bind is the sum of recon_sizes - see zip_update_txt_counters - which is
    // expected to be the sum of txt_len. However, this NOT the size of the original file which is stored in
    // z_file->txt_data_so_far_bind_0.
    ASSERTW (all_txt_len == z_file->txt_data_so_far_bind || flag.data_modified || flag.make_reference, 
             "Hmm... incorrect calculation for %s sizes: total section sizes=%s but file size is %s (diff=%d)", 
             dt_name (z_file->data_type), str_uint_commas (all_txt_len).s, str_uint_commas (z_file->txt_data_so_far_bind).s, 
             (int32_t)(z_file->txt_data_so_far_bind - all_txt_len)); 

    zfile_compress_section_data (evb, SEC_STATS, &z_file->stats_buf);
    zfile_compress_section_data (evb, SEC_STATS, &z_file->STATS_buf);

    buf_free (&count_per_section_buf);
}

void stats_display (void)
{
    Buffer *buf = flag.show_stats == 1 ? &z_file->stats_buf : &z_file->STATS_buf;

    if (!buf_is_alloc (buf)) return; // no stats available

    buf_print (buf , false);

    Section sl = sections_last_sec (SEC_STATS, false) - 1; // first stats section 

    if (z_file->disk_size < (1<<20))  // no need to print this note if total size > 1MB, as the ~2K of overhead is rounded off anyway
        // stats text doesn't include SEC_STATS and SEC_GENOZIP_HEADER - the last 3 sections in the file - since stats text is generated before these sections are compressed
        iprintf ("\nNote: ZIP total file size excludes overhead of %s\n", str_size (z_file->disk_size - sl->offset).s);

    iprint0 ("\n");
}

void stats_read_and_display (void)
{
    Section sl = sections_last_sec (SEC_STATS, true);
    if (!sl) return; // genozip file does not contain stats sections (SEC_STATS was introduced in v 7.0.5)

    // read and uncompress the requested stats section
    zfile_get_global_section (SectionHeader, SEC_STATS, sl - (flag.show_stats==1),
                              flag.show_stats == 1 ? &z_file->stats_buf : &z_file->STATS_buf, "z_file->stats_buf");
    
    stats_display();

    if (exe_type == EXE_GENOCAT) exit_ok(); // if this is genocat - we're done
}

// concatenate txt names of bound files so we can show them all
void stats_add_txt_name (const char *fn)
{
    bufprintf (evb, &z_file->bound_txt_names, "%s%s", z_file->bound_txt_names.len ? " ": "", fn);
    z_file->bound_txt_names.param++; // we store the number of files in param
}
