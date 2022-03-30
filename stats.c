// ------------------------------------------------------------------
//   stats.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
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
#include "gencomp.h"

typedef struct {
    DidIType my_did_i, st_did_i;
    int64_t txt_len, z_size;
    char name[100];
    rom type;
    StrText did_i, words, hash, uncomp_dict, comp_dict, comp_b250, comp_data;
    float pc_of_txt, pc_of_z, pc_dict, pc_in_local, pc_failed_singletons, pc_hash_occupancy;
} StatsByLine;

// store the sizes of dict / b250 / local in zctx->*.param, and of other sections in sbl[st].z_size
static void stats_get_compressed_sizes (StatsByLine *sbl, CompIType comp_i) 
{
    ContextIndex ctx_index[MAX_DICTS];

    // initialize & prepare context index
    for (DidIType did_i=0; did_i < z_file->num_contexts; did_i++) {
        ZCTX(did_i)->b250.len = ZCTX(did_i)->b250.count; // number of b250 words - move to len
        ZCTX(did_i)->b250.count = ZCTX(did_i)->dict.count = ZCTX(did_i)->local.count = 0; 
        ctx_index[did_i] = (ContextIndex){ .did_i = did_i, .dict_id = ZCTX(did_i)->dict_id };
    }
    
    qsort (ctx_index, z_file->num_contexts, sizeof (ContextIndex), sort_by_dict_id);

    for (Section sec = section_next(0); sec; sec = section_next (sec)) {

        if (comp_i != COMP_NONE && comp_i != sec->comp_i) continue;

        if (sec->st != SEC_B250 && sec->st != SEC_DICT && sec->st != SEC_LOCAL && sec->st != SEC_COUNTS) 
            sbl[sec->st].z_size += sec->size;
            
        else {
            DidIType did_i = ctx_get_existing_did_i_do (sec->dict_id, z_file->contexts, z_file->dict_id_to_did_i_map,
                                                        ctx_index, z_file->num_contexts);

            // accumulate z_size for its context in its local/b250/dict.param
            switch (sec->st) {
                case SEC_LOCAL  : ZCTX(did_i)->local.count += sec->size; break;
                case SEC_B250   : ZCTX(did_i)->b250.count  += sec->size; break;
                case SEC_COUNTS :
                case SEC_DICT   : ZCTX(did_i)->dict.count  += sec->size; break;
                default         : break;
            }
        }
    }

    // update allocations of compressed sizes between contexts
    if (DTPZ(stats_reallocate)) DTPZ(stats_reallocate)();
}

static void stats_output_file_metadata (Buffer *buf)
{
    bufprint0 (evb, buf, "\n\n");
    if (txt_file->name) 
        bufprintf (evb, buf, "%s file%s%s: %.*s\n", dt_name (z_file->data_type), 
                   z_file->bound_txt_names.count > 1 ? "s" : "", // param holds the number of txt files
                   flag.pair ? " (paired)" : "",
                   (int)z_file->bound_txt_names.len, z_file->bound_txt_names.data);
    
    if (flag.reference == REF_MAKE_CHAIN) {
        bufprintf (evb, buf, "PRIM reference: %s MD5=%s genozip_version=%u\n", ref_get_filename (prim_ref), digest_display (ref_get_file_md5 (prim_ref)).s, ref_get_genozip_version (prim_ref));
        bufprintf (evb, buf, "LUFT reference: %s MD5=%s genozip_version=%u\n", ref_get_filename (gref), digest_display (ref_get_file_md5 (gref)).s, ref_get_genozip_version (gref));
    }

    else if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE || flag.reference == REF_LIFTOVER) 
        bufprintf (evb, buf, "Reference: %s MD5=%s genozip_version=%u\n", ref_get_filename (gref), digest_display (ref_get_file_md5 (gref)).s, ref_get_genozip_version (gref));

    if (Z_DT(DT_VCF)) 
        bufprintf (evb, buf, "Samples: %u   ", vcf_header_get_num_samples());

    uint32_t num_used_ctxs=0;
    for (ContextP ctx=ZCTX(0); ctx < ZCTX(z_file->num_contexts); ctx++)
        if (ctx->nodes.len || ctx->txt_len) num_used_ctxs++;

    if (Z_DT(DT_SAM) || Z_DT(DT_BAM))
        bufprintf (evb, buf, "%ss: %s (Prim: %s Supp/Secn: %s)  Contexts: %u   Vblocks: %u x %u MB  Sections: %u\n", 
                   DTPZ (line_name), str_int_commas (z_file->num_lines).s, str_int_commas (gencomp_get_num_lines (SAM_COMP_PRIM)).s, 
                   str_int_commas (gencomp_get_num_lines (SAM_COMP_DEPN)).s, num_used_ctxs, 
                   z_file->num_vbs, (uint32_t)(segconf.vb_size >> 20), (uint32_t)z_file->section_list_buf.len);

    else if (z_is_dvcf)
        bufprintf (evb, buf, "%ss: %s (Prim-only: %s Luft-only: %s)  Contexts: %u   Vblocks: %u x %u MB  Sections: %u\n", 
                   DTPZ (line_name), str_int_commas (z_file->num_lines).s, str_int_commas (gencomp_get_num_lines (VCF_COMP_PRIM_ONLY)).s, 
                   str_int_commas (gencomp_get_num_lines (VCF_COMP_LUFT_ONLY)).s, num_used_ctxs, 
                   z_file->num_vbs, (uint32_t)(segconf.vb_size >> 20), (uint32_t)z_file->section_list_buf.len);

    else
        bufprintf (evb, buf, "%ss: %s   Contexts: %u   Vblocks: %u x %u MB  Sections: %u\n", 
                   DTPZ (line_name), str_int_commas (z_file->num_lines).s, num_used_ctxs, 
                   z_file->num_vbs, (uint32_t)(segconf.vb_size >> 20), (uint32_t)z_file->section_list_buf.len);

    if (Z_DT(DT_KRAKEN)) {
        int64_t dominant_taxid_count;
        rom dominant_taxid = ctx_get_snip_with_largest_count (KRAKEN_TAXID, &dominant_taxid_count);

        if (dominant_taxid_count != -1)
            bufprintf (evb, buf, "Dominant TaxID: %s  %s: %s (%-5.2f%%)\n", dominant_taxid, DTPZ (line_name),
                       str_int_commas (dominant_taxid_count).s, 100.0 * (float)dominant_taxid_count / (float)z_file->num_lines); 
        else
            bufprint0 (evb, buf, "Dominant TaxID: No dominant species\n"); 
    }  
    
    else if (kraken_is_loaded) 
        bufprint0 (evb, buf, "Features: Per-line taxonomy ID data\n");

    else if (Z_DT(DT_CHAIN) && flag.reference == REF_MAKE_CHAIN && !segconf.chain_mismatches_ref)
        bufprint0 (evb, buf, "Features: Chain file suitable for use with genozip --chain\n");

    if (Z_DT(DT_VCF)) {
        if (z_is_dvcf)  
            bufprintf (evb, buf, "Features: Dual-coordinates: Main VBs: %u Prim-only VBs: %u Luft-only VBs: %u\n", 
                       sections_get_num_vbs(VCF_COMP_MAIN), sections_get_num_vbs(VCF_COMP_PRIM_ONLY), sections_get_num_vbs(VCF_COMP_LUFT_ONLY));
    }

    if (Z_DT(DT_SAM) || Z_DT(DT_BAM)) {
        if (segconf.sam_is_sorted)   
            bufprint0 (evb, buf, "Sorting: Sorted by POS\n");
        
        if (segconf.sam_is_collated) 
            bufprint0 (evb, buf, "Sorting: Collated by QNAME\n");
        
        if (!segconf.sam_is_sorted && !segconf.sam_is_collated) 
            bufprint0 (evb, buf, "Sorting: Not sorted or collated\n");
        
        if (z_has_gencomp) 
            bufprintf (evb, buf, "Main VBs: %u Prim VBs: %u Depn VBs: %u\n", 
                       sections_get_num_vbs(SAM_COMP_MAIN), sections_get_num_vbs(SAM_COMP_PRIM), sections_get_num_vbs(SAM_COMP_DEPN));
    }

    if (Z_DT(DT_FASTA))
        bufprintf (evb, buf, "Sequence type: %s\n", segconf.seq_type==SQT_AMINO ? "Amino acids" : "Nucleotide bases");

    if (segconf.qname_flavor) 
        bufprintf (evb, buf, "Read name style: %s%s%s\n", 
                   qf_name(segconf.qname_flavor), segconf.qname_flavor2 ? " + " : "", segconf.qname_flavor2 ? qf_name(segconf.qname_flavor2) : "");

    bufprintf (evb, buf, "Genozip version: %s %s\nDate compressed: %s\n", 
               GENOZIP_CODE_VERSION, arch_get_distribution(), str_time().s);

    bufprint0 (evb, buf, "Command line: ");
    buf_add_string (evb, buf, flags_command_line()); // careful not to use bufprintf with command_line as it can exceed the maximum length in bufprintf

    if (license_has_details())
        bufprintf (evb, buf, "\n%s\n", license_get_one_line());
}

static DESCENDING_SORTER (stats_sort_by_z_size, StatsByLine, z_size)

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

void stats_set_consolidation (VBlockP vb, DidIType parent, unsigned num_deps, ...)
{
    va_list args;
    va_start (args, num_deps);

    for (unsigned d=0; d < num_deps; d++) {
        DidIType dep = (DidIType)va_arg (args, int); 
        if (dep != parent) 
            CTX(dep)->st_did_i = parent;
    }

    CTX(parent)->is_stats_parent = true;
    
    va_end (args);
}

void stats_set_consolidation_(VBlockP vb, DidIType parent, unsigned num_deps, ContextP *dep_ctxs)
{
    for (unsigned d=0; d < num_deps; d++)
        if (dep_ctxs[d]->did_i != parent) 
            dep_ctxs[d]->st_did_i = parent;

    CTX(parent)->is_stats_parent = true;
}

void stats_set_consolidationN (VBlockP vb, DidIType parent, DidIType first_dep, unsigned num_deps)
{
    for (ContextP ctx=CTX(first_dep); ctx < CTX(first_dep + num_deps); ctx++)
        if (ctx->did_i != parent) 
            ctx->st_did_i = parent;

    CTX(parent)->is_stats_parent = true;
}


static void stats_consolidate_non_ctx (StatsByLine *sbl, unsigned num_stats, rom consolidated_name, 
                                       unsigned num_deps, ...)
{
    va_list args;
    va_start (args, num_deps);

    rom deps[num_deps];
    for (unsigned d=0; d < num_deps; d++) 
        deps[d] = va_arg (args, rom );

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
    bufprintf (evb, &z_file->stats_buf, "NAME                   GENOZIP      %%       TXT      %%   RATIO\n%s", "");

    for (uint32_t i=0; i < num_stats; i++, s++)
        if (s->z_size)
            bufprintf (evb, &z_file->stats_buf, "%-20.20s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                       s->name, 
                       str_size (s->z_size).s, s->pc_of_z, // z size and % of total z that is in this line
                       str_size ((float)s->txt_len).s, s->pc_of_txt, // txt size and % of total txt which is in this line
                       (float)s->txt_len / (float)s->z_size); // ratio z vs txt

    if (src_comp_ratio != 1 && flag.show_stats != STATS_SHORT_COMP)
        bufprintf (evb, &z_file->stats_buf, 
                   "GENOZIP vs %-9s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                   codec_name (codec),
                   str_size (all_z_size).s, all_pc_of_z, // total z size and sum of all % of z (should be 100)
                   str_size (all_txt_len_0 / src_comp_ratio).s, all_pc_of_txt, // total txt fize and ratio z vs txt
                   all_comp_ratio / src_comp_ratio);
    
    // note: no point showing this per component in SAM and DVCF, bc all_txt_len_0 is fully accounted for in the MAIN component and it is 0 in the others
    if (!(z_has_gencomp && flag.show_stats == STATS_SHORT_COMP))
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
    if (flag.show_stats != STATS_SHORT_COMP && flag.show_stats != STATS_LONG_COMP) {
        bufprintf (evb, &z_file->STATS_buf, "\nSystem info: OS=%s cores=%u endianity=%s\n", 
                arch_get_os(), arch_get_num_cores(), arch_get_endianity());
        bufprintf (evb, &z_file->STATS_buf, "\nSections (sorted by %% of genozip file):%s\n", "");
    }

    bufprintf (evb, &z_file->STATS_buf, "did_i Name              Parent            #Words  Snips-(%% of #Words)    Hash-table    uncomp      comp      comp      comp      comp       txt    comp   %% of   %% of  %s\n", "");
    bufprintf (evb, &z_file->STATS_buf, "                                         in file   Dict  Local FailSton   Size Occp      dict      dict      b250     local     TOTAL             ratio    txt    zip%s\n", "");

    for (uint32_t i=0; i < num_stats; i++, s++)
        if (s->z_size)
            bufprintf (evb, &z_file->STATS_buf, "%-2.2s    %-17.17s %-17.17s %6s  %4.*f%%  %4.*f%%  %4.*f%% %8s %3.0f%% %9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
                       s->did_i.s, s->name, s->type, s->words.s, 
                       PC (s->pc_dict), s->pc_dict, PC(s->pc_in_local), s->pc_in_local, PC(s->pc_failed_singletons), s->pc_failed_singletons, 
                       s->hash.s, s->pc_hash_occupancy, // Up to here - these don't appear in the total
                       s->uncomp_dict.s, s->comp_dict.s, s->comp_b250.s, s->comp_data.s, str_size (s->z_size).s, 
                       str_size ((float)s->txt_len).s, (float)s->txt_len / (float)s->z_size, s->pc_of_txt, s->pc_of_z);

    // note: no point showing this per component in SAM and DVCF, bc all_txt_len_0 is fully accounted for in the MAIN component and it is 0 in the others
    if (!(z_has_gencomp && flag.show_stats == STATS_LONG_COMP))
        bufprintf (evb, &z_file->STATS_buf, "TOTAL                                                                               "
                "%9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
                str_size (all_uncomp_dict).s, str_size (all_comp_dict).s,  str_size (all_comp_b250).s, 
                str_size (all_comp_data).s,   str_size (all_z_size).s, str_size (all_txt_len).s, 
                all_comp_ratio, all_pc_of_txt, all_pc_of_z);
}

// generate the stats text - all sections except genozip header and the two stats sections 
void stats_generate (CompIType comp_i) // specific section, or COMP_NONE if for the entire file
{
    if (comp_i == COMP_NONE) {
        stats_output_file_metadata (&z_file->stats_buf);
        buf_copy (evb, &z_file->STATS_buf, &z_file->stats_buf, char,0,0, "z_file->STATS_buf");
    }

    int64_t all_comp_dict=0, all_uncomp_dict=0, all_comp_b250=0, all_comp_data=0, all_z_size=0, all_txt_len=0;

    // prepare data
    #define NUM_SBL (NUM_SEC_TYPES + z_file->num_contexts + 2) // 2 for consolidated groups
    static Buffer sbl_buf = EMPTY_BUFFER;
    ARRAY_alloc (StatsByLine, sbl, NUM_SBL, true, sbl_buf, evb, "sbl");

    #define ST_NAME(st) (&st_name(st)[4]) // cut off "SEC_" 

    stats_get_compressed_sizes (sbl, comp_i); // ctx sizes in ctx->local/dict/ctx.param, and other sections in sbl[st]->z_size.
    
    // use watermarks, to calculate amount attributable to each component
    static uint64_t z_size_watermark=0, txt_size_water_mark=0, txt_size_0_water_mark=0;
    if (comp_i == COMP_NONE || comp_i == COMP_MAIN) z_size_watermark = txt_size_water_mark = txt_size_0_water_mark = 0; // initialize

    uint64_t z_size   = z_file->disk_so_far - z_size_watermark;
    z_size_watermark  = z_file->disk_so_far;

    uint64_t txt_size = z_file->txt_data_so_far_bind - txt_size_water_mark;
    txt_size_water_mark = z_file->txt_data_so_far_bind;

    uint64_t txt_size_0 = z_file->txt_data_so_far_bind_0 - txt_size_0_water_mark;
    txt_size_0_water_mark = z_file->txt_data_so_far_bind_0;

    // non-contexts
    StatsByLine *s = sbl;
    for (SectionType st=0; st < NUM_SEC_TYPES; st++, s++) { 

        if (st == SEC_DICT || st == SEC_B250 || st == SEC_LOCAL || st == SEC_COUNTS) continue; // these are covered by individual contexts

        s->txt_len    = st == SEC_TXT_HEADER  ? z_file->txt_txtheader_so_far_bind : 0; // note: this excludes generated headers for DVCF and SAM/BAM
        s->type       = "Global";
        s->my_did_i   = s->st_did_i = DID_I_NONE;
        s->did_i.s[0] = s->words.s[0] = s->hash.s[0] = s->uncomp_dict.s[0] = s->comp_dict.s[0] = '-';
        s->pc_of_txt  = txt_size ? 100.0 * (float)s->txt_len / (float)txt_size : 0;
        s->pc_of_z    = z_size   ? 100.0 * (float)s->z_size  / (float)z_size   : 0;
        strcpy (s->name, ST_NAME (st)); 

        all_z_size    += s->z_size;
        all_comp_data += s->z_size;
        all_txt_len   += s->txt_len;
    }

    z_file->txt_txtheader_so_far_bind = 0; // reset (in case of showing components)
    
    // contexts
    for (DidIType did_i=0; did_i < z_file->num_contexts; did_i++) { 
        ContextP ctx = ZCTX(did_i);

        s->z_size = ctx->dict.count + ctx->b250.count + ctx->local.count;

        if (!ctx->b250.count && !ctx->txt_len && !ctx->b250.len && !ctx->is_stats_parent && !s->z_size) 
            continue;

        s->txt_len = ctx->txt_len;
        
        all_comp_dict   += ctx->dict.count;
        all_uncomp_dict += ctx->dict.len;
        all_comp_b250   += ctx->b250.count;
        all_comp_data   += ctx->local.count;
        all_z_size      += s->z_size;
        all_txt_len     += s->txt_len;

        if (Z_DT(DT_VCF) && dict_id_type(ctx->dict_id))
            sprintf (s->name, "%s/%s", dtype_name_z(ctx->dict_id), ctx->tag_name);
        else 
            strcpy (s->name, ctx->tag_name);

        s->type  = ctx->st_did_i != DID_I_NONE ? ZCTX(ctx->st_did_i)->tag_name : ctx->tag_name;

        // note: each VB contributes local.len contains its b250.count if it has it, and local_num_len if not 
        float n_words = ctx->local.len;

        s->my_did_i             = ctx->did_i;
        s->st_did_i             = ctx->st_did_i;
        s->did_i                = str_int_commas ((uint64_t)ctx->did_i); 
        s->words                = str_uint_commas_limit (n_words, 99999);
        s->pc_dict              = !n_words ? 0 : 100.0 * (float)MIN_(ctx->nodes.len, n_words) / n_words; // MIN_ is a workaround - not sure why nodes.len sometimes exceeds the dictionary words on the file (eg in TOPLEVEL)
        s->pc_in_local          = !n_words ? 0 : 100.0 * (float)ctx->local_num_words / n_words;
        s->pc_failed_singletons = !ctx->b250.count ? 0 : 100.0 * (float)ctx->num_failed_singletons / (float)ctx->b250.len;
        s->pc_hash_occupancy    = !ctx->global_hash_prime   ? 0 : 100.0 * (float)(ctx->nodes.len + ctx->ol_nodes.len) / (float)ctx->global_hash_prime;
        s->hash                 = str_size (ctx->global_hash_prime);
        s->uncomp_dict          = str_size (ctx->dict.len);
        s->comp_dict            = str_size (ctx->dict.count);
        s->comp_b250            = str_size (ctx->b250.count);
        s->comp_data            = str_size (ctx->local.count);
        s->pc_of_txt            = txt_size ? 100.0 * (float)s->txt_len / (float)txt_size : 0;
        s->pc_of_z              = z_size   ? 100.0 * (float)s->z_size  / (float)z_size   : 0;
    
        ctx->b250.count = ctx->dict.count = ctx->local.count = ctx->txt_len = 0; // reset

        s++; // increment only if it has some data, otherwise already continued
    }

    unsigned num_stats = s - sbl;

    // note: for txt size and compression ratio in the TOTAL line (all_comp_ratio) we use txt_data_so_far_bind_0 
    // (the original txt data size) and not all_txt_size (size after ZIP modifications like --optimize). 
    // Therefore, in case of ZIP-modified txt, the sum of the (modified) fields in the TXT column will NOT equal the
    // TOTAL in the TXT column. That's ok.
    float all_comp_ratio = (float)txt_size_0 /* without modifications */ / (float)all_z_size;
    float all_pc_of_txt  = txt_size ? 100.0 * (float)all_txt_len / (float)txt_size : 0 /* with modifications */;
    float all_pc_of_z    = z_size   ? 100.0 * (float)all_z_size  / (float)z_size   : 0;

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

    stats_consolidate_non_ctx (sbl, num_stats, "Other", 22 + (DTPZ(txt_header_required) == HDR_NONE), "E1L", "E2L", "EOL", 
                               "SAMPLES", "AUX", TOPLEVEL, "ToPLUFT", "TOP2BAM", "TOP2FQ", "TOP2FQEX", "TOP2VCF", "TOP2HASH", 
                               "LINEMETA", "CONTIG", "COORDS", "SAGROUP", "SAALN",
                               ST_NAME (SEC_RANDOM_ACCESS), ST_NAME (SEC_DICT_ID_ALIASES), ST_NAME (SEC_RECON_PLAN),
                               ST_NAME (SEC_VB_HEADER), ST_NAME (SEC_BGZF), ST_NAME(SEC_TXT_HEADER)/*must be last*/);
    
    ASSERTW (all_txt_len == txt_size || flag.make_reference, // all_txt_len=0 in make-ref as there are no contexts
             "Expecting all_txt_len=%"PRId64" == txt_size=%"PRId64, all_txt_len, txt_size);

    // short form stats from --stats    
    qsort (sbl, num_stats, sizeof (sbl[0]), stats_sort_by_z_size);  // re-sort after consolidation

    // source compression, eg BGZF, against txt before any modifications
    float src_comp_ratio = (float)txt_size_0 / 
                           (float)((flag.bind != BIND_FQ_PAIR && txt_file->disk_size) ? txt_file->disk_size : z_file->txt_disk_so_far_bind); 

    stats_output_stats (sbl, num_stats, src_comp_ratio, txt_file->source_codec, all_txt_len, txt_size_0, all_z_size, all_pc_of_txt, all_pc_of_z, all_comp_ratio);
    
    // if we're showing stats of a single components - output it now
    if (comp_i != COMP_NONE) {
        iprintf ("\n\nComponent=%s:\n", comp_name (comp_i));
        buf_print (flag.show_stats == STATS_SHORT_COMP ? &z_file->stats_buf : &z_file->STATS_buf, false);
        buf_free (z_file->stats_buf);
        buf_free (z_file->STATS_buf);
    }

    // note: we use txt_data_so_far_bind is the sum of recon_sizes - see zip_update_txt_counters - which is
    // expected to be the sum of txt_len. However, this NOT the size of the original file which is stored in
    // z_file->txt_data_so_far_bind_0.
    ASSERTW (comp_i != COMP_NONE || all_txt_len == txt_size || flag.data_modified || flag.make_reference, 
             "Hmm... incorrect calculation for %s sizes: total section sizes=%s but file size is %s (diff=%d)", 
             dt_name (z_file->data_type), str_int_commas (all_txt_len).s, str_int_commas (txt_size).s, 
             (int32_t)(txt_size - all_txt_len)); 

    if (comp_i == COMP_NONE) {
        zfile_compress_section_data (evb, SEC_STATS, &z_file->stats_buf);
        zfile_compress_section_data (evb, SEC_STATS, &z_file->STATS_buf);
    }

    buf_free (sbl_buf);
}

void stats_display (void)
{
    Buffer *buf = flag.show_stats == 1 ? &z_file->stats_buf : &z_file->STATS_buf;

    if (!buf_is_alloc (buf)) return; // no stats available

    buf_print (buf , false);

    Section sec = sections_last_sec (SEC_STATS, false) - 1; // first stats section 

    if (z_file->disk_size < (1<<20))  // no need to print this note if total size > 1MB, as the ~2K of overhead is rounded off anyway
        // stats text doesn't include SEC_STATS and SEC_GENOZIP_HEADER - the last 3 sections in the file - since stats text is generated before these sections are compressed
        iprintf ("\nNote: ZIP total file size excludes overhead of %s\n", str_size (z_file->disk_size - sec->offset).s);

    iprint0 ("\n");
}

void stats_read_and_display (void)
{
    Section sec = sections_last_sec (SEC_STATS, true);
    if (!sec) {
        iprint0 ("No stats available for this file.\n"); // possibly the file was compressed with --cstats or --CSTATS
        return; // genozip file does not contain stats sections (SEC_STATS was introduced in v 7.0.5)
    }

    // read and uncompress the requested stats section
    zfile_get_global_section (SectionHeader, sec - (flag.show_stats==1),
                              flag.show_stats == 1 ? &z_file->stats_buf : &z_file->STATS_buf, "z_file->stats_buf");
    
    stats_display();

    if (exe_type == EXE_GENOCAT) exit_ok(); // if this is genocat - we're done
}

// concatenate txt names of bound files so we can show them all
void stats_add_txt_name (rom fn)
{
    bufprintf (evb, &z_file->bound_txt_names, "%s%s", z_file->bound_txt_names.len ? " ": "", fn);
    z_file->bound_txt_names.count++; // number of files
}
