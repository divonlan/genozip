// ------------------------------------------------------------------
//   stats.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <stdarg.h>
#include <errno.h>
#include "stats.h"
#include "file.h"
#include "txtheader.h"
#include "zfile.h"
#include "arch.h"
#include "codec.h"
#include "license.h"
#include "qname.h"
#include "gencomp.h"
#include "url.h"
#include "crypt.h"
#include "tar.h"
#include "contigs.h"
#include "mgzip.h"
#include "huffman.h"
#include "sorter.h"
#include "fastq.h"

#define SHORT_HEADER "NAME                   GENOZIP      %       TXT      %   RATIO\n"

#define LONG_HEADER "did_i Name              Parent           ----------- #Words -----------     Hash-table        uncomp      comp      comp    uncomp      comp      comp       txt    comp   % of   % of\n" \
                    "                                         in file   Dict  Local FailSton   Size Occp SOccp     dict      dict      b250     local     local     TOTAL             ratio    txt    zip\n"

static Buffer stats={}, STATS={};
Buffer sbl_buf={}, features={}, exceptions={}, internals={};
Buffer stats_programs = {}; // data-type specific programs (eg @PG for SAM/BAM FORMAT/INFO tags for VCF)
int64_t all_txt_len=0;
float src_comp_ratio=0, all_comp_ratio=0;

void stats_add_one_program (STRp(prog_name))
{
    stats_programs.name = "stats_programs"; // initialize name (note: destroyed between files, so needs to be reinitialized)

    buf_append (evb, stats_programs, char, prog_name, prog_name_len, "stats_programs");

    BNXTc (stats_programs) = ';'; // note: buf_append allocates one extra character
}

rom stats_find_in_programs (rom signature)
{
    if (!stats_programs.len) return false;

    SAFE_NULB(stats_programs);
    rom found = strstr (B1STc(stats_programs), signature);
    SAFE_RESTORE;

    return found;
}

void stats_remove_data_after_program_name (rom program)
{
    rom prog = stats_find_in_programs (program);

    if (!program || (prog != B1STc (stats_programs) && prog[-1] != '\n')) return;

    rom sep = strpbrk (prog, " \t;"); // note: programs always end with a ;
    if (*sep == ';') return; // no data after program name 

    rom end = strchr (sep, ';');

    buf_remove (stats_programs, char, BNUM (stats_programs, sep), end - sep);
}

// substitute ; and , with their UTF-8 equivalents, replace '\' with '\\'
static StrText1K stats_subs_seps_in_name (rom name)
{
    ASSERT0 (sizeof (StrText) > MAX_TAG_LEN - 10, "bad string size");

    StrText1K s = {};
    char *next = s.s;

    // similar escaping as in segconf_get_qual_histo
    while (*name) {
        switch (*name) {
            case ','  : memcpy (next, "⸲",    STRLEN("⸲"));     next += STRLEN("⸲");     break; // Unicode "Turned Comma"
            case ';'  : memcpy (next, "；",   STRLEN("；"));   next += STRLEN("；");    break; // Unicode "Full-width" semicolon
            case '\\' : memcpy (next, "\\\\", STRLEN("\\\\")); next += STRLEN("\\\\"); break; // a backslash must be escaped (can exist in channel context names of long muteces)
            case '"'  : memcpy (next, "\\\"", STRLEN("\\\"")); next += STRLEN("\\\""); break;
            default   : *next++ = *name;
        }
        name++;
    }

    return s;
}

// calculate exceptions before consolidating stats
static void stats_calc_hash_occ (StatsByLine *sbl, unsigned num_stats)
{
    int need_sep=0;

    for (StatsByLine *s=sbl, *after=sbl + num_stats; s < after; s++) {
        ContextP zctx = (s->my_did_i != DID_NONE) ? ZCTX(s->my_did_i) : NULL;

        if (s->pc_hash_occupancy > 100 || // > 100%
            (zctx && zctx->nodes.len > 1000000 && z_file->num_lines / zctx->nodes.len < 16)) { // more than 10% of num_lines (and at least 1M nodes)
             
            // in case of an over-populated hash table, we send the first 3 and last 3 words in the dictionary, which will help debugging the issue
            bufprintf (evb, &exceptions, "%s%s,%s,%s,%u%%", 
                       need_sep++ ? ";" : "", s->name, s->type, s->hash.s, (int)s->pc_hash_occupancy);
            
            uint32_t n_words = zctx->nodes.len32; // note: this can be a low number despite pc_hash_occupancy being large - if words ended up as singletons
            WordIndex words[NUM_COLLECTED_WORDS] = { 0, 1, 2, n_words-3, n_words-2, n_words-1 }; // first three and last threewords in the the dictionary of this field
            for (int i=0; i < MIN_(NUM_COLLECTED_WORDS, n_words); i++) {
                STR(snip);
                ctx_get_z_snip_ex (zctx, words[i], pSTRa(snip));
                
                char *s = str_snip_ex (DT_NONE, STRa(snip), false).s; 
                int s_len = strlen (s);
                
                #define MAX_LEN_EXECP_SNIP 100 // limit chars per snip
                if (s_len > MAX_LEN_EXECP_SNIP) {  
                    s[MAX_LEN_EXECP_SNIP] = 0; 
                    s_len = MAX_LEN_EXECP_SNIP;
                }

                bufprintf (evb, &exceptions, ",%s", str_replace_letter (STRa(s), ',', -127)); 
            }
        }
    }

    // we send the first 6 qnames of unrecognized QNAME flavor or in case of existing by unrecognized QNAME2 flavor
    for (QType q=QNAME1; q < NUM_QTYPES; q++) {
        if (z_file->n_1st_flav_qnames[q]) {

            bufprintf (evb, &exceptions, "%s%s,,", need_sep++ ? ";" : "", qtype_name(q));

            for (int i=0; i < z_file->n_1st_flav_qnames[q]; i++)
                if (z_file->unk_flav_qnames[q][i][0])
                    bufprintf (evb, &exceptions, ",%s", stats_subs_seps_in_name (z_file->unk_flav_qnames[q][i]).s); 
        }
    }

    for (int id_i=0; id_i < NUM_UNK_ID_CTXS && z_file->unk_ids_tag_name[id_i][0]; id_i++) {
        bufprintf (evb, &exceptions, "%s%s,,,", need_sep++ ? ";" : "", z_file->unk_ids_tag_name[id_i]);

        for (int i=0; i < NUM_COLLECTED_WORDS && z_file->unk_ids[id_i][i][0]; i++)
            bufprintf (evb, &exceptions, ",%s", str_replace_letter (z_file->unk_ids[id_i][i], strlen(z_file->unk_ids[id_i][i]), ',', -127)); 
    }

    if (z_file->sam_malformed_XA[0]) 
        bufprintf (evb, &exceptions, "%sBAD_XA,%s", (need_sep++ ? ";" : ""), str_replace_letter (z_file->sam_malformed_XA, strlen(z_file->sam_malformed_XA), ',', -127)); 

    // nonbio linkers, but unrecognized format
    int n_linkers = fastq_nonbio_get_n_linkers();
    if (n_linkers && !segconf.nonbio_type) {
        bufprintf (evb, &exceptions, "%sNonbio_Linkers,,", need_sep++ ? ";" : "");

        for (int i=0; i < n_linkers + 1/*umi*/; i++)
            bufprintf (evb, &exceptions, ",%s", fastq_nonbio_get_linker_for_stats(i).s); 
    }
}

// stats of contexts and sections contribution to Z
static void stats_prepare_internals (StatsByLine *sbl, unsigned num_stats, uint64_t all_z_size)
{
    #define MAX_INTERNALS 128 
    if (!all_z_size) return;
    
    int need_sep = 0;

    // note: sbl is sorted by z_size
    int count = 0;
    for (StatsByLine *after=sbl + num_stats; sbl < after && count < MAX_INTERNALS; sbl++) 
        if (sbl->z_size) {
            bufprintf (evb, &internals, "%s%s,%s,%.1f%%,", 
                       need_sep++ ? ";" : "", stats_subs_seps_in_name (sbl->name).s, sbl->type, sbl->pc_of_z);

            if (sbl->my_did_i != DID_NONE) {
                ContextP zctx = ZCTX(sbl->my_did_i);

                bufprintf (evb, &internals, "%s,%.1f%%,",
                           sbl->lcodec, percent (zctx->local.count, all_z_size));  // local

                bufprintf (evb, &internals, "%s,%.1f%%,",
                           sbl->bcodec, percent (zctx->b250.count, all_z_size));  // b250

                bufprintf (evb, &internals, "%s,%.1f%%,%u,",
                           sbl->dcodec, percent (zctx->dict.count, all_z_size), zctx->nodes.len32); // dict
            }

            count++; // count only those with z_size
        }
}

// store the sizes of dict / b250 / local in zctx->*.param, and of other sections in sbl[st].z_size
static void stats_get_compressed_sizes (StatsByLine *sbl) 
{
    for_zctx {
        zctx->b250.len   = zctx->b250.count; // number of b250 words - move to len
        zctx->b250.count = zctx->dict.count = zctx->local.count = 0; 
    }
    
    ctx_create_ctx_index (evb, &z_file->ca);    

    for (Section sec = section_next(0); sec; sec = section_next (sec)) {

        if (flag.show_stats_comp_i == COMP_NONE && !IS_DICTED_SEC (sec->st))
            sbl[sec->st].z_size += sec->size;
            
        else if (flag.show_stats_comp_i != COMP_NONE && flag.show_stats_comp_i == sec->comp_i && 
            (sec->st==SEC_VB_HEADER || sec->st==SEC_TXT_HEADER || sec->st==SEC_RECON_PLAN))
            sbl[sec->st].z_size += sec->size;
            
        else if ((flag.show_stats_comp_i == COMP_NONE || flag.show_stats_comp_i == sec->comp_i) && IS_DICTED_SEC (sec->st)) {
            Did did_i = ctx_get_existing_did_i_do (sec->dict_id, &z_file->ca);

            ASSERT (did_i != DID_NONE, "Cannot find zctx for %s with dict_id=%s\n", st_name(sec->st), dis_dict_id (sec->dict_id).s);

            // accumulate z_size for its context in its local/b250/dict.param
            switch (sec->st) {
                case SEC_LOCAL    : ZCTX(did_i)->local.count += sec->size; break;
                case SEC_B250     : ZCTX(did_i)->b250.count  += sec->size; break;
                case SEC_COUNTS   : // include z_size of COUNTS, SUBDICTS, HUFFMAN in dict column of stats
                case SEC_SUBDICTS :
                case SEC_HUFFMAN  :
                case SEC_DICT     : ZCTX(did_i)->dict.count  += sec->size; break;
                default           : break;
            }
        }
    }

    // update allocations of compressed sizes between contexts
    if (DTPZ(stats_reallocate)) DTPZ(stats_reallocate)();
}

static void stats_output_file_metadata (void)
{
    #define FEATURE0(cond,stats_str,features_str)                       \
    if (cond) {                                                         \
        bufprint0 (evb, &stats, stats_str "\n");                        \
        bufprint0 (evb, &features, features_str ";");                   \
    }

    #define FEATURE(cond,stats_format,features_format,...)              \
    if (cond) {                                                         \
        bufprintf (evb, &stats, stats_format "\n", __VA_ARGS__);        \
        bufprintf (evb, &features, features_format ";", __VA_ARGS__);   \
    }

    bufprintf (evb, &features, "VBs=%u X %s;", z_file->num_vbs, str_size (segconf.vb_size).s);
    if (!Z_DT(GNRIC) && !flag.make_reference) {
        bufprintf (evb, &features, "num_lines=%"PRIu64";", z_file->num_lines);

        if (z_file->num_components > 1) {
            bufprint0 (evb, &features, "comp_num_lines=");
            for (CompIType comp_i=0; comp_i < z_file->num_components; comp_i++)
                bufprintf (evb, &features, "%"PRId64"+", z_file->comp_num_lines[comp_i]); 

            features.len32 -= 1; // remove final "+"
            bufprint0 (evb, &features, ";");
        }
    }
        
    bufprint0 (evb, &stats, "\n\n");
    if (txt_file->name) 
        bufprintf (evb, &stats, "%s%s%s file%s%s: %.*s\n", 
                   z_dt_name(),
                   flag.deep ? "/" : "",
                   flag.deep ? (segconf.fasta_as_fastq ? "FASTA" : "FASTQ") : "",
                   s_or_nil(z_file->bound_txt_names.count), // param holds the number of txt files
                   (flag.deep && flag.pair)?" (deep & paired)" : flag.deep?" (deep)" : flag.pair?" (paired)" : "", 
                   (int)z_file->bound_txt_names.len, z_file->bound_txt_names.data);
    
    if (IS_REF_EXTERNAL || IS_REF_EXT_STORE) 
        bufprintf (evb, &stats, "Reference: %s %s=%s ref_genozip_ver=%s\n", 
                   ref_get_filename(),  digest_alg_name (ref_get_genome_digest_alg()), 
                   digest_display_(ref_get_genome_digest(), ref_get_genome_digest_alg()).s, 
                   STRver(ref_get_genozip_ver()).s);

    uint32_t num_used_ctxs=0;
    for_zctx_that (zctx->nodes.len || zctx->txt_len) num_used_ctxs++;

    #define REPORT_VBs ({                                                                                   \
        bufprintf (evb, &stats, "%ss: %s   Contexts: %u   Vblocks: %u x %s   Sections: %u\n",               \
                   DTPZ (line_name), str_int_commas (z_file->num_lines).s, num_used_ctxs,                   \
                   z_file->num_vbs, TXT_IS_VB_SIZE_BY_MGZIP         ? "(var-length)"                        \
                                  : (segconf.vb_size % (1 MB) != 0) ? str_int_commas (segconf.vb_size).s    \
                                  :                                   str_size (segconf.vb_size).s, z_file->section_list.len32); })

    #define REPORT_QNAME                                                                                    \
        FEATURE (z_file->num_lines, "Read name style: %s%s%s%s", "QNAME=%s%s%s%s",                          \
                 segconf_qf_name (QNAME1),                                                                  \
                 cond_str(segconf.qname_flavor[QNAME2], "+", segconf_qf_name (QNAME2)), /* no space surrounding the '+' as expected by batch_qname_flavors */ \
                 cond_str(segconf.qname_flavor[QLINE3], "+", segconf_qf_name (QLINE3)),                     \
                 cond_str(segconf.qname_flavor[QEMBED], "+", segconf_qf_name (QEMBED)));

    switch (z_file->data_type) {
        case DT_SAM:
        case DT_BAM: {
            uint64_t num_alignments = z_file->comp_num_lines[SAM_COMP_MAIN] + z_file->comp_num_lines[SAM_COMP_PRIM] + z_file->comp_num_lines[SAM_COMP_DEPN]; // excluding Deep FQ components
            unsigned num_fq_files = MAX_(0, (int)z_file->num_components - SAM_COMP_FQ00); 
            uint64_t num_fq_reads = 0;
            for (CompIType comp_i=SAM_COMP_FQ00; comp_i < z_file->num_components; comp_i++)
                num_fq_reads += z_file->comp_num_lines[comp_i];

            double deep_pc = percent (z_file->deep_stats[NDP_DEEPABLE] + z_file->deep_stats[NDP_DEEPABLE_TRIM], z_file->deep_stats[NDP_FQ_READS]);

            if (z_has_gencomp) 
                bufprintf (evb, &stats, "%s %ss: %s (in Prim VBs: %s in Depn VBs: %s)\n", 
                           z_dt_name(),
                           DTPZ (line_name), str_int_commas (num_alignments).s, str_int_commas (z_file->comp_num_lines[SAM_COMP_PRIM]).s, 
                           str_int_commas (z_file->comp_num_lines[SAM_COMP_DEPN]).s);
            else
                bufprintf (evb, &stats, "%s %ss: %s\n", 
                           z_dt_name(), DTPZ (line_name), str_int_commas (num_alignments).s);

            if (flag.deep)
                bufprintf (evb, &stats, "FASTQ %ss: %s in %u FASTQ files\n", 
                           dt_props[DT_FASTQ].line_name, str_int_commas (z_file->deep_stats[NDP_FQ_READS]).s, num_fq_files);

            bufprintf (evb, &stats, "Contexts: %u  Vblocks: %u x %s  Sections: %s\n", 
                       num_used_ctxs, z_file->num_vbs, (segconf.vb_size % (1 MB) != 0) ? str_int_commas (segconf.vb_size).s : str_size (segconf.vb_size).s, str_int_commas (z_file->section_list.len).s);

            uint32_t num_deep_fq_vbs = 0;
            if (flag.deep) 
                for (CompIType comp_i = SAM_COMP_FQ00; comp_i < z_file->num_components; comp_i++)
                    num_deep_fq_vbs += sections_get_num_vbs (comp_i);

            if (z_has_gencomp) {
                bufprintf (evb, &stats, "Main VBs: %u Prim VBs: %u Depn VBs: %u", 
                           sections_get_num_vbs(SAM_COMP_MAIN), sections_get_num_vbs(SAM_COMP_PRIM), sections_get_num_vbs(SAM_COMP_DEPN));

                if (flag.deep) 
                    bufprintf (evb, &stats, " FASTQ VBs: %u", num_deep_fq_vbs);

                bufprint0 (evb, &stats, "\n");
            }
            else if (flag.deep) 
                bufprintf (evb, &stats, "%s VBs: %u FASTQ VBs: %u\n", z_dt_name(), sections_get_num_vbs(SAM_COMP_MAIN), num_deep_fq_vbs);
            else
                REPORT_VBs;

            if (sam_num_header_contigs()) bufprintf (evb, &features, "hdr_contigs=%u (%"PRIu64");", sam_num_header_contigs(), contigs_get_nbases (sam_hdr_contigs));
            if (IS_REF_LOADED_ZIP) bufprintf (evb, &features, "ref_contigs=%u (%"PRIu64");", ref_get_ctgs()->contigs.len32, contigs_get_nbases (ref_get_ctgs()));
            if (flag.deep) {
                bufprintf (evb, &features, "deep=%.1f%%;", deep_pc);                
                bufprintf (evb, &features, "deep_qtype=%s;", qtype_name (segconf.deep_qtype));
                bufprintf (evb, &features, "deep_trim=%s;", segconf_deep_trimming_name());
                bufprintf (evb, &features, "deep_qual=%s;", TF (!segconf.deep_no_qual));
                
                if (segconf.deep_N_sam_score && segconf.deep_N_fq_score) 
                    bufprintf (evb, &features, "deep_N_qual=%s%c→%c;", 
                            segconf.deep_N_fq_score == '\'' ? "'" : "", // escape a ' at the beginning of a cell
                            segconf.deep_N_fq_score, segconf.deep_N_sam_score);

                if (segconf.is_interleaved) bufprint0 (evb, &features, "interleaved=True;"); // fastq is interleaved
            }

            if (num_alignments) {
                FEATURE0 (segconf.is_sorted && !segconf.sam_is_unmapped, "Sorting: Sorted by POS", "sorted_by=POS");        
                FEATURE0 (segconf.is_sorted && segconf.sam_is_unmapped, "Sorting: Unmapped", "sorted_by=Unmapped");        
                FEATURE0 (segconf.is_collated, "Sorting: Collated by QNAME", "sorted_by=QNAME");
                FEATURE0 (!segconf.is_sorted && !segconf.is_collated, "Sorting: Not sorted or collated", "sorted_by=NONE");
                FEATURE0 (segconf.multiseq, "Multiseq", "multiseq=True");

                if (segconf.is_minimap2 || segconf.CIGAR_has_eqx)    bufprintf (evb, &features, "CIGAR_has_eqx=%s;",    TF(segconf.CIGAR_has_eqx));
                if (segconf.is_minimap2 || segconf.SA_NM_by_CIGAR_X) bufprintf (evb, &features, "SA_NM_by_CIGAR_X=%s;", TF(segconf.SA_NM_by_CIGAR_X));
                if (segconf.is_minimap2 || segconf.SA_CIGAR_abbreviated==yes) bufprintf (evb, &features, "SA_CIGAR_abbrev=%s;",  YN(segconf.SA_CIGAR_abbreviated));
            }
                        
            FEATURE (true, "Aligner: %s", "mapper=%s", segconf_sam_mapper_name()); 

            if (ZCTX(SAM_QUAL)->qual_codec == CODEC_DOMQ && z_file->num_lines)
                bufprintf (evb, &features, "QUAL=DOMQ%s (DIVR:%.f%%);", segconf.sam_has_xcons ? ".xcons" : "",
                           percent (z_file->divr_lines, z_file->num_lines));                
            else
                bufprintf (evb, &features, "QUAL=%s;", !segconf.nontrivial_qual ? "Trivial" : ZCTX(SAM_QUAL)->qual_codec != CODEC_UNKNOWN ? codec_name (ZCTX(SAM_QUAL)->qual_codec) : codec_name (ZCTX(SAM_QUAL)->lcodec));

            bufprintf (evb, &features, "QUAL_histo=%s;", segconf_get_qual_histo(QHT_QUAL).s);
            
            if (segconf_has(OPTION_OQ_Z)) {
                bufprintf (evb, &features, "OQ=%s;", (ZCTX(OPTION_OQ_Z)->qual_codec != CODEC_UNKNOWN) ? codec_name (ZCTX(OPTION_OQ_Z)->qual_codec) : codec_name (ZCTX(OPTION_OQ_Z)->lcodec));
                bufprintf (evb, &features, "OQ_histo=%s;", segconf_get_qual_histo(QHT_OQ).s);
            }

            bufprintf (evb, &features, "smux_max_stdv=%2.1f%% '%s';",
                       segconf.smux_max_stdv * 100, 
                       segconf.smux_max_stdv_q=='='?"≐" : segconf.smux_max_stdv_q==';'?"；" : char_to_printable_json (segconf.smux_max_stdv_q).s); // Unicode ；and ≐ (in UTF-8) to avoid breaking spreadsheet

            FEATURE0 (segconf.sam_bisulfite, "Feature: Bisulfite", "Bisulfite");
            FEATURE0 (segconf.has_10xGen, "Feature: 10xGenomics_tags", "has_10xGen");
            FEATURE0 (segconf.has_Parse, "Feature: Parse_tags", "has_Parse");

            if (segconf.sam_ms_type) {
                rom names[] = ms_type_NAME;
                bufprintf (evb, &features, "ms:i_type=%s;", names[segconf.sam_ms_type]);
            }

            if (segconf.sam_XG_inc_S != unknown) 
                bufprintf (evb, &features, "XG_include_S=%s;", TF(segconf.sam_XG_inc_S));

            if (segconf_has(OPTION_RG_Z))
                bufprintf (evb, &features, "RG_method=%s;", RG_method_name (segconf.RG_method));

            if (num_alignments) {
                double mate_line_pc   = percent (z_file->mate_line_count, num_alignments);
                double saggy_near_pc  = percent (z_file->saggy_near_count, num_alignments);
                double depn_far_pc    = percent (z_file->depn_far_count, num_alignments);
                double secondary_pc   = percent (z_file->secondary_count, num_alignments);
                double supp_pc        = percent (z_file->supplementary_count, num_alignments);
                double far_of_depn_pc = percent (z_file->depn_far_count, z_file->comp_num_lines[SAM_COMP_DEPN]);
                #define PREC(f) (((f) && ((f)<10 || (f)>95)) ? (1 + ((f)<1 || (f)>99.5)) : 0)
                FEATURE(true, "Buddying: sag_type=%s mate=%.*f%% saggy_near=%.*f%% depn_far=%.*f%% depn_far/num_DEPN=%.*f%% sec=%.*f%% supp=%.*f%%",
                        "sag_type=%s;mate=%.*f%%;saggy_near=%.*f%%;depn_far=%.*f%%;depn_far/num_DEPN=%.*f%%;secondary=%.*f%%;suppl=%.*f%%",
                        sag_type_name (segconf.sag_type), PREC(mate_line_pc), mate_line_pc, PREC(saggy_near_pc), saggy_near_pc, PREC(depn_far_pc), depn_far_pc, 
                        PREC(far_of_depn_pc), far_of_depn_pc, PREC(secondary_pc), secondary_pc, PREC(supp_pc), supp_pc);

                if (flag.deep)
                    bufprintf (evb, &stats, "Deep: fq_reads_deeped=%.1f%% qtype=%s no_qual=%s trim=%s%s%s\n", 
                               deep_pc, 
                               qtype_name (segconf.deep_qtype), TF(segconf.deep_no_qual), 
                               segconf_deep_trimming_name(),
                               cond_str (segconf.deep_N_sam_score && segconf.deep_N_fq_score, " N_qual=", char_to_printable_json (segconf.deep_N_fq_score).s),
                               cond_str (segconf.deep_N_sam_score && segconf.deep_N_fq_score, "→", char_to_printable_json (segconf.deep_N_sam_score).s));

                if (flag.deep)
                    FEATURE (z_file->num_lines, "SAM qname: %s%s", "SAM_Qname=%s%s", 
                             segconf_qf_name (QSAM), cond_str(segconf.deep_sam_qname_flavor[1], "+", segconf_qf_name(QSAM2)))

                REPORT_QNAME;
                FEATURE (segconf.tech, "Sequencer: %s", "sequencer=%s", tech_name(segconf.tech));
                
                bufprintf (evb, &features, "tlen_pred=%.1f%%;", percent (z_file->sam_num_tlen_pred, num_alignments));

                if (z_file->sam_num_seq_by_aln) // seg SEQ vs internal or external reference according to SAM alignment 
                    bufprintf (evb, &features, "seq_by_sam_aln=%.1f%%;", percent (z_file->sam_num_seq_by_aln, num_alignments));

                if (z_file->sam_num_by_prim)    // seg SEQ vs PRIM VB or vs saggy line
                    bufprintf (evb, &features, "seq_by_prim=%.1f%%;", percent (z_file->sam_num_by_prim, num_alignments)); 

                if (z_file->sam_num_aligned)    // seg SEQ vs external reference using our aligner
                    bufprintf (evb, &features, "seq_by_aligner (perfect)=%.1f%% (%.1f%%);", percent (z_file->sam_num_aligned, num_alignments), percent (z_file->sam_num_aligned_perfect, num_alignments)); // report even if num_aligned=0 (i.e. wrong reference)           

                if (z_file->sam_num_verbatim)   // seg SEQ by storing verbatim
                    bufprintf (evb, &features, "seq_by_verbatim=%.1f%%;", percent (z_file->sam_num_verbatim, num_alignments));

                if (flag.deep && (z_file->deep_stats[NDP_SAM_DUP] || z_file->deep_stats[NDP_SAM_DUP_TRIM]))
                    bufprintf (evb, &features, "sam_alns_deep_hash_contn=%"PRIu64";", z_file->deep_stats[NDP_SAM_DUP] + z_file->deep_stats[NDP_SAM_DUP_TRIM]);
            }

            if (num_fq_reads) {
                if (flag.deep) 
                    bufprintf (evb, &features, "fq_deep=%.1f%%;", percent (z_file->deep_stats[NDP_DEEPABLE] + z_file->deep_stats[NDP_DEEPABLE_TRIM], num_fq_reads));

                if (z_file->fq_num_aligned)
                    bufprintf (evb, &features, "fq_aligner_ok (perfect|spliced|per&spl)=%.1f%% (%.1f%% | %.1f%% | %.1f%%);", 
                               percent (z_file->fq_num_aligned, num_fq_reads), 
                               percent (z_file->fq_num_aligned_perfect, num_fq_reads),  // this includes all perfect matches - spliced or not
                               percent (z_file->fq_num_aligned_spliced, num_fq_reads),  // this includes all spliced matches - perfect or not
                               percent (z_file->fq_num_perfect_spliced, num_fq_reads)); // report even if num_aligned=0 (i.e. wrong reference)           

                if (z_file->fq_num_verbatim)
                    bufprintf (evb, &features, "fq_verbatim=%.1f%%;", percent (z_file->fq_num_verbatim, num_fq_reads));
            }

            if (segconf.tech != segconf.tech_by_RG) {
                if (segconf.tech_by_RG)
                    bufprintf (evb, &features, "tech_by_RG=%s;", tech_name (segconf.tech_by_RG));
                else if (z_file->tech_by_RG_unidentified[0])
                    bufprintf (evb, &features, "@RG_PL=%s;", z_file->tech_by_RG_unidentified);
            }

            break;
        }

        case DT_FASTQ:
            REPORT_VBs;
            REPORT_QNAME;
            FEATURE (segconf.tech, "Sequencer: %s", "sequencer=%s", tech_name(segconf.tech));
            FEATURE0 (FAF, "FASTA-as-FASTQ", "FASTA-as-FASTQ=True");
            FEATURE0 (segconf.multiseq, "Multiseq", "multiseq=True");
            FEATURE0 (segconf.is_interleaved, "Interleaved", "interleaved=True");
            FEATURE (segconf.nonbio_type, "Nonbio: %s", "nonbio_type=%s", fastq_nonbio_type_name());
            
            double bamass_pc = percent (z_file->deep_stats[NDP_DEEPABLE] + z_file->deep_stats[NDP_DEEPABLE_TRIM], z_file->deep_stats[NDP_FQ_READS]);

            if (flag.bam_assist) {
                bufprintf (evb, &stats, "Bamass: fq_reads_deeped=%.1f%% qtype=%s trim=%s\n", 
                           bamass_pc, qtype_name (segconf.deep_qtype),
                           segconf_deep_trimming_name());

                bufprintf (evb, &features, "bamass_qtype=%s;", qtype_name (segconf.deep_qtype));
                bufprintf (evb, &features, "bamass_trim=%s;", segconf_deep_trimming_name());
                bufprintf (evb, &features, "bamass_cigar_trims=%s;", (rom[])TRIM_IS_NAMES[segconf.bamass_trims]);
            }

            if (IS_REF_LOADED_ZIP) {
                bufprintf (evb, &features, "ref_ncontigs=%u;", ref_get_ctgs()->contigs.len32);
                bufprintf (evb, &features, "ref_nbases=%"PRIu64";", contigs_get_nbases (ref_get_ctgs()));
            }

            if (bamass_pc) 
                bufprintf (evb, &features, "bamass=%.1f%%;", bamass_pc);

            if (z_file->fq_num_aligned) 
                bufprintf (evb, &features, "aligned (perfect|spliced|per&spl)=%.1f%% (%.1f%% | %.1f%% | %.1f%%);",
                           percent (z_file->fq_num_aligned,         z_file->num_lines), 
                           percent (z_file->fq_num_aligned_perfect, z_file->num_lines),  // this includes all perfect matches - spliced or not
                           percent (z_file->fq_num_aligned_spliced, z_file->num_lines),  // this includes all spliced matches - perfect or not
                           percent (z_file->fq_num_perfect_spliced, z_file->num_lines)); // argument ignored if !fq_num_nonbio
            if (z_file->fq_num_monochar) 
                bufprintf (evb, &features, "monochar=%.1f%%;", percent (z_file->fq_num_monochar, z_file->num_lines));

            if (z_file->fq_num_empty_read) 
                bufprintf (evb, &features, "empty_read=%.1f%%;", percent (z_file->fq_num_empty_read, z_file->num_lines));

            if (z_file->fq_num_nonbio) 
                bufprintf (evb, &features, "nonbio=%.1f%%;", percent (z_file->fq_num_nonbio, z_file->num_lines));

            if (z_file->fq_num_verbatim)
                bufprintf (evb, &features, "verbatim=%.1f%%;", percent (z_file->fq_num_verbatim, z_file->num_lines));

            if (!FAF) {
                if (ZCTX(SAM_QUAL)->qual_codec == CODEC_DOMQ && z_file->num_lines)
                    bufprintf (evb, &features, "QUAL=DOMQ (DIVR:%.f%%);", percent (z_file->divr_lines, z_file->num_lines));                
                else
                    bufprintf (evb, &features, "QUAL=%s;", !segconf.nontrivial_qual ? "Trivial" : ZCTX(SAM_QUAL)->qual_codec != CODEC_UNKNOWN ? codec_name (ZCTX(SAM_QUAL)->qual_codec) : codec_name (ZCTX(SAM_QUAL)->lcodec));
                
                bufprintf (evb, &features, "Qual_histo=%s;", segconf_get_qual_histo(QHT_QUAL).s);

                bufprintf (evb, &features, "smux_max_stdv=%2.1f%% '%s';",
                        segconf.smux_max_stdv * 100.0, 
                        segconf.smux_max_stdv_q=='='?"≐" : segconf.smux_max_stdv_q==';'?"；" : char_to_printable_json (segconf.smux_max_stdv_q).s); // Unicode ；and ≐ (in UTF-8) to avoid breaking spreadsheet

                if (segconf.r1_or_r2 && !flag.pair) 
                    bufprintf (evb, &features, "R1_or_R2=R%d;", (segconf.r1_or_r2 == PAIR_R1) ? 1 : 2);
            }
            break;

        case DT_VCF:
        case DT_BCF:
            bufprintf (evb, &stats, "Samples: %u   ", vcf_header_get_num_samples()); //  no newline
            bufprintf (evb, &features, "num_samples=%u;", vcf_header_get_num_samples());
            if (vcf_header_get_num_contigs()) bufprintf (evb, &features, "hdr_contigs=%u (%"PRIu64");", vcf_header_get_num_contigs(), vcf_header_get_nbases());
            if (IS_REF_LOADED_ZIP) bufprintf (evb, &features, "ref_contigs=%u (%"PRIu64");", ref_get_ctgs()->contigs.len32, contigs_get_nbases (ref_get_ctgs()));

            if (z_file->max_ploidy != 2) 
                bufprintf (evb, &features, "ploidy=%u;", z_file->max_ploidy);

            if (z_file->mate_line_count && z_file->num_lines) {
                double pc = percent (z_file->mate_line_count, z_file->num_lines);
                bufprintf (evb, &stats, "Mated: %.*f%%   ", PREC(pc), pc); //  no newline
                bufprintf (evb, &features, "mated=%.*f%%;", PREC(pc), pc);
            }

            if (z_file->vcf_num_samples_copied && z_file->num_lines) {
                double pc = percent (z_file->vcf_num_samples_copied, z_file->num_lines * vcf_header_get_num_samples());
                bufprintf (evb, &features, "samples_copied=%.*f%%;", PREC(pc), pc);
            }
            else
                bufprintf (evb, &features, "samples_copied=%s;", segconf.vcf_sample_copy ? "None" : "Disabled");

            if (segconf.vcf_is_gvcf) 
                bufprint0 (evb, &features, "gvcf=true;");

            bufprintf (evb, &features, "QUAL_method=%s;", VCF_QUAL_method_name (segconf.vcf_QUAL_method));

            bufprintf (evb, &features, "INFO_method=%s;", VCF_INFO_method_name (segconf.vcf_INFO_method));

            if (segconf_has(INFO_DP))
                bufprintf (evb, &features, "INFO_DP_method=%s;", INFO_DP_method_name (segconf.INFO_DP_method));

            if (segconf_has(FORMAT_DP))
                bufprintf (evb, &features, "FMT_DP_method=%s;", FMT_DP_method_name (segconf.FMT_DP_method));

            if (segconf_has(FORMAT_GQ))
                bufprintf (evb, &features, "GQ_method=%s;", FMT_GQ_method_name (segconf.FMT_GQ_method));

            if (segconf_has(FORMAT_PL))
                bufprintf (evb, &features, "PL_method=%s;",  segconf.PL_mux_by_DP==yes ? "DosageXDP" : "Dosage");

            if (segconf.FMT_GP_content)
                bufprintf (evb, &features, "GP_content=%s;",  FMT_GP_content_name (segconf.FMT_GP_content));

            if (segconf.vcf_is_freebayes)
                bufprintf (evb, &features, "RO_AO_method=%s;",  FMT_ROAO_method_name (segconf.FMT_RO_AO_method));
            
            REPORT_VBs;

            break;
        
        case DT_FASTA:
            FEATURE0 (segconf.fasta_seq_type==SQT_AMINO, "Sequence type: Amino acids",      "Amino_acids");
            FEATURE0 (segconf.fasta_seq_type==SQT_NUKE, "Sequence type: Nucleotide bases", "Nucleotide_bases");
            bufprint0 (evb, &features, "FASTA-as-FASTQ=False;");
            FEATURE0 (segconf.multiseq, "Multiseq", "multiseq");
            FEATURE (true, "Sequences: %"PRIu64, "num_sequences=%"PRIu64, z_file->num_sequences);
            REPORT_VBs;
            break;

        case DT_GFF:
            REPORT_VBs;
            FEATURE (true, "GFF version: %d", "GFF_version=%d", segconf.gff_version); 
            FEATURE (segconf.has_embedded_fasta, "FASTA Sequences: %"PRIu64, "num_fasta_sequences=%"PRIu64, z_file->num_sequences);
            break;
            
        case DT_BED:
            if (z_file->num_lines) {
                FEATURE (true, "Columns: %u", "columns=%u", segconf.bed_num_flds);
                FEATURE0 (segconf.is_sorted, "Sorting: Sorted", "Sorted");        
            }
            break;
            
        case DT_REF: 
            bufprintf (evb, &stats, "%s of genome (in-memory): %s\n", digest_name(), digest_display (z_file->digest).s); 
            FEATURE (true, "Contigs: %u (%"PRIu64")", "ref_contigs=%u (%"PRIu64")", ref_contigs_get_num_contigs(), ref_contigs_get_genome_nbases());
            break;

        case DT_GNRIC:
            REPORT_VBs;
            bufprintf (evb, &stats, "Vblocks: %u x %u MB  Sections: %u\n", 
                       z_file->num_vbs, (uint32_t)(segconf.vb_size >> 20), z_file->section_list.len32);
    
            bufprintf (evb, &features, "magic=%s;extension=%s;", generic_get_magic().s, generic_get_ext().s);
            break;

        default:
            REPORT_VBs;
    }
    
    if (!flag.make_reference && z_file->num_lines) {
        bufprintf (evb, &features, "segconf.line_len=%u;", segconf.line_len); 
        if      (segconf.std_seq_lR2) bufprintf (evb, &features, "segconf.std_seq_len=%u+%u;", segconf.std_seq_len, segconf.std_seq_lR2); 
        else if (segconf.std_seq_len) bufprintf (evb, &features, "segconf.std_seq_len=%u;", segconf.std_seq_len); 
    }

    if (stats_programs.len) {
        if (*BLSTc(stats_programs) == ';') stats_programs.len--; // remove final ';'
     
        bufprintf (evb, &stats, "Programs: %.*s\n", STRfb(stats_programs));
    }

    if (flag.optimize)
        bufprintf (evb, &stats, "Fields optimized: %s\n", segconf_get_optimizations().s);

    bufprintf (evb, &stats, "Genozip version: %s %s\nDate compressed: %s\n", 
               STRver(code_version()).s, get_distribution(), str_time().s);

    bufprint0 (evb, &stats, "Command line: ");
    buf_append_string (evb, &stats, flags_command_line()); // careful not to use bufprintf with command_line as it can exceed the maximum length in bufprintf

    bufprintf (evb, &stats, "\n%s\n", license_get_one_line());
    
    *BAFTc(features) = '\0';
}

static DESCENDING_SORTER (stats_sort_by_z_size, StatsByLine, z_size)

static void stats_consolidate_ctxs (StatsByLine *sbl, unsigned num_stats)
{
    for (unsigned parent=0; parent < num_stats; parent++) 

        // case: we might consolidate to this context
        if (sbl[parent].my_did_i != DID_NONE && 
            (sbl[parent].st_did_i == DID_NONE || sbl[parent].st_did_i == sbl[parent].my_did_i))  {
            
            for (unsigned child=0; child < num_stats; child++) {

                if (child != parent && sbl[parent].my_did_i  == sbl[child].st_did_i) {
                    sbl[parent].txt_len   += sbl[child].txt_len;
                    sbl[parent].z_size    += sbl[child].z_size;  
                    sbl[parent].pc_of_txt += sbl[child].pc_of_txt;
                    sbl[parent].pc_of_z   += sbl[child].pc_of_z;  

                    if (flag.debug_stats)
                        iprintf ("Consolidated %s (did=%u) (txt_len=%"PRIu64" z_size=%"PRIu64") into %s (did=%u) (AFTER: txt_len=%"PRIu64" z_size=%"PRIu64")\n",
                                 sbl[child].name, sbl[child].my_did_i, sbl[child].txt_len, sbl[child].z_size, 
                                 sbl[parent].name, sbl[parent].my_did_i, sbl[parent].txt_len, sbl[parent].z_size);

                    sbl[child] = (StatsByLine){ .my_did_i = DID_NONE, .st_did_i = DID_NONE };
                }
            }

            if (!strcmp (sbl[parent].name, "SQBITMAP")) strcpy (sbl[parent].name, "SEQ"); // rename
        }
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
    for (StatsByLine *s=sbl, *after=sbl + num_stats; s < after; s++) 
        if (s->name[0] && !strcmp (consolidated_name, s->name)) {
            survivor = s;
            break;
        }

    for (StatsByLine *s=sbl, *after=sbl + num_stats; s < after; s++) {
        
        if (!s->name[0]) continue; // unused entry

        for (unsigned d=0; d < num_deps; d++) 
            if (!strcmp (deps[d], s->name)) {
                if (!survivor) {
                    survivor = s; // first found gets to be the survivor
                }
                else {
                    survivor->txt_len   += s->txt_len;
                    survivor->z_size    += s->z_size;  
                    survivor->pc_of_txt += s->pc_of_txt;
                    survivor->pc_of_z   += s->pc_of_z; 
                    strcpy (survivor->name, consolidated_name); // rename only if at least one was consolidated

                    if (flag.debug_stats)
                        iprintf ("Consolidated %s (txt_len=%"PRIu64" z_size=%"PRIu64") into %s (AFTER: txt_len=%"PRIu64" z_size=%"PRIu64")\n",
                                 s->name, s->txt_len, s->z_size, survivor->name, survivor->txt_len, survivor->z_size);

                    *s = (StatsByLine){}; 
                }
                break;
            } 
    }

    va_end (args);
}

static void stats_output_stats (StatsByLine *s, unsigned num_stats, float src_comp_ratio,
                                int64_t all_txt_len, int64_t all_txt_len_0, int64_t all_z_size, float all_pc_of_txt, float all_pc_of_z, float all_comp_ratio)
{
    bufprintf (evb, &stats, "\nSections (sorted by %% of genozip file):%s\n", "");
    bufprintf (evb, &stats, "%s", SHORT_HEADER);

    // if -W or -w appear multiple times, we print the stats header to stderr to allow use with "| grep"
    if (flag.show_stats == STATS_SHORT_GREP) { 
        fprintf (stderr, "%s", SHORT_HEADER);
        fflush (stderr);
    }

    for (uint32_t i=0; i < num_stats; i++, s++)
        if (s->z_size || (s->my_did_i != DID_NONE && segconf_optimize (s->my_did_i) && s->txt_len))
            bufprintf (evb, &stats, "%-20.20s %9s %5.1f%% %9s %5.1f%% %6.*fX\n", 
                       s->name, 
                       str_size (s->z_size).s, s->pc_of_z, // z size and % of total z that is in this line
                       str_size ((float)s->txt_len).s, s->pc_of_txt, // txt size and % of total txt which is in this line
                       ((float)s->txt_len / (float)s->z_size) < 1000 ? 1 : 0, // precision
                       (float)s->txt_len / (float)s->z_size); // ratio z vs txt

    if (src_comp_ratio != 1 && flag.show_stats_comp_i == COMP_NONE) {

        // display codec name if same source codec for all components, or "DISK_SIZE" if not
        rom source_code_name = codec_name (z_file->comp_src_codec[0]);

        for (CompIType comp_i=1; comp_i < z_file->num_components; comp_i++)
            if (z_file->comp_src_codec[comp_i] != z_file->comp_src_codec[0])
                source_code_name = "DISK_SIZE";

        bufprintf (evb, &stats, 
                   "GENOZIP vs %-9s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                   source_code_name,
                   str_size (all_z_size).s, all_pc_of_z, // total z size and sum of all % of z (should be 100)
                   str_size (all_txt_len_0 / src_comp_ratio).s, all_pc_of_txt, // total txt fize and ratio z vs txt
                   all_comp_ratio / src_comp_ratio);
    }

    bufprintf (evb, &stats, 
                "%-20s %9s %5.1f%% %9s %5.1f%% %6.1fX\n", 
                src_comp_ratio != 1 ? "GENOZIP vs TXT" : "TOTAL",
                str_size (all_z_size).s, all_pc_of_z, // total z size and sum of all % of z (should be 100)
                str_size (all_txt_len_0).s, all_pc_of_txt, // total txt fize and ratio z vs txt
                all_comp_ratio);
}

static void stats_output_STATS (StatsByLine *s, unsigned num_stats,
                                int64_t all_txt_len, int64_t all_uncomp_dict, int64_t all_comp_dict, int64_t all_comp_b250, int64_t all_comp_local, 
                                int64_t all_z_size, float all_pc_of_txt, float all_pc_of_z, float all_comp_ratio)
{
#define PC(pc) ((pc==0 || pc>=10) ? 0 : (pc<1 ? 2:1))

    // add diagnostic info
    if (flag.show_stats_comp_i == COMP_NONE) {
        bufprintf (evb, &STATS, "\nSystem info: OS=%s cores=%u RAM=%1.1f GB\n", 
                   arch_get_os(), arch_get_num_cores(), arch_get_physical_mem_size());
        bufprintf (evb, &STATS, "\nSections (sorted by %% of genozip file):%s\n", "");
    }

    bufprintf (evb, &STATS, "%s", LONG_HEADER);
    // note: bufprint0 requires % to display a %, while fprintf requires %%.
    
    // if -W or -w appear multiple times, we print the stats header to stderr to allow use with "| grep"
    if (flag.show_stats == STATS_LONG_GREP) 
        fprintf (stderr, "%s", LONG_HEADER);

    for (uint32_t i=0; i < num_stats; i++, s++)
        if (s->z_size || (s->my_did_i != DID_NONE && segconf_optimize (s->my_did_i) && s->txt_len))
            bufprintf (evb, &STATS, "%-2.2s    %-17.17s %-17.17s %6s %6s %6s %6s %8s %3.0f%% %3.0f%% %9s %9s %9s %9s %9s %9s %9s %6.*fX %5.1f%% %5.1f%%\n", 
                       s->did_i.s, s->name, s->type, s->words.s, 
                       s->dict_words.s, s->local_words.s, s->failed_ston_words.s, 
                       (s->is_dict_alias ? "ALIAS  " : s->hash.s), s->pc_hash_occupancy, s->pc_ston_hash_occup, // Up to here - these don't appear in the total
                       s->uncomp_dict.s, s->comp_dict.s, s->comp_b250.s, s->uncomp_local.s, s->comp_local.s, str_size (s->z_size).s, 
                       str_size ((float)s->txt_len).s, 
                       ((float)s->txt_len / (float)s->z_size) < 1000 ? 1 : 0, // precision
                        (float)s->txt_len / (float)s->z_size, // ratio z vs txt
                       s->pc_of_txt, s->pc_of_z);

    // note: no point showing this per component in SAM and DVCF, bc all_txt_len_0 is fully accounted for in the MAIN component and it is 0 in the others
    if (!(z_has_gencomp && flag.show_stats_comp_i != COMP_NONE))
        bufprintf (evb, &STATS, "TOTAL                                                                               "
                "%9s %9s %9s %9s %9s %9s %6.1fX %5.1f%% %5.1f%%\n", 
                str_size (all_uncomp_dict).s, str_size (all_comp_dict).s, str_size (all_comp_b250).s, 
                str_size (all_comp_local).s,  str_size (all_z_size).s,    str_size (all_txt_len).s, 
                all_comp_ratio, all_pc_of_txt, all_pc_of_z);
}

// generate the stats text - all sections except genozip header and the two stats sections 
void stats_generate (void) // specific section, or COMP_NONE if for the entire file
{
    // initial allocation
    buf_alloc (evb, &stats,      0, 8 KB, char, 1, "stats");
    buf_alloc (evb, &STATS,      0, 8 KB, char, 1, "stats");
    buf_alloc (evb, &features,   0, 1 KB, char, 1, "stats");
    buf_alloc (evb, &exceptions, 0, 1 KB, char, 1, "stats");
    buf_alloc (evb, &internals,  0, 8 KB, char, 1, "stats");

    if (flag.show_stats_comp_i == COMP_NONE) {
        stats_output_file_metadata();
        buf_copy (evb, &STATS, &stats, char,0,0, "stats");
    }

    int64_t all_comp_dict=0, all_uncomp_dict=0, all_comp_b250=0, all_comp_local=0, all_z_size=0;

    // prepare data
    #define NUM_SBL (NUM_SEC_TYPES + z_file->ca.num_contexts + 2) // 2 for consolidated groups
    ARRAY_alloc (StatsByLine, sbl, NUM_SBL, true, sbl_buf, evb, "stats");

    #define ST_NAME(st) (&st_name(st)[4]) // cut off "SEC_" 

    stats_get_compressed_sizes (sbl); // ctx sizes in ctx->local/dict/ctx.param, and other sections in sbl[st]->z_size.
    
    uint64_t z_size     = flag.show_stats_comp_i==COMP_NONE ? z_file->disk_so_far            : z_file->disk_so_far_comp           [flag.show_stats_comp_i];
    uint64_t txt_size   = flag.show_stats_comp_i==COMP_NONE ? z_file->txt_data_so_far_bind   : z_file->txt_data_so_far_bind_comp  [flag.show_stats_comp_i];
    uint64_t txt_size_0 = flag.show_stats_comp_i==COMP_NONE ? z_file->txt_data_so_far_bind_0 : z_file->txt_data_so_far_bind_0_comp[flag.show_stats_comp_i];
    StatsByLine *s = sbl;

    for (SectionType st=0; st < NUM_SEC_TYPES; st++, s++) { 

        if (ST(DICT) || ST(B250) || ST(LOCAL) || ST(COUNTS)) continue; // these are covered by individual contexts

        s->txt_len    = ST(TXT_HEADER) ? z_file->header_size : 0; 
        s->type       = (ST(REFERENCE) || ST(REF_IS_SET) || ST(REF_CONTIGS) || ST(CHROM2REF_MAP) || ST(REF_IUPACS)) ? "SEQUENCE" 
                      : (ST(RANDOM_ACCESS) || ST(REF_RAND_ACC))                                                     ? "RandomAccessIndex"
                      :                                                                                               "Other"; // note: some contexts appear as "Other" in --stats, but in --STATS their parent is themself, not "Other"
        s->my_did_i   = s->st_did_i = DID_NONE;
        s->did_i.s[0] = s->words.s[0] = s->hash.s[0] = s->uncomp_dict.s[0] = s->comp_dict.s[0] = '-';
        s->pc_of_txt  = percent (s->txt_len, txt_size);
        s->pc_of_z    = percent (s->z_size, z_size);
        strcpy (s->name, ST_NAME (st)); 

        all_z_size     += s->z_size;
        all_comp_local += s->z_size;
        all_txt_len    += s->txt_len;
    }

    z_file->header_size = 0; // reset (in case of showing components)
    
    // contexts
    for_zctx {    
        s->z_size = zctx->dict.count + zctx->b250.count + zctx->local.count;

        if (!zctx->b250.count && !zctx->txt_len && !zctx->b250.len && !zctx->is_stats_parent && !zctx->txt_shrinkage && 
            !s->z_size) 
            continue;

        s->txt_len = zctx->txt_len + zctx->txt_shrinkage; // original txt_len before modifications
        
        all_comp_dict   += zctx->dict.count;
        all_uncomp_dict += zctx->dict.len;
        all_comp_b250   += zctx->b250.count;
        all_comp_local  += zctx->local.count;
        all_txt_len     += zctx->txt_len; // sum of txt after modifications
        all_z_size      += s->z_size;

        if ((Z_DT(VCF) || Z_DT(GFF)) && dict_id_type(zctx->dict_id) != DTYPE_FIELD)
            snprintf (s->name, sizeof (s->name), "%s/%.*s", dtype_name_z(zctx->dict_id), MAX_TAG_LEN, zctx->tag_name);
        else  
            strcpy (s->name, zctx->tag_name);

        // parent
        s->type  = zctx->dict_id.num                 == _SAM_SQBITMAP ? "SEQUENCE"
                 : zctx->st_did_i == DID_NONE                         ? zctx->tag_name
                 : ZCTX(zctx->st_did_i)->dict_id.num == _SAM_SQBITMAP ? "SEQUENCE"
                 :                                                      ZCTX(zctx->st_did_i)->tag_name;

        // note: each VB contributes local.len contains its b250.count if it has it, and local_num_len if not 
        float n_words = zctx->word_list.count; // fields segged of this type in the file : set in ctx_update_stats()

        s->my_did_i           = zctx->did_i;
        s->st_did_i           = zctx->st_did_i;
        s->did_i              = str_int_commas ((uint64_t)zctx->did_i); 
        s->words              = str_uint_commas_limit (n_words, 99999);
        s->dict_words         = str_uint_commas_limit (MIN_(zctx->nodes.len, n_words), 99999); // MIN_ is a workaround - not sure why nodes.len sometimes exceeds the dictionary words on the file (eg in TOPLEVEL)
        s->local_words        = str_uint_commas_limit (zctx->local_num_words, 99999);
        s->failed_ston_words  = str_uint_commas_limit (zctx->num_failed_singletons, 99999);
        s->pc_hash_occupancy  = percent (zctx->nodes.len, zctx->global_hash.len32);
        s->pc_ston_hash_occup = percent (zctx->ston_ents.len, zctx->global_hash.len32);
        s->hash               = str_size (zctx->global_hash.len32);
        s->uncomp_dict        = str_size (zctx->dict.len);
        s->comp_dict          = str_size (zctx->dict.count);
        s->comp_b250          = str_size (zctx->b250.count);
        s->uncomp_local       = str_size (zctx->local.len); // set in ctx_update_stats()
        s->comp_local         = str_size (zctx->local.count);
        s->pc_of_txt          = percent (s->txt_len, txt_size);
        s->pc_of_z            = percent (s->z_size, z_size);
        s->dcodec             = codec_name (zctx->dcodec);
        s->bcodec             = codec_name (zctx->bcodec);
        s->lcodec             = codec_name (zctx->lcodec_non_inherited ? zctx->lcodec_non_inherited : zctx->lcodec);
        s->is_dict_alias      = (zctx->did_i != zctx->dict_did_i);
        
        s++; // increment only if it has some data, otherwise already continued
    }

    sbl_buf.len32 = s - sbl;

    // note: for txt size and compression ratio in the TOTAL line (all_comp_ratio) we use txt_data_so_far_bind_0 
    // (the original txt data size) and not all_txt_size (size after ZIP modifications like --optimize). 
    // Therefore, in case of ZIP-modified txt, the sum of the (modified) fields in the TXT column will NOT equal the
    // TOTAL in the TXT column. That's ok.
    all_comp_ratio      = (float)txt_size_0 /* without modifications */ / (float)all_z_size;
    float all_pc_of_txt = percent (all_txt_len, txt_size/* with modifications */);
    float all_pc_of_z   = percent (all_z_size, z_size);

    // long form stats from --STATS    
    qsort (sbl, sbl_buf.len32, sizeof (sbl[0]), stats_sort_by_z_size);  // sort by compressed size

    stats_prepare_internals (sbl, sbl_buf.len32, all_z_size);

    stats_output_STATS (sbl, sbl_buf.len32, 
                        all_txt_len, all_uncomp_dict, all_comp_dict, all_comp_b250, all_comp_local, all_z_size, all_pc_of_txt, all_pc_of_z, all_comp_ratio);

    // calculate hash occupancy before consolidating stats
    stats_calc_hash_occ (sbl, sbl_buf.len32);
    
    // consolidates stats of child contexts into the parent one
    stats_consolidate_ctxs (sbl, sbl_buf.len32);
    
    stats_consolidate_non_ctx (sbl, sbl_buf.len32, 
                               IS_REF_EXT_STORE ? "Reference" : "SEQ", // when compressing SAM/FASTQ with REF_EXT_STORE, account for the reference in its own "Parent"
                               5, ST_NAME (SEC_REFERENCE), ST_NAME (SEC_REF_IS_SET), 
                               ST_NAME (SEC_REF_CONTIGS), ST_NAME (SEC_CHROM2REF_MAP),
                               ST_NAME (SEC_REF_IUPACS));

    stats_consolidate_non_ctx (sbl, sbl_buf.len32, "Other", 20 + (DTPZ(txt_header_required) == HDR_NONE), "E1L", "E2L", "EOL", 
                               "SAMPLES", "AUX", TOPLEVEL, "TOP2BAM", "TOP2NONE", "TOP2VCF", "BAM_BIN", "LINEMETA", "CONTIG", "SAG", "SAALN",
                               ST_NAME (SEC_DICT_ID_ALIASES), ST_NAME (SEC_RECON_PLAN), ST_NAME (SEC_GENCOMP),
                               ST_NAME (SEC_VB_HEADER), ST_NAME (SEC_GZ_ISIZES), ST_NAME (SEC_GZ_DIGESTS), ST_NAME(SEC_TXT_HEADER)/*must be last*/);
        
    stats_consolidate_non_ctx (sbl, sbl_buf.len32, "RandomAccessIndex", 2, ST_NAME (SEC_RANDOM_ACCESS), ST_NAME (SEC_REF_RAND_ACC));
    
    ASSERTW (all_txt_len == txt_size || flag.make_reference, // all_txt_len=0 in make-ref as there are no contexts
             "Expecting all_txt_len=Σ(ctx.txt_len)=%"PRId64" == txt_size%s=%"PRId64" (diff=%"PRId64")", 
             all_txt_len, (segconf.zip_txt_modified ? "(as modified)" : ""), txt_size, (int64_t)txt_size - all_txt_len);

    // consolidated stats for --stats    
    qsort (sbl, sbl_buf.len32, sizeof (sbl[0]), stats_sort_by_z_size);  // re-sort after consolidation

    // source compression, eg BGZF, against txt before any modifications
    src_comp_ratio = (float)txt_size_0 / z_file->txt_file_disk_sizes_sum; 

    stats_output_stats (sbl, sbl_buf.len32, src_comp_ratio, all_txt_len, txt_size_0, all_z_size, all_pc_of_txt, all_pc_of_z, all_comp_ratio);
    
    // if we're showing stats of a single components - output it now
    if (flag.show_stats_comp_i != COMP_NONE) {
        iprintf ("\n\nComponent=%s:\n", comp_name (flag.show_stats_comp_i));
        buf_print (ABS((int)flag.show_stats) == STATS_SHORT ? &stats : &STATS, false);
        buf_free (stats);
        buf_free (STATS);
    }

    // note: we use txt_data_so_far_bind is the sum of recon_sizes - see zip_update_txt_counters - which is
    // expected to be the sum of txt_len. However, this NOT the size of the original file which is stored in
    // z_file->txt_data_so_far_bind_0.
    ASSERT (!flag.debug_or_test || flag.show_stats_comp_i != COMP_NONE || all_txt_len == txt_size || segconf.zip_txt_modified || flag.make_reference, 
            "Hmm... incorrect calculation for %s sizes: total section sizes=%s but file size is %s (diff=%d). "
            _TIP"Use --debug-recon-size for debugging this", 
            z_dt_name(), str_int_commas (all_txt_len).s, str_int_commas (txt_size).s, 
            (int32_t)(txt_size - all_txt_len)); 

    if (flag.show_stats_comp_i == COMP_NONE) { 
        zfile_compress_section_data (evb, SEC_STATS, &stats);
        zfile_compress_section_data (evb, SEC_STATS, &STATS);

        // store stats overhead, for stats_display (because when stats_display runs, section list won't be accessible since it is already converted to SectionEntFileFormat)
        Section sec = sections_last_sec (SEC_NONE, HARD_FAIL) - 1; // first stats section     
        stats.count = z_file->disk_so_far - sec->offset;
    }

    if (flag.show_stats_comp_i != COMP_NONE) 
        iprint0 ("\nNote: Components stats don't include global sections like SEC_DICT, SEC_REFERENCE etc\n");

    if (!(flag.show_stats && z_file->z_closes_after_me)) { // case: we're not going to display these stats - free now
        buf_destroy (stats);
        buf_destroy (STATS);
    }
}

static void stats_display (void)
{
    BufferP buf = ABS((int)flag.show_stats) == 1 ? &stats : &STATS;

    if (!buf_is_alloc (buf)) return;  // no stats available

    buf_print (buf , false);

    if (stats.count && z_file->disk_so_far < 1 MB && IS_ZIP)  // no need to print this note if z size > 1MB, as the 2-4KB of overhead is rounded off anyway
        // stats text doesn't include SEC_STATS and SEC_GENOZIP_HEADER - the last 3 sections in the file - since stats text is generated before these sections are compressed
        iprintf ("\nNote: ZIP total file size excludes overhead of %s\n", str_size (stats.count).s);

    iprint_newline();

    buf_destroy (stats);
    buf_destroy (STATS);
}

void stats_read_and_display (void)
{
    Section sec = sections_last_sec (SEC_STATS, SOFT_FAIL);
    if (!sec) {
        iprint0 ("No stats available for this file.\n"); // possibly the file was compressed with --cstats or --CSTATS
        return; // genozip file does not contain stats sections (SEC_STATS was introduced in v 7.0.5)
    }

    // if -W or -w appear multiple times, we print the stats header to stderr to allow use with "| grep"
    if (flag.show_stats == STATS_SHORT_GREP) 
        fprintf (stderr, "%s", SHORT_HEADER);

    else if (flag.show_stats == STATS_LONG_GREP) 
        fprintf (stderr, "%s", LONG_HEADER);

    // read and uncompress the requested stats section
    zfile_get_global_section (SectionHeader, sec - (ABS((int)flag.show_stats)==1),
                              ABS((int)flag.show_stats) == 1 ? &stats : &STATS, "stats");
    
    stats_display();
}

// concatenate txt names of bound files so we can show them all
void stats_add_txt_name (rom fn)
{
    bufprintf (evb, &z_file->bound_txt_names, "%s%s", z_file->bound_txt_names.len ? " ": "", fn);
    z_file->bound_txt_names.count++; // number of files
}

void stats_finalize (void)
{
    if (flag.show_stats) stats_display();

    buf_destroy (stats);
    buf_destroy (STATS);
    buf_destroy (features);
    buf_destroy (exceptions);
    buf_destroy (internals);
    buf_destroy (STATS);
    buf_destroy (sbl_buf);
    buf_destroy (stats_programs);

    all_txt_len = 0;
    src_comp_ratio = all_comp_ratio = 0;    
}

// called if --show-seg-summary
void stats_show_seg_summary (void)
{
    ASSERTNOTNULL (z_file);

    uint64_t num_alignments=0, num_fq_reads=0, num_samples=0;
    
    if (Z_DT(SAM) || Z_DT(BAM)) {
        num_alignments = z_file->comp_num_lines[SAM_COMP_MAIN] + z_file->comp_num_lines[SAM_COMP_PRIM] + z_file->comp_num_lines[SAM_COMP_DEPN]; // excluding Deep FQ components

        for (CompIType comp_i=SAM_COMP_FQ00; comp_i < z_file->num_components; comp_i++)
            num_fq_reads += z_file->comp_num_lines[comp_i];
    }

    else if (Z_DT(FASTQ))
        num_fq_reads = z_file->num_lines;

    else if (Z_DT(VCF))
        num_samples = z_file->num_lines * vcf_header_get_num_samples();

    iprint0 ("\nSeg summary:\n");

    #define SUMMARIZE(x, total, level)        if (z_file->x && (total)) iprintf ("%.*s" #x "=%"PRIu64" (%1.1f%%)\n", level*2, "          ",         z_file->x, percent (z_file->x, total));
    #define SUMMARIZE_(x, total, level, name) if (z_file->x && (total)) iprintf ("%.*s%s=%"PRIu64" (%1.1f%%)\n",    level*2, "          ", (name), z_file->x, percent (z_file->x, total));
    if (num_alignments) iprintf ("sam_num_alignments=%"PRIu64" (100%%)\n", num_alignments);
    SUMMARIZE(sam_num_seq_by_aln,       num_alignments, 1);
    SUMMARIZE(sam_num_aligned,          num_alignments, 1);
    SUMMARIZE(sam_num_aligned_perfect,  num_alignments, 2);
    SUMMARIZE(sam_num_verbatim,         num_alignments, 1);
    SUMMARIZE(sam_num_by_prim,          num_alignments, 1);
    SUMMARIZE(sam_num_tlen_pred,        num_alignments, 1);
    
    if (num_fq_reads) iprintf ("fq_num_reads=%"PRIu64" (100%%)\n", num_fq_reads);
    SUMMARIZE_(deep_stats[NDP_DEEPABLE],      num_fq_reads, 1, (flag.bam_assist ? "fq_num_bamass" : "fq_num_deep"));
    SUMMARIZE_(deep_stats[NDP_DEEPABLE_TRIM], num_fq_reads, 1, (flag.bam_assist ? "fq_num_bamass_trim" : "fq_num_deep_trim"));
    SUMMARIZE(fq_num_nonbio,            num_fq_reads, 1);
    SUMMARIZE(fq_num_aligned,           num_fq_reads, 1);
    SUMMARIZE(fq_num_aligned_perfect,   num_fq_reads, 2);
    SUMMARIZE(fq_num_perfect_r2,        num_fq_reads, 3);
    SUMMARIZE(fq_num_r2_gpos_Δ,         num_fq_reads, 2);
    SUMMARIZE(fq_num_aligned_spliced,   num_fq_reads, 2);
    SUMMARIZE(fq_num_perfect_spliced,   num_fq_reads, 3);
    SUMMARIZE(fq_num_monochar,          num_fq_reads, 1);
    SUMMARIZE(fq_num_verbatim,          num_fq_reads, 1);
    SUMMARIZE(fq_num_empty_read,        num_fq_reads, 1);
    
    SUMMARIZE(vcf_num_samples_copied,   num_samples,  1);

    // TO DO: add deep, bamass summaries
}
