// ------------------------------------------------------------------
//   sam_piz.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "sam_private.h"
#include "regions.h"
#include "aligner.h"
#include "coverage.h"
#include "bases_filter.h"
#include "lookback.h"
#include "qname.h"
#include "writer.h"
#include "codec.h"
#include "qname_filter.h"
#include "libdeflate_1.19/libdeflate.h"

static void seq_filter_destroy (void); // forward

void sam_piz_xtra_line_data (VBlockP vb_)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    if (SAM_PIZ_HAS_SAG) {
        iprintf ("grp_i=%u", ZGRP_I(vb->sag));
        if (IS_SAG_SA) iprintf (" aln_i=%"PRIu64" (%d within group)", ZALN_I(vb->sa_aln), (int)(ZALN_I(vb->sa_aln) - vb->sag->first_aln_i));
        iprint0 ("\n");
        
        sam_show_sag_one_grp (ZGRP_I(VB_SAM->sag));
    }
}

void sam_piz_genozip_header (ConstSectionHeaderGenozipHeaderP header)
{
    if (VER(14)) {
        segconf.sam_seq_len           = BGEN32 (header->sam.segconf_seq_len); 
        segconf.seq_len_to_cm         = header->sam.segconf_seq_len_cm;
        segconf.sam_ms_type           = header->sam.segconf_ms_type;
        segconf.has_MD_or_NM          = header->sam.segconf_has_MD_or_NM;
        segconf.sam_bisulfite         = header->sam.segconf_bisulfite;
        segconf.is_paired             = header->sam.segconf_is_paired;
        segconf.sag_type              = header->sam.segconf_sag_type;
        segconf.sag_has_AS            = header->sam.segconf_sag_has_AS;
        segconf.pysam_qual            = header->sam.segconf_pysam_qual;
        segconf.has_10xGen            = header->sam.segconf_10xGen;
        segconf.SA_HtoS               = header->sam.segconf_SA_HtoS;
        segconf.is_sorted             = header->sam.segconf_is_sorted;
        segconf.is_collated           = header->sam.segconf_is_collated;
        segconf.seq_len_dict_id       = header->sam.segconf_seq_len_dict_id; 
        segconf.MD_NM_by_unconverted  = header->sam.segconf_MD_NM_by_un;
        segconf.sam_predict_meth_call = header->sam.segconf_predict_meth;
    }

    if (VER(15)) {
        segconf.deep_qtype            = header->sam.segconf_deep_qname1 ? QNAME1
                                      : header->sam.segconf_deep_qname2 ? QNAME2
                                      :                                   QNONE;
        segconf.deep_no_qual          = header->sam.segconf_deep_no_qual;
        segconf.deep_N_fq_score       = header->sam.segconf_deep_N_fq_score;
        segconf.use_insertion_ctxs    = header->sam.segconf_use_ins_ctxs;
        segconf.est_sam_factor        = (double)header->sam.segconf_sam_factor / (double)SAM_FACTOR_MULT;
        segconf.MAPQ_use_xq           = header->sam.segconf_MAPQ_use_xq;
        segconf.has[OPTION_MQ_i]      = header->sam.segconf_has_MQ;
    }
}

// main thread: called for each txt file, after reading global area, before reading txt header
bool sam_piz_initialize (CompIType comp_i)
{
    if (flag.deep && (comp_i == SAM_COMP_MAIN || comp_i == SAM_COMP_NONE/*single component or interleave*/)) 
        sam_piz_deep_initialize();

    return true;
}

void sam_piz_header_init (void)
{
    if (flag.qnames_file)
        qname_filter_initialize_from_file (flag.qnames_file); // note: uses SectionHeaderTxtHeader.flav_props to canonize qnames

    if (flag.qnames_opt)
        qname_filter_initialize_from_opt (flag.qnames_opt); 
}

bool sam_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header)
{
    if (IS_PRIM(vb))
        VB_SAM->plsg_i = sam_piz_get_plsg_i (vb->vblock_i);

    if (flag.deep)
        sam_piz_deep_init_vb (VB_SAM, header);

    return true; // all good
}

// PIZ compute thread: called after uncompressing contexts and before reconstructing
void sam_piz_recon_init (VBlockP vb)
{
    // we can proceed with reconstructing a PRIM or DEPN vb, only after SA Groups is loaded - busy-wait for it
    while ((IS_PRIM(vb) || IS_DEPN(vb)) && !sam_is_sag_loaded())
        usleep (250000); // 250 ms

    buf_alloc_zero (vb, &CTX(SAM_CIGAR)->cigar_anal_history, 0, vb->lines.len, CigarAnalItem, 0, "cigar_anal_history"); // initialize to exactly one per line.
}

// PIZ: piz_after_recon callback: called by the compute thread from piz_reconstruct_one_vb. order of VBs is arbitrary
void sam_piz_after_recon (VBlockP vb)
{
}

// PIZ: main thread: piz_process_recon callback: usually called in order of VBs, but out-of-order if --test with no writer
void sam_piz_process_recon (VBlockP vb)
{
    if (flag.collect_coverage)
        coverage_add_one_vb (vb);
    
    if (flag.deep)
        sam_piz_deep_grab_deep_ents (VB_SAM);
}

// PIZ: called by main thread after each SAM/BAM txt file reconstruction is done. Note: not called if
// txt_file is filtered out with --fastq, --R1 etc
void sam_piz_finalize (bool is_last_z_file)
{
    sam_header_finalize();

    // if next, we're going to reconstruct deep FASTQ components we can free some memory now
    if (flag.deep) {
        buf_destroy (z_file->sag_grps);
        buf_destroy (z_file->sag_grps_index);
        buf_destroy (z_file->sag_alns);
        buf_destroy (z_file->sag_qnames);
        buf_destroy (z_file->sag_depn_index);
        buf_destroy (z_file->sag_cigars);
        buf_destroy (z_file->sag_seq);
        buf_destroy (z_file->sag_qual);

        vb_dehoard_memory (true);
    }

    if (is_last_z_file && flag.debug_valgrind)
        seq_filter_destroy();
}

// returns true if section is to be skipped reading / uncompressing
IS_SKIP (sam_piz_is_skip_section)
{
    #define SKIP return true
    #define KEEP return false
    #define SKIPIF(cond)  ({ if (cond) SKIP; })
    #define KEEPIF(cond)  ({ if (cond) KEEP; })
    #define SKIPIFF(cond) ({ if (cond) SKIP; else KEEP; })
    #define KEEPIFF(cond) ({ if (cond) KEEP; else SKIP;})

    // case FASTQ components in a Deep file: use FASTQ skip fucntion instead
    if (comp_i >= SAM_COMP_FQ00 && comp_i != COMP_NONE)
        return fastq_piz_is_skip_section (st, comp_i, dict_id, f8, purpose);

    SectionFlags f = { .flags = f8 };
    bool is_dict = (ST(DICT));
    
    bool preproc = (purpose == SKIP_PURPOSE_PREPROC) && (ST(B250) || ST(LOCAL)); // when loading SA, don't skip the needed B250/LOCAL
    bool dict_needed_for_preproc = (is_dict && z_file->z_flags.has_gencomp);  // when loading a SEC_DICT in a file that has gencomp, don't skip dicts needed for loading SA

    // select representative dict_id

    // case: QNAME subfields: treat as _SAM_Q1NAME
    if (dict_id.num == 0) {}

    else if (dict_id_is_qname_sf(dict_id)) dict_id.num = _SAM_Q1NAME; 

    // case: a mux channel of an OPTION - consider its parent instead. Eg. "M0C:Z0", "M1C:Z1" -> "MC:Z"
    else if (dict_id.id[3]==':' && dict_id.id[1] == dict_id.id[5] && !dict_id.id[6])
        dict_id = (DictId){ .id = { dict_id.id[0], dict_id.id[2], ':', dict_id.id[4] } };

    // case: a mux channel or sub-context of a SAM field - consider its parent instead
    else if (dict_id_is_field (dict_id)) {
        if (dict_id.id[3] == '-') // eg Q03-QUAL (PACB) or CCC-QUAL (SMUX)
            dict_id.num = _SAM_QUAL;

        else if (IS_DIGIT (dict_id.id[1]) && 
                 dict_id.id[3] != '_') // exclude eg C0Y_DOMQ : subfields of QT/CY/BZ/OQ/2Y/QX/TQ are incorrectly of type FIELD, but kept that way for backcomp
            switch (dict_id.id[2]) {
                case 'A': dict_id.num = _SAM_MAPQ  ; break; // eg M0APQ0
                case 'N': dict_id.num = _SAM_PNEXT ; break; // eg P0NEXT
                case 'L': dict_id.num = _SAM_FLAG  ; break; // eg F0LAG0
                case 'O': dict_id.num = _SAM_POS   ; break; // eg P0OS0
                default: ABORT ("unrecognized FIELD dict_id %s", dis_dict_id(dict_id).s);
            }
    }

    bool is_prim = (comp_i == SAM_COMP_PRIM) && !preproc;
    bool is_main = (comp_i == SAM_COMP_MAIN);
    bool cov  = flag.collect_coverage;
    bool cnt  = flag.count && !flag.grep; // we skip if we're only counting, but not also if based on --count, but not if --grep, because grepping requires full reconstruction
    bool deep_fq = flag.deep && (comp_i <= SAM_COMP_DEPN) && (flag.one_component-1 >= SAM_COMP_FQ00); // deep: --R1 or --R2
    #define has_sa (ZCTX(OPTION_SA_Z)->z_data_exists > 0)
    #define is_aux dict_id_is_aux_sf(dict_id) /* #define so calculated only when (rarely) needed*/

    switch (dict_id.num) {
        case _SAM_SQBITMAP : 
            SKIPIFF ((cov || cnt) && !flag.bases && 
                     !has_sa); // if this file has SA:Z: needed by MD:Z which is needed by NM:i which is needed by SA:Z
            
        case _SAM_NONREF   : case _SAM_NONREF_X : case _SAM_GPOS     : case _SAM_STRAND :
        case _SAM_SEQMIS_A : case _SAM_SEQMIS_C : case _SAM_SEQMIS_G : case _SAM_SEQMIS_T : 
        case _SAM_SEQINS_A : case _SAM_SEQINS_C : case _SAM_SEQINS_G : case _SAM_SEQINS_T : 
            SKIPIF (is_prim); // in PRIM, we skip sections that we used for loading the SA Groups in sam_piz_load_sags, but not needed for reconstruction
                           // (during PRIM SA Group loading, skip function is temporarily changed to sam_plsg_only). see also: sam_load_groups_add_grps
            SKIPIFF ((cov || cnt) && !flag.bases && !has_sa);

        case _SAM_QUAL  : case _SAM_DOMQRUNS  : case _SAM_QUALMPLX  : case _SAM_DIVRQUAL  :
        case _SAM_CQUAL : case _SAM_CDOMQRUNS : case _SAM_CQUALMPLX : case _SAM_CDIVRQUAL : 
            SKIPIF (is_prim);                                         
            SKIPIFF (cov || cnt);

        case _OPTION_np_i : case _OPTION_ec_f : // np is needed for reconstruting QUAL (PACB codec), and ec is needed to rereconstruct np
        case _OPTION_iq_sq_dq : case SAM_QUAL_PACBIO_DIFF: // iq_sq_dq, if it exists, is combined with DIFF to construct QUAL
            SKIPIFF (cov || cnt);

        case _SAM_TLEN    :
        case _SAM_BAM_BIN :
        case _SAM_EOL     :
            SKIPIF (deep_fq);
            // fallthrough

        case _SAM_QUALSA  :
            SKIPIFF (preproc || cov || (cnt && !(flag.bases && (OUT_DT(BAM) || OUT_DT(CRAM)))));

        case _SAM_Q1NAME : case _SAM_QNAMESA :
            KEEPIF (preproc || dict_needed_for_preproc || (cnt && flag.bases && (OUT_DT(BAM) || OUT_DT(CRAM)))); // if output is BAM we need the entire BAM record to correctly analyze the SEQ for IUPAC, as it is a structure.
            SKIPIFF (is_prim);                                         

        case _OPTION_SA_Z :
        case _OPTION_SA_MAIN :
            KEEPIF (IS_SAG_SA && preproc && ST(LOCAL));
            SKIPIF (IS_SAG_SA && is_prim && ST(LOCAL));
            KEEP; // need to reconstruct fields (RNAME, POS etc) against saggy
                     
        case _OPTION_SA_RNAME  : case _OPTION_SA_POS : case _OPTION_SA_NM : case _OPTION_SA_CIGAR :
        case _OPTION_SA_STRAND : case _OPTION_SA_MAPQ : 
            KEEPIF (IS_SAG_SA && (preproc || dict_needed_for_preproc));
            KEEPIF (is_main || is_dict); // need to reconstruct from prim line
            SKIPIFF (IS_SAG_SA && is_prim);    
            
        case _OPTION_AS_i  : // we don't skip AS in preprocessing unless it is entirely skipped
            SKIPIF (preproc && !segconf.sag_has_AS);
            SKIPIF (deep_fq);
            SKIPIFF ((cov || cnt) && is_aux);
  
        // data stored in SAGs - needed for reconstruction, and also for preproccessing
        case _OPTION_NH_i: 
        case _OPTION_CR_Z: case _OPTION_CR_Z_X: case _OPTION_CB_Z: case _OPTION_CB_ARR: case _OPTION_CB_SUFFIX:
        case _OPTION_RX_Z: case _OPTION_RX_Z_X: case _OPTION_BX_Z:  
        case _OPTION_BC_Z: case _OPTION_BC_ARR:
        case _OPTION_CY_Z: case _OPTION_CY_ARR: case _OPTION_CY_DIVRQUAL: case _OPTION_CY_DOMQRUNS: case _OPTION_CY_QUALMPLX:
        case _OPTION_QT_Z: case _OPTION_QT_ARR: case _OPTION_QT_DIVRQUAL: case _OPTION_QT_DOMQRUNS: case _OPTION_QT_QUALMPLX:
        case _OPTION_QX_Z:                      case _OPTION_QX_DIVRQUAL: case _OPTION_QX_DOMQRUNS: case _OPTION_QX_QUALMPLX:
            SKIPIFF (deep_fq || cov || cnt);

        case _SAM_FQ_AUX   : KEEPIFF (OUT_DT(FASTQ) || flag.collect_coverage);
        
        case _SAM_FLAG     : KEEP; // needed for demultiplexing by has_prim
        case _SAM_BUDDY    : KEEP; // always needed (if any of these are needed: QNAME, FLAG, MAPQ, CIGAR...)
        case _SAM_QNAME    : KEEP; // always needed as it is used to determine whether this line has a buddy
        case _SAM_RNAME    : KEEP;
        case _SAM_SAG      : KEEP;
        case _SAM_SAALN    : KEEP;
        case _SAM_AUX      : KEEP; // needed in preproc for container_peek_get_idxs and codec_pacb_piz_get_np

        case _OPTION_MQ_i  : // needed for MAPQ reconstruction (mate copy)
        case _OPTION_xq_i  : // needed for MAPQ reconstruction
        case _SAM_MAPQ     : // required for SA:Z reconstruction, which is required for RNAME, POS etc if segged against saggy
        case _OPTION_NM_i  : // required for SA:Z reconstruction
        case _OPTION_MD_Z  : // required NM:i reconstruction
            if (!has_sa) goto other; // not needed if this file has no SA:Z
            SKIPIFF (preproc); 
            
        case _OPTION_MC_Z  :          
        case _SAM_CIGAR    : SKIPIFF (preproc && IS_SAG_SA);
        case _SAM_RNEXT    : // Required by RNAME and PNEXT
        case _SAM_PNEXT    : // PNEXT is required by POS
        case _SAM_POS      : SKIPIFF (preproc && IS_SAG_SA);
        case _SAM_TOPLEVEL : KEEPIF (flag.deep && is_dict); // needed to reconstruct the FASTQ components
                             SKIPIFF (preproc || OUT_DT(BAM) || OUT_DT(CRAM) || OUT_DT(FASTQ));
        case _SAM_TOP2BAM  : SKIPIFF (preproc || OUT_DT(SAM) || OUT_DT(FASTQ));
        case 0             : KEEPIFF (ST(VB_HEADER) || !preproc);
        
        default            : other :
            SKIPIF ((cov || cnt) && is_aux && !(is_dict && dict_id.num == _OPTION_OC_Z/*dict-alias of MC0_Z needed to reconstruct CIGAR*/));
            SKIPIF (deep_fq && is_aux && !is_dict); // Dictionaries might be needed for FASTQ AUX fields
            SKIPIF (deep_fq && is_dict && is_aux && !f.dictionary.deep_fastq); // Deep --R1/2: we don't need SAM-only AUX dicts
            SKIPIF (is_dict && !flag.deep && is_aux && f.dictionary.deep_fastq && !f.dictionary.deep_sam);   // Deep file, but we are not reconstructed FQ (hence !flag.deep) - we don't need FASTQ-only AUX dicts
            SKIPIFF (preproc);
    }

    #undef KEEP
    #undef SKIP
}

//-----------------------------
// --seqs-file filter
//-----------------------------

typedef struct {
    uint32_t hash; // hash of seq
    uint32_t seq_len    : 31;
    uint32_t is_alloced : 1;
    rom seq;
} SeqFilterItem;
static Buffer seqs_filter = {};

static ASCENDING_SORTER (seq_filter_sort_by_hash, SeqFilterItem, hash)
static BINARY_SEARCHER (find_seq_in_filter, SeqFilterItem, uint32_t, hash, false)

// initialize seqs_filter and flag.seq_filter from command line
void seq_filter_initialize (rom filename)
{
    flag.seq_filter = (filename[0] == '^') ? -1 : 1;

    if (flag.seq_filter == -1) filename++; // negative filter

    file_split_lines (filename, "seqs_file", VERIFY_ASCII);

    // 2 entries per requested seq - forward and revcomp
    ARRAY_alloc (SeqFilterItem, seq, n_lines * 2, true, seqs_filter, evb, "seqs_filter");
    
    // forward entries 
    for (int i=0; i < n_lines; i++) 
        seq[i] = (SeqFilterItem){ .seq     = lines[i], // pointer into data buffer defined in file_split_lines
                                  .seq_len = line_lens[i],
                                  .hash    = crc32 (0, STRi(line,i)) };

    // revcomp entries
    for (int i=0; i < n_lines; i++) {
        seq[n_lines + i].seq        = MALLOC (line_lens[i]);
        seq[n_lines + i].seq_len    = line_lens[i];
        str_revcomp ((char *)seq[n_lines + i].seq, STRi(line, i));
        seq[n_lines + i].hash       = crc32 (0, STRa(seq[n_lines + i].seq));
        seq[n_lines + i].is_alloced = true;
    }

    qsort (STRb(seqs_filter), sizeof(SeqFilterItem), seq_filter_sort_by_hash);

    // remove dups (note: only consecutive dups are removed - if 2+ seqs have the same hash, dups might be interleaved and not removed)
    int next = 1;
    for (int i=1; i < n_lines * 2; i++)
        if ((seq[i].hash != seq[i-1].hash || seq[i].seq_len != seq[i-1].seq_len ||
            memcmp (seq[i].seq, seq[i-1].seq, seq[i].seq_len))) { // not dup
            
            if (i != next) seq[next] = seq[i];
            next++;
        }
    seqs_filter.len32 = next;

    // display sorted list of (hash,seq) pairs (uncomment for debugging)
    // for_buf (SeqFilterItem, e, seqs_filter) iprintf ("%u %s\n", e->hash, e->seq);
}

static void seq_filter_destroy (void)
{
    for_buf (SeqFilterItem, seq, seqs_filter)
        if (seq->is_alloced) FREE (seq->seq);
}

static bool sam_piz_line_survives_seq_filter (STRp(seq))
{
    uint32_t hash = crc32 (0, STRa(seq));
    SeqFilterItem *ent = binary_search (find_seq_in_filter, SeqFilterItem, seqs_filter, hash);

    bool found = false;

    for (; ent && ent < BAFT (SeqFilterItem, seqs_filter) && ent->hash == hash; ent++)
        if (str_issame (seq, ent->seq)) {
            found = true;
            break;
        }

    return (flag.seq_filter == 1) ? found : !found; // positive or negative filter
}

// set --FLAG filtering from command line argument
void sam_set_FLAG_filter (rom optarg)
{
    #define FLAG_ERR "Bad argument of --FLAG: \"%s\". It should be one of + - ^ (+:INCLUDE_IF_ALL ; -:INCLUDE_IF_NONE ; ^:EXCLUDE_IF_ALL) followed by a decimal or hexadecimal integer (eg 0x1c). These values (and their prefixes) are also accepted in lieu of a number: MULTI, ALIGNED, UNMAPPED, NUNMAPPED, REVCOMP, NREVCOMP, FIRST, LAST, SECONDARY, FILTERED, DUPLICATE, SUPPLEMENTARY"
    switch (optarg[0]) {
        case '+': flag.sam_flag_filter = SAM_FLAG_INCLUDE_IF_ALL  ; break;
        case '-': flag.sam_flag_filter = SAM_FLAG_INCLUDE_IF_NONE ; break;
        case '^': flag.sam_flag_filter = SAM_FLAG_EXCLUDE_IF_ALL  ; break;
        default : ABORTINP (FLAG_ERR, optarg);
    }

    STR(value) = strlen (optarg)-1;
    value = &optarg[1];
    
    if (str_get_int_range_allow_hex16 (STRa(value), 1, 65535, &flag.FLAG)) {} // done
    else if (!strncmp (value, "MULTI",         value_len)) flag.FLAG = SAM_FLAG_MULTI_SEG;
    else if (!strncmp (value, "ALIGNED",       value_len)) flag.FLAG = SAM_FLAG_IS_ALIGNED;
    else if (!strncmp (value, "UNMAPPED",      value_len)) flag.FLAG = SAM_FLAG_UNMAPPED;
    else if (!strncmp (value, "NUNMAPPED",     value_len)) flag.FLAG = SAM_FLAG_NEXT_UNMAPPED;
    else if (!strncmp (value, "REVCOMP",       value_len)) flag.FLAG = SAM_FLAG_REV_COMP;
    else if (!strncmp (value, "NREVCOMP",      value_len)) flag.FLAG = SAM_FLAG_NEXT_REV_COMP;
    else if (!strncmp (value, "FIRST",         value_len)) flag.FLAG = SAM_FLAG_IS_FIRST;
    else if (!strncmp (value, "LAST",          value_len)) flag.FLAG = SAM_FLAG_IS_LAST;
    else if (!strncmp (value, "SECONDARY",     value_len)) flag.FLAG = SAM_FLAG_SECONDARY;
    else if (!strncmp (value, "FILTERED",      value_len)) flag.FLAG = SAM_FLAG_FILTERED;
    else if (!strncmp (value, "DUPLICATE",     value_len)) flag.FLAG = SAM_FLAG_DUPLICATE;
    else if (!strncmp (value, "SUPPLEMENTARY", value_len)) flag.FLAG = SAM_FLAG_SUPPLEMENTARY;
    else ABORTINP (FLAG_ERR, optarg);
}

// set --MAPQ filtering from command line argument
void sam_set_MAPQ_filter (rom optarg)
{
    #define MAPQ_ERR "Bad argument of --MAPQ: \"%s\". It should be a number 0-255 (INCLUDE lines with MAPQ of at least this) or ^ (eg ^1) (EXCLUDE lines with MAPQ of at least this)"
    
    if (optarg[0] == '^') {
        flag.sam_mapq_filter = SAM_MAPQ_EXCLUDE_IF_AT_LEAST;
        optarg++;
    }
    else 
        flag.sam_mapq_filter = SAM_MAPQ_INCLUDE_IF_AT_LEAST;

    ASSERT (str_get_int_range8 (optarg, 0, 0, 255, &flag.MAPQ), MAPQ_ERR, optarg);
}

static inline void sam_piz_update_coverage (VBlockP vb, const uint16_t sam_flag, uint32_t soft_clip)
{
    ARRAY (uint64_t, read_count, vb->read_count);
    ARRAY (uint64_t, coverage, vb->coverage);
    uint64_t *coverage_special   = BAFT64 (vb->coverage)   - NUM_COVER_TYPES;
    uint64_t *read_count_special = BAFT64 (vb->read_count) - NUM_COVER_TYPES;
    WordIndex chrom_index = vb->last_index(SAM_RNAME);

    if (chrom_index == WORD_INDEX_NONE ||
             sam_flag & SAM_FLAG_UNMAPPED)       { coverage_special[CVR_UNMAPPED]      += vb->seq_len; read_count_special[CVR_UNMAPPED]     ++; }
    else if (sam_flag & SAM_FLAG_FILTERED)       { coverage_special[CVR_FAILED]        += vb->seq_len; read_count_special[CVR_FAILED]       ++; }
    else if (sam_flag & SAM_FLAG_DUPLICATE)      { coverage_special[CVR_DUPLICATE]     += vb->seq_len; read_count_special[CVR_DUPLICATE]    ++; }
    else if (sam_flag & SAM_FLAG_SECONDARY)      { coverage_special[CVR_SECONDARY]     += vb->seq_len; read_count_special[CVR_SECONDARY]    ++; }
    else if (sam_flag & SAM_FLAG_SUPPLEMENTARY)  { coverage_special[CVR_SUPPLEMENTARY] += vb->seq_len; read_count_special[CVR_SUPPLEMENTARY]++; }
    else {
        coverage_special[CVR_SOFT_CLIP] += soft_clip;
        coverage[chrom_index] += vb->seq_len - soft_clip;
        read_count[chrom_index]++;
    }
}

// Case 1: BIN is set to SPECIAL, we will set new_value here to -1 and wait for CIGAR to calculate it, 
//         as we need vb->ref_consumed - sam_cigar_special_CIGAR will update the reconstruced value
// Case 2: BIN is an textual integer snip - its BIN.last_value will be set as normal and transltor will reconstruct it
SPECIAL_RECONSTRUCTOR (bam_piz_special_BIN)
{
    ctx->semaphore = true; // signal to sam_cigar_special_CIGAR to calculate
    return NO_NEW_VALUE;
}

// note of float reconstruction:
// When compressing SAM, floats are stored as textual string in earlier version, reconstruced 
//    natively for SAM and via sam_piz_sam2bam_FLOAT for BAM. This way, the correct (as in the SAM file) number of textual digits is reconstructed.
// When compressing BAM, floats are stored as 32-bit binaries (in local since v15 and encoded as uint32 in the snip up to v14), 
//    They are reconstructed, either as textual for SAM or binary for BAM via bam_piz_special_FLOAT. This way, the BAM binary float is reconstructed precisely.
SPECIAL_RECONSTRUCTOR (bam_piz_special_FLOAT)
{
    union { uint32_t i; float f; } machine_en; // number in machine endianity

    // case: since v15, float is *may* stored in local
    if (!snip_len) 
        machine_en.f = NEXTLOCAL(float, ctx);
    
    // case: up to v14, float was stored in the snip, as a binary-identical uint32_t
    else {
        int64_t n;
        ASSERT (str_get_int (snip, snip_len, &n), "failed to read integer in %s", ctx->tag_name);
        
        machine_en.i = (uint32_t)n;
    }

    if (!reconstruct) goto finish;

    // binary reconstruction in little endian - BAM format
    if (OUT_DT(BAM) || OUT_DT(CRAM)) {
        uint32_t n32_lten = LTEN32 (machine_en.i); // little endian (per BAM format spec)
        RECONSTRUCT (&n32_lten, sizeof (uint32_t)); // in binary - float and uint32 are the same
    }

    // textual reconstruction - SAM format 
    else { 
        #define NUM_SIGNIFICANT_DIGITS 6 // 6 significant digits, as samtools does
        
        // calculate digits before and after the decimal point
        double log_f = log10 (machine_en.f >= 0 ? machine_en.f : -machine_en.f);
        unsigned int_digits = (log_f >= 0) + (unsigned)log_f;
        unsigned dec_digits = MAX_(0, NUM_SIGNIFICANT_DIGITS - int_digits);
        
        // reconstruct number with exactly NUM_SIGNIFICANT_DIGITS digits
        snprintf (BAFTtxt, Rtxt, "%.*f", dec_digits, machine_en.f); 
        unsigned len = strlen (BAFTtxt); 
        Ltxt += len;

        // remove trailing decimal zeros:  "5.500"->"5.5" ; "5.0000"->"5" ; "50"->"50"
        if (dec_digits) {
            unsigned trailing_zeros=0;
            for (int i=Ltxt-1; i >= Ltxt-dec_digits; i--)
                if (*Btxt (i) == '0') 
                    trailing_zeros++;
                else
                    break;
            
            Ltxt -= (dec_digits==trailing_zeros) ? dec_digits+1 : trailing_zeros;
        }
    }

finish:
    new_value->f = (double)machine_en.f;
    return true; // have new value
}

//-----------------------------------------------------------------
// Translator functions for reconstructing SAM data into BAM format
//-----------------------------------------------------------------

// output the word_index of RNAME, which is verified in ref_contigs_get_ref_chrom during seg
// to be the same as the reference id 
TRANSLATOR_FUNC (sam_piz_sam2bam_RNAME)
{
    STR0(snip);
    ctx_get_snip_by_word_index (ctx, ctx->last_value.i, snip);

    // if it is '*', reconstruct -1
    if (IS_ASTERISK(snip))
        RECONSTRUCT_BIN32 (-1);

    // if its RNEXT and =, emit the last index of RNAME
    else if (ctx->dict_id.num == _SAM_RNEXT && IS_EQUAL_SIGN (snip)) 
        RECONSTRUCT_BIN32 (CTX(SAM_RNAME)->last_value.i);

    // otherwise - output the word_index which was stored here because of flags.store=STORE_INDEX set in seg 
    else     
        RECONSTRUCT_BIN32 (ctx->last_value.i); 
    
    return 0;
}

// output, in binary form, POS-1 as BAM uses 0-based POS
TRANSLATOR_FUNC (sam_piz_sam2bam_POS)
{
    RECONSTRUCT_BIN32 (ctx->last_value.i - 1);
    return 0;
}

// translate AUX SAM->BAM - called as translator-only item on within the Aux reconstruction
// fix prefix eg MX:i: -> MXs
TRANSLATOR_FUNC (sam_piz_sam2bam_AUX_SELF)
{
    ContainerP con = (ContainerP)recon;

    // if this translator is called due to SAM->FASTQEXT, cancel translation for the AUX container and its items
    if (OUT_DT(FASTQ)) {
        con->no_translation = true; // turn off translation for this container and its items
        return 0; 
    }

    if (recon_len == -1) return 0; // no Aux data in this alignment

    char *prefixes_before = &recon[con_sizeof (*con)] + 2; // +2 to skip the empty prefixes of container wide, and item[0]
    char *prefixes_after = prefixes_before;

    uint32_t num_items = con_nitems (*con);

    for (unsigned i=1; i < num_items; i++, prefixes_before+=6, prefixes_after+=4) {
        prefixes_after[0] = prefixes_before[0]; // tag[0] 
        prefixes_after[1] = prefixes_before[1]; // tag[1]
        prefixes_after[2] = prefixes_before[3]; // type
        prefixes_after[3] = SNIP_CONTAINER; // end of prefix

        // a SAM 'i' translate to one of several BAM types using the translator code
        // that may be 0->6 (NONE to SAM2BAM_LTEN_U32)
        if (prefixes_after[2] == 'i') 
            prefixes_after[2] = "\0cCsSiI"[con->items[i].translator];
    }

    return -2 * (num_items-1); // change in prefixes_len
}

// translate AUX SAM->BAM - called after Aux reconstruction is done
TRANSLATOR_FUNC (sam_piz_sam2bam_AUX)
{
    /* up to v11 we used this translator to set alignment.block_size. due to container logic change, it
       has now moved to sam_piz_container_cb. we keep this translator function for backward
       compatability with v11 bam.genozip files */
    return 0;
}

// note of float reconstruction:
// When compressing SAM, floats are stored as a textual string, reconstruced natively for SAM and via sam_piz_sam2bam_FLOAT for BAM.
//    Done this way so when reconstructing SAM, the correct number of textual digits is reconstructed.
// When compressing BAM, floats are stored as 32-bit binaries, encoded as uint32, and stringified to a snip. They are reconstructed,
//    either as textual for SAM or binary for BAM via bam_piz_special_FLOAT. Done this way so BAM binary float is reconstructd precisely.
TRANSLATOR_FUNC (sam_piz_sam2bam_FLOAT)
{
    union {
        float f; // 32 bit float
        uint32_t i;
    } value;
    
    ASSERT0 (sizeof (value)==4, "expecting value to be 32 bits"); // should never happen

    value.f = (float)ctx->last_value.f;
    RECONSTRUCT_BIN32 (value.i);

    return 0;
}

// remove the comma from the prefix that contains the type, eg "i,"->"i"
TRANSLATOR_FUNC (sam_piz_sam2bam_ARRAY_SELF_1) // single context, multi repeat array
{
    ContainerP con = (ContainerP)recon;
    char *prefixes = &recon[con_sizeof (*con)];

    // remove the ',' from the prefix, and terminate with CON_PX_SEP_SHOW_REPEATS - this will cause
    // the number of repeats (in LTEN32) to be outputed after the prefix
    prefixes[1] = CON_PX_SEP_SHOW_REPEATS; // prefixes is now { type, CON_PX_SEP_SHOW_REPEATS }
    
    return con->repeats ? -1 : 0; // change in prefixes length (no change in length if repeats=0, because no comman in SAM eg "ML:B:C" - empty array)
}

// remove the comma from the prefix that contains the type, eg "i,"->"i"
TRANSLATOR_FUNC (sam_piz_sam2bam_ARRAY_SELF_M) // multiple context, single repeat array (15.0.35)
{
    ContainerP con = (ContainerP)recon;
    char *prefixes = &recon[con_sizeof (*con)];

    // remove the ',' from the prefix, and terminate with CON_PX_SEP_SHOW_N_ITEMS - this will cause
    // the number of items (in LTEN32) to be outputed after the prefix
    prefixes[1] = CON_PX_SEP_SHOW_N_ITEMS; // prefixes is now { type, CON_PX_SEP_SHOW_N_ITEMS }
    
    return -1; 
}

//------------------------------------------------------------------------------------
// Translator and filter functions for reconstructing SAM / BAM data into FASTQ format
//------------------------------------------------------------------------------------

TXTHEADER_TRANSLATOR (txtheader_sam2fq)
{
    txtheader_buf->len = 0; // fastq has no header
}

// filtering during reconstruction: called by container_reconstruct for each sam alignment (repeat)
CONTAINER_CALLBACK (sam_piz_container_cb)
{    
    if (is_top_level) {
        #define DROP_LINE(reason) ({ vb->drop_curr_line = (reason); goto dropped; })

        if (flag.add_line_numbers && TXT_DT(SAM)) {
            Ltxt -= 1 + (*(BLSTtxt-1) == '\r'); // remove \n or \r\n
            Ltxt += snprintf (BAFTtxt, Rtxt, "\tVB:Z:%s\n", LN_NAME);
        }

        // case SAM to BAM translation: set alignment.block_size (was in sam_piz_sam2bam_AUX until v11)
        if (dict_id.num == _SAM_TOP2BAM) { 
            BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)Btxt (vb->line_start);
            PUT_UINT32_(alignment, block_size, LTEN32 (Ltxt - vb->line_start - sizeof (uint32_t))); // block_size doesn't include the block_size field itself
        }
        
        // --FLAG
        if (flag.sam_flag_filter) {

            uint16_t this_sam_flag = (uint16_t)vb->last_int (SAM_FLAG);
        
            bool all_flags_set = (this_sam_flag & flag.FLAG) == flag.FLAG;
            bool no_flags_set  = (this_sam_flag & flag.FLAG) == 0;
            if ((flag.sam_flag_filter == SAM_FLAG_INCLUDE_IF_ALL  && !all_flags_set)
            ||  (flag.sam_flag_filter == SAM_FLAG_INCLUDE_IF_NONE && !no_flags_set)
            ||  (flag.sam_flag_filter == SAM_FLAG_EXCLUDE_IF_ALL  &&  all_flags_set))
                DROP_LINE ("FLAG");
        }

        // --MAPQ
        if (flag.sam_mapq_filter) {     
            uint8_t this_mapq = (uint8_t)vb->last_int (SAM_MAPQ);
        
            if ((flag.sam_mapq_filter == SAM_MAPQ_INCLUDE_IF_AT_LEAST && this_mapq < flag.MAPQ) ||
                (flag.sam_mapq_filter == SAM_MAPQ_EXCLUDE_IF_AT_LEAST && this_mapq >= flag.MAPQ))
                DROP_LINE ("MAPQ");
        }

        // --bases
        if (flag.bases && 
            !(TXT_DT(BAM) ? iupac_is_included_bam (last_txt (vb, SAM_SQBITMAP), ((BAMAlignmentFixed *)recon)->l_seq)
                          : iupac_is_included_ascii (STRlst (SAM_SQBITMAP))))
            DROP_LINE ("bases");
        
        // --qnames and --qnames-file
        if (flag.qname_filter && !qname_filter_does_line_survive (STRlst (IS_PRIM(vb) ? SAM_QNAMESA : SAM_QNAME))) // works also for FASTQ as SAM_QNAME==FASTQ_QNAME
            DROP_LINE ("qname_filter");

        // --seqs-file
        if (flag.seq_filter && !sam_piz_line_survives_seq_filter (STRlst (SAM_SQBITMAP))) // works also for FASTQ as SAM_SQBITMAP==FASTQ_SQBITMAP
            DROP_LINE ("seq_filter");

        // count coverage, if needed    
        if (flag.show_coverage)
            sam_piz_update_coverage (vb, vb->last_int(SAM_FLAG), VB_SAM->soft_clip[0] + VB_SAM->soft_clip[1]);

        if (flag.idxstats) {
            if (vb->last_int(SAM_FLAG) & SAM_FLAG_UNMAPPED)   
                (*B64 (vb->unmapped_read_count, vb->last_index(SAM_RNAME)))++;
            else
                (*B64 (vb->read_count, vb->last_index(SAM_RNAME)))++;
        }

        dropped: {}
    }
}

bool sam_piz_filter_up_to_v13_stuff (VBlockP vb, DictId dict_id, int item, bool *filter_ret_value)
{
    // BAM: set buddy at the beginning of each line, as QNAME is reconstructed much later
    // For MAIN and DEPN: buddy, if existing, will have QNAME=SNIP_COPY_BUDDY (see sam_seg_QNAME)
    // note: we always load buddy, to prevent a situation when in some lines it is consumed
    // and other lines, which have buddy, it is not consumed because the field that consumes it is skipped
    if (dict_id.num == _SAM_TOP2BAM && item == 0 && 
        CTX(SAM_QNAME)->b250.len32) { // might be 0 in special cases, like flag.count
        STR(snip);
        PEEK_SNIP(SAM_QNAME);
        if (snip_len && *snip == v13_SNIP_COPY_BUDDY)
            sam_piz_set_buddy_v13(vb);
        *filter_ret_value = true;
        return true;
    }

    // collect_coverage: set buddy_line_i here, since we don't reconstruct QNAME
    // note: we always load buddy, to prevent a situation when in some lines it is consumed
    // and other lines, which have buddy, it is not consumed because the field that consumes it is skipped
    else if (dict_id.num == _SAM_QNAME && flag.collect_coverage) {
        STR(snip);
        LOAD_SNIP(SAM_QNAME);
        if (snip_len && *snip == v13_SNIP_COPY_BUDDY)
            sam_piz_set_buddy_v13(vb);
        *filter_ret_value = false; // don't reconstruct QNAME
        return true;
    }

    // XA:Z: insert RNAME, STRAND and POS lookbacks - moved to sam_piz_container_cb in v14
    else if (dict_id.num == _OPTION_XA_Z && item == 4) {  // importantly, after RNAME and STRAND_POS
        sam_piz_XA_field_insert_lookback_v13 (vb);
        *filter_ret_value = true;
        return true;
    }

    else    
        return false; // nothing was handled here
} 

// filter is called before reconstruction of a repeat or an item, and returns false if item should 
// not be reconstructed. contexts are not consumed.
CONTAINER_FILTER_FUNC (sam_piz_filter)
{
    bool v13_ret_value;
    if (!VER(14) && sam_piz_filter_up_to_v13_stuff (vb, dict_id, item, &v13_ret_value))
        return v13_ret_value;

    else if (dict_id.num == _SAM_TOP2NONE) { // bug 857
        // if (con->items[item].dict_id.num == _SAM_SQBITMAP && (flag.header_only_fast || flag.qual_only)) 
        //     return false; // don't reconstruct SEQ

        // if (con->items[item].dict_id.num == _SAM_QUAL     && (flag.header_only_fast || flag.seq_only)) 
        //     return false; // don't reconstruct QUAL

        // if (con->items[item].dict_id.num == _SAM_QNAME    && (flag.qual_only || flag.seq_only))
        //     return false; // don't reconstruct QNAME
    }

    // collect_coverage: rather than reconstructing optional, reconstruct SAM_FQ_AUX that just consumes MC:Z if it exists
    else if (dict_id.num == _SAM_AUX) {
        if (flag.collect_coverage) {
            if      (CTX(SAM_FQ_AUX    )->is_loaded) reconstruct_from_ctx (vb, SAM_FQ_AUX,     0, false); // filter_repeats is set in the AUX container since v14
            else if (CTX(SAM_FQ_AUX_OLD)->is_loaded) reconstruct_from_ctx (vb, SAM_FQ_AUX_OLD, 0, false);
            else ABORT0 ("Neither SAM_FQ_AUX or SAM_FQ_AUX_OLD are loaded");

            return false; // don't reconstruct AUX
        }

        else
            VB_SAM->aux_con = con;
    }

    return true; // go ahead and reconstruct
}

bool sam_piz_line_has_aux_field (VBlockSAMP vb, DictId dict_id)
{
    ASSPIZ0 (vb->aux_con, "this function can only be called while reconstructing the AUX container");

    for_con (vb->aux_con)
        if (item->dict_id.num == dict_id.num) return true;

    return false;
}

// v14: De-multiplex by has_mate
SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_BY_MATE)
{
    return reconstruct_demultiplex (vb, ctx, STRa(snip), sam_has_mate, new_value, reconstruct);
}

// v14: De-multiplex by has_buddy (which can be mate or prim)
SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_BY_BUDDY)
{
    return reconstruct_demultiplex (vb, ctx, STRa(snip), sam_has_mate || sam_has_saggy, new_value, reconstruct);
}

// v14: De-multiplex by has_mate and has_prim
SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_BY_MATE_PRIM)
{
    // note: when reconstructing POS or MAPQ in BAM, FLAG is not known yet, so we peek it here
    if (!ctx_has_value_in_line_(vb, CTX(SAM_FLAG)))
        ctx_set_last_value (vb, CTX(SAM_FLAG), reconstruct_peek (vb, CTX(SAM_FLAG), 0, 0));

    int channel_i = sam_has_mate?1 : sam_has_prim?2 : 0;

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

//---------------------------------------------
// Consuming data stored by a "historical" line
//---------------------------------------------

// used starting v14
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_SET_BUDDY)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    ContextP buddy_ctx = CTX(SAM_BUDDY);

    // for prim VB we set sag here, as all lines have a sag. for depn, we will set as we encounter lines that use sage
    if (IS_PRIM(vb) && !vb->preprocessing)
        sam_piz_set_sag (VB_SAM); 

    // case 1: when reconstructing prim, we don't normally consume QNAME (we consume QNAMESA instead), so we consume it here
    // case 2: when loading SAG, we set buddy multiple times - as each field is loaded separately (iterators are reset after loading each field)
    //         we refrain to consume when loading the QNAME field, as QNAME itself we will be reconstructed
    //         when preprocessing, we take "reconstruct" to mean "consume qname"
    if (IS_PRIM(vb) && (!vb->preprocessing || (vb->preprocessing && reconstruct))) 
        LOAD_SNIP (SAM_QNAME);
    else
        PEEK_SNIP (SAM_QNAME);

    if (snip_len == 3 && snip[1] == SAM_SPECIAL_COPY_BUDDY) { // has buddy

        BuddyType bt = snip[2] - '0';

        // note: if both mate and saggy are available, the first one in BUDDY.local is mate

        if (bt & BUDDY_MATE) { // has mate
            int32_t num_lines_back = reconstruct_from_local_int (VB, buddy_ctx, 0, RECON_OFF);
            vb->mate_line_i = vb->line_i - num_lines_back;
        }

        if (bt & BUDDY_SAGGY) { // has saggy
            int32_t num_lines_back = reconstruct_from_local_int (VB, buddy_ctx, 0, RECON_OFF);
            vb->saggy_line_i  = vb->line_i - num_lines_back;
            vb->saggy_is_prim = !sam_is_depn ((SamFlags){ .value = history64(SAM_FLAG, vb->saggy_line_i)});
        }

        if (bt && flag.show_buddy) {
            if      (bt == BUDDY_EITHER) iprintf ("%s: mate_line_i=%u saggy_line_i=%u%s\n", LN_NAME, VB_SAM->mate_line_i, VB_SAM->saggy_line_i, vb->preprocessing ? " (preprocessing)" : "");
            else if (bt == BUDDY_MATE)   iprintf ("%s: mate_line_i=%u%s\n",  LN_NAME, VB_SAM->mate_line_i, vb->preprocessing ? " (preprocessing)" : "");
            else if (bt == BUDDY_SAGGY)  iprintf ("%s: saggy_line_i=%u%s\n", LN_NAME, VB_SAM->saggy_line_i, vb->preprocessing ? " (preprocessing)" : "");
        }
    }

    return NO_NEW_VALUE;
}

// used up to v13 - existance of a mate was conveyed by a QNAME snip that was { SNIP_COPY_BUDDY, SNIP_COPY_BUDDY },
// and the mate line itself was conveyed via SAM_BUDDY.local
void sam_piz_set_buddy_v13 (VBlockP vb)
{
    if (VB_SAM->mate_line_i != NO_LINE) return; // already set

    ContextP buddy_ctx = CTX(SAM_BUDDY);

    int32_t num_lines_back = reconstruct_from_local_int (vb, buddy_ctx, 0, RECON_OFF);

    // a bug that existed 12.0.41-13.0.1 (bug 367): we stored buddy in machine endianty instead of BGEN32.
    // v13 starting 13.0.2 set SAM_BUDDY.local.prm8[0] so can detect buggy files by local.prm8[0]=0 and convert it back to machine endianity.
    if (!buddy_ctx->local.prm8[0])    
        num_lines_back = BGEN32 ((uint32_t)num_lines_back);

    VB_SAM->mate_line_i = vb->line_i - num_lines_back; // convert value passed (distance in lines to buddy) to 0-based buddy_line_i

    if (flag.show_buddy)
        iprintf ("%s: mate_line_i=%u\n", LN_NAME, VB_SAM->mate_line_i);
        
    ASSPIZ (VB_SAM->mate_line_i != NO_LINE, "Expecting vb->mate_line_i=%d to be non-negative. num_lines_back=%d buddy_ctx->local.prm8[0]=%d", 
            VB_SAM->mate_line_i, num_lines_back, buddy_ctx->local.prm8[0]);
}

void sam_reconstruct_from_buddy_get_textual_snip (VBlockSAMP vb, ContextP ctx, BuddyType bt, pSTRp(snip))
{
    LineIType buddy_line_i = (bt == BUDDY_MATE) ? vb->mate_line_i : vb->saggy_line_i;

    recon_history_get_historical_snip (VB, ctx, buddy_line_i, snip, snip_len);
}

// Copy from buddy: buddy is data that appears on a specific "buddy line", in this context or another one. Not all lines need
// Note of difference vs. lookback: with buddy, not all lines need to have the data (eg MC:Z), so the line number is constant,
// but if we had have used lookback, the lookback value would have been different between different fields.  
// two options for snip:
// 1. single character - buddy_type    (up to v13 - no snip, implied BUDDY_MATE)
// 2. other_ctx followed by buddy type (up to v13 - only other_ctx, implied BUDDY_MATE)
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_COPY_BUDDY)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    BuddyType buddy_type = VER(14) ? (snip[snip_len-1] - '0') : BUDDY_MATE;
    ContextP base_ctx = ctx;

    // up to v13: set buddy if needed and not already set
    if (!VER(14) && snip_len==1 && *snip == v13_SNIP_COPY_BUDDY) {
        sam_piz_set_buddy_v13 (VB); 
        snip++;
        snip_len--;
    }

    // optional: base context is different than ctx
    if (snip_len > 1) 
        base_ctx = reconstruct_special_get_base_ctx (VB, ctx, pSTRa(snip));

    LineIType buddy_line_i = buddy_type == BUDDY_MATE   ? vb->mate_line_i
                           : buddy_type == BUDDY_SAGGY  ? vb->saggy_line_i
                           : vb->mate_line_i != NO_LINE ? ({ buddy_type = BUDDY_MATE;  vb->mate_line_i;  })  // BUDDY_EITHER and mate is available - take it (mate before saggy! seg may rely on this order)
                           :                              ({ buddy_type = BUDDY_SAGGY; vb->saggy_line_i; }); // BUDDY_EITHER - mate it not available, so saggy must be

    ASSPIZ (buddy_line_i != NO_LINE, "expecting buddy of type %s when copying from %s's %s to our %s but there is none", 
            buddy_type_name(buddy_type), buddy_type_name (buddy_type), base_ctx->tag_name, ctx->tag_name);

    // note: a non-existant STORE_INT buddy is taken as 0 (bc in seg, dl->field is 0 is non-existent and might be segged as buddy)
    ASSPIZ (ctx->flags.store == STORE_INT || buddy_line_i < base_ctx->history.len32, "history not set for %s, perhaps seg forgot to set store_per_line? (buddy_type=%s buddy_line_i=%d history.len=%u)", 
            base_ctx->tag_name, buddy_type_name (buddy_type), buddy_line_i, base_ctx->history.len32);

    // case: numeric value 
    if (ctx->flags.store == STORE_INT) {
        new_value->i = (buddy_line_i < base_ctx->history.len32) ? *B64(base_ctx->history, buddy_line_i) : 0; // note: a non-existant STORE_INT buddy is taken as 0
        if (reconstruct) RECONSTRUCT_INT (new_value->i);
        return HAS_NEW_VALUE;
    }

    // case: word index (use if we always store word_index. to mix word_index with other method, use LookupDict)
    else if (ctx->flags.store == STORE_INDEX) {
        new_value->i = *B(WordIndex, base_ctx->history, buddy_line_i);

        if (reconstruct) {
            ctx_get_snip_by_word_index (ctx, new_value->i, snip);
            RECONSTRUCT_snip;
        }
            
        return HAS_NEW_VALUE;
    }

    // case: textual value
    else {
        if (reconstruct) {
            sam_reconstruct_from_buddy_get_textual_snip (vb, base_ctx, buddy_type, pSTRa(snip));
            RECONSTRUCT_snip;

            if (ctx->did_i == SAM_QNAME && buddy_type == BUDDY_MATE && segconf.flav_prop[QNAME1].is_mated)
                *BLSTtxt = (*BLSTtxt=='1') ? '2' : '1';

            ctx_set_last_value (VB, ctx, (int64_t)buddy_line_i); // this context has STORE_NONE so this can't conflict with anything
        }

        return NO_NEW_VALUE; 
    }
}

// invoked from TOP2FQ and TOP2NONE (but not TOP2FQEX, bc it reconstructs AUX) to consume some AUX fields if they exists 
// in this line, in case this line - AUX fields that might be needed to reconstruct required fields in 
// subsequent lines 
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_FASTQ_CONSUME_AUX)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    ContainerPeekItem peek_items[] = { 
        { _OPTION_MC_Z, -1 }, { _OPTION_SA_Z, -1 }, { _OPTION_ec_f, -1 }, { _OPTION_np_i, -1 }, 
        { _OPTION_iq_Z, -1 }, { _OPTION_dq_Z, -1 }, { _OPTION_sq_Z, -1 } };

    vb->aux_con = container_peek_get_idxs (VB, CTX(SAM_AUX), ARRAY_LEN(peek_items), peek_items, true);

    // if this line has an MC:Z field store it directly history - might be needed for subsequent line CIGAR
    if (CTX(OPTION_MC_Z)->flags.store_per_line && CTX(OPTION_MC_Z)->is_loaded && // MC:Z is buddied
        peek_items[0].idx != -1) // line has MC:Z field
        reconstruct_to_history (VB, CTX(OPTION_MC_Z));

    // likewise for SA:Z - might be needed for subsequent line POS 
    if (CTX(OPTION_SA_Z)->flags.store_per_line && CTX(OPTION_SA_Z)->is_loaded &&
        peek_items[1].idx != -1) // line has SA:Z field
        reconstruct_to_history (VB, CTX(OPTION_SA_Z));

    // np:i and ec:f are PacBio tags used to reconstruct QUAL,iq,dq,sq - we just need the last_value, no need to store history
    if (CTX(OPTION_ec_f)->is_loaded && peek_items[2].idx != -1) // line has ec:f field
        reconstruct_from_ctx (VB, OPTION_ec_f, 0, false); 

    if (CTX(OPTION_np_i)->is_loaded && peek_items[3].idx != -1) // line has np:i field
        reconstruct_from_ctx (VB, OPTION_np_i, 0, false); 

    return NO_NEW_VALUE; 
}

bool sam_is_last_flags_rev_comp (VBlockP vb)
{
    return last_flags.rev_comp;
}
