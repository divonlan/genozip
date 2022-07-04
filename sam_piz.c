// ------------------------------------------------------------------
//   sam_piz.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <math.h>
#include "sam_private.h"
#include "seg.h"
#include "context.h"
#include "piz.h"
#include "reconstruct.h"
#include "strings.h"
#include "dict_id.h"
#include "codec.h" // must be included before reference.h
#include "reference.h"
#include "regions.h"
#include "aligner.h"
#include "file.h"
#include "container.h"
#include "coverage.h"
#include "bases_filter.h"
#include "endianness.h"
#include "lookback.h"
#include "qname.h"
#include "writer.h"

void sam_piz_xtra_line_data (VBlockP vb_)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    if (SAM_PIZ_HAS_SA_GROUP) {
        iprintf ("grp_i=%u aln_i=%"PRIu64" (%d within group)\n", ZGRP_I(vb->sag), ZALN_I(vb->sa_aln), (int)(ZALN_I(vb->sa_aln) - vb->sag->first_aln_i));
        sam_show_sag_one_grp (ZGRP_I(VB_SAM->sag));
    }
}

void sam_piz_genozip_header (const SectionHeaderGenozipHeader *header)
{
    if (VER(14)) {
        segconf.sam_seq_len           = BGEN32 (header->sam.segconf_seq_len); 
        segconf.seq_len_to_cm         = header->sam.segconf_seq_len_cm;
        segconf.sam_ms_type           = header->sam.segconf_ms_type;
        segconf.has_MD_or_NM          = header->sam.segconf_has_MD_or_NM;
        segconf.sam_bisulfite         = header->sam.segconf_bisulfite;
        segconf.is_paired             = header->sam.segconf_is_paired;
        segconf.sag_type              = header->sam.segconf_sag_type;
        segconf.sag_has_AS       = header->sam.segconf_sag_has_AS;
        segconf.AS_is_2ref_consumed   = header->sam.segconf_AS_is_2refc;
        segconf.pysam_qual            = header->sam.segconf_pysam_qual;
        segconf.qname_seq_len_dict_id = header->sam.segconf_seq_len_dict_id; 
    }
}

// main thread: is it possible that genocat of this file will re-order lines
bool sam_piz_maybe_reorder_lines (void)
{
    return z_file->z_flags.has_gencomp; // has PRIM or DEPN components
}

bool sam_piz_init_vb (VBlockP vb, const SectionHeaderVbHeader *header, uint32_t *txt_data_so_far_single_0_increment)
{
    if (vb->comp_i == SAM_COMP_PRIM)
        VB_SAM->plsg_i = sam_piz_get_plsg_i (vb->vblock_i);

    return true; // all good
}

// PIZ compute thread: called after uncompressing contexts and before reconstructing
void sam_piz_recon_init (VBlockP vb)
{
    // we can proceed with reconstructing a PRIM or DEPN vb, only after SA Groups is loaded - busy-wait for it
    while ((sam_is_prim_vb || sam_is_depn_vb) && !sam_is_SA_Groups_loaded())
        usleep (250000); // 250 ms

    if (CTX(OPTION_XA_Z)->dict.len) { // this file has XA
        lookback_init (vb, CTX(OPTION_XA_LOOKBACK), CTX (OPTION_XA_RNAME),  STORE_INDEX);
        lookback_init (vb, CTX(OPTION_XA_LOOKBACK), CTX (OPTION_XA_STRAND), STORE_INDEX);
        lookback_init (vb, CTX(OPTION_XA_LOOKBACK), CTX (OPTION_XA_POS),    STORE_INT);
    }

    buf_alloc_zero (vb, &CTX(SAM_CIGAR)->ref_consumed_history, 0, vb->lines.len, uint32_t, 0, "ref_consumed_history"); // initialize to exactly one per line.
}

// PIZ: piz_after_recon callback: called by the compute thread from piz_reconstruct_one_vb. order of VBs is arbitrary
void sam_piz_after_recon (VBlockP vb)
{
}

void sam_piz_finalize (void)
{
    sam_header_finalize();
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

    // if this is a mux channel of an OPTION - consider its parent instead. Eg. "M0C:Z0", "M1C:Z1" -> "MC:i"
    if (dict_id.id[3]==':' && dict_id.id[1] == dict_id.id[5]&& !dict_id.id[6])
        dict_id = (DictId){ .id = { dict_id.id[0], dict_id.id[2], ':', dict_id.id[4] } };

    uint64_t dnum = dict_id.num;
    bool preproc = (purpose == SKIP_PURPOSE_PREPROC) && (st == SEC_B250 || st == SEC_LOCAL); // when loading SA, don't skip the needed B250/LOCAL
    bool dict_needed_for_preproc = (st == SEC_DICT && z_file->z_flags.has_gencomp);  // when loading a SEC_DICT in a file that has gencomp, don't skip dicts needed for loading SA

    if (dict_id_is_qname_sf(dict_id)) dnum = _SAM_Q1NAME; // treat all QNAME subfields as _SAM_Q1NAME
    bool prim = (comp_i == SAM_COMP_PRIM) && !preproc;
    bool main = (comp_i == SAM_COMP_MAIN);
    bool cov  = flag.collect_coverage;
    bool cnt  = flag.count && !flag.grep; // we skip based on --count, but not if --grep, because grepping requires full reconstruction
    
    switch (dnum) {
        case _SAM_SQBITMAP : 
            SKIPIFF ((cov || cnt) && !flag.bases);
            KEEPIFF (st == SEC_B250);
            
        case _SAM_NONREF   : case _SAM_NONREF_X : case _SAM_GPOS     : case _SAM_STRAND :
        case _SAM_SEQMIS_A : case _SAM_SEQMIS_C : case _SAM_SEQMIS_G : case _SAM_SEQMIS_T : 
            SKIPIF (prim); // in PRIM, we skip sections that we used for loading the SA Groups in sam_piz_load_sags, but not needed for reconstruction
                           // (during PRIM SA Group loading, skip function is temporarily changed to sam_plsg_only). see also: sam_load_groups_add_grps
            SKIPIFF ((cov || cnt) && !flag.bases);

        case _SAM_QUAL : case _SAM_DOMQRUNS : case _SAM_QUALMPLX : case _SAM_DIVRQUAL :
            SKIPIF (prim);                                         
            SKIPIFF (cov || cnt);
        
        case _SAM_QUALSA  :
        case _SAM_TLEN    :
        case _SAM_EOL     :
        case _SAM_BAM_BIN :
            SKIPIFF (preproc || cov || (cnt && !(flag.bases && flag.out_dt == DT_BAM)));

        case _SAM_RNEXT   :
            SKIPIFF ((preproc && IS_SAG_SA) || (cnt && !flag.bases)); // needed to reconstruct RNAME from a mate

        case _SAM_Q1NAME : case _SAM_QNAMESA :
            KEEPIF (preproc || dict_needed_for_preproc || (cnt && flag.bases && flag.out_dt == DT_BAM)); // if output is BAM we need the entire BAM record to correctly analyze the SEQ for IUPAC, as it is a structure.
            SKIPIFF (prim);                                         

        case _OPTION_SA_Z :
        case _OPTION_SA_MAIN :
            KEEPIF (preproc && st == SEC_LOCAL);
            SKIPIF (prim && st == SEC_LOCAL);
            SKIPIFF (cnt && !flag.bases);
                     
        case _OPTION_SA_RNAME  :    
        case _OPTION_SA_POS    :    
        case _OPTION_SA_CIGAR  :
        case _OPTION_SA_STRAND : 
        case _OPTION_SA_NM     :
            KEEPIF (preproc || dict_needed_for_preproc);
            SKIPIF (cnt && !flag.bases); 
            KEEPIF (main || st == SEC_DICT); // need to reconstruct from prim line
            SKIPIFF (prim);    

        case _OPTION_SA_MAPQ : 
            KEEPIF (preproc || dict_needed_for_preproc || main || st == SEC_DICT);
            SKIPIF (cnt && !flag.sam_mapq_filter);
            KEEPIF (main || st == SEC_DICT);
            SKIPIFF (prim);                                         
        
        case _SAM_FLAG     : 
            KEEPIF (preproc || dict_needed_for_preproc);            
            SKIPIFF (cnt && !flag.sam_flag_filter && !flag.bases);
            
        case _OPTION_NH_i  : 
            SKIPIFF (preproc || cnt || cov);

        case _OPTION_AS_i  : // we don't skip AS in preprocessing unless it is entirely skipped
            SKIPIF (preproc && !segconf.sag_has_AS);
            SKIPIFF ((cov || cnt) && dict_id_is_aux_sf(dict_id));
  
        case _OPTION_NH_PRIM :
            KEEPIFF (preproc || dict_needed_for_preproc);        

        case _OPTION_CR_Z: case _OPTION_CB_Z: case _OPTION_CR_CB:
        case _OPTION_CY_Z: case _OPTION_CY_ARR: case _OPTION_CY_DIVRQUAL: case _OPTION_CY_DOMQRUNS: case _OPTION_CY_QUALMPLX:
        case _OPTION_UR_Z: // note: UB is an alias so no data
        case _OPTION_UY_Z: case _OPTION_UY_DIVRQUAL: case _OPTION_UY_DOMQRUNS: case _OPTION_UY_QUALMPLX:
        case _OPTION_gx_Z:
        case _OPTION_gn_Z:
        case _OPTION_GX_Z:
        case _OPTION_GN_Z:  
            SKIPIFF (cov || cnt);

        case _SAM_MC_Z     : KEEPIFF (flag.out_dt == DT_FASTQ);
        
        case _SAM_BUDDY    : KEEP; // always needed (if any of these are needed: QNAME, FLAG, MAPQ, CIGAR...)
        case _SAM_QNAME    : KEEP; // always needed as it is used to determine whether this line has a buddy
        case _SAM_RNAME    : KEEP;
        case _SAM_SAG  : KEEP;
        case _SAM_SAALN    : KEEP;
        case _SAM_CIGAR    : SKIPIFF ((preproc && IS_SAG_SA) || (cnt && !flag.bases));
        case _OPTION_MC_Z  : SKIPIFF (preproc || (cnt && !flag.bases));
        case _SAM_AUX      : SKIPIFF (preproc && (!segconf.sag_has_AS && !IS_SAG_SOLO));
        case _SAM_MAPQ     : 
        case _OPTION_MQ_i  : SKIPIFF (preproc || ((cov || cnt) && !flag.bases && !flag.sam_mapq_filter));
        case _SAM_PNEXT    : case _SAM_P0NEXT : case _SAM_P1NEXT : case _SAM_P2NEXT : case _SAM_P3NEXT :
        case _SAM_POS      : SKIPIFF ((preproc && IS_SAG_SA) || (cnt && !flag.regions && !flag.bases));
        case _SAM_TAXID    : SKIPIFF (preproc || flag.kraken_taxid == TAXID_NONE);
        case _SAM_TOPLEVEL : SKIPIFF (preproc || flag.out_dt == DT_BAM || flag.out_dt == DT_FASTQ);
        case _SAM_TOP2BAM  : SKIPIFF (preproc || flag.out_dt == DT_SAM || flag.out_dt == DT_FASTQ);
        case _SAM_TOP2FQ   : SKIPIFF (preproc || flag.out_dt == DT_SAM || flag.out_dt == DT_BAM || flag.extended_translation);
        case _SAM_TOP2FQEX : SKIPIFF (preproc || flag.out_dt == DT_SAM || flag.out_dt == DT_BAM || !flag.extended_translation);
        case 0             : KEEPIFF (st == SEC_VB_HEADER || !preproc);
        
        default            : 
            SKIPIF ((cov || cnt) && dict_id_is_aux_sf(dict_id));
            SKIPIFF (preproc);
    }

    #undef KEEP
    #undef SKIP
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
// When compressing SAM, floats are stored as a textual string, reconstruced natively for SAM and via sam_piz_sam2bam_FLOAT for BAM.
//    Done this way so when reconstructing SAM, the correct number of textual digits is reconstructed.
// When compressing BAM, floats are stored as 32-bit binaries, encoded as uint32, and stringified to a snip. They are reconstructed,
//    either as textual for SAM or binary for BAM via bam_piz_special_FLOAT. Done this way so BAM binary float is reconstructed precisely.
SPECIAL_RECONSTRUCTOR (bam_piz_special_FLOAT)
{
    int64_t n;
    ASSERT (str_get_int (snip, snip_len, &n), "failed to read integer in %s", ctx->tag_name);
    
    union {
        uint32_t i;
        float f;
    } machine_en = { .i = (uint32_t)n };

    if (!reconstruct) goto finish;

    // binary reconstruction in little endian - BAM format
    if (flag.out_dt == DT_BAM) {
        uint32_t n32_lten = LTEN32 (machine_en.i); // little endian (BAM format)
        RECONSTRUCT (&n32_lten, sizeof (uint32_t)); // in binary - the float and uint32 are the same
    }

    // textual reconstruction - SAM format 
    else { 
        #define NUM_SIGNIFICANT_DIGITS 6 // 6 significant digits, as samtools does
        
        // calculate digits before and after the decimal point
        double log_f = log10 (machine_en.f >= 0 ? machine_en.f : -machine_en.f);
        unsigned int_digits = (log_f >= 0) + (unsigned)log_f;
        unsigned dec_digits = MAX_(0, NUM_SIGNIFICANT_DIGITS - int_digits);
        
        // reconstruct number with exactly NUM_SIGNIFICANT_DIGITS digits
        sprintf (BAFTtxt, "%.*f", dec_digits, machine_en.f); 
        unsigned len = strlen (BAFTtxt); 
        vb->txt_data.len += len;

        // remove trailing decimal zeros:  "5.500"->"5.5" ; "5.0000"->"5" ; "50"->"50"
        if (dec_digits) {
            unsigned trailing_zeros=0;
            for (int i=vb->txt_data.len-1; i >= vb->txt_data.len-dec_digits; i--)
                if (*Bc (vb->txt_data, i) == '0') 
                    trailing_zeros++;
                else
                    break;
            
            vb->txt_data.len -= (dec_digits==trailing_zeros) ? dec_digits+1 : trailing_zeros;
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
    if (snip_len == 1 && *snip == '*') 
        RECONSTRUCT_BIN32 (-1);

    // if its RNEXT and =, emit the last index of RNAME
    else if (ctx->dict_id.num == _SAM_RNEXT && snip_len == 1 && *snip == '=') 
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
    if (flag.out_dt == DT_FASTQ) {
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
TRANSLATOR_FUNC (sam_piz_sam2bam_ARRAY_SELF)
{
    ContainerP con = (ContainerP)recon;
    char *prefixes = &recon[con_sizeof (*con)];

    // remove the ',' from the prefix, and terminate with CON_PX_SEP_SHOW_REPEATS - this will cause
    // the number of repeats (in LTEN32) to be outputed after the prefix
    prefixes[1] = CON_PX_SEP_SHOW_REPEATS; // prefixes is now { type, CON_PX_SEP_SHOW_REPEATS }
    
    return -1; // change in prefixes length
}

//------------------------------------------------------------------------------------
// Translator and filter functions for reconstructing SAM / BAM data into FASTQ format
//------------------------------------------------------------------------------------

TXTHEADER_TRANSLATOR (txtheader_sam2fq)
{
    txtheader_buf->len = 0; // fastq has no header
}

// filtering during reconstruction: called by container_reconstruct_do for each sam alignment (repeat)
CONTAINER_CALLBACK (sam_piz_container_cb)
{
    if (is_top_level) {

        if (flag.add_line_numbers && TXT_DT(DT_SAM)) {
            vb->txt_data.len32 -= 1 + (*(BLSTtxt-1) == '\r'); // remove \n or \r\n
            vb->txt_data.len32 += sprintf (BAFTtxt, "\tVB:Z:%s\n", LN_NAME);
        }

        // case SAM to BAM translation: set alignment.block_size (was in sam_piz_sam2bam_AUX until v11)
        if (dict_id.num == _SAM_TOP2BAM) { 
            BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)Bc (vb->txt_data, vb->line_start);
            alignment->block_size = vb->txt_data.len - vb->line_start - sizeof (uint32_t); // block_size doesn't include the block_size field itself
            alignment->block_size = LTEN32 (alignment->block_size);
        }
        
        // case SAM to FASTQ translation: drop line if this is not a primary alignment (don't show its secondary or supplamentary alignments)
        else if ((dict_id.num == _SAM_TOP2FQ || dict_id.num == _SAM_TOP2FQEX)
        && ((uint16_t)vb->last_int(SAM_FLAG) & (SAM_FLAG_SECONDARY | SAM_FLAG_SUPPLEMENTARY)))
            vb->drop_curr_line = "not_primary";

        // --taxid: filter out by Kraken taxid (SAM, BAM, FASTQ)
        if (flag.kraken_taxid != TAXID_NONE && !vb->drop_curr_line 
        && (   (kraken_is_loaded  && !kraken_is_included_loaded (vb, last_txt(vb, SAM_QNAME), vb->last_txt_len (SAM_QNAME)))// +1 in case of FASTQ to skip "@"
            || (!kraken_is_loaded && !kraken_is_included_stored (vb, SAM_TAXID, !flag.collect_coverage && !flag.count)))) 
            vb->drop_curr_line = "taxid";

        // --FLAG
        if (flag.sam_flag_filter && !vb->drop_curr_line ) {

            uint16_t this_sam_flag = (uint16_t)vb->last_int (SAM_FLAG);
        
            bool all_flags_set = (this_sam_flag & flag.FLAG) == flag.FLAG;
            bool no_flags_set  = (this_sam_flag & flag.FLAG) == 0;
            if ((flag.sam_flag_filter == SAM_FLAG_INCLUDE_IF_ALL  && !all_flags_set)
            ||  (flag.sam_flag_filter == SAM_FLAG_INCLUDE_IF_NONE && !no_flags_set)
            ||  (flag.sam_flag_filter == SAM_FLAG_EXCLUDE_IF_ALL  &&  all_flags_set))
                vb->drop_curr_line = "FLAG";
        }

        // --MAPQ
        if (flag.sam_mapq_filter && !vb->drop_curr_line) {
            
            if (dict_id.num == _SAM_TOP2FQ || dict_id.num == _SAM_TOP2FQEX)
                reconstruct_from_ctx (vb, SAM_MAPQ, 0, false); // when translating to FASTQ, MAPQ is normally not reconstructed

            uint8_t this_mapq = (uint8_t)vb->last_int (SAM_MAPQ);
        
            if ((flag.sam_mapq_filter == SAM_MAPQ_INCLUDE_IF_AT_LEAST && this_mapq < flag.MAPQ) ||
                (flag.sam_mapq_filter == SAM_MAPQ_EXCLUDE_IF_AT_LEAST && this_mapq >= flag.MAPQ))
                vb->drop_curr_line = "MAPQ";
        }

        // --bases
        if (flag.bases && !vb->drop_curr_line && 
            !(TXT_DT(DT_BAM) ? iupac_is_included_bam   (last_txt (vb, SAM_SQBITMAP), ((BAMAlignmentFixed *)recon)->l_seq)
                             : iupac_is_included_ascii (last_txt (vb, SAM_SQBITMAP), vb->last_txt_len (SAM_SQBITMAP))))
            vb->drop_curr_line = "bases";
        
        // count coverage, if needed    
        if ((flag.show_sex || flag.show_coverage) && is_top_level && !vb->drop_curr_line)
            sam_piz_update_coverage (vb, vb->last_int(SAM_FLAG), VB_SAM->soft_clip[0] + VB_SAM->soft_clip[1]);

        if (flag.idxstats && !vb->drop_curr_line) {
            if (vb->last_int(SAM_FLAG) & SAM_FLAG_UNMAPPED)   
                (*B64 (vb->unmapped_read_count, vb->last_index(SAM_RNAME)))++;
            else
                (*B64 (vb->read_count, vb->last_index(SAM_RNAME)))++;
        }
    }
}

// filter is called before reconstruction of a repeat or an item, and returns false if item should 
// not be reconstructed. contexts are not consumed.
CONTAINER_FILTER_FUNC (sam_piz_filter)
{
    // BAM: set buddy at the beginning of each line, as QNAME is reconstructed much later
    // For MAIN and DEPN: buddy, if existing, will have QNAME=SNIP_COPY_BUDDY (see sam_seg_QNAME)
    // note: we always load buddy, to prevent a situation when in some lines it is consumed
    // and other lines, which have buddy, it is not consumed because the field that consumes it is skipped
    if (!sam_is_prim_vb && dict_id.num == _SAM_TOP2BAM && item == 0 && 
        CTX(SAM_QNAME)->b250.len32) { // might be 0 in special cases, like flag.count
        STR(snip);
        PEEK_SNIP(SAM_QNAME);
        if (snip_len && *snip == SNIP_COPY_BUDDY)
            reconstruct_set_buddy(vb);
    }
 
    // BAM of a primary VB: sets sag, and then sets buddy if indicated by the sag.
    // note: it will be so, if SAM_QNAME which was loaded to the sag, had a buddy.
    else if (sam_is_prim_vb && dict_id.num == _SAM_TOP2BAM && item == 0)
        sam_piz_set_sag (VB_SAM); 

    // collect_coverage: set buddy_line_i here, since we don't reconstruct QNAME
    // note: we always load buddy, to prevent a situation when in some lines it is consumed
    // and other lines, which have buddy, it is not consumed because the field that consumes it is skipped
    else if (dict_id.num == _SAM_QNAME && flag.collect_coverage) {
        STR(snip);
        LOAD_SNIP(SAM_QNAME);
        if (snip_len && *snip == SNIP_COPY_BUDDY)
            reconstruct_set_buddy(vb);
        return false; // don't reconstruct QNAME
    }

    // XA: insert RNAME, STRAND and POS lookbacks
    else if (dict_id.num == _OPTION_XA_Z && item == 4)  // importantly, after RNAME and STRAND_POS
        sam_piz_XA_field_insert_lookback (vb);
    
    // collect_coverage: rather than reconstructing optional, reconstruct SAM_MC_Z that just consumes MC:Z if it exists
    else if (dict_id.num == _SAM_AUX && flag.collect_coverage) {
        reconstruct_from_ctx (vb, SAM_MC_Z, 0, false);
        return false; // don't reconstruct AUX
    }

    return true; // go ahead and reconstruct
}

// emit 1 if (FLAGS & 0x40) or 2 of (FLAGS & 0x80) except if --FLAG is specified too (--FLAG can be
// used to get only R1 or R2)
TRANSLATOR_FUNC (sam_piz_sam2fastq_FLAG)
{
    if (flag.sam_flag_filter) return 0;
    
    uint16_t sam_flag = (uint16_t)vb->last_int(SAM_FLAG);

    if (sam_flag & (SAM_FLAG_IS_FIRST | SAM_FLAG_IS_LAST)) {
        
        recon -= item_prefix_len + 1; // move to before prefix and previous field's separator (\n for regular fq or \t for extended)
        memmove (recon+2, recon, recon_len + item_prefix_len + 1); // make room for /1 or /2 
        recon[0] = '/';
        recon[1] = (sam_flag & SAM_FLAG_IS_FIRST) ? '1' : '2';
        vb->txt_data.len += 2; 
    }

    return 0;
}

// v14: De-multiplex by has_mate
SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_BY_MATE)
{
    return reconstruct_demultiplex (vb, ctx, STRa(snip), piz_has_mate, new_value, reconstruct);
}

// v14: De-multiplex by has_buddy (which can be mate or prim)
SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_BY_BUDDY)
{
    return reconstruct_demultiplex (vb, ctx, STRa(snip), piz_has_buddy, new_value, reconstruct);
}

// v14: De-multiplex by has_mate and has_prim
SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_BY_MATE_PRIM)
{
    // note: when reconstructing POS in BAM, FLAG is not known yet, so we peek it here
    if (!ctx_has_value_in_line_(vb, CTX(SAM_FLAG)))
        ctx_set_last_value (vb, CTX(SAM_FLAG), reconstruct_peek (vb, CTX(SAM_FLAG), 0, 0));

    int channel_i = piz_has_mate?1 : piz_has_real_prim?2 : 0;

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

