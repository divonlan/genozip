// ------------------------------------------------------------------
//   sam_sag_zip.c
//   Copyright (C) 2022-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "sam_private.h"
#include "chrom.h"
#include "writer.h"
#include "compressor.h"
#include "qname.h"
#include "zfile.h"
#include "huffman.h"
#include "sorter.h"

typedef struct {
    VBIType vb_i;              // MAIN VB vblock_i
    uint64_t first_gc_line[2]; // first PRIM or DEPN gencomp line that originates from this MAIN VB - 0-based line within the PRIM or DEPN component
    uint32_t num_gc_lines[2];  // number of PRIM or DEPN gencomp lines that originate from this MAIN VB
} SamMainVbInfo;

typedef struct {
    VBIType vb_i;              // This PRIM or DEPN VB vblock_i
    uint64_t first_line;       // 0-based line within this comp
    uint32_t num_lines;
    uint32_t lines_consumed;
} SamGcVbInfo;

const SoloProp solo_props[NUM_SOLO_TAGS] = SOLO_PROPS;

// called by main thread after reading a VB - called from sam_zip_init_vb
void sam_sag_zip_init_vb (VBlockSAMP vb)
{
    if (IS_MAIN(vb) || !Ltxt) return;

    // PRIM or DEPN - add to vb_info
    BufferP vb_info = &z_file->vb_info[vb->comp_i]; 
    
    SamGcVbInfo new_info = {
        .vb_i       = vb->vblock_i,
        .first_line = (vb_info->len == 0) ? 0 : (BLST(SamGcVbInfo, *vb_info)->first_line + BLST(SamGcVbInfo, *vb_info)->num_lines),
        .num_lines  = vb->lines.len
    };

    buf_append_one (z_file->vb_info[vb->comp_i], new_info);

    if (flag_debug_gencomp)
        iprintf ("SetVbInfo %s: first_line=%"PRIu64" num_lines=%u\n", VB_NAME, new_info.first_line, new_info.num_lines);
}

// called by compute thread holding mutex, after absorbing a MAIN VB (in order of aborption)
void sam_add_main_vb_info (VBlockP vb,
                           uint64_t prim_first_line, uint32_t prim_num_lines, 
                           uint64_t depn_first_line, uint32_t depn_num_lines)
{    
    VB_SAM->main_vb_info_i = z_file->vb_info[SAM_COMP_MAIN].len32;

    SamMainVbInfo new_info = {
        .vb_i             = vb->vblock_i,
        .first_gc_line[0] = prim_first_line, .num_gc_lines[0] = prim_num_lines,
        .first_gc_line[1] = depn_first_line, .num_gc_lines[1] = depn_num_lines,
    };

    // Note: we are adding main in VBs in order of their aborption, 
    // therefore prim_first_line and depn_first_line are each in order (i.e. monotonically increasing)
    buf_append_one (z_file->vb_info[SAM_COMP_MAIN], new_info);

    if (flag_debug_gencomp)
        iprintf ("SetVbInfo %s: prim_lines={first=%"PRIu64" n=%u} depn_lines={first=%"PRIu64" n=%u}\n", 
                 VB_NAME, new_info.first_gc_line[0], new_info.num_gc_lines[0], new_info.first_gc_line[1], new_info.num_gc_lines[1]);
}                           

// main thread: called from sam_zip_after_segconf
void sam_set_sag_type (void)
{
    // --no-gencomp prohibits gencomp
    if (flag.no_gencomp)
        segconf.sag_type = SAG_NONE;

    // --fast, --low-memory or a non-sorted file prohibit gencomp, but may be overrided with --force-gencomp
    else if ((!segconf.is_sorted || flag.fast || flag.low_memory || flag_has_head) && !flag.force_gencomp)
        segconf.sag_type = SAG_NONE;

    else if (MP(LONGRANGER))
        segconf.sag_type = SAG_NONE; // TO DO - new SAG_BY_LONGRANGER SamGcVbInfo : longranger has non-standard SA:Z format

    else if (MP(NOVOALIGN))
        segconf.sag_type = SAG_BY_NH; // NovoAlign may have both NH and SA, we go by NH

    // SAG_BY_NH: STAR and other aligners with NH:i/HI:i for secondary alignments.
    // Note that STAR uses NH:i for secondary and SA:Z for supplamentary. We choose to gencomp for
    // secondaries as they are usually more numerous that supplamentaries (bug 1148).
    // we identify a PRIM as a line that is (1) not supp/secondary (2) has NH >= 2
    else if (segconf.HI_has_two_plus && segconf.has[OPTION_NH_i])
        segconf.sag_type = segconf.has_barcodes ? SAG_BY_SOLO : SAG_BY_NH; 

    // SAG_BY_SA: the preferd Sag method - by SA:Z
    // we identify a PRIM line - if it (1) has SA:Z and (2) supp/secondary flags are clear
    else if (segconf.sam_has_SA_Z || segconf.has[OPTION_SA_Z]) 
        segconf.sag_type = SAG_BY_SA;

    // SAG_BY_CC: cases where there is NH, but no HI or SA, and all lines in a SAG have NH>=2, and, 
    // except the last line in the SAG, also have CC and CP.
    // Note: we don't test for CC/CP here, because only lines with NH>=2 have them, and that might not have been encountered yet
    else if (segconf.has[OPTION_NH_i] && !segconf.has[OPTION_HI_i] && !segconf.has[OPTION_SA_Z])
        segconf.sag_type = SAG_BY_CC;

    // SAG_BY_FLAG: cases (eg BLASR) where there might be dependent lines, but no auxilliary fields to support them
    // We treat all non-depn mapped lines as PRIM (and warn the user of memory consumption)
    else if ((MP(BLASR) || segconf.sam_has_depn || flag.force_gencomp) && 
             (flag.best || flag.force_gencomp) &&  // too slow for normal mode
             !txt_file->redirected && !txt_file->is_remote) { // conditions for using sam_sag_by_flag_scan_for_depn
        sam_sag_by_flag_scan_for_depn ();
        segconf.sag_type = !z_file->sag_depn_index.len ? SAG_NONE
                         : segconf.has[OPTION_SA_Z]    ? SAG_BY_SA  // scan can detect SA:Z not previously detected by segconf
                         :                               SAG_BY_FLAG;
    }

    else 
        segconf.sag_type = SAG_NONE;      

    if (segconf.sag_type) {
        sam_sa_prim_initialize_ingest(); // the PRIM component is compressed (out-of-band) at the same time as MAIN
        gencomp_initialize (SAM_COMP_PRIM, GCT_OOB); 
        gencomp_initialize (SAM_COMP_DEPN, GCT_DEPN); 

        buf_alloc (evb, &z_file->vb_info[SAM_COMP_MAIN], 0, 10000, SamMainVbInfo, 0, "z_file->vb_info");
        buf_alloc (evb, &z_file->vb_info[SAM_COMP_PRIM], 0, 3000,  SamGcVbInfo,   0, "z_file->vb_info");
        buf_alloc (evb, &z_file->vb_info[SAM_COMP_DEPN], 0, 3000,  SamGcVbInfo,   0, "z_file->vb_info");
    }
}

void sam_seg_gc_initialize (VBlockSAMP vb)
{
    // DEPN stuff
    if (IS_DEPN(vb)) {
        CTX(SAM_SAG)->ltype = sizeof(SAGroup)==4 ? LT_UINT32 : LT_UINT64;
        CTX(SAM_SAALN)->ltype   = LT_UINT16; // index of alignment with SA Group
    }

    // PRIM stuff 
    else if (IS_PRIM(vb)) { 
        CTX(SAM_SAALN)->ltype   = LT_UINT16; // index of alignment with SA Group
        CTX(SAM_FLAG)->no_stons = true;        

        ctx_set_dyn_int (VB, OPTION_SA_Z, DID_EOL);    // we store num_alns in local

        ctx_set_store (VB, STORE_INDEX, OPTION_SA_RNAME, OPTION_SA_STRAND, DID_EOL);
        ctx_set_store (VB, STORE_INT, OPTION_SA_POS, OPTION_SA_MAPQ, OPTION_SA_NM, OPTION_NH_i, DID_EOL);
    }
}

// Note: in MAIN vbs, failure means we will keep this line in MAIN and not move it to PRIM
// In PRIM, we abort because this is not not expected as MAIN component should not have moved this line to PRIM
#define FAILIF(condition, format, ...) \
  ( { if (condition) { \
          if (flag.debug_sag) iprintf ("%s: " format "\n", LN_NAME, __VA_ARGS__); \
          if (IS_MAIN(vb)) return false; \
          else { progress_newline(); fprintf (stderr, "%s: Error in %s:%u: Failed PRIM line because ", LN_NAME, __FUNCLINE); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "%s", report_support_if_unexpected()); fflush (stderr); exit_on_error(true); }} \
    } )

// Call in seg PRIM line (both in MAIN and PRIM vb): 
// MAIN: test that line is a valid prim before moving it to the PRIM gencomp file
// PRIM: add SA Group (not alignments) to the data structure in VB 
bool sam_seg_prim_add_sag (VBlockSAMP vb, ZipDataLineSAMP dl, uint16_t num_alns/* inc primary aln*/, bool is_bam)
{
    rom seq = vb->textual_seq_str;
    uint32_t seq_len = dl->SEQ.len;

    ASSERT (seq, "%s: seq=NULL", LN_NAME);

    // MAIN: check that values are within limits defined in Sag (no need to check in PRIM as we already checked in MAIN)
    if (IS_MAIN(vb)) {
        FAILIF (!sam_might_have_saggies_in_other_VBs (vb, dl, num_alns), "all sag alignments are contained in this VB%s", "");
        FAILIF (dl->QNAME_len > SAM_MAX_QNAME_LEN, "dl->QNAME_len=%u > %u", dl->QNAME_len, SAM_MAX_QNAME_LEN);
        FAILIF (seq_len==1 && *seq == '*', "SEQ=\"*\"%s", ""); // we haven't segged seq yet, so vb->seq_missing is not yet set
        FAILIF (seq_len > MAX_SA_SEQ_LEN, "seq_len=%u > %u", seq_len, MAX_SA_SEQ_LEN);
        FAILIF (!dl->POS, "POS=0%s", ""); // unaligned
        FAILIF (num_alns > MAX_SA_NUM_ALNS, "%s=%u > %u", (IS_SAG_NH || IS_SAG_SOLO || IS_SAG_CC) ? "NH:i" : "num_alns", num_alns, MAX_SA_NUM_ALNS);
        FAILIF (dl->hard_clip[0] || dl->hard_clip[1], "has_hard_clips%s", ""); // primary cannot have hard clips
        
        uint32_t bad_i;
        FAILIF (!str_is_ACGT (STRa(seq), &bad_i), 
                "RNAME=\"%.*s\" POS=%d. Found '%c' in base_i=%u in SEQ is not A,C,G or T. SEQ=\"%.*s\"", 
                STRf(vb->chrom_name), dl->POS, seq[bad_i], bad_i, STRf(seq));
    }

    // PRIM: actually add the group
    else if (IS_PRIM(vb)) { 
        buf_alloc (vb, &vb->sag_grps, 1, 1000, Sag, CTX_GROWTH, "sag_grps");

        // add SA group
        BNXT (Sag, vb->sag_grps) = (Sag){
            .first_aln_i = vb->sag_alns.len,
            .num_alns    = num_alns,
            .qname       = BNUMtxt (dl_qname(dl)),
            .qname_len   = dl->QNAME_len,
            .qual        = dl->QUAL.index,
            .no_qual     = vb->qual_missing,
            .revcomp     = dl->FLAG.rev_comp,
            .multi_segs  = dl->FLAG.multi_segs,
            .is_first    = dl->FLAG.is_first,
            .is_last     = dl->FLAG.is_last,
            .seq         = dl->SEQ.index,
            .seq_len     = seq_len,
            .as          = CAP_SA_AS (dl->AS) // [0,65535] capped at 0 and 65535
        };
    }

    return true;
}

// MAIN: called in sam_seg_is_gc_line to test that line is a valid prim before moving it to the PRIM gencomp file
// PRIM: call in seg of seg of SA:Z add SA Group (including alignments) to the data structure in VB
// return n_alns if successful and 0 if not
int32_t sam_seg_prim_add_sag_SA (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(sa), int64_t this_nm, bool is_bam)
{
    str_split (sa, sa_len, 0, ';', aln, false); // n_alns is 1+number of repeats in SA, because of terminating ';'. so it is the total number of alignments, including the primary
    FAILIF (n_alns < 2 || n_alns > MAX_SA_NUM_ALNS, "n_alns=%d is not in [2,%u]", n_alns, MAX_SA_NUM_ALNS); // we cannot add this SA - either invalid or over the max number of alignments

    if (!sam_seg_prim_add_sag (vb, dl, n_alns, is_bam)) return 0; 

    // check that values are within limits defined in SAAln
    FAILIF (dl->POS  > MAX_SA_POS, "POS=%d > MAX_SA_POS=%u", dl->POS, MAX_SA_POS);
    FAILIF (dl->MAPQ > MAX_SA_MAPQ, "MAPQ=%u > MAX_SA_MAPQ=%u", dl->MAPQ, MAX_SA_MAPQ);
    FAILIF (this_nm  > MAX_SA_NM, "NM=%"PRId64" > MAX_SA_NM=%u", this_nm, MAX_SA_NM);
    FAILIF (vb->cigar_missing, "No CIGAR%s", "");
    FAILIF (vb->hard_clip[0], "hard_clip[LEFT]=%u > 0", vb->hard_clip[0]); // we can't add primary lines with a CIGAR with H as we need the full sequence and qual 
    FAILIF (vb->hard_clip[1], "hard_clip[RIGHT]=%u > 0", vb->hard_clip[1]);

    if (IS_PRIM(vb)) { // in PRIM, we actually add the group, in MAIN, we are just testing
        buf_alloc (vb, &vb->sag_alns, n_alns /*+1 for primary aln*/, 64, SAAln, CTX_GROWTH, "sag_alns");

        // in Seg, we add the Adler32 of CIGAR to save memory, as we just need to verify it, not reconstruct it
        uint32_t textual_cigar_len = (is_bam ? vb->textual_cigar.len32 : dl->CIGAR.len);
        rom textual_cigar = is_bam ? B1STc(vb->textual_cigar) : Btxt (dl->CIGAR.index);
        
        FAILIF (textual_cigar_len > MAX_SA_CIGAR_LEN, "CIGAR.len=%u > MAX_SA_CIGAR_LEN=%u", textual_cigar_len, MAX_SA_CIGAR_LEN);

        // add primary alignment as first alignment in SA group
        BNXT (SAAln, vb->sag_alns) = (SAAln){
            .rname           = vb->chrom_node_index,
            .pos             = dl->POS,
            .revcomp         = dl->FLAG.rev_comp,
            .cigar.signature = cigar_sign (vb, dl, STRa(textual_cigar)),
            .mapq            = dl->MAPQ,
            .nm              = this_nm
        };

        // chew non-dict cigar if this is the first prim vb (for use in piz - only non-dict cigars get nico-compressed)
        if (textual_cigar_len > MAX_CIGAR_LEN_IN_DICT && z_file->SA_CIGAR_chewing_vb_i == vb->vblock_i)
            nico_chew_one_cigar (OPTION_SA_CIGAR, CIG(vb->binary_cigar));

        vb->sag_alns.count += n_alns;
    }

    // test or add all dependent alignments
    for (uint32_t i=0; i < n_alns-1; i++) {
        
        // get items - rname, pos, revcomp, CIGAR, mapQ, NM
        str_split (alns[i], aln_lens[i], NUM_SA_ITEMS, ',', item, true);
        FAILIF (n_items != NUM_SA_ITEMS, "in SA alignment %u - n_items=%u != NUM_SA_ITEMS=%u", i, n_items, NUM_SA_ITEMS);

        FAILIF (aln_lens[i] > MAX_SA_CIGAR_LEN, "SA_CIGAR.len=%u > MAX_SA_CIGAR_LEN=%u", aln_lens[i], MAX_SA_CIGAR_LEN);

        // rname (pre-populated from sam header)
        WordIndex rname = ctx_get_ol_node_index_by_snip (VB, CTX(SAM_RNAME), STRi(item,SA_RNAME));
        FAILIF (rname == WORD_INDEX_NONE, "in SA alignment %u - rname==WORD_INDEX_NONE", i);

        // protect rname from removal by ctx_shorten_unused_dict_words if it is unused in the main RNAME field. 
        ctx_protect_from_removal (VB, CTX(SAM_RNAME), rname); 

        // pos, mapq, nm - get integers and verify limits
        int64_t pos, mapq, nm;
        FAILIF (!str_get_int_range64 (STRi(item,SA_POS ), 0, MAX_SA_POS,  &pos ), "in SA alignment %u - pos=%"PRId64" ∉ [0, %u]", i, pos, MAX_SA_POS);
        FAILIF (!str_get_int_range64 (STRi(item,SA_MAPQ), 0, MAX_SA_MAPQ, &mapq), "in SA alignment %u - mapq=%"PRId64" ∉ [0, %u]", i, mapq, MAX_SA_MAPQ);
        FAILIF (!str_get_int_range64 (STRi(item,SA_NM  ), 0, MAX_SA_NM,   &nm  ), "in SA alignment %u - nm=%"PRId64" ∉ [0, %u]", i, nm, MAX_SA_NM);
 
        // revcomp
        bool revcomp = *items[SA_STRAND] == '-';
        FAILIF (item_lens[SA_STRAND] != 1 || (!revcomp && *items[SA_STRAND] != '+'), "in SA alignment %u - STRAND is not '+' or '-'", i);

        if (IS_PRIM(vb)) { // in PRIM, we actually add the group, in MAIN, we are just testing

            if (!segconf.SA_CIGAR_can_have_H && sam_cigar_has_H (STRi(item,SA_CIGAR)))
                segconf.SA_CIGAR_can_have_H = true; // no worries about thread safety, this just gets set to true while segging MAIN and is inspected only when segging MAIN is over

            uint64_t cigar_sig = cigar_sign (vb, NULL, STRi(item,SA_CIGAR));

            BNXT (SAAln, vb->sag_alns) = (SAAln){ 
                .rname           = rname, 
                .pos             = pos, 
                .revcomp         = revcomp, 
                .cigar.signature = cigar_sig,
                .mapq            = mapq,
                .nm              = nm
            };

            if (z_file->SA_CIGAR_chewing_vb_i == vb->vblock_i)
                nico_chew_one_textual_cigar (VB, OPTION_SA_CIGAR, STRi(item,SA_CIGAR));
        }
    }

    return n_alns;
}       

// PRIM: call in seg of seg of NH:i add SA Group to the data structure in VB 
void sam_seg_prim_add_sag_NH (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t nh)
{
    sam_seg_prim_add_sag (vb, dl, nh, IS_BAM_ZIP);
    vb->sag_alns.count += nh;
}

// PRIM: call in seg of seg of NH:i add SA Group to the data structure in VB 
void sam_seg_prim_add_sag_CC (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t nh)
{
    sam_seg_prim_add_sag (vb, dl, nh, IS_BAM_ZIP);

    buf_alloc (vb, &vb->sag_alns, 1, 64, CCAln, CTX_GROWTH, "sag_alns");
    BNXT (CCAln, vb->sag_alns) = (CCAln){ 
        // we don't enter new RNAMES (not in the sam header / reference file / committed by previous VBs) as we don't have a word_index yet.
        .rname = (dl->RNAME < CTX(SAM_RNAME)->ol_nodes.len32) ? dl->RNAME : WORD_INDEX_NONE, 
        .pos = dl->POS 
    };

    vb->sag_alns.count += nh;
}

// PRIM: called to ingest a segged VB data into solo "alignments" 
void sam_seg_prim_add_sag_SOLO (VBlockSAMP vb, ZipDataLineSAMP dl)
{
    sam_seg_prim_add_sag (vb, dl, dl->NH, IS_BAM_ZIP);

    // TO DO: consider the tradeoff - this way, we don't allocate memory and add in each line (like sam_seg_prim_add_sag_CC)
    // BUT, we do hundreds of thousands of small memcpy's during ingest - with VBs serialized on a mutex, instead of just
    // one big memcpy
}

// SA:Z was encountered for the first time (perhaps multiple threads concurrently) and abbreviation is possible
static void sam_seg_MAIN_determine_abbreviated (VBlockSAMP vb, bool is_bam)
{
    // one of the threads gets the mutex 
    mutex_lock (z_file->test_abbrev_mutex);

    // case: not already determined by another concurrent thread
    if (segconf.SA_CIGAR_abbreviated == unknown) 
        segconf.SA_CIGAR_abbreviated = sam_test_SA_CIGAR_abbreviated (STRauxZ (SA_Z, is_bam));

    mutex_unlock (z_file->test_abbrev_mutex);
}

// Seg MAIN: returns true if this line is a Primary or Dependent line of a supplementary/secondary group - 
// and should be moved to a generated component. Called by compute thread in seg of Normal VB.
bool sam_seg_is_gc_line (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(alignment), bool is_bam)
{
    START_TIMER;

    SamNMType NM;
    int32_t n_alns; 
    SamComponentType comp_i = COMP_MAIN; // generated component
    
    // unmapped or otherwise lacking alignment - not gencomp
    if (!segconf.sag_type || vb->chrom_node_index == WORD_INDEX_NONE || // RNAME='*' or headerless SAM (we need header-contigs for sam_sa_add_sa_group)
        dl->FLAG.unmapped || vb->cigar_missing || !dl->POS)
            goto done; 

    bool has_hards = (vb->hard_clip[0] > 0 || vb->hard_clip[1] > 0);
    bool has_softs = (vb->soft_clip[0] > 0 || vb->soft_clip[1] > 0);
    
    switch (segconf.sag_type) {
        case SAG_BY_SA:
            if (has_hards && has_softs)
                goto done; // we don't support adding an alignment with both soft and hard clips to an SA-based sag

            if (has(SA_Z) && segconf.SA_CIGAR_abbreviated == unknown)
                sam_seg_MAIN_determine_abbreviated (vb, is_bam); // SA_CIGAR_abbreviated was not seg in segconf, we set it now

            if (has(SA_Z) && sam_line_is_depn(dl)) {
                comp_i = SAM_COMP_DEPN;

                if (vb->qname_count.len32)
                    n_alns = str_count_char (STRauxZ (SA_Z, is_bam), ';') + 1; // +1 for this aln
            
                if (has_hards)
                   segconf.depn_CIGAR_can_have_H = true; // no worries about thread safety, this just gets set to true while segging MAIN and is inspected only when segging MAIN is over
            }

            // case: we determine the nm subfield is number of X bases in CIGAR  
            else if (segconf.SA_NM_by_CIGAR_X && has(SA_Z) &&
                     (n_alns = sam_seg_prim_add_sag_SA (vb, dl, STRauxZ (SA_Z, is_bam), vb->mismatch_bases_by_CIGAR, is_bam)))  // testing to see if we can successfully add a sag based on SA
                comp_i = SAM_COMP_PRIM;

            // case: standard SA:Z            
            else if (has(SA_Z) && has(NM_i) &&
                     sam_seg_get_aux_int (vb, vb->idx_NM_i, &NM, is_bam, MIN_NM_i, MAX_NM_i, SOFT_FAIL) &&
                     (n_alns = sam_seg_prim_add_sag_SA (vb, dl, STRauxZ (SA_Z, is_bam), NM, is_bam)))  // testing to see if we can successfully add a sag based on SA
                comp_i = SAM_COMP_PRIM;

            break;
        
        case SAG_BY_NH: 
        case SAG_BY_SOLO: 
            // case: IH=1 - no gencomp
            if (sam_seg_peek_int_field (vb, OPTION_IH_i, vb->idx_IH_i, 1, 1, false, NULL))
                goto done;

            if (has(NH_i) && sam_seg_get_aux_int (vb, vb->idx_NH_i, &n_alns, is_bam, 2/*at least*/, MAX_HI_NH, SOFT_FAIL)) {
                
                if (sam_line_is_depn(dl)) 
                    comp_i = SAM_COMP_DEPN;

                else if (sam_seg_prim_add_sag (vb, dl, n_alns, is_bam)) // testing to see if we can successfully add a sag based on NH
                    comp_i = SAM_COMP_PRIM;
            }
            break;
            
        case SAG_BY_CC:
            if (has(NH_i) && sam_seg_get_aux_int (vb, vb->idx_NH_i, &n_alns, is_bam, 2, MAX_HI_NH, SOFT_FAIL)) {  // not out of range, i.e. at least 2

                if (has(CC_Z) && has(CP_i))
                    comp_i = SAM_COMP_DEPN;

                else if (!has(CC_Z) && !has(CP_i) && !has(SA_Z) && !has(HI_i) && !has_hards && 
                         sam_seg_prim_add_sag (vb, dl, n_alns, is_bam)) // testing to see if we can successfully add a sag based on NH
                    comp_i = SAM_COMP_PRIM;
                }
            break;

        case SAG_BY_FLAG: // note: we only possibly set this for BAM (not SAM) files
            if (sam_line_is_depn(dl))
                comp_i = SAM_COMP_DEPN;

            else if (sam_seg_prim_add_sag (vb, dl, 0, is_bam))
                comp_i = SAM_COMP_PRIM;
            
            // note: if the line is prim, is stays in the VB. if its MAIN (see above), it moves to a gencomp.
            break;
        
        default: ABORT ("Invalid sag_type=%u", segconf.sag_type);
    }

    // check if all alignments of this sag are contained in this VB - in which case we will 
    // seg them vs saggy_line_i rather than using gencomp (note: for SAM_COMP_PRIM, we already tested in sam_seg_prim_add_sag)
    if (comp_i == SAM_COMP_DEPN && !sam_might_have_saggies_in_other_VBs (vb, dl, n_alns)) 
        comp_i = COMP_MAIN;

    if (comp_i != COMP_MAIN) {

        // store location where this gc line should be inserted    
        gencomp_seg_add_line (VB, comp_i, STRa(alignment));
        
        vb->line_i--;
        vb->recon_size -= alignment_len;
        vb->txt_size   -= alignment_len;
        vb->num_gc_lines++;

        if      (comp_i == SAM_COMP_PRIM) vb->seg_found_prim_line = true;
        else if (comp_i == SAM_COMP_DEPN) vb->seg_found_depn_line = true;
    }

done:        
    #define last_plan BLST(VbPlanItem, vb->vb_plan)

    if (vb->vb_plan.len32 && vb->vb_plan.prev_comp_i == comp_i && last_plan->n_lines < MAX_PLAN_ITEM_LINES)
        last_plan->n_lines++;
    
    else {
         // case: we need a 0-lines transition entry to modify prev_comp_i to be different from comp_i
        if (comp_i == vb->vb_plan.prev_comp_i) {
            CompIType transition_comp_i = (comp_i==SAM_COMP_MAIN ? SAM_COMP_PRIM : SAM_COMP_MAIN);
            int16_t transition_comp = VB_PLAN_COMP_ZIP (transition_comp_i);
            buf_append_one (vb->vb_plan, ((VbPlanItem){ .comp=transition_comp, .n_lines = 0 })); 
            
            vb->vb_plan.prev_comp_i = transition_comp_i;
        }

        int16_t comp = VB_PLAN_COMP_ZIP(comp_i);
        buf_append_one (vb->vb_plan, ((VbPlanItem){ .comp=comp, .n_lines = 1})); // initial allocation in sam_seg_initialize
        vb->vb_plan.prev_comp_i = comp_i;
    }

    COPY_TIMER (sam_seg_is_gc_line);
    return comp_i != COMP_MAIN; // true if line moves to generated component
}

// ------------------
// QNAME
// ------------------

typedef struct { uint32_t qname_hash, grp_i; } SAGroupIndexEntry; 

static BINARY_SEARCHER (sam_sa_binary_search_for_qname_hash, SAGroupIndexEntry, int64_t, qname_hash, false, IfNotExact_ReturnNULL);

// ZIP DEPN: find group index with this_qname in z_file->sag_grps, and if there are several - return the first
static Sag *sam_sa_get_first_group_by_qname_hash (VBlockSAMP vb, STRp(this_qname), bool is_last, int64_t *grp_index_i, 
                                                  uint32_t *this_qname_hash) // out
{
    // search for a group with qname in z_file->sa_qname
    *this_qname_hash = qname_calc_hash (QNAME1, COMP_NONE, this_qname, this_qname_len, is_last, false, CRC32, NULL);

    // search for index entry by qname_hash 
    SAGroupIndexEntry *index_entry = binary_search (sam_sa_binary_search_for_qname_hash, SAGroupIndexEntry, z_file->sag_grps_index, *this_qname_hash);
    *grp_index_i = index_entry ? BNUM (z_file->sag_grps_index, index_entry) : -1;

    const SAGroupIndexEntry *index_ent = B(SAGroupIndexEntry, z_file->sag_grps_index, *grp_index_i); // invalid pointer if grp_index_i==-1, that's ok    
    return (*grp_index_i >= 0) ? B(Sag, z_file->sag_grps, index_ent->grp_i) : NULL; 
}

// ZIP DEPN: if there are more groups in z_file->sag_grps with the qname adlers, return the next group index
static Sag *sam_sa_get_next_group_by_qname_hash (VBlockSAMP vb, int64_t *grp_index_i)
{
    const SAGroupIndexEntry *index_ent = B(SAGroupIndexEntry, z_file->sag_grps_index, *grp_index_i);

    if (*grp_index_i < z_file->sag_grps_index.len-1 && index_ent->qname_hash == (index_ent+1)->qname_hash) {
        (*grp_index_i)++;
        return B(Sag, z_file->sag_grps, (index_ent+1)->grp_i);
    }

    return NULL; // not found
}

// Seg DEPN: check if the seq from this DEPN line matches a SA Group, considering hard-clips in its CIGAR
// returns number of mismatches
static bool sam_seg_depn_is_subseq_of_prim (VBlockSAMP vb, bytes depn_textual_seq, uint32_t depn_seq_len,
                                            bool xstrand, // depn_seq is opposite strand vs primary
                                            const Sag *g, bool is_bam)
{
    Bits *prim_seq = (BitsP)&z_file->sag_seq; // ACGT format
    uint64_t start_p = g->seq + vb->hard_clip[0];

    if (!xstrand)
        for (uint64_t i=0; i < depn_seq_len; i++) {
            uint8_t p = bits_get2 (prim_seq, (start_p + i) * 2);

            switch (depn_textual_seq[i]) {
                case 'A': if (p != 0) return false; break; // if prim base is not 'A'(=0) then the sequences don't match
                case 'C': if (p != 1) return false; break; 
                case 'G': if (p != 2) return false; break; 
                case 'T': if (p != 3) return false; break; 
                default :             return false; // depn seq has a IUPAC "base" - we can't seg it against prim
            }
        }

    else {
        uint64_t start_p = g->seq + vb->hard_clip[1];

        for (uint64_t i=0; i < depn_seq_len; i++) {
            uint8_t p = bits_get2 (prim_seq, (start_p + i) * 2);

            switch (depn_textual_seq[depn_seq_len - i - 1]) {
                case 'T': if (p != 0) return false; break; // since xstrand depn seq is T, we expect prim seq to be A(=0)
                case 'G': if (p != 1) return false; break; 
                case 'C': if (p != 2) return false; break; 
                case 'A': if (p != 3) return false; break; 
                default :             return false; // depn seq has a IUPAC "base" - we can't seg it against prim
            }
        }
    }

    return true; 
}

// --------------------------------------------------
// Seg Alignments (main fields + SA-field alignments)
// --------------------------------------------------

// Seg DEPN line: verify that this line is a member of group (identical alignments) + return aln_i of this line
static inline bool sam_seg_depn_find_SA_aln (VBlockSAMP vb, const Sag *g,
                                             uint32_t n_my_alns, const SAAln *my_alns,
                                             uint16_t *my_aln_i, // out - relative to group
                                             uint16_t *prim_aln_index_in_SA_Z) // out: index of PRIM alignment in my SA:Z (1-based)
{
    const SAAln *grp_alns = B(SAAln, z_file->sag_alns, g->first_aln_i);
    *prim_aln_index_in_SA_Z = 0;

    // the group's primary alignmet (grp_alns[0]) usually appears in SA[0] in the depns lines, but
    // not always (e.g. in pbmm2, it may appear in another position - but a consistent one across all depn alignments)
    for (uint32_t my_SA_aln_i=1; my_SA_aln_i < n_my_alns; my_SA_aln_i++) // iterate on this depn line's SA:Z alignments (hence starting from 1, as [0] is this line main-field alignment)
        if (!memcmp (&grp_alns[0], &my_alns[my_SA_aln_i], sizeof(SAAln))) {
            *prim_aln_index_in_SA_Z = my_SA_aln_i; // most commonly this is 1 (except for pbmm2)
            break;
        }

    if (*prim_aln_index_in_SA_Z == 0) 
        return false; // the primary alignment is missing from my SA:Z

    uint32_t my_SA_aln_i=1; // index into my_alns (starting from [1], becaues [0] if the main-fields DEPN alignment)
    for (uint64_t grp_aln_i=1; grp_aln_i < g->num_alns; grp_aln_i++) { // skip 0, iterate on depn alignments of this group
        if (my_SA_aln_i == *prim_aln_index_in_SA_Z) 
            my_SA_aln_i++; // skip primary alignment in SA:Z - we already handled it above

        // case: group alignment i matches one of the alignments in my SA:Z
        // note: alignments in SA:Z are expected to be in the same order as in the prim alignment (i.e. the group)
        if (my_SA_aln_i < n_my_alns && !memcmp (&grp_alns[grp_aln_i], &my_alns [my_SA_aln_i], sizeof(SAAln))) {
            my_SA_aln_i++;
            continue;
        }

        // case: group alignment i matches my alignment (i.e. this line)
        // note: we don't allow a depn line to have my_aln_i=0 (i.e. the primary alignment), bc we make an assumption (eg in sam_reconstruct_main_cigar_from_sag) that a DEPN line cannot have the primary alignment (aln_i=0) in its main fields
        if (grp_aln_i && !memcmp (&grp_alns[grp_aln_i], &my_alns[0], sizeof (SAAln))) 
            *my_aln_i = grp_aln_i;

        // case: group alignment doesn't match this line or any of the alignment on SA:Z of this line
        else
            return false; 
    }    

    // return true;
    return my_SA_aln_i == n_my_alns; // all good if we verified all SA alignments (i.e. exactly one match of my alignment)
}

// parse my SA field, and build my alignments based on my main fields data + SA data
static inline bool sam_sa_seg_depn_get_my_SA_alns (VBlockSAMP vb, ZipDataLineSAMP dl,
                                                   WordIndex my_rname_contig, PosType32 my_pos, uint8_t my_mapq, 
                                                   STRp(my_cigar), int64_t my_nm, bool my_revcomp, 
                                                   STRp(SA), // 0,0 if no SA (eg STAR alignment) 
                                                   uint32_t n_my_alns, 
                                                   SAAln *my_alns) // out 
{
    str_split (SA, SA_len, n_my_alns, ';', SA_aln, true);

    my_alns[0] = (SAAln){ .cigar.signature = cigar_sign (vb, dl, STRa(my_cigar)),
                          .mapq            = my_mapq,
                          .nm              = my_nm,
                          .pos             = my_pos,
                          .revcomp         = my_revcomp,
                          .rname           = my_rname_contig };

    for (uint32_t aln_i=1; aln_i < n_my_alns; aln_i++) { 

        // "rname,pos,strand,CIGAR,mapQ,NM"
        str_split (SA_alns[aln_i-1], SA_aln_lens[aln_i-1], NUM_SA_ITEMS, ',', SA_item, true); 
        if (!n_SA_items) return false;

        // check syntax
        int32_t SA_mapq=0, SA_nm=0, SA_pos=0; 
        if (!str_get_int_range32 (STRi(SA_item, SA_POS),  0, MAX_SA_POS,  &SA_pos)  ||  // pos
            !str_get_int_range32 (STRi(SA_item, SA_NM),   0, MAX_SA_NM,   &SA_nm)   ||  // nm
            !str_get_int_range32 (STRi(SA_item, SA_MAPQ), 0, MAX_SA_MAPQ, &SA_mapq) ||  // mapq
            SA_item_lens[SA_STRAND] != 1                                            || 
            (*SA_items[SA_STRAND] != '-' && *SA_items[SA_STRAND] != '+'))               // strand
            return false;

        my_alns[aln_i] = (SAAln){
            .rname           = ctx_get_ol_node_index_by_snip (VB, CTX(SAM_RNAME), STRi(SA_item, SA_RNAME)),
            .cigar.signature = cigar_sign (vb, NULL, STRi(SA_item, SA_CIGAR)),
            .mapq            = SA_mapq,
            .nm              = SA_nm,
            .pos             = SA_pos,
            .revcomp         = *SA_items[SA_STRAND] == '-',
        };

        if (my_alns[aln_i].rname == WORD_INDEX_NONE) return false;            
    }

    return true; // all good
}

// Seg of DEPN by SA:Z tag: find matching SAGroup and seg SAM_SAG. We identify the group by matching (qname,rname,pos,cigar,revcomp,nm) 
// and matching flags and SEQ
// sets vb->sag (NULL if not found).
static void sam_sa_seg_depn_find_sagroup_SAtag (VBlockSAMP vb, ZipDataLineSAMP dl, 
                                                STRp(textual_cigar), rom textual_seq, bool is_bam)
{
    #define DONE ({ sam_cigar_restore_H (htos); return; })

    vb->sag = NULL; // initialize to "not found"
    vb->sa_aln = NULL;

    bool revcomp = dl->FLAG.rev_comp;
    uint32_t seq_len = dl->SEQ.len;

    buf_alloc (vb, &CTX(SAM_SAG)->local,   1, vb->lines.len, SAGroup,  1, CTX_TAG_LOCAL);
    buf_alloc (vb, &CTX(SAM_SAALN)->local, 1, vb->lines.len, uint16_t, 1, CTX_TAG_LOCAL);

    int64_t grp_index_i=-1;
    uint32_t qname_hash;
    const Sag *g = sam_sa_get_first_group_by_qname_hash (vb, STRqname(dl), dl->FLAG.is_last, &grp_index_i, &qname_hash);
    if (!g) return; // no PRIM with this qname

    // temporarily replace H with S if needed
    HtoS htos = (segconf.SA_HtoS==yes) ? sam_cigar_H_to_S (vb, (char*)STRa(textual_cigar), false) : (HtoS){};

    // if this is BAM, and we have an odd number of bases, the final seq value in the BAM file of the "missing" base
    // we can't encode last "base" other than 0 (see also bug 531)
    if (is_bam && (seq_len&1) && (*B8 (vb->txt_data, dl->SEQ.index + seq_len/2) & 0xf)) DONE;

    ctx_set_encountered (VB, CTX(OPTION_SA_Z));
    uint32_t n_my_alns = str_count_char (STRauxZ(SA_Z, is_bam), ';') + 1; // +1 for alignment of my main fields

    SamNMType my_nm = -1; // -1 means no nm
    SAAln my_alns[n_my_alns];

    // get alignments of this DEPN line ([0]=main field [1...]=SA alignments)
    if (segconf.SA_NM_by_CIGAR_X) my_nm = vb->mismatch_bases_by_CIGAR;
    else if (has(NM_i)) sam_seg_get_aux_int (vb, vb->idx_NM_i, &my_nm, is_bam, MIN_NM_i, MAX_NM_i, HARD_FAIL);

    // populate my_alns: this depn line's alignment in my_alns[0], and the alignments in SA:Z following in my_alns[]
    if (!sam_sa_seg_depn_get_my_SA_alns (vb, dl, vb->chrom_node_index, dl->POS, dl->MAPQ, STRa(textual_cigar), my_nm, revcomp, 
                                         STRauxZ(SA_Z, is_bam), n_my_alns, my_alns))
        DONE; // cannot make sense of the SA alignments

    // iterate on all groups with matching QNAME hash
    do {
        uint16_t my_aln_i=0; // relative to group

        // case: get alignment where SA and my alignment perfectly match the group 
        if (dl->FLAG.is_first   == g->is_first                                   && // note: order tests in the condition is optimized to fail fast with minimal effort
            g->seq_len          == vb->hard_clip[0] + seq_len + vb->hard_clip[1] &&
            dl->FLAG.is_last    == g->is_last                                    &&
            dl->FLAG.multi_segs == g->multi_segs                                 &&
            n_my_alns           == g->num_alns                                   && 
            huffman_issame (SAM_QNAME, GRP_QNAME(g), g->qname_len, STRqname(dl)) &&
            sam_seg_depn_find_SA_aln (vb, g, n_my_alns, my_alns, &my_aln_i, &vb->prim_aln_index_in_SA_Z)     &&
            sam_seg_depn_is_subseq_of_prim (vb, (uint8_t*)textual_seq, dl->SEQ.len, (revcomp != g->revcomp), g, is_bam)) { // this will fail if SEQ has non-ACGT or if DEPN sequence invalidly does not match the PRIM sequence (observed in the wild)

            // found - seg into local
            BNXT (SAGroup, CTX(SAM_SAG)->local) = ZGRP_I(g);
            BNXT16 (CTX(SAM_SAALN)->local) = my_aln_i;
            
            vb->sag = g;
            vb->sa_aln = has(SA_Z) ? B(SAAln, z_file->sag_alns, g->first_aln_i + my_aln_i) : NULL;
            vb->depn_far_count++; // for stats

            if (flag.show_depn || flag.show_buddy)
                 iprintf ("%s: %.*s grp=%u aln=%s qname_hash=%08x\n", LN_NAME, STRfQNAME, ZGRP_I(g), sam_show_sag_one_SA_aln (VB, vb->sag, vb->sa_aln).s, qname_hash);                
                        
            break;
        }

    } while ((g = sam_sa_get_next_group_by_qname_hash (vb, &grp_index_i)));

    DONE;
}

// Seg of DEPN by NH:Z tag: find matching SAGroup and seg SAM_SAG. We identify the group by matching (qname,rname,pos,cigar,revcomp,nm) 
// and matching flags and SEQ
// sets vb->sag (NULL if not found).
static void sam_sa_seg_depn_find_sagroup_noSA (VBlockSAMP vb, ZipDataLineSAMP dl, rom textual_seq, bool is_bam)
{
    vb->sag = NULL; // initialize to "not found"

    bool revcomp = dl->FLAG.rev_comp;
    uint32_t seq_len = dl->SEQ.len;

    buf_alloc (vb, &CTX(SAM_SAG)->local, 1, vb->lines.len, SAGroup, 1, CTX_TAG_LOCAL);

    int64_t grp_index_i=-1;
    uint32_t qname_hash;
    const Sag *g = sam_sa_get_first_group_by_qname_hash (vb, STRqname(dl), dl->FLAG.is_last, &grp_index_i, &qname_hash);

    PosType32 cp = -1;
    STR0(cc); cc="";
    int32_t hi=-1; // stays -1 if the line has no HI:i
    if (flag.show_depn) {
        if (has(HI_i)) sam_seg_get_aux_int (vb, vb->idx_HI_i, &hi, is_bam, 1, 0x7fffffff, SOFT_FAIL);
        if (has(CP_i)) sam_seg_get_aux_int (vb, vb->idx_CP_i, &cp, is_bam, 0, MAX_POS_SAM, SOFT_FAIL);
        if (has(CC_Z)) sam_seg_get_aux_Z (vb, vb->idx_CC_Z, pSTRa(cc), is_bam);
    }

    if (!g) {
        if (flag.show_depn) iprintf ("vb=%u FAIL:NO_QNAME_MATCH QNAME=\"%.*s\"(%08x) HI=%d CC=\"%.*s\" CP=%d\n", vb->vblock_i, STRfQNAME, qname_hash, hi, STRf(cc), cp);
        return; // no PRIM with this qname
    }

    // if this is BAM, and we have an odd number of bases, the final seq value in the BAM file of the "missing" base
    // we can't encode last "base" other than 0 (see also bug 531)
    if (is_bam && (seq_len&1) && (*B8 (vb->txt_data, dl->SEQ.index + seq_len/2) & 0xf)) {
        if (flag.show_depn) iprintf ("vb=%u FAIL:ODD_BASE_NON0 QNAME=\"%.*s\"(%08x) HI=%d CC=\"%.*s\" CP=%d\n", vb->vblock_i, STRfQNAME, qname_hash, hi, STRf(cc), cp);
        return;
    }
    
    int32_t nh;
    if ((IS_SAG_NH || IS_SAG_SOLO || IS_SAG_CC) && !sam_seg_get_aux_int (vb, vb->idx_NH_i, &nh, is_bam, 1, 0x7fffffff, SOFT_FAIL)) {
        if (flag.show_depn) iprintf ("vb=%u FAIL:NO_VALID_NH QNAME=\"%.*s\"(%08x) HI=%d CC=\"%.*s\" CP=%d\n", vb->vblock_i, STRfQNAME, qname_hash, hi, STRf(cc), cp);
        return; // missing or invalid NH:i in depn line
    }

    // iterate on all groups with matching QNAME hash
    do {
        if ((IS_SAG_FLAG ||  nh == g->num_alns)                                  && // note: order tests in the condition is optimized to fail fast with minimal effort
            dl->FLAG.is_first   == g->is_first                                   && 
            g->seq_len          == vb->hard_clip[0] + seq_len + vb->hard_clip[1] &&
            dl->FLAG.is_last    == g->is_last                                    &&
            dl->FLAG.multi_segs == g->multi_segs                                 &&
            huffman_issame (SAM_QNAME, GRP_QNAME(g), g->qname_len, STRqname(dl)) &&
            sam_seg_depn_is_subseq_of_prim (vb, (uint8_t*)textual_seq, dl->SEQ.len, (revcomp != g->revcomp), g, is_bam)) {

            // found - seg into local
            BNXT (SAGroup, CTX(SAM_SAG)->local) = ZGRP_I(g);
            vb->sag = g;
            vb->depn_far_count++; // for stats
            
            if (IS_SAG_CC)
                vb->cc_aln = B(CCAln, z_file->sag_alns, ZGRP_I(g));

            else if (IS_SAG_SOLO)
                vb->solo_aln = B(SoloAln, z_file->sag_alns, ZGRP_I(g));

            if (flag.show_depn || flag.show_buddy) {
                iprintf ("%s: %.*s grp=%u qname_hash=%08x", LN_NAME, STRfQNAME, ZGRP_I(g), qname_hash);
                if (has(HI_i)) iprintf (" HI=%d\n", hi);
                else if (has(CC_Z) && has(CP_i)) iprintf (" CC=\"%.*s\" CP=%d\n", STRf(cc), cp);
                else iprint_newline();
            }
            return; // done
        }

    } while ((g = sam_sa_get_next_group_by_qname_hash (vb, &grp_index_i)));

    if (flag.show_depn) iprintf ("vb=%u FAIL:ALN_MISMATCH QNAME=\"%.*s\"(%08x) HI=%d CC=\"%.*s\" CP=%d\n", 
                                 vb->vblock_i, STRfQNAME, qname_hash, hi, STRf(cc), cp);
}

// Seg compute VB: called when segging PRIM/DEPN VBs
void sam_seg_sag_stuff (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(textual_cigar), rom textual_seq, bool is_bam)
{
    START_TIMER;

    // in Dependent component - try to find which SAGroup this line belongs to, and seg it to SAM_SAG and SAM_SAALN
    // (if successful, this will set vb->sag/sa_aln)
    if (IS_DEPN(vb) && IS_SAG_SA) 
        sam_sa_seg_depn_find_sagroup_SAtag (vb, dl, STRa(textual_cigar), textual_seq, is_bam);

    else if (IS_DEPN(vb) && !IS_SAG_SA) 
        sam_sa_seg_depn_find_sagroup_noSA (vb, dl, textual_seq, is_bam);

    else if (IS_PRIM(vb) && IS_SAG_SA && has(SA_Z)) 
        ctx_set_encountered (VB, CTX(OPTION_SA_Z)); // note: for DEPN, this is done in sam_sa_seg_depn_find_sagroup_SAtag

    COPY_TIMER (sam_seg_sag_stuff);
}

void sam_seg_against_sa_group (VBlockSAMP vb, ContextP ctx, uint32_t add_bytes)
{
    seg_special0 (VB, SAM_SPECIAL_pull_from_sag, ctx, add_bytes);
}

// seg a DEPN line against the SA data in z_file
void sam_seg_against_sa_group_int (VBlockSAMP vb, ContextP ctx, int64_t parameter, uint32_t add_bytes)
{
    SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_pull_from_sag, parameter);
    seg_by_ctx (VB, STRa(snip), ctx, add_bytes);
}

//-------------------
// Zip - VB_HEADER
//-------------------

// Main thread, PRIM VB. Set sam_prim fields of VB_HEADER. Callback from zfile_compress_vb_header
void sam_zip_set_vb_header_specific (VBlockP vb_, SectionHeaderVbHeaderP vb_header)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    if (IS_PRIM(vb)) {
        uint32_t total_seq_len=0;
        for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) 
            total_seq_len += DATA_LINE (line_i)->SEQ.len;   // length in bases, not bytes

        vb_header->sam_prim_seq_len         = BGEN32 (total_seq_len);
        vb_header->sam_prim_comp_qual_len   = BGEN32 (vb->comp_qual_len); 
        vb_header->sam_prim_comp_qname_len  = BGEN32 (CTX(SAM_QNAME)->huffman.comp_len);
        vb_header->sam_prim_num_sag_alns    = BGEN32 ((uint32_t)VB_SAM->sag_alns.count); 
        vb_header->sam_prim_first_grp_i     = BGEN32 (vb->first_grp_i);
        if (IS_SAG_SA)
            vb_header->sam_prim_comp_cigars_len = BGEN32 (vb->comp_cigars_len); 
        else if (IS_SAG_SOLO)
            vb_header->sam_prim_solo_data_len   = BGEN32 (vb->comp_solo_data_len); 
    }

    // In MAIN VBs with gencomp lines, we compress their vb_plan (=dt_specific_vb_header_payload) as the payload of the VB_HEADER section
    if (IS_MAIN(vb) && vb->gencomp_lines.len) {
        BGEN_u16_buf (&vb->vb_plan, NULL);
        
        vb->vb_plan.len *= sizeof (VbPlanItem);
        vb_header->data_uncompressed_len    = BGEN32 (vb->vb_plan.len32);
    }
}

static ASCENDING_SORTER (sort_main_vb_info, SamMainVbInfo, vb_i)
static BINARY_SEARCHER (find_gc_vb_by_line, SamGcVbInfo, uint64_t, first_line, true, IfNotExact_ReturnLower)

// simulate a reconstruction and count max concurrent loaded VBs
uint32_t sam_zip_calculate_max_conc_writing_vbs (void)
{
    START_TIMER;

    // sort MAIN VBs as they were originally added out-of-order
    qsort (STRb(z_file->vb_info[COMP_MAIN]), sizeof(SamMainVbInfo), sort_main_vb_info);

    // byte-map of set when a VB is accessed
    ASSERTNOTINUSE (evb->scratch);
    ARRAY_alloc (bool, simulate_vb_is_loaded, z_file->num_vbs+1, true, evb->scratch, evb, "scratch");
    int32_t conc_vbs = 1; // 1 for the currently reconstrucing MAIN vb
    int32_t conc_vbs_high_watermark = 1; // maximum reached
    
    for_buf (SamMainVbInfo, info, z_file->vb_info[COMP_MAIN]) {
        simulate_vb_is_loaded[info->vb_i] = true;

        for (CompIType comp_i=SAM_COMP_PRIM; comp_i <= SAM_COMP_DEPN; comp_i++) {
            SamGcVbInfo *first_gc_vb_info = binary_search (find_gc_vb_by_line, SamGcVbInfo, z_file->vb_info[comp_i], info->first_gc_line[comp_i-1]);
            SamGcVbInfo *last_gc_vb_info  = binary_search (find_gc_vb_by_line, SamGcVbInfo, z_file->vb_info[comp_i], info->first_gc_line[comp_i-1] + info->num_gc_lines[comp_i-1] - 1);

            if (!first_gc_vb_info) continue; // this file has no lines of this component
            ASSERTNOTNULL (last_gc_vb_info); // if we have a first line, we surely have a last line...
            
            uint32_t remaining_gc_lines = info->num_gc_lines[comp_i-1]; // Number of gencomp lines that originate from this MAIN VB that are of comp_i (PRIM or DEPN)

            for (SamGcVbInfo *info = first_gc_vb_info; info <= last_gc_vb_info; info++) {
                if (!simulate_vb_is_loaded[info->vb_i]) {
                    simulate_vb_is_loaded[info->vb_i] = true;
                    conc_vbs++;
                    if (conc_vbs > conc_vbs_high_watermark) conc_vbs_high_watermark = conc_vbs;
                }

                uint32_t lines_consumed = MIN_(remaining_gc_lines, info->num_lines - info->lines_consumed);
                info->lines_consumed += lines_consumed;

                if (info->num_lines == info->lines_consumed) 
                    conc_vbs--;
            }
        }
    }

    // sanity
    for (VBIType vb_i=1; vb_i <= z_file->vb_info[COMP_MAIN].len + z_file->vb_info[SAM_COMP_PRIM].len + z_file->vb_info[SAM_COMP_DEPN].len/*exc. deep FASTQ*/; vb_i++)
        ASSERT (simulate_vb_is_loaded[vb_i], "gencomp vb_i=%u was not consumed in conc_writing_vbs simulation", vb_i);

    for (CompIType comp_i=SAM_COMP_PRIM; comp_i <= SAM_COMP_DEPN; comp_i++) 
        for_buf (SamGcVbInfo, gc, z_file->vb_info[comp_i])
            ASSERT (gc->lines_consumed == gc->num_lines, "Expecting lines_consumed=%u == num_lines=%u for vb=%s/%u",
                    gc->lines_consumed, gc->num_lines, comp_name (comp_i), gc->vb_i);

    buf_free (evb->scratch);

    if (flag_show_memory || flag_debug_gencomp) 
        iprintf ("\nconc_writing_vbs=%d (consumes memory in piz, not threads)\n\n", conc_vbs_high_watermark);

    // catch values that are non-sensical. work-around: use --no-gencomp
    ASSERTINRANGX (conc_vbs_high_watermark, 0, MAX_CONC_WRITING_VBS);
    
    COPY_TIMER_EVB (sam_zip_calculate_max_conc_writing_vbs);
    return conc_vbs_high_watermark; // note: this may change slightly between executions, because MAIN VBs are absorbed at non-deterministic order, and as a result the order of PRIM and DEPN lines is different between executions
}

// prepare and output the SEC_GENCOMP section
void sam_zip_compress_sec_gencomp (void)
{
    ASSERTNOTINUSE (evb->scratch);

    buf_alloc (evb, &evb->scratch, 0, z_file->vb_info[SAM_COMP_MAIN].len, GencompSecItem, 0, "scratch");

    // note: vb_info[SAM_COMP_MAIN] is in the order of main VB absorption, there the prim/depn lines
    // refered in it are in order. This allows us to store the number of lines only, without their start.
    for_buf (SamMainVbInfo, info, z_file->vb_info[SAM_COMP_MAIN])
        BNXT (GencompSecItem, evb->scratch) = (GencompSecItem){ 
            .vb_i            = BGEN32 (info->vb_i),
            .num_gc_lines[0] = BGEN32 (info->num_gc_lines[0]), // PRIM
            .num_gc_lines[1] = BGEN32 (info->num_gc_lines[1])  // DEPN
        };

    evb->scratch.len *= sizeof (GencompSecItem);
    zfile_compress_section_data (evb, SEC_GENCOMP, &evb->scratch);

    buf_free (evb->scratch);
}

// Callback from stats_get_compressed_sizes: in PRIM VBs, we seg main field CIGAR into OPTION_SA_CIGAR.
// We correct the z_data counts to be shown, so that the right amount of OPTION_SA_CIGAR data is shown as SAM_CIGAR data.
static void sam_stats_reaccount_one (Did main_ctx_did_i, Did sa_ctx_did_i)
{    
    ContextP main_ctx = ZCTX(main_ctx_did_i);
    ContextP sa_ctx   = ZCTX(sa_ctx_did_i);

    // sa_ctx->counts.count - non-strand: txt_len of the PRIM data segged into sa_ctx (NOT added to sa_ctx->txt_len)
    //                        strand: 1 byte per Flag.revcomp segged into SA_STRAND, so comparable to "+" or "-" normally segged into SA_STRAND
    // sa_ctx->txt_len      - the SA data added to sa_ctx
    double pc_main = (sa_ctx->txt_len || sa_ctx->counts.count) 
        ? sa_ctx->counts.count / (sa_ctx->txt_len + sa_ctx->counts.count) // allocated % of z_data of OPTION_SA_*  that is actually due to the main field
        : 1; // move all (added "0" to NM case)

    // move counts of z_data from SA Cigar to main field CIGAR, based on txt_len proportions
    uint64_t amount_to_move = pc_main * sa_ctx->dict.count;
    main_ctx->dict.count  += amount_to_move; // note: stats counts the z_data to be reported in param
    sa_ctx->dict.count    -= amount_to_move;

    amount_to_move = pc_main * sa_ctx->b250.count;
    main_ctx->b250.count  += amount_to_move;
    sa_ctx->b250.count    -= amount_to_move;

    amount_to_move = pc_main * sa_ctx->local.count;
    main_ctx->local.count += amount_to_move;
    sa_ctx->local.count   -= amount_to_move;
}

void sam_stats_reallocate (void)
{    
    sam_stats_reaccount_one (SAM_CIGAR, OPTION_SA_CIGAR); // counts.param accumulated sam_cigar_seg_prim_cigar
    sam_stats_reaccount_one (SAM_RNAME, OPTION_SA_RNAME);
    sam_stats_reaccount_one (SAM_MAPQ,  OPTION_SA_MAPQ);
    sam_stats_reaccount_one (SAM_POS,   OPTION_SA_POS);    
    sam_stats_reaccount_one (SAM_FLAG,  OPTION_SA_STRAND);    

    // Note: if all PRIM lines either have SA+NM or don't, then the math will be correct. In the unlikely case of mixed SA:Z and NH:i PRIM lines, the math will be off.
    if (ZCTX(OPTION_NM_i)->txt_len) // we have NM:i in this file
        sam_stats_reaccount_one (OPTION_NM_i, OPTION_SA_NM);
    
    // case: no NM:i in this field, and hence also no SA:Z - all OPTION_SA_NM data is due to "0" added in sam_seg_sag_stuff, as we will account for it as 
    // belonging to SAM_SAG
    else 
        sam_stats_reaccount_one (SAM_SAG, OPTION_SA_NM); 
}

rom sag_type_name (SagType sagt)
{
    return IN_RANGE (sagt, 0, NUM_SAG_TYPES) ? (rom[])SAM_SAG_TYPE_NAMES[sagt] : "InvalidSagType";
}

