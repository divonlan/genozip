// ------------------------------------------------------------------
//   sam_sag_zip.c
//   Copyright (C) 2022-2022 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "chrom.h"
#include "context.h"
#include "profiler.h"
#include "codec.h"
#include "bits.h"
#include "writer.h"
#include "sections.h"
#include "compressor.h"
#include "bits.h"
#include "flags.h"
#include "recon_plan_io.h"
#include "libdeflate/libdeflate.h"

typedef struct {
    VBIType vb_i;
    uint32_t num_lines;
    uint32_t txt_len; // txt_data.len in this VB - limited by VB size which is at most 2GB
} SamGcVbInfo;

#define SAM_GC_UPDATE_PRIM 0xfffffffe
#define SAM_GC_UPDATE_DEPN 0xffffffff

const SoloProp solo_props[NUM_SOLO_TAGS] = SOLO_PROPS;

//------------------------
// Zip - After Reading VB
//------------------------

// called by main thread after reading a VB - callback of zip_init_vb
void sam_sag_zip_init_vb (VBlockP vb)
{
    if (vb->comp_i == SAM_COMP_MAIN || !vb->txt_data.len) return;

    // PRIM or DEPN - add to vb_info
    buf_alloc (evb, &z_file->vb_info[vb->comp_i-1], 1, 20, SamGcVbInfo, 2, "z_file->vb_info");
    BNXT (SamGcVbInfo, z_file->vb_info[vb->comp_i-1]) = (SamGcVbInfo){
        .vb_i      = vb->vblock_i,
        .num_lines = vb->lines.len,
        .txt_len   = vb->txt_data.len
    };
}

//-------------------
// Zip - Seg
//-------------------

// this is SAM/BAM's zip_after_segconf callback
void sam_set_sag_type (void)
{
    if (flag.no_gencomp || flag.fast || (!segconf.is_sorted && !flag.force_gencomp))
        segconf.sag_type = SAG_NONE;

    else if (MP(LONGRANGER))
        segconf.sag_type = SAG_NONE; // TO DO - new SAG_BY_LONGRANGER type

    // SAG_BY_SA: the preferd Sag method - by SA:Z
    // we identify a PRIM line - if it (1) has SA:Z and (2) supp/secondary flags are clear
    else if (sam_has_SA_Z() || (segconf.has[OPTION_SA_Z])) // longranger has non-standard SA:Z format
        segconf.sag_type = SAG_BY_SA;

    // SAG_BY_NH: STAR and other aligners with NH:i/HI:i instead of SA:Z
    // we identify a PRIM as a line that is (1) not supp/secondary (2) has NH >= 2
    else if ((MP(STAR) || (segconf.has[OPTION_NH_i] && segconf.has[OPTION_HI_i])) && !segconf.has[OPTION_SA_Z])
        segconf.sag_type = segconf.has_barcodes ? SAG_BY_SOLO : SAG_BY_NH; 

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
        segconf.sag_type = SAG_BY_FLAG;
    }

    else 
        segconf.sag_type = SAG_NONE;      

    if (segconf.sag_type) {
        sam_sa_prim_initialize_ingest(); // the PRIM component is compressed (out-of-band) at the same time as MAIN
        gencomp_initialize (SAM_COMP_PRIM, GCT_OOB); 
        gencomp_initialize (SAM_COMP_DEPN, GCT_DEPN); 
    }
}

void sam_seg_gc_initialize (VBlockSAMP vb)
{
    // DEPN stuff
    if (sam_is_depn_vb) {
        CTX(SAM_SAG)->ltype = sizeof(SAGroup)==4 ? LT_UINT32 : LT_UINT64;
        CTX(SAM_SAALN)->ltype   = LT_UINT16; // index of alignment with SA Group
    }

    // PRIM stuff 
    else if (sam_is_prim_vb) { 
        CTX(SAM_SAALN)->ltype   = LT_UINT16;   // index of alignment with SA Group
        CTX(OPTION_SA_Z)->ltype = LT_UINT8;    // we store num_alns in local
        CTX(SAM_FLAG)->no_stons = true;        
        
        ctx_set_store (VB, STORE_INDEX, OPTION_SA_RNAME, OPTION_SA_STRAND, DID_EOL);
        ctx_set_store (VB, STORE_INT, OPTION_SA_POS, OPTION_SA_MAPQ, OPTION_SA_NM, OPTION_NH_i, DID_EOL);
    }
}

// Note: in MAIN vbs, failure means we will keep this line in MAIN and not move it to PRIM
// In PRIM, we abort because this is not not expected as MAIN component should not have moved this line to PRIM
#define FAILIF(condition, format, ...) \
  ( { if (condition) { \
          if (flag.debug_sag) iprintf ("%s: " format "\n", LN_NAME, __VA_ARGS__); \
          if (sam_is_main_vb) return false; \
          else { progress_newline(); fprintf (stderr, "%s: Error in %s:%u: Failed PRIM line because ", LN_NAME, __FUNCLINE); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, SUPPORT); fflush (stderr); exit_on_error(true); }} \
    } )

// Call in seg PRIM line (both in MAIN and PRIM vb): 
// MAIN: test that line is a valid prim before moving it to the PRIM gencomp file
// PRIM: add SA Group (not alignments) to the data structure in VB 
bool sam_seg_prim_add_sag (VBlockSAMP vb, ZipDataLineSAM *dl, uint16_t num_alns/* inc primary aln*/, bool is_bam)
{
    rom textual_seq = is_bam ? B1STc(vb->textual_seq) : Btxt (dl->SEQ.index);
    uint32_t seq_len = is_bam ? vb->textual_seq.len32 : dl->SEQ.len;

    ASSERTNOTNULL (textual_seq);

    // MAIN: check that values are within limits defined in Sag (no need to check in PRIM as we already checked in MAIN)
    if (sam_is_main_vb) {
        FAILIF (!sam_might_have_saggies_in_other_VBs (vb, dl, num_alns), "all sag alignments are contained in this VB%s", "");
        FAILIF (dl->QNAME.len > MAX_SA_QNAME_LEN, "dl->QNAME.len=%u > %u", dl->QNAME.len, MAX_SA_QNAME_LEN);
        FAILIF (seq_len==1 && *textual_seq == '*', "SEQ=\"*\"%s", ""); // we haven't segged seq yet, so vb->seq_missing is not yet set
        FAILIF (seq_len > MAX_SA_SEQ_LEN, "seq_len=%u > %u", seq_len, MAX_SA_SEQ_LEN);
        FAILIF (!dl->POS, "POS=0%s", ""); // unaligned
        FAILIF (num_alns > MAX_SA_NUM_ALNS, "%s=%u > %u", (IS_SAG_NH || IS_SAG_SOLO || IS_SAG_CC) ? "NH:i" : "num_alns", num_alns, MAX_SA_NUM_ALNS);
        FAILIF (dl->hard_clip[0] || dl->hard_clip[1], "has_hard_clips%s", ""); // primary cannot have hard clips
        
        uint32_t bad_i;
        FAILIF (!str_is_only_ACGT (textual_seq, seq_len, &bad_i), 
                "RNAME=\"%.*s\" POS=%d. Found '%c' in base_i=%u in SEQ is not A,C,G or T. SEQ=\"%.*s\"", 
                vb->chrom_name_len, vb->chrom_name, dl->POS, textual_seq[bad_i], bad_i, seq_len, textual_seq);
    }

    // PRIM: actually add the group
    else if (sam_is_prim_vb) { 
        buf_alloc (vb, &vb->sag_grps, 1, 1000, Sag, CTX_GROWTH, "sag_grps");

        // add SA group
        BNXT (Sag, vb->sag_grps) = (Sag){
            .first_aln_i = vb->sag_alns.len,
            .num_alns    = num_alns,
            .qname       = dl->QNAME.index,
            .qname_len   = dl->QNAME.len,
            .qual        = dl->QUAL.index,
            .no_qual     = vb->qual_missing,
            .revcomp     = dl->FLAG.rev_comp,
            .multi_segs  = dl->FLAG.multi_segs,
            .is_first    = dl->FLAG.is_first,
            .is_last     = dl->FLAG.is_last,
            .seq         = dl->SEQ.index,
            .seq_len     = seq_len,
            .as          = CAP_SA_AS (dl->AS) // [0,255] capped at 0 and 255
        };
    }

    return true;
}

// MAIN: called in sam_seg_is_gc_line to test that line is a valid prim before moving it to the PRIM gencomp file
// PRIM: call in seg of seg of SA:Z add SA Group (including alignments) to the data structure in VB
// return n_alns if successful and 0 if not
int32_t sam_seg_prim_add_sag_SA (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(sa), int64_t this_nm, bool is_bam)
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

    if (sam_is_prim_vb) { // in PRIM, we actually add the group, in MAIN, we are just testing
        buf_alloc (vb, &vb->sag_alns, n_alns /*+1 for primary aln*/, 64, SAAln, CTX_GROWTH, "sag_alns");

        // in Seg, we add the Adler32 of CIGAR to save memory, as we just need to verify it, not reconstruct it
        uint32_t cigar_len = (is_bam ? vb->textual_cigar.len32 : dl->CIGAR.len);
        CigarSignature cigar_sig = cigar_sign (is_bam ? B1STc(vb->textual_cigar) : Btxt(dl->CIGAR.index), cigar_len);
        
        // add primary alignment as first alignment in SA group
        BNXT (SAAln, vb->sag_alns) = (SAAln){
            .rname           = vb->chrom_node_index,
            .pos             = dl->POS,
            .revcomp         = dl->FLAG.rev_comp,
            .cigar.signature = cigar_sig,
            .mapq            = dl->MAPQ,
            .nm              = this_nm
        };

        // in BAM, add textual CIGAR to vb->sag_cigars, as txt_data CIGAR is binary
        if (is_bam)
            buf_add_buf (vb, &vb->sa_prim_cigars, &vb->textual_cigar, char, "sa_prim_cigars");

        vb->sag_alns.count += n_alns;
    }

    // test or add all dependent alignments
    for (uint32_t i=0; i < n_alns-1; i++) {
        
        // get items - rname, pos, revcomp, CIGAR, mapQ, NM
        str_split (alns[i], aln_lens[i], NUM_SA_ITEMS, ',', item, true);
        FAILIF (n_items != NUM_SA_ITEMS, "in SA alignment %u - n_items=%u != NUM_SA_ITEMS=%u", i, n_items, NUM_SA_ITEMS);

        // rname (pre-populated from sam header)
        WordIndex rname = ctx_get_ol_node_index_by_snip (VB, CTX(SAM_RNAME), STRi(item,SA_RNAME));
        FAILIF (rname == WORD_INDEX_NONE, "in SA alignment %u - rname==WORD_INDEX_NONE", i);

        // protect rname from removal by ctx_shorten_unused_dict_words if it is unused in the main RNAME field. 
        ctx_protect_from_removal (VB, CTX(SAM_RNAME), rname); 

        // pos, mapq, nm - get integers and verify limits
        int64_t pos, mapq, nm;
        FAILIF (!str_get_int_range64 (STRi(item,SA_POS ), 0, MAX_SA_POS,  &pos ), "in SA alignment %u - pos=%"PRId64" out of range [0, %u]", i, pos, MAX_SA_POS);
        FAILIF (!str_get_int_range64 (STRi(item,SA_MAPQ), 0, MAX_SA_MAPQ, &mapq), "in SA alignment %u - mapq=%"PRId64" out of range [0, %u]", i, mapq, MAX_SA_MAPQ);
        FAILIF (!str_get_int_range64 (STRi(item,SA_NM  ), 0, MAX_SA_NM,   &nm  ), "in SA alignment %u - nm=%"PRId64" out of range [0, %u]", i, nm, MAX_SA_NM);
 
        // revcomp
        bool revcomp = *items[SA_STRAND] == '-';
        FAILIF (item_lens[SA_STRAND] != 1 || (!revcomp && *items[SA_STRAND] != '+'), "in SA alignment %u - STRAND is not '+' or '-'", i);

        if (sam_is_prim_vb) { // in PRIM, we actually add the group, in MAIN, we are just testing

            if (!segconf.SA_CIGAR_can_have_H && sam_cigar_has_H (STRi(item,SA_CIGAR)))
                segconf.SA_CIGAR_can_have_H = true; // no worries about thread safety, this just gets set to true while segging MAIN and is inspected only when segging MAIN is over

            CigarSignature cigar_sig = cigar_sign (STRi(item,SA_CIGAR));

            BNXT (SAAln, vb->sag_alns) = (SAAln){ 
                .rname           = rname, 
                .pos             = pos, 
                .revcomp         = revcomp, 
                .cigar.signature = cigar_sig,
                .mapq            = mapq,
                .nm              = nm
            };
        }
    }

    return n_alns;
}       

// PRIM: call in seg of seg of NH:i add SA Group to the data structure in VB 
void sam_seg_prim_add_sag_NH (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t nh)
{
    sam_seg_prim_add_sag (vb, dl, nh, IS_BAM_ZIP);
    vb->sag_alns.count += nh;
}

// PRIM: call in seg of seg of NH:i add SA Group to the data structure in VB 
void sam_seg_prim_add_sag_CC (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t nh)
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

// PRIM: called to ingested a segged VB data into solo "alignments" 
void sam_seg_prim_add_sag_SOLO (VBlockSAMP vb, ZipDataLineSAM *dl)
{
    sam_seg_prim_add_sag (vb, dl, dl->NH, IS_BAM_ZIP);

    // TO DO: consider the tradeoff - this way, we don't allocate memory and add in each line (like sam_seg_prim_add_sag_CC)
    // BUT, we do hundreds of thousands of small memcpy's during ingest - with VBs serialized on a mutex, instead of just
    // one big memcpy
}

// Seg MAIN: returns true if this line is a Primary or Dependent line of a supplementary/secondary group - 
// and should be moved to a generated component. Called by compute thread in seg of Normal VB.
bool sam_seg_is_gc_line (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(alignment), bool is_bam)
{
    START_TIMER;

    SamNMType NM;
    int32_t n_alns; 
    SamComponentType comp_i = COMP_NONE; // generated component
    
    // unmapped or otherwise lacking alignment - not gencomp
    if (!segconf.sag_type || vb->chrom_node_index == WORD_INDEX_NONE || // RNAME='*' or headerless SAM (we need header-contigs for sam_sa_add_sa_group)
        dl->FLAG.unmapped || vb->cigar_missing || !dl->POS)
            goto done; 

    switch (segconf.sag_type) {
        case SAG_BY_SA:
            if ((vb->hard_clip[0]>0 || vb->hard_clip[1]>0) && (vb->soft_clip[0]>0 || vb->soft_clip[1]>0))
                goto done; // we don't support adding an alignment with both soft and hard clips to an SA-based sag

            if (has_SA && sam_line_is_depn(dl)) {
                comp_i = SAM_COMP_DEPN;

                if (vb->qname_count.len32)
                    n_alns = str_count_char (STRauxZ (SA_Z, is_bam), ';') + 1; // +1 for this aln
            
                if (dl->hard_clip[0] || dl->hard_clip[1])
                   segconf.depn_CIGAR_can_have_H = true; // no worries about thread safety, this just gets set to true while segging MAIN and is inspected only when segging MAIN is over
            }
            
            else if (has_SA && has_NM &&
                     sam_seg_get_aux_int (vb, vb->idx_NM_i, &NM, is_bam, MIN_NM_i, MAX_NM_i, true) &&
                     (n_alns = sam_seg_prim_add_sag_SA (vb, dl, STRauxZ (SA_Z, is_bam), NM, is_bam)))  // testing to see if we can successfully add a sag based on SA
                comp_i = SAM_COMP_PRIM;
            break;
        
        case SAG_BY_NH: 
        case SAG_BY_SOLO: 
            if (has_NH && sam_seg_get_aux_int (vb, vb->idx_NH_i, &n_alns, is_bam, 2/*at least*/, MAX_HI_NH, true)) {
                
                if (sam_line_is_depn(dl)) 
                    comp_i = SAM_COMP_DEPN;

                else if (sam_seg_prim_add_sag (vb, dl, n_alns, is_bam)) // testing to see if we can successfully add a sag based on NH
                    comp_i = SAM_COMP_PRIM;
            }
            break;
            
        case SAG_BY_CC:
            if (has_NH && sam_seg_get_aux_int (vb, vb->idx_NH_i, &n_alns, is_bam, 2, MAX_HI_NH, true)) {  // not out of range, i.e. at least 2

                if (has_CC && has_CP)
                    comp_i = SAM_COMP_DEPN;

                else if (!has_CC && !has_CP && !has_SA && !has_HI && !dl->hard_clip[0] && !dl->hard_clip[1] && 
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
        comp_i = COMP_NONE;

    if (comp_i != COMP_NONE) {

        // store location where this gc line should be inserted    
        gencomp_seg_add_line (VB, comp_i, STRa(alignment));
        
        vb->line_i--;
        vb->recon_size -= alignment_len;
        vb->txt_size   -= alignment_len;

        if      (comp_i == SAM_COMP_PRIM) vb->seg_found_prim_line = true;
        else if (comp_i == SAM_COMP_DEPN) vb->seg_found_depn_line = true;
    }

done:    
    COPY_TIMER (sam_seg_is_gc_line);
    
    return comp_i != COMP_NONE; // true if line moves to generated component
}

// ------------------
// QNAME
// ------------------

typedef struct __attribute__ ((__packed__)) { uint32_t qname_hash, grp_i; } SAGroupIndexEntry; 

// ZIP DEPN: search for index entry by qname_hash 
static int64_t sam_sa_binary_search_for_qname_hash (const SAGroupIndexEntry *index, uint64_t this_qname_hash, int64_t first, int64_t last)
{
    if (first > last) return -1; // not found

    int64_t mid = (first + last) / 2;

    int64_t cmp = (int64_t)index[mid].qname_hash - (int64_t)this_qname_hash;
    if      (cmp < 0) return sam_sa_binary_search_for_qname_hash (index, this_qname_hash, mid+1, last);
    else if (cmp > 0) return sam_sa_binary_search_for_qname_hash (index, this_qname_hash, first, mid-1);

    // mid contains the this_qname_hash - scan the index backwards for the first matching entry as there might be several different qnames with the same hash
    while (mid >= 1 && index[mid-1].qname_hash == this_qname_hash) mid--;

    return mid;
}

// ZIP DEPN: find group index with this_qname in z_file->sag_grps, and if there are several - return the first
static Sag *sam_sa_get_first_group_by_qname_hash (VBlockSAMP vb, STRp(this_qname), bool is_last, int64_t *grp_index_i, 
                                                          uint32_t *this_qname_hash) // out
{
    // search for a group with qname in z_file->sa_qname
    *this_qname_hash = QNAME_HASH (this_qname, this_qname_len, is_last);
    *grp_index_i = sam_sa_binary_search_for_qname_hash (B1ST (SAGroupIndexEntry, z_file->sag_gps_index), *this_qname_hash, 0, z_file->sag_gps_index.len-1);

    const SAGroupIndexEntry *index_ent = B(SAGroupIndexEntry, z_file->sag_gps_index, *grp_index_i); // invalid pointer if grp_index_i==-1, that's ok    
    return (*grp_index_i >= 0) ? B(Sag, z_file->sag_grps, index_ent->grp_i) : NULL; 
}

// ZIP DEPN: if there are more groups in z_file->sag_grps with the qname adlers, return the next group index
static Sag *sam_sa_get_next_group_by_qname_hash (VBlockSAMP vb, int64_t *grp_index_i)
{
    const SAGroupIndexEntry *index_ent = B(SAGroupIndexEntry, z_file->sag_gps_index, *grp_index_i);

    if (*grp_index_i < z_file->sag_gps_index.len-1 && index_ent->qname_hash == (index_ent+1)->qname_hash) {
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
                                             uint16_t *my_aln_i) // out - relative to group
{
    uint32_t SA_aln_i=1; // index into my_alns
    for (uint64_t grp_aln_i=0; grp_aln_i < g->num_alns; grp_aln_i++) { 

        const SAAln *grp_aln = B(SAAln, z_file->sag_alns, g->first_aln_i + grp_aln_i);
        const SAAln *SA_aln = &my_alns[SA_aln_i];

        // case: group alignment matches the SA alignment
        if (SA_aln_i < n_my_alns && !memcmp (grp_aln, SA_aln, sizeof(SAAln))) {
            SA_aln_i++;
            continue;
        }

        // case: group alignment matches my alignment
        if (grp_aln_i && !memcmp (grp_aln, &my_alns[0], sizeof (SAAln))) // note: we don't allow a depn line to have my_aln_i=0, bc we make an assumption (eg in sam_reconstruct_main_cigar_from_sag) that a DEPN line cannot have the primary alignment (aln_i=0) in its main fields
            *my_aln_i = grp_aln_i;

        // case: group alignment doesn't match the SA alignment and doesn't match my alignment
        else
            return false; 
    }    

    return SA_aln_i == n_my_alns; // all good if we verified all SA alignments (i.e. exactly one match of my alignment)
}

// parse my SA field, and build my alignments based on my main fields data + SA data
static inline bool sam_sa_seg_depn_get_my_SA_alns (VBlockSAMP vb,
                                                   WordIndex my_rname_contig, SamPosType my_pos, uint8_t my_mapq, 
                                                   STRp(my_cigar), int64_t my_nm, bool my_revcomp, 
                                                   STRp(SA), // 0,0 if no SA (eg STAR alignment) 
                                                   uint32_t n_my_alns, 
                                                   SAAln *my_alns) // out 
{
    str_split (SA, SA_len, n_my_alns, ';', SA_aln, true);

    my_alns[0] = (SAAln){ .cigar.signature = cigar_sign (STRa(my_cigar)),
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
        int32_t SA_mapq, SA_nm, SA_pos; 
        if (n_SA_items &&
            (   !str_get_int_range32 (STRi(SA_item, SA_POS),  0, MAX_SA_POS,  &SA_pos)  ||  // pos
                !str_get_int_range32 (STRi(SA_item, SA_NM),   0, MAX_SA_NM,   &SA_nm)   ||  // nm
                !str_get_int_range32 (STRi(SA_item, SA_MAPQ), 0, MAX_SA_MAPQ, &SA_mapq) ||  // mapq
                SA_item_lens[SA_STRAND] != 1                                            || 
                (*SA_items[SA_STRAND] != '-' && *SA_items[SA_STRAND] != '+')))              // strand
            return false;

        my_alns[aln_i] = (SAAln){
            .rname           = ctx_get_ol_node_index_by_snip (VB, CTX(SAM_RNAME), STRi(SA_item, SA_RNAME)),
            .cigar.signature = cigar_sign (STRi(SA_item, SA_CIGAR)),
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
static void sam_sa_seg_depn_find_sagroup_SAtag (VBlockSAMP vb, ZipDataLineSAM *dl, 
                                                STRp(textual_cigar), rom textual_seq, bool is_bam)
{
    #define DONE ({ sam_cigar_restore_H (htos); return; })

    vb->sag = NULL; // initialize to "not found"
    vb->sa_aln = NULL;

    bool revcomp = dl->FLAG.rev_comp;
    uint32_t seq_len = dl->SEQ.len;

    buf_alloc (vb, &CTX(SAM_SAG)->local, 1, vb->lines.len, SAGroup, 1, "contexts->local");
    buf_alloc (vb, &CTX(SAM_SAALN)->local,   1, vb->lines.len, uint16_t,  1, "contexts->local");

    int64_t grp_index_i=-1;
    uint32_t qname_hash;
    const Sag *g = sam_sa_get_first_group_by_qname_hash (vb, STRtxtw(dl->QNAME), dl->FLAG.is_last, &grp_index_i, &qname_hash);
    if (!g) return; // no PRIM with this qname

    // temporarily replace H with S if needed
    HtoS htos = (segconf.SA_HtoS==yes) ? sam_cigar_H_to_S (vb, (char*)STRa(textual_cigar), false) : (HtoS){};

    // if this is BAM, and we have an odd number of bases, the final seq value in the BAM file of the "missing" base
    // we can't encode last "base" other than 0 (see also bug 531)
    if (is_bam && (seq_len&1) && (*B8 (vb->txt_data, dl->SEQ.index + seq_len/2) & 0xf)) DONE;

    // we can't encode depn seq bases other than A,C,G,T as we operate in 2bit
    if (!str_is_only_ACGT (textual_seq, seq_len, NULL)) DONE;

    ctx_set_encountered (VB, CTX(OPTION_SA_Z));
    uint32_t n_my_alns = str_count_char (STRauxZ(SA_Z, is_bam), ';') + 1; // +1 for alignment of my main fields

    SamNMType my_nm = -1; // -1 means no nm
    SAAln my_alns[n_my_alns];

    // get alignments of this DEPN line ([0]=main field [1...]=SA alignments)
    if (has_NM) sam_seg_get_aux_int (vb, vb->idx_NM_i, &my_nm, is_bam, MIN_NM_i, MAX_NM_i, false);

    if (!sam_sa_seg_depn_get_my_SA_alns (vb, vb->chrom_node_index, dl->POS, dl->MAPQ, STRa(textual_cigar), my_nm, revcomp, 
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
            str_issame_(STRtxtw(dl->QNAME), GRP_QNAME(g), g->qname_len)          && 
            sam_seg_depn_find_SA_aln (vb, g, n_my_alns, my_alns, &my_aln_i)) {

            // found - seg into local
            BNXT (SAGroup, CTX(SAM_SAG)->local) = ZGRP_I(g);
            BNXT16 (CTX(SAM_SAALN)->local) = my_aln_i;
            
            vb->sag = g;
            vb->sa_aln = has_SA ? B(SAAln, z_file->sag_alns, g->first_aln_i + my_aln_i) : NULL;
            vb->prim_far_count++; // for stats

            if (flag.show_depn || flag.show_buddy)
                 iprintf ("%s: %.*s grp=%u aln=%s qname_hash=%08x\n", LN_NAME, STRfw(dl->QNAME), ZGRP_I(g), sam_show_sag_one_aln (vb->sag, vb->sa_aln).s, qname_hash);                
                        
            break;
        }

    } while ((g = sam_sa_get_next_group_by_qname_hash (vb, &grp_index_i)));

    DONE;
}

// Seg of DEPN by NH:Z tag: find matching SAGroup and seg SAM_SAG. We identify the group by matching (qname,rname,pos,cigar,revcomp,nm) 
// and matching flags and SEQ
// sets vb->sag (NULL if not found).
static void sam_sa_seg_depn_find_sagroup_noSA (VBlockSAMP vb, ZipDataLineSAM *dl, rom textual_seq, bool is_bam)
{
    vb->sag = NULL; // initialize to "not found"

    bool revcomp = dl->FLAG.rev_comp;
    uint32_t seq_len = dl->SEQ.len;

    buf_alloc (vb, &CTX(SAM_SAG)->local, 1, vb->lines.len, SAGroup, 1, "contexts->local");

    int64_t grp_index_i=-1;
    uint32_t qname_hash;
    const Sag *g = sam_sa_get_first_group_by_qname_hash (vb, STRtxtw(dl->QNAME), dl->FLAG.is_last, &grp_index_i, &qname_hash);

    SamPosType cp = -1;
    STR0(cc); cc="";
    int32_t hi=-1; // stays -1 if the line has no HI:i
    if (flag.show_depn) {
        if (has_HI) sam_seg_get_aux_int (vb, vb->idx_HI_i, &hi, is_bam, 1, 0x7fffffff, true);
        if (has_CP) sam_seg_get_aux_int (vb, vb->idx_CP_i, &cp, is_bam, 0, MAX_POS_SAM, true);
        if (has_CC) sam_seg_get_aux_Z (vb, vb->idx_CC_Z, pSTRa(cc), is_bam);
    }

    if (!g) {
        if (flag.show_depn) iprintf ("vb=%u FAIL:NO_QNAME_MATCH QNAME=\"%.*s\"(%08x) HI=%d CC=\"%.*s\" CP=%d\n", vb->vblock_i, STRfw(dl->QNAME), qname_hash, hi, STRf(cc), cp);
        return; // no PRIM with this qname
    }

    // if this is BAM, and we have an odd number of bases, the final seq value in the BAM file of the "missing" base
    // we can't encode last "base" other than 0 (see also bug 531)
    if (is_bam && (seq_len&1) && (*B8 (vb->txt_data, dl->SEQ.index + seq_len/2) & 0xf)) {
        if (flag.show_depn) iprintf ("vb=%u FAIL:ODD_BASE_NON0 QNAME=\"%.*s\"(%08x) HI=%d CC=\"%.*s\" CP=%d\n", vb->vblock_i, STRfw(dl->QNAME), qname_hash, hi, STRf(cc), cp);
        return;
    }
    
    int32_t nh;
    if ((IS_SAG_NH || IS_SAG_SOLO || IS_SAG_CC) && !sam_seg_get_aux_int (vb, vb->idx_NH_i, &nh, is_bam, 1, 0x7fffffff, true)) {
        if (flag.show_depn) iprintf ("vb=%u FAIL:NO_VALID_NH QNAME=\"%.*s\"(%08x) HI=%d CC=\"%.*s\" CP=%d\n", vb->vblock_i, STRfw(dl->QNAME), qname_hash, hi, STRf(cc), cp);
        return; // missing or invalid NH:i in depn line
    }

    // iterate on all groups with matching QNAME hash
    do {
        if ((IS_SAG_FLAG ||  nh == g->num_alns)  && // note: order tests in the condition is optimized to fail fast with minimal effort
            dl->FLAG.is_first   == g->is_first   && 
            g->seq_len          == vb->hard_clip[0] + seq_len + vb->hard_clip[1] &&
            dl->FLAG.is_last    == g->is_last    &&
            dl->FLAG.multi_segs == g->multi_segs &&
            str_issame_(STRtxtw(dl->QNAME), GRP_QNAME(g), g->qname_len) && 
            sam_seg_depn_is_subseq_of_prim (vb, (uint8_t*)textual_seq, dl->SEQ.len, (revcomp != g->revcomp), g, is_bam)) {

            // found - seg into local
            BNXT (SAGroup, CTX(SAM_SAG)->local) = ZGRP_I(g);
            vb->sag = g;
            vb->prim_far_count++; // for stats
            
            if (IS_SAG_CC)
                vb->cc_aln = B(CCAln, z_file->sag_alns, ZGRP_I(g));

            else if (IS_SAG_SOLO)
                vb->solo_aln = B(SoloAln, z_file->sag_alns, ZGRP_I(g));

            if (flag.show_depn || flag.show_buddy) {
                iprintf ("%s: %.*s grp=%u qname_hash=%08x", LN_NAME, STRfw(dl->QNAME), ZGRP_I(g), qname_hash);
                if (has_HI) iprintf (" HI=%d\n", hi);
                else if (has_CC && has_CP) iprintf (" CC=\"%.*s\" CP=%d\n", STRf(cc), cp);
                else iprint0 ("\n");
            }
            return; // done
        }

    } while ((g = sam_sa_get_next_group_by_qname_hash (vb, &grp_index_i)));

    if (flag.show_depn) iprintf ("vb=%u FAIL:ALN_MISMATCH QNAME=\"%.*s\"(%08x) HI=%d CC=\"%.*s\" CP=%d\n", 
                                 vb->vblock_i, STRfw(dl->QNAME), qname_hash, hi, STRf(cc), cp);
}

// Seg compute VB: called when segging PRIM/DEPN VBs
void sam_seg_sag_stuff (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(textual_cigar), rom textual_seq, bool is_bam)
{
    START_TIMER;

    // in Dependent component - try to find which SAGroup this line belongs to, and seg it to SAM_SAG and SAM_SAALN
    // (if successful, this will set vb->sag/sa_aln)
    if (sam_is_depn_vb && IS_SAG_SA) 
        sam_sa_seg_depn_find_sagroup_SAtag (vb, dl, STRa(textual_cigar), textual_seq, is_bam);

    else if (sam_is_depn_vb && !IS_SAG_SA) 
        sam_sa_seg_depn_find_sagroup_noSA (vb, dl, textual_seq, is_bam);

    else if (sam_is_prim_vb && IS_SAG_SA && has_SA) 
        ctx_set_encountered (VB, CTX(OPTION_SA_Z)); // for DEPN, this is done in sam_sa_seg_depn_find_sagroup_SAtag

    COPY_TIMER (sam_seg_sag_stuff);
}

void sam_seg_against_sa_group (VBlockSAMP vb, ContextP ctx, uint32_t add_bytes)
{
    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SAG }, 2, ctx, add_bytes);
}

// seg a DEPN line against the SA data in z_file
void sam_seg_against_sa_group_int (VBlockSAMP vb, ContextP ctx, int64_t parameter, uint32_t add_bytes)
{
    SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_SAG, parameter);
    seg_by_ctx (VB, STRa(snip), ctx, add_bytes);
}

//---------------------
// Zip - After Compute
//---------------------

// called in the main thread after_compute - VBs might be out of order
void sam_zip_gc_after_compute_main (VBlockSAMP vb)
{
    BufferP recon_plan = &txt_file->recon_plan;
    BufferP recon_plan_index = &txt_file->recon_plan_index;
    ARRAY (GencompLineIEntry, gc_lines, vb->gencomp_lines);

    buf_alloc_zero (evb, recon_plan_index, 0, MAX_(1000, vb->vblock_i+1), BufWord, 2, "txt_file->recon_plan_index");
    recon_plan_index->len32 = MAX_(recon_plan_index->len32, vb->vblock_i+1);

    uint64_t recon_plan_vb_start = recon_plan->len;
    
    // we create a "Full VB" recon plan if VB has no Supplementary/Secondary or Dependent lines 
    // Other VBs might have them. If no VB has them, we will get rid of the recon_plan in sam_zip_generate_recon_plan
    if (!gc_lines_len) {
        buf_alloc (evb, recon_plan, 1, 100000, ReconPlanItem, 2, "txt_file->recon_plan");

        BNXT (ReconPlanItem, *recon_plan) = (ReconPlanItem){
            .flavor = PLAN_FULL_VB,
            .vb_i   = vb->vblock_i
        };
    }

    // case: we removed some gencomp lines from this VB, we now create a recon plan to insert them back.
    // the actual line location in the gencomp vb will be updated later
    else {         
        buf_alloc (evb, recon_plan, gc_lines_len * 2 + 2, 100000, ReconPlanItem, 2, "txt_file->recon_plan");
        buf_alloc (evb, &txt_file->line_info[0], gc_lines_len, 100000, uint32_t, 2, "txt_file->line_info");
        buf_alloc (evb, &txt_file->line_info[1], gc_lines_len, 100000, uint32_t, 2, "txt_file->line_info");

        uint32_t normal_line_i=0;
    
        for (uint32_t gc_line_i=0 ; gc_line_i < gc_lines_len; gc_line_i++) {
            GencompLineIEntry gc_line = gc_lines[gc_line_i];
            
            // insert normal lines before the next gc line
            if (gc_line.line_i > normal_line_i) {

                BNXT (ReconPlanItem, *recon_plan) = (ReconPlanItem){
                    .vb_i       = vb->vblock_i,
                    .start_line = normal_line_i,
                    .num_lines  = gc_line.line_i - normal_line_i
                };
 
                normal_line_i = gc_line.line_i; // next normal line will be after this gc line and possibly subsequent ones that are marked as coming before it
            }

            // insert gc lines - a bunch of them that have the same component and are consecutive
            // vb_i within the gencomp components and start_line will be updated later as we don't know them yet
            BNXT (ReconPlanItem, *recon_plan) = 
                (ReconPlanItem){ .vb_i = (gc_line.comp_i == SAM_COMP_PRIM) ? SAM_GC_UPDATE_PRIM : SAM_GC_UPDATE_DEPN };
          
            // store line lengths, to be used later to calculate vb_info
            BNXT32 (txt_file->line_info[gc_line.comp_i-1]) = gc_line.line_len;
        }

        // insert final normal lines
        if (normal_line_i < vb->lines.len32) 
            BNXT (ReconPlanItem, *recon_plan) = (ReconPlanItem){
                .vb_i       = vb->vblock_i,
                .start_line = normal_line_i,
                .num_lines  = vb->lines.len32 - normal_line_i
            };

        // insert end-of-VB (note: we insert after all lines, not just normal lines, since it is at this point that
        // writer_main_loop calculates the digest)
        BNXT (ReconPlanItem, *recon_plan) = (ReconPlanItem){ 
            .vb_i   = vb->vblock_i,
            .flavor = PLAN_END_OF_VB,
        }; 
    }

    *B(BufWord, *recon_plan_index, vb->vblock_i) = (BufWord){ .index = recon_plan_vb_start, 
                                                              .len   = recon_plan->len - recon_plan_vb_start };
}

//-------------------
// Zip - VB_HEADER
//-------------------

// Main thread, PRIM VB. Set sam_prim fields of VB_HEADER. Callback from zfile_compress_vb_header
void sam_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeader *vb_header)
{
    if (sam_is_prim_vb) {
        uint32_t total_seq_len=0, total_qname_len=0;
        for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) {
            ZipDataLineSAM *dl = DATA_LINE (line_i);
            total_seq_len   += dl->SEQ.len;   // length in bases, not bytes
            total_qname_len += dl->QNAME.len;
        }

        vb_header->sam_prim_seq_len          = BGEN32 (total_seq_len);
        vb_header->sam_prim_comp_qual_len    = BGEN32 (VB_SAM->comp_qual_len);
        vb_header->sam_prim_qname_len        = BGEN32 (total_qname_len);
        vb_header->sam_prim_num_sag_alns     = BGEN32 ((uint32_t)VB_SAM->sag_alns.count); 
        vb_header->sam_prim_first_grp_i      = BGEN32 (VB_SAM->first_grp_i);
        vb_header->sam_prim_comp_cigars_len  = BGEN32 (VB_SAM->comp_cigars_len); // note: both the header and VB field are a union with solo_data_len 
    }
}

//------------------
// Zip - RECON_PLAN
//------------------

// callback function of compress_recon_plan, called from zip_one_file
static bool sam_zip_recon_plan_full_vb_only (void)
{
    ARRAY (ReconPlanItem, recon_plan, txt_file->recon_plan);

    for (uint32_t i=0; i < recon_plan_len; i++)
        if (recon_plan[i].flavor != PLAN_FULL_VB) return false;

    return true;
}

// ZIP main thread: merge consecutive recon plan entries of the same VB to reduce the recon plan size
static void sam_zip_recon_plan_optimize_gc_lines (void)
{
    ARRAY (ReconPlanItem, recon_plan, txt_file->recon_plan);

    uint32_t last_i = 0;
    for (uint32_t i=1; i < recon_plan_len; i++) 
        // range items - merge if same vb_i
        if (recon_plan[i]     .flavor == PLAN_RANGE 
         && recon_plan[last_i].flavor == PLAN_RANGE 
         && recon_plan[i].vb_i == recon_plan[last_i].vb_i) 
            recon_plan[last_i].num_lines++;
                
        else 
            recon_plan[++last_i] = recon_plan[i];

    txt_file->recon_plan.len = last_i + 1;
}

// ZIP main thread: compressing recon plan after DEPN: update recon plan entries of PRIM and DEPN with their details
static uint32_t sam_zip_recon_plan_add_gc_lines (void)
{
    START_TIMER;

    #define UPDATE_WRITERS_END_VB(vb_i) ({  \
        if (vb_in_use[vb_i]) {              \
            curr_conc_writers--;            \
            vb_in_use[vb_i] = false;        \
        } })

    #define UPDATE_WRITERS(vb_i) ({         \
        if (!vb_in_use[vb_i]) {             \
            vb_in_use[vb_i] = true;         \
            curr_conc_writers++;            \
            if (curr_conc_writers > max_conc_writers) max_conc_writers = curr_conc_writers; \
        } })

    SamGcVbInfo *vb_info[2] =             { B1ST (SamGcVbInfo, z_file->vb_info[0]),   // PRIM 
                                            B1ST (SamGcVbInfo, z_file->vb_info[1]) }; // DEPN

    const SamGcVbInfo *after_vb_info[2] = { BAFT (SamGcVbInfo, z_file->vb_info[0]),   // PRIM
                                            BAFT (SamGcVbInfo, z_file->vb_info[1]) }; // DEPN

    uint32_t gc_vb_line_i[2] = { 0, 0 }; // PRIM/DEPN line within current VB
    uint32_t curr_conc_writers=0, max_conc_writers=0;

    // byte-map of set when a VB is accessed
    ASSERTNOTINUSE (evb->scratch);
    ARRAY_alloc (bool, vb_in_use, z_file->num_vbs, true, evb->scratch, evb, "scratch");

    for (uint64_t i=0; i < txt_file->recon_plan.len; i++) { // note: recon_plan may get extended within the loop with INSERBtxtAFTER

        ReconPlanItem *pi = B(ReconPlanItem, txt_file->recon_plan, i);
        VBIType vb_i = pi->vb_i;

        // case: not a depn or prim plan item
        if (vb_i < SAM_GC_UPDATE_PRIM) {
            if (pi->flavor == PLAN_END_OF_VB) 
                UPDATE_WRITERS_END_VB (vb_i);

            else if (pi->flavor == PLAN_FULL_VB) {
                UPDATE_WRITERS (vb_i);
                UPDATE_WRITERS_END_VB (vb_i);
            }

            else 
                UPDATE_WRITERS (vb_i);
            
            continue; // not a gc item
        }

        bool is_depn = (vb_i == SAM_GC_UPDATE_DEPN);

        ASSERTNOTNULL (vb_info[is_depn]);
        ASSERT (vb_info[is_depn] < after_vb_info[is_depn], "vb_info[%u] out if bounds", is_depn);

        vb_i = vb_info[is_depn]->vb_i; // now we know what the true gencomp vb_i is

        *pi = (ReconPlanItem) {
            .vb_i       = vb_i,  // Re-assemble a line from a prim/depn: vb=vb_i
            .start_line = gc_vb_line_i[is_depn]++, // line within prim/depn vb
            .num_lines  = 1
        };

        UPDATE_WRITERS (vb_i);

        // last line of this PRIM/DEPN VB - move to next VB
        if (! (--vb_info[is_depn]->num_lines)) {
            
            // add PLAN_END_OF_VB after last line of the prim/depn VB
            *INSERBtxtAFTER (ReconPlanItem, txt_file->recon_plan, i) = (ReconPlanItem){
                .flavor = PLAN_END_OF_VB,
                .vb_i   = vb_i,
            };

            UPDATE_WRITERS_END_VB (vb_i);
            
            i++; // skip inserted PLAN_END_OF_VB

            vb_info[is_depn]++; // move to next VB
            gc_vb_line_i[is_depn] = 0;
        }
    }

    buf_free (evb->scratch);

    // merge consecutive GC lines from the same GC & VB
    sam_zip_recon_plan_optimize_gc_lines();

    if (flag.show_memory) iprintf ("\nconcurrent_writer_vblocks (consumes memory in decompression, not threads)=%u\n\n", max_conc_writers);

    COPY_TIMER_VB (evb, sam_zip_recon_plan_add_gc_lines);

    return max_conc_writers;
}

// ZIP main thread
void sam_zip_generate_recon_plan (void)
{
    START_TIMER;

    // case: recon_plan is just FULL_VBs, meaning we have no PRIM or DEPN lines. 
    // This file will have just one (main) component with no recon_plan
    if (sam_zip_recon_plan_full_vb_only()) 
        flag.bind = BIND_NONE; // no PRIM or DEPN lines found - single-component SAM/BAM without a recon_plan

    // case: MAIN component (not all full VBs - we have some PRIM and/or DEPN lines) - plan that incorporates everything 
    else {
        // recon_plan may have VBs out of order - fix that now
        recon_plan_sort_by_vb (txt_file); 

        // update recon plan with gencomp lines
        uint32_t conc_writing_vbs = sam_zip_recon_plan_add_gc_lines();
        
        // output the SEC_RECON_PLAN section
        recon_plan_compress (conc_writing_vbs, false);
    }
    
    COPY_TIMER_VB (evb, generate_recon_plan);
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
    sam_stats_reaccount_one (SAM_CIGAR,   OPTION_SA_CIGAR); // counts.param accumulated sam_cigar_seg_prim_cigar
    sam_stats_reaccount_one (SAM_RNAME,   OPTION_SA_RNAME);
    sam_stats_reaccount_one (SAM_MAPQ,    OPTION_SA_MAPQ);
    sam_stats_reaccount_one (SAM_POS,     OPTION_SA_POS);    
    sam_stats_reaccount_one (SAM_FLAG,    OPTION_SA_STRAND);    

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
    rom names[] = SAM_SAG_TYPE_NAMES;

    if (sagt < 0 || sagt >= NUM_SAG_TYPES) return "InvalidSagType";
    else                                   return names[sagt];
}
