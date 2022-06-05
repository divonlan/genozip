// ------------------------------------------------------------------
//   sam_gc_zip.c
//   Copyright (C) 2022-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

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
#include "bit_array.h"
#include "writer.h"
#include "sections.h"
#include "codec.h"
#include "compressor.h"
#include "bit_array.h"
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

//------------------------
// Zip - After Reading VB
//------------------------

// called by main thread after reading a VB - callback of zip_init_vb
void sam_gc_zip_init_vb (VBlockP vb)
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

void sam_seg_gc_initialize (VBlockSAMP vb)
{
    // DEPN stuff
    if (sam_is_depn_vb) {
        CTX(SAM_SAGROUP)->ltype = sizeof(SAGroup)==4 ? LT_UINT32 : LT_UINT64;
        CTX(SAM_SAALN)->ltype   = LT_UINT16; // index of alignment with SA Group
    }

    // PRIM stuff 
    else if (sam_is_prim_vb) { 
        CTX(SAM_SAALN)->ltype              = LT_UINT16;   // index of alignment with SA Group
        CTX(OPTION_SA_Z)->ltype            = LT_UINT8;    // we store num_alns in local
        
        // set store, consumed by sam_piz_prim_add_Alns
        CTX(OPTION_SA_RNAME)->flags.store  = STORE_INDEX;
        CTX(OPTION_SA_STRAND)->flags.store = STORE_INDEX; // 1 for - and 0 for +
        CTX(OPTION_SA_POS)->flags.store    = STORE_INT;
        CTX(OPTION_SA_MAPQ)->flags.store   = STORE_INT;
        CTX(OPTION_SA_NM)->flags.store     = STORE_INT;

        CTX(SAM_FLAG)->no_stons            = true;        // 
    }
}

// Note: in MAIN vbs, failure means we will keep this line in MAIN and not move it to PRIM
// In PRIM, we abort because this is not not expected as MAIN component should not have moved this line to PRIM
#define FAILIF(condition, format, ...) \
  ( { if (condition) { \
          if (flag.debug_sa) iprintf ("%s: " format "\n", LN_NAME, __VA_ARGS__); \
          if (sam_is_main_vb) return false; \
          else { progress_newline(); fprintf (stderr, "%s: Error in %s:%u: Failed PRIM line because ", LN_NAME, __FUNCLINE); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, SUPPORT); fflush (stderr); exit_on_error(true); }} \
    } )

// Call in seg PRIM line (both in MAIN and PRIM vb): 
// MAIN: test that line is a valid prim before moving it to the PRIM gencomp file
// PRIM: add SA Group (not alignments) to the data structure in VB 
static bool sam_seg_prim_add_sa_group (VBlockSAMP vb, ZipDataLineSAM *dl, uint16_t num_alns/* inc primary aln*/, bool is_bam)
{
    rom textual_seq = is_bam ? B1STc(vb->textual_seq) : Btxt (dl->SEQ.index);
    uint32_t seq_len = is_bam ? vb->textual_seq.len32 : dl->SEQ.len;

    ASSERTNOTNULL (textual_seq);

    // MAIN: check that values are within limits defined in SAGroupType (no need to check in PRIM as we already checked in MAIN)
    if (sam_is_main_vb) {
        FAILIF (dl->QNAME.len > MAX_SA_QNAME_LEN, "dl->QNAME.len=%u > %u", dl->QNAME.len, MAX_SA_QNAME_LEN);
        FAILIF (seq_len==1 && *textual_seq == '*', "SEQ=\"*\"%s", ""); // we haven't segged seq yet, so vb->seq_missing is not yet set
        FAILIF (seq_len > MAX_SA_SEQ_LEN, "seq_len=%u > %u", seq_len, MAX_SA_SEQ_LEN);
        FAILIF (!dl->POS, "POS=0%s", ""); // unaligned

        uint32_t bad_i;
        FAILIF (!str_is_only_ACGT (textual_seq, seq_len, &bad_i), 
                "RNAME=\"%.*s\" POS=%d. Found '%c' in base_i=%u in SEQ is not A,C,G or T. SEQ=\"%.*s\"", 
                vb->chrom_name_len, vb->chrom_name, dl->POS, textual_seq[bad_i], bad_i, seq_len, textual_seq);
    }

    // PRIM: actually add the group
    else if (sam_is_prim_vb) { 
        buf_alloc (vb, &vb->sa_groups, 1, 32, SAGroupType, CTX_GROWTH, "sa_groups");

        // add SA group
        BNXT (SAGroupType, vb->sa_groups) = (SAGroupType){
            .first_aln_i    = vb->sa_alns.len,
            .num_alns       = num_alns,
            .qname          = dl->QNAME.index,
            .qname_len      = dl->QNAME.len,
            .qual           = dl->QUAL.index,
            .no_qual        = vb->qual_missing,
            .revcomp        = dl->FLAG.bits.rev_comp,
            .multi_segments = dl->FLAG.bits.multi_segments,
            .is_first       = dl->FLAG.bits.is_first,
            .is_last        = dl->FLAG.bits.is_last,
            .seq            = dl->SEQ.index,
            .seq_len        = seq_len,
        };
    }

    return true;
}

// Call in seg PRIM line (both in MAIN and PRIM vb): 
// MAIN: test that line is a valid prim before moving it to the PRIM gencomp file
// PRIM: add SA Group (including alignments) to the data structure in VB (called from sam_SA_Z_seg)
bool sam_seg_prim_add_sa_group_SA (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(sa), int64_t this_nm, bool is_bam)
{
    str_split (sa, sa_len, 0, ';', aln, false); // n_alns is 1+number of repeats in SA, because of terminating ';'. so it is the total number of alignments, including the primary
    FAILIF (n_alns < 2 || n_alns > MAX_SA_NUM_ALNS, "n_alns=%d is not in [2,%u]", n_alns, MAX_SA_NUM_ALNS); // we cannot add this SA - either invalid or over the max number of alignments

    if (!sam_seg_prim_add_sa_group (vb, dl, n_alns, is_bam)) return false; 

    // check that values are within limits defined in SAAlnType
    FAILIF (dl->POS  > MAX_SA_POS, "POS=%d > MAX_SA_POS=%u", dl->POS, MAX_SA_POS);
    FAILIF (dl->MAPQ > MAX_SA_MAPQ, "MAPQ=%u > MAX_SA_MAPQ=%u", dl->MAPQ, MAX_SA_MAPQ);
    FAILIF (this_nm  > MAX_SA_NM, "NM=%"PRId64" > MAX_SA_NM=%u", this_nm, MAX_SA_NM);
    FAILIF (vb->cigar_missing, "No CIGAR%s", "");
    FAILIF (vb->hard_clip[0], "hard_clip[LEFT]=%u > 0", vb->hard_clip[0]); // we can't add primary lines with a CIGAR with H as we need the full sequence and qual 
    FAILIF (vb->hard_clip[1], "hard_clip[RIGHT]=%u > 0", vb->hard_clip[1]);

    if (sam_is_prim_vb) { // in PRIM, we actually add the group, in MAIN, we are just testing
        buf_alloc (vb, &vb->sa_alns, n_alns /*+1 for primary aln*/, 64, SAAlnType, CTX_GROWTH, "sa_alns");

        // in Seg, we add the Adler32 of CIGAR to save memory, as we just need to verify it, not reconstruct it
        uint32_t cigar_len = (is_bam ? vb->textual_cigar.len32 : dl->CIGAR.len);
        CigarSignature cigar_sig = cigar_sign (is_bam ? B1STc(vb->textual_cigar) : Btxt(dl->CIGAR.index), cigar_len);

        // add primary alignment as first alignment in SA group
        BNXT (SAAlnType, vb->sa_alns) = (SAAlnType){
            .rname           = vb->chrom_node_index,
            .pos             = dl->POS,
            .revcomp         = dl->FLAG.bits.rev_comp,
            .cigar.signature = cigar_sig,
            .mapq            = dl->MAPQ,
            .nm              = this_nm
        };

        // in BAM, add textual CIGAR to vb->sa_cigars, as txt_data CIGAR is binary
        if (is_bam)
            buf_add_buf (vb, &vb->sa_prim_cigars, &vb->textual_cigar, char, "sa_prim_cigars");
    }

    // add all dependent alignments
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

            CigarSignature cigar_sig = cigar_sign (STRi(item,SA_CIGAR));

            BNXT (SAAlnType, vb->sa_alns) = (SAAlnType){ 
                .rname           = rname, 
                .pos             = pos, 
                .revcomp         = revcomp, 
                .cigar.signature = cigar_sig,
                .mapq            = mapq,
                .nm              = nm
            };
        }
    }

    return true;
}       

// Seg MAIN: returns true if this line is a Primary or Dependent line of a supplementary/secondary group - 
// and should be moved to a generated component. Called by compute thread in seg of Normal VB.
bool sam_seg_is_gc_line (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(alignment), STRps(aux), bool is_bam)
{
    if (vb->chrom_node_index == WORD_INDEX_NONE) return false; // RNAME='*' or headerless SAM (we need header-contigs for sam_sa_add_sa_group)

    int32_t HI, NH;
    SamNMType NM; 
    SamComponentType comp_i = 0;;
    
    // dependent - always has secondary or supplementary flag
    if (dl->FLAG.bits.secondary || dl->FLAG.bits.supplementary) 
        comp_i = SAM_COMP_DEPN;

    // primary - no sec/sup flag, identify by SA, and very all the fields
    else if (vb->idx_SA_Z != -1 && vb->idx_NM_i != -1 &&
             sam_seg_get_aux_int (vb, STRi(aux, vb->idx_NM_i), &NM, is_bam, MIN_NM_i, MAX_NM_i, true) &&
             sam_seg_prim_add_sa_group_SA (vb, dl, STRauxZ (SA_Z, is_bam), NM, is_bam))  // testing to see if we can successfully add an SA Group based on SA
        comp_i = SAM_COMP_PRIM;

    // primary - no sec/sup flag, identify by NH>1 && HI==1
    // note: some files have NH but no HI. This usually have no supplementary alignments in the file.
    else if (vb->idx_NH_i != -1 && vb->idx_HI_i != -1 &&
             sam_seg_get_aux_int (vb, STRi(aux, vb->idx_NH_i), &NH, is_bam, 2, MAX_HI_NH, true) &&
             sam_seg_get_aux_int (vb, STRi(aux, vb->idx_HI_i), &HI, is_bam, 1, 1, true) && 
             sam_seg_prim_add_sa_group (vb, dl, 0, is_bam)) // testing to see if we can successfully add an SA Group based on NH/HI
        comp_i = SAM_COMP_PRIM;

    else
        return false; // not a PRIM or DEPN line

    // store location where this gc line should be inserted    
    gencomp_seg_add_line (VB, comp_i, STRa(alignment));
    
    vb->line_i--;
    vb->recon_size -= alignment_len;
    vb->txt_size   -= alignment_len;

    return true;
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

// ZIP DEPN: find group index with this_qname in z_file->sa_groups, and if there are several - return the first
static SAGroupType *sam_sa_get_first_group_by_qname_hash (VBlockSAMP vb, STRp(this_qname), bool is_last, int64_t *grp_index_i)
{
    // search for a group with qname in z_file->sa_qname
    uint32_t this_qname_hash = QNAME_HASH (this_qname, this_qname_len, is_last);
    *grp_index_i = sam_sa_binary_search_for_qname_hash (B1ST (SAGroupIndexEntry, z_file->sa_groups_index), this_qname_hash, 0, z_file->sa_groups_index.len-1);

    const SAGroupIndexEntry *index_ent = B(SAGroupIndexEntry, z_file->sa_groups_index, *grp_index_i); // invalid pointer if grp_index_i==-1, that's ok    
    return (*grp_index_i >= 0) ? B(SAGroupType, z_file->sa_groups, index_ent->grp_i) : NULL; 
}

// ZIP DEPN: if there are more groups in z_file->sa_groups with the qname adlers, return the next group index
static SAGroupType *sam_sa_get_next_group_by_qname_hash (VBlockSAMP vb, int64_t *grp_index_i)
{
    const SAGroupIndexEntry *index_ent = B(SAGroupIndexEntry, z_file->sa_groups_index, *grp_index_i);

    if (*grp_index_i < z_file->sa_groups_index.len-1 && index_ent->qname_hash == (index_ent+1)->qname_hash) {
        (*grp_index_i)++;
        return B(SAGroupType, z_file->sa_groups, (index_ent+1)->grp_i);
    }

    return NULL; // not found
}

// Seg DEPN: check if the seq from this DEPN line matches a SA Group, considering hard-clips in its CIGAR
static bool sam_seg_depn_is_subseq_of_prim (VBlockSAMP vb, bytes depn_textual_seq, uint32_t depn_seq_len,
                                            bool xstrand, // depn_seq is opposite strand vs primary
                                            const SAGroupType *g,
                                            bool is_bam)
{
    // seq - DEPN may have hard-clips - check the DEPN sequence as a sub-sequence of PRIM
    if (g->seq_len != vb->hard_clip[0] + depn_seq_len + vb->hard_clip[1]) return false;
    
    BitArray *prim_seq = buf_get_bitarray (&z_file->sa_seq); // ACGT format
    uint32_t start_p = g->seq + vb->hard_clip[0];

    if (!xstrand)
        for (uint32_t i=0; i < depn_seq_len; i++) {
            uint8_t p = bit_array_get2 (prim_seq, (start_p + i) * 2);

            switch (depn_textual_seq[i]) {
                case 'A': if (p != 0) return false; break; // if prim base is not 'A'(=0) then the sequences don't match
                case 'C': if (p != 1) return false; break; 
                case 'G': if (p != 2) return false; break; 
                case 'T': if (p != 3) return false; break; 
                default : return false; // depn seq has a non-ACGT "base" - we can't seg it against prim
            }
        }

    else {
        uint32_t start_p = g->seq + vb->hard_clip[1];

        for (uint32_t i=0; i < depn_seq_len; i++) {
            uint8_t p = bit_array_get2 (prim_seq, (start_p + i) * 2);

            switch (depn_textual_seq[depn_seq_len - i - 1]) {
                case 'T': if (p != 0) return false; break; // since xstrand depn seq is T, we expect prim seq to be A(=0)
                case 'G': if (p != 1) return false; break; 
                case 'C': if (p != 2) return false; break; 
                case 'A': if (p != 3) return false; break; 
                default : return false; // depn seq has a non-ACGT "base" - we can't seg it against prim
            }
        }
    }

    return true;
}

// --------------------------------------------------
// Seg Alignments (main fields + SA-field alignments)
// --------------------------------------------------

// Seg DEPN line: verify that this line is a member of group (identical alignments) + return aln_i of this line
static inline bool sam_seg_depn_find_aln (VBlockSAMP vb, const SAGroupType *g,
                                          uint32_t n_my_alns, const SAAlnType *my_alns,
                                          uint16_t *my_aln_i) // out - relative to group
{
    uint32_t SA_aln_i=1; // index into my_alns
    for (uint32_t grp_aln_i=0; grp_aln_i < g->num_alns; grp_aln_i++) { 

        const SAAlnType *grp_aln = B(SAAlnType, z_file->sa_alns, g->first_aln_i + grp_aln_i);
        const SAAlnType *SA_aln = &my_alns[SA_aln_i];

        // case: group alignment matches the SA alignment
        if (SA_aln_i < n_my_alns && !memcmp (grp_aln, SA_aln, sizeof(SAAlnType))) {
            SA_aln_i++;
            continue;
        }

        // case: group alignment matches my alignment
        if (grp_aln_i && !memcmp (grp_aln, &my_alns[0], sizeof (SAAlnType))) // note: we don't allow a depn line to have my_aln_i=0, bc we make an assumption (eg in sam_reconstruct_main_cigar_from_SA_Group) that a DEPN line cannot have the primary alignment (aln_i=0) in its main fields
            *my_aln_i = grp_aln_i;

        // case: group alignment doesn't match the SA alignment and doesn't match my alignment
        else
            return false; 
    }    

    return SA_aln_i == n_my_alns; // all good if we verified all SA alignments (i.e. exactly one match of my alignment)
}

// parse my SA field, and build my alignments based on my main fields data + SA data
static inline bool sam_sa_seg_depn_get_my_alns (VBlockSAMP vb,
                                                WordIndex my_rname_contig, SamPosType my_pos, uint8_t my_mapq, 
                                                STRp(my_cigar), int64_t my_nm, bool my_revcomp, 
                                                STRp(SA), // 0,0 if no SA (eg STAR alignment) 
                                                uint32_t n_my_alns, 
                                                SAAlnType *my_alns) // out 
{
    str_split (SA, SA_len, n_my_alns, ';', SA_aln, true);

    my_alns[0] = (SAAlnType){ .cigar.signature = cigar_sign (STRa(my_cigar)),
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

        my_alns[aln_i] = (SAAlnType){
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

// Seg of DEPN: find matching SAGroup and seg SAM_SAGROUP. We identify the group by matching (qname,rname,pos,cigar,revcomp,nm) 
// sets vb->sa_grp (NULL if not found).
static void sam_sa_seg_depn_find_sagroup (VBlockSAMP vb, ZipDataLineSAM *dl, STRps(aux), 
                                          STRp(textual_cigar), rom textual_seq, bool is_bam)
{
    #define DONE ({ sam_cigar_restore_H (htos); return; })

    vb->sa_grp = NULL; // initialize to "not found"
    vb->sa_aln = NULL;

    bool revcomp = dl->FLAG.bits.rev_comp;
    uint32_t seq_len = dl->SEQ.len;

    buf_alloc (vb, &CTX(SAM_SAGROUP)->local, 1, vb->lines.len, SAGroup, 1, "contexts->local");
    buf_alloc (vb, &CTX(SAM_SAALN)->local,   1, vb->lines.len, uint16_t,  1, "contexts->local");

    int64_t grp_index_i=-1;
    const SAGroupType *g = sam_sa_get_first_group_by_qname_hash (vb, STRtxtw(dl->QNAME), dl->FLAG.bits.is_last, &grp_index_i);
    if (!g) return; // no PRIM with this qname

    // temporarily replace H with S
    HtoS htos = sam_cigar_H_to_S (vb, (char*)STRa(textual_cigar));

    // if this is BAM, and we have an odd number of bases, the final seq value in the BAM file of the "missing" base
    // we can't encode last "base" other than 0 (see also bug 531)
    if (is_bam && (seq_len&1) && (*B8 (vb->txt_data, dl->SEQ.index + seq_len/2) & 0xf)) DONE;

    // we can't encode depn seq bases other than A,C,G,T as we operate in 2bit
    if (!str_is_only_ACGT (textual_seq, seq_len, NULL)) DONE;

    uint32_t n_my_alns=1; // initialize: only main-fields alignment (eg in the case of Star)

    if (has_SA) {
        ctx_set_encountered (VB, CTX(OPTION_SA_Z));
        n_my_alns = str_count_char (STRauxZ(SA_Z, is_bam), ';') + 1; // +1 for alignment of my main fields
    }

    SamNMType my_nm = -1; // -1 means no nm
    SAAlnType my_alns[n_my_alns];

    // get alignments of this DEPN line ([0]=main field [1...]=SA alignments)
    if (has_SA) {
        if (has_NM) sam_seg_get_aux_int (vb, STRi(aux, vb->idx_NM_i), &my_nm, is_bam, MIN_NM_i, MAX_NM_i, false);

        if (!sam_sa_seg_depn_get_my_alns (vb, vb->chrom_node_index, dl->POS, dl->MAPQ, STRa(textual_cigar), my_nm, revcomp, 
                                          STRauxZ(SA_Z, is_bam), n_my_alns, my_alns))
            DONE; // cannot make sense of the SA alignments
    }

    // iterate on all groups with matching QNAME hash
    do {
        uint16_t my_aln_i=0; // relative to group

        // case: we have SA - get alignment where SA and my alignment perfectly match the group 
        // for SEQ - we tolerate a diff in SA alignments, but if no SA - it has to be a perfect match (modulo revcomp) 
        if (dl->FLAG.bits.is_first == g->is_first && // note: order tests in the condition is optimized to fail fast with minimal effort
            g->seq_len == vb->hard_clip[0] + seq_len + vb->hard_clip[1] &&
            vb->qual_missing == g->no_qual &&
            dl->FLAG.bits.is_last  == g->is_last &&
            dl->FLAG.bits.multi_segments == g->multi_segments &&
            str_issame_(STRtxtw(dl->QNAME), GRP_QNAME(g), g->qname_len) && 
            (!has_SA || (n_my_alns == g->num_alns && sam_seg_depn_find_aln (vb, g, n_my_alns, my_alns, &my_aln_i))) &&
            (has_SA || sam_seg_depn_is_subseq_of_prim (vb, (uint8_t*)textual_seq, dl->SEQ.len, (revcomp != g->revcomp), g, is_bam))) {

            // found - seg into local
            BNXT (SAGroup, CTX(SAM_SAGROUP)->local) = ZGRP_I(g);
            BNXT16 (CTX(SAM_SAALN)->local) = my_aln_i;
            
            vb->sa_grp = g;
            vb->sa_aln = has_SA ? B(SAAlnType, z_file->sa_alns, g->first_aln_i + my_aln_i) : NULL;

            if (flag.show_depn)
                 iprintf ("vb=%u grp=%u aln=%s\n", vb->vblock_i, ZGRP_I(g), 
                          has_SA ? sam_show_sa_one_aln (vb->sa_grp, vb->sa_aln).s : "N/A");                
            
            break;
        }

    } while ((g = sam_sa_get_next_group_by_qname_hash (vb, &grp_index_i)));

    DONE;
}

// Seg compute VB: called when segging PRIM/DEPN VBs
void sam_seg_sa_group_stuff (VBlockSAMP vb, ZipDataLineSAM *dl, 
                             STRps(aux), STRp(textual_cigar), rom textual_seq, bool is_bam)
{
    // in Dependent component - try to find which SAGroup this line belongs to, and seg it to SAM_SAGROUP and SAM_SAALN
    // (if successful, this will set vb->sa_grp/sa_aln)
    if (sam_is_depn_vb) 
        sam_sa_seg_depn_find_sagroup (vb, dl, STRas(aux), STRa(textual_cigar), textual_seq, is_bam);

    else if (sam_is_prim_vb) {
        if (has_SA) 
            ctx_set_encountered (VB, CTX(OPTION_SA_Z)); // for DEPN, this is done in sam_sa_seg_depn_find_sagroup

        // primary line without SA:Z (i.e. defined by NH:i/HI:i) - seg "line has only one alignment" 
        else { 
            // we seg the main fields into the single (primary) alignment for non-SA lines
            uint8_t num_alns_8b = 1; 
            seg_add_to_local_nonresizeable (VB, CTX(OPTION_SA_Z), &num_alns_8b, false, 0);

            // case PRIM line with no SA and no NM: we still need NM, which is used to load the primary alignment in piz. 
            // If this line has NM, it will be segged in sam_seg_NM_field. If it doesn't, we seg it as 0 here. In case of a file
            // that doesn't have NM at all (the likely case - if one line is lacking, likely all lines are) - this will just be an all-the-same context.
            // note: we don't increment OPTION_SA_NM->counts.param, as the NM:i didn't contribute to SA:Z... probably we don't even have an NM:i field in the file. We will address this in sam_stats_reallocate.
            if (!dl->NM_len) 
                seg_by_did_i (VB, "0", 1, OPTION_SA_NM, 0); 
        }
    }
}

void sam_seg_against_sa_group (VBlockSAMP vb, ContextP ctx, uint32_t add_bytes)
{
    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SAGROUP }, 2, ctx, add_bytes);
}

// seg a DEPN line against the SA data in z_file
void sam_seg_against_sa_group_int (VBlockSAMP vb, ContextP ctx, int64_t parameter, uint32_t add_bytes)
{
    SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_SAGROUP, parameter);
    seg_by_ctx (VB, STRa(snip), ctx, add_bytes);
}

void sam_seg_against_sa_group_bool (VBlockSAMP vb, ContextP ctx, bool parameter, uint32_t add_bytes)
{
    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SAGROUP, '0'+parameter }, 3, ctx, add_bytes);
}

//---------------------
// Zip - After Compute
//---------------------

// called in the main thread after_compute, in order of VBs - MAIN component
void sam_zip_gc_after_compute_main (VBlockSAMP vb)
{
    BufferP recon_plan = &txt_file->recon_plan;
    ARRAY (GencompLineIEntry, gc_lines, vb->gencomp_lines);

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
        if (normal_line_i < vb->lines.len) 
            BNXT (ReconPlanItem, *recon_plan) = (ReconPlanItem){
                .vb_i       = vb->vblock_i,
                .start_line = normal_line_i,
                .num_lines  = vb->lines.len - normal_line_i
            };

        // insert end-of-VB (note: we insert after all lines, not just normal lines, since it is at this point that
        // writer_main_loop calculates the digest)
        BNXT (ReconPlanItem, *recon_plan) = (ReconPlanItem){ 
            .vb_i   = vb->vblock_i,
            .flavor = PLAN_END_OF_VB,
        }; 
    }
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
        vb_header->sam_prim_num_alns         = BGEN32 (CTX(OPTION_SA_CIGAR)->b250.len32); // this is total number of alignments, because we add a CIGAR for each alignment including the primary alignment
        vb_header->sam_prim_first_grp_i      = BGEN32 (VB_SAM->first_grp_i);
        vb_header->sam_prim_comp_cigars_len  = BGEN32 (VB_SAM->comp_cigars_len);
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
        // update recon plan with gencomp lines
        uint32_t conc_writing_vbs = sam_zip_recon_plan_add_gc_lines();
        
        // output the SEC_RECON_PLAN section
        recon_plan_compress (conc_writing_vbs, false);
    }
    
    COPY_TIMER_VB (evb, generate_recon_plan);
}

// Callback from stats_get_compressed_sizes: in PRIM VBs, we seg main field CIGAR into OPTION_SA_CIGAR.
// We correct the z_data counts to be shown, so that the right amount of OPTION_SA_CIGAR data is shown as SAM_CIGAR data.
static void sam_stats_reaccount_one (DidIType main_ctx_did_i, DidIType sa_ctx_did_i)
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

    // Note: if all PRIM lines either have SA+NM or don't, then the math will be correct. In the unlikely case of mixed SA:Z and NH/HI PRIM lines, the math will be off.
    if (ZCTX(OPTION_NM_i)->txt_len) // we have NM:i in this file
        sam_stats_reaccount_one (OPTION_NM_i, OPTION_SA_NM);
    
    // case: no NM:i in this field, and hence also no SA:Z - all OPTION_SA_NM data is due to "0" added in sam_seg_sa_group_stuff, as we will account for it as 
    // belonging to SAM_SAGROUP
    else 
        sam_stats_reaccount_one (SAM_SAGROUP, OPTION_SA_NM); 
}
