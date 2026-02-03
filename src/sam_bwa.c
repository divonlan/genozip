// ------------------------------------------------------------------
//   sam_bwa.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "chrom.h"
#include "lookback.h"
#include "b250.h"
#include "zip_dyn_int.h"

// -------------------------------------------------------
// XA:Z "Alternative alignments" (a BWA, gem3 feature)
// XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
// Example XA:Z:chr9,-60942781,150M,0;chr9,-42212061,150M,0;chr9,-61218415,150M,0;chr9,+66963977,150M,1;
// See: http://bio-bwa.sourceforge.net/bwa.shtml
// -------------------------------------------------------

// Lookup buffer
#define lookback_value last_value.i

static void sam_seg_BWA_XA_initialize (VBlockSAMP vb)
{
    ContextP rname_ctx    = CTX(OPTION_XA_RNAME);
    ContextP strand_ctx   = CTX(OPTION_XA_STRAND);
    ContextP pos_ctx      = CTX(OPTION_XA_POS);
    ContextP lookback_ctx = CTX(OPTION_XA_LOOKBACK); // invalid if lookback_did_i=DID_NONE, that's ok
 
    dyn_int_init_ctx (VB, pos_ctx, 0);

    if (!segconf_running && segconf.is_sorted) {

        // note: we need to allocate lookback even if reps_per_line=0, lest an XA shows up despite not being in segconf
        rname_ctx->no_stons         = true;  // as we store by index
        strand_ctx->no_stons        = true;
        dyn_int_init_ctx (VB, lookback_ctx, 0);
        lookback_ctx->local_param   = true;
        lookback_ctx->local.prm8[0] = lookback_size_to_local_param (1024);
        lookback_ctx->local_always  = (lookback_ctx->local.param != 0); // no need for a SEC_LOCAL section if the parameter is 0 (which is the default anyway)

        lookback_init (VB, lookback_ctx, rname_ctx,  STORE_INDEX); // lookback_ctx->local.param must be set before
        lookback_init (VB, lookback_ctx, strand_ctx, STORE_INDEX);
        lookback_init (VB, lookback_ctx, pos_ctx,    STORE_INT);

        // create strand nodes (nodes will be deleted in sam_zip_after_vbs if not used)
        ctx_create_node (VB, OPTION_XA_STRAND, cSTR("-"));  // word_index=0
        ctx_create_node (VB, OPTION_XA_STRAND, cSTR("+"));  // word_index=1
        ctx_create_node (VB, OPTION_XA_STRAND, cSTR("-C")); // word_index=2. these 4 are used by gem3 mapper with --bisulfite-conversion
        ctx_create_node (VB, OPTION_XA_STRAND, cSTR("+C")); // word_index=3
        ctx_create_node (VB, OPTION_XA_STRAND, cSTR("-G")); // word_index=4
        ctx_create_node (VB, OPTION_XA_STRAND, cSTR("+G")); // word_index=5
    }

    // case: we're not going to use lookback - create an all-the-same XA_LOOKBACK dict
    else
        ctx_create_node (VB, OPTION_XA_LOOKBACK, "0", 1);   // word_index=0

    CTX(OPTION_XA_Z)->is_initialized = true;
}

#define SET_XA(old,new) ({ thool expected = (old); \
                           __atomic_compare_exchange_n (&segconf.sam_has_BWA_XA_Z, &expected, (new), false, __ATOMIC_RELAXED, __ATOMIC_RELAXED); })

static bool sam_seg_verify_BWA_XA (VBlockSAMP vb, STRp(xa))
{
    str_split (xa, xa_len, 0, ';', rep, false); // verify at least one semicolon
    if (n_reps < 2) 
        SET_XA(unknown, no); // set only if not already set by another thread

    else { // test first repeat
        if (segconf.sam_semcol_in_contig && !memchr (reps[0], ',', rep_lens[0])) 
            rep_lens[0] += rep_lens[1] + 1; // if first rep is a half-contig, merge it with second rep
        
        str_split (reps[0], rep_lens[0], 4, ',', item, true);

        // verify 4 items, and POS is a -/+ followed by an integer. note: GEM3 XA format will be a "no". TO DO: recognize GEM3 format 
        if (n_items != 4 || item_lens[1] < 2           || // 4 items, item [1] has at least 2 characters
            (items[1][0] != '+' && items[1][0] != '-') || // item [1] must start with a + or -
            !str_is_int (&items[1][1], item_lens[1]-1) || // rest of item [1] must be an integer (POS)
            !sam_is_cigar (STRi(item,2), false)        || // item [2] must be a valid textual CIGAR
            !str_is_int (STRi(item,3)))                   // item [3] must be an integer (NM)
            
            SET_XA(unknown, no);
        else
            SET_XA(unknown, yes);
    }

    return (segconf.sam_has_BWA_XA_Z == yes); // as set by this thread, or perhaps another thread that beat us to it
}

static bool sam_seg_BWA_XA_strand_cb (VBlockP vb, ContextP ctx, STRp(strand), uint32_t rep)
{
    if (strand_len != 1 || (*strand != '+' && *strand != '-'))  
        return false; // invalid XA format - expecting pos to begin with the strand

    // seg normally for now, we might change it in sam_seg_BWA_XA_pos
    WordIndex strand_index = (*strand == '+');
    seg_known_node_index (vb, CTX (OPTION_XA_STRAND), strand_index, 1); 

    return true; // segged successfully
}

static bool sam_seg_BWA_XA_pos_cb (VBlockP vb, ContextP ctx, STRp(pos_str), uint32_t rep)
{
    #define MAX_POS_DISTANCE 10000 // the smaller it is, the more we search for a better XA - the delta-pos will be better, but lookback worse. 10K works well.
    START_TIMER;

    ContextP rname_ctx  = CTX (OPTION_XA_RNAME);
    ContextP pos_ctx    = CTX (OPTION_XA_POS);
    ContextP strand_ctx = CTX (OPTION_XA_STRAND);
    ContextP lb_ctx     = CTX (OPTION_XA_LOOKBACK);

    // look back for a node with this index and a similar POS - we use word_index to store the original rname_node_index, pos
    WordIndex rname_index = b250_seg_get_last (rname_ctx);
    int64_t lookback = 0;
    PosType32 pos = -MAX_POS_DISTANCE; // initial to "invalid pos" - value chosen so we can store it in poses, in lieu of an invalid non-integer value, without a future pos being considered close to it

    WordIndex strand_wi = b250_seg_get_last (strand_ctx);

    if (str_get_int_range32 (STRa(pos_str), 0, MAX_POS_SAM, &pos)) {

        int64_t iterator = -1;
        while ((lookback = lookback_get_next (vb, lb_ctx, rname_ctx, rname_index, &iterator))) {

            PosType32 lookback_pos = lookback_get_value (vb, lb_ctx, pos_ctx, lookback).i;

            // case: we found a lookback - same rname and close enough pos
            if (ABS (pos-lookback_pos) < MAX_POS_DISTANCE) {
                //    if (ABS (pos-lookback_pos) <= ABS(CTX(SAM_TLEN)->last_value.i)) { <-- better POS deltas but bigger index - its a wash
                // replace rname with lookback
                b250_seg_remove_last (vb, rname_ctx, rname_index);
                seg_by_ctx (vb, STRa(XA_lookback_snip), rname_ctx, 0); // add_bytes=0 bc we already added them when we segged rname the first time

                // possibly replace strand with lookback
                WordIndex lookback_strand_index = lookback_get_index (vb, lb_ctx, strand_ctx, lookback);
                if (strand_wi == lookback_strand_index) {
                    b250_seg_remove_last (vb, strand_ctx, lookback_strand_index);
                    seg_by_ctx (vb, STRa(XA_lookback_snip), strand_ctx, 0); // add_bytes=0 bc we already added them when we segged strand the first time
                }

                // seg pos as a delta
                SNIP(48);
                memcpy (snip, XA_lookback_snip, XA_lookback_snip_len);
                snip_len = XA_lookback_snip_len + str_int (pos - lookback_pos, &snip[XA_lookback_snip_len]);
                seg_by_ctx (vb, STRa(snip), pos_ctx, pos_str_len);
                
                break;
            }
        }
    }

    if (!lookback) 
        seg_integer_or_not (vb, pos_ctx, STRa(pos_str), pos_str_len);

    dyn_int_append (vb, lb_ctx, lookback, 0);

    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_RNAME,  false, (int64_t)rname_index);
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_POS,    false, (int64_t)pos);
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_STRAND, false, (int64_t)strand_wi);

    COPY_TIMER (sam_seg_BWA_XA_pos);        
    return true; // segged successfully
}

// used in case of non-sorted files - no lookback as we are not expecting near XAs to be similar
static bool sam_seg_SA_no_lookback_cb (VBlockP vb, ContextP ctx, STRp(snip), uint32_t rep)
{
    seg_known_node_index (vb, ctx, 0, 0); // node created in sam_seg_BWA_XA_initialize
    return true; // segged successfully
}


// Called after splitting in case contig names may contain a ; - in which case we need to correct over-splitting
// Example: XA:Z:ScxkALA_1000;HRSCAF=1302,+101183,26S19M104S,0;ScxkALA_1805;HRSCAF=2562,+35201024,17S19M111S,0
// Original split: [0]=ScxkALA_1000  [1]=HRSCAF=1302,+101183,26S19M104S,0  [2]=ScxkALA_1805  [3]=HRSCAF=2562,+35201024,17S19M111S,0
// Correct to:     [0]=ScxkALA_1000;HRSCAF=1302,+101183,26S19M104S,0  [1]=ScxkALA_1805;HRSCAF=2562,+35201024,17S19M111S,0
void sam_seg_correct_for_semcol_in_contig (uint32_t *n_repeats, rom *repeats, uint32_t *repeat_lens)
{
    for (uint32_t i=0; i < *n_repeats-1; i++)  // -1 bc last repeat cannot be a half contig
        if (!memchr (repeats[i], ',', repeat_lens[i])) { // a repeat containing no commas is a "half contig" (eg "ScxkALA_1000")
            // extend to include the semicolon and the subsequent repeat
            repeat_lens[i] += 1 + repeat_lens[i+1]; 
            
            // eliminate the subsequent repeat
            memmove (&repeats[i+1],     &repeats[i+2],     (*n_repeats - i - 2) * sizeof (rom));
            memmove (&repeat_lens[i+1], &repeat_lens[i+2], (*n_repeats - i - 2) * sizeof (uint32_t));
            (*n_repeats)--;
        }
}

void sam_seg_BWA_XA_Z (VBlockSAMP vb, STRp(xa), unsigned add_bytes)
{
    START_TIMER;

    // case: still "unknown". set to "yes" or "no" and proceed accordingly. this can happen in any VB.
    if (segconf.sam_has_BWA_XA_Z == unknown && !sam_seg_verify_BWA_XA (vb, STRa(xa))) {
        DO_ONCE memcpy (segconf.sam_malformed_XA, xa, MIN_(xa_len, sizeof (segconf.sam_malformed_XA)));

        seg_by_did (VB, STRa(xa), OPTION_XA_Z, add_bytes); // fallback;        
        return;
    }

    // case: first encounter with XA in this VB, possibly after this thread or another thread set "sam_has_BWA_XA_Z" to "yes" in sam_seg_verify_BWA_XA 
    if (!CTX(OPTION_XA_Z)->is_initialized)
        sam_seg_BWA_XA_initialize (vb);

    static const MediumContainer container_XA = {
        .repeats      = 0, 
        .nitems_lo    = 6, 
        .repsep       = {';'}, // including last item
        .items        = { { .dict_id = { _OPTION_XA_LOOKBACK }, .separator = { CI0_INVISIBLE, CI1_LOOKBACK } }, 
                          { .dict_id = { _OPTION_XA_RNAME    }, .separator = { ',',           CI1_LOOKBACK } }, 
                          { .dict_id = { _OPTION_XA_STRAND   }, .separator = { CI0_DIGIT,     CI1_LOOKBACK } },
                          { .dict_id = { _OPTION_XA_POS      }, .separator = { ',',           CI1_LOOKBACK } },
                          { .dict_id = { _OPTION_XA_CIGAR    }, .separator = { ','}                          },
                          { .dict_id = { _OPTION_XA_NM       },                                              } }  };

    // case: for a collated (or otherwise unsorted) file, we can lookup against our mate, but if we have no
    // mate, we just seg without lookback, because XA:Z in near lines are not expected to be similar (also: in segconf.running)
    bool use_lb = segconf.is_sorted && !segconf_running;

    SegCallback callbacks_no_lb[6] = { sam_seg_SA_no_lookback_cb, chrom_seg_cb, 0, seg_integer_or_not_cb, sam_seg_0A_cigar_cb, 0 };

    SegCallback callbacks[6] = { seg_do_nothing_cb, chrom_seg_cb, 
                                 (MP(GEM3) && segconf.sam_bisulfite) ? sam_seg_gem3_XA_strand_cb : sam_seg_BWA_XA_strand_cb, 
                                 sam_seg_BWA_XA_pos_cb, sam_seg_0A_cigar_cb, 0 };

    int32_t repeats = seg_array_of_struct (VB, CTX(OPTION_XA_Z), container_XA, STRa(xa), 
                                           use_lb ? callbacks : callbacks_no_lb, 
                                           segconf.sam_semcol_in_contig ? sam_seg_correct_for_semcol_in_contig : NULL,
                                           add_bytes);

    // case: we failed to seg as a container - flush lookbacks (rare condition, and complicated to rollback given the round-robin and unlimited repeats)
    if (use_lb && repeats == -1) {
        SET_XA (yes, no); // this file's XA:Z is not bwa format after all 
        lookback_flush (VB, &container_XA);
    }

    COPY_TIMER(sam_seg_BWA_XA_Z);
}

// PIZ up to v13: this is called for XA that are a container, but not for invalid XA that are segged as a simple snip
void sam_piz_XA_field_insert_lookback_v13 (VBlockP vb)
{
    if (!CTX(OPTION_XA_LOOKBACK)->is_initialized) {
        lookback_init (vb, CTX(OPTION_XA_LOOKBACK), CTX (OPTION_XA_RNAME),  (CTX(OPTION_XA_RNAME) ->flags.store = STORE_INDEX));
        lookback_init (vb, CTX(OPTION_XA_LOOKBACK), CTX (OPTION_XA_STRAND), (CTX(OPTION_XA_STRAND)->flags.store = STORE_INDEX));
        lookback_init (vb, CTX(OPTION_XA_LOOKBACK), CTX (OPTION_XA_POS),    (CTX(OPTION_XA_POS)   ->flags.store = STORE_INT  ));
    }

    // copy last_value to lookback buffer
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_RNAME,  true, (int64_t)0);
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_STRAND, true, (int64_t)0);
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_POS,    true, (int64_t)0);
}

// ----------------------------------------------------------------------------------------------------------
// XC:i (bwa) Undocumented: usually seq_len minus the final soft-clip (right if forward and left if rev-comp) 
// ----------------------------------------------------------------------------------------------------------
void sam_seg_BWA_XC_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t XC, unsigned add_bytes)
{
    decl_ctx (OPTION_XC_i);
    int64_t prediction = dl->SEQ.len - vb->soft_clip[!dl->FLAG.rev_comp || dl->FLAG.unmapped];

    if (XC == prediction) 
        seg_special0 (VB, SAM_SPECIAL_BWA_XC, ctx, add_bytes);

    // case: prediction failed: cases: 1. not bwa 2. rev_comp != next_rev_comp
    else if (ABS(XC - prediction) < 1000) {  // if its not to big - seg a delta to avoid creating a local section for rare cases
        SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_BWA_XC, XC - prediction);
        seg_by_ctx (VB, STRa(snip), ctx, add_bytes);
    }

    else {
        seg_integer (VB, ctx, XC, true, add_bytes);
    }

    ctx_set_encountered (VB, ctx); // needed for segconf xcons_std_line_histogram building to work
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_BWA_XC)
{
    int64_t delta=0;
    if (snip_len)  // we have a delta
        ASSPIZ (str_get_int (STRa(snip), &delta), "Invalid delta=\"%.*s\"", STRf(snip));

    new_value->i = vb->seq_len - VB_SAM->soft_clip[!last_flags.rev_comp || last_flags.unmapped] + delta;
    
    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// XT:A (bwa) Unique/Repeat/N/Mate-sw
// ----------------------------------------------------------------------------------------------------------
void sam_seg_BWA_XT_A (VBlockSAMP vb, char XT, unsigned add_bytes)
{
    int64_t X0 = has(X0_i) ? sam_seg_get_aux_int_(vb, X0_i) : 0;
    
    // predict based on the existance and value of X0
    char prediction = !has(X0_i) ? 'M'  // Mate Smith-Waterman used for mapping. Note: could also be 'N', but we set our prediction to 'M' as it is very much more common
                    : X0 == 1    ? 'U'  // Unique mapping
                    :              'R'; // Repeat (i.e. not unique)

    if (prediction == XT)
        seg_special0 (VB, SAM_SPECIAL_BWA_XT, CTX(OPTION_XT_A), add_bytes);

    else
        seg_by_did (VB, (char[]){ XT }, 1, OPTION_XT_A, add_bytes);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_BWA_XT)
{
    if (reconstruct) {
        char prediction = !curr_container_has (vb, _OPTION_X0_i)   ? 'M'
                        : reconstruct_peek (vb, CTX(OPTION_X0_i), 0, 0).i == 1 ? 'U'
                        :                                                        'R';
        RECONSTRUCT1 (prediction);
    }

    return NO_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// X1:A (bwa) Number of suboptimal hits found by BWA
// ----------------------------------------------------------------------------------------------------------
void sam_seg_BWA_X1_i (VBlockSAMP vb, int64_t X1, unsigned add_bytes)
{
    if (!has(X0_i)) goto fallback; // need X0 to calculate prediction

    int64_t prediction = 0;

    // in some files, where XA:Z includes sub-optimal alignments, this is usually true: 
    // (X0 + X1 = 1 + XA.repeats), and (X0 >= 1)
    if (has(XA_Z)) {
        STR(xa);
        sam_seg_get_aux_Z (vb, vb->idx_XA_Z, pSTRa(xa), IS_BAM_ZIP);

        int xa_alns = str_count_char (STRa(xa), ';');

        prediction = xa_alns + 1 - sam_seg_get_aux_int_(vb, X0_i);;
    }
    
    if (prediction == X1)
        seg_special0 (VB, SAM_SPECIAL_BWA_X1, CTX(OPTION_X1_i), add_bytes);

    // note: rarely, X1 can be very large (tens of thousands)
    else 
        fallback:
        seg_integer (VB, CTX(OPTION_X1_i), X1, true, add_bytes);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_BWA_X1)
{
    new_value->i = 0; // default if no XA

    if (curr_container_has (vb, _OPTION_XA_Z)) {
        int64_t XA_repeats = container_peek_repeats (vb, CTX(OPTION_XA_Z), ';');
        int64_t X0 = reconstruct_peek (vb, CTX(OPTION_X0_i), 0, 0).i;
        
        new_value->i = XA_repeats + 1 - X0;
    }

    if (reconstruct) 
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------
// Suboptimal alignment score : XS:i in bwa, dragen, bsbolt, tmap and bowtie2 ; ZS:i in hisat2
// ----------------------------------------------------------------------------------------------

static inline int sam_XS_get_mux_channel (uint8_t MAPQ)
{
    switch (MAPQ) {
        case 0         : return 0;
        case 1  ... 10 : return 1;
        case 11 ... 59 : return 2;
        default        : return 3;
    }    
}

// usually XS:i, but ZS:i in hisat2
void sam_seg_BWA_XS_i (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, int64_t xs, unsigned add_bytes)
{
    START_TIMER;

    // "Suboptimal alignment score" - multiplex by MAPQ and (sometimes) delta vs AS
    if (has(AS_i) && ABS(xs) < 10000) {

        int channel_i = sam_XS_get_mux_channel (dl->MAPQ);
        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, OPTION_XS_i, (MultiplexerP)&vb->mux_XS, channel_i);

        if (channel_i == 3 ||
            segconf.AS_is_2ref_consumed) { // when AS:i is inflated x2 we observe that XS:i is not very near AS:i or AS:i/2 - we're better of segging in local

            seg_integer (VB, channel_ctx, xs, true, add_bytes);
        }

        else if (xs && dl->AS) { // XS can be before or after AS
            CTX(OPTION_AS_i)->last_value.i = dl->AS;
            seg_delta_vs_other_localN (VB, channel_ctx, CTX(OPTION_AS_i), xs, -1, add_bytes);
        }

        else if (!xs)
            seg_by_ctx (VB, "0", 1, channel_ctx, add_bytes);

        else
            seg_integer_as_snip_do (VB, channel_ctx, xs, add_bytes);

        seg_by_did (VB, STRa(vb->mux_XS.snip), did_i, 0);
    }
    
    // delta doesn't make sense - store as snip
    else
        seg_integer (VB, CTX(did_i), xs, true, add_bytes);

    COPY_TIMER (sam_seg_BWA_XS_i);
}

// v14: De-multiplex XS by MAPQ
SPECIAL_RECONSTRUCTOR (sam_piz_special_BWA_XS)
{
    int channel_i = sam_XS_get_mux_channel (CTX(SAM_MAPQ)->last_value.i);
    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}
