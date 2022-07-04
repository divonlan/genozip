// ------------------------------------------------------------------
//   sam_xa.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

// -------------------------------------------------------
// XA:Z "Alternative alignments" (a BWA, gem3 feature)
//
// XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
// Example XA:Z:chr9,-60942781,150M,0;chr9,-42212061,150M,0;chr9,-61218415,150M,0;chr9,+66963977,150M,1;
// See: http://bio-bwa.sourceforge.net/bwa.shtml
//
// gem3 also uses the same XA:Z but when used with --bisulfite produce an XB:C component of POS:
// Example of gem with bisulftie: XA:Z:chr1,+G10006,103M1I36M1D11M,9;
//
// Note that a different XA was observed in IonXpress data: "XA:Z:map4-1". This will be just segged as fallback
// -------------------------------------------------------

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "chrom.h"
#include "codec.h"
#include "lookback.h"

// Lookup buffer
#define lookback_buf zip_lookback_buf // we store the previous rname, pos, strand in their ctx->zip_lookback_buf buffer
#define lookback_value last_value.i

void sam_seg_BWA_XA_initialize (VBlockSAMP vb)
{
    ContextP rname_ctx    = CTX(OPTION_XA_RNAME);
    ContextP strand_ctx   = CTX(OPTION_XA_STRAND);
    ContextP pos_ctx      = CTX(OPTION_XA_POS);
    ContextP cigar_ctx    = CTX(OPTION_XA_CIGAR);
    ContextP lookback_ctx = CTX(OPTION_XA_LOOKBACK); // invalid if lookback_did_i=DID_NONE, that's ok
 
    // note: we need to allocate lookback even if reps_per_line=0, lest an XA shows up despite not being in segconf
    if (segconf.is_sorted) {
        rname_ctx->no_stons              = true;  // as we store by index
        strand_ctx->no_stons             = true;
        lookback_ctx->flags.store        = STORE_INT;
        lookback_ctx->ltype              = LT_DYN_INT;
        lookback_ctx->local_param        = true;
        lookback_ctx->local.prm8[0]      = lookback_size_to_local_param (1024);
        lookback_ctx->local_always       = (lookback_ctx->local.param != 0); // no need for a SEC_LOCAL section if the parameter is 0 (which is the default anyway)

        lookback_init (VB, lookback_ctx, rname_ctx,  STORE_INDEX); // lookback_ctx->local.param must be set before
        lookback_init (VB, lookback_ctx, strand_ctx, STORE_INDEX);
        lookback_init (VB, lookback_ctx, pos_ctx,    STORE_INT);
    }

    // create strand nodes (nodes will be deleted in sam_seg_finalize if not used)
    ctx_create_node (VB, OPTION_XA_STRAND, cSTR("-"));  // word_index=0
    ctx_create_node (VB, OPTION_XA_STRAND, cSTR("+"));  // word_index=1
    ctx_create_node (VB, OPTION_XA_STRAND, cSTR("-C")); // word_index=2. these 4 are used by gem3 mapper with --bisulfite-conversion
    ctx_create_node (VB, OPTION_XA_STRAND, cSTR("+C")); // word_index=3
    ctx_create_node (VB, OPTION_XA_STRAND, cSTR("-G")); // word_index=4
    ctx_create_node (VB, OPTION_XA_STRAND, cSTR("+G")); // word_index=5
    strand_ctx->no_vb1_sort = true; // keep them in this ^ order
}

// Seg lookback callback for POS item of XA - seg the POS relative to lookback pos, and replace the already segged RNAME with a lookback if there is one.
static void sam_seg_XA_pos (VBlockP vb, STRp(pos_str), uint32_t rep)
{
    #define MAX_POS_DISTANCE 10000 // the smaller it is, the more we search for a better XA - the delta-pos will be better, but lookback worse. 10K works well.
    START_TIMER;

    ContextP rname_ctx = CTX (OPTION_XA_RNAME);
    ContextP pos_ctx   = CTX (OPTION_XA_POS);
    ContextP lb_ctx    = CTX (OPTION_XA_LOOKBACK);

    // look back for a node with this index and a similar POS - we use word_index to store the original rname_node_index, pos
    WordIndex rname_index = LASTb250(rname_ctx);
    int64_t lookback = 0;
    SamPosType pos = -MAX_POS_DISTANCE; // initial to "invalid pos" - value chosen so we can storeit in poses, in lieu of an invalid non-integer value, without a future pos being considered close to it

    if (!segconf.running && str_get_int_range32 (STRa(pos_str), 0, MAX_POS_SAM, &pos)) {

        int64_t iterator = -1;
        while ((lookback = lookback_get_next (vb, lb_ctx, rname_ctx, rname_index, &iterator))) {

            SamPosType lookback_pos = lookback_get_value (vb, lb_ctx, pos_ctx, lookback).i;

            // case: we found a lookback - same rname and close enough pos
            if (ABS (pos-lookback_pos) < MAX_POS_DISTANCE) {
            //    if (ABS (pos-lookback_pos) <= ABS(CTX(SAM_TLEN)->last_value.i)) { <-- better POS deltas but bigger index - its a wash
                // replace rname with lookback
                rname_ctx->b250.len--;
                ctx_decrement_count (vb, rname_ctx, rname_index);
                
                seg_by_ctx (vb, STRa(XA_lookback_snip), rname_ctx, 0); // add_bytes=0 bc we already added them when we segged rname the first time

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

    seg_add_to_local_resizable (vb, CTX(OPTION_XA_LOOKBACK), lookback, 0);

    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_RNAME, false, (ValueType){.i = rname_index }, true);
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_POS,   false, (ValueType){.i = pos }, false);
    
    CTX(OPTION_XA_Z)->lookback_value = lookback; // for use when segging STRAND

    COPY_TIMER (sam_seg_XA_pos);
}

// Seg lookback callback for STRAND item of XA
static void sam_seg_XA_strand (VBlockP vb, WordIndex strand_index, unsigned add_bytes)
{
    ContextP strand_ctx = CTX (OPTION_XA_STRAND);
    int64_t lookback = CTX(OPTION_XA_Z)->lookback_value; // calculated in sam_seg_XA_pos

    if (lookback && lookback_get_index (vb, CTX(OPTION_XA_LOOKBACK), strand_ctx, lookback) == strand_index) 
        seg_by_ctx (vb, STRa(XA_lookback_snip), strand_ctx, add_bytes);
    else
        seg_known_node_index (vb, strand_ctx, strand_index, add_bytes); 

    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_STRAND, false, (ValueType){ .i = strand_index }, true);
}

// split the pos strand-pos string, eg "-10000" to strand "-" and pos "10000"
static bool seg_XA_strand_pos_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep)
{
    char strand = field[0];
    
    if (field_len < 2 || (strand != '+' && strand != '-'))  
        return false; // invalid XA format - expecting pos to begin with the strand

    if (segconf.is_sorted && !segconf.running) {
        sam_seg_XA_pos (vb, &field[1], field_len-1, rep);
        sam_seg_XA_strand (vb, strand == '+', 1); // index is 0 if '-' and 1 if '+' (set in sam_seg_0X_initialize)
    }
    // case: for a collated (or otherwise unsorted) file, we just seg normally (also: in segconf.running)
    else {
        seg_by_did (VB, field, 1, OPTION_XA_STRAND, 1);
        seg_integer_or_not (vb, CTX(OPTION_XA_POS), &field[1], field_len-1, field_len-1);
    }
    seg_by_ctx (VB, xa_strand_pos_snip, xa_strand_pos_snip_len, ctx, 0); // pre-created constant container

    return true; // segged successfully
}

// split the pos strand-XB-pos string, eg "-C10000" to strand "-" XB "C" and pos "10000"
static bool seg_seg_gem3_XA_strand_XB_pos_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep)
{
    char strand = field[0];
    char xb     = field[1];

    if (field_len < 3 || (strand != '+' && strand != '-') || (xb != 'C' && xb != 'G'))  
        return false; // invalid XA format - expecting pos to begin with the strand and XB

    if (segconf.is_sorted && !segconf.running) {
        sam_seg_XA_pos (vb, &field[2], field_len-2, rep);
        sam_seg_XA_strand (vb, 2 + ((strand=='+') | (xb=='G') << 1), 2); // word_indices set in sam_seg_BWA_XA_initialize
    }
    
    // case: for a collated (or otherwise unsorted) file, we just seg normally (also: in segconf.running)
    else {
        seg_by_did (VB, field, 2, OPTION_XA_STRAND, 2);
        seg_integer_or_not (vb, CTX(OPTION_XA_POS), &field[2], field_len-2, field_len-2);
    }
    seg_by_ctx (VB, xa_strand_pos_snip, xa_strand_pos_snip_len, ctx, 0); // pre-created constant container

    return true; // segged successfully
}

static bool seg_XA_lookup_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep)
{
    return true; // "segged successfully" - do nothing, we will seg it in seg_XA_strand_pos_cb
}

void sam_seg_BWA_XA_Z (VBlockSAMP vb, STRp(xa), unsigned add_bytes)
{
    static const MediumContainer container_XA = {
        .repeats      = 0, 
        .nitems_lo    = 5, 
        .filter_items = true, 
        .repsep       = {';'}, // including last item
        .items        = { { .dict_id = { _OPTION_XA_LOOKBACK   }, .separator = { CI0_INVISIBLE } }, 
                          { .dict_id = { _OPTION_XA_RNAME      }, .separator = {','} }, 
                          { .dict_id = { _OPTION_XA_STRAND_POS }, .separator = {','} },
                          { .dict_id = { _OPTION_XA_CIGAR      }, .separator = {','} }, // we don't mix the prirmary CIGAR field as the primary has a SNIP_SPECIAL
                          { .dict_id = { _OPTION_XA_NM         },                    } }  };

    SegCallback callbacks[5] = { [0]=seg_XA_lookup_cb, [1]=chrom_seg_cb, [3]=sam_seg_0A_cigar_cb,
                                 [2]=(MP(GEM3) && segconf.sam_bisulfite) ? seg_seg_gem3_XA_strand_XB_pos_cb : seg_XA_strand_pos_cb };
    int32_t repeats = seg_array_of_struct (VB, CTX(OPTION_XA_Z), container_XA, STRa(xa), callbacks, add_bytes);

    // case: we failed to seg as a container - flush lookbacks (rare condition, and complicated to rollback given the round-robin and unlimited repeats)
    if (repeats == -1) {
        lookback_flush (VB, CTX(OPTION_XA_RNAME));
        lookback_flush (VB, CTX(OPTION_XA_POS));
        lookback_flush (VB, CTX(OPTION_XA_STRAND));
    }
}

// ---------------
// PIZ
// ---------------

// this is called for XA that are a container, but not for invalid XA that are segged as a simple snip
void sam_piz_XA_field_insert_lookback (VBlockP vb)
{
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_RNAME,  true, (ValueType){ .i = 0 }, true);
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_STRAND, true, (ValueType){ .i = 0 }, true);
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_POS,    true, (ValueType){ .i = 0 }, false);
}
