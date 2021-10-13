// ------------------------------------------------------------------
//   sam_fields.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "sam_private.h"
#include "seg.h"
#include "codec.h"
#include "piz.h"
#include "reconstruct.h"
#include "context.h"
#include "container.h"
#include "chrom.h"
#include "optimize.h"
#include "strings.h"
#include "segconf.h"
#include "profiler.h"
#include "lookback.h"

static const StoreType optional_field_store_flag[256] = {
    ['c']=STORE_INT, ['C']=STORE_INT, 
    ['s']=STORE_INT, ['S']=STORE_INT,
    ['i']=STORE_INT, ['I']=STORE_INT,
    ['f']=STORE_FLOAT
};

const char optional_sep_by_type[2][256] = { { // compressing from SAM
        ['c']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['C']=CI_NATIVE_NEXT | CI_TRANS_NOR, // reconstruct number and \t separator is SAM, and don't reconstruct anything if BAM (reconstruction will be done by translator)
        ['s']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['S']=CI_NATIVE_NEXT | CI_TRANS_NOR, // -"-
        ['i']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['I']=CI_NATIVE_NEXT | CI_TRANS_NOR, // -"-
        ['f']=CI_NATIVE_NEXT | CI_TRANS_NOR,                                      // compressing SAM - a float is stored as text, and when piz with translate to BAM - is not reconstructed, instead - translated
        ['Z']=CI_NATIVE_NEXT | CI_TRANS_NUL, ['H']=CI_NATIVE_NEXT | CI_TRANS_NUL, // reconstruct text and then \t separator if SAM and \0 if BAM 
        ['A']=CI_NATIVE_NEXT,                                                     // reconstruct character and then \t separator if SAM and no separator for BAM
        ['B']=CI_NATIVE_NEXT                                                      // reconstruct array and then \t separator if SAM and no separator for BAM
}, 
{ // compressing from BAM
        ['c']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['C']=CI_NATIVE_NEXT | CI_TRANS_NOR, // reconstruct number and \t separator is SAM, and don't reconstruct anything if BAM (reconstruction will be done by translator)
        ['s']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['S']=CI_NATIVE_NEXT | CI_TRANS_NOR, // -"-
        ['i']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['I']=CI_NATIVE_NEXT | CI_TRANS_NOR, // -"-
        ['f']=CI_NATIVE_NEXT,                                                     // compressing SAM - a float is stored as a SPECIAL, and the special reconstructor handles the SAM and BAM reconstructing
        ['Z']=CI_NATIVE_NEXT | CI_TRANS_NUL, ['H']=CI_NATIVE_NEXT | CI_TRANS_NUL, // reconstruct text and then \t separator if SAM and \0 if BAM 
        ['A']=CI_NATIVE_NEXT,                                                     // reconstruct character and then \t separator if SAM and no separator for BAM
        ['B']=CI_NATIVE_NEXT                                                      // reconstruct array and then \t separator if SAM and no separator for BAM
} };

//--------------
// FLAG
//--------------

void sam_seg_FLAG (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(flag_str)/* optional, otherwise in dl */, unsigned add_bytes)
{
    if (flag_str_len) 
        ASSSEG (str_get_int_range16 (STRa(flag_str), 0, SAM_MAX_FLAG, &dl->FLAG.value), flag_str, "invalid FLAG field: %.*s", STRf(flag_str));

    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // note: an invalid pointer if buddy_line_i is -1
    
    // case: we can retrieve the FLAG from this line's buddy
    #define SAME_AS_BUDDY_FLAGS (SAM_FLAG_MULTI_SEGMENTS | SAM_FLAG_IS_ALIGNED | SAM_FLAG_SECONDARY | SAM_FLAG_FAILED_FILTERS | \
                                 SAM_FLAG_DUPLICATE | SAM_FLAG_SUPPLEMENTARY)
    if (segconf.sam_is_sorted && vb->buddy_line_i != -1 &&
        (dl->FLAG.value & SAME_AS_BUDDY_FLAGS) == (buddy_dl->FLAG.value & SAME_AS_BUDDY_FLAGS) &&
        dl->FLAG.bits.unmapped       == buddy_dl->FLAG.bits.next_unmapped &&
        dl->FLAG.bits.next_unmapped  == buddy_dl->FLAG.bits.unmapped &&
        dl->FLAG.bits.rev_comp       == buddy_dl->FLAG.bits.next_rev_comp &&
        dl->FLAG.bits.next_rev_comp  == buddy_dl->FLAG.bits.rev_comp &&
        dl->FLAG.bits.is_first       == buddy_dl->FLAG.bits.is_last &&
        dl->FLAG.bits.is_last        == buddy_dl->FLAG.bits.is_first)
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_COPY_BUDDY_FLAG }, 2, SAM_FLAG, add_bytes); // added 12.0.41
    
    // case: normal snip
    else {
        seg_integer (vb, SAM_FLAG, dl->FLAG.value, false);
        CTX(SAM_FLAG)->txt_len += add_bytes; 
    }
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_COPY_BUDDY_FLAG)
{
    ASSPIZ0 (vb->buddy_line_i >= 0, "No buddy line is set for the current line");

    SamFlags flag = { .value = *ENT (int64_t, ctx->history, vb->buddy_line_i) }; 
    SWAPbit (flag.bits.unmapped, flag.bits.next_unmapped);
    SWAPbit (flag.bits.rev_comp, flag.bits.next_rev_comp);
    SWAPbit (flag.bits.is_first, flag.bits.is_last);
    
    new_value->i = flag.value;
    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return true; // new value
}

//--------------
// TLEN
//--------------

// TLEN - 4 cases: 
// 1. case: a non-zero value that is the negative of its mate in a sorted file - a COPY_BUDDY_TLEN special
// 2. case: a non-zero value that is the negative of the previous line (usually a mate in a collated file) - a SNIP_DELTA & "-" (= value negation)
// 3. case: tlen>0 and pnext_pos_delta>0 and seq_len>0 tlen is stored as SNIP_SPECIAL & tlen-pnext_pos_delta-seq_len
// 4. otherwise: stored as is
void sam_seg_TLEN (VBlockSAM *vb, ZipDataLineSAM *dl, 
                         const char *tlen, unsigned tlen_len, // option 1
                         int64_t tlen_value, // option 2
                         PosType pnext_pos_delta, int32_t cigar_seq_len)
{
    Context *ctx = CTX(SAM_TLEN);

    if (tlen) { // option 1
        ASSSEG0 (tlen_len, tlen, "empty TLEN");

        bool is_int = str_get_int (tlen, tlen_len, &tlen_value); // note: tlen_value remains 0 if not a valid integer
        ASSSEG (is_int, tlen, "expecting TLEN to be an integer, but found \"%.*s\"", tlen_len, tlen);
    }

    unsigned add_bytes = IS_BAM ? sizeof (uint32_t) : tlen_len + 1;

    // case 1: tlen is minus its buddy in a sorted file - mate is likely close but not necessarily adjacent
    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // note: an invalid pointer if buddy_line_i is -1
    if (segconf.sam_is_sorted && tlen_value && vb->buddy_line_i != -1 && buddy_dl->TLEN == -tlen_value) 
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_COPY_BUDDY_TLEN }, 2, ctx, add_bytes); // added 12.0.41

    // case 2: tlen is minus previous tlen - usually adjacent mates in a collated file
    else if (tlen_value && tlen_value == -ctx->last_value.i) 
        seg_by_ctx (VB, (char[]){ SNIP_SELF_DELTA, '-' }, 2, ctx, add_bytes);
    
    // case 3:
    else if (tlen_value > 0 && pnext_pos_delta > 0 && cigar_seq_len > 0) {
        char tlen_by_calc[30] = { SNIP_SPECIAL, SAM_SPECIAL_TLEN };
        unsigned tlen_by_calc_len = str_int (tlen_value - pnext_pos_delta - (int64_t)cigar_seq_len, &tlen_by_calc[2]);
        seg_by_ctx (VB, tlen_by_calc, tlen_by_calc_len + 2, ctx, add_bytes);
    }
    
    // case 4: tlen=ref_consumed (happens eg in pacbio)
    else if (tlen_value == vb->ref_consumed) 
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED }, 2, ctx, add_bytes);
    
    // case default: add as is (option 1)
    else if (tlen)
        seg_by_ctx (VB, tlen, tlen_len, ctx, add_bytes);

    // case default: add as is (option 2)
    else {
        char snip[20];
        unsigned snip_len = str_int (tlen_value, snip);
        seg_by_ctx (VB, snip, snip_len, ctx, add_bytes);
    }

    ctx->last_value.i = tlen_value;
    dl->TLEN = tlen_value;
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_TLEN)
{
    ASSERT0 (snip_len, "snip_len=0");

    int32_t tlen_by_calc = atoi (snip);
    int32_t tlen_val = tlen_by_calc + CTX(SAM_PNEXT)->last_delta + vb->seq_len;

    new_value->i = tlen_val;

    if (reconstruct) RECONSTRUCT_INT (tlen_val);

    return true; // new value
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_COPY_BUDDY_TLEN)
{
    ASSPIZ0 (vb->buddy_line_i >= 0, "No buddy line is set for the current line");

    new_value->i = -*ENT (int64_t, ctx->history, vb->buddy_line_i); // minus the buddy
    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return true; // new value
}

// -------------------------------------------------------------------------------------------------------------------------------------------
// U2:Z "Phred probability of the 2nd call being wrong conditional on the best being wrong" (https://samtools.github.io/hts-specs/SAMtags.pdf)
// -------------------------------------------------------------------------------------------------------------------------------------------

// callback function for compress to get data of one line
void sam_zip_U2 (VBlock *vb, uint64_t vb_line_i, char **line_u2_data,  uint32_t *line_u2_len, uint32_t maximum_len) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    *line_u2_len = MIN_(maximum_len, dl->U2.snip_len);

    if (!line_u2_data) return; // only lengths were requested

    *line_u2_data = ENT (char, vb->txt_data, dl->U2.char_index);

    if (flag.optimize_QUAL)
        optimize_phred_quality_string (*line_u2_data, *line_u2_len);
}

// ---------
// BD and BI
// ---------

static void sam_seg_BD_BI_field (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(field), DictId dict_id, unsigned add_bytes)
{
    bool is_bi = (dict_id.num == _OPTION_BI_Z);
    Context *this_ctx  = is_bi ? CTX(OPTION_BI_Z) : CTX (OPTION_BD_Z);

    if (field_len != dl->seq_len) {
        seg_by_ctx (VB, field, field_len, this_ctx, field_len);
        return;
    }
    
    dl->BD_BI[is_bi] = WORD_IN_TXT_DATA (field);

    CTX(OPTION_BD_BI)->txt_len += add_bytes; 

    if (!dl->BD_BI[!is_bi].char_index) // the first of BD and BI increments local.len, so it is incremented even if just one of BD/BI appears
        CTX(OPTION_BD_BI)->local.len += field_len * 2;

    seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, SAM_SPECIAL_BDBI }), 2, this_ctx, 0);
}

// callback function for compress to get BD_BI data of one line: this is an
// interlaced line containing a character from BD followed by a character from BI - since these two fields are correlated
void sam_zip_BD_BI (VBlock *vb_, uint64_t vb_line_i, 
                    char **line_data, uint32_t *line_len,  // out 
                    uint32_t maximum_len)
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    
    const char *bd = dl->BD_BI[0].char_index ? ENT (char, vb->txt_data, dl->BD_BI[0].char_index) : NULL;
    const char *bi = dl->BD_BI[1].char_index ? ENT (char, vb->txt_data, dl->BD_BI[1].char_index) : NULL;
    
    if (!bd && !bi) return; // no BD or BI on this line

    ASSERT (bd && bi, "A line has one of the BD:Z/BI:Z pair - Genozip can only compress lines that have either both BD:Z and BI:Z or neither. vb=%u vb_line_i=%"PRIu64, 
            vb->vblock_i, vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_len  = MIN_(maximum_len, dl->seq_len * 2);

    if (!line_data) return; // only length was requested

    buf_alloc (vb, &vb->bd_bi_line, 0, dl->seq_len * 2, uint8_t, 2, "bd_bi_line");

    // calculate character-wise delta
    for (unsigned i=0; i < dl->seq_len; i++) {
        *ENT (uint8_t, vb->bd_bi_line, i*2    ) = bd[i];
        *ENT (uint8_t, vb->bd_bi_line, i*2 + 1) = bi[i] - (bd ? bd[i] : 0);
    }

    *line_data = FIRSTENT (char, vb->bd_bi_line);
}   

// BD and BI - reconstruct from BD_BI context which contains interlaced BD and BI data. 
SPECIAL_RECONSTRUCTOR (sam_piz_special_BD_BI)
{
    if (!vb->seq_len || !reconstruct) goto done;

    Context *bdbi_ctx = CTX(OPTION_BD_BI);

    // note: bd and bi use their own next_local to retrieve data from bdbi_ctx. the actual index
    // in bdbi_ctx.local is calculated given the interlacing
    ASSERT (ctx->next_local + vb->seq_len * 2 <= bdbi_ctx->local.len, "Error reading txt_line=%"PRIu64": unexpected end of %s data", vb->line_i, dis_dict_id (ctx->dict_id).s);

    char *dst        = AFTERENT (char, vb->txt_data);
    const char *src  = ENT (char, bdbi_ctx->local, ctx->next_local * 2);
    uint32_t seq_len = vb->seq_len; // automatic var for effeciency

    if (ctx->dict_id.num == _OPTION_BD_Z)
        for (uint32_t i=0; i < seq_len; i++, src+=2, dst++) *dst = *src;
    else
        for (uint32_t i=0; i < seq_len; i++, src+=2, dst++) *dst = *src + *(src+1);
    
    vb->txt_data.len += vb->seq_len;    
    ctx->next_local  += vb->seq_len;

done:
    return false; // no new value
}

// -------------------------------------------------------
// XA:Z "Alternative alignments" (a BWA-backtrack feature)
// -------------------------------------------------------

static bool sam_seg_0A_cigar_cb (VBlockP vb, ContextP ctx, STRp (cigar), uint32_t repeat)
{
    // complicated CIGARs, likely to be relatively uncommon, are better off in local - anything more than eg 112M39S 
    // note: we set no_stons=true in sam_seg_initialize_0X so we can use local for this rather than singletons
    if (cigar_len > 7 && !flag.best) {
        seg_add_to_local_text (vb, ctx, STRa(cigar), cigar_len);
        seg_simple_lookup (vb, ctx, 0);
    }
    else 
        seg_by_ctx (vb, STRa(cigar), ctx, cigar_len);

    return true;
}

// Lookup buffer
#define lookback_buf ctx_specific_buf // we store the previous rname, pos, strand in their ctx->ctx_specific_buf buffer
#define lookback_value last_value.i
#define MAX_LOOKUPS (1 << XA_LOOKBACK_DEPTH_BITS)

// Seg lookback callback for POS item of XA - seg the POS relative to lookback pos, and replace the already segged RNAME with a lookback if there is one.
static void sam_seg_XA_pos (VBlockP vb, STRp(pos_str), uint32_t rep)
{
    #define MAX_POS_DISTANCE 10000 // the smaller it is, the more we search for a better XA - the delta-pos will be better, but lookback worse. 10K works well.
    START_TIMER;

    ContextP rname_ctx = CTX (OPTION_XA_RNAME);
    ContextP pos_ctx   = CTX (OPTION_XA_POS);

    // look back for a node with this index and a similar POS - we use word_index to store the original rname_node_index, pos
    WordIndex rname_index = *LASTENT (WordIndex, rname_ctx->b250);
    int64_t lookback = 0;
    PosType pos = -MAX_POS_DISTANCE; // initial to "invalid pos" - value chosen so we can storeit in poses, in lieu of an invalid non-integer value, without a future pos being considered close to it

    if (!segconf.running && str_get_int (STRa(pos_str), &pos)) {

        int64_t iterator = -1;
        while ((lookback = lookback_get_next (vb, rname_ctx, rname_index, &iterator))) {

            PosType lookback_pos = lookback_get_value (vb, pos_ctx, lookback);
//printf ("xxx lookback_pos=%"PRId64" lookback=%u iterator=%"PRId64"\n", lookback_pos, lookback, iterator);
            // case: we found a lookback - same rname and close enough pos
            if (ABS (pos-lookback_pos) < MAX_POS_DISTANCE) {
            //    if (ABS (pos-lookback_pos) <= ABS(CTX(SAM_TLEN)->last_value.i)) { <-- better POS deltas but bigger index - its a wash

                // replace rname with lookback
                rname_ctx->b250.len--;
                ctx_decrement_count (vb, rname_ctx, rname_index);
                
                seg_by_ctx (vb, STRa(xa_lookback_snip), rname_ctx, 0); // add_bytes=0 bc we already added them when we segged rname the first time

                // seg pos as a delta
                SNIP(48);
                memcpy (snip, xa_lookback_snip, xa_lookback_snip_len);
                snip_len = xa_lookback_snip_len + str_int (pos - lookback_pos, &snip[xa_lookback_snip_len]);
                seg_by_ctx (vb, STRa(snip), pos_ctx, pos_str_len);
                pos_ctx->numeric_only = false; // not all numeric (cancel the setting by seg_integer_or_not)
//if (!segconf.running) printf ("xxx rep=%u Lookback: \"%.*s\" pos=%u prev_pos=%u delta=%d\n", rep, STRf(pos_str), pos, lookback_pos, pos - lookback_pos);     

                break;
            }
        }
    }

    if (!lookback) {
//if (!segconf.running)printf ("xxx rep=%u No lookback: %.*s\n", rep, STRf(pos_str));     
        seg_integer_or_not (vb, pos_ctx, STRa(pos_str), pos_str_len);
    }

    seg_add_to_local_uint32 (vb, CTX(OPTION_XA_LOOKBACK), lookback, 0);

    lookback_insert (vb, OPTION_XA_RNAME, rname_index, true);
    lookback_insert (vb, OPTION_XA_POS, pos, false);
    
    CTX(OPTION_XA_Z)->lookback_value = lookback; // for use when segging STRAND

    COPY_TIMER (sam_seg_XA_pos);
}

// Seg lookback callback for STRAND item of XA
static void sam_seg_XA_strand (VBlockP vb, WordIndex strand_index)
{
    ContextP strand_ctx = CTX (OPTION_XA_STRAND);
    int64_t lookback = CTX(OPTION_XA_Z)->lookback_value; // calculated in sam_seg_XA_pos

    if (lookback && lookback_get_index (vb, strand_ctx, lookback) == strand_index) 
        seg_by_ctx (vb, STRa(xa_lookback_snip), strand_ctx, 1);
    else
        seg_known_node_index (vb, strand_ctx, strand_index, 1); 

    lookback_insert (vb, OPTION_XA_STRAND, strand_index, true);
}

// split the pos strand-pos string, eg "-10000" to strand "-" and pos "10000"
static bool seg_XA_strand_pos_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep)
{
    if (field_len < 2 || (field[0] != '+' && field[0] != '-'))  
        return false; // invalid XA format - expecting pos to begin with the strand

    if (segconf.sam_is_sorted && !segconf.running) {
        sam_seg_XA_pos (vb, &field[1], field_len-1, rep);
        sam_seg_XA_strand (vb, *field == '+'); // index is 0 if '-' and 1 if '+' (set in sam_seg_initialize_0X)
    }
    // case: for a collated (or otherwise unsorted) file, we just seg normally (also: in segconf.running)
    else {
        seg_by_did_i (VB, field, 1, OPTION_XA_STRAND, 1);
        seg_integer_or_not (vb, CTX(OPTION_XA_POS), &field[1], field_len-1, field_len-1);
    }
    seg_by_ctx (VB, xa_strand_pos_snip, xa_strand_pos_snip_len, ctx, 0); // pre-created constant container

    return true; // segged successfully
}

static bool seg_XA_lookup_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep)
{
    return true; // "segged successfully" - do nothing, we will seg it in seg_XA_strand_pos_cb
}

// XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
// Example XA:Z:chr9,-60942781,150M,0;chr9,-42212061,150M,0;chr9,-61218415,150M,0;chr9,+66963977,150M,1;
// See: http://bio-bwa.sourceforge.net/bwa.shtml
// Note that a different XA was observed in IonXpress data: "XA:Z:map4-1". This will be just segged as fallback
static void sam_seg_BWA_XA_field (VBlockSAM *vb, STRp(xa))
{
    static const MediumContainer container_XA = {
        .repeats      = 0, 
        .nitems_lo    = 5, 
        .filter_items = true, 
        .repsep       = {';'}, // including last item
        .items        = { { .dict_id = { _OPTION_XA_LOOKBACK   }, .separator = { CI_INVISIBLE } }, 
                          { .dict_id = { _OPTION_XA_RNAME      }, .separator = {','} }, 
                          { .dict_id = { _OPTION_XA_STRAND_POS }, .separator = {','} },
                          { .dict_id = { _OPTION_XA_CIGAR      }, .separator = {','} }, // we don't mix the prirmary CIGAR field as the primary has a SNIP_SPECIAL
                          { .dict_id = { _OPTION_XA_NM         },                    } }  };

    static const SegCallback callbacks[5] = { [0]=seg_XA_lookup_cb, [1]=chrom_seg_cb, [2]=seg_XA_strand_pos_cb, [3]=sam_seg_0A_cigar_cb };
    int32_t repeats = seg_array_of_struct (VB, CTX(OPTION_XA_Z), container_XA, STRa(xa), callbacks);
    CTX(OPTION_XA_Z)->txt_len++; // +1 for '\t' in SAM or '\0' in BAM 

    // case: we failed to seg as a container - flush lookbacks (rare condition, and complicated to rollback given the round-robin and unlimited repeats)
    if (repeats == -1) {
        lookback_flush (VB, CTX(OPTION_XA_RNAME));
        lookback_flush (VB, CTX(OPTION_XA_POS));
        lookback_flush (VB, CTX(OPTION_XA_STRAND));
    }
    else if (segconf.running) 
        segconf.XA_reps += repeats;
}

static int sam_seg_which_XA (STRp(xa))
{
    unsigned semicolons;
    if ((semicolons = str_count_char (STRa(xa), ';')) >= 1 && str_count_char (STRa(xa), ',') == semicolons * 3) 
        return XA_BWA;

    if (memchr (xa, '-', xa_len))
        return XA_IONTORRENT;

    else
        return XA_UNKNOWN;
}

// this is called for XA that are a container, but not for invalid XA that are segged as a simple snip
void sam_piz_XA_field_insert_lookback (VBlockP vb)
{
    lookback_insert (vb, OPTION_XA_RNAME,  TAKE_LAST_VALUE, true);
    lookback_insert (vb, OPTION_XA_STRAND, TAKE_LAST_VALUE, true);
    lookback_insert (vb, OPTION_XA_POS,    TAKE_LAST_VALUE, false);
}

// ---------------------------------------------------------
// SA:Z "Other canonical alignments in a chimeric alignment"
// ---------------------------------------------------------

// OA and SA format is: (rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ . in OA - NM is optional (but its , is not)
// Example SA:Z:chr13,52863337,-,56S25M70S,0,0;chr6,145915118,+,97S24M30S,0,0;chr18,64524943,-,13S22M116S,0,0;chr7,56198174,-,20M131S,0,0;chr7,87594501,+,34S20M97S,0,0;chr4,12193416,+,58S19M74S,0,0;
// See: https://samtools.github.io/hts-specs/SAMtags.pdf
// note: even though SA, OA, XA contain similar fields amongst each other and similar to the primary fields,
// the values of subsequent lines tend to be similar for each one of them seperately, so we maintain separate contexts

static void sam_seg_SA_field (VBlockSAM *vb, STRp(field))
{
    static const MediumContainer container_SA = { .nitems_lo = 6,      
                                                  .repsep    = { ';' }, // including on last repeat    
                                                  .items     = { { .dict_id = { _OPTION_SA_RNAME  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_POS    }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_STRAND }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_CIGAR  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_MAPQ   }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_NM     },                  } } };

    SegCallback callbacks[6] = { [0]=chrom_seg_cb, [1]=seg_pos_field_cb, [3]=sam_seg_0A_cigar_cb };
     
    int32_t repeats = seg_array_of_struct (VB, CTX(OPTION_SA_Z), container_SA, field, field_len, callbacks);

    CTX(OPTION_SA_Z)->txt_len++; // 1 for \t in SAM and \0 in BAM 

    if (segconf.running && repeats > 0) segconf.SA_reps += repeats;
}

// -------------------------
// OA:Z "Original alignment"
// -------------------------

static void sam_seg_OA_field (VBlockSAM *vb, STRp(field))
{
    static const MediumContainer container_OA = { .nitems_lo = 6,          
                                                  .repsep    = { ';' }, // including on last repeat    
                                                  .items     = { { .dict_id = { _OPTION_OA_RNAME  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_POS    }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_STRAND }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_CIGAR  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_MAPQ   }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_NM     },                    } } };

    SegCallback callbacks[6] = { [0]=chrom_seg_cb, [1]=seg_pos_field_cb, [3]=sam_seg_0A_cigar_cb };
     
    int32_t repeats = seg_array_of_struct (VB, CTX(OPTION_OA_Z), container_OA, field, field_len, callbacks);

    CTX(OPTION_OA_Z)->txt_len++; // 1 for \t in SAM and \0 in BAM 

    if (segconf.running && repeats > 0) segconf.SA_reps += repeats;
}

// ----------------------------
// NM:i "Number of differences"
// ----------------------------

// Two variations:
// 1) Integer NM per SAM specification https://samtools.github.io/hts-specs/SAMtags.pdf: "Number of differences (mismatches plus inserted and deleted bases) 
// between the sequence and reference, counting only (case-insensitive) A, C, G and T bases in sequence and reference as potential matches, with everything
// else being a mismatch. Note this means that ambiguity codes in both sequence and reference that match each other, such as ‘N’ in both, or compatible 
// codes such as ‘A’ and ‘R’, are still counted as mismatches. The special sequence base ‘=’ will always be considered to be a match, even if the reference 
// is ambiguous at that point. Alignment reference skips, padding, soft and hard clipping (‘N’, ‘P’, ‘S’ and ‘H’ CIGAR operations) do not count as mismatches,
// but insertions and deletions count as one mismatch per base."
// 2) Binary NM: 0 if sequence fully matches the reference when aligning according to CIGAR, 1 is not.
static void sam_seg_NM_field (VBlockSAM *vb, STRp(field), unsigned add_bytes)
{
    int32_t NM;
    if (!str_get_int_range32 (STRa(field), 0, 10000000, &NM)) goto fallback;
    
    if (segconf.running && NM > 1) segconf.NM_is_integer = true; // we found evidence of integer NM

    ctx_set_last_value (VB, CTX (OPTION_NM_i), ((LastValueType){ .i = NM }));

    if (segconf.NM_is_integer && NM == vb->mismatch_bases)
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'i'}, 3, OPTION_NM_i, add_bytes); 

    else if (!segconf.NM_is_integer && (NM > 0) == (vb->mismatch_bases > 0))
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'b'}, 3, OPTION_NM_i, add_bytes); 

    else
        fallback:
        seg_by_did_i (VB, field, field_len, OPTION_NM_i, add_bytes); 
}

SPECIAL_RECONSTRUCTOR (bam_piz_special_NM)
{
    if (*snip == 'i') 
        new_value->i = ((VBlockSAMP)vb)->mismatch_bases;

    else if (*snip == 'b')
        new_value->i = ((VBlockSAMP)vb)->mismatch_bases > 0;

    else 
        ASSPIZ (false, "unrecognized opcode '%c'", *snip);

    if (reconstruct) // will be false if BAM, reconstruction is done by translator based on new_value set here
        RECONSTRUCT_INT (new_value->i);

    return true; // has new value
}

// -------------------------------------------------------------------------------------------------------------------
// XM:i Case 1: BWA: "Number of mismatches in the alignment" 
//      Case 2: IonTorrent TMAP: "The target length, that is, the number of reference bases spanned by the alignment." 
// -------------------------------------------------------------------------------------------------------------------

static void sam_seg_XM_field (VBlockSAM *vb, STRp(xm_str), unsigned add_bytes)
{
    int64_t xm;
    ContextP NM_ctx;
    if (!str_get_int (STRa(xm_str), &xm))
        goto fallback;
    
    // check for BWA case - XM is similar to NM - in our test file > 99% identical to NM.
    if (ctx_has_value_in_line (vb, _OPTION_NM_i, &NM_ctx) && xm == NM_ctx->last_value.i) 
        seg_by_did_i (VB, STRa(XM_snip), OPTION_XM_i, add_bytes); // copy from NM
    
    // check IonTorrent TMAP case - XM is supposed to be ref_consumed
    else if (xm == vb->ref_consumed)                              
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED }, 2, OPTION_XM_i, add_bytes);
        
    else
        fallback:
        seg_by_did_i (VB, STRa(xm_str), OPTION_XM_i, add_bytes);  // normal seg
}

// ----------------------------------------------------------------------------------------------
// AS:i "Alignment score generated by aligner" (https://samtools.github.io/hts-specs/SAMtags.pdf)
// ----------------------------------------------------------------------------------------------

// AS has a value set (at least as set by BWA and IonTorrent TMAP) of at most vb->ref_consumed, and often equal to it. we modify
// it to be new_value=(value-ref_consumed) 
static inline void sam_seg_AS_field (VBlockSAM *vb, STRp(as_str), unsigned add_bytes)
{
    int32_t as;    
    if (str_get_int_range32 (STRa(as_str), 0, 100000, &as)) {
        ctx_set_last_value (VB, CTX (OPTION_AS_i), ((LastValueType){ .i = as }));

        // store a special snip with delta
        char new_snip[20] = { SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED };
        unsigned delta_len = str_int ((int32_t)vb->ref_consumed-as, &new_snip[2]);

        seg_by_ctx (VB, new_snip, delta_len+2, CTX (OPTION_AS_i), add_bytes); 
    }

    // not possible - just store unmodified
    else
        seg_by_ctx (VB, STRa(as_str), CTX (OPTION_AS_i), add_bytes); 
}

// reconstruct seq_len or (seq_len-snip)
// Note: This is used by AS:i fields in files compressed up to 12.0.37
SPECIAL_RECONSTRUCTOR (sam_piz_special_SEQ_LEN)
{
    new_value->i = (int32_t)vb->seq_len - atoi (snip); // seq_len if snip=""
    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return true; // new value
}

// reconstruct ref_consumed or (ref_consumed-snip)
SPECIAL_RECONSTRUCTOR (sam_piz_special_REF_CONSUMED)
{
    new_value->i = (int32_t)((VBlockSAMP)vb)->ref_consumed - atoi (snip); // ref_consumed if snip=""
    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return true; // new value
}

// ----------------------------------------------------------------------------------------------
// XS:i 
// ----------------------------------------------------------------------------------------------

// often its the same as AS:i
static inline void sam_seg_XS_field (VBlockSAM *vb, STRp(xs_str), unsigned add_bytes)
{
    int32_t xs = 0;
    if (ctx_has_value_in_line_(VB, CTX(OPTION_AS_i)) && str_get_int_range32 (STRa(xs_str), 0, 10000, &xs) && 
        xs == CTX(OPTION_AS_i)->last_value.i) 
        seg_by_did_i (VB, STRa(XS_snip), OPTION_XS_i, add_bytes);

    else if (xs) {
        // store a special snip with delta
        char new_snip[20] = { SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED };
        unsigned delta_len = str_int (vb->ref_consumed-xs, &new_snip[2]);

        seg_by_did_i (VB, new_snip, delta_len+2, OPTION_XS_i, add_bytes); 
    }
    else
        seg_by_did_i (VB, STRa(xs_str),  OPTION_XS_i, add_bytes); 
}

// MQ:i is often very similar to MAPQ
static inline void sam_seg_MQ_field (VBlockSAM *vb, ZipDataLineSAM *dl, 
                                     const char *snip, unsigned snip_len, unsigned add_bytes)
{
    int32_t mq, mapq;
    if (str_get_int_range32 (STRa(snip), 0, 1000, &mq) &&
        str_get_int_range32 (last_txt(vb, SAM_MAPQ), vb->last_txt_len(SAM_MAPQ), 0, 1000, &mapq)) {
        
        int32_t delta = mq - mapq;

        SNIP(100);
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _SAM_MAPQ, true, delta, snip);
        seg_by_ctx (VB, STRa(snip), CTX(OPTION_MQ_i), add_bytes);
    }
    else
        seg_by_ctx (VB, snip, snip_len, CTX (OPTION_MQ_i), add_bytes); 
}

// mc:i: (output of bamsormadup and other biobambam tools - mc in small letters) 
// appears to be a pos value usually close to PNEXT, but it is -1 is POS=PNEXT.
// from bamsort manual: "adddupmarksupport=<0|1>: add information required for streaming duplicate marking in the aux fields MS and MC.
// Input is assumed to be collated by query name. This option is ignored unless fixmates=1. By default it is disabled."
// https://github.com/gt1/biobambam2/blob/master/src/programs/bamsort.cpp says: "biobambam used MC as a mate coordinate tag which now has a clash
// with the official SAM format spec.  New biobambam version uses mc."
// ms="MateBaseScore" - sum all characters in QUAL of the mate, where the value of each character is its ASCII minus 33 (i.e. the Phred score)
// mc="MateCoordinate"
static inline void sam_seg_mc_field (VBlockSAM *vb, DictId dict_id, 
                                     const char *snip, unsigned snip_len, unsigned add_bytes)
{
    uint8_t mc_did_i = ctx_get_ctx (vb, dict_id)->did_i;
    
    // if snip is "-1", store as simple snip
    if (snip_len == 2 && snip[0] == '-' && snip[1] == '1')
        seg_by_did_i (VB, snip, snip_len, mc_did_i, add_bytes);
    
    // delta vs PNEXT
    else
        seg_pos_field (VB, mc_did_i, SAM_PNEXT, SPF_BAD_SNIPS_TOO, 0, snip, snip_len, 0, add_bytes);
}

// possibly a special snip for copying RG from our buddy
static inline void sam_seg_RG_field (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(rg), unsigned add_bytes)
{
    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1

    if (segconf.sam_buddy_RG && vb->buddy_line_i != -1 && 
        buddy_dl->RG.snip_len == rg_len && !memcmp (rg, ENT (char, vb->txt_data, buddy_dl->RG.char_index), rg_len)) 

        seg_by_did_i (VB, (char[]){ SNIP_COPY_BUDDY }, 1, OPTION_RG_Z, add_bytes);
    else
        seg_by_did_i (VB, STRa(rg), OPTION_RG_Z, add_bytes);    

    dl->RG = WORD_IN_TXT_DATA(rg);
}

// optimization for Ion Torrent flow signal (ZM) - negative values become zero, positives are rounded to the nearest 10
static void sam_optimize_ZM (const char **snip, unsigned *snip_len, char *new_str)
{
    char *after;
    int number = strtoul (*snip, &after, 10);

    if ((unsigned)(after - *snip) > 0) {
        if (number >= 0) number = ((number + 5) / 10) * 10;
        else             number = 0;

        *snip_len = str_int (number, new_str);
        *snip = new_str;
    }    
}

// E2 - SEQ data. Currently broken. To do: fix.
/*static void sam_seg_E2_field (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(field), unsigned add_bytes)
{
    ASSSEG0 (dl->seq_len, field, "E2 tag without a SEQ"); 
    ASSINP (field_len == dl->seq_len, 
            "Error in %s: Expecting E2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
            txt_name, dl->seq_len, field_len, field_len, field);

    PosType this_pos = vb->last_int(SAM_POS);

    sam_seg_SEQ (vb, OPTION_E2_Z, (char *)field, field_len, this_pos, vb->last_cigar, vb->ref_consumed, vb->ref_and_seq_consumed, 0, field_len, // remove const bc SEQ data is actually going to be modified
                        vb->last_cigar, add_bytes); 
}*/

// U2 - QUAL data (note: U2 doesn't have a context - it shares with QUAL)
static void sam_seg_U2_field (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(field), unsigned add_bytes)
{
    ASSSEG0 (dl->seq_len, field, "U2 tag without a SEQ"); 
    ASSINP (field_len == dl->seq_len, 
            "Error in %s: Expecting U2 data to be of length %u as indicated by CIGAR, but it is %u. U2=%.*s",
            txt_name, dl->seq_len, field_len, field_len, field);

    dl->U2 = WORD_IN_TXT_DATA (field);
    CTX(OPTION_U2_Z)->txt_len   += add_bytes;
    CTX(OPTION_U2_Z)->local.len += field_len;
}

static inline unsigned sam_seg_optional_add_bytes (char type, unsigned value_len, bool is_bam)
{
    if (is_bam)
        switch (type) {
            case 'c': case 'C': case 'A': return 1;
            case 's': case 'S':           return 2;
            case 'i': case 'I': case 'f': return 4;
            case 'Z': case 'H':           return value_len + 1; // +1 for \0
            default : return 0;
        }
    else // SAM
        return value_len + 1; // +1 for \t
}

// an array - all elements go into a single item context, multiple repeats
static void sam_seg_array_field (VBlock *vb, DictId dict_id, const char *value, unsigned value_len)
{   
    // get optimization function, if there is one
    SegOptimize optimize = NULL;
    if (flag.optimize_ZM && dict_id.num == _OPTION_ZM_B && value_len > 3 && value[0] == 's')  // XM:B:s,
        optimize = sam_optimize_ZM;

    // prepare array container - a single item, with number of repeats of array element. array type is stored as a prefix
    Context *container_ctx     = ctx_get_ctx (vb, dict_id);

    SmallContainer con = { .nitems_lo = 2, 
                           .drop_final_item_sep_of_final_repeat = true, // TODO - get rid of this flag and move to making the seperators to be repeat seperators as they should have been, using drop_final_repeat_sep and obsoleting this flag 
                           .repsep    = {0,0}, 
                           .items     = { { .translator = SAM2BAM_ARRAY_SELF  },  // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM 
                                          { .separator  = {0, ','}            } } // item[1] is actual array item
                         };
    
    char prefixes[] = { CON_PREFIX_SEP, value[0], ',', CON_PREFIX_SEP }; // prefix contains type eg "i,"
    
    const char *str = value + 2;      // remove type and comma
    int str_len = (int)value_len - 2; // must be int, not unsigned, for the for loop

    // prepare context where array elements will go in
    char arr_dict_id_str[8]   = "XX_ARRAY";
    arr_dict_id_str[0]        = FLIP_CASE (dict_id.id[0]);
    arr_dict_id_str[1]        = FLIP_CASE (dict_id.id[1]);
    con.items[1].dict_id      = dict_id_make (arr_dict_id_str, 8, DTYPE_SAM_OPTIONAL);
    con.items[1].translator   = optional_field_translator ((uint8_t)value[0]); // instructions on how to transform array items if reconstructing as BAM (value[0] is the subtype of the array)
    con.items[1].separator[0] = optional_sep_by_type[IS_BAM][(uint8_t)value[0]];
    
    Context *element_ctx      = ctx_get_ctx (vb, con.items[1].dict_id);
    element_ctx->st_did_i     = container_ctx->did_i;
    element_ctx->flags.store  = optional_field_store_flag[(uint8_t)value[0]];

    for (con.repeats=0; con.repeats < CONTAINER_MAX_REPEATS && str_len > 0; con.repeats++) { // str_len will be -1 after last number

        const char *snip = str;
        for (; str_len && *str != ','; str++, str_len--) {};

        unsigned number_len = (unsigned)(str - snip);
        unsigned snip_len   = number_len; // might be changed by optimize
             
        char new_number_str[30];
        if (optimize && snip_len < 25)
            optimize (&snip, &snip_len, new_number_str);

        seg_by_ctx (VB, snip, snip_len, element_ctx, IS_BAM ? 0 : number_len+1);
        
        str_len--; // skip comma
        str++;
    }

    ASSSEG (con.repeats < CONTAINER_MAX_REPEATS, value, "array has too many elements, more than %u", CONTAINER_MAX_REPEATS);

    // add bytes here in case of BAM - all to main field
    unsigned container_add_bytes=0;
    if (IS_BAM) {
        unsigned add_bytes_per_repeat = sam_seg_optional_add_bytes (value[0], 0, true);
        container_add_bytes = add_bytes_per_repeat * con.repeats + 4 /* count */ + 1 /* type */ ;
    }
    else 
        container_add_bytes = 2; // type - eg "i,"

    container_seg (vb, container_ctx, (ContainerP)&con, prefixes, sizeof(prefixes), container_add_bytes);
}

// process an optional subfield, that looks something like MX:Z:abcdefg. We use "MX" for the field name, and
// the data is abcdefg. The full name "MX:Z:" is stored as part of the OPTIONAL dictionary entry
DictId sam_seg_optional_field (VBlockSAM *vb, ZipDataLineSAM *dl, bool is_bam, 
                               const char *tag, char bam_type, const char *value, unsigned value_len)
{
    char sam_type = sam_seg_bam_type_to_sam_type (bam_type);
    char dict_name[4] = { tag[0], tag[1], ':', sam_type };
    DictId dict_id = dict_id_make (dict_name, 4, DTYPE_SAM_OPTIONAL);

    unsigned add_bytes = sam_seg_optional_add_bytes (bam_type, value_len, is_bam);

    switch (dict_id.num) {

        case _OPTION_SA_Z: sam_seg_SA_field (vb, STRa(value)); break;

        case _OPTION_OA_Z: sam_seg_OA_field (vb, STRa(value)); break;

        case _OPTION_XA_Z: 
            if (segconf.running && !segconf.has_XA) segconf.has_XA = sam_seg_which_XA (STRa(value));

            if (segconf.has_XA == XA_BWA)
                sam_seg_BWA_XA_field (vb, STRa(value)); 
            else
                goto fallback;
            break;
        
        case _OPTION_MC_Z: sam_cigar_seg_MC (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_MD_Z: sam_md_seg (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_NM_i: sam_seg_NM_field (vb, STRa(value), add_bytes); break;

        case _OPTION_BD_Z:
        case _OPTION_BI_Z: sam_seg_BD_BI_field (vb, dl, STRa(value), dict_id, add_bytes); break;
        
        case _OPTION_AS_i: sam_seg_AS_field (vb, STRa(value), add_bytes); break;

        case _OPTION_XS_i: sam_seg_XS_field (vb, STRa(value), add_bytes); break;

        case _OPTION_XM_i: sam_seg_XM_field (vb, STRa(value), add_bytes); break;

        case _OPTION_MQ_i: sam_seg_MQ_field (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_mc_i: sam_seg_mc_field (vb, dict_id, STRa(value), add_bytes); break;

        case _OPTION_RG_Z: sam_seg_RG_field (vb, dl, STRa(value), add_bytes); break;

        // TX:i: - we seg this as a primary field SAM_TAX_ID
        case _OPTION_TX_i: seg_by_did_i (VB, taxid_redirection_snip, taxid_redirection_snip_len, OPTION_TX_i, add_bytes); break;

        //case _OPTION_E2: sam_seg_E2_field (vb, dl, STRa(value), add_bytes); // BROKEN. To do: fix.

        case _OPTION_U2_Z: sam_seg_U2_field (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_Z5_i: if (value_len && IS_DIGIT(value[0])) seg_pos_field (VB, OPTION_Z5_i, SAM_PNEXT, 0, 0, STRa(value), 0, add_bytes); 
                           else goto fallback; 
                           break;

        default:
            fallback:
            // integer
            if (sam_type == 'i')
                seg_integer_or_not (VB, ctx_get_ctx (vb, dict_id), STRa(value), add_bytes);

            // Numeric array array
            else if (bam_type == 'B') 
                sam_seg_array_field (VB, dict_id, STRa(value));

            // All other subfields - normal snips in their own dictionary
            else        
                seg_by_dict_id (VB, STRa(value), dict_id, add_bytes); 
    }

    // integer and float fields need to be STORE_INT/FLOAT to be reconstructable as BAM
    if (optional_field_store_flag[(uint8_t)sam_type]) {
        Context *ctx;
        if ((ctx = ECTX (dict_id)))
            ctx->flags.store = optional_field_store_flag[(uint8_t)sam_type];
    }
 
    return dict_id;
}
