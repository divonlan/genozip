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

#define DICT_ID_ARRAY(dict_id) (DictId){ .id = { (dict_id).id[0], (dict_id).id[1], '_','A','R','R','A','Y' } } // DTYPE_2

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
        seg_integer_as_text (vb, SAM_FLAG, dl->FLAG.value, false);
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

// -------------------------------------------------------------------------------------------------------------------------------------------
// U2:Z "Phred probability of the 2nd call being wrong conditional on the best being wrong" (https://samtools.github.io/hts-specs/SAMtags.pdf)
// -------------------------------------------------------------------------------------------------------------------------------------------

// callback function for compress to get data of one line
COMPRESSOR_CALLBACK (sam_zip_U2)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    *line_data_len = MIN_(maximum_size, dl->U2.snip_len);

    if (!line_data) return; // only lengths were requested

    *line_data = ENT (char, vb->txt_data, dl->U2.char_index);

    if (flag.optimize_QUAL)
        optimize_phred_quality_string (*line_data, *line_data_len);
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
COMPRESSOR_CALLBACK (sam_zip_BD_BI)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    
    const char *bd = dl->BD_BI[0].char_index ? ENT (char, vb->txt_data, dl->BD_BI[0].char_index) : NULL;
    const char *bi = dl->BD_BI[1].char_index ? ENT (char, vb->txt_data, dl->BD_BI[1].char_index) : NULL;
    
    if (!bd && !bi) return; // no BD or BI on this line

    ASSERT (bd && bi, "A line has one of the BD:Z/BI:Z pair - Genozip can only compress lines that have either both BD:Z and BI:Z or neither. vb=%u vb_line_i=%"PRIu64, 
            vb->vblock_i, vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_data_len  = MIN_(maximum_size, dl->seq_len * 2);

    if (!line_data) return; // only length was requested

    buf_alloc (vb, &VB_SAM->bd_bi_line, 0, dl->seq_len * 2, uint8_t, 2, "bd_bi_line");

    // calculate character-wise delta
    for (unsigned i=0; i < dl->seq_len; i++) {
        *ENT (uint8_t, VB_SAM->bd_bi_line, i*2    ) = bd[i];
        *ENT (uint8_t, VB_SAM->bd_bi_line, i*2 + 1) = bi[i] - (bd ? bd[i] : 0);
    }

    *line_data = FIRSTENT (char, VB_SAM->bd_bi_line);
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
    if (cigar_len > 7) {
        seg_add_to_local_text (vb, ctx, STRa(cigar), cigar_len);
        seg_simple_lookup (vb, ctx, 0);
    }
    else 
        seg_by_ctx (vb, STRa(cigar), ctx, cigar_len);

    return true;
}

// Lookup buffer
#define lookback_buf zip_lookback_buf // we store the previous rname, pos, strand in their ctx->zip_lookback_buf buffer
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
    WordIndex rname_index = LASTb250(rname_ctx);
    int64_t lookback = 0;
    PosType pos = -MAX_POS_DISTANCE; // initial to "invalid pos" - value chosen so we can storeit in poses, in lieu of an invalid non-integer value, without a future pos being considered close to it

    if (!segconf.running && str_get_int (STRa(pos_str), &pos)) {

        int64_t iterator = -1;
        while ((lookback = lookback_get_next (vb, rname_ctx, rname_index, &iterator))) {

            PosType lookback_pos = lookback_get_value (vb, pos_ctx, lookback);

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
        seg_by_ctx (vb, STRa(XA_lookback_snip), strand_ctx, 1);
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
     
    seg_array_of_struct (VB, CTX(OPTION_SA_Z), container_SA, field, field_len, callbacks);

    CTX(OPTION_SA_Z)->txt_len++; // 1 for \t in SAM and \0 in BAM 
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
     
    seg_array_of_struct (VB, CTX(OPTION_OA_Z), container_OA, field, field_len, callbacks);

    CTX(OPTION_OA_Z)->txt_len++; // 1 for \t in SAM and \0 in BAM 
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
static void sam_seg_NM_field (VBlockSAM *vb, ValueType NM, unsigned add_bytes)
{
    if (segconf.running && NM.i > 1) segconf.NM_is_integer = true; // we found evidence of integer NM

    ContextP ctx = CTX (OPTION_NM_i);
    ctx_set_last_value (VB, ctx, NM);

    if (segconf.NM_is_integer && NM.i == vb->mismatch_bases)
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'i'}, 3, OPTION_NM_i, add_bytes); 

    else if (!segconf.NM_is_integer && (NM.i > 0) == (vb->mismatch_bases > 0))
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'b'}, 3, OPTION_NM_i, add_bytes); 

    else 
        seg_integer (VB, ctx, NM.i, true, add_bytes);
}

SPECIAL_RECONSTRUCTOR (bam_piz_special_NM)
{
    if (*snip == 'i') 
        new_value->i = VB_SAM->mismatch_bases;

    else if (*snip == 'b')
        new_value->i = VB_SAM->mismatch_bases > 0;

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

static void sam_seg_XM_field (VBlockSAM *vb, ValueType XM, unsigned add_bytes)
{
    ContextP NM_ctx;
    
    // check for BWA case - XM is similar to NM - in our test file > 99% identical to NM.
    if (ctx_has_value_in_line (vb, _OPTION_NM_i, &NM_ctx) && XM.i == NM_ctx->last_value.i) 
        seg_by_did_i (VB, STRa(XM_snip), OPTION_XM_i, add_bytes); // copy from NM
    
    // check IonTorrent TMAP case - XM is supposed to be ref_consumed
    else if (XM.i == vb->ref_consumed)                              
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED }, 2, OPTION_XM_i, add_bytes);
        
    else
        seg_integer (VB, CTX(OPTION_XM_i), XM.i, true, add_bytes);
}

// ----------------------------------------------------------------------------------------------
// AS:i "Alignment score generated by aligner" (https://samtools.github.io/hts-specs/SAMtags.pdf)
// ----------------------------------------------------------------------------------------------

// AS has a value set (at least as set by BWA and IonTorrent TMAP) of at most vb->ref_consumed, and often equal to it. we modify
// it to be new_value=(value-ref_consumed) 
static inline void sam_seg_AS_field (VBlockSAM *vb, ValueType AS, unsigned add_bytes)
{
    ctx_set_last_value (VB, CTX (OPTION_AS_i), AS);

    // store a special snip with delta
    char new_snip[20] = { SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED };
    unsigned delta_len = str_int ((int32_t)vb->ref_consumed-AS.i, &new_snip[2]);

    seg_by_ctx (VB, new_snip, delta_len+2, CTX (OPTION_AS_i), add_bytes); 
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
    new_value->i = (int32_t)VB_SAM->ref_consumed - atoi (snip); // ref_consumed if snip=""
    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return true; // new value
}

// ----------------------------------------------------------------------------------------------
// XS:i Suboptimal alignment score 
// ----------------------------------------------------------------------------------------------

// alignment score of second best alignment - delta vs AS
static inline void sam_seg_XS_field (VBlockSAM *vb, ValueType XS, unsigned add_bytes)
{
    if (ctx_has_value_in_line_(VB, CTX(OPTION_AS_i)) && XS.i >= 0 && XS.i < 10000) {
        ctx_set_last_value (VB, CTX (OPTION_XS_i), XS); // needed for seg_delta_vs_other_do
        seg_delta_vs_other_do (VB, CTX(OPTION_XS_i), CTX(OPTION_AS_i), NULL, 0, -1, add_bytes);
    }
    else
        seg_integer_as_text_do (VB, CTX(OPTION_XS_i), XS.i, add_bytes); // seg as text to not prevent singletons
}

// MQ:i Mapping quality of the mate/next segment
// Seg against buddy if we have one, or else against MAPQ as it is often very similar
static inline void sam_seg_MQ_field (VBlockSAM *vb, ZipDataLineSAM *dl, ValueType MQ, unsigned add_bytes)
{
    if (segconf.running) segconf.has_MQ = true;

    dl->MQ = MQ.i;
    
    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1

    int64_t mapq = ctx_has_value_in_line_(vb, CTX(SAM_MAPQ)) ? CTX(SAM_MAPQ)->last_value.i : -1;

    if (!segconf.running && vb->buddy_line_i != -1 && MQ.i != mapq && MQ.i == buddy_dl->MAPQ)
        seg_by_did_i (VB, STRa(MQ_buddy_snip), OPTION_MQ_i, add_bytes); // copy MAPQ from earlier-line buddy 

    else if (mapq >= 0) {
        int64_t delta = MQ.i - mapq;

        SNIP(100);
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _SAM_MAPQ, true, delta, snip);
        seg_by_did_i (VB, STRa(snip), OPTION_MQ_i, add_bytes);
    }
    else
        seg_integer_as_text_do (VB, CTX(OPTION_MQ_i), MQ.i, add_bytes); // seg as text so to not prevent singletons
}

// mc:i: (output of bamsormadup and other biobambam tools - mc in small letters) 
// appears to be a pos value usually close to PNEXT, but it is -1 is POS=PNEXT.
// from bamsort manual: "adddupmarksupport=<0|1>: add information required for streaming duplicate marking in the aux fields MS and MC.
// Input is assumed to be collated by query name. This option is ignored unless fixmates=1. By default it is disabled."
// https://github.com/gt1/biobambam2/blob/master/src/programs/bamsort.cpp says: "biobambam used MC as a mate coordinate tag which now has a clash
// with the official SAM format spec.  New biobambam version uses mc."
// ms="MateBaseScore" - sum all characters in QUAL of the mate, where the value of each character is its ASCII minus 33 (i.e. the Phred score)
// mc="MateCoordinate"
static inline void sam_seg_mc_field (VBlockSAM *vb, ValueType mc, unsigned add_bytes)
{
    // if snip is "-1", store as simple snip
    if (mc.i == -1)
        seg_by_did_i (VB, "-1", 2, OPTION_mc_i, add_bytes);
    
    // delta vs PNEXT
    else
        seg_pos_field (VB, OPTION_mc_i, SAM_PNEXT, SPF_BAD_SNIPS_TOO, 0, 0, 0, mc.i, add_bytes);
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
static void sam_optimize_ZM (VBlockSAMP vb, ContextP ctx, void *cb_param, ValueType n, uint32_t repeat, unsigned add_bytes)
{
    ctx->dynamic_size_local = true;
    
    if (n.i >= 0) n.i = ((n.i + 5) / 10) * 10;
    else          n.i = 0;

    seg_integer (VB, ctx, n.i, false, add_bytes);
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
    if (!is_bam || type=='Z' || type=='H')
        return value_len + 1; // +1 for \0 (BAM Z/H) or \t (SAM)

    else
        return aux_width[(uint8_t)type]; // BAM

    // note: this will return 0 for 'B'
}

static inline SmallContainer *sam_seg_array_field_get_con (VBlockSAMP vb, Context *con_ctx, uint8_t type, bool has_callback,
                                                           ContextP *elem_ctx) // out
{
    // case: cached with correct type
    if (con_ctx->con_cache.param == type) {
        SmallContainer *con = FIRSTENT (SmallContainer, con_ctx->con_cache);
        *elem_ctx = ctx_get_ctx (vb, con->items[1].dict_id);
        return con; // note: we return a pointer into con_cache- sam_seg_array_field may modify the number of repeats. that's ok.
    }

    buf_alloc (vb, &con_ctx->con_cache, 0, 1, SmallContainer, 0, "con_cache");

    SmallContainer *con = FIRSTENT (SmallContainer, con_ctx->con_cache);
    *con = (SmallContainer){ .nitems_lo = 2, 
                             .drop_final_item_sep_of_final_repeat = true, // TODO - get rid of this flag and move to making the seperators to be repeat seperators as they should have been, using drop_final_repeat_sep and obsoleting this flag 
                             .repsep    = {0,0}, 
                             .items     = { { .translator = SAM2BAM_ARRAY_SELF  },  // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM 
                                            { .separator  = {0, ','}            } } // item[1] is actual array item
                           };            
    
    // prepare context where array elements will go in
    con->items[1].dict_id      = DICT_ID_ARRAY (con_ctx->dict_id);
    con->items[1].translator   = optional_field_translator (type); // instructions on how to transform array items if reconstructing as BAM (array[0] is the subtype of the array)
    con->items[1].separator[0] = optional_sep_by_type[IS_BAM][type];
    
    *elem_ctx = ctx_get_ctx (vb, con->items[1].dict_id);
    (*elem_ctx)->st_did_i = con_ctx->did_i;
    
    StoreType store_type = optional_field_store_flag[type];
    (*elem_ctx)->flags.store = store_type;

    ASSERT (store_type, "Invalid type \"%c\" in array of %s. vb=%u line=%"PRIu64, type, con_ctx->tag_name, vb->vblock_i, vb->line_i);

    if (store_type == STORE_INT && !has_callback) 
        (*elem_ctx)->dynamic_size_local = true; // Note: fields with callback should set ltype or dynamic in the callback if needed

    con_ctx->con_cache.param = type;

    return con;
}

static inline int64_t sam_seg_array_get_int (uint8_t type, const char *array)
{
    switch (type) {
        case 'I': return (uint32_t)GET_UINT32(array);
        case 'i': return (int32_t) GET_UINT32(array);
        case 'S': return (uint16_t)GET_UINT16(array);
        case 's': return (int16_t) GET_UINT16(array);
        case 'C': return (uint8_t) GET_UINT8 (array);
        case 'c': return (int8_t)  GET_UINT8 (array);
        default : ABORT_R ("Type is not int '%c'(%u)", type, type); 
    }
}

// an array - all elements go into a single item context, multiple repeats. items are segged as dynamic integers or floats, or a callback is called to seg them.
typedef void (*ArrayItemCallback) (VBlockSAMP vb, ContextP ctx, void *cb_param, ValueType n, uint32_t repeat, unsigned add_bytes);
static void sam_seg_array_field (VBlockSAMP vb, DictId dict_id, uint8_t type, 
                                 const char *array, int/*signed*/ array_len, // SAM: comma separated array ; BAM : arrays original width and machine endianity
                                 ArrayItemCallback callback, void *cb_param) // optional - call back for each item to seg the item
{   
    // prepare array container - a single item, with number of repeats of array element. array type is stored as a prefix
    Context *con_ctx = ctx_get_ctx (vb, dict_id), *elem_ctx;
    SmallContainer *con = sam_seg_array_field_get_con (vb, con_ctx, type, !!callback, &elem_ctx);

    int width = aux_width[type];
    bool is_int = (elem_ctx->flags.store == STORE_INT);

    ASSERT (is_int || (type=='f' && !callback), "Type not supported for SAM/BAM arrays '%c'(%u) in ctx=%s vb=%u line=%"PRIu64,
            type, type, con_ctx->tag_name, vb->vblock_i, vb->line_i);

    if (IS_BAM && !callback) {
        con->repeats = array_len;
        elem_ctx->txt_len += width * array_len;
        buf_alloc (vb, &elem_ctx->local, array_len, array_len * 50, int64_t, 2, "contexts->local"); // careful not * line.len - we can get OOM
        
        for (int i=0; i < array_len; i++, array += width) 
            if (is_int)
                NEXTENT (int64_t, elem_ctx->local) = sam_seg_array_get_int (type, array);
            else { // f
                char snip[16] = { SNIP_SPECIAL, SAM_SPECIAL_FLOAT };
                unsigned snip_len = 2 + str_int (GET_UINT32(array), &snip[2]);
                seg_by_ctx (VB, STRa(snip), elem_ctx, 0); // TODO: seg in local and update SPECIAL to get from local if no snip
            }
    }

    else if (IS_BAM && callback) {
        con->repeats = array_len;
        for (int i=0; i < array_len; i++, array += width) { 
            ValueType value = { .i = sam_seg_array_get_int (type, array) };
            callback (vb, elem_ctx, cb_param, value, i, width);
        }
    }

    else { // SAM

        // note: we're not using str_split on array, because the number of elements can be very large (eg one per base in PacBio ip:B) - possibly stack overflow
        for (con->repeats=0; con->repeats < CONTAINER_MAX_REPEATS && array_len > 0; con->repeats++) { // str_len will be -1 after last number

            const char *snip = array;
            for (; array_len && *array != ','; array++, array_len--) {};

            unsigned snip_len = (unsigned)(array - snip);

            if (elem_ctx->flags.store == STORE_INT) { 
                ValueType n;
                ASSERT (str_get_int (STRa(snip), &n.i), "Invalid array: \"%.*s\", expecting an integer in an array element of %s. vb=%u line=%"PRIu64, 
                        snip_len, snip, con_ctx->tag_name, vb->vblock_i, vb->line_i);
                
                if (callback)
                    callback (vb, elem_ctx, cb_param, n, con->repeats, snip_len+1);
                else
                    seg_integer (VB, elem_ctx, n.i, false, snip_len+1);
            }
            else {
                // TODO: change to seg_float_or_not - need to implement reconstruction, see comment in seg_float_or_not
                //seg_float_or_not (vb, elem_ctx, STRa(snip), snip_len+1);
                seg_by_ctx (VB, STRa(snip), elem_ctx, snip_len+1);
            }
            
            array_len--; // skip comma
            array++;
        }
    }

    ASSERT (con->repeats < CONTAINER_MAX_REPEATS, "array has too many elements, more than %u", CONTAINER_MAX_REPEATS);

    // add bytes here in case of BAM - all to main field
    unsigned container_add_bytes = IS_BAM ? 4/*count*/ + 1/*type*/ : 2/*type - eg "i,"*/;
    container_seg (vb, con_ctx, (ContainerP)con, ((char[]){ CON_PX_SEP, type, ',', CON_PX_SEP }), 4, container_add_bytes);
}

// process an optional subfield, that looks something like MX:Z:abcdefg. We use "MX" for the field name, and
// the data is abcdefg. The full name "MX:Z:" is stored as part of the OPTIONAL dictionary entry
DictId sam_seg_optional_field (VBlockSAM *vb, ZipDataLineSAM *dl, bool is_bam, 
                               const char *tag, char bam_type, char array_subtype, 
                               STRp(value), ValueType numeric) // two options 
{
    char sam_type = sam_seg_bam_type_to_sam_type (bam_type);
    char dict_name[4] = { tag[0], tag[1], ':', sam_type };
    DictId dict_id = dict_id_make (dict_name, 4, DTYPE_SAM_OPTIONAL); // match dict_id as declared in #pragma GENDICT

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

        case _OPTION_NM_i: sam_seg_NM_field (vb, numeric, add_bytes); break;

        case _OPTION_BD_Z:
        case _OPTION_BI_Z: sam_seg_BD_BI_field (vb, dl, STRa(value), dict_id, add_bytes); break;
        
        case _OPTION_AS_i: sam_seg_AS_field (vb, numeric, add_bytes); break;

        case _OPTION_XS_i: sam_seg_XS_field (vb, numeric, add_bytes); break;

        case _OPTION_XM_i: sam_seg_XM_field (vb, numeric, add_bytes); break;

        case _OPTION_MQ_i: sam_seg_MQ_field (vb, dl, numeric, add_bytes); break;

        case _OPTION_mc_i: sam_seg_mc_field (vb, numeric, add_bytes); break;

        case _OPTION_ms_i: sam_seg_ms_field (vb, numeric, add_bytes); break;

        case _OPTION_RG_Z: sam_seg_RG_field (vb, dl, STRa(value), add_bytes); break;

        // tx:i: - we seg this as a primary field SAM_TAX_ID
        case _OPTION_tx_i: seg_by_did_i (VB, taxid_redirection_snip, taxid_redirection_snip_len, OPTION_tx_i, add_bytes); break;

        //case _OPTION_E2: sam_seg_E2_field (vb, dl, STRa(value), add_bytes); // BROKEN. To do: fix.

        case _OPTION_U2_Z: sam_seg_U2_field (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_Z5_i: seg_pos_field (VB, OPTION_Z5_i, SAM_PNEXT, 0, 0, 0, 0, numeric.i, add_bytes); break;

        case _OPTION_ZM_B: sam_seg_array_field (vb, (DictId)_OPTION_ZM_B, array_subtype, STRa(value), 
                                                (segconf.tech == TECH_IONTORR && flag.optimize_ZM) ? sam_optimize_ZM : NULL, 0);
                           break;

        default: fallback:
            
            // all types of integer
            if (sam_type == 'i') {
                ContextP ctx = ctx_get_ctx (vb, dict_id);
                ctx->dynamic_size_local = true;
                ctx->flags.store = STORE_INT; // needs this be reconstructable as BAM
                seg_integer (VB, ctx, numeric.i, true, add_bytes);
            }

            else if (sam_type == 'f') {
                ContextP ctx = ctx_get_ctx (vb, dict_id);
                ctx->flags.store = STORE_FLOAT; // needs this be reconstructable as BAM

                if (is_bam) {
                    char snip[16] = { SNIP_SPECIAL, SAM_SPECIAL_FLOAT };
                    unsigned snip_len = 2 + str_int (numeric.i, &snip[2]); // we emit the textual integer that is binary-identical to the float, when converted to binary for
                    seg_by_ctx (VB, STRa(snip), ctx, add_bytes); 
                }
                else
                    seg_by_ctx (VB, STRa(value), ctx, add_bytes); 
            }

            // Numeric array
            else if (sam_type == 'B') 
                sam_seg_array_field (vb, dict_id, array_subtype, STRa(value), NULL, NULL);

            // Z,H,A - normal snips in their own dictionary
            else        
                seg_by_dict_id (VB, STRa(value), dict_id, add_bytes); 
    }
    
    return dict_id;
}
