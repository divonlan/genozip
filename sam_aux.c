// ------------------------------------------------------------------
//   sam_aux.c
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

static char taxid_redirection_snip[100], xa_strand_pos_snip[100];
static unsigned taxid_redirection_snip_len, xa_strand_pos_snip_len;

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
        ['Z']=CI_NATIVE_NEXT | CI_TRANS_NUL, ['H']=CI_NATIVE_NEXT | CI_TRANS_NUL, // reconstruct text and then \t seperator if SAM and \0 if BAM 
        ['A']=CI_NATIVE_NEXT,                                                     // reconstruct character and then \t seperator if SAM and no seperator for BAM
        ['B']=CI_NATIVE_NEXT                                                      // reconstruct array and then \t seperator if SAM and no seperator for BAM
}, 
{ // compressing from BAM
        ['c']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['C']=CI_NATIVE_NEXT | CI_TRANS_NOR, // reconstruct number and \t separator is SAM, and don't reconstruct anything if BAM (reconstruction will be done by translator)
        ['s']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['S']=CI_NATIVE_NEXT | CI_TRANS_NOR, // -"-
        ['i']=CI_NATIVE_NEXT | CI_TRANS_NOR, ['I']=CI_NATIVE_NEXT | CI_TRANS_NOR, // -"-
        ['f']=CI_NATIVE_NEXT,                                                     // compressing SAM - a float is stored as a SPECIAL, and the special reconstructor handles the SAM and BAM reconstructing
        ['Z']=CI_NATIVE_NEXT | CI_TRANS_NUL, ['H']=CI_NATIVE_NEXT | CI_TRANS_NUL, // reconstruct text and then \t seperator if SAM and \0 if BAM 
        ['A']=CI_NATIVE_NEXT,                                                     // reconstruct character and then \t seperator if SAM and no seperator for BAM
        ['B']=CI_NATIVE_NEXT                                                      // reconstruct array and then \t seperator if SAM and no seperator for BAM
} };

void sam_aux_zip_initialize (void)
{
    taxid_redirection_snip_len = sizeof (taxid_redirection_snip);
    seg_prepare_snip_other (SNIP_REDIRECTION, _SAM_TAXID, false, 0, 
                            taxid_redirection_snip, &taxid_redirection_snip_len);

    static SmallContainer xa_strand_pos_con = { .repeats=1, .nitems_lo=2, .items = { { { _OPTION_XA_STRAND } }, { { _OPTION_XA_POS } } } };
    xa_strand_pos_snip_len = sizeof (xa_strand_pos_snip);
    container_prepare_snip ((ConstContainerP)&xa_strand_pos_con, 0, 0, xa_strand_pos_snip, &xa_strand_pos_snip_len); 
}

// -------------------------------------------------------------------------------------------------------------------------------------------
// U2:Z "Phred probability of the 2nd call being wrong conditional on the best being wrong" (https://samtools.github.io/hts-specs/SAMtags.pdf)
// -------------------------------------------------------------------------------------------------------------------------------------------

// callback function for compress to get data of one line
void sam_zip_U2 (VBlock *vb, uint64_t vb_line_i, char **line_u2_data,  uint32_t *line_u2_len, uint32_t maximum_len) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    *line_u2_len = MIN_(maximum_len, dl->u2_data_len);

    if (!line_u2_data) return; // only lengths were requested

    *line_u2_data = ENT (char, vb->txt_data, dl->u2_data_start);

    if (flag.optimize_QUAL)
        optimize_phred_quality_string (*line_u2_data, *line_u2_len);
}

// ---------
// BD and BI
// ---------

static void sam_seg_BD_BI_field (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(field), DictId dict_id, unsigned add_bytes)
{
    bool is_bi = (dict_id.num == _OPTION_BI);
    Context *this_ctx  = is_bi ? CTX(OPTION_BI) : CTX (OPTION_BD);

    if (field_len != dl->seq_len) {
        seg_by_ctx (VB, field, field_len, this_ctx, field_len);
        return;
    }
    
    dl->bdbi_data_start[is_bi] = ENTNUM (vb->txt_data, field);

    CTX(OPTION_BD_BI)->txt_len += add_bytes; 

    if (!dl->bdbi_data_start[!is_bi]) // the first of BD and BI increments local.len, so it is incremented even if just one of BD/BI appears
        CTX(OPTION_BD_BI)->local.len += field_len * 2;

    seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, SAM_SPECIAL_BDBI }), 2, this_ctx, 0);
}

// callback function for compress to get BD_BI data of one line: this is an
// interlaced line containing a character from BD followed by a character from BI - since these two fields are correlated
// note: if only one of BD or BI exists, the missing data in the interlaced string will be 0 (this should is not expected to ever happen)
void sam_zip_BD_BI (VBlock *vb_, uint64_t vb_line_i, 
                    char **line_data, uint32_t *line_len,  // out 
                    uint32_t maximum_len)
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    
    const char *bd = dl->bdbi_data_start[0] ? ENT (char, vb->txt_data, dl->bdbi_data_start[0]) : NULL;
    const char *bi = dl->bdbi_data_start[1] ? ENT (char, vb->txt_data, dl->bdbi_data_start[1]) : NULL;
    
    if (!bd && !bi) return; // no BD or BI on this line

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_len  = MIN_(maximum_len, dl->seq_len * 2);

    if (!line_data) return; // only length was requested

    buf_alloc (vb, &vb->bd_bi_line, 0, dl->seq_len * 2, uint8_t, 2, "bd_bi_line");

    // calculate character-wise delta
    for (unsigned i=0; i < dl->seq_len; i++) {
        *ENT (uint8_t, vb->bd_bi_line, i*2    ) = bd ? bd[i] : 0;
        *ENT (uint8_t, vb->bd_bi_line, i*2 + 1) = bi ? bi[i] - (bd ? bd[i] : 0) : 0;
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

    if (ctx->dict_id.num == _OPTION_BD)
        for (uint32_t i=0; i < seq_len; i++, src+=2, dst++) *dst = *src;
    else
        for (uint32_t i=0; i < seq_len; i++, src+=2, dst++) *dst = *src + *(src+1);
    
    vb->txt_data.len += vb->seq_len;    
    ctx->next_local  += vb->seq_len;

done:
    return false; // no new value
}

// ----------------------------------------------------------------------------------------------
// OA:Z "Original alignment"
// SA:Z "Other canonical alignments in a chimeric alignment"
// XA:Z
// ----------------------------------------------------------------------------------------------

// OA and SA format is: (rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ . in OA - NM is optional (but its , is not)
// Example SA:Z:chr13,52863337,-,56S25M70S,0,0;chr6,145915118,+,97S24M30S,0,0;chr18,64524943,-,13S22M116S,0,0;chr7,56198174,-,20M131S,0,0;chr7,87594501,+,34S20M97S,0,0;chr4,12193416,+,58S19M74S,0,0;
// See: https://samtools.github.io/hts-specs/SAMtags.pdf
// note: even though SA, OA, XA contain similar fields amongst each other and similar to the primary fields,
// the values of subsequent lines tend to be similar for each one of them seperately, so we maintain separate contexts
static void sam_seg_SA_field (VBlockSAM *vb, STRp(field))
{
    static const MediumContainer container_SA = { .nitems_lo = 6,      
                                                  .repsep    = { ';' }, // including on last repeat    
                                                  .items     = { { .dict_id = { _OPTION_SA_RNAME  }, .seperator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_POS    }, .seperator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_STRAND }, .seperator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_CIGAR  }, .seperator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_MAPQ   }, .seperator = {','} },  
                                                                 { .dict_id = { _OPTION_SA_NM     },                  } } };

    SegCallback callbacks[6] = { [0]=chrom_seg_cb, [1]=seg_pos_field_cb, [3]=seg_add_to_local_text_cb };
    seg_array_of_struct (VB, CTX(OPTION_SA), container_SA, field, field_len, callbacks);
    CTX(OPTION_SA)->txt_len++; // 1 for \t in SAM and \0 in BAM 
}

static void sam_seg_OA_field (VBlockSAM *vb, STRp(field))
{
    static const MediumContainer container_OA = { .nitems_lo = 6,          
                                                  .repsep    = { ';' }, // including on last repeat    
                                                  .items     = { { .dict_id = { _OPTION_OA_RNAME  }, .seperator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_POS    }, .seperator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_STRAND }, .seperator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_CIGAR  }, .seperator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_MAPQ   }, .seperator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_NM     },                    } } };

    SegCallback callbacks[6] = { [0]=chrom_seg_cb, [1]=seg_pos_field_cb, [3]=seg_add_to_local_text_cb };
    seg_array_of_struct (VB, CTX(OPTION_OA), container_OA, field, field_len, callbacks);
    CTX(OPTION_OA)->txt_len++; // 1 for \t in SAM and \0 in BAM 
}

// split the pos strand-pos string, eg "-10000" to strand "-" and pos "10000"
static void seg_xa_strand_pos_cb (VBlockP vb, ContextP ctx, STRp(field))
{
    if (field_len < 2 || (field[0] != '+' && field[0] != '-'))  // invalid format - expecting pos to begin with the strand
        seg_by_ctx (VB, field, field_len, ctx, field_len);

    else {
        seg_by_did_i (VB, field, 1, OPTION_XA_STRAND, 1);
        seg_integer_or_not (vb, CTX(OPTION_XA_POS), &field[1], field_len-1, field_len-1);
        seg_by_ctx (VB, xa_strand_pos_snip, xa_strand_pos_snip_len, ctx, 0); // pre-created constant container
    }
}

static void sam_seg_XA_field (VBlockSAM *vb, STRp(field))
{
    // XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
    // Example XA:Z:chr9,-60942781,150M,0;chr9,-42212061,150M,0;chr9,-61218415,150M,0;chr9,+66963977,150M,1;
    // See: http://bio-bwa.sourceforge.net/bwa.shtml
    static const MediumContainer container_XA = {
        .repeats     = 0, 
        .nitems_lo   = 4, 
        .repsep      = {';'}, // including last item
        .items       = { { .dict_id = { _OPTION_XA_RNAME      }, .seperator = {','} }, // note: optional fields are DTYPE_2, in which short ids are left as-is, so we can skip dict_id_make
                         { .dict_id = { _OPTION_XA_STRAND_POS }, .seperator = {','} },
                         { .dict_id = { _OPTION_XA_CIGAR      }, .seperator = {','} }, // we don't mix the primary as the primary has a SNIP_SPECIAL
                         { .dict_id = { _OPTION_XA_NM         },                    } }  };

    SegCallback callbacks[4] = { [0]=chrom_seg_cb, [1]=seg_xa_strand_pos_cb, [2]=seg_add_to_local_text_cb };
    seg_array_of_struct (VB, CTX(OPTION_XA), container_XA, field, field_len, callbacks);
    CTX(OPTION_XA)->txt_len++; // 1 for \t in SAM and \0 in BAM 
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

    if (segconf.NM_is_integer && NM == vb->mismatch_bases)
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'i'}, 3, OPTION_NM, add_bytes); 

    else if (!segconf.NM_is_integer && (NM > 0) == (vb->mismatch_bases > 0))
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'b'}, 3, OPTION_NM, add_bytes); 

    else
        fallback:
        seg_by_did_i (VB, field, field_len, OPTION_NM, add_bytes); 
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


// ----------------------------------------------------------------------------------------------
// AS:i "Alignment score generated by aligner" (https://samtools.github.io/hts-specs/SAMtags.pdf)
// XS:i
// ----------------------------------------------------------------------------------------------

// AS has a value (at least as set by BWA) of at most seq_len, and often equal to it. we modify
// it to be new_value=(value-seq_len) 
static inline void sam_seg_AS_field (VBlockSAM *vb, ZipDataLineSAM *dl, 
                                     const char *snip, unsigned snip_len, unsigned add_bytes)
{
    int32_t as;
    bool positive_delta = str_get_int_range32 (STRa(snip), 0, 10000, &as);
    
    if (positive_delta) {
        ctx_set_last_value (VB, CTX (OPTION_AS), ((LastValueType){ .i = as }));
        if (dl->seq_len < as) positive_delta=false;
    }

    // if possible, store a special snip with the positive delta
    if (positive_delta) {
        char new_snip[20] = { SNIP_SPECIAL, SAM_SPECIAL_AS };
        unsigned delta_len = str_int (dl->seq_len-as, &new_snip[2]);

        seg_by_ctx (VB, new_snip, delta_len+2, CTX (OPTION_AS), add_bytes); 
    }

    // not possible - just store unmodified
    else
        seg_by_ctx (VB, snip, snip_len, CTX (OPTION_AS), add_bytes); 
}
/* doesn't help - variability vs AS is the same as independent variability 
static inline void sam_seg_XS_field (VBlockSAM *vb, ZipDataLineSAM *dl, 
                                     const char *snip, unsigned snip_len, unsigned add_bytes)
{
    int32_t xs;
    if (ctx_encountered_in_line_(VB, CTX(OPTION_AS)) && str_get_int_range32 (STRa(snip), 0, 10000, &xs)) {
        int32_t delta = xs - CTX(OPTION_AS)->last_value.i;

        SNIP(100);
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _OPTION_AS, true, delta, snip, &snip_len);
        seg_by_ctx (VB, STRa(snip), CTX(OPTION_XS), add_bytes);
    }
    else
        seg_by_ctx (VB, snip, snip_len, CTX (OPTION_AS), add_bytes); 
}*/

// MQ:i is often very similar to MAPQ
static inline void sam_seg_MQ_field (VBlockSAM *vb, ZipDataLineSAM *dl, 
                                     const char *snip, unsigned snip_len, unsigned add_bytes)
{
    int32_t mq, mapq;
    if (str_get_int_range32 (STRa(snip), 0, 1000, &mq) &&
        str_get_int_range32 (last_txt(vb, SAM_MAPQ), vb->last_txt_len(SAM_MAPQ), 0, 1000, &mapq)) {
        
        int32_t delta = mq - mapq;

        SNIP(100);
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _SAM_MAPQ, true, delta, snip, &snip_len);
        seg_by_ctx (VB, STRa(snip), CTX(OPTION_MQ), add_bytes);
    }
    else
        seg_by_ctx (VB, snip, snip_len, CTX (OPTION_MQ), add_bytes); 
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_AS)
{
    new_value->i = vb->seq_len - atoi (snip);
    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return true; // new value
}

// mc:i: (output of bamsormadup and other biobambam tools - mc in small letters) 
// appears to be a pos value usually close to PNEXT, but it is -1 is POS=PNEXT.
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

    sam_seg_seq_field (vb, OPTION_E2, (char *)field, field_len, this_pos, vb->last_cigar, 0, field_len, // remove const bc SEQ data is actually going to be modified
                        vb->last_cigar, add_bytes); 
}*/

// U2 - QUAL data (note: U2 doesn't have a context - it shares with QUAL)
static void sam_seg_U2_field (VBlockSAM *vb, ZipDataLineSAM *dl, STRp(field), unsigned add_bytes)
{
    ASSSEG0 (dl->seq_len, field, "U2 tag without a SEQ"); 
    ASSINP (field_len == dl->seq_len, 
            "Error in %s: Expecting U2 data to be of length %u as indicated by CIGAR, but it is %u. U2=%.*s",
            txt_name, dl->seq_len, field_len, field_len, field);

    dl->u2_data_start = field - vb->txt_data.data;
    dl->u2_data_len   = field_len;
    CTX(OPTION_U2)->txt_len   += add_bytes;
    CTX(OPTION_U2)->local.len += field_len;
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
    if (flag.optimize_ZM && dict_id.num == _OPTION_ZM && value_len > 3 && value[0] == 's')  // XM:B:s,
        optimize = sam_optimize_ZM;

    // prepare array container - a single item, with number of repeats of array element. array type is stored as a prefix
    Context *container_ctx     = ctx_get_ctx (vb, dict_id);

    SmallContainer con = { .nitems_lo = 2, 
                           .drop_final_item_sep_of_final_repeat = true, // TODO - get rid of this flag and move to making the seperators to be repeat seperators as they should have been, using drop_final_repeat_sep and obsoleting this flag 
                           .repsep    = {0,0}, 
                           .items     = { { .translator = SAM2BAM_ARRAY_SELF  },  // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM 
                                          { .seperator  = {0, ','}            } } // item[1] is actual array item
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
    con.items[1].seperator[0] = optional_sep_by_type[IS_BAM][(uint8_t)value[0]];
    
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

        case _OPTION_SA: sam_seg_SA_field (vb, value, value_len); break;

        case _OPTION_OA: sam_seg_OA_field (vb, value, value_len); break;

        case _OPTION_XA: sam_seg_XA_field (vb, value, value_len); break;

        // fields containing CIGAR format data - aliases of _OPTION_CIGAR (not the main CIGAR field that all snips have SNIP_SPECIAL)
        // MC: "Mate Cigar", added by eg https://manpages.debian.org/unstable/biobambam2/bamsort.1.en.html  
        case _OPTION_MC:
        case _OPTION_OC: seg_by_did_i (VB, value, value_len, OPTION_CIGAR, add_bytes); break;
        
        case _OPTION_MD: sam_md_seg (vb, dl, value, value_len, add_bytes); break;

        case _OPTION_NM: sam_seg_NM_field (vb, value, value_len, add_bytes); break;

        case _OPTION_BD:
        case _OPTION_BI: sam_seg_BD_BI_field (vb, dl, value, value_len, dict_id, add_bytes); break;
        
        case _OPTION_AS: sam_seg_AS_field (vb, dl, value, value_len, add_bytes); break;

        //case _OPTION_XS: sam_seg_XS_field (vb, dl, value, value_len, add_bytes); break;

        case _OPTION_MQ: sam_seg_MQ_field (vb, dl, value, value_len, add_bytes); break;

        case _OPTION_mc: sam_seg_mc_field (vb, dict_id, value, value_len, add_bytes); break;

        // TX:i: - we seg this as a primary field SAM_TAX_ID
        case _OPTION_TX: seg_by_did_i (VB, taxid_redirection_snip, taxid_redirection_snip_len, OPTION_TX, add_bytes); break;

        //case _OPTION_E2: sam_seg_E2_field (vb, dl, value, value_len, add_bytes); // BROKEN. To do: fix.

        case _OPTION_U2: sam_seg_U2_field (vb, dl, value, value_len, add_bytes); break;

        default:
            // Numeric array array
            if (bam_type == 'B') 
                sam_seg_array_field (VB, dict_id, value, value_len);

            // All other subfields - normal snips in their own dictionary
            else        
                seg_by_dict_id (VB, value, value_len, dict_id, add_bytes); 
    }

    // integer and float fields need to be STORE_INT/FLOAT to be reconstructable as BAM
    if (optional_field_store_flag[(uint8_t)sam_type]) {
        Context *ctx;
        if ((ctx = ECTX (dict_id)))
            ctx->flags.store = optional_field_store_flag[(uint8_t)sam_type];
    }
 
    return dict_id;
}
