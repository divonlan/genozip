// ------------------------------------------------------------------
//   sam_optimize.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "sam_private.h"

static rom sam_add_seq (VBlockSAMP vb, rom line_start, uint32_t remaining)
{
    ASSERT0 (!flag.optimize, "combining --optimize and --add-seq is not yet supported");

    bool has_13;
    rom next_line = line_start;
    str_split_by_tab (next_line, remaining, MAX_FIELDS + AUX, &has_13, false, true); // also advances next_line to next line
    
    ASSSEG (n_flds >= 10, "%s: Bad SAM file: alignment expected to have at least 10 fields (lacking SEQ), but found only %u", LN_NAME, n_flds);

    buf_alloc (vb, &vb->optimized_line, 0, (next_line - line_start) + (fld_lens[SEQ] + 1), char, 0, "optimized_line"); // x2+1000 is plenty for types of modifications we have so far.
    
    // the start of the line: until TLEN inc. the following \t
    char *next = mempcpy (B1STc(vb->optimized_line), line_start, flds[SEQ] - line_start); // initialize to exact copy of fields 1-10

    // generate SEQ
    for (int i=0; i < fld_lens[SEQ]; i++) // note: since there is no SEQ field, [SEQ] contains the QUAL data
        *next++ = 'A';

    CTX(SAM_SQBITMAP)->txt_shrinkage -= fld_lens[SEQ] + 1;

    // the end of the line: from QUAL inc. the preceding \t
    next = mempcpy (next, flds[SEQ] - 1/*\t*/, next_line - (flds[SEQ] - 1) + has_13);
    
    vb->optimized_line.len32 = BNUM(vb->optimized_line, next);
    return next_line;
}

void sam_segconf_finalize_optimizations (void)
{
    // optimize QUAL and other tags containing base qualities unless already binned (8 is the number of bins in Illimina: https://sapac.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf)
    if (segconf_get_num_qual_scores(QHT_QUAL) > 8) {
        segconf.optimize[SAM_QUAL] = true; 
        segconf.optimize[OPTION_QT_Z] = segconf.has[OPTION_QT_Z]; 
        segconf.optimize[OPTION_CY_Z] = segconf.has[OPTION_CY_Z]; 
        segconf.optimize[OPTION_BZ_Z] = segconf.has[OPTION_BZ_Z]; 
        segconf.optimize[OPTION_UY_Z] = segconf.has[OPTION_UY_Z] && segconf.has_10xGen; 
        segconf.optimize[OPTION_QX_Z] = segconf.has[OPTION_QX_Z] && segconf.has_10xGen; 
        segconf.optimize[OPTION_sQ_Z] = segconf.has[OPTION_sQ_Z] && segconf.has_10xGen; 
        segconf.optimize[OPTION_2Y_Z] = segconf.has[OPTION_2Y_Z] && segconf.has_10xGen; 
        segconf.optimize[OPTION_GY_Z] = segconf.has[OPTION_GY_Z] && segconf.has_10xGen; 
        segconf.optimize[OPTION_fq_Z] = segconf.has[OPTION_fq_Z] && segconf.has_10xGen; 
    }

    segconf.optimize[OPTION_ZM_B_s] = segconf.has[OPTION_ZM_B_s] && ((MP(TMAP/*mapped file*/) || MP(TORRENT_BC/*unmapped file*/)));

    // set float optimizations (note: all new contexts discovered by segconf were already added to z_file->contexts)
    for (Did did_i=SAM_FIRST_OPTIONAL_DID; did_i < z_file->num_contexts; did_i++) {
        #define ID(i,c) (ZCTX(did_i)->dict_id.id[i] == (c))
        if (segconf.has[did_i] && ((ID(2,':') && ID(3,'f') && ID(4,'\0')) || 
                                   (ID(2,':') && ID(3,'B') && ID(4,':') && ID(5,'f') && ID(6,'\0'))))
            segconf.optimize[did_i] = true;
    }
}

// SAM / BAM / FASTQ
// change the quality scores to be in a small number of bins, similar to Illumina: https://sapac.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf
// This is assuming standard Sanger format of Phred scores between 0 and 93 encoded in ASCII 33->126
// See here: https://pythonhosted.org/OBITools/fastq.html
// 0,1,2  -  unchanged : note regarding 2: according to the documentation, Illumina bins 2->6 but it does not modify scores of an 'N' base. In practice, in Illumina data, a 2 ('#') can only be a score of a N base, so it effectively also never changes.
// 3–9   -> 6
// 10–19 -> 15
// 20–24 -> 22
// 25–29 -> 27
// 30–34 -> 33
// 35–39 -> 37
// ≥ 40  -> 40
char *optimize_phred_quality_string (STRp(qual), char *out, bool is_bam,
                                     bool keep_underscore) 
{
#   define P(x) ((x)+33)
    static uint8_t qual_bins[256] = {
        // non Phread ASCII 0-32 - unchanged
        0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32, 
        
        // Phread values 0-93 -> ASCII 33-126
        P(0),  P(1),  P(2),  P(6),  P(6),  // Phred 0-4: same as Illumina
        P(6),  P(6),  P(6),  P(6),  P(6),  // Phred 5-9: same as Illumina
        P(15), P(15), P(15), P(15), P(15), // Phred 10-14: same as Illumina
        P(15), P(15), P(15), P(15), P(15), // Phred 14-19: same as Illumina
        P(22), P(22), P(22), P(22), P(22), // Phred 20-24: same as Illumina
        P(27), P(27), P(27), P(27), P(27), // Phred 25-29: same as Illumina
        P(33), P(33), P(33), P(33), P(33), // Phred 30-34: same as Illumina
        P(37), P(37), P(37), P(37), P(37), // Phred 35-39: same as Illumina
        P(42), P(42), P(42), P(42), P(42), // Phred 40-44
        P(47), P(47), P(47), P(47), P(47), // Phred 45-49
        P(52), P(52), P(52), P(52), P(52), // Phred 50-54
        P(57), P(57), P(57), P(57), P(57), // Phred 59-59
        P(62), P(62), P(62), P(62), P(62), // Phred 60-64
        P(67), P(67), P(67), P(67), P(67), // Phred 65-69
        P(72), P(72), P(72), P(72), P(72), // Phred 70-74
        P(77), P(77), P(77), P(77), P(77), // Phred 74-79
        P(82), P(82), P(82), P(82), P(82), // Phred 80-84
        P(87), P(87), P(87), P(87), P(87), // Phred 85-89
        P(91), P(91), P(91),               // Phred 90-92
        P(93),                             // Phred 93: keep maximum value as is (common in Pac Bio)

        // non Phread ASCII 127-255 - unchanged
        127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,
        157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,
        187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,
        217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,
        247,248,249,250,251,252,253,254,255
    };
#   undef P

    // do the mapping
    const uint8_t *after = (uint8_t *)qual + qual_len;
    if (is_bam)
        for (uint8_t *q=(uint8_t*)qual, *o=(uint8_t*)out; q < after; q++, o++) 
            *o = (char)qual_bins[(*q + 33) & 0xff] - 33;

    else  // SAM / FASTQ
        for (uint8_t *q=(uint8_t*)qual, *o=(uint8_t*)out; q < after; q++, o++) 
            *o = (char)qual_bins[*q];

    // revert changes to '_' if needed. separate loop for relatively rare case, to not burden the main loop
    // note: in STARsolo barcode fields with multiple parts, '_' is used as a separator (despite SAM spec recommending a space separator)
    if (keep_underscore) 
        for (uint8_t *q=(uint8_t*)qual, *o=(uint8_t*)out; q < after; q++, o++) 
            if (*q == (is_bam ? ('_'-33) : '_')) *o = *q;

    return out + qual_len;
}

// optimization for Ion Torrent flow signal (ZM:B:s) - negative values become zero, positives are rounded to the nearest 10
static void bam_optimize_TMAP_ZM (const BamArray_s *in, BamArray_s *out)
{
    for (uint32_t i=0; i < in->array_len; i++)
        out->elem[i] = (in->elem[i] >= 0) ? LTEN16 (((LTEN16(in->elem[i]) + 5) / 10) * 10) : 0;
}

static char *sam_optimize_TMAP_ZM (VBlockSAMP vb, STRp(in), char *out)
{
    char *save_out = out;
    out = mempcpy (out, in, 6);  // "ZM:B:s"

    str_split_ints (in+7, in_len-7, 0, ',', elem, false);

    for (uint32_t i=0; i < n_elems; i++) {
        *out++ = ',';
        // note: a out element can be 1 digit character more than the in: "99"->"100"
        out += str_int_ex ((elems[i] >= 0) ? ((elems[i] + 5) / 10) * 10 : 0, out, false);
    }

    CTX(OPTION_ZM_B_s)->txt_shrinkage += in_len - (out - save_out);
    return out;
}

static inline void optimize_float (const void *in, void *out) // non-aligned address of float in Little Endian
{
    uint32_t f32 = LTEN32 (GET_UINT32 (in)); // in native machine endianity
    // to preserve 10 significant bits we can get rid of the other 13 bits of the 23-bit mantissa, 
    // see: https://www.mimosa.org/ieee-floating-point-format/
    // TO DO: it would be better to round-up if the first dropped bit is 1, but it is not easy to do
    // as the exponent might change too. If we do so, we can make-out on more bit and achieve equivalent accuracy
    f32 &= 0b11111111111111111110000000000000;
    PUT_UINT32 (out, f32); // back to BAM's Little Endian
}

rom sam_zip_modify (VBlockP vb_, rom line_start, uint32_t remaining)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    
    if (flag.add_seq)   
        return sam_add_seq (vb, line_start, remaining);

    bool has_13;
    rom next_line = line_start;
    str_split_by_tab (next_line, remaining, MAX_FIELDS + AUX, &has_13, false, true); // also advances next_line to next line
    
    ASSSEG (n_flds >= 11, "%s: (sam_zip_modify) Bad SAM file: alignment expected to have at least 11 fields, but found only %u", LN_NAME, n_flds);

    uint32_t n_auxs = n_flds - AUX;
    rom *auxs = &flds[AUX]; // note: pointers to data on the stack
    uint32_t *aux_lens = &fld_lens[AUX];

    buf_alloc (vb, &vb->optimized_line, 0, (next_line - line_start) * 2 + 1000, char, 0, "optimized_line"); // x2+1000 is plenty for types of modifications we have so far.
    char *next = mempcpy (B1STc(vb->optimized_line), line_start, flds[QUAL] - line_start); // initialize to exact copy of fields 1-10

    next = segconf.optimize[SAM_QUAL] && !(TXT_DT(SAM) && str_issame_(STRfld(QUAL), "*", 1))
           ? optimize_phred_quality_string (STRfld(QUAL), next, false, false) 
           : mempcpy (next, flds[QUAL], fld_lens[QUAL]);
  
    *next++ = n_auxs ? '\t' : '\n'; // drop \r and drop terminal tab, if there were any

    bool terminal_tab = n_auxs && aux_lens[n_auxs-1] == 0;
    if (terminal_tab) n_auxs--;
    ASSSEG0 (!terminal_tab || !has_13, "line ends with \\t\\r\\n: this is not currently supported by Genozip");

    for (uint32_t f=0; f < n_auxs; f++) {
        DictId dict_id = auxs[f][3] != 'B' ? (DictId)DICT_ID_MAKE2_4(((char[]){ auxs[f][0], auxs[f][1], ':', auxs[f][3] }))
                                           : (DictId)DICT_ID_MAKE2_6(((char[]){ auxs[f][0], auxs[f][1], ':', auxs[f][3], ':', auxs[f][5] }));
        ContextP ctx;
        #define IF(cond) if ((ctx = ECTX(dict_id)) && segconf.optimize[ctx->did_i] && (cond)) 
        #define CASE(x) case x: if (!(ctx = ECTX(dict_id)) || !segconf.optimize[ctx->did_i]) goto fallback; 
        
        IF (auxs[f][3] == 'f') {
            next = mempcpy (next, auxs[f], 5);
            next = optimize_float_3_sig_dig (VB, ctx, auxs[f]+5, aux_lens[f]-5, next);
        }

        else IF (auxs[f][3] == 'B' && auxs[f][5] == 'f') {
            str_split (auxs[f]+7, aux_lens[f]-7, 0, ',', item, false);
            if (!n_items) goto fallback;

            next = mempcpy (next, auxs[f], 7);
            for (uint32_t i=0; i < n_items; i++) {
                next = optimize_float_3_sig_dig (VB, ctx, STRi(item,i), next);
                *next++ = ',';
            }
            next--; // remove final comma
        }
        
        else switch (dict_id.num) {    
            case _OPTION_QT_Z: case _OPTION_CY_Z: case _OPTION_BZ_Z: case _OPTION_UY_Z:
            case _OPTION_QX_Z: case _OPTION_sQ_Z: case _OPTION_2Y_Z: case _OPTION_GY_Z:
            CASE (_OPTION_fq_Z)
                next = mempcpy (next, auxs[f], 5);
                next = optimize_phred_quality_string (auxs[f]+5, aux_lens[f]-5, next, false, segconf.star_solo);
                break;

            CASE(_OPTION_ZM_B_s)
                next = sam_optimize_TMAP_ZM (vb, STRi(aux,f), next);
                break;

            default: fallback:
                next = mempcpy (next, auxs[f], aux_lens[f]);
            #undef CASE
        }
        

        *next++ = (f < n_auxs-1) ? '\t' : '\n'; // drop \r and drop terminal tab, if there were any
    }

    CTX(SAM_EOL)->txt_shrinkage += has_13;

    vb->optimized_line.len32 = BNUM(vb->optimized_line, next);
    return next_line;
}

rom bam_zip_modify (VBlockP vb_, rom line_start, uint32_t remaining)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    BAMAlignmentFixed *aln = (BAMAlignmentFixed *)line_start;

    rom after = line_start + LTEN32 (aln->block_size) + 4; // +4 for block_size field itself

    buf_alloc_exact (vb, vb->optimized_line, after - line_start, char, "optimized_line"); 
    memcpy (B1STc(vb->optimized_line), line_start, after - line_start); // initialize to exact copy
    
    uint32_t l_seq = LTEN32 (aln->l_seq);
    rom qual = aln->read_name + aln->l_read_name + LTEN16 (aln->n_cigar_op) * sizeof(uint32_t) + (l_seq+1)/2;

    if (segconf.optimize[SAM_QUAL] && l_seq && (uint8_t)qual[0] != 0xff) // in case SEQ is present but QUAL is omitted, all qual is 0xff
        optimize_phred_quality_string (qual, l_seq, Bc(vb->optimized_line, qual - line_start), true, false); 

    STR_ARRAY (aux,MAX_FIELDS) = bam_split_aux (vb, line_start, qual + l_seq, after, auxs, aux_lens);

    for (uint32_t f=0; f < n_auxs; f++) {
        DictId dict_id = auxs[f][2] != 'B' ? (DictId)DICT_ID_MAKE2_4(((char[]){ auxs[f][0], auxs[f][1], ':', auxs[f][2] }))
                                           : (DictId)DICT_ID_MAKE2_6(((char[]){ auxs[f][0], auxs[f][1], ':', auxs[f][2], ':', auxs[f][3] }));

        ContextP ctx;
        #define NEED_OPT ((ctx = ECTX(dict_id)) && segconf.optimize[ctx->did_i])
        #define CASE(x) case x: if (!NEED_OPT) break; 

        // optimize float fields
        if (auxs[f][2] == 'f') {  
            if (NEED_OPT) 
                optimize_float (auxs[f] + 3, Bc(vb->optimized_line, (auxs[f] + 3) - line_start));
        }

        // optimization float arrays
        else if (auxs[f][2] == 'B' && auxs[f][3] == 'f') {  
            if (NEED_OPT) {
                BamArray_f *in  = (BamArray_f *)(auxs[f]); // note: not endianity safe! TO DO: support Big Endian
                BamArray_f *out = (BamArray_f *)Bc(vb->optimized_line, (rom)in - line_start);
                for (uint32_t i=0; i < in->array_len; i++)
                    optimize_float (&in->elem[i], &out->elem[i]);
            }
        }

        else switch (dict_id.num) {    
            case _OPTION_QT_Z: case _OPTION_CY_Z: case _OPTION_BZ_Z: case _OPTION_UY_Z:
            case _OPTION_QX_Z: case _OPTION_sQ_Z: case _OPTION_2Y_Z: case _OPTION_GY_Z:
            CASE (_OPTION_fq_Z)
                optimize_phred_quality_string (auxs[f]+3, aux_lens[f]-3, Bc(vb->optimized_line, (auxs[f]+3) - line_start), true, segconf.star_solo); 
                break;

            CASE(_OPTION_ZM_B_s)
                const BamArray_s *in  = (BamArray_s *)(auxs[f]);
                BamArray_s *out = (BamArray_s *)Bc(vb->optimized_line, (rom)in - line_start);
                bam_optimize_TMAP_ZM (in, out);
                break;

            default: break;
            #undef CASE
        }
    }

    return after;
}
