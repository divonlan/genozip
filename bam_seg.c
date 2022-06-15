// ------------------------------------------------------------------
//   sam_bam.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "profiler.h"
#include "digest.h"
#include "buffer.h"
#include "vblock.h"
#include "txtfile.h"
#include "file.h"
#include "endianness.h"
#include "sam_private.h"
#include "seg.h"
#include "strings.h"
#include "random_access.h"
#include "dict_id.h"
#include "codec.h"
#include "flags.h"
#include "profiler.h"
#include "context.h"
#include "kraken.h"
#include "segconf.h"
#include "qname.h"
#include "libdeflate/libdeflate.h"

void bam_seg_initialize (VBlockP vb)
{
    sam_seg_initialize (vb);

    // we store special quality values at the end of txt_data, as they need to be somewhere in txt_data
    buf_alloc (vb, &vb->txt_data, 2, 0, char, 0, 0); // add 2 character after the end of txt_data
    BAFTtxt[0] = '*'; // QUAL_MISSING_STANDARD;
    BAFTtxt[1] = 127; // QUAL_MISSING_PYSAM;

    if (!segconf.running && line_textual_cigars_used)
        buf_alloc (vb, &VB_SAM->line_textual_cigars, 0, segconf.sam_cigar_len * vb->lines.len32 / (segconf.is_long_reads ? 4 : 1),/*divide in case sam_cigar_len is not representative*/
                   char, CTX_GROWTH, "line_textual_cigars");
}

static int32_t bam_unconsumed_scan_forwards (VBlockP vb)
{
    ARRAY (char, txt, vb->txt_data);
    
    if (txt_len < sizeof (BAMAlignmentFixed)) return -1; // this VB doesn't not even contain one single full alignment

    uint32_t aln_size=0, i;
    for (i=0 ; i < txt_len-3; i += aln_size) 
        aln_size = LTEN32 ((BAMAlignmentFixed *)&txt[i])->block_size + 4;

    if (aln_size > txt_len) 
        return -1; // this VB doesn't not even contain one single full alignment

    else if (i == txt_len)
        return 0;  // will will consume all data - nothing to pass to next VB
    
    else 
        return aln_size - (i - txt_len); // we pass the data of the final, partial, alignment to the next VB
}

static int32_t bam_unconsumed_scan_backwards (VBlockP vb, uint32_t first_i, int32_t *i)
{
    *i = MIN_(*i, vb->txt_data.len - sizeof(BAMAlignmentFixed));

    // find the first alignment in the data (going backwards) that is entirely in the data - 
    // we identify and alignment by l_read_name and read_name
    for (; *i >= (int32_t)first_i; (*i)--) {
        const BAMAlignmentFixed *aln = (const BAMAlignmentFixed *)Bc (vb->txt_data, *i);

        uint32_t block_size = LTEN32 (aln->block_size);
        if (block_size > 100000000) continue; // quick short-circuit - more than 100M for one alignment - clearly wrong

        uint32_t l_seq      = LTEN32 (aln->l_seq);
        uint16_t n_cigar_op = LTEN16 (aln->n_cigar_op);

        // test to see block_size makes sense
        if ((uint64_t)*i + (uint64_t)block_size + 4 > (uint64_t)vb->txt_data.len || // 64 bit arith to catch block_size=-1 that will overflow in 32b
            block_size + 4 < sizeof (BAMAlignmentFixed) + 4*n_cigar_op  + aln->l_read_name + l_seq + (l_seq+1)/2)
            continue;

        // test to see l_read_name makes sense
        if (LTEN32 (aln->l_read_name) < 2 ||
            &aln->read_name[aln->l_read_name] > BAFTtxt) continue;

        // test pos
        int32_t pos = LTEN32 (aln->pos);
        if (pos < -1) continue;

        // test read_name    
        if (aln->read_name[aln->l_read_name-1] != 0 || // nul-terminated
            !str_is_in_range (aln->read_name, aln->l_read_name-1, '!', '~')) continue;  // all printable ascii (per SAM spec)

        // final test option 1: test l_seq vs seq_len implied by cigar
        if (aln->l_seq && n_cigar_op) {
            uint32_t seq_len_by_cigar=0;
            uint32_t *cigar = (uint32_t *)((uint8_t *)(aln+1) + aln->l_read_name);
            for (uint16_t cigar_op_i=0; cigar_op_i < n_cigar_op; cigar_op_i++) {
                uint8_t cigar_op = *(uint8_t *)&cigar[cigar_op_i] & 0xf; // LSB by Little Endian - take 4 LSb
                uint32_t op_len = cigar[cigar_op_i] >> 4;
                if (cigar_lookup_bam[cigar_op] & CIGAR_CONSUMES_QUERY) seq_len_by_cigar += op_len; 
            }
            if (l_seq != seq_len_by_cigar) continue;
        }
        
        // final test option 2: in case we know the qname flavor
        else if (segconf.qname_flavor) {
            if (qname_test_flavor (aln->read_name, aln->l_read_name-1, segconf.qname_flavor, 0, 0)) // not qname_flavor
                continue;
        }

        // we don't have a final test - skip
        else continue;

        // Note: we don't use add aln->bin calculation because in some files we've seen data that doesn't
        // agree with our formula. see comment in bam_reg2bin

        // all tests passed - this is indeed an alignment
        return vb->txt_data.len - (*i + LTEN32 (aln->block_size) + 4); // everything after this alignment is "unconsumed"
    }

    return -1; // we can't find any alignment - need more data (lower first_i)    
}  

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
// if first_i > 0, we attempt to heuristically detect the start of a BAM alignment.
int32_t bam_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    ASSERT (*i >= 0 && *i < vb->txt_data.len, "*i=%d is out of range [0,%"PRIu64"]", *i, vb->txt_data.len);

    int32_t result;

    // if we have all the data - search forward - faster in BAM as we have aln.block_size
    if (first_i == 0)
        result = bam_unconsumed_scan_forwards (vb);

    // stringent -either CIGAR needs to match seq_len, or qname needs to match flavor
    else
        result = bam_unconsumed_scan_backwards (vb, first_i, i); 

    return result; // if -1 - we will be called again with more data
}

uint32_t bam_split_aux (VBlockSAMP vb, rom aux, rom after_aux, rom *auxs, uint32_t *aux_lens)
{
    uint32_t n_auxs = 0;
    while (aux < after_aux) {
        auxs[n_auxs] = aux;

        static unsigned const size[256] = { ['A']=1, ['c']=1, ['C']=1, ['s']=2, ['S']=2, ['i']=4, ['I']=4, ['f']=4 };

        if (aux[2] == 'Z' || aux[2] == 'H') {
            SAFE_NUL (after_aux);
            aux_lens[n_auxs] = strlen (aux+3) + 4; // add tag[2], type and \0
            SAFE_RESTORE;
        }
        else if (aux[2] == 'B') 
            aux_lens[n_auxs] = 8 + GET_UINT32 (aux+4) * size[(int)aux[3]]; 
        
        else if (size[(int)aux[2]])
            aux_lens[n_auxs] = 3 + size[(int)aux[2]];
        else
            ABORT ("%s: Unrecognized aux type '%c' (ASCII %u)", LN_NAME, aux[2], aux[2]);
        
        aux += aux_lens[n_auxs++];
    }

    ASSERT (aux == after_aux, "%s: overflow while parsing auxilliary fields", LN_NAME);

    return n_auxs;
}

void bam_seg_BIN (VBlockSAMP vb, ZipDataLineSAM *dl, uint16_t bin /* used only in bam */, bool is_bam)
{
    SamPosType this_pos = dl->POS;
    SamPosType last_pos = dl->FLAG.bits.unmapped ? this_pos : (this_pos + vb->ref_consumed - 1);
    uint16_t reg2bin = bam_reg2bin (this_pos, last_pos); // zero-based, half-closed half-open [start,end)

    if (!is_bam || (last_pos <= MAX_POS_SAM && reg2bin == bin))
        seg_by_did_i (VB, ((char []){ SNIP_SPECIAL, SAM_SPECIAL_BIN }), 2, SAM_BAM_BIN, is_bam ? sizeof (uint16_t) : 0);
    
    else {
#ifdef DEBUG // we show this warning only in DEBUG because I found actual files that have edge cases that don't work with our formula (no harm though)
        WARN_ONCE ("FYI: %s: bad bin value in: this_pos=%d ref_consumed=%u flag=%u last_pos=%d: bin=%u but reg2bin=%u. No harm. This warning will not be shown again for this file.",
                    LN_NAME, this_pos, vb->ref_consumed, dl->FLAG.value, last_pos, bin, reg2bin);
#endif
        seg_integer_as_text (vb, SAM_BAM_BIN, bin, is_bam);
        CTX(SAM_BAM_BIN)->flags.store = STORE_INT;
    }
}

static inline void bam_seg_ref_id (VBlockSAMP vb, ZipDataLineSAM *dl, DidIType did_i, int32_t ref_id, int32_t compare_to_ref_i)
{
    ASSERT (ref_id >= -1 && ref_id < (int32_t)sam_hdr_contigs->contigs.len, "%s: encountered ref_id=%d but header has only %u contigs",
            LN_NAME, ref_id, sam_hdr_contigs->contigs.len32);

    // get snip and snip_len
    STR0(snip);
    if (ref_id >= 0) {
        if (ref_id == compare_to_ref_i) {
            snip = "=";
            snip_len = 1;            
        }
        else 
            snip = contigs_get_name (sam_hdr_contigs, ref_id, &snip_len);
    }
    else { 
        snip = "*";
        snip_len = 1;
    }

    if (did_i == SAM_RNAME) sam_seg_RNAME (vb, dl, STRa(snip), true, sizeof (int32_t));
    else                    sam_seg_RNEXT (vb, dl, STRa(snip),       sizeof (int32_t));
}

// Rewrite the QUAL field - add +33 to Phred scores to make them ASCII
static inline bool bam_rewrite_qual (uint8_t *qual, uint32_t l_seq)
{
    if (qual[0] == 0xff) return false; // in case SEQ is present but QUAL is omitted, all qual is 0xFF

    for (uint32_t i=0; i < l_seq; i++)
        qual[i] += 33;

    return true;
}

static inline QualMissingType bam_get_missing_qual_type (VBlockSAMP vb, bytes qual, uint32_t l_seq)
{
    if (l_seq <= 1) return QUAL_MISSING_STANDARD;

    uint8_t filler = qual[1]; // SAM spec - needs to be 0xff. Genozip also supports - as created by pysam - 0x00.
    ASSERT (filler == 0 || filler == 0xff, "%s: Non-compliant QUAL field. l_seq=%u filler=%u (only 0 and 0xff are accepted)",
            LN_NAME, l_seq, filler);

    for (uint32_t i=2; i < l_seq; i++)
        ASSERT (qual[i] == filler, "%s: Non-compliant QUAL field. l_seq=%u filler=%u qual[%u]=%u",
                LN_NAME, l_seq, filler, i, qual[i]);
    
    return filler ? QUAL_MISSING_STANDARD : QUAL_MISSING_PYSAM;
}

void bam_get_one_aux (VBlockSAMP vb, STRp(aux),
                      rom *tag, char *type, char *array_subtype, // out
                      pSTRp(value), ValueType *numeric) // out - one of these depending on the type
{
    *tag  = aux;
    *type = aux[2]; // c, C, s, S, i, I, f, A, Z, H or B
    aux += 3;
    *array_subtype = 0;
    
    switch (*type) {
        // in case of an numeric type, we pass the value as a ValueType
        case 'i': *value_len = 4; numeric->i = (int32_t)LTEN32 (GET_UINT32 (aux)); break;
        case 'I': 
        case 'f': *value_len = 4; numeric->i = LTEN32 (GET_UINT32 (aux));          break; // note: uint32 and float are binary-identical so this effectively sets value->f            
        case 's': *value_len = 2; numeric->i = (int16_t)LTEN16 (GET_UINT16 (aux)); break;
        case 'S': *value_len = 2; numeric->i = LTEN16 (GET_UINT16 (aux));          break;
        case 'c': *value_len = 1; numeric->i = (int8_t)*aux;                       break;
        case 'C': *value_len = 1; numeric->i = (uint8_t)*aux;                      break;
        case 'Z': 
        case 'H': *value_len = aux_len - 4; *value = aux;                          break; // value_len excludes the terminating \0
        case 'A': *value_len = 1; *value = aux;                                    break;

        // in case of a numerical value we pass the data as is, in machine endianity
        case 'B':
            *array_subtype = *aux++; // type of elements of array
            *value_len = GET_UINT32(aux);  // number of elements, not number of bytes
            *value = aux + sizeof (uint32_t);
            
            // switch from BAM's little endian to machine endianity
            if (!flag.is_lten) 
                switch (*array_subtype) {
                    case 'i': case 'I': case 'f': 
                        for (int i=0; i < *value_len; i++) ((uint32_t*)(*value))[i] = LTEN32(((uint32_t*)(*value))[i]);
                        break;

                    case 's': case 'S': 
                        for (int i=0; i < *value_len; i++) ((uint32_t*)(*value))[i] = LTEN16(((uint32_t*)(*value))[i]);
                        break;

                    default: break;
                }

            ASSERT (aux_width[(uint8_t)*array_subtype], "%s: Invalid array type: '%c' for field \"%c%c\"", 
                    LN_NAME, *array_subtype, (*tag)[0], (*tag)[1]);

            break;

        default:
            ABORT ("%s: Invalid field type: '%c' for field \"%c%c\"", LN_NAME, *type, (*tag)[0], (*tag)[1]);
    }
}

static int64_t bam_get_one_aux_int (VBlockSAMP vb, STRp(aux))
{
    rom tag                  __attribute__((unused));
    char type, array_subtype __attribute__((unused));
    STR (value)              __attribute__((unused));
    ValueType numeric;

    bam_get_one_aux (vb, STRa(aux), &tag, &type, &array_subtype, pSTRa(value), &numeric);
    return numeric.i;
}

rom bam_seg_txt_line (VBlockP vb_, rom alignment /* BAM terminology for one line */,
                      uint32_t remaining_txt_len, bool *has_13_unused)   
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    sam_reset_line (VB);

    // case: --show - only print BAM lines, no segging
    if (flag.show_bam) return bam_show_line (vb, alignment, remaining_txt_len);

    ZipDataLineSAM *dl = DATA_LINE (vb->line_i);
    rom next_field = alignment;
    // BAM alignment fixed-length fields 
    uint32_t block_size = NEXT_UINT32;

    WordIndex prev_line_chrom = vb->chrom_node_index;

    // a non-sensical block_size might indicate an false-positive identification of a BAM alignment in bam_unconsumed
    ASSERT (block_size + 4 >= sizeof (BAMAlignmentFixed) && block_size + 4 <= remaining_txt_len, 
            "%s: (block_size+4)=%u is out of range - too small, or goes beyond end of txt data: remaining_txt_len=%u",
            LN_NAME, block_size+4, remaining_txt_len);

    rom after = alignment + block_size + sizeof (uint32_t);

    // note: in BAM, we calculate here instead of in seg_all_data_lines, bc we rewrite qual
    if (flag.debug_lines) 
        vb->debug_line_hash = adler32 (2, alignment, after - alignment); 

    dl->RNAME = vb->chrom_node_index = (int32_t)NEXT_UINT32; // corresponding to CHROMs in the BAM header. -1 in BAM means '*' (no RNAME) - which luckily is WORD_INDEX_NONE.    
    
    // expecting all contigs to be defined in SAM header, and hence in ol_nodes/ol_dict 
    ASSERT (vb->chrom_node_index >= -1 && vb->chrom_node_index < (int32_t)CTX(SAM_RNAME)->ol_nodes.len32, 
            "%s: RNAME=%d out of range [-1,%d]", LN_NAME, vb->chrom_node_index, (int)CTX(SAM_RNAME)->ol_nodes.len32-1);

    dl->POS              = 1 + (int32_t)NEXT_UINT32; // pos in BAM is 0 based, -1 for unknown 
    uint8_t l_read_name  = NEXT_UINT8;               // QNAME length
    dl->MAPQ             = NEXT_UINT8;
    uint16_t bin         = NEXT_UINT16;
    uint16_t n_cigar_op  = NEXT_UINT16;
    dl->FLAG.value       = NEXT_UINT16;              // not to be confused with our global var "flag"
    uint32_t l_seq       = NEXT_UINT32;              // note: we stick with the same logic as SAM for consistency - dl->SEQ.len is determined by CIGAR 
    dl->RNEXT            = (int32_t)NEXT_UINT32;     // corresponding to CHROMs in the BAM header
    SamPosType next_pos  = 1 + (int32_t)NEXT_UINT32; // pos in BAM is 0 based, -1 for unknown
    SamTlenType tlen     = (SamTlenType)NEXT_UINT32;
    rom read_name        = next_field;
    dl->QNAME            = (TxtWord){ .index = BNUMtxt (read_name), .len = l_read_name-1 }; // -1 don't count \0
    rom cigar            = read_name + l_read_name;
    bytes seq   = (uint8_t *)cigar + sizeof(uint32_t)*n_cigar_op;
    dl->SEQ.index        = BNUMtxt (seq);
    rom qual             = (rom)seq + (l_seq+1)/2;
    dl->QUAL             = (TxtWord){ .index = BNUMtxt (qual), .len = l_seq };
    rom aux              = qual + l_seq;
    
    vb->RNEXT_is_equal = (dl->RNAME == dl->RNEXT);
 
    // split auxillary fields
    rom auxs[MAX_FIELDS]; 
    uint32_t aux_lens[MAX_FIELDS];
    uint32_t n_auxs = bam_split_aux (vb, aux, after, auxs, aux_lens);
    sam_seg_index_aux (vb, STRas(aux));
    
    if (vb->chrom_node_index != WORD_INDEX_NONE) 
        ctx_get_vb_snip_ex (CTX(SAM_RNAME), vb->chrom_node_index, pSTRa(vb->chrom_name));

    else {
        vb->chrom_name = "";
        vb->chrom_name_len = 0;
    }

    // convert BAM seq format to SAM
    buf_alloc (vb, &vb->textual_seq, 0, l_seq+1 /* +1 for last half-byte */, char, 1.5, "textual_seq");
    bam_seq_to_sam (vb, seq, l_seq, false, true, &vb->textual_seq);

    // if this is a secondary / supplamentary read (aka Dependent) or a read that has an associated sec/sup read 
    // (aka Primary) - move the line to the appropriate component and skip it here (no segging done yet)
    if (vb->check_for_gc && sam_seg_is_gc_line (vb, dl, alignment, after - alignment, STRas(aux), true)) {
        vb->debug_line_hash_skip = true;
        goto done;
    }

    // case: to biopsy_line: we just needed to pass sam_seg_is_gc_line and we're done
    if (flag.biopsy_line.line_i != NO_LINE && sam_seg_test_biopsy_line (VB, alignment, block_size + 4)) 
        goto done;  

    // analyze (but not seg yet) cigar
    buf_add_more_(vb, &vb->binary_cigar, BamCigarOp, cigar, n_cigar_op, "binary_cigar");
    sam_cigar_binary_to_textual (vb, n_cigar_op, (uint32_t*)cigar, &vb->textual_cigar); // re-write BAM format CIGAR as SAM textual format in vb->textual_cigar

    uint32_t seq_len;
    bam_seg_cigar_analyze (vb, &seq_len);
    dl->SEQ.len = seq_len; // do it this way to avoid compiler warning

    // SEQ - calculate diff vs. reference (denovo or loaded)
    ASSERT (dl->SEQ.len == l_seq || (vb->textual_cigar.len == 1 && *B1STc(vb->textual_cigar) == '*') || !l_seq, 
            "seq_len implied by CIGAR=%s is %u, but actual SEQ length is %u, SEQ=%.*s", 
            B1STc(vb->textual_cigar), dl->SEQ.len, l_seq, l_seq, B1STc(vb->textual_seq));

    // convert BAM QUAL format to SAM
    if (!l_seq || !bam_rewrite_qual ((uint8_t *)qual, l_seq)) // add 33 to Phred scores to make them ASCII
        vb->qual_missing = bam_get_missing_qual_type (vb, (uint8_t *)qual, l_seq);

    if (has_NM) 
        dl->NM_len = sam_seg_get_aux_int (vb, STRi(aux, vb->idx_NM_i), &dl->NM, true, MIN_NM_i, MAX_NM_i, false);

    if (!sam_is_main_vb)
        sam_seg_sa_group_stuff (vb, dl , STRas(aux), STRb(vb->textual_cigar), B1STc(vb->textual_seq), true);

    // seg QNAME first, as it will find the buddy
    sam_seg_QNAME (vb, dl, read_name, l_read_name-1, 2); // QNAME. account for \0 and l_read_name

    bam_seg_ref_id (vb, dl, SAM_RNAME, vb->chrom_node_index, -1); // ref_id (RNAME)

    // note: pos can have a value even if ref_id=-1 (RNAME="*") - this happens if a SAM with a RNAME that is not in the header is converted to BAM with samtools
    sam_seg_POS (vb, dl, prev_line_chrom, sizeof (uint32_t)); // POS
    
    if (vb->chrom_node_index >= 0) sam_seg_verify_RNAME_POS (vb, NULL, dl->POS);

    sam_seg_MAPQ (vb, dl, sizeof (uint8_t));

    sam_seg_FLAG (vb, dl, sizeof (uint16_t));
    
    bam_seg_ref_id (vb, dl, SAM_RNEXT, dl->RNEXT, vb->chrom_node_index); // RNEXT

    sam_seg_PNEXT (vb, dl, 0, 0, next_pos, sizeof (uint32_t));

    // Segment BIN after we've gathered bin, flags, pos and vb->ref_confumed (and before sam_seg_SEQ which ruins vb->ref_consumed)
    bam_seg_BIN (vb, dl, bin, true);

    // we search forward for MD:Z now, XG:Z as we will need it for SEQ if it exists
    if (has_MD && !segconf.running)
        sam_seg_MD_Z_analyze (vb, dl, STRauxZ(MD_Z,true), dl->POS);

    if (has_XG)
        sam_seg_XG_Z_analyze (vb, dl, STRauxZ(XG_Z,true), dl->POS);

    // analyzing X0 - needed for segging XT:A and X1:i
    if (has_X0)
        ctx_set_last_value (VB, CTX(OPTION_X0_i), bam_get_one_aux_int (vb, STRi(aux, vb->idx_X0_i))); 

    // analyzing XA - needed for segging X1:i - set to number of repeats
    if (has_XA && has_X1)
        ctx_set_last_value (VB, CTX(OPTION_XA_Z), (int64_t)str_count_char (auxs[vb->idx_XA_Z]+3, aux_lens[vb->idx_XA_Z]-3, ';')); 

    sam_seg_SEQ (vb, dl, STRb(vb->textual_seq), (l_seq+1)/2 + sizeof (uint32_t) /* account for l_seq and seq fields */);

    // QUAL
    if (!vb->qual_missing) // case we have both SEQ and QUAL
        sam_seg_QUAL (vb, dl, qual, l_seq, l_seq /* account for qual field */ );

    else { // cases 1. were both SEQ and QUAL are '*' (seq_len=0) and 2. SEQ exists, QUAL not (bam_rewrite_qual returns false)
        *(char *)alignment = (vb->qual_missing == QUAL_MISSING_STANDARD) ? '*' : 127; // overwrite as we need it somewhere in txt_data
        dl->QUAL = (TxtWord){ .index=vb->txt_data.len + (vb->qual_missing == QUAL_MISSING_PYSAM), // special values set in bam_seg_initialize 
                              .len = 1 }; // index of "*" (or 127)

        sam_seg_QUAL (vb, dl, alignment, 1, l_seq /* account of l_seq 0xff */);
        
        vb->qual_codec_no_longr = true; // we cannot compress QUAL with CODEC_LONGR in this case
    }

    // finally we can segment the textual CIGAR now (including if n_cigar_op=0)
    sam_cigar_seg_binary (vb, dl, l_seq, cigar, n_cigar_op);

    // AUX fields - up to MAX_FIELDS of them
    sam_seg_aux_all (vb, dl, STRas(aux));
    
    // TLEN - must be after AUX as we might need data from MC:Z
    sam_seg_TLEN (vb, dl, 0, 0, tlen, vb->chrom_node_index == dl->RNEXT); // TLEN

done: {}

    buf_free (vb->textual_cigar);
    buf_free (vb->textual_seq);

    return after;
}
